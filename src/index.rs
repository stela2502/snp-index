//! Binned SNP index with `FeatureIndex` support.
//!
//! Design rules:
//! - `SnpIndex::loci` is the only canonical SNP storage.
//! - Bins only store locus ids into `loci`.
//! - Chromosome names are resolved through a fuzzy lookup table.
//! - Read/SNP overlap returns `ObservedSnp`, not pre-classified monster structs.

use crate::locus::SnpLocus;
use crate::read::{AlignedRead, ObservedBase};
use crate::vcf::{SnpVcfReader, VcfReadOptions};
use anyhow::{anyhow, Context, Result};
use scdata::feature_index::FeatureIndex;
use std::collections::HashMap;
use std::path::Path;
use crate::RawSnpRecord;
use crate::Strand;
use std::fmt;
use crate::ReadOpKind::*;

pub const DEFAULT_BIN_WIDTH: usize=10_000;

/// Binned SNP index.
///
/// SNPs are stored once in [`SnpIndex::loci`].
/// Bins store only ids into that vector.
#[derive(Debug, Clone)]
pub struct SnpIndex {
    pub chr_info: Vec<ChrInfo>,
    pub chr_name_to_id: HashMap<String, usize>,
    pub bin_width: u32,

    /// Canonical SNP table.
    ///
    /// Invariant:
    /// `loci[id].id == id`
    pub loci: Vec<SnpLocus>,

    /// Flat bin vector.
    ///
    /// Global bin:
    /// `chr_info[chr_id].bin_offset + local_bin`
    pub bins: Vec<SnpLocusBin>,

    /// Feature name / VCF id -> feature id.
    pub name_to_id: HashMap<String, u64>,
}

/// Per-chromosome metadata.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ChrInfo {
    pub name: String,
    pub length: u32,
    pub bin_offset: usize,
    pub bin_count: usize,
}

/// One genomic bin.
///
/// This intentionally does not know how to search.
/// Search belongs to [`SnpIndex`].
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct SnpLocusBin {
    /// Indices into `SnpIndex::loci`.
    ///
    /// Invariant:
    /// sorted by `(loci[id].pos0, loci[id].id)`.
    pub locus_ids: Vec<usize>,
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct SnpReadMatch<'a> {
    pub reference: Vec<ObservedSnp<'a>>,
    pub alternate: Vec<ObservedSnp<'a>>,
    pub other: Vec<ObservedSnp<'a>>,
}

impl<'a> SnpReadMatch<'a> {
    pub fn push(&mut self, hit: ObservedSnp<'a>) {
        if hit.is_ref() {
            self.reference.push(hit);
        } else if hit.is_alt() {
            self.alternate.push(hit);
        } else {
            self.other.push(hit);
        }
    }

    pub fn is_empty(&self) -> bool {
        self.reference.is_empty() && self.alternate.is_empty() && self.other.is_empty()
    }

    pub fn has_reference(&self) -> bool {
        !self.reference.is_empty()
    }

    pub fn has_alternate(&self) -> bool {
        !self.alternate.is_empty()
    }

    pub fn reference_ids(&self) -> impl Iterator<Item = u64> + '_ {
        self.reference.iter().map(|hit| hit.feature_id())
    }

    pub fn alternate_ids(&self) -> impl Iterator<Item = u64> + '_ {
        self.alternate.iter().map(|hit| hit.feature_id())
    }

    pub fn other_ids(&self) -> impl Iterator<Item = u64> + '_ {
        self.other.iter().map(|hit| hit.feature_id())
    }
}

/// SNP observed in one aligned read.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ObservedSnp<'a> {
    pub locus: &'a SnpLocus,
    pub observed: ObservedBase,
}

impl<'a> ObservedSnp<'a> {
    pub fn is_ref(&self) -> bool {
        self.observed.base.to_ascii_uppercase() == self.locus.reference.to_ascii_uppercase()
    }

    pub fn is_alt(&self) -> bool {
        let base = self.observed.base.to_ascii_uppercase();

        self.locus
            .alternates
            .iter()
            .any(|alt| base == alt.to_ascii_uppercase())
    }

    pub fn is_other(&self) -> bool {
        !self.is_ref() && !self.is_alt()
    }

    pub fn feature_id(&self) -> u64 {
        self.locus.id as u64
    }
}

impl SnpIndex {
    /// Build a SNP index from already parsed loci.
    ///
    /// `chr_names` and `chr_lengths` define the chromosome id space.
    /// Each `SnpLocus::chr_id` must refer to this space.
    pub fn new(
        chr_names: &[String],
        chr_lengths: &[u32],
        mut raw_loci: Vec<RawSnpRecord>,
        bin_width: u32,
    ) -> Result<Self> {
        Self::validate_genome_inputs(&chr_names, &chr_lengths)?;
        Self::validate_bin_width(bin_width)?;

        let chr_name_to_id = Self::build_chr_map(&chr_names);
        let chr_info = Self::build_chr_info(chr_names, chr_lengths, bin_width)?;

        Self::validate_and_canonicalize_raw_loci(&chr_info, &mut raw_loci)?;

        raw_loci.sort_by(|a, b| {
            (a.chr_id, a.pos0, a.name.as_str(), a.vcf_id.as_str())
                .cmp(&(b.chr_id, b.pos0, b.name.as_str(), b.vcf_id.as_str()))
        });

        let loci: Vec<SnpLocus> = raw_loci
            .into_iter()
            .enumerate()
            .map(|(id, raw)| SnpLocus {
                id,
                chr_id: raw.chr_id,
                pos0: raw.pos0,
                reference: raw.reference.to_ascii_uppercase(),
                alternates: raw
                    .alternates
                    .into_iter()
                    .map(|b| b.to_ascii_uppercase())
                    .collect(),
                name: raw.name.clone(),
                vcf_id: raw.vcf_id.clone(),
            })
            .collect();

        let mut bins = Self::build_empty_bins(&chr_info);
        Self::fill_bins(&chr_info, bin_width, &loci, &mut bins)?;
        Self::sort_bins(&loci, &mut bins);

        let name_to_id = Self::build_name_to_id(&loci);

        Ok(Self {
            chr_info,
            chr_name_to_id,
            bin_width,
            loci,
            bins,
            name_to_id,
        })
    }

    /// Build the index directly from a VCF path.
    ///
    /// Adapt the two `SnpVcfReader` calls if your reader API uses different names.
    pub fn from_vcf_path<P: AsRef<Path>>(
        path: P,
        chr_names: Vec<String>,
        chr_lengths: Vec<u32>,
        bin_width: u32,
        options: &VcfReadOptions,
    ) -> Result<Self> {
        let chr_name_to_id = Self::build_chr_map(&chr_names);

        let raw_loci =  SnpVcfReader::read_path(path, &chr_name_to_id, options)
            .context("failed to read SNP loci from VCF")?;

        Self::new(&chr_names, &chr_lengths, raw_loci, bin_width)
    }

    /// Fuzzy chromosome name lookup.
    pub fn chr_id(&self, chr_name: &str) -> Option<usize> {
        self.chr_name_to_id.get(chr_name).copied()
    }

    pub fn chr_name(&self, chr_id: usize) -> Option<&str> {
        self.chr_info.get(chr_id).map(|chr| chr.name.as_str())
    }

    /// Return global bin id for one chromosome position.
    pub fn global_bin_for_pos(&self, chr_id: usize, pos0: u32) -> Option<usize> {
        let chr = self.chr_info.get(chr_id)?;

        if pos0 >= chr.length {
            return None;
        }

        let local_bin = (pos0 / self.bin_width) as usize;

        if local_bin >= chr.bin_count {
            return None;
        }

        Some(chr.bin_offset + local_bin)
    }

    pub fn global_bin_for_chr_name_pos(&self, chr_name: &str, pos0: u32) -> Option<usize> {
        let chr_id = self.chr_id(chr_name)?;
        self.global_bin_for_pos(chr_id, pos0)
    }

    pub fn locus_ids_in_global_bin(&self, global_bin: usize) -> Option<&[usize]> {
        self.bins
            .get(global_bin)
            .map(|bin| bin.locus_ids.as_slice())
    }

    /// SNPs in one exact position.
    pub fn snps_at_pos(&self, chr_id: usize, pos0: u32) -> SnpPosIter<'_> {
        let Some(global_bin) = self.global_bin_for_pos(chr_id, pos0) else {
            return SnpPosIter::empty(self);
        };

        let Some(bin) = self.bins.get(global_bin) else {
            return SnpPosIter::empty(self);
        };

        let ids = bin.locus_ids.as_slice();

        let start = ids.partition_point(|&id| self.loci[id].pos0 < pos0);
        let end = ids.partition_point(|&id| self.loci[id].pos0 <= pos0);

        SnpPosIter {
            index: self,
            ids: &ids[start..end],
            offset: 0,
        }
    }

    pub fn snps_at_chr_name_pos(
        &self,
        chr_name: &str,
        pos0: u32,
    ) -> SnpPosIter<'_> {
        let Some(chr_id) = self.chr_id(chr_name) else {
            return SnpPosIter::empty(self);
        };

        self.snps_at_pos(chr_id, pos0)
    }

    /// SNPs in closed-open interval `[start0, end0)`.
    pub fn snps_in_range(
        &self,
        chr_id: usize,
        start0: u32,
        end0: u32,
    ) -> impl Iterator<Item = &SnpLocus> {
        if end0 <= start0 {
            return SnpRangeIter::empty(self);
        }

        let Some(chr) = self.chr_info.get(chr_id) else {
            return SnpRangeIter::empty(self);
        };

        if start0 >= chr.length {
            return SnpRangeIter::empty(self);
        }

        let end0 = end0.min(chr.length);

        let first_bin = (start0 / self.bin_width) as usize;
        let last_bin = ((end0 - 1) / self.bin_width) as usize;

        SnpRangeIter {
            index: self,
            start0,
            end0,
            current_global_bin: chr.bin_offset + first_bin,
            last_global_bin: chr.bin_offset + last_bin,
            current_ids: &[],
            current_offset: 0,
        }
    }


    /// Return SNP loci actually observed by a read.
    ///
    /// This does not classify ref/alt/other globally.
    /// Each returned [`ObservedSnp`] can classify itself.
    pub fn observed_snps<'a>(
        &'a self,
        read: &'a AlignedRead,
        min_baseq: u8,
    ) -> Vec<ObservedSnp<'a>> {
        let mut out = Vec::new();

        let Some((start0, end0)) = read.ref_span() else {
            return out;
        };

        for locus in self.snps_in_range(read.chr_id, start0, end0) {
            let Some(mut observed) = read.base_at_ref_pos(locus.pos0) else {
                continue;
            };

            if let Some(q) = observed.qual {
                if q < min_baseq {
                    continue;
                }
            }

            out.push(ObservedSnp { locus, observed });
        }

        out
    }

    /// Return a tupel of SNP_ids:
    ///
    /// ( ref, alt, other )
    pub fn get_ref_alt_other_ids_for_read( &self, read: &AlignedRead, min_baseq: u8,) -> (Vec<u64>, Vec<u64>, Vec<u64>){
        let obs = self.observed_snps( &read, min_baseq );
        let mut reference = Vec::with_capacity(obs.len());
        let mut alternate =  Vec::with_capacity(obs.len());
        let mut other =     Vec::with_capacity(obs.len());

        for obs_snp in obs {
            let feature_id = obs_snp.feature_id();
            if obs_snp.is_ref() {
                reference.push(feature_id);
            } else if obs_snp.is_alt() {
                alternate.push(feature_id);
            } else {
                other.push(feature_id);
            }
        }
        (reference, alternate, other)
    }

    fn complement(b: u8) -> u8 {
        match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'N',
        }
    }
    /// Build a fuzzy chromosome name map.
    ///
    /// Examples:
    /// - `chr1` maps also from `1`
    /// - `1` maps also from `chr1`
    /// - `MT`, `M`, `chrM` are treated as aliases
    pub fn build_chr_map(chr_names: &[String]) -> HashMap<String, usize> {
        let mut map = HashMap::with_capacity(chr_names.len() * 4);

        for (chr_id, name) in chr_names.iter().enumerate() {
            map.entry(name.clone()).or_insert(chr_id);

            if let Some(no_chr) = name.strip_prefix("chr") {
                map.entry(no_chr.to_string()).or_insert(chr_id);
            } else {
                map.entry(format!("chr{name}")).or_insert(chr_id);
            }

            match name.as_str() {
                "MT" => {
                    map.entry("chrM".to_string()).or_insert(chr_id);
                    map.entry("M".to_string()).or_insert(chr_id);
                }
                "chrM" => {
                    map.entry("MT".to_string()).or_insert(chr_id);
                    map.entry("M".to_string()).or_insert(chr_id);
                }
                "M" => {
                    map.entry("MT".to_string()).or_insert(chr_id);
                    map.entry("chrM".to_string()).or_insert(chr_id);
                }
                _ => {}
            }
        }

        map
    }

    fn validate_genome_inputs(chr_names: &[String], chr_lengths: &[u32]) -> Result<()> {
        if chr_names.len() != chr_lengths.len() {
            return Err(anyhow!(
                "chr_names length ({}) != chr_lengths length ({})",
                chr_names.len(),
                chr_lengths.len()
            ));
        }

        if chr_names.is_empty() {
            return Err(anyhow!("cannot build SnpIndex without chromosomes"));
        }

        for (i, name) in chr_names.iter().enumerate() {
            if name.is_empty() {
                return Err(anyhow!("chromosome name at index {i} is empty"));
            }
        }

        for (name, length) in chr_names.iter().zip(chr_lengths.iter()) {
            if *length == 0 {
                return Err(anyhow!("chromosome '{name}' has length 0"));
            }
        }

        Ok(())
    }

    fn validate_bin_width(bin_width: u32) -> Result<()> {
        if bin_width == 0 {
            return Err(anyhow!("bin_width must be > 0"));
        }

        Ok(())
    }

    fn build_chr_info(
        chr_names: &[String],
        chr_lengths: &[u32],
        bin_width: u32,
    ) -> Result<Vec<ChrInfo>> {
        let mut out = Vec::with_capacity(chr_names.len());
        let mut bin_offset = 0usize;

        for (name, length) in chr_names.into_iter().zip(chr_lengths.into_iter()) {
            let bin_count = length.div_ceil(bin_width) as usize;

            out.push(ChrInfo {
                name: name.to_string(),
                length: *length,
                bin_offset,
                bin_count,
            });

            bin_offset = bin_offset
                .checked_add(bin_count)
                .ok_or_else(|| anyhow!("too many SNP bins; usize overflow"))?;
        }

        Ok(out)
    }

    fn build_empty_bins(chr_info: &[ChrInfo]) -> Vec<SnpLocusBin> {
        let total_bins = chr_info
            .last()
            .map(|chr| chr.bin_offset + chr.bin_count)
            .unwrap_or(0);

        vec![SnpLocusBin::default(); total_bins]
    }

    pub fn match_read<'a>(
        &'a self,
        read: &'a AlignedRead,
        min_baseq: u8,
    ) -> SnpReadMatch<'a> {
        let mut out = SnpReadMatch::default();

        for hit in self.observed_snps(read, min_baseq) {
            out.push(hit);
        }

        out
    }

    fn validate_and_canonicalize_raw_loci(
        chr_info: &[ChrInfo],
        loci: &mut [RawSnpRecord],
    ) -> Result<()> {
        loci.sort_by(|a, b| {
            (a.chr_id, a.pos0, a.name.as_str(), a.vcf_id.as_str())
                .cmp(&(b.chr_id, b.pos0, b.name.as_str(), b.vcf_id.as_str()))
        });

        for locus in loci.iter_mut() {
            let Some(chr) = chr_info.get(locus.chr_id) else {
                return Err(anyhow!(
                    "SNP '{}' has invalid chr_id {}",
                    locus.name,
                    locus.chr_id
                ));
            };

            if locus.pos0 >= chr.length {
                return Err(anyhow!(
                    "SNP '{}' position {} is outside chromosome '{}' length {}",
                    locus.name,
                    locus.pos0,
                    chr.name,
                    chr.length
                ));
            }

            if locus.alternates.is_empty() {
                return Err(anyhow!("SNP '{}' has no alternate alleles", locus.name));
            }

            locus.reference = locus.reference.to_ascii_uppercase();

            for alt in &mut locus.alternates {
                *alt = alt.to_ascii_uppercase();
            }

            if locus.name.is_empty() {
                locus.name = Self::make_default_locus_name(chr, locus);
            }

        }

        Ok(())
    }

    fn make_default_locus_name(chr: &ChrInfo, locus: &RawSnpRecord) -> String {
        let pos1 = locus.pos0 + 1;
        let reference = locus.reference as char;

        let alt = locus
            .alternates
            .iter()
            .map(|base| (*base as char).to_string())
            .collect::<Vec<_>>()
            .join(",");

        format!("{}:{}:{}/{}", chr.name, pos1, reference, alt)
    }

    fn fill_bins(
        chr_info: &[ChrInfo],
        bin_width: u32,
        loci: &[SnpLocus],
        bins: &mut [SnpLocusBin],
    ) -> Result<()> {
        for locus in loci {
            let chr = &chr_info[locus.chr_id];
            let local_bin = (locus.pos0 / bin_width) as usize;
            let global_bin = chr.bin_offset + local_bin;

            let Some(bin) = bins.get_mut(global_bin) else {
                return Err(anyhow!(
                    "internal error: global bin {} does not exist for SNP '{}'",
                    global_bin,
                    locus.name
                ));
            };

            bin.locus_ids.push(locus.id);
        }

        Ok(())
    }

    fn sort_bins(loci: &[SnpLocus], bins: &mut [SnpLocusBin]) {
        for bin in bins {
            bin.locus_ids
                .sort_by_key(|&id| (loci[id].pos0, loci[id].id));
        }
    }

    fn build_name_to_id(loci: &[SnpLocus]) -> HashMap<String, u64> {
        let mut out = HashMap::with_capacity(loci.len() * 2);

        for locus in loci {
            let id = locus.id as u64;

            out.entry(locus.name.clone()).or_insert(id);

            if !locus.vcf_id.is_empty() && locus.vcf_id != "." {
                out.entry(locus.vcf_id.clone()).or_insert(id);
            }
        }

        out
    }
}

impl FeatureIndex for SnpIndex {
    fn feature_name(&self, feature_id: u64) -> &str {
        &self.loci[feature_id as usize].name
    }

    fn feature_id(&self, name: &str) -> Option<u64> {
        self.name_to_id.get(name).copied()
    }

    fn to_10x_feature_line(&self, feature_id: u64) -> String {
        let locus = &self.loci[feature_id as usize];
        let chr = &self.chr_info[locus.chr_id];

        let pos1 = locus.pos0 + 1;
        let reference = locus.reference as char;

        let alt = locus
            .alternates
            .iter()
            .map(|base| (*base as char).to_string())
            .collect::<Vec<_>>()
            .join(",");

        let feature_name = format!("{}:{}:{}/{}", chr.name, pos1, reference, alt);

        format!("{}\t{}\tSNP", locus.name, feature_name)
    }

    fn ordered_feature_ids(&self) -> Vec<u64> {
        self.loci.iter().map(|locus| locus.id as u64).collect()
    }
}

pub struct SnpPosIter<'a> {
    index: &'a SnpIndex,
    ids: &'a [usize],
    offset: usize,
}

impl<'a> SnpPosIter<'a> {
    fn empty(index: &'a SnpIndex) -> Self {
        Self {
            index,
            ids: &[],
            offset: 0,
        }
    }
}

impl<'a> Iterator for SnpPosIter<'a> {
    type Item = &'a SnpLocus;

    fn next(&mut self) -> Option<Self::Item> {
        let id = *self.ids.get(self.offset)?;
        self.offset += 1;
        Some(&self.index.loci[id])
    }
}

pub struct SnpRangeIter<'a> {
    index: &'a SnpIndex,
    start0: u32,
    end0: u32,
    current_global_bin: usize,
    last_global_bin: usize,
    current_ids: &'a [usize],
    current_offset: usize,
}

impl<'a> SnpRangeIter<'a> {
    fn empty(index: &'a SnpIndex) -> Self {
        Self {
            index,
            start0: 0,
            end0: 0,
            current_global_bin: 1,
            last_global_bin: 0,
            current_ids: &[],
            current_offset: 0,
        }
    }
}

impl<'a> Iterator for SnpRangeIter<'a> {
    type Item = &'a SnpLocus;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            while let Some(&id) = self.current_ids.get(self.current_offset) {
                self.current_offset += 1;

                let locus = &self.index.loci[id];

                if locus.pos0 >= self.start0 && locus.pos0 < self.end0 {
                    return Some(locus);
                }
            }

            if self.current_global_bin > self.last_global_bin {
                return None;
            }

            let bin = &self.index.bins[self.current_global_bin];
            self.current_ids = bin.locus_ids.as_slice();
            self.current_offset = 0;
            self.current_global_bin += 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::read::{AlignedRead, ReadOpKind, Strand};
    use crate::Genome;


    /*fn locus(
        chr_id: usize,
        pos0: u32,
        reference: u8,
        alternates: &[u8],
        name: &str,
        vcf_id: &str,
    ) -> RawSnpRecord {
        RawSnpRecord {
            chr_id,
            pos0,
            reference,
            alternates: alternates.to_vec(),
            name: name.to_string(),
            vcf_id: vcf_id.to_string(),
        }
    }*/

    fn test_index() -> SnpIndex {
        SnpIndex::new(
            &["chrA".to_string(), "chrB".to_string(), "MT".to_string()],
            &[1_000, 500, 100],
            vec![
                RawSnpRecord::new( 0, 125, b'a', b"T", "rs_geneA_125", "rsA"),
                RawSnpRecord::new( 0, 445, b'C', b"G", "rs_geneB_445", "rsB"),
                RawSnpRecord::new( 0, 545, b'G', b"T", "rs_geneC_545", "rsC"),
                RawSnpRecord::new( 0, 835, b'C', b"A", "rs_geneD_835", "rsD"),
                RawSnpRecord::new( 1, 25, b'T', b"C", "rs_chrB_25", "rsE"),
                RawSnpRecord::new( 2, 5, b'A', b"G", "rs_mt_5", "rsMT"),
            ],
            100,
        )
        .unwrap()
    }

    #[test]
    fn builds_chr_info_and_bins() {
        let index = test_index();

        assert_eq!(index.chr_info[0].name, "chrA");
        assert_eq!(index.chr_info[0].bin_offset, 0);
        assert_eq!(index.chr_info[0].bin_count, 10);

        assert_eq!(index.chr_info[1].name, "chrB");
        assert_eq!(index.chr_info[1].bin_offset, 10);
        assert_eq!(index.chr_info[1].bin_count, 5);

        assert_eq!(index.chr_info[2].name, "MT");
        assert_eq!(index.chr_info[2].bin_offset, 15);
        assert_eq!(index.chr_info[2].bin_count, 1);

        assert_eq!(index.bins.len(), 16);
    }

    #[test]
    fn fuzzy_chr_lookup_works() {
        let index = test_index();

        assert_eq!(index.chr_id("chrA"), Some(0));
        assert_eq!(index.chr_id("A"), Some(0));

        assert_eq!(index.chr_id("chrB"), Some(1));
        assert_eq!(index.chr_id("B"), Some(1));

        assert_eq!(index.chr_id("MT"), Some(2));
        assert_eq!(index.chr_id("chrM"), Some(2));
        assert_eq!(index.chr_id("M"), Some(2));

        assert_eq!(index.chr_id("does_not_exist"), None);
    }

    #[test]
    fn global_bin_for_pos_works() {
        let index = test_index();

        assert_eq!(index.global_bin_for_pos(0, 0), Some(0));
        assert_eq!(index.global_bin_for_pos(0, 99), Some(0));
        assert_eq!(index.global_bin_for_pos(0, 100), Some(1));
        assert_eq!(index.global_bin_for_pos(0, 999), Some(9));
        assert_eq!(index.global_bin_for_pos(0, 1_000), None);

        assert_eq!(index.global_bin_for_pos(1, 25), Some(10));
        assert_eq!(index.global_bin_for_pos(2, 5), Some(15));
    }

    #[test]
    fn loci_are_canonicalized_and_sorted() {
        let index = test_index();

        for (expected_id, locus) in index.loci.iter().enumerate() {
            assert_eq!(locus.id, expected_id);
        }

        let names: Vec<_> = index.loci.iter().map(|l| l.name.as_str()).collect();

        assert_eq!(
            names,
            vec![
                "rs_geneA_125",
                "rs_geneB_445",
                "rs_geneC_545",
                "rs_geneD_835",
                "rs_chrB_25",
                "rs_mt_5",
            ]
        );

        assert_eq!(index.loci[0].reference, b'A');
        assert_eq!(index.loci[0].alternates, vec![b'T']);
    }

    #[test]
    fn snps_at_pos_finds_exact_hits() {
        let index = test_index();

        let hits: Vec<_> = index.snps_at_pos(0, 125).map(|s| s.name.as_str()).collect();
        assert_eq!(hits, vec!["rs_geneA_125"]);

        let hits: Vec<_> = index.snps_at_pos(0, 126).collect();
        assert!(hits.is_empty());

        let hits: Vec<_> = index
            .snps_at_chr_name_pos("A", 445)
            .map(|s| s.name.as_str())
            .collect();

        assert_eq!(hits, vec!["rs_geneB_445"]);
    }

    #[test]
    fn snps_at_pos_allows_multiple_snps_at_same_position() {
        let index = SnpIndex::new(
            &["chrA".to_string()],
            &[1_000],
            vec![
                RawSnpRecord::new(0, 125, b'A', b"T", "snp1", "vcf1"),
                RawSnpRecord::new(0, 125, b'A', b"G", "snp2", "vcf2"),
            ],
            100,
        )
        .unwrap();

        let hits: Vec<_> = index.snps_at_pos(0, 125).map(|s| s.name.as_str()).collect();

        assert_eq!(hits, vec!["snp1", "snp2"]);
    }

    #[test]
    fn snps_in_range_crosses_bins() {
        let index = test_index();

        let hits: Vec<_> = index
            .snps_in_range(0, 100, 600)
            .map(|s| s.name.as_str())
            .collect();

        assert_eq!(hits, vec!["rs_geneA_125", "rs_geneB_445", "rs_geneC_545"]);
    }

    #[test]
    fn feature_index_works() {
        let index = test_index();

        assert_eq!(index.feature_id("rs_geneA_125"), Some(0));
        assert_eq!(index.feature_id("rsA"), Some(0));
        assert_eq!(index.feature_name(0), "rs_geneA_125");

        assert_eq!(
            index.to_10x_feature_line(0),
            "rs_geneA_125\tchrA:126:A/T\tSNP"
        );

        assert_eq!(index.ordered_feature_ids(), vec![0, 1, 2, 3, 4, 5]);
    }

    #[test]
    fn observed_snps_reports_ref_alt_and_other_cleanly() {
        let index = SnpIndex::new(
            &["chrA".to_string()],
            &[1_000],
            vec![
                RawSnpRecord::new(0, 100, b'A', b"T", "ref_hit", "rs1"),
                RawSnpRecord::new(0, 101, b'C', b"T", "alt_hit", "rs2"),
                RawSnpRecord::new(0, 102, b'G', b"A", "other_hit", "rs3"),
            ],
            100,
        )
        .unwrap();

        let read = AlignedRead::new(
            0,
            Strand::Plus,
            100,
            b"ATC".to_vec(),
            Some(vec![30, 30, 30]),
            vec![(ReadOpKind::Match, 3)],
        );

        let hits = index.observed_snps(&read, 20);

        assert_eq!(hits.len(), 3);

        assert_eq!(hits[0].locus.name, "ref_hit");
        assert!(hits[0].is_ref());
        assert!(!hits[0].is_alt());

        assert_eq!(hits[1].locus.name, "alt_hit");
        assert!(!hits[1].is_ref());
        assert!(hits[1].is_alt());

        assert_eq!(hits[2].locus.name, "other_hit");
        assert!(hits[2].is_other());
    }

    #[test]
    fn observed_snps_respects_base_quality() {
        let index = SnpIndex::new(
            &["chrA".to_string()],
            &[1_000],
            vec![RawSnpRecord::new(0, 100, b'A', b"T", "lowq", "rs1")],
            100,
        )
        .unwrap();

        let read = AlignedRead::new(
            0,
            Strand::Plus,
            100,
            b"A".to_vec(),
            Some(vec![10]),
            vec![(ReadOpKind::Match, 1)],
        );

        assert!(index.observed_snps(&read, 20).is_empty());
        assert_eq!(index.observed_snps(&read, 10).len(), 1);
    }

    #[test]
    fn constructor_rejects_bad_inputs() {
        assert!(SnpIndex::new(&[], &[], vec![], 100).is_err());

        assert!(SnpIndex::new(
            &["chrA".to_string()],
            &[1_000],
            vec![],
            0,
        )
        .is_err());

        assert!(SnpIndex::new(
            &["chrA".to_string()],
            &[1_000],
            vec![RawSnpRecord::new( 0, 1_000, b'A', b"T", "bad", "rsBad")],
            100,
        )
        .is_err());
    }

    #[test]
    fn real_world_stress_test() {
        let genome= Genome::from_fasta("tests/data/chr17_slice.fa.gz").unwrap();
        let mut read = AlignedRead::new(
            0,                  // e.g. 0 for chr17 slice
            Strand::Minus,         // flag 16
            1,               // POS-1
            b"GTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTAGTATTACTGAATTTATTTGTATTCACTTTTTATGTTCAAGACAAGAAAACTCAACTGAAAGAGGAGTATATTAAACCTTTCTTTTCTCCCATACCATGGCCCATTCATTTATACTTTCAATTTGCCCGACTTCCTCCATTAGGGCTCATCCCAAATTCCATTCCACTTCCTCCATTCCCCCTTTTACATTCACAACCCTTGTGCAACATTCCTTTTCTCTCAAGTCAGACCAGACCACATGTTTCAGTGCACATGGATCCCATGTACTCCTGTCTCTTCATACCACTGCTTACTACACGACGCTCTTCCGACCATAGTGAGCAGGTGACGGCACAGGCCATTTTTTTTTTTTTTTTTTTTTTTTTTTCCCAGCTAATTTTGTATTTTTAGTAGAGATGGGGTTTCATCTTGTTGGTCAGGCTAGTCTCCAACTCCTGACCTCAAGTGATCCGCCGGCCTTGGCCTCCCAGAGTGTCGGGATTACAGGCGTGAGCCACCGCGCCCAGCTGCAATTCTCAATCCAGCAATCTATCAGCTCACAAAGAAGGGTTTGGCAGAGTCAGAATGGATCTGAAGCAGCACAGCAGAAATAAAAACAGAGCTAGGGACATAATCTGAGGAAAAGGTTGGGATTGAAAGGGTGGAGGATGGTTGCGTAAAAGGGAGATTTTCTGTTTACTTGAGATCTGCAGAGTAGGGTCCATGTGTCCTCGTCTGTAGGTGAACACACGTGCACACACAAATGGTAACTGTGTGTGATGATGGACGTGCTGATTAATTTGATTGTGAATGGTTTCACAATATACCCCCATATCAAATCATCACACTGCATACCTTGAATGTATGCAATTTTGGGGGATGAATTATGCCCAACAAAGGAAATAAAAAAAATTAAAAAAAAAAAAAAGGAGTCATCAAGGTGACCTGTCTACAGAGGCAAGGACAGGGACTGAGCTTCAGGAGCTCTAGTTTGCCTGTTGGGTAGGGGACAGATGTTTAAGTTAAAAGTCTCTGAAAGAGGTGAAATTCTGGATCTCCTGGGAAGAGTGTTTGGCATTCCCTGGAGTGAGCCCTGCTCCCCCCTGGGCCCTTCTGGAGCCTGGGCATCCTTGAGTTCCAAGGCTGATTCAGCTCTCTGGAACATCTCGAAGCGCTCACGCCCACGGATCTGAAGGGTGAGAGATTCTCCATCCAGCGGTTTCTTCTTTGGCTGGGGAGAGGAGCTGGTGTTGTTGGGCAGTGCTGCTTAGCTCCTGGGGCAGCTCTGTGGTGAGGCTCCCCTTTCTTGCGGACATCTTCTTCCCTCTGTGTGCCTCGGGCCCCTCCCAGGACAGGCACATGCACCTCAAAGCTGTTCCGTCCCAGTACATGACCACTGGAGTCTTCCAGTGTGATGATGGTGAGGATGGGCCTCCGGTTCATGCCGCCCATGCAGGAACTGTTACACATGTAGTTGTAGTGGATGGTGACAGAGTCAGAGCCAACCTCAGGCGGCTCATAGGGCACCACCACACGATGTCGAAGAGTGTTTCTGTCATCCAAATACTCCACACGGAAATTTCCTTCCACGTGCAGAGCCTAGGAGGGCGAGGACGGCCATCTACAGTGCTCACACACGGCGGGGCAGCGCCTCACAGCTCCAGGCCGTGACTGCTTGGGATGGCCATGGCGCGGACGGGCGGGCGCCGGCGGGCGTGGAATCAACCCACAGCTGCACAGGGCAGGTCTTGGCCAGTTGGCAAACATCTTGTTGAGGGCAGGGCAGTATGTGGACATACTTGGCTGTCCCAGAATGCAAGAGAAGCCACGGGAAACCGGAGCCGCCCTGGTAGGTTTTCTGGGAGGGATAGAAGATGACAGGGGCCAGGAGGGGGCTGGTGCAGGGGCCGCCGGTGTAGGAGCTGCTGGTGCAGGGGCCACGGGGAGAGTCCCTCCCGCATTTGTAGTTTCATCCTGACCTGGGTCTTCACCC".to_vec(),
            Some(b"=>?BBC;74444555554.*.*(((('%%&),8?AECAA@@@EEGGDC97789@@ABFGGE>?=FFHB<:770++-)&(())**+222.-0227=?A999::?>>@@AC@BA@?@<55-)&&&&&&&,00@@@@>A@A65432344463222122266679KHIIGHD<;;:;FDCCCDCD11,,2,+**+*++)'''()+)((()67=;=85-)))*DAAABBHJIEDACDABA??A@ACAA=4;744321+++++++211110/.**,,,,?A@7777)((''((((,+**('''('%%%%%&%%%%%&(>?C>><<;<==;;::553334,,,,,++('&'(*/011159:;?@A@BBCDIHJG<;;;:<<<=7777899887665444320/246;>??GGHHEN@@>@NKSQSSKF554548889KHRPMFCCBB=:99;;BIKKSCBBBDSFKFEDGNLJLIJKSIQHDJFSNKOGI@@@@B>>>=<@@DHF22223;(((()>DCHSJDCABCCGOKSEMHEEEEEDB7732211101.)))+----.66777766?>>EGCA=700001FGIC3212122212222;;===;99:00/0/++++*.---+++1444971234?=>;;8776232222''&&':;.7<//73223110'&&')*13.-00+++1130-++,++34211001277<?>?A85543200*('&))''''(,-)(342(''(-/.44442///14B=;65445557677<532322222679:CGIKNSJNHMEOKJILQSMIJICCCBBC876665....>??:::<>AB;84346654477???>BEIFG@?=AGFHJNIDDEEEHHHJEDCCBBGJDC>>:669:<:78:51,++*),,++-.,+(()))---.,21,*++.3561034661-++*)(&'''$$%&')-...79EEADE=;::656778,,,++7884433577:445668<:1198@AKIEKJRJSKSKNJSG9B=<;*311//0004332122.-*'))&''(+,,-,,,-,/*)'''&'))(+--.+,,0*).**/45;;<??55667?-,8))).?BCEGJAB?<=<;;;;)'''&-//-.+(&&&,++*/0,8:;:;.---34/))))*,*((***+*++,,,)&&&&))//66111.-33444446<=022221,+''(()6./-)(&&&&&&'),----1----.CBE8789;9<=;>??36<::<;=@?@>86888>((('')/./41--,*)))*++)(01-../013***(*--./013334HEA:39774(&&%%%%%&')*0+'&)*(((+0/*(()(&&&&&'((53550198><<,+'%&'(/./.---0.0001::::;B><;:100+*+++<<<=>@BCEGIJGIDFHHHAA@@@@.-,**+0020*,.459==;00../569:43333>>ACEGGGHEISRJIDDAAAGDJGFDDB>977776)(((*)))003478;;<=?BABB975552211155?A@<=>?@DB<,,,,,/..-'&%&&)*)***)'')))))*&44554111151(&&%&&&&&&&,/,'$$#&&&&'(*)*)())''&$##$$##$$&%$$####'&)('(&((,*(('(**+2354343,)('&&('''%%&&%%%%%''(''%&%%%((()*.2210.(((&%('&$&(''((''(()**&&&&&+')')***,....45===?@?A97:87101122678999>>669+)/017:<=?><;;>,9<96)*('&%%$%%,&###$(&'',.89;>:8***+++*)((,,-./-+***,)./032...//////154++,,,----,,-,(*++,--'''''97788=AJLLKECE@@AA@?<43322--++*)+)))))))*2;;<64213++++.237---:0:+*437,0(&%%&%$$%%%%%%%)&%%$%%&'&'%%%&%%%$%%(**8::;?CB?>555".to_vec()),
            vec![
                (SoftClip, 400),
                (Match, 249),
                (Ins, 1),
                (Match, 79),
                (Ins, 1),
                (Match, 18),
                (Ins, 1),
                (Match, 17),
                (Del, 2),
                (Match, 55),
                (Del, 2),
                (Match, 80),
                (Del, 1),
                (Match, 10),
                (Del, 3),
                (Match, 5),
                (Del, 1),
                (Match, 28),
                (Del, 1),
                (Match, 73),
                (Ins, 1),
                (Match, 38),
                (Ins, 2),
                (Match, 35),
                (RefSkip, 3183),
                (Match, 33),
                (Ins, 1),
                (Match, 28),
                (Del, 1),
                (Match, 14),
                (Ins, 1),
                (Match, 31),
                (RefSkip, 2819),
                (Match, 74),
                (RefSkip, 92),
                (Match, 2),
                (Del, 1),
                (Match, 5),
                (Del, 2),
                (Match, 3),
                (Del, 1),
                (Match, 3),
                (Del, 1),
                (Match, 10),
                (Ins, 1),
                (Match, 35),
                (Ins, 1),
                (Match, 12),
                (Ins, 3),
                (Match, 19),
                (Del, 4),
                (Match, 39),
                (RefSkip, 343),
                (Match, 110),
                (RefSkip, 568),
                (Match, 84),
                (Del, 3),
                (Match, 7),
                (Del, 1),
                (Match, 8),
                (Del, 1),
                (Match, 6),
                (Ins, 1),
                (Match, 3),
                (RefSkip, 81),
                (Del, 3),
                (Match, 8),
                (Del, 1),
                (Match, 1),
                (Del, 1),
                (Match, 7),
                (Ins, 4),
                (Match, 5),
                (Del, 1),
                (Match, 16),
                (Del, 1),
                (Match, 5),
                (Del, 3),
                (Match, 2),
                (Del, 1),
                (Match, 1),
                (Del, 1),
                (Match, 14),
                (Del, 1),
                (Match, 19),
                (Ins, 2),
                (Match, 9),
                (Del, 1),
                (Match, 3),
                (Del, 2),
                (Match, 50),
                (Del, 1),
                (Match, 27),
                (RefSkip, 757),
                (Match, 4),
                (Del, 3),
                (Match, 1),
                (Del, 2),
                (Match, 25),
                (Ins, 2),
                (Match, 7),
                (Del, 3),
                (Match, 2),
                (Ins, 1),
                (Match, 33),
                (Del, 1),
                (Match, 86),
                (Ins, 1),
                (Match, 1),
                (Del, 2),
                (Match, 11),
                (Del, 1),
                (Match, 1),
                (Del, 1),
                (Match, 27),
                (SoftClip, 3),
            ],
        );
        let slice_start_1based = 7_666_729u32;
        let snps = vec![
            RawSnpRecord::new(0, 7_675_994 - slice_start_1based, b'G', b"C", "TP53_c375G>C", "TP53_c375G>C"),
            RawSnpRecord::new(0, 7_674_894 - slice_start_1based, b'C', b"T", "TP53_c637C>T","TP53_c637C>T"),
            RawSnpRecord::new(0, 7_674_953 - slice_start_1based, b'A', b"T", "TP53_c578A>T","TP53_c637C>T"),
        ];

        let pos0 =  7_674_953 - slice_start_1based;
        let obs = read.base_at_ref_pos(pos0).unwrap();
        assert_eq!(obs.base, b'A', "base at {}='{}' could be 'A' (in the read))",
            pos0,
            obs.base as char
        );

        read.seq[obs.read_pos as usize] = b'T';

        let obs = read.base_at_ref_pos(pos0).unwrap();
        assert_eq!(obs.base, b'T', "base at {}='{}' not changed to 'T' (in the read))",
            pos0,
            obs.base as char
        );

        let index = SnpIndex::new(
            &genome.chr_names(),
            &genome.chr_lengths(),
            snps,
            1_000_000,
        );

        assert!(
            index.is_ok(),
            "SNP index should build for chr17 slice-local SNP coordinates: {:?}",
            index.err()
        );

        let index = index.unwrap();

        assert!(read.base_at_ref_pos(7_675_994 - slice_start_1based).is_none());

        assert!(read
            .base_at_ref_pos(7_674_894 - slice_start_1based)
            .is_some());

        assert!(read
            .base_at_ref_pos(7_674_953 - slice_start_1based)
            .is_some());

        let (ref_ids, alt_ids, other_ids) =
            index.get_ref_alt_other_ids_for_read(&read, 0);

        assert_eq!(ref_ids, vec![0u64]);
        assert_eq!(alt_ids, vec![1u64]);
        assert!(other_ids.is_empty());
    }
}


impl fmt::Display for SnpIndex {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let total_bins = self.bins.len();
        let non_empty_bins = self
            .bins
            .iter()
            .filter(|bin| !bin.locus_ids.is_empty())
            .count();

        writeln!(
            f,
            "SnpIndex: {} SNPs in {} chromosomes; bin_width={}; bins={}/{} non-empty",
            self.loci.len(),
            self.chr_info.len(),
            self.bin_width,
            non_empty_bins,
            total_bins
        )?;

        if let Some(first) = self.loci.first() {
            let chr = self.chr_name(first.chr_id).unwrap_or("?");
            let alt = first
                .alternates
                .iter()
                .map(|b| (*b as char).to_string())
                .collect::<Vec<_>>()
                .join(",");

            writeln!(
                f,
                "first SNP: {}:{} {}>{} ({})",
                chr,
                first.pos0 + 1,
                first.reference as char,
                alt,
                first.name
            )?;
        }

        Ok(())
    }
}