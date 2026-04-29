//! SNP index construction.
//!
//! This module owns the final `SnpIndex` structure.
//! It stores SNP loci sorted by chromosome and position, plus a flat CSR-style
//! bin index for fast genomic lookup.

use crate::locus::SnpLocus;
use crate::vcf::{SnpVcfReader, VcfReadOptions};
use anyhow::{Result, anyhow};
use scdata::feature_index::FeatureIndex;
use std::collections::HashMap;
use std::path::Path;

use std::fmt;

use crate::AlignedRead;

/// Default genomic bin width.
///
/// This is intentionally fairly coarse. SNPs are point features, so the
/// optimal value is usually determined by read length and SNP density.
pub const DEFAULT_BIN_WIDTH: u32 = 16_384;

/// Result of matching a read against SNP loci.
///
/// For each overlapped SNP:
/// - matching REF base → added to `ref_ids`
/// - matching ALT base → added to `alt_ids`
///
/// Notes:
/// - Each SNP appears in at most one list.
/// - SNPs not covered or not matching REF/ALT are ignored.
/// - Stores SNP *ids* (not positions).
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct SnpReadMatch {
    pub ref_ids: Vec<u32>,
    pub alt_ids: Vec<u32>,
}

/// A flat, chromosome-aware SNP index.
///
/// Bins are flattened into a single global array:
///
/// `global_bin = chr_bin_offsets[chr_id] + local_bin`
///
/// `bin_starts` uses a CSR-style layout:
///
/// SNPs in bin `b` are stored in:
///
/// `loci[bin_starts[b]..bin_starts[b + 1]]`
#[derive(Debug, Clone)]
pub struct SnpIndex {
    pub chr_names: Vec<String>,
    pub chr_lengths: Vec<u32>,
    pub chr_bin_offsets: Vec<usize>,
    pub chr_bin_counts: Vec<usize>,
    pub bin_width: u32,

    /// Loci sorted by `(chr_id, pos0, name)`, then reordered by bin.
    ///
    /// Within each bin, loci remain sorted by `(chr_id, pos0, name)`.
    pub loci: Vec<SnpLocus>,

    /// CSR-style bin starts.
    ///
    /// Length is `number_of_global_bins + 1`.
    pub bin_starts: Vec<usize>,

    /// Feature name to feature id lookup.
    pub name_to_id: HashMap<String, u64>,
}

impl FeatureIndex for SnpIndex {
    /// id to name translation
    fn feature_name(&self, feature_id: u64) -> &str {
        &self.loci[feature_id as usize].name
    }

    /// name to id translation
    fn feature_id(&self, name: &str) -> Option<u64> {
        self.name_to_id.get(name).copied()
    }

    /// One line in features.tsv.
    ///
    /// 10x convention is:
    /// feature_id<TAB>feature_name<TAB>feature_type
    fn to_10x_feature_line(&self, feature_id: u64) -> String {
        let locus = &self.loci[feature_id as usize];

        format!("{}\t{}\tSNP", locus.vcf_id, locus.name)
    }

    /// Export SNPs in feature-id order.
    fn ordered_feature_ids(&self) -> Vec<u64> {
        (0..self.loci.len() as u64).collect()
    }
}

impl SnpIndex {
    /// Build an SNP index from already parsed loci.
    ///
    /// `chr_names` and `chr_lengths` should normally come from the BAM header.
    /// Their order defines the chromosome id space used by this index.
    pub fn new(
        chr_names: Vec<String>,
        chr_lengths: Vec<u32>,
        loci: Vec<SnpLocus>,
        bin_width: u32,
    ) -> Result<Self> {
        Self::validate_genome(&chr_names, &chr_lengths)?;
        Self::validate_bin_width(bin_width)?;

        let chr_bin_counts = Self::build_chr_bin_counts(&chr_lengths, bin_width);
        let chr_bin_offsets = Self::build_chr_bin_offsets(&chr_bin_counts);
        let total_bins = Self::total_bins(&chr_bin_counts);

        let mut index = Self {
            chr_names,
            chr_lengths,
            chr_bin_offsets,
            chr_bin_counts,
            bin_width,
            loci: Vec::new(),
            bin_starts: vec![0; total_bins + 1],
            name_to_id: HashMap::new(),
        };

        index.rebuild_from_loci(loci)?;

        Ok(index)
    }

    /// Match an aligned read against indexed SNP loci.
    ///
    /// For every SNP covered by the read, the observed read base is compared
    /// against the locus reference and alternate alleles.
    pub fn match_read(&self, read: &AlignedRead, min_baseq: u8) -> SnpReadMatch {
        let Some((start0, end0)) = read.ref_span() else {
            return SnpReadMatch::default();
        };

        let Some(bin_range) = self.global_bins_for_span(read.chr_id, start0, end0) else {
            return SnpReadMatch::default();
        };

        let mut hits = SnpReadMatch::default();

        for global_bin in bin_range {
            let Some(loci) = self.loci_in_global_bin(global_bin) else {
                continue;
            };

            for locus in loci {
                if locus.chr_id != read.chr_id {
                    continue;
                }

                if locus.pos0 < start0 || locus.pos0 >= end0 {
                    continue;
                }

                let Some(obs) = read.base_at_ref_pos(locus.pos0) else {
                    continue;
                };

                if let Some(q) = obs.qual
                    && q < min_baseq
                {
                    continue;
                }

                if locus.is_reference_base(obs.base) {
                    hits.ref_ids.push(locus.id as u32);
                } else if locus.is_alternate_base(obs.base) {
                    hits.alt_ids.push(locus.id as u32);
                }
            }
        }

        hits
    }

    /// Build an SNP index directly from a VCF/BCF file.
    ///
    /// The `chr_names` and `chr_lengths` should come from the BAM header.
    /// VCF contigs are mapped into this chromosome id space by name.
    pub fn from_vcf_path<P: AsRef<Path>>(
        path: P,
        chr_names: Vec<String>,
        chr_lengths: Vec<u32>,
        bin_width: u32,
        options: &VcfReadOptions,
    ) -> Result<Self> {
        let chr_map = Self::build_chr_map(&chr_names);
        let raw = SnpVcfReader::read_path(path, &chr_map, options)?;
        let loci = SnpLocus::from_raw_records(raw);

        Self::new(chr_names, chr_lengths, loci, bin_width)
    }

    /// Rebuild bins and name lookup from a locus vector.
    ///
    /// Loci are sorted and reassigned ids so feature ids are deterministic.
    pub fn rebuild_from_loci(&mut self, loci: Vec<SnpLocus>) -> Result<()> {
        let mut loci = Self::sort_and_reassign_loci(loci);
        Self::validate_loci_against_genome(&loci, &self.chr_lengths)?;

        let old_to_new =
            Self::build_loci_order_by_bin(&loci, &self.chr_bin_offsets, self.bin_width)?;

        loci = Self::reorder_loci_by_indices(loci, &old_to_new);
        Self::reassign_locus_ids_in_place(&mut loci);

        self.bin_starts = Self::build_bin_starts(
            &loci,
            self.chr_bin_offsets.len(),
            &self.chr_bin_offsets,
            &self.chr_bin_counts,
            self.bin_width,
        )?;

        self.name_to_id = Self::build_name_to_id(&loci)?;
        self.loci = loci;

        Ok(())
    }

    /// Return the number of SNP loci.
    pub fn len(&self) -> usize {
        self.loci.len()
    }

    /// Return true if the index contains no SNP loci.
    pub fn is_empty(&self) -> bool {
        self.loci.is_empty()
    }

    /// Return the number of global bins.
    pub fn n_bins(&self) -> usize {
        self.bin_starts.len().saturating_sub(1)
    }

    /// Return the global bin for `chr_id` and `pos0`.
    pub fn global_bin_for_pos(&self, chr_id: usize, pos0: u32) -> Option<usize> {
        if chr_id >= self.chr_lengths.len() {
            return None;
        }

        if pos0 >= self.chr_lengths[chr_id] {
            return None;
        }

        let local_bin = (pos0 / self.bin_width) as usize;

        if local_bin >= self.chr_bin_counts[chr_id] {
            return None;
        }

        Some(self.chr_bin_offsets[chr_id] + local_bin)
    }

    /// Return the loci slice belonging to a global bin.
    pub fn loci_in_global_bin(&self, global_bin: usize) -> Option<&[SnpLocus]> {
        if global_bin + 1 >= self.bin_starts.len() {
            return None;
        }

        let start = self.bin_starts[global_bin];
        let end = self.bin_starts[global_bin + 1];

        Some(&self.loci[start..end])
    }

    /// Return the loci slice for a chromosome-local bin.
    pub fn loci_in_chr_bin(&self, chr_id: usize, local_bin: usize) -> Option<&[SnpLocus]> {
        if chr_id >= self.chr_bin_offsets.len() {
            return None;
        }

        if local_bin >= self.chr_bin_counts[chr_id] {
            return None;
        }

        self.loci_in_global_bin(self.chr_bin_offsets[chr_id] + local_bin)
    }

    /// Return all global bins overlapped by a genomic half-open interval.
    ///
    /// The interval is `[start0, end0)`.
    pub fn global_bins_for_span(
        &self,
        chr_id: usize,
        start0: u32,
        end0: u32,
    ) -> Option<std::ops::RangeInclusive<usize>> {
        if chr_id >= self.chr_lengths.len() {
            return None;
        }

        if start0 >= end0 {
            return None;
        }

        let chr_len = self.chr_lengths[chr_id];
        if start0 >= chr_len {
            return None;
        }

        let clipped_end0 = end0.min(chr_len);
        let start_bin = (start0 / self.bin_width) as usize;
        let end_bin = ((clipped_end0 - 1) / self.bin_width) as usize;

        let offset = self.chr_bin_offsets[chr_id];

        Some((offset + start_bin)..=(offset + end_bin))
    }

    /// Build a chromosome name to chromosome id lookup.
    pub fn build_chr_map(chr_names: &[String]) -> HashMap<String, usize> {
        let mut map = HashMap::with_capacity(chr_names.len() * 3);

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

    /// Check that chromosome names and lengths describe the same genome.
    pub fn validate_genome(chr_names: &[String], chr_lengths: &[u32]) -> Result<()> {
        if chr_names.len() != chr_lengths.len() {
            return Err(anyhow!(
                "chromosome name count ({}) does not match chromosome length count ({})",
                chr_names.len(),
                chr_lengths.len()
            ));
        }

        if chr_names.is_empty() {
            return Err(anyhow!("genome must contain at least one chromosome"));
        }

        for (idx, name) in chr_names.iter().enumerate() {
            if name.is_empty() {
                return Err(anyhow!("chromosome name at index {idx} is empty"));
            }
        }

        for (idx, length) in chr_lengths.iter().enumerate() {
            if *length == 0 {
                return Err(anyhow!("chromosome length at index {idx} is zero"));
            }
        }

        Ok(())
    }

    /// Check that bin width is valid.
    pub fn validate_bin_width(bin_width: u32) -> Result<()> {
        if bin_width == 0 {
            Err(anyhow!("bin_width must be greater than zero"))
        } else {
            Ok(())
        }
    }

    /// Compute the number of bins for each chromosome.
    pub fn build_chr_bin_counts(chr_lengths: &[u32], bin_width: u32) -> Vec<usize> {
        chr_lengths
            .iter()
            .map(|len| Self::ceil_div_u32(*len, bin_width) as usize)
            .collect()
    }

    /// Compute chromosome-level global bin offsets from per-chromosome bin counts.
    pub fn build_chr_bin_offsets(chr_bin_counts: &[usize]) -> Vec<usize> {
        let mut offsets = Vec::with_capacity(chr_bin_counts.len());
        let mut running = 0usize;

        for count in chr_bin_counts {
            offsets.push(running);
            running += *count;
        }

        offsets
    }

    /// Return the total number of bins.
    pub fn total_bins(chr_bin_counts: &[usize]) -> usize {
        chr_bin_counts.iter().sum()
    }

    /// Sort loci by genomic order and reassign ids.
    pub fn sort_and_reassign_loci(mut loci: Vec<SnpLocus>) -> Vec<SnpLocus> {
        loci.sort_by(|a, b| {
            a.chr_id
                .cmp(&b.chr_id)
                .then_with(|| a.pos0.cmp(&b.pos0))
                .then_with(|| a.name.cmp(&b.name))
        });

        Self::reassign_locus_ids_in_place(&mut loci);

        loci
    }

    /// Reassign locus ids to match their current vector positions.
    pub fn reassign_locus_ids_in_place(loci: &mut [SnpLocus]) {
        for (id, locus) in loci.iter_mut().enumerate() {
            locus.id = id;
        }
    }

    /// Validate loci against chromosome ids and chromosome lengths.
    pub fn validate_loci_against_genome(loci: &[SnpLocus], chr_lengths: &[u32]) -> Result<()> {
        for locus in loci {
            if locus.chr_id >= chr_lengths.len() {
                return Err(anyhow!(
                    "locus {} has chr_id {}, but genome has only {} chromosomes",
                    locus.name,
                    locus.chr_id,
                    chr_lengths.len()
                ));
            }

            if locus.pos0 >= chr_lengths[locus.chr_id] {
                return Err(anyhow!(
                    "locus {} at chr_id {} pos0 {} is outside chromosome length {}",
                    locus.name,
                    locus.chr_id,
                    locus.pos0,
                    chr_lengths[locus.chr_id]
                ));
            }
        }

        Ok(())
    }

    /// Build locus vector order by global bin.
    ///
    /// Returns old locus indices ordered by the bin they belong to.
    pub fn build_loci_order_by_bin(
        loci: &[SnpLocus],
        chr_bin_offsets: &[usize],
        bin_width: u32,
    ) -> Result<Vec<usize>> {
        let mut order: Vec<usize> = (0..loci.len()).collect();

        order.sort_by_key(|old_idx| {
            let locus = &loci[*old_idx];
            chr_bin_offsets[locus.chr_id] + (locus.pos0 / bin_width) as usize
        });

        Ok(order)
    }

    /// Reorder loci according to an index vector.
    pub fn reorder_loci_by_indices(loci: Vec<SnpLocus>, order: &[usize]) -> Vec<SnpLocus> {
        order.iter().map(|old_idx| loci[*old_idx].clone()).collect()
    }

    /// Build CSR-style bin starts from loci already ordered by bin.
    pub fn build_bin_starts(
        loci: &[SnpLocus],
        n_chr: usize,
        chr_bin_offsets: &[usize],
        chr_bin_counts: &[usize],
        bin_width: u32,
    ) -> Result<Vec<usize>> {
        if chr_bin_offsets.len() != n_chr || chr_bin_counts.len() != n_chr {
            return Err(anyhow!("chromosome bin metadata length mismatch"));
        }

        let total_bins = Self::total_bins(chr_bin_counts);
        let mut bin_counts = vec![0usize; total_bins];

        for locus in loci {
            let global_bin = chr_bin_offsets[locus.chr_id] + (locus.pos0 / bin_width) as usize;

            if global_bin >= total_bins {
                return Err(anyhow!(
                    "locus {} maps to invalid global bin {}",
                    locus.name,
                    global_bin
                ));
            }

            bin_counts[global_bin] += 1;
        }

        let mut bin_starts = vec![0usize; total_bins + 1];

        for i in 0..total_bins {
            bin_starts[i + 1] = bin_starts[i] + bin_counts[i];
        }

        Ok(bin_starts)
    }

    /// Build a feature-name lookup table.
    pub fn build_name_to_id(loci: &[SnpLocus]) -> Result<HashMap<String, u64>> {
        let mut map = HashMap::with_capacity(loci.len());

        for locus in loci {
            let previous = map.insert(locus.name.clone(), locus.id as u64);

            if previous.is_some() {
                return Err(anyhow!("duplicate SNP feature name: {}", locus.name));
            }
        }

        Ok(map)
    }

    /// Integer ceiling division for positive denominator.
    pub fn ceil_div_u32(value: u32, denominator: u32) -> u32 {
        value.div_ceil(denominator)
    }
}




