//! VCF/BCF parsing for `snp-index`.
//!
//! This module reads variant files and converts accepted SNV records into a
//! small intermediate representation. It does not build bins and it does not
//! inspect BAM reads.

use anyhow::{Context, Result};
use rust_htslib::bcf;
use rust_htslib::bcf::Read as BcfRead;
use std::collections::HashMap;
use std::path::Path;

/// One raw SNP record parsed from a VCF/BCF file.
///
/// Coordinates are already converted to 0-based Rust/internal coordinates.
/// `chr_id` must refer to the chromosome id space provided by the caller,
/// typically derived from the BAM header.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RawSnpRecord {
    pub chr_id: usize,
    pub pos0: u32,
    pub reference: u8,
    pub alternates: Vec<u8>,
    pub name: String,
}

/// Options controlling which VCF/BCF records are accepted.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct VcfReadOptions {
    /// Keep only records with FILTER=PASS or no FILTER.
    pub pass_only: bool,

    /// Keep only simple SNVs: REF length 1 and all ALT alleles length 1.
    pub snps_only: bool,

    /// Keep only canonical A/C/G/T REF and ALT bases.
    pub acgt_only: bool,
}

impl Default for VcfReadOptions {
    /// Return strict defaults suitable for single-cell SNP counting.
    fn default() -> Self {
        Self {
            pass_only: true,
            snps_only: true,
            acgt_only: true,
        }
    }
}

/// Reader/converter for VCF/BCF SNP files.
pub struct SnpVcfReader;

impl SnpVcfReader {
    /// Read a VCF/BCF file and return accepted SNP records.
    ///
    /// `chr_map` translates VCF chromosome names into your internal chromosome
    /// ids. In your case this should normally come from the BAM header.
    ///
    /// Supported by htslib:
    /// - `.vcf`
    /// - `.vcf.gz` if bgzip-compressed
    /// - `.bcf`
    pub fn read_path<P: AsRef<Path>>(
        path: P,
        chr_map: &HashMap<String, usize>,
        options: &VcfReadOptions,
    ) -> Result<Vec<RawSnpRecord>> {
        let path_ref = path.as_ref();

        let mut reader = bcf::Reader::from_path(path_ref)
            .with_context(|| format!("failed to open VCF/BCF file {}", path_ref.display()))?;

        let header = reader.header().clone();
        let mut records = Vec::new();

        for rec_result in reader.records() {
            let record = rec_result.with_context(|| {
                format!("failed to read VCF/BCF record from {}", path_ref.display())
            })?;

            if let Some(raw) = Self::convert_record(&record, &header, chr_map, options)? {
                records.push(raw);
            }
        }

        Ok(records)
    }

    /// Convert one htslib BCF record into a raw SNP record.
    ///
    /// Returns `Ok(None)` when the record should simply be skipped.
    pub fn convert_record(
        record: &bcf::Record,
        header: &bcf::header::HeaderView,
        chr_map: &HashMap<String, usize>,
        options: &VcfReadOptions,
    ) -> Result<Option<RawSnpRecord>> {
        if options.pass_only && !Self::record_is_pass(record) {
            return Ok(None);
        }

        let rid = match record.rid() {
            Some(rid) => rid,
            None => return Ok(None),
        };

        let chr_name = std::str::from_utf8(header.rid2name(rid)?)
            .context("VCF/BCF chromosome name is not valid UTF-8")?;

        let chr_id = match chr_map.get(chr_name) {
            Some(chr_id) => *chr_id,
            None => return Ok(None),
        };

        let pos0 = match Self::record_pos0(record) {
            Some(pos0) => pos0,
            None => return Ok(None),
        };

        let alleles = record.alleles();
        let (reference, alternates) = match Self::parse_alleles(&alleles, options) {
            Some(parsed) => parsed,
            None => return Ok(None),
        };

        let name = Self::record_name(record, chr_name, pos0, reference, &alternates);

        Ok(Some(RawSnpRecord {
            chr_id,
            pos0,
            reference,
            alternates,
            name,
        }))
    }

    /// Return the 0-based record position, or `None` for invalid negative positions.
    pub fn record_pos0(record: &bcf::Record) -> Option<u32> {
        let pos = record.pos();

        if pos < 0 {
            None
        } else {
            Some(pos as u32)
        }
    }

    /// Parse REF and ALT alleles according to the supplied options.
    ///
    /// This currently focuses on SNVs. With default options:
    /// - REF must have length 1
    /// - every ALT must have length 1
    /// - all bases must be A/C/G/T
    pub fn parse_alleles(
        alleles: &[&[u8]],
        options: &VcfReadOptions,
    ) -> Option<(u8, Vec<u8>)> {
        if alleles.len() < 2 {
            return None;
        }

        let reference = alleles[0];

        if options.snps_only && reference.len() != 1 {
            return None;
        }

        if reference.is_empty() {
            return None;
        }

        let reference_base = reference[0].to_ascii_uppercase();

        if options.acgt_only && !Self::is_acgt(reference_base) {
            return None;
        }

        let mut alternates = Vec::with_capacity(alleles.len() - 1);

        for alt in &alleles[1..] {
            if alt.is_empty() {
                continue;
            }

            if options.snps_only && alt.len() != 1 {
                return None;
            }

            let alt_base = alt[0].to_ascii_uppercase();

            if options.acgt_only && !Self::is_acgt(alt_base) {
                return None;
            }

            alternates.push(alt_base);
        }

        if alternates.is_empty() {
            None
        } else {
            Some((reference_base, alternates))
        }
    }

    /// Check whether a VCF/BCF record passes filtering.
    ///
    /// Records with explicit `PASS` are accepted.
    /// Records with no FILTER values are also accepted.
    pub fn record_is_pass(record: &bcf::Record) -> bool {
        record.has_filter(&b"PASS"[..]) || record.filters().next().is_none()
    }

    /// Build a stable feature name for a SNP.
    ///
    /// If the VCF ID field is present and not `.`, that ID is used.
    /// Otherwise a coordinate-based name is generated using VCF-style 1-based
    /// coordinates:
    ///
    /// `chr17:7673803:C>T`
    pub fn record_name(
        record: &bcf::Record,
        chr_name: &str,
        pos0: u32,
        reference: u8,
        alternates: &[u8],
    ) -> String {
        let id = record.id();

        if !id.is_empty() && id != b"."
            && let Ok(id_string) = std::str::from_utf8(&id) {
                return id_string.to_string();
            }

        Self::fallback_name(chr_name, pos0, reference, alternates)
    }

    /// Build a coordinate/allele fallback name.
    pub fn fallback_name(
        chr_name: &str,
        pos0: u32,
        reference: u8,
        alternates: &[u8],
    ) -> String {
        let pos1 = pos0 + 1;
        let alt_string = alternates.iter().map(|b| *b as char).collect::<String>();

        format!("{chr_name}:{pos1}:{}>{alt_string}", reference as char)
    }

    /// Return true for canonical DNA bases.
    pub fn is_acgt(base: u8) -> bool {
        matches!(base, b'A' | b'C' | b'G' | b'T')
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_options_are_strict() {
        let options = VcfReadOptions::default();

        assert!(options.pass_only);
        assert!(options.snps_only);
        assert!(options.acgt_only);
    }

    #[test]
    fn is_acgt_accepts_only_canonical_uppercase_bases() {
        assert!(SnpVcfReader::is_acgt(b'A'));
        assert!(SnpVcfReader::is_acgt(b'C'));
        assert!(SnpVcfReader::is_acgt(b'G'));
        assert!(SnpVcfReader::is_acgt(b'T'));

        assert!(!SnpVcfReader::is_acgt(b'N'));
        assert!(!SnpVcfReader::is_acgt(b'a'));
        assert!(!SnpVcfReader::is_acgt(b'-'));
    }

    #[test]
    fn parse_alleles_accepts_simple_snv() {
        let options = VcfReadOptions::default();
        let alleles: Vec<&[u8]> = vec![b"C", b"T"];

        let parsed = SnpVcfReader::parse_alleles(&alleles, &options).unwrap();

        assert_eq!(parsed.0, b'C');
        assert_eq!(parsed.1, vec![b'T']);
    }

    #[test]
    fn parse_alleles_uppercases_bases() {
        let options = VcfReadOptions::default();
        let alleles: Vec<&[u8]> = vec![b"c", b"t"];

        let parsed = SnpVcfReader::parse_alleles(&alleles, &options).unwrap();

        assert_eq!(parsed.0, b'C');
        assert_eq!(parsed.1, vec![b'T']);
    }

    #[test]
    fn parse_alleles_accepts_multiallelic_snv() {
        let options = VcfReadOptions::default();
        let alleles: Vec<&[u8]> = vec![b"G", b"A", b"T"];

        let parsed = SnpVcfReader::parse_alleles(&alleles, &options).unwrap();

        assert_eq!(parsed.0, b'G');
        assert_eq!(parsed.1, vec![b'A', b'T']);
    }

    #[test]
    fn parse_alleles_rejects_indel_with_default_options() {
        let options = VcfReadOptions::default();
        let alleles: Vec<&[u8]> = vec![b"A", b"AT"];

        assert!(SnpVcfReader::parse_alleles(&alleles, &options).is_none());
    }

    #[test]
    fn parse_alleles_rejects_non_acgt_with_default_options() {
        let options = VcfReadOptions::default();
        let alleles: Vec<&[u8]> = vec![b"N", b"T"];

        assert!(SnpVcfReader::parse_alleles(&alleles, &options).is_none());
    }

    #[test]
    fn parse_alleles_can_allow_non_snp_lengths() {
        let options = VcfReadOptions {
            pass_only: true,
            snps_only: false,
            acgt_only: true,
        };

        let alleles: Vec<&[u8]> = vec![b"A", b"AT"];

        let parsed = SnpVcfReader::parse_alleles(&alleles, &options).unwrap();

        assert_eq!(parsed.0, b'A');
        assert_eq!(parsed.1, vec![b'A']);
    }

    #[test]
    fn fallback_name_uses_one_based_position() {
        let name = SnpVcfReader::fallback_name("chr17", 7673802, b'C', b"T");

        assert_eq!(name, "chr17:7673803:C>T");
    }

    #[test]
    fn fallback_name_handles_multiallelic_alt_string() {
        let name = SnpVcfReader::fallback_name("chr1", 99, b'G', b"AT");

        assert_eq!(name, "chr1:100:G>AT");
    }
}
