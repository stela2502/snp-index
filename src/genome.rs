//! Reference genome access.
//!
//! This module provides a small FASTA-backed genome object used for
//! sequence-aware read refinement and SNP validation.

use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// In-memory reference genome.
///
/// Chromosome ids are assigned in FASTA order.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Genome {
    pub chr_names: Vec<String>,
    pub seqs: Vec<Vec<u8>>,
    pub name_to_id: HashMap<String, usize>,
}

impl Genome {
    /// Read a plain or gzip-compressed FASTA file into memory.
    ///
    /// Supports:
    /// - `.fa` / `.fasta`
    /// - `.fa.gz` / `.fasta.gz`
    pub fn from_fasta<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path_ref = path.as_ref();
        let reader = Self::open_fasta_reader(path_ref)?;

        let mut records: Vec<(String, Vec<u8>)> = Vec::new();
        let mut current_name: Option<String> = None;
        let mut current_seq: Vec<u8> = Vec::new();

        for line_result in reader.lines() {
            let line = line_result.with_context(|| {
                format!("failed to read FASTA line from {}", path_ref.display())
            })?;

            if let Some(rest) = line.strip_prefix('>') {
                if let Some(name) = current_name.take() {
                    records.push((name, std::mem::take(&mut current_seq)));
                }

                let name = rest.split_whitespace().next().unwrap_or("").to_string();

                if name.is_empty() {
                    anyhow::bail!("empty FASTA record name");
                }

                current_name = Some(name);
            } else if !line.trim().is_empty() {
                current_seq.extend(line.trim().as_bytes());
            }
        }

        if let Some(name) = current_name.take() {
            records.push((name, current_seq));
        }

        if records.is_empty() {
            anyhow::bail!("FASTA file contains no records: {}", path_ref.display());
        }

        Self::new(records)
    }

    /// Open plain or gzip-compressed FASTA as a buffered reader.
    fn open_fasta_reader<P: AsRef<Path>>(path: P) -> Result<Box<dyn BufRead>> {
        let path_ref = path.as_ref();

        let file = File::open(path_ref)
            .with_context(|| format!("failed to open FASTA file {}", path_ref.display()))?;

        let is_gz = path_ref.extension().map(|ext| ext == "gz").unwrap_or(false);

        if is_gz {
            let decoder = MultiGzDecoder::new(file);
            Ok(Box::new(BufReader::new(decoder)))
        } else {
            Ok(Box::new(BufReader::new(file)))
        }
    }
    /// Create a genome directly from named sequences.
    ///
    /// Useful for tests and small synthetic references.
    pub fn new(records: Vec<(String, Vec<u8>)>) -> Result<Self> {
        let mut chr_names = Vec::with_capacity(records.len());
        let mut seqs = Vec::with_capacity(records.len());
        let mut name_to_id = HashMap::with_capacity(records.len());

        for (name, seq) in records {
            if name.is_empty() {
                anyhow::bail!("empty chromosome name");
            }

            if name_to_id.contains_key(&name) {
                anyhow::bail!("duplicate chromosome name: {name}");
            }

            let chr_id = chr_names.len();
            name_to_id.insert(name.clone(), chr_id);
            chr_names.push(name);
            seqs.push(Self::uppercase_sequence(&seq));
        }

        Ok(Self {
            chr_names,
            seqs,
            name_to_id,
        })
    }

    /// Return chromosome id for a chromosome name.
    pub fn chr_id(&self, name: &str) -> Option<usize> {
        self.name_to_id.get(name).copied()
    }

    /// Return chromosome name for an id.
    pub fn chr_name(&self, chr_id: usize) -> Option<&str> {
        self.chr_names.get(chr_id).map(|s| s.as_str())
    }

    /// Return chromosome length.
    pub fn chr_len(&self, chr_id: usize) -> Option<u32> {
        self.seqs.get(chr_id).map(|seq| seq.len() as u32)
    }

    /// Return all chromosome lengths in chromosome id order.
    pub fn chr_lengths(&self) -> Vec<u32> {
        self.seqs.iter().map(|seq| seq.len() as u32).collect()
    }

    /// Return all chromosome names in chromosome id order.
    pub fn chr_names(&self) -> Vec<String> {
        self.chr_names.clone()
    }

    /// Return reference base at `chr_id:pos0`.
    pub fn base(&self, chr_id: usize, pos0: u32) -> Option<u8> {
        self.seqs
            .get(chr_id)
            .and_then(|seq| seq.get(pos0 as usize))
            .copied()
            .map(|b| b.to_ascii_uppercase())
    }

    /// Return reference slice `[start0, end0)`.
    pub fn slice(&self, chr_id: usize, start0: u32, end0: u32) -> Option<&[u8]> {
        if start0 > end0 {
            return None;
        }

        self.seqs
            .get(chr_id)
            .and_then(|seq| seq.get(start0 as usize..end0 as usize))
    }

    /// Return true if `base` matches the reference at `chr_id:pos0`.
    pub fn base_matches(&self, chr_id: usize, pos0: u32, base: u8) -> bool {
        self.base(chr_id, pos0)
            .map(|ref_base| ref_base == base.to_ascii_uppercase())
            .unwrap_or(false)
    }

    /// Parse a FASTA/FASTQ record id into a chromosome name.
    pub fn record_name(id: &[u8]) -> Result<String> {
        let id_str = std::str::from_utf8(id).context("FASTA record id is not valid UTF-8")?;

        let name = id_str.split_whitespace().next().unwrap_or("").to_string();

        if name.is_empty() {
            anyhow::bail!("empty FASTA record id");
        }

        Ok(name)
    }

    /// Uppercase a sequence.
    pub fn uppercase_sequence(seq: &[u8]) -> Vec<u8> {
        seq.iter().map(|b| b.to_ascii_uppercase()).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::AlignedRead;
    use crate::ReadOpKind;
    use crate::Strand;


    #[test]
    fn new_builds_genome_and_uppercases_sequences() {
        let genome = Genome::new(vec![
            ("chr1".to_string(), b"acgtN".to_vec()),
            ("chr2".to_string(), b"ttcc".to_vec()),
        ])
        .unwrap();

        assert_eq!(genome.chr_id("chr1"), Some(0));
        assert_eq!(genome.chr_id("chr2"), Some(1));
        assert_eq!(genome.chr_name(0), Some("chr1"));
        assert_eq!(genome.chr_len(0), Some(5));
        assert_eq!(genome.chr_lengths(), vec![5, 4]);

        assert_eq!(genome.base(0, 0), Some(b'A'));
        assert_eq!(genome.base(0, 1), Some(b'C'));
        assert_eq!(genome.base(0, 4), Some(b'N'));
    }

    #[test]
    fn new_rejects_duplicate_names() {
        let result = Genome::new(vec![
            ("chr1".to_string(), b"AAAA".to_vec()),
            ("chr1".to_string(), b"CCCC".to_vec()),
        ]);

        assert!(result.is_err());
    }

    #[test]
    fn base_returns_none_outside_reference() {
        let genome = Genome::new(vec![("chr1".to_string(), b"ACGT".to_vec())]).unwrap();

        assert_eq!(genome.base(0, 4), None);
        assert_eq!(genome.base(1, 0), None);
    }

    #[test]
    fn slice_returns_half_open_reference_interval() {
        let genome = Genome::new(vec![("chr1".to_string(), b"ACGTACGT".to_vec())]).unwrap();

        assert_eq!(genome.slice(0, 2, 6), Some(&b"GTAC"[..]));
        assert_eq!(genome.slice(0, 6, 2), None);
        assert_eq!(genome.slice(0, 0, 99), None);
    }

    #[test]
    fn base_matches_is_case_insensitive_for_query_base() {
        let genome = Genome::new(vec![("chr1".to_string(), b"ACGT".to_vec())]).unwrap();

        assert!(genome.base_matches(0, 1, b'C'));
        assert!(genome.base_matches(0, 1, b'c'));
        assert!(!genome.base_matches(0, 1, b'A'));
        assert!(!genome.base_matches(0, 99, b'A'));
    }

    #[test]
    fn record_name_uses_first_whitespace_token() {
        let name = Genome::record_name(b"chr1 some description here").unwrap();

        assert_eq!(name, "chr1");
    }

    #[test]
    fn record_name_rejects_empty_id() {
        assert!(Genome::record_name(b"").is_err());
        assert!(Genome::record_name(b"   ").is_err());
    }

    fn make_read(
        seq: &[u8],
        pairs: &[(ReadOpKind, u32)],
    ) -> AlignedRead {
        AlignedRead::new(
            0,
            Strand::Plus,
            0,
            seq,
            Some(vec![30; seq.len()]),
            pairs,
        )
    }

    #[test]
    fn base_at_ref_pos_simple() {
        let read = make_read(
            b"ACGT",
            &[(ReadOpKind::Match, 4)],
        );

        // Assuming default pos0 = 0
        assert_eq!(read.base_at_ref_pos(0).unwrap().base, b'A');
        assert_eq!(read.base_at_ref_pos(1).unwrap().base, b'C');
    }

    #[test]
    fn base_at_ref_pos_softclip() {
        let read = make_read(
            b"NNNNACGT",
            &[
                (ReadOpKind::SoftClip, 4),
                (ReadOpKind::Match, 4),
            ],
        );

        // first aligned base should be 'A'
        assert_eq!(read.base_at_ref_pos(0).unwrap().base, b'A');
    }

    #[test]
    fn base_at_ref_pos_splice() {
        let read = make_read(
            b"AAAACCCC",
            &[
                (ReadOpKind::Match, 4),
                (ReadOpKind::RefSkip, 100),
                (ReadOpKind::Match, 4),
            ],
        );

        assert_eq!(read.base_at_ref_pos(0).unwrap().base, b'A');
        assert_eq!(read.base_at_ref_pos(3).unwrap().base, b'A');

        // skipped region → no base
        assert!(read.base_at_ref_pos(50).is_none());

        // second block starts at 104
        assert_eq!(read.base_at_ref_pos(104).unwrap().base, b'C');
    } 
    
    #[test]
    fn base_at_ref_pos_insertion() {
        let read = make_read(
            b"AACGT",
            &[
                (ReadOpKind::Match, 2),
                (ReadOpKind::Ins, 1),
                (ReadOpKind::Match, 2),
            ],
        );

        assert_eq!(read.base_at_ref_pos(0).unwrap().base, b'A');
        assert_eq!(read.base_at_ref_pos(1).unwrap().base, b'A');
        assert_eq!(read.base_at_ref_pos(2).unwrap().base, b'G');
    }      
}
