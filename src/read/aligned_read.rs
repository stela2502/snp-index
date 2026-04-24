//! Aligned read representation used for SNP matching and refinement.
//!
//! This struct provides a normalized, sequence-aware view of an alignment,
//! independent of the original input format. While it can be constructed from
//! BAM records (`rust-htslib`), it is designed to act as a stable intermediate
//! representation for downstream processing (e.g. SNP matching, refinement).

use rust_htslib::bam::record::{Cigar, CigarStringView};
use rust_htslib::bam::Record;
use gtf_splice_index::types::RefBlock;

/// Read strand/orientation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Plus,
    Minus,
    Unknown,
}

/// CIGAR-like operation kind.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReadOpKind {
    Match,
    Equal,
    Diff,
    Ins,
    Del,
    RefSkip,
    SoftClip,
    HardClip,
    Pad,
}

/// One alignment operation with explicit reference/read coordinates.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ReadOp {
    pub kind: ReadOpKind,
    pub len: u32,
    pub ref_start0: u32,
    pub read_start: u32,
}

/// A base observed in the read at a genomic reference position.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ObservedBase {
    pub base: u8,
    pub qual: Option<u8>,
    pub read_pos: u32,
}

/// BAM-independent aligned read.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AlignedRead {
    pub chr_id: usize,
    pub strand: Strand,
    pub ref_start0: u32,
    pub seq: Vec<u8>,
    pub qual: Option<Vec<u8>>,
    pub ops: Vec<ReadOp>,
    finalized: bool,
}

impl ReadOpKind {
    /// Return true if this operation consumes reference coordinates.
    pub fn consumes_ref(&self) -> bool {
        matches!(
            self,
            Self::Match | Self::Equal | Self::Diff | Self::Del | Self::RefSkip
        )
    }

    /// Return true if this operation consumes read/query coordinates.
    pub fn consumes_read(&self) -> bool {
        matches!(
            self,
            Self::Match | Self::Equal | Self::Diff | Self::Ins | Self::SoftClip
        )
    }

    /// Return true if this operation has aligned read bases.
    pub fn aligned_bases(&self) -> bool {
        matches!(self, Self::Match | Self::Equal | Self::Diff)
    }
}

impl ReadOp {
    /// Create a positioned read operation.
    pub fn new(kind: ReadOpKind, len: u32, ref_start0: u32, read_start: u32) -> Self {
        Self {
            kind,
            len,
            ref_start0,
            read_start,
        }
    }

    /// Return the exclusive reference end coordinate.
    pub fn ref_end0(&self) -> u32 {
        if self.kind.consumes_ref() {
            self.ref_start0 + self.len
        } else {
            self.ref_start0
        }
    }

    /// Return the exclusive read/query end coordinate.
    pub fn read_end(&self) -> u32 {
        if self.kind.consumes_read() {
            self.read_start + self.len
        } else {
            self.read_start
        }
    }

    /// Return true if this operation can yield a base for a reference position.
    pub fn can_observe_reference_base(&self) -> bool {
        self.kind.aligned_bases()
    }

    /// Return true if this operation covers `pos0` on the reference.
    pub fn contains_ref_pos(&self, pos0: u32) -> bool {
        self.kind.consumes_ref() && self.ref_start0 <= pos0 && pos0 < self.ref_end0()
    }

    /// Map reference coordinate to read coordinate for aligned-base operations.
    pub fn read_pos_for_ref_pos(&self, pos0: u32) -> Option<u32> {
        if !self.can_observe_reference_base() || !self.contains_ref_pos(pos0) {
            return None;
        }

        Some(self.read_start + (pos0 - self.ref_start0))
    }
}

impl ObservedBase {
    /// Create a new observed base.
    pub fn new(base: u8, qual: Option<u8>, read_pos: u32) -> Self {
        Self {
            base: base.to_ascii_uppercase(),
            qual,
            read_pos,
        }
    }
}

impl AlignedRead {
    /// Create a new aligned read from raw pieces.
    ///
    /// `ops_input` contains only `(ReadOpKind, len)`. This constructor walks the
    /// operations and fills reference/read coordinates.
    pub fn new(
        chr_id: usize,
        strand: Strand,
        ref_start0: u32,
        seq: Vec<u8>,
        qual: Option<Vec<u8>>,
        ops_input: Vec<(ReadOpKind, u32)>,
    ) -> Self {
        let ops = Self::build_positioned_ops(ref_start0, &ops_input);

        Self {
            chr_id,
            strand,
            ref_start0,
            seq: Self::uppercase_seq(seq),
            qual,
            ops,
            finalized: false,
        }
    }

    /// Build positioned operations from plain operation/length pairs.
    pub fn build_positioned_ops(ref_start0: u32, ops_input: &[(ReadOpKind, u32)]) -> Vec<ReadOp> {
        let mut ref_pos = ref_start0;
        let mut read_pos = 0u32;
        let mut ops = Vec::with_capacity(ops_input.len());

        for (kind, len) in ops_input {
            ops.push(ReadOp::new(*kind, *len, ref_pos, read_pos));

            if kind.consumes_ref() {
                ref_pos += *len;
            }

            if kind.consumes_read() {
                read_pos += *len;
            }
        }

        ops
    }

    /// Build an `AlignedRead` from a BAM record.
    ///
    /// `chr_id` must already be mapped into the same chromosome id space used by
    /// `SnpIndex`.
    pub fn from_record(record: &Record, chr_id: usize) -> Self {
        let strand = if record.is_reverse() {
            Strand::Minus
        } else {
            Strand::Plus
        };

        Self::new(
            chr_id,
            strand,
            record.pos() as u32,
            record.seq().as_bytes(),
            Some(record.qual().to_vec()),
            Self::cigar_to_read_ops(&record.cigar()),
        )
    }

    /// Convert aligned read bases into genomic reference blocks.
    ///
    /// Only operations that contain aligned read bases are emitted:
    /// `Match`, `Equal`, and `Diff`.
    ///
    /// Deletions and reference skips consume reference coordinates but do not
    /// produce read-supported blocks.
    pub fn ref_blocks(&self) -> Vec<RefBlock> {
        let mut blocks: Vec<RefBlock> = Vec::new();

        for op in &self.ops {
            if !op.kind.aligned_bases() {
                continue;
            }

            let new_block = RefBlock {
                start: op.ref_start0,
                end: op.ref_end0(),
            };

            match blocks.last_mut() {
                Some(last) if new_block.start <= last.end => {
                    // Merge adjacent or overlapping (robust, even if overlap shouldn't happen)
                    last.end = last.end.max(new_block.end);
                }
                _ => blocks.push(new_block),
            }
        }

        blocks
    }

    /// Convert a BAM CIGAR into `AlignedRead` operation pairs.
    fn cigar_to_read_ops(cigar: &CigarStringView) -> Vec<(ReadOpKind, u32)> {
        cigar
            .iter()
            .map(|op| Self::cigar_op_to_read_op(*op))
            .collect()
    }

    /// Convert one BAM CIGAR op into a `ReadOpKind`.
    pub fn cigar_op_to_read_op(op: Cigar) -> (ReadOpKind, u32) {
        match op {
            Cigar::Match(len) => (ReadOpKind::Match, len),
            Cigar::Ins(len) => (ReadOpKind::Ins, len),
            Cigar::Del(len) => (ReadOpKind::Del, len),
            Cigar::RefSkip(len) => (ReadOpKind::RefSkip, len),
            Cigar::SoftClip(len) => (ReadOpKind::SoftClip, len),
            Cigar::HardClip(len) => (ReadOpKind::HardClip, len),
            Cigar::Pad(len) => (ReadOpKind::Pad, len),
            Cigar::Equal(len) => (ReadOpKind::Equal, len),
            Cigar::Diff(len) => (ReadOpKind::Diff, len),
        }
    }

    /// Validate the read and mark it finalized.
    ///
    /// Later refinement steps can also call this after changing operations.
    pub fn finalize(&mut self) -> Result<(), String> {
        self.validate()?;
        self.finalized = true;
        Ok(())
    }

    /// Return whether this read has been finalized.
    pub fn is_finalized(&self) -> bool {
        self.finalized
    }

    /// Validate sequence, quality, and operation consistency.
    pub fn validate(&self) -> Result<(), String> {
        if let Some(qual) = &self.qual {
            if qual.len() != self.seq.len() {
                return Err(format!(
                    "quality length ({}) does not match sequence length ({})",
                    qual.len(),
                    self.seq.len()
                ));
            }
        }

        let expected = self.read_len_from_ops() as usize;
        if expected != self.seq.len() {
            return Err(format!(
                "read ops consume {expected} read bases, but sequence length is {}",
                self.seq.len()
            ));
        }

        Ok(())
    }

    /// Return read length implied by operations.
    pub fn read_len_from_ops(&self) -> u32 {
        self.ops.last().map(|op| op.read_end()).unwrap_or(0)
    }

    /// Return full reference span `[start0, end0)`.
    pub fn ref_span(&self) -> Option<(u32, u32)> {
        let start = self
            .ops
            .iter()
            .filter(|op| op.kind.consumes_ref())
            .map(|op| op.ref_start0)
            .min()?;

        let end = self
            .ops
            .iter()
            .filter(|op| op.kind.consumes_ref())
            .map(|op| op.ref_end0())
            .max()?;

        Some((start, end))
    }

    /// Return observed base at reference position `pos0`.
    ///
    /// Returns `None` for deletions, ref-skips, insertions, clips, pads, and
    /// positions outside the alignment.
    pub fn base_at_ref_pos(&self, pos0: u32) -> Option<ObservedBase> {
        for op in &self.ops {
            let read_pos = match op.read_pos_for_ref_pos(pos0) {
                Some(read_pos) => read_pos,
                None => continue,
            };

            let base = *self.seq.get(read_pos as usize)?;
            let qual = self
                .qual
                .as_ref()
                .and_then(|q| q.get(read_pos as usize).copied());

            return Some(ObservedBase::new(base, qual, read_pos));
        }

        None
    }

    /// Replace operations from plain operation/length pairs.
    ///
    /// This is useful for refinement modules that rewrite CIGAR-like structure.
    pub fn replace_ops(&mut self, ops_input: Vec<(ReadOpKind, u32)>) {
        self.ops = Self::build_positioned_ops(self.ref_start0, &ops_input);
        self.finalized = false;
    }

    /// Return operations as plain `(kind, len)` pairs.
    pub fn op_pairs(&self) -> Vec<(ReadOpKind, u32)> {
        self.ops.iter().map(|op| (op.kind, op.len)).collect()
    }

    /// Uppercase sequence bases.
    pub fn uppercase_seq(seq: Vec<u8>) -> Vec<u8> {
        seq.into_iter().map(|b| b.to_ascii_uppercase()).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use gtf_splice_index::types::RefBlock;


    impl AlignedRead {
        fn simple_read() -> Self {
            Self::new(
                0,
                Strand::Plus,
                100,
                b"ACGTACGTAA".to_vec(),
                Some(vec![30; 10]),
                vec![(ReadOpKind::Match, 10)],
            )
        }
    }

    #[test]
    fn build_positioned_ops_tracks_coordinates() {
        let ops = AlignedRead::build_positioned_ops(
            100,
            &[
                (ReadOpKind::SoftClip, 5),
                (ReadOpKind::Match, 10),
                (ReadOpKind::Ins, 2),
                (ReadOpKind::RefSkip, 100),
                (ReadOpKind::Match, 15),
            ],
        );

        assert_eq!(ops[0], ReadOp::new(ReadOpKind::SoftClip, 5, 100, 0));
        assert_eq!(ops[1], ReadOp::new(ReadOpKind::Match, 10, 100, 5));
        assert_eq!(ops[2], ReadOp::new(ReadOpKind::Ins, 2, 110, 15));
        assert_eq!(ops[3], ReadOp::new(ReadOpKind::RefSkip, 100, 110, 17));
        assert_eq!(ops[4], ReadOp::new(ReadOpKind::Match, 15, 210, 17));
    }

    #[test]
    fn validate_accepts_consistent_read() {
        let mut read = AlignedRead::simple_read();

        assert!(read.finalize().is_ok());
        assert!(read.is_finalized());
    }

    #[test]
    fn validate_rejects_quality_length_mismatch() {
        let read = AlignedRead::new(
            0,
            Strand::Plus,
            100,
            b"ACGT".to_vec(),
            Some(vec![30; 3]),
            vec![(ReadOpKind::Match, 4)],
        );

        assert!(read.validate().is_err());
    }

    #[test]
    fn validate_rejects_sequence_length_mismatch() {
        let read = AlignedRead::new(
            0,
            Strand::Plus,
            100,
            b"ACGT".to_vec(),
            None,
            vec![(ReadOpKind::Match, 3)],
        );

        assert!(read.validate().is_err());
    }

    #[test]
    fn ref_span_includes_skips_and_deletions() {
        let read = AlignedRead::new(
            0,
            Strand::Plus,
            100,
            b"AAAAATTTTT".to_vec(),
            None,
            vec![
                (ReadOpKind::SoftClip, 5),
                (ReadOpKind::Match, 5),
                (ReadOpKind::RefSkip, 100),
                (ReadOpKind::Match, 5),
            ],
        );

        assert_eq!(read.ref_span(), Some((100, 210)));
    }

    #[test]
    fn base_at_ref_pos_maps_match_positions() {
        let read = AlignedRead::simple_read();

        let obs = read.base_at_ref_pos(100).unwrap();
        assert_eq!(obs.base, b'A');
        assert_eq!(obs.qual, Some(30));
        assert_eq!(obs.read_pos, 0);

        let obs = read.base_at_ref_pos(103).unwrap();
        assert_eq!(obs.base, b'T');
        assert_eq!(obs.read_pos, 3);

        assert!(read.base_at_ref_pos(99).is_none());
        assert!(read.base_at_ref_pos(110).is_none());
    }

    #[test]
    fn base_at_ref_pos_skips_deletions_and_refskips() {
        let read = AlignedRead::new(
            0,
            Strand::Plus,
            100,
            b"AAAAACCCCC".to_vec(),
            Some(vec![40; 10]),
            vec![
                (ReadOpKind::Match, 5),
                (ReadOpKind::Del, 3),
                (ReadOpKind::RefSkip, 100),
                (ReadOpKind::Match, 5),
            ],
        );

        assert!(read.base_at_ref_pos(102).is_some());
        assert!(read.base_at_ref_pos(105).is_none());
        assert!(read.base_at_ref_pos(108).is_none());
        assert!(read.base_at_ref_pos(208).is_some());
    }

    #[test]
    fn replace_ops_rebuilds_coordinates() {
        let mut read = AlignedRead::simple_read();

        read.replace_ops(vec![
            (ReadOpKind::Match, 5),
            (ReadOpKind::RefSkip, 100),
            (ReadOpKind::Match, 5),
        ]);

        assert_eq!(read.ref_span(), Some((100, 210)));
        assert!(!read.is_finalized());
    }

    #[test]
    fn ref_blocks_merges_adjacent_aligned_ops() {
        let read = AlignedRead::new(
            0,
            Strand::Plus,
            100,
            b"CCCCCCCCCCTGGGGGGGGGG".to_vec(),
            Some(vec![30; 21]),
            vec![
                (ReadOpKind::Match, 10),
                (ReadOpKind::Diff, 1),
                (ReadOpKind::Match, 10),
            ],
        );

        assert_eq!(
            read.ref_blocks(),
            vec![RefBlock {
                start: 100,
                end: 121,
            }]
        );
    }

    #[test]
    fn ref_blocks_does_not_merge_across_refskip_or_deletion() {
        let read = AlignedRead::new(
            0,
            Strand::Plus,
            100,
            b"CCCCCCCCCCGGGGGGGGGG".to_vec(),
            Some(vec![30; 20]),
            vec![
                (ReadOpKind::Match, 10),
                (ReadOpKind::RefSkip, 100),
                (ReadOpKind::Match, 5),
                (ReadOpKind::Del, 3),
                (ReadOpKind::Match, 5),
            ],
        );

        assert_eq!(
            read.ref_blocks(),
            vec![
                RefBlock {
                    start: 100,
                    end: 110,
                },
                RefBlock {
                    start: 210,
                    end: 215,
                },
                RefBlock {
                    start: 218,
                    end: 223,
                },
            ]
        );
    }
}
