//! Genome-aware read refinement.
//!
//! Refinement must be sequence-aware. A mismatch near a splice junction may be
//! either a true mutation or a mapping artifact. We only rewrite it when the
//! read base matches the reference base at the proposed corrected position.

use crate::genome::Genome;
use crate::read::aligned_read::{AlignedRead, ReadOpKind};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RefineOptions {
    pub max_diff_before_refskip: u32,
    pub max_diff_after_refskip: u32,
    pub merge_adjacent_ops: bool,
}

impl Default for RefineOptions {
    fn default() -> Self {
        Self {
            max_diff_before_refskip: 1,
            max_diff_after_refskip: 1,
            merge_adjacent_ops: true,
        }
    }
}

impl AlignedRead {
    /// Refine this read against a reference genome.
    ///
    /// This is the safe production refinement entry point.
    pub fn refine_against_genome(&mut self, genome: &Genome, options: RefineOptions) {
        let pairs = self.op_pairs();

        let pairs = self.absorb_diff_before_refskip_against_genome(&pairs, genome, options);
        //let pairs = self.absorb_diff_after_refskip_against_genome(&pairs, genome, options);

        let pairs = if options.merge_adjacent_ops {
            Self::merge_adjacent_op_pairs(&pairs)
        } else {
            pairs
        };

        self.replace_ops(pairs);
    }

    /// Convenience wrapper using default options.
    pub fn refine_against_genome_default(&mut self, genome: &Genome) {
        self.refine_against_genome(genome, RefineOptions::default());
    }

    /// Absorb `M X N M` into `M N M` only if the X base matches the first
    /// reference base after the skipped region.
    pub fn absorb_diff_before_refskip_against_genome(
        &self,
        pairs: &[(ReadOpKind, u32)],
        genome: &Genome,
        options: RefineOptions,
    ) -> Vec<(ReadOpKind, u32)> {
        if options.max_diff_before_refskip == 0 || pairs.len() < 4 {
            return pairs.to_vec();
        }

        let ops = Self::build_positioned_ops(self.ref_start0, pairs);
        let mut out = Vec::with_capacity(pairs.len());
        let mut i = 0usize;

        while i < ops.len() {
            if i + 3 < ops.len()
                && Self::is_aligned_kind(ops[i].kind)
                && ops[i + 1].kind == ReadOpKind::Diff
                && ops[i + 1].len <= options.max_diff_before_refskip
                && ops[i + 2].kind == ReadOpKind::RefSkip
                && Self::is_aligned_kind(ops[i + 3].kind)
                && self.diff_matches_reference_at_target(genome, &ops[i + 1], ops[i + 3].ref_start0)
            {
                out.push((ops[i].kind, ops[i].len));
                out.push((ops[i + 2].kind, ops[i + 2].len));
                out.push((ops[i + 3].kind, ops[i + 3].len + ops[i + 1].len));
                i += 4;
            } else {
                out.push((ops[i].kind, ops[i].len));
                i += 1;
            }
        }

        out
    }

    /// Absorb `M N X M` into `M N M` only if the X base matches the last
    /// reference base before the skipped region.
    pub fn absorb_diff_after_refskip_against_genome(
        &self,
        pairs: &[(ReadOpKind, u32)],
        genome: &Genome,
        options: RefineOptions,
    ) -> Vec<(ReadOpKind, u32)> {
        if options.max_diff_after_refskip == 0 || pairs.len() < 4 {
            return pairs.to_vec();
        }

        let ops = Self::build_positioned_ops(self.ref_start0, pairs);
        let mut out = Vec::with_capacity(pairs.len());
        let mut i = 0usize;

        while i < ops.len() {
            if i + 3 < ops.len()
                && Self::is_aligned_kind(ops[i].kind)
                && ops[i + 1].kind == ReadOpKind::RefSkip
                && ops[i + 2].kind == ReadOpKind::Diff
                && ops[i + 2].len <= options.max_diff_after_refskip
                && Self::is_aligned_kind(ops[i + 3].kind)
                && self.diff_matches_reference_at_target(genome, &ops[i + 2], ops[i].ref_end0())
            {
                out.push((ops[i].kind, ops[i].len + ops[i + 2].len));
                out.push((ops[i + 1].kind, ops[i + 1].len));
                out.push((ops[i + 3].kind, ops[i + 3].len));
                i += 4;
            } else {
                out.push((ops[i].kind, ops[i].len));
                i += 1;
            }
        }

        out
    }

    /// Check whether every read base in a Diff operation matches the reference
    /// at the proposed corrected target position.
    pub fn diff_matches_reference_at_target(
        &self,
        genome: &Genome,
        diff_op: &crate::read::aligned_read::ReadOp,
        target_ref_start0: u32,
    ) -> bool {
        for offset in 0..diff_op.len {
            let read_pos = diff_op.read_start + offset;
            let target_pos = target_ref_start0 + offset;

            let Some(read_base) = self.seq.get(read_pos as usize).copied() else {
                return false;
            };

            if !genome.base_matches(self.chr_id, target_pos, read_base) {
                return false;
            }
        }

        true
    }

    pub fn merge_adjacent_op_pairs(pairs: &[(ReadOpKind, u32)]) -> Vec<(ReadOpKind, u32)> {
        let mut out: Vec<(ReadOpKind, u32)> = Vec::with_capacity(pairs.len());

        for pair in pairs {
            if pair.1 == 0 {
                continue;
            }

            match out.last_mut() {
                Some(last) if last.0 == pair.0 => last.1 += pair.1,
                _ => out.push(*pair),
            }
        }

        out
    }

    pub fn is_aligned_kind(kind: ReadOpKind) -> bool {
        matches!(
            kind,
            ReadOpKind::Match | ReadOpKind::Equal | ReadOpKind::Diff
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::read::aligned_read::Strand;

    impl Genome {
        fn synthetic_splice_genome() -> Self {
            // Coordinates:
            // 100..110 = first exon, C
            // 110..210 = intron, A
            // 210..226 = second exon, G
            let mut seq = vec![b'N'; 300];

            for b in &mut seq[100..110] {
                *b = b'C';
            }
            for b in &mut seq[110..210] {
                *b = b'A';
            }
            for b in &mut seq[210..226] {
                *b = b'G';
            }

            Genome::new(vec![("chr1".to_string(), seq)]).unwrap()
        }
    }

    impl AlignedRead {
        fn make_read(seq: &[u8], pairs: Vec<(ReadOpKind, u32)>) -> Self {
            Self::new(
                0,
                Strand::Plus,
                100,
                seq.to_vec(),
                Some(vec![30; seq.len()]),
                pairs,
            )
        }
    }

    #[test]
    fn real_internal_mutation_is_not_refined_without_refskip() {
        let genome = Genome::synthetic_splice_genome();

        let mut read = AlignedRead::make_read(
            b"CCCCCCCCCCTCCCC",
            vec![
                (ReadOpKind::Match, 10),
                (ReadOpKind::Diff, 1),
                (ReadOpKind::Match, 4),
            ],
        );

        read.refine_against_genome(&genome, RefineOptions::default());

        assert_eq!(
            read.op_pairs(),
            vec![
                (ReadOpKind::Match, 10),
                (ReadOpKind::Diff, 1),
                (ReadOpKind::Match, 4),
            ]
        );

        let obs = read.base_at_ref_pos(110).unwrap();
        assert_eq!(obs.base, b'T');
        assert_eq!(obs.read_pos, 10);
    }

    #[test]
    fn diff_before_refskip_is_preserved_if_it_is_real_mutation_not_second_exon_base() {
        let genome = Genome::synthetic_splice_genome();

        let mut read = AlignedRead::make_read(
            b"CCCCCCCCCCTGGGGGGGGGGGGGGG",
            vec![
                (ReadOpKind::Match, 10),
                (ReadOpKind::Diff, 1),
                (ReadOpKind::RefSkip, 100),
                (ReadOpKind::Match, 15),
            ],
        );

        read.refine_against_genome(&genome, RefineOptions::default());

        assert_eq!(
            read.op_pairs(),
            vec![
                (ReadOpKind::Match, 10),
                (ReadOpKind::Diff, 1),
                (ReadOpKind::RefSkip, 100),
                (ReadOpKind::Match, 15),
            ]
        );

        let mutation = read.base_at_ref_pos(110).unwrap();
        assert_eq!(mutation.base, b'T');
        assert_eq!(mutation.read_pos, 10);
    }

    #[test]
    fn diff_after_refskip_is_preserved_if_it_is_not_left_exon_reference() {
        let genome = Genome::synthetic_splice_genome();

        let mut read = AlignedRead::make_read(
            b"CCCCCCCCCCTGGGGGGGGGGGGGGG",
            vec![
                (ReadOpKind::Match, 10),
                (ReadOpKind::RefSkip, 100),
                (ReadOpKind::Diff, 1),
                (ReadOpKind::Match, 15),
            ],
        );

        read.refine_against_genome(&genome, RefineOptions::default());

        assert_eq!(
            read.op_pairs(),
            vec![
                (ReadOpKind::Match, 10),
                (ReadOpKind::RefSkip, 100),
                (ReadOpKind::Diff, 1),
                (ReadOpKind::Match, 15),
            ]
        );

        let mutation = read.base_at_ref_pos(210).unwrap();
        assert_eq!(mutation.base, b'T');
        assert_eq!(mutation.read_pos, 10);
    }

    #[test]
    fn larger_diff_is_not_refined_by_default() {
        let genome = Genome::synthetic_splice_genome();

        let mut read = AlignedRead::make_read(
            b"CCCCCCCCCCGGGGGGGGGGGGGGGGG",
            vec![
                (ReadOpKind::Match, 10),
                (ReadOpKind::Diff, 2),
                (ReadOpKind::RefSkip, 100),
                (ReadOpKind::Match, 15),
            ],
        );

        read.refine_against_genome(&genome, RefineOptions::default());

        assert_eq!(
            read.op_pairs(),
            vec![
                (ReadOpKind::Match, 10),
                (ReadOpKind::Diff, 2),
                (ReadOpKind::RefSkip, 100),
                (ReadOpKind::Match, 15),
            ]
        );
    }

    #[test]
    fn larger_diff_can_be_refined_when_configured_and_reference_matches() {
        let genome = Genome::synthetic_splice_genome();

        let mut read = AlignedRead::make_read(
            b"CCCCCCCCCCGGGGGGGGGGGGGGGGG",
            vec![
                (ReadOpKind::Match, 10),
                (ReadOpKind::Diff, 2),
                (ReadOpKind::RefSkip, 100),
                (ReadOpKind::Match, 15),
            ],
        );

        read.refine_against_genome(
            &genome,
            RefineOptions {
                max_diff_before_refskip: 2,
                max_diff_after_refskip: 1,
                merge_adjacent_ops: true,
            },
        );

        assert_eq!(
            read.op_pairs(),
            vec![
                (ReadOpKind::Match, 10),
                (ReadOpKind::RefSkip, 100),
                (ReadOpKind::Match, 17),
            ]
        );

        assert_eq!(read.base_at_ref_pos(210).unwrap().base, b'G');
        assert_eq!(read.base_at_ref_pos(211).unwrap().base, b'G');
    }

    #[test]
    fn merge_adjacent_op_pairs_collapses_neighbors_and_removes_zero_length_ops() {
        let pairs = vec![
            (ReadOpKind::Match, 5),
            (ReadOpKind::Match, 3),
            (ReadOpKind::Diff, 0),
            (ReadOpKind::RefSkip, 100),
            (ReadOpKind::Match, 2),
            (ReadOpKind::Match, 4),
        ];

        let merged = AlignedRead::merge_adjacent_op_pairs(&pairs);

        assert_eq!(
            merged,
            vec![
                (ReadOpKind::Match, 8),
                (ReadOpKind::RefSkip, 100),
                (ReadOpKind::Match, 6),
            ]
        );
    }
}
