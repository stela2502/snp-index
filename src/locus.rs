//! Final SNP locus objects.
//!
//! This module converts raw VCF/BCF records into stable indexed SNP features.
//! It does not build bins and it does not inspect reads.
use crate::vcf::RawSnpRecord;

/// One indexed SNP locus.
///
/// `id` is the feature id used in sparse matrices.
/// `chr_id` refers to the chromosome id space supplied by the caller,
/// normally derived from the BAM header.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SnpLocus {
    pub id: usize,
    pub chr_id: usize,
    pub pos0: u32,
    pub reference: u8,
    pub alternates: Vec<u8>,
    pub name: String,
    pub vcf_id: String,
}
impl SnpLocus {
    /// Create a new SNP locus.
    pub fn new(
        id: usize,
        chr_id: usize,
        pos0: u32,
        reference: u8,
        alternates: &[u8],
        name: &str,
        vcf_id: &str
    ) -> Self {
        Self {
            id,
            chr_id,
            pos0,
            reference,
            alternates: alternates.to_vec(),
            name: name.to_string(),
            vcf_id: vcf_id.to_string(),
        }
    }

    /// Convert one raw SNP record into an indexed SNP locus.
    pub fn from_raw(id: usize, raw: RawSnpRecord) -> Self {
        Self::new(
            id,
            raw.chr_id,
            raw.pos0,
            raw.reference,
            &raw.alternates,
            &raw.name,
            &raw.vcf_id,
        )
    }

    /// Convert raw SNP records into indexed loci.
    ///
    /// Loci are sorted by `(chr_id, pos0, name)` before ids are assigned.
    /// This gives deterministic feature ids independent of VCF input order.
    pub fn from_raw_records(mut records: Vec<RawSnpRecord>) -> Vec<Self> {
        records.sort_by(|a, b| {
            a.chr_id
                .cmp(&b.chr_id)
                .then_with(|| a.pos0.cmp(&b.pos0))
                .then_with(|| a.name.cmp(&b.name))
        });

        records
            .into_iter()
            .enumerate()
            .map(|(id, raw)| Self::from_raw(id, raw))
            .collect()
    }

    /// Return true if this locus has the given chromosome and position.
    pub fn is_at(&self, chr_id: usize, pos0: u32) -> bool {
        self.chr_id == chr_id && self.pos0 == pos0
    }

    /// Return true if `base` equals the reference allele.
    pub fn is_reference_base(&self, base: u8) -> bool {
        self.reference == base.to_ascii_uppercase()
    }

    /// Return true if `base` matches any alternative allele.
    pub fn is_alternate_base(&self, base: u8) -> bool {
        let base = base.to_ascii_uppercase();
        self.alternates.contains(&base)
    }

    /// Return the 1-based genomic position, useful for human-readable output.
    pub fn pos1(&self) -> u32 {
        self.pos0 + 1
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    impl RawSnpRecord {
        fn test_record(chr_id: usize, pos0: u32, name: &str) -> Self {
            Self {
                chr_id,
                pos0,
                reference: b'C',
                alternates: vec![b'T'],
                name: name.to_string(),
                vcf_id: format!("vcf_{name}"),
            }
        }
    }

    #[test]
    fn new_creates_locus() {
        let locus = SnpLocus::new(
            7,
            1,
            100,
            b'A',
            b"G",
            "chr2:101:A>G",
            "rs_test_7",
        );

        assert_eq!(locus.id, 7);
        assert_eq!(locus.chr_id, 1);
        assert_eq!(locus.pos0, 100);
        assert_eq!(locus.reference, b'A');
        assert_eq!(locus.alternates, vec![b'G']);
        assert_eq!(locus.name, "chr2:101:A>G");
        assert_eq!(locus.vcf_id, "rs_test_7");
    }

    #[test]
    fn from_raw_assigns_id_and_preserves_vcf_id() {
        let raw = RawSnpRecord::test_record(0, 42, "snp_a");

        let locus = SnpLocus::from_raw(3, raw);

        assert_eq!(locus.id, 3);
        assert_eq!(locus.chr_id, 0);
        assert_eq!(locus.pos0, 42);
        assert_eq!(locus.name, "snp_a");
        assert_eq!(locus.vcf_id, "vcf_snp_a");
    }

    #[test]
    fn from_raw_records_sorts_before_assigning_ids() {
        let records = vec![
            RawSnpRecord::test_record(1, 50, "c"),
            RawSnpRecord::test_record(0, 20, "b"),
            RawSnpRecord::test_record(0, 10, "a"),
        ];

        let loci = SnpLocus::from_raw_records(records);

        assert_eq!(loci.len(), 3);

        assert_eq!(loci[0].id, 0);
        assert_eq!(loci[0].chr_id, 0);
        assert_eq!(loci[0].pos0, 10);
        assert_eq!(loci[0].name, "a");
        assert_eq!(loci[0].vcf_id, "vcf_a");

        assert_eq!(loci[1].id, 1);
        assert_eq!(loci[1].chr_id, 0);
        assert_eq!(loci[1].pos0, 20);
        assert_eq!(loci[1].name, "b");
        assert_eq!(loci[1].vcf_id, "vcf_b");

        assert_eq!(loci[2].id, 2);
        assert_eq!(loci[2].chr_id, 1);
        assert_eq!(loci[2].pos0, 50);
        assert_eq!(loci[2].name, "c");
        assert_eq!(loci[2].vcf_id, "vcf_c");
    }

    #[test]
    fn from_raw_records_sorts_same_position_by_name() {
        let records = vec![
            RawSnpRecord::test_record(0, 10, "b"),
            RawSnpRecord::test_record(0, 10, "a"),
        ];

        let loci = SnpLocus::from_raw_records(records);

        assert_eq!(loci[0].id, 0);
        assert_eq!(loci[0].name, "a");
        assert_eq!(loci[0].vcf_id, "vcf_a");

        assert_eq!(loci[1].id, 1);
        assert_eq!(loci[1].name, "b");
        assert_eq!(loci[1].vcf_id, "vcf_b");
    }

    #[test]
    fn is_at_checks_chr_and_position() {
        let locus = SnpLocus::from_raw(0, RawSnpRecord::test_record(2, 123, "snp"));

        assert!(locus.is_at(2, 123));
        assert!(!locus.is_at(2, 124));
        assert!(!locus.is_at(1, 123));
    }

    #[test]
    fn allele_checks_are_case_insensitive() {
        let locus = SnpLocus::new(
            0,
            0,
            10,
            b'C',
            b"AT",
            "snp",
            "rs_snp",
        );

        assert!(locus.is_reference_base(b'C'));
        assert!(locus.is_reference_base(b'c'));

        assert!(locus.is_alternate_base(b'A'));
        assert!(locus.is_alternate_base(b'a'));
        assert!(locus.is_alternate_base(b'T'));
        assert!(locus.is_alternate_base(b't'));

        assert!(!locus.is_reference_base(b'A'));
        assert!(!locus.is_alternate_base(b'G'));
    }

    #[test]
    fn pos1_returns_one_based_position() {
        let locus = SnpLocus::from_raw(0, RawSnpRecord::test_record(0, 99, "snp"));

        assert_eq!(locus.pos1(), 100);
    }
}