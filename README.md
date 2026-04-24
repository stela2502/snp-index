[![Rust](https://github.com/stela2502/snp-index/actions/workflows/rust.yml/badge.svg)](https://github.com/stela2502/snp-index/actions/workflows/rust.yml)

# snp-index

A lightweight, high-performance Rust crate for SNP-aware read matching and sparse matrix generation.

## Overview

`snp-index` provides the core building blocks to:

* Load SNPs from a VCF file into a fast lookup index
* Convert BAM records into a normalized `AlignedRead`
* Optionally refine alignments using a reference genome
* Match reads against SNP loci (REF vs ALT)
* Aggregate results into **cell × SNP sparse matrices**
* Export results in **10x-compatible format**

---

## Key Features

* **Fast SNP lookup** via binned index
* **BAM → AlignedRead conversion built-in**
* **Sequence-aware refinement** (genome-based)
* **UMI-aware deduplication** via `scdata`
* **Spliced read support** (`N` / ref-skip)
* **10x-compatible output**

---

## Architecture

```id="pipeline"
BAM → AlignedRead::from_record()
FASTA → Genome (optional, for refinement)
VCF → SnpIndex

AlignedRead
    ↓ refine_against_genome()
SnpIndex.match_read()
    ↓
SnpReadMatch
    ↓
Scdata (cell × SNP × UMI)
    ↓
write_sparse()
```

---

## Example Workflow

```rust id="example"
use snp_index::{
    AlignedRead, Genome, RefineOptions, SnpIndex,
};
use scdata::cell_data::GeneUmiHash;

// Load inputs
let genome = Genome::from_fasta("genome.fa")?;
let snp_index = SnpIndex::from_vcf_path(...)?;

// Iterate BAM (pseudo-code)
for rec in bam.records() {
    let rec = rec?;

    // convert BAM → AlignedRead
    let mut read = AlignedRead::from_record(&rec, chr_id);

    // optional but recommended
    read.refine_against_genome(&genome, RefineOptions::default());

    let hits = snp_index.match_read(&read, 20);

    for snp_id in hits.ref_ids {
        scdata_ref.try_insert(
            &cell,
            GeneUmiHash(snp_id, umi),
            1.0,
            &mut report,
        );
    }

    for snp_id in hits.alt_ids {
        scdata_alt.try_insert(
            &cell,
            GeneUmiHash(snp_id, umi),
            1.0,
            &mut report,
        );
    }
}

// export 10x-style matrices
scdata_ref.write_sparse("ref_out", &snp_index)?;
scdata_alt.write_sparse("alt_out", &snp_index)?;
```

---

## Core Data Structures

### `AlignedRead`

A normalized alignment representation:

* sequence + qualities
* CIGAR-like operations
* explicit reference coordinates
* independent of BAM after construction

```rust id="aligned_read_example"
let read = AlignedRead::new(
    chr_id,
    strand,
    start0,
    seq,
    qual,
    vec![(ReadOpKind::Match, 50)],
);
```

---

### `SnpIndex`

* stores SNP loci
* fast lookup via genomic bins
* implements `FeatureIndex` for `scdata`

---

### `SnpReadMatch`

Result of matching a read:

```rust id="snp_match"
pub struct SnpReadMatch {
    pub ref_ids: Vec<u32>,
    pub alt_ids: Vec<u32>,
}
```

---

## Refinement (Important)

Refinement fixes splice-edge artifacts such as:

```text
10M1X100N15M → 10M100N16M
```

**Only if the read base matches the reference genome at the corrected position.**

This prevents destroying real mutations.

---

## Known Limitations

### Allele conflicts per UMI

If the same `(cell, SNP, UMI)` supports both REF and ALT:

* current behavior: effectively first-observation / duplicate suppression
* no explicit conflict resolution

Possible future strategies:

* majority vote per UMI
* discard conflicting UMIs
* separate “ambiguous” matrix
* quality-weighted decisions

---

## Design Philosophy

* Keep SNP matching **fast and pure**
* Decouple logic from file formats
* Make everything **testable with synthetic data**
* Push statistical decisions downstream

---

## Future Extensions

* multi-allelic SNP tracking
* UMI-level consensus refinement
* conflict-aware allele assignment
* SNP coverage / depth metrics
* performance optimizations (SIMD)

---

## Status

Experimental but fully functional.
Includes full integration tests (VCF + genome + reads → sparse output).

---

## Author

Stefan Lang
Lund University – Bioinformatics Core Facility

