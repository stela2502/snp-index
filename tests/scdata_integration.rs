//scdata_integration.rs
use anyhow::Result;
use scdata::cell_data::GeneUmiHash;
use scdata::{MatrixValueType, Scdata};
use snp_index::{AlignedRead, Genome, ReadOpKind, RefineOptions, SnpIndex, Strand, VcfReadOptions};

use std::fs;
use std::path::PathBuf;

#[test]
fn snp_matches_can_be_inserted_into_scdata_and_exported() -> Result<()> {
    let tmp = std::env::temp_dir().join(format!(
        "snp_index_scdata_export_test_{}",
        std::process::id()
    ));

    fs::create_dir_all(&tmp)?;

    let fasta_path = tmp.join("genome.fa");
    let vcf_path = tmp.join("snps.vcf");
    let ref_out = tmp.join("ref_matrix");
    let alt_out = tmp.join("alt_matrix");

    write_test_genome(&fasta_path)?;
    write_test_vcf(&vcf_path)?;

    let genome = Genome::from_fasta(&fasta_path)?;

    let snp_index = SnpIndex::from_vcf_path(
        &vcf_path,
        genome.chr_names(),
        genome.chr_lengths(),
        50,
        &VcfReadOptions::default(),
    )?;

    let mut scdata_ref = Scdata::new(1, MatrixValueType::Real);
    let mut scdata_alt = Scdata::new(1, MatrixValueType::Real);

    // TODO: replace this with your real MappingInfo constructor/import.
    let mut report = mapping_info::MappingInfo::new(None, 0.0, usize::MAX);

    let reads = make_test_reads();

    for (i, mut read) in reads.into_iter().enumerate() {
        read.refine_against_genome(&genome, RefineOptions::default());

        let snp_match = snp_index.match_read(&read, 20);

        let cell_id = (i as u64) + 1;
        let umi = 0_u64;

        for hit in &snp_match.reference {
            let snp_id = hit.feature_id();
            scdata_ref.try_insert(&cell_id, GeneUmiHash(snp_id as u64, umi), 1.0, &mut report);
        }
        
        for hit in &snp_match.alternate {
            let snp_id = hit.feature_id();
            scdata_alt.try_insert(&cell_id, GeneUmiHash(snp_id as u64, umi), 1.0, &mut report);
        }
    }

    // TODO: adjust min_cell_counts if your Scdata export needs more cells/counts.
    scdata_ref.finalize_for_export(1, &snp_index);
    scdata_alt.finalize_for_export(1, &snp_index);

    scdata_ref
        .write_sparse(&ref_out, &snp_index)
        .map_err(anyhow::Error::msg)?;
    scdata_alt
        .write_sparse(&alt_out, &snp_index)
        .map_err(anyhow::Error::msg)?;

    assert!(ref_out.exists());
    assert!(alt_out.exists());

    let ref_features = read_gz_to_string(ref_out.join("features.tsv.gz"))?;
    let alt_features = read_gz_to_string(alt_out.join("features.tsv.gz"))?;

    let ref_matrix = read_gz_to_string(ref_out.join("matrix.mtx.gz"))?;
    let alt_matrix = read_gz_to_string(alt_out.join("matrix.mtx.gz"))?;

    assert!(ref_features.contains("snp_100_A_G"));
    assert!(ref_features.contains("snp_150_C_T"));
    assert!(ref_features.contains("snp_210_G_A"));

    assert!(alt_features.contains("snp_100_A_G"));
    assert!(alt_features.contains("snp_150_C_T"));
    assert!(alt_features.contains("snp_210_G_A"));

    println!("REF features:\n{ref_features}");
    println!("ALT features:\n{alt_features}");
    println!("REF matrix:\n{ref_matrix}");
    println!("ALT matrix:\n{alt_matrix}");

    fs::remove_dir_all(&tmp).ok();

    Ok(())
}

fn write_test_genome(path: &PathBuf) -> Result<()> {
    let mut chr1 = vec![b'N'; 300];

    // SNP positions:
    // pos0 100: REF A, ALT G
    // pos0 150: REF C, ALT T
    // pos0 210: REF G, ALT A

    chr1[100] = b'A';
    chr1[150] = b'C';
    chr1[210] = b'G';

    // Extra exon/intron-like context for spliced/refinement tests.
    for b in &mut chr1[90..110] {
        *b = b'C';
    }
    chr1[100] = b'A';

    for b in &mut chr1[110..210] {
        *b = b'A';
    }

    for b in &mut chr1[210..230] {
        *b = b'G';
    }

    let fasta = format!(
        ">chr1 test chromosome\n{}\n",
        String::from_utf8(chr1).unwrap()
    );

    fs::write(path, fasta)?;

    Ok(())
}

fn write_test_vcf(path: &PathBuf) -> Result<()> {
    fs::write(
        path,
        "\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=300>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t101\tsnp_100_A_G\tA\tG\t.\tPASS\t.
chr1\t151\tsnp_150_C_T\tC\tT\t.\tPASS\t.
chr1\t211\tsnp_210_G_A\tG\tA\t.\tPASS\t.
",
    )?;

    Ok(())
}

fn make_test_reads() -> Vec<AlignedRead> {
    vec![
        // 1: covers SNP 100 as REF A
        make_linear_read(96, 20, &[(100, b'A')]),
        // 2: covers SNP 100 as ALT G
        make_linear_read(96, 20, &[(100, b'G')]),
        // 3: covers SNP 150 as REF C
        make_linear_read(146, 20, &[(150, b'C')]),
        // 4: covers SNP 150 as ALT T
        make_linear_read(146, 20, &[(150, b'T')]),
        // 5: covers SNP 210 as REF G
        make_linear_read(206, 20, &[(210, b'G')]),
        // 6: covers SNP 210 as ALT A
        make_linear_read(206, 20, &[(210, b'A')]),
        // 7: covers SNP 100 and 150: REF at 100, ALT at 150
        make_linear_read(96, 80, &[(100, b'A'), (150, b'T')]),
        // 8: covers SNP 100 and 150: ALT at 100, REF at 150
        make_linear_read(96, 80, &[(100, b'G'), (150, b'C')]),
        // 9: spliced read, second block covers SNP 210 as REF G
        make_spliced_read_around_210(b'G'),
        // 10: spliced read, second block covers SNP 210 as ALT A
        make_spliced_read_around_210(b'A'),
    ]
}

fn make_linear_read(start0: u32, len: u32, edits: &[(u32, u8)]) -> AlignedRead {
    let mut seq = vec![b'N'; len as usize];

    for (pos0, base) in edits {
        let read_pos = (*pos0 - start0) as usize;
        seq[read_pos] = *base;
    }

    AlignedRead::new(
        0,
        Strand::Plus,
        start0,
        seq,
        Some(vec![30; len as usize]),
        vec![(ReadOpKind::Match, len)],
    )
}

fn make_spliced_read_around_210(base_at_210: u8) -> AlignedRead {
    // Reference:
    // first block: 100..110
    // skip:        110..210
    // second block:210..220
    //
    // Read positions:
    // 0..10  first block
    // 10..20 second block
    // read_pos 10 maps to ref pos 210

    let mut seq = vec![b'C'; 20];
    seq[10] = base_at_210;

    for b in &mut seq[11..20] {
        *b = b'G';
    }

    AlignedRead::new(
        0,
        Strand::Plus,
        100,
        seq,
        Some(vec![30; 20]),
        vec![
            (ReadOpKind::Match, 10),
            (ReadOpKind::RefSkip, 100),
            (ReadOpKind::Match, 10),
        ],
    )
}

use flate2::read::GzDecoder;
use std::io::Read;

fn read_gz_to_string(path: impl AsRef<std::path::Path>) -> Result<String> {
    let file = std::fs::File::open(path)?;
    let mut decoder = GzDecoder::new(file);

    let mut text = String::new();
    decoder.read_to_string(&mut text)?;
    Ok(text)
}
