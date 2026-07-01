//! Value-exact compatibility with vcftools 0.1.17 `--TajimaD`.
//!
//! Every `expected` string below was captured verbatim from vcftools 0.1.17 on
//! crafted inputs and frozen as a constant. No oracle binary, interpreter, or
//! filesystem is touched at test time — inputs feed `compute_tajima_d_reader`
//! straight from in-memory bytes.
//!
//! Cases deliberately span the edges the original golden set omitted: missing
//! genotypes (`./.`), half-calls (`0/.`), phased genotypes, bare-`.` haploid
//! sites, all-missing sites, empty/gap and trailing windows (`nan`),
//! monomorphic and multiallelic sites, `ALT=.`, single/two/many-sample files,
//! low frequencies, and a `GT:DP` FORMAT.

use rsomics_vcf_tajima_d::{compute_tajima_d_reader, header};

fn render(vcf: &str, window_size: u64) -> String {
    let rows = compute_tajima_d_reader(vcf.as_bytes(), window_size).expect("compute");
    let mut out = header().to_string();
    for row in &rows {
        out.push_str(&row.to_text());
    }
    out
}

/// (label, window_size, vcf_input, expected_vcftools_output)
const CASES: &[(&str, u64, &str, &str)] = &[
    (
        "missing_gt",
        1000,
        "##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\nchr1\t100\t.\tA\tT\t60\tPASS\t.\tGT\t0/0\t0/1\t./.\nchr1\t200\t.\tG\tC\t60\tPASS\t.\tGT\t0/1\t1/1\t0/0\n",
        "CHROM\tBIN_START\tN_SNPS\tTajimaD\nchr1\t0\t2\t0.941775\n",
    ),
    (
        "half_and_phased",
        1000,
        "##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\nchr1\t100\t.\tA\tT\t60\tPASS\t.\tGT\t0/0\t0/1\t0/.\nchr1\t200\t.\tG\tC\t60\tPASS\t.\tGT\t0|1\t1|1\t0|0\n",
        "CHROM\tBIN_START\tN_SNPS\tTajimaD\nchr1\t0\t2\t0.584729\n",
    ),
    (
        "haploid_skipped",
        1000,
        "##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\nchr1\t100\t.\tA\tT\t60\tPASS\t.\tGT\t0\t1\t0/1\nchr1\t200\t.\tG\tC\t60\tPASS\t.\tGT\t0/1\t1/1\t0/0\n",
        "CHROM\tBIN_START\tN_SNPS\tTajimaD\nchr1\t0\t1\t1.4451\n",
    ),
    (
        "bare_dot_skipped",
        1000,
        "##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\nchr1\t100\t.\tA\tT\t60\tPASS\t.\tGT\t0/1\t.\t1/1\nchr1\t200\t.\tG\tC\t60\tPASS\t.\tGT\t0/1\t1/1\t0/0\n",
        "CHROM\tBIN_START\tN_SNPS\tTajimaD\nchr1\t0\t1\t1.4451\n",
    ),
    (
        "all_missing_leading",
        1000,
        "##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\nchr1\t100\t.\tA\tT\t60\tPASS\t.\tGT\t./.\t./.\t./.\nchr1\t3200\t.\tG\tC\t60\tPASS\t.\tGT\t0/1\t1/1\t0/0\n",
        "CHROM\tBIN_START\tN_SNPS\tTajimaD\nchr1\t3000\t1\t1.4451\n",
    ),
    (
        "multiallelic_skipped",
        1000,
        "##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\nchr1\t100\t.\tA\tT,C\t60\tPASS\t.\tGT\t0/1\t1/2\t0/0\nchr1\t200\t.\tG\tC\t60\tPASS\t.\tGT\t0/1\t1/1\t0/0\n",
        "CHROM\tBIN_START\tN_SNPS\tTajimaD\nchr1\t0\t1\t1.4451\n",
    ),
    (
        "alt_dot_skipped",
        1000,
        "##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\nchr1\t100\t.\tA\t.\t60\tPASS\t.\tGT\t0/0\t0/0\t0/0\nchr1\t200\t.\tG\tC\t60\tPASS\t.\tGT\t0/1\t1/1\t0/0\n",
        "CHROM\tBIN_START\tN_SNPS\tTajimaD\nchr1\t0\t1\t1.4451\n",
    ),
    (
        "gap_windows",
        1000,
        "##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\nchr1\t100\t.\tA\tT\t60\tPASS\t.\tGT\t0/1\t1/1\t0/0\nchr1\t3100\t.\tG\tC\t60\tPASS\t.\tGT\t0/1\t1/1\t0/0\n",
        "CHROM\tBIN_START\tN_SNPS\tTajimaD\nchr1\t0\t1\t1.4451\nchr1\t1000\t0\tnan\nchr1\t2000\t0\tnan\nchr1\t3000\t1\t1.4451\n",
    ),
    (
        "trailing_mono",
        1000,
        "##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\nchr1\t100\t.\tA\tT\t60\tPASS\t.\tGT\t0/1\t1/1\t0/0\nchr1\t3100\t.\tG\tC\t60\tPASS\t.\tGT\t0/0\t0/0\t0/0\n",
        "CHROM\tBIN_START\tN_SNPS\tTajimaD\nchr1\t0\t1\t1.4451\nchr1\t1000\t0\tnan\nchr1\t2000\t0\tnan\nchr1\t3000\t0\tnan\n",
    ),
    (
        "trailing_all_missing",
        1000,
        "##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\nchr1\t100\t.\tA\tT\t60\tPASS\t.\tGT\t0/1\t1/1\t0/0\nchr1\t3100\t.\tG\tC\t60\tPASS\t.\tGT\t./.\t./.\t./.\n",
        "CHROM\tBIN_START\tN_SNPS\tTajimaD\nchr1\t0\t1\t1.4451\nchr1\t1000\t0\tnan\nchr1\t2000\t0\tnan\nchr1\t3000\t0\tnan\n",
    ),
    (
        "leading_mono",
        1000,
        "##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\nchr1\t100\t.\tA\tT\t60\tPASS\t.\tGT\t0/0\t0/0\t0/0\nchr1\t3100\t.\tG\tC\t60\tPASS\t.\tGT\t0/1\t1/1\t0/0\n",
        "CHROM\tBIN_START\tN_SNPS\tTajimaD\nchr1\t3000\t1\t1.4451\n",
    ),
    (
        "single_sample",
        1000,
        "##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\nchr1\t100\t.\tA\tT\t60\tPASS\t.\tGT\t0/1\nchr1\t200\t.\tG\tC\t60\tPASS\t.\tGT\t1/1\n",
        "CHROM\tBIN_START\tN_SNPS\tTajimaD\nchr1\t0\t1\tnan\n",
    ),
    (
        "two_sample",
        1000,
        "##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\nchr1\t100\t.\tA\tT\t60\tPASS\t.\tGT\t0/1\t0/1\nchr1\t200\t.\tG\tC\t60\tPASS\t.\tGT\t0/0\t0/1\n",
        "CHROM\tBIN_START\tN_SNPS\tTajimaD\nchr1\t0\t2\t0.59158\n",
    ),
    (
        "small_freq",
        1000,
        "##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\tS5\tS6\tS7\tS8\tS9\tS10\nchr1\t100\t.\tA\tT\t60\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/1\nchr1\t200\t.\tG\tC\t60\tPASS\t.\tGT\t0/1\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\n",
        "CHROM\tBIN_START\tN_SNPS\tTajimaD\nchr1\t0\t2\t-1.51284\n",
    ),
    (
        "pos_exact_multiple",
        100,
        "##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\nchr1\t100\t.\tA\tT\t60\tPASS\t.\tGT\t0/1\t1/1\t0/0\nchr1\t200\t.\tG\tC\t60\tPASS\t.\tGT\t0/1\t1/1\t0/0\n",
        "CHROM\tBIN_START\tN_SNPS\tTajimaD\nchr1\t100\t1\t1.4451\nchr1\t200\t1\t1.4451\n",
    ),
    (
        "multi_chrom",
        1000,
        "##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\nchr1\t100\t.\tA\tT\t60\tPASS\t.\tGT\t0/1\t1/1\t0/0\t0/1\nchr1\t1500\t.\tG\tC\t60\tPASS\t.\tGT\t0/1\t1/1\t0/0\t1/1\nchr2\t500\t.\tA\tG\t60\tPASS\t.\tGT\t1/1\t0/1\t0/0\t0/0\nchr2\t4500\t.\tT\tA\t60\tPASS\t.\tGT\t0/1\t0/1\t0/1\t0/0\n",
        "CHROM\tBIN_START\tN_SNPS\tTajimaD\nchr1\t0\t1\t1.44416\nchr1\t1000\t1\t1.1665\nchr2\t0\t1\t1.1665\nchr2\t1000\t0\tnan\nchr2\t2000\t0\tnan\nchr2\t3000\t0\tnan\nchr2\t4000\t1\t1.1665\n",
    ),
    (
        "format_gt_dp",
        1000,
        "##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\nchr1\t100\t.\tA\tT\t60\tPASS\t.\tGT:DP\t0|1:30\t1|1:25\t0|0:40\t0/.:10\nchr1\t250\t.\tG\tC\t60\tPASS\t.\tGT:DP\t0/1:30\t./.:0\t1/1:22\t0|1:18\n",
        "CHROM\tBIN_START\tN_SNPS\tTajimaD\nchr1\t0\t2\t1.43081\n",
    ),
    (
        "default_window",
        10000,
        "##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\tS5\nchr1\t500\t.\tA\tT\t60\tPASS\t.\tGT\t0/1\t1/1\t0/0\t0/1\t1/1\nchr1\t8000\t.\tG\tC\t60\tPASS\t.\tGT\t0/0\t0/1\t0/1\t1/1\t0/1\nchr1\t25000\t.\tC\tA\t60\tPASS\t.\tGT\t0/1\t0/0\t1/1\t0/1\t0/0\n",
        "CHROM\tBIN_START\tN_SNPS\tTajimaD\nchr1\t0\t2\t1.74286\nchr1\t10000\t0\tnan\nchr1\t20000\t1\t1.30268\n",
    ),
    (
        "twenty_sample_lowfreq",
        1000,
        "##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\ts4\ts5\ts6\ts7\ts8\ts9\ts10\ts11\ts12\ts13\ts14\ts15\ts16\ts17\ts18\ts19\ts20\nchr1\t100\t.\tA\tT\t60\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/1\nchr1\t300\t.\tG\tC\t60\tPASS\t.\tGT\t0/1\t0/1\t0/1\t0/1\t0/1\t0/1\t0/1\t0/1\t0/1\t0/1\t0/1\t0/1\t0/1\t0/1\t0/1\t0/1\t0/1\t0/1\t0/1\t0/1\nchr1\t700\t.\tC\tG\t60\tPASS\t.\tGT\t0/0\t0/0\t0/1\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\n",
        "CHROM\tBIN_START\tN_SNPS\tTajimaD\nchr1\t0\t3\t-0.285802\n",
    ),
];

// Each case's `expected` was captured from vcftools 0.1.17 `--TajimaD`.
#[test]
fn matches_reference_goldens() {
    for (label, window_size, vcf, expected) in CASES {
        let got = render(vcf, *window_size);
        assert_eq!(&got, expected, "case `{label}` diverges from the oracle");
    }
}

/// vcftools aborts with a fatal error on a polyploid genotype
/// ("Polyploidy found, and not supported by vcftools"); we mirror that by
/// refusing the file rather than emitting a fabricated statistic.
#[test]
fn polyploid_is_rejected() {
    let vcf = "##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\nchr1\t100\t.\tA\tT\t60\tPASS\t.\tGT\t0/0/1\t0/1\t1/1\n";
    assert!(compute_tajima_d_reader(vcf.as_bytes(), 1000).is_err());
}
