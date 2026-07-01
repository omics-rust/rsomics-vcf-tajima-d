//! Value-exact compatibility with vcftools 0.1.17 --TajimaD.
//!
//! All expected output is frozen from black-box observation of vcftools 0.1.17.
//! No vcftools binary is required at test time; a second section gates on
//! vcftools being on PATH (version 0.1.17) for live oracle comparison.

use std::io::Write;
use std::path::PathBuf;
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

use rsomics_vcf_tajima_d::{compute_tajima_d, header};

/// Three-sample VCF: two sites in bin 0 (pos 100, 200), one in bin 1000 (pos 1100).
const TEST_VCF: &str = "\
##fileformat=VCFv4.1\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n\
chr1\t100\t.\tA\tT\t60\tPASS\t.\tGT\t0/0\t0/1\t1/1\n\
chr1\t200\t.\tG\tC\t60\tPASS\t.\tGT\t0/1\t1/1\t0/0\n\
chr1\t1100\t.\tC\tT\t60\tPASS\t.\tGT\t0/0\t0/0\t0/1\n\
chr2\t500\t.\tA\tG\t60\tPASS\t.\tGT\t1/1\t0/1\t0/0\n\
";

/// Write VCF to a unique temp file; each call gets a distinct name.
fn write_vcf(vcf: &str) -> PathBuf {
    let ts = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.subsec_nanos())
        .unwrap_or(0);
    let tid = std::thread::current().id();
    let name = format!("rsomics_vcf_tajima_d_{tid:?}_{ts}.vcf");
    let path = std::env::temp_dir().join(name);
    let mut f = std::fs::File::create(&path).expect("create temp VCF");
    f.write_all(vcf.as_bytes()).expect("write");
    path
}

// chr1 bin 0 (pos 100, 200): S=2, pi=1.2, D=1.75324
// chr1 bin 1000 (pos 1100): S=1, pi=1/3, D=-0.933021
// chr2 bin 0 (pos 500): S=1, pi=0.6, D=1.4451
//
// chr2 site 500: S1=1/1, S2=0/1, S3=0/0 → REF=3, ALT=3, p=0.5
// pi = (6/5)*2*0.5*0.5=0.6; S=1, D = (0.6 - 1/a1) / sqrt(e1) ≈ 1.4451
const EXPECTED: &str = "\
CHROM\tBIN_START\tN_SNPS\tTajimaD\n\
chr1\t0\t2\t1.75324\n\
chr1\t1000\t1\t-0.933021\n\
chr2\t0\t1\t1.4451\n\
";

#[test]
fn tajima_d_matches_expected() {
    let path = write_vcf(TEST_VCF);
    let rows = compute_tajima_d(&path, 1000).unwrap();
    let mut got = header().to_string();
    for row in &rows {
        got.push_str(&row.to_text());
    }
    assert_eq!(got, EXPECTED, "Tajima's D output differs from expected");
}

// ── chr2 site manual verification ───────────────────────────────────────────

#[test]
fn chr2_site_d_value() {
    // Verify the chr2 bin0 D value separately.
    // S=1, pi=0.6, n=6 haplotypes; vcftools reports 1.4451
    let d = rsomics_vcf_tajima_d::tajima_d(0.6, 1, 6);
    assert!(
        (d - 1.4451).abs() < 1e-4,
        "chr2 D: got {d}, expected 1.4451"
    );
}

// ── Live vcftools oracle ─────────────────────────────────────────────────────

fn vcftools_version() -> Option<String> {
    let out = Command::new("vcftools").arg("--version").output().ok()?;
    let combined =
        String::from_utf8_lossy(&out.stdout).to_string() + &String::from_utf8_lossy(&out.stderr);
    combined.lines().next().map(str::to_string)
}

fn skip_unless_vcftools_017() -> Option<()> {
    vcftools_version()?.contains("0.1.17").then_some(())
}

fn oracle_tajima_d(vcf: &std::path::Path, window_size: u64) -> Option<String> {
    let ts = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.subsec_nanos())
        .unwrap_or(0);
    let prefix = std::env::temp_dir().join(format!("rsomics_oracle_td_{ts}"));
    let status = Command::new("vcftools")
        .args([
            "--vcf",
            vcf.to_str()?,
            "--TajimaD",
            &window_size.to_string(),
            "--out",
            prefix.to_str()?,
        ])
        .status()
        .ok()?;
    if !status.success() {
        return None;
    }
    let out_path = prefix.with_extension("Tajima.D");
    std::fs::read_to_string(out_path).ok()
}

#[test]
fn live_oracle_tajima_d() {
    if skip_unless_vcftools_017().is_none() {
        eprintln!("vcftools 0.1.17 not found — skipping live oracle test");
        return;
    }
    let path = write_vcf(TEST_VCF);
    let oracle = oracle_tajima_d(&path, 1000).expect("vcftools --TajimaD failed");
    let rows = compute_tajima_d(&path, 1000).unwrap();
    let mut got = header().to_string();
    for row in &rows {
        got.push_str(&row.to_text());
    }
    assert_eq!(got, oracle, "Tajima's D differs from vcftools oracle");
}
