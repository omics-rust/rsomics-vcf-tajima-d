//! Tajima's D statistic in fixed-width windows from VCF.
//!
//! Reproduces vcftools `--TajimaD` exactly, including its whole-file constant
//! sample size:
//! - `n` = 2 × (number of samples), fixed for the entire file. Every variance
//!   coefficient (a1, a2, b1, b2, c1, c2, e1, e2) is derived from this constant
//!   n, never from a per-site called-chromosome count.
//! - Only biallelic (single concrete ALT), fully diploid sites participate. A
//!   site with any haploid genotype is skipped; a polyploid genotype aborts.
//! - Per site: `p = ref_count / non_missing_chr`; a site contributes to a bin
//!   only when `0 < p < 1`, adding `p(1−p)` and incrementing the SNP count.
//! - Per window: `π = 2·Σp(1−p)·n/(n−1)`, `θ_W = S/a1`,
//!   `D = (π − θ_W) / sqrt(e1·S + e2·S·(S−1))`.
//!
//! Binning: `BIN_START = floor(pos/W)·W`, computed as `(pos·(1/W)) as u64` to
//! match vcftools' floating-point truncation. Once the first SNP-bearing bin of
//! a chromosome is emitted, every subsequent bin up to the last biallelic-diploid
//! site is emitted too — empty windows in that span carry `N_SNPS=0` and `nan`.

use std::collections::HashMap;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

use flate2::read::MultiGzDecoder;
use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

#[derive(Debug, Clone, Serialize)]
pub struct TajimaDRow {
    pub chrom: String,
    pub bin_start: u64,
    pub n_snps: u64,
    pub tajima_d: f64,
}

impl TajimaDRow {
    pub fn to_text(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\n",
            self.chrom,
            self.bin_start,
            self.n_snps,
            format_g(self.tajima_d),
        )
    }
}

/// Render `x` exactly as C `printf("%g", x)` (default precision 6), matching the
/// C++ ostream default float format vcftools writes with.
pub fn format_g(x: f64) -> String {
    if x.is_nan() {
        return "nan".to_string();
    }
    if x.is_infinite() {
        return if x < 0.0 { "-inf" } else { "inf" }.to_string();
    }
    format_g_prec(x, 6)
}

fn format_g_prec(x: f64, prec: usize) -> String {
    if x == 0.0 {
        return if x.is_sign_negative() { "-0" } else { "0" }.to_string();
    }

    let sig = prec.max(1);
    let sci = format!("{:.*e}", sig - 1, x);
    let epos = sci.find('e').unwrap();
    let exp: i32 = sci[epos + 1..].parse().unwrap();

    if exp < -4 || exp >= sig as i32 {
        let mantissa = sci[..epos].trim_end_matches('0').trim_end_matches('.');
        format!(
            "{mantissa}e{}{:02}",
            if exp < 0 { '-' } else { '+' },
            exp.abs()
        )
    } else {
        let decimals = (sig as i32 - 1 - exp).max(0) as usize;
        let s = format!("{x:.decimals$}");
        if s.contains('.') {
            s.trim_end_matches('0').trim_end_matches('.').to_string()
        } else {
            s
        }
    }
}

/// Whole-file Tajima variance coefficients, derived once from the constant
/// chromosome count `n = 2 × samples`.
struct Coeffs {
    a1: f64,
    e1: f64,
    e2: f64,
    n: f64,
}

impl Coeffs {
    fn new(n: u64) -> Self {
        let nf = n as f64;
        let a1: f64 = (1..n).map(|i| 1.0 / i as f64).sum();
        let a2: f64 = (1..n).map(|i| 1.0 / (i as f64 * i as f64)).sum();
        let b1 = (nf + 1.0) / 3.0 / (nf - 1.0);
        let b2 = 2.0 * (nf * nf + nf + 3.0) / 9.0 / nf / (nf - 1.0);
        let c1 = b1 - 1.0 / a1;
        let c2 = b2 - (nf + 2.0) / (a1 * nf) + a2 / a1 / a1;
        let e1 = c1 / a1;
        let e2 = c2 / (a1 * a1 + a2);
        Self { a1, e1, e2, n: nf }
    }

    fn tajima_d(&self, sum_pq: f64, s: u64) -> f64 {
        if s == 0 {
            return f64::NAN;
        }
        let sf = s as f64;
        let pi = 2.0 * sum_pq * self.n / (self.n - 1.0);
        let tw = sf / self.a1;
        let var = self.e1 * sf + self.e2 * sf * (sf - 1.0);
        (pi - tw) / var.sqrt()
    }
}

/// Standalone Tajima's D for the given accumulated Σp(1−p), segregating-site
/// count, and constant chromosome count `n`.
pub fn tajima_d(sum_pq: f64, s: u64, n: u64) -> f64 {
    if n < 2 {
        return f64::NAN;
    }
    Coeffs::new(n).tajima_d(sum_pq, s)
}

/// Outcome of scanning one biallelic site's genotypes.
enum SiteScan {
    /// Site is not fully diploid (a haploid genotype present) — skip entirely.
    NotDiploid,
    /// Site is diploid: reference-allele count and non-missing chromosome count.
    Diploid { ref_count: u64, non_missing: u64 },
    /// A polyploid genotype — vcftools aborts here.
    Polyploid,
}

/// Locate the GT sub-field index within a FORMAT string.
fn gt_index(format: &str) -> usize {
    format.split(':').position(|k| k == "GT").unwrap_or(0)
}

/// Scan a site's genotype fields, mirroring vcftools' allele accounting.
fn scan_site(gt_fields: &[&str], gt_idx: usize) -> SiteScan {
    let mut ref_count = 0u64;
    let mut non_missing = 0u64;
    for field in gt_fields {
        let gt = field.split(':').nth(gt_idx).unwrap_or("");
        let sep_count = gt.bytes().filter(|&b| b == b'/' || b == b'|').count();
        match sep_count {
            0 => return SiteScan::NotDiploid,
            1 => {}
            _ => return SiteScan::Polyploid,
        }
        for allele in gt.split(['/', '|']) {
            if allele == "." {
                continue;
            }
            non_missing += 1;
            if allele == "0" {
                ref_count += 1;
            }
        }
    }
    SiteScan::Diploid {
        ref_count,
        non_missing,
    }
}

/// Per-chromosome bins, indexed by window index; each is (S, Σp(1−p)).
struct Bins {
    chrom_order: Vec<String>,
    prev_chrom: Option<String>,
    map: HashMap<String, Vec<(u64, f64)>>,
}

impl Bins {
    fn new() -> Self {
        Self {
            chrom_order: Vec::new(),
            prev_chrom: None,
            map: HashMap::new(),
        }
    }

    fn touch(&mut self, chrom: &str, idx: usize) -> &mut (u64, f64) {
        let bins = self.map.entry(chrom.to_string()).or_default();
        if idx >= bins.len() {
            bins.resize(idx + 1, (0, 0.0));
        }
        if self.prev_chrom.as_deref() != Some(chrom) {
            self.chrom_order.push(chrom.to_string());
            self.prev_chrom = Some(chrom.to_string());
        }
        &mut bins[idx]
    }

    fn into_rows(self, coeffs: &Coeffs, window_size: u64) -> Vec<TajimaDRow> {
        let mut rows = Vec::new();
        for chrom in &self.chrom_order {
            let bins = &self.map[chrom];
            let mut started = false;
            for (s, (n_snps, sum_pq)) in bins.iter().enumerate() {
                if *n_snps > 0 {
                    started = true;
                }
                if started {
                    rows.push(TajimaDRow {
                        chrom: chrom.clone(),
                        bin_start: s as u64 * window_size,
                        n_snps: *n_snps,
                        tajima_d: coeffs.tajima_d(*sum_pq, *n_snps),
                    });
                }
            }
        }
        rows
    }
}

fn open_reader(path: &Path) -> Result<Box<dyn Read>> {
    let file = std::fs::File::open(path).map_err(|e| {
        RsomicsError::Io(std::io::Error::new(
            e.kind(),
            format!("cannot open {}: {e}", path.display()),
        ))
    })?;
    let is_gz = path
        .extension()
        .is_some_and(|e| e.eq_ignore_ascii_case("gz"));
    Ok(if is_gz {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    })
}

const FIRST_SAMPLE: usize = 9;
const COL_CHROM: usize = 0;
const COL_POS: usize = 1;
const COL_ALT: usize = 4;
const COL_FORMAT: usize = 8;

pub fn compute_tajima_d(path: &Path, window_size: u64) -> Result<Vec<TajimaDRow>> {
    let reader = open_reader(path)?;
    compute_tajima_d_reader(reader, window_size)
}

/// Core computation over any byte source (plain-text VCF). `compute_tajima_d`
/// wraps this after transparently opening plain or gzip files.
pub fn compute_tajima_d_reader<R: Read>(reader: R, window_size: u64) -> Result<Vec<TajimaDRow>> {
    if window_size == 0 {
        return Err(RsomicsError::InvalidInput(
            "window size must be positive".into(),
        ));
    }
    let lines = BufReader::new(reader).lines();

    let mut n_samples: Option<usize> = None;
    let mut coeffs: Option<Coeffs> = None;
    let mut bins = Bins::new();
    let inv_w = 1.0 / window_size as f64;

    for line in lines {
        let line = line?;
        let line = line.trim_end_matches('\r');
        if line.starts_with("##") {
            continue;
        }
        if line.starts_with('#') {
            let cols = line.split('\t').count();
            let samples = cols.saturating_sub(FIRST_SAMPLE);
            if samples * 2 < 2 {
                return Err(RsomicsError::InvalidInput(
                    "VCF has fewer than two chromosomes (need at least one sample)".into(),
                ));
            }
            n_samples = Some(samples);
            coeffs = Some(Coeffs::new(samples as u64 * 2));
            continue;
        }
        if n_samples.is_none() {
            return Err(RsomicsError::InvalidInput(
                "VCF missing #CHROM header line".into(),
            ));
        }

        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() <= FIRST_SAMPLE {
            continue;
        }
        let alt = cols[COL_ALT];
        if alt == "." || alt.contains(',') {
            continue;
        }
        let pos: u64 = match cols[COL_POS].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let chrom = cols[COL_CHROM];
        let gt_idx = gt_index(cols[COL_FORMAT]);
        let (ref_count, non_missing) = match scan_site(&cols[FIRST_SAMPLE..], gt_idx) {
            SiteScan::NotDiploid => continue,
            SiteScan::Polyploid => {
                return Err(RsomicsError::InvalidInput(format!(
                    "polyploidy found, not supported: {chrom}:{pos}"
                )));
            }
            SiteScan::Diploid {
                ref_count,
                non_missing,
            } => (ref_count, non_missing),
        };

        let idx = (pos as f64 * inv_w) as usize;
        let p = ref_count as f64 / non_missing as f64;
        let bin = bins.touch(chrom, idx);
        if p > 0.0 && p < 1.0 {
            bin.0 += 1;
            bin.1 += p * (1.0 - p);
        }
    }

    let coeffs = coeffs
        .ok_or_else(|| RsomicsError::InvalidInput("VCF missing #CHROM header line".into()))?;
    Ok(bins.into_rows(&coeffs, window_size))
}

pub fn header() -> &'static str {
    "CHROM\tBIN_START\tN_SNPS\tTajimaD\n"
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn format_g_matches_printf() {
        assert_eq!(format_g(0.0), "0");
        assert_eq!(format_g(f64::NAN), "nan");
        assert_eq!(format_g(f64::INFINITY), "inf");
        assert_eq!(format_g(f64::NEG_INFINITY), "-inf");
        assert_eq!(format_g(1.4450969118871448), "1.4451");
        assert_eq!(format_g(0.9417745826447361), "0.941775");
        assert_eq!(format_g(-1.51284), "-1.51284");
        assert_eq!(format_g(0.0709896), "0.0709896");
    }

    #[test]
    fn tajima_d_single_sample_is_nan() {
        // n = 2 → e1 = e2 = 0 → variance 0 → nan.
        assert!(tajima_d(0.25, 1, 2).is_nan());
    }

    #[test]
    fn tajima_d_zero_sites_is_nan() {
        assert!(tajima_d(0.0, 0, 6).is_nan());
    }
}
