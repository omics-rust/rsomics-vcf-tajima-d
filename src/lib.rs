//! Tajima's D statistic in fixed-width windows from VCF.
//!
//! Implements vcftools `--TajimaD` exactly:
//! - Non-overlapping windows of `window_size` bp per chromosome.
//! - BIN_START = floor(pos / window_size) × window_size (0-based; sites at exact multiples of W go into the higher bin).
//! - Tajima's π per site: `(2n/(2n−1)) × 2p(1−p)` (only polymorphic sites).
//! - Watterson's θ_W = S / a₁ where a₁ = Σ(1/i) for i=1..2n−1, 2n = haplotype count.
//! - D = (π_sum − θ_W) / sqrt(e₁·S + e₂·S·(S−1))
//!   using Tajima (1989) variance coefficients e₁, e₂.
//!
//! Windows with fewer than 2 segregating sites return NaN for D (matching vcftools).

use std::collections::HashMap;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

use flate2::read::MultiGzDecoder;
use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

// ── Output row ───────────────────────────────────────────────────────────────

#[derive(Debug, Clone, Serialize)]
pub struct TajimaDRow {
    pub chrom: String,
    pub bin_start: u64,
    pub n_snps: u64,
    pub tajima_d: f64,
}

impl TajimaDRow {
    /// Render as the tab-separated line vcftools emits.
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

/// Format like vcftools `%g` (6 significant figures, no trailing zeros/dot).
pub fn format_g(x: f64) -> String {
    if x.is_nan() {
        return "nan".to_string();
    }
    if x.is_infinite() {
        return if x > 0.0 {
            "inf".to_string()
        } else {
            "-inf".to_string()
        };
    }
    if x == 0.0 {
        return "0".to_string();
    }
    let mag = x.abs().log10().floor() as i32;
    let dec = (5 - mag).max(0) as usize;
    let s = format!("{x:.dec$}");
    let s = s.trim_end_matches('0').trim_end_matches('.');
    s.to_string()
}

// ── Tajima coefficients for a given haplotype count ─────────────────────────

/// Precomputed Tajima (1989) variance coefficients for a given 2n.
struct TajimaCoeffs {
    e1: f64,
    e2: f64,
}

impl TajimaCoeffs {
    fn new(two_n: u64) -> Self {
        let n = two_n as f64;
        let a1: f64 = (1..two_n).map(|i| 1.0 / i as f64).sum();
        let a2: f64 = (1..two_n).map(|i| 1.0 / (i as f64).powi(2)).sum();
        let b1 = (n + 1.0) / (3.0 * (n - 1.0));
        let b2 = 2.0 * (n * n + n + 3.0) / (9.0 * n * (n - 1.0));
        let c1 = b1 - 1.0 / a1;
        let c2 = b2 - (n + 2.0) / (a1 * n) + a2 / (a1 * a1);
        let e1 = c1 / a1;
        let e2 = c2 / (a1 * a1 + a2);
        Self { e1, e2 }
    }
}

/// Compute Tajima's D from π_sum, S segregating sites, and 2n haplotypes.
///
/// Returns NaN when variance is zero (S=0 or S=1 with e2·0=0 → e1·S).
pub fn tajima_d(pi_sum: f64, s: u64, two_n: u64) -> f64 {
    if s == 0 || two_n < 2 {
        return f64::NAN;
    }
    let a1: f64 = (1..two_n).map(|i| 1.0 / i as f64).sum();
    let theta_w = s as f64 / a1;
    let coeffs = TajimaCoeffs::new(two_n);
    let v = coeffs.e1 * s as f64 + coeffs.e2 * s as f64 * (s as f64 - 1.0);
    if v <= 0.0 {
        return f64::NAN;
    }
    (pi_sum - theta_w) / v.sqrt()
}

// ── Window accumulator ───────────────────────────────────────────────────────

/// Per-bin accumulator: (pi_sum, segregating_sites, haplotype_count_of_first_site).
/// vcftools uses the haplotype count from the first site in each window.
struct BinAcc {
    pi_sum: f64,
    n_snps: u64,
    two_n: u64,
}

struct WindowAcc {
    window_size: u64,
    cur_chrom: String,
    rows: Vec<TajimaDRow>,
    bins: HashMap<u64, BinAcc>,
    max_bin: Option<u64>,
}

impl WindowAcc {
    fn new(window_size: u64) -> Self {
        Self {
            window_size,
            cur_chrom: String::new(),
            rows: Vec::new(),
            bins: HashMap::new(),
            max_bin: None,
        }
    }

    fn bin_of(&self, pos: u64) -> u64 {
        // vcftools: BIN_START = floor(pos/W)*W  (0-indexed bins at multiples of W)
        (pos / self.window_size) * self.window_size
    }

    fn flush_chrom(&mut self) {
        if self.cur_chrom.is_empty() {
            return;
        }
        let Some(max_bin) = self.max_bin else {
            return;
        };
        // Emit only bins that had at least one SNP.
        let mut bin_keys: Vec<u64> = self.bins.keys().copied().collect();
        bin_keys.sort_unstable();
        for k in bin_keys {
            let acc = &self.bins[&k];
            let d = tajima_d(acc.pi_sum, acc.n_snps, acc.two_n);
            self.rows.push(TajimaDRow {
                chrom: self.cur_chrom.clone(),
                bin_start: k,
                n_snps: acc.n_snps,
                tajima_d: d,
            });
        }
        let _ = max_bin;
        self.bins.clear();
        self.max_bin = None;
    }

    fn push_site(&mut self, chrom: &str, pos: u64, pi_site: f64, two_n: u64) {
        if chrom != self.cur_chrom {
            self.flush_chrom();
            self.cur_chrom = chrom.to_string();
        }
        let k = self.bin_of(pos);
        let e = self.bins.entry(k).or_insert(BinAcc {
            pi_sum: 0.0,
            n_snps: 0,
            two_n,
        });
        e.pi_sum += pi_site;
        e.n_snps += 1;
        self.max_bin = Some(self.max_bin.map_or(k, |m| m.max(k)));
    }

    fn finish(mut self) -> Vec<TajimaDRow> {
        self.flush_chrom();
        self.rows
    }
}

// ── Per-site π and haplotype count ───────────────────────────────────────────

/// Compute per-site π (Tajima's estimator) and haplotype count from genotype fields.
/// Returns `None` when fewer than 2 called haplotypes.
/// Returns `Some((0.0, two_n))` for monomorphic (all-REF or all-ALT) sites.
pub fn site_pi_with_n(gt_fields: &[&str]) -> Option<(f64, u64)> {
    let mut n_ref: u64 = 0;
    let mut n_alt: u64 = 0;
    for field in gt_fields {
        let gt = if let Some(c) = field.find(':') {
            &field[..c]
        } else {
            field
        };
        for a in gt.split(['/', '|']) {
            match a {
                "0" => n_ref += 1,
                "." => {}
                _ if a.parse::<u64>().is_ok_and(|v| v > 0) => n_alt += 1,
                _ => {}
            }
        }
    }
    let two_n = n_ref + n_alt;
    if two_n < 2 {
        return None;
    }
    let p = n_alt as f64 / two_n as f64;
    let q = 1.0 - p;
    let pi = (two_n as f64 / (two_n - 1) as f64) * 2.0 * p * q;
    Some((pi, two_n))
}

// ── VCF scanner ─────────────────────────────────────────────────────────────

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

pub fn compute_tajima_d(path: &Path, window_size: u64) -> Result<Vec<TajimaDRow>> {
    let reader = open_reader(path)?;
    let mut lines = BufReader::new(reader).lines();
    let mut acc = WindowAcc::new(window_size);
    let mut found_chrom = false;

    for line in lines.by_ref() {
        let line = line?;
        let line = line.trim_end_matches('\r');
        if line.starts_with("##") {
            continue;
        }
        if line.starts_with('#') {
            found_chrom = true;
            continue;
        }
        if !found_chrom {
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
        let gt_fields = &cols[FIRST_SAMPLE..];
        if let Some((pi, two_n)) = site_pi_with_n(gt_fields) {
            // Only polymorphic sites (pi > 0) contribute to N_SNPS.
            if pi > 0.0 {
                acc.push_site(chrom, pos, pi, two_n);
            }
        }
    }

    if !found_chrom {
        return Err(RsomicsError::InvalidInput(
            "VCF missing #CHROM header line".into(),
        ));
    }
    Ok(acc.finish())
}

pub fn header() -> &'static str {
    "CHROM\tBIN_START\tN_SNPS\tTajimaD\n"
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn format_g_zero() {
        assert_eq!(format_g(0.0), "0");
    }

    #[test]
    fn format_g_nan() {
        assert_eq!(format_g(f64::NAN), "nan");
    }

    #[test]
    fn tajima_d_window1() {
        // n=3 samples → 2n=6; S=2; pi_sum=0.6+0.6=1.2
        let d = tajima_d(1.2, 2, 6);
        // Expected: 1.75324
        assert!((d - 1.75324).abs() < 1e-5, "d={d}");
    }

    #[test]
    fn tajima_d_window2() {
        // n=3 samples → 2n=6; S=1; pi = (6/5)*2*(5/6)*(1/6) = 1/3
        let pi = (6.0_f64 / 5.0) * 2.0 * (5.0 / 6.0) * (1.0 / 6.0);
        let d = tajima_d(pi, 1, 6);
        // Expected: -0.933021
        assert!((d - (-0.933021)).abs() < 1e-6, "d={d}");
    }

    #[test]
    fn tajima_d_zero_snps_is_nan() {
        let d = tajima_d(0.0, 0, 6);
        assert!(d.is_nan());
    }
}
