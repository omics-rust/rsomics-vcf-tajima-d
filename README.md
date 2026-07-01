# rsomics-vcf-tajima-d

Tajima's D statistic in fixed-width windows from a VCF file.

Output is byte-identical to `vcftools --TajimaD`.

## Install

```
cargo install rsomics-vcf-tajima-d
```

## Usage

```
rsomics-vcf-tajima-d [OPTIONS] <VCF>

Arguments:
  <VCF>  Input VCF (plain or .gz)

Options:
      --window-size <BP>  Window size in bp [default: 10000]
      --json              Emit JSON envelope instead of tab-delimited text
  -h, --help              Print help
  -V, --version           Print version
```

Output columns: `CHROM`, `BIN_START`, `N_SNPS`, `TajimaD` — identical to vcftools.

## Performance

On a 100k-variant synthetic VCF (3 chromosomes, 3 samples, window 10 kbp):

| Tool | Time (mean ± σ) |
|------|----------------|
| rsomics-vcf-tajima-d 0.1.0 | 43.7 ms ± 3.0 ms |
| vcftools 0.1.17 --TajimaD 10000 | 131.0 ms ± 6.0 ms |

**3.00× faster** (mini_m2 aarch64-apple-darwin, 1 thread, hyperfine 10 runs).

## Origin

This crate is an independent Rust reimplementation based on:
- Tajima F (1989). *Statistical method for testing the neutral mutation hypothesis by DNA polymorphism.* Genetics 123(3):585–595.
- vcftools 0.1.17 black-box behavior (Danecek et al. 2011, Bioinformatics 27(15):2156–2158).

The Tajima's D variance coefficients (e₁, e₂) follow Tajima (1989) Table 2.
BIN_START uses 0-based coordinates aligned to multiples of the window size.

License: MIT OR Apache-2.0.
Upstream credit: vcftools <https://vcftools.github.io/> (GPL-3.0).
