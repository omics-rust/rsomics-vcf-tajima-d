use criterion::{Criterion, criterion_group, criterion_main};
use std::path::PathBuf;

fn bench_tajima_d(c: &mut Criterion) {
    let fixture = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../rsomics-fixtures/vcf-tajima-d/bench_100k.vcf");
    if !fixture.exists() {
        eprintln!("fixture not found: {}", fixture.display());
        return;
    }
    c.bench_function("tajima_d_100k", |b| {
        b.iter(|| rsomics_vcf_tajima_d::compute_tajima_d(&fixture, 10000).unwrap());
    });
}

criterion_group!(benches, bench_tajima_d);
criterion_main!(benches);
