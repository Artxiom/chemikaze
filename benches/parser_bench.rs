use criterion::{
    BenchmarkId, Criterion, SamplingMode, Throughput, black_box, criterion_group, criterion_main,
};
use std::fs;

// Import the modules we need to benchmark
mod atom_counts {
    include!("../src/main/rust/atom_counts.rs");
}
mod errors {
    include!("../src/main/rust/errors.rs");
}
mod periodic_table {
    include!("../src/main/rust/periodic_table.rs");
}
mod mf_parser {
    include!("../src/main/rust/mf_parser.rs");
}

use mf_parser::MfParser;

fn bench_single_mf_parsing(c: &mut Criterion) {
    let test_cases = vec![
        ("simple", "H2O"),
        ("medium", "C67H132N8O3"),
        ("complex", "[(2H2O.NaCl)3S.N]2-"),
        ("parentheses", "(C(OH)2)2P"),
        ("leading_coeff", "2H2O"),
        ("dots", "NH3.2CH3"),
    ];

    let mut group = c.benchmark_group("single_mf_parsing");

    for (name, mf) in test_cases {
        group.bench_with_input(BenchmarkId::new("new_parser_each", name), &mf, |b, mf| {
            b.iter(|| black_box(MfParser::parse_single(black_box(*mf)).unwrap()));
        });

        group.bench_with_input(BenchmarkId::new("reusable_parser", name), &mf, |b, mf| {
            let mut parser = MfParser::new();
            b.iter(|| black_box(parser.parse(black_box(*mf)).unwrap()))
        });
    }
    group.finish();
}

fn bench_full_dataset(c: &mut Criterion) {
    // Load the full test dataset
    let content =
        fs::read_to_string("src/test/resources/MFs.csv").expect("Could not read test file");
    let lines: Vec<&str> = content.lines().collect();
    let repeat = 50usize;

    let mut group = c.benchmark_group("full_dataset");
    group.sample_size(10); // Fewer samples since this is expensive
    group.throughput(Throughput::Elements((lines.len() * repeat) as u64));
    group.sampling_mode(SamplingMode::Flat);

    group.bench_function("reusable_parser_full", |b| {
        let mut parser = MfParser::new();
        b.iter(|| {
            for _ in 0..repeat {
                for &mf in &lines {
                    black_box(parser.parse(black_box(mf)).unwrap());
                }
            }
        })
    });

    group.finish();
}

criterion_group!(benches, bench_single_mf_parsing, bench_full_dataset,);
criterion_main!(benches);
