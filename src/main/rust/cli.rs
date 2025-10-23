use std::time::*;
use std::{env, fs};

mod atom_counts;
mod errors;
mod mf_parser;
mod periodic_table;

use mf_parser::MfParser;

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        panic!("[ERROR] Pass filename as a parameter!");
    }
    let filepath = &args[1];
    let repeats = 50;
    let content = fs::read_to_string(filepath).expect(&format!("Couldn't read {filepath}"));
    let lines: Vec<&[u8]> = content.split("\n").map(|l| l.as_bytes()).collect();
    let mf_cnt = repeats * lines.len();
    let mut parser = MfParser::new();

    let start = Instant::now();
    parse_mfs(&mut parser, &lines, repeats);
    let elapsed = start.elapsed();
    println!(
        "[RUST BENCHMARK] {mf_cnt} MFs in {:.2?} ({}MF/s)",
        elapsed,
        (mf_cnt as f64 / elapsed.as_secs_f64()) as u32
    );
}

fn parse_mfs(parser: &mut MfParser, mfs: &Vec<&[u8]>, n: usize) -> u32 {
    let mut hcount: u32 = 0;
    for _ in 0..n {
        for mf in mfs {
            hcount += parser.parse(*mf).unwrap().counts[0];
        }
    }
    hcount // return something so that this isn't optimized out
}
