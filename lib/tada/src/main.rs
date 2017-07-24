//
// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
//

//! # tada, assembler components in Rust
//! tada implements a front-end DeBruijn graph construction using the MSP strategy.

#![crate_name = "tada"]
#![allow(dead_code)]

extern crate byteorder;
extern crate bincode;
extern crate rustc_serialize;
extern crate bio;
extern crate flate2;
extern crate docopt;
extern crate rand;
extern crate chrono;
extern crate time;
extern crate itertools;
extern crate linked_hash_map;
extern crate csv;
extern crate martian;
extern crate rayon;
extern crate crossbeam;
extern crate libc;
extern crate bit_set;
extern crate shardio;
extern crate tempfile;
extern crate fastq_10x;
extern crate nix;
extern crate backtrace;
extern crate lz4;
extern crate smallvec;
extern crate rustc_version_runtime;

#[macro_use]
extern crate log;
extern crate env_logger;

#[cfg(test)]
#[macro_use]
extern crate quickcheck;

pub mod cmd_msp;
pub mod cmd_shard_asm;
pub mod cmd_main_asm;
pub mod cmd_graph_stats;

pub mod kmer;
pub mod msp;
pub mod debruijn;
pub mod fx;
pub mod bitenc;
pub mod bwt;
pub mod load_graph;
pub mod external;
pub mod multifastq;
pub mod utils;
pub mod sim_tests;
pub mod min_hash;
pub mod scaf_graph;
pub mod cmd_sort_fastq;

use docopt::Docopt;
use std::path::{Path, PathBuf};
use std::env;
use log::{LogRecord, LogLevelFilter};
use env_logger::LogBuilder;
use std::collections::HashMap;
use chrono::*;
use martian::MartianStage;
use nix::sys::signal;


const USAGE: &'static str = "
Top-down assembler

Usage:
  tada exp <graph> [<edge-db>]
  tada msp <whitelist> <permutation> <output> <files>...
  tada shard-asm [--min-kmer-obs=<mk>] <output> <files>...
  tada stats <graph> <bcs> <stats>
  tada fasta <graph> <fasta>
  tada gfa <graph> <gfa>
  tada bv <graph> <output>
  tada bcmat <graph> <bcs> <mm-file>
  tada scaf-graph <graph> <bcs> <mm-file> <min-ctg> <min-bcs> <max-bcs>
  tada martian <adapter>...
  tada ingest [--num-threads=<n>] <whitelist> <chunk-info> <output>
  tada (-h | --help)
  tada --version

Options:
  --num-threads=<n>    Number of threads to use. By default use one thread per CPU [default: 0]
  --min-kmer-obs=<mk>  Minimum number of kmer observations required [default: 3]
  -h --help            Show this screen.
  --version            Show version.
";

#[derive(Debug, RustcDecodable)]
struct Args {
    arg_whitelist: String,
    arg_permutation: String,
    arg_output: String,
    arg_files: Vec<String>,
    arg_edge_db: Option<String>,

    flag_min_kmer_obs: u32,
    flag_num_threads: usize,
    flag_version: bool,

    arg_graph: String,
    arg_bcs: String,
    arg_stats: String,
    arg_gfa: String,
    arg_mm_file: String,
    arg_fasta: String,
    arg_min_bcs: usize,
    arg_min_ctg: usize,
    arg_max_bcs: usize,
    arg_chunk_info: String,

    cmd_exp: bool,
    cmd_msp: bool,
    cmd_shard_asm: bool,
    cmd_stats: bool,
    cmd_gfa: bool,
    cmd_bcmat: bool,
    cmd_fasta: bool,
    cmd_scaf_graph: bool,
    cmd_bv: bool,
    cmd_ingest: bool,

    // Martian interface mode
    cmd_martian: bool,
    arg_adapter: Vec<String>,
}

// Load git build version from env
fn get_version() -> String {
    let v = include_str!(concat!(env!("OUT_DIR"), "/version.txt"));
    let mut version = String::from(v);
    let nl = version.len() - 1;
    version.truncate(nl);
    version
}

fn main() {

    let args: Args = Docopt::new(USAGE)
                         .and_then(|d| d.decode())
                         .unwrap_or_else(|e| e.exit());
    println!("{:?}", args);
    info!("{:?}", args);
    setup_sig();

    if args.flag_version {
        println!("Git version of repo: {}", get_version());
        println!("rustc version {}", rustc_version_runtime::version());
        return;
    }

    // Martian does it's own logging -- don't configure if we're being called by martian
    if !args.cmd_martian {

        // setup logging
        let format = |record: &LogRecord| {
            let time = Local::now().format("%Y-%m-%d %H:%M:%S").to_string();
            format!("[{}] {} - {}", time, record.level(), record.args())
        };

        let mut builder = LogBuilder::new();
        builder.format(format).filter(None, LogLevelFilter::Info);

        if env::var("RUST_LOG").is_ok() {
           builder.parse(&env::var("RUST_LOG").unwrap());
        }

        builder.init().unwrap();
    }

    // Log version info
    info!("Git version of repo: {}", get_version());
    info!("rustc version {}", rustc_version_runtime::version());

    if args.cmd_msp {
        let mut paths: Vec<PathBuf> = Vec::new();
        for p in args.arg_files {
            paths.push(PathBuf::from(p));
        }

        cmd_msp::main_msp(
            7,
            Path::new(&args.arg_whitelist), paths,
            Path::new(&args.arg_permutation),
            Path::new(&args.arg_output));

    } else if args.cmd_shard_asm {

        let mut paths: Vec<PathBuf> = Vec::new();
        for p in args.arg_files {
            paths.push(PathBuf::from(p));
        }
        cmd_shard_asm::main_shard_asm(args.flag_min_kmer_obs, 0, 1, paths, Path::new(&args.arg_output), Path::new("FIXME"), 8);

    } else if args.cmd_martian {

        let mut stage_registry: HashMap<String, Box<MartianStage>> = HashMap::new();
        stage_registry.insert("msp".to_string(), Box::new(cmd_msp::MspMartian));
        stage_registry.insert("shard-asm".to_string(), Box::new(cmd_shard_asm::ShardAsmMartian));
        stage_registry.insert("main-asm".to_string(), Box::new(cmd_main_asm::MainAsmMartian));
        stage_registry.insert("main-asm-sn".to_string(), Box::new(cmd_main_asm::MainAsmMartianSupernova));
        stage_registry.insert("bucket-bcs".to_string(), Box::new(cmd_sort_fastq::BucketBcsMartian));
        stage_registry.insert("sort-bcs".to_string(), Box::new(cmd_sort_fastq::SortBcsMartian));

        // Run the built-in martian adapter
        martian::martian_main(args.arg_adapter, stage_registry);

    } else if args.cmd_exp {
        debruijn::do_explore(args.arg_graph, args.arg_edge_db);
    } else if args.cmd_stats {
        cmd_graph_stats::main_graph_stats(
            Path::new(&args.arg_graph),
            Path::new(&args.arg_bcs),
            Path::new(&args.arg_stats),
        );
    } else if args.cmd_gfa {
        cmd_graph_stats::main_write_gfa(
            Path::new(&args.arg_graph),
            Path::new(&args.arg_gfa),
        );
    } else if args.cmd_bcmat {
        cmd_graph_stats::write_graph_bcs_matrix(
            Path::new(&args.arg_graph),
            Path::new(&args.arg_bcs),
            Path::new(&args.arg_mm_file),
        );
    } else if args.cmd_fasta {
        cmd_graph_stats::write_graph_fasta(
            Path::new(&args.arg_graph),
            Path::new(&args.arg_fasta),
        );
    } else if args.cmd_scaf_graph {
            scaf_graph::write_scaf_graph(
                Path::new(&args.arg_graph),
                Path::new(&args.arg_fasta),
                Path::new(&args.arg_mm_file),
                args.arg_min_ctg, args.arg_min_bcs, args.arg_max_bcs);
    } else if args.cmd_bv {
        cmd_graph_stats::write_basevector_graph(
            Path::new(&args.arg_graph),
            Path::new(&args.arg_output),
        )
    } else if args.cmd_ingest {
        let bc_counts = args.arg_output.clone() + ".bcs";
        
        cmd_sort_fastq::go(
            Path::new(&args.arg_whitelist),
            &args.arg_chunk_info,
            Path::new(&args.arg_output),
            Path::new(&bc_counts),
            args.flag_num_threads,
        );
    }
}


extern fn handle_sigint(_:i32) {
        println!("caught signal - Interrupted!");
            let bt = backtrace::Backtrace::new();
                println!("bt:\n{:?}", bt);
}

fn setup_sig() {
      let sig_action = signal::SigAction::new(signal::SigHandler::Handler(handle_sigint),
                                              signal::SaFlags::empty(),
                                              signal::SigSet::empty());
      unsafe { signal::sigaction(signal::SIGQUIT, &sig_action).unwrap(); }
}
