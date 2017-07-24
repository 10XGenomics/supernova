//
// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
//

use utils;
use std::path::Path;
use utils::MultiVec;

use std::fs::File;
use std::io::{BufWriter, Write};

use kmer;
use kmer::dir::Dir;
use debruijn::{self, VEdge, TempGraph};
use csv;


#[derive(RustcEncodable)]
struct GraphNode {
    id: usize,
    len: usize,
    num_bcs: usize,
    exts_left: u8,
    exts_right: u8,
    sequence: String,
}


pub fn main_graph_stats(graph: &Path, bcs: &Path, stats: &Path)
{
    let graph : TempGraph = utils::read_obj(graph).expect("read graph");
    let bcs : MultiVec<u32> = utils::read_obj(bcs).expect("read bcs");

    println!("bcs: start_pos: {:?}, items: {:?}, vec_len: {:?}", &bcs.start_pos[0..10], &bcs.items[0..40], &bcs.vec_len[0..10]);

    let mut wtr = csv::Writer::from_file(stats).expect("open csv").delimiter(b'\t');

    for e in graph.iter()
    {
        let record = GraphNode {
            id: e.id,
            len: e.sequence.length,
            num_bcs: bcs.vec_len[e.id] as usize,
            exts_left: e.exts.num_exts_l(),
            exts_right: e.exts.num_exts_r(),
            sequence: e.sequence.to_string(),
        };

        wtr.encode(record).expect("csv write error");
    }
}


pub fn write_vedge_to_gfa(w: &mut Write, e: VEdge) {
    writeln!(w, "S\t{}\t{}", e.id, e.sequence.to_dna_string()).unwrap();

    for (target, dir, _) in e.l_exts {
        if target > e.id as usize {
            let to_dir = match dir { Dir::Left => "+", Dir::Right => "-" };
            writeln!(w, "L\t{}\t{}\t{}\t{}\t{}M", e.id, "-", target, to_dir, kmer::K-1).unwrap();
        }
    }

    for (target, dir, _) in e.r_exts {
        if target > e.id as usize {
            let to_dir = match dir { Dir::Left => "+", Dir::Right => "-" };
            writeln!(w, "L\t{}\t{}\t{}\t{}\t{}M", e.id, "+", target, to_dir, kmer::K-1).unwrap();
        }
    }
}

pub fn main_write_gfa(graph: &Path, gfa_out: &Path)
{
    let graph : TempGraph = utils::read_obj(graph).expect("read graph");
    let edge_db = debruijn::EdgeDb::new(&graph);

    let mut wtr = File::create(gfa_out).unwrap();

    writeln!(wtr, "H\tVN:Z:tadaV1").unwrap();

    for e in 0 .. graph.start_pos.len() {
        let e = debruijn::make_vedge(&graph, &edge_db, e);
        write_vedge_to_gfa(&mut wtr, e);
    }
}



pub fn write_graph_bcs_matrix(_graph: &Path, _bcs: &Path, mm_file: &Path) {
    let graph : TempGraph = utils::read_obj(_graph).expect("read graph");
    let bcs : MultiVec<u32> = utils::read_obj(_bcs).expect("read bcs");

    let max_bc_id = *bcs.items.iter().max().unwrap();
    let num_edges = graph.start_pos.len();

    let mut _wtr = File::create(mm_file).unwrap();
    let mut wtr = BufWriter::new(_wtr);

    // Main header
    writeln!(wtr, "%%MatrixMarket {} {} {} {}", "matrix", "coordinate", "pattern", "general").unwrap();
    
    // Comment
    writeln!(wtr, "% edge -> bc matrix for graph: {:?}", _graph).unwrap();

    // size line
    writeln!(wtr, "{} {} {}", num_edges, max_bc_id+1, bcs.items.len()).unwrap();

    for edge_id in 0 .. graph.start_pos.len() {
        let bc_slice = bcs.get_slice(edge_id);

        for bc_id in bc_slice {
            writeln!(wtr, "{} {}", edge_id + 1, bc_id + 1).unwrap();
        }
    }
}

pub fn write_graph_fasta(_graph: &Path, fasta: &Path) {

    let graph: TempGraph = utils::read_obj(_graph).expect("can't read graph");
    let edge_db = debruijn::EdgeDb::new(&graph);

    let mut _wtr = File::create(fasta).unwrap();
    let mut wtr = BufWriter::new(_wtr);

    for edge_id in 0 .. graph.start_pos.len() {
        let e = debruijn::make_vedge(&graph, &edge_db, edge_id);
        writeln!(wtr, ">{}", e.id).unwrap();
        writeln!(wtr, "{}", e.sequence.to_dna_string()).unwrap();
    }
}


pub fn write_basevector_graph(_graph: &Path, bv_path: &Path) {

    let graph: TempGraph = utils::read_obj(_graph).expect("can't read graph");

    let mut _wtr = File::create(bv_path).unwrap();
    let mut wtr = BufWriter::new(_wtr);
    graph.write_to_sn_format(&mut wtr);
}
