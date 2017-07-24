//
// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
//

use utils::MultiVec;
use debruijn::TempGraph;
use min_hash::quick_jaccard;
use std::cmp::{min};
use itertools::Itertools;
use utils;
use std::io::BufWriter;
use std::fs::File;
use std::path::Path;
use std::io::Write;

pub fn estimate_distance(intersection: usize, union: usize, s1: usize, s2: usize, total_diversity: f32) -> f32 {

    let expected_ovl = s1 as f32 / total_diversity * s2 as f32;

    let exp_d = f32::max(1.0, intersection as f32 - expected_ovl) *  union as f32 / s1 as f32 / s2 as f32;
    -exp_d.ln()
}



pub fn find_ctg_overlaps(to_consider: &Vec<usize>, bcs: MultiVec<u32>, total_diversity: f32) -> Vec<(usize, usize, f32)> {

    let mut result = Vec::new();

    for i1 in 0..to_consider.len() {
        for i2 in (i1+1)..to_consider.len() {

            let slc1 = bcs.get_slice(to_consider[i1]);
            let slc2 = bcs.get_slice(to_consider[i2]);
            let (_, intersect, union) = quick_jaccard(slc1, slc2);

            let est_dist = estimate_distance(intersect, union, slc1.len(), slc2.len(), total_diversity);
            if est_dist < 2.0 {
                result.push((to_consider[i1], to_consider[i2], est_dist));
            }
        }
    }

    result
}

pub fn build_bc_scaffold_graph(
    g: TempGraph, bcs: MultiVec<u32>, 
    max_links: usize, 
    min_ctg: usize,
    max_bcs: usize,
    min_bcs: usize) -> Vec<(usize, usize, f32)> {

    // Pick longest 
    let mut to_consider = Vec::new();
    for e in g.iter() {
        if e.sequence.len() > min_ctg && bcs.get_slice(e.id).len() > min_bcs && bcs.get_slice(e.id).len() < max_bcs {
            to_consider.push(e.id);
        }
    }
    
    println!("nodes to consider: {}", to_consider.len());

    let overlaps = find_ctg_overlaps(&to_consider, bcs, 1.5e6);

    let mut best_overlaps = Vec::new();

    let mut sink_buf = Vec::new();
    for (_, sinks) in &overlaps.into_iter().group_by(|x| x.0) {
        sink_buf.clear();
        sink_buf.extend(sinks);
        sink_buf.sort_by(|a, b| a.2.partial_cmp(&b.2).unwrap());

        for i in 0 .. min(sink_buf.len(), max_links) {
            best_overlaps.push(sink_buf[i]);
        }
    }

    best_overlaps
}
    


pub fn write_scaf_graph(_graph: &Path, _bcs: &Path, mm_file: &Path, min_ctg: usize, min_bcs: usize, max_bcs: usize) {
    let graph : TempGraph = utils::read_obj(_graph).expect("read graph");
    let bcs : MultiVec<u32> = utils::read_obj(_bcs).expect("read bcs");

    let ovl = build_bc_scaffold_graph(graph, bcs, 5, min_ctg, min_bcs, max_bcs);

    let mut _wtr = File::create(mm_file).unwrap();
    let mut wtr = BufWriter::new(_wtr);

    for (i,j,v) in ovl {
        writeln!(wtr, "{}, {}, {}", i, j, v).unwrap();
    }
}
    

    /*
    let mut scaf_graph = Graph::<u32, f32>::new();

    let mut nodeIdxs = HashMap::new();
    for e in to_consider {
        let node = scaf_graph.add_node(e as u32);
        nodeIdxs.insert(e, node);
    }

    let mut sink_buf = Vec::new();
    for (src, sinks) in &overlaps.into_iter().group_by(|x| x.0) {
        sink_buf.clear();
        sink_buf.extend(sinks);
        sink_buf.sort_by_key(|x| x.2 as i32);

        for i in 0 .. min(sink_buf.len(), max_links) {
            let (src, sink, pval) = sink_buf[i];
            scaf_graph.add_edge(nodeIdxs[&src], nodeIdxs[&sink], pval);
        }
    }

    scaf_graph
    }
    */
