//
// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
//

//! Final assembly command
//!
//! Assemble short contigs into larger sequences, terminated by branching extensions or hard stops
//! Compute the set of barcodes that touched each sequence in the graph.

use std::fs::File;
use std::io::BufWriter;
use std::path::{Path, PathBuf};
use utils;
use debruijn;
use std::iter;
use utils::MultiVec;
use kmer::{Lmer, Exts, K};
use martian::{JsonDict, MartianStage};
use rustc_serialize::json::Json;
use std::collections::BTreeMap;
use std::fs;
use std::string;


pub fn main_main_asm(shard_asms: Vec<PathBuf>, shard_bcs: Vec<PathBuf>, out_graph: &Path, out_edge_bcs: Option<&Path>, out_hair: Option<&Path>, censor_sedges: Option<&Path>, write_bv: bool)
{
    // Generate a graph from the shard assemblies
    let (mut graph, n_sedges) =
    {
        
        let censored_sedges: Option<Vec<u32>> = censor_sedges.map(|f| {
            info!("Loading censored edges");
            utils::read_obj(f).expect("loading censor_sedges")
        });


        let mut sedges : Vec<(Lmer, Exts)> = Vec::new();
        for p in shard_asms {
            let shard_sedges: Vec<(Lmer, Exts)> = utils::read_obj(&p).expect("loading shard");
            sedges.extend(shard_sedges);
        }

        let n_sedges = sedges.len();
        info!("Loaded {} sedges", n_sedges);

        info!("Building full graph");
        let (graph, _) = debruijn::build_edges(sedges, censored_sedges);
        info!("Built graph. Got {} edges", graph.len());

        (graph, n_sedges)
    };

    // Generate a list of sedges involved in hairy edges.  We will drop these
    out_hair.map(|f| {
        info!("Finding hairy edges");
        let mut hair_sedges = Vec::new();
        let hair = graph.find_hairy_edges(2*K);

        for edge_id in hair {
            let e = graph.get(edge_id);
            for sedge_idx in e.sedges.unwrap() {
                hair_sedges.push(*sedge_idx);
            }
        }

        info!("Writing hair");
        utils::write_obj(&hair_sedges, f).expect("write failed");
    });

    let mut sedge_to_edge : Vec<isize> = iter::repeat(-1).take(n_sedges).collect();
    for (edge_idx, edge) in Iterator::enumerate(graph.iter()) {
        for sedge_idx in edge.sedges.unwrap()
        {
            sedge_to_edge[*sedge_idx as usize] = edge_idx as isize;
        }
    }

    // Dont save the sedge list of the graph
    graph.edge_sedges = None;
    info!("Writing graph");

    // Write to SN 'basevector' format if requested
    if write_bv {
        let mut _wtr = File::create(out_graph).unwrap();
        let mut wtr = BufWriter::new(_wtr);
        graph.write_to_sn_format(&mut wtr);
    } else {
        utils::write_obj(&graph, out_graph).expect("write failed");
    }

    if out_edge_bcs.is_some() {

        info!("Collating BCs for edges");
        // Keep the BC list for each edge in an array of HashSets
        let mut bcs_sets : Vec<Vec<u32>> = Vec::with_capacity(graph.len());
        for _ in 0..graph.len() {
            bcs_sets.push(Vec::new());
        }

        // Iterate over shards, add the BCS for each sedge to the set for each edge
        let mut shard_idx_base = 0;
        for (i, shard_file) in shard_bcs.iter().enumerate() {
            if i % 100 == 0 { info!("Coalating BC for shard: {}", i); }

            let shard_bcs : MultiVec<u32> = utils::read_obj(&shard_file).expect("read");

            for sedge_idx in 0..shard_bcs.len() {
                let global_sedge_idx = sedge_idx + shard_idx_base;
                let edge_idx = sedge_to_edge[global_sedge_idx];

                if edge_idx >= 0 {
                    let mut bcs_set = bcs_sets.get_mut(edge_idx as usize).unwrap();
                    

                    if bcs_set.len() < 20000 {
                        bcs_set.extend(shard_bcs.get_slice(sedge_idx));
                    }
                }
            }

            // Occasionally compact-ify the bcs_sets
            if i % 64 == 0
            {
                info!("Compacting BC sets");
                // Occasionally compact-ify the bcs_sets
                for bcs_set in bcs_sets.iter_mut()
                {
                    bcs_set.sort();
                    bcs_set.dedup();
                }
            }

            shard_idx_base += shard_bcs.len();
        }

        for bcs_set in bcs_sets.iter_mut()
        {
            bcs_set.sort();
            bcs_set.dedup();
        }

        // Pull out per-edge BC lists into a MultiVec
        let mut edge_bcs = MultiVec::new();

        for bc_set in bcs_sets {
            edge_bcs.add(bc_set);
        }

        info!("writing edge bcs");
        utils::write_obj(&edge_bcs, out_edge_bcs.unwrap()).expect("write failed");
    }
}


pub struct MainAsmMartian;

impl MartianStage for MainAsmMartian {
    fn split(&self, _: JsonDict) -> JsonDict {
        panic!("non-splitting stage");
    }

    fn main(&self, args: JsonDict, outs: JsonDict) -> JsonDict {
        let sedge_asm: Vec<PathBuf> = args["sedge_asm"].as_array().unwrap().iter().map(|v| PathBuf::from(v.as_string().unwrap())).collect();
        let sedge_bcs: Vec<PathBuf> = args["sedge_bcs"].as_array().unwrap().iter().map(|v| PathBuf::from(v.as_string().unwrap())).collect();

        let out_asm = Path::new({outs.get("asm_graph").unwrap().as_string().unwrap().clone()});
        let out_bcs = Path::new({outs.get("asm_bcs").unwrap().as_string().unwrap().clone()});
        let censor_file = PathBuf::from({outs.get("asm_graph").unwrap().as_string().unwrap().clone().to_string() + ".hair"});

        main_main_asm(sedge_asm.clone(), sedge_bcs.clone(), &out_asm, None, Some(&censor_file), None, false);
        main_main_asm(sedge_asm, sedge_bcs, &out_asm, Some(out_bcs), None,  Some(&censor_file), false);
        outs.clone()
    }

    fn join(&self, _: JsonDict, _: JsonDict, _: Vec<JsonDict>, _: Vec<JsonDict>) -> JsonDict {
        panic!("Non-splitting stage");
    }
}


pub struct MainAsmMartianSupernova;

impl MartianStage for MainAsmMartianSupernova {
    fn split(&self, _: JsonDict) -> JsonDict {
        let mut bt = BTreeMap::new();
        bt.insert("__mem_gb".to_string(), Json::F64(64.0));
        let obj = Json::Object(bt);
        let mut chunks = Vec::new();
        chunks.push(obj);
        let chunk_array = Json::Array(chunks);
        let mut cc =  BTreeMap::new();
        cc.insert("chunks".to_string(), chunk_array);
        cc
    }

    fn main(&self, args: JsonDict, outs: JsonDict) -> JsonDict {
        let sedge_asm: Vec<PathBuf> = args["sedge_asm"].as_array().unwrap().iter().map(|v| PathBuf::from(v.as_string().unwrap())).collect();
        let sedge_bcs: Vec<PathBuf> = args["sedge_bcs"].as_array().unwrap().iter().map(|v| PathBuf::from(v.as_string().unwrap())).collect();

        let out_asm = Path::new({outs.get("asm_graph").unwrap().as_string().unwrap().clone()});
        main_main_asm(sedge_asm, sedge_bcs, &out_asm, None, None, None, true);
        outs.clone()
    }

    fn join(&self, _: JsonDict, outs: JsonDict, _: Vec<JsonDict>, chunk_outs: Vec<JsonDict>) -> JsonDict {
        let chunk_name = chunk_outs[0].get("asm_graph").unwrap().as_string().unwrap();
        let outs_name = outs.get("asm_graph").unwrap().as_string().unwrap();
        fs::rename( chunk_name, outs_name ).unwrap();
        let mut obj = BTreeMap::new();
        obj.insert( "asm_graph".to_string(), Json::String( outs_name.to_string() ) );
        obj
    }
}
