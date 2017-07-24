//
// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
//

//! Command for assembling MSP shard Kmers into short DeBruijn graph snippets
//!
//! Enumerates all Kmers observed in the MSP shard and the barcodes observed on each.  Filters
//! the kmers by an arbitrary rule.  Then labels each Kmer with all the observed extensions from it.
//! Note, the extension in this phase will be a super-set of to Kmers that pass the filtering.
//! This will be resolved later.

use std::path::{Path, PathBuf};
use rustc_serialize::json::Json;
use std::collections::{HashMap, BTreeMap};
use std::fs;
use std::hash::{Hash, BuildHasher};
use std::env;
use std::mem;
use rayon::prelude::*;
use rayon;

use martian::{JsonDict, MartianStage};
use utils;
use debruijn;
use shardio::shard::ShardReaderSet;
use kmer;

pub fn est_gb_hashmap<K: Eq + Hash, V, S: BuildHasher>(hash_map: HashMap<K,V,S>) -> f32 {
    11.0/10.0 * ((hash_map.capacity() * (mem::size_of::<K>() + mem::size_of::<V>() + mem::size_of::<u64>())) as f32) / 1e9
}

pub fn est_gb_vec<V>(vec: Vec<V>) -> f32 {
    (vec.capacity() * ( mem::size_of::<V>())) as f32 / 1e9
}


pub fn main_shard_asm(min_kmer_obs: u32, chunk_id: usize, total_chunks: usize, shard_chunks: Vec<PathBuf>, sedge_asm_out: &Path, sedge_bcs_out: &Path, num_threads: usize) {

    let reader = ShardReaderSet::open(&shard_chunks, kmer::BspSer{});
    let my_shards: Vec<(usize, usize)> = (0..reader.num_shards()).filter(|x| x % total_chunks == chunk_id).enumerate().collect();
    let my_total_shards = my_shards.len();

    if num_threads > 0 {
        println!("setup rayon with {} threads", num_threads);
        let cfg = rayon::Configuration::new().set_num_threads(num_threads);
        rayon::initialize(cfg);
    }

    // process each shard, and collect results
    let mut results = Vec::new();
    my_shards.into_par_iter().map(|(shard_count, shard_id)| {
        info!("Processing shard {} of {}", shard_count, my_total_shards);
        let mut bsps = Vec::new();
        reader.read_shard(shard_id, &mut bsps);
        let (valid_kmers, all_kmers) = utils::process_kmer_shard_opt(&bsps, min_kmer_obs, false);
        let kmer_extensions = debruijn::kmer_extensions(&valid_kmers, &all_kmers, &bsps);
        let sedges = debruijn::build_sedges(&kmer_extensions);

        info!("Shard:{}. BPSs:{}, valid_kmers:{}, all_kmers:{}", shard_count, est_gb_vec(bsps), est_gb_hashmap(valid_kmers), est_gb_hashmap(all_kmers));
        //let mut bc_vec = Vec::new();
        let mut sedge_bcs: utils::MultiVec<u32> = utils::MultiVec::new();

        info!("Skipping BC counting");
        /*
        info!("Getting barcodes for sedges");
        for &(sedge, _) in sedges.iter() {
            debruijn::barcodes_for_sedge(&mut bc_vec, sedge, &valid_kmers);
            sedge_bcs.add(bc_vec.drain(..));
        }
        */

        (sedges, sedge_bcs)
    }).collect_into(&mut results);


    let mut all_sedges = Vec::new();
    let mut all_sedge_bcs = utils::MultiVec::new();

    for (mut sedges, sedge_bcs) in results {
        all_sedges.extend(sedges.drain(..));

        for i in 0..sedge_bcs.start_pos.len() {
            all_sedge_bcs.add_slice(sedge_bcs.get_slice(i));
        }
    }

    info!("Got sedges: {}", all_sedges.len());

    info!("Writing sedges: {:?}", sedge_asm_out);
    utils::write_obj(&all_sedges, sedge_asm_out).expect("write failed");

    info!("Writing sedge bcs: {:?}", sedge_bcs_out);
    utils::write_obj(&all_sedge_bcs, sedge_bcs_out).expect("write failed");
}


pub struct ShardAsmMartian;

impl MartianStage for ShardAsmMartian {
    fn split(&self, _: JsonDict) -> JsonDict {

        let total_chunks = 256;

        let mut chunks = Vec::new();
        for chunk_id in 0..total_chunks
        {
            let mut bt = BTreeMap::new();
            bt.insert("chunk_id".to_string(), Json::I64(chunk_id));
            bt.insert("total_chunks".to_string(), Json::I64(total_chunks));
            bt.insert("__threads".to_string(), Json::I64(4));
            bt.insert("__mem_gb".to_string(), Json::F64(36.0));
            chunks.push(Json::Object(bt));
        }

        let chunk_array = Json::Array(chunks);
        let mut cc =  BTreeMap::new();
        cc.insert("chunks".to_string(), chunk_array);
        cc
    }

    fn main(&self, args: JsonDict, outs: JsonDict) -> JsonDict {

        let min_kmer_obs = args["min_kmer_obs"].as_u64().unwrap() as u32;
        let shard_chunks = args["chunks"].as_array().unwrap().iter().map(|v| PathBuf::from(v.as_string().unwrap())).collect();

        let chunk_id = args["chunk_id"].as_u64().unwrap() as usize;
        let total_chunks = args["total_chunks"].as_u64().unwrap() as usize;
        let threads = args["__threads"].as_u64().unwrap() as usize;

        let sedge_asm = Path::new({outs.get("sedge_asm").unwrap().as_string().unwrap().clone()});
        let sedge_bcs = Path::new({outs.get("sedge_bcs").unwrap().as_string().unwrap().clone()});


        main_shard_asm(min_kmer_obs, chunk_id, total_chunks, shard_chunks, sedge_asm, sedge_bcs, threads);
        outs.clone()
    }

    fn join(&self, _: JsonDict, _: JsonDict, _: Vec<JsonDict>, chunk_outs: Vec<JsonDict>) -> JsonDict {

            let mut sedge_asms: Vec<Json> = Vec::new();
            let mut sedge_bcss : Vec<Json> = Vec::new();

            let pwd = env::current_dir().unwrap();

            for (idx, c) in chunk_outs.into_iter().enumerate() {

                //let _old_name =
                let old_name = c.get("sedge_asm").unwrap().as_string().unwrap();
                let mut new_name = pwd.clone();
                new_name.push(format!("chunk{}.sedge_asm", idx));
                fs::rename(old_name, new_name.clone()).unwrap();
                sedge_asms.push(Json::String(new_name.to_str().unwrap().to_string()));

                let old_name = c.get("sedge_bcs").unwrap().as_string().unwrap();
                let mut new_name = pwd.clone();
                new_name.push(format!("chunk{}.sedge_bcs", idx));
                fs::rename(old_name, new_name.clone()).unwrap();
                sedge_bcss.push(Json::String(new_name.to_str().unwrap().to_string()));
            }

            let mut obj = BTreeMap::new();
            obj.insert("sedge_asm".to_string(), Json::Array(sedge_asms));
            obj.insert("sedge_bcs".to_string(), Json::Array(sedge_bcss));
            obj
    }
}
