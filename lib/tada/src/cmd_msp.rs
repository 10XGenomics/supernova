//
// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
//

//! Command to generate barcoded MSP substrings from FASTQ data
//!
//! For each input read, compute the MSP substrings. Each substring is at most 2K-P long, and
//! is stored in a fixed-size struct Bsp, along with the barcode id, read id, and read position,
//! and an Exts bitfield indicating the neighboring bases (if any).

#![allow(dead_code)]

use std::path::{Path, PathBuf};
use std::collections::{BTreeMap, HashMap};
use std::env;

use utils::BcIndexer;
use shardio::shard::{ShardWriteManager, ShardSender};
use multifastq::{MultiFastqIter, InputRead, Fastq};
use msp::simple_scan;
use martian::{JsonDict, MartianStage};
use rustc_serialize::json::{Json};
use crossbeam;
use std::mem;

use kmer::{self, Pmer, Bsp, BspSer, base_to_bits};
use utils;

const K: usize = kmer::K;
const P: usize = 8;
const MEM_ALLOC_GB: usize = 8;

/// Entry point for an MSP run. bc_wl is the path to the BC whitelist file, which will be used to
/// convert BC sequences to ints. fastq is set of pseudo-FASTQ files created by the _SORT_FASTQ_BY_BARCODE
/// pipeline. Permutation is the permutation of p-mers to use, serialized into the give path.
/// The best permutation to use is determined by msp_get_pmer_permutation below.
#[inline(never)]
pub fn main_msp(trim_min_qual: u8, bc_wl: &Path, fastq: Vec<PathBuf>, permutation: &Path, out_path: &Path) {

    let bc_idx = &BcIndexer::new(bc_wl);
    let permutation : &Vec<usize> = &utils::read_obj(permutation).expect("couldn't read permutation");

    let mem = (MEM_ALLOC_GB - 1) * (10 as usize).pow(9);
    let nshards = 8192;
    let bsp_size = mem::size_of::<Bsp>();
    let buf_per_shard = mem / nshards / bsp_size * 2;       // the *2 here is an awful fudge
    info!("MSP mem config: mem: {}, bsp_size: {}, Shards: {}, Buffered items per shard: {}", mem,
            bsp_size,  nshards, buf_per_shard);

    let shard = ShardWriteManager::new(out_path, buf_per_shard, nshards, 2, BspSer{});

    let mut nreads : u64 = 0;


    crossbeam::scope(|scope| {

        let mut handles = Vec::new();

        for f in fastq.chunks(2) {
            let chnk = Vec::from(f);
            let sender = shard.get_sender();

            let mut handle = scope.spawn(move ||
            {
                let fq_src = chnk.iter().flat_map(|p| MultiFastqIter::new(&p, &bc_idx));
                let nr = process_fastq_set(trim_min_qual, K, P, &permutation, fq_src, sender);
                return nr;
            } );
            handles.push(handle);
        }

        while let Some(handle) = handles.pop() {
            nreads += handle.join();
        }

//        handles.iter().flat_map( |h| nreads += h.join() );

    });
    info!( "MSP read count: {}", nreads);
}


/// Load a sampling of reads and determine the distribution of pmers.  We will redefine the order of
/// pmers to make the most common by the lexicographically highest, which will gives more uniform
/// chunk sizes.
pub fn msp_get_pmer_permutation(fastq: Vec<PathBuf>, bc_wl: &Path, permutation_path: &Path) {
    let bc_idx = &BcIndexer::new(bc_wl);
    let fq_src_small = fastq.iter().flat_map(|p| MultiFastqIter::new(&p, &bc_idx).take(25000));
    let pmer_counts = count_pmers_fastq(fq_src_small, P);
    let permutation = pmer_inverse_freq_permutation(P, &pmer_counts);

    info!("Writing pmer permutation: {:?}", permutation_path);
    utils::write_obj(&permutation, permutation_path).expect("write failed");
}


/// Do the MSP splitting of sequence seq, and pack the MSP partitions into Bsp structs,
/// and send these struct out to disk.
#[inline(never)]
fn msp_read(k: usize,
            p: usize,
            partition: u32,
            permutation: &Vec<usize>,
            read: u16,
            seq: &[u8],
            shard_sender: &mut ShardSender<Bsp>) {

    // Can't do anything with strings shorter than K
    if seq.len() < k { return; }

    let msp_parts = simple_scan(k, p, seq, permutation);

    for (msp_key, _, start_pos, slice_length) in msp_parts {
        assert!(slice_length >= k);
        let b = Bsp::new(seq, start_pos, slice_length, partition, read, msp_key as u16);

        if (b.len() as usize) < K || (b.len() as usize) != b.sequence.len()
        {
            println!("bad bsp: {:?}", b);
            panic!("bad bsp: {:?}", b);
        }

        shard_sender.send(b);
    }
}

/// Find the longest prefix of read such that all bases in the final kmer have qv >= min_qual
/// Should implement https://github.com/10XDev/supernova/blob/master/lib/assembly/src/paths/long/BuildReadQGraph48.cc#L70
fn find_trim_len(fq: &Fastq, min_qual: u8) -> usize {

    let qvs = fq.qual.clone().into_bytes();
    let mut good = 0;

    for i in (0..qvs.len()).rev() {
        if qvs[i]-33 < min_qual {
            good = 0;
        } else {
            good += 1;
            if good == K {
                return i + K;
            }
        }
    }

    return 0;
}

fn process_fastq_set<T>(min_qual: u8,
                        k: usize,
                        p: usize,
                        permutation: &Vec<usize>,
                        input: T,
                        mut shard_sender: ShardSender<Bsp>) -> u64
    where T: IntoIterator<Item = InputRead>
{
    let mut iter = input.into_iter();
    let mut nreads : u64 = 0;
    loop {
        match iter.next() {
            Some(inp) => {
                // R1 is numbered with an even number, R2 is the next value
                let trim_len1 = find_trim_len(&inp.r1, min_qual);
                let r1b: Vec<u8> = inp.r1.read.into_bytes().into_iter().map(base_to_bits).collect();
                let r1_slice = &r1b[..trim_len1];
                msp_read(k,
                         p,
                         inp.bc_id,
                         permutation,
                         inp.read_id,
                         r1_slice,
                         &mut shard_sender);

                let trim_len2 = find_trim_len(&inp.r2, min_qual);
                let r2b: Vec<u8> = inp.r2.read.into_bytes().into_iter().map(base_to_bits).collect();
                let r2_slice = &r2b[..trim_len2];
                msp_read(k,
                         p,
                         inp.bc_id,
                         permutation,
                         inp.read_id + 1,
                         r2_slice,
                         &mut shard_sender);
                nreads += 2;
            }
            None => break,
        }
    }

    return nreads;
}

fn pmer_inverse_freq_permutation(p: usize, counts: &HashMap<Pmer, usize>) -> Vec<usize> {
    let mut freq_vec: Vec<(usize, Pmer)> = Vec::new();

    for pmer in Pmer::enumerate(p) {
        match counts.get(&pmer) {
            Some(count) => freq_vec.push((*count, pmer)),
            None => freq_vec.push((0, pmer)),
        }
    }

    freq_vec.sort();

    let mut perm: Vec<(usize, usize)> = freq_vec.into_iter().enumerate().map(|(idx, (_, pmer))| (pmer.value(), idx)).collect();
    perm.sort();
    perm.into_iter().map(|(_, idx)| idx).collect()
}


fn count_pmers_fastq<T: Iterator<Item = InputRead>>(mut input: T,
                                                    p: usize)
                                                    -> HashMap<Pmer, usize> {
    let mut counter = HashMap::new();
    loop {
        match input.next() {
            Some(inp) => {

                // R1 is numbered with an even number, R2 is the next value
                let r1b: Vec<u8> = inp.r1.read.into_bytes().into_iter().map(base_to_bits).collect();
                let r1_slice = &r1b[..];
                let pmers = Pmer::pmers_from_string(r1_slice, p);

                for pmer in pmers {
                    let (min_pmer, _) = pmer.min_rc();
                    let e = counter.entry(min_pmer).or_insert_with(|| 0);
                    *e = *e + 1;
                }

                let r2b: Vec<u8> = inp.r2.read.into_bytes().into_iter().map(base_to_bits).collect();
                let r2_slice = &r2b[..];
                let pmers = Pmer::pmers_from_string(r2_slice, p);

                for pmer in pmers {
                    let (min_pmer, _) = pmer.min_rc();
                    let e = counter.entry(min_pmer).or_insert_with(|| 0);
                    *e = *e + 1;
                }

            }
            None => break,
        }
    }

    counter
}

pub struct MspMartian;


impl MartianStage for MspMartian {
    fn split(&self, args: JsonDict) -> JsonDict {

        let fastqs : Vec<String> = args["fastqs"].as_array().unwrap().into_iter().map(|x| x.as_string().unwrap().to_string()).collect();
        //let fastqs_fofn = args["fastqs"].as_string().unwrap();
        //let fastqs = utils::read_fofn(Path::new(fastqs_fofn));

        let _whitelist = args["barcode_whitelist"].as_string().unwrap();
        let whitelist  = Path::new(_whitelist);

        let mut permutation_file = env::current_dir().unwrap();
        permutation_file.push("permutation.perm");

        // Compute the permutation
        let paths = fastqs.clone().iter().map(PathBuf::from).collect();
        msp_get_pmer_permutation(paths, &whitelist, &permutation_file);

        let mut chunks : Vec<Json> = Vec::new();
        for in_files in fastqs.chunks(8)
        {
            let fastqs = in_files.iter().map(|c| Json::String(c.clone())).collect();

            let mut bt = BTreeMap::new();
            bt.insert("chunk".to_string(), Json::Array(fastqs));
            let p2s = permutation_file.clone().into_os_string().into_string().unwrap();
            bt.insert("permutation".to_string(), Json::String(p2s));
            bt.insert("__mem_gb".to_string(), Json::F64(MEM_ALLOC_GB as f64));
            bt.insert("__threads".to_string(), Json::I64(4));
            chunks.push(Json::Object(bt));
        }


        let chunk_array = Json::Array(chunks);
        let mut cc =  BTreeMap::new();
        cc.insert("chunks".to_string(), chunk_array);
        cc
    }

    fn main(&self, args: JsonDict, outs: JsonDict) -> JsonDict {
        let chunk = args["chunk"].as_array().unwrap().iter().map(|v| PathBuf::from(v.as_string().unwrap())).collect();
        //let chunk  = vec![PathBuf::from(_chunk)];

        let _whitelist = args["barcode_whitelist"].as_string().unwrap();
        let whitelist  = Path::new(_whitelist);

        let _permutation = args["permutation"].as_string().unwrap();
        let permutation = Path::new(_permutation);

        let trim_min_qual = args["trim_min_qual"].as_u64().unwrap() as u8;

        let _out_fn = {outs.get("chunks").unwrap().as_string().unwrap().clone()};
        let out_fn = Path::new(_out_fn);

        let files = main_msp(trim_min_qual, &whitelist, chunk, &permutation, &out_fn);
        let mut final_outs = outs.clone();
        final_outs
    }

    fn join(&self, _: JsonDict, _: JsonDict, _: Vec<JsonDict>, chunk_outs: Vec<JsonDict>) -> JsonDict {

            let mut chunk_vec: Vec<Json> = Vec::new();

            for c in chunk_outs {
                let chunk = c.get("chunks").unwrap();
                chunk_vec.push(chunk.clone());
            }

            let mut obj = BTreeMap::new();
            obj.insert("chunks".to_string(), Json::Array(chunk_vec));
            obj
    }
}


#[cfg(test)]
mod tests {
    use multifastq::Fastq;
    use kmer::K;

    #[test]
    fn test_qv_trim_read() {
        let seq = "TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAC".to_string().into_bytes();
        let quals = "FFFFFFFFIFFFFFFFFFFIIIFFBFIFFFFIFBFFIBFIFBBFFIFFIFFFFFFFFFFFBBBBBBBBBB07BB7BB<BBBBBBBBBB".to_string().into_bytes();
        

        for i in 0..quals.len() {
            // put in a bad qv and make sure we get the right length
            let mut myquals = quals.clone();
            myquals[i] = 34 as u8;

            let fq = Fastq {
                read: String::from_utf8(seq.clone()).unwrap(),
                qual: String::from_utf8(myquals).unwrap(),
            };

            let trim_length = super::find_trim_len(&fq, 10);

            let expected_length = if i < quals.len() - K { quals.len() } else if i >= K { i } else { 0 };
            assert_eq!(expected_length, trim_length);
        }
    }
}
