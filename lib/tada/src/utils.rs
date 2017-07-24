//
// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
//

//! Various useful methods (could be cleaned up)

use std::fs::File;
use std::path::{Path};
use std::collections::HashMap;
use std::io::{BufReader, BufWriter};
use std::io::BufRead;
use kmer::{Kmer, Bsp};
use smallvec::SmallVec;
use std::mem;
use std::fmt::Debug;
use std::str::FromStr;
use std::io;

use bincode;
use rustc_serialize::{Encodable, Decodable};
use bincode::rustc_serialize::{EncodingResult, DecodingResult, encode_into, decode_from};

use itertools::{Itertools};
use fx::FxHashMap;


pub fn write_obj<T: Encodable, P: AsRef<Path> + Debug>(g: &T, filename: P) -> EncodingResult<()> {
    let f = match File::create(&filename) {
        Err(err) => panic!("couldn't create file {:?}: {}", filename, err),
        Ok(f) => f,
    };
    let mut writer = BufWriter::new(f);
    encode_into(&g, &mut writer, bincode::SizeLimit::Infinite)
}

pub fn read_obj<T: Decodable, P: AsRef<Path> + Debug>(filename: P) -> DecodingResult<T> {
    let f = match File::open(&filename) {
        Err(err) => panic!("couldn't open file {:?}: {}", filename, err),
        Ok(f) => f,
    };
    let mut reader = BufReader::new(f);
    decode_from(&mut reader, bincode::SizeLimit::Infinite)
}




#[derive(Clone, RustcEncodable, RustcDecodable)]
pub struct  MultiVec<T> {
    pub items: Vec<T>,
    pub start_pos: Vec<usize>,
    pub vec_len: Vec<u32>,
}

impl<T : Clone> MultiVec<T> {
    pub fn new() -> MultiVec<T> {
        MultiVec {
            items: Vec::new(),
            start_pos: Vec::new(),
            vec_len: Vec::new(),
        }
    }

    pub fn add<S: IntoIterator<Item = T>>(&mut self, items: S) {
        let start_pos = self.items.len();
        self.start_pos.push(start_pos);

        let mut n = 0;
        for i in items {
            self.items.push(i);
            n += 1;
        }

        self.vec_len.push(n);
    }

    pub fn add_slice(&mut self, items: &[T]) {
        let start_pos = self.items.len();
        self.start_pos.push(start_pos);

        let mut n = 0;
        for i in items {
            self.items.push(i.clone());
            n += 1;
        }

        self.vec_len.push(n);
    }

    pub fn len(&self) -> usize {
        return self.start_pos.len();
    }

    pub fn get_slice(&self, sub_vec: usize) -> &[T]
    {
        &self.items[(self.start_pos[sub_vec])..(self.start_pos[sub_vec] + self.vec_len[sub_vec] as usize)]
    }
}


pub struct BcIndexer {
    bc_map: HashMap<String, u32>,
    num_bcs: u32,
}

impl BcIndexer {
    pub fn new<P:AsRef<Path>>(filename: P) -> BcIndexer {
        let f = File::open(filename.as_ref()).unwrap();
        Self::from_reader(f)
    }

    pub fn from_reader<R: io::Read>(reader: R) -> BcIndexer {
        let br = BufReader::new(reader);

        let mut bc_map = HashMap::new();

        let mut i = 0u32;
        for l in br.lines() {
            bc_map.insert(l.unwrap(), i);
            i += 1;
        }

        BcIndexer {
            bc_map: bc_map,
            num_bcs: i 
        }
    }

    /// Parse a BC string of the form 'ACGATA-1', or 'ACGAT' into a 
    /// numerical bc id (if the whole sequence or the subsequence before the first dash)
    /// and a gem group ID.
    pub fn get_bc_parts(&self, bc: &str) -> (Option<u32>, Option<u8>) {
        if bc.contains('-') {
            let mut parts = bc.split('-');
            let bc_seq = parts.next().unwrap();
            let bc_id = self.bc_map.get(bc_seq).cloned();
            let gem_group_str = parts.next().unwrap();
            let gem_group = u8::from_str(gem_group_str).expect("invalid gem group string");

            (bc_id, Some(gem_group))
        } else {
            (self.bc_map.get(bc).cloned(), None)
        }
    }

    /// Return a BC id for a BC string.  
    /// bc_id = bc_number + (gem_group - 1) * num_bcs + 1
    /// BC ID 0 is reserved, so the bc_is offset by 1.
    /// The GEM group increments the id by num_bcs.
    /// BC sequences lacking a gem group suffix are assumed to have gem group 1.
    pub fn get_bc_id(&self, bc: &str) -> Option<u32> {
        let (bc, gg) = self.get_bc_parts(bc);

        match (bc, gg) {
            (None, _) => None,
            (Some(b), gg) => {
                let g = gg.unwrap_or(1);
                let base = ((g as u32) - 1).checked_mul(self.num_bcs).expect("too many gem groups - BC id overflowed");
                let r = base + b + 1;
                Some(r)
            },
        }
    }
}


pub fn load_barcode_whitelist2<P: AsRef<Path>>(filename: P) -> HashMap<Vec<u8>, u32> {
    let f = File::open(filename).expect("couldn't load barcode whitelist");
    let br = BufReader::new(f);

    let mut bc_map = HashMap::new();

    let mut i = 1u32;
    for l in br.lines() {
        bc_map.insert(l.unwrap().into_bytes(), i);
        i += 1;
    }

    bc_map
}

pub fn load_barcode_whitelist_inverse<P: AsRef<Path>>(filename: P) -> HashMap<u32, Vec<u8>> {
    let f = File::open(filename).unwrap();
    let br = BufReader::new(f);

    let mut bc_map = HashMap::new();

    let mut i = 1u32;
    for l in br.lines() {
        bc_map.insert(i, l.unwrap().into_bytes());
        i += 1;
    }

    bc_map
}

/// Read a shard and determine the valid kmers
#[inline(never)]
pub fn process_kmer_shard2(bsps: &Vec<Bsp>) -> FxHashMap<Kmer, Vec<(u32, u16)>> {
    let mut kmer_dict = FxHashMap();

    for b in bsps {

        //if (b.length as usize) < K || (b.length as usize) != b.sequence.length() || b.pos + b.length > 114
        //{
        //    println!("bad bsp: {:?}", b);
        //    panic!("bad bsp: {:?}", b);
        //}

        let pri = (b.partition, b.read);
        for k in b.kmers() {
            let e = kmer_dict.entry(k.min_rc()).or_insert_with(|| Vec::new());
            (*e).push(pri)
        }
    }

    // filter out kmers we don't like
    let mut final_kmers = FxHashMap();
    let mut observed_kmers = 0;

    for (kmer, pris) in kmer_dict {
        observed_kmers += 1;
        let mut multiple_bcs = false;
        let first_bc = pris[0].0;

        for &(bc_id, _) in pris.iter() {
            if bc_id != first_bc {
                multiple_bcs = true;
                break;
            }
        }

        if multiple_bcs {
            final_kmers.insert(kmer, pris);
        }
    }

    println!("Kmers observed: {}. Kmers accepted: {}", observed_kmers, final_kmers.len());
    final_kmers
}

/// Read a shard and determine the valid kmers
#[inline(never)]
pub fn process_kmer_shard_old(bsps: &Vec<Bsp>, min_kmer_obs: u32) -> (FxHashMap<Kmer, Vec<u32>>, FxHashMap<Kmer, u32>) {

    let mut kmer_buckets : Vec<Vec<(Kmer, u32)>> = Vec::new();
    for _ in 0..256 {
        kmer_buckets.push(Vec::new());
    }

    info!("Enumerating kmers...");
    for b in bsps {
        for k in b.iter_kmers() {
            let min_kmer = k.min_rc();
            let bucket = min_kmer.bucket();
            kmer_buckets[bucket as usize].push((min_kmer, b.partition));
        }
    }

    let mut all_kmer_count: FxHashMap<Kmer, u32> = FxHashMap();
    let mut final_kmers: FxHashMap<Kmer, Vec<u32>> = FxHashMap();
    let mut unique_kmers = 0;

    let mut total_kmers = 0;
    for bucket in kmer_buckets.iter() {
        total_kmers += bucket.len();
    }

    info!("Validating kmers...");
    for mut kmer_vec in kmer_buckets {

        kmer_vec.sort();

        for (kmer, group) in &kmer_vec.into_iter().group_by(|elt| elt.0) {
            let mut pris = SmallVec::<[u32; 8]>::new();
            unique_kmers += 1;

            let mut num_bcs = 0;
            let mut last_bc = 0;
            let mut num_obs = 0;
            for (_, bc) in group {
                if bc != last_bc {
                    num_bcs += 1;
                    last_bc = bc;
                    pris.push(bc);
                }
                num_obs += 1
            }

            if false {
                let test = String::from("TACACTGACACCTCACACGGCAGGGTATTCCAACAGACCTGCAGCTGA").into_bytes();
                let test_kmer = Kmer::kmers_from_string(&test[..])[0]; 
                let test_kmer_rc = test_kmer.rc();

                if kmer == test_kmer {
                    info!("NEIL: fw num_obs={} num_bcs={}", num_obs, num_bcs);
                }

                if kmer == test_kmer_rc {
                    info!("NEIL: rc num_obs={} num_bcs={}", num_obs, num_bcs);
                }
            }

            all_kmer_count.insert(kmer, num_obs);
            if num_bcs > 1 && num_obs >= min_kmer_obs {
                let mut vec = Vec::with_capacity(pris.len());
                vec.extend(pris);
                final_kmers.insert(kmer, vec);
            }
        }
    }

    info!("Total BSPs: {}, Total kmers observed: {}, Unique Kmers observed: {}. Kmers accepted: {}", bsps.len(), total_kmers, unique_kmers, final_kmers.len());
    (final_kmers, all_kmer_count)
}


/// Read a shard and determine the valid kmers
/// Low memory implementation that should consume < 4G of temporary memory
/// To reduce memory consumption, set track_bcs to false to forget about BC lists.
#[inline(never)]
pub fn process_kmer_shard_opt(bsps: &Vec<Bsp>, min_kmer_obs: u32, track_bcs: bool) -> (FxHashMap<Kmer, Vec<u32>>, FxHashMap<Kmer, u32>) {

    let mut all_kmer_count: FxHashMap<Kmer, u32> = FxHashMap();
    let mut final_kmers: FxHashMap<Kmer, Vec<u32>> = FxHashMap();
    let mut unique_kmers = 0;
    let mut total_kmers = 0;

    // Estimate memory consumed by Kmer vectors, and set iteration count appropriately
    let expected_kmers = bsps.len() * 30;
    let kmer_mem = expected_kmers * mem::size_of::<(Kmer, u32)>();
    let max_mem = 4 * (10 as usize).pow(9);
    let slices = kmer_mem / max_mem + 1;
    let sz = 256 / slices + 1;
    
    let mut bucket_ranges = Vec::new();
    let mut start = 0;
    while start < 256 {
        bucket_ranges.push(start..start+sz);
        start += sz;
    }

    if bucket_ranges.len() > 1 {
        info!("processing {} BSPs in {} passes. Bucket ranges: {:?}", bsps.len(), bucket_ranges.len(), bucket_ranges);
    }


    for bucket_range in bucket_ranges {

        let mut kmer_buckets : Vec<Vec<(Kmer, u32)>> = Vec::new();
        for _ in 0..256 {
            kmer_buckets.push(Vec::new());
        }

        info!("Enumerating kmers...");
        for b in bsps {
            for k in b.iter_kmers() {
                let min_kmer = k.min_rc();
                let bucket = min_kmer.bucket() as usize;

                if bucket >= bucket_range.start && bucket < bucket_range.end {
                    kmer_buckets[bucket].push((min_kmer, b.partition));
                }
            }
        }
        
        for bucket in kmer_buckets.iter() {
            total_kmers += bucket.len();
        }

        info!("Validating kmers...");
        for mut kmer_vec in kmer_buckets {

            kmer_vec.sort();

            for (kmer, group) in &kmer_vec.into_iter().group_by(|elt| elt.0) {
                let mut pris = SmallVec::<[u32; 8]>::new();
                unique_kmers += 1;

                let mut num_bcs = 0;
                let mut last_bc = 0;
                let mut num_obs = 0;
                for (_, bc) in group {
                    if bc != last_bc {
                        num_bcs += 1;
                        last_bc = bc;
                        pris.push(bc);
                    }
                    num_obs += 1
                }

                all_kmer_count.insert(kmer, num_obs);
                if num_bcs > 1 && num_obs >= min_kmer_obs {
                    let mut vec = Vec::with_capacity(pris.len());
                    vec.extend(pris);
                    if track_bcs {
                        final_kmers.insert(kmer, vec);
                    } else {
                        final_kmers.insert(kmer, vec![]);
                    }
                }
            }
        }
    }

    info!("Total BSPs: {}, Total kmers observed: {}, Unique Kmers observed: {}. Kmers accepted: {}", bsps.len(), total_kmers, unique_kmers, final_kmers.len());
    (final_kmers, all_kmer_count)
}


#[inline(never)]
pub fn process_kmer_shard(bsps: &Vec<Bsp>, min_kmer_obs: u32) -> (FxHashMap<Kmer, Vec<u32>>, FxHashMap<Kmer, u32>) {
    process_kmer_shard_opt(bsps, min_kmer_obs, true)
}




pub fn read_fofn(fofn: &Path) -> Vec<String> {
    let f = File::open(fofn).unwrap();
    BufReader::new(f).lines().map(Result::unwrap).collect()
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::io;

    const WL_FILE: &'static [u8] = b"ACGTA
ACGTC
ACGTG
ACGTT";

    #[test]
    fn test_bc_indexed() {
        let index = BcIndexer::from_reader(WL_FILE);

        assert_eq!(index.get_bc_id("ACGTA-1"), Some(1));
        assert_eq!(index.get_bc_id("ACGTA"), Some(1));
        assert_eq!(index.get_bc_id("ACGTA-2"), Some(5));

        assert_eq!(index.get_bc_id("ACGTT-1"), Some(4));
        assert_eq!(index.get_bc_id("ACGTT-2"), Some(8));

        assert_eq!(index.get_bc_id("AACG-1"), None);
        assert_eq!(index.get_bc_id("AACG"), None);
    }
}
