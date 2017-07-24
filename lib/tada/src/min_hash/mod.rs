//
// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
//

use fx::FxHasher;
use std::hash::Hash;
use std::marker::PhantomData;
use std::hash::Hasher;
use std::cmp::min;
use std::u64::MAX;
use std::collections::{HashMap, HashSet};
use std::cmp::{Ordering, Ord};
use std::iter::FromIterator;

struct MinHasher<T:Hash>
{
    num_hashes: usize,
    hash_seeds: Vec<u64>,
    phantom: PhantomData<T>,
}

impl<T:Hash> MinHasher<T>
{
    fn new(n: usize) -> MinHasher<T>
    {
        MinHasher
        {
            num_hashes: n,
            hash_seeds: (0u64..(n as u64)).collect(),
            phantom: PhantomData,
        }
    }

    fn min_hash_signature(&self, values: &[T], signature: &mut Vec<u64>)
    {
        signature.clear();
        for seed in self.hash_seeds.iter()
        {
            let mut hasher = FxHasher::new(*seed as usize);
            let mut min_hash : u64 = MAX;
            for v in values.iter()
            {
                v.hash(&mut hasher);
                let hash_value = hasher.finish();
                min_hash = min(min_hash, hash_value);
            }

            signature.push(min_hash);
        }
    }
}

struct MinHashIndex<K, T:Hash>
{
    hasher: MinHasher<T>,
    items: Vec<K>,
    hash_dicts: Vec<HashMap<u64, Vec<usize>>>,
    buf: Vec<u64>,
}

impl<K, T:Hash> MinHashIndex<K, T>
{
    pub fn new(num_hashes: usize) -> MinHashIndex<K,T>
    {
        let hash_dicts = (0..num_hashes).map(|_| { HashMap::new() }).collect();

        MinHashIndex
        {
            hasher: MinHasher::new(num_hashes),
            items: Vec::new(),
            hash_dicts: hash_dicts,
            buf: Vec::with_capacity(num_hashes),
        }
    }

    pub fn insert(&mut self, key: K, set: &[T])
    {
        // Insert item
        let id = self.items.len();
        self.items.push(key);

        // Min-hash signature
        self.hasher.min_hash_signature(&set, &mut self.buf);

        for (mut h, sig_item) in self.hash_dicts.iter_mut().zip(self.buf.iter())
        {
            let mut item_vec = h.entry(*sig_item).or_insert_with(||{Vec::new()});
            item_vec.push(id);
        }
    }

    /// Query the min-hash index for items similar to the qset signature
    /// Returns items from the index with expected Jaccard similarity greater
    /// than min_similarity.
    pub fn query(&mut self, min_similarity: f64, query_set: &Vec<T>) -> Vec<(f64, &K)>
    {
        let mut hit_counts = HashMap::new();

        // Compute the signature of the query set
        self.hasher.min_hash_signature(query_set, &mut self.buf);

        for (h, sig_item) in self.hash_dicts.iter().zip(self.buf.iter())
        {
            for i in h.get(&sig_item).unwrap_or(&Vec::new())
            {
                let counts = hit_counts.entry(*i).or_insert(0);
                *counts = *counts + 1;
            }
        }

        let mut res = Vec::new();
        for (item_id,hits) in hit_counts.iter()
        {
            let jacc = (*hits as f64) / (self.hasher.num_hashes as f64);
            if jacc > min_similarity
            {
                res.push((jacc, self.items.get(*item_id).unwrap()));
            }
        }

        res
    }
}

/// Compute the jaccard index, intersection size and union size of two sorted, deduped arrays
pub fn quick_jaccard<T: Eq + Ord>(s1: &[T], s2: &[T]) -> (f32, usize, usize) {

    let mut p1 = 0;
    let mut p2 = 0;
    let mut intersection = 0;
    let mut union = 0;

    while p1 < s1.len() && p2 < s2.len() {

        match s1[p1].cmp(&s2[p2]) {
            Ordering::Greater => {
                p2 += 1;
                union += 1;
            },
            Ordering::Less => {
                p1 += 1;
                union += 1;
            },
            Ordering::Equal =>
            {
                p1 += 1;
                p2 += 1;
                intersection += 1;
                union += 1;
            }
        }
    }

    union += (s1.len() - p1) + (s2.len() - p2);
    let jacc = intersection as f32 / union as f32;
    (jacc, intersection, union)
}


pub fn jaccard<T: Hash + Eq + Clone + Copy, I1: Iterator<Item=T>, I2: Iterator<Item=T>>(s1: I1, s2: I2) -> (f32, usize, usize) {

    let s1: HashSet<T> = HashSet::from_iter(s1);
    let s2: HashSet<T> = HashSet::from_iter(s2);

    let intersection = s1.intersection(&s2).count();
    let union = s1.union(&s2).count();
    let jacc = intersection as f32 / union as f32;

    (jacc, intersection, union)
}







#[cfg(test)]
mod tests {
    use super::*;
    use rand::{self, Rng};
    use csv;

   quickcheck! {
        fn test_jacc(s1: Vec<u32>, s2: Vec<u32>) -> bool {

            let mut s1 = s1;
            let mut s2 = s2;

            let j1 = jaccard(s1.iter(), s2.iter());

            s1.sort();
            s1.dedup();

            s2.sort();
            s2.dedup();

            let j2 = quick_jaccard(&s1, &s2);
            j1 == j2
        }
   }   



    pub fn rand_set(n: usize, m: usize) -> Vec<usize> {
        let mut r = rand::thread_rng();
        (0..n).map(|_| (r.next_u64()as usize) % m).collect()
    }


    pub fn test_hashes(sz1: usize, sz2: usize, max_items: usize, num_hashes: usize) -> (f32, f32) {
        let set1 = rand_set(sz1, max_items);
        let set2 = rand_set(sz2, max_items);

        let hashes = super::MinHasher::new(num_hashes);

        let mut sig1 = Vec::new();
        hashes.min_hash_signature(&set1, &mut sig1);

        let mut sig2 = Vec::new();
        hashes.min_hash_signature(&set2, &mut sig2);

        let n_match = sig1.iter().zip(sig2.iter()).filter(|&(x1, x2)| x1 == x2).count();

        let jacc = jaccard(sig1.iter(), sig2.iter());

        (jacc.0, n_match as f32 / num_hashes as f32)
    }


    #[test]
    pub fn test_min_hash() {

        let mut wtr = csv::Writer::from_file("stats.csv").expect("open csv").delimiter(b'\t');
        let mut rng = rand::thread_rng();

        for _ in 0 .. 1000 {
            let n1 = rng.gen_range(10, 100);
            let n2 = rng.gen_range(10, 100);
            let c = rng.gen_range(10, 300);
            let r = test_hashes(n1, n2, c, 500);
            wtr.encode(r).expect("csv write error");
        }
    }
}
