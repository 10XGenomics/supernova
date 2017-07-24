// Copyright 2014 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! A fixed-width bit encoding implementation. This allows to store a sequence of values over
//! a reduced alphabet by packing them bit-encoded into a sequence of bytes.

use std::fmt;
use kmer;
use std::hash::{Hash, Hasher};
use kmer::{K, Kmer, reverse_by_twos, Dir};

const BLOCK_BITS: usize = 64;
const WIDTH: usize = 2;

fn mask() -> u64 {
    (1 << WIDTH) - 1
}

/// A sequence of bitencoded values. This implementation is specialized to 2-bit alphabets
#[derive(Clone, PartialEq, Eq, RustcEncodable, RustcDecodable)]
pub struct BitEnc {
    storage: Vec<u64>,
    mask: u64,
    len: usize,
}


impl BitEnc {
    /// Create a new instance with a given encoding width (e.g. width=2 for using two bits per value).
    pub fn new() -> BitEnc {
        BitEnc {
            storage: Vec::new(),
            mask: mask(),
            len: 0,
        }
    }

    /// Create a new instance with a given capacity.
    pub fn with_capacity(n: usize) -> Self {
        BitEnc {
            storage: Vec::with_capacity(n * WIDTH / 64),
            mask: mask(),
            len: 0,
        }
    }

    pub fn from_dna_string(dna: &str) -> BitEnc {
        let mut bitenc = BitEnc {
            storage: Vec::new(),
            mask: mask(),
            len: 0,
        };
        for c in dna.chars() {
            bitenc.push(kmer::base_to_bits(c as u8));
        }
        bitenc
    }

    pub fn from_bytes(bytes: &Vec<u8>) -> BitEnc {
        let mut bitenc = BitEnc {
            storage: Vec::new(),
            mask: mask(),
            len: 0,
        };

        for b in bytes.iter() {
            bitenc.push(*b)
        }

        bitenc
    }

    pub fn to_dna_string(&self) -> String {
        let mut dna: String = String::new();
        for v in self.iter() {
            dna.push(kmer::bits_to_base(v));
        }
        dna
    }

    pub fn to_ascii_vec(&self) -> Vec<u8> {
        let mut res = Vec::new();
        for v in self.iter() {
            res.push(kmer::bits_to_ascii(v));
        }
        res
    }

    /// Append a value.
    pub fn push(&mut self, value: u8) {
        let (block, bit) = self.addr(self.len);
        if bit == 0 {
            self.storage.push(0);
        }
        self.set_by_addr(block, bit, value);
        self.len += 1;
    }

    /// Push values read from a byte array.
    ///
    /// # Arguments
    /// `bytes`: byte array to read values from
    /// `seq_length`: how many values to read from the byte array. Note that this
    /// is number of values not number of elements of the byte array.
    pub fn push_bytes(&mut self, bytes: &Vec<u8>, seq_length: usize) {
        assert!(seq_length <= bytes.len() * 8 / WIDTH,
                "Number of elements to push exceeds array length");

        for i in 0..seq_length {
            let byte_index = (i * WIDTH) / 8;
            let byte_slot = (i * WIDTH) % 8;

            let v = bytes[byte_index];
            let bits = (v >> byte_slot) & (self.mask as u8);

            self.push(bits);
        }
    }

    /// Append `n` times the given value.
    pub fn push_values(&mut self, mut n: usize, value: u8) {
        {
            // fill the last block
            let (block, mut bit) = self.addr(self.len);
            if bit > 0 {
                // TODO use step_by once it has been stabilized: for bit in (bit..64).step_by(self.width) {
                while bit <= 64 {
                    self.set_by_addr(block, bit, value);
                    n -= 1;
                    bit = bit + WIDTH
                }
            }
        }

        // pack the value into a block
        let mut value_block = 0;
        {
            let mut v = value as u64;
            for _ in 0..(64/WIDTH) {
                value_block |= v;
                v <<= WIDTH;
            }
        }

        // push as many value blocks as needed
        let i = self.len + n;
        let (block, bit) = self.addr(i);
        for _ in self.storage.len()..block {
            self.storage.push(value_block);
        }

        if bit > 0 {
            // add the remaining values to a final block
            self.storage.push(value_block >> (64 - bit));
        }

        self.len = i;
    }

    /// Set the value as position `i`.
    pub fn set(&mut self, i: usize, value: u8) {
        let (block, bit) = self.addr(i);
        self.set_by_addr(block, bit, value);
    }

    /// Get the value at position `i`.
    pub fn get(&self, i: usize) -> Option<u8> {
        if i >= self.len {
            None
        } else {
            let (block, bit) = self.addr(i);
            Some(self.get_by_addr(block, bit))
        }
    }

    /// Iterate over stored values (values will be unpacked into bytes).
    pub fn iter(&self) -> BitEncIter {
        BitEncIter {
            bitenc: self,
            i: 0,
        }
    }

    /// Clear the sequence.
    pub fn clear(&mut self) {
        self.storage.clear();
        self.len = 0;
    }

    fn get_by_addr(&self, block: usize, bit: usize) -> u8 {
        ((self.storage[block] >> bit) & self.mask) as u8
    }

    fn set_by_addr(&mut self, block: usize, bit: usize, value: u8) {
        let mask = self.mask << bit;
        self.storage[block] |= mask;
        self.storage[block] ^= mask;
        self.storage[block] |= (value as u64 & self.mask) << bit;
    }

    fn addr(&self, i: usize) -> (usize, usize) {
        let k = i * WIDTH;
        (k / BLOCK_BITS, k % BLOCK_BITS)
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn width(&self) -> usize {
        WIDTH
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Get the prefix of a given size.
    pub fn pref(&self, k: usize) -> BitEnc {
        assert!(k <= self.len, "Prefix size exceeds number of elements.");
        if k == 0 {
            return BitEnc::new();
        }
        // Get position of k-th element (element with index k - 1)
        let (block, bitpos) = self.addr(k - 1);
        // We need to keep the first block-1 elements of the storage vector
        // plus the first bitpos positions of the block-th element.
        let mut new_storage: Vec<u64> = self.storage.iter().take(block + 1).cloned().collect();
        // zero-out all bits after bitpos
        if bitpos + WIDTH < 64 {
            new_storage[block] &= (1 << (bitpos + WIDTH)) - 1;
        }
        BitEnc {
            storage: new_storage,
            mask: self.mask,
            len: k,
        }
    }

    pub fn suf(&self, k: usize) -> BitEnc {
        assert!(k <= self.len, "Suffix size exceeds number of elements.");
        // There's a faster way but involves shifting bits down from "higher" blocks...
        // let (block, bitpos) = self.addr(self.len - k);
        let values: Vec<u8> = self.iter().skip(self.len - k).collect();
        let mut bitenc = BitEnc::new();
        for v in values {
            bitenc.push(v);
        }
        bitenc
    }

    pub fn reverse(&self) -> BitEnc {
        let values: Vec<u8> = self.iter().collect();
        let mut bitenc = BitEnc::new();
        for v in values.iter().rev() {
            bitenc.push(*v);
        }
        bitenc
    }

    pub fn rc(&self) -> BitEnc {
        let mut bitenc = BitEnc::new();
        for i in (0..self.len()).rev() {
            let v = 3 - self.get(i).expect("oob");
            bitenc.push(v);
        }
        bitenc
    }

    // pub fn complement(&self) -> BitEnc {
    //    assert!(self.width == 2, "Complement only supported for 2bit encodings.");
    //    let values: Vec<u32> = Vec::with_capacity(self.len());
    //    for i, v in self.storage.iter() {
    //        values[i] = v;
    //    }
    //    values[values.len() - 1] =
    // }

    pub fn get_kmer_slow(&self, pos: usize) -> kmer::Kmer {
        // superseded by get_kmer

        let mut kmer = kmer::Kmer::empty();
        for (j, i) in (pos..(pos + kmer::K)).enumerate() {
            let v = self.get(i).expect("out of bounds");
            kmer = kmer.set(j, v);
        }
        kmer
    }

    pub fn get_kmer(&self, pos: usize) -> kmer::Kmer {
        let (block, bit) = self.addr(pos);

        if bit == 0 {
            let b0 = reverse_by_twos(self.storage[block]);
            let b1 = reverse_by_twos(self.storage[block+1]);
            kmer::Kmer { storage: [b0, b1 & kmer::Kmer::lower_mask()] }
        } else if bit <= (BLOCK_BITS - 2*(K-BLOCK_BITS/2)) {
            let up_shift = bit;
            let down_shift = 64-bit;
            let b0 = reverse_by_twos(self.storage[block]);
            let b1 = reverse_by_twos(self.storage[block+1]);

            let r0 = b0 << up_shift | b1 >> down_shift;
            let r1 = (b1 << up_shift) & kmer::Kmer::lower_mask();
            kmer::Kmer { storage: [r0, r1 & kmer::Kmer::lower_mask()] }
        } else {

            let up_shift = bit;
            let down_shift = 64-bit;
            let b0 = reverse_by_twos(self.storage[block]);
            let b1 = reverse_by_twos(self.storage[block+1]);
            let b2 = reverse_by_twos(self.storage[block+2]);

            let r0 = b0 << up_shift | b1 >> down_shift;
            let r1 = (b1 << up_shift | b2 >> down_shift) & kmer::Kmer::lower_mask();
            kmer::Kmer { storage: [r0, r1 & kmer::Kmer::lower_mask()] }
        }
    }


    pub fn kmers(&self) -> Vec<kmer::Kmer> {

        let mut vec = Vec::with_capacity(self.len() - kmer::K + 1);
        let first_kmer = self.get_kmer(0);
        vec.push(first_kmer);
        for i in 0..(self.len() - kmer::K) {
            let mut prev_kmer: kmer::Kmer = vec[vec.len() - 1].clone();
            prev_kmer = prev_kmer.extend_right(self.get(i + kmer::K).unwrap());
            vec.push(prev_kmer);
        }
        vec
    }

    pub fn first_kmer(&self) -> kmer::Kmer {
        self.get_kmer(0)
    }

    pub fn last_kmer(&self) -> kmer::Kmer {
        self.get_kmer(self.len() - kmer::K)
    }
}

impl Hash for BitEnc {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.storage.hash(state);
    }
}

impl fmt::Debug for BitEnc {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut s = String::new();
        for pos in 0..self.len() {
            s.push(kmer::bits_to_base(self.get(pos).expect("pos")))
        }

        write!(f, "{}", s)
    }
}


/// Iterator over values of a bitencoded sequence (values will be unpacked into bytes).
pub struct BitEncIter<'a> {
    bitenc: &'a BitEnc,
    i: usize,
}

impl<'a> Iterator for BitEncIter<'a> {
    type Item = u8;

    fn next(&mut self) -> Option<u8> {
        let value = self.bitenc.get(self.i);
        self.i += 1;
        value
    }
}



pub struct BitEncSlice<'a> {
    pub bitenc: &'a BitEnc,
    pub start: usize,
    pub length: usize,
}

impl<'a> BitEncSlice<'a> {
    pub fn len(&self) -> usize {
        self.length
    }

    pub fn kmers(&self) -> Vec<Kmer> {
        let mut kmers = Vec::new();

        let mut k0 = self.bitenc.get_kmer(self.start);
        kmers.push(k0);

        for i in 1..(self.length - K + 1) {
            let next_base = self.bitenc.get(self.start + K - 1 + i).expect("overflow");
            k0 = k0.extend_right(next_base);
            kmers.push(k0);
        }

        kmers
    }

    pub fn first_kmer(&self) -> Kmer {
        self.bitenc.get_kmer(self.start)
    }

    pub fn last_kmer(&self) -> Kmer {
        self.bitenc.get_kmer(self.start + self.length - K)
    }

    pub fn term_kmer(&self, dir: Dir) -> Kmer {
        match dir {
            Dir::Left => self.first_kmer(),
            Dir::Right => self.last_kmer(),
        }
    }

    pub fn ascii(&self) -> Vec<u8> {
        let mut v = Vec::new();
        for pos in self.start..(self.start + self.length) {
            v.push(kmer::bits_to_ascii(self.bitenc.get(pos).expect("pos")));
        }
        v
    }

    pub fn to_dna_string(&self) -> String {
        let mut dna: String = String::new();
        for pos in self.start..(self.start + self.length) {
            dna.push(kmer::bits_to_base(self.bitenc.get(pos).expect("pos")));
        }
        dna
    }

    pub fn to_string(&self) -> String {
        String::from_utf8(self.ascii()).unwrap()
    }

    pub fn to_bitenc(&self) -> BitEnc {
        let mut be = BitEnc::with_capacity(self.length);
        for pos in self.start..(self.start + self.length) {
            be.push(self.bitenc.get(pos).expect("pos"));
        }

        be
    }

    pub fn get(&self, pos: usize) -> u8 {
        self.bitenc.get(self.start + pos).unwrap()
    }
}

impl<'a> fmt::Debug for BitEncSlice<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut s = String::new();
        for pos in self.start..(self.start + self.length) {
            s.push(kmer::bits_to_base(self.bitenc.get(pos).expect("pos")))
        }

        write!(f, "{}", s)
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use kmer;

    // use bio::data_structures::bwt::{bwt};
    // use bio::data_structures::suffix_array::{suffix_array};
    // use bio::alphabets::{dna, Alphabet};
    // use bio::data_structures::fmindex::FMIndex;

    #[test]
    fn test_bitenc() {
        let mut bitenc = BitEnc::new();
        bitenc.push(0);
        bitenc.push(2);
        bitenc.push(1);
        let mut values: Vec<u8> = bitenc.iter().collect();
        assert_eq!(values, [0, 2, 1]);
        bitenc.set(1, 3);
        values = bitenc.iter().collect();
        assert_eq!(values, [0, 3, 1]);
    }

    #[test]
    fn test_push_values() {
        let mut bitenc = BitEnc::new();
        bitenc.push_values(32, 0);
        assert_eq!(bitenc.storage, [0]);
    }

    #[test]
    fn test_push_bytes() {
        let in_values: Vec<u8> = vec![2, 20];

        let mut bitenc = BitEnc::new();
        bitenc.push_bytes(&in_values, 8);
        // Contents should be 00010100 00000010
        let values: Vec<u8> = bitenc.iter().collect();
        assert_eq!(values, [2, 0, 0, 0, 0, 1, 1, 0]);
        assert_eq!(bitenc.storage, [5122]);

        let mut bitenc = BitEnc::new();
        bitenc.push_bytes(&in_values, 2);
        // Contents should be 00000010
        let values: Vec<u8> = bitenc.iter().collect();
        assert_eq!(values, [2, 0]);
        assert_eq!(bitenc.storage, [2]);
    }

    #[test]
    fn test_from_dna_string() {
        let dna = "ACGTACGT";
        let bitenc = BitEnc::from_dna_string(dna);
        let values: Vec<u8> = bitenc.iter().collect();

        assert_eq!(bitenc.len, 8);
        assert_eq!(values, [0, 1, 2, 3, 0, 1, 2, 3]);

        let dna_cp = bitenc.to_dna_string();
        assert_eq!(dna, dna_cp);
    }

    #[test]
    fn test_pref() {
        let in_values: Vec<u8> = vec![2, 20];
        let mut bitenc = BitEnc::new();
        bitenc.push_bytes(&in_values, 8);
        // Contents should be 00010100 00000010

        let pref_bitenc = bitenc.pref(0);
        assert_eq!(pref_bitenc.storage.len(), 0);

        let pref_bitenc = bitenc.pref(8);
        assert_eq!(pref_bitenc, bitenc);

        let pref_bitenc = bitenc.pref(4);
        assert_eq!(pref_bitenc.storage, [2]);

        let pref_bitenc = bitenc.pref(6);
        assert_eq!(pref_bitenc.storage, [1026]);
        let values: Vec<u8> = pref_bitenc.iter().collect();
        assert_eq!(values, [2, 0, 0, 0, 0, 1]);

        bitenc.push_bytes(&in_values, 8);
        bitenc.push_bytes(&in_values, 8);

        let pref_bitenc = bitenc.pref(17);
        let values: Vec<u8> = pref_bitenc.iter().collect();
        assert_eq!(values, [2, 0, 0, 0, 0, 1, 1, 0, 2, 0, 0, 0, 0, 1, 1, 0, 2]);
    }

    #[test]
    fn test_suf() {
        let in_values: Vec<u8> = vec![2, 20];
        let mut bitenc = BitEnc::new();
        bitenc.push_bytes(&in_values, 8);
        // Contents should be 00010100 00000010

        let suf_bitenc = bitenc.suf(0);
        assert_eq!(suf_bitenc.storage.len(), 0);

        let suf_bitenc = bitenc.suf(8);
        assert_eq!(suf_bitenc, bitenc);

        let suf_bitenc = bitenc.suf(4);
        assert_eq!(suf_bitenc.storage, [20]);

        // 000101000000 64+256
        let suf_bitenc = bitenc.suf(6);
        assert_eq!(suf_bitenc.storage, [320]);
        let values: Vec<u8> = suf_bitenc.iter().collect();
        assert_eq!(values, [0, 0, 0, 1, 1, 0]);

        bitenc.push_bytes(&in_values, 8);
        bitenc.push_bytes(&in_values, 8);

        let suf_bitenc = bitenc.suf(17);
        let values: Vec<u8> = suf_bitenc.iter().collect();
        assert_eq!(values, [0, 2, 0, 0, 0, 0, 1, 1, 0, 2, 0, 0, 0, 0, 1, 1, 0]);
    }

    #[test]
    fn test_reverse() {
        let in_values: Vec<u8> = vec![2, 20];
        let mut bitenc = BitEnc::new();
        let rev_bitenc = bitenc.reverse();
        assert_eq!(bitenc, rev_bitenc);

        bitenc.push_bytes(&in_values, 8);
        // Contents should be 00010100 00000010

        let rev_bitenc = bitenc.reverse();
        let values: Vec<u8> = rev_bitenc.iter().collect();
        assert_eq!(values, [0, 1, 1, 0, 0, 0, 0, 2]);
    }

    #[test]
    fn test_kmers() {
        let dna = "TGCATTAGAAAACTCCTTGCCTGTCAGCCCGACAGGTAGAAACTCATTAATCCACACATTGA".to_string() +
                  "CTCTATTTCAGGTAAATATGACGTCAACTCCTGCATGTTGAAGGCAGTGAGTGGCTGAAACAGCATCAAGGCGTGAAGGC";
        let bitenc = BitEnc::from_dna_string(&dna);
        let kmers = bitenc.kmers();

        for i in 0..(dna.len() - kmer::K + 1) {
            assert_eq!(kmers[i].to_string(), &dna[i..(i + kmer::K)]);
        }

        let last_kmer = bitenc.last_kmer();
        assert_eq!(last_kmer.to_string(), &dna[(dna.len() - kmer::K)..]);

        for (idx, &k) in kmers.iter().enumerate() {
            assert_eq!(k, bitenc.get_kmer(idx));
        }

    }
}
