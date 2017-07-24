//
// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
//

//! Read FASTQs write to shardio buckets by barcode while counting barcodes. 
//! Reads without a correct barcode go to a special 'unbarcoded' bucket.
//! First-pass read subsampling is applied, and bc subsampling is applied.
//! Then go back over the unbarcoded bucket & correct barcodes, using the 
//!! BC counts.  Finally sort each bucket and write out FASTH files.

#![allow(dead_code)]

use std::path::{Path, PathBuf};
use std::collections::{BTreeMap, HashMap};
use std::io::{Read,  Write, BufWriter};
use std::hash::{Hash, Hasher};
use std::fs::rename;
use std::cmp::min;
use std::fs::File;
use std::boxed::Box;

use bincode;
use bincode::rustc_serialize::{encode_into, decode_from};

use utils::{load_barcode_whitelist2, load_barcode_whitelist_inverse};
use shardio::shard::{ShardWriteManager, ShardSender, Shardable, Serializer, ShardReaderSet};

use martian::{self, JsonDict, MartianStage};
use rustc_serialize::json::{self, Json, ParserError};
use rustc_serialize::Decodable;

use libc;
use lz4;
use flate2::write::GzEncoder;
use flate2::Compression;

use rayon::prelude::*;
use rayon;
use rayon::par_iter::reduce::{reduce, ReduceOp, SumOp};

use rand::XorShiftRng;
use rand::Rng;
use rand::SeedableRng;

use fastq_10x::*;
use fastq_10x::barcode::BarcodeValidator;

use fx;
use utils;


#[derive(RustcEncodable, RustcDecodable, PartialOrd, Ord, Eq, PartialEq)]
pub struct CrReadPair {
    gem_group: u8,
    corrected_bc: Option<u32>,

    read_group: String,
    bc_length: u8,
    trim_r1: u8,
    data: Vec<u8>,

    // Offsets of the 7 fields
    offsets: [u16; 8]
}

/// An efficient container for the components of a paired-end, barcoded Chromium Genome read.
/// We assume the the 10x barcode is at the beginning of R1
impl CrReadPair {
    pub fn fld(&self, field: usize) -> &[u8] {
        &self.data[self.offsets[field] as usize .. self.offsets[field+1] as usize]
    }

    /// FASTQ read header
    pub fn head(&self) -> &[u8] {
        self.fld(0)
    }

    /// Full raw R1 sequence
    pub fn r1_seq_raw(&self) -> &[u8] {
        self.fld(1)
    }

    /// Full raw R1 QVs
    pub fn r1_qual_raw(&self) -> &[u8] {
        self.fld(2)
    }

    /// Full R2 sequence
    pub fn r2_seq(&self) -> &[u8] {
        self.fld(3)
    }

    /// Full R2 QVs
    pub fn r2_qual(&self) -> &[u8] {
        self.fld(4)
    }

    /// Sample index (I1) sequence
    pub fn si_seq(&self) -> &[u8] {
        self.fld(5)
    }

    /// Sample index (I1) QVs
    pub fn si_qual(&self) -> &[u8] {
        self.fld(6)
    }

    /// Raw, uncorrected barcode sequence
    pub fn raw_bc_seq(&self) -> &[u8] {
        &self.r1_seq_raw()[0..self.bc_length as usize]
    }

    /// Raw barcode QVs
    pub fn raw_bc_qual(&self) -> &[u8] {
        &self.r1_qual_raw()[0..self.bc_length as usize]
    }

    /// Bases trimmed after the 10x BC, before the start of bases used from R1
    pub fn r1_trim_seq(&self) -> &[u8] {
        &self.r1_seq_raw()[self.bc_length as usize.. (self.bc_length+self.trim_r1) as usize]
    }

    /// QVs trimmed after the 10x BC, before the start of bases used from R1
    pub fn r1_trim_qual(&self) -> &[u8] {
        &self.r1_qual_raw()[self.bc_length as usize .. (self.bc_length+self.trim_r1) as usize]
    }

    /// Usable R1 bases after removal of BC and trimming
    pub fn r1_seq(&self) -> &[u8] {
        &self.r1_seq_raw()[(self.bc_length+self.trim_r1) as usize..]
    }

    /// Usable R1 bases after removal of BC and trimming
    pub fn r1_qual(&self) -> &[u8] {
        &self.r1_qual_raw()[(self.bc_length+self.trim_r1) as usize..]
    }

    /// Read group that read came from
    pub fn read_group(&self) -> &[u8] {
        self.read_group.as_bytes()
    }
}


impl Shardable for CrReadPair {
    fn shard(&self) -> usize {
        match self.corrected_bc {
            Some(bc) => bc as usize,
            None => {
                let mut h = fx::FxHasher::default();
                self.head().hash(&mut h);
                h.finish() as usize
            }
        }
    }
}

fn writeln<W: Write>(w: &mut W, bytes: &[u8]) {
    w.write(bytes).unwrap();
    w.write(b"\n").unwrap();
}

impl CrReadPair {
    
    pub fn write_to_fasth<W: Write>(&self, w: &mut W, inv_wl: &HashMap<u32, Vec<u8>>) {
        w.write(&self.head()).unwrap();
        w.write(b" ").unwrap();
        writeln(w, &self.read_group());

        writeln(w, &self.r1_seq());
        writeln(w, &self.r1_qual());

        writeln(w, &self.r2_seq());
        writeln(w, &self.r2_qual());

        // Special handling of barcode line.  Uncorrected is entered at bare sequence.
        // Corrected is entered as <corrected>-<gem_group>,raw
        match self.corrected_bc {
            Some(id) => {
                let ref bc = inv_wl[&id];
                w.write(&bc);
                write!(w, "-{},", self.gem_group);
                w.write(self.raw_bc_seq());
                w.write(b"\n");
            },
            None => writeln(w, &self.raw_bc_seq()),
        }
        writeln(w, &self.raw_bc_qual());

        writeln(w, &self.si_seq());
        writeln(w, &self.si_qual());
    }
}

#[derive(Clone)]
pub struct CrReadPairSer { }

impl Serializer<CrReadPair> for CrReadPairSer {
    fn serialize(&self, items: &Vec<CrReadPair>, buf: &mut Vec<u8>) {
        let mut encoder = lz4::EncoderBuilder::new().build(buf).unwrap();
        //let mut encoder = zstd::stream::Encoder::new(buf, 1).unwrap();
        encode_into(items, &mut encoder, bincode::SizeLimit::Infinite).unwrap();
        encoder.finish();
    }

    fn deserialize(&self, buf: &mut Vec<u8>, data: &mut Vec<CrReadPair>) {
        //let mut decoder = zstd::stream::Decoder::new(buf.as_slice()).unwrap();
        let mut decoder = lz4::Decoder::new(buf.as_slice()).unwrap();
        let r: Vec<CrReadPair> = decode_from(&mut decoder, bincode::SizeLimit::Infinite).unwrap();
        data.extend(r);
    }
}


struct FqProcessor<'a> {
    bc_length: usize,
    trim_length: usize,
    gem_group: u8,
    read_group: String,
    bc_map: &'a HashMap<Vec<u8>, u32>,
    bc_counts: HashMap<u32, u32>,
    no_bc_count: u32
}

impl<'a> FqProcessor<'a> {

    pub fn process_fastq<I: Iterator<Item=RawReadSet>>(&mut self, reads: I, sender: &mut ShardSender<CrReadPair>) {
        for raw_read in reads {
            let cr_read_pair = self.process_read(raw_read);

            match cr_read_pair.corrected_bc {
                Some(bc) => {
                    let entry = self.bc_counts.entry(bc).or_insert(0);
                    *entry += 1;
                },
                None => self.no_bc_count += 1,
            }

            sender.send(cr_read_pair);
        }
    }

    pub fn try_correct_bc(&self, raw_bc: &[u8]) -> Option<u32> {
        self.bc_map.get(raw_bc).map(|x| x.clone())
    }

    pub fn process_read(&self, raw: RawReadSet) -> CrReadPair {

        let (r1, r2, si) = raw;
        let mut data = Vec::new();
        let mut offsets = [0_u16;8];
        let mut pos = 0;

        {
            let mut add = |slc: &[u8]| {
                data.extend_from_slice(slc);
                offsets[pos+1] = data.len() as u16;
                pos += 1;
            };

            // FIXME -- trim to first space
            add(&r1.0); // head
            add(&r1.1); // r1 seq
            add(&r1.2); // r1 qual

            add(&r2.1); // r2 seq
            add(&r2.2); // r2 qual

            match si {
                Some(si_fq) => {
                    add(&si_fq.1); // si seq
                    add(&si_fq.2); // si qual
                },
                None => {
                    add(&[]);
                    add(&[])
                },
            };
        }

        let mut rp = 
            CrReadPair {
                data: data,
                offsets: offsets,
                bc_length: self.bc_length as u8,
                trim_r1: self.trim_length as u8,

                gem_group: self.gem_group,
                read_group: self.read_group.clone(),
                corrected_bc: None,
            };

        rp.corrected_bc = self.try_correct_bc(rp.raw_bc_seq());
        rp
    }
}


/// Returns a number that is unique to the calling thread.
///
/// Calling this function twice from the same thread will return the same
/// number. Calling this function from a different thread will return a
/// different number.
#[inline]
pub fn get_thread_id() -> usize {
    get_internal()
}

#[cfg(unix)]
#[inline]
fn get_internal() -> usize {
    unsafe { libc::pthread_self() as usize }
}


fn init_threads(num_threads: usize) {
    if num_threads > 0 {
        println!("setup rayon with {} threads", num_threads);
        let cfg = rayon::Configuration::new().set_num_threads(num_threads);
        rayon::initialize(cfg).expect("Couldn't initialize rayon");
    }
}


struct CountReduce { }

impl<K> ReduceOp<HashMap<K,u32>> for CountReduce where K: Hash + Eq {
    fn start_value(&self) -> HashMap<K,u32> {
        HashMap::new()
    }

    fn reduce(&self, mut v1: HashMap<K, u32>, mut v2: HashMap<K, u32>) -> HashMap<K, u32> {
        for (k,v) in v2.drain() {
            let counter = v1.entry(k).or_insert(0);
            *counter += v;
        }

        v1
    }
}

/// Entry point for an MSP run. bc_wl is the path to the BC whitelist file, which will be used to
/// convert BC sequences to ints. fastq is set of pseudo-FASTQ files created by the _SORT_FASTQ_BY_BARCODE
/// pipeline. Permutation is the permutation of p-mers to use, serialized into the give path.
/// The best permutation to use is determined by msp_get_pmer_permutation below.

// Count 'passed' reads
// Subsample at requested rate
// BC subsample at requested rate

const NUM_READS_BUFFER: usize = 2048;

#[inline(never)]
pub fn main_bucket_bcs(
    bc_wl_path: &Path, 
    fastq: Vec<FqInputChunk>,
    read_path: &Path, 
    no_bc_read_path: &Path,
    bc_count_path: &Path,
    trim_length: usize,
    num_threads: usize) {

    init_threads(num_threads);
    let ref bc_wl = load_barcode_whitelist2(bc_wl_path);

    let ref bc_shards = ShardWriteManager::new(read_path, 2048, NUM_READS_BUFFER, 2, CrReadPairSer{});
    let ref no_bc_shards = ShardWriteManager::new(no_bc_read_path, 2048, NUM_READS_BUFFER, 1, CrReadPairSer{});

    let fq_w_senders: Vec<(FqInputChunk, ShardSender<CrReadPair>, ShardSender<CrReadPair>)> = 
    fastq.into_iter().map(|fq| (fq, bc_shards.get_sender(), no_bc_shards.get_sender())).collect();

    let par_iter = fq_w_senders.into_par_iter().map(|(fq, bc_sender, no_bc_sender)| {
        info!("start proc: {}, on thread: {}", fq.read1, get_thread_id()); 
        let counts = process_fastq(bc_wl, fq.clone(), trim_length, bc_sender, no_bc_sender);
        info!("end proc: {}", fq.read1);
        counts
    });

    let bc_counts = reduce(par_iter, &CountReduce{});
    utils::write_obj(&bc_counts, bc_count_path).expect("error writing BC counts");
}


pub fn main_sort_bcs(
    read_subsample_rate: Option<f64>,
    read_shards: Vec<String>,
    inv_whitelist: &HashMap<u32, Vec<u8>>,
    chunk_id: usize,
    total_chunks: usize,
    num_output_files: usize,
    out_dir: &Path,
    num_threads: usize) -> (Vec<PathBuf>, usize) {

    // output_file_idx % num_chunks == chunk_id --> set of output files to be generated by this chunk
    // shard % num_output_files == output_file_idx --> set of reads in output_file_idx
    
    init_threads(num_threads);
    let read_subsample_rate = read_subsample_rate.unwrap_or(1.0);

    // BC correction is done prior to this step & folded into read_shards
    let reader = ShardReaderSet::open(&read_shards,CrReadPairSer{});

    let my_outputs = (0..num_output_files).filter(|x| x % total_chunks == chunk_id);
    let my_output_files: Vec<(usize, PathBuf)> = my_outputs.map(|x| (x, out_dir.join(format!("chunk-{}.fasth.gz", x)))).collect();

    let par_iter = my_output_files.clone().into_par_iter().map(|(output_file_idx, out_path)| {

        // Open FASTH writer
        let wa = File::create(&out_path).unwrap();
        let wb = GzEncoder::new(wa, Compression::Fast);
        let mut w = BufWriter::new(wb);

        info!("Outputing file {}: {:?}", output_file_idx, out_path);

        // Seed based on the shard ID.  Still deterministic, but avoid correlations in sampling between shards
        let mut rand = XorShiftRng::from_seed([0,0,1,output_file_idx as u32]);

        let mut reads = Vec::new();
        let mut emitted_reads = 0;

        for shard in 0 .. reader.num_shards() {
            if shard % num_output_files == output_file_idx {
                reads.clear();
                reader.read_shard(shard, &mut reads);
                reads.sort();

                for r in reads.iter() {
                    // Final layer of read subsampling, to get within Poisson sample of requested number of reads
                    if rand.next_f64() < read_subsample_rate {
                        r.write_to_fasth(&mut w, inv_whitelist);
                        emitted_reads += 1;
                    }
                }
            }
        }

        emitted_reads
    });

    let total_reads = reduce(par_iter, &SumOp);
    (my_output_files.into_iter().map(|(_,f)| f).collect(), total_reads)
}

pub fn process_fastq<'a>(
    wl: &'a HashMap<Vec<u8>, u32>, 
    chunk: FqInputChunk,
    trim_length: usize,
    mut bc_sender: ShardSender<CrReadPair>,
    mut no_bc_sender: ShardSender<CrReadPair>) -> HashMap<(u8, u32), u32> {

    // Counter for BCs
    let mut counts = HashMap::new();

    // Always subsample deterministically
    let mut rand = XorShiftRng::new_unseeded();
    // below: 1.0-r "inverts" the sense of the fraction so that values near
    // 1.0 don't overflow.  Values near zero are a problem.
    let bc_subsample_thresh = chunk.bc_subsample_rate.map_or(0,
                                |r| {   assert!( r > 0.0, "zero barcode fraction passed in");
                                        ((1.0-r) * u64::max_value() as f64) as u64 } );
    let read_subsample_rate = chunk.subsample_rate.unwrap_or(1.0);

    let fq_proc = FqProcessor {
        bc_length: chunk.bc_length,
        trim_length: trim_length,
        gem_group: chunk.gem_group as u8,
        read_group: chunk.read_group,
        bc_map: wl,
        bc_counts: HashMap::new(),
        no_bc_count: 0,
    };

    let iter: Box<Iterator<Item=RawReadSet>> =
        match chunk.read2 {
            Some(r2) => open_fastq_pair_iter(chunk.read1, r2, chunk.sample_index),
            None => open_interleaved_fastq_pair_iter(chunk.read1, chunk.sample_index),
        };

    for r in iter {

        // Read subsampling
        if rand.next_f64() < read_subsample_rate {

            let read_pair = fq_proc.process_read(r);
            match read_pair.corrected_bc {
                Some(bc) => {
                    // barcode subsampling
                    if fx::hash(&bc) >= bc_subsample_thresh {
                        // count reads on each (gem_group, bc)
                        let c = counts.entry((read_pair.gem_group, bc)).or_insert(0);
                        *c += 1;
                        bc_sender.send(read_pair)
                    }
                },

                
                None => {
                    // don't bc subsample these reads yet -- that will
                    // happen in the 2nd pass in 'process_unbarcoded'
                    // count non-barcoded reads in the 0 bin
                    let c = counts.entry((read_pair.gem_group, 0)).or_insert(0);
                    *c += 1;
                    no_bc_sender.send(read_pair)
                },
            }
        }
    }

    counts
}

pub fn process_unbarcoded<'a>(
    bc_corrected_reads_fn: &str,
    unbarcoded_read_shards: Vec<String>,
    barcode_validator: &BarcodeValidator<'a>,
    bc_subsample_rate: Option<f64>) -> HashMap<(u8, u32), u32> {

    let ref writer = ShardWriteManager::new(Path::new(bc_corrected_reads_fn), 2048, NUM_READS_BUFFER, 2, CrReadPairSer{});

    // reader
    let reader = ShardReaderSet::open(&unbarcoded_read_shards, CrReadPairSer{});
    // below: 1.0-r "inverts" the sense of the fraction so that values near
    // 1.0 don't overflow.  Values near zero are a problem.
    let bc_subsample_thresh = bc_subsample_rate.map_or(0,
                                |r| {   assert!( r > 0.0, "zero barcode fraction passed in");
                                        ((1.0-r) * u64::max_value() as f64) as u64 } );

    let iter = 
        (0..reader.num_shards()).
        map(|x| (x, writer.get_sender())).collect::<Vec<_>>().
        into_par_iter().
        map(|(x, mut read_sender)| {

        // Counter for BCs
        let mut counts = HashMap::new();
        let mut data = Vec::new(); 
        reader.read_shard(x, &mut data); 
        data.sort(); // need sort to preserve determinism

        for mut read_pair in data {

            // Correct the BC, if possible
            read_pair.corrected_bc = 
                barcode_validator.correct_barcode(
                    read_pair.gem_group, 
                    read_pair.raw_bc_seq(), 
                    read_pair.raw_bc_qual());

            // Read subsampling has already occured in process_fastq -- don't repeat here
            match read_pair.corrected_bc {
                Some(bc) => {
                    // barcode subsampling
                    if fx::hash(&bc) >= bc_subsample_thresh {
                        // count reads on each (gem_group, bc)
                        let c = counts.entry((read_pair.gem_group, bc)).or_insert(0);
                        *c += 1;
                        read_sender.send(read_pair)
                    }
                },

                None => {
                    // apply the barcode subsample factor to non-barcoded reads, to keep the balance correct
                    if fx::hash(&read_pair.head()) >= bc_subsample_thresh {
                        let c = counts.entry((read_pair.gem_group, 0)).or_insert(0);
                        *c += 1;
                        read_sender.send(read_pair)
                    }
                },
            }
        }

        counts
    });
    
    reduce(iter, &CountReduce{})
}


#[derive(RustcDecodable, RustcEncodable, Debug, Clone)]
pub struct FqInputChunk {
    /// Does the barcode read need to be RC'd
    barcode_reverse_complement: bool,

    /// Which read is the BC found in?
    bc_in_read: usize,

    /// Length of the barcode
    bc_length: usize,

    /// GEM group of the read
    gem_group: usize,

    /// FASTQ file with R1, or interleaved R1 & R2 is read2 is None
    read1: String,

    /// FASTQ file with R2, unless interleaved with read1
    read2: Option<String>,

    /// Read group of this set of reads
    read_group: String,

    /// FASTQ file with SI reads
    sample_index: Option<String>,

    /// subsample rate
    subsample_rate: Option<f64>,

    /// bc_subsample rate
    bc_subsample_rate: Option<f64>,
}


fn read_json(path: &str) -> Result<Json, ParserError> {
    let mut f = try!(File::open(path));
    let mut buf = String::new();
    try!(f.read_to_string(&mut buf));
    Json::from_str(&buf)
}

pub fn go(bc_whitelist: &Path, chunk_args: &str, output: &Path, bc_counts: &Path, num_threads: usize)  {
        let args = read_json(chunk_args).unwrap();

        let mut decoder = json::Decoder::new(args["chunks"].clone());
        let chunks: Vec<FqInputChunk> = Decodable::decode(&mut decoder).unwrap();

        for f in chunks.iter() {
            println!("proc: {}", f.read1);
        }

        //[main_bucket_bcs(bc_whitelist, chunks, output, bc_counts, 0, num_threads);
    }

#[derive(Default)]
pub struct BucketBcsMartian;

impl MartianStage for BucketBcsMartian {
    fn split(&self, args: JsonDict) -> JsonDict {

        let mut decoder = json::Decoder::new(args["chunks"].clone());
        let chunks: Vec<FqInputChunk> = Decodable::decode(&mut decoder).unwrap();

        let mut stage_defs: Vec<Json> = Vec::new();
        for (idx, _) in chunks.iter().enumerate() {

            let mut bt = BTreeMap::new();
            bt.insert("__mem_gb".to_string(), Json::F64(8.0));
            bt.insert("__threads".to_string(), Json::I64(3));
            bt.insert("which".to_string(), Json::I64(idx as i64));
            stage_defs.push(Json::Object(bt));
        }

        let join_def: Json = Json::from_str( r#"
                {
                  "__mem_gb": 12.0,
                  "__threads": 4
                }
        "# ).unwrap();

        let mut cc = BTreeMap::new();
        cc.insert("chunks".to_string(), Json::Array(stage_defs));
        cc.insert("join".to_string(), join_def );
        cc
    }

    fn main(&self, args: JsonDict, outs: JsonDict) -> JsonDict {

        let mut decoder = json::Decoder::new(args["chunks"].clone());
        let chunks: Vec<FqInputChunk> = Decodable::decode(&mut decoder).unwrap();

        let which_chunk = args["which"].as_i64().unwrap() as usize;
        let chunk = chunks[which_chunk].clone();

        let nthreads = args["__threads"].as_i64().unwrap() as usize;
        let trim_length = args["trim_length"].as_i64().unwrap() as usize;

        let _whitelist = args["barcode_whitelist_path"].as_string().unwrap();
        let whitelist  = Path::new(_whitelist);

        let read_fn = Path::new({outs.get("read_buckets").unwrap().as_string().unwrap().clone()});
        let no_bc_read_fn = Path::new({outs.get("no_bc_read_buckets").unwrap().as_string().unwrap().clone()});
 
        let _bc_count_fn = {outs.get("bc_counts").unwrap().as_string().unwrap().clone()};
        let bc_count_fn = Path::new(_bc_count_fn);

        main_bucket_bcs(&whitelist, vec![chunk], &read_fn, &no_bc_read_fn, &bc_count_fn, trim_length, nthreads);
        outs.clone()
    }

    fn join(&self, args: JsonDict, mut outs: JsonDict, _: Vec<JsonDict>, chunk_outs: Vec<JsonDict>) -> JsonDict {

            let mut read_buckets: Vec<Json> = Vec::new();
            let mut no_bc_read_buckets: Vec<String> = Vec::new();

            let cr = CountReduce{};
            let mut bc_counts = cr.start_value();            

            for c in chunk_outs {
                // This is the output data -- pass to outs
                read_buckets.push(c.get("read_buckets").unwrap().clone());

                // These are the unbarcoded reads -- try to correct barcodes and writer to a new read bucket
                no_bc_read_buckets.push(c["no_bc_read_buckets"].as_string().unwrap().to_string());

                // Load and aggregate the BC counts from this chunk
                let chunk_counts: HashMap<(u8, u32), u32> = utils::read_obj(c["bc_counts"].as_string().unwrap()).unwrap();
                bc_counts = cr.reduce(bc_counts, chunk_counts);
            }

            let mut bc_reads = 0;
            let mut no_bc_reads = 0;
            for (&(_, bc),v) in bc_counts.iter() {
                if bc == 0 { no_bc_reads += *v } else { bc_reads += *v }
            }
            info!("First pass: {} barcoded reads, {} unbarcoded reads", bc_reads, no_bc_reads);

            // Output the corrected read bucket
            let bc_corrected_read_bucket_fn = outs["read_buckets"].as_string().unwrap().to_string();
            read_buckets.push(Json::String(bc_corrected_read_bucket_fn.clone()));

            let bc_wl_path = args["barcode_whitelist_path"].as_string().unwrap();
            let bc_wl = load_barcode_whitelist2(bc_wl_path);

            // Get the subsample rates from the first input chunk (they are assumed to be constant)
            let mut decoder = json::Decoder::new(args["chunks"].clone());
            let chunks: Vec<FqInputChunk> = Decodable::decode(&mut decoder).unwrap();

            let bc_subsample_rate = chunks[0].bc_subsample_rate;

            // Correct barcodes & count resulting reads
            let bc_corrected_counts =  {
                let validator = BarcodeValidator {
                    max_expected_barcode_errors: args["max_expected_barcode_errors"].as_f64().unwrap(),
                    bc_confidence_threshold: args["bc_confidence_threshold"].as_f64().unwrap(),
                    whitelist: &bc_wl,
                    bc_counts: &bc_counts,
                };

                process_unbarcoded(
                    &bc_corrected_read_bucket_fn,
                    no_bc_read_buckets,
                    &validator,
                    bc_subsample_rate)
            };

            let mut corr_bc_reads = 0;
            let mut corr_no_bc_reads = 0;
            for (&(_, bc),v) in bc_corrected_counts.iter() {
                if bc == 0 { corr_no_bc_reads += *v } else { corr_bc_reads += *v }
            }
            info!("Correction of unbarcoded: {} barcoded reads, {} unbarcoded reads", corr_bc_reads, corr_no_bc_reads);

            // The unbarcoded reads were all counted in (gem_group, 0). Remove these counts, then add the new ones
            for gg in 0_u8 .. 255_u8 { bc_counts.remove(&(gg,0)); }
            bc_counts = cr.reduce(bc_counts, bc_corrected_counts);            

            utils::write_obj(&bc_counts, outs["bc_counts"].as_string().unwrap());


            let mut total_reads = 0_i64;
            for v in bc_counts.values() {
                total_reads += *v as i64;
            }

            let requested_read_pairs = args["requested_read_pairs"].as_f64().map(|x| x as i64);
            let final_subsample_rate = requested_read_pairs.map(|r| r as f64 / total_reads as f64);

            info!("Total read pairs: {}, Requested read pairs: {:?}, Final subsample_rate: {:?}", total_reads, requested_read_pairs, final_subsample_rate);

            match final_subsample_rate {
                Some(f) => outs.insert("final_subsample_rate".to_string(), Json::F64(f)),
                None => outs.insert("final_subsample_rate".to_string(), Json::Null)
            };

            outs.insert("total_reads".to_string(), Json::I64(total_reads));
            outs.insert("read_buckets".to_string(), Json::Array(read_buckets));
            outs.insert("no_bc_read_buckets".to_string(), Json::Array(vec![]));
            outs
    }
}


pub struct SortBcsMartian;


#[derive(RustcDecodable, RustcEncodable, Debug)]
struct SortArgs {
    read_buckets: Vec<String>,
    subsample_rate: Option<f64>,
    barcode_whitelist_path: String,
}

#[derive(RustcDecodable, RustcEncodable, Debug)]
struct SortSplit {
    __threads: usize,
    __mem_gb: f64,
    chunk_id: usize,
    total_chunks: usize,
}

#[derive(RustcDecodable, RustcEncodable, Debug)]
struct SortOuts {
    reads: Vec<String>,
    total_reads: f64,
}


const NUM_OUT_FILES: usize = 128;
const CHUNK_THREADS: usize = 4;

impl MartianStage for SortBcsMartian {
    fn split(&self, _: JsonDict) -> JsonDict {
        
        let mut stage_defs = Vec::new();

        let nchunks = 16;
        for i in 0..nchunks {
            let mut chnk = JsonDict::new();
            chnk.insert("__threads".to_string(), Json::I64(3));
            chnk.insert("__mem_gb".to_string(), Json::F64(8.0));
            chnk.insert("chunk_id".to_string(), Json::I64(i));
            chnk.insert("total_chunks".to_string(), Json::I64(nchunks));
            stage_defs.push(Json::Object(chnk));

        }
        
        let mut defs = JsonDict::new();
        defs.insert("chunks".to_string(), Json::Array(stage_defs));
        defs
    }

    fn main(&self, json_args: JsonDict, outs: JsonDict) -> JsonDict {
        // Martian packs both the original args and the chunk-specific split args
        // into the args json. Pull them out sepeartely.

        let args: SortArgs = martian::obj_decode(json_args.clone());
        let split: SortSplit = martian::obj_decode(json_args);

        // Load id -> bc seq map
        let inverse_whitelist = load_barcode_whitelist_inverse(args.barcode_whitelist_path);

        let out_path = outs["reads"].as_string().unwrap().clone();
        let mut path = PathBuf::from(out_path);
        path.pop(); // This is where we will write to

        let (out_files, total_reads) = main_sort_bcs(
            args.subsample_rate, args.read_buckets,
            &inverse_whitelist,
            split.chunk_id, split.total_chunks, 
            NUM_OUT_FILES,  &path, CHUNK_THREADS);
        
        let out_strings = out_files.into_iter().map(|x| x.to_str().unwrap().to_string()).collect();

        martian::encode_to_json(&SortOuts{ reads: out_strings, total_reads: total_reads as f64 }).as_object().unwrap().clone()
    }

    fn join(&self, _: JsonDict, outs: JsonDict, chunk_defs: Vec<JsonDict>, chunk_outs: Vec<JsonDict>) -> JsonDict {

        let _path = outs["reads"].as_string().unwrap().clone();
        let mut out_path = PathBuf::from(_path);
        out_path.pop();

        let mut out_files = Vec::new();

        let mut total_reads = 0.0;
        for (_, _chnk_out) in chunk_defs.into_iter().zip(chunk_outs.into_iter()) {
            let chnk_out: SortOuts = martian::obj_decode(_chnk_out);
            total_reads += chnk_out.total_reads;

            for _src_path in chnk_out.reads {
                let src_path = PathBuf::from(_src_path);
                let src_name = src_path.file_name().unwrap();
                let dest_path = out_path.join(src_name);

                rename(&src_path, &dest_path).unwrap();
                out_files.push(dest_path.to_str().unwrap().to_string());
            }
        }

        info!("Total read pairs emitted: {}", total_reads as u64);        
        let new_outs = SortOuts { reads: out_files, total_reads: total_reads };
        martian::encode_to_json(&new_outs).as_object().unwrap().clone()
    }
}
