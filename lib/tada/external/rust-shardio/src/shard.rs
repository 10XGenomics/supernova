//
// Copyright (c) 2017 10x Genomics, Inc. All rights reserved.
//

use std::fs::File;
use std::io::{Result, Error, Seek, SeekFrom};
use std::os::unix::io::{RawFd, AsRawFd};
use std::collections::HashMap;

use std::sync::{Arc, Mutex};
use std::sync::mpsc::sync_channel;
use std::sync::mpsc::{SyncSender, Receiver};
use std::path::{Path};

use std::thread;
use std::thread::JoinHandle;
use std::marker::PhantomData;

use bincode::rustc_serialize::{encode_into, decode_from};
use bincode;
use byteorder::{BigEndian, WriteBytesExt, ReadBytesExt};

use libc::{pread, pwrite, c_void, off_t, size_t, ssize_t};

fn err(e: ssize_t) -> Result<usize> {
    if e == -1 as ssize_t {
        Err(Error::last_os_error())
    } else {
        Ok(e as usize)
    }
}


fn read_at(fd: &RawFd, pos: u64, buf: &mut [u8]) -> Result<usize> {
    err(unsafe {
        pread(*fd, buf.as_mut_ptr() as *mut c_void, buf.len() as size_t, pos as off_t)
    })
}


fn write_at(fd: &RawFd, pos: u64, buf: &[u8]) -> Result<usize> {
    err(unsafe {
        pwrite(*fd, buf.as_ptr() as *const c_void, buf.len() as size_t, pos as off_t)
    })
}


pub trait Serializer<T>: Clone + Send {
    fn serialize(&self, items: &Vec<T>, buf: &mut Vec<u8>);
    fn deserialize(&self, buf: &mut Vec<u8>, data: &mut Vec<T>);
}


#[derive(Clone, RustcEncodable, RustcDecodable, Debug, PartialEq, Eq)]
struct  ShardRecord {
    offset: usize,
    shard: usize,
    block_size: usize,
    n_items: usize
}


pub struct FileRegionState {
    // Current start position of next chunk
    cursor: usize,

    // Record of chunks written
    regions: Vec<ShardRecord>,
}

pub struct FileRegionManager {
    state: Mutex<FileRegionState>
}

impl FileRegionManager {
    pub fn new() -> FileRegionManager {
        let state = FileRegionState {
            cursor: 4096,
            regions: Vec::new(),
        };

        FileRegionManager {
            state: Mutex::new(state)
        }
    }

    pub fn register_write(&self, bucket: usize, block_size: usize, n_items: usize) -> usize
    {
        let mut state = self.state.lock().unwrap();
        let cur_offset = state.cursor;
        let reg = 
            ShardRecord {
                offset: state.cursor,
                shard: bucket,
                block_size: block_size,
                n_items: n_items };

        state.regions.push(reg);
        state.cursor = state.cursor + block_size;

        cur_offset
    }
}

struct FileChunkWriter {
    file: RawFd,
    region_manager: Arc<FileRegionManager>
}

impl FileChunkWriter {
    fn write(&mut self, bucket: usize, n_items: usize, data: &Vec<u8>) {
        let offset = self.region_manager.register_write(bucket, data.len(), n_items);
        write_at(&mut self.file, offset as u64, data).unwrap();
    }
}


pub trait Shardable {
    fn shard(&self) -> usize;
}

struct ShardWriterThread<T, S> where T: Sync + Send + Shardable, S: Serializer<T> {
    thread_id: usize,
    total_shards: usize,
    thread_bits: usize,
    serializer: S,

    writer: FileChunkWriter,
    rx: Receiver<Option<Vec<T>>>,
    item_buffers: Vec<Vec<T>>,
    write_buffer: Vec<u8>,
}

impl<T, S> ShardWriterThread<T, S> where T: 'static + Sync + Send + Shardable + , S: 'static + Serializer<T> {
    fn new(
        buffer_size: usize,
        thread_id: usize,
        total_shards: usize,
        thread_bits: usize,
        serializer: S,
        writer: FileChunkWriter,
        rx: Receiver<Option<Vec<T>>>) -> ShardWriterThread<T, S>{

        let local_shards = total_shards >> thread_bits;

        let mut item_buffers = Vec::new();
        for _ in 0 .. local_shards {
            item_buffers.push(Vec::with_capacity(buffer_size))
        }

        ShardWriterThread {
            thread_id: thread_id,
            total_shards: total_shards,
            thread_bits: thread_bits,
            serializer: serializer,

            writer: writer,
            rx: rx,
            item_buffers: item_buffers,
            write_buffer: Vec::new(),
        }
    }

    fn process(&mut self) {

        loop {
            match self.rx.recv() {
                Ok(Some(v)) => {
                    for item in v {
                        self.add(item);
                    }
                },

                Ok(None) => {
                    self.flush();
                    break;
                },

                Err(_) => break,
            }
        }
    }

    fn add(&mut self, item: T) {
        let key = item.shard();

        let shard = key % self.total_shards;
        let local_shard = shard >> self.thread_bits;

        let write = {
            let buf = self.item_buffers.get_mut(local_shard).unwrap();
            buf.push(item);
            buf.len() == buf.capacity()
        };

        if write
        {
            self.write_item_buffer(local_shard);
        }
    }

    fn write_item_buffer(&mut self, local_shard: usize)
    {
        let main_shard = local_shard << self.thread_bits | (self.thread_id as usize);
        let items = self.item_buffers.get_mut(local_shard).unwrap();

        self.write_buffer.clear();
        self.serializer.serialize(items, &mut self.write_buffer);
        self.writer.write(main_shard, items.len(), &self.write_buffer);
        items.clear();
    }

    fn flush(&mut self)
    {
        for i in 0 .. self.item_buffers.len()
        {
            if self.item_buffers[i].len() > 0 {
                self.write_item_buffer(i)
            }
        }
    }
}


pub struct ShardWriteManager<T: 'static + Sync + Send + Shardable> {
    total_shards: usize,
    handles: Vec<JoinHandle<()>>,
    txs: Vec<SyncSender<Option<Vec<T>>>>,
    region_manager: Arc<FileRegionManager>,
    file: File,
}

impl<T> ShardWriteManager<T> where T: 'static + Sync + Send + Shardable {
    pub fn new<S: 'static + Serializer<T>>(path: &Path, per_shard_buffer_size: usize, num_shards: usize, thread_bits: usize, serializer: S) -> ShardWriteManager<T> {
        let mut txs = Vec::new();
        let mut handles = Vec::new();

        let file = File::create(path).unwrap();

        let regions = FileRegionManager::new();
        let arc_regions = Arc::new(regions);
        let num_threads = 1 << thread_bits;

        for thread_id in (0..num_threads).into_iter() {
            // Create communication channel
            let (tx, rx) = sync_channel::<Option<Vec<T>>>(10);
            txs.push(tx);

            // Copy the file handle
            let writer = FileChunkWriter {
                file: file.as_raw_fd(),
                region_manager: arc_regions.clone(),
            };

            let mut thread = ShardWriterThread::new(
                per_shard_buffer_size,
                thread_id,
                num_shards,
                thread_bits,
                serializer.clone(),
                writer,
                rx);

            let handle = thread::spawn(move || { thread.process() });
            handles.push(handle);
        }

        ShardWriteManager {
            total_shards: num_shards,
            region_manager: arc_regions,
            handles: handles,
            txs: txs,
            file: file,
        }
    }

    pub fn num_threads(&self) -> usize {
        self.txs.len()
    }

    pub fn get_sender(&self) -> ShardSender<T>
    {
        let sender: ShardSender<T> = ShardSender::new(&self);
        sender
    }

    /// Write out the shard positioning data
    fn write_index_block(&mut self) {

        let ref _regs = *self.region_manager;
        let regs = _regs.state.lock().unwrap();
        let mut buf = Vec::new();

        encode_into(&regs.regions, &mut buf, bincode::SizeLimit::Infinite).unwrap();

        let index_block_position = regs.cursor;
        let index_block_size = buf.len();

        write_at(&self.file.as_raw_fd(), index_block_position as u64, buf.as_slice()).unwrap();

        self.file.seek(SeekFrom::Start((index_block_position + index_block_size) as u64)).unwrap();
        self.file.write_u64::<BigEndian>(self.total_shards as u64).unwrap();
        self.file.write_u64::<BigEndian>(index_block_position as u64).unwrap();
        self.file.write_u64::<BigEndian>(index_block_size as u64).unwrap();
    }

    /// Shutdown the writing
    fn finish(&mut self) {
        for tx in self.txs.iter() {
            tx.send(None).unwrap();
        }

        for t in self.handles.drain(..) {
            t.join().unwrap();
        }

        self.write_index_block();
    }
}


impl<T> Drop for ShardWriteManager<T>  where T: 'static + Sync + Send + Shardable
{
    fn drop(&mut self) {
        self.finish();
    }
}



pub struct ShardSender<T: Sync + Send + Shardable> {
    tx_channels: Vec<SyncSender<Option<Vec<T>>>>,
    buffers: Vec<Vec<T>>,
    buf_size: usize,
    thread_shards: usize,
}

impl<T: Sync + Send + Shardable> ShardSender<T> {
    fn new(manager: &ShardWriteManager<T>) -> ShardSender<T> {
        let mut new_txs = Vec::new();
        for t in manager.txs.iter() {
            new_txs.push(t.clone())
        }

        let n = manager.txs.len();

        let mut buffers = Vec::with_capacity(n);
        for _ in 0..n {
            buffers.push(Vec::with_capacity(256));
        }

        ShardSender{
            tx_channels: new_txs,
            buffers: buffers,
            buf_size: 256,
            thread_shards: n,
        }
    }

    pub fn send(&mut self, item: T) {
        let shard_idx = item.shard() % self.thread_shards;
        let send = {
            let buf = self.buffers.get_mut(shard_idx).unwrap();
            buf.push(item);
            buf.len() == self.buf_size
        };

        if send {
            self.buffers.push(Vec::with_capacity(self.buf_size));
            let send_buf = self.buffers.swap_remove(shard_idx);

            let out_ch = self.tx_channels.get(shard_idx).unwrap();
            out_ch.send(Some(send_buf)).unwrap();
        }
    }

    pub fn finished(&mut self) {
        for (idx, buf) in self.buffers.drain(..).enumerate() {
            let out_ch = self.tx_channels.get(idx).unwrap();
            out_ch.send(Some(buf)).unwrap();
        }
    }
}

impl<T: Sync + Send + Shardable> Drop for ShardSender<T> {
    fn drop(&mut self) {
        self.finished();
    }
}

pub struct ShardReader<T, S> {
    serializer: S,
    file: File,
    num_shards: usize,
    index: HashMap<usize, Vec<ShardRecord>>,
    phantom: PhantomData<T>,
}

impl<T, S> ShardReader<T, S> where S: Serializer<T> {
    pub fn open<P: AsRef<Path>>(path: P, serializer: S) -> ShardReader<T, S> {
        let mut f = File::open(path).unwrap();

        let (num_shards, index_rows) = Self::read_index_block(&mut f);
        let mut index: HashMap<usize, Vec<ShardRecord>> = HashMap::new();

        for rec in index_rows {
            let shard_recs = index.entry(rec.shard).or_insert_with(|| Vec::new());
            shard_recs.push(rec);
        }

        ShardReader {
            serializer: serializer,
            file: f,
            num_shards: num_shards,
            index: index,
            phantom: PhantomData,
        }
    }

    /// Write out the shard positioning data
    fn read_index_block(file: &mut File) -> (usize, Vec<ShardRecord>) {
    
        let _ = file.seek(SeekFrom::End(-24)).unwrap();
        let num_shards = file.read_u64::<BigEndian>().unwrap() as usize;
        let index_block_position = file.read_u64::<BigEndian>().unwrap();
        let _ = file.read_u64::<BigEndian>().unwrap();
        file.seek(SeekFrom::Start(index_block_position as u64)).unwrap();
        let regs = decode_from::<File, Vec<ShardRecord>>(file, bincode::SizeLimit::Infinite).unwrap();

        (num_shards, regs)
    }

    pub fn read_shard_buf(&self, shard: usize, data: &mut Vec<T>, buf: &mut Vec<u8>) {
        match self.index.get(&shard) {
            Some(recs) => {
                for rec in recs.iter() {
                    buf.resize(rec.block_size, 0);
                    let read_len = read_at(&self.file.as_raw_fd(), rec.offset as u64, buf.as_mut_slice()).unwrap();
                    assert_eq!(read_len, rec.block_size);
                    self.serializer.deserialize(buf, data);
                }
            },

            None => (),
        }
    }

    pub fn read_shard(&self, shard: usize) -> Vec<T> {
        let mut buf = Vec::new();
        let mut data = Vec::new();

        self.read_shard_buf(shard, &mut data, &mut buf);
        data
    }

    pub fn shard_len(&self, shard: usize) -> usize {
        match self.index.get(&shard) {
            Some(shard_idx) => {
                shard_idx.iter().map(|x| x.n_items).sum()
            }, 
            None => 0,
        }
    }

    pub fn num_shards(&self) -> usize {
        *self.index.keys().max().unwrap_or(&0) + 1
    }
}

pub struct ShardReaderSet<T,S> {
    readers: Vec<ShardReader<T,S>>
}

impl<T, S> ShardReaderSet<T, S> where S: Serializer<T> {
    pub fn open<P: AsRef<Path>>(shard_files: &Vec<P>, serializer: S) -> ShardReaderSet<T,S> {
        let mut readers = Vec::new();

        for p in shard_files {
            let reader = ShardReader::open(p, serializer.clone());
            readers.push(reader);
        }

        ShardReaderSet{ 
            readers: readers
        }
    }

    pub fn read_shard(&self, shard: usize, data: &mut Vec<T>) {
        let mut buf = Vec::new();
        for r in self.readers.iter() {
            r.read_shard_buf(shard, data, &mut buf)
        }
    }

    pub fn num_shards(&self) -> usize {
        self.readers[0].num_shards
    }

    pub fn shard_len(&self, shard: usize) -> usize {
        self.readers.iter().map(|r| r.shard_len(shard)).sum()
    }
}



#[cfg(test)]
mod shard_tests {
    use bincode::rustc_serialize::{encode_into, decode};
    use bincode;
    use tempfile;

    use super::*;

    #[derive(Copy, Clone, Eq, PartialEq, RustcEncodable, RustcDecodable, Debug)]
    struct T1 {
        a: u64,
        b: u32,
        c: u16,
        d: u8,
    }

    impl Shardable for T1 {
        fn shard(&self) -> usize {
            self.a as usize
        }
    }

    #[derive(Clone)]
    struct T1Ser {}

    impl Serializer<T1> for T1Ser {
        fn serialize(&self, items: &Vec<T1>, buf: &mut Vec<u8>) {
            encode_into(items, buf, bincode::SizeLimit::Infinite).unwrap();
        }

        fn deserialize(&self, buf: &mut Vec<u8>, data: &mut Vec<T1>) {
            let mut buf_slice = buf.as_mut_slice();
            let r: Vec<T1> = decode(&mut buf_slice).unwrap();
            data.extend(r);
        }
    }

    #[test]
    fn test_shard_round_trip() {

        let ns = 16;
        // Test different numbers of threads
        check_round_trip(1024, ns, 0, 1<<4);
        check_round_trip(4096, ns, 0, 1<<8);
        check_round_trip(128, ns, 0, 1<<12);
        check_round_trip(50, ns, 0, 1<<16);


        // Test different numbers of threads
        check_round_trip(100, 1024, 0, 2<<16);
        check_round_trip(200, 1024, 1, 2<<16);
        check_round_trip(400, 1024, 2, 2<<16);
        check_round_trip(800, 1024, 3, 2<<16);

        // Test different numbers of shards
        check_round_trip(1000, 1 << 4, 2, 2<<16);
        check_round_trip(2000, 1 << 8, 2, 2<<16);
        check_round_trip(20, 1 << 12, 2, 2<<16);
        check_round_trip(100, 1 << 16, 2, 2<<16);
    }


    #[test]
    fn test_shard_round_trip_big() {

        // Test different numbers of threads
        check_round_trip(1024, 8192, 4, 2<<16);
        check_round_trip(2048, 1024, 1, 2<<16);
        check_round_trip(4096, 1024, 2, 2<<16);
        check_round_trip(8192, 1024, 3, 2<<16);

        // Test different numbers of shards
        check_round_trip(8192, 1 << 4, 2, 2<<16);
        check_round_trip(4096, 1 << 8, 2, 2<<16);
        check_round_trip(2048, 1 << 12, 2, 2<<16);
        check_round_trip(1024, 1 << 16, 2, 2<<16);
    }


    fn check_round_trip(shard_buf_size: usize, n_shards: usize, thread_bits: usize, n_items: usize) {
        
        println!("test round trip: n_shards: {}, thread_bits: {}, n_items: {}", n_shards, thread_bits, n_items);

        let ser = T1Ser{};

        let tmp = tempfile::NamedTempFile::new().unwrap();

        // Write and close file
        let true_items = {
            let manager = ShardWriteManager::new(tmp.path(), shard_buf_size, n_shards, thread_bits, ser.clone());
            let mut true_items = Vec::new();

            // Sender must be closed
            {
                let mut sender = manager.get_sender();

                for i in 0..n_items {
                    let tt = T1 {
                        a: (i/2) as u64,
                        b: i as u32,
                        c: (i * 2) as u16,
                        d: i as u8,
                    };
                    sender.send(tt);
                    true_items.push(tt);
                }
            }
            true_items
        };

        // Open finished file
        let reader = ShardReader::open(tmp.path(), ser.clone());

        let mut all_items = Vec::new();

        for i in 0..reader.num_shards() {
            let items = reader.read_shard(i);
            assert_eq!(reader.shard_len(i), items.len());
            all_items.extend(items);
        }

        all_items.sort_by_key(|x| x.a);

        if !(true_items == all_items) {
            println!("true len: {:?}", true_items.len());
            println!("round trip len: {:?}", all_items.len());
            assert!(false);
        }
    }
}
