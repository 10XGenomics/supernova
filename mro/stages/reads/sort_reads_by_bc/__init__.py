#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# Sort BAM file by sorting buckets, then concatenating bucket files
#

import pysam
import tenkit.bam as tk_bam
import tenkit.bio_io as tk_io
import tenkit.cache as tk_cache
import heapq
import martian

__MRO__ = """
stage SORT_BY_BC(
    in  map    bc_buckets,
    in  map    non_bc_buckets,
    in  bam    possorted_bam,
    out int    total_reads,
    out bam,
    src py     "stages/reads/sort_reads_by_bc",
) split using (
    in  int    index,
    in  string prefix,
    in  bam[]  bucket,
)
"""

def split(args):
    chunk_defs = []
    for i, buckets in enumerate([args.non_bc_buckets, args.bc_buckets]):
        chunk_defs += [{'prefix': prefix, 'bucket': bucket, 'index': i} for prefix, bucket in buckets.iteritems()]
    return {'chunks': chunk_defs}

def bc_sort_key(read):
    return (tk_io.get_read_barcode(read), read.qname)

def main(args, outs):
    outs.coerce_strings()
    bam_in = tk_bam.create_bam_infile(args.bucket[0])
    bam_out, _ = tk_bam.create_bam_outfile(outs.default, None, None, template=bam_in, pgs=tk_bam.make_pg_header(martian.get_pipelines_version(), "sort_reads_by_bc"))
    bam_in.close()

    outs.total_reads = merge_by_key(args.bucket, bc_sort_key, bam_out)
    bam_out.close()

def merge_by_key(bam_filenames, key_func, bam_out):
    file_cache = tk_cache.FileHandleCache(mode='rb', open_func=pysam.Samfile)
    total_reads = 0
    heap  = []

    for bam_filename in bam_filenames:
        try:
            bam = file_cache.get(bam_filename)
            first_read = bam.next()
            heapq.heappush(heap, (key_func(first_read), first_read, bam_filename))
        except StopIteration:
            pass

    while len(heap) > 0:
        # Get the minimum item and write it to the bam.
        key, read, bam_filename = heapq.heappop(heap)
        bam = file_cache.get(bam_filename)
        bam_out.write(read)
        total_reads += 1

        # Get the next read from the source bam we just wrote from
        # If that BAM file is out of reads, then we leave that one out
        try:
            next_read = bam.next()
            heapq.heappush(heap, (key_func(next_read), next_read, bam_filename))
        except StopIteration:
            pass

    return total_reads

def join(args, outs, chunk_defs, chunk_outs):
    chunk_lists = [[], []]
    outs.total_reads = 0
    for chunk in zip(chunk_defs, chunk_outs):
        index = chunk[0].index
        chunk_lists[index].append(chunk)
        outs.total_reads += chunk[1].total_reads

    # Sanity check vs. position-sorted BAM
    with tk_bam.create_bam_infile(args.possorted_bam) as possorted_bam_in:
        assert possorted_bam_in.unmapped + possorted_bam_in.mapped == outs.total_reads

    buckets = []
    for chunks in chunk_lists:
        chunks = sorted(chunks, key=lambda chunk: chunk[0].prefix)
        buckets += [chunk[1].default for chunk in chunks]
    tk_bam.concatenate(outs.default, buckets)
