#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# Put BAM file into buckets by read qname
#
import os.path
import tenkit.bam as tk_bam
import tenkit.bio_io as tk_io
import martian

__MRO__ = """
stage BUCKET_BY_QNAME(
    in  string[] qnames,
    in  bam      input,
    out map      buckets,
    out bam,
    src py       "stages/reads/bucket_reads_by_qname",
) split using (
    in  string   chunk_start,
    in  string   chunk_end,
)
"""

def split(args):
    bam_in = tk_bam.create_bam_infile(args.input)
    chunk_defs = tk_bam.chunk_bam_records(bam_in, chunk_bound_key=None, chunk_size_gb=0.5)
    return {'chunks': chunk_defs}

def bc_sort_key(read):
    return (tk_io.get_read_barcode(read), read.qname)

def main(args, outs):
    chunk_start = args.chunk_start
    chunk_end = args.chunk_end

    bam_in = tk_bam.create_bam_infile(args.input)
    reads = list(tk_bam.read_bam_chunk(bam_in, (chunk_start, chunk_end)))

    tmp_dir = os.path.dirname(outs.default)
    bams_out = {}
    outs.buckets = {}
    buckets = {}
    for qname in args.qnames:
        filename = os.path.join(tmp_dir, "qname_%s.bam" % qname)
        bam_out, _ = tk_bam.create_bam_outfile(filename, None, None, template=bam_in)
        bams_out[qname] = bam_out
        outs.buckets[qname] = filename
        buckets[qname] = []

    qname_ranges = zip(args.qnames, args.qnames[1:])
    for r in reads:
        qname = None
        for qnames in qname_ranges:
            if qnames[0] <= r.qname and r.qname <= qnames[1]:
                qname = qnames[0]
                break
        if qname is None:
            qname = args.qnames[-1]
        buckets[qname].append(r)

    for qname, bucket in buckets.iteritems():
        bucket.sort(key=bc_sort_key)
        bam_out = bams_out[qname]
        for r in bucket:
            bam_out.write(r)
        bam_out.close()

def join(args, outs, chunk_defs, chunk_outs):
    outs.buckets = {}
    for out in chunk_outs:
        for qname, filename in out.buckets.iteritems():
            if qname not in outs.buckets:
                outs.buckets[qname] = []
            outs.buckets[qname].append(filename)
