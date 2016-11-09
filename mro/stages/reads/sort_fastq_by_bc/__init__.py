#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# Sort BAM file by sorting buckets, then concatenating bucket files
#

import subprocess
import json
import os
__MRO__ = """
stage SORT_FASTQ_BY_BC(
    in  map    bc_buckets,
    out fastq.gz[] reads,
    out int[] gem_group,
    out bool[] valid_bc,
    out string[] prefix,
    src py     "stages/reads/sort_reads_by_bc",
) split using (
    in  int    index,
    in  string prefix,
    in  string gem_group,
    in  map[]  bucket,
)
"""
def split(args):
    with open(args.bc_buckets) as buckets:
        bc_buckets = json.load(buckets)
    chunk_defs = []
    for gem_group, gem_group_bucket in bc_buckets.iteritems():
        for prefix, prefix_bucket in gem_group_bucket.iteritems():
            chunk_defs += [{'gem_group': gem_group, 'prefix': prefix, 'bucket': prefix_bucket, '__mem_gb': 3.0}]
    return {'chunks': chunk_defs}

def main(args, outs):
    output_dir = os.path.dirname(os.path.realpath(outs.reads))
    files = {}
    for index, file_map in enumerate(args.bucket):
        files[str(index)] = file_map
    barcode_bucket_json = outs.reads+"file_map.json"
    with open(barcode_bucket_json,'w') as json_file:
        json.dump(files, json_file)

    subprocess.check_call(['sort_fastq_by_bc',
                           '-prefix='+args.prefix,
                           '-barcode_bucket_json='+barcode_bucket_json,
                           '-output_dir='+output_dir])
    outs.reads = output_dir + "/reads.fastq.gz"
    outs.gem_group = args.gem_group
    outs.valid_bc = not "N" in args.prefix
    outs.prefix = args.prefix

def join(args, outs, chunk_defs, chunk_outs):
    outs.reads = [chunk_out.reads for chunk_out in chunk_outs]
    outs.gem_group = [int(chunk_out.gem_group) for chunk_out in chunk_outs]
    outs.prefix = [chunk_out.prefix for chunk_out in chunk_outs]
    outs.valid_bc = [chunk_out.valid_bc for chunk_out in chunk_outs]
