#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# Group FASTQs by 10X barcode.
#
__MRO__ = """
stage BUCKET_FASTQ_BY_BC(
    in  map[]    chunk,
    in  int      buckets,
    in  string   barcode_whitelist,
    in  float    max_expected_barcode_errors,
    in  float    bc_confidence_threshold,
    in  json     barcode_counts,
    out json     file_map,
    out int      gem_group,
    out fastq.gz,
    src py       "stages/reads/bucket_fastq_by_bc",
) split using (
    in  fastq.gz read_file,
    in  fastq.gz barcode_file,
    in  fastq.gz sample_index_file,
    in  int      gem_group,
)
"""

import subprocess
import os
import json
from tenkit.constants import BARCODE_LOCATION

import martian

def split(args):
    chunk = args.chunk
    chunks = [{'chunk':x, '__mem_gb':3.0} for x in chunk]
    return {'chunks': chunks}

def join(args, outs, chunk_defs, chunk_outs):
    buckets = {}
    directories = [os.path.dirname(os.path.realpath(chunk_out.default)) for chunk_out in chunk_outs]
    for directory in directories:
        for fastq in os.listdir(directory):
            chunk = int(fastq.split("_")[1])
            gem_group = int(fastq.split("_")[2][:-9]) #strip off .fastq.gz
            gem_group_file_map = buckets.setdefault(gem_group,{})
            chunk_list = gem_group_file_map.setdefault(chunk,[])
            chunk_list.append({"R":directory+"/"+fastq})
    with open(outs.file_map, 'w') as out_file:
        json.dump(buckets, out_file)

def main(args, outs):
    chunk = args.chunk
    assert(chunk['reads_interleaved'])
    if not chunk['reads_interleaved'] and (chunk['read1'] is None or chunk['read2'] is None):
        martian.throw("must supply a read1 and read2 when reads_interleave == False")
    output_dir = os.path.dirname(os.path.realpath(outs.default))
    if args.barcode_whitelist:
        barcode_whitelist = BARCODE_LOCATION + "/" + args.barcode_whitelist + ".txt"
    else:
        barcode_whitelist = "none"
    gem_group = chunk["gem_group"] or 1
    barcode_read = chunk["barcode"] or "none"
    barcode_counts = args.barcode_counts or "none"
    sample_index = chunk["sample_index"] or "none"
    read_group_string = chunk["read_group"] or "none"

    subprocess.check_call(['bucket_fastq_by_bc',
                           '-reads='+chunk["read1"],
                           '-read_group_string='+read_group_string,
                           '-barcodes='+barcode_read,
                           '-barcodeCounts='+barcode_counts,
                           '-bcConfidenceThreshold='+str(args.bc_confidence_threshold),
                           '-output_directory='+output_dir,
                           '-sample_index_reads='+sample_index,
                           '-gem_group='+str(gem_group),
                           '-barcode_whitelist='+barcode_whitelist,
                           '-interleaved='+str(chunk['reads_interleaved']),
                           '-max_expected_barcode_errors='+str(args.max_expected_barcode_errors),
                           '-buckets='+str(args.buckets)])
