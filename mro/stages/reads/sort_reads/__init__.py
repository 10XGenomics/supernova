#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# Sort BAM file
#
import os.path
import resource
import shutil
import tenkit.bam as tk_bam
import tenkit.read_filter

__MRO__ = """
stage SORT_BY_POS(
    in  bam input,
    out bam,
    src py  "stages/reads/sort_reads",
) split using (
    in  bam chunk_input,
)
"""

def split(args):
    chunk_defs = [{'chunk_input':x} for x in args.input if os.path.isfile(x)]
    return {'chunks': chunk_defs, 'join': {'__threads': 6}}

def main(args, outs):
    main_sort_reads(args, outs)


def main_sort_reads(args, outs):
    args.coerce_strings()
    bam_prefix, ext = os.path.splitext(outs.default)
    perfect_read_count = 0
    bam = tk_bam.create_bam_infile(str(args.chunk_input))
    tk_bam.sort(str(args.chunk_input), str(bam_prefix))
    while True:
        try:
            read = bam.next()
            if tenkit.read_filter.stringent_read_filter(read, True):
                perfect_read_count += 1
        except StopIteration:
            break
    outs.perfect_read_count = perfect_read_count


def merge(input_bams, output_bam, threads=1):
    ''' Merge the sorted bam chunks hierarchically to conserve open file handles '''
    soft, _ = resource.getrlimit(resource.RLIMIT_NOFILE)
    soft -= 100

    tmp_dir = os.path.dirname(output_bam)
    while len(input_bams) > 1:
        new_bams = []
        for i in range(0, len(input_bams), soft):
            bam_chunk = input_bams[i:i+soft]
            if len(bam_chunk) > 1:
                new_bam = os.path.join(tmp_dir, "%d-%d.bam" % (i, len(input_bams)))
                tk_bam.merge(new_bam, bam_chunk, threads)
                new_bams.append(new_bam)
            else:
                new_bams.append(input_bams[i])
        input_bams = new_bams
    shutil.move(input_bams[0], output_bam)

def join(args, outs, chunk_defs, chunk_outs):
    outs.coerce_strings()
    input_bams = [str(chunk.default) for chunk in chunk_outs]
    merge(input_bams, outs.default, args.__threads)
    tk_bam.index(outs.default)
    outs.perfect_read_count = sum([chunk.perfect_read_count for chunk in chunk_outs])
