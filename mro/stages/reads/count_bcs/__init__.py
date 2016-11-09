#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# Put the sample index barcode and 10X barcode into read tags in a BAM file.
# Screen each 10X barcode against the given barcode whitelist
#
import gzip
import numpy as np
import tenkit.safe_json
import json
import tenkit.seq as tk_seq
import tenkit.fasta as tk_fasta

import martian

__MRO__ = """
stage COUNT_BCS(
    in  string barcode_whitelist,
    in  map[]  chunks,
    out json   bc_counts,
    src py     "stages/reads/count_bcs",
) split using (
    in  map chunk,
)
"""

def split(args):
    ''' Attach BCS to each chunk of the input files '''

    chunk_defs = [{'chunk':c} for c in args.chunks]
    return {'chunks': chunk_defs}


def join(args, outs, chunk_defs, chunk_outs):

    valid_counts = [c.bc_counts for c in chunk_outs if c.bc_counts is not None]

    # No counts if there's no whitelist or actual couts
    if args.barcode_whitelist is None or len(valid_counts) == 0:
        outs.bc_counts = None
        return

    result = {}

    for (c_out, c_def) in zip(chunk_outs, chunk_defs):
        gem_group = c_def.chunk['gem_group']
        if c_out.bc_counts is None:
            continue

        with open(c_out.bc_counts) as f:
            r = json.load(f)

        gg_result = result.setdefault(gem_group, {'bad_bc_count':0, 'bc_counts':None})

        gg_result['bad_bc_count'] += r['bad_bc_count']

        if gg_result['bc_counts'] is None:
            gg_result['bc_counts'] = np.array(r['bc_counts'], dtype=np.int32)
        else:
            gg_result['bc_counts'] += np.array(r['bc_counts'], dtype=np.int32)


    for gg in result.keys():
        rgg = result[gg]
        rgg['bc_error_rate'] = float(rgg['bad_bc_count']) / (float(rgg['bad_bc_count'] + rgg['bc_counts'].sum()))

    with open(outs.bc_counts, 'w') as f:
        tenkit.safe_json.dump_numpy(result, f)

def main(args, outs):
    """ Attaches barcodes. Attaches raw barcode to RAW_BC tag and filters those to form set of PROCESSES_BARCODES """

    # Bail out if there's no barcodes or whitelist
    if args.barcode_whitelist is None or args.chunk['barcode'] is None:
        outs.bc_counts = None
        return

    def open_maybe_gzip(fn):
        if fn[-2:] == "gz":
            return gzip.open(fn)
        else:
            return open(fn)

    barcode_whitelist = sorted(list(tk_seq.load_barcode_whitelist(args.barcode_whitelist)))
    bc_idx = { bc:idx for (idx,bc) in enumerate(barcode_whitelist) }
    bc_counts = np.zeros(len(barcode_whitelist), dtype=np.int32)
    bad_count = 0

    barcode_file = open_maybe_gzip(args.chunk['barcode'])
    bc_iterator = tk_fasta.read_generator_fastq(barcode_file)

    for (bc_read, raw_bc_seq, raw_bc_qual) in bc_iterator:
        idx = bc_idx.get(raw_bc_seq)

        if idx is not None:
            bc_counts[idx] += 1
        else:
            bad_count += 1

    # Write BC count array and bad count to pickle
    result = {}
    result['bad_bc_count'] = bad_count
    result['bc_counts'] = list(bc_counts)

    with open(outs.bc_counts, 'w') as bc_counts_out:
        tenkit.safe_json.dump_numpy(result, bc_counts_out)
