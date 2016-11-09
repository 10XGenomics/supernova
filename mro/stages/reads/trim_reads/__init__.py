#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# Consume multiple fastq files and spit out a fastq file with specified number of leading/trailing bases removed from each read
#
import os.path
import tenkit.fasta as tk_fasta
import tenkit.seq as tk_seq
import tenkit.stats as tk_stats
import itertools
import martian
import subprocess
import random
import numpy as np
import tenkit.safe_json
import json

from tenkit.constants import WHITELIST_TO_LOT_MAP

__MRO__ = """
stage TRIM_READS(
    in  map[]    chunks,
    int string   barcode_whitelist,
    in  int      max_read_num,
    in  int      read1_trim_length            "this is the trim length",
    in  int      read2_trim_length,
    out map[]    chunks,
    out fastq    placeholder,
    out json     bc_counts,
    out json     lot_info,
    src py       "stages/reads/trim_reads",
) split using (
    in map chunk,
)
"""

def split(args):
    mem_gb = 3
    chunks = [{'chunk': c, "__mem_gb": mem_gb} for c in args.chunks]
    return {'chunks': chunks}

def join(args, outs, chunk_defs, chunk_outs):
    final_chunks = []

    for cl in chunk_outs:
        final_chunks.extend(cl.chunks)

    outs.chunks = final_chunks
    valid_counts = [c.bc_counts for c in chunk_outs if c.bc_counts is not None]

    # No counts if there's no whitelist or actual counts
    if args.barcode_whitelist is None or len(valid_counts) == 0:
        outs.bc_counts = None
        outs.lot_info = None
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

    total_counts = 0
    total_errors = 0
    for gg in result.keys():
        rgg = result[gg]
        rgg['bc_error_rate'] = tk_stats.robust_divide(float(rgg['bad_bc_count']), float(rgg['bad_bc_count'] + rgg['bc_counts'].sum()))
        total_counts += float(rgg['bad_bc_count'] + rgg['bc_counts'].sum())
        total_errors += float(rgg['bad_bc_count'])

    # Hardcoded bail-out if the BC-correct rate is extremely high
    bc_error_rate = total_errors / total_counts
    if bc_error_rate > 0.97:
        martian.exit("Extremely high rate of incorrect barcodes observed (%.2f %%). Check that input is 10x Chromium data, and that there are no missing cycles in the first 16bp of Read 1. Please note that Supernova does not support GemCode data." % (bc_error_rate * 100.0))


    # possibly do lot detection
    lot_detection = {}
    lot_map = WHITELIST_TO_LOT_MAP.get(args.barcode_whitelist)
    if lot_map is not None:
        # get BC counts histogram
        # for now, just sum over all gem groups
        bc_seq = sorted(list(tk_seq.load_barcode_whitelist(args.barcode_whitelist)))
        bc_cts = np.sum([ggr['bc_counts'] for ggr in result.values()], axis=0)
        bc_hist = {seq: cts for seq, cts in zip(bc_seq, bc_cts)}

        (gelbead_lot, gelbead_lot_confidence, gelbead_lot_counts) = identify_gelbead_lot(bc_hist, lot_map)
        # only report on lots with nonzero counts
        gelbead_lot_counts_nonzero = {lot: count for lot, count in gelbead_lot_counts.items() if count > 0}

        lot_detection['gelbead_lot'] = gelbead_lot
        lot_detection['gelbead_lot_confidence'] = gelbead_lot_confidence
        lot_detection['gelbead_lot_counts'] = gelbead_lot_counts_nonzero

        martian.log_info("Gelbead lot detected: %s, reason (if None): %s" % (gelbead_lot, gelbead_lot_confidence))

    with open(outs.lot_info, 'w') as f:
        tenkit.safe_json.dump_numpy(lot_detection, f, pretty=True)

    with open(outs.bc_counts, 'w') as f:
        tenkit.safe_json.dump_numpy(result, f)

def openfq(fn):
    if fn[-2:] == 'gz':
        proc = subprocess.Popen(["gunzip", "--stdout", fn], stdout=subprocess.PIPE)
        return proc.stdout
    else:
        return open(fn, 'r')

def identify_gelbead_lot(bc_hist, lot_to_bcs, min_frac=0.95, min_counts=1000):
    # flip from {lot: [bcs]} to {bc: lot} for easier lookup
    bc_to_lot = {bc: lot for lot, bc_list in lot_to_bcs.items() for bc in bc_list}
    lot_counts = {lot: 0 for lot in lot_to_bcs}
    part_a_len = len(bc_to_lot.keys()[0]) # infer length of part A

    # count the total BC observations for each lot
    for bc, count in bc_hist.items():
        part_a = bc[:part_a_len]
        lot = bc_to_lot.get(part_a)
        if lot is not None:
            lot_counts[lot] += count

    # determine best match
    best_lot = max(lot_counts, key=lambda lot: lot_counts[lot])
    best_counts = lot_counts[best_lot]
    total_counts = sum(lot_counts.values())

    best_frac = 1.0 * best_counts / total_counts if total_counts > 0 else 0.0

    if best_frac >= min_frac and total_counts >= min_counts:
        result_lot = best_lot
        result_conf = "confident"
    elif total_counts < min_counts:
        result_lot = None
        result_conf = "insufficient data"
    else:
        result_lot = None
        result_conf = "ambiguous"

    return (result_lot, result_conf, lot_counts)

def main(args, outs):
    """ Trim the reads in a series of fasta files """

    # Set a fixed random seed to eliminate noise in metrics
    random.seed(0)

    chunk = args.chunk
    interleaved = chunk['reads_interleaved']
    have_read2 = chunk['read2'] is not None
    paired = interleaved or have_read2

    read1_trim = args.read1_trim_length
    read2_trim = args.read2_trim_length

    subsample_rate=chunk['subsample_rate']

    # BC config -- BC come from separate fastq, or are embedded in R1 or R2
    have_barcode = False
    bc_in_read1 = False
    bc_in_read2 = False
    bc_in_fastq = False

    # If we have bc in read, use that & ignore a separate BC read
    if chunk.get('bc_in_read', None) is not None and chunk.get('bc_length', 0) > 0:
        have_barcode = True
        bc_length = chunk['bc_length']
        if chunk['bc_in_read'] == 1:
            bc_in_read1 = True
            read1_trim += bc_length
        elif chunk['bc_in_read']== 2:
            bc_in_read2 = True
            read2_trim += bc_length
        else:
            martian.exit("bc_in_read configuration incorrect -- read must be 1 or 2")

    # Otherwise use the BC file
    elif chunk['barcode'] is not None:
        have_barcode = True
        bc_in_fastq = True

    have_sample_index = chunk['sample_index'] is not None

    output_directory = os.path.dirname(os.path.realpath(outs.placeholder))
    max_read_num = args.max_read_num

    # counter for sub-chunked files
    file_number = 1

    # open the available read files and make the appropriate iterators
    if interleaved:
        read_in = openfq(chunk['read1'])
        read_iter = tk_fasta.read_generator_fastq(read_in, paired_end=True)
    else:
        if have_read2:
            read1_in = openfq(chunk['read1'])
            read1_iter = tk_fasta.read_generator_fastq(read1_in)

            read2_in = openfq(chunk['read2'])
            read2_iter = tk_fasta.read_generator_fastq(read2_in)

            read_iter = itertools.imap(lambda x,y: (x[0],x[1],x[2],y[0],y[1],y[2]), read1_iter, read2_iter)
        else:
            read1_in = openfq(chunk['read1'])
            read_iter = tk_fasta.read_generator_fastq(read1_in)


    # open read file
    read_name = output_directory + "/read"+str(file_number)+".fastq"
    read_names = [read_name]
    out_read_fastq = open(read_name, 'w')

    # Bail out if there's no barcodes or whitelist
    if args.barcode_whitelist is None:
        outs.bc_counts = None
        bc_idx = None
    else:
        barcode_whitelist = sorted(list(tk_seq.load_barcode_whitelist(args.barcode_whitelist)))
        bc_idx = { bc:idx for (idx,bc) in enumerate(barcode_whitelist) }
        bc_counts = np.zeros(len(barcode_whitelist), dtype=np.int32)
        bad_count = 0

    # open barcode file if there is one
    if have_barcode:
        bc_name = output_directory + "/BC"+str(file_number)+".fastq"
        out_bc_fastq = open(bc_name,'w')
        bc_names = [bc_name]
        if bc_in_fastq:
            bc_in = openfq(chunk['barcode'])
            bc_iter = tk_fasta.read_generator_fastq(bc_in)
        elif bc_in_read1 or bc_in_read2:
            # BC in read -- have output file but no input file
            bc_iter = itertools.repeat(None)
    else:
        bc_iter = itertools.repeat(None)
        bc_names = [None]
        outs.bc_counts = None


   # open sample_index file if there is one
    if have_sample_index:
        si_name = output_directory + "/SI"+str(file_number)+".fastq"
        out_si_fastq = open(si_name,'w')
        si_in = openfq(chunk['sample_index'])
        si_iter = tk_fasta.read_generator_fastq(si_in)
        si_names = [si_name]
    else:
        si_iter = itertools.repeat(None)
        si_names = [None]

    # loop through reads
    read_num = 0
    for read, barcode_read, sample_index_read in itertools.izip(read_iter, bc_iter, si_iter):
        if read_num > 0 and random.random() > subsample_rate:
            continue

        if paired:
            (name1, seq1, qual1, name2, seq2, qual2) = read
        else:
            (name1, seq1, qual1) = read

        new_seq1 = seq1[read1_trim:]
        new_qual1 = qual1[read1_trim:]
        if paired:
            new_seq2 = seq2[read2_trim:]
            new_qual2 = qual2[read2_trim:]

        # Get BC sequence out of the read, for BC-in-read schemes
        if bc_in_read1:
            barcode_read = (name1, seq1[:bc_length], qual1[:bc_length])

        if bc_in_read2:
            barcode_read = (name2, seq2[:bc_length], qual2[:bc_length])

        read_num += 1
        if read_num > max_read_num:
            read_num = 1
            file_number += 1
            read_name = output_directory + "/read"+str(file_number)+".fastq"
            out_read_fastq.close()
            out_read_fastq = open(read_name,'w')
            read_names.append(read_name)

            if have_barcode:
                bc_name = output_directory + "/BC"+str(file_number)+".fastq"
                out_bc_fastq.close()
                out_bc_fastq = open(bc_name,'w')
                bc_names.append(bc_name)
            else:
                bc_names.append(None)

            if have_sample_index:
                si_name = output_directory + "/SI"+str(file_number)+".fastq"
                out_si_fastq.close()
                out_si_fastq = open(si_name,'w')
                si_names.append(si_name)
            else:
                si_names.append(None)

        if have_barcode:
            barcode_seq = barcode_read[1]
            barcode_qual = barcode_read[2]
            if chunk['barcode_reverse_complement']:
                barcode_seq = tk_seq.get_rev_comp(barcode_seq)
                barcode_qual = barcode_qual[::-1] # obscure way to reverse string
            if bc_idx is not None:
                idx = bc_idx.get(barcode_seq)
                if idx is not None:
                    bc_counts[idx] += 1
                else:
                    bad_count += 1

            tk_fasta.write_read_fastq(out_bc_fastq, barcode_read[0], barcode_seq, barcode_qual)
        if have_sample_index:
            tk_fasta.write_read_fastq(out_si_fastq, sample_index_read[0], sample_index_read[1], sample_index_read[2])

        tk_fasta.write_read_fastq(out_read_fastq, name1, new_seq1, new_qual1)
        if paired:
            tk_fasta.write_read_fastq(out_read_fastq, name2, new_seq2, new_qual2)

    if have_barcode:
        out_bc_fastq.close()
        # Only emit BC counts if we had a whitelist
        if outs.bc_counts is not None:
            result = {}
            result['bad_bc_count'] = bad_count
            result['bc_counts'] = list(bc_counts)
            with open(outs.bc_counts, 'w') as bc_counts_out:
                tenkit.safe_json.dump_numpy(result, bc_counts_out)
    if have_sample_index:
        out_si_fastq.close()
    out_read_fastq.close()

    chunks = []
    for (r,bc,si) in zip(read_names, bc_names, si_names):
        new_chunk = {
            'read1': r,
            'read2': None,
            'barcode': bc,
            'sample_index': si,
            'barcode_reverse_complement': False,
            'reads_interleaved': have_read2 or interleaved,
            'gem_group': chunk['gem_group'],
            'read_group': chunk['read_group']
        }
        chunks.append(new_chunk)

    outs.chunks = chunks
