#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# Setup read chunks.
#
import os.path
import tenkit.fasta as tk_fasta
import itertools
import tenkit.seq as tk_seq
import tenkit.stats as tk_stats
import tenkit.preflight as tk_preflight
import tenkit.safe_json
import tenkit.bam as tk_bam
from tenkit.constants import BARCODE_LOCATION
import martian
import gzip
import re

__MRO__ = """
stage SETUP_CHUNKS(
    in  string    sample_id          "id of the sample",
    in  map[]     sample_def         "list of dictionary specifying input data",
    in  string    input_mode         "configuration of the input fastqs",
    in  string    barcode_whitelist,
    in  map       downsample         "map specifies either subsample_rate (float) or gigabases (int)",
    out txt       barcode_whitelist_path,
    out map[]     chunks             "map has barcode, barcode_reverse_complement, sample_index, read1, read2, gem_group, and read_group fields",
    out string[]  read_groups        "list of strings representing read groups",
    out json      downsample_info    "info about downsampling"
    out txt       barcode_whitelist_path,
    out int       requested_read_pairs,
    src py        "stages/reads/setup_chunks",
)
"""

def main(args, outs):
    """Combine reads from multiple input FASTQ files, and potentially trim.
       Demultiplex outputs a series of FASTQ files with filenames of the form:
       read-[RA|I1|I2]_si-AGTAACGT_lane-001_chunk_001.fastq[.gz].
    """

    def check_key(n, dict_in, name, tys):
        if not dict_in.has_key(name):
            martian.exit("Entry %d in sample_def missing required field: %s" % (n, name))

        if not (type(dict_in[name]) in tys):
            martian.exit("Entry %d in sample_def for '%s' has incorrect type -- expecting %s, got %s" % (n, name, str(tys), type(dict_in[name])))


    # Check for self-consistent gem_group settings in the sample_def entries
    gem_groups = [x['gem_group'] for x in args.sample_def]
    all_null = all([x is None for x in gem_groups])
    all_int = all([type(x) is int for x in gem_groups])

    if not (all_null or all_int):
        martian.exit("Inconsistent gem_group tags. Please specify all gem_group tags as null, or all gem_group tags with an integer")

    # If all gem_groups are set to null, then set them all to 1
    if all_null:
        for sample_item in args.sample_def:
            sample_item['gem_group'] = 1

    # Predicted input bases
    total_seq_bases = 0

    # verify input mode upfront
    if args.input_mode not in ["BCL_PROCESSOR", "ILMN_BCL2FASTQ"]:
        martian.throw("Unrecognized input_mode: %s" % args.input_mode)

    for (idx, sample_item) in enumerate(args.sample_def):
        # validate fields
        check_key(idx, sample_item, "read_path", [str, unicode])
        check_key(idx, sample_item, "lanes",  [list, type(None)])
        check_key(idx, sample_item, "gem_group", [int, type(None)])
        if args.input_mode == "BCL_PROCESSOR":
            check_key(idx, sample_item, "sample_indices", [list, type(None)])
        elif args.input_mode == "ILMN_BCL2FASTQ":
            check_key(idx, sample_item, "sample_names", [list, type(None)])

    interleaved_read_type = "RA"

    chunks = []
    read_groups = set()

    for read_chunk in args.sample_def:
        # Each sample_def entry can have a separate pre-applied downsampling rate
        # We adjust the estimated data in that chunk to account for this
        # subsampling
        chunk_subsample_rate = read_chunk.get('subsample_rate', 1.0)

        bc_in_read = {}
        if read_chunk.has_key('bc_in_read'):
            if read_chunk['bc_in_read'] is not None:
                bc_in_read['bc_in_read'] = read_chunk['bc_in_read']
                bc_in_read['bc_length'] = read_chunk['bc_length']

        path = read_chunk['read_path']
        lanes = read_chunk['lanes']
        gem_group = read_chunk['gem_group']
        unbarcoded = read_chunk.get('unbarcoded')
        sample_id = args.sample_id
        library_id = read_chunk.get('library', 'MissingLibrary')

        # split on BCL_PROCESSOR / ILMN_BCL2FASTQ
        # the main difference is that BCL_PROCESSOR uses interleaved reads and labels FASTQs by sample index;
        # whereas ILMN_BCL2FASTQ uses R1/R2 and labels by sample name

        if args.input_mode == "BCL_PROCESSOR":
            sample_index_strings, msg = tk_preflight.check_sample_indices(read_chunk)
            if sample_index_strings is None:
                martian.exit(msg)

            sample_seq_bases = 0
            find_func = tk_fasta.find_input_fastq_files_10x_preprocess
            for sample_index in sample_index_strings:
                # process interleaved reads
                reads = find_func(path, interleaved_read_type, sample_index, lanes)
                for read in reads:
                    _, predicted_seq_bases, read_length = fastq_data_estimate(read)
                    sample_seq_bases += predicted_seq_bases

            sample_seq_bases = chunk_subsample_rate * sample_seq_bases
            bp_per_read_pair = 2*read_length

            martian.log_info("Input data: Predict %f GB from %s. (%d bp per read pair)" % (float(sample_seq_bases)/1e9, path, bp_per_read_pair))
            total_seq_bases += sample_seq_bases

            for sample_index in sample_index_strings:
                reads = find_func(path, interleaved_read_type, sample_index, lanes)
                # TODO confirm that this works with cellranger
                si_read, bc_read = ("I1", "I2")
                if 'barcode_read' in read_chunk and read_chunk['barcode_read'] == 'I1':
                    si_read, bc_read = ("I2", "I1")
                sis = find_func(path, si_read, sample_index, lanes)

                # allow empty sample index case if all reads in lane are same sample
                if sis is None or sis == []:
                    sis = [None] * len(reads)

                if not unbarcoded:
                    barcodes = find_func(path, bc_read, sample_index, lanes)
                    if len(barcodes) == 0:
                        barcodes = [None] * len(reads)
                else:
                    barcodes = [None] * len(reads)

                # calculate chunks
                for r,b,si in zip(reads, barcodes, sis):
                    (flowcell, lane) = get_run_data(r)
                    rg_string = tk_bam.pack_rg_string(sample_id, library_id, gem_group, flowcell, lane)
                    new_chunk = {
                        'read1': r, 'read2': None, 'reads_interleaved': True, 'barcode': b,
                        'sample_index': si, 'barcode_reverse_complement': False, 'gem_group': gem_group,
                        'subsample_rate': chunk_subsample_rate, 'read_group': rg_string
                    }
                    new_chunk.update(bc_in_read)
                    chunks.append(new_chunk)
                    read_groups.add(rg_string)

        elif args.input_mode == "ILMN_BCL2FASTQ":
            sample_names = read_chunk['sample_names']

            sample_seq_bases = 0
            find_func = tk_fasta.find_input_fastq_files_bcl2fastq_demult
            for sample_name in sample_names:
                # process read 1
                reads = find_func(path, "R1", sample_name, lanes)
                for read in reads:
                    _, predicted_seq_bases, read_length1 = fastq_data_estimate(read)
                    sample_seq_bases += predicted_seq_bases
                # process read 2
                reads = find_func(path, "R2", sample_name, lanes)
                for read in reads:
                    _, predicted_seq_bases, read_length2 = fastq_data_estimate(read)
                    sample_seq_bases += predicted_seq_bases

            sample_seq_bases = chunk_subsample_rate * sample_seq_bases
            bp_per_read_pair = read_length1 + read_length2

            martian.log_info("Input data: Predict %f GB from %s. (%d bp per read pair)" % (float(sample_seq_bases)/1e9, path, bp_per_read_pair))
            total_seq_bases += sample_seq_bases

            for sample_name in sample_names:
                r1_reads = find_func(path, "R1", sample_name, lanes)
                r2_reads = find_func(path, "R2", sample_name, lanes)

                # TODO confirm that this works with cellranger
                si_read, bc_read = ("I1", "I2")
                if 'barcode_read' in read_chunk and read_chunk['barcode_read'] == 'I1':
                    si_read, bc_read = ("I2", "I1")
                sis = find_func(path, si_read, sample_name, lanes)

                # allow empty sample index case if all reads in lane are same sample
                if sis is None or sis == []:
                    sis = [None] * len(r1_reads)

                # in Chromium chemistry... there shouldn't be separate barcode reads...
                if not unbarcoded:
                    barcodes = find_func(path, bc_read, sample_name, lanes)
                    if len(barcodes) == 0:
                        barcodes = [None] * len(r1_reads)
                else:
                    barcodes = [None] * len(r1_reads)

                # again, with Chromium, the barcodes should be an array of Nones, but
                # just in case...
                if not (len(r1_reads) == len(r2_reads) == len(barcodes)):
                    martian.log_info("Read 1 files: %s" % str(r1_reads))
                    martian.log_info("Read 2 files: %s" % str(r2_reads))
                    martian.log_info("Barcode files: %s" % str(barcodes))
                    martian.exit("Read1, Read2, and Barcode files are mismatched. Exiting pipline")

                # calculate chunks
                for r1,r2,b,si in zip(r1_reads, r2_reads, barcodes, sis):
                    (flowcell, lane) = get_run_data(r1)
                    rg_string = tk_bam.pack_rg_string(sample_id, library_id, gem_group, flowcell, lane)
                    new_chunk = {
                        'read1': r1, 'read2': r2, 'reads_interleaved': False, 'barcode': b,
                        'sample_index': si, 'barcode_reverse_complement': False, 'gem_group': gem_group,
                        'subsample_rate': chunk_subsample_rate, 'read_group': rg_string
                    }
                    new_chunk.update(bc_in_read)
                    chunks.append(new_chunk)
                    read_groups.add(rg_string)

    martian.log_info("Input data: Predict %f total GB" % (float(total_seq_bases)/1e9))

    if len(chunks) == 0:
        martian.exit("No input FASTQs were found for the requested parameters.")


    #
    # Downsampling setup
    #

    # The total available input raw gigabases of input data (est_gb), and the base pairs per read pair (bp_per_read_pair)
    # are estimated above.
    (est_gb, bp_per_read_pair) = (float(total_seq_bases)/1e9, bp_per_read_pair)

    downsample = args.downsample if args.downsample is not None else {}

    # Possible BC subsampling -- try to get the requested amount of data _after_ bc subsampling
    est_gb_post_bc = est_gb * downsample.get("bc_subsample_rate", 1.0)

    # Aim high to ensure that we won't be left with too few reads
    fudge_factor = 1.25

    downsample_succeeded = True

    if downsample.has_key("gigabases"):
        read_sample_rate = min(1.0, fudge_factor * downsample['gigabases'] / est_gb_post_bc)
        requested_read_pairs = int(1e9 * downsample['gigabases'] / bp_per_read_pair)
        downsample_succeeded = downsample['gigabases'] > est_gb_post_bc

    elif downsample.has_key("target_reads"):
        requested_read_pairs = int(downsample['target_reads'] / 2)
        est_read_pair_post_bc = 1e9 * est_gb_post_bc / bp_per_read_pair
        read_sample_rate = min(1.0, fudge_factor * requested_read_pairs / est_read_pair_post_bc)
        downsample_succeeded = requested_read_pairs > est_read_pair_post_bc

    elif downsample.has_key("subsample_rate"):
        read_sample_rate = min(1.0, downsample["subsample_rate"] / downsample.get("bc_subsample_rate", 1.0))
        requested_read_pairs = None

    else:
        read_sample_rate = 1.0
        requested_read_pairs = None


    martian.log_info("Downsampling request: %s" % str(downsample))
    martian.log_info("Base pairs per read pair: %s" % bp_per_read_pair)
    martian.log_info("Estimated Input: %.2f GB, Initial Downsample Rate: %.3f. Requested total reads: %s" % (est_gb, read_sample_rate, str(requested_read_pairs)))

    # Copy over the per-chunk subsample rates
    if read_sample_rate is not None:
        for chunk in chunks:
            chunk['subsample_rate'] = chunk.get('subsample_rate', 1.0) * read_sample_rate
            if downsample.has_key("bc_subsample_rate"):
                chunk["bc_subsample_rate"] = downsample["bc_subsample_rate"]

    outs.requested_read_pairs = requested_read_pairs

    martian.log_info("Input reads: %s" % str(chunks))
    outs.chunks = chunks
    outs.read_groups = [rg for rg in read_groups]

    downsample_info = {}
    downsample_info['available_gb'] = est_gb
    downsample_info['requested_gb'] = downsample.get('gigabases', None)
    downsample_info['requested_rate'] = read_sample_rate
    downsample_info['post_downsample_gb'] = float(requested_read_pairs * bp_per_read_pair) / 1e9 if requested_read_pairs is not None else None
    downsample_info['downsample_succeeded'] = downsample_succeeded

    with open(outs.downsample_info, 'w') as downsample_out:
        tenkit.safe_json.dump_numpy(downsample_info, downsample_out)

    check_fastqs(outs.chunks)

    # Give out full path to BC whitelist
    if args.barcode_whitelist:
        outs.barcode_whitelist_path = BARCODE_LOCATION + "/" + args.barcode_whitelist + ".txt"
    else:
        outs.barcode_whitelist_path = None


def check_fastqs(chunks):
    keys = ['read1', 'read2', 'barcode', 'sample_index']

    def check_fastq(fastq):
        # Check if fastq is readable
        if not os.access(fastq, os.R_OK):
            martian.exit("Do not have file read permission for FASTQ file: %s" % fastq)

        # Check if fastq is gzipped
        gzip_suffix = '.gz'
        is_gzip_fastq = True
        try:
            with gzip.open(fastq) as f:
                f.read(1)
        except:
            is_gzip_fastq = False

        if is_gzip_fastq and not fastq.endswith(gzip_suffix):
            martian.exit("Input FASTQ file is gzipped but filename does not have %s suffix: %s" % (fastq, gzip_suffix))
        if not is_gzip_fastq and fastq.endswith(gzip_suffix):
            martian.exit("Input FASTQ file is not gzipped but filename has %s suffix: %s" % (fastq, gzip_suffix))

    for chunk in chunks:
        for key in keys:
            fastq = chunk.get(key)
            if fastq is not None:
                check_fastq(fastq)

def infer_barcode_reverse_complement(barcode_whitelist, barcode_files):
    rc_valid_count = 0
    reg_valid_count = 0
    if barcode_whitelist:
        barcode_rc = []
        for barcode_file in barcode_files:
            read_num = 0

            if barcode_file[-3:] == ".gz":
                barcode_open_file = gzip.open(barcode_file)
            else:
                barcode_open_file = open(barcode_file, 'r')
            read_iter = tk_fasta.read_generator_fastq(barcode_open_file)
            for (name, seq, qual) in read_iter:
                if seq in barcode_whitelist:
                    reg_valid_count += 1
                if tk_seq.get_rev_comp(seq) in barcode_whitelist:
                    rc_valid_count += 1
                if read_num > 1000:
                    break
                read_num += 1

            if tk_stats.robust_divide(float(rc_valid_count), float(rc_valid_count + reg_valid_count)) > 0.75:
                barcode_rc.append(True)
            else:
                barcode_rc.append(False)
            barcode_open_file.close()
        return barcode_rc
    else:
        return [False] * len(barcode_files)


import struct

def get_bsize(f, subfield_start, extra_len):
    (id1, id2, length) = struct.unpack("BBH", f.read(4))

    if id1 == 66 and id2 == 67 and length == 2:
        bsize = struct.unpack("H", f.read(2))
        return bsize[0]
    #elif f.tell() - subfield_start < extra_len:
    #    get_bsize(f, subfield_start, extra_len)
    else:
        raise Exception("BSIZE field not found -- this a valid BAM file?")

def parse_bgzf_header(f):
    cur_pos = f.tell()
    header_fmt = "BBBBIBBH"

    d = f.read(12)

    # We are at EOF when read returns an empty string
    if d == '':
        return None

    header = struct.unpack(header_fmt, d)

    # Check for a valid gzip header
    if header[0] != 31 or header[1] != 139:
        raise Exception("Not a valid gzip header")

    xlen = header[7]
    bsize = get_bsize(f, f.tell(), xlen)

    next_pos = cur_pos + bsize + 1
    f.seek(next_pos)
    return next_pos

def parse_gzip_sz(fn):
    sz = os.path.getsize(fn)
    f = open(fn, 'r')
    f.seek(sz-4)

    d = f.read(4)
    sz = struct.unpack("I", d)[0]
    return sz

'''
# TODO
def guess_uncompressed_data(fn):
    is_bgzf = True
    try:
        f = open(fn, "r")
        n = parse_bgzf_header(f)
    except:
        is_bgzf = False

    if is_bgzf:
        # Estimate bgzf size
        pass
    else:
        pass
'''

def estimate_gzip_uncompressed_size(fn):
    ''' Estimate the uncompressed size of a non-BGZF gzip file, using the best estimate from the size mod 2^32 encoded in the gzip chunk '''
    file_sz = os.path.getsize(fn)
    reader = gzip.open(fn)

    # Read through 50k reads or ~4% of the file, whichever is greater, to estimate the compression ratio
    num_to_read = max(50000, 0.04 * (file_sz * 4 / 200))

    lens = [len(l) for l in itertools.islice(reader, num_to_read)]
    uncompressed_sz = sum(lens)
    compressed_sz = reader.myfileobj.tell()

    # Rough guess of the uncompressed size
    estimated_full_sz = 1.01 * float(uncompressed_sz) / compressed_sz * file_sz
    uncompressed_full_sz_mod = parse_gzip_sz(fn)

    min_err = 1.0e12
    true_size = 0
    max_units = max(int(float(estimated_full_sz) * 2 / (2**32)), 2)

    for true_size_opt in [ uncompressed_full_sz_mod + 2**32 * i for i in range(max_units + 1) ]:
        err = abs(estimated_full_sz - true_size_opt)
        if err < min_err:
            min_err = err
            true_size = true_size_opt

    assert(min_err <= 2**31)
    return (true_size, estimated_full_sz)


def fastq_data_estimate(fn, num_reads = 5000):
    # Open reader
    if fn[-2:] == 'gz':
        reader = gzip.open(fn)
        is_gz = True
    else:
        reader = open(fn, 'r')
        is_gz = False

    gen = tk_fasta.read_generator_fastq(reader)
    rds = itertools.islice(gen, num_reads)

    input_lens = [(len(header) + len(r) + len(qual) + 4, len(r)) for (header,r,qual) in rds]
    total_seq_len = sum(x[1] for x in input_lens)
    total_data_len = sum(x[0] for x in input_lens)
    file_sz = os.path.getsize(fn)

    read_length = total_seq_len / len(input_lens)

    if is_gz:
        (uncomp_size, predicted_sz) = estimate_gzip_uncompressed_size(fn)
    else:
        uncomp_size = file_sz
        predicted_sz = file_sz

    read_yield = float(len(input_lens)) / total_data_len
    seq_yield = float(total_seq_len) / total_data_len
    predicted_reads = read_yield * uncomp_size
    predicted_seq = seq_yield * uncomp_size

    # Log estimate of downsampling
    gzip_sz = parse_gzip_sz(fn)
    martian.log_info("Estimates for: %s" % fn)
    dbg_str =  "compressed_size: %.2f, predicted_size: %.2f, predicted_size_mod: %.2f, gzip_size_mod: %.2f, gzip_predicted_size: %.2f" % (float(file_sz)/1e9, float(predicted_sz)/1e9, float(predicted_sz % 2**32)/1e9, float(gzip_sz)/1e9, float(uncomp_size)/1e9)
    martian.log_info(dbg_str)

    return (predicted_reads, predicted_seq, read_length)

def get_run_data(fn):
    """ Parse flowcell + lane from the first FASTQ record.
    NOTE: we don't check whether there are multiple FC / lanes in this file.
    """
    if fn[-2:] == 'gz':
        reader = gzip.open(fn)
    else:
        reader = open(fn, 'r')

    gen = tk_fasta.read_generator_fastq(reader)

    try:
        (name, seq, qual) = gen.next()
        (flowcell, lane) = re.split(':', name)[2:4]
        return (flowcell, lane)
    except StopIteration:
        # empty fastq
        martian.exit("FASTQ is empty: %s" % fn)
