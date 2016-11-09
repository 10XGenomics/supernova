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
    out map[]     chunks             "map has barcode, barcode_reverse_complement, sample_index, read1, read2, gem_group, and read_group fields",
    out string[]  read_groups        "list of strings representing read groups",
    out json      downsample_info    "info about downsampling"
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

    global_subsample_rate = 1.0
    downsample_gigabases = False
    downsample_reads      = False
    if args.downsample is not None:
        ## make sure that exactly one downsampling option is specified
        options_supplied=0
        for subsample_key in ["gigabases", "subsample_rate", "target_reads"]:
            if args.downsample.get(subsample_key, None) is not None:
                options_supplied += 1
        assert( options_supplied == 1 )
        ##
        if 'subsample_rate' in args.downsample and args.downsample['subsample_rate'] is not None:
            global_subsample_rate = args.downsample['subsample_rate']
            assert( global_subsample_rate <= 1.0 )
        elif 'target_reads' in args.downsample and args.downsample['target_reads'] is not None:
            downsample_reads = True
        else:
            downsample_gigabases = True

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
    total_seq_reads = 0

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
        # Check if subsample_rate exists in sample_def
        if 'subsample_rate' in read_chunk.keys():
            subsample_rate = global_subsample_rate * read_chunk['subsample_rate']
        else:
            subsample_rate = global_subsample_rate

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
        library_id = read_chunk.get('library_id', 'MissingLibrary')

        # split on BCL_PROCESSOR / ILMN_BCL2FASTQ
        # the main difference is that BCL_PROCESSOR uses interleaved reads and labels FASTQs by sample index;
        # whereas ILMN_BCL2FASTQ uses R1/R2 and labels by sample name

        if args.input_mode == "BCL_PROCESSOR":
            sample_index_strings, msg = tk_preflight.check_sample_indices(read_chunk)
            if sample_index_strings is None:
                martian.exit(msg)

            sample_seq_bases = 0
            sample_seq_reads = 0
            find_func = tk_fasta.find_input_fastq_files_10x_preprocess
            for sample_index in sample_index_strings:
                # process interleaved reads
                reads = find_func(path, interleaved_read_type, sample_index, lanes)
                for read in reads:
                    predicted_seq_reads, predicted_seq_bases = fastq_data_estimate(read)
                    sample_seq_bases += predicted_seq_bases
                    sample_seq_reads += predicted_seq_reads

            martian.log_info("Input data: Predict %f GB from %s" % (float(sample_seq_bases)/1e9, path))
            total_seq_bases += sample_seq_bases
            total_seq_reads += sample_seq_reads

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
                    rg_string = ':'.join([sample_id, library_id, str(gem_group), flowcell, lane])
                    new_chunk = {
                        'read1': r, 'read2': None, 'reads_interleaved': True, 'barcode': b, 
                        'sample_index': si, 'barcode_reverse_complement': False, 'gem_group': gem_group,
                        'subsample_rate': subsample_rate, 'read_group': rg_string
                    }
                    new_chunk.update(bc_in_read)
                    chunks.append(new_chunk)
                    read_groups.add(rg_string)

        elif args.input_mode == "ILMN_BCL2FASTQ":
            sample_names = read_chunk['sample_names']

            sample_seq_bases = 0
            sample_seq_reads = 0
            find_func = tk_fasta.find_input_fastq_files_bcl2fastq_demult
            for sample_name in sample_names:
                # process read 1
                reads = find_func(path, "R1", sample_name, lanes)
                for read in reads:
                    predicted_seq_reads, predicted_seq_bases = fastq_data_estimate(read)
                    sample_seq_bases += predicted_seq_bases
                    sample_seq_reads += predicted_seq_reads
                # process read 2
                reads = find_func(path, "R2", sample_name, lanes)
                for read in reads:
                    predicted_seq_reads, predicted_seq_bases = fastq_data_estimate(read)
                    sample_seq_bases += predicted_seq_bases
                    sample_seq_reads += predicted_seq_reads

            martian.log_info("Input data: Predict %f GB from %s" % (float(sample_seq_bases)/1e9, path))
            total_seq_bases += sample_seq_bases
            total_seq_reads += sample_seq_reads

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
                    rg_string = ':'.join([sample_id, library_id, str(gem_group), flowcell, lane])
                    new_chunk = {
                        'read1': r1, 'read2': r2, 'reads_interleaved': False, 'barcode': b,
                        'sample_index': si, 'barcode_reverse_complement': False, 'gem_group': gem_group,
                        'subsample_rate': subsample_rate, 'read_group': rg_string
                    }
                    new_chunk.update(bc_in_read)
                    chunks.append(new_chunk)
                    read_groups.add(rg_string)

    martian.log_info("Input data: Predict %f total GB" % (float(total_seq_bases)/1e9))
    martian.log_info("            Predict %d total reads" % total_seq_reads)

    if len(chunks) == 0:
        martian.exit("No input FASTQs were found for the requested parameters.")

    if downsample_gigabases and args.downsample['gigabases'] is not None:
        # Calculate global downsample rate
        global_subsample_rate = min(1.0, float(args.downsample['gigabases'])*1e9 / float(total_seq_bases))
        martian.log_info("Input data downsampling: Requested: %.2f GB, Estimated Input: %.2f GB, Downsample Rate: %.3f" \
         % (float(args.downsample['gigabases']), float(total_seq_bases)/1e9, global_subsample_rate))

        for chunk in chunks:
            chunk['subsample_rate'] = chunk['subsample_rate'] * global_subsample_rate
    elif downsample_reads:
        global_subsample_rate = min(1.0, float(args.downsample['target_reads'])/float(total_seq_reads))
        martian.log_info("Input data downsampling: Requested: %.2f M reads, Estimated Input: %.2f M reads, Downsample Rate: %.3f" \
         % (float(args.downsample['target_reads'])/1e6, float(total_seq_reads)/1e6, global_subsample_rate))

        for chunk in chunks:
            chunk['subsample_rate'] = chunk['subsample_rate'] * global_subsample_rate



    martian.log_info("Input reads: %s" % str(chunks))
    outs.chunks = chunks
    outs.read_groups = [rg for rg in read_groups]

    # log info about input vs requested GB
    # first, set defaults
    available_gb = float(total_seq_bases)/1e9
    requested_gb = None
    available_reads = total_seq_reads
    requested_reads = None
    requested_rate = None
    post_downsample_gb = requested_gb
    downsample_succeeded = True

    if args.downsample is not None and args.downsample.get('gigabases') is not None:
        requested_gb = float(args.downsample['gigabases'])
        post_downsample_gb = min(available_gb, requested_gb)
        if available_gb < requested_gb:
            martian.log_info("Downsample requested more GB than was available; will not downsample.")
            downsample_succeeded = False

    elif args.downsample is not None and args.downsample.get('subsample_rate') is not None:
        requested_rate = float(args.downsample['subsample_rate'])
        post_downsample_gb = available_gb * requested_rate

    elif args.downsample is not None and args.downsample.get('target_reads') is not None:
        requested_reads = float(args.downsample['target_reads'])


    downsample_info = {}
    downsample_info['available_gb'] = available_gb
    downsample_info['requested_gb'] = requested_gb
    downsample_info['available_reads'] = available_reads
    downsample_info['requested_reads'] = requested_reads
    downsample_info['requested_rate'] = requested_rate
    downsample_info['post_downsample_gb'] = post_downsample_gb
    downsample_info['downsample_succeeded'] = downsample_succeeded

    with open(outs.downsample_info, 'w') as downsample_out:
        tenkit.safe_json.dump_numpy(downsample_info, downsample_out)

    check_fastqs(outs.chunks)

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

    lens = [len(l) for l in itertools.islice(reader, 20000)]
    uncompressed_sz = sum(lens)
    compressed_sz = reader.myfileobj.tell()

    # Rough guess of the uncompressed size
    estimated_full_sz = 1.1 * float(uncompressed_sz) / compressed_sz * file_sz
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
    return true_size


def fastq_data_estimate(fn, num_reads = 8000):
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

    if is_gz:
        #file_len = reader.myfileobj.tell()
        uncomp_size = estimate_gzip_uncompressed_size(fn) * 0.8        # STUPID FUDGE FACTOR TO GET >= REQUESTED AMT
    else:
        #file_len = data_len
        uncomp_size = file_sz

    read_yield = float(len(input_lens)) / total_data_len
    seq_yield = float(total_seq_len) / total_data_len
    predicted_reads = read_yield * uncomp_size
    predicted_seq = seq_yield * uncomp_size

    # For debugging
    #predicted_sz = float(total_data_len) / file_len * file_sz
    #gzip_sz = parse_gzip_sz(fn)
    #print "comp: %.2f, pred: %.2f, pred_mod2: %.2f, gzip_mod2: %.2f, gzip_est: %.2f" % (float(file_sz)/1e9, float(predicted_sz)/1e9, float(predicted_sz % 2**32)/1e9, float(gzip_sz)/1e9, float(uncomp_gzip_est2)/1e9)

    return (predicted_reads, predicted_seq)

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
