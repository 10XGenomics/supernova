#!/usr/bin/env python
#
# Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
#

import os.path
import tenkit.fasta as tk_fasta
import tenkit.preflight as tk_preflight
import martian
import gzip
import re
import socket
import tenkit.supernova as tk_sn

def check_key(n, dict_in, name, tys):
    if not dict_in.has_key(name):
        martian.exit("Entry %d in sample_def missing required field: %s" % (n, name))

    if not (type(dict_in[name]) in tys):
        martian.exit("Entry %d in sample_def for '%s' has incorrect type -- expecting %s, got %s" % (n, name, str(tys), type(dict_in[name])))



### main preflight checks for SN
SUPERNOVA_MEMORY_MB=512*10**3

def main(args, outs):
    hostname = socket.gethostname()
    tk_preflight.record_package_versions()
    
    ## no barcode whitelist
    if args.barcode_whitelist is None:
        martian.exit("No barcode whitelist specified.")
    
    ## there must be a barcode in each sample
    ## and it should be 16 bases long
    ## and it should be on read 1 or read 2
    for sd in args.sample_def:
        if sd.get("bc_length", 0) != 16 or sd.get("bc_in_read", 3) not in [1,2]:
            martian.exit("Barcode must be 16 bases and on read1 or read2.")


    print "Checking FASTQ folder..."
    for sample_def in args.sample_def:
        read_path = sample_def["read_path"]
        if not read_path:
            martian.exit("Must specify a read_path containing FASTQs.")
        if not read_path.startswith('/'):
            martian.exit("Specified FASTQ folder must be an absolute path: %s" % read_path)
        if not os.path.exists(read_path):
            martian.exit("On machine: %s, specified FASTQ folder does not exist: %s" % (hostname, read_path))
        if not os.access(read_path, os.X_OK):
            martian.exit("On machine: %s, supernova does not have permission to open FASTQ folder: %s" % (hostname, read_path))
        if not os.listdir(read_path):
            martian.exit("Specified FASTQ folder is empty: " + read_path)

        library_id = sample_def.get("library_id")
        if library_id is not None:
            if not re.match("^[\w-]+$", library_id):
                martian.exit("Library name may only contain letters, numbers, underscores, and dashes: " + library_id)

        lanes = sample_def["lanes"]
        if lanes is not None:
            for lane in lanes:
                if not is_int(lane):
                    martian.exit("Lanes must be a comma-separated list of numbers.")

    # Open file handles limit
    ok, msg = tk_preflight.check_open_fh()
    if not ok:
        martian.exit(msg)

    ## compile a list of fastq files
    fastq_files = []
    if args.input_mode == "BCL_PROCESSOR":
        # Validate the sample_def fields are correct
        for (idx, sample_item) in enumerate(args.sample_def):
            # validate
            check_key(idx, sample_item, "sample_indices", [list, type(None)])
            check_key(idx, sample_item, "read_path", [str, unicode])
            check_key(idx, sample_item, "lanes",  [list, type(None)])

        main_read_type = "RA"
        find_func = tk_fasta.find_input_fastq_files_10x_preprocess

        for read_chunk in args.sample_def:
            sample_index_strings, msg = tk_preflight.check_sample_indices(read_chunk)
            if sample_index_strings is None:
                martian.exit(msg)

            path = read_chunk['read_path']
            lanes = read_chunk['lanes']
            
            for sample_index in sample_index_strings:
                reads = find_func(path, main_read_type, sample_index, lanes)
                fastq_files.extend(reads)
    elif args.input_mode == "ILMN_BCL2FASTQ":
        # Validate the sample_def fields are correct
        for (idx, sample_item) in enumerate(args.sample_def):
            # validate
            check_key(idx, sample_item, "read_path", [str, unicode])
            check_key(idx, sample_item, "lanes",  [list, type(None)])
            check_key(idx, sample_item, "sample_names", [list, type(None)])

        find_func = tk_fasta.find_input_fastq_files_bcl2fastq_demult

        for read_chunk in args.sample_def:
            sample_names = read_chunk['sample_names']
            path = read_chunk['read_path']
            lanes = read_chunk['lanes']

            for sample_name in sample_names:
                reads = find_func(path, "R1", sample_name, lanes)
                fastq_files.extend(reads)
                reads = find_func(path, "R3", sample_name, lanes)
                fastq_files.extend(reads)
    else:
        martian.throw("Unrecognized input_mode: %s" % args.input_mode)

    ## if we found nothing then break
    if len(fastq_files) == 0:
        martian.exit("No input FASTQs were found with the requested lanes and sample indices.")
    
    ## make sure they are okay first
    check_fastqs(fastq_files)
    
    total_reads = 0.0
    global_avg  = 0.0
    num_files   = 0
    for fn in fastq_files:
        reads_fn, avg_read_len_fn = estimate_read_count_and_length(fn, num_reads=1000)
        total_reads += reads_fn
        global_avg  += avg_read_len_fn
        num_files   += 1
    global_avg = global_avg/num_files
    martian.log_info("Estimated read length = %.1f, Estimated total read input = %.1f"% (global_avg, total_reads))

    PreflightAlert = tk_sn.AlertLogger(stage="preflight")
    PreflightAlert.issue("mean_read_length", global_avg)
     
    ### if the user has supplied a target number of reads
    ### use that in the memory check.
    #if args.downsample is not None:
    #    target_reads = args.downsample.get("target_reads", None)
    #    if target_reads is not None and target_reads < total_reads:
    #        total_reads = target_reads
    #
    ### now do a read count estimate
    ### and perform a memory check
    #mem = psutil.virtual_memory().available/1e6 ## in MB
    #mem_reqd = SUPERNOVA_MEMORY_MB/1.2e9*total_reads
    #if mem < 0:
    #    martian.log_warn("couldn't check available memory")
    #else:
    #    martian.log_info("available mem = %.2f, mem required (estimate for %d reads) =  %.2f " % (mem, total_reads, mem_reqd))
    #    if mem < mem_reqd:
    #        martian.exit("Insufficient memory to run Supernova on the supplied input: %.1f MB available < %.1f MB required. If you would like to proceed anyway turn off the preflight check." % (mem, mem_reqd) )


def check_fastqs(fastq_files):
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

    for fastq in fastq_files:
        check_fastq(fastq)

def estimate_read_count_and_length(fn, num_reads=1000):
    ''' 
    Estimate the number of reads AND the average read length
    in the fastq file fn by only reading in the first
    num_reads (default 1000) reads.
    '''
    # Open reader
    if fn[-2:] == 'gz':
        reader = gzip.open(fn)
        is_gz = True
    else:
        reader = open(fn, 'r')
        is_gz = False
    ## first compute the average read length
    avg_read_length = 0.0    
    gen = tk_fasta.read_generator_fastq(reader)        
    rec_count = 0
    for (header, r, qual) in gen:
        avg_read_length += len(r)
        rec_count += 1
        if rec_count == num_reads:
            break
    avg_read_length = avg_read_length/rec_count
    if is_gz:
        file_len = reader.myfileobj.tell() 
    else:
        file_len = reader.tell()
    ## total file size 
    file_sz = os.path.getsize(fn)
    total_reads_est= float(num_reads)/file_len*file_sz
    return (total_reads_est, avg_read_length)

def is_int(s):
    try:
        int(s)
    except ValueError:
        return False
    return True


 
