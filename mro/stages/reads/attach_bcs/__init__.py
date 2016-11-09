#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# Put the sample index barcode and 10X barcode into read tags in a BAM file.
# Screen each 10X barcode against the given barcode whitelist
#
import gzip
import itertools
import random
import numpy as np
import array
import json
import tenkit.seq as tk_seq
import tenkit.bam as tk_bam
import tenkit.fasta as tk_fasta
import tenkit.bio_io as tk_io
import tenkit.read_filter
from tenkit.constants import PROCESSED_BARCODE_TAG, RAW_BARCODE_TAG, RAW_BARCODE_QUAL_TAG, SAMPLE_INDEX_TAG, SAMPLE_INDEX_QUAL_TAG, ILLUMINA_QUAL_OFFSET


import martian

__MRO__ = """
 stage ATTACH_BCS(
     in  string barcode_whitelist,
     in  bam    align,
     in  map[]  chunks,
     in  bool   paired_end,
     in  bool   exclude_non_bc_reads,
     in  float  bc_confidence_threshold,
     in  json   bc_counts,
     out bam    output,
     out int    perfect_read_count,
     src py     "stages/reads/attach_bcs",
 ) split using (
     in  bam    align_chunk,
     in  map    chunk,
 )
"""

def split(args):
    ''' Attach BCS to each chunk of the input files '''

    chunk_defs = [{'chunk': fastq_chunk, 'align_chunk': aln_chunk} for (fastq_chunk, aln_chunk) in zip(args.chunks, args.align)]
    return {'chunks': chunk_defs}



def join(args, outs, chunk_defs, chunk_outs):
    ''' Pass through the BAM files for sorting '''
    outs.output = [chunk.output for chunk in chunk_outs]
    outs.perfect_read_count = sum([chunk.perfect_read_count for chunk in chunk_outs])


def main(args, outs):
    """ Attaches barcodes. Attaches raw barcode to RAW_BC tag and filters those to form set of PROCESSES_BARCODES """

    chunk = args.chunk

    #subsample_rate = 1.0
    #if args.subsample_rate is not None:
    #    subsample_rate = args.subsample_rate

    bam_in = tk_bam.create_bam_infile(args.align_chunk)
    bam_out, tids = tk_bam.create_bam_outfile(outs.output, None, None, template=bam_in, pgs=tk_bam.make_pg_header(martian.get_pipelines_version(), "attach_bcs"))

    if args.barcode_whitelist is None or args.bc_counts is None:
        # If there's no whitelist or counts then all high quality BC reads get allowed.
        barcode_whitelist = None
        wl_idxs = None
        bc_dist = None
    else:
        barcode_whitelist = tk_seq.load_barcode_whitelist(args.barcode_whitelist)

        # Load the bc counts for this GEM group
        counts = json.load(open(args.bc_counts, 'r'))
        counts = counts[str(chunk['gem_group'])]['bc_counts']

        # Prior distribution over barcodes, with pseudo-count
        bc_dist = np.array(counts, dtype=np.float) + 1.0
        bc_dist = bc_dist / bc_dist.sum()
        wl_idxs = { bc:idx for (idx,bc) in enumerate(sorted(list(barcode_whitelist))) }

    # set random seed to get deterministic subsampling
    random.seed(0)

    def open_maybe_gzip(fn):
        if fn[-2:] == "gz":
            return gzip.open(fn)
        else:
            return open(fn)

    if chunk['barcode']:
        processed_barcode_iter = get_raw_processed_barcodes(open_maybe_gzip(chunk['barcode']), barcode_whitelist, args.bc_confidence_threshold, chunk['gem_group'], chunk['barcode_reverse_complement'], wl_idxs, bc_dist)
        require_barcode_for_stringent = True
    else:
        processed_barcode_iter = itertools.repeat(None)
        require_barcode_for_stringent = False

    if chunk['sample_index']:
        sample_index_iter = tk_fasta.read_generator_fastq(open_maybe_gzip(chunk['sample_index']))
    else:
        sample_index_iter = itertools.repeat(None)

    iters = itertools.izip(processed_barcode_iter, sample_index_iter)

    # First read
    read = bam_in.next()

    # Number of perfect reads -- used to compute down-sampling rates in mark_duplicates
    perfect_read_count = 0

    # Due to secondary alignments, we must apply the tags to all
    # reads with the same cluster name.
    for (barcode_info, sample_index_info) in iters:
        tags = []
        read_name = None

        if read is None:
            break

        if barcode_info:
            (bc_read_name, raw_bc_seq, processed_bc_seq, raw_bc_qual) = barcode_info
            tags.append((RAW_BARCODE_TAG, raw_bc_seq))
            tags.append((RAW_BARCODE_QUAL_TAG, raw_bc_qual))
            if processed_bc_seq is not None:
                tags.append((PROCESSED_BARCODE_TAG, processed_bc_seq))
            read_name = bc_read_name.split()[0]


        if sample_index_info:
            (si_read_name, seq, qual) = sample_index_info
            tags.append((SAMPLE_INDEX_TAG, seq))
            tags.append((SAMPLE_INDEX_QUAL_TAG, qual))

            if read_name != None:
                if si_read_name.split()[0] != read_name:
                    martian.log_info("mismatch: si_read_name: %s, bam_read_name: %s" % (si_read_name, read_name))
                assert(si_read_name.split()[0] == read_name)
            else:
                read_name = si_read_name.split()[0]

        reads_attached = 0
        #emit_read_pair = random.random() < subsample_rate
        emit_read_pair = True

        while read.qname == read_name or read_name == None:
            if len(tags) > 0:
                existing_tags = read.tags
                existing_tags.extend(tags)
                read.tags = existing_tags

            reads_attached += 1
            if not (read_name is None):
                assert(read.qname == read_name)

            if emit_read_pair:
                # Count the perfect reads -- will be used when subsampling in dedup
                if tenkit.read_filter.stringent_read_filter(read, require_barcode_for_stringent):
                    perfect_read_count += 1

                if args.exclude_non_bc_reads:
                    if not(tk_io.get_read_barcode(read) is None):
                        bam_out.write(read)
                else:
                    bam_out.write(read)

            try:
                read = bam_in.next()

            except StopIteration:
                read = None
                break

        # We may have more than 2 reads is there was a
        # secondary alignment, but less than 2 means
        # something went wrong
        assert(reads_attached >= 2)


    outs.perfect_read_count = perfect_read_count
    bam_out.close()


def get_raw_processed_barcodes(barcode_file, barcode_whitelist, bc_confidence_threshold, gem_group, barcodes_reverse_complement, wl_idxs, wl_dist):
    """ Stream the barcodes and the 'processed' barcode """
    bc_iterator = tk_fasta.read_generator_fastq(barcode_file)

    gem_group_str = "-" + str(gem_group)

    for (name, seq, qual) in bc_iterator:
        if barcodes_reverse_complement:
            seq = tk_seq.get_rev_comp(seq)
            qual = qual[::-1] #reverse qual string
        # Check for valid bc sequences
        if barcode_whitelist is None:
            # No whitelist case -- attach BC if there are no Ns
            if not ('N' in seq):
                processed_bc = seq + gem_group_str
                yield (name, seq, processed_bc, qual)
            else:
                yield (name, seq, None, qual)
        else:
            # whitelist case -- attach bc if posterior probability of best
            # BC sequence exceeds the confidence threshold
            bc_seq = handle_10x_barcode(bc_confidence_threshold, seq, qual, wl_idxs, wl_dist)
            if bc_seq is None:
                yield (name, seq, None, qual)
            else:
                processed_bc = bc_seq + gem_group_str
                yield (name, seq, processed_bc, qual)

def handle_10x_barcode(bc_confidence_threshold, seq, qual, wl_idxs, wl_dist):
    '''Compute the MAP whitelist BC sequence using the read QV * bc distribution
       model. If the MAP probability exceeds the bc_confidence_threshold,
       return that sequence, otherwise return None.  In the case the BC is on
       the whitelist, consider Hamming distance=2 neighbors, if not, consider
       Hamming distance=1 neighbors'''

    if seq in wl_idxs:
        return check_correct_bc(bc_confidence_threshold, seq, qual, wl_idxs, wl_dist)
    else:
        return correct_bc_error(bc_confidence_threshold, seq, qual, wl_idxs, wl_dist)


def correct_bc_error(bc_confidence_threshold, seq, qual, wl_idxs, wl_dist):
    '''Attempt to correct an incorrect BC sequence by computing
       the probability that a Hamming distance=1 BC generated
       the observed sequence, accounting for the prior distribution
       of the whitelist barcodes (wl_dist), and the QV of the base
       that must have been incorrect'''

    # QV values
    qvs = np.fromstring(qual, dtype=np.byte) - ILLUMINA_QUAL_OFFSET

    # Char array of read
    a = array.array('c', seq)

    # Likelihood of candidates
    wl_cand = []
    likelihoods = []

    # Enumerate Hamming distance 1 sequences - if a sequence
    # is on the whitelist, compute it's likelihood.
    for pos in range(len(a)):
        existing = a[pos]
        for c in ['A', 'C', 'G', 'T']:
            if c == existing:
                continue
            a[pos] = c
            test_str = a.tostring()

            idx = wl_idxs.get(test_str)
            if idx is not None:
                # prior probability of this BC
                p_bc = wl_dist[idx]

                # probability of the base error
                edit_qv = min(33.0, float(qvs[pos]))
                p_edit = 10.0**(-edit_qv/10.0)
                wl_cand.append(test_str)
                likelihoods.append(p_bc * p_edit)

        a[pos] = existing

    posterior = np.array(likelihoods)
    posterior /= posterior.sum()
    if len(posterior) > 1:
        pmax = posterior.max()
        if pmax > bc_confidence_threshold:
            return wl_cand[np.argmax(posterior)]

    return None


def check_correct_bc(bc_confidence_threshold, seq, qual, wl_idxs, wl_dist):
    '''Attempt to correct an incorrect BC sequence by computing
       the probability that a Hamming distance=1 BC generated
       the observed sequence, accounting for the prior distribution
       of the whitelist barcodes (wl_dist), and the QV of the base
       that must have been incorrect'''

    # QV values
    qvs = np.fromstring(qual, dtype=np.byte) - ILLUMINA_QUAL_OFFSET

    # Only examine the barcode if there is a QV <= 24
    if (qvs > 24).all():
        return seq

    # Char array of read
    a = array.array('c', seq)

    # Likelihood of candidates
    wl_cand = []
    likelihoods = []

    # include the un-edited case here -- this assumes the pErrs are small
    original_bc = a.tostring()
    wl_cand.append(original_bc)
    likelihoods.append(wl_dist[wl_idxs[original_bc]])

    # Enumerate Hamming distance 1 sequences - if a sequence
    # is on the whitelist, compute it's likelihood.
    for pos in range(len(a)):
        existing = a[pos]
        for c in ['A', 'C', 'G', 'T']:
            if c == existing:
                continue
            a[pos] = c

            for pos2 in range(pos+1, len(a)):
                existing2 = a[pos2]
                for c2 in ['A', 'C', 'G', 'T']:
                    if c2 == existing2:
                        continue
                    a[pos2] = c2
                    test_str = a.tostring()

                    idx = wl_idxs.get(test_str)
                    if idx is not None:
                        # prior probability of this BC
                        p_bc = wl_dist[idx]

                        # probability of the base errors
                        edit_qv1 = min(33.0, max(3.0, float(qvs[pos])-1.0))
                        edit_qv2 = min(33.0, max(3.0, float(qvs[pos2])-1.0))
                        p_edit = 10.0**(-edit_qv1/10.0) * 10.0**(-edit_qv2/10.0)
                        likelihoods.append(p_bc * p_edit)

                        wl_cand.append(test_str)

                a[pos2] = existing2
        a[pos] = existing

    posterior = np.array(likelihoods)
    posterior /= posterior.sum()
    if len(posterior) > 1:
        pmax = posterior.max()
        if pmax > bc_confidence_threshold:
            return wl_cand[np.argmax(posterior)]

    return None
