#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# Runs alignment code
#
import itertools
import tenkit.align
import tenkit.fasta as tk_fasta
import tenkit.reference
import martian
import math

__MRO__ = """
stage ALIGN(
    in  map[] chunks,
    in  string aligner,
    in  string aligner_method,
    in  string reference_path,
    in  string read_group_sample,
    in  int    num_threads,
    out bam,
    src py     "stages/reads/align_reads",
) split using (
    in map chunk,
)
"""

# set a min read length -- BWA-MEM won't work below this
MIN_READ_LENGTH = 25

def split(args):
    '''We just align each chunk independently -- joining will happen in the join step of SORT_READS'''

    # Pull some reads from fastq files -- bail out if it's less than 25bp
    fastq_tests = [x['read1'] for x in args.chunks]

    for fastq_test in fastq_tests:
        with open(fastq_test) as in_file:
            reader = tk_fasta.read_generator_fastq(in_file)
            for name, read, qual in itertools.islice(reader, 10):
                continue
                if len(read) < MIN_READ_LENGTH:
                    martian.alarm("BWA-MEM can't handle reads <25bp -- reads will be unmapped.")
                    continue

    # estimated amount of memory needed to process genome is 2x(num gigabases)+4GB
    reference_pyfasta = tenkit.reference.open_reference(args.reference_path)
    reference_bases = sum(len(reference_pyfasta[contig]) for contig in reference_pyfasta)
    base_mem_in_gb = int(math.ceil(2*reference_bases / (1024.0**3)))

    mem_in_gb = base_mem_in_gb + 4
    chunks = [{'chunk':x, '__threads': args.num_threads, '__mem_gb': mem_in_gb} for x in args.chunks]
    return {'chunks': chunks}

def main(args, outs):
    chunk = args.chunk

    if not chunk['reads_interleaved'] and (chunk['read1'] is None or chunk['read2'] is None):
        martian.throw("must supply a read1 and read2 when reads_interleave == False")

    if chunk['reads_interleaved']:
        reads = chunk['read1']
    else:
        reads = [chunk['read1']]
        if chunk['read2'] is not None:
            reads.append(chunk['read2'])

    a = tenkit.align.Aligner(reads, outs.default)
    aligner = args.aligner

    ref_fasta = tenkit.reference.get_fasta(args.reference_path)
    a.output_alignment(aligner=aligner, aligner_params={'ref_fasta':ref_fasta, 'algorithm': args.aligner_method}, num_threads=args.__threads, sample=args.read_group_sample)

def join(args, outs, chunk_defs, chunk_outs):
    outs.default = [chunk.default for chunk in chunk_outs]
