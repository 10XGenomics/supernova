#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#

# Test attach bcs
import pysam
import tenkit.test as tk_test
import tenkit.fasta as tk_fasta
from .. import *

import martian

IN_FASTQ = tk_test.in_path("test_bwa.fastq")
OUT_BAM = tk_test.out_path("test_aligner.bam")

class TestFunctions(tk_test.UnitTestBase):
    def setUp(self):
        pass

    def test_align(self):
        args = { 'chunk_input': IN_FASTQ, 'aligner':'bwa', 'aligner_method': 'MEM', 'reference_path':'hg19', '__threads': 1, 'reads_interleaved': True,
                'read_group_sample': 'test_sample'}
        outs = { 'default': OUT_BAM }

        args = martian.Record(args)
        outs = martian.Record(outs)

        main(args, outs)

        # Ensure each read has a barcode
        out_bam = pysam.Samfile(OUT_BAM)
        bam_reads = list(out_bam)

        fq_file = open(IN_FASTQ)
        fq_reads = list(tk_fasta.read_generator_fastq(fq_file, paired_end=False))

        self.assertEqual(len(bam_reads), len(fq_reads))
