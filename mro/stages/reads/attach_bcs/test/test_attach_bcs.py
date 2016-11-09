#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#

# Test attach bcs
import pysam
from itertools import groupby
import tenkit.test as tk_test
import tenkit.seq as tk_seq
from .. import *
from tenkit.constants import PROCESSED_BARCODE_TAG, RAW_BARCODE_TAG, SAMPLE_INDEX_TAG

import martian

IN_BAM = tk_test.in_path('attach_bcs/alignment_with_secondary.bam')
IN_I1 = tk_test.in_path('attach_bcs/phix_I1.fastq')
IN_I2 = tk_test.in_path('attach_bcs/phix_I2.fastq')
IN_WHITELIST = '737K-april-2014'

OUT_BAM = tk_test.out_path('test_attach_bcs.bam')

class TestFunctions(tk_test.UnitTestBase):
    def setUp(self):
        pass

    def test_attach_bcs(self):
        #  --align_input alignment_output.bam --barcode_input phix_I2.fastq --output test2.out --complete ~/c --stats ~/s
        args = {
            'barcode_whitelist' : IN_WHITELIST,
            'align_chunk' : IN_BAM,
            'barcode_chunk' : IN_I2,
            'sample_index_chunk' : IN_I1,
            'gem_group' : None,
            'paired_end' : True,
            'exclude_non_bc_reads' : False,
            'max_expected_bc_error': 0.75,
            'subsample_rate' : 1.0,
        }
        outs = { 'output': OUT_BAM }

        args = martian.Record(args)
        outs = martian.Record(outs)

        main(args, outs)

        # Get the barcodes
        barcode_whitelist = tk_seq.load_barcode_whitelist(IN_WHITELIST)

        # Ensure each read has a barcode
        out_bam = pysam.Samfile(OUT_BAM)
        for r in out_bam:
            tag_dict = { k:v for (k,v) in r.tags }
            tag_names = [ k for (k,v) in r.tags ]
            self.assertTrue(RAW_BARCODE_TAG in tag_names)

            if tag_dict[RAW_BARCODE_TAG] in barcode_whitelist:
                self.assertTrue(PROCESSED_BARCODE_TAG in tag_names)

            self.assertTrue(SAMPLE_INDEX_TAG in tag_names)


        # Make sure we put out the full BAM file
        out_len = len([ x for x in pysam.Samfile(OUT_BAM)])
        in_len  = len([ x for x in pysam.Samfile(IN_BAM)])
        self.assertEqual(out_len, in_len)


        def get_bc(r):
            tags = { k:v for (k,v) in r.tags }
            return tags[RAW_BARCODE_TAG]

        # Ensure each read pair has the same barcode
        out_bam = pysam.Samfile(OUT_BAM)
        reads = [ x for x in out_bam ]

        for (grp, reads) in groupby(reads, lambda x: x.qname):
            bcs = set(tk_io.get_read_barcode(r) for r in reads)
            self.assertEqual(len(bcs), 1)



if __name__ == '__main__':
    tk_test.run_tests()
