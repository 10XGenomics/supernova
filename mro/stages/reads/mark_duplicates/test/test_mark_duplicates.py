#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#
import tenkit.test as tk_test
from .. import *
import tenkit.bio_io as tk_io
import tenkit.seq as tk_seq
import pysam
import itertools
import numpy as np
import tenkit.constants
import tenkit.read_filter

import martian

martian.test_initialize("")

# Patrick Marks
# Simple test of deduper

IN_BAM = tk_test.in_path('test_analyze_bias.bam')
OUT_BAM = tk_test.out_path('test_dedup_out.bam')
OUT_JSON = tk_test.out_path( 'test_dedup_summary.json')
IN_BAM_BIG = tk_test.in_path("test_mark_duplicates.bam")

class TestFunctions(tk_test.UnitTestBase):
    def setUp(self):
        pass

    def test_dedup(self):
        tenkit.constants.DUPLICATE_SUBSAMPLE_COVERAGES = [0.00001, 0.0001]
        args = martian.Record({ 'input': IN_BAM, 'estimated_coverage': 100.0, 'perfect_read_count': 1000, 'chunk_start':None, 'chunk_end':None })
        outs = martian.Record({ 'output': OUT_BAM, 'duplicate_summary': OUT_JSON })
        main_mark_duplicates(args, outs)

        out_bam = pysam.Samfile(OUT_BAM)
        dups = [ x.is_duplicate for x in out_bam ]

        self.assertEqual(dups, [ False, True, False, False, True, False ])


    def test_big_dedup(self):
        tenkit.constants.DUPLICATE_SUBSAMPLE_COVERAGES = [0.000003, 0.000015]
        args = martian.Record({ 'input': IN_BAM_BIG, 'estimated_coverage':100.0, 'perfect_read_count': 100000, 'chunk_start': None, 'chunk_end': None })
        outs = martian.Record({ 'output': OUT_BAM, 'duplicate_summary': OUT_JSON })
        main_mark_duplicates(args, outs)

        out_bam = pysam.Samfile(OUT_BAM)
        out_reads = list(out_bam)

        in_bam = pysam.Samfile(IN_BAM_BIG)
        in_reads = list(in_bam)

        # Check we haven't lost any reads
        self.assertEqual(len(out_reads), len(in_reads))

        def read_tuple(r):
            bc = tk_io.get_read_barcode(r)
            return (bc, r.tid, r.pos, r.mrnm, r.mpos, r.is_reverse, r.is_read1)
            #return (bc, r.is_read1, r.is_reverse, r.tid, r.pos, r.mrnm, r.mpos)

        def mark_duplicates(read_set):

            # Re-run the dup analysis manually
            read_tups = [(read_tuple(r), r) for r in read_set]
            read_tups.sort(key = lambda x: x[0])

            groups = itertools.groupby(read_tups, lambda x: x[0])

            for (k, reads) in groups:
                rl = list(reads)
                rl[0][1].is_duplicate = False

                for i in range(1, len(rl)):
                    rl[i][1].is_duplicate = True


        mark_duplicates(in_reads)

        # Make sure our 'all-reads' analysis matches the code
        out_dup_marks = np.array([ r.is_duplicate for r in out_reads if (not r.is_unmapped) and (not r.mate_is_unmapped)])
        test_dup_marks = np.array([ r.is_duplicate for r in in_reads if (not r.is_unmapped) and (not r.mate_is_unmapped)])

        print "len(start_bam): %d  -- len(out_bam): %d" % (len(out_dup_marks), len(test_dup_marks))
        eq = (out_dup_marks == test_dup_marks).all()

        print "mean dups code: %f" % out_dup_marks.mean()
        print "mean dups test: %f" % test_dup_marks.mean()

        self.assertTrue(eq)

        # Read the molecule count histogram and verify
        count_hist = json.load(file(OUT_JSON))['no_filter_full_use_bcs']

        dups = sum([ (int(times_observed) - 1) * n for (times_observed, n) in count_hist.items() ])
        total_reads = sum([ int(times_observed) * n for (times_observed, n) in count_hist.items() ])
        summary_dup_rate = float(dups) / total_reads

        mapped_in_reads = np.array([r.is_duplicate for r in in_reads if not(r.is_unmapped or r.mate_is_unmapped) and tk_io.get_read_barcode(r) is not None ])
        self.assertEqual(summary_dup_rate, mapped_in_reads.mean())


        # Get the perfect reads, mark dups and compare stats
        perfect_reads = [x for x in in_reads if tenkit.read_filter.stringent_read_filter(x, True)]
        mark_duplicates(perfect_reads)

        # Read the molecule count histogram and verify -- perfect reads
        count_hist = json.load(file(OUT_JSON))['full_use_bcs']

        dups = sum([ (int(times_observed) - 1) * n for (times_observed, n) in count_hist.items() ])
        total_reads = sum([ int(times_observed) * n for (times_observed, n) in count_hist.items() ])
        summary_dup_rate = float(dups) / total_reads

        mapped_in_reads = np.array([r.is_duplicate for r in perfect_reads])
        self.assertEqual(summary_dup_rate, mapped_in_reads.mean())






