#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#
import os
import tenkit.test as tk_test
import tenkit.fasta as tk_fasta
import martian
from .. import *

martian.test_initialize(tk_test.out_path(""))

IN_PREFIX = tk_test.in_path('combine_and_trim_reads')
IN_RA_ALL = [ tk_test.in_path('combine_and_trim_reads/read-RA_si-AAAA_lane-1_chunk-1.fastq'),
              tk_test.in_path('combine_and_trim_reads/read-RA_si-AAAN_lane-1_chunk-1.fastq'),
              tk_test.in_path('combine_and_trim_reads/read-RA_si-CCCC_lane-1_chunk-1.fastq') ]

class TestFunctions(tk_test.UnitTestBase):

    def test_setup_chunks(self):

        args = martian.Record({
            'input_mode': 'BCL_PROCESSOR',
            'sample_def': [ { 'read_path': IN_PREFIX,
                              'sample_indices': ["AAAA", "CCCC"],
                              'lanes': None,
                              'gem_group': None, } ],
            'barcode_whitelist':"737K-april-2014",
            }
            )
        outs = martian.Record({})

        main(args, outs)
        print outs
        self.assertTrue(len(outs.chunks) == 3)

    def test_ilmn_bclfastq_mode(self):
        args = martian.Record({
            'input_mode': "ILMN_BCL2FASTQ",
            'sample_def': [{
                            "gem_group": None,
                            "lanes": None,
                            "read_path": "/mnt/projects/pat/bcl_direct/t2/p1/a",
                            "samples": [
                                                "aa"
                                            ]
                            },
                            {
                            "gem_group": None,
                            "lanes": None,
                            "read_path": "/mnt/projects/pat/bcl_direct/t2/p2/b",
                            "samples": [
                                           "bb"
                                       ]
                            }],
            'barcode_whitelist': "737K-april-2014",
        }
        )
        outs = martian.Record({})
        if os.path.exists("/mnt/projects/pat/bcl_direct/t2/p1/a"):
            main(args,outs)
            self.assertTrue(outs.barcodes == [
                        "/mnt/projects/pat/bcl_direct/t2/p1/a/aa_S1_L001_R2_001.fastq.gz",
                        "/mnt/projects/pat/bcl_direct/t2/p1/a/aa_S1_L002_R2_001.fastq.gz",
                        "/mnt/projects/pat/bcl_direct/t2/p1/a/aa_S1_L003_R2_001.fastq.gz",
                        "/mnt/projects/pat/bcl_direct/t2/p1/a/aa_S1_L004_R2_001.fastq.gz",
                        "/mnt/projects/pat/bcl_direct/t2/p1/a/aa_S1_L005_R2_001.fastq.gz",
                        "/mnt/projects/pat/bcl_direct/t2/p1/a/aa_S1_L006_R2_001.fastq.gz",
                        "/mnt/projects/pat/bcl_direct/t2/p1/a/aa_S1_L007_R2_001.fastq.gz",
                        "/mnt/projects/pat/bcl_direct/t2/p1/a/aa_S1_L008_R2_001.fastq.gz",
                        "/mnt/projects/pat/bcl_direct/t2/p2/b/bb_S2_L001_R2_001.fastq.gz",
                        "/mnt/projects/pat/bcl_direct/t2/p2/b/bb_S2_L002_R2_001.fastq.gz",
                        "/mnt/projects/pat/bcl_direct/t2/p2/b/bb_S2_L003_R2_001.fastq.gz",
                        "/mnt/projects/pat/bcl_direct/t2/p2/b/bb_S2_L004_R2_001.fastq.gz",
                        "/mnt/projects/pat/bcl_direct/t2/p2/b/bb_S2_L005_R2_001.fastq.gz",
                        "/mnt/projects/pat/bcl_direct/t2/p2/b/bb_S2_L006_R2_001.fastq.gz",
                        "/mnt/projects/pat/bcl_direct/t2/p2/b/bb_S2_L007_R2_001.fastq.gz",
                        "/mnt/projects/pat/bcl_direct/t2/p2/b/bb_S2_L008_R2_001.fastq.gz"
                    ])
            self.assertTrue(outs.barcodes_reverse_complement[0] == True)
            self.assertTrue(outs.reads_interleaved == False)
            self.assertTrue(outs.sample_indices == [
                        "/mnt/projects/pat/bcl_direct/t2/p1/a/aa_S1_L001_I1_001.fastq.gz",
                        "/mnt/projects/pat/bcl_direct/t2/p1/a/aa_S1_L002_I1_001.fastq.gz",
                        "/mnt/projects/pat/bcl_direct/t2/p1/a/aa_S1_L003_I1_001.fastq.gz",
                        "/mnt/projects/pat/bcl_direct/t2/p1/a/aa_S1_L004_I1_001.fastq.gz",
                        "/mnt/projects/pat/bcl_direct/t2/p1/a/aa_S1_L005_I1_001.fastq.gz",
                        "/mnt/projects/pat/bcl_direct/t2/p1/a/aa_S1_L006_I1_001.fastq.gz",
                        "/mnt/projects/pat/bcl_direct/t2/p1/a/aa_S1_L007_I1_001.fastq.gz",
                        "/mnt/projects/pat/bcl_direct/t2/p1/a/aa_S1_L008_I1_001.fastq.gz",
                        "/mnt/projects/pat/bcl_direct/t2/p2/b/bb_S2_L001_I1_001.fastq.gz",
                        "/mnt/projects/pat/bcl_direct/t2/p2/b/bb_S2_L002_I1_001.fastq.gz",
                        "/mnt/projects/pat/bcl_direct/t2/p2/b/bb_S2_L003_I1_001.fastq.gz",
                        "/mnt/projects/pat/bcl_direct/t2/p2/b/bb_S2_L004_I1_001.fastq.gz",
                        "/mnt/projects/pat/bcl_direct/t2/p2/b/bb_S2_L005_I1_001.fastq.gz",
                        "/mnt/projects/pat/bcl_direct/t2/p2/b/bb_S2_L006_I1_001.fastq.gz",
                        "/mnt/projects/pat/bcl_direct/t2/p2/b/bb_S2_L007_I1_001.fastq.gz",
                        "/mnt/projects/pat/bcl_direct/t2/p2/b/bb_S2_L008_I1_001.fastq.gz"
                    ])
            self.assertTrue(outs.is_read1[0] == True and outs.is_read1[1] == False)

    def test_ilmn_bcl2fastq_t1(self):
        args = martian.Record({
            'input_mode': "ILMN_BCL2FASTQ",
            'sample_def': [{
                        "gem_group": None,
                        "lanes": None,
                        "read_path": "/mnt/projects/pat/bcl_direct/t1/p1",
                        "samples": [
                                        "a","b","c"
                                    ]
                    }],
            'barcode_whitelist': "737K-april-2014",
        }
        )
        outs = martian.Record({})
        if os.path.exists("/mnt/projects/pat/bcl_direct/t1/p1"):
            main(args,outs)
            self.assertTrue(outs.barcodes == [
                        "/mnt/projects/pat/bcl_direct/t1/p1/a_S1_L001_R2_001.fastq.gz",
                        "/mnt/projects/pat/bcl_direct/t1/p1/a_S1_L002_R2_001.fastq.gz",
                        "/mnt/projects/pat/bcl_direct/t1/p1/a_S1_L003_R2_001.fastq.gz",
                        "/mnt/projects/pat/bcl_direct/t1/p1/a_S1_L004_R2_001.fastq.gz",
                        "/mnt/projects/pat/bcl_direct/t1/p1/a_S1_L005_R2_001.fastq.gz",
                        "/mnt/projects/pat/bcl_direct/t1/p1/a_S1_L006_R2_001.fastq.gz",
                        "/mnt/projects/pat/bcl_direct/t1/p1/a_S1_L007_R2_001.fastq.gz",
                        "/mnt/projects/pat/bcl_direct/t1/p1/a_S1_L008_R2_001.fastq.gz",
                        "/mnt/projects/pat/bcl_direct/t1/p1/b_S2_L001_R2_001.fastq.gz",
                        "/mnt/projects/pat/bcl_direct/t1/p1/b_S2_L002_R2_001.fastq.gz",
                        "/mnt/projects/pat/bcl_direct/t1/p1/b_S2_L003_R2_001.fastq.gz",
                        "/mnt/projects/pat/bcl_direct/t1/p1/b_S2_L004_R2_001.fastq.gz",
                        "/mnt/projects/pat/bcl_direct/t1/p1/b_S2_L005_R2_001.fastq.gz",
                        "/mnt/projects/pat/bcl_direct/t1/p1/b_S2_L006_R2_001.fastq.gz",
                        "/mnt/projects/pat/bcl_direct/t1/p1/b_S2_L007_R2_001.fastq.gz",
                        "/mnt/projects/pat/bcl_direct/t1/p1/b_S2_L008_R2_001.fastq.gz",
                        "/mnt/projects/pat/bcl_direct/t1/p1/c_S3_L001_R2_001.fastq.gz",
                        "/mnt/projects/pat/bcl_direct/t1/p1/c_S3_L002_R2_001.fastq.gz",
                        "/mnt/projects/pat/bcl_direct/t1/p1/c_S3_L003_R2_001.fastq.gz",
                        "/mnt/projects/pat/bcl_direct/t1/p1/c_S3_L004_R2_001.fastq.gz",
                        "/mnt/projects/pat/bcl_direct/t1/p1/c_S3_L005_R2_001.fastq.gz",
                        "/mnt/projects/pat/bcl_direct/t1/p1/c_S3_L006_R2_001.fastq.gz",
                        "/mnt/projects/pat/bcl_direct/t1/p1/c_S3_L007_R2_001.fastq.gz",
                        "/mnt/projects/pat/bcl_direct/t1/p1/c_S3_L008_R2_001.fastq.gz"
                    ])
            self.assertTrue(outs.barcodes_reverse_complement[0] == True)
            self.assertTrue(outs.reads_interleaved == False)
            self.assertTrue(outs.sample_indices[5] == None)
            self.assertTrue(outs.is_read1[0] == True and outs.is_read1[1] == False)
            self.assertTrue(outs.reads == [
                    "/mnt/projects/pat/bcl_direct/t1/p1/a_S1_L001_R1_001.fastq.gz",
                    "/mnt/projects/pat/bcl_direct/t1/p1/a_S1_L001_R3_001.fastq.gz",
                    "/mnt/projects/pat/bcl_direct/t1/p1/a_S1_L002_R1_001.fastq.gz",
                    "/mnt/projects/pat/bcl_direct/t1/p1/a_S1_L002_R3_001.fastq.gz",
                    "/mnt/projects/pat/bcl_direct/t1/p1/a_S1_L003_R1_001.fastq.gz",
                    "/mnt/projects/pat/bcl_direct/t1/p1/a_S1_L003_R3_001.fastq.gz",
                    "/mnt/projects/pat/bcl_direct/t1/p1/a_S1_L004_R1_001.fastq.gz",
                    "/mnt/projects/pat/bcl_direct/t1/p1/a_S1_L004_R3_001.fastq.gz",
                    "/mnt/projects/pat/bcl_direct/t1/p1/a_S1_L005_R1_001.fastq.gz",
                    "/mnt/projects/pat/bcl_direct/t1/p1/a_S1_L005_R3_001.fastq.gz",
                    "/mnt/projects/pat/bcl_direct/t1/p1/a_S1_L006_R1_001.fastq.gz",
                    "/mnt/projects/pat/bcl_direct/t1/p1/a_S1_L006_R3_001.fastq.gz",
                    "/mnt/projects/pat/bcl_direct/t1/p1/a_S1_L007_R1_001.fastq.gz",
                    "/mnt/projects/pat/bcl_direct/t1/p1/a_S1_L007_R3_001.fastq.gz",
                    "/mnt/projects/pat/bcl_direct/t1/p1/a_S1_L008_R1_001.fastq.gz",
                    "/mnt/projects/pat/bcl_direct/t1/p1/a_S1_L008_R3_001.fastq.gz",
                    "/mnt/projects/pat/bcl_direct/t1/p1/b_S2_L001_R1_001.fastq.gz",
                    "/mnt/projects/pat/bcl_direct/t1/p1/b_S2_L001_R3_001.fastq.gz",
                    "/mnt/projects/pat/bcl_direct/t1/p1/b_S2_L002_R1_001.fastq.gz",
                    "/mnt/projects/pat/bcl_direct/t1/p1/b_S2_L002_R3_001.fastq.gz",
                    "/mnt/projects/pat/bcl_direct/t1/p1/b_S2_L003_R1_001.fastq.gz",
                    "/mnt/projects/pat/bcl_direct/t1/p1/b_S2_L003_R3_001.fastq.gz",
                    "/mnt/projects/pat/bcl_direct/t1/p1/b_S2_L004_R1_001.fastq.gz",
                    "/mnt/projects/pat/bcl_direct/t1/p1/b_S2_L004_R3_001.fastq.gz",
                    "/mnt/projects/pat/bcl_direct/t1/p1/b_S2_L005_R1_001.fastq.gz",
                    "/mnt/projects/pat/bcl_direct/t1/p1/b_S2_L005_R3_001.fastq.gz",
                    "/mnt/projects/pat/bcl_direct/t1/p1/b_S2_L006_R1_001.fastq.gz",
                    "/mnt/projects/pat/bcl_direct/t1/p1/b_S2_L006_R3_001.fastq.gz",
                    "/mnt/projects/pat/bcl_direct/t1/p1/b_S2_L007_R1_001.fastq.gz",
                    "/mnt/projects/pat/bcl_direct/t1/p1/b_S2_L007_R3_001.fastq.gz",
                    "/mnt/projects/pat/bcl_direct/t1/p1/b_S2_L008_R1_001.fastq.gz",
                    "/mnt/projects/pat/bcl_direct/t1/p1/b_S2_L008_R3_001.fastq.gz",
                    "/mnt/projects/pat/bcl_direct/t1/p1/c_S3_L001_R1_001.fastq.gz",
                    "/mnt/projects/pat/bcl_direct/t1/p1/c_S3_L001_R3_001.fastq.gz",
                    "/mnt/projects/pat/bcl_direct/t1/p1/c_S3_L002_R1_001.fastq.gz",
                    "/mnt/projects/pat/bcl_direct/t1/p1/c_S3_L002_R3_001.fastq.gz",
                    "/mnt/projects/pat/bcl_direct/t1/p1/c_S3_L003_R1_001.fastq.gz",
                    "/mnt/projects/pat/bcl_direct/t1/p1/c_S3_L003_R3_001.fastq.gz",
                    "/mnt/projects/pat/bcl_direct/t1/p1/c_S3_L004_R1_001.fastq.gz",
                    "/mnt/projects/pat/bcl_direct/t1/p1/c_S3_L004_R3_001.fastq.gz",
                    "/mnt/projects/pat/bcl_direct/t1/p1/c_S3_L005_R1_001.fastq.gz",
                    "/mnt/projects/pat/bcl_direct/t1/p1/c_S3_L005_R3_001.fastq.gz",
                    "/mnt/projects/pat/bcl_direct/t1/p1/c_S3_L006_R1_001.fastq.gz",
                    "/mnt/projects/pat/bcl_direct/t1/p1/c_S3_L006_R3_001.fastq.gz",
                    "/mnt/projects/pat/bcl_direct/t1/p1/c_S3_L007_R1_001.fastq.gz",
                    "/mnt/projects/pat/bcl_direct/t1/p1/c_S3_L007_R3_001.fastq.gz",
                    "/mnt/projects/pat/bcl_direct/t1/p1/c_S3_L008_R1_001.fastq.gz",
                    "/mnt/projects/pat/bcl_direct/t1/p1/c_S3_L008_R3_001.fastq.gz"
                ])

    def test_ilmn_bcl2fastq_v184(self):
        args = martian.Record({
            'input_mode': "ILMN_BCL2FASTQ",
            'sample_def': [
                        {
                                        "gem_group": None,
                                        "lanes": None,
                                        "read_path": "/mnt/projects/pat/bcl_direct/v184_t1/Project_proj1/Sample_a",
                                        "samples": [
                                                            "a"
                                                        ]
                                    }
                        ],
            'barcode_whitelist': "737K-april-2014",
        }
        )
        outs = martian.Record({})
        if os.path.exists("/mnt/projects/pat/bcl_direct/v184_t1/Project_proj1/Sample_a"):
            main(args,outs)
            self.assertTrue(outs.barcodes[0] == "/mnt/projects/pat/bcl_direct/v184_t1/Project_proj1/Sample_a/a_CGCTTCAA_L001_R2_001.fastq.gz")
            self.assertTrue(outs.barcodes_reverse_complement[0] == False)
            self.assertTrue(outs.reads_interleaved == False)
            self.assertTrue(outs.sample_indices[10] == None)
            self.assertTrue(outs.is_read1[0] == True and outs.is_read1[1] == False)

    def test_ilmn_bcl2fastq_miseq(self):
        args = martian.Record({
            'input_mode': "ILMN_BCL2FASTQ",
            'sample_def': [{
                            "gem_group": None,
                            "lanes": None,
                            "read_path": "/mnt/projects/pat/bcl_direct/t1_miseq/p1",
                            "samples":  [
                                                "9968"
                                        ]
                            }],
            'barcode_whitelist': "737K-april-2014",
        }
        )
        outs = martian.Record({})
        if os.path.exists("/mnt/projects/pat/bcl_direct/t1_miseq/p1"):
            main(args,outs)
            self.assertTrue(outs.barcodes == [
            "/mnt/projects/pat/bcl_direct/t1_miseq/p1/9968_S3_L001_R2_001.fastq.gz"
            ])
            self.assertTrue(outs.barcodes_reverse_complement[0] == False)
            self.assertTrue(outs.reads_interleaved == False)
            self.assertTrue(outs.sample_indices == [None])
            self.assertTrue(outs.is_read1[0] == True and outs.is_read1[1] == False)
            self.assertTrue(outs.reads == [
                        "/mnt/projects/pat/bcl_direct/t1_miseq/p1/9968_S3_L001_R1_001.fastq.gz",
                        "/mnt/projects/pat/bcl_direct/t1_miseq/p1/9968_S3_L001_R3_001.fastq.gz"
                    ])


