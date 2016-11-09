#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# Mark PCR duplicates in a BAM file
#
import random
import json
import itertools
import math
import tenkit.stats as tk_stats
import tenkit.dict_utils
import tenkit.bio_io as tk_io
import tenkit.bam as tk_bam
import tenkit.lane as tk_lane
from tenkit.read_filter import stringent_read_filter
from tenkit.constants import DUPLICATE_SUBSAMPLE_COVERAGES
import tenkit.coverage

import martian

__MRO__ = """
tage MARK_DUPLICATES(
    in  bam     input,
    in  int     perfect_read_count,
    in  bed     targets_file,
    out bam     output,
    out bam.bai index,
    out json    duplicate_summary,
    src py      "stages/reads/mark_duplicates",
) split using (
    in  float   estimated_coverage,
    in  map     lane_map,
    in  string  chunk_start,
    in  string  chunk_end,
)
"""

# For optical duplicate detection
OPTICAL_DUPLICATE_DISTANCE=100

# For diffusion duplicate detection, max distance over which diffusion is expected
MAX_DIFFUSION_DUP_DISTANCE=25e3

def chunk_bound_func(read):
    if not read.is_unmapped:
        return (read.tid, read.pos)
    else:
        return None

def split(args):
    # Chunk bam to get 1GB per chunk
    bam_in = tk_bam.create_bam_infile(args.input)
    chunk_defs = tk_bam.chunk_bam_records(bam_in, chunk_bound_func, chunk_size_gb=0.75)

    if args.targets_file is None:
        targets_path = None
    else:
        targets_path = args.targets_file

    estimated_coverage = tenkit.coverage.estimate_mean_coverage(targets_path, bam_in, lambda x: stringent_read_filter(x, False))
    for chunk in chunk_defs:
        chunk['estimated_coverage'] = estimated_coverage

    lane_coord_sys = tk_lane.LaneCoordinateSystem()

    # Reopen BAM for estimating tile extents
    bam_in = tk_bam.create_bam_infile(args.input)
    lane_coord_sys.estimate_tile_extents(bam_in)
    for chunk in chunk_defs:
        chunk['lane_map'] = lane_coord_sys.to_dict()

    return {'chunks': chunk_defs}

def main(args, outs):
    main_mark_duplicates(args, outs)


# Generator utilities -- move to tenkit?
def consumer(func):
    ''' decorator for initializing a generator consumer function '''
    def start(*args,**kwargs):
        c = func(*args,**kwargs)
        c.next()
        return c
    return start

def broadcast(source, consumers):
    ''' send each item in the source generator to each consumer in the list '''
    for item in source:
        for c in consumers:
            c.send(item)

    for c in consumers:
        c.close()


def main_mark_duplicates(args, outs):
    """
    Mark exact duplicate reads in the BAM file. Duplicates have the same read1 start site and read2 start site
    """

    lane_coord_sys = tk_lane.LaneCoordinateSystem.from_dict(args.lane_map)

    args.coerce_strings()
    outs.coerce_strings()

    bam_in = tk_bam.create_bam_infile(args.input)
    bam_out, tids = tk_bam.create_bam_outfile(outs.output, None, None, template=bam_in, pgs=tk_bam.make_pg_header(martian.get_pipelines_version(), "mark_duplicates"))

    # Determine whether the BAM has 10x barcodes
    bam_in.reset()
    has_barcodes = [tk_io.read_has_barcode(x) for x in itertools.islice(bam_in, 1000)]
    have_barcodes = (float(sum(has_barcodes)) / len(has_barcodes)) > 0.1

    # We do the subsampling to achieve the desired coverage on _perfect reads_, as
    # defined by tenkit.read_filter.stringent_read_filter.  This is tallied in ATTACH_BCS,
    # and passed into the perfect_read_count argument.  We will fail if it's not supplied.
    total_coverage = args.estimated_coverage

    # Set a fixed random seed to eliminate noise in metrics
    random.seed(0)

    sampling_rates = []
    for sample_cov in DUPLICATE_SUBSAMPLE_COVERAGES:
        rate = tk_stats.robust_divide(float(sample_cov), total_coverage)
        sampling_rates.append((rate, sample_cov))

    # All read duplicate marking - these dup decisions are written to bam_out
    # the output bam has BC aware dup marking if available.
    # Ensure the summary key indicates what kind of dup marking was actually performed.
    if have_barcodes:
        no_filter_dups_bcs =    DupSummary(False, 1.0, True,  "no_filter_full_use_bcs", lane_coord_sys, bam_out)
        no_filter_dups_no_bcs = DupSummary(False, 1.0, False, "no_filter_full_ignore_bcs", lane_coord_sys, write_to_stdout=False)
    else:
        no_filter_dups_bcs =    DupSummary(False, 1.0, True,  "no_filter_full_use_bcs", lane_coord_sys)
        no_filter_dups_no_bcs = DupSummary(False, 1.0, False, "no_filter_full_ignore_bcs", lane_coord_sys, bam_out, write_to_stdout=False)


    # Dup marking on all perfect reads
    full_dups_bcs = DupSummary(True, 1.0, True, "full_use_bcs", lane_coord_sys)
    full_dups_no_bcs = DupSummary(True, 1.0, False, "full_ignore_bcs", lane_coord_sys)

    # Make a battery of duplicate summaries at different coverages, with and w/o
    # barcode splitting
    split_options = [True, False]

    dup_sums = [full_dups_bcs, full_dups_no_bcs, no_filter_dups_bcs, no_filter_dups_no_bcs]

    # Duplicate marking on perfect reads subsampled to the requested coverage
    for (sr, cov) in sampling_rates:
        for split_bc in split_options:
            description = "cov_" + str(cov) + ('_use_bcs' if split_bc else '_ignore_bcs')
            dup_sums.append(DupSummary(True, sr, split_bc, description, lane_coord_sys))

    # Now broadcast the selected reads to the summarizers
    # We can't do the points the require a sample_rate > 1.0 so, skip those.
    # If we don't have barcodes, don't run the set that are split by barcode.
    consumers = [x.read_consumer() for x in dup_sums if x.sample_rate <= 1.0 and ((not x.split_bcs) or have_barcodes)]

    source = tk_bam.read_bam_chunk(bam_in, (args.chunk_start, args.chunk_end))
    broadcast(source, consumers)

    # We close the BAM
    bam_out.close()
    # Note - the indexing happens in join

    # Package up the summaries:
    dup_results = {}
    for x in dup_sums:
        (dups, optical_dups, diff_dups) = x.result
        desc = x.description
        dup_results[desc] = dups
        optical_desc = "optical_" + desc
        dup_results[optical_desc] = optical_dups
        diff_desc = "diffusion_" + desc
        dup_results[diff_desc] = diff_dups

    if outs.duplicate_summary:
        f = open(outs.duplicate_summary, 'w')
        json.dump(dup_results, f)
        f.close()




class DupSummary:
    def __init__(self, perfect_read_filter, sample_rate, split_bcs, description, lane_coordinate_system, output_bam=None, write_to_stdout=False):
        ''' Summarize dups at a given subsampling rate, and barcode
            splitting policy.  If an open output_bam pysam.Samfile
            is passed, dups will be marked and reads will be written
            to output_bam '''

        self.perfect_read_filter = perfect_read_filter
        self.require_barcode = perfect_read_filter and split_bcs
        self.sample_rate = sample_rate
        self.split_bcs = split_bcs
        self.description = description
        self.output_bam = output_bam
        self.write_to_stdout = write_to_stdout

        self.lane_coordinate_system = lane_coordinate_system

        # default return value -- will be used if the sample_rate > 1.0
        self.result = (None, None, None)


    def count_dups_by_distance(self, reads):
        '''Count number of nearby duplicates in a set of reads.  A pair is counted as 1'''
        # Get (flowcell, lane, surface, swath, tile, x, y) tuples for each read
        read_locs = []
        for (key, read, idx) in reads:
            read_loc = tk_lane.extract_read_position(read)
            if read_loc is not None:
                read_locs.append((read_loc, read))

        # Sort by flowcell_lane
        def flowcell_lane(read_loc):
            return "%s_%s" % (read_loc[0].flowcell, read_loc[0].lane)
        read_locs.sort(key=flowcell_lane)
        lane_groups = itertools.groupby(read_locs, flowcell_lane)

        opt_dups_found = 0   # really close dupes
        diff_dups_found = 0  # somewhat close dupes

        # Measure distances between all pairs in a lane
        for (lane, lane_reads) in lane_groups:
            lane_reads = list(lane_reads)

            layout = self.lane_coordinate_system.get_layout_for_read_loc(lane_reads[0][0])
            test_dups = layout.has_diffusion_duplicates(MAX_DIFFUSION_DUP_DISTANCE)

            if len(lane_reads) > 100:
                martian.log_info("Got dup cluster of size: %d" % len(lane_reads))
                first_read = lane_reads[0][1]
                martian.log_info("tid: %d, pos: %d, mapq: %d, seq: %s" % (first_read.tid, first_read.pos, first_read.mapq, first_read.seq))

            opt_dups = set()
            diff_dups = set()
            dump = []
            cmp_reads = min(100, len(lane_reads))
            for i in range(cmp_reads):
                loc1, read1 = lane_reads[i]
                lane_loc1 = self.lane_coordinate_system.convert_to_lane_coords(loc1)

                for j in range(i+1, len(lane_reads)):
                    loc2, read2 = lane_reads[j]
                    lane_loc2 = self.lane_coordinate_system.convert_to_lane_coords(loc2)

                    dist = math.sqrt((lane_loc1[0]-lane_loc2[0])**2 + (lane_loc1[1] - lane_loc2[1])**2)
                    if test_dups and dist < MAX_DIFFUSION_DUP_DISTANCE:
                        diff_dups.add(j)
                        if self.write_to_stdout and j not in diff_dups:
                            dump.append(("%d\t"+("%d\t"*14)) % (dist,
                                                                loc1.surface, loc1.swath, loc1.tile, loc1.x, loc1.y, lane_loc1[0], lane_loc1[1],
                                                                loc2.surface, loc2.swath, loc2.tile, loc2.x, loc2.y, lane_loc2[0], lane_loc2[1]))

                    if dist < OPTICAL_DUPLICATE_DISTANCE:
                        opt_dups.add(j)

            if self.write_to_stdout and len(diff_dups) >= 2:
                for x in dump:
                    print ("%d\t%s" % (len(diff_dups), x))

            diff_dups_found += len(diff_dups)
            opt_dups_found += len(opt_dups)

        return (opt_dups_found, diff_dups_found)


    @consumer
    def read_consumer(self):

        def process_read_block(reads, count_hist, optical_dup_hist, diffusion_dup_hist):
            ''' dedup a block of reads, then write to BAM in original order '''
            read_tuples = []
            for read in reads:
                # Don't consider unmapped reads
                if read.is_unmapped:
                    continue

                bc_sequence = None

                if self.split_bcs:
                    bc_sequence = tk_io.get_read_barcode(read)

                read_tup = (bc_sequence, read.is_read1, read.is_reverse, read.tid, read.pos, read.mrnm, read.mpos)

                # Include the read and the original index so that we can get back to the original order
                read_tuples.append((read_tup, read, len(read_tuples)))

            # Sort by the read tuple -- need to do this for the groupby to get all the common items together
            read_tuples.sort(key=lambda x: x[0])

            # Both reads in a pair must be mapped, otherwise we drop the pair
            mapped_pair_tuples = itertools.ifilter(lambda x: (not x[1].is_unmapped) and (not x[1].mate_is_unmapped), read_tuples)

            # Group the reads by the read_tuple
            dup_groups = itertools.groupby(mapped_pair_tuples, lambda x: x[0])

            for (key, dup_group) in dup_groups:
                # Note how many dups we have
                dup_group = list(dup_group)
                n_dups = len(dup_group)

                if n_dups > 1:
                    optical_dups, diffusion_dups = self.count_dups_by_distance(dup_group)
                else:
                    optical_dups = 0
                    diffusion_dups = 0

                # diffusion dups encompass optical dups, if
                non_proximal_dups = n_dups - max(diffusion_dups, optical_dups)

                # If we are splitting on bcs, then only counts stats for read groups with BCs
                group_bc = key[0]
                if not self.split_bcs or group_bc is not None:
                    count_hist[non_proximal_dups] = count_hist.setdefault(non_proximal_dups,0) + 1
                    if optical_dups > 0:
                        optical_dup_hist['count'] = optical_dup_hist.setdefault('count',0) + optical_dups
                        diffusion_dup_hist['count'] = diffusion_dup_hist.setdefault('count',0) + diffusion_dups

                # Mark dups
                if self.output_bam:
                    # mark the first read in the set as not-a-dup
                    dup_group[0][1].is_duplicate = False

                    for i in range(1, n_dups):
                        dup_group[i][1].is_duplicate = True

            if self.output_bam:
                for read in reads:
                    self.output_bam.write(read)

            # done processing block of reads
            return


        # bam is only sorted by chrom and pos.  We pull in all the reads with common chrom and pos, and dedup
        # them as a group
        current_bam_key = (-1, -1)
        current_reads = []

        # we will track the histogram of the reads/molecule
        molecule_count_hist = {}
        optical_dup_hist = {'count':0}
        diffusion_dup_hist = {'count':0}

        try:
            while True:
                # accept the next read
                read = (yield)

                # Sample to the requested subsample rate
                if self.sample_rate < 1.0:
                    if random.random() > self.sample_rate:
                        continue

                # Apply the perfect read filter - the perfect read filter requires barcodes if they are
                # being split on
                if self.perfect_read_filter and not stringent_read_filter(read, self.require_barcode):
                    continue

                new_bam_key = (read.tid, read.pos)

                # If the dup group gets extremely large we can run out of memory.
                # Process things in groups of 500K to prevent memory blow-up
                # May cause us to miss a few dups, but it doesn't really matter in these crazy regions
                if new_bam_key != current_bam_key or len(current_reads) > 500000:
                    # Otherwise process our block of reads
                    process_reads = current_reads
                    current_reads = [read]
                    current_bam_key = new_bam_key

                    if len(process_reads) > 0:
                        process_read_block(process_reads, molecule_count_hist, optical_dup_hist, diffusion_dup_hist)
                else:
                    # accumulate block of reads with same start position
                    current_reads.append(read)

                # Short circuit the bolus of reads
                # that have no mapping info --
                # they pass straight through
                if current_bam_key == (-1, -1) or current_bam_key == (-1, 0):
                    current_reads = []
                    if self.output_bam:
                        self.output_bam.write(read)


        except GeneratorExit:
            # Finish up final batch
            process_read_block(current_reads, molecule_count_hist, optical_dup_hist, diffusion_dup_hist)

            # save the result for later use
            self.result = (molecule_count_hist, optical_dup_hist, diffusion_dup_hist)
            return

def join(args, outs, chunk_defs, chunk_outs):
    outs.coerce_strings()

    # combine the duplicate summary counts
    dup_summaries = [json.load(open(out.duplicate_summary)) for out in chunk_outs]
    combined_dups = reduce(lambda x,y: tenkit.dict_utils.add_dicts(x,y,2), dup_summaries)
    combined_dups['read_counts'] = {}
    combined_dups['read_counts']['perfect_read_count'] = args.perfect_read_count

    f = open(outs.duplicate_summary, 'w')
    json.dump(combined_dups, f)
    f.close()

    # combine & index the chunks of the BAM
    tk_bam.concatenate(outs.output, [c.output for c in chunk_outs])
    tk_bam.index(outs.output)
    outs.index = outs.output + '.bai'
