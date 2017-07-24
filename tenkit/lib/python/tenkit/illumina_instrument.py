#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#
# Code to detect illumina instrument from read names

import gzip
import re
import martian

# dictionary of instrument id regex: [platform(s)]
InstrumentIDs = {"HWI-M[0-9]{4}$" : ["MiSeq"],
        "HWUSI" : ["Genome Analyzer IIx"],
        "M[0-9]{5}$" : ["MiSeq"],
        "HWI-C[0-9]{5}$" : ["HiSeq 1500"],
        "C[0-9]{5}$" : ["HiSeq 1500"],
        "HWI-D[0-9]{5}$" : ["HiSeq 2500"],
        "D[0-9]{5}$" : ["HiSeq 2500"],
        "J[0-9]{5}$" : ["HiSeq 3000"],
        "K[0-9]{5}$" : ["HiSeq 3000","HiSeq 4000"],
        "E[0-9]{5}$" : ["HiSeq X"],
        "NB[0-9]{6}$": ["NextSeq"],
        "NS[0-9]{6}$" : ["NextSeq"],
        "MN[0-9]{5}$" : ["MiniSeq"]}

# dictionary of flow cell id regex: ([platform(s)], flow cell version and yeild)
FCIDs = {"C[A-Z,0-9]{4}ANXX$" : (["HiSeq 1500", "HiSeq 2000", "HiSeq 2500"], "High Output (8-lane) v4 flow cell"),
         "C[A-Z,0-9]{4}ACXX$" : (["HiSeq 1000", "HiSeq 1500", "HiSeq 2000", "HiSeq 2500"], "High Output (8-lane) v3 flow cell"),
         "H[A-Z,0-9]{4}ADXX$" : (["HiSeq 1500", "HiSeq 2500"], "Rapid Run (2-lane) v1 flow cell"),
         "H[A-Z,0-9]{4}BCXX$" : (["HiSeq 1500", "HiSeq 2500"], "Rapid Run (2-lane) v2 flow cell"),
         "H[A-Z,0-9]{4}BCXY$" : (["HiSeq 1500", "HiSeq 2500"], "Rapid Run (2-lane) v2 flow cell"),
         "H[A-Z,0-9]{4}BBXX$" : (["HiSeq 4000"], "(8-lane) v1 flow cell"),
         "H[A-Z,0-9]{4}BBXY$" : (["HiSeq 4000"], "(8-lane) v1 flow cell"),
         "H[A-Z,0-9]{4}CCXX$" : (["HiSeq X"], "(8-lane) flow cell"),
         "H[A-Z,0-9]{4}CCXY$" : (["HiSeq X"], "(8-lane) flow cell"),
         "H[A-Z,0-9]{4}ALXX$" : (["HiSeq X"], "(8-lane) flow cell"),
         "H[A-Z,0-9]{4}BGXX$" : (["NextSeq"], "High output flow cell"),
         "H[A-Z,0-9]{4}BGXY$" : (["NextSeq"], "High output flow cell"),
         "H[A-Z,0-9]{4}BGX2$" : (["NextSeq"], "High output flow cell"),
         "H[A-Z,0-9]{4}AFXX$" : (["NextSeq"], "Mid output flow cell"),
         "A[A-Z,0-9]{4}$" : (["MiSeq"], "MiSeq flow cell"),
         "B[A-Z,0-9]{4}$" : (["MiSeq"], "MiSeq flow cell"),
         "D[A-Z,0-9]{4}$" : (["MiSeq"], "MiSeq nano flow cell"),
         "G[A-Z,0-9]{4}$" : (["MiSeq"], "MiSeq micro flow cell"),
         "H[A-Z,0-9]{4}DMXX$" : (["NovaSeq"], "S2 flow cell")}


SUPERNOVA_PLATFORM_BLACKLIST = ["HiSeq 3000","HiSeq 4000","HiSeq 3000/4000"]

_upgrade_set1 = set(["HiSeq 2000", "HiSeq 2500"])
_upgrade_set2 = set(["HiSeq 1500", "HiSeq 2500"])
_upgrade_set3 = set(["HiSeq 3000", "HiSeq 4000"])
_upgrade_set4 = set(["HiSeq 1000", "HiSeq 1500"])
_upgrade_set5 = set(["HiSeq 1000", "HiSeq 2000"])

fail_msg = "Cannot determine sequencing platform"
success_msg_template = "(likelihood: {})"
null_template = "{}"

# do intersection of lists
def intersect(a, b):
    return list(set(a) & set(b))

def union(a, b):
    return list(set(a) | set(b))

# extract ids from reads
def parse_readhead(head):
    fields = head.strip("\n").split(":")

    # if ill-formatted/modified non-standard header, return cry-face
    if len(fields) < 3:
        return -1,-1
    iid = fields[0][1:]
    fcid = fields[2]
    return iid, fcid

# infer sequencer from ids from single fastq
def infer_sequencer(iid, fcid):
    seq_by_iid = []
    for key in InstrumentIDs:
        if re.search(key,iid): 
            seq_by_iid += InstrumentIDs[key]

    seq_by_fcid = []
    for key in FCIDs:
        if re.search(key,fcid):
            seq_by_fcid += FCIDs[key][0]

    sequencers = []
   
    # if both empty
    if not seq_by_iid and not seq_by_fcid:
        return sequencers, "fail"
    
    # if one non-empty
    if not seq_by_iid:
        return seq_by_fcid, "likely"
    if not seq_by_fcid:
        return seq_by_iid, "likely"

    # if neither empty
    sequencers = intersect(seq_by_iid, seq_by_fcid)
    if sequencers:
        return sequencers, "high"
    # this should not happen, but if both ids indicate different sequencers..
    else: 
        sequencers = union(seq_by_iid, seq_by_fcid)
        return sequencers, "uncertain"

# process the flag and detected sequencer(s) for single fastq
def infer_sequencer_with_message(iid, fcid):
    sequencers, flag = infer_sequencer(iid, fcid)
    if not sequencers:
        return [""], fail_msg 

    if flag == "high":
        msg_template = null_template
    else:
        msg_template = success_msg_template

    if set(sequencers) <= _upgrade_set1:
        return ["HiSeq2000/2500"], msg_template.format(flag)
    if set(sequencers) <= _upgrade_set2:
        return ["HiSeq1500/2500"], msg_template.format(flag)
    if set(sequencers) <= _upgrade_set3:
        return ["HiSeq3000/4000"], msg_template.format(flag)
    return sequencers, msg_template.format(flag)

def test_sequencer_detection():
    Samples = ["@ST-E00314:132:HLCJTCCXX:6:2206:31213:47966 1:N:0", 
               "@D00209:258:CACDKANXX:6:2216:1260:1978 1:N:0:CGCAGTT",
               "@D00209:258:CACDKANXX:6:2216:1586:1970 1:N:0:GAGCAAG"]

    seqrs = set()
    for head in Samples:
        iid, fcid = parse_readhead(head) 
        seqr, msg = infer_sequencer_with_message(iid,fcid)
        for sr in seqr:
            signal = (sr,msg)
        seqrs.add(signal)

    print seqrs

def sequencer_detection_message(fastq_files):
    seqrs = set()
    # accumulate (sequencer, status) set
    for fastq in fastq_files:
        with gzip.open(fastq) as f:
            head = ""
            line = f.readline()
            if len(line)>0:
                if line[0]=="@":
                    head = line
                else:
                    martian.exit("Incorrectly formatted first read in FASTQ file: %s" % fastq)

        iid, fcid = parse_readhead(head)
        seqr, msg = infer_sequencer_with_message(iid,fcid)
        for sr in seqr:
            signal = (sr,msg)
        seqrs.add(signal)

    # get a list of sequencing platforms
    platforms = set()
    for platform, _ in seqrs:
        platforms.add(platform)
    sequencers = list(platforms)

    # if no sequencer detected at all
    message = ""
    fails = 0
    for platform, status in seqrs:
        if status == fail_msg:
            fails += 1
    if fails == len(seqrs):
        message = "could not detect the sequencing platform(s) used to generate the input FASTQ files"
        return message, sequencers

    # if partial or no detection failures
    if fails > 0:
        message = "could not detect the sequencing platform used to generate some of the input FASTQ files, "
    message += "detected the following sequencing platforms- "
    for platform, status in seqrs:
        if status != fail_msg:
            message += platform + " " + status + ", "
    message = message.strip(", ")
    return message, sequencers

