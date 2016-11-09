///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef PREFAB_H
#define PREFAB_H

#include "CoreTools.h"
#include "paths/long/CreateGenome.h"
#include "paths/long/Friends.h"
#include "paths/long/LongProtoTools.h"

// Select a set of reads for initial processing.  Currently this returns all
// reads.

// SelectInitialReads.  Select a set of reads that should provide 'adequate' 
// coverage of the genome, for initial processing.
// 
// input:  rid = totality of reads to be processed
// output: rid1 = subset of rid1, to be processed initially
//
// Initial implementation: returns a random 20% subset of the reads.

void SelectInitialReads( const vec<int>& rid, const IAndOsVec& F, vec<int>& rid1,
     const long_heuristics& heur, const long_logging_control& log_control );

// CorrectSomeMoreReads.  After initial reads have been processed, correct some
// more of them by somehow using the results of the initial processing.
//
// input:  rid = totality of reads to be processed
//         rid1 = subset of rid1, processed initially
//         (corrected, cid) = results from processing of rid1
// output: (corrected2, cid2) = results of this step
//
// Initial implementation: does nothing!!

void CorrectSomeMoreReads( 
     const vecbasevector& reads,            // reads
     const vec<int>& rid, const vec<int>& rid1, const vec< vec<int> >& friend_ids1,
     const IAndOsVec& F, const VecEFasta& corrected, const vec<int>& cid,
     VecEFasta& corrected2, vec<int>& cid2, const long_heuristics& heur,
     const long_logging_control& log_control, const long_logging& logc,
     const ref_data& ref );

#endif
