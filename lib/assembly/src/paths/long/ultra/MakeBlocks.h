///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeBlocks. This takes as input K, some reads, and some kmer aligns of the 
// reads with id >= 1 to read 0.  It produces as output "blocks", sequences that
// are supposed to be true sequences for the genome underlying the read.  These may
// overlap by up to K-1.
// (Now produces a system of threaded blocks.)

#ifndef MAKE_BLOCKS_H
#define MAKE_BLOCKS_H

#include "Basevector.h"
#include "CoreTools.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/ultra/ThreadedBlocks.h"

template<int K> void MakeBlocks( 
     const int ec,                            // read id
     const vecbasevector& reads,
     const vec< vec< pair<int,int> > >& a,    // kmer aligns to read 0
     threaded_blocks& tb,                     // output
     const basevector& gread0,                // true genomic sequence for read 0
                                              // and flanking region
     const vec<basevector>& gkmers,
     const int gid,
     const int gstart,
     ostream& out,                            // for logging
     const Bool FILTER_BAD_BLOCKS,
     const long_logging_control& log_control,
     const long_logging& logc );

#endif
