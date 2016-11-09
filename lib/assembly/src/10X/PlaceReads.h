// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

#ifndef TENX_PLACE_READS_H
#define TENX_PLACE_READS_H

#include "CoreTools.h"
#include "math/HoInterval.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "10X/paths/ReadPathVecX.h"

// Align a path to the supergraph.

void Align(
     const digraphE<vec<int>>& D, // supergraph
     const vec<int>& to_left,     // edge to left
     const vec<int>& to_right,    // edge to right
     const vec<int>& x,           // input path
     // input/output:
     ho_interval& xpos,           // aligned start/stop positions on x
     vec<int>& d,                 // path in D
     int& dstart,                 // index of start edge on first edge in d
     int& dstop                   // index+1 of stop edge on last edge in d
          );

void Align2(
     const digraphE<vec<int>>& D, // supergraph
     const vec<int>& to_left,     // edge to left
     const vec<int>& to_right,    // edge to right
     const vec<int>& x,           // input path
     // input/output:
     ho_interval& xpos,           // aligned start/stop positions on x
     vec<int>& d,                 // path in D
     int& dstart,                 // index of start edge on first edge in d
     int& dstop                   // index+1 of stop edge on last edge in d
          );

// PlaceReads: only find unique placements.

void PlaceReads( const HyperBasevectorX& hb, const ReadPathVec& paths, 
     const vec<Bool>& dup, const digraphE<vec<int>>& D, ReadPathVec& dpaths,
     const Bool verbose, const Bool single, const Bool align2 = False );

// PlaceReadsSmart.  Use barcode localization to enhance output of PlaceReads.
// There are some limitations of the current code:
//
// 1. We require, in effect, three het sites to localize a barcode to a phased 
//    region.  Three is probably too high a bar.
//
// 2. There will be barcodes having parallel but ambiguous placements on the two
//    arms of a megabubble.  It would probably be useful to rescue these placements,
//    but tag them as having two homes.
//
// 3. One sees occasional instances where a barcode has been placed twice on a line,
//    with separation roughly in the range of 100-200 kb, presumably representing
//    very long molecules having low coverage.  Probably the reads in the middle
//    could be rescued.

void PlaceReadsSmart( const HyperBasevectorX& hb, const ReadPathVec& paths,
     const vec<Bool>& dup, const digraphE<vec<int>>& D, const vec<int>& dinv,
     ReadPathVec& dpaths, const vec<vec<vec<vec<int>>>>& dlines, 
     const vec<int64_t>& bci, const Bool verbose, const int btest = -1,
     const Bool align2 = False );

void PlaceReads( const HyperBasevectorX& hb, const ReadPathVecX& paths, 
     const vec<Bool>& dup, const digraphE<vec<int>>& D, ReadPathVec& dpaths,
     const Bool verbose, const Bool single, const Bool align2 = False );

void PlaceReadsSmart( const HyperBasevectorX& hb, const ReadPathVecX& paths,
     const vec<Bool>& dup, const digraphE<vec<int>>& D, const vec<int>& dinv,
     ReadPathVec& dpaths, const vec<vec<vec<vec<int>>>>& dlines, 
     const vec<int64_t>& bci, const Bool verbose, const int btest = -1,
     const Bool align2 = False );

#endif
