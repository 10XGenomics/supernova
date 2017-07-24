// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

#ifndef TENX_CLEAN_THE_H
#define TENX_CLEAN_THE_H

#include "CoreTools.h"
#include "feudal/PQVec.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "10X/paths/ReadPathVecX.h"

void BuildLineGraph( const digraphE<vec<int>>& D, const vec<int>& dinv,
     const vec<vec<vec<vec<int>>>>& dlines, digraphE<int>& D2 );

void SnipFlipSquares( digraphE<vec<int>>& D, vec<int>& dinv,
     const ReadPathVec& dpaths, const Bool verbose = False );

void BarcodeJoin( const HyperBasevectorX& hb, const vec<int>& inv,
     digraphE<vec<int>>& D, vec<int>& dinv,
     ReadPathVec& dpaths, const vec<Bool>& dup, const ReadPathVecX& pathsx, 
     const vec<Bool>& pmask, const vec<int64_t>& bci, 
     const vec< triple<int,int,int> >& qept, const VecIntVec& ebcx,
     const Bool break_trinads = False );

void BigPullInv( const HyperBasevectorX& hb, const vec<int>& inv,
     digraphE<vec<int>>& D, vec<int>& dinv, const int MIN_LEN );

void DiploidFilter( const HyperBasevectorX& hb, digraphE<vec<int>>& D, 
     vec<int>& dinv );

void RemoveVerySmallComponents( const HyperBasevectorX& hb, 
     digraphE<vec<int>>& D, vec<int>& dinv, const int MIN_COMP );

void DeleteSomeRedundantGaps(
     const HyperBasevectorX& hb, digraphE<vec<int>>& D, vec<int>& dinv, Bool verbose = False );

void HangBeGone( const HyperBasevectorX& hb, digraphE<vec<int>>& D, vec<int>& dinv,
     const int MIN_RATIO2 = 20, const int MAX_KILL2 = 700, 
     const int verbosity = 1, const Bool NEW = False );

void RemoveRedundantCells( digraphE<vec<int>>& D, vec<int>& dinv );

void DeleteParallelGaps( digraphE<vec<int>>& D, vec<int>& dinv, 
     Bool verbose = False );

void ZapPairGaps( digraphE<vec<int>>& D, vec<int>& dinv,
     Bool verbose=False );

void KillPairGapsAtBranches( digraphE<vec<int>>& D, vec<int>& dinv, 
     Bool verbose = False );

void DeletePairGapsInCells( digraphE<vec<int>>& D, vec<int>& dinv,
     vec<vec<vec<vec<int>>>> & dlines, vec<int> & dels, Bool verbose = False );

// Are we around a cell or sequence gap?

Bool AroundGap( const digraphE<vec<int>> & D, const vec<int> & to_left,
     const vec<int> & to_right, const int & d );

// Wherever there is a branch, with support >=10 to 0, liberate the weak
// branch, so long as it contains some unique kmers.
     
void KillZeroSupportBranches( const HyperBasevectorX & hb, const vec<int> & inv,
     digraphE<vec<int>> & D, vec<int> & dinv, const ReadPathVec & dpaths, 
     const VecULongVec & dpaths_index, vec<Bool> & to_delete, const Bool verbose );

// Delete zero support edges of a graph based on the following criterion:
// if an edge between vertices v and w has zero weight, AND
// there is an alternate path with support between v and w, DELETE.
// support: < 0 for edges that shouldn't be included in alt paths
//        : = 0 for edges that are considered for deletion
//        : > 0 for edges that are to count towards alt paths
// edges that are to be deleted are marked so in to_del

template <class F>
void DeleteZeroSupportEdges( const digraphE<F> & G, const vec<int> & ginv,
     const vec<int> & to_left, const vec<int> & to_right,
     const vec<int> & support, vec<Bool> & to_del );

// If there is a pair gap between two vertices, and there is a path
// through the graph that contains at least one regular edge
// then delete the pair gap

void PairGapKiller( digraphE<vec<int>> & D, vec<int> & dinv, const int MAX_DEPTH );

void MarkEdgesWithZeroReads(digraphE<vec<int>> & D, vec<int> & dinv,
     const vec<int> to_left, const vec<int> to_right,
     const VecULongVec & dpaths_index, vec<Bool> & to_delete );

// Delete zero support bridges
void MarkWeakBridges( const digraphE<vec<int>> & D, const vec<int> & dinv,
     const vec<int> & to_left, const vec<int> & to_right,
     const VecULongVec & dpaths_index, vec<Bool> & to_delete );

// Delete zero support lines that are hanging
void MarkDuplicateHangingLines( const HyperBasevectorX & hb, 
     const digraphE<vec<int>> & D, const vec<int> & dinv,
     const vec<vec<vec<vec<int>>>> & dlines, const VecULongVec & dpaths_index,
     vec<Bool> & to_delete );

// If a bc gap is in the canonical topology between two edges
// and there is sufficient overlap between the flanks then
// just glue them together

void DeleteRedundantBCGaps( const HyperBasevectorX & hb, 
     digraphE<vec<int>> & D, vec<int> & dinv );

// Mark inversion artifacts for deletion

void MarkInversionArtifacts( const HyperBasevectorX & hb, const vec<int> & inv,
     digraphE<vec<int>> & D, vec<int> & dinv, const vec<vec<vec<vec<int>>>> & dlines,
     const VecULongVec & dpaths_index, vec<Bool> & to_delete,
     const Bool verbose = False );

// Compute the number of read-pair inserts that go across pairs of lines
// stored in a map<pair<int,int>,int> data structure
// (line1, line2) -> # of read pair inserts

void GetInsertsConnectingLines( const vec<int> & dinv, const vec<int> & linv,
     const ReadPathVec & dpaths, const VecULongVec & dpaths_index, const vec<int> & dtol,
     map<pair<int,int>,int> & lsupp );

// Just a validation of D, dinv and dpaths

void SanityCheck( const digraphE<vec<int>> & D, const vec<int> & dinv,
     const ReadPathVec & dpaths );

// Pull apart 2-in/1-middle/2-out topologies using read pairs

void ReadPairPullApart( digraphE<vec<int>> & D, vec<int> & dinv,
     const ReadPathVec & dpaths, const VecULongVec & dpaths_index );

// If there are sufficiently many read-pairs between two lines,
// and those two lines don't have a similar relationship with
// other lines, then connect them with a pair gap

void HookLinesWithReadPairs( const HyperBasevectorX & hb, 
     digraphE<vec<int>> & D, vec<int> & dinv, const vec<vec<vec<vec<int>>>> & dlines,
     const ReadPathVec & dpaths, VecULongVec & dpaths_index );

// If there is an edge with zero support then delete it if there is an 
// alternate path that contains at least one edge with > 0 reads.
// This is a little dangerous, simplifies but potentially deletes true
// sequence

void DeleteZeroSupportEdgesWithAltPath( digraphE<vec<int>> & D, vec<int> & dinv,
     vec<Bool> & to_delete, const VecULongVec & dpaths_index );

// Master clean up function

void CleanTheAssembly( const HyperBasevectorX& hb, const vec<int>& inv,
     digraphE<vec<int>>& D, vec<int>& dinv, ReadPathVec& dpaths,
     const vec<Bool>& dup,
     const ReadPathVecX& pathsx, const vec<Bool>& pmask, const vec<int64_t>& bci,
     const vec< triple<int,int,int> >& qept, const VecIntVec& ebcx,
     const String& udir, const String& suffix,
     const Bool intermediates, const String INSTANCE, const int substart, 
     const int substop, const Bool LOWMEMB, const Bool LOOPS2_VERBOSE,
     const Bool ZIPPER, const Bool MORE );

#endif
