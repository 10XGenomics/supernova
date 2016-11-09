// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

#ifndef TENX_SUPER_H
#define TENX_SUPER_H

#include "CoreTools.h"
#include "graph/Digraph.h"
#include "math/HoInterval.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "10X/IntIndex.h"

#include "10X/paths/ReadPathVecX.h"

void KillMisassembledCells( const HyperBasevectorX& hb, const vec<Bool>& dup,
     const vec<int64_t>& bci, const ReadPathVecX& paths,
     digraphE<vec<int>>& D, vec<int>& dinv, ReadPathVec& dpaths,
     IntIndex& dpaths_index, vec<int32_t>& bc,
     vec<vec<vec<vec<int>>>>& dlines, vec<int>& dels, vec <double> &fhist,
     vec <pair<float, int>> & lr,
     const int BC_REQUIRE, int BC_FLANK, int BC_IGNORE, const Bool verbose );

void FixMisassemblies( const HyperBasevectorX& hb, const vec<Bool>& dup,
     const vec<int64_t>& bci, const ReadPathVecX& paths,
     digraphE<vec<int>>& D, vec<int>& dinv, ReadPathVec& dpaths );

void FlattenSomeBubbles( const HyperBasevectorX& hb, const vec<Bool>& dup,
     const ReadPathVecX& paths, digraphE<vec<int>>& D, vec<int>& dinv, 
     ReadPathVec& dpaths, const int MAX_DELTA, const double MIN_RATIO,
     const int MAX_DEL );

void Cleaner( const HyperBasevectorX& hb, const vec<int>& inv, 
     const ReadPathVecX& paths, const vec<Bool>& dup, digraphE<vec<int>>& D, 
     vec<int>& dinv, ReadPathVec& dpaths, const Bool verbose );

void BucketLines( const vec<vec<vec<vec<int>>>>& dlines, const vec<int>& llens,
     vec<vec<int>>& buckets, const int min_len = 0 );

void ZapInversionBubbles( const digraphE<vec<int>>& D, const vec<int>& dinv,
     vec<int>& dels );

void ZapMegaInversionBubbles( digraphE<vec<int>>& D, const vec<int>& dinv );

void KillMisassembledCellsAlt( const HyperBasevectorX& hb, digraphE<vec<int>>& D,
     vec<int>& dinv, const VecIntVec& ebcx, vec<vec<vec<vec<int>>>>& dlines,
     vec<int>& dels );

void DelWeak3( digraphE< vec<int> >& D, vec<int>& dinv, const ReadPathVec& dpaths,
     const IntIndex& dpaths_index, vec<int>& dels, 
     const int MAX_WEAK3_LOSE, const int MIN_WEAK3_WIN, const int MIN_WEAK3_RATIO,
     const Bool verbose );

void DelWeak4( digraphE< vec<int> >& D, vec<int>& dinv, const ReadPathVec& dpaths,
     const IntIndex& dpaths_index, vec<int>& dels, 
     const int MAX_WEAK4_LOSE, const int MAX_WEAK4_LOSE_TOTAL, 
     const int MIN_WEAK4_WIN, const int MIN_WEAK4_RATIO, const Bool verbose );

void DelWeak5( digraphE< vec<int> >& D, vec<int>& dinv, const ReadPathVec& dpaths,
     const IntIndex& dpaths_index, vec<int>& dels, 
     const int MAX_WEAK4_LOSE, const int MAX_WEAK4_LOSE_TOTAL, 
     const int MIN_WEAK4_WIN, const int MIN_WEAK4_RATIO, const Bool verbose );

void ComputeMult( const HyperBasevectorX& hb, const digraphE<vec<int>>& D,
     vec<int>& mult );

void LineCN( const vec<int>& kmers, const MasterVec<SerfVec<pair<int,int>>>& lbp,
     const digraphE<vec<int>>& D, const vec<vec<vec<vec<int>>>>& dlines,
     const vec<int>& llens, vec<double>& COV );

// Splay vertices at the ends of long lines.

void Splay( digraphE<vec<int>>& D, const vec<int>& dinv,
     const vec<vec<vec<vec<int>>>>& dlines, const vec<int>& llens,
     const int MIN_SPLAY );

void GetLineLengths( const HyperBasevectorX& hb, const digraphE<vec<int>>& D,
     const vec<vec<vec<vec<int>>>>& dlines, vec<int>& llens );
void GetLineLengths( const vec<int>& kmers, const digraphE<vec<int>>& D,
     const vec<vec<vec<vec<int>>>>& dlines, vec<int>& llens );
void GetLineLengths( const vec<int>& dlens,
     const vec<vec<vec<vec<int>>>>& dlines, vec<int>& llens );

void FindCompoundHangs( const digraphE<vec<int>>& D, const vec<int>& dinv,
     const vec<int>& lens, const vec<int>& dfw, vec<int>& dels,
     const int MAX_TINY, const double MIN_RATIO, const Bool verbose );

void PullApart( const vec<int>& inv, digraphE<vec<int>>& D, vec<int>& dinv, 
     const ReadPathVec& dpaths, vec<int>& dels, Bool verbose );

// Zipper: look for cases where two non-gap superedges emanate from the same vertex
// and agree on their first base edge.  In most cases, but not all, it is possible
// to "zipper up", eliminating the common base edge.

void Zipper( digraphE<vec<int>>& D, vec<int>& dinv );

void KillLowUnique( const HyperBasevectorX& hb, digraphE<vec<int>>& D,
     vec<int>& dels, const Bool verbose );

void KillLowUniqueFrac( const HyperBasevectorX& hb, digraphE<vec<int>>& D,
     vec<int>& dels, const double MIN_UNIQ_FRAC = 0.1, const Bool verbose = True );

void RemoveDuff( const HyperBasevectorX& hb, digraphE<vec<int>>& D,
     const vec<int>& dinv, vec<int>& dels, const Bool verbose );

void RemoveDuff2( const HyperBasevectorX& hb, digraphE<vec<int>>& D,
     vec<int>& dinv );

void SimpleHangs( const HyperBasevectorX& hb, digraphE<vec<int>>& D,
     const vec<int>& dinv, vec<int>& dels, const int MAX_KILL,
     const double MIN_RATIO, const Bool verbose, const Bool single );

void Cleaner( const HyperBasevectorX& hb, const vec<int>& inv, 
     const ReadPathVec& paths, const vec<Bool>& dup, digraphE<vec<int>>& D, 
     vec<int>& dinv, ReadPathVec& dpaths, const Bool verbose );

void CleanupCore( digraphE<vec<int>>& D, vec<int>& dinv );

// Validate: validate a supergraph assembly.  Hypothetically we could require that
// the assembly is zippered up, but only if we exclude those cases that can't be
// zippered, as defined in the Zipper code.

void Validate( const HyperBasevectorX& hb, const vec<int>& inv, 
     const digraphE<vec<int>>& D, const vec<int>& dinv );

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

void Emanate( digraphE<vec<int>>& D, vec<int>& dinv, const Bool verbose );
void Emanate2( const vec<int>& inv, digraphE<vec<int>>& D, vec<int>& dinv, 
     const Bool verbose );

void RemoveUnneededVertices( digraphE<vec<int>>& D, vec<int>& dinv );

void LineProx( const HyperBasevectorX& hb, const vec<int>& inv, 
     const VecIntVec& ebcx, const digraphE<vec<int>>& D, const vec<int>& dinv,
     const vec<vec<vec<vec<int>>>>& dlines,
     const vec< triple<int,int,int> >& qept, vec< vec< pair<int,int> > >& lhood );

int LineN50( const HyperBasevectorX& hb, const digraphE<vec<int>>& D,
     const vec<vec<vec<vec<int>>>>& dlines );

void KillInversionArtifacts( const digraphE<vec<int>>& D, const vec<int>& dinv,
     const ReadPathVec& dpaths, const IntIndex& dpaths_index, 
     const vec<int64_t>& bid, vec<int>& dels, const int MAX_CAN_INS_DEL );

int64_t CheckSum( const HyperBasevectorX& hb, const vec<int>& inv,
     const digraphE<vec<int>>& D, const vec<int>& dinv );

// MakeMerges: M = { ( (d1,p1), (d2,p2), n ) } where 
// interval [p1,p1+n) on superedge d1 matches interval [p2,p2+n) on superedge d2.
// It looks like merge set should be symmetric with respect to the involution,
// but not checked carefully.

int MakeMerges( const vec< triple< pair<int,int>, pair<int,int>, int > >& M,
     digraphE<vec<int>>& D, vec<int>& dinv, const Bool verbose = True );

// Try to merge along relatively short overlaps that are proximate in the graph.

void MergeShortOverlaps( HyperBasevectorX& hb, vec<int>& inv,
     digraphE<vec<int>>& D, vec<int>& dinv, const int LOOK_MERGE = 400,
     const Bool allow_two = False );

void LineLevelPullApart( const vec<int32_t>& bc, const HyperBasevectorX& hb,
     const vec<int>& inv, const ReadPathVec& paths, const vec<Bool>& dup,
     digraphE<vec<int>>& D, vec<int>& dinv );

void FlattenSomeBubbles( const HyperBasevectorX& hb, const vec<Bool>& dup,
     const ReadPathVec& paths, digraphE<vec<int>>& D, vec<int>& dinv, 
     ReadPathVec& dpaths, const int MAX_DELTA, const double MIN_RATIO,
     const int MAX_DEL );

#endif
