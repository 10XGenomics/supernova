// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

#ifndef MAKE_LOCALS_TOOLS_H
#define MAKE_LOCALS_TOOLS_H

#include "CoreTools.h"
#include "feudal/PQVec.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "10X/paths/ReadPathVecX.h"

void FindLinesFixed( const digraphE<vec<int>>& D, const vec<int>& dinv,
     vec<vec<vec<vec<int>>>>& dlines, const int max_cell_paths,
     const int max_cell_depth, const Bool verbose = False );

/*
template<class T> 
*/
     void MakeLocalAssemblies( const HyperBasevectorX& hb,
     const vec<int>& inv, const vec<int64_t>& bci, 
     // T& pathsx, 
     ReadPathVecX& pathsx,
     const vec<Bool>& dup,
     VirtualMasterVec<basevector>& bases, VirtualMasterVec<PQVec>& quals,
     const vec<vec<int>>& bs, vec< digraphE<vec<int>> >& DL,
     vec< vec<int> >& DLI, vec< vec<int> >& BCLIST,
     vec<vec<vec<int>>>& BPATHS, const int MAX_BARCODES, const int MIN_LINE_TO_WALK,
     const String& udir, const String& suffix, 
     const String& DIR, const MasterVec< SerfVec<triple<int,int,int> > >& alignsb,
     Bool C2G=False, Bool STORE_LOCALS = False, const Bool PRINT_DETAILS = False,
     const String& INSTANCE = "undefined", const Bool NEW_STUFF = False,
     const Bool CANON = False );

void AddPairGaps( const HyperBasevectorX& hb, const vec<int> & inv,
     digraphE<vec<int>>& D, vec<int>& dinv, const ReadPathVec & dpaths, 
     Bool verbose=False );

void MergeBarcodeSets( vec<vec<int>>& bs, const vec<int64_t>& bci, 
     const int MAX_BARCODES );

void FillPairGaps( const HyperBasevectorX& hb, const vec<int>& inv,
     digraphE<vec<int>>& D, vec<int>& dinv, const ReadPathVecX& pathsx,
     const vec<int64_t>& bci, const vec<Bool>& dup );

void FBSCore( const String& DIR, const HyperBasevectorX& hb,
     const digraphE<vec<int>>& D, const vec<int>& dinv,
     const vec<vec<vec<vec<int>>>>& dlines,
     const int MIN_LINE_TO_WALK, const int MIN_KMERS, const int MAX_BARCODES,
     const int MAX_VIEW, vec<vec<int>>& bs, vec<int>& bsl, 
     const Bool FIX_BUG = False );

void FormBarcodeSets1( const String& DIR, const HyperBasevectorX& hb,
     const digraphE<vec<int>>& D, const vec<int>& dinv,
     const vec<vec<vec<vec<int>>>>& dlines, const vec<int64_t>& bci,
     const int MIN_LINE_TO_WALK, const int MIN_KMERS, const int MAX_BARCODES,
     const int MAX_VIEW, vec<vec<int>>& bs, 
     const MasterVec< SerfVec<triple<int,int,int> > >& alignsbf );

void FormBarcodeSets2( const String& DIR, const HyperBasevectorX& hb,
     const vec<int>& inv, const digraphE<vec<int>>& D, const vec<int>& dinv,
     const vec<vec<vec<vec<int>>>>& dlines, const vec<int64_t>& bci,
     const int MIN_LINE_TO_WALK, const int MIN_KMERS, const int MAX_BARCODES,
     const int MAX_VIEW, vec<vec<int>>& bs, const vec<int>& fin,
     const String& REGION, const vec<int>& linelist,
     const MasterVec< SerfVec<triple<int,int,int> > >& alignsb,
     const MasterVec< SerfVec<triple<int,int,int> > >& alignsbf,
     const Bool PRINT_DETAILS );

void EvalVersusFinished( const HyperBasevectorX& hb, const vec<int>& inv,
     const digraphE<vec<int>>& D, const vec<int>& dinv, const String& FIN,
     const vec<int>& fin, const String& DIR, const String& INDIR, const Bool VERB1 );

void PrintSummaryStats( const HyperBasevectorX& hb,
     const digraphE<vec<int>>& D, const vec<int>& dinv, ostream & out );

void ExtractSubsetOfGlobal( const HyperBasevectorX& hb, const digraphE<vec<int>>& D,
     const vec<int>& dinv,
     const vec<vec<vec<vec<int>>>>& dlines, const vec<pair<int,int>>& tol,
     const MasterVec< SerfVec<triple<int,int,int> > >& alignsb,
     const MasterVec< SerfVec<triple<int,int,int> > >& alignsbf,
     const vec<int>& fin, const String& REGION, const Bool PRINT_FETCH,
     vec<Bool>& keep );

void ParseRegion( const String& REGION, int& g_R, int& start_R, int& stop_R );

void PlaceReadsMasked2( const HyperBasevectorX& hb, const digraphE<vec<int>>& D,
     const vec<Bool>& dup, const ReadPathVecX& pathsx,
     const vec<Bool>& pmask, MasterVec<IntVec>& dpaths );

void PlaceReadsMasked( const HyperBasevectorX& hb, const digraphE<vec<int>>& D,
     const vec<Bool>& dup, const ReadPathVecX& pathsx,
     const vec<Bool>& pmask, ReadPathVec& dpaths );

void PlaceLinkedReadsMasked( const HyperBasevectorX & hb, const vec<int> & inv,
     const digraphE<vec<int>> & D, const vec<int> & dinv, 
     const vec<vec<vec<vec<int>>>> & dlines, const vec<Bool> & dup,
     const vec<int64_t> & bci, const ReadPathVecX & pathsx, ReadPathVec & dpaths,
     const vec<Bool> & pmask, const int MAX_PASSES, const Bool verbose );

void LBC2GBC( const vec<vec<int>>& BCLIST, const int A, const vec<int>& bcseq);

int GE2LA( const vec<digraphE<vec<int>>>& DL, const int64_t gE);

int GE2LA( const vec<int64_t>& dci, const int64_t gE);

#endif
