// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

#ifndef TENX_LOCAL_H
#define TENX_LOCAL_H

#include "CoreTools.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "10X/paths/ReadPathVecX.h"

void GlueAssemblies( const HyperBasevectorX& hb, digraphE<vec<int>>& D, 
     vec<int>& dinv, const vec<int> mult,
     const int MIN_LINE_TO_WALK, const int MAX_MULT, const Bool verbose = True );

void BigPull( const HyperBasevectorX& hb, const vec<int>& inv,
     digraphE<vec<int>>& D, vec<int>& dinv, const ReadPathVec& paths,
     const vec<Bool>& dup, const Bool verbose = False );

void FindRelevantPaths ( digraphE<vec<int>> & D, vec<int> & dinv, ReadPathVecX & pathsx,
     VecULongVec & dpaths_index, vec<int64_t> & bci, 
     vec <Bool> & keep, vec < vec<int> > & bs, vec <Bool> & pmask );

#endif
