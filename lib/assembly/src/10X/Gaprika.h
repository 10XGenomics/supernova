// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

#ifndef TENX_GAPRIKA_H
#define TENX_GAPRIKA_H

#include "CoreTools.h"
#include "10X/astats/RefAlign.h"

void Gaprika( const vec<int>& kmers, digraphE<vec<int>>& D, const vec<int>& dinv,
     const vec<vec<vec<vec<int>>>>& dlines, 
     const MasterVec<SerfVec<pair<int,int>>>& lbpx, const int MAX_GAP,
     const MasterVec<SerfVec<refalign>>& galigns, 
     const int VERBOSITY, const Bool EVAL );

#endif
