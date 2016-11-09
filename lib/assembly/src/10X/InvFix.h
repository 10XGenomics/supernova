// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

#ifndef TENX_INVFIX_H
#define TENX_INVFIX_H

#include "CoreTools.h"
#include "paths/HyperBasevector.h"
#include "10X/astats/RefAlign.h"

void InvFix( const HyperBasevectorX& hb, const vec<int>& inv,
     const vec<int>& kmers, digraphE<vec<int>>& D, const vec<int>& dinv,
     vec<vec<vec<vec<int>>>>& dlines,
     const MasterVec<SerfVec<pair<int,int>>>& lbpx,
     const MasterVec<SerfVec<refalign>>& galigns );

#endif
