// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

#ifndef TENX_PULL_APART_H
#define TENX_PULL_APART_H

#include "CoreTools.h"
#include "graph/Digraph.h"
#include "paths/long/ReadPath.h"
#include "10X/IntIndex.h"

void PullApartInversions( digraphE<vec<int>>& D, vec<int>& dinv );

void PullApart( const vec<int>& inv, digraphE<vec<int>>& D, vec<int>& dinv, 
     const ReadPathVec& dpaths, vec<int>& dels, Bool verbose );

#endif
