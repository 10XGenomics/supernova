// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

#ifndef TENX_MEASURE_GAPS_H
#define TENX_MEASURE_GAPS_H

#include "CoreTools.h"
#include "graph/Digraph.h"

void MeasureGaps( const int K, const vec<int>& kmers, const vec<int>& inv,
     const digraphE<vec<int>>& D, const vec<vec<vec<vec<int>>>>& dlines,
     const vec<int>& llens,
     const MasterVec< SerfVec<triple<int,int,int> > >& alignsb,
     vec< pair<int,int> >& gaps, const Bool verbose );

#endif
