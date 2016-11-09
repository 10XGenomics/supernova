// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

#ifndef TENX_VIEW_H
#define TENX_VIEW_H

#include "CoreTools.h"
#include "graph/Digraph.h"
#include "math/HoInterval.h"

// View: find genomic locations of a line.  Note:
// 1. This only looks at non-bubble parts of the line.
// 2. It only looks at aligns to 1,...,22,X,Y.

template<class VA> void View( const int l, const int K, const vec<int>& kmers, 
     const vec<int>& inv, const digraphE<vec<int>>& D,
     const vec<vec<vec<vec<int>>>>& dlines,
     const vec<vec< pair<int,int> >>& linelocs,
     VA&& alignsb, SerfVec< quad<int,Bool,ho_interval,ho_interval> >& view );

#endif
