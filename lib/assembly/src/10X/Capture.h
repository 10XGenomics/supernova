// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

#ifndef TENX_CAPTURE_H
#define TENX_CAPTURE_H

#include "CoreTools.h"
#include "graph/Digraph.h"
#include "paths/HyperBasevector.h"

void CaptureLoops( const HyperBasevectorX& hb, const vec<int>& inv,
     digraphE<vec<int>>& D, vec<int>& dinv, const Bool verbose, const Bool single );

void CaptureSimpleLoops( digraphE<vec<int>>& D, vec<int>& dinv, vec<int>& dels,
     const Bool verbose, const Bool single );

void CaptureMultiLoops( digraphE<vec<int>>& D, vec<int>& dinv, vec<int>& dels,
     const Bool verbose, const Bool single );

void CaptureMessyLoops( const HyperBasevectorX& hb, const vec<int>& inv,
     digraphE<vec<int>>& D, vec<int>& dinv, vec<int>& dels, 
     const Bool allow_point = False, const int LONG_LINE = 10000,
     const int MAX_EDGE_IN_LOOP = 2000 );

void CaptureMessyLoops2( const HyperBasevectorX& hb, const vec<int>& inv,
     digraphE<vec<int>>& D, vec<int>& dinv, const int MAX_EDGE_IN_LOOP = 2000, 
     const int END_SEARCH = 10, const int MAX_MESS = 20,
     const Bool verbose = False );

void CaptureCanonicalLoops( digraphE<vec<int>>& D, vec<int>& dinv, vec<int>& dels,
     const Bool verbose, const Bool single );

#endif
