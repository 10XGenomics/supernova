// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

#ifndef TENX_CAPTURE_H
#define TENX_CAPTURE_H

#include "CoreTools.h"
#include "graph/Digraph.h"
#include "paths/HyperBasevector.h"

void CaptureMultiLoops( digraphE<vec<int>>& D, vec<int>& dinv, vec<int>& dels,
     const Bool verbose, const Bool single );

void CaptureMessyLoops( const HyperBasevectorX& hb, const vec<int>& inv,
     digraphE<vec<int>>& D, vec<int>& dinv, vec<int>& dels, 
     const Bool allow_point = False, const int LONG_LINE = 10000 );

void CaptureSimpleLoops( digraphE<vec<int>>& D, vec<int>& dinv, vec<int>& dels,
     const Bool verbose, const Bool single );

void CaptureCanonicalLoops( digraphE<vec<int>>& D, vec<int>& dinv, vec<int>& dels,
     const Bool verbose, const Bool single );

#endif
