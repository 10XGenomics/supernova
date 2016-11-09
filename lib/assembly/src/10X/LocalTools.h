// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

#ifndef TENX_LOCAL_H
#define TENX_LOCAL_H

#include "CoreTools.h"
#include "paths/HyperBasevector.h"

void GlueAssemblies( const HyperBasevectorX& hb, digraphE<vec<int>>& D, 
     vec<int>& dinv, const vec<int> mult,
     const int MIN_LINE_TO_WALK, const int MAX_MULT, const Bool verbose = True );

#endif
