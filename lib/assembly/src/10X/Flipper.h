// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

#ifndef TENX_FLIPPER_H
#define TENX_FLIPPER_H

#include "CoreTools.h"
#include "paths/HyperBasevector.h"
#include "10X/IntIndex.h"

void Flipper( const HyperBasevectorX& hb, const vec<int>& inv, const vec<Bool>& dup,
     const vec<Bool>& bad, const vec<int32_t>& bc, digraphE<vec<int>>& D, 
     vec<int>& dinv, const vec<vec<vec<vec<int>>>>& dlines,
     const IntIndex& dpaths_index, vec<String>& report, const int test_line = -1 );

#endif
