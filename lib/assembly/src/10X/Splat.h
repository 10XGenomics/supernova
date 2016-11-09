// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

// Patch in original gap closures.

#ifndef TENX_SPLAT_H
#define TENX_SPLAT_H

#include "CoreTools.h"
#include "paths/HyperBasevector.h"

void Splat( const HyperBasevectorX& hb, const vec<int>& inv,
     const vecbasevector& closures, digraphE<vec<int>>& D, vec<int>& dinv );

#endif
