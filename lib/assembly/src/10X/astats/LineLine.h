// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

#ifndef TENX_ASTATS_LINE_LINE_H
#define TENX_ASTATS_LINE_LINE_H

#include "CoreTools.h"

void FindLineLines(
     // inputs:
     const digraphE<vec<int>>& D, const vec<int>& dinv,
     const vec<vec<vec<vec<int>>>>& dlines,
     // outputs:
     vec<vec<vec<vec<int>>>>& dlines2,
     vec<int>* linv2p = nullptr );

void GetLineLineLengths( const vec<int>& llens,
     const vec<vec<vec<vec<int>>>>& dlines2, vec<int>& lens2 );

#endif
