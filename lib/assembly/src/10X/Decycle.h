// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

#ifndef TENX_DECYCLE_H
#define TENX_DECYCLE_H

#include "CoreTools.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "10X/Gap.h"
#include "10X/Super.h"

void Decycle( const HyperBasevectorX& hb, const vec<int>& inv,
     digraphE<vec<int>>& D, vec<int>& dinv, const ReadPathVecX& paths );

#endif
