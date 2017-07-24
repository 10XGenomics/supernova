// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

#ifndef TENX_PROJECT_RUNWAY_H
#define TENX_PROJECT_RUNWAY_H

#include "CoreTools.h"
#include "ParallelVecUtilities.h"
#include "math/HoInterval.h"
#include "paths/HyperBasevector.h"
#include "10X/astats/ProjectRunway.h"

void GetMaxPerfs(

     // inputs:
     const vec<int>& gs, const HyperBasevectorX& hb, const vec<int>& kmers,
     const digraphE<vec<int>>& D, const vec<int>& to_left, const vec<int>& to_right,
     const vec<int>& dlens, const vecbasevector& tigs, const vecbasevector& genomef,
     const MasterVec< SerfVec< triple< ho_interval, int, int > > >& flocs0,
     const int minlen, const int minlen_force,
     // output:
     vec<vec<vec<triple<ho_interval,int,int>>>>& X );

#endif
