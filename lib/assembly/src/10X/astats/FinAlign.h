// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

#ifndef TENX_ASTATS_FIN_ALIGN_H
#define TENX_ASTATS_FIN_ALIGN_H

#include "CoreTools.h"
#include "math/HoInterval.h"
#include "10X/astats/RefLookup.h"

String FinAlign( const String& suffix, const vec<int>& gs, const Bool VERBOSE, 
     const String& INDIR, const HyperBasevectorX& hb, vec<int> inv, 
     vecbasevector tigs, digraphE<vec<int>> D, vec<int> dinv,
     vec<vec<vec<vec<int>>>> dlines, const vecbasevector& genomef,
     const MasterVec< SerfVec< triple< ho_interval, int, int > > >& flocs,
     vec<vec<int>>& statslist, const Bool write_files = True,
     const Bool PS = False // Compute the N50 perfect stretch
     );

#endif
