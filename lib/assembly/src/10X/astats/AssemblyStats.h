// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

#ifndef TENX_ASSEMBLY_STATS_H
#define TENX_ASSEMBLY_STATS_H

#include "CoreTools.h"
#include "graph/Digraph.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "math/HoInterval.h"
#include "10X/paths/ReadPathVecX.h"

void ReportAssemblyStats( const vec<int64_t>& bci, const vecbasevector& genome,
     const vec< pair<int,ho_interval> >& ambint,
     const HyperBasevectorX& hb, const vec<int>& inv,
     digraphE<vec<int>> D, const vec<int> dinv, 
     const vec<vec<vec<vec<int>>>>& lines,
     MasterVec< SerfVec<triple<int,int,int> > >& alignsb,
     const vecbasevector& G, ostream& out, const String& DIR, const String& OUTDIR );

#endif
