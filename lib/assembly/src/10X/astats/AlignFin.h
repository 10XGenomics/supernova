// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

#ifndef TENX_ALIGN_FIN
#define TENX_ALIGN_FIN

#include "CoreTools.h"
#include "PackAlign.h"
#include "paths/HyperBasevector.h"

void AlignFin( 

     // inputs:

     const HyperBasevectorX& hb, vec<int> inv, digraphE<vec<int>> D, vec<int> dinv,
     const vec<vec<vec<vec<int>>>>& dlines, 
     const vecbasevector& genome,  // finished sequence
     MasterVec< SerfVec<triple<int,int,int> > > galignsb,  // align->finished

     // outputs:

     vec< vec< vec< triple< vec<int>, align, int > > > >& Matches,
     double& errw,     // weighted error rate
     int64_t& N50PS,   // N50 perfect stretch
     String& report,

     // control:

     const Bool longest_only, // kill all but longest alignment
     const Bool stats = True, // compute N50 perfect stretch and weighted error rate
     const vec<int>& targets = vec<int>( ) // use only these finished ids

     );

#endif
