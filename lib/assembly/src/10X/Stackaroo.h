// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

// Try to close gaps in supergraph assembly.

#ifndef TENX_STACKAROO_H
#define TENX_STACKAROO_H

#include "CoreTools.h"
#include "feudal/PQVec.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"

template <class VP, class VPI>
void Stackaroo( 

     // inputs:

     VirtualMasterVec<basevector> bases,
     VirtualMasterVec<PQVec> quals, const HyperBasevectorX& hb,
     const vec<int>& inv, VP & xpaths,
     VPI & xpaths_index, 
     const int MODE,

     // inputs and outputs:

     digraphE<vec<int>>& D, vec<int>& dinv, 

     // control over gap set:

     const String& S, String& R, const Bool ALL,

     // logging:

     const Bool VERBOSE, const int VERBOSITY, const HyperBasevectorX& ddn,
     const Bool VISUAL_ABBR, const Bool DIRECT, const Bool SHOW_MERGED_STACKS
     );

#endif
