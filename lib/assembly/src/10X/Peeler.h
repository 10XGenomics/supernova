// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

#ifndef TENX_PEELER_H
#define TENX_PEELER_H

#include "CoreTools.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "10X/Super.h"

void Peeler( const HyperBasevectorX& hb, const vec<int>& inv,
     digraphE<vec<int>>& D, vec<int>& dinv, 
     ReadPathVec& dpaths, // should be current on input, not current on exit
     MasterVec<ULongVec>& dpaths_index, // ditto
     ostream& rout // deep logging
     );

#endif
