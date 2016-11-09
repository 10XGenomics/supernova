// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

#ifndef TENX_RESCUE_H
#define TENX_RESCUE_H

#include "Basevector.h"
#include "CoreTools.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"

void Rescue( 

     // For experimentation, otherwise zero:

     const int64_t START, 
     
     // Inputs:

     const vecbasevector& bases,
     const vec<int32_t>& bc, const vec<DataSet>& datasets,
     const HyperBasevectorX& hb, const vec<int>& inv, const ReadPathVec& paths,

     // Output:

     vec<basevector>& patches );

#endif
