// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

#ifndef TENX_LAWNMOWER_H
#define TENX_LAWNMOWER_H

#include "CoreTools.h"
#include "feudal/PQVec.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "10X/paths/ReadPathVecX.h"
#include "feudal/ObjectManager.h"

template< class VB, class VQ, class VP > void MowLawn(

     // Inputs, read-only:

     HyperBasevectorX& hb, vec<int>& inv, VB& bases, VQ& quals,
     VP& paths, String paths_index_file, const vec<int32_t>& bc,
     vec<Bool>& dup, double interdup,

     // Output:

     vec<int>& dels,

     // For tiny experiments, a list of edges that might be deleted:

     const vec<int>& = vec<int>( ) );

template< class VB, class VQ > void MowLawn(

     // Inputs, read-only:

     HyperBasevectorX& hb, vec<int>& inv, VB& bases, ObjectManager<VQ>& quals_om,
     ReadPathVecX& paths, String paths_index_file, const vec<int32_t>& bc,
     vec<Bool>& dup, double interdup,

     // Output:

     vec<int>& dels,

     // For tiny experiments, a list of edges that might be deleted:

     const vec<int>& = vec<int>( ) );

#endif
