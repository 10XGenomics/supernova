// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

#ifndef TENX_CLOSOMATIC_H
#define TENX_CLOSOMATIC_H

#include "CoreTools.h"
#include "feudal/PQVec.h"
#include "math/Functions.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "paths/long/ReadStack.h"
#include "10X/DfTools.h"
#include <unordered_map>

void FindEdgePairs( const HyperBasevectorX& hb, const vec<int>& inv,
     MasterVec<ReadPath>& paths, String pi_file,
     const vec<Bool>& bad,
     vec< pair<int,int> >& pairs, const vec<DataSet>& datasets,
     const vec<int32_t>& bc, const Bool one_good );

template< class VB, class VQ, class VP, class VPI >
void CloseGap( const HyperBasevectorX& hb, const vec<int>& inv, 
     VB& bases, VQ& quals, VP& paths, VPI& paths_index,
     const int e1, const int e2, vec<basevector>& closures, Bool verbose, 
     const int pi, const int max_width );

template< class VB, class VQ, class VP, class VPI >
void CloseGap2( const HyperBasevectorX& hb, const vec<int>& inv, 
     VB& bases, VQ& quals, VP& paths, VPI& paths_index,
     const int e1, const int e2, vec<basevector>& closures, Bool verbose, 
     const int pi, const int max_width );

#endif
