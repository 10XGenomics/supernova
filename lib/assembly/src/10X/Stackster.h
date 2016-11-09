// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

#ifndef TENX_STACKSTER_H
#define TENX_STACKSTER_H

#include "CoreTools.h"
#include "feudal/PQVec.h"
#include "paths/long/ReadPath.h"
#include "paths/long/ReadStack.h"

template< class VB, class VQ, class VP, class VPI >
void Stackster( const int e1, const int e2, const vec<basevector>& edges,
     VB const& bases, VQ const& quals, const int K, const vec<DataSet>& datasets,
     const vec<int>& kmers, const vec<int>& inv, const vec<Bool>& dup,
     VP const& paths, VPI const& paths_index, vec<basevector>& closures, 
     vec<int>& trim, const int VERBOSITY, const Bool ALT, const Bool EXP,
     const vec< pair<int64_t,Bool> >& idsfw2 = vec< pair<int64_t,Bool> >( ) );

#endif
