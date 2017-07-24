// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

#ifndef TENX_BASE_FIN_LOOKUP_H
#define TENX_BASE_FIN_LOOKUP_H

#include "MainTools.h"
#include "math/HoInterval.h"
#include "paths/HyperBasevector.h"

template<int K> void BaseFinLookup( const HyperBasevectorX& hb, 
     const vecbasevector& G,
     MasterVec< SerfVec< triple< ho_interval, int, int > > >& X );

template<int K> void BaseFinLookupSup( const HyperBasevectorX& hb, 
     const vecbasevector& tigs, const vecbasevector& G,
     MasterVec< SerfVec< triple< ho_interval, int, int > > >& X );

#endif
