// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

#ifndef SCAFFOLD_H
#define SCAFFOLD_H

#include "CoreTools.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "10X/DfTools.h"

void Scaffold( const HyperBasevectorX& hb, const vec<int>& inv,
     const VecIntVec& ebcx, digraphE<vec<int>>& D, vec<int>& dinv, 
     const ReadPathVec& dpaths, const vec<int64_t>& bid, 
     const vec<DataSet>& datasets, const Bool verbose, String& link_report,
     const Bool single );

#endif
