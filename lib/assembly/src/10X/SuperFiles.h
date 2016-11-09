// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

// Generate full set of super-assembly files.  Starts with only the super-assembly
// data structures D and dinv.

#ifndef TENX_SUPER_FILESX_H
#define TENX_SUPER_FILESX_H

#include "CoreTools.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "10X/paths/ReadPathVecX.h"

void SuperFiles( const HyperBasevectorX& hb, const vec<int>& inv,
     const vec<int64_t>& bci,
     const vec<int32_t>& bc,
     const ReadPathVecX& paths, const vec<Bool>& dup,
     const digraphE<vec<int>>& D, const vec<int>& dinv,
     ReadPathVec& dpaths,              // recomputed!
     vec<vec<vec<vec<int>>>>& dlines,  // recomputed!
     vec<double>& COV,                 // recomputed!
     const VecIntVec& ebcx, const vec< triple<int,int,int> >& qept,
     const vecbasevector& genome,
     MasterVec< SerfVec<triple<int,int,int> > >& alignsb,
     const String& DIR, const String& WRITE_SUB );

// generate unique read counts for dpaths
void GenerateDpathsCounts(ReadPathVec const& dpaths, VecULongVec const& dpaths_index, 
          vec<int> const& dinv, vec<uint16_t>& dpaths_counts);

#endif
