// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
#include "Basevector.h"
#include "Intvector.h"
#include "MainTools.h"
#include "VecUtilities.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "ParallelVecUtilities.h"
#include "10X/paths/ReadPathVecX.h"

void writePathsIndex( ReadPathVecX & paths, const HyperBasevectorX& hb,
        vec<int> & inv, String base_dir, String pi_file,
        String counts_file, const int chunks, bool verbose);

void writePathsIndex( ReadPathVec & paths,
        vec<int> & inv, String base_dir, String pi_file,
        String counts_file, const int chunks, bool verbose);

void computeEdgeToBarcodeX(const ReadPathVecX & paths,const HyperBasevectorX & hb, 
        const vec <int32_t> & bc, const vec<int> & inv, 
        const vec<int64_t> & bci, VecIntVec & ebcx, bool verbose );

