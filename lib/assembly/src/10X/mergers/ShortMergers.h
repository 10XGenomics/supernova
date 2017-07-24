// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

#ifndef TENX_MERGERS_SHORTMERGERS_H
#define TENX_MERGERS_SHORTMERGERS_H

#include "10X/Super.h"
#include "ParallelVecUtilities.h"
#include "10X/DfTools.h"
#include "10X/mergers/GlueGraphs.h"

void ExploreRightToDepth( const digraphE<vec<int>> & D, const vec<int> & to_right,
     const vec<Bool> & used, const int & v, vec<pair<uint8_t, int>> & de,
     const uint8_t DEPTH );

void MergeShortOverlapsLowMem( const HyperBasevectorX& hb, const vec<int>& inv, 
     digraphE<vec<int>>& D, vec<int>& dinv, const int LOOK_MERGE, const int LOOK, 
     const Bool allow_two, const Bool verbose, const Bool debug = False );

void ZipperRecursive( const HyperBasevectorX & hb, const vec<int> & inv, 
     digraphE<vec<int>> & D, vec<int> & dinv, const String DEBUG_DIR,
     const int MAX_DEPTH = 50, const int MAX_PASSES = 2, const int MIN_ZIPS = -1,
     const Bool verbose = False, const Bool gg_verbose = False);

void ZipperCombo( const HyperBasevectorX & hb, const vec<int> & inv, 
     digraphE<vec<int>> & D, vec<int> & dinv, const String DEBUG_DIR = ".", 
     const Bool verbose = False );

#endif
