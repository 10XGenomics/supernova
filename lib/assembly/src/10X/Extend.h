// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

#ifndef TENX_EXTEND_H
#define TENX_EXTEND_H

#include "10X/Extend.h"
#include "feudal/PQVec.h"
#include "graph/DigraphTemplate.h"
#include "paths/HyperBasevector.h"
#include "paths/long/large/GapToyTools.h"
#include "10X/paths/ReadPathVecX.h"

template <class VQ>
void ExtendPathsNew( const HyperBasevector& hb, const vec<int>& inv,
     const vecbasevector& bases, const VQ& qualsx, ReadPathVecX& paths,
     const Bool BACK_EXTEND );
#endif
