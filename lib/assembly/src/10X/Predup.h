// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

// Pre-deduplication.

#ifndef TENX_PREDUP_H
#define TENX_PREDUP_H

#include "Basevector.h"
#include "CoreTools.h"
#include "feudal/PQVec.h"

void Predup( vecbasevector& bases, VecPQVec& quals, vec<int64_t>& bci );

#endif
