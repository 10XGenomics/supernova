// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

#ifndef TENX_REF_LOOKUP_H
#define TENX_REF_LOOKUP_H

#include "CoreTools.h"
#include "paths/HyperBasevector.h"

void FetchFinished( const String& SAMPLE, vecbasevector& G, int N = 0,
     const int seed = 0 );

String SampleGenome( const String& SAMPLE );
#endif
