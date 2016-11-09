// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

// Align an assembly to a genome reference sequence.  

#ifndef TENX_GENOME_ALIGN_H
#define TENX_GENOME_ALIGN_H

#include "MainTools.h"
#include "paths/HyperBasevector.h"

void AlignToGenomeCore( const vecbasevector& tigs, const vecbasevector& genome,
     vec< vec< pair<int,int> > >& hits, const int K2 );

template<int K> void GenomeAlign( const HyperBasevectorX& hb, const vec<int>& inv,
     const vecbasevector& genome, 
     MasterVec< SerfVec<triple<int,int,int> > >& alignsb );

#endif
