// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

#ifndef TENX_CLOSURES_TO_GRAPH_H
#define TENX_CLOSURES_TO_GRAPH_H

#include "10X/mergers/MergerStructs.h"
#include "10X/mergers/MergerHelper.h"
#include "10X/mergers/NucleateGraph.h"
#include "10X/mergers/NucleateGraph_HP.h"

// Convert a digraphE<int> into a digraphE<vec<int>> in which unneeded vertices
// have been removed.

void Vectorify(

     // input: pointer to digraphE<int> -- GETS DELETED MIDSTREAM!!

     digraphE<int>* D0p,

     // input and output:

     vec<int>& dinv,

     // output:

     digraphE<vec<int>>& D,

     // control:

     const Bool verbose, const Bool single, const Bool use_inv );

template <class T>
void ClosuresToGraph( const HyperBasevectorX& hb, const vec<int>& inv, 
     vec<vec<T>>& all_closures, digraphE<vec<T>>& D, vec<int>& dinv,
     const Bool verbose, String dir="", Bool Canon = False);

template <class T>
void ClosuresToGraph_HP( const HyperBasevectorX& hb, const vec<int>& inv, 
     vec<vec<T>>& all_closures, digraphE<vec<T>>& D, vec<int>& dinv,
     const Bool verbose, String dir="", Bool Canon = False);


#endif
