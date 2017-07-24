// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

#ifndef TENX_GRAPH_EDITORS_H
#define TENX_GRAPH_EDITORS_H

#include "10X/mergers/MergerStructs.h"

void CallGraphEditors(digraphE<vec<unsigned char>>& D, vec<int>& dinv, 
        vec<vec<unsigned char>>& all_closures, const vec<int64_t>& cinv,
        vec<vec<quad<int,uint16_t,int,uint16_t>>>& r_supp, const Bool verbose=False);
 
template<class T>
void RemoveSimpleHangs(digraphE<vec<T>>& D, vec<int>& dinv, vec<int>& dels);

template <class T>
void BabyZipper(digraphE<vec<T>>& D, vec<int>& dinv, T nE, const int max_pass, const int min_zips );

#endif
