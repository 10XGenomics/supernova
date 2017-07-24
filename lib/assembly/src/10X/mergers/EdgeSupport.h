// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

#ifndef TENX_EDGESUPP_H
#define TENX_EDGESUPP_H

#include "10X/mergers/MergerStructs.h"

template <class EL>
void MergeSupport(vec<vec<quad<int,EL,int,EL>>>& support,
        const vec<int>& picked, vec<quad<int,EL,int,EL>>& merged_el);

template <class EL>
void SplitSupport(vec<vec<quad<int,EL,int,EL>>>& split_support, 
        const vec<int>& markers, const vec<quad<int,EL,int,EL>>& el);

template <class EL>
void DeleteEdges(vec<vec<quad<int,EL,int,EL>>>& supp, const vec<int>& to_delete);

template <class EL>
void DeleteEdges(vec<vec<quad<int,EL,int,EL>>>& supp, const vec<Bool>& to_delete);

template <class EL>
void EraseEdges(vec<vec<quad<int,EL,int,EL>>>& supp, const vec<int>& to_delete);

template <class EL>
void EraseEdges(vec<vec<quad<int,EL,int,EL>>>& supp, const vec<Bool>& to_delete);

template <class EL>
void CleanupSupport(vec<quad<int,EL,int,EL>>& supp);

template <class EL>
void CleanupSupport(vec<vec<quad<int,EL,int,EL>>>& supp);

template <class EL>
void AppendSupport(vec<quad<int,EL,int,EL>>& supp, const vec<quad<int,EL,int,EL>>& add);

template <class EL, class T>
void RCSupport(vec<quad<int,EL,int,EL>>& supp, vec<quad<int,EL,int,EL>>& rcsupp,
        const vec<int64_t>& cinv, const vec<vec<T>>& cl, int Ds);

template <class EL, class T>
void RCSupport(vec<quad<int,EL,int,EL>>& supp, const vec<int64_t>& cinv, 
        const vec<vec<T>>& cl, int Ds);

template <class EL>
void GetReadPathsAndIndex(vec<vec<quad<int,EL,int,EL>>>& r_supp,
    ReadPathVec& paths, vec<vec<int>>& paths_index, int nreads);

#endif
