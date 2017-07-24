// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>

#include "10X/mergers/EdgeSupport.h"
#include "10X/mergers/NicePrints.h"
#include "10X/Super.h"
#include "10X/mergers/ExtremeAsm.h"


// merger ops between many elements into one
// this isn't necessarily well define considering that underlying structure
template <class EL>
void MergeSupport(vec<vec<quad<int,EL,int,EL>>>& support, 
        const vec<int>& picked, vec<quad<int,EL,int,EL>>& merged_el)
{
    ForceAssertLe(picked.size(),support.size());
    merged_el.clear();

    // aggregate
    for(auto p : picked)
        if(p < support.size())
            merged_el.append(support[p]);

    CleanupSupport<EL>(merged_el);
}

template void MergeSupport(vec<vec<quad<int,uint16_t,int,uint16_t>>>& support,
        const vec<int>& picked, vec<quad<int,uint16_t,int,uint16_t>>& mergerd_el);

// split ops from one element at various markers on the superedge
// this isn't necessarily well define around loops
template <class EL>
void SplitSupport(vec<vec<quad<int,EL,int,EL>>>& split_support, 
        const vec<int>& markers, const vec<quad<int,EL,int,EL>>& el)
{
   split_support.clear();
   split_support.resize(markers.size()+1); // n dividers create n+1 parts
   auto el1 = el;
   for(int m = 0; m < markers.size(); m++){
        for(auto& q : el1){
            if(q.third < markers[m] && (q.third+q.fourth >= markers[m])){ 
                auto q1 = q;
                q1.fourth = markers[m]-q.third;
                split_support[m].push_back(q1);
                q.third = q1.third+q1.fourth;
                q.fourth = q.fourth - q1.fourth;
            }
        }
   }
}

template void SplitSupport(vec<vec<quad<int,uint16_t,int,uint16_t>>>& split_support,
        const vec<int>& markers, const vec<quad<int,uint16_t,int,uint16_t>>& el);

template <class EL>
void EraseEdges(vec<vec<quad<int,EL,int,EL>>>& supp, const vec<Bool>& to_delete)
{
    if(to_delete.size()==0)
        return;
    // sanitize
    ForceAssertEq(to_delete.size(),supp.size());
    EraseIf(supp,to_delete);
}

template void EraseEdges(vec<vec<quad<int,uint16_t,int,uint16_t>>>& supp, const vec<Bool>& to_delete);


template <class EL>
void EraseEdges(vec<vec<quad<int,EL,int,EL>>>& supp, const vec<int>& to_delete)
{
    if(to_delete.size()==0)
        return;

    // sanitize
    auto to_del = to_delete;
    ParallelUniqueSort(to_del);
    ForceAssertLe(to_del.size(),supp.size());
    ForceAssertLe(Max(to_del),supp.size()-1);

    vec<Bool> remove(supp.size(),False);
    for(auto d: to_del)
        remove[d]=True;
    EraseIf(supp,remove);
}

template void EraseEdges(vec<vec<quad<int,uint16_t,int,uint16_t>>>& supp, const vec<int>& to_delete);

template <class EL>
void DeleteEdges(vec<vec<quad<int,EL,int,EL>>>& supp, const vec<int>& to_delete)
{
    if(to_delete.size()==0)
        return;

    // sanitize
    auto to_del = to_delete;
    ParallelUniqueSort(to_del);
    ForceAssertLe(to_del.size(),supp.size());
    ForceAssertLe(Max(to_del),supp.size()-1);

    #pragma omp parallel for
    for(int i = 0; i < to_del.size(); i++)
        supp[to_del[i]].clear();

}

template void DeleteEdges(vec<vec<quad<int,uint16_t,int,uint16_t>>>& supp, 
        const vec<int>& to_delete);

template <class EL>
void DeleteEdges(vec<vec<quad<int,EL,int,EL>>>& supp, const vec<Bool>& to_delete)
{
    if(to_delete.size()==0)
        return;

    // sanitize
    ForceAssertEq(to_delete.size(),supp.size());
    #pragma omp parallel for
    for(int i = 0; i < to_delete.size(); i++)
        if(to_delete[i])
            supp[i].clear();

}

template void DeleteEdges(vec<vec<quad<int,uint16_t,int,uint16_t>>>& supp, 
        const vec<Bool>& to_delete);

template <class EL>
void CleanupSupport(vec<vec<quad<int,EL,int,EL>>>& supp)
{
    #pragma omp parallel for
    for(int i = 0; i < supp.isize(); i++)
        CleanupSupport(supp[i]);
}
template void CleanupSupport(vec<vec<quad<int,uint16_t,int,uint16_t>>>& supp);

template <class EL>
void CleanupSupport(vec<quad<int,EL,int,EL>>& merged_el)
{
    // sort by reads and positions
    // don't uniq sort since appending a loop to itself will increase length of read
    ParallelUniqueSort(merged_el);
    if(merged_el.isize()<2) return;
    auto el = merged_el; merged_el.clear();
    merged_el.push_back(el[0]);

    for(int i = 1; i < el.isize(); i++){
        auto& last = merged_el.back();
        if(el[i].first!=last.first){ // new read
            merged_el.push_back(el[i]);
            continue;
        }// else try merging
        //make sure connectedness, also avoids merging in self loops
        if(el[i].second==last.second+last.fourth){
            if(el[i].third==last.third+last.fourth){ 
                last.fourth+=el[i].fourth;
                continue;
            }
        }
        // else it is a new segment
        merged_el.push_back(el[i]);
    }
    Sort(merged_el);
}

template void CleanupSupport(vec<quad<int,uint16_t,int,uint16_t>>& supp);

template <class EL>
void GetReadPathsAndIndex(vec<vec<quad<int,EL,int,EL>>>& r_supp,
    ReadPathVec& paths, vec<vec<int>>& paths_index, int nreads){
    paths_index.clear();
    paths_index.resize(r_supp.size());
    // get index
    #pragma omp parallel for
    for(int i = 0; i < r_supp.isize(); i++){
        paths_index[i].reserve(r_supp[i].size());
        // store first
        for(auto& q: r_supp[i])
            paths_index[i].push_back(q.first);
        UniqueSort(paths_index[i]);
        UniqueSort(r_supp[i]);
    }
    
    vec<quad<int,int,int,int>> matches;
    for(int i = 0; i < r_supp.isize(); i++)
        for(auto q: r_supp[i])
            matches.push(q.first,q.second,q.third,i);

    ParallelUniqueSort(matches); // uniq just to be safe form dups
    // read off path and offset
    int rd = -1;
    paths.resize(nreads);
    for(int i = 0; i < matches.isize(); i++){
        if(rd!=matches[i].first){
            rd = matches[i].first;
            paths[rd].setOffset(matches[i].third);
        }
        paths[rd].push_back(matches[i].fourth);
    }
}

template void GetReadPathsAndIndex(
        vec<vec<quad<int,uint16_t,int,uint16_t>>>& r_supp,
        ReadPathVec& paths, vec<vec<int>>& paths_index, int nreads);

template <class EL>
void AppendSupport(vec<quad<int,EL,int,EL>>& supp, const vec<quad<int,EL,int,EL>>& add)
{
    // aggregate
    supp.append(add);
    CleanupSupport<EL>(supp);
}

template void AppendSupport(vec<quad<int,uint16_t,int,uint16_t>>& supp, 
        const vec<quad<int,uint16_t,int,uint16_t>>& add);

template <class EL, class T>
void RCSupport(vec<quad<int,EL,int,EL>>& supp, vec<quad<int,EL,int,EL>>& rcsupp,
        const vec<int64_t>& cinv, const vec<vec<T>>& cl, int Ds)
{
   for(auto q: supp){
        auto rcq = q;
        rcq.first = cinv[rcq.first];
        rcq.second = cl[rcq.first].size() - q.fourth - q.second;
        rcq.third = Ds - (q.fourth%Ds) - q.third;  // the % is to take care of loops
        rcsupp.push_back(rcq);
   }
   ParallelUniqueSort(rcsupp);
}

template <class EL, class T>
void RCSupport(vec<quad<int,EL,int,EL>>& supp, const vec<int64_t>& cinv,
        const vec<vec<T>>& cl, int Ds)
{
    auto rcsupp = supp;
    rcsupp.clear();
    RCSupport<EL>(supp,rcsupp,cinv,cl,Ds);
}

template void RCSupport(vec<quad<int,uint16_t,int,uint16_t>>& supp, 
        vec<quad<int,uint16_t,int,uint16_t>>& rcsupp, const vec<int64_t>& cinv,
        const vec<vec<uchar>>& cl, int Ds);

template void RCSupport(vec<quad<int,uint16_t,int,uint16_t>>& supp,
        const vec<int64_t>& cinv, const vec<vec<uchar>>& cl, int Ds);


