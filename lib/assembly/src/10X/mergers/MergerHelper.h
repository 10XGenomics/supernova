// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
#ifndef TENX_MERGER_HELPER_H
#define TENX_MERGER_HELPER_H

#include "10X/mergers/MergerStructs.h"
#include "10X/mergers/NicePrints.h"

template <class T>
void InvertClosures(const digraphE<vec<T>>& D, vec<vec<int> >& ci, 
        const T nE);

template <class T>
void LeftEdges(const digraphE<vec<T>>& D, int64_t c, 
        const vec<int>& to_lefts, vec<int>& ds);

template <class T>
void RightEdges(const digraphE<vec<T>>& D, int64_t c, 
        const vec<int>& to_rights, vec<int>& ds);

template <class EL>
std::size_t checksum(const vec<etuple<EL>> & V);
std::size_t checksum(const vec<vtuple> & V);
std::size_t checksum(const vec<int>& V);

template <class T>
void InvertClosures(const digraphE<vec<T>>& D, vec<vec<int> >& ci, 
        const T nE);

template <class T>
void InvertClosures(const vec<vec<T>>& all_closures,vec<vec<int> >& ci,const T nE);

template <class T>
void ComputeClosureIndex(const vec<int>& inv, vec<int64_t>& cinv, 
        vec<vec<int>>& ci, T hbE, const vec<vec<T>>& all_closures, 
        const Bool verbose);

class GlueElements{
    // use only with primitive data types 
    public:
    template <class EL>
    static void JoinEs(vec<etuple<EL>>& e_inst, const vec<int64_t>& estart, 
            int64_t idx1, int64_t idx2);

    template <class EL>
    static void JoinEs(vec<etuple<EL>>& e_inst, const vec<int64_t>& estart, 
        vec<int>& osize, int64_t idx1, int64_t idx2);

    static void JoinVs(vec<vtuple>& v_inst, const int idx1, const int idx2);

};

// Help Flatten D's 
class All_Closures{
    public:
    template <class T>
    static int elem(const digraphE<vec<T>>& D, int c, int m);

    template <class T>
    static int isize(const digraphE<vec<T>>& D, int c);

    template <class T>
    static unsigned int size(const digraphE<vec<T>>& D, int c);
        
    static unsigned int size(const vec<int64_t>& estart, int c);

    static int isize(const vec<int64_t>& estart, int c);
    
    template <class T>
    static vec<T> obj(const digraphE<vec<T>>& D, int c);

    template <class T>
    static int last(const digraphE<vec<T>>& D, int c);

    template <class T>
    static int first(const digraphE<vec<T>>& D, int c);
};

#endif
