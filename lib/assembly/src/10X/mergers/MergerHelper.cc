// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

#include "10X/mergers/MergerHelper.h"
#include <omp.h>

template <class T>
void SwapE(T& a, T& b){
    T c = a;
    a = b;
    b = c;
}

template void SwapE(int& a, int& b);
template void SwapE(uint16_t& a, uint16_t& b);
template void SwapE(uchar& a, uchar& b);

template <class EL>
void GlueElements::JoinEs(vec<etuple<EL>>& e_inst, const vec<int64_t>& estart, 
        int64_t idx1, int64_t idx2){
    if(e_inst[idx1].eoID==e_inst[idx2].eoID) return;
    SwapE(e_inst[idx1].c,e_inst[idx2].c); 
    SwapE(e_inst[idx1].m,e_inst[idx2].m);
    int64_t idx = estart[e_inst[idx1].c] + e_inst[idx1].m;
    while(e_inst[idx].eoID!=e_inst[idx1].eoID){
        e_inst[idx].eoID = e_inst[idx1].eoID;
        idx = estart[e_inst[idx].c] + e_inst[idx].m;
    }
}
template void GlueElements::JoinEs(vec<etuple<uint16_t>>& e_inst, const vec<int64_t>& estart, 
        int64_t idx1, int64_t idx2);

template <class EL>
void GlueElements::JoinEs(vec<etuple<EL>>& e_inst, const vec<int64_t>& estart, 
        vec<int>& osize, int64_t idx1, int64_t idx2){
    if(e_inst[idx1].eoID==e_inst[idx2].eoID) return;
    SwapE(e_inst[idx1].c,e_inst[idx2].c); 
    SwapE(e_inst[idx1].m,e_inst[idx2].m);
    int s1 = osize[idx1], s2 = osize[idx2];
    int finsize = s1+s2;
    osize[idx1] = finsize;
    osize[idx2] = finsize;
    if(s1>=s2){
        int64_t idx = estart[e_inst[idx1].c] + e_inst[idx1].m;
        while(e_inst[idx].eoID!=e_inst[idx1].eoID){
            e_inst[idx].eoID = e_inst[idx1].eoID;
            osize[idx] = finsize;
            idx = estart[e_inst[idx].c] + e_inst[idx].m;
        }
    }else{
        int64_t idx = estart[e_inst[idx2].c] + e_inst[idx2].m;
        while(e_inst[idx].eoID!=e_inst[idx2].eoID){
            e_inst[idx].eoID = e_inst[idx2].eoID;
            osize[idx] = finsize;
            idx = estart[e_inst[idx].c] + e_inst[idx].m;
        }
    }
}

template void GlueElements::JoinEs(vec<etuple<uint16_t>>& e_inst, const vec<int64_t>& estart, 
        vec<int>& osize, int64_t idx1, int64_t idx2);


void GlueElements::JoinVs(vec<vtuple>& v_inst, const int idx1, const int idx2){
    if(v_inst[idx1].voID == v_inst[idx2].voID) return; // nothing to do
    SwapE(v_inst[idx1].nv,v_inst[idx2].nv);
    int idx = v_inst[idx1].nv; 
    while(v_inst[idx].voID!=v_inst[idx1].voID){
        v_inst[idx].voID = v_inst[idx1].voID;
        idx = v_inst[idx].nv; 
    }
}


// Help Flatten D's 
template <class T>
int All_Closures::elem(const digraphE<vec<T>>& D, int c, int m){
    return D.O(c)[m];
}
template int All_Closures::elem(const digraphE<vec<uchar>>& D, int c, int m);
template int All_Closures::elem(const digraphE<vec<int>>& D, int c, int m);

template <class T>
int All_Closures::isize(const digraphE<vec<T>>& D, int c){
    return D.O(c).isize();
}
template int All_Closures::isize(const digraphE<vec<uchar>>& D, int c);
template int All_Closures::isize(const digraphE<vec<int>>& D, int c);

template <class T>
unsigned int All_Closures::size(const digraphE<vec<T>>& D, int c){
    return D.O(c).size();
}
template unsigned int All_Closures::size(const digraphE<vec<uchar>>& D, int c);
template unsigned int All_Closures::size(const digraphE<vec<int>>& D, int c);
    
unsigned int All_Closures::size(const vec<int64_t>& estart, int c){
    return (unsigned int) estart[c+1]-estart[c];
}

int All_Closures::isize(const vec<int64_t>& estart, int c){
    return (int) estart[c+1]-estart[c];
}

template <class T>
vec<T> All_Closures::obj(const digraphE<vec<T>>& D, int c){
    return D.O(c);
}
template vec<uchar> All_Closures::obj(const digraphE<vec<uchar>>& D, int c);
template vec<int> All_Closures::obj(const digraphE<vec<int>>& D, int c);

template <class T>
int All_Closures::last(const digraphE<vec<T>>& D, int c){
    return D.O(c).back();
}
template int All_Closures::last(const digraphE<vec<uchar>>& D, int c);
template int All_Closures::last(const digraphE<vec<int>>& D, int c);

template <class T>
int All_Closures::first(const digraphE<vec<T>>& D, int c){
    return D.O(c).front();
}
template int All_Closures::first(const digraphE<vec<uchar>>& D, int c);
template int All_Closures::first(const digraphE<vec<int>>& D, int c);

template <class T>
void LeftEdges(const digraphE<vec<T>>& D, int64_t c, 
        const vec<int>& to_lefts, vec<int>& ds){
        ds.clear();
        int v = to_lefts[c];
        for (int j = 0; j<D.To(v).isize(); j++)
            ds.push_back(D.ITo(v,j));
}

template void LeftEdges(const digraphE<vec<uchar>>& D, int64_t c, 
        const vec<int>& to_lefts, vec<int>& ds);
template void LeftEdges(const digraphE<vec<int>>& D, int64_t c, 
        const vec<int>& to_lefts, vec<int>& ds);

template <class T>
void RightEdges(const digraphE<vec<T>>& D, int64_t c, 
        const vec<int>& to_rights, vec<int>& ds){
        ds.clear();
        int w = to_rights[c];
        for (int j = 0; j<D.From(w).isize(); j++)
            ds.push_back(D.IFrom(w,j));
}

template void RightEdges(const digraphE<vec<uchar>>& D, int64_t c, 
        const vec<int>& to_rights, vec<int>& ds);
template void RightEdges(const digraphE<vec<int>>& D, int64_t c, 
        const vec<int>& to_rights, vec<int>& ds);

template <class EL>
std::size_t checksum(const vec<etuple<EL>> & V) {
    std::size_t seed = 4;
    for(auto& q: V){
        seed ^= std::hash<int>()(q.c) + 0x9e3779b9 + 
            (seed << 6) + (seed >> 2);
        seed ^= std::hash<EL>()(q.m) + 0x9e3779b9 + 
            (seed << 6) + (seed >> 2);
        seed ^= std::hash<unsigned int>()(q.eoID) + 0x9e3779b9 + 
            (seed << 6) + (seed >> 2);
        seed ^= std::hash<char>()(q.type) + 0x9e3779b9 + 
            (seed << 6) + (seed >> 2);
    }
    return seed;
}
std::size_t checksum(const vec<etuple<uint16_t>> & V);

std::size_t checksum(const vec<vtuple> & V) {
    std::size_t seed = 4;
    for(auto& q: V){
        seed ^= std::hash<int64_t>()(q.v) + 0x9e3779b9 + 
            (seed << 6) + (seed >> 2);
        seed ^= std::hash<int>()(q.nv) + 0x9e3779b9 + 
            (seed << 6) + (seed >> 2);
        seed ^= std::hash<int>()(q.voID) + 0x9e3779b9 + 
            (seed << 6) + (seed >> 2);
        seed ^= std::hash<char>()(q.type) + 0x9e3779b9 + 
            (seed << 6) + (seed >> 2);
    }
    return seed;
}

std::size_t checksum(const vec<int>& V){
    std::size_t seed = 4;
    for(auto& e: V){
        seed ^= std::hash<int>()(e) + 0x9e3779b9 + 
            (seed << 6) + (seed >> 2);
    }
    return seed;
}


template <class T>
void InvertClosures(const digraphE<vec<T>>& D, vec<vec<int> >& ci, 
        const T nE){
    ci.clear();
    ci.resize(nE); 
    int64_t N = D.E();

    for ( int i = 0; i < N; i++ )
        for ( int j = 0; j < All_Closures::isize<T>(D,i) ; j++ )
            ci[ All_Closures::elem<T>(D,i,j) ].push_back(i);

    #pragma omp parallel for schedule( dynamic, 1000 )
    for ( T e = 0; e < nE; e++ )
    {   UniqueSort( ci[e] );}
}

template void InvertClosures(const digraphE<vec<int>>& D, vec<vec<int> >& ci, 
        const int nE);
template void InvertClosures(const digraphE<vec<uchar>>& D, vec<vec<int> >& ci, 
        const uchar nE);

template <class T>
void InvertClosures(const vec<vec<T>>& all_closures, vec<vec<int> >& ci, const T nE){
    ci.resize(nE);
    for ( int64_t i = 0; i < (int64_t) all_closures.size(); i++ )
        for ( int j = 0; j < all_closures[i].isize( ); j++ ){
            ci[ all_closures[i][j] ].push_back(i);
        }

    #pragma omp parallel for schedule( dynamic, 1000 )
    for ( T e = 0; e < nE; e++ )
    {   UniqueSort( ci[e] );}
}
template void InvertClosures(const vec<vec<int>>& all_closures, vec<vec<int>>& ci,
        const int nE);
template void InvertClosures(const vec<vec<uchar>>& all_closures,
        vec<vec<int>>& ci, const uchar nE);

template <class T>
void ComputeClosureIndex(const vec<int>& inv, vec<int64_t>& cinv, 
        vec<vec<int>>& ci, T hbE, const vec<vec<T>>& all_closures, 
        const Bool verbose)
{
    if (verbose) PRINTDEETS("indexing closures + uniquesorting");

    ci.resize(hbE);
    InvertClosures<T>(all_closures,ci,hbE);

    // Force match symmetry.
    if (verbose) PRINTDEETS("computing involution for all_closures");
    
    int64_t N = all_closures.size();
    cinv.resize(N,-1);
    #pragma omp parallel for
    for ( int64_t i = 0; i < N; i++ ){
        const vec<T>& x = all_closures[i];
        vec<T> rx( x.size( ) );
        for ( int j = 0; j < x.isize( ); j++ )
            rx[ x.isize( ) - j - 1 ] = inv[ x[j] ];
        int64_t ip = BinPosition( all_closures, rx );
        ForceAssertGe( ip, 0 );
        cinv[i] = ip;
    }

}
template void ComputeClosureIndex(const vec<int>& inv, vec<int64_t>& cinv, 
        vec<vec<int>>& ci, uchar hbE, 
        const vec<vec<uchar>>& all_closures, const Bool verbose);

template void ComputeClosureIndex(const vec<int>& inv, vec<int64_t>& cinv, 
        vec<vec<int>>& ci, int hbE, const vec<vec<int>>& all_closures, 
        const Bool verbose);
