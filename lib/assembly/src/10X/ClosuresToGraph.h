// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

#ifndef TENX_CLOSURES_TO_GRAPH_H
#define TENX_CLOSURES_TO_GRAPH_H

#include "CoreTools.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "10X/ClosuresToGraph.h"
#include <unordered_set>
#include <unordered_map>
#include "ParallelVecUtilities.h"
#include "graph/Digraph.h"

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

void ClosuresToGraph( const HyperBasevectorX& hb, const vec<int>& inv, 
     vec<vec<int>>& all_closures, digraphE<vec<int>>& D, vec<int>& dinv,
     const Bool verbose, String dir="" );

class etuple{
    public:
    int c;
    unsigned short m;
    unsigned int eoID;
    char type;

    etuple(){};

    etuple(const int ci, const unsigned short mi, const unsigned int eoIDi,
            const char typei){
                c = ci; m = mi; eoID = eoIDi; type = typei; }

    etuple(const etuple& source){
                c = source.c; m = source.m; eoID = source.eoID; type = source.type; }

    etuple& operator=(const etuple& source){
                if (this == &source)
                return *this;
                c = source.c; m = source.m; eoID = source.eoID; type = source.type;
                return *this;}

    bool operator==(const etuple& source){
                if(this == &source)
                    return true;
                if((c==source.c) && (m==source.m) && 
                        (eoID==source.eoID) && (type==source.type))
                    return true;
                return false;
    }
};

TRIVIALLY_SERIALIZABLE(etuple);

class vtuple{
    public:
    int64_t v;
    int nv;
    int voID;
    char type;

    vtuple(){};

    vtuple(const int64_t v1, const int nv1, const int voID1,
            const char type1){
                v = v1; nv = nv1; voID = voID1; type = type1; }

    vtuple(const vtuple& source){
                v = source.v; nv = source.nv; voID = source.voID; type = source.type; }

    vtuple& operator=(const vtuple& source){
                if (this == &source)
                return *this;
                v = source.v; nv = source.nv; voID = source.voID; type = source.type;
                return *this;}

    bool operator==(const vtuple& source){
                if(this == &source)
                    return true;
                if(v == source.v && nv == source.nv && voID == source.voID && type == source.type)
                    return true;
                return false;
    }
};

TRIVIALLY_SERIALIZABLE(vtuple);

namespace std{
    template <>
    struct hash<vec<int>>{ 
        public:
        std::size_t operator()(const vec<int> & vec) const {
            std::size_t seed = vec.size();
            for(auto& v : vec) {
                seed ^= std::hash<int>()(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            }
            return seed;
        }
    };

    template <>
    struct hash<vec<int64_t>>{ 
        public:
        std::size_t operator()(const vec<int64_t> & vec) const {
            std::size_t seed = vec.size();
            for(auto& v : vec) {
                seed ^= std::hash<int64_t>()(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            }
            return seed;
        }
    };
}

class om {

     public:

     om( ) { }
     om( const int c2, const short start1, const short start2, const short len ) :
          c2(c2), start1(start1), start2(start2), len(len) { }

     int c2;
     short start1, start2;
     short len;

     short Offset( ) const { return start1 - start2; }

     short Start1( ) const { return start1; }
     short Start2( ) const { return start2; }
     short Stop1( ) const { return start1 + len; }
     short Stop2( ) const { return start2 + len; }
     short Len() const {return len;}

     void Extend( const vec<int>& C1, const vec<int>& C2 )
     {    int n1 = C1.size( ), n2 = C2.size( );
          while( start1 > 0 && start2 > 0 && C1[start1-1] == C2[start2-1] )
          {    start1--;
               start2--;
               len++;    }
          while( start1 + len < n1 && start2 + len < n2
               && C1[start1+len] == C2[start2+len] )
          {    len++;    }    }

     void Validate( const vec<int>& C1, const vec<int>& C2 )
     {    for ( int i = 0; i < len; i++ )
          {    if ( C1[ Start1( ) + i ] != C2[ Start2( ) + i ] )
               {    cout << "\nInvalid om." << endl;
                    TracebackThisProcess( );    }    }    }

     friend Bool operator==( const om& o1, const om& o2 )
     {    return o1.c2 == o2.c2 && o1.start1 == o2.start1 && o1.start2 == o2.start2 
               && o1.len == o2.len;    }

};

TRIVIALLY_SERIALIZABLE(om);

#endif
