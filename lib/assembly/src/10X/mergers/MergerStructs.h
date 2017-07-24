// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

#ifndef TENX_MERGER_STRUCTS_H
#define TENX_MERGER_STRUCTS_H

#include "CoreTools.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include <unordered_set>
#include <unordered_map>
#include "ParallelVecUtilities.h"
#include "graph/Digraph.h"
#include "10X/paths/ReadPathVecX.h"
#include "10X/BuildLocal.h"
#include "10X/Super.h"
#include "feudal/PQVec.h"
#include "graph/DigraphTemplate.h"
#include "kmers/KmerRecord.h"
#include "paths/long/large/Lines.h"
#include "10X/DfTools.h"
#include "10X/Gap.h"
#include "10X/Heuristics.h"
#include "10X/Scaffold.h"
#include "10X/SecretOps.h"
#include "paths/long/large/GapToyTools.h"
#include "10X/PlaceReads.h"
#include "paths/RemodelGapTools.h"


typedef triple<pair<int,int>, pair<int,int>, int> nptuple;

template <class EL>
class etuple{
    public:
    int c;
    EL m;
    unsigned int eoID;
    char type;

    etuple(){};

    etuple(const int ci, const EL mi, const unsigned int eoIDi,
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

TRIVIALLY_SERIALIZABLE(etuple<int16_t>);
TRIVIALLY_SERIALIZABLE(etuple<int>);

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

struct qtuple{
    int d1;
    int a1;
    int d2;
    int a2;

    qtuple(){};

    qtuple(const int D1, const int A1, const int D2, const int A2){
        d1 = D1, d2 = D2, a1 = A1, a2 = A2;}

    qtuple(const qtuple& source){
        d1 = source.d1, d2 = source.d2, a1 = source.a1, a2 = source.a2;}

    qtuple& operator=(const qtuple& source){
        if(this == & source)
            return *this;
        d1 = source.d1, d2 = source.d2, a1 = source.a1, a2 = source.a2;
        return *this;
    }

    void swap(){
        int t2 = a1;
        int h2 = d1;
        a1 = a2;
        d1 = d2;
        a2 = t2;
        d2 = h2;
    }

    bool operator==(const qtuple& source)const{
        if(this==&source)
            return true;
        if(d1 == source.d1 && d2 == source.d2 && a1 == source.a1 && a2 == source.a2 )
            return true;
        return false;
    }

    bool operator<(const qtuple& source)const{
        if(a1 < source.a1)
            return true;
        if(d1 < source.d1 && a1 == source.a1)
            return true;
        if(d1 == source.d1 && a1==source.a1 && a2 < source.a2)
            return true;
        if(d1 == source.d1 && a1==source.a1 && d2<source.d2 && a2 == source.a2)
            return true;
        return false;
    }

};

TRIVIALLY_SERIALIZABLE(qtuple);

struct ptuple{
    int d1;
    int p1;
    int d2;
    int p2;
    int len;

    ptuple(){};

    ptuple(const int D1, const int P1, const int D2, const int P2, const int L){
        d1 = D1, d2 = D2, p1 = P1, p2 = P2; len = L;}

    ptuple(const ptuple& source){
        d1 = source.d1, d2 = source.d2, p1 = source.p1, p2 = source.p2; len = source.len;}

    ptuple& operator=(const ptuple& source){
        if(this == & source)
            return *this;
        d1 = source.d1, d2 = source.d2, p1 = source.p1, p2 = source.p2; len = source.len;
        return *this;
    }

    void swap(){
        int t2 = p1;
        int h2 = d1;
        p1 = p2;
        d1 = d2;
        p2 = t2;
        d2 = h2;
    }

    bool operator==(const ptuple& source)const{
        if(this==&source)
            return true;
        if(d1 == source.d1 && d2 == source.d2 && 
                p1 == source.p1 && p2 == source.p2 && len == source.len)
            return true;
        return false;
    }

    bool operator<(const ptuple& source)const{
        if(d1 < source.d1)
            return true;
        if(d1 == source.d1 && p1 < source.p1)
            return true;
        if(d1 == source.d1 && p1==source.p1 && d2 < source.d2)
            return true;
        if(d1 == source.d1 && p1==source.p1 && d2==source.d2 && p2 < source.p2)
            return true;
        if(d1 == source.d1 && p1==source.p1 && d2==source.d2 && p2 == source.p2 
                && len < source.len)
            return true;
        return false;
    }

};

TRIVIALLY_SERIALIZABLE(ptuple);

template <class EL>
class om {

     public:

     om( ) { }
     om( const int c2, const EL start1, const EL start2, const EL len ) :
          c2(c2), start1(start1), start2(start2), len(len) { }

     int c2;
     EL start1, start2;
     EL len;

     int Offset( ) const { return ((int)start1) - ((int)start2); }

     EL Start1( ) const { return start1; }
     EL Start2( ) const { return start2; }
     EL Stop1( ) const { return start1 + len; }
     EL Stop2( ) const { return start2 + len; }
     EL Len() const {return len;}

     template <class T>
     void Extend( const vec<T>& C1, const vec<T>& C2 )
     {    int n1 = C1.size( ), n2 = C2.size( );
          while( start1 > 0 && start2 > 0 && C1[start1-1] == C2[start2-1] )
          {    start1--;
               start2--;
               len++;    }
          while( start1 + len < n1 && start2 + len < n2
               && C1[start1+len] == C2[start2+len] )
          {    len++;    }    }

     template <class T>
     void Validate( const vec<T>& C1, const vec<T>& C2 )
     {    for ( int i = 0; i < len; i++ )
          {    if ( C1[ Start1( ) + i ] != C2[ Start2( ) + i ] )
               {    cout << "\nInvalid om." << endl;
                    TracebackThisProcess( );    }    }    }

    bool operator==(const om& source){
         if(this==&source)
             return true;
         if(c2 == source.c2 && start1 == source.start1 && start2 == source.start2 
                 && len == source.len )
             return true;
         return false;
    }

    bool operator<(const om& source)const{
         if(c2 < source.c2)
             return true;
         if(c2 == source.c2 && start1 < source.start1)
             return true;
         if(c2 == source.c2 && start1 == source.start1 && start2 < source.start2)
             return true;
         if(c2 == source.c2 && start1 == source.start1 && start2==source.start2 
                 && len < source.len)
             return true;
         return false;
     }

     friend Bool operator==( const om& o1, const om& o2 )
     {    return o1.c2 == o2.c2 && o1.start1 == o2.start1 && o1.start2 == o2.start2 
               && o1.len == o2.len;    }

};

TRIVIALLY_SERIALIZABLE(om<uint16_t>);
TRIVIALLY_SERIALIZABLE(om<int>);

class olap {

     public:

     olap( ) { }
     olap( const int start1, const int start2, const int len ) :
          start1(start1), start2(start2), len(len) { }

     int start1, start2;
     int len;

     int Offset( ) const { return start1 - start2; }

     int Start1( ) const { return start1; }
     int Start2( ) const { return start2; }
     int Stop1( ) const { return start1 + len; }
     int Stop2( ) const { return start2 + len; }
     int Len() const {return len;}

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
               {    cout << "\nInvalid olap." << endl;
                    TracebackThisProcess( );    }    }    }

     friend Bool operator==( const olap& o1, const olap& o2 )
     {    return o1.start1 == o2.start1 && o1.start2 == o2.start2 
               && o1.len == o2.len;    }

     bool operator==(const olap& source){
         if(this==&source)
             return true;
         if(start1 == source.start1 && start2 == source.start2 && len == source.len )
             return true;
         return false;
     }

     bool operator<(const olap& source)const{
         if(start1 < source.start1)
             return true;
         if(start1 == source.start1 && start2 < source.start2)
             return true;
         if(start1 == source.start1 && start2==source.start2 && len < source.len)
             return true;
         return false;
     }
};

TRIVIALLY_SERIALIZABLE(olap);

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

    template<>
    struct hash<qtuple>{
        public:
        std::size_t operator()(const qtuple & q) const{
            std::size_t seed = 4;
            seed ^= std::hash<int>()(q.d1) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            seed ^= std::hash<int>()(q.a1) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            seed ^= std::hash<int>()(q.d2) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            seed ^= std::hash<int>()(q.a2) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            return seed;
        }
    };
    template<>
    struct hash<pair<int64_t,int64_t>>{
        public:
        std::size_t operator()(const pair<int64_t,int64_t> & q) const{
            std::size_t seed = 4;
            seed ^= std::hash<int64_t>()(q.first) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            seed ^= std::hash<int64_t>()(q.second) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            return seed;
        }
    };
    template<>
    struct hash<pair<int,int>>{
        public:
        std::size_t operator()(const pair<int,int> & q) const{
            std::size_t seed = 4;
            seed ^= std::hash<int>()(q.first) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            seed ^= std::hash<int>()(q.second) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            return seed;
            /* return (((std::size_t)q.first+q.second)*(q.first+q.second+1)/2 + q.second); */
        }
    };
    template<>
    struct hash<pair<uchar,int>>{
        public:
        std::size_t operator()(const pair<uchar,int> & q) const{
            std::size_t seed = 4;
            seed ^= std::hash<uchar>()(q.first) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            seed ^= std::hash<int>()(q.second) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            return seed;
            /* return (((std::size_t)q.first+q.second)*(q.first+q.second+1)/2 + q.second); */
        }
    };
};

#endif
