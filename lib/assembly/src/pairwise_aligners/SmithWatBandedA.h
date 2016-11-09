///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef SMITHWATBANDEDA
#define SMITHWATBANDEDA

#include "Alignment.h"
#include "Basevector.h"

// Warning.  An object of type X contains the error count at a given position in
// the dynamic programming matrix.  Below, we use X = unsigned char, which could
// easily overflow.  Use of a larger type is recommended but carries with it a
// performance penalty.

// Note also the extremely ugly inclusion of the dummy argument x, which is needed,
// for unknown reasons.

template<class X> float SmithWatBandedA2( const basevector& S, 
     const basevector& T, int offset, int bandwidth, align& a, int& errors, 
     ostream* log = 0, int mis=2, int ins=3, int del=3 );

inline float SmithWatBandedA( const basevector& S, const basevector& T, 
     int offset, int bandwidth, align& a, int& errors, ostream* log = 0, 
     int mis=2, int gap=3  )
{    
     return SmithWatBandedA2<unsigned char>( 
          S, T, offset, bandwidth, a, errors, log, mis, gap, gap  );    }

// invoking run() will perform the same function as SmithWatBandedA, as of r48232
// the point is that each instance retains its storage and thus memory-bound behavior
// is lessened, if each thread uses its own instance.
class SmithWatBandedAEngine{
public:
    SmithWatBandedAEngine(size_t a, size_t b): m_from(a,2*b+4),m_x(a+1),m_s(a+1) { };
    float run( const basevector& S, const basevector& T
             , int offset, int bandwidth, align& a, int& errors
             , ostream* log = 0, int mis=2, int gap=3  ){
        return run2A( S, T, offset, bandwidth, a, errors, log, mis, gap, gap  );
    };
private:
    float run2A( const basevector& S, const basevector& T
               , int offset, int bandwidth, align& a, int& errors
               , ostream* log = 0, int mis=2, int ins=3, int del=3 );
    class matrix_t{
    public:
        matrix_t():storage(),n1(0),n2(0){};
        matrix_t(size_t a,size_t b):storage(a*b,0),n1(a),n2(b){};
        unsigned char& operator()(size_t a, size_t b){return storage[a*n2+b];};
        const unsigned char& operator()(size_t a, size_t b)const {return storage[a*n2+b];};
        void reset(size_t a, size_t b, unsigned char c){
            storage.clear();
            n1=a;
            n2=b;
            storage.resize(n1*n2,c);
        }
    private:
        vec<unsigned char> storage;
        size_t n1,n2;
    };
    matrix_t m_from;
    vec<int> m_x;
    vec<char> m_s;
};

inline float SmithWatBandedA( const basevector& S, const basevector& T, 
     int offset, int bandwidth, alignment& a, ostream *log = 0 )
{    align temp; int errors;
     float result = SmithWatBandedA( S, T, offset, bandwidth, temp, errors, log );
     a.Set( temp, errors );
     return result;    }

#endif
