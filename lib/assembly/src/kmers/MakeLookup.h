///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef MAKE_LOOKUP_H
#define MAKE_LOOKUP_H

#include "Basevector.h"
#include "CoreTools.h"
#include "kmers/KmerRecord.h"
#include "system/SortInPlace.h"

template<int K> void MakeKmerLookup0Pre( const vecbasevector& unibases,
     const String& prefix, vec< triple<kmer<K>,int,int> >& kmers_plus );

template<int K> void MakeKmerLookup2x( const vecbasevector& unibases,
     vec< triple<kmer<K>,int,int> >& kmers_plus )
{    vec<int64_t> starts;
     starts.push_back(0);
     for ( size_t i = 0; i < unibases.size( ); i++ )
     {    const basevector& u = unibases[i];
          starts.push_back( starts.back( ) + Max( 0, u.isize( ) - K + 1 ) );    }
     kmers_plus.resize( starts.back( ) );
     #pragma omp parallel for
     for ( size_t i = 0; i < unibases.size( ); i++ )
     {    const basevector& u = unibases[i];
          kmer<K> x, xrc;
          for ( int j = 0; j <= u.isize( ) - K; j++ )
          {    int64_t r = starts[i] + j;
               x.SetToSubOf( u, j );
               xrc = x;
               xrc.ReverseComplement( );
               Bool fw = ( x < xrc );
               kmers_plus[r].first = ( fw ? x : xrc );
               kmers_plus[r].second = i;
               kmers_plus[r].third = ( fw ? j : -j-1 );    }    }
     sortInPlaceParallel( kmers_plus.begin( ), kmers_plus.end( ) );    }

template<int K> void MakeKmerLookup0x( const vecbasevector& unibases,
     vec< triple<kmer<K>,int,int> >& kmers_plus )
{    vec<int64_t> starts;
     starts.push_back(0);
     for ( size_t i = 0; i < unibases.size( ); i++ )
     {    const basevector& u = unibases[i];
          starts.push_back( starts.back( ) + Max( 0, u.isize( ) - K + 1 ) );    }
     kmers_plus.resize( starts.back( ) );
     #pragma omp parallel for
     for ( size_t i = 0; i < unibases.size( ); i++ )
     {    const basevector& u = unibases[i];
          kmer<K> x;
          for ( int j = 0; j <= u.isize( ) - K; j++ )
          {    int64_t r = starts[i] + j;
               x.SetToSubOf( u, j ); 
               kmers_plus[r].first = x;
               kmers_plus[r].second = i; 
               kmers_plus[r].third = j;    }    }
     cout << Date( ) << ": sorting kmers" << endl;
     sortInPlaceParallel( kmers_plus.begin( ), kmers_plus.end( ) );    }

#endif
