// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

// BaseFinLookup: find all perfect matches of finished sequence to base graph edges.
// Output: X[g] = sorted ( (interval on finished sequence g, e, estart ) )}.

#include <omp.h>

#include "MainTools.h"
#include "ParallelVecUtilities.h"
#include "kmers/KmerRecord.h"
#include "math/HoInterval.h"
#include "paths/HyperBasevector.h"
#include "paths/RemodelGapTools.h"
#include "10X/astats/BaseFinLookup.h"

template<int K> void BaseFinLookup( const HyperBasevectorX& hb, 
     const vecbasevector& G,
     MasterVec< SerfVec< triple< ho_interval, int, int > > >& X )
{
     ForceAssertEq( K, hb.K( ) );
     X.resize( G.size( ) );
     vec< triple<kmer<K>,int,int> > kmers_plus;
     MakeKmerLookup0( G, kmers_plus );
     const int bucket_size = 1000000;
     vec<vec<int>> buckets = { { } };
     int nk = 0;
     for ( int e = 0; e < hb.E( ); e++ )
     {    if ( nk + hb.Kmers(e) <= bucket_size )
          {    buckets.back( ).push_back(e);
               nk += hb.Kmers(e);    }
          else
          {    buckets.push_back( {e} );
               nk = hb.Kmers(e);    }   }
     vec< quad<int,int,int,int> > locs; // { (g,gpos,e,epos) }
     int nthreads = omp_get_max_threads( ), cur = 0;
     #pragma omp parallel for
     for ( int t = 0; t < nthreads; t++ )
     {    vec< quad<int,int,int,int> > locst;
          kmer<K> y;
          while(1)
          {    int b;
               Bool done = False;
               #pragma omp critical
               {    if ( cur == buckets.isize( ) ) done = True;
                    else b = cur++;    }
               if (done) break;
               const vec<int>& x = buckets[b];
               for ( int j = 0; j < x.isize( ); j++ )
               {    int e = x[j];
                    for ( int l = 0; l <= hb.O(e).isize( ) - K; l++ )
                    {    y.SetToSubOf( hb.O(e), l );
                         int64_t low = LowerBound1( kmers_plus, y );
                         for ( int64_t m = low; m < kmers_plus.jsize( ); m++ )
                         {    if ( kmers_plus[m].first != y ) break;
                              locst.push( kmers_plus[m].second, e,
                                   kmers_plus[m].third, l );    
                                   }    }    }    }
          #pragma omp critical
          {    locs.append(locst);    }    }
     ParallelSort(locs);
     for ( int64_t i = 0; i < locs.jsize( ); i++ )
     {    int64_t j;
          for ( j = i + 1; j < locs.jsize( ); j++ )
          {    if ( locs[j].first != locs[i].first ) break;
               if ( locs[j].third != locs[j-1].third + 1 ) break;
               if ( locs[j].second != locs[j-1].second ) break;
               if ( locs[j].fourth != locs[j-1].fourth + 1 ) break;    }
          X[ locs[i].first ].push_back( make_triple(
               ho_interval( locs[i].third, locs[j-1].third + K ),
               locs[i].second, locs[i].fourth ) );
          i = j - 1;    }    }





template<int K> void BaseFinLookupSup( const HyperBasevectorX& hb, 
     const vecbasevector& tigs, const vecbasevector& G,
     MasterVec< SerfVec< triple< ho_interval, int, int > > >& X )
{
     ForceAssertEq( K, hb.K( ) );
     X.resize( G.size( ) );
     vec< triple<kmer<K>,int,int> > kmers_plus;
     MakeKmerLookup0( G, kmers_plus );
     const int bucket_size = 1000000;
     vec<vec<int>> buckets = { { } };
     int nk = 0;
     for ( int e = hb.E( ); e < (int) tigs.size( ); e++ )
     {    if ( nk + tigs[e].isize( ) - K + 1 <= bucket_size )
          {    buckets.back( ).push_back(e);
               nk += tigs[e].isize( ) - K + 1;    }
          else
          {    buckets.push_back( {e} );
               nk = tigs[e].isize( ) - K + 1;    }   }
     vec< quad<int,int,int,int> > locs; // { (g,gpos,e,epos) }
     int nthreads = omp_get_max_threads( ), cur = 0;
     #pragma omp parallel for
     for ( int t = 0; t < nthreads; t++ )
     {    vec< quad<int,int,int,int> > locst;
          kmer<K> y;
          while(1)
          {    int b;
               Bool done = False;
               #pragma omp critical
               {    if ( cur == buckets.isize( ) ) done = True;
                    else b = cur++;    }
               if (done) break;
               const vec<int>& x = buckets[b];
               for ( int j = 0; j < x.isize( ); j++ )
               {    int e = x[j];
                    for ( int l = 0; l <= tigs[e].isize( ) - K; l++ )
                    {    y.SetToSubOf( tigs[e], l );
                         int64_t low = LowerBound1( kmers_plus, y );
                         for ( int64_t m = low; m < kmers_plus.jsize( ); m++ )
                         {    if ( kmers_plus[m].first != y ) break;
                              locst.push( kmers_plus[m].second, e,
                                   kmers_plus[m].third, l );    
                                   }    }    }    }
          #pragma omp critical
          {    locs.append(locst);    }    }
     ParallelSort(locs);
     for ( int64_t i = 0; i < locs.jsize( ); i++ )
     {    int64_t j;
          for ( j = i + 1; j < locs.jsize( ); j++ )
          {    if ( locs[j].first != locs[i].first ) break;
               if ( locs[j].third != locs[j-1].third + 1 ) break;
               if ( locs[j].second != locs[j-1].second ) break;
               if ( locs[j].fourth != locs[j-1].fourth + 1 ) break;    }
          X[ locs[i].first ].push_back( make_triple(
               ho_interval( locs[i].third, locs[j-1].third + K ),
               locs[i].second, locs[i].fourth ) );
          i = j - 1;    }    }












typedef triple<ho_interval,int,int> tripto;
extern template class SmallVec<tripto,MempoolAllocator<tripto> >;
extern template class OuterVec< SerfVec<tripto> >;
#include "feudal/SmallVecDefs.h"
#include "feudal/OuterVecDefs.h"
template class SmallVec<tripto,MempoolAllocator<tripto> >;
template class OuterVec< SerfVec<tripto> >;

template void BaseFinLookup<48>( const HyperBasevectorX& hb,
     const vecbasevector& G,
     MasterVec< SerfVec< triple< ho_interval, int, int > > >& X );

template void BaseFinLookupSup<48>( const HyperBasevectorX& hb,
     const vecbasevector& tigs, const vecbasevector& G,
     MasterVec< SerfVec< triple< ho_interval, int, int > > >& X );
