///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "CoreTools.h"
#include "PackAlign.h"
#include "lookup/LookAlign.h"
#include "pairwise_aligners/ClusterAligner.h"
#include "pairwise_aligners/SmithWatBandedA.h"
#include "paths/LongReadTools.h"

class kmer_cluster {

     public:

     kmer_cluster( ) { }
     kmer_cluster( const int n, const int t, const int offset_low,
          const int offset_high, const Bool fw ) : n(n), t(t), 
          offset_low(offset_low), offset_high(offset_high), fw(fw) { }

     int n;                        // hits in cluster
     int t;                        // target contig
     int offset_low, offset_high;  // offset range
     Bool fw;                      // orientation

     friend Bool operator<( const kmer_cluster& c1, const kmer_cluster& c2 )
     {    return c1.n < c2.n;    }

};

void ClusterAligner( basevector q, const vecbasevector& G, const int K, 
     const VecIntPairVec& Glocs, vec<look_align>& aligns, const Bool FW_ONLY,
     const int BW_ADD, const int MIN_CLUSTER, const int MAX_OFFSET_DIFF,
     const int MISMATCH_PENALTY, const int GAP_PENALTY )
{
     // Heuristics.

     const double WINNING_EDGE = 2.5;

     // Align.

     aligns.clear( );
     align a;
     vec<int> errs;
     int npasses = ( FW_ONLY ? 1 : 2 );
     for ( int pass = 1; pass <= npasses; pass++ )
     {    if ( pass == 2 ) q.ReverseComplement( );
          vec< triple<int,int,int> > locs;
          vec<kmer_cluster> clusters;
          int N = 0;
          for ( int j = 0; j <= q.isize( ) - K; j++ )
          {    int n = KmerId( q, K, j );
               N += Glocs[n].size( );    }
          locs.reserve(N);
          for ( int j = 0; j <= q.isize( ) - K; j++ )
          {    int n = KmerId( q, K, j );
               if ( Glocs[n].size( ) > 1000 ) continue;
               for ( unsigned z = 0; z < Glocs[n].size( ); z++ )
               {    locs.push( Glocs[n][z].first, 
                         Glocs[n][z].second - j, j );    }    }
          Sort(locs);
          for ( int i = 0; i < locs.isize( ); i++ )
          {    int j;
               for ( j = i + 1; j < locs.isize( ); j++ )
               {    if ( locs[j].first != locs[i].first ) break;
                    if ( locs[j].second - locs[i].second > MAX_OFFSET_DIFF )
                         break;    }
               vec<Bool> hit( q.isize( ), False );
               for ( int k = i; k < j; k++ )
                    hit[ locs[k].third ] = True;
               clusters.push( Sum(hit), locs[i].first, 
                    locs[i].second, locs[j-1].second, pass == 1 );
               i = j - 1;    }
          if ( clusters.nonempty( ) )
          {    ReverseSort(clusters);
               int high = clusters[0].n;
               int low = int( floor( high - WINNING_EDGE * sqrt(double(high)) ) );
               low = Max( low, MIN_CLUSTER );
               for ( int j = 0; j < clusters.isize( ); j++ )
               {    if ( clusters[j].n < low ) break;
                    int offset_low = clusters[j].offset_low;
                    int offset_high = clusters[j].offset_high;
                    int offset = (offset_low + offset_high) / 2;
                    int bandwidth = (offset_high - offset_low) / 2;
                    int errsx;
                    int g = clusters[j].t;
                    SmithWatBandedA( q, G[g], -offset, bandwidth + BW_ADD, a, errsx,
                         0, MISMATCH_PENALTY, GAP_PENALTY );
                    vector<int> mgg;
                    mgg = a.MutationsGap1Gap2( q, G[g] );
                    int mutations = mgg[0], indels = mgg[1] + mgg[2];
                    // PRINT7( clusters[j].n, pass, g, offset, bandwidth, mutations, indels ); // XXXX
                    aligns.push( 0, g, q.size( ), G[g].size( ), !clusters[j].fw,
                         a, clusters[j].n, mutations, indels );
                    errs.push_back( mutations + indels );    }    }    }
     if ( aligns.empty( ) ) return;
     SortSync( errs, aligns );
     int n;
     for ( n = 0; n < errs.isize( ); n++ )
          if ( errs[n] > errs[0] + WINNING_EDGE * sqrt(double(errs[0])) ) break;
     aligns.resize(n);    }
