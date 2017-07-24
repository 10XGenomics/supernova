// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>

#include "CoreTools.h"
#include "Equiv.h"
#include "math/HoInterval.h"
#include "paths/HyperBasevector.h"
#include "10X/astats/BaseFinLookup.h"
#include "10X/astats/ProjectRunway.h"

// Find maximal perfect stretches in super graph of at least a given size.  If
// two perfect stretches have the same endpoints on the graph and on the genome,
// we arbitrarily pick one.

// Note that use of flocs isn't quite right in cases where edges have been
// trimmed back because of sequence gap edges.

void GetMaxPerfs( 

     // inputs:

     const vec<int>& gs, const HyperBasevectorX& hb, const vec<int>& kmers,
     const digraphE<vec<int>>& D, const vec<int>& to_left, const vec<int>& to_right,
     const vec<int>& dlens, const vecbasevector& tigs, const vecbasevector& genomef,
     const MasterVec< SerfVec< triple< ho_interval, int, int > > >& flocs0,
     const int minlen, const int minlen_force,

     // output:
     // one entry per element of gs
     // inside, one entry per path
     // path = ( ( interval on finished seq, d, start on d ) )

     vec<vec<vec<triple<ho_interval,int,int>>>>& X )
{
     // Sanity checks on input.

     for ( auto g : gs ) ForceAssertLt( g, flocs0.size( ) );

     const int K = hb.K( );
     const int E = hb.E( );

     // Initialize output.

     X.clear_and_resize( gs.size( ) );

     // Align extra "base" edges.

     MasterVec< SerfVec< triple< ho_interval, int, int > > > flocs = flocs0, flocs2;
     BaseFinLookupSup<48>( hb, tigs, genomef, flocs2 );
     for ( int g = 0; g < (int) genomef.size( ); g++ ) 
     for ( int j = 0; j < (int) flocs2[g].size( ); j++ )
          flocs[g].push_back( flocs2[g][j] );

     // Go through reference tigs.

     int cur = 0, nthreads = omp_get_max_threads( );
     #pragma omp parallel for schedule(dynamic, 1)
     for ( int nt = 0; nt < nthreads; nt++ )
     {    
          // Data structures shared across thread.

          vec<Bool> used( kmers.size( ), False ), mused( kmers.size( ), False );
          vec<Bool> to_delete;
          vec< quad< int, int, int, int > > matches;
          vec<vec<int>> paths;

          // Loop through clones.

          while(1)
          {
               // Get next clone.

               int g = -1, gi = -1;
               #pragma omp critical
               {    if ( cur < gs.isize( ) ) 
                    {    gi = cur;
                         g = gs[cur++];    }    }
               if ( g < 0 ) break;

               // Find the relevant superedges, filtering a bit.
     
               vec<int> x;
               for ( int i = 0; i < (int) flocs[g].size( ); i++ ) 
                    x.push_back( flocs[g][i].second );
               UniqueSort(x);
               vec<int> ds;
               for ( auto e : x ) mused[e] = True;
               for ( int d = 0; d < D.E( ); d++ )
               {    if ( D.O(d)[0] < 0 ) continue;
                    for ( int i = 0; i < D.O(d).isize( ); i++ )
                    {    int e = D.O(d)[i];
                         if ( mused[e] )
                         {    if ( i > 0 && i < D.O(d).isize( ) - 1 )
                              {    int nk = kmers[e];
                                   Bool ok = True;
                                   for ( int j = i - 1; j >= 0; j-- )
                                   {    if ( !mused[ D.O(d)[j] ] ) break;
                                        if ( j == 0 ) ok = False;
                                        nk += kmers[ D.O(d)[j] ];    }
                                   for ( int j = i + 1; j < D.O(d).isize( ); j++ )
                                   {    if ( !mused[ D.O(d)[j] ] ) break;
                                        if ( j == D.O(d).isize( ) - 1 ) ok = False;
                                        nk += kmers[ D.O(d)[j] ];    }
                                   if ( ok && nk + K - 1 < minlen ) continue;    }
                              ds.push_back(d);
                              break;    }    }    }
               for ( auto e : x ) mused[e] = False;
     
               // Index X[g].
     
               vec< pair<int,int> > xin( flocs[g].size( ) );
               for ( int i = 0; i < (int) flocs[g].size( ); i++ )
                    xin[i] = make_pair( flocs[g][i].second, i );
               Sort(xin);
               for ( auto x : xin ) used[ x.first ] = True;
     
               // Find matches.
     
               matches.clear( );
               for ( int i = 0; i < ds.isize( ); i++ )
               {    int d = ds[i], dpos = 0;
                    const vec<int>& x = D.O(d);
                    for ( auto e : x )
                    {    if ( !used[e] ) 
                         {    dpos += kmers[e];
                              continue;    }
                         int low = LowerBound1( xin, e ); 
                         int high = UpperBound1( xin, e );
                         for ( int j = low; j < high; j++ )
                         {    int l = xin[j].second;
                              const ho_interval& h = flocs[g][l].first;
                              int epos = flocs[g][l].third;
                              int dstart = dpos + epos;
                              matches.push( d, dstart - h.Start( ), 
                                   h.Start( ), h.Stop( ) );    }
                         dpos += kmers[e];    }    }
               for ( auto x : xin ) used[ x.first ] = False;
               Sort(matches);
     
               // Combine matches within a single superedge.
     
               to_delete.resize_and_set( matches.size( ), False );
               for ( int j = 1; j < matches.isize( ); j++ )
               {    if ( matches[j].first != matches[j-1].first ) continue;
                    if ( matches[j].second != matches[j-1].second ) continue;
                    if ( matches[j].third != matches[j-1].fourth - K + 1 ) continue;
                    matches[j].third = matches[j-1].third;
                    to_delete[j-1] = True;    }
               EraseIf( matches, to_delete );
               for ( auto& x : matches )
               //                  hstart   d        hstop     dstart
                    x = make_quad( x.third, x.first, x.fourth, x.second+x.third );
               Sort(matches);
     
               // Filter.
     
               to_delete.resize_and_set( matches.size( ), False );
               for ( int j = 0; j < matches.isize( ); j++ )
               {    int d = matches[j].second, dstart = matches[j].fourth;
                    int hstart = matches[j].first, hstop = matches[j].third;
                    if ( dstart == 0 ) continue;
                    if ( dstart + hstop - hstart == dlens[d] + K - 1 ) continue;
                    if ( hstop - hstart < minlen ) to_delete[j] = True;    }
               EraseIf( matches, to_delete );
     
               // Create an acyclic digraph G, whose edges are matches, and for which 
               // edges are connected if they are connected in the super graph and
               // match up on the genome.
     
               int NM = matches.size( );
               vec<int> mi( NM, vec<int>::IDENTITY );
               digraphE<int> H( mi, digraphE<int>::EDGES_SEPARATE );
               equiv_rel e( 2*NM );
               vec< triple<int,int,int> > ends;
               for ( int i = 0; i < matches.isize( ); i++ )
               {    const quad<int,int,int,int>& x = matches[i];
                    int hstart = x.first, hstop = x.third;
                    int d = x.second, dstart = x.fourth;
                    int dstop = dstart + hstop - hstart;
                    int v = to_left[d], w = to_right[d];
                    if ( dstart == 0 ) ends.push( v, hstart, 2*i );
                    if ( dstop == dlens[d] + K - 1 )
                         ends.push( w, hstop-(K-1), 2*i+1 );    }
               Sort(ends);
               for ( int j = 0; j < ends.isize( ); j++ )
               {    int k;
                    for ( k = j + 1; k < ends.isize( ); k++ )
                    {    if ( ends[k].first != ends[j].first ) break;
                         if ( ends[k].second != ends[j].second ) break;    }
                    for ( int l = j; l < k - 1; l++ )
                         e.Join( ends[l].third, ends[l+1].third );
                    j = k - 1;    }
               digraphE<int> G( H, e );
     
               // In G, For each source and each sink to which it is connected,
               // arbitarily pick one path.
     
               paths.clear( );
               for ( int v = 0; v < G.N( ); v++ )
               {    if ( !G.Source(v) ) continue;
                    vec<int> from_v;
                    G.GetSuccessors1( v, from_v );
                    for ( auto w : from_v )
                    {    if ( !G.Sink(w) ) continue;
                         vec<int> to_w, mid;
                         G.GetPredecessors1( w, to_w );
                         mid = Intersection( from_v, to_w );
                         int y = v;
                         vec<int> p;
                         while( y != w )
                         {    for ( int j = 0; j < G.From(y).isize( ); j++ )
                              {    if ( BinMember( mid, G.From(y)[j] ) )
                                   {    p.push_back( G.IFrom(y,j) );
                                        y = G.From(y)[j];
                                        break;    }    }    }
                         paths.push_back(p);    }    }
     
               // Delete short paths.
     
               to_delete.resize_and_set( paths.size( ), False );
               for ( int i = 0; i < paths.isize( ); i++ )
               {    const vec<int>& p = paths[i];
                    int hstart = matches[ G.O( p.front( ) ) ].first;
                    int hstop = matches[ G.O( p.back( ) ) ].third;
                    if ( hstop - hstart < minlen ) to_delete[i] = True;    }
               EraseIf( paths, to_delete );
     
               // Delete dominated paths.
     
               to_delete.resize_and_set( paths.size( ), False );
               vec<ho_interval> cov( paths.size( ) );
               for ( int i = 0; i < paths.isize( ); i++ )
               {    const vec<int>& p = paths[i];
                    int hstart = matches[ G.O( p.front( ) ) ].first;
                    int hstop = matches[ G.O( p.back( ) ) ].third;
                    cov[i] = ho_interval( hstart, hstop );    }
               vec<int> ids( paths.size( ), vec<int>::IDENTITY ); 
               vec<int> stop( cov.size( ) );
               SortSync( cov, ids );
               for ( int i = 0; i < cov.isize( ); i++ )
               {    if ( i == 0 ) stop[i] = cov[0].Stop( );
                    else stop[i] = Max( cov[i].Stop( ), stop[i-1] );    }
               for ( int i = 0; i < cov.isize( ); i++ )
               {    for ( int j = i - 1; j >= 0; j-- )
                    {    if ( stop[j] < stop[i] ) break;
                         if ( cov[j].Length( ) > cov[i].Length( )
                              && cov[i].Length( ) < minlen_force )
                         {    to_delete[ ids[i] ] = True;
                              break;    }    }
                    for ( int j = i + 1; j < cov.isize( ); j++ )
                    {    if ( cov[j].Start( ) > cov[i].Start( ) ) break;
                         if ( cov[j].Length( ) > cov[i].Length( ) 
                              && cov[i].Length( ) < minlen_force )
                         {    to_delete[ ids[i] ] = True;
                              break;    }    }    }
               EraseIf( paths, to_delete );
     
               // Pack up the results.
     
               X[gi].resize( paths.size( ) );
               for ( int i = 0; i < paths.isize( ); i++ )
               {    const vec<int>& p = paths[i];
                    X[gi][i].clear( );
                    for ( auto u : p )
                    {    int j = G.O(u);
                         X[gi][i].push( ho_interval( matches[j].first, 
                              matches[j].third ), matches[j].second, 
                              matches[j].fourth );    }    }    }    }    }
