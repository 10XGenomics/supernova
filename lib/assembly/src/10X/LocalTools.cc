// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"
#include "ParallelVecUtilities.h"
#include "paths/HyperBasevector.h"
#include "10X/LocalTools.h"
#include "10X/Super.h"

void GlueAssemblies( const HyperBasevectorX& hb, digraphE<vec<int>>& D, 
     vec<int>& dinv, const vec<int> mult,
     const int MIN_LINE_TO_WALK, const int MAX_MULT, const Bool verbose )
{
     // Go through multiple join steps.

     for ( int pass = 1; ; pass++ )
     {
          // Gather long edges in the graph.

          if (verbose) cout << "\n" << Date( ) << ": pass = " << pass << endl;
          if (verbose) cout << Date( ) << ": getting edge lengths" << endl;
          vec<int> dlens( D.E( ), 0 ), dmm( D.E( ), 1000000000 );
          #pragma omp parallel for
          for ( int d = 0; d < D.E( ); d++ )
          {    if ( D.O(d)[0] < 0 ) continue;
               for ( int j = 0; j < D.O(d).isize( ); j++ )
               {    dlens[d] += hb.Kmers( D.O(d)[j] );
                    dmm[d] = Min( dmm[d], mult[ D.O(d)[j] ] );    }    }
          if (verbose) cout << Date( ) << ": processing long edges" << endl;
          vec<triple<int,int,int>> xlocs;
          for ( int d = 0; d < D.E( ); d++ )
          {    if ( dlens[d] < MIN_LINE_TO_WALK ) continue;
               if ( dmm[d] > MAX_MULT ) continue;
               for ( int j = 0; j < D.O(d).isize( ); j++ )
               {    int e = D.O(d)[j];
                    if ( mult[e] > 1 ) xlocs.push( e, d, j );    }    }
          if (verbose) cout << Date( ) << ": sorting" << endl;
          ParallelSort(xlocs);
          vec<int64_t> xstarts( hb.E( ) + 1, -1 );
          xstarts[0] = 0;
          for ( int64_t j = xlocs.jsize( ) - 1; j >= 0; j-- )
               xstarts[ xlocs[j].first ] = j;
          xstarts[ hb.E( ) ] = xlocs.size( );
          for ( int64_t j = hb.E( ) - 1; j >= 0; j-- )
               if ( xstarts[j] < 0 ) xstarts[j] = xstarts[j+1];

          // Find long matches.  This finds instances of two edges (d1, d2) that 
          // share an identical stretch of >= MIN_LINE_TO_WALK kmers.  Once we find 
          // an instance for an edge d1, we move on to the next d1.

          if (verbose) cout << Date( ) << ": looking for long matches" << endl;
          vec< triple< pair<int,int>, pair<int,int>, int > > M;
          #pragma omp parallel for schedule(dynamic,10000)
          for ( int d1 = 0; d1 < D.E( ); d1++ )
          {    if ( dlens[d1] < MIN_LINE_TO_WALK ) continue;
               if ( dmm[d1] > MAX_MULT ) continue;
               for ( int j1 = 0; j1 < D.O(d1).isize( ); j1++ )
               {    int e = D.O(d1)[j1], jx, s0 = 0;
                    for ( jx = j1; jx < D.O(d1).isize( ); jx++ )
                    {    s0 += hb.Kmers( D.O(d1)[jx] );
                         if ( s0 >= MIN_LINE_TO_WALK ) break;    }
                    if ( s0 < MIN_LINE_TO_WALK ) break;
                    int64_t low = xstarts[e], high = xstarts[e+1];
                    for ( int64_t i = low; i < high; i++ )
                    {    int d2 = xlocs[i].second, j2 = xlocs[i].third;
                         if ( d2 == d1 && j2 == j1 ) continue;
                         if ( j1 > 0 && j2 > 0 && D.O(d1)[j1-1] == D.O(d2)[j2-1] ) 
                              continue;
                         if ( j2 - j1 + jx >= D.O(d2).isize( ) ) continue;
                         int rd1 = dinv[d1], rd2 = dinv[d2];
                         if ( !IsUnique( vec<int>{ d1, d2, rd1, rd2 } ) ) continue;
                         Bool fail = False;
                         for ( int l1 = jx; l1 > j1; l1-- )
                         {    int l2 = j2 - j1 + l1;
                              if ( D.O(d1)[l1] != D.O(d2)[l2] )
                              {    fail = True;
                                   break;    }    }
                         if (fail) continue;
                         int s = 0, l1;
                         for ( l1 = jx + 1; l1 < D.O(d1).isize( ); l1++ )
                         {    int l2 = j2 - j1 + l1;
                              if ( l2 >= D.O(d2).isize( ) ) break;
                              if ( D.O(d1)[l1] != D.O(d2)[l2] ) break;
                              s += hb.Kmers( D.O(d1)[l1] );    }
                         int min_mult = 1000000000;
                         for ( int m = j1; m < l1; m++ )
                              min_mult = Min( min_mult, mult[ D.O(d1)[m] ] );
                         if ( min_mult > MAX_MULT ) continue;
                         #pragma omp critical
                         {    M.push( 
                                   make_pair(d1,j1), make_pair(d2,j2), l1-j1 );    }
                         goto next_d1;    }    }
                    next_d1: continue;    }
          if (verbose) cout << Date( ) << ": sorting matches" << endl;
          ParallelSort(M);

          // Make merges.

          if ( MakeMerges( M, D, dinv, verbose ) == 0 ) break;    }    }
