// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

// Patch in original gap closures.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"
#include "ParallelVecUtilities.h"
#include "kmers/MakeLookup.h"
#include "paths/HyperBasevector.h"
#include "paths/RemodelGapTools.h"
#include "10X/Capture.h"
#include "10X/Gap.h"
#include "10X/Splat.h"
#include "10X/Super.h"

void Splat( const HyperBasevectorX& hb, const vec<int>& inv,
     const vecbasevector& closures, digraphE<vec<int>>& D, vec<int>& dinv )
{
     // Check K.

     const int K = 48;
     ForceAssertEq( K, hb.K( ) );

     // Path the closures.

     vec<vec<int>> cpaths;
     {    cout << Date( ) << ": making lookup table, mem usage = " 
               << MemUsageGBString( ) << endl;
          vec< triple<kmer<K>,int,int> > kmers_plus;
          // MakeKmerLookup0( hb.Edges( ), kmers_plus );
          MakeKmerLookup0x( hb.Edges( ), kmers_plus ); // slower but less mem
          cout << Date( ) << ": pathing " << ToStringAddCommas( closures.size( ) )
               << " closures" << endl;
          cpaths.resize( 2 * closures.size( ) );
          #pragma omp parallel for schedule( dynamic, 1000 )
          for ( int c = 0; c < (int) closures.size( ); c++ )
          {    const basevector& C = closures[c];
               if ( C.isize( ) < K ) continue;
               kmer<K> x;
               Bool fail = False;
               vec< pair<int,int> > X( C.isize( ) - K + 1 );
               for ( int j = 0; j <= C.isize( ) - K; j++ )
               {    x.SetToSubOf( C, j );
                    int64_t p = BinPosition1( kmers_plus, x );
                    if ( p < 0 )
                    {    fail = True;
                         break;    }
                    X[j] = make_pair( 
                         kmers_plus[p].second, kmers_plus[p].third );    }
               if (fail) continue;
               vec<int> es = { X[0].first };
               for ( int j = 1; j < X.isize( ); j++ )
               {    int e = X[j].first, ep = X[j-1].first;
                    if ( e == ep && X[j].second == X[j-1].second + 1 );
                    else if ( X[j].second == 0 && X[j-1].second == hb.Kmers(ep) - 1
                         && hb.ToRight(ep) == hb.ToLeft(e) )
                    {    es.push_back(e);    }
                    else
                    {    fail = True;
                         break;   }   }
               if (fail) continue;
               cpaths[2*c] = es;
               es.ReverseMe( );
               for ( int i = 0; i < es.isize( ); i++ )
                    es[i] = inv[ es[i] ];
               cpaths[2*c+1] = es;    }    }
     cout << Date( ) << ": sorting, mem usage = " << MemUsageGBString( ) << endl;
     ParallelUniqueSort(cpaths);
     cout << Date( ) << ": have " << ToStringAddCommas( cpaths.size( ) ) 
          << " closure paths" << endl;

     // Index the closure paths.

     cout << Date( ) << ": indexing closure paths" << endl;
     vec<vec<pair<int,int>>> pos( hb.E( ) );
     for ( int i = 0; i < cpaths.isize( ); i++ )
     for ( int j = 0; j < cpaths[i].isize( ); j++ )
          pos[ cpaths[i][j] ].push( i, j );

     // Find gap closures.

     cout << Date( ) << ": finding gap closures" << endl;
     const int MAX_BACK = 100;
     const int MAX_PATHS = 4;
     vec< quad< int, vec<int>, vec<int>, vec<vec<int>> > > edits;
     vec<int> to_left,to_right;
     D.ToLeft(to_left), D.ToRight(to_right);
     #pragma omp parallel for schedule( dynamic, 1000 )
     for ( int d = 0; d < D.E( ); d++ )
     {    if ( !IsPairGap( D.O(d) ) || dinv[d] < d ) continue;
          int v = to_left[d], w = to_right[d];
          if ( !D.From(v).solo( ) || !D.To(v).solo( ) ) continue;
          if ( !D.From(w).solo( ) || !D.To(w).solo( ) ) continue;
          int d1 = D.ITo( v, 0 ), d2 = D.IFrom( w, 0 );
          const vec<int> &x1 = D.O(d1), &x2 = D.O(d2);
          if ( x1[0] < 0 || x2[0] < 0 ) continue;

          // Find the number of edges on the left and right that add up to
          // MAX_BACK kmers.

          vec<int> play(2, 0);
          for ( int j = 0; j < 2; j++ )
          {    int n = 0;
               if ( j == 0 )
               {    for ( int i = x1.isize( ) - 1; i >= 0; i-- )
                    {    n += hb.Kmers( x1[i] );
                         play[j]++;
                         if ( n >= MAX_BACK ) break;    }    }
               else
               {    for ( int i = 0; i < x2.isize( ); i++ )
                    {    n += hb.Kmers( x2[i] );
                         play[j]++;
                         if ( n >= MAX_BACK ) break;    }    }    }

          // Look for patches.

          vec<vec<int>> Z;
          for ( int i1 = x1.isize( ) - play[0]; i1 < x1.isize( ); i1++ )
          {    int e1 = x1[i1];
               for ( int j1 = 0; j1 < pos[e1].isize( ); j1++ )
               {    int p1 = pos[e1][j1].first, k1 = pos[e1][j1].second;
                    if ( i1 > x1.isize( ) - play[0] && k1 > 0 ) continue;
                    for ( int i2 = 0; i2 < play[1]; i2++ )
                    {    int e2 = x2[i2];
                         for ( int j2 = 0; j2 < pos[e2].isize( ); j2++ )
                         {    int p2 = pos[e2][j2].first, k2 = pos[e2][j2].second;
                              if ( p1 != p2 || k1 > k2 ) continue;
                              int p = p1;
                              if ( i2 < play[1] - 1 && k2 < cpaths[p].isize( ) - 1 )
                                   continue;
     
                              // Now cpaths[p][k1] = x1[i1]
                              // and cpaths[p][k2] = x2[i2].
     
                              vec<int> z;
                              for ( int m = x1.isize( ) - play[0]; m < i1; m++ )
                                   z.push_back( x1[m] );
                              for ( int m = k1; m <= k2; m++ )
                                   z.push_back( cpaths[p][m] );
                              for ( int m = i2 + 1; m < play[1]; m++ )
                                   z.push_back( x2[m] );
                              Z.push_back(z);    }    }    }    }
          UniqueSort(Z);
          if ( Z.empty( ) || Z.isize( ) > MAX_PATHS ) continue;
          vec<int> left, right;
          for ( int i = x1.isize( ) - play[0]; i < x1.isize( ); i++ )
               left.push_back( x1[i] );
          for ( int i = 0; i < play[1]; i++ )
               right.push_back( x2[i] );
          #pragma omp critical
          {    edits.push( d, left, right, Z );    }    }
     Sort(edits);
     cout << Date( ) << ": found " << ToStringAddCommas( edits.size( ) ) 
          << " edits" << endl;

     // Edit the graph.

     cout << Date( ) << ": editing the graph" << endl;
     vec<int> dels;
     for ( int i = 0; i < edits.isize( ); i++ )
     {    int d = edits[i].first; 
          vec<int> &left = edits[i].second, &right = edits[i].third;
          vec<vec<int>>& Z = edits[i].fourth;
          int v = to_left[d], w = to_right[d];
          int d1 = D.ITo( v, 0 ), d2 = D.IFrom( w, 0 );
          int rd = dinv[d];
          int rv = to_left[rd], rw = to_right[rd];
          int rd2 = D.ITo( rv, 0 ), rd1 = D.IFrom( rw, 0 );
          if ( !IsUnique( vec<int>{ d1, d2, rd1, rd2 } ) ) continue;
          int n1 = D.O(d1).size( ), n2 = D.O(d2).size( );
          if ( left.isize( ) > n1 ) continue;
          if ( !D.O(d1).Contains( left, n1 - left.isize( ) ) ) continue;
          if ( !D.O(d2).Contains( right, 0 ) ) continue;
          D.OMutable(d1).resize( n1 - left.size( ) );
          D.OMutable(d2).SetToSubOf( D.O(d2), right.size( ), n2 - right.size( ) );
          D.OMutable(rd1) = D.O(d1), D.OMutable(rd2) = D.O(d2);
          D.OMutable(rd1).ReverseMe( ), D.OMutable(rd2).ReverseMe( );
          for ( int j = 0; j < D.O(rd1).isize( ); j++ )
               D.OMutable(rd1)[j] = inv[ D.O(rd1)[j] ];
          for ( int j = 0; j < D.O(rd2).isize( ); j++ )
               D.OMutable(rd2)[j] = inv[ D.O(rd2)[j] ];
          for ( int j = 0; j < Z.isize( ); j++ )
          {    int N = D.E( );
               D.AddEdge( v, w, Z[j] );
               to_left.push_back(v), to_right.push_back(w);
               Z[j].ReverseMe( );
               for ( int l = 0; l < Z[j].isize( ); l++ )
                    Z[j][l] = inv[ Z[j][l] ];
               D.AddEdge( rv, rw, Z[j] );
               to_left.push_back(rv), to_right.push_back(rw);
               dinv.push_back( N+1, N );    }
          dels.push_back( d, dinv[d] );    }

     // Remove empty edges just now created.

     cout << Date( ) << ": removing empty edges" << endl;
     for ( int d = 0; d < D.E( ); d++ )
     {    if ( D.O(d).empty( ) )
          {    int v = to_left[d], w = to_right[d];
               dels.push_back(d);
               if ( v != w )
                    D.TransferEdgesWithUpdate( v, w, to_left, to_right );    }    }
     D.DeleteEdges(dels);
     CleanupCore( D, dinv );

     // Clean up.

     cout << Date( ) << ": cleaning up after patching" << endl;
     const int MAX_KILL = 100;
     const double MIN_RATIO = 10;
     Bool verbose = True;
     Bool single = False;
     dels.clear( );
     SimpleHangs( hb, D, dinv, dels, MAX_KILL, MIN_RATIO, verbose, single );
     D.DeleteEdges(dels);
     Zipper( D, dinv );
     dels.clear( );
     CaptureMessyLoops( hb, inv, D, dinv, dels, True );
     D.DeleteEdges(dels);
     RemoveUnneededVertices( D, dinv );
     CleanupCore( D, dinv );
     Validate( hb, inv, D, dinv );    }
