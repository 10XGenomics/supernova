// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"
#include "graph/DigraphTemplate.h"
#include "paths/long/ReadPath.h"
#include "paths/long/large/Lines.h"
#include "10X/Heuristics.h"
#include "10X/IntIndex.h"
#include "10X/PullApart.h"

void PullApartInversions( digraphE<vec<int>>& D, vec<int>& dinv )
{    cout << Date( ) << ": pulling apart inversions" << endl;
     vec<vec<vec<vec<int>>>> dlines;
     FindLines( D, dinv, dlines, MAX_CELL_PATHS, MAX_CELL_DEPTH );
     vec<int> to_left, to_right, linv;
     D.ToLeft(to_left), D.ToRight(to_right);
     LineInv( dlines, dinv, linv );
     vec<int> tol( D.E( ), -1 );
     for ( int i = 0; i < dlines.isize( ); i++ )
     for ( int j = 0; j < dlines[i].isize( ); j++ )
     for ( int k = 0; k < dlines[i][j].isize( ); k++ )
     for ( int l = 0; l < dlines[i][j][k].isize( ); l++ )
          tol[ dlines[i][j][k][l] ] = i;
     vec<int> mods;
     #pragma omp parallel for schedule(dynamic, 10000)
     for ( int m = 0; m < dlines.isize( ); m++ )
     {    if ( linv[m] != m ) continue;
          int d1 = dlines[m][0][0][0];
          int v1 = to_left[d1];
          if ( D.To(v1).size( ) != 2 ) continue;
          #pragma omp critical
          {    mods.push_back(m);    }    }
     Sort(mods);
     int count = 0;
     for ( int mi = 0; mi < mods.isize( ); mi++ )
     {    int m = mods[mi];
          int d1 = dlines[m].front( )[0][0], d2 = dlines[m].back( )[0][0];
          int v1 = to_left[d1], v2 = to_right[d2];
          int f1 = D.ITo(v1,0), f2 = D.ITo(v1,1);
          int g1 = dinv[f1], g2 = dinv[f2];
          int l1 = tol[f1], l2 = tol[f2];
          int r1 = linv[l1], r2 = linv[l2];
          vec<int> es;
          for ( int j = 0; j < dlines[m].isize( ); j++ )
          for ( int k = 0; k < dlines[m][j].isize( ); k++ )
          for ( int l = 0; l < dlines[m][j][k].isize( ); l++ )
               es.push_back( dlines[m][j][k][l] );
          UniqueSort(es);
          count++;
          digraphE<vec<int>> M( digraphE<vec<int>>::COMPLETE_SUBGRAPH_EDGES,
               D, es, to_left, to_right );
          vec<int> mto_left, mto_right;
          M.ToLeft(mto_left), M.ToRight(mto_right);
          int ml = mto_left[ BinPosition( es, d1 ) ];
          int mr = mto_right[ BinPosition( es, d2 ) ];
          int N = D.N( ), E = D.E( );
          for ( int d = 0; d < M.E( ); d++ )
          {    mto_left[d] += N, mto_right[d] += N;    }

          for ( int i = 0; i < es.isize( ); i++ )
          {    int e = es[i];
               dinv.push_back( dinv[e] );
               dinv[e] = E + BinPosition( es, dinv[e] );    }
          D.Append(M);
          to_left.append(mto_left), to_right.append(mto_right);
          D.GiveEdgeNewToVx( f2, v1, N + ml );
          D.GiveEdgeNewFromVx( g1, v2, N + mr );    }
     cout << Date( ) << ": pulled apart " << count << " inversions" << endl;    }

Bool SupportSplit( const int d1, const int d2, const int f1, const int f2,
     const vec<int>& dinv, const ReadPathVec& dpaths, const IntIndex& dpaths_index )
{
     int sup11 = 0, sup22 = 0, sup12 = 0, sup21 = 0;
     vec<vec<int>> X;
     for ( int pass = 1; pass <= 2; pass++ )
     {    const int d = ( pass == 1 ? d1 : d2 );
          int rd = dinv[d];
          vec<int64_t> used;
          for ( int j = 0; j < dpaths_index.Count(d); j++ )
          {    int64_t id1 = dpaths_index.Val( d, j );
               int64_t id2 = ( id1 % 2 == 0 ? id1+1 : id1-1 );
               used.push_back(id1/2);
               vec<int> x;
               for ( int j = 0; j < (int) dpaths[id1].size( ); j++ )
                    x.push_back( dpaths[id1][j] );
               for ( int j = 0; j < (int) dpaths[id2].size( ); j++ )
                    x.push_back( dinv[ dpaths[id2][j] ] );
               X.push_back(x);    }
          UniqueSort(used);
          for ( int j = 0; j < dpaths_index.Count(rd); j++ )
          {    int64_t id1 = dpaths_index.Val( rd, j );
               if ( BinMember( used, id1/2 ) ) continue;
               vec<int> x;
               for ( int j = 0; j < (int) dpaths[id1].size( ); j++ )
                    x.push_back( dinv[ dpaths[id1][j] ] );
               X.push_back(x);    }    }
     for ( int i = 0; i < X.isize( ); i++ )
     {    const vec<int>& x = X[i];
          int p1 = Position( x, d1 ), p2 = Position( x, d2 );
          if ( p1 >= 0 && p2 >= 0 ) continue;
          if ( p1 >= 0 )
          {    for ( int j = p1 + 1; j < x.isize( ); j++ )
               {    if ( x[j] == f1 )
                    {    sup11++;
                         break;    }    }
               for ( int j = p1 + 1; j < x.isize( ); j++ )
               {    if ( x[j] == f2 )
                    {    sup12++;
                         break;    }    }    }
          if ( p2 >= 0 )
          {    for ( int j = p2 + 1; j < x.isize( ); j++ )
               {    if ( x[j] == f1 )
                    {    sup21++;
                         break;    }    }
               for ( int j = p2 + 1; j < x.isize( ); j++ )
               {    if ( x[j] == f2 )
                    {    sup22++;
                         break;    }    }    }    }
     if ( sup11 < 5 || sup22 < 5 ) return False;
     if ( sup11 + sup22 < 5 * ( sup12 + sup21 ) ) return False;
     if ( sup12 + sup21 > 4 ) return False;
     return True;    }

void PullApart( const vec<int>& inv, digraphE<vec<int>>& D, vec<int>& dinv, 
     const ReadPathVec& dpaths, vec<int>& dels, Bool verbose )
{    // Seems a little dubious that we're not recomputing dpaths.
     IntIndex dpaths_index( dpaths, D.E( ), verbose ); // THIS IS A DUPLICATION!!!!!
     vec<int> to_left, to_right;
     D.ToLeft(to_left), D.ToRight(to_right);
     if (verbose) cout << Date( ) << ": making canonical pullaparts" << endl;
     vec<vec<int>> pulls, pulls2;
     #pragma omp parallel for
     for ( int v = 0; v < D.N( ); v++ )
     {    if ( D.To(v).size( ) != 2 || D.From(v).size( ) != 1 ) continue;
          int w = D.From(v)[0];
          if ( D.To(w).size( ) != 1 || D.From(w).size( ) != 2 ) continue;
          int d1 = D.ITo(v,0), d2 = D.ITo(v,1), e = D.IFrom(v,0);
          for ( int p = 0; p < 2; p++ )
          {    int j1 = p, j2 = 1 - p;
               int f1 = D.IFrom(w,j1), f2 = D.IFrom(w,j2);
               int re = dinv[e];
               if ( !IsUnique( vec<int> { to_left[e], to_right[e],
                    to_left[re], to_right[re] } ) )
               {    continue;    }
               if ( !SupportSplit( d1, d2, f1, f2, dinv, dpaths, dpaths_index ) )
                    continue;
               #pragma omp critical
               {    pulls.push( vec<int>( { d1, d2, e, f1, f2 } ) );    }    
               break;    }    }
     Sort(pulls);
     #pragma omp parallel for
     for ( int v = 0; v < D.N( ); v++ )
     {    if ( D.To(v).size( ) != 2 || D.From(v).size( ) != 2 ) continue;
          int d1 = D.ITo(v,0), d2 = D.ITo(v,1);
          for ( int p = 0; p < 2; p++ )
          {    int j1 = p, j2 = 1 - p;
               int f1 = D.IFrom(v,j1), f2 = D.IFrom(v,j2);
               int rd1 = dinv[d1], rd2 = dinv[d2], rf1 = dinv[f1], rf2 = dinv[f2];
               int rv = to_left[rd1];
               if ( v == rv ) continue;
               if ( !SupportSplit( d1, d2, f1, f2, dinv, dpaths, dpaths_index ) )
                    continue;
               #pragma omp critical
               {    pulls2.push( vec<int>( { d1, d2, f1, f2 } ) );    }    
               break;    }    }
     Sort(pulls2);
     int pullcount = 0;
     vec<Bool> touchedv( D.N( ), False );
     for ( int i = 0; i < pulls.isize( ); i++ )
     {    int d1 = pulls[i][0], d2 = pulls[i][1], e = pulls[i][2];
          int f1 = pulls[i][3], f2 = pulls[i][4];
          int v1 = to_left[d1], w1 = to_right[f1];
          int v2 = to_left[d2], w2 = to_right[f2];
          int rv1 = to_right[dinv[d1]], rw1 = to_left[dinv[f1]];
          int rv2 = to_right[dinv[d2]], rw2 = to_left[dinv[f2]];
          int re = dinv[e], rd2 = dinv[d2], rf2 = dinv[f2];
          int a = to_left[e], b = to_right[e];
          if ( touchedv[a] || touchedv[b] ) continue;
          int ra = to_left[re], rb = to_right[e];
          touchedv[a] = touchedv[b] = True;
          touchedv[ra] = touchedv[rb] = True;
          pullcount++;
          int N = D.N( ), E = D.E( );
          D.AddVertices(4);
          D.AddEdgeWithUpdate( N, N+1, D.O(e), to_left, to_right );
          D.AddEdgeWithUpdate( N+2, N+3, D.O(re), to_left, to_right );
          dinv.push_back( E+1, E );
          D.GiveEdgeNewToVx( d2, to_right[d2], N );
          to_right[d2] = N;
          D.GiveEdgeNewFromVx( f2, to_left[f2], N+1 );
          to_left[f2] = N+1;
          D.GiveEdgeNewToVx( rf2, to_right[rf2], N+2 );
          to_right[rf2] = N+2;
          D.GiveEdgeNewFromVx( rd2, to_left[rd2], N+3 );
          to_left[rd2] = N+3;
          touchedv.push_back( True, True, True, True );    }
     if (verbose) cout << Date( ) << ": made " << pullcount << " pullaparts" << endl;
     int pullcount2 = 0;
     for ( int i = 0; i < pulls2.isize( ); i++ )
     {    int d1 = pulls2[i][0], d2 = pulls2[i][1];
          int f1 = pulls2[i][2], f2 = pulls2[i][3];
          if ( touchedv[ to_right[d2] ] ) continue;
          int rd1 = dinv[d1], rd2 = dinv[d2], rf1 = dinv[f1], rf2 = dinv[f2];
          touchedv[ to_right[d2] ] = True;
          touchedv[ to_left[rd2] ] = True;
          pullcount2++;
          int N = D.N( );
          D.AddVertices(2);
          D.GiveEdgeNewToVx( d2, to_right[d2], N );
          to_right[d2] = N;
          D.GiveEdgeNewFromVx( f2, to_left[f2], N );
          to_left[f2] = N;
          D.GiveEdgeNewFromVx( rd2, to_left[rd2], N+1 );
          to_left[rd2] = N+1;
          D.GiveEdgeNewToVx( rf2, to_right[rf2], N+1 );
          to_right[rf2] = N+1;
          touchedv.push_back( True, True );    }
     if (verbose) 
     {    cout << Date( ) << ": made " << pullcount2 << " pullaparts, type 2" 
               << endl;    }
     D.DeleteEdges(dels);    }
