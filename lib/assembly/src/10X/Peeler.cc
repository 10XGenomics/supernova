// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

// Peel off lines to make direct connections between them.  For each line of 
// length at least 1 kb (MIN_LINE), find read and read pair links forward from its 
// end.  [TO EXPAND]

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"
#include "graph/Digraph.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "paths/long/large/Lines.h"
#include "10X/Gap.h"
#include "10X/Heuristics.h"
#include "10X/MakeLocalsTools.h"
#include "10X/Peeler.h"
#include "10X/Super.h"
#include "10X/paths/ReadPathVecX.h"

void Peeler( const HyperBasevectorX& hb, const vec<int>& inv,
     digraphE<vec<int>>& D, vec<int>& dinv,
     ReadPathVec& dpaths, // should be current on input, not current on exit
     MasterVec<ULongVec>& dpaths_index, // ditto
     ostream& rout // deep logging
     )
{    
     cout << Date( ) << ", entering Peeler, mem = " << MemUsageGBString( ) << endl;

     // Define heuristics.

     const int MIN_LINE = 1000;
     const int MAX_DIST = 500;
     const int MIN_WIN = 5;

     // Find lines.

     vec<vec<vec<vec<int>>>> dlines;
     FindLines( D, dinv, dlines, MAX_CELL_PATHS, MAX_CELL_DEPTH );

     // Compute ancillary info.  Then form lines into a graph.  Copied from 
     // LineLine.cc, should be separate function.  Change: in this versions, the 
     // edges are the line lengths.

     cout << Date( ) << ": computing ancillaries" << endl;
     vec<int> llens, tol, linv, to_left, to_right;
     GetLineLengths( hb, D, dlines, llens );
     GetTol( hb, dlines, tol );
     LineInv( dlines, dinv, linv );
     D.ToLeft(to_left), D.ToRight(to_right);
     vec<int> verts( 2 * dlines.size( ) );
     #pragma omp parallel for
     for ( int i = 0; i < dlines.isize( ); i++ )
     {    int v = to_left[ dlines[i].front( )[0][0] ];
          int w = to_right[ dlines[i].back( )[0][0] ];
          verts[2*i] = v, verts[2*i+1] = w;    }
     UniqueSort(verts);
     int N = verts.size( );
     vec<vec<int>> from(N), to(N), from_edge_obj(N), to_edge_obj(N);
     for ( int i = 0; i < dlines.isize( ); i++ )
     {    int v = BinPosition( verts, to_left[ dlines[i].front( )[0][0] ] );
          int w = BinPosition( verts, to_right[ dlines[i].back( )[0][0] ] );
          from[v].push_back(w), to[w].push_back(v);
          from_edge_obj[v].push_back(i), to_edge_obj[w].push_back(i);    }
     #pragma omp parallel for
     for ( int v = 0; v < N; v++ )
     {    SortSync( from[v], from_edge_obj[v] );
          SortSync( to[v], to_edge_obj[v] );    }
     digraphE<int> D2( from, to, llens, to_edge_obj, from_edge_obj );
     vec<int> to_left2, to_right2;
     D2.ToLeft(to_left2), D2.ToRight(to_right2);

     // Go through the lines.

     cout << Date( ) << ": traversing lines, mem = " << MemUsageGBString( ) << endl;
     int naccepted = 0;
     vec<int> linksto( dlines.size( ), -1 );
     vec<int> linksto_n( dlines.size( ), 0 );
     vec<vec<int>> between( dlines.size( ) );
     for ( int li = 0; li < dlines.isize( ); li++ )
     {    
          if ( llens[li] < MIN_LINE ) continue;
          int d = dlines[li].back( )[0][0];
          if ( D.O(d)[0] < 0 ) continue;

          // Start output.

          ostringstream out;

          // Find out where the line reaches forward to, using the graph.

          const int MAX_REACH = 1000;
          int v = to_right[d];
          vec<pair<int,int>> reach;
          for ( int j = 0; j < D.From(v).isize( ); j++ )
          {    int d = D.IFrom(v)[j];
               if ( D.O(d)[0] < 0 ) continue;
               int l = tol[d];
               reach.push( 0, l );    }
          vec<Bool> to_delete( reach.size( ), False );
          for ( int i = 0; i < reach.isize( ); i++ )
          {    if ( STRONG_VALIDATE && reach.isize( ) > MAX_REACH ) break;
               int dist = reach[i].first, l = reach[i].second;
               int w = to_right[ dlines[l].back( )[0][0] ];
               for ( int k = 0; k < D.From(w).isize( ); k++ )
               {    int d = D.IFrom( w, k );
                    if ( D.O(d)[0] < 0 ) continue;
                    int l2 = tol[d], dist2 = dist + llens[l];
                    if ( dist2 > MAX_DIST ) continue;
                    Bool found = False;
                    for ( int j = 0; j < reach.isize( ); j++ )
                    {    if ( reach[j].second == l2 )
                         {    if ( reach[j].first > dist2 ) to_delete[j] = True;
                              else found = True;
                              break;    }    }
                    if (found) continue;
                    reach.push( dist2, l2 );
                    if ( STRONG_VALIDATE && reach.isize( ) > MAX_REACH ) break;
                    to_delete.push_back(False);    }    }
          if ( STRONG_VALIDATE && reach.isize( ) > MAX_REACH ) continue;
          EraseIf( reach, to_delete );
          Sort(reach);
          int maxlen = 0;
          for ( auto& x : reach ) maxlen = Max( maxlen, llens[x.second] );
          if ( maxlen < MIN_LINE ) continue;
          out << "\n======================================================="
               << "=============================\n";
          out << "\nREACHING FORWARD FROM LINE " << li << endl;
          out << "using graph:\n";
          for ( int i = 0; i < reach.isize( ); i++ )
          {    int dist = reach[i].first, l = reach[i].second;
               out << "at dist " << dist << ", line " << l
                    << ", length = " << llens[l] << endl;    }
          vec<int> range;
          for ( auto& x : reach ) range.push_back( x.second );
          Sort(range);
          
          // Figure out where line reaches forward to using reads.

          out << "using reads:\n";
          vec< pair<int,int64_t> > ls;
          int rd = dinv[d];
          for ( int i = 0; i < (int) dpaths_index[d].size( ); i++ )
          {    int64_t id1 = dpaths_index[d][i];
               int64_t id2 = ( id1 % 2 == 0 ? id1+1 : id1-1 );
               const ReadPath &p1 = dpaths[id1], &p2 = dpaths[id2];
               for ( int j = 0; j < (int) p1.size( ); j++ )
               {    if ( p1[j] == d )
                    {    for ( int k = j + 1; k < (int) p1.size( ); k++ )
                         {    if ( D.O( p1[k] )[0] < 0 ) continue;
                              int l = tol[ p1[k] ];
                              ls.push( l, id1/2 );    }
                         break;    }    }
               for ( int j = 0; j < (int) p2.size( ); j++ )
               {    if ( D.O( p2[j] )[0] < 0 ) continue;
                    int l = tol[ p2[j] ];
                    ls.push( linv[l], id1/2 );    }    }
          for ( int i = 0; i < (int) dpaths_index[rd].size( ); i++ )
          {    int64_t id = dpaths_index[rd][i];
               const ReadPath& p = dpaths[id];
               for ( int j = 0; j < (int) p.size( ); j++ )
               {    if ( p[j] == rd )
                    {    for ( int k = 0; k < j; k++ )
                         {    if ( D.O( dinv[ p[k] ] )[0] < 0 ) continue;
                              int l = tol[ dinv[ p[k] ] ];
                              if ( STRONG_VALIDATE 
                                   && D.O( dlines[l].front( )[0][0] )[0] < 0 )
                              {    continue;    }
                              ls.push( l, id/2 );    }
                         break;    }    }    }
          UniqueSort(ls);
          vec<pair<int,int>> lx;
          for ( int i = 0; i < ls.isize( ); i++ )
          {    int j;
               for ( j = i + 1; j < ls.isize( ); j++ )
                    if ( ls[j].first != ls[i].first ) break;
               int l = ls[i].first;
               lx.push( j - i, l );
               i = j - 1;    }
          ReverseSort(lx);

          // Excluded certain links.

          to_delete.resize_and_set( lx.size( ), True );
          for ( int i = 0; i < lx.isize( ); i++ )
          {    int n = lx[i].first, l = lx[i].second;
               if ( n == 1 ) continue;
               Bool out_of_range = !BinMember( range, l );
               if ( n == 2 && out_of_range ) continue;
               if ( ( l == li || l == linv[li] ) && out_of_range ) continue;
               to_delete[i] = False;    }
          EraseIf( lx, to_delete );

          // Print evidence.

          for ( int i = 0; i < lx.isize( ); i++ )
          {    int n = lx[i].first, l = lx[i].second;
               out << "[" << n << "] " << l
                    << " (len=" << llens[l] << ",inv=" << linv[l] << ")";
               if ( !BinMember( range, l ) ) out << " [out of range]";
               if ( l == li ) out << " [this]";
               if ( l == linv[li] ) out << " [invthis]";
               out << endl;    }

          // Test for simple win.  Suppose that the longest edge that is linked to
          // has at least MIN_WIN links, and that it has length at least MIN_LINE,
          // and that all other links are in range and between.  Then we accept 
          // the link.

          int best = 0, best_id = -1;
          Bool out_of_range = False;
          for ( int i = 0; i < lx.isize( ); i++ )
          {    int l = lx[i].second;
               if ( !BinMember( range, l ) ) out_of_range = True;
               if ( llens[l] > best )
               {    best = llens[l], best_id = i;    }    }
          Bool accepted1 = False;
          if ( !out_of_range && best >= MIN_LINE && lx[best_id].first >= MIN_WIN )
          {    int l = lx[best_id].second;

               // Find all the lines between li and l.

               const int MAX_PATHS = 100;
               const int MAX_LOOPS = 1000;
               const int MAX_COPIES = 2;
               vec<vec<int>> paths;
               Bool aol = 
                    D2.AllPathsLengthRange( to_right2[li], to_left2[l], 0, MAX_DIST,
                    to_left2, to_right2, paths, MAX_PATHS, MAX_LOOPS, False,
                    MAX_COPIES );
               if ( !aol || paths.empty( ) ) continue;
               vec<int> bet = Contents(paths);

               // Check that all links are between.

               Bool OK = True;
               for ( int i = 0; i < lx.isize( ); i++ )
               {    int l2 = lx[i].second;
                    if ( l2 == l || l2 == li ) continue;
                    if ( !BinMember( bet, l2 ) ) OK = False;    }
               if ( l == li || l == linv[li] ) OK = False;

               // OK now maybe accept.

               if (OK) 
               {    out << "\nACCEPT LINK TO " << l << endl;
                    out << "total lines between = " << between.size( ) << endl;
                    rout << out.str( );
                    accepted1 = True;    
                    linksto[li] = l;
                    #pragma omp critical
                    {    naccepted++;    
                         linksto_n[l]++;    }
                    between[li] = bet;    }    }
          if ( !accepted1 ) out << out.str( );    }

     // Remove asymmetric links.

     cout << Date( ) << ": reviewing links" << endl;
     int pass = 0;
     for ( int l1 = 0; l1 < dlines.isize( ); l1++ )
     {    int l2 = linksto[l1];
          Bool OK = True;
          if ( l2 < 0 || linksto_n[l2] > 1 ) OK = False;
          if (OK)
          {    int rl1 = linv[l1], rl2 = linv[l2];
               if ( linksto[rl2] != rl1 ) OK = False;
               if (OK)
               {    vec<int> x = between[rl2], y;
                    for ( auto l : x ) y.push_back( linv[l] );
                    UniqueSort(y);
                    if ( y != between[l1] ) OK = False;
                    if (OK) pass++;    }    }
          if ( !OK ) 
          {    linksto[l1] = -1, between[l1].clear( );    }    }

     // Gather data for edits.

     cout << Date( ) << ": gathering data for edits" << endl;
     vec<vec<int>> ds( dlines.size( ) );
     for ( int l1 = 0; l1 < dlines.isize( ); l1++ )
     {    int l2 = linksto[l1];
          if ( l2 < 0 ) continue;
          int rl1 = linv[l1], rl2 = linv[l2];
          if ( make_pair( rl2, rl1 ) < make_pair( l1, l2 ) ) continue;
          int N = D.N( );
          int d1 = dlines[l1].back( )[0][0], d2 = dlines[l2].front( )[0][0];
          int rd1 = dinv[d1], rd2 = dinv[d2];
          for ( auto l : between[l1] )
          {    const vec<vec<vec<int>>>& L = dlines[l];
               for ( int j = 0; j < L.isize( ); j++ )
               {    const vec<vec<int>>& M = L[j];
                    ds[l1].append( Contents(M) );
                    if ( M.solo( ) && M[0].empty( ) )
                    {    int g = -1;
                         if ( j > 0 && L[j-1].nonempty( )
                              && L[j-1][0].nonempty( ) )
                         {    int d = L[j-1][0][0];
                              int v = to_right[d];
                              if ( v >= 0 ) g = D.IFrom(v,0);    }
                         if ( g >= 0 ) ds[l1].push_back(g);    }    }    }
          Sort( ds[l1] );    }

     // Introduce edits.

     cout << Date( ) << ": editing graph" << endl;
     for ( int l1 = 0; l1 < dlines.isize( ); l1++ )
     {    int l2 = linksto[l1];
          if ( l2 < 0 ) continue;
          int rl1 = linv[l1], rl2 = linv[l2];
          if ( make_pair( rl2, rl1 ) < make_pair( l1, l2 ) ) continue;
          int d1 = dlines[l1].back( )[0][0], d2 = dlines[l2].front( )[0][0];
          int rd1 = dinv[d1], rd2 = dinv[d2];
          if ( !IsUnique( vec<int>{ d1, d2, rd1, rd2 } ) ) continue;
          int E = D.E( ), n = ds[l1].size( );
          vec<int> rds;
          for ( auto d : ds[l1] ) rds.push_back( dinv[d] );
          digraphE<vec<int>> M( digraphE<vec<int>>::COMPLETE_SUBGRAPH_EDGES,
               D, ds[l1], to_left, to_right );

          // This is stupid, should just merge adjacent cells:

          /*
          if (STRONG_VALIDATE)
          {    Bool bad = False;
               for ( int v = 0; v < M.N( ); v++ )
               {    if ( M.To(v).solo( ) && M.From(v).solo( ) )
                    {    int x = M.To(v)[0], y = M.From(v)[0];
                         if ( !IsUnique( vec<int>{ x, v, y } ) ) continue;
                         int d1 = M.ITo(v,0), d2 = M.IFrom(v,0);
                         if ( !IsCell( M.O(d1) ) || !IsCell( M.O(d2) ) ) continue;
                         bad = True;    }    }
               if (bad) continue;    }
          */

          D.AppendWithUpdate( M, to_left, to_right );
          digraphE<vec<int>> RM( digraphE<vec<int>>::COMPLETE_SUBGRAPH_EDGES,
               D, rds, to_left, to_right );
          D.AppendWithUpdate( RM, to_left, to_right );
          dinv.resize( E + 2*n );
          for ( int i = 0; i < n; i++ )
          {    dinv[ E + i ] = E + n + i, dinv[ E + n + i ] = E + i;    }
          if ( ds[l1].empty( ) )
          {    int N = D.N( );
               D.AddVertices(2);
               D.GiveEdgeNewToVxWithUpdate( d1, to_right[d1], N, to_right );
               D.GiveEdgeNewFromVxWithUpdate( d2, to_left[d2], N, to_left );
               D.GiveEdgeNewToVxWithUpdate( rd2, to_right[rd2], N+1, to_right );
               D.GiveEdgeNewFromVxWithUpdate( rd1, to_left[rd1], N+1, to_left );    }
          else
          {    int v1 = -1, v2 = -1, rv1 = -1, rv2 = -1;
               for ( int i = 0; i < ds[l1].isize( ); i++ )
               {    int d = ds[l1][i];
                    if ( to_left[d] == to_right[d1] )
                    {    v1 = to_left[E+i];
                         break;    }    }
               for ( int i = 0; i < ds[l1].isize( ); i++ )
               {    int d = ds[l1][i];
                    if ( to_right[d] == to_left[d2] )
                    {    v2 = to_right[E+i];
                         break;    }    }
               for ( int i = 0; i < rds.isize( ); i++ )
               {    int d = rds[i];
                    if ( to_left[d] == to_right[rd2] )
                    {    rv2 = to_left[E+n+i];
                         break;    }    }
               for ( int i = 0; i < rds.isize( ); i++ )
               {    int d = rds[i];
                    if ( to_right[d] == to_left[rd1] )
                    {    rv1 = to_right[E+n+i];
                         break;    }    }
               D.GiveEdgeNewToVxWithUpdate( d1, to_right[d1], v1, to_right );
               D.GiveEdgeNewFromVxWithUpdate( d2, to_left[d2], v2, to_left );
               D.GiveEdgeNewToVxWithUpdate( rd2, to_right[rd2], rv2, to_right );
               D.GiveEdgeNewFromVxWithUpdate( rd1, to_left[rd1], rv1, to_left );    
                    }    }

     // Clean up.

     cout << Date( ) << ": cleaning up" << endl;
     RemoveUnneededVertices( D, dinv );
     CleanupCore( D, dinv );    }
