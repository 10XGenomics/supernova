// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>

#include "CoreTools.h"
#include "graph/DigraphTemplate.h"
#include "paths/HyperBasevector.h"
#include "paths/long/large/Lines.h"
#include "10X/Gap.h"
#include "10X/Heuristics.h"
#include "10X/Super.h"

// Capture multiple loops at a vertex.

void CaptureMultiLoops( digraphE<vec<int>>& D, vec<int>& dinv, vec<int>& dels,
     const Bool verbose, const Bool single )
{    vec<int> to_left;
     D.ToRight(to_left);
     vec< triple<int,int,int> > loops;
     int nthreads = ( single ? 1 : omp_get_max_threads( ) );
     vec< triple<int,int,int> > loopers;
     #pragma omp parallel for num_threads(nthreads)
     for ( int v = 0; v < D.N( ); v++ )
     {    int nloops = 0;
          for ( int i = 0; i < D.From(v).isize( ); i++ )
               if ( D.From(v)[i] == v ) nloops++;
          if ( nloops < 2 ) continue;
          if ( D.From(v).isize( ) != nloops + 1 ) continue;
          if ( D.To(v).isize( ) != nloops + 1 ) continue;
          vec<int> loops;
          int x = -1, y = -1;
          for ( int i = 0; i < D.To(v).isize( ); i++ )
               if ( D.To(v)[i] != v ) x = D.ITo( v, i );
          for ( int i = 0; i < D.From(v).isize( ); i++ )
          {    if ( D.From(v)[i] == v ) loops.push_back( D.IFrom( v, i ) );
               else y = D.IFrom( v, i );    }
          if ( make_pair( dinv[y], dinv[x] ) < make_pair( x, y ) ) continue;
          if ( !IsUnique( vec<int>{ x, y, dinv[x], dinv[y] } ) ) continue;
          Bool gap = False;
          for ( int i = 0; i < loops.isize( ); i++ )
               if ( D.O( loops[i] )[0] < 0 ) gap = True;
          if (gap) continue;
          #pragma omp critical
          {    loopers.push( v, x, y );    }    }
     Sort(loopers);
     for ( int i = 0; i < loopers.isize( ); i++ )
     {    int v = loopers[i].first, x = loopers[i].second, y = loopers[i].third;
          int rx = dinv[x], ry = dinv[y];
          int rv = to_left[ry];
          vec<int> loops;
          for ( int i = 0; i < D.From(v).isize( ); i++ )
               if ( D.From(v)[i] == v ) loops.push_back( D.IFrom( v, i ) );

          // Encode as cells.

          vec<int> zero( loops.size( ), 0 );
          vec<int> z( loops.size( ), vec<int>::IDENTITY );
          vec<vec<int>> from({zero}), to({zero}); 
          vec<vec<int>> from_edge_obj({z}), to_edge_obj({z});
          vec<vec<int>> edges;
          for ( int j = 0; j < loops.isize( ); j++ )
               edges.push_back( D.O( loops[j] ) );
          digraphE<vec<int>> G1( from, to, edges, to_edge_obj, from_edge_obj );
          for ( auto& d : loops ) d = dinv[d];
          edges.clear( );
          for ( int j = 0; j < loops.isize( ); j++ )
               edges.push_back( D.O( loops[j] ) );
          digraphE<vec<int>> G2( from, to, edges, to_edge_obj, from_edge_obj );
          cell c1( G1, 0, 0 ), c2( G2, 0, 0 );
          vec<int> x1, x2;
          c1.CellEncode(x1), c2.CellEncode(x2);

          // Edit to:
          //
          // --------x---------> v -------E------> N -------y------>
          //
          // <-------rx-------- rv <-----E+1------ N+1 <----ry------

          for ( int i = 0; i < D.To(v).isize( ); i++ )
               if ( D.To(v)[i] == v ) dels.push_back( D.ITo(v,i) );
          for ( int i = 0; i < D.To(rv).isize( ); i++ )
               if ( D.To(rv)[i] == rv ) dels.push_back( D.ITo(rv,i) );
          int N = D.N( ), E = D.E( );
          D.AddVertices(2);
          D.AddEdge( v, N, x1 );
          D.AddEdge( N+1, rv, x2 );
          dinv.push_back( E+1, E );
          D.GiveEdgeNewFromVx( y, v, N );
          D.GiveEdgeNewToVx( ry, rv, N+1 );    }
     if (verbose)
     {    cout << Date( ) << ": captured " << loopers.size( ) 
               << " multiloops" << endl;    }    }

void CaptureMessyLoops( const HyperBasevectorX& hb, const vec<int>& inv,
     digraphE<vec<int>>& D, vec<int>& dinv, vec<int>& dels, 
     const Bool allow_point, const int LONG_LINE )
{
     cout << Date( ) << ": capturing messy loops" << endl;
     int mess = 0;
     vec<vec<vec<vec<int>>>> dlines;
     FindLines( D, dinv, dlines, MAX_CELL_PATHS, MAX_CELL_DEPTH );
     vec<int> to_left, to_right;
     D.ToLeft(to_left), D.ToRight(to_right);
     vec<Bool> used;
     D.Used(used);
     const int END_SEARCH = 10;
     const int MAX_MESS = 20;
     const int MAX_EDGE_IN_LOOP = 2000;
     vec<int> lens( D.E( ), 0 );
     // vec<bool> visited( D.E( ), false);
     // compute all the super edge lengths
     #pragma omp parallel for
     for ( int d = 0; d < D.E( ); d++ )
     {    if ( D.O(d)[0] < 0 ) continue;
          for ( int j = 0; j < D.O(d).isize( ); j++ )
               lens[d] += hb.Kmers( D.O(d)[j] );    }
     vec<int> llens;
     GetLineLengths( hb, D, dlines, llens );
     vec< pair<int,int> > long_left, long_right;
     // long_left is a list of (start vertex of line i, i)
     // for lines that are at least LONG_LINE in length.
     // long_right (end vertex of line i, i) and length condition.
     for ( int i = 0; i < dlines.isize( ); i++ )
     {    if ( llens[i] < LONG_LINE ) continue;
          const vec<vec<vec<int>>>& L = dlines[i];
          int d1 = L.front( )[0][0], d2 = L.back( )[0][0];
          if ( used[d1] ) long_left.push( to_left[d1], i );
          if ( used[d2] ) long_right.push( to_right[d2], i );    }
     Sort(long_left), Sort(long_right);
     cout << Date( ) << ": main mess" << endl;
     for ( int i = 0; i < long_right.isize( ); i++ )
     {    int L1 = long_right[i].second, L2 = -1;
          int v = long_right[i].first, w = -1;
          // vs is the list of vertices that can be reached
          // FROM v, by exploring up to depth END_SEARCH.
          // vs lies inside the MESSY loop. And OUTSIDE L1.
          vec<int> vs = {v};
          for ( int d = 0; d < END_SEARCH; d++ )
          {    int n = vs.size( );
               for ( int j = 0; j < n; j++ )
                    vs.append( D.From( vs[j] ) );
               UniqueSort(vs);    }
          // In our exploration, did we reach the start of another line?
          // If not, then move on. Otherwise, this defines L2 and w.
          for ( int j = 0; j < vs.isize( ); j++ )
          {    int p = BinPosition1( long_left, vs[j] );
               if ( p >= 0 )
               {    w = vs[j];
                    L2 = long_left[p].second;
                    break;    }    }
          if ( w < 0 ) continue;
          int d1 = dlines[L1].back( )[0][0], d2 = dlines[L2].front( )[0][0];
          int rv = to_right[ dinv[d2] ], rw = to_left[ dinv[d1] ];
          int rd1 = dinv[d1], rd2 = dinv[d2];
          // Let's do this only once and not repeat for rev comp
          if ( make_pair( rd2, rd1 ) <= make_pair( d1, d2 ) ) continue;
          
          // make sure that we are not colliding with the rev comp
          vs = {v,w};
          if ( allow_point && v == w )
          {    if ( v == rv ) continue;
               vs = {v};    }
          else if ( !IsUnique( v, w, rv, rw ) ) continue;

          // Find the vertices in the loop thing.
          // L1--edge d1-->v:(MESS):w--edge d2-->L2
          // Collect up to MAX_MESS+2 vertices into vs that are all in the mess.

          while( vs.isize( ) <= MAX_MESS + 2 )
          {    int n = vs.size( );
               for ( int j = 0; j < n; j++ )
               {    int t = vs[j];
                    for ( int l = 0; l < D.From(t).isize( ); l++ )
                         if ( D.IFrom(t,l) != d2 ) vs.push_back( D.From(t)[l] );
                    for ( int l = 0; l < D.To(t).isize( ); l++ )
                         if ( D.ITo(t,l) != d1 ) vs.push_back( D.To(t)[l] );    }
               UniqueSort(vs);
               if ( vs.isize( ) == n ) break;    }
          // this will help us break out of the loop if we need to
          Bool bad = False;
          // skip if any edge in vs is a source or sink.
          for ( int j = 0; j < vs.isize( ); j++ )
               if ( D.Source( vs[j] ) || D.Sink( vs[j] ) ) bad = True;
          if ( vs.isize( ) > MAX_MESS + 2 || bad ) continue;

          // Let ds be the set of all edges that are part of vs, except d1 and d2.
          // All these are in the mess.

          vec<int> ds;
          for ( int j = 0; j < vs.isize( ); j++ )
          {    int t = vs[j];
               for ( int l = 0; l < D.From(t).isize( ); l++ )
               {    int d = D.IFrom(t,l);
                    if ( d != d2 ) ds.push_back(d);    }
               for ( int l = 0; l < D.To(t).isize( ); l++ )
               {    int d = D.ITo(t,l);
                    if ( d != d1 ) ds.push_back(d);    }    }
          UniqueSort(ds);
          if ( ds.empty( ) ) continue;

          // If we have already dealt with these edges before, continue.
          // since edges are sorted, we won't hit -4 edge first.

          /*
          int nds = ds.isize();
          for ( int j = 0; j < nds; j++) 
          {    if ( visited[ ds[j] ] || visited[ dinv[ ds[j] ] ] ) 
               {    bad=true;
                    break;    }    }
          if (bad) continue;
          */

          // We have captured a subset vs of vertices involved in the mess,
          // and edges ds associated with these. 
          // We want to now make sure that ALL non-d1 edges going INTO v, and 
          // ALL non-d2 edges emanating from w are marked for deletion.
          // That way when we delete these edges, the mess will become a 
          // disconnected component.
          // if not, we cannot delete, and replace with a -4 gap edge.

          for ( int j = 0; j < D.To(v).isize( ); j++ )
          {    int d = D.ITo(v,j);
               if ( d != d1 && !BinMember( ds, d ) ) bad = True;    }
          for ( int j = 0; j < D.From(w).isize( ); j++ )
          {    int d = D.IFrom(w,j);
               if ( d != d2 && !BinMember( ds, d ) ) bad = True;    }
          if (bad) continue;

          // Don't do anything is the mess contains an edge longer than 
          // MAX_EDGE_IN_LOOP.  Ditto if we have a BarcodeOnlyGap.

          for ( int j = 0; j < ds.isize( ); j++ )
          {    if ( ds[j] >= lens.isize( ) )
               {    bad = True;
                    break;    }
               if ( lens[ ds[j] ] > MAX_EDGE_IN_LOOP ) bad = True;
               if ( IsBarcodeOnlyGap( D.O( ds[j] ) ) ) bad = True;    }
          if (bad) continue;

          // OK good to go.  Expand any cell edges in ds.

          int nds = ds.size( );
          vec<Bool> to_delete( ds.size( ), False );
          for ( int j = 0; j < nds; j++ )
          {    int d = ds[j];
               if ( IsCell( D.O(d) ) )
               {    to_delete[j] = True;
                    int E = D.E( );
                    ReinsertLoop( d, hb, inv, D, dinv, to_left, to_right );
                    if ( dinv[d] != d )
                    {    for ( int f = E; f < E + (D.E( ) - E )/2; f++ )
                         {    ds.push_back(f);
                              to_delete.push_back(False);    }    }
                    else
                    {    for ( int f = E; f < D.E( ); f++ )
                         {    ds.push_back(f);    
                              to_delete.push_back(False);    }    }    }    }
          EraseIf( ds, to_delete );

          // Update vertices, as expansion may have invalidated them.

          v = to_right[d1], w = to_left[d2];
          rv = to_right[ dinv[d2] ], rw = to_left[ dinv[d1] ];

          // Form the ds edges into a cell.  Ditto for rc.

          digraphE<vec<int>> G( digraphE<vec<int>>::COMPLETE_SUBGRAPH_EDGES,
               D, ds, to_left, to_right );
          vec<int> gto_left, gto_right;
          G.ToLeft(gto_left), G.ToRight(gto_right);
          int cv = -1, cw = -1;
          for ( int pass = 1; pass <= 2; pass++ )
          for ( int j = 0; j < ds.isize( ); j++ )
          {    int d = ds[j];
               if ( to_left[d] == ( pass == 1 ? v : w ) )
               {    ( pass == 1 ? cv : cw ) = gto_left[j];
                    break;    }
               if ( to_right[d] == ( pass == 1 ? v : w ) )
               {    ( pass == 1 ? cv : cw ) = gto_right[j];
                    break;    }    }
          ForceAssertGe( cv, 0 ), ForceAssertGe( cw, 0 );
          ForceAssert( ds.nonempty( ) );
          vec<int> rds;
          for ( auto d : ds ) rds.push_back( dinv[d] );
          Sort(rds);
          cell c1( G, cv, cw );
          G.Reverse( );
          for ( int g = 0; g < G.E( ); g++ )
               G.OMutable(g) = D.O( dinv[ ds[g] ] );
          cell c2( G, cw, cv );
          vec<int> x1, x2;
          c1.CellEncode(x1), c2.CellEncode(x2);
          
          // Add the rev comp edges to ds.

          nds = ds.isize( );
          for ( int j = 0; j < nds; j++ ) ds.push_back( dinv[ ds[j] ] );

          // Make the edit.

          // for ( auto d : ds ) visited[d] = True;
          mess++;
          dels.append(ds);
          int N = D.N( ), E = D.E( );
          if ( v != w )
          {    D.AddEdgeWithUpdate( v, w, x1, to_left, to_right );
               D.AddEdgeWithUpdate( rv, rw, x2, to_left, to_right );    }
          else
          {    D.AddVertices(2);
               D.GiveEdgeNewFromVx( d2, v, N );
               to_left[d2] = N;
               D.AddEdgeWithUpdate( v, N, x1, to_left, to_right );
               D.GiveEdgeNewToVx( rd2, rv, N+1 );
               to_right[rd2] = N+1;
               D.AddEdgeWithUpdate( N+1, rv, x2, to_left, to_right );    }
          dinv.push_back( E + 1, E );    }
     cout << Date( ) << ": captured " << 2*mess << " messy loops, comprising " 
          << dels.size( ) << " edges" << endl;
     Validate( hb, inv, D, dinv );    }

// Simple loops have the following topology:
// Vertices u,v,w. Edge d from u to v. Edge f from v to w. Edge e from v to v.
// This function removes the loop e. A new edge is added with a gap label.
// f is transferred to the new vertex. 
// Symmetrically implemented on reverse complement.

void CaptureSimpleLoops( digraphE<vec<int>>& D, vec<int>& dinv, vec<int>& dels,
     const Bool verbose, const Bool single )
{    
    vec< triple<int,int,int> > loops;
    int nthreads = ( single ? 1 : omp_get_max_threads( ) );
    #pragma omp parallel for num_threads(nthreads)
    for ( int v = 0; v < D.N( ); v++ )
    {   
        // Canonically define vertices u,v,w and edges d,e,f
        // according to definition above.
        // skip if the vertex v does not have this topology.

        if ( D.From(v).size( ) != 2 || D.To(v).size( ) != 2 ) continue;
        int u = -1, w = -1, d, e, f;
        if ( D.To(v)[0] == v ) {
            u = D.To(v)[1];
            d = D.ITo(v,1);
            e = D.ITo(v,0);    }
        else if ( D.To(v)[1] == v ) {
            u = D.To(v)[0];
            d = D.ITo(v,0);
            e = D.ITo(v,1);    }
        else continue;
        if ( D.From(v)[0] == v )  {
            w = D.From(v)[1];
            f = D.IFrom(v,1);    }
        else if ( D.From(v)[1] == v ) {
            w = D.From(v)[0];
            f = D.IFrom(v,0);    }
        else continue;
        if ( !IsUnique( u, v, w ) ) continue;
        if ( D.O(e)[0] < 0 ) continue;
        int rd = dinv[d], re = dinv[e], rf = dinv[f];
        // let's not get into weird situations with reverse complements.
        if ( !IsUnique( vec<int>{ d, e, f, rd, re, rf } ) ) continue;
        #pragma omp critical
        {    loops.push( e, f, v );    }    }
     Sort(loops);

     // now we have a sorted list of loops to deal with.
     // add in the new vertices and edges, and
     // add in the edges to delete in del.

     int N = D.N( ), E = D.E( );
     D.AddVertices( loops.size( ) );
     for ( int i = 0; i < loops.isize( ); i++ ) {
        int e = loops[i].first, f = loops[i].second, v = loops[i].third;
        int re = dinv[e];
        int ri = BinPosition1( loops, re );
        ForceAssertGe( ri, 0 );

        // Encode as cell.

        vec<vec<int>> from = {{0}}, to = {{0}};
        vec<vec<int>> from_edge_obj = {{0}}, to_edge_obj = {{0}};
        vec<vec<int>> edges = { D.O(e) };
        digraphE<vec<int>> G( from, to, edges, to_edge_obj, from_edge_obj );
        cell c( G, 0, 0 );
        vec<int> x;
        c.CellEncode(x);

        // Now edit the graph.

        D.AddEdge( v, N + i, x );
        D.GiveEdgeNewFromVx( f, v, N + i );
        // modify involution. note that e != re by design.
        dinv.push_back( E + ri );
        dels.push_back(e);    }
     if (verbose) {
        cout << Date( ) << ": captured " << loops.size( ) 
             << " simple loops" << endl;    }    }

// Canonical loops.  Replace:
//
//    (x) -------> (y) ---d1--> (z) -------> (w)
//                     <--d2---
//
//    by
//
//    (x) -------> (y) --gap--> (z) -------> (w).
          
void CaptureCanonicalLoops( digraphE<vec<int>>& D, vec<int>& dinv, vec<int>& dels,
     const Bool verbose, const Bool single )
{    vec<int> to_left, to_right;
     D.ToLeft(to_left), D.ToRight(to_right);
     vec< pair<int,int> > canloops;
     int nthreads = ( single ? 1 : omp_get_max_threads( ) );
     #pragma omp parallel for schedule(dynamic, 10000), num_threads(nthreads)
     for ( int v = 0; v < D.N( ); v++ )
     {    if ( D.To(v).size( ) != 2 || D.From(v).size( ) != 1 ) continue;
          int w = D.From(v)[0];
          if ( D.From(w).size( ) != 2 || D.To(w).size( ) != 1 ) continue;
          int d1 = D.ITo(w,0), d2 = -1;
          if ( D.From(w)[0] == v ) d2 = D.IFrom(w,0);
          else if ( D.From(w)[1] == v ) d2 = D.IFrom(w,1);
          else continue;
          vec<int> x = D.To(v);
          x.append( D.From(w) );
          UniqueSort(x);
          if ( x.size( ) != 4 || D.O(d1)[0] < 0 || D.O(d2)[0] < 0 ) continue;
          int rd1 = dinv[d1], rd2 = dinv[d2];
          if ( !IsUnique( vec<int>{ d1, d2, rd1, rd2 } ) ) continue;
          #pragma omp critical
          {    canloops.push( d1, d2 );    }    }
     if (verbose)
     {    cout << Date( ) << ": captured " << canloops.size( ) 
               << " canonical loops" << endl;    }
     Sort(canloops);
     int N = D.E( );
     for ( int i = 0; i < canloops.isize( ); i++ )
     {    int d1 = canloops[i].first, d2 = canloops[i].second;
          int rd1 = dinv[d1], rd2 = dinv[d2];
          int v = to_left[d1], w = to_right[d1];

          // Encode as cell.

          vec<vec<int>> from = { {1}, {0} }, to = { {1}, {0} };
          vec<vec<int>> from_edge_obj = { {0}, {1} }, to_edge_obj = { {1}, {0} };
          vec<vec<int>> edges = { D.O(d1), D.O(d2) };

          digraphE<vec<int>> G( from, to, edges, to_edge_obj, from_edge_obj );
          cell c( G, 0, 1 );
          vec<int> x;
          c.CellEncode(x);

          // Edit graph.

          D.AddEdge( v, w, x );
          dels.push_back( d1, d2 );
          int p = BinPosition( canloops, make_pair( rd1, rd2 ) );
          ForceAssertGe( p, 0 );
          dinv.push_back( N + p );    }    }
