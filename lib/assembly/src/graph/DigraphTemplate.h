///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// This file contains some template functions from Digraph.h.  They are here in
// a separate file so that these functions do not have to be inlined, thereby 
// allowing for reduction of compilation time and executable size (in principle).
//
// See Digraph.h for notes about usage of this file.
// In particular, do not include this file to resolve link errors.
// The digraph-derived template classes are explicitly instantiated in a single
// module for each template parameter.  This is typically in the .cc file
// associated with the .h file that defines the template parameter.  So, for
// example, all the explicit instantiations of methods of digraphE<KmerPath>
// are declared in KmerPath.cc.  To resolve link errors, find the right place to
// explicitly instantiate the missing method, and add it to the list of explicit
// instantiations at the end of that file.

#ifndef DIGRAPH_TEMPLATE_H
#define DIGRAPH_TEMPLATE_H

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"
#include "Equiv.h"
#include "FeudalMimic.h"
#include "Set.h"
#include "VecUtilities.h"
#include "graph/Digraph.h"
#include <cstddef>

template<class E> vec<int> digraphE<E>:: EdgesBoundedBy( const int e1, const int e2,
     const vec<int>& to_left, const vec<int>& to_right ) const
{    int v = to_right[e1], w = to_left[e2];
     vec<int> edges, verts;
     set<int> edgesx, vertsx;
     edges.push_back( e1, e2 );
     edgesx.insert(e1), edgesx.insert(e2);
     verts.push_back( v, w );
     vertsx.insert(v), vertsx.insert(w);
     for ( int i = 0; i < verts.isize( ); i++ )
     {    int x = verts[i];
          for ( int j = 0; j < To(x).isize( ); j++ )
          {    int y = To(x)[j];
               if ( Member( vertsx, y ) ) continue;
               int e = EdgeObjectIndexByIndexTo( x, j );
               if ( e == e1 || e == e2 ) continue;
               verts.push_back(y);
               vertsx.insert(y);
               if ( Member( edgesx, e ) ) continue;
               edges.push_back(e);
               edgesx.insert(e);    }
          for ( int j = 0; j < From(x).isize( ); j++ )
          {    int y = From(x)[j];
               if ( Member( vertsx, y ) ) continue;
               int e = EdgeObjectIndexByIndexFrom( x, j );
               if ( e == e1 || e == e2 ) continue;
               verts.push_back(y);
               vertsx.insert(y);
               if ( Member( edgesx, e ) ) continue;
               edges.push_back(e);
               edgesx.insert(e);    }    }
     return edges;    }

template<class E> void digraphE<E>::InitialEdges( vec<int>& v ) const
{    v.clear( );
     for ( int x = 0; x < N( ); x++ )
     {    if ( To(x).empty( ) )
          {    for ( int j = 0; j < From(x).isize( ); j++ )
                    v.push_back( EdgeObjectIndexByIndexFrom( x, j ) );    }    }    }

template<class E> void digraphE<E>::TerminalEdges( vec<int>& v ) const
{    v.clear( );
     for ( int x = 0; x < N( ); x++ )
     {    if ( From(x).empty( ) )
          {    for ( int j = 0; j < To(x).isize( ); j++ )
                    v.push_back( EdgeObjectIndexByIndexTo( x, j ) );    }    }    }

template<class E> void digraphE<E>::Initialize(
     const ConstructorType2 constructor_type, const digraphE& g, const vec<int>& ed,
     const vec<int>& to_left, const vec<int>& to_right )
{    ForceAssertEq( (int) constructor_type, (int) COMPLETE_SUBGRAPH_EDGES );
     edges_.resize( ed.size( ) );
     for ( int i = 0; i < ed.isize( ); i++ )
          edges_[i] = g.EdgeObject( ed[i] );
     vec<int> verts;
     for ( int i = 0; i < ed.isize( ); i++ )
          verts.push_back( to_left[ ed[i] ], to_right[ ed[i] ] );
     UniqueSort(verts);
     int N = verts.size( );
     from_.resize(N), to_.resize(N);
     from_edge_obj_.resize(N), to_edge_obj_.resize(N);
     for ( int i = 0; i < ed.isize( ); i++ )
     {    int e = ed[i];
          int v = to_left[e], w = to_right[e];
          int iv = BinPosition( verts, v ), iw = BinPosition( verts, w );
          from_[iv].push_back(iw);
          from_edge_obj_[iv].push_back(i);
          to_[iw].push_back(iv);
          to_edge_obj_[iw].push_back(i);    }
     for ( int v = 0; v < this->N( ); v++ )
     {    SortSync( from_[v], from_edge_obj_[v] );
          SortSync( to_[v], to_edge_obj_[v] );    }    }

template<class E> digraphE<E>::digraphE( 
     const ConstructorType2 constructor_type, const digraphE& g, const vec<int>& ed,
     const vec<int>& to_left, const vec<int>& to_right )
{    Initialize( constructor_type, g, ed, to_left, to_right );    }

template<class EE> void digraphE<EE>::Initialize(
     const ConstructorType2 constructor_type, const digraphEX<EE>& g, 
     const vec<int>& ed )
{    ForceAssertEq( (int) constructor_type, (int) COMPLETE_SUBGRAPH_EDGES );
     edges_.resize( ed.size( ) );
     for ( int i = 0; i < ed.isize( ); i++ )
          edges_[i] = g.EdgeObject( ed[i] );
     vec<int> verts;
     for ( int i = 0; i < ed.isize( ); i++ )
          verts.push_back( g.ToLeft( ed[i] ), g.ToRight( ed[i] ) );
     UniqueSort(verts);
     int N = verts.size( );
     from_.resize(N), to_.resize(N);
     from_edge_obj_.resize(N), to_edge_obj_.resize(N);
     for ( int i = 0; i < ed.isize( ); i++ )
     {    int e = ed[i];
          int v = g.ToLeft(e), w = g.ToRight(e);
          int iv = BinPosition( verts, v ), iw = BinPosition( verts, w );
          from_[iv].push_back(iw);
          from_edge_obj_[iv].push_back(i);
          to_[iw].push_back(iv);
          to_edge_obj_[iw].push_back(i);    }
     for ( int v = 0; v < this->N( ); v++ )
     {    SortSync( from_[v], from_edge_obj_[v] );
          SortSync( to_[v], to_edge_obj_[v] );    }    }

template<class EE> digraphE<EE>::digraphE( const ConstructorType2 constructor_type, 
     const digraphEX<EE>& g, const vec<int>& ed )
{    Initialize( constructor_type, g, ed );    }

template<class E> void digraphE<E>::Initialize( 
     const ConstructorType1 constructor_type, const digraphE& g, const vec<int>& v )
{    ForceAssertEq( (int) constructor_type, (int) COMPLETE_SUBGRAPH );
     from_.resize( v.size( ) ), to_.resize( v.size( ) );
     from_edge_obj_.resize( v.size( ) ), to_edge_obj_.resize( v.size( ) );
     int edgecount = 0;
     vec<int> vsorted(v), vindex( v.size( ), vec<int>::IDENTITY );
     SortSync( vsorted, vindex );
     for ( int i = 0; i < v.isize( ); i++ )
     {    int x = v[i];
          for ( int j = 0; j < g.From(x).isize( ); j++ )
          {    int y = g.From(x)[j];
               int p2 = BinPosition( vsorted, y );
               if ( p2 < 0 ) continue;
               int i2 = vindex[p2];
               from_[i].push_back(i2);    
               to_[i2].push_back(i);
               from_edge_obj_[i].push_back(edgecount);
               to_edge_obj_[i2].push_back(edgecount);
               ++edgecount;    }    }
     edges_.reserve(edgecount);
     for ( int i = 0; i < v.isize( ); i++ )
     {    int x = v[i];
          for ( int j = 0; j < g.From(x).isize( ); j++ )
          {    int y = g.From(x)[j];
               int p2 = BinPosition( vsorted, y );
               if ( p2 < 0 ) continue;
               int i2 = vindex[p2];
               edges_.push_back( g.EdgeObjectByIndexFrom( x, j ) );    }    }
     for ( int i = 0; i < v.isize( ); i++ )
     {    SortSync( from_[i], from_edge_obj_[i] );
          SortSync( to_[i], to_edge_obj_[i] );    }    }

template<class E> digraphE<E>::digraphE( 
     const ConstructorType1 constructor_type, const digraphE& g, const vec<int>& v )
{    Initialize( constructor_type, g, v );    }

template<class E> vec<int> digraphE<E>::EdgesConnectedTo( const vec<int>& v ) const
{    vec<int> G = VerticesConnectedTo(v), e;
     for ( int x = 0; x < G.isize( ); x++ )
     {    for ( int j = 0; j < From( G[x] ).isize( ); j++ )
               e.push_back( EdgeObjectIndexByIndexFrom( G[x], j ) );
          for ( int j = 0; j < To( G[x] ).isize( ); j++ )
               e.push_back( EdgeObjectIndexByIndexTo( G[x], j ) );    }
     UniqueSort(e);
     return e;    }


template<class E> digraphE<E> digraphE<E>::Subgraph( const vec<int>& v ) const
{    digraphE result;
     result.from_.resize( v.size( ) );
     result.to_.resize( v.size( ) );
     result.from_edge_obj_.resize( v.size( ) );
     result.to_edge_obj_.resize( v.size( ) );
     int edgecount = 0;
     vec<int> vsorted(v), vindex( v.size( ), vec<int>::IDENTITY );
     SortSync( vsorted, vindex );
     for ( int i = 0; i < v.isize( ); i++ )
     {    int x = v[i];
          for ( int j = 0; j < From(x).isize( ); j++ )
          {    int y = From(x)[j];
               int p2 = BinPosition( vsorted, y );
               if ( p2 < 0 ) continue;
               int i2 = vindex[p2];
               result.from_[i].push_back(i2);
               result.to_[i2].push_back(i);
               result.from_edge_obj_[i].push_back(edgecount);
               result.to_edge_obj_[i2].push_back(edgecount);
               ++edgecount;    }    }
     result.edges_.reserve(edgecount);
     for ( int i = 0; i < v.isize( ); i++ )
     {    int x = v[i];
          for ( int j = 0; j < From(x).isize( ); j++ )
          {    int y = From(x)[j];
               int p2 = BinPosition( vsorted, y );
               if ( p2 < 0 ) continue;
               int i2 = vindex[p2];
               result.edges_.push_back( EdgeObjectByIndexFrom( x, j ) );  }  }
     for ( int i = 0; i < v.isize( ); i++ )
     {    SortSync( result.from_[i], result.from_edge_obj_[i] );
          SortSync( result.to_[i], result.to_edge_obj_[i] );    }
     return result; }

template<class F> digraphE<F>::digraphE( const ConstructorName cname,
     const digraphE& g, const equiv_rel& ee )
{    ForceAssert( cname == FROM_EDGE_EQUIV );
     ForceAssertEq( ee.Size( ), g.E( ) );
     vec<int> to_left, to_right, ereps;
     g.ToLeft(to_left), g.ToRight(to_right);
     ee.OrbitRepsAlt(ereps);
     equiv_rel ev( g.N( ) );
     vec<F> edges( ereps.size( ) );
     for ( int r = 0; r < ereps.isize( ); r++ )
     {    vec<int> o;
          ee.Orbit( ereps[r], o );
          for ( int j1 = 0; j1 < o.isize( ) - 1; j1++ )
          {    int j2 = j1 + 1;
               ForceAssert( g.O( o[j1] ) == g.O( o[j2] ) );
               ev.Join( to_left[ o[j1] ], to_left[ o[j2] ] );
               ev.Join( to_right[ o[j1] ], to_right[ o[j2] ] );    }
          edges[r] = g.O( o[0] );    }
     vec<int> vreps;
     ev.OrbitRepsAlt(vreps);
     int N = vreps.size( );
     vec<vec<int>> from(N), to(N), from_edge_obj(N), to_edge_obj(N);
     for ( int v = 0; v < N; v++ )
     {    vec<int> o;
          ev.Orbit( vreps[v], o );
          vec< pair<int,int> > we;
          for ( int j = 0; j < o.isize( ); j++ )
          {    int vp = o[j];
               for ( int k = 0; k < g.From(vp).isize( ); k++ )
               {    int wp = g.From(vp)[k], ep = g.IFrom( vp, k );
                    int w = BinPosition( vreps, ev.ClassId(wp) );
                    int e = BinPosition( ereps, ee.ClassId(ep) );
                    we.push( w, e );    }    }    
          UniqueSort(we);
          for ( int j = 0; j < we.isize( ); j++ )
          {    int w = we[j].first, e = we[j].second;
               from[v].push_back(w), to[w].push_back(v);
               from_edge_obj[v].push_back(e);
               to_edge_obj[w].push_back(e);     }     }
     for ( int v = 0; v < N; v++ )
     {    SortSync( from[v], from_edge_obj[v] );
          SortSync( to[v], to_edge_obj[v] );    }
     Initialize( from, to, edges, to_edge_obj, from_edge_obj );    }

template<class E> digraphE<E>::digraphE( const ConstructorName cname,
     const digraphE& g, const vec< vec<int> >& C )
{    ForceAssert( cname == FROM_SUBS );
     int nedges = 0;
     for ( int i = 0; i < C.isize( ); i++ )
          nedges += C[i].size( );
     edges_.reserve(nedges);
     vec<int> to_left, to_right;
     g.ToLeft(to_left), g.ToRight(to_right);
     for ( int i = 0; i < C.isize( ); i++ )
     {    for ( int j = 0; j < C[i].isize( ); j++ )
               edges_.push_back( g.EdgeObject( C[i][j] ) );    }
     for ( int pass = 1; pass <= 2; pass++ )
     {    int nverts = 0, nedges = 0;
          for ( int i = 0; i < C.isize( ); i++ )
          {    vec<int> verts;
               for ( int j = 0; j < C[i].isize( ); j++ )
                    verts.push_back( to_left[ C[i][j] ], to_right[ C[i][j] ] );
               UniqueSort(verts);
               if ( pass == 2 )
               {    for ( int j = 0; j < C[i].isize( ); j++ )
                    {    int v = BinPosition( verts, to_left[ C[i][j] ] );
                         int w = BinPosition( verts, to_right[ C[i][j] ] );
                         from_[ nverts + v ].push_back( nverts + w );
                         to_[ nverts + w ].push_back( nverts + v );
                         from_edge_obj_[ nverts + v ].push_back(nedges + j);
                         to_edge_obj_[ nverts + w ].push_back(nedges + j);    }    }
               nverts += verts.size( );
               nedges += C[i].size( );    }
          if ( pass == 1 )
          {    from_.resize(nverts), to_.resize(nverts);
               from_edge_obj_.resize(nverts), to_edge_obj_.resize(nverts);    }    }
     for ( int v = 0; v < N( ); v++ )
     {    SortSync( from_[v], from_edge_obj_[v] );
          SortSync( to_[v], to_edge_obj_[v] );    }    }

template<class E> digraphE<E>::digraphE( const digraphE& g, int n )
{    equiv_rel e;
     g.ComponentRelation(e);
     vec<int> reps, o;
     e.OrbitRepsAlt(reps);
     ForceAssertLt( n, reps.isize( ) );
     e.Orbit( reps[n], o );
     *this = g.Subgraph(o); }

template<class E> void digraphE<E>::Initialize( const int n ){
  edges_.resize(n);
  from_.resize(n);
  to_.resize(n);
  from_edge_obj_.resize(n);
  to_edge_obj_.resize(n);     
}

template<class F> void digraphE<F>::EdgeEquivConstructor( 
     const vec<F>& edges, const equiv_rel& e )
{    edges_ = edges;
     int ne = edges.size( );
     vec<int> reps;
     e.OrbitReps(reps);
     int nv = 2 * reps.isize( );
     to_edge_obj_.resize(nv);
     from_edge_obj_.resize(nv);
     to_.resize(nv);
     from_.resize(nv);
     for ( int i = 0; i < reps.isize( ); i++ )
     {    vec<int> o;
          e.Orbit( reps[i], o );
          for ( int j = 0; j < o.isize( ); j++ )
          {    from_[ 2*i ].push_back( 2*i + 1 );
               from_edge_obj_[ 2*i ].push_back( o[j] );
               to_[ 2*i + 1 ].push_back( 2*i );
               to_edge_obj_[ 2*i + 1 ].push_back( o[j] );    }    }    }

template<class F> digraphE<F>::digraphE( const vec<F>& edges, const equiv_rel& e )
{    EdgeEquivConstructor( edges, e );    }

template<class F> digraphE<F>::digraphE( const vec<F>& edges, 
     const ConstructorBehavior constructor_type )
{    edges_ = edges;
     int ne = edges.size( );
     int nv = ( constructor_type == EDGES_SEPARATE ? ne * 2 : ne + 1 );
     to_edge_obj_.resize(nv);
     from_edge_obj_.resize(nv);
     to_.resize(nv);
     from_.resize(nv);
     if ( constructor_type == EDGES_SEPARATE )
     {    for ( int i = 0; i < ne; i++ )
          {    from_[ 2*i ].push_back( 2*i + 1 );
               from_edge_obj_[ 2*i ].push_back(i);
               to_[ 2*i + 1 ].push_back( 2*i );
               to_edge_obj_[ 2*i + 1 ].push_back(i);    }    }
     else if ( constructor_type == EDGES_IN_LINE )
     {    for ( int i = 0; i <= ne; i++ )
          {    if ( i < ne )
               {    from_[i].push_back(i+1);
                    from_edge_obj_[i].push_back(i);    }
               if ( i > 0 ) 
               {    to_[i].push_back(i-1);
                    to_edge_obj_[i].push_back(i-1);    }    }    }
     else ForceAssert( 0 == 1 );    }

template<class F> void digraphE<F>::Used( vec<Bool>& used ) const
{    used.resize_and_set( edges_.size( ), False );
     for ( int i = 0; i < N( ); i++ )
     {    for ( int j = 0; j < to_edge_obj_[i].isize( ); j++ )
               used[ to_edge_obj_[i][j] ] = True;    }    }

template<class F> int digraphE<F>::UsedCount( ) const
{    vec<Bool> used;
     Used(used);
     return Sum(used);    }

template<class F> void digraphEX<F>::Used( vec<Bool>& used ) const
{    used.resize_and_set( edges_.size( ), False );
     for ( int i = 0; i < N( ); i++ )
     {    for ( size_t j = 0; j < to_edge_obj_[i].size( ); j++ )
               used[ to_edge_obj_[i][j] ] = True;    }    }

template<class F> int digraphEX<F>::UsedCount( ) const
{    vec<Bool> used;
     Used(used);
     return Sum(used);    }

template<class F> void digraphE<F>::JoinEdges( int x, const F& e )
{    if ( from_[x].size( ) != 1 || to_[x].size( ) != 1 )
     {    cout << "Problem in JoinEdges.\n";
          PRINT(x);
          cout << "edges in = " << printSeq( ToEdgeObj(x) ) << endl;
          cout << "edges out = " << printSeq( FromEdgeObj(x) ) << endl;    }
     ForceAssert( from_[x].size( ) == 1 && to_[x].size( ) == 1 );
     int v = to_[x][0], w = from_[x][0];
     ForceAssert( x != v || x != w );
     from_[x].clear( ), from_edge_obj_[x].clear( );
     to_[x].clear( ), to_edge_obj_[x].clear( );
     for ( int i = 0; i < from_[v].isize( ); i++ )
     {    if ( from_[v][i] == x )
          {    from_[v].erase( from_[v].begin( ) + i );
               from_edge_obj_[v].erase( from_edge_obj_[v].begin( ) + i );
               break;    }    }
     for ( int i = 0; i < to_[w].isize( ); i++ )
     {    if ( to_[w][i] == x )
          {    to_[w].erase( to_[w].begin( ) + i );
               to_edge_obj_[w].erase( to_edge_obj_[w].begin( ) + i );
               break;    }    }
     AddEdge( v, w, e );    }

template<class F> void digraphE<F>::RemoveUnneededVertices( )
{    for ( int i = 0; i < N( ); i++ )
     {    if ( From(i).size( ) == 1 && To(i).size( ) == 1 && From(i)[0] != i )
          {    F p = EdgeObjectByIndexTo( i, 0 );
               p.append( EdgeObjectByIndexFrom( i, 0 ) );
               JoinEdges( i, p );    }    }
     RemoveEdgelessVertices( );    }

// Input is a set of vertices v.  Each v must be located at the opening of
// a bubble, with exactly two edges that lead to the same successor w:
//        _-_
//  --> v     w -->
//        -_-
template<class E> void digraphE<E>::PopBubbles( const vec<int> & bubble_vs )
{
  vec<int> bubble_edges;
  bubble_edges.reserve( bubble_vs.size() );
  
  for ( int i = 0; i < bubble_vs.isize(); i++ ) {
    int v = bubble_vs[i];
    ForceAssertEq( from_[v].size(), 2u );
    ForceAssertEq( from_[v][0], from_[v][1] );
    
    // Choose one of the edges that make up this bubble, and delete it.
    // Arbitrarily, we choose the higher-indexed path.
    bubble_edges.push_back( from_edge_obj_[v][1] );
  }
  
  DeleteEdges( bubble_edges );
  // Combine edges.  For bubbles v->w in which v had only one predecessor
  // and/or w had only one successor, this will combine the remaining edge in
  // the bubble with the edge leading to/from the bubble.
  RemoveUnneededVertices( );
  // Clear out edges that have been removed from the graph.
  RemoveDeadEdgeObjects( );
}

// Input is a set of vertices v.  Each v must be located at the opening of
// a bubble, with two or more edges that lead to the same successor w:
//         _-_
//        -   -
//  --> v ----- w -->
//        _   _
//         -_-
template<class E> void digraphE<E>::PopHyperBubbles( const vec<int> & bubble_vs )
{
  vec<int> bubble_edges;
  bubble_edges.reserve( bubble_vs.size() );
  
  for ( int i = 0; i < bubble_vs.isize(); i++ ) {
    int v = bubble_vs[i];
    ForceAssertGe( from_[v].size(), 2u );
    ForceAssertEq( Min(from_[v]), Max(from_[v]) );
    
    // Choose one of the edges that make up this bubble, and delete it.
    // Arbitrarily, we choose the higher-indexed path.
    for ( size_t ib = 1; ib < from_edge_obj_[v].size(); ib++ )
      bubble_edges.push_back( from_edge_obj_[v][ib] );
  }
  
  DeleteEdges( bubble_edges );
  // Combine edges.  For bubbles v->w in which v had only one predecessor
  // and/or w had only one successor, this will combine the remaining edge in
  // the bubble with the edge leading to/from the bubble.
  RemoveUnneededVertices( );
  // Clear out edges that have been removed from the graph.
  RemoveDeadEdgeObjects( );
}

template<class E> void digraphE<E>::RemoveEdgelessVertices( 
     const vec<int>& to_remove )
{    vec<Bool> remove( N( ), False );
     for ( int i = 0; i < to_remove.isize( ); i++ )
          remove[ to_remove[i] ] = True;
     vec<int> new_vertex_id( N( ), -1 );
     int id = 0;
     for ( int i = 0; i < N( ); i++ )
     {    if ( remove[i] )
          {    ForceAssert( from_[i].empty( ) );
               ForceAssert( to_[i].empty( ) );    }
          else
          {    new_vertex_id[i] = id;
               ++id;    }    }
     EraseIf( from_, remove ), EraseIf( from_edge_obj_, remove );
     EraseIf( to_, remove ), EraseIf( to_edge_obj_, remove );
     for ( int i = 0; i < N( ); i++ )
     {    for ( int j = 0; j < from_[i].isize( ); j++ )
               from_[i][j] = new_vertex_id[ from_[i][j] ];
          for ( int j = 0; j < to_[i].isize( ); j++ )
               to_[i][j] = new_vertex_id[ to_[i][j] ];    }    }

template<class E> void digraphE<E>::RemoveEdgelessVertices( )
{    vec<int> to_remove;
     for ( int i = 0; i < N( ); i++ )
          if ( from_[i].empty( ) && to_[i].empty( ) ) to_remove.push_back(i);
     RemoveEdgelessVertices(to_remove);    }

template<class V> void digraphV<V>::RemoveEdgelessVertices( )
{    vec<int> new_vertex_id( N( ), -1 );
     int id = 0;
     vec<Bool> remove( N( ), False );
     for ( int i = 0; i < N( ); i++ )
     {    if ( from_[i].empty( ) && to_[i].empty( ) ) remove[i] = True;
          else
          {    new_vertex_id[i] = id;
               ++id;    }    }
     EraseIf( from_, remove ), EraseIf( to_, remove ), EraseIf( verts_, remove );
     for ( int i = 0; i < N( ); i++ )
     {    for ( int j = 0; j < from_[i].isize( ); j++ )
               from_[i][j] = new_vertex_id[ from_[i][j] ];
          for ( int j = 0; j < to_[i].isize( ); j++ )
               to_[i][j] = new_vertex_id[ to_[i][j] ];    }    }

template<class E> void digraphE<E>::Reverse( )
{    for ( int i = 0; i < N( ); i++ )
     {    swap( from_[i], to_[i] );
          swap( from_edge_obj_[i], to_edge_obj_[i] );    }    }

template<class E> void digraphE<E>::ReverseComponent( int x )
{    equiv_rel e( N( ) );
     for ( int v = 0; v < N( ); v++ )
     {    for ( int i = 0; i < from_[v].isize( ); i++ )
          {    int w = from_[v][i];
               e.Join( v, w );    }    }
     vec<int> o;
     e.Orbit( x, o );
     for ( int j = 0; j < o.isize( ); j++ )
     {    int i = o[j];
          swap( from_[i], to_[i] );
          swap( from_edge_obj_[i], to_edge_obj_[i] );    }    }

template<class E> void digraphE<E>::ReorderVertices( const vec<int>& new_order )
{    ForceAssertEq( new_order.isize( ), N( ) );
     vec<int> order_new( N( ) );
     for ( int i = 0; i < N( ); i++ )
          order_new[ new_order[i] ] = i;
     PermuteVec( from_, order_new );
     PermuteVec( from_edge_obj_, order_new );
     PermuteVec( to_, order_new );
     PermuteVec( to_edge_obj_, order_new );
     for ( int v = 0; v < N( ); v++ )
     {    for ( int j = 0; j < from_[v].isize( ); j++ )
               from_[v][j] = order_new[ from_[v][j] ];
          for ( int j = 0; j < to_[v].isize( ); j++ )
               to_[v][j] = order_new[ to_[v][j] ];
          SortSync( from_[v], from_edge_obj_[v] );
          SortSync( to_[v], to_edge_obj_[v] );    }    }

template<class E> void digraphE<E>::ReorderComponents( const vec<int>& new_order )
{    equiv_rel e( N( ) );
     for ( int v = 0; v < N( ); v++ )
     {    for ( int i = 0; i < from_[v].isize( ); i++ )
          {    int w = from_[v][i];
               e.Join( v, w );    }    }
     vec<int> reps;
     for ( int v = 0; v < N( ); v++ )
          if ( e.Representative(v) ) reps.push_back(v);
     ForceAssertEq( new_order.size( ), reps.size( ) );
     vec<int> new_vertex_order;
     for ( int i = 0; i < reps.isize( ); i++ )
     {    int v = reps[ new_order[i] ];
          vec<int> o;
          e.Orbit( v, o );
          new_vertex_order.append(o);    }
     ReorderVertices(new_vertex_order);    }

template<class E> void
digraphE<E>::ComponentEdges( vec< vec<edge_t> >& edges ) const
{
  vec<vec<int> > vertices;
  Components( vertices );
  int n = vertices.isize( );
  
  edges.resize( 0 );
  edges.resize( n );
  for ( int i = 0; i < n; i++ ) {
    for ( int j = 0; j < vertices[i].isize( ); j++ )
      edges[i].append( FromEdgeObj( vertices[i][j] ) );
    UniqueSort( edges[i] );
  }
}

template<class F> void digraphE<F>::Append( const digraphE<F>& D )
{    int nedges = edges_.size( );
     edges_.append( D.edges_ );
     int nvertices = from_.size( );
     from_.append( D.from_ );
     to_.append( D.to_ );
     from_edge_obj_.append( D.from_edge_obj_ );
     to_edge_obj_.append( D.to_edge_obj_ );
     for ( int i = nvertices; i < N( ); i++ )
     {    for ( int j = 0; j < from_[i].isize( ); j++ )
          {    from_[i][j] += nvertices;
               from_edge_obj_[i][j] += nedges;    }
          for ( int j = 0; j < to_[i].isize( ); j++ )
          {    to_[i][j] += nvertices;
               to_edge_obj_[i][j] += nedges;    }    }    }

template<class F> int digraphE<F>::AddEdgeWithUpdate( 
     const int v, const int w, const F& e, vec<int>& to_left, vec<int>& to_right )
{    to_left.push_back(v), to_right.push_back(w);
     return AddEdge( v, w, e );    }

template<class F> void digraphE<F>::AppendWithUpdate( 
     const digraphE<F>& D, vec<int>& to_left, vec<int>& to_right )
{    to_left.resize( to_left.isize( ) + D.E( ), -1 );
     to_right.resize( to_right.isize( ) + D.E( ), -1 );
     for ( int v = 0; v < D.N( ); v++ )
     for ( int j = 0; j < D.From(v).isize( ); j++ )
     {    int w = D.From(v)[j], e = D.IFrom(v,j);
          to_left[ e + E( ) ] = v + N( );
          to_right[ e + E( ) ] = w + N( );    }
     Append(D);    }

template<class F> void digraphE<F>::SplitEdge( int v, int j, const F& e1, const F& e2 )
{    int n = N( );
     int ne = edges_.size( );
     edges_.push_back( e1, e2 );
     int w = from_[v][j];
     int we = from_edge_obj_[v][j];
     int i = InputFromOutputTo( v, j );
     from_[v].erase( from_[v].begin( ) + j );
     from_edge_obj_[v].erase( from_edge_obj_[v].begin( ) + j );
     to_[w].erase( to_[w].begin( ) + i );
     to_edge_obj_[w].erase( to_edge_obj_[w].begin( ) + i );
     from_[v].push_back(n), from_edge_obj_[v].push_back(ne);
     vec<int> nfrom, nto;
     vec<int> nfrom_edge_obj, nto_edge_obj;
     nfrom.push_back(w), nfrom_edge_obj.push_back(ne+1);
     nto.push_back(v), nto_edge_obj.push_back(ne);
     from_.push_back(nfrom), to_.push_back(nto);
     from_edge_obj_.push_back(nfrom_edge_obj);
     to_edge_obj_.push_back(nto_edge_obj);
     for ( int u = 0; u < to_[w].isize( ); u++ )
     {    if ( to_[w][u] == v && we == to_edge_obj_[w][u] )
          {    to_.erase( to_.begin( ) + u );
               to_edge_obj_.erase( to_edge_obj_.begin( ) + u );
               break;    }    }
     to_[w].push_back(n), to_edge_obj_[w].push_back(ne+1);    }

template<class F> void digraphE<F>::Glue( const EmbeddedSubPath<F>& a,
     const EmbeddedSubPath<F>& b, const vec<int>& EE, const vec<int>& FF, 
     const digraphE<F>& c )
{    
     // Sanity check.

     ForceAssertGe( a.NVertices( ), 2 ); ForceAssertGe( b.NVertices( ), 2 );
     ForceAssert( !HasSharedEdge(a, b) );
     ForceAssertEq( EE.isize( ), a.NVertices( ) ); 
     ForceAssertEq( FF.isize( ), b.NVertices( ) );
     vec<int> Esort = EE, Fsort = FF;
     Sort(Esort), Sort(Fsort);
     ForceAssert( Esort.UniqueOrdered( ) );
     ForceAssert( Fsort.UniqueOrdered( ) );
     ForceAssertEq( EE.front( ), 0 ); ForceAssertEq( EE.back( ), c.N( ) - 1 );
     ForceAssertEq( FF.front( ), 0 ); ForceAssertEq( FF.back( ), c.N( ) - 1 );

     // Delete edges appearing in a and b.

     for ( int i = 0; i < a.NVertices( ) - 1; i++ )
     {    int v = a.Vertex(i), w = a.Vertex(i+1);
          int e = a.EdgeObjectIndexAbs(i);
          int ef = EdgeObjectIndexToFromIndex( v, e );
          int et = InputFromOutputTo( v, ef );
          from_[v].erase( from_[v].begin( ) + ef );
          from_edge_obj_[v].erase( from_edge_obj_[v].begin( ) + ef );
          to_[w].erase( to_[w].begin( ) + et );
          to_edge_obj_[w].erase( to_edge_obj_[w].begin( ) + et );    }
     for ( int i = 0; i < b.NVertices( ) - 1; i++ )
     {    int v = b.Vertex(i), w = b.Vertex(i+1);
          int e = b.EdgeObjectIndexAbs(i);
          int ef = EdgeObjectIndexToFromIndex( v, e );
          int et = InputFromOutputTo( v, ef );
          from_[v].erase( from_[v].begin( ) + ef );
          from_edge_obj_[v].erase( from_edge_obj_[v].begin( ) + ef );
          to_[w].erase( to_[w].begin( ) + et );
          to_edge_obj_[w].erase( to_edge_obj_[w].begin( ) + et );    }

     // Attach c.

     int nvertices = N( );
     Append(c);
     for ( int i = 0; i < a.NVertices( ); i++ )
          TransferEdges( a.Vertex(i), EE[i] + nvertices );
     for ( int i = 0; i < b.NVertices( ); i++ )
          TransferEdges( b.Vertex(i), FF[i] + nvertices );    

     // If data implies that some vertices in c should be identified, do so.

     vec< vec<int> > sources( c.N( ) );
     for ( int i = 0; i < a.NVertices( ); i++ )
          sources[ EE[i] ].push_back( a.Vertex(i) );
     for ( int i = 0; i < b.NVertices( ); i++ )
          sources[ FF[i] ].push_back( b.Vertex(i) );
     for ( int i = 0; i < c.N( ); i++ )
          Sort( sources[i] );
     for ( int i1 = 0; i1 < c.N( ); i1++ )
     {    for ( int i2 = i1 + 1; i2 < c.N( ); i2++ )
          {    if ( Meet( sources[i1], sources[i2] ) )
                    TransferEdges( i1 + nvertices, i2 + nvertices );    }    }    }

template<class E> void digraphE<E>::TransferEdges( int v, int w, 
     const Bool enter_only )
{    ForceAssert( v != w );

     // Change edges v --> v to edges w --> w.

     if ( !enter_only )
     {
     vec<Bool> remove_from_v;
     remove_from_v.resize_and_set( from_[v].size( ), False );
     for ( int i = 0; i < from_[v].isize( ); i++ )
     {    if ( from_[v][i] == v )
          {    from_[w].push_back(w);
               from_edge_obj_[w].push_back( from_edge_obj_[v][i] );
               to_[w].push_back(w);
               to_edge_obj_[w].push_back( from_edge_obj_[v][i] );
               remove_from_v[i] = True;
               int j = InputFromOutputTo( v, i );
               to_[v].erase( to_[v].begin( ) + j );
               to_edge_obj_[v].erase( to_edge_obj_[v].begin( ) + j );    }    }
     EraseIf( from_[v], remove_from_v );
     EraseIf( from_edge_obj_[v], remove_from_v );
     SortSync( from_[w], from_edge_obj_[w] );
     SortSync( to_[w], to_edge_obj_[w] );
     }

     // Change edges u --> v to edges u --> w.
     
     for ( int i = 0; i < to_[v].isize( ); i++ )
     {    int u = to_[v][i];
          int j = InputToOutputFrom( v, i );
          from_[u][j] = w;
          SortSync( from_[u], from_edge_obj_[u] );    }

     // Change edges v --> x to edges w --> x.

     // if ( !enter_only )
     {    for ( int i = 0; i < from_[v].isize( ); i++ )
          {    int x = from_[v][i];
               int j = InputFromOutputTo( v, i );
               if ( !enter_only ) to_[x][j] = w;
               else to_[x][j] = v;
               SortSync( to_[x], to_edge_obj_[x] );    }    }

     // Do the rest.

     if ( !enter_only )
     {    from_[w].append( from_[v] );
          from_edge_obj_[w].append( from_edge_obj_[v] );    }
     SortSync( from_[w], from_edge_obj_[w] );
     to_[w].append( to_[v] );
     to_edge_obj_[w].append( to_edge_obj_[v] );
     SortSync( to_[w], to_edge_obj_[w] );
     to_[v].clear( ), to_edge_obj_[v].clear( );    
     if ( !enter_only ) { from_[v].clear( ), from_edge_obj_[v].clear( ); }    }

template<class E> void digraphE<E>::TransferEdgesWithUpdate( int v, int w, 
          vec<int>& to_left, vec<int>& to_right, const Bool enter_only )
{    vec<int> x = FromEdgeObj(v), y = ToEdgeObj(v);
     TransferEdges( v, w, enter_only );
     if ( !enter_only )
     {    for ( int i = 0; i < x.isize( ); i++ )
               to_left[ x[i] ] = w;    }
     for ( int i = 0; i < y.isize( ); i++ )
          to_right[ y[i] ] = w;    }

template<class E> void digraphE<E>::RemoveDuplicateEdges( )
{    for ( int v = 0; v < N( ); v++ )
     {    vec<Bool> remove;
          remove.resize_and_set( from_[v].size( ), False );
          for ( int j = 0; j < from_[v].isize( ); j++ )
          {    int k;
               for ( k = j + 1; k < from_[v].isize( ); k++ )
                    if ( from_[v][k] != from_[v][j] ) break;
               for ( int u1 = j; u1 < k; u1++ )
               {    if ( remove[u1] ) continue;
                    for ( int u2 = u1 + 1; u2 < k; u2++ )
                    {    if ( edges_[ from_edge_obj_[v][u1] ]
                                   == edges_[ from_edge_obj_[v][u2] ] )
                         {    remove[u2] = True;    }    }    }
               j = k - 1;    }
          for ( int i = 0; i < remove.isize( ); i++ )
          {    if ( remove[i] )
               {    int w = from_[v][i];
	            int j = InputFromOutputTo( v, i );
                    to_[w].erase( to_[w].begin( ) + j );
                    to_edge_obj_[w].erase( to_edge_obj_[w].begin( ) + j );    }    }
          EraseIf( from_[v], remove );
          EraseIf( from_edge_obj_[v], remove );    }    }

template<class E> void digraphE<E>::DeleteEdgesAtVertex( int v )
{    for ( int i = 0; i < from_[v].isize( ); i++ )
     {    int w = from_[v][i];
          int j = InputFromOutputTo( v, i );
          if ( v == w ) continue;
          to_[w].erase( to_[w].begin( ) + j );
          to_edge_obj_[w].erase( to_edge_obj_[w].begin( ) + j );    }
     for ( int i = 0; i < to_[v].isize( ); i++ )
     {    int w = to_[v][i];
          int j = InputToOutputFrom( v, i );
          if ( v == w ) continue;
          from_[w].erase( from_[w].begin( ) + j );
          from_edge_obj_[w].erase( from_edge_obj_[w].begin( ) + j );    }
     from_[v].clear( ), from_edge_obj_[v].clear( );
     to_[v].clear( ), to_edge_obj_[v].clear( );    }
                    
template<class E> vec<int> digraphE<E>::RemoveDeadEdgeObjects( )
{    vec<Bool> used;
     Used(used);
     int count = 0;
     vec<int> to_new_id( edges_.size( ), -1 );
     for ( int i = 0; i < edges_.isize( ); i++ )
     {    if ( used[i] )
          {    if ( count != i ) edges_[count] = edges_[i];
               to_new_id[i] = count;
               ++count;    }    }
     edges_.resize(count);
     for ( int v = 0; v < N( ); v++ )
     {    for ( int i = 0; i < from_[v].isize( ); i++ )
               from_edge_obj_[v][i] = to_new_id[ from_edge_obj_[v][i] ];
          for ( int i = 0; i < to_[v].isize( ); i++ )
               to_edge_obj_[v][i] = to_new_id[ to_edge_obj_[v][i] ];    }

     return to_new_id;
}

template<class E> Bool digraphE<E>::TestValid( const Bool exit ) const
{    if ( !digraph(*this).TestValid( ) ) return False;
     if ( from_edge_obj_.size( ) != to_edge_obj_.size( ) )
         DIGRAPH_INVALID( "sizes of from_edge_obj_ and to_edge_obj_ are different", exit );
     if ( from_.size( ) != from_edge_obj_.size( ) )
         DIGRAPH_INVALID( "sizes of from_ and from_edge_obj_ are different", exit );
     for ( int v = 0; v < N( ); v++ )
     {    if ( from_[v].size( ) != from_edge_obj_[v].size( ) )
          {    DIGRAPH_INVALID( "sizes of from_[" << v << "] and "
                    << "from_edge_obj_[" << v << "] are different", exit );    }    }
     for ( int v = 0; v < N( ); v++ )
     {    if ( to_[v].size( ) != to_edge_obj_[v].size( ) )
          {    DIGRAPH_INVALID( "sizes of to_[" << v << "] and "
                    << "to_edge_obj_[" << v << "] are different", exit );    }    }
     for ( int v = 0; v < N( ); v++ )
     {    for ( int j = 0; j < from_[v].isize( ); j++ )
          {    int w = from_[v][j];
               int ei = from_edge_obj_[v][j];
               if ( ei < 0 || ei >= EdgeObjectCount( ) )
                    DIGRAPH_INVALID( "Illegal from_edge_obj value.", exit );
               Bool found = False;
               for ( int r = 0; r < to_[w].isize( ); r++ )
                    if ( to_[w][r] == v && to_edge_obj_[w][r] == ei ) found = True;
               if ( !found )
               {    DIGRAPH_INVALID( "There is an edge from " << v << " to " << w
                         << " in from_[" << v 
                         << "], but not in to_[" << w << "].", exit );    }    }    }
     for ( int v = 0; v < N( ); v++ )
     {    for ( int j = 0; j < to_[v].isize( ); j++ )
          {    int w = to_[v][j];
               int ei = to_edge_obj_[v][j];
               if ( ei < 0 || ei >= EdgeObjectCount( ) )
                    DIGRAPH_INVALID( "Illegal to_edge_obj value.", exit );
               Bool found = False;
               for ( int r = 0; r < from_[w].isize( ); r++ )
               {    if ( from_[w][r] == v && from_edge_obj_[w][r] == ei ) 
                         found = True;    }
               if ( !found )
               {    DIGRAPH_INVALID( "There is an edge from " << v << " to " << w
                         << " in to_[" << v << "], but not in from_[" << w 
                         << "].", exit );    }    }    }
     return True;    }

template<class E> Bool digraphE<E>::TestValidParallel( ) const
{    if ( !digraph(*this).TestValidParallel( ) ) return False;
     if ( from_edge_obj_.size( ) != to_edge_obj_.size( ) )
         DIGRAPH_INVALIDP( "sizes of from_edge_obj_ and to_edge_obj_ are different"  );
     if ( from_.size( ) != from_edge_obj_.size( ) )
         DIGRAPH_INVALIDP( "sizes of from_ and from_edge_obj_ are different" );
     #pragma omp parallel for
     for ( int v = 0; v < N( ); v++ )
     {    if ( from_[v].size( ) != from_edge_obj_[v].size( ) )
          {    DIGRAPH_INVALIDP( "sizes of from_[" << v << "] and "
                    << "from_edge_obj_[" << v << "] are different" );    }    }
     #pragma omp parallel for
     for ( int v = 0; v < N( ); v++ )
     {    if ( to_[v].size( ) != to_edge_obj_[v].size( ) )
          {    DIGRAPH_INVALIDP( "sizes of to_[" << v << "] and "
                    << "to_edge_obj_[" << v << "] are different" );    }    }
     #pragma omp parallel for
     for ( int v = 0; v < N( ); v++ )
     {    for ( int j = 0; j < from_[v].isize( ); j++ )
          {    int w = from_[v][j];
               int ei = from_edge_obj_[v][j];
               if ( ei < 0 || ei >= EdgeObjectCount( ) )
                    DIGRAPH_INVALIDP( "Illegal from_edge_obj value." );
               Bool found = False;
               for ( int r = 0; r < to_[w].isize( ); r++ )
                    if ( to_[w][r] == v && to_edge_obj_[w][r] == ei ) found = True;
               if ( !found )
               {    DIGRAPH_INVALIDP( "There is an edge from " << v << " to " << w
                         << " in from_[" << v 
                         << "], but not in to_[" << w << "]." );    }    }    }
     #pragma omp parallel for
     for ( int v = 0; v < N( ); v++ )
     {    for ( int j = 0; j < to_[v].isize( ); j++ )
          {    int w = to_[v][j];
               int ei = to_edge_obj_[v][j];
               if ( ei < 0 || ei >= EdgeObjectCount( ) )
                    DIGRAPH_INVALIDP( "Illegal to_edge_obj value." );
               Bool found = False;
               for ( int r = 0; r < from_[w].isize( ); r++ )
               {    if ( from_[w][r] == v && from_edge_obj_[w][r] == ei ) 
                         found = True;    }
               if ( !found )
               {    DIGRAPH_INVALIDP( "There is an edge from " << v << " to " << w
                         << " in to_[" << v << "], but not in from_[" << w 
                         << "]." );    }    }    }
     return True;    }

template<class F>
void digraphE<F>::Initialize( const vec< vec<int> >& from, const vec< vec<int> >& to,
          const vec<F>& edges, const vec< vec<int> >& to_edge_obj,
          const vec< vec<int> >& from_edge_obj, const Bool allow_unused_edges )
{    digraph::Initialize( from, to );
     edges_ = edges; 
     to_edge_obj_ = to_edge_obj;
     from_edge_obj_ = from_edge_obj;
     int N = from.size( );
     ForceAssertEq( N, to_edge_obj.isize( ) );
     ForceAssertEq( N, from_edge_obj.isize( ) );
     vec<int> used( edges.size( ), 0 );
     for ( int i = 0; i < N; i++ ) 
     {    ForceAssertEq( to_edge_obj[i].size( ), to[i].size( ) );
          ForceAssertEq( from_edge_obj[i].size( ), from[i].size( ) );
          for ( int j = 0; j < to_edge_obj[i].isize( ); j++ )
          {    int o = to_edge_obj[i][j];
               ForceAssertGe( o, 0 );
               ForceAssertLt( o, edges.isize( ) );
               ++used[o];
               int w = i, v = to_[i][j];
               int wf = BinPosition( from[v], w );
               // The following assert won't do what we want if there are multiple
               // edges between two given vertices (in which case wf doesn't
               // make sense).
               // ForceAssertEq( o, from_edge_obj[v][wf] );    
                    }    }
     for ( int i = 0; i < used.isize( ); i++ )
     {    if ( used[i] > 1 || ( !allow_unused_edges && used[i] == 0 ) )
          {    cout << "Edge " << i << " is used " << used[i]
                    << " times, whereas it should be used exactly once.\n";    }
          if (allow_unused_edges) ForceAssertLe( used[i], 1 );
          else ForceAssertEq( used[i], 1 );    }    }

template<class F>
digraphE<F>::digraphE( const vec< vec<int> >& from, const vec< vec<int> >& to,
          const vec<F>& edges, const vec< vec<int> >& to_edge_obj,
          const vec< vec<int> >& from_edge_obj, const Bool allow_unused_edges )
      : digraph(from, to) // redundant with initialize?
{    Initialize( from, to, edges, to_edge_obj, from_edge_obj, 
          allow_unused_edges );    }

template<class V>
void digraphV<V>::Initialize( const vec< vec<int> >& from, const vec< vec<int> >& to,
          const vec<V>& verts )
{    digraph::Initialize( from, to );
     verts_ = verts;
     ForceAssertEq( N( ), verts.isize( ) );    }

template<class V>
digraphV<V>::digraphV( const vec< vec<int> >& from, const vec< vec<int> >& to,
          const vec<V>& verts ) 
     : digraph(from, to) // redundant with initialize?
{    Initialize( from, to, verts );    }

template<class V, class E>
void digraphVE<V,E>::Initialize( const vec< vec<int> >& from, 
     const vec< vec<int> >& to, const vec<V>& verts, const vec<E>& edges, 
     const vec< vec<int> >& to_edge_obj, const vec< vec<int> >& from_edge_obj )
{    digraphE<E>::Initialize( from, to, edges, to_edge_obj, from_edge_obj );
     verts_ = verts;
     ForceAssertEq( from.size( ), verts.size( ) );    }

template<class V, class E>
digraphVE<V,E>::digraphVE( const vec< vec<int> >& from, 
     const vec< vec<int> >& to, const vec<V>& verts, const vec<E>& edges, 
     const vec< vec<int> >& to_edge_obj, const vec< vec<int> >& from_edge_obj )
     : digraphE<E>( from, to, edges, to_edge_obj, from_edge_obj ) // redundant??
{    Initialize( from, to, verts, edges, to_edge_obj, from_edge_obj );    }

template<class V, class E>
digraphVE<V,E>::digraphVE( const digraphE<E>& G, const vec<V>& verts )
     : digraphE<E>(G)
{    verts_ = verts;
     ForceAssertEq( G.N( ), verts.isize( ) );    }

template<class E> Bool digraphE<E>::IsComplete( 
     const vec<int>& vertices, const vec<int>& edges ) const
{    ForceAssert( vertices.UniqueOrdered( ) );
     ForceAssert( edges.UniqueOrdered( ) );
     for ( int u = 0; u < vertices.isize( ); u++ )
     {    int v = vertices[u];
          for ( int j = 0; j < From(v).isize( ); j++ )
          {    int w = From(v)[j];
               if ( !BinMember( vertices, w ) ) return False;
               int e = EdgeObjectIndexByIndexFrom( v, j );
               if ( !BinMember( edges, e ) ) return False;    }
          for ( int j = 0; j < To(v).isize( ); j++ )
          {    int w = To(v)[j];
               if ( !BinMember( vertices, w ) ) return False;
               int e = EdgeObjectIndexByIndexTo( v, j );
               if ( !BinMember( edges, e ) ) return False;    }    }
     return True;    }

template<class E> void digraphE<E>::DualComponentRelation( 
     equiv_rel& e, const vec<Bool>& exclude ) const
{    e.Initialize( EdgeObjectCount( ) );
     for ( int v = 0; v < N( ); v++ )
     {    for ( int j1 = 0; j1 < To(v).isize( ); j1++ )
          {    int e1 = EdgeObjectIndexByIndexTo( v, j1 );
               if ( exclude.nonempty( ) && exclude[e1] ) continue;
               for ( int j2 = 0; j2 < From(v).isize( ); j2++ )
               {    int e2 = EdgeObjectIndexByIndexFrom( v, j2 );
                    if ( exclude.nonempty( ) && exclude[e2] ) continue;
                    e.Join(e1, e2);    }     }    }    }

template<class E> void digraphE<E>::Initialize( 
     const digraphE& g, const equiv_rel& e )
{    edges_ = g.edges_;
     vec<int> reps;
     e.OrbitRepsAlt(reps);
     int nreps = reps.size( );
     vec<int> to_reps( g.N( ) );
     for ( int i = 0; i < nreps; i++ )
     {    vec<int> o;
          e.Orbit( reps[i], o );
          for ( int j = 0; j < o.isize( ); j++ )
               to_reps[ o[j] ] = i;    }
     from_.resize(nreps), to_.resize(nreps);
     from_edge_obj_.resize(nreps), to_edge_obj_.resize(nreps);
     int nedges = g.EdgeObjectCount( );
     vec<int> to_left_vertex(nedges, -1), to_right_vertex(nedges, -1);
     for ( int w = 0; w < g.N( ); w++ )
     {    for ( int j = 0; j < g.To(w).isize( ); j++ )
          {    int m = g.EdgeObjectIndexByIndexTo( w, j );
               int v = g.To(w)[j];
               to_left_vertex[m] = v, to_right_vertex[m] = w;    }    }
     for ( int m = 0; m < nedges; m++ )
     {    if ( to_left_vertex[m] < 0 || to_right_vertex[m] < 0 ) continue;
          int v = to_reps[ to_left_vertex[m] ];
          int w = to_reps[ to_right_vertex[m] ];
          from_[v].push_back(w), to_[w].push_back(v);
          from_edge_obj_[v].push_back(m), to_edge_obj_[w].push_back(m);    }
     for ( int v = 0; v < N( ); v++ )
     {    SortSync( from_[v], from_edge_obj_[v] );
          SortSync( to_[v], to_edge_obj_[v] );    }    }

template<class E> digraphE<E>::digraphE( 
     const digraphE& g, const equiv_rel& e ) : edges_( g.edges_ )
{    Initialize( g, e );    }

template<class E> digraphE<E>::digraphE( const vec<digraphE>& g )
{    Initialize(g);    }

template<class E> void digraphE<E>::Initialize( const vec<digraphE>& g )
{    for ( int i = 0; i < g.isize( ); i++ )
          Append( g[i] );    }

template<class F> void digraphE<F>::Initialize( const vec<digraphE>& g,
     const vec< pair< pair<int,int>, pair<int,int> > >& joins )
{    digraphE<F> G(g);
     equiv_rel e( G.N( ) );
     vec<int> start( g.isize( ) );
     start[0] = 0;
     for ( int i = 1; i < g.isize( ); i++ )
          start[i] = start[i-1] + g[i-1].N( );
     for ( int i = 0; i < joins.isize( ); i++ )
     {    int v = start[ joins[i].first.first ] + joins[i].first.second;
          int w = start[ joins[i].second.first ] + joins[i].second.second;
          e.Join( v, w );    }
     Initialize( G, e );    }

template<class F> digraphE<F>::digraphE( const vec<digraphE>& g,
     const vec< pair< pair<int,int>, pair<int,int> > >& joins )
{    Initialize( g, joins );    }


template<class F> void digraphE<F>::Initialize( const digraph& g, const vec<F>& edges ){    
  int nedges = g.N();
  ForceAssertEq( nedges, edges.isize() );
  equiv_rel e( 2*nedges );
  for ( int v = 0; v < nedges; v++ ){    
    for ( size_t j = 0; j < g.From(v).size( ); j++ ){    
      int w = g.From(v)[j];
      e.Join( 2*v + 1, 2*w );    
    }    
  }
  vec<int> reps;
  e.OrbitRepsAlt(reps);
 
  int N = reps.size( );
  vec< vec<int> > from(N), to(N);
  vec< vec<int> > from_edge_obj(N), to_edge_obj(N);
  for ( int i = 0; i < edges.isize( ); i++ ){    
    int x = BinPosition( reps, e.ClassId( 2*i ) );
    int y = BinPosition( reps, e.ClassId( 2*i + 1 ) );
    from[x].push_back(y), to[y].push_back(x);    
    from_edge_obj[x].push_back(i), to_edge_obj[y].push_back(i);    
  }
  
  for ( int i = 0; i < N; i++ ){
    SortSync( from[i], from_edge_obj[i] );
    SortSync( to[i], to_edge_obj[i] );
  }    

  Initialize( from, to, edges, to_edge_obj, from_edge_obj );
}

template<class F> digraphE<F>::digraphE( const digraph& g, const vec<F>& edges ){
  Initialize( g, edges );
}


template<class E> digraphE<E>::digraphE( const digraph& g ){
  vec<int> edges( g.N(), vec<int>::IDENTITY );
  Initialize( g, edges );
}


template<class F> Bool digraphE<F>::ThisClose( int v, int w, F d ) const
{    if ( d < 0 ) return False;
     if ( v == w ) return True;
     set< pair<int,F> > unprocessed, processed;

     unprocessed.insert( make_pair( v, 0 ) );
     while( !unprocessed.empty( ) )
     {    int x = unprocessed.begin( )->first;
          F dx = unprocessed.begin( )->second;
          typename set< pair<int,F> >::iterator u 
               = processed.lower_bound( make_pair( x, 0 ) );
          unprocessed.erase( unprocessed.begin( ) );
          if ( u != processed.end( ) && u->first == x )
          {    if ( u->second <= dx ) continue;
               processed.erase(u);    }
          processed.insert( make_pair( x, dx ) );
          for ( int j = 0; j < From(x).isize( ); j++ )
          {    int y = From(x)[j];
               F dy = dx + EdgeObjectByIndexFrom( x, j );
               if ( dy > d ) continue;
               if ( y == w ) return True;
               typename set< pair<int,F> >::iterator p 
                    = processed.lower_bound( make_pair( y, 0 ) );
               if ( p != processed.end( ) && p->first == y )
               {    if ( p->second <= dy ) continue;
                    processed.erase(p);    }
               typename set< pair<int,F> >::iterator u 
                    = unprocessed.lower_bound( make_pair( y, 0 ) );
               if ( u != unprocessed.end( ) && u->first == y )
               {    if ( u->second <= dy ) continue;
                    unprocessed.erase(u);    }
               unprocessed.insert( make_pair( y, dy ) );    }    }
     return False;    }

template<class E> void digraphE<E>::ToLeft( vec<int>& to_left ) const
{    to_left.resize( EdgeObjectCount( ), -1 );
     for ( int i = 0; i < N( ); i++ )
     {    for ( int j = 0; j < From(i).isize( ); j++ )
          {    int e = EdgeObjectIndexByIndexFrom( i, j );
               to_left[e] = i;    }    }    }

template<class E> void digraphE<E>::ToRight( vec<int>& to_right ) const
{    to_right.resize( EdgeObjectCount( ), -1 );
     for ( int i = 0; i < N( ); i++ )
     {    for ( int j = 0; j < To(i).isize( ); j++ )
          {    int e = EdgeObjectIndexByIndexTo( i, j );
               to_right[e] = i;    }    }    }

template<class E> void digraphE<E>::ToLeftParallel( vec<int>& to_left ) const
{    to_left.resize( EdgeObjectCount( ) );
     #pragma omp parallel for
     for ( int i = 0; i < N( ); i++ )
     {    for ( int j = 0; j < From(i).isize( ); j++ )
          {    int e = EdgeObjectIndexByIndexFrom( i, j );
               to_left[e] = i;    }    }    }

template<class E> void digraphE<E>::ToRightParallel( vec<int>& to_right ) const
{    to_right.resize( EdgeObjectCount( ) );
     #pragma omp parallel for
     for ( int i = 0; i < N( ); i++ )
     {    for ( int j = 0; j < To(i).isize( ); j++ )
          {    int e = EdgeObjectIndexByIndexTo( i, j );
               to_right[e] = i;    }    }    }

template<class E> void digraphE<E>::ValidateLR( 
     const vec<int>& to_left, const vec<int>& to_right ) const
{    vec<int> to_leftx, to_rightx;
     ToLeft(to_leftx), ToRight(to_rightx);
     ForceAssert( to_left == to_leftx );
     ForceAssert( to_right == to_rightx );    }

template<class E> void digraphE<E>::ValidateLRParallel( 
     const vec<int>& to_left, const vec<int>& to_right ) const
{    ForceAssertEq( to_left.isize( ), E( ) );
     ForceAssertEq( to_right.isize( ), E( ) );
     #pragma omp parallel for
     for ( int v = 0; v < N( ); v++ )
     for ( int j = 0; j < From(v).isize( ); j++ )
     {    int w = From(v)[j], e = IFrom(v,j);
          ForceAssertEq( to_left[e], v );
          ForceAssertEq( to_right[e], w );    }    }

template<class F>
void digraphE<F>::GetSuccessors( const vec<int>& v, vec< pair<int,F> >& from_v )
{    set< pair<int,F> > check, fromv;

     for ( int i = 0; i < v.isize( ); i++ )
          check.insert( make_pair( v[i], 0 ) );
     while( !check.empty( ) )
     {    int x = check.begin( )->first;
          F dx = check.begin( )->second;
          typename set< pair<int,F> >::iterator u 
               = fromv.lower_bound( make_pair( x, 0 ) );
          check.erase( check.begin( ) );
          if ( u != fromv.end( ) && u->first == x )
          {    if ( u->second <= x ) continue;
               fromv.erase(u);    }
          fromv.insert( make_pair( x, dx ) );
          for ( int i = 0; i < From(x).isize( ); i++ )
          {    int y = From(x)[i]; 
               F dy = dx + EdgeObjectByIndexFrom( x, i );
               typename set< pair<int,F> >::iterator a 
                    = check.lower_bound( make_pair( y, 0 ) );
               if ( a != check.end( ) && a->first == y )
               {    if ( a->second <= dy ) continue;
                    check.erase(a);    }
               typename set< pair<int,F> >::iterator b 
                    = fromv.lower_bound( make_pair( y, 0 ) );
               if ( b != fromv.end( ) && b->first == y )
               {    if ( b->second <= dy ) continue;
                    fromv.erase(b);    }
               check.insert( make_pair( y, dy ) );    }    }
     from_v.clear( );
     for ( typename set< pair<int,F> >::iterator i = fromv.begin( ); 
          i != fromv.end( ); ++i )
     {    from_v.push_back(*i);    }    }

template<class F>
void PrintEdge( const int v, const int w, const int ei, const vec<double>& lengths,
     const vec<Bool>* dashed, const vec<String>* edge_color, const int tiny_top,
     const typename digraphE<F>::edge_label_info eli, ostream& out )
{
     float wd = 0.1; // this value not used
     String color, label;
     Bool bold = False;
     Bool is_dashed = False;
     if ( dashed != NULL ) is_dashed = (*dashed)[ei];
     double len = lengths[ei];
     if ( len < tiny_top ) 
     {    color = "gray";
	  if ( v == w ) label = ToString( len, 0 );
	  wd = 1.0;    }
     else if ( len >= tiny_top && len < 1000.0 ) 
     {    color = "black";
	  wd = 2.0;    }
     else if ( len >= 1000.0 && len < 10000.0 ) 
     {    color = "red";
          wd = 4.0;
          label = ToString( len/1000.0, 1 ) + " kb";    }
     else 
     {    color = "magenta";
	  bold = True;
	  wd = 8.0;
	  label = ToString( len/1000.0, 0 ) + " kb";    }
     if ( edge_color != NULL && (*edge_color)[ei] != "" ) color = (*edge_color)[ei];
     out << v << " -> " << w << " [minlen=" << wd << ",color=" << color;
     if ( color == "brown" ) out << ",penwidth=4";
     if (is_dashed) out << ",style=dashed";
     else if (bold) out << ",style=bold";
     if ( eli.edge_id_names != NULL ) 
     {    if ( label == "" ) label = (*eli.edge_id_names)[ei];
          else
          {    label = (*eli.edge_id_names)[ei] + " (" + label + ")";    }    }
     else if ( eli.label_edges ) 
     {    if ( label == "" ) 
               label = ( eli.edge_labels_base_alpha ? BaseAlpha(ei) : ToString(ei) );
          else
          {    label = ( eli.edge_labels_base_alpha ? BaseAlpha(ei) 
                    : ToString(ei) ) + " (" + label + ")";    }    }
     if ( eli.label_edges_extra ) label += " " + (*eli.label_edges_extra)[ei];
     if ( label != "" ) out << ",label=\"" << label << "\"";    }

template<class F>
void PrintEdge2( const int v, const int w, const int ei, const vec<double>& lengths,
     const vec<Bool>* dashed, const vec<String>* edge_color, const int tiny_top,
     const typename digraphEX<F>::edge_label_info eli, ostream& out )
{
     float wd = 0.1; // this value not used
     String color, label;
     Bool bold = False;
     Bool is_dashed = False;
     if ( dashed != NULL ) is_dashed = (*dashed)[ei];
     double len = lengths[ei];
     if ( len < tiny_top ) 
     {    color = "gray";
	  if ( v == w ) label = ToString( len, 0 );
	  wd = 1.0;    }
     else if ( len >= tiny_top && len < 1000.0 ) 
     {    color = "black";
	  wd = 2.0;    }
     else if ( len >= 1000.0 && len < 10000.0 ) 
     {    color = "red";
          wd = 4.0;
          label = ToString( len/1000.0, 1 ) + " kb";    }
     else 
     {    color = "magenta";
	  bold = True;
	  wd = 8.0;
	  label = ToString( len/1000.0, 0 ) + " kb";    }
     if ( edge_color != NULL && (*edge_color)[ei] != "" ) color = (*edge_color)[ei];
     out << v << " -> " << w << " [minlen=" << wd << ",color=" << color;
     if ( color == "brown" ) out << ",penwidth=4";
     if (is_dashed) out << ",style=dashed";
     else if (bold) out << ",style=bold";
     if ( eli.edge_id_names != NULL ) 
     {    if ( label == "" ) label = (*eli.edge_id_names)[ei];
          else
          {    label = (*eli.edge_id_names)[ei] + " (" + label + ")";    }    }
     else if ( eli.label_edges ) 
     {    if ( label == "" ) 
               label = ( eli.edge_labels_base_alpha ? BaseAlpha(ei) : ToString(ei) );
          else
          {    label = ( eli.edge_labels_base_alpha ? BaseAlpha(ei) 
                    : ToString(ei) ) + " (" + label + ")";    }    }
     if ( eli.label_edges_extra ) label += " " + (*eli.label_edges_extra)[ei];
     if ( label != "" ) out << ",label=\"" << label << "\"";    }

template<class E> void FindLeftMostVertex( const digraphE<E>& G,
     const vec<double>& lengths, const vec<int>& o, const vec<Bool>* invisible, 
     int& leftv )
{
    // Restrict attention to visible vertices.

    vec<int> oo;
    for ( int i1 = 0; i1 < o.isize( ); i1++ ) 
    {    int v = o[i1];
         if ( invisible != NULL )
         {    Bool owned = False;
              for ( int j = 0; j < G.From(v).isize( ); j++ ) 
              {    int e = G.EdgeObjectIndexByIndexFrom( v, j );
                   if ( !(*invisible)[e] ) owned = True;    }
              for ( int j = 0; j < G.To(v).isize( ); j++ ) 
              {    int e = G.EdgeObjectIndexByIndexTo( v, j );
                   if ( !(*invisible)[e] ) owned = True;    }
              if ( !owned ) continue;    }
         oo.push_back(v);    }
    Sort(oo);

    vec<float> pos( oo.size( ) );
    vec<Bool> placed( oo.size( ), False );
    pos[0] = 0.0, placed[0] = True;

    // Note that the following block of code is very slow on large components.
    // Performance may be OK now.

    while( Sum(placed) < oo.isize( ) ) 
    {    Bool progress = False;
         for ( int i1 = 0; i1 < oo.isize( ); i1++ ) 
         {    int v = oo[i1];
	      for ( int j = 0; j < G.From(v).isize( ); j++ ) 
              {    int w = G.From(v)[j];
	           int i2 = BinPosition( oo, w );
                   if ( i2 < 0 ) continue;
	           if ( !( placed[i1] ^ placed[i2] ) ) continue;
                   progress = True;
	           edge_t e = G.EdgeObjectIndexByIndexFrom( v, j );
	           if ( placed[i1] ) pos[i2] = pos[i1] + lengths[e];
	           else pos[i1] = pos[i2] - lengths[e];
	           placed[i1] = placed[i2] = True;    }    }
         if ( !progress ) break;    }
    
    float left = Min(pos);
    int leftj = 0;
    for ( leftj = 0; leftj < pos.isize( ); leftj++ )
      if ( pos[leftj] == left ) break;
    leftv = oo[leftj];    }

template<class E> void FindLeftMostVertex( const digraphEX<E>& G,
     const vec<double>& lengths, const vec<int>& o, const vec<Bool>* invisible, 
     int& leftv )
{
    // Restrict attention to visible vertices.

    vec<int> oo;
    for ( int i1 = 0; i1 < o.isize( ); i1++ ) 
    {    int v = o[i1];
         if ( invisible != NULL )
         {    Bool owned = False;
              for ( int j = 0; j < (int) G.From(v).size( ); j++ ) 
              {    int e = G.IFrom( v, j );
                   if ( !(*invisible)[e] ) owned = True;    }
              for ( int j = 0; j < (int) G.To(v).size( ); j++ ) 
              {    int e = G.ITo( v, j );
                   if ( !(*invisible)[e] ) owned = True;    }
              if ( !owned ) continue;    }
         oo.push_back(v);    }
    Sort(oo);

    vec<float> pos( oo.size( ) );
    vec<Bool> placed( oo.size( ), False );
    pos[0] = 0.0, placed[0] = True;

    // Note that the following block of code is very slow on large components.
    // Performance may be OK now.

    while( Sum(placed) < oo.isize( ) ) 
    {    Bool progress = False;
         for ( int i1 = 0; i1 < oo.isize( ); i1++ ) 
         {    int v = oo[i1];
	      for ( int j = 0; j < (int) G.From(v).size( ); j++ ) 
              {    int w = G.From(v)[j];
	           int i2 = BinPosition( oo, w );
                   if ( i2 < 0 ) continue;
	           if ( !( placed[i1] ^ placed[i2] ) ) continue;
                   progress = True;
	           edge_t e = G.IFrom( v, j );
	           if ( placed[i1] ) pos[i2] = pos[i1] + lengths[e];
	           else pos[i1] = pos[i2] - lengths[e];
	           placed[i1] = placed[i2] = True;    }    }
         if ( !progress ) break;    }
    
    float left = Min(pos);
    int leftj = 0;
    for ( leftj = 0; leftj < pos.isize( ); leftj++ )
      if ( pos[leftj] == left ) break;
    leftv = oo[leftj];    }

template<class E> void LabelTransitionVertices( const digraphE<E>& G,
     const int v, const vec<Bool>* invisible, ostream& out )
{    int vis_count = 0, invis_count = 0;
     for ( int j = 0; j < G.From(v).isize( ); j++ )
     {    int ei = G.EdgeObjectIndexByIndexFrom( v, j );
          if ( (*invisible)[ei] ) invis_count++;
          else vis_count++;    }
     for ( int j = 0; j < G.To(v).isize( ); j++ )
     {    int ei = G.EdgeObjectIndexByIndexTo( v, j );
          if ( (*invisible)[ei] ) invis_count++;
          else vis_count++;    }
     if ( vis_count > 0 && invis_count > 0 ) out << v << " [color=red];\n";    }

template<class E> void LabelTransitionVertices( const digraphEX<E>& G,
     const int v, const vec<Bool>* invisible, ostream& out )
{    int vis_count = 0, invis_count = 0;
     for ( int j = 0; j < (int) G.From(v).size( ); j++ )
     {    int ei = G.IFrom( v, j );
          if ( (*invisible)[ei] ) invis_count++;
          else vis_count++;    }
     for ( int j = 0; j < (int) G.To(v).size( ); j++ )
     {    int ei = G.ITo( v, j );
          if ( (*invisible)[ei] ) invis_count++;
          else vis_count++;    }
     if ( vis_count > 0 && invis_count > 0 ) out << v << " [color=red];\n";    }

template<class E> void CreateContigLabels( const vec<vec<int>>& components,
     const vec<String>* label_contigs_extra, vec<String>& contig_labels0,
     vec<String>& contig_labels )
{    vec<int> label_distance;
     if ( label_contigs_extra ) contig_labels0 = *label_contigs_extra;
     else 
     {    contig_labels0.resize( components.size( ) );
          for (size_t ii=0; ii<components.size( ); ii++)
               contig_labels0[ii] = "contig " + ToString( ii );    }
     label_distance.resize( contig_labels0.size( ), 0 );
     for (int ii=0; ii<(int)contig_labels0.size( ); ii++)
          label_distance[ii] = 1 + (int)( contig_labels0[ii].size( ) / 2 );
     for ( int i = 0; i < contig_labels0.isize( ); i++ )
     {    contig_labels.push_back( ",taillabel=\"" + contig_labels0[i] 
               + "\",labelangle=180," + "weight=10000," + "labeldistance=" 
               + ToString(label_distance[i]) + ",labelfontsize=18," 
               + "labelfontname=\"Times-Bold\"" );    }    }

template<class E>
void DotHeader( const Bool label_contigs, const Bool label_vertices,
     const String layout, const double fontsize, const double scale, 
     const double aspect, ostream& out )
{    out << "digraph G {\n\n";
     if (label_vertices)
     {    out << "node [width=" << scale * 0.1 << ",height=" << scale * 0.1 
               << ",fontsize=12,shape=plaintext];\n";    }
     else 
     {    out << "node [width=" << scale * 0.1
               << ",height=" << scale * 0.1 << ",fontsize=10,shape=point];\n";    }
     out << "edge [fontsize=" << fontsize 
          << ",penwidth=" << scale * 1.0
          << ",arrowsize=" << scale * 1.0
          << ",fontname=Arial];\n";
     if (label_contigs) out << "margin=1.0;\n";
     out << "rankdir=LR;\n";
     out << "labeljust=l;\n";
     out << "margin=0;\n";
     if ( layout != "" ) out << "layout=" << layout << ";\n";
     if ( aspect >= 0 ) out << "ratio=" << aspect << ";\n";    }

template<class F> void
digraphE<F>::PrettyDOT( ostream& out, const vec<double>& lengths,
     const edge_label_info eli, Bool label_contigs, Bool label_vertices,
     const vec<int>* componentsToPrint, const vec<String> *label_contigs_extra,
     const vec<int> *verticesToPrint, const vec<Bool>* dashed,
     const vec<Bool>* invisible, const vec<String>* edge_color,
     const vec<int>* pen_widths, const String layout, const double tiny_top,
     const double fontsize, const double scale, const double aspect ) const
{
     // Define components and those that are selected.

     vec< vec<int> > components;
     if ( invisible == NULL ) Components(components);
     else
     {    vec<int> to_left, to_right;
          ToLeft(to_left), ToRight(to_right);
          vec<Bool> invis( N( ), True );
          for ( int e = 0; e < EdgeObjectCount( ); e++ )
          {    if ( !(*invisible)[e] )
                    invis[ to_left[e] ] = invis[ to_right[e] ] = False;    }
          Components( components, &invis );    }
     vec<int> select;
     if (componentsToPrint) select = *componentsToPrint;
     else 
     {    select.reserve( components.size( ) );
          for (int ii=0; ii<(int)components.size( ); ii++) 
               select.push_back( ii );    }

     // Set up output and contig labels.

     DotHeader<F>( 
          label_contigs, label_vertices, layout, fontsize, scale, aspect, out );
     vec<String> contig_labels0, contig_labels;
     if (label_contigs) 
     {    CreateContigLabels<F>( components, label_contigs_extra, 
               contig_labels0, contig_labels );    }

     // Define vertices to skip.

     vec<bool> skip_vtx;
     if (verticesToPrint) 
     {    skip_vtx.resize( this->N( ), true );
          for (size_t ii=0; ii<verticesToPrint->size( ); ii++)
               skip_vtx[ (*verticesToPrint)[ii] ] = false;    }
  
      // Print the contigs.  We put each contig in its own cluster (the
      // subgraph's name MUST start with "cluster" for this to have any effect).

      for ( int sel_id = select.isize( ) - 1; sel_id >= 0; sel_id-- ) 
      {    int i = select[sel_id];
           vec<int> &o = components[i];
           if ( invisible != NULL )
           {    int vis_count = 0;
                for ( int vi = 0; vi < o.isize( ); vi++ ) 
                {    int v = o[vi];
                     if ( verticesToPrint && skip_vtx[v] ) continue;
                     for ( int j = 0; j < From(v).isize( ); j++ )
                     {    int ei = EdgeObjectIndexByIndexFrom( v, j );
                          if ( !(*invisible)[ei] ) vis_count++;    }    }
                if ( vis_count == 0 ) continue;    }

          out << "\nsubgraph cluster" << i << " {\n";
          out << "color=white;\n";
          if ( label_contigs && label_contigs_extra )
          {    out << "label=\"" << contig_labels0[i] << "\"," << "fontsize=18,"
	            << "fontname=\"Times-Bold\"\n";    }
    
          // Find "leftmost" vertex in graph.

          Sort(o);
          int leftv;
          FindLeftMostVertex( *this, lengths, o, invisible, leftv );

          // Print component.
    
          for ( int vi = 0; vi < o.isize( ); vi++ ) 
          {    int v = o[vi];
               if ( verticesToPrint && skip_vtx[v] ) continue;
               if (label_vertices)
	       {    out << v << " [label=" << "\"" << v << "\"" 
	                 << ",fontcolor=black];\n";    }

               // If some edges touching a vertex are invisible, and some are 
               // visible, make vertex red.  Note incompatibility with 
               // label_vertices.

               else if ( invisible != NULL ) 
               {    LabelTransitionVertices( *this, v, invisible, out );    }

               for ( int j = 0; j < From(v).isize( ); j++ ) 
               {    int ei = EdgeObjectIndexByIndexFrom( v, j );
                    if ( invisible != NULL && (*invisible)[ei] ) continue;
	            int w = From(v)[j];
                    PrintEdge<F>( 
                         v, w, ei, lengths, dashed, edge_color, tiny_top, eli, out );
	            if ( label_contigs && v == leftv && j == 0 
                         && !label_contigs_extra )
	            {    out << contig_labels[i];    }
                    if ( pen_widths != NULL && (*pen_widths)[ei] > 0 )
                         out << ",penwidth=" << (*pen_widths)[ei];
	            out << "];\n";    }    }
          out << "}\n";    }
     out << "\n}" << endl;
     out << "#done" << endl;    }

template<class F> void
digraphEX<F>::PrettyDOT( ostream& out, const vec<double>& lengths,
     const edge_label_info eli, Bool label_contigs, Bool label_vertices,
     const vec<int>* componentsToPrint, const vec<String> *label_contigs_extra,
     const vec<int> *verticesToPrint, const vec<Bool>* dashed,
     const vec<Bool>* invisible, const vec<String>* edge_color,
     const vec<int>* pen_widths, const String layout, const double tiny_top,
     const double fontsize, const double scale, const double aspect ) const
{
     // Define components and those that are selected.

     vec< vec<int> > components;
     if ( invisible == NULL ) Components(components);
     else
     {    vec<Bool> invis( N( ), True );
          for ( int e = 0; e < E( ); e++ )
          {    if ( !(*invisible)[e] )
                    invis[ ToLeft(e) ] = invis[ ToRight(e) ] = False;    }
          Components( components, &invis );    }
     vec<int> select;
     if (componentsToPrint) select = *componentsToPrint;
     else 
     {    select.reserve( components.size( ) );
          for (int ii=0; ii<(int)components.size( ); ii++) 
               select.push_back( ii );    }

     // Set up output and contig labels.

     DotHeader<F>( 
          label_contigs, label_vertices, layout, fontsize, scale, aspect, out );
     vec<String> contig_labels0, contig_labels;
     if (label_contigs) 
     {    CreateContigLabels<F>( components, label_contigs_extra, 
               contig_labels0, contig_labels );    }

     // Define vertices to skip.

     vec<bool> skip_vtx;
     if (verticesToPrint) 
     {    skip_vtx.resize( this->N( ), true );
          for (size_t ii=0; ii<verticesToPrint->size( ); ii++)
               skip_vtx[ (*verticesToPrint)[ii] ] = false;    }
  
      // Print the contigs.  We put each contig in its own cluster (the
      // subgraph's name MUST start with "cluster" for this to have any effect).

      for ( int sel_id = select.isize( ) - 1; sel_id >= 0; sel_id-- ) 
      {    int i = select[sel_id];
           vec<int> &o = components[i];
           if ( invisible != NULL )
           {    int vis_count = 0;
                for ( int vi = 0; vi < o.isize( ); vi++ ) 
                {    int v = o[vi];
                     if ( verticesToPrint && skip_vtx[v] ) continue;
                     for ( int j = 0; j < (int) From(v).size( ); j++ )
                     {    int ei = IFrom( v, j );
                          if ( !(*invisible)[ei] ) vis_count++;    }    }
                if ( vis_count == 0 ) continue;    }

          out << "\nsubgraph cluster" << i << " {\n";
          out << "color=white;\n";
          if ( label_contigs && label_contigs_extra )
          {    out << "label=\"" << contig_labels0[i] << "\"," << "fontsize=18,"
	            << "fontname=\"Times-Bold\"\n";    }
    
          // Find "leftmost" vertex in graph.

          Sort(o);
          int leftv;
          FindLeftMostVertex( *this, lengths, o, invisible, leftv );

          // Print component.
    
          for ( int vi = 0; vi < o.isize( ); vi++ ) 
          {    int v = o[vi];
               if ( verticesToPrint && skip_vtx[v] ) continue;
               if (label_vertices)
	       {    out << v << " [label=" << "\"" << v << "\"" 
	                 << ",fontcolor=black];\n";    }

               // If some edges touching a vertex are invisible, and some are 
               // visible, make vertex red.  Note incompatibility with 
               // label_vertices.

               else if ( invisible != NULL ) 
               {    LabelTransitionVertices( *this, v, invisible, out );    }

               for ( int j = 0; j < (int) From(v).size( ); j++ ) 
               {    int ei = IFrom( v, j );
                    if ( invisible != NULL && (*invisible)[ei] ) continue;
	            int w = From(v)[j];
                    PrintEdge2<F>( 
                         v, w, ei, lengths, dashed, edge_color, tiny_top, eli, out );
	            if ( label_contigs && v == leftv && j == 0 
                         && !label_contigs_extra )
	            {    out << contig_labels[i];    }
                    if ( pen_widths != NULL && (*pen_widths)[ei] > 0 )
                         out << ",penwidth=" << (*pen_widths)[ei];
	            out << "];\n";    }    }
          out << "}\n";    }
     out << "\n}" << endl;
     out << "#done" << endl;    }

// Method: DumpGraphML
// Output the digraph structure in a textual format that can be easily
// read without reference to our code base.

template<class E>
void
digraphE<E>::DumpGraphML( const String& graphMLFileName ) const
{
  vec< vec< String > > edgeLabels( N() );
  for ( int v = 0; v < N( ); v++ ) {
    for ( int j = 0; j < From(v).isize( ); j++ ) {
      int w = From(v)[j];
      edgeLabels[ v ].push_back( BaseAlpha( EdgeObjectIndexByIndexFrom( v, j ) ) );
    }
  }
  
  Ofstream( grml, graphMLFileName );
  WriteGraphML( grml, edgeLabels );
}

template<class E> void digraphE<E>::ComponentsE( vec< vec<int> >& comp ) const
{    comp.clear( );
     equiv_rel e( N( ) );
     for ( int v = 0; v < N( ); v++ )
     {    for ( int j = 0; j < From(v).isize( ); j++ )
               e.Join( v, From(v)[j] );    }
     for ( int x = 0; x < N( ); x++ )
     {    if ( e.Representative(x) )
          {    vec<int> o;
               e.Orbit( x, o );
               Sort(o);
               vec<int> C;
               for ( int i = 0; i < o.isize( ); i++ )
               {    int v = o[i];
                    for ( int j = 0; j < From(v).isize( ); j++ )
                         C.push_back( EdgeObjectIndexByIndexFrom( v, j ) );    }
               comp.push_back(C);    }    }    }

template<class E> void digraphE<E>::ComponentsEFast( vec< vec<int> >& comp ) const
{    comp.clear( );
     Components(comp);
     for ( int64_t j = 0; j < comp.jsize( ); j++ )
     {    vec<int>& o = comp[j];
          vec<int> C;
          for ( int i = 0; i < o.isize( ); i++ )
          {    int v = o[i];
               for ( int j = 0; j < From(v).isize( ); j++ )
                    C.push_back( EdgeObjectIndexByIndexFrom( v, j ) );    }
          o = C;    }    }

template<class E> void LongestPath( const digraphE<E>& G, int (E::*len)( ) const,
     vec<int>& a_longest_path )
{    vec<int> D;
     const int infinity = 2000000000;
     DistancesToEnd( G, len, infinity, True, D );
     int M = 0, v = 0;
     for ( int x = 0; x < G.N( ); x++ )
          if ( D[x] > M ) { M = D[x], v = x; }
     a_longest_path.clear( );
     while( G.From(v).nonempty( ) )
     {    for ( int j = 0; j < G.From(v).isize( ); j++ )
          {    int w = G.From(v)[j];
               if ( D[w] == D[v] - ((G.EdgeObjectByIndexFrom( v, j )).*len)( ) )
               {    a_longest_path.push_back( G.EdgeObjectIndexByIndexFrom( v, j ) );
                    v = w;
                    break;    }    }    }    }

template<class E> void DistancesToEndArr( const digraphE<E>& G,
     vec<int> const& edgeLens, const int max_dist, const Bool fw, vec<int>& D )
{
     // Let D(v) be the maximum length of a path starting at v, to be computed.
     // Define initial values for D(v) to be 'undefined', except for sinks, 
     // which are zero.

     D.resize_and_set( G.N( ), -1 );
     for ( int v = 0; v < G.N( ); v++ )
     {    if ( fw && G.Sink(v) ) D[v] = 0;
          if ( !fw && G.Source(v) ) D[v] = 0;    }

     // Initialize vertices to process.

     vec<Bool> to_process( G.N( ), False );
     vec<int> to_processx;
     for ( int v = 0; v < G.N( ); v++ )
     {    if ( (fw && G.Sink(v)) || (!fw && G.Source(v)) ) 
          {    to_process[v] = True, to_processx.push_back(v);    }    }

     // Now compute D.  Uncomputed values are set to 'infinity'.

     while( to_processx.nonempty( ) )
     {    int v = to_processx.back( );
          to_processx.pop_back( );
          to_process[v] = False;
          for ( int j = 0; j < (fw ? G.To(v) : G.From(v) ).isize( ); j++ )
          {    int w = ( fw ? G.To(v) : G.From(v) )[j];
               if ( D[w] >= max_dist ) continue;
               int edgeId = ( fw ? G.EdgeObjectIndexByIndexTo(v, j)
                    : G.EdgeObjectIndexByIndexFrom(v, j) );
               int Dw_new = edgeLens[edgeId] + D[v];
               if ( Dw_new > D[w] ) 
               {    D[w] = Dw_new;
                    if ( !to_process[w] )
                    {    to_process[w] = True; 
                         to_processx.push_back(w);    }    }    }    }
     for ( int v = 0; v < G.N( ); v++ )
          if ( D[v] < 0 ) D[v] = max_dist;    }

template<class E> void RemoveHangingEnds( digraphE<E>& G, 
     int (E::*len)( ) const, const int max_del, const double min_ratio )
{
     // Track hanging ends.

     vec<Bool> hanging( G.EdgeObjectCount( ), False );

     // Define the maximum length that we care about.

     const int max_dist = int( ceil( double(max_del) * min_ratio ) );

     // Go through two passes (forward and reverse).

     for ( int pass = 1; pass <= 2; pass++ )
     {
          // Compute distances to end.

          vec<int> D;
          DistancesToEnd( G, len, max_dist, pass == 1, D );

          // Identify hanging ends.

          for ( int v = 0; v < G.N( ); v++ )
          {    const vec<int>& V = ( pass == 1 ? G.From(v) : G.To(v) );
               vec<int> d( V.size( ) ); 
               vec<int> id( V.size( ), vec<int>::IDENTITY );
               for ( int j = 0; j < V.isize( ); j++ )
               {    d[j] = ((pass == 1 
                         ? G.EdgeObjectByIndexFrom(v,j) : G.EdgeObjectByIndexTo(v,j))
                         .*len)( ) + D[ V[j] ];    }
               ReverseSortSync( d, id );
               for ( int j = 1; j < d.isize( ); j++ )
               {    if ( d[j] <= max_del && d[0] >= d[j] * min_ratio )
                    {    hanging[ ( pass == 1 
                              ? G.EdgeObjectIndexByIndexFrom( v, id[j] )
                              : G.EdgeObjectIndexByIndexTo( v, id[j] ) ) ] 
                              = True;    }    }    }    }

     // Remove hanging ends.

     vec<int> to_delete;
     for ( int i = 0; i < G.EdgeObjectCount( ); i++ )
          if ( hanging[i] ) to_delete.push_back(i);
     G.DeleteEdges(to_delete);    }


// Remove short hanging ends.  Look for
//
//                 x
//                 |
//                 e
//                 |
//        u --c--> v --d--> w
//
// where x is a source or sink, e is short (and can go either way), whereas
// c and d are long.  Works for T = HyperKmerPath and T = HyperFastavector.

template<class T> void RemoveHangingEnds2( T& h,const int max_del,
     const double min_ratio )
{

  for ( int x = 0; x < h.N( ); x++ ) {
    
    // Check that basic assumptions are satisfied, including length(e) <= 5kb.
    
    int v, c, d, e;
    if ( h.Source(x) && h.From(x).size( ) == 1 ) {
      v = h.From(x)[0];
      e = h.EdgeObjectIndexByIndexFrom( x, 0 );
    } else if ( h.Sink(x) && h.To(x).size( ) == 1 ) {
      v = h.To(x)[0];
      e = h.EdgeObjectIndexByIndexTo( x, 0 );
    } else 
      continue;

    if ( h.EdgeLengthKmers(e) > max_del ) continue;

    if ( h.Source(x) ) {
      if ( !( h.From(v).size( ) == 1 && h.To(v).size( ) == 2 ) ) continue;
      d = h.EdgeObjectIndexByIndexFrom( v, 0 );
      c = h.EdgeObjectIndexByIndexTo( v, 0 );
      if ( c == e ) c = h.EdgeObjectIndexByIndexTo( v, 1 );
    } else {
      if ( !( h.From(v).size( ) == 2 && h.To(v).size( ) == 1 ) ) continue;
      c = h.EdgeObjectIndexByIndexTo( v, 0 );
      d = h.EdgeObjectIndexByIndexFrom( v, 0 );
      if ( d == e ) d = h.EdgeObjectIndexByIndexFrom( v, 1 );
    }

    // We require that there is an edge "competing with e", that is at least
    // 20 times longer.
    
    static vec<int> v_only(1), to_v, from_v;
    v_only[0] = v;
    int max_competitor = 0;
    if ( h.Source(x) ) {
      h.digraph::GetPredecessors( v_only, to_v );
      for ( int j = 0; j < to_v.isize( ); j++ ) {
	int z = to_v[j];
	for ( int i = 0; i < h.To(z).isize( ); i++ ) {
	  int e = h.EdgeObjectIndexByIndexTo( z, i );
	  max_competitor = Max( max_competitor, h.EdgeLengthKmers(e) );
	}
      }
    } else {
      h.digraph::GetSuccessors( v_only, from_v );
      for ( int j = 0; j < from_v.isize( ); j++ ) {
	int z = from_v[j];
	for ( int i = 0; i < h.From(z).isize( ); i++ ) {
	  int e = h.EdgeObjectIndexByIndexFrom( z, i );
	  max_competitor = Max( max_competitor, h.EdgeLengthKmers(e) );
	}
      }
    }

    if ( min_ratio * h.EdgeLengthKmers(e) > max_competitor ) continue;
    
    // Edit the graph.
    
    if ( h.Source(x) ) h.DeleteEdgeFrom( x, 0 );
    else h.DeleteEdgeTo( x, 0 );

  }
}


// Find the indices of all edges e that form self-loops, i.e., e goes from v -> v.
template<class E> vec<int> digraphE<E>::SelfLoops( ) const
{
  vec<int> to_left, to_right;
  ToLeft( to_left );
  ToRight( to_right );
  vec<Bool> used;
  Used( used );
  
  vec<int> self_loops;
  for ( int i = 0; i < EdgeObjectCount(); i++ )
    if ( to_left[i] == to_right[i] && used[i] )
      self_loops.push_back( i );
  
  return self_loops;
}




template<class E> void digraphE<E>::LoopSubgraph( vec<int>& loop_edges ) const
{    loop_edges.clear( );
     vec< vec<int> > SCC;
     StronglyConnectedComponents(SCC);
     for ( int i = 0; i < SCC.isize( ); i++ )
     {    const vec<int>& V = SCC[i];
          for ( int r = 0; r < V.isize( ); r++ )
          {    int v = V[r];
               for ( int j = 0; j < From(v).isize( ); j++ )
               {    if ( BinMember( V, From(v)[j] ) )
                    {    loop_edges.push_back( EdgeObjectIndexByIndexFrom( 
                              v, j ) );    }    }    }    }
     Sort(loop_edges);    }

template<class E> void digraphE<E>::SplayVertex( const int v )
{    int n = N( );
     AddVertices( To(v).size( ) );
     for ( int j = To(v).isize( ) - 1; j >= 0; j-- )
          GiveEdgeNewToVx( EdgeObjectIndexByIndexTo( v, j ), v, n + j );
     n = N( );
     AddVertices( From(v).size( ) );
     for ( int j = From(v).isize( ) - 1; j >= 0; j-- )
     {    GiveEdgeNewFromVx( EdgeObjectIndexByIndexFrom( v, j ), 
               v, n + j );    }    }

template<class E> void digraphE<E>::SplayVertexWithUpdate( const int v,
     vec<int>& to_left, vec<int>& to_right )
{    int n = N( );
     AddVertices( To(v).size( ) );
     for ( int j = To(v).isize( ) - 1; j >= 0; j-- )
     {    int e = ITo( v, j );
          GiveEdgeNewToVx( e, v, n + j );
          to_right[e] = n + j;    }
     n = N( );
     AddVertices( From(v).size( ) );
     for ( int j = From(v).isize( ) - 1; j >= 0; j-- )
     {    int e = IFrom( v, j );
          GiveEdgeNewFromVx( e, v, n + j );
          to_left[e] = n + j;    }    }

template<class E> void digraphE<E>::LiberateEdge( 
     const int e, const int v, const int w )
{    int j = EdgeObjectIndexToFromIndex( v, e );
     DeleteEdgeFrom( v, j );
     SplayVertex(v), SplayVertex(w);    }

template<class E> void digraphE<E>::GiveEdgeNewFromVx
( int edge_id, int old_from_v, int new_from_v ) {
       int i = Position( from_edge_obj_[old_from_v], edge_id );
       ForceAssert( i != -1 );
       int w = from_[old_from_v][i];
       int j = Position( to_edge_obj_[w],edge_id );
       ForceAssert( j != -1 );
       to_[w][j] = new_from_v;
       from_[old_from_v].erase( from_[old_from_v].begin() + i );
       from_edge_obj_[old_from_v].erase( from_edge_obj_[old_from_v].begin() + i );
       from_[new_from_v].push_back(w);
       from_edge_obj_[new_from_v].push_back(edge_id);
       SortSync( to_[w], to_edge_obj_[w] );
       SortSync( from_[new_from_v], from_edge_obj_[new_from_v] );
     }

template<class E> void digraphE<E>::GiveEdgeNewToVx
( int edge_id, int old_to_w, int new_to_w ) {
       int j = Position( to_edge_obj_[old_to_w], edge_id );
       if ( j == -1 )
       {    cout << "Edge does not stop at the claimed vertex." << endl;
            TracebackThisProcess( );    }
       int v = to_[old_to_w][j];
       int i = Position( from_edge_obj_[v],edge_id );
       ForceAssert( i != -1 );
       from_[v][i] = new_to_w;
       to_[old_to_w].erase( to_[old_to_w].begin() + j );
       to_edge_obj_[old_to_w].erase( to_edge_obj_[old_to_w].begin() + j );
       to_[new_to_w].push_back(v);
       to_edge_obj_[new_to_w].push_back(edge_id);
       SortSync( from_[v], from_edge_obj_[v] );
       SortSync( to_[new_to_w], to_edge_obj_[new_to_w] );
     }

template<class F> int digraphE<F>::AddEdge( int v, int w, const F& e )
{    int n = EdgeObjectCount( );
     edges_.push_back(e);
     int i = upper_bound( from_[v].begin(), from_[v].end(), w ) - from_[v].begin();
     from_[v].insert( from_[v].begin()+i, w );
     from_edge_obj_[v].insert( from_edge_obj_[v].begin()+i, n );
     int j = upper_bound( to_[w].begin(), to_[w].end(), v ) - to_[w].begin();
     to_[w].insert( to_[w].begin()+j, v );
     to_edge_obj_[w].insert( to_edge_obj_[w].begin()+j, n );
     return n;
}

template<class E>
Bool digraphE<E>::EdgePaths( const vec<int>& left, const vec<int>& right,
     const int v, const int w, vec< vec<int> >& paths, const int max_copies, 
     const int max_paths, const int max_iterations, 
     const Bool allow_terminal_loops ) const
{    
     // Pretest to determine if the computation will explode.  This only works if
     // max_copies is not set.

     if ( max_copies < 0 && ( max_paths >= 0 || max_iterations >= 0 ) )
     {    vec<int> subs;
          int path_count = 0;
          for ( int i = 0; i < From(v).isize( ); i++ )
          {    int e = IFrom( v, i );
               subs.push_back(e);    }
          int iterations = 0;
          while( subs.nonempty( ) )
          {    if ( max_iterations > 0 && ++iterations > max_iterations ) 
                    return False;
               int p = subs.back( );
               subs.pop_back( );
               int x = right[p];
               if ( x == w ) 
               {    if ( max_paths >= 0 && ++path_count > max_paths ) 
                         return False;    }
               if ( x == w && !allow_terminal_loops ) continue;
               for ( int j = 0; j < From(x).isize( ); j++ )
               {    int e = IFrom( x, j );
                    subs.push_back(e);    }    }    }

     // Now do the computation for real.

     vec< vec<int> > subs;
     paths.clear( );
     for ( int i = 0; i < From(v).isize( ); i++ )
     {    int e = IFrom( v, i );
          vec<int> one = {e};
          subs.push_back(one);    }
     int iterations = 0;
     while( subs.nonempty( ) )
     {    if ( max_iterations > 0 && ++iterations > max_iterations ) return False;
          vec<int> p = subs.back( );
          subs.pop_back( );
          int x = right[ p.back( ) ];
          if ( x == w ) 
          {    paths.push_back(p);
               if ( max_paths >= 0 && paths.isize( ) > max_paths ) return False;    }
          if ( x == w && !allow_terminal_loops ) continue;
          for ( int j = 0; j < From(x).isize( ); j++ )
          {    int e = IFrom( x, j );
               vec<int> pp(p);
               pp.push_back(e);
               if ( max_copies >= 0 )
               {    vec<int> pps(pp);
                    Sort(pps);
                    Bool fail = False;
                    for ( int r = 0; r < pps.isize( ); r++ )
                    {    int s = pps.NextDiff(r);
                         if ( s - r > max_copies )
                         {    fail = True;
                              break;    }
                         r = s - 1;    }
                    if (fail) continue;    }
               subs.push_back(pp);    }    }
     return True;    }

template<class E>
Bool digraphE<E>::EdgePathsLim( const vec<int>& left, const vec<int>& right,
     const int v, const int w, const int e_verboten, vec< vec<int> >& paths, 
     const int max_copies, const int max_paths, const int max_iterations ) const
{    
     vec< vec<int> > subs;
     paths.clear( );
     vec<int> one(1);
     for ( int i = 0; i < From(v).isize( ); i++ )
     {    int e = IFrom( v, i );
          if ( e == e_verboten ) continue;
          one[0] = e;
          subs.push_back(one);    }
     int iterations = 0;
     vec<int> p, pp, pps;
     while( subs.nonempty( ) )
     {    if ( max_iterations > 0 && ++iterations > max_iterations ) return False;
          p = subs.back( );
          subs.pop_back( );
          int x = right[ p.back( ) ];
          if ( x == w ) 
          {    paths.push_back(p);
               if ( max_paths >= 0 && paths.isize( ) > max_paths ) return False;    }
          for ( int j = 0; j < From(x).isize( ); j++ )
          {    int e = IFrom( x, j );
               if ( e == e_verboten ) continue;
               pp = p;
               pp.push_back(e);
               pps = pp;
               Sort(pps);
               Bool fail = False;
               for ( int r = 0; r < pps.isize( ); r++ )
               {    int s = pps.NextDiff(r);
                    if ( s - r > max_copies )
                    {    fail = True;
                         break;    }
                    r = s - 1;    }
               if ( !fail ) subs.push_back(pp);    }    }
     return True;    }

template<class E> Bool digraphEX<E>::EdgePaths( 
     const int v, const int w, vec< vec<int> >& paths, const int max_copies, 
     const int max_paths, const int max_iterations,
     const Bool allow_terminal_loops ) const
{    
     // Pretest to determine if the computation will explode.  This only works if
     // max_copies is not set.

     if ( max_copies < 0 && ( max_paths >= 0 || max_iterations >= 0 ) )
     {    vec<int> subs;
          int path_count = 0;
          for ( int i = 0; i < (int) From(v).size( ); i++ )
          {    int e = IFrom( v, i );
               subs.push_back(e);    }
          int iterations = 0;
          while( subs.nonempty( ) )
          {    if ( max_iterations > 0 && ++iterations > max_iterations ) 
                    return False;
               int p = subs.back( );
               subs.pop_back( );
               int x = ToRight(p);
               if ( x == w ) 
               {    if ( max_paths >= 0 && ++path_count > max_paths ) 
                         return False;    }
               if ( x == w && !allow_terminal_loops ) continue;
               for ( int j = 0; j < (int) From(x).size( ); j++ )
               {    int e = IFrom( x, j );
                    subs.push_back(e);    }    }    }

     // Now do the computation for real.

     vec< vec<int> > subs;
     paths.clear( );
     vec<int> one(1);
     for ( int i = 0; i < (int) From(v).size( ); i++ )
     {    int e = IFrom( v, i );
          one[0] = e;
          subs.push_back(one);    }
     int iterations = 0;
     vec<int> p, pp, pps;
     while( subs.nonempty( ) )
     {    if ( max_iterations > 0 && ++iterations > max_iterations ) return False;
          p = subs.back( );
          subs.pop_back( );
          int x = ToRight( p.back( ) );
          if ( x == w ) 
          {    paths.push_back(p);
               if ( max_paths >= 0 && paths.isize( ) > max_paths ) return False;    }
          if ( x == w && !allow_terminal_loops ) continue;
          for ( int j = 0; j < (int) From(x).size( ); j++ )
          {    int e = IFrom( x, j );
               pp = p;
               pp.push_back(e);
               if ( max_copies >= 0 )
               {    pps = pp;
                    Sort(pps);
                    Bool fail = False;
                    for ( int r = 0; r < pps.isize( ); r++ )
                    {    int s = pps.NextDiff(r);
                         if ( s - r > max_copies )
                         {    fail = True;
                              break;    }
                         r = s - 1;    }
                    if (fail) continue;    }
               subs.push_back(pp);    }    }
     return True;    }

template<class E>
Bool digraphE<E>::EdgePaths( const int v, const int w, vec< vec<int> >& paths,
     const int max_copies, const int max_paths, const int max_iterations,
     const Bool allow_terminal_loops ) const
{    vec<int> left, right;
     ToLeft(left), ToRight(right);
     return EdgePaths( left, right, v, w, paths, max_copies, max_paths, 
          max_iterations, allow_terminal_loops );    }

template<class E> void digraphE<E>::DeleteEdgeTo( int w, int j )
{    int v = to_[w][j];
     int i = InputToOutputFrom( w, j );
     to_[w].erase( to_[w].begin( ) + j );
     to_edge_obj_[w].erase( to_edge_obj_[w].begin( ) + j );
     from_[v].erase( from_[v].begin( ) + i );
     from_edge_obj_[v].erase( from_edge_obj_[v].begin( ) + i );    }

template<class E> void digraphE<E>::DeleteEdgeFrom( int v, int j )
{    int w = from_[v][j];
     int i = InputFromOutputTo( v, j );
     from_[v].erase( from_[v].begin( ) + j );
     from_edge_obj_[v].erase( from_edge_obj_[v].begin( ) + j );
     to_[w].erase( to_[w].begin( ) + i );
     to_edge_obj_[w].erase( to_edge_obj_[w].begin( ) + i );    }

template<class E> void digraphE<E>::DeleteEdgesTo( int w, const vec<int>& js )
{    for ( int l = js.isize( ) - 1; l >= 0; l-- )
          DeleteEdgeTo( w, js[l] );    }

template<class E> void digraphE<E>::DeleteEdgesFrom( int v, const vec<int>& js )
{    for ( int l = js.isize( ) - 1; l >= 0; l-- )
          DeleteEdgeFrom( v, js[l] );    }

template<class E> 
vec<int> digraphE<E>::EdgesBetween( const int v, const int w ) const
{    vec<int> b;
     for ( int i = 0; i < From(v).isize( ); i++ )
     {    if ( From(v)[i] == w )
               b.push_back( EdgeObjectIndexByIndexFrom( v, i ) );    }
     return b;    }

template<class E> vec<int> digraphE<E>::EdgesBetween( const vec<int>& v ) const
{    vec<int> b;
     for ( int j = 0; j < v.isize( ); j++ )
     {    for ( int i = 0; i < From(v[j]).isize( ); i++ )
          {    if ( BinMember( v, From(v[j])[i] ) )
                    b.push_back( EdgeObjectIndexByIndexFrom( v[j], i ) );    }    }
     Sort(b);
     return b;    }

template<class F> 
vec<F> digraphE<F>::EdgeObjectsBetween( const int v, const int w ) const
{    vec<F> b;
     for ( int i = 0; i < From(v).isize( ); i++ )
     {    if ( From(v)[i] == w )
               b.push_back( EdgeObjectByIndexFrom( v, i ) );    }
     return b;    }

template<class E> int digraphE<E>::InputToOutputFrom( int w, int i ) const
{    int v = to_[w][i];
     int ei = to_edge_obj_[w][i];
     for ( int j = 0; j < from_[v].isize( ); j++ )
          if ( from_edge_obj_[v][j] == ei ) return j;
     cout << "\nInputToOutputFrom failed." << endl;
     TracebackThisProcess( );
     return -1;    }

template<class E> int digraphE<E>::InputFromOutputTo( int w, int i ) const
{    int v = from_[w][i];
     int ei = from_edge_obj_[w][i];
     for ( int j = 0; j < to_[v].isize( ); j++ )
          if ( to_edge_obj_[v][j] == ei ) return j;
     cout << "\nInputFromOutputTo failed." << endl;
     TracebackThisProcess( );
     return -1;    }

template<class F> void digraphE<F>::ChangeEdgeObjectFrom( int v, int i, const F& e )
{    int ne = edges_.size( );
     edges_.push_back(e);
     int w = From(v)[i];
     int j = InputFromOutputTo( v, i );
     from_edge_obj_[v][i] = ne;
     to_edge_obj_[w][j] = ne;    }

template<class F> F digraphE<F>::MinEdge( int v, int w )
{    F m = 0;
     Bool first = True;
     for ( int j = 0; j < From(v).isize( ); j++ )
     {    if ( From(v)[j] != w ) continue;
          if (first) m = EdgeObjectByIndexFrom( v, j );
          else m = Min( m, EdgeObjectByIndexFrom( v, j ) );
          first = False;    }
     ForceAssert( !first );
     return m;    }

template<class F> F digraphE<F>::MaxEdge( int v, int w )
{    F M = 0;
     Bool first = True;
     for ( int j = 0; j < From(v).isize( ); j++ )
     {    if ( From(v)[j] != w ) continue;
          if (first) M = EdgeObjectByIndexFrom( v, j );
          else M = Max( M, EdgeObjectByIndexFrom( v, j ) );
          first = False;    }
     ForceAssert( !first );
     return M;    }

template<class F> void digraphE<F>::AddVertices( int nadd ) 
{    int nvert = N( );
     from_.resize( nvert + nadd );
     to_.resize( nvert + nadd );
     from_edge_obj_.resize( nvert + nadd );
     to_edge_obj_.resize( nvert + nadd );    }

template<class F> void digraphE<F>::DeleteEdges( const vec<int>& to_delete )
{    vec<int> to_delete_local;
     if ( !to_delete.UniqueOrdered( ) )
     {    to_delete_local = to_delete;
          UniqueSort(to_delete_local);    }
     const vec<int>& tod 
          = ( to_delete_local.nonempty( ) ? to_delete_local: to_delete );
     for ( int v = 0; v < N( ); v++ )
     {    for ( int j = From(v).isize( ) - 1; j >= 0; j-- )
          {    int e = EdgeObjectIndexByIndexFrom( v, j );
               if ( BinMember( tod, e ) ) DeleteEdgeFrom( v, j );    }    }    }

template<class F> void digraphE<F>::DeleteEdgesParallel( const vec<Bool>& to_delete )
{    ForceAssertEq(to_delete.isize(), EdgeObjectCount());
     #pragma omp parallel for
     for ( int v = 0; v < N( ); v++ )
     {    for ( int j = From(v).isize( ) - 1; j >= 0; j-- )
          {    int e = from_edge_obj_[v][j];
               if (to_delete[e]) {
                 from_[v].erase( from_[v].begin( ) + j );
                 from_edge_obj_[v].erase( from_edge_obj_[v].begin( ) + j ); } }
          for ( int j = To(v).isize( ) - 1; j >= 0; j-- )
          {    int e = to_edge_obj_[v][j];
               if (to_delete[e]) {
                 to_[v].erase( to_[v].begin( ) + j );
                 to_edge_obj_[v].erase( to_edge_obj_[v].begin( ) + j ); }  }  }  }

template<class F> void digraphE<F>::DeleteEdgesParallel( const vec<int>& to_delete )
{    vec<Bool> to_deletex( E( ), False );
     for ( int i = 0; i < to_delete.isize( ); i++ )
          to_deletex[ to_delete[i] ] = True;
     DeleteEdgesParallel(to_deletex);    }

template<class F> void digraphE<F>::DeleteEdges( const vec<int>& to_delete,
     const vec<int>& to_left )
{    vec<int> to_delete_local;
     if ( !to_delete.UniqueOrdered( ) )
     {    to_delete_local = to_delete;
          UniqueSort(to_delete_local);    }
     const vec<int>& tod 
          = ( to_delete_local.nonempty( ) ? to_delete_local: to_delete );
     vec<int> vs;
     for ( int i = 0; i < to_delete.isize( ); i++ )
          vs.push_back( to_left[ to_delete[i] ] );
     UniqueSort(vs);
     for ( int i = 0; i < vs.isize( ); i++ )
     {    int v = vs[i];
          for ( int j = From(v).isize( ) - 1; j >= 0; j-- )
          {    int e = EdgeObjectIndexByIndexFrom( v, j );
               if ( BinMember( tod, e ) ) DeleteEdgeFrom( v, j );    }    }    }

template<class F> void digraphE<F>::DeleteEdgesWithUpdate( 
     const vec<int>& to_delete, vec<int>& to_left, vec<int>& to_right )
{    DeleteEdges( to_delete, to_left );
     for ( int i = 0; i < to_delete.isize( ); i++ )
     {    to_left[ to_delete[i] ] = -1;
          to_right[ to_delete[i] ] = -1;    }    }

template<class F> int digraphE<F>::EdgeObjectIndexByIndexTo( int v, int j ) const
{    CheckGoodVertex(v);
     AssertGe( j, 0 );
     AssertLt( j, to_edge_obj_[v].isize( ) );
     return to_edge_obj_[v][j];    }
   
template<class F> int digraphE<F>::EdgeObjectIndexToFromIndex( int v, int e ) const
{    AssertGe( v, 0 );
     AssertLt( v, from_edge_obj_.isize( ) );
     for ( int i = 0; i < from_edge_obj_[v].isize( ); i++ )
          if ( from_edge_obj_[v][i] == e ) return i;
     return -1;    }

template<class F> int digraphE<F>::EdgeObjectIndexToToIndex( int v, int e ) const
{    AssertGe( v, 0 );
     AssertLt( v, to_edge_obj_.isize( ) );
     for ( int i = 0; i < to_edge_obj_[v].isize( ); i++ )
          if ( to_edge_obj_[v][i] == e ) return i;
     return -1;    }

template<class F> bool operator!=( const digraphE<F>& g1, const digraphE<F>& g2 )
{ return !(g1==g2); }

template<class F> bool EqualExceptEdgeObjectOrder( 
     const digraphE<F>& g1, const digraphE<F>& g2 )
{
    if ( static_cast<digraph const&>(g1) != static_cast<digraph const&>(g2) )
        return false;

    // digraphs are the same, now check edge objects
    typedef vec<int> V;
    typedef V::const_iterator VI;
    typedef vec<V> VV;
    typedef VV::const_iterator VVI;
    VV const& vv1 = g1.FromEdgeObj();
    VV const& vv2 = g2.FromEdgeObj();
    if ( vv1.size() != vv2.size() )
        return false;

    VVI oE(vv1.end());
    for ( VVI o1(vv1.begin()), o2(vv2.begin()); o1 != oE; ++o1, ++o2 )
    {
        if ( o1->size() != o2->size() )
            return false;

        VI iE(o1->end());
        for ( VI i1(o1->begin()), i2(o2->begin()); i1 != iE; ++i1, ++i2 )
            if ( !(g1.EdgeObject(*i1) == g2.EdgeObject(*i2)) )
                return false;
    }
    return true;
}

template<class E> bool operator==( const digraphE<E>& g1, const digraphE<E>& g2 )
{   if ( !EqualExceptEdgeObjectOrder( g1, g2 ) ) return false;
    return g1.Edges( ) == g2.Edges( );   }

template<class E> 
     void Compare( ostream& out, const digraphE<E>& g1, const digraphE<E>& g2 )
{    if ( g1.N( ) != g2.N( ) )
          cout << "first graph has " << g1.N( ) << " vertices but "
               << "second graph has " << g2.N( ) << "\n";
     if ( g1.From( ) != g2.From( ) ) cout << "from_ not the same\n";
     if ( g1.To( ) != g2.To( ) ) cout << "to_ not the same\n";
     if ( g1.Edges( ) != g2.Edges( ) ) cout << "edges_ not the same\n";
     if ( g1.ToEdgeObj( ) != g2.ToEdgeObj( ) )
          cout << "to_edge_obj_ not the same\n";
     if ( g1.FromEdgeObj( ) != g2.FromEdgeObj( ) )
          cout << "from_edge_obj_ not the same\n";
     if ( g1 != g2 ) cout << "DIGRAPHS ARE NOT EQUAL\n";
     return;    }

template<class E> void digraphE<E>::Clear( )
{    from_.clear( ), to_.clear( );
     from_edge_obj_.clear( ), to_edge_obj_.clear( );
     edges_.clear( );    }

template<class E> const E& digraphE<E>::EdgeObject( int i ) const
{    AssertGe( i, 0 );
     AssertLt( i, edges_.isize( ) );
     return edges_[i];    }

template<class E> E& digraphE<E>::EdgeObjectMutable( int i )
{    AssertGe( i, 0 );
     AssertLt( i, edges_.isize( ) );
     return edges_[i];    }

template<class V> const V& digraphV<V>::Vert( int v ) const
{    AssertGe( v, 0 );
     AssertLt( v, N( ) );
     return verts_[v];    }

template<class V> V& digraphV<V>::VertMutable( int v )
{    AssertGe( v, 0 );
     AssertLt( v, N( ) );
     return verts_[v];    }

template<class V, class E> const V& digraphVE<V,E>::Vert( int v ) const
{    AssertGe( v, 0 );
     AssertLt( v, N( ) );
     return verts_[v];    }

template<class V, class E> V& digraphVE<V,E>::VertMutable( int v )
{    AssertGe( v, 0 );
     AssertLt( v, N( ) );
     return verts_[v];    }

template<class V> void digraphV<V>::DeleteVertex( const int v )
{    int n = N( );
     AssertGe( v, 0 );
     AssertLt( v, n );
     DeleteEdgesAtVertex(v);
     verts_.erase( verts_.begin( ) + v );
     from_.erase( from_.begin( ) + v );
     to_.erase( to_.begin( ) + v );
     for ( int x = 0; x < n - 1; x++ )
     {    for ( int j = 0; j < From(x).isize( ); j++ )
               if ( From(x)[j] >= v ) FromMutable(x)[j]--;
          for ( int j = 0; j < To(x).isize( ); j++ )
               if ( To(x)[j] >= v ) ToMutable(x)[j]--;    }    }

template<class V> void digraphV<V>::DeleteVertices( const vec<int>& v )
{    for ( int m = v.isize( ) - 1; m >= 0; m-- )
          DeleteVertex( v[m] );    }

template<class V> int digraphV<V>::AddVertex( const V& v )
{    verts_.push_back(v);
     from_.resize( from_.size( ) + 1 );
     to_.resize( to_.size( ) + 1 );
     return verts_.size() - 1;    }

template<class V, class E> void digraphVE<V,E>::AddVertex( const V& v )
{    verts_.push_back(v);
     this->FromMutable( ).resize( this->From( ).size( ) + 1 );
     this->ToMutable( ).resize( this->To( ).size( ) + 1 );
     this->FromEdgeObjMutable( ).resize( this->FromEdgeObj( ).size( ) + 1 );
     this->ToEdgeObjMutable( ).resize( this->ToEdgeObj( ).size( ) + 1 );    }

template<class V, class E> void digraphVE<V,E>::RemoveVertices( const vec<int>& to_remove )
{
    vec<Bool> to_delete(verts_.size(),False );
    for( auto entry: to_remove){
        to_delete[entry]=True;
        digraphE<E>::DeleteEdgesAtVertex(entry);
    }
    digraphE<E>::RemoveEdgelessVertices(to_remove);
    EraseIf(verts_,to_delete);
}

template<class E> 
vec<int> digraphE<E>::EdgesSomewhereBetween( const int v, const int w ) const
{    vec<int> answer, after_v, before_w, both;
     GetSuccessors1( v, after_v ), GetPredecessors1( w, before_w );
     Intersection( after_v, before_w, both );
     for ( int l = 0; l < both.isize( ); l++ )
     {    int s = both[l];
          for ( int j = 0; j < From(s).isize( ); j++ )
          {    int t = From(s)[j];
               if ( BinMember( both, t ) ) 
                    answer.append( EdgesBetween( s, t ) );    }    }
     UniqueSort(answer);
     return answer;    }

template<class E>
void digraphE<E>::writeBinary( BinaryWriter& writer ) const
{
    digraph::writeBinary(writer);
    writer.write(from_edge_obj_);
    writer.write(to_edge_obj_);
    writer.write(edges_);   }

template<class E>
void digraphE<E>::readBinary( BinaryReader& reader )
{
    digraph::readBinary(reader);
    reader.read(&from_edge_obj_);
    reader.read(&to_edge_obj_);
    reader.read(&edges_);   }

template<class F> void digraphEX<F>::writeBinary( BinaryWriter& writer ) const
{   digraphX::writeBinary(writer);
    writer.write(from_edge_obj_);
    writer.write(to_edge_obj_);
    writer.write(edges_);
    writer.write(to_left_);
    writer.write(to_right_);     }

template<class F> void digraphEX<F>::readBinary( BinaryReader& reader )
{   digraphX::readBinary(reader);
    reader.read(&from_edge_obj_);
    reader.read(&to_edge_obj_);
    reader.read(&edges_);
    reader.read(&to_left_);
    reader.read(&to_right_);    }

template<class V>
void digraphV<V>::writeBinary( BinaryWriter& writer ) const
{
    digraph::writeBinary(writer);
    writer.write(verts_);  }

template<class V>
void digraphV<V>::readBinary( BinaryReader& reader )
{   
    digraph::readBinary(reader);
    reader.read( &verts_ );    }

template<class V, class E>
void digraphVE<V,E>::writeBinary( BinaryWriter& writer ) const
{    digraphE<E>::writeBinary(writer);
     writer.write(verts_);  }

template<class V, class E>
void digraphVE<V,E>::readBinary( BinaryReader& reader )
{    digraphE<E>::readBinary(reader);
     reader.read( &verts_ );    }

template<class E> void EmbeddedSubPath<E>::TestValid( ) const
{    ForceAssertEq( e_.isize( ), a_.isize( ) - 1 );
     for ( int u = 0; u < a_.isize( ) - 1; u++ )
     {    const vec<int>& fr = D_->From( a_[u] );
          ForceAssertGe( e_[u], 0 );
          ForceAssertLt( e_[u], fr.isize( ) );
          ForceAssertEq( fr[ e_[u] ], a_[u+1] );
          ForceAssertEq( D_->EdgeObjectIndexByIndexFrom( a_[u], e_[u] ),
               esafe_[u] );    }    }

template<class E> void EmbeddedSubPath<E>::Repair( )
{    for ( int u = 0; u < e_.isize( ); u++ )
     {    if ( D_->EdgeObjectIndexByIndexFrom( a_[u], e_[u] ) != esafe_[u] )
               e_[u] = D_->EdgeObjectIndexToFromIndex(
                    a_[u], esafe_[u] );    }    }

template<class E> void DistancesToEnd3( const digraphE<E>& G,
     int (E::*len)( ) const, const int max_dist, const Bool fw, vec<int>& D,
     vec<Bool>& complete, const int max_paths )
{
     D.resize( G.N( ), 0 );
     complete.resize( G.N( ) );
     #pragma omp parallel for
     for ( int v = 0; v < D.isize( ); v++ )
     {    vec< pair< vec<int>, int > > 
               paths( { make_pair( vec<int>({v}), 0 ) } );
          while( paths.isize( ) <= max_paths )
          {    vec< pair< vec<int>, int > > paths2;
               for ( const auto& p : paths )
               {    int x = p.first.back( );
                    vec< pair<int,int> > ext;
                    for ( int j = 0; 
                         j < ( fw ? G.From(x).isize( ) : G.To(x).isize( ) ); j++ )
                    {    int y = ( fw ? G.From(x)[j] : G.To(x)[j] );
                         if ( Member( p.first, y ) ) continue;
                         int e = ( fw ? G.EdgeObjectIndexByIndexFrom( x, j )
                              : G.EdgeObjectIndexByIndexTo( x, j ) );
                         int l = (G.EdgeObject(e).*len)( );
                         ext.push( y, l );    }
                    ReverseSort(ext);
                    for ( int i = 0; i < ext.isize( ); i++ )
                    {    int j;
                         for ( j = i + 1; j < ext.isize( ); j++ )
                              if ( ext[j].first != ext[i].first ) break;
                         auto q(p);
                         q.first.push_back( ext[i].first );
                         q.second += ext[i].second;
                         paths2.push_back(q);
                         i = j - 1;   }
                    if ( ext.empty( ) ) paths2.push_back(p);    }
               if ( paths2 == paths ) break;
               paths = paths2;    }
          complete[v] = ( paths.isize( ) <= max_paths );
          for ( int i = 0; i < paths.isize( ); i++ )
               D[v] = Max( D[v], paths[i].second );    }    }

template<class E> void RemoveHangingEnds3( digraphE<E>& G, 
     int (E::*len)( ) const, const int max_del, const double min_ratio,
     const int max_paths )
{
     // Track hanging ends.

     vec<Bool> hanging( G.EdgeObjectCount( ), False );

     // Define the maximum length that we care about.

     const int max_dist = int( ceil( double(max_del) * min_ratio ) );

     // Go through two passes (forward and reverse).

     for ( int pass = 1; pass <= 2; pass++ )
     {
          // Compute distances to end.

          vec<int> D;
          vec<Bool> complete;
          DistancesToEnd3( G, len, max_dist, pass == 1, D, complete, max_paths );

          // Identify hanging ends.

          #pragma omp parallel for
          for ( int v = 0; v < G.N( ); v++ )
          {    const vec<int>& V = ( pass == 1 ? G.From(v) : G.To(v) );
               vec<int> d( V.size( ) ), id( V.size( ), vec<int>::IDENTITY );
               vec<Bool> c( V.size( ) );
               for ( int j = 0; j < V.isize( ); j++ )
               {    d[j] = ((pass == 1 
                         ? G.EdgeObjectByIndexFrom(v,j) : G.EdgeObjectByIndexTo(v,j))
                         .*len)( ) + D[ V[j] ];
                    c[j] = complete[ V[j] ];    }
               ReverseSortSync( d, c, id );
               for ( int j = 1; j < d.isize( ); j++ )
               {    if ( d[j] <= max_del && d[0] >= d[j] * min_ratio && c[j] )
                    {    hanging[ ( pass == 1 
                              ? G.EdgeObjectIndexByIndexFrom( v, id[j] )
                              : G.EdgeObjectIndexByIndexTo( v, id[j] ) ) ] 
                              = True;    }    }    }    }

     // Remove hanging ends.

     vec<int> to_delete;
     for ( int i = 0; i < G.EdgeObjectCount( ); i++ )
          if ( hanging[i] ) to_delete.push_back(i);
     G.DeleteEdges(to_delete);    }

template<class E> int64_t digraphE<E>::CheckSum( ) const
{    int64_t x = 0;
     for ( int v = 0; v < N( ); v++ )
     {    for ( int j = 0; j < From(v).isize( ); j++ )
               x += ( v + 1 ) * (j + 1 ) * ( From(v)[j] + 1 );
          for ( int j = 0; j < To(v).isize( ); j++ )
               x += ( v + 1 ) * (j + 1 ) * ( To(v)[j] + 1 );
          for ( int j = 0; j < FromEdgeObj(v).isize( ); j++ )
               x += ( v + 1 ) * (j + 1 ) * ( FromEdgeObj(v)[j] + 1 );
          for ( int j = 0; j < ToEdgeObj(v).isize( ); j++ )
               x += ( v + 1 ) * (j + 1 ) * ( ToEdgeObj(v)[j] + 1 );    }
     return x;    }

template<class E> int64_t digraphEX<E>::CheckSum( ) const
{    int64_t x = 0;
     for ( int v = 0; v < N( ); v++ )
     {    for ( int j = 0; j < From(v).size( ); j++ )
               x += ( v + 1 ) * (j + 1 ) * ( From(v)[j] + 1 );
          for ( int j = 0; j < To(v).size( ); j++ )
               x += ( v + 1 ) * (j + 1 ) * ( To(v)[j] + 1 );
          for ( int j = 0; j < From(v).size( ); j++ )
               x += ( v + 1 ) * (j + 1 ) * ( IFrom(v,j) + 1 );
          for ( int j = 0; j < To(v).size( ); j++ )
               x += ( v + 1 ) * (j + 1 ) * ( ITo(v,j) + 1 );    }
     return x;    }

template<class E> const E& digraphEX<E>::EdgeObject( int i ) const
{    AssertGe( i, 0 );
     AssertLt( i, (int) edges_.size( ) );
     return edges_[i];    }

template<class F> digraphEX<F>::digraphEX( const digraphE<F>& G )
{    from_.resize( G.N( ) );
     to_.resize( G.N( ) );
     from_edge_obj_.resize( G.N( ) );
     to_edge_obj_.resize( G.N( ) );
     for ( int i = 0; i < G.N( ); i++ )
     {    from_[i].resize( G.From(i).size( ) );
          for ( int j = 0; j < G.From(i).isize( ); j++ )
               from_[i][j] = G.From(i)[j];
          to_[i].resize( G.To(i).size( ) );
          for ( int j = 0; j < G.To(i).isize( ); j++ )
               to_[i][j] = G.To(i)[j];
          from_edge_obj_[i].resize( G.From(i).size( ) );
          for ( int j = 0; j < G.From(i).isize( ); j++ )
               from_edge_obj_[i][j] = G.IFrom( i, j );
          to_edge_obj_[i].resize( G.To(i).size( ) );
          for ( int j = 0; j < G.To(i).isize( ); j++ )
               to_edge_obj_[i][j] = G.ITo( i, j );    }
     edges_.resize( G.EdgeObjectCount( ) );
     for ( int e = 0; e < G.EdgeObjectCount( ); e++ )
          edges_[e] = G.EdgeObject(e);
     to_left_.resize( E( ) ), to_right_.resize( E( ) );
     for ( int v = 0; v < N( ); v++ )
     {    for ( int j = 0; j < (int) From(v).size( ); j++ )
               to_left_[ IFrom( v, j ) ] = v;
          for ( int j = 0; j < (int) To(v).size( ); j++ )
               to_right_[ ITo( v, j ) ] = v;    }    }

template<class F> digraphE<F> digraphEX<F>::AsDigraphE() const {
    vec< vec<int> > from(N()), to(N());
    vec< vec<int> > to_edge_obj(N()), from_edge_obj(N());
    vec<F> edges(E());
    for ( int i = 0; i < N( ); i++ ) {
	from[i].resize( From(i).size( ) );
	for ( size_t j = 0; j < From(i).size( ); j++ )
	    from[i][j] = From(i)[j];
	to[i].resize( To(i).size( ) );
	for ( size_t j = 0; j < To(i).size( ); j++ )
	    to[i][j] = To(i)[j];
	from_edge_obj[i].resize( From(i).size( ) );
	for ( size_t j = 0; j < From(i).size( ); j++ )
	    from_edge_obj[i][j] = IFrom( i, j );
	to_edge_obj[i].resize( To(i).size( ) );
	for ( size_t j = 0; j < To(i).size( ); j++ )
	    to_edge_obj[i][j] = ITo( i, j );    
    }
    for ( int e = 0; e < E(); e++ )
	edges[e] = EdgeObject(e);
    
    return digraphE<F>(from, to, edges, to_edge_obj, from_edge_obj, true);
}

template<class V> void digraphV<V>::DeleteEdgeTo( int w, int j )
{    int v = to_[w][j];
     int i = InputToOutputFrom( w, j );
     to_[w].erase( to_[w].begin( ) + j );
     from_[v].erase( from_[v].begin( ) + i );    }

template<class V> void digraphV<V>::DeleteEdgeFrom( int v, int j )
{    int w = from_[v][j];
     int i = InputFromOutputTo( v, j );
     from_[v].erase( from_[v].begin( ) + j );
     to_[w].erase( to_[w].begin( ) + i );    }

#endif
