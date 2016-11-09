///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// =================================================================================
// PLEASE READ: THERE ARE OTHER FILES YOU MAY WANT TO LOOK AT:
// - DigraphPaths.h
// - FindCells.h
// Please add to this list if you add new files having digraph functions.
// =================================================================================

// This file contains the definition of class digraph, which represents a finite
// directed graph, digraphV<V>, which has vertex objects from class V, and
// digraphE<E>, which has edge objects from class E.  There is also a class 
// digraphVE<V,E>.
//
// The file Digraph.cc contains some of the definitions for digraph functions.
//
// The file DigraphTemplate.h contains some of the templatized definitions for
// digraphE functions.  Putting these in a separate file makes this file more
// readable, but has the additional advantage of reducing the inlining of these
// functions, thereby reducing compilation time and executable size (in principle).
//
// DANGEROUS THINGS
// - The fewer places DigraphTemplate.h is included, the better.  In particular,
//   as with all template instantiations, if somehow a template gets instantiated
//   twice (perhaps because DigraphTemplate.h was included in two locations, leading
//   to duplicate instantiations), the linker may handle this badly, choosing
//   randomly among the multiple instantiations and potentially injecting
//   unanticipated link-time dependencies on unrelated code.  It's a bear to
//   sort this out.  So the rule is to explicitly instantiate needed methods for
//   a particular template parameter in one, and only one place.
// - All the digraph classes are based on representation of vertex and edge indices
//   as ints.  This would be problematic if a graph were to contain > ~2 x 10^9
//   vertices or edges.

#ifndef DIGRAPH_H
#define DIGRAPH_H

#include "Bitvector.h"
#include "CoreTools.h"
#include "Equiv.h"
#include "Intvector.h"
#include "VecUtilities.h"
#include "feudal/BinaryStream.h"
#include "math/Functions.h"
#include "math/Matrix.h"
#include "system/TraceVal.h"
#include <cstddef>

typedef pair<int, int> VertexPair;

typedef int vrtx_t;
typedef int edge_t;

/**
   Class: digraph

   A digraph may have multiple edges between two given vertices, but no edges
   from a vertex to itself.

   A digraph is stored as two vec< vec<int> >'s, "to" and "from".  For each 
   vertex v, to[v] is the *sorted* list of vertices which have an edge to v, and 
   from[v] is the *sorted* list of vertices which have an edge from v.
*/

// =================================================================================
// ============================ DIGRAPH CLASS ======================================
// =================================================================================

class digraph
{
     public:

     digraph( ) { }
     digraph( const vec< vec<int> >& from, const vec< vec<int> >& to );
     digraph( const matrix<Bool>& );
     explicit digraph( int n ) { Initialize( n ); }

     // "Reducing" constructors
     // Reduce g to the vertices in v

     digraph( const digraph& g, const vec<int>& v );

     // Reduce g to the vertex/vertices in component #i

     digraph( const digraph& g, const int i );
     
     void Initialize( const vec< vec<int> >& from, const vec< vec<int> >& to );
     void Initialize( const digraph& g, const vec<int>& v );
     void Initialize( const matrix< Bool >& );

     // Create an empty graph with n nodes.

     void Initialize( int64_t n ) {
       Clear();
       from_.resize( n );
       to_.resize( n );
     }

     Bool TestValid( const Bool exit = True ) const; // incomplete test
     Bool TestValidParallel( ) const; // incomplete test

     void Clear( )
     {    from_.clear( ), to_.clear( );    }

     int64_t N( ) const { return from_.size( ); } // number of vertices

     void CheckGoodVertex( int64_t v ) const
     {    AssertGe( v, 0 );
          AssertLt( v, N( ) );    }

     const vec< vec<int> >& From( ) const { return from_; }
     const vec< vec<int> >& To( ) const { return to_; }
     vec< vec<int> >& FromMutable( ) { return from_; }
     vec< vec<int> >& ToMutable( ) { return to_; }

     const vec<int>& From( int64_t v ) const { CheckGoodVertex(v); return from_[v]; }
     const vec<int>& To( int64_t v ) const { CheckGoodVertex(v); return to_[v]; }
     int FromSize( int64_t v ) const { CheckGoodVertex(v); return from_[v].size( ); }
     int ToSize( int64_t v ) const { CheckGoodVertex(v); return to_[v].size( ); }

     vec<int>& FromMutable( int64_t v ) { CheckGoodVertex(v); return from_[v]; }
     vec<int>& ToMutable( int64_t v ) { CheckGoodVertex(v); return to_[v]; }

     // Create new vertices.

     void AddVertices( int64_t nadd ) {
       int64_t nvert = N( );
       from_.resize( nvert + nadd );
       to_  .resize( nvert + nadd );
     }

     // Add an edge from v to w.

     void AddEdge( int64_t v, int64_t w );
  
     Bool HasEdge( int64_t v, int64_t w ) const;

     Bool Source( int64_t v ) const { return To(v).empty( ); }
     Bool Sink( int64_t v ) const { return From(v).empty( ); }

     // ============================================================================
     // ==================== METHODS TO FIND A SUBGRAPH ============================
     // ============================================================================

     // All these methods are const and return a set of vertices or edges in the
     // graph.

     // Find the sources or sinks in a given graph.

     void Sources( vec<int>& v ) const;
     void Sinks( vec<int>& v ) const;

     // Given the complete subgraph defined by a set of vertices S, find its sources
     // and sinks.  The set S must be sorted.

     void SubgraphSources( const vec<int>& S, vec<int>& v ) const;
     void SubgraphSinks( const vec<int>& S, vec<int>& v ) const;

     // GetPredecessors: find all vertices which have a directed path to a vertex
     // in v.  This includes the vertices in v by definition.  Return a sorted list 
     // to_v.  GetSuccessors: go the other way.  
     // See also another version of GetSuccessors for class digraphE.

     void GetPredecessors( const vec<int>& v, vec<int>& to_v ) const;
     void GetSuccessors( const vec<int>& v, vec<int>& from_v ) const;
     void GetPredecessors1( const int v, vec<int>& to_v ) const;
     void GetSuccessors1( const int v, vec<int>& from_v ) const;

     // Return the connected component containing a given vertex.  Returns a sorted
     // list of vertices.

     vec<int> ComponentOf( const int v ) const;

     // VerticesConnectedTo.  Return all vertices that are connected directly or
     // indirectly to a vertex in a given list (which need not be ordered and may 
     // contain duplicates).  Return a sorted list of vertices.  If provided a 
     // single vertex as input, this is the same as ComponentOf.  

     vec<int> VerticesConnectedTo( const vec<int>& v ) const;

     // ============================================================================
     // ======================== END ABOVE SECTION =================================
     // ============================================================================

     // Reverse the graph.

     void Reverse( );

     // Components: find the connected components.  Each component is a sorted list
     // of vertices.  If invisible vertices are provided, they are ignored.
     // ComponentsAlt orders the components differently.
     
     void Components( vec< vec<int> >& comp, const vec<Bool>* invisible = NULL ) 
          const;
     void ComponentsAlt( vec< vec<int> >& comp ) const;
     size_t NComponents() const;

     // Determine if a graph has a path of nonzero length from v to v.

     Bool LoopAt( const int v ) const;

     // Return the strongly connected components of a graph.  The code is an
     // algorithm of Tarjan.  We borrowed pseudocode from Wikipedia.  The answer
     // is a sorted vector of sorted vectors.

     void StronglyConnectedComponents( vec< vec<int> >& SCC ) const;

     // Determine if the connected component defined by vertices sub has a directed
     // cycle.  Don't call this with anything other than a connected component.

     Bool HasCycle( const vec<int>& sub ) const;

     // Determine if graph is acyclic.  Can be very slow.  See CyclicCore, below.

     Bool Acyclic( ) const;

     // CyclicCore: return a subset of vertices that define a subgraph having no 
     // sources and sinks, and which is empty iff the graph is acyclic.  This 
     // should be faster than Acyclic.  The cyclic core should be the union of all
     // vertices that appear in cycles.

     void CyclicCore( vec<int>& core ) const;

     // ComponentRelation: return equivalence relation corresponding to connected 
     // components.  Note that this can be very slow on a large graph.  The function
     // "Components" can return the same information, faster.

     void ComponentRelation( equiv_rel& e ) const;

     // Return number of connected components in a graph.

     int ConnectedComponents( ) const;

     // Return a sorted list of all vertices whose removal would increase the 
     // number of components in the graph.  Complexity = O( sum( n^2 log(n) ) ),
     // where n ranges over the component sizes.

     void CutPoints( vec<int>& cuts ) const;

     // Given an edge referred to by To, find the corresponding edge in From,
     // and the other way.

     int InputToOutputFrom( int w, int i ) const;
     int InputFromOutputTo( int w, int i ) const;

     // Delete an edge.  Indices of edges <j are untouched, while those >j decrease 
     // by 1.  In particular, adding an edge and then deleting it is guaranteed to 
     // leave the indices of other edges unchanged.

     void DeleteEdgeTo( int w, int j );
     void DeleteEdgeFrom( int v, int j );

     // Remove all edges entering or exiting a given vertex.

     void DeleteEdgesAtVertex( int v );

     // AllPaths: find all paths between two vertices which contain no duplicated
     // vertices.  Zero length paths (i.e. paths from a vertex to the same vertex)
     // are included.  Return answer as a list of lists, with each inner list a 
     // sequence of vertices.  This last can then be expanded out to reflect
     // alternative edges between adjacent vertices.
     //
     // If w = -1, instead find all paths from v to a sink.
     // If v = -1, instead find all paths from a source to w.
     //
     // If maxpaths >= 0, fail, returning False if greater than that many paths or 
     // partial paths found.
     //
     // If allow_self_loop and v = w (and >= 0), then loops from v to v are
     // allowed (but v cannot appear more than twice).
     //
     // Note that the behavior of this code may not be as expected.  Thus in this 
     // case
     //                    y        x
     //      * ------> * <----- * -----> *
     //                  -----> 
     //
     // AllPaths( -1, -1, ... ) will not use edge y.

     Bool AllPaths( int v, int w, vec< vec<int> >& paths, int maxpaths = -1,
          const Bool allow_self_loop = False, const int maxpushes = -1 ) const;

     // TransferEdges( v, w ): move all edges entering and exiting v, so that they
     // enter and exit w instead, leaving v without any edges entering and exiting 
     // it.  If enter_only, only transfer entering edges.
     // Also for digraphE.

     void TransferEdges( int v, int w, const Bool enter_only = False );

     // Find simple lines in a graph: these are the unbranched paths in it,
     // returned here as sequences of vertices.  The starting point of a circle is
     // somewhat accidental.  For circles the first vertex equals the last vertex.

     void FindSimpleLines( vec<vec<int>>& lines ) const;

     // Similar but slightly different definition:

     void FindSimpleLines2( vec<vec<int>>& lines ) const;

     // ============================================================================
     // ======================= PRINTERS AND WRITERS ===============================
     // ============================================================================

     // Create a representation of a given graph in DOT:
     // http://www.graphviz.org/doc/info/lang.html
     // The edge_colors and edge_labels arguments should be in bijective 
     // correspondence with from_.

     void DOT( ostream& out ) const;
     void DOT( ostream& out, const vec<String>& vertex_colors ) const;
     void DOT( ostream& out, const vec<String>& vertex_colors,
          const vec< vec<String> >& edge_colors ) const;
     void DOT( ostream& out, const vec< vec<String> >& edge_labels ) const;
     void DOT( ostream& out, const vec< vec<String> >& edge_labels,
          const vec<String>& vertex_colors ) const;
     void DOT( ostream& out, const vec< vec<String> >& edge_labels,
          const vec<String>& vertex_colors, 
          const vec<vec<String>>& edge_colors, const vec<vec<String>>& edge_attrs,
          const vec<vec<String>>& legends, const vec<String>& legend_colors ) const;
     void PrettyDOT( ostream& out, Bool label_contigs = True, 
		     Bool label_vertices = False,
		     const vec<int>* componentsToPrint = NULL, 
		     const vec<String> *label_contigs_extra = NULL,
		     const vec<int> *verticesToPrint = NULL ) const;

     // DOT_vl: generate a graph with vertex labels in circles.  There is an 
     // optional argument layout, which could be set to circo if you want 
     // circular layout.  There is another optional argument legend, which is 
     // rendered inside a square box, but at a random place (not necessarily at
     // the top).  Lines within the box are left-right centered.  The version with
     // legends is similar: you can have more than one.
     // vshape: vertex shape, default = ellipse

     void DOT_vl( ostream& out, const vec<String> & vertex_labels,
          const String& layout = "", 
          const vec<String>& legend = vec<String>( ), 
          const String& color = "", const String& vshape = "" ) const;

     void DOT_vl( ostream& out, const vec<String> & vertex_labels,
          const String& layout = "", 
          const vec< vec<String> >& legends = vec< vec<String> >( ),
          const vec<String>& colors = vec<String>( ),
          const vec< vec<String> >& edge_labels = vec< vec<String> >( ),
          const vec<String>& vertex_colors = vec<String>( ),
          const String& vshape = "" ) const;

     void DOT_el( ostream& out, const String& layout = "", 
          const vec< vec<String> >& edge_labels = vec< vec<String> >( ) ) const;

     // Print the subgraph defined by the given vertices vs, which need
     // not be sorted.

     void DOT_vl( ostream& out, const vec<String> & vertex_labels, 
          vec<int> vs, const vec<String>& vertex_colors = vec<String>( ) ) const;

     // Create a representation of a given graph in GraphML:
     // http://graphml.graphdrawing.org
     // The edge_labels arguments should be in bijective 
     // correspondence with from_.

     void WriteGraphML( ostream& out, const vec< vec<String> >& edge_labels ) const;

     void writeBinary( BinaryWriter& writer ) const
     { writer.write(from_); }

     void readBinary( BinaryReader& reader );

     // ============================================================================
     // ======================== END ABOVE SECTION =================================
     // ============================================================================

     /// ReplaceWithTransitiveClosure()
     /// Replace this graph with its transitive closure.  That is, if there is a
     /// directed path from u to w in the original graph, there is a direct edge
     /// (u,w) in the new graph.

     void ReplaceWithTransitiveClosure();

     static size_t externalSizeof() { return 0; }

     friend bool operator==( digraph const& g1, digraph const& g2 )
     { return g1.from_ == g2.from_; }

     friend bool operator!=( digraph const& g1, digraph const& g2 )
     { return g1.from_ != g2.from_; }

     protected:

     vec< vec<int> > from_, to_;

};

SELF_SERIALIZABLE(digraph);

// =================================================================================
// ============================ DIGRAPHX CLASS =====================================
// =================================================================================

class digraphX
{
     public:
     
     digraphX( ) { }

     int N( ) const { return from_.size( ); } // number of vertices

     void CheckGoodVertex( int v ) const
     {    AssertGe( v, 0 );
          AssertLt( v, N( ) );    }

     const SerfVec<int>& From( int v ) const { CheckGoodVertex(v); return from_[v]; }
     const SerfVec<int>& To( int v ) const { CheckGoodVertex(v); return to_[v]; }

     void Components( vec< vec<int> >& comp, const vec<Bool>* invisible = NULL ) 
          const;

     void writeBinary( BinaryWriter& writer ) const
     { writer.write(from_);
       writer.write(to_); }

     void readBinary( BinaryReader& reader )
     { reader.read(&from_);
       reader.read(&to_); }

     protected:

     VecIntVec from_, to_;

};

SELF_SERIALIZABLE(digraphX);

// =================================================================================
// ============================ DIGRAPHV CLASS =====================================
// =================================================================================

template<class V> class digraphV : public digraph {

     public:

     Bool TestValid( const Bool exit = True ) const; // incomplete test

     digraphV( ) { }

     digraphV( const vec< vec<int> >& from, const vec< vec<int> >& to,
          const vec<V>& verts );
     void Initialize( const vec< vec<int> >& from, const vec< vec<int> >& to,
          const vec<V>& verts );

     const V& Vert( int v ) const;
     V& VertMutable( int v );

     int AddVertex( const V& v );

     // DeleteVertex: delete a vertex and all edges incident upon it.
     // Deleting a vertex renumbers everything.

     void DeleteVertex( const int v );

     // DeleteVertices take a unique-sorted vector as input.

     void DeleteVertices( const vec<int>& v );

     void RemoveEdgelessVertices( );

     // Delete an edge.  Indices of edges <j are untouched, while those >j decrease 
     // by 1.  In particular, adding an edge and then deleting it is guaranteed to 
     // leave the indices of other edges unchanged.

     void DeleteEdgeTo( int w, int j );
     void DeleteEdgeFrom( int v, int j );

     void writeBinary( BinaryWriter& writer ) const;
     void readBinary( BinaryReader& reader );

     private:

     vec<V> verts_;

};
extern template class digraphV< vec<int> >;

template <class V>
struct Serializability< digraphV<V> >
{ typedef SelfSerializable type; };


// =================================================================================

template<class F> class digraphEX;       // forward declaration
template<class E> class EmbeddedSubPath; // forward declaration

// =================================================================================
// ============================ DIGRAPHE CLASS =====================================
// =================================================================================

template<class E> class digraphE;

/**
   Class Template: digraphE

   A digraphE<E> is stored as a digraph, plus a vec<E> edges_, plus indexing vectors
   (vec< vec<int> >'s to_edge_obj_ and from_edge_obj_) such that 
   if to_[w][i] = v, then to_edge_obj_[w][i] is the index in edges_ of the edge
   object for v --> w, and similarly for from_edge_obj_.

   The vec<E> edges_ is allowed to have unused entries (although not initially) and 
   there is a RemoveDeadEdgeObjects( ) function to eliminate them.

   Note that there are a bunch of member functions given here for digraphE<E>,
   that could also be implemented for class digraph.
*/

template<class F> class digraphE : public digraph {

     public:

     Bool TestValid( const Bool exit = True ) const; // incomplete test
     Bool TestValidParallel( ) const; // incomplete test

     enum ConstructorBehavior { EDGES_SEPARATE, EDGES_IN_LINE };
     enum ConstructorType1 { COMPLETE_SUBGRAPH };
     enum ConstructorType2 { COMPLETE_SUBGRAPH_EDGES };

     // Constructor 1: build the empty digraph.

     digraphE( ) { }

     // Constructor 2: build digraph from arbitrary digraph data.
     // (Note the INSANE order of the arguments: first from before to, then
     // to_edge_obj before from_edge_obj.)

     digraphE( const vec< vec<int> >& from, const vec< vec<int> >& to,
          const vec<F>& edges, const vec< vec<int> >& to_edge_obj,
          const vec< vec<int> >& from_edge_obj, 
          const Bool allow_unused_edges = False );
     void Initialize( const vec< vec<int> >& from, const vec< vec<int> >& to,
          const vec<F>& edges, const vec< vec<int> >& to_edge_obj,
          const vec< vec<int> >& from_edge_obj,
          const Bool allow_unused_edges = False );

     // Constructor 3a: given a collection of edge objects, create a graph having
     // one edge and two vertices for each of the edge objects.
     // [ invoked by: digraphE( const vec<F>& edges, EDGES_SEPARATE ); ]
     // Constructor 3b: given a collection of edge objects, create a line graph out
     // of the edges.
     // [ invoked by: digraphE( const vec<F>& edges, EDGES_IN_LINE ); ]

     enum ConstructorName {
       FROM_EDGE_EQUIV, /* 8b */
       FROM_SUBS /* 11 */
     };

     digraphE( const vec<F>& edges, const ConstructorBehavior constructor_type  );

     // Constructor 4: given a collection of edge objects, and an equivalence 
     // relation on them, build a graph having two vertices per equivalence class,
     // with one edge between those two vertices for each member of the equivalence
     // class.

     digraphE( const vec<F>& edges, const equiv_rel& e );

     // Same as above, but not a constructor:

     void EdgeEquivConstructor( const vec<F>& edges, const equiv_rel& e );

     // Constructor 5a: given a digraph, and given a list of vertex indices,
     // create the digraph having those vertices (with indices starting at 0,
     // but in the given order), and having all the edges that were between those
     // vertices.  Thus this is a "complete subgraph" constructor.  The order of the
     // edges is first by vertex v, and then by the order within From(v).

     digraphE( const ConstructorType1 constructor_type, const digraphE& g, 
          const vec<int>& v ); // call with first arg = COMPLETE_SUBGRAPH
     void Initialize( const ConstructorType1 constructor_type, const digraphE& g, 
          const vec<int>& v ); // call with first arg = COMPLETE_SUBGRAPH

     // Constructor 5b: given a digraph, and list of edge indices (not necessarily
     // sorted or unique),
     // create the digraph having just those edges (in the given order), 
     // and whose vertices consist of those vertices appearing as an end of one of 
     // the edges, in the same order and under the same equivalence relation implied 
     // by the original graph.  This is a sort of complete subgraph constructor.

     digraphE( const ConstructorType2 constructor_type, const digraphE& g, 
          const vec<int>& ed, const vec<int>& to_left, const vec<int>& to_right ); 
          // call with first arg = COMPLETE_SUBGRAPH_EDGES

     void Initialize( const ConstructorType2 constructor_type, 
          const digraphE& g, const vec<int>& ed, const vec<int>& to_left, 
          const vec<int>& to_right ); 
          // call with first arg = COMPLETE_SUBGRAPH_EDGES

     // Same, from digraphEX.

     digraphE( const ConstructorType2 constructor_type, const digraphEX<F>& g, 
          const vec<int>& ed ); 
          // call with first arg = COMPLETE_SUBGRAPH_EDGES

     void Initialize( const ConstructorType2 constructor_type, 
          const digraphEX<F>& g, const vec<int>& ed ); 
          // call with first arg = COMPLETE_SUBGRAPH_EDGES

     // Constructor 6: extract the nth connected component from another digraph.

     digraphE( const digraphE& g, int n );

     // Constructor 7: from another digraph and replacement edge objects.

     template<class G> digraphE( const digraphE<G>& g, const vec<F>& edges )
     {    from_ = g.From( );
          to_ = g.To( );
          from_edge_obj_ = g.FromEdgeObj( );
          to_edge_obj_ = g.ToEdgeObj( );
          edges_ = edges;    }

     template<class G> void Initialize ( const digraphE<G>& g, const vec<F>& edges )
     {    from_ = g.From( );
          to_ = g.To( );
          from_edge_obj_ = g.FromEdgeObj( );
          to_edge_obj_ = g.ToEdgeObj( );
          edges_ = edges;    }

     // Constructor 8: from a given digraph and an equivalence relation on the
     // vertices.  The vertices of the constructed digraph are the equivalence
     // classes of vertices of the given graph.
     
     digraphE( const digraphE& g, const equiv_rel& e );
     void Initialize( const digraphE& g, const equiv_rel& e );

     // Constructor 8b: from a given digraph and an equivalence relation on its
     // edges.  The edges of the constructed digraph are the equivalence classes
     // of edges in the given graph.  We assume that the values of the edges in an
     // equivalence class are the same.  A vertex equivalence relation follows from
     // the edge equivalence relation.

     digraphE( const ConstructorName cname /* FROM_EDGE_EQUIV */, const digraphE& g, 
          const equiv_rel& e );

     // Constructor 9: form the disjoint union of a collection of digraphs.

     explicit digraphE( const vec<digraphE>& g );
     void Initialize( const vec<digraphE>& g );

     // Constructor 10: from a collection of digraphs and a set of identifications
     // between vertices in their disjoint union, each of which is specified as
     // ( (g1,v1), (g2,v2) ) where g1, g2 refer to graphs and v1, v2 refer to
     // vertices on those graphs.

     digraphE( const vec<digraphE>& g, 
          const vec< pair< pair<int,int>, pair<int,int> > >& joins );
     void Initialize( const vec<digraphE>& g, 
          const vec< pair< pair<int,int>, pair<int,int> > >& joins );

     // Constructor 11: from a given digraph and a collection of subsets of its
     // edges, each of which is given the induced subgraph structure, which are then
     // merged into a disjoint union.  The numbering of edge objects in the new
     // graph is the obvious order C[0][0], C[0][1], ..., C[1][0], C[1][1], ... .

     digraphE( const ConstructorName cname /* FROM_SUBS */, const digraphE& g, 
          const vec< vec<int> >& C );

     // Constructor 12: from a digraph and the collection of edges corresponding
     // to digraph vertices

     digraphE( const digraph& g, const vec<F>& edges );
     void Initialize( const digraph& g, const vec<F>& edges ); 

     // Constructor 13: from a digraph (edges will correspond to vertex numbers 
     // in digraph) 

     explicit digraphE( const digraph& g );

     // Initialize 
     void Initialize( const int n );

     // Delete the entire graph.

     void Clear( );

     // Return a checksum.  This does NOT take into account the edge objects
     // themselves.

     int64_t CheckSum( ) const;

     // Given a list of vertex indices,
     // create the digraph having those vertices (with indices starting at 0,
     // but in the given order), and having all the edges that were between those
     // vertices.  Thus this is a "complete subgraph" constructor.  The order of the
     // edges is first by vertex v, and then by the order within From(v).

     digraphE Subgraph( const vec<int>& v ) const;

     // Return number of edge objects.

     int EdgeObjectCount( ) const { return edges_.size( ); }
     int E( ) const { return edges_.size( ); }

     // A bunch of ways to access edge objects:

     const F& EdgeObject( int i ) const;
     const F& O( int i ) const { return EdgeObject(i); }
     F& OMutable( int i ) { return EdgeObjectMutable(i); }

     F& EdgeObjectMutable( int i );

     const F& EdgeObjectByIndexFrom( int v, int j ) const
     {    CheckGoodVertex(v);
          AssertGe( j, 0 );
          AssertLt( j, from_edge_obj_[v].isize( ) );
          return edges_[ from_edge_obj_[v][j] ];    }

     F& EdgeObjectByIndexFromMutable( int v, int j )
     {    CheckGoodVertex(v);
          AssertGe( j, 0 );
          AssertLt( j, from_edge_obj_[v].isize( ) );
          return edges_[ from_edge_obj_[v][j] ];    }

     int EdgeObjectIndexByIndexFrom( int v, int j ) const
     {    CheckGoodVertex(v);
          AssertGe( j, 0 );
          AssertLt( j, from_edge_obj_[v].isize( ) );
          return from_edge_obj_[v][j];    }

     int IFrom( int v, int j ) const { return EdgeObjectIndexByIndexFrom( v, j ); }

     const F& OFrom( int v, int j ) const 
     {    return EdgeObject( from_edge_obj_[v][j] );    }

     const F& OTo( int v, int j ) const 
     {    return EdgeObject( to_edge_obj_[v][j] );    }

     int IGoes( int v, int w ) const 
     {    for ( int j = 0; j < From(v).isize( ); j++ )
               if ( From(v)[j] == w ) return from_edge_obj_[v][j];
          ForceAssert( 0 == 1 );
          return -1;    }

     const F& OGoes( int v, int w ) const 
     {    for ( int j = 0; j < From(v).isize( ); j++ )
               if ( From(v)[j] == w ) return EdgeObject( from_edge_obj_[v][j] );
          ForceAssert( 0 == 1 );
          return EdgeObject(0);   } // totally invalid but never executed

     const F& EdgeObjectByIndexTo( int v, int j ) const
     {    CheckGoodVertex(v);
          AssertGe( j, 0 );
          AssertLt( j, to_edge_obj_[v].isize( ) );
          return edges_[ to_edge_obj_[v][j] ];    }

     F& EdgeObjectByIndexToMutable( int v, int j )
     {    CheckGoodVertex(v);
          AssertGe( j, 0 );
          AssertLt( j, to_edge_obj_[v].isize( ) );
          return edges_[ to_edge_obj_[v][j] ];    }

     int ITo( int v, int j ) const { return EdgeObjectIndexByIndexTo( v, j ); }

     int EdgeObjectIndexByIndexTo( int v, int j ) const;

     // The following two functions return -1 if they can't find e.
   
     int EdgeObjectIndexToFromIndex( int v, int e ) const;
     int EdgeObjectIndexToToIndex( int v, int e ) const;

     const vec<int>& FromEdgeObj( int v ) const { return from_edge_obj_[v]; }
     const vec<int>& IFrom( int v ) const { return from_edge_obj_[v]; }
     const vec<int>& ToEdgeObj( int v ) const { return to_edge_obj_[v]; }
     const vec<int>& ITo( int v ) const { return to_edge_obj_[v]; }
     vec<int>& FromEdgeObjMutable( int v ) { return from_edge_obj_[v]; }
     vec<int>& ToEdgeObjMutable( int v ) { return to_edge_obj_[v]; }

     const vec< vec<int> >& FromEdgeObj( ) const { return from_edge_obj_; }
     const vec< vec<int> >& ToEdgeObj( ) const { return to_edge_obj_; }
     vec< vec<int> >& FromEdgeObjMutable( ) { return from_edge_obj_; }
     vec< vec<int> >& ToEdgeObjMutable( ) { return to_edge_obj_; }

     const vec<F>& Edges( ) const { return edges_; }
     vec<F>& EdgesMutable( ) { return edges_; }

     

     // ============================================================================
     // ==================== METHODS TO FIND A SUBGRAPH ============================
     // ============================================================================

     // All these methods are const and return a set of vertices or edges in the
     // graph.

     // Find indices of all edges between two vertices.

     vec<int> EdgesBetween( const int v, const int w ) const;

     // Return sorted list of indices of all edges that are between two vertices, 
     // possibly indirectly.

     vec<int> EdgesSomewhereBetween( const int v, const int w ) const;

     // Return sorted list of indices of all edges that lie in the part of the
     // graph bounded by two edges.

     vec<int> EdgesBoundedBy( const int e1, const int e2, 
          const vec<int>& to_left, const vec<int>& to_right ) const;

     // Find indices of all edges between a set v of vertices.  We assume that
     // v is sorted.

     vec<int> EdgesBetween( const vec<int>& v ) const;

     // Find all edges between two vertices.

     vec<F> EdgeObjectsBetween( const int v, const int w ) const;
     
     // Find the 'loop subgraph' of a given graph.  This consists of all edges
     // are part of a loop.  We return a sorted list of the edge indices.  

     void LoopSubgraph( vec<int>& loop_edges ) const;

     // EdgesConnectedTo.  Return all edges that are connected directly or
     // indirectly to a vertex in a given list (which need not be ordered and may 
     // contain duplicates).  Return a sorted list of edge indices.  

     vec<int> EdgesConnectedTo( const vec<int>& v ) const;

     // ============================================================================
     // ================= DELETE OR ADD VERTICES OR EDGES ==========================
     // ============================================================================

     // Add vertices.

     void AddVertices( int nadd );

     // Add an edge from v to w.  Indices of edges from v to any x<=w are untouched,
     // while ones from v to x>w are incremented by 1; likewise for w's edges.
     // Returns the edge number added.

     int AddEdge( const int v, const int w, const F& e );

     // Delete an edge.  Indices of edges < j are untouched, while those > j 
     // decrease by 1.  In particular, adding an edge and then deleting it is
     // guaranteed to leave the indices of other edges unchanged.

     void DeleteEdgeTo( int w, int j );
     void DeleteEdgeFrom( int v, int j );

     // Delete edges entering or exiting a vertex.  Takes argument js, which should
     // be unique-sorted.

     void DeleteEdgesTo( int w, const vec<int>& js );
     void DeleteEdgesFrom( int v, const vec<int>& js );

     // Two functions to delete a bunch of edges.  For both functions, the input
     // list does not have to be sorted, and may contain duplicates.  The first
     // Delete function is O(N), where N is the number of vertices in the graph.
     // The second Delete function is O(to_delete).

     void DeleteEdges( const vec<int>& to_delete );
     void DeleteEdges( const vec<int>& to_delete, const vec<int>& to_left );
     // Fast edge deleting with flags set for each edge directly.
     void DeleteEdgesParallel( const vec<Bool>& to_delete );
     void DeleteEdgesParallel( const vec<int>& to_delete );

     // Remove all edges entering or exiting a given vertex.

     void DeleteEdgesAtVertex( int v );

     // If two edges have the same stop and start, and are same as F-objects,
     // delete one.  This does not actually delete the objects.

     void RemoveDuplicateEdges( );

     // Eliminate unused edge objects.  Renumbers the edges and returns the
     // translation table of old edge ids to new edge ids.

     vec<int> RemoveDeadEdgeObjects( );

     // Remove vertices having no edges coming in or going out, or only specified
     // vertices having this same property.

     void RemoveEdgelessVertices( );
     void RemoveEdgelessVertices( const vec<int>& to_remove );

     // ============================================================================
     // ==============VARIOUS OPERATIONS WITH LEFT/RIGHT UPDATING===================
     // ============================================================================

     void DeleteEdgesWithUpdate( 
          const vec<int>& to_delete, vec<int>& to_left, vec<int>& to_right );
     void TransferEdgesWithUpdate( int v, int w, 
          vec<int>& to_left, vec<int>& to_right, const Bool enter_only = False );
     void AppendWithUpdate( 
          const digraphE<F>& D, vec<int>& to_left, vec<int>& to_right );
     int AddEdgeWithUpdate( const int v, const int w, const F& e, 
          vec<int>& to_left, vec<int>& to_right );
     void SplayVertexWithUpdate( 
          const int v, vec<int>& to_left, vec<int>& to_right );

     // ============================================================================
     // ======================== END ABOVE SECTION =================================
     // ============================================================================

     // InitialEdges: return sorted list of indices of all edges that start from a 
     // source.  TerminalEdges: return sorted list of indices of all edges that 
     // stop at a sink.

     void InitialEdges( vec<int>& v ) const;
     void TerminalEdges( vec<int>& v ) const;

     // GetSuccessors: find all vertices which have a directed path from a vertex
     // in v, and return the minimum distance to each, as a sorted list.

     void GetSuccessors( const vec<int>& v, vec< pair<int,F> >& from_v );

     // Given an edge referred to by To, find the corresponding edge in From,
     // and the other way.

     int InputToOutputFrom( int w, int i ) const;
     int InputFromOutputTo( int w, int i ) const;

     // Change an edge object.

     void SetEdgeObject( int i, const F& e )
     {    edges_[i] = e;    }

     void ChangeEdgeObjectFrom( int v, int i, const F& e );

     // Move an edge's endpoint, without changing the edge object itself.

     void GiveEdgeNewFromVx( int edge_id, int old_from_v, int new_from_v );
     void GiveEdgeNewToVx( int edge_id, int old_to_w, int new_to_w );

     // MinEdge: find the minimum length of an edge from v to w.  Assert if there is
     // no edge.  This only makes sense if Min(F,F) is defined.

     F MinEdge( int v, int w );

     // MaxEdge: find the maximum length of an edge from v to w.  Assert if there is
     // no edge.  This only makes sense if Max(F,F) is defined.

     F MaxEdge( int v, int w );

     // Return equivalence relation e on the edges generated by the rule that if 
     // edge A follows edge B, then they are equivalent.  The classes of e are thus
     // the components of the directed line graph associated to the given graph.
     // If "exclude" is specified, it should have one entry per edge.  The excluded
     // edges are not used in creating the equivalence relation.

     void DualComponentRelation( equiv_rel& e, const vec<Bool>& exclude ) const;

     // ToLeft, ToRight: create vectors that map edge indices to vertex indices.

     void ToLeft( vec<int>& to_left ) const;
     void ToRight( vec<int>& to_right ) const;
     void ToLeftParallel( vec<int>& to_left ) const;
     void ToRightParallel( vec<int>& to_right ) const;
     void ValidateLR( const vec<int>& to_left, const vec<int>& to_right ) const;
     void ValidateLRParallel( 
          const vec<int>& to_left, const vec<int>& to_right ) const;

     // Use Dijkstra's algorithm or the Bellman-Ford algorithm to find the shortest 
     // path in a digraphE<T> (T = int or double) between the 'start' vertex and 
     // the 'stop' vertex, where the length of the path is defined to be the sum of 
     // the edge values.  Negative edge values are allowed but negative cycles will 
     // cause the code to assert.  The handling of ties is not specified here.  
     // Return the list of vertices that defines the path.  In the event of 
     // failure, return an empty path.
     //
     // Implementation based on pseudocode in Wikipedia.  
     //
     // Note that there is a faster algorithm for acyclic graphs, see 
     // Cormen et al. p. 536.
     // 
     // Compare DistancesToEnd.

     void ShortestPath( const int start, const int stop, vec<int>& path ) const;
     
     // SplayVertex.  Replace a given vertex v by To(v).size( ) + From(v).size( )
     // vertices and move the edges entering and exiting v to these vertices.
     // This does not actually touch the edge objects.

     void SplayVertex( const int v );

     // LiberateEdge.  Remove a given edge e: v --> w, and spread out all the
     // edges that touch v or w, so that the two vertices v and w are replaced
     // by From(v) + To(v) + From(w) + To(w) - 2 edges.  (The - 2 is there to
     // account for e's contribution.)

     void LiberateEdge( const int e, const int v, const int w );

     // ComponentsE: find the connected components.  Each component is a list 
     // of edges.  And a version that appears to be stupdendously faster.
     
     void ComponentsE( vec< vec<int> >& comp ) const;
     void ComponentsEFast( vec< vec<int> >& comp ) const; // comp[i] NOT sorted

     // ThisClose: determine if there is a path from v to w whose sum of edge
     // objects is <= d.  This assumes that edge objects can be added and that
     // <= makes sense on them.  In the special case v = w, return True (so long
     // as d >= 0).

     Bool ThisClose( int v, int w, F d ) const;

     // Find the indices of all edges e that form self-loops, i.e.,
     // e goes from v -> v.

     vec<int> SelfLoops( ) const;

     // Determine if the given vertices and edges (sorted lists) comprise the 
     // union of zero or more connected components.

     Bool IsComplete( const vec<int>& vertices, const vec<int>& edges ) const;

     // Determine which edge objects are used.

     void Used( vec<Bool>& used ) const;

     // Return count of number of used edge objects.

     int UsedCount( ) const;

     // Reverse the graph or the connected component containing a given vertex of 
     // it.  The current version does nothing to the edge objects.

     void Reverse( );
     void ReverseComponent( int v );

     // Reverse a vertex.

     void ReverseVertex( const int v )
     {    swap( from_[v], to_[v] );
          swap( from_edge_obj_[v], to_edge_obj_[v] );    }

     // Change order of vertices.

     void ReorderVertices( const vec<int>& new_order );

     // Change order of components.

     void ReorderComponents( const vec<int>& new_order );

     // Find the edges in each connected component.

     void ComponentEdges( vec< vec<int> >& edges ) const;

     // Distance: determine the lengths of all directed paths from vertex v to
     // vertex w, which are no longer than max_dist.  Return answer as D.
     // This is only implemented for the case where F = int.

     void Distance( int v, int w, int max_dist, vec<int>& D ) const;

     // EdgePaths: find all paths from vertex v to w, as lists of edges.  Cycles at 
     // w will not be seen unless you specifically allow terminal loops.  This is
     // asymmetrical and it would have been better to make this the default.  The
     // second version precomputes to_left and to_right, and should be called if 
     // one is invoking EdgePaths many times for the same graph, as in that case it
     // is much faster.
     //
     // WARNING.  The option allow_terminal_loops is probably broken.  It's actually
     // not clear what should be done, if in path exploration one arrives at w,
     // and one would like to allow possible looping from w to w.  Clearly one could
     // allow simple loops (single edges from w to w).  But if one allows multi-edge
     // loops, then that opens up the possibility of exploring way beyond w in the
     // graph.
     //
     // WARNING.  Calling EdgePaths with v = w is probably not a good idea.  The 
     // algorithm will just keep walking beyond w.

     Bool EdgePaths( const int v, const int w, vec< vec<int> >& paths,
          const int max_copies = -1, const int max_paths = -1,
          const int max_iterations = -1, 
          const Bool allow_terminal_loops = False ) const;
     Bool EdgePaths( const vec<int>& to_left, const vec<int>& to_right, const int v,
          const int w, vec< vec<int> >& paths, const int max_copies = -1, 
          const int max_paths = -1, const int max_iterations = -1,
          const Bool allow_terminal_loops = False ) const;

     // EdgePathsLim: find all paths from vertex v to w, as lists of edges, but
     // never use edge e_verboten.  The nature application for this is where e is 
     // an edge emanating from w that one does not want to use, but where there may
     // be cycles between v and w.  We require specification of max_copies, and
     // allow all manner of loops.

     Bool EdgePathsLim( const vec<int>& to_left, const vec<int>& to_right, 
          const int v, const int w, const int e_verboten, vec< vec<int> >& paths, 
          const int max_copies, const int max_paths = -1, 
          const int max_iterations = -1 ) const;

     // AllPathsFixedLength: find all paths from v to w having length L.  This
     // is only implemented for the case where F = int.

     void AllPathsFixedLength( int v, int w, int L, 
          vec< vec<int> >& paths ) const;

     // AllPathsLengthRange: find all paths from v to w having length between L1
     // and L2, as lists of edges.  This is only implemented for the case where 
     // F = int.  If max_paths is specified and more than that many paths are found,
     // return False.  If max_loops is specified and the code loops more than that 
     // number of times, return False.

     Bool AllPathsLengthRange( int v, int w, int L1, int L2, 
          const vec<int>& to_right, vec< vec<int> >& paths, 
          int max_paths = 0, int max_loops = 0, const Bool no_dups = False ) const;

     // AllPathsLengthRangeAlt: same as AllPathsLengthRange, but processes as FIFO
     // instead of LIFO.

     Bool AllPathsLengthRangeAlt( int v, int w, int L1, int L2, 
          const vec<int>& to_right, vec< vec<int> >& paths, 
          int max_paths = 0, int max_loops = 0, const Bool no_dups = False,
          const Bool eq_ok = False, const int max_partials = 0 ) const;

     // SplitEdge( v, j, e1, e2 ): Replace edge from_[v][j] from v to w by a pair 
     // of edges 
     //
     //   e1    e2
     // v --> n --> w,
     //
     // where n is a new vertex [numbered N( ), if N( ) is called before SplitEdge].
     // This pushes e1 and then e2 onto edges_, and leaves an unused entry in 
     // edges_ where the original edge was.

     void SplitEdge( int v, int j, const F& e1, const F& e2 );

     // JoinEdges( x, e ): Suppose vertex x has exactly one edge entering it (from
     // a vertex v) and exactly one edge exiting it (to a vertex w, v != w).  
     // Remove both edges and substitute a single edge from v to w, leaving x as a 
     // vertex having no edges coming in or going out, and also leaving unused 
     // entries in edges_.

     void JoinEdges( int x, const F& e );

     // RemoveUnneededVertices: Apply JoinEdges to all applicable vertices x.  This 
     // assumes that class F has a member function append( const F& d ), used here 
     // to compute the edge e in JoinEdges( x, e ).  This function removes edgeless 
     // vertices, but will leave unused entries in edges_.  If you want to remove
     // them, call RemoveDeadEdgeObjects.

     void RemoveUnneededVertices( );

     // PopBubbles.
     // Input is a set of vertices v.  Each v must be located at the opening of
     // a bubble, with exactly two edges that lead to the same successor w:
     //        _-_
     //  --> v     w -->
     //        -_-
     // Re-route all edges entering v so that they enter w instead.  Then remove v
     // and the bubble edges from the graph.  Last, contract any redundant edges by 
     // a call to RemoveUnneededVertices.

     void PopBubbles( const vec<int> & bubble_vs );

     // PopHyperBubbles.
     // Input is a set of vertices v.  Each v must be located at the opening of
     // a bubble, with two or more edges that lead to the same successor w:
     //         _-_
     //        -   -
     //  --> v ----- w -->
     //        _   _
     //         -_-
     // Re-route all edges entering v so that they enter w instead.  Then remove v
     // and the bubble edges from the graph.  Last, contract any redundant edges by 
     // a call to RemoveUnneededVertices.

     void PopHyperBubbles( const vec<int> & bubble_vs );

     // Append(D): append another digraph, forming the disjoint union.  The
     // new vertices from D are numbered after the existing vertices.

     void Append( const digraphE<F>& D );

     // TransferEdges( v, w ): move all edges entering and exiting v, so that they
     // enter and exit w instead, leaving v without any edges entering and exiting 
     // it.  If enter_only, only transfer entering edges.

     void TransferEdges( int v, int w, const Bool enter_only = False );

     // ContractEdgeTo, ContractEdgeFrom: contract the specified edge.
     // Afterwards that edge no longer exists, and all other edges of
     // the far-end vertex are transfered to v.

     void ContractEdgeTo( int w, int j ) 
     {    int v = to_[w][j];
          DeleteEdgeTo( w, j );
	  if ( v != w ) TransferEdges( v, w );    }
     void ContractEdgeFrom( int v, int j ) 
     {    int w = from_[v][j];
          DeleteEdgeFrom( v, j );
	  if ( v != w ) TransferEdges( w, v );    }

     // Glue.  Given two sequences of edges in the graph, given by subpaths a and b,
     //
     //    e0       em-1            f0      fn-1
     // a0 --> ...  --> am,      b0 --> ... --> bn
     //
     // such that {e0,...,em-1} does not meet {f0,...,fn-1},
     //
     // plus a third sequence of edges not in the graph
     //
     //    f0       fp-1  
     // c0 --> ...  --> cp,  
     // 
     // and two injective maps:
     //
     // EE: {0,...m} --> {0,...p},     FF: {0,...,n} --> {0,...,p},
     //
     // such that EE(0) = 0, EE(m) = p, FF(0) = 0, and FF(n) = p,
     //
     // modify the graph by replacing the "a" and "b" sequences by the "c"
     // sequence, using the maps EE and FF to transfer edges involving any of
     // a0,...,am,b0,...,bn.  However, these vertices are left in the graph as
     // vertices having no edges entering or exiting them.
     //
     // Notes on inputs:
     // The digraph c gives c1,...,cp and the intervening edges.

     void Glue( const EmbeddedSubPath<F>& a, const EmbeddedSubPath<F>& b,
          const vec<int>& EE, const vec<int>& FF, const digraphE<F>& c );

     static size_t externalSizeof() { return 0; }

     // ============================================================================
     // ======================= PRINTERS AND WRITERS ===============================
     // ============================================================================

     // DOT_vel. This is for a digraphE<int>.  Print in DOT with given vertex 
     // labels in circles and edges labeled with the integers they carry.

     void DOT_vel( ostream& out, const vec<String> & vertex_labels );

     // DOT0. This is for a digraphE<int> or digraphE<String>.  Print in DOT 
     // with the edges labeld by the integers or Strings they carry.

     void DOT0( ostream& out );

     // Edge label handler for PrettyDOT, below.

     class edge_label_info {
          public:
          enum ConstructorBehavior { DEFAULT, EXPLICIT, DIRECT };
          edge_label_info( const ConstructorBehavior ctype )
          {    ForceAssert( ctype == DEFAULT );
               label_edges = False;
               edge_labels_base_alpha = False;
               edge_id_names = NULL;
               label_edges_extra = NULL;    }
          edge_label_info( const ConstructorBehavior ctype,
               const Bool label_edges, const Bool edge_labels_base_alpha,
               const vec<String>* label_edges_extra )
               : label_edges(label_edges), 
               edge_labels_base_alpha(edge_labels_base_alpha),
               edge_id_names(NULL), label_edges_extra(label_edges_extra)
          {    ForceAssert( ctype == EXPLICIT );    }
          edge_label_info( const ConstructorBehavior ctype,
               const vec<String>* edge_id_names )
               : label_edges(True), edge_labels_base_alpha(False), 
               edge_id_names(edge_id_names), label_edges_extra(NULL)
          {    ForceAssert( ctype == DIRECT );    }
          Bool label_edges;                     // should label include the edge id?
          Bool edge_labels_base_alpha;          // and if so should it be base-alpha?
          const vec<String>* edge_id_names;     // use this alt edge id
          const vec<String>* label_edges_extra; // tack this onto edge labels
     };

     // PrettyDOT: Like digraph::DOT, but prettier...
     // -- Components are labeled (as "contigs", since this is a holdover from
     //    HyperKmerPath and HyperBasevector)
     // -- Longer edges (as indicated by the lengths vector) are drawn longer
     // -- Edges are colored according to their length:
     //    < 100: gray
     //    100-1000: black
     //    1000-10000: red
     //    > 10000: magenta
     //
     // Edges may be declared 'dashed' (and if so shown that way), or 'invisible' 
     // (and if so not shown at all).  So as to see the boundaries of the visible
     // part of a graph, vertices touching both visible and invisible edges are 
     // shown in red.
     //
     // The 'edge_color' argument (if specified) allows for overriding of the 
     // default color.  Entries for which the color is "" don't override.
     //
     // The 'pen_widths' argument (if specified) allows for overriding of the 
     // default penwidth of the edge.  Entries for which the value <= 0 don't 
     // override.
   
     void PrettyDOT( ostream& out, const vec<double>& lengths,
          const edge_label_info = edge_label_info(edge_label_info::DEFAULT),
          Bool label_contigs = True, Bool label_vertices = False,
	  const vec<int>* componentsToPrint = NULL,
	  const vec<String> *label_contigs_extra = 0,
	  const vec<int>* verticesToPrint = NULL, const vec<Bool>* dashed = NULL,
          const vec<Bool>* invisible = NULL,
          const vec<String>* edge_color = NULL,
          const vec<int>* pen_widths = NULL,
          const String layout = "", const double tiny_top = 100.0,
          const double fontsize = 12, const double scale = 1,
          const double aspect = -1 ) const;
  
     // Method: DumpGraphML
     // Output the digraph structure in a textual format that can be easily
     // read without reference to our code base.

     void DumpGraphML( const String& graphMLFileName ) const;

     void writeBinary( BinaryWriter& writer ) const;
     void readBinary( BinaryReader& reader );

     // ============================================================================
     // ======================= PRIVATE DEFINITIONS ================================
     // ============================================================================

     private:

     vec<F> edges_;
     vec< vec<int> > to_edge_obj_, from_edge_obj_;

};
template <> void digraphE<int>::DOT_vel( ostream& out, 
const vec<String> & vertex_labels );
template <> void digraphE<int>::DOT0( ostream& out );
template <> void digraphE<String>::DOT0( ostream& out );
extern template class digraphE<int>;
extern template class digraphE<VertexPair>;

template <class F> struct Serializability< digraphE<F> >
{ typedef SelfSerializable type; };

// Routines to compare two digraphs.  The == operator checks for complete
// formal equality, whereas EqualExceptEdgeObjectOrder still returns true
// even if the edge objects are not in the same order, so long as 
// FromEdgeObj points to the same objects.

template<class F> bool operator!=( const digraphE<F>& g1, const digraphE<F>& g2 );
template<class F> bool operator==( const digraphE<F>& g1, const digraphE<F>& g2 );
template<class F> 
     void Compare( ostream& out, const digraphE<F>& g1, const digraphE<F>& g2 );

template<class F> bool EqualExceptEdgeObjectOrder( 
     const digraphE<F>& g1, const digraphE<F>& g2 );

// DistancesToEndArr.  For each vertex v, let D(v) be the maximum (possibly infinity)
// of the lengths over all paths that start with v, where length is defined by
// a vector of edge lengths.  Stop counting as soon as length has reached max_dist.  If
// fw = True, do as just stated.  If fw = False, go in reverse direction.

template <class F>
void DistancesToEndArr( digraphE<F> const& G, vec<int> const& edgeLens,
                            int const max_dist, Bool const fw, vec<int>& D );

// same as above, but edge length is defined by a functor
template <class F, class Func>
void DistancesToEndFunc( digraphE<F> const& G, Func func, int const max_dist,
                            Bool const fw, vec<int>& D )
{ vec<int> edgeLens; edgeLens.reserve(G.E());
  for ( F const& edge : G.Edges() ) edgeLens.push_back(func(edge));
  DistancesToEndArr(G,edgeLens,max_dist,fw,D); }

// same as above, but edge length is defined by a member function of the edge type
template<class F> void DistancesToEnd( const digraphE<F>& G,
     int (F::*len)( ) const, const int max_dist, const Bool fw, vec<int>& D )
{ DistancesToEndFunc(G,[len](F const& edge){return (edge.*len)();},max_dist,fw,D); }

// LongestPath.  For a given nonempty acyclic graph, find a longest path from
// a source to a sink.  Do not call this on a graph with cycles, as it may
// either assert or run forever.  The answer is a sequence of edge indices.

template<class F> void LongestPath( const digraphE<F>& G, int (F::*len)( ) const, 
     vec<int>& a_longest_path );

// RemoveHangingEnds.  For each edge e, let M(e) be the maximum (possibly
// infinity) of the lengths over all paths that start with e, where length
// is defined by the function len.  If edges ei, ej, ... emanate from a vertex 
// and M(ei) <= max_del and M(ej)/M(ei) >= min_ratio, delete ei.
//
// Unfortunately this can do "bad" things.  For example in the following graph
//
//                    y        x
//        ------> * <----- * -----> *
//                  -----> 
//
// edge x can be deleted using edge y.

template<class F> void RemoveHangingEnds( digraphE<F>& G,
     int (F::*len)( ) const, const int max_del, const double min_ratio );

// Remove short hanging ends.  Look for
//
//                 x
//                 |
//                 e
//                 |
//        u --c--> v --d--> w
//
// where x is a source or sink, e is short (and can go either way), whereas
// c and d are long.  Works for T = HyperKmerPath, HyperFastavector, or
// HyperBasevector.

template<class T> void RemoveHangingEnds2( T& h, const int max_del,
     const double min_ratio );

// RemoveHangingEnds3 is a lot like RemoveHangingEnds, but in measuring distances
// to end, doesn't count cycles (or more precisely, doesn't count them more than
// once), and hence does not have the defect that RemoveHangingEnds has.  However, 
// for efficiency reasons, RemoveHangingEnds3 needs to be given a bound on the 
// number of paths it will look at, leading off from a given vertex.  Therefore it
// may fail to identify a hanging end because a computation has been truncated.

template<class F> void RemoveHangingEnds3( digraphE<F>& G,
     int (F::*len)( ) const, const int max_del, const double min_ratio,
     const int max_paths );

// =================================================================================
// ============================ DIGRAPHEX CLASS ====================================
// =================================================================================

template<class F> class digraphEX : public digraphX {

     public:

     digraphEX( ) { }

     explicit digraphEX( const digraphE<F>& G );

     // Why isn't this a constructor in digraphE? - because if you do that
     // you need to instantiate digraphEX for all digraphE instantiations.

     digraphE<F> AsDigraphE() const;

     int IFrom( int v, int j ) const
     {    CheckGoodVertex(v);
          AssertGe( j, 0 );
          AssertLt( j, (int) from_edge_obj_[v].size( ) );
          return from_edge_obj_[v][j];    }
     const SerfVec<int>& IFrom( int v ) const
     {    return from_edge_obj_[v];     }

     const F& OFrom( int v, int j ) const 
     {    return EdgeObject( from_edge_obj_[v][j] );    }

     const F& OTo( int v, int j ) const 
     {    return EdgeObject( to_edge_obj_[v][j] );    }

     int IGoes( int v, int w ) const 
     {    for ( int j = 0; j < From(v).isize( ); j++ )
               if ( From(v)[j] == w ) return from_edge_obj_[v][j];
          ForceAssert( 0 == 1 );
          return -1;    }

     // Abut(f,g): does edge target(f) = source(g)?

     Bool Abut( const int f, const int g ) const
     {    return ToRight(f) == ToLeft(g);    }

     const F& OGoes( int v, int w ) const 
     {    for ( int j = 0; j < From(v).isize( ); j++ )
               if ( From(v)[j] == w ) return EdgeObject( from_edge_obj_[v][j] );
          ForceAssert( 0 == 1 );
          return EdgeObject(0);   } // totally invalid but never executed

     int ITo( int v, int j ) const
     {    CheckGoodVertex(v);
          AssertGe( j, 0 );
          AssertLt( j, (int) to_edge_obj_[v].size( ) );
          return to_edge_obj_[v][j];    }
     const SerfVec<int>& ITo( int v ) const
     {    return to_edge_obj_[v];     }

     int E( ) const { return edges_.size( ); }

     int64_t CheckSum() const;

     const MasterVec<F>& Edges( ) const { return edges_; }

     const F& EdgeObject( int i ) const;
     const F& O( int i ) const { return EdgeObject(i); }
     F& OMutable( int i ) { return EdgeObject(i); }

     int ToLeft( int e ) const { return to_left_[e]; }
     int ToRight( int e ) const { return to_right_[e]; }

     Bool SinkEdge( const int e ) const { return From( ToRight(e) ).size( ) == 0; }
     Bool SourceEdge( const int e ) const { return To( ToLeft(e) ).size( ) == 0; }

     vec<int>& to_left() { return to_left_; }
     vec<int>& to_right() { return to_right_; }

     Bool EdgePaths( const int v, const int w, vec< vec<int> >& paths,
          const int max_copies = -1, const int max_paths = -1,
          const int max_iterations = -1, 
          const Bool allow_terminal_loops = False ) const;

     void Used( vec<Bool>& used ) const;
     int UsedCount() const;

     class edge_label_info {
          public:
          enum ConstructorBehavior { DEFAULT, EXPLICIT, DIRECT };
          edge_label_info( const ConstructorBehavior ctype )
          {    ForceAssert( ctype == DEFAULT );
               label_edges = False;
               edge_labels_base_alpha = False;
               edge_id_names = NULL;
               label_edges_extra = NULL;    }
          edge_label_info( const ConstructorBehavior ctype,
               const Bool label_edges, const Bool edge_labels_base_alpha,
               const vec<String>* label_edges_extra )
               : label_edges(label_edges), 
               edge_labels_base_alpha(edge_labels_base_alpha),
               edge_id_names(NULL), label_edges_extra(label_edges_extra)
          {    ForceAssert( ctype == EXPLICIT );    }
          edge_label_info( const ConstructorBehavior ctype,
               const vec<String>* edge_id_names )
               : label_edges(True), edge_labels_base_alpha(False), 
               edge_id_names(edge_id_names), label_edges_extra(NULL)
          {    ForceAssert( ctype == DIRECT );    }
          Bool label_edges;                     // should label include the edge id?
          Bool edge_labels_base_alpha;          // and if so should it be base-alpha?
          const vec<String>* edge_id_names;     // use this alt edge id
          const vec<String>* label_edges_extra; // tack this onto edge labels
     };

     void PrettyDOT( ostream& out, const vec<double>& lengths,
          const edge_label_info = edge_label_info(edge_label_info::DEFAULT),
          Bool label_contigs = True, Bool label_vertices = False,
	  const vec<int>* componentsToPrint = NULL,
	  const vec<String> *label_contigs_extra = 0,
	  const vec<int>* verticesToPrint = NULL, const vec<Bool>* dashed = NULL,
          const vec<Bool>* invisible = NULL,
          const vec<String>* edge_color = NULL,
          const vec<int>* pen_widths = NULL,
          const String layout = "", const double tiny_top = 100.0,
          const double fontsize = 12, const double scale = 1,
          const double ratio = -1 ) const;

     void writeBinary( BinaryWriter& writer ) const;
     void readBinary( BinaryReader& reader );

     private:

     MasterVec<F> edges_;
     VecIntVec to_edge_obj_, from_edge_obj_;
     vec<int> to_left_, to_right_;

};

template <class E> struct Serializability< digraphEX<E> >
{ typedef SelfSerializable type; };

// =================================================================================
// ============================ DIGRAPHEV1 CLASS ===================================
// =================================================================================

// Virtual digraphE class, supporting From but not To.

template<class T> class digraphE_V1 {

     public:

     digraphE_V1( const int N, const function< vec<pair<int,T>> (int)>& from )
     {    N_ = N;
          from_ = from;    }

     int N( ) const { return N_; }

     vec< pair<int,T> > FromVec( const int v ) const
     {    return from_(v);    }

     // Shortest path assumes nonnegative edge values.

     T ShortestPath( const int start, const int stop, vec<int>& path ) const;
     
     private:

     int N_;
     function< vec< pair<int,T> > (int) > from_;

};

// =================================================================================
// ============================ DIGRAPHVE CLASS ====================================
// =================================================================================

template<class V, class E> class digraphVE : public digraphE<E> {

     public:

     digraphVE( ) { }

     digraphVE( const vec< vec<int> >& from, const vec< vec<int> >& to,
          const vec<V>& verts, const vec<E>& edges, 
          const vec< vec<int> >& to_edge_obj,
          const vec< vec<int> >& from_edge_obj );

     void Initialize( const vec< vec<int> >& from, const vec< vec<int> >& to,
          const vec<V>& verts, const vec<E>& edges, 
          const vec< vec<int> >& to_edge_obj,
          const vec< vec<int> >& from_edge_obj );

     digraphVE( const digraphE<E>& G, const vec<V>& verts );

     int N( ) const { return verts_.size( ); }

     const V& Vert( int v ) const;
     V& VertMutable( int v );

     void AddVertex( const V& v );

     void RemoveVertices( const vec<int>& to_remove );

     void writeBinary( BinaryWriter& writer ) const;
     void readBinary( BinaryReader& reader );

     private:

     vec<V> verts_;

};
extern template class digraphVE<int,int>;

template <class V, class E> struct Serializability< digraphVE<V,E> >
{ typedef SelfSerializable type; };

// =================================================================================
// ============================ EMBEDDEDSUBPATH CLASS ==============================
// =================================================================================

// Class Template: EmbeddedSubPath
//
// An EmbeddedSubPath is a directed path within a given <digraphE>.  It keeps the
// address of the digraph, so should not be used in situations where that address
// might change.
//
// Constructor: takes as input a sequence of consecutive edges
//
//    e0       em-1
// a0 --> ...  --> am,
//
// where the edge from ai to ai+1 is given by from_[ a[i] ][ e[i] ].  

template<class E> class EmbeddedSubPath {

     public:

     void TestValid( ) const;

     // Repair e_ entries which are wrong because digraph has been edited.

     void Repair( );

     EmbeddedSubPath( ) { D_ = 0; }

     EmbeddedSubPath( const digraphE<E>& D, const vec<int>& a, const vec<int>& e )
          : D_(&D), a_(a), e_(e)
     {    esafe_.resize( e.size( ) );
          for ( int i = 0; i < e.isize( ); i++ )
               esafe_[i] = D.EdgeObjectIndexByIndexFrom( a[i], e[i] );
          TestValid( );    }

     // The following constructor takes a vector of vertices.  There must be a
     // unique edge from each vertex to the next.

     EmbeddedSubPath( const digraphE<E>& D, const vec<int>& a ) : D_(&D), a_(a)
     {    ForceAssertGt( a.size( ), 0 );
          e_.resize( a.isize( ) - 1 ), esafe_.resize( a.isize( ) - 1 );
          for ( int i = 0; i < e_.isize( ); i++ )
          {    int v = a[i], w = a[i+1];
               Bool found = False;
               for ( int j = 0; j < D.From(v).isize( ); j++ )
               {    if ( D.From(v)[j] == w )
                    {    ForceAssert( !found );
                         found = True;
                         e_[i] = j;
                         esafe_[i] = D.EdgeObjectIndexByIndexFrom( v, j );    }    }
               if ( !found )
               {    cout << "No edge found from vertex " << v << " to vertex "
                         << w << ".\n";
                    ForceAssert( 0 == 1 );    }    }    }

     int NVertices( ) const { return a_.size( ); }
     int NEdges( ) const { return e_.size( ); }

     int Vertex( int i ) const { return a_[i]; }
     int FirstVertex( ) const { return a_.front( ); }
     int LastVertex( ) const { return a_.back( ); }

     bool IsFromSource( ) const { return D_->Source(FirstVertex()); }
     bool IsToSink( ) const { return D_->Sink(LastVertex()); }

     const E& EdgeObject( int i ) const
     {    AssertGe( i, 0 );
          AssertLt( i, e_.isize( ) );
          return D_->EdgeObjectByIndexFrom( a_[i], e_[i] );    }

     int EdgeObjectIndexAbs( int i ) const
     {    return esafe_[i];    }

     int EdgeObjectIndex( int i ) const
     {    AssertGe( i, 0 );
          AssertLt( i, e_.isize( ) );
          return D_->EdgeObjectIndexByIndexFrom( a_[i], e_[i] );    }

     int EdgeObjectFromIndex( int i ) const
     {    AssertGe( i, 0 );
          AssertLt( i, e_.isize( ) );
          return e_[i];    }

     void SetVertex( int i, int v ) { a_[i] = v; }
     void SetEdge( int i, int e ) 
     {    e_[i] = e; 
          esafe_[i] = D_->EdgeObjectIndexByIndexFrom( a_[i], e );    }

     // Add vertex and edge to left or right.

     void Prepend( int a, int e )
     {    a_.push_front(a);
          e_.push_front(e);
          esafe_.push_front( D_->EdgeObjectIndexByIndexFrom(a,e) );    }
     void Append( int a, int e )
     {    esafe_.push_back( D_->EdgeObjectIndexByIndexFrom(a_.back(),e) );
          a_.push_back(a);
          e_.push_back(e);    }
  
     // Return true if this edge already exists in this EmbeddedSubPath.
     // If true, calling Prepend( a, e ) or Append( a, e ) would produce a loop
     // and/or a duplicated edge.
 
     bool Contains( int edge_id ) 
     {    for ( int i = 0; i < e_.isize( ); i++ )
               if ( edge_id == esafe_[i] ) return true;
          return false;     }

     // Return true if there are any edges in common between p1 and p2.

     friend Bool HasSharedEdge( const EmbeddedSubPath& p1, 
          const EmbeddedSubPath& p2 )
     {    vec<int> edges1, edges2;
          edges1.reserve( p1.NEdges( ) );
          edges2.reserve( p2.NEdges( ) );
          for ( int i = 0; i < p1.NEdges( ); i++ )
          {    edges1.push_back( p1.D_-> 
                    EdgeObjectIndexByIndexFrom( p1.a_[i], p1.e_[i] ) );    }
          for ( int i = 0; i < p2.NEdges( ); i++ )
          {    edges2.push_back( p2.D_-> 
                    EdgeObjectIndexByIndexFrom( p2.a_[i], p2.e_[i] ) );    }
          Sort(edges1), Sort(edges2);
          return Meet( edges1, edges2 );    }

     friend bool operator== ( const EmbeddedSubPath& lhs, 
          const EmbeddedSubPath& rhs ) 
     {    return ( lhs.D_ == rhs.D_ && lhs.a_ == rhs.a_ && lhs.e_ == rhs.e_ &&
                lhs.esafe_ == rhs.esafe_ );    }

     private:

     const digraphE<E>* D_;
     vec<int> a_;
     vec<int> e_;
     vec<int> esafe_;

};

#define DIGRAPH_INVALID(REASON,EXIT) \
  { std::cout << "\nDigraph is invalid.\n" << REASON << "\nAbort." << std::endl; \
  if (EXIT) TracebackThisProcess(); \
  else return False; }

#define DIGRAPH_INVALIDP(REASON) \
  { std::cout << "\nDigraph is invalid.\n" << REASON << "\nAbort." << std::endl; \
  TracebackThisProcess();  }

void PrintStandardDOTHeader( ostream& out );

template<class E> void FindLeftMostVertex( const digraphE<E>& G,
     const vec<double>& lengths, const vec<int>& o, const vec<Bool>* invisible,
     int& leftv );

void PrintDotLegends( 
     ostream& out, const vec<vec<String>>& legends, const vec<String>& colors );

#endif
