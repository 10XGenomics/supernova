///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Define class HyperBasevector, which is a directed graph whose edges are
// basevectors, and other ancillary classes.

#ifndef HYPER_BASEVECTOR_H
#define HYPER_BASEVECTOR_H

#include "Basevector.h"
#include "Equiv.h"
#include "Qualvector.h"
#include "feudal/BinaryStream.h"
#include "graph/Digraph.h"
#include "paths/HyperKmerPath.h"
#include "paths/KmerBaseBroker.h"
#include "paths/KmerPath.h"

// Class: HyperBasevector
// 
// A HyperBasevector is a <kmer numbering>-independent representation of 
// a HyperKmerPath.  However, it is not independent of K.

class HyperBasevectorX;  // forward declaration

class HyperBasevector;

extern template class digraphE<basevector>;
class HyperBasevector : public digraphE<basevector> {

     public:
     using Superclass = digraphE<basevector>;

     HyperBasevector( ) { K_ = 0; }
     HyperBasevector( int K ) : K_(K) { }

     // Constructor from a HyperKmerPath having no gaps:

     HyperBasevector( const HyperKmerPath& h, const KmerBaseBroker& kbb );
     HyperBasevector( const HyperKmerPath& h, const KmerBaseBrokerBig& kbb );

     void Initialize( const HyperKmerPath& h, const KmerBaseBroker& kbb );

     // Constructor: given a collection of basevectors, create a graph having
     // one edge and two vertices for each of the edge objects

    HyperBasevector( int K, const vec<basevector>& p )
          : digraphE<basevector>( p, EDGES_SEPARATE )
     {    K_ = K;    }

     // Constructor: form the disjoint union of some HyperBasevectors.
     
     HyperBasevector( int K, const vec<HyperBasevector>& v )
     {    SetK(K);   
          SetToDisjointUnionOf(v);    }

     explicit HyperBasevector( const HyperBasevectorX& hbx );

     void Initialize( int K, const vec< vec<int> >& from, 
          const vec< vec<int> >& to, const vec<basevector>& edges, 
          const vec< vec<int> >& from_edge_obj, const vec< vec<int> >& to_edge_obj )
     {    K_ = K;
          from_ = from;
          to_ = to;
          FromEdgeObjMutable( ) = from_edge_obj;
          ToEdgeObjMutable( ) = to_edge_obj;
          EdgesMutable( ) = edges;    }

     // Constructor: given a HyperBasevector, and given a list of vertex indices,
     // create the HyperBasevector having those vertices (with indices starting at 0,
     // but in the given order), and having all the edges that were between those
     // vertices.  Thus this is a "complete subgraph" constructor.

     HyperBasevector( const HyperBasevector& h, const vec<int>& v )
          : digraphE<basevector>(h.Subgraph(v))
     {    K_ = h.K( );    }

     explicit HyperBasevector( const String& filename );

     // Constructor: from K, a digraphE g, and some edge objects.  This ignores 
     // the edge objects of g and puts in the new edge objects.

     template<class T>
     HyperBasevector( int K, const digraphE<T>& g, const vec<basevector>& edges )
          : digraphE<basevector>( g, edges )
     {    K_ = K;    }

     // Check that the edge adjacencies in this HyperBasevector make sense.
     // If they don't, this causes a FatalErr.

     Bool TestValid( const Bool exit = True ) const;

     // Return a checksum.

     int64_t CheckSum( ) const;
  
     // SetToDisjointUnionOf: clear a given HyperBasevector and set it to the 
     // disjoint union of a given collection of HyperBasevectors.  If clear is set to
     // false, then we don't clear first and form the disjoint union WITH the current HyperBasevector.

     void SetToDisjointUnionOf( const vec<HyperBasevector>& v, const bool clear = true );

     // AppendDisjointUnionOf: efficiently form the disjoint union of v WITH the current HyperBasevector.
     // Calls SetToDisjointUnionOf to do the dirty work.

     void AppendDisjointUnionOf( const vec<HyperBasevector>& v ) { SetToDisjointUnionOf( v, false ); }

     // EdgePathToBases.  Given a sequence of edge ids, return the associated
     // basevector.

     basevector EdgePathToBases( const vec<int>& e ) const;

     int K( ) const { return K_; }
     void SetK( int K ) { K_ = K; }

     // Cat: return concatenation of edges, with K-1 bases trimmed between.

     basevector Cat( const int e1, const int e2 ) const;
     basevector Cat( const int e1, const int e2, const int e3 ) const;
     basevector Cat( const int e1, const int e2, const int e3, const int e4 ) const;
     basevector Cat( const int e1, const int e2, const int e3, const int e4,
          const int e5 ) const;
     basevector Cat( const vec<int>& e ) const;

     // Get involution of a HyperBasevector.  If it doesn't exist you'll get
     // garbage.  Parallel.

     void Involution( vec<int>& inv );

     void LowerK( int newK );

     // ReduceK is like LowerK but smarter.  It is consistent with an involution.

     void ReduceK( int newK );

     void RemoveUnneededVertices( );

     int EdgeLength( int e ) const { return EdgeObject(e).size( ); }
     int EdgeLengthBases( int e ) const { return EdgeObject(e).size( ); }
     int Bases( int e ) const { return EdgeObject(e).size( ); }

     int EdgeLengthKmers( int e ) const { return EdgeObject(e).isize( ) - K( ) + 1; }
     int Kmers( int e ) const { return EdgeObject(e).isize( ) - K( ) + 1; }
     int KmerSum( const vec<int>& p ) const
     {    if ( p.empty( ) || p[0] < 0 ) return 0;
          int sum = 0;
          for ( auto e : p ) sum += Kmers(e);
          return sum;    }

     // TotalEdgeLength: return the total number of bases.

     longlong TotalEdgeLength() const {
       longlong length = 0;
       for ( int e = 0; e < EdgeObjectCount(); e++ )
	 length += EdgeObject(e).size();
       return length;
     }

     // Reverse entire graph.

     void Reverse( );

     // If two components are reverse complements of each other, delete one of them.

     void DeleteReverseComplementComponents( );

     // Generate a vecbasevector consisting of the edges.  Note that these are
     // not in the natural order, but are instead ordered by component.

     void GenerateVecbasevector( vecbvec &bases ) const;

     void RemoveSmallComponents( int min_kmers );

     // FindIsomorphicComponents2.  First set comp1 to the component relation on
     // the vertices of the HyperBasevector.  Then extract representatives of the
     // orbits (using OrbitRepsAlt), and set comp2 equal to an equivalence relation
     // on these corresponding to the isomophism of components in the 
     // HyperBasevector.
     
     void FindIsomorphicComponents2( equiv_rel& comp1, equiv_rel& comp2 ) const;

     // Write this HyperBasevector to file.
     // default behavior for FASTA is to name edges edge_A, edge_B, etc.
     // if edge_alpha == False, they'll be written out as 0, 1, 2

     void DumpFastb( const String& fn ) const;
     void DumpFasta( const String& fn, Bool edge_alpha = True ) const;

     void PrintSummaryDOT0w( ostream& out, Bool label_contigs = True,
			     Bool label_vertices = False, Bool label_edges = False,
			     const vec<int>* componentsToPrint = NULL,
                             const Bool edge_labels_base_alpha = False,
			     const vec<String> *label_edges_extra = NULL,
			     const vec<String> *label_contigs_extra = NULL,
			     const vec<int> *verticesToPrint = NULL,
			     const vec<String> *edge_color = NULL ) const;

     void writeBinary( BinaryWriter& writer ) const;
     void readBinary( BinaryReader& reader );
     static size_t externalSizeof() { return 0; }

     friend Bool operator==( const HyperBasevector& h1, const HyperBasevector& h2 );

     friend ostream& operator<<( ostream& out, const HyperBasevector& h )
     {    for ( int v = 0; v < h.N( ); v++ )
          {    for ( size_t j = 0; j < h.From(v).size(); j++ )
               {    int e = h.EdgeObjectIndexByIndexFrom( v, j );
                    int w = h.From(v)[j];
                    h.EdgeObject(e).Print( out, BaseAlpha(e) + " [vert_" 
                         + ToString(v) + "-->vert_" + ToString(w) 
                         + "]" );    }    }
          return out;    }

     private:

     int K_;

};

class HyperBasevectorX : public digraphEX<basevector> {

     public:
     using Superclass = digraphEX<basevector>;

     HyperBasevectorX( ) { }
     HyperBasevectorX( int K ) : K_(K) { }

     explicit HyperBasevectorX( const HyperBasevector& hb )
          : digraphEX<basevector>( (const digraphE<basevector>&) hb )
     {    K_ = hb.K( );    }

     int K( ) const { return K_; }

     int Bases( int e ) const { return EdgeObject(e).size( ); }
     int Kmers( int e ) const { return EdgeObject(e).isize( ) - K( ) + 1; }

     basevector Cat( const vec<int>& e ) const;

     int64_t CheckSum( ) const;

     void writeBinary( BinaryWriter& writer ) const;
     void readBinary( BinaryReader& reader );

     private:

     int K_;

};

// ComponentsAreIsomorphic: see HyperKmerPath.h.

bool ComponentsAreIsomorphic( const HyperBasevector& hkp1, int seed_vx_1,
			      const HyperBasevector& hkp2, int seed_vx_2,
			      vec<int>* p_isomorphism = NULL );

SELF_SERIALIZABLE(HyperBasevector);
SELF_SERIALIZABLE(HyperBasevectorX);

#endif
