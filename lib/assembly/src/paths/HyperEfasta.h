///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Define class HyperEfasta, which is a directed graph whose edges are
// efastas, with overlaps of fixed size K.

#ifndef HYPER_EFASTA_H
#define HYPER_EFASTA_H

#include "Basevector.h"
#include "efasta/EfastaTools.h"
#include "feudal/BinaryStream.h"
#include "graph/Digraph.h"
#include "paths/HyperBasevector.h"

extern template class digraphE<efasta>;

class HyperEfasta : public digraphE<efasta> {

     public:

     HyperEfasta( ) { }
     HyperEfasta( int K ) : K_(K) { }
     HyperEfasta( const HyperBasevector& h );
     void Initialize( const HyperBasevector& h );

     int K( ) const { return K_; }
     void SetK( int K ) { K_ = K; }

     int EdgeLengthBases( int e ) const { return EdgeObject(e).Length1( ); }
     int EdgeLengthKmers( int e ) const 
     { return EdgeObject(e).Length1( ) - K( ) + 1; }

     void RemoveUnneededVertices( );

     void writeBinary( BinaryWriter& writer ) const;
     void readBinary( BinaryReader& reader );
     static size_t externalSizeof() { return 0; }

     void PrintSummaryDOT0w( ostream& out, Bool label_contigs = True,
                             Bool label_vertices = False, Bool label_edges = False,
                             const vec<int>* componentsToPrint = NULL,
                             const Bool edge_labels_base_alpha = False,
                             const vec<String> *label_edges_extra = NULL,
                             const vec<String> *label_contigs_extra = NULL,
                             const vec<int> *verticesToPrint = NULL ) const;

     friend ostream& operator<<( ostream& out, const HyperEfasta& h )
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



vec<String> SpecialEfastaTag2Graph( digraphE<basevector>& out, efasta in);


SELF_SERIALIZABLE(HyperEfasta);

#endif
