// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

// Define gap types.

#ifndef TENX_GAP_H
#define TENX_GAP_H

#include "CoreTools.h"
#include "paths/HyperBasevector.h"

// We represent a "supergraph" assembly D as a graph whose edges are vectors of
// integers.  Each such path represents a path in another fixed graph hb.
// However certain edges are special "gap" edges.  Gap edges are denoted by 
// vectors whose first entry is negative.
//
// Categorization of special edges (which will expand over time).

// 1. Read pair captured gaps.
//    representation: (-1).

inline Bool IsPairGap( const vec<int>& x ) { return x[0] == -1; }

// 2. Barcode-only gaps.
//    representation: (-2) or (-2, gap) where gap is the predicted gap size.

inline Bool IsBarcodeOnlyGap( const vec<int>& x ) { return x[0] == -2; }

// 3. Sequence of length n >= K.  The interpretation of this is that it should 
//    overlap abutting edges by K-1 bases.  However these are to first be trimmed by 
//    ltrim and rtrim bases.
//    representation (-3, ltrim, rtrim, n, encoding of basevector).
//
//    Notes.
//    1. A sequence gap edge must be abutted on the left by exactly one edge,
//       and likewise on the right.
//    2. We allow multiple parallel sequence edges between two vertices.
//       They must all have the same ltrim and rtrim values.

inline Bool IsSequence( const vec<int>& x ) { return x[0] == -3; }

void SeqToGap( const int ltrim, const int rtrim, const basevector& x, vec<int>& y );
void GapToSeq( const vec<int>& y, int& ltrim, int& rtrim, basevector& x );

// 4. General cells.
//    representation (-4, stuff encoding the cell).
//    This replaces a pretty arbitrary cell.  It is intended for use with cycles,
//    but some cyclic cases have snuck through.  Perhaps need to modify 
//    CaptureMessyLoops to exclude them.

inline Bool IsCell( const vec<int>& x ) { return x[0] == -4; }

void ReinsertLoop( const int d, const HyperBasevectorX& hb, const vec<int>& inv,
     digraphE<vec<int>>& D, vec<int>& dinv, vec<int>& to_left, vec<int>& to_right );

void ReinsertLoopsMap( digraphE<vec<int>> const& D, vec<int> const& dinv, 
     std::map<int,int>& edge_map );

void ReinsertLoops( const HyperBasevectorX& hb, const vec<int>& inv,
     digraphE<vec<int>>& D, vec<int>& dinv );

// this is a copy of ReinsertLoops above that is to be used by
// ValidateMakeFasta. It enables tracking of cell index.
void ReinsertLoopsWithTracking( const HyperBasevectorX& hb, const vec<int>& inv,
     digraphE<vec<int>>& D, vec<int>& dinv , vec <int> & ctracker );

void Encode( const digraphE<vec<int>>& D, vec<int>& x );
void Decode( const vec<int>& x, int& pos, digraphE<vec<int>>& D );

// A cell represents an 'abstracted' subgraph of a supergraph, that has defined
// entry and exit vertices (possibly equal).

class cell {

     public:

     cell( ) { }

     cell( const digraphE<vec<int>>& G, const int v, const int w )
          : G_(G), v_(v), w_(w) { }

     const digraphE<vec<int>>& G( ) const { return G_; }
     digraphE<vec<int>>& GMutable( ) { return G_; }
     int Left( ) const { return v_; }
     int Right( ) const { return w_; }

     // Conversion to vec<int> and back.

     void CellEncode( vec<int>& x ) const;
     void CellDecode( const vec<int>& x );

     // FindPath.  Find a path through a cell, which is represented as a sequence
     // of edges p in G.  Tries to find a path that goes through as many edges as 
     // possible, but this is not guaranteed.  With extremely low frequency (not 
     // observed), can return path of length zero, which should be treated as a gap.

     void FindPath( vec<int>& p ) const;

     private:

     digraphE<vec<int>> G_;    // the subgraph
     int v_, w_;               // entry and exit vertices in the subgraph

};

// This does not test everything that should be tested.

void ValidateGapEdges( const HyperBasevectorX& hb, const vec<int>& inv,
     const digraphE<vec<int>>& D, const vec<int>& dinv );

// Munch.  Take as input a supergraph assembly D and a base edge vector tigs.
// Modify these in an unholy fashion by, for each sequence gap edge:
// * push back onto tigs, creating a new "base" graph edge, which is not really
//   a base graph edge;
// * trim abutting supergraph edges according to the instructions associated to
//   the seq graph edge -- this will sometimes involve creating new "base" graph
//   edges in tigs;
// * edit edges in D accordingly, converting the seq graph edges to edges of the
//   form {e} where e is one of the newly created pseudo-base-graph edges, and
//   similarly making modifications to trimmed edges.
// This process can fail in rare cases if the trimming instructions are 
// contradictory.  The number of failures is reported.
//
// Code largely borrowed from SuperToSeqGraph.
//
// If we wanted to subsequently modify hb itself, this could would be a good place
// to start.

void Munch( digraphE<vec<int>>& D, const vec<int>& dinv, vecbasevector& tigs, 
     vec<int>& inv, const int HBK );

#endif
