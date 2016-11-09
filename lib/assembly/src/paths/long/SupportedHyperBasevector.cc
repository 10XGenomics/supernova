///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"
#include "ParallelVecUtilities.h"
#include "Set.h"
#include "VecUtilities.h"
#include "fastg/FastgGraph.h"
#include "math/Functions.h"
#include "paths/BigMapTools.h"
#include "paths/HyperBasevector.h"
#include "paths/long/CleanEfasta.h"
#include "paths/long/DigraphFromWords.h"
#include "paths/long/LargeKDispatcher.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/MakeKmerStuff.h"
#include "paths/long/PairInfo.h"
#include "paths/long/SupportedHyperBasevector.h"

namespace { // open anonymous namespace

vec<int> GetNextEdges( const HyperBasevector& hb, const vec< vec<int> >& paths,
     const vec<fix64_6>& counts_fw, const vec<fix64_6>& counts_rc, 
     const vec< vec< pair<int,int> > >& paths_index,
     const int median_read, const vec<int>& p, const fix64_6 keep_score, 
     const int win_ratio, const int verbosity, ostringstream& out )
{    
     // First find the edges that follow p, from an overlap of a path with p.

     vec<int> nexts;
     int b = p.back( );
     ForceAssertGe( b, 0 );
     for ( int k = 0; k < paths_index[b].isize( ); k++ )
     {    int id2 = paths_index[b][k].first, pos2 = paths_index[b][k].second;
          if ( pos2 == paths[id2].isize( ) - 1 ) continue;
          int n = paths[id2][pos2+1];
          if ( n < 0 || Member( nexts, n ) ) continue;
          Bool mismatch = False;
          for ( int p2 = pos2 - 1; p2 >= 0; p2-- )
          {    int p1 = p.isize( ) - 1 - ( pos2 - p2 );
               if ( p1 < 0 ) break;
               if ( p[p1] != paths[id2][p2] )
               {    mismatch = True;
                    break;    }    }
          if ( !mismatch ) nexts.push_back( paths[id2][pos2+1] );    }
     Sort(nexts);
     if ( nexts.size( ) > 1 )
     {    
          // Try to eliminate some nexts using the nonexistence criterion.

          if ( verbosity >= 2 )
          {    out << "- checking " << printSeq( p, " " );
               out << " { " << printSeq( nexts, " " ) << " }\n";    }
          for ( int s = 0; s < p.isize( ); s++ )
          {    vec<fix64_6> score( nexts.size( ), 0 );
               int len = 0;
               for ( int j = s + 1; j < p.isize( ); j++ )
               {    ForceAssertGe( p[j], 0 );
                    len += hb.EdgeLengthKmers( p[j] );    }
               if ( len > median_read ) continue;
               vec<int> p1, offsets;
               p1.SetToSubOf( p, s, p.isize( ) - s );
               for ( int l = 0; l < nexts.isize( ); l++ )
               {    vec<int> p2(p1);
                    p2.push_back( nexts[l] );

                    // Compute the number of reads that contain p2.

                    vec<int> ids;
                    int b = p2.back( );
                    for ( int k = 0; k < paths_index[b].isize( ); k++ )
                    {    int id3 = paths_index[b][k].first; 
                         int pos3 = paths_index[b][k].second;
                         if ( pos3 < p2.isize( ) - 1 ) continue;
                         Bool mismatch = False;
                         for ( int p3 = pos3 - 1; p3 >= 0; p3-- )
                         {    int p2x = p2.isize( ) - 1 - ( pos3 - p3 );
                              if ( p2x < 0 ) break;
                              if ( p2[p2x] != paths[id3][p3] )
                              {    mismatch = True;
                                   break;    }    }
                         if (mismatch) continue;
                         ids.push_back(id3);    }
                    UniqueSort(ids);
                    for ( int i = 0; i < ids.isize( ); i++ )
                         score[l] += counts_fw[ ids[i] ] + counts_rc[ ids[i] ];    }
               if ( verbosity >= 3 )
               {    PRINT2_TO( out, s, len );
                    out << "score = " << printSeq( score, " " ) << "\n";    }
               vec<Bool> to_delete( nexts.size( ), False );
               for ( int l = 0; l < nexts.isize( ); l++ )
               {    if ( score[l] < keep_score 
                         && Max(score) >= win_ratio * score[l] ) 
                    {    to_delete[l] = True;     }    }
               EraseIf( nexts, to_delete );    }    }
     return nexts;    }

vec<int> Betweens( const HyperBasevector& hb, const vec<int>& to_right, 
     const vec<int>& q, const bool debug = false )
{    vec<int> d, v;
     d.push_back( q.front( ), q.back( ) );
     digraphE<basevector> G1(hb);
     G1.DeleteEdges(d);
     for ( int j = 0; j < q.isize( ) - 1; j++ )
          v.push_back( to_right[ q[j] ] );

     if ( debug ) {
	  cout << "BETWEENS: v=" << printSeq(v) << endl;
     }
     return G1.EdgesConnectedTo(v);    }

vec<int> Betweens2(const HyperBasevector& hb, const vec<int>& to_left, const vec<int>& to_right,
	  const vec<int>& q, const bool debug = false)
{
     vec<int> v;
     typedef pair<int, int > VertexPair;
     vec<VertexPair > v_verboten(2);		// edges between any VertexPair here will be not be explored
						// this is how we isolate the box.  Size is hardcoded.

     // all interior vertices on v
     for (int j = 0; j < q.isize() - 1; j++)
	  v.push_back(to_right[q[j]]);

     // connections from the exterior to the interior are on the "verboten" list
     // we've hard-coded the size as the tests below are hard-coded for speed.
     v_verboten[0] = VertexPair( to_left[ q.front() ], to_right[ q.front() ] );
     v_verboten[1] = VertexPair( to_left[ q.back() ], to_right[ q.back() ] );

     // find all vertices that can be reached from vertices in v, without crossing one of the
     // verboten edges.  This is a modified version of VerticesConnectedTo.
     vec<int> v_neighbors( v );
     UniqueSort(v_neighbors);
     while ( 1 ) {
	  size_t n1 = v_neighbors.size();
	  for ( size_t i = 0; i < n1; i++ ) {

	       const vec<int>& from = hb.From( v_neighbors[i] );	// from this vertex to others
	       for ( size_t j = 0; j < from.size(); ++j ) {
		   VertexPair test( v_neighbors[i], from[j] );
		   if ( test != v_verboten[0] && test != v_verboten[1] )
			 v_neighbors.push_back( from[j] );
	       }
	       const vec<int>& to = hb.To( v_neighbors[i] );		// to this vertex from others
	       for ( size_t j = 0; j < to.size(); ++j ) {
		    VertexPair test( to[j], v_neighbors[i] );
		    if ( test != v_verboten[0] && test != v_verboten[1] )
			 v_neighbors.push_back( to[j] );
	       }
	  }
	  UniqueSort(v_neighbors);
	  if (v_neighbors.size() == n1)
	       break;

     }

     // debugging cruft
     if ( debug ) {
	  cout << "q=" << printSeq(q) << endl;
	  cout << "to_left[q]=" << to_left[q[0]] << "," << to_left[q[1]] << endl;
	  cout << "to_right[q]=" << to_right[q[0]] << "," << to_right[q[1]] << endl;
	  cout << "v_neighbors=" << printSeq(v_neighbors) << endl;
	  for ( size_t i = 0; i < v_neighbors.size(); ++i ) {
	       cout << "vneighbors[" << i << "].From()=" << printSeq(hb.From(v_neighbors[i])) << endl;
	       cout << "vneighbors[" << i << "].To()=" << printSeq(hb.To(v_neighbors[i])) << endl;
	  }
	  vec<int> dels = Betweens( hb, to_right, q, true );
	  cout << "dels = " << printSeq(dels) << endl;
     }

     // Again, a modified version of EdgesConnectedTo: find all edges connected to a list of vertices, excluding
     // edges between verboten pairs.
     vec<int> e;
     for ( size_t i = 0; i < v_neighbors.size(); i++ ) {
	  for ( size_t j = 0; j < hb.From( v_neighbors[i] ).size(); j++ ) {
	       int x = hb.EdgeObjectIndexByIndexFrom( v_neighbors[i], j );
	       VertexPair test(v_neighbors[i], hb.From(v_neighbors[i])[j] );
	       if ( test != v_verboten[0] && test != v_verboten[1]  )
		    e.push_back(x);
	  }
	  for ( size_t j = 0; j < hb.To( v_neighbors[i] ).size(); j++ ) {
	       int x = hb.EdgeObjectIndexByIndexTo( v_neighbors[i], j );
	       VertexPair test( hb.To(v_neighbors[i])[j], v_neighbors[i] );
	       if ( test != v_verboten[0] && test != v_verboten[1] )
		    e.push_back(x);
	  }
     }
     UniqueSort(e);
     return e;
}

class seq_place {

     public:

     seq_place( ) { }
     seq_place( const vec<int>& x, const int g, const Bool fw, const int start,
          const int stop ) : x(x), g(g), fw(fw), start(start), stop(stop) { }

     vec<int> x;    // seq is s[x][0], s[x[1]], ...
     int g;         // placed at G[g]
     Bool fw;       // placement orientation
     int start;     // start of placement
     int stop;      // stop of placement

};

} // close anonymous namespace

Bool Overlap( const vec<int>& v, const vec<int>& w, const int o )
{    for ( int i = 0; i < v.isize( ); i++ )
     {    int j = i - o;
          if ( j < 0 || j >= w.isize( ) ) continue;
          if ( v[i] != w[j] ) return False;    }
     return True;    }

Bool PullApartProcessVertex( SupportedHyperBasevector& shb, vec<int>& to_left, 
     vec<int>& to_right, vec< vec< pair<int,int> > >& paths_index, 
     const int v, const int w, const int pass, const double min_weight_split, 
     const long_logging& logc )
{    int id1 = ( pass == 1 ? 0 : 1 ), id2 = ( pass == 1 ? 1 : 0 );
     int x1 = shb.EdgeObjectIndexByIndexTo( v, id1 );
     int x2 = shb.EdgeObjectIndexByIndexTo( v, id2 );
     int r = shb.EdgeObjectIndexByIndexFrom( v, 0 );
     int y1 = shb.EdgeObjectIndexByIndexFrom( w, 0 );
     int y2 = shb.EdgeObjectIndexByIndexFrom( w, 1 );
     if ( x1 == y1 || x1 == y2 || x2 == y1 || x2 == y2 ) return False;

     fix64_6 weight_11 = 0, weight_12 = 0, weight_21 = 0, weight_22 = 0;
     vec<int> p11, p12, p21, p22;
     p11.push_back(x1,r,y1), p12.push_back(x1,r,y2);
     p21.push_back(x2,r,y1), p22.push_back(x2,r,y2);

     vec<int> p11_paths, p12_paths, p21_paths, p22_paths;
     #pragma omp parallel for
     for ( int l = 0; l < paths_index[x1].isize( ); l++ )
     {    int i = paths_index[x1][l].first, j = paths_index[x1][l].second;
          const vec<int>& p = shb.Path(i);
          if ( p.Contains( p11, j ) )
          {
               #pragma omp critical
               {    p11_paths.push_back(i);    }    }
          if ( p.Contains( p12, j ) )
          {
               #pragma omp critical
               {    p12_paths.push_back(i);    }    }    }
     #pragma omp parallel for
     for ( int l = 0; l < paths_index[x2].isize( ); l++ )
     {    int i = paths_index[x2][l].first, j = paths_index[x2][l].second;
          const vec<int>& p = shb.Path(i);
          if ( p.Contains( p21, j ) )
          {
               #pragma omp critical
               {    p21_paths.push_back(i);    }    }
          if ( p.Contains( p22, j ) )
          {
               #pragma omp critical
               {    p22_paths.push_back(i);    }    }    }
     UniqueSort(p11_paths), UniqueSort(p12_paths);
     UniqueSort(p21_paths), UniqueSort(p22_paths);
     for ( int l = 0; l < p11_paths.isize( ); l++ )
     {    int i = p11_paths[l];
          weight_11 += shb.Weight(i);    }
     for ( int l = 0; l < p12_paths.isize( ); l++ )
     {    int i = p12_paths[l];
          weight_12 += shb.Weight(i);    }
     for ( int l = 0; l < p21_paths.isize( ); l++ )
     {    int i = p21_paths[l];
          weight_21 += shb.Weight(i);    }
     for ( int l = 0; l < p22_paths.isize( ); l++ )
     {    int i = p22_paths[l];
          weight_22 += shb.Weight(i);    }

     const int min_weight_split_low = 2;
     Bool OK = False;
     if ( weight_11 >= min_weight_split_low && weight_22 >= min_weight_split_low
          && weight_12 == 0 && weight_21 == 0 )
     {    OK = True;    }
     if ( weight_11 >= min_weight_split && weight_22 >= min_weight_split
          // && weight_12 <= 1 && weight_21 <= 1
          && weight_12 + weight_21 <= 2
          && shb.EdgeLengthKmers(r) <= shb.MedianCorrectedReadLengthFudge( ) )
     {    OK = True;    }
     if ( weight_11 >= min_weight_split/2
          && weight_22 >= min_weight_split/2
          && weight_12 + weight_21 < 2
          && shb.EdgeLengthKmers(r) <= shb.MedianCorrectedReadLengthFudge( ) )
     {    OK = True;    } 
     if ( !OK ) return False;

     // Find images of edges under involution.

     int rx1 = shb.Inv(x1), ry1 = shb.Inv(y1), rr = shb.Inv(r);
     int rx2 = shb.Inv(x2), ry2 = shb.Inv(y2);

     // Test for first special case.
     //
     // a1 --x1-->            --y1=ry2-->               --rx2--> c1
     //            v --r--> w             rv --rr--> rw
     // a2 --x2-->            --y2=ry1-->               --rx1--> c2
     //
     // This pulls apart to x1,r,y1,rr,rx2
     //                     x2,r,y2,rr,rx1,
     // which are rc to each other.

     Bool special1 = ( rx1 >= 0 && y1 == ry2 && r != rr
          && rx1 != x1 && rx1 != x2 && rx2 != x1 && rx2 != x2 );
     if (special1)
     {    if ( logc.PULL_APART_DEBUG ) cout << "\n";
          if ( logc.verb[ "PULL_APART" ] >= 1 || logc.PULL_APART_DEBUG )
          {    cout << "pulling part (special1) " << x1 << "," << r << "," << y1 
                    << " from " << x2 << "," << r << "," << y2 << endl;    }
          if ( logc.PULL_APART_DEBUG ) 
          {    PRINT8( x1, rx1, x2, rx2, r, rr, y1, y2 );
               PRINT4( v, w, id1, id2 );    }
          const basevector &X1 = shb.EdgeObject(x1), &RX1 = shb.EdgeObject(rx1); 
          const basevector &X2 = shb.EdgeObject(x2); 
          const basevector &RX2 = shb.EdgeObject(rx2);
          const basevector &R = shb.EdgeObject(r); 
          const basevector &RR = shb.EdgeObject(rr);
          const basevector &Y1 = shb.EdgeObject(y1), &Y2 = shb.EdgeObject(y2);
          basevector Z1 = shb.Cat( x1, r, y1, rr, rx2 );
          basevector Z2 = shb.Cat( x2, r, y2, rr, rx1 );
          int a1 = shb.To(v)[id1], a2 = shb.To(v)[id2];
          int rv = to_left[rr], rw = to_right[rr];
          shb.DeleteEdgesAtVertex(v), shb.DeleteEdgesAtVertex(w);
          shb.DeleteEdgesAtVertex(rw);
          int c1 = to_right[rx2], c2 = to_right[rx1];
          int z1 = shb.EdgeObjectCount( ), z2 = shb.EdgeObjectCount( ) + 1;
          shb.InvMutable( ).push_back( z2, z1 );
          shb.AddEdge(a1,c1,Z1), shb.AddEdge(a2,c2,Z2);
          to_left.push_back(a1,a2), to_right.push_back(c1,c2);
          vec< vec<int> > e(2);
          e[0].push_back( x1, r, y1, rr, rx2 );
          e[1].push_back( x2, r, y2, rr, rx1 );
          vec<int> f;
          f.push_back( z1, z2 );
          shb.TransformPaths( e, f, paths_index );
          return True;    }

     // Find images of edges under involution (different labeling).

     rx1 = shb.Inv(y1), ry1 = shb.Inv(x1), rr = shb.Inv(r),
     rx2 = shb.Inv(y2), ry2 = shb.Inv(x2);

     // Check for consistency with involution.

     Bool eq = False;
     if ( rx1 >= 0 )
     {    Bool left_eq = ( rx1 == x1 && rx2 == x2 ) || ( rx1 == x2 && rx2 == x1 );
          Bool right_eq = ( ry1 == y1 && ry2 == y2 ) || ( ry1 == y2 && ry2 == y1 );
          eq = left_eq && right_eq && rr == r;
          vec<int> s, t;
          if ( !eq )
          {    s.push_back( x1, x2, r, y1, y2 );
               t.push_back( rx1, rx2, rr, ry1, ry2 );
               UniqueSort(s), UniqueSort(t);    }
          if ( !eq && Meet(s,t) ) return False;    }

     // Two passes.

     vec< vec<int> > e;
     vec<int> f;
     if ( logc.PULL_APART_DEBUG ) 
     {    cout << "\n";
          PRINT6( x1, rx1, r, rr, y1, ry1 );
          PRINT4( x2, rx2, y2, ry2 );
          PRINT4( v, w, id1, id2 );
          PRINT4( to_left[x1], to_right[x1], to_left[x2], to_right[x2] );
          PRINT2( to_left[r], to_right[r] );
          PRINT4( to_left[y1], to_right[y1], to_left[y2], to_right[y2] );
          PRINT4( to_left[rx1], to_right[rx1], to_left[rx2], to_right[rx2] );
          PRINT2( to_left[rr], to_right[rr] );
          PRINT4( to_left[ry1], to_right[ry1], to_left[ry2], to_right[ry2] );    }
     int a1_1 = -1, a1_2 = -1, a2_1 = -1, a2_2 = -1;
     int b1_1 = -1, b1_2 = -1, b2_1 = -1, b2_2 = -1;
     basevector Z1_1, Z1_2, Z2_1, Z2_2;
     int v_1 = -1, v_2 = -1, w_1 = -1, w_2 = -1;
     int npasses = ( !eq ? 2 : 1 );
     for ( int xpass = 1; xpass <= npasses; xpass++ )
     {    if ( xpass == 2 && rx1 < 0 ) continue;
          int px1, px2, py1, py2, pr;
          if ( xpass == 1 )
          {    px1 = x1, px2 = x2, py1 = y1, py2 = y2, pr = r;    }
          else
          {    px1 = rx1, px2 = rx2, py1 = ry1, py2 = ry2;
               pr = rr;    }
          int v = to_left[pr], w = to_right[pr];
          if ( logc.PULL_APART_DEBUG ) PRINT3( xpass, v, w );
          int id1 = ( shb.EdgeObjectIndexByIndexTo(v,0) == px1 ? 0 : 1 );
          int id2 = 1 - id1;
          int id1b
               = ( shb.EdgeObjectIndexByIndexFrom(w,0) == py1 ? 0 : 1 );
          if ( xpass == 1 )
          {    a1_1 = shb.To(v)[id1], a2_1 = shb.To(v)[id2];    }
          else
          {    a1_2 = shb.To(v)[id1], a2_2 = shb.To(v)[id2];    }
          if ( xpass == 1 )
          {    b1_1 = shb.From(w)[id1b], b2_1 = shb.From(w)[1-id1b];    }
          else
          {    b1_2 = shb.From(w)[id1b], b2_2 = shb.From(w)[1-id1b];    }
          ( xpass == 1 ? Z1_1 : Z1_2 ) = shb.Cat( px1, pr, py1 );
          ( xpass == 1 ? Z2_1 : Z2_2 ) = shb.Cat( px2, pr, py2 );
          ( xpass == 1 ? v_1 : v_2 ) = v;
          ( xpass == 1 ? w_1 : w_2 ) = w;     }

     // Check for an illegal condition.

     if ( eq && ReverseComplement(Z1_1) != Z2_1 ) 
     {    if ( logc.PULL_APART_DEBUG ) cout << "illegal, bailing\n";
          return False;    }

     for ( int xpass = 1; xpass <= npasses; xpass++ )
     {    if ( xpass == 2 && rx1 < 0 ) continue;
          int px1, px2, py1, py2, pr;
          if ( xpass == 1 )
          {    px1 = x1, px2 = x2, py1 = y1, py2 = y2, pr = r;    }
          else
          {    px1 = rx1, px2 = rx2, py1 = ry1, py2 = ry2, pr = rr;    }
          vec<int> p11, p22;
          p11.push_back(px1,pr,py1), p22.push_back(px2,pr,py2);
          int v = ( xpass == 1 ? v_1 : v_2 );
          int w = ( xpass == 1 ? w_1 : w_2 );

          // Announce.

          if ( logc.verb[ "PULL_APART" ] >= 1 || logc.PULL_APART_DEBUG )
          {    cout << "pulling part " << px1 << "," << pr << "," << py1 
                    << " from " << px2 << "," << pr << "," << py2 << endl;    }

          // Edit.

          basevector Z1 = ( xpass == 1 ? Z1_1 : Z1_2 );
          basevector Z2 = ( xpass == 1 ? Z2_1 : Z2_2 );
          int a1 = ( xpass == 1 ? a1_1 : a1_2 );
          int a2 = ( xpass == 1 ? a2_1 : a2_2 );
          int b1 = ( xpass == 1 ? b1_1 : b1_2 );
          int b2 = ( xpass == 1 ? b2_1 : b2_2 );

          if ( logc.PULL_APART_DEBUG ) 
               cout << "deleting edges at " << v << " and " << w << endl;
          shb.DeleteEdgesAtVertex(v), shb.DeleteEdgesAtVertex(w);
          int z1 = shb.EdgeObjectCount( ), z2 = shb.EdgeObjectCount( ) + 1;
          if ( logc.PULL_APART_DEBUG ) 
          {    cout << "adding edge from " << a1 << " to " << b1 << endl;
               cout << "adding edge from " << a2 << " to " << b2 << endl;
               PRINT2( Z1.size( ), Z2.size( ) );    }
          shb.AddEdge(a1,b1,Z1), shb.AddEdge(a2,b2,Z2);
          if (eq) shb.InvMutable( ).push_back( z2, z1 );
          else if ( rx1 < 0 ) shb.InvMutable( ).push_back( -1, -1 );
          else if ( xpass == 1 ) 
          {    shb.InvMutable( ).push_back( z1 + 2, z2 + 2 );    }
          else if ( xpass == 2 ) 
          {    shb.InvMutable( ).push_back( z1 - 2, z2 - 2 );    }
          to_left.push_back(a1,a2), to_right.push_back(b1,b2);

          e.push_back( p11, p22 );
          f.push_back( z1, z2 );     }
     shb.TransformPaths( e, f, paths_index );
     return True;     }

void SupportedHyperBasevector::PullApart( const double min_weight_split,
     const long_logging& logc )
{    if (logc.STATUS_LOGGING) cout << Date( ) << ": pulling apart" << endl;
     while(1)
     {    double clock1 = WallClockTime( );
          vec<int> to_left, to_right;
          ToLeft(to_left), ToRight(to_right);

          vec< vec< pair<int,int> > > paths_index( EdgeObjectCount( ) );
          for ( int i = 0; i < Paths( ).isize( ); i++ )
	  for ( int j = 0; j < Path(i).isize( ); j++ ) 
	       paths_index[ Path(i,j) ].push(i,j);

          Bool changed = False;
          if ( logc.verb[ "PULL_APART" ] >= 1 )
               cout << Date( ) << ": start round of pulling apart" << endl;
          for ( int v = 0; v < N( ); v++ )
          {    if ( To(v).size( ) != 2 || From(v).size( ) != 1 ) continue;
               int w = From(v)[0];
               if ( To(w).size( ) != 1 || From(w).size( ) != 2 ) continue;
               for ( int pass = 1; pass <= 2; pass++ )
               {    if ( PullApartProcessVertex( *this, to_left, to_right, 
                         paths_index, v, w, pass, min_weight_split, logc ) )
                    {    changed = True;
                         if ( logc.PULL_APART_DEBUG ) 
                         {    TestValid( logc, False );
                              vec<int> to_leftx, to_rightx;
                              ToLeft(to_leftx), ToRight(to_rightx);
                              ForceAssertEq( to_left.size( ), to_leftx.size( ) );
                              ForceAssertEq( to_right.size( ), to_rightx.size( ) );
                              vec<Bool> used( EdgeObjectCount( ), False );
                              Used(used);
                              for ( int e = 0; e < EdgeObjectCount( ); e++ )
                              {    if ( !used[e] ) continue;
                                   if ( to_left[e] != to_leftx[e] )
                                   {    cout << "to_left wrong at edge " << e
                                             << endl;
                                        PRINT2( to_left[e], to_leftx[e] );
                                        TracebackThisProcess( );    }
                                   if ( to_right[e] != to_rightx[e] )
                                   {    cout << "to_right wrong at edge " << e
                                             << endl;
                                        PRINT2( to_right[e], to_rightx[e] );
                                        TracebackThisProcess( );    }    }    }
                         break;    }    }    }
          REPORT_TIME( clock1, "used in pull apart main loop" );
          if ( !changed ) break;
          double clock2 = WallClockTime( );

          vec<Bool> to_delete( NPaths( ), False );
          for ( int i = 0; i < NPaths( ); i++ )
               if ( Path(i).empty( ) ) to_delete[i] = True;
          EraseIf( PathsMutable( ), to_delete );
          EraseIf( WeightsFwMutable( ), to_delete );
          EraseIf( WeightsRcMutable( ), to_delete );

          to_delete.resize_and_set( NPairs( ), False );
          for ( int i = 0; i < NPairs( ); i++ )
          {    if ( PairLeft(i).empty( ) || PairRight(i).empty( ) ) 
                    to_delete[i] = True;    }
          EraseIf( PairsMutable( ), to_delete );
          EraseIf( PairDataMutable( ), to_delete );

          UniqueOrderPaths( );    
          REPORT_TIME( clock2, "used in pull apart tail" );
          RemoveDeadEdgeObjects( );    }
     TestValid(logc);
     DeleteReverseComplementComponents(logc);    }

void SupportedHyperBasevector::UnwindAssembly( 
     const long_logging_control& log_control, const long_logging& logc )
{
     double clock1 = WallClockTime();
     int verbosity = logc.verb["UNWIND"];
     if (verbosity >= 1)
	  cout << Date() << ": entering UnwindAssembly" << endl;
     if (verbosity >= 3) {
	  cout << "\nEdges in assembly:\n";
	  for (int e = 0; e < EdgeObjectCount(); e++) {
	       cout << e << " [" << EdgeLengthKmers(e) << " kmers]";
	       if (logc.UNWIND_TRUTH) {
		    vec<placementy> p = log_control.FindGenomicPlacements(
			      EdgeObject(e));
		    if (p.nonempty()) {
			 cout << "    [";
			 for (int k = 0; k < p.isize(); k++) {
			      cout << " " << p[k].g << "." << p[k].pos << "-"
					<< p[k].Pos << ":"
					<< (p[k].fw ? "fw" : "rc");
			 }
			 cout << " ]";
		    }
	       }
	       cout << "\n";
	  }
	  cout << "\nRead paths:\n";
	  for (int i = 0; i < Paths().isize(); i++) {
	       cout << "[" << i << ",weight=" << Weight(i) << "] "
			 << printSeq(Path(i), " ") << "\n";
	  }
     }

     // Define heuristics.

     const int max_paths = 10;
     const int max_len = 20;
     const fix64_6 keep_score(2001, 1000);
     const int win_ratio = 5;
     const double min_protect = 5;

     // Compute the multiplicity of each edge in hb.  Create an index to the paths.

     vec<fix64_6> mult(EdgeObjectCount(), 0);
     vec<vec<pair<int, int> > > paths_index(EdgeObjectCount());
     for (int i = 0; i < Paths().isize(); i++)
	  for (int j = 0; j < Path(i).isize(); j++) {
	       if (Path(i, j) < 0)
		    continue;
	       mult[Path(i, j)] += Weight(i);
	       paths_index[Path(i, j)].push(i, j);
	  }

     // Look forward from each edge.  
     // 
     // In this process, we reject a subpath a,b1,...,bn,c if it is not seen in the
     // reads at least keep_score (2) times and the total kmer count in the bi is 
     // <= the median read length, provided that for some c', we see a,b1,...,bn,c' 
     // at least win_ratio (5) times more often.
     //
     // This defines a vector recs, whose entries are pairs (v,w), where v is a
     // vec< vec<int> > and w is a vec<int>.  Here v consists of a sequence of
     // paths that represent replacements for some part of the assembly.  All
     // have the same start and the same end.  Between them lie a collection of
     // edges w.

     vec<int> to_left, to_right;
     ToLeft(to_left), ToRight(to_right);
     vec<pair<vec<vec<int> >, vec<int> > > recs;
     if (verbosity >= 1)
	  cout << Date() << ": start main loop" << endl;
     REPORT_TIME( clock1, "used in unwind setup");
     if (logc.STATUS_LOGGING) cout << Date( ) << ": begin unwind main" << endl;
#pragma omp parallel for schedule(dynamic, 1)
     for (int e = 0; e < EdgeObjectCount(); e++) {
	  double clock = WallClockTime();
	  ostringstream out;
	  double clock1 = WallClockTime();
	  if (verbosity >= 2)
	       out << "\n" << Date() << ": looking forward from edge " << e
			 << endl;

	  // We create ps, a list of paths starting at e.  As soon as a path in ps
	  // gets to have length max_len, we move it to ps_term.  As soon as we
	  // have more than max_paths paths, we stop extending.

	  vec<vec<int> > ps, ps_term;
	  vec<int> p;
	  p.push_back(e);
	  ps.push_back(p);
	  while (ps.nonempty() && ps.isize() + ps_term.isize() <= max_paths) {
	       int min_size = 1000000000, best = -1;
	       for (int j = 0; j < ps.isize(); j++) {
		    if (ps[j].isize() < min_size) {
			 best = j;
			 min_size = ps[j].size();
		    }
	       }
	       vec<int> p = ps[best];
	       ps.erase(ps.begin() + best);
	       if (p.isize() == max_len)
		    ps_term.push_back(p);
	       else {
		    vec<int> nexts = GetNextEdges(*this, Paths(), WeightsFw(),
			      WeightsRc(), paths_index,
			      MedianCorrectedReadLengthFudge(), p, keep_score,
			      win_ratio, verbosity, out);
		    if (nexts.empty())
			 ps_term.push_back(p);
		    else {
			 for (int j = 0; j < nexts.isize(); j++) {
			      vec<int> q(p);
			      q.push_back(nexts[j]);
			      ps.push_back(q);
			 }
		    }
	       }
	  }
	  ps.append(ps_term);
	  if (verbosity >= 2) {
	       out << "paths seen:\n";
	       for (int j = 0; j < ps.isize(); j++)
		    out << "- " << printSeq(ps[j], " ") << "\n";
	  }
	  REPORT_TIMEX( clock1, "used in unwind main 1", out);

          // Test to see if e appears twice.  Note that we should also check later 
          // to see if the terminating edge appears twice.

          Bool twice = False;
          for (int j = 0; j < ps.isize(); j++)
          {    int count = 0;
               for ( int l = 0; l < ps[j].isize( ); l++ )
                    if ( ps[j][l] == e ) count++;
               if ( count >= 2 ) twice = True;    }
          if (twice)
	  {    if (verbosity >= 2)
		    out << "nothing to do\n";
               continue;    }

	  // Now truncate the ps on the right if needed.  To define this truncation
	  // criterion, suppose that a given p is represented thusly:
	  //
	  //              v1 --p1--> v2 ... vn --pn--> vn+1
	  //
	  // and G is the graph obtained from hb by first removing p1 and pn,
	  // then forming the smallest complete subgraph that contains v2,...,vn.
	  // Now consider the edges in G that are not contained in p.  If the
	  // multiplicity of such an edge is >= min_protect (5), then the p is
	  // regarded as invalid.  However truncation of p on the right may fix
	  // the problem.
	  //
	  // Note: rewritten, a little different, see code.

	  Bool found_good = False;
	  for (int s = ps[0].isize(); s > 1; s--) {
	       double clocka = WallClockTime();
	       vec<vec<int> > psx(ps);
	       psx[0].resize(s);
	       if (psx[0].CountValue(psx[0].back()) != 1)
		    continue;
	       Bool bad = False;
	       for (int j = 1; j < psx.isize(); j++) {
		    int m = Position(psx[j], psx[0].back());
		    if (m < 0) {
			 bad = True;
			 break;
		    }
		    psx[j].resize(m + 1);
	       }
	       REPORT_TIMEX( clocka, "used in unwind main a", out);
	       if (bad)
		    continue;
	       double clockb = WallClockTime();

	       // Define truncated path set.

	       vec<int> v, qe;
	       if (verbosity >= 2) {
		    out << "testing truncated path set:\n";
		    for (int j = 0; j < psx.isize(); j++)
			 out << "- " << printSeq(psx[j], " ") << "\n";
	       }
	       for (int j = 0; j < psx.isize(); j++) {
		    const vec<int>& q = psx[j];
		    for (int j = 0; j < q.isize() - 1; j++) {
			 v.push_back(to_right[q[j]]);
			 if (j > 0)
			      qe.push_back(q[j]);
		    }
	       }
	       UniqueSort(qe);

	       // Define edges not to be used = first and last edges of psx[0].

	       int edge1 = psx[0].front(), edge2 = psx[0].back();
	       int vertex1 = to_left[ psx[0].front() ], vertex2 = to_right[ psx[0].back() ];

	       // Find dels = all edges that are connected to v, without passing
	       // through edge1 or edge2, and test for protected edges.  For
	       // efficiency, we test as we build dels.

	       vec<Bool> evil(EdgeObjectCount(), False);
	       for (int e = 0; e < EdgeObjectCount(); e++) {
		    if (!BinMember(qe, e) && mult[e] >= min_protect
			      && e != edge1 && e != edge2) {
			 evil[e] = True;
		    }
	       }
	       vec<int> dels0;
	       set<int> dels;
	       REPORT_TIMEX( clockb, "used in unwind main b", out);
	       double clockc = WallClockTime();
	       for (int j = 0; j < v.isize(); j++) {
		    int x = v[j];
		    for (int l = 0; l < To(x).isize(); l++) {
			 int e = EdgeObjectIndexByIndexTo(x, l);
			 if (e == edge1 || e == edge2)
			      continue;
			 if (evil[e]) {
			      bad = True;
			      break;
			 }
			 dels0.push_back(e);
			 dels.insert(e);
		    }
		    if (bad)
			 break;
		    for (int l = 0; l < From(x).isize(); l++) {
			 int e = EdgeObjectIndexByIndexFrom(x, l);
			 if (e == edge1 || e == edge2)
			      continue;
			 if (evil[e]) {
			      bad = True;
			      break;
			 }
			 dels0.push_back(e);
			 dels.insert(e);
		    }
		    if (bad)
			 break;
	       }
	       REPORT_TIMEX( clockc, "used in unwind main c", out);
	       if (bad)
		    continue;
	       double clockd = WallClockTime();
	       UniqueSort(dels0);
	       while (dels0.nonempty()) {
		    int e = dels0.back();
		    dels0.pop_back();
		    int v = to_left[e], w = to_right[e];
		    if ( v == vertex1 || v == vertex2 || w == vertex1 || w == vertex2 ) {
			 bad = True;
			 break;
		    }
		    for (int pass = 1; pass <= 2; pass++) {
			 int y = (pass == 1 ? v : w);
			 for (int j = 0; j < To(y).isize(); j++) {
			      int m = EdgeObjectIndexByIndexTo(y, j);
			      if (Member(dels, m))
				   continue;
			      if (m == edge1 || m == edge2)
				   continue;
			      if (evil[m]) {
				   bad = True;
				   break;
			      }
			      dels0.push_back(m);
			      dels.insert(m);
			 }
			 if (bad)
			      break;
			 for (int j = 0; j < From(y).isize(); j++) {
			      int m = EdgeObjectIndexByIndexFrom(y, j);
			      if (Member(dels, m))
				   continue;
			      if (m == edge1 || m == edge2)
				   continue;
			      if (evil[m]) {
				   bad = True;
				   break;
			      }
			      dels0.push_back(m);
			      dels.insert(m);
			 }
			 if (bad)
			      break;
		    }
		    if (bad)
			 break;
	       }
	       REPORT_TIMEX( clockd, "used in unwind main d", out);
	       if (bad)
		    continue;
	       double clocke = WallClockTime();

	       // Now test for compatibility with involution.

	       if (InvDef(edge1)) {
		    int redge1 = Inv(edge2), redge2 = Inv(edge1);
		    set<int> rdels;
		    for (set<int>::iterator i = dels.begin(); i != dels.end();
			      i++) {
			 rdels.insert(Inv(*i));
		    }
		    set<int> delsp(dels), rdelsp(rdels);
		    delsp.insert(edge1), delsp.insert(edge2);
		    rdelsp.insert(redge1), rdelsp.insert(redge2);
		    int over = 0;
		    for (set<int>::iterator i = delsp.begin(); i != delsp.end();
			      i++) {
			 if (Member(rdelsp, *i))
			      over++;
		    }
		    if (over > 0) {
			 if (edge1 != redge1 || edge2 != redge2
				   || dels != rdels)
			      continue;
		    }
	       }

	       // OK, everything is hunky-dory.

	       if (verbosity >= 2) {
		    out << "v = " << printSeq(v, " ") << "\n"
			      << "connected to ";
		    if (dels.size() <= 100)
			 out << printSeq(dels, " ") << "\n";
		    else
			 out << "[" << dels.size() << " edges]\n";
	       }
	       ps = psx;
	       found_good = True;
	       REPORT_TIMEX( clocke, "used in unwind main e", out);
	       break;
	  }
	  double clock3 = WallClockTime();

	  // Save recommendation.

	  UniqueSort(ps);

	  if (ps[0].solo() || !found_good) {
	       if (verbosity >= 2)
		    out << "nothing to do\n";
	  } else {
//	       vec<int> dels = Betweens(*this, to_right, ps[0]);
	       vec<int> dels = Betweens2( *this, to_left, to_right, ps[0] );

	       // we will get back no edges in the following case:
	       //
	       //         o | --> o --> | o
	       //
	       // as long as the central vertex connects to nothing else.
	       // Therefore there is nothing to do.
	       if ( dels.size() == 0 )
		    continue;

	       if (verbosity >= 1) {
		    out << "[" << recs.size() + 1 << "] recommend replacing "
			      << ps[0].front() << " {" << printSeq(dels, ",")
			      << "} " << ps[0].back() << " by";
		    for (int j = 0; j < ps.isize(); j++) {
			 if (j > 0)
			      out << " or";
			 out << " " << printSeq(ps[j], " ");
		    }
		    out << "\n";
	       }


	       if (!InvDef(ps[0].front())) {
#pragma omp critical
		    recs.push(ps, dels);
	       } else {
		    pair<vec<vec<int> >, vec<int> > recp = pair<vec<vec<int> >,
			      vec<int> >(ps, dels);
		    pair<vec<vec<int> >, vec<int> > recq;
		    for (int j = 0; j < recp.first.isize(); j++) {
			 vec<int> v;
			 for (int k = 0; k < recp.first[j].isize(); k++)
			      v.push_back(Inv(recp.first[j][k]));
			 v.ReverseMe();
			 recq.first.push_back(v);
		    }
		    Sort(recq.first);
		    for (int j = 0; j < recp.second.isize(); j++)
			 recq.second.push_back(Inv(recp.second[j]));
		    Sort(recq.second);
#pragma omp critical
		    {
			 recs.push_back(recp);
			 recs.push_back(recq);
		    }
	       }

	  }
	  if (verbosity >= 2)
	       out << "done, time used = " << TimeSince(clock) << endl;
	  REPORT_TIMEX( clock3, "used in unwind main 3", out);
	  if (verbosity >= 1 || logc.PRINT_TIME_USED) {
#pragma omp critical
	       {
		    cout << out.str();
	       }
	  }
     }

     // Intersect recommendations.

     double iclock = WallClockTime();
     if (logc.STATUS_LOGGING) 
          cout << Date( ) << ": intersecting recommendations" << endl;
     UniqueSort(recs);
     {
	  vec<vec<int> > content(recs.size());
	  vec<int> fronts, backs;
	  for (int j = 0; j < recs.isize(); j++) {
	       fronts.push_back(recs[j].first[0].front());
	       backs.push_back(recs[j].first[0].back());
	  }
	  for (int j = 0; j < recs.isize(); j++) {
	       content[j].push_back(fronts[j], backs[j]);
	       content[j].append(recs[j].second);
	       UniqueSort(content[j]);
	  }
	  SortSync(content, recs);
	  for (int i = 0; i < content.isize(); i++) {
	       int j = content.NextDiff(i);
	       vec<vec<vec<int> > > firsts;
	       for (int k = i; k < j; k++)
		    firsts.push_back(recs[k].first);
	       vec<vec<int> > com;
	       Intersection(firsts, com);
	       if (com.nonempty()) {
		    for (int k = i; k < j; k++)
			 recs[k].first = com;
	       }
	       i = j - 1;
	  }
     }
     UniqueSort(recs);
     REPORT_TIME( iclock, "used intersecting");

     // Process the recommendations.  Delete subsets, and if there's a conflict,
     // delete the smaller one.

     double tclock = WallClockTime();
     if (verbosity >= 1) cout << Date() << ": processing recommendations" << endl;
     UniqueSort(recs);
     vec<Bool> to_delete(recs.size(), False);
     vec<vec<int> > content(recs.size());
     vec<int> fronts, backs;
     for (int j = 0; j < recs.isize(); j++) {
	  fronts.push_back(recs[j].first[0].front());
	  backs.push_back(recs[j].first[0].back());
     }
     for (int j = 0; j < recs.isize(); j++) {
	  content[j].push_back(fronts[j], backs[j]);
	  content[j].append(recs[j].second);
	  UniqueSort(content[j]);
     }
     for (int i1 = 0; i1 < recs.isize(); i1++)
	  for (int i2 = 0; i2 < recs.isize(); i2++) {
	       if (to_delete[i1] || to_delete[i2])
		    continue;
	       if (i1 == i2)
		    continue;
	       if (fronts[i1] == fronts[i2] && backs[i1] == backs[i2]
			 && content[i1] == content[i2]) {
		    if (i2 > i1)
			 to_delete[i2] = True;
	       } else if (BinSubset(content[i2], content[i1]))
		    to_delete[i2] = True;
	       else {
		    vec<int> s1 = content[i1], s2 = content[i2];
		    if (Meet(s1, s2)) {
			 if (s1.size() > s2.size()
				   || (s1.size() == s2.size() && i1 > i2)) {
			      if (verbosity >= 1) {
				   cout << "conflict, deleting recommendation "
					     << i2 + 1 << "\n";
			      }
			      to_delete[i2] = True;
			 }
		    }
	       }
	  }

     for (int i1 = 0; i1 < recs.isize(); i1++)
	  for (int i2 = 0; i2 < recs.isize(); i2++) {
	       if (i1 == i2)
		    continue;
	       if (to_delete[i1] || to_delete[i2])
		    continue;
	       if (fronts[i1] != fronts[i2] || backs[i1] != backs[i2])
		    continue;
	       pair<vec<vec<int> >, vec<int> > recp1 = recs[i1], recp2 =
			 recs[i2];
	       if (recp1.second != recp2.second)
		    continue;

	       if (recp1.first.size() < recp2.first.size())
		    to_delete[i2] = True;
	       if (recp1.first.size() == recp2.first.size() && i1 < i2)
		    to_delete[i2] = True;

	       if (InvDef(fronts[i1])) {
		    pair<vec<vec<int> >, vec<int> > recq1, recq2;
		    for (int j = 0; j < recp1.first.isize(); j++) {
			 vec<int> v;
			 for (int k = 0; k < recp1.first[j].isize(); k++)
			      v.push_back(Inv(recp1.first[j][k]));
			 v.ReverseMe();
			 recq1.first.push_back(v);
		    }
		    Sort(recq1.first);
		    for (int j = 0; j < recp1.second.isize(); j++)
			 recq1.second.push_back(Inv(recp1.second[j]));
		    Sort(recq1.second);
		    for (int j = 0; j < recp2.first.isize(); j++) {
			 vec<int> v;
			 for (int k = 0; k < recp2.first[j].isize(); k++)
			      v.push_back(Inv(recp2.first[j][k]));
			 v.ReverseMe();
			 recq2.first.push_back(v);
		    }
		    Sort(recq2.first);
		    for (int j = 0; j < recp2.second.isize(); j++)
			 recq2.second.push_back(Inv(recp2.second[j]));
		    Sort(recq2.second);
		    int j1 = BinPosition(recs, recq1), j2 = BinPosition(recs,
			      recq2);
		    if (j2 >= 0 && to_delete[i2])
			 to_delete[j2] = True;
	       }
	  }

     EraseIf(recs, to_delete);
     REPORT_TIME( tclock, "used dealing with recs");

     // Find the boxes defined by each recommendation.

     double zclock = WallClockTime();		// ******* "unwind tail" starts here
     vec<triple<int, int, vec<int> > > box;
     vec<int> box_id(recs.size(), vec<int>::IDENTITY);
     for (int i = 0; i < recs.isize(); i++) {
	  box.push(recs[i].first[0].front(), recs[i].first[0].back(),
		    recs[i].second);
     }
     SortSync(box, box_id);

     // Define involution of recs.

     if ( verbosity >= 1 ) cout << ": defining involution of recs" << endl;
     vec<int> rinv(recs.size(), -1);
     for (int i = 0; i < box.isize(); i++) {
	  triple<int, int, vec<int> > b;
	  b.first = Inv(box[i].second), b.second = Inv(box[i].first);
	  for (int j = 0; j < box[i].third.isize(); j++)
	       b.third.push_back(Inv(box[i].third[j]));
	  Sort(b.third);
	  int p = BinPosition(box, b);
	  if (p >= 0)
	       rinv[i] = box_id[p];
     }

     // Use involution to improve recs.

     for (int i = 0; i < recs.isize(); i++) {
	  int r = rinv[i];
	  if (r < 0)
	       continue;
	  vec<vec<int> > x;
	  for (int j = 0; j < recs[i].first.isize(); j++) {
	       vec<int> v;
	       for (int k = 0; k < recs[i].first[j].isize(); k++)
		    v.push_back(Inv(recs[i].first[j][k]));
	       v.ReverseMe();
	       x.push_back(v);
	  }
	  Sort(x);
	  if (x != recs[r].first && BinSubset(x, recs[r].first))
	       recs[r].first = x;
     }

     // Make sure that what we're doing is consistent with the involution.

     vec<int> to_inv;
     for (int pass = 1; pass <= 2; pass++) {
	  to_inv.clear();
	  to_inv.resize(recs.size(), -1);
	  vec<Bool> to_delete2(recs.size(), False);
	  for (int i = 0; i < recs.isize(); i++) {
	       pair<vec<vec<int> >, vec<int> >& recp = recs[i];
	       if (!InvDef(recp.first[0].front()))
		    continue;
	       pair<vec<vec<int> >, vec<int> > recq;
	       for (int j = 0; j < recp.first.isize(); j++) {
		    vec<int> v;
		    for (int k = 0; k < recp.first[j].isize(); k++)
			 v.push_back(Inv(recp.first[j][k]));
		    v.ReverseMe();
		    recq.first.push_back(v);
	       }
	       Sort(recq.first);
	       for (int j = 0; j < recp.second.isize(); j++)
		    recq.second.push_back(Inv(recp.second[j]));
	       Sort(recq.second);
	       to_inv[i] = BinPosition(recs, recq);
	       if (pass == 1 && to_inv[i] < 0)
		    to_delete2[i] = True;
	  }
	  if (pass == 1)
	       EraseIf(recs, to_delete2);
     }

     // List surviving recommendations.

     if (verbosity >= 1) {
	  cout << "\ninvolution:\n";
	  for (int e = 0; e < EdgeObjectCount(); e++)
	       if (InvDef(e))
		    cout << e << " --> " << Inv(e) << "\n";
	  cout << "\nsurviving recommendations:\n";
	  for (int i = 0; i < recs.isize(); i++) {
	       Bool placed = False;
	       cout << "\n[" << i + 1 << "] replace "
			 << recs[i].first[0].front() << " {"
			 << printSeq(recs[i].second, ",") << "} "
			 << recs[i].first[0].back() << " by\n  ";
	       for (int j = 0; j < recs[i].first.isize(); j++) {
		    if (j > 0)
			 cout << "\nor";
		    cout << " " << printSeq(recs[i].first[j], " ");
		    if (logc.UNWIND_TRUTH) {
			 basevector b;
			 for (int l = 0; l < recs[i].first[j].isize(); l++) {
			      if (l == 0)
				   b = EdgeObject(recs[i].first[j][l]);
			      else
				   b = TrimCat(K(), b,
					EdgeObject(recs[i].first[j][l]));
			 }
			 vec<placementy> p = log_control.FindGenomicPlacements(
				   b);
			 cout << "    [";
			 if (p.nonempty()) {
			      placed = True;
			      for (int k = 0; k < p.isize(); k++) {
				   cout << " " << p[k].g << "." << p[k].pos
					     << "-" << p[k].Pos << ":"
					     << (p[k].fw ? "fw" : "rc");
			      }
			 } else {
			      cout << " N/A ";
			 }
			 cout << " ]";
		    }
	       }
	       cout << "\n";
	       if (logc.UNWIND_TRUTH && !placed) {
		    vec<int> s;
		    s.push_back(recs[i].first[0].front());
		    s.append(recs[i].second);
		    s.push_back(recs[i].first[0].back());
		    vec<vec<placementy> > places(s.size());
		    for (int j = 0; j < s.isize(); j++) {
			 places[j] = log_control.FindGenomicPlacements(
				   EdgeObject(s[j]));
		    }
		    vec<seq_place> splacesp, splacesf;
		    for (int j = 0; j < places[0].isize(); j++) {
			 vec<int> x;
			 x.push_back(0);
			 const placementy& p = places[0][j];
			 splacesp.push(x, p.g, p.fw, p.pos, p.Pos);
		    }
		    while (splacesp.nonempty()) {
			 seq_place p = splacesp.back();
			 splacesp.pop_back();
			 int u = s[p.x.back()];
			 int ur = to_right[u];
			 for (int r = 0; r < From(ur).isize(); r++) {
			      int v = EdgeObjectIndexByIndexFrom(ur, r);
			      int l = Position(s, v);
			      if (l < 0)
				   continue;
			      for (int m = 0; m < places[l].isize(); m++) {
				   const placementy& q = places[l][m];
				   if (q.g != p.g || q.fw != p.fw)
					continue;
				   if ((p.fw && p.stop - q.pos == K() - 1)
					     || (!p.fw
						       && q.Pos - p.start
								 == K() - 1)) {
					vec<int> x(p.x);
					x.push_back(l);
					int start = (p.fw ? p.start : q.pos);
					int stop = (p.fw ? q.Pos : p.stop);
					seq_place pn(x, p.g, p.fw, start, stop);
					if (l == s.isize() - 1)
					     splacesf.push_back(pn);
					else
					     splacesp.push_back(pn);
				   }
			      }
			 }
		    }
		    for (int j = 0; j < splacesf.isize(); j++) {
			 const seq_place& p = splacesf[j];
			 cout << "Note:";
			 for (int l = 0; l < p.x.isize(); l++)
			      cout << " " << s[p.x[l]];
			 cout << " placed at " << p.g << "." << p.start << "-"
				   << p.stop << ":" << (p.fw ? "fw" : "rc")
				   << ".\n";
		    }
	       }
	  }
	  cout << "\n";
     }

     // Now instantiate the recommendations.

     if (logc.STATUS_LOGGING)
     {    cout << Date( ) << ": instantiating " << recs.size( ) 
               << " recommendations 1" << endl;    }
     vec< digraphE<int> > Ds( recs.size( ) );
     vec< vec<int> > ss( recs.size( ) );
     vec<int> t1s( recs.size( ) ), t2s( recs.size( ) );
     vec<int> s1s( recs.size( ) ), s2s( recs.size( ) );
     vec< vec<int> > edge_maps( recs.size( ) );
     vec<Bool> ok( recs.size( ), False );
     #pragma omp parallel for
     for (int i = 0; i < recs.isize(); i++) 
     {    digraphE<int>& D = Ds[i];
          vec<int>& s = ss[i];
          int &s1 = s1s[i], &s2 = s2s[i], &t1 = t1s[i], &t2 = t2s[i];
          vec<int>& edge_map = edge_maps[i];
          int ri = to_inv[i];
	  if ( ri > i ) continue;
	  vec< vec<int> > w = recs[i].first;
	  Sort(w);
	  DigraphFromWords( w, D );
	  s1 = recs[i].first[0].front(), s2 = recs[i].first[0].back();
	  s.push_back(s1, s2);
	  s.append(recs[i].second);

	  // Find the unique source and sink of D.

	  vec<int> sources, sinks;
	  D.Sources(sources), D.Sinks(sinks);
	  ForceAssert( sources.solo( ) && sinks.solo( ));
	  t1 = sources[0], t2 = sinks[0];

	  // Set up for symmetric-identical case.  This can fail.

	  edge_map.resize( D.EdgeObjectCount(), -1 );
	  if (ri == i) 
          {    // We need to find an involution of D that is compatible with
	       // the given involution.

	       vec<int> to_left, to_right;
	       D.ToLeft(to_left), D.ToRight(to_right);
	       vec<int> vert_map(D.N(), -1);
	       while (1) 
               {    Bool progress = False;
		    for (int j = 0; j < D.EdgeObjectCount(); j++) 
                    {    if (edge_map[j] >= 0) continue;
			 int e = D.EdgeObject(j);
			 int v = to_left[j], w = to_right[j];
			 vec<int> places;
			 for (int l = 0; l < D.EdgeObjectCount(); l++) 
                         {    if (D.EdgeObject(l) == Inv(e)) 
                              {    if (edge_map[l] >= 0 && edge_map[l] != j)
					continue;
				   int rv = to_left[l], rw = to_right[l];
				   if (vert_map[v] >= 0 && vert_map[v] != rw)
					continue;
				   if (vert_map[w] >= 0 && vert_map[w] != rv)
					continue;
				   if (vert_map[rv] >= 0 && vert_map[rv] != w)
					continue;
				   if (vert_map[rw] >= 0 && vert_map[rw] != v)
					continue;
				   places.push_back(l);    }    }
			 if (places.solo()) 
                         {    progress = True;
			      int l = places[0];
			      edge_map[j] = l;
			      edge_map[l] = j;
			      int rv = to_left[l], rw = to_right[l];
			      vert_map[v] = rw, vert_map[w] = rv;
			      vert_map[rv] = w, vert_map[rw] = v;    }    }
		    if (!progress) break;    }

	       // Check for success.  Give up otherwise.

	       Bool fail = False;
	       for (int j = 0; j < D.EdgeObjectCount(); j++)
	       {    if (edge_map[j] < 0) fail = True;    }
	       if (fail) continue;    }
          ok[i] = True;    }

     if (logc.STATUS_LOGGING)
          cout << Date( ) << ": instantiating recommendations 2" << endl;
     for (int i = 0; i < recs.isize(); i++) 
     {    if ( !ok[i] ) continue;
          const digraphE<int>& D = Ds[i];
          const vec<int>& s = ss[i];
          const int &s1 = s1s[i], &s2 = s2s[i], &t1 = t1s[i], &t2 = t2s[i];
          const vec<int>& edge_map = edge_maps[i];
          int ri = to_inv[i];

	  // Convert D into HyperBasevector Db.

	  vec<basevector> edgesb;
	  for (int j = 0; j < D.EdgeObjectCount(); j++)
	       edgesb.push_back( EdgeObject(D.EdgeObject(j)) );
	  HyperBasevector Db( K(), D, edgesb );

	  // Now edit the assembly.

	  DeleteEdges(s);

	  vec<HyperBasevector> hs;
	  hs.push_back(Db);
	  int nverts = N();
	  int nedges0 = EdgeObjectCount();
	  AppendDisjointUnionOf(hs);
	  TransferEdges(nverts + t1, to_left[s1]);
	  TransferEdges(nverts + t2, to_right[s2]);

	  ToLeft(to_left), ToRight(to_right); /* probably not necessary */
	  InvMutable().resize(EdgeObjectCount(), -1);

	  // Do the symmetric-identical case.

	  if (ri == i) 
          {    for (int j = 0; j < D.EdgeObjectCount(); j++)
		    InvMutable(nedges0 + j) = nedges0 + edge_map[j];
	       continue;    }

	  // Do the symmetric-nonidentical case.

	  if (ri < 0) continue;
	  int rs1 = Inv(s2), rs2 = Inv(s1);
	  vec<int> rs;
	  for (int j = 0; j < s.isize(); j++)
	       rs.push_back(Inv(s[j]));

	  DeleteEdges(rs);
	  vec<HyperBasevector> rhs;
	  HyperBasevector rDb(Db);
	  rDb.Reverse();
	  rhs.push_back(rDb);
	  nverts = N();
	  AppendDisjointUnionOf(rhs);
	  TransferEdges( nverts + t2, to_left[rs1]);
	  TransferEdges( nverts + t1, to_right[rs2]);

	  ToLeft(to_left), ToRight(to_right); /* probably not necessary */
	  InvMutable().resize(EdgeObjectCount(), -1);
	  for (int j = 0; j < D.EdgeObjectCount(); j++) {
	       InvMutable(nedges0 + j) = nedges0 + D.EdgeObjectCount() + j;
	       InvMutable(nedges0 + D.EdgeObjectCount() + j) = nedges0 + j;
	  }
     }


     // Clean up.

     if (logc.STATUS_LOGGING) cout << Date( ) << ": unwind cleaning up" << endl;
     vec<triple<int, int, int> > merges;
     RemoveUnneededVertices0(merges);
     RemoveEdgelessVertices();
     REPORT_TIME( zclock, "used in unwind tail");
     RemoveDeadEdgeObjects0();
     if (verbosity >= 1) cout << Date() << ": leaving UnwindAssembly" << endl;    }
