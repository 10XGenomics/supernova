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
#include "lookup/LookAlign.h"
#include "pairwise_aligners/ClusterAligner.h"
#include "paths/long/EvalByReads.h"
#include "paths/long/LargeKDispatcher.h"
#include "paths/long/MakeKmerStuff.h"
#include "paths/long/RefTrace.h"
#include "paths/long/SupportedHyperBasevector.h"
#include "random/Bernoulli.h"
#include "reporting/PerfStat.h"
#include "util/NullOStream.h"
#include "paths/long/ReadOriginTracker.h"

namespace { // open anonymous namespace

void FixPath( vec<int>& p, const vec< triple<int,int,int> >& merges,
     const vec< vec<int> >& merges_index, const HyperBasevector& hb_orig, 
     int& left_add, int& right_add )
{    left_add = 0, right_add = 0;
     while(1)
     {    Bool changed = False;
          vec<int> mids;
          for ( int l = 0; l < p.isize( ); l++ )
               if ( p[l] >= 0 ) mids.append( merges_index[ p[l] ] );
          UniqueSort(mids);
          for ( int mj = 0; mj < mids.isize( ); mj++ )
          {    int j = mids[mj];
               int e1 = merges[j].first, e2 = merges[j].second;
               ForceAssert( e1 != e2 );
               int enew = merges[j].third;
               for ( int l = 0; l < p.isize( ); l++ )
               {    if ( l < p.isize( ) - 1 && p[l] == e1 && p[l+1] == e2 )
                    {    p[l] = enew;
                         for ( int m = l+2; m < p.isize( ); m++ )
                              p[m-1] = p[m];
                         p.pop_back( );
                         changed = True;    }
                    else if ( l == p.isize( ) - 1 && p[l] == e1 )
                    {    right_add += hb_orig.EdgeLengthKmers(e2);
                         p[l] = enew;
                         changed = True;    }
                    else if ( l == 0 && p[l] == e2 )
                    {    left_add += hb_orig.EdgeLengthKmers(e1);
                         p[l] = enew;
                         changed = True;    }
                    if (changed) break;    }
               if (changed) break;    }
          if ( !changed ) break;    }    }

} // close anonymous namespace

void SupportedHyperBasevector::RemoveUnneededVertices0( 
     vec< triple<int,int,int> >& merges )
{    vec<int> to_left, to_right;
     ToLeft(to_left), ToRight(to_right);

     for ( int i = 0; i < N( ); i++ )
     {    if ( From(i).size( ) == 1 && To(i).size( ) == 1 && From(i)[0] != i )
          {    int e1 = EdgeObjectIndexByIndexTo( i, 0 );
               int e2 = EdgeObjectIndexByIndexFrom( i, 0 );
               basevector p = Cat( e1, e2 );
               int enew = EdgeObjectCount( );
               int re1 = Inv(e1), re2 = Inv(e2);
               ForceAssert( ( re1 < 0 && re2 < 0 ) || ( re1 >= 0 && re2 >= 0 ) );
               int v = To(i)[0], w = From(i)[0];
               Bool loop = ( v == w && From(v).solo( ) && To(v).solo( ) );
               // v --e1--> i --e2--> w
               merges.push( e1, e2, enew );
               JoinEdges( i, p );
               to_left.push_back(v), to_right.push_back(w);
               if ( re1 < 0 ) InvMutable( ).push_back(-1);
               else if ( re1 == e2 && re2 == e1 )
               {    // * --e1=re2--> * --e2=re1--> *
                    InvMutable( ).push_back(enew);    }
               else if ( re2 == e2 && re1 != e1 )
               {    // * --e1--> * --e2=re2--> * --re1--> *
                    int enew2 = EdgeObjectCount( );
                    merges.push( enew, re1, enew2 );
                    basevector p2 = TrimCat( K( ), p, EdgeObject(re1) );
                    to_left.push_back(v), to_right.push_back( to_right[re1] );
                    JoinEdges( w, p2 );
                    InvMutable( ).push_back( -1, enew2 );    }
               else if ( re1 == e1 && re2 != e2 )
               {    // * --re2--> * --e1=re1--> * --e2--> *
                    int enew2 = EdgeObjectCount( );
                    merges.push( re2, enew, enew2 );
                    basevector p2 = TrimCat( K( ), EdgeObject(re2), p );
                    to_left.push_back( to_left[re2] ), to_right.push_back(w);
                    JoinEdges( v, p2 );
                    InvMutable( ).push_back( -1, enew2 );    }
               else if ( re1 == e1 && re2 == e2 )
               {    if (loop) InvMutable( ).push_back(-1);
                    else
                    {    // not sure if this can happen
                         ForceAssert( 0 == 1 );    }    }
               else
               {    // e1, e2, re1, re2 all different
                    int renew = EdgeObjectCount( );
                    basevector rp = Cat( re2, re1 );
                    merges.push( re2, re1, renew );
                    int ri = to_right[re2];
                    JoinEdges( ri, rp );
                    int rv = to_left[re2], rw = to_right[re1];
                    to_left.push_back(rv), to_right.push_back(rw);
                    InvMutable( ).push_back(renew, enew);    }    }    }    }

void SupportedHyperBasevector::RemoveUnneededVertices( )
{    vec< triple<int,int,int> > merges;
     RemoveUnneededVertices0(merges);
     vec< vec<int> > merges_index( EdgeObjectCount( ) );
     for ( int i = 0; i < merges.isize( ); i++ )
     {    merges_index[ merges[i].first ].push_back(i);
          merges_index[ merges[i].second ].push_back(i);    }
     #pragma omp parallel for
     for ( int i = 0; i < NPaths( ); i++ )
     {    int left_add, right_add;
          FixPath( PathMutable(i), merges, merges_index, *this,
               left_add, right_add );    }
     #pragma omp parallel for
     for ( int i = 0; i < NPairs( ); i++ )
     {    int left_add1, right_add1, left_add2, right_add2;
          FixPath( PairLeftMutable(i), merges, merges_index, *this,
               left_add1, right_add1 );
          FixPath( PairRightMutable(i), merges, merges_index, *this,
               left_add2, right_add2 );
          AddTrim( i, right_add1 + left_add2 );    }
     UniqueOrderPaths( );
     RemoveEdgelessVertices( );    }

void SupportedHyperBasevector::DeleteUnusedPaths( )
{    vec<Bool> used, to_delete( NPaths( ), False );
     Used(used);
     for ( int i = 0; i < NPaths( ); i++ )
     {    for ( int j = 0; j < Path(i).isize( ); j++ )
               if ( Path(i,j) >= 0 && !used[ Path(i,j) ] ) to_delete[i] = True;    }
     EraseIf( PathsMutable( ), to_delete );
     EraseIf( WeightsFwMutable( ), to_delete );    
     EraseIf( WeightsRcMutable( ), to_delete );    
     // The following if is temporary - until origins fully implemented.
     if ( to_delete.size( ) == WeightsFwOrigin( ).size( ) )
     {    EraseIf( WeightsFwOriginMutable( ), to_delete );    
          EraseIf( WeightsRcOriginMutable( ), to_delete );    }
     to_delete.resize_and_set( NPairs( ), False );
     for ( int i = 0; i < NPairs( ); i++ )
     {    vec<int> &p1 = PairLeftMutable(i), &p2 = PairRightMutable(i);
          for ( int pass = 1; pass <= 2; pass++ )
          {    vec<int>& p = ( pass == 1 ? p1 : p2 );
               for ( int j = 0; j < p.isize( ); j++ )
                    if ( p[j] >= 0 && !used[ p[j] ] ) to_delete[i] = True;    }    }
     EraseIf( PairsMutable( ), to_delete );
     EraseIf( PairDataMutable( ), to_delete );    }

void SupportedHyperBasevector::RemoveDeadEdgeObjects0( )
{    vec<Bool> used;
     Used(used);
     vec<int> to_new_id( used.size( ), -1 );
     {    int count = 0;
          for ( int i = 0; i < used.isize( ); i++ )
               if ( used[i] ) to_new_id[i] = count++;    }
     vec<int> inv2;
     for ( int i = 0; i < EdgeObjectCount( ); i++ )
     {    if ( !used[i] ) continue;
          if ( !InvDef(i) ) inv2.push_back(-1);
          else inv2.push_back( to_new_id[ Inv(i) ] );    }
     InvMutable( ) = inv2;
     HyperBasevector::RemoveDeadEdgeObjects( );    }

void SupportedHyperBasevector::RemoveDeadEdgeObjects( )
{    vec<Bool> used;
     Used(used);
     vec<int> to_new_id( used.size( ), -1 );
     {    int count = 0;
          for ( int i = 0; i < used.isize( ); i++ )
               if ( used[i] ) to_new_id[i] = count++;    }
     vec<int> inv2;
     for ( int i = 0; i < EdgeObjectCount( ); i++ )
     {    if ( !used[i] ) continue;
          if ( !InvDef(i) ) inv2.push_back(-1);
          else inv2.push_back( to_new_id[ Inv(i) ] );    }
     InvMutable( ) = inv2;

     vec<Bool> to_delete( NPaths( ), False );
     for ( int i = 0; i < NPaths( ); i++ )
     {    vec<int>& p = PathMutable(i);
          for ( int j = 0; j < Path(i).isize( ); j++ )
          {    if ( Path(i,j) >= 0 )
               {    int n = to_new_id[ Path(i,j) ];
                    if ( n < 0 ) to_delete[i] = True;
                    else PathMutable(i)[j] = n;    }    }    }
     EraseIf( PathsMutable( ), to_delete );
     // Temporary if.
     if ( WeightsFwOrigin( ).size( ) == WeightsFw( ).size( ) )
     {    EraseIf( WeightsFwOriginMutable( ), to_delete );
          EraseIf( WeightsRcOriginMutable( ), to_delete );    }
     EraseIf( WeightsFwMutable( ), to_delete );
     EraseIf( WeightsRcMutable( ), to_delete );

     to_delete.resize_and_set( NPairs( ), False );
     for ( int i = 0; i < NPairs( ); i++ )
     {    vec<int> &p1 = PairLeftMutable(i), &p2 = PairRightMutable(i);
          for ( int pass = 1; pass <= 2; pass++ )
          {    vec<int>& p = ( pass == 1 ? p1 : p2 );
               for ( int j = 0; j < p.isize( ); j++ )
               {    if ( p[j] >= 0 )
                    {    int n = to_new_id[ p[j] ];
                         if ( n < 0 ) to_delete[i] = True;
                         else p[j] = n;    }    }    }    }
     EraseIf( PairsMutable( ), to_delete );
     EraseIf( PairDataMutable( ), to_delete );

     HyperBasevector::RemoveDeadEdgeObjects( );    }

void SupportedHyperBasevector::UniqueOrderPaths( )
{    
     // Temporary if.
     if ( WeightsFwOrigin( ).size( ) == WeightsFw( ).size( ) )
     {    SortSync( PathsMutable( ), WeightsFwMutable( ), WeightsRcMutable( ),
               WeightsFwOriginMutable( ), WeightsRcOriginMutable( ) );    }
     else SortSync( PathsMutable( ), WeightsFwMutable( ), WeightsRcMutable( ) );
     vec<Bool> to_delete( NPaths( ), False );
     for ( int i = 0; i < NPaths( ); i++ )
     {    int j = Paths( ).NextDiff(i);
          fix64_6 cfw = 0.0, crc = 0.0;
          for ( int k = i; k < j; k++ )
               cfw += WeightFw(k);
          WeightFwMutable(i) = cfw;
          for ( int k = i; k < j; k++ )
               crc += WeightRc(k);
          WeightRcMutable(i) = crc;
          // Temporary if.
          if ( WeightsFwOrigin( ).size( ) == WeightsFw( ).size( ) )
          {    for ( int k = i+1; k < j; k++ )
               {    WeightFwOriginMutable(i).append(
                         WeightFwOriginMutable(k) );
                    WeightRcOriginMutable(i).append(
                         WeightRcOriginMutable(k) );    }
               Sort( WeightFwOriginMutable(i) );
               Sort( WeightRcOriginMutable(i) );    }
          for ( int k = i+1; k < j; k++ )
               to_delete[k] = True;
          if ( cfw + crc == 0 ) to_delete[i] = True;
          i = j - 1;   }
     EraseIf( PathsMutable( ), to_delete );
     // Temporary if.
     if ( WeightsFwOrigin( ).size( ) == WeightsFw( ).size( ) )
     {    EraseIf( WeightsFwOriginMutable( ), to_delete );
          EraseIf( WeightsRcOriginMutable( ), to_delete );    }
     EraseIf( WeightsFwMutable( ), to_delete );
     EraseIf( WeightsRcMutable( ), to_delete );

     // Now do pairs.

     SortSync( PairsMutable( ), PairDataMutable( ) );
     to_delete.resize_and_set( NPairs( ), False );
     for ( int i = 0; i < NPairs( ); i++ )
     {    vec<int> &p1 = PairLeftMutable(i), &p2 = PairRightMutable(i);
          int j = Pairs( ).NextDiff(i);
          vec<pair_point> x;
          for ( int k = i; k < j; k++ )
               x.append( PairData(k) );
          Sort(x);
          PairDataMutable(i) = x;
          for ( int k = i+1; k < j; k++ )
               to_delete[k] = True;
          if ( PairData(i).empty( ) ) to_delete[i] = True;
          i = j - 1;   }
     EraseIf( PairsMutable( ), to_delete );
     EraseIf( PairDataMutable( ), to_delete );    }

void SupportedHyperBasevector::KillWeakExits2( const long_heuristics& heur,
     const long_logging& logc )
{    
     // Assign weights to vertices.  Every path containing the vertex internally
     // contributes to the vertex's weight.

     double clock = WallClockTime( );
     if (logc.STATUS_LOGGING) cout << Date( ) << ": killing weak exits 2" << endl;
     
     // Run in both orientations.

     const fix64_6 weak = 2;
     const int strong_mult = 10;
     const fix64_6 min_save = 10;
     vec<int> dels;
     for ( int pass = 1; pass <= 2; pass++ )
     {    Reverse( );

          // Assign entering and exiting weights to edges.  Each path containing
          // the edge and something after contributes to the entering weight,
          // whereas each path containing the edge and something before contributes
          // to the exiting weight.

          vec<fix64_6> enter_wt( EdgeObjectCount( ), 0 );
          vec<fix64_6> exit_wt( EdgeObjectCount( ), 0 );
          vec<fix64_6> wt( EdgeObjectCount( ), 0 );
          for ( int i = 0; i < NPaths( ); i++ )
          {    const vec<int>& p = Path(i);
               for ( int j = 0; j < p.isize( ); j++ )
               {    wt[ p[j] ] += Weight(i);
                    if ( j > 0 ) exit_wt[ p[j] ] += Weight(i);
                    if ( j < p.isize( ) - 1 ) 
                         enter_wt[ p[j] ] += Weight(i);    }    }

          // Explore each edge as an entering point.

          for ( int v = 0; v < N( ); v++ )
          for ( int i = 0; i < To(v).isize( ); i++ )
          {    int e = EdgeObjectIndexByIndexTo( v, i );

               // Find a neighborhood that starts from v.  This neighborhood 
               // consists of all vertices within distance max_dist of v, but if 
               // there are more than max_see such vertices, the computation is 
               // aborted.
     
               const int max_dist = 4;
               const int max_see = 20;
               vec< pair<int,int> > see; // (vertex, distance)
               see.push( v, 0 );
               Bool fail = False;
               while(1)
               {    Bool progress = False;
                    for ( int s = 0; s < see.isize( ); s++ )
                    {    int w = see[s].first, d = see[s].second;
                         if ( d == max_dist ) continue;
                         for ( int j = 0; j < From(w).isize( ); j++ )
                         {    int x = From(w)[j];
                              Bool found = False;
                              for ( int t = 0; t < see.isize( ); t++ )
                              {    if ( see[t].first == x )
                                   {    found = True;
                                        if ( d + 1 < see[t].second )
                                        {    see[t].second = d + 1;
                                             progress = True;    }
                                        break;    }    }
                              if ( !found )
                              {    see.push( x, d + 1 );
                                   progress = True;    
                                   if ( see.isize( ) > max_see )
                                   {    fail = True;
                                        break;    }    }    }
                         if (fail) break;    }
                    if ( !progress || fail ) break;    }
               if (fail) continue;

               // Get the vertices in 'see'.

               vec<int> vs;
               for ( int s = 0; s < see.isize( ); s++ )
                    vs.push_back( see[s].first );
               Sort(vs);

               // Look for a weak exiting edge v1 --e1--> and a strong exiting edge
               // v2 --e2-->.  We also require that v is strong relative to v1.

               for ( int j1 = 0; j1 < vs.isize( ); j1++ )
               for ( int i1 = 0; i1 < From( vs[j1] ).isize( ); i1++ )
               {    int v1 = vs[j1];
                    int e1 = EdgeObjectIndexByIndexFrom( v1, i1 );
                    if ( exit_wt[e1] > weak ) continue;
                    if ( v1 == v ) continue; // not sure if needed
                    if ( enter_wt[e] < strong_mult * exit_wt[e1] ) continue;
                    if ( enter_wt[e] < strong_mult ) continue;

                    for ( int j2 = 0; j2 < vs.isize( ); j2++ )
                    for ( int i2 = 0; i2 < From( vs[j2] ).isize( ); i2++ )
                    {    int v2 = vs[j2];
                         int e2 = EdgeObjectIndexByIndexFrom( v2, i2 );
                         if ( exit_wt[e2] < strong_mult * exit_wt[e1] ) continue;
                         if ( exit_wt[e2] < strong_mult ) continue;
                         if ( v2 == v || v2 == v1 ) continue; // not sure if needed

                         // Now we have candidates v1, v2.

                         if ( logc.verb[ "KILL_WEAK_EXITS" ] >= 2 )
                         {    cout << "\nKillWeakExits2 "
                                   << ( pass == 1 ? "rc" : "fw" ) << ", looking at "
                                   << "e = " << e << "[enter_wt=" << enter_wt[e] 
                                   << "], e1 = " << e1 << "[exit_wt=" << exit_wt[e1]
                                   << "]" << ", e2 = " << e2 
                                   << "[exit_wt=" << exit_wt[e2] << "]" << endl;    }
     
                         // Determine if e --> ... --> {e1,e2} is bounding.
     
                         vec<Bool> vis( vs.size( ), False ); 
                         vec<Bool> checked( vs.size( ), False );
                         for ( int l = 0; l < vs.isize( ); l++ )
                         {    if ( vs[l] == v ) vis[l] = True;
                              if ( vs[l] == v1 || vs[l] == v2 ) vis[l] = True;    }
                         Bool fail = False;
                         while(1)
                         {    Bool progress = False;
                              for ( int j = 0; j < vis.isize( ); j++ )
                              {    if ( !vis[j] || checked[j] ) continue;
                                   checked[j] = True;
                                   int w = vs[j];
                                   for ( int l = 0; l < From(w).isize( ); l++ )
                                   {    if ( w == v1 && l == i1 ) continue;
                                        if ( w == v2 && l == i2 ) continue;
                                        int x = From(w)[l];
                                        int p = BinPosition( vs, x );
                                        if ( p < 0 )
                                        {    fail = True;
                                             break;    }
                                        if ( !vis[p] ) progress = True;
                                        vis[p] = True;    }
                                   if (fail) break;    }
                              if ( !progress ) break;    }
                         if (fail) continue;
                         for ( int j = 0; j < vis.isize( ); j++ )
                         {    if ( !vis[j] ) continue;
                              int w = vs[j];
                              for ( int l = 0; l < To(w).isize( ); l++ )
                              {    if ( w == v && l == i ) continue;
                                   int p = BinPosition( vs, To(w)[l] );
                                   if ( p < 0 || !vis[p] )
                                   {    fail = True;
                                        break;    }    }
                              if (fail) break;    }
                         if (fail) continue;
                         int p = BinPosition( vs, To(v)[i] );
                         if ( p >= 0 && vis[p] ) continue;
                         int p1 = BinPosition( vs, From(v1)[i1] );
                         if ( p1 >= 0 && vis[p1] ) continue;
                         int p2 = BinPosition( vs, From(v2)[i2] );
                         if ( p2 >= 0 && vis[p2] ) continue;
                         if ( logc.verb[ "KILL_WEAK_EXITS" ] >= 2 )
                              cout << "bounding\n";

                         // Don't delete high-weight edges.

                         if ( wt[e1] >= min_save ) continue;

                         // Edge e1 is bad.  Mark it for deletion.
     
                         vec<int> delsx;
                         delsx.push_back(e1);
                         dels.append(delsx);
                         /*
                         for ( int j = 0; j < delsx.isize( ); j++ )
                         {    if ( InvDef( delsx[j] ) ) 
                                   dels.push_back( Inv( delsx[j] ) );    }
                         */

                         if ( logc.verb[ "KILL_WEAK_EXITS" ] >= 1 )
                         {    cout << "\nKillWeakExits2, success, "
                                   << "e = " << e << "[enter_wt=" << enter_wt[e] 
                                   << "], e1 = " << e1 << "[exit_wt=" << exit_wt[e1]
                                   << "]" << ", e2 = " << e2 
                                   << "[exit_wt=" << exit_wt[e2] << "]" << endl;    }
                                   }    }    }    }

     // Delete edges.
                    
     UniqueSort(dels);
     if (logc.STATUS_LOGGING)
          cout << Date( ) << ": deleting " << dels.size( ) << " edges" << endl;
     DeleteEdges(dels);
     RemoveDeadEdgeObjects( );
     RemoveEdgelessVertices( );

     // Handle 'degenerate' case missed by above.  This is the case where there
     // is a vertex with one edge coming in and two going out.

     if ( heur.KILL_WEAK_EXITS_SOLO )
     {    int sever_count = 0;
          for ( int pass = 1; pass <= 2; pass++ )
          {    Reverse( );
               vec<vec<pair<int, int> > > paths_index( EdgeObjectCount( ) );
               for ( int i = 0; i < Paths( ).isize( ); i++ )
	       for ( int j = 0; j < Path(i).isize( ); j++ ) 
                    paths_index[Path(i,j)].push(i,j);
               for ( int v = 0; v < N( ); v++ )
               {    if ( To(v).size( ) != 1 || From(v).size( ) != 2 ) continue;
                    int e = EdgeObjectIndexByIndexTo( v, 0 );
                    vec<int> f(2);
                    f[0] = EdgeObjectIndexByIndexFrom( v, 0 );
                    f[1] = EdgeObjectIndexByIndexFrom( v, 1 );
                    vec<fix64_6> w(2, 0), s(2, 0);
                    for ( int i = 0; i < paths_index[e].isize( ); i++ )
                    {    int id = paths_index[e][i].first;
                         int j = paths_index[e][i].second;
                         for ( int m = 0; m < 2; m++ )
                              s[m] += Weight(id);
                         if ( j == Path(id).isize( ) - 1 ) continue;
                         for ( int m = 0; m < 2; m++ )
                              if ( Path(id)[j+1] == f[m] ) w[m] += Weight(id);    }
                    if ( w[0] < w[1] )
                    {    swap( f[0], f[1] );
                         swap( w[0], w[1] );
                         swap( s[0], s[1] );    }
                    const int solo_weak = 4;
                    if ( w[1] > weak ) continue;
                    if ( heur.KILL_WEAK_EXITS_SOLO_CAREFUL )
                    {    if ( Inv(f[0]) != f[1] && s[1] > solo_weak 
                              && From(v)[0] == From(v)[1] )
                         {    continue;    }    }
                    if ( w[0] < strong_mult * Max( w[1], fix64_6(1,2) ) ) continue;
                    if ( logc.verb[ "KILL_WEAK_EXITS" ] >= 1 )
                    {    cout << "KillWeakExits2, solo, severing connection from "
                              << "edge " << e << " to edge " << f[1] << endl;    }
                    sever_count++;
                    AddVertices(1);
                    GiveEdgeNewFromVx( f[1], v, N( ) - 1 );    }    }
          if (logc.STATUS_LOGGING)
          {    cout << Date( ) << ": " << sever_count << " connections severed" 
                    << endl;    }    }

     // Delete illegal paths.

     vec<int> to_left, to_right;
     ToLeft(to_left), ToRight(to_right);
     vec<Bool> to_delete( NPaths( ), False );
     for ( int j = 0; j < NPaths( ); j++ )
     {    for ( int i = 0; i < Path(j).isize( ) - 1; i++ )
          {    if ( to_right[ Path(j,i) ] != to_left[ Path(j,i+1) ] )
                    to_delete[j] = True;    }    }
     EraseIf( PathsMutable( ), to_delete );
     EraseIf( WeightsFwMutable( ), to_delete );
     EraseIf( WeightsRcMutable( ), to_delete );
     // Temporary if.
     if ( WeightsFwOrigin( ).size( ) == WeightsFw( ).size( ) )
     {    EraseIf( WeightsFwOriginMutable( ), to_delete );
          EraseIf( WeightsRcOriginMutable( ), to_delete );    }
     vec<Bool> p_to_delete( NPairs( ), False );
     for ( int i = 0; i < NPairs( ); i++ )
     {    for ( int j = 0; j < PairLeft(i).isize( ) - 1; j++ )
          {    if ( to_right[ PairLeft(i,j) ] != to_left[ PairLeft(i,j+1) ] )
                    p_to_delete[i] = True;    }
          for ( int j = 0; j < PairRight(i).isize( ) - 1; j++ )
          {    if ( to_right[ PairRight(i,j) ] != to_left[ PairRight(i,j+1) ] )
                    p_to_delete[i] = True;    }    }
     EraseIf( PairsMutable( ), p_to_delete );
     EraseIf( PairDataMutable( ), p_to_delete );

     // Finish up.

     RemoveUnneededVertices( );
     RemoveEdgelessVertices( );
     REPORT_TIME( clock, "used killing weak exits" );
     RemoveDeadEdgeObjects( );
     TestValid(logc);    }

namespace
{

template<int K> void OrientCore( SupportedHyperBasevector& shb, 
     const vecbasevector& genome )
{    vec< triple<kmer<K>,int,int> > kmers_plus;
     MakeKmerLookup1( genome, kmers_plus );
     int fw = 0, rc = 0;
     vec< kmer<K> > kmers( kmers_plus.size( ) );
     for ( int64_t i = 0; i < kmers_plus.jsize( ); i++ )
          kmers[i] = kmers_plus[i].first;
     vec<int> to_left, to_right;
     shb.ToLeft(to_left), shb.ToRight(to_right);
     vec< vec<int> > comp;
     shb.ComponentsE(comp);
     vec<Bool> reversed( shb.EdgeObjectCount( ), False );
     for ( int c = 0; c < comp.isize( ); c++ )
     {    if ( comp[c].empty( ) ) continue;
          if ( shb.Inv( comp[c][0] ) >= 0 ) continue;
          int fw = 0, rc = 0;
          for ( int i = 0; i < comp[c].isize( ); i++ )
          {    int e = comp[c][i];
               const basevector& E = shb.EdgeObject(e);
               vec<Bool> fwv( E.isize( ) - K + 1, False ); 
               vec<Bool> rcv( E.isize( ) - K + 1, False );
               #pragma omp parallel for
               for ( int j = 0; j <= E.isize( ) - K; j++ )
               {    kmer<K> x;
                    x.SetToSubOf( E, j );
                    if ( BinMember( kmers, x ) ) fwv[j] = True;
                    x.ReverseComplement( );
                    if ( BinMember( kmers, x ) ) rcv[j] = True;    }
               fw += Sum(fwv);
               rc += Sum(rcv);    }
          if ( rc > fw )
          {    int v = to_left[ comp[c][0] ];
               for ( int i = 0; i < comp[c].isize( ); i++ )
               {    int e = comp[c][i];
                    shb.EdgeObjectMutable(e).ReverseComplement( );
                    reversed[e] = True;    }
               shb.ReverseComponent(v);    }    }
     for ( int i = 0; i < shb.NPaths( ); i++ )
          if ( reversed[ shb.Path(i)[0] ] ) shb.PathMutable(i).ReverseMe( );

     // Reverse pairs.  Note that if a pair goes from one component to another,
     // and only one of the components is reversed, then the pair no longer makes
     // sense.  Not sure what to do about that.

     for ( int i = 0; i < shb.NPairs( ); i++ )
     {    Bool rev1 = reversed[ shb.PairLeft(i)[0] ];
          Bool rev2 = reversed[ shb.PairRight(i)[0] ];
          if (rev1) shb.PairLeftMutable(i).ReverseMe( );
          if (rev2) shb.PairRightMutable(i).ReverseMe( );
          if (rev1)
          {    shb.PairMutable(i) = make_pair(
                    shb.PairRight(i), shb.PairLeft(i) );    }    }    }

template <int K>
struct OrientCoreFunctor
{
    void operator()( SupportedHyperBasevector& shb,
                        const vecbasevector& genome )
    { OrientCore<K>(shb,genome); }
};

}

void OrientToReference( SupportedHyperBasevector& shb, const vecbasevector& genome,
     const long_logging& logc )
{    double clock = WallClockTime( );
     BigK::dispatch<OrientCoreFunctor>(shb.K(),shb,genome);
     shb.UniqueOrderPaths( );
     shb.TestValid(logc);
     REPORT_TIME( clock, "used orienting to reference" );    }

void AssessAssembly( const String& SAMPLE, const SupportedHyperBasevector& shb, 
     const HyperEfasta& he, const vec<Bool>& hide, const String& TMP, 
     const ref_data& ref, const String& HUMAN_CONTROLS,
     const long_logging& logc, const uint NUM_THREADS, RefTraceControl RTCtrl )
{
     const vecbasevector& G = ref.G;
     const vec<HyperBasevector>& GH = ref.GH; 
     const vec<bool>& is_circular = ref.is_circular;

     if (logc.STATUS_LOGGING) DATE_MSG( "assessing assembly" );
     if ( logc.COUNT_COV > 0 ) CountCov( shb, TMP, logc.COUNT_COV );
     int assembly_count = 0, reference_count = 0;
     if ( logc.READ_EVAL == "True" )
     {    double eclock = WallClockTime( );
          Bool print_a = False;
          vecbasevector bases( TMP + "/frag_reads_orig.fastb" );
          vecqualvector quals( TMP + "/frag_reads_orig.qualb" );
          vec<basevector> R;
          for ( int g = 0; g < (int) G.size( ); g++ )
               R.push_back( G[g] );
          HyperBasevector hb_R( shb.K( ), R );
          EvalByReads( shb, hb_R, bases, quals, assembly_count, reference_count, 
               print_a, logc.PRINT_FAVORING_REF );
          REPORT_TIME( eclock, "used evaluating by reads" );    }
     cout << "\n=================================================================="
          << "==================\n\n";
     cout << "SUMMARY STATS  --  dexter longread stats shown as -[...]-\n\n";
     if ( SAMPLE != "unknown" && logc.REFTRACE == "True" )
     {    RTCtrl.ReadySampleLookup( );
          ReadOriginTracker read_tracker(RTCtrl);
          RefTraceHeuristics rth = RefTraceHeuristics( );
          Bool fix_bug = False;
          NullOStream nullOS;
          if ( HUMAN_CONTROLS != "" )
          {    fix_bug = True;
               int np = 3;
               vec<int> penalty(np), gaps(np), meta_events(np);
               std::vector<ostringstream> out(np);
               for ( int pass = 0; pass < np; pass++ )
               {    if ( pass == 0 )
                    {    rth.min_group_frac = 0.1;
                         rth.min_group_save = 200;    }
                    if ( pass == 1 )
                    {    rth.max_offset_diff = 30;
                         rth.max_error_rate = 0.31;
                         rth.offset_add = 5;
                         rth.min_group_frac = 0.1;
                         rth.max_twiddle = 5;    }
                    if ( pass == 2 )
                    {    rth.max_offset_diff = 350;
                         rth.max_error_rate = 0.31;
                         rth.offset_add = 5;
                         rth.min_group_frac = 0.75;
                         rth.max_twiddle = 120;    }
                    RefTraceResults res = RefTrace( ref, shb, shb.Inv( ),
                            logc.verb[ "REFTRACE" ], logc, out[pass], rth, "",
                            fix_bug);
                    //RefTraceResults res = RefTrace( ref, shb, shb.Inv( ),
                    //        logc.verb[ "REFTRACE" ], logc, out[pass], nullOS,
                    //        rth, "", fix_bug, RTCtrl, &read_tracker);
                    penalty[pass] = res.penalty;
                    gaps[pass] = res.gaps;
                    meta_events[pass] = res.meta_events;    }
               int pix = -1;
               for ( int pi = 0; pi < np; pi++ )
               {    if ( penalty[pi] == Min(penalty) )
                    {    pix = pi;
                         cout << out[pi].str( );
                         if ( logc.SHOW_REFTRACE_EVENTS && logc.PERF_STATS )
                         {    PerfStat::log( ) << std::fixed << setprecision(0)
                                   << PerfStat( "error_meta_events", 
                                   "error meta-events", meta_events[pi] );
                              PerfStat::log( ) << std::fixed << setprecision(0)
                                   << PerfStat( "gaps", "gaps", gaps[pi] );    }
                         break;    }    }    
               for ( int pi = 0; pi < np; pi++ )
               {    if ( pi == pix ) continue;
                    istringstream in( out[pi].str( ) );
                    String line;
                    while(1)
                    {    getline( in, line );
                         if ( in.fail( ) ) break;
                         if ( line.Contains( " seconds " ) 
                              || line.Contains( " minutes " )
                              || line.Contains( " hours " ) 
                              || line.Contains( " days " ) )
                         {    cout << line << "\n";    }    }    }    }
          else {
            RefTraceAndCallVaraint( ref, shb, shb.Inv( ), logc.verb[ "REFTRACE"],
                    logc, cout, nullOS, rth, "", fix_bug, RTCtrl, &read_tracker); }
     }
     int hidden = 0;
     for ( int e = 0; e < shb.EdgeObjectCount( ); e++ )
          if ( shb.Inv(e) >= 0 && shb.Inv(e) < e ) hidden++;
     cout << shb.EdgeObjectCount( ) << " edges (" << shb.EdgeObjectCount( ) - hidden
          << " visible, " << hidden << " hidden)\n";
     cout << he.EdgeObjectCount( ) << " edges in efasta" 
          << " (" << he.EdgeObjectCount( ) - Sum(hide) << " visible, "
          << Sum(hide) << " hidden)" << endl;
     if ( !shb.Acyclic( ) ) cout << "assembly has a cycle\n";
     cout << "assembly has " << shb.NComponents( ) << " components\n";
     if ( G.size( ) > 0 )
     {    int64_t gsize = 0;
          for ( int g = 0; g < (int) G.size( ); g++ )
               gsize += G[g].size( );
          int excess_edges = shb.EdgeObjectCount( ) - (int) G.size( );
          double excess_edges_per = 1000000.0 * double(excess_edges) / double(gsize);
          cout << "-[" << setiosflags(ios::fixed) << setprecision(0) 
               << excess_edges_per 
               << resetiosflags(ios::fixed) << "]- excess edges per Mb" 
               << " (total excess edges = " << excess_edges << ")" << endl;
          if ( logc.PERF_STATS )
          {    PerfStat::log( ) << std::fixed << setprecision(0) 
                    << PerfStat( "excess_edges", "excess edges per Mb",
                    excess_edges_per );    }    }
     if ( SAMPLE != "unknown" ) ReportExcessEdges( he, G, logc.PERF_STATS );
     if ( logc.READ_EVAL == "True" )
     {    double eclock = WallClockTime( );
          int64_t gsize = 0;
          for ( int g = 0; g < (int) G.size( ); g++ )
               gsize += G[g].size( );
          double bad_reads_per = 1000000.0 * double(reference_count) / double(gsize);
          cout << "-[" << setiosflags(ios::fixed) << setprecision(1) << bad_reads_per
               << resetiosflags(ios::fixed) << "]- reads favoring reference per Mb" 
               << " (favors ref = " << reference_count << ")\n";
          if ( logc.PERF_STATS )
          {    PerfStat::log( ) << std::fixed << setprecision(1) << PerfStat( 
                    "bad_reads_per", "bad reads per Mb", bad_reads_per );    }
          vec<int> sources, sinks;
          shb.Sources(sources), shb.Sinks(sinks);
          int expected = 0;
          for ( int g = 0; g < (int) G.size( ); g++ )
               if ( !is_circular[g] ) expected += 2;
          int excess_ends = sources.isize( ) + sinks.isize( ) - expected;
          double excess_ends_per = double(excess_ends)/(double(gsize)/1000000.0);
          cout << "-[" << setiosflags(ios::fixed) << setprecision(2) 
               << excess_ends_per << resetiosflags(ios::fixed) 
               << "]- excess ends per Mb (excess ends = " << excess_ends << ")\n";
          if ( logc.PERF_STATS )
          {    PerfStat::log( ) << std::fixed << setprecision(2) 
                    << PerfStat( "excess_ends_per", "excess ends per Mb",
                    excess_ends_per );    }
          REPORT_TIME( eclock, "used evaluating by reads tail" );    }    }

void CountCov( const SupportedHyperBasevector& shb, const String& TMP,
     const int gp1 )
{
     // Determine file count.

     int fcount = 0;
     for ( fcount = 0; ; fcount++ )
          if ( !IsRegularFile( TMP + "/" + ToString(fcount) + ".fastb" ) ) break;

     // Set up counts.

     vec< vec<int> > cov( shb.EdgeObjectCount( ), vec<int>(2,0) );

     // Build data structures.

     const int L = 12;
     HyperBasevector hb_fw(shb), hb_rc(shb);
     hb_rc.Reverse( );
     vec<int> to_right_fw, to_right_rc;
     hb_fw.ToRight(to_right_fw), hb_rc.ToRight(to_right_rc);
     vecbasevector x_fw, x_rc;
     for ( int i = 0; i < hb_fw.EdgeObjectCount( ); i++ )
          x_fw.push_back( hb_fw.EdgeObject(i) );
     for ( int i = 0; i < hb_rc.EdgeObjectCount( ); i++ )
          x_rc.push_back( hb_rc.EdgeObject(i) );
     VecIntPairVec locs_fw, locs_rc;
     CreateGlocs( x_fw, L, locs_fw );
     CreateGlocs( x_rc, L, locs_rc );

     // Heuristics.

     const double prox = 20 * 1000;
     const int max_qual = 100 * 1000;
     const int min_pos_evidence = 8;

     // Go through the files.

     for ( int f = 0; f < fcount; f++ )
     {    vecbasevector bases( TMP + "/" +ToString(f) + ".fastb" );
          vecqualvector quals( TMP + "/" +ToString(f) + ".qualb" );
          int gid = ( f < gp1 ? 0 : 1 );

          // Align reads.

          #pragma omp parallel for
          for ( int64_t id = 0; id < (int64_t) bases.size( ); id++ )
          {    int n = KmerId( bases[id], L, 0 );
               vec<read_place> places;
               const int infinity = 1000000000;
               int qual_sum = infinity;
               const int min_qual = 1;
               FindPlaces( bases[id], quals[id], n, hb_fw, hb_rc, to_right_fw,
                    to_right_rc, locs_fw, locs_rc, places, qual_sum, min_qual,
                    prox );

               vec<Bool> to_delete( places.size( ), False );
               for ( int i = 0; i < places.isize( ); i++ )
                    if ( places[i].Qsum( ) > max_qual ) to_delete[i] = True;
               EraseIf( places, to_delete );
               Bool OK = False;
               if ( places.solo( ) ) OK = True;
               if ( places.size( ) == 2 && shb.InvDef( places[0].E( ).front( ) ) )
               {    vec<int> e = places[0].E( );
                    for ( int j = 0; j < e.isize( ); j++ )
                         e[j] = shb.Inv( e[j] );
                    // Doing this two ways, not sure which is right.
                    if ( e == places[1].E( ) ) OK = True;
                    e.ReverseMe( );
                    if ( e == places[1].E( ) ) OK = True;    }

               if (OK)
               {    for ( int i = 0; i < places.isize( ); i++ )
                    {    const read_place& p = places[i];
                         #pragma omp critical
                         {    for ( int j = 0; j < p.N( ); j++ )
                                   cov[ p.E(j) ][gid]++;    }    }    }    }    }

     // Announce results.

     cout << "\ncoverage:\n";
     for ( int e = 0; e < shb.EdgeObjectCount( ); e++ )
     {    if ( cov[e][0] < min_pos_evidence ) continue;
          if ( cov[e][1] > 0 ) continue;
          cout << e << " --> " << printSeq( cov[e] ) << "\n";    }    }

void SupportedHyperBasevector::RemoveSmallComponents2( const long_logging& logc,
     const int max_small_comp )
{    const int max_iterations = 1000;
     vec<int> e_to_delete;
     vec< vec<int> > comps;
     Components(comps);
     #pragma omp parallel for
     for ( size_t i = 0; i < comps.size( ); i++ )
     {    const vec<int>& o = comps[i];
          if ( HasCycle(o) ) continue;
          vec<int> sources, sinks;
          for ( int j = 0; j < o.isize( ); j++ )
          {    if ( Source( o[j] ) ) sources.push_back( o[j] );
               if ( Sink( o[j] ) ) sinks.push_back( o[j] );    }
          Bool failed = False;
          int max_path = 0;
          for ( int j1 = 0; j1 < sources.isize( ); j1++ )
          for ( int j2 = 0; j2 < sinks.isize( ); j2++ )
          {    vec<vec<int>> paths;
               // could upgrade to faster version that uses to_left, to_right
               Bool ok = EdgePaths( sources[j1], sinks[j2], paths, 
                    -1, -1, max_iterations );
               if ( !ok )
               {    failed = True;
                    break;    }
               for ( int l = 0; l < paths.isize( ); l++ )
               {    int mp = 0;
                    for ( int r = 0; r < paths[l].isize( ); r++ )
                         mp += EdgeLengthKmers( paths[l][r] );
                    max_path = Max( mp, max_path );    }    }
          if ( !failed && max_path <= max_small_comp )
          {    for ( size_t j = 0; j < o.size( ); j++ )
               {    int v = o[j];
                    #pragma omp critical
                    {    for ( size_t t = 0; t < From(v).size( ); t++ )
                         {    e_to_delete.push_back( EdgeObjectIndexByIndexFrom( 
                                   v, t ) );    }    }    }    }    }
     DeleteEdges(e_to_delete);
     TruncatePaths(logc);
     RemoveDeadEdgeObjects( );
     RemoveEdgelessVertices( );
     FixWeights(logc);
     TestValid(logc);    }

void SupportedHyperBasevector::RemoveSmallMostlyAcyclicComponents( 
     const long_logging& logc, const int max_small_comp )
{
     double clock = WallClockTime( );
     vec<int> e_to_delete;
     vec< vec<int> > comps;
     Components(comps);
     #pragma omp parallel for
     for ( size_t i = 0; i < comps.size( ); i++ )
     {    const vec<int>& o = comps[i];
          int nkmers = 0;
          int max_kmers = 0;
          for ( size_t j = 0; j < o.size( ); j++ )
          {    int v = o[j];
               for ( size_t t = 0; t < From(v).size( ); t++ )
               {    int n = EdgeObjectByIndexFrom( v, t ).size( ) - K( ) + 1;
                    nkmers += n;    
                    max_kmers = Max( max_kmers, n );    }    }
          const int min_cyclic_kmers_toss = 500;
          if ( nkmers > max_small_comp - K( ) + 1 
               || ( HasCycle(o) 
                    && ( nkmers > min_cyclic_kmers_toss || max_kmers > K( ) ) ) )
          {    continue;    }
          for ( size_t j = 0; j < o.size( ); j++ )
          {    int v = o[j];
               #pragma omp critical
               {    for ( size_t t = 0; t < From(v).size( ); t++ )
                    {    e_to_delete.push_back( EdgeObjectIndexByIndexFrom( 
                              v, t ) );    }    }    }    }
     DeleteEdges(e_to_delete);
     REPORT_TIME( clock, "removing mostly acyclic components" );
     TruncatePaths(logc);
     RemoveDeadEdgeObjects( );
     RemoveEdgelessVertices( );
     FixWeights(logc);
     TestValid(logc);    }

void SupportedHyperBasevector::RemovePathsWithoutReverseComplements()
{
     vec<Bool> delete_or_not( NPaths(), false);
     for ( int i1 = 0; i1 < NPaths( ); i1++ )
     {    const vec<int>& p1 = Path(i1);
          if ( !InvDef( p1[0] ) ) continue;
          vec<int> p2;
          for ( int j = 0; j < p1.isize( ); j++ ){
               if ( p1[j] < 0 ) p2.push_back( p1[j] );
               else             p2.push_back( Inv( p1[j] ) );
          }
          p2.ReverseMe( );
          int i2 = BinPosition( Paths( ), p2 );
          if ( i2 < 0 ) { delete_or_not[i1]=true; }
     }
     EraseIf(PathsMutable(),delete_or_not);
     EraseIf(WeightsRcMutable(),delete_or_not);
     EraseIf(WeightsFwMutable(),delete_or_not);
     EraseIf(WeightsFwOriginMutable(),delete_or_not);
     EraseIf(WeightsRcOriginMutable(),delete_or_not);
}
void SupportedHyperBasevector::RemovePairsWithoutReverseComplements()
{
    vec<Bool> delete_or_not( NPairs(), false);
    for ( int i1 = 0; i1 < NPairs( ); i1++ )
    {    const vec<int> &p1 = PairLeft(i1), &q1 = PairRight(i1);
         if ( !InvDef( p1[0] ) || !InvDef( q1[0] ) ) continue;
         vec<int> p2, q2;
         for ( int j = 0; j < p1.isize( ); j++ )
         {     if ( p1[j] < 0 ) p2.push_back( p1[j] );
              else p2.push_back( Inv( p1[j] ) );    }
         for ( int j = 0; j < q1.isize( ); j++ )
         {     if ( q1[j] < 0 ) q2.push_back( q1[j] );
              else q2.push_back( Inv( q1[j] ) );    }
         p2.ReverseMe( ), q2.ReverseMe( );
         int i2 = BinPosition( Pairs( ), make_pair( q2, p2 ) );
         if ( i2 < 0 )
         {    delete_or_not[i1]=true; }
    }
     EraseIf(PairsMutable(),delete_or_not);
     EraseIf(PairDataMutable(),delete_or_not);
}

void ReportAssemblyStats( const SupportedHyperBasevector& shb )
{    cout << Date( ) << ": initial assembly has K = " << shb.K( ) 
          << " and " << shb.EdgeObjectCount( ) << " edges" << endl;
     vec<int> nuni, nuni2;
     for ( int64_t id = 0; id < (int64_t) shb.NPaths( ); id++ )
     {    nuni.push_back( shb.Path(id).size( ) );
          if ( shb.Weight(id) >= 2 ) nuni2.push_back( shb.Path(id).size( ) );    }
     Sort(nuni), Sort(nuni2);
     DATE_MSG( "median path has " << Median(nuni) << " edges in it" );
     if ( nuni2.nonempty( ) )
     {    DATE_MSG( "median nonsolo path has " << Median(nuni2) 
               << " edges in it" );    }
     else DATE_MSG ("every unipath is solo");    }

void DumpEfastaAssembly( const SupportedHyperBasevector& shb, 
     const HyperEfasta& he, const vec<int>& inv, 
     const vec<Bool>& hide, const String& OUT_HEAD,
     const long_logging& logc )
{    double clock = WallClockTime( );
     Ofstream( dout, OUT_HEAD + ".dot" );
     const Bool DOT_LABEL_CONTIGS = True;
     const Bool DOT_LABEL_VERTICES = False;
     vec<double> lengths( he.EdgeObjectCount( ) );
     for ( int i = 0; i < he.EdgeObjectCount( ); i++ )
          lengths[i] = he.EdgeLengthKmers(i);
     vec<String> edge_id_names( he.EdgeObjectCount( ) );
     for ( int i = 0; i < he.EdgeObjectCount( ); i++ )
     {    if ( inv[i] < 0 ) edge_id_names[i] = ToString(i);
          else edge_id_names[i] = ToString(i) + "=" + ToString( inv[i] ) + "'";    }
     he.PrettyDOT( dout, lengths, HyperEfasta::edge_label_info( 
          HyperEfasta::edge_label_info::DIRECT, &edge_id_names ), 
          DOT_LABEL_CONTIGS, DOT_LABEL_VERTICES, NULL, NULL, NULL, &hide );
     Ofstream( eout, OUT_HEAD + ".efasta" );
     for ( int v = 0; v < he.N( ); v++ )
     {    for ( size_t j = 0; j < he.From(v).size(); j++ )
          {    int e = he.EdgeObjectIndexByIndexFrom( v, j );
               int w = he.From(v)[j];
               he.EdgeObject(e).Print( eout, ToString(e) + " [vert_" + ToString(v) 
                    + "-->vert_" + ToString(w) + "]" );    }    }
     BinaryWriter::writeFile( OUT_HEAD + ".hbv", (HyperBasevector) shb );
     REPORT_TIME( clock, "used dumping efasta" );    }

// Detect variants by finding all bubbles and alignments reads back to the
// branches. Currently code is copied from DivineBubbles.

void SupportedHyperBasevector::DetectVarients( const vecbasevector& bases,
        const vecqualvector& quals,
             const long_logging& logc, vec<VariantSignature>* snp_bubbles,
             vec<int>* snp_edge_list) const
{
     // Find the bubbles.

     double clock1 = WallClockTime( );
     Bool verbose = logc.verb[ "DETECT_VARIENTS" ];
     vec< vec<int> > bubbles;
     for ( int v = 0; v < N( ); v++ )
     {    if ( To(v).size( ) != 1 || From(v).size( ) != 2 ) continue;
          if ( From(v)[0] != From(v)[1] ) continue;
          int w = From(v)[0];
          if ( To(w).size( ) != 2 || From(w).size( ) != 1 ) continue;
          int x = To(v)[0], y = From(w)[0];
          vec<int> all, b;
          all.push_back( x, v, w, y );
          UniqueSort(all);
          if ( all.size( ) != 4 ) continue;
          int e = EdgeObjectIndexByIndexTo( v, 0 );
          int f1 = EdgeObjectIndexByIndexFrom( v, 0 );
          int f2 = EdgeObjectIndexByIndexFrom( v, 1 );
          if ( f1 > f2 ) swap( f1, f2 );
          int g = EdgeObjectIndexByIndexFrom( w, 0 );
          b.push_back( e, f1, f2, g );
          bubbles.push_back(b);    }

     // Build data structures.

     const int L = 12;
     HyperBasevector hb_fw(*this), hb_rc(*this);
     hb_rc.Reverse( );
     vec<int> to_right_fw, to_right_rc;
     hb_fw.ToRight(to_right_fw), hb_rc.ToRight(to_right_rc);
     vecbasevector x_fw, x_rc;
     for ( int i = 0; i < hb_fw.EdgeObjectCount( ); i++ )
          x_fw.push_back( hb_fw.EdgeObject(i) );
     for ( int i = 0; i < hb_rc.EdgeObjectCount( ); i++ )
          x_rc.push_back( hb_rc.EdgeObject(i) );
     VecIntPairVec locs_fw, locs_rc;
     CreateGlocs( x_fw, L, locs_fw );
     CreateGlocs( x_rc, L, locs_rc );

     // Go through the reads.

     vec< vec<fix64_6> > votes_bub( bubbles.size( ), vec<fix64_6>(4,0) );
     vec< triple< int, int, pair<fix64_6,int> > > adj;
     vec<fix64_6> support( EdgeObjectCount( ), 0 );
     #pragma omp parallel for
     for ( int64_t id = 0; id < (int64_t) bases.size( ); id++ )
     {    int n = KmerId( bases[id], L, 0 );
          vec<read_place> places;
          const int infinity = 1000000000;
          int qual_sum = infinity;
          FindPlaces( bases[id], quals[id], n, hb_fw, hb_rc, to_right_fw, 
               to_right_rc, locs_fw, locs_rc, places, qual_sum );
          int np = places.size( );

          #pragma omp critical
          {    for ( int i = 0; i < np; i++ )
               {    for ( int j = 0; j < places[i].N( ); j++ )
                         support[ places[i].E(j) ] += fix64_6(1,np);
                    for ( int j = 0; j < places[i].N( ) - 1; j++ )
                    {    if ( places[i].Fw( ) )
                         {    adj.push( places[i].E(j), places[i].E(j+1), 
                                   make_pair( fix64_6(1,np), qual_sum ) );    }
                         else
                         {    adj.push( places[i].E(j+1), places[i].E(j), 
                                   make_pair( fix64_6(1,np), qual_sum ) );    }
                                   }    }    }

          #pragma omp critical
          for ( int i = 0; i < np; i++ )
          {    vec<int> bids;
               for ( int j = 0; j < bubbles.isize( ); j++ )
               {    if ( Member( places[i].E( ), bubbles[j][1] ) )
                    {    if ( places[i].Fw( ) ) votes_bub[j][0] += fix64_6(1,np);
                         else votes_bub[j][1] += fix64_6(1,np);    }
                    if ( Member( places[i].E( ), bubbles[j][2] ) )
                    {    if ( places[i].Fw( ) ) votes_bub[j][2] += fix64_6(1,np);
                         else votes_bub[j][3] += fix64_6(1,np);    }
                    if ( Member( places[i].E( ), bubbles[j][1] )
                         || Member( places[i].E( ), bubbles[j][2] ) )
                    {     bids.push_back(j);    }    }
               if ( bids.empty( ) ) continue;  }    }
     
     // Test bubbles.

     const double max_asym_rarity = 0.00001;
     const int min_to_save = 10;
     if (verbose) cout << "\nbubbles:\n";
     for ( int i = 0; i < bubbles.isize( ); i++ )
     {    double f1 = votes_bub[i][0].ToDouble( ), r1 = votes_bub[i][1].ToDouble( );
          double f2 = votes_bub[i][2].ToDouble( ), r2 = votes_bub[i][3].ToDouble( );
          int e1 = 1, e2 = 2;
          if ( f2 + r2 > f1 + r1 || ( f2 + r2 == f1 + r1 && f2 > f1 ) )
          {    swap( e1, e2 );
               swap( f1, f2 );
               swap( r1, r2 );    }
          if ( f2 > r2 || ( f2 == r2 && f1 > r1 ) )
          {    swap( f1, r1 ); 
               swap( f2, r2 );    }
          long double p = Min( 0.5, f1/(f1+r1) ) / 2;
          int n = int(floor(f1+r1+f2+r2));
          double q = -1;
          if ( n > 0 && n <= 10000 ) 
          {    q = BinomialSum( n, int(ceil(f2)), p );
               //if ( q < max_asym_rarity && f2 + r2 < min_to_save )
               //{  // } unlikely bubbles are deleted
              //int branch1 = (e1 == 1 ? bubbles[i][1] : bubbles[i][2]);
              //int branch2 = (e1 == 1 ? bubbles[i][2] : bubbles[i][1]);
              if ( q > max_asym_rarity ) {
                  int branch1 = bubbles[i][1];
                  int branch2 = bubbles[i][2];
                  snp_bubbles->push( x_fw[branch1], x_fw[branch2] );
                  if ( snp_edge_list != NULL ) {
                      snp_edge_list->push_back( branch1 );
                      snp_edge_list->push_back( branch2 );
                  }
              }
          }
          if (verbose)
          {    cout << "[" << i << "] " << bubbles[i][0] << " --> {" 
                    << bubbles[i][1] << "," << bubbles[i][2] << "} --> " 
                    << bubbles[i][3] << "; vote: " << setiosflags(ios::fixed) 
                    << setprecision(1) << f1 << "+" << r1 << " vs " << f2 << "+" 
                    << r2 << resetiosflags(ios::fixed);
               if ( q >= 0 ) cout << ", surprise(" << e2 << ") = " << q;
               cout << endl;    }    }

     if (logc.STATUS_LOGGING)
          cout << Date( ) << ": found " << snp_bubbles->size( ) << " bubbles" << endl;
     REPORT_TIME( clock1, "used detecting variants" );    }

void TraceEdges( const SupportedHyperBasevector& shb, const String& TRACE_EDGES,
     const vecbasevector& bases, const vecqualvector& quals )
{    if ( TRACE_EDGES == "" ) return;
     cout << "\nplacements of reads on trace edges:\n";
     const int infinity = 1000000000;
     vec<int> trace, interesting;
     ParseIntSet( TRACE_EDGES, trace );
     const int L = 12;
     HyperBasevector hb_fw(shb), hb_rc(shb);
     hb_rc.Reverse( );
     vec<int> to_left_fw, to_left_rc;
     hb_fw.ToLeft(to_left_fw), hb_rc.ToLeft(to_left_rc);
     vec<int> to_right_fw, to_right_rc;
     hb_fw.ToRight(to_right_fw), hb_rc.ToRight(to_right_rc);

     // Have we been given two edges to trace and do they have the
     // same start and end points?

     Bool bubble = False;
     if ( trace.size( ) == 2 && to_left_fw[ trace[0] ] == to_left_fw[ trace[1] ]
          && to_right_fw[ trace[0] ] == to_right_fw[ trace[1] ] )
     {    bubble = True;    }

     vecbasevector x_fw, x_rc;
     for ( int i = 0; i < hb_fw.EdgeObjectCount( ); i++ )
          x_fw.push_back( hb_fw.EdgeObject(i) );
     for ( int i = 0; i < hb_rc.EdgeObjectCount( ); i++ )
          x_rc.push_back( hb_rc.EdgeObject(i) );
     VecIntPairVec locs_fw, locs_rc;
     CreateGlocs( x_fw, L, locs_fw );
     CreateGlocs( x_rc, L, locs_rc );
     vec< vec< vec<fix64_6> > > count( trace.size( ) );
     for ( int i = 0; i < trace.isize( ); i++ )
          count[i].resize( shb.EdgeObject( trace[i] ).size( ), vec<fix64_6>(4,0) );
     #pragma omp parallel for
     for ( int64_t id = 0; id < (int64_t) bases.size( ); id++ )
     {    int n = KmerId( bases[id], L, 0 );
          vec<read_place> places;
          int qual_sum = infinity;
          FindPlaces( bases[id], quals[id], n, hb_fw, hb_rc, to_right_fw, 
               to_right_rc, locs_fw, locs_rc, places, qual_sum );
          int np = places.size( );
          for ( int i = 0; i < np; i++ )
          {    vec<int> e = places[i].E( );
               UniqueSort(e);
               if ( Meet( trace, e ) )
               {    
                    #pragma omp critical
                    {    const read_place& p = places[i];
                         const HyperBasevector& hb = ( p.Fw( ) ? hb_fw : hb_rc );
                         const basevector& b = bases[id];
                         const qualvector& q = quals[id];
                         vec<Bool> match( b.size( ) );
                         const int flank = 5;
                         for ( int pass = 1; pass <= 2; pass++ )
                         {    int ei = 0, pos = p.P( );
                              for ( int l = 0; l < b.isize( ); l++ )
                              {    int x = BinPosition( trace, p.E(ei) );
                                   const basevector& E = hb.EdgeObject( p.E(ei) );
                                   if ( pass == 1 ) match[l] = ( b[l] == E[pos] );
                                   if ( x >= 0 && q[l] > 2 && pass == 2 )
                                   {    int n = E.size( );

                                        Bool flanked = False, mismatch = False;
                                        if ( l >= flank && b.isize( ) - l > flank )
                                        {    for ( int m = 1; m <= flank; m++ )
                                             {    if ( !match[l-m] || !match[l+m] )
                                                       mismatch = True;    }
                                             flanked = !mismatch;    }

                                        if (flanked)
                                        {
                                        if ( p.Fw( ) )
                                        {    count[x][pos][ b[l] ] 
                                                  += fix64_6( q[l], np );    }
                                        else
                                        {    count[x][n-pos-1][ 3 - b[l] ]
                                                  += fix64_6( q[l], np );    }    
                                             }
                                        }

                                   pos++;
                                   if ( pos == hb.EdgeObject( p.E(ei) ).isize( ) )
                                   {    ei++;
                                        if ( ei == p.N( ) ) break;
                                        pos = hb.K( ) - 1;    }    }    }    }
                    int id1 = id/2 * 2;
                    int id2 = id1 + 1;
                    #pragma omp critical
                    {    interesting.push_back( id1, id2 );    }
                    break;    }    }    }
     for ( int i = 0; i < trace.isize( ); i++ )
     {    cout << "\ncoverage of edge " << trace[i] << endl;
          int n = count[i].size( );
          vec<Bool> interesting(n);
          for ( int j = 0; j < n; j++ )
          {    int bigs = 0;
               for ( int k = 0; k < 4; k++ )
                    if ( count[i][j][k] >= 50 ) bigs++;
               interesting[j] = ( bigs >= 2);    }
          vec<Bool> near_interesting( n, False );
          for ( int j = 0; j < n; j++ )
          {    if ( interesting[j] )
               {    for ( int k = -10; k <= 10; k++ )
                    {    if ( j+k >= 0 && j+k < n )
                              near_interesting[j+k] = True;    }    }    }
          int last_printed = -1;
          for ( int j = 0; j < n; j++ )
          {    if ( !near_interesting[j] ) continue;
               if ( j > 0 && last_printed < 0 ) cout << "\n...\n";
               else if ( j > last_printed + 1 ) cout << "...\n\n...\n";
               last_printed = j;
               cout << j << " -->";
               int bigs = 0;
               for ( int k = 0; k < 4; k++ )
               {    if ( count[i][j][k] >= 50 ) bigs++;
                    if ( count[i][j][k] == 0 ) continue;
                    cout << " " << as_base(k) << "[" << count[i][j][k] << "]";    }
               if ( bigs >= 2 ) cout << " *****";
               cout << "\n";    }
          if ( last_printed < n - 1 ) cout << "...\n";    }
     cout << "\n";
     UniqueSort(interesting);
     for ( int j = 0; j < interesting.isize( ); j++ )
     {    int id = interesting[j];
          int n = KmerId( bases[id], L, 0 );
          vec<read_place> places;
          int qual_sum = infinity;
          FindPlaces( bases[id], quals[id], n, hb_fw, hb_rc, to_right_fw, 
               to_right_rc, locs_fw, locs_rc, places, qual_sum );

          // In the case of a bubble, create places for the 'alternate' placement,
          // if not already computed.  This doesn't properly handle the case where
          // the read goes through the bubble more than once.

          if (bubble)
          {    int np = places.size( );
               for ( int i = 0; i < np; i++ )
               {    for ( int m = 0; m < 2; m++ )
                    {    int p = Position( places[i].E( ), trace[m] );
                         if ( p < 0 ) continue;
                         read_place x = places[i];
                         x.EMutable( )[p] = trace[1-m];
                         Bool found = False;
                         for ( int l = 0; l < places.isize( ); l++ )
                         {    if ( places[l].E( ) == x.E( )
                                   && places[l].P( ) == x.P( )
                                   && places[l].Fw( ) == x.Fw( ) )
                              {    found = True;
                                   break;    }    }
                         if (found) continue;
                         const int min_qual = 3;
                         x.ComputeQsum( bases[id], quals[id],
                              ( x.Fw( ) ? hb_fw : hb_rc ), min_qual );
                         places.push_back(x);    }    }    }

          if ( id % 2 == 0 ) cout << "\n";
          int np = places.size( );
          {    for ( int i = 0; i < np; i++ )
                    cout << id << "." << i+1 << ": " << places[i] << endl;    }    }
     cout << "\n";    }

void AnalyzeAssembly( const SupportedHyperBasevector& shb, const vecbasevector& G,
     const int LG, const VecIntPairVec& Glocs )
{    cout << "\nAlignment of assembly to reference:\n\n";
     for ( int e = 0; e < shb.EdgeObjectCount( ); e++ )
     {    vec<look_align> aligns;
          ClusterAligner( shb.EdgeObject(e), G, LG, Glocs, aligns );
          for ( int j = 0; j < aligns.isize( ); j++ )
          {    look_align& la = aligns[j];
               la.query_id = e;
               int g = la.target_id;
               la.PrintReadableBrief( cout, "edge_" + ToString(e), ToString(g) );
               basevector query = shb.EdgeObject(e);
               la.PrintVisual( cout, query, G[g] );    
               if (la.rc1) query.ReverseComplement( );
               la.a.PrintMutations( query, G[g], cout );    
               cout << "\n";    }    }    }
