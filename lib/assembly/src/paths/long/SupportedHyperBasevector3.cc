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
#include "lookup/LookAlign.h"
#include "math/Functions.h"
#include "pairwise_aligners/ClusterAligner.h"
#include "pairwise_aligners/SmithWatFree.h"
#include "paths/long/CreateGenome.h"
#include "paths/long/EvalByReads.h"
#include "paths/long/SupportedHyperBasevector.h"
#include "paths/long/DiscovarTools.h"
#include "random/Bernoulli.h"

void RemoveNegatives( vec<int>& p )
{    vec<Bool> to_delete( p.size( ), False );
     for ( int j = 0; j < p.isize( ); j++ )
          if ( p[j] == -1 ) to_delete[j] = True;
     EraseIf( p, to_delete );    }

Bool Compati( const vec<int>& X, const vec<int>& v, const int r, const int pos )
{    Bool compati = False;
     for ( int rr = 1; rr < v.isize( ) - 1; rr++ )
     {    if ( v[rr] != v[r] ) continue;
          Bool compatj = True;
          for ( int s = rr + 1; s < v.isize( ); s++ )
          {    int p = pos + s - rr;
               if ( p >= X.isize( ) ) break;
               if ( v[s] != X[p] )
               {    compatj = False;
                    break;    }    }
          for ( int s = rr - 1; s >= 0; s-- )
          {    int p = pos + s - rr;
               if ( p < 0 ) break;
               if ( v[s] != X[p] )
               {    compatj = False;
                    break;    }    }
          if (compatj)
          {    compati = True;
               break;    }    }
     return compati;    }

#if 0
void SupportedHyperBasevector::PullApart2( const double min_weight_split,
     const long_logging& logc )
{
     // Logging stuff.

     if (logc.STATUS_LOGGING) cout << Date( ) << ": starting PullApart2" << endl;
     int verbosity = logc.verb[ "PULL_APART2" ];
     double clock = WallClockTime( );

     // Iterate until no improvement.

     while(1)
     {    Bool progress = False;
          // ---->

     // Set up indices.

     vec<int> to_left, to_right;
     ToLeft(to_left), ToRight(to_right);
     vec< vec< pair<int,int> > > paths_index( EdgeObjectCount( ) );
     for (int i = 0; i < Paths( ).isize(); i++)    
     for (int j = 0; j < Path(i).isize(); j++)
          paths_index[ Path(i,j) ].push(i, j);

     // Keep track of proposed joins.

     vec< triple< vec<int>, vec<int>, vec<int> > > joins;
     vec< vec< vec<int> > > paths1, paths2;

     // Start with an edge e1, ending in a vertex v.

     int count = 0;
     #pragma omp parallel for schedule(dynamic, 1)
     for ( int e1 = 0; e1 < EdgeObjectCount( ); e1++ )
     {    
          if (logc.STATUS_LOGGING)
          {
               #pragma omp critical
               {    int dots1 = (100*count)/EdgeObjectCount( );
                    count++;
                    int dots2 = (100*count)/EdgeObjectCount( );
                    if ( dots2 > dots1 )
                    {    for ( int j = dots1; j < dots2; j++ )
                         {    cout << ".";
                              if ( (j+1) % 50 == 0 ) cout << "\n";
                              else if ( (j+1) % 10 == 0 ) cout << " ";    }
                         flush(cout);    }    }    }

          int v = to_right[e1];

          // Form an edge neighborhood of v.

          const int radius = 4;
          vec< pair<int,int> > ed;
          ed.push( e1, 0 );
          while(1)
          {    Bool progress = False;
               for ( int i = 1; i < ed.isize( ); i++ )
               {    if ( ed[i].second == radius ) continue;
                    int w = to_left[ ed[i].first ];
                    for ( int j = 0; j < To(w).isize( ); j++ )
                    {    int e = EdgeObjectIndexByIndexTo( w, j );
                         int d = ed[i].second + 1;
                         Bool found = False;
                         for ( int l = 0; l < ed.isize( ); l++ )
                         {    if ( ed[l].first == e )
                              {    found = True;
                                   if ( ed[l].second > d )
                                   {    ed[l].second = d;
                                        progress = True;    }    }    }
                         if ( !found )
                         {    ed.push( e, d );
                              progress = True;    }    }    }
               for ( int i = 0; i < ed.isize( ); i++ )
               {    if ( ed[i].second == radius ) continue;
                    int w = to_right[ ed[i].first ];
                    for ( int j = 0; j < From(w).isize( ); j++ )
                    {    int e = EdgeObjectIndexByIndexFrom( w, j );
                         int d = ed[i].second + 1;
                         Bool found = False;
                         for ( int l = 0; l < ed.isize( ); l++ )
                         {    if ( ed[l].first == e )
                              {    found = True;
                                   if ( ed[l].second > d )
                                   {    ed[l].second = d;
                                        progress = True;    }    }    }
                         if ( !found )
                         {    ed.push( e, d );
                              progress = True;    }    }    }
               if ( !progress ) break;    }

          // Define (e1,f1) and (e2,f2) pairs.

          vec< pair<int,int> > s1, s2;
          {    vec<int> x, s;
               for ( int l = 0; l < ed.isize( ); l++ )
                    x.push_back( ed[l].first );
               UniqueSort(x);
               s.push_back(e1);
               for ( int l = 0; l < radius; l++ )
               {    int n = s.size( );
                    for ( int m = 0; m < n; m++ )
                         s.append( FromEdgeObj( to_right[ s[m] ] ) );
                    UniqueSort(s);    }
               for ( int l = 0; l < s.isize( ); l++ )
                    s1.push( e1, s[l] );
               for ( int l = 0; l < x.isize( ); l++ )
               {    int e2 = x[l];
                    vec<int> t;
                    t.push_back(e2);
                    for ( int k = 0; k < radius; k++ )
                    {    int n = t.size( );
                         for ( int m = 0; m < n; m++ )
                              t.append( FromEdgeObj( to_right[ t[m] ] ) );
                         UniqueSort(t);    }
                    t = Intersection( t, x );
                    for ( auto m : t ) s2.push( e2, m );    }    }

          // Define e2, f1, f2.
          
          for ( int me = 0; me < s1.isize( ); me++ )
          for ( int mf = 0; mf < s2.isize( ); mf++ )
          {    int f1 = s1[me].second, e2 = s2[mf].first, f2 = s2[mf].second;
               vec<int> all;
               all.push_back( e1, e2, f1, f2 );
               UniqueSort(all);
               if ( all.size( ) != 4 ) continue;

               // Find edges reachable from e1 and e2, and bridge paths.

               vec<int> reach;
               Bool found_f1 = False, found_f2 = False;
               vec< vec<int> > bridge11, bridge22, bridge12, bridge21;
               Bool bad = False;
               for ( int z = 1; z <= 2; z++ )
               {    int e = ( z == 1 ? e1 : e2 );
                    for ( int i = 0; i < paths_index[e].isize( ); i++ )
                    {    int id = paths_index[e][i].first; 
                         int p = paths_index[e][i].second;
                         for ( int j = p + 1; j < Path(id).isize( ); j++ )
                         {    int g = Path(id)[j];
                              if ( e == e1 && g == e2 ) bad = True;
                              if ( e == e2 && g == e1 ) bad = True;
                              if ( e == e1 && g == f1 )
                              {    vec<int> b;
                                   for ( int l = p; l <= j; l++ )
                                        b.push_back( Path(id)[l] );
                                   bridge11.push_back(b);    }
                              if ( e == e2 && g == f2 )
                              {    vec<int> b;
                                   for ( int l = p; l <= j; l++ )
                                        b.push_back( Path(id)[l] );
                                   bridge22.push_back(b);    }
                              if ( e == e1 && g == f2 )
                              {    vec<int> b;
                                   for ( int l = p; l <= j; l++ )
                                        b.push_back( Path(id)[l] );
                                   bridge12.push_back(b);    }
                              if ( e == e2 && g == f1 )
                              {    vec<int> b;
                                   for ( int l = p; l <= j; l++ )
                                        b.push_back( Path(id)[l] );
                                   bridge21.push_back(b);    }
                              if ( g == f1 ) 
                              {    found_f1 = True;
                                   for ( int k = j + 1; k < Path(id).isize( ); k++ )
                                   {    if ( e == e1 && Path(id)[k] == e2 ) 
                                             bad = True;
                                        if ( e == e2 && Path(id)[k] == e1 ) 
                                             bad = True;    }
                                   break;    }
                              if ( g == f2 ) 
                              {    found_f2 = True;
                                   for ( int k = j + 1; k < Path(id).isize( ); k++ )
                                   {    if ( e == e1 && Path(id)[k] == e2 ) 
                                             bad = True;
                                        if ( e == e2 && Path(id)[k] == e1 ) 
                                             bad = True;    }
                                   break;    }
                              reach.push_back(g);    }    }    }
               if ( bad || !found_f1 || !found_f2 ) continue;
               UniqueSort(reach), UniqueSort(bridge11), UniqueSort(bridge22);
               if ( BinMember( reach, e1 ) || BinMember( reach, e2 ) ) continue;

               // Check for paths that are consistent with no bridge.

               vec< vec<int> > all_bridges;
               all_bridges.append(bridge11);
               all_bridges.append(bridge22);
               all_bridges.append(bridge12);
               all_bridges.append(bridge21);
               for ( int z = 1; z <= 2; z++ )
               {    int e = ( z == 1 ? e1 : e2 );
                    for ( int i = 0; i < paths_index[e].isize( ); i++ )
                    {    int id = paths_index[e][i].first; 
                         int p = paths_index[e][i].second;
                         Bool term = False;
                         for ( int j = p + 1; j < Path(id).isize( ); j++ )
                         {    if ( Path(id)[j] == f1 || Path(id)[j] == f2 )
                                   term = True;    }
                         if (term) continue;
                         vec<int> b;
                         for ( int j = p; j < Path(id).isize( ); j++ )
                              b.push_back( Path(id)[j] );
                         Bool normal = False;
                         for ( int i = 0; i < all_bridges.isize( ); i++ )
                              if ( all_bridges[i].Contains( b, 0 ) ) normal = True;
                         if ( !normal) bad = True;    }    }
               if (bad) continue;

               // See if {e1,e2}{reach}{f1,f2} makes sense.  Note that we're 
               // requiring all intermediate edges to be reachable by paths 
               // starting at e1 or e2, which is probably too stringent.

               for ( int i = 0; i < reach.isize( ); i++ )
               {    int x = to_left[ reach[i] ], y = to_right[ reach[i] ];
                    for ( int j = 0; j < To(x).isize( ); j++ )
                    {    int g = EdgeObjectIndexByIndexTo( x, j );
                         if ( g != e1 && g != e2 && !BinMember( reach, g ) )
                              bad = True;    }
                    for ( int j = 0; j < From(y).isize( ); j++ )
                    {    int g = EdgeObjectIndexByIndexFrom( y, j );
                         if ( g != f1 && g != f2 && !BinMember( reach, g ) )
                              bad = True;    }    }
               for ( int i = 0; i < 2; i++ )
               {    int x = ( i == 0 ? to_left[f1] : to_left[f2] );
                    for ( int j = 0; j < To(x).isize( ); j++ )
                    {    int g = EdgeObjectIndexByIndexTo( x, j );
                         if ( g != e1 && g != e2 && !BinMember( reach, g ) )
                              bad = True;    }
                    int y = ( i == 0 ? to_right[e1] : to_right[e2] );
                    for ( int j = 0; j < From(y).isize( ); j++ )
                    {    int g = EdgeObjectIndexByIndexFrom( y, j );
                         if ( g != f1 && g != f2 && !BinMember( reach, g ) )
                              bad = True;    }    }
               if (bad) continue;

               // The reach edges must comprise a connected subgraph.

               if ( reach.nonempty( ) )
               {    vec<int> r;
                    r.push_back( reach[0] );
                    for ( int i = 0; i < r.isize( ); i++ )
                    {    int v = to_left[ r[i] ], w = to_right[ r[i] ];
                         for ( int pass = 1; pass <= 2; pass++ )
                         {    int x = ( pass == 1 ? v : w );
                              for ( int j = 0; j < From(x).isize( ); j++ )
                              {    int e = EdgeObjectIndexByIndexFrom( x, j );
                                   if ( BinMember( reach, e ) && !Member( r, e ) )
                                        r.push_back(e);    }
                              for ( int j = 0; j < To(x).isize( ); j++ )
                              {    int e = EdgeObjectIndexByIndexTo( x, j );
                                   if ( BinMember( reach, e ) && !Member( r, e ) )
                                        r.push_back(e);    }    }    }
                    if ( r.size( ) < reach.size( ) ) continue;    }

               // Compute maximum bridge length.

               int B = 0;
               for ( int i = 0; i < bridge11.isize( ); i++ )
               {    int b = 0;
                    for ( int j = 1; j < bridge11[i].isize( ) - 1; j++ )
                         b += EdgeLengthKmers( bridge11[i][j] );
                    B = Max( B, b );    }
               for ( int i = 0; i < bridge22.isize( ); i++ )
               {    int b = 0;
                    for ( int j = 1; j < bridge22[i].isize( ) - 1; j++ )
                         b += EdgeLengthKmers( bridge22[i][j] );
                    B = Max( B, b );    }

               // See if cross-paths satisfy hypotheses.

               fix64_6 weight_11 = 0, weight_12 = 0, weight_21 = 0, weight_22 = 0;
               for ( int i = 0; i < paths_index[e1].isize( ); i++ )
               {    int id = paths_index[e1][i].first, p = paths_index[e1][i].second;
                    for ( int j = p + 1; j < Path(id).isize( ); j++ )
                    {    if ( Path(id)[j] == f1 ) weight_11 += Weight(id);
                         if ( Path(id)[j] == f2 ) weight_12 += Weight(id);    }    }
               for ( int i = 0; i < paths_index[e2].isize( ); i++ )
               {    int id = paths_index[e2][i].first, p = paths_index[e2][i].second;
                    for ( int j = p + 1; j < Path(id).isize( ); j++ )
                    {    if ( Path(id)[j] == f1 ) weight_21 += Weight(id);
                         if ( Path(id)[j] == f2 ) weight_22 += Weight(id);    }    }
               const int min_weight_split_low = 2;
               Bool OK = False;
               if ( weight_11 >= min_weight_split_low
                    && weight_22 >= min_weight_split_low
                    && weight_12 == 0 && weight_21 == 0 )
               {    OK = True;    } 
               if ( weight_11 >= min_weight_split && weight_22 >= min_weight_split
                    && weight_12 + weight_21 <= 2
                    && B <= MedianCorrectedReadLengthFudge( ) )
               {    OK = True;    }
               if ( weight_11 >= min_weight_split/2 
                    && weight_22 >= min_weight_split/2
                    && weight_12 + weight_21 < 2
                    && B <= MedianCorrectedReadLengthFudge( ) )
               {    OK = True;    }
               if ( ( OK && verbosity >= 2 ) || verbosity >= 3 )
               {    
                    #pragma omp critical
                    {    cout << "\n";
                         PRINT4( e1, e2, f1, f2 );
                         cout << "reach = " << printSeq(reach) << endl;
                         cout << "bridge11:\n";
                         for ( int i = 0; i < bridge11.isize( ); i++ )
                         {    cout << "[" << i+1 << "] " << printSeq( bridge11[i] ) 
                                   << "\n";    }
                         cout << "bridge22:\n";
                         for ( int i = 0; i < bridge22.isize( ); i++ )
                         {    cout << "[" << i+1 << "] " << printSeq( bridge22[i] ) 
                                   << "\n";    }
                         if ( verbosity >= 3 && !OK ) 
                         {    PRINT4( weight_11, weight_22, weight_12, weight_21 );
                              cout << "rejecting" << endl;    }    }    }
               if ( !OK ) continue;
     
               // Save join.

               vec<int> e12, f12;
               e12.push_back(e1,e2), f12.push_back(f1,f2);
               #pragma omp critical
               {    joins.push( e12, reach, f12 );
                    paths1.push_back(bridge11); 
                    paths2.push_back(bridge22);    }    }    }

     // Process joins.

     ParallelSortSync( joins, paths1, paths2 );
     if (logc.STATUS_LOGGING)
     {    cout << Date( ) << ": processing " << joins.size( ) 
               << " potential joins" << endl;    }
     if ( verbosity >= 2 ) cout << "\n";
     vec<Bool> touched( EdgeObjectCount( ), False );
     for ( int i = 0; i < joins.isize( ); i++ )
     {    Bool overlap = False;
          for ( int j = 0; j < joins[i].first.isize( ); j++ )
               if ( touched[ joins[i].first[j] ] ) overlap = True;
          for ( int j = 0; j < joins[i].second.isize( ); j++ )
               if ( touched[ joins[i].second[j] ] ) overlap = True;
          for ( int j = 0; j < joins[i].third.isize( ); j++ )
               if ( touched[ joins[i].third[j] ] ) overlap = True;
          if (overlap) continue;

          vec< triple< vec<int>, vec<int>, vec<int> > > proc;
          vec< vec< vec<int> > > proc1, proc2;
          proc.push_back( joins[i] );
          proc1.push_back( paths1[i] ), proc2.push_back( paths2[i] );

          if ( Inv( joins[i].first[0] ) >= 0 )
          {    vec<int> a, b, c;
               a.push_back( Inv( joins[i].third[0] ), Inv( joins[i].third[1] ) );
               for ( int j = 0; j < joins[i].second.isize( ); j++ )
                    b.push_back( Inv( joins[i].second[j] ) );
               Sort(b);
               c.push_back( Inv( joins[i].first[0] ), Inv( joins[i].first[1] ) );
               vec<int> all1, all2;
               all1.append( joins[i].first );
               all1.append( joins[i].second );
               all1.append( joins[i].third );
               all2.append(a), all2.append(b), all2.append(c);
               Sort(all1), Sort(all2);
               if ( Meet( all1, all2 ) ) continue;
               for ( int j = 0; j < all2.isize( ); j++ )
                    if ( touched[ all2[j] ] ) overlap = True;
               if (overlap) continue;
               proc.push( a, b, c );   
               vec< vec<int> > p1 = paths1[i], p2 = paths2[i];
               for ( int j = 0; j < p1.isize( ); j++ )
               {    p1[j].ReverseMe( );
                    for ( int l = 0; l < p1[j].isize( ); l++ )
                         p1[j][l] = Inv( p1[j][l] );   }
               for ( int j = 0; j < p2.isize( ); j++ )
               {    p2[j].ReverseMe( );
                    for ( int l = 0; l < p2[j].isize( ); l++ )
                         p2[j][l] = Inv( p2[j][l] );   }
               proc1.push_back(p1), proc2.push_back(p2);    }

          vec<fix64_6> weight( EdgeObjectCount( ), 0 );
          for ( int i = 0; i < NPaths( ); i++ )
          for ( int j = 0; j < Path(i).isize( ); j++ )
               weight[ Path(i)[j] ] += Weight(i);
          vec<int> count(2, 0);
          for ( int p = 0; p < proc.isize( ); p++ )
          {    if ( verbosity >= 2 )
               {    int e1 = proc[p].first[0], e2 = proc[p].first[1];
                    int f1 = proc[p].third[0], f2 = proc[p].third[1];
                    cout << "joining: ";
                    PRINT4(e1,e2,f1,f2);    }
               for ( int j = 0; j < proc[p].first.isize( ); j++ )
                    touched[ proc[p].first[j] ] = True;
               for ( int j = 0; j < proc[p].second.isize( ); j++ )
                    touched[ proc[p].second[j] ] = True;
               for ( int j = 0; j < proc[p].third.isize( ); j++ )
                    touched[ proc[p].third[j] ] = True;
               vec<int> dels;
               dels.append( proc[p].first );
               dels.append( proc[p].second );
               dels.append( proc[p].third );
               const int max_del_weight = 4;
               Bool bad = False;
               for ( int i = 0; i < dels.isize( ); i++ )
               {    Bool used = False;
                    for ( int j = 0; j < proc1[p].isize( ); j++ )
                         if ( Member( proc1[p][j], dels[i] ) ) used = True;
                    for ( int j = 0; j < proc2[p].isize( ); j++ )
                         if ( Member( proc2[p][j], dels[i] ) ) used = True;
                    if (used) continue;
                    if ( weight[ dels[i] ] > max_del_weight ) bad = True;    }
               if (bad)
               {    if ( verbosity >= 2 ) cout << "aborting join" << endl;
                    continue;    }
               progress = True;
               int v1 = to_left[ proc[p].first[0] ]; 
               int v2 = to_left[ proc[p].first[1] ];
               int w1 = to_right[ proc[p].third[0] ]; 
               int w2 = to_right[ proc[p].third[1] ];
               vec< vec<int> > e( proc1[p] );
               e.append( proc2[p] );
               vec<int> f;
               for ( int j = 0; j < proc1[p].isize( ); j++ )
               {    f.push_back( EdgeObjectCount( ) );
                    int rid = -1;
                    if ( p == 1 ) rid = EdgeObjectCount( ) - count[0];
                    InvMutable( ).push_back(rid);
                    if ( p == 1 )
                    {    InvMutable( EdgeObjectCount( ) - count[0] ) 
                              = EdgeObjectCount( );    }
                    AddEdge( v1, w1, Cat( proc1[p][j] ) );
                    count[p]++;    }
               for ( int j = 0; j < proc2[p].isize( ); j++ )
               {    f.push_back( EdgeObjectCount( ) );
                    int rid = -1;
                    if ( p == 1 ) rid = EdgeObjectCount( ) - count[0];
                    InvMutable( ).push_back(rid);
                    if ( p == 1 )
                    {    InvMutable( EdgeObjectCount( ) - count[0] ) 
                              = EdgeObjectCount( );    }
                    AddEdge( v2, w2, Cat( proc2[p][j] ) );
                    count[p]++;    }
               DeleteEdges(dels);
               TransformPaths( e, f );    }    }
     if ( verbosity >= 2 ) cout << "\n";

     // Clean up.

     UniqueOrderPaths( );
     RemoveEdgelessVertices( );
     REPORT_TIME( clock, "used in PullApart2" );
     RemoveDeadEdgeObjects( );
     TestValid(logc);
     DeleteReverseComplementComponents(logc);    

     if ( !progress ) break;    }    }
#endif

// Note that we're not yet testing pairs for symmetry, and it's not 100% clear
// what the tests would be.

void SupportedHyperBasevectorPrivate::TestValid( const long_logging& logc,
     const Bool test_paths ) const
{    double clock = WallClockTime( );

     // Test for empty.

     if ( EdgeObjectCount( ) == 0 ) DiscovarTools::ExitAssemblyEmpty( );

     // Test HyperBasevector structure.

     HyperBasevector::TestValid( );
     vec<int> to_left, to_right;
     ToLeft(to_left), ToRight(to_right);

     // Test involution.

     if ( inv_.isize( ) != EdgeObjectCount( ) )
     {    cout << "\n";
          PRINT2( EdgeObjectCount( ), inv_.size( ) );
          cout << "SupportedHyperBasevector is invalid.\n"
               << "Involution has wrong size.\n" << "Abort." << endl;
          TracebackThisProcess( );    }
     for ( int e = 0; e < EdgeObjectCount( ); e++ )
     {    if ( !InvDef(e) ) continue;
          if ( Inv(e) < 0 || Inv(e) >= EdgeObjectCount( ) )
          {    cout << "\n";
               PRINT3( e, Inv(e), EdgeObjectCount( ) );
               cout << "SupportedHyperBasevector is invalid.\n"
                    << "Illegal involution value.\n" << "Abort." << endl;
               TracebackThisProcess( );    }
          basevector b = EdgeObject(e);
          b.ReverseComplement( );
          if ( b != EdgeObject( Inv(e) ) )
          {    cout << "\n";
               int re = Inv(e);
               PRINT4( e, re, b.size( ), EdgeObject(re).size( ) );
               cout << "SupportedHyperBasevector is invalid.\n"
                    << "Involution value not rc.\n" << "Abort." << endl;
               TracebackThisProcess( );    }
          if ( Inv(Inv(e)) != e )
          {    cout << "\nSupportedHyperBasevector is invalid.\n"
                    << "Involution is not an involution.\n" << "Abort." << endl;
               TracebackThisProcess( );    }    }
     for ( int v = 0; v < N( ); v++ )
     {    for ( int i1 = 0; i1 < To(v).isize( ); i1++ )
          for ( int i2 = 0; i2 < From(v).isize( ); i2++ )
          {    int e1 = EdgeObjectIndexByIndexTo( v, i1 );
               int e2 = EdgeObjectIndexByIndexFrom( v, i2 );
               if ( InvDef(e1) != InvDef(e2) )
               {    cout << "\n";
                    PRINT3( e1, EdgeLengthKmers(e1), Inv(e1) );
                    PRINT3( e2, EdgeLengthKmers(e2), Inv(e2) );
                    cout << "SupportedHyperBasevector is invalid.\n"
                         << "Some edges in the same component have the inversion "
                         << "defined and some do not.\n" 
                         << "Abort." << endl;
                    TracebackThisProcess( );    }    }    }
     for ( int v = 0; v < N( ); v++ )
     {    for ( int i1 = 0; i1 < To(v).isize( ); i1++ )
          for ( int i2 = 0; i2 < From(v).isize( ); i2++ )
          {    int e1 = EdgeObjectIndexByIndexTo( v, i1 );
               int e2 = EdgeObjectIndexByIndexFrom( v, i2 );
               int re1 = Inv(e1), re2 = Inv(e2);
               if ( InvDef(e1) && InvDef(e2) && to_right[re2] != to_left[re1] )
               {    cout << "\nSupportedHyperBasevector is invalid.\n"
                         << "Involution does not preserve graph structure.\n" 
                         << "Abort." << endl;
                    TracebackThisProcess( );    }    }    }
     // Note sure the following can happen:
     for ( int v = 0; v < N( ); v++ )
     {    for ( int i = 0; i < From(v).isize( ); i++ )
          {    int w = From(v)[i];
               int e = EdgeObjectIndexByIndexFrom( v, i );
               if ( !InvDef(e) ) continue;
               int re = Inv(e);
               int rv = to_left[re], rw = to_right[re];
               if ( From(rv).size( ) != To(w).size( )
                    || To(rw).size( ) != From(v).size( ) )
               {    cout << "\nSupportedHyperBasevector is invalid.\n"
                         << "Graph structure is asymmetric.\n" << "Abort." << endl;
                    TracebackThisProcess( );    }    }    }    

     // Return if we're not testing paths.

     if ( !test_paths )
     {    REPORT_TIME( clock, "used testing valid" );
          return;    }

     for ( int j = 0; j < NPaths( ); j++ )
     {    if ( Path(j).empty( ) )
          {    cout << "\nSupportedHyperBasevector is invalid.\n"
                    << "Path is empty.\n" << "Abort." << endl;
               TracebackThisProcess( );    }
          for ( int i = 0; i < Path(j).isize( ); i++ )
          {    if ( Path(j,i) < 0 )
               {    cout << "\nSupportedHyperBasevector is invalid.\n"
                         << "Illegal negative path value.\n" << "Abort." << endl;
                    TracebackThisProcess( );    }    }
          for ( int i = 0; i < Path(j).isize( ) - 1; i++ )
          {    if ( Path(j,i) < 0 || Path(j,i+1) < 0 ) continue;
               if ( to_right[ Path(j,i) ] != to_left[ Path(j,i+1) ] )
               {    cout << "\nSupportedHyperBasevector is invalid.\n"
                         << "Non-adjacent edges in path " << j << ": " << Path(j,i) 
                         << ", " << Path(j,i+1) << "\n" << "Abort." << endl;
                    TracebackThisProcess( );    }    }    }    

     if ( Paths( ).size( ) != WeightsFw( ).size( ) )
     {    cout << "\nSupportedHyperBasevector is invalid.\n"
               << "Paths and weights have different sizes.\n" << "Abort." << endl;
          TracebackThisProcess( );    }
     if ( WeightsFw( ).size( ) != WeightsRc( ).size( ) )
     {    cout << "\nSupportedHyperBasevector is invalid.\n"
               << "Fw and rc weights have different sizes.\n" << "Abort." << endl;
          TracebackThisProcess( );    }

     // Test for no support.

     if ( NPaths( ) == 0 )
     {    cout << "\nSupportedHyperBasevector is invalid.\n";
          cout << "There are no paths.\n" << "Abort." << endl;
          DiscovarTools::ExitPathsEmpty( ); }
//          TracebackThisProcess( );    }

     // Test paths for unique ordering.

     if ( !Paths( ).UniqueOrdered( ) )
     {    cout << "\nSupportedHyperBasevector is invalid.\n"
               << "Paths are not unique-ordered.\n" << "Abort." << endl;
               TracebackThisProcess( );    }

     // Test paths and weights for symmetry.

     for ( int i1 = 0; i1 < NPaths( ); i1++ )
     {    const vec<int>& p1 = Path(i1);
          if ( p1[0] < 0 || !InvDef( p1[0] ) ) continue;
          vec<int> p2;
          for ( int j = 0; j < p1.isize( ); j++ )
               p2.push_back( Inv( p1[j] ) );
          p2.ReverseMe( );
          int i2 = BinPosition( Paths( ), p2 );
          if ( i2 < 0 )
          {    cout << "\nPath(" << i1 << ") = " << printSeq( Path(i1) ) << "\n"
                    << "SupportedHyperBasevector is invalid.\n"
                    << "The reverse complement of path " << i1 << " is missing.\n" 
                    << "Abort." << endl;
                    TracebackThisProcess( );    }
          if ( WeightFw(i1) != WeightRc(i2) || WeightRc(i1) != WeightFw(i2) )
          {    cout << "\n";
               PRINT2( WeightFw(i1), WeightRc(i2) );
               PRINT2( WeightRc(i1), WeightFw(i2) );
               cout << "SupportedHyperBasevector is invalid.\n"
                    << "Weights are asymmetric.\n" << "Abort." << endl;
               TracebackThisProcess( );    }    }

     // Test pairs.

     if ( Pairs( ).size( ) != PairData( ).size( ) )
     {    cout << "SupportedHyperBasevector is invalid.\n"
               << "Pairs and PairData have different sizes.\n" << "Abort." << endl;
          TracebackThisProcess( );    }
     if ( !Pairs( ).UniqueOrdered( ) )
     {    cout << "\nSupportedHyperBasevector is invalid.\n"
               << "Pairs are not unique-ordered.\n" << "Abort." << endl;
          TracebackThisProcess( );    }
     for ( int i = 0; i < NPairs( ); i++ )
     {    if ( PairLeft(i).empty( ) || PairRight(i).empty( ) )
          {    cout << "\nSupportedHyperBasevector is invalid.\n"
                    << "Pair side is empty.\n" << "Abort." << endl;
               TracebackThisProcess( );    }    }
     for ( int j = 0; j < NPairs( ); j++ )
     {    for ( int i = 0; i < PairLeft(j).isize( ); i++ )
          {    if ( PairLeft(j,i) < 0 || PairLeft(j,i) >= EdgeObjectCount( ) )
               {    cout << "\nSupportedHyperBasevector is invalid.\n"
                         << "Illegal edge id in left side of pair.\n" 
                         << "Abort." << endl;
                    TracebackThisProcess( );    }    }
          for ( int i = 0; i < PairRight(j).isize( ); i++ )
          {    if ( PairRight(j,i) < 0 || PairRight(j,i) >= EdgeObjectCount( ) )
               {    cout << "\nSupportedHyperBasevector is invalid.\n"
                         << "Illegal edge id in right side of pair.\n" 
                         << "Abort." << endl;
                    TracebackThisProcess( );    }    }
          for ( int i = 0; i < PairLeft(j).isize( ) - 1; i++ )
          {    if ( PairLeft(j,i) < 0 || PairLeft(j,i+1) < 0 ) continue;
               if ( to_right[ PairLeft(j,i) ] != to_left[ PairLeft(j,i+1) ] )
               {    cout << "\nSupportedHyperBasevector is invalid.\n"
                         << "Non-adjacent edges in left side of pair " << j << ": " 
                         << PairLeft(j,i) << ", " << PairLeft(j,i+1) << "\n" 
                         << "Abort." << endl;
                    TracebackThisProcess( );    }    }
          for ( int i = 0; i < PairRight(j).isize( ) - 1; i++ )
          {    if ( PairRight(j,i) < 0 || PairRight(j,i+1) < 0 ) continue;
               if ( to_right[ PairRight(j,i) ] != to_left[ PairRight(j,i+1) ] )
               {    cout << "\nSupportedHyperBasevector is invalid.\n"
                         << "Non-adjacent edges in right side of pair " << j << ": " 
                         << PairRight(j,i) << ", " << PairRight(j,i+1) << "\n" 
                         << "Abort." << endl;
                    TracebackThisProcess( );    }    }
          if ( PairData(j).empty( ) )
          {    cout << "\nSupportedHyperBasevector is invalid.\n"
                    << "No pair data." << endl << "Abort." << endl;
               TracebackThisProcess( );    }
          for ( int i = 0; i < PairData(j).isize( ); i++ )
          {    const pair_point& p = PairData(j,i);
               if ( p.Trim( ) < -1000000000 || p.Trim( ) > 1000000000 )
               {    cout << "\nSupportedHyperBasevector is invalid.\n"
                         << "Insane value " << p.Trim( ) << " for trim." 
                         << endl << "Abort." << endl;
                    TracebackThisProcess( );    }
               if ( p.Lib( ) < 0 )
               {    cout << "\nSupportedHyperBasevector is invalid.\n"
                         << "Negative library id." << endl << "Abort." << endl;
                    TracebackThisProcess( );    }
               if ( p.Weight( ) < 0 )
               {    cout << "\nSupportedHyperBasevector is invalid.\n"
                         << "Pair has negative weight." << endl << "Abort." << endl;
                    TracebackThisProcess( );    }    }    }

     REPORT_TIME( clock, "used testing valid" );    }

void SupportedHyperBasevector::DeleteReverseComplementComponents( const long_logging& logc , const int64_t iDirSortFactor)
{    double clock = WallClockTime( );
     vec< vec<int> > components;
     ComponentEdges(components);
     
     if(iDirSortFactor==0){
         Sort(components);
     }
     else{
         vec<int> to_left, to_right;
         ToLeft(to_left), ToRight(to_right);
         vec<std::tuple<int64_t,int64_t,vec<int>>> order;
         order.reserve(components.size());
         
         for(const auto& entry: components){
             int64_t diff=0,bdiff=0;
             for(const auto& edge: entry){
                 if( Source(to_left[edge])){
                     bdiff-=iDirSortFactor*EdgeObject(edge).size();
                     diff-=iDirSortFactor;
                 }
                 if( Sink(to_right[edge])){
                     bdiff+=iDirSortFactor*EdgeObject(edge).size();
                     diff+=iDirSortFactor;
                 }
             }
             order.emplace_back(diff,bdiff,entry);
         }
         Sort(order);
         
         components.clear();
         for(const auto&entry: order){ components.push_back(std::get<2>(entry)); }
     }
     
     vec<int> rc_to_delete;
     for ( size_t i = 0; i < components.size(); i++ )
     {    vec<int> rc;
          for ( size_t j = 0; j < components[i].size(); j++ )
               rc.push_back( Inv( components[i][j] ) );
          Sort(rc);
          int p = (iDirSortFactor==0)?BinPosition( components, rc ):Position( components, rc );
          if ( p > (int)i )
          {    for ( size_t j = 0; j < components[p].size(); j++ )
                    rc_to_delete.push_back( components[p][j] );
               for ( size_t j = 0; j < components[i].size(); j++ )
                    InvMutable( components[i][j] ) = -1;    }    }
     Sort(rc_to_delete);
     if ( logc.verb[ "DELETE_RC" ] >= 1 )
     {    cout << Date( ) << ": rc component deletion -- deleting " 
               << rc_to_delete.size( ) << " edges" << endl;    }
     DeleteEdges(rc_to_delete);
     DeleteUnusedPaths( );
     REPORT_TIME( clock, "used deleting reverse complement components" );
     double eclock = WallClockTime( );
     RemoveEdgelessVertices( );
     RemoveUnneededVertices( );
     REPORT_TIME( eclock, "used cleaning up after deleting rc components" );
     RemoveDeadEdgeObjects( );
     TestValid(logc);    }

void SupportedHyperBasevector::RemoveSmallComponents( const int K )
{    HyperBasevector::RemoveSmallComponents(K);
     DeleteUnusedPaths( );    }

void SupportedHyperBasevector::Reverse( )
{    HyperBasevector::Reverse( );
     for ( int i = 0; i < NPaths( ); i++ )
          PathMutable(i).ReverseMe( );
     for ( int i = 0; i < NPairs( ); i++ )
     {    PairLeftMutable(i).ReverseMe( );
          PairRightMutable(i).ReverseMe( );
          PairMutable(i) = make_pair( PairRight(i), PairLeft(i) );    }    }

void SupportedHyperBasevector::ChunkEdges( const long_logging& logc )
{    
     // Define heuristics.

     double clock1 = WallClockTime( );
     cout << Date( ) << ": chunking edges" << endl;
     const double min_chunk = 5.0;

     // Set up indices.

     vec< vec< pair<int,int> > > paths_index( EdgeObjectCount( ) );
     for ( int i = 0; i < Paths( ).isize( ); i++ )
     for ( int j = 0; j < Path(i).isize( ); j++ )
          if ( Path(i,j) >= 0 ) paths_index[ Path(i,j) ].push( i, j );
     vec< vec< pair<int,int> > > lefts_index( EdgeObjectCount( ) );
     for ( int i = 0; i < Pairs( ).isize( ); i++ )
     for ( int j = 0; j < PairLeft(i).isize( ); j++ )
          if ( PairLeft(i,j) >= 0 ) lefts_index[ PairLeft(i,j) ].push( i, j );
     vec< vec< pair<int,int> > > rights_index( EdgeObjectCount( ) );
     for ( int i = 0; i < Pairs( ).isize( ); i++ )
     for ( int j = 0; j < PairRight(i).isize( ); j++ )
          if ( PairRight(i,j) >= 0 ) rights_index[ PairRight(i,j) ].push( i, j );
     vec<int> to_left, to_right;
     ToLeft(to_left), ToRight(to_right);

     // Go through the edges.

     vec< vec<int> > chunks;
     for ( int x1 = 0; x1 < EdgeObjectCount( ); x1++ )
     {    
          // Find eligible extensions of x1, and take the longest one.

          vec< vec<int> > exts;
          vec<fix64_6> extsw;
          for ( int i = 0; i < paths_index[x1].isize( ); i++ )
          {    int id = paths_index[x1][i].first, pos = paths_index[x1][i].second;
               if ( pos == Path(id).isize( ) - 1 ) continue;
               vec<int> v;
               int nkmers = 0;
               for ( int j = pos; j < Path(id).isize( ); j++ )
               {    if ( j-1 > pos ) nkmers += EdgeLengthKmers( Path(id,j-1) );
                    if ( nkmers > MedianCorrectedReadLengthFudge( ) ) break;
                    v.push_back( Path(id,j) );
                    if ( j > pos )
                    {    exts.push_back(v); 
                         extsw.push_back( Weight(id) );    }    }    }
          SortSync( exts, extsw );
          vec<Bool> to_delete( exts.size( ), False );
          for ( int i = 0; i < exts.isize( ); i++ )
          {    int j = exts.NextDiff(i);
               for ( int k = i + 1; k < j; k++ )
               {    extsw[i] += extsw[k];
                    to_delete[k] = True;    }
               i = j - 1;    }
          EraseIf( exts, to_delete ), EraseIf( extsw, to_delete );
          vec<int> len;
          for ( int i = 0; i < exts.isize( ); i++ )
               len.push_back( exts[i].size( ) );
          SortSync( len, exts, extsw );
          int fb;
          for ( fb = 0; fb < exts.isize( ) - 1; fb++ )
               if ( len[fb] == len[fb+1] || extsw[fb] < min_chunk ) break;
          exts.resize(fb);
          vec<Bool> to_delete2( exts.size( ), False );
          for ( int i = 0; i < exts.isize( ); i++ )
          {    int xn = exts[i].back( );
               for ( int j = 0; j < paths_index[xn].isize( ); j++ )
               {    int id = paths_index[xn][j].first; 
                    int pos = paths_index[xn][j].second;
                    Bool incomp = False;
                    for ( int k = pos - 1; k >= 0; k-- )
                    {    int e = exts[i].isize( ) - 1 - (pos-k);
                         if ( e < 0 ) break;
                         if ( exts[i][e] != Path(id,k) )
                         {    incomp = True;
                              break;    }    }
                    if (incomp) 
                    {    to_delete2[i] = True;
                         break;    }    }    }
          EraseIf( exts, to_delete2 );
          if ( exts.nonempty( ) ) chunks.push_back( exts.back( ) );    }

     // Keep the largest chunks that don't overlap other chunks.  Also we require
     // symmetry with respect to the involution.

     vec<int> chunklen;
     for ( int i = 0; i < chunks.isize( ); i++ )
          chunklen.push_back( chunks[i].size( ) );
     ReverseSortSync( chunklen, chunks );
     vec<Bool> to_delete( chunks.size( ), False ), used( EdgeObjectCount( ), False );
     for ( int i = 0; i < chunks.isize( ); i++ )
     {    Bool conflict = False;
          for ( int j = 0; j < chunks[i].isize( ); j++ )
               if ( used[ chunks[i][j] ] ) conflict = True;
          if (conflict)
          {    to_delete[i] = True;
               continue;    }
          for ( int j = 0; j < chunks[i].isize( ); j++ )
               used[ chunks[i][j] ] = True;    }
     EraseIf( chunks, to_delete );
     Sort(chunks);
     vec<Bool> to_delete2( chunks.size( ), False );
     for ( int i = 0; i < chunks.isize( ); i++ )
     {    if ( !InvDef( chunks[i][0] ) ) continue;
          vec<int> x;
          for ( int j = 0; j < chunks[i].isize( ); j++ )
               x.push_back( Inv( chunks[i][j] ) );
          x.ReverseMe( );
          if ( !BinMember( chunks, x ) ) to_delete2[i] = True;    }
     EraseIf( chunks, to_delete2 );

     // Define involution of chunks.

     vec<int> inv2( chunks.size( ), -1 );
     for ( int i = 0; i < chunks.isize( ); i++ )
     {    if ( !InvDef( chunks[i][0] ) ) continue;
          vec<int> x;
          for ( int j = 0; j < chunks[i].isize( ); j++ )
               x.push_back( Inv( chunks[i][j] ) );
          x.ReverseMe( );
          inv2[i] = BinPosition( chunks, x );    }

     // Process the chunks.

     Bool verbose = ( logc.verb[ "CHUNK_EDGES" ] >= 1 );
     int nedges = EdgeObjectCount( );
     vec<Bool> deleted( nedges, False );
     for ( int ci = 0; ci < chunks.isize( ); ci++ )
     {    const vec<int>& v = chunks[ci];
          if (verbose) 
          {    cout << "processing " << printSeq(v);
               if ( InvDef( v[0] ) )
               {    vec<int> rv;
                    for ( int j = 0; j < v.isize( ); j++ )
                         rv.push_back( Inv(v[j]) );
                    rv.ReverseMe( );
                    cout << " (inv = " << printSeq(rv) << ")";    }
               cout << endl;    }

          // Edit.

          int x1 = v.front( ), xn = v.back( );
          AddEdge( to_left[x1], to_right[xn], Cat(v) );
          to_left.push_back( to_left[x1] ), to_right.push_back( to_right[xn] );
          if ( inv2[ci] < 0 ) InvMutable( ).push_back(-1);
          else InvMutable( ).push_back( inv2[ci] + nedges );
          vec<int> dels, edited, edited_left, edited_right;
          dels.push_back( x1, xn );

          // Apply condition (C).

          for ( int i = 0; i < paths_index[x1].isize( ); i++ )
          {    int id = paths_index[x1][i].first, pos1 = paths_index[x1][i].second;
               if ( Path(id,pos1) != x1 ) continue;
               if (verbose) cout << id << ": 1 replacing " << printSeq( Path(id) );
               int pos2;
               for ( pos2 = pos1 + 1; pos2 < Path(id).isize( ); pos2++ )
               {    if ( Path(id,pos2) == xn )
                    {    pos2++;
                         break;    }    }
               PathMutable(id)[pos1] = EdgeObjectCount( ) - 1;
               for ( int j = pos1 + 1; j < pos2; j++ )
                    PathMutable(id)[j] = -1;
               if (verbose) cout << " by " << printSeq( Path(id) ) << endl;
               edited.push_back(id);    }
          for ( int i = 0; i < paths_index[xn].isize( ); i++ )
          {    int id = paths_index[xn][i].first, pos1 = paths_index[xn][i].second;
               if ( Path(id,pos1) != xn ) continue;
               if (verbose) cout << id << ": 2 replacing " << printSeq( Path(id) );
               int pos2;
               for ( pos2 = pos1 - 1; pos2 >= 0; pos2-- )
               {    if ( Path(id,pos2) == x1 )
                    {    pos2--;
                         break;    }    }
               PathMutable(id)[pos1] = EdgeObjectCount( ) - 1;
               for ( int j = pos1 - 1; j > pos2; j-- )
                    PathMutable(id)[j] = -1;
               if (verbose) cout << " by " << printSeq( Path(id) ) << endl;
               edited.push_back(id);    }
          for ( int i = 0; i < paths_index[xn].isize( ); i++ )
          {    int id = paths_index[xn][i].first, pos1 = paths_index[xn][i].second;
               if ( Path(id,pos1) != xn ) continue;
               int pos2;
               for ( pos2 = pos1 - 1; pos2 >= 0; pos2-- )
               {    if ( Path(id,pos2) == x1 )
                    {    pos2--;
                         break;    }    }
               PathMutable(id)[pos1] = EdgeObjectCount( ) - 1;
               for ( int j = pos1 - 1; j > pos2; j-- )
                    PathMutable(id)[j] = -1;
               edited.push_back(id);    }
          for ( int i = 0; i < lefts_index[x1].isize( ); i++ )
          {    int id = lefts_index[x1][i].first, pos1 = lefts_index[x1][i].second;
               if ( PairLeft(id,pos1) != x1 ) continue;
               int pos2;
               for ( pos2 = pos1 + 1; pos2 < PairLeft(id).isize( ); pos2++ )
               {    if ( PairLeft(id,pos2) == xn )
                    {    pos2++;
                         break;    }    }
               PairLeftMutable(id)[pos1] = EdgeObjectCount( ) - 1;
               for ( int j = pos1 + 1; j < pos2; j++ )
                    PairLeftMutable(id)[j] = -1;
               int add = 0;
               for ( int j = pos2; j < pos1 + v.isize( ); j++ )
                    add += EdgeLengthKmers( v[ j - pos1 ] );
               AddTrim( id, add );
               edited_left.push_back(id);    }

          for ( int i = 0; i < lefts_index[xn].isize( ); i++ )
          {    int id = lefts_index[xn][i].first, pos1 = lefts_index[xn][i].second;
               if ( PairLeft(id,pos1) != xn ) continue;
               int pos2;
               for ( pos2 = pos1 - 1; pos2 >= 0; pos2-- )
               {    if ( PairLeft(id,pos2) == x1 )
                    {    pos2--;
                         break;    }    }
               PairLeftMutable(id)[pos1] = EdgeObjectCount( ) - 1;
               for ( int j = pos1 - 1; j > pos2; j-- )
                    PairLeftMutable(id)[j] = -1;
               edited_left.push_back(id);    }

          for ( int i = 0; i < rights_index[x1].isize( ); i++ )
          {    int id = rights_index[x1][i].first, pos1 = rights_index[x1][i].second;
               if ( PairRight(id,pos1) != x1 ) continue;
               int pos2;
               for ( pos2 = pos1 + 1; pos2 < PairRight(id).isize( ); pos2++ )
               {    if ( PairRight(id,pos2) == xn )
                    {    pos2++;
                         break;    }    }
               PairRightMutable(id)[pos1] = EdgeObjectCount( ) - 1;
               for ( int j = pos1 + 1; j < pos2; j++ )
                    PairRightMutable(id)[j] = -1;
               edited_right.push_back(id);    }

          for ( int i = 0; i < rights_index[xn].isize( ); i++ )
          {    int id = rights_index[xn][i].first, pos1 = rights_index[xn][i].second;
               if ( PairRight(id,pos1) != xn ) continue;
               int pos2;
               for ( pos2 = pos1 - 1; pos2 >= 0; pos2-- )
               {    if ( PairRight(id,pos2) == x1 )
                    {    pos2--;
                         break;    }    }
               PairRightMutable(id)[pos1] = EdgeObjectCount( ) - 1;
               for ( int j = pos1 - 1; j > pos2; j-- )
                    PairRightMutable(id)[j] = -1;
               int add = 0;
               for ( int j = pos1 - v.isize( ) + 1; j < 0; j++ )
                    add += EdgeLengthKmers( v[-j-1] );
               AddTrim( id, add );
               edited_right.push_back(id);    }

          // Apply condition (D).

          for ( int r = 1; r < v.isize( ) - 1; r++ )
          {    int xi = v[r];
               if (verbose) cout << "looking at xi = " << xi << endl;
               Bool compat = True;
               for ( int i = 0; i < paths_index[xi].isize( ); i++ )
               {    int id = paths_index[xi][i].first; 
                    int pos = paths_index[xi][i].second;
                    if ( Path(id,pos) != xi ) continue;
                    Bool compati = Compati( Path(id), v, r, pos );
                    if ( !compati ) 
                    {    compat = False;
                         break;    }    }
               for ( int i = 0; i < lefts_index[xi].isize( ); i++ )
               {    int id = lefts_index[xi][i].first; 
                    int pos = lefts_index[xi][i].second;
                    if ( PairLeft(id,pos) != xi ) continue;
                    Bool compati = Compati( PairLeft(id), v, r, pos );
                    if ( !compati ) 
                    {    compat = False;
                         break;    }    }
               for ( int i = 0; i < rights_index[xi].isize( ); i++ )
               {    int id = rights_index[xi][i].first; 
                    int pos = rights_index[xi][i].second;
                    if ( PairRight(id,pos) != xi ) continue;
                    Bool compati = Compati( PairRight(id), v, r, pos );
                    if ( !compati ) 
                    {    compat = False;
                         break;    }    }
               if ( !compat ) continue;

               for ( int i = 0; i < paths_index[xi].isize( ); i++ )
               {    int id = paths_index[xi][i].first; 
                    int pos = paths_index[xi][i].second; 
                    if ( Path(id,pos) != xi ) continue;
                    if (verbose) 
                         cout << id << ": 3 replacing " << printSeq( Path(id) );
                    PathMutable(id,0) = EdgeObjectCount( ) - 1;
                    for ( int j = 1; j < Path(id).isize( ); j++ )
                         PathMutable(id,j) = -1;
                    if (verbose) cout << " by " << printSeq( Path(id) ) << endl;
                    edited.push_back(id);    }

               for ( int i = 0; i < lefts_index[xi].isize( ); i++ )
               {    int id = lefts_index[xi][i].first; 
                    int pos = lefts_index[xi][i].second; 
                    if ( PairLeft(id,pos) != xi ) continue;
                    PairLeftMutable(id,0) = EdgeObjectCount( ) - 1;
                    for ( int j = 1; j < PairLeft(id).isize( ); j++ )
                         PairLeftMutable(id,j) = -1;
                    int start = Position( v, PairLeft(id) ), add = 0;
                    ForceAssertGe( start, 0 );
                    int stop = start + v.isize( );
                    for ( int j = stop; j < v.isize( ) - 1; j++ )
                         add += EdgeLengthKmers( v[j] );
                    AddTrim( id, add );
                    edited_left.push_back(id);    }

               for ( int i = 0; i < rights_index[xi].isize( ); i++ )
               {    int id = rights_index[xi][i].first; 
                    int pos = rights_index[xi][i].second; 
                    if ( PairRight(id,pos) != xi ) continue;
                    PairRightMutable(id,0) = EdgeObjectCount( ) - 1;
                    for ( int j = 1; j < PairRight(id).isize( ); j++ )
                         PairRightMutable(id,j) = -1;
                    int start = Position( v, PairRight(id) ), add = 0;
                    ForceAssertGe( start, 0 );
                    for ( int j = 0; j < start; j++ )
                         add += EdgeLengthKmers( v[j] );
                    AddTrim( id, add );
                    edited_right.push_back(id);    }

               dels.push_back(xi);    }

          // Apply condition (E).

          for ( int i = 1; i < v.isize( ) - 2; i++ )
          for ( int j = i+1; j < v.isize( ); j++ )
          {    int xi = v[i], xj = v[j];
               vec<int> w;
               for ( int m = i; m <= j; m++ )
                    w.push_back( v[m] );

               // First find all positions of xi,...,xj within x2,...,xn-1.

               vec<int> starts;
               for ( int l = 1; l < v.isize( ) - (j-i); l++ )
                    if ( v.Contains( w, l ) ) starts.push_back(l);
               ForceAssertGe( starts.isize( ), 1 ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXX

               // Now find all paths containing xi,...,xj, and test for 
               // compatibility.

               Bool compat = True;

               for ( int l = 0; l < paths_index[xi].isize( ); l++ )
               {    int id = paths_index[xi][l].first; 
                    int pos = paths_index[xi][l].second;
                    if ( !Path(id).Contains( w, pos ) ) continue;
                    Bool compatl = False;
                    for ( int m = 0; m < starts.isize( ); m++ )
                    {    if ( w.Contains( Path(id), starts[m] ) )
                         {    compatl = True;
                              break;    }    }
                    if ( !compatl )
                    {    compat = False;
                         break;    }    }

               for ( int l = 0; l < lefts_index[xi].isize( ); l++ )
               {    int id = lefts_index[xi][l].first; 
                    int pos = lefts_index[xi][l].second;
                    if ( !PairLeft(id).Contains( w, pos ) ) continue;
                    Bool compatl = False;
                    for ( int m = 0; m < starts.isize( ); m++ )
                    {    if ( w.Contains( PairLeft(id), starts[m] ) )
                         {    compatl = True;
                              break;    }    }
                    if ( !compatl )
                    {    compat = False;
                         break;    }    }

               for ( int l = 0; l < rights_index[xi].isize( ); l++ )
               {    int id = rights_index[xi][l].first; 
                    int pos = rights_index[xi][l].second;
                    if ( !PairRight(id).Contains( w, pos ) ) continue;
                    Bool compatl = False;
                    for ( int m = 0; m < starts.isize( ); m++ )
                    {    if ( w.Contains( PairRight(id), starts[m] ) )
                         {    compatl = True;
                              break;    }    }
                    if ( !compatl )
                    {    compat = False;
                         break;    }    }

               if ( !compat ) continue;

               for ( int l = 0; l < paths_index[xi].isize( ); l++ )
               {    int id = paths_index[xi][l].first; 
                    int pos = paths_index[xi][l].second;
                    if ( !Path(id).Contains( w, pos ) ) continue;
                    if (verbose) 
                         cout << id << ": 4 replacing " << printSeq( Path(id) );
                    PathMutable(id,0) = EdgeObjectCount( ) - 1;
                    for ( int j = 1; j < Path(id).isize( ); j++ )
                         PathMutable(id,j) = -1;    
                    if (verbose) cout << " by " << printSeq( Path(id) ) << endl;
                    edited.push_back(id);    }

               for ( int l = 0; l < lefts_index[xi].isize( ); l++ )
               {    int id = lefts_index[xi][l].first; 
                    int pos = lefts_index[xi][l].second;
                    if ( !PairLeft(id).Contains( w, pos ) ) continue;
                    PairLeftMutable(id,0) = EdgeObjectCount( ) - 1;
                    for ( int m = 1; m < PairLeft(id).isize( ); m++ )
                         PairLeftMutable(id,m) = -1;    
                    int add = 0;
                    for ( int k = j + 1; k < v.isize( ); k++ )
                         add += EdgeLengthKmers( v[k] );
                    AddTrim( id, add );
                    edited_left.push_back(id);    }

               for ( int l = 0; l < rights_index[xi].isize( ); l++ )
               {    int id = rights_index[xi][l].first; 
                    int pos = rights_index[xi][l].second;
                    if ( !PairRight(id).Contains( w, pos ) ) continue;
                    PairRightMutable(id,0) = EdgeObjectCount( ) - 1;
                    for ( int m = 1; m < PairRight(id).isize( ); m++ )
                         PairRightMutable(id,m) = -1;    
                    int add = 0;
                    for ( int k = 0; k < i; k++ )
                         add += EdgeLengthKmers( v[k] );
                    AddTrim( id, add );
                    edited_right.push_back(id);    }    }

          // Finish up.

          DeleteEdges(dels);
          for ( int j = 0; j < dels.isize( ); j++ )
               deleted[ dels[j] ] = True;

          UniqueSort(edited);
          for ( int j = 0; j < edited.isize( ); j++ )
               RemoveNegatives( PathMutable( edited[j] ) );
          for ( int i = 0; i < v.isize( ); i++ )
          {    vec<Bool> to_delete( paths_index[ v[i] ].size( ), False );
               for ( int j = 0; j < paths_index[ v[i] ].isize( ); j++ )
               {    if ( BinMember( edited, paths_index[ v[i] ][j].first ) )
                         to_delete[j] = True;    }
               EraseIf( paths_index[ v[i] ], to_delete );    }

          UniqueSort(edited_left);
          for ( int j = 0; j < edited_left.isize( ); j++ )
               RemoveNegatives( PairLeftMutable( edited_left[j] ) );
          for ( int i = 0; i < v.isize( ); i++ )
          {    vec<Bool> to_delete( lefts_index[ v[i] ].size( ), False );
               for ( int j = 0; j < lefts_index[ v[i] ].isize( ); j++ )
               {    if ( BinMember( edited_left, lefts_index[ v[i] ][j].first ) )
                         to_delete[j] = True;    }
               EraseIf( lefts_index[ v[i] ], to_delete );    }

          UniqueSort(edited_right);
          for ( int j = 0; j < edited_right.isize( ); j++ )
               RemoveNegatives( PairRightMutable( edited_right[j] ) );
          for ( int i = 0; i < v.isize( ); i++ )
          {    vec<Bool> to_delete( rights_index[ v[i] ].size( ), False );
               for ( int j = 0; j < rights_index[ v[i] ].isize( ); j++ )
               {    if ( BinMember( edited_right, rights_index[ v[i] ][j].first ) )
                         to_delete[j] = True;    }
               EraseIf( rights_index[ v[i] ], to_delete );    }

          vec<int> vps(v);
          vps.push_back( EdgeObjectCount( ) - 1 );
          UniqueSort(vps);
          paths_index.resize( EdgeObjectCount( ) );
          lefts_index.resize( EdgeObjectCount( ) );
          rights_index.resize( EdgeObjectCount( ) );
          vec<int> touched(vps);

          for ( int j = 0; j < edited.isize( ); j++ )
          {    int id = edited[j];
               const vec<int>& p = Path(id);
               for ( int l = 0; l < p.isize( ); l++ )
               {    int y = p[l];
                    if ( BinMember( vps, y ) ) continue;
                    touched.push_back(y);
                    vec<Bool> to_delete( paths_index[y].isize( ), False );
                    for ( int m = 0; m < paths_index[y].isize( ); m++ )
                    {    if ( paths_index[y][m].first == id )
                              to_delete[m] = True;    }
                    EraseIf( paths_index[y], to_delete );    }
               for ( int l = 0; l < p.isize( ); l++ )
                    paths_index[ p[l] ].push( id, l );    }

          for ( int j = 0; j < edited_left.isize( ); j++ )
          {    int id = edited_left[j];
               const vec<int>& p = PairLeft(id);
               for ( int l = 0; l < p.isize( ); l++ )
               {    int y = p[l];
                    if ( BinMember( vps, y ) ) continue;
                    touched.push_back(y);
                    vec<Bool> to_delete( lefts_index[y].isize( ), False );
                    for ( int m = 0; m < lefts_index[y].isize( ); m++ )
                    {    if ( lefts_index[y][m].first == id )
                              to_delete[m] = True;    }
                    EraseIf( lefts_index[y], to_delete );    }
               for ( int l = 0; l < p.isize( ); l++ )
                    lefts_index[ p[l] ].push( id, l );    }

          for ( int j = 0; j < edited_right.isize( ); j++ )
          {    int id = edited_right[j];
               const vec<int>& p = PairRight(id);
               for ( int l = 0; l < p.isize( ); l++ )
               {    int y = p[l];
                    if ( BinMember( vps, y ) ) continue;
                    touched.push_back(y);
                    vec<Bool> to_delete( rights_index[y].isize( ), False );
                    for ( int m = 0; m < rights_index[y].isize( ); m++ )
                    {    if ( rights_index[y][m].first == id )
                              to_delete[m] = True;    }
                    EraseIf( rights_index[y], to_delete );    }
               for ( int l = 0; l < p.isize( ); l++ )
                    rights_index[ p[l] ].push( id, l );    }

          UniqueSort(touched);
          for ( int i = 0; i < touched.isize( ); i++ )
          {    Sort( paths_index[ touched[i] ] );
               Sort( lefts_index[ touched[i] ] );
               Sort( rights_index[ touched[i] ] );    }
          for ( int i = 0; i < dels.isize( ); i++ )
          {    to_left[ dels[i] ] = -1, to_right[ dels[i] ] = -1;    }    }

     // Clean up.

     RemoveEdgelessVertices( );
     RemoveUnneededVertices( );
     REPORT_TIME( clock1, "used in ChunkEdges 1" );
     RemoveDeadEdgeObjects( );
     double clock2 = WallClockTime( );
     UniqueOrderPaths( );
     cout << Date( ) << ": after chunking there are "
          << EdgeObjectCount( ) << " edges" << endl;
     REPORT_TIME( clock2, "used in ChunkEdges 2" );
     TestValid(logc);    }
