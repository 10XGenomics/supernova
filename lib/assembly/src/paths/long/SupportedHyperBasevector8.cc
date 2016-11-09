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
#include "paths/long/SupportedHyperBasevector.h"

// Hookup.
// Design:
// Look for edge sequences e = e1,...,en such that enough read paths contain the
// sequence, such that all paths containing e1 are consistent with, and likewise
// with en.  Remove e1 and en from the graph and substitute the concatenation of e.
// The actual behavior is more complicated.

void SupportedHyperBasevector::Hookup( const long_logging& logc )
{
     int verbosity = logc.verb[ "HOOKUP" ];

     if ( verbosity >= 1 )
     {    cout << EdgeObjectCount( ) << " edges and " << NComponents( )
               << " components" << endl;
          vec<int> len;
          for ( int e = 0; e < EdgeObjectCount( ); e++ )
               len.push_back( EdgeLengthKmers(e) );
          Sort(len);
          cout << "N50 edge = " << N50(len) << endl;    }

     // Define heuristics.

     const fix64_6 min_wt = 5;

     // Run two passes.

     for ( int pass = 1; pass <= 2; pass++ )
     {
          // Compute unique right and left extensions.
     
          int ne = EdgeObjectCount( );
          vec< vec<int> > right(ne), left(ne);
          for ( int i = 0; i < NPaths( ); i++ )
          {    const vec<int>& p = Path(i);
               for ( int j = 0; j < p.isize( ) - 1; j++ )
               {    int e = p[j];
                    if ( !( right[e].solo( ) && right[e][0] < 0 ) )
                    {    Bool conflict = False;
                         for ( int k = j + 1; k < p.isize( ); k++ )
                         {    if ( k - j - 1 >= right[e].isize( ) ) break;
                              if ( p[k] != right[e][ k - j - 1 ] )
                              {    conflict = True;
                                   break;    }    }
                         if (conflict)
                         {    right[e].resize(1);
                              right[e][0] = -1;    }
                         else if ( p.isize( ) - j - 1 > right[e].isize( ) )
                         {    right[e].reserve( p.isize( ) - j - 1 );
                              for ( int k = right[e].isize( ) + j + 1; k < p.isize( ); 
                                   k++ )
                              {    right[e].push_back( p[k] );    }    }    }    }
               for ( int j = 1; j < p.isize( ); j++ )
               {    int e = p[j];
                    if ( !( left[e].solo( ) && left[e][0] < 0 ) )
                    {    Bool conflict = False;
                         for ( int k = j - 1; k >= 0; k-- )
                         {    if ( (j-1) - k >= left[e].isize( ) ) break;
                              if ( p[k] != left[e][ (j-1) - k ] )
                              {    conflict = True;
                                   break;    }    }
                         if (conflict)
                         {    left[e].resize(1);
                              left[e][0] = -1;    }
                         else if ( j > left[e].isize( ) )
                         {    left[e].reserve(j);
                              for ( int k = j - 1 - left[e].isize( ); k >= 0; k-- )
                                   left[e].push_back( p[k] );    }    }    }    }

          // Truncate so that ends go back.

          for ( int e = 0; e < ne; e++ )
          {    while(1)
               {    Bool progress = False;
                    if ( right[e].nonempty( ) && right[e][0] >= 0 )
                    {    int f = right[e].back( );
                         if ( left[f].empty( ) || left[f][0] < 0 )
                         {    right[e].pop_back( );
                              progress = True;    }    }
                    if ( left[e].nonempty( ) && left[e][0] >= 0 )
                    {    int f = left[e].back( );
                         if ( right[f].empty( ) || right[f][0] < 0 )
                         {    left[e].pop_back( );
                              progress = True;    }    }
                    if ( !progress ) break;    }    }
     
          // Truncate until multiplicities high enough.
     
          vec< vec<fix64_6> > rmult(ne), lmult(ne);
          for ( int e = 0; e < ne; e++ )
          {    if ( right[e].nonempty( ) && right[e][0] >= 0 )
                    rmult[e].resize( right[e].size( ), 0 );
               if ( left[e].nonempty( ) && left[e][0] >= 0 )
                    lmult[e].resize( left[e].size( ), 0 );    }
          for ( int i = 0; i < NPaths( ); i++ )
          {    const vec<int>& p = Path(i);
               for ( int j = 0; j < p.isize( ); j++ )
               {    int e = p[j];
                    if ( right[e].nonempty( ) && right[e][0] >= 0 )
                    {    int rext = Min( right[e].isize( ), p.isize( ) - j - 1 );
                         for ( int k = 0; k < rext; k++ )
                              rmult[e][k] += Weight(i);    }
                    if ( left[e].nonempty( ) && left[e][0] >= 0 )
                    {    int lext = Min( left[e].isize( ), j );
                         for ( int k = 0; k < lext; k++ )
                              lmult[e][k] += Weight(i);    }    }    }
          for ( int e = 0; e < ne; e++ )
          {    while ( right[e].nonempty( ) && right[e][0] >= 0 
                    && rmult[e].back( ) < min_wt )
               {    right[e].pop_back( ), rmult[e].pop_back( );    }
               while ( left[e].nonempty( ) && left[e][0] >= 0 
                    && lmult[e].back( ) < min_wt )
               {    left[e].pop_back( ), lmult[e].pop_back( );    }    }

          // Truncate so that ends go back.

          for ( int e = 0; e < ne; e++ )
          {    while(1)
               {    Bool progress = False;
                    if ( right[e].nonempty( ) && right[e][0] >= 0 )
                    {    int f = right[e].back( );
                         if ( left[f].empty( ) || left[f][0] < 0 )
                         {    right[e].pop_back( );
                              progress = True;    }    }
                    if ( left[e].nonempty( ) && left[e][0] >= 0 )
                    {    int f = left[e].back( );
                         if ( right[f].empty( ) || right[f][0] < 0 )
                         {    left[e].pop_back( );
                              progress = True;    }    }
                    if ( !progress ) break;    }    }
     
          // Find reciprocal extensions having high enough weight.
     
          int count = 0;
          vec< vec<int> > joins;
          for ( int e = 0; e < ne; e++ )
          {    if ( right[e].empty( ) || right[e][0] < 0 ) continue;
               int f = right[e].back( );
               if ( left[f].empty( ) || left[f][0] < 0 ) continue;
               vec<int> x = {e};
               x.append( right[e] );
               joins.push_back(x);
               // cout << "[" << count++ << "] " << printSeq(x) << " (" 
               //      << EdgeLengthKmers(e);
               // for ( int j = 0; j < right[e].isize( ); j++ )
               //      cout << "," << EdgeLengthKmers( right[e][j] );
               // cout << ")" << endl;    
                    }
          for ( int e = 0; e < ne; e++ )
          {    if ( left[e].empty( ) || left[e][0] < 0 ) continue;
               int f = left[e].back( );
               if ( right[f].empty( ) || right[f][0] < 0 ) continue;
               vec<int> x = {e};
               x.append( left[e] );
               x.ReverseMe( );
               joins.push_back(x);    }
          UniqueSort(joins);

          // Find paths through joins.

          digraphE< vec<int> > G;
          vec< vec<int> > paths;
          {    vec<vec<int>> from(ne), to(ne), from_edge_obj(ne), to_edge_obj(ne);
               vec<vec<int>> edges;
               for ( int i = 0; i < joins.isize( ); i++ )
               {    const vec<int>& x = joins[i];
                    from[ x.front( ) ].push_back( x.back( ) );
                    from_edge_obj[ x.front( ) ].push_back( edges.size( ) );
                    to[ x.back( ) ].push_back( x.front( ) );
                    to_edge_obj[ x.back( ) ].push_back( edges.size( ) );
                    vec<int> e;
                    for ( int j = 1; j < x.isize( ) - 1; j++ )
                         e.push_back( x[j] );
                    edges.push_back(e);    }
               for ( int e = 0; e < ne; e++ )
               {    SortSync( from[e], from_edge_obj[e] );
                    SortSync( to[e], to_edge_obj[e] );    }
               G.Initialize( from, to, edges, to_edge_obj, from_edge_obj );
               G.AllPaths( -1, -1, paths );
               vec<Bool> to_delete( paths.size( ), False );
               for ( int i = 0; i < paths.isize( ); i++ )
                    if ( paths[i].solo( ) ) to_delete[i] = True;
               EraseIf( paths, to_delete );    }

          // Reconstruct joins from the paths.

          joins.clear_and_resize( paths.size( ) );
          vec< vec<int> > ends( paths.size( ) );
          for ( int i = 0; i < paths.isize( ); i++ )
          {    const vec<int>& p = paths[i];
               joins[i].push_back( p[0] );
               for ( int j = 0; j < p.isize( ) - 1; j++ )
               {    for ( int k = 0; k < G.From( p[j] ).isize( ); k++ )
                    {    if ( G.From( p[j] )[k] == p[j+1] )
                         {    joins[i].append( 
                                   G.EdgeObjectByIndexFrom( p[j], k ) );    }    }
                    joins[i].push_back( p[j+1] );    }
               ends[i] = p;    }
     
          // Fix overlaps.

          if ( verbosity >= 3 ) cout << "\n";
          vec< vec< pair<int,int> > > users(ne);
          for ( int i = 0; i < joins.isize( ); i++ )
          for ( int j = 0; j < joins[i].isize( ); j++ )
               users[ joins[i][j] ].push( i, j );
          for ( int e = 0; e < ne; e++ )
          for ( int i1 = 0; i1 < users[e].isize( ); i1++ )
          for ( int i2 = 0; i2 < users[e].isize( ); i2++ )
          {    if ( i2 == i1 ) continue;
               int j1 = users[e][i1].first, j2 = users[e][i2].first;
               vec<int> &x1 = joins[j1], &x2 = joins[j2];
               if ( x1.empty( ) || x2.empty( ) ) continue;
               int p1 = users[e][i1].second, p2 = users[e][i2].second;
               Bool mismatch = False;
               for ( int l1 = 0; l1 < x1.isize( ); l1++ )
               {    int l2 = l1 + p2 - p1;
                    if ( l2 < 0 || l2 >= x2.isize( ) ) continue;
                    if ( x1[l1] != x2[l2] ) mismatch = True;    }
               if ( !mismatch )
               {    if ( verbosity >= 3 )
                    {    cout << "see overlap between " << printSeq(x1)
                              << " and " << printSeq(x2) << endl;    }
                    if ( x1.size( ) >= x2.size( ) ) x2.resize(0);
                    else x1.resize(0);    }    }

          // Remove joins that behave strangely with respect to the inversion.

          SortSync( joins, ends );
          vec<Bool> dead( joins.size( ), False );
          for ( int i = 0; i < joins.isize( ); i++ )
          {    if ( joins[i].empty( ) || Inv( joins[i][0] ) < 0 ) continue;
               vec<int> x;
               for ( int j = joins[i].isize( ) - 1; j >= 0; j-- )
                    x.push_back( Inv( joins[i][j] ) );
               int p = BinPosition( joins, x );
               if ( p < 0 )
               {    if ( verbosity >= 3 )
                    {    cout << "problem with join" << endl;
                         cout << "join = " << printSeq(joins[i]) << endl;
                         cout << "inv = " << printSeq(x) << endl;    }
                    dead[i] = True;    }
               else
               {    vec<int> y = joins[i];
                    Sort(x), Sort(y);
                    if ( Meet( x, y ) ) dead[i] = dead[p] = True;    
                    else
                    {    vec<int> e;
                         for ( int j = 0; j < ends[i].isize( ); j++ )
                              e.push_back( Inv( ends[i][j] ) );
                         ends[p] = e;    }    }    }
     
          // Delete dead joins.

          for ( int i = 0; i < joins.isize( ); i++ )
               if ( joins[i].empty( ) ) dead[i] = True;
          EraseIf( joins, dead ), EraseIf( ends, dead );

          // List joins.

          if ( verbosity >= 2 )
          {    cout << "\nFINAL JOINS:\n";
               for ( int i = 0; i < joins.isize( ); i++ )
               {    const vec<int>& x = joins[i];
                    cout << "[" << i+1 << "] " << printSeq(x) << " (" 
                         << EdgeLengthKmers( x[0] );
                    for ( int j = 1; j < x.isize( ); j++ )
                         cout << "," << EdgeLengthKmers( x[j] );
                    cout << ") ends = " << printSeq( ends[i] ) << endl;    }    }
     
          // Edit graph.

          if ( verbosity >= 1 ) cout << "\n" << Date( ) << ": editing graph" << endl;
          vec<int> to_left, to_right;
          ToLeft(to_left), ToRight(to_right);
          for ( int i = 0; i < joins.isize( ); i++ )
          {    const vec<int>& x = joins[i];
               AddEdge( to_left[ x.front( ) ], to_right[ x.back( ) ], Cat(x) );
               if ( Inv( x[0] ) < 0 ) InvMutable( ).push_back(-1);
               else
               {    vec<int> x;
                    for ( int j = joins[i].isize( ) - 1; j >= 0; j-- )
                         x.push_back( Inv( joins[i][j] ) );
                    int p = BinPosition( joins, x );
                    ForceAssert( p >= 0 );
                    InvMutable( ).push_back( ne + p );    }
               DeleteEdges( ends[i], to_left );    }
     
          // Transform paths.

          if ( verbosity >= 1 ) cout << Date( ) << ": transforming paths" << endl;
          vec<int> news;
          for ( int i = 0; i < joins.isize( ); i++ )
               news.push_back( ne + i );
          TransformPaths( joins, news );
          UniqueOrderPaths( );

          // Clean up.

          if ( verbosity >= 1 ) cout << Date( ) << ": cleaning up" << endl;
          RemoveUnneededVertices( );
          RemoveEdgelessVertices( );
          RemoveDeadEdgeObjects( );
          if ( verbosity >= 1 ) cout << Date( ) << ": testing" << endl;
          TestValid(logc);
          DeleteReverseComplementComponents(logc);
          const double junk_ratio = 10.0;
          const int max_del = 1000;
          long_heuristics heur( "" );
          TrimHangingEnds( max_del, junk_ratio, heur, logc );
          DeleteWeakEdges(logc);
          RemoveSmallMostlyAcyclicComponents(logc);
          if ( verbosity >= 1 ) cout << Date( ) << ": done" << endl;    

          // Print some stats.

          if ( verbosity >= 1 )
          {    cout << EdgeObjectCount( ) << " edges and " << NComponents( )
                    << " components" << endl;
               {    vec<int> len;
                    for ( int e = 0; e < EdgeObjectCount( ); e++ )
                         len.push_back( EdgeLengthKmers(e) );
                    Sort(len);
                    cout << "N50 edge = " << N50(len) << endl;    }
               vec< vec<int> > comp;
               ComponentsE(comp);
               vec<int> compsizes;
               for ( int i = 0; i < comp.isize( ); i++ )
               {    const vec<int>& c = comp[i];
                    int s = 0;
                    for ( int j = 0; j < c.isize( ); j++ )
                         s += EdgeLengthKmers( c[j] );
                    compsizes.push_back(s);    }
               ReverseSort(compsizes);
               cout << "component sizes: " 
                    << printSeq(compsizes) << endl;    }    }    }
