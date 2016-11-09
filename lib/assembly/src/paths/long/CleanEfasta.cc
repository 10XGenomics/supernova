///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "Set.h"
#include "VecUtilities.h"
#include "efasta/EfastaTools.h"
#include "pairwise_aligners/SmithWatBandedA.h"
#include "paths/HyperEfasta.h"
#include "paths/long/CleanEfasta.h"
#include "paths/long/Heuristics.h"
#include "paths/long/Logging.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/SupportedHyperBasevector.h"
#include "paths/long/Variants.h"

void MakeEfastaAssembly( HyperEfasta& he, vec<int>& inv, vec<Bool>& hide,
     const vec<VariantSignature>& snp_bubbles, const long_heuristics& heur,
     const long_logging& logc )
{    double clock = WallClockTime( );
     Bool have_inv = False;
     for ( int i = 0; i < inv.isize( ); i++ )
          if ( inv[i] >= 0 ) have_inv = True;
     if (have_inv) CleanEfasta( he, inv, logc );
     else
     {    Reduce( he, logc.verb[ "REDUCE" ], logc );
          inv.resize( he.EdgeObjectCount( ), -1 );    }
     REPORT_TIME( clock, "used cleansing efasta" );
     FlagEdgesForHiding( he, inv, hide, logc );
     if (heur.DETECT_VARIANTS) MarkVariants( he, snp_bubbles, logc );    }

void CleanEfasta( HyperEfasta& he, vec<int>& inv, const long_logging& logc )
{    int K = he.K( );
     vec<int> to_left, to_right;
     he.ToLeft(to_left), he.ToRight(to_right);
     for ( int v = 0; v < he.N( ); v++ )
     {    if ( he.From(v).size( ) != 2 ) continue;
          int w = he.From(v)[0];
          if ( he.From(v)[1] != w ) continue;
          if ( !he.To(v).solo( ) || !he.From(w).solo( ) ) continue;
          int x = he.To(v)[0], y = he.From(w)[0];
          const efasta& e1 = he.EdgeObjectByIndexFrom( v, 0 );
          const efasta& e2 = he.EdgeObjectByIndexFrom( v, 1 );
          if ( e1.Contains( "{" ) || e2.Contains( "{" ) ) continue;
          basevector b1(e1), b2(e2);
          const int bw = 20;
          const int max_errs = 10;
          align a;
          int e;
          SmithWatBandedA( b1, b2, 0, bw, a, e, 0, 1, 1 );
          e += a.pos1( ) + a.pos2( )
               + b1.isize( ) - a.Pos1( ) + b2.isize( ) - a.Pos2( );
          // PRINT(e);
          if ( e > max_errs ) continue;
          vec<basevector> b;
          b.push_back( b1, b2 );
          // PrintVisualAlignment( True, cout, b1, b2, a );
          efasta left = he.EdgeObjectByIndexTo( v, 0 );
          efasta right = he.EdgeObjectByIndexFrom( w, 0 );
          Bool bad = False;
          for ( int j = 0; j < K - 1; j++ )
          {    if ( right[j] == '{' ) bad = True;
               if ( left[ left.isize( ) - 1 - j ] == '}' ) bad = True;    }
          if (bad) continue;
          efasta be(b);
          left.resize( left.isize( ) - (K-1) );
          right = right.substr( K-1, right.isize( ) - (K-1) );
          efasta enew = left + be + right;

          int p = he.EdgeObjectIndexByIndexTo( v, 0 );
          int x1 = he.EdgeObjectIndexByIndexFrom( v, 0 );
          int x2 = he.EdgeObjectIndexByIndexFrom( v, 1 );
          int q = he.EdgeObjectIndexByIndexFrom( w, 0 );
          int rp = inv[p], rx1 = inv[x1], rx2 = inv[x2], rq = inv[q];
          vec<int> s1, s2;
          s1.push_back( p, x1, x2, q );
          s2.push_back( rp, rx1, rx2, rq );
          UniqueSort(s1), UniqueSort(s2);
          vec<int> dels(s1);

          if ( rp < 0 )
          {    he.DeleteEdges(dels);
               inv.push_back(-1);
               he.AddEdge( x, y, enew );
               to_left.push_back(x), to_right.push_back(y);    }

          else if ( rp == q && s1 == s2 )
          {    he.DeleteEdges(dels);
               inv.push_back( he.EdgeObjectCount( ) );
               he.AddEdge( x, y, enew );
               to_left.push_back(x), to_right.push_back(y);    }

          else if ( Meet( s1, s2 ) )
          {    continue;    }

          else
          {    dels.append(s2);
               int rx = to_left[rq], ry = to_right[rp];
               he.DeleteEdges(dels);
               inv.push_back( he.EdgeObjectCount( ) + 1 );
               inv.push_back( he.EdgeObjectCount( ) );
               he.AddEdge( x, y, enew );
               to_left.push_back(x), to_right.push_back(y);
               efasta renew(enew);
               renew.ReverseComplement( );
               he.AddEdge( rx, ry, renew );
               to_left.push_back(rx), to_right.push_back(ry);    }    }

     // Remove unneeded vertices.

     for ( int i = 0; i < he.N( ); i++ )
     {    if ( he.From(i).size( ) == 1 && he.To(i).size( ) == 1 
               && he.From(i)[0] != i )
          {    efasta p = he.EdgeObjectByIndexTo( i, 0 );
               p.resize( p.size( ) - K + 1 );
               p.append( he.EdgeObjectByIndexFrom( i, 0 ) );
               int e1 = he.EdgeObjectIndexByIndexTo( i, 0 );
               int e2 = he.EdgeObjectIndexByIndexFrom( i, 0 );
               int enew = he.EdgeObjectCount( );
               int v = he.To(i)[0], w = he.From(i)[0];
               he.JoinEdges( i, p );
               to_left.push_back(v), to_right.push_back(w);
               if ( inv[e1] < 0 ) inv.push_back(-1);
               else
               {    ForceAssertGe( inv[e2], 0 );
                    int re1 = inv[e1], re2 = inv[e2];
                    if ( re1 != e1 && re1 != e2 )
                    {    ForceAssert( re2 != e1 && re2 != e2 );
                         int renew = he.EdgeObjectCount( );
                         efasta rp = he.EdgeObject(re2);
                         rp.resize( rp.isize( ) - K + 1 );
                         rp.append( he.EdgeObject(re1) );
                         int ri = to_right[re2];
                         he.JoinEdges( ri, rp );
                         int rv = to_left[re2], rw = to_right[re1];
                         to_left.push_back(rv), to_right.push_back(rw);
                         inv.push_back(renew, enew);    }    }    }    }

     // Remove dead edge objects.

     vec<Bool> used;
     he.Used(used);
     vec<int> to_new_id( used.size( ), -1 );
     {    int count = 0;
          for ( int i = 0; i < used.isize( ); i++ )
               if ( used[i] ) to_new_id[i] = count++;    }
     vec<int> inv2;
     for ( int i = 0; i < he.EdgeObjectCount( ); i++ )
     {    if ( !used[i] ) continue;
          if ( inv[i] == -1 ) inv2.push_back(-1);
          else inv2.push_back( to_new_id[ inv[i] ] );    }
     inv = inv2;
     he.RemoveDeadEdgeObjects( );
     he.RemoveEdgelessVertices( );    }

template<class H> void FlagEdgesForHiding( 
     const H& he, const vec<int>& inv, vec<Bool>& hide, const long_logging& logc )
{    double clock = WallClockTime( );
     vec<int> to_left, to_right;
     he.ToLeft(to_left), he.ToRight(to_right);

     // We break edges into five categories:
     // - keep_auto: not mapped by involution, or palindrome
     // - discard: to be hidden
     // - to_decide: edges not yet processed
     // - keep_look: edges to keep, but whose neighbors have not yet been looked at
     // - keep_looked: edges to keep whose neighbors have been examined.
     
     set<int> keep_auto, to_decide, keep_looked, keep_look, discard;

     // Define keep_auto, and initial to_decide.

     for ( int e = 0; e < he.EdgeObjectCount( ); e++ )
     {    if ( inv[e] < 0 || inv[e] == e ) keep_auto.insert(e);
          else to_decide.insert(e);    }

     // Main loop.

     while( keep_look.size( ) > 0 || to_decide.size( ) > 0 )
     {    if ( keep_look.size( ) > 0 )
          {    int e = *keep_look.begin( );
               keep_look.erase(e);
               keep_looked.insert(e);
               vec<int> nhood;
               int v = to_left[e], w = to_right[e];
               nhood.append( he.FromEdgeObj(v) );
               nhood.append( he.ToEdgeObj(v) );
               nhood.append( he.FromEdgeObj(w) );
               nhood.append( he.ToEdgeObj(w) );
               UniqueSort(nhood);
               for ( int j = 0; j < nhood.isize( ); j++ )
               {    int n = nhood[j];
                    if ( !Member( to_decide, n ) ) continue;
                    keep_look.insert(n), discard.insert( inv[n] );
                    to_decide.erase(n), to_decide.erase( inv[n] );    }    }
          else
          {    int e = *to_decide.begin( );
               to_decide.erase(e), to_decide.erase( inv[e] );
               keep_look.insert(e), discard.insert( inv[e] );    }    }

     // Report results.

     hide.resize( he.EdgeObjectCount( ), False );
     for ( set<int>::iterator i = discard.begin( ); i != discard.end( ); i++ )
          hide[*i] = True;
     REPORT_TIME( clock, "used flagging edges for hiding" );    }

template void FlagEdgesForHiding(
     const HyperEfasta& he, const vec<int>& inv, vec<Bool>& hide, 
     const long_logging& logc );

template void FlagEdgesForHiding(
     const SupportedHyperBasevector& he, const vec<int>& inv, vec<Bool>& hide,
     const long_logging& logc );

void CollapseBubbles( HyperEfasta& he )
{    const int max_brackets = 4;
     while(1)
     {    Bool changed = False;
          for ( int v = 0; v < he.N( ); v++ )
          {    for ( int l = 0; l < he.From(v).isize( ); l++ )
               {    int w = he.From(v)[l];
                    if ( w == v ) continue;
                    vec<int> ms, es;
                    for ( int m = 0; m < he.From(v).isize( ); m++ )
                    {    if ( he.From(v)[m] == w ) 
                         {    ms.push_back(m);
                              es.push_back( 
                                   he.EdgeObjectIndexByIndexFrom( v, m ) );    }    }
                    if ( es.solo( ) ) continue;

                    // Now we have two or more edges from vertex v to vertex w.

                    VecEFasta x;
                    for ( int m = 0; m < es.isize( ); m++ )
                         x.push_back( he.EdgeObject( es[m] ) );
                    int brackets = 0;
                    for ( size_t j = 0; j != x.size( ); j++ )
                    for ( int i = 0; i < x[j].isize( ); i++ )
                         if ( x[j][i] == '{' ) brackets++;
                    if ( brackets > max_brackets ) continue;
                    efasta e(x);
                    for ( int r = ms.isize( ) - 1; r >= 1; r-- )
                         he.DeleteEdgeFrom( v, ms[r] );
                    he.EdgeObjectByIndexFromMutable( v, ms[0] ) = e;
                    changed = True;
                    break;    }    }
          if ( !changed ) break;
          he.RemoveUnneededVertices( );    }
     he.RemoveDeadEdgeObjects( );
     he.RemoveEdgelessVertices( );    }

void GetCells( const HyperEfasta& he, vec<vec<int>>& cells )
{    cells.clear( );

     // Go through the vertices.

     for ( int v = 0; v < he.N( ); v++ )
     {    if ( he.From(v).solo( ) ) continue;

          // Now we have a vertex v, with at least two edges emanating from it.

          vec<int> suc, pre, between, betweenp;
          he.GetSuccessors1( v, suc );
          suc.EraseValue(v);
          for ( int wi = 0; wi < suc.isize( ); wi++ )
          {    int w = suc[wi];
               if ( he.To(w).solo( ) ) continue;

               // Now we have a vertex w which is not v, but follows v.  Moreover w
               // has at least two edges leading into it.  We posit v and w as
               // (possibly) bounding a cell.

               he.GetPredecessors1( w, pre );
               pre.EraseValue(w);
               between = Intersection( suc, pre );

               // Now between is the set of all vertices between v and w,
               // exclusive of them.

               betweenp = between;
               betweenp.push_back( v, w );
               Sort(betweenp);
               Bool bad = False;
               for ( int j = 0; j < between.isize( ); j++ )
               {    int x = between[j];
                    if ( !BinSubset( he.To(x), betweenp ) ) bad = True;
                    if ( !BinSubset( he.From(x), betweenp ) ) bad = True;    }
               if ( !BinSubset( he.From(v), betweenp ) ) bad = True;
               if ( !BinSubset( he.To(w), betweenp ) ) bad = True;
               if (bad) continue;

               // Now we know that v and w bound a cell.

               vec<int> cell;
               cell.push_back(v);
               cell.append(between);
               cell.push_back(w);
               cells.push_back(cell);    }    }    }

void Reduce( HyperEfasta& he, const int verbosity, const long_logging& logc )
{
     // First collapse bubbles.

     if ( verbosity >= 1 )
     {    cout << Date( ) << ": have " << he.EdgeObjectCount( ) << " edges, "
               << "collapsing bubbles" << endl;    }
     CollapseBubbles(he);

     // Now collapse more complex cells.

     if ( verbosity >= 1 ) 
     {    cout << Date( ) << ": have " << he.EdgeObjectCount( ) << " edges, "
               << "looking for cells to collapse" << endl;    }
     const int max_paths = 50;
     while(1)
     {    Bool progress = False;
          vec< vec<int> > cells;

          // Print.

          if ( verbosity >= 3 )
          {    cout << "\nat top of Reduce loop, graph is:\n";
               for ( int v = 0; v < he.N( ); v++ )
               for ( int j = 0; j < he.From(v).isize( ); j++ )
               {    int w = he.From(v)[j];
                    int e = he.EdgeObjectIndexByIndexFrom( v, j );
                    cout << v << " --(" << e << ",kmers=" 
                         << he.EdgeLengthKmers(e) 
                         << ", amb_events = " << he.EdgeObject(e).AmbEventCount( )
                         << ")--> " << w << endl;    }
               cout << endl;    }

          // Get cells.

          GetCells( he, cells );
          if ( verbosity >= 1 ) 
               cout << Date( ) << ": found " << cells.size( ) << " cells" << endl;

          // Go through the cells.

          for ( int i = 0; i < cells.isize( ); i++ )
          {    const vec<int>& cell = cells[i];
               if ( verbosity >= 2 )
                    cout << "examining cell " << printSeq(cell) << endl;
               if ( verbosity >= 3 )
               {    cout << "\ngraph is:\n";
                    for ( int v = 0; v < he.N( ); v++ )
                    for ( int j = 0; j < he.From(v).isize( ); j++ )
                    {    int w = he.From(v)[j];
                         int e = he.EdgeObjectIndexByIndexFrom( v, j );
                         cout << v << " --(" << e << ",kmers=" 
                              << he.EdgeLengthKmers(e) 
                              << ", amb_events = " 
                              << he.EdgeObject(e).AmbEventCount( )
                              << ")--> " << w << endl;    }
                    cout << endl;    }

               // Test for loop at w.  This is in effect working around a bug in
               // EdgePaths, which incorrectly ignores loops at w.

               if ( he.LoopAt( cell.back( ) ) ) continue;

               vec< vec<int> > paths;
               const int max_iterations = 1000;
               // note could upgrade to faster version that uses to_left, to_right
               Bool OK = he.EdgePaths( cell.front( ), cell.back( ), paths, 
                         -1, max_paths, max_iterations );
               if ( verbosity >= 1 ) PRINT2( int(OK), paths.size( ) );
               if ( !OK ) continue;
               if ( paths.isize( ) > max_paths ) continue;
               if ( paths.size( ) <= 1 ) continue;
               Bool bad = False;
               vec<int> edges;
               for ( int j = 0; j < paths.isize( ); j++ )
               for ( int k = 0; k < paths[j].isize( ); k++ )
                    edges.push_back( paths[j][k] );
               UniqueSort(edges);
               if ( verbosity >= 2 )
                    cout << "see edges " << printSeq(edges) << endl;
               VecEFasta bpaths;
               Bool fail = False;
               int K = he.K( );
               for ( int j = 0; j < paths.isize( ); j++ )
               {    if (fail) break;
                    efasta b = he.EdgeObject( paths[j][0] );
                    for ( int k = 1; k < paths[j].isize( ); k++ )
                    {    efasta c = he.EdgeObject( paths[j][k] );
                         if ( b.isize( ) < K - 1 )
                         {    fail = True;
                              break;    }
                         for ( int i = 0; i < K - 1; i++ )
                         {    if ( b[ b.isize( ) - 1 - i ] == '}' )
                              {    fail = True;
                                   break;    }    }
                         if ( !fail ) b.resize( b.isize( ) - ( K - 1 ) );
                         else
                         {    fail = False;
                              if ( c.isize( ) >= K )
                              {    for ( int i = 0; i < K - 1; i++ )
                                        if ( c[i] == '{' ) fail = True;    }
                              if ( !fail ) 
                                   c = c.substr( K-1, c.isize( ) - (K-1) );    }
                         if (fail) break;
                         b += c;    }
                    bpaths.push_back(b);    }
               if (fail) continue;
               int bracks = 0;
               for ( size_t i = 0; i != bpaths.size( ); i++ )
               for ( int j = 0; j < bpaths[i].isize( ); j++ )
                    if ( bpaths[i][j] == '{' ) bracks++;
               const int max_bracks = 10;
               if ( bracks > max_bracks ) continue;
               if ( verbosity >= 2 ) cout << "combining" << endl;
               efasta e(bpaths);
               he.DeleteEdges(edges);
               he.AddEdge( cell.front( ), cell.back( ), e );
               progress = True;
               he.RemoveUnneededVertices( );
               he.RemoveDeadEdgeObjects( );
               he.RemoveEdgelessVertices( );
               break;    }
          if ( !progress ) break;    }

     // Collapse bubbles again.

     CollapseBubbles(he);    

     // Find small unresolved but isolated cyclic edge clumps and 'box' them up in
     // a fastg-like contruction (no longer efasta).  Currently very restrictive.
     //
     // First define heuristics.

     const int min_long = 1000;
     const int max_short = 450;
     const int max_dist = 4;
     const int max_int = 15;

     // Repeat until no further progress.

     while(1)
     {    Bool improved = False;
          Bool deflower_verbose = ( logc.verb[ "DEFLOWER" ] >= 1 );
          if (deflower_verbose)
          {    cout << "\n";
               for ( int v = 0; v < he.N( ); v++ )
               for ( int j = 0; j < he.From(v).isize( ); j++ )
               {    int w = he.From(v)[j];
                    int e = he.EdgeObjectIndexByIndexFrom( v, j );
                    int len = he.EdgeLengthKmers(e);
                    PRINT4( v, w, e, len );    }
               cout << "\n";    }

          // Look for a long edge that will function as a left boundary.

          int K = he.K( );
          for ( int u = 0; u < he.N( ); u++ )
          for ( int j = 0; j < he.From(u).isize( ); j++ )
          {    int v = he.From(u)[j];
               int e1 = he.EdgeObjectIndexByIndexFrom( u, j );
               int n1 = he.EdgeObject(e1).isize( ) - K + 1;
               if ( n1 < min_long ) continue;
               if (deflower_verbose) PRINT(e1);
               
               // The left edge cannot have a bracket near its right end.
     
               Bool bad = False;
               for ( int s = 0; s < K - 1; s++ )
               {    if ( he.EdgeObject(e1)[ 
                         he.EdgeObject(e1).isize( ) - s - 1 ] == '}' )
                    {    bad = True;     }    }
               if (bad) continue;

               // Look for all vertices within distance d of v.

               vec< pair<int,int> > vd;
               vd.push( v, 0 );
               while(1)
               {    Bool progress = False;
                    for ( int i = 0; i < vd.isize( ); i++ )
                    {    int x = vd[i].first, d = vd[i].second;
                         if ( d == max_dist ) continue;
                         for ( int l = 0; l < he.From(x).isize( ); l++ )
                         {    int y = he.From(x)[l];
                              Bool found = False;
                              for ( int m = 0; m < vd.isize( ); m++ )
                              {    if ( vd[m].first == y )
                                   {    if ( d + 1 < vd[m].second )
                                        {    vd[m].second = d + 1;  
                                             progress = True;    }
                                        found = True;    }    }
                              if ( !found )
                              {    vd.push( y, d + 1 );
                                   progress = True;    }    }    }
                    if ( !progress ) break;    }
     
               // Now look for a long edge that could function as a right boundary.
     
               Bool edited = False;
               for ( int i = 0; i < vd.isize( ); i++ )
               {    if (edited) break;
                    int w = vd[i].first;
                    for ( int l = 0; l < he.From(w).isize( ); l++ )
                    {    int e2 = he.EdgeObjectIndexByIndexFrom( w, l );
                         if ( e2 == e1 ) continue;
                         int n2 = he.EdgeObject(e2).isize( ) - K + 1;
                         if ( n2 < min_long ) continue;
                         if (deflower_verbose) PRINT(e2);
     
                         // Now we have long edges that might function as left/right 
                         // boundaries.  Define what's between them.
          
                         vec<int> vint, eint;
                         vint.push_back( v, w );
                         //Bool get_back = False;
                         for ( int m = 0; m < vint.isize( ); m++ )
                         {    int x = vint[m];
                              for ( int r = 0; r < he.From(x).isize( ); r++ )
                              {    int y = he.From(x)[r];
                                   //if ( y == v ) get_back = True;
                                   int e = he.EdgeObjectIndexByIndexFrom( x, r );
     
                                   // Don't follow the right edge.
     
                                   if ( e == e2 ) continue;
     
                                   // Save vertex and edge.
     
                                   if (deflower_verbose)
                                        cout << "using edge " << e << endl;
                                   eint.push_back(e);
                                   if ( !Member( vint, y ) ) vint.push_back(y);
                                   if ( vint.isize( ) > max_int ) break;    }    }
                         if (deflower_verbose) PRINT( vint.size( ) );
     
                         // We have to be able to get back to v, and there can't be
                         // too many vertices.
     
                         if ( vint.isize( ) > max_int /* || !get_back */ ) continue;
     
                         // Going backwards can't take us out of vint.  Keep adding 
                         // to edge list.
     
                         Bool bad = False;
                         for ( int m = 0; m < vint.isize( ); m++ )
                         {    int x = vint[m];
                              for ( int r = 0; r < he.To(x).isize( ); r++ )
                              {    int e = he.EdgeObjectIndexByIndexTo( x, r );
                                   if ( e == e1 ) continue;
                                   eint.push_back(e);
                                   int y = he.To(x)[r];
                                   if ( !Member( vint, y ) ) bad = True;    }    }
                         if (deflower_verbose) PRINT(int(bad));
                         if (bad) continue;
     
                         // Make sure that edges are short.
     
                         UniqueSort(eint);
                         if (deflower_verbose) 
                              cout << "eint = {" << printSeq(eint) << "}" << endl;
                         Bool OK = True;
                         for ( int r = 0; r < eint.isize( ); r++ )
                         {    if ( he.EdgeLengthKmers( eint[r] ) > max_short ) 
                                   OK = False;    }
                         if (deflower_verbose) PRINT(int(OK));
                         if ( !OK ) continue;

                         // Create pseudo-fastg for the intermediate stuff.

                         if (deflower_verbose) cout << "deflowering" << endl;
                         if ( vint.size( ) > 2 )
                              swap( vint[1], vint[ vint.isize( ) - 1 ] );
                         String mid = "{";
                         for ( int m = 0; m < vint.isize( ); m++ )
                         {    int x = vint[m];
                              for ( int r = 0; r < he.From(x).isize( ); r++ )
                              {    int y = he.From(x)[r];
                                   int e = he.EdgeObjectIndexByIndexFrom( x, r );
                                   if ( e == e2 ) continue;
                                   int n = Position( vint, y );
                                   vec<basevector> b;
                                   he.EdgeObject(e).ExpandTo(b);
                                   for ( int z = 0; z < b.isize( ); z++ )
                                   {    if ( mid.size( ) > 1 ) mid += ",";
                                        mid += "(" + ToString(m) + "," 
                                        + ToString(n) + ",";
                                        for ( int s = 0; 
                                             s < b[z].isize( ) - ( K - 1 ); s++ )
                                        {    mid += as_base( b[z][s] );    }    
                                        mid += ")";
                                   }
                              }
                         }
                         mid += "}";

                         // Edit the assembly.
     
                         String left = he.EdgeObject(e1), right = he.EdgeObject(e2);
                         left.resize( left.isize( ) - ( K - 1 ) );
                         he.AddEdge( u, he.From(w)[l], left + mid + right );
                         vec<int> dels = eint;
                         dels.push_back( e1, e2 );
                         he.DeleteEdges(dels);
                         edited = True;
                         improved = True;
                         break;    }    }    }
          if ( !improved ) break;    }

     // Clean up.

     he.RemoveDeadEdgeObjects( );
     he.RemoveEdgelessVertices( );    }
