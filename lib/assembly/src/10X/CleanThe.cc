// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"
#include "Equiv.h"
#include "ParallelVecUtilities.h"
#include "Set.h"
#include "TokenizeString.h"
#include "feudal/PQVec.h"
#include "graph/DigraphTemplate.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "paths/long/large/GapToyTools.h"
#include "paths/long/large/Lines.h"
#include "10X/Capture.h"
#include "10X/CleanThe.h"
#include "10X/Gap.h"
#include "10X/Heuristics.h"
#include "10X/LineLine.h"
#include "10X/LineOO.h"
#include "10X/LocalTools.h"
#include "10X/MakeLocalsTools.h"
#include "10X/Peeler.h"
#include "10X/PlaceReads.h"
#include "10X/PullApart.h"
#include "10X/Scaffold.h"
#include "10X/SecretOps.h"
#include "10X/Super.h"
#include "10X/astats/FinAlign.h"
#include "10X/astats/RefLookup.h"
#include "10X/DfTools.h"
#include "10X/mergers/ShortMergers.h"
#include "10X/LineGraphOps.h"
#include "10X/ThreadedLogger.h"

void RepairShop( const HyperBasevectorX& hb, const vec<int>& inv, 
     digraphE<vec<int>>& D, vec<int>& dinv )
{    
     cout << Date( ) << ": welcome to the repair shop" << endl;

     // Remove certain unneeded vertices at cell.cell.  This creates cells having 
     // unneeded vertices and edges.  Would require more patience to clean up.
     // Also punting on cells within cells (shouldn't be allowed) and sequence gap
     // edges (which should be allowed).

     while(1)
     {    vec<int> to_left, to_right, dels;
          D.ToLeft(to_left), D.ToRight(to_right);
          int repairs = 0;
          for ( int v = 0; v < D.N( ); v++ )
          {    if ( D.To(v).solo( ) && D.From(v).solo( ) )
               {    int x = D.To(v)[0], y = D.From(v)[0];
                    if ( !IsUnique( vec<int>{ x, v, y } ) ) continue;
                    int d1 = D.ITo(v,0), d2 = D.IFrom(v,0);
                    if ( !IsCell( D.O(d1) ) || !IsCell( D.O(d2) ) ) continue;
                    int rd1 = dinv[d1], rd2 = dinv[d2];
                    if ( !IsUnique( vec<int>{ d1, d2, rd1, rd2 } ) ) continue;
                    cell c1, c2;
                    c1.CellDecode( D.O(d1) ), c2.CellDecode( D.O(d2) );
                    digraphE<vec<int>> G = c1.G( );
                    G.Append( c2.G( ) );
                    Bool ok = True;
                    {    for ( int g = 0; g < G.E( ); g++ )
                         {    if ( IsCell( G.O(g) ) || IsSequence( G.O(g) ) )
                                   ok = False;    }    }
                    if ( !ok ) continue;
                    repairs++;
                    int n1 = c1.G( ).N( );
                    int ll = c1.Left( ), lr = c1.Right( );
                    int rl = c2.Left( ) + n1, rr = c2.Right( ) + n1;
                    G.TransferEdges( lr, rl );
                    cell c( G, ll, rr );
                    digraphE<vec<int>> RG(G);
                    RG.Reverse( );
                    for ( auto& x : RG.EdgesMutable( ) )
                         if ( x[0] >= 0 ) for ( auto& e : x ) e = inv[e];
                    cell rc( RG, rr, ll );
                    int ry = to_right[rd1], rx = to_left[rd2];
                    vec<int> z, rz;
                    c.CellEncode(z), rc.CellEncode(rz);
                    dels.push_back( d1, d2, rd1, rd2 );
                    int E = D.E( );
                    D.AddEdgeWithUpdate( x, y, z, to_left, to_right );
                    D.AddEdgeWithUpdate( rx, ry, rz, to_left, to_right );
                    dinv.push_back(E+1);
                    dinv.push_back(E);    }    }
          if ( repairs > 0 )
          {    D.DeleteEdges(dels);
               RemoveUnneededVertices( D, dinv );
               CleanupCore( D, dinv );
               cout << Date( ) << ": made " << repairs << " cell.cell repairs" 
                    << endl;    }
          else break;    }    }

// Caution: the following could put either the line ids or line lengths in as
// the edge objects.  Make sure its does what you need.

void BuildLineGraph( const digraphE<vec<int>>& D, const vec<int>& dinv,
     const vec<vec<vec<vec<int>>>>& dlines, digraphE<int>& D2 )
{    vec<int> to_left, to_right;
     D.ToLeft(to_left), D.ToRight(to_right);
     vec<int> verts( 2 * dlines.size( ) );
     #pragma omp parallel for
     for ( int i = 0; i < dlines.isize( ); i++ )
     {    int v = to_left[ dlines[i].front( )[0][0] ];
          int w = to_right[ dlines[i].back( )[0][0] ];
          verts[2*i] = v, verts[2*i+1] = w;    }
     UniqueSort(verts);
     int N = verts.size( );
     vec<vec<int>> from(N), to(N), from_edge_obj(N), to_edge_obj(N);
     for ( int i = 0; i < dlines.isize( ); i++ )
     {    int v = BinPosition( verts, to_left[ dlines[i].front( )[0][0] ] );
          int w = BinPosition( verts, to_right[ dlines[i].back( )[0][0] ] );
          from[v].push_back(w), to[w].push_back(v);
          from_edge_obj[v].push_back(i), to_edge_obj[w].push_back(i);    }
     #pragma omp parallel for
     for ( int v = 0; v < N; v++ )
     {    SortSync( from[v], from_edge_obj[v] );
          SortSync( to[v], to_edge_obj[v] );    }
     vec<int> edges( dlines.size( ), vec<int>::IDENTITY );
     D2.Initialize( from, to, edges, to_edge_obj, from_edge_obj );    }

void SnipFlipSquares( digraphE<vec<int>>& D, vec<int>& dinv,
     const ReadPathVec& dpaths, const Bool verbose )
{
     // Heuristics.

     const int MAX_CAN_INS_DEL = 5;
     const int MIN_CAN_INS_RATIO = 4;

     // Set up.

     cout << Date( ) << ": start to snip flip squares" << endl << endl;
     vec<vec<vec<vec<int>>>> dlines;
     FindLinesFixed( D, dinv, dlines, 
          MAX_CELL_PATHS_EVALUATION, MAX_CELL_DEPTH_EVALUATION );
     digraphE<int> D2;
     BuildLineGraph( D, dinv, dlines, D2 );
     vec<int> linv;
     LineInv( dlines, dinv, linv );
     IntIndex dpaths_index( dpaths, D.E( ) );

     // Traverse vertices of line graph and identify edges to be deleted.

     vec<int> dels;
     cout << Date( ) << ": identify edges to be deleted" << endl;
     for ( int v = 0; v < D2.N( ); v++ )
     for ( int m = 0; m < D2.From(v).isize( ); m++ )
     {    
          // Suppose that two lines exit v, and one is a single nongap edge, 
          // leading to w.  Suppose also that one line enters v, that one other line 
          // enters w, and that these two lines are inverse to each other.

          if ( D2.From(v).size( ) != 2 ) continue;
          if ( !dlines[ D2.IFrom(v,m) ].solo( ) ) continue;
          int l1 = D2.IFrom(v,m), l2 = D2.IFrom(v,1-m);
          int w = D2.From(v)[m];
          if ( !D2.To(v).solo( ) || D2.To(w).size( ) != 2 ) continue;
          int il = 0;
          if ( D2.ITo(w,il) == l1 ) il = 1;

          int l3 = D2.ITo(v,0), l4 = D2.ITo(w,il);
          if ( l4 != linv[l3] ) continue;
          int d1 = dlines[l1][0][0][0], d2 = dlines[l2][0][0][0];
          if ( D.O(d1)[0] < 0 ) continue;
          int g = dlines[l3].back( )[0][0];

          // If d1 is very weakly supported, kill it.

          int n1 = 0, n2 = 0;
          for ( int mpass = 1; mpass <= 2; mpass++ )
          for ( int pass = 1; pass <= 2; pass++ )
          {    int f;
               if ( mpass == 1 ) f = ( pass == 1 ? d1 : dinv[d1] );
               else f = ( pass == 1 ? d2 : dinv[d2] );
               for ( int j = 0; j < dpaths_index.Count(f); j++ )
               {    int64_t id = dpaths_index.Val( f, j );
                    const ReadPath& p = dpaths[id];
                    for ( int l = 0; l < (int) p.size( ) - 1; l++ )
                    {    if ( pass == 1 )
                         {    if ( p[l] == g && p[l+1] == f )
                              {    ( mpass == 1 ? n1 : n2 )++;
                                   break;    }    }
                         else
                         {    if ( p[l] == f && p[l+1] == dinv[g] )
                              {   ( mpass == 1 ? n1 : n2 )++;
                                   break;    }    }    }    }    }
          if (verbose) PRINT3( d1, n1, n2 );
          if ( n1 > MAX_CAN_INS_DEL ) continue;
          if ( n2 == 0 || n2 < MIN_CAN_INS_RATIO * n1 ) continue;
          if (verbose) 
               cout << "deleting edges " << d1 << " and " << dinv[d1] << endl;
          dels.push_back( d1, dinv[d1] );    }

     // Delete the edges.

     UniqueSort(dels);
     cout << Date( ) << ": snipping " << dels.size( ) << " edges" << endl;
     D.DeleteEdges(dels);
     RemoveUnneededVertices( D, dinv );
     CleanupCore( D, dinv );    }

void BarcodeJoin( const HyperBasevectorX& hb, const vec<int>& inv,
     digraphE<vec<int>>& D, vec<int>& dinv,
     ReadPathVec& dpaths, const vec<Bool>& dup, const ReadPathVecX& pathsx, 
     const vec<Bool>& pmask, const vec<int64_t>& bci, 
     const vec< triple<int,int,int> >& qept, const VecIntVec& ebcx,
     const Bool break_trinads )
{
     // Heuristics.

     const int MIN_BIG = 10000;
     const int MIN_LEN = 4000;  // min length of line in nhood
     const double MAX_CN_DIFF = 0.25;
     const int MIN_LEFT_IGNORE = 100;
     const int MIN_ADVANTAGE = 100;
     const int MAX_DEPTH = 25;

     // Make ancillary stuff.

     vec<vec<vec<vec<int>>>> dlines;
     FindLinesFixed( 
          D, dinv, dlines, MAX_CELL_PATHS_EVALUATION, MAX_CELL_DEPTH_EVALUATION );
     vec<int> tol;
     MakeTol( D, dlines, tol );
     cout << Date( ) << ": expanding barcode index" << endl;
     vec<int32_t> bc( bci.back( ), -1 );
     #pragma omp parallel for
     for ( int b = 0; b < bci.isize( ) - 1; b++ )
     {    int64_t start = bci[b], stop = bci[b+1];
          for ( int64_t j = start; j < stop; j++ )
               bc[j] = b;    }

     // Compute line copy number.

     ThreadedLogger log;
     double clock = WallClockTime( );
     MasterVec<SerfVec<pair<int,int>>> lbpx;
     vec<double> cov;
     cout << Date( ) << ": BEGIN COPY NUMBER CALCULATION BLOCK" << endl;
     vec<int> kmers( hb.E( ) ), llens;
     #pragma omp parallel for
     for ( int e = 0; e < kmers.isize( ); e++ )
          kmers[e] = hb.Kmers(e);
     vec<vec<pair<int,int>>> lbp;
     PlaceReadsMasked( hb, D, dup, pathsx, pmask, dpaths );
     {    IntIndex dpaths_index( dpaths, D.E( ) );
          BarcodePos( bc, hb, D, dinv, dpaths, dpaths_index, dlines, lbp, 0 );    }
     {    for ( auto x : lbp )
          {    SerfVec<pair<int,int>> y( x.size( ) );
               for ( int j = 0; j < x.isize( ); j++ )
                    y[j] = x[j];
               lbpx.push_back(y);    }
          GetLineLengths( hb, D, dlines, llens );
          LineCN( kmers, lbpx, D, dlines, llens, cov );    }
     cout << Date( ) << ": END COPY NUMBER CALCULATION BLOCK" << endl;
     cout << Date( ) << ": total time used computing copy number = "
          << TimeSince(clock) << endl;

     // Compute proximate lines.

     vec< vec< pair<int,int> > > lhood;
     LineProx( hb, inv, ebcx, D, dinv, dlines, qept, lhood );

     // Build lineline graph.

     digraphE<int> D2;
     BuildLineGraph( D, dinv, dlines, D2 );
     vec<int> linv, to_left, to_right, to_left2, to_right2;
     LineInv( dlines, dinv, linv );
     D.ToLeft(to_left), D.ToRight(to_right);
     D2.ToLeft(to_left2), D2.ToRight(to_right2);

     // Find unambiguous joins.

     Bool verbose = False;
     cout << Date( ) << ": looking for links, mem = " << MemUsageGBString( ) << endl;
     if (verbose) cout << endl;
     vec< pair<int,int> > links;

     int nthreads = omp_get_max_threads( );
     vec<Bool> processed( dlines.size( ), False );
     #pragma omp parallel for schedule(dynamic,1)
     for ( int bi = 0; bi < nthreads; bi++ )
     {
          // Go through the lines.

          int L = 0;
          while(1)
          {    
               // Find next unprocessed line.

               if ( L == dlines.isize( ) ) break;
               if ( processed[L] )
               {    L++;
                    continue;    }
               else
               {
                    Bool next = False;
                    #pragma omp critical
                    {    if ( processed[L] )
                         {    L++;
                              next = True;    }
                         else processed[L] = True;    }
                    if (next) continue;    }

               // Process it.

               if ( llens[L] < MIN_BIG ) continue;
               vec<pair<int,int>> LH;
               for ( int i = 0; i < lhood[L].isize( ); i++ )
               {    int L2 = lhood[L][i].second;
                    if ( L2 == L || L2 == linv[L] ) continue;
                    if ( llens[L2] >= MIN_LEN ) LH.push_back( lhood[L][i] );    }
     
               // Define right reach of L.

               vec<int> reach, vs = { to_right[ dlines[L].back( )[0][0] ] };
               int dp = 0, vsp = 0;
               Bool fail = False;
               for ( ; dp < MAX_DEPTH; dp++ )
               {    int n = vs.size( );
                    if ( vsp == n ) break;
                    for ( ; vsp < n; vsp++ )
                    {    int v = vs[vsp];
                         for ( int j = 0; j < D.From(v).isize( ); j++ )
                         {    int L2 = tol[ D.IFrom(v,j) ];
                              if ( L2 < 0 ) // don't know if this can happen
                              {    fail = True;
                                   break;    }
                              if ( llens[L2] >= MIN_LEN ) reach.push_back(L2);
                              else 
                              {    vs.push_back( to_right[ 
                                        dlines[L2].back( )[0][0] ] );   
                                        }    }    }    }
               if ( dp == MAX_DEPTH || fail ) reach.clear( );
               UniqueSort(reach);

               // Now search.

               Bool confused = False;
               vec<int> X;
               vec<Bool> good;
               for ( int i = 0; i < LH.isize( ); i++ )
               {    int L2 = LH[i].second;
     
                    // Ignore if L2 too short.
     
                    if ( llens[L2] < MIN_LEN ) continue;
     
                    // Compute scores.
                    
                    vec<double> scores;
                    vec< triple<int,int,int> > M;
                    scores.push_back( ScoreOrder( { L2, L }, lbp, llens, M ) );
                    scores.push_back( ScoreOrder( { linv[L2], L }, lbp, llens, M ) );
                    scores.push_back( ScoreOrder( { L, L2 }, lbp, llens, M ) );
                    scores.push_back( ScoreOrder( { L, linv[L2] }, lbp, llens, M ) );
     
                    // If we're pretty sure L2 goes on the left, ignore, and if we
                    // can't tell, give up.
     
                    double left_adv = 
                         Min( scores[2], scores[3] ) - Min( scores[0], scores[1] );
                    if ( reach.nonempty( ) && !BinMember( reach, L2 )
                         && left_adv < MIN_LEFT_IGNORE 
                         && left_adv > -MIN_LEFT_IGNORE )
                    {    continue;    }
                    if ( left_adv >= MIN_LEFT_IGNORE ) continue;
                    if ( left_adv > -MIN_LEFT_IGNORE )
                    {    confused = True;
                         break;    }

                    // Find winner.

                    vec<int> ids( 4, vec<int>::IDENTITY );
                    SortSync( scores, ids );
                    double win = scores[1] - scores[0];
     
                    // Save.

                    X.push_back( ids[0]==2 ? L2 : linv[L2] );
                    good.push_back( win >= MIN_ADVANTAGE );    }
               if (confused) continue;
               if ( X.empty( ) ) continue;
     
               // Find leftmost L2.
     
               if ( X.size( ) > 1 )
               {    for ( int j2 = 0; j2 < X.isize( ); j2++ )
                    {    int L2 = X[j2];
                         Bool confused = False;
                         for ( int j3 = 0; j3 < X.isize( ); j3++ )
                         {    int L3 = X[j3];
                              if ( L3 == L2 ) continue;
                              vec<double> scores;
                              vec< triple<int,int,int> > M;
                              scores.push_back( ScoreOrder( 
                                   { L3, L2 }, lbp, llens, M ) );
                              scores.push_back( ScoreOrder( 
                                   { linv[L3], L2 }, lbp, llens, M ) );
                              scores.push_back( ScoreOrder( 
                                   { L2, L3 }, lbp, llens, M ) );
                              scores.push_back( 
                                   ScoreOrder( { L2, linv[L3] }, lbp, llens, M ) );
                              double left_adv = Min( scores[2], scores[3] ) 
                                   - Min( scores[0], scores[1] );
                              if ( reach.nonempty( ) && !BinMember( reach, L3 ) &&
                                   left_adv >= -MIN_LEFT_IGNORE && left_adv <= 0 )
                              {    continue;    }
                              if ( left_adv >= -MIN_LEFT_IGNORE )
                              {    confused = True;
                                   break;    }    }
                         if ( !confused )
                         {    if ( good[j2] ) X = {L2};
                              break;    }    }   }
               if ( X.size( ) > 1 ) continue;

               // Save if link is to something long and CN is close.

               int L2 = X[0];
               if ( llens[L2] >= MIN_BIG && Abs( cov[L] - cov[L2] ) < MAX_CN_DIFF )
               {
                    #pragma omp critical
                    {    links.push( L, L2 );    }    }    }    }

     Sort(links);

     // Break trinads.

     if (break_trinads)
     {    for ( auto x : links )
          for ( int i = 0; i < links.isize( ); i++ )
          {    int L1 = links[i].first, L2 = links[i].second;
               if ( to_right2[L1] != to_left2[L2] ) continue;
               int v = to_right2[L1];
               if ( D.To(v).size( ) != 2 || !D.From(v).solo( ) ) continue;
               int L3 = -1;
               for ( int j = 0; j < 2; j++ )
                    if ( D2.ITo(v,j) != L1 ) L3 = D2.ITo(v,j);
               int d = dlines[L3].back( )[0][0];
               int rd = dinv[d];
               if ( d == rd || to_right[d] == to_left[rd] ) continue;
               int N = D.N( );
               D.AddVertices(2);
               D.GiveEdgeNewToVx( d, to_right[d], N );
               D.GiveEdgeNewFromVx( rd, to_left[rd], N+1 );    }
          RemoveUnneededVertices( D, dinv );
          CleanupCore( D, dinv );
          Validate( hb, inv, D, dinv );
          return;    }

     // Remove asymmetric links.

     vec<Bool> to_delete( links.size( ), False );
     for ( auto x : links )
     for ( int i = 0; i < links.isize( ); i++ )
     {    int L1 = links[i].first, L2 = links[i].second;
          if ( !BinMember( links, make_pair( linv[L2], linv[L1] ) ) )
               to_delete[i] = True;    }
     EraseIf( links, to_delete );

     // Announce.

     if (verbose)
     {    for ( auto x : links )
          {    int L1 = x.first, L2 = x.second;
               {    cout << "recommend linking from L" << L1 << " to L" << L2
                         << endl;    }    }
          cout << endl;    }

     // Make the links (or more precisely, some of them).

     int nlinks = 0;
     for ( int i = 0; i < links.isize( ); i++ )
     {    int L1 = links[i].first, L2 = links[i].second;
          int RL1 = linv[L1], RL2 = linv[L2];
          if ( !IsUnique( vec<int>{ L1, L2, RL1, RL2 } ) ) continue;
          if ( make_pair( RL2, RL1 ) < make_pair( L1, L2 ) ) continue;
          int d1 = dlines[L1].back( )[0][0], d2 = dlines[L2].front( )[0][0];
          int v = to_right[d1], w = to_left[d2];

          // Type 1.  Connect simple dead ends.

          if ( D.From(v).empty( ) && D.To(w).empty( )
               && D.To(v).solo( ) && D.From(w).solo( ) )
          {    if (verbose)
                    cout << "\ntype 1 linking from L" << L1 << " to L" << L2 << endl;
               nlinks++;
               int E = D.E( );
               int rd1 = dlines[RL1].front( )[0][0], rd2 = dlines[RL2].back( )[0][0];
               int rv = to_left[rd1], rw = to_right[rd2];
               D.AddEdgeWithUpdate( v, w, vec<int>{-2}, to_left, to_right );
               D.AddEdgeWithUpdate( rw, rv, vec<int>{-2}, to_left, to_right );
               dinv.push_back( E+1, E );
               continue;    }

          // Type 2.  Connect if there's something in between.

          vec<pair<int,int>> nhood; // (line, dist from right end of L1)
          map<int,int> distx;
          nhood = { make_pair(L1,0) };
          distx[L1] = 0;
          int x = to_right2[L1], y = to_left2[L2];
          for ( int j = 0; j < nhood.isize( ); j++ )
          {    int L = nhood[j].first, dist = nhood[j].second;
               int z = to_right2[L];
               for ( int m = 0; m < D2.From(z).isize( ); m++ )
               {    int LP = D2.IFrom(z,m);
                    int distp = dist + llens[LP];
                    if ( distp >= MIN_BIG ) continue;
                    if ( distx.find(LP) == distx.end( ) )
                         nhood.push(LP,distp), distx[LP] = distp;
                    else
                    {    int dist_old = distx[LP];
                         if ( distp < dist_old )
                         {    nhood.push( LP, distp );
                              distx[LP] = distp;    }    }    }    }
          vec<int> ls;
          for ( auto x : nhood ) if ( x.first != L1 ) ls.push_back( x.first );
          UniqueSort(ls);
          vec<int> lsr;
          set<int> lsrx;
          for ( int j = 0; j < D2.To(y).isize( ); j++ )
          {    int L = D2.ITo(y,j);
               if ( BinMember( ls, L ) ) 
               {    lsr.push_back(L);
                    lsrx.insert(L);    }    }
          for ( int j = 0; j < lsr.isize( ); j++ )
          {    int z = to_left2[ lsr[j] ];
               for ( int m = 0; m < D2.To(z).isize( ); m++ )
               {    int L = D2.ITo(z,m);
                    if ( !BinMember( ls, L ) || Member( lsrx, L ) ) continue;
                    lsr.push_back(L), lsrx.insert(L);    }    }
          if ( lsr.empty( ) && y != x ) continue;
          UniqueSort(lsr);
          const int MAX_INTERMEDIATES = 100;
          if ( lsr.isize( ) > MAX_INTERMEDIATES ) continue;
          if (verbose)
          {    cout << "\ntype 2 linking from L" << L1 << " to L" << L2 << endl;
               cout << "ls = " << printSeq(ls) << endl;
               PRINT3( lsr.size( ), x, y );
               cout << "final intermediate lines = " << printSeq(lsr) << endl;    }
          nlinks++;
          vec<int> em, emr;
          for ( auto L : lsr )
          {    const vec<vec<vec<int>>>& X = dlines[L];
               em.append( Contents(X) );
               for ( int j = 1; j < X.isize( ); j += 2 )
               {    if ( X[j].solo( ) && X[j][0].empty( ) )
                    {    int v = to_right[ X[j-1][0][0] ];
                         em.push_back( D.IFrom(v,0) );    }    }    }
          UniqueSort(em); 
          for ( auto d : em ) emr.push_back( dinv[d] );
          int E = D.E( ), n = em.size( );
          digraphE<vec<int>> M( digraphE<vec<int>>::COMPLETE_SUBGRAPH_EDGES,
               D, em, to_left, to_right );
          D.AppendWithUpdate( M, to_left, to_right );
          digraphE<vec<int>> RM( digraphE<vec<int>>::COMPLETE_SUBGRAPH_EDGES,
               D, emr, to_left, to_right );
          D.AppendWithUpdate( RM, to_left, to_right );
          dinv.resize( E + 2*n );
          for ( int i = 0; i < n; i++ )
          {    dinv[ E + i ] = E + n + i, dinv[ E + n + i ] = E + i;    }
          int rd1 = dinv[d1], rd2 = dinv[d2];
          if ( n == 0 )
          {    int N = D.N( );
               D.AddVertices(2);
               D.GiveEdgeNewToVxWithUpdate( d1, to_right[d1], N, to_right );
               D.GiveEdgeNewFromVxWithUpdate( d2, to_left[d2], N, to_left );
               D.GiveEdgeNewToVxWithUpdate( rd2, to_right[rd2], N+1, to_right );
               D.GiveEdgeNewFromVxWithUpdate( rd1, to_left[rd1], N+1, to_left );    }
          else
          {    int v1 = -1, v2 = -1, rv1 = -1, rv2 = -1;
               for ( int i = 0; i < em.isize( ); i++ )
               {    int d = em[i];
                    if ( to_left[d] == to_right[d1] )
                    {    v1 = to_left[E+i];
                         break;    }    }
               for ( int i = 0; i < em.isize( ); i++ )
               {    int d = em[i];
                    if ( to_right[d] == to_left[d2] )
                    {    v2 = to_right[E+i];
                         break;    }    }
               for ( int i = 0; i < emr.isize( ); i++ )
               {    int d = emr[i];
                    if ( to_left[d] == to_right[rd2] )
                    {    rv2 = to_left[E+n+i];
                         break;    }    }
               for ( int i = 0; i < emr.isize( ); i++ )
               {    int d = emr[i];
                    if ( to_right[d] == to_left[rd1] )
                    {    rv1 = to_right[E+n+i];
                         break;    }    }
               D.GiveEdgeNewToVxWithUpdate( d1, to_right[d1], v1, to_right );
               D.GiveEdgeNewFromVxWithUpdate( d2, to_left[d2], v2, to_left );
               D.GiveEdgeNewToVxWithUpdate( rd2, to_right[rd2], rv2, to_right );
               D.GiveEdgeNewFromVxWithUpdate( rd1, to_left[rd1], rv1, to_left );    
                    }    }
     cout << Date( ) << ": made " << nlinks << " links" << endl;
     RemoveUnneededVertices( D, dinv );
     CleanupCore( D, dinv );
     Validate( hb, inv, D, dinv );    }

// Look for a line have two lines entering on the left and two exiting on the
// right, so that the lefts are inv to the rights.  Make two copies of the line
// and separate.

void BigPullInv( const HyperBasevectorX& hb, const vec<int>& inv,
     digraphE<vec<int>>& D, vec<int>& dinv, const int MIN_LEN )
{    cout << Date( ) << ": starting big pull inv" << endl;
     const int MAX_EDGES = 100;
     vec<vec<vec<vec<int>>>> dlines;
     while(1)
     {    int lpulls = 0;
          FindLines( D, dinv, dlines, MAX_CELL_PATHS, MAX_CELL_DEPTH);
          vec<int> llens, to_left, to_right, tol;
          GetLineLengths( hb, D, dlines, llens );
          D.ToLeft(to_left), D.ToRight(to_right);
          MakeTol( D, dlines, tol );
          vec<Bool> touched( dlines.size( ), False );
          cout << Date( ) << ": start loop" << endl; // XXXXXXXXXXXXXXXXXXXXXXXXXXXX
          int E0 = D.E( );
          for ( int m = 0; m < dlines.isize( ); m++ )
          {    const vec<vec<vec<int>>>& L = dlines[m];
               int u1 = L.front( )[0][0], u2 = L.back( )[0][0];
               if ( dinv[u1] != u2 ) continue;
               int v = to_left[u1], w = to_right[u2];
               if ( D.To(v).size( ) != 2 || !D.From(v).solo( ) ) continue;
               int d1 = D.ITo(v,0), d2 = D.ITo(v,1);
               int f1 = dinv[d1], f2 = dinv[d2];
               if ( u1 >= E0 || u2 >= E0 || d1 >= E0 || d2 >= E0 || f1 >= E0 
                    || f2 >= E0 )
               {    continue;    }
               int a1 = tol[d1], b1 = tol[d2], a2 = tol[f1], b2 = tol[f2];
               if ( a1 < 0 || b1 < 0 || a2 < 0 || b2 < 0 ) continue;
               if ( llens[a1] < MIN_LEN || llens[b1] < MIN_LEN ) continue;
               if ( !IsUnique( vec<int>{ m, a1, a2, b1, b2 } ) ) continue;
               vec<int> em = Contents( dlines[m] ), mto_left, mto_right;
               for ( int j = 1; j < L.isize( ); j += 2 )
               {    if ( L[j].solo( ) && L[j][0].empty( ) )
                    {    int v = to_right[ L[j-1][0][0] ];
                         em.push_back( D.IFrom(v,0) );    }    }
               Sort(em);
               if ( em.isize( ) > MAX_EDGES ) continue;
               if ( touched[m] ) continue;
               if ( touched[a1] || touched[a2] || touched[b1] || touched[b2] ) 
                    continue;
               // cout << endl;
               // PRINT6( d1, d2, f1, f2, u1, u2 ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXX
               // PRINT5( m, a1, a2, b1, b2 ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
               touched[m] = True;
               touched[a1] = touched[a2] = touched[b1] = touched[b2] = True;
               digraphE<vec<int>> M( digraphE<vec<int>>::COMPLETE_SUBGRAPH_EDGES,
                    D, em, to_left, to_right );
               M.ToLeft(mto_left), M.ToRight(mto_right);
               int N = D.N( ), E = D.E( );
               D.AppendWithUpdate( M, to_left, to_right );
               int p1 = BinPosition( em, u1 ), p2 = BinPosition( em, u2 );
               D.GiveEdgeNewToVxWithUpdate(
                    d1, to_right[d1], N + mto_left[p1], to_right );
               D.GiveEdgeNewFromVxWithUpdate(
                    f2, to_left[f2], N + mto_right[p2], to_left );
               vec<int> eminv( em.size( ) );
               for ( int i = 0; i < em.isize( ); i++ )
                    eminv[i] = BinPosition( em, dinv[ em[i] ] );
               for ( int i = 0; i < em.isize( ); i++ )
               {    dinv[ em[i] ] = E + eminv[i];
                    dinv.push_back( em[ eminv[i] ] );    }
               // Validate( hb, inv, D, dinv ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
               lpulls++;    }
          RemoveUnneededVertices( D, dinv );    
          CleanupCore( D, dinv );
          cout << Date( ) << ": made " << lpulls << " big pull inv edits" << endl;
          if ( lpulls == 0 ) break;    }    }

// Remove lines that fail a diploid filter.

void DiploidFilter( const HyperBasevectorX& hb, digraphE<vec<int>>& D, 
     vec<int>& dinv )
{
     // Remove lines that fail a diploid filter.

     const int MIN_DIP = 1000;
     const int EXT_DEPTH = 15;
     vec<vec<vec<vec<int>>>> dlines;
     FindLines( D, dinv, dlines, MAX_CELL_PATHS, MAX_CELL_DEPTH );
     vec<int> toll( D.E( ), -1 );
     for ( int i = 0; i < dlines.isize( ); i++ )
     for ( int j = 0; j < dlines[i].isize( ); j++ )
     for ( int k = 0; k < dlines[i][j].isize( ); k++ )
     for ( int l = 0; l < dlines[i][j][k].isize( ); l++ )
          toll[ dlines[i][j][k][l] ] = i;
     vec<int> dlens( D.E( ), 0 );
     for ( int d = 0; d < D.E( ); d++ )
     {    if ( D.O(d)[0] < 0 ) continue;
          for ( auto e : D.O(d) ) dlens[d] += hb.Kmers(e);    }
     vec<int> linv;
     LineInv( dlines, dinv, linv );
     vec<Bool> keepl( dlines.size( ), False );
     for ( int l = 0; l < dlines.isize( ); l++ )
     {    const vec<vec<vec<int>>>& L = dlines[l];
          int n = 0;
          for ( int j = 0; j < L.isize( ); j++ )
          {    if ( L[j].size( ) > 2 )
               {    n = 0;
                    continue;    }
               vec<int> lens( L[j].size( ) );
               for ( int k = 0; k < L[j].isize( ); k++ )
                    for ( auto d : L[j][k] ) lens[k] += dlens[d];
               Sort(lens);
               if ( lens.nonempty( ) ) n += Median(lens);
               if ( n >= MIN_DIP ) break;    }
          if ( n >= MIN_DIP ) 
               keepl[l] = keepl[ linv[l] ] = True;    }
     vec<int> keep( D.E( ), False );
     for ( int d = 0; d < D.E( ); d++ )
          if ( D.O(d)[0] < 0 || keepl[ toll[d] ] ) keep[d] = True;
     vec<int> to_left, to_right;
     D.ToLeft(to_left), D.ToRight(to_right);
     for ( int l = 0; l < dlines.isize( ); l++ )
     {    if ( !keepl[l] ) continue;
          const vec<vec<vec<int>>>& L = dlines[l];
          vec<int> ds = { L.front( )[0][0], L.back( )[0][0] };
          for ( int p = 0; p < EXT_DEPTH; p++ )
          {    int nds = ds.size( );
               for ( int i = 0; i < nds; i++ )
               {    int d = ds[i];
                    for ( int pass = 1; pass <= 2; pass++ )
                    {    int v = (pass == 1 ? to_left[d] : to_right[d]);
                         for (int j = 0; j < D.From(v).isize( ); j++)
                              ds.push_back( D.IFrom(v,j) );
                         for (int j = 0; j < D.To(v).isize( ); j++)
                              ds.push_back( D.ITo(v,j) );    }    }
               UniqueSort(ds);    }
          for ( auto d : ds ) keep[d] = True;    }
     vec<int> dels;
     for ( int d = 0; d < D.E( ); d++ ) if ( !keep[d] ) dels.push_back(d);
     D.DeleteEdges(dels);
     RemoveUnneededVertices( D, dinv );    
     CleanupCore( D, dinv );    }

// Remove loops that are identical to another path.

void RemoveSuperfluousLoops( const HyperBasevectorX& hb, 
     digraphE<vec<int>>& D, vec<int>& dinv )
{    cout << Date( ) << ": remove loops identical to another path" << endl;
     vec<Bool> deleted( D.E( ), False );
     vec<int> dels, to_left, to_right;
     D.ToLeft(to_left), D.ToRight(to_right);
     for ( int d = 0; d < D.E( ); d++ )
     {    if ( D.O(d)[0] < 0 ) continue;
          int v = to_left[d];
          if ( to_right[d] != v ) continue;
          int n = 0;
          for ( auto e : D.O(d) ) n += hb.Kmers(e);

          // Now we have edge d, which is a loop from v to v.

          vec<vec<int>> partials = { { } };
          int loops = 0;
          const int MAX_LOOPS = 1000;
          while( partials.nonempty( ) && ++loops <= MAX_LOOPS )
          {    vec<int> p = partials.back( ), x;
               partials.pop_back( );
               for ( auto f : p ) x.append( D.O(f) );
               if ( x == D.O(d) && to_right[ p.back( ) ] == v )
               {    deleted[d] = deleted[ dinv[d] ] = True;
                    dels.push_back( d, dinv[d] );
                    break;    }
               int vn = ( p.empty( ) ? v : to_right[ p.back( ) ] );
               for ( int j = 0; j < D.From(vn).isize( ); j++ )
               {    int f = D.IFrom( vn, j );
                    if ( D.O(f)[0] < 0 ) continue;
                    if ( f == d || deleted[f] ) continue;
                    if ( x.size( ) + D.O(f).size( ) > D.O(d).size( ) ) continue;
                    vec<int> p2 = p;
                    p2.push_back(f);
                    int nn = 0;
                    for ( auto e : p2 ) nn += hb.Kmers(e);
                    if ( nn > n ) continue;
                    partials.push_back(p2);    }    }    }
     UniqueSort(dels);
     cout << Date( ) << ": deleting " << dels.size( ) << " loops" << endl;
     D.DeleteEdges(dels);
     RemoveUnneededVertices( D, dinv );
     CleanupCore( D, dinv );    }

void RemoveVerySmallComponents( const HyperBasevectorX& hb, 
     digraphE<vec<int>>& D, vec<int>& dinv, const int MIN_COMP_SIZE )
{    vec<vec<int>> comp;
     D.ComponentsEFast(comp);
     vec<int> dels;
     for ( auto& x : comp )
     {    int64_t n = 0;
          for ( auto d : x )
          {    if ( D.O(d)[0] >= 0 ) for ( auto e : D.O(d) ) n += hb.Kmers(e);
               else if ( IsSequence( D.O(d) ) )
               {    int ltrim, rtrim;
                    basevector b;
                    GapToSeq( D.O(d), ltrim, rtrim, b );
                    n += b.isize( ) - hb.K( ) + 1;    }
               else if ( IsCell( D.O(d) ) )
               {    cell c;
                    c.CellDecode( D.O(d) );
                    const digraphE<vec<int>>& C = c.G( );
                    for ( int l = 0; l < C.E( ); l++ )
                    {    const vec<int>& z = C.O(l);
                         if ( z[0] < 0 ) continue;
                         for ( auto e : z ) n += hb.Kmers(e);    }    }    }
          if ( n < MIN_COMP_SIZE )
          {    for ( auto d : x ) dels.push_back( d, dinv[d] );    }    }
     D.DeleteEdges(dels);
     RemoveUnneededVertices( D, dinv );    
     CleanupCore( D, dinv );    }

// Delete pair gaps inside a cell if there is a gapless path through it

void DeletePairGapsInCells( digraphE<vec<int>>& D, vec<int>& dinv,
     vec<vec<vec<vec<int>>>> & dlines, vec<int> & dels, Bool verbose )
{
     vec<int> dels2;
     cout << Date( ) << ": deleting some gap edges" << endl;
     for ( int i = 0; i < dlines.isize( ); i++ )
     {    const vec<vec<vec<int>>>& L = dlines[i];
          for ( int j = 1; j < L.isize( ); j += 2 )
          {    Bool gapless = False;
               for ( int k = 0; k < L[j].isize( ); k++ )
               {    Bool gap = False;
                    for ( int l = 0; l < L[j][k].isize( ); l++ )
                    {    int d = L[j][k][l];
                         if ( IsPairGap( D.O(d) ) || 
                              IsBarcodeOnlyGap( D.O(d) ) )
                              gap = True;    }
                    if ( !gap ) gapless = True;    }
               if (gapless)
               {    for ( int k = 0; k < L[j].isize( ); k++ )
                    for ( int l = 0; l < L[j][k].isize( ); l++ )
                         if ( IsPairGap( D.O( L[j][k][l] ) ) )
                              dels2.push_back( L[j][k][l] );    }    }    }
     dels.append( dels2 );
     if ( verbose ) {
          cout << Date ( ) << ": deleted " << dels.size() << " pair gaps"
               << " in cells" << endl;  }    }

// THESE ARE SUPPORT FUNCTIONS THAT DELETE PAIR GAPS
// ALL COLLECTED IN ONE PLACE

// Rank paths through a cell based on "badness" = ( num gaps, -length )
// so low gaps, long edges => low badness
// keep all edges that are in paths with zero gaps
// if there is at least one gap,
//   keep paths with equal badness
//   if there are two paths, one that is more bad, delete the more bad one

void DeleteSomeRedundantGaps(
     const HyperBasevectorX& hb, digraphE<vec<int>>& D, vec<int>& dinv,
     Bool verbose )
{    vec<vec<vec<vec<int>>>> dlines;
     FindLines( D, dinv, dlines, MAX_CELL_PATHS, MAX_CELL_DEPTH );
     vec<int> dels;
     for ( int i = 0; i < dlines.isize( ); i++ )
     {    const vec<vec<vec<int>>>& L = dlines[i];
          for ( int j = 1; j < L.isize( ); j += 2 )
          {    vec<vec<int>> M = L[j];
               vec<pair<int,int>> badness( M.size( ), make_pair(0,0) );
               for ( int k = 0; k < M.isize( ); k++ )
               {    for ( int l = 0; l < M[k].isize( ); l++ )
                    {    int d = M[k][l];
                         if ( IsPairGap( D.O(d) ) ) badness[k].first++;
                         else if ( D.O(d)[0] >= 0 )
                         {    for ( auto e : D.O(d) )
                              {   badness[k].second 
                                        -= hb.Kmers(e);    }    }    }    }
               SortSync( badness, M );
               if ( M.nonempty( ) )
               {    vec<int> keep = M[0];
                    for ( int j = 1; j < M.isize( ); j++ )
                    {    if ( badness[j] > badness[j-1] 
                              && badness[j].first >= 1 )
                         {    UniqueSort(keep);
                              for ( int k = j; k < M.isize( ); k++ )
                              for ( auto d : M[k] )
                              {    if ( !BinMember( keep, d ) )
                                        dels.push_back( d, dinv[d] );    }
                              break;    }
                         else keep.append( M[j] );    }    }    }    }
     
     if ( verbose ) {
          cout << Date( ) << ": deleted " << dels.size() 
               << " redundant pair gaps" << endl;  }
     D.DeleteEdges(dels);
     RemoveUnneededVertices( D, dinv );
     CleanupCore( D, dinv );    }

// If we have two or more pair gaps between the same vertices,
// i.e., they are parallel, then only retain one

void DeleteParallelGaps( digraphE<vec<int>>& D, vec<int>& dinv, 
     Bool verbose )
{    vec<Bool> deleted( D.E( ), False );
     vec<int> dels, to_left, to_right;
     D.ToLeft(to_left), D.ToRight(to_right);
     for ( int d = 0; d < D.E( ); d++ )
     {    if ( !IsPairGap( D.O(d) ) || deleted[d] ) continue;
          int v = to_left[d], w = to_right[d];
          for ( int j = 0; j < D.From(v).isize( ); j++ )
          {    if ( D.From(v)[j] != w ) continue;
               int d2 = D.IFrom(v,j);
               if ( d2 == d || deleted[d2] ) continue;
               if ( !IsPairGap( D.O(d2) ) ) continue;
               deleted[d2] = deleted[ dinv[d2] ] = True;
               dels.push_back( d2, dinv[d2] );    }    }
     
     if ( verbose ) {
          cout << Date( ) << ": deleted " << dels.size() 
               << " parallel gaps" << endl;  }
     D.DeleteEdges(dels);
     RemoveUnneededVertices( D, dinv );
     CleanupCore( D, dinv );    }

// Delete pair gaps v--->w where you can reach w from v
// by exploring up to depth MAX_DEPTH in the fwd direction
// from v.

void ZapPairGaps( digraphE<vec<int>>& D, vec<int>& dinv,
     Bool verbose )
{    const int MAX_DEPTH = 6;
     const int NPASSES = 2;
     for ( int pass = 1; pass <= NPASSES; pass++ )
     {    vec<int> dels;
          #pragma omp parallel for
          for ( int v = 0; v < D.N( ); v++ )
          for ( int j = 0; j < D.From(v).isize( ); j++ )
          {    int d = D.IFrom(v,j);
               if ( !IsPairGap( D.O(d) ) ) continue;
               int w = D.From(v)[j];
               if ( v == w ) continue;
               vec<int> vs = {v};
               for ( int p = 0; p < MAX_DEPTH; p++ )
               {    int nvs = vs.size( );
                    for ( int h = 0; h < nvs; h++ )
                    {    int z = vs[h];
                         for ( int l = 0; l < D.From(z).isize( ); l++ )
                         {    int g = D.IFrom(z,l);
                              if ( D.O(g)[0] >= 0 ) 
                                   vs.push_back( D.From(z)[l] );    }    }
                    UniqueSort(vs);
                    if ( BinMember( vs, w ) )
                    {    
                         #pragma omp critical
                         {    dels.push_back( d, dinv[d] );    }
                         break;    }    }    }
          if ( verbose ) {
               cout << Date ( ) << ": zapped " << dels.size() 
                    << " pair gaps" << endl; }
          D.DeleteEdges(dels);
          RemoveUnneededVertices( D, dinv );
          CleanupCore( D, dinv );    }    }

// If a pair gap and a sequence-ful edge emerge from the same vertex
// kill the pair gap

void KillPairGapsAtBranches( digraphE<vec<int>>& D, vec<int>& dinv,
     Bool verbose )
{    vec<int> dels;
     for ( int v = 0; v < D.N( ); v++ )
     {    for ( int j1 = 0; j1 < D.From(v).isize( ); j1++ )
          {    int d1 = D.IFrom(v,j1);
               if ( D.O(d1)[0] < 0 ) continue;
               for ( int j2 = 0; j2 < D.From(v).isize( ); j2++ )
               {    int d2 = D.IFrom(v,j2);
                    if ( IsPairGap( D.O(d2) ) )
                         dels.push_back( d2, dinv[d2] );    }    }    }
     if ( verbose ) {
          cout << Date( ) << ": killed " << dels.size() << " pair gaps"
               << " at branches" << endl;    }
     D.DeleteEdges(dels);
     RemoveUnneededVertices( D, dinv );
     CleanupCore( D, dinv );    }

/******************************************************************************/

void HangBeGone( const HyperBasevectorX& hb, digraphE<vec<int>>& D, vec<int>& dinv,
     const int MIN_RATIO2, const int MAX_KILL2, const int verbosity, const Bool NEW )
{    if ( verbosity >= 1 ) cout << Date( ) << ": removing hanging ends" << endl;
     vec<int> lens( D.E( ), 0 ), dfw;
     #pragma omp parallel for
     for ( int e = 0; e < D.E( ); e++ )
     {    if ( D.O(e)[0] < 0 ) continue;
          for ( int j = 0; j < D.O(e).isize( ); j++ )
               lens[e] += hb.Kmers( D.O(e)[j] );    }
     const int MIN_RATIO1 = 5;
     const int MAX_KILL1 = 4;
     if ( NEW )
          MaxDistanceToEndArr( D, lens, MAX_KILL2 * MIN_RATIO2, True, dfw );
     else
          DistancesToEndArr( D, lens, MAX_KILL2 * MIN_RATIO2, True, dfw );

     vec<int> dels;
     #pragma omp parallel for
     for ( int v = 0; v < D.N( ); v++ )
     {    for ( int j1 = 0; j1 < D.From(v).isize( ); j1++ )
          for ( int j2 = 0; j2 < D.From(v).isize( ); j2++ )
          {    if ( j2 == j1 ) continue;
               int w1 = D.From(v)[j1], w2 = D.From(v)[j2];
               int d1 = D.IFrom(v,j1), d2 = D.IFrom(v,j2);
               int len1 = lens[d1] + dfw[w1], len2 = lens[d2] + dfw[w2];
               if ( ( len1 <= MAX_KILL1 && len2 >= MIN_RATIO1 * len1
                    && len2 >= MAX_KILL1 )
                    || ( len1 <= MAX_KILL2 && len2 >= MIN_RATIO2 * len1
                    && len2 >= MAX_KILL2 ) )
               {
                    #pragma omp critical
                    {    dels.push_back( d1, dinv[d1] );    
                         if ( verbosity >= 2 )
                         {    cout << "deleting edges " << d1 << " and "
                                   << dinv[d1] << endl;    }    }    }    }    }
     D.DeleteEdges(dels);
     RemoveUnneededVertices( D, dinv );
     CleanupCore( D, dinv );    }

void RemoveRedundantCells( digraphE<vec<int>>& D, vec<int>& dinv )
{    vec<int> dels;
     vec<Bool> deleted( D.E( ), False );
     for ( int v = 0; v < D.N( ); v++ )
     for ( int j1 = 0; j1 < D.From(v).isize( ); j1++ )
     for ( int j2 = 0; j2 < D.From(v).isize( ); j2++ )
     {    if ( j2 == j1 ) continue;
          if ( D.From(v)[j1] != D.From(v)[j2] ) continue;
          int d1 = D.IFrom(v)[j1], d2 = D.IFrom(v)[j2];
          if ( deleted[d1] || deleted[d2] ) continue;
          if ( !IsCell( D.O(d1) ) || !IsCell( D.O(d2) ) ) continue;
          if ( D.O(d1) != D.O(d2) ) continue;
          dels.push_back( d2, dinv[d2] );
          deleted[d2] = deleted[ dinv[d2] ] = True;    }
     D.DeleteEdges(dels);
     RemoveUnneededVertices( D, dinv );
     CleanupCore( D, dinv );    }


// More cleanup functions below
//
//

// Are we around a cell or sequence gap?

Bool AroundGap( const digraphE<vec<int>> & D, const vec<int> & to_left,
     const vec<int> & to_right, const int & d )
{
     for ( int vpass = 1; vpass <= 2; vpass++ ) {
          const int v = (vpass == 1 ? to_left[d] : to_right[d] );
          for ( int pass = 1; pass <=2; pass++ ) {
               const auto & es = ( pass == 1 ? D.IFrom(v) : D.ITo(v) );
               for ( auto & d2 : es )
                    if ( IsCell( D.O(d2) ) || IsSequence( D.O(d2) ) )
                         return True;
          }
     }
     return False;
}

// Wherever there is a branch, with support >=10 to 0, liberate the weak
// branch, so long as it contains some unique kmers.
     
void KillZeroSupportBranches( const HyperBasevectorX & hb, const vec<int> & inv,
     digraphE<vec<int>> & D, vec<int> & dinv, const ReadPathVec & dpaths, 
     const VecULongVec & dpaths_index, vec<Bool> & to_delete, const Bool verbose )
{
     vec<int> to_left, to_right, mult;
     D.ToLeftParallel(to_left), D.ToRightParallel(to_right);
     ComputeMult( hb, D, mult );
     vec<pair<Bool,Bool>> liberate(D.E(), make_pair(False,False));
     for ( int v = 0; v < D.N( ); v++ )
     {    if ( !D.To(v).solo( ) || !D.From(v).duo( ) ) continue;
          int g = D.ITo(v,0), f1 = D.IFrom(v,0), f2 = D.IFrom(v,1);
          if ( D.O(g)[0] < 0 ) continue;
          
          if ( AroundGap(D, to_left, to_right, g ) ) continue;
          if ( AroundGap(D, to_left, to_right, f1 ) ) continue;
          if ( AroundGap(D, to_left, to_right, f2 ) ) continue;
          
          for ( int z = 0; z < 2; z++ )
          {    int d1 = D.IFrom(v,z), d2 = D.IFrom(v,1-z);
               if ( D.O(d1)[0] < 0 || D.O(d2)[0] < 0 ) continue;
               int rd1 = dinv[d1], rd2 = dinv[d2];
               Bool solo1 = False;
               for ( auto e : D.O(d1) ) if ( mult[e] == 1 ) solo1 = True;
               if ( !solo1 ) continue;
               int n1 = 0, n2 = 0;
               for ( int mpass = 1; mpass <= 2; mpass++ )
               for ( int pass = 1; pass <= 2; pass++ )
               {    int f;
                    if ( mpass == 1 ) f = ( pass == 1 ? d1 : dinv[d1] );
                    else f = ( pass == 1 ? d2 : dinv[d2] );
                    for ( auto & id : dpaths_index[f])
                    {    const ReadPath& p = dpaths[id];
                         for ( int l = 0; l < (int) p.size( ) - 1; l++ )
                         {    if ( pass == 1 )
                              {    if ( p[l] == g && p[l+1] == f )
                                   {    ( mpass == 1 ? n1 : n2 )++;
                                        break;    }    }
                              else
                              {    if ( p[l] == f && p[l+1] == dinv[g] )
                                   {   ( mpass == 1 ? n1 : n2 )++;
                                        break;    }    }    }    }    }
               if ( n1 > 0 ) continue;
               if ( n2 < 5 || n2 < 10 * n1 ) continue;
               liberate[d1].first = True;
               liberate[rd1].second = True;
               break;    }    }
     for ( int d = 0; d < D.E( ); d++ ) {
          if ( to_delete[d] ) continue;
          const Bool from = liberate[d].first, to = liberate[d].second;
          if ( from && to )
               to_delete[d] = True;
          else if ( from ) {
               if ( verbose )
                    cout << "liberating " << d << " on left" << endl;
               int N = D.N( );
               D.AddVertices(1);
               D.GiveEdgeNewFromVx( d, to_left[d], N );
          } else if ( to ) {
               if ( verbose )
                    cout << "liberating " << d << " on right" << endl;
               int N = D.N( );
               D.AddVertices(1);
               D.GiveEdgeNewToVx( d, to_right[d], N );
          }
     }
}

// Delete zero support edges of a graph based on the following criterion:
// if an edge between vertices v and w has zero weight, AND
// there is an alternate path with support between v and w, DELETE.
// support: < 0 for edges that shouldn't be included in alt paths
//        : = 0 for edges that are considered for deletion
//        : > 0 for edges that are to count towards alt paths
// edges that are to be deleted are marked so in to_del

template <class F>
void DeleteZeroSupportEdges( const digraphE<F> & G, const vec<int> & ginv,
     const vec<int> & to_left, const vec<int> & to_right,
     const vec<int> & support, vec<Bool> & to_del )
{
     ForceAssertEq( support.isize(), G.E() );
     ForceAssertEq( to_del.isize(), G.E() );

     const int MAX_DEPTH=10;
     vec<int> expl;
     #pragma omp parallel for schedule (dynamic,1) firstprivate(expl)
     for ( int d = 0; d < G.E(); d++ ) {
          if ( to_del[d] ) continue;
          if ( support[d] != 0 ) continue;

          const int v = to_left[d], target = to_right[d];
          
          expl.clear();
          set<int> all = {v};
          set<int> vset = {v};
          Bool alt = False;
          for ( int depth = 0; depth < MAX_DEPTH; depth++ ) {
               const int start = expl.size( );
               for ( auto & vp : vset ) {
                    for ( int i = 0; i != G.From(vp).isize(); i++ ) {
                         const int edge = G.IFrom( vp, i );
                         if ( edge == d ) continue;
                         if ( support[edge] > 0 ) expl.push_back( edge );
                    }
               }
               vset.clear( );
               
               if ( expl.isize() == start ) break;

               for ( int i = start; i != expl.isize( ); i++ ) {
                    const int w = to_right[expl[i]];
                    if ( w == target ) {
                         alt = True;
                         goto found;
                    }
                    if ( all.count( w ) ) continue;
                    vset.insert( w );
                    all.insert( w );
               }
          }
found:    const int rd = ginv[d];
          if ( alt && !(to_left[rd] == v && to_right[rd] == target) ) {
               #pragma omp critical
               {    to_del[d] = True;
                    to_del[rd] = True;  }
          }
         
     }
}

// If there is a pair gap between two vertices, and there is a path
// through the graph that contains at least one regular edge
// then delete the pair gap

void PairGapKiller( digraphE<vec<int>> & D, vec<int> & dinv, const int MAX_DEPTH )
{
     cout << Date( ) << ": killing pair gaps in multiple passes" << endl;
     vec<Bool> to_delete;
     vec<Bool> skip( D.E(), False );
     vec<Bool> noalt( D.E(), False );
     vec<int> to_right;
     D.ToRightParallel( to_right );
     int pass = 1;
     int dels = 0;
     // make multiple passes
     while ( True )
     {
          pass++;
          to_delete.clear();
          to_delete.resize( D.E(), False );
          skip.clear();
          skip.resize( D.E(), False );
          vec<int> expl;
          #pragma omp parallel for firstprivate(expl)
          for ( int v = 0; v < D.N(); v++ ) {
               for ( auto & d : D.IFrom(v) ) {
                    if ( !IsPairGap( D.O(d) ) ) continue;
                    if ( noalt[d] ) continue;

                    // can we delete this pair gap
                    
                    // if it's a loop then delete
                    const int target = to_right[d];
                    if ( target == v ) {
                         to_delete[d] = True;
                         break;
                    }
                    // stores all edges explored EXCEPT d
                    expl.clear();
                    set<int> all = {v};
                    set<int> vset = {v};
                    Bool alt = False;
                    for ( int depth = 0; depth < MAX_DEPTH; depth++ ) {
                         const int start = expl.size( );
                         for ( auto & vp : vset ) {
                              for ( int i = 0; i != D.From(vp).isize(); i++ ) {
                                   const int edge = D.IFrom( vp, i );
                                   if ( edge != d )
                                        expl.push_back( edge );
                              }
                         }
                         vset.clear( );
                         
                         if ( expl.isize() == start ) break;

                         for ( int i = start; i != expl.isize( ); i++ ) {
                              const int w = to_right[expl[i]];
                              if ( w == target ) {
                                   alt = True;
                                   goto found;
                              }
                              if ( all.count( w ) ) continue;
                              vset.insert( w );
                              all.insert( w );
                         }
                    }
found:              if ( alt ) {
                         for ( auto & d2 : expl )
                              skip[d2] = True;
                         to_delete[d] = True;
                         break;
                    } else
                         noalt[d] = True;
               }
          }

          
          Bool noop = True;
          #pragma omp parallel for
          for ( int d = 0; d < D.E(); d++ ) {
               const int rd = dinv[d];
               if ( rd < d ) continue;
               const Bool del = to_delete[d] || to_delete[rd];
               const Bool save = skip[d] || skip[rd];
               if ( save ) {
                    to_delete[d] = False;
                    to_delete[rd] = False;
               } else if ( del ) {
                    noop = False;
                    to_delete[d] = True;
                    to_delete[rd] = True;
                    dels += 1 + int( d != rd );
               }
          }
          if ( noop ) break;
          D.DeleteEdgesParallel( to_delete );
     }
     cout << Date( ) << ": deleted " << ToStringAddCommas( dels ) << " pair gaps"
          << " in " << pass << " passes" << endl;
     RemoveUnneededVertices( D, dinv );
     CleanupCore( D, dinv );
     RemoveUnneededVertices( D, dinv );
     CleanupCore( D, dinv );
}

void MarkEdgesWithZeroReads(digraphE<vec<int>> & D, vec<int> & dinv,
     const vec<int> to_left, const vec<int> to_right,
     const VecULongVec & dpaths_index, vec<Bool> & to_delete )
{
     cout << Date( ) << ": deleting zero support edges" << endl;
     // define support
     const int THRESHOLD = 1;
     vec<int> support( D.E( ), -1 );
     for ( int d = 0; d < D.E(); d++ ) {
          if ( D.O(d)[0] < 0 ) continue;
          const int s = dpaths_index[d].size() + dpaths_index[dinv[d]].size();
          if ( s >= THRESHOLD )
               support[d] = 1;
          else
               support[d] = 0;
     }
     DeleteZeroSupportEdges( D, dinv, to_left, to_right, support, to_delete );
     int numdels = 0;
     for ( auto & isdel : to_delete )
          if ( isdel ) numdels++;
     
     cout << Date( ) << ": marked " << ToStringAddCommas(numdels) 
          << " edges for deletion" << endl;
 
}

// Delete zero support bridges
void MarkWeakBridges( const digraphE<vec<int>> & D, const vec<int> & dinv,
     const vec<int> & to_left, const vec<int> & to_right,
     const VecULongVec & dpaths_index, vec<Bool> & to_delete )
{          
     cout << Date( ) << ": deleting zero support bridges" << endl;

     // zero support edges are 
     const int RELIABLE = 10;
     const int DELETE = 1;
     vec<int> support( D.E( ), -1 );
     #pragma omp parallel for schedule(dynamic,1)
     for ( int d = 0; d < D.E(); d++ ) {
          if ( IsPairGap(D.O(d)) )
               support[d] = 0;
          if ( D.O(d)[0] < 0 ) continue;
          const int s = dpaths_index[d].size() + dpaths_index[dinv[d]].size();
          if ( s <= DELETE )
               support[d] = 0;
          else {
               vec<int64_t> a;
               a.reserve(s);
               for ( auto & id : dpaths_index[d] )
                    a.push_back( id );
               for ( auto & id : dpaths_index[dinv[d]] )
                    a.push_back( id );
               UniqueSort(a);
               if ( a.size() >= RELIABLE )
                    support[d]=1;
               else
                    support[d]=-1; 
          }
     }
     #pragma omp parallel for
     for ( int d = 0; d < D.E( ); d++ ) {
          if ( support[d] != 0 ) continue;
          Bool alt = False;
          for ( int vpass = 1; vpass <= 2; vpass++ ) {
          for ( int pass = 1; pass <=2; pass++ ) {
               const int v = ( vpass == 1 ? to_left[d] : to_right[d] );
               const int num = ( pass == 1 ? D.From(v).size() : D.To(v).size() );
               alt=False;
               for ( int i = 0; i != num; i++ ) {
                    const int d2 = ( pass == 1 ? D.IFrom(v,i) : D.ITo(v,i) );
                    if ( d2 == d ) continue;
                    if ( to_delete[d2] ) continue;
                    if ( support[d2] == 1 )
                         alt = True;
               }
               if ( !alt )
                    goto skipout;
          }
          }
skipout:       if ( alt )
               to_delete[d] = True;
     }
     int numdels = 0;
     for ( auto & isdel : to_delete )
          if ( isdel ) numdels++;
     cout << Date( ) << ": marked " << ToStringAddCommas(numdels) 
          << " edges for deletion" << endl;
}

// Delete zero support lines that are hanging
void MarkDuplicateHangingLines( const HyperBasevectorX & hb, 
     const digraphE<vec<int>> & D, const vec<int> & dinv,
     const vec<vec<vec<vec<int>>>> & dlines, const VecULongVec & dpaths_index,
     vec<Bool> & to_delete )
{ 
     cout << Date( ) << ": building line graph" << endl;
     digraphE<int> LG;
     BuildLineGraph( D, dinv, dlines, LG );
     
     vec<int> llens;
     cout << Date( ) << ": computing line lengths" << endl;
     GetLineLengths( hb, D, dlines, llens );

     cout << Date( ) << ": computing line coverage" << endl;
     vec<double> lcov ( LG.E(), 0.0 );
     #pragma omp parallel for
     for ( int l = 0; l < LG.E(); l++ ) {
          const auto & L = dlines[l];
          
          double supp = 0;
          for ( int j = 0; j < L.isize( ); j++ )
               for ( auto & path : L[j] )
                    for ( auto & d : path ) {
                         supp += dpaths_index[d].size();
                         supp += dpaths_index[dinv[d]].size();
                    }
          supp *= 140;
          supp /= Max( 1, llens[l] );
          lcov[l] = supp;
     }
     cout << Date( ) << ": computing base edge content" << endl;
     vec<set<int>> tol( hb.E() );
     for ( int l = 0; l < LG.E(); l++ ) {
          const auto & L = dlines[l];
          for ( int j = 0; j < L.isize( ); j++ )
               for ( auto & path : L[j] )
                    for ( auto & d : path ) {
                         if ( D.O(d)[0] < 0 ) continue;
                         for ( auto & e : D.O(d) )
                              tol[e].insert( l );
                    }
     }
     
     const double MIN_LCOV = 1;
     const double MAX_SHARE = 0.95;
     const int MIN_LEN = 5000;
     vec<float>shared;
     vec<int>ldels;
     int lowcovlines=0;
     #pragma omp parallel for schedule(dynamic,1) firstprivate(shared)
     for ( int v = 0; v < LG.N(); v++ ) {
          if ( LG.To(v).nonempty() ) continue;
          if ( LG.From(v).size() != 1 ) continue;
          const int l = LG.IFrom(v,0);
          if ( lcov[l] > MIN_LCOV ) continue;
          if ( llens[l] > MIN_LEN ) continue;

          lowcovlines++;
          shared.clear();
          shared.resize( LG.E(), 0 );
          const auto & L = dlines[l];
          set<int> eset;
          for ( int j = 0; j < L.isize( ); j++ ) {
               eset.clear();
               for ( auto & path : L[j] )
               for ( auto & d : path ) {
                    if ( D.O(d)[0] < 0 ) continue;
                    for ( auto & e : D.O(d) ) {
                         eset.insert(e);
                    }
               }
          }
          int numbase = eset.size();
          for ( auto & e : eset )
               for ( auto & l2 : tol[e] ) {
                    if ( l2 == l ) continue;
                    shared[l2]+=1;
               }

          float sharefrac = float(Max(shared))/float(Max(1,numbase));
          //PRINT5( l, llens[l], lcov[l], numbase, sharefrac );
          ForceAssertLe( sharefrac, 1.0 );
          if ( sharefrac > MAX_SHARE ) {
               #pragma omp critical
               {    ldels.push_back( l );    }
          }
     }
     cout << Date( ) << ": found " << lowcovlines << " low coverage hanging lines"
          << endl;

     for ( auto & l : ldels ) {
          const auto & L = dlines[l];
          for ( int j = 0; j < L.isize( ); j++ )
          for ( auto & path : L[j] )
          for ( auto & d : path ) {
               to_delete[d] = True;
               to_delete[dinv[d]] = True;
          }
     }
}
 

// If a bc gap is in the canonical topology between two edges
// and there is sufficient overlap between the flanks then
// just glue them together

void DeleteRedundantBCGaps( const HyperBasevectorX & hb, 
     digraphE<vec<int>> & D, vec<int> & dinv )
{
     cout << Date( ) << ": delete unneeded bc gaps" << endl;
     vec<int> to_left, to_right;
     D.ToLeftParallel( to_left );
     D.ToRightParallel( to_right );
     vec<triple<int,int8_t,int>> edp;
     vec<vec<int>> match, score;
     vec<vec<Bool>> backtrack;
     const int MIN_GLUE = 100;
     vec< triple< pair<int,int>, pair<int,int>, int > > M;
     vec<int> dels;
     const int initE = D.E();
     vec<Bool> touched (D.N(), False);
     #pragma omp parallel for ordered
     for ( int d = 0; d < initE; d++ ) {
          if ( !IsBarcodeOnlyGap( D.O(d) ) ) continue;
          if( dinv[d] < d ) continue;
          const int v = to_left[d], w = to_right[d];
          if ( D.To(v).size() != 1 || D.From(w).size() != 1 )
               continue;
          const int l = D.ITo(v,0), r = D.IFrom(w,0);
          const auto & vl = D.O(l), & vr = D.O(r);
          if ( vl[0] < 0  || vr[0] < 0 )
               continue;
          edp.clear();
          for ( int i = 0; i < D.O(l).isize(); i++ )
               edp.push( vl[i], 0, i );
          for ( int i = 0; i < D.O(r).isize(); i++ )
               edp.push( vr[i], 1, i );
          Sort( edp );

          match.clear();
          score.clear();
          backtrack.clear();
          match.resize( vl.size(), vec<int>( vr.size(), 0 ) );
          score.resize( vl.size()+1, vec<int>( vr.size()+1, 0 ) );
          backtrack.resize( vl.size()+1, vec<Bool>( vr.size()+1, False ) );
          for ( int i = 0; i < edp.isize(); i++ ) {
               int j = i+1;
               for (; j < edp.isize(); j++ )
                    if ( edp[i].first != edp[j].first )
                         break;
               if ( j-i == 1 ) continue;
               vec<vec<int>> x(2);
               for ( int k = i; k < j; k++ )
                    x[int(edp[k].second)].push_back( edp[k].third );
               for ( auto & x1 : x[0] )
               for ( auto & y1 : x[1] )
                    match[x1][y1] = hb.Kmers( edp[i].first );
               i = j-1;
          }
          
          for ( int i1 = 0; i1 < vl.isize(); i1++ ) {
          for ( int i2 = 0; i2 < vr.isize(); i2++ ) {
               if ( match[i1][i2] > 0 ) {
                    score[i1+1][i2+1] = score[i1][i2] + match[i1][i2];
                    backtrack[i1+1][i2+1] = True;
               }
          }
          }
          pair<int,int> best(-1,-1);
          int max = 0;
          for ( int i1 = 0; i1 < score.isize(); i1++ ) {
          for ( int i2 = 0; i2 < score[i1].isize(); i2++ ) {
               if ( score[i1][i2] > max ) {
                    max = score[i1][i2];
                    best = make_pair( i1, i2 );
               }
          }    }
          if ( max < MIN_GLUE ) continue;

          int i1 = best.first, i2 = best.second;
          while ( i1 >= 0 && i2 >= 0 && backtrack[i1][i2] ) {
               i1--;
               i2--;
          }
          const int len = best.first - i1;

          // check for seq gap
          Bool invalid = False;
          for ( int pass = 1; pass <= 2; pass++ ) {
               const int vert = ( pass == 1 ? to_right[r] : to_left[l] );
               const auto & edges = ( pass == 1 ? D.IFrom(vert) : D.ITo(vert) );
               if ( edges.size() == 1 ) {
                    const int g = edges[0];
                    if ( IsSequence(D.O(g)) ) {
                         int ltrim, rtrim;
                         basevector seq;
                         GapToSeq( D.O(g), ltrim, rtrim, seq );
                         int trim = ( pass == 1 ? ltrim : rtrim );
                         if ( trim == 0 ) continue;
                         invalid = True;
                    }
               }
          }

          if ( invalid )
               continue;
          
          vec<int> alledges({l, dinv[l], r, dinv[r], d, dinv[d]});
          if ( !IsUnique( alledges ) ) continue;
          const int lmost = to_left[l], rmost = to_right[r];
          const int ilmost = to_left[dinv[r]], irmost = to_right[dinv[l]];
          vec<int> affected({lmost, rmost, ilmost, irmost});
          if ( !IsUnique( affected ) ) continue;
          
          #pragma omp ordered
          {
               Bool skip = False;
               for ( auto & vertex : affected ) {
                    for ( auto & cvert : D.To(vertex) )
                         if ( touched[cvert] )
                              skip = True;
                    for ( auto & cvert : D.From(vertex) )
                         if ( touched[cvert] )
                              skip = True;
               }
               if ( !skip ) {
                    for ( auto & vertex : affected ) {
                         for ( auto & cvert : D.To(vertex) )
                              touched[cvert] = True;
                         for ( auto & cvert : D.From(vertex) )
                              touched[cvert] = True;
                    }    
                    // Make the edit
                    
                    dels.append( alledges );
                    vec<int> common, rcommon;
                    const int llen = vl.size(), rlen = vr.size();
                    const int N = D.N();
                    int E = D.E();
                    common.SetToSubOf( D.O(l), i1, len );
                    rcommon.SetToSubOf( D.O(dinv[l]), llen-i1-len, len );
                    D.AddVertices(4);
                    D.AddEdgeWithUpdate( N, N+1, common, to_left, to_right );
                    D.AddEdgeWithUpdate( N+2, N+3, rcommon, to_left, to_right );
                    dinv.push_back( E+1 );
                    dinv.push_back( E );
                    E += 2;

                    if ( i1 > 0 ) {
                         vec<int> l1, rl1;
                         l1.SetToSubOf( D.O(l), 0, i1 );
                         rl1.SetToSubOf( D.O(dinv[l]), llen-i1, i1 );
                         ForceAssert(l1.nonempty());
                         ForceAssert(rl1.nonempty());
                         D.AddEdgeWithUpdate( lmost, N, l1, to_left, to_right );
                         D.AddEdgeWithUpdate( N+3, irmost, rl1, to_left, to_right );
                         dinv.push_back( E+1 );
                         dinv.push_back( E );
                         E += 2;
                    } else {
                         D.TransferEdgesWithUpdate( lmost, N, to_left,
                                                            to_right, False );
                         D.TransferEdgesWithUpdate( irmost,N+3, to_left,
                                                            to_right, False );
                    }
                    
                    if ( i1+len < llen ) {
                         vec<int> l3, rl3;
                         l3.SetToSubOf( D.O(l), i1+len, llen-i1-len );
                         rl3.SetToSubOf( D.O(dinv[l]), 0, llen-i1-len );
                         ForceAssert(l3.nonempty());
                         ForceAssert(rl3.nonempty());
                         const int V = D.N();
                         D.AddVertices(2);
                         D.AddEdgeWithUpdate( N+1, V, l3, to_left, to_right );
                         D.AddEdgeWithUpdate( V+1, N+2, rl3, to_left, to_right );
                         dinv.push_back( E+1 );
                         dinv.push_back( E );
                         E += 2;
                    }

                    if ( i2 > 0 ) {
                         vec<int> r1, rr1;
                         r1.SetToSubOf( D.O(r), 0, i2 );
                         rr1.SetToSubOf( D.O(dinv[r]), rlen-i2, i2 );
                         ForceAssert( r1.nonempty() );
                         ForceAssert( rr1.nonempty() );
                         const int V = D.N();
                         D.AddVertices(2);
                         D.AddEdgeWithUpdate( V, N, r1, to_left, to_right );
                         D.AddEdgeWithUpdate( N+3, V+1, rr1, to_left, to_right );
                         dinv.push_back( E+1 );
                         dinv.push_back( E );
                         E += 2;
                    }

                    if ( i2+len < rlen ) {
                         vec<int> r3, rr3;
                         r3.SetToSubOf( D.O(r), i2+len, rlen-i2-len );
                         rr3.SetToSubOf( D.O(dinv[r]), 0, rlen-i2-len );
                         ForceAssert( r3.nonempty() );
                         ForceAssert( rr3.nonempty() );
                         D.AddEdgeWithUpdate( N+1, rmost, r3, to_left, to_right );
                         D.AddEdgeWithUpdate( ilmost, N+2, rr3, to_left, to_right );
                         dinv.push_back( E+1 );
                         dinv.push_back( E );
                         E += 2;
                    } else {
                         D.TransferEdgesWithUpdate( rmost, N+1, to_left,
                                                       to_right, False );
                         D.TransferEdgesWithUpdate( ilmost, N+2, to_left,
                                                       to_right, False );
                    }
               }
          } // end ordered block 
     }
     cout << Date( ) << ": deleted " << dels.size()/3 << " bc gaps" << endl;
     D.DeleteEdgesParallel(dels);
     RemoveUnneededVertices( D, dinv );
     CleanupCore( D, dinv );
}

// Mark inversion artifacts for deletion

void MarkInversionArtifacts( const HyperBasevectorX & hb, const vec<int> & inv,
     digraphE<vec<int>> & D, vec<int> & dinv, const vec<vec<vec<vec<int>>>> & dlines,
     const VecULongVec & dpaths_index, vec<Bool> & to_delete,
     const Bool verbose )
{
     cout << Date( ) << ": build data structures" << endl;
     digraphE<int> LG;
     BuildLineGraph( D, dinv, dlines, LG );
     
     vec<int> lg_to_left, lg_to_right;
     LG.ToLeft( lg_to_left );
     LG.ToRight( lg_to_right );
     vec<vec<int>> tpaths, bpaths;
     vec<Bool> ldel( dlines.size(), False );
     const int MIN_READS = 5;

     cout << Date( ) << ": finding inversions" << endl;
     #pragma omp parallel for firstprivate(tpaths,bpaths)
     for ( int lv = 0; lv < LG.N(); lv++ ) {
          if ( LG.From(lv).size() != 2 || LG.To(lv).size() != 1 )
               continue;
          const int top = LG.IFrom(lv,0), bottom = LG.IFrom(lv,1);
          const int rlv = lg_to_right[top];
          if ( lg_to_right[bottom] != rlv )
               continue;
          if ( LG.From(rlv).size() != 1 )
               continue;
          
          // in the canonical topology
          // --->o---line bubble---->o--->
          
          const int nblocks = dlines[top].size();
          if ( dlines[bottom].isize() != nblocks )
               continue;
          const int lflank = dlines[LG.ITo(lv,0)].back()[0][0];
          const int rflank = dlines[LG.IFrom(rlv,0)][0][0][0];
          Bool inversion = True;
          for ( int bi = 0; bi < nblocks; bi++ ) {
               if (dlines[top][bi].size() != dlines[bottom][nblocks-1-bi].size()) {
                    inversion = False;
                    break;
               }
               const int npaths = dlines[top][bi].size();    
               
               /*
               // no gaps in cells
               Bool gaps = False;
               for ( int pass = 1; pass <= 2; pass++ ) {
                    const auto & paths = ( pass == 1 ? dlines[top][bi] : 
                                             dlines[bottom][nblocks-1-bi] );
                    for ( int pi = 0; pi < npaths; pi++ ) {
                         const auto & path = paths[pi];
                         for ( auto & d : path )
                              if ( D.O(d)[0] < 0 )
                                   gaps = True;
                    }
               }
               if ( gaps ) {
                    inversion = False;
                    break;
               }
               */

               tpaths.clear();
               bpaths.clear();
               tpaths.resize( npaths );
               bpaths.resize( npaths );
               
               for ( int pass = 1; pass <= 2; pass++ ) {
                    for ( int pi = 0; pi < npaths; pi++ ) {
                         const vec<int> & in = ( pass == 1 ? dlines[top][bi][pi] : 
                              dlines[bottom][nblocks-1-bi][pi] );
                         auto & out = ( pass == 1 ? tpaths[pi] : bpaths[pi] );
                         for ( auto & d : in )
                              out.append( (D.O(d)[0] < 0 ?
                                   vec<int>{D.O(d)[0]} : D.O(d)) );
                    }
               }

               for ( auto & path : bpaths ) {
                    for ( int pi = 0; pi < path.isize(); pi++ ) 
                         path[pi] = (path[pi] >= 0 ? inv[path[pi]] : path[pi] );
                    path.ReverseMe();
               }
               Sort( tpaths );
               Sort( bpaths );
               
               if ( tpaths != bpaths ) {
                    inversion = False;
                    break;
               }
          }
          if ( !inversion )
               continue;
          
          // get the super edges that are terminal on all lines at the
          // left and right line vertex
          const int ltop = dlines[top][0][0][0];
          const int rtop = dlines[top].back()[0][0];
          const int lbot = dlines[bottom][0][0][0];
          const int rbot = dlines[bottom].back()[0][0];
          vec<vec<int>> esets;
          esets.push_back({lflank, ltop, lbot});
          esets.push_back({rflank, rtop, rbot});
          vec<vec<int64_t>> pids;
          vec<int> choices(2, -1);
          for ( int esi = 0; esi < 2; esi++ ) {
               const auto & es = esets[esi];
               pids.clear();
               pids.resize(3);
               for ( int i = 0; i < 3; i++ ) {
                    for ( auto & id : dpaths_index[es[i]] )
                         pids[i].push_back( id/2 );
                    for ( auto & id : dpaths_index[dinv[es[i]]] )
                         pids[i].push_back( id/2 );
                    UniqueSort( pids[i] );
               }
               vec<double> score(2,0);
               score[0] = Max((double)MeetSize( pids[0], pids[1] ), .1);
               score[1] = Max((double)MeetSize( pids[0], pids[2] ), .1);
               if ( Min(score) < MIN_READS && Max(score) >= 2 * MIN_READS )
                    choices[esi] = ( score[0] > score[1] ? bottom : top );
          }
          if ( choices[0] == choices[1] && choices[0] != -1 ) {
               ldel[choices[0]] = True; 
               if ( verbose ) {
                    #pragma omp critical
                    {
                         cout << "--" << endl;
                         PRINT2(top,bottom);
                         PRINT2(lflank,rflank);
                         PRINT2( ltop, lbot );
                         cout << "deleting " << choices[0] << endl;
                    }
               }
          }
     }
     cout << Date( ) << ": marking edges for deletion" << endl;
     for ( int l = 0; l < LG.E(); l++ ) {
          if ( ldel[l] ) {
               for ( auto & blocks : dlines[l] ) {
               for ( auto & paths : blocks ) {
               for ( auto & d : paths ) {
                    to_delete[d] = True;
                    to_delete[dinv[d]] = True;
               }    }    }
          }
     }
}

// Compute the number of read-pair inserts that go across pairs of lines
// stored in a map<pair<int,int>,int> data structure
// (line1, line2) -> # of read pair inserts

void GetInsertsConnectingLines( const vec<int> & dinv, const vec<int> & linv,
     const ReadPathVec & dpaths, const VecULongVec & dpaths_index, const vec<int> & dtol,
     map<pair<int,int>,int> & lsupp )
{
     cout << Date( ) << ": computing insert support" << endl;
     vec<int> insert;
     #pragma omp parallel for firstprivate(insert) schedule(dynamic,1000)
     for ( size_t id = 0; id < dpaths.size(); id +=2  ) {
          if (dpaths[id].size() + dpaths[id+1].size() < 2)
               continue;
          insert.clear();
          for ( auto & x : dpaths[id] ) {
               const int l = dtol[x];
               if ( insert.nonempty() && insert.back() == l )
                    continue;
               insert.push_back( l );
          }
          int s = dpaths[id+1].size();
          for ( int i = 0; i < s; i++ ) {
               const int l = dtol[dinv[dpaths[id+1][s-i-1]]];
               if ( insert.nonempty() && insert.back() == l )
                    continue;
               insert.push_back( l );
          }
          if ( insert.size() < 2 ) continue; 
          for ( int pass = 1; pass <= 2; pass++ ) {
               for ( int i = 0; i < insert.isize()-1; i++ ) {
                    const int l1 = (pass == 1 ? insert[i] : insert[i+1]);
                    const int l2 = (pass == 1 ? linv[insert[i+1]] : 
                                                linv[insert[i]]);
                    const pair<int,int> p(l1,l2);
                    //if ( l1 < 0 || l2 < 0 )
                    //     continue;
                    #pragma omp critical
                    {    lsupp[p]++;    }
               }
          }
     }
}

// Pull apart 2-in/1-middle/2-out topologies using read pairs

void ReadPairPullApart( digraphE<vec<int>> & D, vec<int> & dinv,
     const ReadPathVec & dpaths, const VecULongVec & dpaths_index )
{
     vec<int> to_left, to_right;
     D.ToLeftParallel( to_left ), D.ToRightParallel( to_right );
     vec<quad<int,int,int,int>> pull;
     #pragma omp parallel for
     for ( int d = 0; d < D.E(); d++ ) {
          if ( dinv[d] >= d ) continue;
          const int l = to_left[d], r = to_right[d];
          if ( D.From(l).size() != 1 ) continue;
          if ( D.To(r).size() != 1 ) continue;
          if ( D.From(r).size() != 2 ) continue;
          if ( D.To(l).size() != 2 ) continue;
          
          Bool skip = False;
          for ( auto & x : D.ITo(l) ) {
               if ( dpaths_index[x].size() + dpaths_index[dinv[x]].size() < 20 )
                    skip = True;
          }
          for ( auto & x : D.IFrom(r) ) {
               if ( dpaths_index[x].size() + dpaths_index[dinv[x]].size() < 20 )
                    skip = True;
          }
          if ( skip ) continue;

          map<pair<int,int>,int> conn;
          for ( auto & x : D.ITo(l) )
          for ( auto & y : D.IFrom(r) )
               conn[make_pair(x,y)] = 0;
          for ( auto & id : dpaths_index[d] ) {
               const ReadPath & rp = dpaths[id];
               int pos = 0;
               for ( ; pos < int(rp.size()); pos++ )
                    if ( rp[pos] == d )
                         break;
               if ( pos == int(rp.size()) ) {
                    #pragma omp critical
                    {
                    //weird if this were to happen
                    PRINT(pos);
                    PRINT(id);
                    for ( auto & x : rp )
                         cout << x << " ";
                    cout << endl;
                    FatalErr( "dpaths_index and dpaths out of sync" );
                    }
               }
               if ( pos == 0 ) continue;
               if ( pos == int(rp.size()-1) ) continue;
               const pair<int,int> p ( rp[pos-1], rp[pos+1] );
               conn[p]++;
          }
          // repeat for rc edge
          for ( auto & id : dpaths_index[dinv[d]] ) {
               const ReadPath & rp = dpaths[id];
               int pos = 0;
               for ( ; pos < int(rp.size()); pos++ )
                    if ( rp[pos] == dinv[d] )
                         break;
               if ( pos == int(rp.size()) ) {
                    #pragma omp critical
                    {
                    //weird if this were to happen
                    PRINT(pos);
                    PRINT(id);
                    for ( auto & x : rp )
                         cout << x << " ";
                    cout << endl;
                    FatalErr( "dpaths_index and dpaths out of sync" );
                    }
               }
               if ( pos == 0 ) continue;
               if ( pos == int(rp.size()-1) ) continue;
               if ( rp[pos+1] >= dinv.isize() || rp[pos-1] >= dinv.isize() ) {
                    #pragma omp critical
                    {
                    cout << "rp = ";
                    for ( auto & x : rp )
                         cout << x << " ";
                    cout << endl;
                    PRINT3( D.E(), pos, dinv.size() );
                    cout << "ERROR!" << endl;
                    Scram(1);
                    }
               }
               const pair<int,int> p ( dinv[rp[pos+1]], dinv[rp[pos-1]] );
               conn[p]++;
          }
          
          vec<triple<int,int,int>> sef;
          for ( auto & x : conn )
               sef.push( x.second, x.first.first, x.first.second );
          ReverseSort( sef );

          int max = sef[1].first, min = sef[2].first;
          if ( min > 5 && max < 20 * min ) continue;
          if ( max < 10 ) continue;
          if ( max < 5*min ) continue;

          if (!IsUnique( sef[0].second, sef[0].third, sef[1].second, sef[1].third ))
               continue;
          #pragma omp critical
          {    pull.push( sef[0].second, sef[0].third,
                          sef[1].second, sef[1].third ); }
     }
     cout << Date( ) << ": making edits" << endl;
     int edits = 0;
     for ( auto & x : pull ) {
          //if ( toc[d] < 0 ) continue;
          
          /*
          cout << "CLONE: " << toc[d] << endl;
          cout << d << endl;
          cout << x.first << "(" << toc[x.first] << ") -> ";
          cout << x.second << "(" << toc[x.second] << ")" << endl;
          cout << x.third << "(" << toc[x.third] << ") -> ";
          cout << x.fourth << "(" << toc[x.fourth] << ")" << endl;
          */

          // Make edits
          const int d = D.IFrom(to_right[x.first], 0);
          const int le1 = x.first, re1 = x.second;
          const int le2 = x.third, re2 = x.fourth;
          const int l = to_left[d], r = to_right[d];
          const int il = to_left[dinv[d]], ir = to_right[dinv[d]];
          const int N = D.N();
          
          edits++;
          D.AddVertices(4);
          
          D.GiveEdgeNewToVxWithUpdate( le1, l, N, to_right );
          D.GiveEdgeNewFromVxWithUpdate( re1, r, N, to_left );
          D.GiveEdgeNewToVxWithUpdate( le2, l, N+1, to_right );
          D.GiveEdgeNewFromVxWithUpdate( re2, r, N+1, to_left );
          vec<int> &eo1 = D.OMutable ( le1 ), & eo2 = D.OMutable( le2 );
          eo1.append( D.O( d ) );
          eo2.append( D.O( d ) );
          D.DeleteEdgeFrom( l, 0 );
          to_left[d] = -1, to_right[d] = -1;

          D.GiveEdgeNewToVxWithUpdate( dinv[re1], il, N+2, to_right );
          D.GiveEdgeNewFromVxWithUpdate( dinv[le1], ir, N+2, to_left );
          D.GiveEdgeNewToVxWithUpdate( dinv[re2], il, N+3, to_right );
          D.GiveEdgeNewFromVxWithUpdate( dinv[le2], ir, N+3, to_left );
          vec<int> &reo1 = D.OMutable( dinv[le1] ), & reo2 = D.OMutable( dinv[le2] );
          vec<int> v1 = D.O(dinv[d]), v2 = D.O(dinv[d]);
          v1.append( D.O( dinv[le1] ) );
          v2.append( D.O( dinv[le2] ) );
          reo1 = v1;
          reo2 = v2;
          D.DeleteEdgeFrom( il, 0 );
          to_left[dinv[d]] = -1, to_right[dinv[d]] = -1;
/*
          // Update dpaths and dpaths_index
          for ( int pass = 1; pass <= 2; pass++ ) {
               const int x = (pass == 1 ? d : dinv[d] );
               for ( auto & id : dpaths_index[x] ) {
                    ReadPath & rp = dpaths[id];
                    const int n = rp.size();
                    int pos = 0;
                    for ( ; pos < n; pos++ )
                         if ( rp[pos] == x )
                              break;
                    for ( int j = pos; j < n-1; j++ )
                         rp[j] = rp[j+1];
                    rp.resize(n-1);
               }
               dpaths_index[x].clear();
          }
*/
     }
     cout << Date( ) << ": made " << edits << " edits" << endl;
}


// If there are sufficiently many read-pairs between two lines,
// and those two lines don't have a similar relationship with
// other lines, then connect them with a pair gap

void HookLinesWithReadPairs( const HyperBasevectorX & hb, 
     digraphE<vec<int>> & D, vec<int> & dinv, const vec<vec<vec<vec<int>>>> & dlines,
     const ReadPathVec & dpaths, VecULongVec & dpaths_index )
{
     vec<int> to_left, to_right;
     D.ToLeftParallel( to_left ), D.ToRightParallel( to_right );

     vec<int> linv;
     LineInv( dlines, dinv, linv );
     
     vec<int> dlens;
     ComputeDlens( hb, D, dlens );

     cout << Date( ) << ": indexing lines" << endl;
     vec<int> dtol( D.E(), -1 );
     //vec<int> ltoc( dlines.size(), -1 );
     for ( int l = 0; l < dlines.isize(); l++ ) {
          const auto & L = dlines[l];
          for ( int j = 0; j < L.isize( ); j++ )
               for ( auto & path : L[j] )
                    for ( auto & d : path ) {
                         dtol[d] = l;
                         //if ( ltoc[l] < 0 && toc[d] >= 0 )
                         //     ltoc[l] = toc[d];
                    }
     }
     
     map<pair<int,int>, int> lsupp;
     GetInsertsConnectingLines( dinv, linv, dpaths, dpaths_index, dtol, lsupp );

     vec<vec<int>> lconn ( dlines.size() );
     for ( auto & ls : lsupp ) {
          const int l1 = ls.first.first;
          const int l2 = ls.first.second;
          if ( ls.second < 5 ) continue;
          if ( l1 == l2  || l1 == linv[l2] ) continue;
          const int left = dlines[l1].back()[0].back();
          const int right = dlines[l2][0][0][0];
          if ( dlens[left] < 100 || dlens[right] < 100 ) continue;
          const Bool l1dead = D.From(to_right[left]).empty();
          const Bool l2dead = D.To(to_left[right]).empty();
          if ( l1dead && l2dead ) {
               set<int64_t> pids;
               int common = 0;
               for ( auto & id : dpaths_index[left] )
                    pids.insert( id/2 );
               for ( auto & id : dpaths_index[dinv[right]] ) {
                    if ( pids.count( id/2 ) )
                         common++;
               }
               if ( common < 5 ) continue;
               //cout << l1 << " <-> " << l2 << " = " << ls.second << endl;
               lconn[l1].push_back(l2);
          }
     }
     for ( auto & ls : lconn ) {
          if ( ls.nonempty() )
               Sort( ls );
     }
     vec<pair<int,int>> connect;
     vec<Bool> done( dlines.size(), False );
     for ( int l1 = 0; l1 < lconn.isize(); l1++ ) {
          if ( done[l1] ) continue;
          if (lconn[l1].size() != 1) continue;
          const int l2 = lconn[l1][0];
          if (lconn[linv[l2]].size() != 1) continue;
          if (lconn[linv[l2]][0] != linv[l1]) continue;
          const int left = dlines[l1].back()[0].back();
          const int right = dlines[l2][0][0][0];
          /*
          cout << "CONNECT: " << l1;
          if ( ltoc[l1] >= 0 )
               cout << "(CLONE:" << ltoc[l1] << ") ";
          cout << " --> " << l2;
          if ( ltoc[l2] >= 0 )
               cout << "(CLONE:" << ltoc[l2] << ")";
          cout << endl;
          */
          const int v1 = to_right[left], v2 = to_left[right];
          const int E = D.E();
          D.AddEdgeWithUpdate( v1, v2, {-1}, to_left, to_right );
          D.AddEdgeWithUpdate( to_right[dinv[right]], to_left[dinv[left]],
               {-1}, to_left, to_right );
          dinv.push_back( E+1 );
          dinv.push_back( E );
          done[l1] = True;
          done[linv[l2]] = True;
     }
     // update dpaths index
     dpaths_index.resize( D.E() );
}


// If there is an edge with zero support then delete it if there is an 
// alternate path that contains at least one edge with > 0 reads.
// This is a little dangerous, simplifies but potentially deletes true
// sequence

void DeleteZeroSupportEdgesWithAltPath( digraphE<vec<int>> & D, vec<int> & dinv,
     vec<Bool> & to_delete, const VecULongVec & dpaths_index )
{
     ForceAssertEq(to_delete.size(), D.E() );

     const int MAX_DEPTH = 10;
     
     vec<Bool> to_save ( D.E(), False );
     vec<int> expl;
     map<int,int> backtrack;
     int pass = 1;
     int tot = 0;
     while (True) {
     
     vec<int> to_left, to_right;
     D.ToLeftParallel( to_left );
     D.ToRightParallel( to_right );
     
     int dels = 0;
     cout << Date( ) << ": pass " << pass << endl;
     pass++;
     #pragma omp parallel for ordered firstprivate(expl, backtrack)
     for ( int d = 0; d < D.E(); d++ ) {
          if ( D.O(d)[0] < 0 ) continue;

          if ( dpaths_index[d].size() + dpaths_index[dinv[d]].size() > 0 )
               continue;

          const int v = to_left[d], target = to_right[d];
          
          // stores all edges explored EXCEPT d
          expl.clear();
          backtrack.clear();
          set<int> all = {v};
          set<int> vset = {v};
          vec<int> thepath;
          Bool alt = False;
          for ( int depth = 0; depth < MAX_DEPTH; depth++ ) {
               const int start = expl.size( );
               for ( auto & vp : vset ) {
                    for ( int i = 0; i != D.From(vp).isize(); i++ ) {
                         const int edge = D.IFrom( vp, i );
                         if ( IsPairGap( D.O(edge) ) ) continue;
                         if ( IsBarcodeOnlyGap( D.O(edge) ) ) continue;
                         if ( edge != d )
                              expl.push_back( edge );
                    }
               }
               vset.clear( );
               
               if ( expl.isize() == start ) break;

               for ( int i = start; i != expl.isize( ); i++ ) {
                    const int w = to_right[expl[i]];
                    if ( w == target ) {
                         thepath.push_back( expl[i] );
                         int vx = to_left[ expl[i] ];
                         while ( vx != v ) {
                              const int backedge = backtrack[vx];
                              vx = to_left[backedge];
                         }
                         for ( auto & x : thepath ) {
                              if ( dpaths_index[x].size() || 
                                        dpaths_index[dinv[x]].size() )
                                   alt = True;
                         }
                         break;
                    }
                    if ( all.count( w ) ) continue;
                    backtrack[w] = expl[i];
                    vset.insert( w );
                    all.insert( w );
               }
               if ( alt )
                    break;
          }

          if ( !alt ) continue;

          #pragma omp ordered
          {
               if ( !to_save[d] ) {
                    Bool mod = False;
                    for ( auto & d : thepath ) {
                         if ( to_delete[d] ) {
                              mod = True;
                              break;
                         }
                    }
                    if ( !mod ) {
                         to_delete[d] = True;
                         to_delete[dinv[d]] = True;
                         dels += 1 + int(d != dinv[d]);
                         for ( auto & d : thepath ) {
                              to_save[d] = True;
                              to_save[dinv[d]] = True;
                         }
                    }
               }
          }
     }
     tot += dels;
     if ( dels == 0 ) {
          break;
     }

     cout << Date( ) << ": identified " << dels << " redundant edges" << endl;
     D.DeleteEdgesParallel( to_delete );
     
     to_delete.clear();
     to_save.clear();
     to_delete.resize(D.E(), False );
     to_save.resize(D.E(), False );
     
     } // END OF MULTIPLE PASSES
     cout << Date( ) << ": deleted " << ToStringAddCommas(tot) << " edges" << endl;
}


//
//
//

void CleanTheAssembly( const HyperBasevectorX& hb, const vec<int>& inv,
     digraphE<vec<int>>& D, vec<int>& dinv, ReadPathVec& dpaths,
     const vec<Bool>& dup,
     const ReadPathVecX& pathsx, const vec<Bool>& pmask, const vec<int64_t>& bci,
     const vec< triple<int,int,int> >& qept, const VecIntVec& ebcx,
     const String& udir, const String& suffix,
     const Bool intermediates, const String INSTANCE, const int substart, 
     const int substop, const Bool LOWMEMB, const Bool LOOPS2_VERBOSE,
     const Bool ZIPPER, const Bool MORE )
{
     // Debugging control.

     auto WriteMe = [&]( int n )
     {    cout << Date( ) << ": ***** checksum for clean" << n << " = " 
               << D.CheckSum( ) << " *****" << endl;
          cout << Date( ) << ": mem = " << MemUsageGBString( )
               << ", peak = " << PeakMemUsageGBString( ) << endl;
          RepairShop( hb, inv, D, dinv );
          if (intermediates)
          {    cout << Date( ) << ": ***** WRITING clean" << n << " *****" << endl;
               {    String OUTDIR = udir + "/" + INSTANCE + "/clean" + ToString(n);
                    Mkdir777(OUTDIR);
                    BinaryWriter::writeFile( OUTDIR + "/a.sup", D );
                    BinaryWriter::writeFile( 
                         OUTDIR + "/a.sup.inv", dinv );    }    }
          Validate( hb, inv, D, dinv );    
          if ( n == substop )
          {    cout << endl << Date( ) << ": done, peak mem = " 
                    << PeakMemUsageGBString( ) << endl << endl;
               Scram(0);    }    };
     #define WM(s) WriteMe(s); start##s:

     // Heuristics.

     const int MIN_RATIO = 5;
     const int MAX_KILL = 4;
     const int MAX_KILLX = 2500;
     const double MIN_RATIOX = 20;
     const int LOOK_MERGE = 250;
     const int LOOK = 6; // raising to 8 may be very slow
     const int MESS_LONG_LINE = 2000;
     const int MAX_CAN_INS_DEL = 5;
     const int MIN_CAN_INS_RATIO = 4;

     // Declarations.

     vec<int> dels, dfw, lens, mult, linv, to_left, to_right;
     vec<vec<vec<vec<int>>>> dlines, dlines2;
     Bool del_verbose = True;
     vec<Bool> deleted;

     // Expand barcode index.

     cout << Date( ) << ": expanding barcode index" << endl;
     vec<int32_t> bc( bci.back( ), -1 );
     #pragma omp parallel for
     for ( int b = 0; b < bci.isize( ) - 1; b++ )
     {    int64_t start = bci[b], stop = bci[b+1];
          for ( int64_t j = start; j < stop; j++ )
               bc[j] = b;    }

     // Create a vector of integers, one for each read, such that
     // "having two" nonzero elements is enough.

     DataSet d;
     d.dt = ReadDataType::BAR_10X;
     d.start = 0;
     vec<DataSet> datasets = {d};
     vec<int64_t> bid( pathsx.size( ) );
     for ( int64_t id = 0; id < (int64_t) pathsx.size( ); id++ )
     {    int di;
          for ( di = 0; di < datasets.isize( ); di++ )
               if ( id < datasets[di].start ) break;
          const ReadDataType& dtype = datasets[di-1].dt;
          if ( dtype == ReadDataType::BAR_10X )
               bid[id] = (int64_t) pathsx.size( ) + bc[id] + 1;
          else if ( dtype == ReadDataType::UNBAR_10X ) bid[id] = 0;
          else if ( dtype == ReadDataType::PCR_FREE ) bid[id] = id + 1;    }

     // Jump.

     if ( substart == 1 ) goto start1;
     if ( substart == 2 ) goto start2;
     if ( substart == 3 ) goto start3;
     if ( substart == 4 ) goto start4;
     if ( substart == 5 ) goto start5;
     if ( substart == 6 ) goto start6;
     if ( substart == 7 ) goto start7;
     if ( substart == 8 ) goto start8;
     if ( substart == 9 ) goto start9;
     if ( substart == 10 ) goto start10;
     if ( substart == 11 ) goto start11;
     if ( substart == 18 ) goto start18;
     if ( substart == 20 ) goto start20;
     if ( substart == 21 ) goto start21;
     if ( substart == 22 ) goto start22;

     // Remove hanging ends.
     
     cout << Date( ) << ": removing hanging ends, peak mem = " 
          << PeakMemUsageGBString( ) << endl;
     lens.resize( D.E( ), 0 );
     #pragma omp parallel for
     for ( int e = 0; e < D.E( ); e++ )
     {    if ( D.O(e)[0] < 0 ) continue;
          for ( int j = 0; j < D.O(e).isize( ); j++ )
               lens[e] += hb.Kmers( D.O(e)[j] );    }
     DistancesToEndArr( D, lens, MAX_KILL * MIN_RATIO, True, dfw );
     #pragma omp parallel for
     for ( int v = 0; v < D.N( ); v++ )
     {    for ( int j1 = 0; j1 < D.From(v).isize( ); j1++ )
          for ( int j2 = 0; j2 < D.From(v).isize( ); j2++ )
          {    if ( j2 == j1 ) continue;
               int w1 = D.From(v)[j1], w2 = D.From(v)[j2];
               int d1 = D.IFrom(v,j1), d2 = D.IFrom(v,j2);
               int len1 = lens[d1] + dfw[w1], len2 = lens[d2] + dfw[w2];
               if ( len1 <= MAX_KILL && len2 >= MIN_RATIO * len1
                    && len2 >= MAX_KILL )
               {    
                    #pragma omp critical
                    {    dels.push_back( d1, dinv[d1] );    }    }    }    }

     // Clean unnecessary pair gaps from cells.

     FindLines( D, dinv, dlines, MAX_CELL_PATHS, MAX_CELL_DEPTH );
     DeletePairGapsInCells( D, dinv, dlines, dels, del_verbose );
     D.DeleteEdges(dels);     
     
     // Remove hanging ends.

     HangBeGone( hb, D, dinv );

     // Delete parallel gaps.

     DeleteParallelGaps( D, dinv, del_verbose );
     WM(1);

     // Delete some redundant gaps.

     DeleteSomeRedundantGaps( hb, D, dinv, del_verbose );
     
     // Capture simple loops.

     ReinsertLoops( hb, inv, D, dinv );
     RemoveSuperfluousLoops( hb, D, dinv );
     dels.clear( );
     CaptureLoops( hb, inv, D, dinv, False, False );
     WM(2);

     // Remove duplicate edges.

     deleted.resize( D.E( ), False );
     dels.clear( );
     for ( int v = 0; v < D.N( ); v++ )
     for ( int j1 = 0; j1 < D.From(v).isize( ); j1++ )
     {    int d1 = D.IFrom(v,j1);
          if ( deleted[d1] ) continue;
          for ( int j2 = j1 + 1; j2 < D.From(v).isize( ); j2++ )
          {    int d2 = D.IFrom(v,j2);
               if ( D.O(d2) != D.O(d1) ) continue;
               if ( D.From(v)[j1] != D.From(v)[j2] ) continue;
               deleted[d2] = deleted[ dinv[d2] ] = True;
               dels.push_back( d2, dinv[d2] );    }    }
     D.DeleteEdges(dels);
     RemoveUnneededVertices( D, dinv );
     CleanupCore( D, dinv );    
     WM(3);

     // Remove redundant cells.

     RemoveRedundantCells( D, dinv );
     WM(4);

     // Delete compound hanging ends.

     lens.resize_and_set( D.E( ), 0 );
     #pragma omp parallel for
     for ( int e = 0; e < D.E( ); e++ )
     {    if ( D.O(e)[0] < 0 ) continue;
          for ( int j = 0; j < D.O(e).isize( ); j++ )
               lens[e] += hb.Kmers( D.O(e)[j] );    }
     DistancesToEndArr( D, lens, MAX_KILLX * MIN_RATIOX, True, dfw );
     dels.clear( );
     FindCompoundHangs( D, dinv, lens, dfw, dels, MAX_KILLX, MIN_RATIOX, False );
     D.DeleteEdges(dels);
     RemoveUnneededVertices( D, dinv );
     CleanupCore( D, dinv );    
     WM(5);
     
     // Delete more redundant pair gaps.
     // Delete pair gaps at branches.  This is a big hammer.
     
     cout << Date( ) << ": start pair gap cleanup, peak mem = "
          << PeakMemUsageGBString( ) << endl;
     ZapPairGaps( D, dinv, del_verbose );
     {    
          // Weaker form of killing pair gaps at branches.

          PlaceReadsMasked( hb, D, dup, pathsx, pmask, dpaths );
          IntIndex dpaths_index( dpaths, D.E( ) );
          dels.clear( );
          const int MIN_SUPPORT = 2;
          for ( int v = 0; v < D.N( ); v++ )
          {    for ( int j1 = 0; j1 < D.From(v).isize( ); j1++ )
               {    int d1 = D.IFrom(v,j1);
                    if ( D.O(d1)[0] < 0 ) continue;
                    for ( int j2 = 0; j2 < D.From(v).isize( ); j2++ )
                    {    int d2 = D.IFrom(v,j2);
                         if ( IsPairGap( D.O(d2) ) )
                         {    int w = D.From(v)[j2];
                              vec<int> bcs;
                              for ( int k1 = 0; k1 < D.To(v).isize( ); k1++ )
                              for ( int k2 = 0; k2 < D.From(w).isize( ); k2++ )
                              {    int f = D.ITo(v,k1), g = D.IFrom(w,k2);
                                   if ( D.O(f)[0] < 0 || D.O(g)[0] < 0 ) continue;
                                   for ( int l = 0; l < dpaths_index.Count(f); l++ )
                                   {    int64_t id1 = dpaths_index.Val(f,l);
                                        int64_t id2 = (id1 % 2 == 0 ? id1+1 : id1-1);
                                        int b = bc[id1];
                                        if ( b == 0 ) continue;
                                        const ReadPath& p = dpaths[id2];
                                        for ( auto d : p )
                                        {    if ( dinv[d] == g )
                                             {    bcs.push_back(b);
                                                  break;    }    }    }    }
                              UniqueSort(bcs);
                              if ( bcs.isize( ) < MIN_SUPPORT )
                              {    dels.push_back( 
                                        d2, dinv[d2] );    }    }    }    }    }
          cout << Date( ) << ": killing " << dels.size( ) << " pair gaps"
               << " at branches" << endl;
          D.DeleteEdges(dels);
          RemoveUnneededVertices( D, dinv );
          CleanupCore( D, dinv );
          Validate( hb, inv, D, dinv );    }

     // Kill pair gap loops.

     {    dels.clear( );
          for ( int v = 0; v < D.N( ); v++ )
          for ( int j = 0; j < D.From(v).isize( ); j++ )
          {    int d = D.IFrom(v,j);
               if ( IsPairGap( D.O(d) ) && D.From(v)[j] == v ) 
                    dels.push_back(d);    }
          D.DeleteEdges(dels);
          RemoveUnneededVertices( D, dinv );
          CleanupCore( D, dinv );    }
     Validate( hb, inv, D, dinv );

     // Capture some loops.

     dels.clear( );
     CaptureMessyLoops( hb, inv, D, dinv, dels, True );
     D.DeleteEdges(dels);
     RemoveUnneededVertices( D, dinv );
     CleanupCore( D, dinv );    
     dels.clear( );
     CaptureCanonicalLoops( D, dinv, dels, False, False );
     D.DeleteEdges(dels);
     RemoveUnneededVertices( D, dinv );
     CleanupCore( D, dinv );    
     dels.clear( );
     CaptureSimpleLoops( D, dinv, dels, False, False );
     D.DeleteEdges(dels);
     RemoveUnneededVertices( D, dinv );
     CleanupCore( D, dinv );    
     WM(6);

     // Do more of same stuff done above.

     
     cout << Date( ) << ": merging short overlaps, in six passes" << endl;
     if ( LOWMEMB ) {
          for ( int pass = 1; pass <= 6; pass++ ) {
               MergeShortOverlapsLowMem(hb, inv, D, dinv, LOOK_MERGE, LOOK, False, False);
               if ( ZIPPER )
                    ZipperFast( D, dinv, True );
               else
                    Zipper( D, dinv, False, True );
          }
     } else {
          for ( int pass = 1; pass <= 6; pass++ )
               MergeShortOverlaps(hb, inv, D, dinv, LOOK_MERGE, LOOK, False, False);
     }

     CleanupCore( D, dinv );
     WM(7);
     DeleteSomeRedundantGaps( hb, D, dinv, del_verbose );
     HangBeGone( hb, D, dinv );
     DeleteParallelGaps( D, dinv, del_verbose );
     ZapPairGaps( D, dinv, del_verbose );
     WM(8);
     
     // Remove loops that are identical to another path.

     RemoveSuperfluousLoops( hb, D, dinv );

     // Remove very small components.

     cout << Date( ) << ": remove very small components, peak mem = "
          << PeakMemUsageGBString( ) << endl;
     RemoveVerySmallComponents( hb, D, dinv, 200 );
     WM(9);
     
     // Capture messy loops.

     CaptureMessyLoops2( hb, inv, D, dinv, MESS_LONG_LINE, 10, 20, LOOPS2_VERBOSE ); 
     WM(10);

     // Kill low unique.

     dels.clear( );
     KillLowUnique( hb, D, dels, False );
     D.DeleteEdges(dels);
     RemoveUnneededVertices( D, dinv );
     CleanupCore( D, dinv );

     // Delete stuff having weak support.

     cout << Date( ) << ": deleting stuff with weak support, peak mem = " 
          << PeakMemUsageGBString( ) << endl;
     PlaceReadsMasked( hb, D, dup, pathsx, pmask, dpaths );
     {
     IntIndex dpaths_index( dpaths, D.E( ) );
     ComputeMult( hb, D, mult );
     dels.clear( );
     for ( int v = 0; v < D.N( ); v++ )
     {    int n = D.From(v).size( );
          if ( n < 2 ) continue;
          vec<int> count(n), ids( n, vec<int>::IDENTITY );
          for ( int j = 0; j < n; j++ )
          {    int d = D.IFrom(v,j);
               count[j] 
                    = dpaths_index.Count(d) + dpaths_index.Count(dinv[d]);    }
          ReverseSortSync( count, ids );
          const int MAX_DEL = 2;
          const int MIN_RAT = 10;
          if ( count[0] >= 1 && count[1] <= MAX_DEL 
               && count[0] >= MIN_RAT * count[1] )
          {    for ( int j = 1; j < n; j++ )
               {    int d = D.IFrom( v, ids[j] );
                    Bool uniq = False;
                    if ( D.O(d)[0] >= 0 )
                         for ( auto e : D.O(d) ) if (mult[e] == 1) uniq = True;
                    if (uniq) dels.push_back( d, dinv[d] );    }    }    }
     }
     D.DeleteEdges(dels);
     RemoveUnneededVertices( D, dinv );
     CleanupCore( D, dinv );
     WM(11);

     // Peel off lines.

     // String OUTDIR = INDIR + ".res3";
     // Mkdir777(OUTDIR);
     // Ofstream( rout, OUTDIR + "/Peel.log" );
     {
     ostringstream rout;
     PlaceReadsMasked( hb, D, dup, pathsx, pmask, dpaths );
     {    VecULongVec dpaths_index;
          invert( dpaths, dpaths_index, D.E( ) );
          Peeler( hb, inv, D, dinv, dpaths, dpaths_index, rout );    }
     }
     WM(12);

     // Place reads again.

     cout << Date( ) << ": placing reads" << endl;
     PlaceReadsMasked( hb, D, dup, pathsx, pmask, dpaths );

     // Pull apart.

     cout << Date( ) << ": pulling apart" << endl;
     dels.clear( );
     PullApart( inv, D, dinv, dpaths, dels, False, 2 );
     D.DeleteEdges(dels);
     RemoveUnneededVertices( D, dinv );
     CleanupCore( D, dinv );
     Validate( hb, inv, D, dinv );

     // Capture messy loops.

     dels.clear( );
     CaptureMessyLoops( hb, inv, D, dinv, dels, True, MESS_LONG_LINE );
     D.DeleteEdges(dels);
     RemoveUnneededVertices( D, dinv );
     CleanupCore( D, dinv );

     // Kill inversion artifacts.  Note that it's possible that we could dispense
     // with this by relaxing the parameters in the same way in the local 
     // assemblies.

     PlaceReadsMasked( hb, D, dup, pathsx, pmask, dpaths );
     dels.clear( );
     {    IntIndex dpaths_index( dpaths, D.E( ) );
          Bool advanced = False;
          KillInversionArtifacts( hb, D, dinv, dpaths, dpaths_index, bid, dels,
               MAX_CAN_INS_DEL, MIN_CAN_INS_RATIO, advanced, False );    }
     D.DeleteEdges(dels);
     RemoveUnneededVertices( D, dinv );
     CleanupCore( D, dinv );

     // Fix a rare type of inversion artifact, where an inversion bubble has two
     // branches entering on each side.  We copy the bubble.

     for ( int v = 0; v < D.N( ); v++ )
     {    if ( D.From(v).size( ) != 2 || D.From(v)[0] != D.From(v)[1] ) continue;
          int w = D.From(v)[0], d1 = D.IFrom(v,0), d2 = D.IFrom(v,1);
          if ( D.O(d1)[0] < 0 || dinv[d1] != d2 ) continue;
          if ( D.To(v).size( ) != 2 || D.From(w).size( ) != 2 ) continue;
          int f1 = D.ITo(v,0), f2 = D.ITo(v,1), g1 = D.IFrom(w,0), g2 = D.IFrom(w,1);
          if ( g1 == dinv[f1] ) swap( g1, g2 );
          int N = D.N( );
          D.AddVertices(2);
          D.GiveEdgeNewToVx( f2, v, N );
          D.GiveEdgeNewFromVx( g2, w, N+1 );
          int x1 = D.AddEdge( N, N+1, D.O(d1) );
          int x2 = D.AddEdge( N, N+1, D.O(d2) );
          dinv.push_back( d2, d1 );
          dinv[d1] = x2, dinv[d2] = x1;    }

     // Capture messy loops.

     dels.clear( );
     CaptureMessyLoops( hb, inv, D, dinv, dels, True, MESS_LONG_LINE );
     D.DeleteEdges(dels);
     RemoveUnneededVertices( D, dinv );
     CleanupCore( D, dinv );    
     Validate( hb, inv, D, dinv );

     // Explode vertices at palindromes.

     /*
     cout << Date( ) << ": exploding palindromes" << endl;
     vec<int> to_left, to_right;
     D.ToLeft(to_left), D.ToRight(to_right);
     int np = 0;
     for ( int d = 0; d < D.E( ); d++ )
     {    if ( d == dinv[d] )
          {    if ( D.To( to_left[d] ).empty( ) ) continue;
               np++;
               int v = to_left[d];
               D.SplayVertexWithUpdate( v, to_left, to_right );
               int w = to_right[d];
               D.SplayVertexWithUpdate( w, to_left, to_right );    }    }
     cout << Date( ) << ": exploded " << np << " palindromes" << endl;
     */

     // Kill inversion artifacts.  Note that it's possible that we could dispense
     // with this by relaxing the parameters in the same way in the local 
     // assemblies.  Also note that we do nearly the same thing above, perhaps
     // unnecessarily.

     PlaceReadsMasked( hb, D, dup, pathsx, pmask, dpaths );
     dels.clear( );
     {    IntIndex dpaths_index( dpaths, D.E( ) );
          Bool advanced = True;
          KillInversionArtifacts( hb, D, dinv, dpaths, dpaths_index, bid, dels,
               MAX_CAN_INS_DEL, MIN_CAN_INS_RATIO, advanced, False );    }
     D.DeleteEdges(dels);
     RemoveUnneededVertices( D, dinv );
     CleanupCore( D, dinv );
     WM(13);

     // Pull apart.

     BigPullInv( hb, inv, D, dinv, 10000 );

     // Delete compound hanging ends.

     lens.resize_and_set( D.E( ), 0 );
     #pragma omp parallel for
     for ( int e = 0; e < D.E( ); e++ )
     {    if ( D.O(e)[0] < 0 ) continue;
          for ( int j = 0; j < D.O(e).isize( ); j++ )
               lens[e] += hb.Kmers( D.O(e)[j] );    }
     DistancesToEndArr( D, lens, MAX_KILLX * MIN_RATIOX, True, dfw );
     dels.clear( );
     FindCompoundHangs( D, dinv, lens, dfw, dels, MAX_KILLX, MIN_RATIOX, False );
     D.DeleteEdges(dels);
     RemoveUnneededVertices( D, dinv );
     CleanupCore( D, dinv );    

     // Remove very small components.

     cout << Date( ) << ": remove very small components, peak mem = "
          << PeakMemUsageGBString( ) << endl;
     RemoveVerySmallComponents( hb, D, dinv, 300 );
     WM(14);

     // Join using barcodes.

     BarcodeJoin( hb, inv, D, dinv, dpaths, dup, pathsx, pmask, bci, qept, ebcx );
     WM(15);

     // Snip flip squares.

     SnipFlipSquares( D, dinv, dpaths );
     Validate( hb, inv, D, dinv );
     WM(16);

     // Capture messy loops.

     CaptureMessyLoops2( hb, inv, D, dinv, 25000, 10, 20, LOOPS2_VERBOSE );

     // Pull apart.

     BigPullInv( hb, inv, D, dinv, 5000 );

     // Join using barcodes.

     BarcodeJoin( hb, inv, D, dinv, dpaths, dup, pathsx, pmask, bci, qept, ebcx );

     // Delete gap edges in certain special and very rare configurations.

     FindLinesFixed( D, dinv, dlines, 
          MAX_CELL_PATHS_EVALUATION, MAX_CELL_DEPTH_EVALUATION );
     FindLineLines( D, dinv, dlines, dlines2 );
     LineInv( dlines, dinv, linv );
     D.ToLeft(to_left), D.ToRight(to_right);
     {
     vec<int> verts( 2 * dlines.size( ) );
     #pragma omp parallel for
     for ( int i = 0; i < dlines.isize( ); i++ )
     {    int v = to_left[ dlines[i].front( )[0][0] ];
          int w = to_right[ dlines[i].back( )[0][0] ];
          verts[2*i] = v, verts[2*i+1] = w;    }
     UniqueSort(verts);
     int N = verts.size( );
     vec<vec<int>> from(N), to(N), from_edge_obj(N), to_edge_obj(N);
     for ( int i = 0; i < dlines.isize( ); i++ )
     {    int v = BinPosition( verts, to_left[ dlines[i].front( )[0][0] ] );
          int w = BinPosition( verts, to_right[ dlines[i].back( )[0][0] ] );
          from[v].push_back(w), to[w].push_back(v);
          from_edge_obj[v].push_back(i), to_edge_obj[w].push_back(i);    }
     #pragma omp parallel for
     for ( int v = 0; v < N; v++ )
     {    SortSync( from[v], from_edge_obj[v] );
          SortSync( to[v], to_edge_obj[v] );    }
     vec<int> edges( dlines.size( ), vec<int>::IDENTITY );
     digraphE<int> D2( from, to, edges, to_edge_obj, from_edge_obj );
     vec<int> to_left2, to_right2;
     D2.ToLeft(to_left2), D2.ToRight(to_right2);
     dels.clear( );
     for ( int v = 0; v < D2.N( ); v++ )
     {    if ( D2.From(v).size( ) != 2 ) continue;
          int l1 = D2.IFrom(v,0), l2 = D2.IFrom(v,1);
          if ( !dlines[l2].solo( ) || !IsPairGap( D.O( dlines[l2][0][0][0] ) ) )
               swap( l1, l2 );
          if ( !dlines[l2].solo( ) || !IsPairGap( D.O( dlines[l2][0][0][0] ) ) )
               continue;
          int d2 = dlines[l2][0][0][0];
          int w = to_right2[l1], r = to_right2[l2];
          if ( !D2.From(r).solo( ) ) continue;
          int l3 = D2.IFrom(r,0);
          if ( linv[l1] != l3 ) continue;
          int s = to_right2[l3];
          Bool gap2 = False;
          for ( int j = 0; j < D2.From(w).isize( ); j++ )
               if ( D2.From(w)[j] == s && D2.OFrom(w,j) == linv[l2] ) gap2 = True;
          if (gap2) dels.push_back( d2, dinv[d2] );    }
     D.DeleteEdges(dels);
     RemoveUnneededVertices( D, dinv );
     CleanupCore( D, dinv );
     }

     // Remove very small components.

     cout << Date( ) << ": remove very small components, peak mem = "
          << PeakMemUsageGBString( ) << endl;
     RemoveVerySmallComponents( hb, D, dinv, 300 );

     // Remove duff at ends of long lines including cells.

     RemoveDuff2( hb, D, dinv );

     // Remove hanging ends.

     HangBeGone( hb, D, dinv );

     // Make barcode-only joins.

     for ( int pass = 1; pass <= 3; pass++ )
     {    BarcodeJoin( hb, inv, D, dinv, dpaths, dup, 
               pathsx, pmask, bci, qept, ebcx );    }
     WM(17);

     // Fix misassemblies.

     cout << Date( ) << ": fixing misassemblies" << endl;
     PlaceReadsMasked( hb, D, dup, pathsx, pmask, dpaths );
     FixMisassemblies( hb, dup, bci, pathsx, D, dinv, dpaths );
     CleanupCore( D, dinv );

     // Explode unidirectional vertices.

     for ( int v = 0; v < D.N( ); v++ )
     {    if ( ( D.From(v).empty( ) && D.To(v).size( ) > 1 )
               || ( D.From(v).size( ) > 1 && D.To(v).empty( ) ) )
          {    D.SplayVertex(v);    }    }
     Validate( hb, inv, D, dinv );

     // Do stuff like done above.

     RemoveDuff2( hb, D, dinv, 2500, 100 );
     Validate( hb, inv, D, dinv );
     CaptureMessyLoops2( hb, inv, D, dinv, 25000, 10, 20, LOOPS2_VERBOSE );
     Validate( hb, inv, D, dinv );
     for ( int pass = 1; pass <= 3; pass++ )
          BarcodeJoin(hb, inv, D, dinv, dpaths, dup, pathsx, pmask, bci, qept, ebcx);
     Validate( hb, inv, D, dinv );
     CaptureMessyLoops2( hb, inv, D, dinv, 25000, 10, 20, LOOPS2_VERBOSE );
     HangBeGone( hb, D, dinv, 10, 1000 );
     RemoveDuff2( hb, D, dinv, 3000, 100 );
     CaptureMessyLoops2( hb, inv, D, dinv, 25000, 25, 100, LOOPS2_VERBOSE );
     BarcodeJoin( hb, inv, D, dinv, dpaths, dup, pathsx, pmask, bci, qept, ebcx );

     // Appears to kill 0.5% of coverage, but raises N50:
     
     // RemoveDuffx( hb, D, dinv, 3000, 100, True );

     // A bunch of interesting stuff, commented out:

     /*

     // Kill low unique.

     vec<int> dels;
     KillLowUnique( hb, D, dels, False );
     // KillLowUniqueFrac( hb, D, dels );
     D.DeleteEdges(dels);
     RemoveUnneededVertices( D, dinv );
     CleanupCore( D, dinv );

     // Describe lines that got deleted (TEMP).
     
     vec<int> ls;
     for ( auto d : dels )
     {    int l = tol[d];
          if ( l >= 0 ) ls.push_back(l);    }
     UniqueSort(ls);
     cout << "\ndeleted lines: " << printSeq(ls) << endl << endl;

     // Wherever there is a branch, with support >=10 to 0, liberate the weak
     // branch, so long as it contains some unique kmers.

     PlaceReadsMasked( hb, D, dup, pathsx, pmask, dpaths );
     IntIndex dpaths_index( dpaths, D.E( ) );
     vec<int> to_left, to_right, mult;
     D.ToLeft(to_left), D.ToRight(to_right);
     ComputeMult( hb, D, mult );
     for ( int v = 0; v < D.N( ); v++ )
     {    if ( !D.To(v).solo( ) || !D.From(v).duo( ) ) continue;
          int g = D.ITo(v,0);
          if ( D.O(g)[0] < 0 ) continue;
          for ( int z = 0; z < 2; z++ )
          {    int d1 = D.IFrom(v,z), d2 = D.IFrom(v,1-z);
               if ( D.O(d1)[0] < 0 || D.O(d2)[0] < 0 ) continue;
               int rd1 = dinv[d1], rd2 = dinv[d2];
               Bool solo1 = False;
               for ( auto e : D.O(d1) ) if ( mult[e] == 1 ) solo1 = True;
               if ( !solo1 ) continue;
               int n1 = 0, n2 = 0;
               for ( int mpass = 1; mpass <= 2; mpass++ )
               for ( int pass = 1; pass <= 2; pass++ )
               {    int f;
                    if ( mpass == 1 ) f = ( pass == 1 ? d1 : dinv[d1] );
                    else f = ( pass == 1 ? d2 : dinv[d2] );
                    for ( int j = 0; j < dpaths_index.Count(f); j++ )
                    {    int64_t id = dpaths_index.Val( f, j );
                         const ReadPath& p = dpaths[id];
                         for ( int l = 0; l < (int) p.size( ) - 1; l++ )
                         {    if ( pass == 1 )
                              {    if ( p[l] == g && p[l+1] == f )
                                   {    ( mpass == 1 ? n1 : n2 )++;
                                        break;    }    }
                              else
                              {    if ( p[l] == f && p[l+1] == dinv[g] )
                                   {   ( mpass == 1 ? n1 : n2 )++;
                                        break;    }    }    }    }    }
               if ( n1 > 0 ) continue;
               if ( n2 < 10 || n2 < 10 * n1 ) continue;
               int N = D.N( );
               D.AddVertices(2);
               D.GiveEdgeNewFromVx( d1, v, N );
               D.GiveEdgeNewToVx( rd1, to_right[rd1], N+1 );
               Bool verbose = True;
               if (verbose)
               {    cout << "liberating " << d1 << " on left" << endl;
                    cout << "liberating " << rd1 << " on right" << endl;    }
               break;    }    }
     RemoveUnneededVertices( D, dinv );
     CleanupCore( D, dinv );
     HangBeGone( hb, D, dinv );

     */

     WM(18);
     
     // Number of passes to use while placing reads, each extra pass
     // after the first takes about 3 mins.

     if ( MORE ){
          const int MAX_PASSES=1;
          PlaceReadsMasked( hb, D, dup, pathsx, pmask, dpaths );
          
          VecULongVec dpaths_index;
          cout << Date( ) << ": inverting dpaths" << endl;
          invert ( dpaths, dpaths_index, D.E( ) );

          // STAGE 19
          // kill off zero weight components
          // Next kill pair gaps that are not informative
          // BC gaps that have identical sequence on both ends.
          // Zipper to kill some bubbles with identical edges
          // Remove some hanging ends
          
          cout << Date( ) << ": deleting components with zero support" << endl;
          KillZeroWeightComponents( hb, D, dinv, dpaths_index, True );
          
          to_left.clear();
          to_right.clear();
          D.ToLeftParallel( to_left );
          D.ToRightParallel( to_right );
          PairGapKiller( D, dinv, 10 );
          
          DeleteRedundantBCGaps( hb, D, dinv );
          
          ZipperFast( D, dinv, False );

          cout << Date( ) << ": deleting some hanging ends" << endl;
          {
               const int MIN_RATIO = 15;
               const int MAX_KILL  = 4000;
               const int verbosity = 1;
               const Bool NEW = True;
               HangBeGone( hb, D, dinv, MIN_RATIO, MAX_KILL, verbosity, NEW );
          }
     }

     WM(19);

     // STAGE 20:
     
     // Clean up based on read placement
     // 1. some inversion artifacts
     // 2. hanging lines
     // 3. zero support edges that can be replaced with a path with >= 1 support.
     if ( MORE ) {
          VecULongVec dpaths_index;
          
          // Update data structures
          cout << Date( ) << ": finding lines" << endl;
          FindLines( D, dinv, dlines, MAX_CELL_PATHS, MAX_CELL_DEPTH );
          
          const int MAX_PASSES=1;
          PlaceLinkedReadsMasked( hb, inv, D, dinv, dlines, dup, bci, pathsx, dpaths,
               pmask, MAX_PASSES, True );
          
          cout << Date( ) << ": inverting dpaths" << endl;
          invert ( dpaths, dpaths_index, D.E( ) );

          cout << Date( ) << ": finding lines" << endl;
          FindLines( D, dinv, dlines, 2, MAX_CELL_DEPTH );

          vec<Bool> to_delete( D.E(), False );
          MarkInversionArtifacts( hb, inv, D, dinv, dlines, dpaths_index, 
                                   to_delete, False );
          
          MarkDuplicateHangingLines( hb, D, dinv, dlines, dpaths_index, to_delete );
          
          int tot = 0;
          for ( auto & isdel : to_delete )
               if ( isdel ) tot++;
          cout << Date( ) << ": marked " << ToStringAddCommas(tot)
               << " edges for deletion" << endl;
 
          DeleteZeroSupportEdgesWithAltPath( D, dinv, to_delete, dpaths_index );
          
          cout << Date( ) << ": cleaning up" << endl;
          RemoveUnneededVertices( D, dinv );
          CleanupCore( D, dinv );
          RemoveUnneededVertices( D, dinv );
          CleanupCore( D, dinv );

     }

     // Update data structures
     
     WM(20);
     
     // STAGE 21: 
     // connect lines with pair gaps when possible
     // liberate branches with zero support 
     // delete weak bridges
     // kill disconnected components with zero weight
     // delete hanging ends
     if ( MORE ) {
          VecULongVec dpaths_index;

          cout << Date( ) << ": finding lines" << endl;
          FindLines( D, dinv, dlines, MAX_CELL_PATHS, MAX_CELL_DEPTH );

          PlaceLinkedReads( hb, inv, D, dinv, dlines, dup, bci, pathsx, dpaths,
               3, True );
          
          cout << Date( ) << ": inverting dpaths" << endl;
          invert ( dpaths, dpaths_index, D.E( ) );

          HookLinesWithReadPairs( hb, D, dinv, dlines, dpaths, dpaths_index );

          vec<Bool> to_delete( D.E(), False );
          cout << Date( ) << ": liberating zero support branches" << endl;
          KillZeroSupportBranches( hb, inv, D, dinv, dpaths, dpaths_index,
               to_delete, False );

          to_left.clear();
          to_right.clear();
          D.ToLeftParallel( to_left ), D.ToRightParallel( to_right );
          MarkWeakBridges( D, dinv, to_left, to_right, dpaths_index, to_delete );

          D.DeleteEdgesParallel( to_delete );

          cout << Date( ) << ": deleting some hanging ends" << endl;
          const int MIN_RATIO = 20;
          const int MAX_KILL  = 1000;
          const int verbosity = 1;
          const Bool NEW = True;
          HangBeGone( hb, D, dinv, MIN_RATIO, MAX_KILL, verbosity, NEW );
     }
     

     WM(21);

     // STAGE 22:
     // Read-pair pull apart
     if ( MORE ) {
          VecULongVec dpaths_index;
          
          // Update data structures
          
          cout << Date( ) << ": finding lines" << endl;
          FindLines( D, dinv, dlines, MAX_CELL_PATHS, MAX_CELL_DEPTH );
          
          cout << Date( ) << ": USING PLACELINKEDREADS" << endl;
          const int MAX_PASSES=1;
          PlaceLinkedReadsMasked( hb, inv, D, dinv, dlines, dup, bci,
               pathsx, dpaths, pmask, MAX_PASSES, True );
     
          cout << Date( ) << ": inverting dpaths" << endl;
          invert ( dpaths, dpaths_index, D.E( ) );

          ReadPairPullApart( D, dinv, dpaths, dpaths_index );
          RemoveUnneededVertices( D, dinv );
          CleanupCore( D, dinv );
     }

     WM(22);
     
     // STAGE 23:
     // Inversion bubble cleanup and tangle cleanup
     // (initially very targeted w.r.t topologies)
     //
     if ( MORE ) {
          cout << Date( ) << ": finding lines" << endl;
          FindLines( D, dinv, dlines, MAX_CELL_PATHS, MAX_CELL_DEPTH );

          PlaceReadsMasked( hb, D, dup, pathsx, pmask, dpaths );
          VecULongVec dpaths_index;
          invert ( dpaths, dpaths_index, D.E( ) );

          LineGraphOps lgo( hb, inv, pathsx, dup, bci, qept, ebcx, 
                              D, dinv, dlines, dpaths, dpaths_index);
          auto pivots = lgo.FindInvPivots();
          lgo.ScoreInvAndResolve(pivots);
     }
     
     WM(23);
}
