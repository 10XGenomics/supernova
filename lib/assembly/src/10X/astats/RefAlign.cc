// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

// RefAlign.  Align a supernova assembly to a reference sequence.  Intended for 
// aligning a human assembly to hg19.
//
// Probably robust to unused edges in D.

// Mr.X:
// SUMMARY
// kmers        count  %aligned  %multiply
// 1 - 249  4,060,378     74.95       5.29
// 250-499    196,234     75.42      18.44
// 500-999    164,383     84.63       9.01
// 1000 -   1,109,403     97.60       1.18
// weakly compatible abutting solo placements: 93.64%
// kmers aligned uniquely and consistently: 96%
// errs per kb = 3.36
// time used = 6.44 minutes, peak mem = 356.81 GB

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"
#include "PrintAlignment.h"
#include "math/HoInterval.h"
#include "pairwise_aligners/AlignFromMutmersAndSW.h"
#include "pairwise_aligners/Mutmer.h"
#include "pairwise_aligners/SmithWatFree.h"
#include "pairwise_aligners/SmithWatAffine.h"
#include "paths/HyperBasevector.h"
#include "paths/long/MakeKmerStuff.h"
#include "random/Random.h"
#include "10X/astats/RefAlign.h"

template<int K> void RefAlignCore(

     // Edge to be aligned.

     const vec<int>& X,

     // Assembly info.

     const int HBK, 
     const vecbasevector& tigs,
     const vec<int>& inv,
     const digraphE<vec<int>>& D,
     const vec<int>& dinv,
     const vec<vec<vec<vec<int>>>>& dlines,
     const vecbasevector& genome,
     const MasterVec< SerfVec<triple<int,int,int> > >& alignsb,

     // Output and output control.

     SerfVec<refalign>& galigns,
     const int verbosity,
     String& report,

     // Scratch data carried around to reduce memory allocation overhead.

     RefAlignScratch<K>& s,

     // Use only this g.

     const int g_use,

     // More control.

     const Bool allow_second_best,
     const Bool optb
)
{
     // Heuristics.

     const int MAX_ABS_OFFSET = 1000;
     const int MAX_TO_END = 300;
     const int MAX_ERRORS = 2000;
     const double END_STRETCH = 3.0;
     const int MIN_DELTA = -100;
     const int MAX_DELTA = 2500;
     const int MAX_ADD = 2000;
     const int MAX_LOCS_DIFF = 100;
     const int BW_ADD = 10;

     // Sanity check input.

     if ( X.nonempty( ) && X[0] < 0 )
     {    
          #pragma omp critical
          {    cout << "\nInternal error.  RefAlignCore appears to have been passed "
                    << "a gap edge.\n" << endl;
               TracebackThisProcess( );
               Scram(0);    }    }
     for ( auto e : X )
     {    if ( e >= inv.isize( ) )
          {    
               #pragma omp critical
               {    cout << "\nInternal error.  RefAlignCore appears to have been "
                         << "passed a non-current inv.\n" << endl;
                    TracebackThisProcess( );
                    Scram(0);    }    }
          if ( e < 0 )
          {    
               #pragma omp critical
               {    cout << "\nInternal error.  RefAlignCore appears to have been "
                         << "passed a scrambled X, found the value " << e << ".\n" 
                         << endl;
                    TracebackThisProcess( );
                    Scram(0);    }    }    }

     // Define some functions.

     auto Cat = [&]( const vec<int>& e )
     {    basevector x = tigs[ e[0] ];
          for ( int j = 1; j < e.isize( ); j++ )
               x = TrimCat( HBK, x, tigs[ e[j] ] );
          return x;    };
     auto Squarex = [&]( const int64_t t ){ return t*t; };

     // Localize scratch variables.

     vec<vec<triple<int,int,int>>>& locs = s.locs;
     vec< triple< pair<int,int>, triple< int, Bool, int >, String > >& 
          valigns = s.valigns;
     vec< pair<ho_interval,ho_interval> >& M = s.M;
     vec< triple<kmer<K>,int,int> >& kmers_plus = s.kmers_plus;
     basevector& B = s.B;
     basevector& G = s.G;
     vecbasevector& blob = s.blob;
     vec<Bool>& to_delete = s.to_delete;
     vec<int>& path = s.path;

     // Start reporting.

     ostringstream out;

     // Proceed.

     Bool punted = False;
     for ( int j = 0; j < 2; j++ ) locs[j].clear( );
     valigns.clear( ), galigns.clear( );
     vec<int> RX(X);
     RX.ReverseMe( );
     for ( int i = 0; i < RX.isize( ); i++ )
          RX[i] = inv[ RX[i] ];
     for ( int ps = 0; ps < 2; ps++ )
     {    const vec<int>& Y = ( ps == 0 ? X : RX );
          int dlensf = 0;
          for ( int i = 0; i < Y.isize( ); i++ )
               dlensf += tigs[ Y[i] ].isize( ) - HBK + 1;
          int offset = 0;
          for ( int i = 0; i < Y.isize( ); i++ )
          {    int e = Y[i];
               if ( e < (int) alignsb.size( ) )
               {    for ( int j = 0; j < (int) alignsb[e].size( ); j++ )
                    {    int g = alignsb[e][j].first;
                         if ( g_use >= 0 && g != g_use ) continue;
                         int start = alignsb[e][j].second + offset;
                         int stop = alignsb[e][j].third 
                              + offset - tigs[e].isize( ) + dlensf + HBK - 1;
                         locs[ps].push( g, start, stop );    }    }
               offset -= tigs[e].isize( ) - HBK + 1;    }
          UniqueSort(locs[ps]);

          // Merge overlapping intervals.

          vec<Bool> to_delete( locs[ps].size( ), False );
          for ( int i = 0; i < locs[ps].isize( ); i++ )
          {    int j;
               for ( j = i + 1; j < locs[ps].isize( ); j++ )
               {    if ( locs[ps][j].first != locs[ps][i].first ) break;
                    if ( locs[ps][j].second > locs[ps][j-1].third ) break;
                    to_delete[j-1] = True;
                    locs[ps][j].second = locs[ps][j-1].second;
                    locs[ps][j].third 
                         = Max( locs[ps][j-1].third, locs[ps][j].third );    }
               i = j -1;    }
          EraseIf( locs[ps], to_delete );

          // Process locations.

          for ( int i = 0; i < locs[ps].isize( ); i++ )
          {    ostringstream xout;
               if ( verbosity >= 1 )
               {    xout << ( ps == 0 ? "fw" : "rc" ) << " locs:\n";
                    xout << "[" << i+1 << "] " << locs[ps][i].first << ":"
                         << locs[ps][i].second << "-" << locs[ps][i].third 
                         << "\n";    }
               int j, g = locs[ps][i].first; 
               int start = locs[ps][i].second, stop = locs[ps][i].third;
               for ( j = i + 1; j < locs[ps].isize( ); j++ )
               {    if ( locs[ps][j].first != locs[ps][i].first ) break;
                    if ( locs[ps][j].second - locs[ps][j-1].second 
                         > MAX_LOCS_DIFF ) 
                    {    break;    }
                    stop = Max( stop, locs[ps][j].third );    }
               start = Max( 0, start ), stop = Min( stop, genome[g].isize( ) );
               if ( verbosity >= 2 ) PRINT2_TO( xout, start, stop );
               if ( start >= stop ) 
               {    if ( verbosity >= 2 ) out << xout.str( );
                    i = j - 1;
                    continue;    }
               B = Cat(Y);
               G.SetToSubOf( genome[g], start, stop - start );
               if ( verbosity >= 1 )
               {    xout << "B has size " << B.size( ) << ", G has size "
                         << G.size( ) << endl;    }

               // Two passes to allow enlargement of G as needed.

               for ( int pass = 1; pass <= 2; pass++ )
               {    M.clear( );

                    // Find matches.

                    {    blob[0] = B, blob[1] = G;
                         MakeKmerLookup3( blob, kmers_plus );
                         for ( int i = 0; i < kmers_plus.isize( ); i++ )
                         {    int j;
                              for ( j = i + 1; j < kmers_plus.isize( ); j++ )
                              {    if ( kmers_plus[j].first != kmers_plus[i].first )
                                        break;    }
                              for ( int k1 = i; k1 < j; k1++ )
                              for ( int k2 = i; k2 < j; k2++ )
                              {    if ( kmers_plus[k1].second != 0 ) continue;
                                   if ( kmers_plus[k2].second != 1 ) continue;
                                   int pos1 = kmers_plus[k1].third;
                                   int pos2 = kmers_plus[k2].third;
                                   if ( pos1 > 0 && pos2 > 0 
                                        && B[pos1-1] == G[pos2-1] )
                                   {    continue;    }
                                   int p1 = pos1+K, p2 = pos2+K;
                                   for ( ; p1 < B.isize( ); p1++, p2++ )
                                   {    if ( p2 == G.isize( ) ) break;
                                        if ( B[p1] != G[p2] ) break;    }
                                   M.push( ho_interval(pos1, p1), 
                                        ho_interval(pos2, p2) );    }
                              i = j - 1;    }    }

               /*
               if ( verbosity >= 1 )
               {    xout << "initial matches:\n";
                    for ( int i = 0; i < M.isize( ); i++ )
                    {    xout << "[" << i+1 << "] match A[" 
                              << M[i].first << ") to " << "G[" << M[i].second
                              << "), l = " << M[i].first.Length( ) 
                              << endl;    }    }
               */

                    // Delete subsumed matches.  Only checking on assembly side,
                    // could do genome side too.
     
                    Sort(M);
                    to_delete.resize_and_set( M.size( ), False );
                    for ( int a = 0; a < M.isize( ); a++ )
                    for ( int b = 0; b < M.isize( ); b++ )
                    {    if (  !to_delete[b]
                              && ProperSubset( M[a].first, M[b].first ) )
                         {    to_delete[a] = True;    }    }
                    EraseIf( M, to_delete );

                    // Enlarge G.

                    if ( pass == 1 )
                    {    int u1 = 1000000000, u2 = 1000000000;
                         for ( int i = 0; i < M.isize( ); i++ )
                         {    u1 = Min( u1, M[i].second.Start( ) 
                                   - M[i].first.Start( ) );
                              u2 = Min( u2, ( G.isize( ) - M[i].second.Stop( ) )
                                   - ( B.isize( ) - M[i].first.Stop( ) ) );    }
                         if ( u1 >= -MAX_ADD && u2 >= -MAX_ADD )
                         {    if ( u1 < 0 ) start += u1;
                              if ( u2 < 0 ) stop -= u2;
                              start = Max( 0, start );
                              stop = Min( genome[g].isize( ), stop );
                              G.SetToSubOf( genome[g], start, stop - start );
                              if ( verbosity >= 1 )
                              {    xout << "now B has size " << B.size( ) 
                                        << ", G has size " << G.size( ) 
                                        << endl;    }    }
                         i = j - 1;
                         if ( verbosity >= 2 ) out << xout.str( );
                         continue;    }    }

               // Trim overlapping matches.

               for ( int i = 0; i < M.isize( ) - 1; i++ )
               {    int over = 0;
                    if ( M[i].first.Stop( ) > M[i+1].first.Start( ) )
                         over = M[i].first.Stop( ) - M[i+1].first.Start( );
                    if ( M[i].second.Stop( ) > M[i+1].second.Start( ) )
                    {    over = Max( over, M[i].second.Stop( ) 
                              - M[i+1].second.Start( ) );    }
                    if ( over < M[i].first.Length( ) 
                              && over < M[i].second.Length( )
                              && over < M[i+1].first.Length( ) 
                              && over < M[i+1].second.Length( ) )
                    {    M[i].first.AddToStop(-over);
                         M[i].second.AddToStop(-over);    }    }

               // Delete subsize match blocks.

               to_delete.resize_and_set( M.size( ), False );
               for ( int i = 0; i < M.isize( ); i++ )
                    if ( M[i].first.Length( ) < K ) to_delete[i] = True;
               EraseIf( M, to_delete );

               // Merge across SNPs.  The number 8 comes from mismatch
               // penalty = 3, gap opening penalty = 12: (12+12)/3 = 8.
               //
               // BUG.
               // Note that nerrors is computed incorrectly, messed up when
               // we do shortest path.

               int nerrors = 0;
               to_delete.resize_and_set( M.size( ), False );
               for ( int i = 0; i < M.isize( ) - 1; i++ )
               {    int gap1 = M[i+1].first.Start( ) - M[i].first.Stop( );
                    int gap2 = M[i+1].second.Start( ) - M[i].second.Stop( );
                    if ( gap1 > 0 && gap1 == gap2 )
                    {    int mis = 0;
                         const int max_mis = 8;
                         for ( int j = 0; j < gap1; j++ )
                         {    if ( B[ M[i].first.Stop( ) + j ]
                                   != G[ M[i].second.Stop( ) + j ] )
                              {    mis++;
                                   if ( mis > max_mis ) break;    }    }
                         if ( mis <= max_mis )
                         {    to_delete[i] = True;
                              nerrors += mis;
                              M[i+1].first.SetStart( M[i].first.Start( ) );
                              M[i+1].second.SetStart( 
                                   M[i].second.Start( ) );    }    }    }
               EraseIf( M, to_delete );

               /*
               if ( verbosity >= 1 )
               {    xout << "post merging matches:\n";
                    for ( int i = 0; i < M.isize( ); i++ )
                    {    xout << "[" << i+1 << "] match A[" 
                              << M[i].first << ") to " << "G[" << M[i].second
                              << "), l = " << M[i].first.Length( ) << endl;    }    }
               */

               // Find a shortest path through the match blocks.  Let each match
               // block be a vertex in a digraph, plus add begin and end vertices.

               int N = M.size( ) + 2;
               vec<vec<int>> from(N), to(N), from_edge_obj(N), to_edge_obj(N);
               vec<double> edges;
               for ( int i1 = 0; i1 < M.isize( ); i1++ )
               for ( int i2 = 0; i2 < M.isize( ); i2++ )
               {    if ( i1 == i2 ) continue;
                    from[i1].push_back(i2), to[i2].push_back(i1);
                    from_edge_obj[i1].push_back( edges.size( ) );
                    to_edge_obj[i2].push_back( edges.size( ) );
                    auto square = [&]( double x ){ return x*x; };
                    double dist = sqrt(
                         Squarex( M[i2].first.Start( ) - M[i1].first.Stop( ) ) +
                         Squarex( M[i2].second.Start( ) - M[i1].second.Stop( ) ) );
                    edges.push_back(dist);   }
               for ( int i = 0; i < M.isize( ); i++ )
               {    from[N-2].push_back(i), to[i].push_back(N-2);
                    from_edge_obj[N-2].push_back( edges.size( ) );
                    to_edge_obj[i].push_back( edges.size( ) );
                    double sum = Squarex( M[i].first.Start( ) );
                    if ( !optb ) sum += Squarex( M[i].second.Start( ) );
                    edges.push_back( sqrt(sum) );    }
               for ( int i = 0; i < M.isize( ); i++ )
               {    from[i].push_back(N-1), to[N-1].push_back(i);
                    from_edge_obj[i].push_back( edges.size( ) );
                    to_edge_obj[N-1].push_back( edges.size( ) );
                    double sum = Squarex( B.isize( ) - M[i].first.Stop( ) );
                    if ( !optb )
                         sum += Squarex( G.isize( ) - M[i].second.Stop( ) );
                    edges.push_back( sqrt(sum) );    }
               for ( int i = 0; i < N; i++ )
               {    SortSync( from[i], from_edge_obj[i] );
                    SortSync( to[i], to_edge_obj[i] );    }
               digraphE<double> H(from, to, edges, to_edge_obj, from_edge_obj);
               H.ShortestPath( N-2, N-1, path );
               {    to_delete.resize_and_set( M.size( ), True );
                    for ( int i = 1; i < path.isize( ) - 1; i++ )
                         to_delete[ path[i] ] = False;
                    EraseIf( M, to_delete );    }

               /*
               if ( verbosity >= 1 )
               {    xout << "shortest path matches:\n";
                    for ( int i = 0; i < M.isize( ); i++ )
                    {    xout << "[" << i+1 << "] match A[" 
                              << M[i].first << ") to " << "G[" << M[i].second
                              << "), l = " << M[i].first.Length( ) << endl;    }    }
               */

               // Trim overlapping matches (again), then delete subsize match blocks.

               for ( int i = 0; i < M.isize( ) - 1; i++ )
               {    int over = 0;
                    if ( M[i].first.Stop( ) > M[i+1].first.Start( ) )
                         over = M[i].first.Stop( ) - M[i+1].first.Start( );
                    if ( M[i].second.Stop( ) > M[i+1].second.Start( ) )
                    {    over = Max( over, M[i].second.Stop( ) 
                              - M[i+1].second.Start( ) );    }
                    if ( over < M[i].first.Length( ) 
                              && over < M[i].second.Length( )
                              && over < M[i+1].first.Length( ) 
                              && over < M[i+1].second.Length( ) )
                    {    M[i].first.AddToStop(-over);
                         M[i].second.AddToStop(-over);    }    }

               // Print matches.

               if ( verbosity >= 1 )
               {    xout << "matches:\n";
                    for ( int i = 0; i < M.isize( ); i++ )
                    {    xout << "[" << i+1 << "] match A[" 
                              << M[i].first << ") to " << "G[" << M[i].second
                              << "), l = " << M[i].first.Length( ) << endl;    }    }

               // Punt if we don't have a match block near the end.

               Bool punt = ( M.empty( ) || 
                    ( M.front( ).first.Start( ) > MAX_TO_END
                         && M.front( ).second.Start( ) > MAX_TO_END )
                    || ( B.isize( ) - M.back( ).first.Stop( ) > MAX_TO_END
                         && G.isize( ) - M.back( ).second.Stop( ) > MAX_TO_END ) );
               if (punt) 
               {    punted = True;
                    if ( verbosity >= 2 ) 
                    {    out << "punting\n" << xout.str( );
                         PRINT_TO( out, M.front( ).first.Start( ) );
                         PRINT_TO( out, M.front( ).second.Start( ) );
                         PRINT2_TO( out, B.size( ), M.back( ).first.Stop( ) );
                         PRINT2_TO( out, G.size( ), M.back( ).second.Stop( ) );    }
                    i = j - 1;
                    continue;    }

               // Do the easy cases.

               align a;
               Bool easy = True;
               for ( int i = 0; i < M.isize( ) - 1; i++ )
               {    int delta1 = M[i+1].first.Start( ) - M[i].first.Stop( );
                    int delta2 = M[i+1].second.Start( ) - M[i].second.Stop( );
                    if ( delta1 == 0 || delta2 == 0 ) continue;
                    if ( delta1 < MIN_DELTA || delta2 < MIN_DELTA 
                         || delta1 > MAX_DELTA || delta2 > MAX_DELTA )
                    {    easy = False;    
                         if ( verbosity >= 2 )
                         {    xout << "Not easy." << endl;
                              PRINT3_TO( 
                                   xout, delta1, delta2, MAX_DELTA );    }    }    }
               if (easy)
               {    shortvector<mutmer> m( M.size( ) );
                    for ( int l = 0; l < M.isize( ); l++ )
                    {    int e = 0; // actually knowable in advance
                         for ( int r1 = M[l].first.Start( ); 
                              r1 < M[l].first.Stop( ); r1++ )
                         {    int r2 = r1 - M[l].first.Start( )
                                   + M[l].second.Start( );
                              if ( B[r1] != G[r2] ) e++;    }
                         m(l).SetFrom( M[l].first.Start( ), M[l].second.Start( ), 
                              M[l].first.Length( ), e );    }
                    nerrors = 0;
                    a.CreateFromMutmersAndSW( K, m, B, G, MAX_ERRORS,
                         END_STRETCH, nerrors, True );    }
               else 
               {    i = j - 1;
                    if ( verbosity >= 2 ) out << xout.str( );
                    continue;    }

               // Print and save.
     
               if ( verbosity >= 1 ) PrintVisualAlignmentClean(True, xout, B, G, a);
               if ( ( a.pos1( ) > 0 && a.pos2( ) > 0 )
                    || ( a.Pos1( ) < B.isize( ) && a.Pos2( ) < G.isize( ) ) )
               // if ( !a.FullLength( B.size( ) ) )
               {    if ( verbosity >= 2 )
                    {    xout << "Not full length!" << endl;
                         out << xout.str( );    }
                    i = j - 1;
                    continue;    }
               int unaligned = a.pos1( ) + B.isize( ) - a.Pos1( );
               valigns.push( make_pair(unaligned, nerrors), 
                    make_triple( g, ps, start ), xout.str( ) );
               a.AddToStartOnTarget(start);

               // If alignments is one base short, fix it.  I think this
               // same fix may be somewhere else.  Note not adjusting nerrors.
               // Also note that the PrintVisual... above doesn't see this change.

               if ( a.Pos1( ) == B.isize( ) - 1 && a.Pos2( ) < genome[g].isize( ) )
                    a.AddToLength( a.Nblocks( ) - 1, 1 );

               // Save.

               galigns.push_back( refalign( ( ps == 0 ? g+1 : -g-1 ), a ) );

               // Advance.

               i = j - 1;    }    }
     SortSync( valigns, galigns );
     if ( punted && valigns.empty( ) ) out << "\nPUNTED" << endl << endl;
     else if ( valigns.empty( ) ) out << "\nUNALIGNED" << endl << endl;
     int bests = 0;
     for ( int i = 0; i < valigns.isize( ); i++ )
     {    if ( i > 0 && valigns[i].first.first > valigns[i-1].first.first )
          {    bests++;
               if ( !allow_second_best || bests == 2 )
               {    galigns.resize(i);
                    break;    }    }
          if ( i > 0 && valigns[i].first.second > 2 * valigns[i-1].first.second + 5 )
          {    galigns.resize(i);
               break;    }
          int g = valigns[i].second.first;
          String chr;
          if ( g == 22 ) chr = "X";
          else if ( g == 23 ) chr = "Y";
          else chr = ToString(g+1);
          out << "\n" << ( valigns[i].second.second == 0 ? "+" : "-" )
               << chr << ":" << valigns[i].second.third << endl;
          out << valigns[i].third;    }
     report += out.str( );    }

void RefAlign(
     // inputs:
     const vecbasevector& genome, const int HBK, const vecbasevector& tigs,
     const vec<int>& inv, const digraphE<vec<int>>& D, 
     const vec<int>& dinv, const vec<vec<vec<vec<int>>>>& dlines,
     const MasterVec< SerfVec<triple<int,int,int> > >& alignsb,
     // output:
     MasterVec<SerfVec<refalign>>& galigns,
     // control:
     const int verbosity, const vec<int>& test_edges, const Bool allow_second_best )
{
     // Set up.

     double clock = WallClockTime( );
     vec<int> dlens( D.E( ), 0 );
     vec<Bool> used;
     D.Used(used);
     #pragma omp parallel for
     for ( int e = 0; e < D.E( ); e++ )
     {    if ( !used[e] || D.O(e)[0] < 0 ) continue;
          for ( int j = 0; j < D.O(e).isize( ); j++ )
               dlens[e] += tigs[ D.O(e)[j] ].isize( ) - HBK + 1;    }
     vec<int> to_left, to_right;
     D.ToLeft(to_left), D.ToRight(to_right);

     // Heuristics.

     const int K = 40;
     const int MIN_LEN = 25;
     const double MAX_BUB_DIFF = 0.5;

     // Define cat function, depends on HBK and tigs.

     auto Cat = [&]( const vec<int>& e )
     {    basevector x = tigs[ e[0] ];
          for ( int j = 1; j < e.isize( ); j++ )
               x = TrimCat( HBK, x, tigs[ e[j] ] );
          return x;    };

     // Select edges.  We don't have to align identical edges or reverse
     // complement edges.

     cout << Date( ) << ": selecting edges" << endl;
     vec<vec<int>> all;
     vec<int> ds;
     if ( test_edges.nonempty( ) )
     {    for ( int i = 0; i < test_edges.isize( ); i++ )
          {    int d = test_edges[i];
               if ( dlens[d] >= MIN_LEN )
               {    all.push_back( D.O(d) );
                    ds.push_back(d);    }    }    }
     else 
     {    for ( int d = 0; d < D.E( ); d++ )
          {    if ( dlens[d] >= MIN_LEN ) 
               {    if ( D.O(d) <= D.O( dinv[d] ) )
                    {    all.push_back( D.O(d) );
                         ds.push_back(d);    }    }    }    }
     ParallelUniqueSortSync( all, ds );

     // Output data structure galigns = {(chr,start,stop,align)},
     // where chr = 1+g, and negative if rc.

     galigns.clear( );
     galigns.resize( D.E( ) );

     // Align edges to ref.

     cout << Date( ) << ": aligning edges to ref" << endl;
     if ( verbosity >= 1 ) cout << endl;
     vec<String> reports( ds.size( ) );
     const int mismatch_penalty = 3;
     const int gap_open_penalty = 12;
     const int gap_extend_penalty = 1;
     const int batch = 1;
     #pragma omp parallel for schedule(dynamic, 1)
     for ( int bi = 0; bi < ds.isize( ); bi += batch )
     {    RefAlignScratch<K> scr;
          for ( int s = bi; s < Min( bi + batch, ds.isize( ) ); s++ )
          {    int d = ds[s];
               if ( test_edges.nonempty( ) && !Member( test_edges,d ) ) continue;
               ostringstream out;
               out << "========================================================"
                    << "========================\n";
               out << "\n[" << s+1 << "] d = " << d << "[l=" << dlens[d] << "]" 
                    << endl;
               reports[s] = out.str( );
               RefAlignCore<K>( D.O(d), HBK, tigs, inv, D, dinv, dlines, genome, 
                    alignsb, galigns[d], verbosity, reports[s], scr, -1,
                    allow_second_best );

               // Print summary line.

               if ( verbosity >= 1 && galigns[d].nonempty( ) )
               {    const refalign& r = galigns[d][0];
                    int chr = r.chr;
                    if ( chr >= 0 )
                    {    const align& a = r.a;
                         const basevector B = Cat( D.O(d) );
                         const basevector& G = genome[chr-1];
                         int p1 = a.pos1( ), p2 = a.pos2( );
                         int errs = 0;
                         for ( int j = 0; j < a.Nblocks( ); j++ )
                         {    if ( a.Gaps(j) > 0 ) p2 += a.Gaps(j);
                              if ( a.Gaps(j) < 0 ) p1 -= a.Gaps(j);
                              if ( a.Gaps(j) != 0 )
                              {    errs += gap_open_penalty 
                                        + gap_extend_penalty * Abs( a.Gaps(j) );    }
                              for ( int x = 0; x < a.Lengths(j); x++ )
                              {    if ( B[p1] != G[p2] ) errs += mismatch_penalty;
                                   ++p1; ++p2;    }    }
                         reports[s] += "ALIGNED " + ToString(d)
                              + " with " + ToString(errs) + " errors\n\n";    }
                              }    }    }

     // Fill in alignments.

     cout << Date( ) << ": filling in alignments" << endl;
     #pragma omp parallel for
     for ( int d = 0; d < D.E( ); d++ )
     {    if ( dlens[d] < MIN_LEN || BinMember( ds, d ) ) continue;
          {    if ( D.O(d) <= D.O( dinv[d] ) )
               {    int p = BinPosition( all, D.O(d) );
                    if ( p >= 0 ) galigns[d] = galigns[ ds[p] ];    }
               else
               {    int p = BinPosition( all, D.O( dinv[d] ) );
                    if ( p >= 0 )
                    {    galigns[d] = galigns[ ds[p] ];
                         for ( int j = 0; j < (int) galigns[d].size( ); j++ )
                              galigns[d][j].chr *= -1;    }    }    }    }

     // Destroy all since it isn't used beyond this
     Destroy( all );

     // Try to fill in bubbles.

     cout << Date( ) << ": filling in bubbles" << endl;
     #pragma omp parallel for
     for ( int d = 0; d < D.E( ); d++ )
     {    if ( !used[d] || D.O(d)[0] < 0 ) continue;
          int v = to_left[d], w = to_right[d];
          if ( !D.To(v).solo( ) || !D.From(w).solo( ) ) continue;
          int d1 = D.ITo(v,0), d2 = D.IFrom(w,0);
          if ( galigns[d1].size( ) != 1 || galigns[d2].size( ) != 1 ) continue;
          int chr1 = galigns[d1][0].chr, chr2 = galigns[d2][0].chr;
          if ( chr1 < 0 || chr1 != chr2 ) continue;
          if ( galigns[d].size( ) > 1 ) continue;
          if ( galigns[d].size( ) == 1 && galigns[d][0].chr == chr1
               && galigns[d1][0].a.Pos2( ) - galigns[d][0].a.pos2( ) == HBK-1
               && galigns[d][0].a.Pos2( ) - galigns[d2][0].a.pos2( ) == HBK-1 )
          {    continue;    }
          int rd = dinv[d], rd1 = dinv[d1], rd2 = dinv[d2];
          if ( rd == d ) continue;
          if ( galigns[rd1].size( ) != 1 || galigns[rd2].size( ) != 1 ) continue;
          int rchr1 = galigns[rd1][0].chr, rchr2 = galigns[rd2][0].chr;
          if ( rchr1 != -chr1 || rchr1 != rchr2 ) continue;
          int g = chr1 - 1, K = HBK;
          const align &a1 = galigns[d1][0].a, &a2 = galigns[d2][0].a;
          int gstart = a1.Pos2( ) - (K-1), gstop = a2.pos2( ) + (K-1);
          if ( gstart < 0 || gstop >= genome[g].isize( ) ) continue;
          int glen = gstop - gstart, alen = dlens[d] + (K-1);
          if ( glen <= 0 || Abs(glen-alen) > alen * MAX_BUB_DIFF ) continue;
          basevector B = Cat( D.O(d) ), G( genome[g], gstart, gstop - gstart );
          align a;
          int nerrors;
          int offset =  (alen - glen) / 2, bandwidth = Abs(glen-alen) + 2;
          SmithWatAffineBanded( B, G, offset, bandwidth, a, nerrors );
          if ( verbosity >= 1 ) 
          {    
               #pragma omp critical
               {    cout << "bubble attempted, ";
                    PRINT5( 
                         d, a.pos1( ), a.Pos1( ), a.pos2( ), a.Pos2( ) );    }    }
          if ( a.pos1( ) != 0 || a.Pos1( ) != B.isize( ) ||
               a.pos2( ) != 0 || a.Pos2( ) != G.isize( ) )
          {    int best_loc;
               alignment al;
               if ( B.size( ) <= G.size( ) )
               {    SmithWatFree( B, G, best_loc, al, true, true, 2, 3, 1000000 );
                    a = al;    }
               else
               {    SmithWatFree( G, B, best_loc, al, true, true, 2, 3, 1000000 );
                    a = al;
                    a.Flip( );    }    }
          if ( a.pos1( ) != 0 || a.Pos1( ) != B.isize( ) ) continue;
          if ( a.pos2( ) != 0 || a.Pos2( ) != G.isize( ) ) continue;
          a.AddToStartOnTarget(gstart);
          if ( verbosity >= 1 ) 
          {    
               #pragma omp critical
               {    cout << "bubble filled, ";
                    PRINT3( d, a.pos2( ), a.Pos2( ) );    }    }
          galigns[d].clear( ), galigns[rd].clear( );
          galigns[d].push_back( refalign( chr1, a ) );
          galigns[rd].push_back( refalign( -chr1, a ) );    }
     cout << Date( ) << ": done" << endl;

     // Try to fill in cells, and chains of cells.

     cout << Date( ) << ": filling in cells" << endl;
     const int MAX_CELL_DIAMETER = 1000;
     vec< triple< int, int, ho_interval > > locs;
     vec<Bool> locating( D.E( ), False );
     for ( int dl = 0; dl < dlines.isize( ); dl++ )
     {    const vec<vec<vec<int>>>& L = dlines[dl];

          // Define anchored endpoints of cells.

          vec<pair<int,int>> left( L.size( ), make_pair(-1,-1) );
          vec<pair<int,int>> right( L.size( ), make_pair(-1,-1) );
          for ( int j = 0; j < L.isize( ); j++ )
          {    Bool left_nonsolo = False, right_nonsolo = False;
               vec<pair<int,int>> lefts, rights;
               for ( int k = 0; k < L[j].isize( ); k++ )
               {    if ( L[j][k].empty( ) ) continue;
                    int d = L[j][k].front( );
                    if ( !galigns[d].solo( ) ) left_nonsolo = True;
                    else
                    {    const refalign& r = galigns[d][0];
                         lefts.push( r.chr, r.a.pos2( ) );    }
                    d = L[j][k].back( );
                    if ( !galigns[d].solo( ) ) right_nonsolo = True;
                    else
                    {    const refalign& r = galigns[d][0];
                         rights.push( r.chr, r.a.Pos2( ) );    }    }
               UniqueSort(lefts), UniqueSort(rights);
               if ( !left_nonsolo && lefts.solo( ) && lefts[0].first >= 0 ) 
                    left[j] = lefts[0];
               if ( !right_nonsolo && rights.solo( ) && rights[0].first >= 0 ) 
                    right[j] = rights[0];    }

          // Look for edges that might be aligned.

          for ( int j = 0; j < L.isize( ); j++ )
          {    const vec<vec<int>>& C = L[j];
               pair<int,int> left_bound = left[j], right_bound = right[j];
               Bool left_is_left = True, right_is_right = True;
               int jl = j, jr = j;
               while( left_bound.first < 0 )
               {    jl--;
                    if ( jl < 0 ) break;
                    if ( right[jl].first >= 0 )
                    {    left_bound = right[jl];
                         left_is_left = False;
                         break;    }
                    if ( left[jl].first >= 0 )
                    {    left_bound = left[jl];
                         break;    }    }
               while( right_bound.first < 0 )
               {    jr++;
                    if ( jr == L.isize( ) ) break;
                    if ( left[jr].first >= 0 )
                    {    right_bound = left[jr];
                         right_is_right = False;
                         break;    }
                    if ( right[jr].first >= 0 )
                    {    right_bound = right[jr];
                         break;    }    }
               int chr = left_bound.first;
               if ( chr < 0 || chr != right_bound.first ) continue;
               int start = left_bound.second;
               if ( !left_is_left ) start -= (HBK-1);
               int stop = right_bound.second;
               if ( !right_is_right) stop += (HBK-1);
               int sep = stop - start;
               if ( sep <= 0 || sep > MAX_CELL_DIAMETER ) continue;
               vec<int> ds = Contents(C);
               for ( auto d : ds )
               {    if ( galigns[d].size( ) > 0 || D.O(d)[0] < 0 ) continue;
                    locs.push( d, chr, ho_interval(start,stop) );
                    locating[d] = True;    }    }    }
     #pragma omp parallel for
     for ( int u = 0; u < locs.isize( ); u++ )
     {    int d = locs[u].first, chr = locs[u].second;
          int start = locs[u].third.Start( ), stop = locs[u].third.Stop( );
          basevector G( genome[chr-1], start, stop - start );
          int rd = dinv[d];
          if ( locating[rd] || galigns[rd].size( ) > 0 ) continue;
          basevector B = Cat( D.O(d) );
          alignment al;
          Bool left_fixed = False, right_fixed = False;
          int v = to_left[d], w = to_right[d];
          if ( D.To(v).solo( ) )
          {    int f = D.ITo(v,0);
               if ( galigns[f].solo( ) ) left_fixed = True;    }
          if ( D.From(w).solo( ) )
          {    int f = D.IFrom(w,0);
               if ( galigns[f].solo( ) ) right_fixed = True;    }
          int best_loc;
          if ( B.isize( ) <= G.isize( ) )
          {    SmithWatFree( B, G, best_loc, al, left_fixed, right_fixed );
               if ( al.pos1( ) > 0 || al.Pos1( ) < B.isize( ) ) continue;
               start += al.pos2( );
               stop -= ( G.isize( ) - al.Pos2( ) );    }
          else
          {    SmithWatFree( G, B, best_loc, al );
               if ( al.pos2( ) > 0 || al.Pos2( ) < B.isize( ) ) continue;
               start += al.pos1( );
               stop -= ( G.isize( ) - al.Pos1( ) );    }
          if ( start >= stop ) continue;
          basevector G2( genome[chr-1], start, stop - start );
          SmithWatAffine( B, G2, al, left_fixed, right_fixed );
          if ( al.pos1( ) > 0 || al.Pos1( ) < B.isize( ) ) continue;
          align a(al);
          a.AddToStartOnTarget(start);
          if ( verbosity >= 1 ) 
          {    
               #pragma omp critical
               {    cout << "cell filled, ";
                    PRINT3( d, a.pos2( ), a.Pos2( ) );    }    }
          galigns[d].push_back( refalign( chr, a ) );
          galigns[rd].push_back( refalign( -chr, a ) );    }

     // Output reports.

     if ( verbosity >= 1 )
     {    cout << endl;
          for ( int s = 0; s < ds.isize( ); s++ ) cout << reports[s];
          cout << "========================================================"
               << "========================\n\n";    }

     // Report stats.

     if ( verbosity >= 1 )
     {    int bins = 4;
          vec<int> total(bins,0), aligned(bins,0), multiply(bins,0);
          for ( int d = 0; d < D.E( ); d++ )
          {    if ( !used[d] || D.O(d)[0] < 0 ) continue;
               if ( dlens[d] < 250 )
               {    total[0]++;
                    if ( galigns[d].size( ) > 0 ) aligned[0]++;
                    if ( galigns[d].size( ) > 1 ) multiply[0]++;    }
               else if ( dlens[d] < 500 )
               {    total[1]++;
                    if ( galigns[d].size( ) > 0 ) aligned[1]++;
                    if ( galigns[d].size( ) > 1 ) multiply[1]++;    }
               else if ( dlens[d] < 1000 )
               {    total[2]++;
                    if ( galigns[d].size( ) > 0 ) aligned[2]++;
                    if ( galigns[d].size( ) > 1 ) multiply[2]++;    }
               else
               {    total[3]++;
                    if ( galigns[d].size( ) > 0 ) aligned[3]++;
                    if ( galigns[d].size( ) > 1 ) multiply[3]++;    }    }
          auto per = [&]( double n, double d )
          {    ostringstream out;
               out << fixed << setprecision(2) << setw(5) << right << 100*n/d;
               return out.str( );    };
          vec<vec<String>> rows;
          cout << "\n// SUMMARY OF REFERENCE ALIGNMENTS" << endl;
          vec<String> titles 
               = { "// 1 - 249", "// 250-499", "// 500-999", "// 1000 -" };
          rows.push( vec<String>( 
               { "// kmers", "count", "%aligned", "%multiply" } ) );
          for ( int j = 0; j < bins; j++ )
          {    rows.push( vec<String>( { titles[j], ToStringAddCommas(total[j]),
                    per(aligned[j],total[j]), per(multiply[j],total[j]) } ) );    }
          PrintTabular( cout, rows, 2, "lrrr" );

          // Test for compatibility.  Note that we only test to see if the reference
          // intervals overlap sanely.  We do not look at the alignments themselves,
          // beyond the intervals that they define.

          int looked = 0, compat = 0;
          vec< pair<int,int> > incompats;
          vec<Bool> bad( D.E( ), False );
          for ( int d1 = 0; d1 < D.E( ); d1++ )
          {    if ( galigns[d1].size( ) != 1 ) continue;
               int v = to_right[d1];
               for ( int j = 0; j < D.From(v).isize( ); j++ )
               {    int d2 = D.IFrom(v,j);
                    if ( galigns[d2].size( ) != 1 ) continue;
                    int chr1 = galigns[d1][0].chr, chr2 = galigns[d2][0].chr;
                    align a1 = galigns[d1][0].a, a2 = galigns[d2][0].a;
                    int f1 = d1, f2 = d2;
                    Bool compatible = False;
                    looked++;
                    if ( chr1 == chr2 )
                    {    if ( chr1 < 0 )
                         {    swap( a1, a2 );
                              swap( f1, f2 );
                              f1 = dinv[f1], f2 = dinv[f2];    }
                         if ( a1.Pos2( ) - a2.pos2( ) == HBK - 1 )
                              compatible = True;    }
                    if (compatible) compat++;
                    else 
                    {    incompats.push( d1, d2 );    
                         bad[d1] = bad[d2] = True;    }    }    }
          cout << "// weakly compatible abutting solo placements: " 
               << per(compat, looked) << "%" << endl;
          vec<int> unaligned;
          for ( int d = 0; d < D.E( ); d++ )
          {    if ( !used[d] || D.O(d)[0] < 0 ) continue;
               if ( galigns[d].size( ) > 0 ) continue;
               unaligned.push_back(d);
               bad[d] = True;    }
          int64_t total_kmers = 0, good = 0;
          for ( int d = 0; d < D.E( ); d++ )
          {    total_kmers += dlens[d];
               if ( !bad[d] ) good += dlens[d];    }
          cout << "// kmers aligned uniquely and consistently: "
               << PERCENT_RATIO( 3, good, total_kmers ) << endl;

          // Compute alignment error rate.

          int64_t total_bases = 0, total_errs = 0;
          #pragma omp parallel for schedule(dynamic, 1000)
          for ( int d = 0; d < D.E( ); d++ )
          {    if ( galigns[d].empty( ) ) continue;
               const refalign& r = galigns[d][0];
               int chr = r.chr;
               if ( chr < 0 ) continue;
               const align& a = r.a;
               const basevector B = Cat( D.O(d) );
               const basevector& G = genome[chr-1];
               int p1 = a.pos1( ), p2 = a.pos2( );
               int errs = 0;
               for ( int j = 0; j < a.Nblocks( ); j++ )
               {    if ( a.Gaps(j) > 0 ) p2 += a.Gaps(j);
                    if ( a.Gaps(j) < 0 ) p1 -= a.Gaps(j);
                    if ( a.Gaps(j) != 0 )
                    {    errs += gap_open_penalty 
                              + gap_extend_penalty * Abs( a.Gaps(j) );    }
                    for ( int x = 0; x < a.Lengths(j); x++ )
                    {    if ( B[p1] != G[p2] ) errs += mismatch_penalty;
                         ++p1; ++p2;    }    }
               #pragma omp critical
               {    total_errs += errs;
                    total_bases += B.size( );    }    }
          ostringstream out;
          out << fixed << setprecision(2) 
               << ( 1000.0 * total_errs ) / (mismatch_penalty * total_bases) << endl;
          cout << "// errs per kb = "  << out.str( );    }

     // Done.

     if ( verbosity >= 1 )
     {    cout << "// time used = " << TimeSince(clock)
               << ", peak mem = " << PeakMemUsageGBString( ) 
               << endl << endl;    }    }

template void RefAlignCore(
     const vec<int>& X,
     const int HBK,
     const vecbasevector& tigs,
     const vec<int>& inv,
     const digraphE<vec<int>>& D,
     const vec<int>& dinv,
     const vec<vec<vec<vec<int>>>>& dlines,
     const vecbasevector& genome,
     const MasterVec< SerfVec<triple<int,int,int> > >& alignsb,
     SerfVec<refalign>& galigns,
     const int verbosity,
     String& report,
     RefAlignScratch<40>& s,
     const int g_use,
     const Bool allow_second_best,
     const Bool optb );

