// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

// This code aligns finished sequence to an assembly and reports results.  The
// method is as follows:
//
// 1. First find perfect matches of the finished sequence to the assembly.  We
// exclude short perfect matches, and "dominated" matches, unless the dominated
// match is long.  Where there are two perfect matches that have the same start
// and stop points in the assembly and on the reference, we pick one arbitrarily.
//
// 2. We then trim back the ends of the perfect matches.  The goal is to connect
// up the perfect matches, but the ends can overlap, so they have to be trimmed.
//
// 3. We find paths from one perfect path to another, then align, allowing for
// captured gaps in the middle.  We save the best alignment and its score.
//
// 4. We find the shortest path, extending from end to end on the finished sequence,
// containing the perfect paths and the imperfect paths between them, as well as
// uncaptured gaps.
//
// 5. We then build a composite alignment by gluing together the perfect and 
// imperfect alignments.
//
// 6. A visual alignment and statistics are then read off the composite alignment.
//
// WARNING.  This module creates files o.report, o.fullreport and o.finpaths
// each time your run it.  You might accidentally overwrite these files using
// your favorite value for F.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"
#include "PackAlign.h"
#include "PrintAlignment.h"
#include "math/HoInterval.h"
#include "pairwise_aligners/SmithWatAffine.h"
#include "paths/long/large/Lines.h"
#include "10X/Gap.h"
#include "10X/Heuristics.h"
#include "10X/astats/FinAlign.h"
#include "10X/astats/ProjectRunway.h"
#include "random/Shuffle.h"
#include "10X/DfTools.h"

// Alignment structures.  An malign represents a general alignment of a clone to
// an assembly.

class aalign {

     public:

     vec<int> p;        // path through assembly
     ho_interval gpos;  // interval on reference matching
                        // (above goofy, should just roll into a.pos2)
     align a;           // alignment of assembly to reference
     int score;         // alignment score

};

class malign {

     public:

     malign( ) { }
     vec<aalign> aligns;  // alignments of paths
     vec<vec<int>> gaps;  // gap "path" or empty if uncaptured

};

// Version of affine Smith-Waterman alignment that captures left and right gaps.

int SmithWatAffineEndFix( const basevector& S, const basevector& T, align& al,
     bool penalize_left_gap = true, bool penalize_right_gap = true,
     const int mismatch_penalty = 3, const int gap_open_penalty = 12,
     const int gap_extend_penalty = 1 )
{    alignment a;
     int score = SmithWatAffine( S, T, a, penalize_left_gap, penalize_right_gap,
          mismatch_penalty, gap_open_penalty, gap_extend_penalty );
     al = a;
     if (penalize_left_gap)
     {    int left_hang = al.pos2( );
          if ( left_hang > 0 )
          {    int n = al.Nblocks( );
               al.SetNblocks( n + 1 );
               for ( int j = n - 1; j >= 0; j-- )
               {    al.SetGap( j+1, al.Gaps(j) );
                    al.SetLength( j+1, al.Lengths(j) );    }
               al.SetGap( 0, left_hang );
               al.SetLength( 0, 0 );
               al.Setpos2(0);    }
          left_hang = al.pos1( );
          if ( left_hang > 0 )
          {    int n = al.Nblocks( );
               al.SetNblocks( n + 1 );
               for ( int j = n - 1; j >= 0; j-- )
               {    al.SetGap( j+1, al.Gaps(j) );
                    al.SetLength( j+1, al.Lengths(j) );    }
               al.SetGap( 0, -left_hang );
               al.SetLength( 0, 0 );
               al.Setpos1(0);    }    }
     if (penalize_right_gap)
     {    int right_hang = T.isize( ) - al.Pos2( );
          if ( right_hang > 0 )
          {    int n = al.Nblocks( );
               al.SetNblocks( n + 1 );
               al.SetGap( n, right_hang );
               al.SetLength( n, 0 );    }
          right_hang = S.isize( ) - al.Pos1( );
          if ( right_hang > 0 )
          {    int n = al.Nblocks( );
               al.SetNblocks( n + 1 );
               al.SetGap( n, -right_hang );
               al.SetLength( n, 0 );    }    }
     return score;    }

int BuildGapConnection( malign& ma, const int d1, const int d2, 
     const vecbasevector& tigs,
     const vec<int>& left, const vec<int>& right, const int aleft, 
     const int aright, const int g, const int gstop1, const int gright,
     const vec<int>& gpath, const HyperBasevectorX& hb, const digraphE<vec<int>>& D,
     const vecbasevector& genomef,
     const int MISMATCH_PENALTY, const int GAP_OPEN_PENALTY,
     const int GAP_EXTEND_PENALTY, const int GAP_SLACK, const int MIN_ADD,
     const Bool VERBOSE, ostringstream& out )
{    
     int K = hb.K( );

     // Define a function.

     auto Cat = [&]( const vec<int>& e )
     {    basevector x = tigs[ e[0] ];
          for ( int j = 1; j < e.isize( ); j++ )
               x = TrimCat( K, x, tigs[ e[j] ] );
          return x;    };

     basevector A, G, b1 = Cat( D.O(d1) );
     A.SetToSubOf( b1, aleft, b1.isize( ) - aleft );
     for ( int m = 0; m < left.isize( ); m++ )
     {    basevector c = Cat( D.O( left[m] ) );
          c.SetToSubOf( c, K-1, c.isize( ) - (K-1) );
          A.append(c);    }
     int gleft1 = gstop1;
     int minlen = Max( GAP_SLACK, A.isize( ) + MIN_ADD );
     int gleft2 = Min( gstop1 + minlen, genomef[g].isize( ) );
     G.SetToSubOf( genomef[g], gleft1, gleft2 - gleft1 );
                              
     // Affine align A to G, anchored on the left, with right ends free.

     aalign aa1;
     aa1.p = {d1};
     aa1.p.append(left);
     int lscore;
     if ( A.size( ) == 0 ) // align suppressed later
     {    lscore = 0;
          avector<int> gaps, lengths;
          aa1.a.Set( 0, 0, gaps, lengths );    }
     else if ( G.size( ) == 0 )
     {    lscore = GAP_OPEN_PENALTY + ( A.size( ) - 1 ) * GAP_EXTEND_PENALTY;    
          avector<int> gaps(1), lengths(1);
          gaps(0) = A.size( ), lengths(0) = 0;
          aa1.a.Set( 0, A.size( ), gaps, lengths );    }
     else
     {    align a;
          if (VERBOSE)
          {    out << "aligning " << G.ToString( ) << " to " << A.ToString( ) 
                    << endl;    }
          // note weird call to A, G to bypass bug (?)
          lscore = SmithWatAffineEndFix( A, G, a, true, false, MISMATCH_PENALTY, 
               GAP_OPEN_PENALTY, GAP_EXTEND_PENALTY );
          a.Flip( );
          aa1.a = a;    }
     aa1.gpos.Set( gleft1, gleft1 + aa1.a.Pos1( ) );
     aa1.score = lscore;
     aa1.a.AddToPos2(aleft);
     if (VERBOSE) out << "left score = " << lscore/10.0 << endl;
     A = Cat( D.O(d1) );
     for ( int m = 0; m < left.isize( ); m++ )
     {    basevector c = Cat( D.O( left[m] ) );
          c.SetToSubOf( c, K-1, c.isize( ) - (K-1) );
          A.append(c);    }
     if (VERBOSE) PRINT_TO( out, A.size( ) );

     // Now do the right side.

     A.clear( );
     for ( int m = 0; m < right.isize( ); m++ )
     {    A.append( Cat( D.O( right[m] ) ) );
          A.resize( A.isize( ) - (K-1) );    }
     A.append( basevector( Cat( D.O(d2) ), 0, aright ) );
     minlen = Max( GAP_SLACK, A.isize( ) + MIN_ADD );
     int gright1 = Max( gright - minlen, 0 ), gright2 = gright;
     if (VERBOSE) PRINT2_TO( out, gright1, gright2 );
     G.SetToSubOf( genomef[g], gright1, gright2 - gright1 );

     // Affine align A to G, anchored on the right, with left ends free.

     aalign aa2;
     aa2.p = right;
     aa2.p.push_back(d2);
     int rscore;
     if ( A.size( ) == 0 ) // align suppressed later
     {    rscore = 0;
          avector<int> gaps, lengths;
          aa2.a.Set( 0, 0, gaps, lengths );    }
     else if ( G.size( ) == 0 )
     {    rscore = GAP_OPEN_PENALTY + ( A.size( ) - 1 ) * GAP_EXTEND_PENALTY;
          avector<int> gaps(1), lengths(1);
          gaps(0) = A.size( ), lengths(0) = 0;
          aa2.a.Set( 0, A.size( ), gaps, lengths );    }
     else
     {    align a;
          if (VERBOSE)
          {    out << "aligning " << G.ToString( ) << " to " << A.ToString( ) 
                    << endl;    }
          // note weird call to A, G to bypass bug (?)
          rscore = SmithWatAffineEndFix( A, G, a, false, true, MISMATCH_PENALTY, 
               GAP_OPEN_PENALTY, GAP_EXTEND_PENALTY );    
          a.Flip( );
          if (VERBOSE)
          {    out << "initial alignment" << endl;
               PRINT2_TO( out, a.pos1( ), a.pos2( ) );
               PrintVisualAlignmentClean( True, out, G, A, a );    }
          aa2.a = a;    }
     aa2.gpos.Set( gright1 + aa2.a.pos1( ), gright2 );
     aa2.a.Setpos1(0);
     if (VERBOSE) 
     {    PRINT2_TO( out, aa2.gpos.Start( ), aa2.gpos.Stop( ) );
          PRINT2_TO( out, aa2.a.pos2( ), aa2.a.Pos2( ) );
          out << "right score = " << rscore/10.0 << endl;    }
     aa2.score = rscore;

     // Package aligns.

     ma.aligns = {aa1,aa2};
     ma.gaps = {gpath};
     return lscore + rscore;    }

String FinAlign( const String& suffix, const vec<int>& gs, const Bool VERBOSE, 
     const String& INDIR, const HyperBasevectorX& hb, vec<int> inv, 
     vecbasevector tigs, digraphE<vec<int>> D, vec<int> dinv,
     vec<vec<vec<vec<int>>>> dlines, const vecbasevector& genomef,
     const MasterVec< SerfVec< triple< ho_interval, int, int > > >& flocs,
     vec<vec<int>>& statslist, const Bool write_files, const Bool PS )
{    
     // Make lines if needed.

     if ( dlines.empty( ) ) 
          FindLines( D, dinv, dlines, MAX_CELL_PATHS, MAX_CELL_DEPTH );

     // Reinsert loops, then munch.

     const int K = 48;
     ForceAssertEq( K, hb.K( ) );
     int nd = D.E( ), TE = tigs.size( );
     ReinsertLoops( hb, inv, D, dinv );
     cout << Date( ) << ": munching" << endl;
     Munch( D, dinv, tigs, inv, K );

     // Compute ancillary stuff.

     cout << Date( ) << ": computing ancillary stuff" << endl;
     vec<int> kmers( tigs.size( ) ), to_left, to_right, dlens( D.E( ), 0 );
     for ( int e = 0; e < (int) tigs.size( ); e++ )
          kmers[e] = tigs[e].isize( ) - hb.K( ) + 1;
     D.ToLeft(to_left), D.ToRight(to_right);
     #pragma omp parallel for
     for ( int e = 0; e < D.E( ); e++ )
     {    if ( D.O(e)[0] < 0 ) continue;
          for ( int j = 0; j < D.O(e).isize( ); j++ )
               dlens[e] += kmers[ D.O(e)[j] ];    }
     digraphE<int> M( D, dlens );
     vec<int> tol( D.E( ), -1 );
     for ( int i = 0; i < dlines.isize( ); i++ )
     for ( int j = 0; j < dlines[i].isize( ); j++ )
     for ( int k = 0; k < dlines[i][j].isize( ); k++ )
     for ( int l = 0; l < dlines[i][j][k].isize( ); l++ )
          tol[ dlines[i][j][k][l] ] = i;

     // Penalty structure.  After division by ten, penalties are as follows.
     // A mismatch counts as 1, as does a gap open.  A gap extension counts as 0.5,
     // but in the final scoring is computed as 0.1.  A captured gap counts as 100.
     // An uncaptured gap counts as 1000, plus the penalty for the associated gap.  
     // The penalty for a gap in alignment at the end is the number of bases in 
     // the gap.

     const int MISMATCH_PENALTY = 10;
     const int GAP_OPEN_PENALTY = 10;
     const int GAP_EXTEND_PENALTY = 5;
     const int GAP_EXTEND_PENALTY_SHOW = 1;
     const int CAP_PENALTY = 1000;
     const int UNCAP_PENALTY = 10000;
     const int END_PENALTY = 10;

     // Other heuristics.

     const int MINLEN = 150;
     const int MINLEN_FORCE = 500;
     const int MAX_OVERLAP = 250;
     const int MAX_GAP = 350;
     const int MAX_PATHS = 250;
     const int MAX_ITERATIONS = 2000;
     const int MAX_COPIES = 2;
     const int MAX_TRIM = 150;
     const int GAP_SLACK = 200;
     const int MAX_LENGTH = 1000;
     const int MIN_ADD = 100;
     const int MIN_UNCAP = -100;
     const int MAX_UNCAP = 200;
     const int MAX_MOVE = 2;

     // Compute maxperfs.  This gives us a collection of perfect matches of the
     // clone sequence to the assembly.  Note that this will see cycles but not
     // sequence gap edges.

     cout << Date( ) << ": computing maxperfs" << endl;
     double clock = WallClockTime( );
     ostringstream out;
     vec<String> vars;
     vec<vec<vec<triple<ho_interval,int,int>>>> X;
     GetMaxPerfs( gs, hb, kmers, D, to_left, to_right, dlens, tigs,
          genomef, flocs, MINLEN, MINLEN_FORCE, X );
     cout << Date( ) << ": done, time used = " << TimeSince(clock) << endl;
     
     // Data structures for N50_perf and weighted error rate

     vec<vec<ho_interval>> perfs( gs.size( ) );
     vec< pair<int,int> > werr( gs.size( ), make_pair(0,0) );

     // Go through the finished clones.

     int64_t total_score = 0, total_caps = 0, total_uncaps = 0, total_endgaps = 0;
     int64_t total_cov = 0, total_total = 0, total_edges = 0, total_lines = 0;
     vec<quad<int,int,int,int>> stats( gs.size( ), make_quad(0,0,0,0) );
     vec<pair<int,int>> covstats( gs.size( ), make_pair(0,0) );
     vec<pair<int,int>> counts( gs.size( ), make_pair(0,0) );
     int anoms = 0, fails = 0, problems = 0;
     vec<String> reports( gs.size( ) ), fullreports( gs.size( ) );
     vec<vec<int>> allpaths( gs.size( ) );
     vec<malign> allmaligns( gs.size( ) );
     vec<vec<vec<int>>> allpaths_flat( gs.size( ) );
     statslist.clear( );
     statslist.resize( gs.size( ) );
     #pragma omp parallel for schedule(dynamic, 1)
     for ( int i = 0; i < gs.isize( ); i++ )
     {    int g = gs[i];

          // Announce.

          ostringstream out;
          out << "\nCLONE " << g << endl << endl;
          for ( int j = 0; j < X[i].isize( ); j++ )
          {    vec<int> p;
               for ( auto m : X[i][j] ) p.push_back( m.second );
               int gstart = X[i][j].front( ).first.Start( );
               int gstop = X[i][j].back( ).first.Stop( );
               out << "[" << j+1 << "] " << gstart << "-" << gstop;
               if ( gstop == genomef[g].isize( ) ) out << " = end";
               out << ": ";
               if ( p.size( ) <= 3 ) out << printSeq(p) << endl;
               else out << p.front( ) << ".." << p.back( ) << endl;    }
          out << endl;

          // Trim.  We remove up to MAX_TRIM kmers from the ends of each perfect
          // placement, unless it goes up to an end of the clone.

          vec<vec<triple<ho_interval,int,int>>> xleft( X[i].size( ) );
          vec<vec<triple<ho_interval,int,int>>> xright( X[i].size( ) );
          for ( int j = 0; j < X[i].isize( ); j++ )
          {    int left = 0, right = 0;
               while( X[i][j].size( ) > 1 )
               {    int gstart = X[i][j].front( ).first.Start( );
                    if ( gstart == 0 ) break;
                    int d = X[i][j].front( ).second, dstart = X[i][j].front( ).third;
                    left += dlens[d] - dstart;
                    if ( left > MAX_TRIM ) break;
                    xleft[j].push_back( X[i][j].front( ) );
                    X[i][j].pop_front( );    }
               while( X[i][j].size( ) > 1 )
               {    int gstop = X[i][j].back( ).first.Stop( );
                    if ( gstop == genomef[g].isize( ) ) break;
                    int d = X[i][j].back( ).second, dstart = X[i][j].back( ).third;
                    right += dlens[d] - dstart;
                    if ( right > MAX_TRIM ) break;
                    xright[j].push_front( X[i][j].back( ) );
                    X[i][j].pop_back( );    }    }

          // Extra trimming step.

          for ( int j1 = 0; j1 < X[i].isize( ); j1++ )
          for ( int j2 = 0; j2 < X[i].isize( ); j2++ )
          {    if ( j1 == j2 ) continue;
               int gstop1 = X[i][j1].back( ).first.Stop( );
               int gstart2 = X[i][j2].front( ).first.Start( );
               int gap = gstart2 - gstop1;
               if ( gap > MAX_GAP || -gap > MAX_OVERLAP ) continue;
               int d1 = X[i][j1].back( ).second, d2 = X[i][j2].front( ).second;
               if ( d1 != d2 )
               {    if ( to_left[d1] == to_left[d2] || to_right[d1] == to_right[d2] )
                    {    if ( !X[i][j1].solo( ) && !X[i][j2].solo( ) )
                         {    xright[j1].push_front( X[i][j1].back( ) );
                              X[i][j1].pop_back( ); 
                              xleft[j2].push_back( X[i][j2].front( ) );
                              X[i][j2].pop_front( );    }    }    }    }

          // Extra extra trimming step.

          for ( int j1 = 0; j1 < X[i].isize( ); j1++ )
          for ( int j2 = 0; j2 < X[i].isize( ); j2++ )
          {    if ( j1 == j2 ) continue;
               if ( X[i][j1].empty( ) || X[i][j2].empty( ) ) continue;
               vec<int> v1, v2;
               v1.push_back( to_left[ X[i][j1][0].second ] );
               for ( auto x : X[i][j1] ) v1.push_back( to_right[ x.second ] );
               v2.push_back( to_left[ X[i][j2][0].second ] );
               for ( auto x : X[i][j2] ) v2.push_back( to_right[ x.second ] );
               vec<int> s1(v1), s2(v2);
               UniqueSort(s1), UniqueSort(s2);
               vec<int> I = Intersection( s1, s2 );
               if ( !I.solo( ) ) continue;
               vec<int> p1, p2;
               for ( int l = 0; l < v1.isize( ); l++ )
                    if ( v1[l] == I[0] ) p1.push_back(l);
               for ( int l = 0; l < v2.isize( ); l++ )
                    if ( v2[l] == I[0] ) p2.push_back(l);
               if ( !p1.solo( ) || !p2.solo( ) ) continue;
               if ( p1.isize( ) - p1[0] > MAX_MOVE + 1 ) continue;
               if ( p2[0] > MAX_MOVE ) continue;
               if ( p1[0] == 0 || p2[0] == v2.isize( ) - 1 ) continue;
               int len = X[i][j1].back( ).first.Length( );
               int d = X[i][j1].back( ).second, dstart = X[i][j1].back( ).third;
               int skip1 = dlens[d] + K - 1 - dstart - len;
               int skip2 = X[i][j2].front( ).third;
               int trim1 = -skip1, trim2 = -skip2;
               for ( int l = p1[0]; l < v1.isize( ) - 1; l++ )
                    trim1 += dlens[ X[i][j1][l].second ];
               for ( int l = 0; l < p2[0]; l++ ) 
                    trim2 += dlens[ X[i][j2][l].second ];
               if ( trim1 > MAX_TRIM || trim2 > MAX_TRIM ) continue;
               for ( int l = p1[0]; l < v1.isize( ) - 1; l++ ) 
               {    xright[j1].push_front( X[i][j1].back( ) );
                    X[i][j1].pop_back( );    }
               for ( int l = 0; l < p2[0]; l++ ) 
               {    xleft[j2].push_back( X[i][j2].front( ) );
                    X[i][j2].pop_front( );    }    }

          UniqueSortSync( X[i], xleft, xright );

          // Define a function.

          auto Cat = [&]( const vec<int>& e )
          {    basevector x = tigs[ e[0] ];
               for ( int j = 1; j < e.isize( ); j++ )
                    x = TrimCat( K, x, tigs[ e[j] ] );
               return x;    };

          // Find connections.

          vec<quad<int,int,int,pair<int,int>>> connections;
          vec<malign> maligns;
          if (VERBOSE) out << "Finding pairs that might be bridged:\n";
          for ( int j1 = 0; j1 < X[i].isize( ); j1++ )
          for ( int j2 = 0; j2 < X[i].isize( ); j2++ )
          {    if ( j1 == j2 ) continue;
               int gstart1 = X[i][j1].front( ).first.Start( );
               int gstop1 = X[i][j1].back( ).first.Stop( );
               int gstart2 = X[i][j2].front( ).first.Start( );
               int gstop2 = X[i][j2].back( ).first.Stop( );
               int gap = gstart2 - gstop1;
               if ( gap > MAX_GAP || -gap > MAX_OVERLAP ) continue;
               vec<int> p1, p2;
               for ( auto m : X[i][j1] ) p1.push_back( m.second );
               for ( auto m : X[i][j2] ) p2.push_back( m.second );
               int d1 = p1.back( ), d2 = p2.front( );
               triple<ho_interval,int,int> z1 = X[i][j1].back( );
               triple<ho_interval,int,int> z2 = X[i][j2].front( );
               int gleft = z1.first.Stop( ), gright = z2.first.Start( );
               int aleft = z1.third + z1.first.Length( ), aright = z2.third;
               if (VERBOSE)
               {    out << "\n";
                    out << j1+1 << " ==> " << j2+1 << ": ";
                    PRINT3_TO( out, gstop1, gstart2, gap );
                    out << "p1 = ";
                    if ( p1.size( ) <= 3 ) out << printSeq(p1) << endl;
                    else out << p1.front( ) << ".." << p1.back( ) << endl;
                    out << "p2 = ";
                    if ( p2.size( ) <= 3 ) out << printSeq(p2) << endl;
                    else out << p2.front( ) << ".." << p2.back( ) << endl;
                    PRINT5_TO( out, p1.back( ), gleft, gright, aleft, aright );
                    PRINT3_TO( 
                         out, z1.third, z1.first.Start( ), z1.first.Stop( ) );    }
               basevector G, A;
               align al;

               // Aligner and other functions.

               auto Align = [&]( )
               {    if (VERBOSE)
                    {    DPRINT2_TO( out, G.size( ), A.size( ) );
                         out << "aligning " << G.ToString( ) << " to "
                              << A.ToString( ) << endl;    }
                    int score;
                    if ( G.size( ) == 0 && A.size( ) == 0 )
                    {    score = 0;
                         avector<int> gaps, lengths;
                         if (VERBOSE) out << "adding empty align" << endl;
                         al.Set( 0, 0, gaps, lengths );    }
                    else if ( G.size( ) == 0 && A.size( ) > 0 )
                    {    score = GAP_OPEN_PENALTY 
                              + ( A.size( ) - 1 ) * GAP_EXTEND_PENALTY;
                         avector<int> gaps(1), lengths(1);
                         gaps(0) = A.size( ), lengths(0) = 0;
                         al.Set( 0, 0, gaps, lengths );    }
                    else if ( A.size( ) == 0 && G.size( ) > 0 )
                    {    score = GAP_OPEN_PENALTY 
                              + ( G.size( ) - 1 ) * GAP_EXTEND_PENALTY;
                         avector<int> gaps(1), lengths(1);
                         gaps(0) = -G.size( ), lengths(0) = 0;
                         al.Set( 0, 0, gaps, lengths );    }
                    else
                    {    score = SmithWatAffineEndFix( G, A, al, true, true,
                              MISMATCH_PENALTY, GAP_OPEN_PENALTY, 
                              GAP_EXTEND_PENALTY );    }
                    double dscore = score/10.0;
                    if (VERBOSE) DPRINT_TO( out, dscore );
                    return score;    };
               // incomplete validation!!
               auto Validate = [&]( const aalign& aa )
               {    if ( aa.a.Nblocks( ) == 0 ) return;
                    if ( aa.p.empty( ) )
                    {    out << "Validate FAILS: empty p" << endl;
                         fails++;
                         return;    }
                    basevector b = Cat( D.O( aa.p[0] ) );
                    for ( int j = 1; j < aa.p.isize( ); j++ )
                         b = TrimCat( K, b, Cat( D.O( aa.p[j] ) ) );
                    basevector G;
                    G.SetToSubOf( genomef[g], aa.gpos.Start( ), aa.gpos.Length( ) );
                    int score = 0;
                    const align& a = aa.a;
                    if ( a.pos1( ) > 0 )
                    {    score += GAP_OPEN_PENALTY
                              + ( a.pos1( ) - 1 ) * GAP_EXTEND_PENALTY;    }
                    if ( a.Pos1( ) < G.isize( ) )
                    {    score += GAP_OPEN_PENALTY
                              + ( ( G.isize( ) - a.Pos1( ) ) - 1 ) 
                              * GAP_EXTEND_PENALTY;    }
                    int p1 = a.pos1( ), p2 = a.pos2( );
                    for ( int j = 0; j < a.Nblocks( ); j++ ) 
                    {    if ( a.Gaps(j) > 0 )
                         {    int penalty = GAP_OPEN_PENALTY
                                   + ( a.Gaps(j) - 1 ) * GAP_EXTEND_PENALTY;
                              score += penalty;
                              p2 += a.Gaps(j);    }
                         if ( a.Gaps(j) < 0 )
                         {    int penalty = GAP_OPEN_PENALTY
                                   + ( -a.Gaps(j) - 1 ) * GAP_EXTEND_PENALTY;
                              score += penalty;
                              p1 -= a.Gaps(j);    }
                         for ( int x = 0; x < a.Lengths(j); x++ ) 
                         {    if ( p1 < 0 || p1 >= G.isize( )
                                   || p2 < 0 || p2 >= A.isize( ) )
                              {    out << "\nValidate FAILS: invalid base "
                                        << "position" << endl;
                                   fails++;
                                   out << "p = " << printSeq(aa.p) << endl;
                                   PRINT3_TO( out, p2, G.size( ), A.size( ) );
                                   PRINT2_TO( 
                                        out, aa.gpos.Start( ), aa.gpos.Length( ) );
                                   PRINT2_TO( out, aa.a.pos1( ), aa.a.pos2( ) );
                                   for ( int j = 0; j < aa.a.Nblocks( ); j++ )
                                   {    out << "[" << j+1 << "] gap = " 
                                             << aa.a.Gaps(j) << ", length = " 
                                             << aa.a.Lengths(j) << endl;    }
                                   PRINT2_TO( out, score, aa.score );
                                   out << endl;
                                   return;    }
                              if ( G[p1] != A[p2] ) score += MISMATCH_PENALTY;
                              ++p1;
                              ++p2;    }    }
                    if ( score != aa.score )
                    {    out << "\nValidate FAILS: score does not agree" << endl;
                         fails++;
                         PRINT2_TO( out, G.size( ), A.size( ) );
                         PRINT2_TO( out, aa.gpos.Start( ), aa.gpos.Length( ) );
                         PRINT2_TO( out, aa.a.pos1( ), aa.a.pos2( ) );
                         for ( int j = 0; j < aa.a.Nblocks( ); j++ )
                         {    out << "[" << j+1 << "] gap = " << aa.a.Gaps(j)
                                   << ", length = " << aa.a.Lengths(j) << endl;    }
                         PRINT2_TO( out, score, aa.score );
                         out << endl;
                         return;    }    };

               // Idea should be to construct a digraph whose vertices are the
               // perfect stretches, and whose edges are connections, with penalties
               // computed via what's done below.  Then find all paths through
               // the graph from a source to a sink.  Each has start/stop on
               // genome and an error penalty.

               if ( p1.back( ) == p2.front( ) ) 
               {    if (VERBOSE) out << "shared edge\n";

                    // Define the exact things to be aligned.

                    // on g = [ z1.first.Stop( ), z2.second.Start( ) )
                    // path = { p1.back( ) } = { p2.front( ) }
                    // on path = [ z1.third + z1.first.Length( ), z2.third )
                    // but both may have negative length!

                    // Deal with negative intervals.

                    if ( gleft > gright )
                    {    if ( aleft - (gleft-gright) < 0 )
                         {    out << "Oops, need different method." << endl;
                              continue;    }
                         aleft -= ( gleft-gright ), gleft = gright;    }
                    if ( aleft > aright )
                    {    if ( gleft - (aleft-aright) < 0 )
                         {    out << "Gotta do something different." << endl;
                              continue;    }
                         gleft -= (aleft-aright), aleft = aright;    }

                    // Now g:[gleft,gright) may be aligned to
                    // d: [aleft,aright)
                    // where d = p1.back( ).

                    // Seems like we're aligning K-1 bases too many.
                    // Actually K-1 should be added, so there's a bug somewhere.

                    if (VERBOSE) PRINT4_TO( out, gleft, gright, aleft, aright );
                    G.SetToSubOf( genomef[g], gleft, gright - gleft + K - 1 );
                    int d = p1.back( );
                    A = Cat( D.O(d) );
                    A.SetToSubOf( A, aleft, aright - aleft + K - 1 );
                    int score = Align( );
                    aalign aa;
                    aa.p = {d};
                    aa.gpos = ho_interval( gleft, gright + K - 1 );
                    aa.a = al;
                    aa.a.AddToPos2(aleft);
                    aa.score = score;
                    A = Cat( D.O(d) );
                    Validate(aa); // YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY
                    malign ma;
                    ma.aligns = {aa};
                    connections.push( j1, j2, score, make_pair(0,0) );
                    maligns.push_back(ma);    }
               else
               {    int v = to_right[ p1.back( ) ], w = to_left[ p2.front( ) ];
                    vec<vec<int>> paths;
                    if ( v == w ) 
                    {    paths.push_back( { } );
                         if (VERBOSE) out << "WOW, SAME VERTEX\n";    }
                    else
                    {    if (VERBOSE)
                         {    out << "FINDING PATHS FROM " << p1.back( )
                                   << " TO " << p2.front( ) << endl;    }
                         Bool OK = M.AllPathsLengthRange( v, w, 0, MAX_LENGTH, 
                              to_left, to_right, paths, MAX_PATHS, MAX_ITERATIONS,
                              False, MAX_COPIES );
                         if ( !OK ) 
                         {    if (VERBOSE) out << "FAILED!\n";
                              continue;    }
                         else if ( paths.empty( ) ) 
                         {    if (VERBOSE) out << "NO PATH!\n";
                              continue;    }    }
                    if (VERBOSE) 
                         out << "Found " << paths.size( ) << " paths." << endl;
                    int min_score = 1000000000;
                    malign ma;
                    Bool have_gapfree = False;
                    for ( int l = 0; l < paths.isize( ); l++ )
                    {    const vec<int>& p = paths[l];
                         vec<int> gaps;
                         for ( int i = 0; i < p.isize( ); i++ )
                         {    int d = p[i];
                              if ( D.O(d)[0] < 0 ) gaps.push_back(i);    }
                         if ( gaps.empty( ) ) have_gapfree = True;    }
                    for ( int l = 0; l < paths.isize( ); l++ )
                    {    const vec<int>& p = paths[l];

                         // Handle gaps.  Note that if there is more than one gap,
                         // we treat the entire section (from the first gap to the
                         // last one) as a single gap.

                         vec<int> gaps;
                         for ( int i = 0; i < p.isize( ); i++ )
                         {    int d = p[i];
                              if ( D.O(d)[0] < 0 ) gaps.push_back(i);    }
                         if ( gaps.nonempty( ) )
                         {    if (have_gapfree) continue;
                              if (VERBOSE)
                              {    out << "Punting path with " << gaps.size( ) 
                                        << " gaps.\n";    }
                              vec<int> left, right;
                              for ( int m = 0; m < gaps.front( ); m++ )
                                   left.push_back( p[m] );
                              for ( int m = gaps.back( ) + 1; m < p.isize( ); m++ )
                                   right.push_back( p[m] );
                              vec<int> gpath;
                              for ( int m = gaps.front( ); m <= gaps.back( ); m++ )
                                   gpath.push_back( p[m] );
                              malign mac;
                              int score = BuildGapConnection( mac, d1, d2, tigs, 
                                   left, right,  aleft, aright, g, gstop1, gright, 
                                   gpath, hb, D, genomef, MISMATCH_PENALTY, 
                                   GAP_OPEN_PENALTY, GAP_EXTEND_PENALTY, GAP_SLACK, 
                                   MIN_ADD, VERBOSE, out );
                              if ( score < min_score )
                              {    min_score = score;
                                   ma = mac;    }
                              continue;    }
                         // Huh, how could paths be empty??
                         if ( paths.nonempty( ) && !have_gapfree )
                         {    connections.push( j1, j2, 
                                   ma.aligns[0].score + ma.aligns[1].score,
                                   make_pair(1,0) );
                              maligns.push_back(ma);
                              continue;    }

                         // Now do the case of a gapfree path.

                         if (VERBOSE) out << "ALIGNING PATH" << endl;

                         // Deal with negative intervals.

                         if ( gleft > gright )
                         {    if ( aleft - (gleft-gright) < 0 )
                              {    out << "Oops, need different method - 3." 
                                        << endl;
                                   continue;    }
                              aleft -= ( gleft-gright ), gleft = gright;    }

                         // Make sure enough room.

                         if ( dlens[d1] + K - 1 - aleft < K - 1 )
                         {    int delta = aleft - dlens[d1];
                              if ( delta > aleft )
                              {    out << "This shouldn't have happened."
                                        << " (2)" << endl;
                                   continue;    }
                              aleft -= delta, gleft -= delta;    }

                         // Define the exact things to be aligned.

                         // on g = [ z1.first.Stop( ), z2.second.Start( ) )
                         // path = p1.back( ), p, p2.front( )
                         // on path from z1.third + z1.first.Length( ) 
                         // to z2.third

                         // Align.

                         int dstop1 = aleft, dstart2 = aright;
                         G.SetToSubOf( genomef[g], gleft, gright - gleft );
                         basevector A1 = Cat( D.O(d1) ), A2 = Cat( D.O(d2) );
                         A1.SetToSubOf( A1, dstop1, A1.isize( ) - dstop1 );
                         A2.SetToSubOf( A2, 0, dstart2 );
                         A = A1;
                         for ( auto d : p ) A = TrimCat( K, A, Cat( D.O(d) ) );
                         A = TrimCat( K, A, A2 );
                         int score = Align( );    
                         if ( score < min_score )
                         {    min_score = score;
                              aalign aa;
                              aa.p = {d1};
                              aa.p.append(p);
                              aa.p.push_back( {d2} );
                              aa.gpos = ho_interval( gleft, gright );
                              aa.a = al;
                              aa.a.AddToPos2(dstop1);
                              aa.score = score;

                              // Testing for score > 0 below is a workaround for a
                              // bug in assemblies which allows a seq gap edge
                              // to have length K-1 bases.  These are cases where
                              // the seq gap edges should not be there at all, and
                              // rather the two abutting edges should be directly
                              // joined.

                              if ( score > 0 ) ma.aligns = {aa};

                              A1 = Cat( D.O(d1) );
                              A = A1;
                              for ( auto d : p ) A = TrimCat( K, A, Cat( D.O(d) ) );
                              A = TrimCat( K, A, A2 );
                              Validate(aa); // YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY
                              ma.gaps.clear( );    }    }
                    if ( paths.nonempty( ) && !have_gapfree )
                    {    connections.push( j1, j2, 
                              ma.aligns[0].score + ma.aligns[1].score, 
                              make_pair(1,0) );
                         maligns.push_back(ma);    }
                    else if ( min_score < 1000000000 )
                    {    connections.push( j1, j2, min_score, make_pair(0,0) );
                         maligns.push_back(ma);    }    }    }

          // Announce.

          out << "\nSUMMARY FOR " << g << endl << endl;
          for ( int j = 0; j < X[i].isize( ); j++ )
          {    vec<int> p;
               for ( auto m : X[i][j] ) p.push_back( m.second );
               int gstart = X[i][j].front( ).first.Start( );
               int gstop = X[i][j].back( ).first.Stop( );
               out << "[" << j+1 << "] " << gstart << "-" << gstop;
               if ( gstop == genomef[g].isize( ) ) out << " = end";
               out << ": ";
               if ( p.size( ) <= 3 ) out << printSeq(p) << endl;
               else out << p.front( ) << ".." << p.back( ) << endl;    }

          // Find likely uncaptured gaps and extend the alignment of edges that 
          // flank them.

          for ( int v1 = 0; v1 < X[i].isize( ); v1++ )
          for ( int v2 = 0; v2 < X[i].isize( ); v2++ )
          {    if ( v1 == v2 ) continue;
               Bool owned = False;
               for ( auto m : connections )
               {    if ( m.first == v1 && m.second == v2 )
                    {    owned = True;
                         break;    }    }
               if (owned) continue;
               int start = X[i][v1].back( ).first.Stop( );
               int stop = X[i][v2].front( ).first.Start( );
               if ( stop - start < MIN_UNCAP || stop - start > MAX_UNCAP ) continue;
               triple<ho_interval,int,int> z1 = X[i][v1].back( );
               triple<ho_interval,int,int> z2 = X[i][v2].front( );
               int d1 = z1.second, d2 = z2.second;
               int aleft = z1.third + z1.first.Length( ), aright = z2.third;
               vec<int> left, right, gpath;
               int gleft = z1.first.Stop( ), gright = z2.first.Start( );
               if ( gleft > gright )
               {    if ( aleft - (gleft-gright) < 0 ) continue;
                    aleft -= ( gleft-gright ), gleft = gright;    }
               if ( dlens[d1] + K - 1 - aleft > MAX_UNCAP || aright > MAX_UNCAP ) 
                    continue;
               if (VERBOSE)
               {    out << endl << Date( ) << ": building gap connection from [" 
                         << v1+1 << "] to [" << v2+1 << "]" << endl;
                    PRINT2_TO( out, start, stop );    }
               malign ma;
               int score = BuildGapConnection( ma, d1, d2, tigs, left, right, aleft, 
                    aright, g, gleft, gright, gpath, hb, D, genomef,
                    MISMATCH_PENALTY, GAP_OPEN_PENALTY, GAP_EXTEND_PENALTY, 
                    GAP_SLACK, MIN_ADD, VERBOSE, out );
               connections.push( v1, v2, score, make_pair(0,1) );
               maligns.push(ma);    }

          // Put back some trims.

          vec<Bool> ext_left( X[i].size( ), False );
          vec<Bool> ext_right( X[i].size( ), False );
          for ( int i = 0; i < connections.isize( ); i++ )
          {    const quad<int,int,int,pair<int,int>>& m = connections[i];
               int v = m.first, w = m.second;
               ext_right[v] = True, ext_left[w] = True;    }
          for ( int j = 0; j < X[i].isize( ); j++ )
          {    if ( !ext_left[j] ) X[i][j].prepend( xleft[j] );
               if ( !ext_right[j] ) X[i][j].append( xright[j] );    }

          // Start to build graph.

          int N = X[i].size( ) + 2;
          vec<vec<int>> from(N), to(N), from_edge_obj(N), to_edge_obj(N);
          vec<int> edges;
          vec<triple<int,int,int>> edgesf;
          vec<malign> edgesm;

          // Add regular edges.

          for ( int j = 0; j < connections.isize( ); j++ )
          {    const quad<int,int,int,pair<int,int>>& m = connections[j];
               int v = m.first, w = m.second;
               from[v].push_back(w), to[w].push_back(v);
               from_edge_obj[v].push_back( edges.size( ) );
               to_edge_obj[w].push_back( edges.size( ) );
               int penalty = m.third + m.fourth.first * CAP_PENALTY;
               if ( m.fourth.second > 0 )
               {    penalty += UNCAP_PENALTY;
                    int gstart = X[i][v].back( ).first.Stop( );
                    int gstop = X[i][w].front( ).first.Start( );
                    int gap = Abs( gstop - gstart );
                    if ( gap > 0 )
                         penalty += GAP_OPEN_PENALTY + (gap-1) + GAP_EXTEND_PENALTY;
                              }
               edges.push(penalty);
               edgesf.push( m.third, m.fourth.first, m.fourth.second );
               edgesm.push( maligns[j] );    }

          // Add end edges.

          int vstart = X[i].size( ), vstop = X[i].size( ) + 1;
          for ( int v = 0; v < X[i].isize( ); v++ )
          {    from[vstart].push_back(v), to[v].push_back(vstart);
               from_edge_obj[vstart].push_back( edges.size( ) );
               to_edge_obj[v].push_back( edges.size( ) );
               int gstart = X[i][v].front( ).first.Start( );
               int gstop = X[i][v].back( ).first.Stop( );
               edges.push_back(END_PENALTY*gstart);
               edgesf.push( END_PENALTY*gstart, 0, 0 );
               edgesm.push( malign( ) );
               from[v].push_back(vstop), to[vstop].push_back(v);
               from_edge_obj[v].push_back( edges.size( ) );
               to_edge_obj[vstop].push_back( edges.size( ) );
               int dist_to_end = genomef[g].isize( ) - gstop;
               edges.push_back(END_PENALTY*dist_to_end);
               edgesf.push( END_PENALTY*dist_to_end, 0, 0 );
               edgesm.push( malign( ) );    }

          // Add new uncaptured gap edges.

          for ( int v1 = 0; v1 < X[i].isize( ); v1++ )
          for ( int v2 = 0; v2 < X[i].isize( ); v2++ )
          {    if ( v1 == v2 ) continue;
               Bool owned = False;
               for ( auto m : connections )
               {    if ( m.first == v1 && m.second == v2 )
                    {    owned = True;
                         break;    }    }
               if (owned) continue;
               from[v1].push_back(v2), to[v2].push_back(v1);
               from_edge_obj[v1].push_back( edges.size( ) );
               to_edge_obj[v2].push_back( edges.size( ) );
               int gstart = X[i][v1].back( ).first.Stop( );
               int gstop = X[i][v2].front( ).first.Start( );
               int gap = Abs( gstop - gstart ), penalty = 0;
               if ( gap > 0 )
                    penalty = GAP_OPEN_PENALTY + (gap-1) + GAP_EXTEND_PENALTY;
               edges.push_back( penalty + UNCAP_PENALTY );
               edgesf.push( 0, 0, 1 );
               edgesm.push( malign( ) );    }

          // Now build the graph and find the shortest path through it.

          for ( int v = 0; v < N; v++ )
          {    SortSync( from[v], from_edge_obj[v] );
               SortSync( to[v], to_edge_obj[v] );    }
          digraphE<int> G( from, to, edges, to_edge_obj, from_edge_obj );
          vec<int> p;
          G.ShortestPath( vstart, vstop, p );

          // Build alignment stack.

          vec<malign> M;
          for ( int l = 1; l <= p.isize( ) - 2; l++ )
          {    aalign aa;
               const vec<triple<ho_interval,int,int>>& P = X[i][p[l]];
               for ( auto m : X[i][p[l]] ) aa.p.push_back( m.second );
               int gstart = P.front( ).first.Start( );
               int gstop = P.back( ).first.Stop( );
               aa.gpos = ho_interval( gstart, gstop );
               align a;
               avector<int> gaps(1), lengths(1);
               gaps(0) = 0, lengths(0) = gstop - gstart;
               a.Set( 0, gstop - gstart, gaps, lengths );
               a.Setpos2( P.front( ).third );
               aa.a = a;
               aa.score = 0;
               malign ma;
               ma.aligns = {aa};
               M.push(ma);
               if ( l < p.isize( ) - 2 )
               {    int e = -1, v = p[l], w = p[l+1], min_score = 1000000000;
                    for ( int j = 0; j < G.From(v).isize( ); j++ )
                    {    if ( G.From(v)[j] != w ) continue;
                         int x = G.IFrom( v, j );
                         int score = G.O(x);
                         if ( score < min_score )
                         {    min_score = score;
                              e = x;    }    }
                    M.push_back( edgesm[e] );    }    }

          // Remove overlaps.

          remove_overlaps:
          for ( int i = 1; i < M.isize( ); i++ )
          {    if ( M[i-1].aligns.empty( ) || M[i].aligns.empty( ) ) continue;
               aalign &aa1 = M[i-1].aligns.back( ), &aa2 = M[i].aligns.front( );
               if ( aa1.a.Nblocks( ) == 0 || aa2.a.Nblocks( ) == 0 ) continue;
               int overlap = aa1.gpos.Stop( ) - aa2.gpos.Start( );
               if ( overlap > 0 )
               {    if ( overlap > aa1.a.Lengths( aa1.a.Nblocks( ) - 1 ) )
                    {    if ( overlap <= aa2.a.Lengths(0) )
                         {    aa2.gpos.AddToStart(overlap);
                              aa2.a.AddToLength( 0, -overlap );
                              aa2.a.AddToPos2(overlap);
                              continue;    }
                         else 
                         {    out << "overlap problem" << endl;
                              PRINT6_TO( out, g, aa1.gpos.Stop( ), aa2.gpos.Start( ),
                                   overlap, aa1.a.Lengths( aa1.a.Nblocks( ) - 1 ),
                                   aa2.a.Lengths(0) );
                              problems++;    
                              // GRUESOME WORKAROUND:
                              M.erase( M.begin( ) + i-1 );
                              goto remove_overlaps;    }    }
                    aa1.gpos.AddToStop(-overlap);
                    aa1.a.AddToLength( aa1.a.Nblocks( ) - 1, -overlap );    }    }

          // Announce and merge to form one malign.

          if (VERBOSE) out << "\nALIGNMENT STACK\n";
          int last_gstop = -1, last_a = -1, last_astop = -1;
          malign mal;
          aalign cur;
          Bool started = False, broken = False;
          for ( int j = 0; j < M.isize( ); j++ )
          {    int count = 0;
               if ( M[j].aligns.empty( ) ) 
               {    if (VERBOSE) out << "[" << j+1 << ".1" << "] uncaptured gap\n";
                    last_gstop = -1;    
                    mal.aligns.push_back(cur);
                    mal.gaps.push_back( { } );
                    started = False;    }
               for ( int k = 0; k < M[j].aligns.isize( ); k++ )
               {    const aalign& aa = M[j].aligns[k];
                    if ( aa.a.Nblocks( ) > 0 )
                    {    if ( last_gstop < 0 ) 
                         {    cur = aa;
                              started = True;    }
                         else
                         {    if ( aa.gpos.Start( ) != last_gstop )
                              {    out << "\nANOMALY 1\n";
                                   PRINT4_TO( out, 
                                        j, k, aa.gpos.Start( ), last_gstop );
                                   out << endl;
                                   anoms++;
                                   broken = True;
                                   break;    }
                              Bool aok = False;
                              int a = aa.p.front( ), apos = aa.a.pos2( );
                              if ( a == last_a && apos == last_astop ) aok = True;
                              else
                              {    int v = to_right[last_a];
                                   for ( int m = 0; m < D.From(v).isize( ); m++ )
                                   {    if ( D.IFrom(v)[m] == a )
                                        {    if ( last_astop - dlens[last_a]
                                                  == apos )
                                             {    aok = True;    }    }    }    }
                              if ( !aok ) 
                              {    out << "\nANOMALY 2\n\n";
                                   PRINT7_TO( out, j, k, aa.gpos.Start( ), a, last_a,
                                        apos, last_astop );
                                   anoms++;
                                   broken = True;
                                   break;    }
                              if ( !started ) 
                              {    cur = aa;
                                   started = True;    }
                              else if ( a == last_a && apos == last_astop )
                              {    cur.gpos.SetStop( aa.gpos.Stop( ) );
                                   for ( int l = 1; l < aa.p.isize( ); l++ )
                                        cur.p.push_back( aa.p[l] );
                                   int nb1 = cur.a.Nblocks( ), nb2 = aa.a.Nblocks( );
                                   cur.a.SetNblocks( nb1 + nb2 );
                                   for ( int l = 0; l < nb2; l++ )
                                   {    cur.a.SetGap( nb1 + l, aa.a.Gaps(l) );
                                        cur.a.SetLength( nb1 + l, 
                                             aa.a.Lengths(l) );    }
                                   cur.score += aa.score;    }
                              else
                              {    cur.gpos.SetStop( aa.gpos.Stop( ) );
                                   cur.p.append( aa.p );
                                   int nb1 = cur.a.Nblocks( ), nb2 = aa.a.Nblocks( );
                                   cur.a.SetNblocks( nb1 + nb2 );
                                   for ( int l = 0; l < nb2; l++ )
                                   {    cur.a.SetGap( nb1 + l, aa.a.Gaps(l) );
                                        cur.a.SetLength( nb1 + l, 
                                             aa.a.Lengths(l) );    }
                                   cur.score += aa.score;    }    }
                         last_gstop = aa.gpos.Stop( );    
                         int apos = aa.a.pos2( ) + aa.a.extent2( );
                         const vec<int>& p = aa.p;
                         for ( int m = 0; m < p.isize( ) - 1; m++ )
                              apos -= dlens[ p[m] ];
                         last_a = p.back( );
                         last_astop = apos;    }
                    if (broken) break;
                    if (VERBOSE)
                    {    if ( aa.a.Nblocks( ) > 0 )
                         {    out << "[" << j+1 << "." << ++count << "] ";
                              out << aa.a.pos2( ) << ":(" << printSeq(aa.p) 
                                   << ") --> [" << aa.gpos.Start( ) << "," 
                                   << aa.gpos.Stop( ) << ")" << ", score = " 
                                   << aa.score/10.0 << endl;    }    }
                    if ( k < M[j].aligns.isize( ) - 1 )
                    {    if (VERBOSE)
                         {    out << "[" << j+1 << "." << ++count << "] gap = " 
                                   << printSeq( M[j].gaps[k] ) << endl;    }
                         last_gstop = -1;    
                         mal.aligns.push_back(cur);
                         mal.gaps.push_back( M[j].gaps[k] );    
                         started = False;    }    }    }
          stats[i] = make_quad( 0, 0, 0, 0 );
          covstats[i] = make_pair( 0, genomef[g].size( ) );
          #pragma omp critical
          {    total_total += genomef[g].size( );    }
          if ( M.empty( ) )
          {    out << "Nothing aligned." << endl;
               fullreports[i] = out.str( );
               continue;    }
          mal.aligns.push_back(cur);
          if (broken) 
          {    fullreports[i] = out.str( );
               continue;    }

          // Compute coverage.

          vec<ho_interval> cov;
          for ( int l = 0; l < mal.aligns.isize( ); l++ )
          {    int start = mal.aligns[l].gpos.Start( );
               const align& a = mal.aligns[l].a;
               int p1 = a.pos1( );
               for ( int j = 0; j < a.Nblocks( ); j++ ) 
               {    if ( a.Gaps(j) < 0 ) p1 -= a.Gaps(j);
                    if ( a.Lengths(j) < 0 ) continue; // ZZZZZZZZZZZZZZZZZZZZZZZZZZZ
                    cov.push( start + p1, start + p1 + a.Lengths(j) );
                    p1 += a.Lengths(j);    }    }

          // Compute number of edges and lines, save paths.

          vec<int> P;
          for ( int l = 0; l < mal.aligns.isize( ); l++ ) 
          {    P.append( mal.aligns[l].p );
               allpaths_flat[i].push_back( mal.aligns[l].p );    }
          int nedges = P.size( ), nlines = 1;
          vec<int> ls;
          for ( auto d : P ) if ( d < D.E( ) && tol[d] >= 0 ) ls.push_back( tol[d] );
          for ( int l = 1; l < ls.isize( ); l++ ) if ( ls[l] != ls[l-1] ) nlines++;

          // Keep mal.

          allmaligns[i] = mal;

          // Create allpaths.

          for ( int l = 0; l < mal.aligns.isize( ); l++ ) 
          {    allpaths[i].append( mal.aligns[l].p );
               if ( l < mal.aligns.isize( ) - 1 ) 
                    allpaths[i].append( mal.gaps[l] );    }

          // Announce path.  This incorrectly computes the number of uncaptured
          // gaps, not sure why, so computed a second time.

          out << "\nPATH " << "[" << p[1] + 1 << "] --> [" 
               << p[ p.isize( ) - 2 ] + 1 << "]:" << endl;
          {    int score = 0, caps = 0, uncaps = 0, endgaps = 0;
               for ( int l = 1; l < p.isize( ) - 2; l++ )
               {    out << "[" << p[l] + 1 << "]";
                    int e = -1, v = p[l], w = p[l+1], min_score = 1000000000;
                    for ( int j = 0; j < G.From(v).isize( ); j++ )
                    {    if ( G.From(v)[j] != w ) continue;
                         int x = G.IFrom( v, j );
                         int score = G.O(x);
                         if ( score < min_score )
                         {    min_score = score;
                              e = x;    }    }
                    int s = edgesf[e].first, c = edgesf[e].second;
                    int u = edgesf[e].third;

                    // Recompute score.

                    const malign& m = edgesm[e];
                    s = 0;
                    for ( auto& aa : m.aligns )
                    {    const align& a = aa.a;
                         basevector G( 
                              genomef[g], aa.gpos.Start( ), aa.gpos.Length( ) );
                         basevector A = Cat( D.O( aa.p[0] ) );
                         for ( int j = 1; j < aa.p.isize( ); j++ )
                              A = TrimCat( K, A, Cat( D.O( aa.p[j] ) ) );
                         int p1 = a.pos1( ), p2 = a.pos2( );
                         for ( int j = 0; j < a.Nblocks( ); j++ ) 
                         {    if ( a.Gaps(j) != 0 )
                              {    int g = Abs( a.Gaps(j) );
                                   s += GAP_OPEN_PENALTY
                                        + GAP_EXTEND_PENALTY_SHOW * (g-1);    }
                              if ( a.Gaps(j) > 0 ) p2 += a.Gaps(j);
                              if ( a.Gaps(j) < 0 ) p1 -= a.Gaps(j);
                              for ( int x = 0; x < a.Lengths(j); x++ )
                              {    if ( G[p1] != A[p2] ) s += MISMATCH_PENALTY;
                                   ++p1;
                                   ++p2;    }    }    }

                    // Print and save.

                    out << " --(" << s/10.0 << "," << c << "," << u 
                         << ";" << edges[e] << ")--> ";
                    score += s, caps += c, uncaps += u;    }
               if ( X[i][p[1]].front( ).first.Start( ) > 0 ) endgaps++;
               if ( X[i][ p[ p.isize( ) - 2 ] ].back( ).first.Stop( ) 
                    < genomef[g].isize( ) )
               {    endgaps++;    }
               out << "[" << p[ p.isize( ) - 2 ] + 1 << "]" << endl;

               // Recompute uncaptured gaps, see comment above.

               uncaps = 0;
               for ( int l = 0; l < mal.aligns.isize( ); l++ )
               {    if ( l < mal.gaps.isize( ) )
                         if ( mal.gaps[l].empty( ) ) uncaps++;    }

               out << "\nSTATS FOR " << g << ": edges = " << nedges << ", lines = "
                    << nlines << ", errs = " << score/10.0 
                    << ", caps = " << caps << ", uncaps = " << uncaps 
                    << ", endgaps = " << endgaps << ", cov = " 
                    << PERCENT_RATIO( 4, TotalCovered(cov), genomef[g].isize( ) )
                    << endl;
               stats[i] = make_quad( score, caps, uncaps, endgaps );
               covstats[i] = make_pair( TotalCovered(cov), genomef[g].size( ) );
               counts[i] = make_pair( nedges, nlines );
               #pragma omp critical
               {    total_score += score, total_caps += caps;
                    total_uncaps += uncaps, total_endgaps += endgaps;
                    total_cov += TotalCovered(cov);
                    total_edges += nedges, total_lines += nlines;    }    }
          

          // Print alignment.

          ostringstream yout;
          yout << "\n===================================================="
               << "============================\n";
          yout << "\nALIGNMENT TO FINISHED SEQUENCE " << g << "\n";
          ostringstream xout;
          xout << "{" << setiosflags(ios::fixed) << setprecision(2)
               << round( ( 10000.0 * covstats[i].first ) 
                    / covstats[i].second ) / 100.0 << resetiosflags(ios::fixed) 
               << "," << stats[i].first/10.0 << "," << stats[i].second << ","
               << stats[i].third << "," << stats[i].fourth << "}";
          out << "stats = " << xout.str( ) << endl;
          
          const int MINPERF = 400;
          for ( int l = 0; l < mal.aligns.isize( ); l++ )
          {    aalign& aa = mal.aligns[l];
               align& a = aa.a;
               basevector G( genomef[g], aa.gpos.Start( ), aa.gpos.Length( ) );
               basevector A = Cat( D.O( aa.p[0] ) );
               for ( int j = 1; j < aa.p.isize( ); j++ )
                    A = TrimCat( K, A, Cat( D.O( aa.p[j] ) ) );
               a.Flip( );
               yout << "\n[" << aa.gpos.Start( ) << "," << aa.gpos.Stop( ) << "], ";
               yout << "path = ";
               {    if ( aa.p.size( ) <= 3 ) yout << printSeq(aa.p) << endl;
                    else 
                    {    yout << aa.p.front( ) << ".." << aa.p.back( ) 
                              << " [n=" << aa.p.size( ) << "]" << endl;    }    }
               if ( a.Pos1( ) > A.isize( ) || a.Pos2( ) > G.isize( ) )
               {    
                    #pragma omp critical
                    {    cout << "Warning: alignment overflow on clone " << g 
                              << "." << endl;    }    }
               PrintVisualAlignmentClean( True, yout, A, G, a );
               
               // Compute perfect match blocks
               // and the contribution to the weighted error rate
               // only if there is a sufficiently long match

               if ( PS && a.MaxPerf(A,G) >= MINPERF ) {
                    
                    vec<ho_interval> pi;
                    a.PerfectIntervals2( A, G, pi );

                    for ( int k = 0; k < pi.isize( ); k++ )
                    {    ho_interval h = pi[k];
                         int s = perfs[i].empty( ) ? 0 : perfs[i].back( ).Stop( );
                         if ( h.Start( ) > s )
                         {    for ( int m = s; m < h.Start( ); m++ )
                                   perfs[i].push( m, m+1 );    }
                         int start = Max( s, h.Start( ) ), stop = h.Stop( );
                         if ( stop > start ) 
                              perfs[i].push( start, stop );    }
                             
                    werr[i].first  += mal.aligns[l].score;
                    werr[i].second += a.extent2();
               }

               if ( l < mal.gaps.isize( ) )
               {    int start = mal.aligns[l].gpos.Stop( );
                    int stop = mal.aligns[l+1].gpos.Start( );
                    int len = stop - start;
                    if ( mal.gaps[l].empty( ) )
                    {    yout << "UNCAPTURED gap from " << start << " to "
                              << stop << ", len = " << len << endl;    }
                    else 
                    {    yout << "captured gap from " << start << " to "
                              << stop << ", len = " << len 
                              << ", path = {" << printSeq( mal.gaps[l] ) << "}"
                              << endl;    }    }    }
          reports[i] = yout.str( );
          yout << "===================================================="
               << "============================\n";
          out << yout.str( );
          fullreports[i] = out.str( );  }
     
     // Compute the N50 perfect stretch and weighted error rate
          
     if ( PS )
          {    cout << Date( ) << ": computing perfect stretch" << endl;
               // Compute weighted error rate.

               int64_t num = 0, den = 0;
               for ( int g = 0; g < gs.isize( ); g++ )
               {    num += werr[g].first;
                    den += werr[g].second;    }
               double err_rate = 0.0;
               if (den > 0)
                    err_rate = (num*100.0) / den;
               
               /*
               for ( int i = 0; i < gs.isize(); i++ ) {
                    cout << "CLONE " << gs[i] << endl;
                    vec<int> len;
                    for ( auto & x : perfs[i] )
                         len.push_back( x.Length() );
                    Sort(len);
                    PRINT5( perfs[i].size(), Min(len), Max(len), Mean(len), Median(len) );
               }*/
               for ( int i = 0; i < gs.isize(); i++ ) {
                    const int g = gs[i];
                    if ( perfs[g].nonempty( ) )
                    {    while( perfs[g].back( ).Stop( ) < genomef[g].isize( ) )
                         {    perfs[g].push( perfs[g].back( ).Stop( ),
                                   perfs[g].back( ).Stop( ) + 1 );    }    }
                    else
                    {    for ( int s = 0; s < genomef[g].isize( ); s++ )
                              perfs[g].push( s, s+1 );    }
               } 
               // Compute N50, allowing for random changing of the Fosmids.
          
               const int npasses = 1000;
               vec<int> N50s(npasses);
               vec<vec<int>> R(npasses);
               for ( int pass = 0; pass < npasses; pass++ )
               {    R[pass] = vec<int>( gs.size( ), vec<int>::IDENTITY );
                    Shuffle( R[pass].begin( ), R[pass].end( ) );    }
               cout << Date( ) << ": computing N50" << endl;
               #pragma omp parallel for
               for ( int pass = 0; pass < npasses; pass++ )
               {    vec<int> q;
                    for ( int gg = 0; gg < (int) gs.size( ); gg++ )
                    {    int g = R[pass][gg];
                         if ( gg > 0 )
                         {    q.back( ) += perfs[g].front( ).Length( );
                              for ( int j = 1; j < perfs[g].isize( ); j++ )
                                   q.push_back( perfs[g][j].Length( ) );    }
                         else
                         {    for ( int j = 0; j < perfs[g].isize( ); j++ )
                                   q.push_back( perfs[g][j].Length( ) );    }    }
                    N50s[pass] = N50(q);    }
               int N50PS = int(round( Mean(N50s)));
               cout << Date( ) << ": new evaluation ";
               PRINT(N50PS);
               cout << Date( ) << ": new evaluation ";
               PRINT(err_rate);
               //StatLogger::log("perf_N50", N50PS, "N50 perfect stretch size" );
               //StatLogger::log("weighted_error_rate", err_rate, "weighted error rate" );
          }


     // Print stats summary.

     ostringstream qout;
     qout << "\nSUMMARY\n" << "edges: " << total_edges << ", lines: " 
          << total_lines << ", cov: " << setiosflags(ios::fixed) 
          << setprecision(2) << round( ( 10000.0 * total_cov ) 
               / total_total ) / 100.0 << resetiosflags(ios::fixed) 
          << ", errs: " << setiosflags(ios::fixed) 
          << setprecision(1) << total_score/10.0 << resetiosflags(ios::fixed) 
          << ", caps: " << total_caps << ", uncaps: " << total_uncaps 
          << ", endgaps: " << total_endgaps << endl;
     String sout_str = qout.str( );
     qout << endl;
     ostringstream yout;
     yout << "{" << total_edges << "," << total_lines << ","
          << setiosflags(ios::fixed) 
          << setprecision(2) << round( ( 10000.0 * total_cov ) 
               / total_total ) / 100.0 << resetiosflags(ios::fixed) 
          << "," << setiosflags(ios::fixed) 
          << setprecision(1) << total_score/10.0 << resetiosflags(ios::fixed) 
          << "," << total_caps << "," << total_uncaps << "," << total_endgaps << "}";
     if ( anoms + fails + problems > 0 ) 
     {    qout << "SOMETHING FUNNY HAPPENED." << endl;
          PRINT3_TO( qout, anoms, fails, problems );    }
     cout << qout.str( );

     // Open files.

     if ( !write_files ) return yout.str( );
     Ofstream( wout, INDIR + "/o.report." + suffix );
     for ( int i = 0; i < gs.isize( ); i++ ) wout << reports[i];
     Ofstream( zout, INDIR + "/o.fullreport." + suffix );
     for ( int i = 0; i < gs.isize( ); i++ ) zout << fullreports[i];
     Ofstream( pout, INDIR + "/o.finpaths." + suffix );
     for ( int i = 0; i < gs.isize( ); i++ )
          pout << gs[i] << ": " << printSeq( allpaths[i] ) << endl;
     Ofstream( sout, INDIR + "/o.summary." + suffix );
     Ofstream( fout, INDIR + "/o.finseqs." + suffix );
     Ofstream( vout, INDIR + "/o.findata." + suffix );
     for ( int i = 0; i < gs.isize( ); i++ )
     {    const malign& mal = allmaligns[i];
          for ( int l = 0; l < mal.aligns.isize( ); l++ ) 
          {    const aalign& a = mal.aligns[l];
               // output = ( clone, clone_start, clone_stop, path_start, path_stop,
               //            assembly path )
               vout << i << " " << a.gpos.Start( ) << " " << a.gpos.Stop( )
                    << " " << a.a.pos1( ) << " " << a.a.Pos1( )
                    << " " << printSeq(a.p) << endl;    }    }

     auto Cat = [&]( const vec<int>& e )
     {    basevector x = tigs[ e[0] ];
          for ( int j = 1; j < e.isize( ); j++ )
               x = TrimCat( K, x, tigs[ e[j] ] );
          return x;    };

     for ( int i = 0; i < gs.isize( ); i++ )
     for ( int j = 0; j < allpaths_flat[i].isize( ); j++ )
     {    const vec<int>& p = allpaths_flat[i][j];
          basevector b = Cat( D.O( p[0] ) );
          for ( int j = 1; j < p.isize( ); j++ )
               b = TrimCat( K, b, Cat( D.O( p[j] ) ) );
          b.Print( fout, ToString(gs[i]) + "." + ToString(j) );    }

     // Print stats table.

     vec<vec<String>> rows;
     rows.push_back( vec<String> { "fin", "edges", "lines",
          "cov%", "errs", "cap", "uncap", "endgap" } );
     for ( int i = 0; i < gs.isize( ); i++ )
     {    int g = gs[i];
          vec<String> row;
          ostringstream xout, out;
          xout << "{" << setiosflags(ios::fixed) << setprecision(2)
               << round( ( 10000.0 * covstats[i].first ) 
                    / covstats[i].second ) / 100.0 << resetiosflags(ios::fixed) 
               << "," << stats[i].first/10.0 << "," << stats[i].second << ","
               << stats[i].third << "," << stats[i].fourth << "}";
          out << setiosflags(ios::fixed) << setprecision(2)
               << round( ( 10000.0 * covstats[i].first ) 
                    / covstats[i].second ) / 100.0 << resetiosflags(ios::fixed);
          row.push_back( ToString(g), ToString(counts[i].first), 
               ToString(counts[i].second), out.str( ), 
               ToString( stats[i].first/10.0 ), ToString(stats[i].second), 
               ToString(stats[i].third), ToString(stats[i].fourth) );
          rows.push_back(row);
          statslist[i].push_back( counts[i].first, counts[i].second,
               covstats[i].first, covstats[i].second, stats[i].first,
               stats[i].second, stats[i].third, stats[i].fourth );    }
     wout << endl;
     zout << endl;
     PrintTabular( wout, rows, 2, "rrrrrrrr" );
     PrintTabular( zout, rows, 2, "rrrrrrrr" );
     PrintTabular( sout, rows, 2, "rrrrrrrr" );

     // Print stats summary.

     sout << sout_str;
     qout << "= " << yout.str( ) << endl << endl;
     wout << qout.str( ), zout << qout.str( );
     // cout << "\nused " << TimeSince(clock) << endl << endl;
     cout << "For details, please see o.report." << endl << endl;
     return yout.str( );    }

typedef triple<ho_interval,int,int> tripto;
extern template class SmallVec<tripto,MempoolAllocator<tripto> >;
extern template class OuterVec< SerfVec<tripto> >;
#include "feudal/SmallVecDefs.h"
#include "feudal/OuterVecDefs.h"
template class SmallVec<tripto,MempoolAllocator<tripto> >;
template class OuterVec< SerfVec<tripto> >;
