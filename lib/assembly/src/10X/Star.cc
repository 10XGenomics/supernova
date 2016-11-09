// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>

#include "CoreTools.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "paths/long/large/Lines.h"
#include "10X/Heuristics.h"
#include "10X/IntIndex.h"
#include "10X/LineOO.h"
#include "10X/PlaceReads.h"
#include "10X/Star.h"
#include "10X/Super.h"

void Star( const HyperBasevectorX& hb, const vec<int>& inv, const vec<Bool>& dup,
     const vec<int32_t>& bc, const ReadPathVecX& paths, 
     digraphE<vec<int>>& D, vec<int>& dinv, 
     ReadPathVec& dpaths, const VecIntVec& ebcx, 
     const vec< triple<int,int,int> >& qept, 
     const double MIN_ADVANTAGE, const Bool OO, 
     const Bool DJANGO, const Bool verbose )
{
     cout << Date( ) << ": begin Star" << endl;

     // Heuristics.

     const int MIN_STAR = 5000;
     const double MAX_CN_DIFF = 0.5;
     const int MAX_VIEW = 10;
     int MAX_RIGHTS = 6;
     if (OO) MAX_RIGHTS = 4;
     const int MIN_BAR_TO = 2000;
     const int BC_VIEW = 50000;

     // Place reads.

     cout << Date( ) << ": placing reads" << endl;
     PlaceReads( hb, paths, dup, D, dpaths, True, False );
     cout << Date( ) << ": making index" << endl;
     IntIndex dpaths_index( dpaths, D.E( ) );

     // Find lines.

     cout << Date( ) << ": finding lines" << endl;
     vec<vec<vec<vec<int>>>> dlines;
     FindLines( D, dinv, dlines, MAX_CELL_PATHS, MAX_CELL_DEPTH, True );

     // Compute kmers.

     vec<int> kmers( hb.E( ) );
     #pragma omp parallel for
     for ( int e = 0; e < hb.E( ); e++ )
          kmers[e] = hb.Kmers(e);

     // Compute coverage.

     cout << Date( ) << ": mult before CN" << endl;
     vec<int> mult, llens;
     ComputeMult( hb, D, mult );
     GetLineLengths( hb, D, dlines, llens );
     vec<double> cov;
     {    vec<vec<pair<int,int>>> lbp;
          BarcodePos( bc, hb, D, dinv, dpaths, dpaths_index, dlines, lbp, 0 );
          MasterVec<SerfVec<pair<int,int>>> lbpx;
          for ( auto x : lbp )
          {    SerfVec<pair<int,int>> y( x.size( ) );
               for ( int j = 0; j < x.isize( ); j++ )
                    y[j] = x[j];
               lbpx.push_back(y);    }
          LineCN( kmers, lbpx, D, dlines, llens, cov );    }

     // Start passes.

     int sp = 0;
     while(1)
     {    cout << Date( ) << ": start star pass " << ++sp << endl;

          // Introduce barcode-only gaps.  First step: find line neighbors.

          cout << Date( ) << ": get line ancillary data" << endl;
          vec<int> llens, to_left, to_right, linv;
          vec< vec< pair<int,int> > > lhood;
          GetLineLengths( hb, D, dlines, llens );
          LineProx( hb, inv, ebcx, D, dinv, dlines, qept, lhood );
          D.ToLeft(to_left), D.ToRight(to_right);
          LineInv( dlines, dinv, linv );

          // Get barcode positions.

          cout << Date( ) << ": getting barcode positions" << endl;
          vec<vec<pair<int,int>>> lbp;
          BarcodePos( bc, hb, D, dinv, dpaths, dpaths_index, dlines, lbp, BC_VIEW );
     
          // Go through long lines and find their nearest right neighbors.

          cout << Date( ) << ": begin star stage" << endl;
          vec< pair<int,int> > stars;
          vec<int> l1s;
          for ( int L1 = 0; L1 < dlines.isize( ); L1++ )
          {    if ( llens[L1] < MIN_STAR ) continue;
               int v = to_right[ dlines[L1].back( )[0][0] ];
               if ( v < 0 ) continue;
               if ( D.From(v).nonempty( ) || !D.To(v).solo( ) ) continue;
               l1s.push_back(L1);    }
          vec<int> l1slen( l1s.isize( ) );
          for ( int i = 0; i < l1s.isize( ); i++ )
               l1slen[i] = llens[ l1s[i] ];
          ReverseSortSync( l1slen, l1s );
          cout << Date( ) << ": start star loop on " << l1s.size( )
               << " L1 values" << endl;
          map< vec<int>, double > memory;
          double clock = WallClockTime( );
          const int batch = Max( 1, l1s.isize( )/500 );
          #pragma omp parallel for schedule(dynamic, batch)
          for ( int li = 0; li < l1s.isize( ); li++ )
          {    int L1 = l1s[li];
               vec< triple<int,int,int> > M;

               // Set up logging.

               ostringstream* outp;
               if (verbose) outp = new ostringstream;
               auto DumpOut = [&]
               {    if (verbose)
                    {
                         #pragma omp critical
                         {    cout << outp->str( );    }    
                         delete outp;    }    };

               // Find right neighbors of L1.

               vec<int> rights;
               for ( int i = 0; i < Min( MAX_VIEW, lhood[L1].isize( ) ); i++ )
               {    int L2 = lhood[L1][i].second;
                    if ( AbsDiff( cov[L1], cov[L2] ) > MAX_CN_DIFF ) continue;
                    if ( llens[L2] < MIN_BAR_TO ) continue;
                    vec<double> scores;
                    scores.push_back( MemoryScoreOrder( 
                         { L2, L1 }, lbp, llens, M, memory ) );
                    scores.push_back( MemoryScoreOrder( 
                         { linv[L2], L1 }, lbp, llens, M, memory ) );
                    scores.push_back( MemoryScoreOrder( 
                         { L1, L2 }, lbp, llens, M, memory ) );
                    scores.push_back( MemoryScoreOrder( 
                         { L1, linv[L2] }, lbp, llens, M, memory ) );
                    vec<int> ids( 4, vec<int>::IDENTITY );
                    SortSync( scores, ids );
                    double ad = scores[1] - scores[0];
                    if ( ad < MIN_ADVANTAGE ) continue;
                    if ( ids[0] <= 1 ) continue;
                    if ( ids[0] == 2 ) rights.push_back(L2);
                    else rights.push_back( linv[L2] );    }

               // Cap number of rights.  Go with the biggest lines.

               if ( rights.isize( ) > MAX_RIGHTS )
               {    vec<int> lens( rights.size( ) );
                    for ( int i = 0; i < rights.isize( ); i++ )
                         lens[i] = llens[ rights[i] ];
                    ReverseSortSync( lens, rights );
                    rights.resize(MAX_RIGHTS);    }
               if (verbose)
               {    (*outp) << "\nlooking right from L" << L1 << endl;
                    (*outp) << "see " << printSeq(rights) << endl;    }
               if ( rights.empty( ) )
               {    DumpOut( );
                    continue;    }
     
               // Find leftmost right.

               int R = -1;
               if ( rights.size( ) == 1 ) R = rights[0];
               else if ( !DJANGO )
               {    vec<int> brights;
                    double ad;
                    if (OO)
                    {    MemoryOrderAndOrientN( rights, lbp, llens, 
                              linv, M, memory, ad, brights );    }
                    else MemoryOrderN( rights, lbp, llens, M, memory, ad, brights );
                    if (verbose)
                    {    (*outp) << "\nlooking right from L" << L1 << endl
                              << "see " << printSeq(rights) << endl
                              << "winner = L" << brights[0]
                              << " [" << ad << "]" << endl;    }
                    if ( ad < MIN_ADVANTAGE )
                    {    DumpOut( );
                         continue;    }
                    R = brights[0];    }
               else
               {    vec<int> lens( rights.size( ) );
                    for ( int i = 0; i < rights.isize( ); i++ )
                         lens[i] = llens[ rights[i] ];
                    ReverseSortSync( lens, rights );
                    while( rights.size( ) > 1 )
                    {    vec<int> brights;
                         double ad;
                         if (OO)
                         {    MemoryOrderAndOrientN( rights, lbp, llens, 
                                   linv, M, memory, ad, brights );    }
                         else 
                         {    MemoryOrderN( rights, lbp, llens, 
                                   M, memory, ad, brights );    }
                         if (verbose)
                         {    (*outp) << "\nlooking right from L" << L1 << endl
                                   << "see " << printSeq(rights) << endl
                                   << "winner = L" << brights[0]
                                   << " [" << ad << "]" << endl;    }
                         if ( ad >= MIN_ADVANTAGE )
                         {    R = brights[0];
                              break;    }
                         rights.pop_back( );    }    }

               // Check for done and otherwise save.

               if ( R == -1 )
               {    DumpOut( );
                    continue;    }
               int w = to_left[ dlines[R].front( )[0][0] ]; 
               if ( w < 0 )
               {    DumpOut( );
                    continue;    }
               if ( D.To(w).nonempty( ) || !D.From(w).solo( ) )
               {    DumpOut( );
                    continue;    }
               if (verbose) (*outp) << "saving " << L1 << " --> " << R << endl;
               DumpOut( );
               #pragma omp critical
               {    stars.push( L1, R );    }    }
     
          // Introduce joins.
     
          cout << Date( ) << ": done, time used = " << TimeSince(clock) << endl;
          Sort(stars);
          vec<pair<int,int>> stars2;
          for ( int i = 0; i < stars.isize( ); i++ )
          {    int L = stars[i].first, R = stars[i].second;
               if (verbose) cout << "\nexamining " << L << " --> " << R << endl;
               int ri = BinPosition( stars, make_pair( linv[R], linv[L] ) );
               if ( !DJANGO && ri < i ) continue;
               int d1 = dlines[L].back( )[0][0], d2 = dlines[R].front( )[0][0];
               int v = to_right[d1], w = to_left[d2];
               int rv = to_right[ dinv[d2] ], rw = to_left[ dinv[d1] ];
               if (DJANGO)
               {    if ( D.From(v).nonempty( ) || D.From(rv).nonempty( )
                         || D.To(w).nonempty( ) || D.To(rw).nonempty( ) )
                    {    continue;    }    }
               dinv.push_back( D.E( ) + 1, D.E( ) );
               D.AddEdge( v, w, {-2} );
               D.AddEdge( rv, rw, {-2} );
               stars2.push( L, R );
               stars2.push( linv[R], linv[L] );
               if (verbose)
               {    cout << "linking from L" << L << " to L" << R << endl;
                    cout << "linking from L" << linv[R] << " to L" << linv[L]
                         << endl;    }    }
          cout << Date( ) << ": made " << stars2.size( ) << " star joins" << endl;
          if ( stars2.empty( ) ) break;
          
          // Find lines of lines.

          vec<vec<int>> linelines;
          vec<int> ll_inv;
          vec<Bool> touched( dlines.size( ), False );
          {    Sort(stars2);
               vec<pair<int,int>> stars2r( stars2.size( ) );
               for ( int i = 0; i < stars2.isize( ); i++ )
                    stars2r[i] = make_pair( stars2[i].second, stars2[i].first );
               Sort(stars2r);
               for ( int i = 0; i < stars2.isize( ); i++ )
               {    int L = stars2[i].first, R = stars2[i].second;
                    if ( touched[L] || touched[R] ) continue;
                    vec<int> x = {L,R};
                    touched[L] = touched[R] = True;
                    Bool circle = False;
                    while(1)
                    {    int p = BinPosition1( stars2, x.back( ) );
                         if ( p < 0 ) break;
                         int R = stars2[p].second;
                         if ( R == x.front( ) ) 
                         {    circle = True;
                              break;    }
                         x.push_back(R);
                         touched[R] = True;    }
                    if ( !circle )
                    {    while(1)
                         {    int p = BinPosition1( stars2r, x.front( ) );
                              if ( p < 0 ) break;
                              int L = stars2r[p].second;
                              // if ( L == x.back( ) ) 
                              // {    circle = True;
                              //      break;    }
                              x.push_front(L);
                              touched[L] = True;    }    }
                    vec<int> y;
                    for ( int j = x.isize( ) - 1; j >= 0; j-- )
                    {    y.push_back( linv[ x[j] ] );
                         touched[ linv[ x[j] ] ] = True;    }
                    vec<int> xs(x), ys(y);
                    Sort(xs), Sort(ys);
                    if (circle)
                    {    x.push_back( x[0] ), y.push_back( y[0] );    }
                    int nll = linelines.size( );
                    linelines.push_back(x);
                    if ( ys != xs ) 
                    {    linelines.push_back(y);
                         ll_inv.push_back( nll+1, nll );    }
                    else ll_inv.push_back(nll);    }    }

          // Update lines and parallel data structures.

          cout << Date( ) << ": updating lines" << endl;
          vec<vec<int>> gap = { { } };
          for ( int i = 0; i < linelines.isize( ); i++ )
          {    const vec<int>& x = linelines[i];
               vec<vec<vec<int>>> y = dlines[ x[0] ];
               for ( int j = 1; j < x.isize( ); j++ )
               {    y.push_back(gap);
                    if ( x[j] == x[0] ) y.push_back( dlines[ x[j] ][0] );
                    else y.append( dlines[ x[j] ] );    }
               if ( ll_inv[i] == i - 1 )
               {    y = dlines.back( );
                    y.ReverseMe( );
                    for ( int j = 0; j < y.isize( ); j++ )
                    for ( int k = 0; k < y[j].isize( ); k++ )
                    for ( int l = 0; l < y[j][k].isize( ); l++ )
                         y[j][k][l] = dinv[ y[j][k][l] ];    }
               dlines.push_back(y);
               int64_t lensum = 0;
               double cnsum = 0.0;
               for ( int j = 0; j < x.isize( ); j++ )
               {    lensum += llens[ x[j] ];
                    cnsum += llens[ x[j] ] * cov[ x[j] ];    }
               cov.push_back( cnsum/lensum );
               touched.push_back(False);    }
          EraseIf( dlines, touched );
          EraseIf( cov, touched );    }

     // Clean up.

     cout << Date( ) << ": cleaning up at end of Star" << endl;
     CleanupCore( D, dinv );    }

/*

// Star2
// 1. given L, find its lefts and rights
// 2. take its largest right R
// 3. get lefts and rights of R
// 4. must have L in lefts(R)
// 5. any member of rights(L) >= 10 kb must be present in R+lefts(R)+rights(R)
// 6. any member of lefts(R) >= 10 kb must be present in L+lefts(L)+rights(L)
// 7. pool all, drop any lefts of L and rights of R
// 8. drop tiny stuff until <= 5
// 9. O&O, dropping tiny stuff until confident
// 10. save
// 11. traverse largest by Min(L,R) and don’t join if touched

void Star2( const HyperBasevectorX& hb, const vec<int>& inv, const vec<Bool>& dup,
     const vec<int32_t>& bc, const ReadPathVecX& paths, 
     digraphE<vec<int>>& D, vec<int>& dinv, 
     ReadPathVec& dpaths, const VecIntVec& ebcx, 
     const vec< triple<int,int,int> >& qept, 
     const double MIN_ADVANTAGE, const Bool OO, 
     const Bool DJANGO, const Bool verbose )
{
     cout << Date( ) << ": begin Star" << endl;

     // Heuristics.

     const int MIN_STAR = 5000;
     const double MAX_CN_DIFF = 0.5;
     const int MAX_VIEW = 10;
     int MAX_RIGHTS = 6;
     if (OO) MAX_RIGHTS = 4;
     const int MIN_BAR_TO = 2000;
     const int BC_VIEW = 50000;

     // Place reads.

     cout << Date( ) << ": placing reads" << endl;
     PlaceReads( hb, paths, dup, D, dpaths, True, False );
     cout << Date( ) << ": making index" << endl;
     IntIndex dpaths_index( dpaths, D.E( ) );

     // Find lines.

     cout << Date( ) << ": finding lines" << endl;
     vec<vec<vec<vec<int>>>> dlines;
     FindLines( D, dinv, dlines, MAX_CELL_PATHS, MAX_CELL_DEPTH, True );

     // Compute kmers.

     vec<int> kmers( hb.E( ) );
     #pragma omp parallel for
     for ( int e = 0; e < hb.E( ); e++ )
          kmers[e] = hb.Kmers(e);

     // Compute coverage.

     cout << Date( ) << ": mult before CN" << endl;
     vec<int> mult, llens;
     ComputeMult( hb, D, mult );
     GetLineLengths( hb, D, dlines, llens );
     vec<double> cov;
     {    vec<vec<pair<int,int>>> lbp;
          BarcodePos( bc, hb, D, dinv, dpaths, dpaths_index, dlines, lbp, 0 );
          MasterVec<SerfVec<pair<int,int>>> lbpx;
          for ( auto x : lbp )
          {    SerfVec<pair<int,int>> y( x.size( ) );
               for ( int j = 0; j < x.isize( ); j++ )
                    y[j] = x[j];
               lbpx.push_back(y);    }
          LineCN( kmers, lbpx, D, dlines, llens, cov );    }

     // Start passes.

     int sp = 0;
     while(1)
     {    cout << Date( ) << ": start star pass " << ++sp << endl;

          // Introduce barcode-only gaps.  First step: find line neighbors.

          cout << Date( ) << ": get line ancillary data" << endl;
          vec<int> llens, to_left, to_right, linv;
          vec< vec< pair<int,int> > > lhood;
          GetLineLengths( hb, D, dlines, llens );
          LineProx( hb, inv, ebcx, D, dinv, dlines, qept, lhood );
          D.ToLeft(to_left), D.ToRight(to_right);
          LineInv( dlines, dinv, linv );

          // Get barcode positions.

          cout << Date( ) << ": getting barcode positions" << endl;
          vec<vec<pair<int,int>>> lbp;
          BarcodePos( bc, hb, D, dinv, dpaths, dpaths_index, dlines, lbp, BC_VIEW );
     
          // Find long lines that terminate at their right end.

          cout << Date( ) << ": begin star stage" << endl;
          vec< pair<int,int> > stars;
          vec<int> l1s;
          for ( int L1 = 0; L1 < dlines.isize( ); L1++ )
          {    if ( llens[L1] < MIN_STAR ) continue;
               int v = to_right[ dlines[L1].back( )[0][0] ];
               if ( v < 0 ) continue;
               if ( D.From(v).nonempty( ) || !D.To(v).solo( ) ) continue;
               l1s.push_back(L1);    }
          vec<int> l1slen( l1s.isize( ) );
          for ( int i = 0; i < l1s.isize( ); i++ )
               l1slen[i] = llens[ l1s[i] ];
     
          // Go through long lines and find their nearest left and right neighbors.

          ReverseSortSync( l1slen, l1s );
          cout << Date( ) << ": start star loop on " << l1s.size( )
               << " L1 values" << endl;
          map< vec<int>, double > memory;
          double clock = WallClockTime( );
          const int batch = Max( 1, l1s.isize( )/500 );
          #pragma omp parallel for schedule(dynamic, batch)
          for ( int li = 0; li < l1s.isize( ); li++ )
          {    int L1 = l1s[li];
               vec< triple<int,int,int> > M;

               // Set up logging.

               ostringstream* outp;
               if (verbose) outp = new ostringstream;
               auto DumpOut = [&]
               {    if (verbose)
                    {
                         #pragma omp critical
                         {    cout << outp->str( );    }    
                         delete outp;    }    };

               // Find left and right neighbors of L1.

               vec<int> lefts, rights;
               for ( int i = 0; i < Min( MAX_VIEW, lhood[L1].isize( ) ); i++ )
               {    int L2 = lhood[L1][i].second;
                    if ( AbsDiff( cov[L1], cov[L2] ) > MAX_CN_DIFF ) continue;
                    if ( llens[L2] < MIN_BAR_TO ) continue;
                    vec<double> scores;
                    scores.push_back( MemoryScoreOrder( 
                         { L2, L1 }, lbp, llens, M, memory ) );
                    scores.push_back( MemoryScoreOrder( 
                         { linv[L2], L1 }, lbp, llens, M, memory ) );
                    scores.push_back( MemoryScoreOrder( 
                         { L1, L2 }, lbp, llens, M, memory ) );
                    scores.push_back( MemoryScoreOrder( 
                         { L1, linv[L2] }, lbp, llens, M, memory ) );
                    vec<int> ids( 4, vec<int>::IDENTITY );
                    SortSync( scores, ids );
                    double ad = scores[1] - scores[0];
                    if ( ad < MIN_ADVANTAGE ) continue;
                    if ( ids[0] == 0 ) lefts.push_back(L2);
                    else if ( ids[0] == 1 ) lefts.push_back( linv[L2] );
                    else if ( ids[0] == 2 ) rights.push_back(L2);
                    else rights.push_back( linv[L2] );    }

               // Find the largest right neighbor L2.

               if (verbose)
               {    (*outp) << "\nlooking right from L" << L1 << endl;
                    (*outp) << "see " << printSeq(rights) << endl;    }
               if ( rights.empty( ) )
               {    DumpOut( );
                    continue;    }
               int maxlen = -1, L2 = -1;
               for ( auto l : right )
               {    if ( llens[l] > maxlen )
                    {    maxlen = llens[l];
                         L2 = l;    }    }
               if (verbose) (*outp) << "largest right = " << L2 << endl;

====================================================================================

               // Cap number of rights.  Go with the biggest lines.

               if ( rights.isize( ) > MAX_RIGHTS )
               {    vec<int> lens( rights.size( ) );
                    for ( int i = 0; i < rights.isize( ); i++ )
                         lens[i] = llens[ rights[i] ];
                    ReverseSortSync( lens, rights );
                    rights.resize(MAX_RIGHTS);    }
               if (verbose)
               {    (*outp) << "\nlooking right from L" << L1 << endl;
                    (*outp) << "see " << printSeq(rights) << endl;    }
               if ( rights.empty( ) )
               {    DumpOut( );
                    continue;    }
     
               // Find leftmost right.

               int R = -1;
               if ( rights.size( ) == 1 ) R = rights[0];
               else if ( !DJANGO )
               {    vec<int> brights;
                    double ad;
                    if (OO)
                    {    MemoryOrderAndOrientN( rights, lbp, llens, 
                              linv, M, memory, ad, brights );    }
                    else MemoryOrderN( rights, lbp, llens, M, memory, ad, brights );
                    if (verbose)
                    {    (*outp) << "\nlooking right from L" << L1 << endl
                              << "see " << printSeq(rights) << endl
                              << "winner = L" << brights[0]
                              << " [" << ad << "]" << endl;    }
                    if ( ad < MIN_ADVANTAGE )
                    {    DumpOut( );
                         continue;    }
                    R = brights[0];    }
               else
               {    vec<int> lens( rights.size( ) );
                    for ( int i = 0; i < rights.isize( ); i++ )
                         lens[i] = llens[ rights[i] ];
                    ReverseSortSync( lens, rights );
                    while( rights.size( ) > 1 )
                    {    vec<int> brights;
                         double ad;
                         if (OO)
                         {    MemoryOrderAndOrientN( rights, lbp, llens, 
                                   linv, M, memory, ad, brights );    }
                         else 
                         {    MemoryOrderN( rights, lbp, llens, 
                                   M, memory, ad, brights );    }
                         if (verbose)
                         {    (*outp) << "\nlooking right from L" << L1 << endl
                                   << "see " << printSeq(rights) << endl
                                   << "winner = L" << brights[0]
                                   << " [" << ad << "]" << endl;    }
                         if ( ad >= MIN_ADVANTAGE )
                         {    R = brights[0];
                              break;    }
                         rights.pop_back( );    }    }

               // Check for done and otherwise save.

               if ( R == -1 )
               {    DumpOut( );
                    continue;    }
               int w = to_left[ dlines[R].front( )[0][0] ]; 
               if ( w < 0 )
               {    DumpOut( );
                    continue;    }
               if ( D.To(w).nonempty( ) || !D.From(w).solo( ) )
               {    DumpOut( );
                    continue;    }
               if (verbose) (*outp) << "saving " << L1 << " --> " << R << endl;
               DumpOut( );
               #pragma omp critical
               {    stars.push( L1, R );    }    }
     
          // Introduce joins.
     
          cout << Date( ) << ": done, time used = " << TimeSince(clock) << endl;
          Sort(stars);
          vec<pair<int,int>> stars2;
          for ( int i = 0; i < stars.isize( ); i++ )
          {    int L = stars[i].first, R = stars[i].second;
               if (verbose) cout << "\nexamining " << L << " --> " << R << endl;
               int ri = BinPosition( stars, make_pair( linv[R], linv[L] ) );
               if ( !DJANGO && ri < i ) continue;
               int d1 = dlines[L].back( )[0][0], d2 = dlines[R].front( )[0][0];
               int v = to_right[d1], w = to_left[d2];
               int rv = to_right[ dinv[d2] ], rw = to_left[ dinv[d1] ];
               if (DJANGO)
               {    if ( D.From(v).nonempty( ) || D.From(rv).nonempty( )
                         || D.To(w).nonempty( ) || D.To(rw).nonempty( ) )
                    {    continue;    }    }
               dinv.push_back( D.E( ) + 1, D.E( ) );
               D.AddEdge( v, w, {-2} );
               D.AddEdge( rv, rw, {-2} );
               stars2.push( L, R );
               stars2.push( linv[R], linv[L] );
               if (verbose)
               {    cout << "linking from L" << L << " to L" << R << endl;
                    cout << "linking from L" << linv[R] << " to L" << linv[L]
                         << endl;    }    }
          cout << Date( ) << ": made " << stars2.size( ) << " star joins" << endl;
          if ( stars2.empty( ) ) break;
          
          // Find lines of lines.

          vec<vec<int>> linelines;
          vec<int> ll_inv;
          vec<Bool> touched( dlines.size( ), False );
          {    Sort(stars2);
               vec<pair<int,int>> stars2r( stars2.size( ) );
               for ( int i = 0; i < stars2.isize( ); i++ )
                    stars2r[i] = make_pair( stars2[i].second, stars2[i].first );
               Sort(stars2r);
               for ( int i = 0; i < stars2.isize( ); i++ )
               {    int L = stars2[i].first, R = stars2[i].second;
                    if ( touched[L] || touched[R] ) continue;
                    vec<int> x = {L,R};
                    touched[L] = touched[R] = True;
                    Bool circle = False;
                    while(1)
                    {    int p = BinPosition1( stars2, x.back( ) );
                         if ( p < 0 ) break;
                         int R = stars2[p].second;
                         if ( R == x.front( ) ) 
                         {    circle = True;
                              break;    }
                         x.push_back(R);
                         touched[R] = True;    }
                    if ( !circle )
                    {    while(1)
                         {    int p = BinPosition1( stars2r, x.front( ) );
                              if ( p < 0 ) break;
                              int L = stars2r[p].second;
                              // if ( L == x.back( ) ) 
                              // {    circle = True;
                              //      break;    }
                              x.push_front(L);
                              touched[L] = True;    }    }
                    vec<int> y;
                    for ( int j = x.isize( ) - 1; j >= 0; j-- )
                    {    y.push_back( linv[ x[j] ] );
                         touched[ linv[ x[j] ] ] = True;    }
                    vec<int> xs(x), ys(y);
                    Sort(xs), Sort(ys);
                    if (circle)
                    {    x.push_back( x[0] ), y.push_back( y[0] );    }
                    int nll = linelines.size( );
                    linelines.push_back(x);
                    if ( ys != xs ) 
                    {    linelines.push_back(y);
                         ll_inv.push_back( nll+1, nll );    }
                    else ll_inv.push_back(nll);    }    }

          // Update lines and parallel data structures.

          cout << Date( ) << ": updating lines" << endl;
          vec<vec<int>> gap = { { } };
          for ( int i = 0; i < linelines.isize( ); i++ )
          {    const vec<int>& x = linelines[i];
               vec<vec<vec<int>>> y = dlines[ x[0] ];
               for ( int j = 1; j < x.isize( ); j++ )
               {    y.push_back(gap);
                    if ( x[j] == x[0] ) y.push_back( dlines[ x[j] ][0] );
                    else y.append( dlines[ x[j] ] );    }
               if ( ll_inv[i] == i - 1 )
               {    y = dlines.back( );
                    y.ReverseMe( );
                    for ( int j = 0; j < y.isize( ); j++ )
                    for ( int k = 0; k < y[j].isize( ); k++ )
                    for ( int l = 0; l < y[j][k].isize( ); l++ )
                         y[j][k][l] = dinv[ y[j][k][l] ];    }
               dlines.push_back(y);
               int64_t lensum = 0;
               double cnsum = 0.0;
               for ( int j = 0; j < x.isize( ); j++ )
               {    lensum += llens[ x[j] ];
                    cnsum += llens[ x[j] ] * cov[ x[j] ];    }
               cov.push_back( cnsum/lensum );
               touched.push_back(False);    }
          EraseIf( dlines, touched );
          EraseIf( cov, touched );    }

     // Clean up.

     cout << Date( ) << ": cleaning up at end of Star" << endl;
     CleanupCore( D, dinv );    }

*/
