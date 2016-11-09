// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

// Estimate sizes of barcode-only gaps.
//
// results for Mr.X: mean error = 2168, median err = 447
//
// updated:
// errs = 10204, mean error = 1976, median err = 447
// fraction of gap estimates off by >= 50 kb = 0.304%

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"
#include "10X/Gap.h"
#include "10X/Gaprika.h"
#include "10X/Super.h"
#include "10X/astats/RefAlign.h"

void Gaprika( const vec<int>& kmers, digraphE<vec<int>>& D, const vec<int>& dinv,
     const vec<vec<vec<vec<int>>>>& dlines, 
     const MasterVec<SerfVec<pair<int,int>>>& lbpx, const int MAX_GAP,
     const MasterVec<SerfVec<refalign>>& galigns, 
     const int VERBOSITY, const Bool EVAL )
{
     // Heuristics.

     const int WINDOW = 10000;
     const int MULT = 50;
     const int GAP_DELTA = 5000;
     const int MIN_GAP = 400;
     const int MIN_POINTS = 2;

     // Set up.

     vec<int> to_left, to_right;
     D.ToLeft(to_left), D.ToRight(to_right);

     // Find line lengths and barcode-only gap positions.

     double clock = WallClockTime( );
     cout << Date( ) << ": getting line data" << endl;
     vec<int> xpos( dlines.size( ), 0 );
     vec<vec<int>> bc_gap_pos( dlines.size( ) );
     vec<vec<pair<int,int>>> bc_gap_loc( dlines.size( ) );
     int print_count = 0;
     int changes = 0;
     int wayoff = 0, wayoff_total = 0;
     #pragma omp parallel for schedule(dynamic,10)
     for ( int l = 0; l < dlines.isize( ); l++ )
     {    const vec<vec<vec<int>>>& L = dlines[l];
          for ( int m = 0; m < L.isize( ); m++ )
          {    const vec<vec<int>>& M = L[m];
               vec<int> xlens;
               if ( M.solo( ) && M[0].empty( ) )
               {    if ( m == 0 ) continue; // shouldn't happen
                    if ( L[m-1].empty( ) ) continue; // shouldn't happen
                    int v = to_right[ L[m-1][0].back( ) ];
                    if ( !D.From(v).solo( ) ) continue; // shouldn't happen
                    int d = D.IFrom(v,0);
                    if ( IsBarcodeOnlyGap( D.O(d) ) )
                    {    bc_gap_pos[l].push_back(xpos[l]);
                         bc_gap_loc[l].push( d, xpos[l] );    }    }
               for ( int k = 0; k < M.isize( ); k++ )
               {    int len = 0;
                    for ( auto d : M[k] )
                    {    if ( D.O(d)[0] < 0 )
                         {    if ( IsBarcodeOnlyGap( D.O(d) ) )
                                   bc_gap_pos[l].push_back(xpos[l]);
                              continue;    }
                         for ( auto e : D.O(d) ) len += kmers[e];    }
                    xlens.push_back(len);    }
               Sort(xlens);
               if ( xlens.nonempty( ) ) xpos[l] += Median(xlens);    }
          UniqueSort(bc_gap_pos[l]);    }

     // Range over gap sizes.

     vec<vec<int>> buckets;
     BucketLines( dlines, xpos, buckets, 2*WINDOW );
     cout << Date( ) << ": starting loop" << endl;
     if ( VERBOSITY >= 1 ) cout << endl;
     vec<pair<int,double>> gb;
     for ( int gap = 0; gap <= MAX_GAP; gap += GAP_DELTA )
     {    vec<int> nbridges;
          vec<double> fracs;
          #pragma omp parallel for schedule(dynamic,1)
          for ( int bu = 0; bu < buckets.isize( ); bu++ )
          {    vec< pair<int,int>> pb;
               vec<int> bcount, lefts, rights, lr;
               vec<double> fr;
               for ( auto l : buckets[bu] )
               {    if ( xpos[l] < 2*WINDOW + gap ) continue;
                    const vec<vec<vec<int>>>& L = dlines[l];
                    pb.clear( ), bcount.clear( );
                    fr.clear( );
                    for ( auto x : lbpx[l] ) pb.push( x.second, x.first );
                    Sort(pb);
                    for ( int i = WINDOW; i <= xpos[l] - WINDOW - gap; 
                         i += WINDOW * MULT )
                    {    int left1 = i - WINDOW, left2 = i;
                         int right1 = i + gap, right2 = i + WINDOW + gap;
                         Bool at_gap = False;
                         for ( auto b : bc_gap_pos[l] ) 
                              if ( b >= left1 && b <= right2 ) at_gap = True;
                         if (at_gap) continue;
                         int low = LowerBound1( pb, left1 ); 
                         lefts.clear( ), rights.clear( );
                         for ( int j = low; j < pb.isize( ); j++ )
                         {    int p = pb[j].first;
                              if ( p > right2 ) break;
                              int b = pb[j].second;
                              if ( p >= left1 && p < left2 ) lefts.push_back(b);
                              if ( p >= right1 && p < right2 ) 
                                   rights.push_back(b);    }

                         Sort(lefts), Sort(rights);
                         vec<Bool> del_left( lefts.size( ), False );
                         for ( int j = 0; j < lefts.isize( ); j++ )
                         {    int k;
                              for ( k = j + 1; k < lefts.isize( ); k++ )
                                   if ( lefts[k] != lefts[j] ) break;
                              if ( k - j < MIN_POINTS )
                                   for ( int l = j; l < k; l++ ) del_left[l] = True;
                              j = k - 1;    }
                         EraseIf( lefts, del_left );
                         vec<Bool> del_right( rights.size( ), False );
                         for ( int j = 0; j < rights.isize( ); j++ )
                         {    int k;
                              for ( k = j + 1; k < rights.isize( ); k++ )
                                   if ( rights[k] != rights[j] ) break;
                              if ( k - j < MIN_POINTS )
                                   for ( int l = j; l < k; l++ ) del_right[l] = True;
                              j = k - 1;    }
                         EraseIf( rights, del_right );

                         UniqueSort(lefts), UniqueSort(rights);
                         lr = lefts;
                         lr.append(rights);
                         UniqueSort(lr);
                         int bridges = MeetSize( lefts, rights );
                         if ( lr.empty( ) ) continue;
                         double bridge_frac = double(bridges) / lr.size( );
                         bcount.push_back(bridges);    
                         fr.push_back(bridge_frac);    }
                    #pragma omp critical
                    {    nbridges.append(bcount);    
                         fracs.append(fr);    }    }    }
          Sort(fracs);
          int count = fracs.size( );
          if ( count > 0 )
          {    double m = Mean(fracs);
               double dev = StdDev( fracs, m );
               if ( VERBOSITY >= 1 ) PRINT5( gap, count, m, dev, Median(fracs) );    
               gb.push( gap, m );    }    }
     if ( VERBOSITY >= 1 ) cout << endl;
     cout << Date( ) << ": done, time used = " << TimeSince(clock) << endl;

     // Now estimate sizes of bc-only gaps.

     cout << Date( ) << ": estimating gap sizes" << endl;
     if ( VERBOSITY >= 1 ) cout << endl;
     vec<int> errs;
     #pragma omp parallel for schedule(dynamic,10)
     for ( int l = 0; l < dlines.isize( ); l++ )
     {    const vec<vec<vec<int>>>& L = dlines[l];
          vec< pair<int,int>> pb;
          for ( auto x : lbpx[l] ) pb.push( x.second, x.first );
          Sort(pb);
          vec<int> lr;
          for ( auto x : bc_gap_loc[l] )
          {    int d = x.first, pos = x.second;
               vec<int> bcount, lefts, rights;
               vec<double> fracs;
               int left1 = pos - WINDOW, left2 = pos;
               int right1 = pos, right2 = pos + WINDOW;
               if ( left1 < 0 || right2 > xpos[l] ) continue; // FAIL
               Bool at_gap = False;
               for ( auto b : bc_gap_pos[l] ) 
                    if ( b != pos && b >= left1 && b <= right2 ) at_gap = True;
               if (at_gap) continue; // FAIL
               int low = LowerBound1( pb, left1 ); 
               for ( int j = low; j < pb.isize( ); j++ )
               {    int p = pb[j].first;
                    if ( p > right2 ) break;
                    int b = pb[j].second;
                    if ( p >= left1 && p < left2 ) lefts.push_back(b);
                    if ( p >= right1 && p < right2 ) rights.push_back(b);    }

               Sort(lefts), Sort(rights);
               vec<Bool> del_left( lefts.size( ), False );
               for ( int j = 0; j < lefts.isize( ); j++ )
               {    int k;
                    for ( k = j + 1; k < lefts.isize( ); k++ )
                         if ( lefts[k] != lefts[j] ) break;
                    if ( k - j < MIN_POINTS )
                         for ( int l = j; l < k; l++ ) del_left[l] = True;
                    j = k - 1;    }
               EraseIf( lefts, del_left );
               vec<Bool> del_right( rights.size( ), False );
               for ( int j = 0; j < rights.isize( ); j++ )
               {    int k;
                    for ( k = j + 1; k < rights.isize( ); k++ )
                         if ( rights[k] != rights[j] ) break;
                    if ( k - j < MIN_POINTS )
                         for ( int l = j; l < k; l++ ) del_right[l] = True;
                    j = k - 1;    }
               EraseIf( rights, del_right );

               UniqueSort(lefts), UniqueSort(rights);
               lr = lefts;
               lr.append(rights);
               UniqueSort(lr);
               int bridges = MeetSize( lefts, rights );
               if ( lr.empty( ) ) continue;
               double bridge_frac = double(bridges) / lr.size( );
               if ( VERBOSITY >= 2 )
               {
                    #pragma omp critical
                    {    PRINT6( l, d, lefts.size( ), rights.size( ),
                              bridges, bridge_frac );    }    }

               // If linking is too weak, don't set gap.  In some cases these
               // instances would be due to misassemblies.

               if ( gb.nonempty( ) && bridge_frac < gb.back( ).second/2 ) continue;

               int best_gap = -1;
               double score = 1000000000;
               for ( int i = 0; i < gb.isize( ); i++ )
               {    int gap = gb[i].first;
                    double count = gb[i].second;
                    if ( Abs(count-bridge_frac) < score )
                    {    score = Abs(count-bridge_frac);
                         best_gap = gap;    }    }
               if ( best_gap < 0 ) continue;
               best_gap = Max( best_gap, MIN_GAP );

               // Set the gap.

               int rd = dinv[d];
               if ( d < rd ) 
               {    ForceAssert( D.O(d).solo( ) );
                    D.OMutable(d).push_back(best_gap);
                    D.OMutable(rd).push_back(best_gap);
                    #pragma omp critical
                    {    changes += 2;    }    }

               // Assess gap accuracy.

               if (EVAL)
               {    int v = to_left[d], w = to_right[d];
                    int d1 = -1, d2 = -1;
                    if ( D.To(v).solo( ) ) d1 = D.ITo(v,0);
                    if ( D.From(w).solo( ) ) d2 = D.IFrom(w,0);
                    const int MAX_PRINT = 10;
                    #pragma omp critical
                    {    if ( d1 >= 0 && d2 >= 0 && galigns[d1].solo( )
                              && galigns[d2].solo( ) 
                              && galigns[d1][0].chr == galigns[d2][0].chr )
                         {    int chr = galigns[d1][0].chr;
                              const align& a1 = galigns[d1][0].a; 
                              const align& a2 = galigns[d2][0].a;
                              if ( chr > 0 )
                              {    int ref_gap = a2.pos2( ) - a1.Pos2( );
                                   if ( Abs(ref_gap) < 300000 )
                                   {    int err = Abs( best_gap - ref_gap );
                                        wayoff_total++;
                                        if ( err >= 50000 ) wayoff++;
                                        if ( randomx() % 100 == 0 
                                             && ++print_count <= MAX_PRINT )
                                        {    cout << "l = " << l << ", d = " << d
                                                  << ", predicted gap = " << best_gap
                                                  << ", ref gap = " << ref_gap 
                                                  << ", err = " << err << endl;    }
                                        errs.push_back(err);    
                                             }    }    }    }    }    }    }
     if (EVAL)
     {    Sort(errs);
          int mean_err = int(round(Mean(errs)));
          cout << endl << "errs = " << errs.size( ) << ", mean error = " << mean_err
               << ", median err = " << Median(errs) << endl;
          cout << "fraction of gap estimates off by >= 50 kb = "
               << PERCENT_RATIO( 3, wayoff, wayoff_total ) << endl << endl;    }
     cout << Date( ) << ": really done, " << changes << " gap values set" 
          << endl;    }
