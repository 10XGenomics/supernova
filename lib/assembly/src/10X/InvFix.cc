// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

// Flip parts of lines that seem inverted.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"
#include "paths/HyperBasevector.h"
#include "paths/long/large/Lines.h"
#include "10X/Gap.h"
#include "10X/Heuristics.h"
#include "10X/InvFix.h"
#include "10X/Super.h"
#include "10X/astats/RefAlign.h"

void InvFix( const HyperBasevectorX& hb, const vec<int>& inv,
     const vec<int>& kmers, digraphE<vec<int>>& D, const vec<int>& dinv,
     vec<vec<vec<vec<int>>>>& dlines,
     const MasterVec<SerfVec<pair<int,int>>>& lbpx,
     const MasterVec<SerfVec<refalign>>& galigns )
{
     // Set up.

     vec<int> to_left, to_right, linv;
     D.ToLeft(to_left), D.ToRight(to_right);
     LineInv( dlines, dinv, linv );
     const Bool verbose = False;

     // Count inversions versus reference.

     auto InvCount = [&]( )
     {    cout << Date( ) << ": counting inversions versus reference" << endl;
          int invcount = 0;
          #pragma omp parallel for
          for ( int l = 0; l < dlines.isize( ); l++ )
          {    const vec<vec<vec<int>>>& L = dlines[l];
               vec<refalign> R;
               for ( int m = 0; m < L.isize( ); m += 2 )
               {    int d = L[m][0][0];
                    if ( galigns[d].size( ) == 1 ) R.push_back( galigns[d][0] );    }
               for ( int i = 1; i < R.isize( ); i++ )
               {    const refalign &r1 = R[i-1], &r2 = R[i];
                    if ( Abs(r1.chr) == Abs(r2.chr) && r1.chr != r2.chr ) 
                    {
                         #pragma omp critical
                         {    invcount++;    }    }    }    }
          DPRINT(invcount);    };
     if ( galigns.size( ) > 0 ) InvCount( );
     
     // Find line lengths and barcode-only gap positions.

     cout << Date( ) << ": getting line data" << endl;
     vec<int> xpos( dlines.size( ), 0 );
     vec<vec<int>> bc_gap_pos( dlines.size( ) ), bc_gap_loc( dlines.size( ) );
     vec<vec<int>> bc_gap_id( dlines.size( ) );
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
                         bc_gap_loc[l].push(m);
                         bc_gap_id[l].push(d);    }    }
               for ( int k = 0; k < M.isize( ); k++ )
               {    int len = 0;
                    for ( auto d : M[k] )
                    {    if ( D.O(d)[0] < 0 ) continue;
                         for ( auto e : D.O(d) ) len += kmers[e];    }
                    xlens.push_back(len);    }
               Sort(xlens);
               if ( xlens.nonempty( ) ) xpos[l] += Median(xlens);    }    }

     // Scan for pairs of barcode-only gaps.

     cout << Date( ) << ": scanning" << endl;
     const int WINDOW = 10000;
     #pragma omp parallel for schedule(dynamic,10)
     for ( int l = 0; l < dlines.isize( ); l++ )
     {    if ( linv[l] <= l ) continue;
          const vec<vec<vec<int>>>& L = dlines[l];
          vec<pair<int,int>> pb;
          for ( auto x : lbpx[l] ) pb.push( x.second, x.first );
          Sort(pb);

          // Define scoring function.

          auto Score = [&]( const int j1, const int j2 )
          {    int start = bc_gap_pos[l][j1], stop = bc_gap_pos[l][j2];
               vec<int> left1fw, right1fw, left2fw, right2fw;
               int low = LowerBound1( pb, start - WINDOW );
               int high = UpperBound1( pb, stop + WINDOW );
               for ( int i = low; i < high; i++ )
               {    int p = pb[i].first, b = pb[i].second;
                    if ( p < start && p >= start - WINDOW ) left1fw.push_back(b);
                    if ( p >= start && p < start + (stop-start)/2 )
                         right1fw.push_back(b);
                    if ( p < stop && p >= stop - (stop-start)/2 )
                         left2fw.push_back(b);
                    if ( p >= stop && p < stop + WINDOW ) right2fw.push_back(b);    }
               UniqueSort(left1fw), UniqueSort(right1fw);
               UniqueSort(left2fw), UniqueSort(right2fw);
               int n1 = MeetSize(left1fw, right1fw) + MeetSize(left2fw, right2fw);
               int n2 = MeetSize(left1fw, left2fw) + MeetSize(right1fw, right2fw);
               return n2 - n1;    };

          // Look for inversions.

          for ( int j1 = 0; j1 < bc_gap_pos[l].isize( ) - 1; j1++ )
          {    int add = 1;
               int j2 = j1 + add;
               int score = Score( j1, j2 );
               if ( score > 0 )
               {    int m1 = bc_gap_loc[l][j1], m2 = bc_gap_loc[l][j2];
                    int d1 = bc_gap_id[l][j1], d2 = bc_gap_id[l][j2];
                    int rd1 = dinv[d1], rd2 = dinv[d2];
                    int v1 = to_left[d1], w1 = to_right[d2];
                    if ( !D.To(v1).solo( ) || !D.From(w1).solo( ) ) continue;
                    int v2 = to_left[rd2], w2 = to_right[rd1];
                    int f1 = D.ITo(v1,0), g1 = D.IFrom(w1,0);
                    int f2 = D.ITo(v2,0), g2 = D.IFrom(w2,0);
                    D.GiveEdgeNewToVx( f1, v1, v2 );
                    D.GiveEdgeNewToVx( f2, v2, v1 );
                    D.GiveEdgeNewFromVx( g1, w1, w2 );
                    D.GiveEdgeNewFromVx( g2, w2, w1 );
                    if (verbose)
                    {
                         #pragma omp critical
                         {    cout << "L=" << l << "." << m1 << "-" << m2 
                                   << ", add = " << add
                                   << ", score = " << score << endl;    }    }
     
                    // Update pb.

                    int start = bc_gap_pos[l][j1], stop = bc_gap_pos[l][j2];
                    int low = LowerBound1(pb, start), high = UpperBound1(pb, stop);
                    for ( int j = low; j < high; j++ )
                         pb[j].first = stop - start - pb[j].first;
                    for ( int j = low; j < low + (high-low)/2; j++ )
                         swap( pb[j], pb[ high - 1 - j ] );

                    // Advance.

                    j1 = j2 + 1;
                    while( j1 < bc_gap_pos[l].isize( ) - 1 
                         && bc_gap_pos[l][j1] - bc_gap_pos[l][j2] < WINDOW ) 
                    {    j1++;    }
                    j1--;    }    }    }

     // Regenerate lines and look for inversions versus reference.

     cout << Date( ) << ": regenerating lines" << endl;
     FindLines( D, dinv, dlines, MAX_CELL_PATHS, MAX_CELL_DEPTH );
     if ( galigns.size( ) > 0 ) InvCount( );    }
