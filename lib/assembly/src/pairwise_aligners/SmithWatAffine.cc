///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// SmithWatAffine( S, T )
//
// Not optimized

#include "Basevector.h"
#include "math/Functions.h"
#include "math/Array.h"

#include "Alignment.h"
#include "PackAlign.h"
#include "ShortVector.h"
#include "pairwise_aligners/SmithWatAffine.h"
#include "PrintAlignment.h"

#include "kmers/naif_kmer/Kmers.h"
#include "kmers/naif_kmer/NaifKmerizer.h"
#include "kmers/naif_kmer/KernelKmerStorer.h"
#include "kmers/naif_kmer/KmerMap.h"

#define FIX_RIGHT_GAP 1
//#undef FIX_RIGHT_GAP

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

namespace {
const int Infinity = 100000000;

// Unmatched ends in source read are not penalized
// Unmatched ends in target read are optionally penalized
unsigned int SmithWatAffineBandedCore( const basevector& S, const basevector& T,
			     align* a, int offset, int bandwidth,
                             bool penalize_left_gap, bool penalize_right_gap,
                             const int mismatch_penalty,
                             const int gap_open_penalty,
                             const int gap_extend_penalty )
{
     ForceAssertGt( S.size(), 0u );
     ForceAssertGt( T.size(), 0u );

     int n = S.size( ), N = T.size( );

     avector<char> s, t;
     s.resize(n);
     for (int i = 0; i < n; i++)
          s(i) = S[i];
     t.resize(N);
     for (int i = 0; i < N; i++)
          t(i) = T[i];

     const int Infinity = 100000000;
     int best_score = Infinity;
     RecArray<unsigned int> score_x(n+1, N+1, Infinity);
     RecArray<unsigned int> score_y(n+1, N+1, Infinity);
     RecArray<unsigned int> score_z(n+1, N+1, Infinity);

     RecArray<unsigned char> x_from(n+1, N+1, 's');
     RecArray<unsigned char> y_from(n+1, N+1, 's');
     RecArray<unsigned char> z_from(n+1, N+1, 's');

     // y and z are for the gap extending in target and source, respectively
     // x is for base substitution
     score_x[0][0] = 0;
     score_y[0][0] = Infinity;
     score_z[0][0] = Infinity;
     for (int i = 1; i <= n; i++ ) {
         if ( 0 < i - offset - bandwidth || 0  > i - offset + bandwidth)
             continue;
         score_z[i][0] = 0;
     }
     if (!penalize_left_gap) {
         for (int j = 1; j <= N; j++) {
             if (j < 0 - offset - bandwidth || j > 0 - offset + bandwidth)
                 continue;
             score_y[0][j] = 0;
         }
     }
     for (int i = 1; i <= n; i++ )
     {    int low = Max( 1, i - offset - bandwidth );
          int high = Min( N, i - offset + bandwidth );
          for (int j = low; j <= high; j++ )
	  {
               // if (j < i - offset - bandwidth || j > i - offset + bandwidth)
               //     continue;

               unsigned int x_x = score_x[i-1][j-1] + mismatch_penalty * ( s(i-1) != t(j-1) );
	       unsigned int x_y = score_y[i-1][j-1] + mismatch_penalty * ( s(i-1) != t(j-1) );
	       unsigned int x_z = score_z[i-1][j-1] + mismatch_penalty * ( s(i-1) != t(j-1) );
	       unsigned int y_x = score_x[i][j-1] + (penalize_right_gap || i != n ? gap_open_penalty : 0);
	       unsigned int y_y = score_y[i][j-1] + (penalize_right_gap || i != n ? gap_extend_penalty : 0);
	       unsigned int y_z = Infinity; //score_z[i][j-1] + gap_open_penalty;
	       unsigned int z_x = score_x[i-1][j] + (j != N ? gap_open_penalty : 0);
	       unsigned int z_y = Infinity; //score_y[i-1][j] + gap_open_penalty;
	       unsigned int z_z = score_z[i-1][j] + (j != N ? gap_extend_penalty: 0);

	       score_x[i][j] = Min( Min( x_x, x_y ), x_z );
	       score_y[i][j] = Min( Min( y_x, y_y ), y_z );
	       score_z[i][j] = Min( Min( z_x, z_y ), z_z );
               if ((int)score_x[i][j] > Infinity) score_x[i][j] = Infinity;

	       if ( x_x <= x_y )
	       {    if ( x_x <= x_z ) x_from[i][j] = 'x';
	            else x_from[i][j] =  'z';    }
	       else
	       {    if ( x_y <= x_z ) x_from[i][j] = 'y';
	            else x_from[i][j] =  'z';    }

	       if ( y_x <= y_y )
	       {    if ( y_x <= y_z ) y_from[i][j] = 'x';
	            else y_from[i][j] =  'z';    }
	       else
	       {    if ( y_y <= y_z ) y_from[i][j] = 'y';
	            else y_from[i][j] =  'z';    }

	       if ( z_x <= z_y )
	       {    if ( z_x <= z_z ) z_from[i][j] = 'x';
	            else z_from[i][j] =  'z';    }
	       else
	       {    if ( z_y <= z_z ) z_from[i][j] = 'y';
	            else z_from[i][j] =  'z';    }    }    }

     best_score = Min( score_x[n][N], Min( score_y[n][N], score_z[n][N] ) );

     // Find the end point of the best matching trace
     int ii = n, jj = N;
     for(int k = n; k >= 0; k--){
         int best_score_k = Min( score_x[k][N], Min( score_y[k][N], score_z[k][N] ) );
         if (best_score_k < best_score) {
             best_score = best_score_k;
             ii = k;
             jj = N;
         }
     }
     if (!penalize_right_gap) {
     for(int k = N; k >= 0; k--){
         int best_score_k = Min( score_x[n][k], Min( score_y[n][k], score_z[n][k] ) );
         if (best_score_k < best_score) {
             best_score = best_score_k;
             ii = n;
             jj = k;
         }
     }
     }
     if (best_score == Infinity) return best_score;

     //cout << "ii= " << ii << endl;
     //cout << "jj= " << jj << endl;
     //cout << "best_score= " << best_score << endl;
     //vec< vec<unsigned char> > *from;
     RecArray<unsigned char> *from;
     if ( score_x[ii][jj] <= score_y[ii][jj] )
     {    if ( score_x[ii][jj] <= score_z[ii][jj] ) from = &x_from;
          else from = &z_from;    }
     else
     {    if ( score_y[ii][jj] <= score_z[ii][jj] ) from = &y_from;
          else from = &z_from;    }
     int i = ii;
     int j = jj;

     int lcount = 0, g1count = 0, g2count = 0;
     int last_length = 0;
     avector<int> gaps(0), lengths(0);
     while(1)
     {
          unsigned char dir = (*from)[i][j];
	  //cout << dir;
          if ( from == &x_from )
          {    if ( g1count > 0 )
               {    if ( last_length > 0 )
		    {    gaps.Prepend( g1count );
		         lengths.Prepend( last_length );    }
	            g1count = 0;    }
               if ( g2count > 0 )
               {    if ( last_length > 0 )
	            {    gaps.Prepend( -g2count );
                         lengths.Prepend( last_length );    }
                    g2count = 0;    }
               ++lcount;
               --i;
               --j;   }
          else if ( from == &z_from )  // gap on long sequence
          {    if ( lcount > 0 )
               {    last_length = lcount;
                    lcount = 0;    }
               ForceAssert( g1count == 0 );
               ++g2count;
               --i;    }
          else                           // gap on short sequence
          {    if ( lcount > 0 )
               {    last_length = lcount;
                    lcount = 0;    }
               ForceAssert( g2count == 0 );
               ++g1count;
               --j;    }

	  if ( dir == 'x') from = &x_from;
	  else if ( dir == 'y') from = &y_from;
	  else from = &z_from;

	  if( (*from)[i][j] == 's' ) break;

    }

     if ( g1count != 0 ) gaps.Prepend( g1count );
     else if ( g2count != 0 ) gaps.Prepend( -g2count );
     else gaps.Prepend(0);

     lengths.Prepend( lcount );

     int pos1 = i;
     int pos2 = j;

     if ( gaps(0) < 0 )
     {   pos2 -= gaps(0);
         gaps(0) = 0;    }

     if ( gaps(0) > 0 )
     {   pos1 += gaps(0);
         gaps(0) = 0;    }

     *a = align( pos1, pos2, gaps, lengths );
     return best_score;
}

unsigned int SmithWatAffineBandedCoreFast( const basevector& S, const basevector& T,
			     align* a, int offset, int bandwidth,
                             bool penalize_left_gap, bool penalize_right_gap,
                             const int mismatch_penalty,
                             const int gap_open_penalty,
                             const int gap_extend_penalty )
{
     ForceAssertGt( S.size(), 0u );
     ForceAssertGt( T.size(), 0u );

     int n = S.size( ), N = T.size( );

     avector<char> s, t;
     s.resize(n);
     for (int i = 0; i < n; i++)
          s(i) = S[i];
     t.resize(N);
     for (int i = 0; i < N; i++)
          t(i) = T[i];

     int best_score = Infinity;
     typedef BandedArray<unsigned int, Infinity> ScoringArrayT;
     ScoringArrayT score_x(n+1, N+1, offset, bandwidth);
     ScoringArrayT score_y(n+1, N+1, offset, bandwidth);
     ScoringArrayT score_z(n+1, N+1, offset, bandwidth);

     typedef BandedArray<unsigned char, 's'> DirArrayT;
     DirArrayT x_from(n+1, N+1, offset, bandwidth);
     DirArrayT y_from(n+1, N+1, offset, bandwidth);
     DirArrayT z_from(n+1, N+1, offset, bandwidth);

     // y and z are for the gap extending in target and source, respectively
     // x is for base substitution
     //score_x[0][0] = 0;
     //score_y[0][0] = Infinity;
     //score_z[0][0] = Infinity;
     // valid only if (j >= i - offset - bandwidth && j <= i - offset + bandwidth)
     if ( 0 >= 0 - offset - bandwidth && 0 <= 0 - offset + bandwidth)
         score_x.Mutable(0,0) = 0;
     for (int i = 1; i <= n; i++ ) {
         if ( 0 < i - offset - bandwidth || 0  > i - offset + bandwidth)
             continue;
         score_z.Mutable(i,0) = 0;
     }
     if (!penalize_left_gap) {
         for (int j = 1; j <= N; j++) {
             if (j < 0 - offset - bandwidth || j > 0 - offset + bandwidth)
                 continue;
             score_y.Mutable(0,j) = 0;
         }
     }
     double time = WallClockTime();
     for (int i = 1; i <= n; i++ )
     {
          const int left = i - offset - bandwidth;
          const int right = i - offset + bandwidth;
          int low = Max( 1, left );
          int high = Min( N, right );
          for (int j = low; j <= high; j++ )
	  {
               //if (j < i - offset - bandwidth || j > i - offset + bandwidth) continue;
               // if (j < left || j > right) continue;
               int mismatch_score = (s(i-1) == t(j-1) ? 0 : mismatch_penalty);
               unsigned int x_x = score_x[i-1][j-1] + mismatch_score;
	       unsigned int x_y = score_y[i-1][j-1] + mismatch_score;
	       unsigned int x_z = score_z[i-1][j-1] + mismatch_score;
	       unsigned int y_x = score_x[i][j-1] + (penalize_right_gap || i != n ? gap_open_penalty : 0);
	       unsigned int y_y = score_y[i][j-1] + (penalize_right_gap || i != n ? gap_extend_penalty : 0);
	       unsigned int y_z = Infinity; //score_z[i][j-1] + gap_open_penalty;
	       unsigned int z_x = score_x[i-1][j] + (j != N ? gap_open_penalty : 0);
	       unsigned int z_y = Infinity; //score_y[i-1][j] + gap_open_penalty;
	       unsigned int z_z = score_z[i-1][j] + (j != N ? gap_extend_penalty: 0);

	       score_x.Mutable(i,j) = Min( Min( x_x, x_y ), x_z );
	       score_y.Mutable(i,j) = Min( Min( y_x, y_y ), y_z );
	       score_z.Mutable(i,j) = Min( Min( z_x, z_y ), z_z );
               if ((int)score_x[i][j] > Infinity) {
                   cout << "Renormalized at " << " i= " << i << " j= " << j
                      << " score_x[i][j]= " << score_x[i][j] << endl;
                   score_x.Mutable(i,j) = Infinity;
               }
	       if ( x_x <= x_y )
	       {    if ( x_x <= x_z ) x_from.Mutable(i,j) = 'x';
	            else x_from.Mutable(i,j) =  'z';    }
	       else
	       {    if ( x_y <= x_z ) x_from.Mutable(i,j) = 'y';
	            else x_from.Mutable(i,j) =  'z';    }

	       if ( y_x <= y_y )
	       {    if ( y_x <= y_z ) y_from.Mutable(i,j) = 'x';
	            else y_from.Mutable(i,j) =  'z';    }
	       else
	       {    if ( y_y <= y_z ) y_from.Mutable(i,j) = 'y';
	            else y_from.Mutable(i,j) =  'z';    }

	       if ( z_x <= z_y )
	       {    if ( z_x <= z_z ) z_from.Mutable(i,j) = 'x';
	            else z_from.Mutable(i,j) =  'z';    }
	       else
	       {    if ( z_y <= z_z ) z_from.Mutable(i,j) = 'y';
	            else z_from.Mutable(i,j) =  'z';    }    }    }

     best_score = Min( score_x[n][N], Min( score_y[n][N], score_z[n][N] ) );

     // Find the end point of the best matching trace
     int ii = n, jj = N;
     for(int k = n; k >= 0; k--){
         int best_score_k = Min( score_x[k][N], Min( score_y[k][N], score_z[k][N] ) );
         if (best_score_k < best_score) {
             best_score = best_score_k;
             ii = k;
             jj = N;
         }
     }
     if (!penalize_right_gap) {
     for(int k = N; k >= 0; k--){
         int best_score_k = Min( score_x[n][k], Min( score_y[n][k], score_z[n][k] ) );
         if (best_score_k < best_score) {
             best_score = best_score_k;
             ii = n;
             jj = k;
         }
     }
     }
     if (best_score == Infinity) return best_score;

     DirArrayT *from;
     if ( score_x[ii][jj] <= score_y[ii][jj] )
     {    if ( score_x[ii][jj] <= score_z[ii][jj] ) from = &x_from;
          else from = &z_from;    }
     else
     {    if ( score_y[ii][jj] <= score_z[ii][jj] ) from = &y_from;
          else from = &z_from;    }
     int i = ii;
     int j = jj;

     int lcount = 0, g1count = 0, g2count = 0;
     int last_length = 0;
     avector<int> gaps(0), lengths(0);
     while(1)
     {
          unsigned char dir = (*from)[i][j];
	  //cout << dir;
          if ( from == &x_from )
          {    if ( g1count > 0 )
               {    if ( last_length > 0 )
		    {    gaps.Prepend( g1count );
		         lengths.Prepend( last_length );    }
	            g1count = 0;    }
               if ( g2count > 0 )
               {    if ( last_length > 0 )
	            {    gaps.Prepend( -g2count );
                         lengths.Prepend( last_length );    }
                    g2count = 0;    }
               ++lcount;
               --i;
               --j;   }
          else if ( from == &z_from )  // gap on long sequence
          {    if ( lcount > 0 )
               {    last_length = lcount;
                    lcount = 0;    }
               ForceAssert( g1count == 0 );
               ++g2count;
               --i;    }
          else                           // gap on short sequence
          {    if ( lcount > 0 )
               {    last_length = lcount;
                    lcount = 0;    }
               ForceAssert( g2count == 0 );
               ++g1count;
               --j;    }

	  if ( dir == 'x') from = &x_from;
	  else if ( dir == 'y') from = &y_from;
	  else from = &z_from;

	  if( (*from)[i][j] == 's' ) break;

    }

     if ( g1count != 0 ) gaps.Prepend( g1count );
     else if ( g2count != 0 ) gaps.Prepend( -g2count );
     else gaps.Prepend(0);

     lengths.Prepend( lcount );

     int pos1 = i;
     int pos2 = j;

     if ( gaps(0) < 0 )
     {   pos2 -= gaps(0);
         gaps(0) = 0;    }

     if ( gaps(0) > 0 )
     {   pos1 += gaps(0);
         gaps(0) = 0;    }

     *a = align( pos1, pos2, gaps, lengths );
     return best_score;
}

} // end anonymous namespace

unsigned int SmithWatAffine( const basevector& S, const basevector& T,
			     alignment& a,
			     bool penalize_left_gap,
			     bool penalize_right_gap,
                             const int mismatch_penalty,
                             const int gap_open_penalty,
                             const int gap_extend_penalty )
{
     ForceAssertGt( S.size(), 0u );
     ForceAssertGt( T.size(), 0u );

     unsigned int n = S.size( ), N = T.size( );

     //     ForceAssertLe( n, N );

     avector<char> s, t;
     s.resize(n);
     for ( unsigned int i = 0; i < n; i++ )
          s(i) = S[i];
     t.resize(N);
     for ( unsigned int i = 0; i < N; i++ )
          t(i) = T[i];

     int best_score = Infinity;
     vec< vec<unsigned int> > score_x;
     vec< vec<unsigned int> > score_y;
     vec< vec<unsigned int> > score_z;

     vec< vec<unsigned char> > x_from;
     vec< vec<unsigned char> > y_from;
     vec< vec<unsigned char> > z_from;

     score_x.resize( n+1 );
     score_y.resize( n+1 );
     score_z.resize( n+1 );
     x_from.resize( n+1 );
     y_from.resize( n+1 );
     z_from.resize( n+1 );

     for ( unsigned int i = 0; i <= n; ++i )
     {
         score_x[i].resize( N+1 );
         score_y[i].resize( N+1 );
         score_z[i].resize( N+1 );
         x_from[i].resize( N+1 );
         y_from[i].resize( N+1 );
         z_from[i].resize( N+1 );
     }

     score_x[0][0] = 0;
     score_y[0][0] = Infinity;
     score_z[0][0] = Infinity;
     x_from[0][0] = 's';
     y_from[0][0] = 's';
     z_from[0][0] = 's';

     for ( unsigned int i = 1; i <= n; i++ )
     {    score_x[i][0] = Infinity;
	  score_y[i][0] = Infinity;
          score_z[i][0] = gap_open_penalty + gap_extend_penalty * (i-1);
	  x_from[i][0] = 's';
	  y_from[i][0] = 's';
	  z_from[i][0] = 's';   }

     for ( unsigned int j = 1; j <= N; j++)
       {  score_x[0][j] = Infinity;
          score_y[0][j] = (penalize_left_gap ? gap_open_penalty + gap_extend_penalty * (j-1) : 0);
	  score_z[0][j] = Infinity;
	  x_from[0][j] = 's';
	  y_from[0][j] = 's';
	  z_from[0][j] = 's';   }

     for ( unsigned int i = 1; i <= n; i++ )
     {   for ( unsigned int j = 1; j <= N; j++ )
	  {    unsigned int x_x = score_x[i-1][j-1] + mismatch_penalty * ( s(i-1) != t(j-1) );
	       unsigned int x_y = score_y[i-1][j-1] + mismatch_penalty * ( s(i-1) != t(j-1) );
	       unsigned int x_z = score_z[i-1][j-1] + mismatch_penalty * ( s(i-1) != t(j-1) );
	       unsigned int y_x = score_x[i][j-1] + (i != n || penalize_right_gap ? gap_open_penalty : 0);
	       unsigned int y_y = score_y[i][j-1] + (i != n || penalize_right_gap ? gap_extend_penalty : 0);
	       unsigned int y_z = Infinity; //score_z[i][j-1] + gap_open_penalty;
	       unsigned int z_x = score_x[i-1][j] + gap_open_penalty;
	       unsigned int z_y = Infinity; //score_y[i-1][j] + gap_open_penalty;
	       unsigned int z_z = score_z[i-1][j] + gap_extend_penalty;

	       score_x[i][j] = Min( Min( x_x, x_y ), x_z );
	       score_y[i][j] = Min( Min( y_x, y_y ), y_z );
	       score_z[i][j] = Min( Min( z_x, z_y ), z_z );

	       if ( x_x <= x_y )
	       {    if ( x_x <= x_z ) x_from[i][j] = 'x';
	            else x_from[i][j] =  'z';    }
	       else
	       {    if ( x_y <= x_z ) x_from[i][j] = 'y';
	            else x_from[i][j] =  'z';    }

	       if ( y_x <= y_y )
	       {    if ( y_x <= y_z ) y_from[i][j] = 'x';
	            else y_from[i][j] =  'z';    }
	       else
	       {    if ( y_y <= y_z ) y_from[i][j] = 'y';
	            else y_from[i][j] =  'z';    }

	       if ( z_x <= z_y )
	       {    if ( z_x <= z_z ) z_from[i][j] = 'x';
	            else z_from[i][j] =  'z';    }
	       else
	       {    if ( z_y <= z_z ) z_from[i][j] = 'y';
	            else z_from[i][j] =  'z';    }    }    }

     best_score = Min( score_x[n][N], Min( score_y[n][N], score_z[n][N] ) );
     int ii = n;
     int jj = N;

#ifdef  FIX_RIGHT_GAP
    int right = Min(n,N);
    if (!penalize_right_gap) {
        for(int k = N; k >= right; k--){
            int best_score_k = Min( score_x[n][k], Min( score_y[n][k], score_z[n][k] ) );
            if (best_score_k < best_score) {
                best_score = best_score_k;
                ii = n;
                jj = k;
            }
        }
    }
#endif

     vec< vec<unsigned char> > *from;
     if ( score_x[n][N] <= score_y[n][N] )
     {    if ( score_x[n][N] <= score_z[n][N] ) from = &x_from;
          else from = &z_from;    }
     else
     {    if ( score_y[n][N] <= score_z[n][N] ) from = &y_from;
          else from = &z_from;    }

     int i = ii;
     int j = jj;
     int lcount = 0, g1count = 0, g2count = 0;
     int last_length = 0;
     avector<int> gaps(0), lengths(0);
     while(1)
     {
          unsigned char dir = (*from)[i][j];
	  //cout << dir;
          if ( from == &x_from )
          {    if ( g1count > 0 )
               {    if ( last_length > 0 )
		    {    gaps.Prepend( g1count );
		         lengths.Prepend( last_length );    }
	            g1count = 0;    }
               if ( g2count > 0 )
               {    if ( last_length > 0 )
	            {    gaps.Prepend( -g2count );
                         lengths.Prepend( last_length );    }
                    g2count = 0;    }
               ++lcount;
               --i;
               --j;   }
          else if ( from == &z_from )  // gap on long sequence
          {    if ( lcount > 0 )
               {    last_length = lcount;
                    lcount = 0;    }
               ForceAssert( g1count == 0 );
               ++g2count;
               --i;    }
          else                           // gap on short sequence
          {    if ( lcount > 0 )
               {    last_length = lcount;
                    lcount = 0;    }
               ForceAssert( g2count == 0 );
               ++g1count;
               --j;    }

	  if ( dir == 'x') from = &x_from;
	  else if ( dir == 'y') from = &y_from;
	  else from = &z_from;

	  if( (*from)[i][j] == 's' ) break;

    }

     //cout << "\n";

     if ( g1count != 0 ) gaps.Prepend( g1count );
     else if ( g2count != 0 ) gaps.Prepend( -g2count );
     else gaps.Prepend(0);

     lengths.Prepend( lcount );

     int pos1 = i;
     int pos2 = j;

     if ( gaps(0) < 0 )
     {   pos2 -= gaps(0);
         gaps(0) = 0;    }

     if ( gaps(0) > 0 )
     {   pos1 += gaps(0);
         gaps(0) = 0;    }

     int errors = best_score;
     a = alignment( pos1, pos2, errors, gaps, lengths );

     return best_score;    }


unsigned int SmithWatAffineParallel(const basevector & S, const basevector & T,
                            alignment & a,
                            bool penalize_left_gap,
                            bool penalize_right_gap,
                            const int mismatch_penalty,
                            const int gap_open_penalty,
                            const int gap_extend_penalty
                            )
{

    ForceAssertGt(S.size(), 0u);
    ForceAssertGt(T.size(), 0u);

    unsigned int n = S.size(), N = T.size();

    //     ForceAssertLe( n, N );

    avector < char >s, t;
    s.resize(n);
    for (unsigned int i = 0; i < n; i++)
        s(i) = S[i];
    t.resize(N);
    for (unsigned int i = 0; i < N; i++)
        t(i) = T[i];

    const int Infinity = 100000000;
    int best_score = Infinity;

    RecArray<int> score_x( n+1, N+1 );
    RecArray<int> score_y( n+1, N+1 );
    RecArray<int> score_z( n+1, N+1 );

    RecArray<unsigned char> x_from( n+1, N+1 );
    RecArray<unsigned char> y_from( n+1, N+1 );
    RecArray<unsigned char> z_from( n+1, N+1 );

    score_x[0][0] = 0;
    score_y[0][0] = Infinity;
    score_z[0][0] = Infinity;
    x_from[0][0] = 's';
    y_from[0][0] = 's';
    z_from[0][0] = 's';

    for (unsigned int i = 1; i <= n; i++) {
        score_x[i][0] = Infinity;
        score_y[i][0] = Infinity;
        score_z[i][0] = gap_open_penalty + gap_extend_penalty * i;
        x_from[i][0] = 's';
        y_from[i][0] = 's';
        z_from[i][0] = 's';
    }

    // Zero penalty for unmatched bases on left side of T
    for (unsigned int j = 1; j <= N; j++) {
        score_x[0][j] = Infinity;
        score_y[0][j] = (penalize_left_gap ? gap_open_penalty +
                gap_extend_penalty * j :0);
        score_z[0][j] = Infinity;
        x_from[0][j] = 's';
        y_from[0][j] = 's';
        z_from[0][j] = 's';
    }

    // parallelize the smith-waterman calculation.
    // The idea is to divide the n by N matrix (score[1..n][1..N]) into many (BxB) blocks.
    // Every block will only depend on the block on it's left, up, and up-left
    // blocks, and all the block(i,j) with same i + j can be calculated parallelized.
    //
    // block size = 200 is chosen heuristically, which seems most efficient
    //
    const int block_size = 200;
    int n_block_s = n / block_size;
    if ( n_block_s * block_size < (int)n ) n_block_s += 1;
    int n_block_t = N / block_size;
    if ( n_block_t * block_size < (int)N ) n_block_t += 1;
    //cout <<  Date( ) << ": start the main loop" << endl;
    //cout <<  Date( ) << ": time used = " << TimeSince(clock) << endl;
    for ( int level = 0; level < n_block_s + n_block_t; ++level ) {
    // number of blocks in the diagonal that are independent to each other
    // start and stop of the block index
    // i_block_s in [0, n_block_s)
    // i_block_t in [0, n_block_t)
    // where i_block_s + i_block_t = level

    // in other words i_block_s is in ( level - n_block_t, level ] && [0, n_block_s)

    #pragma omp parallel for
    for ( int i_block_s = max( 0, level - n_block_t + 1 );
          i_block_s < min( n_block_s, level + 1 ); i_block_s++ ) {
        int i_block_t = level - i_block_s;
        ForceAssertGe( i_block_t , 0 );
        ForceAssertLt( i_block_t , n_block_t );

        int istart = i_block_s * block_size + 1;
        int istop = min( istart + block_size, (int)n + 1 );
        int jstart = i_block_t * block_size + 1;
        int jstop = min( jstart + block_size, (int)N + 1 );


    for ( int i = istart; i < istop; i++) {
        for ( int j = jstart; j < jstop; j++) {
            unsigned int x_x =
                score_x[i - 1][j - 1] + mismatch_penalty * (s(i - 1) !=
                                                            t(j - 1));
            unsigned int x_y =
                score_y[i - 1][j - 1] + mismatch_penalty * (s(i - 1) !=
                                                            t(j - 1));
            unsigned int x_z =
                score_z[i - 1][j - 1] + mismatch_penalty * (s(i - 1) !=
                                                            t(j - 1));
            //unsigned int y_x = score_x[i][j - 1] + gap_open_penalty;
            //unsigned int y_y = score_y[i][j - 1] + gap_extend_penalty;
            unsigned int y_x = score_x[i][j-1] +
                (penalize_right_gap || i != (int)n ? gap_open_penalty : 0);
            unsigned int y_y = score_y[i][j-1] +
                (penalize_right_gap || i != (int)n ? gap_extend_penalty : 0);

            unsigned int y_z = Infinity;        //score_z[i][j-1] + gap_open_penalty;
            unsigned int z_x = score_x[i - 1][j] + gap_open_penalty;
            unsigned int z_y = Infinity;        //score_y[i-1][j] + gap_open_penalty;
            unsigned int z_z = score_z[i - 1][j] + gap_extend_penalty;

            score_x[i][j] = Min(Min(x_x, x_y), x_z);
            score_y[i][j] = Min(Min(y_x, y_y), y_z);
            score_z[i][j] = Min(Min(z_x, z_y), z_z);

            if (x_x <= x_y) {
                if (x_x <= x_z)
                    x_from[i][j] = 'x';
                else
                    x_from[i][j] = 'z';
            }
            else {
                if (x_y <= x_z)
                    x_from[i][j] = 'y';
                else
                    x_from[i][j] = 'z';
            }

            if (y_x <= y_y) {
                if (y_x <= y_z)
                    y_from[i][j] = 'x';
                else
                    y_from[i][j] = 'z';
            }
            else {
                if (y_y <= y_z)
                    y_from[i][j] = 'y';
                else
                    y_from[i][j] = 'z';
            }

            if (z_x <= z_y) {
                if (z_x <= z_z)
                    z_from[i][j] = 'x';
                else
                    z_from[i][j] = 'z';
            }
            else {
                if (z_y <= z_z)
                    z_from[i][j] = 'y';
                else
                    z_from[i][j] = 'z';
            }
        }
    }



    // end parallelization
    }
    }

    int ii = n, jj = N;
    best_score = Min(score_x[n][N], Min(score_y[n][N], score_z[n][N]));
#ifdef FIX_RIGHT_GAP
    int right = Min(n,N);
#else
    int right = 0;
#endif
    if (!penalize_right_gap) {
        for(int k = N; k >= right; k--){
            int best_score_k = Min( score_x[n][k], Min( score_y[n][k], score_z[n][k] ) );
            if (best_score_k < best_score) {
                best_score = best_score_k;
                ii = n;
                jj = k;
            }
        }
    }

    //vec < vec < unsigned char > >*from;
    RecArray<unsigned char> *from;
    if (score_x[ii][jj] <= score_y[ii][jj]) {
        if (score_x[ii][jj] <= score_z[ii][jj])
            from = &x_from;
        else
            from = &z_from;
    }
    else {
        if (score_y[ii][jj] <= score_z[ii][jj])
            from = &y_from;
        else
            from = &z_from;
    }

    int i = ii;
    int j = jj;
    int lcount = 0, g1count = 0, g2count = 0;
    int last_length = 0;
    avector < int >gaps(0), lengths(0);
    while (1) {
        unsigned char dir = (*from)[i][j];
        //cout << dir;
        if (from == &x_from) {
            if (g1count > 0) {
                if (last_length > 0) {
                    gaps.Prepend(g1count);
                    lengths.Prepend(last_length);
                }
                g1count = 0;
            }
            if (g2count > 0) {
                if (last_length > 0) {
                    gaps.Prepend(-g2count);
                    lengths.Prepend(last_length);
                }
                g2count = 0;
            }
            ++lcount;
            --i;
            --j;
        }
        else if (from == &z_from)       // gap on long sequence
        {
            if (lcount > 0) {
                last_length = lcount;
                lcount = 0;
            }
            ForceAssert(g1count == 0);
            ++g2count;
            --i;
        }
        else                    // gap on short sequence
        {
            if (lcount > 0) {
                last_length = lcount;
                lcount = 0;
            }
            ForceAssert(g2count == 0);
            ++g1count;
            --j;
        }

        if (dir == 'x')
            from = &x_from;
        else if (dir == 'y')
            from = &y_from;
        else
            from = &z_from;

        if ((*from)[i][j] == 's')
            break;

    }

    //cout << "\n";

    if (g1count != 0)
        gaps.Prepend(g1count);
    else if (g2count != 0)
        gaps.Prepend(-g2count);
    else
        gaps.Prepend(0);

    lengths.Prepend(lcount);

    int pos1 = i;
    int pos2 = j;

    if (gaps(0) < 0) {
        pos2 -= gaps(0);
        gaps(0) = 0;
    }

    if (gaps(0) > 0) {
        pos1 += gaps(0);
        gaps(0) = 0;
    }

    int errors = best_score;
    a = alignment(pos1, pos2, errors, gaps, lengths);

    return best_score;
}

// Perform a banded alignment of sequence S to T, given that the relative
// offset between S and T.
//
// Offset is the relative shift of target with respect to source sequence
// Unmatched ends on both S and T are not penalized. Meaningful alignments
// outputs are guaranteed only if reasonable offset and bandwidth are given.
unsigned int SmithWatAffineBanded( const basevector& S, const basevector& T,
                             int offset, int bandwidth,
			     align& a, int& error,
                             const int mismatch_penalty,
                             const int gap_open_penalty,
                             const int gap_extend_penalty )
{
    ForceAssertGt(S.isize(), 0);
    ForceAssertGt(T.isize(), 0);
    ForceAssertGt(bandwidth, 0);

    int t_start = Max(0, -offset - bandwidth);
    int t_stop = Min(T.isize(), S.isize() - offset + bandwidth);
    if (t_start >= t_stop) {
        a = alignment();
        return Infinity;
    }
    basevector T2(T, t_start, t_stop - t_start);
    // new offset = s_pos - (t_pos - t_start) = offset + t_start
    int score = SmithWatAffineBandedCoreFast(S, T2, &a, t_start + offset, bandwidth,
            false, false, mismatch_penalty, gap_open_penalty, gap_extend_penalty);
    error = a.Errors(S, T2);
    a.Setpos2(a.pos2()  + t_start);
    if ( a.pos1() < 0 || a.pos2() < 0 || a.Pos1() > (int)S.size()
            || a.Pos2() > (int)T.size() ) {
        a = align();
        score = Infinity;
    }
    return score;
}


// SmithWatAffineSuper - long edge alignment
//
// The algorithm kmerizes the reads, S and T, and then calls
// SmithWatAffine on the sequences between perfect matches.  Perfect
// matches are defined as runs of uniquely agreeing (and properly
// ordered) kmers between the two sequences.
//
// two things could be improved:
// 1. we currently extend via kmerization, which at the default K=501 is
// more expensive than it needs to be -- we could just count matching
// bases, given a kmer seed.
//
// 2. we break a kmer run at the point of a disagreement, which puts the
// disagreement at the edge of the sequence which gets passed to
// Smith-Waterman.  To put it another way: we only pad the sequence with
// matching bases on one side of the event.
//

struct perf_interval {
    size_t pos1;
    size_t pos2;
    size_t len;

    struct sort_by_pos1 {
        bool operator() ( perf_interval const& a, perf_interval const& b ) const {
            if ( a.pos1 < b.pos1 ) return true;
            else if ( a.pos1 == b.pos1 && a.pos2 < b.pos2 ) return true;
            else return false;
        }
    };

    struct sort_by_pos2 {
        bool operator() ( perf_interval const& a, perf_interval const& b ) const {
            if  ( a.pos2 < b.pos2 ) return true;
            else if ( a.pos2 == b.pos2 && a.pos1 < b.pos1 ) return true;
            else return false;
        }
    };
};


static std::ostream& operator<< (std::ostream &os, const perf_interval& p) {
    os << "pos1=" << p.pos1 << ", pos2=" << p.pos2 << ", len=" << p.len;
    return os;
};

// filter_perfect_intervals, score_perfect_intervals -- by design, these intervals are
// non-overlapping and non-crossing with respect to pos2, however pos1
// may cross over.  This code removes that.

static size_t score_perfect_intervals( vec<perf_interval>& perf_intervals, size_t skip = 0 )
{
    ForceAssertGe( perf_intervals.size(), 1U );

    size_t badness = 0;

    perf_interval last = perf_intervals[0];
    for (size_t i = 1; i < perf_intervals.size(); ++i ) {
        if ( i == skip ) continue;
        perf_interval curr = perf_intervals[i];
        if ( curr.pos1 < last.pos1 + last.len || curr.pos2 < last.pos2 + last.len )
            badness++;

        last = curr;
    }

    return badness;
}

static void filter_perfect_intervals( vec<perf_interval>& input, int verbose = 0 )
{
    ForceAssertGe( input.size(), 1U );

    perf_interval last = input[0];
    vector<perf_interval> output;

    // okay, so we assume that the endpoints are not removable...
    for ( size_t i = 1; i < input.size(); ++i ) {
        auto const& curr = input[i];

        if ( curr.pos1 < last.pos1 + last.len || curr.pos2 < last.pos2 + last.len ) {
            // okay, we have a conflict between two adjacent intervals,
            // the current and the last.  let's decide which to drop
            bool drop_last = false;     // this handles the i==1 case, drop curr

            if ( i == input.size()-1 ) {  // can't drop curr, so we drop last
                if (verbose >= 2 ) cout << "case 1" << endl;
                drop_last = true;
            } else if ( i > 1 ) {         // in the middle, decide which to drop
                size_t score_drop_last = score_perfect_intervals( input, i-1 );
                size_t score_drop_curr = score_perfect_intervals( input, i );
                drop_last = score_drop_last < score_drop_curr;
                if ( verbose >= 2 ) {
                    cout << "case 2" << endl;
                    cout << "last=" << last;
                    cout << ", curr=" << curr << endl;
                    PRINT3( score_drop_last, score_drop_curr, drop_last );
                }
            } else
               if ( verbose >= 2 ) cout << "case 3" << endl;

            // okay, so if we're dropping 'last', we don't push back and
            // we make curr the new last.  If we're dropping 'curr',
            // then we keep 'last' for the next iteration.
            if ( drop_last ) {
                last=curr;
            }

        } else {        // no conflict, so we push back the last one
            output.push_back( last );
            last = curr;
        }
    }
    output.push_back( input.back() );

    std::swap(input, output);
}

// okay, try to drop seeds that imply a massive gap
// run AFTER removing non-affine
static void filter_perfect_intervals3( vec<perf_interval>& input, bool penalize_left_gap, bool penalize_right_gap )
{
    vec<Bool> delete_if( input.size(), False );

    size_t start = penalize_left_gap ? 1 : 2;
    size_t end = input.size() - 2;      // we're looking left, so the right gap doesn't matter

    size_t j = start - 1;
    size_t pad = 10U;
    for ( size_t i = start; i <= end; i++ ) {
        ForceAssertGe( input[i].pos1, input[j].pos1 );
        ForceAssertGe( input[i].pos2, input[j].pos2 );
        size_t len1 = max( input[i].pos1 - input[j].pos1, pad );
        size_t len2 = max( input[i].pos2 - input[j].pos2, pad );

        if ( max( len1, len2 ) > 500 &&
                ( static_cast<double>(len1)/len2 > (3./2.) ||
                  static_cast<double>(len1)/len2 < (2./3.) ) ) {
            delete_if[ i ] = True;
        } else
            j=i;
    }

    EraseIf( input, delete_if );
}

// toss non-affine relationships
static void filter_perfect_intervals2( vec<perf_interval>& input )
{
    // so the intervals come sorted by pos1 by design
    ForceAssert( std::is_sorted(
              input.begin(), input.end(),
              perf_interval::sort_by_pos1() ) );
    // let's sort by pos2 and see what happens
    vec<perf_interval> input_tmp( input );
    vec<size_t> ident( input_tmp.size(), vec<size_t>::IDENTITY );
    SortSync( input_tmp, ident, perf_interval::sort_by_pos2() );

    vec<size_t> rejects_pos1;
    size_t j = 0;
    size_t i = 0;
    while ( j < input.size() ) {
        if ( ident[i] == j ) {
            i++;
            j++;
        } else {
            if ( Member( rejects_pos1, j ) ) {
                j++;
            } else {
                rejects_pos1.push_back( ident[i] );
                i++;
            }
        }
    }

    vec<Bool> delete_if( input.size(), False);
    for ( size_t pos1 : rejects_pos1 )
        delete_if[pos1] = True;

    EraseIf( input, delete_if );
}


static void normalize_gaps_lengths( size_t& pos1, size_t& pos2, avector<int>& gaps, avector<int>& lengths )
{
    ForceAssertEq(gaps.length, lengths.length);

    std::vector<int> new_gaps, new_lengths;

    // convert initial gap to an offset, only if it doesn't make a
    // conflicting offset
    if ( gaps.length > 0 ) {
        if ( gaps(0) > 0 && pos1 == 0 ) {
            pos2 += gaps(0);
            gaps(0)=0;
        } else if ( gaps(0) < 0 && pos2 == 0 ) {
            pos1 -= gaps(0);
            gaps(0)=0;
        }
    }

    for (size_t i = 0; i < gaps.length-1; ++i ) {
        if ( gaps(i) == 0 && lengths(i) == 0 ) {
            // skip no-op
            continue;
        } else if ( gaps(i+1) == 0 ) {
            // forward propagate
            lengths(i+1) += lengths(i);
            gaps(i+1) += gaps(i);
        } else {
            // emit
            new_gaps.push_back( gaps(i) );
            new_lengths.push_back( lengths(i) );
        }
    }

    // only add final entry if there's a length and not just a final gap
    if ( lengths(lengths.length-1) > 0 ) {
        new_gaps.push_back( gaps(gaps.length-1) );
        new_lengths.push_back( lengths(lengths.length-1) );
    }


    gaps.Setsize( new_gaps.size() );
    lengths.Setsize( new_lengths.size() );

    ForceAssertEq( new_gaps.size(), new_lengths.size() );
    for ( size_t i = 0; i < new_gaps.size(); ++i ) {
        gaps(i) = new_gaps[i];
        lengths(i) = new_lengths[i];
    }
}


static void find_perfect_seeds( const basevector& S, const basevector& T, int K,
        bool penalize_left_gap, bool penalize_right_gap, vec<perf_interval>& perf_intervals, int verbose = 0 )
{
    // seeding threshold is the minimum size problem (in bases (S.size()*T.size())
    // for which we will seed.  For smaller problems, say a 4k edge versus a 40k
    // reference, we will pass on to the "subaligner".  This is generally just
    // falling back to the usual aligner
    const size_t seeding_threshold = 40000*5000;

    if ( size_t(S.size())*size_t(T.size()) < seeding_threshold ) {   // non-structured short-circuit for short edges
        perf_intervals.clear();
        perf_intervals.push_back({0,0,0});
        perf_intervals.push_back({S.size(),T.size(),0});
        return;
    }

    ForceAssertGt(min(S.size(), T.size()), static_cast<size_t>(K)+10u);

    vec<perf_interval> new_perf_intervals;

    size_t NUM_THREADS = 1;		// boundNumThreads(0);

    // kmerize the "reference"
    ForceAssertLt(K,504);
    typedef Kmer504 Kmer_t;
    typedef KmerKmerBVLoc<Kmer_t> KmerRec_t;
    vec<KmerRec_t> kmer_vec;


    vecbasevector Tv(1,T);

    // calculate the kmer db

    KernelKmerBVLocStorer<KmerRec_t> storer( Tv, K, &kmer_vec);
    naif_kmerize( &storer, NUM_THREADS, 0 );

    // kmer map hash into the kmer db
    KmerMap<KmerRec_t> kmer_map( kmer_vec );

    typedef KmerKmerFreq<Kmer_t> KmerFreqRec_t;
    vec<KmerFreqRec_t> kmer_freq_vec;

#if 0
    vecbasevector Sv(1,S);
    KernelKmerStorer<KmerFreqRec_t> freq_storer(Sv, K, &kmer_freq_vec);
    naif_kmerize( &freq_storer, NUM_THREADS, 0 );

    KmerMap<KmerFreqRec_t> kmer_freq_map( kmer_freq_vec );
#endif

    // process S kmer-by-kmer:
    // at the start of a match we record the start position, and at the
    // end of the match we add the perfect interval to the list.
    SubKmers<basevector, Kmer_t> subs( K, S );
    size_t match_start_pos2 = 0U;
    size_t match_start_pos1 = 0U;
    size_t match_end_pos1 = 0U;
    bool in_match = false;

    new_perf_intervals.clear();
    new_perf_intervals.push_back({0,0,0});          // left edge is a zero-length, perfect interval

    if ( verbose >= 2 )
        cout << "SmithWatAffineSuper: align len " << S.size() << " to len " << T.size() << endl;

    while ( 1 ) {

        KmerRec_t& hit = kmer_map( subs.canonical() );

//        KmerFreqRec_t& hit2 = kmer_freq_map( subs.canonical() );
//        ForceAssert(hit2.is_valid_kmer() );

        // first establish whether this is a valid, unique hit in the
        // reference:
        // 1. valid kmer
        // 2. unique hit
        // 3. correct direction
        // 4. is in the position we expect
        bool valid = \
//            hit2.freq() == 1 &&
            subs.not_done() &&
            hit.is_valid_kmer() &&
            hit.locs().size() == 1 &&
            subs.is_canonical_fw() == hit.locs()[0].is_fw();

        if ( valid ) {
            // we shouldn't ever have this problem as the test for a
            // unique hit should cover it
            auto const& loc = hit.locs()[0];
            if ( in_match && loc.ib() != match_end_pos1 ) valid = false;
        }


        // now deal with it
        if ( in_match && !valid ) {

            // ENDING a perfect match here
            match_end_pos1 += K-1;
            in_match = false;
            if ( verbose >= 2 ) {
                cout << "perfect unique match from pos1=["
                         <<  match_start_pos1 << "," <<
                          match_end_pos1 << ") to " << endl;
                cout << "\tpos2=[" << match_start_pos2 << "," <<
                          subs.index_start()+K-1 << ")" << endl;
            }
            // push back result here
            if ( match_end_pos1 - match_start_pos1 > static_cast<unsigned int>(K) )
                new_perf_intervals.push_back( { match_start_pos2, match_start_pos1,
                    match_end_pos1 - match_start_pos1 } );

        } else if ( !in_match && valid ) {

            // STARTING a perfect match here
            auto const& loc = hit.locs()[0];
            match_start_pos2 = subs.index_start();
            match_start_pos1 = loc.ib();
            match_end_pos1 = match_start_pos1;
            in_match = true;
        }

        if ( !subs.not_done() )
            break;

        match_end_pos1++;
        subs.next();    // next kmer
    }

    // end with zero-length perfect interval
    new_perf_intervals.push_back( { S.size(), T.size(), 0 } );


    if ( verbose >= 2 ) {
        cout << "BEFORE backing off...." << endl;
        for (auto const& pi : new_perf_intervals )
            cout << "pos1=" << pi.pos1 << ", pos2=" << pi.pos2 << ", len=" << pi.len << endl;
    }

    // we trim back the perfect seeds to provide some context and breathing
    // room for adjancet alignments.  This is to avoid unintended consequences
    // e.g. we had a situation where it was better for a deletion in a homopolymer
    // to occur *within* the perfect seed, rather than next to it
    //
    // if the pefect seed runs right up to the edge of the sequence, on either
    // sequence, we do not trim there to avoid introducing an extra call to
    // the S-W aligner.
    perf_intervals.clear();
    const int trim_size = 15;
    for ( auto const& pi : new_perf_intervals ) {
	if ( pi.len == 0 )
	    perf_intervals.push_back( pi );
	else  {
	    size_t trim_left = trim_size;
	    size_t trim_right = trim_size;
	    if ( pi.pos1 == 0 || pi.pos2 == 0 ) trim_left = 0;
	    if ( pi.pos1 + pi.len == S.size() || pi.pos2 + pi.len == T.size() ) trim_right = 0;
	    if ( pi.len > trim_left + trim_right )
		perf_intervals.push_back( {pi.pos1 + trim_left, pi.pos2 + trim_left, pi.len - trim_left - trim_right});
	}
    }

    // okay, now remove any non-affine relationships implied by the
    // intervals.  This assumes that they've been created in-order for
    // pos2.
    if ( verbose >= 2 ) {
        cout << "BEFORE filtering...." << endl;
        for (auto const& pi : perf_intervals )
            cout << "pos1=" << pi.pos1 << ", pos2=" << pi.pos2 << ", len=" << pi.len << endl;
    }

    filter_perfect_intervals2( perf_intervals );

    if ( verbose >= 2 ) {
        cout << "AFTER filtering NON-AFFINE..." << endl;
        for (auto const& pi : perf_intervals )
            cout << "pos1=" << pi.pos1 << ", pos2=" << pi.pos2 << ", len=" << pi.len << endl;
    }

    filter_perfect_intervals3( perf_intervals, penalize_left_gap, penalize_right_gap );

    if ( verbose >= 2 ) {
        cout << "AFTER filtering large implied gaps ..." << endl;
        for (auto const& pi : perf_intervals )
            cout << "pos1=" << pi.pos1 << ", pos2=" << pi.pos2 << ", len=" << pi.len << endl;
    }

    filter_perfect_intervals( perf_intervals );

}


// convert an offset into an initial gap
// the point here is that we need to add a bunch of alignments together
// and while you can add gaps and lengths, there's only one pos1/pos2
// pair per alignment
static void convert_offset_initial_gap( alignment& al )
{
    avector<int> gaps, lengths;
    int pos1, pos2, errors;
    al.Unpack( pos1, pos2, errors, gaps, lengths );
    if ( pos1 == 0 )
    { if ( pos2 )
      { ForceAssertGt(gaps.length,0u);
        gaps(0) += pos2; } }
    else if ( pos2 == 0 )
    { if ( pos1 )
      { ForceAssertGt(gaps.length,0u);
        gaps(0) -= pos1; } }
    else
      FatalErr("neither pos1 nor pos2 is 0");
    al.Set( 0, 0, errors, gaps, lengths );
}


static void set_pos1pos2( alignment& al, int pos1, int pos2 )
{
    avector<int> gaps, lengths;
    int dummy1, dummy2, errors;
    al.Unpack( dummy1, dummy2, errors, gaps, lengths );
    al.Set( pos1, pos2, errors, gaps, lengths );
}


// add an initial gap to an existing alignment and score the penalty
static void alignment_add_gap( alignment& al, int gap, int gap_open_penalty, int gap_extend_penalty )
{
    avector<int> gaps, lengths;
    int pos1, pos2, errors;
    al.Unpack( pos1, pos2, errors, gaps, lengths );
    gaps(0) += gap;
    errors += gap_open_penalty + std::abs( gap-1 ) * gap_extend_penalty;
    al.Set( pos1, pos2, errors, gaps, lengths );
}


static void dumpit ( const alignment& a )
{
    avector<int> gaps, lengths;
    int errors;
    int pos1;
    int pos2;

    a.Unpack( pos1, pos2, errors, gaps, lengths );
    PRINT3( pos1, pos2, errors );
    cout << "\tGAPS:";
    for ( size_t i = 0; i < gaps.length; i++ )
        cout << " " << gaps(i);
    cout << endl << "\tLENS:";
    for ( size_t i = 0; i < lengths.length; i++ )
        cout << " " << lengths(i);
    cout << endl;
}



unsigned int SmithWatAffineSuper( const basevector& S, const basevector& T,
			     alignment& a, int K,
			     bool penalize_left_gap,
			     bool penalize_right_gap,
                             const int verbose,
                             const int mismatch_penalty,
                             const int gap_open_penalty,
                             const int gap_extend_penalty,
                             SubAlignFuncPtr subalign,
                             const int64_t max_product )
{
    float heuristic_size_min = 10.;        // require at least heuristic_size_min*K size for each

    ForceAssertGt(K, 0);
    if ( K % 2 == 0 ) K += 1;   // palindrome avoidance


    if ( verbose >= 1 ) {
        cout << Date() << ": SmithWatAffineSuper align S.len=" << S.size()  <<
            " T.len=" << T.size() << ", K=" << K << ", penalize L,R=" <<
            penalize_left_gap << "," << penalize_right_gap << endl;
    }

    // find seed intervals of extended perfect matches
    // and filter them for undesirable properties (e.g. non-affine
    // relationships)
    vec<perf_interval> perf_intervals;
    find_perfect_seeds( S, T, K, penalize_left_gap, penalize_right_gap, perf_intervals, verbose );

    // debugging
    if ( verbose >= 2 ) {
        cout << "AFTER filtering overlaps...." << endl;
        for (auto const& pi : perf_intervals )
            cout << "pos1=" << pi.pos1 << ", pos2=" << pi.pos2 << ", len=" << pi.len << endl;
    }

    // loop through intervals
    size_t new_pos1 = 0, new_pos2 = 0;

    vec<alignment> new_aligns;
    vec<int> extra_gaps;


    for ( size_t pi = 1; pi < perf_intervals.size(); pi++ ) {

	bool const final = ( pi == perf_intervals.size() - 1);

        auto const& last = perf_intervals[pi-1];
        auto const& cur = perf_intervals[pi];
        if ( verbose >= 2 ) {
            cout << "last " << pi-1 << ": pos1=" << last.pos1 << " pos2=" << last.pos2 <<
                " len=" << last.len << endl;
            cout << "cur " << pi << ": pos1=" << cur.pos1 << " pos2=" << cur.pos2 <<
                " len=" << cur.len << endl;
        }
        // set up subsequences S' & T' between the start of the current
        // interval and the end of the previous interval
        basevector Spr, Tpr;
        size_t start1 = last.pos1 + last.len;
        size_t len1 = cur.pos1 - last.pos1 - last.len;
        size_t start2 = last.pos2 + last.len;
        size_t len2 = cur.pos2 - last.pos2 - last.len;

        int next_gap = 0;

        if ( len1 > 0 || len2 > 0 ) {
            alignment tmp_align;
	    if (len1 > 0 && len2 > 0) {
		if (verbose >= 2)
		    cout << "SetToSubOf S: " << last.pos1 + last.len << " len "
			    << cur.pos1 - last.pos1 - last.len << endl;
		Spr.SetToSubOf(S, last.pos1 + last.len,
			cur.pos1 - last.pos1 - last.len);
		if (verbose >= 2)
		    cout << "SetToSubOf T: " << last.pos2 + last.len << " len "
			    << cur.pos2 - last.pos2 - last.len << endl;
		Tpr.SetToSubOf(T, last.pos2 + last.len,
			cur.pos2 - last.pos2 - last.len);

		// call SmithWatAffine on the subsequences
		bool t_left_gap = true, t_right_gap = true;
		if (pi == 1)
		    t_left_gap = penalize_left_gap;
		if (pi == perf_intervals.size() - 1)
		    t_right_gap = penalize_right_gap;
		if (verbose >= 2)
		    cout << "sw-align on ref [" << last.pos1 + last.len << ","
			    << cur.pos1 << ")  sizes " << Spr.size() << " to "
			    << Tpr.size() << endl << "L,R gap penalty="
			    << t_left_gap << "," << t_right_gap << endl;
                if ( (int64_t) Spr.size( ) * (int64_t) Tpr.size( ) > max_product )
                     return SmithWatAffineSuper_FAIL;
//		unsigned int tmp_score = SmithWatAffineParallel2(Spr, Tpr,
		unsigned int tmp_score = subalign(Spr, Tpr,
			tmp_align, t_left_gap, t_right_gap, mismatch_penalty,
			gap_open_penalty, gap_extend_penalty);
		if (verbose >= 2) {
		    cout << "\tscore=" << tmp_score << endl;
		    PRINT4(tmp_align.pos1(), tmp_align.Pos1(), tmp_align.pos2(),
			    tmp_align.Pos2());
		}
		if (verbose >= 3)
		    PrintVisualAlignment(True, cout, Spr, Tpr, tmp_align);

//		if ( tmp_align.extent1() == 0 && tmp_align.extent2() == 0 ) {
//		    tmp_align = alignment(0, 0, 0, avector<int>(1,0), avector<int>(1,len1) );
//		}

	    } else if (len1 == 0 && len2 != 0) {
		if (pi == 1 && !penalize_left_gap) {
		    // introduce an offset only
		    tmp_align = alignment(0, len2, 0, avector<int>(1, 0),
			    avector<int>(1, 0));
		} else if (!final || penalize_right_gap) {
		    // introduce a gap
		    int tmp_errors = gap_open_penalty
			    + (len2-1) * gap_extend_penalty;
		    tmp_align = alignment(0, 0, tmp_errors,
			    avector<int>(1, len2), avector<int>(1, 0));
		}
	    } else if (len1 != 0 && len2 == 0) {
		// introduce a gap
		int tmp_errors = gap_open_penalty + (len1-1) * gap_extend_penalty;
		tmp_align = alignment(0, 0, tmp_errors, avector<int>(1, -len1),
			avector<int>(1, 0));
	    }

            // fix up the front of the alignment:
            // convert the offset to an initial gap, unless it's the
            // first alignment and we're not penalizing a left gap
            ForceAssert(tmp_align.pos1() == 0 || tmp_align.pos2() == 0);
            if ( pi != 1 || penalize_left_gap ) {
                convert_offset_initial_gap( tmp_align );
            } else {
                // pi == 1 && !penalize_left_gap
                new_pos1 = tmp_align.pos1();
                new_pos2 = tmp_align.pos2();
            }

            // fix up the tail of the alignment:
            // add extensions if the alignment didn't cover the whole
            // space.  Since the gap comes first, we add it on to the
            // perfect interval later, if it's needed.
            ForceAssertEq( next_gap, 0 );
            if (  pi != perf_intervals.size() - 1 ) {
                int ilen1 = static_cast<int>(len1);
                int ilen2 = static_cast<int>(len2);
                int tmp_Pos1 = tmp_align.Pos1();
                int tmp_Pos2 = tmp_align.Pos2();
                if ( tmp_Pos1 <  ilen1 ) {
                    ForceAssertEq( tmp_Pos2, ilen2 );
                    next_gap = -( ilen1 - tmp_Pos1 );
                }  else if ( tmp_Pos2 < ilen2 ) {
                    ForceAssertEq( tmp_Pos1, ilen1 );
                    next_gap = ilen2 - tmp_Pos2;
                }

                if ( next_gap && verbose >= 2 )
                    cout << "NIW: adding extension: " << next_gap << endl;
            }

            // okay, now we squirrel away the current position on S and
            // T in pos1 and pos2 for later use
            set_pos1pos2( tmp_align, cur.pos1, cur.pos2);

            // add SW alignment results to the list
            new_aligns.push_back( tmp_align );
            extra_gaps.push_back( next_gap );
        }

        // add the perfect interval (and any carry-over gap from the end
        // of the SW alignment above)
        int next_gap_abs = abs( next_gap-1 );
        int next_gap_errors = next_gap_abs ? ( gap_open_penalty+gap_extend_penalty*next_gap_abs ) : 0;
        new_aligns.push_back( alignment( cur.pos1+cur.len, cur.pos2+cur.len, next_gap_errors,
                                avector<int>(1, 0),
                                avector<int>(1, cur.len) ) );
        extra_gaps.push_back(0);
    }

#if 0
    //NIW
    cout << "ALIGNS:" << endl;
    for ( auto const& align : new_aligns )
	align.Print(cout);
    cout << "EXTRA GAPS:" << endl;
    for ( auto const& gap : extra_gaps )
	cout << gap << endl;
    cout << endl;
#endif


#if 0
    // We calculate the cost of the left and right sides of the
    // alignment versus just *dumping* them for a gap. We walk forwards
    // for the left side and backwards for the right side looking for
    // the last alignment where the cumulative errors drops below the
    // cumulative cost of a gap

    if ( 0 && !penalize_left_gap  ) {
        size_t left = 0;
        int total_errors = 0;
        int proposed_gap = 0;
        int proposed_errors = 0;
        int proposed_pos1 = 0;
        int proposed_pos2 = 0;
        for ( size_t i = 0; i < new_aligns.size(); ++i ) {
            dumpit(new_aligns[i]);
            total_errors += new_aligns[i].Errors();
            int len1 = new_aligns[i].pos1();
            int len2 = new_aligns[i].pos2();
            int gap_cost_one = (len1?gap_open_penalty:0)+len1*gap_extend_penalty;
            int gap_cost_two = (len2?gap_open_penalty:0)+len2*gap_extend_penalty;
            PRINT3(gap_cost_one, gap_cost_two, total_errors);
            if ( gap_cost_one <= gap_cost_two && gap_cost_one < total_errors ) {
                left = i+1;
                proposed_gap = -new_aligns[i].pos1();
                total_errors = proposed_errors = gap_cost_one;
                proposed_pos1 = 0;
                proposed_pos2 = new_aligns[i].pos2();
            } else if ( gap_cost_two < total_errors ) {
                left = i+1;
                proposed_gap = new_aligns[i].pos2();
                total_errors = proposed_errors = gap_cost_two;
                proposed_pos1 = new_aligns[i].pos1();
                proposed_pos2 = 0;
            }
            PRINT6(i, left, proposed_gap, proposed_errors, proposed_pos1, proposed_pos2);
        }

        if ( left > 0 ) {
            vec<Bool> erase( new_aligns.size(), False );
            for (size_t i = 1; i < left; ++i )
                erase[i] = True;
            new_aligns[0] = alignment( 0, 0, proposed_errors,
                        avector<int>(1, proposed_gap), avector<int>(1, 0) );
            extra_gaps[0] = 0;
            EraseIf( new_aligns, erase);
            EraseIf( extra_gaps, erase);
            new_pos1 = proposed_pos1;
            new_pos2 = proposed_pos2;
        }
    }

    if ( 0 && !penalize_right_gap  ) {
        PRINT4( new_aligns.back().pos1(), new_aligns.back().Pos1(), new_aligns.back().pos2(), new_aligns.back().Pos2() );
        size_t right = new_aligns.size()-1;
        int total_errors = 0;
        int proposed_errors = 0;
        int zero1 = new_aligns.back().pos1();
        int zero2 = new_aligns.back().pos2();
        for ( size_t i = right; i > 0; --i ) {
            cout << endl;
            dumpit(new_aligns[i]);
            total_errors += new_aligns[i].Errors();
            int len1 = zero1 - new_aligns[i-1].pos1();
            int len2 = zero2 - new_aligns[i-1].pos2();
            int gap_cost_one = (len1?gap_open_penalty:0)+len1*gap_extend_penalty;
            int gap_cost_two = (len2?gap_open_penalty:0)+len2*gap_extend_penalty;
            PRINT3(gap_cost_one, gap_cost_two, total_errors);
            if ( gap_cost_one <= gap_cost_two && gap_cost_one < total_errors ) {
                right = i-1;
                total_errors = proposed_errors = gap_cost_one;
            } else if ( gap_cost_two < total_errors ) {
                right = i-1;
                total_errors = proposed_errors = gap_cost_two;
            }
            PRINT3(i, right, proposed_errors);
        }

        if ( right > 0 ) {
            vec<Bool> erase( new_aligns.size(), False );
            for (size_t i = new_aligns.size()-1; i > right; --i )
                erase[i] = True;
            EraseIf( new_aligns, erase);
            EraseIf( extra_gaps, erase);

            int gap1 = zero1 - new_aligns.back().pos1();
            int gap2 = zero2 - new_aligns.back().pos2();

            if ( gap1 < gap2 )
                new_aligns.push_back( alignment( 0, 0, proposed_errors, avector<int>(1, -gap1 ), avector<int>(1,0)));
            else
                new_aligns.push_back( alignment( 0, 0, proposed_errors, avector<int>(1, gap2 ), avector<int>(1,0)));
            PRINT2(gap1, gap2);
        }
    }
#endif


    // okay, bust out the alignments into gaps and lengths
    avector<int> new_align_gaps;
    avector<int> new_align_lengths;
    int new_align_errors = 0;

    for ( size_t i = 0; i < new_aligns.size(); ++i ) {
        auto const& tmp_align = new_aligns[i];

        avector<int> tmp_gaps, tmp_lengths;
        int tmp_pos1, tmp_pos2, tmp_errors;

        tmp_align.Unpack( tmp_pos1, tmp_pos2, tmp_errors, tmp_gaps, tmp_lengths );

        // tack on add-on gap from last alignment
        if ( i > 0 ) {
            int gap_len = abs( extra_gaps[i-1] );
            if ( gap_len )
                tmp_gaps(0) += extra_gaps[i-1];
            tmp_errors += ( gap_len ? gap_open_penalty : 0 ) + (gap_len-1)*gap_extend_penalty;
        }

        // okay, now tack on the gaps and lengths and errors
        new_align_errors += tmp_errors;
        ForceAssertEq( tmp_gaps.length, tmp_lengths.length );
        for ( size_t i = 0; i < tmp_gaps.length; ++i ) {
            if ( verbose >= 2 ) {
                cout << "NIW: adding gap: " << tmp_gaps(i) << endl;
                cout << "NIW: adding length: " << tmp_lengths(i) << endl;
            }
            new_align_gaps.Append( tmp_gaps(i) );
            new_align_lengths.Append( tmp_lengths(i) );
        }

    }



#if 0
        // before we add the pefect interval, we can calculate the cost
        // of *everything* that we've done so far compared with the cost
        // of simply introducing a gap here and incrementing the offset.
        if ( !penalize_left_gap && cur.len > 0 ) {
            int gap_cost_one = gap_open_penalty+cur.pos1*gap_extend_penalty;
            int gap_cost_two = gap_open_penalty+cur.pos2*gap_extend_penalty;
            if ( verbose >= 2 )
                PRINT3(gap_cost_one, gap_cost_two, new_align_errors );

            if ( gap_cost_one <= gap_cost_two && gap_cost_one < new_align_errors ) {
                if ( verbose >= 2 )
                    cout << "introducing offset on two and gap on one" << endl;
                new_align_lengths.Setsize(1);
                new_align_lengths(0) = 0;
                new_align_gaps.Setsize(1);
                new_align_gaps(0) = -cur.pos1;
                new_align_errors= gap_cost_one;
                new_pos1 = 0;
                new_pos2 = cur.pos2;
            } else if (  gap_cost_two < new_align_errors ) {
                if ( verbose >= 2 )
                    cout << "introducing offset on one and gap on two" << endl;
                new_align_lengths.Setsize(1);
                new_align_lengths(0) = 0;
                new_align_gaps.Setsize(1);
                new_align_gaps(0) = cur.pos2;
                new_align_errors= gap_cost_two;
                new_pos1 = cur.pos1;
                new_pos2 = 0;
            }
        }

        // add the perfect interval, if it's of non-zero length
        if ( cur.len > 0 ) {
                new_align_lengths.Append(cur.len);
                new_align_gaps.Append(0);
        }

        ForceAssertEq( new_align_gaps.length, new_align_lengths.length );
    } // foreach perf_interval
#endif

    if ( verbose >= 2 ) {
        PRINT3( new_pos1, new_pos2, new_align_errors );
        cout << "before NORMALIZATION...." << endl;
        cout << "GAPS LENGTHS" << endl;
        for ( size_t i = 0; i < new_align_gaps.length; ++i )
            cout << new_align_gaps(i) << " " << new_align_lengths(i) << endl;
    }

    normalize_gaps_lengths( new_pos1, new_pos2, new_align_gaps, new_align_lengths );

    if ( verbose >= 2 ) {
        PRINT3( new_pos1, new_pos2, new_align_errors );
        cout << "after NORMALIZATION...." << endl;
        cout << "GAPS LENGTHS" << endl;
        for ( size_t i = 0; i < new_align_gaps.length; ++i )
            cout << new_align_gaps(i) << " " << new_align_lengths(i) << endl;
    }


    // construct the resulting alignment
    a.Set(new_pos1, new_pos2, new_align_errors, new_align_gaps,
            new_align_lengths, new_align_lengths.length );

    if ( verbose >= 3 ) {
        PrintVisualAlignment( True, cout, S, T, a );
    }

    if ( verbose >= 1 ) {
        cout << Date() << ": SmithWatAffineSuper END OF ALIGNMENT S.size()=" << S.size() << " vs. T.size()=" << T.size() <<
        ", penalize L,R=" << penalize_left_gap << "," << penalize_right_gap << endl;
        PRINT5(a.pos1(), a.Pos1(), a.pos2(), a.Pos2(), new_align_errors);
    }

    if ( verbose >= 4 ) {
        alignment tmp_align(a);

        static size_t do_compare = 0;
        if ( ++do_compare > 0 ) {
            cout << "ORIGINAL GAPS LENGTHS" << endl;
            for ( size_t i = 0; i < new_align_gaps.length; ++i )
                cout << new_align_gaps(i) << " " << new_align_lengths(i) << endl;
            cout << endl;

            cout << "RUNNING COMPARISON " << do_compare << ": " << endl;
            SmithWatAffine( S, T, tmp_align, penalize_left_gap, penalize_right_gap );
            cout << "COMPARISON: " << endl;
            PrintVisualAlignment( True, cout, S, T, tmp_align );
            int pos1, pos2, errors;
            int Pos1, Pos2;
            Pos1 = tmp_align.Pos1();
            Pos2 = tmp_align.Pos2();
            avector<int> gaps, lengths;
            tmp_align.Unpack(pos1, pos2, errors, gaps, lengths);
            PRINT5( pos1, Pos1, pos2, Pos2, errors );
            for ( size_t i = 0; i < gaps.length; ++i )
                cout << gaps(i) << " " << lengths(i) << endl;
        }
    }

    return new_align_errors;
}
unsigned int SmithWatAffineParallel2(const basevector & S, const basevector & T,
                            alignment & a,
                            bool penalize_left_gap,
                            bool penalize_right_gap,
                            const int mismatch_penalty,
                            const int gap_open_penalty,
                            const int gap_extend_penalty
                            )
{

    ForceAssertGt(S.size(), 0u);
    ForceAssertGt(T.size(), 0u);

    unsigned int n = S.size(), N = T.size();

    if( S.size() > 120000 && T.size() > 120000){
        std::cout << "WARNING: SWA aligning two sequences of length "
                  << S.size() << " and " << T.size()
                  << std::endl;
    }

    //     ForceAssertLe( n, N );

    avector < char >s, t;
    s.resize(n);
    for (unsigned int i = 0; i < n; i++)
        s(i) = S[i];
    t.resize(N);
    for (unsigned int i = 0; i < N; i++)
        t(i) = T[i];

    const int Infinity = 100000000;
    int best_score = Infinity;
/*
    RecArray<int> score_x( n+1, N+1 );
    RecArray<int> score_y( n+1, N+1 );
    RecArray<int> score_z( n+1, N+1 );

    RecArray<unsigned char> x_from( n+1, N+1 );
    RecArray<unsigned char> y_from( n+1, N+1 );
    RecArray<unsigned char> z_from( n+1, N+1 );
    */

    struct elem_t{
        int score[3];
        unsigned char from[3];
    };
    RecArray<elem_t> matrix(n+1,N+1);

    matrix[0][0].score[0] = 0;
    matrix[0][0].score[1] = Infinity;
    matrix[0][0].score[2] = Infinity;
    matrix[0][0].from[0] = 's';
    matrix[0][0].from[1] = 's';
    matrix[0][0].from[2] = 's';

    for (unsigned int i = 1; i <= n; i++) {
        matrix[i][0].score[0] = Infinity;
        matrix[i][0].score[1] = Infinity;
        matrix[i][0].score[2] = gap_open_penalty + gap_extend_penalty * i;
        matrix[i][0].from[0] = 's';
        matrix[i][0].from[1] = 's';
        matrix[i][0].from[2] = 's';
    }

    // Zero penalty for unmatched bases on left side of T
    for (unsigned int j = 1; j <= N; j++) {
        matrix[0][j].score[0] = Infinity;
        matrix[0][j].score[1] = (penalize_left_gap ? gap_open_penalty +
                gap_extend_penalty * j :0);
        matrix[0][j].score[2] = Infinity;
        matrix[0][j].from[0] = 's';
        matrix[0][j].from[1] = 's';
        matrix[0][j].from[2] = 's';
    }

    // parallelize the smith-waterman calculation.
    // The idea is to divide the n by N matrix (score[1..n][1..N]) into many (BxB) blocks.
    // Every block will only depend on the block on it's left, up, and up-left
    // blocks, and all the block(i,j) with same i + j can be calculated parallelized.
    //
    // block size = 200 is chosen heuristically, which seems most efficient
    //
    const int block_size = 128;
    int n_block_s = n / block_size;
    if ( n_block_s * block_size < (int)n ) n_block_s += 1;
    int n_block_t = N / block_size;
    if ( n_block_t * block_size < (int)N ) n_block_t += 1;
    //cout <<  Date( ) << ": start the main loop" << endl;
    //cout <<  Date( ) << ": time used = " << TimeSince(clock) << endl;
    #pragma omp parallel
    {
    for ( int level = 0; level < n_block_s + n_block_t; ++level ) {
    // number of blocks in the diagonal that are independent to each other
    // start and stop of the block index
    // i_block_s in [0, n_block_s)
    // i_block_t in [0, n_block_t)
    // where i_block_s + i_block_t = level

    // in other words i_block_s is in ( level - n_block_t, level ] && [0, n_block_s)

    #pragma omp for schedule(dynamic,1)
    for ( int i_block_s = max( 0, level - n_block_t + 1 );
          i_block_s < min( n_block_s, level + 1 ); i_block_s++ ) {
        int i_block_t = level - i_block_s;
        ForceAssertGe( i_block_t , 0 );
        ForceAssertLt( i_block_t , n_block_t );

        int istart = i_block_s * block_size + 1;
        int istop = min( istart + block_size, (int)n + 1 );
        int jstart = i_block_t * block_size + 1;
        int jstop = min( jstart + block_size, (int)N + 1 );


    for ( int i = istart; i < istop; i++) {
        for ( int j = jstart; j < jstop; j++) {
            unsigned int x_x =
                matrix[i - 1][j - 1].score[0] + mismatch_penalty * (s(i - 1) !=
                                                            t(j - 1));
            unsigned int x_y =
                matrix[i - 1][j - 1].score[1] + mismatch_penalty * (s(i - 1) !=
                                                            t(j - 1));
            unsigned int x_z =
                matrix[i - 1][j - 1].score[2] + mismatch_penalty * (s(i - 1) !=
                                                            t(j - 1));
            //unsigned int y_x = score_x[i][j - 1] + gap_open_penalty;
            //unsigned int y_y = score_y[i][j - 1] + gap_extend_penalty;
            unsigned int y_x = matrix[i][j-1].score[0] +
                (penalize_right_gap || i != (int)n ? gap_open_penalty : 0);
            unsigned int y_y = matrix[i][j-1].score[1] +
                (penalize_right_gap || i != (int)n ? gap_extend_penalty : 0);

            unsigned int y_z = Infinity;        //score_z[i][j-1] + gap_open_penalty;
            unsigned int z_x = matrix[i - 1][j].score[0] + gap_open_penalty;
            unsigned int z_y = Infinity;        //score_y[i-1][j] + gap_open_penalty;
            unsigned int z_z = matrix[i - 1][j].score[2] + gap_extend_penalty;

            matrix[i][j].score[0] = Min(Min(x_x, x_y), x_z);
            matrix[i][j].score[1] = Min(Min(y_x, y_y), y_z);
            matrix[i][j].score[2] = Min(Min(z_x, z_y), z_z);

            if (x_x <= x_y) {
                if (x_x <= x_z)
                    matrix[i][j].from[0] = 'x';
                else
                    matrix[i][j].from[0] = 'z';
            }
            else {
                if (x_y <= x_z)
                    matrix[i][j].from[0] = 'y';
                else
                    matrix[i][j].from[0] = 'z';
            }

            if (y_x <= y_y) {
                if (y_x <= y_z)
                    matrix[i][j].from[1] = 'x';
                else
                    matrix[i][j].from[1] = 'z';
            }
            else {
                if (y_y <= y_z)
                    matrix[i][j].from[1] = 'y';
                else
                    matrix[i][j].from[1] = 'z';
            }

            if (z_x <= z_y) {
                if (z_x <= z_z)
                    matrix[i][j].from[2] = 'x';
                else
                    matrix[i][j].from[2] = 'z';
            }
            else {
                if (z_y <= z_z)
                    matrix[i][j].from[2] = 'y';
                else
                    matrix[i][j].from[2] = 'z';
            }
        }
    }



    // end parallelization
    }
    }
    }//omp parallel

    int ii = n, jj = N;
    best_score = Min(matrix[n][N].score[0], Min(matrix[n][N].score[1], matrix[n][N].score[2]));
#ifdef FIX_RIGHT_GAP
    int right = Min(n,N);
#else
    int right = 0;
#endif
    if (!penalize_right_gap) {
        for(int k = N; k >= right; k--){
            int best_score_k = Min( matrix[n][k].score[0], Min( matrix[n][k].score[1], matrix[n][k].score[2] ) );
            if (best_score_k < best_score) {
                best_score = best_score_k;
                ii = n;
                jj = k;
            }
        }
    }

    //vec < vec < unsigned char > >*from;
//    RecArray<unsigned char> *from;
    int from_idx=-1;
    if (matrix[ii][jj].score[0] <= matrix[ii][jj].score[1]) {
        if (matrix[ii][jj].score[0] <= matrix[ii][jj].score[2])
            from_idx=0;//from = &x_from;
        else
            from_idx=2;//from = &z_from;
    }
    else {
        if (matrix[ii][jj].score[1] <= matrix[ii][jj].score[2])
            from_idx=1;//from = &y_from;
        else
            from_idx=2;//from = &z_from;
    }

    int i = ii;
    int j = jj;
    int lcount = 0, g1count = 0, g2count = 0;
    int last_length = 0;
    avector < int >gaps(0), lengths(0);
    while (1) {
        unsigned char dir = matrix[i][j].from[from_idx];//(*from)[i][j];
        //cout << dir;
        if (from_idx==0/*from == &x_from*/) {
            if (g1count > 0) {
                if (last_length > 0) {
                    gaps.Prepend(g1count);
                    lengths.Prepend(last_length);
                }
                g1count = 0;
            }
            if (g2count > 0) {
                if (last_length > 0) {
                    gaps.Prepend(-g2count);
                    lengths.Prepend(last_length);
                }
                g2count = 0;
            }
            ++lcount;
            --i;
            --j;
        }
        else if (from_idx==2/*from == &z_from*/)       // gap on long sequence
        {
            if (lcount > 0) {
                last_length = lcount;
                lcount = 0;
            }
            ForceAssert(g1count == 0);
            ++g2count;
            --i;
        }
        else                    // gap on short sequence
        {
            if (lcount > 0) {
                last_length = lcount;
                lcount = 0;
            }
            ForceAssert(g2count == 0);
            ++g1count;
            --j;
        }

        if (dir == 'x')
            from_idx=0;//from = &x_from;
        else if (dir == 'y')
            from_idx=1;//from = &y_from;
        else
            from_idx=2;//from = &z_from;

        if (matrix[i][j].from[from_idx]/*(*from)[i][j]*/ == 's')
            break;

    }

    //cout << "\n";

    if (g1count != 0)
        gaps.Prepend(g1count);
    else if (g2count != 0)
        gaps.Prepend(-g2count);
    else
        gaps.Prepend(0);

    lengths.Prepend(lcount);

    int pos1 = i;
    int pos2 = j;

    if (gaps(0) < 0) {
        pos2 -= gaps(0);
        gaps(0) = 0;
    }

    if (gaps(0) > 0) {
        pos1 += gaps(0);
        gaps(0) = 0;
    }

    int errors = best_score;
    a = alignment(pos1, pos2, errors, gaps, lengths);

    return best_score;
}

