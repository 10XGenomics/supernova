///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "CoreTools.h"
#include "Qualvector.h"
#include "paths/BigMapTools.h"
#include "paths/long/CorrectByStack.h"
#include "paths/long/MakeKmerStuff.h"
#include "paths/long/ReadStack.h"

void CorrectByStack( StackBaseVecVec& call, StackQualVecVec& callq,
     basevector& bases_new, qualvector& quals_new, int& trim_to,
     vec< triple<int,int,Bool> >& offset_id_rc2, const Bool verbose )
{

     // Heuristics.

     const int min_win = 50;
     const int min_win_ratio = 10;
     const int max_lose = 100;

     // Set up.

     int n = call.size( ), k = call[0].size( );
     trim_to = k;
     const int id1 = offset_id_rc2[0].second; // for verbose output
     if (verbose) cout << "\n";

     // Identify reads that have a Q30 mismatch with the founder.  These are
     // suspected of being false friends.

     const int critical_q = 30;
     vec<Bool> suspect( n, False );
     for ( int j = 1; j < n; j++ )
     {    for ( int c = 0; c < k; c++ )
          {    if ( call[j][c] != call[0][c]
                    && callq[j][c] >= critical_q && callq[0][c] >= critical_q )
               {    suspect[j] = True;
                    break;    }    }    }

     // Raise quality scores.  For each 11 base window, suppose first that the
     // founder has not been edited on the window, except possibly for raising 
     // quality scores, and has score < 30 at the middle base.  Suppose in addition 
     // that there are at least 3 reads that agree with the founder, have not been 
     // edited on the window, except possibly for editing quality scores, have score 
     // > 30 at the middle base, and have no Q30 difference with the founder.
     // Then change the score of the middle base on the given read to 30.  Note that
     // this does not affect bases produced by Correct1, unless you run it a 
     // second time.
     //
     // Amelioration.  Don't raise q score if there is a viable alternate 
     // hypothesis.  This is a sequence in the window that agrees with the founder,
     // except for the middle base, and which occurs at least 3 times in reads 
     // that have no Q30 difference with the founder.

     const int rwindow = 11;
     const int min_agree = 3;
     vec<int> raises;
     for ( int c = 0; c <= k - rwindow; c++ )
     {    if ( callq[0][ c + rwindow/2 ] >= critical_q ) continue;

          /*
          Bool edited = False;
          for ( int l = 0; l < rwindow; l++ )
               if ( callq[0][c+l] == 0 ) edited = True;
          if (edited) continue;
          */

          if ( callq[0][ c + rwindow/2 ] == 0 ) continue;

          int support = 0;
          for ( int j = 1; j < n; j++ )
          {    if ( suspect[j] || callq[j][ c + rwindow/2 ] < critical_q ) continue;
               Bool bad = False;
               for ( int l = 0; l < rwindow; l++ )
               {    if ( callq[j][c+l] <= 0 || call[j][c+l] != call[0][c+l] )
                    {    bad = True;
                         break;    }    }
               if (bad) continue;
               if ( ++support == min_agree ) break;    }
          if ( support < min_agree ) continue;

          int alts[4]{0,0,0,0};
          for ( int j = 1; j < n; j++ )
          {    if ( suspect[j] ) continue;
               Bool bad = False;
               for ( int l = 0; l < rwindow; l++ )
               {    if ( l == rwindow/2 ) continue;
                    if ( callq[j][c+l] <= 0 || call[j][c+l] != call[0][c+l] )
                    {    bad = True;
                         break;    }    }
               if (bad) continue;
               char x = call[j][ c + rwindow/2 ];
               if ( x != call[0][ c + rwindow/2 ] ) alts[int(x)]++;    }
          if ( *std::max_element(alts,alts+4) >= min_agree ) continue;

          if (verbose) raises.push_back( c + rwindow/2 );
          // callq[0][ c + rwindow/2 ] = critical_q;    }
          quals_new[ c + rwindow/2 ] = critical_q;    }
     if ( verbose && raises.nonempty( ) )
     {
          #pragma omp critical
          {    cout << "raising q score for read " << id1 << ", bases ";
               for ( int i = 0; i < raises.isize( ); i++ )
               {    if ( i > 0 ) cout << ", ";
                    int j;
                    for ( j = i + 1; j < raises.isize( ); j++ )
                         if ( raises[j] - raises[j-1] != 1 ) break;
                    if ( j - i == 1 ) cout << raises[i];
                    else cout << raises[i] << "-" << raises[j-1];
                    i = j - 1;    }
               cout << endl;    }    }

     // Find motifs and use them to delete putatively false friends.

     vec<Bool> to_delete( n, False );
     int width = 10;
     int min_mult = 10;
     for ( int i = 0; i <= k - width; i += width )
     {    if ( i + width > k ) break;
          vec< vec<char> > X;
          vec<int> ids;
          for ( int j = 0; j < n; j++ )
          {    Bool bad = False;
               for ( int l = i; l < i + width; l++ )
                    if ( call[j][l] == ' ' ) bad = True;
               if (bad) continue;
               vec<char> x(width);
               for ( int l = i; l < i + width; l++ )
                    x[l-i] = as_base( call[j][l] );
               X.push_back(x);
               ids.push_back(j);    }
          SortSync( X, ids );
          Bool found = False;
          vec<int> bigs;
          int this_one = -1;
          for ( int r = 0; r < X.isize( ); r++ )
          {    int s = X.NextDiff(r);
               if ( s - r >= min_mult )
               {    Bool agree = True;
                    for ( int l = 0; l < width; l++ )
                    {    if ( X[r][l] != as_base( call[0][i+l] ) ) 
                              agree = False;    }
                    if (agree) this_one = bigs.size( );
                    bigs.push_back(r);    }
               r = s - 1;    }    
          if ( this_one >= 0 )
          {    for ( int b = 0; b < bigs.isize( ); b++ )
               {    if ( b == this_one ) continue;
                    Bool hq_diff = False;
                    for ( int l = 0; l < width; l++ )
                    {    if ( X[ bigs[b] ][l] != X[this_one][l] )
                         {    if ( callq[0][i+l] >= 20 )
                                   hq_diff = True;    }    }
                    if ( !hq_diff ) continue;
                    for ( int d = bigs[b]; d < X.isize( ); d++ )
                    {    if ( X[d] != X[ bigs[b] ] ) break;
                         to_delete[ ids[d] ] = True;    }    }    }    }
     EraseIf( offset_id_rc2, to_delete );
     call.EraseIf( to_delete );
     callq.EraseIf( to_delete );
     EraseIf( suspect, to_delete );
     n = call.size( );

     // Remove a friend if it has a Q30+/Q30+ mismatch with the founder.


     suspect.resize_and_set( n, False );
     for ( int j = 1; j < n; j++ )
     {    for ( int c = 0; c < k; c++ )
          {    if ( call[j][c] != call[0][c]
                    && callq[j][c] >= critical_q && callq[0][c] >= critical_q )
               {    suspect[j] = True;
                    break;    }    }    }



     EraseIf( offset_id_rc2, suspect );
     call.EraseIf( suspect );
     callq.EraseIf( suspect );
     n = call.size( );

     // Remove friends having inadequate glue to the founder.

     const int min_glue = 20;
     const int max_homopolymer_in_glue = 10;
     to_delete.resize_and_set( n, False );
     for ( int j = 1; j < n; j++ )
     {    vec< vec<char> > glue, glue2;
          for ( int r = 0; r < k; r++ )
          {    if ( call[j][r] != call[0][r] ) continue;
               int s;
               for ( s = r + 1; s < k; s++ )
                    if ( call[j][s] != call[0][s] ) break;
               vec<char> g;
               if ( s - r >= min_glue )
               {    for ( int l = r; l < s; l++ )
                         g.push_back( call[j][l] );    }
               glue.push_back(g);
               r = s;    }
          for ( int i = 0; i < glue.isize( ); i++ )
          {    const vec<char>& g = glue[i];
               vec<char> gnew;
               int hcount = 1;
               for ( int l = 0; l < g.isize( ); l++ )
               {    if ( gnew.nonempty( ) && g[l] == gnew.back( ) ) hcount++;
                    else hcount = 1;
                    if ( hcount <= max_homopolymer_in_glue )
                         gnew.push_back( g[l] );    }
               glue2.push_back(gnew);    }
          int M = 0;
          for ( int i = 0; i < glue2.isize( ); i++ )
               M = Max( M, glue2[i].isize( ) );
          if ( M < min_glue ) to_delete[j] = True;    }
     EraseIf( offset_id_rc2, to_delete );
     call.EraseIf( to_delete );
     callq.EraseIf( to_delete );
     n = call.size( );

     // Trace read.

     if (verbose) TraceRead( call, callq, offset_id_rc2 );

     // Define multiplicity.  Note that this is neutered, because although
     // it makes sense, there is of yet no evidence that it helps.

     vec<int> mult(n);
     vec<int> where( n, vec<int>::IDENTITY );
     vec<int> ids(n);
     for ( int i = 0; i < n; i++ )
          ids[i] = offset_id_rc2[i].second;
     SortSync( ids, where );
     for ( int i = 0; i < n; i++ )
     {    int j = ids.NextDiff(i);
          for ( int k = i; k < j; k++ )
               mult[ where[k] ] = j - i;
          i = j - 1;    }
     for ( int i = 0; i < n; i++ )
          mult[i] = 1;   // NEUTERING!!

     // Go through the columns.

     for ( int i = 0; i < k; i++ )
     {    
          // Compute quality score sum for each base.  We count Q2 bases as
          // next to nothing.

          vec<double> sum( 4, 0 ); 
          vec<int> ids( 4, vec<int>::IDENTITY );
          for ( int j = 0; j < n; j++ )
          {    double q = callq[j][i];
               if ( q <= 2 ) q = Min( q, 0.2 );
               if ( callq[j][i] >= 0 ) sum[ call[j][i] ] += q/mult[j];    }
          ReverseSortSync( sum, ids );

          int winner = ids[0];
          vec< pair<int,int> > competitors;
          vec<double> top( 4, 0 );
          for ( int j = 0; j < n; j++ )
          {    int b = call[j][i];
               if ( b != ' ' ) 
                    top[b] = Max( top[b], callq[j][i] / double( mult[j] ) );    }

          // Drop the top score for competitors.

          for ( int j = 1; j < 4; j++ )
               sum[j] -= top[ ids[j] ];

          Bool OK = False;
          if ( sum[0] >= min_win && sum[0] >= min_win_ratio * sum[1]
               && sum[1] <= max_lose )
          {    OK = True;    }
          if (OK) 
          {    if ( call[0][i] != winner )
               {    if (verbose)
                    {    
                         #pragma omp critical
                         {    cout << "trace read " << id1 << ", changing base "
                                   << i << " from " << as_base( call[0][i] ) 
                                   << " to " << as_base(winner) << endl;    }    }
                    bases_new.Set( i, winner );
                    quals_new[i] = 0;    }    }
          else 
          {    trim_to = i;
               if (verbose)
               {    
                    #pragma omp critical
                    {    cout << "\ntrimming trace read " << id1 << " to " << i 
                              << " bases" << endl;    }    }
               break;    }    }    }

void TraceRead( StackBaseVecVec& call, StackQualVecVec& callq,
     vec< triple<int,int,Bool> >& offset_id_rc2 )
{
     // Compute provisional consensus.

     int n = call.size( ), k = call[0].size( );
     vec<char> con(k);
     for ( int c = 0; c < k; c++ )
     {    vec<int> count( 4, 0 ), ids( 4, vec<int>::IDENTITY );
          for ( int j = 0; j < n; j++ )
          {    char x = call[j][c];
               int q = callq[j][c];
               if ( x != ' ' ) count[x] += q;    }
          ReverseSortSync( count, ids );
          con[c] = ids[0];    }

     // Set up for quality score display.

     vec< vec<char> > callq1( n, vec<char>( k, ' ' ) );
     vec< vec<char> > callq2( n, vec<char>( k, ' ' ) );
     for ( int j = 0; j < n; j++ )
     {    for ( int c = 0; c < k; c++ )
          {    int qa = callq[j][c] / 10, qb = callq[j][c] % 10;
               if ( callq[j][c] >= 0 )
               {    callq1[j][c] = ( qa == 0 ? ' ' : '0' + qa );
                    callq2[j][c] = '0' + qb;    }    }    }

     // Report glue.

     /*
     #pragma omp critical
     {    for ( int j = 1; j < n; j++ )
          {    cout << "\nglue for friend " << j << ", id = " 
                    << offset_id_rc2[j].second << " "
                    << ( offset_id_rc2[j].third ? "rc" : "fw" )
                    << ", offset = " << offset_id_rc2[j].first << "\n";
               for ( int r = 0; r < k; r++ )
               {    if ( call[j][r] != call[0][r] ) continue;
                    int s;
                    for ( s = r + 1; s < k; s++ )
                         if ( call[j][s] != call[0][s] ) break;
                    if ( s - r >= 20 )
                    {    for ( int l = r; l < s; l++ )
                              cout << as_base( call[j][l] );
                         cout << "\n";    }
                    r = s;    }    }    }
     */

     // Display windows.

     const int id1 = offset_id_rc2[0].second;
     #pragma omp critical
     {    cout << "\nwindows for read " << id1 << endl;
          const int w = 70;
          for ( int cstart = 0; cstart < k; cstart += w )
          {    int c1 = cstart, c2 = Min( cstart + w, k );
               int wid = cstart/w + 1;
               cout << "\n=================================================="
                    << "==================================" << endl;
               cout << "\n\n" << "[" << wid << "] window from " << c1 
                    << " to " << c2 << "\n\n\n";
          
               // Windows are sorted by quality = number of matches.
          
               vec<int> matches( n, 0 ), ids( n, vec<int>::IDENTITY );
               matches[0] = 1000000000;
               for ( int j = 1; j < n; j++ )
               {    for ( int c = c1; c < c2; c++ )
                    {    if ( call[0][c] == ' ' || call[j][c] == ' ' ) continue;
                         if ( con[c] == call[j][c] ) matches[j]++;    }    }
               ReverseSortSync( matches, ids);

               for ( int jm = 0; jm < n; jm++ )
               {    int j = ids[jm];

                    vec<String> tag(4);
                    tag[0] = "     id=" + ToString( offset_id_rc2[j].second );
                    tag[1] = String("     ") 
                         + ( offset_id_rc2[j].third ? "rc" : "fw" );
                    tag[2] = "     " + ToString( offset_id_rc2[j].first );
                    tag[3] = "     window [" + ToString(wid) + "]";
                    int tt = 0;

                    Bool present = False, nonblank = False;
                    for ( int c = c1; c < c2; c++ )
                         if ( call[j][c] != ' ' ) present = True;
                    if ( !present ) continue;

                    for ( int c = c1; c < c2; c++ )
                         cout << ( call[j][c] == ' ' ? ' ' : '-' );
                    cout << "\n";

                    for ( int c = c1; c < c2; c++ )
                         if ( callq1[j][c] != ' ' ) nonblank = True;
                    if (nonblank)
                    {    for ( int c = c1; c < c2; c++ )
                              cout << callq1[j][c];
                         if ( tt < 4 ) cout << tag[tt++];
                         cout << "\n";    }
                    for ( int c = c1; c < c2; c++ )
                         cout << callq2[j][c];
                    if ( tt < 4 ) cout << tag[tt++];
                    cout << "\n";
                    for ( int c = c1; c < c2; c++ )
                         cout << ( call[j][c] == ' ' ? ' ' : ':' );
                    if ( tt < 4 ) cout << tag[tt++];
                    cout << "\n";

                    int ms = 0;
                    for ( int pass = 1; pass <= 2; pass++ )
                    {    if ( pass == 2 && ms == 0 ) break;
                         for ( int c = c1; c < c2; c++ )
                         {    char m;  
                              if ( call[0][c] == ' ' ) m = ' ';
                              else if ( call[j][c] == ' ' ) m = ' ';
                              else if ( con[c] == call[j][c] ) m = ' ';
                              else m = '*';
                              if ( m == '*' ) ms++;    
                              if ( pass == 2 ) cout << m;    }
                         if ( pass == 2 ) 
                         {    if ( tt < 4 ) cout << tag[tt++];
                              cout << "\n";    }    }
                    for ( int c = c1; c < c2; c++ )
                    {    if ( call[j][c] == ' ' ) cout << ' '; // ?????????????
                         else cout << as_base( call[j][c] );    }
                    if ( tt < 4 ) cout << tag[tt++];
                    cout << "\n";
                    for ( int c = c1; c < c2; c++ )
                         cout << ( call[j][c] == ' ' ? ' ' : '-' );
                    cout << "\n\n\n";    }    }

          // Analyze.
     
          for ( int c = 0; c < k; c++ )
          {    vec<char> base;
               vec<int> qual;
               vec<int> total( 4, 0 ), id( 4, vec<int>::IDENTITY );
               for ( int j = 0; j < n; j++ )
               {    base.push_back( call[j][c] );
                    qual.push_back( callq[j][c] );    
                    if ( call[j][c] != ' ' )
                    {    total[ call[j][c] ] += callq[j][c];    }    }

               vec< pair<int,char> > qb;
               vec<int> ids, origin;
               for ( int j = 0; j < base.isize( ); j++ )
               {    qb.push( qual[j], base[j] );
                    ids.push( offset_id_rc2[j].second );
                    origin.push_back(j);    }
               ReverseSortSync( qb, ids, origin );

               // Call a position 'squeaky clean' if it is supported by 3 Q30+
               // bases and there is nothing else except Q2s.  In that case
               // we're done and don't print.

               int top_competitor = 0, tc_pos = -1, tc_id = -1;
               for ( int j = 1; j < qb.isize( ); j++ )
               {    if ( qb[j].second != qb[0].second )
                    {    top_competitor = qb[j].first;
                         tc_pos = j;
                         tc_id = ids[j];
                         break;    }    }
               if ( top_competitor <= 2 && tc_pos >= 3 && qb[2].first >= 30 ) 
                    continue;

               cout << "\ncolumn " << c << "\n";
               for ( int j = 0; j < base.isize( ); j++ )
               {    if ( j > 0 ) cout << " ";
                    int k;
                    for ( k = j + 1; k < base.isize( ); k++ )
                         if ( qb[k] != qb[j] ) break;
                    if ( qb[j].first >= 0 )
                    {    cout << as_base(qb[j].second) << "[" << qb[j].first << "]";
                         if ( k - j > 1 ) cout << "*" << k-j;    }
                    j = k - 1;    }
               ReverseSortSync( total, id );
               cout << "  --  ";
               for ( int j = 0; j < 4; j++ )
               {    if ( total[j] == 0 ) break;
                    cout << as_base( id[j] ) << "[" << total[j] << "]";    }
               cout << "\n";
               
               // Assess Q30 competitors.
          
               cout << "winning base = " << as_base( id[0] ) << endl;
               vec<int> competitors;
               for ( int j = 0; j < base.isize( ); j++ )
               {    if ( qb[j].first >= 20 && qb[j].second != id[0] )
                    {    cout << "competitor = " << ids[j] << ", base = " 
                              << as_base(qb[j].second) << ", qual = " 
                              << qb[j].first << endl;
                    competitors.push_back(j);    }    }    }    }    }
