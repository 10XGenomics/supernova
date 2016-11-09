///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"
#include "ParallelVecUtilities.h"
#include "VecUtilities.h"
#include "math/Functions.h"
#include "paths/HyperBasevector.h"
#include "paths/long/CreateGenome.h"
#include "paths/long/DigraphFromWords.h"
#include "paths/long/EvalByReads.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/SupportedHyperBasevector.h"

void Snorkle( const vec<int>& w, const vec< vec<int> >& dindex,
     const vec<int>& to_left, const vec<int>& to_right, vec< vec<int> >& paths )
{    vec< vec<int> > places( w.size( ) );
     for ( int l = 0; l < w.isize( ); l++ )
     for ( int m = 0; m < dindex[ w[l] ].isize( ); m++ )
          places[l].push_back( dindex[ w[l] ][m] );
     paths.clear( );
     paths.push_back( vec<int>( ) );
     for ( int l = 0; l < w.isize( ); l++ )
     {    vec< vec<int> > paths2;
          for ( int p = 0; p < paths.isize( ); p++ )
          {    const vec<int>& q = paths[p];
               for ( int r = 0; r < places[l].isize( ); r++ )
               {    vec<int> x(q);
                    if ( q.nonempty( ) 
                         && to_right[ q.back( ) ] != to_left[ places[l][r] ] )
                    {    continue;    }
                    x.push_back( places[l][r] );
                    paths2.push_back(x);    }    }
          paths = paths2;    }    }

void SupportedHyperBasevector::WordifyAlt2( const long_logging_control& log_control, 
     const long_logging& logc )
{
     // Set up index to paths.

     double clock = WallClockTime( );
     vec< vec<int> > w = Paths( );
     int wmax = -1;
     for ( int i = 0; i < w.isize( ); i++ )
     for ( int j = 0; j < w[i].isize( ); j++ )
     {    ForceAssertGe( w[i][j], 0 );
          wmax = Max( wmax, w[i][j] );    }
     vec< vec< pair<int,int> > > windex( wmax + 1 );
     for ( int i = 0; i < w.isize( ); i++ )
     for ( int j = 0; j < w[i].isize( ); j++ )
               windex[ w[i][j] ].push( i, j );

     // Note trace ids.

     vec<int> trace;
     ParseIntSet( "{" + logc.WORDIFY_TRACE + "}", trace );

     // Find overlaps between paths.

     if (logc.STATUS_LOGGING) cout << Date( ) << ": wordify finding overlaps" << endl;
     vec< vec< triple<int,int,int> > > over( w.size( ) );
     #pragma omp parallel for schedule( dynamic, 1 )
     for ( int i1 = 0; i1 < w.isize( ); i1++ )
     {    vec< triple<int,int,int> > x;
          for ( int j1 = 0; j1 < w[i1].isize( ); j1++ )
          {    for ( int l = 0; l < windex[ w[i1][j1] ].isize( ); l++ )
               {    int i2 = windex[ w[i1][j1] ][l].first;
                    int j2 = windex[ w[i1][j1] ][l].second;
                    int o = j1 - j2;
                    int ov = IntervalOverlap(
                         0, w[i1].isize( ), o, o + w[i2].isize( ) );
                    if ( i1 == i2 && o == 0 ) continue;
                    if ( Overlap( w[i1], w[i2], o ) ) x.push( i2, o, ov );    }    }
          UniqueSort(x);
          over[i1] = x;    }

     // Find overlap segments.  Note that we could choose to generalize the notion 
     // of segment, so that it is not restricted to those defined by overlaps.

     if (logc.STATUS_LOGGING)
          cout << Date( ) << ": wordify finding overlap segments" << endl;
     vec< vec<int> > segs;
     const int batches = 1000;
     int batch = Max( 1, NPaths( ) / batches );
     #pragma omp parallel for schedule( dynamic, 1 )
     for ( int bi1 = 0; bi1 < NPaths( ); bi1 += batch )
     {    vec< vec<int> > segsi;
          for ( int i1 = bi1; i1 < Min( bi1 + batch, NPaths( ) ); i1++ )
          {    for ( int j1 = 0; j1 < over[i1].isize( ); j1++ )
               {    int i2 = over[i1][j1].first, offset = over[i1][j1].second;
                    int start1 = 0, stop1 = Path(i1).size( );
                    int start2 = offset, stop2 = offset + Path(i2).isize( );
                    int low = Max( start1, start2 ), high = Min( stop1, stop2 );
                    vec<int> s;
                    for ( int l = low; l < high; l++ )
                         s.push_back( Path(i1)[l] );
                    segsi.push_back(s);    }    }
          UniqueSort(segsi);
          #pragma omp critical
          {    segs.append(segsi);    }    }
     if (logc.STATUS_LOGGING)
          cout << Date( ) << ": wordify sorting overlap segments" << endl;
     ParallelUniqueSort(segs);
     if (logc.STATUS_LOGGING)
          cout << Date( ) << ": found " << segs.size( ) << " segments" << endl;

     // Find the paths containing each segment.

     if (logc.STATUS_LOGGING) cout << Date( ) << ": wordify finding paths" << endl;
     vec< vec< pair<int,int> > > slocs( segs.size( ) );
     #pragma omp parallel for
     for ( int s = 0; s < segs.isize( ); s++ )
     {    int v = segs[s][0];
          for ( int u = 0; u < windex[v].isize( ); u++ )
          {    int id = windex[v][u].first, pos = windex[v][u].second;
               if ( pos + segs[s].isize( ) > Path(id).isize( ) ) continue;
               Bool mismatch = False;
               for ( int l = 1; l < segs[s].isize( ); l++ )
               {    if ( segs[s][l] != Path(id)[ pos + l ] )
                    {    mismatch = True;
                         break;    }    }
               if ( !mismatch ) slocs[s].push( id, pos );    }    }

     // Set up data structure to record paths destined for deletion.  This
     // is deducible from start but convenient to record separately.
     //
     // Note the issue of whether paths should be deleted, or just truncated.

     vec<Bool> delp( NPaths( ), False );
     vec< vec<Bool> > over_del( w.size( ) );
     for ( int i = 0; i < w.isize( ); i++ )
          over_del[i].resize( over[i].size( ), False );

     // Go through the overlap segments.

     if (logc.STATUS_LOGGING) cout << Date( ) << ": begin wordify main loop" << endl;
     #pragma omp parallel for schedule( dynamic, 1 )
     for ( int s = 0; s < segs.isize( ); s++ )
     {    int nkmers = 0;
          int w = segs[s].size( );
          for ( int l = 0; l < w; l++ )
               nkmers += EdgeLengthKmers( segs[s][l] );
          int left = 0, right = 0;
          for ( int j = 0; j < slocs[s].isize( ); j++ )
          {    left = Max( left, slocs[s][j].second );
               int id = slocs[s][j].first;
               right = Max( right, Path(id).isize( ) 
                    - slocs[s][j].second - segs[s].isize( ) );    }

          // Form matrix.

          vec< vec<int> > M;
          vec<fix64_6> W;
          for ( int j = 0; j < slocs[s].isize( ); j++ )
          {    vec<int> row;
               for ( int l = 0; l < left - slocs[s][j].second; l++ )
                    row.push_back( -1 );
               int id = slocs[s][j].first;
               for ( int l = 0; l < Path(id).isize( ); l++ )
                    row.push_back( Path(id)[l] );
               int pad = right - ( Path(id).isize( ) - slocs[s][j].second - w );
               for ( int l = 0; l < pad; l++ )
                    row.push_back( -1 );
               M.push_back(row);
               W.push_back( Weight(id) );    }

          // Test total coverage.

          const int min_mult = 10;
          const int max_weak = 2;
          fix64_6 sw = 0;
          for ( int j = 0; j < W.isize( ); j++ )
               sw += W[j];
          if ( sw < min_mult ) continue;

          // Don't bother unless some path crosses the segment.

          Bool crossed = False;
          for ( int j = 0; j < slocs[s].isize( ); j++ )
          {    int id = slocs[s][j].first;
               if ( slocs[s][j].second > 0 && Path(id).isize( )
                    - slocs[s][j].second - w > 0 )
               {    crossed = True;    }    }
          if ( !crossed ) continue;

          // Don't bother if the data are consistent.

          Bool inconsistent = False;
          int mc = left + segs[s].isize( ) + right;
          for ( int i = 0; i < mc; i++ )
          {    vec<int> x;
               for ( int j = 0; j < M.isize( ); j++ )
                    if ( M[j][i] >= 0 ) x.push_back( M[j][i] );
               UniqueSort(x);
               if ( !x.solo( ) ) inconsistent = True;    }
          if ( !inconsistent ) continue;

          // Determine the number of kmers that the 'core' segment could be 
          // extended on either side.  These sets could be made smaller, which
          // would make the subsequent computations more efficient.

          vec<int> left_ext, right_ext;
          int cols = left + w + right;
          for ( int i = 0; i < M.isize( ); i++ )
          {    int lext = 1, rext = 1;
               for ( int j = left - 1; j >= 0; j-- )
               {    if ( M[i][j] < 0 ) break;
                    left_ext.push_back(lext);
                    lext += EdgeLengthKmers( M[i][j] );    }
               for ( int j = left + w; j < cols; j++ )
               {    if ( M[i][j] < 0 ) break;
                    right_ext.push_back(rext);
                    rext += EdgeLengthKmers( M[i][j] );    }    }
          UniqueSort(left_ext), UniqueSort(right_ext);

          // Now traverse all possible extensions.  

          int64_t joint_exts = left_ext.size( ) * right_ext.size( );
          const int64_t max_joint_exts = 1000;
          if ( joint_exts > max_joint_exts ) continue;
          ostringstream out;
          if ( logc.verb[ "WORDIFY" ] >= 2 )
          {    out << "\n=========================================================="
                    << "======================\n";    }
          for ( int e1 = 0; e1 < left_ext.isize( ); e1++ )
          for ( int e2 = 0; e2 < right_ext.isize( ); e2++ )
          {
                // Build truncated matrix.

               if ( logc.verb[ "WORDIFY" ] >= 2 )
               {    out << "\nseg = " << printSeq(segs[s]) << ", left_ext = " 
                         << left_ext[e1]
                         << ", right_ext = " << right_ext[e2] << endl;    }
               vec< vec<int> > MT(M);
               int mrows = M.size( );
               vec<Bool> left_full( mrows, False ), right_full( mrows, False );
               vec< pair<Bool,fix64_6> > WT(mrows);
               for ( int i = 0; i < MT.isize( ); i++ )
               {    int total = 0;
                    for ( int j = left - 1; j >= 0; j-- )
                    {    if ( left_full[i] ) MT[i][j] = -1;
                         if ( MT[i][j] >= 0 )
                         {    total += EdgeLengthKmers( MT[i][j] );
                              if ( total >= left_ext[e1] ) 
                                   left_full[i] = True;    }    }
                    total = 0;
                    for ( int j = left + w; j < cols; j++ )
                    {    if ( right_full[i] ) MT[i][j] = -1;
                         if ( MT[i][j] >= 0 )
                         {    total += EdgeLengthKmers( MT[i][j] );
                              if ( total >= right_ext[e2] ) 
                                   right_full[i] = True;    }    }
                    WT[i] = make_pair( left_full[i] && right_full[i], W[i] );    }
               SortSync( MT, WT, left_full, right_full );
               vec<Bool> to_delete( MT.size( ), False );
               for ( int i = 0; i < MT.isize( ); i++ )
               {    int j = MT.NextDiff(i);
                    for ( int k = i + 1; k < j; k++ )
                    {    WT[i].second += WT[k].second;
                         to_delete[k] = True;    }
                    i = j - 1;    }
               EraseIf( MT, to_delete ), EraseIf( WT, to_delete );
               EraseIf( left_full, to_delete), EraseIf( right_full, to_delete );
               ReverseSortSync( WT, MT, left_full, right_full );

               // See if minimum threshold is passed.

               if ( !( WT.nonempty( ) && WT[0].first && WT[0].second >= min_mult ) )
                    continue;

               // Index nonzero elements of MT.

               int nrows = MT.size( );
               vec<int> start(nrows), stop(nrows);
               for ( int i = 0; i < nrows; i++ )
               {    int p = 0;
                    while( MT[i][p] < 0 ) p++;
                    start[i] = p;
                    while( p < cols && MT[i][p] >= 0 ) p++;
                    stop[i] = p;    }

               // Add in potential joins, assigning weight 0.

               for ( int i1 = 0; i1 < nrows; i1++ )
               {    if ( !left_full[i1] || MT[i1][left+w] >= 0 ) continue;
                    for ( int i2 = 0; i2 < nrows; i2++ )
                    {    if ( !right_full[i2] || MT[i2][left-1] >= 0 ) continue;
                         vec<int> join( left + w + right, -1 );
                         for ( int l = start[i1]; l < stop[i1]; l++ )
                              join[l] = MT[i1][l];
                         for ( int l = start[i2]; l < stop[i2]; l++ )
                              join[l] = MT[i2][l];
                         if ( Member( MT, join ) ) continue;
                         MT.push_back(join), WT.push( True, 0 );
                         left_full.push_back(True); 
                         right_full.push_back(True);    }    }
               nrows = MT.size( );

               // Index nonzero elements of MT (again).

               start.resize(nrows), stop.resize(nrows);
               for ( int i = 0; i < nrows; i++ )
               {    int p = 0;
                    while( MT[i][p] < 0 ) p++;
                    start[i] = p;
                    while( p < cols && MT[i][p] >= 0 ) p++;
                    stop[i] = p;    }

               // Define left ext and right ext for each row.

               vec<int> lext(nrows,0), rext(nrows,0);
               for ( int i = 0; i < nrows; i++ )
               {    for ( int j = left - 1; j >= 0; j-- )
                    {    if ( j == 0 || MT[i][j-1] < 0 )
                         {    lext[i]++;
                              break;    }
                         else lext[i] += EdgeLengthKmers( MT[i][j] );    }
                    for ( int j = left + w; j < cols; j++ )
                    {    if ( j == left + w + right - 1 || MT[i][j+1] < 0 )
                         {    rext[i]++;
                              break;    }
                         else rext[i] += EdgeLengthKmers( MT[i][j] );    }    }

               // Define keepers.

               vec<Bool> keep( nrows, False );
               for ( int i = 0; i < nrows; i++ )
               {    if ( WT[i].first && ( WT[i].second > max_weak 
                         || WT[i].second * min_mult > WT[0].second ) )
                    {    keep[i] = True;    }    }

               // Define followers.

               vec<Bool> follow( nrows, False );
               for ( int i = 0; i < nrows; i++ )
               {    if ( WT[i].first ) continue;
                    for ( int j = 0; j < MT.isize( ); j++ )
                    {    if ( !keep[j] ) continue;
                         Bool OK = True;
                         for ( int l = start[i]; l < stop[i]; l++ )
                         {    if ( MT[i][l] != MT[j][l] )
                              {    OK = False;
                                   break;    }    }
                         if (OK) 
                         {    follow[i] = True;
                              break;    }    }    }

               // Define losers.

               vec<Bool> lose( nrows, False );
               for ( int i = 0; i < nrows; i++ )
               {    if ( !WT[i].first || keep[i] || follow[i] ) continue;
                    fix64_6 total = WT[i].second;
                    for ( int j = 0; j < MT.isize( ); j++ )
                    {    if ( WT[j].first || follow[j] ) continue;
                         Bool OK = True;
                         for ( int l = start[j]; l < stop[j]; l++ )
                         {    if ( MT[i][l] != MT[j][l] )
                              {    OK = False;
                                   break;    }    }
                         if (OK) total += WT[j].second;    }
                    fix64_6 control = 0;
                    for ( int j = 0; j < nrows; j++ )
                    {    if ( !WT[j].first ) continue;
                         if ( lext[j] < lext[i] || rext[j] < rext[i] ) continue;
                         control = Max( control, WT[j].second );    }
                    if ( total <= max_weak && total * min_mult <= control )
                         lose[i] = True;    }

               // Print.

               if ( logc.verb[ "WORDIFY" ] >= 3 )
               {    out << "\nMT,WT:\n";
                    for ( int i = 0; i < nrows; i++ )
                    {    out << printSeq( MT[i] ) << " [" << WT[i].second << "]";
                         if ( keep[i] ) out << " (keep)";
                         if ( follow[i] ) out << " (follow)";
                         if ( lose[i] ) out << " (lose)";
                         out << "\n";    }    }
               if ( logc.verb[ "WORDIFY" ] >= 3 )
               {    for ( int i = 0; i < nrows; i++ )
                    {    if ( lose[i] )
                         {    out << "losing ";
                              for ( int j = start[i]; j < stop[i]; j++ )
                              {    if ( j > start[i] ) out << ",";
                                   out << MT[i][j];    }
                              out << endl;    }    }    }

               // Record deletions.

               for ( int i = 0; i < nrows; i++ )
               {    if ( !lose[i] ) continue;
                    for ( int j = 0; j < mrows; j++ )
                    {    Bool OK = True;
                         for ( int l = start[i]; l < stop[i]; l++ )
                         {    if ( MT[i][l] != M[j][l] )
                              {    OK = False;
                                   break;    }    }
                         if ( OK && logc.verb[ "WORDIFY" ] >= 3 )
                              out << "should delete " << printSeq(M[j]) << endl;
                         if (OK) delp[ slocs[s][j].first ] = True;    }
                    for ( int j1 = 0; j1 < mrows; j1++ )
                    {    if ( M[j1][left+w] >= 0 ) continue;
                         for ( int j2 = 0; j2 < mrows; j2++ )
                         {    if ( M[j2][left-1] >= 0 ) continue;
                              Bool OK = True;
                              for ( int l = start[i]; l < stop[i]; l++ )
                              {    if ( ( M[j1][l] >= 0 && MT[i][l] != M[j1][l] )
                                      || ( M[j2][l] >= 0 && MT[i][l] != M[j2][l] )
                                      || ( M[j1][l] < 0 && M[j2][l] < 0 ) )
                                   {    OK = False;
                                        break;    }    }
                              if ( !OK ) continue;
                              if ( logc.verb[ "WORDIFY" ] >= 3 )
                              {    out << "should delete join of " 
                                        << printSeq(M[j1]) << " and " 
                                        << printSeq(M[j2]) << endl;    }
                              int start1 = 0, start2 = 0;
                              for ( int l = 0; l < cols; l++ )
                              {    if ( M[j1][l] >= 0 ) 
                                   {    start1 = l;
                                        break;    }   }
                              for ( int l = 0; l < cols; l++ )
                              {    if ( M[j2][l] >= 0 ) 
                                   {    start2 = l;
                                        break;    }   }
                              int offset = start2 - start1;
                              int v1 = slocs[s][j1].first, v2 = slocs[s][j2].first;
                              for ( int l = 0; l < over[v1].isize( ); l++ )
                              {    if ( over[v1][l].first != v2 ) continue;
                                   if ( over[v1][l].second != offset ) continue;
                                   over_del[v1][l] = True;
                                   // out << "deleting overlap" << endl;    
                                        }
                              for ( int l = 0; l < over[v2].isize( ); l++ )
                              {    if ( over[v2][l].first != v1 ) continue;
                                   int offsetr = -offset;
                                   if ( over[v2][l].second != offsetr ) continue;
                                   over_del[v2][l] = True;
                                   // out << "deleting overlap" << endl;    
                                        }    }    }    }    }
                              
          // Print.

          if ( logc.verb[ "WORDIFY" ] >= 1 )
          {    vec< vec<String> > rows;
               vec<String> row;
               for ( int l = 0; l < left; l++ )
                    row.push_back( "" );
               for ( int l = 0; l < segs[s].isize( ); l++ )
                    row.push_back( ToString( segs[s][l] ) );
               rows.push_back(row);
               for ( int j = 0; j < slocs[s].isize( ); j++ )
               {    vec<String> row;
                    for ( int l = 0; l < left - slocs[s][j].second; l++ )
                         row.push_back( "" );
                    int id = slocs[s][j].first;
                    for ( int l = 0; l < Path(id).isize( ); l++ )
                         row.push_back( ToString( Path(id)[l] ) );
                    int pad = right - ( Path(id).isize( ) - slocs[s][j].second - w );
                    for ( int l = 0; l < pad; l++ )
                         row.push_back( "" );
                    double w = Weight(id).ToDouble( );
                    double wr = int(round(10*w)) / 10.0;
                    String wrs = ToString(wr);
                    if ( !wrs.Contains( "." ) ) wrs = wrs + ".0";
                    row.push_back( "   " + wrs );
                    rows.push_back(row);    }
               out << "\nsegment " << " = " << printSeq( segs[s] ) 
                    << " (kmers=" << nkmers << ")" << endl;
               PrintTabular( out, rows, 2, "rrrrrrrrrr" );    }

          /*
          // =======================================================================
          // Start of an alternate approach, on hold for now.
          vec<int> lefts, rights;
          for ( int i = 0; i < M.isize( ); i++ )
          {    if ( M[i][left-1] >= 0 ) lefts.push_back( M[i][left-1] );
               if ( M[i][left+w] >= 0 ) rights.push_back( M[i][left+w] );    }
          UniqueSort(lefts), UniqueSort(rights);
          if ( logc.verb[ "WORDIFY" ] >= 1 )
          {    out << "\nlefts = " << printSeq(lefts) << endl;
               out << "rights = " << printSeq(rights) << endl;    }
          vec<fix64_6> counts;
          for ( int j1 = 0; j1 < lefts.isize( ); j1++ )
          for ( int j2 = 0; j2 < rights.isize( ); j2++ )
          {    fix64_6 count = 0;
               for ( int i = 0; i < M.isize( ); i++ )
               {    if ( M[i][left-1] == lefts[j1] && M[i][left+w] == rights[j2] )
                         count += W[i];    }
               if ( logc.verb[ "WORDIFY" ] >= 1 )
               {    out << lefts[j1] << "|" << rights[j2] << " --> " 
                         << count << endl;    }
               counts.push_back(count);    }
          Sort(counts);
          vec<fix64_6> cleft( lefts.size( ), 0 ), cright( rights.size( ), 0 );
          for ( int i = 0; i < M.isize( ); i++ )
          {    if ( M[i][left-1] >= 0 ) 
                    cleft[ BinPosition( lefts, M[i][left-1] ) ] += W[i];
               if ( M[i][left+w] >= 0 ) 
                    cright[ BinPosition( rights, M[i][left+w] ) ] += W[i];    }
          const int see_low = 1;
          const int see_high = 5;
          vec< pair<int,int> > dels;
          if ( counts.nonempty( ) && counts[0] <= see_low 
               && counts.back( ) >= see_high )
          {    Bool ok = True;
               for ( int i = 0; i < counts.isize( ); i++ )
                    if ( counts[i] > see_low && counts[i] < see_high ) ok = False;
               if (ok)
               {    for ( int j1 = 0; j1 < lefts.isize( ); j1++ )
                    for ( int j2 = 0; j2 < rights.isize( ); j2++ )
                    {    fix64_6 count = 0;
                         for ( int i = 0; i < M.isize( ); i++ )
                         {    if ( M[i][left-1] == lefts[j1] 
                                   && M[i][left+w] == rights[j2] )
                              {    count += W[i];    }    }
                         if ( count <= see_low ) dels.push( j1, j2 );    }
                    vec<Bool> to_delete( dels.size( ), False );

                    for ( int r = 0; r < lefts.isize( ); r++ )
                    {    if ( cleft[r] > see_low )
                         {    int d = 0;
                              for ( int q = 0; q < dels.isize( ); q++ )
                                   if ( dels[q].first == r ) d++;
                              if ( d == rights.isize( ) )
                              {    for ( int q = 0; q < dels.isize( ); q++ )
                                   {    if ( dels[q].first == r )
                                             to_delete[q] = True;    }    }    }    }
                    for ( int r = 0; r < rights.isize( ); r++ )
                    {    if ( cright[r] > see_low )
                         {    int d = 0;
                              for ( int q = 0; q < dels.isize( ); q++ )
                                   if ( dels[q].second == r ) d++;
                              if ( d == lefts.isize( ) )
                              {    for ( int q = 0; q < dels.isize( ); q++ )
                                   {    if ( dels[q].second == r )
                                             to_delete[q] = True;    }    }    }    }
                    EraseIf( dels, to_delete );

                    for ( int q = 0; q < dels.isize( ); q++ )
                    {    int j1 = dels[q].first, j2 = dels[q].second;
                         if ( logc.verb[ "WORDIFY" ] >= 1 )
                         {    out << "recommend deleting "
                                   << lefts[j1] << "|" << rights[j2] << endl;    }
                         for ( int i = 0; i < M.isize( ); i++ )
                         {    if ( M[i][left-1] == lefts[j1] 
                                   && M[i][left+w] == rights[j2] )
                              {    if ( logc.verb[ "WORDIFY" ] >= 1 )
                                   {    out << "xshould delete " << printSeq(M[i])
                                             << endl;    }
                                   delp[ slocs[s][i].first ] = True;    }    }
                         for ( int k1 = 0; k1 < M.isize( ); k1++ )
                         {    if ( M[k1][left+w] >= 0 ) continue;
                              if ( M[k1][left-1] != lefts[j1] ) continue;
                              for ( int k2 = 0; k2 < M.isize( ); k2++ )
                              {    if ( M[k2][left-1] >= 0 ) continue;
                                   if ( M[k2][left+w] != rights[j2] ) continue;
                                   if ( logc.verb[ "WORDIFY" ] >= 1 )
                                   {    out << "xshould delete join of " 
                                             << printSeq(M[k1]) << " and " 
                                             << printSeq(M[k2]) << endl;    }
                                   int start1 = 0, start2 = 0;
                                   for ( int l = 0; l < cols; l++ )
                                   {    if ( M[k1][l] >= 0 ) 
                                        {    start1 = l;
                                             break;    }   }
                                   for ( int l = 0; l < cols; l++ )
                                   {    if ( M[k2][l] >= 0 ) 
                                        {    start2 = l;
                                             break;    }   }
                                   int offset = start2 - start1;
                                   int v1 = slocs[s][k1].first; 
                                   int v2 = slocs[s][k2].first;
                                   for ( int l = 0; l < over[v1].isize( ); l++ )
                                   {    if ( over[v1][l].first != v2 ) continue;
                                        if ( over[v1][l].second != offset ) 
                                             continue;
                                        over_del[v1][l] = True;    }
                                   for ( int l = 0; l < over[v2].isize( ); l++ )
                                   {    if ( over[v2][l].first != v1 ) continue;
                                        int offsetr = -offset;
                                        if ( over[v2][l].second != offsetr ) 
                                             continue;
                                        over_del[v2][l] = True;    }    }    }
                                   }    }    }
          // =======================================================================
          */

          Bool print = False;
          if ( trace.empty( ) ) print = True;
          else
          {    for ( int j = 0; j < segs[s].isize( ); j++ )
                    if ( BinMember( trace, segs[s][j] ) ) print = True;    }
          if (print)
          {    
               #pragma omp critical
               {    cout << out.str( );    }    }    }

     // Forcing symmetry.  May not be necessary.

     vec<int> inv = Inv( );
     for ( int i = 0; i < w.isize( ); i++ )
     for ( int j = 0; j < over[i].isize( ); j++ )
     {    if ( !over_del[i][j] ) continue;
          if ( inv[ w[i][0] ] < 0 ) continue;
          int ip = over[i][j].first;
          vec<int> p( w[i] );
          p.ReverseMe( );
          for ( int l = 0; l < p.isize( ); l++ )
               p[l] = inv[ p[l] ];
          int i2 = BinPosition( w, p );
          ForceAssertGe( i2, 0 );
          vec<int> pp( w[ip] );
          pp.ReverseMe( );
          for ( int l = 0; l < pp.isize( ); l++ )
               pp[l] = inv[ pp[l] ];
          int ip2 = BinPosition( w, pp );
          ForceAssertGe( ip2, 0 );
          for ( int j2 = 0; j2 < over[i2].isize( ); j2++ )
          {    if ( over[i2][j2].first == ip2
                    && over[i2][j2].second == over[i][j].second )
               {    over_del[i2][j2] = True;    }    }    }

     // Given word i, consider all the words j that contain it properly.  Suppose 
     // that every legal extension of i is contained in some extension of some j.  
     // Then we delete i.
     //
     // Note that when we check for containment, we should have computed an offset
     // for this in advance, and check for containment relative to that offset.
     //
     // TURNED OFF BECAUSE TOO SLOW IN MOUSE.

     /*
     double zclock = WallClockTime( ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     for ( int i = 0; i < w.isize( ); i++ )
     {    vec< pair<int,int> > con;
          // cout << "\nassessing " << printSeq(w[i]) << endl; // XXXXXXXXXXXXXXXXXX
          for ( int u = 0; u < windex[ w[i][0] ].isize( ); u++ )
          {    int j = windex[ w[i][0] ][u].first, o = windex[ w[i][0] ][u].second;
               if ( j == i ) continue;
               // cout << "contained in " << printSeq(w[j]) << endl; // XXXXXXXXXXXX
               if ( w[j].Contains( w[i], o ) ) con.push( j, o );    }
          if ( con.empty( ) ) continue;
          vec< vec<int> > joins2;
          for ( int r = 0; r < con.isize( ); r++ )
          {    int j = con[r].first, o = con[r].second;
               for ( int n = 0; n < over[j].isize( ); n++ )
               {    int l = over[j][n].first, q = over[j][n].second;
                    vec<int> join2;
                    const vec<int> &v1 = w[j], &v2 = w[l];
                    int m1 = v1.size( ), m2 = v2.size( );
                    if ( q >= 0 )
                    {    join2.resize( Max( m1, q + m2 ) );
                         for ( int l = 0; l < m1; l++ )
                              join2[l] = v1[l];
                         for ( int l = 0; l < m2; l++ )
                              join2[q+l] = v2[l];    }
                    else
                    {    join2.resize( Max( m2, -q + m1 ) );
                         for ( int l = 0; l < m1; l++ )
                              join2[-q+l] = v1[l];
                         for ( int l = 0; l < m2; l++ )
                              join2[l] = v2[l];    }
                    // cout << "join2 is " << printSeq(join2) << endl; // XXXXXXXXXX
                    joins2.push_back(join2);    }    }
          Bool valid = True;
          for ( int m = 0; m < over[i].isize( ); m++ )
          {    if ( over_del[i][m] ) continue;
               int k = over[i][m].first, p = over[i][m].second;
               vec<int> join;
               const vec<int> &w1 = w[i], &w2 = w[k];
               // cout << "looking at overlap with " << printSeq(w[k]) << endl; //XXX
               int n1 = w1.size( ), n2 = w2.size( );
               if ( p >= 0 )
               {    join.resize( Max( n1, p + n2 ) );
                    for ( int l = 0; l < n1; l++ )
                         join[l] = w1[l];
                    for ( int l = 0; l < n2; l++ )
                         join[p+l] = w2[l];    }
               else
               {    join.resize( Max( n2, -p + n1 ) );
                    for ( int l = 0; l < n1; l++ )
                         join[-p+l] = w1[l];
                    for ( int l = 0; l < n2; l++ )
                         join[l] = w2[l];    }
               // cout << "join is " << printSeq(join) << endl; // XXXXXXXXXXXXXXXXX
               Bool found = False;
               for ( int s = 0; s < joins2.isize( ); s++ )
               {    if ( joins2[s].Contains(join) )
                    {    found = True;
                         break;    }    }
               if ( !found )
               {    valid = False;
                    break;    }    }
          if (valid)
          {    delp[i] = True;
               // cout << "killing " << printSeq(w[i]) << endl; // XXXXXXXXXXXXXXXXX
               if ( Inv( w[i][0] ) >= 0 )
               {    vec<int> wr;
                    for ( int j = w[i].isize( ) - 1; j >= 0; j-- )
                         wr.push_back( Inv( w[i][j] ) );
                    int g = BinPosition( w, wr );
                    if ( g >= 0 ) delp[g] = True;    }    }    }
     if (logc.STATUS_LOGGING)
          cout << TimeSince(zclock) << " used in word containment" << endl;
     */

     // Erase dead overlaps.

     for ( int i = 0; i < w.isize( ); i++ )
          EraseIf( over[i], over_del[i] );

     // Delete paths.

     vec<int> to_new( w.size( ), -1 );
     int ct = 0;
     for ( int i = 0; i < w.isize( ); i++ )
     {    if ( delp[i] ) continue;
          to_new[i] = ct++;    }
     for ( int i1 = 0; i1 < w.isize( ); i1++ )
     {    vec<Bool> del( over[i1].size( ), False );
          for ( int j = 0; j < over[i1].isize( ); j++ )
          {    int& i2 = over[i1][j].first;
               if ( to_new[i2] < 0 ) del[j] = True;
               else i2 = to_new[i2];    }
          EraseIf( over[i1], del );    }
     EraseIf( w, delp );
     EraseIf( over, delp );
     // cout << "\nWORDS:\n"; // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     // for ( int i = 0; i < w.isize( ); i++ ) // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     //      cout << "[" << i << "] " << printSeq(w[i]) << endl; // XXXXXXXXXXXXXXXX

     // Define overx.

     vec< vec< pair<int,int> > > overx( over.size( ) );
     for ( int i = 0; i < w.isize( ); i++ )
     {    overx[i].resize( over[i].size( ) );
          for ( int j = 0; j < over[i].isize( ); j++ )
          {    overx[i][j].first = over[i][j].first;
               overx[i][j].second = over[i][j].second;    }    }

     // Delete subsumption overlaps.

     /*
     for ( int i1 = 0; i1 < overx.isize( ); i1++ )
     {    vec<Bool> to_delete( overx[i1].size( ), True );
          for ( int j = 0; j < overx[i1].isize( ); j++ )
          {    int i2 = overx[i1][j].first;
               int o = overx[i1][j].second;
               const vec<int> &w1 = w[i1], &w2 = w[i2];
               int n1 = w1.size( ), n2 = w2.size( );
               if ( o >= 0 && o + n2 <= n1 ) continue;
               if ( o <= 0 && o + n2 >= n1 ) continue;
               to_delete[j] = False;    }
          EraseIf( overx[i1], to_delete );    }
     */
     
     if ( logc.verb[ "WORDIFY" ] >= 4 )
     {    cout << "\nOVERLAPS:\n";
          for ( int i1 = 0; i1 < w.isize( ); i1++ )
          for ( int j = 0; j < overx[i1].isize( ); j++ )
          {    int i2 = overx[i1][j].first;
               int o = overx[i1][j].second;
               const vec<int> &w1 = w[i1], &w2 = w[i2];
               int n1 = w1.size( ), n2 = w2.size( );
               // Don't print overlaps that are "subsumptions".
               if ( o >= 0 && o + n2 <= n1 ) continue;
               if ( o <= 0 && o + n2 >= n1 ) continue;
               vec<int> w12;
               if ( o >= 0 )
               {    w12.resize( o + n2 );
                    for ( int l = 0; l < n1; l++ )
                         w12[l] = w1[l];
                    for ( int l = 0; l < n2; l++ )
                         w12[o+l] = w2[l];    }
               else
               {    w12.resize( -o + n1 );
                    for ( int l = 0; l < n1; l++ )
                         w12[-o+l] = w1[l];
                    for ( int l = 0; l < n2; l++ )
                         w12[l] = w2[l];    }
               cout << printSeq(w[i1]) << " with " << printSeq(w[i2]) 
                    << " yields " << printSeq(w12) << endl;    }    }

     // Create digraph.

     digraphE< vec<int> > D;
     WordsToDigraphAlt2( w, overx, D, inv, log_control, logc );

     // Print.

     if ( logc.verb[ "WORDIFY" ] >= 2 )
     {    cout << "\ndigraph:\n";
          digraphE< vec<int> > DX(D);
          DX.RemoveUnneededVertices( );
          int count = 0;
          for ( int v = 0; v < DX.N( ); v++ )
          for ( int j = 0; j < DX.From(v).isize( ); j++ )
          {    int w = DX.From(v)[j];
               const vec<int>& e = DX.EdgeObjectByIndexFrom( v, j );
               cout << "[" << count++ << "] " << v
                    << "---(" << printSeq(e) << ")---> " << w << "\n";    }
          cout << "\n";    }

     // Make new graph.

     vec<basevector> xedges( D.EdgeObjectCount( ) );
     for ( int e = 0; e < D.EdgeObjectCount( ); e++ )
     {    const vec<int>& x = D.EdgeObject(e);
          xedges[e] = EdgeObject( x[0] );
          for ( int j = 1; j < x.isize( ); j++ )
               xedges[e] = TrimCat( K( ), xedges[e], EdgeObject( x[j] ) );    }
     int edges_orig = EdgeObjectCount( );
     HyperBasevector hb( K( ), D, xedges );
     (HyperBasevector&)(*this) = hb;
     InvMutable( ) = inv;

     // Translate paths.

     vec<int> to_left, to_right;
     D.ToLeft(to_left), D.ToRight(to_right);
     vec< vec<int> > dindex(edges_orig);
     for ( int j = 0; j < D.EdgeObjectCount( ); j++ )
          dindex[ D.EdgeObject(j)[0] ].push_back(j);

     vec< vec<int> > pathsx;
     vec<fix64_6> weightsx_fw, weightsx_rc;
     if ( logc.verb[ "WORDIFY" ] >= 1 )
          cout << "\ntranslating " << NPaths( ) << " paths" << endl;
     #pragma omp parallel for
     for ( int z = 0; z < NPaths( ); z++ )
     {    const vec<int>& w = Path(z);
          vec< vec<int> > paths;
          Snorkle( w, dindex, to_left, to_right, paths );
          if ( logc.verb[ "WORDIFY" ] >= 3 )
          {
               #pragma omp critical
               {    cout << "\n" << Date( ) << ": z = " << z << ", Path(z) = " 
                         << printSeq(w) << endl << paths.size( ) << " paths found" 
                         << endl;    }    }
          #pragma omp critical
          {    pathsx.append(paths);
               for ( int j = 0; j < paths.isize( ); j++ )
               {    weightsx_fw.push_back( WeightFw(z) / paths.size( ) );
                    weightsx_rc.push_back( 
                         WeightRc(z) / paths.size( ) );    }    }    }
     PathsMutable( ) = pathsx; 
     WeightsFwMutable( ) = weightsx_fw;
     WeightsRcMutable( ) = weightsx_rc;
     SortSync( PathsMutable( ), WeightsFwMutable( ), WeightsRcMutable( ) );
     
     vec< pair< vec<int>, vec<int> > > pairsx;
     vec< vec<pair_point> > pairdatax;
     #pragma omp parallel for
     for ( int z = 0; z < NPairs( ); z++ )
     {    const vec<int> &w1 = PairLeft(z), &w2 = PairRight(z);
          vec< vec<int> > pairs1, pairs2;
          Snorkle( w1, dindex, to_left, to_right, pairs1 );
          Snorkle( w2, dindex, to_left, to_right, pairs2 );
          vec< pair< vec<int>, vec<int> > > pairs;
          vec<pair_point> pairs_data;
          for ( int j1 = 0; j1 < pairs1.isize( ); j1++ )
          for ( int j2 = 0; j2 < pairs2.isize( ); j2++ )
               pairs.push( pairs1[j1], pairs2[j2] );
          int n = pairs1.size( ) * pairs2.size( );
          vec<pair_point> pd = PairData(z);

          // Test for n > 0 to avoid dividing by zero.  Perhaps the n = 0 case is
          // indicative of a bug.

          if ( n > 0 )
          {    for ( int i = 0; i < pd.isize( ); i++ )
                    pd[i].WeightMutable( ) /= n;    }
          #pragma omp critical
          {    pairsx.append(pairs);
               for ( int i = 0; i < n; i++ )
                    pairdatax.push_back(pd);    }    }
     PairsMutable( ) = pairsx, PairDataMutable( ) = pairdatax;
     SortSync( PairsMutable( ), PairDataMutable( ) );

     // Eliminate duplicate edges.  I'm not sure how this can happen.

     vec<int> trans( EdgeObjectCount( ), vec<int>::IDENTITY );
     for ( int v = 0; v < N( ); v++ )
     {    for ( int j1 = 0; j1 < From(v).isize( ); j1++ )
          for ( int j2 = j1 + 1; j2 < From(v).isize( ); j2++ )
          {    if ( From(v)[j1] != From(v)[j2] ) continue;
               int e1 = EdgeObjectIndexByIndexFrom( v, j1 );
               int e2 = EdgeObjectIndexByIndexFrom( v, j2 );
               if ( EdgeObject(e1) != EdgeObject(e2) ) continue;
               if ( trans[e1] != e1 || trans[e2] != e2 ) continue;
               trans[e2] = e1;
               // If the following fails, need to tweak involution.
               ForceAssert( e2 != Inv(e1) );
               if ( Inv(e2) >= 0 ) trans[ Inv(e2) ] = Inv(e1);    }    }
     vec<int> dels;
     for ( int e = 0; e < EdgeObjectCount( ); e++ )
          if ( trans[e] != e ) dels.push_back(e);
     for ( int i = 0; i < NPaths( ); i++ )
     {    for ( int j = 0; j < Path(i).isize( ); j++ )
               PathMutable(i,j) = trans[ Path(i,j) ];    }
     for ( int i = 0; i < NPairs( ); i++ )
     {    for ( int j = 0; j < PairLeft(i).isize( ); j++ )
               PairLeftMutable(i,j) = trans[ PairLeft(i,j) ];
          for ( int j = 0; j < PairRight(i).isize( ); j++ )
               PairRightMutable(i,j) = trans[ PairRight(i,j) ];    }
     DeleteEdges(dels);

     // Clean up.

     REPORT_TIME( clock, "used in WordifyAlt2" );
     RemoveUnneededVertices( );
     RemoveDeadEdgeObjects( );
     double clock5 = WallClockTime( );
     UniqueOrderPaths( );
     REPORT_TIME( clock5, "used in WordifyAlt2 - tail2" );

     // Force symmetry of weights.  This is due to small arithmetic differences.
     // Might be better to use rational numbers to represent floats.

     FixWeights(logc);    }
