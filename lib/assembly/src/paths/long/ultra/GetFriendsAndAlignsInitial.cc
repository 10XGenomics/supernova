///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "CoreTools.h"
#include "kmers/KmerRecord.h"
#include "math/Combinatorics.h"
#include "math/Functions.h"
#include "pairwise_aligners/SmithWatBandedA.h"
#include "paths/long/ultra/ConsensusScoreModel.h"
#include "paths/long/Friends.h"
#include "paths/long/KmerAlign.h"
#include "paths/long/Logging.h"
#include "paths/long/ultra/FounderAlignment.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/MakeKmerStuff.h"
#include "random/Bernoulli.h"

#include "PrintAlignment.h"

#include <omp.h>
// MakeDepend: library OMP

namespace { // open anonymous namespace

template<int K> void MakeKmerLookup0Single( const vecbasevector& unibases,
     vec< triple<kmer<K>,int,int> >& kmers_plus, 
     const long_logging_control& log_control, const long_logging& logc, ostream& out )
{    double clock1 = WallClockTime( );
     vec<int64_t> starts;
     starts.push_back(0);
     for ( size_t i = 0; i < unibases.size( ); i++ )
     {    const basevector& u = unibases[i];
          starts.push_back( starts.back( ) + Max( 0, u.isize( ) - K + 1 ) );    }
     kmers_plus.resize( starts.back( ) );
     for ( size_t i = 0; i < unibases.size( ); i++ )
     {    const basevector& u = unibases[i];
          for ( int j = 0; j <= u.isize( ) - K; j++ )
          {    int64_t r = starts[i] + j;
               kmers_plus[r].first.SetToSubOf( u, j );
               kmers_plus[r].second = i; 
               kmers_plus[r].third = j;    }    }
     REPORT_TIMEX( clock1, "used in kmer pile creation", out );
     double clock2 = WallClockTime( );
     Sort(kmers_plus);
     REPORT_TIMEX( clock2, "used in kmer pile sorting", out );    }

} // close anonymous namespace

template<int K> void GetFriendsAndAlignsInitial( 
     const IAndOsVec& F,                   // given friends
     const vecbasevector& reads,            // reads
     int ec,                                // where founder was
     vecbasevector& gang,                   // the gang; read 0 = founder
     vec<int> * p_friends,               // optionally output friend ids in the gang
     vec< vec< pair<int,int> > >& a,        // kmer aligns of gang to read 0
     ostream& out,                          // for logging
     const ConsensusScoreModel& error_model,
     const long_heuristics& heur,
     const long_logging_control& log_control, const long_logging& logc )
{
     // Use given friends.  In this initial version we go ahead and run the
     // rest of the code as before.  This is unlikely to be ideal.

     vec< triple<kmer<K>,int,int> > kmers_plus;
     vec<int> friend_ids;
     int ec_orig = ec;
     double gclock = WallClockTime( );
     friend_ids.push_back(ec);
     gang.resize(0);
     gang.reserve( F[ec].size( ) + 1 );
     gang.push_back( reads[ec] );
     if ( F[ec].empty( ) && logc.verb[ "ULTRA" ] >= 1 ) 
          out << "Warning: this read has no friends!\n";
     for ( int i = 0; i < (int) F[ec].size( ); i++ )
     {    const IdAndOrientation& x = F[ec][i];
          int id = x.getId( );
          friend_ids.push_back(id);
          if ( !x.isRC( ) ) gang.push_back(reads[id]);
          else
          {    basevector b = reads[id];
               b.ReverseComplement( );
               gang.push_back(b);    }    }
     REPORT_TIMEX( gclock, "used defining gang", out );
     MakeKmerLookup0Single( gang, kmers_plus, log_control, logc, out );

     // Compute kmer matches between reads and read ec (which is now read 0).

     int N = gang.size( );
     double match_clock = WallClockTime( );
     vec< vec< pair<int,int> > > offsets(N);
     for ( int64_t i = 0; i < (int64_t) kmers_plus.size( ); i++ )
     {    int64_t j, z1, z2;
          for ( j = i + 1; j < (int64_t) kmers_plus.size( ); j++ )
               if ( kmers_plus[j].first != kmers_plus[i].first ) break;
          for ( z1 = i; z1 < j; z1++ )
               if ( kmers_plus[z1].second == 0 ) break;
          for ( z2 = z1; z2 < j; z2++ )
               if ( kmers_plus[z2].second != 0 ) break;
          for ( int64_t z = z1; z < z2; z++ )
          {    int pos1 = kmers_plus[z].third;
               for ( int64_t k = i; k < j; k++ )
               {    int id2 = kmers_plus[k].second; 
                    if ( id2 == 0 ) continue;
                    int pos2 = kmers_plus[k].third; // position on read id2
                         offsets[id2].push( pos1-pos2, pos1 );    }    }
          i = j - 1;    }
     for ( int id = 1; id < N; id++ )
          Sort( offsets[id] );
     if ( logc.PRINT_TIME_USED && omp_get_thread_num( ) == 0 )
          out << TimeSince(match_clock) << " used computing kmer matches" << endl;

     // Find the maximum number mh of hits.  For each read achieving at least
     // mh/2 hits, find the median offset.  Do a banded Smith-Waterman versus
     // read ec.  Build kmer aligns to read ec.

     double align_clock = WallClockTime( );
     vec<align> aligns( gang.size( ) );
     vec<Bool> accepted( N, False );
     vec<int> errs(N, -1);
     vec<ho_interval> ext1(N);
     {    int mh = 0;
          for ( int id = 0; id < N; id++ )
               if ( id != 0 ) mh = Max( mh, offsets[id].isize( ) );
          vec<int> pos1_count( N );
          for ( int id = 0; id < N; id++ )
          {    if ( id == 0 ) continue;
               vec<int> pos1;
               for ( int j = 0; j < offsets[id].isize( ); j++ )
                    pos1.push_back( offsets[id][j].second );
               UniqueSort(pos1);
               pos1_count[id] = pos1.size( );    }
          int max_pos1_count = 0;
          for ( int id = 0; id < N; id++ )
          {    if ( id != 0 ) 
                    max_pos1_count = Max( max_pos1_count, pos1_count[id] );    }
          if ( logc.verb[ "FRIEND" ] >= 1 ) PRINT(mh);
          a.clear_and_resize(N);
          for ( int id = 0; id < N; id++ )
          {    if ( id == 0 ) continue;
               if ( pos1_count[id] >= max_pos1_count/2 && offsets[id].nonempty( ) )
               {    accepted[id] = True;
                    // Using kmer-align method instead of smith-waterman
                    if (heur.USE_KMER_ALIGN_METHOD) 
                    {    KmerAlign( offsets[id], a[id], 
                              logc.verb[ "FRIEND" ] ); 
                         continue;    }   

                    // Form offsets into groups, breaking whenever there is a
                    // > max_sep = 20 separation in offsets.

                    vec<int> ostarts;
                    ostarts.push_back(0);
                    const int max_sep = 20;
                    for ( int j = 0; j < offsets[id].isize( ) - 1; j++ )
                    {    if ( offsets[id][j+1].first > offsets[id][j].first + 20 )
                              ostarts.push_back(j+1);    }
                    ostarts.push_back( offsets[id].size( ) );

                    // Compute the size of each group, as measured by its number
                    // of distinct rpos2 values.  Find the largest group.

                    int gp_max = 0, gp_best = -1;
                    for ( int j = 0; j < ostarts.isize( ) - 1; j++ )
                    {    vec<int> rpos2;
                         for ( int l = ostarts[j]; l < ostarts[j+1]; l++ )
                              rpos2.push_back( offsets[id][l].second );
                         UniqueSort(rpos2);
                         if ( rpos2.isize( ) > gp_max )
                         {    gp_max = rpos2.size( );
                              gp_best = j;    }    }

                    // Set offset to the median within the winning group, and
                    // set bandwidth to span the group, but no more than 
                    // max_bandwidth = 200;

                    int offset, bandwidth;
                    if ( !heur.NEW_SMITH_WAT )
                    {    bandwidth = 200;
                         offset = offsets[id][ offsets[id].isize( )/2 ].first;    }
                    else
                    {    int start = ostarts[gp_best], stop = ostarts[gp_best+1];
                         int mid = start + (stop-start)/2;
                         offset = offsets[id][mid].first;
                         const int bw_add = 12;
                         bandwidth = Max( offset - offsets[id][start].first,
                              offsets[id][stop-1].first - offset ) + bw_add;
                         const int max_bandwidth = 200;
                         bandwidth = Min( bandwidth, max_bandwidth );    }

                    if ( bandwidth < 0 ) // SHOULD NOT HAPPEN!!!!!!!!!!!!!!!!!!!!!!!
                         bandwidth = 0;

                    // Align.

                    align x;
                    int errors;
                    SmithWatBandedA( gang[0], gang[id], offset, bandwidth, x,
                         errors, 0, 1, 1 );    
                    errs[id] = errors;
                    ext1[id] = x.Extent1( );
                    aligns[id] = x;
                    vec<ho_interval> p1, p2;
                    x.PerfectIntervals1( gang[0], gang[id], p1 );
                    x.PerfectIntervals2( gang[0], gang[id], p2 );    
                    for ( int j = 0; j < p1.isize( ); j++ )
                    {    const ho_interval &h1 = p1[j], &h2 = p2[j];
                         for ( int l = h1.Start( ); l <= h1.Stop( ) - K; l++ )
                         {    a[id].push( l + h2.Start( ) - h1.Start( ), 
                                   l );    }    }    }    }    }
     for ( int p = 0; p <= gang[0].isize( ) - K; p++ )
          a[0].push( p, p );
     if ( logc.PRINT_TIME_USED && omp_get_thread_num( ) == 0 )
          out << TimeSince(align_clock) << " used aligning" << endl;

     // Filter out reads based on error rate.  Note that this doesn't actually
     // remove the reads yet (which might be better).

     double err_clock = WallClockTime( );
     vec<Bool> evil( gang.size( ), False );
     if (heur.ERR_FILTER)
     {    vec<Bool> e_to_delete( N, False );
          vec<int> eids( N, vec<int>::IDENTITY );
          for ( int i = 0; i < N; i++ )
               if ( errs[i] < 0 ) e_to_delete[i] = True;
          EraseIf( errs, e_to_delete );
          EraseIf( ext1, e_to_delete );
          EraseIf( eids, e_to_delete );
          vec<double> erate( errs.size( ) );
          for ( int j = 0; j < errs.isize( ); j++ )
               erate[j] = double(errs[j]) / double( ext1[j].Length( ) );
          SortSync( erate, errs, ext1, eids );

          // Walk through errs until the founder is covered.

          vec<ho_interval> cov;
          int ep;
          for ( ep = 0; ep < ext1.isize( ); ep++ )
          {    cov.push_back( ext1[ep] );
               if ( TotalCovered(cov) == gang[0].isize( ) ) break;    }
          if ( ep < ext1.isize( ) )
          {    
               // Filter out reads that have too high an error rate.

               const double dev_mult = 3.0;
               double add = dev_mult * erate[ep] * sqrt(errs[ep]) / errs[ep];
               for ( int j = 0; j < errs.isize( ); j++ )
               {    if ( erate[j] > erate[ep] + add )
                         evil[ eids[j] ] = True;    }    }    }
     if ( logc.PRINT_TIME_USED && omp_get_thread_num( ) == 0 )
          out << TimeSince(err_clock) << " used filtering based on errors" << endl;

     // Filter out reads that appear not to belong.

     double cclock = WallClockTime( );
     vec<double> evil_score( gang.size( ), 1 );
     VecUCharVec multi;
     vec<Bool> funkyc( gang[0].size( ), False );
     if ( heur.COLUMN_FILTER && !heur.USE_KMER_ALIGN_METHOD )
     {
          // Create substitution-only multiple alignment.

          unsigned char EMPTY_CELL = 5;
          int nrows = gang.size( ), ncols = gang[0].size();
          if ( !heur.COLUMN_FILTER_SUBST_ONLY )
          {    Scorer scorer( error_model.GetSubRate( ),
                    error_model.GetDelRate( ), error_model.GetInsRate( ) );
               AlignFriendsToFounder( gang, 0, aligns, scorer, &multi );
               ncols = multi[0].size();
               funkyc.resize(ncols, False);    }
          else
          {
              multi.resize(nrows);
              for ( int id = 0; id < nrows; id++ )
              {    multi[id].resize( ncols, EMPTY_CELL );
                   if ( id == 0 )
                   {    for ( int j = 0; j < ncols; j++ )
                             multi[id][j] = gang[0][j];    }
                   else
                   {    const align& a = aligns[id];
                        int p1 = a.pos1( ), p2 = a.pos2( );
                        for ( int j = 0; j < a.Nblocks( ); j++ )
                        {    if ( a.Gaps(j) > 0 ) p2 += a.Gaps(j);
                             if ( a.Gaps(j) < 0 ) p1 -= a.Gaps(j);
                             for ( int x = 0; x < a.Lengths(j); x++ )
                             {    multi[id][p1] = gang[id][p2];
                                  ++p1; ++p2;    }    }    }    }
          }

          // Look for funky columns.
     
          // For haploid genome, we may want the threshold to be strict so that
          // to reduce false negatives. However it is necessary to increase the
          // threshold for polymorphic genome.
          // double eps = 1.0e-40;

          const double eps = 1.0e-20;
          VecUCharVec F(nrows);
          for ( int c = 0; c < ncols; c++ )
          {    Bool funky = True;
               vec<int> x(5, 0); // (A,C,G,T,-)
               for ( int r = 0; r < nrows; r++ )
                    if ( multi[r][c] != EMPTY_CELL ) x[ multi[r][c] ]++;
               if ( logc.verb[ "FRIEND" ] >= 2 )
               {    out << "\n";
                    PRINT6_TO( out, c, x[0], x[1], x[2], x[3], x[4] );    }
               int N = Sum(x);
               double worst = 0.0;
               for ( int j = 0; j < 5; j++ )
               {    Bool funkyj = False;
                    double w = 1.0;
                    for ( int k = 0; k < 5; k++ )
                    {    int n = x[k];
                         if ( k == j || n == 0 ) continue;
                         double p;
                         if ( j < 4 )
                         {    if ( k < 4 ) p = error_model.GetSubRate( )/3.0;
                              else p = error_model.GetDelRate( );    }
                         else p = error_model.GetInsRate( )/3.0;
                         double q = BinomialProbCum( 1.0-p, N, N-n );
                         if ( q < w ) w = q;
                         if ( q <= eps )
                         {    funkyj = True;
                              // break;    
                              }    }
                    if ( w > worst ) worst = w;
                    if ( !funkyj )
                    {    funky = False;
                         break;    }    }
               if ( !funky ) continue;
               funkyc[c] = True;
               if ( logc.verb[ "FRIEND" ] >= 1 )
               {    out << "column " << c << " looks funky, score = " 
                         << worst << "\n";    }
               for ( int r = 0; r < nrows; r++ )
                    F[r].push_back( multi[r][c] );    }

          // Look for reads that appear not to belong.

          const double pd = 0.05;
          const double eps2 = 0.000001; // 10^-6
          for ( int j = 0; j < nrows; j++ )
          {    if ( j == 0 ) continue;
               if ( !accepted[j] ) continue;
               int total = 0, diffs = 0;
               for ( size_t c = 0; c < F[j].size( ); c++ )
               {    if ( F[j][c] == EMPTY_CELL ) continue;
                    total++;
                    if ( F[j][c] != F[0][c] ) diffs++;    }
               if ( diffs == 0 ) continue;
               double q = 1.0 - BinomialSum( total, diffs-1, pd );
               evil_score[j] = q;
               if ( q <= eps2 )
               {    evil[j] = True;    }    }    }
     if ( logc.PRINT_TIME_USED && omp_get_thread_num( ) == 0 )
          out << TimeSince(cclock) << " used in column filtering" << endl;

     // Delete reads for which number of alignments is less than half the max.
     // Print gang.

     double delete_clock = WallClockTime( );
     int M = 0;
     for ( int id = 1; id < N; id++ )
          M = Max( M, a[id].isize( ) );
     vec<Bool> to_delete( N, False );
     vec< vec<String> > gang_rows;
     vec< pair<int,int> > gang_start;
     vec<int> colsp;
     const vec<ref_loc>& locs = *log_control.readlocs;
     if (logc.PRINT_GANG_LOCS)
     {    ParseIntSet( logc.COLS_TO_PRINT, colsp );
          out << "\ngang locations for read " << ec_orig << ":\n";
          if ( !multi.empty() && colsp.nonempty( ) && heur.COLUMN_FILTER )
          {    vec<String> row;
               String s;
               for ( int j = 0; j < colsp.isize( ); j++ )
               {    if ( funkyc[ colsp[j] ] ) s += "X";
                    else s += ".";   }
               row.push_back( "", "", "", "", s );
               gang_rows.push_back(row);
               gang_start.push( -1, -1 );    }
          int fid = friend_ids[0];
          vec<String> row;
          row.push_back( "[0]", ToString(fid) );
          row.push_back( ToString(locs[fid].id) + "." + ToString(locs[fid].start)
               + "-" + ToString(locs[fid].stop) );
          row.push_back( "***** FOUNDER *****" );
          if ( !multi.empty() && colsp.nonempty( ) )
          {    String s;
               for ( int j = 0; j < colsp.isize( ); j++ )
               {    ForceAssert( colsp[j] >= 0 && colsp[j] < (int)multi[0].size( ) );
                    s += multi[0][ colsp[j] ];    }
               row.push_back(s);    }
          gang_rows.push_back(row);
          gang_start.push( locs[fid].id, locs[fid].start );    }
     int count = 0, true_friends = 0, false_friends = 0;
     for ( int id = 1; id < N; id++ )
     {    if ( a[id].isize( ) < M/2 || a[id].empty( ) /* || evil[id] */ ) 
               to_delete[id] = True;
          else if (logc.PRINT_GANG_LOCS)
          {    int fid = friend_ids[id];
               vec<String> row;
               row.push_back( "[" + ToString(++count) + "]" );
               row.push_back( ToString(fid) );
               row.push_back( ToString(locs[fid].id) + "." 
                    + ToString(locs[fid].start)
                    + "-" + ToString(locs[fid].stop) );
               if ( !evil[id] )
               {    if ( locs[fid].id == locs[ friend_ids[0] ].id
                         && IntervalOverlap( locs[fid].start, locs[fid].stop,
                         (int) locs[ friend_ids[0] ].start, 
                         (int) locs[ friend_ids[0] ].stop ) > 0 )
                    {    true_friends++;    }
                    else false_friends++;    }
               if (heur.COLUMN_FILTER)
               {    if ( evil[id] ) row.push_back( "(evil)" );
                    else row.push_back( "(" + ToString( evil_score[id] ) + ")" );   }
               if ( !multi.empty() && colsp.nonempty( ) )
               {    String s;
                    for ( int j = 0; j < colsp.isize( ); j++ )
                    {    ForceAssert( colsp[j] >= 0 
                              && colsp[j] < (int)multi[id].size( ) );
                         s += multi[id][ colsp[j] ];    }
                    row.push_back(s);    }
               gang_rows.push_back(row);
               gang_start.push( locs[fid].id, locs[fid].start );    }
          if ( evil[id] ) to_delete[id] = True;    }
     if (logc.PRINT_GANG_LOCS) 
     {    SortSync( gang_start, gang_rows );
          PrintTabular( out, gang_rows, 2, "lrll" );
          out << "\nfound " << true_friends << " true friends" 
               << " and " << false_friends << " false friends" << endl;    }
     gang.EraseIf(to_delete);
     EraseIf( a, to_delete );
     EraseIf( friend_ids, to_delete );
     if ( p_friends != 0 )  *p_friends = friend_ids;
     if ( logc.PRINT_TIME_USED && omp_get_thread_num( ) == 0 )
          out << TimeSince(delete_clock) << " used deleting reads" << endl;    }

template void GetFriendsAndAlignsInitial<16>( const IAndOsVec& F,
     const vecbasevector& reads, int ec, vecbasevector& gang, vec<int> *p_friend_ids,
     vec< vec< pair<int,int> > >& a, ostream& out, const ConsensusScoreModel&,
     const long_heuristics& heur, const long_logging_control&,
     const long_logging& logc );

template void GetFriendsAndAlignsInitial<20>( const IAndOsVec& F,
     const vecbasevector& reads, int ec, vecbasevector& gang, vec<int> *p_friend_ids,
     vec< vec< pair<int,int> > >& a, ostream& out, const ConsensusScoreModel&,
     const long_heuristics& heur, const long_logging_control&,
     const long_logging& logc );

template void GetFriendsAndAlignsInitial<24>( const IAndOsVec& F,
     const vecbasevector& reads, int ec, vecbasevector& gang, vec<int> *p_friend_ids,
     vec< vec< pair<int,int> > >& a, ostream& out, const ConsensusScoreModel&,
     const long_heuristics& heur, const long_logging_control&,
     const long_logging& logc );

template void GetFriendsAndAlignsInitial<28>( const IAndOsVec& F,
     const vecbasevector& reads, int ec, vecbasevector& gang, vec<int> *p_friend_ids,
     vec< vec< pair<int,int> > >& a, ostream& out, const ConsensusScoreModel&,
     const long_heuristics& heur, const long_logging_control&,
     const long_logging& logc );
