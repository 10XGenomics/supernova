///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>

#include "Basevector.h"
#include "CoreTools.h"
#include "FastIfstream.h"
#include "PairsManager.h"
#include "ParseSet.h"
#include "ParallelVecUtilities.h"
#include "TokenizeString.h"
#include "efasta/EfastaTools.h"
#include "kmers/KmerRecord.h"
#include "reporting/PerfStat.h"
#include "paths/long/CreateGenome.h"
#include "paths/long/Friends.h"
#include "paths/long/ultra/GetFriendsAndAlignsInitial.h"
#include "paths/long/LargeKDispatcher.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/MakeKmerStuff.h"
#include "paths/long/SupportedHyperBasevector.h"
#include "paths/long/DiscovarTools.h"
#include "paths/long/ultra/ConsensusScoreModel.h"
#include "paths/long/ultra/MakeBlocks.h"
#include "paths/long/ultra/Prefab.h"
#include "paths/long/ultra/ThreadedBlocks.h"
#include "random/Shuffle.h"
#include "reporting/PerfStat.h"
#include "system/ParsedArgs.h"

template<int K> void BuildCorrectedReads( const vecbasevector& reads,
     const IAndOsVec& F, const vec<int>& rid, VecEFasta& corrected,
     vec<int>& cid, const ConsensusScoreModel& error_model,
     const long_heuristics& heur, const long_logging_control& log_control,
     const long_logging& logc, const ref_data& ref, const int NUM_THREADS )
{
     int N = reads.size( );
     vec<String> reports(N);

     int npasses = ( heur.PREFAB ? 2 : 1 );
     vec<int> rid1;
     vec< vec<int> > friend_ids; // also record actual friends found during pass1
     for ( int pass = 1; pass <= npasses; pass++ )
     {
          // Do prefab if requested.

          vec<int> ridx;
          if ( heur.PREFAB ) 
          {    if ( pass == 1 ) 
               {    double clock = WallClockTime( );
                    SelectInitialReads( rid, F, rid1, heur, log_control );
                    cout << Date() <<  ": Selected  " << rid1.size() 
                         << " reads for correction" << endl;
                    ridx = rid1;    
                    friend_ids.resize( ridx.size() );
                    REPORT_TIME( clock, "used in initial read selection" );    }
               if ( pass == 2 )
               {    double clock = WallClockTime( );
                    VecEFasta corrected2;
                    vec<int> cid2;
                    size_t n_pre_corrected = corrected.size();
                    if ( n_pre_corrected * 10 < reads.size() ) 
                    {    cout << "Warning: than 10% of the reads are pre-corrected" 
                              << endl;    }
                    CorrectSomeMoreReads( reads, rid, rid1, friend_ids, F, corrected,
                         cid, corrected2, cid2, heur, log_control, logc, ref );
                    vec<int> processed(rid1);
                    processed.append(cid2);
                    UniqueSort(processed);
                    for ( int i = 0; i < rid.isize( ); i++ )
                    {    if ( !BinMember( processed, rid[i] ) )
                              ridx.push_back( rid[i] );    }
                    corrected.append(corrected2.begin(),corrected2.end());
                    cid.append(cid2);    
                    REPORT_TIME( clock, "used correcting more reads" );    }    }
          else ridx = rid;

          // Define batches.

          double bclock = WallClockTime( );
#if 0
// here's a completely deterministic version
          typedef std::pair<unsigned,unsigned> XId; // friendCount, ridx index
          vec<XId> xIds;
          typedef vec<XId>::iterator XItr;
          xIds.reserve(ridx.size());
          size_t totalCount = 0;
          for ( size_t idx = 0; idx != ridx.size(); ++idx )
          {
              unsigned friendCount = F[ridx[idx]].size();
              totalCount += friendCount;
              xIds.push_back( XId(friendCount,idx) );
          }
          std::sort(xIds.begin(),xIds.end(),std::greater<XId>());
          Shuffle(xIds.begin(),xIds.begin()+xIds.size()/2);

          size_t const N_BATCHES = 1000;
          size_t batchSize = totalCount / N_BATCHES;
          vec<int> batch;
          vec< vec<int> > batches;
          batches.reserve(N_BATCHES);
          XItr itr(xIds.begin());
          XItr end(xIds.end());
          while ( itr != end )
          {
              size_t curCount = 0;
              batch.clear();
              do
              {
                  batch.push_back(itr->second);
                  curCount += itr->first;
              }
              while ( ++itr != end && curCount < batchSize );
              batches.push_back(batch);
          }
#else
          vec<int> ids( ridx.size( ), vec<int>::IDENTITY );
          vec<int> friend_count( ridx.size( ) );
          for ( int xi = 0; xi < ridx.isize( ); xi++ )
          {    int ec = ridx[xi];
               friend_count[xi] = F[ec].size( );    }
          ReverseSortSync( friend_count, ids );
          vec<int> zz( ridx.size( ), vec<int>::IDENTITY );
          random_shuffle( 
               zz.begin( ), zz.begin( ) + ( zz.end( ) - zz.begin( ) ) / 2 );
          vec<int> friend_count2( ridx.size( ) ), ids2( ridx.size( ) );
          for ( int i = 0; i < ridx.isize( ); i++ )
          {    friend_count2[i] = friend_count[ zz[i] ];
               ids2[i] = ids[ zz[i] ];    }
          friend_count = friend_count2;
          ids = ids2;
          int64_t total_friends = Sum(friend_count);
          int64_t target_batches = 1000;
          int64_t target_friends = total_friends / target_batches;
          // PRINT(target_friends); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          vec< vec<int> > batches;
          for ( int j = 0; j < ids.isize( ); j++ )
          {    int f = 0, k;
               for ( k = j + 1; k < ids.isize( ); k++ )
               {    if ( f >= target_friends ) break;
                    f += friend_count[k];    }
               vec<int> b;
               for ( int l = j; l < k; l++ )
                    b.push_back( ids[l] );
               batches.push_back(b);
               j = k - 1;    }
#endif
          REPORT_TIME( bclock, "used defining batches" );
          /*
          if ( pass == 1 ) // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          {    for ( int j = 0; j < ids.isize( ); j++ ) // XXXXXXXXXXXXXXXXXXXXXXXXX
                    PRINT2( j, friend_count[j] ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
               for ( int j = 0; j < batches.isize( ); j++ ) // XXXXXXXXXXXXXXXXXXXXX
                    PRINT2( j, batches[j].size( ) ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXX
                    } // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
               // Scram();    } // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          */

          // Process reads.

          double mclock = WallClockTime( );
          cout << Date( ) << ": starting main loop" << endl;
          int count = 0;
          #pragma omp parallel for if(!logc.DIRECT) schedule(dynamic, 1)
          for ( int bi = 0; bi < batches.isize( ); bi++ )
          {    
               double biclock = -1; 
               if (logc.BATCH_TRACKING)
               {    biclock = WallClockTime( );
                    #pragma omp critical
                    {    cout << "\n" << Date( ) << ": pass = " << pass 
                              << ", begin batch " << bi << endl;
                         for ( int j = 0; j < batches[bi].isize( ); j++ )
                         {    int xi = batches[bi][j];
                              int ec = ridx[xi];
                              PRINT4( j, xi, ec, F[ec].size( ) );    }    }    }

               for ( int j = 0; j < batches[bi].isize( ); j++ )
               {    int xi = batches[bi][j];
                    int ec = ridx[xi];
     
                    // Logging.
     
                    ostringstream out;
                    ostream& xout = ( !logc.DIRECT ? out : cout );
                    if ( logc.verb[ "ULTRA" ] >= 1 )
                         xout << "\ncorrecting read " << ec << " of " << N << endl;
                    if (logc.ANNOUNCE)
                    {
                         #pragma omp critical
                         {    cout << "begin read " << ec << endl;    }    }
     
                    // Get friends and kmer alignments of the gang to read 0.
     
                    vecbasevector gang;
                    vec< vec< pair<int,int> > > a;
                    vec<int>* p_friends 
                         = ( heur.PREFAB &&  pass == 1 ? &friend_ids[xi] : 0 );
                    GetFriendsAndAlignsInitial<K>( F, reads, ec, gang, p_friends,
                         a, xout, error_model, heur, log_control, logc );
                    if ( logc.verb[ "ULTRA" ] >= 1 )
                         xout << "found " << gang.size( ) - 1 << " friends" << endl;
     
                    // Define the 'true' sequence of read plus flanks, for 
                    // evaluation.
     
                    int gid = -1, gstart = -1;
                    vec<basevector> gkmers;
                    basevector g;
                    if (logc.PRINT_READ) 
                         reads[ec].Print( xout, "read_" + ToString(ec) );
	            const vec<ref_loc>& locs = *log_control.readlocs;
                    if ( locs.nonempty( ) )
                    {    double gclock = WallClockTime( );
                         const int flank = 100;
                         gid = locs[ec].id;
                         gstart = Max( 
                              0, locs[ec].start - flank );
                         int gstop = Min( (*(log_control.G))[gid].isize( ), 
                              locs[ec].stop + flank );
                         g = basevector( 
                         (*(log_control.G))[gid], gstart, gstop-gstart );
                         if ( locs[ec].rc2 ) g.ReverseComplement( );
                         if (logc.PRINT_READ) g.Print( xout, "truth_" + ToString(ec) );
                         for ( int p = 0; p <= g.isize( ) - K; p++ )
                              gkmers.push( g, p, K );
                         Sort(gkmers);
                         if ( logc.PRINT_TIME_USED && omp_get_thread_num( ) == 0 )
                         {    xout << TimeSince(gclock) 
                                   << " used computing genomic kmers"
                                   << endl;    }    }
     
                    // Make threaded blocks, then build "corrected read", and trim K
                    // bases off its ends.
     
                    threaded_blocks tb;
                    MakeBlocks<K>( ec, gang, a, tb, g, gkmers, gid, gstart, xout, 
                         heur.FILTER_BAD_BLOCKS, log_control, logc );
                    double cclock = WallClockTime( );
                    efasta r = tb.MakeCorrectedRead(
                         error_model, out, heur, log_control, logc );
                    double tclock = WallClockTime( );
                    int left = 0, right = 0;
                    for ( int j = 0; j < r.isize( ); j++ )
                    {    if ( r[j] == '{' ) break;
                         left++;    }
                    for ( int j = r.isize( ) - 1; j >= 0; j-- )
                    {    if ( r[j] == '}' ) break;
                         right++;    }
                    if ( left > K ) r = r.substr( K, r.isize( ) - K );
                    if ( right > K ) r.resize( r.isize( ) - K );
                    if ( logc.PRINT_TIME_USED && omp_get_thread_num( ) == 0 )
                    {    xout << TimeSince(tclock) << " used trimming corrected read"
                         << endl;    }
                    if ( r.size( ) > 0 ) 
                    {    
                         #pragma omp critical
                         {    corrected.push_back(r);    
                              cid.push_back(ec);    }    }    

                    // Logging.
     
                    reports[ec] = out.str( );
                    if (logc.ANNOUNCE)
                    {
                         #pragma omp critical
                         {    cout << "end read " << ec << endl;    }    }    }

               if (logc.BATCH_TRACKING)
               {    
                    #pragma omp critical
                    {    cout << "\n" << Date( ) << ": pass = " << pass 
                              << ", end batch " << bi 
                              << ", time used = " << TimeSince(biclock) << endl;
                         for ( int j = 0; j < batches[bi].isize( ); j++ )
                         {    int xi = batches[bi][j];
                              int ec = ridx[xi];
                              PRINT4( j, xi, ec, F[ec].size( ) );    }    }    }

               if ( !logc.ANNOUNCE && !logc.DIRECT )
               {
                    #pragma omp critical
                    {    int dots1 = (100*count)/batches.isize( );
                         count++;
                         int dots2 = (100*count)/batches.isize( );
                         if ( dots2 > dots1 )
                         {    for ( int j = dots1; j < dots2; j++ )
                              {    cout << ".";
                                   if ( (j+1) % 50 == 0 ) cout << "\n";
                                   else if ( (j+1) % 10 == 0 ) cout << " ";    }
                              flush(cout);    }    }    }    }
          cout << Date( ) << ": main loop complete, time used = " 
               << TimeSince(mclock) << endl;

          // Print report.

          for ( int ec = 0; ec < N; ec++ )
               cout << reports[ec];    }    }

template void BuildCorrectedReads<16>( const vecbasevector& reads,
     const IAndOsVec& F, const vec<int>& rid,
     VecEFasta& corrected, vec<int>& cid, const ConsensusScoreModel& error_model,
     const long_heuristics& heur, const long_logging_control& log_control,
     const long_logging& logc, const ref_data& ref, const int NUM_THREADS );

template void BuildCorrectedReads<20>( const vecbasevector& reads,
     const IAndOsVec& F, const vec<int>& rid,
     VecEFasta& corrected, vec<int>& cid, const ConsensusScoreModel& error_model,
     const long_heuristics& heur, const long_logging_control& log_control,
     const long_logging& logc, const ref_data& ref, const int NUM_THREADS );

template void BuildCorrectedReads<24>( const vecbasevector& reads,
     const IAndOsVec& F, const vec<int>& rid,
     VecEFasta& corrected, vec<int>& cid, const ConsensusScoreModel& error_model,
     const long_heuristics& heur, const long_logging_control& log_control,
     const long_logging& logc, const ref_data& ref, const int NUM_THREADS );

template void BuildCorrectedReads<28>( const vecbasevector& reads,
     const IAndOsVec& F, const vec<int>& rid,
     VecEFasta& corrected, vec<int>& cid, const ConsensusScoreModel& error_model,
     const long_heuristics& heur, const long_logging_control& log_control,
     const long_logging& logc, const ref_data& ref, const int NUM_THREADS );

void ReportPeakMem( const String msg )
{    cout << Date( ) << ": " << msg << ( msg != "" ? ", " : "" ) 
          << "peak mem = " << setiosflags(ios::fixed)
          << setprecision(1) << PeakMemUsageBytes( ) / 1000000000.0
          << resetiosflags(ios::fixed) << " GB" << endl;    }

int SelectK2( const VecEFasta& corrected, const double K2frac,
     const long_logging& logc, const long_heuristics& heur )
{    int K2 = -1;
     double hclock = WallClockTime( );
     vec<int> lens;
     for ( size_t i = 0; i < corrected.size( ); i++ )
     {    if ( corrected[i].size( ) == 0 ) continue;
          lens.push_back( corrected[i].Length1( ) );    }
     Sort(lens);
     if ( lens.empty( ) ) DiscovarTools::ExitNoCorrectedReads( );
     int med = Median(lens);
     double target_length = K2frac * double(med);
     if (logc.MIN_LOGGING)
     {    cout << Date( ) << ": " << ToStringAddCommas( lens.size( ) ) 
               << " corrected/closed pairs, having median length " << med 
               << endl;    }
     double min_err = 1000000;
     // iterate over allowable K values
     for ( auto itr=BigK::begin(),end=BigK::end(); itr != end; ++itr )
     {    double err = Abs( target_length - *itr );
          if ( err < min_err ) 
          {    min_err = err;
               K2 = *itr;    }    }
     if ( K2 == -1 )
         FatalErr("Unable to identify a suitable value for K2.");
     if (logc.STATUS_LOGGING) cout << Date( ) << ": using K2 = " << K2 << endl;
     REPORT_TIME( hclock, "used choosing K2" );
     return K2; }

void ReportExcessEdges( const HyperEfasta& he, const vecbasevector& G,
     const Bool PERF_STATS )
{    int64_t gsize = 0;
     for ( int g = 0; g < (int) G.size( ); g++ )
          gsize += G[g].size( );
     int excess_edges = he.EdgeObjectCount( ) - (int) G.size( );
     double excess_edges_per = 1000000.0 * double(excess_edges) / double(gsize);
     cout << "-[" << setiosflags(ios::fixed) << setprecision(0) << excess_edges_per
          << resetiosflags(ios::fixed) << "]- excess edges per Mb in EFASTA" 
          << " (total excess edges in EFASTA = " << excess_edges << ")" << endl;
     if (PERF_STATS)
     {    PerfStat::log( ) << std::fixed << setprecision(0)
               << PerfStat( "excess_edges_efasta", 
               "excess edges per Mb in EFASTA", excess_edges_per );    }    }

void Done( const double clock )
{    cout << "\n" << Date( ) << ": done, time used = " << TimeSince(clock)
          << ", peak mem used = " << setiosflags(ios::fixed) << setprecision(1)
          << PeakMemUsageBytes( ) / 1000000000.0 << resetiosflags(ios::fixed)
          << " GB" << endl;
     Scram();    }

void DumpSimLocs( const vecbasevector& reads, const vec<ref_loc>& readlocs )
{    vec< triple<int,ho_interval,int> > x;
     for ( int id = 0; id < (int) reads.size( ); id++ )
     {    x.push( readlocs[id].id, ho_interval( readlocs[id].start, 
               readlocs[id].stop ), id );    }
     ParallelSort(x);
     cout << "\nsimulated read locs:\n";
     for ( int id = 0; id < (int) reads.size( ); id++ )
     {    if ( x[id].first >= 0 )
          {    cout << x[id].third << "   " << x[id].first << "."
                    << x[id].second << "\n";    }    }    }

int StartStep( const String& IN_SHBV, const int START_STEP )
{    int start_step = 2;
     if ( IN_SHBV != "" )
     {    String start = IN_SHBV;
          if ( start.Contains( ".shbv", -1 ) )
          {    start = start.RevBefore( ".shbv" );
               if ( start.Contains( "." ) )
               {    start = start.RevAfter( "." );
                    if ( start.IsInt( ) )
                         start_step = start.Int( ) + 1;    }    }    }
     if ( START_STEP > 1 ) start_step = START_STEP;
     return start_step;    }

void TestSampleAndReads( const String& SAMPLE, const String& READS,     
     long_logging& logc )
{    if ( SAMPLE == "unknown" && READS == "" )
          FAIL_MSG( "If SAMPLE=unknown, you have to specify READS." );
     if ( SAMPLE == "hpool1" && READS != "#picard" )
          FAIL_MSG( "Currently, SAMPLE=hpool1 only works with READS=#picard." );
     if ( SAMPLE == "hpool2" && READS != "#picard" )
          FAIL_MSG( "Currently, SAMPLE=hpool2 only works with READS=#picard." );
     if ( SAMPLE == "hpool3" && READS != "#picard" )
          FAIL_MSG( "Currently, SAMPLE=hpool3 only works with READS=#picard." );
     if ( SAMPLE == "unknown" && logc.EVAL_CORRECTED )
          FAIL_MSG( "If EVAL_CORRECTED=True, you can't have SAMPLE=unknown." );
     if ( SAMPLE == "unknown" && logc.PRINT_EDITS )
          FAIL_MSG( "If PRINT_EDITS=True, you can't have SAMPLE=unknown." );
     if ( logc.VALIDATE >= 0 && READS != "#picard" )
          FAIL_MSG( "VALIDATE can only be used with READS=#picard." );
     if ( logc.READ_EVAL == "True" && READS != "#picard" )
          FAIL_MSG( "READ_EVAL only works when READS == #picard." );
     if ( ( SAMPLE == "human" || SAMPLE == "hpool2" || SAMPLE == "hpool3" ) 
          && READS == "#picard" )
     {    if ( logc.REFTRACE == "" ) logc.REFTRACE = "False";
          if ( logc.READ_EVAL == "" ) logc.READ_EVAL = "True";    }
     if ( logc.REFTRACE == "" ) logc.REFTRACE = "True";
     if ( logc.READ_EVAL == "" ) logc.READ_EVAL = "False";    }

void PrintPerformanceStats( const double clock, const long_logging& logc )
{    double hours = double( WallClockTime( ) - clock ) / 3600.0;
     if (logc.PERF_STATS)
     {    PerfStat::log( ) << std::fixed << setprecision(2)
               << PerfStat( "etime_h", "elapsed time in hours", hours );
          PerfStat::log( ) << std::fixed << setprecision(1) << PerfStat(
               "mem_gigabytes", "memory usage peak in GB", PeakMemUsageGB( ) );    }
     cout << "\n" << Date( ) << ": done, time used = " << TimeSince(clock)
          << ", peak mem used = " << setiosflags(ios::fixed) << setprecision(1)
          << PeakMemUsageGB( ) << resetiosflags(ios::fixed) << " GB" << endl;    }

void PrintCorrectedReadStats( const VecEFasta& corrected )
{    int n_reads_amb = 0;
     for ( size_t id = 0; id < corrected.size( ); id++ )
          if ( corrected[id].Contains( "{" ) ) n_reads_amb++;
     cout << Date( ) << ": "
          << PERCENT_RATIO( 3, n_reads_amb, (long)corrected.size( ) )
          << " of corrected reads contain ambiguities" << endl;    }

void PrintTail( const parsed_args& command, const String& X, const String& X_actual,
     const double clock, const long_logging& logc )
{    if (logc.STATUS_LOGGING)
     {    String tc;
          SlashFold( command.TheCommand( ), tc );
          cout << "\n" << tc;    }
     if ( X.Contains( "random" ) ) cout << "using X=" << X_actual << endl;
     PrintPerformanceStats( clock, logc );
     cout << "\n================================================================="
          << "===================\n" << endl;    }

void LoadEfastaLong( const String& fn, VecEFasta& corrected,
     vec<int>& cid, vec<pairing_info>& cpartner, const long_logging& logc )
{    double clock = WallClockTime( );
     vec<String> headers;
     LoadEfastaIntoStrings( fn, corrected, headers, True );
     for ( int i = 0; i < headers.isize( ); i++ )
          cid.push_back( headers[i].After( "corrected_" ).Int( ) );
     cpartner.resize( corrected.size(), pairing_info(0,-1,-1) );
     REPORT_TIME( clock, "used loading efasta" );    }
