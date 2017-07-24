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
#include "paths/long/LargeKDispatcher.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/MakeKmerStuff.h"
#include "paths/long/SupportedHyperBasevector.h"
#include "paths/long/DiscovarTools.h"
#include "random/Shuffle.h"
#include "reporting/PerfStat.h"
#include "system/ParsedArgs.h"

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
