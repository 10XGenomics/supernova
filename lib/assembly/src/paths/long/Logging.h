///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef LONGPROTO_LOGGING_H
#define LONGPROTO_LOGGING_H

#include "CoreTools.h"
#include "ParseSet.h"

class simple_map2 {

     public:

     vec< pair<String,int> > x;

     int operator[]( const String& s ) const
     {    for ( int j = 0; j < x.isize( ); j++ )
               if ( x[j].first == s ) return x[j].second;
          return 0;    }

};

#define LoggingBool( L, DEFAULT, EXPLANATION )                   \
     if ( String(#DEFAULT) == "True" ) L = True;                 \
     else if ( String(#DEFAULT) == "False" ) L = False;          \
     else FatalErr( "Illegal value for " << #L << "." );         \
     for ( int i = 0; i < logging.isize( ); i++ )                \
     {    const String& x = logging[i];                          \
          if ( !x.Contains( "=" ) )                              \
          {    FatalErr( "Illegal logging " << x << "." );    }  \
          if ( x.Before( "=" ) == #L )                           \
          {    used[i] = True;                                   \
               if ( x.After( "=" ) == "True" ) L = True;         \
               else if ( x.After( "=" ) == "False" ) L = False;  \
               else                                              \
               {    FatalErr( "Illegal value for logging for "                 \
                         << #L << ": " << x.After( "=" ) << "." );    }        \
               break;    }    }

#define LoggingInt( L, DEFAULT, EXPLANATION )                       \
     if ( String(#DEFAULT).IsInt( ) ) L = String(#DEFAULT).Int( );  \
     else FatalErr( "Illegal value for " << #L << "." );            \
     for ( int i = 0; i < logging.isize( ); i++ )                   \
     {    const String& x = logging[i];                             \
          if ( !x.Contains( "=" ) )                                 \
          {    FatalErr( "Illegal logging " << x << "." );    }     \
          if ( x.Before( "=" ) == #L )                              \
          {    used[i] = True;                                      \
               if ( x.After( "=" ).IsInt( ) )                       \
                    L = x.After( "=" ).Int( );                      \
               else                                                 \
               {    FatalErr( "Illegal value for logging for "                 \
                         << #L << ": " << x.After( "=" ) << "." );    }        \
               break;    }    }

#define LoggingDouble( L, DEFAULT, EXPLANATION )                          \
     if ( String(#DEFAULT).IsDouble( ) ) L = String(#DEFAULT).Double( );  \
     else FatalErr( "Illegal value for " << #L << "." );                  \
     for ( int i = 0; i < logging.isize( ); i++ )                         \
     {    const String& x = logging[i];                                   \
          if ( !x.Contains( "=" ) )                                       \
          {    FatalErr( "Illegal logging " << x << "." );    }           \
          if ( x.Before( "=" ) == #L )                                    \
          {    used[i] = True;                                            \
               if ( x.After( "=" ).IsDouble( ) )                          \
                    L = x.After( "=" ).Double( );                         \
               else                                                       \
               {    FatalErr( "Illegal value for logging for "                 \
                         << #L << ": " << x.After( "=" ) << "." );    }        \
               break;    }    }

#define LoggingString( L, DEFAULT, EXPLANATION )                 \
     for ( int i = 0; i < logging.isize( ); i++ )                \
     {    const String& x = logging[i];                          \
          if ( !x.Contains( "=" ) )                              \
          {    FatalErr( "Illegal logging " << x << "." );    }  \
          if ( x.Before( "=" ) == #L )                           \
          {    used[i] = True;                                   \
               L = x.After( "=" );                               \
               break;    }    }

class long_logging {

     public:

     long_logging( const String& L = "", const String& VERB = "" )
     {    vec<String> logging;
          ParseStringSet( L, logging );
          vec<Bool> used( logging.size( ), False );
          vec<String> m;
          ParseStringSet( VERB, m );
          for ( int j = 0; j < m.isize( ); j++ )
          {    if ( !m[j].Contains( ":" ) )
               {    cout << "\nI'm sorry, the VERB option " << m[j]
                         << " does not contain a colon, as required." << endl;
                    exit(1);    }
               verb.x.push( m[j].Before( ":" ),
                    m[j].After( ":" ).Int( ) );    }

// Definition of logging:

LoggingBool( PRINT_BLOCKS, False, "print blocks" );

LoggingBool( PRINT_THREADS, False, "print threads" );

LoggingBool( PRINT_GANG_LOCS, False, "print locations of gang" );

LoggingBool( PRINT_EDITS, False, "print edits implied by thread alignments" );

LoggingBool( PRINT_READ, False, "print read and true genome sequence" );

LoggingBool( PRINT_BLOCK_ERROR_RATES, False,
     "print error rate for each block for each friend, and aggregate" );

LoggingBool( TRACE_READS, False, "trace path of reads on assembly" );

LoggingBool( EVAL_CORRECTED, False, "evaluate corrected reads" );

LoggingBool( EVAL_UNCORRECTED, False, "evaluate uncorrected reads" );

LoggingBool( ANALYZE, False, "analyze assembly" );

LoggingBool( UNWIND_TRUTH, False,
     "use truth data to provide some debugging data for UNWIND" );

LoggingBool( DUMP_SIM_LOCS, False, "print out simulated reads locs" );

LoggingBool( TREAT_AS_UNKNOWN, False,
     "even if genome sequence is specified, ignore it" );

LoggingBool( BATCH_TRACKING, False,
     "track timing of error correction batches" );

LoggingBool( PRINT_FAVORING_REF, False,
     "in assessment, print alignments of reads favoring reference" );

LoggingBool( DUMP_PRE, False, 
     "dump pre-corrected reads as frag_reads_pre.{fastb,qualb}" );

LoggingBool( DIRECT, False, "send main output to cout directly, runs as single "
     "thread intended for debugging of crashes, but ANNOUNCE is faster" );

LoggingBool( SHOW_REFTRACE_EVENTS, False,
     "in RefTrace, print catalog of error events" );

LoggingBool( REFTRACE_VARIANTS, False,
     "in RefTrace, find variants by aligning assembly to reference" );

LoggingBool( REFTRACE_VARIANTS_DC, False,
     "call variants using disconected region" );

LoggingBool( REFTRACE_VARIANTS_LIMIT_K1_EFFORT, True,
     "give up on aligning a read to the K=1 bubble graph if long, exponentially scaled run time is detected." );

LoggingBool( REFTRACE_VARIANTS_MINOR_PHASE_RECONSTRUCTION, True,
     "Carefully use phasing information to process variant probabilites." );

LoggingBool( REFTRACE_VARIANTS_SUMMARY, False,
     "in RefTrace, print summary statistics" );

LoggingInt( SHOW_READS_ON_VARIANT, -1,
     "print reads support the variant with given id" );

LoggingBool( DETECT_VARIANT_ON_SINGLE_EDGES, False,
     "Calculate probability of variant on single assembly edges separately" );

LoggingBool( IGNORE_VARIANT_ON_REPETITIVE_EDGE, True,
     "Ignore variants located on short repetitive edges" );

LoggingBool( ANNOUNCE, False,
     "announce start and stop for each read, for debugging of crashes" );

LoggingBool( USE_GENOME_FOR_DUMP, True,
     "use genome to generate .glocs and .gpaths when dumping files" );

LoggingBool( PERF_STATS, False, "generate readable performance stats" );

LoggingBool( PRINT_TIME_USED, False,
     "print time used for each step; useful for performance analysis "
     "and also for debugging where the code crashed" );

LoggingInt( VALIDATE, -1,
     "if set, use original reads to try to validate differences between "
     "reference and this assembly edge; only works when READS=#picard" );

LoggingInt( VALIDATE_TRACE_ID, -1, "trace this read id when running VALIDATE" );

LoggingBool( VALIDATE_DISPLAY_ALL, False,
     "report on status of all edits, not just ones that fail" );

LoggingInt( SHOW_MISSING, 0,
     "print this many kmer clumps that are missing from corrected reads" );

LoggingInt( EVAL_SAMPLE, 50000, "sample size for read evaluation" );

LoggingInt( MIN_ERRS_TO_PRINT, -1,
     "if nonnegative, print alignments of corrected reads to reference "
     "if they have at least this many errors" );

LoggingString( READ_EVAL, "",
     "False (default) or True, to use reads to assess assembly; "
     "default for SAMPLE=human READS=#picard is True" );

LoggingString( REFTRACE, "", "True (default) or False, to evaluate assembly by "
     "tracing reference through it; default for SAMPLE=human READS=#picard is "
     "False, unless HUMAN_CONTROLS is specified" );

LoggingBool( EVAL_ASSEMBLY, False, "Evaluate the assembly accuracy "
      "using an experimental method derived from RefTrace " );

LoggingString( WORDIFY_TRACE, "",
     "comma-separated list of edge ids to be traced by Wordify" );

LoggingString( TRACE_IDS, "", "trace these reads; ParseIntSet format" );

LoggingString( TRACE_PIDS, "", "trace these read pair ids; ParseIntSet format" );

LoggingString( COLS_TO_PRINT, "",
     "show these columns (ParseIntSet) when printing gang" );

LoggingString( TRACE_EDGES, "", "show where reads land on these edges" );

LoggingString( TRACE_EDGES0, "", 
     "show where reads land on these edges on initial assembly" );

LoggingBool( KEEP_LOCS, False,
     "if specified and READS=#picard, keep locations of reads from bam file" );

LoggingBool( PRINT_INITIAL_PAIRS, False, "print info on initial pairs" );

LoggingBool( PULL_APART_DEBUG, False, "debug mode for PullApart" );

LoggingString( LAYOUT, "", "layout option for dot, for example neato" );

LoggingInt( COUNT_COV, 0, "count coverage across datasets" );

LoggingBool( PRINT_BAMS, False, "print paths of bam files that are loaded" );

LoggingBool( STATUS_LOGGING, True, "generate cryptic logging that reports on the "
     "status of intermediate calculations" );

LoggingBool( MIN_LOGGING, True, "generate at least minimal logging" );

LoggingBool( MAKE_FASTG, False, "generates fastg output when OUT_INT_HEAD is ON");

     // Test for illegal arguments.

     Bool fail = False;
     for ( int i = 0; i < used.isize( ); i++ )
     {    if ( !used[i] )
          {    cout << "\nIllegal logging " << logging[i] << "." << endl;
               fail = True;    }    }
     if (fail)
     {    cout << "Abort." << endl;
          _exit(1);    }

     // Special argument checking.

     if ( READ_EVAL != "" && READ_EVAL != "True" && READ_EVAL != "False" )
          FatalErr( "You can only set READ_EVAL to True or False." );
     if ( REFTRACE != "" && REFTRACE != "True" && REFTRACE != "False" )
          FatalErr( "You can only set REFTRACE to True or False." );
     if ( DIRECT && verb[ "ULTRA" ] == 0 )
     {    FatalErr( "Setting DIRECT=True and VERB(ULTRA)=0 doesn't really make "
               "sense.\nPlease set VERB(ULTRA)=1 if that's what you want.");    }
     if (REFTRACE_VARIANTS)
          REFTRACE = "True"; 

}

     simple_map2 verb;
     Bool PRINT_BLOCKS;
     Bool PRINT_THREADS;
     Bool PRINT_GANG_LOCS;
     Bool PRINT_EDITS;
     Bool PRINT_READ;
     Bool PRINT_BLOCK_ERROR_RATES;
     Bool TRACE_READS;
     Bool EVAL_CORRECTED;
     Bool EVAL_UNCORRECTED;
     Bool ANALYZE;
     Bool UNWIND_TRUTH;
     Bool DUMP_SIM_LOCS;
     Bool TREAT_AS_UNKNOWN;
     Bool BATCH_TRACKING;
     Bool PRINT_FAVORING_REF;
     Bool DUMP_PRE;
     Bool DIRECT;
     Bool SHOW_REFTRACE_EVENTS;
     Bool REFTRACE_VARIANTS;
     Bool REFTRACE_VARIANTS_DC;
     Bool REFTRACE_VARIANTS_LIMIT_K1_EFFORT;
     Bool REFTRACE_VARIANTS_SUMMARY;
     Bool REFTRACE_VARIANTS_MINOR_PHASE_RECONSTRUCTION;
     int  SHOW_READS_ON_VARIANT;
     Bool DETECT_VARIANT_ON_SINGLE_EDGES;
     Bool IGNORE_VARIANT_ON_REPETITIVE_EDGE;
     Bool ANNOUNCE;
     Bool VALIDATE_DISPLAY_ALL;
     Bool PERF_STATS;
     Bool PRINT_TIME_USED;
     int SHOW_MISSING;
     int VALIDATE;
     int VALIDATE_TRACE_ID;
     int EVAL_SAMPLE;
     int MIN_ERRS_TO_PRINT;
     String READ_EVAL;
     String REFTRACE;
     Bool EVAL_ASSEMBLY;
     String WORDIFY_TRACE;
     String TRACE_IDS;
     String TRACE_PIDS;
     String COLS_TO_PRINT;
     String TRACE_EDGES;
     String TRACE_EDGES0;
     Bool KEEP_LOCS;
     Bool PRINT_INITIAL_PAIRS;
     Bool PULL_APART_DEBUG;
     String LAYOUT;
     int COUNT_COV;
     Bool PRINT_BAMS;
     Bool STATUS_LOGGING;
     Bool MIN_LOGGING;
     Bool MAKE_FASTG;
     Bool USE_GENOME_FOR_DUMP;

};

#endif
