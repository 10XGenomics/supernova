///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef LONGPROTO_HEURISTICS_H
#define LONGPROTO_HEURISTICS_H

#include "CoreTools.h"
#include "ParseSet.h"

#define HeuristicBool( HEUR, DEFAULT, EXPLANATION )                     \
     if ( String(#DEFAULT) == "True" ) HEUR = True;                     \
     else if ( String(#DEFAULT) == "False" ) HEUR = False;              \
     else FatalErr( "Illegal value for " << #HEUR << "." );             \
     for ( int i = 0; i < heuristics.isize( ); i++ )                    \
     {    const String& x = heuristics[i];                              \
          if ( !x.Contains( "=" ) )                                     \
          {    FatalErr( "Illegal bool heuristic " << x << "." );    }  \
          if ( x.Before( "=" ) == #HEUR )                               \
          {    used[i] = True;                                          \
               if ( x.After( "=" ) == "True" ) HEUR = True;             \
               else if ( x.After( "=" ) == "False" ) HEUR = False;      \
               else                                                     \
               {    FatalErr( "Illegal value for heuristic "            \
                         << x.After( "=" ) << "." );    }               \
               break;    }    }

#define HeuristicInt( HEUR, DEFAULT, EXPLANATION )                     \
     if ( String(#DEFAULT).IsInt( ) ) HEUR = String(#DEFAULT).Int( );  \
     else FatalErr( "Illegal value for " << #HEUR << "." );            \
     for ( int i = 0; i < heuristics.isize( ); i++ )                   \
     {    const String& x = heuristics[i];                             \
          if ( !x.Contains( "=" ) )                                    \
          {    FatalErr( "Illegal int heuristic " << x << "." );    }  \
          if ( x.Before( "=" ) == #HEUR )                              \
          {    used[i] = True;                                         \
               if ( x.After( "=" ).IsInt( ) )                          \
                    HEUR = x.After( "=" ).Int( );                      \
               else                                                    \
               {    FatalErr( "Illegal value for heuristic "           \
                         << x.After( "=" ) << "." );    }              \
               break;    }    }

#define HeuristicDouble( HEUR, DEFAULT, EXPLANATION )                        \
     if ( String(#DEFAULT).IsDouble( ) ) HEUR = String(#DEFAULT).Double( );  \
     else FatalErr( "Illegal value for " << #HEUR << "." );                  \
     for ( int i = 0; i < heuristics.isize( ); i++ )                         \
     {    const String& x = heuristics[i];                                   \
          if ( !x.Contains( "=" ) )                                          \
          {    FatalErr( "Illegal double heuristic " << x << "." );    }     \
          if ( x.Before( "=" ) == #HEUR )                                    \
          {    used[i] = True;                                               \
               if ( x.After( "=" ).IsDouble( ) )                             \
                    HEUR = x.After( "=" ).Double( );                         \
               else                                                          \
               {    FatalErr( "Illegal value for heuristic "                 \
                         << x.After( "=" ) << "." );    }                    \
               break;    }    }

#define HeuristicString( HEUR, DEFAULT, EXPLANATION )                        \
     HEUR = DEFAULT;                                                         \
     for ( int i = 0; i < heuristics.isize( ); i++ )                         \
     {    const String& x = heuristics[i];                                   \
          if ( !x.Contains( "=" ) )                                          \
          {    FatalErr( "Illegal string heuristic " << x << "." );    }     \
          if ( x.Before( "=" ) == #HEUR )                                    \
          {    used[i] = True;                                               \
               HEUR = x.After( "=" );                                        \
               break;    }    }

class long_heuristics {

     public:

     long_heuristics( const String& HEURISTICS )
     {    vec<String> heuristics;
          ParseStringSet( HEURISTICS, heuristics );
          vec<Bool> used( heuristics.size( ), False );

// Definition of heuristics:
HeuristicBool( USE_BROKEN_FF, True, "use known broken friend-finding" );

HeuristicBool( USE_KMER_ALIGN_METHOD, False, "use KmerAlign(...)" );

HeuristicBool( USE_MULTIPLE_ALIGNER, False,
     "use the multiple aligner for thread consensus" );

HeuristicBool( FILTER_BAD_BLOCKS, True,
     "remove reads having blocks that align with too many errors" );

HeuristicBool( COLUMN_FILTER, True, 
     "identify and remove false friends based on columns in multiple alignment" );

HeuristicBool( COLUMN_FILTER_SUBST_ONLY, True, "identify and remove "
     "false friends based on columns in substitutions-only multiple alignment" );

HeuristicBool( LENGTHEN_HOMOPOLYMERS, True,
     "if thread consensus includes a homopolymer of length n >= 20 bases, "
     "adjoin an alternative of length n+1" );

HeuristicBool( LENGTHEN_HOMOPOLYMERS2, True,
     "try again with LENGTHEN_HOMOPOLYMERS if n >= 40" );

HeuristicBool( LENGTHEN_DINUKES, True,
     "if thread consensus includes a dinucleotide of length n >= 30 bases, "
     "adjoin an alternative of length n+2" );

HeuristicBool( IMPROVE_HOMOPOLYMERS, True,
     "delete unlikely homopolymers length generated by LENGTHEN_HOMOPOLYMERS" );

HeuristicBool( LC_CAREFUL, False,
     "be more careful when deleting low coverage edges" );

HeuristicBool( NEW_SMITH_WAT, True, "use new method for aligning reads to reads" );

HeuristicBool( ERR_FILTER, True, "filter out reads based on too many errors" );

HeuristicBool( PREFAB, True,
     "method based on initial correction of a small fraction of reads" );

HeuristicBool( CHUNK_EDGES, False, "chunk edges" );

HeuristicBool( WORDIFY, True, "wordify" );

HeuristicBool( MAKE_LOOPS, True, "make loops" );

HeuristicBool( GULP, True, "gulp edges" );

HeuristicBool( CP_CONDENSE_HOMOPOLYMERS, True,
     "in CorrectPairs1, attempt to condense homopolymers" );

HeuristicBool( CP_RAISE_ZERO, False,
     "in CorrectPairs1, passed to StrongConsensus1" );

HeuristicBool( CP2, True, "run second pass of CorrectPairs1" );

HeuristicBool( USE_EM, False, "use EM algorithm in Correct1Pre" );

HeuristicBool( USE_PF_ONLY, False, "use only PF reads" );

HeuristicBool( FIX_ASYMMETRY_BUG, True, "fix bug in testing for asymmetric data" );

HeuristicBool( KILL_WEAK_EXITS_SOLO, True,
     "run single-vertex version of KillWeakExits2" );

HeuristicBool( KILL_WEAK_EXITS_SOLO_CAREFUL, True, 
     "applies to KILL_WEAK_EXIT_SOLO" );

HeuristicBool( HQ_DIFF_WINDOW, True, "call HighQualDiffWindow in Correct1Pre" );

HeuristicBool( PULL_APART2, True, "run PullApart2" );

HeuristicBool( PERTURB_TRANSLATIONS, False, "perturb read translations; "
     "use only with Illumina data" );

HeuristicBool( NEW_LC_FILT, False, "use new low-coverage filter" );

HeuristicBool( EXTRA_STEPS, False, "do some extra graph simplification steps" );

HeuristicBool( WEAK_LOOPS, True, "remove weakly supported loops" );

HeuristicBool( DIVINE_STRONG, False, "stronger form of DivineBubbles" );

HeuristicBool( PACBIO_PATCH, False, 
        "use additional pacbio reads to resolve fosmid assembly graph" );

HeuristicBool( CORRECT_PACBIO_PATCH, True, 
        "validate and correct pacbio patches in fosmid assembly" );

HeuristicBool( DETECT_VARIANTS, False, "mark SNP bubble edges in efasta" );

HeuristicBool( PRECORRECT_OLD_NEW, False, "run PreCorrectOldNew" );

HeuristicBool( PRECORRECT_ALT1, False,
     "for Correct1, use only precorrection to get initial read set" );

HeuristicBool( CORRECT_PAIRS, True, "do correction at pair level" );

HeuristicBool( FILL_PAIRS_ALT, False, "use new pair filling code" );

HeuristicInt( READ_SAMPLE_SIZE, -1, "sample size for error "
     "model estimation, if non-positive then then it is computed internally" );

HeuristicInt( CP_MINQ_FLOOR, 10,
     "in CorrectPairs1, floor for minimum quality of closure consensus" );

HeuristicInt( CP_MIN_GLUE, 30,
     "in CorrectPairs1, floor for minimum glue to hold joint stack together");

HeuristicInt( DELTA_MIS, 0, "passed to readstack::GetOffsets1" );

HeuristicInt( MAX_STACK, 10000,
     "maximum stack size to process in error correction" );

HeuristicInt( FF_MAKE_ALIGN_IMPL, 1,
     "algorithm to use in friend finding process" );

HeuristicInt( FF_MIN_FREQ, 2,
     "minimum kmer frequency in friend finding process" );

HeuristicInt( FF_MAX_FREQ, 1000,
     "maximum kmer frequency in friend finding process" );

HeuristicInt( FF_MIN_QUAL, 20,
     "minimum kmer quality score in friend finding process" );

HeuristicInt( FF_COVERAGE, 30,
     "anticipated coverage in friend finding process" );

HeuristicBool( FF_DOWN_SAMPLE, False,
     "downsample in friend finding process" );

HeuristicInt( FF_VERBOSITY, 0,
     "log verbosity in friend finding process" );

HeuristicInt( K, 20,
     "minimum overlap for blocks; currently 16, 20, 24 and 28 are supported" );

HeuristicDouble( CP_MAX_QDIFF, 20.0,
     "in CorrectPairs1, how much lower pair qual can be" );

HeuristicDouble( K2frac, 0.22,
     "K2 is set to this fraction of the median corrected read length, "
     "times an internally determined variable K2frac_mult" );

HeuristicInt( K2_FORCE, -1, "force K2 to be this value" );

HeuristicInt( K2_FLOOR, -1, "minimum value for K2" );

HeuristicBool( CORRECT, True, "if set to False, don't correct reads" );

HeuristicString( PRECORRECT_SEQ, "24,40",
     "sequence of K values to use for Correct1Pre" );

HeuristicBool( REMOVE_VECTOR, True,  
     "attempt to remove vector from edited reads if READS=#picard" );

HeuristicBool( ORIENT_TO_REFERENCE, True,
     "reverse assembly if needed to give it the same orientation as the reference" );

HeuristicBool( REMAP_READS, False,
     "after generating initial assembly, remap the original reads to the assembly, "
     "and where this leads to a unique and different placement of a pair, move it "
     "to the new location" );

HeuristicBool( DIVINE_BRANCHES, False,
     "use a quality score distribution to compare branches" );

HeuristicBool( DIVINE_BRANCHES2, False,
     "use a quality score distribution to compare branches, alt version" );

HeuristicBool( POP_BUBBLES, True, "simplify graph by popping bubbles" );

// Satellite deletion.

HeuristicBool( DELETE_SATELLITE, False, 
     "delete read pairs carrying human satellite DNA" );
HeuristicString( SATELLITE_TARGETS, "all", 
     "targets for DELETE_SATELLITE, parenthesized list" );
HeuristicDouble( MAX_ALPHA_SCORE, 20.0,
     "maximum score to delete read as alpha satellite using DELETE_SATELLITE" );

HeuristicBool( DELETE_Q2, False, "delete Q2 tails of reads" );

HeuristicBool( KILL_HOMOPOLYMER_EDGES, False, 
     "delete graph edges that look too much like a very long homopolymer" );

HeuristicBool( INJECT_REF, True, "inject reference into assembly" );

HeuristicBool( REMOVE_FOSMID_VECTOR, True, "remove Fosmid vector" );

HeuristicBool( DELETE_WEAK, False, "run DeleteWeakEdges" );

HeuristicBool( HOOKUP, False, "run Hookup" );
HeuristicBool( DIVINE_SINGLE_MUTATION, False, "deal with bubble with single mutation" );

HeuristicInt( DIVINE_MAX_LOCS, -1, 
          "maximum number of locations for a kmer in DivineBubbles" );

HeuristicInt( DIVINE_K, 12, "K value for DivineBubbles" );

HeuristicString( REQUIRE_EDGE_MATCH, "",
     "delete edges not having a kmer match to this fastb file" );

     Bool fail = False;
     for ( int i = 0; i < used.isize( ); i++ )
     {    if ( !used[i] )
          {    cout << "\nUnrecognized heuristic " << heuristics[i] << "." << endl;
               fail = True;    }    }
     if (fail)
     {    cout << "Abort." << endl;
          _exit(1);    }

}

     int K;
     double K2frac;
     int K2_FORCE;
     int K2_FLOOR;
     Bool USE_BROKEN_FF;
     Bool USE_KMER_ALIGN_METHOD;
     Bool COLUMN_FILTER;
     Bool COLUMN_FILTER_SUBST_ONLY;
     Bool FILTER_BAD_BLOCKS;
     Bool USE_MULTIPLE_ALIGNER;
     Bool LENGTHEN_HOMOPOLYMERS;
     Bool LENGTHEN_HOMOPOLYMERS2;
     Bool IMPROVE_HOMOPOLYMERS;
     Bool LENGTHEN_DINUKES;
     int READ_SAMPLE_SIZE;
     Bool NEW_SMITH_WAT;
     Bool ERR_FILTER;
     Bool PREFAB;
     Bool CHUNK_EDGES;
     Bool WORDIFY;
     Bool MAKE_LOOPS;
     Bool GULP;
     Bool CP_CONDENSE_HOMOPOLYMERS;
     int CP_MINQ_FLOOR;
     double CP_MAX_QDIFF;
     int CP_MIN_GLUE;
     Bool CP_RAISE_ZERO;
     Bool CP2;
     Bool USE_EM;
     int DELTA_MIS;
     int FF_MAKE_ALIGN_IMPL;
     int FF_MIN_QUAL;
     int FF_MIN_FREQ;
     int FF_MAX_FREQ;
     int FF_COVERAGE;
     Bool FF_DOWN_SAMPLE;
     int FF_VERBOSITY;
     int MAX_STACK;
     Bool FIX_ASYMMETRY_BUG;
     Bool KILL_WEAK_EXITS_SOLO;
     Bool KILL_WEAK_EXITS_SOLO_CAREFUL;
     Bool HQ_DIFF_WINDOW;
     Bool PULL_APART2;
     Bool WORDIFY_ALT;
     Bool PERTURB_TRANSLATIONS;
     Bool NEW_LC_FILT;
     Bool EXTRA_STEPS;
     Bool WEAK_LOOPS;
     Bool DIVINE_STRONG;
     Bool PACBIO_PATCH;
     Bool CORRECT_PACBIO_PATCH;
     Bool DETECT_VARIANTS;
     Bool PRECORRECT_OLD_NEW;
     Bool PRECORRECT_ALT1;
     Bool CORRECT_PAIRS;
     Bool FILL_PAIRS_ALT;
     Bool CORRECT;
     String PRECORRECT_SEQ;
     Bool REMOVE_VECTOR;
     Bool ORIENT_TO_REFERENCE;
     Bool REMAP_READS;
     Bool DIVINE_BRANCHES;
     Bool DIVINE_BRANCHES2;
     Bool POP_BUBBLES;
     Bool DELETE_SATELLITE;
     String SATELLITE_TARGETS;
     double MAX_ALPHA_SCORE;
     Bool DELETE_Q2;
     Bool KILL_HOMOPOLYMER_EDGES;
     Bool INJECT_REF;
     Bool REMOVE_FOSMID_VECTOR;
     Bool DELETE_WEAK;
     Bool HOOKUP;
     Bool DIVINE_SINGLE_MUTATION;
     int DIVINE_MAX_LOCS;
     int DIVINE_K;
     String REQUIRE_EDGE_MATCH;
     Bool LC_CAREFUL;
     Bool USE_PF_ONLY;

};

#endif
