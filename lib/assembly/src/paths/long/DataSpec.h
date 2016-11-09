///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef DATA_SPEC_H
#define DATA_SPEC_H

#include "CoreTools.h"
#include "ParseSet.h"

#define DataSpecBool( L, DEFAULT, EXPLANATION )                   \
     if ( String(#DEFAULT) == "True" ) L = True;                  \
     else if ( String(#DEFAULT) == "False" ) L = False;           \
     else FatalErr( "Illegal value for " << #L << "." );          \
     for ( int i = 0; i < spec.isize( ); i++ )                    \
     {    const String& x = spec[i];                              \
          if ( !x.Contains( "=" ) )                               \
          {    FatalErr( "Illegal spec " << x << "." );    }      \
          if ( x.Before( "=" ) == #L )                            \
          {    used[i] = True;                                    \
               if ( x.After( "=" ) == "True" ) L = True;          \
               else if ( x.After( "=" ) == "False" ) L = False;   \
               else                                               \
               {    FatalErr( "Illegal value for spec "           \
                         << x.After( "=" ) << "." );    }         \
               break;    }    }

#define DataSpecInt( L, DEFAULT, EXPLANATION )                      \
     if ( String(#DEFAULT).IsInt( ) ) L = String(#DEFAULT).Int( );  \
     else FatalErr( "Illegal value for " << #L << "." );            \
     for ( int i = 0; i < spec.isize( ); i++ )                      \
     {    const String& x = spec[i];                                \
          if ( !x.Contains( "=" ) )                                 \
          {    FatalErr( "Illegal spec " << x << "." );    }        \
          if ( x.Before( "=" ) == #L )                              \
          {    used[i] = True;                                      \
               if ( x.After( "=" ).IsInt( ) )                       \
                    L = x.After( "=" ).Int( );                      \
               else                                                 \
               {    FatalErr( "Illegal value for spec "             \
                         << x.After( "=" ) << "." );    }           \
               break;    }    }

#define DataSpecDouble( L, DEFAULT, EXPLANATION )                         \
     if ( String(#DEFAULT).IsDouble( ) ) L = String(#DEFAULT).Double( );  \
     else FatalErr( "Illegal value for " << #L << "." );                  \
     for ( int i = 0; i < spec.isize( ); i++ )                            \
     {    const String& x = spec[i];                                      \
          if ( !x.Contains( "=" ) )                                       \
          {    FatalErr( "Illegal spec " << x << "." );    }              \
          if ( x.Before( "=" ) == #L )                                    \
          {    used[i] = True;                                            \
               if ( x.After( "=" ).IsDouble( ) )                          \
                    L = x.After( "=" ).Double( );                         \
               else                                                       \
               {    FatalErr( "Illegal value for spec "                   \
                         << x.After( "=" ) << "." );    }                 \
               break;    }    }

#define DataSpecString( L, DEFAULT, EXPLANATION )                \
     for ( int i = 0; i < spec.isize( ); i++ )                   \
     {    const String& x = spec[i];                             \
          if ( !x.Contains( "=" ) )                              \
          {    FatalErr( "Illegal spec " << x << "." );    }     \
          if ( x.Before( "=" ) == #L )                           \
          {    used[i] = True;                                   \
               L = x.After( "=" );                               \
               break;    }    }

class long_data_spec {

     public:

     long_data_spec( const String& L )
     {    vec<String> spec;
          ParseStringSet( L, spec );
          vec<Bool> used( spec.size( ), False );

// Definition of spec:

DataSpecBool( HPOOL_ORIG, False,
     "if True, use original rather than realigned Fosmid reads" );

DataSpecBool( PACBIO_MORE, False,
     "for NA12878 Fosmid pools, use full set of PacBio data" );

DataSpecBool( SHUFFLE_INPUT, True, "randomize input read order" );

DataSpecBool( REMOVE_FINISHED, False,
     "for NA12878 Fosmid pools, remove reads from 'finished' Fosmids" );

DataSpecInt( HPOOL_COVERAGE_CAP, 150, "default COVERAGE_CAP for hpool2 and hpool3" );

DataSpecInt( SELECT_SEED, 666,
     "random seed for selecting reads from full picard datasets" );

DataSpecDouble( SELECT_FRAC, -1, 
     "fraction of full picard dataset to use; default is 1.0 except for "
     "hpool1, which has default 0.125; only one of SELECT_FRAC and "
     "COVERAGE_CAP may be used" );

DataSpecDouble( COVERAGE_CAP, -1, "if set to a positive "
     "value, and nominal coverage exceeds this value, cut coverage to this "
     "value; default for hpool2 and hpool3 is 150; only one of SELECT_FRAC and "
     "COVERAGE_CAP may be used" );

DataSpecString( HUMAN_CONTROLS, "",
     "if specified, a list of integers referring to known NA12878 Fosmids, "
     "or 'all'; this causes reads from those regions to be spiked in, and "
     "evaluation to be done relative to the corresponding haploid Fosmid "
     "references; note that if X is not specified, then in this context "
     "it is interpreted as empty, rather than the entire genome." );

DataSpecBool( NEW_Q2, False, "use reprocessed bams without Q2s if available" );

DataSpecBool( NEW_Q2_NEW_ALIGN, False, 
     "use reprocessed bams without Q2s if available, with new alignments" );

     Bool fail = False;
     for ( int i = 0; i < used.isize( ); i++ )
     {    if ( !used[i] )
          {    cout << "\nIllegal spec " << spec[i] << "." << endl;
               fail = True;    }    }
     if (fail)
     {    cout << "Abort." << endl;
          _exit(1);    }

     if ( SELECT_FRAC >= 0 && COVERAGE_CAP >= 0 )
          FatalErr( "Only one of SELECT_FRAC and COVERAGE_CAP may be specified.");

}

     Bool HPOOL_ORIG;
     Bool PACBIO_MORE;
     Bool REMOVE_FINISHED;
     Bool SHUFFLE_INPUT;
     int HPOOL_COVERAGE_CAP;
     int SELECT_SEED;
     double SELECT_FRAC;
     double COVERAGE_CAP;
     String HUMAN_CONTROLS;
     Bool NEW_Q2;
     Bool NEW_Q2_NEW_ALIGN;

};

#endif
