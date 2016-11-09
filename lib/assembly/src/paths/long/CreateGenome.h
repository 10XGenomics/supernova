///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef CREATE_GENOME_H
#define CREATE_GENOME_H

#include "Basevector.h"
#include "CoreTools.h"
#include "IntPairVec.h"
#include "paths/HyperBasevector.h"
#include "paths/long/Logging.h"
#include "paths/long/DiscovarTools.h"

class ref_data {

     public:

     ref_data( ) : LG(12) { }

     vecbasevector G, G3, G3plus;
     vec<int> Gplus_ext;
     vec<HyperBasevector> GH;
     vec<double> ploidy;
     vec<bool> is_circular;
     VecIntPairVec Glocs, G3locs, G3pluslocs;
     int LG;
};

class RefTraceControl;

void ParseRegions( const String& X, vec<String>& regions );

// note that this function is duplicated and modded to
// CreateGenome_Discovar, as a quick fix of memory/speed-related issue.
// Any changes here should also be made there

void CreateGenome( const String& IN_GENOME, const String& SAMPLE, const String& X,
     const String& HUMAN_CONTROLS, String& X_actual, ref_data& ref,
     const double GENOME_SUB_PERCENT,RefTraceControl& RTCtrl );

// note that this function is duplicated and modded from
// CreateGenome, as a quick fix of memory/speed-related issue.
// Any changes here should also be made there

void CreateGenome_Discovar( ref_data& ref, const DiscovarTools::DiscovarRefTraceControl& vControl );

void CreateGlocs( const vecbasevector& G, unsigned const LG, VecIntPairVec& Glocs );

// note that this function is duplicated and modded slightly to
// BuildRefData_discovar, as a quick fix of memory/speed-related issue.
// Any changes here should also be made there

inline void BuildRefData( const String& IN_GENOME, const String& SAMPLE, 
     const String& X, const String& HUMAN_CONTROLS, String& X_actual, ref_data& ref,
     const double GENOME_SUB_PERCENT, const long_logging& logc,RefTraceControl& RTCtrl )
{
     CreateGenome( IN_GENOME, SAMPLE, X, HUMAN_CONTROLS, X_actual, ref,
          GENOME_SUB_PERCENT,RTCtrl );
     ref.G3 = ref.G;
     for ( int g = 0; g < (int) ref.G.size( ); g++ )
     {    if ( ref.is_circular[g] )
          {    ref.G3[g] = Cat( ref.G[g], ref.G[g], ref.G[g] );
               ref.G3plus[g] 
                    = Cat( ref.G3plus[g], ref.G3plus[g], ref.G3plus[g] );    }    }
     CreateGlocs( ref.G, ref.LG, ref.Glocs );
     CreateGlocs( ref.G3, ref.LG, ref.G3locs );
     CreateGlocs( ref.G3plus, ref.LG, ref.G3pluslocs );    }

// note that this function is duplicated and modded slightly from
// BuildRefData, as a quick fix of memory/speed-related issue.
// Any changes to BuildRefData should also be made here

inline void BuildRefData_Discovar( ref_data& ref, const DiscovarTools::DiscovarRefTraceControl& vControl )
{
     CreateGenome_Discovar( ref, vControl );
     ref.G3 = ref.G;
     for ( int g = 0; g < (int) ref.G.size( ); g++ )
     {    if ( ref.is_circular[g] )
          {    ref.G3[g] = Cat( ref.G[g], ref.G[g], ref.G[g] );
               ref.G3plus[g]
                    = Cat( ref.G3plus[g], ref.G3plus[g], ref.G3plus[g] );    }    }
     CreateGlocs( ref.G, ref.LG, ref.Glocs );
     CreateGlocs( ref.G3, ref.LG, ref.G3locs );
     CreateGlocs( ref.G3plus, ref.LG, ref.G3pluslocs );    }

#endif
