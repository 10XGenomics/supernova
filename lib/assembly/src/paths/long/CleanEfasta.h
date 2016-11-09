///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef CLEAN_EFASTA_H
#define CLEAN_EFASTA_H

#include "CoreTools.h"
#include "paths/HyperEfasta.h"
#include "paths/long/Heuristics.h"
#include "paths/long/Logging.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/Variants.h"

void GetCells( const HyperEfasta& he, vec<vec<int>>& cells );

void MakeEfastaAssembly( HyperEfasta& he, vec<int>& inv, vec<Bool>& hide,
     const vec<VariantSignature>& snp_bubbles, const long_heuristics& heur,
     const long_logging& logc );

void CleanEfasta( HyperEfasta& he, vec<int>& inv, const long_logging& logc );

// FlagEdgesForHiding( he,... ) -- it is intended that he is a HyperBasevector-type
// object.

template<class H> void FlagEdgesForHiding( 
     const H& he, const vec<int>& inv, vec<Bool>& hide,
     const long_logging& logc );

void CollapseBubbles( HyperEfasta& he );

void Reduce( HyperEfasta& he, const int verbosity, const long_logging& logc );

#endif
