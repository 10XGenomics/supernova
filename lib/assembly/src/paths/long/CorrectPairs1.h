///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef CORRECT_PAIRS1_H
#define CORRECT_PAIRS1_H

#include "Basevector.h"
#include "CoreTools.h"
#include "PairsManager.h"
#include "Qualvector.h"
#include "efasta/EfastaTools.h"
#include "paths/long/LongProtoTools.h"

void CorrectPairs1( String const& tmpDir, const int K, const int max_freq, vecbasevector& bases,
     vecqualvector& quals, const PairsManager& pairs, const vec<Bool>& to_edit, 
     const vec<int>& trace_ids, const long_heuristics& heur, 
     const long_logging_control& log_control, const long_logging& logc,
     VecEFasta& corrected );

#endif
