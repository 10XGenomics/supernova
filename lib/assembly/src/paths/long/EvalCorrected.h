///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Assess corrected reads.  For the moment we assume that all reads are
// forward.  This is a preliminary method that counts 100-mers.
//
// Warning: this writes files corrected.fasta and corrected.aligns.  It would
// be better not to write any files at all.  Also better not to take a lookup
// table as input.

#ifndef EVAL_CORRECTED_H
#define EVAL_CORRECTED_H

#include "Basevector.h"
#include "CoreTools.h"
#include "IntPairVec.h"
#include "efasta/EfastaTools.h"
#include "paths/long/CreateGenome.h"
#include "paths/long/LongProtoTools.h"

void EvalCorrected( 
     const VecEFasta& corrected0,                // corrected reads
     const vec<int>& cid,                        // ids of corrected reads
     const ref_data& ref,
     const long_logging_control& log_control,
     const long_logging& logc );

#endif
