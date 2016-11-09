///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Make the corrected reads into a HyperKmerPath assembly.

#ifndef LONG_HYPER_H
#define LONG_HYPER_H

#include "CoreTools.h"
#include "efasta/EfastaTools.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/PairInfo.h"
#include "paths/long/SupportedHyperBasevector.h"

Bool LongHyper( const String& READS, const VecEFasta& correctede, 
     const vec<pairing_info>& cpartner, SupportedHyperBasevector& shb,
     const long_heuristics& heur,
     const long_logging_control& log_control, const long_logging& logc, 
     const LongProtoTmpDirManager& tmp_mgr, bool useOldLRPMethod );

#endif
