///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef LOAD_CORRECT_CORE_H
#define LOAD_CORRECT_CORE_H

#include "Basevector.h"
#include "CoreTools.h"
#include "Qualvector.h"
#include "efasta/EfastaTools.h"
#include "paths/long/DataSpec.h"
#include "paths/long/Logging.h"
#include "paths/long/LongProtoTools.h"

void CapQualityScores( vecqualvector& cquals, const vec<Bool>& done );

// zero all quality scores associated with corrections (for reads in file)
void ZeroCorrectedQuals( String const& readsFile, vecbvec const& creads,
                            vecqvec* pQuals );
// zero all quality scores associated with corrections (for reads already in memory)
void ZeroCorrectedQuals( vecbasevector const& readsFile, vecbvec const& creads,
                            vecqvec* pQuals );

void SamIAm( const int i, const String& getsam, const String& TMP,
                bool keepLocs = false, String const& dexterLibs = "",
                const Bool PF_ONLY = False, const Bool KEEP_NAMES = False );

void MergeReadSets( const vec<String>& heads, const String& TMP, 
     const long_logging& logc, const Bool KEEP_NAMES = False );

void SelectRandom( const String& TMP, const double SELECT_FRAC, 
     const long_logging& logc, const long_data_spec& spec );

void CorrectionSuite( LongProtoTmpDirManager& tmp_mgr, const long_heuristics& heur,
     const long_logging& logc, const long_logging_control& log_control,
     vecbasevector& creads, VecEFasta& corrected, vec<int>& cid, 
     vec<pairing_info>& cpartner, const uint NUM_THREADS, const String& EXIT, 
     const double clock, bool useOldLRPMethod );

void CorrectionSuite( const String& TMP, const long_heuristics& heur,
     const long_logging& logc, const long_logging_control& log_control,
     vecbasevector& creads, VecEFasta& corrected, vec<int>& cid,
     vec<pairing_info>& cpartner, const uint NUM_THREADS, const String& EXIT,
     const double clock, bool useOldLRPMethod );

void DefinePairingInfo( const LongProtoTmpDirManager& tmp_mgr, const vecbasevector& creads,
     const vec<Bool>& to_delete, vec<int>& cid, VecEFasta& corrected,
     vec<pairing_info>& cpartner, const long_logging& logc );

#endif
