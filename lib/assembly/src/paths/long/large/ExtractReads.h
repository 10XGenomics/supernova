///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef EXTRACT_READS_H
#define EXTRACT_READS_H

#include "Basevector.h"
#include "CoreTools.h"
#include "feudal/PQVec.h"

void ExtractReads( const String& sample, const String& species, String reads,
     String& SELECT_FRAC, const int READS_TO_USE, const vec<String>& regions, 
     const String& tmp_dir1, const String& work_dir, const Bool all, 
     const Bool USE_PF_ONLY, const Bool KEEP_NAMES, vec<String>& subsam_names, 
     vec<int64_t>& subsam_starts, vecbvec* pReads, ObjectManager<VecPQVec>& quals );

#endif
