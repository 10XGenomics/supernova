///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef MERGE_READ_SETS_CORE_H
#define MERGE_READ_SETS_CORE_H

#include "CoreTools.h"

bool CheckFileSetExists(const vec<String>& heads, const String& suffix, bool verbose = false); 

// Merge the feudal files at <heads_in>.suffix into a single file at
// <head_out>.suffix.

void MergeFeudal(const String& head_out, const vec<String>& heads_in, const String& suffix, const vec<uint64_t> & sizes_check); 

// Merge the alignment files at <heads_in>.qltout into a single file at
// <head_out>.qltout.  The dataset sizes are important, because qltout files
// do not know how many reads they contain.

void MergeQltout(const String& head_out, const vec<String>& heads_in, const vec<uint64_t> & sizes);
  
#endif
