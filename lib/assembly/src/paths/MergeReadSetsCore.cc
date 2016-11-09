///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "lookup/LookAlign.h"

bool CheckFileSetExists(const vec<String>& heads, const String& suffix, bool verbose ) {
  bool found_all = true;
  for (size_t i = 0; i < heads.size(); ++i)
    if (IsRegularFile(heads[i] + suffix) == false) {
      if (verbose) 
	cout << "Could not find file: " << heads[i] + suffix << endl;
      found_all = false;
    }
  return found_all;
}

// Merge the feudal files at <heads_in>.suffix into a single file at
// <head_out>.suffix.

void MergeFeudal(const String& head_out, const vec<String>& heads_in, const String& suffix, const vec<uint64_t> & sizes_check) {
  vec<String> full_in(heads_in.size());
  for (size_t i = 0; i < heads_in.size(); ++i) {
    full_in[i] = heads_in[i] + suffix;
    ForceAssertEq(MastervecFileObjectCount(full_in[i]), sizes_check[i]);
  }
  Remove(head_out + suffix);
  MergeMastervecs(head_out + suffix, full_in);
}

// Merge the alignment files at <heads_in>.qltout into a single file at
// <head_out>.qltout.  The dataset sizes are important, because qltout files
// do not know how many reads they contain.

void MergeQltout(const String& head_out, const vec<String>& heads_in, const vec<uint64_t> & sizes) {
  
  size_t offset = 0;
  Ofstream(out, head_out + ".qltout");
  
  // Loop over all of the input files.
  for (size_t i = 0 ; i < heads_in.size(); i++) {
    if (IsRegularFile(heads_in[i] + ".mapq"))
         Cp(heads_in[i] + ".mapq", head_out + ".mapq", True);
    
    // Load LookAligns from this input file.
    vec<look_align> aligns;
    LoadLookAligns( heads_in[i] + ".qltout", aligns, /* ignore bads */ True );
    
    // Apply an offset to each of these LookAligns' read IDs, to account for the
    // merging.  Then write to the output file.
    for (size_t j = 0; j < aligns.size(); j++) {
      ForceAssertLt((uint64_t)aligns[j].query_id, sizes[i]);
      aligns[j].query_id += offset;
      aligns[j].PrintParseable(out);
      aligns[j].PrintReadableBrief(out);
    }
    
    offset += sizes[i];

  }
}
