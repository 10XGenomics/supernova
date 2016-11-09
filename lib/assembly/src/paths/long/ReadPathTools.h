///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// ReadPath Utilities

#ifndef READPATHTOOLS_H_
#define READPATHTOOLS_H_

#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "10X/paths/ReadPathVecX.h"


// Validates a single read path against an HBV. Returns false if the path is invalid.

bool ValidateReadPath(const HyperBasevector& hbv, const vec<int>& to_left,
		      const vec<int>& to_right, const int offset,
		      const vec<int>& edge_list, String& message,
		      const int read_length = 0);

// Validates all read paths against an HBV. Returns false if any invalid paths are found.

bool ValidateAllReadPaths(const HyperBasevector& hbv, const ReadPathVec& readpaths );
bool ValidateAllReadPaths(const HyperBasevector& hb, const ReadPathVecX& readpaths );
bool ValidateAllReadPaths(const HyperBasevector& hb, const HyperBasevectorX& hbx, const ReadPathVecX& readpaths );

// Displays a path through the HBV given a list of edges and the start base on
// the first edge. If supplied with the corresponding read and quality scores
// these are displayed too.

void DisplayReadPath(ostream &out, const HyperBasevector& hbv,
                     const vec<int>& to_left, const vec<int>& to_right,
                     const int offset, const vec<int>& edge_list,
                     const BaseVec& seq = BaseVec(), const QualVec& quals = QualVec(),
                     bool show_edges = false, bool show_ruler = false,
                     bool show_seq = true, bool show_quals = true,
                     bool show_mismatch = true, bool show_overlap = true);

void DisplayReadPath(ostream &out, const vecbvec& graph_edges, int K,
		     const int offset, const vec<int>& edge_list,
		     const BaseVec& seq, const QualVec& quals,
                     bool show_edges = false, bool show_ruler = false,
                     bool show_seq = true, bool show_quals = true,
                     bool show_mismatch = true, bool show_overlap = true);



#endif /* READPATHTOOLS_H_ */
