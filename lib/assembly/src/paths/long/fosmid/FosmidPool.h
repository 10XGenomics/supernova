///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// For tracking detailed information about a Fosmid pool via a special file.  For 
// now, hardwired for one pool.  A line in the file consists of blank-separated 
// fields:
// 1 = index of region (not used, for convenience);
// 2 = coordinates (c:x-y, where c is in {1,2,...,22,X});
// 3 = changes, a comma-separated list of the following entity types
//     * junction, looks like seq1|seq2, where seq1 is 'vector' and seq2 is 
//                 'genomic';
//       action: when located on a given edge, seq1 and bases to the left of it are
//               deleted
//     * break, looks like seq1$seq2, where both seq1 and seq2 are genomic;
//       action: break at the junction
//     * edit, looks like seq1-->seq2, a correction to the Fosmid sequence.
//       action: when located on a given edge, seq1 is replaced by seq2;
//     * delete, written DELETED, means that clone entry is to be deleted from set;
//       this change is not processed but serves to alert the reader of the file.
//     * low-coverage, written LOW_COVERAGE, means that clone has low coverage and
//       should be ignored; this change is not processed but serves to alert the 
//       reader of the file.
// Within a line, C++-style comments of the form /*...*/ are recognized and removed.

#ifndef FOSMID_POOL_H
#define FOSMID_POOL_H

#include "CoreTools.h"
#include "FastIfstream.h"

void ParseFosmidPoolMetainfo( vec<String>& regions, 
     vec< vec< pair<String,String> > >& junctions,
     vec< vec< pair<String,String> > >& breaks,
     vec< vec< pair<String,String> > >& edits );

vec<String> GetFinishedFosmidFiles( );

#endif
