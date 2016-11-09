///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef FILL_PAIRS_H
#define FILL_PAIRS_H

#include "Basevector.h"
#include "PairsManager.h"

// FillPairs.  First truncate given reads when the multiplicity of a kmer
// drops below min_freq.  Path these truncated reads, and in cases where the 
// first and last kmer of a pair lie on a single unipath, fill in the sequence
// between them.
//
// This should be rewritten from scratch.

void FillPairs( const vecbvec& bases, const PairsManager& pairs,
                    const int min_freq, vecbvec& filled, bool newMethod,
                    bool useOldLRPMethod );

#endif
