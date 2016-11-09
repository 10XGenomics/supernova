///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file FounderAlignment.h
 * \author tsharpe
 * \date May 24, 2012
 *
 * \brief
 */
#ifndef PATHS_LONG_FOUNDERALIGNMENT_H_
#define PATHS_LONG_FOUNDERALIGNMENT_H_

#include "Basevector.h"
#include "Charvector.h"
#include "PackAlign.h"
#include "Vec.h"
#include "paths/long/ultra/MultipleAligner.h"
#include <cstddef>

typedef vec<align> VecAlign;

void AlignFriendsToFounder(
        // the reads
        vecbvec const& friends,

        // index of founder read against which alignments were made
        size_t founderIdx,

        // alignments of founder to read for all reads except founderIdx
        VecAlign const& alignments,

        // error model for alignment
        Scorer const& scorer,

        // Multiple alignment is output as a rectangular matrix of small
        // integers.  There are six codes.  The first four are for A=0, C=1,
        // G=2, and T=3.  There is also a code for a pad (=4) marking a
        // deletion, and one (=5) for an initial or final space where the friend
        // doesn't overlap the founder.
        VecUCharVec* pMultipleAlignment
);

#endif /* PATHS_LONG_FOUNDERALIGNMENT_H_ */
