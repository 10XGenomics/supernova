///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file ReadPather.cc
 * \author tsharpe
 * \date Dec 13, 2011
 *
 * \brief
 */
#include "kmers/ReadPather.h"

unsigned const EdgeID::NULLVAL;
size_t const KDef::MAX_OFFSET;

#include "feudal/SmallVecDefs.h"
template class SmallVec< UnipathEvidence, MempoolAllocator<UnipathEvidence> >;

#include "feudal/OuterVecDefs.h"
template class OuterVec<UnipathEvidenceVec>;
