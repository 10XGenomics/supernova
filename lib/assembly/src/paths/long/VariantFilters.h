///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
#ifndef VARIANT_FILTERS_H
#define VARIANT_FILTERS_H
#include "CoreTools.h"

class VariantCallGroup;
class HyperBasevector;

void RemoveRepetitiveEdges(vec<VariantCallGroup>& groups, 
        const vec<size_t>& ref_index, 
        const vec<size_t>& ref_shift, 
        const String& ref_m100_file,
        const HyperBasevector& hbp,
        const vec<pair<int,Bool>>& hbp_to_hb,
        int verbosity = 0);

#endif
