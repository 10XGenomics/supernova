///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "graph/DigraphTemplate.h"
#include "paths/long/KmerCount.h"

template void RemoveHangingEnds<kmer_count>(digraphE<kmer_count>&, 
     int (kmer_count::*)() const, int, double);

template void RemoveHangingEnds3<kmer_count>(digraphE<kmer_count>&, 
     int (kmer_count::*)() const, int, double, int);

template void digraphE<kmer_count>::Used(vec<Bool>&) const;
template void DistancesToEndArr( const digraphE<kmer_count>&, vec<int> const&, int const, Bool const, vec<int>& );
