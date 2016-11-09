///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file LongReadsToPaths.h
 * \author tsharpe
 * \date Aug 21, 2012
 *
 * \brief
 */
#ifndef PATHS_LONG_LONGREADSTOPATHS_H_
#define PATHS_LONG_LONGREADSTOPATHS_H_

#include "Basevector.h"
#include "Vec.h"
#include "paths/HyperBasevector.h"
#include "paths/HyperKmerPath.h"
#include "paths/KmerPath.h"
#include "paths/KmerPathInterval.h"

void LongReadsToPaths( vecbvec const& reads,
                            unsigned K, unsigned coverage,
                            unsigned logLevel, bool useOldLRPMethod,
                            HyperBasevector* pHBV,
                            HyperKmerPath* pHKP=nullptr,
                            vecKmerPath* pPaths=nullptr,
                            vecKmerPath* pPathsRC=nullptr,
                            vec<big_tagged_rpint>* pPathsDB=nullptr );

#endif /* PATHS_LONG_LONGREADSTOPATHS_H_ */
