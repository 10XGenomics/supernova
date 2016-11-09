///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * BigKPather.h
 *
 *  Created on: Dec 12, 2013
 *      Author: tsharpe
 */

#ifndef BIGKPATHER_H_
#define BIGKPATHER_H_

#include "Basevector.h"
#include "paths/HyperBasevector.h"
#include "paths/HyperKmerPath.h"
#include "paths/KmerPath.h"
#include "paths/long/ReadPath.h"

// If the vecbvec is actually a bunch of reads, then set lowCoverage to false.
// Set the lowCoverage bool to true if you're building an HBV from edges of some
// other graph, for example, to tell the code that the kmers are going to be
// mostly unique already and that it doesn't make sense to do a Map/Reduce
// cycle.
void buildBigKHBVFromReads( unsigned K, vecbvec const& reads, unsigned coverage,
                                HyperBasevector* pHBV,
                                ReadPathVec* pReadPaths=nullptr,
                                HyperKmerPath* pHKP=nullptr,
                                vecKmerPath* pKmerPaths=nullptr );

#endif /* BIGKPATHER_H_ */
