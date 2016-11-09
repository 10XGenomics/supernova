///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * BuildReadQGraph60.h
 *
 *  Created on: Jan 22, 2014
 *      Author: tsharpe
 */

#ifndef PATHS_LONG_BUILDREADQGRAPH_60_H_
#define PATHS_LONG_BUILDREADQGRAPH_60_H_

#include "String.h"
#include "feudal/ObjectManager.h"
#include "feudal/PQVec.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"

void buildReadQGraph60( vecbvec const& reads, ObjectManager<VecPQVec>& quals,
                        bool doFillGaps, bool doJoinOverlaps,
                        unsigned minQual, unsigned minFreq,
                        double minFreq2Fract, unsigned maxGapSize,
                        String const& refFasta,
       		        bool useNewAligner, bool repathUnpathed,
                        HyperBasevector* pHBV, ReadPathVec* pPaths,
						float const meanMemFrac = 0.9,
                        bool const VERBOSE = False );

#endif /* PATHS_LONG_BUILDREADQGRAPH_60_H_ */
