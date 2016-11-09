///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * BuildReadQGraph48.h
 *
 *  Created on: Jan 22, 2014
 *      Author: tsharpe
 */

#ifndef PATHS_LONG_BUILDREADQGRAPH_48_H_
#define PATHS_LONG_BUILDREADQGRAPH_48_H_

#include "String.h"
#include "feudal/ObjectManager.h"
#include "feudal/PQVec.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"

void buildReadQGraph48(String const& work_dir, String const& read_head, 
                        std::string const mspFilename,
                        vecbvec& reads, ObjectManager<VecPQVec>& quals,
                        bool doFillGaps, bool doJoinOverlaps,
                        unsigned minQual, unsigned minFreq,
                        int64_t ignBcBelow,
                        unsigned minBC, vec<int32_t> const* bcp,
                        double minFreq2Fract, unsigned maxGapSize,
                        String const& refFasta,
        		         bool useNewAligner, bool repathUnpathed,
                        HyperBasevector* pHBV, ReadPathVec* pPaths,
				    float const meanMemFrac = 0.9,
                        bool const VERBOSE = False );

#endif /* PATHS_LONG_BUILDREADQGRAPH_48_H_ */
