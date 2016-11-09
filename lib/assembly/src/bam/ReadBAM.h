///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * ReadBAM.h
 *
 *  Created on: Aug 29, 2014
 *      Author: tsharpe
 */

#ifndef BAM_READBAM_H_
#define BAM_READBAM_H_

#include "Basevector.h"
#include "VecString.h"
#include "feudal/PQVec.h"

class BAMReader
{
public:
    explicit BAMReader( bool pfOnly=false,
                        bool uniquifyNames=true,
                        double selectFrac=1.,
                        size_t readsToUse=~0ul )
    : mPFOnly(pfOnly), mUniquifyNames(uniquifyNames),
      mSelectFrac(selectFrac), mReadsToUse(readsToUse) {}

    // Appends sequence, quals, and (optionally) names of the paired reads
    // in the given BAM to the given vecvecs.
    void readBAM( String const& bamFile,
                    vecbvec* pVBV, VecPQVec* pVPQV,
                    vecString* pReadNames=nullptr,
				set<String>* pTagsSelect=nullptr,
	   			vecString* pReadTags=nullptr );

private:
    bool mPFOnly;
    bool mUniquifyNames;
    double mSelectFrac;
    size_t mReadsToUse;
};

#endif /* BAM_READBAM_H_ */
