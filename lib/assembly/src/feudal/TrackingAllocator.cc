///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * TrackingAllocator.cc
 *
 *  Created on: May 20, 2013
 *      Author: tsharpe
 */

#include "feudal/TrackingAllocator.h"
#include "String.h"
#include <iostream>

size_t MemUse::gMinReportSize = 256ul*1024*1024;

void MemUse::report()
{
    std::cout << "MemUse by " << mType << ": "
              << ToStringAddCommas(mTotalAllocated)
              << " bytes in " << ToStringAddCommas(mNAllocs)
              << " allocs, inuse=" << ToStringAddCommas(mInUse)
              << ", max use=" << ToStringAddCommas(mMaxUsed)
              << std::endl;
}
