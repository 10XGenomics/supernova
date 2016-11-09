///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "paths/long/ReadOriginTracker.h"
#include "paths/long/RefTrace.h"
#include "paths/long/RefTraceControl.h"
#include "Basevector.h"
#include "Qualvector.h"
#include "PairsManager.h"

ReadOriginTracker::ReadOriginTracker(const RefTraceControl& refc) 
    : ref_control_(refc) 
{
    const PairsManager& pm = refc.GetPM();
    ForceAssertEq((size_t)pm.nReads(), refc.Reads().size());
    read_dirs_.assign(pm.nReads(), UNKNOWN);
    for (size_t i = 0; i < pm.nPairs(); i++) {
        read_dirs_[pm.ID1(i)] = PLUS;
        read_dirs_[pm.ID2(i)] = MINUS;
    }
}

vec<String> ReadOriginTracker::getSampleList() const 
{
    return ref_control_.getSampleList();
}

int ReadOriginTracker::getSampleID(size_t read_ID) const
{
    return ref_control_.getSampleID(read_ID);
}

String ReadOriginTracker::getSampleName(size_t read_ID) const
{
    return ref_control_.getSampleName(read_ID);
}
       
const vecbasevector& ReadOriginTracker::Reads() const
{
    return ref_control_.Reads();
}

const vecqualvector& ReadOriginTracker::Quals() const
{
    return ref_control_.Quals();
}
