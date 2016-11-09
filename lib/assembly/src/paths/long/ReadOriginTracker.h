///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef READ_ORIGIN_TRACKER_H
#define READ_ORIGIN_TRACKER_H

#include "CoreTools.h"
#include "Basevector.h"
#include "Qualvector.h"
class RefTraceControl;

class ReadOriginTracker {
public:
    ReadOriginTracker (const RefTraceControl& refc);
    virtual ~ReadOriginTracker () {};

    vec<String> getSampleList() const;

    int getSampleID(size_t read_ID) const;

    String getSampleName(size_t read_ID) const;

    const vecbasevector& Reads() const;
    const vecqualvector& Quals() const;

    // return the direction of the read +1 or -1. 0 for unknown
    enum READ_DIR {PLUS, MINUS, UNKNOWN};
    READ_DIR ReadDirection(size_t read_id) const { return read_dirs_[read_id]; };
       
private:
    const RefTraceControl& ref_control_;
    vec<READ_DIR> read_dirs_;
};

#endif
