///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
#include "paths/long/RefTraceControl.h"
#include "Basevector.h"
#include "Qualvector.h"

void RefTraceControl::ReadySampleLookup(){ //ready the pair manager
    if (read_head != "" && reads.empty()) {
        reads.ReadAll(read_head + ".fastb");
        quals.ReadAll(read_head + ".qualb");
        pm.Read(read_head + ".pairs");
        pm.makeCache();
        pm.printLibraryStats(std::cout);
        bPairsManagerInitialized = true;
        ForceAssertEq(reads.size(), quals.size());
        ForceAssertEq(reads.size(), pm.nReads());
        cout << Date() << ": Loaded " << reads.size() << " reads from " << 
            read_head << endl;
    }
}

vec<String> RefTraceControl::getSampleList()const{ //get a list of library names
    if( bPairsManagerInitialized) return pm.getLibraryNames();
    else                          return vec<String>();
};

int RefTraceControl::getSampleID(size_t read_ID)const { //get sample index from a read id
    int out= -1;
    if(bPairsManagerInitialized) { 
        int pair_id = pm.getPairID(read_ID);
        if (pair_id != -1)
            out= pm.libraryID(pair_id); 
    }
    return out;
}

String RefTraceControl::getSampleName(size_t read_ID)const { //get sample name from a read id
    String ans = "UNSPECIFIED";
    if(bPairsManagerInitialized){ 
        int pair_id = pm.getPairID(read_ID);
        if (pair_id != -1)
            ans = pm.libraryName(pair_id);
    }
    return ans;
};

