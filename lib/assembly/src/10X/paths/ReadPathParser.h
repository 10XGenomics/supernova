// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

#ifndef READPATHPARSER_H_
#define READPATHPARSER_H_

#include "10X/paths/ReadPathDefs.h"

class RPParser{

    public:

    /* P S E U D O   C O N S T R U C T O R S */

    // zips a readpath into a readpathX
    void LLzip(uCharArray& bitseq, const ReadPath& rp, const HyperBasevectorX& hb);
 
    void LLzip(uCharArray& bitseq, const vec<int>& rp, const int offset, const HyperBasevectorX& hb);

    /* A C C E S S O R S: U N C O M P R E S S O R S */

    void LLunzip(const uCharArray& bitseq, ReadPath& rp, const HyperBasevectorX& hb, const int64_t index = 0)const;
    
    void LLunzip(const uCharArray& bitseq, vec<int>& rp, int& offset, const HyperBasevectorX& hb, const int64_t index = 0) const;


    /* E N C R Y P T I O N   A L G O R I T H M S */

    // standard 2bit scheme OR huffman scheme
    void LLencodeBranchId(uCharArray& bitseq, const unsigned int& id, unsigned int& sub_idx, int64_t& idx, CODING_SCHEME sch=FIXED_WIDTH);

    int LLdecodeBranchId(const uCharArray& bitseq, unsigned int& sub_idx, int64_t& idx, CODING_SCHEME sch=FIXED_WIDTH) const;


    /* A C C E S S O R S  &  M O D I F I E R S */

    int LLgetOffset(const uCharArray& bitseq, const int64_t idx = 0) const;

    int LLgetFirstEdge(const uCharArray& bitseq, const int64_t idx = 0) const;

    int LLgetNumEdges(const uCharArray& bitseq, const int64_t idx = 0) const; 

    int LLgetNumBytes(const uCharArray& bitseq) const;

    /* /1* S T D   O/P *1/ */

    // well formatted output, use only with RPX
    void LLprintMe(const uCharArray& bitseq, const HyperBasevectorX& hb, const std::string messg="", const int64_t index=0) const;

};

#endif /* READPATHPARSER_H_ */
