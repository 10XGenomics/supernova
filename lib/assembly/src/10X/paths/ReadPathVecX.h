// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

#ifndef READPATHVECX_H_
#define READPATHVECX_H_

#include "10X/paths/ReadPathX.h"
#include "10X/paths/ReadPathParser.h"

class ReadPathVecX : private RPParser
{

public:

    /* I T E R A T O R S */

    uchararr_iterator begin();

    uchararr_iterator end();

    const_uchararr_iterator begin() const;

    const_uchararr_iterator end() const;

    const_uchararr_iterator cbegin() const;

    const_uchararr_iterator cend() const;
 

    /* C O N S T R U C T O R S */

    // default 
    ReadPathVecX();

    // reserve appropriate space using heuristics. expect ~ 10byte meanlength
    // we expect at least a few thousand reads per RPY
    ReadPathVecX(const int64_t numReads); 

    ReadPathVecX(const int64_t sid, const int64_t numReads);

    // The following 2 methods should be used sparingly
    // It appends the zipped sequence to existing one, otherwise creates new
    ReadPathVecX(const ReadPath& rp, const HyperBasevectorX& hb);

    // same as above, accepts vec<int> and offset as argument explicitly
    ReadPathVecX(const vec<int>& edge_list, const int offset, const HyperBasevectorX& hb);

    ReadPathVecX(const ReadPathX rpx);

    ReadPathVecX(const ReadPathVec& paths, const HyperBasevectorX& hb, Bool parallel=true);

    ReadPathVecX(const ReadPathVec& paths, const HyperBasevectorX& hb, const int64_t sid, const int64_t numReads, Bool parallel=true);

    ReadPathVecX(const VReadPathVec& paths, const HyperBasevectorX& hb, Bool parallel=true);

    ReadPathVecX(const VReadPathVec& paths, const HyperBasevectorX& hb, const int64_t sid, const int64_t numReads, Bool parallel=true);

    ReadPathVecX(const String Fname);


    /* P S E U D O - C O N S T R U C T O R S */
    
    // It appends the zipped sequence to existing one
    // ideally should reserve mem calling 2nd constructor prior
    void append(const ReadPath& rp, const HyperBasevectorX& hb);

    // same as above, accepts vec<int> and offset as argument explicitly
    void append(const vec<int>& edge_list, const int offset, const HyperBasevectorX& hb);

    // appends raw compressed seq
    void append(const ReadPathX& rpx);

    // appends chunks of ReadPathVecXs
    void append(const ReadPathVecX& rpvx);

    void append(const ReadPathVec& paths, const HyperBasevectorX& hb);

    void append(const ReadPathVec& paths, const HyperBasevectorX& hb, const int64_t sid, const int64_t numReads);

    void append(const VReadPathVec& paths, const HyperBasevectorX& hb);

    void append(const VReadPathVec& paths, const HyperBasevectorX& hb, const int64_t sid, const int64_t numReads);

    void parallel_append(const ReadPathVec& paths, const HyperBasevectorX& hb);

    void parallel_append(const ReadPathVec& paths, const HyperBasevectorX& hb, const int64_t sid, const int64_t numReads);

    void parallel_append(const VReadPathVec& paths, const HyperBasevectorX& hb);

    void parallel_append(const VReadPathVec& paths, const HyperBasevectorX& hb, const int64_t sid, const int64_t numReads);


    // unpacks a ReadPathVecX into a ReadPath
    void unzip(ReadPath& rp, const HyperBasevectorX& hb, const int64_t read_id) const;

    // unpacks a ReadPathVecX into an {edge list, offset}
    void unzip(vec<int>& edge_list, int& offset, const HyperBasevectorX& hb, const int64_t read_id) const;
    
    void unzip(ReadPathVec& paths, const HyperBasevectorX& hb) const;

    void parallel_unzip(ReadPathVec& paths, const HyperBasevectorX& hb) const;

    void unzipAppend(ReadPathVec& paths, const HyperBasevectorX& hb) const;

    void parallel_unzipAppend(ReadPathVec& paths, const HyperBasevectorX& hb) const;
    
    /* C O P Y   &   A S S I G N M E N T */

    ReadPathVecX(const ReadPathVecX& source);

    ReadPathVecX& operator=(const ReadPathVecX& source);


    /* D E S T R U C T O R */

    ~ReadPathVecX();

    void destroy();

    /* P S E D O   M O D I F I E R S   &   A C C E S S O R S  F U N C T I O N S */

    int64_t size() const;

    int64_t storageSize() const;

    int64_t SizeSum() const;

    int64_t numReadPaths () const;

    double meanPathLength() const;

    void extractRPX(ReadPathX& rpx, const int64_t read_id);

    // modifiers and accessors for numEdges
    unsigned char getNumEdges(const int64_t read_id) const; 

    // modifiers and accessors for offsets
    int getOffset(const int64_t read_id) const;

    void setOffset(const int64_t read_id, int offset);

    void addOffset(const int64_t read_id, int add);

    unsigned getFirstSkip(const int64_t read_id) const; 

    void setFirstSkip(const int64_t read_id, unsigned firstSkip);


    /* S T D   O/P */

    // You don't want to use this- outputs too much stuff
    friend ostream& operator<<( ostream& os, const ReadPathVecX& rpvx);

    void printReadPath(const HyperBasevectorX& hb, const std::string messg="", const int64_t index=0) const;

    /* M E M A N A G E R */
    
    // use sparingly - only if you know what you're doing!
    void reserve(const int64_t numReads);

    void reserve(const int64_t sid, const int64_t numReads);

    void resize(const int64_t sz);


    /* R E A D  -  W R I T E */

    // write binary format (ignoring '|' delimiter): 
    // SKIP|bin|STID|bin|EDID|bin|ZDTA|bin|ZIDX|bin
    void readBinary(std::string pathname);
    void ReadAll(std::string pathname);

    void writeBinary(std::string pathname) const;
    void WriteAll(std::string pathname) const;

private:

    int64_t skip = 10; // determines the size of the zipindex
    int64_t start_rid; // inclusive
    int64_t next_start_rid; // exclusive
    uCharArray ZippedData;
    vec<int64_t> ZipIndex;

    /* U P D A T E   &   M A I N T A I N   I N D E X */

    // updates the ZipIndex every skip read_ids, updates the next_start_rid even if not explicitly included in the index
    void updateZipIndex();

    // recompile index after appending a block of rpvx
    void recompileZipIndex();
 
    /* S P E C I A L   I N T E R N A L   I T E R A T O R S */
    // iterator to a random rpx
    const_uchararr_iterator accessItr(const int64_t read_id) const; 

    // use sparingly, not protected
    void jumpPtr(const_uchararr_iterator& itr) const; 

    // index to a random rpx
    int64_t accessIdx(const int64_t read_id) const; 

    // use sparingly, not protected
    void jumpIdx(int64_t& idx) const; 

};

SELF_SERIALIZABLE(ReadPathVecX);

inline void Destroy(ReadPathVecX& v)
{ v.destroy(); }

// extern template class OuterVec<ReadPathVecX>;

#endif
/* READPATHVECX_H_ */
