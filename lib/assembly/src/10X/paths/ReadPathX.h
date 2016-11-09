// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

#ifndef READPATHX_H_
#define READPATHX_H_

#include "10X/paths/ReadPathParser.h"


class ReadPathX: private RPParser
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

    // default input
    ReadPathX();

    // complete initializer accepting an existing ReadPath and HyperBasevectorX
    ReadPathX(const ReadPath& rp, const HyperBasevectorX& hb);

    // complete initializer accepting an existing intvec and offset and HyperBasevectorX
    ReadPathX(const vec<int>& edge_list, const int offset, const HyperBasevectorX& hb);

    // initialize directly from a correctly formatted vector
    ReadPathX(const uCharArray& V);

    // initialize directly from a subvector
    ReadPathX(const uCharArray& RPY_subunit, int64_t start, int64_t end);

    // initialize taking iterators
    ReadPathX(const_uchararr_iterator start, const_uchararr_iterator end);


    /* P S E U D O   C O N S T R U C T O R S */

    // zips a readpath into a readpathX
    void zip(const ReadPath& rp, const HyperBasevectorX& hb);
 
    void zip(const vec<int>& rp, const int offset, const HyperBasevectorX& hb);

    // counterparts of constructors
    void assign(const uCharArray& V);

    // initialize directly from a subvector
    void assign(const uCharArray& RPY_subunit, int64_t start, int64_t end);

    // initialize taking iterators
    void assign(const_uchararr_iterator start, const_uchararr_iterator end);


    /* C O P Y   &   A S S I G N M E N T */

    ReadPathX(const ReadPathX& source);

    ReadPathX& operator=(const ReadPathX& source);


    /* A C C E S S O R S: U N C O M P R E S S O R S */

    void unzip(ReadPath& rp, const HyperBasevectorX& hb)const;
    
    void unzip( vec<int>& rp, int& offset, const HyperBasevectorX& hb) const;


    /* D E S T R U C T O R S */

    ~ReadPathX();


    /* A C C E S S O R S  &  M O D I F I E R S */

    int getOffset() const;

    int getFirstEdge() const;

    int getNumEdges() const; 

    int getNumBytes() const;

    void swap( ReadPathX& that);


    /* S T D   O/P */

    // raw formatted output
    void printRaw();

    // well formatted output
    void printMe(const HyperBasevectorX& hb, std::string messg="") const;

    // doesn't decode completely
    friend ostream& operator<<(ostream & os, const ReadPathX& rpx);
 
private:
    uCharArray ZippedData;

};
// TODO - Check how this will actually work
SELF_SERIALIZABLE(ReadPathX);

#endif /* READPATHX_H_ */
