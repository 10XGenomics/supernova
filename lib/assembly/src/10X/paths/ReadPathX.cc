// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

#include "10X/paths/ReadPathX.h"

#define READPATHX_CC_
#ifdef READPATHX_CC_

// TODO: implement similar stuff for HyperBasevector

/* I T E R A T O R S */

/**
 * @brief returns begin iterator to internal storage
 *
 * @return 
 */
uchararr_iterator ReadPathX::begin() { return ZippedData.begin(); }

/**
 * @brief returns begin iterator to internal storage
 *
 * @return 
 */
uchararr_iterator ReadPathX::end() { return ZippedData.end(); }

/**
 * @brief returns a const begin iterator
 *
 * @return 
 */
const_uchararr_iterator ReadPathX::begin() const { return ZippedData.begin(); }

/**
 * @brief returns a const end iterator
 *
 * @return 
 */
const_uchararr_iterator ReadPathX::end() const { return ZippedData.end(); }

/**
 * @brief returns a cbegin iterator
 *
 * @return 
 */
const_uchararr_iterator ReadPathX::cbegin() const { return ZippedData.cbegin(); }

/**
 * @brief returns a cend iterator
 *
 * @return 
 */
const_uchararr_iterator ReadPathX::cend() const { return ZippedData.cend(); }


/* C O N S T R U C T O R S */

/**
 * @brief default constructor
 */
ReadPathX::ReadPathX() {}

/**
 * @brief  complete initializer accepting an existing ReadPath and HyperBasevectorX
 *
 * @param rp: const ReadPath&
 * @param hb: const HyperBasevectorX&
 */
ReadPathX::ReadPathX(const ReadPath& rp, const HyperBasevectorX& hb){
    zip(rp,hb);
}

/**
* @brief complete initializer accepting an existing intvec and offset and HyperBasevectorX
*
* @param edge_list: const vec<int>&
* @param offset: const int
* @param hb: const HyperBasevectorX&
*/
ReadPathX::ReadPathX(const vec<int>& edge_list, const int offset, const HyperBasevectorX& hb){
    zip(edge_list,offset,hb);
}

/**
 * @brief  initialize directly from a correctly formatted vector
 *
 * @param V: const vec<unsigned char>&
 */
ReadPathX::ReadPathX(const uCharArray& V){
    ZippedData.assign(V.begin(),V.end());
}

/**
 * @brief initialize directly from a subvector
 *
 * @param RPY_subunit: const vec<unsigned char>&
 * @param start: int64_t
 * @param end: int64_t
 */
ReadPathX::ReadPathX(const uCharArray& RPY_subunit, int64_t start, int64_t end){
    ZippedData.assign(RPY_subunit.begin()+start,RPY_subunit.begin()+end);
}

/**
 * @brief initialize taking iterators
 *
 * @param start: const iterator to vec<unsigned char>
 * @param end: const iterator to vec<unsigned char>
 */
ReadPathX::ReadPathX(const_uchararr_iterator start, const_uchararr_iterator end){
    ZippedData.assign(start,end);
}


/* P S E U D O   C O N S T R U C T O R S */

/**
 * @brief zips a readpath into a readpathX
 *
 * @param rp: const ReadPath&
 * @param hb: const HyperBasevectorX&
 */
void ReadPathX::zip(const ReadPath& rp, const HyperBasevectorX& hb){

    LLzip(ZippedData,rp,hb);
}

/**
 * @brief zips an edgelist and offset and appends to existing internal storage
 *
 * @param rp: const vec<int>&
 * @param offset: const int
 * @param hb: const HyperBasevectorX&
 */
void ReadPathX::zip(const vec<int>& rp, const int offset, const HyperBasevectorX& hb){
    LLzip(ZippedData,rp,offset,hb);     
}

/**
 * @brief counterparts of constructors
 *
 * @param V: const vec<unsigned char>&
 */
void ReadPathX::assign(const uCharArray& V){
    ZippedData.assign(V.begin(),V.end());
}

/**
 * @brief initialize directly from a subvector
 *
 * @param RPY_subunit: const vec<unsigned char>&
 * @param start: int64_t
 * @param end: int64_t
 */
void ReadPathX::assign(const uCharArray& RPY_subunit, int64_t start, int64_t end){
    ZippedData.assign(RPY_subunit.begin()+start,RPY_subunit.begin()+end);
}

/**
 * @brief initialize taking iterators
 *
 * @param start: const iterator to vec<unsigned char>
 * @param end: const iterator to vec<unsigned char>
 */
// initialize taking iterators
void ReadPathX::assign(const_uchararr_iterator start, const_uchararr_iterator end){
    ZippedData.assign(start,end);
}


/* C O P Y   &   A S S I G N M E N T */

/**
 * @brief copy constr.
 *
 * @param source: const ReadPathX&
 */
ReadPathX::ReadPathX(const ReadPathX& source){
    ZippedData = source.ZippedData;
}

/**
 * @brief assignment operator
 *
 * @param source: const ReadPathX&
 *
 * @return 
 */
ReadPathX& ReadPathX::operator=(const ReadPathX& source){
    if(this == &source)
        return *this;
    ZippedData = source.ZippedData;
        return *this;
}

/* A C C E S S O R S: U N C O M P R E S S O R S */

/**
 * @brief unzip internal data into a ReadPath
 *
 * @param rp: ReadPath&
 * @param hb: const HyperBasevectorX&
 */
void ReadPathX::unzip(ReadPath& rp, const HyperBasevectorX& hb)const{
    LLunzip(ZippedData,rp,hb);
}

/**
 * @brief unzip internal data into an edge list and an offset
 *
 * @param rp: vec<int>&
 * @param offset: int&
 * @param hb: HyperBasevectorX&
 */
void ReadPathX::unzip( vec<int>& rp, int& offset, const HyperBasevectorX& hb) const{
    LLunzip(ZippedData,rp,offset,hb); 
}


/* D E S T R U C T O R S */

/**
 * @brief destructor
 */
ReadPathX::~ReadPathX(){ 
    ZippedData.clear(); 
}


/* A C C E S S O R S  &  M O D I F I E R S */

/**
 * @brief gets offset
 *
 * @return 
 */
int ReadPathX::getOffset()const {
    return LLgetOffset(ZippedData);
}

/**
 * @brief gets first edge
 *
 * @return 
 */
int ReadPathX::getFirstEdge()const {
    return LLgetFirstEdge(ZippedData);
}

/**
 * @brief gets number of edges in edgelist
 *
 * @return 
 */
int ReadPathX::getNumEdges()const {
    return LLgetNumEdges(ZippedData);
}

/**
 * @brief gets number of bytes for compressed storage
 *
 * @return 
 */
int ReadPathX::getNumBytes()const {
    return LLgetNumBytes(ZippedData);
}

/**
 * @brief does a swap of internal storages of two ReadPathXs
 *
 * @param that
 */
void ReadPathX::swap( ReadPathX& that){
    using std::swap;
    swap(ZippedData,that.ZippedData);
}


/* S T D   O/P */

/**
 * @brief raw formatted output
 */
void ReadPathX::printRaw(){
    cout<<"[";
    for(unsigned int i = 0; i<ZippedData.size(); i++)
        cout<<static_cast<int>(ZippedData[i])<<",";
    cout<<"]"<<endl;
}

/**
 * @brief well formatted output
 *
 * @param hb: const HyperBasevectorX&
 * @param messg: std::string
 */
void ReadPathX::printMe(const HyperBasevectorX& hb, std::string messg) const {
    LLprintMe(ZippedData,hb,messg);
    return;
}

/**
 * @brief partially decoded, well formatted output
 *
 * @param os: ostream&
 * @param rpx: const ReadPathX&
 *
 * @return 
 */
ostream& operator<<(ostream & os, const ReadPathX& rpx){
    os<<"[";
    int size = rpx.getNumEdges();
    os<<"offset: "<<rpx.getOffset()<<" | numEdges: "<<size<<" | Edges: ";
    unsigned int idx = 1+sizeof(int16_t);

    if(size==0){
        os<<"<empty>]";
        return os;
    }

    // set edge1
    int edge = rpx.getFirstEdge();
    os<<edge;
    idx += sizeof(uint32_t);

    if(size==1){
        os<<"]";
        return os;
    }

    os<<", (enc:) ";
    for(unsigned int i = 1+sizeof(int16_t)+sizeof(uint32_t); i < rpx.ZippedData.size(); i++){
        os<<static_cast<int>(rpx.ZippedData[i])<<" ";
    }
    os<<"]";
    return os;
}

#endif /* READPATHX_CC_ */
