///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * ReadPath.h
 *
 *  Created on: Dec 11, 2013
 *      Author: tsharpe
 */

#ifndef READPATH_H_
#define READPATH_H_

#include "Intvector.h"

// A description of a graph traversal by some sequence (a read, let's say).
// It's just a vector of edge IDs, but it also tells you how many bases at the
// start of the first edge to skip, and how many bases at the end of the last
// edge to skip.
class ReadPath : public IntVec
{
public:
    ReadPath() :  mOffset(0), mLastSkip(0) {}
    ReadPath( int offset )
    : mOffset(offset), mLastSkip(0) {}

    ReadPath( int offset, const vec<int>& edge_list )
    : mOffset(offset), mLastSkip(0) {
	this->assign(edge_list.begin(), edge_list.end());
    }

    int getOffset() const { return mOffset; }
    void setOffset( int offset ) { mOffset = offset; }
    void addOffset( int add ) { mOffset += add; }

    // FirstSkip was replaced by mOffset, which could be negative.  If mOffset is <= 0,
    // then FirstSkip is zero and the read completely covers the start of the edge.
    // If mOffset is positive, then the start of the read is past the start of the edge.
    unsigned getFirstSkip() const { return (mOffset < 0 ? 0u : static_cast<unsigned>(mOffset)); }
    void setFirstSkip( unsigned firstSkip ) { mOffset = firstSkip; }

    // stuff to make this class feudal
    explicit ReadPath( alloc_type const& a )
    : IntVec(a), mOffset(0), mLastSkip(0) {}

    void swap( ReadPath& that )
    { using std::swap;
      swap(mOffset,that.mOffset);
      swap(mLastSkip,that.mLastSkip);
      swap(static_cast<IntVec&>(*this),static_cast<IntVec&>(that)); }

    void readFeudal( BinaryReader& reader, unsigned long dataLen, void* fixed )
    { reader.read(&mOffset); reader.read(&mLastSkip);
      dataLen -= sizeof(mOffset)+sizeof(mLastSkip);
      static_cast<IntVec*>(this)->readFeudal(reader,dataLen,fixed); }

    void writeFeudal( BinaryWriter& writer, void const** pFixed ) const
    { writer.write(mOffset); writer.write(mLastSkip);
      static_cast<IntVec const*>(this)->writeFeudal(writer,pFixed); }

    void writeBinary( BinaryWriter& writer ) const
    {  writer.write(mOffset); writer.write(mLastSkip);
      static_cast<IntVec const*>(this)->writeBinary(writer); }

    void readBinary( BinaryReader& reader )
    {  reader.read(&mOffset); reader.read(&mLastSkip);
      static_cast<IntVec*>(this)->readBinary(reader); }

    static size_t externalSizeof() { return 0ul; }

    friend ostream& operator<<( ostream& os, ReadPath const& rp )
    {
      os << "(" << rp.getOffset() << ") [" << rp.getFirstSkip() << "] ";
      std::copy( rp.begin(), rp.end(),
	      std::ostream_iterator<unsigned int>(os, " ") );
      os << "[NA]";
      return os;
    }


private:
    int mOffset;
    unsigned mLastSkip;
};
SELF_SERIALIZABLE(ReadPath);

typedef MasterVec<ReadPath> ReadPathVec;
extern template class OuterVec<ReadPath>;

#endif /* READPATH_H_ */
