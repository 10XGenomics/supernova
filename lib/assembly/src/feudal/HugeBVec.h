///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file HugeBVec.h
 * \author tsharpe
 * \date Sep 16, 2010
 *
 * \brief
 */
#ifndef FEUDAL_HUGEBVEC_H_
#define FEUDAL_HUGEBVEC_H_

#include "dna/Bases.h"
#include "feudal/BinaryStream.h"
#include "feudal/Iterator.h"
#include "system/Assert.h"
#include <algorithm>
#include <cstddef>
#include <cstring>
#include <iterator>
#include <ostream>
#include <vector>

/// A bvec that can be as big as all memory.
class HugeBVec
{
public:
    typedef std::vector<unsigned char> container;
    typedef container::value_type value_type;
    typedef container::size_type size_type;
    typedef container::difference_type difference_type;
    typedef std::iterator<std::random_access_iterator_tag,
                          value_type,
                          difference_type,
                          void,
                          value_type> ItrTagBase;
    class const_iterator
    : public ItrTagBase,
      public IteratorBase<const_iterator,size_type,difference_type>
    {
        typedef IteratorBase<const_iterator,size_type,difference_type> BaseT;
    public:
        const_iterator() : mpHBV(0) {}
        const_iterator( HugeBVec const* pHBV, size_type pos )
        : BaseT(pos), mpHBV(pHBV) {}

        // compiler-supplied copying and destructor are OK

        value_type operator*() const { return (*mpHBV)[pos()]; }
        value_type operator[]( difference_type idx ) const
        { return (*mpHBV)[pos()+idx]; }

    private:
        HugeBVec const* mpHBV;
    };

    class const_rc_iterator
    : public ItrTagBase,
      public IteratorBase<const_rc_iterator,size_type,difference_type>
    {
        typedef IteratorBase<const_rc_iterator,size_type,difference_type> BaseT;
    public:
        const_rc_iterator() : mLast(~0ul), mpHBV(0) {}
        const_rc_iterator( HugeBVec const* pHBV, size_type end, size_type pos )
        : BaseT(pos), mLast(end-1), mpHBV(pHBV) {}

        // compiler-supplied copying and destructor are OK

        value_type operator*() const
        { return (*mpHBV)[mLast-pos()] ^ 3; }

        value_type operator[]( difference_type idx ) const
        { return (*mpHBV)[mLast-(pos()+idx)] ^ 3; }

    private:
        size_type mLast;
        HugeBVec const* mpHBV;
    };

    HugeBVec( char const* fileName )
    { BinaryReader::readFile(fileName,this); }

    HugeBVec( size_type size ) : mSize(size), mVec(physSize()) {}

    HugeBVec() : mSize(0) {}

    // compiler-supplied copying and destructor are OK

    size_type size() const { return mSize; }
    size_type capacity() const { return 4*mVec.capacity(); }

    HugeBVec& clear() { mVec.clear(); mSize = 0; return *this; }

    value_type at( size_type idx ) const
    { ForceAssertLt(idx,mSize); return (mVec[idx>>2] >> 2*(~idx&3))&3; }

    value_type operator[]( size_type idx ) const
    { AssertLt(idx,mSize); return (mVec[idx>>2] >> 2*(~idx&3))&3; }

    void set( size_type idx, value_type val )
    { AssertLt(idx,mSize); AssertLt(val,4u);
      unsigned shift = 2*(~idx&3);
      unsigned char& chr = mVec[idx>>2];
      chr ^= (chr ^ (val<<shift)) & (3<<shift); }

    const_iterator begin( size_type idx = 0 ) const
    { return const_iterator(this,idx); }
    const_iterator end() const
    { return const_iterator(this,size()); }

    const_rc_iterator rcbegin() const
    { return const_rc_iterator(this,size(),0); }
    const_rc_iterator rcend() const
    { return const_rc_iterator(this,size(),size()); }

    // Note that these are coded in a kind of funky way to allow you to go
    // backward starting from a known forward position.
    const_rc_iterator rcbegin( size_type end, size_type idx = 0 ) const
    { return const_rc_iterator(this,end,idx); }
    const_rc_iterator rcend( size_type end ) const
    { return const_rc_iterator(this,end,end); }

    HugeBVec& reserve( size_type nBases )
    { mVec.reserve( (nBases+3)/4 ); return *this; }

    HugeBVec& resize( size_type siz )
    { mSize = siz; mVec.resize(physSize()); return *this; }

    HugeBVec& push_back( value_type val )
    { AssertLt(val,4u);
      size_type idx = mSize++;
      if ( physSize() > mVec.capacity() ) mVec.reserve(2*physSize());
      mVec.resize(physSize());
      mVec[idx>>2] |= ((val&3) << 2*(~idx&3));
      return *this; }

    template <class Itr>
    HugeBVec& append( Itr itr, Itr const& end )
    { using std::distance;
      size_type idx = mSize; mSize += distance(itr,end);
      if ( physSize() > mVec.capacity() ) mVec.reserve(2*physSize());
      mVec.resize(physSize());
      value_type* ppp = &mVec[idx>>2];
      value_type vvv = *ppp;
      if ( (idx & 3) ) vvv >>= 2*(-idx&3);
      while ( idx != mSize )
      { value_type val = *itr; AssertLt(val,4u); ++itr;
        vvv = (vvv << 2) | (val & 3);
        if ( !(++idx & 3) ) *ppp++ = vvv; }
      if ( (idx & 3) ) *ppp = vvv << 2*(-idx&3);
      return *this; }

    void writeBinary( BinaryWriter& bw ) const
    { bw.write(mSize);
      if ( mSize )
      { value_type const* ppp = &mVec[0];
        bw.write(ppp,ppp+physSize()); } }

    void readBinary( BinaryReader& br )
    { br.read(&mSize);
      mVec.resize(physSize());
      if ( mSize )
      { value_type* ppp = &mVec[0];
        br.read(ppp,ppp+physSize()); } }

    void swap( HugeBVec& bv )
    { using std::swap; swap(mSize,bv.mSize); mVec.swap(bv.mVec); }

    friend void swap( HugeBVec& bv1, HugeBVec& bv2 )
    { bv1.swap(bv2); }

    static size_t externalSizeof() { return 0; }

    friend bool operator==( HugeBVec const& v1, HugeBVec const& v2 )
    { return v1.mSize == v2.mSize && v1.mVec == v2.mVec; }

    friend bool operator!=( HugeBVec const& v1, HugeBVec const& v2 )
    { return !(v1 == v2); }

    friend std::ostream& operator<<( std::ostream& os, const HugeBVec& hbv )
    { std::transform(hbv.begin(),
                        hbv.end(),
                        std::ostream_iterator<char>(os),
                        BaseToCharMapper());
      return os; }

    class Builder
    {
    public:
        Builder( char const* fileName )
        : mWriter(fileName), mSize(0), mByte(0)
        { mSeek = mWriter.tell(); mWriter.write(mSize); }

        ~Builder() { if ( mWriter.isOpen() ) close(); }

        void close()
        { if ( mSize & 3 )
          { mByte <<= 2*(-mSize&3); mWriter.write(mByte); }
          mWriter.seek(mSeek); mWriter.write(mSize); mWriter.close(); }

        size_type size() const { return mSize; }

        void push_back( value_type val )
        { AssertLt(val,4u); mByte <<= 2; mByte |= val&3;
          if ( !(++mSize & 3) ) mWriter.write(mByte); }

        template <class Itr>
        void append( Itr itr, Itr const& end )
        { while ( itr != end ) { push_back(*itr); ++itr; } }

    private:
        Builder( Builder const& ); // unimplemented -- no copying
        Builder& operator=( Builder const& ); // unimplemented -- no copying

        BinaryWriter mWriter;
        size_type mSize;
        value_type mByte;
        size_t mSeek;
    };

private:
    size_t physSize() const { return (mSize+3)/4; }

    size_type mSize;
    container mVec;
};

SELF_SERIALIZABLE(HugeBVec);

#endif /* FEUDAL_HUGEBVEC_H_ */
