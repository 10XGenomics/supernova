///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * PQVec.h
 *
 *  Created on: Aug 22, 2014
 *      Author: tsharpe
 */

#ifndef PQVEC_H_
#define PQVEC_H_

#include "Qualvector.h"
#include "feudal/BinaryStream.h"
#include "feudal/MasterVec.h"
#include "feudal/Mempool.h"
#include "system/Assert.h"
#include "system/WorklistN.h"
#include <algorithm>
#include <cstddef>
#include <cstring>
#include <iterator>
#include <memory>
#include <numeric>
#include <vector>

// helper class to do the buffer encoding and decoding
class PQVecEncoder
{
public:
    using byte = unsigned char;

    PQVecEncoder() : mpQV(nullptr) {}
    PQVecEncoder( PQVecEncoder const& )=delete;
    PQVecEncoder( qvec const& qv ) { init(qv); }
    PQVecEncoder& operator=( PQVecEncoder const& )=delete;

    void init( qvec const& qv );

    size_t size() const
    { return std::accumulate(mBlocks.begin(),mBlocks.end(),1ul,
           []( size_t acc, Block const& blk ) { return acc+blk.size(); }); }

    byte* encode( byte* pBuf ) const;

    static void decode( byte const* pqBuf, byte* pQs );

private:
    struct Block
    { Block( byte nQs, byte bits, byte minQ )
      : mNQs(nQs), mBits(bits), mMinQ(minQ) {}
      unsigned size() const { return blockSize(mNQs,mBits); }
      static unsigned blockSize( unsigned nQs, unsigned nBits )
      { return (nQs*nBits+17+7)>>3; }
      byte mNQs; byte mBits; byte mMinQ; };

    std::vector<Block> mBlocks;
    std::vector<unsigned> mCosts;
    qvec const* mpQV;
};

// a compressed qvec.
// the only way you can make one is with a qvec (or by copying).
// the only way you can use one is to turn it into a qvec.
template <class Alloc=MempoolAllocator<unsigned char>>
class PQVecA
{
    using byte = unsigned char;
public:
    using value_type = byte;
    using size_type = unsigned;
    using allocator_type = Alloc;

    // these are the only methods that are interesting
    explicit PQVecA( qvec const& qv )
    { new (&allocator()) Alloc;
      PQVecEncoder enc(qv);
      byte* buf = alloc(enc.size());
      enc.encode(buf); }

    void unpack( qvec* pQV ) const
    { pQV->clear();
      size_type nQs = vSize();
      if ( !nQs ) return;
      pQV->resize(nQs);
      PQVecEncoder::decode(data(),&pQV->front()); }

    operator qvec() const { qvec qv; unpack(&qv); return qv; }

    // all the rest of this crap is boilerplate
    PQVecA() { new (&allocator()) Alloc; }

    explicit PQVecA( Alloc const& alloc ) { new (&allocator()) Alloc(alloc); }

    PQVecA( PQVecA const& that ) { new (&allocator()) Alloc; *this = that; }

    PQVecA( PQVecA&& that )
    { new (&allocator()) Alloc;
      if ( allocator()==that.allocator() )
      { using std::swap; swap(mData,that.mData); }
      else *this = that; }

    ~PQVecA() { clear(); allocator().~Alloc(); }

    PQVecA& operator=( PQVecA const& that )
    { if ( this != &that )
      { clear(); size_type sz = that.size();
        if ( sz ) memcpy(alloc(sz),that.data(),sz); }
      return *this; }

    PQVecA& operator=( PQVecA&& that )
    { if ( allocator()==that.allocator() )
      { using std::swap; swap(mData,that.mData); }
      else *this = that;
      return *this; }

    // number of bytes in compressed representation
    size_type size() const
    { byte const* buf = data(); if ( !buf ) return 0;
      size_t nQs;
      while ( (nQs = *buf++) ) buf += PQVecEncoder::Block::blockSize(nQs,*buf&7)-1;
      return buf-data(); }

    // number of bytes in original, uncompressed representation
    size_type vSize() const
    { byte const* buf = data(); if ( !buf ) return 0;
      size_t nQs, result = 0;
      while ( (nQs = *buf++) )
      { result += nQs; buf += PQVecEncoder::Block::blockSize(nQs,*buf&7)-1; }
      return result; }

    PQVecA& clear()
    { byte* buf = data(); if ( buf ) allocator().deallocate(buf,size());
      mData &= ~PTRMASK;
      return *this; }

    void swap( PQVecA& that )
    { if ( allocator()==that.allocator() )
      { using std::swap; swap(mData,that.mData); }
      else
      { byte* thisBuf = nullptr;
        size_type sz = that.size();
        if ( sz )
        { thisBuf = allocator().allocate(sz); memcpy(thisBuf,that.data(),sz); }
        byte* thatBuf = nullptr;
        sz = size();
        if ( sz )
        { thatBuf = that.allocator().allocate(sz); memcpy(thatBuf,data(),sz); }
        clear().setData(thisBuf);
        that.clear().setData(thatBuf); } }

    size_type allocSize() const { return size(); }
    void readFeudal( BinaryReader& reader, size_t sz, void* )
    { clear(); byte* buf = alloc(sz); reader.read(buf,buf+sz); }
    void writeFeudal( BinaryWriter& writer, void const** ) const
    { byte const* buf = data(); if ( buf ) writer.write(buf,buf+size()); }
    void writeBinary( BinaryWriter& writer ) const
    { size_type sz = size(); writer.write(sz);
      if ( sz ) { byte const* buf=data(); writer.write(buf,buf+sz); } }
    void readBinary( BinaryReader& reader )
    { clear(); size_type sz; reader.read(&sz);
      if ( sz ) { byte* buf = alloc(sz); reader.read(buf,buf+sz); } }

    static size_t externalSizeof() { return 0; }
    static unsigned fixedDataLen() { return 0; }
    static size_type interpretSize( void*, size_t sz ) { return sz; }

    Alloc const& get_allocator() const
    { return reinterpret_cast<Alloc const*>(&mData+1)[-1]; }

    // man, you'd really have to be careful using this
    byte* setData( byte* buf )
    { Assert(!data());
      mData |= reinterpret_cast<size_t>(buf)&PTRMASK; return buf; }

private:
    static_assert(sizeof(byte*)==sizeof(size_t),"Weird pointer size.");
    static_assert(sizeof(Alloc)<=size_t(2),"Allocator too big.");

    static size_t const PTRMASK = 0xffffffffffff; // lowest 48 bits

    byte* data() { return reinterpret_cast<byte*>(mData&PTRMASK); }
    byte const* data() const { return reinterpret_cast<byte*>(mData&PTRMASK); }
    byte* alloc( size_type sz )
    { return setData(allocator().allocate(sz)); }
    Alloc& allocator() { return reinterpret_cast<Alloc*>(&mData+1)[-1]; }

    size_t mData = 0;
};

template <class Alloc>
struct Serializability< PQVecA<Alloc> >
{ typedef SelfSerializable type; };

template <class Alloc>
void swap( PQVecA<Alloc>& v1, PQVecA<Alloc>& v2 ) { v1.swap(v2); }

using PQVec = PQVecA<>;
using VecPQVec = MasterVec<PQVec>;
extern template class OuterVec<PQVec>;

template <class Itr> // Itr is a random-access iterator over const qvec's
void convertCopy( Itr beg, Itr end, VecPQVec::iterator oItr )
{ if ( beg != end )
  { unsigned char* buf = nullptr;
    size_t remain = 0;
    PQVec::allocator_type alloc = oItr->get_allocator();
    size_t maxChunkSz = alloc.getMaxEnchunkableSize();
    size_t maxUncompressed = (5*maxChunkSz+1)/2;
    PQVecEncoder enc;
    while ( beg != end )
    { enc.init(*beg); ++beg;
      size_t need = enc.size();
      if ( need > remain )
      { if ( remain ) alloc.deallocate(buf,remain);
        remain = 0;
        for ( auto itr=beg; itr != end; ++itr )
            if ( (remain += itr->size()) > maxUncompressed )
                break;
        remain = std::accumulate(beg,end,0ul,
                     []( size_t val, qvec const& qv ){ return val+qv.size(); });
        remain = std::max(need,2*remain/5);
        remain = std::min(remain,maxChunkSz);
        buf = alloc.allocate(remain); }
      oItr->clear().setData(buf); ++oItr;
      buf = enc.encode(buf); remain -= need; }
    if ( remain ) alloc.deallocate(buf,remain); } }

template <class Itr> // Itr is a random-access iterator over const qvec's
void convertAppendParallel( Itr beg, Itr end, VecPQVec& vpqv )
{
    size_t const BATCH_SIZE = 100000ul;
    size_t nnn = end-beg;
    vpqv.resize(vpqv.size()+nnn);
    auto oItr = vpqv.end()-nnn;
    size_t nBatches = (nnn+BATCH_SIZE-1)/BATCH_SIZE;
    parallelFor(0ul,nBatches,
            [beg,BATCH_SIZE,nnn,oItr]( size_t batchId ) mutable
            { size_t off1 = batchId*BATCH_SIZE;
              size_t off2 = std::min(nnn,off1+BATCH_SIZE);
              convertCopy(beg+off1,beg+off2,oItr+off1); });
}

template <class Itr> // Itr is a random-access iterator over const qvec's
void convertAssignParallel( Itr beg, Itr end, VecPQVec& vpqv )
{
    vpqv.clear();
    convertAppendParallel(beg,end,vpqv);
}
#endif /* PQVEC_H_ */
