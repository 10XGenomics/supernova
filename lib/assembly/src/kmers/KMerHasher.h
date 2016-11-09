///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * KMerHasher.h
 *
 *  Created on: Dec 10, 2013
 *      Author: tsharpe
 */

#ifndef KMERS_KMERHASHER_H_
#define KMERS_KMERHASHER_H_

#include <cstdint>

template <unsigned K>
class BuzHasher
{
    constexpr uint64_t rotate( uint64_t val )
    { return (val<<(K&0x3F))|(val>>(-K&0x3F)); }

public:
    typedef unsigned char Call;

    BuzHasher()
    : mBits0{0x4a2b3ce2251bb6e3,0xbcd60b1e3af4a91d,
             0x6381f301d5e0512c,0x957cc4fdca0f4ed2},
      mBitsK{rotate(0x4a2b3ce2251bb6e3),rotate(0xbcd60b1e3af4a91d),
             rotate(0x6381f301d5e0512c),rotate(0x957cc4fdca0f4ed2)}
    {}

    template <class Itr>
    uint64_t hash( Itr itr )
    { return hashFront(itr); }

    template <class Itr>
    uint64_t step( Itr itr, uint64_t prevHash )
    { return stepFrontFwd(prevHash,itr[-1],itr[K-1]); }


    template <class Itr>
    uint64_t hashFront( Itr itr )
    { uint64_t hash = 0;
      unsigned nnn = K;
      while ( nnn-- )
      { hash = (hash << 1) | (hash >> 63);
        hash ^= mBits0[*itr]; ++itr; }
      return hash; }

    template <class Itr>
    uint64_t hashRear( Itr itr )
    { uint64_t hash = 0;
      unsigned nnn = K;
      while ( nnn-- )
      { hash = (hash << 1) | (hash >> 63);
        hash ^= mBits0[*--itr^3]; }
      return hash; }

    uint64_t stepFrontFwd( uint64_t hash, Call pred, Call succ )
    { hash = (hash << 1) | (hash >> 63);
      hash ^= mBitsK[pred];
      hash ^= mBits0[succ];
      return hash; }

    uint64_t stepRearFwd( uint64_t hash, Call pred, Call succ )
    { hash ^= mBitsK[succ^3];
      hash ^= mBits0[pred^3];
      hash = (hash >> 1) | (hash << 63);
      return hash; }

    uint64_t stepFrontRev( uint64_t hash, Call pred, Call succ )
    { hash ^= mBitsK[pred];
      hash ^= mBits0[succ];
      hash = (hash >> 1) | (hash << 63);
      return hash; }

    uint64_t stepRearRev( uint64_t hash, Call pred, Call succ )
    { hash = (hash << 1) | (hash >> 63);
      hash ^= mBits0[pred^3];
      hash ^= mBitsK[succ^3];
      return hash; }

private:
    uint64_t const mBits0[4];
    uint64_t const mBitsK[4];
};

// A hash function for a fixed-length string of base calls.
// Has the miraculous property that the RC hashes to the same value as the
// forward strand.
//
// It does this by computing two BuzHash's:  one on the first half of the
// string, and one on the RC of the 2nd half, and xor-ing them together.
//
// Additionally, being based on BuzHash means that this is a bi-directional
// rolling hash: successive (previous) kmers can be hashed very quickly knowing
// only the hash value of the previous (successor) kmer and just four
// base calls, regardless of the size of K.  I.e., calculating the next or
// previous kmer's hash value is a constant-time operation once an initial
// hash value is calculated.
// Class is templated on K, but methods are templated to take any random-access
// iterator that returns base codes.
template <unsigned K>
class KMerHasher
{
public:
    // compute a stand-alone hash value for some iterator
    // note that this functor version operates on a const object and doesn't
    // prime the pump for rolling operation.
    template <class Itr> // Itr is a random-access iterator over base codes
    size_t operator()( Itr itr ) const
    { return gBH.hashFront(itr) ^ gBH.hashRear(itr+K); }

    // ready ourselves for a rolling hash computation.
    template <class Itr>
    size_t hash( Itr itr )
    { mHashF = gBH.hashFront(itr);
      mHashR = gBH.hashRear(itr+K);
      return mHashF ^ mHashR; }

    // assumes that the last call to this hasher was on --itr.
    template <class Itr>
    size_t stepF( Itr itr )
    { mHashF = gBH.stepFrontFwd(mHashF,itr[-1],itr[EFRONT-1]);
      mHashR = gBH.stepRearFwd(mHashR,itr[SREAR-1],itr[K-1]);
      return mHashF ^ mHashR; }

    // assumes that the last call to this hasher was on ++itr.
    template <class Itr>
    size_t stepR( Itr itr )
    { mHashF = gBH.stepFrontRev(mHashF,itr[0],itr[EFRONT]);
      mHashR = gBH.stepRearRev(mHashR,itr[SREAR],itr[K]);
      return mHashF ^ mHashR; }

private:
    static unsigned const HALFK = (K+1)/2;
    static unsigned const EFRONT = HALFK;
    static unsigned const SREAR = K - HALFK;

    uint64_t mHashF;
    uint64_t mHashR;
    static BuzHasher<HALFK> gBH;
};

template <unsigned K>
BuzHasher<KMerHasher<K>::HALFK> KMerHasher<K>::gBH;

#endif /* KMERS_KMERHASHER_H_ */
