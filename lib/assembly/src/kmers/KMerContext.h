///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file KMerContext.h
 * \author tsharpe
 * \date Aug 14, 2012
 *
 * \brief
 */
#ifndef KMERS_KMERCONTEXT_H_
#define KMERS_KMERCONTEXT_H_

#include "dna/Bases.h"
#include "system/Assert.h"
#include <ostream>

// keeps track of a kmer's predecessor and successor base codes
class KMerContext
{
public:
    KMerContext() : mVal(0) {}
    KMerContext( unsigned char predCode, unsigned char succCode )
    : mVal(pred2Val(predCode)|succ2Val(succCode)) {}

    friend KMerContext operator|( KMerContext const& c1, KMerContext const& c2 )
    { return KMerContext(c1.mVal|c2.mVal); }

    KMerContext& operator|=( KMerContext const& context )
    { mVal |= context.mVal; return *this; }

    // bit 0 is set if A is a predecessor, bit1=C, bit2=G, bit3=T
    unsigned char getPredecessors() const { return mVal >> 4; }

    // base code of the single predecessor (error if there's not exactly 1)
    unsigned char getSinglePredecessor() const
    { return bits2Val(mVal>>4); }

    bool isPredecessor( unsigned char predCode ) const
    { return pred2Val(predCode) & mVal; }

    void setPredecessor( unsigned char predCode )
    { *this |= finalContext(predCode); }

    void removePredecessor( unsigned char predCode )
    { mVal &= ~pred2Val(predCode); }

    unsigned getPredecessorCount() const
    { return sideCount(getPredecessors()); }

    // bit 0 is set if A is a successor, bit1=C, bit2=G, bit3=T
    unsigned char getSuccessors() const { return mVal & 0xfu; }

    // base code of the single predecessor (error if there's not exactly 1)
    unsigned char getSingleSuccessor() const
    { return bits2Val(mVal&0xfu); }

    bool isSuccessor( unsigned char succCode ) const
    { return succ2Val(succCode) & mVal; }

    void setSuccessor( unsigned char succCode )
    { *this |= initialContext(succCode); }

    void removeSuccessor( unsigned char succCode )
    { mVal &= ~succ2Val(succCode); }

    unsigned getSuccessorCount() const
    { return sideCount(getSuccessors()); }

    KMerContext rc() const
    { return KMerContext(gRCVals[mVal]); }

    bool isSingleContext() const
    { return getSuccessorCount()==1 && getPredecessorCount()==1; }

    friend bool operator==( KMerContext kc1, KMerContext kc2 )
    { return kc1.mVal == kc2.mVal; }

    friend bool operator!=( KMerContext kc1, KMerContext kc2 )
    { return kc1.mVal != kc2.mVal; }

    friend std::ostream& operator<<( std::ostream& os, KMerContext kc )
    { return os << GeneralizedBase::bits2Char(kc.getPredecessors())
		<< "^"
                << GeneralizedBase::bits2Char(kc.getSuccessors()); }

    static KMerContext initialContext( unsigned char succCode )
    { return KMerContext(succ2Val(succCode)); }

    static KMerContext finalContext( unsigned char predCode )
    { return KMerContext(pred2Val(predCode)); }

    static KMerContext NNContext()
    { return KMerContext(0xff); }

private:
    explicit KMerContext( unsigned char val ) : mVal(val) {}

    static unsigned char pred2Val( unsigned char predCode )
    { return Base::val2Bits(predCode) << 4; }

    static unsigned char succ2Val( unsigned char succCode )
    { return Base::val2Bits(succCode); }

    static unsigned sideCount( unsigned char bits )
    { AssertLe(bits,0xfu); return gSideCounts[bits&0xful]; }

    static unsigned char bits2Val( unsigned char bits )
    { unsigned result = gBits2Val[bits];
      AssertLt(result,4u);
      return result; }

    unsigned char mVal;
    static unsigned char gRCVals[];
    static unsigned gSideCounts[];
    static unsigned gBits2Val[];
};

#endif /* KMERS_KMERCONTEXT_H_ */
