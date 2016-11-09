///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * BigKMer.h
 *
 *  Created on: Dec 9, 2013
 *      Author: tsharpe
 */

#ifndef BIGKMER_H_
#define BIGKMER_H_

#include "Basevector.h"
#include "dna/CanonicalForm.h"
#include "feudal/Iterator.h"
#include "kmers/KMerContext.h"
#include <ostream>

template <unsigned K>
class BigKMer
{
public:
    typedef unsigned char Call;
    struct hasher
    {
        typedef BigKMer argument_type;
        size_t operator()( BigKMer const& bigKMer ) const
        { return bigKMer.mHash; }
    };

    class iterator
    :  public std::iterator<std::random_access_iterator_tag,Call,int,void,Call>,
       public IteratorBase<iterator,unsigned,int>
    {
    public:
        iterator( BigKMer const* pBigKMer, unsigned offset )
        : IteratorBase<iterator,unsigned,int>(offset), mpBigKMer(pBigKMer) {}

        Call operator*() const
        { return (*mpBigKMer)[this->pos()]; }

        Call operator[]( int idx ) const
        { return (*mpBigKMer)[this->pos()+idx]; }

    private:
        BigKMer const* mpBigKMer;
    };

    BigKMer() : mpBV(nullptr), mOffset(0), mRC(false), mAssigned(false),
                mHash(0)
    {}

    BigKMer( bvec const& bv, size_t hash,
                KMerContext context=KMerContext(), unsigned offset = 0 )
    : mpBV(&bv), mOffset(offset), mRC(false), mAssigned(false),
      mContext(context), mHash(hash)
    {}

    bvec const& getBV() const { return *mpBV; }
    unsigned getOffset() const { return mOffset; }
    bool isRC() const { return mRC; }

    KMerContext getContext() const { return mContext; }
    void addContext( KMerContext context ) { mContext |= context; }

    bool isUnassigned() const { return !mAssigned; }
    void setAssigned() { mAssigned = true; }
    void setUnassigned() { mAssigned = false; }

    iterator begin() const { return iterator(this,0); }
    iterator end() const { return iterator(this,K); }

    Call operator[]( unsigned offset ) const
    { offset += mOffset;
      return mRC ? (*mpBV)[mpBV->size()-offset-1]^3 : (*mpBV)[offset]; }

    Call front() const { return (*this)[0]; }
    Call back() const { return (*this)[K-1]; }

    void updateLocation( BigKMer const& kmer )
    { mpBV = kmer.mpBV; mOffset = kmer.mOffset; mRC = kmer.mRC;
      AssertEq(mHash,kmer.mHash); }

    BigKMer rc() const
    { unsigned offset = mpBV->size()-mOffset-K;
      return BigKMer(mpBV,offset,!mRC,mAssigned,mContext.rc(),mHash); }

    CanonicalForm getForm() const
    { return CF<K>::getForm(begin()); }

    bool isFwd() const { return getForm()==CanonicalForm::FWD; }
    bool isRev() const { return getForm()==CanonicalForm::REV; }
    bool isPalindrome() const { return getForm()==CanonicalForm::PALINDROME; }

    BigKMer& successor( size_t hash )
    { ++mOffset; mHash = hash; return *this; }

    BigKMer& successor( size_t hash, KMerContext context )
    { successor(hash); mContext = context; return *this; }

    friend bool identical( BigKMer const& bigKmer1, BigKMer const& bigKmer2 )
    { return bigKmer1.mpBV == bigKmer2.mpBV &&
             bigKmer1.mOffset == bigKmer2.mOffset &&
             bigKmer1.mRC == bigKmer2.mRC; }

    friend bool operator==( BigKMer const& bigKmer1, BigKMer const& bigKmer2 )
    { return bigKmer1.mHash==bigKmer2.mHash &&
            (identical(bigKmer1,bigKmer2) ||
             std::equal(bigKmer1.begin(),bigKmer1.end(),bigKmer2.begin())); }

    friend bool operator<( BigKMer const& bigKmer1, BigKMer const& bigKmer2 )
    { auto itr1=bigKmer1.begin(), end=bigKmer1.end();
      for ( auto itr2=bigKmer2.begin(); itr1 != end; ++itr1,++itr2 )
      { if ( *itr1 < *itr2 ) return true;
        if ( *itr1 > *itr2 ) return false; }
      return false; }

    friend ostream& operator<<( ostream& os, BigKMer const& kmer )
    { PrintBasesIter(os,kmer.begin(),kmer.end(),K+1); return os; }

 private:
    BigKMer( bvec const* pBV, unsigned offset, bool rc, bool assigned,
                    KMerContext context, size_t hash )
    : mpBV(pBV), mOffset(offset), mRC(rc), mAssigned(assigned),
      mContext(context), mHash(hash)
    {}

    bvec const* mpBV;
    unsigned mOffset;
    bool mRC;
    bool mAssigned;
    KMerContext mContext;
    size_t mHash;
};


#endif /* BIGKMER_H_ */
