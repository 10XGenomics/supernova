///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file FounderAlignment.cc
 * \author tsharpe
 * \date May 24, 2012
 *
 * \brief
 */

#include "paths/long/ultra/FounderAlignment.h"
#include "paths/long/ultra/MultipleAligner.h"
#include "system/Assert.h"
#include <algorithm>
#include <vector>

namespace
{

class Insert
{
public:
    Insert( size_t readIdx, unsigned start, unsigned len )
    : mReadIdx(readIdx), mStart(start), mLen(len) {}

    // compiler-supplied destructor and copying is OK

    size_t readIdx() const { return mReadIdx; }
    unsigned start() const { return mStart; }
    unsigned end() const { return mStart+mLen; }
    unsigned len() const { return mLen; }

    struct LenCompare // comparator to sort on length ascending
    { bool operator()( Insert const& i1, Insert const& i2 )
      { return i1.mLen < i2.mLen; } };

    struct IdxCompare // comparator to sort on read index ascending
    { bool operator()( Insert const& i1, Insert const& i2 )
      { return i1.mReadIdx < i2.mReadIdx; } };

    size_t mReadIdx;
    unsigned mStart;
    unsigned mLen;
};

typedef std::vector<Insert> VecInsert;

struct Mapper : public std::unary_function<unsigned char,char>
{ char operator()( unsigned char val ) const { return "ACGT- "[val]; } };

void print( VecUCharVec const& vucv )
{
    static SpinLockedData gLock;

    SpinLocker locker(gLock);
    typedef VecUCharVec::const_iterator Itr;
    for ( Itr itr(vucv.begin()), end(vucv.end()); itr != end; ++itr )
    {
        std::transform(itr->begin(),itr->end(),
                        std::ostream_iterator<char>(std::cout),Mapper());
        std::cout << std::endl;
    }
}

}
void AlignFriendsToFounder( vecbvec const& friends, size_t founderIdx,
                   VecAlign const& alignments, Scorer const& scorer,
                   VecUCharVec* pMultipleAlignment )
{
    unsigned char const GAP_CODE = 4;
    unsigned char const SPACE_CODE = 5;
    ForceAssertEq(friends.size(),alignments.size());
    ForceAssertLt(founderIdx,friends.size());
    pMultipleAlignment->clear().resize(friends.size());
    typedef vecbvec::const_iterator Itr;
    Itr beg(friends.begin());
    Itr end(friends.end());
    Itr fndr(friends.begin(founderIdx));
    size_t fndrSz = fndr->size();
    VecInsert* inserts = new VecInsert[fndrSz+1];
    VecAlign::const_iterator aItr(alignments.begin());
    VecUCharVec::iterator rItr(pMultipleAlignment->begin());
    for ( Itr itr(beg); itr != end; ++itr,++aItr,++rItr )
    {
        rItr->reserve(3*fndrSz/2);
        if ( itr == fndr )
            std::copy(fndr->begin(),fndr->end(),std::back_inserter(*rItr));
        else
        {
            align const& aln = *aItr;
            ForceAssertGe(aln.StartOnQuery(),0);
            ForceAssertGe(aln.StartOnTarget(),0);
            if ( aln.StartOnQuery() )
                rItr->resize(aln.StartOnQuery(),SPACE_CODE);

            bvec::const_iterator rd(itr->begin(aln.StartOnTarget()));
            for ( int blk = 0; blk != aln.Nblocks(); ++blk )
            {
                int gapSize = aln.Gaps(blk);
                if ( gapSize < 0 )
                    rItr->append(-gapSize,GAP_CODE);
                else if ( gapSize > 0 )
                {
                    size_t rIdx = itr-beg;
                    size_t fIdx = rItr->size();
                    inserts[fIdx].push_back(Insert(rIdx,rd.pos(),gapSize));
                    rd += gapSize;
                }
                int len = aln.Lengths(blk);
                ForceAssertGe(len,0);
                std::copy(rd,rd+len,std::back_inserter(*rItr));
                rd += len;
            }
            ForceAssertLe(rItr->size(),fndrSz);
            rItr->resize(fndrSz,SPACE_CODE);
        }
    }
    if ( fndrSz )
    {
        size_t col = fndrSz;
        size_t endRow = pMultipleAlignment->size();
        UCharVec nullIns;
        UCharVec emptyIns;
        while ( --col )
        {
            VecInsert& vIns = inserts[col];
            if ( !vIns.size() )
                continue;  // NON-STRUCTURED!

            typedef VecInsert::iterator IItr;
            IItr beg(vIns.begin());
            IItr end(vIns.end());
            std::sort(beg,end,Insert::LenCompare());
            MultipleAligner ma(scorer,vIns[0].mLen);
            for ( IItr itr(beg); itr != end; ++itr )
            {
                bvec const& rd = friends[itr->readIdx()];
                ma.addRead(rd.begin(itr->start()),rd.begin(itr->end()));
            }
            size_t consLen = ma.getConsensusLen();
            nullIns.resize(consLen,GAP_CODE);
            emptyIns.resize(consLen,SPACE_CODE);
            std::sort(beg,end,Insert::IdxCompare());
            for ( size_t row = 0; row != endRow; ++row )
            {
                UCharVec& ucv = (*pMultipleAlignment)[row];
                UCharVec::iterator insPtr(ucv.begin(col));
                if ( beg == end || beg->readIdx() != row )
                {
                    if ( insPtr[-1] == SPACE_CODE || insPtr[0] == SPACE_CODE )
                        ucv.insert(insPtr,emptyIns.begin(),emptyIns.end());
                    else
                        ucv.insert(insPtr,nullIns.begin(),nullIns.end());
                }
                else if ( beg->len() == consLen )
                {
                    bvec const& rd = friends[beg->readIdx()];
                    ucv.insert(insPtr,rd.begin(beg->start()),
                                rd.begin(beg->end()));
                    ++beg;
                }
                else
                {
                    bvec const& rd = friends[beg->readIdx()];
                    gvec ins = ma.getAlignment(rd.begin(beg->start()),
                                                rd.begin(beg->end()));
                    ucv.insert(insPtr,ins.begin(),ins.begin()+consLen);
                    ++beg;
                }
            }
            ForceAssert(beg==end);
        }
    }
    delete [] inserts;
}
