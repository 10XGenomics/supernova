///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * PQVec.cc
 *
 *  Created on: Aug 25, 2014
 *      Author: tsharpe
 */

#include "feudal/PQVec.h"
#include "math/PowerOf2.h"

void PQVecEncoder::init( qvec const& qv )
{
    mpQV = &qv;
    mBlocks.clear();
    mCosts.clear();
    mCosts.reserve(qv.size()+1);
    mCosts.push_back(1); // cost of an empty qv
    auto beg = qv.begin();
    auto end = qv.end();
    auto itr = beg;
    while ( itr != end )
    {
        unsigned char const MAX_Q = 63;
        if ( *itr > MAX_Q )
        {   cout << "\nYour input reads are funny.  I found a quality score of "
                 << unsigned(*itr) << ".\nThe maximum value that I allow is "
                 << unsigned(MAX_Q) << ".\n" << endl;
            Scram(1);    }

        auto iCost = mCosts.end();
        unsigned minVal = std::min(63u,unsigned(*itr));
        unsigned maxVal = *itr;
        unsigned bits = PowerOf2::ceilLg2(maxVal+1u-minVal);
        unsigned prevCost = *--iCost;
        unsigned nQs = 1;
        unsigned bestCost = prevCost + Block::blockSize(nQs,bits);
        Block best(1,bits,minVal);
        auto itr2 = itr;
        ++itr;
        while ( itr2 != beg && nQs < 255 )
        {
            unsigned val = *--itr2;
            if ( val > maxVal )
                maxVal = val;
            if ( val < minVal )
                minVal = val;
            bits = PowerOf2::ceilLg2(maxVal+1u-minVal);
            prevCost = *--iCost;
            unsigned curCost = prevCost + Block::blockSize(++nQs,bits);
            if ( curCost < bestCost )
            {
                bestCost = curCost;
                best = Block(nQs,bits,minVal);
            }
        }
        mCosts.push_back(bestCost);
        unsigned toRemove = best.mNQs - 1;
        if ( !toRemove )
            mBlocks.push_back(best);
        else
        {
            Assert(!mBlocks.empty());
            while ( toRemove > mBlocks.back().mNQs )
            {
                toRemove -= mBlocks.back().mNQs;
                mBlocks.pop_back();
                Assert(!mBlocks.empty());
            }
            if ( toRemove == mBlocks.back().mNQs )
                mBlocks.back() = best;
            else
            {
                mBlocks.back().mNQs -= toRemove;
                mBlocks.push_back(best);
            }
        }
    }
}

PQVecEncoder::byte* PQVecEncoder::encode( byte* pBuf ) const
{
    Assert(mpQV);
    auto itr = mpQV->begin();
    for ( Block const& block : mBlocks )
    {
        uint64_t nQs = block.mNQs;
        uint64_t nBits = block.mBits;
        uint64_t minQ = block.mMinQ;
        *pBuf++ = nQs;
        uint64_t bits = nBits;
        bits |= minQ << 3;
        *pBuf++ = bits;
        bits >>= 8;
        if ( !nBits )
        {
            *pBuf++ = bits;
            itr += nQs;
        }
        else
        {
            uint64_t off = 1;
            while ( nQs-- )
            {
                uint64_t val = *itr - minQ;
                ++itr;
                bits |= val << off;
                if ( (off += nBits) >= 8 )
                {
                    *pBuf++ = bits;
                    off -= 8;
                    bits >>= 8;
                }
            }
            if ( off )
                *pBuf++ = bits;
        }
    }
    *pBuf++ = 0;
    return pBuf;
}

void PQVecEncoder::decode( byte const* pqBuf, byte* pQs )
{
    uint64_t addr = reinterpret_cast<uint64_t>(pqBuf);
    uint64_t* buf = reinterpret_cast<uint64_t*>(addr&~7);
    uint64_t bits = *buf++;
    addr = (addr & 7) << 3;
    uint64_t remain = 64 - addr;
    bits >>= addr;
    uint64_t nQs;
    while ( (nQs = bits&0xff) )
    {
        bits >>= 8;
        if ( !(remain -= 8) )
        {
            bits = *buf++;
            remain = 64;
        }
        uint64_t nBits = bits & 0x07;
        bits >>= 3;
        uint64_t minQ = bits & 0x3f;
        bits >>= 6;
        if ( remain < 9 )
        {
            bits = *buf++;
            minQ |= (bits & 1) << 5;
            bits >>= 1;
            remain += 64;
        }
        remain -= 9;
        if ( !nBits )
            while ( nQs-- )
                *pQs++ = minQ;
        else
        {
            uint64_t mask = (1ul<<nBits)-1ul;
            while ( nQs-- )
            {
                uint64_t val = bits;
                uint64_t used = nBits;
                if ( remain < nBits )
                {
                    bits = *buf++;
                    val |= bits << remain;
                    used -= remain;
                    remain += 64;
                }
                remain -= nBits;
                bits >>= used;
                *pQs++ = minQ + (val & mask);
            }
        }
        bits >>= remain & 7;
        remain &= ~7ul;
        if ( !remain )
        {
            bits = *buf++;
            remain = 64;
        }
    }
}

#include "feudal/OuterVecDefs.h"
template class OuterVec<PQVec>;
