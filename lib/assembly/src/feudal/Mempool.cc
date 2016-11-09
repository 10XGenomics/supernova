///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * \file Mempool.cc
 * \author tsharpe
 * \date Jul 6, 2009
 *
 * \brief Simple memory sub-allocator.
 */
#include "feudal/Mempool.h"
#include "feudal/TrackingAllocator.h"
#include <iostream>

using std::cout;
using std::endl;

MempoolFinder* MempoolFinder::gpInstance;

void* Mempool::allocate( size_t siz, size_t alignmentReq )
{
    if ( tooBig(siz) )
    {
#ifdef TRACK_MEMUSE
        mpMemUse->alloc(siz);
#endif
        return new char[siz];
    }

    SpinLocker locker(*this);
    void* result;
    if ( !mpChunk || !(result = mpChunk->allocate(siz,alignmentReq)) )
    {
        if ( mpPreallocatedChunk )
        {
            mpChunk = mpPreallocatedChunk;
            mpPreallocatedChunk = 0;
        }
        else
        {
            void* ppp = new char[mChunkSize+sizeof(Chunk)];
            mpChunk = new (ppp) Chunk(mpChunk,mChunkSize);
#ifdef TRACK_MEMUSE
            mpMemUse->alloc(mChunkSize);
#endif
            mTotalSize += mChunkSize;
            mFreeSize += mChunkSize;
        }
        result = mpChunk->allocate(siz,alignmentReq);
    }

    mFreeSize -= siz;
    AssertLe(mFreeSize,mTotalSize);
    return result;
}

void Mempool::free( void* ppp, size_t siz )
{
    if ( tooBig(siz) )
    {
        delete [] static_cast<char*>(ppp);
#ifdef TRACK_MEMUSE
        mpMemUse->free(siz);
#endif
        return;
    }

    Chunk* pPre = 0;
    Chunk* pChunk = 0;

    if ( true )
    {
        SpinLocker locker(*this);
        mFreeSize += siz;
        mpChunk->free(ppp,siz);
        AssertLe(mFreeSize,mTotalSize);
        if ( mFreeSize >= mTotalSize )
        {
            pPre = mpPreallocatedChunk;
            pChunk = mpChunk;
            mpPreallocatedChunk = mpChunk = 0;
#ifdef TRACK_MEMUSE
            mpMemUse->free(mTotalSize);
#endif
            mTotalSize = mFreeSize = 0;
        }
    }

    if ( pPre )
    {
        reportUnusedPreallocation(pPre);
        killChunkChain(pPre);
    }
    else if ( pChunk )
    {
        killChunkChain(pChunk);
    }
}

void Mempool::preAllocate( size_t nBytes )
{
    char* ppp = new char[nBytes+sizeof(Chunk)];
#ifdef TRACK_MEMUSE
    mpMemUse->alloc(nBytes);
#endif

    SpinLocker locker(*this);
    if ( mpPreallocatedChunk )
    {
        delete [] ppp;
#ifdef TRACK_MEMUSE
        mpMemUse->free(nBytes);
#endif
    }
    else
    {
        mpPreallocatedChunk = new (ppp) Chunk(mpChunk,nBytes);
        mTotalSize += nBytes;
        mFreeSize += nBytes;
    }
}

void Mempool::killChunkChain( Chunk* pChunk )
{
    if ( pChunk->mpNext )
        killChunkChain(pChunk->mpNext);
    pChunk->~Chunk();
    delete [] reinterpret_cast<char*>(pChunk);
}

void Mempool::reportUnusedPreallocation( Chunk* pChunk )
{
    cout << "Warning: Mempool has an unused " << pChunk->size() <<
            "-byte preallocation." << endl;
}

void MempoolFinder::leakReport() const
{
    unsigned int trouble = 0;
    Mempool const* end = mMempools + N_POOLS;
    for ( Mempool const* pPool = mMempools; pPool != end; ++pPool )
    {
        size_t nUsed = pPool->bytesInUse();
        if ( nUsed )
            cout << "Warning: Mempool " << pPool-mMempools
                 << " is reclaiming " << nUsed << " leaked bytes." << endl;
        if ( pPool->nRefs() )
            trouble += 1;
    }
    if ( trouble )
        cout << "Warning: There were " << trouble <<
        " leaked mempools which were cleaned up at the end of the program run."
        << endl;
}
