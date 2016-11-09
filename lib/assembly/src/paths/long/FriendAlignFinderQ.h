///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * FriendAlignFinderQ.h
 *
 *  Created on: Jul 29, 2013
 *      Author: tsharpe
 */

#ifndef PATHS_LONG_FRIENDALIGNFINDERQ_H_
#define PATHS_LONG_FRIENDALIGNFINDERQ_H_

#include "paths/long/FriendAligns.h"
#include "Intvector.h"
#include "dna/Bases.h"
#include "feudal/BinaryStream.h"
#include "feudal/HashSet.h"
#include "feudal/VirtualMasterVec.h"
#include "kmers/KMerContext.h"
#include "kmers/KMer.h"
#include "system/SpinLockedData.h"
#include "system/WorklistN.h"
#include <algorithm>
#include <condition_variable>
#include <iostream>
#include <limits>
#include <vector>
#include <unistd.h>

// FriendContext class to record the occurrence of some kmer at some offset
// within a specified read, whether it occurs there canonically or RC, and the
// context in which it occurs.
class FriendContext
{
public:
    FriendContext()=default;

    FriendContext( size_t readId )
    : mReadId(readId), mOffset(0), mRCOffset(0), mRC(0)
    { AssertLe(readId,std::numeric_limits<unsigned>::max()); }

    FriendContext( size_t readId, unsigned offset, unsigned rcOffset, bool rc,
                                          KMerContext context = KMerContext() )
    : mReadId(readId), mOffset(offset), mRCOffset(rcOffset), mRC(rc),
      mContext(context)
    { ForceAssertLe(readId,std::numeric_limits<unsigned>::max());
      ForceAssertLt(offset,1u<<11);
      ForceAssertLt(rcOffset,1u<<11); }

    size_t readId() const { return mReadId; }
    unsigned offset() const { return mOffset; }
    unsigned rcOffset() const { return mRCOffset; }
    bool isRC() const { return mRC; }

    KMerContext getContext() const { return mContext; }
    void setContext( KMerContext context ) { mContext = context; }

    bool isFriend( FriendContext const& buddy ) const
    { return mContext == KMerContext::NNContext() ||
                buddy.mContext == KMerContext::NNContext() ||
                mContext != buddy.mContext; }

    Friend makeFriend( FriendContext const& buddy ) const
    { bool rc = mRC != buddy.mRC;
      int offset = mOffset - (rc ? buddy.mRCOffset : buddy.mOffset);
      return Friend(buddy.mReadId,offset,rc); }

    friend bool operator==( FriendContext const& f1, FriendContext const& f2 )
    { return f1.mReadId == f2.mReadId &&
            f1.mOffset == f2.mOffset && f1.mRC == f2.mRC; }

    friend bool operator<( FriendContext const& f1, FriendContext const& f2 )
    { if ( f1.mReadId < f2.mReadId ) return true;
      if ( f1.mReadId == f2.mReadId )
      { if ( f1.mOffset < f2.mOffset ) return true;
        if ( f1.mOffset == f2.mOffset && f1.mRC < f2.mRC ) return true; }
      return false; }

    friend ostream& operator<<( ostream& os, FriendContext const& fc )
    { if ( fc.mRC ) os << '~';
      return os << fc.mReadId << '[' << fc.mOffset << ']'
              << GeneralizedBase::bits2Char(fc.mContext.getPredecessors())
              << GeneralizedBase::bits2Char(fc.mContext.getSuccessors()); }
private:
    unsigned mReadId;
    unsigned mOffset : 11;
    unsigned mRCOffset : 11;
    unsigned mRC : 1;
    KMerContext mContext;
};
TRIVIALLY_SERIALIZABLE(FriendContext);
typedef SerfVec<FriendContext> Circle;

class FriendAlignFinderQ : public FriendAlignerImpl
{
public:
    FriendAlignFinderQ( String const& friendsCache )
    : mCircles(friendsCache), mIdx(idxName(friendsCache))
    {}

    void getAligns( size_t readId, Friends* pFriends )
    { pFriends->clear();
      ULongVec const& circVec = mIdx[readId];
      for ( size_t circleId : circVec )
        resolveCircle(readId,mCircles[circleId],pFriends);
      uniqueSortFriends(*pFriends); }

    static String idxName( String const& friendsCache )
    { return friendsCache.ReplaceExtension(".friends",".friendsidx"); }

    template <class ReadC, class QualC>
    static void storeFriends(  ReadC const& reads, QualC const& quals,
                               String const& friendsCache,
                               unsigned FF_MIN_QUAL, unsigned FF_COVERAGE,
                               unsigned K,
                               unsigned FF_MIN_FREQ, unsigned FF_MAX_FREQ,
                               int FF_VERBOSITY );
private:

    // What we put into the dictionary: a KMer, how many times we saw it,
    // the context in which we saw it, and the form in which the kmer was
    // initially presented (in particular, whether it's a palindrome or not).
    template <unsigned K>
    class HashEnt1 : public KMer<K>
    {
    public:
        HashEnt1() : mCount(0), mForm(CanonicalForm::FWD) {}
        HashEnt1( bvec const& bv, unsigned pos )
        : KMer<K>(bv.begin(pos)), mCount(0)
        { if ( !pos || pos+K==bv.size() )
            mContext = KMerContext::NNContext();
          else
            mContext = KMerContext(bv[pos-1],bv[pos+K]);
          mForm = this->getCanonicalForm();
          if ( mForm==CanonicalForm::REV )
          { this->rc(); mContext = mContext.rc(); } }

        HashEnt1& operator=( KMer<K> const& kmer )
        { static_cast<KMer<K>&>(*this) = kmer;
          mCount = 0; mContext = KMerContext();
          return *this; }

        // callback for hashset's apply function
        void operator()( HashEnt1 const& constEnt )
        { HashEnt1& ent = const_cast<HashEnt1&>(constEnt);
          ent.mCount += 1; ent.mContext |= mContext; }

        unsigned getCount() const { return mCount; }
        KMerContext getContext() const { return mContext; }
        void setContext( KMerContext context ) { mContext = context; }
        bool isRC() const { return mForm==CanonicalForm::REV; }
        bool isPalindrome() const { return mForm==CanonicalForm::PALINDROME; }

    private:
        unsigned mCount;
        KMerContext mContext;
        CanonicalForm mForm;
    };

    template <unsigned K>
    using HS1=HashSet<HashEnt1<K>,typename KMer<K>::Hasher,std::equal_to<KMer<K>>>;

    // Kmerization procedure to kmerize and enter into the dictionary each
    // high-quality kmer with appropriate multiplicity observed in a batch of reads.
    template <unsigned K, class ReadC, class QualC>
    class Proc1
    {
    public:
        Proc1( ReadC const& reads, QualC const& quals,
               unsigned batchSz, unsigned minQual, HS1<K>* pHS )
        : mReads(reads), mQuals(quals), mHS(*pHS), mBatchSz(batchSz),
          mMinQual(minQual) {}

        void operator()( size_t batchId )
        { size_t startId = std::min(mReads.size(),batchId*mBatchSz);
          size_t endId = std::min(mReads.size(),startId+mBatchSz);
          auto qItr = mQuals.begin(startId);
          auto bEnd = mReads.begin(endId);
          for ( auto bItr = mReads.begin(startId); bItr != bEnd; ++bItr, ++qItr )
          { qvec const& qv = *qItr;
            if ( qv.size() < K ) continue;
            bvec const& bv = *bItr;
            unsigned goodQuals = 0;
            unsigned off = 1u-K;
            bool first = true;
            HashEnt1<K> lastEnt;
            auto qvEnd = qv.end();
            for ( auto qvItr=qv.begin(); qvItr!=qvEnd; ++qvItr, ++off )
            { if ( *qvItr < mMinQual ) goodQuals = 0;
              else if ( ++goodQuals >= K )
              { HashEnt1<K> ent(bv,off);
                if ( ent.isPalindrome() ) continue;
                if ( !first ) mHS.apply(lastEnt,lastEnt);
                else { first = false; ent.setContext(KMerContext::NNContext()); }
                lastEnt = ent; } }
            if ( !first )
            { lastEnt.setContext(KMerContext::NNContext());
              mHS.apply(lastEnt,lastEnt); } } }

    private:
        ReadC mReads;
        QualC mQuals;
        HS1<K>& mHS;
        unsigned mBatchSz;
        unsigned mMinQual;
    };

    // Element type for a 2nd dictionary, which maps kmers to circle-of-friends IDs
    template <unsigned K>
    class HashEnt2 : public KMer<K>
    {
    public:
        HashEnt2() : mId(0) {}
        HashEnt2( KMer<K> const& kmer, size_t id ) : KMer<K>(kmer), mId(id) {}

        size_t getId() const { return mId; }
        void setId( size_t id ) { mId = id; }

    private:
        size_t mId;
    };

    template <unsigned K>
    using HS2=HashSet<HashEnt2<K>,typename KMer<K>::Hasher,std::equal_to<KMer<K>>>;

    // Procedure to re-kmerize the reads, and build the circles of friends that
    // share a high-quality kmer of appropriate multiplicity.  Builds the reverse
    // index at the same time.
    template <unsigned K, class ReadC, class QualC>
    class Proc2
    {
    public:
        Proc2( ReadC const& reads, QualC const& quals,
               HS2<K> const& hs, unsigned batchSz, unsigned minQual,
               MasterVec<Circle>* pCircles, std::vector<SpinLockedData>* pLocks,
               IncrementalWriter<ULongVec>* pCircIdx, std::mutex* pMutex,
               std::condition_variable* pCondVar )
        : mReads(reads), mQuals(quals), mHS(hs),
          mBatchSz(batchSz), mMinQual(minQual),
          mCircles(*pCircles), mLocks(*pLocks), mCircIdx(*pCircIdx),
          mMutex(*pMutex), mCondVar(*pCondVar) {}

        void operator()( size_t batchId )
        { size_t startId = std::min(mReads.size(),batchId*mBatchSz);
          size_t endId = std::min(mReads.size(),startId+mBatchSz);
          size_t rId = startId;
          auto qItr = mQuals.begin(rId);
          auto bEnd = mReads.begin(endId);
          HashEnt2<K> const* pEnt;
          VecULongVec circIdx;
          circIdx.reserve(endId-startId);
          ULongVec circleIds;
          for ( auto bItr = mReads.begin(rId); bItr != bEnd; ++bItr, ++qItr, ++rId )
          { qvec const& qv = *qItr;
            if ( qv.size() < K ) continue;
            bvec const& bv = *bItr;
            unsigned rcOff = bv.size()-K;
            unsigned goodQuals = 0;
            unsigned off = 1u-K;
            bool first = true;
            FriendContext fr;
            size_t circleId = 0;
            auto qvEnd = qv.end();
            for ( auto qvItr=qv.begin(); qvItr!=qvEnd; ++qvItr, ++off )
            { if ( *qvItr < mMinQual ) goodQuals = 0;
              else if ( ++goodQuals >= K )
              { HashEnt1<K> ent(bv,off);
                if ( (pEnt = mHS.lookup(ent)) )
                { if ( !first )
                  { SpinLocker lock(mLocks[circleId]);
                    mCircles[circleId].push_back(fr); }
                  else
                  { first = false; ent.setContext(KMerContext::NNContext()); }
                  fr = FriendContext(rId,off,rcOff-off,ent.isRC(),ent.getContext());
                  circleId = pEnt->getId();
                  circleIds.push_back(circleId); } } }
            if ( !first )
            { fr.setContext(KMerContext::NNContext());
              SpinLocker lock(mLocks[circleId]);
              mCircles[circleId].push_back(fr); }
            circIdx.push_back(circleIds); circleIds.clear(); }
          std::unique_lock<std::mutex> lock(mMutex);
          mCondVar.wait(lock,
               [this,startId](){return mCircIdx.getNElements()==startId;});
          mCircIdx.add(circIdx.begin(),circIdx.end());
          ForceAssertEq(mCircIdx.getNElements(),endId);
          mCondVar.notify_all(); }

    private:
        ReadC mReads;
        QualC mQuals;
        HS2<K> const& mHS;
        unsigned mBatchSz;
        unsigned mMinQual;
        MasterVec<Circle>& mCircles;
        std::vector<SpinLockedData>& mLocks;
        IncrementalWriter<ULongVec>& mCircIdx;
        std::mutex& mMutex;
        std::condition_variable& mCondVar;
    };

    static void resolveCircle( size_t rId, Circle const& circle, Friends* pFriends )
    { using std::lower_bound;
      auto sItr = lower_bound(circle.begin(),circle.end(),FriendContext(rId));
      auto sEnd = lower_bound(circle.begin(),circle.end(),FriendContext(rId+1));
      for ( auto self = sItr; self != sEnd; ++self )
      { for ( auto itr = circle.begin(); itr != sItr; ++itr )
          if ( self->isFriend(*itr) )
            pFriends->push_back(self->makeFriend(*itr));
        for ( auto itr=sEnd, end=circle.end(); itr != end; ++itr )
          if ( self->isFriend(*itr) )
            pFriends->push_back(self->makeFriend(*itr)); } }

    static void uniqueSortFriends( Friends& friends )
    { std::sort(friends.begin(),friends.end());
      friends.erase(std::unique(friends.begin(),friends.end()),friends.end()); }

    // A procedure to examine the reverse index (all the circles of friends in which
    // a read participates), and discard the redundant circles.  The first and last
    // circles are always preserved (the last one is swapped into the 2nd position),
    // but circles of friends implied by kmers in the middle of the read often don't
    // add any new friends and these can be safely discarded from the reverse index.
    class Proc3
    {
    public:
        Proc3( MasterVec<Circle> const& circles,
                VirtualMasterVec<ULongVec>& circIdx,
                size_t batchSz, size_t nReads,
                IncrementalWriter<ULongVec>* pCircIdxWriter, std::mutex* pMutex,
                std::condition_variable* pCondVar )
        : mCircles(circles), mCircIdx(circIdx), mBatchSz(batchSz), mNReads(nReads),
          mCircIdxWriter(*pCircIdxWriter), mMutex(*pMutex), mCondVar(*pCondVar)
        {}

        void operator()( size_t batchId )
        { size_t startId = std::min(mNReads,batchId*mBatchSz);
          size_t endId = std::min(mNReads,startId+mBatchSz);
          size_t rId = startId;
          Friends friends;
          VecULongVec circIdx;
          circIdx.reserve(endId-startId);
          while ( rId != endId )
          { ULongVec circVec = mCircIdx[rId];
            if ( circVec.size() > 2 )
            { using std::swap; swap(circVec[1],circVec.back());
              resolveCircle(rId,mCircles[circVec[0]],&friends);
              resolveCircle(rId,mCircles[circVec[1]],&friends);
              uniqueSortFriends(friends);
              size_t circSz = friends.size();
              auto itr = circVec.end();
              auto end = circVec.begin(1);
              while ( --itr != end )
              { resolveCircle(rId,mCircles[*itr],&friends);
                uniqueSortFriends(friends);
                if ( circSz == friends.size() ) circVec.erase(itr);
                else circSz = friends.size(); }
              friends.clear(); }
            circIdx.push_back(circVec);
            ++rId; }
          std::unique_lock<std::mutex> lock(mMutex);
          mCondVar.wait(lock,
              [this,startId](){return mCircIdxWriter.getNElements()==startId;});
          mCircIdxWriter.add(circIdx.begin(),circIdx.end());
          ForceAssertEq(mCircIdxWriter.getNElements(),endId);
          mCondVar.notify_all(); }

    private:
        MasterVec<Circle> const& mCircles;
        VirtualMasterVec<ULongVec> mCircIdx;
        size_t mBatchSz;
        size_t mNReads;
        IncrementalWriter<ULongVec>& mCircIdxWriter;
        std::mutex& mMutex;
        std::condition_variable& mCondVar;
    };

    template <unsigned K, class ReadC, class QualC>
    static void findFriends( ReadC const& reads, QualC const& quals,
                                String const& friendsCache,
                                unsigned minQ, unsigned coverage,
                                unsigned minFreq, unsigned maxFreq, int verbosity )
    {
        size_t kmerCount = reads.getKmerCount(K)/coverage;
        if ( verbosity )
            std::cout << Date() << ": Estimated dict size = "<<kmerCount<<std::endl;
        HS1<K>* pHS1 = new HS1<K>(kmerCount);

        if ( verbosity )
            std::cout << Date() << ": parsing" << std::endl;
        unsigned nThreads = getConfiguredNumThreads();
        unsigned nBatches = 5*nThreads;
        size_t batchSz = (reads.size()+nBatches-1)/nBatches;
        Proc1<K,ReadC,QualC> proc1(reads,quals,batchSz,minQ,pHS1);
        parallelFor(0u,nBatches,proc1,nThreads);

        if ( verbosity )
            std::cout << Date() << ": Actual dict size = "<<pHS1->size()<<std::endl;
        String dictFile = friendsCache.ReplaceExtension(".friends",".dict");
        if ( true )
        {
            size_t gLo = 0;
            size_t gHi = 0;
            size_t gContext = 0;
            BinaryIteratingWriter<std::vector<HashEnt2<K>>> bw(dictFile);
            kmerCount = 0;
            for ( auto const& hhs : *pHS1 )
                for ( auto const& ent : hhs )
                    if ( ent.getCount() < minFreq ) ++gLo;
                    else if ( ent.getCount() > maxFreq ) ++gHi;
                    else if ( ent.getContext().isSingleContext() ) ++gContext;
                    else
                    {
                         // store count in ent's id member, for now
                         bw.write(HashEnt2<K>(ent,ent.getCount()));
                         kmerCount += 1;
                    }
            if ( verbosity )
                std::cout << "gLo=" << gLo << " gHi=" << gHi << " gContext="
                            << gContext << std::endl;
            bw.close();
            delete pHS1;
            pHS1 = nullptr;
        }

        if ( verbosity )
            std::cout << Date() << ": Good kmer dict size = "<<kmerCount<<std::endl;
        HS2<K> hs2(kmerCount);
        MasterVec<Circle> circles(kmerCount);
        if ( true )
        {
            BinaryIteratingReader<std::vector<HashEnt2<K>>> br(dictFile);
            size_t idx = 0;
            HashEnt2<K> ent;
            while ( br.remaining() )
            {
                br.next(&ent);
                // id member is actually the count
                circles[idx].reserve(ent.getId());
                // reset id member to id
                ent.setId(idx++);
                hs2.insertUniqueValueNoLocking(std::move(ent));
            }
            unlink(dictFile.c_str());
        }

        if ( verbosity )
            std::cout << Date() << ": building circles" << std::endl;
        String const& idxName = FriendAlignFinderQ::idxName(friendsCache);
        String const& tmpIdxName = idxName + ".tmp";
        if ( true )
        {
            std::vector<SpinLockedData> locks(kmerCount);
            std::mutex mtx;
            std::condition_variable cndVar;
            IncrementalWriter<ULongVec> circIdxWriter(tmpIdxName,reads.size());
            Proc2<K,ReadC,QualC> proc2(reads,quals,hs2,batchSz,minQ,&circles,
                                        &locks,&circIdxWriter,&mtx,&cndVar);
            parallelFor(0u,nBatches,proc2,nThreads);
            circIdxWriter.close();
        }

        if ( verbosity )
            std::cout << Date() << ": sorting circles" << std::endl;
        size_t circleBatchSize = (circles.size()+nThreads-1)/nThreads;
        parallelForBatch(0ul,circles.size(),circleBatchSize,
                [&circles]( size_t idx )
                { Circle& circle = circles[idx];
                  std::sort(circle.begin(),circle.end()); },nThreads);

        if ( verbosity )
            std::cout << Date() << ": compressing circles index" << std::endl;
        if ( true )
        {
            VirtualMasterVec<ULongVec> circIdx(tmpIdxName);
            std::mutex mtx;
            std::condition_variable cndVar;
            IncrementalWriter<ULongVec> circIdxWriter(idxName,reads.size());
            Proc3 proc3(circles,circIdx,batchSz,reads.size(),
                                                &circIdxWriter,&mtx,&cndVar);
            parallelFor(0u,nBatches,proc3,nThreads);
            circIdxWriter.close();
            unlink(tmpIdxName.c_str());
        }

        if ( verbosity )
            std::cout << Date() << ": writing circles" << std::endl;
        circles.WriteAll(friendsCache);

        if ( verbosity )
            std::cout << Date() << ": done" << std::endl;
    }

    // this will eventually need to be a VirtualMasterVec, but it's not
    // thread-safe, and getAligns needs to be thread-safe.
    MasterVec<Circle> mCircles; // reads in a circle
    VecULongVec mIdx; // reverse index: circles for a read
};

template <unsigned K>
struct Serializability<FriendAlignFinderQ::HashEnt2<K>>
{ typedef TriviallySerializable type; };

template <class ReadC, class QualC>
void FriendAlignFinderQ::storeFriends( ReadC const& reads,
                                        QualC const& quals,
                                        String const& cache,
                                        unsigned minQ,
                                        unsigned cvrg,
                                        unsigned K,
                                        unsigned minFreq,
                                        unsigned maxFreq,
                                        int verbosity )
{
    switch ( K )
    {
    case 12u:
      findFriends<12u>(reads,quals,cache,minQ,cvrg,minFreq,maxFreq,verbosity);
      break;
    case 16u:
      findFriends<16u>(reads,quals,cache,minQ,cvrg,minFreq,maxFreq,verbosity);
      break;
    case 24u:
      findFriends<24u>(reads,quals,cache,minQ,cvrg,minFreq,maxFreq,verbosity);
      break;
    case 28u:
      findFriends<28u>(reads,quals,cache,minQ,cvrg,minFreq,maxFreq,verbosity);
      break;
    case 40u:
      findFriends<40u>(reads,quals,cache,minQ,cvrg,minFreq,maxFreq,verbosity);
      break;
    case 60u:
      findFriends<60u>(reads,quals,cache,minQ,cvrg,minFreq,maxFreq,verbosity);
      break;
    case 80u:
      findFriends<80u>(reads,quals,cache,minQ,cvrg,minFreq,maxFreq,verbosity);
      break;
    default:
      FatalErr("K=" << K << " not supported by FriendAlignFinderQ.");
    }
}

#endif /* PATHS_LONG_FRIENDALIGNFINDERQ_H_ */
