///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * MapReduceEngine.h
 *
 *  Created on: Aug 15, 2013
 *      Author: tsharpe
 */

#ifndef MAPREDUCEENGINE_H_
#define MAPREDUCEENGINE_H_

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"

#include "system/SysConf.h"
#include "system/System.h"
#include "system/Thread.h"
#include <algorithm>
#include <condition_variable>
#include <cstdlib>
#include <functional>
#include <memory>
#include <mutex>
#include <omp.h>
#include <thread>
#include <vector>
#include <utility>
#include "10X/Martian.h"

// Itr is an iterator over the input data type:  e.g., an iterator over bvec's.
//    Don't ignore the possibility for your Impl class to hold the actual input
//    data as a vector or vectors, in which case Itr can be an index number.
// Key is the type to which the input is mapped:  e.g., a KMer<K>.
//    Key must be hashable by a functor, Hash, taking a Key and returning a
//    size_t.  Key must be sortable by a comparator, Comp, which by default will
//    be std::less<Key>.  It must be movable, or at least copyable.
// Impl is the type that does the mapping and reducing:  e.g., a kmer spectrum
//    calculator.  It must be copy-constructable, and must manage any storage
//    associated with the output from reduction in a thread-safe manner.
//    It must have two methods:
//    void map( Itr, OItr ), which maps the dereferenced Itr to
//      any number of Key's, which it stuffs into an output iterator, OItr.
//    void reduce( Key* key1, Key* key2 ), which does whatever it wants
//      to with the information that the range key1..key2 has equal keys.
// Here's a sample implementation for producing a kmer spectrum from a vecbvec:
//    unsigned const K = 25u;
//    typedef KMer<K> Kmer;
//    typedef std::vector<std::atomic<size_t>> Spectrum;
//    typedef vecbvec::const_iterator Itr;
//
//    class KSpecImpl
//    {
//    public:
//      KSpecImpl( Spectrum* pOut )
//      : mpOut(pOut)
//      { mBins.resize(mpOut->size(),0ul); }
//
//      ~KSpecImpl()
//      { auto oItr = mpOut->begin();
//        for ( auto itr=mBins.begin(),end=mBins.end(); itr != end; ++itr,++oItr )
//          oItr->fetch_add(*itr); }
//
//      template <class OItr>
//      void map( Itr itr, OItr oItr )
//      { Kmer::kmerize(itr->cbegin(),itr->cend(),oItr); }
//
//      void reduce( Kmer const* beg, Kmer const* end )
//      { size_t count = end-beg;
//        if ( --count >= mBins.size() ) count = mBins.size()-1;
//        mBins[count] += 1; }
//
//      Kmer* overflow( Kmer* beg, Kmer* end )
//      { return beg+std::min(size_t(end-beg),mBins.size()); }
//
//    private:
//      Spectrum* mpOut;
//      std::vector<size_t> mBins;
//    };
//
//    typedef MapReduceEngine<KSpecImpl,Kmer,Kmer::Hasher> MRE;
//    Spectrum spectrum(1000);
//    KSpecImpl impl(&spectrum);
//    MRE mre(impl);
//    vecbvec reads(READS);
//    mre.process(reads.getKmerCount(K),reads.cbegin(),reads.cend(),true);
//
template <class Impl, class Key, class Hash, class Comp=std::less<Key>>
class MapReduceEngine
{
    struct KeyBounds
    { Key* mpStart; Key* mpCur; Key* mpEnd; };

    class KeyBuf
    {
    public:
        static size_t const EMPTY = ~0ul;

        KeyBuf( size_t nBatches, size_t nKeysPerBatch )
        : mBounds(new KeyBounds[nBatches]), mEnd(mBounds+nBatches)
        { std::allocator<Key> alloc;
          Key* pKeys;
          if ( WeAreTrackingSomeMemory( ) )
          {    
               #pragma omp critical
               {    double memoryFrac = 0.9;
                    int64_t maxMem = MemAvailable(memoryFrac);
                    cout << Date( ) << ": KeyBuf constructor sees " 
                         << ToStringAddCommas(maxMem) << " bytes available" << endl;
                    cout << Date( ) << ": /proc/meminfo says "
                         << ToStringAddCommas( MemReallyAvailable( ) )
                         << " bytes available" << endl;
                    int64_t req = int64_t(sizeof(Key)) 
                         * int64_t(nBatches) * int64_t(nKeysPerBatch);
                    cout << Date( ) << ": requesting " 
                         << ToStringAddCommas(req) << " bytes" << endl;
                    pKeys = alloc.allocate(nBatches*nKeysPerBatch);
                    maxMem = MemAvailable(memoryFrac);
                    cout << Date( ) << ": now see " << ToStringAddCommas(maxMem) 
                         << " bytes available" << endl;
                    cout << Date( ) << ": /proc/meminfo says "
                         << ToStringAddCommas( MemReallyAvailable( ) )
                         << " bytes available" << endl;    }    }
          else pKeys = alloc.allocate(nBatches*nKeysPerBatch);
          for ( KeyBounds* itr = mBounds; itr != mEnd; ++itr )
          { itr->mpStart = itr->mpCur = pKeys;
            pKeys += nKeysPerBatch; itr->mpEnd = pKeys; } }

        KeyBuf( KeyBuf const& ) = delete;
        KeyBuf& operator=( KeyBuf const& ) = delete;

        ~KeyBuf()
        { for ( KeyBounds* itr=mBounds; itr != mEnd; ++itr )
            while ( itr->mpCur != itr->mpStart )
              (--itr->mpCur)->~Key();
          std::allocator<Key> alloc;
          Key* start = mBounds[0].mpStart;
          alloc.deallocate(start,mEnd[-1].mpEnd-start);
          delete [] mBounds; }

        KeyBounds& getBounds( size_t batch )
        { return mBounds[batch]; }

        bool addKey( Key const& key, size_t batch )
        { KeyBounds& kb = mBounds[batch];
          if ( kb.mpCur == kb.mpEnd )
            return false;
          new (kb.mpCur++) Key(key);
          return true; }

        size_t getMaxBatchSize()
        { size_t result = 0;
          for ( KeyBounds* itr=mBounds; itr != mEnd; ++itr )
              result = std::max(result,size_t(itr->mpCur-itr->mpStart));
          return result; }

        void move( KeyBounds& from, size_t batch )
        { KeyBounds& to = mBounds[batch];
          ForceAssertEq(to.mpStart,to.mpCur);
          Key* end = from.mpCur;
          for ( Key* itr = from.mpStart; itr != end; ++itr )
          { new (to.mpCur++) Key(std::move(*itr)); itr->~Key(); }
          from.mpCur = from.mpStart; }

        Key* begin() const { return mBounds->mpStart; }
        Key* compress()
        { KeyBounds* itr = mBounds;
          Key* beg = itr->mpStart;
          Key* end = itr->mpCur;
          itr->mpCur = itr->mpStart;
          KeyBounds* last = mEnd-1;
          ForceAssertEq(last->mpStart,last->mpCur);
          while ( ++itr != last )
          { if ( itr->mpStart == end )
            { end = itr->mpCur; itr->mpCur = itr->mpStart; continue; }
            for ( Key* pKey=itr->mpStart; pKey != itr->mpCur; ++pKey )
            { new (end++) Key(std::move(*pKey)); pKey->~Key(); }
            itr->mpCur = itr->mpStart; }
          return end; }

    private:
        KeyBounds* mBounds;
        KeyBounds* mEnd;
    };

    enum class Cycle { INIT, MAP, SWIZZLE, REDUCE, EXIT };

    class Status
    {
    public:
        Status( size_t nThreads, size_t nPasses, size_t nKeysPerBatch )
        : mNThreads(nThreads), mNPasses(nPasses), mNKeysPerBatch(nKeysPerBatch),
          mNTotalBatches(mNThreads*mNPasses), mPass(0), mOK(true),
          mpKeyBufs(new KeyBuf*[nThreads]), mDoneCount(0),
          mCycle(Cycle::INIT), mSwizDone(0), mSwizGen(0), mNOverflows(0)
        {}

        ~Status()
        { delete [] mpKeyBufs; }

        Status( Status const& ) = delete;
        Status& operator=( Status const& ) = delete;

        size_t getNThreads() const { return mNThreads; }
        size_t getNPasses() const { return mNPasses; }
        size_t getNTotalBatches() const { return mNTotalBatches; }
        size_t getNKeysPerBatch() const { return mNKeysPerBatch; }

        size_t getNOverflows() const { return mNOverflows; }
        void incrementNOverflows() { mNOverflows += 1; }

        size_t getPass() const { return mPass; }
        void setPass( size_t pass ) { mPass = pass; }

        KeyBuf& getKeyBuf( size_t thread ) { return *mpKeyBufs[thread]; }
        void setKeyBuf( size_t thread, KeyBuf* pKeyBuf )
        { mpKeyBufs[thread] = pKeyBuf; }

        void fail() { mOK = false; }
        bool OK() const { return mOK; }

        // notify everyone waiting for the next cycle to begin
        bool setCycle( Cycle cycle )
        { if ( !mOK && cycle != Cycle::EXIT ) return false;
          waitForAllDone();
          std::unique_lock<std::mutex> lock(mCycleMtx);
          mCycle = cycle;
          mCycleCV.notify_all();
          return true; }

        // wait until it's time to do a map cycle
        bool waitForMap()
        { done();
          std::unique_lock<std::mutex> lock(mCycleMtx);
          while ( mCycle != Cycle::MAP && mCycle != Cycle::EXIT )
              mCycleCV.wait(lock);
          return mCycle == Cycle::MAP; }

        // wait until it's time to do a swizzle cycle
        bool waitForSwizzle()
        { done();
          std::unique_lock<std::mutex> lock(mCycleMtx);
          while ( mCycle != Cycle::SWIZZLE && mCycle != Cycle::EXIT )
              mCycleCV.wait(lock);
          return mCycle == Cycle::SWIZZLE; }

        // wait until it's time to do a reduce cycle
        bool waitForReduce()
        { done();
          std::unique_lock<std::mutex> lock(mCycleMtx);
          while ( mCycle != Cycle::REDUCE && mCycle != Cycle::EXIT )
              mCycleCV.wait(lock);
          return mCycle == Cycle::REDUCE; }

        // thread barrier -- wait until everyone is done with the current
        // round of swizzling
        void swizDone()
        { std::unique_lock<std::mutex> lock(mSwizMtx);
          size_t swizGen = mSwizGen;
          if ( ++mSwizDone == mNThreads )
          { mSwizGen += 1; mSwizDone = 0; mSwizCV.notify_all(); }
          while ( swizGen == mSwizGen )
              mSwizCV.wait(lock); }

    private:
        void done()
        { std::unique_lock<std::mutex> lock(mDoneMtx);
          if ( ++mDoneCount == mNThreads ) mDoneCV.notify_one(); }

        void waitForAllDone()
        { std::unique_lock<std::mutex> lock(mDoneMtx);
          while ( mDoneCount != mNThreads )
              mDoneCV.wait(lock);
          mDoneCount = 0; }

        size_t mNThreads;
        size_t mNPasses;
        size_t mNKeysPerBatch;
        size_t mNTotalBatches;
        size_t mPass;
        bool mOK;
        KeyBuf** mpKeyBufs;
        std::mutex mDoneMtx;
        std::condition_variable mDoneCV;
        size_t mDoneCount;
        std::mutex mCycleMtx;
        std::condition_variable mCycleCV;
        Cycle mCycle;
        std::mutex mSwizMtx;
        std::condition_variable mSwizCV;
        size_t mSwizDone;
        size_t mSwizGen;
        std::atomic_size_t mNOverflows;
    };

    class OItr
    {
    public:
        OItr( MapReduceEngine& mre, Status& status, KeyBuf* pKeyBuf )
        : mMRE(mre), mStatus(status), mpKeyBuf(pKeyBuf) {}

        OItr& operator*() { return *this; }
        OItr& operator++() { return *this; }
        OItr& operator++(int) { return *this; }

        Key const& operator=( Key const& key )
        { size_t batch = mMRE.mHasher(key)%mStatus.getNTotalBatches();
          size_t nThreads = mStatus.getNThreads();
          if ( batch/nThreads == mStatus.getPass() )
          { size_t bin = batch%nThreads;
            if ( !mpKeyBuf->addKey(key,bin) )
            { mStatus.incrementNOverflows();
              KeyBounds& kb = mpKeyBuf->getBounds(bin);
              kb.mpCur = mMRE.overflow(kb.mpStart,kb.mpCur);
              if ( !mpKeyBuf->addKey(key,bin) )
                mStatus.fail(); } }
          return key; }

    private:
        MapReduceEngine& mMRE;
        Status& mStatus;
        KeyBuf* mpKeyBuf;
    };

    template <class Itr>
    class Client
    {
    public:
        Client( size_t thread, MapReduceEngine& mre, Status& status,
                    Itr beg, Itr end )
        : mThread(thread), mMRE(mre), mStatus(status)
        { size_t nTotInputs = end - beg;
          size_t nThreads = mStatus.getNThreads();
          size_t nInputsPerThread = (nTotInputs+nThreads-1)/nThreads;
          mBeg = beg+std::min(nTotInputs,mThread*nInputsPerThread);
          mEnd = beg+std::min(nTotInputs,(mThread+1)*nInputsPerThread); }

        void operator()()
        {
            size_t nThreads = mStatus.getNThreads();
            KeyBuf kb(nThreads+1,mStatus.getNKeysPerBatch());
            size_t nSwizzles = (nThreads - 1)/2;
            mStatus.setKeyBuf(mThread,&kb);

            while ( mStatus.waitForMap() )
            {
                // this must be inside loop because it reads the pass #.
                OItr oItr(mMRE,mStatus,&kb);

                for ( Itr itr = mBeg; mStatus.OK() && itr != mEnd; ++itr )
                    mMRE.mImpl.map(itr,oItr);

                if ( !mStatus.waitForSwizzle() )
                    break;
                for ( size_t swizzle = 0; swizzle != nSwizzles; ++swizzle )
                {
                    size_t batch = (mThread+swizzle+1) % nThreads;
                    kb.move( kb.getBounds(batch), nThreads );
                    kb.move( mStatus.getKeyBuf(batch).getBounds(mThread), batch );
                    mStatus.swizDone();

                    batch = (mThread+nThreads-swizzle-1) % nThreads;
                    KeyBuf& kb2 = mStatus.getKeyBuf(batch);
                    kb.move( kb2.getBounds(nThreads), batch );
                    mStatus.swizDone();
                }
                if ( !(nThreads & 1) )
                {
                    size_t batch = (mThread+nSwizzles+1) % nThreads;
                    kb.move( kb.getBounds(batch), nThreads );
                    mStatus.swizDone();

                    KeyBuf& kb2 = mStatus.getKeyBuf(batch);
                    kb.move( kb2.getBounds(nThreads), batch );
                    mStatus.swizDone();
                }

                if ( !mStatus.waitForReduce() )
                    break;
                mMRE.reduce(kb.begin(),kb.compress());
            }
        }

    private:
        size_t mThread;
        MapReduceEngine mMRE;
        Status& mStatus;
        Itr mBeg;
        Itr mEnd;
    };

public:
    MapReduceEngine( Impl const& impl=Impl(),
                        Hash const& hasher=Hash(),
                        Comp const& comparator=Comp() )
    : mImpl(impl), mHasher(hasher), mComparator(comparator), mFailed(false)
    {}

    ~MapReduceEngine()
    { if ( mFailed ) FatalErr("Map/Reduce operation has failed."); }

    enum class VERBOSITY { SILENT, QUIET, NOISY };

    // nKs must be an upper bound on the number of Ks produced by mapping the
    // entire input set.  We assume that largish subsets of the input are
    // linearish in their production of Ks.
    template <class Itr>
    bool run( size_t nKs, Itr beg, Itr end,
                VERBOSITY verbose=VERBOSITY::SILENT,
			 double memoryFrac = .9,
			 double meanHashUsage=.8 )
    {
	ForceAssertLe( meanHashUsage, 1.0 );
	ForceAssertGt( meanHashUsage, 0.0 );

        if ( inParallelSection() )
        {
            runSingleThreaded(nKs,beg,end);
            return true;
        }

        omp_set_num_threads(1);

        // these three consts could be made into an arguments, if necessary
        //size_t const minPasses = 0;

        size_t maxMem;
        if ( WeAreTrackingSomeMemory( ) )
        {
               #pragma omp critical
               {    maxMem = MemAvailable(memoryFrac);
                    cout << Date( ) << ": MapReduceEngine::run sees " 
                         << ToStringAddCommas(maxMem) << " bytes available" << endl;
                    cout << Date( ) << ": /proc/meminfo says "
                         << ToStringAddCommas( MemReallyAvailable( ) )
                         << " bytes available" << endl << endl;    }    }
        else maxMem = MemAvailable(memoryFrac);

        size_t const nThreads = getConfiguredNumThreads();

        size_t maxKs = meanHashUsage*nThreads/(nThreads+1)*maxMem/sizeof(Key);
        if ( !maxKs ) {
            PRINT4(nThreads,maxMem,maxKs,meanHashUsage);
            Martian::exit("Insufficient memory.");
        }

        size_t nPasses = (nKs+maxKs-1)/maxKs;
        if ( nPasses == 0 )
            FatalErr("No work to do. The calling code should watch for this case");
        size_t const MAX_PASSES = 500;
        if ( nPasses > MAX_PASSES )
        {   cout << "\nWell this is unfortunate.  It seems that there is "
                 << "unsufficient memory.\n"
                 << "We would need " << nPasses << " passes, but the maximum "
                 << "allowed here is " << MAX_PASSES << endl;
            cout << "Giving up." << endl << endl;
            _exit(1);    }

        //if ( nPasses < minPasses )
        //    nPasses = minPasses;
        if ( verbose != VERBOSITY::SILENT )
            std::cout << Date( ) << ": MapReduce needs " << nPasses << " passes." << std::endl;

        size_t meanKsPerBatch = nKs/nPasses/nThreads/nThreads;
        size_t nKsPerBatch = meanKsPerBatch/meanHashUsage;
        size_t memUsed = nKsPerBatch*nThreads*(nThreads+1)*sizeof(Key);
        nKsPerBatch = nKsPerBatch*std::min(1.*maxMem/memUsed,2.);
        size_t const MIN_BATCH_SIZE = 100000;
        if ( nKsPerBatch < MIN_BATCH_SIZE )
            nKsPerBatch = MIN_BATCH_SIZE;
        if ( verbose != VERBOSITY::SILENT )
            std::cout << "Expect " << meanKsPerBatch << " keys per batch.\n"
            "Provide " << nKsPerBatch << " keys per batch." << std::endl;

        Status status(nThreads,nPasses,nKsPerBatch);
        std::thread* threads = new std::thread[nThreads];
        for ( size_t thread = 0; thread != nThreads; ++thread )
            threads[thread] = std::thread(Client<Itr>(thread,*this,status,beg,end));

        size_t pass;
        for ( pass = 0; pass != nPasses; ++pass )
        {
            status.setPass(pass);
            if ( !status.setCycle(Cycle::MAP) ) break;
            if ( verbose != VERBOSITY::SILENT )
            {    if ( pass > 0 && pass % 5 == 0 ) cout << " ";
                 cout << "." << flush;    }
            if ( verbose == VERBOSITY::NOISY )
            {
                if ( pass ) std::cout << ".\n";
                std::cout<<Date()<<" Pass "<<pass+1<<" parse,"<< std::flush;
            }
            if ( !status.setCycle(Cycle::SWIZZLE) )
                break;
            if ( verbose == VERBOSITY::NOISY )
            {
#if 0
                size_t maxBatchSz = 0;
                for ( size_t thread = 0; thread != nThreads; ++thread )
                    maxBatchSz = std::max(maxBatchSz,
                                    status.getKeyBuf(thread).getMaxBatchSize());
                std::cout << ' ' << 1.*maxBatchSz/nKsPerBatch;
#endif
                std::cout << " swizzle," << std::flush;
            }
            if ( !status.setCycle(Cycle::REDUCE) )
                break;
            if ( verbose == VERBOSITY::NOISY )
                std::cout << " reduce" << std::flush;
        }
        cout << endl;
        status.setCycle(Cycle::EXIT);
        if ( verbose == VERBOSITY::NOISY )
            std::cout << ".\n";

        for ( size_t thread = 0; thread != nThreads; ++thread )
            threads[thread].join();

        delete [] threads;

        if ( (mFailed = !status.OK()) )
            std::cout << "Map/Reduce operation has failed at pass "
                        << pass << std::endl;

        if ( verbose != VERBOSITY::SILENT )
        {
            size_t nOverflows = status.getNOverflows();
            if ( nOverflows )
                std::cout << "There were " << nOverflows
                            << " buffer overflows." << std::endl;
        }

        omp_set_num_threads(getConfiguredNumThreads());
        cout << Date( ) << ": MapReduceEngine::run complete" << endl;
        return status.OK();
    }

private:
    template <class Itr>
    void runSingleThreaded( size_t nKs, Itr itr, Itr end )
    { std::vector<Key> keys;
      keys.reserve(nKs);
      auto oItr = std::back_inserter(keys);
      for ( ; itr != end; ++itr )
        mImpl.map(itr,oItr);
      if ( !keys.empty() )
      { Key* beg = &keys.front();
        reduce(beg,beg+keys.size()); } }

    Key* overflow( Key* beg, Key* end )
    { std::sort(beg,end,mComparator);
      Key* oItr = beg;
      while ( beg != end )
      { Key* itr2 = beg;
        while ( ++itr2 != end )
          if ( mComparator(*beg,*itr2) )
            break;
        oItr = moveKeys(beg,mImpl.overflow(beg,itr2),oItr);
        beg = itr2; }
      while ( end != oItr ) (--end)->~Key();
      return oItr; }

    void reduce( Key* beg, Key* end )
    { std::sort(beg,end,mComparator);
      Key* itr = beg;
      while ( itr != end )
      { Key* itr2 = itr;
        while ( ++itr2 != end )
          if ( mComparator(*itr,*itr2) )
            break;
        mImpl.reduce(itr,itr2);
        itr = itr2; }
      while ( end != beg ) (--end)->~Key(); }

    static Key* moveKeys( Key* itr, Key* end, Key* oItr )
    { if ( itr == oItr ) return end;
      while ( itr != end ) *oItr++ = std::move(*itr++);
      return oItr; }

    Impl mImpl;
    Hash mHasher;
    Comp mComparator;
    bool mFailed;
};


#endif /* MAPREDUCEENGINE_H_ */
