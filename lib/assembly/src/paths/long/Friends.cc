///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file Friends.cc
 * \author tsharpe
 * \date Apr 24, 2012
 *
 * \brief
 */
#include "paths/long/Friends.h"
#include "dna/CanonicalForm.h"
#include "feudal/HashSet.h"
#include "feudal/IncrementalWriter.h"
#include "feudal/Mempool.h"
#include "kmers/KMer.h"
#include "system/SysConf.h"
#include "system/WorklistN.h"
#include <algorithm>
#include <condition_variable>
#include <cmath>
#include <mutex>

namespace
{

unsigned const K = 25;
unsigned const MAX_KMER_FREQ = 100;

class KmerFriends : public KMer<K>
{
public:
    KmerFriends( MempoolAllocator<IdAndOrientation> const& alloc )
    : mFriends(alloc) {}

    KmerFriends& operator=( KMer<K> const& kmer )
    { static_cast<KMer<K>&>(*this) = kmer; mFriends.clear(); return *this; }

    bool addFriend( unsigned readId, CanonicalForm form )
    { if ( mFriends.size() == MAX_KMER_FREQ )
      { mFriends.clear(); return false; }
      if ( !mFriends.capacity() ) mFriends.reserve(1);
      else mFriends.reserve(MAX_KMER_FREQ);
      mFriends.push_back(IdAndOrientation(readId,form==CanonicalForm::REV));
      return true; }

    IAndOs::const_iterator fbegin() const { return mFriends.begin(); }
    IAndOs::const_iterator fend() const { return mFriends.end(); }

    class Factory
    {
    public:
        Factory( MempoolAllocator<IdAndOrientation> const& alloc )
        : mAlloc(alloc) {}

        template <class X>
        MempoolAllocator<X> alloc( X* ) const
        { return MempoolAllocator<X>(mAlloc); }

        KmerFriends* create( size_t nnn ) const
        { MempoolAllocator<KmerFriends> alloc(mAlloc);
          KmerFriends* result = alloc.allocate(nnn);
          for ( auto itr=result,end=result+nnn; itr != end; ++itr )
            new (itr) KmerFriends(mAlloc);
          return result; }

        void destroy( KmerFriends* pKF, size_t nnn ) const
        { for ( auto itr=pKF,end=pKF+nnn; itr != end; ++itr )
            itr->~KmerFriends();
          MempoolAllocator<KmerFriends> alloc(mAlloc);
          alloc.deallocate(pKF,nnn); }

    private:
        MempoolAllocator<IdAndOrientation> mAlloc;
    };
private:
    IAndOs mFriends;
};

class ReadGroupFriends
{
public:
    explicit ReadGroupFriends( size_t nnn )
    : mHash(nnn,KMer<K>::Hasher(),std::equal_to<KMer<K>>(),
                                            KmerFriends::Factory(mAlloc)) {}

    void processRange( vecbvec const& reads, size_t startId, size_t endId,
                        IncrementalWriter<IAndOs>& writer )
    { for ( size_t readId = startId; readId != endId; ++readId )
        addRead(reads[readId]);
      addFriends(reads);
      dumpFriends(startId,endId,writer); }

private:
    KmerFriends* lookup( KMer<K> const& kmer )
    { KmerFriends const* pKF = mHash.lookup(kmer);
      return const_cast<KmerFriends*>(pKF); }

    void addRead( bvec const& read )
    { if ( read.size() < K ) return;
      unsigned nKmers = read.size()-K+1;
      KMer<K> kmer(read.begin());
      mHash.add(kmer.isRev()?KMer<K>(kmer).rc():kmer);
      for ( auto itr = read.begin()+K, end = read.end(); itr != end; ++itr )
      { kmer.toSuccessor(*itr);
        mHash.add(kmer.isRev()?KMer<K>(kmer).rc():kmer); } }

    void addFriends( vecbvec const& reads )
    { size_t nReads = reads.size();
      for ( size_t readId = 0; readId != nReads; ++readId )
      { bvec const& read = reads[readId];
        if ( read.size() < K ) continue;
        KMer<K> kmer(read.begin());
        CanonicalForm form = kmer.getCanonicalForm();
        KmerFriends* pKF = lookup(form==CanonicalForm::REV ?
                                    KMer<K>(kmer).rc() : kmer);
        if ( pKF && !pKF->addFriend(readId,form) ) mHash.remove(*pKF);
        for ( auto itr = read.begin()+K, end = read.end(); itr != end; ++itr )
        { kmer.toSuccessor(*itr);
          form = kmer.getCanonicalForm();
          pKF = lookup(form==CanonicalForm::REV ? KMer<K>(kmer).rc() : kmer);
          if ( pKF && !pKF->addFriend(readId,form) ) mHash.remove(*pKF); } } }

    void dumpFriends( size_t startId, size_t endId,
                        IncrementalWriter<IAndOs>& writer )
    { OuterVec<IAndOs,MempoolAllocator<IdAndOrientation>> fVec(mAlloc);
      fVec.resize(endId-startId);
      for ( auto itr=mHash.begin(),end=mHash.end(); itr != end; ++itr )
        for ( KmerFriends const& kf : *itr )
          if ( std::distance(kf.fbegin(),kf.fend()) > 1 )
            for ( auto fItr = std::lower_bound(kf.fbegin(),kf.fend(),
                                 IdAndOrientation(startId)),
                  fEnd = std::lower_bound(fItr,kf.fend(),
                                 IdAndOrientation(endId));
                  fItr != fEnd; ++fItr )
              merge(kf.fbegin(),fItr,kf.fend(),fVec[fItr->getId()-startId]);

      std::unique_lock<std::mutex> lock(gWriteLock);
      gMyTurn.wait(lock,
              [startId,&writer](){ return startId==writer.getNElements();});
      writer.add(fVec.begin(),fVec.end());
      gMyTurn.notify_all(); }

    void merge( IAndOs::const_iterator itr, IAndOs::const_iterator itself,
                IAndOs::const_iterator end, IAndOs& friends )
    { friends.reserve(MAX_KMER_FREQ);
      auto oItr = friends.begin();
      auto oEnd = friends.end();
      bool isRC = itself->isRC();
      while ( itr != itself && oItr != friends.end() )
      { IdAndOrientation tmp(itr->getId(),itr->isRC()!=isRC);
        if ( tmp == *oItr ) ++itr;
        else if ( tmp < *oItr )
        { oItr = friends.insert(oItr,tmp); ++itr; }
        ++oItr; }

      if ( oItr == friends.end() )
      { while ( itr != itself )
        { IdAndOrientation tmp(itr->getId(),itr->isRC()!=isRC);
          friends.push_back(tmp); ++itr; }
        oItr = friends.end(); }

      itr = itself + 1;
      while ( itr != end && oItr != friends.end() )
      { IdAndOrientation tmp(itr->getId(),itr->isRC()!=isRC);
        if ( tmp == *oItr ) ++itr;
        else if ( tmp < *oItr )
        { oItr = friends.insert(oItr,tmp); ++itr; }
        ++oItr; }
      while ( itr != end )
      { IdAndOrientation tmp(itr->getId(),itr->isRC()!=isRC);
        friends.push_back(tmp); ++itr; } }

    MempoolOwner<IdAndOrientation> mAlloc;
    HashSet<KmerFriends,KMer<K>::Hasher,std::equal_to<KMer<K>>,

    KmerFriends::Factory> mHash;
    static std::mutex gWriteLock;
    static std::condition_variable gMyTurn;
};
std::mutex ReadGroupFriends::gWriteLock;
std::condition_variable ReadGroupFriends::gMyTurn;

}

// Here is some badly flawed code that emits lots of false positives.
// It's what we're using until we have time to find out why the correct code
// above changes results.
#include "kmers/ReadPatherDefs.h"
#include <map>
#include <set>
#include <vector>

namespace
{
class FriendFinder
{
public:
    FriendFinder()
    : mZMin(1.5), mMinUnipathLen(5u),
      mMinKmersOnUnipath(3u),
      mMinCommonUnipaths(3u) {}

    /// Find everyone's friends.
    void findFriends(
        IAndOsVec* pFriends, //output
        vecbvec const& reads, //the set of reads to explore for friendship
        unsigned coverage = 60u, //est. coverage, just for sizing spectrum
        unsigned nThreads = 0u,  //number of parallel threads to use
        int verbosity = 1 );

private:
    double mZMin;
    unsigned mMinUnipathLen;
    unsigned mMinKmersOnUnipath;
    unsigned mMinCommonUnipaths;
};


class DictRemovalKmerEater
{
public:
    DictRemovalKmerEater( KmerDict<K>* pDict ) : mpDict(pDict) {}

    // compiler-supplied copying and destructor are OK

    void operator()( KMer<K> const& kmer, KMerContext, size_t, size_t )
    { mpDict->remove(kmer); }

private:
    KmerDict<K>* mpDict;
};

size_t removeShortUnipaths( UnipathGraph<K> const& ug, size_t minLength,
                                KmerDict<K>* pDict )
{
    size_t badUnipathCount = 0;
    typedef HugeBVec::const_iterator BVItr;
    UnipathEdgeVec const& unipaths = ug.getAllEdges();
    size_t nnn = unipaths.size();
    DictRemovalKmerEater eater(pDict);
    for ( size_t idx = 0; idx < nnn; ++idx )
    {
        if ( unipaths[idx].getLength() < minLength )
        {
            badUnipathCount += 1;
            EdgeID unipathID(idx);
            KMer<K>::kmerizeIntoEater(ug.getBases(unipathID),
                                ug.getBasesEnd(unipathID),
                                eater,0);
        }
    }
    return unipaths.size()-badUnipathCount;
}

#if 0
void assessDictionary( vecbvec const& genome, unsigned numThreads,
                            KmerDict<K> const& dict )
{
    std::cout << Date() << ": Checking good kmers against genome."
            << std::endl;
    size_t genomeDictSize = genome.SizeSum();
    KmerDict<K> genomeDict( 5*genomeDictSize/4 );
    genomeDict.process(genome.begin(),genome.end(),numThreads,1);
    std::cout << Date( ) << ": there are " << genomeDict.size()
              << " genomic kmers (expected ~" << genomeDictSize << ")"
              << std::endl;

    size_t const GENOME_MAX_COUNT = 30;
    std::vector<size_t> kCounts(GENOME_MAX_COUNT+1);
    typedef KmerDict<K>::OCItr OCItr;
    typedef KmerDict<K>::ICItr ICItr;
    for ( OCItr oItr(dict.begin()),oEnd(dict.end()); oItr!=oEnd; ++oItr )
    {
        for ( ICItr itr(oItr->begin()), end(oItr->end()); itr!=end; ++itr )
        {
            size_t count = 0;
            KDef const* pKDef = genomeDict.lookup(*itr);
            if ( pKDef )
            {
                count = pKDef->getCount();
                if ( count > GENOME_MAX_COUNT )
                    count = GENOME_MAX_COUNT;
            }
            kCounts[count] += 1;
        }
    }
    std::cout << "True copy numbers of the good kmers:" << std::endl;
    for ( size_t idx = 0; idx < GENOME_MAX_COUNT; ++idx )
        std::cout << idx << "x\t" << kCounts[idx] << '\n';
    std::cout << ">=" << GENOME_MAX_COUNT << "x\t"
              << kCounts[GENOME_MAX_COUNT] << std::endl;
    size_t dictSize = dict.size();
    std::cout << std::setiosflags(ios::fixed) << std::setprecision(2)
              << kCounts[0]*100./dictSize
              << "% of the good kmers are actually errors." << std::endl;
    std::cout << std::setiosflags(ios::fixed) << std::setprecision(1)
              << (dictSize-kCounts[0])*100./genomeDict.size()
              << "% of the genomic kmers are present among the good kmers."
              << std::endl;
}
#endif

typedef std::vector<IdAndOrientation> ReadPath; //id is a unipathId
typedef std::vector<ReadPath> ReadPathVec; //index is a readId

class PathingKmerEater
{
public:
    PathingKmerEater( KmerDict<K> const& dict, unsigned minCount )
    : mDict(dict), mMinCount(minCount), mCurCount(0) {}

    // compiler-supplied copying and destructor are OK

    void operator()( KMer<K> const& kmer, KMerContext, size_t, size_t )
    { KmerDict<K>::Entry const* pEntry = mDict.findEntry(kmer);
      if ( pEntry )
      { IdAndOrientation iAndO(pEntry->getKDef().getEdgeID().val(),kmer.isRev());
        if ( iAndO != mLastIAndO )
        {
            mLastIAndO = iAndO;
            mCurCount = 0;
        }
        if ( ++mCurCount == mMinCount )
            mPath.push_back(iAndO); } }

    ReadPath const& getPath() const { return mPath; }
    void clear()
    { mPath.clear(); mCurCount = 0; mLastIAndO = IdAndOrientation(); }

private:
    KmerDict<K> const& mDict;
    unsigned mMinCount;
    unsigned mCurCount;
    IdAndOrientation mLastIAndO;
    ReadPath mPath;
};

class PathingProc
{
public:
    PathingProc( vecbvec const& reads, KmerDict<K> const& dict,
                unsigned minCount, Dotter& dotter, ReadPathVec* pPaths )
    : mReads(reads), mDotter(dotter), mpPaths(pPaths), mEater(dict,minCount)
    { size_t nBatches = mDotter.getNBatches();
      mBatchSize = (mReads.size()+nBatches-1)/nBatches; }

    // compiler-supplied copying and destructor is OK

    void operator()( size_t );

private:
    vecbvec const& mReads;
    Dotter& mDotter;
    ReadPathVec* mpPaths;
    PathingKmerEater mEater;
    size_t mBatchSize;
};

void PathingProc::operator()( size_t batchNo )
{
    size_t readId = std::min(batchNo*mBatchSize,mReads.size());
    size_t endId = std::min(readId+mBatchSize,mReads.size());
    typedef vecbvec::const_iterator Itr;
    Itr end(mReads.begin(endId));
    for ( Itr itr(mReads.begin(readId)); itr != end; ++itr )
    {
        KMer<K>::kmerizeIntoEater(itr->cbegin(),itr->cend(),mEater,readId);
        (*mpPaths)[readId++] = mEater.getPath();
        mEater.clear();
    }
    mDotter.batchDone();
}

// ReadId and the ordinal position of some edge on that read
class IdAndOrder
{
public:
    IdAndOrder() { mVal = ~0ul; }
    IdAndOrder( size_t id, size_t order )
    : mVal(id<<24)
    { AssertLe(id,MAX_ID); AssertLe(order,MAX_ORDER); mVal |= order; }

    // compiler-supplied copying and destructor are OK

    size_t getId() const { return mVal >> 24; }
    size_t getOrder() const { return mVal & MAX_ORDER; }

    friend bool operator<( IdAndOrder const& val1, IdAndOrder const& val2 )
    { return val1.mVal < val2.mVal; }

    friend bool operator==( IdAndOrder const& val1,
                                IdAndOrder const& val2 )
    { return val1.mVal == val2.mVal; }

    friend bool operator!=( IdAndOrder const& val1,
                                IdAndOrder const& val2 )
    { return val1.mVal != val2.mVal; }

private:
    static size_t const MAX_ID = 0xffffffffff;
    static size_t const MAX_ORDER = 0xffffff;
    unsigned long mVal;
};
size_t const IdAndOrder::MAX_ID;
size_t const IdAndOrder::MAX_ORDER;

typedef std::vector<IdAndOrder> ReadSet; //reads that touch some given unipath
typedef std::vector<ReadSet> ReadSetVec; //indexed by unipathId

class OrderAndCount
{
public:
    OrderAndCount() : mOrder(0), mCount(0) {}

    unsigned getOrder() const { return mOrder; }
    OrderAndCount& setOrder( unsigned order )
    { mOrder = order; return *this; }

    unsigned getCount() const { return mCount; }
    OrderAndCount& incrementCount()
    { mCount += 1;  return *this; }

private:
    unsigned mOrder;
    unsigned mCount;
};

class FriendProc
{
public:
    FriendProc( IAndOsVec& friendsVec, vecbvec const& reads,
            unsigned coverage, double zMin, unsigned minUnipathLen,
            unsigned minKmersOnUnipath, unsigned nThreads, int verbosity );

    void findFriends( unsigned nThreads, unsigned minUnipaths )
    { std::cout << Date() << ": finding friends" << std::endl;
      size_t nReads = mPaths.size();
      mFriendsVec.clear();
      mFriendsVec.resize(nReads);
      mDotter.setNBatches(nReads);
      Proc proc(*this,minUnipaths);
      parallelFor(0ul,nReads,proc,nThreads); }

    IAndOs findReadFriends( size_t readId, unsigned minUnipaths )
    { IAndOs result;
      Proc(*this,minUnipaths).findFriends(readId,result);
      return result; }

private:
    void findUnipaths( vecbvec const& reads, KmerDict<K> const& dict,
                            unsigned minKmersOnUnipath, unsigned nThreads,
                            const int verbosity )
    { std::cout << Date() << ": finding unipaths for each read" << std::endl;
      mPaths.resize(reads.size());
      size_t const batchSize = 100;
      size_t nReads = reads.size();
      size_t nBatches = (nReads+batchSize-1)/batchSize;
      Dotter dotter(nBatches);
      PathingProc proc(reads,dict,minKmersOnUnipath,dotter,&mPaths);
      parallelFor(0ul,nBatches,proc,nThreads); }

    void findReads( size_t nUnipaths, unsigned coverage, const int verbosity )
    {
      if ( verbosity >= 1 )
      {  std::cout << Date() << ": finding reads for each unipath" << std::endl;  }
      mReads.resize(nUnipaths); reserve(mReads,2*coverage);
      mRCReads.resize(nUnipaths); reserve(mRCReads,2*coverage);
      size_t nReads = mPaths.size();
      for ( size_t readId = 0; readId != nReads; ++readId )
      { ReadPath const& rp = mPaths[readId];
        typedef ReadPath::const_iterator RPItr;
        size_t order = 0;
        for ( RPItr itr(rp.begin()), end(rp.end()); itr != end; ++itr )
        { ReadSetVec& rs = itr->isRC() ? mRCReads : mReads;
          rs[itr->getId()].push_back(IdAndOrder(readId,order++)); } } }

    class Proc
    {
    public:
        Proc( FriendProc& fp, unsigned minCommonUnipaths )
        : mFP(fp), mMinCommonUnipaths(minCommonUnipaths) {}

        // compiler-supplied copying and destructor are OK

        void operator()( size_t idx )
        { findFriends(idx,mFP.mFriendsVec[idx]); }

        void findFriends( size_t idx, IAndOs& friends )
        { typedef std::map<size_t,OrderAndCount> FriendsMap;
          FriendsMap friendsMap;
          FriendsMap rcFriendsMap;

          ReadPath const& path = mFP.mPaths[idx];
          typedef ReadPath::const_iterator RPItr;
          for ( RPItr itr(path.begin()), end(path.end()); itr != end; ++itr )
          { ReadSet const* pF = &mFP.mReads[itr->getId()];
            ReadSet const* pR = &mFP.mRCReads[itr->getId()];
            using std::swap; if ( itr->isRC() ) swap(pF,pR);
            typedef ReadSet::const_iterator RSItr;
            for ( RSItr itr2(pF->begin()),end2(pF->end()); itr2!=end2; ++itr2 )
            { if ( itr2->getId() != idx )
              { OrderAndCount& onc = friendsMap[itr2->getId()];
                if ( !onc.getCount() || onc.getOrder() < itr2->getOrder() )
                  onc.setOrder(itr2->getOrder()).incrementCount(); } }
            for ( RSItr itr2(pR->begin()),end2(pR->end()); itr2!=end2; ++itr2 )
            { if ( itr2->getId() != idx )
              { OrderAndCount& onc = rcFriendsMap[itr2->getId()];
                if ( !onc.getCount() || onc.getOrder() > itr2->getOrder() )
                  onc.setOrder(itr2->getOrder()).incrementCount(); } } }

          friends.clear();
          friends.reserve(friendsMap.size()+rcFriendsMap.size());
          typedef FriendsMap::iterator MapItr;
          for ( MapItr itr(friendsMap.begin()), end(friendsMap.end());
                          itr != end; ++itr )
          { if ( itr->second.getCount() >= mMinCommonUnipaths )
              friends.push_back(IdAndOrientation(itr->first,false)); }
          for ( MapItr itr(rcFriendsMap.begin()), end(rcFriendsMap.end());
                          itr != end; ++itr )
          { if ( itr->second.getCount() >= mMinCommonUnipaths )
              friends.push_back(IdAndOrientation(itr->first,true)); }
          mFP.mDotter.batchDone(); }

    private:
        FriendProc& mFP;
        unsigned mMinCommonUnipaths;
    };
    friend class Proc;

    static void reserve( ReadSetVec& rsv, size_t size )
    { typedef ReadSetVec::iterator RSItr;
      for ( RSItr itr(rsv.begin()), end(rsv.end()); itr != end; ++itr )
        itr->reserve(size); }

    ReadPathVec mPaths;
    ReadSetVec mReads;
    ReadSetVec mRCReads;
    IAndOsVec& mFriendsVec;
    Dotter mDotter;
};

FriendProc::FriendProc( IAndOsVec& friendsVec,
                            vecbvec const& reads,
                            unsigned coverage,
                            double zMin,
                            unsigned minUnipathLen,
                            unsigned minKmersOnUnipath,
                            unsigned nThreads,
                            int verbosity )
: mFriendsVec(friendsVec)
{
    size_t allReadsLen = reads.SizeSum();
    size_t dictSize = allReadsLen/2;
    if ( verbosity >= 1 )
         std::cout << Date() << ": creating dictionary" << std::endl;
    KmerDict<K> dict( 5*dictSize/4 );
    dict.process(reads,verbosity,false,nThreads,100);
    if ( verbosity >= 2 )
    {
        std::cout << Date( ) << ": there are " << dict.size()
                << " kmers (expected ~" << dictSize << ")" << std::endl;
    }

    if ( verbosity >= 1 )
         std::cout << Date() << ": building kmer spectrum" << std::endl;
    const int spec_cov_mult = 3;
    Spectrum<K> spectrum(dict,spec_cov_mult*coverage);
    if ( verbosity >= 3 )
    {
        std::cout << "Kmer spectrum:" << std::endl;
        std::cout << spectrum << std::endl;
    }
    size_t peakIdx = spectrum.getFirstPeakIndex();
    size_t troughIdx = spectrum.getFirstTroughIndex();
    if ( verbosity >= 2 )
    {
        std::cout << Date( ) << ": trough has " << spectrum[troughIdx]
                 << " kmers at " << troughIdx << "x\n" << Date()
                 << ": peak has " << spectrum[peakIdx] << " kmers at "
                 << peakIdx << 'x' << std::endl;
    }

    if ( !spectrum.isFirstPeakValid() )
        FatalErr("Unable to find CN=1 peak.");
    if ( !spectrum.isPeakDistinct() )
    {   cout << "\n";
        PRINT3( peakIdx, sqrt(peakIdx), troughIdx );
        cout << "Warning: insufficient peak/trough separation.\n\n";   }

    peakIdx = peakIdx - std::min(peakIdx,size_t(zMin*sqrt(peakIdx)));
    if ( peakIdx <= troughIdx )
    {
        peakIdx = troughIdx + 1;
        std::cout << "Warning: Z too large given the peak to trough distance.  "
                "Truncating above the trough." << std::endl;
    }

    double const MIN_GOOD_FRACTION = .1;
    if ( spectrum.goodFract() < MIN_GOOD_FRACTION )
        FatalErr("There's almost nothing but errors here.");

    if ( verbosity >= 1 )
    {   std::cout << Date() << ": discarding error kmers (less than " << peakIdx
                  << "x)" << std::endl;    }
    dictSize = spectrum.sum(peakIdx);
    dict.clean(KmerDict<K>::BadKmerCountFunctor(peakIdx));
    AssertEq(dictSize,dict.size());
    if ( verbosity >= 2 )
    {
        std::cout << Date( ) << ": there are " << dictSize << " good kmers"
                    << std::endl;
    }

    UnipathGraph<K> ug(dict, verbosity);

    if ( verbosity >= 1 )
    {    std::cout << Date()
                   << ": removing kmers for short unipaths from dictionary"
                   << std::endl;    }
    size_t unipathCount = removeShortUnipaths(ug,minUnipathLen,&dict);
    if ( verbosity >= 2 )
    {
        std::cout << Date( ) << ": there are " << unipathCount
                   << " unipaths with length >= " << minUnipathLen << std::endl;
        dictSize = dict.size();
        std::cout << Date( ) << ": there are now " << dictSize << " good kmers"
                    << std::endl;
    }

    findUnipaths(reads,dict,minKmersOnUnipath,nThreads,verbosity);
    findReads(ug.getNEdges(),coverage,verbosity);
}

void FriendFinder::findFriends( IAndOsVec* pFriendsVec,
                                     vecbvec const& reads,
                                     unsigned coverage, unsigned nThreads,
                                     int verbosity )
{
    nThreads = boundNumThreads(nThreads);
    FriendProc fp(*pFriendsVec,reads,coverage,mZMin,mMinUnipathLen,
                    mMinKmersOnUnipath,nThreads,verbosity);
    fp.findFriends(nThreads,mMinCommonUnipaths);
    std::cout << Date() << ": done" << std::endl;
}

} // end of anonymous namespace


void FindFriends( vecbvec const& reads, String const& friendsFile,
                    bool useBrokenVersion )
{
    if ( useBrokenVersion )
    {
        FriendFinder ff;
        IAndOsVec fv;
        ff.findFriends(&fv,reads);
        fv.WriteAll(friendsFile);
        return;
    }
    size_t nReads = reads.size();
    IncrementalWriter<IAndOs> writer(friendsFile,nReads);
    size_t nKmers = reads.getKmerCount(K);
    size_t needed = 1.5*nKmers*(sizeof(KmerFriends)+
                                MAX_KMER_FREQ/2*sizeof(IdAndOrientation)) +
                nReads*(sizeof(IAndOs)+MAX_KMER_FREQ*sizeof(IdAndOrientation));
    size_t nPasses = ceil(sqrt(needed/(.8*MemAvailable())));
    size_t nThreads = getConfiguredNumThreads();
    size_t nBatches = nPasses * nThreads;
    size_t rgfSize = (nKmers+nBatches-1)/nBatches;
    size_t batchSize = (nReads+nBatches-1)/nBatches;
    parallelFor(0ul,nBatches,
            [&reads,batchSize,rgfSize,&writer]( size_t batchNo )
            { size_t start = std::min(reads.size(),batchNo*batchSize);
              size_t end = std::min(reads.size(),start+batchSize);
              ReadGroupFriends rgf(rgfSize);
              rgf.processRange(reads,start,end,writer); },nThreads);
    writer.close();
}

#include "feudal/SmallVecDefs.h"
template class SmallVec<IdAndOrientation,MempoolAllocator<IdAndOrientation> >;

#include "feudal/OuterVecDefs.h"
template class OuterVec<IAndOs>;
