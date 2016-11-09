///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file ReadPatherDefs.h
 * \author tsharpe
 * \date May 4, 2010
 *
 * \brief
 */
#ifndef KMERS_READ_PATHER_DEFS_H_
#define KMERS_READ_PATHER_DEFS_H_


#define LOG(x) ((void)0)
//#define LOG(x) std::cout << x << std::endl

#include "kmers/ReadPather.h"
#include "Vec.h"
#include "feudal/ChunkDumper.h"
#include "feudal/VirtualMasterVec.h"
#include "system/SpinLockedData.h"
#include "system/SysConf.h"
#include "system/System.h"
#include "system/WorklistN.h"
#include <algorithm>
#include <cmath>
#include <ostream>
#include <vector>
#include <time.h>
#include <unistd.h>

namespace
{
// A worklist processor class to kmerize reads.
// A KmerEater is a functor with the following signature:
// operator()( KMer<K> const&, KMerContext, size_t readId, size_t readOffset );
// A separate eater is used for each thread.
// The Itr class is a random access iterator over bvec's.
template <unsigned K, class KmerEater, class Itr>
class KmerizationProcessor
{
public:
    KmerizationProcessor( Itr const& start, Itr const& stop,
                            KmerEater const& eater, Dotter* pDotter )
    : mReads(start), mBatchSize(1000), mEater(eater), mpDotter(pDotter)
    { using std::distance; mNReads = distance(start,stop);
      if ( pDotter )
      { size_t nBatches = pDotter->getNBatches();
        if ( !nBatches )
          pDotter->setNBatches((mNReads+mBatchSize-1)/mBatchSize);
        else
          mBatchSize = (mNReads+nBatches-1)/nBatches; } }

    // compiler-supplied copying and destructor is OK

    void operator()( size_t batchNo );

    KmerEater const& getEater() const { return mEater; }
    Dotter* getDotter() const { return mpDotter; }

private:
    Itr mReads;
    size_t mNReads;
    size_t mBatchSize;
    KmerEater mEater;
    Dotter* mpDotter;
};

template <unsigned K, class KmerEater, class Itr>
void KmerizationProcessor<K,KmerEater,Itr>::operator()( size_t batchNo )
{
    size_t readId = std::min(batchNo*mBatchSize,mNReads);
    size_t endId = std::min(readId+mBatchSize,mNReads);
    for ( Itr itr(mReads+readId), end(mReads+endId); itr != end; ++itr )
        KMer<K>::kmerizeIntoEater(itr->begin(),itr->end(),mEater,readId++);
    if ( mpDotter ) mpDotter->batchDone();
}

template <unsigned K, class KmerEater>
class VMVProcessor
{
    typedef VirtualMasterVec<bvec>::const_iterator Itr;
public:
    VMVProcessor( VirtualMasterVec<bvec> const& vmv,
                              KmerEater const& eater, Dotter* pDotter )
    : mVMV(vmv), mSubProcessor(mVMV.begin(),mVMV.end(),eater,pDotter) {}

    VMVProcessor( VMVProcessor const& that )
    : mVMV(that.mVMV), mSubProcessor(mVMV.begin(),mVMV.end(),
                                      that.mSubProcessor.getEater(),
                                      that.mSubProcessor.getDotter())
    {}

    // compiler-supplied copy by assignment and destructor is OK

    void operator()( size_t batchNo ) { mSubProcessor(batchNo); }

private:
    VirtualMasterVec<bvec> mVMV;
    KmerizationProcessor<K,KmerEater,Itr> mSubProcessor;
};

// output iterator that increments counts for the kmers in a dictionary
template <unsigned K>
class DictionaryKmerEater
{
public:
    DictionaryKmerEater( KmerDict<K>* pDict ) : mpDict(pDict) {}

    // compiler-supplied copying and destructor are OK

    void operator()( KMer<K> const& kmer, KMerContext context, size_t, size_t )
    { switch ( kmer.getCanonicalForm() )
      { case CanonicalForm::FWD:
          mpDict->applyCanonical(kmer,Functor(context));
          break;
        case CanonicalForm::REV:
          mpDict->applyCanonical(KMer<K>(kmer).rc(),Functor(context.rc()));
          break;
        case CanonicalForm::PALINDROME:
          mpDict->applyCanonical(kmer,Functor(context|=context.rc()));
          break; } }

private:
    class Functor
    {
    public:
        Functor( KMerContext context ) : mContext(context) {}

        // compiler-supplied copying and destructor are OK

        void operator()( KmerDictEntry<K> const& entry ) const
        { KDef& kDef = const_cast<KDef&>(entry.getKDef());
          kDef.incrementCount();
          kDef.orContext(mContext); }

    private:
        KMerContext mContext;
    };

    KmerDict<K>* mpDict;
};

template <unsigned K>
class Spectrum
{
public:
    template <class Itr>
    Spectrum( Itr beg, Itr end )
    : mCounts(end-beg), mTotCount(0)
    { std::copy(beg,end,mCounts.begin());
      init(); }

    Spectrum( KmerDict<K> const& dict, size_t maxCount = 100 )
    : mCounts(maxCount+1), mTotCount(0)
    { typedef typename KmerDict<K>::OCItr OCItr;
      typedef typename KmerDict<K>::ICItr ICItr;
      for ( OCItr oItr(dict.begin()), oEnd(dict.end()); oItr!=oEnd; ++oItr )
      { for ( ICItr itr(oItr->begin()), end(oItr->end()); itr != end; ++itr )
        { size_t count = itr->getKDef().getCount();
          mTotCount += count;
          if ( count > maxCount ) count = maxCount;
          mCounts[count] += 1; } }
      init(); }

    size_t size() const { return mCounts.size(); }
    size_t const& operator[]( size_t idx ) const { return mCounts[idx]; }

    size_t getFirstTroughIndex() const { return mTroughIdx; }
    size_t getFirstPeakIndex() const { return mPeakIdx; }

    bool isFirstPeakValid() const
    { double const ERROR_FRACTION = .9;
      return mPeakIdx > mTroughIdx &&
            mCounts[mPeakIdx] >= mCounts[mTroughIdx]; }

    bool isPeakDistinct() const
    { size_t idx = mPeakIdx - sqrt(mPeakIdx);
      return idx > mTroughIdx;  }

    size_t sum( size_t minIdx = 0, size_t maxIdx = ~0ul ) const
    { typedef std::vector<size_t>::const_iterator VItr;
      VItr itr = mCounts.begin();
      if ( maxIdx > mCounts.size() ) maxIdx = mCounts.size();
      return std::accumulate(itr+minIdx,itr+maxIdx,0ul); }

    size_t weightedSum( size_t minIdx = 1, size_t maxIdx = ~0ul ) const
    { if ( maxIdx > mCounts.size() ) maxIdx = mCounts.size();
      size_t result = 0;
      for ( size_t idx = minIdx; idx < maxIdx; ++idx )
          result += idx*mCounts[idx];
      return result; }

    double goodFract() const { return 1.*weightedSum(mTroughIdx)/mTotCount; }
    double errFract() const { return 1.*weightedSum(1,mTroughIdx)/mTotCount; }
    size_t getNKmers() const { return mTotCount; }

    friend ostream& operator<<( ostream& os, Spectrum const& sp )
    { size_t maxCount = sp.size()-1;
      for ( size_t idx = 1; idx < maxCount; ++idx )
        os << idx << "x\t" << sp[idx] << std::endl;
      os << ">=" << maxCount << "x\t" << sp[maxCount] << std::endl;
      return os; }

private:
    void init()
    { size_t min = mCounts[1];
      size_t max = min;
      mTroughIdx = 1;
      mPeakIdx = mTroughIdx;
      double const MIN_PEAK_TO_TROUGH_RATIO = 1.4;
      size_t const MIN_PEAK_COUNT = 100;
      size_t maxCount = mCounts.size()-1;
      for ( size_t idx = mTroughIdx+1; idx < maxCount; ++idx )
      { if ( mCounts[idx] < min )
        { if ( max >= MIN_PEAK_TO_TROUGH_RATIO*min ) break;
          min = max = mCounts[idx];
          mTroughIdx = idx; mPeakIdx = idx; }
        else if ( mCounts[idx] > max && mCounts[idx] >= MIN_PEAK_COUNT )
        { max = mCounts[idx]; mPeakIdx = idx; } } }

    std::vector<size_t> mCounts;
    size_t mTotCount;
    size_t mTroughIdx;
    size_t mPeakIdx;
};

// creates the files that represent a unipath graph
template <unsigned K>
class GraphBuilder
{
public:
    GraphBuilder( KmerDict<K>& dict )
    : mDict(dict), mStepper(dict), mNComponents(0)
    { mSeq.reserve(K*dict.size()/10); mEdges.reserve(dict.size()/100); }

    // compiler-supplied destructor is OK

    void process( const int verbosity = 0 );

private:
    GraphBuilder( GraphBuilder const& ); // unimplemented -- no copying
    GraphBuilder& operator=( GraphBuilder const& ); // unimplemented--no copying

    typedef typename KmerDict<K>::Entry DictEntry;
    struct EdgeStart
    {
        EdgeStart( DictEntry const& entry, bool isRC, bool isExtension )
        : mpEntry(&entry), mRC(isRC), mExtension(isExtension) {}

        DictEntry const* mpEntry;
        bool mRC;
        bool mExtension;
    };
    typedef std::vector<EdgeStart> EdgeStartVec;

    void buildEdges( EdgeStartVec& toBuild );
    void joinGraph();

    EdgeID fromRef( UnipathEdge const& unipath ) const
    { return EdgeID(&unipath - &mEdges[0]); }

    UnipathEdge& toRef( EdgeID const& edgeID )
    { AssertNot(edgeID.isNull()); return mEdges[edgeID.val()]; }

    void joinSuccessor( KMer<K> const& kmer, KDef* pDef, unsigned predCode,
                            EdgeID const& predID )
    { UnipathEdge& pred = toRef(predID);
      EdgeID const& succID = pDef->getEdgeID();
      UnipathEdge& succ = toRef(succID);
      KmerID succKmerID = succ.getKmerID(pDef->getEdgeOffset());
      unsigned succCode = kmer.back();
      using std::equal;
      if ( succKmerID == succ.getInitialKmerID() &&
           equal(kmer.begin(),kmer.end(),
                   mSeq.begin(succ.getInitialKmerID().val())) )
      { pred.setSuccessor(succCode,succID,false);
        succ.setPredecessor(predCode,predID,false); }
      else
      { AssertEq( succKmerID, succ.getFinalKmerID() );
        Assert( equal(kmer.rcbegin(),kmer.rcend(),
                        mSeq.begin(succ.getFinalKmerID().val())) );
        pred.setSuccessor(succCode,succID,true);
        succ.setSuccessor(GetComplementaryBase(predCode),predID,true); } }

    void joinPredecessor( unsigned succCode, EdgeID const& succID,
                            KMer<K> const& kmer, KDef* pDef )
    { UnipathEdge& succ = toRef(succID);
      EdgeID const& predID = pDef->getEdgeID();
      UnipathEdge& pred = toRef(predID);
      KmerID predKmerID = pred.getKmerID(pDef->getEdgeOffset());
      unsigned predCode = kmer.front();
      using std::equal;
      if ( predKmerID == pred.getFinalKmerID() &&
              equal(kmer.begin(),kmer.end(),
                      mSeq.begin(pred.getFinalKmerID().val())) )
      { pred.setSuccessor(succCode,succID,false);
        succ.setPredecessor(predCode,predID,false); }
      else
      { AssertEq( predKmerID, pred.getInitialKmerID() );
        Assert( equal(kmer.rcbegin(),kmer.rcend(),
                        mSeq.begin(pred.getInitialKmerID().val())) );
        pred.setPredecessor(GetComplementaryBase(succCode),succID,true);
        succ.setPredecessor(predCode,predID,true); } }

    // N50 edge length in kmers
    size_t getN50EdgeLen() const
    { vec<size_t> lengths;
      size_t nnn = mEdges.size();
      lengths.reserve(nnn);
      typedef UnipathEdgeVec::const_iterator Itr;
      for ( Itr itr(mEdges.begin()), end(mEdges.end()); itr != end; ++itr )
        lengths.push_back(itr->getLength());
      std::sort(lengths.begin(),lengths.end());
      return N50(lengths); }

    // median edge length in kmers
    size_t medianEdgeLen() const
    { std::vector<size_t> lengths;
      size_t nnn = mEdges.size();
      lengths.reserve(nnn);
      typedef UnipathEdgeVec::const_iterator Itr;
      for ( Itr itr(mEdges.begin()), end(mEdges.end()); itr != end; ++itr )
        lengths.push_back(itr->getLength());
      std::vector<size_t>::iterator med(lengths.begin()+nnn/2);
      std::nth_element(lengths.begin(),med,lengths.end());
      return *med; }

    KmerDict<K>& mDict;
    KmerStepper<K> mStepper;
    HugeBVec mSeq;
    UnipathEdgeVec mEdges;
    size_t mNComponents;
    friend class UnipathGraph<K>;
};

template <unsigned K>
void GraphBuilder<K>::process( const int verbosity )
{
    if ( verbosity >= 1 ) 
         std::cout << Date() << ": building unipaths" << std::endl;

    typedef typename KmerDict<K>::OCItr OCItr;
    typedef typename KmerDict<K>::ICItr ICItr;
    OCItr oEnd(mDict.cend());
    EdgeStartVec toBuild;

    DictEntry const* entries[4];
    for ( OCItr oItr(mDict.cbegin()); oItr != oEnd; ++oItr )
    {
        for ( ICItr itr(oItr->begin()), end(oItr->end()); itr != end; ++itr )
        {
            DictEntry const& entry = *itr;
            if ( entry.getKDef().isNull() )
            {
                if ( mStepper.getPredecessors(entry,entries) != 1 )
                {
                    LOG("IF: " << static_cast<KMer<K> const&>(entry));
                    toBuild.push_back(EdgeStart(entry,false,false));
                    buildEdges(toBuild);
                    mNComponents += 1;
                }
                else if ( mStepper.getSuccessors(entry,entries) != 1 )
                {
                    LOG("IR: " << static_cast<KMer<K> const&>(entry));
                    toBuild.push_back(EdgeStart(entry,true,false));
                    buildEdges(toBuild);
                    mNComponents += 1;
                }
            }
        }
    }

    // leaves only smooth rings, but we'd better scan for that.  you never know.
    for ( OCItr oItr(mDict.cbegin()); oItr != oEnd; ++oItr )
    {
        for ( ICItr itr(oItr->begin()), end(oItr->end()); itr != end; ++itr )
        {
            DictEntry const& entry = *itr;
            if ( entry.getKDef().isNull() )
            {
                LOG("I1: " << static_cast<KMer<K> const&>(entry));
                toBuild.push_back(EdgeStart(entry,false,false));
                buildEdges(toBuild);
                mNComponents += 1;
            }
        }
    }

    joinGraph();
    if ( verbosity >= 1 )
    {
    std::cout << Date( ) << ": there are " << mEdges.size() << " unipaths in "
              << mNComponents << " connected components\n" << Date( )
              << ": N50 unipath = " << getN50EdgeLen() << " kmers; "
              << "median unipath = " << medianEdgeLen() << " kmers" << std::endl;
    }
}

template <unsigned K>
void GraphBuilder<K>::buildEdges( EdgeStartVec& toBuild )
{
    KMer<K> kmer;

    while ( !toBuild.empty() )
    {
        EdgeStart const& es = toBuild.back();
        KDef* pDef = const_cast<KDef*>(&es.mpEntry->getKDef());
        if ( !pDef->isNull() )
        {
            toBuild.pop_back();
            continue;
        }

        kmer = *es.mpEntry;
        if ( es.mRC ) kmer.rc();
        if ( es.mExtension )
            mSeq.push_back(kmer.back());
        else
            mSeq.append(kmer.begin(),kmer.end());

        EdgeID edgeID(mEdges.size());
        bool isPalindrome = kmer.isPalindrome();
        mEdges.push_back(UnipathEdge(KmerID(mSeq.size()-kmer.size()),
                                     ComponentID(mNComponents),
                                     isPalindrome));
        pDef->set( edgeID, 0 );
        toBuild.pop_back();

        LOG("BF: " << kmer << ' ' << edgeID);

        DictEntry const* entries[4];

        // queue predecessors
        if ( mStepper.getPredecessors(kmer,entries) )
        {
            for ( unsigned predCode = 0; predCode < 4u; ++predCode )
            {
                DictEntry const* pEntry = entries[predCode];
                if ( pEntry && pEntry->getKDef().isNull() )
                {
                    KMer<K>& kmerPred = mStepper.getSteppedKmer();
                    kmerPred.setFront(predCode);
                    LOG("QR: " << KMer<K>(kmerPred).rc() << ' ' << edgeID);
                    bool isRC = kmerPred.getCanonicalForm()!=CanonicalForm::REV;
                    toBuild.push_back(EdgeStart(*pEntry,isRC,false));
                }
            }
        }

        // while we're on the straight and narrow, with just a single successor,
        // try to extend the edge
        unsigned nSuccessors;
        while ( (nSuccessors = mStepper.getSuccessors(kmer,entries)) == 1 )
        {
            // if we're a palindrome we must queue successors
            if ( isPalindrome )
                break;
            nSuccessors = 0; // we'll handle our single successor now

            DictEntry const* pEntry;
            for ( unsigned succCode = 0; succCode < 4u; ++succCode )
            {
                if ( (pEntry = entries[succCode]) )
                {
                    kmer = mStepper.getSteppedKmer();
                    kmer.setBack(succCode);
                    break;
                }
            }

            // if single successor is already built, we're done extending
            if ( !pEntry->getKDef().isNull() )
                break;

            // if successor has multiple predecessors, or is a palindrome
            // start a new edge
            CanonicalForm form = kmer.getCanonicalForm();
            if ( form==CanonicalForm::PALINDROME ||
                    mStepper.getPredecessors(kmer,entries) > 1 )
            {
                LOG("QF: " << kmer << ' ' << edgeID);
                toBuild.push_back(
                        EdgeStart(*pEntry,form==CanonicalForm::REV,true) );
                break;
            }

            // no branching -- extend the edge we're building
            LOG("XF: " << kmer << ' ' << edgeID);
            mSeq.push_back(kmer.back());
            pDef = const_cast<KDef*>(&pEntry->getKDef());
            pDef->set( edgeID, toRef(edgeID).extend() );
        }

        // if there are unhandled successors (multiple successors, or the single
        // successor of a palindrome), queue them up
        if ( nSuccessors )
        {
            for ( unsigned succCode = 0; succCode < 4u; ++succCode )
            {
                DictEntry const* pEntry = entries[succCode];
                if ( pEntry && pEntry->getKDef().isNull() )
                {
                    KMer<K>& kmerSucc = mStepper.getSteppedKmer();
                    kmerSucc.setBack(succCode);
                    LOG("QF: " << kmerSucc << ' ' << edgeID);
                    bool isRC = kmerSucc.getCanonicalForm()==CanonicalForm::REV;
                    toBuild.push_back(EdgeStart(*pEntry,isRC,!--nSuccessors));
                }
            }
        }
    }
}

template <unsigned K>
void GraphBuilder<K>::joinGraph()
{
    typedef UnipathEdgeVec::iterator Itr;
    DictEntry const* entries[4];
    KMer<K> kmer;
    for ( size_t idx = 0; idx < mEdges.size(); ++idx )
    {
        UnipathEdge& edge = mEdges[idx];
        kmer.assign(mSeq.begin(edge.getInitialKmerID().val()));
        unsigned succCode = kmer.back();
        unsigned predCode;
        if ( mStepper.getPredecessors(kmer,entries) )
        {
            EdgeID succID = fromRef(edge);
            for ( predCode = 0; predCode < 4u; ++predCode )
            {
                DictEntry const* pEntry = entries[predCode];
                if ( pEntry && edge.getPredecessor(predCode).isNull() )
                {
                    KDef* pDef = const_cast<KDef*>(&pEntry->getKDef());
                    KMer<K>& kmerPred = mStepper.getSteppedKmer();
                    kmerPred.setFront(predCode);
                    LOG("JP: " << kmerPred << ' ' << succID <<
                            " onto " << pDef->getEdgeID());
                    joinPredecessor(succCode,succID,kmerPred,pDef);
                }
            }
        }
        kmer.assign(mSeq.begin(edge.getFinalKmerID().val()));
        predCode = kmer.front();
        if ( mStepper.getSuccessors(kmer,entries) )
        {
            EdgeID predID = fromRef(edge);
            for ( succCode = 0; succCode < 4u; ++succCode )
            {
                DictEntry const* pEntry = entries[succCode];
                if ( pEntry && edge.getSuccessor(succCode).isNull() )
                {
                    KDef* pDef = const_cast<KDef*>(&pEntry->getKDef());
                    KMer<K>& kmerSucc = mStepper.getSteppedKmer();
                    kmerSucc.setBack(succCode);
                    LOG("JS: " << kmerSucc << ' ' << predID <<
                            " onto " << pDef->getEdgeID());
                    joinSuccessor(kmerSucc,pDef,predCode,predID);
                }
            }
        }
    }
}

template <unsigned K> class PathingProcessor;

// Paths sequence onto UnipathGraph.
template <unsigned K>
class PathBuilder
{
public:
    PathBuilder( String const& infoFile, UnipathGraph<K> const& graph );

    void processReads( VirtualMasterVec<bvec> const& reads,
                            KmerDict<K> const& dict,
                            unsigned nThreads,
                            size_t batchSize = 10000 );

    void writeEvidence( unsigned nThreads );

    // compiler-supplied destructor is OK

private:
    PathBuilder( PathBuilder const& ); // unimplemented -- no copying
    PathBuilder& operator=( PathBuilder const& ); // unimplemented -- no copying

    // the rest of these methods are for use by the friendly PathingProcessor
    friend class PathingProcessor<K>;

    bool isRC( KMer<K> const& kmer, KmerID const& kmerID ) const
    { using std::equal;
      HugeBVec::const_iterator seqItr(mGraph.getBases(kmerID));
      bool result = !equal(kmer.begin(),kmer.end(),seqItr);
      Assert(!result || equal(kmer.rcbegin(),kmer.rcend(),seqItr));
      return result; }

    UnipathEdge const& getEdge( EdgeID const& edgeID ) const
    { return mGraph.getEdge(edgeID); }

    bool nextEdge( unsigned baseCode, EdgeID* pEdgeID, bool isRC ) const
    { return mGraph.nextEdge(baseCode,pEdgeID,isRC); }

    PathID addPathSeq( bvec const& bases )
    { SpinLocker locker(mPathSeqLock);
      PathID pathID(mPathSeq.size());
      mPathSeq.append(bases.begin(),bases.end());
      return pathID; }

    String mInfoFile;
    UnipathGraph<K> const& mGraph;
    HugeBVec::Builder mPathSeq; // indexed by PathID
    SpinLockedData mPathSeqLock;
};

// worklist processor class to kmerize reads
template <unsigned K>
class PathingProcessor
{
public:
    PathingProcessor( VirtualMasterVec<bvec> const& vmv,
                      KmerDict<K> const& dict,
                      PathBuilder<K>& pathBuilder,
                      ChunkDumper<UnipathingVec>* pPathingsDumper,
                      Dotter* pDotter )
    : mVMV(vmv), mDict(dict), mPathBuilder(pathBuilder),
      mpPathingsDumper(pPathingsDumper), mpDotter(pDotter)
    {}

    PathingProcessor( PathingProcessor const& that )
    : mVMV(that.mVMV), mDict(that.mDict), mPathBuilder(that.mPathBuilder),
      mpPathingsDumper(that.mpPathingsDumper), mpDotter(that.mpDotter)
    { mReads.reserve(getBatchSize()); }

    // copy assignment prohibited by references, compiler-supplied destructor OK

    size_t getBatchSize() const
    { size_t result = mpDotter->getNBatches();
      if ( result ) result = (mVMV.size()+result-1)/result;
      return result; }

    void operator()( size_t batchNo )
    { size_t batchSize = getBatchSize();
      size_t readId = std::min(batchNo*batchSize,mVMV.size());
      size_t endId = std::min(readId+batchSize,mVMV.size());
      size_t nnn = endId - readId;
      mPathings.reserve(nnn);
      typedef VirtualMasterVec<bvec>::const_iterator Itr;
      Itr end(mVMV.begin(endId));
      for ( Itr itr(mVMV.begin(readId)); itr != end; ++itr, ++readId )
          path(readId,*itr);
      mpPathingsDumper->dumpChunk(batchNo,mPathings);
      mpDotter->batchDone(); }

private:
    void path( size_t readId, bvec const& read );

    VirtualMasterVec<bvec> mVMV;
    KmerDict<K> const& mDict;
    PathBuilder<K>& mPathBuilder;
    ChunkDumper<UnipathingVec>* mpPathingsDumper;
    Dotter* mpDotter;
    vecbvec mReads;
    bvec mPathSeq;
    UnipathingVec mPathings;

    static SpinLockedData gLock;
};

template <unsigned K> SpinLockedData PathingProcessor<K>::gLock;

template <unsigned K>
void PathingProcessor<K>::path( size_t readId, bvec const& read )
{
    size_t readLen = read.size();
    if ( readLen < K )
    {
        mPathings.push_back(Unipathing());
        return; // LOOKOUT:  early return
    }

    KMer<K> kmer(read.begin());
    KDef const* pDef = mDict.lookup(kmer);
    if ( !pDef )
        FatalErr("Can't find kmer " << kmer << " in dictionary.");

    UnipathEdge const* pEdge = &mPathBuilder.getEdge(pDef->getEdgeID());
    KmerID kmerID = pEdge->getKmerID(pDef->getEdgeOffset());
    bool isRC = mPathBuilder.isRC(kmer,kmerID);
    size_t skipped = isRC ? pEdge->getFinalKmerID().val() - kmerID.val() :
                            kmerID.val() - pEdge->getInitialKmerID().val();
    Unipathing up(pDef->getEdgeID(),isRC,skipped);
    size_t readIdx = K + pEdge->getLength() - skipped - 1;

    if ( readIdx < readLen )
    {
        mPathSeq.clear();
        EdgeID edgeID = pDef->getEdgeID();
        do
        {
            unsigned char baseCode = read[readIdx];
            if ( isRC ) baseCode = GetComplementaryBase(baseCode);
            mPathSeq.push_back(baseCode);
            isRC = mPathBuilder.nextEdge(baseCode,&edgeID,isRC);
            if ( edgeID.isNull() )
                FatalErr("Can't path read " << readId << " " << read.ToString()
                         << " at " << readIdx );
            pEdge = &mPathBuilder.getEdge(edgeID);
            readIdx += pEdge->getLength();
        }
        while ( readIdx < readLen );

        unsigned nBases = mPathSeq.size();
        up.setSegments( nBases+1,
                        nBases <= 8 ?
                           PathID(8*mPathSeq.extractKmer(0,nBases)) :
                           mPathBuilder.addPathSeq(mPathSeq) );
    }

    up.setFinalSkip(readIdx-readLen);
    mPathings.push_back(up);
}

template <unsigned K>
PathBuilder<K>::PathBuilder( String const& infoFile,
                                UnipathGraph<K> const& graph )
: mInfoFile(infoFile), mGraph(graph),
  mPathSeq(PathCollection<K>::getPathseqFilename(infoFile).c_str())
{
    unlink(infoFile.c_str());

    // seed traversal sequence with all 8-base sequences
    bvec bv;
    for ( unsigned idx = 0; idx < 65536; ++idx )
    {
        bv.assignBaseBits(8,&idx);
        mPathSeq.append(bv.begin(),bv.end());
    }
}

template <unsigned K>
void PathBuilder<K>::processReads( VirtualMasterVec<bvec> const& reads,
                                            KmerDict<K> const& dict,
                                            unsigned nThreads,
                                            size_t batchSize )
{
    size_t nReads = reads.size();
    size_t nBatches = (nReads+batchSize-1)/batchSize;
    Dotter dotter(nBatches);
    String pathsFilename = PathCollection<K>::getPathsFilename(mInfoFile);
    ChunkDumper<UnipathingVec> pathingsDumper(pathsFilename.c_str(),
                                              nReads,nBatches);
    PathingProcessor<K> proc(reads,dict,*this,&pathingsDumper,&dotter);
    std::cout << Date() << ": Pathing " << nReads << " reads in " << nBatches
                << " batches of " << proc.getBatchSize() << '.' << std::endl;

    if ( !nThreads )
        nThreads = getConfiguredNumThreads();
    parallelFor(0ul,nBatches,proc,nThreads);
}

template <unsigned K>
class EvidenceProcessor
{
public:
    struct Common : private Dotter
    {
        Common( UnipathGraph<K> const& graph,
                      HugeBVec const& pathSeq,
                      UnipathingVec const& pathings,
                      ChunkDumper<VecUnipathEvidenceVec>& dumper,
                      size_t nBatches, size_t batchSize,
                      std::vector<unsigned> const& evCounts )
        : Dotter(nBatches), mGraph(graph), mPathSeq(pathSeq), mPathings(pathings),
          mDumper(dumper), mBatchSize(batchSize), mEvCounts(evCounts)
        {}

        void dumpChunk( size_t chunk, VecUnipathEvidenceVec& vuev )
        { mDumper.dumpChunk(chunk,vuev); batchDone(); }

        UnipathGraph<K> const& mGraph;
        HugeBVec const& mPathSeq;
        UnipathingVec const& mPathings;
        ChunkDumper<VecUnipathEvidenceVec>& mDumper;
        size_t mBatchSize;
        std::vector<unsigned> const& mEvCounts;
    };

    EvidenceProcessor( Common& common )
    : mCommon(common)
    {}

    EvidenceProcessor( EvidenceProcessor const& that )
    : mCommon(that.mCommon)
    {}

    // copy assignment prohibited by references, compiler-supplied destructor OK

    void operator()( size_t chunk );

private:
    Common& mCommon;
    VecUnipathEvidenceVec mVUEV;
};

template <unsigned K>
void EvidenceProcessor<K>::operator()( size_t chunk )
{
    typedef VecUnipathEvidenceVec::iterator Itr;
    size_t minEdgeId = chunk*mCommon.mBatchSize;
    size_t maxEdgeId = std::min(minEdgeId+mCommon.mBatchSize,
                                mCommon.mGraph.getNEdges());
    mVUEV.clear().resize(maxEdgeId-minEdgeId);
    size_t idx = minEdgeId;
    for ( Itr itr(mVUEV.begin()), end(mVUEV.end()); itr != end; ++itr )
        itr->reserve(mCommon.mEvCounts[idx++]);

    size_t nReads = mCommon.mPathings.size();
    for ( size_t readId = 0; readId < nReads; ++readId )
    {
        ReadID rid(readId);
        Unipathing const& pathing = mCommon.mPathings[readId];
        EdgeID edgeID = pathing.getInitialEdgeID();
        bool isRC = pathing.isInitialEdgeRC();
        size_t edgeId = edgeID.val();
        if ( minEdgeId <= edgeId && edgeId < maxEdgeId )
        {
            UnipathEvidenceVec& ev = mVUEV[edgeId-minEdgeId];
            ev.push_back(UnipathEvidence(rid,0u,isRC));
        }

        unsigned nSegs = pathing.getNSegments();
        if ( nSegs > 1 )
        {
            size_t off = pathing.getPathID().val();
            HugeBVec::const_iterator itr = mCommon.mPathSeq.begin(off-1);
            for ( unsigned segNo = 1; segNo < nSegs; ++segNo )
            {
                isRC = mCommon.mGraph.nextEdge(*++itr,&edgeID,isRC);
                if ( edgeID.isNull() )
                    FatalErr("Invalid successor base in path sequence for read "
                             << readId );
                edgeId = edgeID.val();
                if ( minEdgeId <= edgeId && edgeId < maxEdgeId )
                {
                    UnipathEvidenceVec& ev = mVUEV[edgeId-minEdgeId];
                    ev.push_back(UnipathEvidence(rid,segNo,isRC));
                }
            }
        }
    }
    mCommon.dumpChunk(chunk,mVUEV);
}

template <unsigned K>
void PathBuilder<K>::writeEvidence( unsigned nThreads )
{
    std::cout << Date() << " Creating edge location dictionary." << std::endl;

    mGraph.ungetDict();

    UnipathingVec pathings;
    String pathsFilename = PathCollection<K>::getPathsFilename(mInfoFile);
    BinaryReader::readFile(pathsFilename.c_str(),&pathings);

    mPathSeq.close();
    String pathSeqFilename = PathCollection<K>::getPathseqFilename(mInfoFile);
    HugeBVec pathSeq(pathSeqFilename.c_str());

    size_t nEdges = mGraph.getNEdges();
    std::vector<unsigned> evCounts;
    evCounts.resize(nEdges);

    size_t nReads = pathings.size();
    size_t evTot = 0;
    for ( size_t readId = 0; readId < nReads; ++readId )
    {
        Unipathing const& pathing = pathings[readId];
        if ( pathing.isNull() )
            continue; // Non-structured code!

        EdgeID edgeID = pathing.getInitialEdgeID();
        evCounts[edgeID.val()] += 1;
        evTot += 1;
        unsigned nSegs = pathing.getNSegments();
        if ( nSegs > 1 )
        {
            bool isRC = pathing.isInitialEdgeRC();
            size_t off = pathing.getPathID().val();
            HugeBVec::const_iterator itr = pathSeq.begin(off-1);
            for ( unsigned segNo = 1; segNo < nSegs; ++segNo )
            {
                isRC = mGraph.nextEdge(*++itr,&edgeID,isRC);
                evCounts[edgeID.val()] += 1;
                evTot += 1;
            }
        }
    }

    if ( true )
    {
        size_t bytesPer = 2*((evTot*sizeof(UnipathEvidence)+nEdges-1)/nEdges +
                                    sizeof(UnipathEvidenceVec));
        WorklistParameterizer wp(nEdges,bytesPer,100,nThreads);
        size_t nBatches = wp.getNBatches();
        size_t batchSize = wp.getBatchSize();
        String evFilename = PathCollection<K>::getEvidenceFilename(mInfoFile);
        ChunkDumper<VecUnipathEvidenceVec> dumper(evFilename.c_str(),
                                                  nEdges,nBatches);
        typename EvidenceProcessor<K>::Common common(mGraph, pathSeq, pathings,
                                                     dumper, nBatches,
                                                     batchSize, evCounts);
        EvidenceProcessor<K> proc(common);
        std::cout << "Recording evidence for " << nEdges << " unipaths in "
                    << nBatches << " batches of "
                    << batchSize << '.' << std::endl;
        parallelFor(0ul,nBatches,proc,wp.getNThreads());
    }

    PathInfo pi(pathings.size(),pathSeq.size());
    BinaryWriter::writeFile(mInfoFile, pi);
}

} // end of anonymous namespace

template <unsigned K>
void KmerDict<K>::process( VirtualMasterVec<bvec> const& reads,
                            bool validate, unsigned nThreads, size_t batchSize )
{
    size_t nReads = reads.size();
    size_t nBatches = (nReads+batchSize-1)/batchSize;
    Dotter dotter(nBatches);
    typedef DictionaryKmerEater<K> Eater;
    typedef VMVProcessor<K,Eater> KProc;
    KProc proc(reads,Eater(this),&dotter);
    if ( !nThreads )
        nThreads = getConfiguredNumThreads();
    if ( nThreads > nBatches )
        nThreads = nBatches;
    std::cout << Date() << ": Kmerizing " << nReads << " reads in " << nBatches
                << " batches of " << batchSize << '.' << std::endl;
    parallelFor(0ul,nBatches,proc,nThreads);
    if ( validate )
    {
        size_t nErrs = mKSet.validateBinAssignments();
        if ( nErrs )
            FatalErr("There were " << nErrs
                        << " misplaced dictionary assignments.");
    }
}

template <unsigned K>
void KmerDict<K>::process( vecbvec const& reads, bool verbose,
                            bool validate, unsigned nThreads, size_t batchSize )
{
    size_t nReads = reads.size();
    size_t nBatches = (nReads+batchSize-1)/batchSize;
    Dotter dotter(nBatches);
    typedef DictionaryKmerEater<K> Eater;
    typedef vecbvec::const_iterator Itr;
    typedef KmerizationProcessor<K,Eater,Itr> KProc;
    KProc proc(reads.begin(),reads.end(),Eater(this),verbose?&dotter:0);
    if ( !nThreads )
        nThreads = getConfiguredNumThreads();
    if ( nThreads > nBatches )
        nThreads = nBatches;
    if ( verbose )
        std::cout << Date() << ": Kmerizing " << nReads << " reads in "
                    << nBatches << " batches of " << batchSize << '.'
                    << std::endl;
    parallelFor(0ul,nBatches,proc,nThreads);
    if ( validate )
    {
        size_t nErrs = mKSet.validateBinAssignments();
        if ( nErrs )
            FatalErr("There were " << nErrs
                        << " misplaced dictionary assignments.");
    }
}

template <unsigned K>
UnipathGraph<K>::UnipathGraph( VirtualMasterVec<bvec> const& reads,
                                bool validate, unsigned nThreads,
                                size_t nKmersEstimate, const int verbosity )
: mpDict(new KmerDict<K>(nKmersEstimate)), mMyDict(true)
{
    mpDict->process(reads,validate,nThreads);
    GraphBuilder<K> builder(*mpDict);
    builder.process(verbosity);
    mSeq.swap(builder.mSeq);
    mEdges.swap(builder.mEdges);
    if ( validate )
        validateUnipaths(mSeq,mEdges,*mpDict);
}

template <unsigned K>
UnipathGraph<K>::UnipathGraph( vecbvec const& reads,
                                bool validate, unsigned nThreads,
                                size_t nKmersEstimate, const int verbosity )
: mpDict(new KmerDict<K>(nKmersEstimate)), mMyDict(true)
{
    mpDict->process(reads,verbosity,validate,nThreads);
    GraphBuilder<K> builder(*mpDict);
    builder.process(verbosity);
    mSeq.swap(builder.mSeq);
    mEdges.swap(builder.mEdges);
    if ( validate )
        validateUnipaths(mSeq,mEdges,*mpDict);
}

template <unsigned K>
UnipathGraph<K>::UnipathGraph( KmerDict<K>& dict, const int verbosity )
: mpDict(&dict), mMyDict(false)
{
    GraphBuilder<K> builder(dict);
    builder.process(verbosity);
    mSeq.swap(builder.mSeq);
    mEdges.swap(builder.mEdges);
}

template <unsigned K>
void PathCollection<K>::create( String const& fastb, bool validate,
                             unsigned nThreads, size_t nKmersEst )
{
    VirtualMasterVec<bvec> reads(fastb);

    if ( !nKmersEst )
        nKmersEst = 4*reads.sizeSum()/50; // assume 50x coverage
    UnipathGraph<K> graph(reads,validate,nThreads,nKmersEst);
    graph.write(UnipathGraph<K>::getInfoFilename(fastb));

    PathBuilder<K> pb(getInfoFilename(fastb),graph);
    pb.processReads(reads,graph.getDict(),nThreads);
    pb.writeEvidence(nThreads);
}

#endif // KMERS_READ_PATHER_DEFS_H_
