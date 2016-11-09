///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * BuildReadQGraph40.cc
 *
 *  Created on: Jan 22, 2014
 *      Author: tsharpe
 */

#include "paths/long/BuildReadQGraph40.h"
#include "Basevector.h"
#include "FastaFileset.h"
#include "Intvector.h"
#include "IteratorRange.h"
#include "MapReduceEngine.h"
#include "PairsManager.h"
#include "ParallelVecUtilities.h"
#include "Qualvector.h"
#include "Vec.h"
#include "dna/Bases.h"
#include "feudal/BinaryStream.h"
#include "feudal/VirtualMasterVec.h"
//#include "kmers/BigKPather.h"
#include "kmers/ReadPatherDefs.h"
#include "math/Functions.h"
#include "paths/KmerBaseBroker.h"
#include "paths/UnibaseUtils.h"
#include "paths/long/HBVFromEdges.h"
#include "paths/long/KmerCount.h"
#include "system/SortInPlace.h"
#include "system/SpinLockedData.h"
#include "system/WorklistN.h"
#include <algorithm>
#include <atomic>
#include <fstream>
#include <iostream>
#include <numeric>
#include <utility>
#include <vector>
#include "paths/HyperBasevector.h"
#include "paths/long/ExtendReadPath.h"
#include "paths/long/ShortKmerReadPather.h"

namespace
{
unsigned const K = 40;
typedef KMer<K> Kmer;
typedef KMer<K-1> SubKmer;
typedef KmerDictEntry<K> Entry;
typedef KmerDict<K> Dict;
typedef UnipathGraph<K> Graph;
typedef std::vector<std::atomic_size_t> SpectrumBins;

class GoodLenTailFinder
{
public:
    GoodLenTailFinder( VecPQVec const& quals, unsigned minQual,
                        std::vector<unsigned>* pGoodLens )
    : mQuals(quals), mMinQual(minQual), mGoodLens(*pGoodLens) {}

    void operator()( size_t readId )
    { mQuals[readId].unpack(&mQV);
      unsigned minQual = mMinQual;
      auto itr = mQV.end();
      auto beg = mQV.begin();
      unsigned good = 0;
      while ( itr != beg )
        if ( *--itr < minQual ) good = 0;
        else if ( ++good == K )
        { mGoodLens[readId] = (itr-beg)+K; return; }
      mGoodLens[readId] = 0; }

private:
    VecPQVec const& mQuals;
    unsigned mMinQual;
    std::vector<unsigned>& mGoodLens;
    qvec mQV;
};

inline void summarizeEntries( Entry* e1, Entry* e2 )
{
    KMerContext kc;
    size_t count = 0;
    while ( e2-- != e1 )
    {
        KDef const& kDef = e2->getKDef();
        kc |= kDef.getContext();
        count += std::max(kDef.getCount(),1ul);
    }
    KDef& kDef = e1->getKDef();
    kDef.setContext(kc);
    kDef.setCount(count);
}

class Kmerizer
{
public:
    Kmerizer( vecbvec const& reads, std::vector<unsigned> const& goodLengths,
                unsigned minFreq, Dict* pDict, std::atomic_size_t* pTotKmers )
    : mReads(reads), mGoodLengths(goodLengths), mMinFreq(minFreq),
      mpDict(pDict), mpTotKmers(pTotKmers), mNKmers(0)
    {}

    ~Kmerizer()
    { if ( mpTotKmers ) *mpTotKmers += mNKmers; }

    template <class OItr>
    void map( size_t readId, OItr oItr )
    { unsigned len = mGoodLengths[readId];
      if ( len < K+1 ) return;
      auto beg = mReads[readId].begin(), itr=beg+K, last=beg+(len-1);
      Kmer kkk(beg);
      KMerContext kc = KMerContext::initialContext(*itr);
      *oItr++ = kkk.isRev() ? Entry(Kmer(kkk).rc(),kc.rc()) : Entry(kkk,kc);
      while ( itr != last )
      { unsigned char pred = kkk.front();
        kkk.toSuccessor(*itr); ++itr;
        kc = KMerContext(pred,*itr);
        *oItr++ = kkk.isRev() ? Entry(Kmer(kkk).rc(),kc.rc()) : Entry(kkk,kc); }
      kc = KMerContext::finalContext(kkk.front());
      kkk.toSuccessor(*last);
      *oItr++ = kkk.isRev() ? Entry(Kmer(kkk).rc(),kc.rc()) : Entry(kkk,kc); }

    void reduce( Entry* e1, Entry* e2 )
    { summarizeEntries(e1,e2);
      if ( e1->getKDef().getCount() >= mMinFreq )
      { ++mNKmers;
        if ( mpDict ) mpDict->insertEntry(std::move(*e1)); } }

    Entry* overflow( Entry* e1, Entry* e2 )
    { if ( e2-e1 > 1 ) summarizeEntries(e1,e2); return e1+1; }

private:
    vecbvec const& mReads;
    std::vector<unsigned> const& mGoodLengths;
    unsigned mMinFreq;
    Dict* mpDict;
    std::atomic_size_t* mpTotKmers;
    size_t mNKmers;
};
typedef MapReduceEngine<Kmerizer,Entry,Kmer::Hasher> KMRE;

Dict* createDict( vecbvec const& reads, ObjectManager<VecPQVec>& quals,
                        unsigned minQual, unsigned minFreq, float const mem_frac = 0.9 )
{
    // figure out how much of the read to kmerize by examining quals
    //std::cout << Date() << ": processing quals." << std::endl;
    std::vector<unsigned> goodLens(reads.size());
    parallelForBatch(0ul,reads.size(),100000,
                         GoodLenTailFinder(quals.load(),minQual,&goodLens));
    quals.unload();
    size_t nKmers = std::accumulate(goodLens.begin(),goodLens.end(),0ul);

    if ( nKmers == 0 )
    {    cout << "\nLooks like your input data have almost no good bases.\n"
              << "Giving up.\n" << endl;
         Scram(1);    }

    size_t dictSize;
    if ( true )
    { // count uniq kmers that occur at minFreq or more
        std::atomic_size_t nUniqKmers(0);
        Kmerizer impl(reads,goodLens,minFreq,nullptr,&nUniqKmers);
        KMRE mre(impl);
        if ( !mre.run(nKmers,0ul,reads.size(),
				KMRE::VERBOSITY::QUIET, mem_frac) )
            FatalErr("Failed to count unique kmers.  Out of buffer space.  ");
        dictSize = nUniqKmers;
    }

    // kmerize reads into dictionary
    Dict* pDict = new Dict(dictSize);
    Kmerizer impl(reads,goodLens,minFreq,pDict,nullptr);
    KMRE mre(impl);
    if ( !mre.run(nKmers,0ul,reads.size(),
			    KMRE::VERBOSITY::QUIET, mem_frac) )
        FatalErr("Failed to kmerize.  Out of buffer space.  ");

    // some kmers were discarded because they didn't occur often enough to
    // convince us that they were real -- recompute adjacencies to compensate
    // for these missing kmers.
    if ( minFreq > 1 )
        pDict->recomputeAdjacencies();

    return pDict;
}

class EdgeBuilder
{
public:
    EdgeBuilder( Dict const& dict, vecbvec* pEdges )
    : mDict(dict), mEdges(*pEdges) {}

    void buildEdge( Entry const& entry )
    { if ( isPalindrome(entry) )
        make1KmerEdge(entry);
      else if ( upstreamExtensionPossible(entry) )
      { if ( downstreamExtensionPossible(entry) )
          return;
        extendUpstream(entry); }
      else if ( downstreamExtensionPossible(entry) )
        extendDownstream(entry);
      else
        make1KmerEdge(entry); }

    // not thread-safe
    void simpleCircle( Entry const& entry )
    { Entry const* pFirstEntry = &entry;
      mEdgeSeq.assign(entry.begin(),entry.end());
      mEdgeEntries.push_back(pFirstEntry);
      KMerContext context = entry.getKDef().getContext();
      Kmer kmer(entry);
      while ( true )
      { ForceAssertEq(context.getPredecessorCount(),1u);
        ForceAssertEq(context.getSuccessorCount(),1u);
        unsigned char succCode = context.getSingleSuccessor();
        kmer.toSuccessor(succCode);
        Entry const* pEntry = lookup(kmer,&context);
        if ( pEntry == pFirstEntry )
          break;
        if ( !pEntry->getKDef().isNull() )
        { std::cout << "Failed to close circle.\n";
          for ( auto beg=mEdgeEntries.begin(),end=mEdgeEntries.end();
                  beg != end; ++beg )
            std::cout << *beg << ' ' << **beg << '\n';
          std::cout << pEntry << ' ' << *pEntry << std::endl;
          CRD::exit(1); }
        mEdgeSeq.push_back(succCode);
        mEdgeEntries.push_back(pEntry); }
        canonicalizeCircle();
        addEdge(); }

private:
    void canonicalizeCircle()
    { auto itr = std::min_element(mEdgeEntries.begin(),mEdgeEntries.end(),
                            []( Kmer const* pEnt1, Kmer const* pEnt2 )
                            { return *pEnt1 < *pEnt2; });
      Kmer const& minKmer = **itr;
      size_t idx = itr - mEdgeEntries.begin();
      if ( CF<K>::getForm(mEdgeSeq.begin(idx))==CanonicalForm::REV )
      { mEdgeSeq.ReverseComplement();
        std::reverse(mEdgeEntries.begin(),mEdgeEntries.end());
        idx = mEdgeSeq.size()-idx-K;
        Assert(std::equal(minKmer.begin(),minKmer.end(),mEdgeSeq.begin(idx))); }
      if ( !idx )
          return;
      bvec bv; bv.reserve(mEdgeSeq.size());
      bv.assign(mEdgeSeq.begin(idx),mEdgeSeq.end());
      bv.append(mEdgeSeq.begin(K-1),mEdgeSeq.begin(K+idx-1));
      Assert(std::equal(mEdgeSeq.begin(),mEdgeSeq.begin(K-1),
                          mEdgeSeq.end()-(K-1)));
      mEdgeSeq = bv;
      mEdgeEntries.reserve(mEdgeEntries.size()+idx);
      std::copy(mEdgeEntries.begin(),mEdgeEntries.begin()+idx,
                  std::back_inserter(mEdgeEntries));
      mEdgeEntries.erase(mEdgeEntries.begin(),mEdgeEntries.begin()+idx); }

    bool isPalindrome( Kmer const& kmer )
    { if ( !(K&1) )
        return kmer.isPalindrome();
      SubKmer subKmer(kmer);
      if ( subKmer.isPalindrome() )
        return true;
      subKmer.toSuccessor(kmer.back());
      return subKmer.isPalindrome(); }

    bool upstreamExtensionPossible( Entry const& entry )
    { KMerContext context = entry.getKDef().getContext();
      if ( context.getPredecessorCount() != 1 )
        return false;
      Kmer pred(entry);
      pred.toPredecessor(context.getSinglePredecessor());
      if ( isPalindrome(pred) )
        return false;
      lookup(pred,&context);
      return context.getSuccessorCount() == 1; }

    bool downstreamExtensionPossible( Entry const& entry )
    { KMerContext context = entry.getKDef().getContext();
      if ( context.getSuccessorCount() != 1 )
        return false;
      Kmer succ(entry);
      succ.toSuccessor(context.getSingleSuccessor());
      if ( isPalindrome(succ) )
        return false;
      lookup(succ,&context);
      return context.getPredecessorCount() == 1; }

    void make1KmerEdge( Entry const& entry )
    { mEdgeSeq.assign(entry.begin(),entry.end());
      mEdgeEntries.push_back(&entry);
      addEdge(); }

    void extendUpstream( Entry const& entry )
    { mEdgeSeq.assign(entry.rcbegin(),entry.rcend());
      mEdgeEntries.push_back(&entry);
      extend(Kmer(entry).rc(),entry.getKDef().getContext().rc()); }

    void extendDownstream( Entry const& entry )
    { mEdgeSeq.assign(entry.begin(),entry.end());
       mEdgeEntries.push_back(&entry);
       extend(entry,entry.getKDef().getContext()); }

    void extend( Kmer const& kmer, KMerContext context )
    { Kmer next(kmer);
      while ( context.getSuccessorCount() == 1 )
      { unsigned char succCode = context.getSingleSuccessor();
        next.toSuccessor(succCode);
        if ( isPalindrome(next) )
          break;
        Entry const* pEntry = lookup(next,&context);
        if ( context.getPredecessorCount() != 1 )
          break;
        mEdgeSeq.push_back(succCode);
        mEdgeEntries.push_back(pEntry); }
      switch (mEdgeSeq.getCanonicalForm())
      {case CanonicalForm::PALINDROME:
         ForceAssertEq(mEdgeSeq.size(),K);
        // allow flow-through to FWD case
       case CanonicalForm::FWD:
         addEdge(); break;
       case CanonicalForm::REV:
         mEdgeSeq.clear(); mEdgeEntries.clear(); break; } }

    Entry const* lookup( Kmer const& kmer, KMerContext* pContext )
    { Entry const* result;
      if ( kmer.isRev() )
      { result = mDict.findEntryCanonical(Kmer(kmer).rc());
        ForceAssert(result);
        *pContext = result->getKDef().getContext().rc(); }
      else
      { result = mDict.findEntryCanonical(kmer);
        ForceAssert(result);
        *pContext = result->getKDef().getContext(); }
      return result; }

    void addEdge()
    { static SpinLockedData gLock;
      EdgeID edgeID;
      if ( mEdgeSeq.getCanonicalForm() == CanonicalForm::REV )
      {
          mEdgeSeq.ReverseComplement();
          std::reverse(mEdgeEntries.begin(),mEdgeEntries.end());
      }
      if ( true )
      { SpinLocker lock(gLock);
        edgeID.setVal(mEdges.size());
        mEdges.push_back(mEdgeSeq); }
      unsigned offset = 0;
      bool err = false;
      for ( Entry const* pEnt : mEdgeEntries )
      { KDef const& kDef = pEnt->getKDef();
        if ( kDef.isNull() )
        { Assert(std::equal(pEnt->begin(),pEnt->end(),mEdgeSeq.begin(offset)) ||
                 std::equal(pEnt->rcbegin(),pEnt->rcend(),mEdgeSeq.begin(offset)));
          const_cast<KDef&>(kDef).set(edgeID,offset++); }
        else
        { std::cout << edgeID << ':' << offset++ << ' ' << pEnt
                    << " Already occupied as " << kDef.getEdgeID()
                    << ':' << kDef.getEdgeOffset() << std::endl;
          err = true; } }
        if ( err )
            FatalErr("Having trouble with preoccupied kmers.");
        mEdgeSeq.clear();
        mEdgeEntries.clear(); }

    Dict const& mDict;
    vecbvec& mEdges;
    std::vector<Entry const*> mEdgeEntries;
    bvec mEdgeSeq;
};

void buildEdges( Dict const& dict, vecbvec* pEdges )
{
    EdgeBuilder eb(dict,pEdges);
    dict.parallelForEachHHS(
            [eb]( Dict::Set::HHS const& hhs ) mutable
            { for ( Entry const& entry : hhs )
                if ( entry.getKDef().isNull() )
                  eb.buildEdge(entry); });

    size_t nRegularEdges = pEdges->size();
    size_t nTotalLength = pEdges->SizeSum();
    //std::cout << "There are " << nRegularEdges
    //          << " edges of total length " << nTotalLength << '.' << std::endl;

    // if a kmer isn't marked as being on an edge, it must be a part of a smooth
    // circle: add those edges, too.  simpleCircle method isn't thread-safe, so
    // this part is single-threaded.
    //std::cout << Date() << ": finding smooth circles." << std::endl;
    for ( auto const& hhs : dict )
        for ( auto const& entry : hhs )
            if ( entry.getKDef().isNull() )
                eb.simpleCircle(entry);
    //std::cout << "There are " << pEdges->size()-nRegularEdges
    //          << " circular edges of total length "
    //          << pEdges->SizeSum()-nTotalLength << '.' << std::endl;
}

template <class Itr1, class Itr2>
size_t matchLen( Itr1 itr1, Itr1 end1, Itr2 itr2, Itr2 end2 )
{
    size_t result = 0;

    while ( itr1 != end1 && itr2 != end2 && *itr1 == *itr2 )
    { ++result; ++itr1; ++itr2; }

    return result;
}


template <class Itr1, class Itr2, class ItrW, typename sumType >
size_t fuzzyMatchLen( Itr1 itr1, Itr1 end1, Itr2 itr2, Itr2 end2,
        ItrW itrw, ItrW endw, sumType& maxWeight)
{
    size_t result = 0;

    while ( itr1 != end1 && itr2 != end2 && itrw != endw )
    {
        if ( *itr1 != *itr2 ) maxWeight -= *itrw;

        if ( maxWeight < 0 ) break;             // EARLY EXIT!

        ++result; ++itr1; ++itr2; ++itrw;
    }

    return result;
}

// fuzzyMatchLenBi - bidirectional matcher
//
// note the semantics:
//
// forward is the usual iterator traversal from its initial position
// up to, but not including, the end.
//
// backward first compares *(--itr) and continues backwards up to and
// including the passed in "end" position (which should be something like the
// beginning of the read or edge)
template <class Itr1, class Itr2, class ItrW, typename sumType >
size_t fuzzyMatchLenBi( Itr1 itr1, Itr1 end1, Itr2 itr2, Itr2 end2,
        ItrW itrw, ItrW endw, sumType& maxWeight, bool backward=false)
{
    size_t result = 0;
    int offset  = backward ? -1 : 0 ;
    int incr    = backward ? -1 : 1 ;

    while ( itr1 != end1 && itr2 != end2 && itrw != endw )
    {
        if ( itr1[offset] != itr2[offset] ) maxWeight -= *itrw;

        if ( maxWeight < 0 ) break;             // EARLY EXIT!

        ++result;
        itr1 += incr; itr2 += incr; itrw += incr;
    }

    return result;
}


class EdgeLoc
{
public:
    EdgeLoc() : mEdgeOffset(0) {}
    EdgeLoc( EdgeID const& edgeID, bool rc, unsigned edgeOffset )
    : mEdgeID(edgeID), mEdgeOffset(edgeOffset)
    { if ( rc ) mEdgeOffset = ~edgeOffset; }

    EdgeID const& getEdgeID() const { return mEdgeID; }
    bool isRC() const { return mEdgeOffset < 0; }
    unsigned getEdgeOffset() const
    { return mEdgeOffset < 0 ? ~mEdgeOffset : mEdgeOffset; }

    size_t hash() const { return (mEdgeID.val() << 32) | mEdgeOffset; }

    friend bool operator<( EdgeLoc const& el1, EdgeLoc const& el2 )
    { if ( el1.mEdgeID < el2.mEdgeID ) return true;
      if ( el2.mEdgeID < el1.mEdgeID ) return false;
      return el1.mEdgeOffset < el2.mEdgeOffset; }

    friend ostream& operator<<( ostream& os, EdgeLoc const& el )
    { if ( el.isRC() ) os << '~';
      return os << el.getEdgeID() << '[' << el.getEdgeOffset() << ']'; }

private:
    EdgeID mEdgeID;
    int mEdgeOffset;
};

class PathPart
{
public:
    PathPart() : mLen(0u), mEdgeLen(0) {}
    PathPart( unsigned gapLen ) : mLen(gapLen), mEdgeLen(0) {}
    PathPart( EdgeID const& edgeID, bool rc, unsigned edgeOffset,
                                                unsigned len, unsigned edgeLen )
    : mEdgeLoc(edgeID,rc,edgeOffset), mLen(len), mEdgeLen(edgeLen)
    { AssertLt(getEdgeOffset(),mEdgeLen); }

    EdgeID const& getEdgeID() const { return mEdgeLoc.getEdgeID(); }
    bool isRC() const { return mEdgeLoc.isRC(); }
    unsigned getEdgeOffset() const { return mEdgeLoc.getEdgeOffset(); }

    // if true, then only the length is meaningful
    bool isGap() const { return !mEdgeLen; }
    unsigned getLength() const { return mLen; }
    void incrLength(unsigned len) { mLen += len; }

    unsigned getEndOffset() const { return getEdgeOffset()+mLen; }

    // note: edge length is measured in kmers, not bases
    unsigned getEdgeLen() const { return mEdgeLen; }

    EdgeLoc firstLoc() const
    { return mEdgeLoc; }

    EdgeLoc lastLoc() const
    { EdgeLoc el(getEdgeID(),isRC(),getEdgeOffset()+mLen-1);
      AssertLt(el.getEdgeOffset(),mEdgeLen);
      return el; }

    bool isSameEdge( PathPart const& that ) const
    { return getEdgeID() == that.getEdgeID() && isRC() == that.isRC(); }

    bool isConformingCapturedGap( unsigned maxJitter ) const
    { PathPart const& prev = this[-1];
      PathPart const& next = this[1];
      unsigned graphDist = next.getEdgeOffset()-prev.getEndOffset();
      if ( !prev.isSameEdge(next) )
        graphDist += prev.getEdgeLen();
      // if the gap size is consistent with the graph, return true
      return unsigned(std::abs(int(mLen-graphDist))) <= maxJitter; }

    PathPart rc() const
    { return isGap() ? PathPart(mLen) :
       PathPart(getEdgeID(),!isRC(),mEdgeLen-getEndOffset(),mLen,mEdgeLen); }

    friend ostream& operator<<( ostream& os, PathPart const& ppt )
    { if ( ppt.isGap() ) return os << ppt.getLength() << 'X';
      if ( ppt.isRC() ) os << '~';
      int offset = ppt.getEdgeOffset();
      return os << ppt.getEdgeID()
                  << '[' << offset << '-' << offset+ppt.getLength()
                  << ':' << ppt.getEdgeLen() << ']'; }

private:
    EdgeLoc mEdgeLoc;
    unsigned mLen;
    unsigned mEdgeLen;
};

class Pather
{
public:
    Pather( Dict const& dict, vecbvec const& edges )
    : mDict(dict), mEdges(edges) {}

    std::vector<PathPart> const& path( bvec const& read )
    { mPathParts.clear();

      if ( read.size() < K )
      { mPathParts.emplace_back(read.size());
        return mPathParts; } // EARLY RETURN!

      auto itr = read.begin();
      auto end = read.end()-K+1;
      while ( itr != end )
      { Kmer kmer(itr);
        Entry const* pEnt = mDict.findEntry(kmer);
        if ( !pEnt )
        {
          unsigned gapLen = 1u;
          auto itr2 = itr+K; ++itr;
          auto end2 = read.end();
          while ( itr2 != end2 )
          { kmer.toSuccessor(*itr2); ++itr2;
            if ( (pEnt = mDict.findEntry(kmer)) )
                break;
            ++gapLen; ++itr; }
          mPathParts.emplace_back(gapLen); }
        if ( pEnt )
        { KDef const& kDef = pEnt->getKDef();
          bvec const& edge = mEdges[kDef.getEdgeID().val()];
          int offset = kDef.getEdgeOffset();
          auto eBeg(edge.begin(offset));
          size_t len = 1u;
          bool rc = CF<K>::isRC(itr,eBeg);
          if ( !rc )
            len += matchLen(itr+K,read.end(),eBeg+K,edge.end());
          else
          { offset = edge.size() - offset;
            auto eBegRC(edge.rcbegin(offset));
            len += matchLen(itr+K,read.end(),eBegRC,edge.rcend());
            offset = offset - K; }
          unsigned edgeKmers = edge.size()-K+1;

          mPathParts.emplace_back(kDef.getEdgeID(),rc,offset,len,edgeKmers);
          itr += len;
        } }
      return mPathParts; }

    std::vector<PathPart> const& path_careful(bvec const& read) {
        mPathParts.clear();

        if (read.size() < K) {
            mPathParts.emplace_back(read.size());
            return mPathParts;
        } // EARLY RETURN!

        auto itr = read.begin();
        auto end = read.end() - K + 1;
        while (itr != end) {
            Kmer kmer(itr);
            Entry const* pEnt = mDict.findEntry(kmer);
            if (!pEnt) {
                unsigned gapLen = 1u;
                auto itr2 = itr + K;
                ++itr;
                auto end2 = read.end();
                while (itr2 != end2) {
                    kmer.toSuccessor(*itr2);
                    ++itr2;
                    if ((pEnt = mDict.findEntry(kmer)))
                        break;
                    ++gapLen;
                    ++itr;
                }
                if ( mPathParts.size() && mPathParts.back().isGap() )
                    mPathParts.back().incrLength(gapLen);
                else
                    mPathParts.emplace_back(gapLen);
            }
            if (pEnt) {
                KDef const& kDef = pEnt->getKDef();
                bvec const& edge = mEdges[kDef.getEdgeID().val()];
                int offset = kDef.getEdgeOffset();
                auto eBeg(edge.begin(offset));
                size_t len = 1u;
                bool rc = CF<K>::isRC(itr, eBeg);
                if (!rc)
                    len += matchLen(itr + K, read.end(), eBeg + K, edge.end());
                else {
                    offset = edge.size() - offset;
                    auto eBegRC(edge.rcbegin(offset));
                    len += matchLen(itr + K, read.end(), eBegRC, edge.rcend());
                    offset = offset - K;
                }
                unsigned edgeKmers = edge.size() - K + 1;

                if ( len >= 5 )
                    mPathParts.emplace_back(kDef.getEdgeID(), rc, offset, len, edgeKmers);
                else if (mPathParts.size() && mPathParts.back().isGap())
                    mPathParts.back().incrLength(len);
                else
                    mPathParts.emplace_back(len);

                itr += len;
            }
        }
        return mPathParts;
    }

    bool isJoinable( PathPart const& pp1, PathPart const& pp2 ) const
    { if ( pp1.getEdgeID() == pp2.getEdgeID() ) return true;
      bvec const& e1 = mEdges[pp1.getEdgeID().val()];
      bvec const& e2 = mEdges[pp2.getEdgeID().val()];
      SubKmer k1 = pp1.isRC() ? SubKmer(e1.rcend()-K+1) : SubKmer(e1.end()-K+1);
      SubKmer k2 = pp2.isRC() ? SubKmer(e2.rcend()-K+1) : SubKmer(e2.end()-K+1);
      return k1==k2; }

    std::vector<PathPart> const& path_more( bvec const& read, qvec const& qual,
            bool trace = false ) {
        mPathParts.clear();

        if ( read.size() < K ) {
            mPathParts.emplace_back(read.size());
            return mPathParts;                          // EARLY RETURN!
        }

        auto itr  = read.begin();
        auto qitr = qual.begin();
        auto end = read.end()-K+1;


        while ( itr != end ) {
            int qSumRemain = 21;        // reset per edge
            Kmer kmer(itr);
            Entry const* pEnt = mDict.findEntry(kmer);
            int gapLen = 0;
            if ( !pEnt ) {
                ++gapLen;
                auto itr2 = itr+K; ++itr;
                auto end2 = read.end();
                while ( itr2 != end2 ) {
                    kmer.toSuccessor(*itr2); ++itr2;
                    if ( (pEnt = mDict.findEntry(kmer)) )
                        break;
                    ++gapLen; ++itr;
                }
            }
            // if still !pEnt, then we found no good kmers, so register the gap
            // otherwise we look to extend the match in both directions and possibly
            // adjust the gap later.
            if ( !pEnt ) {
                mPathParts.emplace_back(gapLen);
            } else {
                KDef const& kDef = pEnt->getKDef();
                bvec const& edge = mEdges[kDef.getEdgeID().val()];
                int offset = kDef.getEdgeOffset();
                auto eBeg(edge.begin(offset));
                size_t len = 1u;
                size_t backw = 0u;
                bool rc = CF<K>::isRC(itr,eBeg);
                if ( !rc ) {
                    len += fuzzyMatchLen(itr+K,read.end(),eBeg+K,edge.end(),
                            qitr+K, qual.end(), qSumRemain);
                    backw = fuzzyMatchLenBi( itr, read.begin(), eBeg,
                            edge.begin(), qitr, qual.begin(), qSumRemain, true );
                    offset -= backw;
                    gapLen -= backw;
                } else {
                    offset = edge.size() - offset;
                    auto eBegRC(edge.rcbegin(offset));
                    len += fuzzyMatchLen(itr+K,read.end(),eBegRC,edge.rcend(),
                            qitr+K,qual.end(),qSumRemain);
                    offset = offset - K;

                    backw = fuzzyMatchLenBi( itr, read.begin(), eBegRC-K, edge.rcbegin(),
                            qitr, qual.begin(), qSumRemain, true);
                    offset -= backw;
                    gapLen -= backw;
                }
                if ( gapLen > 0 ) mPathParts.emplace_back(gapLen);
                else if ( gapLen < 0 )
                    cout << "ARRRRRRGH!   gapLen=" << gapLen << endl;
                unsigned edgeKmers = edge.size()-K+1;
                mPathParts.emplace_back(kDef.getEdgeID(),rc,offset,len+backw,edgeKmers);
                itr += len;
            }

        }

        if ( trace )
            for ( auto const& part : mPathParts )
                cout << part << endl;

        return mPathParts;
    }


    friend ostream& operator<<( ostream& os, Pather const& pather )
    { os << rangePrinter(pather.mPathParts.begin(),pather.mPathParts.end()," ");
      return os; }

private:
    Dict const& mDict;
    vecbvec const& mEdges;
    std::vector<PathPart> mPathParts;
};

class GapFiller
{
public:
    GapFiller( vecbvec const& reads, vecbvec const& edges, unsigned maxGapSize,
                    unsigned minFreq, Dict* pDict )
    : mReads(reads), mDict(*pDict), mPather(*pDict,edges),
      mMaxGapSize(maxGapSize), mMinFreq(minFreq) {}

    template <class OItr>
    void map( size_t readId, OItr oItr )
    { bvec const& read = mReads[readId];
      std::vector<PathPart> const& parts = mPather.path(read);
      if ( parts.size() < 3 ) return; // EARLY RETURN!
      auto rItr = read.begin(parts[0].getLength());
      auto end = parts.end()-1;
      // for each path part sandwiched between two others
      for ( auto pPart=parts.begin()+1; pPart != end; ++pPart )
      { // if it's not a gap of appropriate size
        if ( !pPart->isGap() ||
               (mMaxGapSize && pPart->getLength() > mMaxGapSize) ||
               // or if it fits cleanly on the existing graph
               pPart->isConformingCapturedGap(MAX_JITTER) )
        { rItr += pPart->getLength();
          continue; } // just CONTINUE
        // we have a gap of appropriate size, and it has new, interesting kmers
        // that might change the graph structure
        auto itr = rItr - 1;
        Kmer kkk(itr); itr += K;
        update(kkk,KMerContext::initialContext(*itr));
        auto last = itr + pPart->getLength();
        while ( itr != last )
        { unsigned char predCode = kkk.front();
          kkk.toSuccessor(*itr);
          KMerContext kc(predCode,*++itr);
          *oItr++ = kkk.isRev()?Entry(Kmer(kkk).rc(),kc.rc()) : Entry(kkk,kc); }
        KMerContext kc = KMerContext::finalContext(kkk.front());
        kkk.toSuccessor(*itr);
        update(kkk,kc);
        rItr += pPart->getLength(); } }

    void reduce( Entry* e1, Entry* e2 )
    { summarizeEntries(e1,e2);
      if ( e1->getKDef().getCount() >= mMinFreq )
      { mDict.insertEntry(std::move(*e1)); } }

    Entry* overflow( Entry* e1, Entry* e2 )
    { if ( e2-e1 > 1 ) summarizeEntries(e1,e2); return e1+1; }

private:
    void update( Kmer kmer, KMerContext kc )
    { if ( kmer.isRev() ) { kmer.rc(); kc=kc.rc(); }
      mDict.applyCanonical(kmer,
              [kc](Entry const& ent)
              { KDef& kd = const_cast<KDef&>(ent.getKDef());
                kd.setContext(kc|kd.getContext()); }); }

    static unsigned const MAX_JITTER = 1;
    vecbvec const& mReads;
    Dict& mDict;
    Pather mPather;
    unsigned mMaxGapSize;
    unsigned mMinFreq;
};
typedef MapReduceEngine<GapFiller,Entry,Kmer::Hasher> GFMRE;

void fillGaps( vecbvec const& reads, unsigned maxGapSize, unsigned minFreq,
                    vecbvec* pEdges, Dict* pDict )
{
    //std::cout << Date() << ": filling gaps." << std::endl;
    GapFiller gf(reads,*pEdges,maxGapSize,minFreq,pDict);
    GFMRE mre(gf);
    if ( !mre.run(5*reads.size(),0ul,reads.size()) )
        FatalErr("Map/Reduce operation failed in gap filling.");

    //std::cout << "Now the dictionary has " << pDict->size()
    //                << " entries." << std::endl;

    //std::cout << Date() << ": finding edge sequences again." << std::endl;
    pDict->nullEntries();
    pDict->recomputeAdjacencies();
    pEdges->clear();
    buildEdges(*pDict,pEdges);
}

void pathRef( String const& refFasta, Dict const& dict, vecbvec const& edges )
{
    vecbvec ref;
    FastFetchReads(ref,nullptr,refFasta);
    Pather pather(dict,edges);
    for ( bvec const& refTig : ref )
    {
        pather.path(refTig);
        std::cout << "Reference path:" << std::endl;
        std::cout << pather << std::endl;
    }
    size_t edgeId = 0;
    for ( bvec const& edge : edges )
    {
        std::cout << edgeId++ << ": ";
        if ( edge.size() < 2*K )
            std::cout << edge;
        else
            std::cout << bvec(edge.begin(),edge.begin(K)) << "..."
                       << bvec(edge.end()-K,edge.end());
        std::cout << ' ' << edge.size()-K+1 << std::endl;
    }
}

class Join
{
public:
    Join( EdgeLoc const& edgeLoc1, EdgeLoc const& edgeLoc2, unsigned overlap )
    : mEdgeLoc1(edgeLoc1), mEdgeLoc2(edgeLoc2), mOverlap(overlap)
    {}

    EdgeLoc const& getEdgeLoc1() const { return mEdgeLoc1; }
    EdgeLoc const& getEdgeLoc2() const { return mEdgeLoc2; }
    unsigned getOverlap() const { return mOverlap; }

    size_t hash() const
    { return 47 * (47*mEdgeLoc1.hash()+mEdgeLoc2.hash()) + mOverlap; }

    struct Hasher
    { size_t operator()( Join const& join ) { return join.hash(); } };

    friend bool operator<( Join const& join1, Join const& join2 )
    { if ( join1.mEdgeLoc1 < join2.mEdgeLoc1 ) return true;
      if ( join2.mEdgeLoc1 < join1.mEdgeLoc1 ) return false;
      if ( join1.mEdgeLoc2 < join2.mEdgeLoc2 ) return true;
      if ( join2.mEdgeLoc2 < join1.mEdgeLoc2 ) return false;
      return join1.mOverlap < join2.mOverlap; }

    friend ostream& operator<<( ostream& os, Join const& join )
    { return os << join.mEdgeLoc1 << "<-" << join.mOverlap << "->"
                    << join.mEdgeLoc2; }

private:
    EdgeLoc mEdgeLoc1;
    EdgeLoc mEdgeLoc2;
    unsigned mOverlap;
};

class Joiner
{
public:
    Joiner( vecbvec const& reads, vecbvec const& edges, Dict const& dict,
                    unsigned maxGapSize, unsigned minFreq,
                    vecbvec* pFakeReads )
    : mReads(reads), mEdges(edges), mPather(dict,edges),
      mMaxGapSize(maxGapSize), mMinFreq(minFreq), mFakeReads(*pFakeReads)
    { ForceAssertLt(maxGapSize,K-1); }

    template <class OItr>
    void map( size_t readId, OItr oItr )
    { bvec const& read = mReads[readId];
      std::vector<PathPart> const& parts = mPather.path(read);
      if ( parts.size() < 3 ) return; // EARLY RETURN!
      auto end = parts.end()-1;
      for ( auto pPart=parts.begin()+1; pPart != end; ++pPart )
      { // if it's a captured gap of appropriate size
        if ( pPart->isGap() && pPart->getLength() <= mMaxGapSize )
        { PathPart const& prev = pPart[-1];
          PathPart const& next = pPart[1];
          unsigned overlap = K-pPart->getLength()-1;
          if ( next.getEdgeID() < prev.getEdgeID() )
            *oItr++ = Join(next.rc().lastLoc(),prev.rc().firstLoc(),overlap);
          else
            *oItr++ = Join(prev.lastLoc(),next.firstLoc(),overlap); } } }

    void reduce( Join* pJoin1, Join* pJoin2 )
    { Assert(validOverlap(*pJoin1));
      if ( pJoin2-pJoin1 >= mMinFreq )
      { //std::cout << "Joining: " << *pJoin1 << std::endl;
        bvec bv; bv.reserve(2*K);
        append(pJoin1->getEdgeLoc1(),0,&bv);
        append(pJoin1->getEdgeLoc2(),pJoin1->getOverlap(),&bv);
        add(bv); } }

    Join* overflow( Join* pJoin1, Join* pJoin2 )
    { return pJoin1+std::min(unsigned(pJoin2-pJoin1),mMinFreq); }

private:
    bool validOverlap( Join const& join )
    { EdgeLoc const& el1 = join.getEdgeLoc1();
      bvec const& bv1 = mEdges[el1.getEdgeID().val()];
      EdgeLoc const& el2 = join.getEdgeLoc2();
      bvec const& bv2 = mEdges[el2.getEdgeID().val()];
      bool result;
      if ( el1.isRC() )
      { auto end = bv1.rcbegin(el1.getEdgeOffset()+K);
        auto beg = end - join.getOverlap();
        if ( el2.isRC() )
          result = std::equal(beg,end,bv2.rcbegin(el2.getEdgeOffset()));
        else
          result = std::equal(beg,end,bv2.begin(el2.getEdgeOffset())); }
      else
      { auto end = bv1.begin(el1.getEdgeOffset()+K);
        auto beg = end - join.getOverlap();
        if ( el2.isRC() )
          result = std::equal(beg,end,bv2.rcbegin(el2.getEdgeOffset()));
        else
          result = std::equal(beg,end,bv2.begin(el2.getEdgeOffset())); }
      return result; }

    void append( EdgeLoc const& el, unsigned indent, bvec* pBV )
    { bvec const& edge = mEdges[el.getEdgeID().val()];
      unsigned offset = el.getEdgeOffset()+indent;
      unsigned len = K-indent;
      if ( el.isRC() )
      { auto itr = edge.rcbegin(offset);
        pBV->append(itr,itr+len); }
      else
      { auto itr = edge.begin(offset);
        pBV->append(itr,itr+len); } }

    void add( bvec const& bv )
    { static SpinLockedData gLock;
      if ( true )
      { SpinLocker lock(gLock);
        mFakeReads.push_back(bv); } }

    vecbvec const& mReads;
    vecbvec const& mEdges;
    Pather mPather;
    unsigned mMaxGapSize;
    unsigned mMinFreq;
    vecbvec& mFakeReads;
};
typedef MapReduceEngine<Joiner,Join,Join::Hasher> JMRE;

void joinOverlaps( vecbvec const& reads, unsigned maxGapSize, unsigned minFreq,
                        vecbvec* pEdges, Dict* pDict )
{
    //std::cout << Date() << ": joining overlaps." << std::endl;
    vecbvec fakeReads;
    fakeReads.reserve(pEdges->size()/10);
    Joiner joiner(reads,*pEdges,*pDict,maxGapSize,minFreq,&fakeReads);
    JMRE mre(joiner);
    if ( !mre.run(reads.size(),0ul,reads.size()) )
        FatalErr("Map/Reduce operation failed when joining overlaps.");

    if ( fakeReads.size() )
    {
        //std::cout << "Found " << fakeReads.size() << " joins." << std::endl;
        pDict->nullEntries();
        pDict->process(fakeReads);
        //std::cout << "Now the dictionary has " << pDict->size()
        //                    << " entries." << std::endl;
        //std::cout << Date() << ": finding edge sequences again." << std::endl;
        pEdges->clear();
        buildEdges(*pDict,pEdges);
    }
}


class HBVPather
{
public:
    typedef enum {ALGORITHM_ONE, ALGORITHM_TWO} Algorithm;

    HBVPather( vecbvec const& reads, VecPQVec const& quals,
                Dict const& dict, vecbvec const& edges,
                HyperBasevector const& hbv,
                vec<int> const& fwdEdgeXlat, vec<int> const& revEdgeXlat,
                Algorithm alg, ReadPathVec* pPaths, Bool const verbose = False )
    : mReads(reads), mQuals(quals), mHBV(hbv), mFwdEdgeXlat(fwdEdgeXlat),
      mRevEdgeXlat(revEdgeXlat), mPather(dict,edges), mAlgorithm(alg),
      mPaths(*pPaths), mVerbose(verbose),
      mExtender(mHBV,&mToLeft,&mToRight,mVerbose) {
        mHBV.ToLeft(mToLeft);
        mHBV.ToRight(mToRight);
    }


    void operator()( size_t readId )
    {
        switch (mAlgorithm) {
        case ALGORITHM_ONE:
            this->algorithmOne(readId);
            break;
        case ALGORITHM_TWO:
            this->algorithmTwo(readId);
            break;
        default:
            FatalErr("BUG: bad algorithm id passed to HBVPather::path()");
            break;
        }
    }

    size_t startOnRead( std::vector<PathPart> const& parts )
    {
        if ( parts.size() > 1 && parts[0].isGap() && parts[0].getLength() > parts[1].getEdgeOffset() )
            return parts[0].getLength() - parts[1].getEdgeOffset();
        else
            return 0;
    }


    void algorithmTwo( size_t readId )
    {
        bool const debug = mVerbose;
        std::vector<PathPart> parts = mPather.path(mReads[readId]);     // needs to become a forward_list

        if ( debug )  {
            cout << endl;
            cout << readId << " A2: startOnRead=" << startOnRead(parts) << " ";
            if ( parts[0].isGap() ) cout << " gap=" << parts[0].getLength() << endl;
            else cout << endl;
            for ( auto const& part : parts ) cout << part << endl;
            pathPartsToReadPath(parts, mPath);
            PRINT(mPath);
        }


        // convert any seeds on hanging edges to gaps
        std::vector<PathPart> new_parts;
        for ( auto part : parts ) {
            if ( !part.isGap() ) {
               size_t edge_id = pathPartToEdgeID(part);
               size_t vleft  = mToLeft[edge_id];
               size_t vright = mToRight[edge_id];
               if ( mHBV.ToSize(vleft) == 0
                   && mHBV.ToSize(vright) > 1
                   && mHBV.FromSize(vright) > 0
                   && part.getEdgeLen() <= 100 ) {
                       // delete a seed on a hanging edge
                   part = PathPart(part.getLength() );
                   if ( debug )
                       cout << "delete SEED on hanging edge HBV edge_id=" <<
                           edge_id << ", orig seed=" << part.getEdgeID().val() << endl;
               }
            }

            // either add to an existing gap or create a new part (gap or not)
            if ( part.isGap() && new_parts.size() && new_parts.back().isGap() )
                new_parts.back().incrLength( part.getLength() );
            else
                new_parts.push_back(part);
        }
        std::swap(parts,new_parts);

        // if a gap is captured between two edges and things don't make sense, delete the
        // seed, if there are more than one, and subsequent seeds.  This hopefully avoids 
        // the loss of sensitivity of just dropping the seeds.
        if ( parts.size() >= 3 ) {
            size_t seeds = ( parts.begin()->isGap()) ? 0u : 1u;
            for ( auto pPart = parts.begin()+1; pPart != parts.end()-1; ++pPart ) {
                if ( !pPart->isGap() ) { seeds++; continue; }
                if ( !pPart->isConformingCapturedGap(MAX_JITTER) ||
                        !mPather.isJoinable(pPart[-1],pPart[1]) ) {


                    if ( seeds > 1 ) {
                        PathPart tmpPart(pPart[-1].getLength());
                        for ( auto pPart2 = pPart; pPart2 != parts.end(); ++ pPart2 )
                            tmpPart.incrLength( pPart2->getLength() );
                        parts.erase( pPart-1, parts.end() );
                        parts.push_back(tmpPart);
                    } else {
                        for ( auto pPart2 = pPart+1; pPart2 != parts.end(); ++ pPart2 )
                            pPart->incrLength( pPart2->getLength() );
                        parts.erase( pPart+1, parts.end() );
                    }

                    break;
                }
            }
        }


        // if a seed has only taken us <= 5 kmers onto an edge, we back off
        // because we probably don't have enough evidence to commit.  Extension
        // code later will be more careful.
        if ( parts.back().isGap() && parts.size() > 1 ) {
            auto const& last2 = parts[parts.size()-2];
            if ( last2.getEdgeOffset() == 0 && last2.getLength() <= 5 ) {
                auto last = parts.back();
                last.incrLength(last2.getLength());
                parts.pop_back();
                parts.pop_back();
                parts.push_back(last);
            }
        } else if ( !parts.back().isGap() ) {
            auto& last = parts.back();
            if ( last.getEdgeOffset()==0 && last.getLength() <= 5 )  {
                last = PathPart( last.getLength() );
            }
        }

        pathPartsToReadPath(parts, mPath);

        if ( debug ) {
            cout << "PARTS AFTER EDITING: " << endl;
            for ( auto const& part: parts )
                cout << part << endl;
            PRINT(mPath);
        }

        mQuals[readId].unpack(&mQV);
        mExtender.attemptLeftRightExtension(mPath,mReads[readId],mQV);

        if ( debug )
            cout << "id #" << readId << ": new ReadPath=" << mPath << endl;
        mPaths[readId] = mPath;
    }

    void algorithmOne( size_t readId )
    {
        bool const debug = mVerbose;

        std::vector<PathPart> const& parts = mPather.path(mReads[readId]);

        if ( debug )  {
            cout << endl;
            cout << readId << " A1: startOnRead=" << startOnRead(parts) << " ";
            if ( parts[0].isGap() ) cout << " gap=" << parts[0].getLength() << endl;
            else cout << endl;
            for ( auto const& part : parts ) cout << part << endl;
            pathPartsToReadPath(parts, mPath);
            PRINT(mPath);
        }


      // for each path part sandwiched between two others
      // if the part is a captured gap, check that it can be squeezed out
      // without a lot of drama
      if ( parts.size() >= 3 )
      { auto end = parts.end()-1;
        for ( auto pPart=parts.begin()+1; pPart != end; ++pPart )
        { if ( !pPart->isGap() ) continue;
          if ( !pPart->isConformingCapturedGap(MAX_JITTER) ) return;
          if ( !mPather.isJoinable(pPart[-1],pPart[1]) ) return; } }
      mPath.clear();

      pathPartsToReadPath(parts, mPath);

      if ( debug ) {
          cout << "PARTS AFTER EDITING: " << endl;
          for ( auto const& part: parts )
              cout << part << endl;
          PRINT(mPath);
      }

      mPaths[readId] = mPath;
    }

    size_t pathPartToEdgeID( PathPart const& part )
    {
        ForceAssert(!part.isGap());
        size_t idx = part.getEdgeID().val();
        return ( part.isRC() ? mRevEdgeXlat[idx] : mFwdEdgeXlat[idx] );
    }

    void pathPartsToReadPath( std::vector<PathPart> const& parts, ReadPath& path )
    {
      path.clear();
      PathPart const* pLast = 0;
      for ( PathPart const& part : parts )
      { if ( part.isGap() ) continue;
        if ( pLast && pLast->isSameEdge(part) ) continue;
//        if ( part.getEdgeOffset() == 0 && part.getLength() < 5) continue;
        size_t idx = part.getEdgeID().val();
        path.push_back(part.isRC()?mRevEdgeXlat[idx]:mFwdEdgeXlat[idx]);
        pLast = &part; }
      if ( path.empty() )
      { path.setOffset(0); /*path.setLastSkip(0);*/ }
      else
      { PathPart const& firstPart = parts.front();
        if ( !firstPart.isGap() )
          path.setOffset(firstPart.getEdgeOffset());
        else
        {
            int eo1 = parts[1].getEdgeOffset();
            int eo2 = firstPart.getLength();
            path.setOffset(eo1 - eo2 );
        }
        /*
        PathPart const& lastPart = parts.back();

        if ( !lastPart.isGap() )
          path.setLastSkip(lastPart.getEdgeLen()-lastPart.getEndOffset());
        else
        { PathPart const& prevPart = parts[parts.size()-2];
          unsigned skip = prevPart.getEdgeLen()-prevPart.getEndOffset();
          skip -= std::min(skip,lastPart.getLength());
          path.setLastSkip(skip); }
       */
       }
    }

private:
    static unsigned const MAX_JITTER = 3;
    vecbvec const& mReads;
    VecPQVec const& mQuals;
    HyperBasevector const& mHBV;
    vec<int> mToLeft, mToRight;
    vec<int> const& mFwdEdgeXlat;
    vec<int> const& mRevEdgeXlat;
    Pather mPather;
    Algorithm mAlgorithm;
    ReadPathVec& mPaths;
    ReadPath mPath;
    Bool mVerbose;
    ExtendReadPath mExtender;
    qvec mQV;
};

void pathReads( vecbvec const& reads, VecPQVec const& quals,
                Dict const& dict, vecbvec const& edges,
                HyperBasevector const& hbv,
                vec<int> const& fwdEdgeXlat, vec<int> const& revEdgeXlat,
                ReadPathVec* pPaths, Bool const NEW_ALIGNER = False,
                Bool const VERBOSE = False )
{
    pPaths->clear().resize(reads.size());
    auto algorithm = NEW_ALIGNER ? HBVPather::ALGORITHM_TWO : HBVPather::ALGORITHM_ONE;
    HBVPather pather(reads,quals,dict,edges,hbv,fwdEdgeXlat,revEdgeXlat,
              algorithm,pPaths,VERBOSE);
    parallelForBatch(0ul,reads.size(),100000,pather);
}

void repathUnpathedReads(const vecbvec& reads, const VecPQVec& quals,
			 const HyperBasevector& hbv, ReadPathVec& paths) {

    vec<size_t> ids;
    for ( size_t i = 0; i < paths.size(); ++i)
	if (paths[i].empty())
	    ids.push_back(i);

    int num_threads = getConfiguredNumThreads();
    ShortKmerReadPather::FindPaths(reads, quals, hbv, paths, 24, num_threads, ids , true);
}

} // end of anonymous namespace

void buildReadQGraph40( vecbvec const& reads, ObjectManager<VecPQVec>& quals,
                        bool doFillGaps, bool doJoinOverlaps,
                        unsigned minQual, unsigned minFreq,
                        double minFreq2Fract, unsigned maxGapSize,
                        String const& refFasta,
                        bool useNewAligner, bool repathUnpathed,
                        HyperBasevector* pHBV, ReadPathVec* pPaths,
                        float const memFrac,
                        bool const VERBOSE )
{
    //std::cout << Date() << ": loading reads." << std::endl;

    Dict* pDict = createDict(reads,quals,minQual,minFreq,memFrac);

    // figure out the complete base sequence of each edge
    //std::cout << Date() << ": finding edge sequences." << std::endl;
    vecbvec edges;
    edges.reserve(pDict->size()/100);
    buildEdges(*pDict,&edges);

    unsigned minFreq2 = std::max(2u,unsigned(minFreq2Fract*minFreq+.5));

    if ( doFillGaps )
        fillGaps(reads,maxGapSize,minFreq2,&edges,pDict);

    if ( doJoinOverlaps )
        joinOverlaps(reads,K/2,minFreq2,&edges,pDict);

    if ( !refFasta.empty() )
        pathRef(refFasta,*pDict,edges);

    vec<int> fwdEdgeXlat;
    vec<int> revEdgeXlat;
    if ( !pPaths )
    {
        delete pDict;
        pDict = 0;
        buildHBVFromEdges(edges,K,pHBV,&fwdEdgeXlat,&revEdgeXlat);
    }
    else
    {
        buildHBVFromEdges(edges,K,pHBV,&fwdEdgeXlat,&revEdgeXlat);
        pathReads(reads,quals.load(),*pDict,edges,*pHBV,
                    fwdEdgeXlat,revEdgeXlat,pPaths,useNewAligner,VERBOSE);
        delete pDict;
        pDict = 0;

	if (repathUnpathed) 
	    repathUnpathedReads(reads,quals.load(),*pHBV, *pPaths);
    }
}
