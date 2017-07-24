///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * BigKPather.cc
 *
 *  Created on: Dec 13, 2013
 *      Author: tsharpe
 */

#include "kmers/BigKPather.h"
#include "feudal/HashSet.h"
#include "kmers/BigKMer.h"
#include "kmers/KmerRecord.h"
#include "kmers/KMerHasher.h"
#include "paths/long/HBVFromEdges.h"
#include "paths/long/LargeKDispatcher.h"
#include "system/SpinLockedData.h"
#include "system/WorklistN.h"
#include <iostream>
#include "10X/DfTools.h"
#include <unordered_map>
#include <unordered_set>
#include "feudal/SubsetMasterVec.h"
#include "10X/mergers/NicePrints.h"
#include "10X/Gap.h"

size_t checksumV(vec<int>& V)
{
    std::size_t seed = 4;
    for(auto v : V)
        seed ^= std::hash<int>()(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    return seed;
}

size_t checksumV(vec<int64_t>& V)
{
    std::size_t seed = 4;
    for(auto v : V)
        seed ^= std::hash<int64_t>()(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    return seed;
}

size_t checksumRP(ReadPathVec& allx_paths)
{
    std::size_t seed = 4;
    for(auto& rp: allx_paths){
        for(auto& e: rp){
            seed ^= std::hash<int>()(e) + 0x9e3779b9 + 
                (seed << 6) + (seed >> 2);
        } 
            seed ^= std::hash<int>()(rp.getOffset()) + 0x9e3779b9 + 
                (seed << 6) + (seed >> 2);
    }
    return seed;
}

namespace
{

template <unsigned BIGK>
using BigDict = HashSet<BigKMer<BIGK>,typename BigKMer<BIGK>::hasher>;


template <unsigned BIGK>
class BigKMerizer
{
public:
    typedef BigDict<BIGK> BigKDict;
    BigKMerizer( BigKDict* pDict ) : mDict(*pDict) {}

    void kmerize( bvec const& bv )
    { if ( bv.size() < BIGK ) return;
      KMerHasher<BIGK> hasher;
      auto itr = bv.begin();
      size_t hash = hasher.hash(itr);
      if ( bv.size() == BIGK ) { add(BigKMer<BIGK>(bv,hash)); return; }
      BigKMer<BIGK> kmer(bv,hash,KMerContext::initialContext(bv[BIGK]));
      add(kmer);
      auto end = bv.end()-BIGK;
      while ( ++itr != end )
      { hash = hasher.stepF(itr);
        kmer.successor(hash,KMerContext(itr[-1],itr[BIGK]));
        add(kmer); }
      hash = hasher.stepF(itr);
      kmer.successor(hash,KMerContext::finalContext(itr[-1]));
      add(kmer); }

    // update the dictionary with only first kmer
    void kmerize_sleek( bvec const& bv )
    { if ( bv.size() < BIGK ) return;
      KMerHasher<BIGK> hasher;
      auto itr = bv.begin();
      size_t hash = hasher.hash(itr);
      // stop at first kmer, don't add context if bv is size of BIGK
      if ( bv.size() == BIGK ) { add(BigKMer<BIGK>(bv,hash)); return; }
      BigKMer<BIGK> kmer(bv,hash,KMerContext::initialContext(bv[BIGK]));
      add(kmer); }

    // get edge and offset that first kmer of bv lands on
    void kmerize_getLanding( bvec const& bv, int* pStarter,
        vecbvec & pEdges)
    { if ( bv.size() < BIGK ) return;
      auto itr = bv.begin();
      // store into readpath
      KMerHasher<BIGK> hasher;
      BigKMer<BIGK> kmer(bv,hasher(bv.begin()));
      BigKMer<BIGK> entry = lookup(kmer);
      pStarter[0] = entry.getOffset();
      size_t idx = &entry.getBV()-&pEdges[0];
      // ** size_t to int conversion may lose info when nReads > intmax
      pStarter[1] = idx;
      pStarter[2] = entry.isRC();
    }

    void kmerize_getLanding_safe( bvec const& bv, int* pStarter,
        vecbvec & pEdges)
    { if ( bv.size() < BIGK ) return;
      auto itr = bv.begin();
      // store into readpath
      KMerHasher<BIGK> hasher;
      BigKMer<BIGK> kmer(bv,hasher(bv.begin()));
      BigKMer<BIGK> entry; 
      bool success = lookup_safe(kmer,entry);
      if(!success){ // lookup failed
          pStarter[0] = -1;
          pStarter[1] = -1;
          pStarter[2] = 1;
          return;
      }
      pStarter[0] = entry.getOffset();
      size_t idx = &entry.getBV()-&pEdges[0];
      // ** size_t to int conversion may lose info when nReads > intmax
      pStarter[1] = idx;
      pStarter[2] = entry.isRC();
    }

    void findEdgeIdx(const vecbvec& edges, const MasterVec<bvec>& oldEdges, 
            int64_t id, vec<triple<int64_t,int64_t,Bool>>& map)
    {
        if(id>=oldEdges.size()) return;
        const bvec& searchEdge = oldEdges[id];
        KMerHasher<BIGK> hasher;
        //!!! oldEdges should have kmers appearing in edges
        //get first kmer idx and offset
        BigKMer<BIGK> skmer(searchEdge,hasher(searchEdge.begin()));
        BigKMer<BIGK> entry = lookup(skmer); 
        size_t idx = &entry.getBV()-&edges[0];
        map[id].first = id;
        map[id].second = -1;
        map[id].third= entry.isRC();
        if(entry.isRC()){
            if(ReverseComplement(edges[idx])==searchEdge)
                map[id].second = idx;
        }else{
            if(edges[idx]==searchEdge)
                map[id].second = idx;
        }
    }

    template <class OItr>
    void map( vecbvec::const_iterator iItr, OItr oItr )
    { bvec const& bv = *iItr;
      if ( bv.size() < BIGK ) return;
      KMerHasher<BIGK> hasher;
      auto itr = bv.begin();
      size_t hash = hasher.hash(itr);
      if ( bv.size() == BIGK )
      { BigKMer<BIGK> kmer(bv,hash);
        *oItr = kmer.isRev() ? kmer.rc() : kmer; ++oItr;
        return; }
      BigKMer<BIGK> kmer(bv,hash,KMerContext::initialContext(bv[BIGK]));
      *oItr = kmer.isRev() ? kmer.rc() : kmer; ++oItr;
      auto end = bv.end()-BIGK;
      while ( ++itr != end )
      { hash = hasher.stepF(itr);
        kmer.successor(hash,KMerContext(itr[-1],itr[BIGK]));
        *oItr = kmer.isRev() ? kmer.rc() : kmer; ++oItr; }
      hash = hasher.stepF(itr);
      kmer.successor(hash,KMerContext::finalContext(itr[-1]));
      *oItr = kmer.isRev() ? kmer.rc() : kmer; ++oItr; }

    void reduce( BigKMer<BIGK>* key1, BigKMer<BIGK>* key2 )
    { for ( auto itr=key1+1; itr != key2; ++itr )
        key1->addContext(itr->getContext());
      mDict.insertUniqueValue(*key1); }

    // We build the dictionary with BigKMers pointing to reads.  This updates
    // all the BigKMers in the dictionary to point to edges instead.
    void updateDict( bvec const& edge )
    { KMerHasher<BIGK> hasher;
      auto itr = edge.begin();
      BigKMer<BIGK> kmer(edge,hasher.hash(itr));
      update(kmer.isRev()?kmer.rc():kmer);
      auto end = edge.end()-BIGK+1;
      while ( ++itr != end )
      { kmer.successor(hasher.stepF(itr));
        update(kmer.isRev()?kmer.rc():kmer); } }

private:
    void add( BigKMer<BIGK> const& kmer )
    { canonicalAdd( kmer.isRev() ? kmer.rc() : kmer ); }

    BigKMer<BIGK> lookup( BigKMer<BIGK> const& kmer )
    { if ( kmer.isRev() )
      { BigKMer<BIGK> const* pEnt = mDict.lookup(kmer.rc());
        ForceAssert(pEnt);
        return pEnt->rc(); }
      BigKMer<BIGK> const* pEnt = mDict.lookup(kmer);
      ForceAssert(pEnt);
      return *pEnt; }

    bool lookup_safe( BigKMer<BIGK> const& kmer, BigKMer<BIGK>& val)
    { if ( kmer.isRev() )
      { BigKMer<BIGK> const* pEnt = mDict.lookup(kmer.rc());
        if(!pEnt) return false;
        val = pEnt->rc();
        return true;}
      BigKMer<BIGK> const* pEnt = mDict.lookup(kmer);
      if(!pEnt) return false;
      val = *pEnt;
      return true; }

    void update( BigKMer<BIGK> const& kmer )
    { const_cast<BigKMer<BIGK>*>(mDict.lookup(kmer))->updateLocation(kmer); }

    void canonicalAdd( BigKMer<BIGK> const& kmer ) const
    { mDict.apply(kmer,
        [&kmer]( BigKMer<BIGK> const& entry )
        { const_cast<BigKMer<BIGK>&>(entry).addContext(kmer.getContext()); }); }

    BigKDict& mDict;
};

template <unsigned BIGK>
class BigKEdgeBuilder
{
public:
    typedef BigKMer<BIGK> BigKmer;
    typedef BigKMerizer<BIGK> BigKmerizer;
    typedef BigDict<BIGK> BigKDict;

    static void buildEdges( BigKDict const& dict, vecbvec* pEdges )
    {
        BigKEdgeBuilder eb(dict,pEdges);
        dict.parallelForEachHHS(
                [eb]( typename BigKDict::HHS const& hhs ) mutable
                { for ( BigKmer const& entry : hhs )
                    if ( entry.isUnassigned() )
                      eb.buildEdge(entry); });

        size_t nRegularEdges = pEdges->size();
        size_t nTotalLength = pEdges->SizeSum();

        // if a kmer isn't marked as being on an edge, it must be a part of a
        // smooth circle: add those edges, too.  simpleCircle method isn't
        // thread-safe, so this part is single-threaded.
        for ( auto const& hhs : dict )
            for ( auto const& entry : hhs )
                if ( entry.isUnassigned() )
                    eb.simpleCircle(entry);
    }

private:
    BigKEdgeBuilder( BigKDict const& dict, vecbvec* pEdges )
    : mDict(dict), mEdges(*pEdges) {}

    void buildEdge( BigKmer const& entry )
    { if ( isPalindrome(entry) )
        make1KmerEdge(entry);
      else if ( upstreamExtensionPossible(entry) )
      { if ( downstreamExtensionPossible(entry) )
          return;
        BigKmer rc = entry.rc();
        mEdgeSeq.assign(rc.begin(),rc.end());
        mEdgeEntries.push_back(&entry);
        extend(rc.getContext()); }
      else if ( downstreamExtensionPossible(entry) )
      { mEdgeSeq.assign(entry.begin(),entry.end());
        mEdgeEntries.push_back(&entry);
        extend(entry.getContext()); }
      else
        make1KmerEdge(entry); }

    // not thread-safe
    void simpleCircle( BigKmer const& entry )
    { BigKmer const* pFirstEntry = &entry;
      mEdgeSeq.assign(entry.begin(),entry.end());
      mEdgeEntries.push_back(pFirstEntry);
      auto itr = mEdgeSeq.begin();
      KMerHasher<BIGK> hasher;
      BigKmer kmer(mEdgeSeq,hasher.hash(itr));
      KMerContext context = entry.getContext();
      while ( true )
      { ForceAssertEq(context.getPredecessorCount(),1u);
        ForceAssertEq(context.getSuccessorCount(),1u);
        mEdgeSeq.push_back(context.getSingleSuccessor());
        kmer.successor(hasher.stepF(++itr));
        BigKmer const* pEntry = lookup(kmer,&context);
        if ( pEntry == pFirstEntry )
        { mEdgeSeq.pop_back();
          break; }
        if ( !pEntry->isUnassigned() )
        { std::cout << "Failed to close circle.\n";
          for ( auto beg=mEdgeEntries.begin(),end=mEdgeEntries.end();
                  beg != end; ++beg )
            std::cout << *beg << ' ' << **beg;
          std::cout << pEntry << ' ' << *pEntry << std::endl;
          CRD::exit(1); }
        mEdgeEntries.push_back(pEntry); }
      canonicalizeCircle();
      addEdge(); }

    void canonicalizeCircle()
    { auto itr = std::min_element(mEdgeEntries.begin(),mEdgeEntries.end(),
                            []( BigKmer const* pEnt1, BigKmer const* pEnt2 )
                            { return *pEnt1 < *pEnt2; });
      BigKmer const& minKmer = **itr;
      size_t idx = itr - mEdgeEntries.begin();
      if ( !std::equal(minKmer.begin(),minKmer.end(),mEdgeSeq.begin(idx)) )
      { mEdgeSeq.ReverseComplement();
        idx = mEdgeSeq.size()-idx-BIGK;
        Assert(std::equal(minKmer.begin(),minKmer.end(),mEdgeSeq.begin(idx))); }
      bvec bv; bv.reserve(mEdgeSeq.size());
      bv.assign(mEdgeSeq.begin(idx),mEdgeSeq.end());
      bv.append(mEdgeSeq.begin(BIGK-1),mEdgeSeq.begin(BIGK-1+idx));
      Assert(std::equal(mEdgeSeq.begin(),mEdgeSeq.begin(BIGK-1),
                          mEdgeSeq.end()-(BIGK-1)));
      mEdgeSeq = bv; }

    bool isPalindrome( BigKmer const& kmer )
    { auto itr = kmer.begin();
      if ( !(BIGK&1) )
        return CF<BIGK>::isPalindrome(itr);
      if ( CF<BIGK-1>::isPalindrome(itr) )
          return true;
      return CF<BIGK-1>::isPalindrome(++itr); }

    bool upstreamExtensionPossible( BigKmer const& entry )
    { KMerContext context = entry.getContext();
      if ( context.getPredecessorCount() != 1 )
        return false;
      mEdgeSeq.clear().push_back(context.getSinglePredecessor());
      mEdgeSeq.append(entry.begin(),entry.end()-1);
      if ( CF<BIGK>::isPalindrome(mEdgeSeq.begin()) )
        return false;
      BigKmer kmer(mEdgeSeq,KMerHasher<BIGK>()(mEdgeSeq.begin()));
      lookup(kmer,&context);
      return context.getSuccessorCount() == 1; }

    bool downstreamExtensionPossible( BigKmer const& entry )
    { KMerContext context = entry.getContext();
      if ( context.getSuccessorCount() != 1 )
        return false;
      mEdgeSeq.assign(entry.begin()+1,entry.end());
      mEdgeSeq.push_back(context.getSingleSuccessor());
      if ( CF<BIGK>::isPalindrome(mEdgeSeq.begin()) )
        return false;
      BigKmer kmer(mEdgeSeq,KMerHasher<BIGK>()(mEdgeSeq.begin()));
      lookup(kmer,&context);
      return context.getPredecessorCount() == 1; }

    void make1KmerEdge( BigKmer const& entry )
    { mEdgeSeq.assign(entry.begin(),entry.end());
      mEdgeEntries.push_back(&entry);
      addEdge(); }

    void extend( KMerContext context )
    { auto itr = mEdgeSeq.begin();
      KMerHasher<BIGK> hasher;
      BigKmer kmer(mEdgeSeq,hasher.hash(itr));
      while ( context.getSuccessorCount() == 1 )
      { mEdgeSeq.push_back(context.getSingleSuccessor());
        kmer.successor(hasher.stepF(++itr));
        if ( isPalindrome(kmer) )
        { mEdgeSeq.pop_back();
          break; }
        BigKmer const* pEntry = lookup(kmer,&context);
        if ( context.getPredecessorCount() != 1 )
        { mEdgeSeq.pop_back();
          break; }
        mEdgeEntries.push_back(pEntry); }
      switch (mEdgeSeq.getCanonicalForm())
      {case CanonicalForm::PALINDROME:
         ForceAssertEq(mEdgeSeq.size(),BIGK);
        // allow flow-through to FWD case
       case CanonicalForm::FWD:
         addEdge(); break;
       case CanonicalForm::REV:
         mEdgeSeq.clear(); mEdgeEntries.clear(); break; } }

    BigKmer const* lookup( BigKmer const& kmer, KMerContext* pContext )
    { BigKmer const* result;
      if ( kmer.isRev() )
      { result = mDict.lookup(kmer.rc());
        ForceAssert(result);
        *pContext = result->getContext().rc(); }
      else
      { result = mDict.lookup(kmer);
        ForceAssert(result);
        *pContext = result->getContext(); }
      return result; }

    void addEdge()
    { static SpinLockedData gLock;
      size_t edgeId;
      if ( mEdgeSeq.getCanonicalForm() == CanonicalForm::REV )
          mEdgeSeq.ReverseComplement();
      if ( true )
      { SpinLocker lock(gLock);
        edgeId = mEdges.size();
        mEdges.push_back(mEdgeSeq); }
      unsigned offset = 0;
      bool err = false;
      for ( BigKmer const* pEnt : mEdgeEntries )
      { if ( pEnt->isUnassigned() )
          const_cast<BigKmer*>(pEnt)->setAssigned();
        else
        { std::cout << edgeId << ':' << offset++ << ' ' << pEnt
                    << " Already occupied." << std::endl;
          err = true; } }
        if ( err )
            FatalErr("Having trouble with preoccupied kmers.");
        mEdgeSeq.clear();
        mEdgeEntries.clear(); }

    BigKDict const& mDict;
    vecbvec& mEdges;
    std::vector<BigKmer const*> mEdgeEntries;
    bvec mEdgeSeq;
};

template <unsigned BIGK>
class Pather_sleek_safe
{
public:

    Pather_sleek_safe( vecbvec const& reads, HyperBasevector const& HBV,
                vec<int> const& fwdXlat, vec<int> const& revXlat,
                ReadPathVec* pReadPaths, 
                HyperKmerPath const* pHKP, vecKmerPath* pKmerPaths,
                const vec<int>& starters )
    : mReads(reads), mHBV(HBV), mFwdXlat(fwdXlat), mRevXlat(revXlat), 
      mpReadPaths(pReadPaths), mpHKP(pHKP), mpKmerPaths(pKmerPaths), mStarters(starters)
    {HBV.ToRight(mToRight);}

    void operator()( size_t readId )
    { bvec const& read = mReads[readId];
      if ( read.size() < BIGK ) return;
      if(mStarters[3*readId+1]<0 && mStarters[3*readId]<0) return;
      KMerHasher<BIGK> hasher;
      size_t idx = mStarters[3*readId+1];
      unsigned offset = mStarters[3*readId];
      BigKMer<BIGK> startMer(read,hasher(read.begin()));
      mTmpReadPath.setFirstSkip(offset);
      size_t edgeId = (mStarters[3*readId+2])?mRevXlat[idx]:mFwdXlat[idx];
      mTmpReadPath.push_back(edgeId);
      size_t readLenRemaining = read.size();
      size_t edgeLenRemaining = mHBV.O(edgeId).size()-offset;
      while ( readLenRemaining > edgeLenRemaining )
      { readLenRemaining = readLenRemaining-edgeLenRemaining+BIGK-1; 
        unsigned readOffset = read.size()-readLenRemaining; 
        kmer<BIGK> nextReadMer, edgeLeadMer;
        nextReadMer.SetToSubOf(read,readOffset);
        size_t rv = mToRight[edgeId];
        bool found = false;
        for(size_t b = 0ul; b < mHBV.From(rv).size(); b++){
            edgeId = mHBV.IFrom(rv,b);
            edgeLeadMer.SetToSubOf(mHBV.O(edgeId),0);
            if(nextReadMer == edgeLeadMer){
                found=true;
                break;
            }
        }
        if(!found){ mTmpReadPath.clear(); return;}
        mTmpReadPath.push_back(edgeId); 
        edgeLenRemaining = mHBV.O(edgeId).size(); }
      edgeLenRemaining -= readLenRemaining;
      // must align end to end
      vec<int> rpv;
      for(auto e: mTmpReadPath) rpv.push_back(e); 
      bvec rpseq = mHBV.Cat(rpv), mseq;
      int offs = (mTmpReadPath.getOffset()>=0)? 0 : mTmpReadPath.getOffset();
      int start = Max(0,mTmpReadPath.getOffset());
      int sz = Min(read.size()+offs,rpseq.size()-start);
      mseq.SetToSubOf(rpseq,start,sz);
      if(mseq!=read) { mTmpReadPath.clear(); return;}
      if ( mpReadPaths )
        (*mpReadPaths)[readId] = mTmpReadPath;
      if ( mpKmerPaths )
      { buildKmerPath(edgeLenRemaining);
        (*mpKmerPaths)[readId] = mTmpKmerPath; mTmpKmerPath.clear(); }
      mTmpReadPath.clear(); }

private:
    void buildKmerPath( size_t edgeLenRemaining )
    { if ( mTmpReadPath.size() == 1 )
      { KmerPath const& edge = mpHKP->EdgeObject(mTmpReadPath.front());
        mTmpKmerPath.Assign(edge.front().Start()+mTmpReadPath.getFirstSkip(),
                            edge.back().Stop()-edgeLenRemaining); }
      else
      { KmerPath const& edge1 = mpHKP->EdgeObject(mTmpReadPath.front());
        mTmpKmerPath.Assign(edge1.front().Start()+mTmpReadPath.getFirstSkip(),
                            edge1.back().Stop());
        auto end = mTmpReadPath.end()-1;
        for ( auto itr=mTmpReadPath.begin()+1; itr != end; ++itr )
        { KmerPath const& edge = mpHKP->EdgeObject(*itr);
          SerfVec<KmerPathInterval>& tmp = mTmpKmerPath;
          tmp.append(edge.begin(),edge.end()); }
        KmerPath const& edgeN = mpHKP->EdgeObject(*end);
        kmer_id_t stop = edgeN.back().Stop()-edgeLenRemaining;
        mTmpKmerPath.Append(edgeN.front().Start(),stop); } }

    vecbvec const& mReads;
    HyperBasevector const& mHBV;
    vec<int> const& mFwdXlat;
    vec<int> const& mRevXlat;
    ReadPathVec* mpReadPaths;
    HyperKmerPath const* mpHKP;
    vecKmerPath* mpKmerPaths;
    vec<int> const& mStarters;
    ReadPath mTmpReadPath;
    KmerPath mTmpKmerPath;
    // This assumes the number of vertices is type int
    vec<int> mToRight;
};

template <unsigned BIGK>
class Pather_sleek
{
public:

    Pather_sleek( vecbvec const& reads, HyperBasevector const& HBV,
                vec<int> const& fwdXlat, vec<int> const& revXlat,
                ReadPathVec* pReadPaths, 
                HyperKmerPath const* pHKP, vecKmerPath* pKmerPaths,
                const vec<int>& starters )
    : mReads(reads), mHBV(HBV), mFwdXlat(fwdXlat), mRevXlat(revXlat), 
      mpReadPaths(pReadPaths), mpHKP(pHKP), mpKmerPaths(pKmerPaths), mStarters(starters)
    {HBV.ToRight(mToRight);}

    void seed( size_t readId )
    { bvec const& read = mReads[readId];
      if ( read.size() < BIGK ) return;
      KMerHasher<BIGK> hasher;
      size_t idx = mStarters[3*readId+1];
      unsigned offset = mStarters[3*readId];
      BigKMer<BIGK> startMer(read,hasher(read.begin()));
      mTmpReadPath.setFirstSkip(offset);
      size_t edgeId = (mStarters[3*readId+2])?mRevXlat[idx]:mFwdXlat[idx];
      mTmpReadPath.push_back(edgeId); 
      if ( mpReadPaths )
        (*mpReadPaths)[readId] = mTmpReadPath;
      mTmpReadPath.clear();
    }

    void extend( size_t readId )
    { bvec const& read = mReads[readId];
      if ( read.size() < BIGK ) return;
      mTmpReadPath = (*mpReadPaths)[readId];
      size_t edgeId = mTmpReadPath[0];
      size_t readLenRemaining = read.size();
      size_t edgeLenRemaining = mHBV.O(edgeId).size()-mTmpReadPath.getOffset();
      while ( readLenRemaining > edgeLenRemaining )
      { readLenRemaining = readLenRemaining-edgeLenRemaining+BIGK-1; 
        unsigned readOffset = read.size()-readLenRemaining; 
        kmer<BIGK> nextReadMer, edgeLeadMer;
        nextReadMer.SetToSubOf(read,readOffset);
        size_t rv = mToRight[edgeId];
        bool found = false;
        for(size_t b = 0ul; b < mHBV.From(rv).size(); b++){
            edgeId = mHBV.IFrom(rv,b);
            edgeLeadMer.SetToSubOf(mHBV.O(edgeId),0);
            if(nextReadMer == edgeLeadMer){
                found=true;
                break;
            }
        }
        ForceAssertEq(found,true);
        mTmpReadPath.push_back(edgeId); 
        edgeLenRemaining = mHBV.O(edgeId).size(); }
      edgeLenRemaining -= readLenRemaining;
      if ( mpReadPaths )
        (*mpReadPaths)[readId] = mTmpReadPath;
      if ( mpKmerPaths )
      { buildKmerPath(edgeLenRemaining);
        (*mpKmerPaths)[readId] = mTmpKmerPath; mTmpKmerPath.clear(); }
      mTmpReadPath.clear(); }

    void operator()( size_t readId )
    { bvec const& read = mReads[readId];
      if ( read.size() < BIGK ) return;
      KMerHasher<BIGK> hasher;
      size_t idx = mStarters[3*readId+1];
      unsigned offset = mStarters[3*readId];
      BigKMer<BIGK> startMer(read,hasher(read.begin()));
      mTmpReadPath.setFirstSkip(offset);
      size_t edgeId = (mStarters[3*readId+2])?mRevXlat[idx]:mFwdXlat[idx];
      mTmpReadPath.push_back(edgeId);
      size_t readLenRemaining = read.size();
      size_t edgeLenRemaining = mHBV.O(edgeId).size()-offset;
      while ( readLenRemaining > edgeLenRemaining )
      { readLenRemaining = readLenRemaining-edgeLenRemaining+BIGK-1; 
        unsigned readOffset = read.size()-readLenRemaining; 
        kmer<BIGK> nextReadMer, edgeLeadMer;
        nextReadMer.SetToSubOf(read,readOffset);
        size_t rv = mToRight[edgeId];
        bool found = false;
        for(size_t b = 0ul; b < mHBV.From(rv).size(); b++){
            edgeId = mHBV.IFrom(rv,b);
            edgeLeadMer.SetToSubOf(mHBV.O(edgeId),0);
            if(nextReadMer == edgeLeadMer){
                found=true;
                break;
            }
        }
        ForceAssertEq(found,true);
        mTmpReadPath.push_back(edgeId); 
        edgeLenRemaining = mHBV.O(edgeId).size(); }
      edgeLenRemaining -= readLenRemaining;
      if ( mpReadPaths )
        (*mpReadPaths)[readId] = mTmpReadPath;
      if ( mpKmerPaths )
      { buildKmerPath(edgeLenRemaining);
        (*mpKmerPaths)[readId] = mTmpKmerPath; mTmpKmerPath.clear(); }
      mTmpReadPath.clear(); }

private:
    void buildKmerPath( size_t edgeLenRemaining )
    { if ( mTmpReadPath.size() == 1 )
      { KmerPath const& edge = mpHKP->EdgeObject(mTmpReadPath.front());
        mTmpKmerPath.Assign(edge.front().Start()+mTmpReadPath.getFirstSkip(),
                            edge.back().Stop()-edgeLenRemaining); }
      else
      { KmerPath const& edge1 = mpHKP->EdgeObject(mTmpReadPath.front());
        mTmpKmerPath.Assign(edge1.front().Start()+mTmpReadPath.getFirstSkip(),
                            edge1.back().Stop());
        auto end = mTmpReadPath.end()-1;
        for ( auto itr=mTmpReadPath.begin()+1; itr != end; ++itr )
        { KmerPath const& edge = mpHKP->EdgeObject(*itr);
          SerfVec<KmerPathInterval>& tmp = mTmpKmerPath;
          tmp.append(edge.begin(),edge.end()); }
        KmerPath const& edgeN = mpHKP->EdgeObject(*end);
        kmer_id_t stop = edgeN.back().Stop()-edgeLenRemaining;
        mTmpKmerPath.Append(edgeN.front().Start(),stop); } }

    vecbvec const& mReads;
    HyperBasevector const& mHBV;
    vec<int> const& mFwdXlat;
    vec<int> const& mRevXlat;
    ReadPathVec* mpReadPaths;
    HyperKmerPath const* mpHKP;
    vecKmerPath* mpKmerPaths;
    vec<int> const& mStarters;
    ReadPath mTmpReadPath;
    KmerPath mTmpKmerPath;
    // This assumes the number of vertices is type int
    vec<int> mToRight;
};


template <unsigned BIGK>
class Pather
{
public:
    typedef BigDict<BIGK> BigKDict;

    Pather( vecbvec const& reads, BigKDict const& dict, vecbvec const& edges,
                vec<int> const& fwdXlat, vec<int> const& revXlat,
                ReadPathVec* pReadPaths,
                HyperKmerPath const* pHKP, vecKmerPath* pKmerPaths )
    : mReads(reads), mDict(dict), mEdges(edges),
      mFwdXlat(fwdXlat), mRevXlat(revXlat), mpReadPaths(pReadPaths),
      mpHKP(pHKP), mpKmerPaths(pKmerPaths) {}

    void operator()( size_t readId )
    { bvec const& read = mReads[readId];
      if ( read.size() < BIGK ) return;
      KMerHasher<BIGK> hasher;
      BigKMer<BIGK> kmer(read,hasher(read.begin()));
      BigKMer<BIGK> entry = lookup(kmer);
      mTmpReadPath.setFirstSkip(entry.getOffset());
      size_t idx = &entry.getBV()-&mEdges[0];
      size_t edgeId = entry.isRC()?mRevXlat[idx]:mFwdXlat[idx];
      mTmpReadPath.push_back(edgeId);
      size_t readLenRemaining = read.size();
      size_t edgeLenRemaining = entry.getBV().size()-entry.getOffset();
      while ( readLenRemaining > edgeLenRemaining )
      { readLenRemaining = readLenRemaining-edgeLenRemaining+BIGK-1;
        unsigned readOffset = read.size()-readLenRemaining;
        BigKMer<BIGK> nextKmer(read,hasher(read.begin(readOffset)),
                                KMerContext(),readOffset);
        BigKMer<BIGK> nextEntry = lookup(nextKmer);
        ForceAssertEq(nextEntry.getOffset(),0u);
        idx = &nextEntry.getBV()-&mEdges[0];
        edgeId = nextEntry.isRC()?mRevXlat[idx]:mFwdXlat[idx];
        mTmpReadPath.push_back(edgeId);
        edgeLenRemaining = nextEntry.getBV().size(); }
      edgeLenRemaining -= readLenRemaining;
      /* mTmpReadPath.setLastSkip(edgeLenRemaining); */
      if ( mpReadPaths )
        (*mpReadPaths)[readId] = mTmpReadPath;
      if ( mpKmerPaths )
      { buildKmerPath(edgeLenRemaining);
        (*mpKmerPaths)[readId] = mTmpKmerPath; mTmpKmerPath.clear(); }
      mTmpReadPath.clear(); }

private:
    BigKMer<BIGK> lookup( BigKMer<BIGK> const& kmer )
    { if ( kmer.isRev() )
      { BigKMer<BIGK> const* pEnt = mDict.lookup(kmer.rc());
        ForceAssert(pEnt);
        return pEnt->rc(); }
      BigKMer<BIGK> const* pEnt = mDict.lookup(kmer);
      ForceAssert(pEnt);
      return *pEnt; }

    void buildKmerPath( size_t edgeLenRemaining )
    { if ( mTmpReadPath.size() == 1 )
      { KmerPath const& edge = mpHKP->EdgeObject(mTmpReadPath.front());
        mTmpKmerPath.Assign(edge.front().Start()+mTmpReadPath.getFirstSkip(),
                            edge.back().Stop()-edgeLenRemaining); }
      else
      { KmerPath const& edge1 = mpHKP->EdgeObject(mTmpReadPath.front());
        mTmpKmerPath.Assign(edge1.front().Start()+mTmpReadPath.getFirstSkip(),
                            edge1.back().Stop());
        auto end = mTmpReadPath.end()-1;
        for ( auto itr=mTmpReadPath.begin()+1; itr != end; ++itr )
        { KmerPath const& edge = mpHKP->EdgeObject(*itr);
          SerfVec<KmerPathInterval>& tmp = mTmpKmerPath;
          tmp.append(edge.begin(),edge.end()); }
        KmerPath const& edgeN = mpHKP->EdgeObject(*end);
        kmer_id_t stop = edgeN.back().Stop()-edgeLenRemaining;
        mTmpKmerPath.Append(edgeN.front().Start(),stop); } }

    vecbvec const& mReads;
    BigKDict const& mDict;
    vecbvec const& mEdges;
    vec<int> const& mFwdXlat;
    vec<int> const& mRevXlat;
    ReadPathVec* mpReadPaths;
    HyperKmerPath const* mpHKP;
    vecKmerPath* mpKmerPaths;
    ReadPath mTmpReadPath;
    KmerPath mTmpKmerPath;
};

#if 0
void ValidatePaths( vecbvec const& reads, ReadPathVec const& readPaths,
                        HyperBasevector const& hbv )
{
    auto iPath = readPaths.begin();
    for ( bvec const& read : reads )
    {
        ReadPath const& path = *iPath;
        ++iPath;
        if ( path.empty() )
        {
            std::cout << "Empty path.\n";
            continue;
        }
        auto eBeg = hbv.EdgeObject(path.front()).begin()+path.getFirstSkip();
        auto eEnd = hbv.EdgeObject(path.front()).end();
        auto rBeg = read.begin();
        auto rEnd = read.end();
        size_t dist = std::min(eEnd-eBeg,rEnd-rBeg);
        if ( !std::equal(eBeg,eBeg+dist,rBeg) )
        {
            std::cout << "Bad path compare.\n";
            continue;
        }
        if ( eBeg+dist == eEnd )
        {
            rBeg = rBeg+dist-hbv.K()+1;
        }
        else if ( path.size() > 1 )
        {
            std::cout << "Extra path info.\n";
            continue;
        }
        for ( auto itr=path.begin()+1,end=path.end(); itr != end; ++itr )
        {
            eBeg = hbv.EdgeObject(*itr).begin();
            eEnd = hbv.EdgeObject(*itr).end();
            dist = std::min(eEnd-eBeg,rEnd-rBeg);
            if ( !std::equal(eBeg,eBeg+dist,rBeg) )
            {
                std::cout << "Bad path compare.\n";
                break;
            }
            if ( eBeg+dist == eEnd )
                rBeg = rBeg+dist-hbv.K()+1;
            else if ( itr+1 != end )
            {
                std::cout << "Extra path info.\n";
                break;
            }
        }
    }
}
#endif

template <unsigned BIGK>
void readsToHBV_Perturb( vecbvec const& reads, unsigned coverage,
                    HyperBasevector* pHBV, HyperBasevectorX& gHBX, 
                    digraphE<vec<int>>& D, ReadPathVecX& pathsx,
                    vec<HyperBasevector>& hbms, vec<digraphE<vec<int>>>& Dms,
                    String pi_file, vecbvec& bases)
{
    size_t nKmers = reads.getKmerCount(BIGK);
    if ( !nKmers ){
        pHBV->Clear();
        return;
    }

    // declare stuff to protect from scoping
    vecbvec edges;
    vec<int> fwdXlat;
    vec<int> revXlat;
    vec<int> startersCR;
    vec<int> startersUR;
    vec<int> startersCS;
    vec<int> startersCCE;
    vec<int> startersSG;
    vec<int> startersMA;
    vec<size_type> affected_old_edges;
    vec<int64_t> changed_readpaths;
    vec<int64_t> unpathed_readpaths;
    vec<int64_t> changed_superedges;
    vec<pair<int64_t,int64_t>> seqgap_edges;
    vec<pair<int64_t,int>> changed_cells_edges;
    vecbvec CRseq,URseq,CSseq,CCEseq,SGseq,MAseq; 
    ReadPathVec* pCR_RP = nullptr;
    ReadPathVec* pUR_RP = nullptr;
    ReadPathVec* pCS_RP = nullptr;
    ReadPathVec* pCCE_RP = nullptr;
    ReadPathVec* pSG_RP = nullptr;
    ReadPathVec* pMA_RP = nullptr;
    vec<int> readbaselen;
    vec<triple<int64_t,int64_t,Bool>> unaffected_map; // dest, src
    int64_t numNewE = 0;
    Bool DO_UNPATHS = True;
    Bool DO_CELLS = True;
    Bool DO_SEQGAPS = True;
    Bool DO_SUPEREDGES = True;
    Bool DO_READPATHS = True;
    Bool DO_MICROASM = True;

    // this saves memory
    if(DO_UNPATHS)
    {
        PRINTDEETS("get unpathed reads");
        // populate unpathed rps;
        for(int64_t i = 0; i < pathsx.size(); i++)
            if(pathsx.getNumEdges(i)==0)
                unpathed_readpaths.push_back(i);
        URseq.reserve(unpathed_readpaths.size());
        for(int64_t i = 0; i < unpathed_readpaths.size(); i++)
            URseq.push_back(bases[unpathed_readpaths[i]]);
        readbaselen.resize(bases.size());
        for(int64_t i = 0; i < bases.size(); i++)
            readbaselen[i] = bases[i].size();
        Destroy(bases);
    }

    // build edges and hbv
    {
        PRINTDEETS("build big dict");

        BigDict<BIGK> bigDict(nKmers/coverage);

        PRINTDEETS("before kmerizer");

        // kmerize input
        BigKMerizer<BIGK> kmerizer(&bigDict);
        parallelFor(0ul,reads.size(),
                [kmerizer,&reads]( size_t readId ) mutable
                { kmerizer.kmerize(reads[readId]); });

        PRINTDEETS("build edges after kmerizer");

        // build edges
        edges.reserve(bigDict.size()/100);
        BigKEdgeBuilder<BIGK>::buildEdges(bigDict,&edges);
        AssertEq(edges.getKmerCount(BIGK),bigDict.size());

        PRINTDEETS("update dict from new edges");
        // update dict from edges
        parallelForBatch(0ul,edges.size(),100,
                [&edges,&kmerizer]( size_t edgeId )
                { kmerizer.updateDict(edges[edgeId]); });

        {
            PRINTDEETS("find changed old edges");
            auto oldEdges = gHBX.Edges();
            unaffected_map.resize(gHBX.E());
            for(int64_t i = 0; i < edges.size(); i++)
                numNewE += (edges[i]==ReverseComplement(edges[i]))?1:2;

            parallelForBatch(0ul,oldEdges.size(),100,
                [&edges,&oldEdges,&unaffected_map,&kmerizer]( size_t edgeId )
                { kmerizer.findEdgeIdx(edges,oldEdges,edgeId,unaffected_map); });
            
            vec<Bool> todel(gHBX.E(),False);
            for(int64_t i = 0; i < gHBX.E(); i++){
                if(unaffected_map[i].second<0){
                    affected_old_edges.push_back(i);
                    todel[i] = True;
                }
                if(i>=numNewE){
                    affected_old_edges.push_back(i);
                    todel[i] = True;
                }
            }
            EraseIf(unaffected_map,todel);
            ParallelUniqueSort(affected_old_edges);
            ForceAssertEq(unaffected_map.size()+affected_old_edges.size(),gHBX.E());

            cout << "Fraction of old hbv edges that will change = "
                << affected_old_edges.size() << "/" << gHBX.E() << endl;
            cout << "New #E = " << numNewE << " vs Old #E = " << gHBX.E() << endl;
        }
 
        if(DO_READPATHS)
        {
            PRINTDEETS("find changed reads and store their baseseq");
            // populate old affected readpaths;
            SubsetMasterVec<ULongVec> paths_index_sub(pi_file,affected_old_edges);
            for(int64_t i = 0; i < affected_old_edges.size(); i++)
                for(auto r : paths_index_sub[affected_old_edges[i]])
                    changed_readpaths.push_back(r);
            ParallelUniqueSort(changed_readpaths);
            cout << "Fraction of readpaths that will change = " 
                << changed_readpaths.size() << "/" << pathsx.size() << endl;

            // build affected sequences
            CRseq.reserve(changed_readpaths.size());
            int64_t batch = 1000;
            #pragma omp parallel for ordered schedule(dynamic,1)
            for(int64_t start = 0; start < changed_readpaths.size(); start += batch){
                int64_t stop = Min(start+batch,(int64_t)changed_readpaths.size());
                vecbvec subseq;
                subseq.reserve(stop-start);
                ReadPath path; bvec seq, seq1; vec<int> pp;
                for(int64_t i = start; i < stop; i++){
                    pathsx.unzip(path,gHBX,changed_readpaths[i]);
                    pp.clear(); 
                    for(auto e: path) 
                        pp.push_back(e); 
                    seq = gHBX.Cat(pp);
                    // this min-max ensures the actual alignment is considered and only repaths it
                    // this requires updating the offset when it was originally negative
                    int start = Max(0,path.getOffset());
                    int offs = (path.getOffset()>=0)? 0 : path.getOffset();
                    int sz = Min(readbaselen[changed_readpaths[i]]+offs,seq.isize()-start);
                    seq1.SetToSubOf(seq,start,sz);
                    subseq.push_back(seq1);
                }
                #pragma omp ordered
                CRseq.append(subseq.begin(),subseq.end());
            }
            Destroy(readbaselen);
        }

        if(DO_SUPEREDGES)
        { 
            PRINTDEETS("find changed superedges and store their baseseq");
            // populate affected superdges
            vec<Bool> AOE(gHBX.E(),False);
            for(auto aoe: affected_old_edges)
                AOE[aoe] = True;
            for(int d = 0; d < D.E(); d++)
                if(D.O(d)[0]>=0)
                    for(auto e: D.O(d))
                        if(AOE[e])
                            changed_superedges.push_back(d);
            ParallelUniqueSort(changed_superedges);
            cout << "Fraction of superedges that will change = " 
                << changed_superedges.size() << "/" << D.E() << endl; 
            
            // build affected sequences
            CSseq.reserve(changed_superedges.size());
            for(int64_t i = 0; i < changed_superedges.size(); i++)
                CSseq.push_back(gHBX.Cat(D.O(changed_superedges[i])));

        }

        if(DO_CELLS)
        {
            PRINTDEETS("find changed cell gap edges and store their baseseq");
            vec<Bool> AOE(gHBX.E(),False);
            for(auto aoe: affected_old_edges)
                AOE[aoe] = True;
            vec<Bool> todel;
            // Get affected cells and their changed edges
            vec<vec<int>> changed_ce;
            vec<int> changed_cells;
            for(int64_t i = 0; i < D.E(); i++){
                if(D.O(i)[0]>=0) continue;
                if(!IsCell(D.O(i))) continue;
                changed_cells.push_back(i);
                todel.push_back(False);
                changed_ce.resize(changed_cells.size());
                cell c;
                c.CellDecode(D.O(i));
                const digraphE<vec<int>>& cD = c.G();
                for(int d = 0; d < cD.E(); d++)
                    if(cD.O(d)[0]>=0)
                        for(auto e: cD.O(d))
                            if(AOE[e])
                                changed_ce.back().push_back(d);
                if(changed_ce.back().size()==0)
                    todel.back() = True;
            }
            int numCells = changed_cells.size();
            EraseIf(changed_cells,todel);
            EraseIf(changed_ce,todel);
            cout << "Fraction of cell gaps that will change = "
                << changed_cells.size() << "/" << numCells << endl;

            for(int64_t i = 0; i < changed_ce.size(); i++)
                for(int j = 0; j < changed_ce[i].isize(); j++)
                    changed_cells_edges.push(changed_cells[i],changed_ce[i][j]);
            Destroy(changed_cells);
            Destroy(changed_ce);
            ParallelUniqueSort(changed_cells_edges);

            //build affected seq
            CCEseq.reserve(changed_cells_edges.size());
            cell c;
            for(int64_t i = 0; i < changed_cells_edges.size(); i++){
                c.CellDecode(D.O(changed_cells_edges[i].first));
                CCEseq.push_back(gHBX.Cat(c.G().O(changed_cells_edges[i].second)));
            }
        }
        
        if(DO_SEQGAPS)
        {
            PRINTDEETS("find seq gaps and store their seq");
            vec<int> TL, TR;
            D.ToLeft(TL);
            D.ToRight(TR);
            const int K = gHBX.K();
            #pragma omp parallel for
            for(int i = 0; i < D.E(); i++){
                if(!IsSequence(D.O(i))) continue;
                int v = TL[i], w = TR[i];
                if (D.To(v).size()!=1 || D.From(w).size()!=1){
                    cout << "incorrect sequence gap topology" << endl;
                    Scram(0);
                }
                int dleft, dright, ltrim, rtrim; bvec seqgap;
                dleft = D.ITo(v,0); dright = D.IFrom(w,0);
                GapToSeq(D.O(i), ltrim, rtrim, seqgap); 
                // append extra base to left and right of seqgap
                bvec tmp;
                bvec lseq = gHBX.Cat(D.O(dleft));
                bvec rseq = gHBX.Cat(D.O(dright));
                tmp.push_back(lseq[lseq.size()-ltrim-K]);
                tmp.append(seqgap);
                tmp.push_back(rseq[rtrim+K-1]);
                #pragma omp critical
                { 
                    SGseq.push_back(tmp);
                    seqgap_edges.push_back(make_pair(i,-1));
                }
            }
            for(int64_t i = 0; i < D.E(); i++){
                if(!IsCell(D.O(i))) continue;
                cell c;
                c.CellDecode(D.O(i));
                const digraphE<vec<int>>& cD = c.G();
                vec<int> tl,tr; cD.ToLeft(tl); cD.ToRight(tr);
                for(int d = 0; d < cD.E(); d++){
                    if(IsSequence(cD.O(d))){
                        int dleft, dright, ltrim, rtrim; bvec seqgap;
                        dleft = cD.ITo(tl[d],0); dright = cD.IFrom(tr[d],0);
                        GapToSeq(cD.O(d), ltrim, rtrim, seqgap); 
                        // append extra base to left and right of seqgap
                        bvec tmp;
                        bvec lseq = gHBX.Cat(cD.O(dleft));
                        bvec rseq = gHBX.Cat(cD.O(dright));
                        tmp.push_back(lseq[lseq.size()-ltrim-K]);
                        tmp.append(seqgap);
                        tmp.push_back(rseq[rtrim+K-1]);
                        SGseq.push_back(tmp);
                        seqgap_edges.push_back(make_pair(i,d));
                    }
                }
            }
            ParallelSortSync(seqgap_edges, SGseq);
            cout << "Fraction of superedges that are sequence gaps = " 
                << seqgap_edges.size() << "/" << D.E() << endl;
        }

        if(DO_MICROASM){
            PRINTDEETS("find microasm edges and store their seq");
            for(int m = 0 ; m < hbms.isize(); m++){
                for(int d = 0; d < hbms[m].E(); d++){
                    MAseq.push_back(hbms[m].O(d));
                }
            }
            cout << "Number of microassembly edges = " << MAseq.size() << endl;
        }
 
        if(DO_READPATHS){
            PRINTDEETS("get starting kmers of changed reads");
            // get start kmers of input selection and find input edges they land on
            startersCR.resize(3*CRseq.size(),0); // store [offset,idx,isRC] per CR
            parallelFor(0ul,CRseq.size(),
                    [kmerizer,&CRseq,&startersCR,&edges]( size_t readId ) mutable
                    {kmerizer.kmerize_getLanding(CRseq[readId],&startersCR[3*readId],edges); });
        }

        if(DO_SUPEREDGES){
            PRINTDEETS("get starting kmers of changed superedges");
            startersCS.resize(3*CSseq.size(),0); // store [offset,idx,isRC] per CS
            parallelFor(0ul,CSseq.size(),
                    [kmerizer,&CSseq,&startersCS,&edges]( size_t readId ) mutable
                    {kmerizer.kmerize_getLanding(CSseq[readId],&startersCS[3*readId],edges); });
        }

        if(DO_CELLS){
            PRINTDEETS("get starting kmers of changed cell gap edges");
            startersCCE.resize(3*CCEseq.size(),0); // store [offset,idx,isRC] per CCE
            parallelFor(0ul,CCEseq.size(),
                    [kmerizer,&CCEseq,&startersCCE,&edges]( size_t readId ) mutable
                    {kmerizer.kmerize_getLanding(CCEseq[readId],&startersCCE[3*readId],edges); });
        }

        if(DO_UNPATHS){
            PRINTDEETS("get starting kmers of unpathed edges");
            startersUR.resize(3*URseq.size(),0); // store [offset,idx,isRC] per UR
            parallelFor(0ul,URseq.size(),
                    [kmerizer,&URseq,&startersUR,&edges]( size_t readId ) mutable
                    {kmerizer.kmerize_getLanding_safe(URseq[readId],&startersUR[3*readId],edges); });
        }
    
        if(DO_SEQGAPS){
            PRINTDEETS("get starting kmers of seqgap edges");
            startersSG.resize(3*SGseq.size(),0); // store [offset,idx,isRC] per SG
            parallelFor(0ul,SGseq.size(),
                    [kmerizer,&SGseq,&startersSG,&edges]( size_t readId ) mutable
                    {kmerizer.kmerize_getLanding(SGseq[readId],&startersSG[3*readId],edges); });
        }

        if(DO_MICROASM){
            PRINTDEETS("get starting kmers of microasm edges");
            startersMA.resize(3*MAseq.size(),0); // store [offset,idx,isRC] per SG
            parallelFor(0ul,MAseq.size(),
                    [kmerizer,&MAseq,&startersMA,&edges]( size_t readId ) mutable
                    {kmerizer.kmerize_getLanding(MAseq[readId],&startersMA[3*readId],edges); });
        }

        // destroy
        bigDict.clear();
    }
    
    {
        PRINTDEETS("build new graph");
        // build new graph
        buildHBVFromEdges(edges,BIGK,pHBV,&fwdXlat,&revXlat);
        ForceAssertEq(numNewE,pHBV->E());
        Destroy(edges);

        PRINTDEETS("reorder edges");
        // reorder edges and find affected edges
        {
            for(int64_t i = 0; i < unaffected_map.size(); i++)
                    unaffected_map[i].second = (unaffected_map[i].third) ? 
                                        revXlat[unaffected_map[i].second] : 
                                        fwdXlat[unaffected_map[i].second];

            // get permutation: src idx -> dest idx
            vec<int64_t> permutation(pHBV->E());
            vec<Bool> used_src(pHBV->E(),False);
            vec<Bool> used_dest(pHBV->E(),False);

            #pragma omp parallel for
            for(int64_t i = 0; i < unaffected_map.size(); i++){
                permutation[unaffected_map[i].second] = unaffected_map[i].first; // source -> dest
                used_dest[unaffected_map[i].first] = True; // mark used dest
                used_src[unaffected_map[i].second] = True; // mark used src
            }
            Destroy(unaffected_map);

            vec<int64_t> unused_src, unused_dest;
            for(int64_t i = 0; i < used_src.size(); i++)
                if(!used_src[i])
                    unused_src.push_back(i);
            for(int64_t i = 0; i < used_dest.size(); i++)
                if(!used_dest[i])
                    unused_dest.push_back(i);
            ParallelSort(unused_src);
            ParallelSort(unused_dest);
            ForceAssertEq(unused_src.size(),unused_dest.size());

            #pragma omp parallel for
            for(int64_t i = 0; i < unused_src.size(); i++)
                permutation[unused_src[i]] = unused_dest[i];
            Destroy(used_src);
            Destroy(used_dest);
            Destroy(unused_src);
            Destroy(unused_dest);

            // modify hb
            pHBV->ReorderEdges(permutation); 

            // modify fwdXlat, revXlat
            #pragma omp parallel for
            for(int64_t i = 0; i < fwdXlat.size(); i++){
                fwdXlat[i] = permutation[fwdXlat[i]];
                revXlat[i] = permutation[revXlat[i]];
            }
        }

        // do pathing of affected seqs
        HyperKmerPath* pHKP = nullptr;
        vecKmerPath* pKmerPaths = nullptr;
        HyperBasevector& spHBV = *pHBV;

        if(DO_READPATHS)
        {
            PRINTDEETS("repath changed reads");
            pCR_RP = new ReadPathVec;
            pCR_RP->resize(CRseq.size());
            Pather_sleek<BIGK> pather_slk(CRseq,spHBV,fwdXlat,revXlat,
                                            pCR_RP,pHKP,pKmerPaths,startersCR);
            parallelForBatch(0ul,CRseq.size(),1000,pather_slk);
            /* cout << "CR_RP checksum: " << checksumRP(*pCR_RP) << endl; */
            Destroy(startersCR);
            Destroy(CRseq);
        }
        
        if(DO_SUPEREDGES)
        {
            PRINTDEETS("repath changed superedges");
            pCS_RP = new ReadPathVec;
            pCS_RP->resize(CSseq.size());
            Pather_sleek<BIGK> pather_slk(CSseq,spHBV,fwdXlat,revXlat,
                                            pCS_RP,pHKP,pKmerPaths,startersCS);
            parallelForBatch(0ul,CSseq.size(),1000,pather_slk);
            /* cout << "CS_RP checksum: " << checksumRP(*pCS_RP) << endl; */
            Destroy(startersCS);
            Destroy(CSseq);
        }

        if(DO_CELLS)
        {
            PRINTDEETS("repath changed cell gap edges");
            pCCE_RP = new ReadPathVec;
            pCCE_RP->resize(CCEseq.size());
            Pather_sleek<BIGK> pather_slk(CCEseq,spHBV,fwdXlat,revXlat,
                                            pCCE_RP,pHKP,pKmerPaths,startersCCE);
            parallelForBatch(0ul,CCEseq.size(),1000,pather_slk);
            /* cout << "CCE_RP checksum: " << checksumRP(*pCCE_RP) << endl; */
            Destroy(startersCCE);
            Destroy(CCEseq);
        }

        if(DO_UNPATHS)
        {
            PRINTDEETS("repath unpathed reads");
            pUR_RP = new ReadPathVec;
            pUR_RP->resize(URseq.size());
            Pather_sleek_safe<BIGK> pather_slk_safe(URseq,spHBV,fwdXlat,revXlat,
                                            pUR_RP,pHKP,pKmerPaths,startersUR);
            parallelForBatch(0ul,URseq.size(),1000,pather_slk_safe);
            /* cout << "UR_RP checksum: " << checksumRP(*pUR_RP) << endl; */
            Destroy(startersUR);
            Destroy(URseq);
        }

        if(DO_SEQGAPS)
        {
            PRINTDEETS("repath seqgap edges");
            pSG_RP = new ReadPathVec;
            pSG_RP->resize(SGseq.size());
            Pather_sleek<BIGK> pather_slk(SGseq,spHBV,fwdXlat,revXlat,
                                            pSG_RP,pHKP,pKmerPaths,startersSG);
            parallelForBatch(0ul,SGseq.size(),1000,pather_slk);
            /* cout << "SG_RP checksum: " << checksumRP(*pSG_RP) << endl; */
            Destroy(startersSG);
            Destroy(SGseq);
        }

        if(DO_MICROASM)
        {
            PRINTDEETS("repath micro asm edges");
            pMA_RP = new ReadPathVec;
            pMA_RP->resize(MAseq.size());
            Pather_sleek<BIGK> pather_slk(MAseq,spHBV,fwdXlat,revXlat,
                                            pMA_RP,pHKP,pKmerPaths,startersMA);
            parallelForBatch(0ul,MAseq.size(),1000,pather_slk);
            /* cout << "MA_RP checksum: " << checksumRP(*pMA_RP) << endl; */
            Destroy(startersMA);
            Destroy(MAseq);
        }

        Destroy(fwdXlat);
        Destroy(revXlat);
   
        if(DO_READPATHS || DO_UNPATHS)
        {
            PRINTDEETS("update readpaths");
            cout << "Fraction of readpaths that changed = " 
                << pCR_RP->size() << "/" << pathsx.size() << endl;
            //!!! Note expensive conversion
            HyperBasevectorX mHBX = HyperBasevectorX(*pHBV);
            // modify all readpaths in with new encryption/ new readpath
            int batch = 5000;
            ReadPathVecX npaths; npaths.reserve(pathsx.size());
            #pragma omp parallel for ordered schedule(dynamic,1)
            for(int64_t start = 0; start < pathsx.size(); start += batch){
                int64_t stop = Min(start+batch,pathsx.size());
                ReadPath path;
                ReadPathVecX subRPVX;
                subRPVX.reserve(stop-start);
                for(int64_t i = start; i < stop; i++){
                    pathsx.unzip(path, gHBX, i);
                    int offs = path.getOffset();
                    int pos;
                    if(DO_READPATHS){
                        pos = BinPosition(changed_readpaths,i);
                        if(pos>=0){
                            path = (*pCR_RP)[pos];
                            // also need to update the offsets when partially aligned
                            if(offs<0)
                                path.setOffset(offs + path.getOffset());
                        }
                    }
                    if(DO_UNPATHS){
                        pos = BinPosition(unpathed_readpaths,i);
                        if(pos>=0)
                            path = (*pUR_RP)[pos];
                    }
                    subRPVX.append(path,mHBX);
                }
                #pragma omp ordered
                npaths.append(subRPVX);
            }

            int64_t saved = 0;
            for(int64_t i = 0; i < unpathed_readpaths.size(); i++)
                if((*pUR_RP)[i].size()) saved++;
            cout << "Fraction of unpathed readpaths that were saved = " 
                << saved << "/" << pUR_RP->size() << endl;

            pathsx = npaths;
            Destroy(changed_readpaths);
            Destroy(unpathed_readpaths);
            Destroy(*pCR_RP);
            Destroy(*pUR_RP);
            pCR_RP=nullptr;
            pUR_RP=nullptr;
        }

        if(DO_SUPEREDGES)
        {   
            PRINTDEETS("update supergraph");
            cout << "Fraction of superedges that changed = " << pCS_RP->size() << "/" << D.E() << endl;
            // update D
            int offext = 0;
            #pragma omp parallel for
            for(int i = 0; i < changed_superedges.isize(); i++){
                vec<int> newD;
                for(auto e : (*pCS_RP)[i])
                    newD.push_back(e);
                ForceAssertGt(newD.isize(),0);
                if((*pCS_RP)[i].getOffset()!=0){
                    #pragma omp atomic
                    offext++;
                }
                D.OMutable(changed_superedges[i]) = newD;
            }
            cout << "number of remade superedges with nnz offsets on base-assembly = " 
                << offext << "/" << changed_superedges.size() << endl;
            Destroy(changed_superedges);
            Destroy(*pCS_RP);
            pCS_RP=nullptr;
        }

        if(DO_CELLS)
        {
            PRINTDEETS("update cell gaps");
            vec<int64_t> splits; splits.push_back(0);
            for(int64_t i = 1; i < changed_cells_edges.size(); i++)
                if(changed_cells_edges[i].first != changed_cells_edges[i-1].first)
                    splits.push_back(i);
            splits.push_back(changed_cells_edges.size()); 
            #pragma omp parallel for
            for(int s = 0; s < splits.size()-1; s++){
                for(int64_t i = splits[s]; i < splits[s+1]; i++){
                    cell c;
                    c.CellDecode(D.O(changed_cells_edges[i].first));
                    vec<int> pp;
                    for(auto e : (*pCCE_RP)[i]) pp.push_back(e);
                    c.GMutable().OMutable(changed_cells_edges[i].second) = pp;
                    D.OMutable(changed_cells_edges[i].first).clear();
                    c.CellEncode(D.OMutable(changed_cells_edges[i].first));
                }
            }
            Destroy(changed_cells_edges);
            Destroy(*pCCE_RP);
            pCCE_RP=nullptr;
        }

        if(DO_SEQGAPS)
        {
            PRINTDEETS("update seq gaps");
            vec<int> TL, TR;
            D.ToLeft(TL);
            D.ToRight(TR);
            std::unordered_map<int,int> vpair;
            int punt = 0;
            for(int64_t i = 0; i < seqgap_edges.size(); i++){
                if(seqgap_edges[i].second>-1) continue; // normal seqgaps first
                int64_t sg = seqgap_edges[i].first;
                ForceAssertEq(IsSequence(D.O(sg)),True);
                // process each gap
                vec<int> pp;
                for(auto e : (*pSG_RP)[i]) pp.push_back(e);
                ForceAssertGt(pp.isize(),0);
                if(vpair.count(TL[sg])==0){
                    vpair[TL[sg]] = TR[sg];
                    vec<int>& ld = D.OMutable(D.ITo(TL[sg],0));
                    vec<int>& rd = D.OMutable(D.IFrom(TR[sg],0));
                    int ldt = ld.isize()-1, rdt = 0;
                    if(pp.size())
                        while(pp[0]!=ld[ldt] && ldt >= 0) ldt--;
                    ForceAssertGe(ldt,0);
                    if(pp.size())
                        while(pp.back()!=rd[rdt] && rdt < rd.isize()) rdt++;
                    ForceAssertLt(rdt,rd.isize());
                    ld.resize(ldt+1);
                    for(int j = 0; j < rd.size()-rdt; j++)
                        rd[j] = rd[j+rdt];
                    rd.resize(rd.size()-rdt);
                }
                vec<int> qq;
                for(int j = 1; j < pp.isize()-1; j++)
                    qq.push_back(pp[j]);
                D.OMutable(sg).clear();
                if(qq.size()!=0)
                    D.OMutable(sg) = qq;
                else{
                    D.OMutable(sg).push_back(-1); // mark as pair gap
                    punt++;
                }
            }
            for(int64_t i = 0; i < seqgap_edges.size(); i++){
                if(seqgap_edges[i].second<0) continue; // seqs in cells
                vpair.clear();
                int64_t cd = seqgap_edges[i].first;
                int64_t sg = seqgap_edges[i].second;
                cell c;
                c.CellDecode(D.O(cd));
                auto& cD = c.GMutable();
                vec<int> tr, tl;
                cD.ToLeft(tl); cD.ToRight(tr);
                ForceAssertEq(IsSequence(cD.O(sg)),True);
                // process each gap
                vec<int> pp;
                for(auto e : (*pSG_RP)[i]) pp.push_back(e);
                ForceAssertGt(pp.isize(),0);
                if(vpair.count(tl[sg])==0){
                    vpair[tl[sg]] = tr[sg];
                    vec<int>& ld = cD.OMutable(cD.ITo(tl[sg],0));
                    vec<int>& rd = cD.OMutable(cD.IFrom(tr[sg],0));
                    int ldt = ld.isize()-1, rdt = 0;
                    if(pp.size())
                        while(pp[0]!=ld[ldt] && ldt >= 0) ldt--;
                    ForceAssertGe(ldt,0);
                    if(pp.size())
                        while(pp.back()!=rd[rdt] && rdt < rd.isize()) rdt++;
                    ForceAssertLt(rdt,rd.isize());
                    ld.resize(ldt+1);
                    for(int j = 0; j < rd.size()-rdt; j++)
                        rd[j] = rd[j+rdt];
                    rd.resize(rd.size()-rdt);
                }
                vec<int> qq;
                for(int j = 1; j < pp.isize()-1; j++)
                    qq.push_back(pp[j]);
                D.OMutable(cd).clear();
                if(qq.size()!=0)
                    cD.OMutable(sg) = qq;
                else{
                    cD.OMutable(sg).push_back(-1);
                    punt++;
                }
                c.CellEncode(D.OMutable(cd));
            }
            cout << "punted on " << punt << "/" << seqgap_edges.size() << " seq-gaps" << endl;
            Destroy(seqgap_edges);
            Destroy(*pSG_RP);
            pSG_RP = nullptr;
        }

        if(DO_MICROASM)
        {
            PRINTDEETS("update microasm edges");
            int i = 0;
            Dms.clear();
            Dms.resize(hbms.size());
            for(int m = 0 ; m < hbms.isize(); m++){
                vec<vec<int>> newEdges(hbms[m].E());
                for(int d = 0; d < hbms[m].E(); d++){
                    vec<int>& newD = newEdges[d];
                    if(hbms[m].O(d).size()==0)
                       newD.push_back(-1);  
                    else
                        for(auto e : (*pMA_RP)[i])
                            newD.push_back(e);
                    i++;
                }
                Dms[m].Initialize(hbms[m],newEdges);
            }
        }
    }
}

template <unsigned BIGK>
void readsToHBV_Sleek( vecbvec const& reads, unsigned coverage,
                    HyperBasevector* pHBV, ReadPathVec* pReadPaths,
                    HyperKmerPath* pHKP, vecKmerPath* pKmerPaths)
{
    
    size_t nKmers = reads.getKmerCount(BIGK);
    if ( !nKmers )
    {
        // return empty reads
        if ( pReadPaths )
        {
            pReadPaths->clear();
            pReadPaths->resize(reads.size());
        }
        if ( pKmerPaths )
        {
            pKmerPaths->clear();
            pKmerPaths->resize(reads.size());
        }
        pHBV->Clear();
        if ( pHKP ) pHKP->Clear();
        return;
    }

    // declare stuff to protect from scoping
    vecbvec edges;
    vec<int> fwdXlat;
    vec<int> revXlat;
    vec<int> starters;

    // build edges and hbv
    {
        MEM(before_big_dict);

        BigDict<BIGK> bigDict(nKmers/coverage);

        MEM(before_big_kmerizer);

        BigKMerizer<BIGK> kmerizer(&bigDict);

        parallelFor(0ul,reads.size(),
                [kmerizer,&reads]( size_t readId ) mutable
                { kmerizer.kmerize(reads[readId]); });

        MEM(after_kmerizer);

        edges.reserve(bigDict.size()/100);
        BigKEdgeBuilder<BIGK>::buildEdges(bigDict,&edges);
        AssertEq(edges.getKmerCount(BIGK),bigDict.size());

        MEM(after_edges);

        parallelForBatch(0ul,edges.size(),100,
                [&edges,&kmerizer]( size_t edgeId )
                { kmerizer.updateDict(edges[edgeId]); });

        if (!( !pReadPaths && !pKmerPaths ))
        {
            // get start kmers of input reads and find input edges they land on
            starters.resize(3*reads.size(),0); // store [offset,idx,isRC] per read
            parallelFor(0ul,reads.size(),
                    [kmerizer,&reads,&starters,&edges]( size_t readId ) mutable
                    {kmerizer.kmerize_getLanding(reads[readId],&starters[3*readId],edges); });
        }
    
        // destroy
        bigDict.clear();

    }
    
    if ( !pReadPaths && !pKmerPaths )
    {
        buildHBVFromEdges(edges,BIGK,pHBV,&fwdXlat,&revXlat);
        if ( pHKP )
            buildHKPFromHBV(*pHBV,fwdXlat,revXlat,pHKP);
        MEM(after_build);
    }
    else{

        buildHBVFromEdges(edges,BIGK,pHBV,&fwdXlat,&revXlat);

        MEM(after_build2);

        HyperKmerPath hkp;
        if ( pKmerPaths && !pHKP )
            pHKP = &hkp;
        if ( pHKP )
            buildHKPFromHBV(*pHBV,fwdXlat,revXlat,pHKP);

        if ( pReadPaths )
        {
            pReadPaths->clear();
            pReadPaths->resize(reads.size());
        }
        if ( pKmerPaths )
        {
            pKmerPaths->clear();
            pKmerPaths->resize(reads.size());
        }
        HyperBasevector& spHBV = *pHBV;
        Pather_sleek<BIGK> pather_slk(reads,spHBV,fwdXlat,revXlat,
                                        pReadPaths,pHKP,pKmerPaths,starters);

        // Alternate way to path, potentially slower but lower mem
        /* parallelFor(0ul,reads.size(), */
        /*         [pather_slk]( size_t readId ) mutable */
        /*         { pather_slk.seed(readId); }); */

        /* Destroy(starters); */
        /* Destroy(fwdXlat); */
        /* Destroy(revXlat); */

        /* parallelFor(0ul,reads.size(), */
        /*         [pather_slk]( size_t readId ) mutable */
        /*         { pather_slk.extend(readId); }); */

        parallelForBatch(0ul,reads.size(),10000,pather_slk);
        MEM(after_pathing);
    }
}

template <unsigned BIGK>
void readsToHBV( vecbvec const& reads, unsigned coverage,
                    HyperBasevector* pHBV, ReadPathVec* pReadPaths,
                    HyperKmerPath* pHKP, vecKmerPath* pKmerPaths )
{
    if ( pReadPaths )
    {
        pReadPaths->clear();
        pReadPaths->resize(reads.size());
    }
    if ( pKmerPaths )
    {
        pKmerPaths->clear();
        pKmerPaths->resize(reads.size());
    }

    size_t nKmers = reads.getKmerCount(BIGK);
    if ( !nKmers )
    {
        pHBV->Clear();
        if ( pHKP ) pHKP->Clear();
        return;
    }

    MEM(before_big_dict);

    BigDict<BIGK> bigDict(nKmers/coverage);

    MEM(before_big_kmerizer);

    BigKMerizer<BIGK> kmerizer(&bigDict);
    parallelFor(0ul,reads.size(),
            [kmerizer,&reads]( size_t readId ) mutable
            { kmerizer.kmerize(reads[readId]); });

    MEM(after_kmerizer);

    vecbvec edges;
    edges.reserve(bigDict.size()/100);
    BigKEdgeBuilder<BIGK>::buildEdges(bigDict,&edges);
    AssertEq(edges.getKmerCount(BIGK),bigDict.size());
    vec<int> fwdXlat;
    vec<int> revXlat;

    MEM(after_edges);

    if ( !pReadPaths && !pKmerPaths )
    {
        buildHBVFromEdges(edges,BIGK,pHBV,&fwdXlat,&revXlat);
        if ( pHKP )
            buildHKPFromHBV(*pHBV,fwdXlat,revXlat,pHKP);
        MEM(after_build);
    }
    else
    {
        parallelForBatch(0ul,edges.size(),100,
                [&edges,&kmerizer]( size_t edgeId )
                { kmerizer.updateDict(edges[edgeId]); });

        buildHBVFromEdges(edges,BIGK,pHBV,&fwdXlat,&revXlat);

        MEM(after_build2);

        HyperKmerPath hkp;
        if ( pKmerPaths && !pHKP )
            pHKP = &hkp;
        if ( pHKP )
            buildHKPFromHBV(*pHBV,fwdXlat,revXlat,pHKP);

        Pather<BIGK> pather(reads,bigDict,edges,fwdXlat,revXlat,
                                pReadPaths,pHKP,pKmerPaths);
        parallelForBatch(0ul,reads.size(),10000,pather);

        MEM(after_pathing);
    }
}

template <int BIGK>
struct SillyFunctor
{
    void operator()( vecbvec const& reads, unsigned coverage,
            HyperBasevector* pHBV, ReadPathVec* pReadPaths,
            HyperKmerPath* pHKV, vecKmerPath* pKmerPaths)
    {   readsToHBV<BIGK>(reads,coverage,pHBV,pReadPaths,pHKV,pKmerPaths);  }
};

template <int BIGK>
struct SillyFunctor_sleek
{
    void operator()( vecbvec const& reads, unsigned coverage,
            HyperBasevector* pHBV, ReadPathVec* pReadPaths,
            HyperKmerPath* pHKV, vecKmerPath* pKmerPaths)
    {   readsToHBV_Sleek<BIGK>(reads,coverage,pHBV,pReadPaths,pHKV,pKmerPaths); }
};

template <int BIGK>
struct SillyFunctor_perturb
{
    void operator()( vecbvec const& reads, unsigned coverage,
            HyperBasevector* pHBV, HyperBasevectorX& gHBX, 
            digraphE<vec<int>>& D, ReadPathVecX& pathsx, vec<HyperBasevector>& hbms,
            vec<digraphE<vec<int>>>& Dms, String pi_file, vecbvec& bases)
    {   readsToHBV_Perturb<BIGK>(reads,coverage,pHBV,gHBX,D,pathsx,hbms,Dms,pi_file,bases); }
};

}

void buildBigKHBVFromReads_perturb( unsigned K, vecbvec const& reads, unsigned coverage,
                            HyperBasevector* pHBV, HyperBasevectorX& gHBX, 
                            digraphE<vec<int>>& D, ReadPathVecX& pathsx, vec<HyperBasevector>& hbms,
                            vec<digraphE<vec<int>>>& Dms, String pi_file, vecbvec& bases)
{
    BigK::dispatch<SillyFunctor_perturb>(K,reads,coverage,pHBV,gHBX,D,pathsx,hbms,Dms,pi_file,bases);
}

void buildBigKHBVFromReads_sleek( unsigned K, vecbvec const& reads, unsigned coverage,
                            HyperBasevector* pHBV, ReadPathVec* pReadPaths,
                            HyperKmerPath* pHKP, vecKmerPath* pKmerPaths)
{
    BigK::dispatch<SillyFunctor_sleek>(K,reads,coverage,
                                    pHBV,pReadPaths,pHKP,pKmerPaths);
}

void buildBigKHBVFromReads( unsigned K, vecbvec const& reads, unsigned coverage,
                            HyperBasevector* pHBV, ReadPathVec* pReadPaths,
                            HyperKmerPath* pHKP, vecKmerPath* pKmerPaths)
{
    BigK::dispatch<SillyFunctor>(K,reads,coverage,
                                    pHBV,pReadPaths,pHKP,pKmerPaths);
}
