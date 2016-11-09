///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file LongReadsToPaths.cc
 * \author tsharpe
 * \date Aug 21, 2012
 *
 * \brief
 */

#include "paths/long/LongReadsToPaths.h"
#include "paths/long/LargeKDispatcher.h"
#include "kmers/BigKPather.h"
#include "kmers/LongReadPather.h"
#include "system/Assert.h"
#include <algorithm>
#include <iterator>
#include <vector>

namespace
{
void extendPath( KmerPath& path, kmer_id_t kmerStart, kmer_id_t kmerStop )
{
    long const MAXDIFF = KmerPathInterval::maxLength - 1;
    path.reserve(path.size()+(kmerStop-kmerStart+MAXDIFF-1)/MAXDIFF);
    while ( kmerStop-kmerStart > MAXDIFF )
    {
        kmer_id_t stop = kmerStart + MAXDIFF;
        path.AddSegmentNoConcatenate(KmerPathInterval(kmerStart,stop));
        kmerStart = stop+1;
    }
    path.AddSegmentNoConcatenate(KmerPathInterval(kmerStart,kmerStop));
}

template <unsigned K>
class KmerNumberingScheme
{
    // An EdgeInterval holds the starting and ending kmer_id_t for some edge.
    // We can't just use a KmerPathInterval, because those have a size limit too
    //   small to represent long edges.
    // This is the same number as the Long::Graph KmerID, except for palindromes.
    // For palindromic edges, which are necessarily a single kmer, we choose a
    //   kmer_id_t in the palindrome region.
    typedef std::pair<kmer_id_t,kmer_id_t> EdgeInterval;

public:
    KmerNumberingScheme( Long::Graph<K> const& graph )
    : mGraph(graph), mNextPalindrome(first_palindrome)
    { mEdgeIntervals.reserve(graph.getNEdges());
      for ( auto itr(graph.begin()), end(graph.end()); itr != end; ++itr )
      { Long::GraphEdge const& edge(*itr);
        if ( !edge.isPalindrome() )
          mEdgeIntervals.push_back(EdgeInterval(edge.getInitialKmerId(),
                                                edge.getFinalKmerId()));
        else
        { mEdgeIntervals.push_back(EdgeInterval(mNextPalindrome,
                                                mNextPalindrome));
          mNextPalindrome += 1; } } }

    void buildHyperKmerPath( HyperKmerPath* pHKP ) const;

    kmer_id_t getStart( size_t edgeId, bool rc ) const
    { EdgeInterval const& edgeInterval = mEdgeIntervals[edgeId];
      return rc ? reverse_kmer(edgeInterval.second) : edgeInterval.first; }

private:
    void addEdge( size_t edgeId, vecbvec const& ends, std::vector<int>& vIds,
                    int& nextVId, HyperKmerPath* pHKP ) const;

    Long::Graph<K> const& mGraph;
    std::vector<EdgeInterval> mEdgeIntervals;
    kmer_id_t mNextPalindrome;
};

template <unsigned K>
void KmerNumberingScheme<K>::buildHyperKmerPath( HyperKmerPath* pHKP ) const
{
    // Try to make the HKP identical even if the edges don't always occur
    // in the same order in the Long::Graph.
    // Tag-sort to get a canonical order.
    Long::Graph<K> const& graph = mGraph;
    size_t nEdges = graph.getNEdges();
    vec<size_t> edgeOrder(nEdges,vec<size_t>::IDENTITY);
    std::sort(edgeOrder.begin(),edgeOrder.end(),
            [&graph](size_t idx1,size_t idx2)
            { size_t edge1Len = graph.getEdge(idx1).getLength();
              size_t edge2Len = graph.getEdge(idx2).getLength();
              if ( edge1Len > edge2Len ) return true;
              if ( edge1Len < edge2Len ) return false;
              return std::lexicographical_compare(graph.basesBegin(idx1),
                                                  graph.basesEnd(idx1),
                                                  graph.basesBegin(idx2),
                                                  graph.basesEnd(idx2)); });

    // Each vertex is associated with a kmer of length K-1.  Figure out the
    // unique set of those for each end of each edge, and its reverse comp.
    vecbvec ends;
    ends.reserve(4*nEdges);
    bvec tmp(K-1);
    for ( size_t idx = 0; idx != nEdges; ++idx )
    {
        auto itr = graph.basesBegin(idx);
        tmp.assign(itr,itr+K-1);
        ends.push_back(tmp);
        ends.push_back(tmp.ReverseComplement());
        itr = graph.basesEnd(idx);
        tmp.assign(itr-K+1,itr);
        ends.push_back(tmp);
        ends.push_back(tmp.ReverseComplement());
    }
    ends.UniqueSort();
    std::vector<int> vIds(ends.size(),-1);

    // build the HKP
    pHKP->Clear();
    pHKP->SetK(K);
    pHKP->AddVertices(ends.size());
    size_t hkpEdgeCount = 2*nEdges-(mNextPalindrome-first_palindrome);
    pHKP->EdgesMutable().reserve(hkpEdgeCount);

    // Process the canonically-ordered list of edges in order.
    // When we find an edge that hasn't yet been added to the HKP, we add it,
    // and then do a breadth-first traversal of the graph to add the other edges
    // in the same component.
    std::vector<bool> edgeDone(nEdges,false);
    std::list<size_t> edgeIdQueue;
    auto oItr = std::back_inserter(edgeIdQueue);
    int nextVId = 0;
    for ( size_t edgeId : edgeOrder )
    {
        if ( edgeDone[edgeId] )
            continue;
        *oItr++ = edgeId;
        while ( !edgeIdQueue.empty() )
        {
            edgeId = edgeIdQueue.front();
            edgeIdQueue.pop_front();
            if ( !edgeDone[edgeId] )
            {
                edgeDone[edgeId] = true;
                addEdge(edgeId,ends,vIds,nextVId,pHKP);
                graph.getEdge(edgeId).getAllConnections(oItr);
            }
        }
    }
    ForceAssertEq(hkpEdgeCount,pHKP->Edges().size());
    ForceAssertEq(nextVId,pHKP->N());
}

template <unsigned K>
void KmerNumberingScheme<K>::addEdge( size_t edgeId, vecbvec const& ends,
                                        std::vector<int>& vIds, int& nextVId,
                                        HyperKmerPath* pHKP ) const
{
    auto itr = mGraph.basesBegin(edgeId);
    bvec tmp(itr,itr+K-1);
    auto beg = std::lower_bound(ends.begin(),ends.end(),tmp);
    AssertEq(*beg,tmp);
    itr = mGraph.basesEnd(edgeId);
    bvec tmp2(itr-K+1,itr);
    auto end = std::lower_bound(ends.begin(),ends.end(),tmp2);
    AssertEq(*end,tmp2);
    int& v1 = vIds[beg-ends.begin()];
    if ( v1 == -1 ) v1 = nextVId++;
    int& v2 = vIds[end-ends.begin()];
    if ( v2 == -1 ) v2 = nextVId++;
    EdgeInterval const& ei = mEdgeIntervals[edgeId];
    KmerPath path;
    extendPath(path,ei.first,ei.second);
    pHKP->AddEdge(v1,v2,path);
    if ( !is_palindrome(ei.first) )
    {
        tmp2.ReverseComplement();
        beg = std::lower_bound(ends.begin(),ends.end(),tmp2);
        AssertEq(*beg,tmp2);
        tmp.ReverseComplement();
        end = std::lower_bound(ends.begin(),ends.end(),tmp);
        AssertEq(*end,tmp);
        int& v3 = vIds[beg-ends.begin()];
        if ( v3 == -1 ) v3 = nextVId++;
        int& v4 = vIds[end-ends.begin()];
        if ( v4 == -1 ) v4 = nextVId++;
        path.Clear();
        extendPath(path,reverse_kmer(ei.second),reverse_kmer(ei.first));
        pHKP->AddEdge(v3,v4,path);
    }
}


template <unsigned K>
void KImpl( vecbvec const& reads, unsigned coverage, unsigned nThreads,
                vecKmerPath& paths, vecKmerPath& paths_rc,
                HyperKmerPath& hkp )
{
    Long::KmerDict<K> dict( reads.SizeSum(), coverage );
    dict.addReads(reads.begin(),reads.end(),nThreads);

    Long::Graph<K> graph(dict);
    //graph.validate(dict);
    KmerNumberingScheme<K> scheme(graph);
    scheme.buildHyperKmerPath(&hkp);
    paths.clear().reserve(reads.size());
    paths_rc.clear().reserve(reads.size());
    KmerPath path;
    typedef vecbvec::const_iterator Itr;
    for ( Itr itr(reads.begin()), end(reads.end()); itr != end; ++itr )
    {
        size_t rdLen = itr->size();
        if ( rdLen >= K )
        {
            rdLen -= K - 1;
            bool rc = false;
            Long::KmerDictEntry const& entry = dict.lookup(itr->begin(),&rc);
            size_t edgeId = entry.getEdgeID();
            Long::GraphEdge const& edge1 = graph.getEdge(edgeId);
            size_t offset = rc ? edge1.getFinalKmerId()-entry.getKmerID() :
                                    entry.getKmerID()-edge1.getInitialKmerId();
            size_t edgeLen = edge1.getLength() - offset;
            kmer_id_t start = scheme.getStart(edgeId,rc)+offset;
            if ( edgeLen >= rdLen )
                extendPath(path,start,start+rdLen-1);
            else
            {
                extendPath(path,start,start+edgeLen-1);
                bvec::const_iterator pBase = itr->begin()+edgeLen+K-1;
                while ( (rdLen -= edgeLen) )
                {
                    edgeId = graph.nextEdgeId(edgeId,*pBase,&rc);
                    Long::GraphEdge const& edge = graph.getEdge(edgeId);
                    edgeLen = std::min(edge.getLength(),rdLen);
                    if ( edge.isPalindrome() )
                        rc = false;
                    start = scheme.getStart(edgeId,rc);
                    extendPath(path,start,start+edgeLen-1);
                    pBase += edgeLen;
                }
            }
        }
        paths.push_back(path);
        path.ReverseNoConcatenate();
        paths_rc.push_back(path);
        path.clear();
    }

};

template <int K>
struct LRP_Functor
{
    void operator()( vecbvec const& reads, unsigned coverage, unsigned nThreads,
                     vecKmerPath& paths, vecKmerPath& paths_rc,
                     HyperKmerPath& hkp )
    { KImpl<K>(reads,coverage,nThreads,paths,paths_rc,hkp); }
};

}

void LongReadsToPaths( vecbvec const& reads, unsigned k, unsigned coverage,
                        unsigned logLevel, bool useOldLRPMethod,
                        HyperBasevector* pHBV, HyperKmerPath* pHKP,
                        vecKmerPath* pPaths, vecKmerPath* pPathsRC,
                        vec<big_tagged_rpint>* pPathsDB )
{
    if ( pPathsRC && !pPaths )
        FatalErr("Must specify destination for paths if you want pathsRC.");
    if ( pPathsDB && (!pPaths || !pPathsRC) )
        FatalErr("Must specify destination for paths and pathRC if you want a pathsDB");

    if ( useOldLRPMethod )
    {
        HyperKmerPath tmpHKP;
        if ( !pHKP ) pHKP = &tmpHKP;
        vecKmerPath tmpPaths;
        if ( !pPaths ) pPaths = &tmpPaths;
        vecKmerPath tmpPathsRC;
        if ( !pPathsRC ) pPathsRC = &tmpPathsRC;
        Long::gLogLevel = logLevel;
        BigK::dispatch<LRP_Functor>(k,reads,coverage,getConfiguredNumThreads(),
                                    *pPaths,*pPathsRC,*pHKP);
        vec<big_tagged_rpint> tmpPathsdb;
        if ( !pPathsDB ) pPathsDB = &tmpPathsdb;
        CreateDatabase(*pPaths,*pPathsRC,*pPathsDB);
        KmerBaseBrokerBig kbb( k, *pPaths, *pPathsRC, *pPathsDB, reads );
        *pHBV = HyperBasevector( *pHKP, kbb );
    }
    else
    {
        buildBigKHBVFromReads(k,reads,coverage,pHBV,nullptr,pHKP,pPaths);
        if ( pPathsRC )
        {
            pPathsRC->assign(pPaths->begin(),pPaths->end());
            for ( KmerPath& path : *pPathsRC )
                path.ReverseNoConcatenate();
        }
        if ( pPathsDB )
            CreateDatabase(*pPaths,*pPathsRC,*pPathsDB);
    }
}
