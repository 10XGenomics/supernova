///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * HBVFromEdges.cc
 *
 *  Created on: Dec 13, 2013
 *      Author: tsharpe
 */

#include "paths/long/HBVFromEdges.h"
#include "MapReduceEngine.h"
#include "dna/CanonicalForm.h"
#include "feudal/HashSet.h"
#include "math/Hash.h"
#include "paths/KmerPathInterval.h"
#include "system/ID.h"
#include <list>
#include "10X/DfTools.h"      // TODO: REMOVE

namespace
{

class EdgeEnd
{
public:
    EdgeEnd()=default;
    EdgeEnd( bvec const* pBV, bool rc, bool distal, unsigned K )
    : mItr(pBV,distal?pBV->size()-K:0,rc), mK(K), mHash(FNV1a(mItr,mItr+K)) {}

    size_t getHash() const { return mHash; }
    bvec::SwitchHitterIter begin() const { return mItr; }
    bvec::SwitchHitterIter end() const { return mItr+mK; }

    struct Hasher
    { typedef EdgeEnd argument_type;
      size_t operator()( EdgeEnd const& ee ) const { return ee.getHash(); } };

    friend bool operator<( EdgeEnd const& ee1, EdgeEnd const& ee2 )
    { if ( ee1.mHash != ee2.mHash ) return ee1.mHash < ee2.mHash;
      for ( auto i1=ee1.begin(),e1=ee1.end(),i2=ee2.begin(); i1!=e1; ++i1,++i2 )
        if ( *i1 != *i2 ) return *i1 < *i2;
      return false; }

    friend bool operator==( EdgeEnd const& ee1, EdgeEnd const& ee2 )
    { if ( ee1.mHash != ee2.mHash ) return false;
      for ( auto i1=ee1.begin(),e1=ee1.end(),i2=ee2.begin(); i1!=e1; ++i1,++i2 )
        if ( *i1 != *i2 ) return false;
      return true; }

private:
    bvec::SwitchHitterIter mItr;
    unsigned mK;
    size_t mHash;
};

class IandO
{
public:
    IandO()=default;
    IandO( size_t edgeId, bool rc ) : mEdgeID(edgeId), mRC(rc) {}

    int getId() const { return mEdgeID.val(); }
    bool isRC() const { return mRC; }

private:
    ID<5> mEdgeID;
    bool mRC;
};

// edges and orientations connected to a vertex
class Vertex : public EdgeEnd
{
    static int const UNASSIGNED = -1;
    static size_t const MAX_EDGES = 8;

public:
    Vertex() : mId(UNASSIGNED), mNEdges(0) {}
    Vertex( EdgeEnd const& ee ) : EdgeEnd(ee), mId(UNASSIGNED), mNEdges(0) {}

    int getId() const { return mId; }
    int getId( int& nextId )
    { if ( mId == UNASSIGNED ) mId = nextId++;
      return mId; }

    void addEdge( size_t edgeId, bool rc )
    { ForceAssertLt(mNEdges,MAX_EDGES);
      mEdges[mNEdges++] = IandO(edgeId,rc); }

    IandO const* begin() const { return mEdges; }
    IandO const* end() const { return mEdges+mNEdges; }

private:
    int mId;
    size_t mNEdges;
    IandO mEdges[MAX_EDGES];
};
size_t const Vertex::MAX_EDGES;

typedef HashSet<Vertex,EdgeEnd::Hasher,std::equal_to<EdgeEnd>> VertexDict;

struct BVComp
{
    bool operator()( bvec const& bv1, bvec const& bv2 )
    { if ( bv1.size() != bv2.size() ) return bv1.size() > bv2.size();
      return bv1 < bv2; }
};

struct EEComp
{
    bool operator()( EdgeEnd const& ee1, EdgeEnd const& ee2 )
    { auto itr1 = ee1.begin(), itr2 = ee2.begin();
      if ( &itr1.container() != &itr2.container() )
          return BVComp()(itr1.container(),itr2.container());
      if ( itr1.isRC() != itr2.isRC() ) return itr1.isRC() < itr2.isRC();
      return itr1.pos() < itr2.pos(); }
};

class EOComp
{
public:
    EOComp( vecbvec const& edges )
    : mEdges(edges) {}

    bool operator()( size_t idx1, size_t idx2 )
    { return BVComp()(mEdges[idx1],mEdges[idx2]); }

private:
    vecbvec const& mEdges;
};

class VertexDictBuilder
{
public:
    VertexDictBuilder( vecbvec const& edges, unsigned K, VertexDict* pDict )
    : mEdges(edges), mKLO(K-1), mpDict(pDict) {}

    template <class OItr>
    void map( size_t edgeId, OItr oItr )
    { bvec const& edge = mEdges[edgeId];
      *oItr = EdgeEnd(&edge,false,false,mKLO); ++oItr;
      *oItr = EdgeEnd(&edge,false,true,mKLO); ++oItr;
      if ( edge.getCanonicalForm() != CanonicalForm::PALINDROME )
      { *oItr = EdgeEnd(&edge,true,false,mKLO); ++oItr;
        *oItr = EdgeEnd(&edge,true,true,mKLO); ++oItr; } }

    void reduce( EdgeEnd* key1, EdgeEnd* key2 )
    { bvec const* pBV0 = &*mEdges.begin();
      std::sort(key1,key2,EEComp());
      mpDict->apply(*key1,
            [key1,key2,pBV0]( Vertex const& vtx )
            { Vertex& vertex = const_cast<Vertex&>(vtx);
              for ( auto key=key1; key != key2; ++key )
              { auto itr = key->begin();
                vertex.addEdge(&itr.container()-pBV0,itr.isRC()); } }); }

    EdgeEnd* overflow( EdgeEnd* key1, EdgeEnd* key2 )
    { return key2; }

private:
    vecbvec const& mEdges;
    unsigned mKLO;
    VertexDict* mpDict;
};

class HBVBuilder
{
    static int const NULL_EDGEID = -1;
public:
    HBVBuilder( vecbvec const& edges, unsigned K, VertexDict const& dict,
                    HyperBasevector* pHBV, vec<int>* pFwd, vec<int>* pRev )
    : mEdges(edges), mKLO(K-1), mDict(dict),
      mpHBV(pHBV), mFwd(*pFwd), mRev(*pRev), mNextVId(0)
    { size_t nEdges = mEdges.size();
      mpHBV->Clear();
      mpHBV->SetK(K);
      mpHBV->AddVertices(dict.size());
      mpHBV->EdgesMutable().reserve(2*nEdges);
      mFwd.clear(); mFwd.resize(nEdges,NULL_EDGEID);
      mRev.clear(); mRev.resize(nEdges,NULL_EDGEID); }

    ~HBVBuilder()
    { AssertEq(mNextVId,mpHBV->N());
      Assert(std::find(mFwd.begin(),mFwd.end(),NULL_EDGEID)==mFwd.end());
      Assert(std::find(mRev.begin(),mRev.end(),NULL_EDGEID)==mRev.end()); }

    void add( IandO const& iAndO )
    { if ( !isDone(iAndO) )
      { mQueue.push_back(iAndO); processQueue(); } }

private:
    bool isDone( IandO const& iAndO ) const
    { return (iAndO.isRC() ? mRev : mFwd)[iAndO.getId()] != NULL_EDGEID; }

    void processQueue()
    { while ( !mQueue.empty() )
      { IandO iAndO = mQueue.front();
        mQueue.pop_front();
        if ( isDone(iAndO) ) continue;

        bvec const& edge = mEdges[iAndO.getId()];
        bool rc = iAndO.isRC();
        Vertex* pV1 =
                const_cast<Vertex*>(mDict.lookup(EdgeEnd(&edge,rc,false,mKLO)));
        ForceAssert(pV1);
        int v1 = pV1->getId(mNextVId);
        Vertex* pV2 =
                const_cast<Vertex*>(mDict.lookup(EdgeEnd(&edge,rc,true,mKLO)));
        ForceAssert(pV2);
        int v2 = pV2->getId(mNextVId);
        int newEdgeId = mpHBV->EdgeObjectCount();
        mpHBV->AddEdge(v1,v2,edge);
        if ( rc ) mpHBV->EdgeObjectMutable(newEdgeId).ReverseComplement();
        bool isPalindrome = edge.getCanonicalForm()==CanonicalForm::PALINDROME;
        if ( !rc || isPalindrome ) mFwd[iAndO.getId()] = newEdgeId;
        if ( rc || isPalindrome ) mRev[iAndO.getId()] = newEdgeId;
        for ( auto itr=pV1->begin(),end=pV1->end(); itr != end; ++itr )
          if ( !isDone(*itr) )
            mQueue.push_back(*itr);
        for ( auto itr=pV2->begin(),end=pV2->end(); itr != end; ++itr )
          if ( !isDone(*itr) )
            mQueue.push_back(*itr);
      }
    }

    vecbvec const& mEdges;
    unsigned mKLO;
    VertexDict const& mDict;
    HyperBasevector* mpHBV;
    vec<int>& mFwd;
    vec<int>& mRev;
    std::list<IandO> mQueue;
    int mNextVId;
};
int const HBVBuilder::NULL_EDGEID;


}

void buildHBVFromEdges( vecbvec const& edges, unsigned K, HyperBasevector* pHBV,
                             vec<int>* pFwdEdgeXlat, vec<int>* pRevEdgeXlat )
{
    if ( edges.empty() )
    {
        pHBV->Clear();
        pFwdEdgeXlat->clear();
        pRevEdgeXlat->clear();
        return;
    }

    // find all the K-1 overlaps on the edge ends: those sequences mark vertices
    size_t nEdges = edges.size();
    MEM(v0);
    VertexDict dict(2*nEdges);
    if ( true )
    {
        cout << Date( ) << ": vertex dict building" << endl;
        MEM(v1);
        VertexDictBuilder vdb(edges,K,&dict);
        MEM(v2);
        MapReduceEngine<VertexDictBuilder,EdgeEnd,EdgeEnd::Hasher> mre(vdb);
        MEM(v3);
        if ( !mre.run(4*nEdges,0ul,nEdges) )
            FatalErr("Map/Reduce operation failed in finding vertices.");
        MEM(v4);
    }

    // try to make the graph identical even if the edges don't always occur
    // in the same order -- tag sort on edge length descending, then lexically.

    cout << Date( ) << ": sort lexicomatic" << endl;
    vec<size_t> edgeOrder(nEdges,vec<size_t>::IDENTITY);
    std::sort(edgeOrder.begin(),edgeOrder.end(),EOComp(edges));
    MEM(v5)

    cout << Date( ) << ": build it dude" << endl;
    HBVBuilder hbvb(edges,K,dict,pHBV,pFwdEdgeXlat,pRevEdgeXlat);
    MEM(v6);
    for ( size_t edgeId : edgeOrder )
    {
        IandO iAndO(edgeId,false);
        hbvb.add(iAndO);
    }
    MEM(v7);
    for ( size_t edgeId : edgeOrder )
    {
        IandO iAndO(edgeId,true);
        hbvb.add(iAndO);
    }
    MEM(v8);

}


void buildHKPFromHBV( HyperBasevector const& hbv,
                        vec<int> const& fwdEdgeXlat,
                        vec<int> const& revEdgeXlat,
                        HyperKmerPath* pHKP )
{
    unsigned K = hbv.K();
    pHKP->Clear();
    pHKP->SetK(K);
    pHKP->ToMutable() = hbv.To();
    pHKP->FromMutable() = hbv.From();
    pHKP->ToEdgeObjMutable() = hbv.ToEdgeObj();
    pHKP->FromEdgeObjMutable() = hbv.FromEdgeObj();
    pHKP->EdgesMutable().resize(hbv.EdgeObjectCount());
    long nextPalindrome = first_palindrome;
    long nextKmerId = 1;
    auto oItr = pHKP->EdgesMutable().begin();
    vec<bvec> const edges = hbv.Edges();
    for ( auto itr=edges.begin(),end=edges.end(); itr != end; ++itr,++oItr )
    {
        size_t len = itr->size();
        switch ( itr->getCanonicalForm() )
        {
        case CanonicalForm::FWD:
            oItr->Assign(nextKmerId,nextKmerId+len-K);
            nextKmerId += len;
            break;
        case CanonicalForm::PALINDROME:
            oItr->Assign(nextPalindrome,nextPalindrome+len-K);
            nextPalindrome += len;
            break;
        case CanonicalForm::REV:
            break;
        }
    }
    auto rItr = revEdgeXlat.begin();
    auto end = fwdEdgeXlat.end();
    for ( auto itr=fwdEdgeXlat.begin(); itr != end; ++itr,++rItr )
    {
        if ( pHKP->EdgeObject(*itr).empty() )
        {
            KmerPath& edge = pHKP->EdgeObjectMutable(*itr);
            edge = pHKP->EdgeObject(*rItr);
            edge.ReverseNoConcatenate();
        }
        else if ( pHKP->EdgeObject(*rItr).empty() )
        {
            KmerPath& edge = pHKP->EdgeObjectMutable(*rItr);
            edge = pHKP->EdgeObject(*itr);
            edge.ReverseNoConcatenate();
        }
        else
            ForceAssertEq(*itr,*rItr);
    }
}
