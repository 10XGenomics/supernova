///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "paths/long/FriendAligns.h"
#include "paths/long/FriendAlignFinder.h"
#include "paths/long/FriendAlignFinderQ.h"
#include "paths/long/FriendAlignFinderNaif.h"
#include "paths/long/MakeAlignments.h"

namespace
{
class FriendsCache : public FriendAlignerImpl
{
public:
    FriendsCache( vec<simple_align_data> const& aligns, size_t nReads )
    : mFriendsVec(nReads)
    { auto itr = aligns.begin();
      auto end = aligns.end();
      Friends tmp;
      for ( size_t idx = 0; idx != nReads; ++idx )
      { tmp.clear();
        while ( itr != end && static_cast<size_t>(itr->id1) == idx )
          tmp.push_back(Friend(itr->id2,itr->offset,itr->rc2));
        ForceAssert( itr==end || static_cast<size_t>(itr->id1) > idx );
        mFriendsVec[idx] = tmp; } }

    virtual void getAligns( size_t readId, Friends* pFriends )
    { *pFriends = mFriendsVec[readId]; }

private:
    MasterVec<Friends> mFriendsVec;
};

template <class C>
class CShim
{
public:
    CShim( C const* pC ) : mpC(pC) {}
    typename C::size_type size() const { return mpC->size(); }
    typename C::const_iterator begin( size_t idx=0 ) const { return mpC->begin(idx); }
    typename C::const_iterator end() const { return mpC->end(); }
    size_t getKmerCount( unsigned K ) const { return mpC->getKmerCount(K); }
private:
    C const* mpC;
};

}


FriendAligner::FriendAligner( vecbvec const& bases, vecqvec const& quals,
                               vec<Bool> const& toEdit,
                               String const& friendsCache,
                               int MAKE_ALIGN_IMPL, unsigned const K,
                               unsigned const min_freq, unsigned const max_freq,
                               unsigned const min_qual, unsigned const coverage,
                               bool downSample, int verbosity )
{
    if ( MAKE_ALIGN_IMPL == 3 )
    {
        FriendAlignFinderQ::storeFriends( CShim<vecbvec>(&bases),
                                          CShim<vecqvec>(&quals),
                                          friendsCache, min_qual, coverage, K,
                                          min_freq, max_freq, verbosity );
        mpImpl = new FriendAlignFinderQ(friendsCache);
    }
    else if ( MAKE_ALIGN_IMPL == 2 )
    {
        if ( K <= 29 )
            mpImpl = new FriendAlignFinderNaif<Kmer29>(friendsCache,K,bases,max_freq);
        else if ( K <= 60 )
            mpImpl = new FriendAlignFinderNaif<Kmer60>(friendsCache,K,bases,max_freq);
        else if ( K <= 124 )
            mpImpl = new FriendAlignFinderNaif<Kmer124>(friendsCache,K,bases,max_freq);
        else {
            cout << "\nIllegal K value for MakeAlignments." << endl;
            CRD::exit(1);
        }
    }
    else if ( MAKE_ALIGN_IMPL==1 )
    {
        if      ( K == 12 ) mpImpl = new FriendAlignFinder<12>( bases, max_freq, downSample, verbosity );
        else if ( K == 16 ) mpImpl = new FriendAlignFinder<16>( bases, max_freq, downSample, verbosity );
        else if ( K == 24 ) mpImpl = new FriendAlignFinder<24>( bases, max_freq, downSample, verbosity );
        else if ( K == 28 ) mpImpl = new FriendAlignFinder<28>( bases, max_freq, downSample, verbosity );
        else if ( K == 40 ) mpImpl = new FriendAlignFinder<40>( bases, max_freq, downSample, verbosity );
        else if ( K == 60 ) mpImpl = new FriendAlignFinder<60>( bases, max_freq, downSample, verbosity );
        else if ( K == 80 ) mpImpl = new FriendAlignFinder<80>( bases, max_freq, downSample, verbosity );
        else {    
            cout << "\nIllegal K value for MakeAlignments." << endl;
            CRD::exit(1);
        }
    }
    else if ( MAKE_ALIGN_IMPL==0 )
    {
        vec<simple_align_data> aligns_all;
        MakeAlignments( K, max_freq, bases, toEdit, aligns_all, verbosity );
        mpImpl = new FriendsCache(aligns_all,bases.size());
        if ( verbosity ) cout << Date( ) << ": read starts defined" << endl;
    }
    else
        FatalErr("Illegal value for MAKE_ALIGN_EXP (must be 0, 1, 2, or 3).");
}

FriendAlignerImpl::~FriendAlignerImpl()
{}

#include "feudal/SmallVecDefs.h"
template class SmallVec<Friend,MempoolAllocator<Friend>>;

#include "feudal/OuterVecDefs.h"
template class OuterVec<Friends>;
