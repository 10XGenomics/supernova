///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
#ifndef FRIEND_ALIGNS_H
#define FRIEND_ALIGNS_H

#include "Basevector.h"
#include "Qualvector.h"
#include "String.h"
#include "Vec.h"
#include "feudal/OuterVec.h"
#include "feudal/SerfVec.h"
#include <ostream>

/// Describes a read that aligns to an implicit target read.
/// A positive offset implies that the friend read should be placed at that
/// offset onto the target read to get a K-length identity:
/// friend1 xxxxxxxxxxxxxxkmerxxxxxx          offset will be negative
/// friend2   xxxxxxxxxxxxkmerxxxxxxxx        offset will be negative
/// target       xxxxxxxxxkmerxxxxxxxxxxx
/// friend3        xxxxxxxkmerxxxxxxxxxxxxx   offset will be positive
/// friend4          xxxxxkmerxxxxxxxxxxxxxxx offset will be positive
/// If isRC() is true, the friend read should be RC'd before aligning.
class Friend
{
public:
    Friend() : mReadId(0), mOffset(0), mRC(false) {}
    Friend( size_t readId, int offset, bool rc )
    : mReadId(readId), mOffset(offset), mRC(rc)
    { ForceAssertLe(readId,std::numeric_limits<unsigned>::max());
      ForceAssertLe(offset,std::numeric_limits<short>::max()); }

    size_t readId() const { return mReadId; }
    int offset() const { return mOffset; }
    bool isRC() const { return mRC; }

    friend bool operator==( Friend const& f1, Friend const& f2 )
    { return f1.mReadId == f2.mReadId &&
            f1.mOffset == f2.mOffset && f1.mRC == f2.mRC; }

    friend bool operator!=( Friend const& f1, Friend const& f2 )
    { return !(f1 == f2); }

    friend bool operator<( Friend const& f1, Friend const& f2 )
    { if ( f1.mReadId < f2.mReadId ) return true;
      if ( f1.mReadId == f2.mReadId )
      { if ( f1.mOffset < f2.mOffset ) return true;
        if ( f1.mOffset == f2.mOffset && f1.mRC < f2.mRC ) return true; }
      return false; }

    friend ostream& operator<<( ostream& os, Friend const& fr )
    { if ( fr.mRC ) os << '~';
      return os << fr.mReadId << '[' << fr.mOffset << ']'; }

private:
    unsigned mReadId;
    short mOffset;
    bool mRC;
};
typedef SerfVec<Friend> Friends;
extern template class SmallVec<Friend,MempoolAllocator<Friend> >;
extern template class OuterVec<Friends>;

class FriendAlignerImpl
{
public:
    virtual ~FriendAlignerImpl();

    // this must be thread-safe.
    virtual void getAligns( size_t readId, Friends* pFriends ) = 0;
};

class FriendAligner
{
  public:
    FriendAligner( vecbvec const& bases, vecqvec const& quals,
                    vec<Bool> const& toEdit,
                    String const& friendsCache,
                    int MAKE_ALIGN_IMPL, unsigned const K,
                    unsigned const min_freq, unsigned const max_freq,
                    unsigned const min_qual, unsigned const coverage,
                    bool downSample, int verbosity );
    ~FriendAligner() { delete mpImpl; }

    void getAligns( size_t readId, Friends* pFriends )
    { mpImpl->getAligns(readId,pFriends); }

  private:
    FriendAlignerImpl* mpImpl;
};

#endif
