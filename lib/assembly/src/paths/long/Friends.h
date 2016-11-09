///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file Friends.h
 * \author tsharpe
 * \date Apr 24, 2012
 *
 * \brief
 */
#ifndef FRIENDS_H_
#define FRIENDS_H_

#include "Basevector.h"
#include "feudal/BinaryStreamTraits.h"
#include "feudal/SerfVec.h"
#include "feudal/MasterVec.h"
#include "feudal/VirtualMasterVec.h"
#include "String.h"
#include <cstddef>

// ID and orientation
class IdAndOrientation
{
public:
    IdAndOrientation() : mVal(~0ul) {}
    explicit IdAndOrientation( size_t id, bool rc = false )
    : mVal(id<<1) { if ( rc ) mVal |= 1; }

    // compiler-supplied copying and destructor are OK

    size_t getId() const { return mVal>>1; }
    bool isRC() const { return mVal&1; }

    friend bool operator<( IdAndOrientation const& val1,
                                IdAndOrientation const& val2 )
    { return val1.mVal < val2.mVal; }

    friend bool operator==( IdAndOrientation const& val1,
                                IdAndOrientation const& val2 )
    { return val1.mVal == val2.mVal; }

    friend bool operator!=( IdAndOrientation const& val1,
                                IdAndOrientation const& val2 )
    { return val1.mVal != val2.mVal; }

    friend ostream& operator<<( ostream& os, IdAndOrientation const& iAndO )
    { if ( iAndO.isRC() ) os << '~';
      return os << iAndO.getId(); }

private:
    size_t mVal;
};
TRIVIALLY_SERIALIZABLE(IdAndOrientation);

typedef SerfVec<IdAndOrientation> IAndOs; // friendly readIds & orientations
typedef MasterVec<IAndOs> IAndOsVec; // indexed by readId
extern template class SmallVec<IdAndOrientation,MempoolAllocator<IdAndOrientation> >;
extern template class OuterVec<IAndOs>;

/// Write a feudal file of all the friends of the given reads.
/// Friends share a perfect 25-mer.
void FindFriends( vecbvec const& reads, String const& friendsFile,
                    bool useBrokenVersion );

#endif /* FRIENDS_H_ */
