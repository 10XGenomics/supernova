///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file RefLocus.h
 * \author tsharpe
 * \date Oct 18, 2012
 *
 * \brief
 */
#ifndef REFLOCUS_H_
#define REFLOCUS_H_

#include <cstddef>
#include <ostream>
#include <vector>

struct RefLocus
{
    RefLocus( size_t refID, unsigned offset, unsigned length, bool rc = false)
    : mRefID(refID), mOffset(offset), mLength(length), mRC(rc) {}
    RefLocus() : mRefID(0), mOffset(0U), mLength(0U) {};

    size_t mRefID;
    unsigned mOffset;
    unsigned mLength;
    bool mRC;

    friend bool operator< ( const RefLocus& l, const RefLocus& r )
    { if ( l.mRefID != r.mRefID ) return l.mRefID < r.mRefID;
      if ( l.mOffset != r.mOffset ) return l.mOffset < r.mOffset;
      return false; }

    friend ostream& operator<<( ostream& os, RefLocus const& loc )
    { os << "RefLocus( mRefID=" << loc.mRefID << ", mOffset=" << loc.mOffset
         << ", mLength=" << loc.mLength << ", mRC=" << ( loc.mRC?"rc":"fw" )
         << " )";
      return os; }
};
TRIVIALLY_SERIALIZABLE(RefLocus);
typedef std::vector<RefLocus> RefLocusVec;

#endif /* REFLOCUS_H_ */
