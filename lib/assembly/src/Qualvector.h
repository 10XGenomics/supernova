///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/// \file
/// This file defines the typedef "qualvector", which stores quality scores
/// between 0 and 255 as a vector of unsigned chars, and the typedef
/// "vecqualvector", which stores a vector of qualvectors.
/// \ingroup grp_quals

#ifndef QUALVECTOR
#define QUALVECTOR

#include "feudal/SerfVec.h"
#include "feudal/MasterVec.h"
#include "Charvector.h"
#include "String.h"
#include "Vec.h"
#include <algorithm>
#include <ostream>

/// Logical type for quality scores
typedef unsigned char qual_t;

/// Vector of quality scores, for example representing the quality of each base
/// in one read.
typedef UCharVec QualVec;
typedef QualVec qualvector;
typedef QualVec qvec;


/// Vector of vectors of quality scores, for example representing the quality
/// of each base in each read in a set of reads.
typedef VecUCharVec QualVecVec;
typedef QualVecVec vecqualvector;
typedef QualVecVec vecqvec;

typedef OuterVec< OuterVec<qvec,MempoolAllocator<qual_t> >,
                  MempoolOwner<qual_t> > qvec3;
extern template class OuterVec<qvec>;
extern template class OuterVec<qvec,qvec::alloc_type>;
extern template class OuterVec< OuterVec<qvec,qvec::alloc_type>,
                                                    MempoolOwner<qual_t> >;

///Produces fasta format quals, mirrors basevector::Print()
void Print( std::ostream &out, const qualvector &q, const String &name,
            const int scores_per_line = 25 );

/// Returns two strings representing the quality scores stacked vertically
std::pair <String, String> Stacked(const qualvector& quals) ;

/// Writes two strings representing the quality scores stacked vertically
/// e.g. 43,31,20,2,2,2 becomes: 432   
///                              310222
void PrintStacked(std::ostream &out , const qualvector& quals) ;


/// CopyQuals: copy elements from one qualvector to another. If rev_from=True,
/// then copy from reverse(from), starting at from_start on reverse(from).
/// This mirrors CopyBases in Basevector.h.

inline void CopyQuals( const qualvector& from, int from_start, qualvector& to,
     int to_start, int count, bool rev_from = false )
{    if ( !rev_from )
     {    for ( int i = 0; i < count; i++ )
               to[ to_start + i ] = from[ from_start + i ];    }
     else
     {    for ( int i = 0; i < count; i++ )
               to[ to_start + i ]
                    = from[ from.size( ) - (from_start + i) - 1 ];    }    }

// ReadFastaQuals: read quality scores from a fasta quality score file.  If
// ids_to_read is supplied, it should be a sorted list of indices of the records
// to be read.

void ReadFastaQuals( const String& fn, vecqualvector& qual,
                        const vec<int>* ids_to_read = 0 );

/// Returns a new qvec with each quality score replaced by the minimum quality
/// score in the range idx-radius to idx+radius.
inline qvec Squash( qvec const& qv, unsigned radius )
{ qvec result;
  if ( !qv.empty() )
  { qvec::const_pointer beg = &qv[0];
    qvec::const_pointer end = beg + qv.size();
    result.reserve(end-beg);
    for ( qvec::const_pointer itr(beg); itr != end; ++itr )
      result.push_back(*std::min_element(std::max(beg,itr-radius),
                                         std::min(end,itr+radius))); }
  return result; }

#endif
