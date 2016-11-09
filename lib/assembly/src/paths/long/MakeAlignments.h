///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef MAKE_ALIGNMENTS_H
#define MAKE_ALIGNMENTS_H

#include "Basevector.h"
#include "Vec.h"

class simple_align_data {

     public:

     simple_align_data()=default;
     simple_align_data( int id1, int id2, int offset, bool rc2 )
     : id1(id1), id2(id2), offset(offset), rc2(rc2)
     { ForceAssertLt(offset,std::numeric_limits<short>::max());
       ForceAssertGt(offset,std::numeric_limits<short>::min()); }

     int id1;
     int id2;
     short offset;
     bool rc2;

     friend bool operator<( const simple_align_data& a,
                             const simple_align_data& b )
     { if ( a.id1 < b.id1 ) return true;
       if ( a.id1 > b.id1 ) return false;
       if ( a.id2 < b.id2 ) return true;
       if ( a.id2 > b.id2 ) return false;
       if ( a.offset < b.offset ) return true;
       if ( a.offset > b.offset ) return false;
       return a.rc2 < b.rc2;    }

     friend bool operator==( const simple_align_data& a,
                                 const simple_align_data& b )
     { return a.id1 == b.id1 && a.id2 == b.id2 &&
                 a.offset == b.offset && a.rc2 == b.rc2;    }

};

TRIVIALLY_SERIALIZABLE( simple_align_data );

// MakeAlignments.  For each bases[id1] for which is_target[id1] = True, find
// all gap-free alignments (id1,id2,offset,rc2?) that subsume a perfect match
// of length at least K.  Exclude kmers occurring more than max_freq times.
// Output is sorted.

void MakeAlignments( const int K, const int max_freq, const vecbasevector& bases,
     const vec<Bool>& is_target, vec<simple_align_data>& aligns_all,
     int verbosity );

#endif
