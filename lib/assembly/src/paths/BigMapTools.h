///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef BIG_MAP_TOOLS_H
#define BIG_MAP_TOOLS_H

#include "Basevector.h"
#include "CoreTools.h"
#include "IntPairVec.h"
#include "graph/Digraph.h"
#include "paths/LongReadTools.h"

// A placementy represents an alignment of sequence u to sequence g.
// A placementy is allowed to go around the 'end' of a circular reference sequence, 
// in which case it will have Pos < pos.  Also note that in the rc case, 

class placementy {
     public:
     placementy( ) { }
     placementy( const int u, const int g, const int pos, const int Pos, 
          const Bool fw, const int mismatches, const int indels, const int unmapped )
          : u(u), g(g), pos(pos), Pos(Pos), fw(fw), mismatches(mismatches),
          indels(indels), unmapped(unmapped) { }
     int u;
     int g;
     int pos;
     int Pos;
     Bool fw;
     int mismatches;
     int indels;
     int unmapped;

     int Errs( ) const { return mismatches + indels + unmapped; }

     friend Bool operator==( const placementy& p1, const placementy& p2 )
     {    return p1.u == p2.u && p1.g == p2.g && p1.pos == p2.pos
               && p1.Pos == p2.Pos && p1.fw == p2.fw 
               && p1.mismatches == p2.mismatches && p1.indels == p2.indels
               && p1.unmapped == p2.unmapped;    }

     friend Bool operator<( const placementy& p1, const placementy& p2 )
     {    if ( p1.u < p2.u ) return True;
          if ( p1.u > p2.u ) return False;
          if ( p1.g < p2.g ) return True;
          if ( p1.g > p2.g ) return False;
          if ( p1.pos < p2.pos ) return True;
          if ( p1.pos > p2.pos ) return False;
          if ( p1.Pos < p2.Pos ) return True;
          if ( p1.Pos > p2.Pos ) return False;
          return p1.fw < p2.fw;    }

};

vec<placementy> FindGenomicPlacementsY( const int u, basevector b, const int L,
     const vecbasevector& genome, const VecIntPairVec& Glocs );

#endif
