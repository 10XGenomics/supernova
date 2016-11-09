///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "CoreTools.h"
#include "Equiv.h"
#include "ParallelVecUtilities.h"
#include "graph/Digraph.h"
#include "kmers/KmerRecord.h"
#include "paths/BigMapTools.h"

vec<placementy> FindGenomicPlacementsY( const int u, basevector b, const int L,
     const vecbasevector& genome, const VecIntPairVec& Glocs )
{    vec<placementy> places;
     for ( int pass = 1; pass <= 2; pass++ )
     {    if ( pass == 2 ) b.ReverseComplement( );
          if ( b.isize( ) >= L )
          {    int n = KmerId( b, L, 0 );
               for ( unsigned z = 0; z < Glocs[n].size( ); z++ )
               {    const basevector& g = genome[ Glocs[n][z].first ];
                    const int gpos = Glocs[n][z].second;
                    if ( !PerfectMatch( b, g, 0, gpos, b.size( ) ) ) continue;
                    places.push( u, Glocs[n][z].first, gpos, 
                         gpos + b.isize( ), pass == 1, 0, 0, 0 );    }    }    }
     return places;    }
