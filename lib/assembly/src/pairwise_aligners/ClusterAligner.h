///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef CLUSTER_ALIGNER
#define CLUSTER_ALIGNER

// This is an aligner modeled after QueryLookupTable, but designed for
// in-memory alignment.  It aligns a query (q) to a target 'genome' (G).
//
// It takes as input an argument Glocs that can be created as follows:
//
//   Glocs.resize( IPow( 4, K ) ); 
//   for ( size_t i = 0; i < G.size( ); i++ )
//   {    for ( int j = 0; j <= G[i].isize( ) - K; j++ )
//        {    int n = KmerId( G[i], K, j );
//             Glocs[n].push( i, j );    }    }

#include "Basevector.h"
#include "CoreTools.h"
#include "IntPairVec.h"
#include "lookup/LookAlign.h"

void ClusterAligner( basevector q, const vecbasevector& G, const int K, 
     const VecIntPairVec& Glocs, vec<look_align>& aligns,
     const Bool FW_ONLY = False,
     const int BW_ADD = 10,
     const int MIN_CLUSTER = 0,
     const int MAX_OFFSET_DIFF = 500,
     const int MISMATCH_PENALTY = 1,
     const int GAP_PENALTY = 1 );

#endif
