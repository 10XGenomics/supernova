///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
//

#include "feudal/PQVec.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"

#ifndef ShortKmerReadPather_H_
#define ShortKmerReadPather_H_

class ShortKmerReadPather {

public:

    static std::tuple<unsigned, unsigned, unsigned>  // score, mismatch, length
    ScoreEdge( bvec const& bases, qvec const& quals,
	       int offset, bvec const& edge) ;

    static void FindPaths(const vecbasevector& bases, const VecPQVec& quals,
			  const HyperBasevector& hbv, ReadPathVec& paths,
			  const int k, const int num_threads, 
			  const vec<size_t>& ids = vec<size_t>(), bool debug = false) ;

};

#endif /* ShortKmerReadPather */
