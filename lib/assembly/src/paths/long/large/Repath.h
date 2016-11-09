///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef LARGE_REPATH_H
#define LARGE_REPATH_H

#include "Basevector.h"
#include "CoreTools.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"

void Repath( const HyperBasevector& hb, const vecbasevector& edges, 
     const vec<int>& inv, ReadPathVec& paths, const int K, const int K2, 
     const String& BIGKHBV0, const Bool REPATH_TRANSLATE, bool INVERT_PATHS,
     const Bool EXTEND_PATHS );

#endif
