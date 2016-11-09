///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef PRECORRECT_OLD_NEW_H
#define PRECORRECT_OLD_NEW_H

#include "Basevector.h"
#include "CoreTools.h"
#include "Qualvector.h"

void PreCorrectOldNew( vecbasevector* bases, vecqualvector const& quals0,
     const vec<int>& trace_ids );

#endif
