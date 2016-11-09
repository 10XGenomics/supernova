///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef ASSESS_BEST_ALIGN_CORE_H
#define ASSESS_BEST_ALIGN_CORE_H

#include "Basevector.h"
#include "CoreTools.h"
#include "PackAlign.h"

// PushBoundaries: for tandem-repeat indels, push left/right boundaries so that
// they extend up to the ends of the repeat.

void PushBoundaries( basevector b1, basevector b2,
     int& left1, int& right1, int& left2, int& right2 );

void DecomposeAlign( const align& a, const basevector& b1, const basevector& b2,
     vec< pair<int,int> > & P1, vec< pair<int,int> >& P2 );

#endif
