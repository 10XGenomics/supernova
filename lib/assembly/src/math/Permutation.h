///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// This file defines class "Permutation".  It is ONE-BASED, for backward
// compatibility.

#ifndef PERMUTATION_H
#define PERMUTATION_H

#include "CoreTools.h"

class Permutation : public vec<int> {

     public:

     Permutation(int);   // construct the identity Permutation
     Permutation( ) { }

};

#endif
