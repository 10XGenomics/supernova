///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// This is an annoying little class that does nothing more than hold an integer.

#ifndef KMER_COUNT_H
#define KMER_COUNT_H

#include "CoreTools.h"

class kmer_count {

     public:

     kmer_count( ) { }
     kmer_count( const int n ) : n(n) { }

     int N( ) const { return n; }

     int n;

};

#endif
