///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file PreCorrectAlt1.h
 * \author tsharpe
 * \date Nov 29, 2012
 *
 * \brief
 */
#ifndef PATHS_LONG_PRECORRECTALT1_H_
#define PATHS_LONG_PRECORRECTALT1_H_

#include "Basevector.h"

void precorrectAlt1( vecbvec* pReads, unsigned COVERAGE=50,
                        int VERBOSITY=0, unsigned NUM_THREADS=0 );

#endif /* PATHS_LONG_PRECORRECTALT1_H_ */
