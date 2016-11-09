///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file MultipleAligner.cc
 * \author tsharpe
 * \date May 8, 2012
 *
 * \brief
 */

#include "paths/long/ultra/MultipleAligner.h"

char const gvec::gCharSet[] = {'A','C','G','T','-'};
gvec::value_type const Scorer::NVALS;
gvec::value_type const MultipleAligner::ConsScore::NVALS;
