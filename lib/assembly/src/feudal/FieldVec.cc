///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file FieldVec.cc
 * \author tsharpe
 * \date Mar 22, 2012
 *
 * \brief
 */

#include "feudal/FieldVecDefs.h"
template class FieldVec< 1, MempoolAllocator<unsigned char> >;
template class FieldVec< 2, MempoolAllocator<unsigned char> >;
template class FieldVec< 4, MempoolAllocator<unsigned char> >;
