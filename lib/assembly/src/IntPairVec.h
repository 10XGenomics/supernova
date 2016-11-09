///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file IntPairVec.h
 * \author tsharpe
 * \date Jan 24, 2013
 *
 * \brief
 */
#ifndef INTPAIRVEC_H_
#define INTPAIRVEC_H_

#include "feudal/MasterVec.h"
#include "feudal/SerfVec.h"
#include <utility>

typedef std::pair<int,int> IntPair;
typedef SerfVec<IntPair> IntPairVec;
typedef MasterVec<IntPairVec> VecIntPairVec;
extern template class OuterVec<IntPairVec>;

#endif /* INTPAIRVEC_H_ */
