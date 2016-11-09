///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * \file Charvector.h
 * \author jbutler
 * \date June 25, 2004
 *
 * \brief Feudal vectors of chars, uchars, and Bools.
 */
#ifndef CHARVECTOR_H_
#define CHARVECTOR_H_

#include "feudal/MasterVec.h"
#include "feudal/SerfVec.h"

typedef SerfVec<char> CharVec;
typedef MasterVec< CharVec > VecCharVec;
extern template class OuterVec<CharVec>;

typedef SerfVec<unsigned char> UCharVec;
typedef MasterVec< UCharVec > VecUCharVec;
extern template class OuterVec<UCharVec>;

typedef UCharVec Boolvector;
typedef VecUCharVec vecBoolvector;

// Remove newlines from in and save the result in out.
void StripNewlines( const CharVec &in, CharVec &out );

#endif /* CHARVECTOR_H_ */
