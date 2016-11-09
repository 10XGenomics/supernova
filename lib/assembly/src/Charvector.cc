///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * \file Charvector.cc
 * \author tsharpe
 * \date Sep 3, 2009
 *
 * \brief
 */
#include "Charvector.h"

void StripNewlines( const CharVec &in, CharVec &out )
{
    out.reserve(in.size()).clear();
    CharVec::const_iterator end(in.end());
    for ( CharVec::const_iterator itr(in.begin()); itr != end; ++itr )
    {
        char val = *itr;
        if ( val != '\n' )
            out.push_back(val);
    }
}

#include "feudal/SmallVecDefs.h"
#include "feudal/OuterVecDefs.h"
template class SmallVec< char, MempoolAllocator<char> >;
template class OuterVec<CharVec>;
template class SmallVec< unsigned char, MempoolAllocator<unsigned char> >;
template class OuterVec<UCharVec>;
