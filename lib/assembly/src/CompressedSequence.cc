///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CompressedSequence.h"

void CompressedSequence::ReverseComplement()
{
    iterator head = begin();
    iterator tail = end();
    while (head != tail)
    {
        value_type tmp = GeneralizedBase::fromBits(*--tail).complement().bits();
        if (head == tail)
        {
            head.set(tmp);
            break;
        }
        tail.set(GeneralizedBase::fromBits(*head).complement().bits());
        head.set(tmp);
        ++head;
    }
}

void CompressedSequence::asBasevector( basevector &bv, bool allow_x ) const
{
    bv.clear().reserve(size());
    for ( const_iterator itr(begin()), stop(end()); itr != stop; ++itr ) {
	if ( ! allow_x )
	    bv.push_back(GeneralizedBase::bits2Val(*itr));
	else {
	    // Treat 'X' as 'N' to prevent crash
	    GeneralizedBase const& gb = GeneralizedBase::fromBits(*itr);
	    bv.push_back( gb ==  GeneralizedBase::X ? GeneralizedBase::N.random() : gb.random() );
	}
    }
}

void CompressedSequence::getAmbBases( bitvector &bitv ) const
{
    bitv.clear().reserve(size());
    for ( const_iterator itr(begin()), stop(end()); itr != stop; ++itr )
        bitv.push_back(GeneralizedBase::bits2Ambig(*itr));
}

void CompressedSequence::assignChars( char const* begin, char const* end )
{
    clear().reserve(end-begin);
    while ( begin != end )
    {
        char chr = *begin++;
        if ( chr != '*' )
        {
            if ( GeneralizedBase::isGeneralizedBase(chr) )
                push_back(GeneralizedBase::char2Bits(chr));
            else
                push_back(GeneralizedBase::N.bits());
        }
    }
}

#include "feudal/OuterVecDefs.h"
template class OuterVec<CompressedSequence>;
