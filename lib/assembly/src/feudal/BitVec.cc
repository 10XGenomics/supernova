///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "feudal/BitVec.h"

BitVec& BitVec::operator|=( BitVec const& bv )
{
    ForceAssertEq(size(),bv.size());

    const_pointer itr2 = bv.data();
    for ( pointer itr=data(),end=dataEnd(); itr != end; ++itr,++itr2 )
        *itr = *itr | *itr2;

    return *this;
}

BitVec& BitVec::operator&=( BitVec const& bv )
{
    ForceAssertEq(size(),bv.size());

    const_pointer itr2 = bv.data();
    for ( pointer itr=data(),end=dataEnd(); itr != end; ++itr,++itr2 )
        *itr = *itr & *itr2;

    return *this;
}

BitVec& BitVec::operator^=( BitVec const& bv )
{
    ForceAssertEq(size(),bv.size());

    const_pointer itr2 = bv.data();
    for ( pointer itr=data(),end=dataEnd(); itr != end; ++itr,++itr2 )
        *itr = *itr ^ *itr2;

    return *this;
}

BitVec& BitVec::invert()
{
    for ( pointer itr=data(),end=dataEnd(); itr != end; ++itr )
        *itr = ~*itr;
    return *this;
}

BitVec& BitVec::Zero()
{
    for ( pointer itr=data(),end=dataEnd(); itr != end; ++itr )
        *itr = 0;
    return *this;
}

BitVec& BitVec::SetToSubOf( const BitVec& src, size_type start, size_type len )
{
    AssertLe(start,src.size());
    AssertLe(len,src.size()-start);

    if ( size() < len )
        resize(len);

    const_iterator head(src.begin(start));
    const_iterator tail(src.begin(start + len));
    for ( iterator dst(begin()); head != tail; ++head, ++dst )
        dst.set(*head);

    if ( size() > len )
        resize(len);

    return *this;
}

void BitVec::PrintFastaStyle( ostream& out, const String& id ) const
{
    out << '>' << id;
    for ( size_type i = 0; i < size(); ++i )
    {
        out << (i % 40 ? ' ' : '\n');
        out << ((*this)[i] ? '1' : '0');
    }
    out << '\n';
}

float Coverage( const vecbitvector& v )
{
    vecbitvector::size_type covered = 0;
    vecbitvector::size_type total = 0;
    for ( vecbitvector::size_type i = 0; i < v.size(); i++ )
    {
        total += v[i].size();
        covered += v[i].Sum();
    }
    AssertGt(total,0u);
    return static_cast<float> (covered) / total;
}

// Logical Operators

vecbitvector& operator|=( vecbitvector& v1, const vecbitvector& v2 )
{
    AssertEq(v1.size(),v2.size());
    vecbitvector::const_iterator src(v2.begin());
    vecbitvector::iterator stop(v1.end());
    for ( vecbitvector::iterator itr(v1.begin()); itr != stop; ++itr, ++src )
        *itr |= *src;
    return v1;
}

vecbitvector& operator &= (vecbitvector& v1, const vecbitvector& v2)
{
    AssertEq(v1.size(),v2.size());
    vecbitvector::const_iterator src(v2.begin());
    vecbitvector::iterator stop(v1.end());
    for ( vecbitvector::iterator itr(v1.begin()); itr != stop; ++itr, ++src )
        *itr &= *src;
    return v1;
}

vecbitvector& operator ^= (vecbitvector& v1, const vecbitvector& v2)
{
    AssertEq(v1.size(),v2.size());
    vecbitvector::const_iterator src(v2.begin());
    vecbitvector::iterator stop(v1.end());
    for ( vecbitvector::iterator itr(v1.begin()); itr != stop; ++itr, ++src )
        *itr ^= *src;
    return v1;
}

#include "feudal/OuterVecDefs.h"
template class OuterVec<BitVec>;
