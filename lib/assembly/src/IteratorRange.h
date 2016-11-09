///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file IteratorRange.h
 * \author tsharpe
 * \date Jan 25, 2013
 *
 * \brief
 */
#ifndef ITERATORRANGE_H_
#define ITERATORRANGE_H_

#include "feudal/Iterator.h"
#include <ostream>

template <class IntT> // IntT is some integer-like thing
class RangeIterator
: public std::iterator<std::random_access_iterator_tag,IntT,
                          typename std::make_signed<IntT>::type>,
  public IteratorBase<RangeIterator<IntT>,IntT,
                          typename std::make_signed<IntT>::type>
{
public:
    typedef typename std::make_signed<IntT>::type difference_type;

    RangeIterator()=default;
    explicit RangeIterator( IntT pos )
    : IteratorBase<RangeIterator<IntT>,IntT,difference_type>(pos) {}

    // compiler-supplied copying, destruction are OK

    IntT operator*() const { return this->mPos; }

    IntT operator[]( difference_type diff ) const
    { return this->mPos+diff; }
};

template <class IntT>
RangeIterator<IntT> rangeItr( IntT pos )
{ return RangeIterator<IntT>(pos); }

template <class Itr>
class IteratorRange
{
public:
    IteratorRange() = default;
    IteratorRange( Itr const& beg, Itr const& end )
    : mBeg(beg), mEnd(end) {}

    // compiler-supplied copying, moving, and destructor are OK

    Itr begin() const { return mBeg; }
    Itr end() const { return mEnd; }

private:
    Itr mBeg;
    Itr mEnd;
};

template <class Itr>
IteratorRange<Itr> makeIteratorRange( Itr const& beg, Itr const& end )
{ return IteratorRange<Itr>(beg,end); }

template <class Itr>
class IteratorRangePrinter
{
public:
    IteratorRangePrinter( Itr const& beg, Itr const& end, char const* sep )
    : mBeg(beg), mEnd(end), mSep(sep) {}

    // compiler-supplied copying, moving, and destructor are OK

    friend std::ostream& operator<<( std::ostream& os,
                                        IteratorRangePrinter const& irp )
    { Itr itr = irp.mBeg, end = irp.mEnd;
      if ( itr != end )
      { os << *itr;
        while ( ++itr != end ) os << irp.mSep << *itr; }
      return os; }

private:
    Itr mBeg;
    Itr mEnd;
    char const* mSep;
};

template <class Itr>
IteratorRangePrinter<Itr> rangePrinter( Itr const& beg, Itr const& end,
                                                       char const* sep = ", " )
{ return IteratorRangePrinter<Itr>(beg,end,sep); }

#endif /* ITERATORRANGE_H_ */
