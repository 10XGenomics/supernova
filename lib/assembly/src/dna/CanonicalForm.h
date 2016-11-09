///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file CanonicalForm.h
 * \author tsharpe
 * \date Aug 30, 2012
 *
 * \brief
 */
#ifndef DNA_CANONICALFORM_H_
#define DNA_CANONICALFORM_H_

#include <iterator>
#include <ostream>

enum class CanonicalForm : unsigned char { FWD, REV, PALINDROME };

inline CanonicalForm complement( CanonicalForm form )
{ switch ( form )
  { case CanonicalForm::FWD: form = CanonicalForm::REV; break;
    case CanonicalForm::REV: form = CanonicalForm::FWD; break;
    case CanonicalForm::PALINDROME: break; }
  return form; }

/// Get CanonicalForm when the length can only be determined at run-time
/// (like for a bvec).
/// The Itr type must have base codes as its value_type.
/// You'll get a compiler error if Itr isn't a random access iterator.
template <class Itr>
inline CanonicalForm getCanonicalForm( Itr beg, Itr end )
{ using std::distance;
  size_t len = distance(beg,end);
  if ( len & 1 )
    return beg[len/2] & 2 ? CanonicalForm::REV : CanonicalForm::FWD;
  while ( beg != end )
  { auto f = *beg;
    auto r = *--end ^ 3;
    if ( f < r ) return CanonicalForm::FWD;
    if ( r < f ) return CanonicalForm::REV;
    ++beg; }
  return CanonicalForm::PALINDROME; }

inline std::ostream& operator<<( std::ostream& os, CanonicalForm form )
{ return os << "+-|"[int(form)]; }

template <unsigned K>
class CF
{
public:
    /// The Itr type must have base codes as its value_type.
    /// You'll get a compiler error if Itr isn't a random access iterator.
    template <class Itr>
    static CanonicalForm getForm( Itr beg )
    { if ( K&1 ) return beg[K/2] & 2 ? CanonicalForm::REV : CanonicalForm::FWD;
      Itr end(beg+K);
      while ( beg != end )
      { auto f = *beg;
        auto r = *--end ^ 3;
        if ( f < r ) return CanonicalForm::FWD;
        if ( r < f ) return CanonicalForm::REV;
        ++beg; }
      return CanonicalForm::PALINDROME; }

    template <class Itr>
    static bool isFwd( Itr beg ) { return getForm(beg) == CanonicalForm::FWD; }

    template <class Itr>
    static bool isRev( Itr beg ) { return getForm(beg) == CanonicalForm::REV; }

    template <class Itr>
    static bool isPalindrome( Itr beg )
    { return (K&1) ? false : getForm(beg) == CanonicalForm::PALINDROME; }

    /// You've got two kmers, and you know they're either identical, or one is
    /// the RC of the other.  This tells you which case is true.  Note that
    /// you'll get unpredictable garbage as an answer if the kmers are neither
    /// equal nor RC.
    /// The iterators must be random access.
    template <class Itr1, class Itr2>
    static bool isRC( Itr1 itr1, Itr2 itr2 )
    { if ( K&1 ) return itr1[K/2] != itr2[K/2];
      Itr1 end(itr1+K);
      while ( itr1 != end )
      { if ( *itr1 != (*--end ^3) ) return *itr1 != *itr2;
        ++itr1, ++itr2; }
      return false; }
};

#endif /* DNA_CANONICALFORM_H_ */
