///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef STLEXTENSIONS_H
#define STLEXTENSIONS_H

using namespace std;

#include <functional>
#include <vector>
#include <bitset>
#include <algorithm>
#include <iostream>

#include <math.h>
#include "feudal/BinaryStream.h"
#include "system/StaticAssert.h"

/// minimum<T> is a function object.
///
/// If f is an object of class minimum<T> and x and y are objects of
/// class T, then f(x,y) returns a copy of whichever of x or y is lower
/// as defined by T::operator<().
///
/// Example:
///
/// vector<int> V(10);
/// iota( V.begin(), V.end(), 1 );
///
/// copy( V.begin(), V.end(), ostream_iterator<int>( cout, " " ) );
/// cout << endl;
///
/// transform( V.begin(), V.end(), V.begin(), bind2nd(minimum<int>(), 5) );
///
/// copy( V.begin(), V.end(), ostream_iterator<int>( cout, " " ) );
/// cout << endl;
///
/// The preceding code produces the following output:
///
/// 1 2 3 4 5 6 7 8 9 10
/// 1 2 3 4 5 5 5 5 5 5

template <class T>
struct minimum : public binary_function<T, T, T>
{
  T operator() (const T& x, const T& y ) const { return ( x<y ? x : y ); }
};

/// maximum<T> is a function object.
///
/// If f is an object of class maximum<T> and x and y are objects of
/// class T, then f(x,y) returns a copy of whichever of x or y has the
/// highest value as defined by T::operator<().
///
/// Example:
///
/// vector<int> V(10);
/// iota( V.begin(), V.end(), 1 );
///
/// copy( V.begin(), V.end(), ostream_iterator<int>( cout, " " ) );
/// cout << endl;
///
/// transform( V.begin(), V.end(), V.begin(), bind2nd(maximum<int>(), 5) );
///
/// copy( V.begin(), V.end(), ostream_iterator<int>( cout, " " ) );
/// cout << endl;
///
/// The preceding code produces the following output:
///
/// 1 2 3 4 5 6 7 8 9 10
/// 5 5 5 5 5 6 7 8 9 10

template <class T>
struct maximum : public binary_function<T, T, T>
{
  T operator() (const T& x, const T& y ) const { return ( x<y ? y : x ); }
};

/// modulus<float> is a function object that implements the template modulus<T>
/// modulus<double> is a function object that implements the template modulus<T>
///
/// The default modulus<T> defines its operator() as
///
/// { return x % y; }
///
/// floats and doubles do not have a % operator defined.  They require specific
/// math functions to perform modulus operations.  These implementations are provided
/// below.

namespace std {
    template <>
    struct modulus<float> : public binary_function<float, float, float>
    {
        float operator() (const float& x, const float& y) const { return fmodf( x, y ); }
    };

    template <>
    struct modulus<double> : public binary_function<double, double, double>
    {
        double operator() (const double& x, const double& y) const { return fmod( x, y ); }
    };
}

/// address_of<T> and dereference<T> are function objects.
///
/// If f is an object of class address_of<T> and x is an object of
/// class T, then f(x) returns the address of x, i.e. a T*.
///
/// If f is an object of class dereference<T> and x is a pointer to
/// an object of class T, then f(x) returns a copy of *x.
///
/// Example:
///
/// vector<int> V1(10);
/// iota( V1.begin(), V1.end(), 1 );
///
/// copy( V1.begin(), V1.end(), ostream_iterator<int>( cout, " " ) );
/// cout << endl;
///
/// vector<int*> V_ptrs(10);
/// transform( V1.begin(), V1.end(), V_ptrs.begin(), address_of<int>() );
///
/// vector<int> V2(10);
/// transform( V_ptrs.begin(), V_ptrs.end(), dereference<int>() );
///
/// copy( V2.begin(), V2.end(), ostream_iterator<int>( cout, " " ) );
/// cout << endl;
///
/// The preceding code will produce the following output:
///
/// 1 2 3 4 5 6 7 8 9 10
/// 1 2 3 4 5 6 7 8 9 10
///
/// Note that address_of<T> is nearly always inexpensive, while
/// dereference<T> calls the copy constructor of T.

template <class T>
struct address_of : public unary_function<T, T*>
{
  const T* operator() (const T& x ) const { return &x; }
  T* operator() ( T& x ) const { return &x; }
};

template <class T>
struct dereference : public unary_function<T*, T>
{
  T operator() (const T* x ) const { return *x; }
};


// g++ 3.x defines select1st and select2nd in <ext/functional>

#if __GNUC__ > 2
#include <ext/functional>
using __gnu_cxx::select1st;
using __gnu_cxx::select2nd;
#endif

// g++ 3.x defines is_sorted in <ext/algorithm>
// g++ 3.x defines iota in <ext/numeric>

#if __GNUC__ > 2
#include <ext/algorithm>
using __gnu_cxx::is_sorted;
#include <ext/numeric>
using __gnu_cxx::iota;
#else
#if __GNUC__ <= 2
#include <algorithm>
#include <numeric>
#endif
#endif

///returns true if the order is strictly ascending.
template <class _ForwardIter>
bool is_sorted_strict(_ForwardIter __first, _ForwardIter __last)
{
  if (__first == __last)
    return true;

  _ForwardIter __next = __first;
  for (++__next; __next != __last; __first = __next, ++__next) {
    if (*__next <= *__first)
      return false;
  }

  return true;
}

///returns true if the order is strictly ascending.
template <class _ForwardIter, class _StrictWeakOrdering>
bool is_sorted_strict(_ForwardIter __first, _ForwardIter __last,
               _StrictWeakOrdering __comp)
{
  if (__first == __last)
    return true;

  _ForwardIter __next = __first;
  for (++__next; __next != __last; __first = __next, ++__next) {
    if (!(__comp(*__first, *__next) ) )
      return false;
  }

  return true;
}

/// sort_unique combines the STL algorithms sort and unique.  The uniquification
/// part of the algorithm uses iter_swap and operator< to do its dirtywork --
/// the same operations as are used by the sorting part of the algorithm --
/// unlike the STL version of unique, which uses operator== and assignment.
template <class RAItr>
RAItr sort_unique( RAItr start, RAItr const& stop )
{ using std::distance; using std::sort; using std::iter_swap;
  if ( distance(start,stop) < 2 ) return stop;
  sort(start,stop);
  RAItr dest(start);
  for ( ++start; start != stop; ++start )
    if ( *dest < *start && ++dest != start ) iter_swap(dest,start);
  return ++dest; }

template <class RAItr, class Comp>
RAItr sort_unique( RAItr start, RAItr const& stop, Comp comp )
{ using std::distance; using std::sort; using std::iter_swap;
  if ( distance(start,stop) < 2 ) return stop;
  sort(start,stop,comp);
  RAItr dest(start);
  for ( ++start; start != stop; ++start )
    if ( comp(*dest,*start) && ++dest != start ) iter_swap(dest,start);
  return ++dest; }

/// convenience function to uniquely sort a vector-like container that has
/// begin, end, and resize methods.
template <class V>
V& sort_unique( V& container )
{ typename V::iterator start(container.begin());
  container.resize(sort_unique(start,container.end())-start);
  return container; }

/// copyIf conditionally copies elements of a range to an output iterator.
/// copyIf is so tremendously useful, it's hard to understand why it's
/// not in the standard.  [Stroustrup says, "We just forgot."]
/// For notes on this implementation see
/// Effective STL (Meyers, 2001), Item 36, pp.154-156.
template<typename Itr, typename OItr, typename Pred>
inline OItr copyIf( Itr itr, Itr const& end, OItr dest, Pred pred )
{ while ( itr != end )
  { if ( pred(*itr) )
    { *dest = *itr; ++dest; }
    ++itr; }
  return dest; }


/// dereference_compare is a functor which allows comparison of
/// iterators using the operator< of the objects they point to.
///
/// usage example: suppose Foo::operator< exists.  Then:
///   list<Foo> foo_list;
///   vec< list<Foo>::iterator > foo_iters;
///   sort( foo_iters.begin(), foo_iters.end(),
///         dereferenced_compare<list<Foo>::iterator>() );

template<typename Iterator>
struct dereferenced_compare :
  public binary_function<Iterator, Iterator, bool>
{
  bool operator()(Iterator lhs, Iterator rhs) const {
    return( *lhs < *rhs );
  }
};

/*
   Template: update_min

   Update a running minimum: compare a value to the current minimum and update the minimum if the value is smaller.
 */
template <typename NumericType>
void update_min( NumericType& currentMin, const NumericType& val ) {
  if (val < currentMin )
    currentMin = val;
}

/*
   Template: update_max

   Update a running maximum: compare a value to the current maximum and update the maximum if the value is smaller.
 */
template <typename NumericType>
void update_max( NumericType& currentMax, const NumericType& val ) {
  if (val > currentMax )
    currentMax = val;
}

#define ForEach_Mut(x,containerType,c) \
    for ( containerType::iterator x = c.begin(); x != c.end(); x ++ )

#define ForEach(x,containerType,c) \
    for ( containerType::const_iterator x = c.begin() ; x != c.end(); x ++ )

template <class CONTAINER, class T> inline
  bool STLContains( const CONTAINER& c, const T& x ) {
  return c.find(x) != c.end();
}

/**
   Class: ShallowViewOf

   A shallow view of a given class: comparison operators are taken from the class,
   but copy constructor and assignment are the default.

   You would not normally create instances of ShallowViewOf<T>, but you can cast
   a (T *) to a (ShallowViewOf<T> *) before calling std::sort(), to avoid unnecessary
   deep copying.

   The template parameter T must be a model of STL concept LessThanComparable.
*/
template <class T>
class ShallowViewOf {
  char data_[sizeof(T)];
};

template <class T>
inline int operator<( const ShallowViewOf<T>& v1, const ShallowViewOf<T>& v2) {
  return ((T&)v1) < ((T&)v2);
}

template <class T>
inline int operator>( const ShallowViewOf<T>& v1, const ShallowViewOf<T>& v2) {
  return ((T&)v1) > ((T&)v2);
}

template <class T>
inline int operator<=( const ShallowViewOf<T>& v1, const ShallowViewOf<T>& v2) {
  return ((T&)v1) <= ((T&)v2);
}
template <class T>
inline int operator>=( const ShallowViewOf<T>& v1, const ShallowViewOf<T>& v2) {
  return ((T&)v1) >= ((T&)v2);
}
template <class T>
inline int operator==( const ShallowViewOf<T>& v1, const ShallowViewOf<T>& v2) {
  return ((T&)v1) == ((T&)v2);
}
template <class T>
inline int operator!=( const ShallowViewOf<T>& v1, const ShallowViewOf<T>& v2) {
  return ((T&)v1) != ((T&)v2);
}

/**
   Functor: order_ShallowView

   An ordering of ShallowViewOf<T> that just calls the given ordering on T.
*/
template <typename T, typename StrictWeakOrdering>
  struct order_ShallowView: public binary_function< T, T, bool > {
    private:
     StrictWeakOrdering orderOnT_;
    public:
     order_ShallowView<T, StrictWeakOrdering>( StrictWeakOrdering orderOnT ): orderOnT_(orderOnT) { }

    bool operator() ( const ShallowViewOf<T>& v1, const ShallowViewOf<T>& v2 ) const {
      return orderOnT_( (const T&)v1, (const T&)v2 );
    }
  };



/**
   Functor: cmp_functor

   Turn a C++ comparator function into a binary functor that always
   calls this specific comparator function.  Since the comparator function
   is a template parameter of this functor, the compiler may have an easier
   time inlining the function body than if you had passed around
   a pointer to the function itself.
*/
template<class T, bool myComparator(const T&, const T&) >
  struct cmp_functor: public binary_function< T, T, bool > {
    bool operator() ( const T& a, const T& b ) const {
      return myComparator(a,b);
    }
  };

#define COMPARE_BY2(T,F1,F2)  \
 friend bool operator< ( const T& v1, const T& v2 ) { \
   return v1.F1 < v2.F1 ? true : \
     ( v1.F1 > v2.F1 ? false : ( v1.F2 < v2.F2 ) );  \
 }

#define COMPARE_BY3(T,F1,F2,F3)  \
 friend bool operator< ( const T& v1, const T& v2 ) { \
   return  \
     (v1.F1) < (v2.F1) ? true : \
     (v1.F1) > (v2.F1) ? false: \
     (v1.F2) < (v2.F2) ? true : \
     (v1.F2) > (v2.F2) ? false : \
       (v1.F3) < (v2.F3);        \
 }

// Return the canonical version of the two pairs
// (a,b) and (b,a).
template <class _T1, class _T2>
inline pair< _T1, _T2 > CanonPair( const pair< _T1, _T2 >& p ) {
  return p.first <= p.second ? p : make_pair( p.second, p.first );
}

template <class _T1, class _T2>
inline pair< _T1, _T2 > CanonPair( const _T1& v1, const _T2& v2  ) {
  return CanonPair( make_pair( v1, v2 ) );
}

/// Class triple is just like class pair, but has three elements.

template <class _T1, class _T2, class _T3>
struct triple {
  typedef _T1 first_type;
  typedef _T2 second_type;
  typedef _T3 third_type;
  _T1 first;
  _T2 second;
  _T3 third;
  triple() : first(), second(), third() {}
  triple(const _T1& __a, const _T2& __b, const _T3& __c)
       : first(__a), second(__b), third(__c) {}
  template <class _U1, class _U2, class _U3>
  triple(const triple<_U1, _U2, _U3>& __p)
       : first(__p.first), second(__p.second), third(__p.third) {}

  friend int compare( triple const& t1, triple const& t2 )
  {
      int result = compare(t1.first,t2.first);
      if ( !result ) result = compare(t1.second,t2.second);
      if ( !result ) result = compare(t1.third,t2.third);
      return result;
  }


  void writeBinary( BinaryWriter& bw )
  { bw.write(first); bw.write(second); bw.write(third); }
  void readBinary( BinaryReader& br )
  { br.read(&first); br.read(&second); br.read(&third); }
  static size_t externalSizeof()
  { size_t sz1 = BinaryReader::externalSizeof(static_cast<_T1*>(0));
    size_t sz2 = BinaryReader::externalSizeof(static_cast<_T2*>(0));
    size_t sz3 = BinaryReader::externalSizeof(static_cast<_T3*>(0));
    return sz1 && sz2 && sz3 ? sz1+sz2+sz3 : 0UL; }
};

template <class ST1, class ST2, class ST3>
struct TripleSerializability
{ typedef SelfSerializable type; };

template <>
struct TripleSerializability<TriviallySerializable,
                                TriviallySerializable,
                                TriviallySerializable>
{ typedef TriviallySerializable type; };

template <class T1, class T2, class T3>
struct Serializability<triple<T1,T2,T3> >
{ typedef typename Serializability<T1>::type T1ST;
  typedef typename Serializability<T2>::type T2ST;
  typedef typename Serializability<T3>::type T3ST;
  typedef typename TripleSerializability<T1ST,T2ST,T3ST>::type type; };


template <class _T1, class _T2, class _T3>
inline bool operator==(const triple<_T1, _T2, _T3>& __x,
     const triple<_T1, _T2, _T3>& __y)
{
  return __x.first == __y.first && __x.second == __y.second
       && __x.third == __y.third;
}

template <class _T1, class _T2, class _T3>
inline bool operator<(const triple<_T1, _T2, _T3>& __x,
     const triple<_T1, _T2, _T3>& __y)
{
     if ( __x.first < __y.first ) return true;
     if ( __x.first > __y.first ) return false;
     if ( __x.second < __y.second ) return true;
     if ( __x.second > __y.second ) return false;
     if ( __x.third < __y.third ) return true;
     return false;

}

template <class _T1, class _T2, class _T3>
inline bool operator!=(const triple<_T1, _T2, _T3>& __x,
     const triple<_T1, _T2, _T3>& __y) {
  return !(__x == __y);
}

template <class _T1, class _T2, class _T3>
inline bool operator>(const triple<_T1, _T2, _T3>& __x,
     const triple<_T1, _T2, _T3>& __y) {
  return __y < __x;
}

template <class _T1, class _T2, class _T3>
inline bool operator<=(const triple<_T1, _T2, _T3>& __x,
     const triple<_T1, _T2, _T3>& __y) {
  return !(__y < __x);
}

template <class _T1, class _T2, class _T3>
inline bool operator>=(const triple<_T1, _T2, _T3>& __x,
     const triple<_T1, _T2, _T3>& __y) {
  return !(__x < __y);
}

template <class _T1, class _T2, class _T3>
inline triple<_T1, _T2, _T3> make_triple(const _T1& __x, const _T2& __y, const _T3& __z)
{
  return triple<_T1, _T2, _T3>(__x, __y, __z);
}

// LowerBound1: LowerBound on x.first.

template<class T1, class T2> int64_t LowerBound1(
     const vector< pair<T1,T2> >& x, const T1& t )
{    int64_t count = x.size( ), it, step, first = 0;
     while ( count > 0 )
     {    it = first;
          step = count/2;
          it += step;
          if ( x[it].first < t )
          {    first = ++it;
               count -= step + 1;    }
          else count = step;    }
     return first;    }

template<class T1, class T2> int64_t UpperBound1(
     const vector< pair<T1,T2> >& x, const T1& t )
{    int64_t count = x.size( ), it, step, first = 0;
     while ( count > 0 )
     {    it = first;
          step = count/2;
          it += step;
          if ( !( t < x[it].first ) )
          {    first = ++it;
               count -= step + 1;    }
          else count = step;    }
     return first;    }


template<class T1, class T2, class T3> int64_t LowerBound1(
     const vector< triple<T1,T2,T3> >& x, const T1& t )
{    int64_t count = x.size( ), it, step, first = 0;
     while ( count > 0 )
     {    it = first;
          step = count/2;
          it += step;
          if ( x[it].first < t )
          {    first = ++it;
               count -= step + 1;    }
          else count = step;    }
     return first;    }

template<class T1, class T2, class T3> int64_t UpperBound1(
     const vector< triple<T1,T2,T3> >& x, const T1& t )
{    int64_t count = x.size( ), it, step, first = 0;
     while ( count > 0 )
     {    it = first;
          step = count/2;
          it += step;
          if ( !( t < x[it].first ) )
          {    first = ++it;
               count -= step + 1;    }
          else count = step;    }
     return first;    }



/// Class quad is just like class pair, but has four elements.

template <class _T1, class _T2, class _T3, class _T4>
struct quad {
  typedef _T1 first_type;
  typedef _T2 second_type;
  typedef _T3 third_type;
  typedef _T4 fourth_type;
  _T1 first;
  _T2 second;
  _T3 third;
  _T4 fourth;
  quad() : first(), second(), third(), fourth() {}
  quad(const _T1& __a, const _T2& __b, const _T3& __c, const _T4& __d)
       : first(__a), second(__b), third(__c), fourth(__d) {}
  template <class _U1, class _U2, class _U3, class _U4>
  quad(const quad<_U1, _U2, _U3, _U4>& __p)
       : first(__p.first), second(__p.second), third(__p.third), fourth(__p.fourth) {}

  friend int compare( quad const& t1, quad const& t2 )
  {
      int result = compare(t1.first,t2.first);
      if ( !result ) result = compare(t1.second,t2.second);
      if ( !result ) result = compare(t1.third,t2.third);
      if ( !result ) result = compare(t1.fourth,t2.fourth);
      return result;
  }



  void writeBinary( BinaryWriter& bw )
  { bw.write(first); bw.write(second); bw.write(third); bw.write(fourth); }
  void readBinary( BinaryReader& br )
  { br.read(&first); br.read(&second); br.read(&third), br.read(&fourth); }
  static size_t externalSizeof()
  { size_t sz1 = BinaryReader::externalSizeof(static_cast<_T1*>(0));
    size_t sz2 = BinaryReader::externalSizeof(static_cast<_T2*>(0));
    size_t sz3 = BinaryReader::externalSizeof(static_cast<_T3*>(0));
    size_t sz4 = BinaryReader::externalSizeof(static_cast<_T4*>(0));
    return sz1 && sz2 && sz3 && sz4 ? sz1+sz2+sz3+sz4 : 0UL; }
};

template <class ST1, class ST2, class ST3, class ST4>
struct QuadSerializability
{ typedef SelfSerializable type; };

template <>
struct QuadSerializability<TriviallySerializable,
                                TriviallySerializable,
                                TriviallySerializable,
                                TriviallySerializable>
{ typedef TriviallySerializable type; };

template <class T1, class T2, class T3, class T4>
struct Serializability<quad<T1,T2,T3,T4> >
{ typedef typename Serializability<T1>::type T1ST;
  typedef typename Serializability<T2>::type T2ST;
  typedef typename Serializability<T3>::type T3ST;
  typedef typename Serializability<T4>::type T4ST;
  typedef typename QuadSerializability<T1ST,T2ST,T3ST,T4ST>::type type; };


template <class _T1, class _T2, class _T3, class _T4>
inline bool operator==(const quad<_T1, _T2, _T3, _T4>& __x,
     const quad<_T1, _T2, _T3, _T4>& __y)
{
  return __x.first == __y.first && __x.second == __y.second
       && __x.third == __y.third && __x.fourth == __y.fourth;
}

template <class _T1, class _T2, class _T3, class _T4>
inline bool operator<(const quad<_T1, _T2, _T3, _T4>& __x,
     const quad<_T1, _T2, _T3, _T4>& __y)
{
     if ( __x.first < __y.first ) return true;
     if ( __x.first > __y.first ) return false;
     if ( __x.second < __y.second ) return true;
     if ( __x.second > __y.second ) return false;
     if ( __x.third < __y.third ) return true;
     if ( __x.third > __y.third ) return false;
     if ( __x.fourth < __y.fourth ) return true;
     return false;

}

template <class _T1, class _T2, class _T3, class _T4>
inline bool operator!=(const quad<_T1, _T2, _T3, _T4>& __x,
     const quad<_T1, _T2, _T3, _T4>& __y) {
  return !(__x == __y);
}

template <class _T1, class _T2, class _T3, class _T4>
inline bool operator>(const quad<_T1, _T2, _T3, _T4>& __x,
     const quad<_T1, _T2, _T3, _T4>& __y) {
  return __y < __x;
}

template <class _T1, class _T2, class _T3, class _T4>
inline bool operator<=(const quad<_T1, _T2, _T3, _T4>& __x,
     const quad<_T1, _T2, _T3, _T4>& __y) {
  return !(__y < __x);
}

template <class _T1, class _T2, class _T3, class _T4>
inline bool operator>=(const quad<_T1, _T2, _T3, _T4>& __x,
     const quad<_T1, _T2, _T3, _T4>& __y) {
  return !(__x < __y);
}

template <class _T1, class _T2, class _T3, class _T4>
inline quad<_T1, _T2, _T3, _T4> make_quad(const _T1& __x, const _T2& __y, const _T3& __z, const _T4& __w)
{
  return quad<_T1, _T2, _T3, _T4>(__x, __y, __z, __w);
}


template<class T1, class T2, class T3, class T4> int64_t LowerBound1(
     const vector< quad<T1,T2,T3,T4> >& x, const T1& t )
{    int64_t count = x.size( ), it, step, first = 0;
     while ( count > 0 )
     {    it = first;
          step = count/2;
          it += step;
          if ( x[it].first < t )
          {    first = ++it;
               count -= step + 1;    }
          else count = step;    }
     return first;    }

template<class T1, class T2, class T3, class T4> int64_t UpperBound1(
     const vector< quad<T1,T2,T3,T4> >& x, const T1& t )
{    int64_t count = x.size( ), it, step, first = 0;
     while ( count > 0 )
     {    it = first;
          step = count/2;
          it += step;
          if ( !( t < x[it].first ) )
          {    first = ++it;
               count -= step + 1;    }
          else count = step;    }
     return first;    }


// BinPosition1.  Using .first, return the position of an element in a sorted
// vector, else -1.  If the element appears more than once, the position of
// one of its instances is returned.

template<class T1, class T2, class A, class U>
int64_t BinPosition1( const vector<pair<T1,T2>,A>& v, const U& x1 )
{    if ( v.size( ) == 0 ) return -1;
     T1 const& x(x1);
     size_t first = 0, last = v.size( ) - 1, next;
     while(1)
     {    if (first == last) return ( !(x < v[last].first) && !(v[last].first < x) ) ? last : -1;
          next = first + (last - first) / 2;
          if ( x < v[next].first ) last = next;
          else if ( v[next].first < x ) first = next + 1;
          else return next;    }    }

template<class T1, class T2, class T3, class A, class U>
int64_t BinPosition1( const vector<triple<T1,T2,T3>,A>& v, const U& x1 )
{    if ( v.size( ) == 0 ) return -1;
     T1 const& x(x1);
     size_t first = 0, last = v.size( ) - 1, next;
     while(1)
     {    if (first == last) return ( !(x < v[last].first) && !(v[last].first < x) ) ? last : -1;
          next = first + (last - first) / 2;
          if ( x < v[next].first ) last = next;
          else if ( v[next].first < x ) first = next + 1;
          else return next;    }    }

template <class _T1, class _T2, class _T3>
  inline ostream& operator<< ( ostream& out, const triple<_T1, _T2, _T3>& __x ) {
  out << "(" << __x.first << ", " << __x.second << ", " << __x.third << ")";
  return out;
}

template < size_t N >
bool operator< ( const bitset< N >& b1, const bitset< N >& b2 ) {
  for ( size_t i = 0; i < N; i++ )
    if ( !b1.test( i )  &&  b2.test( i ) )
      return true;
  return false;
}

template < size_t N >
bool operator<= ( const bitset< N >& b1, const bitset< N >& b2 ) {
  return b1 < b2  ||  b1 == b2;
}

template < size_t N >
bool operator> ( const bitset< N >& b1, const bitset< N >& b2 ) {
  return !( b1 <= b2 );
}

template < size_t N >
bool operator>= ( const bitset< N >& b1, const bitset< N >& b2 ) {
  return !( b1 < b2 );
}

template <class T>
class LtBySize : public std::binary_function<T,T,bool>
{    public:
     bool operator()( T const& t1, T const& t2 )
     { return t1.size() < t2.size(); }
};

// use this via a printSeqhelper function below
template <class Itr>
class SeqPrinter
{
public:
    SeqPrinter( Itr const& beg, Itr const& end, char const* delim, bool exp )
    : mBeg(beg), mEnd(end), mDelim(delim), mExp(exp) {}

    template <class Delim>
    SeqPrinter( Itr const& beg, Itr const& end, Delim const& delim, bool exp )
    : mBeg(beg), mEnd(end), mDelim(delim.c_str()), mExp(exp) {}

    friend std::ostream& operator<<(std::ostream& os, SeqPrinter const& sp)
    {
      Itr itr(sp.mBeg), end(sp.mEnd);
      if ( !sp.mExp )
      {
      while ( itr != end )
      { os << *itr; if ( ++itr != end ) os << sp.mDelim; }
      return os;
      }
      else
      {    while ( itr != end )
           {    os << *itr;
                int count = 1;
                while( ++itr != end )
                {    itr--;
                     if ( *itr != *(++itr) ) break;
                     count++;    }
                --itr;
                if ( count > 1 ) os << "^" << count;
                if ( ++itr != end ) os << sp.mDelim;     }
           return os;     }
      }

private:
    Itr mBeg;
    Itr mEnd;
    std::string mDelim;
    bool mExp;
};

// Itr is an input iterator
template <class Itr>
SeqPrinter<Itr> printSeq( Itr const& beg, Itr const& end,
                                        char const* delim = "," )
{ return SeqPrinter<Itr>(beg,end,delim,false); }

// Itr is an input iterator, Delim has a c_str method (like string or String)
template <class Itr, class Delim>
SeqPrinter<Itr> printSeq( Itr const& beg, Itr const& end,
                              Delim const& delim,
                              char const*(Delim::*)() const=&Delim::c_str )
{ return SeqPrinter<Itr>(beg,end,delim,false); }

// Cont has a begin and end method, and typedefs a const_iterator type
template <class Cont>
SeqPrinter<typename Cont::const_iterator> printSeq( Cont const& cont,
                                                     char const* delim = "," )
{ typedef typename Cont::const_iterator Itr;
  return SeqPrinter<Itr>(cont.begin(),cont.end(),delim,false); }

// Cont has a begin and end method, and typedefs a const_iterator type
// Itr is an input iterator, Delim has a c_str method (like string or String)
template <class Cont, class Delim>
SeqPrinter<typename Cont::const_iterator> printSeq( Cont const& cont,
                                                           Delim const& delim,
                               char const*(Delim::*)() const=&Delim::c_str )
{ typedef typename Cont::const_iterator Itr;
  return SeqPrinter<Itr>(cont.begin(),cont.end(),delim,false); }

// Exp versions:

template <class Itr> SeqPrinter<Itr> printSeqExp( Itr const& beg, Itr const& end,
                                        char const* delim = "," )
{ return SeqPrinter<Itr>(beg,end,delim,true); }

template <class Itr, class Delim> SeqPrinter<Itr> printSeqExp( Itr const& beg,
     Itr const& end, Delim const& delim,
                              char const*(Delim::*)() const=&Delim::c_str )
{ return SeqPrinter<Itr>(beg,end,delim,true); }

template <class Cont> SeqPrinter<typename Cont::const_iterator> printSeqExp(
     Cont const& cont, char const* delim = "," )
{ typedef typename Cont::const_iterator Itr;
  return SeqPrinter<Itr>(cont.begin(),cont.end(),delim,true); }

template <class Cont, class Delim> SeqPrinter<typename Cont::const_iterator>
     printSeqExp( Cont const& cont, Delim const& delim,
                               char const*(Delim::*)() const=&Delim::c_str )
{ typedef typename Cont::const_iterator Itr;
  return SeqPrinter<Itr>(cont.begin(),cont.end(),delim,true); }

template <class C>
void EraseValue( C& container, typename C::value_type const& value )
{ container.erase(
        std::remove(container.begin(),container.end(),value),
        container.end()); }

template<class X> bool IsUnique( const X& x1, const X& x2, const X& x3 )
{    if ( x1 == x2 || x1 == x3 ) return false;
     if ( x2 == x3 ) return false;
     return true;    }

template<class X> bool IsUnique( const X& x1, const X& x2, const X& x3, const X& x4 )
{    if ( x1 == x2 || x1 == x3 || x1 == x4 ) return false;
     if ( x2 == x3 || x2 == x4 ) return false;
     if ( x3 == x4 ) return false;
     return true;    }

#endif
