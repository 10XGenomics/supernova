///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/// Wraps an STL vector, adds asserts for debugging and useful methods.
/// \class Vec
/// Vec.h defines class vec, which wraps the STL class vector, in such a way
/// that if compiled with NDEBUG on, it is (or should be) the same as a vector
/// (with some added functionality -- see below), but otherwise does run-time 
/// checking of each vector reference to make sure it is in range.
///
/// In debug mode, resize and reserve will fail if you (in effect) ask for more
/// than 100GB (or 1.5 GB on 32-bit systems).
///
/// PLEASE help keep Vec.h tidy by avoiding placing functions that operate on
/// vec in this file - consider using VecUtilities.h instead.
///
/// See Also: VecUtilities.h

#ifndef VEC_H
#define VEC_H

#include <unistd.h>

#include <algorithm>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <memory>
#include <strstream>
#include <vector>
#include <functional>
#include <scoped_allocator>

#include "String.h"
#include "feudal/TrackingAllocator.h"
#include "system/Assert.h"
#include "system/StaticAssert.h"
#include "STLExtensions.h"
#include "system/Types.h"
#include "system/TraceVal.h"
#include "feudal/BinaryStream.h"
#include "system/SortInPlace.h"
#include "system/System.h"
#include "Compare.h"
#include "system/file/FileReader.h"

/////////////////////////////////////////////////////////////////////////////
//
//  vec Class Declaration and Template Definitions
//

template <class T>
using StdVec = std::vector<T,typename DefaultAllocator<T>::type>;

template <class T, class A=typename DefaultAllocator<T>::type> class vec
  : public std::vector<T,A>
{
    typedef std::vector<T,A> BaseT;
 public:
    typedef typename BaseT::size_type size_type;

// ===========================================================================
//
// CONSTRUCTORS
//
// ===========================================================================

  explicit vec( A const& a=A() ) : BaseT(a) {}
  vec( const initializer_list<T>& x, A const& a=A() ) : BaseT(x,a) {}
  explicit vec( size_type n ) : BaseT(n) {}
  vec( size_type n, const T& value, A const& a=A() )
  : BaseT(ValidatedSize(n),value,a) {}

  enum ConstructorBehavior { IDENTITY };
  vec( size_type n, ConstructorBehavior constructor_type, A const& a=A() )
  : BaseT(ValidatedSize(n))
  { for ( size_type i = 0; i < n; i++ ) (*this)[i] = i;    }

  template <class A1>
  vec( const std::vector<T,A1>& v, A const&a=A() )
  : BaseT(v.begin(),v.end(),a) {}

  vec( vec const& v ) : BaseT(v.begin(),v.end(),A()) {}
  vec( vec const& v, A const& a ) : BaseT(v.begin(),v.end(),a) {}

  template<class Itr>
  vec( Itr const& first, Itr const& last, A const& a=A() )
  : BaseT(first,last,a) {}

  template <class A1>
  vec& operator=( std::vector<T,A1> const& v )
  { BaseT::assign(v.begin(),v.end()); return *this; }

  vec& operator=( vec const& v ) = default;

  ///Asserts index within bounds.
  typename BaseT::reference operator[]( size_type i ) {
    AssertLt( i,  BaseT::size() ); // Asserts index within bounds
    return BaseT::operator[](i);   // ... and returns the element
  }
  
  ///Asserts index within bounds.
  typename BaseT::const_reference operator[](size_type i) const {
    AssertLt( i, BaseT::size() );  // Asserts index within bounds
    return BaseT::operator[](i);   // ... and returns the element
  }

  void resize( size_type i, T const& c = T() )
  {  BaseT::resize( ValidatedSize(i), c ); }

  void reserve( size_type i ) { BaseT::reserve(ValidatedSize(i)); }

  typename BaseT::reference front( ) {
    AssertGt( BaseT::size( ), 0u ); // Asserts index within bounds
    return BaseT::front( );   // ... and returns the element
  }

  typename BaseT::const_reference front( ) const {
    AssertGt( BaseT::size( ), 0u ); // Asserts index within bounds
    return BaseT::front( );   // ... and returns the element
  }

  typename BaseT::reference back( ) {
    AssertGt( BaseT::size( ), 0u ); // Asserts index within bounds
    return BaseT::back( );   // ... and returns the element
  }

  typename BaseT::const_reference back( ) const {
    AssertGt( BaseT::size( ), 0u ); // Asserts index within bounds
    return BaseT::back( );   // ... and returns the element
  }

// ===========================================================================
//
// FUNCTIONS TO PUSH ELEMENTS ONTO VECTORS
//
// ===========================================================================

  /// push_front (insert item before first element of vector)

  vec& push_front( const T& t1 )
  {    BaseT::insert( BaseT::begin( ), t1 ); return *this;   }

  /// Generalized push_back, allowing up to 12 items pushed back at a time.

  vec& push_back( const T& t1 )
  {    BaseT::push_back(t1); return *this;   }
  vec& push_back( const T& t1, const T& t2 )
  {    return push_back(t1).push_back(t2);    }
  vec& push_back( const T& t1, const T& t2, const T& t3 )
  {    return push_back(t1).push_back(t2).push_back(t3);    }
  vec& push_back( const T& t1, const T& t2, const T& t3, const T& t4 )
  {    return push_back(t1).push_back(t2).push_back(t3).push_back(t4);    }
  vec& push_back( const T& t1, const T& t2, const T& t3, const T& t4, const T& t5 )
  {    return push_back(t1).push_back(t2).push_back(t3).push_back(t4).push_back(t5); }
  vec& push_back( const T& t1, const T& t2, const T& t3, const T& t4, const T& t5,
       const T& t6 )
  {    return push_back( t1, t2, t3, t4, t5 ).push_back( t6 );     }
  vec& push_back( const T& t1, const T& t2, const T& t3, const T& t4, const T& t5,
       const T& t6, const T& t7 )
  {    return push_back( t1, t2, t3, t4, t5 ).push_back( t6, t7 );     }
  vec& push_back( const T& t1, const T& t2, const T& t3, const T& t4, const T& t5,
       const T& t6, const T& t7, const T& t8 )
  {    return push_back( t1, t2, t3, t4, t5 ).push_back( t6, t7, t8 );     }
  vec& push_back( const T& t1, const T& t2, const T& t3, const T& t4, const T& t5,
       const T& t6, const T& t7, const T& t8, const T& t9 )
  {    return push_back( t1, t2, t3, t4, t5 ).push_back( t6, t7, t8, t9 );     }
  vec& push_back( const T& t1, const T& t2, const T& t3, const T& t4, const T& t5,
       const T& t6, const T& t7, const T& t8, const T& t9, const T& t10 )
  {    return push_back( t1, t2, t3, t4, t5 ).push_back( t6, t7, t8, t9, t10 ); }
  vec& push_back( const T& t1, const T& t2, const T& t3, const T& t4, const T& t5,
       const T& t6, const T& t7, const T& t8, const T& t9, const T& t10,
       const T& t11 )
  {    return push_back( t1, t2, t3, t4, t5 ).push_back( t6, t7, t8, t9, t10, t11 );     }
  vec& push_back( const T& t1, const T& t2, const T& t3, const T& t4, const T& t5,
       const T& t6, const T& t7, const T& t8, const T& t9, const T& t10,
       const T& t11, const T& t12 )
  {    return push_back( t1, t2, t3, t4, t5 ).push_back( t6, t7, t8, t9, t10, t11, t12 );     }

  // push: construct object from arbitrary arguments, then push it back.

  template <class... Args>
  vec& push( Args&&... args )
  { BaseT::emplace_back(std::forward<Args>(args)...); return *this; }

  void push_back_copies( const T& t, const size_type n )
  { BaseT::resize(BaseT::size()+n,t); }

  void pop_front( )
  {    for ( int64_t i = 1; i < this->jsize( ); i++ )
            (*this)[i-1] = (*this)[i];
       this->pop_back( );    }

  template <class U>
  void append( const vec<U>& y ) 
  {    this->insert( this->end( ), y.begin( ), y.end( ) );    }

  template <class U>
  void append( const SerfVec<U>& y ) 
  {    this->insert( this->end( ), y.begin( ), y.end( ) );    }
  
  void append( const vec<T>& y, size_type i, size_type j ) {
    if ( j == y.size( ) ) this->insert( this->end( ), y.begin( ) + i, y.end( ) );
    else this->insert( this->end( ), y.begin( ) + i, y.begin( ) + j );   
  }

  // appends values in y, but only those whose indices are in entries
  // IDX should be either (unsigned) int or longlong depending on the size of y
  template<typename IDX> 
  void append( const vec<T>& y, const vec<IDX>& entries ) {
    AssertLe(y.size(), numeric_limits<IDX>::max());
    const size_type n = this->size( );
    resize( n + entries.size( ) );
    for ( size_type j = 0; j < entries.size( ); j++ )
      (*this)[ n + j ] = y[ entries[j] ];
  }

  void appendFromBinaryFile( String const& filename )
  { BinaryIteratingReader< vec<T> > rdr(filename.c_str());
    reserve(this->size()+rdr.remaining());
    T tmp;
    while ( rdr.next(&tmp) ) push_back(tmp); }

// ===========================================================================
//
// FUNCTIONS TO TEST IF VECTOR HAS AN ATTRIBUTE
//
// ===========================================================================

  bool nonempty( ) const { return ! this->empty( ); }
  bool solo( ) const { return this->size( ) == 1; }

  bool Ordered( ) const
  {    for ( size_type i = 1; i < this->size( ); i++ )
            if ( (*this)[i] < (*this)[i-1] ) return False;
       return True;    }

  bool UniqueOrdered( ) const
  {    for ( size_type i = 1; i < this->size( ); i++ )
            if ( (*this)[i] <= (*this)[i-1] ) return False;
       return True;    }

// ===========================================================================
//
// FUNCTIONS TO EXTRACT PART OF A VECTOR
//
// ===========================================================================

  /// SetToSubOf: Set *this to the len entries of that, starting at pos.
  /// It is OK for that to be the same as this.
  void SetToSubOf( vec const& that, size_type pos, size_type len )
  { AssertLe(pos,that.size());
    AssertLe(len,that.size()-pos);
    if ( this != &that )
    { this->assign(that.begin()+pos,that.begin()+pos+len); }
    else
    { this->resize(pos+len);
      this->erase(this->begin(),this->begin()+pos); } }
  // Implementation note:  the STL carefully accounts for self-assignment, so
  // a simpler implementation would be simply
  // this->assign(that.begin()+pos,that.begin()+pos+len);
  // regardless of whether we are slicing ourselves or someone else.
  // However, the STL implementation appears to assume that the elements of the
  // vector also handle self-assignment correctly, which in the case of our
  // stuff hardly seems like a good assumption.  The resize/erase technique
  // implemented above is almost as efficient, and only assumes that element
  // copying to non-self works correctly.

  inline friend vec SubOf( const vec& x, size_type start, size_type n )  {
    vec<T> s(n);
    copy(x.begin() + start, x.begin() + start + n, s.begin());
    return s;
  }


  // SetToSubOf: Set *this to the entries in x pointed to by indices.
  // IDX should be either (unsigned) int or longlong depending on the size of x
  template<typename IDX>
  void SetToSubOf( const vec<T>& x, const vec<IDX>& indices ) {
    AssertLe(x.size(), numeric_limits<IDX>::max());
    const size_type n = indices.size( );
    resize(n);
    for ( size_type i = 0; i < n; i++ )
      (*this)[i] = x[ indices[i] ];
  }

  // IDX should be either (unsigned) int or longlong depending on the size of x
  template<typename IDX>
  inline friend vec SubOf( const vec& x, const vec<IDX>& indices ) {
    AssertLe(x.size(), numeric_limits<IDX>::max());
    const size_type n = indices.size( );
    vec<T> s(n);
    for ( size_type i = 0; i < n; i++ )
      s[i] = x[ indices[i] ];
    return s;
  }
    
  void SetToRangeOf( const vec<T>& v, size_type i, size_type j ) {
    AssertLe( i, j );
    this->assign(v.begin() + i, v.begin() + j);
  }

  inline friend vec RangeOf( const vec& v, size_type i, size_type j ) {
    AssertLe( i, j );
    vec<T> r(v.begin() + i, v.begin() + j);
    return r;
  }

  int isize( ) const { return this->size( ); }
  int64_t jsize( ) const { return this->size( ); }

  void SetCat( const vec<T>& v1, const vec<T>& v2 ) {
    *this = v1;
    append(v2);
  }

  void clear_and_resize( size_type n ) {
    this->clear( );
    resize(n);
  }

  void resize_and_set( size_type n, const T& x ) {
    this->clear();
    resize(n, x);
  }

  void SetToReverseOf( const vec<T>& v ) {
    resize(v.size());
    reverse_copy(v.begin(), v.end(), this->begin());
  }

  void ReverseMe( ) {
    reverse(this->begin(), this->end());
  }

  bool Contains( const vec<T>& v ) const {
    return (search(this->begin(), this->end(), v.begin(), v.end()) != this->end());
  }

  bool Contains( const vec<T>& v, size_type pos ) const {
    if ( v.size( ) + pos > this->size( ) )
      return False;
    size_type j;
    for ( j = 0; j < v.size( ); j++ )
      if ( (*this)[pos+j] != v[j] )
	break;
    return j == v.size( );
  }


  /// CountValue: count all entries having the given value.
  size_type CountValue( const T& x ) const {
    return count(this->begin(), this->end(), x);
  }

// ===========================================================================
//
// FUNCTIONS TO ERASE ELEMENTS FROM VECTORS - PART 1 (MEMBER FUNCTIONS)
//
// ===========================================================================

  /// Erase: erase range of elements, where range is given by half-open interval.

  void Erase( size_type start, size_type stop ) {
    this->erase( this->begin( ) + start, this->begin( ) + stop );
  }

  /// EraseValue: erase all entries having the given value.
  void EraseValue( const T& x ) {
    this->erase(remove(this->begin(), this->end(), x), this->end());
  }

  /// print values to ostream, separated by sep.
  void Print(ostream & os, const char * sep = " ") const {
    copy(this->begin(), this->end(), ostream_iterator<T>(os, sep));
  }

  /// print values to ostream, separated by sep, with newline at end.
  void Println(ostream & os, const char * sep = " ") const {
    Print(os, sep); os << endl;
  }

  ///Set myself from text stream containing list of values of unknown length.
  void ReadFromTextStream(istream & is) {
    this->clear();
    T t;
    while (true) {
      is >> t;
      if (!is) break;
      push_back(t);
    }
  }

  /// ReadSubset: read selected entries from file written with BinaryWrite.
  /// The type T must have a fixed external size.
  void ReadSubset( String const& filename, vec<int> const& ids,
                          bool append = false )
  { BinaryReader br(filename.c_str());
    size_t nnn;
    br.read(&nnn);
    size_t sz = br.externalSizeof( static_cast<T*>(0) );
    if ( !sz )
    { FatalErr("Can't randomly access the binary file " << filename
                << " which contains variable-sized elements."); }
    if ( !append ) this->clear();
    reserve(this->size()+ids.size());
    T tmp;
    typedef vec<int>::const_iterator Itr;
    for ( Itr itr(ids.begin()), end(ids.end()); itr != end; ++itr )
    { br.seekAndFill( *itr*sz + br.tell(), sz );
      br.read(&tmp);
      push_back(tmp); }  }

  /// ReadRange: read selected entries from file written with BinaryWrite.
  /// The type T must have a fixed external size.
  void ReadRange( String const& filename, size_t from, size_t to )
  { BinaryReader br(filename.c_str());
    size_t nnn;
    br.read(&nnn);
    ForceAssertLe(to,nnn);
    ForceAssertLe(from,to);
    size_t sz = br.externalSizeof( static_cast<T*>(0) );
    if ( !sz )
    { FatalErr("Can't randomly access the binary file " << filename
                << " which contains variable-sized elements."); }
    br.seek( from*sz + br.tell() );
    this->clear();
    reserve(to-from);
    T tmp;
    for ( size_t idx = from; idx != to; ++idx )
    { br.read(&tmp);
      push_back(tmp); }  }


// ===========================================================================
//
// MEMBER FUNCTIONS TO DO ARITHMETIC ON VECTORS.
//
// ===========================================================================

  /// Multiply a vector by a constant of type X.

  template<class X>
  vec & operator*=(const X & x) {
    const size_type S = this->size();
    for (size_type i = 0; i != S; ++i) {
      (*this)[i] = static_cast<T>( x * (*this)[i]);
    }
    return *this;
  }

  /// Divide a vector by a constant of type X.

  template<class X>
  vec & operator/=(const X & x) {
    return this->operator*=(1.0/x);
  }

  /// Add two vectors together.
  /// If vx is longer, add up to this vector's size.
  /// if vx is shorter, add up to vx's size only.

  template<class X> vec & operator+=(const vec<X> & vx) {
    const size_type S = min(this->size(), vx.size());
    for (size_type i = 0; i != S; ++i) {
      (*this)[i] += static_cast<T>(vx[i]);
    }
    return *this;
  }

  //stand-alone operators are implemented in terms of op=
  //See meyers, more effective C++, item 22 for reasons.

  /// Multiply a vector by a constant of type X.

  template<class X> friend vec operator*(const vec & v, const X & x) {
    return vec(v) *= x;
  }

  /// Multiply a vector by a constant of type X.

  template<class X> friend vec operator*(const X & x, const vec & v) {
    return vec(v) *= x;
  }
  /// Divide a vector by a constant of type X.

  template<class X> friend vec operator/(const vec & v, const X & x) {
    return vec(v) /= x;
  }

  /// Add two vectors together.
  /// If vx is longer, add up to this vector's size.
  /// if vx is shorter, add up to vx's size only.

  template<class X> friend vec operator+(const vec & v, const vec<X> & vx) {
    return vec(v) += vx;
  }

  /// NextDiff(i): return index of next element after i that is different from
  /// the ith element.

  inline size_type NextDiff( size_type i ) const
  {    size_type j;
       for ( j = i + 1; j < this->size( ); j++ )
            if ( (*this)[j] != (*this)[i] ) break;
       return j;    }

  friend int compare( vec const& v1, vec const& v2 )
  {   auto itr1 = v1.begin();
      auto itr2 = v2.begin();
      using std::min;
      auto end = itr1 + min(v1.size(),v2.size());
      int result = 0;
      while ( !result && itr1 != end )
      { result = compare(*itr1,*itr2); ++itr1; ++itr2; }
      if ( !result ) result = compare(v1.size(),v2.size());
      return result;
  }

  friend void swap( vec& v1, vec& v2 )
  { v1.swap(v2); }

 private:
  static size_t ValidatedSize( size_t nnn )
  {
#ifndef NDEBUG
      if ( nnn > 1000ul*1000ul*1000ul*1000ul/sizeof(T) )
          FatalErr( "Vector too big: Attempt to resize vec<T> object to " << nnn
                      << ", where sizeof(T) = " << sizeof(T) << ".");
#endif
      return nnn;
  }
};

template<class T> bool IsUnique( vec<T> x )
{    Sort(x);
     return x.UniqueOrdered( );    }

template <class T>
struct Serializability< vec<T> >
{ typedef ExternallySerializable type; };

//
//  End of vec Class Declaration and Template Definitions
//
/////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////
//
// Mutating Functions
//

template<class V>
V Reverse( V const& v )
{ return V( v.rbegin( ), v.rend( ) ); }

template<class T,class A>
void ReverseThis( vec<T,A>& v )
{ v.ReverseMe(); }

template<class T,class A>
void RandomShuffle( vec<T,A>& v )
{ random_shuffle( v.begin( ), v.end( ) ); }

//
// OverlapAppend - Given a sequence like A B C D C D and C D C D E F
// find the largest possible overlap and append the non-overlapping
// portion of the first sequence to the second:
//
//    v1:  A B C D C D
//    v2:      C D C D E F
//
// -> v1: A B C D C D E F
//
// returns the number of overlapping entries.  Note that the sequences
// must be vec-like things of the same thing (e.g. vec<int> and SerfVec<int>)
//
//

template <
     template <typename , typename... > class V1,     // e.g. vec<T,A1>
     template <typename , typename... > class V2,     // e.g. SerfVec<T,A2,A3>
     typename T,
     typename... Args1,
     typename... Args2>
size_t OverlapAppend( V1<T,Args1...>& v1, V2<T,Args2...> const& v2 ) {
     size_t best = 0;
     size_t v1_size = v1.size(), v2_size = v2.size();
     size_t max_overl = std::min(v1_size, v2_size);     // max possible overlap
     // for each possible overlap, from longest to shortest...
     for ( size_t  overl = max_overl;  overl; --overl ) {
          if ( std::equal(v1.end()-overl, v1.end(), v2.begin() ) ) {
               best = overl;
               break;
          }
     }

     std::copy( v2.begin()+best, v2.end(), std::back_inserter(v1));

     return best;
}


/////////////////////////////////////////////////////////////////////////////
//
// Search Functions
//

template<class T,class A>
typename vec<T,A>::size_type LowerBound( const vec<T,A>& v, const T& x )
{ return std::lower_bound( v.begin( ), v.end( ), x ) - v.begin( ); }

template<class T,class A>
typename vec<T,A>::size_type UpperBound( const vec<T,A>& v, const T& x )
{ return std::upper_bound( v.begin( ), v.end( ), x ) - v.begin( ); }

template<class T,class A>
inline bool Member( vec<T,A> const& v, const T& x )
{ return (std::find(v.begin(), v.end(), x) != v.end()); }

/// Return the position of an element in a vector, else -1.
template<class T,class A>
inline long Position( vec<T,A> const& v, const T& x )
{
  using std::distance;
  auto pos = std::find(v.begin(), v.end(), x);
  if (pos != v.end())
    return distance(v.begin(),pos);
  else
    return -1L;
}

// Position( v, w ): return first position of w in v, else -1.

template<class T,class A>
inline int64_t Position( const vec<T,A>& v, const vec<T,A>& w )
{    int64_t pos = std::search(v.begin(),v.end(),w.begin(),w.end()) - v.begin();
     int64_t end = v.jsize( );
     if ( pos == end ) return -1;
     else return pos;    }

/// BinPosition.  Return the position of an element in a sorted vector, else -1.
/// If the element appears more than once, the position of one of its instances
/// is returned.
template<class V,class U>
inline long BinPosition( V const& v, const U& x1 )
{    if ( v.size( ) == 0 ) return -1;
     typename V::const_reference& x(x1);
     typename V::size_type first = 0, last = v.size( ) - 1, next;
     while (1)
     {    if (first == last) return ( !(x < v[last]) && !(v[last] < x) ) ? last : -1;
          next = first + (last - first) / 2;
          if ( x < v[next] ) last = next;
          else if ( v[next] < x ) first = next + 1;
          else return next;    }    }

template<class V,class U>
inline bool BinMember( V const& v, const U& x )
{ return std::binary_search(v.begin(),v.end(),x); }

/// BinSubset: determine if v is a subset of w; assumes w only is sorted and that there
/// is no repetition.
template<class V>
inline bool BinSubset( V const& v, V const& w )
{    for ( typename V::const_reference e : v )
          if ( !BinMember(w,e) ) return False;
     return True;    }

/// Determine if v is a subset of w; assumes that there is no repetition.
template<class T,class A,class A1>
inline bool Subset( const vec<T,A>& v, const vec<T,A1>& w )
{    for ( T const& e : v )
          if ( !Member(w,e) ) return False;
     return True;    }



/////////////////////////////////////////////////////////////////////////////
//
// Numerical Functions
//

template<class T,class A>
size_t SizeSum( vec<T,A> const& v, typename T::value_type* = 0 )
{ size_t sum = 0;
  for ( T const& e : v ) sum += e.size();
  return sum; }

template<class A>
inline bool Nonnegative( const vec<int,A>& v )
{ for ( int e : v ) if ( e < 0 ) return false;
  return true;    }


/////////////////////////////////////////////////////////////////////////////
//
// Destroy - returns the memory used by a vector, more or less
//

// Method given in Stroustrup's book, The C++ Programming Language, Special 
// Edition (2000), p. 457.

template<class T,class A>
inline void Destroy( vec<T,A>& v )
{ vec<T,A> empty; swap(v,empty); }


/////////////////////////////////////////////////////////////////////////////
//
// Sort Functions
//

template<class T,class A>
inline void Sort( vec<T,A>& v )
{
  TRACEVAL_STOP_TRACING_COPIES;
  using std::sort; sort( v.begin( ), v.end( ) );
  TRACEVAL_START_TRACING_COPIES;
}

template<class T,class A,class StrictWeakOrdering>
inline void Sort( vec<T,A>& v, StrictWeakOrdering comp )
{
  TRACEVAL_STOP_TRACING_COPIES;
  using std::sort; sort( v.begin( ), v.end( ), comp );
  TRACEVAL_START_TRACING_COPIES;
}

template<class T,class A>
inline void ReverseSort( vec<T,A>& v )
{
  TRACEVAL_STOP_TRACING_COPIES;
  using std::sort; sort( v.rbegin( ), v.rend( ) );
  TRACEVAL_START_TRACING_COPIES;
}

template<class T,class A,class StrictWeakOrdering>
inline void ReverseSort( vec<T,A>& v, StrictWeakOrdering comp )
{
  TRACEVAL_STOP_TRACING_COPIES;
  using std::sort; sort( v.rbegin( ), v.rend( ), comp );
  TRACEVAL_START_TRACING_COPIES;
}

template <class T,class A,class Comp=std::equal_to<T>>
void Unique( vec<T,A>& v, Comp isEqual=Comp() )
{
    if ( v.size() <= 1 )
        return;

    TRACEVAL_STOP_TRACING_COPIES;

    using std::iter_swap;
    auto dest(v.begin());
    auto end(v.end());
    for ( auto itr(dest+1); itr != end; ++itr )
        if ( !isEqual(*itr,*dest) )
            iter_swap(itr, ++dest);
    v.erase(dest+1,end);

    TRACEVAL_START_TRACING_COPIES;
}

/// Leaves only unique elements in a vector \c v; these elements
/// will be also sorted.
template <class T,class A,class StrictWeakOrdering,class EqualPredicate>
inline void UniqueSort(vec<T,A>& v,StrictWeakOrdering comp,EqualPredicate equal)
{
    if ( v.size() <= 1 )
        return;
    Sort(v, comp);
    Unique(v, equal);
}

/// Leaves only unique elements in a vector \c v; these elements
/// will be also sorted.
template <class T,class A>
inline void UniqueSort(vec<T,A> & v)
{
  UniqueSort(v,std::less<T>(),std::equal_to<T>());
}



/////////////////////////////////////////////////////////////////////////////
//
// Erasing Functions
//

/// EraseIf: wrapper around erase-remove_if idiom.
template<class T,class A>
void EraseIf( vec<T,A>& v, bool (T::*f)( ) const )
{
  TRACEVAL_STOP_TRACING_COPIES;
  v.erase( remove_if( v.begin( ), v.end( ), std::mem_fun_ref(f) ), v.end( ) );
  TRACEVAL_START_TRACING_COPIES;
}

/// Another version of EraseIf: erase v[x] if erase[x] = True.
template<class T,class A,class A1>
void EraseIf( vec<T,A>& v, const vec<Bool,A1>& erase )
{
  TRACEVAL_STOP_TRACING_COPIES;
  typename vec<T,A>::size_type count = 0;
  for ( typename vec<T,A>::size_type i = 0; i < v.size( ); i++ ) {
    if ( ! erase[i] ) {
      if ( count != i ) 
        v[count] = v[i];
      ++count;
    }
  }
  v.resize(count);
  TRACEVAL_START_TRACING_COPIES;
}

/// Another version of EraseIf: erase v[x] if erase[x] = True.
template<class T,class A,class A1>
void EraseUnless( vec<T,A>& v, const vec<Bool,A1>& keep )
{
  TRACEVAL_STOP_TRACING_COPIES;
  typename vec<T,A>::size_type count = 0;
  for ( typename vec<T,A>::size_type i = 0; i < v.size( ); i++ ) {
    if ( keep[i] ) {
      if ( count != i ) 
        v[count] = v[i];
      ++count;
    }
  }
  v.resize(count);
  TRACEVAL_START_TRACING_COPIES;
}

/// EraseTheseIndices: Erase some elements of a vector, as determined by a sorted
/// list of indices.  Not efficiently implemented.
template<class T,class A,class IDX,class A1>
void EraseTheseIndices( vec<T,A>& v, const vec<IDX,A1>& these )
{
  AssertLe(v.size(), static_cast<size_t>(numeric_limits<IDX>::max()));
  TRACEVAL_STOP_TRACING_COPIES;
  typename vec<T,A>::size_type count = 0;
  for ( typename vec<T,A>::size_type i = 0; i < v.size( ); i++ )
    {    if ( !BinMember( these, i ) )
      {    if ( count != i ) v[count] = v[i];
      ++count;    }    }
  v.resize(count);
  TRACEVAL_START_TRACING_COPIES;
}


/////////////////////////////////////////////////////////////////////////////
//
// Input and Output Functions
//

/// IsAsciiVec() returns whether or not filename is a saved vec.  It
/// does this by checking whether the first line of the file is an
/// ASCII representation of a number.  If it encounters a non-digit,
/// non-whitespace character before it finds a newline, it returns
/// false.  Otherwise, it returns true.
bool IsAsciiVec( const String &filename );

/// Get the number of elements in saved ASCII vec.
longlong AsciiVecSize( const String& filename );

inline void AsciiBoolVecReadSubset( const String& filename,
                                         const vec<int>& ids,
                                         vec<Bool>& v )
{    String ns;
     {    Ifstream( in, filename );
          in >> ns;    }
     ForceAssert( ns.IsInt() );
     longlong n = ns.Int();
     int k = ns.size()+1;
     FileReader fr(filename.c_str());
     v.resize( ids.size( ) );
     for ( int i = 0; i < ids.isize( ); i++ )
     {    ForceAssertGe( ids[i], 0 );
          ForceAssertLt( ids[i], n );
          fr.seek( k + ids[i] );
          fr.read( &v[i], 1 );    }    }

/// Return number of elements in binary file of a vec<T> of some sort.
inline size_t BinaryVecNumElements( const String & filename )
{   BinaryReader br(filename.c_str());
    size_t nnn;
    br.read(&nnn);
    return nnn; }

inline size_t BinaryVecElementSize( String const& filename )
{    BinaryReader br(filename.c_str());
     size_t nnn;
     br.read(&nnn);
     size_t dataLen = br.getFilesize() - br.tell();
     ForceAssertEq(dataLen%nnn,0ul);
     return dataLen / nnn; }

#define BREAD2( FILE, TYPE, DATA ) \
    TYPE DATA; BinaryReader::readFile( FILE, &DATA )



// ============================================================================
// ========================================================================

void PrettyPrint( ostream& o, const vec<int>& v, int max_items = 0,
     String terminator = "\n" );

void PrettyPrint( ostream& o, const vec<longlong>& v, int max_items = 0,
     String terminator = "\n" );

void PrettyPrint( ostream& o, const vec<double>& v, int max_items = 0,
     String terminator = "\n" );

void PrettyPrint( ostream& o, const vec<TraceInt>& v, int max_items = 0,
     String terminator = "\n" );


// CompactPrint: print a vector of objects, separated by blanks (default) or
// a user-specified separator.

template<class T> void CompactPrint( ostream& out, const vec<T>& v,
    String separator = " " )
{    for ( typename vec<T>::size_type i = 0; i < v.size( ); i++ )
     {    if ( i > 0 ) out << separator;
          out << v[i];    }    }


/// WriteAppend has an implementation for T = alignment_plus in Alignment.{h,cc}.

template<class T> void WriteAppend( const String& f, const vec<T>& v )
{    ForceAssert( !IsRegularFile( f + ".gz" ) );
     if ( !IsRegularFile(f) )
     {    std::ofstream out(f.c_str());
          out << setfill('0') << setw(15) << v.size() << setfill(' ') << '\n';
          for ( typename vec<T>::size_type i = 0; i <  v.size( ); i++ )
               out << v[i];
          out.close();   }
     else
     {    std::fstream out(f.c_str(),std::ios_base::in|std::ios_base::out);
          size_t n;
          out >> n;
          size_t const max_size_bound = 10000000ul * 100000000ul;
          ForceAssertLt( n+v.size(), max_size_bound );
          out.seekp( 0, std::ios_base::beg );
          out << setfill('0') << setw(15) << n+v.size( ) << setfill(' ') << '\n';
          out.seekp( 0, std::ios_base::end );
          for ( typename vec<T>::size_type i = 0; i <  v.size( ); i++ )
               out << v[i];
          out.close(); }    }

/// a specialized version of WriteAppend for String

template <>
void WriteAppend( const String& f, const vec<String>& v );

template<class T> ostream& operator<<(ostream& s, const vec<T>& v)
{    s << v.size( ) << "\n";
     for ( typename vec<T>::size_type i = 0; i < v.size( ); i++ )
          s << v[i];
     return s;    }

ostream& operator<<(ostream& s, const vec<unsigned short>& v);
ostream& operator<<(ostream& s, const vec<int>& v);
ostream& operator<<(ostream& s, const vec<longlong>& v);
ostream& operator<<(ostream& s, const vec<float>& v);
ostream& operator<<(ostream& s, const vec<double>& v);
ostream& operator<<(ostream& s, const vec<String>& v);

template<class T> istream& operator>>(istream& s, vec<T>& v)
{    typename vec<T>::size_type n;
     s >> n;
     v.resize(n);
     char c;
     s.get(c);
     for ( typename vec<T>::size_type i = 0; i < v.size( ); i++ )
       s >> v[i];   // Breaks cxx
     return s;    }


istream& operator>>(istream& s, vec<String>& v);

/// Print out a matrix, with left-justified entries, and given separation between
/// columns.  (Justification may be changed by supplying an optional argument
/// consisting of a string of l's and r's.)

void PrintTabular( ostream& out, const vec< vec<String> >& rows, int sep,
     String justify = String( ) );

void PrintCSV(ostream& out, const vec< vec<String> >& rows);


#define For_(T,x,v) for( vec< T >::const_iterator x = v.begin(); x != v.end(); ++x )
#define ForMut_(T,x,v) for( vec< T >::iterator x = v.begin(); x != v.end(); ++x )

#endif
