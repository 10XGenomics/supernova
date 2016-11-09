///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/* VecUtilities.h
 *
 * Collection of vec class utilities and helpers. 
 *
 * A grab bag of functions that deal with vecs, moved from Vec.h
 * PLEASE help keep Vec.h tidy by placing any new functions here.
 *
 * See also: Vec.h
 *
 * Contains:
 *
 * MkVec
 * JoinVec
 * UniqueSortAndCount
 * SortSync
 * ReverseSortSync
 * UniqueSortSync
 * SortIndex
 * Intersection
 * Meet
 * WhatPermutation
 * PermuteVec
 * indirect_compare
 */

#ifndef VEC_UTILITIES_H
#define VEC_UTILITIES_H

#include "Vec.h"
#include <cstddef>


/////////////////////////////////////////////////////////////////////////////
//
//  MkVec - creates a new vector from a list of entries
//
/////////////////////////////////////////////////////////////////////////////

template<class T> vec<T> MkVec( const T& v1 ) {
  vec<T> v;
  v.push_back( v1 );
  return v;
}

template<class T> vec<T> MkVec( const T& v1, const T& v2 ) {
  vec<T> v;
  v.push_back( v1, v2 );
  return v;
}

template<class T> vec<T> MkVec( const T& v1, const T& v2, const T& v3 ) {
  vec<T> v;
  v.push_back( v1, v2, v3 );
  return v;
}

template<class T> vec<T> MkVec( const T& v1, const T& v2, const T& v3,
				const T& v4 ) {
  vec<T> v;
  v.push_back( v1, v2, v3, v4 );
  return v;
}

template<class T> vec<T> MkVec( const T& v1, const T& v2, const T& v3,
				const T& v4, const T& v5 ) {
  vec<T> v;
  v.push_back( v1, v2, v3, v4, v5 );
  return v;
}

template<class T> vec<T> MkVec( const T& v1, const T& v2, const T& v3,
				const T& v4, const T& v5, const T& v6 ) {
  vec<T> v;
  v.push_back( v1, v2, v3, v4, v5, v6 );
  return v;
}

template<class T> vec<T> MkVec( const T& v1, const T& v2, const T& v3,
				const T& v4, const T& v5, const T& v6,
				const T& v7 ) {
  vec<T> v;
  v.push_back( v1, v2, v3, v4, v5, v6, v7 );
  return v;
}

template<class T> vec<T> MkVec( const T& v1, const T& v2, const T& v3,
				const T& v4, const T& v5, const T& v6,
				const T& v7, const T& v8 ) {
  vec<T> v;
  v.push_back( v1, v2, v3, v4, v5, v6, v7, v8 );
  return v;
}

template<class T> vec<T> MkVec( const T& v1, const T& v2, const T& v3,
				const T& v4, const T& v5, const T& v6,
				const T& v7, const T& v8, const T& v9 ) {
  vec<T> v;
  v.push_back( v1, v2, v3, v4, v5, v6, v7, v8, v9 );
  return v;
}


/////////////////////////////////////////////////////////////////////////////
//
//  JoinVec - creates a new vector from a list of vectors
//
/////////////////////////////////////////////////////////////////////////////

template<class T> vec<T> JoinVecs( const vec<T>& v1, const vec<T>& v2 ) {
  vec<T> v( v1 );
  v.append( v2 );
  return v;
}

template<class T> vec<T> JoinVecs( const vec<T>& v1, const vec<T>& v2, const vec<T>& v3 ) {
  vec<T> v( v1 );
  v.append( v2 );
  v.append( v3 );
  return v;
}

template<class T> vec<T> JoinVecs( const vec<T>& v1, const vec<T>& v2, const vec<T>& v3,
				   const vec<T>& v4 ) {
  vec<T> v( v1 );
  v.append( v2 );
  v.append( v3 );
  v.append( v4 );
  return v;
}

template<class T> vec<T> JoinVecs( const vec<T>& v1, const vec<T>& v2, const vec<T>& v3,
				   const vec<T>& v4, const vec<T>& v5 ) {
  vec<T> v( v1 );
  v.append( v2 );
  v.append( v3 );
  v.append( v4 );
  v.append( v5 );
  return v;
}

template<class T> vec<T> JoinVecs( const vec<T>& v1, const vec<T>& v2, const vec<T>& v3,
				   const vec<T>& v4, const vec<T>& v5, const vec<T>& v6 ) {
  vec<T> v( v1 );
  v.append( v2 );
  v.append( v3 );
  v.append( v4 );
  v.append( v5 );
  v.append( v6 );
  return v;
}


/////////////////////////////////////////////////////////////////////////////
//
//  UniqueSortAndCount
//
/////////////////////////////////////////////////////////////////////////////

/// After this method is applied to a vector \c v (first argument), this vector
/// will contain only unique elements (as defined by EqualPredicate)
/// found in its original content. These elements will be sorted in ascending order
/// (according to the StrictWeakOrdering). The number of times each 
/// returned unique element v[i] occured in the original content of the vector v will
/// be returned in counts[i] (old content of counts vector, if any, will be destroyed). 

template <class T, class StrictWeakOrdering, class EqualPredicate> 
void UniqueSortAndCount(vec<T> & v, StrictWeakOrdering comp, EqualPredicate equal, 
			vec<unsigned int> & counts) {
  counts.clear();
  if ( v.empty() ) return;

  TRACEVAL_STOP_TRACING_COPIES;
  
  Sort( v, comp );
  
  typename vec<T>::iterator iter = v.begin();
  typename vec<T>::iterator iter_end = v.end();
  typename vec<T>::iterator last_unique_elem_iter = iter; 
  typename vec<T>::iterator move_to = v.begin();

  // we know that v is non-empty, so it's ok to increase iter
  // right away in the next line:
  for ( ++iter; iter != iter_end ; ++iter ) {
    // we follow the standard assumptions of STL here.
    // the latters assume that operator== is defined
    // (but not necessarily !=):
    if ( equal(*iter,*last_unique_elem_iter ) ) {
      continue;
    }
    counts.push_back( distance(last_unique_elem_iter, iter) );
    *(++move_to) = *(last_unique_elem_iter = iter);
  }
  counts.push_back( distance(last_unique_elem_iter, iter) );
  v.erase(++move_to,iter_end);

  TRACEVAL_START_TRACING_COPIES;
}


/// When this method is applied to a vector \c v (first argument), it will be 
/// first shrunk to contain only unique elements (as defined by EqualPredicate)
/// found in its original content. These elements will be also sorted in 
/// ascending order (according to the StrictWeakOrdering) and the number of 
/// times each unique element occured in the original content of the vector v 
/// will be counted. The FilterPredicate will be then applied to each such 
/// count and upon return from this method <v> and <counts> will contain only 
/// unique elements (still sorted) and counts, respectively, such that for each 
/// counts[i] the value of passes(counts[i]) is <true>

template <class T, class StrictWeakOrdering, class EqualPredicate, class FilterPredicate> 
void UniqueSortAndCount(vec<T> & v, StrictWeakOrdering comp, EqualPredicate equal, 
			      vec<unsigned int> & counts, FilterPredicate pass) {
  counts.clear();
  if ( v.empty() ) return;

  TRACEVAL_STOP_TRACING_COPIES;
  
  Sort( v, comp );
  
  typename vec<T>::iterator iter = v.begin();
  typename vec<T>::iterator iter_end = v.end();
  typename vec<T>::iterator last_unique_elem_iter = iter; 
  typename vec<T>::iterator move_to = v.begin();
  --move_to; // point before the first element - in this
             // method we may not need to keep the first element at all!

  // we know that v is non-empty, so it's ok to increase iter
  // right away in the next line:
  for ( ++iter; iter != iter_end ; ++iter ) {
    // we follow the standard assumptions of STL here.
    // the latters assume that operator== is defined
    // (but not necessarily !=):
    if ( equal(*iter,*last_unique_elem_iter ) ) {
      continue; // we are inside a stretch of same element repetions, 
                // keep counting...
    }
    unsigned int cnt = distance(last_unique_elem_iter, iter); 
    if ( pass(cnt) ) {
      counts.push_back(cnt); // record count and store element only if passes()
      *(++move_to) = *(last_unique_elem_iter) ;
    }
    last_unique_elem_iter = iter; // we just encountered new unique element,
                                  // record its position (but not store yet!)
  }
  unsigned int cnt = distance(last_unique_elem_iter, iter_end); 
  if ( pass(cnt) ) { // take care of the last unique elem - not stored yet!
      counts.push_back( cnt );
      *(++move_to) = *(last_unique_elem_iter) ;
  }      
  v.erase(++move_to,iter_end);

  TRACEVAL_START_TRACING_COPIES;
}

/// After this method is applied to a vector \c v (first argument), this vector
/// will contain only unique elements
/// found in its original content. These elements will be sorted in ascending order.
/// The number of times each returned unique element v[i] occured
/// in the original content of the vector v will
/// be returned in counts[i] (old content of counts vector, if any, will be destroyed). 

template <class T> 
inline  void UniqueSortAndCount(vec<T> & v, vec<unsigned int> & counts) {
  UniqueSortAndCount(v,less<T>(),equal_to<T>(),counts);
}

/// When this method is applied to a vector \c v (first argument), it will be 
/// first shrunk to contain only unique elements.
/// These elements will be also sorted in 
/// ascending order and the number of 
/// times each unique element occured in the original content of the vector v 
/// will be counted. The FilterPredicate will be then applied to each such 
/// count and upon return from this method <v> and <counts> will contain only 
/// unique elements (still sorted) and counts, respectively, such that for each 
/// counts[i] the value of passes(counts[i]) is <true>

template <class T, class FilterPredicate> 
inline  void UniqueSortAndCount(vec<T> & v, vec<unsigned int> & counts, FilterPredicate pass) {
  UniqueSortAndCount(v, less<T>(), equal_to<T>(), counts, pass);
}


/////////////////////////////////////////////////////////////////////////////
//
//  SortSync - sort many vectors based on the order of the first
//
/////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
//
//  Helpers for SortSync and SortIndex
//
/////////////////////////////////////////////////////////////////////////////

/**Permute input vector v in place according to permutation.
   Preconditions:
   - v.size() == permutation.size()
   If the permutation contains a -1, the position corresponding to that
   is essentially ignored and ends up in one of the available empty spaces

   This works for all std::vectors and vecs.
*/
template<class V, typename IDX>
void PermuteVec( V& v, const vec<IDX>& permutation )
{
    AssertEq(v.size(), permutation.size());
    // make sure the index type is large enough
    AssertLe(v.size(),static_cast<size_t>(std::numeric_limits<IDX>::max()));
    IDX n = static_cast<IDX> (v.size());
    vec<IDX> o = permutation;
    const bool unsigned_index = (std::numeric_limits<IDX>::is_signed == false);
    // -1 for a signed index, don't care for unsigned
    const IDX minus1 = static_cast<IDX> (-1);
    for ( IDX i = 0; i != n; ++i )
    {
        while ( o[i] != i && (unsigned_index || o[i] != minus1) )
        { // ignore -1 check for unsigned
            using std::swap;
            swap(v[i], v[o[i]]);
            swap(o[i], o[o[i]]);
        }
    }
}

// Instantiate with a vector v, then use in place of operator<.
// It will tell you i<j if v[i]<v[j].  Helper for WhatPermutation.
template<class V, typename C = std::less<typename V::value_type>>
struct indirect_compare : public std::binary_function<size_t, size_t, bool>
{
    indirect_compare( const V& v, C c = C() ) : mV(v), mComp(c) {}

    bool operator()( size_t i, size_t j ) const
    { return mComp(mV[i],mV[j]); }

    const V& mV;
    C mComp;
};

// What permutation would we apply to V to get it sorted?
// NOTE: If v is a vec<int>, this inverts it.
template<class V, typename C, typename IDX>
void WhatPermutation( const V & v, vec<IDX>& permutation, C comp,
                        bool inv = true )
{
    AssertLe(v.size(), static_cast<size_t>(std::numeric_limits<IDX>::max()));
    IDX n = v.size();
    vec<IDX> perm(n);
    for ( IDX i = 0; i < n; i++ )
        perm[i] = i;
    std::sort(perm.begin(),perm.end(),indirect_compare<V,C>(v,comp));
    if ( !inv )
    {
        using std::swap;
        swap(perm, permutation);
    }
    else
    {
        // That's the inverse of the permutation PermuteVec takes.
        permutation.resize(n);
        for ( IDX i = 0; i < n; i++ )
            permutation[perm[i]] = i;
    }
}

template<class V, typename IDX>
void WhatPermutation( const V & v, vec<IDX>& permutation, bool inv = true )
{
    WhatPermutation(v, permutation, less<typename V::value_type>(), inv);
}

// SortSync( vec& v, vec& w ): sort v, moving the elements of w synchronously.

template<class V1, class V2, typename C, typename IDX>
void SortSyncDispatch( V1& v, V2& w, C comp, IDX )
{
    ForceAssertEq( v.size(), w.size() );
    vec<IDX> perm;
    WhatPermutation(v, perm, comp);
    PermuteVec(v, perm);
    PermuteVec(w, perm);
}

template<class V1, class V2, typename C=std::less<typename V1::value_type>>
void SortSync( V1& v, V2& w, C comp= C(),
                decltype(comp(v[0],v[0])) = false, decltype(w.size()) = 0 )
{
    if ( v.size() < std::numeric_limits<unsigned>::max() )
        SortSyncDispatch(v,w,comp,0u);
    else
        SortSyncDispatch(v,w,comp,0ul);
}

template<class V1, class V2>
void ReverseSortSync( V1& v, V2& w )
{
    SortSync(v,w,std::greater<typename V1::value_type>());
}

template<class V1, class V2>
void UniqueSortSync( V1& v, V2& w )
{
    SortSync(v, w, less<typename V1::value_type>());
    typename V1::size_type count = 0;
    for ( typename V1::size_type i = 0; i < v.size(); i++ )
    {
        if ( i > 0 && v[i] == v[i - 1] )
            continue;
        if ( count != i )
        {
            v[count] = v[i];
            w[count] = w[i];
        }
        ++count;
    }
    v.resize(count), w.resize(count);
}

// SortSync( v, w, x )-type functions.

template<class V1, class V2, class V3, typename C, typename IDX>
void SortSyncDispatch( V1& v, V2& w, V3& x, C comp, IDX )
{
    ForceAssertEq( v.size(), w.size() );
    ForceAssertEq( v.size(), x.size() );
    vec<IDX> perm;
    WhatPermutation(v, perm, comp);
    PermuteVec(v, perm);
    PermuteVec(w, perm);
    PermuteVec(x, perm);
}

template<class V1, class V2, class V3, typename C=std::less<typename V1::value_type>>
void SortSync( V1& v, V2& w, V3& x, C comp = C(),
                decltype(comp(v[0],v[0])) = false, decltype(x.size()) = 0 )
{
    if ( v.size() < std::numeric_limits<unsigned>::max() )
        SortSyncDispatch(v,w,x,comp,0u);
    else
        SortSyncDispatch(v,w,x,comp,0ul);
}

template<class V1, class V2, class V3>
void ReverseSortSync( V1& v, V2& w, V3& x )
{
    SortSync(v,w,x,greater<typename V1::value_type>());
}

// SortSync( v, w, x, y )-type functions.

template<class V1, class V2, class V3, class V4, typename C, typename IDX>
void SortSyncDispatch( V1& v, V2& w, V3& x, V4& y, C comp, IDX )
{
    ForceAssertEq( v.size(), w.size() );
    ForceAssertEq( v.size(), x.size() );
    ForceAssertEq( v.size(), y.size() );
    vec<IDX> perm;
    WhatPermutation(v, perm, comp);
    PermuteVec(v, perm);
    PermuteVec(w, perm);
    PermuteVec(x, perm);
    PermuteVec(y, perm);
}

template<class V1, class V2, class V3, class V4, typename C=std::less<typename V1::value_type>>
void SortSync( V1& v, V2& w, V3& x, V4& y, C comp = C(),
                decltype(comp(v[0],v[0])) = false, decltype(y.size()) = 0 )
{
    if ( v.size() <= std::numeric_limits<unsigned>::max() )
        SortSyncDispatch(v,w,x,y,comp,0u);
    else
        SortSyncDispatch(v,w,x,y,comp,0ul);
}

template<class V1, class V2, class V3, class V4, class V5, typename C,
        typename IDX>
void SortSyncDispatch( V1& v, V2& w, V3& x, V4& y, V5& z, C comp, IDX )
{
    ForceAssertEq( v.size(), w.size() );
    ForceAssertEq( v.size(), x.size() );
    ForceAssertEq( v.size(), y.size() );
    ForceAssertEq( v.size(), z.size() );
    vec<IDX> perm;
    WhatPermutation(v, perm, comp);
    PermuteVec(v, perm);
    PermuteVec(w, perm);
    PermuteVec(x, perm);
    PermuteVec(y, perm);
    PermuteVec(z, perm);
}

template<class V1, class V2, class V3, class V4, class V5, typename C=std::less<typename V1::value_type>>
void SortSync( V1& v, V2& w, V3& x, V4& y, V5& z, C comp = C(),
                decltype(comp(v[0],v[0])) = false, decltype(z.size()) = 0 )
{
    if ( v.size() <= std::numeric_limits<unsigned>::max() )
        SortSyncDispatch(v,w,x,y,z,comp,0u);
    else
        SortSyncDispatch(v,w,x,y,z,comp,0ul);
}

template<class V1, class V2, class V3, class V4, class V5, class V6,
        typename C, typename IDX>
void SortSyncDispatch( V1& v, V2& w, V3& x, V4& y, V5& z, V6& u, C comp, IDX )
{
    ForceAssertEq( v.size(), w.size() );
    ForceAssertEq( v.size(), x.size() );
    ForceAssertEq( v.size(), y.size() );
    ForceAssertEq( v.size(), z.size() );
    ForceAssertEq( v.size(), u.size() );
    vec<IDX> perm;
    WhatPermutation(v, perm, comp);
    PermuteVec(v, perm);
    PermuteVec(w, perm);
    PermuteVec(x, perm);
    PermuteVec(y, perm);
    PermuteVec(z, perm);
    PermuteVec(u, perm);
}

template<class V1, class V2, class V3, class V4, class V5, class V6, typename C=std::less<typename V1::value_type>>
void SortSync( V1& v, V2& w, V3& x, V4& y, V5& z, V6& u, C comp = C(),
                decltype(comp(v[0],v[0])) = false, decltype(u.size()) = 0 )
{
    if ( v.size() <= std::numeric_limits<unsigned>::max() )
        SortSyncDispatch(v,w,x,y,z,u,comp,0u);
    else
        SortSyncDispatch(v,w,x,y,z,u,comp,0ul);
}

template<class V1, class V2, class V3, class V4>
void ReverseSortSync( V1& v, V2& w, V3& x, V4& y )
{
    SortSync(v,w,x,y,greater<typename V1::value_type>());
}

template<class V1, class V2, class V3, class V4, class V5>
void ReverseSortSync( V1& v, V2& w, V3& x, V4& y, V5& z )
{
    SortSync(v,w,x,y,z,greater<typename V1::value_type>());
}

template<class V1, class V2, class V3, class V4, class V5, class V6>
void ReverseSortSync( V1& v, V2& w, V3& x, V4& y, V5& z, V6& u ) {
  SortSync(v,w,x,y,z,u,greater<typename V1::value_type>());
}

/////////////////////////////////////////////////////////////////////////////
//
//  SortIndex
//
/////////////////////////////////////////////////////////////////////////////

// SortIndex( v, index ): calculates an index into v such that
// v[index[n]] returns the nth value of a sorted v.
template<class V1, typename IDX, typename C = std::less<typename V1::value_type>>
void SortIndex( const V1& v, vec<IDX>& index, C comp = C() )
{
    AssertLe(v.size(), static_cast<size_t>(std::numeric_limits<IDX>::max()));
    index.resize(v.size());
    iota(index.begin(), index.end(), 0);
    vec<IDX> perm;
    WhatPermutation(v, perm, comp);
    PermuteVec(index, perm);
}


/////////////////////////////////////////////////////////////////////////////
//
//  Intersection & Meet - functions to find shared elements
//
/////////////////////////////////////////////////////////////////////////////

/// Intersection: make pairs of all elements in the first vector that also
/// occur in the second (which can result in duplicates).  Put the pairs
/// in the result. Both vectors must be SORTED.  This version is meant
/// for situations where two objects may compare equal but contain different
/// information, and we want to preserve both copies of that information.

template<class T> void
Intersection( const vec<T>& v1, const vec<T>& v2, vec<pair<T,T> > & result )
{    
  result.clear();
  typename vec<T>::const_iterator v1iter = v1.begin();
  typename vec<T>::const_iterator v2iter = v2.begin();
  while ( v1iter != v1.end() && v2iter != v2.end() )
  {
    if ( *v1iter < *v2iter )
      ++v1iter;
    else if ( *v2iter < *v1iter )
      ++v2iter;
    else
      result.push_back( make_pair(*v1iter++,*v2iter) );
  }
}

/// Intersection: copy all the elements in the first vector that also
/// occur in the second (which can result in duplicates).  Both
/// vectors must be SORTED.  There are two versions: the first accepts
/// as an argument a container for the result, the second simply
/// returns the result.  The elements in the result are taken from the
/// first vector.

template<class T, class R> void
Intersection( const vec<T>& v1, const vec<T>& v2, vec<R>& result )
{    
  result.clear();
  typename vec<T>::const_iterator v1iter = v1.begin();
  typename vec<T>::const_iterator v2iter = v2.begin();
  while ( v1iter != v1.end() && v2iter != v2.end() )
  {
    if ( *v1iter < *v2iter )
      ++v1iter;
    else if ( *v2iter < *v1iter )
      ++v2iter;
    else
      result.push_back( *v1iter++ );
  }
}

/// Intersection: copy all the elements in the first vector that also
/// occur in the second (which can result in duplicates).  Both
/// vectors must be SORTED.  Returns the result.  
/// The elements in the result are taken from the
/// first vector.

template<class T> vec<T> Intersection( const vec<T>& v1, const vec<T>& v2 )
{    
  vec<T> w;
  Intersection( v1, v2, w ); 
  return w;    
}

// Find the intersection of a family of sorted vectors.

template<class T> void Intersection( const vec< vec<T> >& x, vec<T>& y )
{    ForceAssert( x.nonempty( ) );
     y = x[0];
     for ( typename vec<T>::size_type j = 1; j < x.size(); j++ )
     {    vec<T> z;
          Intersection( x[j], y, z );
          y = z;    }    }

// MeetSize: assumes v1 and v2 are sorted.
// Handles duplicates weirdly.

template<class T> int64_t MeetSize( const vec<T>& v1, const vec<T>& v2 )
{ int64_t meet = 0;
  typename vec<T>::const_iterator v1iter = v1.begin();
  typename vec<T>::const_iterator v2iter = v2.begin();
  while ( v1iter != v1.end() && v2iter != v2.end() )
  {
    if ( *v1iter < *v2iter ) ++v1iter;
    else if ( *v2iter < *v1iter ) ++v2iter;
    else
    {
      meet++;
      ++v1iter;
      ++v2iter;
    }
  }
  return meet;
}

/// Meet: determine if two SORTED vectors have an element in common.  

template<class T> Bool Meet( const vec<T>& v1, const vec<T>& v2 )
{
  typename vec<T>::const_iterator v1iter = v1.begin();
  typename vec<T>::const_iterator v2iter = v2.begin();
  while ( v1iter != v1.end() && v2iter != v2.end() )
  {
    if ( *v1iter < *v2iter )
      ++v1iter;
    else if ( *v2iter < *v1iter )
      ++v2iter;
    else
      return true;
  }
  return false;
}

// Determine if a family of unsorted vectors have no elements in common.

template<class T> Bool Disjoint( const vec<vec<T>>& x )
{    vec< pair<T,int> > y;
     for ( int i = 0; i < x.isize( ); i++ )
     for ( int j = 0; j < x[i].isize( ); j++ )
          y.push( x[i][j], i );
     Sort(y);
     for ( int i = 0; i < y.isize( ) - 1; i++ )
     {    if ( y[i].first == y[i+1].first && y[i].second != y[i+1].second )
               return False;    }
     return True;    }

/// Meet2: determine if two UNSORTED vectors have an element in common.  

template<class T> Bool Meet2( vec<T> v1, vec<T> v2 )
{    Sort(v1), Sort(v2);
     return Meet( v1, v2 );    }

/// Invert a double vector of integral type.
/// For example, if you had a vector of edge ids for each read, calling this
/// would return a vector of read ids for each edge.  As a bonus, the read ids
/// are sorted.
template <class VVIn, class VVOut>
void invert( VVIn const& in, VVOut& out, size_t minOutSize=0 )
{
    typedef typename VVIn::value_type::value_type Value;
    Value maxVal = 0;
    Value minVal = 0;
    for ( auto const& vec : in )
        for ( Value val : vec )
            maxVal = std::max(maxVal,val), minVal = std::min(minVal,val);
    ForceAssertGe(minVal,Value(0));
    size_t outSize = std::max(size_t(maxVal+1),minOutSize);

    std::vector<size_t> counts(outSize,0ul);
    for ( auto const& vec : in )
        for ( auto val : vec )
            counts[val] += 1;

    out.clear();
    out.resize(outSize);
    auto oItr = out.begin();
    for ( size_t sz : counts )
        oItr->reserve(sz), ++oItr;

    size_t nnn = in.size();
    for ( size_t idx = 0; idx != nnn; ++idx )
        for ( auto val : in[idx] )
            out[val].push_back(idx);
}

template<class T> vec<T> Flatten( const vec<vec<T>>& x )
{    vec<T> y;
     for ( auto v : x ) y.append(v);
     return y;    }

inline vec<int> Contents( const vec<int>& x )
{    vec<int> y(x);
     UniqueSort(y);
     return y;    }

inline vec<int> Contents( const vec<vec<int>>& x )
{    vec<int> y;
     for ( auto& v : x ) y.append(v);
     UniqueSort(y);
     return y;    }

inline vec<int> Contents( const vec<vec<vec<int>>>& x )
{    vec<int> y;
     for ( auto& v : x ) y.append( Contents(v) );
     UniqueSort(y);
     return y;    }

#endif
