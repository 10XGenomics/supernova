///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// This file is a combination of part of Vec.h and VecUtilities.h, modified to 
// use gcc parallelized algorithm versions.  At present, only sort functions are
// included, and only the parallelized algorithms are the sorts per se.  Probably
// more could be parallelized.

#ifndef PARALLEL_VEC_UTILITIES_H
#define PARALLEL_VEC_UTILITIES_H

// MakeDepend: cflags OMP_FLAGS

#include <parallel/algorithm>
#include <functional>

#include "Vec.h"
#include "VecUtilities.h"

template<class T> void ParallelSort( vec<T>& v )
{
  TRACEVAL_STOP_TRACING_COPIES;
  __gnu_parallel::sort( v.begin( ), v.end( ), std::less<T>() );
  TRACEVAL_START_TRACING_COPIES;
}

template<class T, class StrictWeakOrdering > 
inline void ParallelSort( vec<T>& v, StrictWeakOrdering comp )
{
  TRACEVAL_STOP_TRACING_COPIES;
  __gnu_parallel::sort( v.begin( ), v.end( ), comp );
  TRACEVAL_START_TRACING_COPIES;
}

template<class T> void ParallelReverseSort( vec<T>& v )
{
  TRACEVAL_STOP_TRACING_COPIES;
  __gnu_parallel::sort( v.rbegin( ), v.rend( ), std::less<T>() );
  TRACEVAL_START_TRACING_COPIES;
}

template<class T, class StrictWeakOrdering> 
void ParallelReverseSort( vec<T>& v, StrictWeakOrdering comp )
{
  TRACEVAL_STOP_TRACING_COPIES;
  __gnu_parallel::sort( v.rbegin( ), v.rend( ), comp );
  TRACEVAL_START_TRACING_COPIES;
}

template <class T, class StrictWeakOrdering, class EqualPredicate> 
  void ParallelUniqueSort(vec<T> & v, StrictWeakOrdering comp, EqualPredicate equal) {
  if ( v.empty() ) return;

  TRACEVAL_STOP_TRACING_COPIES;
  
  ParallelSort( v, comp );
  
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
    if ( equal(*iter,*last_unique_elem_iter ) ) continue;
    *(++move_to) = *(last_unique_elem_iter = iter);
  }
  v.erase(++move_to,iter_end);

  TRACEVAL_START_TRACING_COPIES;
}

template <class T> inline void ParallelUniqueSort(vec<T> & v) {
  ParallelUniqueSort(v,less<T>(),equal_to<T>());
}

template <class T, class StrictWeakOrdering, class EqualPredicate> 
void ParallelUniqueSortAndCount(vec<T> & v, StrictWeakOrdering comp, 
     EqualPredicate equal, vec<unsigned int> & counts) {
  counts.clear();
  if ( v.empty() ) return;

  TRACEVAL_STOP_TRACING_COPIES;
  
  ParallelSort( v, comp );
  
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

template <class T, class StrictWeakOrdering, class EqualPredicate, class FilterPredicate> 
void ParallelUniqueSortAndCount(vec<T> & v, StrictWeakOrdering comp, 
     EqualPredicate equal, vec<unsigned int> & counts, FilterPredicate pass) {
  counts.clear();
  if ( v.empty() ) return;

  TRACEVAL_STOP_TRACING_COPIES;
  
  ParallelSort( v, comp );
  
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

template <class T> 
inline  void ParallelUniqueSortAndCount(vec<T> & v, vec<unsigned int> & counts) {
  ParallelUniqueSortAndCount(v,less<T>(),equal_to<T>(),counts);
}

template <class T, class FilterPredicate> 
inline  void ParallelUniqueSortAndCount(vec<T> & v, vec<unsigned int> & counts, 
     FilterPredicate pass) {
  ParallelUniqueSortAndCount(v, less<T>(), equal_to<T>(), counts, pass);
}

template<class V, typename C, typename IDX>
void ParallelWhatPermutation( const V & v, vec<IDX>& permutation, C comp,
                                bool inv = true )
{
    AssertLe(v.size(), static_cast<size_t>(numeric_limits<IDX>::max()));
    IDX n = v.size();
    vec<IDX> perm(n);
    for ( IDX i = 0; i < n; i++ )
        perm[i] = i;
    __gnu_parallel::sort(perm.begin(),perm.end(),indirect_compare<V,C>(v,comp));
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
void ParallelWhatPermutation( const V & v, vec<IDX>& permutation,
                                bool inv = true )
{
    ParallelWhatPermutation(v,permutation,less<typename V::value_type>(),inv);
}

template<class V1, class V2, typename C, typename IDX>
void ParallelSortSyncDispatch( V1& v, V2& w, C comp, IDX )
{
    ForceAssertEq( v.size(), w.size() );
    vec<IDX> perm;
    ParallelWhatPermutation(v, perm, comp);
    PermuteVec(v, perm);
    PermuteVec(w, perm);
}

template<class V1, class V2, typename C=std::less<typename V1::value_type>>
void ParallelSortSync( V1& v, V2& w, C comp= C(),
                decltype(comp(v[0],v[0])) = false, decltype(w.size()) = 0 )
{
    if ( v.size() < std::numeric_limits<unsigned>::max() )
        ParallelSortSyncDispatch(v,w,comp,0u);
    else
        ParallelSortSyncDispatch(v,w,comp,0ul);
}

template<class V1, class V2>
void ParallelReverseSortSync( V1& v, V2& w )
{
    ParallelSortSync(v,w,std::greater<typename V1::value_type>());
}

template<class V1, class V2>
void ParallelUniqueSortSync( V1& v, V2& w )
{
    ParallelSortSync(v, w, less<typename V1::value_type>());
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

template<class V1, class V2, class V3, typename C, typename IDX>
void ParallelSortSyncDispatch( V1& v, V2& w, V3& x, C comp, IDX )
{
    ForceAssertEq( v.size(), w.size() );
    ForceAssertEq( v.size(), x.size() );
    vec<IDX> perm;
    ParallelWhatPermutation(v, perm, comp);
    PermuteVec(v, perm);
    PermuteVec(w, perm);
    PermuteVec(x, perm);
}

template<class V1, class V2, class V3, typename C=std::less<typename V1::value_type>>
void ParallelSortSync( V1& v, V2& w, V3& x, C comp = C(),
                decltype(comp(v[0],v[0])) = false, decltype(x.size()) = 0 )
{
    if ( v.size() < std::numeric_limits<unsigned>::max() )
        ParallelSortSyncDispatch(v,w,x,comp,0u);
    else
        ParallelSortSyncDispatch(v,w,x,comp,0ul);
}

template<class V1, class V2, class V3>
void ParallelReverseSortSync( V1& v, V2& w, V3& x )
{
    ParallelSortSync(v,w,x,greater<typename V1::value_type>());
}

template<class V1, class V2, class V3, class V4, typename C, typename IDX>
void ParallelSortSyncDispatch( V1& v, V2& w, V3& x, V4& y, C comp, IDX )
{
    ForceAssertEq( v.size(), w.size() );
    ForceAssertEq( v.size(), x.size() );
    ForceAssertEq( v.size(), y.size() );
    vec<IDX> perm;
    ParallelWhatPermutation(v, perm, comp);
    PermuteVec(v, perm);
    PermuteVec(w, perm);
    PermuteVec(x, perm);
    PermuteVec(y, perm);
}

template<class V1, class V2, class V3, class V4, typename C=std::less<typename V1::value_type>>
void ParallelSortSync( V1& v, V2& w, V3& x, V4& y, C comp = C(),
                decltype(comp(v[0],v[0])) = false, decltype(y.size()) = 0 )
{
    if ( v.size() <= std::numeric_limits<unsigned>::max() )
        ParallelSortSyncDispatch(v,w,x,y,comp,0u);
    else
        ParallelSortSyncDispatch(v,w,x,y,comp,0ul);
}

template<class V1, class V2, class V3, class V4>
void ParallelReverseSortSync( V1& v, V2& w, V3& x, V4& y )
{
    ParallelSortSync(v,w,x,y,greater<typename V1::value_type>());
}

template<class V1, typename IDX, typename C = std::less<typename V1::value_type>>
void ParallelSortIndex( const V1& v, vec<IDX>& index, C comp = C() )
{
    AssertLe(v.size(), static_cast<size_t>(std::numeric_limits<IDX>::max()));
    index.resize(v.size());
    iota(index.begin(), index.end(), 0);
    vec<IDX> perm;
    ParallelWhatPermutation(v, perm, comp);
    PermuteVec(index, perm);
}

#endif
