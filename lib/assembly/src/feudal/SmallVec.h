///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * \file SmallVec.h
 * \author tsharpe
 * \date Jul 7, 2009
 *
 * \brief Exactly like std::vector, but only 2^32 elements,
 * and a smaller memory footprint.
 */
#ifndef FEUDAL_SMALLVEC_H_
#define FEUDAL_SMALLVEC_H_

#include "Compare.h"
#include "feudal/Iterator.h"
#include "feudal/Oob.h"
#include "feudal/BinaryStream.h"
#include "feudal/IsSizeT.h"
#include "feudal/Generic.h"
#include "system/Assert.h"
#include "system/StaticAssert.h"
#include <algorithm>
#include <cstddef>
#include <initializer_list>
#include <limits>
#include <iostream>

/// Just like a std::vector, but with a smaller memory footprint.
template <class T, class A>
class SmallVec
{
public:
    typedef T value_type;
    typedef T& reference;
    typedef T const& const_reference;
    typedef T* pointer;
    typedef T const* const_pointer;
    typedef unsigned int size_type;
    typedef std::ptrdiff_t difference_type;
    typedef FwdIterator<T> iterator;
    typedef FwdConstIterator<T> const_iterator;
    typedef RevIterator<T> reverse_iterator;
    typedef RevConstIterator<T> const_reverse_iterator;

    SmallVec()
    { allocInit(A()); }

    explicit SmallVec( A const& alloc )
    { allocInit(alloc); }

    explicit SmallVec( size_type sz )
    { allocInit(A());
      reserve(sz);
      init(data()+sz); }

    explicit SmallVec( size_type sz,
                       T const& exemplar,
                       size_type cap = 0,
                       A const& alloc = A() )
    { allocInit(alloc);
      using std::max; reserve(max(sz,cap));
      init(data()+sz,exemplar); }

    template <class Itr>
    SmallVec( Itr first, Itr const& last,
              size_type cap = 0, A const& alloc = A() )
    { allocInit(alloc); assign(first,last,cap); }

    SmallVec( SmallVec const& that )
    { allocInit(A());
      reserve(that.size());
      copy(that.data(),that.dataEnd()); }

    SmallVec( SmallVec&& that )
    { move(that); }

    SmallVec( std::initializer_list<T> list )
    { assign(list.begin(),list.end()); }

    ~SmallVec()
    { UtilizationReporter::gInstance.report(this,8*sizeof(T),size(),capacity(),"SV");
      clear();
      deallocate();
      allocator().~A(); }

    template <class A1>
    SmallVec& operator=( SmallVec<T,A1> const& that )
    { if ( this != &that )
      { clear(); reserve(that.size()); copy(that.data(),that.dataEnd()); }
      return *this; }

    SmallVec& operator=( SmallVec const& that )
    { if ( this != &that )
      { clear(); reserve(that.size()); copy(that.data(),that.dataEnd()); }
      return *this; }

    SmallVec& operator=( SmallVec&& that )
    { if ( this != &that ) move(that);
      return *this; }

    /*
     * iterator
     */

    iterator begin() { return iterator(data()); }
    iterator begin( size_type idx ) { return iterator(data()+idx); }

    const_iterator begin() const { return const_iterator(data()); }
    const_iterator begin( size_type idx ) const
    { return const_iterator(data()+idx); }

    const_iterator cbegin() const { return const_iterator(data()); }
    const_iterator cbegin( size_type idx ) const
    { return const_iterator(data()+idx); }

    iterator end() { return iterator(dataEnd()); }
    const_iterator end() const { return const_iterator(dataEnd()); }
    const_iterator cend() const { return const_iterator(dataEnd()); }

    reverse_iterator rbegin() { return reverse_iterator(0,dataEnd()); }
    reverse_iterator rbegin( size_type idx )
    { return reverse_iterator(idx,dataEnd()); }

    const_reverse_iterator rbegin() const
    { return const_reverse_iterator(0,dataEnd()); }
    const_reverse_iterator rbegin( size_type idx ) const
    { return const_reverse_iterator(idx,dataEnd()); }

    const_reverse_iterator crbegin() const
    { return const_reverse_iterator(0,dataEnd()); }
    const_reverse_iterator crbegin( size_type idx ) const
    { return const_reverse_iterator(idx,dataEnd()); }

    reverse_iterator rend() { return reverse_iterator(size(),dataEnd()); }
    const_reverse_iterator rend() const
    { return const_reverse_iterator(size(),dataEnd()); }
    const_reverse_iterator crend() const
    { return const_reverse_iterator(size(),dataEnd()); }

    /*
     * counts and sizes
     */

    size_type size() const { return mSize; }
    bool nonempty( ) const { return mSize > 0; }
    bool solo( ) const { return mSize == 1; }
    size_type allocSize() const { return mSize; }

    SmallVec& resize( size_type sz )
    { reserve(sz);
      T* posn = data() + sz;
      T* end = dataEnd();
      if ( end > posn ) destroy(posn);
      if ( posn > end ) init(posn);
      return *this; }

    SmallVec& resize( size_type sz, T const& exemplar )
    { reserve(sz);
      T* posn = data() + sz;
      T* end = dataEnd();
      if ( end > posn ) destroy(posn);
      if ( posn > end ) init(posn,exemplar);
      return *this; }

    size_type max_size() const
    { return std::numeric_limits<size_type>::max() - 1000; }
    size_type capacity() const { return mCapacity; }
    bool empty() const { return !size(); }
    bool full() const { return mSize == mCapacity; }

    SmallVec& reserve( size_type sz )
    { if ( sz > capacity() ) realloc(sz); return *this; }


    /*
     * element access
     */

    reference operator[]( size_type idx )
    { AssertLt( idx, size() );
      return data()[idx]; }

    const_reference operator[]( size_type idx ) const
    { AssertLt( idx, size() );
      return data()[idx]; }

    reference at( size_type idx )
    { if ( idx >= size() ) OutOfBoundsReporter::oob("SmallVec",idx,size());
      return data()[idx]; }

    const_reference at( size_type idx ) const
    { if ( idx >= size() ) OutOfBoundsReporter::oob("SmallVec",idx,size());
      return data()[idx]; }

    reference front() { Assert(!empty()); return data()[0]; }
    const_reference front() const { Assert(!empty()); return data()[0]; }
    reference back() { Assert(!empty()); return dataEnd()[-1]; }
    const_reference back() const { Assert(!empty()); return dataEnd()[-1]; }

    /*
     * modifiers
     */
  
    void set(size_type idx, value_type val) { at(idx) = val; }

    SmallVec& assign( size_type sz )
    { return clear().resize(sz); }

    SmallVec& assign( size_type sz, T const& exemplar, size_type cap = 0 )
    { return assignConst(sz,exemplar,cap); }

    template <class Itr>
    SmallVec& assign( Itr first, Itr const& last, size_type cap = 0 )
    { return assignItr(first,last,cap,IsSizeT<Itr>()); }

    SmallVec& append( size_type sz )
    { resize(size()+sz);
      return *this; }

    SmallVec& append( size_type sz, T const& exemplar, size_type cap = 0 )
    { return appendConst(sz,exemplar,cap); }

    template <class Itr>
    SmallVec& append( Itr first, Itr const& last, size_type cap = 0 )
    { return appendItr(first,last,cap,IsSizeT<Itr>()); }

    SmallVec& push_front( T const& val )
    { insert(begin(),val); return *this; }

    SmallVec& push_back( T const& val )
    { return push_back(val,1.6f,64U); }

    SmallVec& push_back( T const& val, float growthFact, size_type growthIncr )
    { AssertGe(growthFact,0.f);
      if ( !full() )
      { init(dataEnd()+1,val); }
      else if ( isSelf(val) )
      { push_back(T(val),growthFact,growthIncr); }
      else
      { size_t sz = static_cast<size_t>(growthFact*size()+growthIncr);
        using std::min; using std::max;
        sz = min(sz,static_cast<size_t>(max_size()));
        realloc(max(static_cast<size_t>(size()+1ul),sz));
        init(dataEnd()+1,val); }
      return *this; }

    iterator erase( iterator pos )
    { T* dest = &*pos;
      Assert(data() <= dest && dest < dataEnd());
      destroy(swap(dest+1,dest));
      return pos; }

    iterator erase( iterator first, iterator last )
    { T* dest = &*first;
      T* start = &*last;
      Assert(data() <= dest && dest <= start && start <= dataEnd());
      destroy(swap(start,dest));
      return first; }

    iterator insert(iterator position, const T& x)
    { size_type offset = (&*position) - data();
      AssertLe(offset,size());
      size_type end = size();
      push_back(x);
      while(end != offset)
      { swap_back_single(data() + end, 1); --end; }
      return begin(offset); }

    void insert(iterator position, size_type n, const T& x)
    { size_type offset = (&*position) - data();
      AssertLe(offset,size());
      size_type end = size();
      append(n, x);
      while((end - offset) >= n)
      { swap_back_group(data() + end, n); end -= n; }
      while(end != offset)
      { swap_back_single(data() + end, n); --end; } }

    template <class InputIterator>
    void insert(iterator position, InputIterator first, InputIterator last)
    { using std::distance;
      size_type offset = (&*position) - data();
      AssertLe(offset,size());
      size_type end = size();
      size_type n = static_cast<size_type>(distance(first, last));
      append(first, last);
      while((end - offset) >= n)
      { swap_back_group(data() + end, n); end -= n; }
      while(end != offset)
      { swap_back_single(data() + end, n); --end; } }

    /// erases each element for which *Itr is true
    template <class Itr>
    SmallVec& eraseIf( Itr begin, Itr const& end )
    { using std::distance; using std::swap;
      AssertEq(static_cast<size_type>(distance(begin,end)),size());
      T* dest = data();
      for ( T* src = dest; begin != end; ++begin, ++src )
      { if ( !*begin ) { if ( src != dest ) swap(*src,*dest); ++dest; } }
      destroy(dest);
      return *this; }

  /// erases each element for which *Itr is false
    template <class Itr>
    SmallVec& eraseUnless( Itr begin, Itr const& end )
    { using std::distance; using std::swap;
      AssertEq(static_cast<size_type>(distance(begin,end)),size());
      T* dest = data();
      for ( T* src = dest; begin != end; ++begin, ++src )
      { if ( *begin ) { if ( src != dest ) swap(*src,*dest); ++dest; } }
      destroy(dest);
      return *this; }

    SmallVec& pop_back()
    { Assert(!empty()); destroy(dataEnd()-1); return *this; }

    SmallVec& clear() { destroy(data()); return *this; }
    SmallVec& shrink_to_fit() { if ( !full() ) realloc(size()); return *this; }

    SmallVec& swap( SmallVec& that )
    { if ( allocator() != that.allocator() ) // if we have to do a slow swap
        exchange(that);
      else
        quickSwap(that);
      return *this; }

    friend int compare( SmallVec const& v1, SmallVec const& v2 )
    { const_pointer itr1 = v1.data();
      const_pointer itr2 = v2.data();
      using std::min; const_pointer end = itr1 + min(v1.size(),v2.size());
      int result = 0;
      while ( !result && itr1 != end ) result = compare(*itr1++,*itr2++);
      if ( !result ) result = compare(v1.size(),v2.size());
      return result; }

    A get_allocator() const { return allocator(); }

    void writeFeudal( BinaryWriter& writer, void const** ppFixed ) const
    { *ppFixed = &mSize;
      writer.write(data(),dataEnd()); }

    void readFeudal( BinaryReader& reader, unsigned long dataLen, void* fixed );

    static unsigned int fixedDataLen()
    { return BinaryReader::externalSizeof(static_cast<T*>(nullptr)) ?
                0U : // for compatibility with existing feudal files of primitive T's
                sizeof(size_type); } // for efficiency with complex T's

    void writeBinary( BinaryWriter& writer ) const
    { writer.write(mSize);
      writer.write(data(),dataEnd()); }

    void readBinary( BinaryReader& reader )
    { size_type sz; reader.read(&sz); resize(sz);
      reader.read(data(),dataEnd()); }

    static size_t externalSizeof() { return 0; }

    static size_type interpretSize( void* pFixed, size_t varDataLen )
    { size_type sz;
      if ( fixedDataLen() )
        memcpy(&sz,pFixed,sizeof(sz));
      else
        sz = varDataLen/BinaryReader::externalSizeof(static_cast<T*>(nullptr));
      return sz; }

protected:
    T* data() { return reinterpret_cast<T*>(mData&MASK); }
    T const* data() const { return reinterpret_cast<T*>(mData&MASK); }
    T* dataEnd() { return data()+size(); }
    T const* dataEnd() const { return data()+size(); }
    void setSize( size_type sz ) { mSize = sz; }

private:
    void setData( void const* data )
    { mData ^= (mData ^ reinterpret_cast<unsigned long>(data)) & MASK; }

    void allocInit( A const& alloc )
    { STATIC_ASSERT(sizeof(A) <= 2u); new (&allocator()) A(alloc); }

    A& allocator() { return reinterpret_cast<A*>(&mData+1)[-1]; }
    A const& allocator() const
    { return reinterpret_cast<A const*>(&mData+1)[-1]; }

    bool isSelf( T const& val ) { return data() <= &val && &val < dataEnd(); }

    void init( T* last );
    void init( T* last, T const& exemplar );
    void destroy( T* first );
    T* swap( T* first, T* dest );

    void swap_back_group(pointer pos, size_type len)
    { using std::swap;
      for(size_type idx = 0; idx < len; ++idx)
        swap(*(idx + pos), *(idx + pos - len)); }

    void swap_back_single(pointer pos, size_type len)
    { using std::swap;
      for(size_type idx = 0; idx < len; ++idx)
        swap(*(pos - 1 + idx), *(pos + idx)); }

    void quickSwap( SmallVec& that )
    { using std::swap;
      swap(mData,that.mData);
      swap(mSize,that.mSize);
      swap(mCapacity,that.mCapacity); }

    void move( SmallVec& that )
    { if ( allocator() != that.allocator() )
        *this = that;
      else
        quickSwap(that); }

    template <class Itr>
    void copy( Itr first, Itr const& last )
    { T* initialDest = dataEnd();
      T* dest = initialDest;
      A& alloc = allocator();
      while ( first != last )
      { alloc.construct(dest++,*first);
        ++first; }
      mSize += dest-initialDest; }

    SmallVec& assignConst( size_type sz, T const& exemplar, size_type cap )
    { using std::max; size_type reserveSize = max(sz,cap);
      if ( reserveSize > capacity() && isSelf(exemplar) )
      { assignConst(sz,T(exemplar),cap); }
      else
      { clear();
        reserve(reserveSize);
        init(data()+sz,exemplar); }
      return *this; }

    template <class Itr>
    SmallVec& assignItr( Itr first, Itr const& last, size_type cap, NoSizeT )
    { Assert(first==last||!isSelf(*first));
      clear();
      using std::distance; size_t dist = distance(first,last);
      using std::max; reserve(max(dist,static_cast<size_t>(cap)));
      copy(first,last);
      return *this; }

    template <class Itr>
    SmallVec& assignItr( Itr first, Itr const& last, size_type cap, YesSizeT )
    { return assignConst(first,last,cap); }

    SmallVec& appendConst( size_type sz, T const& exemplar, size_type cap )
    { using std::max;
      size_t reserveSize = max(size()+static_cast<size_t>(sz),
                               static_cast<size_t>(cap));
      if ( reserveSize > capacity() && isSelf(exemplar) )
      { appendConst(sz,T(exemplar),cap); }
      else
      { reserve(reserveSize);
        init(dataEnd()+sz,exemplar); }
      return *this; }

    template <class Itr>
    SmallVec& appendItr( Itr first, Itr const& last, size_type cap, NoSizeT )
    { Assert(first==last||!isSelf(*first));
      using std::distance;
      size_t newSize = size() + static_cast<size_t>(distance(first,last));
      using std::max; reserve(max(newSize,static_cast<size_t>(cap)));
      copy(first,last);
      return *this; }

    template <class Itr>
    SmallVec& appendItr( Itr first, Itr const& last, size_type cap, YesSizeT )
    { return appendConst(first,last,cap); }

    void realloc( size_t sz );
    void exchange( SmallVec& that );

    void deallocate()
    { if ( data() ) allocator().deallocate(data(),capacity()); }

    size_type mSize = 0;
    size_type mCapacity = 0;
    unsigned long mData = 0;
    static unsigned long const MASK = 0xffffffffffff; // lowest 48 bits

    template <class charT,class Traits> friend class FeudalString;
};

template <class T, class A>
struct Serializability< SmallVec<T,A> >
{ typedef SelfSerializable type; };

template <class T, class A1, class A2>
bool operator==( SmallVec<T,A1> const& v1, SmallVec<T,A2> const& v2 )
{ using std::equal;
  return v1.size()==v2.size() && equal(v1.begin(),v1.end(),v2.begin()); }

template <class T, class A1, class A2>
bool operator!=( SmallVec<T,A1> const& v1, SmallVec<T,A2> const& v2 )
{ using std::equal;
  return v1.size()!=v2.size() || !equal(v1.begin(),v1.end(),v2.begin()); }

template <class T, class A1, class A2>
bool operator<( SmallVec<T,A1> const& v1, SmallVec<T,A2> const& v2 )
{ using std::lexicographical_compare;
  return lexicographical_compare(v1.begin(),v1.end(),v2.begin(),v2.end()); }

template <class T, class A1, class A2>
bool operator>( SmallVec<T,A1> const& v1, SmallVec<T,A2> const& v2 )
{ return (v2 < v1); }

template <class T, class A1, class A2>
bool operator<=( SmallVec<T,A1> const& v1, SmallVec<T,A2> const& v2 )
{ return !(v2 < v1); }

template <class T, class A1, class A2>
bool operator>=( SmallVec<T,A1> const& v1, SmallVec<T,A2> const& v2 )
{ return !(v1 < v2); }

template <class T, class A>
void swap( SmallVec<T,A>& v1, SmallVec<T,A>& v2 )
{ v1.swap(v2); }

template<class T, class A>
inline void Sort( SmallVec<T,A>& v )
{    sort( v.begin( ), v.end( ) );    }

#endif /* FEUDAL_SMALLVEC_H_ */
