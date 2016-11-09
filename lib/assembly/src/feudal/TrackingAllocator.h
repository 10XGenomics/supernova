///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * TrackingAllocator.h
 *
 *  Created on: May 18, 2013
 *      Author: tsharpe
 */

#ifndef TRACKINGALLOCATOR_H_
#define TRACKINGALLOCATOR_H_

#include <scoped_allocator>
#include <memory>

class MemUse
{
public:
    MemUse( char const* type ) : mType(type) {}
    ~MemUse()
    { if ( mTotalAllocated > gMinReportSize ) report(); }

    size_t getNAllocs() const { return mNAllocs; }
    size_t getTotalAllocated() const { return mTotalAllocated; }
    size_t getCurrentInUse() const { return mInUse; }
    size_t getMaxUsed() const { return mMaxUsed; }

    void alloc( size_t nBytes )
    { mNAllocs += 1; mTotalAllocated += nBytes;
      if ( (mInUse += nBytes) > mMaxUsed ) mMaxUsed = mInUse; }

    void free( size_t nBytes )
    { mInUse -= nBytes; }

    void report();

    static void setMinReportSize( size_t minReportSize )
    { gMinReportSize = minReportSize; }

private:
    char const* mType;
    size_t mNAllocs = 0;
    size_t mTotalAllocated = 0;
    size_t mInUse = 0;
    size_t mMaxUsed = 0;

    static size_t gMinReportSize;
};

template <class T>
class TrackingAllocator
 : public std::scoped_allocator_adaptor<std::allocator<T>>
{
    typedef std::scoped_allocator_adaptor<std::allocator<T>> BaseT;
public:
    typedef unsigned long size_type;
    typedef long difference_type;
    typedef T* pointer;
    typedef T const* const_pointer;
    typedef T& reference;
    typedef T const& const_reference;
    typedef T value_type;

    TrackingAllocator() : mMemUse(new MemUse("TrackingAllocator")) {}

    template <class U>
    TrackingAllocator( TrackingAllocator<U> const& that )
    : mMemUse(that.mMemUse) {}

    TrackingAllocator( TrackingAllocator const& )=default;

    TrackingAllocator( TrackingAllocator&& a )
    : mMemUse(a.mMemUse) {} // it's actually a copy

    TrackingAllocator& operator=( TrackingAllocator const& )=default;

    TrackingAllocator& operator=( TrackingAllocator&& a )
    { mMemUse = a.mMemUse; return *this; } // it's actually a copy

    template <class U>
    struct rebind { typedef TrackingAllocator<U> other; };

    pointer allocate( size_type n, void* hint = 0 )
    { mMemUse->alloc(n*sizeof(T)); return BaseT::allocate(n,hint); }

    void deallocate( pointer p, size_type n )
    { mMemUse->free(n*sizeof(T)); BaseT::deallocate(p,n); }

    bool operator==( TrackingAllocator const& that )
    { return mMemUse == that.mMemUse; }

    bool operator!=( TrackingAllocator const& that )
    { return !(*this == that); }

private:
    std::shared_ptr<MemUse> mMemUse;
    template <class U> friend class TrackingAllocator;
};

template <class T>
struct DefaultAllocator
{
#ifndef TRACK_MEMUSE
typedef std::allocator<T> type;
#else
typedef TrackingAllocator<T> type;
#endif
};

#endif /* TRACKINGALLOCATOR_H_ */
