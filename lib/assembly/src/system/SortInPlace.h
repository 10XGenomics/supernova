///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file SortInPlace.h
 * \author tsharpe
 * \date Feb 23, 2010
 *
 * \brief
 */
#ifndef SORTINPLACE_H_
#define SORTINPLACE_H_

#include "Compare.h"
#include "system/LockedData.h"
#include "system/SysConf.h"
#include "system/Worklist.h"
#include <algorithm>
#include <cstddef>
#include <utility>
#include <unistd.h>

// single-threaded, in-place sorter
template <class Itr, class Comp>
class InPlaceSorter
{
public:
    typedef typename std::iterator_traits<Itr>::difference_type diff_t;
    typedef typename std::iterator_traits<Itr>::value_type value_t;

    InPlaceSorter( Comp const& comp = Comp() )
    : mComp(comp) {}

    void sort( Itr const& first, Itr const& last )
    { internalSort(first,last); }

private:
    // choose a pivot -- it's the middle value
    Itr iter_median( Itr a, Itr b, Itr c )
    {
        int order = mComp(*a,*b);
        int order2 = mComp(*b,*c);
        if ( !order || !order2 ) return b;
        if ( order < 0 )
        {
            if ( order2 <= 0 ) return b;
            return mComp(*a,*c) <= 0 ? c : a;
        }
        else
        {
            if ( order2 >= 0 ) return b;
            return mComp(*a,*c) <= 0 ? a : c;
        }
    }

    // sort some range
    void internalSort( Itr first, Itr last )
    {
        using std::iter_swap;
        while ( true )
        {
            // use insertion sort for tiny ranges
            diff_t siz = last - first;
            if ( siz <= INSERTION_SORT_MAX )
            {
                if ( siz > 1 )
                {
                    for ( Itr itr1(first+1); itr1 != last; ++itr1 )
                    {
                        Itr itr2(itr1);
                        while ( itr2 != first )
                        {
                            Itr itr3(itr2--);
                            if ( mComp(*itr3,*itr2) >= 0 )
                                break;

                            iter_swap(itr3,itr2);
                        }
                    }
                }
                break;
            }

            Itr pivot(iter_median(first,first+siz/2,last-1));
            // move up over early elements that are less than the pivot
            // (big win when input is already in order)
            Itr lt(first);
            while ( mComp(*lt,*pivot) < 0 )
                ++lt;
            // ditto for elements greater than the pivot at the end of the range
            Itr gt(last);
            while ( mComp(*pivot,*--gt) < 0 )
                ;
            ++gt;
            // invariants: elements in range [first,lt) are less than the pivot
            //             elements in range [lt,itr) are equal to the pivot
            //             elements in range [itr,gt) are unclassified
            //             elements in range [gt,last) are greater than the pivot

            // move itr from lt to gt, keeping the invariants true
            Itr itr(lt);
            while ( itr != gt )
            {
                int order = mComp(*itr,*pivot);
                if ( order < 0 )
                {
                    iter_swap(lt,itr);
                    if ( pivot == lt )
                        pivot = itr;
                    ++lt;
                    ++itr;
                }
                else if ( order > 0 )
                {
                    iter_swap(itr,--gt);
                    if ( pivot == gt )
                    {
                        pivot = itr;
                        ++itr;
                    }
                }
                else
                {
                    pivot = itr;
                    ++itr;
                }
            }

            if ( lt-first > 1 )
                internalSort(first,lt);

            first = gt;
        }
    }

    Comp mComp;

    static diff_t const INSERTION_SORT_MAX = 8;
};

// same as above, but uses a worklist to handle one or the other part of the
// range once it has been pivoted.  worklist items are sub-ranges that need
// sorting.
template <class Itr, class Comp>
class InPlaceParallelSorter
{
public:
    typedef typename std::iterator_traits<Itr>::difference_type diff_t;
    typedef std::pair<Itr,Itr> Workitem;

    InPlaceParallelSorter( unsigned nThreads, Comp const& comp = Comp() )
    : mComp(comp),
      mWorklist(Processor(this),nThreads) {}

    void sort( Itr const& first, Itr const& last )
    {
        diff_t nToDo = last - first;
        if ( nToDo > 1 )
        {
            mDoneReporter.init(nToDo);
            mWorklist.add(Workitem(first,last));
            mDoneReporter.wait();
        }
    }

private:
    Itr iter_median( Itr a, Itr b, Itr c )
    {
        int order = mComp(*a,*b);
        int order2 = mComp(*b,*c);
        if ( !order || !order2 ) return b;
        if ( order < 0 )
        {
            if ( order2 <= 0 ) return b;
            return mComp(*a,*c) <= 0 ? c : a;
        }
        else
        {
            if ( order2 >= 0 ) return b;
            return mComp(*a,*c) <= 0 ? a : c;
        }
    }

    void internalSort( Itr first, Itr last )
    {
        using std::iter_swap;
        size_t done = 0;
        while ( true )
        {
            diff_t siz = last - first;
            if ( siz <= INSERTION_SORT_MAX )
            {
                if ( siz > 1 )
                {
                    for ( Itr itr1(first+1); itr1 != last; ++itr1 )
                    {
                        Itr itr2(itr1);
                        while ( itr2 != first )
                        {
                            Itr itr3(itr2--);
                            if ( mComp(*itr3,*itr2) >= 0 )
                                break;

                            iter_swap(itr3,itr2);
                        }
                    }
                }
                done += siz;
                break;
            }

            Itr pivot(iter_median(first,first+siz/2,last-1));
            Itr lt(first);
            while ( mComp(*lt,*pivot) < 0 )
                ++lt;
            Itr gt(last);
            while ( mComp(*pivot,*--gt) < 0 )
                ;
            ++gt;
            Itr itr(lt);
            while ( itr != gt )
            {
                int order = mComp(*itr,*pivot);
                if ( order < 0 )
                {
                    iter_swap(lt,itr);
                    if ( pivot == lt )
                        pivot = itr;
                    ++lt;
                    ++itr;
                }
                else if ( order > 0 )
                {
                    iter_swap(itr,--gt);
                    if ( pivot == gt )
                    {
                        pivot = itr;
                        ++itr;
                    }
                }
                else
                {
                    pivot = itr;
                    ++itr;
                }
            }

            done += gt - lt;
            diff_t siz1 = lt - first;
            diff_t siz2 = last - gt;
            if ( siz1 < siz2 )
            {
                if ( siz1 >= PARALLEL_SORT_MIN )
                {
                    mWorklist.add(Workitem(first,lt));
                    siz1 = 0;
                }
            }
            else if ( siz2 >= PARALLEL_SORT_MIN )
            {
                mWorklist.add(Workitem(gt,last));
                siz2 = 0;
            }

            if ( siz1 <= 1 )
                done += siz1;
            else
                internalSort(first,lt);

            if ( siz2 <= 1 )
            {
                done += siz2;
                break;
            }

            first = gt;
        }

        mDoneReporter.done(done);
    }

    class DoneReporter : private LockedData
    {
    public:
        DoneReporter() : mToDo(0), mDone(0), mCondVar(*this) {}

        // compiler-supplied destructor is ok

        void init( size_t toDo )
        { mToDo = toDo; mDone = 0; }

        size_t progress()
        { Locker locker(*this); return mDone; }

        void done( size_t siz )
        { Locker locker(*this);
          if ( (mDone += siz) >= mToDo )
              mCondVar.signal(); }

        void wait()
        { Locker locker(*this);
          if ( mDone < mToDo ) locker.wait(mCondVar); }

    private:
        DoneReporter( DoneReporter const& ); // unimplemented -- no copying
        DoneReporter& operator=( DoneReporter const& ); // unimplemented -- no copying

        size_t mToDo;
        size_t mDone;
        Condition mCondVar;
    };

    class Processor
    {
    public:
        Processor( InPlaceParallelSorter* pSorter ) : mpSorter(pSorter) {}

        // compiler-supplied copying and destructor are OK

        void operator()( Workitem const& wi ) const
        { mpSorter->internalSort(wi.first,wi.second); }

    private:
        InPlaceParallelSorter* mpSorter;
    };
    friend class Processor;

    Comp mComp;
    DoneReporter mDoneReporter;
    Worklist<Workitem,Processor> mWorklist;

    static diff_t const INSERTION_SORT_MAX = 8;
    static diff_t const PARALLEL_SORT_MIN = 1000;
};

template <class Itr, class Comp>
void sortInPlace( Itr const& first, Itr const& last, Comp const& comp )
{
    InPlaceSorter<Itr,Comp> sorter(comp);
    sorter.sort(first,last);
}

template <class Itr>
void sortInPlace( Itr const& first, Itr const& last )
{
    typedef typename std::iterator_traits<Itr>::value_type value_t;
    typedef CompareFunctor<value_t> Comp;
    InPlaceSorter<Itr,Comp> sorter;
    sorter.sort(first,last);
}

template <class Itr, class Comp>
void sortInPlaceParallel( Itr const& first, Itr const& last, Comp const& comp,
                            unsigned nThreads = 0 )
{
    nThreads = boundNumThreads(nThreads);
    if ( nThreads == 1 )
        sortInPlace(first,last,comp);
    else
    {
        InPlaceParallelSorter<Itr,Comp> sorter(nThreads,comp);
        sorter.sort(first,last);
    }
}

template <class Itr>
void sortInPlaceParallel( Itr const& first, Itr const& last,
                            unsigned nThreads = 0 )
{
    nThreads = boundNumThreads(nThreads);
    if ( nThreads == 1 )
        sortInPlace(first,last);
    else
    {
        typedef typename std::iterator_traits<Itr>::value_type value_t;
        typedef CompareFunctor<value_t> Comp;
        InPlaceParallelSorter<Itr,Comp> sorter(nThreads);
        sorter.sort(first,last);
    }
}

#endif /* SORTINPLACE_H_ */
