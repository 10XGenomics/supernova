/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2010) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////
/*
 * \file WorklistN.h
 * \author tsharpe
 * \date Apr. 14, 2010
 *
 * \brief Use multiple threads of execution to process a number of workitems.
 * This is just like the more general Worklist, but the workitems are simply
 * numbered, rather than existing as Workitem objects in a list.
 */
#ifndef SYSTEM_WORKLISTN_H_
#define SYSTEM_WORKLISTN_H_

#include "system/Assert.h"
#include "system/SysConf.h"
#include "system/WorklistUtils.h"
#include <algorithm>
#include <cstddef>
#include <list>
#include <time.h>

/// See the documentation for Worklist.  This works the same way.
template<class Processor, class Index=size_t>
class WorklistN
{
public:
    struct ProgressMonitor
    {
        ProgressMonitor() : mTime(0) {}
        Index mWorkitem;
        time_t mTime;
    };

    /// Create a Worklist that will be processed by nThreads of execution.
    /// If nThreads is 0, no threads will be created, and the entire
    /// Worklist will be processed by the Worklist destructor.  The default
    /// value for nThreads is the number of processors less one (we assume that
    /// the main thread will continue to do work).
    /// Passing a stack size of 0 (it's default value) will cause threads to be
    /// created with a stack as large as main thread's stack (unless the main
    /// thread's stack size is unlimited, in which case a default value of 8Mb
    /// will be used for each thread's stack size).
    WorklistN( Index start,
               Index end,
               Processor const& p,
               size_t nThreads,
               size_t threadStackSize = 0 )
    : mCurWorkitem(start), mLastWorkitem(end),
      mNThreads(calcThreads(nThreads,end-start)), mProcessor(p),
      mProgressMonitors(new ProgressMonitor[mNThreads]),
      mTP(mNThreads,threadStackSize,threadFunc,this)
    { QueueStateManipulator qsm(mQS); qsm.setSize(end-start); }

    WorklistN( WorklistN const& )=delete;
    WorklistN& operator=( WorklistN const& )=delete;

    /// This will wait until all work is done, and all threads have
    /// been terminated.  That might be a long time!
    ~WorklistN()
    { waitForDone();
      delete [] mProgressMonitors; }

    /// Start up a separate thread to monitor thread death.
    /// This thread wakes up every howOftenToCheckInSecs seconds, and replaces
    /// any threads that have died with a new thread.
    /// If you call this more than once, subsequent calls merely adjust the
    /// check interval.
    /// If you call this with howOftenToCheckInSecs == 0, the monitor exits.
    /// If reportRestarts is true, an error message is written to cerr to inform
    /// you of a thread restart.
    void threadDeathMonitor( long howOftenToCheckInSecs,
                             bool reportRestarts = true )
    { mTP.threadDeathMonitor(howOftenToCheckInSecs,reportRestarts); }

    /// Return the number of unprocessed workitems.
    size_t size()
    { QueueStateManipulator qsm(mQS);
      return qsm.getSize(); }

    /// Note:  you don't have to call this -- the destructor will automatically
    /// wait until all the work is complete.  You may call this for
    /// entertainment value, if you'd like.
    void waitForDone()
    { mTP.threadDeathMonitor(0); // shut down the monitor.
      QueueStateManipulator(mQS).setDone(); // tell threads to quit
      mTP.shutdown(); }

    /// Clear the worklist of all unprocessed items.
    /// If you need to quit early, this is the only way to do it.
    /// The destructor is still going to wait for in-progress work to complete.
    void clear()
    { QueueStateManipulator qsm(mQS);
      mCurWorkitem = mLastWorkitem;
      qsm.setSize(0); }

    /// How many threads are there to process the work.
    size_t getNThreads()
    { return mTP.getNThreads(); }

    /// What happening with a thread.
    /// ThreadIdx runs from 0 to getNThreads()-1.
    ProgressMonitor getProgress( size_t threadIdx )
    { QueueStateManipulator qsm(mQS);
      return mProgressMonitors[threadIdx]; }

private:
    void doWork();
    static void* threadFunc( void* );

    static size_t calcThreads( size_t nThreads, size_t nWorkitems )
    { return std::min(size_t(boundNumThreads(nThreads)),nWorkitems); }

    Index mCurWorkitem;
    Index mLastWorkitem;
    size_t mNThreads;
    QueueState mQS;
    Processor mProcessor;
    ProgressMonitor* mProgressMonitors;
    ThreadPool mTP; // this needs to be the last member because we don't
                    // want to start threads until everything else is ready
};

template<class Processor, class Index>
void WorklistN<Processor,Index>::doWork()
{
    size_t idx = mTP.findThreadIndex(pthread_self());
    ProgressMonitor& progressMonitor = mProgressMonitors[idx];
    Processor processor(mProcessor);

    while ( true )
    {
        if ( true ) // empty block just to control lifetime of qsm
        {
            QueueStateManipulator qsm(mQS);
            qsm.waitForWork();
            if ( !qsm.getSize() ) // wait for work only returns when there's work, or when the done flag is set
                break;            // so break if there's no work (because that means we're all done)

            progressMonitor.mWorkitem = mCurWorkitem;
            ++mCurWorkitem;
            qsm.decSize();
        }
        time(&progressMonitor.mTime);
        processor(progressMonitor.mWorkitem);
    }
}

template<class Processor, class Index>
void* WorklistN<Processor,Index>::threadFunc( void* ptr )
{
    reinterpret_cast<WorklistN*>(ptr)->doWork();
    return 0;
}

template <class Proc, class Index>
void parallelFor( Index start, Index end, Proc const& proc,
                    size_t nThreads = getConfiguredNumThreads() )
{
    nThreads = boundNumThreads(nThreads);
    if ( nThreads > 1 )
    {
        WorklistN<Proc> wl(start,end,proc,nThreads);
    }
    else
    {
        Proc prc(proc);
        for ( Index idx = start; idx != end; ++idx )
            prc(idx);
    }
}

template <class Proc, class Index>
void parallelForBatch( Index start, Index stop, size_t batchSize,
                        Proc proc, size_t nThreads = getConfiguredNumThreads(),
                        bool verbose = false )
{
    size_t nBatches = ((stop - start) + batchSize - 1) / batchSize;
    if ( verbose )
    {
        Dotter dotter(nBatches);
        parallelFor( 0ul, nBatches,
            [start,stop,batchSize,proc,&dotter]( size_t batchNo ) mutable
            { auto beg=start+batchNo*batchSize;
              size_t remain = stop-beg;
              auto end=beg+std::min(remain,batchSize);
              while ( beg != end )
              { proc( beg ); ++beg; }
              dotter.batchDone(); }, nThreads );
    }
    else
    {
        parallelFor( 0ul, nBatches,
            [start,stop,batchSize,proc]( size_t batchNo ) mutable
            { auto beg=start+batchNo*batchSize;
              size_t remain = stop-beg;
              auto end=beg+std::min(remain,batchSize);
              while ( beg != end )
              { proc( beg ); ++beg; } }, nThreads );
    }
}

#endif /* SYSTEM_WORKLISTN_H_ */
