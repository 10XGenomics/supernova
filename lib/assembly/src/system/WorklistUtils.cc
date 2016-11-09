///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * \file WorklistUtils.cc
 * \author tsharpe
 * \date Feb 4, 2009
 *
 * \brief A couple of utility classes for the Worklist.
 */
// MakeDepend: library PTHREAD
// MakeDepend: library OMP
#include "system/System.h"
#include "system/WorklistUtils.h"
#include "system/Assert.h"
#include "system/ErrNo.h"
#include "system/SysConf.h"
#include "String.h"
#include <ctime>
#include <csignal>
#include <cstring>
#include <iostream>
#include <iterator>
#include <omp.h>
#include <sched.h>
#include <sstream>
#include <sys/resource.h>
#include <unistd.h>

namespace
{
    pthread_t gMainThread;
    CRD::HOOKFUNC gOldHook;

    void exitHook( int status )
    {
        if ( !pthread_equal(gMainThread,pthread_self()) )
            pthread_exit( reinterpret_cast<void*>(status) );
    }
}

WorklistParameterizer::WorklistParameterizer( size_t nUnits,
                                              size_t bytesPerUnit,
                                              size_t minBatchSize,
                                              unsigned nThreads )
{
    if ( nUnits < minBatchSize )
        minBatchSize = nUnits;

    if ( !nThreads || nThreads > processorsOnline() )
        nThreads = processorsOnline();

    size_t maxSimultaneous = .9*MemAvailable()/bytesPerUnit;
    if ( maxSimultaneous < minBatchSize )
        FatalErr("Not enough memory.");

    size_t maxThreads = maxSimultaneous/minBatchSize;
    if ( nThreads > maxThreads )
        nThreads = maxThreads;

    size_t maxBatchSize = maxSimultaneous/nThreads;
    if ( maxBatchSize > nUnits )
        maxBatchSize = nUnits;

    size_t minBatches = (nUnits+maxBatchSize-1)/maxBatchSize;
    size_t maxBatches = nUnits / minBatchSize;

    if ( maxBatches <= nThreads )
        nThreads = minBatches = maxBatches;
    else
    {
        // round up minBatches to an even divisor of nThreads, subject to
        // remaining less than maxBatches
        size_t nCycles = (minBatches+nThreads-1)/nThreads;
        size_t batches = nCycles * nThreads;
        if ( batches <= maxBatches )
            minBatches = batches;
        else
            minBatches = maxBatches;
    }
    mNBatches = minBatches;
    mNThreads = nThreads;
    mBatchSize = (nUnits+minBatches-1)/minBatches;

    ForceAssertLe(mNThreads*mBatchSize,maxSimultaneous);
}

QueueStateManipulator::~QueueStateManipulator()
{
    size_t size = mQS.getSize();
    bool done = mQS.isDone();

    if ( !mDone && done )
    {
        mQS.mCondNotEmpty.broadcast();
    }
    else if ( size > mSize )
    {
        if ( size-mSize == 1 )
            mQS.mCondNotEmpty.signal();
        else
            mQS.mCondNotEmpty.broadcast();
    }

    if ( mSize && !size )
    {
        mQS.mCondEmpty.broadcast();
    }
}

// populate_cpu_affinity - parse environment variable GOMP_CPU_AFFINITY and populate mGompCpuAffinity
//

static void populate_cpu_affinity( vector<int>& cpu_list )
{
    size_t nprocs = processorsOnline();
    const char *affp= getenv("GOMP_CPU_AFFINITY");
    if ( affp  ) {
	istringstream buff(affp);
	istream_iterator<std::string> item_begin(buff);
	istream_iterator<std::string> item_end;

	bool good=false;
	for ( auto item = item_begin; item != item_end; ++item ) {
	    size_t start = 0, stop = 0, stride = 1;
	    good = false;

	    istringstream sitem(*item);
	    if  (!(sitem >> start)) break;		// expect start
	    if ( sitem.good() ) {
		if ( sitem.peek() != '-' ) break;	// if more, must be dash (-)
		sitem.ignore();
		if (!(sitem >> stop)) break;		// expect stop
		if ( sitem.good() ) {
		    if ( sitem.peek() != ':' ) break;	// if more, must be colon (:)
		    sitem.ignore();
		    if (!(sitem >> stride)) break;	// expect stride
		}
	    }
	    good = true;

	    ForceAssertGt(nprocs,0U);
	    ForceAssertGt(stride,0U);
	    start = min( nprocs-1, start );		// clamp [0,nprocs)
	    stop = min( nprocs-1, stop );

	    size_t i = start;
	    do {
		cpu_list.push_back(i);
		i += stride;
	    } while ( i <= stop );	// test at end, so that if only start is set, we still get one value
	}

	if (!good) {
	    cpu_list.clear();
	    cout << "Warning: apparent nonsense in GOMP_CPU_AFFINITY: " << affp << endl;
	}
    }
}

// setThreadAffinity() -- assigns cpu affinities for threads round-robin (across subsequent calls to this function) based
// on the cpus listed in mGompCpuAffinity.
//
// Does nothing if mGompCpuAffinity is empty.
//
void ThreadPool::setThreadAffinity() {
    if ( mGompCpuAffinity.size() ) {
	size_t procno = mGompCpuAffinity[mNextAffinity];
	mNextAffinity = ( mNextAffinity + 1 ) % mGompCpuAffinity.size();
	cpu_set_t cpuset;
	CPU_ZERO(&cpuset);
	CPU_SET(procno, &cpuset);
	checkThreadOp(pthread_attr_setaffinity_np(&mAttr, sizeof(cpuset), &cpuset), "failed setting CPU affinity for threads: ");
    }
}


ThreadPool::ThreadPool( size_t nThreads, size_t threadStackSize,
                        void* (*threadFunc)(void*), void* threadFuncArg,
                        bool use_gomp_affinity )
: mThreadFunc(threadFunc), mThreadFuncArg(threadFuncArg),
  mNThreads(validateNThreads(nThreads)), mThreads(new pthread_t[mNThreads]),
  mNextAffinity(0), mCondExit(*this), mMonitorInterval(0),
  mReportRestarts(false)
{
    omp_set_num_threads(1);
    checkThreadOp(pthread_attr_init(&mAttr),
                  "ThreadAttr initialization failed: ");

    if ( use_gomp_affinity ) populate_cpu_affinity( mGompCpuAffinity );

    if ( !threadStackSize )
    {
        struct rlimit rl;
        if ( getrlimit(RLIMIT_STACK,&rl) )
        {
            ErrNo err;
            FatalErr("Unable to determine current stack size" << err);
        }

        if ( rl.rlim_cur == RLIM_INFINITY )
            threadStackSize = DEFAULT_THREAD_STACK_SIZE;
        else
            threadStackSize = rl.rlim_cur;
    }
    checkThreadOp(pthread_attr_setstacksize(&mAttr,threadStackSize),
                  "Unable to set stack size for threads: ");

    setThreadAffinity();	// one for the main thread
    pthread_t* pThread = mThreads + mNThreads;
    while ( pThread-- > mThreads )
    {
	setThreadAffinity();
        checkThreadOp(pthread_create(pThread,&mAttr,threadFunc,threadFuncArg),
                      "Thread creation failed: ");
    }
}

ThreadPool::~ThreadPool()
{
    shutdown();
    checkThreadOp(pthread_attr_destroy(&mAttr),
                    "ThreadAttr destruction failed: ");
    omp_set_num_threads(getConfiguredNumThreads());
}

void ThreadPool::threadDeathMonitor( long newInterval,
                                     bool reportRestarts )
{
    if ( mNThreads )
    {
        long oldInterval;
        if ( true )
        {
            Locker ml(*this);

            oldInterval = mMonitorInterval;
            mMonitorInterval = newInterval;
            mReportRestarts = reportRestarts;
        }

        if ( !oldInterval && newInterval )
        {
            gMainThread = pthread_self();
            gOldHook = CRD::installExitHook(&exitHook);
            checkThreadOp(pthread_create(&mMonitorThread,0,threadFunc,this),
                          "Can't create ThreadPool monitor: ");
        }
        else if ( oldInterval && !newInterval )
        {
            mCondExit.signal();
            checkThreadOp(pthread_join(mMonitorThread,0),
                          "ThreadPool monitor join failed: ");
            CRD::installExitHook(gOldHook);
        }
    }
}

size_t ThreadPool::findThreadIndex( pthread_t handle )
{
    Locker ml(*this);
    pthread_t* pThread = mThreads + mNThreads;
    while ( pThread-- > mThreads )
    {
        if ( pthread_equal(handle,*pThread) )
            return pThread - mThreads;
    }
    FatalErr("Unable to find thread's index.");
}

void ThreadPool::shutdown()
{
    threadDeathMonitor(0);

    bool somebodyDied = false;
    pthread_t* pThread = mThreads + mNThreads;
    while ( pThread-- > mThreads )
    {
        void* retVal = 0;
        checkThreadOp(pthread_join(*pThread,&retVal),"Thread join failed: ");
        if ( reinterpret_cast<long>(retVal) )
        {
            std::cout << "Thread #" << (pThread-mThreads)
                        << " died with a return value of "
                        << reinterpret_cast<long>(retVal) << std::endl;
            somebodyDied = true;
        }
    }
    if ( somebodyDied )
        FatalErr("Quitting due to failure of child threads.");

    delete [] mThreads;
    mThreads = 0;
    mNThreads = 0;
}

void ThreadPool::monitorThreads()
{
    Locker ml(*this);
    while ( true )
    {
        ml.timedWait(mCondExit,mMonitorInterval);
        if ( !mMonitorInterval )
            break;

        pthread_t* pThread = mThreads + mNThreads;
        setNextAffinity(0);	// restart the list
        setThreadAffinity();	// one for the main thread
        while ( pThread-- > mThreads )
        {
	    setThreadAffinity();	// jump past each element just in case we re-create the thread
            if ( pthread_kill(*pThread,0) ==  ESRCH )
            {
                void* retVal;
                checkThreadOp(pthread_join(*pThread,&retVal),
                              "Can't join to dead thread: ");
                if ( mReportRestarts )
                    std::cout << "Restarting dead thread #" <<
                                 (pThread-mThreads) <<
                                 " which died with a return value of "
                                 << reinterpret_cast<long>(retVal) << std::endl;
                checkThreadOp(pthread_create(pThread,&mAttr,mThreadFunc,mThreadFuncArg),
                              "Can't create replacement thread: ");
            }
        }
    }
}

size_t ThreadPool::validateNThreads( size_t nThreads )
{
    size_t nProcs = processorsOnline();
    if ( nThreads > 3*nProcs )
        FatalErr("Trying to start a thread pool with " << nThreads <<
                 " threads, but there are only " << nProcs <<
                 " processors, which probably isn't going to work out well.");

    return nThreads;
}

void* ThreadPool::threadFunc( void* ptr )
{
    reinterpret_cast<ThreadPool*>(ptr)->monitorThreads();
    return 0;
}


void ThreadPool::die( int errNo, char const* msg )
{
    ErrNo err(errNo);
    FatalErr(msg << err);
}


/* -------------------------------------------------------------------------------
 * NaiveThreadPool
 *
 *  Associates an idle thread_ID with an item_ID.
 *
 * 2009-12 ribeiro@broadinstitute.org
 */

NaiveThreadPool::NaiveThreadPool(const int num_threads,
                                 const int num_items) :
  num_threads_(num_threads),
  num_items_(num_items)
{
  item_IDs_ = new int[num_threads];
  for (int i = 0; i < num_threads; i++)
    item_IDs_[i] = -1;

  thread_IDs_ = new int[num_items];
  for (int i = 0; i < num_items; i++)
    thread_IDs_[i] = -1;
}

NaiveThreadPool::~NaiveThreadPool()
{
  delete [] item_IDs_;
  delete [] thread_IDs_;
}

// Obtain a thread_ID and associate it with an item_ID
//    item_IDs_[thread_ID] gives the item_ID that thread_ID is processing
//    thread_IDs_[item_ID] gives the thread_ID that processed the item_ID
int NaiveThreadPool::AssignThread(const int item_ID)
{
  Locker lock(*this);
  int thread_ID = 0;

  // search list of thread_IDs for one that is available
  // (i.e., item_IDs_[thread_ID_] == -1)
  while (thread_ID < num_threads_ && item_IDs_[thread_ID] >= 0)
    thread_ID++;

  // if this assert fails something is wrong.
  // if we ask for a thread then there should be at least one available.
  ForceAssertLt(thread_ID, num_threads_);

  item_IDs_[thread_ID] = item_ID;
  thread_IDs_[item_ID] = thread_ID;

  return thread_ID;
}


void NaiveThreadPool::UnassignThread(const int item_ID)
{
  const int & thread_ID = thread_IDs_[item_ID];
  if (thread_ID != -1)
    item_IDs_[thread_ID] = -1;
  else
    std::cout << "trying to dissociate unassociated item_ID." << std::endl;
}








/* -----------------------------------------------------------------------------
 * ThreadBlockLocks
 *
 *  Multithreaded modules sometimes need to write to a shared data structure.
 *  You can mutex_lock everything down to write, but then the rest of the
 *  threads have to wait. The solution is to divide the data structure in blocks
 *  and lock only the block that is currently being written to, letting the
 *  other threads do something else. Also, data should be buffered so that you
 *  are not locking and unlocking all the time.
 *
 * 2009-12 ribeiro@broadinstitute.org
 */
ThreadBlockLocks::ThreadBlockLocks(const int num_blocks,
                                   const int num_threads) :
  num_blocks_(num_blocks),
  num_threads_(num_threads)
{
  thread_IDs_ = new int[num_blocks];
  for (int i = 0; i < num_blocks; i++)
    thread_IDs_[i] = -1;
}

ThreadBlockLocks::~ThreadBlockLocks()
{
  delete [] thread_IDs_;
}

// a debug function
void ThreadBlockLocks::OutputBlocked(const int block_ID, const int thread_ID)
{
  std::cout << "blocker("<< block_ID <<"): ";
  for (int i = 0; i != num_blocks_; i++) {
    const int & tID = thread_IDs_[i];
    if (tID < 0) std::cout << " -" << std::endl;
    else         std::cout << " " << tID << std::endl;
  }
}


// Locks block_ID, updating thread_IDs vec with the blocking thread_ID.
// returns the number of tries.
int ThreadBlockLocks::LockBlock(const int block_ID, const int thread_ID)
{
  int tries = 0;
  while (thread_IDs_[block_ID] != thread_ID) {
    Locker lock(*this);
    if (thread_IDs_[block_ID] < 0) // test for availability
      thread_IDs_[block_ID] = thread_ID;
    tries++;
  }
  return tries;
}

// Locks one block out of several possible in a set,
// updating thread_IDs vec with the blocking thread_ID.
// returns the block_ID.
int ThreadBlockLocks::LockSomeBlock(const std::set<int> & block_IDs,
                                    const int thread_ID)
{
  bool found = false;

  std::set<int>::iterator i_block_ID = block_IDs.begin();
  while (!found) {

    if (thread_IDs_[*i_block_ID] < 0) {   // block is not locked
      Locker lock(*this);
      if (thread_IDs_[*i_block_ID] < 0) { // test for availability again
        thread_IDs_[*i_block_ID] = thread_ID;
        found = true;
      }
    }

    if (!found) {
      i_block_ID++;
      if (i_block_ID == block_IDs.end())
        i_block_ID = block_IDs.begin();
    }
  }
  return *i_block_ID;
}


// Unlocks block_ID.  Should only be done by the thread that locked the block.
// If so, there is no need to use mutex_locks.
void ThreadBlockLocks::UnlockBlock(const int block_ID, const int thread_ID)
{
  ForceAssertEq(thread_IDs_[block_ID], thread_ID);
  thread_IDs_[block_ID] = -1;
}
