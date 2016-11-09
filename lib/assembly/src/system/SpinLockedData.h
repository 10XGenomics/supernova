///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file SpinLockedData.h
 * \author tsharpe
 * \date Nov 16, 2011
 *
 * \brief
 */
#ifndef SYSTEM_SPINLOCKEDDATA_H_
#define SYSTEM_SPINLOCKEDDATA_H_

#include "system/Assert.h"
#include <atomic>
#include <iostream>

#ifdef NDEBUG
/// A spin-lock.
class SpinLockedData
{
public:
    SpinLockedData() : mLockByte(false) {}
    SpinLockedData( SpinLockedData const& )=delete;
    SpinLockedData& operator=( SpinLockedData const& )=delete;
    void lock() { while ( mLockByte.exchange(true) ) {} }
    void unlock() { mLockByte = false; }

private:
    std::atomic_bool mLockByte;
};
#else
/// A spin-lock.
class SpinLockedData
{
public:
    SpinLockedData() : mLockByte(false), mBusyCount(0) {}
    SpinLockedData( SpinLockedData const& )=delete;
    SpinLockedData& operator=( SpinLockedData const& )=delete;

    ~SpinLockedData()
    { AssertEq(bool(mLockByte),false);
      if ( mBusyCount >= 1024 )
        std::cout << "Busy spin lock: " << mBusyCount << std::endl; }

    void lock()
    { size_t count = 0;
      while ( mLockByte.exchange(true) )
      { ++count; }
      mBusyCount += count >> 10; }

    void unlock() { AssertEq(bool(mLockByte),true); mLockByte = false; }

private:
    std::atomic_bool mLockByte;
    unsigned mBusyCount;
};
#endif

/// Something that operates a spin-lock, and never forgets to unlock it.
class SpinLocker
{
public:
    SpinLocker( SpinLockedData& lock ) : mLock(lock) { mLock.lock(); }
    ~SpinLocker() { mLock.unlock(); }

private:
    SpinLocker( SpinLocker const& ); // unimplemented -- no copying
    SpinLocker& operator=( SpinLocker const& ); // unimplemented -- no copying

    SpinLockedData& mLock;
};

#endif /* SYSTEM_SPINLOCKEDDATA_H_ */
