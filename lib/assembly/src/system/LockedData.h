///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * \file LockedData.h
 * \author tsharpe
 * \date Jul 17, 2009
 *
 * \brief
 */

#ifndef SYSTEM_LOCKEDDATA_H_
#define SYSTEM_LOCKEDDATA_H_

#include "system/Assert.h"
#include <mutex>
#include <condition_variable>

typedef std::mutex LockedData;

class Condition : public std::condition_variable
{
public:
    Condition( LockedData& lockedData ) : mMutex(lockedData) {}

    void signal() { notify_one(); }
    void broadcast() { notify_all(); }

    LockedData const* mutex() const { return &mMutex; }

private:
    LockedData& mMutex;
};

class Locker
{
public:
    Locker( LockedData& lockedData ) : mLock(lockedData) {}
    Locker( Locker const& )=delete;

    void wait( Condition& cv )
    { AssertEq(mLock.mutex(),cv.mutex()); cv.wait(mLock); }

    bool timedWait( Condition& cv, long nSecs )
    { AssertEq(mLock.mutex(),cv.mutex());
      std::chrono::seconds duration(nSecs);
      return cv.wait_for(mLock,duration) == std::cv_status::timeout; }

private:
    std::unique_lock<LockedData> mLock;
};

#endif /* SYSTEM_LOCKEDDATA_H_ */
