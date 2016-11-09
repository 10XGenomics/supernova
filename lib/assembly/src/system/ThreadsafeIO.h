///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * \file ThreadsafeIO.h
 * \author tsharpe
 * \date May 14, 2009
 *
 * \brief Provides threadsafe wrappers to control access to an ostream by multiple threads.
 *
 * To use the static instances that wrap cout and cerr, replace cout in your code with ThreadsafeOStreamFetcher::cout()
 * and cerr with ThreadsafeOStreamFetcher::cerr().  You can define a couple of macros to make this easier,
 * if you like, e.g., #define cout ThreadsafeOStreamFetcher::cout()
 */
#ifndef SYSTEM_THREADSAFEIO_H_
#define SYSTEM_THREADSAFEIO_H_

#include "system/SpinLockedData.h"
#include <ostream>

/// Streambuf which locks a mutex and then writes to a wrapped ostream on overflow and sync.
class ThreadsafeStreambuf : public std::streambuf
{
public:
    ThreadsafeStreambuf( std::ostream& os )
    : mOS(os)
    { setp(mBuf,mBuf+sizeof(mBuf)-1); }

    ThreadsafeStreambuf( ThreadsafeStreambuf const& )=delete;

    ~ThreadsafeStreambuf()
    { sync(); }

private:
    int_type overflow( int_type ch );
    int sync();

    std::ostream& mOS;
    char mBuf[4096];
    static SpinLockedData gLock;
};

/// An ostream that wraps another.
/// If each thread uses a separate one of these, then access to the wrapped ostream becomes threadsafe.
class ThreadsafeOStream : public std::ostream
{
public:
    ThreadsafeOStream( std::ostream& os )
    : std::ostream(&mSB), mSB(os)
    {}

    ThreadsafeOStream( ThreadsafeOStream const& )=delete;

private:
    ThreadsafeStreambuf mSB;
};

#endif /* SYSTEM_THREADSAFEIO_H_ */
