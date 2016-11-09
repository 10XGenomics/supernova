///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file SysConf.cc
 * \author tsharpe
 * \date Apr 15, 2010
 *
 * \brief
 */
// MakeDepend: library OMP
#include "system/SysConf.h"
#include "system/Exit.h"
#include <climits>
#include <cstdlib>
#include <iostream>
#include <unistd.h>
#include <malloc.h>
#include <omp.h>

namespace
{

void sysconfErr( char const* attr )
{
    std::cout << "Sysconf unable to determine value of " << attr << std::endl;
    CRD::exit(1);
}

size_t gPagSiz;
size_t gPhysPages;
size_t gProcsOnline;
size_t gClockTicksPerSec;
size_t gMaxHostNameLen;
int gNumThreads;

}

size_t pageSize()
{
    if ( !gPagSiz )
    {
        long result = sysconf(_SC_PAGESIZE);
        if ( result == -1 ) sysconfErr("_SC_PAGESIZE");
        gPagSiz = result;
    }
    return gPagSiz;
}

size_t physicalMemory()
{
    if ( !gPhysPages )
    {
        long result = sysconf(_SC_PHYS_PAGES);
        if ( result == -1 ) sysconfErr("_SC_PHYS_PAGES");
        gPhysPages = result;
    }
    return gPhysPages*pageSize();
}

size_t processorsOnline()
{
    if ( !gProcsOnline )
    {
        long result = sysconf(_SC_NPROCESSORS_ONLN);
        if ( result == -1 ) sysconfErr("_SC_NPROCESSORS_ONLN");
        if ( result < 1 || result > INT_MAX )
        {
            std::cout << "Sysconf's value for the number of processors online ("
                      << result << ") doesn't make sense." << std::endl;
            CRD::exit(1);
        }
        gProcsOnline = result;
    }
    return gProcsOnline;
}

// If the argument is silly (meaning not positive), then we check the
// environment variable OMP_THREAD_LIMIT for inspiration.  if that's not set
// we go with the number of processors online.  we do a final check to make
// sure that, wherever we got the number, that it isn't set to more than the
// number of processors online.
int configNumThreads( int numThreads )
{
    int procsOnline = processorsOnline();
    if ( numThreads < 1 )
    {
        char const* ompEnv = getenv("OMP_THREAD_LIMIT");
        if ( !ompEnv )
            numThreads = procsOnline;
        else
        {
            char* end;
            long ompVal = strtol(ompEnv,&end,10);
            if ( *end || ompVal < 1 || ompVal > INT_MAX )
            {
                std::cout << "Environment variable OMP_THREAD_LIMIT ("
                          << ompEnv
                          << " cannot be parsed as a positive integer."
                          << std::endl;
                CRD::exit(1);
            }
            numThreads = ompVal;
        }
    }
    if ( numThreads > procsOnline )
        numThreads = procsOnline;

    if ( numThreads > 32 )
        numThreads = 32;           // minimal kludge for 1.0

    gNumThreads = numThreads;

#if defined(M_ARENA_MAX)
    // Set M_ARENA_MAX to be a bit more than eight times the number of threads.
    // This value is based on what glibc 2.12.2 would do, eventually, during a run.
    // This sets it up-front.  The seems to eliminate much of the variability between
    // runs probably by eliminating glibc's "adaptive" behavior in malloc, but that
    // has not been demonstrated conclusively.
    size_t arena_max = ((numThreads) * (sizeof(long) == 4 ? 2 : 8)) + 5;
    if ( ! mallopt( M_ARENA_MAX, arena_max ) )
	std::cout << "Warning: failed trying to set mallopt M_ARENA_MAX=" << arena_max << std::endl;
#endif

    return numThreads;
}

int getConfiguredNumThreads()
{
    if ( !gNumThreads )
        gNumThreads = processorsOnline();
    return gNumThreads;
}

int boundNumThreads( int numThreads )
{ if ( inParallelSection() ) numThreads = 1;
  else if ( numThreads < 1 ) numThreads = getConfiguredNumThreads();
  else numThreads = std::min(numThreads,getConfiguredNumThreads());
  return numThreads; }

size_t clockTicksPerSecond()
{
    if ( !gClockTicksPerSec )
    {
        long result = sysconf(_SC_CLK_TCK);
        if ( result == -1 ) sysconfErr("_SC_CLK_TCK");
        gClockTicksPerSec = result;
    }
    return gClockTicksPerSec;
}

size_t maxHostNameLen()
{
    if ( !gMaxHostNameLen )
    {
        long result = sysconf(_SC_HOST_NAME_MAX);
        if ( result == -1 ) sysconfErr("_SC_HOST_NAME_MAX");
        gMaxHostNameLen = result;
    }
    return gMaxHostNameLen;
}
