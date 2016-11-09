///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * Thread.cc
 *
 *  Created on: Jul 31, 2014
 *      Author: tsharpe
 */
// MakeDepend: library OMP
#include "system/Thread.h"
#include <sys/syscall.h>
#include <unistd.h>
#include <omp.h>

int getTid()
{
    return syscall(SYS_gettid);
}

bool isMainThread()
{
    static int gPid = getpid();
    return getTid()==gPid;
}

bool inParallelSection()
{
    return !isMainThread() || omp_in_parallel();
}
