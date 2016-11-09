///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * \file Exit.h
 * \author tsharpe
 * \date Feb 6, 2009
 *
 * \brief A replacement for the ::exit(int) function.
 *
 * Calls ::exit(0) if the arg is zero, calls abort() otherwise.
 * You can change this behavior by installing a hook, which will be
 * called instead.
 */
#ifndef SYSTEM_EXIT_H_
#define SYSTEM_EXIT_H_

namespace CRD
{

void exit( int ) __attribute__((__noreturn__));

typedef void(*HOOKFUNC)(int);
HOOKFUNC installExitHook( HOOKFUNC fHook );

} // end namespace CRD

inline void TracebackThisProcess() { CRD::exit(1); }

#endif /* SYSTEM_EXIT_H_ */
