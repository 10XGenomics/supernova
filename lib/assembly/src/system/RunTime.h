///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef RUNTIME
#define RUNTIME

#include <ostream>
#include <signal.h>
#include "system/Exit.h"
#include "system/Types.h"

typedef void ArachneSignalHandler(int, siginfo_t*, void*);

void arachne_signal_handler( int signal_number, siginfo_t* info, void* context,
    Bool no_ctrlc = False );
void arachne_signal_handler_standard( int signal_number, siginfo_t* info,
    void* context );
void arachne_signal_handler_no_ctrlc_traceback( int signal_number,
    siginfo_t* info, void* context );
void ArachneInterruptHandler(ArachneSignalHandler* func);

void NoDump( ); // Turn off core dumps.

void RunTime( int no_dump=1,
              ArachneSignalHandler* pSigFunc=&arachne_signal_handler_standard );

// RunTimeLax is now just a synonym for RunTime
inline void RunTimeLax( )
{    RunTime();    }

inline void RunTimeNoTraceback( )
{    RunTime( 1, &arachne_signal_handler_no_ctrlc_traceback );    }

void print_signal_handler_message( int signal_number, siginfo_t* info );

void TracebackThisProcess( ostream& out, Bool exit_when_done, Bool minimal );

#endif
