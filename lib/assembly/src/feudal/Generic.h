///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////


#ifndef GENERIC_H_
#define GENERIC_H_

#include "system/ThreadsafeIO.h"
#include <cstddef>
#include <fstream>

class UtilizationReporter
{
public:
    void report( void const* addr, size_t bits, size_t siz, size_t cap,
                    char const* type )
    { if ( mpOS ) writeLine(addr,bits,siz,cap,type); }

    static UtilizationReporter gInstance;

private:
    UtilizationReporter();
    ~UtilizationReporter();
    UtilizationReporter( UtilizationReporter const& )=delete;
    UtilizationReporter& operator=( UtilizationReporter const& )=delete;

    void writeLine( void const* addr, size_t bits, size_t siz, size_t cap,
                    char const* type );

    std::ofstream mLog;
    ThreadsafeOStream* mpOS;
};

#endif // GENERIC_H_
