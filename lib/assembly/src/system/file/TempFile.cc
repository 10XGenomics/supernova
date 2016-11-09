///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file TempFile.cc
 * \author tsharpe
 * \date Feb 22, 2012
 *
 * \brief
 */

#include "system/file/TempFile.h"
#include "system/ErrNo.h"
#include "system/System.h"
#include <algorithm>
#include <cerrno>
#include <cstdlib>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

temp_file::~temp_file()
{
    if ( unlink(mName.c_str()) == -1 )
    {
        ErrNo err;
        FatalErr("Can't unlink temporary file " << mName << err);
    }
}

String temp_file::generateName( char const* path )
{
    String result(path);
    if ( result.size() < 6 ||
            !std::equal(result.end()-6,result.end(),"XXXXXX") )
        result += "XXXXXX";
    int fd = mkstemp(const_cast<char*>(result.c_str()));
    if ( fd == -1 )
    {
        ErrNo err;
        FatalErr("Failed to create temporary file from " << result << err );
    }
    while ( close(fd) == -1 )
    {
        ErrNo err;
        if ( err.val() != EINTR )
            FatalErr("Unable to close temporary file " << result << err);
    }
    if ( chmod(result.c_str(),0664) == -1 )
    {
        ErrNo err;
        FatalErr("Failed to chmod temporary file " << result << err);
    }
    return result;
}
