///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
#include "system/ProcBuf.h"
#include "system/ErrNo.h"
#include "system/System.h"
#include <cerrno>
#include <cstddef>
#include <string>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

namespace
{
size_t const BUFFER_SIZE = 8192;

std::string getName( char const* cmd, std::ios_base::openmode mode )
{
    return std::string("pipe ") +
            ((mode & std::ios_base::in) ? "from" : "to") +
            " the command '" + cmd + "'";
}
}


procbuf::procbuf( const char *command, std::ios_base::openmode mode,
    bool expect_ret_zero ) : mFD(doOpen(command,mode)), mResult(0),
    mPID(0), mFW(mFD,getName(command,mode).c_str()), mCMD(command),
    mExpectRetZero(expect_ret_zero) {}

procbuf::int_type procbuf::overflow( int_type c )
{
    if ( traits_type::eq_int_type(c,traits_type::eof()) )
        return c;

    if ( pbase() == epptr() )
    {
        char* putBuf = new char[BUFFER_SIZE];
        setp(putBuf,putBuf+BUFFER_SIZE);
    }
    else if ( !flush() )
        return traits_type::eof();

    return sputc(c);
}


procbuf::int_type procbuf::underflow()
{
    if ( gptr() < egptr() || fill() )
        return traits_type::to_int_type( *gptr() );
    return traits_type::eof();
}


procbuf::int_type procbuf::pbackfail( int_type c )
{
    // If the character is EOF, then the caller isn't telling us what it is.
    // If we lost the buffer... return failure
    bool bKnowsChar = !traits_type::eq_int_type(c,traits_type::eof());

    // Check for room at the beginning of the buffer
    if ( gptr() == eback() || (!bKnowsChar && gptr() == eback()+1) )
        return traits_type::eof();

    gbump(-1);
    if ( bKnowsChar )
        *gptr() = traits_type::to_char_type(c);
    return traits_type::not_eof(c);
}

procbuf::int_type procbuf::sync()
{
    if ( pptr() != pbase() && !flush() )
        return traits_type::eof();
    return traits_type::not_eof(0);
}

bool procbuf::flush()
{
    if ( pptr() != pbase() )
        mFW.write(pbase(),pptr()-pbase());
    setp(pbase(),epptr());
    return true;
}

bool procbuf::fill()
{
    if ( eback() == egptr() )
    {
        char* putBuf = new char[BUFFER_SIZE+1];
        char* putBufEnd = putBuf + BUFFER_SIZE + 1;
	setg(putBuf,putBufEnd,putBufEnd);
    }
    char* buf = eback()+1;
    size_t len = mFW.readOnce(buf,BUFFER_SIZE);
    setg(eback(),buf,buf+len);
    return len;
}

int procbuf::doOpen( char const* command, std::ios_base::openmode mode )
{
    int pipe_fds[2];
    if ( pipe(pipe_fds) )
    {
        ErrNo err;
        FatalErr("Opening " << getName(command,mode) << " failed" << err);
    }

    pid_t pid = fork();
    if ( !pid )
    {
        if ( (mode & ios_base::in) )
        {
            while ( dup2(pipe_fds[1], 1) == -1 )
            {
                ErrNo err;
                if ( err.val() != EINTR )
                    FatalErr("Can't dup child's pipe fd to stdout.");
            }
        }
        else
        {
            while ( dup2(pipe_fds[0], 0) == -1 )
            {
                ErrNo err;
                if ( err.val() != EINTR )
                    FatalErr("Can't dup child's pipe fd to stdin.");
            }
        }
        while ( ::close(pipe_fds[0]) == -1 )
        {
            ErrNo err;
            if ( err.val() != EINTR )
                FatalErr("Can't close parent's pipe fd 0 in child.");
        }
        while ( ::close(pipe_fds[1]) == -1 )
        {
            ErrNo err;
            if ( err.val() != EINTR )
                FatalErr("Can't close parent's pipe fd 1 in child.");
        }

        execl("/bin/sh", "sh", "-c", command, (char *) 0);
        _exit(-1);
    }

    if ( pid == -1 )
    {
        ErrNo err;
        FatalErr("Forking child to make a pipeline failed" << err);
    }

    mPID = pid;

    int result;
    if ( (mode & ios_base::in) )
    {
        result = pipe_fds[0];
        while ( ::close(pipe_fds[1]) == -1 )
        {
            ErrNo err;
            if ( err.val() != EINTR )
                FatalErr("Can't close child's pipe fd 1 in parent.");
        }
    }
    else
    {
        result = pipe_fds[1];
        while ( ::close(pipe_fds[0]) == -1 )
        {
            ErrNo err;
            if ( err.val() != EINTR )
                FatalErr("Can't close parent's pipe fd 0 in child.");
        }
    }

    return result;
}

int procbuf::doClose()
{
    flush();
    mFW.close();
    while ( ::close(mFD) == -1 )
    {
        ErrNo err;
        if ( err.val() != EINTR )
            FatalErr("Unable to close pipe" << err);
    }
    mFD = -1;

    /* While there is no child complete, or an interrupt, wait: */
    int wstatus;
    while ( waitpid(mPID,&wstatus,0) == -1 )
    {
        ErrNo err;
        if ( err.val() != EINTR )
            FatalErr("Unable to wait for procbuf's child" << err);
    }
    
    if ( WIFEXITED(wstatus) )
    {
        if ( mExpectRetZero && WEXITSTATUS(wstatus) != 0 )
        {
            FatalErr(mCMD << " [" << mPID << "] exited with "
                << WEXITSTATUS(wstatus) << "\n");
        }
    }
    else
    {
        FatalErr(mCMD << " [" << mPID << "] exited abnormally\n");
    }
    return wstatus;
}
