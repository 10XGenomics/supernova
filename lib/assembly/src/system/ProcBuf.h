///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef PROCBUF_H
#define PROCBUF_H

#include "system/file/FileWriter.h"
#include <streambuf>

class procbuf : public std::basic_streambuf<char>
{
  public:
    procbuf( char const* command, std::ios_base::openmode mode,
        bool expect_ret_zero = false );
    procbuf( procbuf const & ) = delete;
    procbuf& operator=( procbuf const& ) = delete;
    ~procbuf() { close(); delete [] eback(); delete [] pbase(); }

    bool is_open() { return mFD != -1; }
    int close()
    { if ( is_open() ) mResult = doClose();
      return mResult; }

    // read as many characters as are immediately available from the pipe
    size_t read( void* buf, size_t len ) { return mFW.readOnce(buf,len); }

  protected:
    int_type overflow( int_type c = traits_type::eof() );
    int_type underflow();
    int_type pbackfail( int_type c );
    int_type sync();

  private:
    bool flush();
    bool fill();
    int doOpen( char const* command, std::ios_base::openmode mode );
    int doClose();

    int mFD;
    int mResult;
    pid_t mPID;
    FileWriter mFW;
    std::string mCMD;
    bool mExpectRetZero;
};

#endif
