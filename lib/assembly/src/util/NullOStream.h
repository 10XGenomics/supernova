///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file NullOStream.h
 * \author tsharpe
 * \date Nov 7, 2012
 *
 * \brief
 */
#ifndef NULLOSTREAM_H_
#define NULLOSTREAM_H_

#include <ostream>

class NullOStream : private std::streambuf, public std::ostream
{
public:
    NullOStream() : std::ostream(this) {}

    // compiler-supplied copying and destructor are OK

protected:
    virtual int overflow( int c )
    { setp(mBuf,mBuf+sizeof(mBuf)); return 0; }

private:
    char mBuf[1024];
};

#endif /* NULLOSTREAM_H_ */
