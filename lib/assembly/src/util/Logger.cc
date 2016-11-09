/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////
/*
 * \file Logger.cc
 * \author tsharpe
 * \date Feb 20, 2009
 *
 * \brief A class for logging error messages.
 */
#include <cstdio>
#include <cstring>
#include <cstdarg>
#include "util/Logger.h"
#include "util/NullOStream.h"

using std::ostream;
using std::streambuf;
using std::string;
using std::map;
using std::endl;

Logger Logger::nullLogger()
{
    static NullOStream gNullOS;
    return Logger(gNullOS);
}

void Logger::log( char const* str ... )
{
    MsgPacket& pkt = mMap[str];
    if ( pkt.getCount() <= mMaxInstances )
    {
        va_list ap;
        va_start(ap,str);
        char buf[BUF_SIZE];
        vsnprintf(buf,BUF_SIZE,str,ap);
        if ( pkt.getCount() == mMaxInstances )
            pkt.setMessage(buf);
        else
            mOS << buf << endl;
        va_end(ap);
    }
    pkt.incrementCount();
}

void Logger::flush()
{
    map<char const*,MsgPacket>::iterator end = mMap.end();
    for ( map<char const*,MsgPacket>::iterator itr = mMap.begin(); itr != end; ++itr )
    {
        itr->second.flush(mOS,mMaxInstances);
    }
    mMap.clear();
}

void Logger::MsgPacket::flush( ostream& os, size_t maxInstances ) const
{
    if ( mOverflow )
    {
        size_t nnn = mCount - maxInstances;
        if ( nnn == 1 )
            os << mMessage << endl;
        else
            os << "\nThere were an additional " << (mCount-maxInstances) << " messages similar to this one:\n" << mMessage << endl;
    }
}
