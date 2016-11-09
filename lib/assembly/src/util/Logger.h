/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////
/*
 * \file Logger.h
 * \author tsharpe
 * \date Feb 20, 2009
 *
 * \brief A class for logging error messages.
 *
 * Meant to avoid the tedium of seeing a zillion similar error messages.
 * It saves up a count of the number of instances of various message types seen,
 * and will simply report the number of instances (with the first few examples)
 * when the number of instances is large.  If there are only a few instances,
 * they're all reported when the Logger is flushed.
 */


#ifndef UTIL_LOGGER_H_
#define UTIL_LOGGER_H_

#include <iostream>
#include <string>
#include <map>

class Logger
{
public:
    /// log to the specified stream
    explicit Logger( std::ostream& os,
                        size_t maxInstances = DFLT_MAX_INSTANCES )
    : mOS(os), mMaxInstances(maxInstances)
    {}

    // compiler-supplied copy constructor is OK
    // assignment won't work because of reference

    /// flushes pending error messages on its way to oblivion
    ~Logger()
    { flush(); }

    /// log an error message.  sprintf semantics.
    /// the format string really, really has to be a static string.
    /// we're using its address to figure out classes of error messages.
    void log( char const* str ... );

    /// flush pending error messages
    void flush();

    class MsgPacket
    {
    public:
        MsgPacket()
        : mCount(0), mOverflow(false)
        {}

        // compiler-supplied copying and destructor are OK

        size_t getCount() const
        { return mCount; }

        size_t incrementCount()
        { return mCount += 1; }

        void setMessage( std::string const& msg )
        { mMessage = msg; mOverflow = true; }

        void flush( std::ostream& os, size_t maxInstances ) const;

    private:
        std::string mMessage; // first overflow message
        size_t mCount;
        bool mOverflow;
    };

    static Logger nullLogger();

private:
    std::ostream& mOS;
    size_t mMaxInstances;
    std::map<char const*,MsgPacket> mMap;

    static size_t const DFLT_MAX_INSTANCES = 3;
    static size_t const BUF_SIZE = 8192;
};

#endif /* UTIL_LOGGER_H_ */
