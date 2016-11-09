///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file FileWriter.h
 * \author tsharpe
 * \date Jan 20, 2012
 *
 * \brief Handle writing to a file descriptor.
 */
#ifndef FILEWRITER_H_
#define FILEWRITER_H_

#include "system/file/FileReader.h"

#include <unistd.h>

class FileWriter : public FileReader
{
public:
    FileWriter() {}

    explicit FileWriter( char const* path, bool atEnd=false )
    : FileReader(doOpen(path,atEnd,false),path)
    { mMyFD = true; if ( atEnd ) seekEnd(); }


    /// Construct from something with a c_str() member (like string or String)
    template <class C>
    explicit FileWriter( C const& path, bool atEnd=false,
                            bool rdOnly = false,
                            char const* (C::*)() const = &C::c_str )
    : FileReader(doOpen(path.c_str(),atEnd,rdOnly),path.c_str())
    { mMyFD = true; if ( atEnd ) seekEnd(); }

    /// NB: This is for situations where the fd isn't in the filesystem (e.g.,
    /// pipes, sockets, etc.).  You still own the fd, and it will NOT be
    /// automatically closed for you.
    FileWriter( int fd, char const* pseudoFilename )
    : FileReader(fd,pseudoFilename) {}

    // copying allowed if (or when) base class allows it
    // compiler-suppied destructor is OK

    FileWriter const& write( void const* buf, size_t len ) const;

    void flush() { if ( mFD != -1 ) fsync(mFD); }

private:
    // the rdOnly silliness is necessary for LookupTable, which sometimes
    // wants to treat a FileWriter solely as a FileReader.
    static int doOpen( char const* path, bool noTrunc, bool rdOnly );
};

#endif /* FILEWRITER_H_ */
