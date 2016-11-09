///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file TempFile.h
 * \author tsharpe
 * \date Feb 22, 2012
 *
 * \brief
 */
#ifndef TEMPFILE_H_
#define TEMPFILE_H_

#include "String.h"
#include <cstdlib>

class temp_file
{
public:
    /// Creates a temporary file, based on a name template that you supply.
    explicit temp_file( char const* path )
    : mName(generateName(path)) {}

    /// Construct from something with a c_str() member (like string or String)
    template <class C>
    explicit temp_file( C const& path,
                            char const*(C::*)() const = &C::c_str )
    : mName(generateName(path.c_str())) {}

    /// Removes the temporary file.
    ~temp_file();

    operator String const&() const { return mName; }
    char const* c_str() const { return mName.c_str(); }

    /// When you don't want the automatic deletion feature.
    /// This creates the temporary file and returns its name.
    static String generateName( char const* path );

    /// Returns the environment variable TMPDIR, or "/tmp" if that's not set
    static String tmpDir()
    { char const* tmpDir = getenv("TMPDIR");
      if ( !tmpDir ) tmpDir = "/tmp";
      return tmpDir; }

private:
    String mName;
};

#endif /* TEMPFILE_H_ */
