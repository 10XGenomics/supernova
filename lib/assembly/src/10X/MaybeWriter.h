// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

#ifndef _maybe_writer_h
#define _maybe_writer_h
class MaybeWriter {
public:
     MaybeWriter( bool actually_write,
               bool quiet = false ) : _write(actually_write), _quiet(quiet) {};

     template <class T>
     void writeFile( char const* filename, T const& obj ) {
          if ( _write ) BinaryWriter::writeFile(filename, obj);
          else if ( !_quiet ) {
               cout << Date() << ": NOT writing " << filename
                    << " due to settings." << endl;
          }
     }

     template <class T>
     void writeFile( String const& filename, T const& obj ) {
          writeFile( filename.c_str(), obj );
     }

private:
     bool const _write;
     bool const _quiet;
};
#endif
