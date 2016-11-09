///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////


#ifndef FILESOUTPUTITERATOR_H
#define FILESOUTPUTITERATOR_H

#include "String.h"
#include "feudal/BinaryStream.h"
#include "system/LockedData.h"
#include <iterator>
#include <vector>

class FilesOutput
{
public:
    FilesOutput( char const* head, char const* ext, size_t maxrecs,
                    size_t firstFileNo = 0 )
    : _head(head), _ext(ext), _maxrecs(maxrecs), _fileno(firstFileNo) {}

    template <class S> // S is some kind of a stringy thing with a c_str method
    explicit FilesOutput( S const& head, S const& ext, size_t maxrecs,
                            size_t firstFileNo = 0,
                            char const*(S::*)() const=&S::c_str )
    : _head(head), _ext(ext), _maxrecs(maxrecs), _fileno(firstFileNo) {}

    ~FilesOutput() { delete _pwriter; }

    FilesOutput( FilesOutput const& )=delete;
    FilesOutput& operator=( FilesOutput const& )=delete;

    template <class T>
    void write( T const& t )
    { if ( _nrecs >= _maxrecs ) close();
      getWriter().write(t); ++_nrecs; }

    template <class T>
    void write( T const* pBeg, T const* pEnd )
    { size_t len = pEnd-pBeg;
      if ( _nrecs+len >= _maxrecs ) close();
      getWriter().write(pBeg,pEnd); _nrecs+=len; }

    void flush() { getWriter().flush(); }

    void close() { delete _pwriter; _pwriter=0; _nrecs=0; }

    template <class T>
    class Iterator  : public std::iterator<std::output_iterator_tag,T>
    {
    public:
        Iterator() : _pOut(nullptr) {} // bit-bucket constructor
        Iterator( FilesOutput* pOut ) : _pOut(pOut) {}

        Iterator& operator++() { return *this; }
        Iterator& operator++(int) { return *this; }
        Iterator& operator*() { return *this; }

        T const& operator=( T const& t )
        { if ( _pOut ) _pOut->write(t);
          return t; }

    private:
        FilesOutput* _pOut;
    };

    LockedData& getLockedData() { return _lock; }

    template <class T>
    Iterator<T> getIterator( T* ) { return Iterator<T>(this); }

    std::vector<String> getFiles() const { return _files; }

private:
    BinaryWriter& getWriter()
    { if ( !_pwriter )
      { String newName = _head + '.' + ToString(_fileno++) + '.' + _ext;
        _files.push_back(newName);
        _pwriter = new BinaryWriter(newName,false); }
      return *_pwriter; }

    String _head;
    String _ext;
    size_t _maxrecs;
    size_t _fileno;
    size_t _nrecs = 0;
    BinaryWriter* _pwriter = nullptr;
    std::vector<String> _files;
    LockedData _lock;
};


#endif	// FILESOUTPUTITERATOR_H
