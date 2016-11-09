///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * \file VirtualMasterVec.h
 * \author tsharpe
 * \date Aug 18, 2009
 *
 * \brief A class that gives you sequential, read-only access to a feudal file.
 *
 * The memory footprint is very small, and sequential access is very efficient,
 * involving minimal I/O.
 *
 * There is a random access capability as well, but it requires a somewhat
 * expensive file seek and buffer purge and reload.  Nonetheless, the random
 * access capability might be quite suitable for sparse access.
 *
 * The template class T must be a type that was written to a feudal file, so it
 * will necessarily have the right interface, namely a copy-constructor and a
 * constructor that takes void*'s to the raw data:
 * T( void* varData, ulong varDataNBytes, void* fixedData, A const& allocator );
 */
#ifndef FEUDAL_VIRTUALMASTERVEC_H_
#define FEUDAL_VIRTUALMASTERVEC_H_

#include "feudal/FeudalFileReader.h"
#include "feudal/Oob.h"
#include <cstddef>

template<class T>
class VirtualMasterVec
{
public:
    typedef T const value_type;
    typedef unsigned long size_type;

    class Itr
    : public std::iterator<std::input_iterator_tag,T,std::ptrdiff_t>
    {
    public:
        Itr( VirtualMasterVec<T>& vec, size_type idx )
        : mVec(vec), mIdx(idx), mLoadedIdx(~0ul) {}

        Itr( Itr const& that )
        : mVec(that.mVec), mIdx(that.mIdx), mLoadedIdx(~0UL) {}

        Itr& operator=( Itr const& that ) { mIdx = that.mIdx; return *this; }

        // compiler-supplied destructor is OK

        T const& operator*()
        { if ( mLoadedIdx != mIdx ) { mVec.load(mIdx,&mRef); mLoadedIdx = mIdx; }
          return mRef; }
        T const* operator->() { return &operator*(); }
        bool operator==( Itr const& that ) const { return mIdx == that.mIdx; }
        bool operator!=( Itr const& that ) const { return mIdx != that.mIdx; }
        Itr& operator++() { mIdx += 1; return *this; }
        Itr operator++(int) { Itr tmp(*this); mIdx += 1; return tmp; }

        bool operator<( Itr const& itr ) { return mIdx < itr.mIdx; }
        std::ptrdiff_t operator-( Itr const& itr ) { return mIdx - itr.mIdx; }
        Itr operator+( std::ptrdiff_t inc )
        { Itr tmp(*this); tmp.mIdx += inc; return tmp; }

    private:
        VirtualMasterVec& mVec;
        size_type mIdx;
        T mRef;
        size_type mLoadedIdx;
    };

    typedef Itr const_iterator;

    VirtualMasterVec( char const* path ) : mFFR(path) {}

    /// Construct from something with a c_str() member (like string or String)
    template <class C>
    explicit VirtualMasterVec( C const& path,
                                    char const*(C::*)() const=&C::c_str )
    : mFFR(path.c_str()) {}

    VirtualMasterVec( VirtualMasterVec<T> const& that ) : mFFR(that.mFFR) {}
    // compiler-supplied destructor is OK
    VirtualMasterVec& operator=( VirtualMasterVec<T> const& that )
    { if ( this != &that ) { mFFR = that.mFFR; } return *this; }

    Itr begin() { return Itr(*this,0); }
    Itr begin( size_type idx ) { return Itr(*this,idx); }
    Itr end() { return Itr(*this,mFFR.getNElements()); }

    T const front() { return obj(0); }
    T const back() { return obj(size()-1); }
    T const operator[]( size_type idx ) const { return obj(idx); }
    T const at( size_type idx ) const
    { if ( idx >= size() )
      { OutOfBoundsReporter::oob("VirtualMasterVec",idx,size()); }
      return obj(idx); }

    void load( size_type idx, T* pT ) const
    { AssertLt(idx,size());
      pT->readFeudal(mFFR.getData(idx),mFFR.getDataLen(idx),
                      mFFR.getFixedData(idx,T::fixedDataLen())); }

    size_type size() const { return mFFR.getNElements(); }
    typename T::size_type eleSize( size_type idx ) const
    { return T::interpretSize(mFFR.getFixedData(idx,T::fixedDataLen()),
                                mFFR.getDataLen(idx)); }

    /// total number of T::value_type's across all elements.
    /// this will blow up if T::value_type doesn't have a fixed external size
    size_t sizeSum() const
    { size_t eleSiz = valSizeof();
      AssertNe(eleSiz,0u); return mFFR.getDataLenTotal()/eleSiz; }

    /// sums element.size()-K+1 for each element where element.size() >= K
    size_t getKmerCount( unsigned K ) const
    { size_t result = 0;
      size_t idx = size();
      while ( idx-- )
      { size_t eleSz = eleSize(idx);
        if ( eleSz >= K ) result += eleSz-K+1; }
      return result; }

    bool empty() const { return size() == 0; }

    size_t getMappedLen() const { return mFFR.getMappedLen(); }

    VirtualMasterVec<T> clone() const { return VirtualMasterVec<T>( this->mFFR.getFilename() ); }

private:
    static size_t valSizeof()
    { return BinaryReader::externalSizeof(
                              static_cast<typename T::value_type*>(nullptr)); }

    T const obj( size_type idx ) const
    { T result; load(idx,&result); return result; }

    mutable FeudalFileReader mFFR;
};

#endif // FEUDAL_VIRTUALMASTERVEC_H_
