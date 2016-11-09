/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////
/*
 * \file MD5.h
 * \author tsharpe
 * \date Dec 8, 2008
 *
 * \brief MD5 hashing
 * Derived from the RSA Data Security, Inc. MD5 Message-Digest Algorithm
 */

#ifndef MD5_H_
#define MD5_H_

#include <cstring>

/**
 * \class MD5
 * \brief Calculates the MD5 message digest.
 *
 * Use it by making an instance, calling update repeatedly to present all your
 * bytes to be hashed, and then calling getDigest or getHexDigest.  You can
 * reuse your instance to calculate a digest of some other, new sequence of
 * bytes by calling reInit, and then update and get(Hex)Digest as before.
 *
 * Sorry about the stateful interface, but it's part of the algorithm:  Once
 * you've called getDigest or getHexDigest, the object becomes "frozen", and
 * further updates have no effect whatever.  You may call reInit to unfreeze
 * the object, but that discards the value of the current digest and starts afresh.
 * The only way to get the value of the digest of some initial portion of your
 * input, and then continue digestion, is to copy the MD5 object.  Not clear why
 * you'd want to do this, but it would work.
 */
class MD5
{
public:
    typedef unsigned int uint;
    typedef unsigned long ulong;

    /// Make a new instance.
    MD5()
      : mEnd(mBuf+sizeof(mBuf))
    {
        reInit();
    }

    MD5( MD5 const& that )
    {
        memcpy(this,&that,sizeof(MD5));
        mEnd = mBuf + sizeof(mBuf);
        mPut = mBuf + (that.mPut - that.mBuf);
    }

    MD5& operator=( MD5 const& that )
    {
        if ( this != &that )
        {
            memcpy(this,&that,sizeof(MD5));
            mEnd = mBuf + sizeof(mBuf);
            mPut = mBuf + (that.mPut - that.mBuf);
        }
        return *this;
    }

    /// Re-initialize the hash.  This allows you to use the same instance
    /// to digest a new byte sequence.  (But it's cheap to just make a
    /// new instance, if you'd prefer.)
    void reInit()
    {
      mCount = 0;
      mState[0] = 0x67452301;
      mState[1] = 0xefcdab89;
      mState[2] = 0x98badcfe;
      mState[3] = 0x10325476;
      mPut = mBuf;
      mFinalized = false;
    }

    /// Add some bytes to be hashed.
    void update( char const* buf, uint len )
    {
      if ( !mFinalized )
      {
        if ( mPut == mBuf )
        {
          updateWholeBuffers(buf,len);
        }
        else
        {
          while ( len-- > 0 )
          {
            *mPut++ = *buf++;
            if ( mPut == mEnd )
            {
              transform(mBuf);
              updateWholeBuffers(buf,len);
            }
          }
        }
      }
    }

    /// Update, adding a single byte.
    void update( char const& chr )
    {
      *mPut++ = chr;
      if ( mPut == mEnd )
      {
        transform(mBuf);
      }
    }

    /// Add some bytes to be hashed.
    template<class Itr> void update( Itr begin, Itr end )
    { while ( begin != end ) { update(*begin++); } }

    /// Returns the MD5 hash.
    /// This is a pointer to precisely 16 bytes of memory internal to this class.  Don't delete it.
    /// This method "freezes" the object so that further updates are ignored.
    char* getDigest();

    /// Returns the MD5 hash as a C-style string of hex digits.
    /// This is a pointer to 32-bytes of digits + a null internal to this class.  Don't delete it.
    /// This method "freezes" the object so that further updates are ignored.
    char* getHexDigest();

    /// Returns the number of bytes hashed.
    /// This method is valid both pre- and post-finalization by getDigest or getHexDigest.
    ulong getCount() const
    { return mCount/8 + (mPut-mBuf); }

private:
    void updateWholeBuffers( char const* buf, uint len );
    void transform( char const* block );

    /* rotates x left n bits. */
    static uint ROL( uint x, uint n ) { return (x << n) | (x >> (32-n)); }

    static uint FF( uint a, uint b, uint c, uint d, uint x, uint s )
    { return ROL(a+((b&c)|(~b&d))+x,s) + b; }

    static uint GG( uint a, uint b, uint c, uint d, uint x, uint s )
    { return ROL(a+((b&d)|(c&~d))+x,s) + b; }

    static uint HH( uint a, uint b, uint c, uint d, uint x, uint s )
    { return ROL(a+((b^c)^d)+x,s) + b; }

    static uint II( uint a, uint b, uint c, uint d, uint x, uint s )
    { return ROL(a+(c^(b|~d))+x,s) + b; }

    uint mState[4]; /* state (ABCD) */
    ulong mCount; /* number of bits hashed, modulo 2^64 */
    char mBuf[64]; /* input buffer */
    char* mPut;
    char* mEnd;
    bool mFinalized;
};

#endif /* MD5_H_ */
