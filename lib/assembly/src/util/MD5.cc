/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////
/*
 * \file MD5.cc
 * \author tsharpe
 * \date Dec 8, 2008
 *
 * \brief
 *
 *
 */
#include "util/MD5.h"

void MD5::updateWholeBuffers( char const* buf, uint len )
{
    if ( !(reinterpret_cast<ulong>(buf) & 3) ) // are we aligned sufficiently to just work out of the input buffer directly?
    {
        while ( len >= sizeof(mBuf) )
        {
            transform(buf);
            buf += sizeof(mBuf);
            len -= sizeof(mBuf);
        }
    }
    else // not aligned sufficiently:  copy input buffer to our own buffer
    {
        while ( len >= sizeof(mBuf) )
        {
            memcpy(mBuf,buf,sizeof(mBuf));
            transform(mBuf);
            buf += sizeof(mBuf);
            len -= sizeof(mBuf);
        }
    }
    if ( len )
    {
        memcpy(mBuf,buf,len);
        mPut = mBuf + len;
    }
}

void MD5::transform ( char const* buf )
{
    enum {S11 = 7, S12 = 12, S13 = 17, S14 = 22,
          S21 = 5, S22 = 9, S23 = 14, S24 = 20,
          S31 = 4, S32 = 11, S33 = 16, S34 = 23,
          S41 = 6, S42 = 10, S43 = 15, S44 = 21};

    uint a = mState[0];
    uint b = mState[1];
    uint c = mState[2];
    uint d = mState[3];
    uint const* x = reinterpret_cast<uint const*>(buf);

    /* Round 1 */
    a = FF(a, b, c, d, x[ 0]+0xd76aa478, S11); /* 1 */
    d = FF(d, a, b, c, x[ 1]+0xe8c7b756, S12); /* 2 */
    c = FF(c, d, a, b, x[ 2]+0x242070db, S13); /* 3 */
    b = FF(b, c, d, a, x[ 3]+0xc1bdceee, S14); /* 4 */
    a = FF(a, b, c, d, x[ 4]+0xf57c0faf, S11); /* 5 */
    d = FF(d, a, b, c, x[ 5]+0x4787c62a, S12); /* 6 */
    c = FF(c, d, a, b, x[ 6]+0xa8304613, S13); /* 7 */
    b = FF(b, c, d, a, x[ 7]+0xfd469501, S14); /* 8 */
    a = FF(a, b, c, d, x[ 8]+0x698098d8, S11); /* 9 */
    d = FF(d, a, b, c, x[ 9]+0x8b44f7af, S12); /* 10 */
    c = FF(c, d, a, b, x[10]+0xffff5bb1, S13); /* 11 */
    b = FF(b, c, d, a, x[11]+0x895cd7be, S14); /* 12 */
    a = FF(a, b, c, d, x[12]+0x6b901122, S11); /* 13 */
    d = FF(d, a, b, c, x[13]+0xfd987193, S12); /* 14 */
    c = FF(c, d, a, b, x[14]+0xa679438e, S13); /* 15 */
    b = FF(b, c, d, a, x[15]+0x49b40821, S14); /* 16 */

    /* Round 2 */
    a = GG(a, b, c, d, x[ 1]+0xf61e2562, S21); /* 17 */
    d = GG(d, a, b, c, x[ 6]+0xc040b340, S22); /* 18 */
    c = GG(c, d, a, b, x[11]+0x265e5a51, S23); /* 19 */
    b = GG(b, c, d, a, x[ 0]+0xe9b6c7aa, S24); /* 20 */
    a = GG(a, b, c, d, x[ 5]+0xd62f105d, S21); /* 21 */
    d = GG(d, a, b, c, x[10]+0x02441453, S22); /* 22 */
    c = GG(c, d, a, b, x[15]+0xd8a1e681, S23); /* 23 */
    b = GG(b, c, d, a, x[ 4]+0xe7d3fbc8, S24); /* 24 */
    a = GG(a, b, c, d, x[ 9]+0x21e1cde6, S21); /* 25 */
    d = GG(d, a, b, c, x[14]+0xc33707d6, S22); /* 26 */
    c = GG(c, d, a, b, x[ 3]+0xf4d50d87, S23); /* 27 */
    b = GG(b, c, d, a, x[ 8]+0x455a14ed, S24); /* 28 */
    a = GG(a, b, c, d, x[13]+0xa9e3e905, S21); /* 29 */
    d = GG(d, a, b, c, x[ 2]+0xfcefa3f8, S22); /* 30 */
    c = GG(c, d, a, b, x[ 7]+0x676f02d9, S23); /* 31 */
    b = GG(b, c, d, a, x[12]+0x8d2a4c8a, S24); /* 32 */

    /* Round 3 */
    a = HH(a, b, c, d, x[ 5]+0xfffa3942, S31); /* 33 */
    d = HH(d, a, b, c, x[ 8]+0x8771f681, S32); /* 34 */
    c = HH(c, d, a, b, x[11]+0x6d9d6122, S33); /* 35 */
    b = HH(b, c, d, a, x[14]+0xfde5380c, S34); /* 36 */
    a = HH(a, b, c, d, x[ 1]+0xa4beea44, S31); /* 37 */
    d = HH(d, a, b, c, x[ 4]+0x4bdecfa9, S32); /* 38 */
    c = HH(c, d, a, b, x[ 7]+0xf6bb4b60, S33); /* 39 */
    b = HH(b, c, d, a, x[10]+0xbebfbc70, S34); /* 40 */
    a = HH(a, b, c, d, x[13]+0x289b7ec6, S31); /* 41 */
    d = HH(d, a, b, c, x[ 0]+0xeaa127fa, S32); /* 42 */
    c = HH(c, d, a, b, x[ 3]+0xd4ef3085, S33); /* 43 */
    b = HH(b, c, d, a, x[ 6]+0x04881d05, S34); /* 44 */
    a = HH(a, b, c, d, x[ 9]+0xd9d4d039, S31); /* 45 */
    d = HH(d, a, b, c, x[12]+0xe6db99e5, S32); /* 46 */
    c = HH(c, d, a, b, x[15]+0x1fa27cf8, S33); /* 47 */
    b = HH(b, c, d, a, x[ 2]+0xc4ac5665, S34); /* 48 */

    /* Round 4 */
    a = II(a, b, c, d, x[ 0]+0xf4292244, S41); /* 49 */
    d = II(d, a, b, c, x[ 7]+0x432aff97, S42); /* 50 */
    c = II(c, d, a, b, x[14]+0xab9423a7, S43); /* 51 */
    b = II(b, c, d, a, x[ 5]+0xfc93a039, S44); /* 52 */
    a = II(a, b, c, d, x[12]+0x655b59c3, S41); /* 53 */
    d = II(d, a, b, c, x[ 3]+0x8f0ccc92, S42); /* 54 */
    c = II(c, d, a, b, x[10]+0xffeff47d, S43); /* 55 */
    b = II(b, c, d, a, x[ 1]+0x85845dd1, S44); /* 56 */
    a = II(a, b, c, d, x[ 8]+0x6fa87e4f, S41); /* 57 */
    d = II(d, a, b, c, x[15]+0xfe2ce6e0, S42); /* 58 */
    c = II(c, d, a, b, x[ 6]+0xa3014314, S43); /* 59 */
    b = II(b, c, d, a, x[13]+0x4e0811a1, S44); /* 60 */
    a = II(a, b, c, d, x[ 4]+0xf7537e82, S41); /* 61 */
    d = II(d, a, b, c, x[11]+0xbd3af235, S42); /* 62 */
    c = II(c, d, a, b, x[ 2]+0x2ad7d2bb, S43); /* 63 */
    b = II(b, c, d, a, x[ 9]+0xeb86d391, S44); /* 64 */

    mState[0] += a;
    mState[1] += b;
    mState[2] += c;
    mState[3] += d;
    mCount += 8*64;
    mPut = mBuf;
}

char* MD5::getDigest()
{
    if ( !mFinalized )
    {
        // padding may alter count.  stash away the true, final count.
        long count = mCount + 8*(mPut - mBuf);

        char pad = static_cast<char>(0x80);
        update(pad);

        pad = 0;
        char* end = mBuf + 56;
        while ( mPut != end )
            update(pad);

        memcpy(end,&count,sizeof(count));
        transform(mBuf);
        memset(mBuf,0,sizeof(mBuf));

        mCount = count;
        mFinalized = true;
    }
    return reinterpret_cast<char*>(mState);
}

char* MD5::getHexDigest()
{
    if ( !mFinalized )
    {
        char* put = mBuf;
        char* digest = getDigest();
        char* end = digest + 16;
        char const alphaBase = 'a' - 10;
        while ( digest != end )
        {
            char val = *digest++;
            char digit = (val >> 4) & 0xF;
            *put++ = digit > 9 ? alphaBase+digit : '0'+digit;
            digit = val & 0xF;
            *put++ = digit > 9 ? alphaBase+digit : '0'+digit;
        }
        *put = 0;
    }
    return mBuf;
}
