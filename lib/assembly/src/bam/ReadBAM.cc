///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * ReadBAM.cc
 *
 *  Created on: Aug 29, 2014
 *      Author: tsharpe
 */

#include "bam/ReadBAM.h"
#include "IteratorRange.h"
#include "ParallelVecUtilities.h"
#include "system/LockedData.h"
#include "system/SortInPlace.h"
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <istream>
#include <thread>
// MakeDepend: library ZLIB
#include <zlib.h>
#include <iostream>

#include "Vec.h"

namespace
{

// image of a block header for a BGZF-compressed file
struct GZipHeader
{

    size_t getBlockLen() const { return mBlockSizeLessOne+1ul; }

    bool isBGZFBlock() const
    { return mID==BGZF_ID && mCompressionMethod==BGZF_CM &&
            mFlags&BGZF_FXTRA && mXLen==BGZF_XLEN &&
            mSI==BGZF_SI && mSLen==BGZF_SLEN; }

    bool hasFName() const { return mFlags&BGZF_FNAME; }
    bool hasComment() const { return mFlags&BGZF_FCMNT; }
    bool hasHdrCRC() const { return mFlags&BGZF_FHCRC; }

    uint16_t mID;
    uint8_t mCompressionMethod;
    uint8_t mFlags;
    uint8_t mModTime[4]; // actually uint32_t, but we don't want any padding
    uint8_t mExtraFlags;
    uint8_t mOpSys;
    uint16_t mXLen;
    uint16_t mSI;
    uint16_t mSLen;
    uint16_t mBlockSizeLessOne;

    static uint16_t const BGZF_ID = 0x8b1f;
    static uint8_t const BGZF_CM = 8;
    static uint8_t const BGZF_FHCRC = 0x02;
    static uint8_t const BGZF_FXTRA = 0x04;
    static uint8_t const BGZF_FNAME = 0x08;
    static uint8_t const BGZF_FCMNT = 0x10;
    static uint16_t const BGZF_XLEN = 6;
    static uint16_t const BGZF_SI = 0x4342;
    static uint16_t const BGZF_SLEN = 2;
};

// image of the end of a GZip block -- the part after the compressed data
struct GZipFooter
{
    GZipFooter( char const* pData )
    { memcpy(this,pData,sizeof(*this)); }

    GZipFooter( z_stream const& zs )
    { uint8_t const* pBuf = zs.next_out-zs.total_out;
      mCRC32 = ::crc32(gEmptyCRC,pBuf,zs.total_out);
      mISize = zs.total_out; }

    friend bool operator!=( GZipFooter const& ftr1, GZipFooter const& ftr2 )
    { return ftr1.mCRC32 != ftr2.mCRC32 || ftr1.mISize != ftr2.mISize; }

    uint32_t mCRC32;
    uint32_t mISize;
    static uint32_t gEmptyCRC;
};

uint32_t GZipFooter::gEmptyCRC = ::crc32(0,nullptr,0);

// image of the header for a BAM file
struct BAMAlignHead
{
    uint32_t mRemainingBlockSize;
    int32_t mRefID;
    int32_t mPos;
    uint8_t mNameLen;
    uint8_t mMapQ;
    uint16_t mBin;
    uint16_t mCigarLen;
    uint16_t mFlags;
    uint32_t mSeqLen;
    int32_t mMateRefID;
    int32_t mMatePos;
    int32_t mTLen;

    // length of the alignment data excluding the header
    uint32_t remainingLen() const
    { return mRemainingBlockSize - sizeof(BAMAlignHead)
                + sizeof(mRemainingBlockSize); }

    bool isFirstRead() const
    { return (mFlags & FLAG_FIRST_SEGMENT) && !(mFlags & FLAG_LAST_SEGMENT); }

    bool isSecondRead() const
    { return (mFlags & FLAG_LAST_SEGMENT) && !(mFlags & FLAG_FIRST_SEGMENT); }

    static uint16_t const FLAG_REVERSED = 0x10;
    static uint16_t const FLAG_FIRST_SEGMENT = 0x40;
    static uint16_t const FLAG_LAST_SEGMENT = 0x80;
    static uint16_t const FLAG_PF = 0x200;
    static uint16_t const FLAG_SECONDARY_ALIGNMENT = 0x900;
};

// auxiliary tags signal the data type of the tag with these characters
// a return of 0 means "variable length"
// a return of -1 means "illegal data type specifier"
inline int getTagLength( char dataType )
{
    int tagLen;
    switch ( dataType )
    {
    case 'A': case 'c': case 'C': tagLen = 1; break;
    case 's': case 'S':           tagLen = 2; break;
    case 'i': case 'I': case 'f': tagLen = 4; break;
    case 'Z': case 'H': case 'B': tagLen = 0; break;
    default:                      tagLen = -1; break;
    }
    return tagLen;
}

char const* RTFM =
".\n"
"Please validate your BAM file with standard tools such as the Picard tool\n"
" ValidateSamFile. To do this, first make sure you have\n"
"java installed on your server and in your path.\nThen find the bin directory"
" for Picard jar files, denoted here picard-bin.\n"
"Finally, execute the command\n"
"java -jar picard-bin/ValidateSamFile.jar my.bam\n"
"where my.bam is the name of your bam file.\n";

#define BAMERR(file,message)  \
     { cout << "\nBAM file " << file << message << RTFM; Scram(1); }

#if 0
// single-threaded version is slower
class BAMbuf : public std::streambuf
{
public:
    BAMbuf( String const& bamFile ) : mFR(bamFile), mBlockNo(0)
    { char* end = &mInfBuf.back(); setg(&mInfBuf.front(),end,end); }

    BAMbuf( BAMbuf const& )=delete;
    BAMbuf& operator=( BAMbuf const& )=delete;

private:
    int_type underflow() override;
    bool readBlock();

    static size_t const FIL_BUF_SIZ = 64*1024ul;
    static size_t const INF_BUF_SIZ = 256*1024ul;

    FileReader mFR;
    size_t mBlockNo;
    std::array<char,FIL_BUF_SIZ> mFilBuf;
    std::array<char,INF_BUF_SIZ> mInfBuf;
};

std::streambuf::int_type BAMbuf::underflow()
{
    while ( gptr() == egptr() )
    {
        *eback() = gptr()[-1];
        if ( !readBlock() )
            return traits_type::eof();
    }
    return traits_type::to_int_type(*gptr());
}

bool BAMbuf::readBlock()
{
    size_t nRead;
    char* filBuf = &mFilBuf.front();
    GZipHeader const& hdr = *reinterpret_cast<GZipHeader*>(filBuf);
    if ( (nRead = mFR.readSome(filBuf,sizeof(hdr))) != sizeof(hdr) )
    {
        if ( !nRead ) return false;
        BAMERR(mFR.getFilename(),
                " is corrupt.  Partial GZIP header at block " << mBlockNo+1);
    }
    ++mBlockNo;
    if ( !hdr.isBGZFBlock() )
        BAMERR(mFR.getFilename(),
                " is uninterpretable as BGZF at block " << mBlockNo);

    char* itr = filBuf + sizeof(hdr);
    char* end = filBuf + hdr.getBlockLen();
    size_t remain = end - itr;
    if ( (nRead = mFR.readSome(itr,remain)) != remain )
        BAMERR(mFR.getFilename()," is truncated in block " << mBlockNo);

    if ( hdr.hasFName() )
        itr += strlen(itr)+1;
    if ( hdr.hasComment() )
        itr += strlen(itr)+1;
    if ( hdr.hasHdrCRC() )
        itr += 2;

    end -= sizeof(GZipFooter);
    if ( itr == end )
        return true; // read an empty block
    else if ( itr > end )
        BAMERR(mFR.getFilename(),
                " has bogus block length at block " << mBlockNo);

    z_stream zs;
    zs.zalloc = nullptr;
    zs.zfree = nullptr;
    zs.opaque = nullptr;
    zs.data_type = Z_BINARY;
    zs.next_in = reinterpret_cast<uint8_t*>(itr);
    zs.avail_in = end-itr;
    zs.next_out = reinterpret_cast<uint8_t*>(&mInfBuf[1]);
    zs.avail_out = INF_BUF_SIZ-1;

    if ( ::inflateInit2(&zs,-15) != Z_OK ||
         ::inflate(&zs,Z_FINISH) != Z_STREAM_END ||
         ::inflateEnd(&zs) != Z_OK ||
         GZipFooter(end) != GZipFooter(zs) ||
         zs.avail_in )
        BAMERR(mFR.getFilename(),
                " can't be unzipped at block " << mBlockNo);

    char* gptr = reinterpret_cast<char*>(zs.next_out-zs.total_out);
    setg( &mInfBuf.front(), gptr, gptr+zs.total_out );
    return true;
}
#endif

// a streambuf that grabs the next decompressed BGZF block on underflow.
// this version uses one thread to unzip while the main thread reads aligns.
// it appears that unzipping is just a bit slower than processing aligns.
class BAMbuf : public std::streambuf
{
public:
    BAMbuf( String const& bamFile )
    : mFR(bamFile), mBlockNo(0), mBufAvailable(mLock), mBufConsumed(mLock),
      mBeg(nullptr), mNext(nullptr), mEnd(nullptr), mQuit(false),
      mUnzipThread([this](){unzipBlocks();})
    { char* end = &mInfBuf[0].back(); setg(&mInfBuf[0].front(),end,end); }

    BAMbuf( BAMbuf const& )=delete;

    ~BAMbuf()
    { if ( true )
      { Locker locker(mLock); mQuit = true; mBeg = nullptr; }
      mBufConsumed.signal(); mUnzipThread.join(); }

    BAMbuf& operator=( BAMbuf const& )=delete;

private:
    int_type underflow() override;
    void unzipBlocks();
    bool offerBuf( char* beg, char* next, char* end )
    {
        Locker locker(mLock);
        while ( mBeg )
            locker.wait(mBufConsumed);
        if ( mQuit )
            return true;
        mNext = next;
        mEnd = end;
        mBeg = beg;
        return false;
    }

    static size_t const FIL_BUF_SIZ = 64*1024ul;
    static size_t const INF_BUF_SIZ = 256*1024ul;

    FileReader mFR;
    size_t mBlockNo;
    std::array<char,FIL_BUF_SIZ> mFilBuf;
    std::array<char,INF_BUF_SIZ> mInfBuf[3];
    LockedData mLock;
    Condition mBufAvailable;
    Condition mBufConsumed;
    char* mBeg;
    char* mNext;
    char* mEnd;
    bool mQuit;
    std::thread mUnzipThread;
};

std::streambuf::int_type BAMbuf::underflow()
{
    if ( gptr() == egptr() )
    {
        char putBackChr = gptr()[-1];
        if ( true )
        {
            Locker locker(mLock);
            while ( !mBeg )
                locker.wait(mBufAvailable);
            if ( mNext == mEnd )
                return traits_type::eof();
            setg(mBeg,mNext,mEnd);
            mBeg = nullptr;
        }
        mBufConsumed.signal();
        *eback() = putBackChr;
    }
    return traits_type::to_int_type(*gptr());
}

void BAMbuf::unzipBlocks()
{
    size_t bufId = 0;
    while ( true )
    {
        char* infBuf = &mInfBuf[bufId].front();
        size_t nRead;
        char* filBuf = &mFilBuf.front();
        GZipHeader const& hdr = *reinterpret_cast<GZipHeader*>(filBuf);
        if ( (nRead = mFR.readSome(filBuf,sizeof(hdr))) != sizeof(hdr) )
        {
            if ( !nRead )
            {
                if ( !offerBuf(infBuf,infBuf,infBuf) )
                    mBufAvailable.signal();
                break;
            }
            BAMERR(mFR.getFilename(),
                    " is corrupt.  Partial GZIP header at block " << mBlockNo+1);
        }
        ++mBlockNo;
        if ( !hdr.isBGZFBlock() )
            BAMERR(mFR.getFilename(),
                    " is uninterpretable as BGZF at block " << mBlockNo);

        char* itr = filBuf + sizeof(hdr);
        char* end = filBuf + hdr.getBlockLen();
        size_t remain = end - itr;
        if ( (nRead = mFR.readSome(itr,remain)) != remain )
            BAMERR(mFR.getFilename()," is truncated in block " << mBlockNo);

        if ( hdr.hasFName() ) itr += strlen(itr)+1;
        if ( hdr.hasComment() ) itr += strlen(itr)+1;
        if ( hdr.hasHdrCRC() ) itr += 2;

        end -= sizeof(GZipFooter);
        if ( itr == end ) continue; // read an empty block
        else if ( itr > end )
            BAMERR(mFR.getFilename(),
                    " has bogus block length at block " << mBlockNo);

        z_stream zs;
        zs.zalloc = nullptr;
        zs.zfree = nullptr;
        zs.opaque = nullptr;
        zs.data_type = Z_BINARY;
        zs.next_in = reinterpret_cast<uint8_t*>(itr);
        zs.avail_in = end-itr;
        zs.next_out = reinterpret_cast<uint8_t*>(infBuf+1);
        zs.avail_out = INF_BUF_SIZ-1;

        if ( ::inflateInit2(&zs,-15) != Z_OK ||
                ::inflate(&zs,Z_FINISH) != Z_STREAM_END ||
                ::inflateEnd(&zs) != Z_OK ||
                GZipFooter(end) != GZipFooter(zs) ||
                zs.avail_in )
            BAMERR(mFR.getFilename(),
                    " can't be unzipped at block " << mBlockNo);

        if ( zs.total_out )
        {
            char* next = reinterpret_cast<char*>(zs.next_out-zs.total_out);
            if ( offerBuf(infBuf,next,next+zs.total_out) ) break;
            mBufAvailable.signal();
            if ( ++bufId == 3 ) bufId = 0;
        }
    }
}

// read over the header and dictionary part of a BAM file, leaving us
// positioned at the alignments section

void skipToAligns( std::istream& is, String const& bamFile )
{
    uint32_t val;
    if ( !is.read(reinterpret_cast<char*>(&val),sizeof(val)) )
        BAMERR(bamFile," is empty");
    if ( val != 0x014d4142 )
        BAMERR(bamFile," lacks a BAM header");
    if ( !is.read(reinterpret_cast<char*>(&val),sizeof(val)) )
        BAMERR(bamFile," header length is truncated");
    if ( !is.ignore(val) )
        BAMERR(bamFile," is truncated in header");
    uint32_t nRefs;
    if ( !is.read(reinterpret_cast<char*>(&nRefs),sizeof(nRefs)) )
        BAMERR(bamFile," is truncated at ref desc count");
    while ( nRefs-- )
    {
        if ( !is.read(reinterpret_cast<char*>(&val),sizeof(val)) )
            BAMERR(bamFile," is truncated in ref desc len");
        if ( !is.ignore(val) )
            BAMERR(bamFile," is truncated at ref desc name");
        if ( !is.read(reinterpret_cast<char*>(&val),sizeof(val)) )
            BAMERR(bamFile," is truncated in ref desc size");
    }
}

void readAligns( std::istream& is, const String& bamFile, 
        MempoolAllocator<char> const& alloc, bool pfOnly, 
        vecbasevector& reads_b, VecPQVec& reads_q, vecString& reads_n,
	   std::set<String>* tag_select,	// null == grab all tags
	   vecString* reads_tp 			// tags here if not nullptr
	   )
{    
     // Buffer for quality score compression in batches.
     // VecPQVec qualsbuf_c;
     const int qbmax = 10000000;
     vec<qvec> qualsbuf;
     for ( int i = 0; i < qbmax; i++ )
          qualsbuf.emplace_back(alloc);
     int qbcount = 0;

     BAMAlignHead alnHd;
     size_t alnNo = 0;
     RNGen rng;
     std::vector<char> nibbleSeq;
     while ( !is.eof() && is.peek() != std::istream::traits_type::eof() )
     {    ++alnNo;
          if ( !is.read(reinterpret_cast<char*>(&alnHd), sizeof(alnHd)) )
               BAMERR(bamFile," is truncated in alignment " << alnNo);
          if ( (pfOnly && !(alnHd.mFlags & BAMAlignHead::FLAG_PF)) ||
               (alnHd.mFlags & BAMAlignHead::FLAG_SECONDARY_ALIGNMENT) )
          {    if ( !is.ignore(alnHd.remainingLen()) )
                    BAMERR(bamFile," is truncated in skipped alignment " << alnNo);
               continue;    }

          reads_b.push_back( basevector(alloc) );
          reads_n.push_back( String(alloc) );
		if ( reads_tp ) reads_tp->push_back( String(alloc) );

          basevector& seq = reads_b.back();
          String& readName = reads_n.back();

          readName.reserve(alnHd.mNameLen + 1);
          readName.resize(alnHd.mNameLen - 1);
          if ( !is.read(&readName.front(), alnHd.mNameLen - 1) || !is.ignore(1) )
               BAMERR(bamFile," is truncated in read name of alignment " << alnNo);
          readName.push_back('.');
          if ( alnHd.isFirstRead() ) readName.push_back('1');
          else if ( alnHd.isSecondRead() ) readName.push_back('2');
          else readName.push_back('3');

          if ( !is.ignore(4 * alnHd.mCigarLen) )
               BAMERR(bamFile," is truncated in cigar of alignment " << alnNo);

          seq.reserve(alnHd.mSeqLen);
          nibbleSeq.resize((alnHd.mSeqLen+1)/2);
          char* nibbles = &nibbleSeq.front();
          if ( !is.read(nibbles,nibbleSeq.size()) )
              BAMERR(bamFile," is truncated in seq of alignment " << alnNo);
          char packedSeq = 0;
          for ( unsigned idx = 0; idx != alnHd.mSeqLen; ++idx )
          {   if ( idx & 1 ) packedSeq <<= 4;
              else packedSeq = *nibbles++;
              char bits = (packedSeq >> 4) & 0x0f;
              if ( !bits )
                  BAMERR(bamFile," has uninterpretable seq data in alignment "
                              << alnNo);
              seq.push_back(GeneralizedBase::fromBits(bits).random(rng));    }
  
        qvec& quals = qualsbuf[qbcount];
        qbcount++;
        quals.resize(alnHd.mSeqLen);
        char* qBuf = reinterpret_cast<char*>(&quals.front());
        if ( !is.read(qBuf,alnHd.mSeqLen) )
            BAMERR(bamFile," is truncated in quals of alignment " << alnNo);

        long auxLen = alnHd.remainingLen() - alnHd.mNameLen
                    - 4*alnHd.mCigarLen - (alnHd.mSeqLen + 1)/2 - alnHd.mSeqLen;
        while ( auxLen > 0 )
        {
            char tag[3];
            if ( !is.read(tag, sizeof(tag)) )
                BAMERR(bamFile," is truncated in tag header for alignment "
                        << alnNo);
            auxLen -= 3;
            int tagLen = getTagLength(tag[2]);
            if ( tagLen == -1 )
                BAMERR(bamFile," has bad data type in tag header for alignment "
                        << alnNo);
            if ( tag[2] == 'B' )
            {
                char dataType;
                uint32_t arrLen;
                if ( !is.get(dataType)
                        || !is.read(reinterpret_cast<char*>(&arrLen),
                                sizeof(arrLen)) )
                    BAMERR(bamFile,
                            " is truncated in B tag header for alignment "
                               << alnNo);
                tagLen = getTagLength(dataType);
                if ( tagLen <= 0 )
                    BAMERR(bamFile,
                            " has bad data type in B tag header for alignment "
                                << alnNo);
                tagLen *= arrLen;
                auxLen -= 5;
            }

            if ( tagLen )
            {
                if ( !is.ignore(tagLen) )
                    BAMERR(bamFile," is truncated in tag data for alignment "
                                << alnNo);
                auxLen -= tagLen;
            }
            else if ( tag[0] == 'O' && tag[1] == 'Q' )
            {
                if ( tag[2] != 'Z' )
                    BAMERR(bamFile," contains OQ tag with non-Z data type "
                                    "for alignment " << alnNo);
                if ( !is.read(qBuf,alnHd.mSeqLen) )
                    BAMERR(bamFile," is truncated in OQ tag"
                                    " data for alignment " << alnNo);
                for ( unsigned char& val : quals )
                    val -= 33;
                char byte;
                if ( !is.get(byte) || byte )
                    BAMERR(bamFile," contains OQ tag with the wrong length "
                                    "for alignment " << alnNo);
                auxLen -= alnHd.mSeqLen + 1;
            }
            else // has to be H or Z tag type
            {    char byte;

			  String stag(2);
			  stag.push_back(tag[0]);
			  stag.push_back(tag[1]);

			  if ( reads_tp && 
					  ( !tag_select || 
					    tag_select->find( stag ) != tag_select->end() 
					    ) ) {
				  auto& readtag = reads_tp->back();
				  readtag.push_back(tag[0]);
				  readtag.push_back(tag[1]);
				  readtag.push_back(':');
				  do
				  {    if ( !is.get(byte) )
						  BAMERR(bamFile," is truncated in null-delimited tag"
                                        " data for alignment " << alnNo);
					  auxLen -= 1;    
					  readtag.push_back(byte);
				  } while ( byte );    
			  } else {

				  do
				  {    if ( !is.get(byte) )
						  BAMERR(bamFile," is truncated in null-delimited tag"
                                        " data for alignment " << alnNo);
					  auxLen -= 1;    
				  } while ( byte );    
			  }
		  }
        }

     if ( auxLen < 0 )
          BAMERR(bamFile," has bogus alignment block len for alignment " << alnNo);

     if ( alnHd.mFlags & BAMAlignHead::FLAG_REVERSED )
     {    seq.ReverseComplement( );
          quals.ReverseMe( );    }

     // Convert a batch of quality scores.

     if ( qbcount == qbmax )
     {    convertAppendParallel( qualsbuf.begin( ), qualsbuf.begin( ) + qbcount, 
               reads_q );
          qbcount = 0;     }

    }

     // Convert final batch of quality scores.

     convertAppendParallel( qualsbuf.begin( ), qualsbuf.begin( ) + qbcount, 
          reads_q );

}

class FunnyIterator
{
public:
    FunnyIterator( qvec const** pVals ) : mpVals(pVals) {}
    qvec const& operator*() { return **mpVals; }
    qvec const* operator->() { return *mpVals; }
    FunnyIterator& operator++() { ++mpVals; return *this; }
    friend FunnyIterator operator+( FunnyIterator const& fi, size_t off )
    { return FunnyIterator(fi.mpVals+off); }
    friend size_t operator-( FunnyIterator const& fi1, FunnyIterator const& fi2 )
    { return fi1.mpVals - fi2.mpVals; }
    friend bool operator!=( FunnyIterator const& fi1, FunnyIterator const& fi2 )
    { return fi1.mpVals != fi2.mpVals; }

private:
    qvec const** mpVals;
};

template <class T> void movePairs( size_t nReads, vecbasevector& reads_b, 
     VecPQVec& reads_q, vecString& reads_n, vecString& reads_t,
	vec<T> const& readIndices, bool uniquifyNames, 
	vecbvec* pVBV, VecPQVec* pVPQV, vecString* pReadNames, vecString* pReadTags )
{
     size_t nPairs = (nReads + 1)/2;

     // populate pReadNames
     if (pReadNames)
     {    pReadNames->reserve(pReadNames->size() + nReads);
          auto prev = readIndices.begin();
          auto end = readIndices.end();
          for ( auto itr = prev+1; itr != end; ++itr,++prev )
          {    String const& name1 = reads_n[*prev];
               String const& name2 = reads_n[*itr];
               if ( name1.size() == name2.size() &&
                    std::equal(name1.begin(),name1.end()-2,name2.begin()) )
               {    if (uniquifyNames)
                    {    pReadNames->push_back(name1);
                         char& last1 = pReadNames->back().back();
                         if ( last1 == '3' ) last1 = '1';
                         pReadNames->push_back(name2);
                         char& last2 = pReadNames->back().back();
                         if ( last2 == '3' ) last2 = '2';    }
                    else
                    {    pReadNames->push_back(String(name1.begin(),name1.end()-2));
                         pReadNames->push_back(String(name2.begin(),name2.end()-2));
                             }    }
               ++prev;
               if ( ++itr == end || !--nPairs ) break;    }    }

     // pull out ids of valid (based on names) pairs
     vec<T> ids;
     {    auto prev = readIndices.begin();
          auto end = readIndices.end();
          for ( auto itr = prev+1; itr != end; ++itr,++prev )
          {    String const& name1 = reads_n[*prev];
               String const& name2 = reads_n[*itr];
               if ( name1.size() == name2.size() &&
                    std::equal(name1.begin(),name1.end()-2,name2.begin()) )
               {    ids.push_back( prev - readIndices.begin( ) );
                    ++prev;
                    if ( ++itr == end || !--nPairs ) break;    }    }    }
     Destroy(reads_n);

     // move reads
     pVBV->reserve(pVBV->size() + nReads);
     for ( int64_t i = 0; i < ids.jsize( ); i++ )
     {    pVBV->push_back( reads_b[ readIndices[ ids[i] ] ] );
          pVBV->push_back( reads_b[ readIndices[ ids[i] + 1 ] ] );    }
     Destroy(reads_b);

     // move quals
     pVPQV->reserve(pVPQV->size() + nReads);
     for ( int64_t i = 0; i < ids.jsize( ); i++ )
     {    pVPQV->push_back( reads_q[ readIndices[ ids[i] ] ] );
          pVPQV->push_back( reads_q[ readIndices[ ids[i] + 1 ] ] );    }
     Destroy(reads_q);    

	// optionally move tags
	if ( pReadTags ) {
		ForceAssertGt( reads_t.size(), 0u );
		pReadTags->reserve(pReadTags->size() + nReads);
		for ( int64_t i = 0; i < ids.jsize( ); i++ )
		{    pReadTags->push_back( reads_t[ readIndices[ ids[i] ] ] );
			pReadTags->push_back( reads_t[ readIndices[ ids[i] + 1 ] ] );    }
		Destroy(reads_q);    
	}
}

} // end of anonymous namespace


void BAMReader::readBAM( String const& bamFile,
                         vecbvec* pVBV, VecPQVec* pVPQV, vecString* pReadNames,
				     std::set<String>* pTagsSelect,
	   				vecString* pReadTags )
{
    cout << Date( ) << ": processing " << bamFile << '.' << endl;
    cout << Date( ) << ": memory in use = " << MemUsageGBString( ) 
          << ", peak = " << PeakMemUsageGBString( ) << endl;
    vecbasevector reads_b;
    VecPQVec reads_q;
    vecString reads_n;
    vecString reads_t;	// tags
    BAMbuf* pBB = new BAMbuf(bamFile);
    {    std::istream is(pBB);
         skipToAligns(is,bamFile);
         MempoolOwner<char> alloc;
         reads_b.reserve(100000000);
         reads_q.reserve(100000000);
         reads_n.reserve(100000000);
	    vecString* reads_tp = &reads_t;
	    if ( pReadTags ) reads_t.reserve(100000000);
	    else reads_tp = nullptr;
         readAligns(is, bamFile, alloc, mPFOnly, reads_b, reads_q, reads_n,
			    pTagsSelect, reads_tp);    
#if 0
	    if ( reads_tp ) {
		    cout << "tag dump" << endl;
		    for ( auto itr = reads_tp->begin(); itr != reads_tp->end(); ++itr )
			    if (itr->size()) cout << "tags: " << *itr << endl;
	    }
#endif
    }
    delete pBB;

    // mSelectFrac and mReadsToUse determine nReads
    size_t nReads = std::min(size_t(mSelectFrac*reads_b.size()),mReadsToUse);
    int64_t length_sum = 0;
    for ( int64_t i = 0; i < (int64_t) reads_b.size( ); i++ )
         length_sum += reads_b[i].size( );
    cout << Date( ) << ": there are " << ToStringAddCommas(nReads) << " reads";
    if ( reads_b.size( ) > 0 )
    {    int64_t len = length_sum / reads_b.size( );
         cout << " of mean length " << len << endl;    }
    else cout << endl;
    cout << Date( ) << ": memory in use = " << MemUsageGBString( ) 
          << ", peak = " << PeakMemUsageGBString( ) << endl;
    nReads = (nReads + 1ul) & ~1ul;
    if ( !nReads ) return;
    if ( (int64_t) reads_b.size( ) < UINT32_MAX )
    {
    	 // readIndices are numeric from 0..reads.size()
    	 vec<uint32_t> readIndices(rangeItr(0ul),rangeItr(reads_b.size()));
    	 // sort readIndices by read name
         ParallelSort( readIndices, [&reads_n]( size_t idx1, size_t idx2 )
              { return reads_n[idx1] < reads_n[idx2]; } );
         cout << Date( ) << ": reads sorted" << endl;
         cout << Date( ) << ": memory in use = " << MemUsageGBString( ) 
               << ", peak = " << PeakMemUsageGBString( ) << endl;
         movePairs( nReads, reads_b, reads_q, reads_n, reads_t, readIndices, 
              mUniquifyNames, pVBV, pVPQV, pReadNames, pReadTags );    }
    else
    {    vec<uint64_t> readIndices(rangeItr(0ul),rangeItr(reads_b.size()));
         ParallelSort( readIndices, [&reads_n]( size_t idx1, size_t idx2 )
              { return reads_n[idx1] < reads_n[idx2]; } );
         cout << Date( ) << ": reads sorted" << endl;
         cout << Date( ) << ": memory in use = " << MemUsageGBString( ) 
               << ", peak = " << PeakMemUsageGBString( ) << endl;
         movePairs( nReads, reads_b, reads_q, reads_n, reads_t, readIndices, 
              mUniquifyNames, pVBV, pVPQV, pReadNames, pReadTags );    }
    cout << Date( ) << ": data stashed in output structures" << std::endl;
    cout << Date( ) << ": memory in use = " << MemUsageGBString( ) 
          << ", peak = " << PeakMemUsageGBString( ) << endl;
}
