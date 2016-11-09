/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////
/*
 * \file SAM.h
 * \author tsharpe
 * \date Jan 14, 2009
 *
 * \brief SAM file reader
 */
#ifndef SAM_H_
#define SAM_H_

#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <iostream>
#include <iterator>
#include <stdlib.h>
#include "system/ProcBuf.h"
#include "util/Logger.h"
#include "util/RefDesc.h"

namespace SAM
{

/// Data about a set of reads
class ReadGroup
{
public:
    typedef std::string string;

    ReadGroup( string const& id, string const& sample, string const& library,
               string const& description, string const& platformUnit,
               uint predictedInsertSize, string const& seqCenter,
               time_t runDate, string const& platform, string const& flowOrder,
               string const& keySeq, string const& programs )
    : mID(id), mSample(sample), mLibrary(library), mDescription(description),
      mPlatformUnit(platformUnit), mPredictedInsertSize(predictedInsertSize),
      mSeqCenter(seqCenter), mRunDate(runDate), mPlatform(platform),
      mFlowOrder(flowOrder), mKeySeq(keySeq), mPrograms(programs)
    {}

    // compiler-provided destructor and copying are ok

    string const& getID() const
    { return mID; }

    string const& getSample() const
    { return mSample; }

    string const& getLibrary() const
    { return mLibrary; }

    string const& getDescription() const
    { return mDescription; }

    string const& getPlatformUnit() const
    { return mPlatformUnit; }

    uint getPredictedInsertSize() const
    { return mPredictedInsertSize; }

    string const& getSeqCenter() const
    { return mSeqCenter; }

    time_t getRunDate() const
    { return mRunDate; }

    string const& getPlatform() const
    { return mPlatform; }

    string const& getFlowOrder() const
    { return mFlowOrder; }

    string const& getKeySeq() const
    { return mKeySeq; }

    string const& getPrograms() const
    { return mPrograms; }

private:
    string mID;
    string mSample;
    string mLibrary;
    string mDescription;
    string mPlatformUnit;
    uint mPredictedInsertSize;
    string mSeqCenter;
    time_t mRunDate;
    string mPlatform;
    string mFlowOrder;
    string mKeySeq;
    string mPrograms;
};

/// Data about who wrote the SAM file
class Program
{
public:
    typedef std::string string;

    Program( string const& id, string const& name, string const& version,
                string const& commandLine, string const& prevID )
    : mID(id), mName(name), mVersion(version), mCommandLine(commandLine),
      mPrevID(prevID)
    {}

    // compiler-provided destructor and copying are ok
    string const& getID() const
    { return mID; }

    string const& getName() const
    { return mName; }

    string const& getVersion() const
    { return mVersion; }

    string const& getCommandLine() const
    { return mCommandLine; }

    string const& getPrevID() const
    { return mPrevID; }

private:
    string mID;
    string mName;
    string mVersion;
    string mCommandLine;
    string mPrevID;
};

/// Data about the SAM file itself
class Header
{
public:
    typedef std::string string;

    Header()
    {}

    Header( string const& name, string const& sortOrder,
                string const& groupOrder )
    : mName(name), mSortOrder(sortOrder), mGroupOrder(groupOrder)
    {}

    // compiler-provided destructor and copying are ok

    string const& getName() const
    { return mName; }

    string const& getSortOrder() const
    { return mSortOrder; }

    string const& getGroupOrder() const
    { return mGroupOrder; }

private:
    string mName;
    string mSortOrder;
    string mGroupOrder;
};

/// A SAM record.
/// Either a line from a SAM file, or an alignment block from a BAM file.
class Record
{
public:
    typedef std::string string;
    typedef unsigned char uchar;
    typedef std::map<string,string> TagMap;

    Record()
    : mFlag(0), mpRefDesc(0), mRefPos(0), mMapQ(0), mpMateRefDesc(0),
      mMateRefPos(0), mISize(0), mpReadGroup(0)
    {}

    Record( string const& qName, uint flag, RefDesc const* pRefDesc,
            uint refPos, uchar mapQ, string const& cigar,
            RefDesc const* pMateRefDesc, uint mateRefPos, int iSize,
            string const& seq, std::vector<uchar> const& quals,
            std::vector<uchar> const& oQuals, ReadGroup const* pReadGroup,
            TagMap const& tags )
    : mQName(qName), mFlag(flag), mpRefDesc(pRefDesc), mRefPos(refPos),
      mMapQ(mapQ), mCigar(cigar), mpMateRefDesc(pMateRefDesc),
      mMateRefPos(mateRefPos), mISize(iSize), mSeq(seq), mQuals(quals),
      mOQuals(oQuals), mpReadGroup(pReadGroup), mTags(tags)
    {}

    // compiler-provided destructor and copying are ok

    string const& getQueryName() const
    { return mQName; }

    uint getFlag() const
    { return mFlag; }

    bool isPaired() const
    { return mFlag&0x0001; }
    bool isProperlyPaired() const
    { return mFlag&0x0002; }
    bool isMapped() const
    { return !(mFlag&0x0004); }
    bool isMateMapped() const
    { return isPaired() && !(mFlag&0x0008); }
    bool isReversed() const
    { return mFlag&0x0010; }
    bool isMateReversed() const
    { return isPaired() && mFlag&0x0020; }
    bool isFirstReadOfPair() const
    { return isPaired() && mFlag&0x0040; }
    bool isSecondReadOfPair() const
    { return isPaired() && mFlag&0x0080; }
    bool isPrimary() const
    { return !(mFlag&0x0100); }
    bool isPF() const
    { return !(mFlag&0x0200); }
    bool isDuplicate() const
    { return mFlag&0x0400; }

    /// this will be the empty string if isMapped() is false
    string const& getRefName() const
    { return mpRefDesc ? mpRefDesc->getName() : EMPTY; }

    /// RefDesc* will be null if isMapped() is false
    RefDesc const* getRefDesc() const
    { return mpRefDesc; }

    /// This will be 0 if isMapped() is false.
    /// First base of reference is base 1.
    uint getRefPos() const
    { return mRefPos; }

    /// This will be 0 if isMapped() is false
    ushort getMapQ() const
    { return mMapQ; }

    /// this will be the empty string if isMapped() is false
    string const& getCigar() const
    { return mCigar; }

    /// this will be the empty string if isMateMapped() is false
    string const& getMateRefName() const
    { return mpMateRefDesc ? mpMateRefDesc->getName() : EMPTY; }

    /// RefDesc* will be null if isMateMapped() is false
    RefDesc const* getMateRefDesc() const
    { return mpMateRefDesc; }

    /// This will be 0 if isMateMapped() is false
    uint getMateRefPos() const
    { return mMateRefPos; }

    /// This isn't really the insert size.
    /// As per the SAM spec, it's the 5'-to-5' pair distance on the reference.
    int getInsertSize() const
    { return mISize; }

    string const& getSequence() const
    { return mSeq; }

    /// These are the actual phred scores. I.e., they are small numbers, not characters.
    /// (The SAM offset has already been subtracted.)
    std::vector<uchar> const& getQualityScores() const
    { return mQuals; }

    /// Get the original quality scores.  (See above.)
    std::vector<uchar> const& getOriginalQualityScores() const
    { return mOQuals; }

    /// Original quality scores if they exist, otherwise recalibrated scores.
    std::vector<uchar> const& getBestQuals() const
    { return mOQuals.size() ? mOQuals : mQuals; }

    ReadGroup const* getReadGroup() const
    { return mpReadGroup; }

    // Note that optional tag names must include the type, e.g. "NM:i" not
    // just "NM".
    bool hasTag( string const& tagName ) const
    { return mTags.find(tagName) != mTags.end(); }

    string const* getTag( string const& tagName ) const
    { TagMap::const_iterator itr = mTags.find(tagName);
      return itr == mTags.end() ? 0 : &itr->second; }

    void clipRecord() {
        std::istringstream cigar(mCigar);
        size_t len;
        char tag;
        bool begin = true;
        std::ostringstream outCigar;
        string outSeq = mSeq;
        uint outRefPos = mRefPos;
        std::vector<uchar> outQuals = mQuals;
        std::vector<uchar> outOQuals = mOQuals;

        while ( cigar.good() ) {
            if (  (cigar >> len) && (cigar >> tag) ) {
                if ( begin && tag == 'H' )
                    outRefPos += len;
                else if ( begin && tag == 'S' ) {
                    outRefPos += len;
                    outSeq = outSeq.substr( len );
                    outQuals.erase( outQuals.begin(), outQuals.begin()+len );
                    if ( mOQuals.size() )
			outOQuals.erase( outOQuals.begin(), outOQuals.begin()+len );
                } else if ( tag == 'S' ) {    // && !begin
                    outSeq = outSeq.substr( 0, outSeq.size() - len );
                    outQuals.erase(outQuals.end() - len, outQuals.end() );
                    if ( mOQuals.size() )
			outOQuals.erase(outOQuals.end() - len, outOQuals.end() );
                } else if ( tag != 'H' ) {
                    begin = false;  // we allow S&H to abut
                    outCigar << len << tag;
                }
            } else if (!cigar.eof()) {
                std::cerr << "clipRecord: failed to parse cigar string: "
                    << mCigar << endl;
                return;         // hi ted - WARNING: NON-STRUCTURED CODE
            }
        }

        mRefPos = outRefPos;
        mCigar = outCigar.str();
        mSeq = outSeq;
        mQuals = outQuals;
        mOQuals = outOQuals;
    }

private:
    string mQName;
    uint mFlag;
    RefDesc const* mpRefDesc;
    uint mRefPos;
    uchar mMapQ;
    string mCigar;
    RefDesc const* mpMateRefDesc;
    uint mMateRefPos;
    int mISize;
    string mSeq;
    std::vector<uchar> mQuals;
    std::vector<uchar> mOQuals;
    ReadGroup const* mpReadGroup;
    TagMap mTags;
    static string EMPTY;
};

/// Source of Records -- a SAM file.
/// This class keeps a file descriptor open, so don't think you can create 1000
///   of these objects simultaneously.
class SAMFile
{
public:
    typedef std::string string;
    typedef std::map<string,RefDesc> RDMap;
    typedef std::map<string,uint> RGMap;

    SAMFile( string const& fileName, Logger const& logger = Logger(std::cerr) );
    SAMFile( std::istream& istrm, Logger const& logger = Logger(std::cerr) );
    virtual ~SAMFile();

    bool hasReferenceDictionary() const { return !mRefDescMap.empty(); }

    /// set the reference dictionary
    /// if you know your SAM file doesn't have the header info to supply a valid
    /// reference dictionary, but you know which dictionary you want,
    /// rather than letting nextRecord synthesize a phony entry for each novel
    /// reference ID, you can supply the correct dictionary using this method.
    /// do this before your first call to nextRecord, or you'll have a mess.
    void setReferenceDictionary( std::vector<RefDesc> const& dict )
    { typedef std::vector<RefDesc>::const_iterator Itr;
      for ( Itr itr(dict.begin()), end(dict.end()); itr != end; ++itr )
          mRefDescMap[itr->getName()] = *itr; }

    size_t getReferenceDictionarySize() const { return mRefDescMap.size(); }

    /// If all the RefDescs point to the same URI, return it.
    /// Otherwise returns an empty string.
    string getUniqueRefURI() const;

    RefDesc const* getRefDesc( string const& sequenceName ) const
    { RDMap::const_iterator pos = mRefDescMap.find(sequenceName);
      return pos == mRefDescMap.end() ? 0 : &pos->second; }

    /// Get all reference descriptors in alphabetical order by name
    std::vector<RefDesc> getRefDescs() const
    { std::vector<RefDesc> result; result.reserve(mRefDescMap.size());
      for ( auto itr=mRefDescMap.begin(),end=mRefDescMap.end();itr!=end;++itr)
          result.push_back(itr->second);
      return result; }

    std::vector<ReadGroup> const& getReadGroups() const
    { return mReadGroups; }

    ReadGroup const* getReadGroup( string const& readGroupID ) const
    { RGMap::const_iterator pos = mReadGroupMap.find(readGroupID);
      return pos == mReadGroupMap.end() ? 0 : &mReadGroups[pos->second]; }

    Header const* getHeader() const
    { return mpHeader; }

    std::vector<Program> const& getPrograms() const
    { return mPrograms; }

    virtual bool nextRecord( Record& record );

    /// where errors in parsing are being recorded
    Logger& getLogger()
    { return mLogger; }

    class iterator : public std::iterator<std::input_iterator_tag,Record>
    {
    public:
        iterator( SAMFile& file ) : mpFile(&file) { next(); }
        iterator() : mpFile(0) {}

        // compiler-supplied copying and destructor are OK

        iterator& operator++() { next(); return *this; }
        iterator& operator++(int) { next(); return *this; }
        Record const& operator*() const { return mRecord; }
        Record const* operator->() const { return &mRecord; }

        friend bool operator==( iterator const& itr1, iterator const& itr2 )
        { return itr1.mpFile == itr2.mpFile; }

        friend bool operator!=( iterator const& itr1, iterator const& itr2 )
        { return itr1.mpFile != itr2.mpFile; }

    private:
        void next() { if ( !mpFile->nextRecord(mRecord) ) mpFile = 0; }

        SAMFile* mpFile;
        Record mRecord;
    };

    iterator begin() { return iterator(*this); }
    iterator end() { return iterator(); }

protected:
    SAMFile( Logger const& logger = Logger(std::cerr) );
    void init();

private:
    SAMFile( SAMFile const& ); // unimplemented -- no copying
    SAMFile const& operator=( SAMFile const& ); // unimplemented -- no copying

    RefDesc const& addRefDesc( string const& name, string const& url, ulong len,
                                string const& hash, string const& species,
                                string const& asmID );
    std::vector<unsigned char> parseQuals( char const* quals, size_t seqSize,
                                            char const* id );

    std::ifstream mOwnedIStream;
protected:
    std::istream* mIStream;
private:
    RDMap mRefDescMap;
    std::vector<ReadGroup> mReadGroups;
    RGMap mReadGroupMap;
    std::vector<Program> mPrograms;
    Header* mpHeader;
    Logger mLogger;
    bool mHasHeaderRefDescs;

    static uint const BUF_LEN = 256*1024;
};

// A SAM::SAMFile that gets its input via Samtools. Despite the name, it can
// handle BAM or SAM files.
class BAMFile : public SAMFile
{
public:
    BAMFile(string const& fileName, string const& region = "",
                bool header_only = false,
                Logger const& logger = Logger(std::cerr) );
    ~BAMFile();
    virtual bool nextRecord(Record& record);
private:
    BAMFile(BAMFile const&); // unimplemented -- no copying
    BAMFile const& operator=(BAMFile const&); // unimplemented -- no copying
    procbuf* bam_pipe;
};


// a couple of extra classes that aren't necessary for parsing SAM
// files, but that may be helpful in interpreting them.  used for
// converting SAM to qltout, for example.


/// an alignment block.
/// honors the usual C conventions for beginnings and endings:
/// all positions are offsets (0-based).
/// ending positions are 1 past the end of the good stuff.
class Block
{
public:
    enum Eat { EAT_REF = 1, EAT_READ = 2 };
    enum BlockType { PADDING=0, READGAP=Block::EAT_REF, REFGAP=Block::EAT_READ,
                        PAIRWISE=Block::EAT_READ|Block::EAT_REF };

    Block( char op, BlockType type, uint readStart, uint refStart, uint len )
    : mReadStart(readStart), mRefStart(refStart), mLen(len), mOp(op),
      mType(type)
    {}

    // compiler-supplied copying and destructor are OK

    char getCigarOp() const
    { return mOp; }

    BlockType getType() const
    { return static_cast<BlockType>(mType); }

    bool isReadEater() const
    { return mType & EAT_READ; }

    bool isRefEater() const
    { return mType & EAT_REF; }

    uint getReadStart() const
    { return mReadStart; }

    uint getReadEnd() const
    { return mReadStart + getReadLen(); }

    uint getReadLen() const
    { return mType&EAT_READ ? mLen : 0; }

    uint getRefStart() const
    { return mRefStart; }

    uint getRefEnd() const
    { return mRefStart + getRefLen(); }

    uint getRefLen() const
    { return mType&EAT_REF ? mLen : 0; }

    uint getLength() const
    { return mLen; }

private:
    Block()
    {}

    uint mReadStart;
    uint mRefStart;
    uint mLen;
    char mOp;
    char mType;

    friend class Alignment; // uses default constructor before assigning valid values
};

/// a class for interpreting cigar strings as a sequence of alignment blocks
/// helper methods for getting padded read and ref sequence
class Alignment
{
public:
    typedef std::string string;

    Alignment( Record const&, Logger& );

    Alignment( Alignment const& alignment )
    : mBlocks(0) { *this = alignment; }

    ~Alignment() { delete [] mBlocks; }

    Alignment& operator=( Alignment const& alignment )
    { if ( this != &alignment )
      { delete [] mBlocks;
        mBlocks = new Block[alignment.mNBlocks];
        std::copy(alignment.begin(),alignment.end(),mBlocks);
        mNBlocks = alignment.mNBlocks;
        mLeftHardClip = alignment.mLeftHardClip;
        mRightHardClip = alignment.mRightHardClip;
        mRightClip = alignment.mRightClip; }
        return *this; }

    uint getReadStart() const { return front().getReadStart(); }

    uint getReadEnd() const { return back().getReadEnd(); }

    uint getRefStart() const { return front().getRefStart(); }

    uint getRefEnd() const { return back().getRefEnd(); }

    /// unaligned portion of the end of the read
    uint getRightClip() const { return mRightClip; }

    /// read length including all clipping
    uint getTotalReadLen() const { return back().getReadEnd()+mRightClip; }

    /// left clip omitted from read's sequence
    uint getLeftHardClip() const { return mLeftHardClip; }

    /// right clip omitted from read's sequence.  no. of bases clipped.
    uint getRightHardClip() const { return mRightHardClip; }

    /// difference between 1st block's readStart and last block's readEnd
    uint getAlignedReadLen() const
    { return mNBlocks ? back().getReadEnd() - front().getReadStart() : 0u; }

    /// difference between 1st block's refStart and last block's refEnd
    uint getAlignedRefLen() const
    { return mNBlocks ? back().getRefEnd() - front().getRefStart() : 0u; }

    /// get read sequence for a block.
    /// this picks out the right substring of the read's sequence for a specified block.
    /// the samSeq arg is the read sequence as it appears in the SAM file.
    /// if the block you request is a READGAP block, the string comes back as all hyphens.
    string getBlockReadSequence( uint idx, string const& samSeq )
    { Block const& block = at(idx);
      return block.isReadEater() ?
              samSeq.substr(block.getReadStart()-mLeftHardClip,block.getLength()) :
              string(block.getLength(),'-'); }

    /// get the complete, padded read sequence.
    /// this is nothing more than the concatenation of strings returned by
    /// getBlockReadSequence over all blocks.
    string getReadSequence( string const& samSeq )
    { string result;
      for ( uint iii = 0; iii < mNBlocks; ++iii )
        result += getBlockReadSequence(iii,samSeq);
      return result; }

    /// get the ref sequence for a block.
    /// this picks out the right substring of the ref sequence for a specified block.
    /// NB: the refSeq arg is *just the part of the reference sequence that's aligned*, i.e.,
    /// the half-open interval from getRefStart() to getRefEnd().
    /// we don't want you passing around 100Mb strings, now do we.
    /// if the block you request is a REFGAP block, the string comes back as all hyphens.
    string getBlockRefSequence( uint idx, string const& refSeq )
    { Block const& block = at(idx);
      return block.isRefEater() ?
              refSeq.substr(block.getRefStart()-getRefStart(),block.getLength()) :
              string(block.getLength(),'-'); }

    /// get the complete, padded reference sequence for the alignment.
    /// this is nothing more than the concatenation of strings returned by
    /// getBlockRefSequence over all blocks.
    /// see note above about the refSeq argument:  it's just the portion of the reference
    /// that we care about for this alignment.
    string getRefSequence( string const& refSeq )
    { string result;
      for ( uint iii = 0; iii < mNBlocks; ++iii )
        result += getBlockRefSequence(iii,refSeq);
      return result; }


    // the rest of this stuff makes us look like a standard (but const) STL container of Blocks

    typedef Block value_type;
    typedef uint size_type;
    typedef Block const* const_iterator;

    Block const* begin( uint idx = 0 ) const { return mBlocks + idx; }

    Block const* end() const { return mBlocks + mNBlocks; }

    Block const& front() const { return mBlocks[0]; }

    Block const& back() const { return mBlocks[mNBlocks-1]; }

    Block const& operator[]( uint idx ) const { return mBlocks[idx]; }

    Block const& at( uint idx ) const
    { if ( idx >= mNBlocks ) throw idx;
      return mBlocks[idx]; }

    uint size() const { return mNBlocks; }

private:
    void recursiveCigarParse( uint block, char const* cigarPart, uint readStart,
                                uint refStart ) throw(char const*);

    Block* mBlocks;
    uint mNBlocks;
    uint mLeftHardClip;
    uint mRightHardClip;
    uint mRightClip;
};

} // end of namespace SAM

#endif /* SAM_H_ */
