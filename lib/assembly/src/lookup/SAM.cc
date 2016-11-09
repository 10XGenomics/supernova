/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////
/*
 * \file SAM.cc
 * \author tsharpe
 * \date Jan 14, 2009
 *
 * \brief SAM file reader
 */
#include "lookup/SAM.h"
#include "system/System.h"
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <limits.h>

using std::string;
using std::vector;
using std::endl;
using std::map;
using std::vector;
using std::pair;
using std::istream;
using std::ostream;
using std::cerr;
typedef unsigned char uchar;

namespace
{
    void tokenize( char* buf, vector<char*>& fields, char separator = '\t' )
    {
        fields.clear();
        char* end = buf + strlen(buf);
        while ( buf <= end )
        {
            fields.push_back(buf);
            char* pTab = strchr(buf,separator);
            if ( !pTab )
            {
                break;
            }
            *pTab = 0;
            buf = pTab + 1;
        }
    }

    class CharPtrComparator
    {
    public:
        bool operator()( char const* str1, char const* str2 ) const
        { return strcmp(str1,str2) < 0; }
    };

    typedef map<char const*,char const*,CharPtrComparator> FieldMap;

    void mapFields( FieldMap& map, char const* headerRowType, vector<char*> const& fields, Logger& logger )
    {
        map.clear();
        vector<char*>::const_iterator end = fields.end();
        for ( vector<char*>::const_iterator itr = ++fields.begin(); itr != end; ++itr )
        {
            char* field = *itr;
            if ( strchr(field,':') != field+2 )
            {
                logger.log("Ignoring field from %s header row.  Can't understand: %s.",headerRowType,field);
            }
            else
            {
                field[2] = 0;
                map[field] = field+3;
            }
        }
    }

    bool parseULong( char const* buf, ulong& val )
    {
        char* ptr;
        val = strtoul(buf,&ptr,10);
        return *buf && !*ptr;
    }

    bool toTime( tm& tm, time_t& time, time_t tzOff = 0 )
    {
        bool result = false;
        time_t tmp = mktime(&tm);
        if ( tmp != -1 )
        {
            time = tmp - tzOff;
            result = true;
        }
        return result;
    }

    // does not handle time-only case, which is excluded in the SAM spec
    // also does not handle year-weeks
    // fractional time parts are ignored
    bool parseDate8601( char const* buf, time_t& time )
    {
        tm tm;
        tm.tm_mon = 0;
        tm.tm_mday = 1;
        tm.tm_hour = 0;
        tm.tm_min = 0;
        tm.tm_sec = 0;
        tm.tm_isdst = 0;

        char* ptr;
        tm.tm_year = strtoul(buf,&ptr,10);
        switch ( ptr-buf )
        {
        case 8:
            tm.tm_mday = tm.tm_year % 100;
            tm.tm_year /= 100;
            // intentional flow-through
        case 6:
            tm.tm_mon = (tm.tm_year % 100) - 1;
            tm.tm_year /= 100;
            // intentional flow-through
        case 4:
            break;
        case 7:
            tm.tm_mday = tm.tm_year % 1000;
            tm.tm_year /= 1000;
            break;
        default:
            return false;
        }

        tm.tm_year -= 1900;

        if ( !*ptr )
            return toTime(tm,time);

        if ( ptr-buf == 4 )
        {
            if ( *ptr != '-' )
                return false;

            buf = ptr+1;
            tm.tm_mon = strtoul(buf,&ptr,10) - 1;
            if ( ptr-buf != 2 )
                return false;

            if ( !*ptr )
                return toTime(tm,time);

            if ( *ptr != '-' )
                return false;

            buf = ptr+1;
            tm.tm_mday = strtoul(buf,&ptr,10);
            if ( ptr-buf != 2 )
                return false;

            if ( !*ptr )
                return toTime(tm,time);
        }

        if ( *ptr != 'T' )
            return false;

        buf = ptr+1;
        tm.tm_hour = strtoul(buf,&ptr,10);
        switch ( ptr-buf )
        {
        case 6:
            tm.tm_sec = tm.tm_hour % 100;
            tm.tm_hour /= 100;
            if ( *ptr == '.' )
            {
                buf = ptr+1;
                strtoul(buf,&ptr,10); // just eat the fractional secs
            }
        case 4:
            tm.tm_min = tm.tm_hour % 100;
            tm.tm_hour /= 100;
        case 2:
            break;
        default:
            return false;
        }

        if ( !*ptr || *ptr == 'Z' )
            return toTime(tm,time);

        if ( *ptr != '+' && *ptr != '-' )
        {
            if ( *ptr != ':' || ptr-buf != 2 )
                return false;

            buf = ptr+1;
            tm.tm_min = strtoul(buf,&ptr,10);
            if ( ptr-buf != 2 )
                return false;

            if ( !*ptr || *ptr == 'Z' )
                return toTime(tm,time);

            if ( *ptr != '+' && *ptr != '-' )
            {
                if ( *ptr != ':' )
                    return false;

                buf = ptr+1;
                tm.tm_sec = strtoul(buf,&ptr,10);
                if ( ptr-buf != 2 )
                    return false;

                if ( *ptr == '.' )
                {
                    buf = ptr+1;
                    strtoul(buf,&ptr,10); // just eat the fractional secs
                }

                if ( !*ptr || *ptr == 'Z' )
                    return toTime(tm,time);

                if ( *ptr != '+' && *ptr != '-' )
                    return false;
            }
        }

        if ( *ptr == '+' || *ptr == '-' )
        {
            buf = ptr;
            time_t tzAdjust = strtol(buf,&ptr,10);
            if ( ptr-buf == 5 )
            {
                if ( *ptr )
                    return false;
                time_t signum = 1;
                if ( tzAdjust < 0 )
                {
                    signum = -1;
                    tzAdjust = -tzAdjust;
                }
                return toTime(tm,time,signum*(tzAdjust/100*60+tzAdjust%100)*60);
            }
            else if ( ptr-buf == 3 )
            {
                if ( !*ptr )
                    return toTime(tm,time,tzAdjust*60*60);

                if ( *ptr == ':' )
                {
                    buf = ptr+1;
                    time_t tzAdjustMins = strtoul(buf,&ptr,10);
                    if ( tzAdjust < 0 )
                        tzAdjustMins = -tzAdjustMins;
                    if ( ptr-buf == 2 && !*ptr )
                        return toTime(tm,time,(tzAdjust*60+tzAdjustMins)*60);
                }
            }
        }

        return false;
    }

    string unrecognizedFields( char const*const* validFields, FieldMap& definedFields )
    {
        string result;
        while ( *validFields )
        {
            definedFields.erase(*validFields++);
        }
        if ( definedFields.size() )
        {
            string space(" ");
            string colon(":");
            FieldMap::iterator end = definedFields.end();
            for ( FieldMap::iterator itr = definedFields.begin(); itr != end; ++itr )
            {
                result += space + itr->first + colon + itr->second;
            }
        }
        return result;
    }

    string strVal( char const* val )
    {
        string result;
        if ( val )
            result = val;
        return result;
    }

    bool checkSeq( char const* seq )
    {
        static bool initialized;
        static bool okVals[256];
        if ( !initialized )
        {
            okVals['a'&0xFF] = true; okVals['c'&0xFF] = true; okVals['g'&0xFF] = true; okVals['t'&0xFF] = true;
            okVals['A'&0xFF] = true; okVals['C'&0xFF] = true; okVals['G'&0xFF] = true; okVals['T'&0xFF] = true;
            okVals['n'&0xFF] = okVals['N'&0xFF] = true;
            okVals['='&0xFF] = true; okVals['.'&0xFF] = true;
            initialized = true;
        }

        bool result = true;
        while ( *seq )
        {
            if ( !okVals[*seq++&0xFF] )
            {
                result = false;
                break;
            }
        }
        return result;
    }


    char const HDRFLD_VN[] = "VN";
    char const HDRFLD_SO[] = "SO";
    char const HDRFLD_GO[] = "GO";
    char const*const HD_HDRFLDS[] = { HDRFLD_VN, HDRFLD_SO, HDRFLD_GO, 0 };
    char const HDRFLD_SN[] = "SN";
    char const HDRFLD_LN[] = "LN";
    char const HDRFLD_AS[] = "AS";
    char const HDRFLD_M5[] = "M5";
    char const HDRFLD_UR[] = "UR";
    char const HDRFLD_SP[] = "SP";
    char const*const SQ_HDRFLDS[] = { HDRFLD_SN, HDRFLD_LN, HDRFLD_AS, HDRFLD_M5, HDRFLD_UR, HDRFLD_SP, 0 };
    char const HDRFLD_ID[] = "ID";
    char const HDRFLD_SM[] = "SM";
    char const HDRFLD_LB[] = "LB";
    char const HDRFLD_DS[] = "DS";
    char const HDRFLD_PU[] = "PU";
    char const HDRFLD_PI[] = "PI";
    char const HDRFLD_CN[] = "CN";
    char const HDRFLD_DT[] = "DT";
    char const HDRFLD_PL[] = "PL";
    char const HDRFLD_FO[] = "FO";
    char const HDRFLD_KS[] = "KS";
    char const HDRFLD_PG[] = "PG";
    char const*const RG_HDRFLDS[] = { HDRFLD_ID, HDRFLD_SM, HDRFLD_LB, HDRFLD_DS, HDRFLD_PU, HDRFLD_PI, HDRFLD_CN, HDRFLD_DT, HDRFLD_PL, HDRFLD_FO, HDRFLD_KS, HDRFLD_PG, 0 };
    char const HDRFLD_CL[] = "CL";
    char const HDRFLD_PN[] = "PN";
    char const HDRFLD_PP[] = "PP";
    char const*const PG_HDRFLDS[] = { HDRFLD_ID, HDRFLD_PN, HDRFLD_CL, HDRFLD_VN, HDRFLD_PP, 0 };
}

namespace SAM
{
string Record::EMPTY("");

SAMFile::SAMFile( string const& fileName, Logger const& logger )
: mOwnedIStream(fileName.c_str()), mIStream(&mOwnedIStream), mpHeader(0),
  mLogger(logger), mHasHeaderRefDescs(false)
{
    if ( !mOwnedIStream.is_open() || !mOwnedIStream.good() )
    {
        mLogger.log("Unable to open SAM file: %s.",fileName.c_str());
        exit(1);
    }

    init();
}

SAMFile::SAMFile( istream& istrm, Logger const& logger )
: mIStream(&istrm), mpHeader(0), mLogger(logger),  mHasHeaderRefDescs(false)
{
    init();
}

SAMFile::SAMFile( Logger const& logger )
: mpHeader(0), mLogger(logger), mHasHeaderRefDescs(false)
{
    // does nothing, goes nowhere
    return;
}

SAMFile::~SAMFile()
{
    delete mpHeader;
    if ( mOwnedIStream.is_open() )
        mOwnedIStream.close();
}

string SAMFile::getUniqueRefURI() const
{
    string result;
    typedef RDMap::const_iterator Itr;
    for ( Itr itr(mRefDescMap.begin()),end(mRefDescMap.end()); itr!=end; ++itr )
    {
        string const& uri = itr->second.getURI();
        if ( result.empty() )
            result = uri;
        else if ( result != uri )
        {
            result.clear();
            break;
        }
    }
    return result;
}

/* Open a BAM or SAM file for reading by piping in Samtools. Assumes that
   BAMs' filenames end with ".bam" and anything else is a SAM. */
BAMFile::BAMFile( string const& fileName, string const& region,
                    bool header_only, Logger const& logger )
: SAMFile(logger)
{
    if (!IsCommandInPath("samtools"))
    {
        std::cerr << "SAM::BAMFile requires that samtools be in the PATH\n";
        std::cerr << "PATH=" << getenv("PATH") << std::endl;
        exit(1);
    }

    std::string samtools_cmd = "samtools view ";
    
    if (header_only)
    {
        samtools_cmd += "-H ";
    }
    else
    {
        samtools_cmd += "-h ";
    }
    
    if (fileName.compare(fileName.size() - 4, 4, ".bam"))
    {
        samtools_cmd += "-S ";
    }
    samtools_cmd += fileName;
    if ( region.size() )
        samtools_cmd += string(" ") + region;
    bam_pipe = new procbuf(samtools_cmd.c_str(), std::ios::in, true);
    mIStream = new std::istream(bam_pipe);
    init();
    return;
}

/* Destroy the machinery for piping in Samtools. */
BAMFile::~BAMFile()
{
    delete mIStream; // SAM::SAMFile isn't responsible for this.
    delete bam_pipe;
    return;
}

/* Fetch the next record, close the pipe (and clean up the Samtools process)
   if we get an EOF. */
bool BAMFile::nextRecord(Record& record)
{
    if (!SAMFile::nextRecord(record))
    {
        bam_pipe->close();
        return false;
    }
    else
    {
        return true;
    }
}

bool SAMFile::nextRecord( Record& record )
{
    char buf[BUF_LEN];
    vector<char*> tabFields;
    while ( mIStream->good() )
    {
        if ( mIStream->getline(buf,sizeof(buf)).gcount() )
        {
            if ( !strlen(buf) ) // silently ignore empty lines
                continue; // NON-STRUCTURED CODE!

            tokenize(buf,tabFields);
            if ( tabFields.size() < 11 )
            {
                mLogger.log("Ignoring alignment for read %s with only %lu fields:",tabFields[0],static_cast<unsigned long>(tabFields.size()));
                continue; // NON-STRUCTURED CODE!
            }

            string qname(tabFields[0]);
            if ( !qname.size() )
            {
                mLogger.log("Warning: read alignment with empty qname.");
            }

            ulong flag;
            if ( !parseULong(tabFields[1],flag) )
            {
                mLogger.log("Ignoring alignment with the incomprehensible flag value: %s.",tabFields[1]);
                continue; // NON-STRUCTURED CODE!
            }
            if ( !(flag&1) && (flag&0xea) )
            {
                mLogger.log("Warning: Unpaired read %s specifies paired-read flag values: %lx",tabFields[0],flag);
            }
            if ( (flag & 0x14) == 0x14 )
            {
                mLogger.log("Warning: Unmapped read %s claims to be on the reverse strand.",tabFields[0]);
            }
            if ( (flag & 0x104) == 0x104 )
            {
                mLogger.log("Warning: Unmapped read %s is marked as not primary.",tabFields[0]);
            }
            if ( (flag & ~0x7FF) )
            {
                mLogger.log("Warning: Read %s has a flag value with unknown bits set: %s..",tabFields[0],tabFields[1]);
            }

            RefDesc const* pRefDesc = 0;
            if ( (flag&4) ) // unmapped
            {
                if ( strcmp(tabFields[2],"*") ) // but rname not *
                {
                    if ( !(flag&1) )
                    {
                        mLogger.log("Warning: Unmapped, unpaired read %s has %s as a reference name.",tabFields[0],tabFields[2]);
                    }
                    else if ( (flag&8) )
                    {
                        mLogger.log("Warning: Unmapped read %s with an unmapped mate has %s as a reference name.",tabFields[0],tabFields[2]);
                    }
                    else if ( strcmp(tabFields[6],"=") && strcmp(tabFields[2],tabFields[6]) )
                    {
                        mLogger.log("Warning: Unmapped read %s has reference name %s, whereas its mapped mate has reference name %s.",tabFields[0],tabFields[2],tabFields[6]);
                    }
                }
            }
            else
            {
                if ( !strcmp(tabFields[2],"*") )
                {
                    mLogger.log("Warning: Mapped read %s has '*' as a reference name.",tabFields[0]);
                }

                string refName(tabFields[2]);
                RDMap::iterator pos = mRefDescMap.find(refName);
                if ( pos != mRefDescMap.end() )
                {
                    pRefDesc = &pos->second;
                }
                else
                {
                    if ( mHasHeaderRefDescs )
                    {
                        mLogger.log("Warning: Mapped read %s refers to reference name %s, which wasn't present in an @SQ header.  We'll fake it.",tabFields[0],tabFields[2]);
                    }
                    string empty;
                    pRefDesc = &addRefDesc(refName,empty,1000000000,empty,empty,empty);
                }
            }

            ulong posn;
            if ( !parseULong(tabFields[3],posn) )
            {
                mLogger.log("Ignoring alignment on read %s with the incomprehensible reference position: %s.",tabFields[0],tabFields[3]);
                continue; // NON-STRUCTURED CODE!
            }
            if ( pRefDesc && posn > pRefDesc->getLength() )
            {
                mLogger.log("Warning: Alignment on read %s has a reference position that exceeds the reference length.",tabFields[0]);
            }
            if ( !(flag&4) && !posn )
            {
                mLogger.log("Ignoring alignment on read %s because it has a reference position of 0.",tabFields[0]);
                continue; // NON-STRUCTURED CODE!
            }
            if ( (flag&4) )
                posn = 0;

            ulong mapQ;
            if ( !parseULong(tabFields[4],mapQ) )
            {
                mLogger.log("Ignoring alignment on read %s with the incomprehensible map quality: %s.",tabFields[0],tabFields[4]);
                continue; // NON-STRUCTURED CODE!
            }
            if ( mapQ > 255 )
            {
                mLogger.log("Warning: Map quality on read %s exceeds 255. Setting it to 255.",tabFields[0]);
            }
            if ( (flag&4) && mapQ )
            {
                mLogger.log("Warning: Unmapped read %s has a non-zero map quality.",tabFields[0]);
            }

            string cigar;
            if ( !*tabFields[5] )
            {
                mLogger.log("Warning: Alignment on read %s has empty cigar.",tabFields[0]);
            }
            else if ( !strcmp(tabFields[5],"*") )
            {
                if ( !(flag&4) )
                {
                    mLogger.log("Warning: Mapped read %s has a cigar string of '*'.",tabFields[0]);
                }
            }
            else
            {
                if ( (flag&4) )
                {
                    mLogger.log("Warning: Unmapped read %s has a cigar string of %s.",tabFields[0],tabFields[5]);
                }
                string(tabFields[5]).swap(cigar);
            }

            RefDesc const* pMateRefDesc = 0;
            if ( !(flag&1) ) // unpaired
            {
                if ( strcmp(tabFields[6],"*") )
                {
                    mLogger.log("Warning: Unpaired read %s has %s as a mate reference name.",tabFields[0],tabFields[6]);
                }
            }
            else if ( !(flag&8) ) // mate mapped
            {
                if ( !strcmp(tabFields[6],"*") )
                {
                    mLogger.log("Warning: Paired read %s with a mapped mate has '*' as the mate reference name.",tabFields[0]);
                }

                string refName(tabFields[6]);
                if ( !strcmp(tabFields[6],"=") )
                {
                    refName = string(tabFields[2]);
                }
                RDMap::iterator pos = mRefDescMap.find(refName);
                if ( pos != mRefDescMap.end() )
                {
                    pMateRefDesc = &pos->second;
                }
                else
                {
                    if ( mHasHeaderRefDescs )
                    {
                        mLogger.log("Warning: Mapped read %s refers to reference name %s, which wasn't present in an @SQ header.  We'll fake it.",tabFields[0],tabFields[2]);
                    }
                    string empty;
                    pMateRefDesc = &addRefDesc(refName,empty,1000000000,empty,empty,empty);
                    if ( !(flag&4) )
                        pRefDesc = &mRefDescMap.find(tabFields[2])->second;
                }
            }
            else if ( (flag&4) ) // read unmapped
            {
                if ( strcmp(tabFields[6],"*") )
                {
                    mLogger.log("Warning: Paired, unmapped read %s with an unmapped mate has %s as the mate reference name.",tabFields[0],tabFields[6]);
                }
            }
            else if ( strcmp(tabFields[6],"*") && strcmp(tabFields[6],"=") && strcmp(tabFields[2],tabFields[6]) )
            {
                mLogger.log("Warning: Paired, mapped read %s with an unmapped mate has %s as the mate reference name.",tabFields[0],tabFields[6]);
            }

            ulong matePosn;
            if ( !parseULong(tabFields[7],matePosn) )
            {
                mLogger.log("Ignoring alignment on read %s with the incomprehensible mate position: %s.",tabFields[0],tabFields[7]);
                continue; // NON-STRUCTURED CODE!
            }
            if ( pMateRefDesc && matePosn > pMateRefDesc->getLength() )
            {
                mLogger.log("Warning: Alignment on read %s has a mate position that exceeds the reference length.",tabFields[0]);
            }
            if ( (flag&1) && !(flag&8) && !matePosn )
            {
                mLogger.log("Warning: Read %s has a mapped mate with a mate reference position of 0.",tabFields[0]);
            }
            if ( !(flag&1) && matePosn )
            {
                mLogger.log("Warning: Unpaired read %s has a mate reference position of %lu.",tabFields[0],matePosn);
                matePosn = 0;
            }
            if ( (flag&8) )
                matePosn = 0;

            char* ptr;
            long isize = strtol(tabFields[8],&ptr,10);
            if ( *ptr )
            {
                mLogger.log("Ignoring alignment on read %s with the incomprehensible isize: %s.",tabFields[0],tabFields[8]);
                continue; // NON-STRUCTURED CODE!
            }
            if ( isize > (1<<29) )
            {
                mLogger.log("Warning: Alignment on read %s has a huge isize that violates the spec: %ld. Resetting it to the max legal value.",tabFields[0],isize);
                isize = 1 << 29;
            }
            if ( isize < -(1<<29) )
            {
                mLogger.log("Warning: Alignment on read %s has a hugely negative isize that violates the spec: %ld. Resetting it to the min legal value.",tabFields[0],isize);
                isize = -(1 << 29);
            }
            if ( (flag&4) )
            {
                if ( isize )
                    mLogger.log("Warning: Unmapped read %s has a non-zero isize.",tabFields[0]);
            }
            else if ( (flag&8) )
            {
                if ( isize )
                    mLogger.log("Warning: Read %s with unmapped mate has a non-zero isize",tabFields[0]);
            }
            else if ( (flag&1) && pRefDesc != pMateRefDesc && isize )
            {
                mLogger.log("Warning: Read %s has a mate mapped to a different reference sequence, but isize is non-zero.",tabFields[0]);
            }

            string seq;
            if ( strcmp(tabFields[9],"*") )
            {
                if ( !checkSeq(tabFields[9]) )
                {
                    mLogger.log("Warning: Read %s has gibberish in the sequence field: %s",tabFields[0],tabFields[9]);
                }
                string(tabFields[9]).swap(seq);
            }

            vector<uchar> quals;
            if ( strcmp(tabFields[10],"*") )
                quals = parseQuals(tabFields[10],seq.size(),tabFields[0]);

            Record::TagMap tags;
            vector<char*>::iterator end = tabFields.end();
            for ( vector<char*>::iterator itr = tabFields.begin()+11;
                itr != end; ++itr )
            {
                char* tag = *itr;
                if ( strlen(tag) < 5 || tag[2] != ':' || tag[4] != ':' )
                {
                    mLogger.log("Ignoring tag %s in alignment for read %s. "
                        "It's illegally formatted.",tag,tabFields[0]);
                }
                else
                {
                    tag[4] = '\0';
                    string name(tag);
                    string value(tag+5);

                    tags.insert(pair<string const, string>(name, value));
                }
            }

            ReadGroup const* pReadGroup = 0;
            Record::TagMap::iterator tagPos = tags.find("RG:Z");
            if ( tagPos != tags.end() )
            {
                string const& groupName = tagPos->second;
                RGMap::iterator grpPos = mReadGroupMap.find(groupName);
                if ( grpPos == mReadGroupMap.end() )
                {
                    mLogger.log("Warning: Alignment for read %s refers to read group %s which doesn't exist.",tabFields[0],groupName.c_str());
                }
                else
                {
                    pReadGroup = &mReadGroups[grpPos->second];
                }
            }

            vector<uchar> oQuals;
            tagPos = tags.find("OQ:Z");
            if ( tagPos != tags.end() )
            {
                string const& oQualStr = tagPos->second;
                oQuals = parseQuals(oQualStr.c_str(),seq.size(),tabFields[0]);
            }
            record = Record(qname,flag,pRefDesc,posn,mapQ,cigar,pMateRefDesc,matePosn,isize,seq,quals,oQuals,pReadGroup,tags);

            return true; // NON-STRUCTURED CODE!
        }
    }
    return false;
}

void SAMFile::init()
{
    char buf[BUF_LEN];
    vector<char*> tabFields;
    FieldMap fields;
    while (mIStream->good() && mIStream->peek() == '@')
    {
        if (mIStream->getline(buf,sizeof(buf)).gcount() &&
            strncmp(buf, "@CO", 3))
        {
            tokenize(buf,tabFields);
            char const* firstField = tabFields[0];
            mapFields(fields,firstField,tabFields,mLogger);

            if ( !strcmp(firstField,"@HD") )
            {
                if ( mpHeader )
                {
                    mLogger.log("Ignoring duplicate @HD header row.");
                    continue; // NON-STRUCTURED CODE!
                }
                if ( fields.find(HDRFLD_VN) == fields.end() )
                {
                    mLogger.log("Ignoring absence of VN field in @HD header row.");
                }
                mpHeader = new Header(strVal(fields[HDRFLD_VN]),strVal(fields[HDRFLD_SO]),strVal(fields[HDRFLD_GO]));

                string badFields(unrecognizedFields(HD_HDRFLDS,fields));
                if ( fields.size() )
                    mLogger.log("The @HD header row had these fields that we didn't understand:%s.",badFields.c_str());
            }
            else if ( !strcmp(firstField,"@SQ") )
            {
                if ( fields.find(HDRFLD_SN) == fields.end() )
                {
                    mLogger.log("Ignoring @SQ header row without a sequence name (SN) field.");
                    continue; // NON-STRUCTURED CODE!
                }
                string name(fields[HDRFLD_SN]);
                if ( fields.find(HDRFLD_LN) == fields.end() )
                {
                    mLogger.log("@SQ header row for %s has no sequence length (LN) field.  Setting length to 1Gb just to try to stumble on.",name.c_str());
                    fields[HDRFLD_LN] = "1000000000";
                }
                ulong seqLen = 0;
                if ( !parseULong(fields[HDRFLD_LN],seqLen) )
                {
                    seqLen = 1000000000;
                }
                if ( mRefDescMap.find(name) != mRefDescMap.end() )
                {
                    mLogger.log("Ignoring duplicate @SQ header row for %s.",name.c_str());
                    continue; // NON-STRUCTURED CODE!
                }

                mHasHeaderRefDescs = true;
                addRefDesc(name,strVal(fields[HDRFLD_UR]),seqLen,strVal(fields[HDRFLD_M5]),strVal(fields[HDRFLD_SP]),strVal(fields[HDRFLD_AS]));

                string badFields(unrecognizedFields(SQ_HDRFLDS,fields));
                if ( fields.size() )
                    mLogger.log("The @SQ header row for %s had these fields that we didn't understand:%s.",name.c_str(),badFields.c_str());
            }
            else if ( !strcmp(firstField,"@RG") )
            {
                if ( fields.find(HDRFLD_ID) == fields.end() )
                {
                    mLogger.log("Ignoring @RG row without an ID field.");
                    continue; // NON-STRUCTURED CODE!
                }
                string id(fields[HDRFLD_ID]);
                if ( fields.find(HDRFLD_SM) == fields.end() )
                {
                    mLogger.log("@RG header row for ID %s has no sample name.",id.c_str());
                }
                ulong insertSize = 0;
                if ( fields.find(HDRFLD_PI) != fields.end() )
                {
                    if ( !parseULong(fields[HDRFLD_PI],insertSize) )
                    {
                        mLogger.log("Couldn't parse the PI field (%s) as an insert size in @RG header row for ID %s.  Setting it to 0.",fields[HDRFLD_PI],id.c_str());
                        insertSize = 0;
                    }
                    if ( insertSize > INT_MAX )
                    {
                        mLogger.log("Insert size is unreasonably large (%lu) in @RG header row for ID %s.  Setting it to 0.",insertSize,id.c_str());
                        insertSize = 0;
                    }
                }
                time_t runDate = 0;
                if ( fields.find(HDRFLD_DT) != fields.end() )
                {
                    if ( !parseDate8601(fields[HDRFLD_DT],runDate) )
                    {
                        mLogger.log("Couldn't parse the DT field (%s) as an ISO 8601 date or date/time in @RG header row for ID %s.  Setting it to the Epoch.",fields[HDRFLD_DT],id.c_str());
                        runDate = 0;
                    }
                }
                if ( mReadGroupMap.find(id) != mReadGroupMap.end() )
                {
                    mLogger.log("Ignoring @RG header row with duplicate read group ID %s.",id.c_str());
                    continue; // NON-STRUCTURED CODE!
                }

                mReadGroupMap[id] = mReadGroups.size();
                mReadGroups.push_back(ReadGroup(id, strVal(fields[HDRFLD_SM]),
                	strVal(fields[HDRFLD_LB]), strVal(fields[HDRFLD_DS]),
                	strVal(fields[HDRFLD_PU]), insertSize,
                	strVal(fields[HDRFLD_CN]), runDate,
                	strVal(fields[HDRFLD_PL]), strVal(fields[HDRFLD_FO]),
                	strVal(fields[HDRFLD_KS]), strVal(fields[HDRFLD_PG])));

                string badFields(unrecognizedFields(RG_HDRFLDS,fields));
                if ( fields.size() )
                    mLogger.log("The @RG header row for ID %s had these fields that we didn't understand:%s.",id.c_str(),badFields.c_str());
            }
            else if ( !strcmp(firstField,"@PG") )
            {
                if ( fields.find(HDRFLD_ID) == fields.end() )
                {
                    mLogger.log("Ignoring lack of ID field in @PG header row.");
                }

                mPrograms.push_back(Program(strVal(fields[HDRFLD_ID]),
                    strVal(fields[HDRFLD_PN]), strVal(fields[HDRFLD_VN]),
                    strVal(fields[HDRFLD_CL]), strVal(fields[HDRFLD_PP])));

                string badFields(unrecognizedFields(PG_HDRFLDS,fields));
                if ( fields.size() )
                    mLogger.log("The @PG header row had these fields that we didn't understand:%s.",badFields.c_str());
            }
            else
            {
                mLogger.log("Ignoring unrecognized header row: %s.",firstField);
            }
        }
    }
}

RefDesc const& SAMFile::addRefDesc( string const& name, string const& url,
                                   ulong len, string const& hash,
                                   string const& species, string const& asmID )
{
    uint id = mRefDescMap.size();
    return mRefDescMap[name] = RefDesc(name,id,url,len,hash,species,asmID);
}

vector<uchar> SAMFile::parseQuals( char const* quals, size_t seqSize, char const* id )
{
    if ( strlen(quals) != seqSize )
    {
        mLogger.log("Warning: Sequence and quality data differ in length for alignment on read %s.",id);
    }

    vector<uchar> result;
    result.reserve(seqSize);

    bool qualTrouble = false;
    while ( *quals )
    {
        if ( (*quals & 0xff) < 33 )
            qualTrouble = true;
        result.push_back(*quals++ - 33);
    }
    if ( qualTrouble )
    {
        mLogger.log("Warning: Read %s has out-of-range quality values.",id);
    }

    return result;
}

Alignment::Alignment( Record const& rec, Logger& logger )
: mBlocks(0), mNBlocks(0), mLeftHardClip(0), mRightHardClip(0), mRightClip(0)
{
    string const& cigar = rec.getCigar();
    if ( cigar.size() )
    {
        try
        {
            recursiveCigarParse(0,cigar.c_str(),0,rec.getRefPos()-1);
            
            if (mNBlocks == 0) 
            {
                // This code handles the special case where all bases were
                // hard or soft-clipped, but we want to produce a valid data 
                // structure that will return the right length information.
                // This is accomplished by inserting a zero-length padding block
                // (which will not consume either read or reference) at the
                // [left hard clip] + [end of sequence] position.
                // Note that this padding block has a blank space as the CIGAR
                // operator, which will hopefully clue anyone in that this
                // isn't really from the CIGAR string.
                mBlocks = new Block[1];
                mBlocks[0] = Block(' ', Block::PADDING, rec.getSequence().size()
                	+ mLeftHardClip, rec.getRefPos() - 1, 0);
                mNBlocks = 1;
            }
        }
        catch ( char const* msg )
        {
            std::cerr << "Fatal error:\n";
            std::cerr << "Read " << rec.getQueryName() << " has a bad CIGAR "
                "string: " << cigar << ".\n";
            exit(1);
        }

        uint expectedSeqLen = getTotalReadLen() - mLeftHardClip;
        uint seqLen = rec.getSequence().size();
        if ( expectedSeqLen != seqLen )
        {
            logger.log("Read %s has a sequence length of %u, but we expected %u from the cigar string: %s.",rec.getQueryName().c_str(),seqLen,expectedSeqLen,cigar.c_str());
        }
    }
}

void Alignment::recursiveCigarParse( uint blockIdx, char const* cigarPart, uint readStart, uint refStart ) throw (char const*)
{
    if ( !cigarPart[0] )
    {
        mNBlocks = blockIdx;
        if (blockIdx > 0)
        {
            mBlocks = new Block[blockIdx];
        }
        return;
    }

    char* ptr;
    ulong len = strtoul(cigarPart,&ptr,10);
    char op = *ptr++;
    Block::BlockType type;

    switch ( op )
    {
    case 'M':
        type = Block::PAIRWISE;
        break;

    case 'I':
        type = Block::REFGAP;
        break;

    case 'N':
    case 'D':
        type = Block::READGAP;
        break;

    case 'P':
        type = Block::PADDING;
        break;

    case 'H':
        if (!*ptr)
        {
            mRightHardClip = len;
        }
        else if (blockIdx == 0)
        {
            mLeftHardClip = len;

            // consume a succeeding soft-clip (so we know that any
            // soft-clips not followed by a hard-clip or a cigar-end
            // are illegal)
            char* ptr2;
            ulong len2 = strtoul(ptr, &ptr2, 10);
            char op2 = *ptr2++;
            if (op2 == 'S')
            {
                ptr = ptr2;
                len += len2;
            }
        }
        else
        {
            throw "Clipping operators cannot be internal to a cigar string.";
        }        
        recursiveCigarParse(blockIdx, ptr, readStart + len, refStart);
        return;
        break;
        
    case 'S':
        if (*ptr && blockIdx != 0)
        {
            // check that a hard-clip follows if this soft-clip is not at
            // the start or end of the cigar (if a hard-clip preceded this
            // soft-clip, it would have already been consumed - see above)
            char *ptr2;
            strtoul(ptr, &ptr2, 10);
            char op2 = *ptr2++;
            if (op2 != 'H')
            {
                throw "Internal soft clips must be preceded or followed by "
                    "hard clips.";
            }
        }
        
        if (blockIdx != 0)
        {
            mRightClip = len;
        }

        // clip operands don't produce an alignment block
        recursiveCigarParse(blockIdx,ptr,readStart+len,refStart);
        return; // NON-STRUCTURED CODE!
        break;

    default:
        throw "Unrecognized operator in cigar string.";
        break;
    }

    Block block(op,type,readStart,refStart,len);
    recursiveCigarParse(blockIdx+1,ptr,block.getReadEnd(),block.getRefEnd());
    mBlocks[blockIdx] = block;
}
}
