/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////
/*
 * \file RefDesc.cc
 * \author tsharpe
 * \date Dec 11, 2008
 *
 * \brief Metadata for reference sequences.
 */
#include "util/RefDesc.h"
#include "util/MD5.h"
#include "dna/Bases.h"
#include "system/Assert.h"
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <cctype>

// stuff used internally, put in namespace to avoid global namespace pollution
namespace
{
    // internal class used to figure out the length and MD5 hash of a reference sequence
    // used when reading a FASTA file to produce the reference dictionary
    class RefDescBuilder
    {
    public:
        typedef std::string string;

        RefDescBuilder( string const& uri, string const& name = ANONYMOUS, string const& species = UNKNOWN, string const& asmID = UNKNOWN )
         : mName(name), mURI(uri), mSpecies(species), mAsmID(asmID)
        {}

        void processSequence( char const* buf )
        { char const* end = buf + strlen(buf);
          while ( buf != end )
          { char chr = *buf++;
            if ( GeneralizedBase::isGeneralizedBase(chr) )
              mHash.update(toupper(chr));
          }
        }

        bool isEmpty() const
        {
            return !mHash.getCount();
        }

        RefDesc getRefDesc( unsigned id ) const
        {
          return RefDesc(mName,id,mURI,mHash.getCount(),MD5(mHash).
                          getHexDigest(),mSpecies,mAsmID);
        }

    private:
        string mName;
        string mURI;
        string mSpecies;
        MD5 mHash;
        string mAsmID;
        static char const gVals[256];
        static string const ANONYMOUS;
        static string const UNKNOWN;
    };

    RefDescBuilder::string const RefDescBuilder::ANONYMOUS("anonymous");
    RefDescBuilder::string const RefDescBuilder::UNKNOWN("unknown");


    // field names in SAM sequence dictionaries
    char const NAME_FIELD[] = "SN:";
    char const LENGTH_FIELD[] = "LN:";
    char const HASH_FIELD[] = "M5:";
    char const URI_FIELD[] = "UR:";
    char const SPECIES_FIELD[] = "SP:";
    char const ASM_ID_FIELD[] = "AS:";
    char const SQ_FIELD[] = "@SQ";
    size_t const BUF_LEN = 8192;
}

RefDesc::RefDesc( string const& samSQLine )
 : mLength(0)
{
    using std::stringstream;
    stringstream iss(samSQLine);
    char fieldBuf[BUF_LEN];
    while ( iss.good() )
    {
        iss.getline(fieldBuf,sizeof(fieldBuf),'\t');
        if ( iss.gcount() && strlen(fieldBuf) )
        {
            if ( !strncmp(fieldBuf,SQ_FIELD,sizeof(SQ_FIELD)-1) )
            {
                // noise
            }
            else if ( !strncmp(fieldBuf,NAME_FIELD,sizeof(NAME_FIELD)-1) )
            {
                mName = fieldBuf + sizeof(NAME_FIELD) - 1;
            }
            else if ( !strncmp(fieldBuf,HASH_FIELD,sizeof(HASH_FIELD)-1) )
            {
                mHash = fieldBuf + sizeof(HASH_FIELD) - 1;
            }
            else if ( !strncmp(fieldBuf,URI_FIELD,sizeof(URI_FIELD)-1) )
            {
                mURI = fieldBuf + sizeof(URI_FIELD) - 1;
            }
            else if ( !strncmp(fieldBuf,SPECIES_FIELD,sizeof(SPECIES_FIELD)-1) )
            {
                mSpecies = fieldBuf + sizeof(SPECIES_FIELD) - 1;
            }
            else if ( !strncmp(fieldBuf,ASM_ID_FIELD,sizeof(ASM_ID_FIELD)-1) )
            {
                mAsmID = fieldBuf + sizeof(ASM_ID_FIELD) - 1;
            }
            else if ( !strncmp(fieldBuf,LENGTH_FIELD,sizeof(LENGTH_FIELD)-1) )
            {
                stringstream lenISS(fieldBuf + sizeof(LENGTH_FIELD) - 1);
                lenISS >> mLength;
                if ( lenISS.fail() || lenISS.bad() || !lenISS.eof() )
                {
                    FatalErr("Stumbled while reading the sequence length field of the sequence dictionary.  The line was: " << samSQLine);
                }
            }
            else
            {
                FatalErr("Stumbled while reading sequence dictionary.  Couldn't understand field: " << fieldBuf << " in the line " << samSQLine);
            }
        }
    }
    if ( !mName.size() || !mLength )
    {
        FatalErr("Need a name (SN:) and a non-zero length (LN:) at minimum.  This line didn't cut it: " << samSQLine);
    }
}

RefDict RefDesc::readFASTA( string const& fastaFileName )
{
    std::ifstream is(fastaFileName.c_str());
    if ( !is )
    { FatalErr("Can't read fasta file " << fastaFileName); }

    string line;
    RefDict vec;
    string uri("file:///"+fastaFileName);
    RefDescBuilder desc(uri);
    while ( std::getline(is,line) )
    {
        if ( !line.empty() )
        {
            if ( line[0] != '>' )
            {
                desc.processSequence(line.c_str());
            }
            else
            {
                if ( !desc.isEmpty() )
                {
                    vec.push_back(desc.getRefDesc(vec.size()));
                }

                string name(line.substr(1));
                name.erase(0,name.find_first_not_of(" \t"));
                name.resize(name.find_last_not_of(" \t")+1); // npos trick
                if ( name.find("gi|") != 0 )
                {
                    desc = RefDescBuilder(uri,name);
                }
                else
                {
                    size_t pos = name.find('|',3);
                    if ( pos == string::npos )
                        pos = name.size();
                    string uri = "gi:" + name.substr(3,pos);
                    pos = name.rfind('|');
                    name = name.substr(pos+1);
                    pos = name.rfind(", complete genome");
                    if ( pos == string::npos )
                    {
                        desc = RefDescBuilder(uri,name);
                    }
                    else
                    {
                        string species(name.substr(0,pos));
                        species.erase(0,species.find_first_not_of(" \t"));
                        species.resize(species.find_last_not_of(" \t")+1); // npos trick
                        desc = RefDescBuilder(uri,"complete genome",species);
                    }
                }
            }
        }
    }
    ForceAssert(is.eof());
    if ( !desc.isEmpty() )
    {
        vec.push_back(desc.getRefDesc(vec.size()));
    }
    return vec;
}

RefDict RefDesc::readDict( string const& dictFileName )
{
    std::ifstream is(dictFileName.c_str());
    if ( !is )
    {
        FatalErr("Can't read dictionary file " << dictFileName);
    }

    RefDict vec;
    readDict(is,vec);

    if ( !is.eof() )
    {
        FatalErr("Trouble of a mysterious nature reading the dictionary file: " + dictFileName);
    }

    is.close();

    return vec;
}

RefDict RefDesc::readDict( istream& is )
{
    RefDict vec;
    readDict(is,vec);
    return vec;
}

bool RefDesc::writeDict( string const& dictFileName,
                                RefDict const& refDescs,
                                bool noHD )
{
    std::ofstream os(dictFileName.c_str());
    if ( !os.is_open() || !os.good() )
    {
        FatalErr("Can't write dictionary file " << dictFileName);
    }
    writeDict( os, refDescs, noHD );

    os.close();

    return !os.fail() && !os.bad();
}

void RefDesc::writeDict( ostream& os, RefDict const& refDescs, bool noHD )
{
    if ( !noHD )
        os << "@HD\tVN:1.0\tSO:unsorted\n";

    RefDict::const_iterator end = refDescs.end();
    for ( RefDict::const_iterator itr = refDescs.begin(); itr != end; ++itr )
    {
        RefDesc const& desc = *itr;
        os << SQ_FIELD << '\t' << NAME_FIELD << desc.getName() << '\t' << LENGTH_FIELD << desc.getLength();
        if ( desc.getAssemblyID().size() )
        {
            os << '\t' << ASM_ID_FIELD << desc.getAssemblyID();
        }
        if ( desc.getHash().size() )
        {
            os << '\t' << HASH_FIELD << desc.getHash();
        }
        if ( desc.getURI().size() )
        {
            os << '\t' << URI_FIELD << desc.getURI();
        }
        if ( desc.getSpecies().size() )
        {
            os << '\t' << SPECIES_FIELD << desc.getSpecies();
        }
        os << '\n';
    }
}

void RefDesc::readDict( istream& is, RefDict& vec )
{
    string line;
    while ( std::getline(is,line) )
    {
        vec.push_back(RefDesc(line));
    }
}
