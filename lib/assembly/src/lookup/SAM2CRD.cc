/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////
/*
 * \file SAM2CRD.cc
 * \author tsharpe
 * \date Feb 12, 2009
 *
 * \brief
 *
 *
 */
#include "lookup/SAM2CRD.h"
#include "dna/Bases.h"
#include "lookup/SAM.h"
#include "math/Functions.h"
#include "random/RNGen.h"
#include "system/System.h"
#include <map>

int const PairsBuilder::DFLT_SEP;
int const PairsBuilder::DFLT_SD;

using std::istream;
using std::pair;
using std::string;
using std::cerr;
using std::map;
using namespace SAM;

namespace
{

std::string IH_TAG("IH:i");
std::string NH_TAG("NH:i");
std::string NM_TAG("NM:i");

enum RecType { UNPAIRED, PAIRED, END1, END2 };

RecType getRecType( Record rec )
{
    RecType result = UNPAIRED;
    if ( rec.isPaired() )
    {
        if ( rec.isFirstReadOfPair() )
            result = END1;
        else if ( rec.isSecondReadOfPair() )
            result = END2;
        else
            result = PAIRED;
    }
    return result;
}

typedef pair<char const*,RecType> IdMapKey;
typedef pair<uint,uint> IdMapVal; // readId, libId

struct IdMapKeyComp
{
    bool operator()( IdMapKey const& key1, IdMapKey const& key2 )
    {
        int strOrder = strcmp(key1.first,key2.first);
        return strOrder < 0 || strOrder == 0 && key1.second < key2.second;
    }
};

typedef MempoolOwner< pair<IdMapKey const,IdMapVal> > IdMapAllocator;
typedef map<IdMapKey,IdMapVal,IdMapKeyComp,IdMapAllocator> IdMap;


}

bool validateReferenceDictionary( SAMFile& sf, vec<String> const& names,
                                    vecbvec const& refSeqs )
{
    bool result = true;
    Logger& lgr = sf.getLogger();
    size_t nNames = names.size();
    if ( sf.getReferenceDictionarySize() != nNames )
    {
        lgr.log("Header has %ld sequences, but we expected %ld.",
                            sf.getReferenceDictionarySize(), nNames);
        result = false;
    }
    for ( size_t idx = 0; idx != nNames; ++idx )
    {
        RefDesc const* pRefDesc = sf.getRefDesc(names[idx]);
        if ( !pRefDesc )
        {
            lgr.log("Header has no SQ for %s.",names[idx].c_str());
            result = false;
        }
        else
        {
            if ( pRefDesc->getId() != idx )
            {
                lgr.log("Header has %s at index %d, but we expected %ld.",
                        names[idx].c_str(), pRefDesc->getId(), idx);
            }
            if ( pRefDesc->getLength() != refSeqs[idx].size() )
            {
                lgr.log("Header thinks %s has length %ld, but we expected %d.",
                        names[idx].c_str(), pRefDesc->getLength(),
                        refSeqs[idx].size());
                result = false;
            }
        }
    }
    return result;
}

void SAM2CRD( SAMFile& sf,
              vecbasevector& seqs, vecqualvector& quals,
              vec<look_align_x>& alns, vec<pairinfo>& readPairs,
              vecString& readnames, vec<Bool>& first_in_pair, vecString& libNames,
              bool mapped_pairs_only,
	      bool keep_duplicates,
              bool useOQ, bool pfOnly, bool clip, bool readnames_plus )
{
    seqs.reserve(10000000);
    quals.reserve(10000000);
    alns.reserve(10000000);
    readPairs.reserve(10000000);
    readnames.reserve(10000000);
    first_in_pair.reserve(10000000);

    uint nextReadId = seqs.size();

    IdMapAllocator alloc;
    IdMapKeyComp comp;
    IdMap idMap(comp,alloc);
    map<string,uint> libNameToLibId;

    IdMap::const_iterator pos;
    Record rec;
    if ( !sf.nextRecord(rec) )
    {
        sf.getLogger().log("No sam records.");
        return; // EARLY RETURN!
    }

    RNGen rnGen;
    do
    {
        // Check for reads that are marked as duplicates or non-PF when
        // we've been asked to do so.
        if ( (!keep_duplicates && rec.isDuplicate()) ||
             (pfOnly && !rec.isPF()) )
            continue; // NON-STRUCTURED CODE!
	if ( mapped_pairs_only && (!rec.isMapped() || !rec.isMateMapped() ) )
	     continue; // NON-STRUCTURED CODE!

        if ( clip )
            rec.clipRecord();

        IdMapKey key(rec.getQueryName().c_str(),getRecType(rec));
        uint readId;

        IdMap::const_iterator mapEnd = idMap.end();
        if ( (pos = idMap.find(key)) != mapEnd )
        {
            readId = pos->second.first;
        }
        else
        {
            if ( key.second == PAIRED )
            {
                key.second = END1;
                if ( (pos = idMap.find(key)) != mapEnd )
                {
                    key.second = END2;
                    if ( (pos = idMap.find(key)) != mapEnd )
                    {
                        sf.getLogger().log("Warning: There are more than two alignments for the paired-read query name %s.",rec.getQueryName().c_str());
                        key.second = PAIRED;
                    }
                }
            }

            readnames.push_back( rec.getQueryName( ) );
            if (readnames_plus) first_in_pair.push_back( rec.isFirstReadOfPair( ) );

            readId = nextReadId++;
            if ( key.second != PAIRED )
            {
                key.first = readnames.back().c_str();
                ReadGroup const* pRG = rec.getReadGroup();
                string libName = pRG ? pRG->getLibrary() : "unknown";
                map<string,uint>::iterator itr = libNameToLibId.find(libName);
                uint libId;
                if ( itr != libNameToLibId.end() )
                    libId = itr->second;
                else
                {
                    libNames.push_back(libName);
                    libId = libNames.size()-1;
                    libNameToLibId[libName] = libId;
                }
                idMap[key] = IdMapVal(readId,libId);
            }

            std::string const& strSeq = rec.getSequence();
            basevector seq( strSeq.begin(), strSeq.end(),
                    [&rnGen]( char base )
                    { return GeneralizedBase::fromChar(base).random(rnGen); });
            std::vector<qual_t> const& qualVec =
                    useOQ ? rec.getBestQuals() : rec.getQualityScores();
            qualvector qual(qualVec.begin(),qualVec.end());
            if ( rec.isReversed() )
            {
                seq.ReverseComplement();
                qual.ReverseMe();
            }
            seqs.push_back(seq);
            quals.push_back(qual);
        }

        if ( !rec.isMapped() )
            continue; // NON-STRUCTURED CODE!

        Alignment aln(rec,sf.getLogger());
        if ( !aln.size() )
        {
            sf.getLogger().log("Alignment for query %s has no blocks.  Skipping alignment.",rec.getQueryName().c_str());
            continue; // NON-STRUCTURED CODE!
        }

        align laAln;
        Convert_Alignment_to_align( aln, laAln );
        uint nIndels = 0;
        for ( int i = 0; i < laAln.Nblocks( ); i++ )
             nIndels += Abs( laAln.Gaps(i) );

        uint nMismatches = 0;
        if ( rec.hasTag(NM_TAG) )
        {
            char* parse_end;
            nMismatches = strtol(rec.getTag(NM_TAG)->c_str(), &parse_end, 10);
            if ( nMismatches < nIndels )
            {
                sf.getLogger().log("Record for query %s has a cigar with %u indels, but NM is only %u.  We'll record 0 mismatches.",rec.getQueryName().c_str(),nIndels,nMismatches);
                nMismatches = 0;
            }
            else
            {
                nMismatches -= nIndels;
            }
        }

        int nHits = 1;
        if ( rec.hasTag(NH_TAG) )
        {
            char* parse_end;
            nHits = strtol(rec.getTag(NH_TAG)->c_str(), &parse_end, 10);
        }
        else if ( rec.hasTag(IH_TAG) )
        {
            char* parse_end;
            nHits = strtol(rec.getTag(IH_TAG)->c_str(), &parse_end, 10);
        }
        alns.push_back(look_align_x( readId,
                                     rec.getRefDesc()->getId(),
                                     rec.getSequence().size(),
                                     rec.getRefDesc()->getLength(),
                                     rec.isReversed(),
                                     laAln,
                                     nHits,
                                     nMismatches,
                                     nIndels,
                                     rec.getMapQ() ));
    }
    while ( sf.nextRecord(rec) );

    IdMap::iterator itr = idMap.begin();
    IdMap::iterator end = idMap.end();
    if ( itr != end )
    {
        IdMap::value_type const* pLast = &*itr;
        while ( ++itr != end )
        {
            IdMap::value_type const& cur = *itr;
            if ( !strcmp(cur.first.first,pLast->first.first) )
            {
                if ( pLast->first.second == UNPAIRED )
                {
                    sf.getLogger().log("Records for query %s were marked as both unpaired and paired.  We assumed that these were actually different reads.",pLast->first.first);
                }
                else
                {
                    if ( cur.second.second != pLast->second.second )
                        sf.getLogger().log("Paired records for query %s came from two different libraries.  Choosing one arbitrarily.",pLast->first.first);
                    readPairs.push_back(pairinfo(pLast->second.first,cur.second.first,pLast->second.second));
                }
            }
            pLast = &cur;
        }
    }
    idMap.clear();
}

void Convert_Alignment_to_align( const Alignment& aln, align& a )
{    avector<int> gaps, lengths;
     bool expectPairwise = false;
     Block const* end = aln.end();
     for ( Block const* pBlock = aln.begin(); pBlock != end; ++pBlock )
     {    switch ( pBlock->getType() )
          {    case Block::PADDING:
                    break;
               case Block::READGAP:
                    if ( expectPairwise ) lengths.Append(0);
                    expectPairwise = true;
                    gaps.Append(pBlock->getLength());
                    break;
               case Block::REFGAP:
                    if ( expectPairwise ) lengths.Append(0);
                    expectPairwise = true;
                    gaps.Append(-static_cast<int>(pBlock->getLength()));
                    break;
               case Block::PAIRWISE:
                    if ( !expectPairwise ) gaps.Append(0);
                    lengths.Append(pBlock->getLength());
                    expectPairwise = false;
                    break;    }    }
     if (expectPairwise) lengths.Append(0);
     a.Set( aln.getReadStart()-aln.getLeftHardClip(),aln.getRefStart(),
          gaps,lengths );    }
