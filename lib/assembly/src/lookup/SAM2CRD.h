/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////
/*
 * \file SAM2CRD.h
 * \author tsharpe
 * \date Feb 12, 2009
 *
 * \brief Reads a sam file and writes time-honored CRD data structures.
 */
#ifndef LOOKUP_SAM2CRD_H_
#define LOOKUP_SAM2CRD_H_

#include "Basevector.h"
#include "Bitvector.h"
#include "PackAlign.h"
#include "PairsManager.h"
#include "Qualvector.h"
#include "RefLocus.h"
#include "Vec.h"
#include "VecString.h"
#include "lookup/LookAlign.h"
#include "lookup/SAM.h"
#include <ostream>
#include <iterator>
#include <map>

// Shell of an output iterator for SAM::Records.
// Contains everything except an assignment operator that says what to do with
// the record.  Sub-classes define this.
class RecOItr : public std::iterator<std::output_iterator_tag,SAM::Record>,
                    public std::unary_function<SAM::Record,bool>
{};

template <class RecOItr1, class RecOItr2>
class RecOItrMux : public RecOItr
{
public:
    RecOItrMux( RecOItr1 const& r1, RecOItr2 const& r2 ) : mR1(r1), mR2(r2) {}

    bool operator()( SAM::Record const& record ) const
    { return mR1(record) && mR2(record); }

    SAM::Record const& operator=( SAM::Record const& record )
    { mR1 = record; mR2 = record; return record; }

    RecOItrMux& operator*() { return *this; }
    RecOItrMux& operator++() { return *this; }
    RecOItrMux& operator++(int) { return *this; }

private:
    RecOItr1 mR1;
    RecOItr2 mR2;
};

// it would be cuter to have this be named operator+, for example, rather than
// mux, but i haven't been able to figure out how to avoid it being overly
// applicable -- we'd need to restrict its instantiation to RecOItr subtypes
// somehow.
template <class RecOItr1, class RecOItr2>
RecOItrMux<RecOItr1,RecOItr2> mux( RecOItr1 const& r1, RecOItr2 const& r2 )
{ return RecOItrMux<RecOItr1,RecOItr2>(r1,r2); }

// SAM::Record output iterator that snags the sequence and puts it into a
// bvec output iterator.  It RC's reversed sequences so that we recover the
// original sequence.
template <class OItr> // OItr is an output iterator for bvecs
class SeqGrabber : public RecOItr
{
public:
    SeqGrabber( OItr const& oItr ) : mOItr(oItr) {}

    bool operator()( SAM::Record const& record ) const { return true; }

    SAM::Record const& operator=( SAM::Record const& record )
    { std::string const& seq = record.getSequence();
      mBVec.assign(seq.begin(),seq.end(),GenCharToRandomBaseMapper());
      if ( record.isReversed() ) mBVec.ReverseComplement();
      *mOItr = mBVec; ++mOItr; return record; }

    SeqGrabber& operator*() { return *this; }
    SeqGrabber& operator++() { return *this; }
    SeqGrabber& operator++(int) { return *this; }

private:
    bvec mBVec;
    OItr mOItr;
};

// SAM::Record output iterator that snags the quality scores and puts them into
// a qvec output iterator.  It reverses the array for reversed records so that
// we recover the original vector of quality scores.
template <class OItr> // OItr is an output iterator for qvecs
class QualGrabber : public RecOItr
{
public:
    QualGrabber( OItr const& oItr ) : mOItr(oItr) {}

    bool operator()( SAM::Record const& record ) const { return true; }

    SAM::Record const& operator=( SAM::Record const& record )
    { std::vector<unsigned char> const& quals = record.getBestQuals();
      mQVec.assign(quals.begin(),quals.end());
      if ( record.isReversed() ) mQVec.ReverseMe();
      *mOItr = mQVec; ++mOItr; return record; }

    QualGrabber& operator*() { return *this; }
    QualGrabber& operator++() { return *this; }
    QualGrabber& operator++(int) { return *this; }

private:
    qvec mQVec;
    OItr mOItr;
};

// SAM::Record output iterator that emits puts a reference locus into an output
// iterator.
template <class OItr> // OItr is an output iterator for a RefLocus.
class RefLocGrabber : public RecOItr
{
public:
    RefLocGrabber( OItr const& oItr, SAM::SAMFile& samFile )
    : mOItr(oItr), mpSAMFile(&samFile) {}

    bool operator()( SAM::Record const& record ) const
    { return record.isMapped() && record.isPrimary() && record.getRefDesc() &&
                record.getRefPos() && record.getCigar().size(); }

    SAM::Record const& operator=( SAM::Record const& record )
    { SAM::Alignment aln(record,mpSAMFile->getLogger());
      unsigned refStart = aln.getRefStart();
      if ( refStart < aln.getReadStart() )
          refStart += record.getRefDesc()->getLength(); // assume circular ref
      unsigned refLen = aln.getReadStart() +
                          aln.getAlignedRefLen() +
                          aln.getRightClip();
      *mOItr = RefLocus(record.getRefDesc()->getId(),
                          refStart-aln.getReadStart(),
                          refLen,
                          record.isReversed());
      ++mOItr; return record; }

    RefLocGrabber& operator*() { return *this; }
    RefLocGrabber& operator++() { return *this; }
    RefLocGrabber& operator++(int) { return *this; }

private:
    OItr mOItr;
    SAM::SAMFile* mpSAMFile;
};

class PairsBuilder : public RecOItr
{
    class MapVal
    {public:
      MapVal() : mReadId(~0ul), mLibId(~0) {}
      MapVal( size_t readId, int libId ) : mReadId(readId), mLibId(libId) {}
      // compiler-supplied copying and destructor are OK
      bool isNull() const { return mLibId == ~0; }
      size_t getReadId() const { return mReadId; }
      int getLibId() const { return mLibId; }
     private:
      size_t mReadId;
      int mLibId; };

    PairsBuilder& operator*() { return *this; }
    PairsBuilder& operator++() { return *this; }
    PairsBuilder& operator++(int) { return *this; }

public:
    PairsBuilder( PairsManager& pm ) : mpPM(&pm), mId(pm.nReads()) {}
    ~PairsBuilder() { if ( mId > mpPM->nReads() ) mpPM->setNReads(mId); }

    bool operator()( SAM::Record const& record ) const { return true; }

    SAM::Record const& operator=( SAM::Record const& record )
    { if ( record.isPaired() )
      { MapVal& val = mMap[record.getQueryName()];
        if ( val.isNull() ) val = MapVal(mId,getLibId(record));
        else
        { AssertEq(getLibId(record),val.getLibId());
          if ( record.isFirstReadOfPair() )
            mpPM->addPairToLib(mId,val.getReadId(),val.getLibId(),true);
          else
            mpPM->addPairToLib(val.getReadId(),mId,val.getLibId(),true);
          mMap.erase(record.getQueryName()); } }
      mId += 1;
      return record; }

private:
    int getLibId( SAM::Record const& record )
    { String libName("unknown");
      SAM::ReadGroup const* pRG = record.getReadGroup();
      if ( pRG ) libName = String(pRG->getLibrary());
      int result = mpPM->libraryID(libName);
      if ( result == -1 )
        result = mpPM->addLibrary(DFLT_SEP,DFLT_SD,libName);
      return result; }

    static int const DFLT_SEP = 30;
    static int const DFLT_SD = 30;

    PairsManager* mpPM;
    size_t mId;
    std::map<std::string,MapVal> mMap;
};

/// add map quality to look_align
class look_align_x : public look_align
{
public:
    look_align_x( int query_id,
            int target_id,
            unsigned int query_length,
            unsigned int target_length,
            Bool rc1,
            const align& aln,
            int nhits,
            nmuts_t mutations,
            int indels,
            unsigned char mapQ )
    : look_align(query_id,target_id,query_length,target_length,rc1,aln,nhits,mutations,indels), mMapQ(mapQ)
    {}

    unsigned char mapQ() const
    { return mMapQ; }

private:
    unsigned char mMapQ;
};

struct pairinfo
{
    pairinfo( size_t id1, size_t id2, size_t lib )
    : readID1(id1), readID2(id2), libraryID(lib) {}

    // compiler-supplied copying and destructor are OK

    size_t readID1; // index into seq or quals or first member of pair
    size_t readID2; // index into seq or quals or first member of pair
    size_t libraryID; // index into libnames
};

/// check to see that the SAM::SAMFile header's reference dictionary has the same
/// names, in the same order, with the same sizes as what we've (probably) read
/// from some FASTA file that we think is associated with the BAM/SAM file.
bool validateReferenceDictionary( SAM::SAMFile& sf, vec<String> const& names,
                                    vecbvec const& refSeqs );

void SAM2CRD( SAM::SAMFile& sf, // sam input file
              vecbasevector& seqs, vecqualvector& quals,
              vec<look_align_x>& alns, vec<pairinfo>& readPairs,
              vecString& readnames, vec<Bool>& first_in_pair, vecString& libNames,
	      bool mapped_pairs_only = false,
              bool keep_duplicates = true, // if false, exclude reads marked as duplicates.
              bool useOQ = false, // use original quality scores
              bool pfOnly = false, 
              bool clip = false,  // rewrite seq/pos/cigar to apply clipping (S&H)
                                  // experimental 
              bool readnames_plus = false // append .1 or .2 to readname
              ); // non-PF reads ignored

/* DO NOT ADD "using namespace SAM;" HERE, PUT IT IN YOUR .CC FILE. */

void Convert_Alignment_to_align( const SAM::Alignment& aln, align& a );

#endif /* LOOKUP_SAM2CRD_H_ */
