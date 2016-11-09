///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file LongReadPather.h
 * \author tsharpe
 * \date Jun 6, 2012
 *
 * \brief
 */
#ifndef KMERS_LONGREADPATHER_H_
#define KMERS_LONGREADPATHER_H_

#include "Basevector.h"
#include "Vec.h"
#include "dna/Bases.h"
#include "dna/CanonicalForm.h"
#include "feudal/BinaryStream.h"
#include "feudal/HashSet.h"
#include "feudal/HugeBVec.h"
#include "kmers/KMerContext.h"
#include "math/Functions.h"
#include "system/Assert.h"
#include "system/ID.h"
#include "system/SpinLockedData.h"
#include "system/SysConf.h"
#include "system/System.h"
#include "system/Worklist.h"
#include "system/file/Directory.h"
#include "system/file/File.h"
#include <algorithm>
#include <cstddef>
#include <functional>
#include <iterator>
#include <ostream>
#include <utility>
#include <vector>

namespace Long
{

extern unsigned gLogLevel;
#define LRP_LOG(L,x) (gLogLevel>=L ? (std::cout << x << std::endl) : std::cout)

typedef ID<5> KmerID;
typedef ID<4> EdgeID;

// this is what's in the dictionary
class KmerDictEntry
{
public:
    // used by the factory in creating HopscotchHashSets
    KmerDictEntry() : mForm(CanonicalForm::FWD) {}

    // used in looking up existing entries
    KmerDictEntry( size_t kmerId, size_t edgeId )
    : mKmerID(kmerId), mEdgeID(edgeId), mForm(CanonicalForm::FWD) {}

    // used while processing the reads into the dictionary
    KmerDictEntry( size_t kmerId, CanonicalForm form )
    : mKmerID(kmerId), mForm(form) {}

    // used while processing the reads into the dictionary
    KmerDictEntry(size_t kmerId, CanonicalForm form, KMerContext context)
    : mKmerID(kmerId), mContext(context), mForm(form)
    { if ( form == CanonicalForm::PALINDROME ) mContext |= context.rc(); }

    // compiler-supplied copying and destructor are OK

    size_t getKmerID() const { return mKmerID.val(); }
    size_t getEdgeID() const { return mEdgeID.val(); }
    KMerContext getContext() const { return mContext; }
    CanonicalForm getForm() const { return mForm; }
    bool isPalindrome() const { return mForm == CanonicalForm::PALINDROME; }

    bool isAssigned() const { return !mEdgeID.isNull(); }
    void assign( size_t kmerId, size_t edgeId, bool rc )
    { mKmerID.setVal(kmerId); mEdgeID.setVal(edgeId);
      if ( rc )
      { mForm = complement(mForm); mContext = mContext.rc(); } }

    void toPrev( CanonicalForm form )
    { mKmerID.setVal(mKmerID.val()-1); mForm = form; }
    void toNext( CanonicalForm form )
    { mKmerID.setVal(mKmerID.val()+1); mForm = form; }

    static bool equalKeys( KmerDictEntry const& e1, KmerDictEntry const& e2 )
    { return e1.mKmerID == e2.mKmerID && e1.mEdgeID == e2.mEdgeID; }

    // folds in any extra context info when we've seen the kmer before in some
    // other read.
    class Functor
    {
    public:
        Functor( KmerDictEntry const& entry ) : mpEntry(&entry) {}

        // compiler-supplied copying and destructor are OK

        void operator()( KmerDictEntry const& entry ) const
        { if ( entry.mKmerID != mpEntry->mKmerID )
            const_cast<KmerDictEntry&>(entry).mContext |=
                    entry.mForm == mpEntry->mForm ?
                            mpEntry->mContext : mpEntry->mContext.rc(); }

    private:
        KmerDictEntry const* mpEntry;
    };
    friend class Functor;

    class Printer
    {
    public:
        Printer( KmerDictEntry const& entry, unsigned K, HugeBVec const& hbv )
        : mEntry(entry), mK(K), mHBV(hbv) {}

        friend ostream& operator<<( ostream& os, Printer const& pr )
        { KmerDictEntry const& entry = pr.mEntry;
          os << "KId(" << entry.getKmerID() << ") EId(" << entry.getEdgeID()
             << ") Context(" << entry.getContext() << ") Form("
             << entry.getForm() << ") ";
          auto itr = pr.mHBV.begin(entry.getKmerID());
          std::transform(itr,itr+pr.mK,std::ostream_iterator<char>(os),
                          BaseToCharMapper());
          return os; }

    private:
        KmerDictEntry const& mEntry;
        unsigned mK;
        HugeBVec const& mHBV;
    };

private:
    KmerID mKmerID;
    EdgeID mEdgeID;
    KMerContext mContext;
    CanonicalForm mForm;
};

// a great big hash set of KmerDictEntry's.
// there are two HugeBVecs.  one is the concatenation of all the reads, and the
// kmerId's from entries that haven't been assigned to an edge are looked up on
// that HugeBVec.  the other is for edges, and kmerId's from entries that are on
// an edge (i.e., after the graph builder has done its thing) are looked up
// relative to that HugeBVec.
template <unsigned K>
class KmerDict
{
public:
    typedef KmerDictEntry Entry;

    // compute the hash of a KmerDictEntry
    class Hasher : public std::unary_function<Entry,size_t>
    {
    public:
        typedef unsigned char BaseCode;
        static unsigned const HALFK = (K+1)/2;

        Hasher( KmerDict const& dict ) : mpDict(&dict) {}

        // compiler-supplied copying and destructor are OK

        size_t operator()( Entry const& entry ) const
        { HugeBVec::const_iterator itr(getBV(entry).begin(entry.getKmerID()));
          return hashFront(itr) ^ hashRear(itr+K); }

        size_t hashFront( HugeBVec::const_iterator itr ) const
        { size_t hash = 0;
          unsigned nnn = HALFK;
          while ( nnn-- )
          { hash = (hash << 1) | (hash >> 63);
            hash ^= gBits0[*itr]; ++itr; }
          return hash; }

        size_t hashRear( HugeBVec::const_iterator itr ) const
        { size_t hash = 0;
          unsigned nnn = HALFK;
          while ( nnn-- )
          { hash = (hash << 1) | (hash >> 63);
            hash ^= gBits0[*--itr^3]; }
          return hash; }

        size_t stepFrontFwd( size_t hash, BaseCode pred, BaseCode succ ) const
        { hash = (hash << 1) | (hash >> 63);
          hash ^= gBits0[succ];
          hash ^= gBitsK[pred];
          return hash; }

        size_t stepRearFwd( size_t hash, BaseCode pred, BaseCode succ ) const
        { hash ^= gBitsK[succ^3];
          hash ^= gBits0[pred^3];
          hash = (hash >> 1) | (hash << 63);
          return hash; }

        size_t stepFrontRev( size_t hash, BaseCode pred, BaseCode succ ) const
        { hash ^= gBitsK[pred];
          hash ^= gBits0[succ];
          hash = (hash >> 1) | (hash << 63);
          return hash; }

        size_t stepRearRev( size_t hash, BaseCode pred, BaseCode succ ) const
        { hash = (hash << 1) | (hash >> 63);
          hash ^= gBits0[pred^3];
          hash ^= gBitsK[succ^3];
          return hash; }

        HugeBVec const& getBV( Entry const& entry ) const
        { return mpDict->getBV(entry); }

    private:
        KmerDict const* mpDict;
        static size_t const gBits0[4];
        static size_t const gBitsK[4];
    };
    friend class Hasher;

    // compute the hash of a KmerDictEntry, using a shortcut when we've just
    // computed the hash of a neighboring entry
    class IncrementalHasher : public Hasher
    {
        static unsigned const EFRONT = Hasher::HALFK;
        static unsigned const SREAR = K-Hasher::HALFK;

    public:
        IncrementalHasher( KmerDict const& dict )
        : Hasher(dict), mpHBV(0), mKmerId(0), mHashF(0), mHashR(0) {}

        // compiler-supplied copying and destructor are OK

        size_t operator()( Entry const& entry )
        { return hash(entry); }

        size_t hash( Entry const& entry )
        { HugeBVec const& hbv = this->getBV(entry);
          size_t kmerId = entry.getKmerID();
          if ( &hbv == mpHBV )
          { long dist = kmerId - mKmerId;
            switch ( dist )
            {case 3: stepF();
             case 2: stepF();
             case 1: return stepF();
             case 0: return mHashF ^ mHashR;
             case -3: stepR();
             case -2: stepR();
             case -1: return stepR(); } }
          mpHBV = &hbv;
          mKmerId = kmerId;
          HugeBVec::const_iterator itr(hbv.begin(kmerId));
          mHashF = this->hashFront(itr);
          mHashR = this->hashRear(itr+K);
          return mHashF ^ mHashR; }

        void clear() { mpHBV = 0; }

        size_t stepF()
        { HugeBVec::const_iterator itr(mpHBV->begin(mKmerId++));
          mHashF = this->stepFrontFwd(mHashF,itr[0],itr[EFRONT]);
          mHashR = this->stepRearFwd(mHashR,itr[SREAR],itr[K]);
          return mHashF ^ mHashR; }

        size_t stepF( unsigned char predecessor,
                        unsigned char successor )
        { HugeBVec::const_iterator itr(mpHBV->begin(mKmerId++));
          mHashF = this->stepFrontFwd(mHashF,predecessor,itr[EFRONT]);
          mHashR = this->stepRearFwd(mHashR,itr[SREAR],successor);
          return mHashF ^ mHashR; }

        size_t stepR()
        { HugeBVec::const_iterator itr(mpHBV->begin(--mKmerId));
          mHashF = this->stepFrontRev(mHashF,itr[0],itr[EFRONT]);
          mHashR = this->stepRearRev(mHashR,itr[SREAR],itr[K]);
          return mHashF ^ mHashR; }

        size_t stepR( unsigned char predecessor,
                        unsigned char successor )
        { HugeBVec::const_iterator itr(mpHBV->begin(--mKmerId));
          mHashF = this->stepFrontRev(mHashF,predecessor,itr[EFRONT]);
          mHashR = this->stepRearRev(mHashR,itr[SREAR],successor);
          return mHashF ^ mHashR; }

    private:
        HugeBVec const* mpHBV;
        size_t mKmerId;
        size_t mHashF;
        size_t mHashR;
    };
    friend class IncrementalHasher;

    // compares the base codes of two kmers for equality
    class Comparator : public std::binary_function<Entry,Entry,bool>
    {
    public:
        Comparator( KmerDict const& dict ) : mpDict(&dict) {}

        // compiler-supplied copying and destructor are OK

        bool operator()( Entry const& e1, Entry const& e2 ) const
        { if ( Entry::equalKeys(e1,e2) ) return true;
          typedef HugeBVec::const_iterator Itr;
          Itr itr(mpDict->getBV(e1).begin(e1.getKmerID()));
          Itr end(itr+K);
          HugeBVec const& hbv2 = mpDict->getBV(e2);
          size_t kmerId2 = e2.getKmerID();
          using std::equal;
          return e1.getForm() == e2.getForm() ?
                  equal(itr,end,hbv2.begin(kmerId2)) :
                  equal(itr,end,hbv2.rcbegin(kmerId2+K)); }

    private:
        KmerDict const* mpDict;
    };
    friend class Comparator;

    class KeyComparator : public std::binary_function<Entry,Entry,bool>
    {
    public:
        bool operator()( Entry const& e1, Entry const& e2 ) const
        { return Entry::equalKeys(e1,e2); }
    };

    // comparator for K-1 bases at the beginning of the kmer,
    // useful in getting predecessors and successors
    class HeadComparator : public std::binary_function<Entry,Entry,bool>
    {
    public:
        HeadComparator( KmerDict const& dict, unsigned char succCode )
        : mpDict(&dict), mSuccCode(succCode) {}

        // compiler-supplied copying and destructor are OK

        bool operator()( Entry const& e1, Entry const& e2 ) const
        { if ( Entry::equalKeys(e1,e2) ) return true;
          HugeBVec::const_iterator i1(mpDict->getBV(e1).begin(e1.getKmerID()));
          HugeBVec const& hbv2 = mpDict->getBV(e2);
          size_t k2 = e2.getKmerID();
          bool result;
          using std::equal;
          if ( e1.getForm() == e2.getForm() )
          { HugeBVec::const_iterator i2(hbv2.begin(k2));
            HugeBVec::const_iterator end(i2+K-1);
            result = mSuccCode == *end && equal(i2,end,i1); }
          else
          { HugeBVec::const_rc_iterator i2(hbv2.rcbegin(k2+K));
            HugeBVec::const_rc_iterator end(i2+K-1);
            result = mSuccCode == *end && equal(i2,end,i1); }
          return result; }

    private:
        KmerDict const* mpDict;
        unsigned char mSuccCode;
    };
    friend class HeadComparator;

    // comparator for K-1 bases at the end of the kmer,
    // useful in getting predecessors and successors
    class TailComparator : public std::binary_function<Entry,Entry,bool>
    {
    public:
        TailComparator( KmerDict const& dict, unsigned char predCode )
        : mpDict(&dict), mPredCode(predCode) {}

        // compiler-supplied copying and destructor are OK

        bool operator()( Entry const& e1, Entry const& e2 ) const
        { if ( Entry::equalKeys(e1,e2) ) return true;
          typedef HugeBVec::const_iterator Itr;
          Itr itr(mpDict->getBV(e1).begin(e1.getKmerID()+1));
          Itr end(itr+K-1);
          HugeBVec const& hbv2 = mpDict->getBV(e2);
          size_t k2 = e2.getKmerID();
          bool result;
          using std::equal;
          if ( e1.getForm() == e2.getForm() )
          { Itr i2(hbv2.begin(k2));
            result = mPredCode == *i2 && equal(itr,end,++i2); }
          else
          { HugeBVec::const_rc_iterator i2(hbv2.rcbegin(k2+K));
            result = mPredCode == *i2 && equal(itr,end,++i2); }
          return result; }

    private:
        KmerDict const* mpDict;
        unsigned char mPredCode;
    };
    friend class TailComparator;

    // parallel processor for read addition
    class Processor
    {
    public:
        Processor( KmerDict& dict ) : mpDict(&dict) {}
        // compiler-supplied copying and destructor are OK

        void operator()( std::pair<size_t,size_t> const& workItem ) const
        { mpDict->processRead(workItem.first,workItem.second); }

    private:
        KmerDict* mpDict;
    };
    friend class Processor;

    typedef HashSet<Entry,Hasher,Comparator> Set;
    typedef typename Set::const_iterator OCItr;
    typedef typename Set::ICItr ICItr;

    KmerDict()
    : mKSet(0,Hasher(*this),Comparator(*this)) {}

    KmerDict( size_t readBases, unsigned coverage )
    : mKSet(5*readBases/4/coverage,Hasher(*this),Comparator(*this))
    { mReadBV.reserve(readBases+1).push_back(0);
      mNewBV.reserve(readBases/coverage).push_back(0); }

    // compiler-supplied destructor is OK

    OCItr begin() const { return mKSet.begin(); }
    OCItr end() const { return mKSet.end(); }

    size_t size() const { return mKSet.size(); }

    template <class Itr> // Itr is an iterator over bvecs
    void addReads( Itr itr, Itr const& end, unsigned nThreads )
    { typedef std::pair<size_t,size_t> Workitem;
      Processor proc(*this);
      Worklist<Workitem,Processor> worklist(proc,nThreads);
      while ( itr != end )
      { bvec const& read = *itr; ++itr;
        size_t readLen = read.size();
        if ( readLen < K )
          LRP_LOG(5,"Ignoring read of length " << readLen << '.');
        else
        { LRP_LOG(5,"Adding read of length " << readLen << '.');
          size_t offset = mReadBV.size();
          ForceAssertLe(offset+readLen,mReadBV.capacity());
          mReadBV.append(read.begin(),read.end());
          worklist.add(Workitem(offset,readLen)); } }
      worklist.waitForDone();
      AssertEq(mKSet.validateBinAssignments(),0ul); }

    // NOT thread-safe!
    template <class Itr>
    Entry const& lookup( Itr const& itr, bool* pRC )
    { size_t len = mReadBV.size();
      mReadBV.append(itr,itr+K);
      Entry tmp(len,CF<K>::getForm(itr));
      Entry const* pEntry = mKSet.lookup(tmp);
      ForceAssert(pEntry);
      *pRC = tmp.getForm() != pEntry->getForm();
      mReadBV.resize(len);
      return *pEntry; }

    Entry const& lookup( size_t kmerId, size_t edgeId ) const
    { Entry tmp(kmerId,edgeId);
      size_t hash = Hasher(*this)(tmp);
      Entry const* pEntry = mKSet.lookup(hash,tmp,KeyComparator());
      ForceAssert(pEntry);
      return *pEntry; }

    size_t assign( Entry const& entry, size_t edgeId, bool rc, bool extend )
    { AssertNot(entry.isAssigned());
      size_t kmerId = entry.getKmerID();
      if ( rc )
      { if ( extend )
          mNewBV.push_back(mReadBV[kmerId] ^ 3);
        else
        { HugeBVec::const_rc_iterator itr(mReadBV.rcbegin(kmerId+K));
          mNewBV.append(itr,itr+K); } }
      else
      { if ( extend )
          mNewBV.push_back(mReadBV[kmerId+K-1]);
        else
        { HugeBVec::const_iterator itr(mReadBV.begin(kmerId));
          mNewBV.append(itr,itr+K); } }
      CanonicalForm form = entry.getForm();
      size_t newKmerId = mNewBV.size()-K;
      const_cast<KmerDictEntry&>(entry).assign(newKmerId,edgeId,rc);
      AssertEq(mKSet.lookup(Entry(kmerId,form)),&entry);
      return newKmerId; }

    Entry const& lookupPredecessor( IncrementalHasher& hasher,
                                    Entry const& entry, unsigned base0,
                                    bool* pRC )
    { typedef HugeBVec::const_iterator Itr;
      Itr itr(getBV(entry).begin(entry.getKmerID()));
      CanonicalForm form;
      if ( K&1 ) form = itr[K/2-1]<2 ? CanonicalForm::FWD : CanonicalForm::REV;
      else
      { unsigned char last = itr[K-2]^3;
        if ( base0 < last ) form = CanonicalForm::FWD;
        else if ( base0 > last ) form = CanonicalForm::REV;
        else form = CF<K-2>::getForm(itr); }
      Entry pred(entry); pred.toPrev(form);
      unsigned char baseN = itr[K-1];
      hasher.hash(entry);
      size_t hash = hasher.stepR(base0,baseN);
      Entry const* pEntry = mKSet.lookup(hash,pred,TailComparator(*this,base0));
      if ( !pEntry )
          FatalErr("Internal error looking up predecessor "
                      << Base::val2Char(base0) << " of entry "
                      << Entry::Printer(entry,K,getBV(entry)));
      hasher.stepF(base0,baseN);
      *pRC = pred.getForm() != pEntry->getForm();
      return *pEntry; }

    Entry const& lookupSuccessor( IncrementalHasher& hasher,
                                  Entry const& entry, unsigned baseN,
                                  bool* pRC )
    { typedef HugeBVec::const_iterator Itr;
      Itr itr(getBV(entry).begin(entry.getKmerID()));
      CanonicalForm form;
      if ( K&1 ) form = itr[K/2+1]<2 ? CanonicalForm::FWD : CanonicalForm::REV;
      else
      { unsigned char first = itr[1]; unsigned char last = baseN ^ 3;
        if ( first < last ) form = CanonicalForm::FWD;
        else if ( first > last ) form = CanonicalForm::REV;
        else form = CF<K-2>::getForm(itr+2); }
      Entry succ(entry); succ.toNext(form);
      unsigned char base0 = *itr;
      hasher.hash(entry);
      size_t hash = hasher.stepF(base0,baseN);
      Entry const* pEntry = mKSet.lookup(hash,succ,HeadComparator(*this,baseN));
      if ( !pEntry )
          FatalErr("Internal error looking up successor "
                      << Base::val2Char(baseN) << " of entry "
                      << Entry::Printer(entry,K,getBV(entry)));
      hasher.stepR(base0,baseN);
      *pRC = succ.getForm() != pEntry->getForm();
      return *pEntry; }

    String describe( Entry const& entry, bool rc = false ) const
    { size_t kmerId = entry.getKmerID();
      String result; result.reserve(80);
      result = ToString(kmerId); result += ':';
      HugeBVec const& hbv = getBV(entry);
      if ( rc )
        result += describeKmer(hbv.rcbegin(kmerId+K));
      else
        result += describeKmer(hbv.begin(kmerId));
      return result; }

    HugeBVec const& getSeq() const { return mNewBV; }

    void assignmentDone()
    { HugeBVec tmp; tmp.swap(mReadBV); }

    void readBinary( BinaryReader& br )
    { br.read(&mNewBV); br.read(&mKSet); br.read(&mReadBV); }

    void writeBinary( BinaryWriter& bw ) const
    { bw.write(mNewBV); bw.write(mKSet); bw.write(mReadBV); }

    static size_t externalSizeof() { return 0; }

private:
    KmerDict( KmerDict const& ); // unimplemented -- no copying
    KmerDict& operator=( KmerDict const& ); // unimplemented -- no copying

    HugeBVec const& getBV( Entry const& entry ) const
    { return entry.isAssigned() ? mNewBV : mReadBV; }

    void processRead( size_t offset, size_t len )
    { if ( len < K ) return;  // EARLY RETURN!
      auto itr(mReadBV.begin(offset));
      IncrementalHasher hasher(*this);
      if ( !(len -= K) )
      { CanonicalForm form = CF<K>::getForm(itr);
        Entry entry(offset,form);
        mKSet.apply(hasher(entry),entry,Entry::Functor(entry)); }
      else
      { KMerContext context = KMerContext::initialContext(itr[K]);
        CanonicalForm form = CF<K>::getForm(itr);
        Entry entry0(offset,form,context);
        mKSet.apply(hasher(entry0),entry0,Entry::Functor(entry0));
        while ( --len )
        { context = KMerContext(*itr,itr[K+1]);
          form = CF<K>::getForm(++itr);
          Entry entry(++offset,form,context);
          mKSet.apply(hasher(entry),entry,Entry::Functor(entry)); }
        context = KMerContext::finalContext(*itr);
        form = CF<K>::getForm(++itr);
        Entry entryN(++offset,form,context);
        mKSet.apply(hasher(entryN),entryN,Entry::Functor(entryN)); } }

    template <class Itr>
    static String describeKmer( Itr start )
    { String result; result.reserve(28);
      for ( Itr itr(start), end(start+10); itr != end; ++itr )
        result.push_back(Base::val2Char(*itr));
      result.append("...");
      if ( K&1 )
      { result.push_back(Base::val2Char(start[K/2]));
        result.append("..."); }
      for ( Itr itr(start+K-10), end(start+K); itr != end; ++itr )
        result.push_back(Base::val2Char(*itr));
      return result; }

    HugeBVec mReadBV;
    HugeBVec mNewBV;
    Set mKSet;
};

template <unsigned K> size_t const KmerDict<K>::Hasher::gBits0[] =
{
    0x4a2b3ce2251bb6e3,
    0xbcd60b1e3af4a91d,
    0x6381f301d5e0512c,
    0x957cc4fdca0f4ed2
};

template <unsigned K> size_t const KmerDict<K>::Hasher::gBitsK[] =
{
     (gBits0[0] << (HALFK&0x3F)) | (gBits0[0] >> (-HALFK&0x3F)),
     (gBits0[1] << (HALFK&0x3F)) | (gBits0[1] >> (-HALFK&0x3F)),
     (gBits0[2] << (HALFK&0x3F)) | (gBits0[2] >> (-HALFK&0x3F)),
     (gBits0[3] << (HALFK&0x3F)) | (gBits0[3] >> (-HALFK&0x3F))
};

class GraphFileNamer
{
public:
    GraphFileNamer( String const& outDir, String const& fastb, unsigned K )
    : mBase(Directory(outDir).file(File(fastb).removeExtension().filename())) {}

    // compiler-supplied copying and destructor are OK

    std::string getInfoFilename() const
    { return mBase.addExtension(".info").toString(); }
    std::string getEdgesFilename() const
    { return mBase.addExtension(".edges").toString(); }
    std::string getDictFilename() const
    { return mBase.addExtension(".dict").toString(); }
    std::string getFASTGFilename() const
    { return mBase.addExtension(".fastg").toString(); }

private:
    File mBase;
};

typedef ID<4> ComponentID;

class GraphEdge
{
public:
    GraphEdge() : mRCFlags(0), mIsPalindrome(false) {}
    GraphEdge( size_t initialKmerID, size_t componentID, bool isPalindrome )
    : mInitialKmerID(initialKmerID), mFinalKmerID(initialKmerID),
      mComponentID(componentID), mRCFlags(0), mIsPalindrome(isPalindrome) {}

    size_t getInitialKmerId() const { return mInitialKmerID.val(); }
    size_t getFinalKmerId() const { return mFinalKmerID.val(); }

    size_t getLength( unsigned K = 1 ) const
    { return mFinalKmerID.val()-mInitialKmerID.val()+K; }

    size_t getComponentID() const { return mComponentID.val(); }

    void extend() { mFinalKmerID.setVal(mFinalKmerID.val()+1); }

    bool hasPredecessor( unsigned baseCode ) const
    { AssertLt(baseCode,4u); return !mConnections[baseCode].isNull(); }
    size_t getPredecessorId( unsigned baseCode ) const
    { AssertLt(baseCode,4u); return mConnections[baseCode].val(); }
    bool isPredecessorRC( unsigned baseCode ) const
    { AssertLt(baseCode,4u); return mRCFlags & (1 << baseCode); }
    void setPredecessor( unsigned baseCode, size_t edgeId, bool rc )
    { AssertLt(baseCode,4u); setConnection(baseCode,edgeId,rc); }

    bool hasSuccessor( unsigned baseCode ) const
    { AssertLt(baseCode,4u); return !mConnections[baseCode+4u].isNull(); }
    size_t getSuccessorId( unsigned baseCode ) const
    { AssertLt(baseCode,4u); return mConnections[baseCode+4u].val(); }
    bool isSuccessorRC( unsigned baseCode ) const
    { AssertLt(baseCode,4u); return mRCFlags & (1 << (baseCode+4u)); }
    void setSuccessor( unsigned baseCode, size_t edgeId, bool rc )
    { AssertLt(baseCode,4u); setConnection(baseCode+4u,edgeId,rc); }

    template <class OItr>
    void getAllConnections( OItr oItr ) const
    { for ( EdgeID const& edgeID : mConnections )
        if ( !edgeID.isNull() ) *oItr++ = edgeID.val(); }

    bool isPalindrome() const { return mIsPalindrome; }

private:
    void setConnection( unsigned baseCode, size_t edgeId, bool rc )
    { if ( !mConnections[baseCode].isNull() )
      { AssertEq(mConnections[baseCode].val(),edgeId); }
      else
      { mConnections[baseCode].setVal(edgeId);
        if ( rc ) mRCFlags |= (1u << baseCode); } }

    KmerID mInitialKmerID;
    KmerID mFinalKmerID;
    EdgeID mConnections[8];
    ComponentID mComponentID;
    unsigned char mRCFlags;
    bool mIsPalindrome;
};


// Data about the graph files stashed in the .info file.
struct GraphInfo
{
    GraphInfo()
    : mVersion(CURRENT_VERSION), nComponents(0), nKmers(0), nEdges(0),
      kSeqSize(0), K(0) {}
    GraphInfo( size_t components, size_t kmers, size_t edges, size_t seqSize,
                    unsigned k )
    : mVersion(CURRENT_VERSION), nComponents(components), nKmers(kmers),
      nEdges(edges), kSeqSize(seqSize), K(k) {}

    // compiler-supplied copying and destructor are OK

    size_t mVersion;
    size_t nComponents;
    size_t nKmers;
    size_t nEdges;
    size_t kSeqSize;
    unsigned K;
    static size_t const CURRENT_VERSION = 1;
};

// creates the files that represent a unipath graph
template <unsigned K>
class GraphBuilder
{
public:
    typedef std::vector<GraphEdge> EdgeVec;

    GraphBuilder( KmerDict<K>& dict, EdgeVec& edges )
    : mDict(dict), mEdges(edges), mNComponents(0)
    {}

    // compiler-supplied destructor is OK

    void buildUnipaths();
    void joinGraph();
    void dumpGraph( GraphFileNamer const& namer ) const;
    void writeFASTG( GraphFileNamer const& namer ) const;
    void logStats() const
    { LRP_LOG(1,"There are " << mEdges.size() << " unipaths in "
              << mNComponents
              << " connected components\nThe N50 unipath length is "
              << getN50EdgeLen() << " kmers"); }

    size_t getNComponents() const { return mNComponents; }

private:
    GraphBuilder( GraphBuilder const& ); // unimplemented -- no copying
    GraphBuilder& operator=( GraphBuilder const& ); // unimplemented--no copying

    struct EdgeStart
    {
        EdgeStart( KmerDictEntry const& entry, bool rc, bool extend )
        : mpEntry(&entry), mRC(rc), mExtend(extend) {}

        KmerDictEntry const* mpEntry;
        bool mRC;
        bool mExtend;
    };
    typedef std::vector<EdgeStart> EdgeStartVec;

    void buildEdges( EdgeStartVec& toBuild );
    void dumpFASTA( std::ostream& out, GraphEdge const& edge ) const;

    // N50 edge length in kmers
    size_t getN50EdgeLen() const
    { vec<size_t> lengths;
      size_t nnn = mEdges.size();
      lengths.reserve(nnn);
      typedef EdgeVec::const_iterator Itr;
      for ( Itr itr(mEdges.begin()), end(mEdges.end()); itr != end; ++itr )
        lengths.push_back(itr->getLength());
      std::sort(lengths.begin(),lengths.end());
      return N50(lengths); }

    KmerDict<K>& mDict;
    EdgeVec& mEdges;
    size_t mNComponents;
};

template <unsigned K>
void GraphBuilder<K>::buildUnipaths()
{
    mEdges.reserve(100000);

    typedef typename KmerDict<K>::OCItr OCItr;
    typedef typename KmerDict<K>::ICItr ICItr;
    EdgeStartVec toBuild;

    OCItr oEnd(mDict.end());
    for ( OCItr oItr(mDict.begin()); oItr != oEnd; ++oItr )
    {
        for ( ICItr itr(oItr->begin()), end(oItr->end()); itr != end; ++itr )
        {
            KmerDictEntry const& entry = *itr;
            if ( !entry.isAssigned() )
            {
                if ( entry.getContext().getPredecessorCount() != 1 )
                {
                    LRP_LOG(4,"IF: " << mDict.describe(entry));
                    toBuild.push_back(EdgeStart(entry,false,false));
                }
                else if ( entry.getContext().getSuccessorCount() != 1 )
                {
                    LRP_LOG(4,"IR: " << mDict.describe(entry,true));
                    toBuild.push_back(EdgeStart(entry,true,false));
                }
                if ( !toBuild.empty() )
                {
                    buildEdges(toBuild);
                    mNComponents += 1;
                }
            }
        }
    }

    // leaves only smooth rings, but we'd better scan for that.  you never know.
    for ( OCItr oItr(mDict.begin()); oItr != oEnd; ++oItr )
    {
        for ( ICItr itr(oItr->begin()), end(oItr->end()); itr != end; ++itr )
        {
            KmerDictEntry const& entry = *itr;
            if ( !entry.isAssigned() )
            {
                LRP_LOG(4,"I1: " << mDict.describe(entry));
                toBuild.push_back(EdgeStart(entry,false,false));
                buildEdges(toBuild);
                mNComponents += 1;
            }
        }
    }
    mDict.assignmentDone();
}

template <unsigned K>
void GraphBuilder<K>::joinGraph()
{
    typename KmerDict<K>::IncrementalHasher hasher(mDict);

    for ( size_t idx = 0; idx < mEdges.size(); ++idx )
    {
        GraphEdge& edge = mEdges[idx];
        KmerDictEntry const& entry1 = mDict.lookup(edge.getInitialKmerId(),idx);
        KMerContext context = entry1.getContext();
        if ( context.getPredecessorCount() )
        {
            for ( unsigned predCode = 0; predCode < 4u; ++predCode )
            {
                if ( context.isPredecessor(predCode) )
                {
                    bool rc = false;
                    KmerDictEntry const& entry =
                            mDict.lookupPredecessor(hasher,entry1,predCode,&rc);
                    size_t edgeId = entry.getEdgeID();
                    AssertEq(entry.getKmerID(),
                              rc?mEdges[edgeId].getInitialKmerId():
                                 mEdges[edgeId].getFinalKmerId());
                    edge.setPredecessor(predCode,edgeId,rc);
                    if ( gLogLevel >= 4 )
                    {
                        if ( rc ) std::cout << '~';
                        std::cout << edgeId << "<-"
                                    << Base::val2Char(predCode) << ' ';
                    }
                }
            }
        }

        if ( gLogLevel >= 4 )
            std::cout << idx << '(' << edge.getLength() << ')';

        KmerDictEntry const& entryN = mDict.lookup(edge.getFinalKmerId(),idx);
        context = entryN.getContext();
        if ( context.getSuccessorCount() )
        {
            for ( unsigned succCode = 0; succCode < 4u; ++succCode )
            {
                if ( context.isSuccessor(succCode) )
                {
                    bool rc = false;
                    KmerDictEntry const& entry =
                            mDict.lookupSuccessor(hasher,entryN,succCode,&rc);
                    size_t edgeId = entry.getEdgeID();
                    AssertEq(entry.getKmerID(),
                              !rc?mEdges[edgeId].getInitialKmerId():
                                  mEdges[edgeId].getFinalKmerId());
                    edge.setSuccessor(succCode,edgeId,rc);
                    if ( gLogLevel >= 4 )
                    {
                        std::cout << ' ' << Base::val2Char(succCode) << "->";
                        if ( rc ) std::cout << '~';
                        std::cout << edgeId;
                    }
                }
            }
        }

        if ( gLogLevel >= 4 )
            std::cout << std::endl;
    }
}

template <unsigned K>
void GraphBuilder<K>::dumpGraph( GraphFileNamer const& namer ) const
{
    unlink(namer.getInfoFilename().c_str());
    BinaryWriter::writeFile(namer.getDictFilename(),mDict);
    BinaryWriter::writeFile(namer.getEdgesFilename(),mEdges);

    size_t seqSize = mDict.getSeq().size();
    GraphInfo gi(mNComponents,mDict.size(),mEdges.size(),seqSize,K);
    BinaryWriter::writeFile(namer.getInfoFilename(),gi);
}

template <unsigned K>
void GraphBuilder<K>::writeFASTG( GraphFileNamer const& namer ) const
{
    std::ofstream out(namer.getFASTGFilename().c_str());
    out << "#begin FASTG version0.62\n";
    size_t connId = 0;
    size_t nnn = mEdges.size();
    for ( size_t idx = 0; idx != nnn; ++idx )
    {
        size_t cID = mEdges[idx].getComponentID();
        size_t end = idx + 1;
        while ( end != nnn && mEdges[end].getComponentID() == cID )
            end += 1;
        if ( end-idx == 1 )
        {
            out << "\n>" << "edge" << idx << '\n';
            dumpFASTA(out,mEdges[idx]);
        }
        else
        {
            out << '\n';
            for ( size_t idx2 = idx; idx2 != end; ++idx2 )
            {
                size_t node = 2*idx2;
                out << '>' << "edge" << idx2 << ':'
                        << node << '-' << node+1 << ";\n";
                dumpFASTA(out,mEdges[idx2]);
            }
            for ( size_t idx2 = idx; idx2 != end; ++idx2 )
            {
                GraphEdge const& edge = mEdges[idx2];
                for ( unsigned base = 0; base < 4u; ++base )
                {
                    if ( edge.hasPredecessor(base) )
                    {
                        size_t edgeId = edge.getPredecessorId(base);
                        if ( edgeId >= idx2 )
                        {
                            out << ">c" << connId++ << ':';
                            if ( edge.isPredecessorRC(base) )
                                out << 2*idx2 << '-' << 2*edgeId;
                            else
                                out << 2*edgeId+1 << '-' << 2*idx2;
                            out << ":connection;\n";
                        }
                    }
                    if ( edge.hasSuccessor(base) )
                    {
                        size_t edgeId = edge.getSuccessorId(base);
                        if ( edgeId >= idx2 )
                        {
                            out << ">c" << connId++ << ':';
                            if ( edge.isSuccessorRC(base) )
                                out << 2*idx2+1 << '-' << 2*edgeId+1;
                            else
                                out << 2*idx2+1 << '-' << 2*edgeId;
                            out << ":connection;\n";
                        }
                    }
                }
            }
        }
        idx = end - 1;
    }
    out << "\n#end FASTG\n";
    out.close();
}

template <unsigned K>
void GraphBuilder<K>::dumpFASTA( std::ostream& out,
                                        GraphEdge const& edge ) const
{
    std::ostream_iterator<char> outItr(out);
    HugeBVec const& hbv = mDict.getSeq();
    typedef HugeBVec::const_iterator Itr;
    Itr itr(hbv.begin(edge.getInitialKmerId()));
    Itr end(hbv.begin(edge.getFinalKmerId()+1));
    while ( end-itr >= 80 )
    {
        std::transform(itr,itr+80,outItr,&Base::val2Char);
        out << '\n';
        itr += 80;
    }
    if ( itr != end )
    {
        std::transform(itr,end,outItr,&Base::val2Char);
        out << '\n';
    }
}

template <unsigned K>
void GraphBuilder<K>::buildEdges( EdgeStartVec& toBuild )
{
    LRP_LOG(2,"Building component " << mNComponents);

    typename KmerDict<K>::IncrementalHasher hasher(mDict);
    while ( !toBuild.empty() )
    {
        EdgeStart const& es = toBuild.back();
        KmerDictEntry const* pEntry = es.mpEntry;
        if ( pEntry->isAssigned() )
        {
            toBuild.pop_back();
            continue; // Non-structured transfer!
        }

        size_t edgeId = mEdges.size();
        size_t newKmerId = mDict.assign(*pEntry,edgeId,es.mRC,es.mExtend);
        toBuild.pop_back();

        LRP_LOG(3,"EE: " << mDict.describe(*pEntry) << ' ' << edgeId);
        mEdges.push_back(GraphEdge(newKmerId,mNComponents,
                                    pEntry->isPalindrome()));
        GraphEdge& edge = mEdges.back();
        KMerContext context = pEntry->getContext();

        // queue predecessors
        if ( context.getPredecessorCount() )
        {
            for ( unsigned predCode = 0; predCode < 4u; ++predCode )
            {
                if ( context.isPredecessor(predCode) )
                {
                    bool rc = false;
                    KmerDictEntry const& entry =
                           mDict.lookupPredecessor(hasher,*pEntry,predCode,&rc);
                    if ( !entry.isAssigned() )
                    {
                        LRP_LOG(4,"QP: " << mDict.describe(entry,!rc)
                                    << ' ' << edgeId);
                        toBuild.push_back(EdgeStart(entry,!rc,false));
                    }
                }
            }
        }

        // while we're on the straight and narrow, with just a single successor,
        // try to extend the edge
        unsigned nSuccessors;
        while ( (nSuccessors = context.getSuccessorCount()) == 1 )
        {
            // if we're a palindrome we must queue successors
            if ( pEntry->isPalindrome() )
                break;

            nSuccessors = 0; // we'll handle our single successor right now

            unsigned succCode = 0;
            while ( !context.isSuccessor(succCode) )
                succCode += 1;
            bool rc = false;
            KmerDictEntry const& entry =
                    mDict.lookupSuccessor(hasher,*pEntry,succCode,&rc);

            // if single successor is already built, we're done extending
            if ( entry.isAssigned() )
                break;

            // if successor has multiple predecessors, or is a palindrome
            // start a new edge
            KMerContext succContext = entry.getContext();
            if ( rc ) succContext = succContext.rc();
            if ( succContext.getPredecessorCount() > 1 || entry.isPalindrome() )
            {
                LRP_LOG(4,"QS: " << mDict.describe(entry,rc) << ' '
                        << edgeId);
                toBuild.push_back(EdgeStart(entry,rc,true));
                break;
            }

            // no branching -- extend the edge we're building
            mDict.assign(entry,edgeId,rc,true);
            LRP_LOG(5,"XE: " << mDict.describe(entry) << ' ' << edgeId);
            edge.extend();

            pEntry = &entry;
            context = entry.getContext();
        }

        // if there are unhandled successors, queue them up
        if ( nSuccessors )
        {
            for ( unsigned succCode = 0; succCode < 4u; ++succCode )
            {
                if ( context.isSuccessor(succCode) )
                {
                    bool rc = false;
                    KmerDictEntry const& entry =
                            mDict.lookupSuccessor(hasher,*pEntry,succCode,&rc);
                    if ( !entry.isAssigned() )
                    {
                        LRP_LOG(4,"QS: " << mDict.describe(entry,rc)
                                    << ' ' << edgeId);
                        toBuild.push_back(EdgeStart(entry,rc,!--nSuccessors));
                        // the last thing pushed (which is the first to be
                        // processed next time around) can be an extension
                    }
                }
            }
        }
        LRP_LOG(3,"Edge " << edgeId << " has length " << edge.getLength());
    }
}


template <unsigned K>
class Graph
{
public:
    typedef typename GraphBuilder<K>::EdgeVec EdgeVec;
    typedef typename EdgeVec::const_iterator const_iterator;

    Graph( GraphFileNamer const& namer )
    {
        GraphInfo gi;
        BinaryReader::readFile(namer.getInfoFilename(),&gi);
        ForceAssertEq(gi.mVersion,GraphInfo::CURRENT_VERSION);
        ForceAssertEq(gi.K,K);
        BinaryReader::readFile(namer.getEdgesFilename(),&mEdges);
        ForceAssertEq(gi.nEdges,mEdges.size());
        BinaryReader br(namer.getDictFilename());
        br.read(&mBases);
        ForceAssertEq(gi.kSeqSize,mBases.size());
        mNComponents = gi.nComponents;
    }

    Graph( KmerDict<K>& dict )
    {
        GraphBuilder<K> gb(dict,mEdges);
        gb.buildUnipaths();
        gb.joinGraph();
        mBases = dict.getSeq();
        mNComponents = gb.getNComponents();
    }

    // compiler-supplied copying and destructor are OK

    size_t getNEdges() const { return mEdges.size(); }
    size_t getNComponents() const { return mNComponents; }

    GraphEdge const& getEdge( size_t edgeId ) const
    { AssertLt(edgeId,mEdges.size()); return mEdges[edgeId]; }

    size_t nextEdgeId( size_t edgeId, unsigned char baseCode, bool* pNC )
    { GraphEdge const& edge = getEdge(edgeId);
      if ( *pNC )
      { baseCode ^= 3;
        ForceAssert(edge.hasPredecessor(baseCode));
        edgeId = edge.getPredecessorId(baseCode);
        *pNC = !edge.isPredecessorRC(baseCode); }
      else
      { ForceAssert(edge.hasSuccessor(baseCode));
        edgeId = edge.getSuccessorId(baseCode);
        *pNC = edge.isSuccessorRC(baseCode); }
      return edgeId; }

    const_iterator begin() const { return mEdges.begin(); }
    const_iterator end() const { return mEdges.end(); }

    HugeBVec::const_iterator basesBegin( size_t edgeId ) const
    { return mBases.begin(getEdge(edgeId).getInitialKmerId()); }

    HugeBVec::const_iterator basesEnd( size_t edgeId ) const
    { return mBases.begin(getEdge(edgeId).getFinalKmerId()+K); }

    void validate( KmerDict<K> const& dict );

    friend ostream& operator<<( ostream& os, Graph const& graph )
    { size_t nEdges = graph.getNEdges();
      for ( size_t edgeId = 0; edgeId != nEdges; ++edgeId )
      { GraphEdge const& edge = graph.getEdge(edgeId);
        for ( unsigned predCode = 0; predCode != 4u; ++predCode )
        { if ( edge.hasPredecessor(predCode) )
          { if ( edge.isPredecessorRC(predCode) ) os << '~';
            os << edge.getPredecessorId(predCode) << "<-["
                    << Base::val2Char(predCode) << "] "; } }
        if ( edge.isPalindrome() ) os << '|' << edgeId << '|';
        else os << edgeId;
        os << '.' << edge.getComponentID();
        os << '(' << edge.getLength() << ')';
        for ( unsigned succCode = 0; succCode != 4u; ++succCode )
        { if ( edge.hasSuccessor(succCode) )
          { os << " [" << Base::val2Char(succCode) << "]->";
            if ( edge.isSuccessorRC(succCode) ) os << '~';
            os << edge.getSuccessorId(succCode); } }
        os << '\n'; }
      return os; }

private:
    EdgeVec mEdges;
    HugeBVec mBases;
    size_t mNComponents;

};

template <unsigned K>
void Graph<K>::validate( KmerDict<K> const& dict )
{
    std::cout << Date() << " Validating unipath graph." << std::endl;
    size_t nErrors = 0;
    size_t nEdges = mEdges.size();
    for ( size_t edgeId = 0; edgeId != nEdges; ++edgeId )
    {
        GraphEdge const& edge = mEdges[edgeId];
        KMerContext context =
                dict.lookup(edge.getInitialKmerId(),edgeId).getContext();
        for ( unsigned predCode = 0; predCode != 4u; ++predCode )
        {
            if ( context.isPredecessor(predCode) )
            {
                if ( !edge.hasPredecessor(predCode) )
                {
                    std::cout << "Edge " << edgeId << " lacks predecessor "
                            << Base::val2Char(predCode) << std::endl;
                    nErrors += 1;
                }
                else
                {
                    size_t edgeId2 = edge.getPredecessorId(predCode);
                    GraphEdge const& target = mEdges[edgeId2];
                    if ( edge.getComponentID() != target.getComponentID() )
                    {
                        std::cout << "Edge " << edgeId << " has its "
                                << Base::val2Char(predCode)
                                << " predecessor in the wrong component."
                                << std::endl;
                        nErrors += 1;
                    }
                }
            }
            else if ( edge.hasPredecessor(predCode) )
            {
                std::cout << "Edge " << edgeId << " has an unexpected "
                        << Base::val2Char(predCode) << " predecessor."
                        << std::endl;
                nErrors += 1;
            }
        }
        context = dict.lookup(edge.getFinalKmerId(),edgeId).getContext();
        for ( unsigned succCode = 0; succCode != 4u; ++succCode )
        {
            if ( context.isSuccessor(succCode) )
            {
                if ( !edge.hasSuccessor(succCode) )
                {
                    std::cout << "Edge " << edgeId << " lacks successor "
                            << Base::val2Char(succCode) << std::endl;
                    nErrors += 1;
                }
                else
                {
                    size_t edgeId2 = edge.getSuccessorId(succCode);
                    GraphEdge const& target = mEdges[edgeId2];
                    if ( edge.getComponentID() != target.getComponentID() )
                    {
                        std::cout << "Edge " << edgeId << " has its "
                                << Base::val2Char(succCode)
                                << " successor in the wrong component."
                                << std::endl;
                        nErrors += 1;
                    }
                }
            }
            else if ( edge.hasSuccessor(succCode) )
            {
                std::cout << "Edge " << edgeId << " has an unexpected "
                        << Base::val2Char(succCode) << " successor."
                        << std::endl;
                nErrors += 1;
            }
        }
    }
    ForceAssertEq(nErrors,0ul);
}

} // end of Long namespace

template <unsigned K> struct Serializability<Long::KmerDict<K> >
{ typedef SelfSerializable type; };

TRIVIALLY_SERIALIZABLE(Long::KmerDictEntry);
TRIVIALLY_SERIALIZABLE(Long::GraphEdge);
TRIVIALLY_SERIALIZABLE(Long::GraphInfo);

#endif /* KMERS_LONGREADPATHER_H_ */
