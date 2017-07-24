///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file ReadPather.h
 * \author tsharpe
 * \date May 3, 2010
 *
 * \brief Kmerization and Unipathing of a mess of sequences.
 *
 */
#ifndef KMERS_READ_PATHER_H_
#define KMERS_READ_PATHER_H_

#include "Basevector.h"
#include "String.h"
#include "feudal/BinaryStream.h"
#include "feudal/HashSet.h"
#include "feudal/HugeBVec.h"
#include "feudal/MasterVec.h"
#include "feudal/SerfVec.h"
#include "feudal/VirtualMasterVec.h"
#include "kmers/KMer.h"
#include "system/Assert.h"
#include "system/ID.h"
#include "system/System.h"
#include "system/WorklistN.h"
#include "system/SpinLockedData.h"
#include <algorithm>
#include <cstddef>
#include <cstring>
#include <ostream>
#include <unistd.h>
#include <vector>

/// An ID for a KMer.
/// Kmers are assigned IDs in such a way that each edge of a unipath graph can
/// be described by a consecutive sequence of KmerIDs.
class KmerID : public ID<5>
{
public:
    KmerID() {}
    explicit KmerID( size_t id ) : ID<5>(id) {}

    // compiler-supplied copying and destructor are ok
};

TRIVIALLY_SERIALIZABLE(KmerID);

/// An ID for an edge of a unipath graph.
class EdgeID
{
public:
    EdgeID() : mID(NULLVAL) {}
    explicit EdgeID( size_t id ) { setVal(id); }

    // compiler-supplied copying and destructor are ok

    size_t val() const { return mID; }
    void setVal( size_t id ) { ForceAssertLe(id,NULLVAL); mID = id; }
    void setNull() { mID = NULLVAL; }
    bool isNull() const { return mID == NULLVAL; }

    friend bool operator==( EdgeID const& eid1, EdgeID const& eid2 )
    { return eid1.mID == eid2.mID; }

    friend bool operator!=( EdgeID const& eid1, EdgeID const& eid2 )
    { return eid1.mID != eid2.mID; }

    friend bool operator<( EdgeID const& eid1, EdgeID const& eid2 )
    { return eid1.mID < eid2.mID; }

    friend int compare( EdgeID const& eid1, EdgeID const& eid2 )
    { return compare(eid1.mID,eid2.mID); }

    friend ostream& operator<<( ostream& os, EdgeID const& eid )
    { return os << eid.mID; }

private:
    static unsigned const NULLVAL = ~0u;
    unsigned mID;
};

TRIVIALLY_SERIALIZABLE(EdgeID);

/// An ID for a connected set of edges in a unipath graph.
class ComponentID : public ID<3>
{
public:
    ComponentID() {}
    explicit ComponentID( size_t id ) : ID<3>(id) {}

    // compiler-supplied copying and destructor are ok
};

TRIVIALLY_SERIALIZABLE(ComponentID);

/// What you get if you look up a KMer that exists in the dictionary.
/// Namely, you get an EdgeID in which it occurs, and the offset of the kmer
/// in question.
class KDef
{
public:
    KDef() : mEdgeOffset(0) {}
    KDef( KMerContext context ) : mEdgeOffset(0), mContext(context) {}

    // compiler-supplied copying and destructor are ok

    EdgeID const& getEdgeID() const { return mEdgeID; }
    unsigned getEdgeOffset() const { return mEdgeOffset; }

    KMerContext getContext() const { return mContext; }
    void orContext( KMerContext context ) { mContext |= context; }
    void setContext( KMerContext context ) { mContext = context; }

    bool isNull() const { return mEdgeID.isNull(); }
    void setNull() { mEdgeID.setNull(); }
    void set( EdgeID const& edgeID, unsigned edgeOffset )
    { Assert(isNull()); ForceAssertLe(edgeOffset,MAX_OFFSET);
      mEdgeID = edgeID; mEdgeOffset = edgeOffset; }

    // the UnipathGraph stores a kmer count in mEdgeOffset during construction
    size_t getCount() const { Assert(isNull()); return mEdgeOffset; }
    void setCount( size_t count )
    { mEdgeOffset = std::min(MAX_OFFSET,count); }
    size_t incrementCount()
    { Assert(isNull()); if ( mEdgeOffset < MAX_OFFSET ) ++mEdgeOffset;
      return mEdgeOffset; }

    friend ostream& operator<<( ostream& os, KDef const& kd )
    { return os << '[' << kd.mEdgeID << '+' << kd.mEdgeOffset << "] "
            << kd.mContext; }

private:
    void checkAssumptions()
    { STATIC_ASSERT(sizeof(KDef)==sizeof(long)); }

    EdgeID mEdgeID;
    unsigned mEdgeOffset : 24;
    KMerContext mContext;
    static size_t const MAX_OFFSET = (1ul<<24)-1ul;
};

TRIVIALLY_SERIALIZABLE(KDef);

class KDefS
{
public:
    KDefS() {}
    KDefS( KMerContext context ) : mContext(context) {}

    // compiler-supplied copying and destructor are ok

    KMerContext getContext() const { return mContext; }
    void orContext( KMerContext context ) { mContext |= context; }
    void setContext( KMerContext context ) { mContext = context; }

    friend ostream& operator<<( ostream& os, KDefS const& kd )
    { return os << kd.mContext; }

private:
    KMerContext mContext;
};

struct BCWrapper {
     explicit BCWrapper(int32_t bc ) : mTempBC(bc) {};
     BCWrapper() {};
     int32_t getTempBC() const { return mTempBC; };
     void setTempBC(int32_t bc) { mTempBC = bc; };
     int32_t mTempBC;
     friend ostream& operator<<( ostream& os, BCWrapper const& ent ) { 
          if ( ent.getTempBC() > -1 ) os << "[" << ent.getTempBC() << "]";
          return os;
     };
};

struct BCEmpty {
     explicit BCEmpty(int32_t bc) {};
     BCEmpty() {};
     int32_t getTempBC() const { FatalErr("BUG: empty base class used for barcode access"); };
     void setTempBC(int32_t) {};
     friend ostream& operator<<( ostream& os, BCEmpty const& ent ) { return os; };
};

template <unsigned K, typename Bcode = BCWrapper >
class KmerDictEntry : public KMer<K>, public Bcode
{
public:
    KmerDictEntry() {}
    KmerDictEntry( KMer<K> const& kmer, KMerContext context, int32_t tempBC = -1 )
    : KMer<K>(kmer), Bcode(tempBC), mKDef(context) {}
    KmerDictEntry& operator=( KMer<K> const& kmer )
    { static_cast<KMer<K>&>(*this) = kmer; mKDef = KDef(); this->setTempBC(-1); return *this; }

    // compiler-supplied copying and destructor are OK

    KDef& getKDef() { return mKDef; }
    KDef const& getKDef() const { return mKDef; }

    friend ostream& operator<<( ostream& os, KmerDictEntry const& ent )
    {
         os << static_cast<KMer<K> const&>(ent) << " " << ent.getKDef() << " " << ent.getTempBC();
         return os;
    }

private:
    KDef     mKDef;
};

template <unsigned K>
struct Serializability< KmerDictEntry<K> >
{ typedef TriviallySerializable type; };



/// A set of KmerDictEntry's, used as a map from KMer onto KDef.
template <unsigned K, typename Bcode = BCWrapper >
class KmerDict
{
public:
    typedef KmerDictEntry<K, Bcode> Entry;
    typedef typename KMer<K>::Hasher Hasher;
    typedef std::equal_to<KMer<K> > Comparator;
    typedef HashSet<Entry,Hasher,Comparator> Set;
    typedef typename Set::const_iterator OCItr;
    typedef typename Set::iterator OItr;
    typedef typename Set::ICItr ICItr;

    KmerDict( size_t dictSize, double maxLoadFactor=.75 )
    : mKSet(dictSize,Hasher(),Comparator(),TFactory<Entry>(),maxLoadFactor) {}

    KmerDict( KmerDict const& ) = delete;
    KmerDict& operator=( KmerDict const& ) = delete;

    Entry const* findEntryCanonical( KMer<K> const& kmer ) const
    { return mKSet.lookup(kmer); }

    /// Returns null pointer if kmer isn't in dictionary.
    Entry const* findEntry( KMer<K> const& kmer ) const
    { return mKSet.lookup( kmer.getCanonicalForm() == CanonicalForm::REV ?
                            KMer<K>(kmer).rc() :
                            kmer ); }

    /// Applies functor to entry, which will be added if not present.
    template <class Func>
    void applyCanonical( KMer<K> const& kmer, Func const& func )
    { mKSet.apply(kmer,func); }

    /// Returns null pointer if kmer isn't in dictionary.
    KDef* lookup( KMer<K> const& kmer )
    { Entry const* pEnt = findEntry(kmer);
      return pEnt ? const_cast<KDef*>(&pEnt->getKDef()) : 0; }

    KDef const* lookup( KMer<K> const& kmer ) const
    { Entry const* pEnt = findEntry(kmer);
      return pEnt ? &pEnt->getKDef() : 0; }

    /// Canonicalizes and inserts kmer, if necessary.
    KDef& operator[]( KMer<K> const& kmer )
    { Entry const* pEnt;
      if ( kmer.getCanonicalForm() != CanonicalForm::REV ) pEnt = &mKSet[kmer];
      else pEnt = &mKSet[KMer<K>(kmer).rc()];
      return const_cast<KDef&>(pEnt->getKDef()); }

    /// Inserts a kmer known to be in canonical form.
    void insertCanonical( KMer<K> const& kmer )
    { mKSet.add(kmer); }

    /// Inserts a canonical kmer, known to be novel, along with its entry info
    /// directly into the dictionary without further ado.
    void insertEntry( Entry const& entry )
    { mKSet.insertUniqueValue(entry); }

    class BadKmerCountFunctor
    {
    public:
      BadKmerCountFunctor( size_t minCount, size_t maxCount = ~0ul )
      : mMinCount(minCount), mMaxCount(maxCount) {}

      bool operator()( Entry const& entry )
      { size_t count = entry.getKDef().getCount();
        return count < mMinCount || count > mMaxCount; }

    private:
      size_t mMinCount; size_t mMaxCount;
    };

    /// Retain only entries with counts in the given range.
    template <class Functor>
    void clean( Functor const& functor )
    { mKSet.remove_if(functor);
      recomputeAdjacencies(); }

    OCItr begin() const { return mKSet.begin(); }
    OCItr end() const { return mKSet.end(); }
    OCItr cbegin() { return mKSet.cbegin(); }
    OCItr cend() { return mKSet.cend(); }
    OItr begin() { return mKSet.begin(); }
    OItr end() { return mKSet.end(); }

    size_t size() const { return mKSet.size(); }

    void remove( KMer<K> const& kmer )
    { mKSet.remove( kmer.getCanonicalForm() == CanonicalForm::REV ?
                            KMer<K>(kmer).rc() : kmer ); }

    void removeCanonical( KMer<K> const& kmer ) { mKSet.remove(kmer); }

    void process( VirtualMasterVec<bvec> const& reads,
                        bool validate = false,
                        unsigned nThreads = getConfiguredNumThreads(),
                        size_t batchSize = 10000 );

    void process( vecbvec const& reads, bool verbose = false,
                        bool validate = false,
                        unsigned nThreads = getConfiguredNumThreads(),
                        size_t batchSize = 10000 );

    void clear() { mKSet.clear(); }

    friend void swap( KmerDict& kd1, KmerDict& kd2 )
    { swap(kd1.mKSet,kd2.mKSet); }

    friend bool operator==( KmerDict const& kd1, KmerDict const& kd2 )
    { return kd1.mKSet == kd2.mKSet; }

    friend bool operator!=( KmerDict const& kd1, KmerDict const& kd2 )
    { return !(kd1==kd2); }

    void writeBinary( BinaryWriter& writer ) const
    { writer.write(mKSet); }

    void readBinary( BinaryReader& reader )
    { reader.read(&mKSet); }

    static size_t externalSizeof() { return 0; }

    template <class Proc>
    void parallelForEachHHS( Proc const& proc ) const
    { mKSet.parallelForEachHHS(proc); }

    void recomputeAdjacencies()
    { parallelForEachHHS(AdjProc(*this)); }

    void nullEntries()
    { parallelForEachHHS(
                []( typename Set::HHS const& hhs )
                { for ( Entry const& entry : hhs )
                    const_cast<Entry&>(entry).getKDef().setNull(); }); }

private:
    class AdjProc
    {
    public:
        AdjProc( KmerDict& dict ) : mDict(dict) {}

        void operator()( typename Set::HHS const& hhs )
        { for ( auto itr=hhs.begin(),end=hhs.end(); itr != end; ++itr )
          { KDef& kDef = const_cast<KDef&>(itr->getKDef());
            KMerContext context = kDef.getContext();
            if ( context.getSuccessors() )
            { KMer<K> kmer(*itr);
              kmer.toSuccessor(0);
              for ( unsigned succCode = 0; succCode < 4u; ++succCode )
              { if ( context.isSuccessor(succCode) )
                { kmer.setBack(succCode);
                  if ( !mDict.findEntry(kmer) )
                    context.removeSuccessor(succCode); } } }
            if ( context.getPredecessors() )
            { KMer<K> kmer(*itr);
              kmer.toPredecessor(0);
              for ( unsigned predCode = 0; predCode < 4u; ++predCode )
              { if ( context.isPredecessor(predCode) )
                { kmer.setFront(predCode);
                  if ( !mDict.findEntry(kmer) )
                    context.removePredecessor(predCode); } } }
            kDef.setContext(context); } }

    private:
        KmerDict& mDict;
    };

    Set mKSet;
};



/// A list of KmerDictEntry's, used as a predecessor to 
/// KmerDict construction
template <unsigned K>
class KmerVec
{
public:
    typedef KmerDictEntry<K> Entry;

    KmerVec() : mNK(0) {};
    KmerVec(size_t exp_size) : mNK(0), mKVec(exp_size) {};
    KmerVec( KmerVec const& ) = delete;
    KmerVec& operator=( KmerVec const& ) = delete;

    typedef typename vec<Entry>::const_iterator OCItr;
    typedef typename vec<Entry>::iterator OItr;

    typedef typename vec<Entry>::size_type size_type;
    typedef typename vec<Entry>::value_type value_type;

#if 0
    /// Inserts a kmer known to be in canonical form.
    void insertCanonical( KMer<K> const& kmer )
    {  
         mLock.lock();
         mKVec.emplace_back(kmer); 
         mLock.unlock();
    }
#endif

    /// Inserts a canonical kmer, known to be novel, along with its entry info
    /// directly into the dictionary without further ado.
    void insertEntry( Entry const&& entry )
    { 
         mLock.lock();
#if 0
         mKVec.emplace_back( std::move(entry ));
         mLock.unlock();
#endif
         size_t newEntry = mNK++;
         if ( mNK > mKVec.size() )
              mKVec.resize(mNK);
         mLock.unlock();
         mKVec[newEntry] = std::move(entry);
    }

    OCItr begin() const { return mKVec.begin(); }
    OCItr end() const { return mKVec.end(); }
    OCItr cbegin() { return mKVec.cbegin(); }
    OCItr cend() { return mKVec.cend(); }
    OItr begin() { return mKVec.begin(); }
    OItr end() { return mKVec.end(); }

    size_t size() const { return mKVec.size(); }

    void clear() { mKVec.clear(); }

    void fit() { mKVec.resize( mNK ); }

    friend void swap( KmerVec& kd1, KmerVec& kd2 )
    { swap(kd1.mKVec,kd2.mKVec); }

    friend bool operator==( KmerVec const& kd1, KmerVec const& kd2 )
    { return kd1.mKVec == kd2.mKVec; }

    friend bool operator!=( KmerVec const& kd1, KmerVec const& kd2 )
    { return !(kd1==kd2); }

    void writeBinary( BinaryWriter& writer ) const
    { 
         if ( mKVec.size() != mNK ) 
              FatalErr("you should call fit() before writeBinary");
         writer.write(mKVec); 
    }

    void readBinary( BinaryReader& reader )
    { 
         reader.read(&mKVec); 
         mNK = mKVec.size();
    }

    static size_t externalSizeof() { return 0; }


private:
    size_t mNK;
    vec<Entry> mKVec;
    SpinLockedData mLock;
};
template <unsigned K> 
struct Serializability< KmerVec<K> >
{ typedef SelfSerializable type; };

template <unsigned K>
struct Serializability< KmerDict<K> >
{ typedef SelfSerializable type; };

/// An unbranched sequence of kmers.
/// Described by a close-ended range of kmer IDs.
/// Also includes info on predecessor and successor edges.
class UnipathEdge
{
public:
    UnipathEdge() : mRCFlags(0), mIsPalindrome(false) {}

    UnipathEdge( KmerID const& kmerID, ComponentID componentID,
                 bool isPalindrome )
    : mRCFlags(0),
      mFirstKmerID(kmerID), mLastKmerID(kmerID),
      mComponentID(componentID), mIsPalindrome(isPalindrome)
    {}

    // compiler-supplied copying and destructor are ok

    KmerID const& getInitialKmerID() const { return mFirstKmerID; }
    size_t getInitialKmerId() const { return mFirstKmerID.val(); }
    KmerID const& getFinalKmerID() const { return mLastKmerID; }
    size_t getFinalKmerId() const { return mLastKmerID.val(); }
    KmerID getKmerID( unsigned offset ) const
    { return KmerID(mFirstKmerID.val()+offset); }

    size_t getLength() const { return mLastKmerID.val()-mFirstKmerID.val()+1; }
    ComponentID const& getComponentID() const { return mComponentID; }

    size_t extend()
    { size_t newVal = mLastKmerID.val()+1;
      mLastKmerID.setVal(newVal);
      return newVal - mFirstKmerID.val(); }

    bool hasPredecessor( unsigned baseCode ) const
    { return !getPredecessor(baseCode).isNull(); }

    EdgeID const& getPredecessor( unsigned baseCode ) const
    { AssertLt(baseCode,4u); return mConnections[baseCode]; }

    size_t getPredecessorId( unsigned baseCode ) const
    { return getPredecessor(baseCode).val(); }

    bool isPredecessorRC( unsigned baseCode ) const
    { AssertLt(baseCode,4u); return (mRCFlags >> baseCode) & 1; }

    void setPredecessor( unsigned baseCode, EdgeID id, bool rc )
    { AssertLt(baseCode,4u);
      if ( !mConnections[baseCode].isNull() )
      { AssertEq(mConnections[baseCode],id); }
      else
      { mConnections[baseCode] = id;
        if ( rc ) mRCFlags |= 1 << baseCode; } }

    bool hasSuccessor( unsigned baseCode ) const
    { return !getSuccessor(baseCode).isNull(); }

    EdgeID const& getSuccessor( unsigned baseCode ) const
    { AssertLt(baseCode,4u); return mConnections[baseCode+4]; }

    size_t getSuccessorId( unsigned baseCode ) const
    { return getSuccessor(baseCode).val(); }

    bool isSuccessorRC( unsigned baseCode ) const
    { AssertLt(baseCode,4u); return (mRCFlags >> (baseCode+4)) & 1; }

    void setSuccessor( unsigned baseCode, EdgeID id, bool rc )
    { AssertLt(baseCode,4u);
      baseCode += 4u;
      if ( !mConnections[baseCode].isNull() )
      { AssertEq(mConnections[baseCode],id); }
      else
      { mConnections[baseCode] = id;
        if ( rc ) mRCFlags |= 1 << baseCode; } }

    bool isPalindrome() const
    { return mIsPalindrome; }

    bool isSource() const
    { return !hasPredecessor(0) && !hasPredecessor(1) &&
                !hasPredecessor(2) && !hasPredecessor(3); }

    bool isSink() const
    { return !hasSuccessor(0) && !hasSuccessor(1) &&
                !hasSuccessor(2) && !hasSuccessor(3); }

    void setPalindrome( bool isPalindrome ) { mIsPalindrome = isPalindrome; }
    friend bool operator==( UnipathEdge const& ue1, UnipathEdge const& ue2 )
    { return !memcmp(&ue1,&ue2,sizeof(UnipathEdge)); }

    friend bool operator!=( UnipathEdge const& ue1, UnipathEdge const& ue2 )
    { return !(ue1 == ue2); }

    friend ostream& operator<<( ostream& os, UnipathEdge const& edge )
    { os << '[' << edge.getInitialKmerID() << '-' << edge.getFinalKmerID() <<
        '@' << edge.getComponentID();
      char sep = '>';
      for ( size_t idx = 0; idx < 8; ++idx )
      { os << sep; sep = ',';
        if ( edge.mRCFlags & (1 << idx) ) os << '~';
        os << edge.mConnections[idx]; }
      return os << ']'; }

private:
    unsigned char mRCFlags;
    KmerID mFirstKmerID;
    KmerID mLastKmerID;
    ComponentID mComponentID;
    EdgeID mConnections[8];
    bool mIsPalindrome;
};

TRIVIALLY_SERIALIZABLE(UnipathEdge);
typedef std::vector<UnipathEdge> UnipathEdgeVec;

// Data about the graph files stashed in the .info file.
struct GraphInfo
{
    size_t nComponents;
    size_t nKmers;
    size_t nEdges;
    size_t kSeqSize;
    unsigned K;
};

TRIVIALLY_SERIALIZABLE(GraphInfo);

/// A bi-directional graph of UnipathEdges.
template <unsigned K>
class UnipathGraph
{
public:
    typedef KmerDict<K> KDict;

    /// Instantiate a graph from its binary files.
    UnipathGraph( String const& graphInfoFilename );

    /// Create a graph from a dictionary of kmers and hit counts.
    UnipathGraph( KDict& dict, const int verbosity = 0 );

    /// Create a graph from a mess of reads.
    UnipathGraph( VirtualMasterVec<bvec> const& reads,
                    bool validate, unsigned nThreads, size_t nKmersEstimate,
                  const int verbosity = 0 );

    /// Create a graph from a mess of reads.
    UnipathGraph( vecbvec const& reads,
                    bool validate, unsigned nThreads, size_t nKmersEstimate,
                  const int verbosity = 0 );

    ~UnipathGraph() { ungetDict(); }

    size_t getNEdges() const { return mEdges.size(); }

    UnipathEdge const& getEdge( EdgeID const& edgeID ) const
    { Assert(edgeID.val()<mEdges.size()); return mEdges[edgeID.val()]; }

    KmerID getKmerID( KDef const& kDef ) const
    { return KmerID(getEdge(kDef.getEdgeID()).getInitialKmerID().val() +
                        kDef.getEdgeOffset()); }

    UnipathEdgeVec const& getAllEdges() const { return mEdges; }

    HugeBVec::const_iterator getBases( KmerID const& kmerID ) const
    { return mSeq.begin(kmerID.val()); }
    HugeBVec::const_iterator getBases( EdgeID const& edgeID ) const
    { return mSeq.begin(getEdge(edgeID).getInitialKmerID().val()); }
    HugeBVec::const_iterator getBases( KDef const& kDef ) const
    { return mSeq.begin(getEdge(kDef.getEdgeID()).getInitialKmerID().val()
                        +kDef.getEdgeOffset()); }
    HugeBVec::const_iterator getBasesEnd( KmerID const& kmerID ) const
    { return mSeq.begin(kmerID.val()+K); }
    HugeBVec::const_iterator getBasesEnd( EdgeID const& edgeID ) const
    { return mSeq.begin(getEdge(edgeID).getFinalKmerID().val()+K); }
    HugeBVec::const_iterator getBasesEnd( KDef const& kDef ) const
    { return mSeq.begin(getEdge(kDef.getEdgeID()).getInitialKmerID().val()
                        +kDef.getEdgeOffset()+K); }

    HugeBVec::const_rc_iterator getRCBases( KmerID const& kmerID ) const
    { return mSeq.rcbegin(kmerID.val()+K); }
    HugeBVec::const_rc_iterator getRCBases( EdgeID const& edgeID ) const
    { return mSeq.rcbegin(getEdge(edgeID).getFinalKmerID().val()+K); }
    HugeBVec::const_rc_iterator getRCBases( KDef const& kDef ) const
    { return mSeq.rcbegin(getEdge(kDef.getEdgeID()).getInitialKmerID().val()
                        +kDef.getEdgeOffset()+K); }
    HugeBVec::const_rc_iterator getRCBasesEnd( KmerID const& kmerID ) const
    { return mSeq.rcbegin(kmerID.val()+K,K); }
    HugeBVec::const_rc_iterator getRCBasesEnd( EdgeID const& edgeID ) const
    { UnipathEdge const& edge = getEdge(edgeID);
      size_t len = edge.getLength()+K-1;
      return mSeq.rcbegin(edge.getFinalKmerID().val()+K,len); }
    HugeBVec::const_rc_iterator getRCBasesEnd( KDef const& kDef ) const
    { return mSeq.rcbegin(getEdge(kDef.getEdgeID()).getInitialKmerID().val()
                        +kDef.getEdgeOffset()+K,K); }

    HugeBVec const& getAllBases() const { return mSeq; }

    bool nextEdge( unsigned baseCode, EdgeID* pEdgeID, bool isRC ) const
    { UnipathEdge const& edge = getEdge(*pEdgeID);
      if ( isRC )
      { *pEdgeID = edge.getPredecessor(baseCode);
        if ( edge.isPredecessorRC(baseCode) ) isRC = false; }
      else
      { *pEdgeID = edge.getSuccessor(baseCode);
        if ( edge.isSuccessorRC(baseCode) ) isRC = true; }
      AssertNot(pEdgeID->isNull());
      return isRC; }

    // loads dictionary (which is huge) lazily
    KDict const& getDict() const
    { if ( !mpDict ) loadDict(); return *mpDict; }

    // unloads dictionary to recover lots of memory
    void ungetDict() const
    { if ( mMyDict ) { delete mpDict; mpDict = 0; } }

    void write( String const& graphInfoFilename )
    { unlink(graphInfoFilename.c_str()); // don't check results, just try it
      BinaryWriter::writeFile(getSeqFilename(graphInfoFilename),mSeq);
      BinaryWriter::writeFile(getDictFilename(graphInfoFilename),getDict());
      BinaryWriter::writeFile(getEdgesFilename(graphInfoFilename),mEdges);
      GraphInfo gi;
      gi.K = K; gi.kSeqSize = mSeq.size();
      gi.nComponents = 0;
      if ( mEdges.size() )
        gi.nComponents = mEdges.back().getComponentID().val() + 1;
      gi.nKmers = getDict().size(); gi.nEdges = mEdges.size();
      BinaryWriter::writeFile(graphInfoFilename,gi); }

    friend bool operator==( UnipathGraph const& ug1, UnipathGraph const& ug2 )
    { return ug1.mSeq == ug2.mSeq &&
             ug1.mEdges == ug2.mEdges; }

    friend bool operator!=( UnipathGraph ug1, UnipathGraph ug2 )
    { return !(ug1 == ug2); }

    static String getInfoFilename( String const& fastb )
    { return fastb.ReplaceExtension(".fastb",".k"+ToString(K)+".ug.info"); }

    static String getSeqFilename( String const& graphInfoFilename )
    { return graphInfoFilename.ReplaceExtension(".info",".seq"); }
    static String getEdgesFilename( String const& graphInfoFilename )
    { return graphInfoFilename.ReplaceExtension(".info",".edges"); }
    static String getDictFilename( String const& graphInfoFilename )
    { return graphInfoFilename.ReplaceExtension(".info",".dict"); }
    static void removeFiles( String const& graphInfoFilename )
    { unlink(graphInfoFilename.c_str());
      unlink(getSeqFilename(graphInfoFilename).c_str());
      unlink(getEdgesFilename(graphInfoFilename).c_str());
      unlink(getDictFilename(graphInfoFilename).c_str()); }

    static void validateUnipaths( HugeBVec const& seq,
                                  UnipathEdgeVec const& edges,
                                  KDict const& dict );

private:
    UnipathGraph( UnipathGraph const& ); // unimplemented -- no copying
    UnipathGraph& operator=( UnipathGraph const& ); // unimplemented -- no copying

    void loadDict() const;

    HugeBVec mSeq;
    UnipathEdgeVec mEdges;
    String mDictFilename;
    mutable KDict* mpDict;
    mutable bool mMyDict;
};

template <unsigned K>
UnipathGraph<K>::UnipathGraph( String const& graphInfoFile )
: mSeq(getSeqFilename(graphInfoFile).c_str()),
  mDictFilename(getDictFilename(graphInfoFile)), mpDict(0), mMyDict(false)
{
    GraphInfo gi;
    BinaryReader::readFile(graphInfoFile.c_str(),&gi);

    if ( gi.K != K )
        FatalErr("K has changed.");

    if ( mSeq.size() != gi.kSeqSize )
        FatalErr("Length of kmer bases sequence (" << mSeq.size() <<
                 ") doesn't jibe with ug.info's value (" << gi.kSeqSize <<
                 ").");

    String edgesFile = getEdgesFilename(graphInfoFile);
    BinaryReader::readFile(edgesFile.c_str(),&mEdges);
    if ( mEdges.size() != gi.nEdges )
        FatalErr("Number of graph edges (" << mEdges.size() <<
                 ") doesn't jibe with ug.info's value (" << gi.nEdges << ").");
}

/// Figures out successors or predecessors to a given kmer.
template <unsigned K>
class KmerStepper
{
public:
    typedef typename KmerDict<K>::Entry DictEntry;

    KmerStepper( KmerDict<K> const& dict )
    : mDict(dict) {}

    // compiler-supplied copy ctor and dtor are OK

    KMer<K>& getSteppedKmer() { return mKMer; }

    unsigned getSuccessors( KMer<K> const& kmer,
                                    DictEntry const** pEntries )
    { unsigned result = 0;
      KMerContext context = getContext(kmer);
      mKMer = kmer;
      mKMer.toSuccessor(0);
      for ( unsigned succCode = 0; succCode < 4u; ++succCode )
      { DictEntry const* pEntry = 0;
        if ( context.isSuccessor(succCode) )
        { mKMer.setBack(succCode);
          pEntry = mDict.findEntry(mKMer); ForceAssert(pEntry);
          result += 1; }
        *pEntries++ = pEntry; }
      AssertEq(result,context.getSuccessorCount());
      return result; }

    unsigned getPredecessors( KMer<K> const& kmer,
                                        DictEntry const** pEntries )
    { unsigned result = 0;
      KMerContext context = getContext(kmer);
      mKMer = kmer;
      mKMer.toPredecessor(0);
      for ( unsigned predCode = 0; predCode < 4u; ++predCode )
      { DictEntry const* pEntry = 0;
        if ( context.isPredecessor(predCode) )
        { mKMer.setFront(predCode);
          pEntry = mDict.findEntry(mKMer); ForceAssert(pEntry);
          result += 1; }
        *pEntries++ = pEntry; }
      AssertEq(result,context.getPredecessorCount());
      return result; }

private:
    KMerContext getContext( KMer<K> const& kmer )
    { KMerContext result;
      DictEntry const* pEntry;
      if ( kmer.getCanonicalForm() == CanonicalForm::REV )
      { pEntry = mDict.findEntryCanonical(KMer<K>(kmer).rc());
        ForceAssert(pEntry);
        result = pEntry->getKDef().getContext().rc(); }
      else
      { pEntry = mDict.findEntryCanonical(kmer);
        ForceAssert(pEntry);
        result = pEntry->getKDef().getContext(); }
      return result; }

    KmerDict<K> const& mDict;
    KMer<K> mKMer;
};

template <unsigned K>
void UnipathGraph<K>::validateUnipaths( HugeBVec const& bv,
                                     UnipathEdgeVec const& edges,
                                     KmerDict<K> const& dict )
{
    std::cout << Date() << " Validating unipaths." << std::endl;
    size_t nErrors = 0;
    typedef UnipathEdgeVec::const_iterator Itr;
    typename KmerDict<K>::Entry const* entries[4];
    KmerStepper<K> stepper(dict);
    KMer<K> kmer;
    using std::equal;
    for ( Itr itr(edges.begin()), end(edges.end()); itr != end; ++itr )
    {
        UnipathEdge const& edge = *itr;
        kmer.assign(bv.begin(edge.getFinalKmerID().val()));
        stepper.getSuccessors(kmer,entries);
        for ( unsigned succCode = 0; succCode < 4u; ++succCode )
        {
            typename KmerDict<K>::Entry const* pEntry = entries[succCode];
            EdgeID edgeID = edge.getSuccessor(succCode);
            if ( !pEntry )
            {
                if ( !edgeID.isNull() )
                {
                    std::cout << "Edge " << (itr-edges.begin()) << " has a " <<
                        Base::val2Char(succCode) << " successor, edge " <<
                        edgeID << ", that wasn't found by getSuccessors." <<
                        std::endl;
                    nErrors += 1;
                }
            }
            else
            {
                KDef const& def = pEntry->getKDef();
                if ( edgeID.isNull() )
                {
                    std::cout << "Edge " << (itr-edges.begin()) <<
                            " should have a " << Base::val2Char(succCode) <<
                            " successor, edge " << def.getEdgeID() <<
                            ", but it's not marked." << std::endl;
                    nErrors += 1;
                }
                else
                {
                    if ( edgeID != def.getEdgeID() )
                    {
                        std::cout << "Edge " << (itr-edges.begin()) <<
                                " should have a " << Base::val2Char(succCode) <<
                                " successor, edge " << def.getEdgeID() <<
                                ", but it's marked as edge " << edgeID << '.'
                                << std::endl;
                        nErrors += 1;
                        edgeID = def.getEdgeID();
                    }
                    UnipathEdge const& succ = edges[edgeID.val()];
                    if ( edge.getComponentID() != succ.getComponentID() )
                    {
                        std::cout << "Edge " << (itr-edges.begin()) << " has a "
                                << Base::val2Char(succCode) << " successor, " <<
                                edgeID << ", in a different component: " <<
                                edge.getComponentID() << " vs. " <<
                                succ.getComponentID() << std::endl;
                        nErrors += 1;
                    }
                    bool isRC = edge.isSuccessorRC(succCode);
                    KmerID kid = isRC ? succ.getFinalKmerID() :
                                        succ.getInitialKmerID();
                    KmerID kid2 = succ.getKmerID(def.getEdgeOffset());
                    if ( kid != kid2 )
                    {
                        std::cout << "Edge " << (itr-edges.begin()) << " has a "
                                << Base::val2Char(succCode) << " successor, " <<
                                edgeID << ", but its adjacent kmerID " << kid <<
                                " looks up as " << kid2 << '.' << std::endl;
                        nErrors += 1;
                    }
                    KMer<K>& kmerSucc = stepper.getSteppedKmer();
                    kmerSucc.setBack(succCode);
                    if ( isRC )
                    {
                        if ( !equal(kmerSucc.rcbegin(),kmerSucc.rcend(),
                                    bv.begin(kid.val())) )
                        {
                            std::cout << "Edge " << (itr-edges.begin()) <<
                                    " has a " << Base::val2Char(succCode) <<
                                    " successor, " << edgeID <<
                           ", but the adjacent kmer has an incorrect sequence."
                                    << std::endl;
                            nErrors += 1;
                        }
                    }
                    else
                    {
                        if ( !equal(kmerSucc.begin(),kmerSucc.end(),
                                    bv.begin(kid.val())) )
                        {
                            std::cout << "Edge " << (itr-edges.begin()) <<
                                    " has a " << Base::val2Char(succCode) <<
                                    " successor, " << edgeID <<
                           ", but the adjacent kmer has an incorrect sequence."
                                    << std::endl;
                            nErrors += 1;
                        }
                    }
                }
            }
        }

        kmer.assign(bv.begin(edge.getInitialKmerID().val()));
        stepper.getPredecessors(kmer,entries);
        for ( unsigned predCode = 0; predCode < 4u; ++predCode )
        {
            typename KmerDict<K>::Entry const* pEntry = entries[predCode];
            EdgeID edgeID = edge.getPredecessor(predCode);
            if ( !pEntry )
            {
                if ( !edgeID.isNull() )
                {
                    std::cout << "Edge " << (itr-edges.begin()) << " has a " <<
                        Base::val2Char(predCode) << " predecessor, edge " <<
                        edgeID << ", that wasn't found by getPredecessors." <<
                        std::endl;
                    nErrors += 1;
                }
            }
            else
            {
                KDef const& def = pEntry->getKDef();
                if ( edgeID.isNull() )
                {
                    std::cout << "Edge " << (itr-edges.begin()) <<
                            " should have a " << Base::val2Char(predCode) <<
                            " predecessor, edge " << def.getEdgeID() <<
                            ", but it's not marked." << std::endl;
                    nErrors += 1;
                }
                else
                {
                    if ( edgeID != def.getEdgeID() )
                    {
                        std::cout << "Edge " << (itr-edges.begin()) <<
                                " should have a " << Base::val2Char(predCode) <<
                                " predecessor, edge " << def.getEdgeID() <<
                                ", but it's marked as edge " << edgeID << '.'
                            << std::endl;
                        nErrors += 1;
                        edgeID = def.getEdgeID();
                    }
                    UnipathEdge const& pred = edges[edgeID.val()];
                    if ( edge.getComponentID() != pred.getComponentID() )
                    {
                        std::cout << "Edge " << (itr-edges.begin()) << " has a "
                                << Base::val2Char(predCode) << " predecessor, "
                                << edgeID << ", in a different component: " <<
                                edge.getComponentID() << " vs. " <<
                                pred.getComponentID() << std::endl;
                        nErrors += 1;
                    }
                    bool isRC = edge.isPredecessorRC(predCode);
                    KmerID kid = isRC ? pred.getInitialKmerID():
                                        pred.getFinalKmerID();
                    KmerID kid2 = pred.getKmerID(def.getEdgeOffset());
                    if ( kid != kid2 )
                    {
                        std::cout << "Edge " << (itr-edges.begin()) << " has a "
                                << Base::val2Char(predCode) << " predecessor, "
                                << edgeID << ", but its adjacent kmerID " <<
                                kid << " looks up as " << kid2 << '.' <<
                                std::endl;
                        nErrors += 1;
                    }
                    KMer<K>& kmerPred = stepper.getSteppedKmer();
                    kmerPred.setFront(predCode);
                    if ( isRC )
                    {
                        if ( !equal(kmerPred.rcbegin(),kmerPred.rcend(),
                                    bv.begin(kid.val())) )
                        {
                            std::cout << "Edge " << (itr-edges.begin()) <<
                                    " has a " << Base::val2Char(predCode) <<
                                    " predecessor, " << edgeID <<
                           ", but the adjacent kmer has an incorrect sequence."
                                    << std::endl;
                            nErrors += 1;
                        }
                    }
                    else
                    {
                        if ( !equal(kmerPred.begin(),kmerPred.end(),
                                    bv.begin(kid.val())) )
                        {
                            std::cout << "Edge " << (itr-edges.begin()) <<
                                    " has a " << Base::val2Char(predCode) <<
                                    " predecessor, " << edgeID <<
                           ", but the adjacent kmer has an incorrect sequence."
                                    << std::endl;
                            nErrors += 1;
                        }
                    }
                }
            }
        }
    }
    ForceAssertEq(nErrors,0ul);
}

template <unsigned K>
void UnipathGraph<K>::loadDict() const
{
    std::cout << Date() << " Loading kmer dictionary." << std::endl;
    mpDict = new KmerDict<K>(0);
    BinaryReader::readFile(mDictFilename,mpDict);
    mMyDict = true;
}

/// An ID for a Read.
class ReadID : public ID<5>
{
public:
    ReadID() {}
    explicit ReadID( size_t id ) : ID<5>(id) {}

    // compiler-supplied copying and destructor are ok
};

TRIVIALLY_SERIALIZABLE(ReadID);

// A path segment index, and a bit for RC or not
class SegID : public ID<3>
{
public:
    SegID() : ID<3>(0ul) {}
    explicit SegID( size_t idx, bool isRC ) : ID<3>(2*idx+isRC) {}

    size_t getIndex() const { return val()>>1; }
    void setIndex( size_t idx ) { setVal(2*idx+isRC()); }
    bool isRC() const { return val()&1; }

    // compiler-supplied copying and destructor are ok
};

TRIVIALLY_SERIALIZABLE(SegID);

/// Describes a read that traverses some unipath graph edge.
/// (Which edge is not specified by this class: that's an external association.)
/// Also tells you which segment of the read's pathing has the edge, and whether
/// it's reverse-complement or not.
class UnipathEvidence
{
public:
    UnipathEvidence() {}
    UnipathEvidence( ReadID const& readID, size_t segID, bool isRC )
    : mReadID(readID), mSegID(segID,isRC) {}

    // compiler-supplied copying and destructor are ok

    ReadID const& getReadID() const { return mReadID; }
    size_t getSegmentID() const { return mSegID.getIndex(); }
    bool isRC() const { return mSegID.isRC(); }

    friend int compare( UnipathEvidence const& ev1, UnipathEvidence const& ev2 )
    { int result = compare(ev1.mReadID,ev2.mReadID);
      if ( !result ) result = compare(ev1.mSegID,ev2.mSegID);
      return result; }

private:
    ReadID mReadID;
    SegID mSegID;
};

typedef SerfVec<UnipathEvidence> UnipathEvidenceVec;
typedef MasterVec<UnipathEvidenceVec> VecUnipathEvidenceVec;
extern template class SmallVec< UnipathEvidence, MempoolAllocator<UnipathEvidence> >;
extern template class OuterVec<UnipathEvidenceVec>;

TRIVIALLY_SERIALIZABLE(UnipathEvidence);

/// An index into a sequence of bases.  Each base shows which turn to take at
/// each node of a UnipathGraph traversed by some sequence of bases (usually a
/// read) placed onto the graph.  This is a very terse way of representing a
/// path.
class PathID : public ID<5>
{
public:
    PathID() {}
    explicit PathID( size_t id ) : ID<5>(id) {}

    // compiler-supplied copying and destructor are ok
};

TRIVIALLY_SERIALIZABLE(PathID);

/// Describes a traversal of a UnipathGraph.
/// Used, for example, to indicate where a read is placed on the graph.
/// This is not convenient to compute with, but is a good storage format.
class Unipathing
{
public:
    Unipathing() {}
    Unipathing( EdgeID const& edgeID, bool isRC, size_t skip )
    : mEdgeID(edgeID), mInitialSkip(skip), mNSegments(1ul,isRC)
    { ForceAssertLt(skip,1ul<<32); }

    // compiler-supplied copying and destructor are ok

    bool isNull() const { return mPathID.isNull(); }
    size_t getNSegments() const { return mNSegments.getIndex(); }
    EdgeID const& getInitialEdgeID() const { return mEdgeID; }
    bool isInitialEdgeRC() const { return mNSegments.isRC(); }
    PathID const& getPathID() const { return mPathID; }
    size_t getInitialSkip() const { return mInitialSkip; }
    size_t getFinalSkip() const { return mFinalSkip; }

    void setSegments( size_t nSegments, PathID const& pathID )
    { mNSegments.setIndex(nSegments); mPathID = pathID; }

    void setFinalSkip( size_t skip )
    { ForceAssertLt(skip,1ul<<32); mFinalSkip = skip; }

    friend std::ostream& operator<<( std::ostream& os,
                                     Unipathing const& pathing )
    { os << "Edge=" << pathing.getInitialEdgeID() <<
            " IsRC=" << pathing.isInitialEdgeRC() <<
            " IniSkip=" << pathing.getInitialSkip() <<
            " FinSkip=" << pathing.getFinalSkip() <<
            " nSegs=" << static_cast<unsigned>(pathing.getNSegments()) <<
            " PathID=" << pathing.getPathID();
      return os; }

private:
    EdgeID mEdgeID; // initial edge ID
    unsigned mInitialSkip;
    unsigned mFinalSkip;
    SegID mNSegments;
    PathID mPathID;
};

TRIVIALLY_SERIALIZABLE(Unipathing);
typedef std::vector<Unipathing> UnipathingVec;

class EdgeDesc
{
public:
    EdgeDesc( EdgeID const& edgeID, CanonicalForm status )
    : mEdgeID(edgeID), mStatus(status) {}

    // compiler-supplied copying and destructor are OK

    EdgeID const& getEdgeID() const { return mEdgeID; }
    CanonicalForm getStatus() const { return mStatus; }

    EdgeDesc operator~() const
    { return EdgeDesc(mEdgeID,complement(mStatus)); }

    friend bool operator<( EdgeDesc const& ed1, EdgeDesc const& ed2 )
    { return compare(ed1,ed2) < 0; }
    friend bool operator<=( EdgeDesc const& ed1, EdgeDesc const& ed2 )
    { return compare(ed1,ed2) <= 0; }
    friend bool operator>( EdgeDesc const& ed1, EdgeDesc const& ed2 )
    { return compare(ed1,ed2) > 0; }
    friend bool operator>=( EdgeDesc const& ed1, EdgeDesc const& ed2 )
    { return compare(ed1,ed2) >= 0; }
    friend bool operator==( EdgeDesc const& ed1, EdgeDesc const& ed2 )
    { return compare(ed1,ed2) == 0; }
    friend bool operator!=( EdgeDesc const& ed1, EdgeDesc const& ed2 )
    { return compare(ed1,ed2) != 0; }
    friend int compare( EdgeDesc const& ed1, EdgeDesc const& ed2 )
    { int result = compare(ed1.mEdgeID,ed2.mEdgeID);
      if ( !result ) result = compare(ed1.mStatus,ed2.mStatus);
      return result; }

    friend ostream& operator<<( ostream& os, EdgeDesc const& ed )
    { if ( ed.mStatus == CanonicalForm::REV ) os << '~';
      return os << ed.mEdgeID; }

private:
    EdgeID mEdgeID;
    CanonicalForm mStatus;
};

/// Describes a unipath graph traversal as a sequence of EdgeDescs.
/// Also has the number of kmers not covered in the initial and final segments.
class EdgeList : public std::vector<EdgeDesc>
{
    typedef std::vector<EdgeDesc> Base;
public:
    EdgeList()
    : mInitialSkip(0), mFinalSkip(0)
    {}

    EdgeList( size_t initialSkip, size_t finalSkip )
    : mInitialSkip(initialSkip), mFinalSkip(finalSkip)
    {}

    // compiler-supplied copying and destructor are OK

    // Reverse-complement in place.
    EdgeList& rc()
    { using std::swap; swap(mInitialSkip,mFinalSkip);
      iterator head = begin(); iterator tail = end();
      while (head != tail)
      { EdgeDesc tmp = ~*--tail;
        if ( head == tail ) { *head = tmp; break; }
        *tail = ~*head; *head = tmp; ++head; }
      return *this; }

    size_t getInitialSkip() const { return mInitialSkip; }
    size_t getFinalSkip() const { return mFinalSkip; }

    friend bool operator<( EdgeList const& el1, EdgeList const& el2 )
    { return compare(el1,el2) < 0; }
    friend bool operator<=( EdgeList const& el1, EdgeList const& el2 )
    { return compare(el1,el2) <= 0; }
    friend bool operator>( EdgeList const& el1, EdgeList const& el2 )
    { return compare(el1,el2) > 0; }
    friend bool operator>=( EdgeList const& el1, EdgeList const& el2 )
    { return compare(el1,el2) >= 0; }
    friend bool operator==( EdgeList const& el1, EdgeList const& el2 )
    { return compare(el1,el2) == 0; }
    friend bool operator!=( EdgeList const& el1, EdgeList const& el2 )
    { return compare(el1,el2) != 0; }
    friend int compare( EdgeList const& el1, EdgeList const& el2 )
    { int result = compare(static_cast<Base const&>(el1),
                            static_cast<Base const&>(el2));
      if ( !result ) result = compare(el1.mInitialSkip,el2.mInitialSkip);
      if ( !result ) result = compare(el1.mFinalSkip,el2.mFinalSkip);
      return result; }

    friend ostream& operator<<( ostream& os, EdgeList const& el )
    { os << '[' << el.getInitialSkip();
      char sep = ']';
      EdgeList::const_iterator end(el.end());
      for ( EdgeList::const_iterator itr(el.begin()); itr != end; ++itr )
      { os << sep << *itr; sep = ':'; }
      os << '[' << el.getFinalSkip() << ']';
      return os; }

    friend void swap( EdgeList& el1, EdgeList& el2 )
    { el1.swap(el2);
      using std::swap;
      swap(el1.mInitialSkip,el2.mInitialSkip);
      swap(el1.mFinalSkip,el2.mFinalSkip); }

private:
    size_t mInitialSkip;
    size_t mFinalSkip;
};

/// Describes how reads align to a UnipathGraph.
template <unsigned K>
class PathCollection
{
public:
    typedef UnipathGraph<K> UGraph;

    /// Create a new graph, and path onto it.
    static void create( String const& fastbFilename,
                        bool validate,
                        unsigned nThreads,
                        size_t nKmersEst = 0 );

    PathCollection( String const& pathInfoFilename,
                    String const& graphInfoFilename );

    // copying prohibited.  compiler-supplied destructor is OK

    UGraph const& getGraph() const { return mGraph; }

    size_t getNReads() const { return mPathings.size(); }

    EdgeList getEdgeList( size_t readId ) const;
    EdgeList getEdgeList( ReadID const& readID ) const
    { return getEdgeList(readID.val()); }

    // number of segments in the path
    size_t getEdgeListSize( size_t readId ) const
    { return mPathings[readId].getNSegments(); }
    size_t getEdgeListSize( ReadID const& readID ) const
    { return mPathings[readID.val()].getNSegments(); }

    // number of kmers in the path
    size_t getEdgeLen( EdgeDesc const& ed ) const
    { return getEdgeLen(ed.getEdgeID()); }
    size_t getEdgeLen( EdgeID const& edgeID ) const
    { return mGraph.getEdge(edgeID).getLength(); }
    size_t getEdgeListLen( EdgeList const& el ) const
    { size_t len = getEdgeListLen(el.begin(),el.end());
      return len - el.getInitialSkip() - el.getFinalSkip(); }
    size_t getEdgeListLen( EdgeList::const_iterator itr,
                           EdgeList::const_iterator end ) const
    { size_t len = 0;
      while ( itr != end )
      { len += getEdgeLen(itr->getEdgeID()); ++itr; }
      return len; }

    EdgeDesc getNextEdgeDesc( EdgeDesc const& ed, unsigned base ) const
    { UnipathEdge const& edge = mGraph.getEdge(ed.getEdgeID());
      EdgeID id;
      CanonicalForm status = CanonicalForm::FWD;
      if ( ed.getStatus() == CanonicalForm::REV )
      { id = edge.getPredecessor(base);
        if ( !edge.isPredecessorRC(base) ) status = CanonicalForm::REV; }
      else
      { id = edge.getSuccessor(base);
        if ( edge.isSuccessorRC(base) ) status = CanonicalForm::REV; }
      if ( !id.isNull() && mGraph.getEdge(id).isPalindrome() )
        status = CanonicalForm::PALINDROME;
      return EdgeDesc(id,status); }

    bvec getBases( EdgeList const& el ) const
    { bvec result; getBases(el,&result); return result; }

    bvec& getBases( EdgeList const& el, bvec* pBV ) const;

    void validate( String const& fastbFilename );

    static String getInfoFilename( String const& fastb )
    { return fastb.ReplaceExtension(".fastb",".k"+ToString(K)+".pc.info"); }

    static String getPathsFilename( String const& infoFile )
    { return infoFile.ReplaceExtension(".info",".paths"); }

    static String getPathseqFilename( String const& infoFile )
    { return infoFile.ReplaceExtension(".info",".seq"); }

    static String getEvidenceFilename( String const& infoFile )
    { return infoFile.ReplaceExtension(".info",".ev"); }

private:
    PathCollection( PathCollection const& ); // unimplemented -- no copying
    PathCollection& operator=( PathCollection const& ); // unimplemented -- no copying

    UGraph mGraph;
    HugeBVec mPathSeq; // indexed by PathID
    UnipathingVec mPathings; // indexed by ReadID
};

// Data about the pathing files stuffed into the .info file.
struct PathInfo
{
    PathInfo() : mVersion(CURRENT_VERSION), mNReads(0), mPathSeqSize(0) {}
    PathInfo( size_t nReads, size_t pathSeqSize )
    : mVersion(CURRENT_VERSION), mNReads(nReads), mPathSeqSize(pathSeqSize) {}

    // compiler-supplied copying and destructor are OK

    size_t mVersion;
    size_t mNReads;
    size_t mPathSeqSize;
    static size_t const CURRENT_VERSION = 1;
};

TRIVIALLY_SERIALIZABLE(PathInfo);

template <unsigned K>
PathCollection<K>::PathCollection( String const& pathInfoFilename,
                                  String const& graphInfoFilename )
: mGraph(graphInfoFilename),
  mPathSeq(getPathseqFilename(pathInfoFilename).c_str())
{
    String pathsFile = getPathsFilename(pathInfoFilename);
    BinaryReader::readFile(pathsFile.c_str(),&mPathings);
    PathInfo pi;
    BinaryReader::readFile(pathInfoFilename.c_str(),&pi);
    if ( PathInfo::CURRENT_VERSION != pi.mVersion )
        FatalErr("PathInfo was written by an earlier version of the code.");
    if ( mPathings.size() != pi.mNReads )
        FatalErr("Number of traversals (" << mPathings.size()
                 << ") doesn't jibe with pathinfo's value ("
                 << pi.mNReads << ").");
    if ( mPathSeq.size() != pi.mPathSeqSize )
        FatalErr("Length of traversal bases sequence (" <<
                 mPathSeq.size() <<
                 ") doesn't jibe with pathinfo's value (" << pi.mPathSeqSize <<
                 ").");
}

template <unsigned K>
EdgeList PathCollection<K>::getEdgeList( size_t readId ) const
{
    Unipathing const& pathing = mPathings[readId];
    EdgeList result(pathing.getInitialSkip(),pathing.getFinalSkip());
    unsigned nSegs = pathing.getNSegments();
    if ( nSegs )
    {
        result.reserve(nSegs);
        EdgeID edgeID = pathing.getInitialEdgeID();
        bool isRC = pathing.isInitialEdgeRC();
        bool isPalindrome = mGraph.getEdge(edgeID).isPalindrome();
        CanonicalForm status = isPalindrome ?
                                CanonicalForm::PALINDROME :
                                isRC ? CanonicalForm::REV : CanonicalForm::FWD;
        result.push_back(EdgeDesc(edgeID,status));
        if ( nSegs > 1 )
        {
            size_t off = pathing.getPathID().val()-1;
            HugeBVec::const_iterator itr(mPathSeq.begin(off));
            while ( --nSegs )
            {
                isRC = mGraph.nextEdge(*++itr,&edgeID,isRC);
                if ( edgeID.isNull() )
                    FatalErr("Invalid successor base in path sequence for read "
                             << readId );
                bool isPalindrome = mGraph.getEdge(edgeID).isPalindrome();
                CanonicalForm status = isPalindrome ?
                                CanonicalForm::PALINDROME :
                                isRC ? CanonicalForm::REV : CanonicalForm::FWD;
                result.push_back(EdgeDesc(edgeID,status));
            }
        }
    }
    return result;
}

template <unsigned K>
bvec& PathCollection<K>::getBases( EdgeList const& el, bvec* pBV ) const
{
    pBV->clear().reserve(getEdgeListLen(el));

    typedef HugeBVec::const_iterator Itr;
    for ( unsigned idx = 0; idx < el.size(); ++idx )
    {
        EdgeDesc const& ed = el[idx];
        UnipathEdge const& edge = mGraph.getEdge(ed.getEdgeID());
        Itr bItr = mGraph.getBases(edge.getInitialKmerID());
        Itr bEnd = mGraph.getBases(edge.getFinalKmerID()) + K;
        if ( ed.getStatus() != CanonicalForm::REV )
        {
            bItr += idx ? (K-1) : el.getInitialSkip();
            if ( idx+1 == el.size() ) bEnd -= el.getFinalSkip();
            pBV->append(bItr,bEnd);
        }
        else
        {
            bEnd -= idx ? (K-1) : el.getInitialSkip();
            if ( idx+1 == el.size() ) bItr += el.getFinalSkip();
            while ( bEnd != bItr )
                pBV->push_back(GetComplementaryBase(*--bEnd));
        }
    }
    return *pBV;
}

template <unsigned K>
void PathCollection<K>::validate( String const& fastbFilename )
{
    std::cout << Date() << " Validating pathings." << std::endl;
    VirtualMasterVec<bvec> vmv(fastbFilename.c_str());
    size_t nnn = vmv.size();
    bvec scratch;
    for ( size_t idx = 0; idx < nnn; ++idx )
    {
        EdgeList el(getEdgeList(idx));
        if ( el.size() )
            ForceAssertEq(vmv[idx],getBases(el,&scratch));
        else
            ForceAssertLt(vmv[idx].size(),K);
    }
}

template <unsigned K>
class PathsWithEvidence : public PathCollection<K>
{
public:
    PathsWithEvidence( String const& pathInfoFilename,
                       String const& graphInfoFilename );

    // copying prohibited by base class.  compiler-supplied destructor is OK

    size_t getNEdges() const { return mEvidence.size(); }

    UnipathEvidenceVec const& getEvidence( EdgeID const& edgeID ) const
    { return mEvidence[edgeID.val()]; }

    struct ReadLoc
    {
        size_t readId;
        size_t readOffset;
    };
    typedef std::vector<ReadLoc> ReadLocVec;
    ReadLocVec getKmerLocations( KMer<K> const& kmer );

private:
    VecUnipathEvidenceVec mEvidence; // indexed by EdgeID
};

template <unsigned K>
PathsWithEvidence<K>::PathsWithEvidence( String const& pathInfoFilename,
                                        String const& graphInfoFilename )
: PathCollection<K>(pathInfoFilename,graphInfoFilename)
{
    String evidenceFile =
            PathCollection<K>::getEvidenceFilename(pathInfoFilename);
    BinaryReader::readFile(evidenceFile.c_str(),&mEvidence);
    size_t nEdges = PathCollection<K>::getGraph().getNEdges();
    if ( mEvidence.size() != nEdges )
        FatalErr("Length of evidence vector (" << mEvidence.size() <<
                 ") doesn't jibe with the number of edges in the graph (" <<
                 nEdges << ").");
}

template <unsigned K>
typename PathsWithEvidence<K>::ReadLocVec
PathsWithEvidence<K>::getKmerLocations( KMer<K> const& kmer )
{
    ReadLocVec result;
    UnipathGraph<K> const& graph = this->getGraph();
    KmerDict<K> const& dict = graph.getDict();
    KDef const* pDef = dict.lookup(kmer);
    if ( pDef )
    {
        unsigned kmerOffset = pDef->getEdgeOffset();
        UnipathEvidenceVec const& evVec = getEvidence(pDef->getEdgeID());
        result.reserve(evVec.size());
        typedef UnipathEvidenceVec::const_iterator Itr;
        for ( Itr itr(evVec.begin()), end(evVec.end()); itr != end; ++itr )
        {
            EdgeList el = this->getEdgeList(itr->getReadID());
            unsigned seg = itr->getSegmentID();
            if ( !seg && kmerOffset < el.getInitialSkip() )
                return -1L;
            if ( seg == el.size()-1 && kmerOffset >= this->getEdgeLen(el.back().getEdgeID())-el.getFinalSkip() )
                return -1L;
            size_t offset = this->getEdgeLen(el.begin(),el.begin()+seg)
                                    - el.getInitialSkip() + kmerOffset;
            result.push_back(ReadLoc(itr->getReadID().val(),offset));
        }
    }
    return result;
}

#endif /* KMERS_READ_PATHER_H_ */
