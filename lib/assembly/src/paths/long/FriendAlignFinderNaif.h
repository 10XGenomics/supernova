///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS


#ifndef FRIEND_ALIGN_FINDER_NAIF_H
#define FRIEND_ALIGN_FINDER_NAIF_H

#include "ParallelVecUtilities.h"
#include "feudal/FilesOutputIterator.h"
#include "paths/long/FriendAligns.h"
#include "kmers/naif_kmer/Kmers.h"
#include "kmers/naif_kmer/NaifKmerizer.h"
#include "dna/Bases.h"
#include <unistd.h>

struct Alignment
{
    Alignment()=default;
    Alignment( uint32_t readId1, uint32_t readId2, int16_t offset,
                bool isRC )
    : mReadId1(readId1), mReadId2(readId2), mOffset(offset), mIsRC(isRC) {}

    uint32_t mReadId1;
    uint32_t mReadId2;
    int16_t mOffset;
    bool mIsRC;

    friend bool operator<( Alignment const& a1, Alignment const& a2 )
    { if ( a1.mReadId1 < a2.mReadId1 ) return true;
      if ( a1.mReadId1 > a2.mReadId1 ) return false;
      if ( a1.mReadId2 < a2.mReadId2 ) return true;
      if ( a1.mReadId2 > a2.mReadId2 ) return false;
      if ( a1.mOffset < a2.mOffset ) return true;
      if ( a1.mOffset > a2.mOffset ) return false;
      if ( a1.mIsRC < a2.mIsRC ) return true;
      return false; }

    friend bool operator==( Alignment const& a1, Alignment const& a2 )
    { return a1.mReadId1 == a2.mReadId1 && a1.mReadId2 == a2.mReadId2
            && a1.mOffset == a2.mOffset && a1.mIsRC == a2.mIsRC; }
};

// functor for comparing based on ReadId2 THEN ReadId1, unlike the default
// comparison within the Alignment class
struct R2StrictWeakOrdering {
    bool operator()( Alignment const& a1, Alignment const& a2 )
    { if ( a1.mReadId2 < a2.mReadId2 ) return true;
      if ( a1.mReadId2 > a2.mReadId2 ) return false;
      if ( a1.mReadId1 < a2.mReadId1 ) return true;
      if ( a1.mReadId1 > a2.mReadId1 ) return false;
      if ( a1.mOffset < a2.mOffset ) return true;
      if ( a1.mOffset > a2.mOffset ) return false;
      if ( a1.mIsRC < a2.mIsRC ) return true;
      return false; }
};

TRIVIALLY_SERIALIZABLE(Alignment);
TRIVIALLY_SERIALIZABLE(Friend);


// ---- KmerLCLoc -- a kmer, its location, and the character to the left (LC)
//
template <class KMER_t>
class KmerLcLoc : public KMER_t, public BVLoc
{
private:
    uint8_t _lc_start : 2;
    uint8_t _lc_stop : 2;
    uint8_t _lc_null : 1;

public:
    KmerLcLoc(const KMER_t & kmer) : KMER_t(kmer), BVLoc(), _lc_start(0),_lc_stop(0),_lc_null(0) {}
    explicit KmerLcLoc(const unsigned K = 0) : KMER_t(K), BVLoc(), _lc_start(0),_lc_stop(0),_lc_null(0)  {}

    friend
    bool operator < (const KmerLcLoc & a, const KmerLcLoc & b)
    {
	// we want kmers group togther, then sorted by LC, then position
	if (static_cast<const KMER_t &>(a) < static_cast<const KMER_t &>(b)) return true;
	if (static_cast<const KMER_t &>(b) < static_cast<const KMER_t &>(a)) return false;
	if ( !a.is_lc_null() &&  b.is_lc_null() ) return true;			// null sorts toward the end
	if (  a.is_lc_null() && !b.is_lc_null() ) return false;
	if ( !a.is_lc_null() && !b.is_lc_null() ) {
	    if ( a.get_lc_start() < b.get_lc_start() ) return true;
	    if ( a.get_lc_start() > b.get_lc_start() ) return false;
	    if ( a.get_lc_stop() < b.get_lc_stop() ) return true;
	    if ( a.get_lc_stop() > b.get_lc_stop() ) return false;
	}
	return (static_cast<const BVLoc &>(a) < static_cast<const BVLoc &>(b));
    }

    void set_lc_null() { _lc_null = true; }
    bool is_lc_null() const { return _lc_null; }

    void set_lc_start( unsigned char start ) { _lc_start = start; }
    void set_lc_stop( unsigned char stop ) { _lc_stop = stop; }

    unsigned char get_lc_start() const { return _lc_start; }
    unsigned char get_lc_stop() const { return _lc_stop; }
};


// ---- KernelFriendFinder -- kmerizes reads into KmerLcLocs
//
template <class ELEM_t>
class KernelFriendFinder
{
private:
    const BaseVecVec   & _bvv;
    const size_t         _K;

    vec<Alignment>	 _aligns;
    LockedData           _lock;     // lock for merge()

    FilesOutput&	 _outfiles;
    size_t               _max_aligns;

    typedef typename ELEM_t::kmer_type  kmer_type;

public:
    typedef KmerLcLoc<kmer_type> rec_type;		// force KmerLcLoc records

    KernelFriendFinder(const BaseVecVec & bvv, 
                   const size_t       K,
                   FilesOutput&		fo,
                   int max_aligns ) :
	_bvv(bvv), 
	_K(K), 
	_aligns(),
	_lock(),
        _outfiles(fo),
        _max_aligns(max_aligns)
	{ ForceAssertGt(max_aligns,0); }

    // copy constructor for temporary kernels
    explicit KernelFriendFinder(const KernelFriendFinder<ELEM_t> & that) :
	_bvv(that._bvv),
	_K(that._K),
	_aligns(),
	_lock(),
	_outfiles(that._outfiles),
        _max_aligns(that._max_aligns)
	{}

    // interface function needed by naif_kmerize()
    size_t K() const { return _K; }

    // interface function needed by naif_kmerize()
    const BaseVecVec & bases() const { return _bvv; }


    // interface function needed by naif_kmerize()
    void parse_base_vec(ParcelBuffer<rec_type> * p_buf, 
		      const size_t ibv)
    {
	ParcelBuffer<rec_type> & parcel_buf = *p_buf;
	SubKmers<BaseVec, kmer_type> kmer_cur(_K, _bvv[ibv]);
	while (kmer_cur.not_done()) {
	    rec_type kmer = kmer_cur.canonical();

	    if (parcel_buf.in_one_parcel(kmer)) {

#if 0
		// DEBUG
		cout << "---------------------------" << endl;
		cout << "KMER=" << kmer.to_string() << endl;
		if ( kmer_cur.is_canonical_fw() ) {
		    size_t start = kmer_cur.index_start();
		    size_t stop = kmer_cur.index_stop();
		    cout << "FORW=";
		    for ( size_t i = start; i < stop; ++i )
			cout <<  Base::val2Char(_bvv[ibv][i]);
		} else {
		    size_t start = kmer_cur.index_start();
		    size_t stop = kmer_cur.index_stop();
		    cout << "BACK=";
		    for ( size_t i = stop-1; i >= start; --i )
			cout <<  Base::val2Char((3^_bvv[ibv][i]));
		}
		cout << endl;
		// END OF DEBUG
#endif

		size_t ib1, ib2;

		// set ib to one base before the fw k-mer (if fw is
		// canonical) or one base before the rc k-mer.
		size_t blen = _bvv[ibv].size();

		if ( kmer_cur.is_canonical_fw() ) {
		    ib1 = kmer_cur.index_start();
		    ib2 = kmer_cur.index_stop();
		    if ( ib1 == 0 || ib2 >= blen ) kmer.set_lc_null();
		    else {
			kmer.set_lc_start( _bvv[ibv][ib1-1] );
			kmer.set_lc_stop(  _bvv[ibv][ib2] );
		    }
		    kmer.set_fw(true);
		} else {
		    ib1 = kmer_cur.index_start_reverse();
		    ib2 = kmer_cur.index_stop_reverse();

		    if ( ib1 == 0 || ib2 >= blen ) kmer.set_lc_null();			// its fine that these limits are for the RC, because we're checking both sides
		    else {
			kmer.set_lc_start( 3^(_bvv[ibv][kmer_cur.index_stop()]) );	// stop is the base just after on fw, or just before on rc, as in this case
			kmer.set_lc_stop(  3^(_bvv[ibv][kmer_cur.index_start()-1]) );
		    }

		    kmer.set_rc(true);
		}

		kmer.set_ibv(ibv);
		kmer.set_ib(ib1);
		parcel_buf.add(kmer);
	    }
	    kmer_cur.next();          
	}
    }

    // validate_align -- ensure that there's at least a K-base match between the sequences using the given alignment
    // adapted from FriendAlignFinder3::ValidateAlign
    bool validate_align( const Alignment& a ) const {
        typedef bvec::iterator Itr;
        bvec read1 = _bvv[a.mReadId1];
        bvec read2 = _bvv[a.mReadId2];
        if ( a.mIsRC ) read2.ReverseComplement();
        Itr it1 = read1.begin();
        Itr it2;
        it2 = read2.begin();
        if ( a.mOffset > 0 )
            it1 += a.mOffset;
        else
            it2 -= a.mOffset;

        Itr end = it1+std::min(read1.end()-it1,read2.end()-it2);
        bool find_match = false;
        while ( it1 + _K <= end ) {
            pair<Itr,Itr> mis_locs = mismatch( it1, end, it2 );
            if ( static_cast<size_t>(mis_locs.first - it1) >= _K ) {	// cast is to quell warnings, and it can't be negative by design
                find_match = true;
                break;
            } else {
                it1 = mis_locs.first + 1;
                it2 = mis_locs.second + 1;
            }
        }
        return find_match;
    }



    // interface function needed by naif_kmerize()
    void summarize(const vec<rec_type> & krecs,  
		   const size_t i0k,
		   const size_t i1k)
    {
//	vec<Alignment> aligns;

#if 0
	Locker lock2( _outfiles.getLockedData() );
	cout << "KMER=" << krecs[i0k].to_string() << endl;
	for (size_t i = i0k; i < i1k; ++i ) {
	    cout << "XXXX=" << krecs[i].to_string() << "|" << krecs[i].get_lc_char() << "|" << krecs[i].is_rc() << endl;
	    size_t ibv = krecs[i].ibv();
	    size_t ib  = krecs[i].ib();
	    BaseVec bv = _bvv[ibv];

	    if ( krecs[i].is_rc() ) bv.ReverseComplement();

	    cout << "KMER ";
	    for ( size_t i = ib; i < ib+_K; ++i )
		cout << Base::val2Char(bv[i]);
	    cout << endl;
	}
	return;
#endif

	if ( i1k-i0k > _max_aligns ) return;
/*
        // dumb way to count how much output we might generate
        // for now, silently leaves if we exceed the max.  The thought is:
        // if there are too many alignments, then this is a pathological k-mer and we can
        // pick our friends more wisely.
	int count = 0;
	for (size_t i = i0k; i != i1k; ++i ) {		// process a single kmer
	    for ( size_t j = i0k; j != i1k; ++j ) {
                if ( i == j ) continue;

	        const rec_type& r1 = krecs[i];
		const rec_type& r2 = krecs[j];

		if ( !r1.is_lc_null() && !r2.is_lc_null()
			&& r1.get_lc_start() == r2.get_lc_start()
			&& r1.get_lc_stop() == r2.get_lc_stop() ) continue;
		if ( ++count > _max_aligns ) return;
	    }
	}
*/
	for (size_t i = i0k; i != i1k; ++i ) {		// process a single kmer
	    for ( size_t j = i+1; j != i1k; ++j ) {

	        const rec_type& r1 = krecs[i];
		const rec_type& r2 = krecs[j];
                if ( r1.ibv() == r2.ibv() ) continue;

		if ( !r1.is_lc_null() && !r2.is_lc_null()
			&& r1.get_lc_start() == r2.get_lc_start()
			&& r1.get_lc_stop() == r2.get_lc_stop() ) continue;

//		cout << "i=" << i << ", j=" << j << ", r1_lc=" << r1.get_lc_char() << ", r2_lc=" << r2.get_lc_char() << endl;

		// okay, we have a valid pair -- emit an alignment for them
		bool rc = r1.is_rc() ^ r2.is_rc();
		int offset;

		if ( !r1.is_rc() )			// the asymmetry is this - only R2 ever gets flipped by the rc flag
		    offset = r1.ib() - r2.ib();
		else
		    offset = r2.ib() - r1.ib();

		Alignment align( r1.ibv(), r2.ibv(), offset, rc );

//		cout << "r1.is_rc=" << r1.is_rc() << ", r2.is_rc=" << r2.is_rc() << ", rc=" << rc << ", offset=" << offset << ", VALID=" << validate_align( align ) << endl;

		_aligns.push_back( align );
	    }
	}

#if 0			// WRITE OUTPUT HERE
	Locker lock( _outfiles.getLockedData() );
#if 0
	size_t count = 0;
	for ( auto itr = aligns.begin(); itr != aligns.end(); ++itr ) {
	    cout << "count=" << ++count << endl;
	    ForceAssert( validate_align( *itr ) );
	}
#endif

	// need to write out aligns here
	auto iter = _outfiles.getIterator<Alignment>(&aligns[0]);
	std::copy( aligns.begin(), aligns.end(), iter );
#endif

    }


    // interface function needed by naif_kmerize()
    void merge(const KernelFriendFinder<ELEM_t> & kernel_tmp,
	       const size_t i_parcel)
    {
	Locker lock(_lock);
	Locker lock2( _outfiles.getLockedData() );
//	cout << endl << "MERGE i_parcel=" << i_parcel << ", numaligns=" << kernel_tmp._aligns.size();
	auto iter = _outfiles.getIterator((Alignment*)nullptr);
	std::copy( kernel_tmp._aligns.begin(), kernel_tmp._aligns.end(), iter );
//	for ( auto inItr = kernel_tmp._aligns.begin(); inItr != kernel_tmp._aligns.end(); ++inItr )
//	    *iter = *inItr;
    }

};


template <class KMER_T>
class FriendAlignFinderNaif : public FriendAlignerImpl {

public:
    FriendAlignFinderNaif( String const& friendsCache, int K, const vecbvec& reads,
                                const int max_freq, size_t NUM_THREADS = 0 )
    : mFriendsFile(goNuts(friendsCache,K,reads,max_freq, NUM_THREADS)), mFriends(mFriendsFile)
    {  unlink(mFriendsFile.c_str());  }

    FriendAlignFinderNaif( const FriendAlignFinderNaif& ) = delete;
    FriendAlignFinderNaif& operator=( const FriendAlignFinderNaif& ) = delete;

    // Find all alignments of one read
    virtual void getAligns( size_t readId, Friends* pFriends )
    { *pFriends = mFriends[readId]; }

private:
    static String goNuts( String const& friendsCache, int K, vecbvec const& reads,
                            int max_freq, size_t NUM_THREADS = 0 )
    {
        String tmpDir(friendsCache.SafeBeforeLast("/"));
        // note about memory usage: we'll double the Alignments in memory and double again to ParallelSort
        FilesOutput fo(tmpDir+"/tmp",String("aligns"),0.2*MemAvailable()/sizeof(Alignment));
        KernelFriendFinder<KmerLcLoc<KMER_T>> finder(reads, K, fo, max_freq);
        naif_kmerize(&finder,NUM_THREADS,1);
        fo.close();
        std::vector<String> files = fo.getFiles();
        size_t nReads = reads.size();
        size_t nFriends = 0U;
        cout << Date() << ": about to convert aligns to friends" << endl;
        if ( files.size() == 1 )
            alignsToFriends(files[0],friendsCache,nReads,nFriends);
        else
        {
            for ( auto itr=files.begin(),end=files.end(); itr != end; ++itr )
            {
                String friendsFile  = itr->ReplaceExtension(".aligns",".friends");
                alignsToFriends(*itr,friendsFile,nReads,nFriends);
                *itr = friendsFile;
            }
            std::vector<VirtualMasterVec<SerfVec<Friend>>*> inputs;
            inputs.reserve(files.size());
            for ( auto itr=files.begin(),end=files.end(); itr != end; ++itr )
                inputs.push_back(new VirtualMasterVec<SerfVec<Friend>>(*itr));
            IncrementalWriter<SerfVec<Friend>> writer(friendsCache,nReads);
            SerfVec<Friend> readFriends;
            for ( size_t idx = 0; idx != nReads; ++idx )
            {
                readFriends.clear();
                for ( auto itr=inputs.begin(),end=inputs.end();itr!=end;++itr )
                {
                    SerfVec<Friend> someFriends( (**itr)[idx] );
                    readFriends.append(someFriends.begin(),someFriends.end());
                }
                writer.add(readFriends);
            }
            writer.close();
            for ( auto itr=inputs.begin(),end=inputs.end(); itr != end; ++itr )
                delete *itr;
            for ( auto itr=files.begin(),end=files.end(); itr != end; ++itr )
                unlink(itr->c_str()) ;
        }
	cout << Date() << ": mean number of friends per read is " << nFriends / nReads << endl;
        cout << Date() << ": done converting aligns to friends" << endl;
        return friendsCache;
    }

    static void alignsToFriends( String const& input, String const& output,
                                    size_t const nReads, size_t& nFriends )
    {
	cout << Date() << ": alignsToFriends called " << input << " -> " << output
		<< " for " << nReads << " reads." << endl;
        vec<Alignment> aligns;
        FileReader fr(input);
        size_t sz = fr.getSize();
        ForceAssertEq(sz%sizeof(Alignment),0ul);
	size_t alloc_size = sz/sizeof(Alignment);
        cout << Date() << ": about to resize for " << alloc_size << endl;
        aligns.resize( 2 * alloc_size );
        cout << Date() << ": reading in " << sz << " bytes" << endl;
        fr.read( &aligns[0], sz );
        fr.close();
        cout << Date() << ": done reading." << endl;
        unlink(input.c_str());
        cout << Date() << ": populating reverse Alignments" << endl;

        for ( size_t i = 0; i < alloc_size; ++i ) {
            const Alignment& in = aligns[i];
            Alignment& out = aligns[i+alloc_size];
            out.mReadId1 = in.mReadId2;
            out.mReadId2 = in.mReadId1;
            out.mIsRC = in.mIsRC;
            // there are multiple reasons to negate the offset which add up to an even number of
            // negations when the Alignment is RC.
            out.mOffset = in.mIsRC ? in.mOffset : -in.mOffset;
        }



        cout << Date() << ": about to parallel sort " << aligns.size() << " aligns" << endl;
        ParallelSort(aligns);
        cout << Date() << ": before unique-erase " << aligns.size() << " aligns" << endl;
        aligns.erase(std::unique(aligns.begin(),aligns.end()),aligns.end());
        cout << Date() << ": after unique-erase " << aligns.size() << " aligns" << endl;

        IncrementalWriter<SerfVec<Friend>> wtr(output,nReads);
        SerfVec<Friend> friends;
        auto itr = aligns.begin();
        auto end = aligns.end();
        for ( size_t idx = 0; idx != nReads; ++idx )
        {
            while ( itr != end )
            {
                if ( itr->mReadId1 != idx )
                    break;
                friends.push_back(Friend(itr->mReadId2,itr->mOffset,itr->mIsRC));
                ++itr;
            }
            wtr.add(friends);
            nFriends += friends.size();
            friends.clear();
        }
        ForceAssert(itr==end);


//        ParallelSort(aligns, R2Comp());

    }

    String mFriendsFile;
    MasterVec<SerfVec<Friend>> mFriends;
};

#endif //FRIEND_ALIGN_FINDER_NAIF_H
