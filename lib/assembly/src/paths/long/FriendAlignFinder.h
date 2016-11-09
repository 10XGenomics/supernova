///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
#ifndef FRIEND_ALIGNS_FINDER3_H
#define FRIEND_ALIGNS_FINDER3_H

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "paths/long/FriendAligns.h"
#include "paths/long/MakeAlignments.h"
#include "kmers/ReadPather.h"
#include "kmers/ReadPatherDefs.h"
#include "Vec.h"
#include "ParallelVecUtilities.h"

// ====================== FriendAlignFinder implementations =============================
template <int K>
class FriendAlignFinder : public FriendAlignerImpl {
public:

    // Explicitly if a read align is actually valid. Used to removed some
    // false positives.
    bool ValidateAlign( const simple_align_data& a ) const {
        typedef bvec::const_iterator Itr;
        bvec const& read1 = mReads[a.id1];
        bvec read2RC;
        bvec const& read2 = a.rc2 ?
                read2RC.ReverseComplement(mReads[a.id2]) :
                mReads[a.id2];
        Itr it1 = read1.begin();
        Itr it2 = read2.begin();
        if ( a.offset > 0 )
        {
            ForceAssertLt(static_cast<unsigned>(a.offset),read1.size());
            it1 += a.offset;
        }
        else
        {
            ForceAssertLt(static_cast<unsigned>(-a.offset),read2.size());
            it2 -= a.offset;
        }
        Itr end = it1 + std::min(std::distance(it1,read1.end()),
                                    std::distance(it2,read2.end()));
        bool find_match = false;
        while ( !find_match ) {
            pair<Itr,Itr> mis_locs = mismatch( it1, end, it2 );
            if ( mis_locs.first - it1 >= K )
                find_match = true;
            if ( mis_locs.first == end )
                break;
            it1 = mis_locs.first + 1;
            it2 = mis_locs.second + 1;
        }
        return find_match;
    }

    // Data structure to store read location on the unipath
    struct ReadLocOnUnipath {
        EdgeID uid;
        short int start;
        unsigned int rid;
        bool rc;

        friend ostream& operator<<( ostream& out, const ReadLocOnUnipath& a ) {
            char dir = ( a.rc ? '-' : '+' );
            out << "loc " << a.rid<< "(" << a.start << "," << ")" << dir << "@" << a.uid.val() ;
            return out;
        }
        friend bool operator<( const ReadLocOnUnipath& l, const ReadLocOnUnipath& r)  {
            if ( l.uid != r.uid ) return l.uid < r.uid;
            return l.start < r.start;
        }
    };

    // We could use more efficient containers. But just use vec of vec for this moment.
    typedef vec<ReadLocOnUnipath> PathSegVec;
    typedef vec<ReadLocOnUnipath> ReadULocVec;

    FriendAlignFinder (const vecbvec& reads, const int max_freq = 1000, Bool use_down_sampling = False, int verb = 1 )
        : mReads(reads), mpDict(NULL), mpGraph(NULL), 
          mCopyNumberMax( max_freq ), mUseDownSampling(use_down_sampling), mVerbose(verb) 
    { Init(); }
    FriendAlignFinder( const FriendAlignFinder& )=delete;
    FriendAlignFinder& operator= ( const FriendAlignFinder& )=delete;

    virtual ~FriendAlignFinder () { delete mpGraph; delete mpDict; }

    // Find all alignments of one read
    virtual void getAligns( size_t readId, Friends* pFriends )
    { vec<ReadLocOnUnipath> locvec;
      PathOneRead( readId, &locvec );
      // Check for bad alignments
      return GetAlignsOneReadUnsorted( readId, locvec, pFriends ); }

private:

    void Init( unsigned int coverage = 5, unsigned int nThreads = 0) {
        // ========= build the kmer dictionary  =========
        if ( mVerbose >= 1 )
             cout << Date() << ": creating dictionary" << endl;
        size_t dictSize = mReads.SizeSum() / coverage;
        mpDict = new KmerDict<K> ( 5*dictSize/4 );
        mpDict->process(mReads,mVerbose,false,nThreads,100);
        if ( mVerbose >= 1 ) {
            cout << Date( ) << ": there are " << mpDict->size()
                    << " kmers (expected ~" << dictSize << ")" << endl;
            ReportMemUsage();
        }
        size_t old_dict_size = mpDict->size();
        //mpDict->clean( typename KmerDict<K>::BadKmerCountFunctor(2, mCopyNumberMax)); 
        // the result seems  better without using mCopyNumberMax at this stage
        mpDict->clean( typename KmerDict<K>::BadKmerCountFunctor(2));
        if ( mVerbose >= 1 ) 
            cout << Date() << ": Cleaning bad kmers, keeping " << ToStringAddCommas(mpDict->size())
                           << "(" << ( mpDict->size() * 100 / old_dict_size ) << "%)" << endl;

        // ========= build the unipath graph    =========
        mpGraph = new UnipathGraph<K>(*mpDict, mVerbose);
        if ( mVerbose >= 1 ) ReportMemUsage();

        // ======== Index the read locs on unipaths
        GenerateReadLocs();
        if ( mVerbose >= 1 ) ReportMemUsage();
    }

    void GenerateReadLocs( ) {
        int64_t total_locs_deleted = 0;
        mULocs.clear_and_resize( mpGraph->getNEdges() );
        #pragma omp parallel for schedule(dynamic, 100)
        for( size_t iread = 0; iread < mReads.size(); iread++ ) {
            vec<ReadLocOnUnipath> locvec;
            PathOneRead( iread, &locvec );
            int num_locs_deleted = 0;
            if ( mUseDownSampling )
                num_locs_deleted = DownSampleLocsOfOneRead( &locvec );
            #pragma omp critical
            {
                total_locs_deleted += num_locs_deleted;
                for ( size_t j = 0; j < locvec.size(); ++j ) 
                    mULocs[ locvec[j].uid.val() ].push_back( locvec[j] );
            }
        }
        #pragma omp parallel for schedule(dynamic, 100)
        for ( size_t i = 0; i < mULocs.size(); ++i ) {
            int ulen = mpGraph->getEdge( EdgeID(i) ).getLength();
            if ( ulen < 5 && mULocs[i].isize() > mCopyNumberMax )
                mULocs[i].clear();
            else
                Sort(mULocs[i]);
        }
        uint64_t total = SizeSum( mULocs );
        if ( mVerbose >= 1 ) 
            cout << Date() << ": Found " << ToStringAddCommas( total ) << " locs" 
                 << " after deleting " << ToStringAddCommas(total_locs_deleted) << endl;
    }

    void GetAlignsOneReadUnsorted( size_t read_id,
                                   const vec<ReadLocOnUnipath>& locvec,
                                   Friends *pFriends ) const {
        std::set<simple_align_data> uniq_aligns;
        for( size_t i = 0; i < locvec.size(); ++i ) {
            int nkmer1 = mReads[read_id].size() - K + 1;
            const ReadLocOnUnipath& loc1 = locvec[i];
            int stop1 = loc1.start + nkmer1;
            const ReadULocVec& ulocvec = mULocs[ loc1.uid.val() ];
            bool isPalindrome = mpGraph->getEdge(loc1.uid).isPalindrome();
            for ( size_t x2 = 0; x2 < ulocvec.size(); ++x2 ) {
                const ReadLocOnUnipath& loc2 = ulocvec[x2];
                int nkmer2 = mReads[loc2.rid].size() - K + 1;
                int stop2 = loc2.start + nkmer2;
                if ( loc2.rid == loc1.rid ) continue;
                if ( stop2 <= loc1.start ) continue;
                if ( loc2.start >= stop1 ) continue;
                {  // for all cases
                    Bool rc = loc2.rc ^ loc1.rc;
                    int offset2 = ( loc1.rc ? stop1 - stop2 : loc2.start - loc1.start );
                    simple_align_data a(loc1.rid, loc2.rid, offset2, rc);
                    uniq_aligns.insert( a );
                    //if ( ! ValidateAlign(a) ) {
                    //    #pragma omp critical 
                    //    {
                    //        cout << "Could not validate alignment ";
                    //        cout << "read1 on " << loc1.start << "," << stop1
                    //             << " read2 on " << loc2.start << "," << stop2 << endl;
                    //        int ulen = mpGraph->getEdge( loc1.uid ).getLength();
                    //        cout << "ulen= " << ulen << endl;
                    //    }
                    //}
                }
                // Special treatmnet of palindrome cases, where the edge consists only
                // one kmer, and both orientation of the kmers are the same and should
                // all be considered!
                if ( isPalindrome ) {
                    Bool rc = loc2.rc ^ loc1.rc ^ 1;
                    int stop2p = - loc2.start + 1;
                    int start2p = - stop2 + 1;
                    int offset2 = ( loc1.rc ? stop1 - stop2p : start2p - loc1.start );
                    simple_align_data a(loc1.rid, loc2.rid, offset2, rc);
                    uniq_aligns.insert( simple_align_data(loc1.rid, loc2.rid, offset2, rc) );
                    //if ( ! ValidateAlign(a) ) {
                    //    #pragma omp critical 
                    //    {
                    //        cout << "Could not validate alignment " << a.rc2 << endl;
                    //        cout << "read1 on " << loc1.start << "," << stop1
                    //             << " read2(palindrom) on " << loc2.start << "," << stop2 
                    //             << " reverted to "         << start2p << "," << stop2p 
                    //             << endl;
                    //        int ulen = mpGraph->getEdge( loc1.uid ).getLength();
                    //        cout << "ulen= " << ulen << endl;
                    //    }
                    //}
                }
            }
        }
        pFriends->clear();
        int n_false_align = 0;
        for( std::set<simple_align_data>::iterator it = uniq_aligns.begin(), end = uniq_aligns.end();
                it != end; it++ ) {
            if ( ValidateAlign( *it ) )
                pFriends->push_back( Friend(it->id2,it->offset,it->rc2) );
            else
                n_false_align++;
        }
    }

    // Pathing provided the read head on the unipaths graph
    void PathOneRead ( size_t read_id, PathSegVec *loc_vec ) const {
        set<ReadLocOnUnipath> locs;
        const bvec& read = mReads[ read_id ];
        int readLen = read.size();
        int nkmers = readLen - K + 1;
        if ( nkmers < 0 ) return; 
        // pathing
        for ( int rpos = 0; rpos < nkmers; rpos++ ) {
            KMer<K> kmer( read.begin() + rpos );
            KDef const* pDef = mpDict->lookup(kmer);
            if ( !pDef ) { continue; }
            EdgeID edgeID = pDef->getEdgeID();
            const UnipathEdge *pEdge = &mpGraph->getEdge(edgeID);
            KmerID kmerID = pEdge->getKmerID( pDef->getEdgeOffset() );
            bool rc = IsRC( kmer,kmerID );
            // number of skipped bases from the unipath 
            int skipped = kmerID.val() - pEdge->getInitialKmerID().val();
            short ustart = ( rc ? skipped - (nkmers-1 - rpos) : skipped - rpos );
            ReadLocOnUnipath the_loc =
                { edgeID, ustart, static_cast<unsigned>(read_id), rc };
            locs.insert( the_loc );
        } 
        (*loc_vec).assign( locs.begin(), locs.end() );
    }

    int DownSampleLocsOfOneRead( PathSegVec *loc_vec ) const {
        size_t nsegs = (*loc_vec).size();
        if ( nsegs < 2 ) return 0;

        int nkmers = mReads[ (*loc_vec)[0].rid ].size() -K + 1;

        vec< pair<int,int> > seg_coverage(nsegs);
        vec<int> seg_lens(nsegs);
        for ( size_t i = 0; i < nsegs; ++i ) {
            int ulen = mpGraph->getEdge( (*loc_vec)[i].uid ).getLength();
            int rstart = -1, rstop = -1;
            if ( ! (*loc_vec)[i].rc ) {
                rstart = max( -(*loc_vec)[i].start, 0 );
                rstop = min( rstart + ulen , nkmers );
            }
            else {
                rstart = max( (*loc_vec)[i].start + nkmers - ulen, 0 );
                rstop = min( rstart + ulen , nkmers );
            }
            seg_coverage[i] = make_pair(rstart, rstop);
            seg_lens[i] = ulen;
        }

       // Select the segments, starting from the largest until every 10-base
       // division in the read has enough coverage. Long unipaths are always
       // kept.
       const int kDivisionSize = 10;
       const int kTargetDivCoverage = 1;
       const int kGoodUnipathLen = 5;  

       vec<Bool> todel( nsegs, true);
       vec<int> seg_indices( nsegs, vec<int>::IDENTITY );
       ReverseSortSync( seg_lens, seg_indices ); 

       vec<int> times_covered( (nkmers-1)/kDivisionSize + 1, 0);
       bool ignore_tail_division = ( times_covered.size() * kDivisionSize - nkmers < 10 ) ;
       for ( size_t i = 0; i < nsegs; ++i ) {
           size_t seg_index = seg_indices[i];
           //// discard redundant segments ( this division it covers all has enough segments )
           //bool is_redundant = true; 
           //for( int j =  seg_coverage[seg_index].first / kDivisionSize; 
           //        j <= (seg_coverage[seg_index].second-1) / kDivisionSize; ++j ) 
           //    if ( times_covered[j] < kTargetDivCoverage ) {
           //        is_redundant = false;
           //        break;
           //    }
           //if ( seg_lens[i] < kGoodUnipathLen && ! is_redundant 
           //        || seg_lens[i] >= kGoodUnipathLen ) {
               for( int j =  seg_coverage[seg_index].first / kDivisionSize; 
                       j <= (seg_coverage[seg_index].second-1) / kDivisionSize; ++j ) 
                   times_covered[j]++;
               todel[seg_index] = false;
           //}
           // Are all divisions covered?
           bool is_well_covered = true;
           size_t div_end = ( ignore_tail_division ? times_covered.size() -1 : times_covered.size() );
           for ( size_t j = 0; j < div_end; ++j ) {
               if ( times_covered[j] < kTargetDivCoverage) {
                   is_well_covered = false;
                   break;
               }
           }
           // exit conditions
           if ( is_well_covered && seg_lens[i] < kGoodUnipathLen ) { break; }
       }
       // return values
       EraseIf( *loc_vec, todel );
       return  nsegs - (*loc_vec).size();
    }

    bool IsRC( KMer<K> const& kmer, KmerID const& kmerID ) const { 
        using std::equal;
        HugeBVec::const_iterator seqItr( mpGraph->getBases( kmerID ) );
        bool result = ! equal( kmer.begin(),kmer.end(),seqItr );
        Assert( !result || equal( kmer.rcbegin(),kmer.rcend(),seqItr ) );
        return result; 
    }

    void ReportMemUsage() {
        cout << Date() << ": Peak memory use = " 
            << PeakMemUsageBytes( ) / 1000000000.0 << resetiosflags(ios::fixed)
            << " GB" << endl;
    }

private:
    const vecbvec          &mReads;
    KmerDict<K>            *mpDict;
    UnipathGraph<K>        *mpGraph;
    vec<ReadULocVec>        mULocs;                  // read path seg on unipaths, indexed by unipaths id
                                                     // sorted by the starting positon on unipath
    int                     mCopyNumberMax;          // Ignore short unipath with high copy number
    Bool                    mUseDownSampling;
    int                     mVerbose;
};

#endif
