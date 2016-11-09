///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>

#include <algorithm>

#include "CoreTools.h"
#include "PrintAlignment.h"
#include "pairwise_aligners/SmithWatBandedA.h"
#include "pairwise_aligners/SmithWatFree.h"
#include "paths/long/CreateGenome.h"
#include "paths/long/EvalCorrected.h"
#include "paths/long/Friends.h"
#include "paths/long/KmerAlign.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/MakeKmerStuff.h"
#include "paths/long/OverlapReads.h"
#include "paths/long/ultra/Prefab.h"
#include "random/Random.h"

//#define PREFAB_EXT_PARA            // allow external parameters to be fed in

void SelectInitialReads( const vec<int>& rid, const IAndOsVec& F, vec<int>& rid1,
     const long_heuristics& heur, const long_logging_control& log_control )
{    srandomx(1234567);
     for ( int i = 0; i < rid.isize( ); i++ )
          if ( randomx( ) % 5 == 0 ) rid1.push_back( rid[i] );    }

void CorrectSomeMoreReads0( const vec<int>& rid, const vec<int>& rid1,
     const IAndOsVec& F, const VecEFasta& corrected, const vec<int>& cid,
     VecEFasta& corrected2, vec<int>& cid2, const long_heuristics& heur,
     const long_logging_control& log_control )
{    return;    }


namespace {  // anonymous namespace for local functions

    // Local data structure to store the mapping of the founder reads to candidate assemblies
    struct AlignRange {
        align a;
        int pos1, Pos1; // positin on the read
        int pos2, Pos2; // position on the assembly
        int n_errors;   // total number of errors
        float error_rate;   // error rate
    };
    bool operator< (const AlignRange& lhs, const AlignRange& rhs) { return lhs.error_rate < rhs.error_rate; }
    ostream& operator<< ( ostream& out, const AlignRange& r) {
        return out << "[" << r.pos1 << ", " << r.Pos1 << ") to [" << r.pos2 << ", " << r.Pos2 << ") n_error= " 
                   << r.n_errors << "(" << r.error_rate << ")";   }

    // Compare only the second element of two pairs
    struct CompSecond {
        bool operator() ( const pair<int,int>& l, const pair<int,int>& r ) {
            return l.second < r.second;
        }
    };

    // Select the seemingly best assembly based on the alignment results
    // ,which is the first qualifying assembly it finds after sort all 
    // assemblies by error rates.
    int SelectAssembly( const vec<AlignRange>& results, double Max_Err_Rate,
           int Min_Assembly_Lengh, const vec<basevector>& assemblies,
           vec<int>* p_alt_sels = 0 ) 
    {
        if ( p_alt_sels != 0 ) p_alt_sels->clear();
        vec<AlignRange> rr = results;
        vec<int> index( results.size(), vec<int>::IDENTITY );
        vec<int> lens( assemblies.size() );
        for ( size_t i = 0; i < assemblies.size(); ++i ) lens[i] = assemblies[i].size();
        SortSync(rr, index, lens);
        int sel = -1;
        for ( size_t i = 0; i < rr.size(); ++i ) {
            //if ( rr[i].pos2 == 0 || rr[i].Pos2 == lens[i] ) continue; // require full range mapping
            if ( rr[i].error_rate > Max_Err_Rate ) continue;
            if ( (rr[i].Pos1 - rr[i].pos1) < Min_Assembly_Lengh ) continue;
            if ( sel == -1 ) sel = index[i];       
            if ( p_alt_sels != 0 ) p_alt_sels->push_back( index[i] );
        }
        return sel;
    }

    // location of a read on the genome
    struct ReadGLoc {
        size_t start, stop, gid;
    };

    ReadGLoc SimRange( const long_logging_control& log_control, int id ) {
        ReadGLoc ans = { (size_t) (*log_control.readlocs)[id].start, 
               (size_t) (*log_control.readlocs)[id].stop, 
               (size_t) (*log_control.readlocs)[id].id };
        return ans;
    }
    basevector ReadTruth( const long_logging_control& log_control, int id) {
        ReadGLoc loc = SimRange( log_control, id );
        return basevector( (*log_control.G)[loc.gid], loc.start, loc.stop - loc.start );
    }
    ostream& operator<< ( ostream& out, const ReadGLoc& gloc ) {
        return out << "[" << gloc.start << "," << gloc.stop << ")@" << gloc.gid;
    }
    void PrintPathGLocs( ostream& out, const vec<ReadOverlap>& path, 
            const vec<int>& creads_filt_ids, const long_logging_control& log_control ) {
        for ( size_t i = 0; i < path.size(); ++i ) {
            if ( i==0 ) 
                out <<  SimRange(log_control,  creads_filt_ids[ path[i].rid1 ] ) << " ";
            out <<  SimRange(log_control, creads_filt_ids[ path[i].rid2 ] ) << " ";
        }
    }

     bool IsTrueFriend( const long_logging_control& log_control, int id1, int id2 )
     {    const vec<ref_loc>& locs = *log_control.readlocs;
          return ( locs[id1].id == locs[id2].id
               && IntervalOverlap( locs[id1].start, locs[id1].stop,
                    locs[id2].start, locs[id2].stop ) > 0 );    }

     // Quickly check if there are errors in the read

    int HasError( const long_logging_control& log_control, int id, const vec<basevector>& reads ) {
        basevector read_truth = ReadTruth( log_control, id );
        for ( size_t i = 0; i < reads.size(); ++i ) 
            if ( search( read_truth.begin(), read_truth.end(), reads[i].begin(), reads[i].end() ) != read_truth.end() ) 
                return 0;
        return 1;
    }
    // Free Smith-Waterman alignment -- to be discarded
    // There is an undesirable behavior that force the ends of the query sequence
    // to be included in the alignment.
    AlignRange SWFAlign( const basevector& s, const basevector& t) 
    {
        align a;
        SmithWatFreeSym( s, t, a, false, false, 1, 1, 0 );
        int error = a.Errors( s, t );
        float error_rate = error * 1.0 / ( a.Pos1() - a.pos1() );
        AlignRange result = { a, a.pos1(), a.Pos1(), a.pos2(), a.Pos2(), error, error_rate };
        return result;
    }

    // Banded Smith-Waterman alignment 
    AlignRange SWAlign( const basevector& s, const basevector& t) 
    {
        align a;
        int error = 0;
        int extra = s.size() / 10;
        int off_ub = 0 + extra;
        int off_lb = - ( (int)t.size() - (int)s.size() ) - extra;
        int offset = ( off_ub + off_lb ) / 2;
        int bandwidth = ( off_ub - off_lb ) / 2;
        SmithWatBandedA( s, t, offset, bandwidth, a, error, 0, 1, 1 );
        float error_rate = error * 1.0 / ( a.Pos1() - a.pos1() );
        AlignRange result = { a, a.pos1(), a.Pos1(), a.pos2(), a.Pos2(), error, error_rate };
        return result;
    }

    void GetKmersAndSort( const int K,  const basevector& u, 
                          vec<basevector>& kmers, vec<int> & locs )
    {   kmers.resize( u.size() - K + 1 );
        locs.resize( u.size() - K + 1 );
        for ( int j = 0; j <= u.isize( ) - K; j++ )
        {   kmers[j].SetToSubOf( u, j, K );
            locs[j] = j;  } 
        SortSync(kmers, locs);  }

    // Banded Smith-Waterman alignment assisted by kmer align ( for offset and bandwidth calculation )
    AlignRange KmerAlign( const int K, const basevector& s, const basevector t, 
                          const vec<basevector>& kmers, const vec<int> & locs )
    {
        vec< pair<int,int> > offsets;
        for ( size_t i = 0; i < t.size()-K+1; ++i ) {
            basevector tk( t, i, K);
            vec<basevector>::const_iterator it 
                = lower_bound( kmers.begin(), kmers.end(), tk );
            while ( it != kmers.end() && *it == tk ) {
                offsets.push( locs[it - kmers.begin()] - i, 
                             locs[it - kmers.begin()] );
                ++it;
            }
        }
        AlignRange result0 = { align(), 0, 0, 0, 0, 100, 100.0 };
        if ( offsets.empty() ) return result0;
        vec< pair<int,int> > offsets2;
        LargestOffsetCluster( offsets, offsets2, 200, 0.2 );
        {  // kmeralign-assisted SmithWaterman
            Sort( offsets2 );
            int off_u = offsets2.back().first;
            int off_l = offsets2.front().first;
            int off = ( off_u + off_l ) / 2;
            int bandwidth = ( off_u - off_l ) / 2;
            align a;
            int error;
            SmithWatBandedA( s, t, off, bandwidth, a, error, 0, 1, 1 );
            float error_rate = error * 1.0 / ( a.Pos1() - a.pos1() );
            AlignRange result = { a, a.pos1(), a.Pos1(), a.pos2(), a.Pos2(), error, error_rate };
            return result;
        }
    }

    template<int K> 
    void MakeKmerLookup0Single( const vecbasevector& unibases, vec< triple<kmer<K>,int,int> >& kmers_plus, 
         const long_logging_control& log_control, const long_logging& logc, ostream& out )
    {    vec<int64_t> starts;
         starts.reserve( unibases.size() + 1 );
         starts.push_back(0);
         for ( size_t i = 0; i < unibases.size( ); i++ )
         {    const basevector& u = unibases[i];
              starts.push_back( starts.back( ) + Max( 0, u.isize( ) - K + 1 ) );    }
         kmers_plus.resize( starts.back( ) );
         for ( size_t i = 0; i < unibases.size( ); i++ )
         {    const basevector& u = unibases[i];
              for ( int j = 0; j <= u.isize( ) - K; j++ )
              {    int64_t r = starts[i] + j;
                   kmers_plus[r].first.SetToSubOf( u, j );
                   kmers_plus[r].second = i; 
                   kmers_plus[r].third = j;    }    }
         Sort(kmers_plus);
    }

    // Mark those friend reads that have at least half of kmer match of the most matching one.
    void FilterFriends( const vec< vec< pair<int,int> > >& offsets, vec<Bool>& accepted ) {
        vec<int> pos1_count( accepted.size() );
        for ( size_t id = 1; id < accepted.size(); id++ ) {    
            vec<int> pos1;
            for ( int j = 0; j < offsets[id].isize( ); j++ )
                pos1.push_back( offsets[id][j].second );
            UniqueSort(pos1);
            pos1_count[id] = pos1.size( );    
        }
        int max_pos1_count = 0;
        for ( size_t id = 0; id < accepted.size(); id++ ) {
            if ( id != 0 ) max_pos1_count = Max( max_pos1_count, pos1_count[id] );    
            if ( pos1_count[id] < max_pos1_count/2 ) continue;
            if ( offsets[id].empty( ) ) continue; 
            accepted[id] = True;
        }
    }

    // Create the list of kmer offsets between reads
    template<int K>
    void CreateOffsets( const vec< triple<kmer<K>,int,int> >& kmers_plus,  
            vec< vec< pair<int,int> > >& offsets ) 
    {
        for ( int64_t i = 0; i < (int64_t) kmers_plus.size( ); i++ ) {    
            int64_t j, z1, z2;
            for ( j = i + 1; j < (int64_t) kmers_plus.size( ); j++ )
                if ( kmers_plus[j].first != kmers_plus[i].first ) break;
            for ( z1 = i; z1 < j; z1++ )
                if ( kmers_plus[z1].second == 0 ) break;
            for ( z2 = z1; z2 < j; z2++ )
                if ( kmers_plus[z2].second != 0 ) break;
            for ( int64_t z = z1; z < z2; z++ ) {    
                int pos1 = kmers_plus[z].third;
                for ( int64_t k = i; k < j; k++ ) {    
                    int id2 = kmers_plus[k].second; 
                    if ( id2 == 0 ) continue;
                    int pos2 = kmers_plus[k].third; // position on read id2
                    offsets[id2].push( pos1-pos2, pos1 );    
                }    
            }
            i = j - 1;    
        }
        for ( size_t id = 1; id < offsets.size(); id++ ) Sort( offsets[id] );
    }

    // If we have several candidate corrected reads, which one is supported more by all
    // friend reads
    void VoteByReads( const vec<basevector>& candidates, const vec<basevector>& reads, 
         vec<int>& votes ) 
    {
        const int K2=14;
        votes.assign( candidates.size(), 0 );
        // kmer-lookup for each candidates
        vec< vec<basevector> > kmers( candidates.size() ); 
        vec< vec<int> >locs( candidates.size() );
        for ( size_t i = 0; i < candidates.size(); ++i ) 
            GetKmersAndSort( K2, candidates[i], kmers[i], locs[i] );
        for ( size_t i = 0; i < reads.size(); ++i ) {
            vec< AlignRange > results;
            vec<int> index;
            for ( size_t j = 0; j < candidates.size(); ++j ) {
                AlignRange ar = KmerAlign( K2, candidates[j], reads[i],
                        kmers[j], locs[j] );
                results.push_back( ar );
                index.push_back(j);
            }
            SortSync(results, index);
            if ( results[0].error_rate > 0.1 ) continue;
            votes[index[0]]++;
        }
    }

    // Given the gang of raw friend reads (threads), find all suggested
    // edits to the founder read (t)
    int FindEdits( const vec<basevector>& threads, const basevector& t, 
         vec< vec<edit0> >& edits, ostream& out ) 
    {   
        edits.clear_and_resize( t.size( ) + 1 );
        const int K2=14;
        vec<basevector> kmers; vec<int> locs;
        GetKmersAndSort( K2, t, kmers, locs );
        out << "Align threads to read " << endl;
        int n_good_reads = 0;
        for ( int i = 0; i < threads.isize( ); i++ ) {    
            AlignRange ar = KmerAlign( K2, t, threads[i], kmers, locs );
            out << i << " " << ar << endl;
            if ( ar.error_rate > 0.1 ) continue;
            n_good_reads++;
            const align& a = ar.a;
            int p1 = a.pos1( ), p2 = a.pos2( );
            for ( int j = 0; j < a.Nblocks( ); j++ ) {    
                if ( a.Gaps(j) > 0 ) {    
                    basevector b( threads[i], p2, a.Gaps(j) );
                    edits[p1].push( INSERTION, b.ToString( ) );
                    p2 += a.Gaps(j);    
                }
                if ( a.Gaps(j) < 0 ) {
                    edits[p1].push( DELETION, -a.Gaps(j) ); 
                    p1 -= a.Gaps(j);    
                }
                for ( int x = 0; x < a.Lengths(j); x++ ) {
                    if ( threads[i][p2] != t[p1] )
                        edits[p1].push( SUBSTITUTION, (char) as_base( threads[i][p2] ) );
                        p1++, p2++;    
                   }    
              }    
         }
         for ( int p = 0; p <= t.isize( ); p++ )
             Sort( edits[p] );
         return n_good_reads;    
    }
    
    // Given all the edits generated by FindEdits, identify signals suggesting
    // errors in the founder reads, and make the correction.
    int  MakeEdits( const vec<basevector>& threads, basevector& t, 
         const vec< vec<edit0> >& edits )
    {    const int Threshold_Vote = threads.isize( ) * 0.5;
         vec< pair<int,edit0> > keepers;
         for ( int p = 0; p <= t.isize( ); p++ )
         {    for ( int i = 0; i < edits[p].isize( ); i++ )
              {    int j = edits[p].NextDiff(i);
                   if ( (j-i) > Threshold_Vote ) keepers.push( p, edits[p][i] );
                   i = j - 1;    }    }
         int nedits = 0;
         for ( int j = keepers.isize( ) - 1; j >= 0; j-- )
         {    if ( j < keepers.isize( ) - 1 && keepers[j].first == keepers[j+1].first )
                   continue;
              nedits++;
              int p = keepers[j].first;
              const edit0& e = keepers[j].second;
              if ( e.etype == INSERTION )
              {    basevector b1( t, 0, p );
                   basevector b2( e.seq );
                   basevector b3( t, p, t.size( ) - p );
                   t = Cat( b1, b2, b3 );    }
              if ( e.etype == DELETION )
              {    basevector b1( t, 0, p );
                   basevector b2( t, p + e.n, t.isize( ) - ( p + e.n ) );
                   t = Cat( b1, b2 );    }
              if ( e.etype == SUBSTITUTION ) t.Set( p, as_char( e.seq[0] ) );    
         } 
         return nedits;
    }

    // Find the relative offset and uncertainty(bandwidth) between two reads given the kmer offset lists.
    void OffsetAndBandwidth( const vec< pair<int,int> >& offsets, int& offset, int& bandwidth ) {
        // Form offsets into groups, breaking whenever there is a
        // > max_sep = 20 separation in offsets.
        vec<int> ostarts;
        ostarts.push_back(0);
        const int max_sep = 20;
        for ( int j = 0; j < offsets.isize( ) - 1; j++ ) 
            if ( offsets[j+1].first > offsets[j].first + 20 )
                ostarts.push_back(j+1);
        ostarts.push_back( offsets.size( ) );
        // Compute the size of each group, as measured by its number
        // of distinct rpos2 values.  Find the largest group.
        int gp_max = 0, gp_best = -1;
        for ( int j = 0; j < ostarts.isize( ) - 1; j++ ) {    
            vec<int> rpos2;
            for ( int l = ostarts[j]; l < ostarts[j+1]; l++ )
                rpos2.push_back( offsets[l].second );
            UniqueSort(rpos2);
            if ( rpos2.isize( ) > gp_max )
            {   gp_max = rpos2.size( );
                gp_best = j;    }    
        }
        // Set offset to the median within the winning group, and
        // set bandwidth to span the group, but no more than 
        // max_bandwidth = 200;
        int start = ostarts[gp_best], stop = ostarts[gp_best+1];
        int mid = start + (stop-start)/2;
        offset = offsets[mid].first;
        const int bw_add = 12;
        bandwidth = Max( offset - offsets[start].first,
                offsets[stop-1].first - offset ) + bw_add;
        const int max_bandwidth = 200;
        //bandwidth = Min( bandwidth, max_bandwidth );
        bandwidth += abs(offset) * 1.1;  // extra bandwidth
    }

    // Find the relative offset and uncertainty(bandwidth) between two reads given the kmer offset lists.
    // Faster (less accurate) method
    void OffsetAndBandwidth2( const vec< pair<int,int> >& offsets, int& offset, int& bandwidth ) {
        const int Extra_Bandwidth = 50;
        vec< pair<int,int> > offsets_filtered;
        GroupAndSelect( offsets, offsets_filtered, 200, 0.2, 0 );
        if ( offsets_filtered.empty() ) {
            offset = 0;
            bandwidth = -1;
            return;
        }
        vec<int> all( offsets_filtered.size() );
        for ( size_t i = 0; i < all.size(); ++i ) 
            all[i] = offsets_filtered[i].first;
        Sort(all);
        offset = Median(all);
        bandwidth = all.back() - all.front();
        bandwidth += Extra_Bandwidth;
    }

    // Sometimes efasta contains {,}, creating redundant edges when expanding. 
    // We want to remove those empty brackets.
    efasta CleanEfasta( const efasta& input) {
        efasta output;
        output.reserve( input.size() );
        efasta::const_iterator it = input.begin();
        efasta::const_iterator it2;
        while( it!= input.end() ) {
            if ( *it == '{' ) {
                // point it2 to '}'
                it2 = it+1;
                bool is_good = false;
                while( *it2 != '}' ) {
                    if ( *it2 != ',' ) is_good = true;
                    it2++;
                }
                if ( is_good ) output.append( it, it2+1 );
                it = it2;
            } 
            else output.push_back(*it);
            ++it;
        }
        return output;
    }
    enum CSMR_MESSAGE { CSMR_OK, CSMR_NO_FRIENDS, CSMR_NO_CR_FRIENDS, CSMR_NO_ASSEMBLY, CSMR_MANY_ASSEMBLY, CSMR_SHORT_ASSEMBLY } ;
}


pair<efasta,CSMR_MESSAGE> CorrectOneReadFast( int ec, const vecbasevector& reads, 
        const IAndOsVec& F,  const vec<int>& known_friends,
        const VecEFasta& corrected, const vec<int>& cid,
        const vec<bool>& is_corrected, const vec<int>& mapping,
        const long_heuristics& heur, const long_logging_control& log_control,
        int VERBOSITY,
        ostream& out ) 
{
#ifdef PREFAB_EXT_PARA
    extern int Csmr_Edits;
    extern int Known_Friends;
    extern int Cheat_Reads;          // Cheating: replace the read-to-be-corrected with true sequence
    extern int Cheat_Friend_Reads;   // Cheating: replace the pre-corrected reads with true sequence
    extern int Cheat_Friend_List;    // Cheating: the friend list is based their true genome locations
#else
    const int Csmr_Edits = 0;         // Make further edits using piled reads
    const int Known_Friends = 1;      // Use only friends known during pre-correction
    const int Cheat_Reads = 0;        // Cheating: replace the read-to-be-corrected with true sequence
    const int Cheat_Friend_Reads = 0; // Cheating: replace the pre-corrected reads with true sequence
    const int Cheat_Friend_List = 0;  // Cheating: the friend list is based their true genome locations
#endif

    // where is the real sequence of the read

    const vec<ref_loc>& rlocs = *log_control.readlocs;
    int sim_gid = rlocs[ec].id; 
    int sim_start = rlocs[ec].start, sim_stop = rlocs[ec].stop;
    Bool sim_rc = rlocs[ec].rc2;
    basevector read_truth( 
         (*log_control.G)[sim_gid], sim_start, sim_stop - sim_start );
    if ( sim_rc ) read_truth.ReverseComplement();
    // The read to be corrected
    basevector read_input( reads[ec] );
    const basevector& read0 = ( Cheat_Reads ? read_truth : read_input );
    // friends of the read
    vec<int> friends;
    vec<bool> friend_rcs;
    if ( Cheat_Friend_List ) {
        for ( size_t fid = 0; fid < reads.size(); ++fid ){
            if ( (int)fid == ec ) continue;
            if( IsTrueFriend( log_control, ec, fid ) ) {
                friends.push_back(fid);
                Bool sim_rc2 = rlocs[fid].rc2;
                friend_rcs.push_back(sim_rc ^ sim_rc2);
            }
        }
    }
    else {
        for ( size_t i = 0; i < known_friends.size(); ++i ) {
            bool rc = ( known_friends[i] < 0 );
            int fid = ( rc ? ~known_friends[i] : known_friends[i] );
            friends.push_back(fid);
            friend_rcs.push_back(rc);
        }
    }
    if ( VERBOSITY >=2 ) {
        out << "Make corrections using " << friends.size() << " friend reads " << endl;
    }
    if ( friends.size() < 2 ) {
        if ( VERBOSITY >=1 ) 
            out << "Cannot find at least two friends" << endl;
        return make_pair( efasta(), CSMR_NO_FRIENDS );
    }
    // Generate the assembly of the only the corrected friend reads
    vec< vec<basevector> > creads_filtered;           // expanded efasta 
    vec< vec<Ambiguity> > creads_ambs;                // ambiguities from the efasta
    vec< int > creads_filt_ids;
    for ( size_t i = 0; i < friends.size(); ++i ) {
        int fid = friends[i];
        if ( ! is_corrected[ fid ] ) continue;
        creads_filt_ids.push_back( fid );
        efasta cread;
        if ( Cheat_Friend_Reads ) {
            int sim_start2 = rlocs[fid].start, sim_stop2 = rlocs[fid].stop;
            int sim_gid2 = rlocs[fid].id;
            Bool sim_rc2 = rlocs[fid].rc2;
            ForceAssertEq( sim_gid2, sim_gid );
            basevector read_truth2( 
                 (*log_control.G)[sim_gid2],  sim_start2, sim_stop2 - sim_start2 );
            if ( sim_rc2 ) read_truth2.ReverseComplement();
            cread = efasta( read_truth2.ToString() );
        }
        else
            cread = corrected[ mapping[fid] ] ;
        if (  friend_rcs[i] ) cread.ReverseComplement();
        creads_filtered.push_back( vec<basevector>() );
        cread.ExpandTo( creads_filtered.back() );
        fastavector fa;
        creads_ambs.resize( creads_ambs.size() + 1 );
        cread.FlattenTo( fa, creads_ambs.back(), False ); // extract the ambiguities
        //if ( creads_ambs.back().size() > 0 ) {
        //    out << "read " << creads_ambs.size()-1 << endl;
        //    out << cread << endl;
        //    out << corrected[ mapping[fid] ] << endl;
        //}
    }
    // Validation: are the friends correct?
    if ( VERBOSITY >= 1 ) {
        unsigned int n_true_friends = 0;
        for ( size_t i = 0; i < creads_filt_ids.size(); ++i ) 
            if( IsTrueFriend( log_control, ec, creads_filt_ids[i] ) ) n_true_friends++;
        if ( n_true_friends < creads_filt_ids.size() && VERBOSITY >=1 ) {
            out << "Wrong friending found at read " << ec << " length= " << read0.size() << endl;
            out << "n_friends= " << creads_filt_ids.size() 
                << ", n_true_friends= " << n_true_friends 
                << endl;
        }
        vec<int> n_errors( creads_filt_ids.size(), 0 ), n_ambs( creads_filt_ids.size(), 0 );
        for ( size_t i = 0; i < creads_filt_ids.size(); ++i ) {
            n_errors[i] = HasError( log_control, creads_filt_ids[i], creads_filtered[i] );
            n_ambs[i] = corrected[ mapping[ creads_filt_ids[i] ] ].AmbEventCount();
        }
        int n_error_reads = n_errors.size() - n_errors.CountValue(0);
        if ( n_error_reads > 0) {
            out << n_error_reads << " out of " << n_errors.size() << " reads has errors: ";
            for ( size_t i = 0; i < n_errors.size(); ++i ) 
                if ( n_errors[i] > 0 ) out << i << " ";
            out << endl;
        }
        int n_amb_reads = n_ambs.size() - n_ambs.CountValue(0);
        if ( n_amb_reads > 0) {
            out << n_amb_reads << " out of " << n_ambs.size() << " reads has ambiguities: ";
            for ( size_t i = 0; i < n_ambs.size(); ++i ) 
                if ( n_ambs[i] > 0 ) out << i << "(" << n_ambs[i] << ") ";
            out << endl;
        }
    }
    // More about the reads
    if ( VERBOSITY >= 3 ) {
        out << "Read0 start " << 0  << " to " << sim_stop - sim_start << "@" << sim_gid << endl;
        out << "Find " << creads_filt_ids.size() << " corrected friend reads" << endl;
        vec<int> starts;
        vec<int> gids;
        vec<int> index;
        for ( size_t i = 0; i < creads_filt_ids.size(); ++i ) {
            int fid = creads_filt_ids[i];
            int sim_start2 = rlocs[fid].start, sim_stop2 = rlocs[fid].stop;
            int sim_gid2 = rlocs[fid].id;
            Bool sim_rc2 = rlocs[fid].rc2;
            int start = sim_start2 - sim_start;
            starts.push( start );
            index.push( i );
            gids.push( sim_gid2 );
        }
        SortSync(starts, gids, index);
        for ( size_t i = 0; i < starts.size(); ++i ) 
            out << "    read " << index[i] << "_0 " << starts[i]
                << " to " << starts[i] + creads_filtered[ index[i] ][0].size() 
                << "@" << gids[i] << endl;
    }
    if ( creads_filtered.size() < 2 ) {
        if ( VERBOSITY >=1 ) 
            out << "Cannot find at least two corrected friends" << endl;
        return make_pair( efasta(), CSMR_NO_CR_FRIENDS );
    }
    // assembly of the reads
    const int Max_Num_Assembly = 5;
    SomeReadsAssembler assembler( creads_filtered );
    assembler.SetMinOverlap(800);
    assembler.SetMaxAssemblyNum( Max_Num_Assembly );
    assembler.SetCombineAmbReads( creads_ambs );
    assembler.InitGraph();
    assembler.Assembly();
    const vec<basevector>& assemblies = assembler.GetAssemblies();
    if ( assemblies.empty() ) {
        if ( VERBOSITY >=1 ) 
            out << "Cannot generate any assembly" << endl;
        return make_pair( efasta(), CSMR_NO_ASSEMBLY );
    }
    if ( assemblies.isize() >= Max_Num_Assembly ) {
        if ( VERBOSITY >=1 ) 
            out << "Generated too many assemblies" << endl;
        return make_pair( efasta(), CSMR_MANY_ASSEMBLY );
    }
    // kmer-lookup of read0
    const int K2=14;
    vec<basevector> kmers; 
    vec<int> locs;
    GetKmersAndSort( K2, read0, kmers, locs );
    vec<basevector> kmers0; vec<int> locs0;
    if ( VERBOSITY >= 1 )
        GetKmersAndSort( K2, read_truth, kmers0, locs0 );
    vec< AlignRange > results;
    vec< AlignRange > results0; // validation
    for ( size_t i = 0; i < assemblies.size(); ++i ) {
        results.push_back( KmerAlign( K2, read0, assemblies[i], kmers, locs ) );
        //results.push_back( SWFAlign( read0, assemblies[i] ) );
        if ( VERBOSITY >= 1 )
            results0.push_back( KmerAlign( K2, read_truth, assemblies[i], kmers0, locs0 ) );
            //results0.push_back( SWFAlign( read_truth, assemblies[i] ) );
    }
    if ( VERBOSITY >=2 ) {
        if ( VERBOSITY >=3  ) 
            assembler.graph_.PrintGraph(out);
        for ( size_t i = 0; i < assemblies.size(); ++i ) {
            out << "assembly " << i << " (" << assemblies[i].size() << ") "
                << assembler.GetAssemblyPath(i) << endl;
            out << " / " << results[i]
                << " / " << results0[i] << endl; 
        }
    }
    int sel = SelectAssembly( results, 0.1, read0.size() * 0.8, assemblies );
    if ( sel == -1 ) {
        if ( VERBOSITY >= 1 ) 
            out << "Warning: read mapping on assembly is poor" << endl;
        return make_pair( efasta(), CSMR_SHORT_ASSEMBLY );
    }
    // Validate the selected assembly path: 
    // Do they cover the whole read and do they come from the same chromosome?
    if ( VERBOSITY >=1 ) {
        int start0 = sim_start;
        int end0 = sim_stop;
        vec<ReadOverlap> path = assembler.GetAssemblyPathRegion(sel, results[sel].pos2, results[sel].Pos2 );
        int rid1 = creads_filt_ids[ path.front().rid1 ];
        int rid2 = creads_filt_ids[ path.back().rid2 ];
        int start = (sim_rc ? SimRange(log_control, rid2 ).start
                            : SimRange(log_control, rid1 ).start );
        int end = (sim_rc ? SimRange(log_control, rid1 ).stop
                          : SimRange(log_control, rid2 ).stop );
        bool path_incomplete = ( start > start0 || end < end0 );
        bool path_invalid  = false;
        for ( size_t i = 0; i < path.size(); ++i ) {
            int rid1 = creads_filt_ids[ path[i].rid1 ];
            int rid2 = creads_filt_ids[ path[i].rid2 ];
            if ( i==0 && (int) SimRange( log_control, rid1 ).gid != sim_gid
                      || (int) SimRange( log_control, rid2 ).gid != sim_gid ) {
                path_invalid = true;
                break;
            }
        }
        if ( path_invalid || path_incomplete ) {
            out << "Read0 " << SimRange(log_control, ec) << endl;
            String ee = "";
            if ( path_incomplete ) ee += "incomplete";
            if ( path_invalid ) {
                if ( ! ee.empty() ) ee += "|";
                ee += "invalid";
            }
            out << "Full assembly " << sel << " size " << assemblies[sel].size() 
                << " " << ee 
                << " mapped range " << start - start0 << " " << end - start0 << endl;
            out << "Map  range: " ;
            PrintPathGLocs( out, path, creads_filt_ids, log_control );
            out << endl;
            for ( size_t i = 0; i < assemblies.size(); ++i ) {
                out << "Assembly " << i << " size " << assemblies[i].size() 
                    << " Path: " << assembler.GetAssemblyPath(i) << endl;
                out << " Locs: ";
                PrintPathGLocs( out, assembler.GetAssemblyPath(i), creads_filt_ids, log_control );
                out << endl;
            }
        }
    }
    const basevector& assembly = assemblies[sel];
    basevector provisional( assembly, results[sel].pos2, 
                           results[sel].Pos2 - results[sel].pos2 );
    if ( VERBOSITY >= 1 ) {
        AlignRange ar = KmerAlign( K2, read_truth, provisional, kmers0, locs0 );
        if ( ar.n_errors > 0 ) {
            out << "Error found in provisional "
                << ar << " from assembly " << sel << endl ;
            for ( size_t i = 0; i < assemblies.size(); ++i ) {
                out << "assembly " << i << " (" << assemblies[i].size()
                    << ") /" << results[i]
                    << " /" << results0[i] << endl;
            }
            if ( VERBOSITY >=2 ) {
                out << "read_truth: " << read_truth.ToString() << endl;
                out << "read_provi: " << provisional.ToString() << endl;
                PrintVisualAlignment( true, out, read_truth, provisional, ar.a  );
                out << "Check accuracy of corrected reads" << endl;
                vec< pair<int,int> > rid_alts;
                vec<ReadOverlap> path = assembler.GetAssemblyPathRegion(sel, results[sel].pos2, results[sel].Pos2 );
                for ( size_t i = 0; i < path.size(); ++i ) {
                    if ( i== 0 )
                        rid_alts.push( path[i].rid1, path[i].alt1 );
                    rid_alts.push( path[i].rid2, path[i].alt2 );
                }
                for ( size_t i = 0; i < rid_alts.size(); ++i ) {
                    basevector query2 = creads_filtered[ rid_alts[i].first ][ rid_alts[i].second ];
                    basevector read_truth2 = ReadTruth( log_control, creads_filt_ids[ rid_alts[i].first ] );
                    AlignRange ar2 = SWAlign(query2,read_truth2);
                    if ( ar2.n_errors > 0 ) {
                        out << "read " << rid_alts[i].first << "_" << rid_alts[i].second << endl;
                        PrintVisualAlignment( true, out, query2, read_truth2, ar2.a  );
                    }
                }
            }
        }
    }
    efasta provisional_efasta = assembler.GetEfastaAssembly(sel, results[sel].pos2, results[sel].Pos2);

    if ( Csmr_Edits == 0 )
        return make_pair( provisional_efasta, CSMR_OK );

    // Additional editing seems bring no added value to the error correction, given that
    // those errors have survived multiple rounds of error corrections.

    // collecting all the uncorrected friend reads
    vec<basevector> other_reads;
    other_reads.push_back( read0 );
    for ( size_t i = 0; i < friends.size(); ++i ) {
        int fid = friends[i];
        other_reads.push_back( reads[fid] );
        if ( friend_rcs[i] ) 
            other_reads.back().ReverseComplement();
    }
    if ( VERBOSITY >=2 ) 
        out << "Use " << other_reads.size() << " original reads for further correction " << endl;
    if ( other_reads.size() < 10 ) {
        if (VERBOSITY>=1 ) out << "Skipped extra editing due to insufficient reads " << other_reads.size() << endl;
        return make_pair( provisional_efasta, CSMR_OK );
    }

    // Editing method 1
    //    If there are alternative assemblies ( eg. from other chromosome ), use the original
    //    reads to vote and make the selection.
    //
    if ( Csmr_Edits == 1 ) { 
        // find alternative assemblies
        int match_size = results[sel].Pos1 - results[sel].pos1;
        vec<int> alt_sels; // other possible assemblies
        SelectAssembly( results, 0.1, match_size * 0.9, assemblies, &alt_sels );
        vec<basevector> candidates;
        vec<int> index;
        for ( size_t i = 0; i < alt_sels.size(); ++i ) {
            int aid = alt_sels[i];
            candidates.push_back( basevector(assemblies[aid], results[aid].pos2, 
                           results[aid].Pos2 - results[aid].pos2) );
            index.push_back(aid);
        }
        UniqueSortSync(candidates, index);
        if ( candidates.solo() ) return make_pair( provisional_efasta, CSMR_OK );

        if ( VERBOSITY >= 2 )
            out << "Check alternative = " << alt_sels.size() << endl;
        // voting
        vec<int> votes( candidates.size(), 0 );
        VoteByReads( candidates, other_reads, votes );
        ReverseSortSync(votes, candidates, index);
        if ( VERBOSITY >=2 ) {
            out << "Votes " << endl;
            for ( size_t i = 0; i < votes.size(); ++i ) 
                out << "candidate " << i << " " << votes[i] 
                    << ( candidates[i] == provisional ? "*": "" ) 
                    << endl;
        }
        if ( VERBOSITY >= 1 ) {
            AlignRange ar = KmerAlign( K2, read_truth, candidates[0], kmers0, locs0 );
            if ( ar.n_errors > 0 ) {
                out << "Error found in corrected " << ar << endl;
                if ( VERBOSITY >=2 ) 
                    PrintVisualAlignment( true, out, read_truth, candidates[0], ar.a  );
            }
        }
        int sel2 = index[0];
        efasta edited_efasta = assembler.GetEfastaAssembly(sel2, results[sel2].pos2, results[sel2].Pos2);
        return make_pair( edited_efasta, CSMR_OK ); 
    }

    // Editing method 2:
    //    Pile up all the friend reads on the provisional one, and make edits if more than half reads suggest so.
    //
    if ( Csmr_Edits == 2) {
        basevector to_correct = provisional;
        vec< vec<edit0> > edits;
        int n_good_reads = FindEdits( other_reads, to_correct, edits, out );
        if ( VERBOSITY >= 2 ) {
            out << "#good reads" << n_good_reads << endl;
            for ( int p = 0; p < edits.isize( ); p++ )
                if (  edits[p].size() >=1 )
                    out << "site " << p << " has " << edits[p].size() << endl;
        }
        int nedits = MakeEdits( other_reads, to_correct, edits );
        if ( nedits == 0 ) 
            return make_pair( provisional_efasta, CSMR_OK );
        if ( VERBOSITY >= 1 && nedits >0 ) {
            out << "Made " << nedits << "edits" << endl;
            AlignRange ar = KmerAlign( K2, read_truth, to_correct, kmers0, locs0 );
            if ( ar.n_errors > 0 ) {
                out << "Error found in corrected " << ar << endl;
                if ( VERBOSITY >=2 ) 
                    PrintVisualAlignment( true, out, read_truth, to_correct, ar.a  );
            }
        }
        return make_pair( efasta(to_correct), CSMR_OK ); 
    }

    return make_pair( provisional_efasta, CSMR_OK );
}

void CorrectSomeMoreReads( 
     const vecbasevector& reads,            // reads
     const vec<int>& rid, const vec<int>& rid1, const vec< vec<int> >& friend_ids1,
     const IAndOsVec& F, const VecEFasta& corrected, const vec<int>& cid,
     VecEFasta& corrected2, vec<int>& cid2, const long_heuristics& heur,
     const long_logging_control& log_control, const long_logging& logc,
     const ref_data& ref )
{    
#ifdef PREFAB_EXT_PARA
    extern int Csmr_Rid;
    extern int Csmr_Verb;
    const int VERBOSITY = Csmr_Verb;
#else
    const int VERBOSITY = 0;
    const int Csmr_Rid = -1;
#endif
    // parameters for EvalCorrected()
    const int MIN_ERRS_TO_PRINT = -1;
    if ( VERBOSITY >=1 ) {
        cout << "Pre-corrected reads accuracy: " << endl;
        EvalCorrected( corrected, cid, ref, log_control, logc );
    }
    // mark which reads are corrected, and mapping to the list
    vec<bool> is_corrected( reads.size(), false );
    vec<int> mapping( reads.size(), -1 );
    for ( size_t i = 0; i < cid.size(); ++i ) {
        is_corrected[ cid[i] ] = true;
        mapping[ cid[i] ] = i;
    }
    // clean up some spurious {,} in input efasta
    VecEFasta corrected_clean( corrected.size() );
    for ( size_t i = 0; i < corrected.size(); ++i ) 
        corrected_clean[i] = CleanEfasta( corrected[i] );
    // find the direction of the friends
    vec< vec<int> > friend_ids2( friend_ids1.size() ); // a copy of friend_ids1, will be sorted
    for ( size_t i = 0; i < friend_ids1.size(); ++i ) {
        int rid = rid1[i];
        if ( ! is_corrected[rid] ) continue;
        friend_ids2[i] = friend_ids1[i];
        Sort( friend_ids2[i] );
        vec<int> signed_friend_ids;  // read id with negative sign meaning aligned in rc direction
        signed_friend_ids.reserve( friend_ids2[i].size() );
        for ( size_t j = 0; j < F[rid].size(); ++j ) {
            int fid = F[rid][j].getId();
            bool rc = F[rid][j].isRC();
            if ( BinMember( friend_ids2[i], fid ) ) 
                signed_friend_ids.push_back( rc ? ~fid : fid );
        }
        swap( friend_ids2[i], signed_friend_ids );
    }
    // collect all the friending information
    vec< vec<int> > known_friends_all( reads.size() );
    for ( size_t i = 0; i < friend_ids2.size(); ++i ) {
        int rid = rid1[i];
        if ( ! is_corrected[rid] ) continue;
        for ( size_t j = 0; j < friend_ids2[i].size(); ++j ) {
            bool rc = ( friend_ids2[i][j] < 0 );
            int fid = ( rc ? ~friend_ids2[i][j] : friend_ids2[i][j] );
            known_friends_all[ fid  ].push_back( rc ? ~rid : rid );
            //known_friends_all[ friend_ids2[i][j] ].push_back( rid );
            // friends of friends
            for ( size_t k = j+1; k < friend_ids2[i].size(); ++k ) {
                bool rc2 = ( friend_ids2[i][k] < 0 );
                int fid2 = ( rc2 ? ~friend_ids2[i][k] : friend_ids2[i][k] );
                known_friends_all[ fid ].push_back( rc ^ rc2 ? ~fid2 : fid2 );
                known_friends_all[ fid2 ].push_back( rc^rc2 ? ~fid : fid );
                //known_friends_all[ friend_ids2[i][j] ].push_back( friend_ids2[i][k] );
                //known_friends_all[ friend_ids2[i][k] ].push_back( friend_ids2[i][j] );
            }
        }
    }
    for ( size_t i = 0; i < known_friends_all.size(); ++i ) 
        UniqueSort( known_friends_all[i] );
    // dispatch the works
    vec<int> todo_list;
    todo_list.reserve( rid.size() );
    for ( size_t ix = 0; ix < rid.size(); ++ix ) 
        if ( ! is_corrected[ rid[ix] ] ) todo_list.push_back( rid[ix] );
    if ( VERBOSITY >= 1 ) {
        cout << "Number of reads to be corrected " << todo_list.size() << endl;
    }
    vec<String> outs( todo_list.size() );     // log output
    int count = 0;
    cout << Date() << ": Correcting " << todo_list.size() 
                   << " reads assisted by " << corrected.size() 
                   << " previously corrected reads." << endl;
    int batch_size = Min( 500, Max( 1, todo_list.isize() / omp_get_max_threads() ) );
    int n_no_cr_friends = 0, n_no_friends = 0, n_no_assembly = 0, 
        n_many_assembly = 0, n_short_assembly = 0, n_other_failures = 0;
    #pragma omp parallel for schedule( dynamic, batch_size )
    for ( size_t ix = 0; ix < todo_list.size(); ++ix ) {
        int ec = todo_list[ix];
        if ( Csmr_Rid != -1 && ec != Csmr_Rid ) continue;
        ostringstream out;
        const vec<int>& known_friends =  known_friends_all[ec];
        pair<efasta, CSMR_MESSAGE> result = CorrectOneReadFast( ec, reads, F, known_friends, 
                corrected_clean, cid, is_corrected, mapping, heur, log_control, VERBOSITY, out );
        //efasta result = CorrectOneRead( ec, reads, F, known_friends, corrected, cid, is_corrected, mapping,
        //        heur, log_control, VERBOSITY, out );
        #pragma omp critical
        {   
            if ( result.second == CSMR_OK ) {
                corrected2.push_back( result.first );
                cid2.push_back( ec );
            } 
            else if ( result.second == CSMR_NO_FRIENDS ) 
                n_no_friends++;
            else if ( result.second == CSMR_NO_CR_FRIENDS ) 
                n_no_cr_friends++;
            else if ( result.second == CSMR_NO_ASSEMBLY ) 
                n_no_assembly++;
            else if ( result.second == CSMR_MANY_ASSEMBLY ) 
                n_many_assembly++;
            else if ( result.second == CSMR_SHORT_ASSEMBLY )
                n_short_assembly++;
            else
                n_other_failures++;
            // progress dotting
            int dots1 = (100*count)/ todo_list.size();
            count++;
            int dots2 = (100*count)/ todo_list.size();
            if ( dots1 < dots2 ) {
                for( int j = dots1; j < dots2; j++) {
                    cout << ".";
                    if ( (j+1) % 50 == 0 ) cout << "\n";
                    else if ( (j+1) % 10 == 0 ) cout << " ";    
                }
                flush(cout);    
            }
        }  
        outs[ix] = out.str();
    }
    // Make sure the output is sorted by the read id
    // to prevent possible confusion of stochastic output results
    SortSync( cid2, corrected2 );
    if ( VERBOSITY >= 1 ) {
        for ( size_t i = 0; i < outs.size(); ++i ) 
            if ( ! outs[i].empty() )
                cout << "================== " << i << " id= " << todo_list[i] << " =============== \n"
                     << outs[i] << endl;
    }
    cout << Date() << ": Corrected " << cid2.size() << " reads." 
        << "\n " << n_no_friends << " attempts failed : not enough friend reads"
        << "\n " << n_no_cr_friends << " attempts failed : not enough corrected friend reads"
        << "\n " << n_no_assembly << " attempts failed : no linear assembly was formed"
        << "\n " << n_many_assembly << " attempts failed : too many linear assembles formed"
        << "\n " << n_short_assembly << " attempts failed : poor read-to-assembly mapping" 
        << "\n " << n_other_failures << " attempts failed : other reasons"
        << endl;
    if ( VERBOSITY >= 1 ) {
        cout << "corrected reads accuracy: " << endl;
        EvalCorrected( corrected2, cid2, ref, log_control, logc );
    }
    // What can be wrong with the initial correction?
    // - Errors in previously correctly reads
    // - Mis-assembly of the reads
    // - Theses reads are not true friends of the uncorrected read
    // - Misplacement of the uncorrected read to the assembly
}
