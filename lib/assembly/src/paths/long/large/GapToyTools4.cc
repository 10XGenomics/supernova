///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <atomic>

#include "CoreTools.h"
#include "ParallelVecUtilities.h"
#include "PrintAlignment.h"
#include "Qualvector.h"
#include "Set.h"
#include "kmers/BigKPather.h"
#include "kmers/LongReadPather.h"
#include "kmers/naif_kmer/KernelPerfectAligner.h"
#include "kmers/naif_kmer/NaifKmerizer.h"
#include "kmers/naif_kmer/Kmers.h"
#include "math/Functions.h"
#include "paths/HyperBasevector.h"
#include "paths/ReadsToPathsCoreX.h"
#include "paths/RemodelGapTools.h"
#include "paths/Unipath.h"
#include "paths/long/Correct1Pre.h"
#include "paths/long/ExtendReadPath.h"
#include "paths/long/FriendAligns.h"
#include "paths/long/HBVFromEdges.h"
#include "paths/long/KmerCount.h"
#include "paths/long/ReadPath.h"
#include "paths/long/ReadPathTools.h"
#include "paths/long/ReadStack.h"
#include "paths/long/large/Lines.h"
#include "paths/long/large/GapToyTools.h"
#include "random/Bernoulli.h"
#include "random/Random.h"
#include "system/WorklistN.h"

namespace {

// GenerateLookups - generate alignments between two sets of edges
//
// given a vecbasevector trans[] where the first hb.EdgeObjectCount() entries are
// either edges from hb or zero-size and the remaining entries are edges from
// hb3 (or zero-size).  Generate 200-mer based alignments (offset-only) between
// hb edges and hb3 edges.  These go into X.
//
// on exit, X is vec of triple( id1, offset, id2 )
//
void GenerateLookups( vecbasevector const& trans, HyperBasevector const& hb,
          vec<triple<int,int,int>>& X, Bool debug = False )
{
     const int K = 200;
     ForceAssertEq(K, hb.K());
     int n_id1 = hb.EdgeObjectCount();
     vec< triple<kmer<K>,int,int> > kmers_plus;
     MakeKmerLookup0( trans, kmers_plus );
     cout << "after kmers_plus creation, peak mem usage = " 
          << PeakMemUsageGBString( ) << endl;

     for ( int64_t i = 0; i < kmers_plus.jsize( ); i++ )
     {    int64_t j;
          for ( j = i + 1; j < kmers_plus.jsize( ); j++ )
               if ( kmers_plus[j].first != kmers_plus[i].first ) break;
          for ( int64_t k1 = i; k1 < j; k1++ )
          for ( int64_t k2 = i; k2 < j; k2++ )
          {    int id1 = kmers_plus[k1].second;
               int id2 = kmers_plus[k2].second - n_id1;
               // do we map an hb edge to an hb3 edge?
               if ( id1 >= n_id1 || id2 < 0 ) continue;
               int offset = kmers_plus[k1].third - kmers_plus[k2].third;
               X.push( id1, offset, id2 );
          }
          i = j - 1;    }
     ParallelUniqueSort(X);
     if ( debug ) {
          for ( auto const& x : X ) {
               int id1 = x.first; int offset = x.second; int id2 = x.third;
               PRINT3(id1,offset,id2);
          }
     }
}


void GenerateLookupsNaif( vecbasevector const& trans, HyperBasevector const& hb,
          vec<triple<int,int,int>>& X, Bool debug = False )
{
     const int K = 200;
     typedef Kmer248 Kmer_t;    // must be >= K
     ForceAssertEq(K, hb.K());
     int n_id1 = hb.EdgeObjectCount();

     size_t NUM_THREADS = configNumThreads(0);

     KmerAligns p;
     KernelAllKmerAligns<Kmer_t> kpa( trans, n_id1, K, &p,
               true /* ignore palindromes */, true /* ignore relative RC */ );
     naif_kmerize(&kpa, NUM_THREADS, debug);

     // push and pop to avoid double-memory or make naif_kmer use the same
     // type
     // -- could just make Unique() also more flexible
     X.clear();
     p.ReverseMe();
     while ( p.size() ) {
          if ( X.size() == 0 ||
                    X.back().first != p.back().first ||
                         X.back().second != p.back().second ||
                              X.back().third  != p.back().third ) {
               X.push_back(p.back());
          }
          p.pop_back();
     }


#if 0
     if ( debug ) cout << "BEFORE X.size()=" << p.size() << endl;
     Unique(p);
     if ( debug ) cout << "AFTER X.size()=" << p.size() << endl;

     X.clear_and_resize(p.size());
     for ( size_t i = 0; i < p.size(); ++i ) X[i] = p[i];
#endif
}

}; // end of anonymous namespace

void BuildAll( vecbasevector& allx, const HyperBasevector& hb, const int64_t extra )
{    size_t allxSize = hb.E( ) + extra;
     for ( auto itr=hb.To().begin(),end=hb.To().end(),itr2=hb.From().begin();
          itr != end; ++itr,++itr2 )
          allxSize += itr->size()*itr2->size();
     allx.reserve(allxSize);
     allx.resize( hb.EdgeObjectCount() );

     // populate allx with used edges from hb

     for ( int e = 0; e < hb.EdgeObjectCount( ); e++ )
          allx[e]=hb.EdgeObject(e);

     // populate allx with pseudo-edges defining a K+1 overlap
     // for every possible vertex crossing in hb

     const int K = hb.K( );
     bvec tmp(K+1);
     for ( int v = 0; v < hb.N( ); v++ )
     {    for ( int i1 = 0; i1 < hb.To(v).isize( ); i1++ )
          {    int e1 = hb.EdgeObjectIndexByIndexTo( v, i1 );
               for ( int i2 = 0; i2 < hb.From(v).isize( ); i2++ )
               {    int e2 = hb.EdgeObjectIndexByIndexFrom( v, i2 );
                    bvec const& x1 = hb.EdgeObject(e1);
                    bvec const& x2 = hb.EdgeObject(e2);
                    if ( x1.empty( ) || x2.empty( ) ) continue;
                    tmp.assign(x1.end()-K,x1.end()).push_back(x2[K-1]);
                    allx.push_back(tmp);    }    }    }    }

// Note that this truncates to length 1.

void TranslatePaths( ReadPathVec& paths2, const HyperBasevector& hb3,
     const vec<vec<int>>& to3, const vec<int>& left3 )
{
     #pragma omp parallel for
     for ( int64_t i = 0; i < (int64_t) paths2.size( ); i++ )
     {    if ( paths2[i].size( ) == 0 ) continue;
          int start = paths2[i].getOffset( ) + left3[ paths2[i][0] ];

          if ( to3[paths2[i][0]].empty( ) )
          {    paths2[i].clear( );
               continue;    }
     
          if ( start < hb3.Bases( to3[ paths2[i][0] ][0] ) )
          {    paths2[i].resize(1);
               paths2[i][0] = to3[ paths2[i][0] ][0];
               paths2[i].setOffset(start);
               continue;    }

          SerfVec<int> p;
          for ( int j = 0; j < (int) paths2[i].size( ); j++ )
          {    if ( to3[ paths2[i][j] ].empty( ) ) break;
               OverlapAppend( p, to3[paths2[i][j]] );    }
          int trim = 0;
          while( start >= hb3.EdgeLengthBases( p[trim] ) )
          {    start -= hb3.EdgeLengthKmers( p[trim] );
               trim++;
               if ( trim == (int) p.size( ) ) break;    }

          if ( trim == (int) p.size( ) )
          {    paths2[i].clear( );
               continue;    }
          paths2[i].resize(1);
          paths2[i][0] = p[trim];
          paths2[i].setOffset(start);    }    }

void AddNewStuff( vecbvec& new_stuff, HyperBasevector& hb, vec<int>& inv2, 
     ReadPathVec& paths2, const vecbasevector& bases, const VecPQVec& quals, 
     const int MIN_GAIN, const vec<int>& TRACE_PATHS, 
     const String& work_dir, const int EXT_MODE )
{
     vec<int> trace_edges;
     if ( TRACE_PATHS.size() ) 
     {    cout << "TRACE_PATHS:" << endl;
          for ( auto const i : TRACE_PATHS )
          {    const ReadPath& p = paths2[i];
               cout << "[" << i  << "] (" << p.getOffset( ) << ") " << printSeq(p)
                    << endl;
               std::copy( p.begin(), p.end(), 
                    std::back_inserter(trace_edges) );    }    }

     double clock1 = WallClockTime( );
     Validate( hb, paths2 );
     const int K = 200;
     ForceAssertEq( K, hb.K( ) );
     vec<Bool> used;
     hb.Used(used);

     HyperBasevector hb3;
     vec<vec<int>> to3( hb.EdgeObjectCount( ) );
     vec<int> left3( hb.EdgeObjectCount( ) );
     {
          ReadPathVec allx_paths;
          vecbasevector allx;
          BuildAll( allx, hb, new_stuff.size( ) );
          // add in "new_stuff" to allx
          allx.append(new_stuff.begin(),new_stuff.end());

          cout << Date( ) << ": building hb2" << endl;
          cout << TimeSince(clock1) << " used in new stuff 1 test" << endl;
          cout << "memory in use now = " << ToStringAddCommas( MemUsageBytes( ) )
               << endl;
          double clock2 = WallClockTime( );
          const int coverage = 4;
          buildBigKHBVFromReads( K, allx, coverage, &hb3, &allx_paths);
          cout << Date( ) << ": back from buildBigKHBVFromReads" << endl;

          // build to3 and left3 from allx_paths

          for ( int i = 0; i < hb.EdgeObjectCount(); ++i ) 
          {    for ( auto const& p : allx_paths[i] )
                    to3[i].push_back( p );
               left3[i] = allx_paths[i].getFirstSkip();    }

          cout << TimeSince(clock2) << " used in new stuff 2 test" << endl;

          if ( trace_edges.size() ) 
          {    cout << "BigKHBV EDGE-PATHS:" << endl;
               for ( auto const& edge : trace_edges )
                    cout << edge << ": " << allx_paths[edge] << endl;    }    }

     cout << "peak mem usage = " << PeakMemUsageGBString( ) << endl;
     double clock3 = WallClockTime( );
     TranslatePaths( paths2, hb3, to3, left3 );
     hb = hb3;
     hb.Involution(inv2);

     // Extend paths.

     Validate( hb, paths2 );
     Bool extend_paths_verbose = False;
     double clock5 = WallClockTime( );
     if (extend_paths_verbose) cout << "\nINCOMPLETE PATHS\n";
     vec<int> to_right;
     hb.ToRight(to_right);
     #pragma omp parallel for
     for ( int64_t i = 0; i < (int64_t) paths2.size( ); i++ )
     {    if ( paths2[i].size( ) > 0 ) paths2[i].resize(1);
          ExtendPath( paths2[i], i, hb, to_right, bases[i], quals.begin()[i],
                  MIN_GAIN, extend_paths_verbose, EXT_MODE );    }
     Validate( hb, paths2 );
     cout << TimeSince(clock5) << " used in new stuff 5" << endl;    }

// ExtendPath.
// - Does not change a path that has a negative start.
// - Considers only paths for which 
//   the read extends beyond the last edge in the path.
// - Considers all extensions of the given path in the graph, however returns
//   making no change if there are more than 100.
// - Considers only those extensions that cover the read.
// - Does nothing if there are no such extensions.
// - Does nothing unless there is an extension whose quality score sum beats the
//   runner up by at least min_gain, but in that case returns this best  extension.

void ExtendPath( ReadPath& p, const int64_t i, const HyperBasevector& hb, 
     const vec<int>& to_right, const bvec& bases,
     const qvec& quals, const int min_gain, const Bool verbose,
     const int mode )
{    if ( p.size( ) == 0 ) return;
     int K = hb.K( );
     int start = p.getOffset( );
     if ( start < 0 ) return;
     int rstop = hb.EdgeLengthBases( p[0] ) - start;
     for ( int j = 1; j < (int) p.size( ); j++ )
          rstop += hb.EdgeLengthKmers( p[j] );
     int ext = bases.isize( ) - rstop;
     if ( ext <= 0 ) return;
     int e = p.back( );
     int v = to_right[e];
     if ( hb.From(v).empty( ) ) return;
     if (verbose)
     {    cout << "\n[" << i << "] (" << start << ") " << printSeq(p)
               << " -- extends " << ext << " bases" << endl;    }
     vec<vec<int>> exts;
     vec<int> exts_len;
     exts.push_back( vec<int>( ) );
     exts_len.push_back(0);
     Bool fail = False;
     const int max_exts = 100;
     for ( int j = 0; j < exts.isize( ); j++ )
     {    if ( j > max_exts )
          {    fail = True;
               if (verbose) cout << "too many extensions" << endl;
               break;    }
          if ( exts_len[j] >= ext ) continue;
          int y;
          if ( exts[j].size( ) > 0 ) y = to_right[ exts[j].back( ) ];
          else y = v;
          for ( int l = 0; l < hb.From(y).isize( ); l++ )
          {    vec<int> e = exts[j];
               int n = hb.EdgeObjectIndexByIndexFrom( y, l );
               e.push_back(n);
               exts.push_back(e);
               exts_len.push_back( exts_len[j] + hb.EdgeLengthKmers(n) );    }    }
     if (fail) return;
     vec<Bool> to_del( exts.size( ), False );
     for ( int j = 0; j < exts.isize( ); j++ )
          if ( exts_len[j] < ext ) to_del[j] = True;
     EraseIf( exts, to_del );
     EraseIf( exts_len, to_del );
     if (verbose) PRINT( exts.size( ) );
     if ( exts.empty( ) ) return;
     int n = bases.size( );
     basevector r( bases, n - ext, ext );
     if (verbose) cout << "<0> " << r.ToString( ) << endl;
     vec<int> qsum( exts.size( ), 0 );
     for ( int j = 0; j < exts.isize( ); j++ )
     {    const vec<int>& e = exts[j];
          basevector b;
          for ( int l = 0; l < e.isize( ); l++ )
          {    b = Cat( b, basevector( hb.EdgeObject( e[l] ), K-1, 
                    hb.EdgeLengthBases( e[l] ) - (K-1) ) );    }
          for ( int m = 0; m < r.isize( ); m++ )
               if ( r[m] != b[m] ) qsum[j] += quals[ n - ext + m ];
          if (verbose)
          {    cout << "<" << j+1 << "> " << b.ToString( ) 
                    << ", qsum = " << qsum[j] << endl;    }    }
     SortSync( qsum, exts );
     if ( mode == 1 )
     {    if ( exts.size( ) >= 2 && qsum[1] - qsum[0] < min_gain ) return;
          for ( int j = 0; j < exts[0].isize( ); j++ )
               p.push_back( exts[0][j] );    }
     else // mode = 2
     {    int m;
          for ( m = 1; m < exts.isize( ); m++ )
               if ( qsum[m] - qsum[0] >= min_gain ) break;
          for ( int j = 0; j < exts[0].isize( ); j++ )
          {    Bool mis = False;
               for ( int l = 1; l < m; l++ )
               {    if ( j >= exts[l].isize( ) || exts[l][j] != exts[0][j] )
                    {    mis = True;
                         break;    }    }
               if (mis) return;
               p.push_back( exts[0][j] );    }    }    }

void ExtendPath2( ReadPath& p, const int64_t i, const HyperBasevector& hb, 
     const vec<int>& to_left, const vec<int>& to_right, const bvec& bases,
     const qvec& quals, const int min_gain, const Bool verbose,
     const int mode )
{    if ( p.size( ) == 0 ) return;
     int K = hb.K( );
     int start = p.getOffset( );
     int adds = 0;
     if (verbose)
          cout << "\n[" << i << "] (" << start << ") " << printSeq(p) << endl;
     if ( start < 0 ) 
     {    if ( mode == 0 ) return;

          // Try extending backwards.

          if (verbose) cout << "-- left extends " << -start << " bases" << endl;
          int e = p.front( );
          int v = to_left[e];
          while( hb.To(v).nonempty( ) && start < 0 )
          {    if ( hb.To(v).solo( ) )
               {    int f = hb.ITo( v, 0 );
                    p.push_front(f);
                    v = to_left[f];
                    start += hb.Kmers(f);
                    p.setOffset( Min( start, 0 ) );
                    adds++;
                    continue;    }
               vec<int> qsum( hb.To(v).size( ), 0 );
               vec<int> ids( hb.To(v).size( ), vec<int>::IDENTITY );
               for ( int j = 0; j < hb.To(v).isize( ); j++ )
               {    int f = hb.ITo( v, j );
                    int nf = hb.EdgeLengthBases(f);
                    for ( int l = K; l <= nf; l++ )
                    {    int rpos = -start - (l-K+1);
                         if ( rpos < 0 ) break;
                         if ( bases[rpos] != hb.EdgeObject(f)[nf-l] )
                              qsum[j] += quals[rpos];    }    }
               SortSync( qsum, ids );
               if ( qsum[0] == 0 && qsum[1] >= min_gain )
               {    int f = hb.ITo( v, ids[0] );
                    p.push_front(f);
                    v = to_left[f];
                    start += hb.Kmers(f);
                    p.setOffset( Min( start, 0 ) );
                    adds++;    }
               else break;    }
          if ( start < 0 ) 
          {    if ( verbose && adds > 0 ) 
                    cout << "extended to " << printSeq(p) << endl;
               return;    }    }

     int rstop = hb.EdgeLengthBases( p[0] ) - start;
     for ( int j = 1; j < (int) p.size( ); j++ )
          rstop += hb.EdgeLengthKmers( p[j] );
     int ext = bases.isize( ) - rstop;
     if ( ext <= 0 ) return;
     int e = p.back( );
     int v = to_right[e];
     if ( hb.From(v).empty( ) ) return;
     if (verbose) cout << "-- right extends " << ext << " bases" << endl;
     while( hb.From(v).nonempty( ) )
     {    if ( hb.From(v).solo( ) )
          {    int e = hb.IFrom( v, 0 );
               p.push_back(e);
               v = to_right[e];
               ext -= hb.Kmers(e);
               adds++;
               if ( ext < 0 ) break;
               continue;    }
          vec<int> qsum( hb.From(v).size( ), 0 );
          vec<int> ids( hb.From(v).size( ), vec<int>::IDENTITY );
          for ( int j = 0; j < hb.From(v).isize( ); j++ )
          {    int f = hb.IFrom( v, j );
               for ( int l = K-1; l < hb.EdgeLengthBases(f); l++ )
               {    int rpos = rstop + l - (K-1);
                    if ( rpos >= bases.isize( ) ) break;
                    if ( bases[rpos] != hb.EdgeObject(f)[l] )
                         qsum[j] += quals[rpos];    }    }
          SortSync( qsum, ids );
          if ( qsum[0] == 0 && qsum[1] >= min_gain )
          {    int e = hb.IFrom( v, ids[0] );
               p.push_back(e);
               v = to_right[e];
               ext -= hb.Kmers(e);
               adds++;
               if ( ext < 0 ) break;    }
          else break;    }
     if ( verbose && adds > 0 ) cout << "extended to " << printSeq(p) << endl;    }

void bubble_logger::bubble_data_t::addSupport(size_t branch, bubble_logger::bubble_data_t::support_t const&weight){
    SpinLocker locker(*lock_ptr);
    branch_supports[branch].push_back(weight);
}

// if one or more edges in rp is part of a bubble, perform gap-free alignment on the path as the alternate path
// and collect the result
// returns true if at any point (qsum of alt path) < (qsum of orig path)
bool bubble_logger::log_read(basevector const&read, qualvector const&qual, ReadPath const&rp, bool bVerbose){
    bool bErr = false;
    for(size_t rr=0;rr<rp.size();++rr){
        const int edge = rp[rr];
        const int other_edge = this->alt(edge);
        if( other_edge >= 0 ){
            ReadPath other_rp=rp;
            other_rp[rr]=other_edge;
            if( rr == 0){
                other_rp.setOffset(   hb_.EdgeObject(other_edge).isize()
                                    - hb_.EdgeObject(edge).size()
                                    + rp.getOffset()
                                  );
            }
            const int q_cur = getQ(read,qual,rp);
            const int q_alt = getQ(read,qual,other_rp);
            if( q_cur > q_alt){
                bErr = true;
                // if(bVerbose) std::cout << "WARNING: read-path of this read is not the lowest error path, q_alt < q_cur " << q_alt << " " << q_cur << std::endl;
                this->addWeight(other_edge,bubble_data_t::support_t(q_alt,q_cur-q_alt));
            }
            else{
                this->addWeight(edge,bubble_data_t::support_t(q_cur,q_alt-q_cur));
            }
        }
    }
    return bErr;
}
bubble_logger::bubble_logger(const HyperBasevector& hb, const vec<int>& inv)
                            :hb_(hb)
                            ,edge_alt_(hb.EdgeObjectCount(),-1)
                            ,edge_bubble_branch_(hb.EdgeObjectCount(),std::make_pair(-1,-1))
                            ,bubble_data_()
{
    // impromptu validation of the involution
    cout << Date() << ": checking involution" << endl;
    for ( int i = 0; i < inv.isize(); ++i )
        if ( i==-1 || i != inv[inv[i]] ) cout << "involution of " << i << " is wrong " << endl;
    cout << Date() << ": done" << endl;

    for(int vv=0;vv<hb.N();++vv){
        // only log bubble with edge topology e1->{b1,b2}->e2, e1 and e2 can be the same
        if(    hb.ToSize(vv) == 1
            && hb.FromSize(vv) == 2
            && hb.From(vv).front() != vv
            && hb.From(vv).front() == hb.From(vv).back()
            && hb.FromSize( hb.From(vv).front() ) == 1
            && hb.From( hb.From(vv).front() ).front() != hb.From(vv).front()
          ){
            const int edge_0 = hb.EdgeObjectIndexByIndexFrom(vv,0);
            const int edge_0_rc = inv[edge_0];
            const int edge_1 = hb.EdgeObjectIndexByIndexFrom(vv,1);
            const int edge_1_rc = inv[edge_1];

            if(   (edge_0_rc <0  && edge_1_rc >=0)
               || (edge_0_rc >=0 && edge_1_rc <0 )){
                // std::cout << "WARNING: inv-data is inconsistent rc("<<edge_0<<")=" << edge_0_rc << " rc("<<edge_1<<")="<<edge_1_rc << std::endl;
                continue;
            }

            if(edge_alt_[edge_0] < 0 || edge_alt_[edge_1]<0){
                if( edge_alt_[edge_0]>=0 || edge_alt_[edge_1]>=0){
                    std::cout << "inconsistent bubble, probably due to inversion properties" << std::cout;
                    continue;
                }

                const int bubble_idx = bubble_data_.size();

                edge_alt_[edge_0]=edge_1;
                edge_bubble_branch_[edge_0] = std::make_pair(bubble_idx,0);
                edge_alt_[edge_1]=edge_0;
                edge_bubble_branch_[edge_1] = std::make_pair(bubble_idx,1);

                Bool rc_available = False;

                if( edge_0_rc >=0 || edge_1_rc >=0) {

                    // bad inv?
                    if( edge_0_rc<0 || edge_1_rc<0) {
                        std::cout << "inconsistent bubble, probably due to inversion properties" << std::cout;

                    } else // no, inv is fine.

                        // naturally exclude self-inverse bubbles...
                        if ( edge_alt_[edge_0_rc] < 0 || edge_alt_[edge_1_rc]<0 ) {

                            if( edge_alt_[edge_0_rc]>=0 || edge_alt_[edge_1_rc]>=0) {
                                std::cout << "inconsistent bubble, probably due to inversion properties" << std::cout;
                            } else {
                                edge_alt_[edge_0_rc]=edge_1_rc;
                                edge_bubble_branch_[edge_0_rc] = std::make_pair(bubble_idx,2);
                                edge_alt_[edge_1_rc]=edge_0_rc;
                                edge_bubble_branch_[edge_1_rc] = std::make_pair(bubble_idx,3);
                                rc_available = True;
                            }
                        }
                }

                if ( rc_available )
                    bubble_data_.push_back(std::move(bubble_data_t(edge_0,edge_1,edge_0_rc,edge_1_rc)));
                else
                    bubble_data_.push_back(std::move(bubble_data_t(edge_0,edge_1)));
            }
            else{
                if( edge_alt_[edge_0]!=edge_1 && edge_alt_[edge_1]!=edge_0){
                    std::cout << "inconsistent bubble, probably due to inversion properties" << std::cout;
                    continue;
                }
            }
        }
    }

    // consistency check
    for ( size_t i = 0; i < edge_bubble_branch_.size(); ++i ) {
        auto const& ebb = edge_bubble_branch_[i];
        if ( ebb.first > -1 && ebb.first >= (int)bubble_data_.size() )
            FatalErr( ToString(i) + ": ebb.first=" +ToString(ebb.first)
                    + ", ebb.second=" + ToString(ebb.second) );
    }
}

//do a gap-free alignment of the read against the graph according to rp
//returns the sum of read-quality-score at the position the sequence mismatches
int bubble_logger::getQ(basevector const&read, qualvector const&qual, ReadPath const&rp, const qualvector::value_type min_q){
    int out=0;
    int bp=0;
    int shift=rp.getOffset();
    if( shift<0 ){
        bp = -shift;
        shift=0;
    }
    for(size_t ee=0;bp<read.isize() && ee<rp.size();++ee){
        basevector const& edge = hb_.EdgeObject(rp[ee]);
        for(size_t ep=shift ; bp < read.isize() && ep < edge.size() ; ++bp , ++ep){
            if( read[bp] != edge[ep] && qual[bp]>=min_q){
                out += qual[bp];
            }
        }
        shift=hb_.K()-1;
    }
    return out;
}
std::ostream& operator<<(std::ostream&os, bubble_logger const& in){
    for(const auto& entry: in.getData()){ os << "\n" << entry << "\n"; }
    return os;
};
std::ostream& operator<<(std::ostream& os, bubble_logger::bubble_data_t const&in){
    const auto& branch_edges = in.getEdges();

    typedef bubble_logger::bubble_data_t::support_t s_t;

    auto sum_branch_branch=[](vec<s_t>const&i){ return std::accumulate(i.begin(),i.end(),0 ,[](int i, s_t const&s){return i+s.branch_branch;}); };

    auto print_branch_branch=[](std::ostream& o, vec<s_t> i){
        if(i.size()==0) return;
        std::sort(i.begin(),i.end(),[](s_t const&L,s_t const&R){return L.branch_branch > R.branch_branch;} );
        o << i.front().branch_branch;
        for(size_t ii=1;ii<i.size();++ii){ o << "," << i[ii].branch_branch; }
    };
    const bool bHasFour=branch_edges.size()==4;

    os << branch_edges[0] << "," << branch_edges[1]<<":";

    os << "{ QSUM=" << sum_branch_branch(in.getSupport(0)) <<":";
    print_branch_branch(os,in.getSupport(0));
    if(bHasFour){
        os << " ; QSUM=" << sum_branch_branch(in.getSupport(2)) <<":";
        print_branch_branch(os,in.getSupport(2));
    }
    os << " },{ ";

    os << "QSUM=" << sum_branch_branch(in.getSupport(1)) <<":";
    print_branch_branch(os,in.getSupport(1));
    if(bHasFour){
        os << " ; QSUM=" << sum_branch_branch(in.getSupport(3)) <<":";
        print_branch_branch(os,in.getSupport(3));
    }
    os << " }";

    return os;
};

namespace{

typedef size_t work_item_t;

class LogBubblesProcessor
{
public:
    LogBubblesProcessor( bubble_logger& logger
                       , vecbasevector const&bases, vecqualvector const&quals, ReadPathVec const& paths2)
                       : logger_(logger)
                       , bases_(bases)
                       , quals_(quals)
                       , paths2_(paths2){}
    void operator()(work_item_t const&item){
         if( logger_.log_read(bases_[item], quals_[item], paths2_[item]) ){
             SpinLocker ml(message_lock);
             if( nWarnings < nMaxVerbose ){
                 // std::cout << "WARNING: read-path of read " << item << " is not the lowest error path, q_alt < q_cur " << std::endl;
                 if(nWarnings+1==nMaxVerbose){
                     // std::cout << "WARNING: read-path warning silenced" << std::endl;
                 }
             }
             ++nWarnings;
         }
    }
    static std::atomic<size_t> nWarnings;
private:
    bubble_logger& logger_;
    vecbasevector const& bases_;
    vecqualvector const& quals_;
    ReadPathVec const&paths2_;
    static const size_t nMaxVerbose=10;
    static SpinLockedData message_lock;

};
std::atomic<size_t> LogBubblesProcessor::nWarnings(0);
SpinLockedData LogBubblesProcessor::message_lock;

void LogBubbles( bubble_logger& logger, HyperBasevector& hb , const vec<int>& inv2
               , const vecbasevector & bases, const VecPQVec& quals, const ReadPathVec& paths2){
     ForceAssertEq( bases.size(), quals.size() );
     ForceAssertEq( bases.size(), paths2.size() );
     ForceAssert(hb.EdgeObjectCount()==inv2.isize());


     //parallelized version -- not turn on at this moment to order the read-path warnings
/*
     parallelFor(0ul,bases.size(),LogBubblesProcessor(logger,bases,quals,paths2));
     if(LogBubblesProcessor::nWarnings>0){
         std::cout << "WARNING: " << LogBubblesProcessor::nWarnings << " suspicious read-paths." << std::endl;
     }
*/

     const size_t nMaxVerbose=10;
     size_t nWarnings=0;
     auto qvItr = quals.begin();
     for(size_t rr=0;rr<bases.size();++rr,++qvItr){
         if( logger.log_read(bases[rr], *qvItr, paths2[rr]) ){
             if( nWarnings < nMaxVerbose ){
                 // std::cout << "WARNING: read-path of read " << rr << " is not the lowest error path, q_alt < q_cur " << std::endl;
                 if( nWarnings+1 == nMaxVerbose){
                     // std::cout << "WARNING: read-path warning silenced" << std::endl;
                 }
             }
             ++nWarnings;
         }
     }
     if(nWarnings>0){
         std::cout << "WARNING: " << nWarnings << " suspicious read-paths." << std::endl;
     }
}


}

void PrintBubbles( std::ostream& os, HyperBasevector& hb , const vec<int>& inv2
               , const vecbasevector & bases, const VecPQVec& quals, const ReadPathVec& paths2){
     bubble_logger logger( hb , inv2 );
     LogBubbles(logger,hb,inv2,bases,quals,paths2);
     os << logger;
}

void PopBubbles( HyperBasevector& hb , const vec<int>& inv2
               , const vecbasevector & bases, const VecPQVec& quals, const ReadPathVec& paths2){

     // set up bubble logger, then go through each read
     bubble_logger logger( hb , inv2 );
     LogBubbles(logger,hb,inv2,bases,quals,paths2);

     auto sum_expected_number_of_reads=[]( vec<bubble_logger::bubble_data_t::support_t>const& branch0
                                         , vec<bubble_logger::bubble_data_t::support_t>const& branch1
                                         ){
         std::pair<double,double> out(0.0,0.0);

         for(const auto&entry:branch0){
             const double p = std::max(0.5 , 1.0 - pow(10.0, - 0.1 * entry.branch_branch));
             out.first += p;
             out.second += 1.0 - p;
         }
         for(const auto&entry:branch1){
             const double p = std::max(0.5 , 1.0 - pow(10.0, - 0.1 * entry.branch_branch));
             out.first += 1.0 - p;
             out.second += p;
         }
         return out;
     };

     // delete edges based on Divine Bubble's method , with the exception that fw/rc components are treated separately using two p-value tests
     vec<int> to_delete;
     for(const auto&data: logger.getData()){
          vec<int> edges = data.getEdges();

          // don't try to pop a bubble where one edge is the RC of the other.
          // the involution should be symmetric, but the world is an
          // imperfect place
          if ( edges.size() == 2 &&
                  ( inv2[edges[0]] == edges[1] || inv2[edges[1]] == edges[0] ) )
              continue;

          auto tmp = sum_expected_number_of_reads(data.getSupport(0),data.getSupport(1));
          double f1 = tmp.first;
          double f2 = tmp.second;
          double r1,r2;
          if(edges.size()==2){ // if there's no seperate RC edges, approximate by splitting the weight
              f1 *= 0.5;
              f2 *= 0.5;
              r1 = f1;
              r2 = f2;
          }
          else if(edges.size()==4){
              auto tmp = sum_expected_number_of_reads(data.getSupport(2),data.getSupport(3));
              r1 = tmp.first;
              r2 = tmp.second;
          }
          else{
              std::cout << "WARNING: bubble logger was not used correctly, PopBubbles did not change the graph." << std::endl;
              return;
          }

          //this is copy-and-pasted from longproto's Divine Bubble
          int shift = 1;
          if ( f2 + r2 > f1 + r1 || ( f2 + r2 == f1 + r1 && f2 > f1 ) )
          {    shift = 0;
               swap( f1, f2 );
               swap( r1, r2 );    }

          const int n_f = int(floor( 2*(f1+f2) ));
          const int n_r = int(floor( 2*(r1+r2) ));

          //heuristics just like those in LongProto's DivineBubbles
          const double max_asym_rarity = 0.00001;
          const double min_to_save = 10;

          if (    ( n_f>0 || n_r > 0)
               && ( n_f==0 || ( BinomialSum( n_f, int(ceil(f2)), 0.25 ) < max_asym_rarity && f2 < min_to_save ) )
               && ( n_r==0 || ( BinomialSum( n_r, int(ceil(r2)), 0.25 ) < max_asym_rarity && r2 < min_to_save ) )
             ) {
                to_delete.push_back( edges[0+shift] );
                if(edges.size()==4) to_delete.push_back(edges[2+shift]);
           }
     }
     hb.DeleteEdges(to_delete);
}

void DeleteFunkyPathPairs( const HyperBasevector& hb, const vec<int>& inv,
     const vecbasevector& bases, ReadPathVec& paths, const Bool verbose )
{
     double clock = WallClockTime( );
     vec<int> to_left, to_right;
     hb.ToLeft(to_left), hb.ToRight(to_right);

     // Define heuristics.

     const int max_homes = 20;
     const int max_frag = 1000;
     const int K = 40;
     const int min_gain = 5;

     // Compute distances to ends, approximately.  This is a bit dangerous.

     // cout << Date( ) << ": estimating distances to end" << endl;
     vec<int> D( hb.N( ), 1000000000 );
     const int max_check = 20;
     vec<int> xc;
     for ( int v = 0; v < hb.N( ); v++ )
     {    if ( hb.From(v).empty( ) ) 
          {    D[v] = 0;
               xc.push_back(v);    }    }
     for ( int p = 0; p < max_check; p++ )
     {    vec<int> xc2;
          for ( int i = 0; i < xc.isize( ); i++ )
          {    int v = xc[i];
               for ( int j = 0; j < hb.To(v).isize( ); j++ )
               {    int w = hb.To( xc[i] )[j];
                    int e = hb.EdgeObjectIndexByIndexTo( v, j );
                    D[w] = Min( D[w], D[v] + hb.EdgeLengthKmers(e) );
                    xc2.push_back(w);    }    }
          xc = xc2;    }

     // Find badly placed pairs.

     // cout << Date( ) << ": finding badly placed pairs" << endl;
     int64_t npids = paths.size( ) / 2;
     vec<Bool> invalid(npids, False);
     #pragma omp parallel for
     for ( int64_t pid = 0; pid < npids; pid++ )
     {    int64_t id1 = 2*pid, id2 = 2*pid + 1;
          ReadPath &p1 = paths[id1], &p2 = paths[id2];
          if ( p1.size( ) == 0 || p2.size( ) == 0 ) continue;
          vec<int> x1, x2;
          for ( int j = 0; j < (int) p1.size( ); j++ )
               x1.push_back( p1[j] );
          for ( int j = ( (int) p2.size( ) ) - 1; j >= 0; j-- )
               x2.push_back( inv[ p2[j] ] );
          Bool equal = ( p1 == p2 );
          if (equal) 
          {    invalid[pid] = True;
               continue;    } // THESE SHOULD BE DELETED.
          int start1 = p1.getOffset( );
          int start2 = hb.EdgeLengthBases( p2.front( ) ) - p2.getOffset( );

          const int min_frag = 50;
          const int max_frag = 1300;

          if ( x1.solo( ) && x1 == x2 
               && start2 - start1 >= min_frag && start2 - start1 <= max_frag )
          {    continue;    }

          int dist_to_end1 = hb.EdgeLengthBases( p1.front( ) ) 
               - p1.getOffset( ) - bases[id1].size( ) + D[ to_right[ p1.back( ) ] ];
          for ( int j = 1; j < (int) p1.size( ); j++ )
               dist_to_end1 -= hb.EdgeLengthKmers( p1[j] );
          int dist_to_end2 = hb.EdgeLengthBases( p2.front( ) ) 
               - p2.getOffset( ) - bases[id2].size( ) + D[ to_right[ p2.back( ) ] ];
          for ( int j = 1; j < (int) p2.size( ); j++ )
               dist_to_end2 -= hb.EdgeLengthKmers( p2[j] );
          if ( dist_to_end1 + dist_to_end2 
                    + bases[id1].isize( ) + bases[id2].isize( ) <= max_frag )
          {    continue;    }

          const int max_exts = 10;
          vec<vec<int>> paths = {x1};
          Bool good = False;
          for ( int e = 0; e <= max_exts; e++ )
          {    for ( int j = 0; j < paths.isize( ); j++ )
               {    const vec<int>& p = paths[j];
                    if ( x2.size( ) <= p.size( ) 
                         && p.Contains( x2, p.size( ) - x2.size( ) ) )
                    {    int s1 = start1;
                         for ( int l = 0; l < p.isize( ) - 1; l++ )
                              s1 -= hb.EdgeLengthKmers( p[l] );
                         if ( start2 - s1 >= min_frag && start2 - s1 <= max_frag )
                         {    good = True;
                              goto done;    }    }    }
               vec<vec<int>> paths2;
               for ( int j = 0; j < paths.isize( ); j++ )
               {    int v = to_right[ paths[j].back( ) ];
                    for ( int l = 0; l < hb.From(v).isize( ); l++ )
                    {    vec<int> x = paths[j];
                         x.push_back( hb.EdgeObjectIndexByIndexFrom( v, l ) );
                         int n = 0;
                         for ( int m = 1; m < x.isize( ) - 1; m++ )
                              n += hb.EdgeLengthKmers( x[m] );
                         if ( n < max_frag ) paths2.push_back(x);    }    }
               paths = paths2;
               if ( paths.empty( ) ) break;
               if ( e == max_exts ) good = True;    } // not good but giving up
          done: if (good) continue;
          if (verbose)
          {    cout << "[" << pid << "] " << start1 << ":" << printSeq(x1) << ".." 
                    << start2 << ":" << printSeq(x2) << endl;    }
          invalid[pid] = True;     }
     if (verbose) cout << "\n";
     PRINT2( Sum(invalid), npids );
     for ( int64_t pid = 0; pid < npids; pid++ )
     {    if ( invalid[pid] )
          {    paths[2*pid].resize(0);
               paths[2*pid+1].resize(0);    }    }
     LogTime( clock, "finding funky pairs" );    }

void ReduceK( const int newK, HyperBasevector& hb, const vec<int>& inv,
     ReadPathVec& paths )
{    int oldK = hb.K( );
     int delta = oldK - newK;
     hb.ReduceK(newK);
     vec<int> to_left, to_right;
     hb.ToLeft(to_left), hb.ToRight(to_right);
     for ( int64_t pi = 0; pi < (int64_t) paths.size( ); pi++ )
     {    ReadPath& p = paths[pi];
          if ( p.size( ) == 0 ) continue;
          int v = to_left[ p[0] ];
          int in = hb.To(v).size( ), out = hb.From(v).size( );
          int trim;
          if ( in == 0 ) trim = 0;
          else if ( in == 1 && out > 1 ) trim = delta;
          else if ( in > 1 && out == 1 ) trim = 0;
          else trim = delta/2;
          p.setOffset( p.getOffset( ) - trim );
          if ( p.getOffset( ) < 0 && hb.To(v).solo( ) )
          {    int e = hb.EdgeObjectIndexByIndexTo( v, 0 );
               p.push_front(e);
               p.setOffset( p.getOffset( ) + hb.EdgeLengthKmers(e) );    }    }    }

void AssayMisassemblies( const HyperBasevector& hbx, const vec<int>& inv,
     const vec< vec< pair<int,int> > >& hits, const vec<vec<vec<vec<int>>>>& linesx,
     const String& final_dir )
{    Ofstream( out, final_dir + "/misassemblies" );
     vec<int> lens;
     GetLineLengths( hbx, linesx, lens );
     const int min_len = 30000;
     int64_t assayed = 0;
     int bads = 0;
     set< pair<int,int> > seen;
     for ( int i = 0; i < linesx.isize( ); i++ )
     {    const vec<vec<vec<int>>>& L = linesx[i];
          if ( lens[i] < min_len ) continue;
          int e1 = -1, e2 = -1, g1 = -1, g2 = -1;
          for ( int j = 0; j < L.isize( ); j += 2 )
          {    int e = L[j][0][0];
               if ( hits[e].size( ) + hits[ inv[e] ].size( ) == 1 )
               {    if ( hits[e].solo( ) ) g1 = hits[e][0].first;
                    else g1 = hits[ inv[e] ][0].first;
                    if ( g1 < 24 )
                    {    e1 = e;
                         break;    }    }    }
          for ( int j = L.isize( ) - 1; j >= 0; j -= 2 )
          {    int e = L[j][0][0];
               if ( hits[e].size( ) + hits[ inv[e] ].size( ) == 1 )
               {    if ( hits[e].solo( ) ) g2 = hits[e][0].first;
                    else g2 = hits[ inv[e] ][0].first;
                    if ( g2 < 24 )
                    {    e2 = e;
                         break;    }    }    }
          if ( e1 < 0 || e2 < 0 ) continue;
          if ( Member( seen, make_pair( inv[e2], inv[e1] ) ) ) continue;
          seen.insert( make_pair( e1, e2 ) );
          assayed += lens[i];
          if ( g1 != g2 ) 
          {    int b1 = -1, b2 = -1;
               for ( int j = L.isize( ) - 1; j >= 0; j -= 2 )
               {    int e = L[j][0][0];
                    for ( int l = 0; l < hits[e].isize( ); l++ )
                    {    if ( hits[e][l].first == g1 )
                         {    b1 = e;
                              break;    }    }
                    for ( int l = 0; l < hits[ inv[e] ].isize( ); l++ )
                    {    if ( hits[ inv[e] ][l].first == g1 )
                         {    b1 = e;
                              break;    }    }
                    if ( b1 >= 0 ) break;    }
               for ( int j = 0; j < L.isize( ); j += 2 )
               {    int e = L[j][0][0];
                    for ( int l = 0; l < hits[e].isize( ); l++ )
                    {    if ( hits[e][l].first == g2 ) 
                         {    b2 = e;
                              break;    }    }
                    for ( int l = 0; l < hits[ inv[e] ].isize( ); l++ )
                    {    if ( hits[ inv[e] ][l].first == g2 ) 
                         {    b2 = e;
                              break;    }    }
                    if ( b2 >= 0 ) break;    }
               out << ++bads << ".  L=" << i << ", len: " << lens[i] 
                    << ", e: " << e1 << " / " << e2 << ", g: " << g1 << " / " << g2
                    << ", b: " << b1 << " / " << b2 << endl;    }    }
     if ( assayed > 0 )
     {    double bad_rate = ( 1000000000.0 * bads ) / assayed;
          cout << "see " << bads << " putative interchromosomal misassemblies, "
               << "rate = " << ToString(bad_rate,1) << " per Gb" << endl;
          PerfStatLogger::log("misassemblies",bads,
                          "Putative interchromosomal misassemblies.");    }    }

String Chr( const int g )
{    if ( g < 22 ) return ToString(g+1);
     if ( g == 22 ) return "X";
     if ( g == 23 ) return "Y";
     return "?";    }

void MakeFinalFasta( const HyperBasevector& hbx, const vec<int>& inv2,
     const vec<vec<vec<vec<int>>>>& linesx, const vec<int>& npairsx,
     const vec<vec<covcount>>& covs, vec< vec< pair<int,int> > > hits,
     const String& final_dir, const String& work_dir, const Bool ALIGN_TO_GENOME )
{    double clock1 = WallClockTime( );
     vec<int> llensx;
     GetLineLengths( hbx, linesx, llensx );
     Ofstream( out, final_dir + "/a.fasta" );
     vec<int> to_left, to_right;
     hbx.ToLeft(to_left), hbx.ToRight(to_right);
     vec<String> head( hbx.E( ) );
     for ( int e = 0; e < hbx.E( ); e++ ) 
     {    int v = to_left[e], w = to_right[e];
          head[e] = ToString(e) + " " + ToString(v) + ":" + ToString(w);    }

     // Add coverage values to headers.

     int ns = covs.size( );
     for ( int e = 0; e < hbx.E( ); e++ )
     {    Bool defined = False;
          int ns = covs.size( );
          for ( int ss = 0; ss < ns; ss++ )
               if ( covs[ss][e].Def( ) ) defined = True;
          if (defined)
          {    head[e] += " cov=";
               for ( int ss = 0; ss < ns; ss++ )
               {    if ( ss > 0 ) head[e] += ",";
                    if ( covs[ss][e].Def( ) )
                         head[e] += ToString( covs[ss][e].Cov( ), 2 ) + "x";
                    else head[e] += "?x";    }    }    }

     if ( ALIGN_TO_GENOME && IsRegularFile( work_dir + "/genome.fastb" ) )
     {    for ( int e = 0; e < hits.isize( ); e++ ) 
          {    int re = inv2[e];
               if ( hits[e].empty( ) && ( re < 0 || hits[re].empty( ) ) ) 
                    continue;
               UniqueSort(hits[e]), UniqueSort(hits[re]);
               head[e] += " loc=";
               int c = 0;
               for ( int j = 0; j < hits[e].isize( ); j++ )
               {    if ( c++ > 0 ) head[e] += ",";
                    head[e] += "+" + Chr(hits[e][j].first)
                         + ":" + ToString(hits[e][j].second);    }
               if ( re < 0 ) continue;
               for ( int j = 0; j < hits[re].isize( ); j++ )
               {    if ( c++ > 0 ) head[e] += ",";
                    head[e] += "-" + Chr(hits[re][j].first)
                         + ":" + ToString(hits[re][j].second);    }    }    }
     cout << TimeSince(clock1) << " using setting up final fasta" << endl;
     double clock2 = WallClockTime( );
     for ( int e = 0; e < hbx.E( ); e++ ) 
     {    out << ">" << head[e] << "\n";
          hbx.EdgeObject(e).Print(out);    }
     cout << TimeSince(clock2) << " using printing final fasta" << endl;    }

void BuildGenomeMap( const HyperBasevector& hb, const vec<int>& inv,
     const vec<vec<vec<vec<int>>>>& lines, const vec<vec<covcount>>& covs,
     const vec< vec< pair<int,int> > >& hits, const vecbitvector& genome_amb,
     const vec<String>& genome_names, const String& final_dir )
{
     // Get line lengths.

     vec<int> lens;
     GetLineLengths( hb, lines, lens );

     // Remove hits beyond chr 1..22,X,Y and nearby hits.

     vec< vec< pair<int,int> > > hitsx(hits);
     for ( int e = 0; e < hits.isize( ); e++ )
     {    vec<Bool> to_del( hits[e].size( ), False );
          for ( int j = 0; j < hits[e].isize( ); j++ )
          {    if ( hits[e][j].first > 23 ) to_del[j] = True;
               if ( j > 0 && hits[e][j].first == hits[e][j-1].first
                    && hits[e][j].second - hits[e][j-1].second < 1000 )
               {    to_del[j] = True;    }    }
          EraseIf( hitsx[e], to_del );    }

     // Locate gaps in reference.

     vec< quad<int,int,int,int> > locs;
     for ( int g = 0; g < 24; g++ )
     {    for ( int j = 0; j < (int) genome_amb[g].size( ); j++ )
          {    if ( !genome_amb[g][j] ) continue;
               int k;
               for ( k = j + 1; k < (int) genome_amb[g].size( ); k++ )
                    if ( !genome_amb[g][k] ) break;
               if ( k - j >= 10000 ) locs.push( g, j, k, -1 );
               j = k - 1;    }    }

     // Locate lines on genome.

     for ( int i = 0; i < lines.isize( ); i++ )
     {    const vec<vec<vec<int>>>& L = lines[i];
          if ( lens[i] < 1000 ) continue;
          pair<int,int> start( -1, -1 ), stop( -1, -1 );
          for ( int j = 0; j < L.isize( ); j += 2 )
          {    int e = L[j][0][0];
               int re = inv[e];
               if ( e >= hitsx.isize( ) ) continue;
               if ( hitsx[e].size( ) + hitsx[re].size( ) != 1 ) continue;
               if ( hitsx[e].size( ) != 1 ) break;
               start = hitsx[e][0];
               break;    }
          for ( int j = L.isize( ) - 1; j >= 0; j -= 2 )
          {    int e = L[j][0][0];
               int re = inv[e];
               if ( e >= hitsx.isize( ) ) continue;
               if ( hitsx[e].size( ) + hitsx[re].size( ) != 1 ) continue;
               if ( hitsx[e].size( ) != 1 ) break;
               stop = hitsx[e][0];
               stop.second += hb.EdgeLengthBases(e);
               break;    }
          if ( start.first < 0 || stop.first < 0 ) continue;
          if ( start.first != stop.first ) continue;
          if ( !( start.second <= stop.second ) ) continue;
          locs.push( start.first, start.second, stop.second, i );    }
     Sort(locs);
     vec<vec<String>> rows;
     int M = 0;
     for ( int i = 0; i < locs.isize( ); i++ )
     {    int chr = locs[i].first, start = locs[i].second, stop = locs[i].third;
          if ( i > 0 && locs[i].first != locs[i-1].first )
          {    vec<String> row;
               rows.push_back(row);
               M = 0;    }
          M = Max( M, stop );
          int l = locs[i].fourth;
          if ( l < 0 )
          {    vec<String> row;
               row.push_back( genome_names[chr], ToStringAddCommas(start), 
                    ToStringAddCommas(stop), "**REF GAP**" );
               rows.push_back(row);
               continue;    }
          int e1 = lines[l].front( )[0][0], e2 = lines[l].back( )[0][0];
          if ( i > 0 && locs[i].first == locs[i-1].first 
               && locs[i].second - M >= 10000 )
          {    vec<String> row;
               int gap = locs[i].second - M;
               String msg = "gap";
               if ( gap >= 50000 ) msg = "GAP";
               row.push_back( 
                    String("--"), "----" + msg + "----", String( "-----------") );
               rows.push_back(row);    }
          vec<String> row;
          row.push_back( genome_names[chr], ToStringAddCommas(start), 
               ToStringAddCommas(stop), "line[" + ToString(l) + "]", 
               ToString(e1) + ".." + ToString(e2), "len=" + ToString(lens[l]) );
          int ns = covs.size( );
          Bool defined = False;
          for ( int ss = 0; ss < ns; ss++ )
               if ( covs[ss][e1].Def( ) ) defined = True;
          if (defined)
          {    String x = "cov=";
               for ( int ss = 0; ss < ns; ss++ )
               {    if ( ss > 0 ) x += ";";
                    if ( covs[ss][e1].Def( ) )
                         x += ToString( covs[ss][e1].Cov( ), 2 );
                    else x += "?";    }
               row.push_back(x);    }
          rows.push_back(row);    }
     Ofstream( out, final_dir + "/a.lines.map" );
     PrintTabular( out, rows, 1, "lrrllll" );    }

String PrintHits( const int e, vec<vec<pair<int,int>>>& hits, 
     const HyperBasevectorX& hb, const vec<int>& inv, 
     const vec<String>& genome_names )
{    
     String answer;
     int re = inv[e];
     for ( int pass = 1; pass <= 2; pass++ )
     {    int f = ( pass == 1 ? e : re );
          Sort( hits[f] );
          for ( int j = 0; j < hits[f].isize( ); j++ )
          {    int k, g = hits[f][j].first, start = hits[f][j].second;
               for ( k = j + 1; k < hits[f].isize( ); k++ )
               {    if ( hits[f][k].first != g ) break;
                    if ( hits[f][k].second - hits[f][k-1].second > 100 ) break;    }
               answer += " (" + ToString( pass == 1 ? "+" : "-" );
               String chr;
               if ( genome_names.empty( ) ) chr = ToString(g+1);
               else chr = genome_names[g];
               int stop = hits[f][k-1].second + hb.Bases(e);
               String xstart = ToStringAddCommas(start);
               String xstop = ToStringAddCommas(stop);
               if ( xstart.size( ) == xstop.size( ) )
               {    for ( int m = 0; m < xstart.isize( ); m++ )
                    {    if ( xstart[m] != xstop[m] )
                         {    xstop = xstop.substr( m, xstop.isize( ) - m );
                              break;    }    }    }
               answer += chr + ":" + xstart + "-" + xstop + ")";
               j = k - 1;    }    }
     return answer;    }

int N50PerfectStretch( const HyperBasevector& hb, const vecbasevector& genome,
     const Bool concatenate )
{    ForceAssert( genome.size( ) > 0 );
     vec< triple< pair<int,int>, pair<int,int>, int > > perfs;
     vec<perf_place> places;
     AlignToGenomePerf( hb, genome, perfs, places );
     if ( places.empty( ) ) return 0;
     vec<vec<perf_place>> PLACES( genome.size( ) );
     for ( int i = 0; i < places.isize( ); i++ )
          PLACES[ places[i].G( ) ].push_back( places[i] );

     if ( !concatenate )
     {    vec<int> lens;
          for ( int i = 0; i < places.isize( ); i++ )
               lens.push_back( places[i].Len( ) );
          Sort(lens);
          return N50(lens);    }

     vec<int> lens, gs, N50s;
     for ( int g = 0; g < (int) genome.size( ); g++ )
          gs.push_back(g);
     int reps = ( genome.size( ) == 1 ? 1 : 100 );
     srandomx(5386927);
     for ( int r = 0; r < reps; r++ )
     {    vec<int> gs2;
          while( gs2.size( ) < gs.size( ) )
          {    int x = randomx( ) % gs.size( );
               if ( !Member( gs2, x ) ) gs2.push_back(x);    }
          vec<perf_place> placesx;
          for ( int m = 0; m < gs2.isize( ); m++ )
               placesx.append( PLACES[ gs2[m] ] );
          for ( int i = 0; i < placesx.isize( ); i++ )
          {    int l = placesx[i].Len( );
               for ( int j = i + 1; j < placesx.isize( ); j++ )
               {    if ( places[j].G( ) != placesx[j-1].G( )
                         && placesx[j-1].Gstop( ) 
                              == genome[ placesx[j-1].G( ) ].isize( )
                         && placesx[j].Gstart( ) == 0 )
                    {    l += placesx[j].Len( );
                         i++;    }
                    else break;    }
               lens.push_back(l);    }
          Sort(lens);
          N50s.push_back( N50(lens) );    }
     Sort(N50s);
     return int(round(Mean(N50s)));    }
