///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"
#include "ParallelVecUtilities.h"
#include "math/Functions.h"
#include "pairwise_aligners/SmithWatFree.h"
#include "paths/long/CreateGenome.h"
#include "paths/long/EvalByReads.h"
#include "paths/long/LargeKDispatcher.h"
#include "paths/long/SupportedHyperBasevector.h"
#include "random/Bernoulli.h"
#include "Map.h"

namespace
{

template<int K> void WingleSupport( const SupportedHyperBasevector& shb,
     const vecbasevector& bases_local, const vecqualvector& quals_local, 
     vec<int>& to_delete, const int verbosity )
{
     cout << Date( ) << ": building kmers_plus, mem in use = "
          << setiosflags(ios::fixed) << setprecision(1) 
          << MemUsage( ) / (1024.0*1024.0) << resetiosflags(ios::fixed) 
          << " GB" << endl;
     vec< triple<kmer<K>,int,int> > kmers_plus;
     {    vec<int64_t> starts;
          starts.push_back(0);
          const size_t nBases=bases_local.size();
          for ( size_t i = 0; i < nBases; i++ )
          {    const basevector& u = bases_local[i];
               starts.push_back( 
                    starts.back( ) + Max( 0, u.isize( ) - K + 1 ) );    }
          cout << Date( ) << ": mem in use = "
               << setiosflags(ios::fixed) << setprecision(1) 
               << MemUsage( ) / (1024.0*1024.0) << resetiosflags(ios::fixed) 
               << " GB, about to allocate "
               << setiosflags(ios::fixed) << setprecision(1) 
               << double( starts.back( ) * sizeof( triple<kmer<K>,int,int> ) )
               / (1024.0*1024.0*1024.0) 
               << resetiosflags(ios::fixed) << " GB" << endl;
          kmers_plus.resize( starts.back( ) );
          #pragma omp parallel for schedule(dynamic,1)
          for ( size_t i = 0; i < nBases; i++ )
          {    const basevector& u = bases_local[i];
               for ( int j = 0; j <= u.isize( ) - K; j++ )
                {    int64_t r = starts[i] + j;
                     kmers_plus[r].first.SetToSubOf( u, j );
                     kmers_plus[r].second = i;
                     kmers_plus[r].third = j;    }    }
          cout << Date( ) << ": sorting kmers_plus, mem in use = " 
               << setiosflags(ios::fixed) << setprecision(1) 
               << MemUsage( ) / (1024.0*1024.0)
               << resetiosflags(ios::fixed) << " GB" << endl;
           sortInPlaceParallel( kmers_plus.begin( ), kmers_plus.end( ) );    }
     cout << Date( ) << ": done creating kmers_plus, mem in use = " 
          << setiosflags(ios::fixed) << setprecision(1) << MemUsage( ) / (1024.0*1024.0)
          << resetiosflags(ios::fixed) << " GB" << endl;
      const int min_mult = 10;
      const int min_win = 100;
      const int min_win_ratio = 8;
      const int max_lose = 100;

      const int poly = 10;
      const basevector chain0(poly,0);
      const basevector chain1(poly,1);
      const basevector chain2(poly,2);
      const basevector chain3(poly,3);
      #pragma omp parallel
      {   vec<int> to_delete_loc;
          basevector E1,E2,t1,t2;
          kmer<K> x1,x2;
          vec< pair<int,int> > locs1, locs2;
          #pragma omp for schedule(dynamic,1) nowait
          for ( int v = 0; v < shb.N( ); v++ )
          {    for ( int pass = 1; pass <= 2; pass++ )
               {    const vec<int>& x = ( pass == 1 ? shb.To(v) : shb.From(v) );
                    for ( int i1 = 0; i1 < x.isize( ); i1++ )
                    for ( int i2 = i1 + 1; i2 < x.isize( ); i2++ )
                    {    int e1, e2;
                         if ( pass == 1 )
                         {    e1 = shb.EdgeObjectIndexByIndexTo( v, i1 );
                              e2 = shb.EdgeObjectIndexByIndexTo( v, i2 );    }
                         else
                         {    e1 = shb.EdgeObjectIndexByIndexFrom( v, i1 );
                              e2 = shb.EdgeObjectIndexByIndexFrom( v, i2 );    }
                         E1 = shb.EdgeObject(e1), E2 = shb.EdgeObject(e2);

                         // Pretest, just the first kmer.

                         x1.SetToSubOf( E1, ( pass == 1 ? E1.isize( ) - K : 0 ) );
                         x2.SetToSubOf( E2, ( pass == 1 ? E2.isize( ) - K : 0 ) );

                         t1.SetToSubOf(
                              E1, ( pass == 1 ? E1.isize( ) - K : K-poly ), poly );
                         t2.SetToSubOf(
                              E2, ( pass == 2 ? E2.isize( ) - K : K-poly ), poly );

                         if ( t1==chain0 || t1==chain1 || t1==chain2 || t1==chain3
                              || t2==chain0 || t2==chain1 || t2==chain2 || t2==chain3 ) 
                         {    continue;    }

                         typedef triple<kmer<K>,int,int> tri;
                         int64_t low1 = std::lower_bound( kmers_plus.begin( ), 
                              kmers_plus.end( ), x1,
                              []( tri a, kmer<K> val ) { return a.first < val; } )
                              - kmers_plus.begin( );
                         int64_t high1 = std::upper_bound( kmers_plus.begin( ), 
                              kmers_plus.end( ), x1,
                              []( kmer<K> val, tri a ) { return val < a.first; } )
                              - kmers_plus.begin( );
                         int64_t low2 = std::lower_bound( kmers_plus.begin( ), 
                              kmers_plus.end( ), x2,
                              []( tri a, kmer<K> val ) { return a.first < val; } )
                              - kmers_plus.begin( );
                         int64_t high2 = std::upper_bound( kmers_plus.begin( ), 
                              kmers_plus.end( ), x2,
                              []( kmer<K> val, tri a ) { return val < a.first; } )
                              - kmers_plus.begin( );

                         int q1 = 0, q2 = 0;
                         for ( auto i = low1; i < high1; i++ )
                         {    int id = kmers_plus[i].second, pos = kmers_plus[i].third;
                              if ( pass == 2 ) pos += K-1;
                              if ( id < (int64_t) quals_local.size( ) )
                                   q1 += quals_local[id][pos];
                              else
                              {    int idx = id - (int64_t) quals_local.size( );
                                   q1 += quals_local[idx]
                                        [ (int) quals_local[idx].size( ) - pos - 1 ];
                                             }    }
                         for ( auto i = low2; i < high2; i++ )
                         {    int id = kmers_plus[i].second, pos = kmers_plus[i].third;
                              if ( pass == 2 ) pos += K-1;
                              if ( id < (int64_t) quals_local.size( ) )
                                   q2 += quals_local[id][pos];
                              else
                              {    int idx = id - (int64_t) quals_local.size( );
                                   q2 += quals_local[idx]
                                        [ (int) quals_local[idx].size( ) - pos - 1 ];
                                             }    }
                         int count1 = high1 - low1, count2 = high2 - low2;
                         if ( q1 > q2 )
                         {    swap( e1, e2 );
                              swap( count1, count2 );
                              swap( q1, q2 );
                              swap( E1, E2 );    }
                         if (q2 < min_mult * q1 || q2 < min_win || q1 > max_lose) 
                              continue;

                         // Now go through K kmers, if possible.

                         int nk = std::min(
                              {K, shb.EdgeLengthKmers(e1), shb.EdgeLengthKmers(e2)} );
                         locs1.clear();
                         locs2.clear();
                         for ( int j = 0; j < nk; j++ )
                         {    x1.SetToSubOf( 
                                   E1, ( pass == 1 ? E1.isize( ) - K  - j: j ) );
                              x2.SetToSubOf( 
                                   E2, ( pass == 1 ? E2.isize( ) - K  - j: j ) );

                              int64_t low1 = std::lower_bound( kmers_plus.begin( ), 
                                   kmers_plus.end( ), x1,
                                   []( tri a, kmer<K> val ) { return a.first < val; } )
                                   - kmers_plus.begin( );
                              int64_t high1 = std::upper_bound( kmers_plus.begin( ), 
                                   kmers_plus.end( ), x1,
                                   []( kmer<K> val, tri a ) { return val < a.first; } )
                                   - kmers_plus.begin( );
                              int64_t low2 = std::lower_bound( kmers_plus.begin( ), 
                                   kmers_plus.end( ), x2,
                                   []( tri a, kmer<K> val ) { return a.first < val; } )
                                   - kmers_plus.begin( );
                              int64_t high2 = std::upper_bound( kmers_plus.begin( ), 
                                   kmers_plus.end( ), x2,
                                   []( kmer<K> val, tri a ) { return val < a.first; } )
                                   - kmers_plus.begin( );

                              for ( auto i = low1; i < high1; i++ )
                              {    int id = kmers_plus[i].second; 
                                   int pos = kmers_plus[i].third;
                                   pos += ( pass == 1 ? j : K-1-j );
                                   locs1.push( id, pos );    }
                              for ( auto i = low2; i < high2; i++ )
                              {    int id = kmers_plus[i].second; 
                                   int pos = kmers_plus[i].third;
                                   pos += ( pass == 1 ? j : K-1-j );
                                   locs2.push( id, pos );    }    }
                         UniqueSort(locs1), UniqueSort(locs2);
                         q1 = 0, q2 = 0;
                         for ( int i = 0; i < locs1.isize( ); i++ )
                         {    int id = locs1[i].first, pos = locs1[i].second;
                              if ( id < (int64_t) quals_local.size( ) )
                                   q1 += quals_local[id][pos];
                              else
                              {    int idx = id - (int64_t) quals_local.size( );
                                   q1 += quals_local[idx]
                                        [ (int) quals_local[idx].size( ) - pos - 1 ];
                                             }    }
                         for ( int i = 0; i < locs2.isize( ); i++ )
                         {    int id = locs2[i].first, pos = locs2[i].second;
                              if ( id < (int64_t) quals_local.size( ) )
                                   q2 += quals_local[id][pos];
                              else
                              {    int idx = id - (int64_t) quals_local.size( );
                                   q2 += quals_local[idx]
                                        [ (int) quals_local[idx].size( ) - pos - 1 ];
                                             }    }
                         count1 = locs1.size( ), count2 = locs2.size( );
                         if ( q1 > q2 )
                         {    swap( e1, e2 );
                              swap( count1, count2 );
                              swap( q1, q2 );    }
                         if ( count2 < min_win_ratio * Max(1,count1) ) continue;
                         if (q2 < min_mult * q1 || q2 < min_win || q1 > max_lose) 
                              continue;
                         if ( verbosity >= 1 )
                         {
                              #pragma omp critical
                              {    PRINT6(e1, e2, count1, count2, q1, q2);    }    }
                         to_delete_loc.push_back(e1);    }    }    }
          #pragma omp critical
          {    to_delete.append(to_delete_loc);    }
          to_delete_loc.clear( ), locs1.clear( ); locs2.clear( );    }    }

template <int K>
struct WingleFunctor
{
    void operator()( const SupportedHyperBasevector& shb,
                        const vecbasevector& bases_local,
                        const vecqualvector& quals_local,
                        vec<int>& to_delete, const int verbosity )
    { WingleSupport<K>(shb,bases_local,quals_local,to_delete,verbosity); }
};
    typedef std::unordered_map<int,int> exact_match_spec_t;
    void ExactSupports( const basevector& b, const qualvector& q
                      , const read_place& rp, const HyperBasevector& hb
                      , const exact_match_spec_t& match_spec
                      , vec<int>& matched_edges
                      , int iFlank=0
                      , const int min_qual=3) {
        iFlank=abs(iFlank);
        matched_edges.clear();
          
        const auto& edge_list=rp.E();
        bool bHasBubbleEdge=false;
        for(size_t ee=0;!bHasBubbleEdge&&ee<edge_list.size();++ee){
            bHasBubbleEdge=match_spec.find(edge_list[ee])!=match_spec.end();
        }
        if(!bHasBubbleEdge) return;
//        int qsum=0;
        int ei = 0, pos = rp.P( );
        
        int tgt_edge=-1;
        int tgt_pos=-1;
        auto itr=match_spec.find(edge_list[ei]);
        if(itr!=match_spec.end()){
            tgt_edge=(*itr).first;
            tgt_pos=(*itr).second;
        }
        for ( int l = 0; l < b.isize( ); l++ ){
            if(tgt_edge>=0 && tgt_edge==edge_list[ei] && tgt_pos==pos){
                const auto& tgt_edge_object=hb.EdgeObject(tgt_edge);
                bool bMatch= l-iFlank>=0 && l+iFlank < b.isize() && pos-iFlank >=0 && pos+iFlank<tgt_edge_object.isize();
                for(int ss=-iFlank ; bMatch && ss <= iFlank; ++ss){
                    bMatch = b[l+ss] == tgt_edge_object[tgt_pos+ss] ;
                }
                if(bMatch){
                    matched_edges.push_back(tgt_edge);
                }
            }
//            if ( b[l] != hb.EdgeObject( edge_list[ei] )[pos] )
//            {    if ( q[l] >= min_qual ) qsum += q[l] * 1000;
//                 else qsum += q[l];    }
            pos++;
            if ( pos == hb.EdgeObject( edge_list[ei] ).isize( ) ){
                ei++;
                itr=match_spec.find(edge_list[ei]);
                if(itr==match_spec.end()){
                    tgt_edge=-1;
                    tgt_pos=-1;
                }
                else{
                    tgt_edge=(*itr).first;
                    tgt_pos=(*itr).second;
                }
                if ( ei == rp.N( ) ) break;
                pos = hb.K( ) - 1;
            }
        }
//        ForceAssert(qsum==rp.Qsum());
    };


} // end of anonymous namespace

void SupportedHyperBasevector::DivineBubbles( const int L, const vecbasevector& bases, 
     const vecqualvector& quals, const long_heuristics& heur, const long_logging& logc )
{
     // New experimental approach to branches, in progress.

     if (logc.STATUS_LOGGING)
          cout << Date( ) << ": begin wingle branch approach" << endl;
     if ( logc.verb[ "WINGLE_BRANCH" ] >= 1 ) cout << "\n";
     vec<int> to_delete;
     {    vecbasevector bases_local;
          bases_local = bases;
          {    vecbasevector bases_rc(bases);
               for ( int64_t id = 0; id < (int64_t) bases.size( ); id++ )
                    bases_rc[id].ReverseComplement( );
               bases_local.Append(bases_rc);    }
          BigK::dispatch<WingleFunctor>(K(),*this, bases_local,quals,
                                          to_delete, logc.verb["WINGLE_BRANCH"]);    }
     if (logc.STATUS_LOGGING)
          cout << Date( ) << ": end wingle branch approach" << endl;

     // Find the bubbles.

     double clock1 = WallClockTime( );
     cout << Date( ) << ": divining bubbles" << endl;
     int verbosity = logc.verb[ "DIVINE_BUBBLES" ];
     vec< vec<int> > bubbles;

     // This is for reverse lookup to speed up place<->bubble intersection in the loop
     // bubbles_lookup[edge] stores a hashmap of bubbles associated with that edge
     // the hashmap coincidentally can be accessed as a list of pair<a,b>,
     // where a is bubble index, and b is 1 (a hit on bubble's index=1), 
     // 2(a hit on bubble's index=2), or 3 (hit on both index 1 and 2)
     // a bit wise & operation with 1 (01) and 2 (10) thus indicates the nature of the 
     // hit.

     vec<StdUnorderedMap<size_t,int> >  bubbles_lookup(EdgeObjectCount());
     for ( int v = 0; v < N( ); v++ )
     {    if ( To(v).size( ) != 1 || From(v).size( ) != 2 ) continue;
          if ( From(v)[0] != From(v)[1] ) continue;
          int w = From(v)[0];
          if ( To(w).size( ) != 2 || From(w).size( ) != 1 ) continue;
          int x = To(v)[0], y = From(w)[0];
          vec<int> all, b;
          all.push_back( x, v, w, y );
          UniqueSort(all);
          if ( all.size( ) != 4 ) continue;
          int e = EdgeObjectIndexByIndexTo( v, 0 );
          int f1 = EdgeObjectIndexByIndexFrom( v, 0 );
          int f2 = EdgeObjectIndexByIndexFrom( v, 1 );
          if ( f1 > f2 ) swap( f1, f2 );
          int g = EdgeObjectIndexByIndexFrom( w, 0 );
          b.push_back( e, f1, f2, g );
          if( size_t(std::max(f1,f2))>= bubbles_lookup.size())
          { bubbles_lookup.resize(2*std::max(f1,f2)); }
          bubbles_lookup[f1][bubbles.size()]+=1;
          bubbles_lookup[f2][bubbles.size()]+=2;
          bubbles.push_back(b);    }

     // Build data structures.

     if ( verbosity >= 1 ) cout << Date( ) << ": building data structures" << endl;
     HyperBasevector hb_fw(*this), hb_rc(*this);
     hb_rc.Reverse( );
     vec<int> to_right_fw, to_right_rc;
     hb_fw.ToRight(to_right_fw), hb_rc.ToRight(to_right_rc);
     vecbasevector x_fw, x_rc;
     for ( int i = 0; i < hb_fw.EdgeObjectCount( ); i++ )
          x_fw.push_back( hb_fw.EdgeObject(i) );
     for ( int i = 0; i < hb_rc.EdgeObjectCount( ); i++ )
          x_rc.push_back( hb_rc.EdgeObject(i) );
     VecIntPairVec locs_fw, locs_rc;
     CreateGlocs( x_fw, L, locs_fw );
     CreateGlocs( x_rc, L, locs_rc );
     if ( heur.DIVINE_MAX_LOCS >= 0 )
     {    for ( int64_t i = 0; i < (int64_t) locs_fw.size( ); i++ )
          {    if ( (int64_t) locs_fw[i].size( ) > heur.DIVINE_MAX_LOCS )
                    locs_fw[i].resize(0);    }
          for ( int64_t i = 0; i < (int64_t) locs_rc.size( ); i++ )
          {    if ( (int64_t) locs_rc[i].size( ) > heur.DIVINE_MAX_LOCS )
                    locs_rc[i].resize(0);    }    }

     // Define data structure adj.  It consists of data
     //      ( e1, e2, (w,q) )
     // arising from a read alignment, where (e1,e2) is a pair of consecutive
     // edges in the alignment, w = 1/n, where n is the number of same-scoring
     // placements of the read, and q is the quality score sum for the placement.  
     // Also define support[e], for each edge e.  It is the sum of w, as above.

     vec< triple< int, int, pair<fix64_6,int> > > adj;
     const size_t nEdgeObjectCount = EdgeObjectCount();
     vec<int64_t> support_raw(nEdgeObjectCount, 0); 
     //vec<fix64_6> support( EdgeObjectCount( ), 0 );

     // Go through the reads.

     const size_t nBubbles=bubbles.size();
     const int64_t nBases=bases.size();
     //this flattens the bubbles.size()-by-4 matrix row-by-row in to a 1D array
     const size_t uVotesBubWidth=4;
     vec<int64_t> votes_bub_raw(bubbles.size()*uVotesBubWidth,0); 
     //vec< vec<fix64_6> > votes_bub( bubbles.size( ), vec<fix64_6>(4,0) );

     //here, the fix64_6 encapsulation is broken to conform to openmp's 
     // atomic type requirement,
     //the following two scales are used for explicity manipulation of int64_t
     //the encapsulation can be recovered by changing fix64_6 appropiately

     const int64_t iFixedScale=1000000;
     const double  dFixedScale=1000000.0;

     vec< vec< pair<int,read_place> > > bplaces( bubbles.size( ) );
     vec< vec<read_place> > PLACES( bases.size( ) );
     if ( verbosity >= 1 ) cout << Date( ) << ": going through the reads" << endl;
     #pragma omp parallel
     {
          vec<read_place> places;

          vec< triple< int, int, pair<fix64_6,int> > > loc_adj;
          vec< vec< pair<int,read_place> > > loc_bplaces( bplaces.size() );
          #pragma omp for schedule(dynamic,1) nowait
          for ( int64_t id = 0; id < nBases; id++ )
          {    places.clear(); //vec<read_place> places;

               if ( verbosity >= 1 && id % 500000 == 0 )
               {    
                    #pragma omp critical
                    {    cout << Date( ) << ": processing read " << id 
                              << endl;    }    }

               if ( bases[id].isize( ) < L ) continue;
               int n = KmerId( bases[id], L, 0 );
               const int infinity = 1000000000;
               int qual_sum = infinity;
               FindPlaces( bases[id], quals[id], n, hb_fw, hb_rc, to_right_fw,
                    to_right_rc, locs_fw, locs_rc, places, qual_sum );
               if (heur.DIVINE_BRANCHES2) PLACES[id] = places;
               if (heur.DIVINE_STRONG)
               {    const int max_qual_sum = 100000;
                    if ( qual_sum > max_qual_sum ) continue;    }
               int np = places.size( );

//               #pragma omp critical
//               {
               for ( int i = 0; i < np; i++ )
               {
                    if ( verbosity >= 3 )
                    {
                         #pragma omp critical
                         {
                              for ( int j = 0; j < places[i].N( ); j++ )
                              {    cout << "read " << id << " supports edge "
                                        << places[i].E(j) << " with weight "
                                        << fix64_6(1,np) << endl;    }
                         }
                    }
                    if(places[i].Fw()){
                         for ( int j = 0; j < places[i].N( ) - 1; j++ )
                         {    loc_adj.push( places[i].E(j), places[i].E(j+1),
                                   make_pair( fix64_6(1,np), qual_sum ) );    }
                         for (auto edge: places[i].E()){
                             #pragma omp atomic
                             support_raw[ edge ] += iFixedScale/np; 
                             //support[ places[i].E(j) ] += fix64_6(1,np);

                             for( auto&  hit: bubbles_lookup[edge]){
                                 if( hit.second& 1){
                                      #pragma omp atomic
                                      votes_bub_raw[hit.first*uVotesBubWidth+0] 
                                           += iFixedScale/np; 
                                      // votes_bub[j][0] += fix64_6(1,np);
                                 }
                                 if( hit.second& 2){
                                      #pragma omp atomic
                                      votes_bub_raw[hit.first*uVotesBubWidth+2] 
                                           += iFixedScale/np; 
                                      // votes_bub[j][2] += fix64_6(1,np);
                                 }
                             }
                         }
                    }
                    else{
                         for ( int j = 0; j < places[i].N( ) - 1; j++ )
                         {    loc_adj.push( places[i].E(j+1), places[i].E(j),
                                   make_pair( fix64_6(1,np), qual_sum ) );    }
                         for (auto edge: places[i].E()){
                             #pragma omp atomic
                             support_raw[ edge ] += iFixedScale/np; 
                             //support[ places[i].E(j) ] += fix64_6(1,np);

                             for( auto&  hit: bubbles_lookup[edge]){
                                 if( hit.second& 1){
                                      #pragma omp atomic
                                      votes_bub_raw[hit.first*uVotesBubWidth+1] 
                                           += iFixedScale/np; 
                                      // votes_bub[j][1] += fix64_6(1,np);
                                 }
                                 if( hit.second& 2){
                                      #pragma omp atomic
                                      votes_bub_raw[hit.first*uVotesBubWidth+3] 
                                           += iFixedScale/np; 
                                      // votes_bub[j][3] += fix64_6(1,np);
                                 }
                             }
                         }
                    }
                    /*
                    #pragma omp critical
                    {    cout << "read " << id << ", placement " << i+1 << " of "
                              << places.size( ) << ": " << places[i] << "; meets bubbles "
                              << printSeq(bids) << endl;    }
                    */
               }

               if ( np == 1 || ( np == 2 && InvDef( places[0].E(0) ) ) )
               {
                    for ( int jp = 0; jp < places.isize( ); jp++ )
                    for (auto edge: places[jp].E()){
                        for( auto&  hit: bubbles_lookup[edge]){
                            if( hit.second& 1){
                                 loc_bplaces[hit.first].push( id, places[jp] );
                            }
                            if( hit.second& 2){
                                 loc_bplaces[hit.first].push( id, places[jp] );
                            }
                        }
                    }
               }
          }//for
          #pragma omp critical
          {
               adj.append(loc_adj);
               for(size_t ii=0;ii<bplaces.size();++ii){
                   bplaces[ii].append(loc_bplaces[ii]);
               }
          }
          loc_adj.clear();
          loc_bplaces.clear();
     }//omp
     if ( verbosity >= 3 ) cout << "\n";
     if ( verbosity >= 2 )
     {    cout << "\nedge support:\n";
          for ( int e = 0; e < EdgeObjectCount( ); e++ ){
               fix64_6 tmp;
               tmp.fromRaw(support_raw[e]);
               cout << e << " --> " << tmp << "\n";
//               cout << e << " --> " << support[e] << "\n";
          }
          cout << "\n";    }

     // Test branches.

     if ( verbosity >= 1 ) cout << Date( ) << ": testing branches" << endl;
     ParallelSort(adj);
     vec<int> start( EdgeObjectCount( ), -1 ), stop( EdgeObjectCount( ), -1 );
     for ( int i = 0; i < adj.isize( ); i++ )
     {    int j;
          for ( j = i + 1; j < adj.isize( ); j++ )
               if ( adj[j].first != adj[i].first ) break;
          start[ adj[i].first ] = i, stop[ adj[i].first ] = j;
          i = j - 1;    }
     for ( int v = 0; v < N( ); v++ )
     {    if ( To(v).size( ) != 1 || From(v).size( ) != 2 ) continue;
          // if ( From(v)[0] == From(v)[1] ) continue;
          int e = EdgeObjectIndexByIndexTo( v, 0 );
          int f1 = EdgeObjectIndexByIndexFrom( v, 0 );
          int f2 = EdgeObjectIndexByIndexFrom( v, 1 );
          if ( support_raw[f2] > support_raw[f1] ) swap(f1, f2);
          if ( support_raw[f2] > 1*iFixedScale ) continue;
          const int64_t i_min_win = 5*iFixedScale; //const fix64_6 min_win = 5;
          /*
          fix64_6 supp1 = 0, supp2 = 0;
          for ( int j = start[e]; j < stop[e]; j++ )
          {    if ( adj[j].second == f1 ) supp1 += adj[j].third.first;
               if ( adj[j].second == f2 ) supp2 += adj[j].third.first;    }
          if ( supp1 >= min_win ) to_delete.push_back(f2);    }
          */
          if ( support_raw[f1] >= i_min_win ) to_delete.push_back(f2);    } //if ( support[f1] >= min_win ) to_delete.push_back(f2);
     for ( int v = 0; v < N( ); v++ )
     {    if ( To(v).size( ) != 2 || From(v).size( ) != 1 ) continue;
          // if ( To(v)[0] == To(v)[1] ) continue;
          int e1 = EdgeObjectIndexByIndexTo( v, 0 );
          int e2 = EdgeObjectIndexByIndexTo( v, 1 );
          int f = EdgeObjectIndexByIndexFrom( v, 0 );
          if ( support_raw[e2] > support_raw[e1] ) swap(e1, e2);
          if ( support_raw[e2] > 1*iFixedScale ) continue;
          const int64_t i_min_win = 5*iFixedScale; //const fix64_6 min_win = 5;
          /*
          fix64_6 supp1 = 0, supp2 = 0;
          for ( int j = start[e1]; j < stop[e1]; j++ )
               supp1 += adj[j].third.first;
          for ( int j = start[e2]; j < stop[e2]; j++ )
               supp2 += adj[j].third.first;
          if ( supp1 >= min_win ) to_delete.push_back(e2);    }
          */
          if ( support_raw[e1] >= i_min_win ) to_delete.push_back(e2);    } //if ( support_raw[e1] >= min_win ) to_delete.push_back(e2);

     // Remove edges having no support.

     if (heur.DIVINE_STRONG)
     {    for ( int e = 0; e < EdgeObjectCount( ); e++ )
          {    if ( support_raw[e] == 0 )
               {    
                    // Don't delete edges that look like they're in a homopolymer
                    // cell.  Probably would be better to explicitly identify such
                    // cells, and avoid all edges in them.

                    const int minh = 20;
                    const basevector& b = EdgeObject(e);
                    int i;
                    for ( i = 1; i < minh; i++ )
                         if ( b[i] != b[0] ) break;
                    if ( i == minh ) continue;
                    for ( i = 1; i < minh; i++ )
                         if ( b[ b.isize( ) - i - 1 ] != b[ b.isize( ) - 1 ] ) break;
                    if ( i == minh ) continue;

                    if ( verbosity >= 1 ) 
                         cout << "deleting edge " << e << ", no support" << endl;
                    to_delete.push_back(e);    }    }    }

     // Test bubbles.

     const double max_asym_rarity = 0.00001;
     const int min_to_save = 10;
     if ( verbosity >= 1 ) cout << "\nbubbles:\n";
     for ( int i = 0; i < bubbles.isize( ); i++ )
     {    double f1 = votes_bub_raw[uVotesBubWidth*i+0]/dFixedScale, r1 = votes_bub_raw[uVotesBubWidth*i+1]/dFixedScale; //double f1 = votes_bub[i][0].ToDouble( ), r1 = votes_bub[i][1].ToDouble( );
          double f2 = votes_bub_raw[uVotesBubWidth*i+2]/dFixedScale, r2 = votes_bub_raw[uVotesBubWidth*i+3]/dFixedScale; //double f2 = votes_bub[i][2].ToDouble( ), r2 = votes_bub[i][3].ToDouble( );
          int e1 = 1, e2 = 2;
          if ( f2 + r2 > f1 + r1 || ( f2 + r2 == f1 + r1 && f2 > f1 ) )
          {    swap( e1, e2 );
               swap( f1, f2 );
               swap( r1, r2 );    }
          if ( f2 > r2 || ( f2 == r2 && f1 > r1 ) )
          {    swap( f1, r1 ); 
               swap( f2, r2 );    }
          long double p = Min( 0.5, f1/(f1+r1) ) / 2;
          int n = int(floor(f1+r1+f2+r2));
          double q = -1;
          if ( n > 0 && n <= 10000 ) 
          {    q = BinomialSum( n, int(ceil(f2)), p );
               if ( q < max_asym_rarity && f2 + r2 < min_to_save )
               {    if ( e1 == 1 ) to_delete.push_back( bubbles[i][2] );
                    else to_delete.push_back( bubbles[i][1] );    }    }
          if ( verbosity >= 1 )
          {    cout << "[" << i << "] " << bubbles[i][0] << " --> {" 
                    << bubbles[i][1] << "," << bubbles[i][2] << "} --> " 
                    << bubbles[i][3] << "; vote: " << setiosflags(ios::fixed) 
                    << setprecision(1) << f1 << "+" << r1 << " vs " << f2 << "+" 
                    << r2 << resetiosflags(ios::fixed);
               if ( q >= 0 ) cout << ", surprise(" << e2 << ") = " << q;
               cout << endl;    }    }
     if ( verbosity >= 1 ) cout << "\n";

     // Test bubbles, approach 2.  Consider a bubble that does not have indels in it.
     // Consider reads that place uniquely on one branch of the bubble, and score
     // their placement on the other branch.  The differences between the scores
     // are assigned to the winning branches.  These differences comprise two
     // distributions, and we check to see if one has a 'much' larger mean than
     // the other, and in that case delete the branch with the smaller mean.  To
     // measure the difference, we assign to each distribution the associated
     // normal distribution that reflects our expectation of its true mean.  Then
     // we take the difference of these two normal distributions.  If the mean of
     // the difference distribution is 3+ standard deviations away from zero, we 
     // deem one distribution to be much larger than the other.
     //
     // Note that the case where a path goes through a bubble twice is incorrectly 
     // handled.

     if ( verbosity >= 1 ) cout << Date( ) << ": bubble testing, approach 2" << endl;
     const double min_dist = 3.0;
     for ( int i = 0; i < bubbles.isize( ); i++ )
     {    int e1 = bubbles[i][1], e2 = bubbles[i][2];
          if ( verbosity >= 1 ) cout << "\nbubble " << i+1 << "\n";

          // Let's not allow indels.

          Bool indel = False;
          if ( EdgeObject(e1).size( ) != EdgeObject(e2).size( ) ) indel = True;
          else
          {    alignment a;
               int best_loc;
               SmithWatFree(EdgeObject(e1), EdgeObject(e2), best_loc, a, true, true);
               vec<int> mgg = a.MutationsGap1Gap2( EdgeObject(e1), EdgeObject(e2) );
               if ( mgg[1] > 0 || mgg[2] > 0 ) indel = True;    }
          if (indel)
          {    if ( verbosity >= 1 ) cout << "\nhas indel, ignoring" << "\n";
               continue;    }

          // Go through the placements.

          vec<int> q1, q2;
          for ( int j = 0; j < bplaces[i].isize( ); j++ )
          {    read_place z1 = bplaces[i][j].second, z2;
               int k1 = Position( z1.E( ), e1 ), k2 = Position( z1.E( ), e2 );
               if ( k1 >= 0 && k2 >= 0 ) continue;

               // Find the alternate placement.

               int id = -1;
               for ( int m = 0; m < 2; m++ )
               {    int p = Position( z1.E( ), bubbles[i][m+1] );
                    if ( p < 0 ) continue;
                    z2 = z1;
                    z2.EMutable( )[p] = bubbles[i][2-m];
                    id = bplaces[i][j].first;
                    const int min_qual = 3;
                    z2.ComputeQsum( bases[id], quals[id],
                         ( z2.Fw( ) ? hb_fw : hb_rc ), min_qual );
                    break;    }

               // Print.

               if ( k1 < 0 ) swap( z1, z2 );
               if ( verbosity >= 1 )
               {    cout << "\nid = " << id << ", place1 = " << z1 << "\n";
                    cout << "id = " << id << ", place2 = " << z2 << "\n";    }
          
               // Save quality score sum differences.

               if ( z1.Qsum( ) < z2.Qsum( ) )
                    q1.push_back( z2.Qsum( ) - z1.Qsum( ) );
               if ( z2.Qsum( ) < z1.Qsum( ) )
                    q2.push_back( z1.Qsum( ) - z2.Qsum( ) );    }

          // Tally.

          if ( q1.empty( ) && q2.empty( ) ) 
          {    if ( verbosity >= 1 ) cout << "\nno support\n";    }
          else
          {    Sort(q1), Sort(q2);
               double m1 = -1, m2 = -1;
               if ( q1.nonempty( ) ) m1 = Mean(q1);
               if ( q2.nonempty( ) ) m2 = Mean(q2);
               int n1 = q1.size( ), n2 = q2.size( );
               if ( verbosity >= 1 )
               {    cout << setiosflags(ios::fixed) << setprecision(1);
                    cout << "\nedge " << e1 << " is supported by " << n1 << " reads";
                    if ( m1 >= 0 )
                    {    cout << ", mean qsum diff = " << setprecision(1)
                              << m1 / 1000.0 << "\n";    }
                    cout << "edge " << e2 << " is supported by " << n2 << " reads";
                    if ( m2 >= 0 )
                    {    cout << ", mean qsum diff = " << setprecision(1)
                              << m2 / 1000.0 << "\n";    }    }
               if ( q1.empty( ) )
               {    if ( verbosity >= 1 ) cout << "deleting edge " << e1 << "\n";
                    to_delete.push_back(e1);
                    continue;    }
               if ( q2.empty( ) )
               {    if ( verbosity >= 1 ) cout << "deleting edge " << e2 << "\n";
                    to_delete.push_back(e2);
                    continue;    }
               double dist = (m1-m2) / sqrt( m1*m1/n1 + m2*m2/n2 );
               if ( verbosity >= 1 )
               {    cout << "dist = " << setprecision(2) << Abs(dist) << "\n";
                    cout << resetiosflags(ios::fixed);    }
               if ( dist <= -min_dist ) 
               {    if ( verbosity >= 1 ) cout << "deleting edge " << e1 << "\n";
                    to_delete.push_back(e1);    }
               if ( dist >= min_dist )
               {    if ( verbosity >= 1 ) cout << "deleting edge " << e2 << "\n";
                    to_delete.push_back(e2);    }    }    }
     if ( verbosity >= 1 ) cout << "\n";

     // Test branches in a manner similar to the above.

     if (heur.DIVINE_BRANCHES)
     {
     const int max_kill2 = 100000;
     const int max_branchiness = 4;
     int bverbosity = 2;
     vec< triple<int,int,int> > branches;
     for ( int v = 0; v < N( ); v++ )
     {    if ( From(v).isize( ) <= max_branchiness )
          {    for ( int j1 = 0; j1 < From(v).isize( ); j1++ )
               for ( int j2 = j1 + 1; j2 < From(v).isize( ); j2++ )
               {    branches.push( EdgeObjectIndexByIndexFrom( v, j1 ),
                         EdgeObjectIndexByIndexFrom( v, j2 ), 0 );    }    }
          if ( To(v).isize( ) <= max_branchiness )
          {    for ( int j1 = 0; j1 < To(v).isize( ); j1++ )
               for ( int j2 = j1 + 1; j2 < To(v).isize( ); j2++ )
               {    branches.push( EdgeObjectIndexByIndexTo( v, j1 ),
                         EdgeObjectIndexByIndexTo( v, j2 ), 1 );    }    }    }
     UniqueSort(branches);
     vec< vec<read_place> > PLACES2( bases.size( ) );
     const double prox = 50 * 1000;
     #pragma omp parallel for
     for ( int64_t id = 0; id < nBases; id++ )
     {    int n = KmerId( bases[id], L, 0 );
          vec<read_place> places;
          const int infinity = 1000000000;
          int qual_sum = infinity;
          const int min_qual = 1;
          FindPlaces( bases[id], quals[id], n, hb_fw, hb_rc, to_right_fw,
               to_right_rc, locs_fw, locs_rc, places, qual_sum, min_qual, prox );

          if ( bverbosity >= 2 )
          {    
               int best = 1000000000;
               {    for ( int j = 0; j < places.isize( ); j++ )
                         best = Min( best, places[j].Qsum( ) );    }
               #pragma omp critical
               {    for ( int j = 0; j < places.isize( ); j++ )
                    {    cout << "read " << id << " placed: " << places[j] 
                              << ( places[j].Qsum( ) == best ? " [BEST]" : "" )
                              << "\n";    }    }    }

          PLACES2[id] = places;    }

     vec< vec< pair<int,double> > > homes( bases.size( ) );
     for ( int64_t id = 0; id < nBases; id++ )
     {    vec<int> v;
          for ( int j = 0; j < PLACES2[id].isize( ); j++ )
          for ( int l = 0; l < PLACES2[id][j].N( ); l++ )
               v.push_back( PLACES2[id][j].E(l) );
          UniqueSort(v);
          for ( int i = 0; i < v.isize( ); i++ )
          {    int e = v[i];
               int q1 = 1000000000, q2 = 1000000000;
               for ( int j = 0; j < PLACES2[id].isize( ); j++ )
               {    if ( Member( PLACES2[id][j].E( ), e ) )
                         q1 = Min( q1, PLACES2[id][j].Qsum( ) );
                    else q2 = Min( q2, PLACES2[id][j].Qsum( ) );    }
               if ( bverbosity >= 2 ) PRINT3( e, id, (q2-q1)/1000.0 );
               homes[id].push( e, q2 - q1 );    }    }

     vec< vec< pair<int,int> > > homes_index( EdgeObjectCount( ) );
     for ( int64_t id = 0; id < nBases; id++ )
     for ( int j = 0; j < homes[id].isize( ); j++ )
          homes_index[ homes[id][j].first ].push( id, homes[id][j].second );

     for ( int i = 0; i < branches.isize( ); i++ )
     {    int e1 = branches[i].first, e2 = branches[i].second;
          int dir = branches[i].third;

          const int qgood = 30;
          vec<int> count(2, 0);

          for ( int r = 0; r < 2; r++ )
          {    int e = ( r == 0 ? e1 : e2 );
               for ( int j = 0; j < homes_index[e].isize( ); j++ )
                    if ( homes_index[e][j].second/1000.0 >= qgood ) count[r]++;    }

          const double min_ratio = 10.0;
          const double max_weak = 2;

          for ( int r = 0; r < 2; r++ )
          {    int e = ( r == 0 ? e1 : e2 ), eop = ( r == 0 ? e2 : e1 );
               if ( count[r] >= min_ratio * Max( 1, count[1-r] )
                    && count[1-r] <= max_weak )
               {    if ( bverbosity >= 1 )
                    {    cout << "\n";
                         PRINT4( e1, count[0], e2, count[1] );
                         cout << "deleting edge " << eop << "\n";        }
                    to_delete.push_back(eop);    }    }    }
     cout << "\n";
     }

     if (heur.DIVINE_BRANCHES2)
     {
     const int min_qual = 1;
     int bverbosity = 2;
     vec< triple<int,int,int> > branches;
     for ( int v = 0; v < N( ); v++ )
     {    for ( int j1 = 0; j1 < From(v).isize( ); j1++ )
          for ( int j2 = j1 + 1; j2 < From(v).isize( ); j2++ )
          {    branches.push( EdgeObjectIndexByIndexFrom( v, j1 ),
                    EdgeObjectIndexByIndexFrom( v, j2 ), 0 );    }
          for ( int j1 = 0; j1 < To(v).isize( ); j1++ )
          for ( int j2 = j1 + 1; j2 < To(v).isize( ); j2++ )
          {    branches.push( EdgeObjectIndexByIndexTo( v, j1 ),
                    EdgeObjectIndexByIndexTo( v, j2 ), 1 );    }    }
     UniqueSort(branches);
     vec< vec< triple<int,int,int> > > places_index( EdgeObjectCount( ) );
     for ( int i = 0; i < PLACES.isize(); i++ )
     for ( int j = 0; j < PLACES[i].isize(); j++ )
     for ( int k = 0; k < PLACES[i][j].E( ).isize(); k++ )
          places_index[ PLACES[i][j].E(k) ].push( i, j, k );
     for ( int i = 0; i < branches.isize( ); i++ )
     {    int e1 = branches[i].first, e2 = branches[i].second;
          int dir = branches[i].third;
          vec<double> q1, q2;
          for ( int r = 0; r < 2; r++ )
          {    
               int e = ( r == 0 ? e1 : e2 );
               int eop = ( r == 0 ? e2 : e1 );

               for ( int l = 0; l < places_index[e].isize( ); l++ )
               {    int id = places_index[e][l].first;
                    int pi = places_index[e][l].second;
                    int pos = places_index[e][l].third;
                    read_place p = PLACES[id][pi];
                    if ( dir == 0 && pos == 0 ) continue;
                    if ( dir == 1 && pos == p.N( ) - 1 ) continue;

                    // Now find the best placements of the read, assuming it uses
                    // eop at the given place.  We assume that the read has been
                    // placed up to the appropriate point.

                    vec<read_place> places_op;
                    int qual_sum_op = 1000000000;
                    if ( ( dir == 0 && p.Fw( ) ) || ( dir == 1 && !p.Fw( ) ) )
                    {    p.EMutable( )[pos] = eop;
                         p.EMutable( ).resize( pos + 1 );
                         if ( p.P( ) >= EdgeLengthBases( p.E(0) ) ) continue;
                         p.ComputeQsum( bases[id], quals[id],
                              ( p.Fw( ) ? hb_fw : hb_rc ), min_qual );
                         ExtendPlacement( p, bases[id], quals[id], hb_fw, hb_rc, 
                              to_right_fw, to_right_rc, places_op, qual_sum_op, 
                              min_qual );    }
                    else
                    {    basevector br( bases[id] );
                         qualvector qr( quals[id] );
                         br.ReverseComplement( );
                         qr.ReverseMe( );
                         p.Reverse( br, *this );
                         int posr = p.N( ) - pos - 1;
                         p.EMutable( )[posr] = eop;
                         p.EMutable( ).resize( posr + 1 );
                         if ( p.P( ) >= EdgeLengthBases( p.E(0) ) ) continue;
                         p.ComputeQsum( br, qr,
                              ( p.Fw( ) ? hb_fw : hb_rc ), min_qual );
                         ExtendPlacement( p, br, qr, hb_fw, hb_rc, to_right_fw, 
                              to_right_rc, places_op, qual_sum_op, min_qual );
                         for ( int j = 0; j < places_op.isize( ); j++ )
                              places_op[j].Reverse( br, *this );    }

                    read_place px = PLACES[id][pi];
                    px.ComputeQsum( bases[id], quals[id],
                         ( px.Fw( ) ? hb_fw : hb_rc ), min_qual );
                    ForceAssert( places_op.nonempty( ) );
                    double delta = places_op[0].Qsum( ) - px.Qsum( );
                    if ( r == 1 ) delta = -delta;
                    if ( delta > 0 ) q1.push_back(delta);
                    else if ( delta < 0 ) q2.push_back(-delta);

                    if (bverbosity >= 2 )
                    {    cout << "\n";
                         cout << "read " << id << "\n";
                         PRINT2( e1, e2 );
                         cout << px << "\n";
                         for ( int j = 0; j < places_op.isize( ); j++ )
                              cout << places_op[j] << "\n";    }    }    }

          Sort(q1), Sort(q2);
          double m1 = -1, m2 = -1;
          if ( q1.nonempty( ) ) m1 = Mean(q1);
          if ( q2.nonempty( ) ) m2 = Mean(q2);
          int n1 = q1.size( ), n2 = q2.size( );
          if ( bverbosity >= 1 )
          {    cout << setiosflags(ios::fixed) << setprecision(1);
               cout << "\nedge " << e1 << " is supported by " << n1 << " reads";
               if ( m1 >= 0 )
               {    cout << ", mean qsum diff = " << setprecision(1)
                         << m1 / 1000.0;    }
               cout << "\n";
               cout << "edge " << e2 << " is supported by " << n2 << " reads";
               if ( m2 >= 0 )
               {    cout << ", mean qsum diff = " << setprecision(1)
                         << m2 / 1000.0;    }    
               cout << "\n";    }    }
     cout << "\n";
     }

     // Clean up graph.  Note that inverse of edges have to be added in.  I'm not
     // sure why these weren't found in the first place.

     int ndels = to_delete.size( );
     for ( int i = 0; i < ndels; i++ )
          if ( Inv( to_delete[i] ) >= 0 ) to_delete.push_back( Inv( to_delete[i] ) );
     DeleteEdges(to_delete);
     DeleteUnusedPaths( );
     REPORT_TIME( clock1, "used divining bubbles" );
     RemoveDeadEdgeObjects( );
     double clock2 = WallClockTime( );
     RemoveEdgelessVertices( );
     RemoveUnneededVertices( );
     REPORT_TIME( clock2, "used divining bubbles tail" );
     RemoveDeadEdgeObjects( );
     FixWeights(logc);
     TestValid(logc);    
     if (heur.DIVINE_STRONG)
     {    const double junk_ratio = 10.0;
          const int max_del = 1000;
          TrimHangingEnds( max_del, junk_ratio, heur, logc );
          RemoveSmallMostlyAcyclicComponents(logc);
          TestValid(logc);    }
     DeleteReverseComplementComponents(logc);    }

void SupportedHyperBasevector::FixWeights( const long_logging& logc )
{    double clock = WallClockTime( );

     for ( int i1 = 0; i1 < NPaths( ); i1++ )
     {    const vec<int>& p1 = Path(i1);
          if ( !InvDef( p1[0] ) ) continue;
          vec<int> p2;
          for ( int j = 0; j < p1.isize( ); j++ )
          {     if ( p1[j] < 0 ) p2.push_back( p1[j] );
               else p2.push_back( Inv( p1[j] ) );    }
          p2.ReverseMe( );
          int i2 = BinPosition( Paths( ), p2 );
          if ( i2 < 0 )
          {    cout << "\nAttempting to fix weights, found asymmetric path." << endl;
               cout << "path = " << printSeq( Path(i1) ) << endl;
               cout << "inv path = " << printSeq(p2) << endl << "Abort." << endl;
               TracebackThisProcess( );    }
          if ( WeightFw(i1) != WeightRc(i2) )
          {    fix64_6 w = Max( WeightFw(i1), WeightRc(i2) );
               if ( w > WeightRc(i2) ) 
               {    WeightRcMutable(i2) = w;
                    // Temporary if.
                    if ( WeightsFwOrigin( ).size( ) == WeightsFw( ).size( ) )
                         WeightRcOriginMutable(i2) = WeightFwOriginMutable(i1);    }
               else 
               {    WeightFwMutable(i1) = w;    
                    // Temporary if.
                    if ( WeightsFwOrigin( ).size( ) == WeightsFw( ).size( ) )
                         WeightFwOriginMutable(i1) = WeightRcOriginMutable(i2);    }
                         }
          if ( WeightFw(i2) != WeightRc(i1) )
          {    fix64_6 w = Max( WeightFw(i2), WeightRc(i1) );
               if ( w > WeightRc(i1) )
               {    WeightRcMutable(i1) = w;    
                    // Temporary if.
                    if ( WeightsFwOrigin( ).size( ) == WeightsFw( ).size( ) )
                         WeightRcOriginMutable(i1) = WeightFwOriginMutable(i2);    }
               else 
               {    WeightFwMutable(i2) = w;    
                    // Temporary if.
                    if ( WeightsFwOrigin( ).size( ) == WeightsFw( ).size( ) )
                         WeightFwOriginMutable(i2) = WeightRcOriginMutable(i1);    }
                         }    }

     for ( int i1 = 0; i1 < NPairs( ); i1++ )
     {    const vec<int> &p1 = PairLeft(i1), &q1 = PairRight(i1);
          if ( !InvDef( p1[0] ) || !InvDef( q1[0] ) ) continue;
          vec<int> p2, q2;
          for ( int j = 0; j < p1.isize( ); j++ )
          {     if ( p1[j] < 0 ) p2.push_back( p1[j] );
               else p2.push_back( Inv( p1[j] ) );    }
          for ( int j = 0; j < q1.isize( ); j++ )
          {     if ( q1[j] < 0 ) q2.push_back( q1[j] );
               else q2.push_back( Inv( q1[j] ) );    }
          p2.ReverseMe( ), q2.ReverseMe( );
          int i2 = BinPosition( Pairs( ), make_pair( q2, p2 ) );
          if ( i2 < 0 )
          {    cout << "\nAttempting to fix weights, found asymmetric pair." << endl;
               cout << "pair = " << printSeq( PairLeft(i1) ) << " ... "
                    << printSeq( PairRight(i1) ) << endl;
               cout << "inv path = " << printSeq(q2) << " ... " 
                    << printSeq(p2) << endl << "Abort." << endl;
               TracebackThisProcess( );    }
          if ( PairData(i1) != PairData(i2) )
          {    PairDataMutable(i2) = PairData(i1);    }    }

     REPORT_TIME( clock, "used fixing weights" );    }

void SupportedHyperBasevector::MakeLoops( const long_logging& logc )
{    double clock = WallClockTime( );
     if (logc.STATUS_LOGGING) cout << Date( ) << ": making loops" << endl;
     vec<int> to_left, to_right;
     ToLeft(to_left), ToRight(to_right);
     for ( int v = 0; v < N( ); v++ )
     {    if ( !From(v).solo( ) || To(v).size( ) != 2 ) continue;
          int w = From(v)[0];
          if ( From(w).size( ) != 2 || !To(w).solo( ) ) continue;
          int i1 = 0, i2 = 1;
          if ( From(w)[i1] != v ) swap( i1, i2 );
          if ( From(w)[i1] != v ) continue;
          int e1 = EdgeObjectIndexByIndexFrom( v, 0 );
          int e2 = EdgeObjectIndexByIndexFrom( w, i1 );
          int e0 = EdgeObjectIndexByIndexTo( v, 0 );
          if ( e0 == e2 ) e0 = EdgeObjectIndexByIndexTo( v, 1 );
          int u = to_left[e0];

          // u --e0--> v --e1--> w --*--> *
          //             <--e2--

          basevector left = Cat( e0, e1 ), loop = Cat( e2, e1 );
          int re0 = Inv(e0), re1 = Inv(e1), re2 = Inv(e2);
          vec<int> dels, rdels;
          dels.push_back( e0, e1, e2 );
          rdels.push_back( re0, re1, re2 );
          Sort(dels), Sort(rdels);
          if ( re0 >= 0 && Meet( dels, rdels ) ) continue;
          if ( logc.verb[ "MAKE_LOOPS" ] >= 1 )
               cout << "making loop from " << e1 << ", " << e2 << endl;
          DeleteEdges(dels);
          int e21 = EdgeObjectCount( ), e01 = EdgeObjectCount( ) + 1;
          int re21 = EdgeObjectCount( ) + 2, re01 = EdgeObjectCount( ) + 3;
          AddEdge(w, w, loop), AddEdge(u, w, left);
          to_left.push_back(w,u), to_right.push_back(w,w);
          if ( re0 < 0 ) InvMutable( ).push_back( -1, -1 );
          else
          {    int ne = EdgeObjectCount( );
               InvMutable( ).push_back( re21, re01, e21, e01 );
               DeleteEdges(rdels);
               int ru = to_right[re0], rv = to_right[re1], rw = to_left[re1];
               basevector rloop(loop), rleft(left);
               rloop.ReverseComplement( ), rleft.ReverseComplement( );
               AddEdge(rw, rw, rloop), AddEdge(rw, ru, rleft);
               to_left.push_back(rw,rw), to_right.push_back(rw,ru);    }
          vec< vec<int> > e(2);
          e[0].push_back( e2, e1 );
          e[1].push_back( e0, e1 );
          vec<int> f;
          f.push_back( e21, e01 );
          TransformPaths( e, f );
          if ( re0 >= 0 )
          {    vec< vec<int> > e(2);
               e[0].push_back( re1, re2 );
               e[1].push_back( re1, re0 );
               vec<int> f;
               f.push_back( re21, re01 );
               TransformPaths( e, f );    }    }
     UniqueOrderPaths( );
     RemoveEdgelessVertices( );
     REPORT_TIME( clock, "used in MakeLoops" );
     RemoveDeadEdgeObjects( );
     if (logc.STATUS_LOGGING)
     {    cout << Date( ) << ": after making loops there are "
               << EdgeObjectCount( ) << " edges" << endl;    }
     TestValid(logc);    }

void SupportedHyperBasevector::DivineSingleMutationBubbles( const int L, const vecbasevector& bases, 
     const vecqualvector& quals, const long_heuristics& heur, const long_logging& logc )
{
     vec<int> to_delete;
    
     typedef vec<StdUnorderedMap<size_t,int> >  bubble_lookup_t;
     
     class report_t{
     private:
         const bubble_lookup_t& lookup;
         vec<double> _num;
         std::set<int> _bubbles_hit;
         double _den;
         size_t _nEntries;
         const size_t nBubbles;
         int _factor;
         bool _bLogged;
         enum{_WIDTH=4};
         report_t();
     public:
         explicit report_t(const bubble_lookup_t& l,const size_t nB):lookup(l),_num(nB*_WIDTH,0),_den(0), _nEntries(0), nBubbles(nB), _factor(std::numeric_limits<int>::max()),_bLogged(false){};
         double prob(size_t ii)const{ return _num[ii]/_den;};
         const vec<double>& numerators()const{ return _num;};
         const std::set<int>& activeBubbles()const{return _bubbles_hit;};
         double denominator()const{return _den;};
         uint64_t nEntries()const{return _nEntries;};
         bool logged()const{return _bLogged;};
         
         void insert(const read_place& rp, const vec<int>& matched_edges){
             const int iQ=rp.Qsum();
             if( _factor == std::numeric_limits<int>::max()){
                 _factor = iQ;
             }
             else if(iQ<_factor){
                 const double dScale = pow(10,0.0001*(iQ-_factor));
                 for(auto&entry:_num){entry*=dScale;}
                 _den*=dScale;
                 _factor=iQ;
             }
             const double pp = pow(10,-0.0001*(iQ-_factor));
             for(const auto& edge: matched_edges){
                 for(const auto& hit:lookup[edge]){
                     _bubbles_hit.insert(hit.first);
                 }
             }
             const size_t offset = !rp.Fw();
             for(const auto&edge:rp.E()){
                 for(const auto& hit:lookup[edge]){
                     _num[hit.first*_WIDTH+offset+hit.second]+=pp;
                     _bLogged=true;
                     /*
                     if(hit.second&1){
                         _num[hit.first*_WIDTH+offset]+=pp;
                         _bLogged=true;
                     }
                     if(hit.second&2){
                         _num[hit.first*_WIDTH+offset+2]+=pp;
                         _bLogged=true;
                     }
                     */
                 }
             }
             _den+=pp;
             ++_nEntries;
         };
         void reset(){
             _bubbles_hit.clear();
             _num.assign(_num.size(),0.0);
             _den=0;
             _nEntries=0;
             _bLogged=false;
             _factor=std::numeric_limits<int>::max();
         }
     };
     
     class read_place_buffer{
     private:
         vec<read_place> _container;
         size_t _nElements;
     public:
         //typedef typename vec<read_place>::value_type                value_type;
         //typedef typename vec<read_place>::reference                 reference;
         //typedef typename vec<read_place>::const_reference           const_reference;
         //typedef typename vec<read_place>::size_type                 size_type;
         read_place_buffer():_container(),_nElements(0){_container.reserve(10000000);};
         
         size_t size()const{return _nElements;};
         bool nonempty()const{return _nElements>0;};
         void reserve(size_t ii){_container.reserve(ii);};
         
         void clear(){_nElements=0;};
         void push_back(const read_place& rp){
             if(_nElements >= _container.size()){
                 _container.resize(2*_container.size()+1);
             }
             _container[_nElements]=rp;
             ++_nElements;
         };
         void pop_back(){ --_nElements; };
         
         read_place& front()                         { return _container[0]; };
         const read_place& front()const              { return _container[0]; };
         read_place& back()                          { return _container[_nElements-1]; };
         const read_place& back()const               { return _container[_nElements-1]; };
         read_place& operator[](size_t ii)           { return _container[ii]; };
         const read_place& operator[](size_t ii)const{ return _container[ii]; }
         
         typename vec<read_place>::iterator begin() noexcept{return _container.begin();};
         typename vec<read_place>::const_iterator begin() const noexcept{return _container.begin();};
         typename vec<read_place>::iterator end() noexcept{return _container.begin()+_nElements;};
         typename vec<read_place>::const_iterator end() const noexcept{return _container.begin()+_nElements;};
     };
     auto loc_FindPlacesInt=[] ( const basevector& b, const qualvector& q, const int n, 
                              const HyperBasevector& hb_fw, const HyperBasevector& hb_rc, 
                              const vec<int>& to_right_fw, const vec<int>& to_right_rc, 
                              const VecIntPairVec& locs_fw, const VecIntPairVec& locs_rc,
                              report_t& report, int& qual_sum, const int min_qual, const double prox,
                              const exact_match_spec_t& exact_match_spec,
                              read_place_buffer&places_part )
     {    
          report.reset();
          vec<int> matched_edges;
          int best_qsum = qual_sum;
          for ( int opass = 1; opass <= 2; opass++ )
          {    const HyperBasevector& hb = ( opass == 1 ? hb_fw : hb_rc );
               const VecIntPairVec& locs = ( opass == 1 ? locs_fw : locs_rc );
               const vec<int>& to_right = ( opass == 1 ? to_right_fw : to_right_rc );
               places_part.clear();
               for ( unsigned j = 0; j < locs[n].size( ); j++ )
               {    
                    // Test to see if we're starting in the last K-1 bases of the 
                    // edge.  Note that we could have built the locs so that this would
                    // be already done.

                    if ( locs[n][j].second >= hb.EdgeLengthKmers( locs[n][j].first ) )
                         continue;

                    vec<int> v(1);
                    v[0] = locs[n][j].first;
                    read_place p( v, locs[n][j].second, opass == 1, 0 );
                    p.ComputeQsum125( b, q, hb, min_qual );
                    places_part.push_back(p);    }
               while( places_part.nonempty( ) )
               {    read_place p = places_part.back( );
                    places_part.pop_back( );
                    if ( p.Qsum( ) > best_qsum + prox ) continue;
                    int ext1 = b.isize( ) - ( hb.EdgeLengthBases( p.E(0) ) - p.P( ) );
                    for ( int j = 1; j < p.N( ); j++ )
                         ext1 -= hb.EdgeLengthKmers( p.E(j) );
                    if ( ext1 <= 0 ) 
                    {    if ( p.Qsum( ) <= best_qsum + prox )
                         {    ExactSupports(b, q, p, hb, exact_match_spec, matched_edges, 5, min_qual);
                              report.insert(p,matched_edges);
                              if ( p.Qsum( ) <= best_qsum ) {
                                   best_qsum = p.Qsum( );    }    }    }
                    else
                    {    int v = to_right[ p.E( p.N( ) - 1 ) ];
                         for ( int j = 0; j < hb.From(v).isize( ); j++ )
                         {    int e = hb.EdgeObjectIndexByIndexFrom( v, j );
                                   places_part.push_back(p);
                                   places_part.back().AddEdge125( e, b, q, hb, min_qual );
                                   if ( places_part.back().Qsum( ) > best_qsum + prox )
                                        places_part.pop_back(); }    }    }    }
          qual_sum = best_qsum;
     };
     
     double clock1 = WallClockTime( );
     cout << Date( ) << ": divining single mutation bubbles" << endl;
     int verbosity = logc.verb[ "DIVINE_BUBBLES" ];
     vec< vec<int> > bubbles;

     // This is for reverse lookup to speed up place<->bubble intersection in the loop
     // bubbles_lookup[edge] stores a hashmap of bubbles associated with that edge
     // the hashmap coincidentally can be accessed as a list of pair<a,b>,
     // where a is bubble index, and b is 1 (a hit on bubble's index=1), 
     // 2(a hit on bubble's index=2), or 3 (hit on both index 1 and 2)
     // a bit wise & operation with 1 (01) and 2 (10) thus indicates the nature of the 
     // hit.

     // Find the bubbles.
     bubble_lookup_t  bubbles_lookup(EdgeObjectCount());
     for ( int v = 0; v < N( ); v++ )
     {    if ( To(v).size( ) != 1 || From(v).size( ) != 2 ) continue;
          if ( From(v)[0] != From(v)[1] ) continue;
          int w = From(v)[0];
          if ( To(w).size( ) != 2 || From(w).size( ) != 1 ) continue;
          int x = To(v)[0], y = From(w)[0];
          vec<int> all, b;
          all.push_back( x, v, w, y );
          UniqueSort(all);
          if ( all.size( ) != 4 ) continue;
          int e = EdgeObjectIndexByIndexTo( v, 0 );
          int f1 = EdgeObjectIndexByIndexFrom( v, 0 );
          int f2 = EdgeObjectIndexByIndexFrom( v, 1 );
          if(EdgeObject(f1).isize() !=  1 + 2 * (K()-1) ) continue;
          if(EdgeObject(f2).isize() !=  1 + 2 * (K()-1) ) continue;
          if ( f1 > f2 ) swap( f1, f2 );
          int g = EdgeObjectIndexByIndexFrom( w, 0 );
          b.push_back( e, f1, f2, g );
          if( size_t(std::max(f1,f2))>= bubbles_lookup.size())
          { bubbles_lookup.resize(2*std::max(f1,f2)); }
//          bubbles_lookup[f1][bubbles.size()]+=1;
//          bubbles_lookup[f2][bubbles.size()]+=2;
          bubbles_lookup[f1][bubbles.size()]=0;
          bubbles_lookup[f2][bubbles.size()]=2;
          bubbles.push_back(b);    }

     // Build data structures.

     if ( verbosity >= 1 ) cout << Date( ) << ": building data structures" << endl;
     HyperBasevector hb_fw(*this), hb_rc(*this);
     hb_rc.Reverse( );
     vec<int> to_right_fw, to_right_rc;
     hb_fw.ToRight(to_right_fw), hb_rc.ToRight(to_right_rc);
     vecbasevector x_fw, x_rc;
     for ( int i = 0; i < hb_fw.EdgeObjectCount( ); i++ )
          x_fw.push_back( hb_fw.EdgeObject(i) );
     for ( int i = 0; i < hb_rc.EdgeObjectCount( ); i++ )
          x_rc.push_back( hb_rc.EdgeObject(i) );
     VecIntPairVec locs_fw, locs_rc;
     CreateGlocs( x_fw, L, locs_fw );
     CreateGlocs( x_rc, L, locs_rc );
     if ( heur.DIVINE_MAX_LOCS >= 0 )
     {    for ( int64_t i = 0; i < (int64_t) locs_fw.size( ); i++ )
          {    if ( (int64_t) locs_fw[i].size( ) > heur.DIVINE_MAX_LOCS )
                    locs_fw[i].resize(0);    }
          for ( int64_t i = 0; i < (int64_t) locs_rc.size( ); i++ )
          {    if ( (int64_t) locs_rc[i].size( ) > heur.DIVINE_MAX_LOCS )
                    locs_rc[i].resize(0);    }    }

     // Define data structure adj.  It consists of data
     //      ( e1, e2, (w,q) )
     // arising from a read alignment, where (e1,e2) is a pair of consecutive
     // edges in the alignment, w = 1/n, where n is the number of same-scoring
     // placements of the read, and q is the quality score sum for the placement.  
     // Also define support[e], for each edge e.  It is the sum of w, as above.

//     vec< triple< int, int, pair<fix64_6,int> > > adj;
     const size_t nEdgeObjectCount = EdgeObjectCount();
//     vec<int64_t> support_raw(nEdgeObjectCount, 0); 
     //vec<fix64_6> support( EdgeObjectCount( ), 0 );

     // Go through the reads.

     //this flattens the bubbles.size()-by-4 matrix row-by-row in to a 1D array
     const size_t uVotesBubWidth=4;
     vec<int64_t> votes_bub_raw(bubbles.size()*uVotesBubWidth,0); 
     //vec< vec<fix64_6> > votes_bub( bubbles.size( ), vec<fix64_6>(4,0) );

     //here, the fix64_6 encapsulation is broken to conform to openmp's 
     // atomic type requirement,
     //the following two scales are used for explicity manipulation of int64_t
     //the encapsulation can be recovered by changing fix64_6 appropiately

     const int64_t iFixedScale=1000000;
     const double  dFixedScale=1000000.0;

     if ( verbosity >= 1 ) cout << Date( ) << ": going through " << bases.size() << " reads" << endl;
     #pragma omp parallel
     {
          vec<int64_t> votes_bub_raw_loc(votes_bub_raw.size(),0);
          read_place_buffer places,places_part;
          
          exact_match_spec_t exact_match_spec;
          for(size_t ee=0;ee<bubbles_lookup.size();++ee){
              if(bubbles_lookup[ee].size()>0){
                  exact_match_spec[ee]=hb_fw.K()-1;
              }
          }
          report_t report(bubbles_lookup,bubbles.size());
          
          #pragma omp for schedule(dynamic,1) nowait
          for (size_t id=0;id<bases.size();++id){
               if ( verbosity >= 1 && id % 500000 == 0 )
               {    
                    #pragma omp critical
                    {    cout << Date( ) << ": processing read " << id 
                              << endl;    }    }
               if ( bases[id].isize( ) < L ) continue;
               int n = KmerId( bases[id], L, 0 );
               const int infinity = 1000000000;
               int qual_sum = infinity;
               loc_FindPlacesInt(bases[id], quals[id], n, hb_fw, hb_rc, to_right_fw, to_right_rc, locs_fw, locs_rc, report, qual_sum, 3, 20*1000, exact_match_spec, places_part );
               if(report.logged()){
                   const double dScale=iFixedScale/report.denominator();
                   const vec<double>& numerators=report.numerators();
                   ForceAssert(numerators.size()==votes_bub_raw_loc.size());
                   ForceAssert(numerators.size()==bubbles.size()*4);
                   
                   for(const auto& bi:report.activeBubbles()){
                       votes_bub_raw_loc[bi*uVotesBubWidth]+=dScale*numerators[bi*uVotesBubWidth];
                       votes_bub_raw_loc[bi*uVotesBubWidth+1]+=dScale*numerators[bi*uVotesBubWidth+1];
                       votes_bub_raw_loc[bi*uVotesBubWidth+2]+=dScale*numerators[bi*uVotesBubWidth+2];
                       votes_bub_raw_loc[bi*uVotesBubWidth+3]+=dScale*numerators[bi*uVotesBubWidth+3];
                   }
               }
          }
          #pragma omp critical
          {
              for(size_t ii=0;ii<votes_bub_raw_loc.size();++ii){
                  votes_bub_raw[ii]+=votes_bub_raw_loc[ii];
              }
          }
          votes_bub_raw_loc.clear();
     }//omp
     if ( verbosity >= 3 ) cout << "\n";
     // Test bubbles.

     const double max_asym_rarity = 0.1;
     const int min_to_save = 9;
     const double large_ratio=4;
     if ( verbosity >= 1 ) cout << "\nbubbles:\n";
     for ( int i = 0; i < bubbles.isize( ); i++ )
     {    double f1 = votes_bub_raw[uVotesBubWidth*i+0]/dFixedScale, r1 = votes_bub_raw[uVotesBubWidth*i+1]/dFixedScale; //double f1 = votes_bub[i][0].ToDouble( ), r1 = votes_bub[i][1].ToDouble( );
          double f2 = votes_bub_raw[uVotesBubWidth*i+2]/dFixedScale, r2 = votes_bub_raw[uVotesBubWidth*i+3]/dFixedScale; //double f2 = votes_bub[i][2].ToDouble( ), r2 = votes_bub[i][3].ToDouble( );
          if (verbosity>=1) std::cout<< "raw " << setiosflags(ios::fixed) << setprecision(1) << f1 << " " <<r1 << " vs " << f2 << " " << r2 <<std::endl <<resetiosflags(ios::fixed);
          int e1 = 1, e2 = 2;
          if ( f2 + r2 > f1 + r1 || ( f2 + r2 == f1 + r1 && f2 > f1 ) )
          {    swap( e1, e2 );
               swap( f1, f2 );
               swap( r1, r2 );    }
          if ( f2 > r2 || ( f2 == r2 && f1 > r1 ) )
          {    swap( f1, r1 ); 
               swap( f2, r2 );    }
          long double p = Min( 0.5, f1/(f1+r1) ) / 2;
          int n = int(floor(f1+r1+f2+r2));
          double q = -1;
          int iDeletedEdge=-1;

          if ( (f2 < 1.0 || r2 < 1.0 ) && f1>=1.0 && r1>=1.0){
              if ( e1 == 1 ) iDeletedEdge = bubbles[i][2] ;
              else           iDeletedEdge = bubbles[i][1] ;
          }
          else if ( (f1 < 1.0 || r1 < 1.0 ) && f2>=1.0 && r2>=1.0){
              if ( e1 == 1 ) iDeletedEdge = bubbles[i][1] ;
              else           iDeletedEdge = bubbles[i][2] ;
          }
          else if ( n > 0 && n <= 10000 ) 
          {    q = BinomialSum( n, int(ceil(f2)), p );
               if ( (q < max_asym_rarity || f1+r1 > large_ratio*(f2+r2)) && f2 + r2 < min_to_save )
               {    if ( e1 == 1 ) iDeletedEdge = bubbles[i][2] ;
                    else           iDeletedEdge = bubbles[i][1] ;
               }    }
          
          if(iDeletedEdge>=0) to_delete.push_back(iDeletedEdge);

          if ( verbosity >= 1 )
          {    cout << "[" << i << "] " << bubbles[i][0] << " --> {" 
                    << bubbles[i][1] << "," << bubbles[i][2] << "} --> " 
                    << bubbles[i][3] << "; vote: " << setiosflags(ios::fixed) 
                    << setprecision(1) << f1 << "+" << r1 << " vs " << f2 << "+" 
                    << r2 << resetiosflags(ios::fixed);
               if ( q >= 0 ) cout << ", surprise(" << e2 << ") = " << q;
               if (iDeletedEdge>=0) cout << " deleted edge " << iDeletedEdge << std::endl;
               cout << endl;    }    }
     if ( verbosity >= 1 ) cout << "\n";

     // Clean up graph.  Note that inverse of edges have to be added in.  I'm not
     // sure why these weren't found in the first place.

     int ndels = to_delete.size( );
     for ( int i = 0; i < ndels; i++ )
          if ( Inv( to_delete[i] ) >= 0 ) to_delete.push_back( Inv( to_delete[i] ) );
     DeleteEdges(to_delete);
     DeleteUnusedPaths( );
     REPORT_TIME( clock1, "used divining single mutation bubbles" );
     RemoveDeadEdgeObjects( );
     double clock2 = WallClockTime( );
     RemoveEdgelessVertices( );
     RemoveUnneededVertices( );
     REPORT_TIME( clock2, "used divining single mutation bubbles tail" );
     RemoveDeadEdgeObjects( );
     FixWeights(logc);
     TestValid(logc);    
     DeleteReverseComplementComponents(logc);    }
