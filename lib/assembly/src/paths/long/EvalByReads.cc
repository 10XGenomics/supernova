///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS
#include "Basevector.h"
#include "CoreTools.h"
#include "PrintAlignment.h"
#include "Qualvector.h"
#include "paths/HyperBasevector.h"
#include "paths/LongReadTools.h"
#include "paths/long/CreateGenome.h"
#include "paths/long/EvalByReads.h"
#include <queue>
#include <omp.h>

namespace {
    struct CompareReadPlaceByQsum {
        bool operator() (const read_place& a, const read_place& b) 
        { return a.Qsum() < b.Qsum(); }
    };
}

void FindPlaces( const basevector& b, const qualvector& q, const int n, 
     const HyperBasevector& hb_fw, const HyperBasevector& hb_rc, 
     const vec<int>& to_right_fw, const vec<int>& to_right_rc, 
     const VecIntPairVec& locs_fw, const VecIntPairVec& locs_rc,
     vec<read_place>& places, int& qual_sum, const int min_qual, const double prox,
     const int max_diff )
{    
     const int infinity = 1000000000;
     int best_qsum = infinity;
     priority_queue<read_place, std::vector<read_place>, CompareReadPlaceByQsum> candidates;
     for ( int opass = 1; opass <= 2; opass++ )
     {    const HyperBasevector& hb = ( opass == 1 ? hb_fw : hb_rc );
          const VecIntPairVec& locs = ( opass == 1 ? locs_fw : locs_rc );
          const vec<int>& to_right = ( opass == 1 ? to_right_fw : to_right_rc );
          vec<read_place> places_part;
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
               p.ComputeQsum( b, q, hb, min_qual );
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
                    {    candidates.push(p);
                         if ( p.Qsum( ) <= best_qsum ) {
                              best_qsum = p.Qsum( );    
                              while (candidates.top().Qsum() > best_qsum + prox)
                                  candidates.pop();    }    }    }
               else
               {    int v = to_right[ p.E( p.N( ) - 1 ) ];
                    for ( int j = 0; j < hb.From(v).isize( ); j++ )
                    {    int e = hb.EdgeObjectIndexByIndexFrom( v, j );
                         read_place r(p);
                         r.AddEdge( e, b, q, hb, min_qual );
                         if ( r.Qsum( ) <= best_qsum + prox )
                              places_part.push_back(r);    }    }    }    }
     places.clear( );
     places.reserve(candidates.size());
     while ( candidates.size() > 0 ) 
     {   places.push_back(candidates.top());
         candidates.pop();   }
     qual_sum = best_qsum;
}

uint64_t SafeFindPlaces( const basevector& b, const qualvector& q, const int n,
     const HyperBasevector& hb_fw, const HyperBasevector& hb_rc,
     const vec<int>& to_right_fw, const vec<int>& to_right_rc,
     const VecIntPairVec& locs_fw, const VecIntPairVec& locs_rc,
     vec<read_place>& places, int& qual_sum,
     uint64_t maxNSteps,
     const int min_qual, const double prox, const int max_diff )
{
     uint64_t n_terminated=0;
     const int infinity = 1000000000;
     int best_qsum = infinity;
     typedef priority_queue<read_place, std::vector<read_place>, CompareReadPlaceByQsum> candidates_t;
     vec<candidates_t> v_candidates;
     for ( int opass = 1; n_terminated==0 && opass <= 2; opass++ )
     {    const HyperBasevector& hb = ( opass == 1 ? hb_fw : hb_rc );
          const VecIntPairVec& locs = ( opass == 1 ? locs_fw : locs_rc );
          const vec<int>& to_right = ( opass == 1 ? to_right_fw : to_right_rc );
          for ( unsigned j = 0; n_terminated==0 && j < locs[n].size( ); j++ )
          {
               // Test to see if we're starting in the last K-1 bases of the
               // edge.  Note that we could have built the locs so that this would
               // be already done.
               vec<read_place> places_part;
               v_candidates.push_back(candidates_t());

               if ( locs[n][j].second >= hb.EdgeLengthKmers( locs[n][j].first ) )
                    continue;

               vec<int> v(1);
               v[0] = locs[n][j].first;
               read_place p( v, locs[n][j].second, opass == 1, 0 );
               p.ComputeQsum( b, q, hb, min_qual );
               places_part.push_back(p);

               uint64_t count=0;
               while( places_part.nonempty( ) )
               {    read_place p = places_part.back( );
                    places_part.pop_back( );
                    if ( p.Qsum( ) > best_qsum + prox ) continue;
 
                    if( ++count > maxNSteps){
                          places_part.clear();
                          v_candidates.pop_back();
                          ++n_terminated;
                          break;
                    }
 
                    int ext1 = b.isize( ) - ( hb.EdgeLengthBases( p.E(0) ) - p.P( ) );
                    for ( int j = 1; j < p.N( ); j++ )
                         ext1 -= hb.EdgeLengthKmers( p.E(j) );
                    if ( ext1 <= 0 )
                    {    if ( p.Qsum( ) <= best_qsum + prox )
                         {    v_candidates.back().push(p);
                              if ( p.Qsum( ) <= best_qsum ) {
                                   best_qsum = p.Qsum( );
                                   while (v_candidates.back().top().Qsum() > best_qsum + prox)
                                       v_candidates.back().pop();    }    }    }
                    else
                    {    int v = to_right[ p.E( p.N( ) - 1 ) ];
                         for ( int j = 0; j < hb.From(v).isize( ); j++ )
                         {    int e = hb.EdgeObjectIndexByIndexFrom( v, j );
                              read_place r(p);
                              r.AddEdge( e, b, q, hb, min_qual );
                              if ( r.Qsum( ) <= best_qsum + prox )
                                   places_part.push_back(r);    }    }
                }
           }
     }
     places.clear( );
     if( n_terminated >0 ){
//         #pragma omp critical
//         {
//             std::cout.precision(3);
//             std::cout << "SafeFindPlaces for a read has been terminated. The allowable maximum number of steps was set to "
//                       << maxNSteps
//                       << std::endl;
//         }
         qual_sum=infinity;
         return n_terminated;
     }
     size_t total_size=0;
     for(const auto& entry: v_candidates){ total_size+=entry.size(); }
     places.reserve(total_size);
     for(auto&entry:v_candidates){
         for( ; entry.size()>0 ; entry.pop()){
            if( entry.top().Qsum() <= best_qsum + prox ){
                places.push_back(entry.top());
            }
         }
     }
     qual_sum = best_qsum;
     return n_terminated;
}

void EvaluateRead( const int id, const HyperBasevector& hb_A_fw,
     const HyperBasevector& hb_A_rc, const HyperBasevector& hb_R_fw, 
     const HyperBasevector& hb_R_rc, const vec<int>& to_right_A_fw, 
     const vec<int>& to_right_A_rc, const vec<int>& to_right_R_fw,
     const vec<int>& to_right_R_rc, const int L, 
     const VecIntPairVec& Alocs_fw,
     const VecIntPairVec& Alocs_rc,
     const VecIntPairVec& Rlocs_fw,
     const VecIntPairVec& Rlocs_rc, const vecbasevector& bases,
     const vecqualvector& quals, vec<int>& qual_sum, int& max_perf, 
     const Bool verbose, vec< vec< triple<int,int,char> > >& edits )
{
     const basevector& b = bases[id];
     const qualvector& q = quals[id];
     const int infinity = 1000000000;
     qual_sum.resize_and_set( 2, infinity );
     max_perf = 0;
     const int min_qual = 3;
     if ( b.isize( ) < L ) return;
     int n = KmerId( b, L, 0 );
     for ( int pass = 0; pass < 2; pass++ )
     {    vec<read_place> places;

          FindPlaces( b, q, n, ( pass == 0 ? hb_A_fw : hb_R_fw ),
               ( pass == 0 ? hb_A_rc : hb_R_rc ),
               ( pass == 0 ? to_right_A_fw : to_right_R_fw ),
               ( pass == 0 ? to_right_A_rc : to_right_R_rc ),
               ( pass == 0 ? Alocs_fw : Rlocs_fw ),
               ( pass == 0 ? Alocs_rc : Rlocs_rc ), places, qual_sum[pass] );

          vec< triple<int,int,char> > edits0;
          for ( int l = 0; l < places.isize( ); l++ )
          {    const read_place& p = places[l];
               int ei = 0, pos = p.P( ), pc = 0, errs = 0;
               vec<char> ref;
               const HyperBasevector& hb = 
                    ( pass == 0 ? ( p.Fw( ) ? hb_A_fw : hb_A_rc ) 
                    : ( p.Fw( ) ? hb_R_fw : hb_R_rc ) );
               for ( int l = 0; l < b.isize( ); l++ )
               {    const basevector& edge = hb.EdgeObject( p.E(ei) );
                    ref.push_back( edge[pos] );
                    if ( b[l] != edge[pos] )
                    {    pc = 0;
                         if ( p.Fw( ) ) edits0.push( p.E(ei), pos, as_base( b[l] ) );
                         else 
                         {    edits0.push( p.E(ei), edge.isize( ) - pos - 1,
                                   as_base( 3 - b[l] ) );    }
                         errs++;    }
                    else pc++;
                    max_perf = Max( max_perf, pc );
                    pos++;
                    if ( pos == hb.EdgeObject( p.E(ei) ).isize( ) )
                    {    ei++;
                         pos = hb.K( ) - 1;    }    }

               if (verbose)
               {    basevector r( ref.size( ) );
                    for ( int i = 0; i < ref.isize( ); i++ )
                         r.Set( i, ref[i] );
                    #pragma omp critical
                    {    cout << "read " << id << " aligns to " 
                              << ( pass == 0 ? "A" : "R" )
                              << " with " << errs << " errors" 
                              << " and qual sum = " << qual_sum[pass] 
                              << ", start = " << p.E(0) << ":" << p.P( ) << endl;
                         cout << "\nalignment of read " << id << endl;
                         align a;
                         a.SetNblocks(1);
                         a.SetGap( 0, 0 );
                         a.SetLength( 0, b.size( ) );
                         a.Setpos1(0), a.Setpos2(0);
                         if ( p.Fw( ) )
                              PrintVisualAlignment( True, cout, b, r, a, q );
                         else
                         {    a.ReverseThis( b.size( ), r.size( ) );
                              basevector brc(b), rrc(r);
                              brc.ReverseComplement( );
                              rrc.ReverseComplement( );
                              qualvector qrc(q);
                              qrc.ReverseMe( );
                              PrintVisualAlignment( True, cout, 
                                   brc, rrc, a, qrc );    }    }    }    }

          UniqueSort(edits0);
          edits[pass].append(edits0);    }    }

void EvalByReads( const HyperBasevector& hb_A, const HyperBasevector& hb_R,
     const vecbasevector& bases, vecqualvector quals, int& assembly_count,
     int& reference_count, const Bool print_a, const Bool print_r )
{
     // Cap quality scores.

     const int cap_radius = 2;
     const int cap_count = 2;
     for ( int r = 0; r < cap_count; r++ )
     {    for ( int64_t id = 0; id < (int64_t) quals.size( ); id++ )
          {    vec<int> q( quals[id].size( ), 1000000 );
               for ( int j = 0; j < (int) quals[id].size( ); j++ )
               {    int start = Max( 0, j - cap_radius );
                    int stop = Min( (int ) quals[id].size( ) - 1, j + cap_radius );
                    for ( int l = start; l <= stop; l++ )
                         q[j] = Min( q[j], (int) quals[id][l] );    }
               for ( int j = 0; j < (int) quals[id].size( ); j++ )
                    quals[id][j] = q[j];    }    }

     // Set up HyperBasevectors.

     HyperBasevector hb_A_fw(hb_A), hb_A_rc(hb_A), hb_R_fw(hb_R), hb_R_rc(hb_R);
     hb_A_rc.Reverse( ), hb_R_rc.Reverse( );
     vec<int> to_right_A_fw, to_right_R_fw;
     vec<int> to_right_A_rc, to_right_R_rc;
     hb_A_fw.ToRight(to_right_A_fw), hb_R_fw.ToRight(to_right_R_fw);
     hb_A_rc.ToRight(to_right_A_rc), hb_R_rc.ToRight(to_right_R_rc);

     // Create indices to facilitate alignment.

     VecIntPairVec Alocs_fw, Alocs_rc, Rlocs_fw, Rlocs_rc;
     vecbasevector Ax_fw, Ax_rc, Rx_fw, Rx_rc;
     for ( int i = 0; i < hb_A_fw.EdgeObjectCount( ); i++ )
          Ax_fw.push_back( hb_A_fw.EdgeObject(i) );
     for ( int i = 0; i < hb_A_rc.EdgeObjectCount( ); i++ )
          Ax_rc.push_back( hb_A_rc.EdgeObject(i) );
     for ( int i = 0; i < hb_R_fw.EdgeObjectCount( ); i++ )
          Rx_fw.push_back( hb_R_fw.EdgeObject(i) );
     for ( int i = 0; i < hb_R_rc.EdgeObjectCount( ); i++ )
          Rx_rc.push_back( hb_R_rc.EdgeObject(i) );
     const int L = 12;
     CreateGlocs( Ax_fw, L, Alocs_fw );
     CreateGlocs( Ax_rc, L, Alocs_rc );
     CreateGlocs( Rx_fw, L, Rlocs_fw );
     CreateGlocs( Rx_rc, L, Rlocs_rc );

     // Align the reads.

     const int max_min_score = 100;
     const int min_perfect_match = 60;  // minimum perfect match for read alignments
     const int min_delta = 30;
     vec<int> favors_assembly, favors_reference;
     vec< vec< triple<int,int,char> > > all_edits(2);
     #pragma omp parallel for
     for ( int64_t id = 0; id < (int64_t) bases.size( ); id++ )
     {    vec<int> qual_sum;
          int max_perf;
          vec< vec< triple<int,int,char> > > edits(2);
          EvaluateRead( id, hb_A_fw, hb_A_rc, hb_R_fw, hb_R_rc, to_right_A_fw,
               to_right_A_rc, to_right_R_fw, to_right_R_rc, L, Alocs_fw, Alocs_rc,
               Rlocs_fw, Rlocs_rc, bases, quals, qual_sum, max_perf, False, edits );
          if ( Abs( qual_sum[0] - qual_sum[1] ) >= min_delta * 1000
               && max_perf >= min_perfect_match
               && Min( qual_sum[0], qual_sum[1] ) <= max_min_score * 1000 )
          {    
               #pragma omp critical
               {    if ( qual_sum[0] < qual_sum[1] ) 
                    {    favors_assembly.push_back(id);
                         all_edits[0].append( edits[1] );    }
                    else 
                    {    favors_reference.push_back(id);
                         all_edits[1].append( edits[0] );    }    }    }    }
     vec< vec< triple<int,int,char> > > all_editsx(2);
     for ( int m = 0; m < 2; m++ )
     {    Sort( all_edits[m] );
          for ( int i = 0; i < all_edits[m].isize( ); i++ )
          {    int j = all_edits[m].NextDiff(i);
               if ( j - i >= 2 ) all_editsx[m].push_back( all_edits[m][i] );
               i = j - 1;    }    }

     // Find signal reads having confirmed edits.

     vec<int> favors_assemblyx, favors_referencex;
     for ( int pass = 0; pass < 2; pass++ )
     {    const vec<int>& fav = ( pass == 0 ? favors_assembly : favors_reference );
          for ( int i = 0; i < fav.isize( ); i++ )
          {    vec<int> qual_sum;
               int max_perf;
               int id = fav[i];
               vec< vec< triple<int,int,char> > > edits(2);
               EvaluateRead( id, hb_A_fw, hb_A_rc, hb_R_fw, hb_R_rc, to_right_A_fw,
                    to_right_A_rc, to_right_R_fw, to_right_R_rc, L, Alocs_fw, 
                    Alocs_rc, Rlocs_fw, Rlocs_rc, bases, quals, qual_sum, max_perf, 
                    False, edits );
               Bool confirmed = False;
               for ( int j = 0; j < edits[1-pass].isize( ); j++ )
               {    if ( BinMember( all_editsx[pass], edits[1-pass][j] ) )
                         confirmed = True;    }
               if (confirmed)
               {    ( pass == 0 ? favors_assemblyx : favors_referencex )
                         .push_back(id);    }    }    }

     // Print alignments of signal reads.

     if (print_r) 
     {    cout << "\nreads favoring reference:\n";
          for ( int i = 0; i < favors_referencex.isize( ); i++ )
          {    vec<int> qual_sum;
               int max_perf;
               int id = favors_referencex[i];
               vec< vec< triple<int,int,char> > > edits(2);
               EvaluateRead( id, hb_A_fw, hb_A_rc, hb_R_fw, hb_R_rc, to_right_A_fw,
                    to_right_A_rc, to_right_R_fw, to_right_R_rc, L, Alocs_fw, 
                    Alocs_rc, Rlocs_fw, Rlocs_rc, bases, quals, qual_sum, max_perf, 
                    True, edits );
               PRINT4( id, qual_sum[0]/1000.0, qual_sum[1]/1000.0, max_perf );
          cout << "\n";     }    }
     if (print_a)
     {    cout << "\nreads favoring assembly:\n";
          {    for ( int i = 0; i < favors_assemblyx.isize( ); i++ )
               {    vec<int> qual_sum;
                    int max_perf;
                    int id = favors_assemblyx[i];
                    vec< vec< triple<int,int,char> > > edits(2);
                    EvaluateRead( id, hb_A_fw, hb_A_rc, hb_R_fw, hb_R_rc, 
                         to_right_A_fw, to_right_A_rc, to_right_R_fw, to_right_R_rc, 
                         L, Alocs_fw, Alocs_rc, Rlocs_fw, Rlocs_rc, bases, quals, 
                         qual_sum, max_perf, True, edits );
                    PRINT4( id, qual_sum[0]/1000.0, qual_sum[1]/1000.0, max_perf );
                    cout << "\n";     }    }    }

     // Report results.
     
     assembly_count = favors_assemblyx.size( );
     reference_count = favors_referencex.size( );    }

void ExtendPlacement( const read_place& p, const basevector& b, const qualvector& q,
     const HyperBasevector& hb_fw, const HyperBasevector& hb_rc, 
     const vec<int>& to_right_fw, const vec<int>& to_right_rc, 
     vec<read_place>& places, int& qual_sum, const int min_qual )
{    
     const int infinity = 1000000000;
     int best_qsum = infinity;
     places.clear( );
     int opass = ( p.Fw( ) ? 1 : 2 );
     const HyperBasevector& hb = ( opass == 1 ? hb_fw : hb_rc );
     const vec<int>& to_right = ( opass == 1 ? to_right_fw : to_right_rc );
     vec<read_place> places_part = {p}, places_done;
     while( places_part.nonempty( ) )
     {    read_place p = places_part.back( );
          places_part.pop_back( );
          if ( p.Qsum( ) > best_qsum ) continue;
          int ext1 = b.isize( ) - ( hb.EdgeLengthBases( p.E(0) ) - p.P( ) );
          for ( int j = 1; j < p.N( ); j++ )
               ext1 -= hb.EdgeLengthKmers( p.E(j) );
          if ( ext1 <= 0 ) 
          {    if ( p.Qsum( ) <= best_qsum )
               {    places_done.push_back(p);
                    best_qsum = p.Qsum( );    }    }
          else
          {    int v = to_right[ p.E( p.N( ) - 1 ) ];
               for ( int j = 0; j < hb.From(v).isize( ); j++ )
               {    int e = hb.EdgeObjectIndexByIndexFrom( v, j );
                    read_place r(p);
                    r.AddEdge( e, b, q, hb, min_qual );
                    if ( r.Qsum( ) <= best_qsum )
                         places_part.push_back(r);    }    }    }
     for ( int j = 0; j < places_done.isize( ); j++ )
     {    const read_place& p = places_done[j];
          int qsum = p.Qsum( );
          if ( qsum == qual_sum ) places.push_back(p);
          if ( qsum < qual_sum )
          {    places.clear( );
               places.push_back(p);
               qual_sum = qsum;    }    }    }

void read_place::writeBinary( BinaryWriter& writer ) const
{    writer.write( e_ );
     writer.write( p_ );
     writer.write( fw_ );
     writer.write( qsum_ );    }

void read_place::readBinary( BinaryReader& reader )
{    reader.read( &e_ );
     reader.read( &p_ );
     reader.read( &fw_ );
     reader.read( &qsum_ );    }
