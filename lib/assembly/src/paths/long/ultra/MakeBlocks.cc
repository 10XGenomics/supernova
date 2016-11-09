///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "CoreTools.h"
#include "Equiv.h"
#include "ParallelVecUtilities.h"
#include "PrintAlignment.h"
#include "graph/Digraph.h"
#include "math/HoInterval.h"
#include "pairwise_aligners/SmithWatBandedA.h"
#include "pairwise_aligners/SmithWatFree.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/ultra/MakeBlocks.h"
#include "paths/long/ultra/ThreadedBlocks.h"
#include "paths/LongReadTools.h"

template<int K> void MakeBlocks( 
     const int ec,                            // read id
     const vecbasevector& reads,
     const vec< vec< pair<int,int> > >& a,    // kmer aligns to read 0
     threaded_blocks& tb,                     // output
     const basevector& gread0,                // true genomic sequence for read 0
                                              // and flanking region
     const vec<basevector>& gkmers,
     const int gid,
     const int gstart,
     ostream& out,                            // for logging
     const Bool FILTER_BAD_BLOCKS,
     const long_logging_control& log_control,
     const long_logging& logc )
{
     // Heuristics.

     const int max_shift = 5;
     const int min_orbit = 3;
     const double win_ratio = 4.0;
     const int min_block = 50;

     // Test input requirements.

     for ( int id = 1; id < (int) reads.size( ); id++ )
     {    ForceAssert( a[id].nonempty( ) );
          for ( int j = 0; j < a[id].isize( ) - 1; j++ )
          {    ForceAssertLt( a[id][j].first, a[id][j+1].first );
               ForceAssertLt( a[id][j].second, a[id][j+1].second );    }    }

     // Define pointed kmers.

     int N = reads.size( );
     double point_clock = WallClockTime( );
     vec< pair<kmer<K>,int> > pkmers0;
     for ( int p = 0; p <= reads[0].isize( ) - K; p++ )
     {    kmer<K> x;
          x.SetToSubOf( reads[0], p );
          pkmers0.push( x, p );    }
     for ( int id = 1; id < N; id++ )
     {    for ( int p2 = 0; p2 < a[id][0].first; p2++ )
          {    int p1 = a[id][0].second - a[id][0].first + p2;
               if ( p1 < 0 || p1 >= reads[0].isize( ) ) continue;
               kmer<K> x;
               x.SetToSubOf( reads[id], p2 );
               pkmers0.push( x, p1 );    }
          for ( int j = 0; j < a[id].isize( ); j++ )
          {    int start2 = a[id][j].first, stop2;
               if ( j < a[id].isize( ) - 1 ) stop2 = a[id][j+1].first;
               else stop2 = reads[id].size( );
               for ( int p2 = start2; p2 < stop2; p2++ )
               {    int p1 = a[id][j].second + p2 - start2;
                    if ( p1 < 0 || p1 >= reads[0].isize( ) ) continue;
                    if ( p2 > reads[id].isize( ) - K ) continue;
                    kmer<K> x;
                    x.SetToSubOf( reads[id], p2 );
                    pkmers0.push( x, p1 );    }    }    }
     REPORT_TIMEX( point_clock, "used making pointed kmers", out );
     double sort_clock = WallClockTime( );
     Sort(pkmers0);
     REPORT_TIMEX( sort_clock, "used in sorting", out );

     // Define an equivalence relation on the pointed kmers, by joining (x,p) to
     // (x,p') if |p-p'| <= 5.

     double eclock = WallClockTime( );
     equiv_rel e( pkmers0.size( ) );
     for ( int j = 1; j < pkmers0.isize( ); j++ )
     {    if ( pkmers0[j].first == pkmers0[j-1].first
               && pkmers0[j].second - pkmers0[j-1].second <= max_shift )
          {    e.Join( j-1, j );    }    }
     REPORT_TIMEX( eclock, "used defining equivalence relation", out );

     // Delete orbits of size < 3.

     double dclock = WallClockTime( );
     vec<int> reps;
     e.OrbitReps(reps);
     vec<Bool> small_orbit( reps.size( ), False );
     for ( int j = 0; j < reps.isize( ); j++ )
     {    vec<int> o;
          e.Orbit( reps[j], o );
          if ( o.isize( ) < min_orbit ) small_orbit[j] = True;    }
     EraseIf( reps, small_orbit );
     int nreps = reps.size( );
     REPORT_TIMEX( dclock, "deleting small orbits", out );

     // Find connections, make a graph.

     double c1clock = WallClockTime( );
     vec< pair<kmer<K>,int> > qkmers0;
     vec<int> rep_id;
     for ( int j = 0; j < reps.isize( ); j++ )
     {    vec<int> o;
          e.Orbit( reps[j], o );
          for ( int l = 0; l < o.isize( ); l++ )
          {    qkmers0.push_back( pkmers0[ o[l] ] );
               rep_id.push_back(j);    }    }
     SortSync( qkmers0, rep_id );
     REPORT_TIMEX( c1clock, "used finding connections 1", out );
     double c2clock = WallClockTime( );
     vec< vec<int> > from(nreps), to(nreps);
     for ( int id = 1; id < N; id++ )
     {    vec< pair<kmer<K>,int> > tkmers0;
          tkmers0.reserve( reads[id].isize( ) - K + 1 );
          kmer<K> x;
          for ( int j = 0; j < a[id].isize( ); j++ )
          {    int start2 = a[id][j].first, stop2;
               if ( j < a[id].isize( ) - 1 ) stop2 = a[id][j+1].first;
               else stop2 = reads[id].size( );
               for ( int p2 = start2; p2 < stop2; p2++ )
               {    int p1 = a[id][j].second + p2 - start2;
                    if ( p2 > reads[id].isize( ) - K ) continue;
                    x.SetToSubOf( reads[id], p2 );
                    tkmers0.push( x, p1 );    }    }
          for ( int j = 0; j < tkmers0.isize( ) - 1; j++ )
          {    int x1 = BinPosition( qkmers0, tkmers0[j] );
               if ( x1 < 0 ) continue;
               int x2 = BinPosition( qkmers0, tkmers0[j+1] );
               if ( x2 < 0 ) continue;
               int v1 = rep_id[x1], v2 = rep_id[x2];
               from[v1].push_back(v2), to[v2].push_back(v1);    }    }
     for ( int v = 0; v < nreps; v++ )
     {    UniqueSort(from[v]), UniqueSort(to[v]);    }
     digraph G( from, to );
     REPORT_TIMEX( c2clock, "used finding connections 2", out );

     // Project onto read 0 and remove dominated kmers.

     double proj_clock = WallClockTime( );
     vec<Bool> bad( G.N( ), False );
     vec<int> start, stop, osize;
     for ( int v = 0; v < nreps; v++ )
     {    vec<int> o, pos;
          e.Orbit( reps[v], o );
          osize.push_back( o.size( ) );
          for ( int l = 0; l < o.isize( ); l++ )
               pos.push_back( pkmers0[ o[l] ].second );
          Sort(pos);
          start.push_back( pos.front( ) );
          stop.push_back( pos.back( ) );    }
     REPORT_TIMEX( proj_clock, "used in projecting", out );
     double bad_clock = WallClockTime( );
     vec<int> ids( nreps, vec<int>::IDENTITY ), max_stop(nreps);
     SortSync( start, stop, ids );
     int ms = 0;
     for ( int j = 0; j < nreps; j++ )
     {    max_stop[j] = ms;
          ms = Max( ms, stop[j] );    }
     for ( int j1 = 0; j1 < nreps; j1++ )
     {    int v1 = ids[j1];
          for ( int j2 = j1 + 1; j2 < nreps; j2++ )
          {    if ( start[j2] > stop[j1] + max_shift ) break;
               int v2 = ids[j2];
               if ( osize[v2] > win_ratio * osize[v1] ) bad[v1] = True;    }
          for ( int j2 = j1 - 1; j2 >= 0; j2-- )
          {    if ( stop[j2] + max_shift < start[j1] ) continue;
               int v2 = ids[j2];
               if ( osize[v2] > win_ratio * osize[v1] ) bad[v1] = True;
               if ( max_stop[j2] + max_shift < start[j1] ) break;    }    }
     REPORT_TIMEX( bad_clock, "used finding bads", out );

     // Define maximal unbranched stretches in the graph.

     double line_clock = WallClockTime( );
     vec< vec<int> > lines;
     vec<Bool> used( nreps, False );
     for ( int v = 0; v < nreps; v++ )
     {    if ( bad[v] ) continue;
          if ( used[v] ) continue;
          vec<int> line;
          line.push_back(v);
          used[v] = True;
          int w = v;
          while(1)
          {    vec<int> tos;
               for ( int j = 0; j < G.To(w).isize( ); j++ )
               {    int x = G.To(w)[j];
                    if ( !bad[x] ) tos.push_back(x);    }
               if ( !tos.solo( ) ) break;
               w = tos[0];
               if ( Member( line, w ) ) break;
               if ( used[w] ) break; // ???
               int froms = 0;
               for ( int j = 0; j < G.From(w).isize( ); j++ )
               {    int x = G.From(w)[j];
                    if ( !bad[x] ) froms++;    }
               if ( froms != 1 ) break;
               line.push_front(w);
               used[w] = True;    }
          while(1)
          {    vec<int> froms;
               for ( int j = 0; j < G.From(w).isize( ); j++ )
               {    int x = G.From(w)[j];
                    if ( !bad[x] ) froms.push_back(x);    }
               if ( !froms.solo( ) ) break;
               w = froms[0];
               int tos = 0;
               for ( int j = 0; j < G.To(w).isize( ); j++ )
               {    int x = G.To(w)[j];
                    if ( !bad[x] ) tos++;    }
               if ( tos != 1 ) break;
               if ( Member( line, w ) ) break;
               if ( used[w] ) break; // ???
               line.push_back(w);
               used[w] = True;    }
          lines.push_back(line);    }
     REPORT_TIMEX( line_clock, "used finding lines", out );

     // Compute and order blocks.

     double comp_clock = WallClockTime( );
     vec< vec<int> > blocksx;
     vec<int> line_id;
     for ( int m = 0; m < lines.isize( ); m++ )
     {    const vec<int>& L = lines[m];
          if ( L.isize( ) < min_block ) continue;
          blocksx.push_back(L);
          line_id.push_back(m);    }
     vec<int> bpos;
     for ( int m = 0; m < blocksx.isize( ); m++ )
     {    const vec<int>& L = blocksx[m];
          vec<int> pos;
          for ( int k = 0; k < L.isize( ); k++ )
          {    int v = L[k];
               vec<int> o;
               e.Orbit( reps[v], o );
               for ( int l = 0; l < o.isize( ); l++ )
                    pos.push_back( pkmers0[ o[l] ].second );    }
          Sort(pos);
          bpos.push_back( pos[ pos.size( )/2 ] );    }
     SortSync( bpos, blocksx, line_id );

     // Translate pkmers0 into pkmers.

     vec< pair<basevector,int> > pkmers( pkmers0.size( ) );
     for ( int j = 0; j < pkmers0.isize( ); j++ )
     {    pkmers0[j].first.GetBasevector( pkmers[j].first );
          pkmers[j].second = pkmers0[j].second;    }
     REPORT_TIMEX( comp_clock, "used computing blocks", out );

     // Two successive blocks could overlap by up to K-1.  Trim their ends so that 
     // they don't overlap and in fact so that there is a small gap between them.
     // The small gap is desirable between we don't want threads to have negative
     // lengths.  Of course when could further trim the ends as needed at the point
     // when we create the threads.

     double trim_clock = WallClockTime( );
     const int target_gap = 3;
     for ( int j = 0; j < blocksx.isize( ) - 1; j++ )
     {    vec<int> &L1 = blocksx[j], &L2 = blocksx[j+1];
          int over;
          const basevector& b1 = pkmers[ reps[ L1.back( ) ] ].first;
          const basevector& b2 = pkmers[ reps[ L2.front( ) ] ].first;
          for ( over = K-1; over >= 1; over-- )
               if ( b1.Overlap( b2, over ) ) break;
          int trim = over + target_gap;
          for ( int j = 0; j < trim; j++ )
          {    if ( L1.size( ) >= L2.size( ) ) L1.resize( L1.isize( ) - 1 );
               else L2.erase( L2.begin( ) );    }    }
     REPORT_TIMEX( trim_clock, "used trimming", out );

     // Compute block positions again.

     double block_clock = WallClockTime( );
     vec<ho_interval> bposh;
     for ( int m = 0; m < blocksx.isize( ); m++ )
     {    const vec<int>& L = blocksx[m];
          vec<int> pos;
          for ( int k = 0; k < L.isize( ); k++ )
          {    int v = L[k];
               vec<int> o;
               e.Orbit( reps[v], o );
               for ( int l = 0; l < o.isize( ); l++ )
                    pos.push_back( pkmers[ o[l] ].second );    }
          Sort(pos);
          bposh.push( pos.front( ), pos.back( ) + K );    }

     // Compute the sequences associated to the blocks.

     vec<basevector> blocks;
     for ( int m = 0; m < blocksx.isize( ); m++ )
     {    const vec<int>& L = blocksx[m];
          basevector seq;
          for ( int k = 0; k < L.isize( ); k++ )
          {    int v = L[k];
               vec<int> o;
               e.Orbit( reps[v], o );
               const basevector& b = pkmers[ o[0] ].first;
               if ( seq.size( ) == 0 ) seq = b;
               else 
               {    seq.resize( seq.isize( ) - (K-1) );
                    seq = Cat( seq, b );    }    }
          blocks.push_back(seq);    }
     REPORT_TIMEX( block_clock, "used block seqs", out );

     // Predict initial positions of blocks on each read.

     double pred_clock = WallClockTime( );
     vec< vec<Bool> > block_defined(N);
     vec< vec<int> > block_start(N), block_stop(N);
     vec< vec<double> > block_err_rate(N);
     for ( int id = 0; id < N; id++ )
     {
          // Create map of positions on read 0 to positions on read id.

          vec<int> to_id( reads[0].size( ) + 1, -1 );
          for ( int p = 0; p < a[id][0].second; p++ )
               to_id[p] = p + a[id][0].first - a[id][0].second;
          for ( int j = 0; j < a[id].isize( ) - 1; j++ )
          {    for ( int p = a[id][j].second; p < a[id][j+1].second; p++ )
                    to_id[p] = p + a[id][j].first - a[id][j].second;    }
          for ( int p = a[id].back( ).second; p <= reads[0].isize( ); p++ )
               to_id[p] = p + a[id].back( ).first - a[id].back( ).second;

          // Now map the blocks.

          block_defined[id].resize( blocks.size( ), False );
          block_start[id].resize( blocks.size( ) );
          block_stop[id].resize( blocks.size( ) );
          block_err_rate[id].resize( blocks.size( ) );
          for ( int m = 0; m < blocksx.isize( ); m++ )
          {    const vec<int>& L = blocksx[m];
               int left = 0, right = 0;
               if ( bposh[m].Start( ) < 0 ) left = -bposh[m].Start( );
               if ( bposh[m].Stop( ) > reads[0].isize( ) )
                    right = bposh[m].Stop( ) - reads[0].isize( );
               int start = left + to_id[ Max( 0, bposh[m].Start( ) ) ];
               int stop = right 
                    + to_id[ Min( reads[0].isize( ), bposh[m].Stop( ) ) ] + K - 1;
               if ( logc.verb[ "ULTRA" ] >= 3 )
               {    out << "\nread " << id << ", block " << m << " --> "
                         << start << "-" << stop << "\n";    }

               // Align block.  For now assume fully on read id.

               if ( start < 0 || stop > reads[id].isize( ) ) continue;
               if ( !( start < stop ) ) continue;
               align a;
               const int add = 5;
               start = start - add;
               stop = stop + add;
               if ( start < 0 ) start = 0;
               if ( stop > reads[id].isize( ) ) stop = reads[id].size( );
               basevector r( reads[id], start, stop - start );

               int errors;
               const int bw_add = 10 + blocks[m].isize( )/50;
               int offset = -( r.isize( ) - blocks[m].isize( ) ) / 2;
               int bandwidth = Max( 0, ( r.isize( ) - blocks[m].isize( ) ) / 2 );
               int errs = SmithWatBandedA( blocks[m], r, offset, bandwidth + bw_add,
                    a, errors, 0, 1, 1 );

               // SmithWatFreeSym( blocks[m], r, a, false, false, 1, 1 );

               int start2 = start + a.pos2( ), stop2 = start + a.Pos2( );
               block_err_rate[id][m] = double(errs) / double( stop2 - start2 );
               if ( logc.verb[ "ULTRA" ] >= 3 )
               {    PRINT2_TO( out, start2, stop2 );
                    PrintVisualAlignment( True, out, blocks[m], 
                         r, a );    }    
               block_defined[id][m] = True;
               block_start[id][m] = start2, block_stop[id][m] = stop2;    }    }
     REPORT_TIMEX( pred_clock, "used predicting positions", out );
     vec< vec<int> > block_start_orig(block_start), block_stop_orig(block_stop);

     // Delete blocks that appear to be off target.

     double fclock = WallClockTime( );
     redelete:
     vec<Bool> blocks_to_delete( blocks.size( ), False );
     for ( int m = 0; m < blocks.isize( ); m++ )
     {    const double e_max = 0.1;
          if ( block_defined[0][m] && block_err_rate[0][m] >= e_max )
               blocks_to_delete[m] = True;    }
     for ( int m = 0; m < blocks.isize( ) - 1; m++ )
     {    if ( block_defined[0][m] && block_defined[0][m+1] )
          {    int over = block_stop[0][m] - block_start[0][m+1];
               if ( over <= K ) continue; // note overlaps < K are ignored

               // Does it appear that one block has a much higher error rate
               // than the other block?

               double e1 = block_err_rate[0][m], e2 = block_err_rate[0][m+1];
               const double e_min = 0.02;
               const double e_mult = 4.0;
               if ( e2 >= e_mult * Max( e1, e_min ) ) 
               {    blocks_to_delete[m+1] = True;
                    if ( logc.verb[ "ULTRA" ] >= 1 )
                    {    out << "deleting block for line " << line_id[m+1] 
                              << "\n";    }    }
               if ( e1 >= e_mult * Max( e2, e_min ) ) 
               {    blocks_to_delete[m] = True;    
                    if ( logc.verb[ "ULTRA" ] >= 1 )
                    {    out << "deleting block for line " << line_id[m] 
                              << "\n";    }    }    }    }
     if ( Sum(blocks_to_delete) > 0 )
     {    for ( int id = 0; id < N; id++ )
          {    EraseIf( block_start[id], blocks_to_delete );
               EraseIf( block_stop[id], blocks_to_delete );
               EraseIf( block_start_orig[id], blocks_to_delete );
               EraseIf( block_stop_orig[id], blocks_to_delete );
               EraseIf( block_err_rate[id], blocks_to_delete );
               EraseIf( block_defined[id], blocks_to_delete );    }
          EraseIf( blocks, blocks_to_delete );
          EraseIf( blocksx, blocks_to_delete );
          EraseIf( bpos, blocks_to_delete );
          EraseIf( bposh, blocks_to_delete );
          goto redelete;    }
     REPORT_TIMEX( fclock, "used off target", out );

     // Print block error rates.
     
     if (logc.PRINT_BLOCK_ERROR_RATES)
     {    vec< vec<String> > rows;
          vec<double> err(N);
          for ( int id = 0; id < N; id++ )
          {    double total_bases = 0, total_errs = 0;
               for ( int m = 0; m < blocks.isize( ); m++ )
               {    if ( !block_defined[id][m] ) continue;
                    int n = block_stop[id][m] - block_start[id][m];
                    total_bases += n;
                    total_errs += block_err_rate[id][m] * double(n);    }
               err[id] = total_errs / total_bases;    }
          out << "\nblock error rates:\n";
          for ( int id = 0; id < N; id++ )
          {    vec<String> row;
               row.push_back( "[" + ToString(id) + "]" );
               ostringstream out1;
               out1 << setiosflags(ios::fixed) << setprecision(1) << 100.0 * err[id] 
                    << resetiosflags(ios::fixed) << "%;";
               row.push_back( out1.str( ) );
               for ( int m = 0; m < blocks.isize( ); m++ )
               {    ostringstream out2;
                    if ( block_defined[id][m] )
                    {    out2 << " " << setiosflags(ios::fixed) << setprecision(1)
                              << 100.0 * block_err_rate[id][m] << "%";    }
                    row.push_back( out2.str( ) );    }
               rows.push_back(row);    }
          PrintTabular( out, rows, 1, String( blocks.size( ) + 1, 'r' ) );    }

     // Kill reads that have too many errors.  Currently this is done using a
     // hard threshold, which can't be right.

     double bclock = WallClockTime( );
     if (FILTER_BAD_BLOCKS)
     {    const double max_block_err_rate = 0.06;
          for ( int id = 1; id < N; id++ )
          {    Bool bad = False;
               for ( int m = 0; m < blocks.isize( ); m++ )
                    if ( block_err_rate[id][m] > max_block_err_rate ) bad = True;
               if (bad)
               {    for ( int m = 0; m < blocks.isize( ); m++ )
                         block_defined[id][m] = False;    }    }    }
     REPORT_TIMEX( bclock, "used filtering bad blocks", out );

     // Prevent illegal block overlap.

     double over_clock = WallClockTime( );
     for ( int m = 0; m < blocks.isize( ) - 1; m++ )
     {    for ( int id = 0; id < (int) reads.size( ); id++ )
          {    if ( block_defined[id][m] && block_defined[id][m+1] )
               {    ForceAssertLe( block_start[id][m+1], reads[id].isize( ) ); // XXX
                    int over = block_stop[id][m] - block_start[id][m+1];
                    if ( over <= 0 ) continue;
                         
                    // Deal with weird failures.

                    // ForceAssertLt( over, blocksx[m].isize( ) );
                    if ( !( over < blocksx[m].isize( ) ) )
                    {    if ( logc.verb[ "ULTRA" ] >= 2 )
                         {    out << "\nfunny condition" << endl;
                              PRINT2_TO( out, over, blocksx[m].size( ) );    }
                         continue;    }

                    // Proceed with normal case.

                    blocks[m].resize( blocks[m].isize( ) - over );
                    blocksx[m].resize( blocksx[m].isize( ) - over );
                    for ( int id = 0; id < (int) reads.size( ); id++ )
                         block_stop[id][m] -= over;    }    }    }
     REPORT_TIMEX( over_clock, "used preventing overlap", out );

     // Predict the position of each block on each read.  Define threaded blocks.

     double thread_clock = WallClockTime( );
     vec< vec<basevector> > threads( reads.size( ) );
     vec< vec<Bool> > thread_defined( reads.size( ) );
     vec<Bool> alive( reads.size( ), False );
     vec<ho_interval> thread_range( reads.size( ) );
     if ( blocks.nonempty( ) )
     {    for ( int id = 0; id < (int) reads.size( ); id++ )
          {    threads[id].resize( blocks.size( ) - 1 );
               thread_defined[id].resize( blocks.size( ) - 1 , False );
               for ( int m = 0; m < blocks.isize( ) - 1; m++ )
               {    if ( block_defined[id][m] && block_defined[id][m+1] )
                    {    
                         // Work around weird failure.

                         if ( !( block_stop[id][m] <= block_start[id][m+1] ) )
                         {    if ( logc.verb[ "ULTRA" ] >= 2 )
                                   out << "working around weird failure" << endl;
                              break;    }

                         // Proceed with normal case.

                         alive[id] = True;
                         threads[id][m] = basevector( reads[id], block_stop[id][m],
                              block_start[id][m+1] - block_stop[id][m] );    
                         thread_defined[id][m] = True;    }    }
               if ( alive[id] )
               {    int thread_start = -1, thread_stop = -1;
                    for ( int j = 0; j < blocks.isize( ); j++ )
                    {    if ( thread_defined[id][j] )
                         {    thread_start = j;
                              break;    }    }
                    for ( int j = blocks.isize( ) - 2; j >= 0; j-- )
                    {    if ( thread_defined[id][j] )
                         {    thread_stop = j+1;
                              break;    }    }
                    thread_range[id] 
                         = ho_interval( thread_start, thread_stop );    }    }    }
     tb = threaded_blocks( blocks, threads, alive, thread_range );
     REPORT_TIMEX( thread_clock, "used making threads", out );

     // Print summary stats.

     double sum_clock = WallClockTime( );
     if ( logc.verb[ "ULTRA" ] >= 1 )
     {    out << "\nSummary stats for read " << ec << ":\n";
          out << "lines = " << lines.size( ) << "\n";
          vec<int> line_sizes;
          for ( int m = 0; m < lines.isize( ); m++ )
               line_sizes.push_back( lines[m].size( ) );
          ReverseSort(line_sizes);
          line_sizes.resize( Min( line_sizes.isize( ), 20 ) );
          out << "top 20 = " << Sum(line_sizes) 
               << " of " << reads[0].size( ) - K + 1 << "\n";
          if ( gkmers.nonempty( ) )
          {    int mixed = 0;
               for ( int m = 0; m < lines.isize( ); m++ )
               {    const vec<int>& L = lines[m];
                    int trues = 0, falses = 0;
                    for ( int k = 0; k < L.isize( ); k++ )
                    {    int v = L[k];
                         vec<int> o;
                         e.Orbit( reps[v], o );
                         const basevector& b = pkmers[ o[0] ].first;
                         if ( BinMember( gkmers, b ) ) trues++;
                         else falses++;    }
                    if ( trues > 0 && falses > 0 ) mixed++;    }
               PRINT_TO( out, mixed );    }    }
     if ( logc.verb[ "ULTRA" ] >= 1 && gread0.isize( ) > 0 )
     {    int false_blocks = 0;
          vec< vec<ho_interval> > block_pos;
          for ( int m = 0; m < blocksx.isize( ); m++ )
          {    const vec<int> L = blocksx[m];
               basevector seq;
               for ( int k = 0; k < L.isize( ); k++ )
               {    int v = L[k];
                    vec<int> o;
                    e.Orbit( reps[v], o );
                    const basevector& b = pkmers[ o[0] ].first;
                    if ( seq.size( ) == 0 ) seq = b;
                    else
                    {    seq.resize( seq.isize( ) - (K-1) );
                         seq = Cat( seq, b );    }    }
               vec<ho_interval> h;
               String gread0s = gread0.ToString( ), seqs = seq.ToString( );
               for ( int p = 0; p <= gread0.isize( ) - seq.isize( ); p++ )
                    if ( gread0s.Contains( seqs, p ) ) h.push( p, p + seq.isize( ) );
               String gread0s_rc;
               StringReverseComplement( gread0s, gread0s_rc );
               for ( int p = 0; p <= gread0.isize( ) - seq.isize( ); p++ )
               {    if ( gread0s_rc.Contains( seqs, p ) ) 
                         h.push( p, p + seq.isize( ) );    }
               if ( h.empty( ) ) false_blocks++;
               block_pos.push_back(h);    }
          out << "false blocks = " << false_blocks << "\n\n";
          vec< vec<String> > rows;
          vec<String> row1, row2;
          row1.push_back( "block", "line", "size" );
          row2.push_back( "",      " id ", "" );
          if ( gkmers.nonempty( ) )
          {    row1.push_back( "actual pos" );
               row2.push_back( "on genome" );    }
          row1.push_back( "pos on",  "overlaps", "error" );
          row2.push_back( "founder", "next",     "  %  " );
          rows.push_back( row1, row2 );
          for ( int m = 0; m < block_pos.isize( ); m++ )
          {    vec<String> row;
               row.push_back( ToString(m), ToString( line_id[m] ) );
               row.push_back( ToString( blocks[m].size( ) ) );
               ostringstream out2, out3, out4, out5;
               if ( gkmers.nonempty( ) )
               {    for ( int l = 0; l < block_pos[m].isize( ); l++ )
                    {    out2 << gid << "." << block_pos[m][l] + gstart;
                         if ( l < block_pos[m].isize( ) - 1 ) out2 << " ";    }
                    row.push_back( out2.str( ) );    }
               if ( block_defined[0][m] )
               {    out3 << block_start_orig[0][m] << "-" << block_stop_orig[0][m];
                    if ( m < block_pos.isize( ) - 1 && block_defined[0][m+1] )
                    {    int gap = block_start_orig[0][m+1] - block_stop_orig[0][m];
                         if ( gap < 0 ) out4 << "+" << -gap;    }
                    out5 << setiosflags(ios::fixed) << setprecision(1) << 100.0 
                         * block_err_rate[0][m] << resetiosflags(ios::fixed);    }
               else
               {    out3 << "---";
                    out5 << "---";    }
               row.push_back( out3.str( ), out4.str( ), out5.str( ) );
               rows.push_back(row);    }
          PrintTabular( out, rows, 2, "rrrlllr" );    }
     REPORT_TIMEX( sum_clock, "used printing summary stats", out );

     // Make blocks and print lines.

     double make_clock = WallClockTime( );
     for ( int m = 0; m < lines.isize( ); m++ )
     {    const vec<int> L = lines[m];
          int trues = 0, falses = 0;
          if ( logc.verb[ "ULTRA" ] >= 2 && gread0.isize( ) > 0 
               && gkmers.nonempty( ) )
          {    for ( int k = 0; k < L.isize( ); k++ )
               {    int v = L[k];
                    vec<int> o;
                    e.Orbit( reps[v], o );
                    const basevector& b = pkmers[ o[0] ].first;
                    if ( BinMember( gkmers, b ) ) trues++;
                    else falses++;    }    }
          if ( logc.verb[ "ULTRA" ] >= 2 )
          {    out << "\n" << m << ". line of length " << L.size( );
               if ( gread0.isize( ) > 0 && gkmers.nonempty( ) )
               {    if ( falses == 0 ) out << " (all true)";
                    else if ( trues == 0 ) out << " (all false)";
                    else out << " (MIXED)";    }
               out << "\n";    }
          basevector seq;
          for ( int k = 0; k < L.isize( ); k++ )
          {    int v = L[k];
               vec<int> o;
               e.Orbit( reps[v], o );
               const basevector& b = pkmers[ o[0] ].first;
               if ( seq.size( ) == 0 ) seq = b;
               else
               {    seq.resize( seq.isize( ) - (K-1) );
                    seq = Cat( seq, b );    }    }
          if ( logc.verb[ "ULTRA" ] >= 2 )
          {    for ( int k = 0; k < L.isize( ); k++ )
               {    if ( falses == 0 && L.isize( ) > 2 )
                    {    if ( k > 0 && k < L.isize( ) - 1 )
                         {    if ( k == 1 )
                              {    out << "...\n";
                                   seq.Print(out);
                                   out << "...\n";    }
                              continue;    }    }
                    int v = L[k];
                    vec<int> o, pos;
                    e.Orbit( reps[v], o );
                    const basevector& b = pkmers[ o[0] ].first;
                    for ( int l = 0; l < o.isize( ); l++ )
                         pos.push_back( pkmers[ o[l] ].second );
                    Sort(pos);
                    out << "[" << v << "] " << b.ToString( ) << " {" 
                         << pos.size( ) << ":";
                    for ( int i = 0; i < pos.isize( ); i++ )
                    {    int j = pos.NextDiff(i);
                         if ( i > 0 ) out << ",";
                         out << pos[i] << "[" << j-i << "]";
                         i = j - 1;    }
                    out << "}";
                    if ( gread0.isize( ) > 0 &&
                         gkmers.nonempty( ) && !BinMember( gkmers, b ) ) 
                    {    out << " (FALSE)";    }
                    if ( k == L.isize( ) - 1 )
                    {    int froms = 0;
                         for ( int j = 0; j < G.From(v).isize( ); j++ )
                         {    int w = G.From(v)[j];
                              if ( !bad[w] ) froms++;    }
                         out << "\n-->" << " <" << froms << ">";
                         for ( int j = 0; j < G.From(v).isize( ); j++ )
                         {    int w = G.From(v)[j];
                              if ( !bad[w] ) out << " " << w;    }    }
                    out << "\n";    }    }    }
     REPORT_TIMEX( make_clock, "used making blocks", out );
     double tclock = WallClockTime( );
     if ( logc.PRINT_THREADS || logc.PRINT_BLOCKS )
     {    vec< vec<basevector> > gap_truth;
          if ( log_control.G->size( ) > 0 ) 
          {    tb.GetGapTruth( *log_control.G, log_control.LG, *log_control.Glocs, 
                    gap_truth );    }
          for ( int g = 0; g < tb.NBlocks( ); g++ )
          {    if (logc.PRINT_BLOCKS) 
               {    out << "\n";
                    tb.Block(g).Print( out, "block_" + ToString(g+1) );    }
               if ( logc.PRINT_THREADS && g < tb.NBlocks( ) - 1 )
               {    tb.PrintThreadSummary( g, out, *log_control.G, log_control.LG, 
                         *log_control.Glocs, gap_truth );    }    }    }
     REPORT_TIMEX( tclock, "used summarizing blocks and threads", out );    }

template void MakeBlocks<16>( const int ec, const vecbasevector& reads,
     const vec< vec< pair<int,int> > >& a, threaded_blocks& tb,
     const basevector& gread0, const vec<basevector>& gkmers, const int gid, 
     const int gstart, ostream& out, const Bool FILTER_BAD_BLOCKS,
     const long_logging_control&, const long_logging& logc );

template void MakeBlocks<20>( const int ec, const vecbasevector& reads,
     const vec< vec< pair<int,int> > >& a, threaded_blocks& tb,
     const basevector& gread0, const vec<basevector>& gkmers, const int gid, 
     const int gstart, ostream& out, const Bool FILTER_BAD_BLOCKS,
     const long_logging_control&, const long_logging& logc );

template void MakeBlocks<24>( const int ec, const vecbasevector& reads,
     const vec< vec< pair<int,int> > >& a, threaded_blocks& tb,
     const basevector& gread0, const vec<basevector>& gkmers, const int gid, 
     const int gstart, ostream& out, const Bool FILTER_BAD_BLOCKS,
     const long_logging_control&, const long_logging& logc );

template void MakeBlocks<28>( const int ec, const vecbasevector& reads,
     const vec< vec< pair<int,int> > >& a, threaded_blocks& tb,
     const basevector& gread0, const vec<basevector>& gkmers, const int gid, 
     const int gstart, ostream& out, const Bool FILTER_BAD_BLOCKS,
     const long_logging_control&, const long_logging& logc );
