///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"
#include "MapReduceEngine.h"
#include "ParallelVecUtilities.h"
#include "Qualvector.h"
#include "Set.h"
#include "feudal/Algorithms.h"
#include "feudal/HashSet.h"
#include "feudal/Mempool.h"
#include "kmers/KMer.h"
#include "kmers/KmerRecord.h"
#include "kmers/MakeLookup.h"
#include "math/Functions.h"
#include "paths/HyperBasevector.h"
#include "paths/RemodelGapTools.h"
#include "paths/long/ReadPath.h"
#include "paths/long/large/GapToyTools.h"
#include "system/RunTime.h"
#include "system/SpinLockedData.h"
#include "system/WorklistN.h"
#include <algorithm>

// AlignToGenome.  Find alignments of assembly edges to the genome.  Currently this
// maps edges to pairs (g,p) consisting of a genome contig g and an inferred start
// position p of the edge on g.
//
// First we find 60-mer matches between the assembly and the genome.  We require
// that the K2-mer occur exactly once in the assembly and once in the genome.  
// Each such match is required to extend perfectly to a total length of 500 bases.
// (Doesn't seem like we're still doing this.)

void AlignToGenome( const HyperBasevector& hb, const vec<int>& inv,
     const vecbasevector& genome, vec< vec< pair<int,int> > >& hits,
     const int K2 )
{
     // Heuristics.

     const int K = 60;
     const int max_gap = 20000;
     const int max_gmult = 4;

     // Computational performance heuristics.

     const int64_t batches = 100;

     // Go through four passes.

     double clock = WallClockTime( );
     int nobj = hb.EdgeObjectCount( ), ngenome = genome.size( );
     hits.clear( );
     hits.resize(nobj);
     for ( int pi = 0; pi < 4; pi++ )
     {    
          // Build lookup table.

          double clock2 = WallClockTime( );
          vec< triple<kmer<K>,int,int> > kmers_plus;
          vec<int64_t> starts;
          starts.reserve( ngenome + nobj );
          starts.push_back(0);
          {    vec<int> counts( ngenome + nobj, 0 );
               #pragma omp parallel for
               for ( int i = 0; i < ngenome; i++ )
               {    const basevector& u = genome[i];
                    counts[i] = 0;
                    for ( int j = 0; j <= u.isize( ) - K2; j++ )
                         if ( u[j] == pi ) counts[i]++;    }
               #pragma omp parallel for
               for ( int i = 0; i < nobj; i++ )
               {    const basevector& u = hb.EdgeObject(i);
                    counts[ genome.size( ) + i ] = 0;
                    for ( int j = 0; j <= u.isize( ) - K2; j++ )
                         if ( u[j] == pi ) counts[ ngenome + i ]++;    }
               for ( int i = 0; i < ngenome + nobj; i++ )
                    starts.push_back( starts.back( ) + counts[i] );    }
          kmers_plus.resize( starts.back( ) );
          #pragma omp parallel for
          for ( int i = 0; i < ngenome; i++ )
          {    const basevector& u = genome[i];
               int count = 0;
               kmer<K> x;
               for ( int j = 0; j <= u.isize( ) - K2; j++ )
               {    if ( u[j] == pi )
                    {    int64_t r = starts[i] + count;
                         x.SetToSubOf( u, j ); 
                         kmers_plus[r].first = x;
                         kmers_plus[r].second = i, kmers_plus[r].third = j;    
                         count++;    }    }    }
          #pragma omp parallel for
          for ( int i = 0; i < nobj; i++ )
          {    const basevector& u = hb.EdgeObject(i);
               int count = 0;
               kmer<K> x;
               for ( int j = 0; j <= u.isize( ) - K2; j++ )
               {    if ( u[j] == pi )
                    {    int64_t r = starts[ngenome+i] + count;
                         x.SetToSubOf( u, j ); 
                         kmers_plus[r].first = x;
                         kmers_plus[r].second = ngenome + i, kmers_plus[r].third = j;
                         count++;    }    }    }
          // cout << TimeSince(clock2) << " used building kmers" << endl;
          double clock2b = WallClockTime( );
          ParallelSort(kmers_plus);
          // cout << TimeSince(clock2b) << " used sorting kmers" << endl;

          // Traverse the kmers.

          double clock3a = WallClockTime( );
          vec<int64_t> bstart(batches+1);
          for ( int64_t i = 0; i <= batches; i++ )
               bstart[i] = ( (int64_t) kmers_plus.size( ) * i ) / batches;
          #pragma omp parallel for schedule(dynamic, 1)
          for ( int64_t i = 1; i < batches; i++ )
          {    int64_t& s = bstart[i];
               while( s > 0 && kmers_plus[s].first == kmers_plus[s-1].first )
               {    s--;    }    }
          // cout << TimeSince(clock3a) << " used setting up traversal" << endl;
          double clock3 = WallClockTime( );
          #pragma omp parallel for schedule(dynamic, 1)
          for ( int64_t bi = 0; bi < batches; bi++ )
          for ( int64_t i = bstart[bi]; i < bstart[bi+1]; i++ )
          {    int64_t j, k;
               for ( j = i + 1; j < bstart[bi+1]; j++ )
                    if ( kmers_plus[j].first != kmers_plus[i].first ) break;
               for ( k = i; k < j; k++ )
                    if ( kmers_plus[k].second >= ngenome ) break;
               if ( k > i && j - k == 1 && k - i <= max_gmult )
               {    for ( int64_t m = i; m < k; m++ )
                    {    int e = kmers_plus[k].second - ngenome;
                         int g = kmers_plus[m].second;
                         int estart = kmers_plus[k].third; 
                         int gstart = kmers_plus[m].third;
                         int offset = gstart - estart;
                         if ( !Member( hits[e], make_pair( g, offset ) ) )
                         {     Bool mismatch = False;
                              for ( int l = K2 - 1; l >= K; l-- )
                              {    if ( hb.EdgeObject(e)[estart+l] 
                                             != genome[g][gstart+l] )
                                   {    mismatch = True;
                                        break;    }    }
                              if ( !mismatch ) 
                              {
                                   #pragma omp critical
                                   {    hits[e].push( 
                                             g, offset );    }    }    }    }    }
               i = j - 1;    }
          // cout << TimeSince(clock3) << " used traversing kmers" << endl;    
          }    

     // Copy alignments from one bubble branch to another, if the other is naked.
     // Copy alignments into the bubble if it makes sense.

     for ( int v = 0; v < hb.N( ); v++ )
     {    if ( hb.To(v).size( ) != 1 || hb.From(v).size( ) != 2 ) continue;
          if ( hb.From(v)[0] != hb.From(v)[1] ) continue;
          int w = hb.From(v)[0];
          if ( hb.To(w).size( ) != 2 || hb.From(w).size( ) != 1 ) continue;
          int e1 = hb.IFrom( v, 0 ), e2 = hb.IFrom( v, 1 );
          int d = hb.ITo( v, 0 );
          int f = hb.IFrom( w, 0 );
          if ( hits[e1].empty( ) && hits[e2].empty( ) && hits[d].solo( )
               && hits[f].solo( ) && hits[d][0].first == hits[f][0].first )
          {    int pos = hits[d][0].second + hb.Bases(d) - ( hb.K( ) - 1 );
               int Pos = hits[f][0].second + ( hb.K( ) - 1 );
               if ( Pos - pos >= hb.K( ) )
               {    hits[e1].push( hits[d][0].first, pos );
                    hits[e2].push( hits[d][0].first, pos );    }    }
          else
          {    if ( hits[e1].empty( ) ) swap( e1, e2 );
               if ( hits[e1].empty( ) || hits[e2].nonempty( ) ) continue;
               hits[e2] = hits[e1];    }    }

     // Adjudicate some alignments.

     for ( int e = 0; e < hb.E( ); e++ )
     {    int re = inv[e];
          if ( hits[e].size( ) + hits[re].size( ) < 2 || e == re ) continue;
          vec<int> errors( hits[e].size( ) + hits[re].size( ), 0 );
          for ( int pass = 1; pass <= 2; pass++ )
          {    int f = ( pass == 1 ? e : re );
               int start = ( pass == 1 ? 0 : hits[e].size( ) );
               for ( int j = 0; j < hits[f].isize( ); j++ )
               {    int g = hits[f][j].first;
                    int gpos = hits[f][j].second;
                    for ( int l = 0; l < hb.Bases(f); l++ )
                    {    if ( gpos+l < 0 || gpos+l >= genome[g].isize( ) ) 
                              errors[start+j]++;
                         else if ( hb.O(f)[l] != genome[g][gpos+l] ) 
                         {    errors[start+j]++;    }    }    }
               vec<int> ids( errors.size( ), vec<int>::IDENTITY );
               SortSync( errors, ids );
               if ( errors[0] == 0 && errors[1] > 0 )
               {    int id = ids[0];
                    if ( id < hits[e].isize( ) )
                    {    hits[e] = { hits[e][id] };
                         hits[re].clear( );    }
                    else
                    {    hits[re] = { hits[re][ id - hits[e].isize( ) ] };
                         hits[e].clear( );    }    }    }    }
     
     LogTime( clock, "aligning to genome" );    }

void AlignToGenomePerf( const HyperBasevector& hb, const vecbasevector& genome, 
     vec< triple< pair<int,int>, pair<int,int>, int > >& perfs,
     vec<perf_place>& places )
{
     // perfs: { (g,gstart), (e,estart), len ) }

     // Heuristics.  We were using K2 = 500, which is ~50% faster.  In principle
     // one could use a larger K then come back to fill in the small stuff.  Could
     // try K = 399.

     const int K2 = 200;

     // Computational performance heuristics.

     const int64_t batches = 100;
     const int K = 60; // must be no larger than K2

     // Go through four passes.

     double clock = WallClockTime( );
     int nobj = hb.EdgeObjectCount( ), ngenome = genome.size( );
     vec< vec< quad<int,int,int,int> > > hits; // (estart,g,gstart,len)
     hits.clear( );
     hits.resize(nobj);
     for ( int pi = 0; pi < 4; pi++ )
     {    
          // Build lookup table.

          double clock2 = WallClockTime( );
          vec< triple<kmer<K>,int,int> > kmers_plus;
          vec<int64_t> starts;
          starts.reserve( ngenome + nobj );
          starts.push_back(0);
          {    vec<int> counts( ngenome + nobj, 0 );
               #pragma omp parallel for
               for ( int i = 0; i < ngenome; i++ )
               {    const basevector& u = genome[i];
                    counts[i] = 0;
                    for ( int j = 0; j <= u.isize( ) - K2; j++ )
                         if ( u[j] == pi ) counts[i]++;    }
               #pragma omp parallel for
               for ( int i = 0; i < nobj; i++ )
               {    const basevector& u = hb.EdgeObject(i);
                    counts[ genome.size( ) + i ] = 0;
                    for ( int j = 0; j <= u.isize( ) - K2; j++ )
                         if ( u[j] == pi ) counts[ ngenome + i ]++;    }
               for ( int i = 0; i < ngenome + nobj; i++ )
                    starts.push_back( starts.back( ) + counts[i] );    }
          kmers_plus.resize( starts.back( ) );
          #pragma omp parallel for
          for ( int i = 0; i < ngenome; i++ )
          {    const basevector& u = genome[i];
               int count = 0;
               kmer<K> x;
               for ( int j = 0; j <= u.isize( ) - K2; j++ )
               {    if ( u[j] == pi )
                    {    int64_t r = starts[i] + count;
                         x.SetToSubOf( u, j ); 
                         kmers_plus[r].first = x;
                         kmers_plus[r].second = i, kmers_plus[r].third = j;    
                         count++;    }    }    }
          #pragma omp parallel for
          for ( int i = 0; i < nobj; i++ )
          {    const basevector& u = hb.EdgeObject(i);
               int count = 0;
               kmer<K> x;
               for ( int j = 0; j <= u.isize( ) - K2; j++ )
               {    if ( u[j] == pi )
                    {    int64_t r = starts[ngenome+i] + count;
                         x.SetToSubOf( u, j ); 
                         kmers_plus[r].first = x;
                         kmers_plus[r].second = ngenome + i, kmers_plus[r].third = j;
                         count++;    }    }    }
          double clock2b = WallClockTime( );
          ParallelSort(kmers_plus);

          // Traverse the kmers.

          double clock3a = WallClockTime( );
          vec<int64_t> bstart(batches+1);
          for ( int64_t i = 0; i <= batches; i++ )
               bstart[i] = ( (int64_t) kmers_plus.size( ) * i ) / batches;
          #pragma omp parallel for schedule(dynamic, 1)
          for ( int64_t i = 1; i < batches; i++ )
          {    int64_t& s = bstart[i];
               while( s > 0 && kmers_plus[s].first == kmers_plus[s-1].first )
               {    s--;    }    }
          double clock3 = WallClockTime( );
          #pragma omp parallel for schedule(dynamic, 1)
          for ( int64_t bi = 0; bi < batches; bi++ )
          for ( int64_t i = bstart[bi]; i < bstart[bi+1]; i++ )
          {    int64_t j, k;
               for ( j = i + 1; j < bstart[bi+1]; j++ )
                    if ( kmers_plus[j].first != kmers_plus[i].first ) break;
               for ( k = i; k < j; k++ )
                    if ( kmers_plus[k].second >= ngenome ) break;
               for ( int64_t m = i; m < k; m++ )
               for ( int64_t x = k; x < j; x++ )
               {    int e = kmers_plus[x].second - ngenome;
                    int g = kmers_plus[m].second;
                    int estart = kmers_plus[x].third, gstart = kmers_plus[m].third;
                    if ( estart-1 >= 0 && gstart-1 >= 0
                         && hb.EdgeObject(e)[estart-1] == genome[g][gstart-1] )
                    {    continue;    }
                    Bool mismatch = False;
                    for ( int l = K2 - 1; l >= K; l-- )
                    {    if ( hb.EdgeObject(e)[estart+l] != genome[g][gstart+l] )
                         {    mismatch = True;
                              break;    }    }
                    if ( !mismatch ) 
                    {    int len;
                         for ( len = K; 
                              estart + len < hb.EdgeObject(e).isize( ); len++ )
                         {    if ( gstart+len == genome[g].isize( )
                                   || hb.EdgeObject(e)[estart+len] 
                                        != genome[g][gstart+len] )
                              {    break;    }    }
                         #pragma omp critical
                         {    hits[e].push( estart, g, gstart, len );    }    }    }
               i = j - 1;    }    }    
     perfs.clear( );
     for ( int e = 0; e < hb.E( ); e++ )
     for ( int j = 0; j < hits[e].isize( ); j++ )
     {    perfs.push( make_pair( hits[e][j].second, hits[e][j].third ),
               make_pair( e, hits[e][j].first ), hits[e][j].fourth );    }

     // ParallelSort(perfs);

     /*
     // for unknown reasons, does not compile
     ParallelSort( perfs, [&](int i1,int i2)
     {    return make_pair( perfs[i1].first, -(perfs[i1].third) ) 
               < make_pair( perfs[i2].first, -(perfs[i2].third) );    } );
     */

     vec< pair< pair<int,int>, int > > keys( perfs.size( ) );
     for ( int i = 0; i < keys.isize( ); i++ )
          keys[i] = make_pair( perfs[i].first, -perfs[i].third );
     ParallelSortSync( keys, perfs );

     // Sloppy patch for an ordering problem.

     /*
     for ( int i = 0; i < (int) perfs.isize( ) - 1; i++ )
     {    if ( perfs[i].first != perfs[i+1].first ) continue;
          if ( perfs[i+1].second.second != 0 ) continue;
          if ( perfs[i+1].third <= perfs[i].third ) continue;
          swap( perfs[i], perfs[i+1] );    }
     */

     // Another sloppy patch for an ordering problem.

     vec<int> to_left, to_right;
     hb.ToLeft(to_left), hb.ToRight(to_right);
     /*
     for ( int i = 1; i < (int) perfs.isize( ) - 1; i++ )
     {    if ( perfs[i-1].first.first != perfs[i].first.first ) continue;
          if ( perfs[i].first.first != perfs[i+1].first.first ) continue;
          if ( perfs[i+1].second.second != 0 ) continue;
          if ( perfs[i].second.second == 0 ) continue;
          if ( to_right[ perfs[i-1].second.first ]
               != to_left[ perfs[i+1].second.first ] )
          {    continue;    }
          swap( perfs[i], perfs[i+1] );    }
     */

     // Chain the alignments.

     places.clear( );
     vec<Bool> used( perfs.size( ), False );
     for ( int i = 0; i < perfs.isize( ); i++ )
     {    if ( used[i] ) continue;
          used[i] = True;
          perf_place p;
          p.g = perfs[i].first.first;
          p.gstart = perfs[i].first.second;
          p.e.push_back( perfs[i].second.first );
          p.estart = perfs[i].second.second;
          p.len = perfs[i].third;
          int j, jl = i;
          for ( j = i + 1; j < perfs.isize( ); j++ )
          {    if ( perfs[j].first.first != p.g ) break;
               int e1 = perfs[jl].second.first, e2 = perfs[j].second.first;
               int len1 = perfs[jl].third, len2 = perfs[j].third;

               if ( perfs[jl].first.second + len1 - ( hb.K( ) - 1 ) 
                    < perfs[j].first.second )
               {    break;    }

               if ( used[j] ) continue;

               if ( perfs[j].first.second + len2 < p.gstart + p.len ) continue;

               if ( perfs[j].first.second + hb.K( ) - 1
                    != perfs[jl].first.second + len1 )
               {    continue;    }
               if ( to_right[e1] != to_left[e2] ) continue;
               if ( perfs[jl].second.second + len1 != hb.EdgeLengthBases(e1) ) 
                    continue;
               if ( perfs[j].second.second != 0 ) continue;

               used[j] = True;
               p.e.push_back( perfs[j].second.first );
               p.len += perfs[j].third - ( hb.K( ) - 1 );
               jl = j;    }
          places.push_back(p);    }
          // i = j - 1;    }

     // Remove subsumed alignments.  Make efficient later.

     vec<Bool> subsumed( places.size( ), False );
     for ( int i1 = 0; i1 < places.isize( ); i1++ )
     for ( int i2 = 0; i2 < places.isize( ); i2++ )
     {    if ( places[i1].G( ) != places[i2].G( ) ) continue;
          if ( places[i1].Gstart( ) <= places[i2].Gstart( )
               && places[i1].Gstop( ) > places[i2].Gstop( ) )
          {    subsumed[i2] = True;    }
          if ( places[i1].Gstart( ) < places[i2].Gstart( )
               && places[i1].Gstop( ) >= places[i2].Gstop( ) )
          {    subsumed[i2] = True;    }    }
     EraseIf( places, subsumed );

     // Try to extend up to the 'chromosome ends'.

     /*
     for ( int i = 0; i < places.isize( ); i++ )
     {    perf_place& p = places[i];
          int g = p.G( ), gstop = p.Gstop( );
          if ( i < places.isize( ) - 1 && g == places[i+1].G( ) ) continue;
          int ext = genome[g].isize( ) - gstop;
          PRINT3( g, ext, gstop ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          if ( ext == 0 ) continue;
          int e = p.E( ).back( ), estop = p.Estop(hb);
          PRINT2( estop, hb.EdgeLengthBases(e) ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          if ( estop < hb.EdgeLengthBases(e) ) continue;
          int v = to_right[e];
          vec<int> fs;
          for ( int j = 0; j < hb.From(v).isize( ); j++ )
          {    int f = hb.IFrom( v, j );
               PRINT(f); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
               if ( K-1 + ext-1 >= hb.EdgeLengthBases(f) ) continue;
               Bool mismatch = False;
               for ( int l = 0; l < ext; l++ )
               {    if ( genome[g][ gstop + l ] != hb.EdgeObject(f)[ K-1 + l ] )
                    {    mismatch = True;
                         cout << "mismatch" << endl; // XXXXXXXXXXXXXXXXXXXXXXXXXXXX
                         break;    }    }
               if ( !mismatch) fs.push_back(f);    }
          if ( fs.solo( ) ) 
          {    p.e.push_back( fs[0] );
               p.len += ext;    }    }
     */

     LogTime( clock, "aligning to genome perf" );    }

void ReroutePaths( const HyperBasevector& hb, const vec<int>& inv,
     ReadPathVec& paths, const vecbasevector& bases, const VecPQVec& quals )
{
     // Create indices.

     double clock = WallClockTime( );
     vec<int> to_left, to_right;
     hb.ToLeft(to_left), hb.ToRight(to_right);

     // Begin reroute.

     const int max_depth = 3;
     const int max_paths = 200;
     const int max_qsum = 100;

     int improveds = 0;
     #pragma omp parallel for schedule(dynamic, 1000)
     for ( int64_t id = 0; id < (int64_t) paths.size( ); id++ )
     {    ReadPath& p = paths[id];
          // Only consider full placements.

          if ( p.size( ) == 0 ) continue;
          if ( p.getOffset( ) < 0 ) continue;
          vec<int> s( p.size( ) );
          s[0] = p.getOffset( );
          for ( int j = 1; j < (int) p.size( ); j++ )
               s[j] = s[j-1] - hb.EdgeLengthKmers( p[j-1] );
          int n = bases[id].size( );
          if ( s.back( ) + n > hb.EdgeLengthBases( p.back( ) ) ) continue;

          // Find possible starts for the read.

          vec< pair<int,int> > starts = { make_pair( p[0], p.getOffset( ) ) };
          set< pair<int,int> > startsx;
          startsx.insert( make_pair( p[0], p.getOffset( ) ) );
          vec<int> depth = {0};
          for ( int i = 0; i < starts.isize( ); i++ )
          {    if ( depth[i] == max_depth ) continue;
               int e = starts[i].first, start = starts[i].second;
               int v = to_left[e], w = to_right[e];
               for ( int j = 0; j < hb.To(v).isize( ); j++ )
               {    int ex = hb.EdgeObjectIndexByIndexTo( v, j );
                    int startx = start + hb.EdgeLengthKmers(ex);
                    if ( !Member( startsx, make_pair( ex, startx ) ) )
                    {    starts.push( ex, startx );
                         startsx.insert( make_pair( ex, startx ) );
                         depth.push_back( depth[i] + 1 );    }    }
               for ( int j = 0; j < hb.From(w).isize( ); j++ )
               {    int ex = hb.EdgeObjectIndexByIndexFrom( w, j );
                    int startx = start - hb.EdgeLengthKmers(e);
                    if ( !Member( startsx, make_pair( ex, startx ) ) )
                    {    starts.push( ex, startx );
                         startsx.insert( make_pair( ex, startx ) );
                         depth.push_back( depth[i] + 1 );    }    }    }

          // Create initial paths.

          vec<ReadPath> ps;
          for ( int i = 0; i < starts.isize( ); i++ )
          {    if ( starts[i].second < 0 
                    || starts[i].second >= hb.EdgeLengthBases( starts[i].first ) )
               {    continue;    }
               ReadPath q;
               q.push_back( starts[i].first );
               q.setOffset( starts[i].second );
               ps.push_back(q);    }

          // Extend the paths.

          vec<Bool> to_delete( ps.size( ), False );
          for ( int i = 0; i < ps.isize( ); i++ )
          {    if ( i >= max_paths ) break;
               vec<int> s( ps[i].size( ) );
               s[0] = ps[i].getOffset( );
               for ( int j = 1; j < (int) ps[i].size( ); j++ )
                    s[j] = s[j-1] - hb.EdgeLengthKmers( ps[i][j-1] );
               int n = bases[id].size( );
               if ( s.back( ) + n <= hb.EdgeLengthBases( ps[i].back( ) ) ) continue;
               to_delete[i] = True;
               int v = to_right[ ps[i].back( ) ];
               for ( int j = 0; j < hb.From(v).isize( ); j++ )
               {    ReadPath r(ps[i]);
                    r.push_back( hb.EdgeObjectIndexByIndexFrom( v, j ) );
                    ps.push_back(r);    
                    to_delete.push_back(False);    }    }
          if ( ps.isize( ) > max_paths ) continue;
          EraseIf( ps, to_delete );

          // Score the paths.

          const int K = hb.K( );
          vec< pair<int,int> > qsum( ps.size( ), make_pair(0,0) );
          const basevector& r = bases[id];
          qvec qv;
          quals[id].unpack(&qv);
          for ( int i = 0; i < ps.isize( ); i++ )
          {    const ReadPath& q = ps[i];
               qsum[i].second = -q.size( );
               int start = q.getOffset( );
               basevector b = hb.EdgeObject( q[0] );
               for ( int l = 1; l < (int) q.size( ); l++ )
               {    b.resize( b.isize( ) - (K-1) );
                    b = Cat( b, hb.EdgeObject( q[l] ) );    }
               for ( int m = 0; m < r.isize( ); m++ )
                    if ( r[m] != b[start+m] ) qsum[i].first += qv[m];    }
          int qorig = qsum[0].first;
          SortSync( qsum, ps );
          for ( int i = 0; i < ps.isize( ); i++ )
               qsum[i].second = -qsum[i].second;
          Bool ok = False;
          for ( int j = 0; j < ps.isize( ); j++ )
          {    if ( ps[j] == p && qsum[j].first == qsum[0].first )
               {    ok = True;    }    }
          if (ok) continue;

          if ( qsum[0].first > max_qsum ) continue;

          #pragma omp critical
          {    improveds++;    }
     
          int ooo = qsum[0].first;
          while( ps.size( ) >= 2 && qsum[0] == qsum[1] )
          {    ps.SetToSubOf( ps, 2, ps.size( ) - 2 );
               qsum.SetToSubOf( qsum, 2, qsum.size( ) - 2 );    }

          vec<Bool> del( ps.size( ), False );
          for ( int j = 1; j < ps.isize( ); j++ )
          {    if ( qsum[j].first > qsum[0].first ) break;
               if ( qsum[j].second < qsum[0].second ) del[j] = True;    }
          EraseIf( ps, del );
          EraseIf( qsum, del );

          if ( ooo < qsum[0].first ) continue;

          Bool verbose = False;
          if (verbose)
          {    cout << "\n[" << id << "] ";
               for ( int j = 0; j < ps.isize( ); j++ )
               {    cout << ps[j].getOffset( ) << ":" << printSeq(ps[j])
                         << " --> " << qsum[j].first;
                    cout << endl;    
                    break;    }    }
          p = ps[0];    }

     cout << improveds << " paths improved by rerouting" << endl;
     LogTime( clock, "rerouting paths" );
     // cout << "\n" << Date( ) << ": done" << endl;    
          }

//    Tamp down hanging ends.
//
//    Under appropriate conditions, replace
//
//    v-----------------e1--------------->w
//    v------------e2--------->x (where no edges emanate from x)
//
//    by
//
//    v--------e1a------------>x----e1b-->w
//    v--------e2a------------>x
//
//    This has the effect of eliminating a hanging end by allowing it to
//    'continue on'.
//
//    Conditions:
//    (a) picture complete except for edges entering v, exiting w
//    (b) e2 matches e1 at end for >= 40 bases
//    (c) e2 mismatches e1 at <= 4 bases
//    (d) |e1| - |e2| + match > K - 1.

void Tamp( HyperBasevector& hb, vec<int>& inv, ReadPathVec& paths,
     const int max_shift )
{
     // Set up.  We build the paths index, but the updates to it are commented
     // out.  Since the bulk of the run time in this module is spent building the
     // paths index, it might be more efficient to pass the paths index to this
     // module, and have it update the index.

     double clock = WallClockTime( );
     VecULongVec paths_index;
     invert( paths, paths_index, hb.EdgeObjectCount( ) );
     int K = hb.K( );
     vec<int> to_left, to_right;
     hb.ToLeft(to_left), hb.ToRight(to_right);

     // Heuristics.

     const int min_match = 40;
     const int max_mismatches = 4;

     // Identify loci to be edited.

     int count = 0;
     vec< triple<int,int,int> > vj;
     vec<int> shift_vj;
     vec<Bool> touched( hb.EdgeObjectCount( ), False );
     for ( int v = 0; v < hb.N( ); v++ )
     {    if ( hb.From(v).size( ) != 2 ) continue;
          for ( int j = 0; j < 2; j++ )
          {    
               // Find candidate e1 and e2.

               int e1 = hb.EdgeObjectIndexByIndexFrom( v, j );
               int e2 = hb.EdgeObjectIndexByIndexFrom( v, 1-j );

               // Test them.

               int n1 = hb.EdgeLengthBases(e1), n2 = hb.EdgeLengthBases(e2);
               if ( n1 <= n2 ) continue;
               int x = to_right[e2], w = to_right[e1];
               if ( !hb.From(x).empty( ) || !hb.To(x).solo( ) ) continue;
               if ( !hb.To(w).solo( ) ) continue;
               const basevector &x1 = hb.EdgeObject(e1), &x2 = hb.EdgeObject(e2);
               if ( !IsUnique( v, x, w ) ) continue;

               int mis = 0, match = 0;
               for ( int l = n2 - 1; l >= 0; l-- )
               {    if ( x1[l] != x2[l] ) 
                    {    mis++;
                         if ( mis > max_mismatches ) break;    }
                    else if ( mis == 0 ) match++;    }
               int shift = 0;
               if ( max_shift == 0 && K - 1 - match < 0 ) continue;
               if ( max_shift == 0 && ( mis > max_mismatches || match < min_match ) )
                    continue;
               if ( max_shift > 0 ) // note does not use max_mismatches!!!!!!!!!!!!!
               {    vec<int> shifts;
                    for ( int s = -max_shift; s <= max_shift; s++ )
                    {    Bool mis = False;
                         for ( int l = n2-1; l > n2-1-min_match; l-- )
                         {    if ( l+s >= x1.isize( ) || x1[l+s] != x2[l] ) 
                              {    mis = True;
                                   break;    }    }
                         if ( !mis ) shifts.push_back(s);    }
                    if ( !shifts.solo( ) ) continue;
                    shift = shifts[0];    
                    match = min_match;    }

               if ( n1 - n2 - shift + match <= K - 1 ) continue;
               int re1 = inv[e1], re2 = inv[e2];
               if ( !IsUnique( e1, e2, re1, re2 ) ) continue; // not sure possible
               if ( touched[e1] || touched[e2] || touched[re1] || touched[re2] )
                    continue;
               touched[e1] = touched[e2] = touched[re1] = touched[re2] = True;
               vj.push( v, j, match );    
               shift_vj.push_back(shift);
               count++;    }    }

     // Edit loci.

     inv.resize( hb.EdgeObjectCount( ) + 4*count );
     // paths_index.reserve( hb.EdgeObjectCount( ) + 4*count );
     for ( int u = 0; u < vj.isize( ); u++ )
     {    int v = vj[u].first, j = vj[u].second, match = vj[u].third;
          int e1 = hb.EdgeObjectIndexByIndexFrom( v, j );
          int e2 = hb.EdgeObjectIndexByIndexFrom( v, 1-j );
          int n1 = hb.EdgeLengthBases(e1), n2 = hb.EdgeLengthBases(e2);
          int x = to_right[e2], w = to_right[e1];
          const basevector &x1 = hb.EdgeObject(e1), &x2 = hb.EdgeObject(e2);
          int re1 = inv[e1], re2 = inv[e2];
          int shift = shift_vj[u];

          // Fix the graph.  Note that we're not checking for interaction between
          // the edit and the rc edit, and this could lead to problems.

          basevector x2a = Cat( x2, basevector( x1, n2 + shift, K - 1 - match ) );
          basevector x1a = basevector( x1, 0, x2a.size( ) + shift );
          basevector x1b = basevector( x1, x2a.isize( ) - (K-1) + shift, 
               x1.isize( ) - ( x2a.isize( ) - (K-1) + shift ) );

          hb.DeleteEdgeFrom( v, j );
          hb.EdgeObjectMutable(e2) = x2a;
          int e1a = hb.AddEdge( v, x, x1a ), e1b = hb.AddEdge( x, w, x1b );

          // Fix the reverse complement instance on the graph.

          x2a.ReverseComplement( );
          x1a.ReverseComplement( ), x1b.ReverseComplement( );
          int rv = to_right[re1], rw = to_left[re1], rx = to_left[re2], rj;
          for ( rj = 0; rj < 2; rj++ )
               if ( hb.EdgeObjectIndexByIndexTo( rv, rj ) == re1 ) break;
          hb.DeleteEdgeTo( rv, rj );
          hb.EdgeObjectMutable(re2) = x2a;
          int re1a = hb.AddEdge( rx, rv, x1a ), re1b = hb.AddEdge( rw, rx, x1b );

          // Update inversion.

          inv[e1] = inv[re1] = -1;
          inv[e1a] = re1a, inv[re1a] = e1a, inv[e1b] = re1b, inv[re1b] = e1b;

          // Update paths and paths_index.

          for ( int l = 0; l < (int) paths_index[e1].size( ); l++ )
          {    ReadPath& p = paths[ paths_index[e1][l] ];
               for ( int m = 0; m < (int) p.size( ); m++ )
               {    if ( p[m] == e1 )
                    {    if ( m > 0 || p.getOffset( ) < hb.EdgeLengthBases(e1a) )
                         {    p[m] = e1a;
                              int p1a = p.getOffset( );
                              for ( int j = 0; j <= m; j++ )
                                   p1a -= hb.EdgeLengthKmers( p[j] );
                              if ( m < (int) p.size( ) - 1 || p1a >= 0 )
                              {    p.insert( p.begin( ) + m + 1 , e1b );
                                   m++;    }    }
                         else
                         {    p[m] = e1b;
                              p.setOffset( p.getOffset( ) 
                                   - hb.EdgeLengthKmers(e1a) );    }    }    }
               // int pi = paths_index[e1][l];
               // paths_index[e1a].push_back(pi);
               // paths_index[e1b].push_back(pi);    
                    }
          for ( int l = 0; l < (int) paths_index[re1].size( ); l++ )
          {    ReadPath& p = paths[ paths_index[re1][l] ];
               for ( int m = 0; m < (int) p.size( ); m++ )
               {    if ( p[m] == re1 )
                    {    if ( m > 0 || p.getOffset( ) < hb.EdgeLengthBases(re1b) )
                         {    p[m] = re1b;
                              int p1b = p.getOffset( );
                              for ( int j = 0; j <= m; j++ )
                                   p1b -= hb.EdgeLengthKmers( p[j] );
                              if ( m < (int) p.size( ) - 1 || p1b >= 0 )
                              {    p.insert( p.begin( ) + m + 1 , re1a );
                                   m++;    }    }
                         else
                         {    p[m] = re1a;
                              p.setOffset( p.getOffset( )
                                   - hb.EdgeLengthKmers(re1b) );    }    }    }
               // int pi = paths_index[re1][l];
               // paths_index[re1a].push_back(pi);
               // paths_index[re1b].push_back(pi);    
                    }
          // paths_index[e1].resize(0), paths_index[re1].resize(0);
          Bool verbose = False;
          if (verbose) PRINT2( e1, e2 );    }

     // Clean up.

     cout << count << " edges tamped down" << endl;
     LogTime( clock, "tamping" );
     Cleanup( hb, inv, paths );
}



void ExtendPairs60( HyperBasevector& hb, vec<int>& inv, ReadPathVec& paths )
{
     // Set up.

     double clock = WallClockTime( );
     vec<int> to_left, to_right;
     hb.ToLeft(to_left), hb.ToRight(to_right);
     const int K = 200;

     // Find pairs.

     vec< pair<int,int> > con;
     vec<Bool> ext( paths.size( ) / 2, False );
     #pragma omp parallel for
     for ( int64_t i = 0; i < (int64_t) paths.size( ); i += 2 )
     {    const ReadPath &p1 = paths[i], &p2 = paths[i+1];
          if ( p1.size( ) == 0 || p2.size( ) == 0 ) continue;
          vec<int> x1, x2;
          for ( int i = 0; i < (int) p1.size( ); i++ )
               x1.push_back( p1[i] ); 
          for ( int i = ( (int) p2.size( ) ) - 1; i >= 0; i-- )
               x2.push_back( inv[ p2[i] ] );
          if ( x1 == x2 ) continue;
          int v = to_right[ x1.back( ) ], w = to_left[ x2.front( ) ];
          if ( hb.From(v).nonempty( ) || hb.To(w).nonempty( ) ) continue;
          ext[i/2] = True;    }
     for ( int64_t i = 0; i < (int64_t) paths.size( ); i += 2 )
     {    if ( !ext[i/2] ) continue;
          const ReadPath &p1 = paths[i], &p2 = paths[i+1];
          vec<int> x1, x2;
          for ( int i = 0; i < (int) p1.size( ); i++ )
               x1.push_back( p1[i] ); 
          for ( int i = ( (int) p2.size( ) ) - 1; i >= 0; i-- )
               x2.push_back( inv[ p2[i] ] );
          con.push( x1.back( ), x2.front( ) );
          con.push( inv[ x2.front( ) ], inv[ x1.back( ) ] );    }

     // Find overlaps.

     Sort(con);
     vec< pair<int,int> > mix1;
     vec<int> mix2;
     for ( int i = 0; i < con.isize( ); i++ )
     {    int j = con.NextDiff(i);
          int e1 = con[i].first, e2 = con[i].second;
          if ( !( make_pair( e1, e2 ) <= make_pair( inv[e2], inv[e1] ) ) ) continue;

          // Test for overlap.

          const int over1 = 10;
          const int over2 = 20;
          vec<int> overs;
          for ( int o = over1; o < hb.K( ); o++ )
          {    if ( hb.EdgeObject(e1).Overlap( hb.EdgeObject(e2), o ) )
                    overs.push_back(o);    }
          if ( !( overs.solo( ) && overs[0] >= over2 ) ) continue;
          int o = overs[0];
          /*
          cout << "[" << j-i << "] " << con[i].first << ".." << con[i].second;
          cout << " over=" << o << " ";
          for ( int j = 0; j < o; j++ )
               cout << as_base( hb.EdgeObject(e2)[j] );
          cout << endl;
          */

          // Add bridge edge, and do the same for the reverse complement.

          int n1 = hb.EdgeLengthBases(e1), n2 = hb.EdgeLengthBases(e2);
          int k = hb.K( );
          basevector b = Cat( basevector( hb.EdgeObject(e1), n1 - (k-1), k-1 ),
               basevector( hb.EdgeObject(e2), o, (k-1) - o ) );
          int v = to_right[e1], w = to_left[e2];
          int x = hb.AddEdge( v, w, b );
          mix1.push( e1, e2 );
          mix2.push_back(x);
          int re1 = inv[e2], re2 = inv[e1];
          if ( re1 != e1 )
          {    b.ReverseComplement( );
               int rx = hb.AddEdge( to_right[re1], to_left[re2], b );
               mix1.push( re1, re2 );
               mix2.push_back(rx);
               inv.push_back( rx, x );    }
          else inv.push_back(x);

          i = j - 1;    }

     // Extend paths.

     SortSync( mix1, mix2 );
     for ( int64_t i = 0; i < (int64_t) paths.size( ); i += 2 )
     {    if ( !ext[i/2] ) continue;
          ReadPath &p1 = paths[i], &p2 = paths[i+1];
          vec<int> x1, x2;
          for ( int i = 0; i < (int) p1.size( ); i++ )
               x1.push_back( p1[i] ); 
          for ( int i = ( (int) p2.size( ) ) - 1; i >= 0; i-- )
               x2.push_back( inv[ p2[i] ] );
          int p = BinPosition( mix1, make_pair( x1.back( ), x2.front( ) ) );
          if ( p < 0 ) continue;
          int n1 = hb.EdgeLengthBases( x1[0] );
          for ( int i = 1; i < x1.isize( ); i++ )
               n1 += hb.EdgeLengthKmers( x1[i] );
          int n2 = hb.EdgeLengthBases( x2[0] );
          for ( int i = 2; i < x2.isize( ); i++ )
               n2 += hb.EdgeLengthKmers( x2[i] );
          if ( n1 < K ) 
          {    p1.push_back( mix2[p] );
               n1 += hb.EdgeLengthKmers( mix2[p] );    }
          if ( n2 < K ) 
          {    p2.push_back( inv[ mix2[p] ] );
               n2 += hb.EdgeLengthKmers( mix2[p] );    }
          for ( int i = 0; i < x2.isize( ); i++ )
          {    if ( n1 >= K ) break;
               p1.push_back( x2[i] );
               n1 += hb.EdgeLengthKmers( x2[i] );    }
          for ( int i = x1.isize( ) - 1; i >= 0; i-- )
          {    if ( n2 >= K ) break;
               p2.push_back( inv[ x1[i] ] );
               n2 += hb.EdgeLengthKmers( x1[i] );    }    }
          
     cout << TimeSince(clock) << " used extending pairs" << endl;    }

void PlacePartners( const HyperBasevector& hb, const vec<int>& inv,
     ReadPathVec& paths, const vecbasevector& bases, const VecPQVec& quals )
{
     // Set up.

     double clock = WallClockTime( );
     cout << Date( ) << ": placing partners" << endl;
     vec<int> to_left, to_right;
     hb.ToLeft(to_left), hb.ToRight(to_right);

     // Define heuristics.

     const int max_homes = 20;
     const int K = 40;
     const int min_gain = 5;

     // Find one-sided placements.

     cout << "\nnpids = " << paths.size( ) / 2 << endl;
     vec<int64_t> ids;
     for ( int64_t pid = 0; pid < (int64_t) paths.size( ) / 2; pid++ )
     {    int64_t id1 = 2*pid, id2 = 2*pid + 1;
          ReadPath &p1 = paths[id1], &p2 = paths[id2];
          if ( !( p1.size( ) == 0 ^ p2.size( ) == 0 ) ) continue;
          if ( p1.size( ) == 0 ) ids.push_back(id1);
          else ids.push_back(id2);    }
     PRINT2( ids.size( ), paths.size( ) );

     // Define edges on which the partner might land.

     cout << Date( ) << ": defining homes" << endl;
     vec<vec<int>> homes( ids.size( ) );
     for ( int i = 0; i < ids.isize( ); i++ )
     {    int64_t id1 = ids[i]; // UNPLACED
          int64_t id2 = ( (id1 % 2 == 0) ? id1 + 1 : id1 - 1 ); // PLACED
          const ReadPath& p2 = paths[id2];
          int e2 = p2[0], start2 = p2.getOffset( );
          int e1 = inv[e2];
          if ( e1 < 0 ) continue;
          int stop1 = hb.EdgeLengthBases(e2) - start2;
          vec< pair<int,int> > exts;
          exts.push( e1, bases[e1].isize( ) - stop1 );
          int nhomes = 1;
          for ( int j = 0; j < exts.isize( ); j++ )
          {    if ( exts[j].second <= 0 ) continue;
               int v = to_left[ exts[j].first ];
               for ( int l = 0; l < hb.To(v).isize( ); l++ )
               {    int e = hb.EdgeObjectIndexByIndexTo( v, l );
                    int len = hb.EdgeLengthKmers(e);
                    Bool superfluous = False;
                    for ( int m = 0; m < exts.isize( ); m++ )
                    {    if ( exts[m].first == e )
                         {    if ( exts[j].second - len <= exts[m].second )
                              {    superfluous = True;
                                   break;    }
                              else 
                              {    exts[m].first = -1;
                                   nhomes--;
                                   break;    }    }    }
                    if (superfluous) continue;
                    exts.push( e, exts[j].second - len );
                    if ( ++nhomes > max_homes ) break;    }
               if ( nhomes > max_homes ) break;    }
          if ( nhomes > max_homes ) continue;
          for ( int j = 0; j < exts.isize( ); j++ )
               if ( exts[j].first >= 0 ) homes[i].push_back( exts[j].first );
          Sort( homes[i] );    }

     // Form kmer lookup table for assembly and unplaced partners, then find homes
     // for partners.

     cout << Date( ) << ": forming all" << endl;
     vec<int> hall;
     for ( int i = 0; i < ids.isize( ); i++ )
          hall.append( homes[i] );
     UniqueSort(hall);
     PRINT2( hall.size( ), hb.EdgeObjectCount( ) );
     vecbasevector all;
     all.reserve( hall.isize( ) + ids.isize( ) );
     for ( int i = 0; i < hall.isize( ); i++ )
          all.push_back( hb.EdgeObject( hall[i] ) );
     for ( int i = 0; i < ids.isize( ); i++ )
          all.push_back( bases[ ids[i] ] );
     vec< vec< pair<int,int> > > hits( ids.size( ) );
     cout << Date( ) << ": forming lookup table" << endl;
     cout << TimeSince(clock) << " used initially" << endl;
     vec<String> prefixes = { "AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA",
               "GC", "GG", "GT", "TA", "TC", "TG", "TT" };
     cout << "memory usage = " << ToStringAddCommas( MemUsageBytes( ) ) << endl;
     // prefixes = { "A", "C", "G", "T" };
     for ( int pass = 0; pass < prefixes.isize( ); pass++ )
     {    double lclock = WallClockTime( );
          vec< triple<kmer<K>,int,int> > kmers_plus;
          MakeKmerLookup0Pre( all, prefixes[pass], kmers_plus );
          cout << TimeSince(lclock) << " used forming lookup table" << endl;

          // Find homes for partners.

          cout << Date( ) << ": finding homes for partners" << endl;
          double fclock = WallClockTime( );
          const int blocks = 1000;
          vec<int64_t> starts( blocks + 1 );
          for ( int j = 0; j <= blocks; j++ )
               starts[j] = ( j * kmers_plus.jsize( ) ) / blocks;
          for ( int j = blocks - 1; j >= 1; j-- )
          {    while( starts[j] > 0 && kmers_plus[ starts[j] ].first
                    == kmers_plus[ starts[j] - 1 ].first )
               {    starts[j]--;    }    }
          #pragma omp parallel for
          for ( int b = 0; b < blocks; b++ )
          {    vec< triple<int,int,int> > x;
               for ( int64_t i = starts[b]; i < starts[b+1]; i++ )
               {    int64_t j;
                    for ( j = i + 1; j < starts[b+1]; j++ )
                         if ( kmers_plus[j].first != kmers_plus[i].first ) break;
                    for ( int64_t k1 = i; k1 < j; k1++ )
                    for ( int64_t k2 = i; k2 < j; k2++ )
                    {    if ( kmers_plus[k1].second >= hall.isize( ) ) continue;
                         if ( kmers_plus[k2].second < hall.isize( ) ) continue;
                         int eid = hall[ kmers_plus[k1].second ];
                         int uid = kmers_plus[k2].second - hall.isize( );
                         if ( !BinMember( homes[uid], eid ) ) continue;
                         int offset = kmers_plus[k1].third - kmers_plus[k2].third;
                         x.push( uid, eid, offset );    }
                    i = j - 1;    }
               #pragma omp critical
               {    for ( int i = 0; i < x.isize( ); i++ )
                    {    hits[ x[i].first ].push( 
                              x[i].second, x[i].third );    }    }    }
          cout << TimeSince(fclock) << " used finding homes" << endl;    }
     double sclock = WallClockTime( );
     for ( int64_t i = 0; i < ids.jsize( ); i++ )
          UniqueSort( hits[i] );
     cout << TimeSince(sclock) << " used sorting" << endl;

     // Try to eliminate negative offsets.

     double tclock = WallClockTime( );
     for ( int64_t i = 0; i < ids.jsize( ); i++ )
     {    while(1)
          {    Bool changed = False;
               for ( int j = 0; j < hits[i].isize( ); j++ )
               {    int v = to_left[ hits[i][j].first ];
                    int start = hits[i][j].second;
                    if ( start < 0 && hb.To(v).nonempty( ) )
                    {    changed = True;
                         for ( int l = 0; l < hb.To(v).isize( ); l++ )
                         {    int e = hb.EdgeObjectIndexByIndexTo( v, l );
                              int len = hb.EdgeLengthKmers(e);
                              if ( l == 0 ) hits[i][j] = make_pair( e, start + len );
                              else 
                              {    hits[i].push( e, 
                                        start + len );    }    }    }    }
               UniqueSort( hits[i] );
               if ( !changed ) break;    }    }

     // If we can shift a singleton left, do so.

     for ( int64_t i = 0; i < ids.jsize( ); i++ )
     {    if ( hits[i].size( ) != 1 ) continue;
          int id = ids[i], e = hits[i][0].first;
          int start = hits[i][0].second;
          if ( start >= hb.K( ) - 1 ) continue;
          int v = to_left[e];
          if ( !hb.To(v).solo( ) ) continue;
          int ep = hb.EdgeObjectIndexByIndexTo( v, 0 );
          hits[i][0] = make_pair( ep, start + hb.EdgeLengthKmers(ep) );    }

     // If we can shift a pair left to a single placement, do so.

     for ( int64_t i = 0; i < ids.jsize( ); i++ )
     {    if ( hits[i].size( ) != 2 ) continue;
          int id = ids[i], e1 = hits[i][0].first, e2 = hits[i][1].first;
          int start1 = hits[i][0].second, start2 = hits[i][1].second;
          if ( to_left[e1] != to_left[e2] || start1 != start2 ) continue;
          if ( start1 > hb.K( ) - 1 ) continue;
          int v = to_left[e1];
          if ( !hb.To(v).solo( ) ) continue;
          int e = hb.EdgeObjectIndexByIndexTo( v, 0 );
          hits[i].resize(1);
          hits[i][0] = make_pair( e, start1 + hb.EdgeLengthKmers(e) );    }

     // In the case where we have placements on both branches of a bubble, score
     // them.

     auto qvItr = quals.begin();
     for ( int64_t i = 0; i < ids.jsize( ); i++ )
     {    if ( hits[i].size( ) != 2 ) continue;
          int id = ids[i], e1 = hits[i][0].first, e2 = hits[i][1].first;
          qvec const& qv = qvItr[id];
          int start1 = hits[i][0].second, start2 = hits[i][1].second;
          if ( start1 < 0 || start2 < 0 ) continue;
          if ( to_right[e1] != to_right[e2] ) continue;
          if ( hb.EdgeLength(e1) - start1 != hb.EdgeLength(e2) - start2 ) continue;
          vec<int> qsum( 2, 0 );
          for ( int l = 0; 
               l < Min( bases[id].isize( ), hb.EdgeLengthBases(e1) - start1 ); l++ )
          {    if ( bases[id][l] != hb.EdgeObject(e1)[start1+l] )
                    qsum[0] += qv[l];
               if ( bases[id][l] != hb.EdgeObject(e2)[start2+l] )
                    qsum[1] += qv[l];    }
          vec<Bool> to_delete( hits[i].size( ), False );
          if ( qsum[0] <= qsum[1] - min_gain ) to_delete[1] = True;
          if ( qsum[1] <= qsum[0] - min_gain ) to_delete[0] = True;
          if ( Sum(to_delete) > 0 ) EraseIf( hits[i], to_delete );    }

     // If we can shift a singleton left, do so.

     for ( int64_t i = 0; i < ids.jsize( ); i++ )
     {    if ( hits[i].size( ) != 1 ) continue;
          int id = ids[i], e = hits[i][0].first;
          int start = hits[i][0].second;
          if ( start >= hb.K( ) - 1 ) continue;
          int v = to_left[e];
          if ( !hb.To(v).solo( ) ) continue;
          int ep = hb.EdgeObjectIndexByIndexTo( v, 0 );
          hits[i][0] = make_pair( ep, start + hb.EdgeLengthKmers(ep) );    }

     // Extend paths.

     int placed = 0;
     for ( int64_t i = 0; i < ids.jsize( ); i++ )
     {    if ( !hits[i].solo( ) ) continue;
          int e = hits[i][0].first;
          int start = hits[i][0].second;
          if ( start >= 0 && start < hb.EdgeLengthBases(e) )
          {    SerfVec<int> x;
               x.push_back(e);
               ReadPath p;
               (SerfVec<int>&) p = x;
               p.setOffset(start);
               const int EXT_MODE = 1;
               ExtendPath( p, ids[i], hb, to_right, bases[i], qvItr[i],
                    min_gain, False, EXT_MODE );
               paths[ ids[i] ] = p;
               placed++;
               // cout << "\nplaces for " << i << " = " << ids[i] << ":\n";
               // cout << p.getOffset( ) << ":" << printSeq(p) << "\n";    
                    }    }

     // cout << "\n";
     PRINT3( placed, ids.size( ), bases.size( ) );

     cout << TimeSince(tclock) << " used in tail" << endl;
     cout << TimeSince(clock) << " used in total placing partners" << endl;    }

namespace
{
unsigned const KLEN = 28;

size_t findInterestingReadIds( HyperBasevector const& hbv,
                                ReadPathVec const& paths,
                                vecbvec const& reads,
                                vec<size_t>* pReadIds )
{
    vec<int> endEdges;
    if ( true )
    {
        int const MAX_DIST = 10000000; // dangerous!
        int const GOOD_DIST = 500;
        vec<int> D;
        int hbk = hbv.K();
        DistancesToEndFunc( hbv, [hbk](bvec const& bv){return bv.size()-hbk+1;},
                                MAX_DIST, True, D );
        hbv.ToRight(endEdges);
        for ( int& vvv : endEdges )
            vvv = (D[vvv] <= GOOD_DIST);
    }

    size_t nReads = reads.size();
    pReadIds->clear();
    pReadIds->reserve(nReads/200);
    size_t nKmers = 0;
    for ( size_t readId=0; readId != nReads; ++readId )
    {
        if ( paths[readId].empty() )
        {
            ReadPath const& matePath = paths[readId^1];
            if ( !matePath.empty() && endEdges[matePath.back()] )
            {
                bvec const& read = reads[readId];
                if ( read.size() >= KLEN )
                {
                    pReadIds->push_back(readId);
                    nKmers += read.size()-KLEN+1;
                }
            }
        }
    }
    return nKmers;
}

class Loc
{
public:
    Loc()=default;
    Loc( unsigned readId, int offset )
    : mReadId(readId), mOffset(offset) {}

    unsigned getReadId() const { return mReadId; }
    int getOffset() const { return mOffset; }

    void nextOffset() { mOffset += 1; }

    friend bool operator<( Loc const& loc1, Loc const& loc2 )
    { if ( loc1.mReadId != loc2.mReadId ) return loc1.mReadId < loc2.mReadId;
      return loc1.mOffset < loc2.mOffset; }
    friend bool operator==( Loc const& loc1, Loc const& loc2 )
    { return loc1.mReadId == loc2.mReadId && loc1.mOffset == loc2.mOffset; }
    friend ostream& operator<<( ostream& os, Loc const& loc )
    { return os << loc.getReadId() << '[' << loc.getOffset() << ']'; }

private:
    unsigned mReadId;
    int mOffset;
};

typedef KMer<KLEN> Kmer;
class KmerLoc : public Kmer
{
public:
    KmerLoc()=default;
    template <class Itr>
    KmerLoc( Itr itr, unsigned readId )
    : Kmer(itr), mLoc(readId,0u) {}

    void nextLoc( unsigned char nextBase )
    { toSuccessor(nextBase); mLoc.nextOffset(); }

    Loc const& getLoc() const { return mLoc; }

private:
    Loc mLoc;
};

class KmerLocs : public Kmer
{
public:
    KmerLocs()=default;

    KmerLocs( Kmer const& kmer, Loc* pLocs, Loc* pEnd )
    : Kmer(kmer), mpLocs(pLocs), mNLocs(pEnd-pLocs), mELocs(0)
    { if ( mNLocs != 1 ) std::sort(pLocs,pEnd);
      else
      { ForceAssertGe(sizeof(Loc*),sizeof(Loc));
        *reinterpret_cast<Loc*>(&mpLocs) = *pLocs; } }

    KmerLocs( KmerLoc const& kloc )
    : Kmer(kloc), mNLocs(1u), mELocs(0)
    { ForceAssertGe(sizeof(Loc*),sizeof(Loc));
      *reinterpret_cast<Loc*>(&mpLocs) = kloc.getLoc(); }

    size_t getNLocs() const { return mNLocs; }

    Loc const* locsBegin() const
    { return mNLocs!=1 ? mpLocs : reinterpret_cast<Loc const*>(&mpLocs); }

    Loc const* locsEnd() const
    { return locsBegin()+mNLocs; }

    unsigned getELocs() const { return mELocs; }
    void setELocs( unsigned eLocs ) { mELocs = eLocs; }

    size_t getTotalLocs() const { return mNLocs+mELocs; }

private:
    Loc* mpLocs;
    unsigned mNLocs;
    unsigned mELocs;
};

typedef HashSet<KmerLocs,Kmer::Hasher,std::equal_to<Kmer>> Dict;

class MREReadProc
{
public:
    MREReadProc( vec<size_t> const& ids, vecbvec const& reads, Mempool& alloc,
                    size_t maxMultiplicity, Dict* pDict )
    : mIds(ids), mReads(reads), mAlloc(alloc),
      mMaxMultiplicity(maxMultiplicity), mDict(*pDict),
      mpNext(nullptr), mpEnd(nullptr)
    {}

    template <class OItr>
    void map( size_t idxId, OItr oItr )
    { size_t readId = mIds[idxId];
      bvec const& read = mReads[readId];
      if ( read.size() < KLEN ) return;
      KmerLoc loc(read.begin(),readId);
      *oItr = loc; ++oItr;
      for ( auto itr=read.begin(KLEN),end=read.end(); itr != end; ++itr )
      { loc.nextLoc(*itr); *oItr = loc; ++oItr; } }

    void reduce( KmerLoc const* beg, KmerLoc const* end )
    { size_t nLocs = end-beg;
      if ( nLocs > mMaxMultiplicity )
        return;
      if ( nLocs == 1 )
        mDict.insertUniqueValue(KmerLocs(*beg));
      else
      { Loc* pLocs = alloc(nLocs);
        KmerLocs kLocs(*beg,pLocs,pLocs+nLocs);
        while ( beg != end ) *pLocs++ = beg++->getLoc();
        mDict.insertUniqueValue(kLocs); } }

    KmerLoc* overflow( KmerLoc* beg, KmerLoc* end )
    { return beg+std::min(size_t(end-beg),mMaxMultiplicity+1); }

private:
    Loc* alloc( size_t nnn )
    { if ( mpNext+nnn > mpEnd )
      { size_t nAlloc = std::max(10000ul,nnn);
        mpNext = static_cast<Loc*>(mAlloc.allocate(nAlloc*sizeof(Loc),4));
        mpEnd = mpNext+nAlloc; }
      Loc* result = mpNext; mpNext += nnn; return result; }

    vec<size_t> const& mIds;
    vecbvec const& mReads;
    Mempool& mAlloc;
    size_t mMaxMultiplicity;
    Dict& mDict;
    Loc* mpNext;
    Loc* mpEnd;
};

typedef MapReduceEngine<MREReadProc,KmerLoc,Kmer::Hasher,std::less<Kmer>> RMRE;

class MREEdgeProc
{
public:
    MREEdgeProc( HyperBasevector const& hbv, Dict* pDict )
    : mHBV(hbv), mDict(*pDict)
    {}

    template <class OItr>
    void map( size_t edgeId, OItr oItr )
    { bvec const& edge = mHBV.EdgeObject(edgeId);
      if ( edge.size() < KLEN ) return;
      Kmer kmer(edge.begin());
      *oItr = kmer; ++oItr;
      for ( auto itr=edge.begin(KLEN),end=edge.end(); itr != end; ++itr )
      { kmer.toSuccessor(*itr); *oItr = kmer; ++oItr; } }

    void reduce( Kmer const* beg, Kmer const* end )
    { bumpCount(*beg,end-beg); }

    Kmer* overflow( Kmer* beg, Kmer* end )
    { bumpCount(*beg,end-beg); return beg; }

private:
    void bumpCount( Kmer const& kmer, size_t count )
    { KmerLocs const* pLocs = mDict.lookup(kmer);
      if ( pLocs )
          const_cast<KmerLocs*>(pLocs)->setELocs(pLocs->getELocs()+count); }

    HyperBasevector const& mHBV;
    Dict& mDict;
};

typedef MapReduceEngine<MREEdgeProc,Kmer,Kmer::Hasher> EMRE;

class EdgeProc
{
    static int const WINDOW = 60;
    static int const MAX_MISMATCHES = 4;
    static unsigned char const TRUSTED_QUAL = 30;
public:
    static int const NOT_AN_EDGE = -1;

    EdgeProc( HyperBasevector const& hbv, Dict const& dict,
                vecbvec const& reads, VecPQVec const& quals,
                ReadPathVec& paths )
    : mHBV(hbv), mDict(dict), mReads(reads), mQuals(quals), mPaths(paths) {}

    void operator()( size_t edgeId )
    { bvec const& edge = mHBV.EdgeObject(edgeId);
      if ( edge.size() < KLEN ) return;
      Kmer kmer(edge.begin());
      int eOffset = 0;
      mLocs.clear();
      KmerLocs const* pLocs;
      if ( (pLocs = mDict.lookup(kmer)) )
        addLocs(*pLocs,eOffset);
      for ( auto itr=edge.begin(KLEN),end=edge.end(); itr != end; ++itr )
      { kmer.toSuccessor(*itr); eOffset += 1;
        if ( (pLocs = mDict.lookup(kmer)) )
          addLocs(*pLocs,eOffset); }
      std::sort(mLocs.begin(),mLocs.end());
      mLocs.erase(std::unique(mLocs.begin(),mLocs.end()),mLocs.end());
      for ( auto const& loc : mLocs )
      { ReadPath& path = mPaths[loc.getReadId()];
        if ( !path.empty() && path[0] == NOT_AN_EDGE )
          continue; // read placement already known not to be unique
        if ( isGood(edgeId,loc) )
        { static SpinLockedData gLock;
          SpinLocker locker(gLock);
          if ( !path.empty() ) path[0] = NOT_AN_EDGE; // mark as not unique
          else
          { path.push_back(edgeId);
            path.setOffset(-loc.getOffset()); } } }
      mLocs.clear(); }

    void cleanAmbiguousPlacements( vec<size_t> const& readIds )
    { for ( size_t readId : readIds )
      { ReadPath& path = mPaths[readId];
        if ( !path.empty() && path[0] == EdgeProc::NOT_AN_EDGE )
        { path.clear(); path.setOffset(0); } } }

private:
    void addLocs( KmerLocs const& locs, int eOffset )
    { auto lItr=locs.locsBegin(), lEnd=locs.locsEnd();
      auto beg=mLocs.begin(),end=mLocs.end();
      while ( beg != end && lItr != lEnd )
      { --lEnd; Loc loc(lEnd->getReadId(),lEnd->getOffset()-eOffset);
        if ( !(loc == *--end) ) { ++lEnd; break; } }
      for ( ; lItr != lEnd; ++lItr )
        mLocs.push_back(Loc(lItr->getReadId(),lItr->getOffset()-eOffset)); }

    // a loc is good if there are no high-quality mismatches, and there's at
    // least one window with few enough mismatches
    bool isGood( size_t edgeId, Loc const& loc )
    { size_t readId = loc.getReadId();
      bvec const& read = mReads[readId];
      mQuals[readId].unpack(&mQVec);
      bvec const& edge = mHBV.EdgeObject(edgeId);
      int offset = -loc.getOffset();
      auto rBeg = read.begin(), rEnd = read.end();
      auto eBeg = edge.begin(), eEnd = edge.end();
      auto qItr = mQVec.begin();
      if ( offset >= 0 ) eBeg += offset;
      else { rBeg -= offset; qItr -= offset; }
      if ( eEnd-eBeg < WINDOW || rEnd-rBeg < WINDOW ) return false;
      auto rItr = rBeg, eItr = eBeg;
      int misMatches = 0;
      for ( auto rWin=rBeg+WINDOW; rItr != rWin; ++rItr,++eItr,++qItr )
        if ( *rItr != *eItr )
        { if ( *qItr >= TRUSTED_QUAL ) return false;
          misMatches += 1; }
      bool good = misMatches <= MAX_MISMATCHES;
      for ( ; rItr != rEnd && eItr != eEnd; ++rItr,++eItr,++qItr,++rBeg,++eBeg )
      { if ( *rItr != *eItr )
        { if ( *qItr >= TRUSTED_QUAL ) return false;
          misMatches += 1; }
        if ( *rBeg != *eBeg ) misMatches -= 1;
        if ( misMatches <= MAX_MISMATCHES ) good = true; }
      return good; }

    HyperBasevector const& mHBV;
    Dict const& mDict;
    vecbvec const& mReads;
    VecPQVec const& mQuals;
    ReadPathVec& mPaths;
    vec<Loc> mLocs;
    qvec mQVec;
};

}

void PartnersToEnds( const HyperBasevector& hbv, ReadPathVec& paths,
                        const vecbasevector& reads, const VecPQVec& quals )
{
    double clock = WallClockTime();

    // find unplaced partners of reads near sinks

    vec<size_t> readIds;
    cout << Date( ) << ": finding interesting reads" << endl;
    cout << Date( ) << ": memory in use = " << MemUsageGBString( )
          << ", peak = " << PeakMemUsageGBString( ) << endl;
    size_t nKmers = findInterestingReadIds(hbv,paths,reads,&readIds);
    size_t nReads = readIds.size();
    if ( nReads == 0 ) return;

    // kmerize those reads, and reduce them into a dictionary of KmerLocs

    cout << Date( ) << ": building dictionary" << endl;
    cout << Date( ) << ": memory in use = " << MemUsageGBString( )
          << ", peak = " << PeakMemUsageGBString( ) << endl;
    Dict* pDict = new Dict(nKmers);
    Mempool locsAlloc;
    size_t const MAX_MULTIPLICITY = 80;
    cout << Date( ) << ": reducing" << endl;
    cout << Date( ) << ": memory in use = " << MemUsageGBString( )
          << ", peak = " << PeakMemUsageGBString( ) << endl;

    {   RMRE rmre(MREReadProc(readIds,reads,locsAlloc,MAX_MULTIPLICITY,pDict));
        rmre.run(nKmers,0ul,nReads,RMRE::VERBOSITY::QUIET);    }

    // kmerize edges, setting the multiplicity for existing dictionary entries

    cout << Date( ) << ": kmerizing" << endl;
    cout << Date( ) << ": memory in use = " << MemUsageGBString( )
          << ", peak = " << PeakMemUsageGBString( ) << endl;
    {   EMRE emre(MREEdgeProc(hbv,pDict));
        auto const& edges = hbv.Edges();
        size_t edgeKmers = kmerCount(edges.begin(),edges.end(),KLEN);
        emre.run(edgeKmers,0ul,size_t(hbv.E()),EMRE::VERBOSITY::QUIET);    }

    // remove dictionary entries having too great a kmer multiplicity

    cout << Date( ) << ": cleaning" << endl;
    cout << Date( ) << ": memory in use = " << MemUsageGBString( )
          << ", peak = " << PeakMemUsageGBString( ) << endl;
    pDict->remove_if([]( KmerLocs const& kLocs )
                     { return kLocs.getTotalLocs() > MAX_MULTIPLICITY; });

    // find a uniquely aligning edge and path the read on that edge

    cout << Date( ) << ": finding uniquely aligning edges" << endl;
    cout << Date( ) << ": memory in use = " << MemUsageGBString( )
          << ", peak = " << PeakMemUsageGBString( ) << endl;
    EdgeProc proc(hbv,*pDict,reads,quals,paths);
    parallelForBatch(0,hbv.E(),100,proc);
    proc.cleanAmbiguousPlacements(readIds);
    delete pDict;
}

void PartnersToEndsOld( const HyperBasevector& hb, ReadPathVec& paths,
                       const vecbasevector& bases, const VecPQVec& quals )
{
     // Set up.

     vec<int> to_left, to_right;
     hb.ToLeft(to_left), hb.ToRight(to_right);
     int64_t npids = bases.size( ) / 2;

     // Heuristics.

     const int K = 28;
     const int64_t max_mult = 80;
     const int batches = 32;
     const int qtop = 30;
     const int min_len = 60;
     const int max_diffs = 4;

     // Compute distance to end for each vertex.

     const int max_dist = 10000000; // dangerous!
     vec<int> D;
     cout << Date( ) << ": computing distances to end" << endl;
     int hbk = hb.K();
     DistancesToEndFunc( hb, [hbk]( bvec const& bv ){ return bv.size()-hbk+1; },
           max_dist, True, D );

     // Identify incompletely placed pairs, with the placed read landing near an end.

     cout << Date( ) << ": identifying unplaced partners" << endl;
     vec<Bool> un( bases.size( ), False );
     const int min_qual = 1000;
     const int max_prox = 500;
     #pragma omp parallel for
     for ( int64_t pid = 0; pid < npids; pid++ )
     {    int64_t id1 = 2*pid, id2 = 2*pid+1;
          const ReadPath &p1 = paths[id1], &p2 = paths[id2];
          if ( p1.size( ) > 0 && p2.size( ) > 0 ) continue;
          if ( p1.size( ) == 0 && p2.size( ) == 0 ) continue;
          if ( p1.size( ) > 0 && D[ to_right[ p1.back( ) ] ] > max_prox ) continue;
          if ( p2.size( ) > 0 && D[ to_right[ p2.back( ) ] ] > max_prox ) continue;
          if ( p1.size( ) == 0 ) un[id1] = True;
          if ( p2.size( ) == 0 ) un[id2] = True;    }

     // Build "trash".

     vecbasevector trashb;
     trashb.reserve( Sum(un) + hb.E( ) );
     cout << Date( ) << ": building trash" << endl;
     int count = 0;
     vec<int> ids;
     for ( int64_t id = 0; id < 2*npids; id++ )
     {    if ( un[id] )
          {    ids.push_back(id);
               trashb.push_back( bases[id] );
               count++;    }    }
     PRINT2( count, npids );
     cout << PERCENT_RATIO( 3, count, npids ) << endl;
     for ( int e = 0; e < hb.E( ); e++ )
          trashb.push_back( hb.EdgeObject(e) );

     // Make kmers.

     cout << Date( ) << ": making kmers" << endl;
     vec< triple<kmer<K>,int,int> > kmers_plus;
     MakeKmerLookup0( trashb, kmers_plus );

     // Find friends.

     cout << Date( ) << ": finding friends" << endl;
     vec<int64_t> bstart(batches+1);
     for ( int64_t i = 0; i <= batches; i++ )
          bstart[i] = ( (int64_t) kmers_plus.size( ) * i ) / batches;
     #pragma omp parallel for schedule(dynamic, 1)
     for ( int64_t i = 1; i < batches; i++ )
     {    int64_t& s = bstart[i];
          while( s > 0 && kmers_plus[s].first == kmers_plus[s-1].first )
          {    s--;    }    }
     ReportPeakMem( );

     vec< vec< triple<int,int,int> > > places(batches);

     cout << Date( ) << ": start loop" << endl;
     #pragma omp parallel for schedule(dynamic, 1)
     for ( int64_t bi = 0; bi < batches; bi++ )
     {    for ( int64_t i = bstart[bi]; i < bstart[bi+1]; i++ )
          {    int64_t j;
               for ( j = i + 1; j < bstart[bi+1]; j++ )
                    if ( kmers_plus[j].first != kmers_plus[i].first ) break;
               if ( j - i >= 2 && j - i <= max_mult )
               {    for ( int64_t k1 = i; k1 < j; k1++ )
                    for ( int64_t k2 = i; k2 < j; k2++ )
                    {    int id1 = kmers_plus[k1].second; // read index
                         if ( id1 >= ids.isize( ) ) continue;
                         int id2 = kmers_plus[k2].second; // edge index
                         if ( id2 < ids.isize( ) ) continue;
                         id2 -= ids.isize( );
                         int pos1 = kmers_plus[k1].third; 
                         int pos2 = kmers_plus[k2].third;
                         int offset = pos2 - pos1;
                         places[bi].push( ids[id1], id2, offset );    }    }
               i = j - 1;    }    
          UniqueSort( places[bi] );    

          // Kill placements having hq diffs or not having a long-enough near-perfect
          // stretch.

          vec<Bool> del( places[bi].size( ), False );
          qvec qv;
          for ( int j = 0; j < places[bi].isize( ); j++ )
          {    int id1 = places[bi][j].first, id2 = places[bi][j].second;
               quals[id1].unpack(&qv);
               int offset = -places[bi][j].third;
               const basevector& rd2 = hb.EdgeObject(id2);
               Bool bad = False, good = False;
               int mis_count = 0;
               vec<Bool> mis;
               for ( int p1 = 0; p1 < bases[id1].isize( ); p1++ )
               {    int p2 = p1 - offset;
                    if ( p2 < 0 || p2 >= rd2.isize( ) ) continue;
                    if ( bases[id1][p1] != rd2[p2] )
                    {    mis.push_back(True);
                         mis_count++;
                         int q = qv[p1];
                         if ( q >= qtop ) 
                         {    bad = True;
                              break;    }    }
                    else mis.push_back(False);
                    if ( mis.isize() > min_len && mis[mis.isize()-min_len-1] )
                         mis_count--;
                    if ( mis.isize( ) >= min_len && mis_count <= max_diffs )
                         good = True;    }
               if ( bad || !good ) del[j] = True;    }
          EraseIf( places[bi], del );    }

     // cout << "destroying" << endl;
     Destroy(kmers_plus), Destroy(trashb);

     // Collate places.

     cout << Date( ) << ": pushing" << endl;
     ReportPeakMem( );
     vec< triple<int,int,int> > PLACES;
     for ( int bi = 0; bi < places.isize( ); bi++ )
          PLACES.append( places[bi] );
     ParallelUniqueSort(PLACES);

     // For now, require unique placement.

     vec<Bool> pdel( PLACES.size( ), False );
     for ( int64_t i = 0; i < PLACES.jsize( ); i++ )
     {    int64_t j;
          for ( j = i + 1; j < PLACES.jsize( ); j++ )
               if ( PLACES[j].first != PLACES[i].first ) break;
          if ( j - i > 1 )
          {    for ( int64_t k = i; k < j; k++ )
                    pdel[k] = True;    }
          i = j - 1;    }
     EraseIf( PLACES, pdel );

     // Report.

     Bool verbose = False;
     if (verbose)
     {    for ( int i = 0; i < PLACES.isize( ); i++ )
          {    cout << "read " << PLACES[i].first << " at "
                    << PLACES[i].second << "." << PLACES[i].third << endl;    }    }
     for ( int i = 0; i < PLACES.isize( ); i++ )
     {    IntVec x = { PLACES[i].second };
          (IntVec&) paths[ PLACES[i].first ] = x;
          paths[ PLACES[i].first ].setOffset( PLACES[i].third );    }
     cout << "\n" << Date( ) << ": done" << endl;    }


double CNIntegerFraction(const HyperBasevector& hb, const vec<vec<covcount>>& covs,
			 const double frac, const int min_edge_size) {
    
    size_t edge_count = 0;  // all edges larger than min_edge_size
    size_t good_count = 0;  // edges within frac of an int

    for ( size_t ss = 0; ss < covs.size(); ss++ ) // for each sample
	for ( int e = 0; e < hb.E(); e++ ) {
	    int len = hb.EdgeLength(e);
	    if (len >= min_edge_size && covs[ss][e].Def( ) ) {
		edge_count++;
		double frac_cov = covs[ss][e].Cov();
		int int_cov = frac_cov + 0.5;
		if ( abs(int_cov - frac_cov) <= frac)
		    good_count++;
	    }
	}
    return static_cast<double>(good_count) / edge_count;
}
