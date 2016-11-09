// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

// Rescue.  Find sequences in reads that start with an assembly kmer and end with 
// an assembly kmer, contain kmers that are not in the assembly, and which appear
// at least twice in the reads (appropriately counted).  The counting is such
// that duplicates (as determined by single read starts) don't count, nor do
// multiple reads from the same barcode, nor do reads from barcode zero.
//
// We ignore pairs for which both edges are less than K kmers long.
//
// Modified to look only for putative missing SNPs.

// Notes:
// 1. If we knew duplicates in advance, we could avoid processing them, thereby
//    reducing both time and memory.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "MainTools.h"
#include "ParallelVecUtilities.h"
#include "kmers/KmerRecord.h"
#include "math/Functions.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "10X/DfTools.h"
#include "10X/Rescue.h"

inline void GetPassIds( const basevector& u, const int L, const int j, 
     const int npasses, int& pass1, int& pass2 )
{    if ( npasses == 4 )
     {    pass1 = u[j];
          pass2 = 3 - u[j+L-1];    }
     else if ( npasses == 16 )
     {    pass1 = 4 * u[j] + u[j+1];
          pass2 = 4 * ( 3 - u[j+L-1] ) + ( 3 - u[j+L-2] );    }
     else if ( npasses == 64 )
     {    pass1 = 16 * u[j] + 4 * u[j+1] + u[j+2];
          pass2 = 16 * ( 3 - u[j+L-1] ) + 4 * ( 3 - u[j+L-2] )
               + ( 3 - u[j+L-3] );    }
     else if ( npasses == 256 )
     {    pass1 = 64 * u[j] + 16 * u[j+1] + 4 * u[j+2] + u[j+3];
          pass2 = 64 * ( 3 - u[j+L-1] ) + 16 * ( 3 - u[j+L-2] )
               + 4 * ( 3 - u[j+L-3] ) + ( 3 - u[j+L-4] );    }
     else
     {    FatalErr("npasses illegal");    }    }

template<int K, class T> void FindEm( const int64_t START,
     const vecbasevector& bases, vec<int64_t>& bid,
     const HyperBasevectorX& hb, const vec<int>& inv,
     const vec<char>& match, const vec<int64_t>& rstart, 
     const ReadPathVec& paths, vec<basevector>& patches )
{
     // Go through 64 passes.  Might be able to use 16, which is faster but
     // requires more memory.

     cout << Date( ) << ": searching reads" << endl;
     Bool verbose = ( START > 0 );
     const vecbasevector& A = hb.Edges( );
     double clock = WallClockTime( ), sclock = 0, bclock = 0, tclock = 0;
     double pclock = 0;
     // vec< pair<T,int> > hits;
     vec< quad<T,int,int,int> > hits;
     vec< triple< kmer<K>, T, int > > kmers;
     const int npasses = 64;
     for ( int pass = 0; pass < npasses; pass++ )
     {    cout << Date( ) << ": starting pass " << pass+1 << " of " 
               << npasses << endl;

          // Index kmers in the reads.

          pclock -= WallClockTime( );
          int64_t N = bases.size( );
          {    vec<int> counts( bases.size( ) + A.size( ), 0 );
               #pragma omp parallel for schedule(dynamic, 10000)
               for ( size_t i = 0; i < bases.size( ); i++ )
               {    if ( bid[i] == 0 ) continue;
                    const basevector& u = bases[i];
                    int64_t s = rstart[i];
                    for ( int j = 0; j <= u.isize( ) - K; j++ )
                    {    if ( match[ (s+j)/8 ] & ( 1 << ( (s+j) % 8 ) ) ) continue;
                         int pass1 = 0, pass2 = 0;
                         GetPassIds( u, K, j, npasses, pass1, pass2 );
                         if ( pass1 < pass || pass2 < pass ) continue;
                         if ( pass1 > pass && pass2 > pass ) continue;
                         counts[i]++;    }    }
               #pragma omp parallel for schedule(dynamic, 10000)
               for ( size_t i = 0; i < A.size( ); i++ )
               {    const basevector& u = A[i];
                    for ( int j = 0; j <= u.isize( ) - K; j++ )
                    {    int pass1 = 0, pass2 = 0;
                         GetPassIds( u, K, j, npasses, pass1, pass2 );
                         if ( pass1 < pass || pass2 < pass ) continue;
                         if ( pass1 > pass && pass2 > pass ) continue;
                         counts[ i + bases.size( ) ]++;    }    }
               vec<int64_t> starts( counts.size( ) + 1 );
               starts[0] = 0;
               for ( size_t i = 0; i < counts.size( ); i++ )
                    starts[i+1] = starts[i] + counts[i];
               kmers.resize( starts.back( ) );
               pclock += WallClockTime( );
               bclock -= WallClockTime( );
               #pragma omp parallel for schedule(dynamic, 10000)
               for ( size_t i = 0; i < counts.size( ); i++ )
               {    if ( i < N && bid[i] == 0 ) continue;
                    const basevector& u 
                         = ( i < bases.size( ) ? bases[i] : A[ i - bases.size( ) ] );
                    int count = 0;
                    int64_t s = -1;
                    if ( i < N ) s = rstart[i];
                    for ( int j = 0; j <= u.isize( ) - K; j++ )
                    {    if ( i < N )
                         {    if ( match[ (s+j)/8 ] & ( 1 << ( (s+j) % 8 ) ) )
                                   continue;    }
                         int pass1 = 0, pass2 = 0;
                         GetPassIds( u, K, j, npasses, pass1, pass2 );
                         if ( pass1 < pass || pass2 < pass ) continue;
                         if ( pass1 > pass && pass2 > pass ) continue;
                         int64_t r = starts[i] + count++;
                         kmers[r].first.SetToSubOf( u, j );
                         kmers[r].second = i, kmers[r].third = j;    }    }    }
          bclock += WallClockTime( );

          // Sort the kmers.  This is presumably the memory high-water mark, because
          // it's not an in-place sort.  People have tried to write in-place 
          // parallel sorts, but they're slower (google the topic).  However it's
          // not impossible that the kmer sorting plus searching could somehow
          // be done in chunks, thus reducing memory.

          tclock -= WallClockTime( );
          ParallelSort(kmers);    
          tclock += WallClockTime( );

          // Look for kmers that are in the edges but not the reads.

          sclock -= WallClockTime( );
          const int batches = 100;
          vec<int64_t> bstarts( batches+1, 0 );
          for ( int i = 0; i <= batches; i++ )
               bstarts[i] = ( i * kmers.jsize( ) ) / batches;
          for ( int i = 0; i < batches; i++ )
          {    if ( bstarts[i] > 0
                    && kmers[ bstarts[i] ].first == kmers[ bstarts[i] - 1 ].first )
               {    bstarts[i]--;    }    }
          #pragma omp parallel for schedule(dynamic, 1)
          for ( int b = 0; b < batches; b++ )
          {    // vec< pair<T,int> > hitsb;
               vec< quad<T,int,int,int> > hitsb;
               for ( int64_t i = bstarts[b]; i < bstarts[b+1]; i++ )
               {    int64_t j, m;
                    for ( j = i + 1; j < kmers.jsize( ); j++ )
                         if ( kmers[j].first != kmers[i].first ) break;
                    for ( m = i; m < j; m++ ) if ( kmers[m].second >= N ) break;
                    if ( m < j )
                    {    for ( int64_t k1 = i; k1 < m; k1++ )
                         {    int64_t id = kmers[k1].second;
                              int rpos = kmers[k1].third;
                              for ( int64_t k2 = m; k2 < j; k2++ )
                              {    int e = kmers[k2].second - N;
                                   int epos = kmers[k2].third;
                                   hitsb.push( id, rpos, e, epos );    }    }    }
                    i = j - 1;    }
               #pragma omp critical
               {    hits.append(hitsb);    }    }
          kmers.clear( );
          sclock += WallClockTime( );    }
     cout << Date( ) << ": " << pclock << " seconds used in setup" << endl;
     cout << Date( ) << ": " << bclock << " seconds used in build" << endl;
     cout << Date( ) << ": " << tclock << " seconds used in sort" << endl;
     cout << Date( ) << ": " << sclock << " seconds used in search" << endl;

     // Collate results.

     cout << Date( ) << ": " << ToStringAddCommas( hits.size( ) ) 
          << " hits before sorting" << endl;
     ParallelSort(hits);
     // ParallelUniqueSort(hits);
     cout << Date( ) << ": " << ToStringAddCommas( hits.size( ) ) 
          << " hits after sorting" << endl;
     int64_t nhits = 0;
     vec< triple< pair<int,int>, pair<int,int>, basevector > >
          bridges; // { ( (e1,e2), (p1,p2), seq ) }

     const int batches = 100;
     vec<int64_t> hstarts( batches+1, 0 );
     for ( int i = 0; i <= batches; i++ )
          hstarts[i] = ( i * hits.jsize( ) ) / batches;
     for ( int i = 0; i < batches; i++ )
     {    if ( hstarts[i] > 0
               && hits[ hstarts[i] ].first == hits[ hstarts[i] - 1 ].first )
          {    hstarts[i]--;    }    }
     vec< pair<int64_t,int> > bridge_locs;
     int ngroups = 0;
     #pragma omp parallel for
     for ( int bi = 0; bi < batches; bi++ )
     {    vec< triple< pair<int,int>, pair<int,int>, basevector > > bridgesb;
          vec< pair<int64_t,int> > bridge_locsb;
          for ( int64_t i = hstarts[bi]; i < hstarts[bi+1]; i++ )
          {    int64_t j, id = hits[i].first;
               for ( j = i + 1; j < hits.jsize( ); j++ )
                    if ( hits[j].first != id ) break;
               vec<int> es;
               for ( int64_t k = i; k < j; k++ )
               {    int e = hits[k].third;
                    es.push_back(e);    }
               UniqueSort(es);
               const ReadPath& p = paths[id];
               {    
                    // cout << "\nread " << id << " = " << p.getOffset( ) << " : "
                    //      << printSeq(p) << " ==> {" << printSeq(es) << "}\n";    
                    const basevector& r = bases[id];
                    vec< triple<int,int,int> > M; // M[rpos] = (e, epos, status)
                    M.resize( r.isize( ) - K + 1, make_triple( 0, 0, 0 ) );
                    for ( int64_t k = i; k < j; k++ )
                    {    int rpos = hits[k].second;
                         M[rpos] = 
                              make_triple( hits[k].third, hits[k].fourth, 1 );    }
                    if ( p.size( ) > 0 )
                    {    int apos, rpos;
                         if ( p.getOffset( ) < 0 )
                         {    rpos = -p.getOffset( );
                              apos = 0;    }
                         else
                         {    rpos = 0;
                              apos = p.getOffset( );    }
                         for ( int j = 0; j < (int) p.size( ); j++ )
                         {    int e = p[j];
                              if ( j > 0 ) apos = 0;
                              while( apos <= A[e].isize( ) - K )
                              {    if ( rpos > r.isize( ) - K ) break;
                                   if ( M[rpos].third == 0 )
                                   {    Bool mismatch = False;
                                        // inefficient:
                                        for ( int u = 0; u < K; u++ )
                                        {    if ( A[e][apos+u] != r[rpos+u] )
                                             {    mismatch = True;
                                                  break;    }    }
                                        if ( !mismatch )
                                        {    M[rpos] = make_triple( 
                                                  e, apos, 2 );    }    }
                                   apos++; rpos++;    }    }    }
                    for ( int x = 0; x < M.isize( ); x++ )
                    {    int y;
                         for ( y = x + 1; y < M.isize( ); y++ )
                         {    if ( M[y].third != M[x].third ) break;
                              if ( M[y].first != M[x].first ) break;    }
                         int rpos1 = x, rpos2 = y - 1;
     
                         // int nb = bridges.size( );
                         if ( x > 0 && M[x-1].third == 2 && M[x].third == 1
                              && hb.ToRight( M[x-1].first ) 
                              != hb.ToLeft( M[x].first ) )
                         {    bridgesb.push( make_pair( M[x-1].first, M[x].first ),
                                   make_pair( M[x-1].second, M[x].second ),
                                   basevector( ) );
                              bridge_locsb.push( id, x-1 );    }
                         else if ( x > 0 && M[x-1].third == 1 && M[x].third == 2
                              && hb.ToRight( M[x-1].first ) 
                              != hb.ToLeft( M[x].first ) )
                         {    bridgesb.push( make_pair( M[x-1].first, M[x].first ),
                                   make_pair( M[x-1].second, M[x].second ),
                                   basevector( ) );
                              bridge_locsb.push( id, x-1 );    }
                         else if ( x > 0 && M[x-1].third == 0 && M[x].third > 0 )
                         {    int w;
                              for ( w = x - 1; w >= 0; w-- )
                                   if ( M[w].third > 0 ) break;
                              if ( w >= 0 )
                              {    basevector b;
                                   for ( int r = w + 1; r < x; r++ )
                                        b.push_back( bases[id][r] );
                                   bridgesb.push( 
                                        make_pair( M[w].first, M[x].first ),
                                        make_pair( M[w].second, M[x].second ), b );
                                   bridge_locsb.push( id, w );     }     }

                         x = y - 1;    }    }
               nhits++;
               i = j - 1;    }
          #pragma omp critical
          {    bridges.append(bridgesb);
               bridge_locs.append(bridge_locsb);    }    }

     cout << endl;

     // Add in reverse complement bridges.  Note super-inefficient 
     // determination of rb.

     cout << Date( ) << ": adding rc bridges" << endl;
     int64_t nb = bridges.size( );
     for ( int64_t i = 0; i < nb; i++ )
     {    int e1 = bridges[i].first.first, e2 = bridges[i].first.second;
          int p1 = bridges[i].second.first, p2 = bridges[i].second.second;
          basevector b = bridges[i].third;
          int rp1 = hb.Bases(e1) - p1 - K, rp2 = hb.Bases(e2) - p2 - K;
          basevector x( hb.O(e1), 0, p1 + 1 );
          x.append(b);
          x.append( basevector( hb.O(e2), p2, hb.Bases(e2) - p2 ) );
          x.ReverseComplement( );
          basevector rb( x, rp2 + 1, b.size( ) );
          bridges.push( make_pair( inv[e2], inv[e1] ), make_pair( rp2, rp1 ), rb );
          bridge_locs.push( 
               bridge_locs[i].first, - 1 - bridge_locs[i].second );    }

     // Analyze bridges and generate patches.

     cout << Date( ) << ": sorting bridges, mem usage = " << MemUsageGBString( )
          << endl;
     SortSync( bridges, bridge_locs );
     cout << Date( ) << ": analyzing bridges, mem usage = " << MemUsageGBString( )
          << endl;
     for ( int64_t i = 0; i < bridges.jsize( ); i++ )
     {    int64_t j;
          for ( j = i + 1; j < bridges.jsize( ); j++ )
          {    if ( bridges[j].first != bridges[i].first ) break;
               if ( bridges[j].second != bridges[i].second ) break;
               if ( bridges[j].third != bridges[i].third ) break;    }

          // Don't use short edges.

          if ( hb.Kmers( bridges[i].first.first ) < K
               && hb.Kmers( bridges[i].first.second ) < K )
          {    i = j - 1;
               continue;    }

          // Disallow homopolymers.  (PROBABLY NOT NEEDED NOW.)

          if ( bridges[i].first.first == bridges[i].first.second )
          {    i = j - 1;
               continue;    }

          // Require at least two reads.

          const int MIN_SUPPORT = 2;
          if ( j - i >= MIN_SUPPORT )
          {    
               // Find duplicates, but just checking one read.  Also check barcodes,
               // which must be different.

               vec<Bool> dup( j-i, False );
               {    vec<int> ids( j-i, vec<int>::IDENTITY ), pos;
                    for ( int64_t k = i; k < j; k++ )
                         pos.push_back( bridge_locs[k].second );
                    SortSync( pos, ids );
                    for ( int64_t k = 1; k < j-i; k++ )
                         if ( pos[k] == pos[k-1] ) dup[ ids[k] ] = True;    }
               vec<int> ids( j-i, vec<int>::IDENTITY ), bar;
               for ( int64_t k = i; k < j; k++ )
                    bar.push_back( bid[ bridge_locs[k].first ] );
               SortSync( bar, ids );
               for ( int64_t k = 1; k < j-i; k++ )
                    if ( bar[k] == bar[k-1] ) dup[ ids[k] ] = True;

               // Show groups.

               if ( j - i - Sum(dup) >= MIN_SUPPORT )
               {    
                    // Require that event "looks" like a SNP.
          
                    int e1 = bridges[i].first.first, e2 = bridges[i].first.second;
                    Bool looks_like = False;
                    if ( hb.ToLeft(e1) == hb.ToLeft(e2)
                         || hb.ToRight(e1) == hb.ToRight(e2) )
                    {
                         for ( int64_t k = i; k < j; k++ )
                         {    if ( dup[k-i] ) continue;
                              int p1 = bridges[k].second.first;
                              int p2 = bridges[k].second.second;
                              basevector b = bridges[k].third;
                              basevector alt( hb.O(e1), 0, p1+1 );
                              alt.append(b);
                              alt.append( basevector( 
                                   hb.O(e2), p2, hb.Bases(e2) - p2 ) );
                              if ( alt.isize( ) == hb.Bases(e2) )
                              {    int diffs = 0;
                                   for ( int l = 0; l < hb.Bases(e2); l++ )
                                   {    if ( alt[l] != hb.O(e2)[l] ) diffs++;    }
                                   if ( diffs == 1 ) 
                                        looks_like = True;    }    }    }
                    if ( !looks_like )
                    {    i = j - 1;
                         continue;    }

                    // Save group.

                    if (verbose) cout << "\nbridge group " << ++ngroups << ":\n";
                    for ( int64_t k = i; k < j; k++ )
                    {    if ( dup[k-i] ) continue;
                         int e1 = bridges[k].first.first;
                         int e2 = bridges[k].first.second;
                         int p1 = bridges[k].second.first;
                         int p2 = bridges[k].second.second;
                         basevector b = bridges[k].third;
                         if (verbose)
                         {    cout << "e1 = " << e1 << ", e2 = " << e2
                                   << ", p1 = " << p1 << ", p2 = " << p2
                                   << ", seq = " << b.ToString( );
                              if (looks_like) cout << " [looks like SNP]";
                              cout << endl;    }
                         int start = Max( 0, p1 - K ), stop = p1;
                         basevector patch( hb.O(e1), start, stop - start + 1 );
                         patch.append(b);
                         patch.append( basevector( hb.O(e2), p2, K ) );
                         patches.push_back(patch);    
                         patch.ReverseComplement( );
                         patches.push_back(patch); // unnecessary?
                         if (verbose) 
                              cout << "patch = " << patch.ToString( ) << endl;
                         break;    }
                    if (verbose)
                    {    for ( int64_t k = i; k < j; k++ )
                         {    if ( dup[k-i] ) continue;
                              int e1 = bridges[k].first.first;
                              int e2 = bridges[k].first.second;
                              int p1 = bridges[k].second.first;
                              int p2 = bridges[k].second.second;
                              basevector b = bridges[k].third;
                              cout << "id = " << START + bridge_locs[k].first << "."
                                   << bridge_locs[k].second 
                                   << endl;    }    }    }    }
          i = j - 1;    }

     // Done.

     cout << "\n" << Date( ) << ": " << ToStringAddCommas(nhits) 
          << " reads are hit" << endl;
     cout << Date( ) << ": done, time used = " << TimeSince(clock)
          << ", peak mem used = " << setiosflags(ios::fixed) << setprecision(1)
          << PeakMemUsageBytes( ) / 1000000000.0 << resetiosflags(ios::fixed)
          << " GB" << endl;     }


void Rescue(
     const int64_t START,
     const vecbasevector& bases,
     const vec<int32_t>& bc,
     const vec<DataSet>& datasets,
     const HyperBasevectorX& hb,
     const vec<int>& inv,
     const ReadPathVec& paths,
     vec<basevector>& patches )
{
     // Start.

     double aclock = WallClockTime( );
     cout << Date( ) << ": total edges = " << ToStringAddCommas( hb.E( ) ) << endl;
     const int K = 48;
     int64_t N = bases.size( );
     Bool verbose = ( START > 0 );

     // Create a vector of integers, one for each read, such that "having two" 
     // nonzero elements is enough.  The value of an "unbarcoded" 10X read is ZERO.

     cout << Date( ) << ": creating bid" << endl;
     vec<int64_t> bid( paths.size( ) );
     #pragma omp parallel for
     for ( int64_t id = 0; id < (int64_t) paths.size( ); id++ )
     {    int di;
          for ( di = 0; di < datasets.isize( ); di++ )
               if ( START+id < datasets[di].start ) break;
          const ReadDataType& dtype = datasets[di-1].dt;
          if ( dtype == ReadDataType::BAR_10X )
               bid[id] = (int64_t) paths.size( ) + bc[START+id] + 1;
          else if ( dtype == ReadDataType::UNBAR_10X ) bid[id] = 0;
          else if ( dtype == ReadDataType::PCR_FREE ) bid[id] = id + 1;

          // Probably not what we want:

          else if ( dtype == ReadDataType::PCR ) bid[id] = id + 1;    }

     /*
     for ( int64_t id = 0; id < (int64_t) paths.size( ); id++ )
          bid[id] = bc[ START + id ];
     */

     // Mark kmer positions in the reads that match their paths.

     cout << Date( ) << ": finding match positions" << endl;
     int64_t nkmers = 0;
     vec<int64_t> rstart( bases.size( ) );
     for ( int64_t id = 0; id < N; id++ )
     {    rstart[id] = nkmers;
          int n = Max( 0, bases[id].isize( ) - (K-1) );
          n = 8 * ((n+7)/8); // add fill bits so two reads don't touch one byte
          nkmers += n;    }
     vec<char> match( nkmers/8, 0 );
     const int batch = 100000;
     #pragma omp parallel for schedule(dynamic, 1)
     for ( int64_t bi = 0; bi < N; bi += batch )
     {    vec<int> bads;
          for ( int64_t id = bi; id < Min( bi + batch, N ); id++ )
          {    if ( bid[id] == 0 ) continue;
               const basevector& r = bases[id];
               const ReadPath& p = paths[id];
               if ( p.size( ) == 0 ) continue;
               int apos, rpos;
               if ( p.getOffset( ) < 0 )
               {    rpos = -p.getOffset( );
                    apos = 0;    }
               else
               {    rpos = 0;
                    apos = p.getOffset( );    }
               bads.clear( );
               for ( int j = 0; j < rpos; j++ )
                    bads.push_back(j);
               for ( int j = 0; j < (int) p.size( ); j++ )
               {    int e = p[j];
                    if ( j > 0 ) apos = K - 1;
                    while( apos < hb.O(e).isize( ) )
                    {    if ( rpos == r.isize( ) ) break;
                         if ( hb.O(e)[apos] != r[rpos] ) bads.push_back(rpos);
                         apos++; rpos++;    }    }
               for ( int j = rpos; j < r.isize( ); j++ )
                    bads.push_back(j);
               int64_t pos = rstart[id];
               int b = 0;
               for ( int j = 0; j <= r.isize( ) - K; j++ )
               {    while ( b < bads.isize( ) && j > bads[b] ) b++;
                    if ( b == bads.isize( ) || bads[b] - j >= K )
                         match[pos/8] |= ( 1 << (pos%8) );
                    pos++;    }    }    }
     int64_t total = 0, total_match = 0;
     for ( int64_t id = 0; id < N; id++ )
          if ( bid[id] > 0 ) total += Max( 0, bases[id].isize( ) - (K-1) );
     for ( int64_t i = 0; i < match.jsize( ); i++ )
     for ( int j = 0; j < 8; j++ )
          if ( match[i] & (1 << j) ) total_match++;
     cout << Date( ) << ": " << PERCENT_RATIO( 3, total_match, total )
          << " of kmers match paths" << endl;

     // Find hits.

     if ( bases.size( ) < 2147483648 ) 
     {    FindEm<K,int32_t>( 
               START, bases, bid, hb, inv, match, rstart, paths, patches );    }
     else 
     {    FindEm<K,int64_t>( 
               START, bases, bid, hb, inv, match, rstart, paths, patches );    }

     // Done.

     cout << Date( ) << ": total time used = " << TimeSince(aclock) << endl;    }

