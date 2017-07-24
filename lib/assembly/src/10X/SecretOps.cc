// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"
#include "Intvector.h"
#include "ParallelVecUtilities.h"
#include "VecUtilities.h"
#include "feudal/PQVec.h"
#include "graph/Digraph.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "random/Random.h"
#include "10X/Closer.h"
#include "10X/Heuristics.h"
#include "10X/SecretOps.h"
#include "system/SortInPlace.h"
#include "10X/DfTools.h"
#include "10X/paths/ReadPathVecX.h"

template< class VB, class VQ> void MarkBads( const HyperBasevectorX& hb,
     VB& bases, VQ& quals, ReadPathVecX& paths, vec<Bool>& bad )
{
     // Heuristics.

     const int MAX_BAD = 5;
     const int MIN_HQ = 30;
     const int MAX_BAD_SUM = 150;

     // Proceed.

     bad.resize_and_set( bases.size( )/2, False );
     const int64_t batch = 10000;
     int64_t N = bases.size( );
     auto quals_clone = quals;
     #pragma omp parallel for schedule(dynamic,1) firstprivate(quals_clone)
     for ( int64_t bi = 0; bi < N; bi += batch )
     {    qualvector q;
          basevector m;
          vec<int> x;
          for ( int64_t id = bi; id < Min( bi + batch, N ); id++ )
          {    ReadPath p; paths.unzip(p,hb,id);
               if ( p.size( ) == 0 ) continue;
               int hqm = 0;
               const basevector& b = bases[id];
               quals_clone[id].unpack(&q);
               x.clear( );
               for ( int l = 0; l < (int) p.size( ); l++ )
                    x.push_back( p[l] );
               m = hb.Cat(x);
               int badsum = 0;
               for ( int l = 0; l < b.isize( ); l++ )
               {    int epos = p.getOffset( ) + l;
                    if ( epos < 0 || epos >= m.isize( ) ) continue;
                    if ( b[l] != m[epos] ) badsum += q[l];    }
               if ( badsum > MAX_BAD_SUM ) bad[id/2] = True;    }    }
     cout << PERCENT_RATIO( 3, Sum(bad), (int64_t) bases.size( )/2 )
          << " of pairs marked bad" << endl;    }

template void MarkBads( const HyperBasevectorX&, VirtualMasterVec<basevector>&,
     VirtualMasterVec<PQVec>&, ReadPathVecX&, vec<Bool>& );

template void MarkBads( const HyperBasevectorX&, MasterVec<basevector>&,
     MasterVec<PQVec>&, ReadPathVecX&, vec<Bool>& );

template void MarkBads( const HyperBasevectorX&, MasterVec<basevector>&,
     VirtualMasterVec<PQVec>&, ReadPathVecX&, vec<Bool>& );

template< class VB, class VQ, class VP > void MarkBads( const HyperBasevectorX& hb,
     VB& bases, VQ& quals, VP& paths, vec<Bool>& bad )
{
     // Heuristics.

     const int MAX_BAD = 5;
     const int MIN_HQ = 30;
     const int MAX_BAD_SUM = 150;

     // Proceed.

     bad.resize_and_set( bases.size( )/2, False );
     const int64_t batch = 10000;
     int64_t N = bases.size( );
     auto quals_clone = quals;
     #pragma omp parallel for schedule(dynamic,1) firstprivate(quals_clone)
     for ( int64_t bi = 0; bi < N; bi += batch )
     {    qualvector q;
          basevector m;
          vec<int> x;
          for ( int64_t id = bi; id < Min( bi + batch, N ); id++ )
          {    const ReadPath& p = paths[id];
               if ( p.size( ) == 0 ) continue;
               int hqm = 0;
               const basevector& b = bases[id];
               quals_clone[id].unpack(&q);
               x.clear( );
               for ( int l = 0; l < (int) p.size( ); l++ )
                    x.push_back( p[l] );
               m = hb.Cat(x);
               int badsum = 0;
               for ( int l = 0; l < b.isize( ); l++ )
               {    int epos = p.getOffset( ) + l;
                    if ( epos < 0 || epos >= m.isize( ) ) continue;
                    if ( b[l] != m[epos] ) badsum += q[l];    }
               if ( badsum > MAX_BAD_SUM ) bad[id/2] = True;    }    }
     cout << PERCENT_RATIO( 3, Sum(bad), (int64_t) bases.size( )/2 )
          << " of pairs marked bad" << endl;    }

template void MarkBads( const HyperBasevectorX&, VirtualMasterVec<basevector>&,
     VirtualMasterVec<PQVec>&, VirtualMasterVec<ReadPath>&, vec<Bool>& );

template void MarkBads( const HyperBasevectorX&, MasterVec<basevector>&,
     MasterVec<PQVec>&, MasterVec<ReadPath>&, vec<Bool>& );

template void MarkBads( const HyperBasevectorX&, MasterVec<basevector>&,
     VirtualMasterVec<PQVec>&, MasterVec<ReadPath>&, vec<Bool>& );

// IntSpaceExts.  Extend to right in integer space.

template< class VP, class VPI > void IntSpaceExts( const int e,
     const HyperBasevectorX& hb, const vec<int>& kmers, const vec<int>& inv,
     VP& paths, VPI& paths_index, const vec<int32_t>& bc,
     const vec<Bool>& dup, const vec<Bool>& bad, const Bool GLOBAL,
     const Bool PRINT_RIGHT_EXTS, const int mu, vec<vec<int>>& ext,
     const Bool verbose, IntSpaceExtsWorkspace& work )
{
     if (verbose) PRINT(e);
     IntSpaceExtsWorkspace& W = work;
     ext.clear( );
     W.ext2.clear( );
     W.epo.clear( );
     for ( int pass = 1; pass <= 2; pass++ )
     {    const int f = ( pass == 1 ? e : inv[e] );
          for ( int i = 0; i < (int) paths_index[f].size( ); i++ )
          {    int64_t id = paths_index[f][i];
               if ( dup[id/2] ) continue;
               if ( bad[id/2] ) continue;
               const ReadPath& p = paths[id];
               int n = p.size( );
               W.y.clear( );
               if ( pass == 1 )
               {    for ( int j = 0; j < n; j++ )
                         W.y.push_back( p[j] );    }
               else
               {    for ( int j = n - 1; j >= 0; j-- )
                         W.y.push_back( inv[p[j]] );    }
               for ( int j = 0; j < n - 1; j++ )
               {    if ( W.y[j] != e ) continue;
                    W.z.clear( );
                    for ( int k = j; k < n; k++ )
                         W.z.push_back( W.y[k] );
                    W.epo.push( bc[id], W.z );    }    }    }
                    // ext.push_back(z);    }    }    }
     UniqueSort(W.epo);

     W.edel.resize_and_set( W.epo.size( ), False );
     for ( int i = 1; i < W.epo.isize( ); i++ )
     {    if ( W.epo[i].first != W.epo[i-1].first ) continue;
          vec<int> &a1 = W.epo[i-1].second, &a2 = W.epo[i].second;
          if ( a2.Contains(a1) ) W.edel[i-1] = True;
          if ( a1.Contains(a2) ) W.edel[i] = True;

          // Suppose the two reads overlap.  Join them.

          vec<int> offsets;
          vec< triple<int,int,int> > m;
          for ( int j = 0; j < a1.isize( ); j++ )
               m.push( a1[j], 1, j );
          for ( int j = 0; j < a2.isize( ); j++ )
               m.push( a2[j], 2, j );
          Sort(m);
          for ( int j = 1; j < m.isize( ); j++ )
          {    if ( m[j].first == m[j-1].first && m[j].second != m[j-1].second )
                    offsets.push_back( m[j-1].third - m[j].third );    }
          UniqueSort(offsets);
          if ( offsets.solo( ) && offsets[0] >= 0 )
          {    vec<int> join;
               int o = offsets[0];
               for ( int l = 0; l < o; l++ )
                    join.push_back( a1[l] );
               join.append(a2);
               W.epo[i-1].second = join;
               W.edel[i] = True;    }    }

     EraseIf( W.epo, W.edel );

     for ( int i = 0; i < W.epo.isize( ); i++ )
          ext.push_back( W.epo[i].second );
     Sort(ext);

     for ( int i = 0; i < ext.isize( ) - 1; i++ )
     {    int k;
          for ( k = 1; k < Min( ext[i].isize( ), ext[i+1].isize( ) ); k++ )
          {    if ( ext[i][k] != ext[i+1][k] ) break;    }
          vec<int> x;
          for ( int j = 0; j < k; j++ )
               x.push_back( ext[i][j] );
          W.ext2.push_back(x);    }
     UniqueSort(W.ext2);
     W.to_delete2.resize_and_set( W.ext2.size( ), False );
     for ( int i = 0; i < W.ext2.isize( ) - 1; i++ )
          if ( Subset( W.ext2[i], W.ext2[i+1] ) ) W.to_delete2[i] = True;
     EraseIf( W.ext2, W.to_delete2 );
     if (verbose)
     {    cout << "\next2:\n";
          for ( int j = 0; j < W.ext2.isize( ); j++ )
               cout << "(" << j+1 << ") " << printSeq( W.ext2[j] ) << endl;
          cout << endl;    }

     // Find inverted partner paths.  Produces a vast and duplicated list of stuff.

     W.pp.clear( );
     for ( int i = 0; i < (int) paths_index[e].size( ); i++ )
     {    int64_t id1 = paths_index[e][i];
          if ( dup[id1/2] ) continue;
          int64_t id2 = ( id1 % 2 == 0 ? id1 + 1 : id1 - 1 );
          const ReadPath& p = paths[id2];
          if ( p.size( ) == 0 ) continue;
          int n = p.size( );
          W.y.clear( );
          for ( int j = n - 1; j >= 0; j-- )
               W.y.push_back( inv[p[j]] );
          W.pp.push_back(W.y);    }

     // Remove edges appearing only once in partner paths.

     vec<int> all;
     for ( int i = 0; i < W.pp.isize( ); i++ )
          all.append( W.pp[i] );
     Sort(all);
     for ( int j = 0; j < all.isize( ); j++ )
     {    int k = all.NextDiff(j);
          if ( k - j == 1 )
          {    for ( int i = 0; i < W.pp.isize( ); i++ )
               {    for ( int m = 0; m < W.pp[i].isize( ); m++ )
                    {    if ( W.pp[i][m] == all[j] )
                         {    W.pp[i].resize(m);
                              goto pp_trimmed;    }    }    }
               pp_trimmed: continue;    }
          j = k - 1;    }
     UniqueSort(W.pp);
     vec<Bool> pdel( W.pp.size( ), False );
     for ( int i = 0; i < W.pp.isize( ); i++ )
          if ( W.pp[i].empty( ) ) pdel[i] = True;
     EraseIf( W.pp, pdel );

     // XXX:
     /*
     cout << "\next2:\n";
     for ( int i = 0; i < W.ext2.isize( ); i++ )
          cout << "[" << i+1 << "] " << printSeq(W.ext2[i]) << endl;
     cout << "\npp:\n";
     for ( int i = 0; i < W.pp.isize( ); i++ )
          cout << "[" << i+1 << "] " << printSeq(W.pp[i]) << endl;
     */

     // Add in partner paths.  Stupendously inefficient.

     while(1)
     {    Bool progress = False;
          W.ext3.clear( );
          for ( int i = 0; i < W.ext2.isize( ); i++ )
          {    W.x = W.ext2[i];
               Bool extended = False;
               for ( int j = 0; j < W.pp.isize( ); j++ )
               {    W.y = W.pp[j];
                    for ( int k = 0; k < W.x.isize( ); k++ )
                    {    if ( W.x[k] != W.y[0] ) continue;
                         int ext = k + W.y.isize( ) - W.x.isize( );
                         if ( ext <= 0 ) continue;
                         Bool mismatch = False;
                         for ( int l = 1; l < W.y.isize( ); l++ )
                         {    if ( k+l >= W.x.isize( ) ) break;
                              if ( W.x[k+l] != W.y[l] )
                              {    mismatch = True;
                                   break;    }    }
                         if (mismatch) continue;
                         Bool present = False;
                         for ( int j = W.y.isize( ) - ext; j < W.y.isize( ); j++ )
                              if ( Member( W.x, W.y[j] ) ) present = True;
                         if (present) continue;
                         extended = True;
                         W.n = W.x;
                         for ( int j = W.y.isize( ) - ext; j < W.y.isize( ); j++ )
                              W.n.push_back( W.y[j] );
                         progress = True;
                         W.ext3.push_back(W.n);    }    }
               if ( !extended ) W.ext3.push_back( W.ext2[i] );    }
          UniqueSort(W.ext3);
          W.to_delete3.resize_and_set( W.ext3.size( ), False );
          for ( int i = 0; i < W.ext3.isize( ) - 1; i++ )
               if ( Subset( W.ext3[i], W.ext3[i+1] ) ) W.to_delete3[i] = True;
          EraseIf( W.ext3, W.to_delete3 );

          // Advance along simple bubble chains.

          if ( !progress )
          {    for ( int i = 0; i < W.ext3.isize( ); i++ )
               {    vec<int>& x = W.ext3[i];
                    int e = x.back( );
                    int v = hb.ToLeft(e), w = hb.ToRight(e);
                    if ( hb.From(w).size( ) == 1 && hb.To(w).size( ) == 2
                         && hb.To(w)[0] == hb.To(w)[1] )
                    {    x.push_back( hb.IFrom( w, 0 ) );    }
                    while(1)
                    {    int e = x.back( );
                         int w = hb.ToRight(e);
                         if ( hb.From(w).size( ) != 2 ) break;
                         if ( hb.To(w).size( ) != 1 ) break;
                         if ( hb.From(w)[0] != hb.From(w)[1] ) break;
                         int f1 = hb.IFrom(w,0), f2 = hb.IFrom(w,1);
                         int y = hb.From(w)[0];
                         if ( hb.To(y).size( ) != 2 ) break;
                         if ( hb.From(y).size( ) != 1 ) break;
                         int g = hb.IFrom( y, 0 );
                         if ( x.size( ) < 2 ) break;
                         vec<int> count( 2, 0 );
                         int d = x[ x.isize( ) - 2 ];
                         for ( int j = 0; j < (int) paths_index[d].size( ); j++ )
                         {    int64_t id = paths_index[d][j];
                              if ( dup[id/2] || bad[id/2] ) continue;
                              const ReadPath& p = paths[id];
                              for ( int l = 0; l < (int) p.size( ) - 2; l++ )
                              {    if ( p[l] == d && p[l+1] == e && p[l+2] == f1 )
                                        count[0]++;
                                   if ( p[l] == d && p[l+1] == e && p[l+2] == f2 )
                                        count[1]++;    }    }
                         int rd = inv[d];
                         for ( int j = 0; j < (int) paths_index[rd].size( ); j++ )
                         {    int64_t id = paths_index[rd][j];
                              if ( dup[id/2] || bad[id/2] ) continue;
                              const ReadPath& p = paths[id];
                              vec<int> x;
                              for ( int l = (int) p.size( ) - 1; l >= 0; l-- )
                                   x.push_back( inv[ p[l] ] );
                              for ( int l = 0; l < x.isize( ) - 2; l++ )
                              {    if ( x[l] == d && x[l+1] == e && x[l+2] == f1 )
                                        count[0]++;
                                   if ( x[l] == d && x[l+1] == e && x[l+2] == f2 )
                                        count[1]++;    }    }
                         if ( count[0] > 0 && count[1] > 0 ) break;
                         if ( count[0] < 2 && count[1] < 2 ) break;
                         if ( count[0] > 0 ) x.push_back( f1, g );
                         else x.push_back( f2, g );    }    }
               UniqueSort(W.ext3);    }

          // If advancing one edge will get us to a long edge, do that.

          if ( !progress )
          {    for ( int i = 0; i < W.ext3.isize( ); i++ )
               {    vec<int>& x = W.ext3[i];
                    int M = 0;
                    for ( int j = 1; j < x.isize( ); j++ )
                         M = Max( M, hb.Kmers( x[j] ) );
                    if ( M >= MIN_LEN ) continue;
                    int w = hb.ToRight( x.back( ) );
                    if ( hb.From(w).size( ) == 1 )
                    {    int f = hb.IFrom( w, 0 );
                         if ( hb.Kmers(f) >= MIN_LEN ) x.push_back(f);    }    }    }

          if ( PRINT_RIGHT_EXTS && !progress )
          {
               // Print.

               cout << "ext:\n";
               for ( int j = 0; j < W.ext3.isize( ); j++ )
               {    vec<int> x = W.ext3[j];
                    if ( mu > 0 )
                    {    x.ReverseMe( );
                         for ( int l = 0; l < x.isize( ); l++ )
                              x[l] = inv[ x[l] ];    }
                    int T = 0, M = 0;
                    for ( int l = 0; l < x.isize( ); l++ )
                    {    int n = kmers[ x[l] ];
                         T += n;
                         if ( mu == 0 && l == 0 ) continue;
                         if ( mu > 0 && l == x.isize( ) - 1 )
                              continue;
                         M = Max( M, n );    }
                    cout << "(" << j+1 << ") " << printSeq(x)
                         << ", total = " << T << ", max = " << M << endl;
                    /*
                    if (verbose)
                    {    basevector b = hb.Cat(x);
                         b.Print( cout, "ext_" + ToString(mu+1)
                              + "_" + ToString(j+1) );    }
                    */
                         }
               cout << endl;    }
          W.ext2 = W.ext3;
          if ( !progress ) break;    }
     ext = W.ext2;    }

template void IntSpaceExts( const int e, const HyperBasevectorX& hb,
     const vec<int>& kmers, const vec<int>& inv, VirtualMasterVec<ReadPath>&,
     VirtualMasterVec<ULongVec>&, const vec<int32_t>&,
     const vec<Bool>& dup, const vec<Bool>& bad,
     const Bool GLOBAL, const Bool PRINT_RIGHT_EXTS,
     const int mu, vec<vec<int>>& ext, const Bool verbose, IntSpaceExtsWorkspace& );

template void IntSpaceExts( const int e, const HyperBasevectorX& hb,
     const vec<int>& kmers, const vec<int>& inv, MasterVec<ReadPath>&,
     MasterVec<ULongVec>&, const vec<int32_t>&, const vec<Bool>& dup,
     const vec<Bool>& bad, const Bool GLOBAL, const Bool PRINT_RIGHT_EXTS,
     const int mu, vec<vec<int>>& ext, const Bool verbose, IntSpaceExtsWorkspace& );

// Note for MarkDups.  An "artifactual duplicate" is a read that appears to be
// part of a duplicate pair and which is base by base and quality score by quality
// score identical to another read.  If present in large numbers, they would
// presumably have to arise from an informatic mixup.

template< class VB, class VQ, class VP >
void MarkDups( VB& bases, VQ& quals, VP& paths, const vec<int32_t>& bc,
     vec<Bool>& dup, double& interdup_rate, const Bool verbose )
{    
     cout << Date( ) << ": MarkDups, version 1" << endl;

     // Heuristics.

     const int BHEAD = 5;

     // Start.

     double clock = WallClockTime( );
     dup.resize_and_set( bases.size( ) / 2, False );

     // Prepare data structure that will allow duplicate identification.

     vec< quad<int,int,int,int64_t> > X( bases.size( ) ); // (e,offset,head2,id1)
     cout << Date( ) << ": making X" << endl;
     #pragma omp parallel for
     for ( int64_t id1 = 0; id1 < (int64_t) bases.size( ); id1++ )
     {    int64_t id2 = ( id1 % 2 == 0 ? id1 + 1 : id1 - 1 );
          const ReadPath& p = paths[id1];
          if ( p.size( ) == 0 ) X[id1] = make_quad( -1, -1, -1, -1 );
          else
          {    int n = 0;
               for ( int j = 0; j < BHEAD; j++ )
               {    n *= 4;
                    n += bases[id2][j];    }
               X[id1] = make_quad( p[0], p.getOffset( ), n, id1 );    }    }

     // Sort, then mark dups.

     cout << Date( ) << ": sorting" << endl;
//     ParallelSort(X);
     sortInPlaceParallel( X.begin(), X.end() );
     cout << Date( ) << ": marking dups" << endl;
//     vec<Bool> dup1( X.size( ), False );
     std::set< int64_t > dup1;
     int64_t ndups = 0, interdups = 0;
     for ( int64_t j = 0; j < X.jsize( ); j++ )
     {    if ( X[j].first < 0 ) continue;
          int64_t k;
          for ( k = j + 1; k < X.jsize( ); k++ )
          {    if ( X[k].first != X[j].first || X[k].second != X[j].second ) break;
               if ( X[k].third != X[j].third ) break;    }
          if ( k - j > 1 )
          {    for ( int64_t l = j; l < k; l++ ) {
                    dup1.insert( X[l].fourth );
              }
               ndups += k - j - 1;
               Bool inter = False;
               int32_t b = bc[ X[j].fourth ];
               for ( int64_t l = j + 1; l < k; l++ )
               {    if ( b == 0 ) b = bc[ X[l].fourth ];
                    else if ( bc[ X[l].fourth ] != b ) inter = True;    }
               if (inter) interdups += k - j - 1;
               if (verbose)
               {    Bool print = False;
                    const int print_freq = 10000;
                    for ( int64_t l = j; l < k; l++ )
                         if ( randomx( ) % print_freq == 0 ) print = True;
                    if (print)
                    {    cout << "\nduplicate group:\n";
                         {    for ( int64_t l = j; l < k; l++ )
                              {    cout << "[" << l-j+1 << "] " << X[l].fourth
                                        << endl;    }    }    }    }    }
          j = k - 1;    }
     if ( ndups > 0 )
          interdup_rate = double(interdups) / double(ndups);
     else
          interdup_rate = 0.0;

     PRINT(dup1.size());

     // Compute qsums.

     cout << Date( ) << ": computing qsums" << endl;
     vec<int> qsum( X.size( ), 0 );
     const int64_t batch = 100000;
     auto quals_clone = quals;
     auto comp = []( quad<int,int,int,int64_t>& t1, quad<int,int,int,int64_t>& t2) {
           int result = compare(t1.fourth,t2.fourth);
           if ( !result ) result = compare(t1.third,t2.third);
           if ( !result ) result = compare(t1.second,t2.second);
           if ( !result ) result = compare(t1.first,t2.first);
           return result;
     };
     cout << Date() << ": ...sorting" << endl;
     sortInPlaceParallel( X.begin(), X.end(), comp);
     cout << Date() << ": ...summing" << endl;
     #pragma omp parallel for schedule(dynamic,1) firstprivate(quals_clone)
     for ( int64_t bi = 0; bi < X.jsize( ); bi += batch )
     {    qualvector q;
          for ( int64_t j = bi; j < Min( bi + batch, X.jsize( ) ); j++ )
          {
               int64_t id1 = X[j].fourth;
               if ( dup1.count( id1 ) == 0 ) continue;
               int64_t id2 = ( id1 % 2 == 0 ? id1 + 1 : id1 - 1 );
               quals_clone[id1].unpack(&q);
               for ( int l = 0; l < (int) q.size( ); l++ )
                    qsum[id1] += q[l];
               quals_clone[id2].unpack(&q);
               for ( int l = 0; l < (int) q.size( ); l++ )
                    qsum[id1] += q[l];    }    }


     // Complete duplicate identification.
     sortInPlaceParallel( X.begin(), X.end() );

     cout << Date( ) << ": finalizing dups" << endl;
     vec< triple<basevector,qualvector,int64_t> > qb;
     vec<Bool> art( bases.size( )/2, False );
     for ( int64_t j = 0; j < X.jsize( ); j++ )
     {    int64_t k;
          for ( k = j + 1; k < X.jsize( ); k++ )
          {    if ( X[k].first != X[j].first || X[k].second != X[j].second ) break;
               if ( X[k].third != X[j].third ) break;    }

          for ( int64_t l = j; l < k; l++ )
          {    if ( l % 100000000 == 0 )
               {    cout << Date( ) << ": start batch " << l/100000000 + 1 << " of "
                         << X.size( )/100000000 << endl;    }    }

          if ( X[j].first < 0 )
          {    j = k - 1;
               continue;    }

          int64_t best = j;
          int q = qsum[X[j].fourth];

          // Check for tie.  Note that we only test this for the winners.

          Bool tie = False;
          for ( int64_t l = j + 1; l < k; l++ )
          {    if ( qsum[X[l].fourth] == q )
               {    tie = True;
                    if ( X[l].fourth < X[best].fourth ) best = l;    }
               else if ( qsum[X[l].fourth] > q )
               {    q = qsum[X[l].fourth];
                    best = l;    }    }
          if (tie)
          {    qb.resize( k - j );
               for ( int64_t l = j; l < k; l++ )
               {    qb[l-j].first = bases[ X[l].fourth ];
                    quals[ X[l].fourth ].unpack( &qb[l-j].second );
                    qb[l-j].third = X[l].fourth/2;    }
               Sort(qb);
               for ( int m = 0; m < (int) qb.isize( ); m++ )
               {    int n;
                    for ( n = m + 1; n < (int) qb.isize( ); n++ )
                    {    if ( qb[n].first != qb[m].first ) break;
                         if ( qb[n].second != qb[m].second ) break;    }
                    for ( int x = m + 1; x < n; x++ )
                         art[ qb[x].third ] = True;
                    m = n - 1;    }    }
          for ( int64_t l = j; l < k; l++ )
               if ( l != best ) dup[ X[l].fourth/2 ] = True;
          j = k - 1;    }
     
     double dup_perc = 100.0 * Sum(dup) / ( bases.size( )/2 );
     cout << setiosflags(ios::fixed) << setprecision(2)
          << dup_perc  << resetiosflags(ios::fixed)
          << "% of pairs appear to be duplicates" << endl;
     StatLogger::log( "dup_perc", dup_perc, "Pct duplicate", true );
     
     cout << setiosflags(ios::fixed) << setprecision(2)
          << ( 100.0 * interdups ) / ndups << resetiosflags(ios::fixed)
          << "% of duplicates involve multiple barcodes" << endl;
     StatLogger::log( "interdup_perc", interdup_rate*100.0, "Pct multi-bc duplicates", false );

     double art_dup_perc = 100.0 * Sum(art) / ( bases.size( )/2 );
     cout << setiosflags(ios::fixed) << setprecision(2)
          << art_dup_perc << resetiosflags(ios::fixed)
          << "% of pairs appear to be artifactual duplicates" << endl;
     StatLogger::log( "art_dup_perc", art_dup_perc, "Pct artifactual dups", false );
     
     if ( double(Sum(art)) / ( bases.size( ) / 2 ) > 0.01 )
     {    cout << "\nWARNING: Too many of your reads appear to be artifactual "
               << "duplicates.\n";
          cout << "It may be that your input files are damaged.\n\n";    }
     cout << Date( ) << ": done marking dups, time used = " << TimeSince(clock)
          << endl;

}

template< class VB, class VQ>
void MarkDups( VB& bases, VQ& quals, ReadPathVecX& paths, const HyperBasevectorX& hb, const vec<int32_t>& bc,
     vec<Bool>& dup, double& interdup_rate, const Bool verbose)
{
     cout << Date( ) << ": MarkDups, version 2" << endl;

     // Heuristics.

     const int BHEAD = 5;

     // Start.

     double clock = WallClockTime( );
     dup.resize_and_set( bases.size( ) / 2, False );

     // Prepare data structure that will allow duplicate identification.

     vec< quad<int,int,int,int64_t> > X( bases.size( ) ); // (e,offset,head2,id1)
     cout << Date( ) << ": making X" << endl;
     #pragma omp parallel for
     for ( int64_t id1 = 0; id1 < (int64_t) bases.size( ); id1++ )
     {    int64_t id2 = ( id1 % 2 == 0 ? id1 + 1 : id1 - 1 );
          ReadPath p; paths.unzip(p,hb,id1); // RPVX
          if ( p.size( ) == 0 ) X[id1] = make_quad( -1, -1, -1, -1 );
          else
          {    int n = 0;
               for ( int j = 0; j < BHEAD; j++ )
               {    n *= 4;
                    n += bases[id2][j];    }
               X[id1] = make_quad( p[0], p.getOffset( ), n, id1 );    }    }

     // Scoping.

     int64_t ndups = 0, interdups = 0;
     vec<int> qsum;
     {
          // Sort, then mark dups.

          cout << Date( ) << ": sorting, mem = " << MemUsageGBString( )
               << ", peak = " << PeakMemUsageGBString( ) << endl;
          sortInPlaceParallel( X.begin(), X.end() );
          cout << Date( ) << ": marking dups, mem = " << MemUsageGBString( )
               << ", peak = " << PeakMemUsageGBString( ) << endl;
          std::set< int64_t > dup1;
          for ( int64_t j = 0; j < X.jsize( ); j++ )
          {    if ( X[j].first < 0 ) continue;
               int64_t k;
               for ( k = j + 1; k < X.jsize( ); k++ )
               {    if ( X[k].first != X[j].first || X[k].second != X[j].second ) 
                         break;
                    if ( X[k].third != X[j].third ) break;    }
               if ( k - j > 1 )
               {    for ( int64_t l = j; l < k; l++ ) {
                         dup1.insert( X[l].fourth );
                   }
                    ndups += k - j - 1;
                    Bool inter = False;
                    int32_t b = bc[ X[j].fourth ];
                    for ( int64_t l = j + 1; l < k; l++ )
                    {    if ( b == 0 ) b = bc[ X[l].fourth ];
                         else if ( bc[ X[l].fourth ] != b ) inter = True;    }
                    if (inter) interdups += k - j - 1;
                    if (verbose)
                    {    Bool print = False;
                         const int print_freq = 10000;
                         for ( int64_t l = j; l < k; l++ )
                              if ( randomx( ) % print_freq == 0 ) print = True;
                         if (print)
                         {    cout << "\nduplicate group:\n";
                              {    for ( int64_t l = j; l < k; l++ )
                                   {    cout << "[" << l-j+1 << "] " << X[l].fourth
                                             << endl;    }    }    }    }    }
               j = k - 1;    }
          interdup_rate = double(interdups) / double(ndups);
          PRINT(dup1.size());
     
          // Compute qsums.

          cout << Date( ) << ": computing qsums, mem = " << MemUsageGBString( )
               << ", peak = " << PeakMemUsageGBString( ) << endl;
          qsum.resize( X.size( ), 0 );
          const int64_t batch = 100000;
          auto quals_clone = quals;
          auto comp = []( 
               quad<int,int,int,int64_t>& t1, quad<int,int,int,int64_t>& t2) {
               int result = compare(t1.fourth,t2.fourth);
               if ( !result ) result = compare(t1.third,t2.third);
               if ( !result ) result = compare(t1.second,t2.second);
               if ( !result ) result = compare(t1.first,t2.first);
               return result;
          };
          cout << Date() << ": ...sorting" << endl;
          sortInPlaceParallel( X.begin(), X.end(), comp);
          cout << Date() << ": ...summing" << endl;
          #pragma omp parallel for schedule(dynamic,1) firstprivate(quals_clone)
          for ( int64_t bi = 0; bi < X.jsize( ); bi += batch )
          {    qualvector q;
               for ( int64_t j = bi; j < Min( bi + batch, X.jsize( ) ); j++ )
               {
                    int64_t id1 = X[j].fourth;
                    if ( dup1.count( id1 ) == 0 ) continue;
                    int64_t id2 = ( id1 % 2 == 0 ? id1 + 1 : id1 - 1 );
                    quals_clone[id1].unpack(&q);
                    for ( int l = 0; l < (int) q.size( ); l++ ) qsum[id1] += q[l];
                    quals_clone[id2].unpack(&q);
                    for ( int l = 0; l < (int) q.size( ); l++ )
                         qsum[id1] += q[l];    }    }    }

     // Complete duplicate identification.

     sortInPlaceParallel( X.begin(), X.end() );
     cout << Date( ) << ": finalizing dups, mem = " << MemUsageGBString( )
          << ", peak = " << PeakMemUsageGBString( ) << endl;
     vec< triple<basevector,qualvector,int64_t> > qb;
     vec<Bool> art( bases.size( )/2, False );
     for ( int64_t j = 0; j < X.jsize( ); j++ )
     {    int64_t k;
          for ( k = j + 1; k < X.jsize( ); k++ )
          {    if ( X[k].first != X[j].first || X[k].second != X[j].second ) break;
               if ( X[k].third != X[j].third ) break;    }

          for ( int64_t l = j; l < k; l++ )
          {    if ( l % 100000000 == 0 )
               {    cout << Date( ) << ": start batch " << l/100000000 + 1 << " of "
                         << X.size( )/100000000 + 1 << endl;    }    }

          if ( X[j].first < 0 )
          {    j = k - 1;
               continue;    }

          int64_t best = j;
          int q = qsum[X[j].fourth];

          // Check for tie.  Note that we only test this for the winners.

          Bool tie = False;
          for ( int64_t l = j + 1; l < k; l++ )
          {    if ( qsum[X[l].fourth] == q )
               {    tie = True;
                    if ( X[l].fourth < X[best].fourth ) best = l;    }
               else if ( qsum[X[l].fourth] > q )
               {    q = qsum[X[l].fourth];
                    best = l;    }    }
          if (tie)
          {    qb.resize( k - j );
               for ( int64_t l = j; l < k; l++ )
               {    qb[l-j].first = bases[ X[l].fourth ];
                    quals[ X[l].fourth ].unpack( &qb[l-j].second );
                    qb[l-j].third = X[l].fourth/2;    }
               Sort(qb);
               for ( int m = 0; m < (int) qb.isize( ); m++ )
               {    int n;
                    for ( n = m + 1; n < (int) qb.isize( ); n++ )
                    {    if ( qb[n].first != qb[m].first ) break;
                         if ( qb[n].second != qb[m].second ) break;    }
                    for ( int x = m + 1; x < n; x++ )
                         art[ qb[x].third ] = True;
                    m = n - 1;    }    }
          for ( int64_t l = j; l < k; l++ )
               if ( l != best ) dup[ X[l].fourth/2 ] = True;
          j = k - 1;    }

     // Compute stats.

     double dup_perc = 100.0 * Sum(dup) / ( bases.size( )/2 );
     cout << setiosflags(ios::fixed) << setprecision(2)
          << dup_perc  << resetiosflags(ios::fixed)
          << "% of pairs appear to be duplicates" << endl;
     StatLogger::log( "dup_perc", dup_perc, "Pct duplicate", true );
     cout << setiosflags(ios::fixed) << setprecision(2)
          << ( 100.0 * interdups ) / ndups << resetiosflags(ios::fixed)
          << "% of duplicates involve multiple barcodes" << endl;
     StatLogger::log( 
          "interdup_perc", interdup_rate*100.0, "Pct multi-bc duplicates", false );
     double art_dup_perc = 100.0 * Sum(art) / ( bases.size( )/2 );
     cout << setiosflags(ios::fixed) << setprecision(2)
          << art_dup_perc << resetiosflags(ios::fixed)
          << "% of pairs appear to be artifactual duplicates" << endl;
     StatLogger::log( "art_dup_perc", art_dup_perc, "Pct artifactual dups", false );
     if ( double(Sum(art)) / ( bases.size( ) / 2 ) > 0.01 )
     {    cout << "\nWARNING: Too many of your reads appear to be artifactual "
               << "duplicates.\n";
          cout << "It may be that your input files are damaged.\n\n";    }
     cout << Date( ) << ": done marking dups, time used = " << TimeSince(clock)
          << endl;    }

template void MarkDups( VirtualMasterVec<basevector>& bases,
     VirtualMasterVec<PQVec>&, ReadPathVecX&, const HyperBasevectorX&, const vec<int32_t>&,
     vec<Bool>&, double&, const Bool );

template void MarkDups( MasterVec<basevector>&, MasterVec<PQVec>&,
     ReadPathVecX&,const HyperBasevectorX&, const vec<int32_t>&, vec<Bool>&, double&, const Bool );

template void MarkDups( MasterVec<basevector>& bases,
     VirtualMasterVec<PQVec>&, ReadPathVecX&,const HyperBasevectorX&, const vec<int32_t>&,
     vec<Bool>&, double&, const Bool);


template void MarkDups( VirtualMasterVec<basevector>& bases,
     VirtualMasterVec<PQVec>&, VirtualMasterVec<ReadPath>&, const vec<int32_t>&,
     vec<Bool>&, double&, const Bool );

template void MarkDups( MasterVec<basevector>&, MasterVec<PQVec>&,
     MasterVec<ReadPath>&, const vec<int32_t>&, vec<Bool>&, double&, const Bool );

template void MarkDups( MasterVec<basevector>& bases,
     VirtualMasterVec<PQVec>&, MasterVec<ReadPath>&, const vec<int32_t>&,
     vec<Bool>&, double&, const Bool );

void AllTinksCore( const HyperBasevectorX& hb, const vec<int>& inv,
     const ReadPathVecX& paths, const vec<int64_t>& bci, const VecIntVec & ebcx,
     vec< triple<int,int,int> >& qept, const int verbosity )
{    if ( verbosity >= 1 ) cout << Date( ) << ": start AllTinks" << endl;

     // Algorithmic heuristics.

     const int min_bc = 4;
     const int min_reads = 100;
     const int max_reads = 10000;

     // Computational performance heuristics.

     const int npasses = 20;
     const int64_t batches = 20;

     // Identify line-like edge groups.  These are simple lines with binary
     // bubbles.  Map all edges in the group to the first edge.

     if ( verbosity >= 1 ) cout << Date( ) << ": defining groups" << endl;
     vec<int> to_first( hb.E( ), vec<int>::IDENTITY );
     vec<vec<int>> groups( hb.E( ) );
     #pragma omp parallel for
     for ( int e = 0; e < hb.E( ); e++ )
     {    Bool starter = True;
          int w = hb.ToLeft(e);
          if ( hb.To(w).size( ) == 2 && hb.From(w).size( ) == 1 )
          {    if ( hb.To(w)[0] == hb.To(w)[1] )
               {    int v = hb.To(w)[0];
                    if ( v != w && hb.To(v).size( ) == 1 && hb.From(v).size( ) == 2 )
                    {    starter = False;    }    }    }
          if (starter)
          {    groups[e] = {e};
               int f = e;
               while(1)
               {    int v = hb.ToRight(f);
                    if ( hb.From(v).size( ) != 2 || hb.To(v).size( ) != 1 ) break;
                    if ( hb.From(v)[0] != hb.From(v)[1] ) break;
                    int w = hb.From(v)[0];
                    if ( hb.To(w).size( ) != 2 || hb.From(w).size( ) != 1 ) break;
                    if ( w == v ) break;
                    int g1 = hb.IFrom(v,0), g2 = hb.IFrom(v,1);
                    f = hb.IFrom(w,0);
                    to_first[g1] = to_first[g2] = to_first[f] = e;
                    groups[e].push_back( g1, g2, f );    }    }    }

     // Define a few constants.

     if ( verbosity >= 1 ) cout << Date( ) << ": defining constants" << endl;
     int NBC = bci.size( ) - 1;
     if ( verbosity >= 1 ) cout << Date( ) << ": have " << ToStringAddCommas((int64_t) paths.size( )) << " reads" << endl;
     int64_t unbar = bci[1];

     // Find good barcodes.

     double clock = WallClockTime( );
     vec<int> gbc;
     for ( int b = 1; b < NBC; b++ )
     {    int64_t nbc = bci[b+1] - bci[b];
          if ( nbc < min_reads || nbc > max_reads ) continue;
          gbc.push_back(b);    }
     if ( verbosity >= 1 ) cout << Date( ) << ": have " << ToStringAddCommas( gbc.size( ) )
          << " qualified barcodes" << endl;

     // Go through ten passes.

     if ( verbosity >= 1 ) cout << Date( ) << ": start main loop" << endl;
     double mclock = WallClockTime( );
     vec<pair<int,int>> ep, qep;
     for ( int d = 0; d < npasses; d++ )
     {
          // Find edge pairs.

          if ( verbosity >= 2 ) cout << Date( ) << ": finding edge pairs" << endl;
          ep.clear( );
          #pragma omp parallel for
          for ( int ig = 0; ig < gbc.isize( ); ig++ )
          {    int b = gbc[ig];
               vec<pair<int,int>> epb;

               // Find all the edges on the barcode.

               vec<int> es;
               for ( int64_t id = bci[b]; id < bci[b+1]; id++ )
               {    
                    ReadPath p; paths.unzip(p,hb,id);
                    for ( int j = 0; j < (int) p.size( ); j++ )
                    {    int e = to_first[ p[j] ], n = 0;
                         for ( int l = 0; l < groups[e].isize( ); l++ )
                              n = Max( n, hb.Kmers( groups[e][l] ) );
                         if ( n < MIN_KMERS2 ) continue;
                         if ( to_first[ inv[e] ] < e ) e = to_first[ inv[e] ];
                         es.push_back(e);    }    }
               UniqueSort(es);

               // Find edge pairs.

               for ( int j1 = 0; j1 < es.isize( ); j1++ )
               {    int e1 = es[j1], n1 = 0;
                    for ( int l = 0; l < groups[e1].isize( ); l++ )
                         n1 = Max( n1, hb.Kmers( groups[e1][l] ) );
                    if ( n1 < MIN_KMERS1 ) continue;
                    if ( n1 % npasses != d ) continue;
                    for ( int j2 = 0; j2 < es.isize( ); j2++ )
                    {    if ( j2 == j1 ) continue;
                         int e2 = es[j2];
                         epb.push( e1, e2 );    }    }
               #pragma omp critical
               {    ep.append(epb);    }    }

          // Sort edge pairs.

          if ( verbosity >= 2 ) cout << Date( ) << ": sorting edge pairs" << endl;
          ParallelSort(ep);

          // Find qualified edge pairs.

          if ( verbosity >= 2 ) cout << Date( ) << ": finding qualified edge pairs" << endl;
          vec<int64_t> starts( batches+1, 0 );
          for ( int i = 0; i <= batches; i++ )
               starts[i] = ( i * ep.jsize( ) ) / batches;
          for ( int i = 0; i < batches; i++ )
          {    if ( starts[i] > 0 && ep[ starts[i] ] == ep[ starts[i] - 1 ] )
                    starts[i]--;    }
          vec< vec< pair<int,int> > > qepb(batches);
          #pragma omp parallel for
          for ( int b = 0; b < batches; b++ )
          {    for ( int64_t i = starts[b]; i < starts[b+1]; i++ )
               {    int64_t j;
                    for ( j = i + 1; j < starts[b+1]; j++ )
                         if ( ep[j] != ep[i] ) break;
                    if ( j - i >= min_bc ) qepb[b].push_back( ep[i] );
                    i = j - 1;    }    }
          for ( int b = 0; b < batches; b++ )
               qep.append( qepb[b] );    }
     if ( verbosity >= 1 ) cout << Date( ) << ": time used in main loop = " << TimeSince(mclock) << endl;

     // Expand edge pairs by groups.  Bad in multiple ways, may need wholistic
     // solution.

     if ( verbosity >= 1 ) cout << Date( ) << ": expanding " << ToStringAddCommas( qep.size( ) )
          << " edge pairs by groups" << endl;
     Destroy(ep);
     {    int64_t N = qep.size( );
          const int batch = 1000000;
          vec<vec<pair<int,int>>> adds( N/batch + 1 );
          #pragma omp parallel for
          for ( int64_t bi = 0; bi < N; bi += batch )
          {    for ( int64_t i = bi; i < Min( bi + batch, N ); i++ )
               {    int e1 = qep[i].first, e2 = qep[i].second;
                    for ( int f1 : groups[e1] )
                    {    if ( hb.Kmers(f1) < MIN_KMERS1 ) continue;
                         int xf1 = Min( f1, inv[f1] );
                         for ( int f2 : groups[e2] )
                         {    if ( hb.Kmers(f2) < MIN_KMERS2 ) continue;
                              adds[bi/batch].push_back(
                                   make_pair( xf1, f2 ) );    }    }    }    }
          qep.clear( );
          if ( verbosity >= 1 ) cout << Date( ) << ": appending" << endl;
          for ( int i = 0; i < adds.isize( ); i++ )
               qep.append( adds[i] );    }
     if ( verbosity >= 1 ) cout << Date( ) << ": now have " << ToStringAddCommas( qep.size( ) )
          << " edge pairs" << endl;

     // Append counts.

     if ( verbosity >= 1 ) cout << Date( ) << ": appending counts" << endl;

     qept.resize( qep.size( ) );
     #pragma omp parallel for
     for ( int64_t i = 0; i < qep.jsize( ); i++ )
     {    int e1 = qep[i].first ,e2 = qep[i].second;
          qept[i].first = e1, qept[i].second = e2;
          // compute intersection of barcode vectors
          // this works because BOTH vectors are sorted
          // and all elements are unique
          auto & b1 = ebcx[e1];
          auto & b2 = ebcx[e2];
          auto p1 = b1.begin();
          auto p2 = b2.begin();
          int common = 0;
          while (p1 != b1.end() && p2 != b2.end() ) {
               if ( *p1 == *p2 ) {
                    common++;
                    p1++;
                    p2++;
               }
               else if ( *p1 < *p2 )
                    p1++;
               else
                    p2++;
          }
          qept[i].third = common;
     }

     // Sort and write.

     if ( verbosity >= 1 ) cout << Date( ) << ": sorting pairs" << endl;
     ParallelSort(qept);
     if ( verbosity >= 1 ) cout << Date( ) << ": " << ToStringAddCommas( qept.size( ) )
          << " pairs of edges connected by several barcodes" << endl;
     if ( verbosity >= 1 ) cout << Date( ) << ": done, time used = " << TimeSince(clock)
          << ", peak mem = " << PeakMemUsageGBString( ) << endl;    }

double SinglesRate( const vec<int>& kmers, const vec<int>& inv,
     const vec<Bool>& dup, const VecULongVec& paths_index, const vec<int64_t>& bci )
{
     // Expand barcode index.

     vec<int32_t> bc( bci.back( ), -1 );
     #pragma omp parallel for
     for ( int b = 0; b < bci.isize( ) - 1; b++ )
     {    int64_t start = bci[b], stop = bci[b+1];
          for ( int64_t j = start; j < stop; j++ )
               bc[j] = b;    }

     // Compute singles rate.

     int64_t n10 = 0, d10 = 0;
     for ( int e = 0; e < kmers.isize( ); e++ )
     {    if ( kmers[e] >= 10000 && kmers[e] <= 12000 )
          {    vec< pair<int,int64_t> > x;
               int re = inv[e];
               for ( int j = 0; j < (int) paths_index[e].size( ); j++ )
               {    int64_t id = paths_index[e][j];
                    if ( !dup[id/2] ) x.push( bc[id], id );    }
               for ( int j = 0; j < (int) paths_index[re].size( ); j++ )
               {    int64_t id = paths_index[re][j];
                    if ( !dup[id/2] ) x.push( bc[id], id );    }
               UniqueSort(x);
               d10 += x.isize( );
               for ( int j = 0; j < x.isize( ); j++ )
               {    if ( x[j].first == 0 ) continue;
                    int k;
                    for ( k = j + 1; k < x.isize( ); k++ )
                         if ( x[k].first != x[j].first ) break;
                    if ( k - j >= 3
                         || ( k - j == 2 && x[j+1].second/2 != x[j].second/2 ) )
                    {    n10 += k - j;    }
                    j = k - 1;    }    }    }
     return 1.0 - ( double(n10) / double(d10) );    }

void MakeClosures( const HyperBasevectorX& hb, const vec<int>& inv,
     ReadPathVecX& paths, VecULongVec& paths_index,
     const vec<Bool>& dup, const vec<Bool>& bad, vec<vec<int>>& all_closures,
     const Bool verbose, String paths_index_file )
{
     // Exactly one of paths_index or paths_index_file MUST be non-zero
     ForceAssert( (paths_index.size() == 0) ^ (paths_index_file.size() == 0) );

     // Make closures.
     
     if (verbose)
     {    cout << Date( ) << ": making closures, mem = " << MemUsageGBString( )
               << ", peak = " << PeakMemUsageGBString( ) << endl;    }
     int64_t opcount = 0;
     if ( paths_index_file.size() > 0 )
     {    cout << Date ( ) << ": loading paths index" << endl;
          paths_index.clear( );
          paths_index.ReadAll ( paths_index_file );

          if (verbose)
          {    cout << Date( ) << ": load complete, mem = " << MemUsageGBString( )
                    << ", peak = " << PeakMemUsageGBString( ) << endl;    }    }

     Closer( hb, inv, paths, paths_index, dup, bad, all_closures, opcount,
          True, ( verbose ? 1 : 0 ) );
     
     if ( paths_index_file.size() > 0 ) {
          MEM(before_destroy_paths_index);
          Destroy( paths_index );
          MEM(after_destroy_paths_index);
     }
     // Double.

     if (verbose)
     {    cout << Date( ) << ": doubling closures, mem = " << MemUsageGBString( )
               << ", peak = " << PeakMemUsageGBString( ) << endl;    }
     int64_t nac = all_closures.size( );
     all_closures.resize( 2*nac );
     #pragma omp parallel for schedule( dynamic, 10000 )
     for ( int64_t i = 0; i < nac; i++ )
     {    vec<int>& x = all_closures[nac+i];
          x = all_closures[i];
          x.ReverseMe( );
          for ( int j = 0; j < x.isize( ); j++ )
               x[j] = inv[ x[j] ];    }
     if (verbose)
     {    cout << Date( ) << ": sorting closures, mem = " << MemUsageGBString( ) 
               << ", peak = " << PeakMemUsageGBString( ) << endl;    }
     UniqueSort(all_closures);

     // Add long edges that were somehow missed.

     if (verbose)
     {    cout << Date( ) << ": adding back long edges" << endl;    }
     int64_t N = all_closures.size( );
     {    vec<Bool> present( hb.E( ), False );
          for ( int64_t i = 0; i < N; i++ )
          for ( int j = 0; j < all_closures[i].isize( ); j++ )
               present[ all_closures[i][j] ] = True;
          for ( int e = 0; e < hb.E( ); e++ )
          {    if ( !present[e] && hb.Kmers(e) >= 200 ) 
                    all_closures.push_back( {e} );     }    }
     if (verbose) cout << Date( ) << ": sorting again" << endl;
     UniqueSort(all_closures);

     // Index closures.

     N = all_closures.size( );
     if (verbose)
     {    int64_t total = 0;
          for ( int64_t i = 0; i < N; i++ )
               total += all_closures[i].size( );
          cout << Date( ) << ": " << ToStringAddCommas(N) << " closures "
               << "having in total " << ToStringAddCommas(total) << " edges" << endl;
          cout << Date( ) << ": indexing closures, peak mem = "
               << PeakMemUsageGBString( ) << endl;    }
     vec<vec<int>> ci( hb.E( ) );
     for ( int64_t i = 0; i < N; i++ )
     for ( int j = 0; j < all_closures[i].isize( ); j++ )
          ci[ all_closures[i][j] ].push_back(i);
     if (verbose) cout << Date( ) << ": uniquesorting" << endl;
     #pragma omp parallel for schedule( dynamic, 1000 )
     for ( int e = 0; e < hb.E( ); e++ )
          UniqueSort( ci[e] );

     // Kill short subset closures.

     if (verbose) cout << Date( ) << ": killing short subset closures" << endl;
     const int MIN_OVER = 200 - (hb.K()-1);
     vec<Bool> sdel( all_closures.size( ), False );
     #pragma omp parallel for schedule( dynamic, 1000 )
     for ( int64_t i = 0; i < all_closures.jsize( ); i++ )
     {    vec<int>& x = all_closures[i];
          int len = 0;
          for ( int j = 0; j < x.isize( ); j++ )
               len += hb.Kmers( x[j] );
          if ( len >= MIN_OVER ) continue;
          for ( int j = 0; j < ci[ x[0] ].isize( ); j++ )
          {    int64_t ip = ci[ x[0] ][j];
               if ( ip == i ) continue;
               if ( all_closures[ip].Contains(x) )
               {    sdel[i] = True;
                    break;    }    }    }
     EraseIf( all_closures, sdel );    }
