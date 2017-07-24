// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>

#include "CoreTools.h"
#include "ParallelVecUtilities.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "10X/Closer.h"
#include "10X/paths/ReadPathVecX.h"

void DefinePairSet( const HyperBasevectorX& hb, const vec<int>& inv, 
     ReadPathVecX& paths, const vec<Bool>& dup, const vec<Bool>& bad,
     vec<int64_t>& ppids, const Bool verbose )
{
     const int64_t batch = 10000;
     const int64_t num_pairs = paths.size()/2;
     vec< pair <vec<int>, vec<int> > > proc;

     // Reserve space.

     if (verbose) cout << Date( ) << ": reserving space" << endl;
     int64_t count = 0;
     #pragma omp parallel for schedule( dynamic, 1 )
     for ( int64_t bi = 0; bi < num_pairs; bi += batch )
     {    int64_t counti = 0;
          const int64_t pid_max = Min( bi + batch, num_pairs );
          for ( int64_t pid = bi; pid < pid_max; pid++ )
          {    if ( dup[pid] || bad[pid] ) continue;
               int64_t id1 = 2*pid, id2 = 2*pid+1;
               ReadPath p1, p2; paths.unzip(p1,hb,id1); paths.unzip(p2,hb,id2);
               if ( p1.size( ) == 0 || p2.size( ) == 0 ) continue;
               int v = hb.ToRight( p1.back( ) ), w = hb.ToLeft( inv[ p2.back( ) ] );
               if ( hb.From(v).size( ) == 0 || hb.To(w).size( ) == 0 ) continue;
               counti++;    }
          #pragma omp critical
          {    count += counti;    }    }
     proc.reserve(2*count), ppids.reserve(count);

     // Now do the computation.

     if (verbose) cout << Date( ) << ": finding pairs" << endl;
     #pragma omp parallel for schedule( dynamic, 1 )
     for ( int64_t bi = 0; bi < num_pairs; bi += batch )
     {    vec< pair <vec<int>, vec<int> > > proci;
          vec<int64_t> ppidsi;
          vec<int> x1, x2;
          const int64_t pid_max = Min( bi + batch, num_pairs );
          for ( int64_t pid = bi; pid < pid_max; pid++ )
          {    if ( dup[pid] || bad[pid] ) continue;
               int64_t id1 = 2*pid, id2 = 2*pid+1;
               ReadPath p1, p2; paths.unzip(p1,hb,id1); paths.unzip(p2,hb,id2);
               if ( p1.size( ) == 0 || p2.size( ) == 0 ) continue;
               x1.clear( ), x2.clear( );
               for ( int i = 0; i < (int) p1.size( ); i++ )
                    x1.push_back( p1[i] );
               for ( int i = (int) p2.size( ) - 1; i >= 0; i-- )
                    x2.push_back( inv[ p2[i] ] );
               int v = hb.ToRight( x1.back( ) ), w = hb.ToLeft( x2.front( ) );
               if ( hb.From(v).size( ) == 0 || hb.To(w).size( ) == 0 ) continue;
               proci.push( x1, x2 );
               ppidsi.push_back(pid);    }
          #pragma omp critical
          {    proc.append(proci); ppids.append(ppidsi);    }    }
     if (verbose) 
     {    cout << Date( ) << ": sorting proc, mem = " << MemUsageGBString( )
               << ", peak = " << PeakMemUsageGBString( ) << endl;    }
     ParallelUniqueSortSync( proc, ppids );
     if (verbose) 
     {    cout << Date( ) << ": sort complete, mem = " << MemUsageGBString( )
               << ", peak = " << PeakMemUsageGBString( ) << endl;    }    }

template<class VPI > void ClosePair( const vec<int>& x1, 
     const vec<int>& x2, const int64_t pid, const HyperBasevectorX& hb, 
     const vec<int>& inv, ReadPathVecX& paths, VPI& paths_index, int64_t& opcounti, 
     vec<vec<int>>& closures, ostringstream& out, const Bool VERBOSE )
{
     // out << "x: " << printSeq(x1) << " .. " << printSeq(x2) << endl;

     // Heuristics.  
     // Note that the MAX_FREQ requirement is problematic, for at least two reasons:
     // 1. A read may not have a kmer of low enough frequency.
     // 2. The behavior depends on coverage.

     const int MAX_LEVEL = 5;
     const int MAX_EXTS = 50;
     const int MAX_MIDDLE = 300;
     const int MIN_OVER = 20;
     const int MAX_FREQ = 1000;
     const int MAX_OPS = 20000;

     // Try for easy closure.

     if ( hb.ToRight( x1.back( ) ) == hb.ToLeft( x2.front( ) ) )
     {    vec<int> x = x1;
          x.append(x2);
          closures.push_back(x);
          // out << "easy closure = " << printSeq(x) << endl;
          // cout << out.str( );
          return;    }
     if ( x1.back( ) == x2.front( ) )
     {    vec<int> x = x1;
          for ( int j = 1; j < x2.isize( ); j++ )
               x.push_back( x2[j] );
          closures.push_back(x);
          // out << "easy closure = " << printSeq(x) << endl;
          // cout << out.str( );
          return;    }
     vec<vec<int>> joins;
     for ( int i1 = 0; i1 < x1.isize( ); i1++ )
     for ( int i2 = 0; i2 < x2.isize( ); i2++ )
     {    if ( x1[i1] == x2[i2] )
          {    Bool mismatch = False;
               for ( int j1 = 0; j1 < x1.isize( ); j1++ )
               {    int j2 = j1 + i2 - i1;
                    if ( j2 < 0 || j2 >= x2.isize( ) ) continue;
                    if ( x1[j1] != x2[j2] ) 
                    {    mismatch = True;
                         break;    }    }
               if (mismatch) continue;
               vec<int> x = x1;
               x.resize( i1 + 1 );
               for ( int j2 = i2 + 1; j2 < x2.isize( ); j2++ )
                    x.push_back( x2[j2] );
               joins.push_back(x);    }    }
     UniqueSort(joins);
     if ( joins.solo( ) )
     {    const vec<int>& x = joins[0];
          closures.push_back(x);
          // out << "easy closure = " << printSeq(x) << endl;
          // cout << out.str( );
          return;    }

     // State what we're doing.

     out << "==================================================================="
          << "=================\n\n";
     out << "EXAMINING PAIR " << ( pid >= 0 ? "+" : "" ) << pid << "\n\n";

     // Set up some handy functions.

     auto wink = [&]( vec<int>& x )
     {    x.ReverseMe( );
          for ( int j = 0; j < x.isize( ); j++ )
               x[j] = inv[ x[j] ];   };
     auto is_match = [&]( const vec<int>& x1, const int p1, const vec<int>& x2,
          const int p2 )
     {    for ( int j1 = 0; j1 < x1.isize( ); j1++ )
          {    int j2 = j1 + p2 - p1;
               if ( j2 < 0 || j2 >= x2.isize( ) ) continue;
               if ( x1[j1] != x2[j2] ) return False;    }
          return True;    };

     // Second approach.  Find the lowest multiplicity edge in the read pair, find
     // all reads that contain it and their partners if pointing in, and see if the
     // pair can be closed using only those pairs.

     /*
     {    vec<int> y1(x1), y2(x2);
          int low1 = 1000000000, low2 = 1000000000;
          int z1 = -1, z2 = -1;
          for ( int j = 0; j < x1.isize( ); j++ )
          {    int n = hb.Kmers( x1[j] );
               if ( n < low1 )
               {    low1 = n;
                    z1 = j;    }    }
          for ( int j = 0; j < x2.isize( ); j++ )
          {    int n = hb.Kmers( x2[j] );
               if ( n < low2 )
               {    low2 = n;
                    z2 = j;    }    }
          if ( low2 < low1 )
          {    swap(low1,low2), swap(z1,z2), swap(y1,y2);
               wink(y1), wink(y2);
               z1 = y1.isize( ) - z1 - 1, z2 = y2.isize( ) - z2 - 1;    }
          int e = y1[z1];
          int re = inv[e];
          vec<vec<int>> xs;
          for ( int j = 0; j < (int) paths_index[e].size( ); j++ )
          {    int64_t id1 = paths_index[e][j];
               int64_t id2 = ( id1 % 2 == 0 ? id1+1 : id1-1 );
               ReadPath p1, p2; paths.unzip(p1,hb,id1); paths.unzip(p2,hb,id2);
               vec<int> x1;
               for ( int l = 0; l < (int) p1.size( ); l++ )
                    x1.push_back( p1[l] );
               Bool match = False;
               for ( int l = 0; l < x1.isize( ); l++ )
               {    if ( x1[l] == e && is_match( y1, z1, x1, l ) )
                    {    match = True;
                         break;    }    }
               if ( !match ) continue;
               xs.push_back(x1);
               if ( p2.size( ) > 0 )
               {    vec<int> x2;
                    for ( int l = 0; l < (int) p2.size( ); l++ )
                         x2.push_back( p2[l] );
                    wink(x2);
                    xs.push_back(x2);    }    }
          for ( int j = 0; j < (int) paths_index[re].size( ); j++ )
          {    int64_t id1 = paths_index[re][j];
               ReadPath p1; paths.unzip(p1,hb,id1);
               vec<int> x1;
               for ( int l = 0; l < (int) p1.size( ); l++ )
                    x1.push_back( p1[l] );
               wink(x1);
               Bool match = False;
               for ( int l = 0; l < x1.isize( ); l++ )
               {    if ( x1[l] == e && is_match( y1, z1, x1, l ) )
                    {    match = True;
                         break;    }    }
               if ( !match ) continue;
               xs.push_back(x1);    }
          UniqueSort(xs);
          out << "would have " << xs.size( ) << " reads in patch set\n" << endl;
          vec< triple<int,int,int> > ps;
          for ( int i = 0; i < xs.isize( ); i++ )
          for ( int j = 0; j < xs[i].isize( ); j++ )
               ps.push( xs[i][j], i, j );
          Sort(ps);
          opcounti += ps.size( );    }
     */

     // Decide if we should flip the pair.  Not necessarily advantageous.

     Bool flipper = False;
     /*
     int nk = 0, n1 = 1000000000, n2 = 1000000000;
     for ( int b = x1.isize( ) - 1; b >= 0; b-- )
     {    int f = x1[b];
          int n = paths_index[f].size( ) + paths_index[ inv[f] ].size( );
          n1 = Min( n, n1 );
          nk += hb.Kmers(f);
          if ( nk >= MIN_OVER ) break;    }
     nk = 0;
     for ( int b = 0; b < x2.isize( ); b++ )
     {    int f = x2[b];
          int n = paths_index[f].size( ) + paths_index[ inv[f] ].size( );
          n2 = Min( n, n2 );
          nk += hb.Kmers(f);
          if ( nk >= MIN_OVER ) break;    }
     if ( n2 < n1 )
     {    flipper = True;
          out << "flipping" << endl;
          swap( x1, x2 );
          x1.ReverseMe( ), x2.ReverseMe( );
          for ( int j = 0; j < x1.isize( ); j++ )
               x1[j] = inv[ x1[j] ];
          for ( int j = 0; j < x2.isize( ); j++ )
               x2[j] = inv[ x2[j] ];    }
     */

     // Now do the general case.

     vec<vec<int>> xs = {x1};
     // out << Date( ) << ": extending" << endl;
     for ( int pass = 1; pass <= MAX_LEVEL; pass++ )
     {    if ( opcounti > MAX_OPS ) break;
          vec<vec<int>> xs2;
          for ( int i = 0; i < xs.isize( ); i++ )
          {    if ( opcounti > MAX_OPS ) break;
               const vec<int>& x = xs[i];

               int bfi;
               int n;
               for ( bfi = x.isize( ) - 1; bfi >= 0; bfi-- )
               {    int f = x[bfi];
                    n = paths_index[f].size( ) + paths_index[ inv[f] ].size( );
                    if ( n < MAX_FREQ ) break;    }
               if ( bfi < 0 ) continue;

               // Implement MIN_OVER requirement.

               int nk = 0, b;
               for ( b = x.isize( ) - 1; b >= 0; b-- )
               {    int f = x[b];
                    int n2 = paths_index[f].size( ) + paths_index[ inv[f] ].size( );
                    if ( n2 < n ) 
                    {    n = n2;
                         bfi = b;    }
                    nk += hb.Kmers(f);
                    if ( nk >= MIN_OVER ) break;    }
               if ( b < 0 ) continue;
     
               int e = x[bfi];
               int re = inv[e];
               for ( int mpass = 1; mpass <= 2; mpass++ )
               {    if ( opcounti > MAX_OPS ) break;
                    int f = ( mpass == 1 ? e : re );
                    int nind = paths_index[f].size( ) // XXXXXXXXXXXXXXXXXXXXXXXXXXX
                         + paths_index[ inv[f] ].size( ); // XXXXXXXXXXXXXXXXXXXXXXX
                    PRINT5_TO( // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                         out, pid, pass, mpass, f, nind ); // XXXXXXXXXXXXXXXXXXXXXX
                    if ( nind > 10000 ) out << "x = " << printSeq(x) << endl; // XXX
                    if ( paths_index[f].size( ) > 100000 ) continue; // ************
                    for ( auto id : paths_index[f] )
                    {    if ( opcounti > MAX_OPS ) break;
                         ReadPath p; paths.unzip(p,hb,id); 
                         if ( hb.SinkEdge( p.back( ) ) ) continue;
                         vec<int> z;
                         if ( mpass == 1 ) for ( auto r : p ) z.push_back(r);
                         else
                         {    for ( int j = (int) p.size( ) - 1; j >= 0; j-- )
                                   z.push_back( inv[ p[j] ] );    }
                         int tail = x.isize( ) - 1 - bfi;
                         for ( int j = 0; j < z.isize( ) - tail; j++ )
                         {    if ( z[j] == e )
                              {    
                                   // bfi  <-->  j
                                   // b    <-->  j + b - bfi

                                   if ( j + b - bfi < 0 ) continue;

                                   Bool mismatch = False;
                                   for ( int m = j + tail; m >= 0; m-- )
                                   {    if ( m == j ) continue;
                                        int xpos = bfi + m - j;
                                        opcounti++;
                                        if ( xpos < 0 ) break;
                                        if ( x[xpos] != z[m] )
                                        {    mismatch = True;
                                             break;    }    }
                                   if (mismatch) continue;
                                   if ( opcounti > MAX_OPS ) break;
                                   vec<int> y = x;
                                   y.resize( bfi + 1 );
                                   for ( int l = j + 1; l < (int) z.size( ); l++ )
                                   {    y.push_back( z[l] );
                                        opcounti++;    }
                                   int nk = 0;
                                   for ( int l = x1.isize( ); 
                                        l < y.isize( ) - 1; l++ )
                                   {    nk += hb.Kmers( y[l] );    }
                                   if ( nk > MAX_MIDDLE ) continue;
                                   xs2.push_back(y);    }    }    }    }    }
          UniqueSort(xs2);
          // if ( P >= 0 ) PRINT2_TO( out, pass, xs2.size( ) );
          if ( xs2.isize( ) > MAX_EXTS ) break;
          xs = xs2;

          // Look for closures.

          for ( int i = 0; i < xs.isize( ); i++ )
          {    const vec<int>& x = xs[i];
               for ( int j = 0; j < x.isize( ); j++ )
               {    if ( x[j] == x1.back( ) )
                    {    Bool mismatch = False;
                         for ( int m = j - 1; m >= 0; m-- )
                         {    int x1pos = x1.isize( ) - 1 - ( j - m );
                              if ( x1pos < 0 ) break;
                              if ( x1[x1pos] != x[m] )
                              {    mismatch = True;
                                   break;    }    }
                         if (mismatch) continue;
                         for ( int k = j + 1; k < x.isize( ); k++ )
                         {    if ( x[k] == x2[0] )
                              {    Bool mismatch = False;
                                   for ( int m = k + 1; m < x.isize( ); m++ )
                                   {    int x2pos = m - k;
                                        if ( x2pos >= x2.isize( ) ) break;
                                        if ( x2[x2pos] != x[m] )
                                        {    mismatch = True;
                                             break;    }    }
                                   if (mismatch) continue;
                                   vec<int> c;
                                   for ( int l = j + 1; l < k; l++ )
                                        c.push_back( x[l] );
                                   closures.push_back(c);    
                                        }    }    }    }    }
          if ( closures.nonempty( ) ) break;    }
     UniqueSort(closures);

     // Print.

     out << "\nclosures of " << ( pid >= 0 ? "+" : "" ) << pid << ":" << endl;
     for ( int j = 0; j < closures.isize( ); j++ )
     {    out << "[" << j+1 << "] " << printSeq(x1) << " | " 
               << printSeq( closures[j] ) << " | " << printSeq(x2) << endl;    }
     out << "\n" << ToStringAddCommas(opcounti) << " operations used"
          << " closing this pair\n" << endl;
     if (VERBOSE)
     {
          #pragma omp critical
          {    cout << out.str( );    }    }

     // Complete closures.

     for ( int j = 0; j < closures.isize( ); j++ )
     {    vec<int> x = x1;
          x.append( closures[j] );
          x.append(x2);
          closures[j] = x;    }    }

template<class VPI > void Closer( const HyperBasevectorX& hb, 
     const vec<int>& inv, ReadPathVecX& paths, VPI& paths_index, const vec<Bool>& dup, 
     const vec<Bool>& bad, vec<vec<int>>& all_closures,
     int64_t& opcount, const Bool GLOBAL, const int verbosity )
{
     // Define pair set.  Note that when we do all reads, we should 
     // canonicalize orientation, and push back both orientations at the end.

     if ( verbosity >= 1 ) 
     {    cout << Date( ) << ": define pair set, mem = " 
               << MemUsageGBString( ) << endl;    }
     vec<int64_t> ppids;
     DefinePairSet( hb, inv, paths, dup, bad, ppids, verbosity >= 1 );

     // Start main loop.

     if ( verbosity >= 1 )
     {    cout << Date( ) << ": start main loop, mem = " << MemUsageGBString( )
               << ", peakmem = " << PeakMemUsageGBString( ) << endl;
          cout << Date( ) << ": running 200 batches" << endl;    }
     double mclock = WallClockTime( );
     if ( !GLOBAL ) omp_set_num_threads(1);
     const int64_t batches = 200;
     int64_t batch = ( ppids.jsize( ) + batches - 1 ) / batches;
     int ndots = 0;
     #pragma omp parallel for schedule( dynamic, 1 )
     for ( int64_t bi = 0; bi < batches; bi++ )
     {    vec<vec<int>> closures, closuresi;
          int64_t opcountj = 0;
          int64_t start = bi*batch, stop = Min( (bi+1)*batch, ppids.jsize( ) );
          for ( int64_t pi = start; pi < stop; pi++ )
          {    ostringstream out;

               // Close pair.

               int64_t opcounti = 0;
               vec<int> x1, x2;
               // out << "x: " << printSeq(x1) << " .. " << printSeq(x2) << endl;
               const int64_t pid = ppids[pi];
               const int64_t id1 = 2*pid, id2 = 2*pid+1;
               ReadPath p1, p2; paths.unzip(p1,hb,id1); paths.unzip(p2,hb,id2);
               for ( int i = 0; i < (int) p1.size( ); i++ )
                    x1.push_back( p1[i] );
               for ( int i = (int) p2.size( ) - 1; i >= 0; i-- )
                    x2.push_back( inv[ p2[i] ] );
 
               closures.clear( );
               ClosePair( x1, x2, pid, hb, inv, paths, paths_index, 
                    opcounti, closures, out, verbosity >= 2 );
               opcountj += opcounti;

               // Save.

               closuresi.append(closures);    }

          // Save.

          #pragma omp critical
          {    all_closures.append(closuresi);
               opcount += opcountj;    }

          // Note progress.

          if ( verbosity >= 1 )
          {
               #pragma omp critical
               {    cout << ".";
                    ndots++;
                    if ( ndots > 0 && ndots % 50 == 0 ) cout << "\n";
                    else if ( ndots > 0 && ndots % 10 == 0 ) cout << " ";
                    flush(cout);    }    }    }
     if ( verbosity >= 1 )
     {    cout << Date( ) << ": main loop done, " << TimeSince(mclock)
               << " used, mem = " << MemUsageGBString( ) << ", peakmem = " 
               << PeakMemUsageGBString( ) << endl;    }    }

template void Closer( const HyperBasevectorX& hb, 
     const vec<int>& inv, ReadPathVecX& paths, 
     MasterVec<ULongVec>& paths_index, const vec<Bool>& dup, 
     const vec<Bool>& bad, vec<vec<int>>& all_closures,
     int64_t& opcount, const Bool GLOBAL, const int verbosity );

template void Closer( const HyperBasevectorX& hb, 
     const vec<int>& inv, ReadPathVecX & paths, 
     VirtualMasterVec<ULongVec>& paths_index, const vec<Bool>& dup, 
     const vec<Bool>& bad, vec<vec<int>>& all_closures,
     int64_t& opcount, const Bool GLOBAL, const int verbosity );
