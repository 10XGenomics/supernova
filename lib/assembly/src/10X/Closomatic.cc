// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "10X/Closomatic.h"
#include "CoreTools.h"
#include "feudal/PQVec.h"
#include "kmers/KmerRecord.h"
#include "math/Functions.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "paths/long/ReadStack.h"
#include "10X/DfTools.h"
#include "feudal/SubsetMasterVec.h"

void FindEdgePairs( const HyperBasevectorX& hb, const vec<int>& inv,
     MasterVec<ReadPath>& paths, String pi_file,
     const vec<Bool>& bad, vec< pair<int,int> >& pairs, 
     const vec<DataSet>& datasets, const vec<int32_t>& bc, const Bool one_good )
{    
     // Start.

     cout << Date( ) << ": start looking for edge pairs" << endl;
     double bclock = WallClockTime( );
     // read in paths_index in batches to save memory
     const int batch = 100000;
     
     // Create a vector of integers, one for each read, such that "having two"
     // is enough.  Note that the way we're treating barcode zero doesn't 
     // necessarily make sense.

     vec<int64_t> bid( paths.size( ) );
     #pragma omp parallel for
     for ( int64_t id = 0; id < (int64_t) paths.size( ); id++ )
     {    int di;
          for ( di = 0; di < datasets.isize( ); di++ )
               if ( id < datasets[di].start ) break;
          const ReadDataType& dtype = datasets[di-1].dt;
          if ( dtype == ReadDataType::BAR_10X ) 
               bid[id] = (int64_t) paths.size( ) + bc[id];
          else if ( dtype == ReadDataType::UNBAR_10X ) bid[id] = paths.size( );
          else if ( dtype == ReadDataType::PCR_FREE ) bid[id] = id;

          // Probably not what we want:

          else if ( dtype == ReadDataType::PCR ) bid[id] = id;    }

     // Heuristics.

     const int MAX_DIST_TO_END = 120;
     const int MIN_LANDING = 100;

     for (int start = 0; start < hb.E( ); start += batch)
     {
     const int stop = Min( hb.E( ), start + batch );
     VecULongVec paths_index_part;
     paths_index_part.ReadRange( pi_file, start, stop);
     for ( int e1 = start; e1 < stop; e1++ )
     {    
          // See if we're close enough to a sink.

          int v = hb.ToRight(e1);
          Bool OK = True;
          for ( int j = 0; j < (int) hb.From(v).size( ); j++ )
          {    int w = hb.From(v)[j];
               if ( hb.From(w).size( ) > 0 || hb.To(w).size( ) > 1 ) 
               {    OK = False;
                    break;    }
               int e = hb.IFrom( v, j );
               if ( hb.Kmers(e) > MAX_DIST_TO_END )
               {    OK = False;
                    break;    }    }
          if ( !OK ) continue;

          // cout << "\n" << Date( ) << ": looking at edge " << e1 << endl;
          vec< pair<int,int> > e2s;
          ULongVec & pi = paths_index_part[e1-start];
          for ( int i = 0; i < (int) pi.size( ); i++ )
          {    int64_t id1 = pi[i];
               const ReadPath& p1 = paths[id1];
               Bool found = False;
               for ( int j = 0; j < (int) p1.size( ); j++ )
                    if ( p1[j] == e1 ) found = True;
               if (found)
               {    int64_t id2 = ( id1 % 2 == 0 ? id1+1 : id1-1 );
                    const ReadPath& p2 = paths[id2];
                    if ( p2.size( ) > 0 )
                    {    int e2 = inv[ p2.back( ) ];
                         if ( e2 != e1 ) e2s.push( e2, bid[id1] );    }    }    }
          UniqueSort(e2s);
          for ( int i = 0; i < e2s.isize( ); i++ )
          {    int j;
               for ( j = i + 1; j < e2s.isize( ); j++ )
                    if ( e2s[j].first != e2s[i].first ) break;
               int e2 = e2s[i].first;
               if ( j - i >= 2 )
               {    
                    // See if we're close enough to a source.

                    int v = hb.ToLeft(e2);
                    Bool OK = True;
                    for ( int j = 0; j < (int) hb.To(v).size( ); j++ )
                    {    int w = hb.To(v)[j];
                         if ( hb.To(w).size( ) > 0 || hb.From(w).size( ) > 1 ) 
                         {    OK = False;
                              break;    }
                         int e = hb.ITo( v, j );
                         if ( hb.Kmers(e) > MAX_DIST_TO_END )
                         {    OK = False;
                              break;    }    }

                    // Save pair.

                    if ( one_good || OK )
                    {    pairs.push( e1, e2 );
                         // cout << "--> " << e2 << "[" << j-i << "]" << endl;    
                              }    }
               i = j - 1;   }   }
     }
     cout << Date( ) << ": done, time used = " << TimeSince(bclock) << endl;
     PRINT( pairs.size( ) );

     // Now look for cases we missed.

     cout << Date( ) << ": looking for one-sided cases" << endl;
     Sort(pairs);
     vec<Bool> seen( hb.E( ), False );
     for ( int i = 0; i < pairs.isize( ); i++ )
          seen[ pairs[i].first ] = True;
     
     for (int start = 0; start < hb.E( ); start += batch)
     {
     const int stop = Min( hb.E( ), start + batch );
     VecULongVec paths_index_part;
     paths_index_part.ReadRange( pi_file, start, stop); 
     for ( int e1 = start; e1 < stop; e1++ )
     {    if ( seen[e1] ) continue;

          // See if we're close enough to a sink.

          int v = hb.ToRight(e1);
          Bool OK = True;
          for ( int j = 0; j < (int) hb.From(v).size( ); j++ )
          {    int w = hb.From(v)[j];
               if ( hb.From(w).size( ) > 0 || hb.To(w).size( ) > 1 ) 
               {    OK = False;
                    break;    }
               int e = hb.IFrom( v, j );
               if ( hb.Kmers(e) > MAX_DIST_TO_END )
               {    OK = False;
                    break;    }    }
          if ( !OK ) continue;

          vec< pair<int,int> > e2s;
          ULongVec & pi = paths_index_part[e1-start];
          for ( int i = 0; i < (int) pi.size( ); i++ )
          {    int64_t id1 = pi[i];
               const ReadPath& p1 = paths[id1];
               Bool found = False;
               for ( int j = 0; j < (int) p1.size( ); j++ )
                    if ( p1[j] == e1 ) found = True;
               if (found)
               {    int64_t id2 = ( id1 % 2 == 0 ? id1+1 : id1-1 );
                    const ReadPath& p2 = paths[id2];
                    if ( p2.size( ) > 0 )
                    {    int e2 = inv[ p2.back( ) ];
                         if ( e2 != e1 ) e2s.push( e2, bid[id1] );    }    }    }
          UniqueSort(e2s);
          for ( int i = 0; i < e2s.isize( ); i++ )
          {    int j;
               for ( j = i + 1; j < e2s.isize( ); j++ )
                    if ( e2s[j].first != e2s[i].first ) break;
               int e2 = e2s[i].first;
               if ( j - i >= 2 && hb.Kmers(e2) >= MIN_LANDING
                    && e2 != e1 && hb.ToRight(e1) != hb.ToLeft(e2) )
               {    pairs.push( e1, e2 );    }
               i = j - 1;    }    }
     }
     // Add in pairs created by a different method.  The idea of the new 
     // approach is that if we cannot find a long extension of an edge, then
     // probably there is a "gap" at the edge (even if edges follow it in the
     // graph).  Here we require that the long extension be obtained using only
     // reads incident upon the edge.
     //
     // To do:
     // 1. Lower MIN_CAND and MIN_RIGHT.
     // 2. Exclude pairs that are closable in the graph.
     // 3. Eliminate the earlier approaches.

     VirtualMasterVec <ULongVec> vpaths_index ( pi_file );

     #pragma omp parallel for firstprivate(vpaths_index)
     for ( int e = 0; e < hb.E( ); e++ )
     {
          // Heuristics.

          const int MIN_CAND = hb.K( ) + 1;
          const int GOOD_EXT = 100;
          const int MIN_RIGHT = 40; // tried lowering to 25, not really better

          // Only start from long edges.

          if ( hb.Kmers(e) < MIN_CAND ) continue;

          // Find all paths visible from e and which should be to the right of e.

          int re = inv[e];
          vec<vec<int>> X;
          vec<int64_t> bs;
          ULongVec pi_e, pi_re;
//          #pragma omp critical
          {
               pi_e = vpaths_index[e];
               pi_re = vpaths_index[re];
          }
          for ( int i = 0; i < (int) pi_e.size( ); i++ )
          {    int64_t id1 = pi_e[i];
               // if ( dup[id1/2] || bad[id1/2] ) continue;
               bs.push_back( bid[id1] );
               const ReadPath& p1 = paths[id1];
               for ( int j = 0; j < (int) p1.size( ); j++ )
               {    if ( p1[j] == e )
                    {    vec<int> x;
                         for ( int l = j; l < (int) p1.size( ); l++ )
                              x.push_back( p1[l] );
                         X.push_back(x);    }    }
               int64_t id2 = ( id1 % 2 == 0 ? id1+1 : id1-1 );
               const ReadPath& p2 = paths[id2];
               if ( p2.size( ) > 0 )
               {    vec<int> x;
                    for ( int l = (int) p2.size( ) - 1; l >= 0; l-- )
                         x.push_back( inv[ p2[l] ] );
                    X.push_back(x);    }    }
          for ( int i = 0; i < (int) pi_re.size( ); i++ )
          {    int64_t id = pi_re[i];
               bs.push_back( bid[id] );
               // if ( dup[id/2] || bad[id/2] ) continue;
               const ReadPath& p = paths[id];
               for ( int j = 0; j < (int) p.size( ); j++ )
               {    if ( p[j] == re )
                    {    vec<int> x;
                         for ( int l = j; l >= 0; l-- )
                              x.push_back( inv[ p[l] ] );
                         X.push_back(x);    }    }    }
          UniqueSort(X);

          // Require at least two barcodes supporting the edge.

          UniqueSort(bs);
          if ( bs.size( ) >= 2 )
          {
               // See if we can find a path extending right from e by GOOD_EXT
               // kmers.  Seems insanely inefficient.

               Bool extended = False;
               vec<vec<int>> exts;
               for ( int i = 0; i < X.isize( ); i++ )
                    if ( X[i][0] == e ) exts.push_back( X[i] );
               while(1)
               {    vec<vec<int>> exts2;
                    for ( int i = 0; i < exts.isize( ); i++ )
                    {    const vec<int>& x = exts[i];
                         int n = 0;
                         for ( int j = 1; j < x.isize( ); j++ )
                              n += hb.Kmers( x[j] );
                         if ( n >= GOOD_EXT )
                         {    extended = True;
                              break;    }    }
                    if (extended) break;
                    for ( int i = 0; i < exts.isize( ); i++ )
                    {    const vec<int>& x = exts[i];
                         int f = x.back( );
                         for ( int j = 0; j < X.isize( ); j++ )
                         {    const vec<int>& y = X[j];
                              for ( int l = 0; l < y.isize( ) - 1; l++ )
                              {    if ( y[l] != f ) continue;
                                   Bool mismatch = False;
                                   for ( int m = 0; m < y.isize( ); m++ )
                                   {    int n = m + x.isize( ) - 1 - l;
                                        if ( n < 0 || n >= x.isize( ) ) continue;
                                        if ( x[n] != y[m] )
                                        {    mismatch = True;
                                             break;    }    }
                                   if (mismatch) continue;
                                   vec<int> x2(x);
                                   for ( int z = l + 1; z < y.isize( ); z++ )
                                        x2.push_back( y[z] );
                                   exts2.push_back(x2);    }    }    }
                    if ( exts2.empty( ) ) break;
                    UniqueSort(exts2);
                    exts = exts2;    }
               if ( !extended ) 
               {    
                    // Find possible rights.

                    vec<int> too_easy, rights;
                    vec< pair<int,int64_t> > can;
                    for ( int i = 0; i < (int) pi_e.size( ); i++ )
                    {    int64_t id1 = pi_e[i];
                         if ( bad[id1/2] ) continue;
                         const ReadPath& p1 = paths[id1];
                         for ( int j = 0; j < (int) p1.size( ); j++ )
                         {    if ( p1[j] == e )
                              {    for ( int l = j+1; l < (int) p1.size( ); l++ )
                                   {    int f = p1[l];
                                        if ( hb.Kmers(f) >= MIN_RIGHT 
                                             && f != e && f != inv[e] )
                                        {    too_easy.push_back(f);    }    }    }
                                                  }
                         int64_t id2 = ( id1 % 2 == 0 ? id1+1 : id1-1 );
                         const ReadPath& p2 = paths[id2];
                         if ( p2.size( ) > 0 )
                         {    for ( int l = (int) p2.size( ) - 1; l >= 0; l-- )
                              {    int f = inv[ p2[l] ];
                                   if ( hb.Kmers(f) >= MIN_RIGHT 
                                        && f != e && f != inv[e] )
                                   {    can.push( f, bid[id1] );    }    }    }    }
                    for ( int i = 0; i < (int) pi_re.size( ); i++ )
                    {    int64_t id = pi_re[i];
                         if ( bad[id/2] ) continue;
                         const ReadPath& p = paths[id];
                         for ( int j = 0; j < (int) p.size( ); j++ )
                         {    if ( p[j] == re )
                              {    for ( int l = j; l >= 0; l-- )
                                   {    int f = inv[ p[l] ];
                                        if ( hb.Kmers(f) >= MIN_RIGHT 
                                             && f != e && f != inv[e] )
                                        {    can.push( f, 
                                                  bid[id] );    }    }    }    }    }
                    UniqueSort(too_easy), UniqueSort(can);
                    for ( int i = 0; i < can.isize( ); i++ )
                    {    int j;
                         for ( j = i + 1; j < can.isize( ); j++ )
                              if ( can[j].first != can[i].first ) break;
                         if ( j - i >= 2 && !BinMember( too_easy, can[i].first ) )
                              rights.push_back( can[i].first );
                         i = j - 1;    }
                    #pragma omp critical
                    {    for ( int i = 0; i < rights.isize( ); i++ )
                              pairs.push( e, rights[i] );    }    }    }    }

     // Done.

     UniqueSort(pairs);
     cout << Date( ) << ": done, time used = " << TimeSince(bclock) << endl;
     PRINT( pairs.size( ) );
     
     // hash pairs
/*     std::string hpairs;
     hpairs.reserve( 4*pairs.size() );
     for ( auto & p : pairs ) {
          hpairs += std::to_string(p.first);
          hpairs += std::to_string(p.second);
     }
     int64_t pair_hash = std::hash<std::string>()(hpairs);
     PRINT(pair_hash);   */ 
     }

template< class VB, class VQ, class VP, class VPI >
void CloseGap( const HyperBasevectorX& hb, const vec<int>& inv, 
     VB& bases, VQ& quals, VP& paths, VPI& paths_index,
     const int e1, const int e2, vec<basevector>& closures, Bool verbose, 
     const int pi, const int max_width )
{    closures.clear( );
     ostringstream out;
     if (verbose) out << "\n";
     if (verbose) PRINT2_TO( out, e1, e2 );

     int v = hb.ToRight(e1), w = hb.ToLeft(e2);

     // Stuff.

     vec<int> es = { e1, inv[e2] };
     vec<readstack> stacks(2);
     for ( int ei = 0; ei < 2; ei++ )
     {    int e = es[ei];

          int ne = hb.Bases(e);
          if (verbose) out << "\n";
          if (verbose) PRINT4_TO( out, pi, e, ne, inv[e] );

          // Find support.
     
          vec< triple<int64_t,int,Bool> > locs;
          for ( int pass = 1; pass <= 2; pass++ )
          {    const int f = ( pass == 1 ? e : inv[e] );
               for ( int i = 0; i < (int) paths_index[f].size( ); i++ )
               {    int64_t id = paths_index[f][i];
                    const ReadPath& p = paths[id];
                    for ( int j = 0; j < (int) p.size( ); j++ )
                    {    if ( p[j] == f )
                         {    int pos = p.getOffset( );
                              for ( int l = 0; l < j; l++ )
                                   pos -= hb.Kmers( p[l] );
                              locs.push( id, pos, pass == 1 );    }    }    }    }

          // Build stack.
     
          int M = 0;
          for ( int i = 0; i < locs.isize( ); i++ )
          {    int64_t id = locs[i].first;
               int pos = locs[i].second;
               Bool fw = locs[i].third;
               if (fw) M = Max( M, pos + bases[id].isize( ) );
               else M = Max( M, ne - pos );    }
          readstack stack( locs.size( ), M );
          for ( int i = 0; i < locs.isize( ); i++ )
          {    int64_t id = locs[i].first;
               qualvector q;
               quals[id].unpack(&q);
               int pos = locs[i].second;
               Bool fw = locs[i].third;
               stack.SetId( i, id );
               stack.SetLen( i, bases[id].size( ) );
               if (fw)
               {    stack.SetOffset( i, pos );
                    for ( int j = 0; j < bases[id].isize( ); j++ )
                    {    int p = pos + j;
                         if ( p >= 0 /* && p < ne */ )
                         {    stack.SetBase( i, p, bases[id][j] );
                              stack.SetQual( i, p, q[j] );    }    }    }
               else 
               {    stack.SetOffset( i, ne - pos - bases[id].isize() );
                    for ( int j = 0; j < bases[id].isize( ); j++ )
                    {    int p = pos + j;
                         // if ( p >= 0 /* && p < ne */ )
                         if ( ne - p - 1 >= 0 )
                         {    stack.SetBase( i, ne - p - 1, 3 - bases[id][j] );
                              stack.SetQual( 
                                   i, ne - p - 1, q[j] );    }    }    }    }

          // Trim stack.

          if ( stack.Cols( ) > max_width )
               stack.Trim( stack.Cols( ) - max_width, stack.Cols( ) );

          // Print stack.

          if ( e == inv[e2] ) stack.Reverse( );
          stacks[ei] = stack;
          if (verbose) stack.Print(out);    }    
     if (verbose) cout << out.str( );

     // Merge stacks.

     int verbosity = 0; // could be 1 or 2
     if (verbose) verbosity = 3;
     int delta_mis = 0;
     vec<int> o = GetOffsets1( stacks[0], stacks[1], verbosity, delta_mis );
     ostringstream out2;
     if (verbose) out2 << "\noffsets = " << printSeq(o) << endl;
     if ( o.size( ) <= 2 )
     {    for ( int j = 0; j < o.size( ); j++ )
          {    readstack s( stacks[0] );
               s.Merge( stacks[1], o[j] );
               basevector f;
               qualvector fragq;
               s.Consensus1( f, fragq );
               if (verbose) out2 << endl;
               if (verbose) f.Print( out2, "consensus" );    
               closures.push_back(f);    }     }
     if (verbose) cout << out2.str( );    }

template void CloseGap( const HyperBasevectorX& hb, const vec<int>& inv,
     MasterVec<basevector>& bases, MasterVec<PQVec>& quals,
     MasterVec<ReadPath>& paths, MasterVec<ULongVec>& paths_index,
     const int e1, const int e2, vec<basevector>& closure, Bool verbose, 
     const int pi, const int max_width );

template void CloseGap( const HyperBasevectorX& hb, const vec<int>& inv,
     MasterVec<basevector>& bases, SubsetMasterVec<PQVec>& quals,
     SubsetMasterVec<ReadPath>& paths, SubsetMasterVec<ULongVec>& paths_index,
     const int e1, const int e2, vec<basevector>& closure, Bool verbose, 
     const int pi, const int max_width );

template void CloseGap( const HyperBasevectorX& hb, const vec<int>& inv,
     VirtualMasterVec<basevector>& bases, VirtualMasterVec<PQVec>& quals,
     VirtualMasterVec<ReadPath>& paths, VirtualMasterVec<ULongVec>& paths_index,
     const int e1, const int e2, vec<basevector>& closure, Bool verbose, 
     const int pi, const int max_width );

template< class VB, class VQ, class VP, class VPI >
void CloseGap2( const HyperBasevectorX& hb, const vec<int>& inv, 
     VB& bases, VQ& quals, VP& paths, VPI& paths_index,
     const int e1, const int e2, vec<basevector>& closures, Bool verbose, 
     const int pi, const int max_width )
{    
     // Setup.

     closures.clear( );
     ostringstream out;
     if (verbose) out << "\n";
     if (verbose) PRINT2_TO( out, e1, e2 );

     // Heuristics.

     const int MIN_EDGE_COUNT = 5;
     const int L = 32;

     // Find read locations and preceeding edges.

     int v = hb.ToRight(e1), w = hb.ToLeft(e2);
     vec<int> es = { e1, inv[e2] };
     vec< vec< triple<int64_t,int,Bool> > > locs(2);
     vec< vec< pair< kmer<L>, int > > > bod(2);
     for ( int ei = 0; ei < 2; ei++ )
     {    int e = es[ei];
          vec< pair<int,int> > prec;
          // prec.push( e, 0 ); // doesn't seem to help
          for ( int pass = 1; pass <= 2; pass++ )
          {    const int f = ( pass == 1 ? e : inv[e] );
               for ( int i = 0; i < (int) paths_index[f].size( ); i++ )
               {    int64_t id = paths_index[f][i];
                    const ReadPath& p = paths[id];
                    int n = p.size( );
                    for ( int j = 0; j < n; j++ )
                    {    if ( p[j] == f )
                         {    if ( pass == 1 )
                              {    int add = 0;
                                   for ( int l = j + 1; l < n; l++ )
                                   {    int g = p[l];
                                        add += hb.Kmers( p[l-1] );
                                        prec.push( g, add );    }    }
                              else
                              {    int add = 0;
                                   for ( int l = j - 1; l >= 0; l-- )
                                   {    int g = p[l];
                                        add += hb.Kmers( p[l+1] );
                                        prec.push( inv[g], add );    }    }
                              int pos = p.getOffset( );
                              for ( int l = 0; l < j; l++ )
                                   pos -= hb.Kmers( p[l] );
                              locs[ei].push( 
                                   id, pos, pass == 1 );    }    }    }    }    
          Sort(prec);
          if ( verbose ) cout << ( ei == 0 ? "LEFT:\n" : "RIGHT:\n" ) << endl;
          for ( int i = 0; i < prec.isize( ); i++ )
          {    int j;
               for ( j = i + 1; j < prec.isize( ); j++ )
                    if ( prec[j] != prec[i] ) break;
               if ( j - i >= MIN_EDGE_COUNT )
               {    int g = prec[i].first;
                    if (verbose) cout << "e = " << ( ei == 0 ? g : inv[g] ) << ", epos = "
                         << prec[i].second << " [" << j-i << "]" << endl;
                    const basevector& E = hb.O(g);
                    for ( int l = 0; l <= E.isize( ) - L; l++ )
                    {    kmer<L> x;
                         x.SetToSubOf( E, l );
                         int pos = prec[i].second + l;
                         bod[ei].push( x, pos );    }    }
               i = j - 1;    }    }

     // Try to add in unplaced partners.

     vec<vec<int64_t>> ids(2);
     for ( int ei = 0; ei < 2; ei++ )
     {    for ( int j = 0; j < locs[ei].isize( ); j++ )
          {    if ( !locs[ei][j].third ) continue;
               ids[ei].push_back( locs[ei][j].first );    }
          UniqueSort( ids[ei] ), UniqueSort( bod[ei] );    }
     for ( int ei = 0; ei < 2; ei++ )
     {    for ( int j = 0; j < ids[ei].isize( ); j++ )
          {    int64_t id1 = ids[ei][j];
               int64_t id2 = ( id1 % 2 == 0 ? id1+1 : id1-1 );
               if ( BinMember( ids[1-ei], id2 ) ) continue;
     
               // Now we know that id2 might belong in locs[1-ei], but we don't
               // know the position.
     
               if (verbose) cout << "ei = " << 1-ei << ", try to locate " << id2 << endl;
               const basevector& b = bases[id2];
               vec<int> finds;
               for ( int l = 0; l <= b.isize( ) - L; l++ )
               {    kmer<L> x;
                    x.SetToSubOf( b, l );
                    int64_t low = LowerBound1( bod[1-ei], x );
                    int64_t high = UpperBound1( bod[1-ei], x );
                    for ( int64_t m = low; m < high; m++ )
                         finds.push_back( bod[1-ei][m].second - l );    }
               UniqueSort(finds);
               if ( finds.nonempty( ) )
                    if ( verbose )cout << id2 << " found at " << printSeq(finds) << endl;
               if ( finds.solo( ) )
                    locs[1-ei].push( id2, finds[0], True );    }    }

     // Stuff.

     vec<readstack> stacks(2);
     for ( int ei = 0; ei < 2; ei++ )
     {    int e = es[ei];
          int ne = hb.Bases(e);
          if (verbose) out << "\n";
          if (verbose) PRINT4_TO( out, pi, e, ne, inv[e] );

          // Build stack.
     
          int M = 0;
          for ( int i = 0; i < locs[ei].isize( ); i++ )
          {    int64_t id = locs[ei][i].first;
               int pos = locs[ei][i].second;
               Bool fw = locs[ei][i].third;
               if (fw) M = Max( M, pos + bases[id].isize( ) );
               else M = Max( M, ne - pos );    }
          readstack stack( locs[ei].size( ), M );
          for ( int i = 0; i < locs[ei].isize( ); i++ )
          {    int64_t id = locs[ei][i].first;
               qualvector q;
               quals[id].unpack(&q);
               int pos = locs[ei][i].second;
               Bool fw = locs[ei][i].third;
               stack.SetId( i, id );
               stack.SetLen( i, bases[id].size( ) );
               if (fw)
               {    stack.SetOffset( i, pos );
                    for ( int j = 0; j < bases[id].isize( ); j++ )
                    {    int p = pos + j;
                         if ( p >= 0 /* && p < ne */ )
                         {    stack.SetBase( i, p, bases[id][j] );
                              stack.SetQual( i, p, q[j] );    }    }    }
               else 
               {    stack.SetOffset( i, ne - pos - bases[id].isize() );
                    for ( int j = 0; j < bases[id].isize( ); j++ )
                    {    int p = pos + j;
                         // if ( p >= 0 /* && p < ne */ )
                         if ( ne - p - 1 >= 0 )
                         {    stack.SetBase( i, ne - p - 1, 3 - bases[id][j] );
                              stack.SetQual( 
                                   i, ne - p - 1, q[j] );    }    }    }    }

          // Trim stack.

          if ( stack.Cols( ) > max_width )
               stack.Trim( stack.Cols( ) - max_width, stack.Cols( ) );

          // Print stack.

          if ( e == inv[e2] ) stack.Reverse( );
          stacks[ei] = stack;
          if (verbose) stack.Print(out);    }    
     if (verbose) cout << out.str( );

     // Merge stacks.

     int verbosity = 0; // could be 1 or 2
     if (verbose) verbosity = 3;
     int delta_mis = 0;
     vec<int> o = GetOffsets1( stacks[0], stacks[1], verbosity, delta_mis );
     ostringstream out2;
     if (verbose) out2 << "\noffsets = " << printSeq(o) << endl;
     if ( o.size( ) <= 2 )
     {    for ( int j = 0; j < o.size( ); j++ )
          {    readstack s( stacks[0] );
               s.Merge( stacks[1], o[j] );
               basevector f;
               qualvector fragq;
               s.Consensus1( f, fragq );
               if (verbose) out2 << endl;
               if (verbose) f.Print( out2, "consensus" );    
               closures.push_back(f);    }     }
     if (verbose) cout << out2.str( );    }

template void CloseGap2( const HyperBasevectorX& hb, const vec<int>& inv,
     MasterVec<basevector>& bases, MasterVec<PQVec>& quals,
     MasterVec<ReadPath>& paths, MasterVec<ULongVec>& paths_index,
     const int e1, const int e2, vec<basevector>& closure, Bool verbose,
     const int pi, const int max_width );

template void CloseGap2( const HyperBasevectorX& hb, const vec<int>& inv,
     MasterVec<basevector>& bases, SubsetMasterVec<PQVec>& quals,
     SubsetMasterVec<ReadPath>& paths, SubsetMasterVec<ULongVec>& paths_index,
     const int e1, const int e2, vec<basevector>& closure, Bool verbose,
     const int pi, const int max_width );

template void CloseGap2( const HyperBasevectorX& hb, const vec<int>& inv,
     VirtualMasterVec<basevector>& bases, VirtualMasterVec<PQVec>& quals,
     VirtualMasterVec<ReadPath>& paths, VirtualMasterVec<ULongVec>& paths_index,
     const int e1, const int e2, vec<basevector>& closure, Bool verbose,
     const int pi, const int max_width );
