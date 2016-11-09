// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

// Build a local assembly, trying to walk from the end of one line to another.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"
#include "ParallelVecUtilities.h"
#include "graph/Digraph.h"
#include "graph/DigraphTemplate.h"
#include "kmers/KmerRecord.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "paths/long/large/Lines.h"
#include "10X/ClosuresToGraph.h"
#include "10X/DfTools.h"
#include "10X/Gap.h"
#include "10X/Heuristics.h"
#include "10X/Scaffold.h"
#include "10X/SecretOps.h"
#include "paths/long/large/GapToyTools.h"
#include "10X/BuildLocal.h"
#include "10X/Super.h"
#include "10X/PlaceReads.h"

// VP = (Virtual)MasterVec<ReadPath>
// VA = (Virtual)MasterVec< SerfVec<triple<int,int,int> > >
// VE = (Virtual)MasterVec<IntVec>

// Note that this creates some global data structures, that should actually be
// passed as arguments.

template<class VA, class VE> void BuildLocal1(

     // Edges to walk between:

     const int s1, const vec<int>& s2s,

     // Global inputs:

     const HyperBasevectorX& hb, const vec<int>& inv, const vec<int64_t>& bci, 
     vec<Bool>& dup, vec<Bool>& bad, ReadPathVecX& paths, VA& alignsb, VE& ebcx, 
     const digraphE<vec<int>>& D, const vec<int>& dinv, 
     const vec<vec<vec<vec<int>>>>& dlines, const ReadPathVec& dpaths, 
     const VecULongVec& dpaths_index, const vec<Bool>& internal,
     const vec<int>& to_left, const vec<int>& mult, const vec<int>& lens,

     // Control:

     const Bool use_rights,

     // Local outputs:

     vec<int32_t>& bcl, HyperBasevector& hbl, HyperBasevectorX& hbxl,
     vec<int>& kmersl, vec<int>& invl, ReadPathVec& pathsl,
     digraphE<vec<int>>& Dl, vec<int>& dinvl, ReadPathVec& dpathsl,
     VecULongVec& dpaths_indexl, vec<Bool>& dupl,
     MasterVec< SerfVec<triple<int,int,int> > >& alignsbl,
     vec<vec<vec<vec<int>>>>&  dlinesl, String& link_report, 
     digraphE<vec<int>>& Dlp,

     // Logging etc.:

     Bool verbose, Bool results_only,
     double& c1, double& c2a, double& c2b, double& c3, double& c4,
     double& c5, double& c6a1, double& c6a2, double& c6a3, double& c6b,
     double& c6c, double& c7,

     // Control:

     const int MAX_READS, const Bool SINGLE )
{
     double clock1a = WallClockTime( );
     #pragma omp critical
     {    c1 -= WallClockTime( );    }
     ForceAssertGe( D.O(s1)[0], 0 );

     // Compute barcode set.

     if (verbose) cout << Date( ) << ": computing barcode set" << endl;
     vec<int> b;
     const int GRAB = 10000;
     const int MAX_BARCODES = 1000;
     for ( int pass = 1; pass <= 2; pass++ )
     {    int MIN_KMERS = ( pass == 1 ? 1 : 10 );
          b.clear( );
          GetBarcodes( s1, GRAB, MAX_BARCODES, MIN_KMERS, 
               hb, D, to_left, lens, mult, ebcx, b );
          if (use_rights)
          {    for ( int l = 0; l < s2s.isize( ); l++ )
               {    int s2 = s2s[l];
                    GetBarcodes( dinv[s2], GRAB, MAX_BARCODES, MIN_KMERS, 
                         hb, D, to_left, lens, mult, ebcx, b );    }    }
          if ( b.isize( ) <= MAX_BARCODES ) break;    }
     #pragma omp critical
     {    c1 += WallClockTime( );    }
     if ( b.isize( ) > MAX_BARCODES ) return;
     if (verbose) cout << Date( ) << ": used " << b.size( ) << " barcodes" << endl;
     ForceAssert( !BinMember( b, 0 ) ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

     // Make sure we're under MAX_READS.

     if ( MAX_READS >= 0 )
     {    vec<int64_t> bs( b.size( ) );
          for ( int i = 0; i < b.isize( ); i++ )
               bs[i] = bci[ b[i] + 1 ] - bci[ b[i] ];
          SortSync( bs, b );
          int64_t sum = Sum(bs);
          while ( sum > MAX_READS )
          {    sum -= bs.back( );
               bs.pop_back( );
               b.pop_back( );    }
          Sort(b);    }

     // Find edges that are supported by at least two barcodes.  Results were
     // slightly worse when we used a threshold of three.  Add in the ends.

     #pragma omp critical
     {    c2a -= WallClockTime( );    }
     if (verbose) cout << Date( ) << ": finding strong edges" << endl;
     const int MIN_BC = 2;
     vec<int> es;
     vec< pair<int,int> > esb;
     {    for ( int i = 0; i < b.isize( ); i++ )
          for ( int64_t id = bci[ b[i] ]; id < bci[ b[i] + 1 ]; id++ )
          {    // Dup test turned on makes results slightly worse.
               if ( dup[id/2] ) continue;
               if ( internal[id] ) continue;
               ReadPath p;
               paths.unzip(p,hb,id);
               for ( int j = 0; j < (int) p.size( ); j++ )
                    esb.push( Min( p[j], inv[ p[j] ] ), b[i] );    }
          UniqueSort(esb);
          for ( int i = 0; i < esb.isize( ); i++ )
          {    int j;
               for ( j = i + 1; j < esb.isize( ); j++ )
                    if ( esb[j].first != esb[i].first ) break;
               if ( j - i >= MIN_BC ) 
                    es.push_back( esb[i].first, inv[ esb[i].first ] );
               i = j - 1;    }
          for ( int j = 0; j < D.O(s1).isize( ); j++ )
               es.push_back( D.O(s1)[j], inv[ D.O(s1)[j] ] );
          for ( int l = 0; l < s2s.isize( ); l++ )
          {    int s2 = s2s[l];

               // Note working around weird problem with lines start with a gap!!!!!

               if ( D.O(s2)[0] < 0 ) continue;

               for ( int j = 0; j < D.O(s2).isize( ); j++ )
                    es.push_back( D.O(s2)[j], inv[ D.O(s2)[j] ] );    }    }
     UniqueSort(es);
     #pragma omp critical
     {    c2a += WallClockTime( );
          c2b -= WallClockTime( );    }

     // Form local HyperBasevector.

     if (verbose) cout << Date( ) << ": forming local HyperBasevector" << endl;
     ( (digraphE<basevector>&) hbl ).Initialize( 
          digraphE<basevector>::COMPLETE_SUBGRAPH_EDGES, hb, es );
     hbl.SetK( hb.K( ) );
     // hbl.TestValid( ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     #pragma omp critical
     {    c2b += WallClockTime( );
          c3 -= WallClockTime( );    }

     // Form local inversion.

     if (verbose) cout << Date( ) << ": forming other local structures" << endl;
     invl.resize( es.size( ) );
     // #pragma omp parallel for
     for ( int i = 0; i < es.isize( ); i++ )
          invl[i] = BinPosition( es, inv[ es[i] ] );
     // for ( int j = 0; j < invl.isize( ); j++ ) // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     //      ForceAssertGe( invl[j], 0 ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

     // Make list of read ids, and their barcodes.  We also include reads that
     // look off the end of s1.

     vec< pair<int64_t,int> > idsb;
     for ( int i = 0; i < b.isize( ); i++ )
     for ( int64_t id = bci[ b[i] ]; id < bci[ b[i] + 1 ]; id++ )
     {    // Dup test makes results slightly worse.
          if ( dup[id/2] ) continue;
          if ( internal[id] ) continue;
          idsb.push( id, b[i] );    }
     int rs1 = dinv[s1], ns1 = 0;
     for ( int i = 0; i < D.O(s1).isize( ); i++ )
          ns1 += hb.Kmers( D.O(s1)[i] );
     for ( int pass = 1; pass <= 2; pass++ )
     {    int s = ( pass == 1 ? s1 : rs1 );
          for ( int i = 0; i < (int) dpaths_index[s].size( ); i++ )
          {    int64_t id = dpaths_index[s][i];
               if ( id >= bci[1] ) continue; // require barcode zero
               int64_t pid = id/2;
               if ( id < (int) dpaths_index[s].size( ) - 1 &&
                    dpaths_index[s][i+1] == 2*pid+1 )
               {    continue;    } // ignore placed pairs
               const ReadPath& p = dpaths[id];
               if ( pass == 1 && ns1 - p.getOffset( ) > 1000 ) continue;
               if ( pass == 2 && p.getOffset( ) > 1000 ) continue;
               idsb.push( 2*pid, 0 );
               idsb.push( 2*pid+1, 0 );    }    }
     UniqueSort(idsb);

     // Form local paths, truncating as needed.  Note that offsets should be 
     // adjusted, but they're not.

     vec<Bool> badl;
     {    vec<vec<int>> X;
          vec<int> offsets;
          for ( int i = 0; i < idsb.isize( ); i++ )
          {    int64_t id = idsb[i].first;
               int b = idsb[i].second;
               int offs;
               vec<int> x;
               paths.unzip(x,offs,hb,id);
               offsets.push_back( offs);
               if ( id % 2 == 0 )
               {    dupl.push_back( dup[id/2] );
                    badl.push_back( bad[id/2] );    }
               bcl.push_back(b);

               vec<vec<int>> xs;
               for ( int i = 0; i < x.isize( ); i++ )
               {    if ( !BinMember( es, x[i] ) ) continue;
                    int j;
                    for ( j = i + 1; j < x.isize( ); j++ )
                         if ( !BinMember( es, x[j] ) ) break;
                    vec<int> y;
                    for ( int k = i; k < j; k++ )
                         y.push_back( x[k] );
                    xs.push_back(y);
                    i = j - 1;    }
               if ( xs.empty( ) ) x.clear( );
               else
               {    vec<int> lens;
                    for ( int i = 0; i < xs.isize( ); i++ )
                    {    int sum = 0;
                         const vec<int>& y = xs[i];
                         for ( int j = 0; j < y.isize( ); j++ )
                              sum += hb.Kmers( y[j] );
                         lens.push_back(sum);    }
                    int M = Max(lens);
                    for ( int j = 0; j < xs.isize( ); j++ )
                    {    if ( lens[j] == M )
                         {    x = xs[j];
                              break;    }    }    }

               /*
               int start = 0, stop = x.size( );
               while( start < stop && !BinMember( es, x[start] ) ) start++;
               while( start < stop && !BinMember( es, x[stop-1] ) ) stop--;
               x.SetToSubOf( x, start, stop - start );
               for ( int j = 0; j < x.isize( ); j++ )
               {    if ( !BinMember( es, x[j] ) )
                    {    x.clear( );
                         break;    }    }
               */

               /*
               cout << "lpid = " << X.size( ) / 2 << ", lid = " << X.size( ) 
                    << ", id = " << id << ", pid = " << id/2 
                    << ", x = " << printSeq(x) << endl;
               */

               X.push_back(x);    }
          for ( int i = 0; i < X.isize( ); i += 2 )
          {    const vec<int> &x1 = X[i], &x2 = X[i+1];
               // if ( x1.empty( ) && x2.empty( ) ) continue;
               ReadPath p1, p2;
               if ( x1.nonempty( ) ) p1.setOffset( offsets[i] );
               if ( x2.nonempty( ) ) p2.setOffset( offsets[i+1] );
               SerfVec<int> y1( x1.size( ) ), y2( x2.size( ) );
               for ( int j = 0; j < (int) x1.isize( ); j++ )
                    y1[j] = BinPosition( es, x1[j] );
               for ( int j = 0; j < (int) x2.isize( ); j++ )
                    y2[j] = BinPosition( es, x2[j] );
               (SerfVec<int>&) p1 = y1;
               (SerfVec<int>&) p2 = y2;
               pathsl.push_back(p1), pathsl.push_back(p2);    }    }
     #pragma omp critical
     {    c3 += WallClockTime( );
          c4 -= WallClockTime( );    }
     if (verbose)
     {    cout << Date( ) << ": s1 = " << s1 << ", " << b.size( ) << " barcodes"
               << ", " << es.size( ) << " strong edges"
               << ", " << pathsl.size( ) << " paths"
               << ", BuildLocal1a time = " << TimeSince(clock1a) << endl;    }
     double clock1b = WallClockTime( );

     // Create a vector of integers, one for each read, such that "having two"
     // nonzero elements is enough.

     DataSet dl;
     dl.dt = ReadDataType::BAR_10X;
     dl.start = 0;
     vec<DataSet> datasetsl = {dl};
     if (verbose) cout << Date( ) << ": creating bidl" << endl;
     vec<int64_t> bidl( pathsl.size( ) );
     int nthreads = ( SINGLE ? 1 : omp_get_max_threads( ) );
     #pragma omp parallel for num_threads(nthreads)
     for ( int64_t id = 0; id < (int64_t) pathsl.size( ); id++ )
     {    int di;
          for ( di = 0; di < datasetsl.isize( ); di++ )
               if ( id < datasetsl[di].start ) break;
          const ReadDataType& dtype = datasetsl[di-1].dt;
          if ( dtype == ReadDataType::BAR_10X ) 
               bidl[id] = (int64_t) pathsl.size( ) + bcl[id] + 1;
          else if ( dtype == ReadDataType::UNBAR_10X ) bidl[id] = 0;
          else if ( dtype == ReadDataType::PCR_FREE ) bidl[id] = id + 1;

          // Probably not what we want:

          else if ( dtype == ReadDataType::PCR ) bidl[id] = id + 1;    }

     // Create temporary paths for s1 and s2.

     ReadPath ps1;
     ps1.setOffset(0);
     int nps = 1;
     for ( auto e : D.O(s1) ) ps1.push_back( BinPosition( es, e ) );
     pathsl.push_back(ps1);
     for ( int l = 0; l < s2s.isize( ); l++ )
     {    int s2 = s2s[l];
          if ( D.O(s2)[0] < 0 ) continue;
          nps++;
          ReadPath ps2;
          ps2.setOffset(0);
          for ( auto e : D.O(s2) ) ps2.push_back( BinPosition( es, e ) );
          pathsl.push_back(ps2);    }

     // Simplify the graph.  Start by removing unneeded vertices.  Note that we 
     // observe cases of circles where no unneeded vertices are removed, and there 
     // is no intersection with the involution.  Just an interesting observation, 
     // nothing to do with this module per se.
     //
     // We keep the original HyperBasevector, so we can translate back to it.

     if (verbose) cout << Date( ) << ": cleaning" << endl;
     HyperBasevector hbls(hbl);
     RemoveUnneededVertices2( hbl, invl, pathsl, False, True );
     CleanupCore( hbl, invl, pathsl );
     // NOTE EXPENSIVE CONVERSION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     hbxl = HyperBasevectorX(hbl);
     #pragma omp critical
     {    c4 += WallClockTime( );
          c5 -= WallClockTime( );    }

     // Map the HyperBasevector edges back to the original.

     if (verbose) cout << Date( ) << ": remapping" << endl;
     const int K = 48;
     ForceAssertEq( K, hbl.K( ) );
     vec<kmer<K>> starts( hbls.E( ) );
     // #pragma omp parallel for
     for ( int e = 0; e < hbls.E( ); e++ )
     {    starts[e].SetToSubOf( hbls.O(e), 0 );    }
     vec<int> sids( hbls.E( ), vec<int>::IDENTITY );
     if (SINGLE) SortSync( starts, sids );
     else ParallelSortSync( starts, sids );
     vec<vec<int>> backto( hbl.E( ) );
     vec<int> sto_right;
     hbls.ToRight(sto_right);
     // #pragma omp parallel for
     for ( int e = 0; e < hbl.E( ); e++ )
     {    kmer<K> y;
          y.SetToSubOf( hbl.O(e), 0 );
          int p = BinPosition( starts, y );
          ForceAssert( p >= 0 );
          vec<int>& x = backto[e];
          x = { sids[p] };
          int n = hbls.Bases( x[0] );
          ForceAssertLe( n, hbl.Bases(e) );
          while( n < hbl.Bases(e) )
          {    int v = sto_right[ x.back( ) ];
               for ( int j = 0; j < hbls.From(v).isize( ); j++ )
               {    int f = hbls.IFrom( v, j );
                    if ( hbls.O(f)[K-1] == hbl.O(e)[n] )
                    {    x.push_back(f);
                         n += hbls.Kmers(f);
                         break;    }    }    }
          if ( n != hbl.Bases(e) ) PRINT(e);
          ForceAssertEq( n, hbl.Bases(e) );    }

     // Take the s1 and s2 paths off pathsl.

     vec<vec<int>> ss;
     for ( int j = 1; j <= nps; j++ )
     {    vec<int> x;
          for ( auto e : pathsl[ pathsl.size( ) - j ] ) x.push_back(e);
          ss.push_back(x);    }
     pathsl.resize( pathsl.size( ) - nps );

     // Form paths index.

     if (verbose) cout << Date( ) << ": building path index" << endl;
     VecULongVec paths_indexl;
     invert( pathsl, paths_indexl, hbl.E( ) );
     #pragma omp critical
     {    c5 += WallClockTime( );
          c6a1 -= WallClockTime( );    }

     // Form closures and convert to graph.  Note that scaffolding does not use
     // the PCR-free reads, which may not make sense.

     if (verbose)
          cout << Date( ) << ": forming closures and converting to graph" << endl;
     vec<vec<int>> all_closuresl;
     // NOTE EXPENSIVE CONVERSION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ReadPathVecX pathslx;
     HyperBasevectorX hblx(hbl);
     pathslx.append(pathsl,hblx);
     MakeClosures( hblx, 
          invl, pathslx, paths_indexl, dupl, badl, all_closuresl, False );
     #pragma omp critical
     {    c6a1 += WallClockTime( );
          c6a2 -= WallClockTime( );    }

     // Add s1 and s2 paths into closures.
     
     all_closuresl.append(ss);
     for ( int m = 0; m < ss.isize( ); m++ )
     {    ss[m].ReverseMe( );
          for ( int j = 0; j < ss[m].isize( ); j++ )
               ss[m][j] = invl[ ss[m][j] ];    }
     all_closuresl.append(ss);
     UniqueSort(all_closuresl);

     // Make graph from closures.

     ClosuresToGraph( hblx, invl, all_closuresl, Dl, dinvl, False );
     #pragma omp critical
     {    c6a2 += WallClockTime( );
          c6a3 -= WallClockTime( );    }
     // NOTE EXPENSIVE CONVERSION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     Cleaner( hblx, invl, pathsl, dupl, Dl, dinvl, dpathsl, False );
     if (verbose)
          cout << Date( ) << ": local graph has " << Dl.E( ) << " edges" << endl;
     if (verbose) cout << Date( ) << ": placing reads" << endl;
     #pragma omp critical
     {    c6a3 += WallClockTime( );
          c6b -= WallClockTime( );    }
     // NOTE EXPENSIVE CONVERSION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     PlaceReads( hblx, pathsl, dupl, Dl, dpathsl, False, SINGLE );
     if (verbose) cout << Date( ) << ": building ebcxl" << endl;
     #pragma omp critical
     {    c6b += WallClockTime( );
          c6c -= WallClockTime( );    }
     VecIntVec ebcxl( hbl.E( ) );
     {    vec<vec<int>> esb2( es.size( ) );
          for ( int i = 0; i < esb.isize( ); i++ )
          {    int j;
               for ( j = i + 1; j < esb.isize( ); j++ )
                    if ( esb[j].first != esb[i].first ) break;
               int e = esb[i].first;
               int p = BinPosition( es, e );
               if ( p >= 0 )
               {    int rp = BinPosition( es, inv[e] );
                    for ( int k = i; k < j; k++ )
                    {    esb2[p].push_back( esb[k].second );
                         esb2[rp].push_back( esb[k].second );    }    }
               i = j - 1;    }
          // #pragma omp parallel for
          for ( int e = 0; e < hbl.E( ); e++ )
          {    vec<int> x;
               for ( int i = 0; i < backto[e].isize( ); i++ )
               {    int f = backto[e][i];
                    for ( auto b : esb2[f] ) x.push_back(b);    }
               UniqueSort(x);
               ebcxl[e].resize( x.size( ) );
               for ( int j = 0; j < x.isize( ); j++ )
                    ebcxl[e][j] = x[j];    }    }
     if ( !results_only )
     {    if (verbose) cout << Date( ) << ": inverting" << endl;
          invert( dpathsl, dpaths_indexl, Dl.E( ) );
          kmersl.resize( hbl.E( ) );
          for ( int i = 0; i < hbl.E( ); i++ )
               kmersl[i] = hbl.Kmers(i);
          if ( alignsbl.size( ) > 0 )
          if (verbose) cout << Date( ) << ": building alignsbl" << endl;
          {    alignsbl.resize( hbl.E( ) );
               for ( int e = 0; e < hbl.E( ); e++ )
               {    vec< triple<int,int,int> > x;
                    int start = 0;
                    for ( int i = 0; i < backto[e].isize( ); i++ )
                    {    int f = backto[e][i];
                         for ( auto b : alignsb[es[f]] ) 
                         {    b.second += start;
                              b.third += start;
                              x.push_back(b);        }
                         start += hbls.Kmers(f);    }
                    UniqueSort(x);
                    alignsbl[e].resize( x.size( ) );
                    for ( int j = 0; j < x.isize( ); j++ )
                         alignsbl[e][j] = x[j];    }    }    }

     // Make scaffolds.  

     if (verbose) cout << Date( ) << ": making scaffolds" << endl;
     for ( int pass = 1; pass <= 2; pass++ )
     {    vec<int> to_left, to_right;
          Dl.ToLeft(to_left), Dl.ToRight(to_right);
          vec<int> lens( Dl.E( ), 0 );
          // #pragma omp parallel for
          for ( int e = 0; e < Dl.E( ); e++ )
          {    if ( Dl.O(e)[0] < 0 ) continue;
               for ( int j = 0; j < Dl.O(e).isize( ); j++ )
                    lens[e] += hbl.Kmers( Dl.O(e)[j] );    }
          if ( pass == 2 )
          {    vec<int> dels;
               for ( int d = 0; d < Dl.E( ); d++ )
               {    int v = to_left[d], w = to_right[d];
                    if ( v < 0 ) continue;
                    if ( Dl.To(v).nonempty( ) || !Dl.From(v).solo( ) ) continue;
                    if ( Dl.From(w).nonempty( ) || !Dl.To(w).solo( ) ) continue;
                    if ( v == w ) continue;
                    if ( lens[d] >= 320 ) continue;
                    dels.push_back(d);    }
               Dl.DeleteEdges(dels);
               CleanupCore( Dl, dinvl );
               RemoveUnneededVertices( Dl, dinvl );
               // NOTE EXPENSIVE CONVERSION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               PlaceReads( HyperBasevectorX(hbl), 
                    pathsl, dupl, Dl, dpathsl, False, SINGLE );    }
          
          // NOTE EXPENSIVE CONVERSION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          Scaffold( HyperBasevectorX(hbl), invl, ebcxl, Dl, dinvl, dpathsl, bidl, 
               datasetsl, False, link_report, SINGLE );    }

     // Note: removing inversion artifacts at this point did not have a big effect
     // on the assembly, and it was not unambiguously positive.

     if ( !results_only )
     {    // NOTE EXPENSIVE CONVERSION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          PlaceReads( HyperBasevectorX(hbl), 
               pathsl, dupl, Dl, dpathsl, False, SINGLE );
          if (verbose) cout << Date( ) << ": indexing paths" << endl;
          invert( dpathsl, dpaths_indexl, Dl.E( ) );    }
     #pragma omp critical
     {    c6c += WallClockTime( );
          c7 -= WallClockTime( );    }

     // Now convert back to original coordinates.  So Dl and Dlp are "identical"
     // assemblies except for the edge space in which they are defined.

     if (verbose) cout << Date( ) << ": converting back to orig" << endl;
     Dlp = Dl;
     for ( int e = 0; e < Dl.E( ); e++ )
     {    const vec<int>& x = Dl.O(e);
          if ( IsCell(x) )
          {    cell c;
               c.CellDecode(x);
               digraphE<vec<int>>& G = c.GMutable( );
               for ( int d = 0; d < G.E( ); d++ )
               {    vec<int> y;
                    for ( auto e : G.O(d) )
                    {    for ( int j = 0; j < backto[e].isize( ); j++ )
                              y.push_back( es[ backto[e][j] ] );    }
                    G.OMutable(d) = y;    }
               vec<int> y;
               c.CellEncode(y);
               Dlp.OMutable(e) = y;    }
          if ( x[0] < 0 ) continue;
          vec<int> y;
          for ( int i = 0; i < x.isize( ); i++ )
          {    int f = x[i];
               for ( int j = 0; j < backto[f].isize( ); j++ )
                    y.push_back( es[ backto[f][j] ] );    }
          Dlp.OMutable(e) = y;    }
     #pragma omp critical
     {    c7 += WallClockTime( );    }
     if (verbose)
     {    cout << Date( ) << ": s1 = " << s1 << ", Dlp.E( ) = " << Dlp.E( )
               << ", BuildLocal1b time = " << TimeSince(clock1b) << endl;    }    }

void BuildLocal2(

     // Edges to walk between:

     const int s1, const int s2,

     // Global inputs:

     const digraphE<vec<int>>& D, const vec<vec<vec<vec<int>>>>& dlines,

     // Local outputs:

     HyperBasevector& hbl,
     vec<int>& invl,
     ReadPathVec& pathsl,
     digraphE<vec<int>>& Dl,
     vec<int>& dinvl,
     ReadPathVec& dpathsl,
     VecULongVec& dpaths_indexl,
     vec<Bool>& dupl,
     vec<vec<vec<vec<int>>>>&  dlinesl,
     digraphE<vec<int>>& Dlp,
     vec<int>& dinvlp,

     // Status:
     // (On Dlp, match beginning of global edge s1 to position p1 on
     // local edge d1, stop of global edges s2 to position p2 on local edge d2.)

     Bool& closed, int& d1, int& p1, int& d2, int& p2,

     // Logging etc.:

     Bool verbose, Bool results_only, const Bool SINGLE, const Bool new_test )
{
     // Determine bounding edges.

     double clock2 = WallClockTime( );

     // Locate ends in local assembly.  This is a first version, and likely
     // too stringent.  First we find the base edges e1 and e2, at the ends of
     // the originating edges s1 and s2 and abutting the gap.  Then we find the 
     // instances of these in the local assembly, given by loc1 and loc2.

     if (verbose) cout << Date( ) << ": locating ends" << endl;
     vec<Bool> used;
     Dlp.Used(used);
     int e1 = D.O(s1).front( ), e2 = D.O(s2).back( );
     if (verbose) DPRINT2( e1, e2 );
     vec< pair<int,int> > loc1, loc2;
     for ( int d = 0; d < Dlp.E( ); d++ )
     {    if ( !used[d] || Dlp.O(d)[0] < 0 ) continue;
          for ( int i = 0; i < Dlp.O(d).isize( ); i++ )
          {    int f = Dlp.O(d)[i];
               if ( f == e1 ) loc1.push( d, i );
               if ( f == e2 ) loc2.push( d, i );    }    }

     vec<int> dto_left, dto_right;
     Dlp.ToLeft(dto_left), Dlp.ToRight(dto_right);
     vec< triple< pair<int,int>, pair<int,int>, vec<int> > > bangs;
     double btime = WallClockTime( );
     for ( int z1 = 0; z1 < loc1.isize( ); z1++ )
     for ( int z2 = 0; z2 < loc2.isize( ); z2++ )
     {    int d1 = loc1[z1].first, d2 = loc2[z2].first;
          int p1 = loc1[z1].second, p2 = loc2[z2].second;
          vec<int> be;
          if ( d1 == d2 ) be = {d1};
          else
          {    int v1 = dto_right[d1], v2 = dto_left[d2];
               be = Dlp.EdgesSomewhereBetween( v1, v2 );
               if ( be.empty( ) ) continue;
               be.push_back( d1, d2 );
               UniqueSort(be);    }
          bangs.push( loc1[z1], loc2[z2], be );   }
     if (verbose)
     {    cout << Date( ) << ": s1 = " << s1 << ", s2 = " << s2 
               << ", Dlp.E( ) = " << Dlp.E( )
               << ", |loc1| = " << loc1.size( ) << ", |loc2| = " << loc2.size( ) 
               << ", |bangs| = " << bangs.size( )
               << ", bangtime = " << TimeSince(btime) << endl;    }

     closed = False;

     if ( bangs.empty( ) )
     {    if (verbose) 
          {    cout << Date( ) << ": no closures found" << endl;
               cout << "locations of left end:\n";
               for ( int j = 0; j < loc1.isize( ); j++ )
               {    cout << "[" << j+1 << "] edge " << loc1[j].first
                              << ", position " << loc1[j].second << endl;    }
               cout << "locations of right end:\n";
               for ( int j = 0; j < loc2.isize( ); j++ )
               {    cout << "[" << j+1 << "] edge " << loc2[j].first
                         << ", position " << loc2[j].second << endl;    }    }    }

     if ( bangs.size( ) > 1 )
     {    if (verbose)
          {    cout << Date( ) << ": found " << bangs.size( ) << " closures"
                    << endl;
               for ( int j = 0; j < bangs.isize( ); j++ )
               {    cout << "[" << j+1 << "] "
                         << "edge " << bangs[j].first.first
                         << "(pos " << bangs[j].first.second << ") --> "
                         << "edge " << bangs[j].second.first
                         << "(pos " << bangs[j].second.second << ")" 
                         << endl;    }    }    }

     if ( bangs.solo( ) )
     {    closed = True;
          d1 = bangs[0].first.first, d2 = bangs[0].second.first;
          p1 = bangs[0].first.second, p2 = bangs[0].second.second;
          vec<int> be = bangs[0].third;
          if (verbose) DPRINT( be.size( ) );

          if (new_test)
          {    int v = dto_left[d1], w = dto_right[d2];
               vec<int> dels;
               for ( int i = 0; i < Dlp.To(v).isize( ); i++ )
               {    int d = Dlp.ITo(v,i);
                    if ( d != d1 && d != d2 ) dels.push_back(d,dinvl[d]);    }
               for ( int i = 0; i < Dlp.From(w).isize( ); i++ )
               {    int d = Dlp.IFrom(w,i);
                    if ( d != d1 && d != d2 ) dels.push_back(d,dinvl[d]);    }
               v = dto_right[d1], w = dto_left[d2];
               for ( int i = 0; i < Dlp.To(v).isize( ); i++ )
               {    int d = Dlp.ITo(v,i);
                    if ( d != d1 && d != d2 ) dels.push_back(d,dinvl[d]);    }
               for ( int i = 0; i < Dlp.From(w).isize( ); i++ )
               {    int d = Dlp.IFrom(w,i);
                    if ( d != d1 && d != d2 ) dels.push_back(d,dinvl[d]);    }
               Dlp.DeleteEdges(dels);
               be = Dlp.EdgesSomewhereBetween( dto_right[d1], dto_left[d2] );
               be.push_back( d1, d2 );
               UniqueSort(be);    }

          // Restrict graph to the stuff between.  Make sure ends are completely
          // clean.

          vec<int> dels;
          for ( int d = 0; d < Dlp.E( ); d++ )
          {    if ( !BinMember( be, d ) && !BinMember( be, dinvl[d] ) ) 
                    dels.push_back( d, dinvl[d] );    }
          int v = dto_left[d1], w = dto_right[d2];
          for ( int j = 0; j < Dlp.From(v).isize( ); j++ )
          {    int d = Dlp.IFrom(v,j);
               if ( d != d1 ) dels.push_back( d, dinvl[d] );    }
          for ( int j = 0; j < Dlp.To(w).isize( ); j++ )
          {    int d = Dlp.ITo(w,j);
               if ( d != d2 ) dels.push_back( d, dinvl[d] );    }

          dinvlp = dinvl;
          Dlp.DeleteEdges(dels);
          CleanupCore( Dlp, dinvlp );
          RemoveUnneededVertices( Dlp, dinvlp );

          if ( !results_only )
          {    Dl.DeleteEdges(dels);
               CleanupCore( Dl, dinvl );
               RemoveUnneededVertices( Dl, dinvl );    }

          if (verbose)
          {    cout << Date( ) << ": restricted graph has " << Dlp.E( )
                    << " edges" << endl;    }
          if ( !results_only )
          {    // NOTE EXPENSIVE CONVERSION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               PlaceReads( HyperBasevectorX(hbl), 
                    pathsl, dupl, Dl, dpathsl, False, SINGLE );
               invert( dpathsl, dpaths_indexl, Dlp.E( ) );    }
     
          // Define splicing into global assembly.
     
          vec<Bool> used;
          Dlp.Used(used);
          d1 = -1, d2 = -1; 
          p1 = 0, p2 = 0;
          for ( int d = 0; d < Dlp.E( ); d++ )
          {    if ( Dlp.O(d)[0] < 0 ) continue;
               if ( !used[d] ) continue;
               p1 = Position( Dlp.O(d), e1 );
               if ( p1 >= 0 ) 
               {    d1 = d;
                    break;    }    }
          for ( int d = 0; d < Dlp.E( ); d++ )
          {    if ( Dlp.O(d)[0] < 0 ) continue;
               if ( !used[d] ) continue;
               p2 = Position( Dlp.O(d), e2 );
               if ( p2 >= 0 ) 
               {    d2 = d;
                    break;    }    }

          // Test for fail.

          if ( d1 < 0 || d2 < 0 ) closed = False;
          else
          {
               // So: match beginning of global edge s1 to position p1 on
               // local edge d1, and likewise for 2.

               // Check for weird cases.

               Bool weird = False;
               if ( d1 == d2 && !( p1 <= p2 ) ) weird = True;
               if ( !new_test && !( p1 <= p2 ) ) weird = True;
               if (weird)
               {    if (verbose) 
                    {    cout << Date( ) 
                              << ": positions p1, p2 don't make sense" << endl;    }
                    closed = False;    }    
               if ( ( d1 == d2 && !IsUnique( vec<int>{d1, dinvlp[d1]} ) )
                    || ( d1 != d2 
                         && !IsUnique( vec<int>{d1, d2, dinvlp[d1], dinvlp[d2]} ) ) )
               {    if (verbose) 
                         cout << Date( ) << ": weirdly dependent d values" << endl;
                    closed = False;    }    }    }

     // Get lines.

     if ( !results_only )
     {    const int MIN_LEN = 0;
          FindLines( Dl, dinvl, dlinesl, MAX_CELL_PATHS, MAX_CELL_DEPTH );
          // NOTE EXPENSIVE CONVERSION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          Validate( HyperBasevectorX(hbl), invl, Dl, dinvl );     }
     if (verbose)
     {    cout << Date( ) << ": s1 = " << s1 << ", s2 = " << s2 
               << ", BuildLocal2 time = " << TimeSince(clock2) << endl;    }    }

template void BuildLocal1( const int s1, const vec<int>& s2s, 
     const HyperBasevectorX& hb, const vec<int>& inv,
     const vec<int64_t>& bci, vec<Bool>& dup, vec<Bool>& bad,
     ReadPathVecX& paths,
     VirtualMasterVec< SerfVec<triple<int,int,int> > >& alignsb,
     MasterVec<IntVec>& ebcx,
     const digraphE<vec<int>>& D, const vec<int>& dinv,
     const vec<vec<vec<vec<int>>>>& dlines,
     const ReadPathVec& dpaths, const VecULongVec& dpaths_index,
     const vec<Bool>& internal,
     const vec<int>& to_left, const vec<int>& mult, const vec<int>& lens,
     const Bool use_rights,
     vec<int32_t>& bcl,
     HyperBasevector& hbl,
     HyperBasevectorX& hbxl, vec<int>& kmersl, vec<int>& invl, ReadPathVec& pathsl,
     digraphE<vec<int>>& Dl, vec<int>& dinvl, 
     ReadPathVec& dpathsl, VecULongVec& dpaths_indexl,
     vec<Bool>& dupl,
     MasterVec< SerfVec<triple<int,int,int> > >& alignsbl,
     vec<vec<vec<vec<int>>>>&  dlinesl,
     String& link_report, digraphE<vec<int>>& Dlp,
     Bool verbose, const Bool results_only,
     double& c1, double& c2a, double& c2b, double& c3, double& c4, 
     double& c5, double& c6a1, double& c6a2, double& c6a3, double& c6b, 
     double& c6c, double& c7, const int MAX_READS, const Bool );

template void BuildLocal1( const int s1, const vec<int>& s2s,
     const HyperBasevectorX& hb, const vec<int>& inv,
     const vec<int64_t>& bci, vec<Bool>& dup, vec<Bool>& bad,
     ReadPathVecX& paths,
     MasterVec< SerfVec<triple<int,int,int> > >& alignsb,
     MasterVec<IntVec>& ebcx,
     const digraphE<vec<int>>& D, const vec<int>& dinv,
     const vec<vec<vec<vec<int>>>>& dlines,
     const ReadPathVec& dpaths, const VecULongVec& dpaths_index,
     const vec<Bool>& internal,
     const vec<int>& to_left, const vec<int>& mult, const vec<int>& lens,
     const Bool use_rights,
     vec<int32_t>& bcl,
     HyperBasevector& hbl,
     HyperBasevectorX& hbxl, vec<int>& kmersl, vec<int>& invl, ReadPathVec& pathsl,
     digraphE<vec<int>>& Dl, vec<int>& dinvl, 
     ReadPathVec& dpathsl, VecULongVec& dpaths_indexl,
     vec<Bool>& dupl,
     MasterVec< SerfVec<triple<int,int,int> > >& alignsbl,
     vec<vec<vec<vec<int>>>>&  dlinesl,
     String& link_report, digraphE<vec<int>>& Dlp,
     Bool verbose, const Bool results_only,
     double& c1, double& c2a, double& c2b, double& c3, double& c4, 
     double& c5, double& c6a1, double& c6a2, double& c6a3, double& c6b, 
     double& c6c, double& c7, const int MAX_READS, const Bool );

template void BuildLocal1( const int s1, const vec<int>& s2s,
     const HyperBasevectorX& hb, const vec<int>& inv,
     const vec<int64_t>& bci, vec<Bool>& dup, vec<Bool>& bad,
     ReadPathVecX& paths,
     VirtualMasterVec< SerfVec<triple<int,int,int> > >& alignsb,
     VirtualMasterVec<IntVec>& ebcx,
     const digraphE<vec<int>>& D, const vec<int>& dinv,
     const vec<vec<vec<vec<int>>>>& dlines,
     const ReadPathVec& dpaths, const VecULongVec& dpaths_index,
     const vec<Bool>& internal,
     const vec<int>& to_left, const vec<int>& mult, const vec<int>& lens,
     const Bool use_rights,
     vec<int32_t>& bcl,
     HyperBasevector& hbl,
     HyperBasevectorX& hbxl, vec<int>& kmersl, vec<int>& invl, ReadPathVec& pathsl,
     digraphE<vec<int>>& Dl, vec<int>& dinvl, 
     ReadPathVec& dpathsl, VecULongVec& dpaths_indexl,
     vec<Bool>& dupl,
     MasterVec< SerfVec<triple<int,int,int> > >& alignsbl,
     vec<vec<vec<vec<int>>>>&  dlinesl,
     String& link_report, digraphE<vec<int>>& Dlp,
     Bool verbose, const Bool results_only,
     double& c1, double& c2a, double& c2b, double& c3, double& c4, 
     double& c5, double& c6a1, double& c6a2, double& c6a3, double& c6b, 
     double& c6c, double& c7, const int MAX_READS, const Bool );

void Surgery( const HyperBasevectorX& hb, const vec<int>& inv,
     digraphE<vec<int>>& D, vec<int>& dinv,
     vec<vec<vec<vec<int>>>>& dlines, vec< pair< int, vec<int> > >& s1s2,
     const vec< vec< triple< int, pair<int,int>, pair<int,int> > > >& closures1,
     const vec< vec< digraphE<vec<int>> > >& closures2,
     const vec< vec< vec<int> > >& closures3, const Bool require_simple,
     const Bool require_simpler, ostream& out )
{
     // Carry out surgery.

     vec<int> to_left, to_right;
     D.ToLeft(to_left), D.ToRight(to_right);
     vec<Bool> edited( D.E( ), False );  // has edge been edited?
     vec<Bool> editedv( D.N( ), False ); // has vertex been edited?
     int nclosed = 0;
     for ( int li = 0; li < s1s2.isize( ); li++ )
     {    int s1 = s1s2[li].first;
          out << "\n";
          PRINT3_TO( out, s1, closures1[li].size( ), int(edited[s1]) );
          if ( !closures1[li].solo( ) ) continue;
          if ( edited[s1] ) continue;
          
          // (On Dlp, match beginning of global edge s1 to position p1 on
          // local edge d1, stop of global edges s2 to position p2 on local edge d2.)
          //
          // If there's just one closure, patch it in. = surgery

          int s2 = closures1[li][0].first;
          PRINT3_TO( out, s1, s2, int(edited[s2]) );
          if ( edited[s2] ) continue;
          int nev = editedv.size( );
          if ( to_left[s1] >= nev || to_left[s2] >= nev ) continue;
          if ( to_right[s1] >= nev || to_right[s2] >= nev ) continue;
          if ( editedv[ to_left[s1] ] || editedv[ to_left[s2] ] ) continue;
          if ( editedv[ to_right[s1] ] | editedv[ to_right[s2] ] ) continue;
          int d1 = closures1[li][0].second.first, d2 = closures1[li][0].third.first;
          int p1 = closures1[li][0].second.second; 
          int p2 = closures1[li][0].third.second;
          PRINT4_TO( out, s1, s2, d1, d2 );
          const digraphE< vec<int> >& Dl = closures2[li][0];
          vec<int> dinvl = closures3[li][0];

          // More tests.  Failures here are likely due to defects in BuildLocal1.

          {    vec<int> xleft, xright;
               Dl.ToLeft(xleft), Dl.ToRight(xright);
               out << "Dl.O(d1) = " << printSeq( Dl.O(d1) ) << "\n";
               out << "Dl.O(d2) = " << printSeq( Dl.O(d2) ) << "\n";
               PRINT_TO( out, Dl.To( xleft[d1] ).size( ) );
               PRINT_TO( out, Dl.From( xleft[d1] ).size( ) );
               PRINT_TO( out, Dl.From( xright[d2] ).size( ) );
               PRINT_TO( out, Dl.To( xright[d2] ).size( ) );
               if ( Dl.To( xleft[d1] ).size( ) != 0 ) continue;
               if ( Dl.From( xleft[d1] ).size( ) != 1 ) continue;
               if ( Dl.From( xright[d2] ).size( ) != 0 ) continue;
               if ( Dl.To( xright[d2] ).size( ) != 1 ) continue;    }

          // Test for simple.  Shouldn't be able to get from d1 to rd1.

          if (require_simpler)
          {    vec<vec<vec<vec<int>>>> dlinesl;
               FindLines( Dl, dinvl, dlinesl, MAX_CELL_PATHS, MAX_CELL_DEPTH );
               PRINT3_TO( out, s1, s2, dlinesl.size( ) );
               if ( dlinesl.size( ) != 2 ) continue;
               if ( dinvl[ dlinesl[0].front( )[0][0] ]
                    != dlinesl[1].back( )[0][0] )
               {    PRINT3_TO( out, s1, s2, "back fail" );
                    continue;    }    }
          else if (require_simple)
          {    int rd1 = dinvl[d1];
               vec<int> dto_left, dto_right;
               Dl.ToLeft(dto_left), Dl.ToRight(dto_right);
               vec<int> x;
               Dl.GetSuccessors1( dto_right[d1], x );
               Sort(x);
               if ( BinMember( x, dto_left[rd1] ) ) 
               {    // cout << "nonsimple" << endl;
                    // continue;    
                         }
               vec<vec<int>> X;
               Dl.ComponentsEFast(X);
               Bool same = False;
               for ( int i = 0; i < X.isize( ); i++ )
               {    Sort( X[i] );
                    if ( BinMember( X[i], d1 ) && BinMember( X[i], rd1 ) )
                    {    same = True;
                         break;    }    }
               if (same)
               {    out << "really not simple" << endl;
                    continue;    }    }
     
          // Proceed with making the patch.

          nclosed++;
          PRINT3_TO( out, s1, s2, "passes" );
          int rs1 = dinv[s1], rs2 = dinv[s2];
          edited[s1] = edited[s2] = edited[rs1] = edited[rs2] = True;
          editedv[ to_left[s1] ] = editedv[ to_right[s1] ] = True;
          editedv[ to_left[s2] ] = editedv[ to_right[s2] ] = True;
          editedv[ to_left[rs1] ] = editedv[ to_right[rs1] ] = True;
          editedv[ to_left[rs2] ] = editedv[ to_right[rs2] ] = True;
          for ( int e = 0; e < dinvl.isize( ); e++ )
               dinvl[e] += D.E( );
          d1 += D.E( ), d2 += D.E( );
          double clock1 = WallClockTime( );
          D.AppendWithUpdate( Dl, to_left, to_right );
          dinv.append(dinvl);
          // Validate( hb, inv, D, dinv ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          out << Date( ) << ": s1 = " << s1 << ", " << "Dl.E( ) = " << Dl.E( )
               << ", " << TimeSince(clock1) << " used appending" << endl;
          int rd1 = dinv[d1], rd2 = dinv[d2];

          // If there's a gap edge in between, delete it.

          int v = to_right[s1], w = to_left[s2];
          if ( D.From(v).solo( ) && D.From(v)[0] == w )
          {    int g = D.IFrom(v,0);
               /*
               ForceAssert( IsBarcodeOnlyGap( D.O(g) ) ); // XXXXXXXXXXXXXXXXXXXXXXX
               #pragma omp critical
               {    cout << "killing gap" << endl;    }
               */
               D.DeleteEdgesWithUpdate( { g, dinv[g] }, to_left, to_right );    }

          D.TransferEdgesWithUpdate( to_left[d1], to_left[s1], to_left, to_right );
          D.TransferEdgesWithUpdate( 
               to_right[d2], to_right[s2], to_left, to_right );
          D.TransferEdgesWithUpdate( to_left[rd2], to_left[rs2], to_left, to_right );
          D.TransferEdgesWithUpdate( 
               to_right[rd1], to_right[rs1], to_left, to_right );
          if ( d1 == d2 )
          {    vec<int> x;
               for ( int j = p1; j <= p2; j++ )
                    x.push_back( D.O(d1)[j] );
               D.OMutable(d1) = x;
               x.ReverseMe( );
               for ( int j = 0; j < x.isize( ); j++ )
                    x[j] = inv[ x[j] ];
               D.OMutable(rd1) = x;    }
          else
          {    vec<int> x1, x2;
               for ( int j = p1; j < D.O(d1).isize( ); j++ )
                    x1.push_back( D.O(d1)[j] );
               ForceAssertEq( x1.front( ), D.O(s1).front( ) );
               for ( int j = 0; j <= p2; j++ )
                    x2.push_back( D.O(d2)[j] );
               ForceAssertEq( x2.back( ), D.O(s2).back( ) );
               D.OMutable(d1) = x1, D.OMutable(d2) = x2;
               x1.ReverseMe( ), x2.ReverseMe( );
               for ( int j = 0; j < x1.isize( ); j++ )
                    x1[j] = inv[ x1[j] ];
               for ( int j = 0; j < x2.isize( ); j++ )
                    x2[j] = inv[ x2[j] ];
               D.OMutable(rd1) = x1, D.OMutable(rd2) = x2;    }
          // Validate( hb, inv, D, dinv );    
               }
     cout << Date( ) << ": total uncaptured gaps closed = " << nclosed << endl;
     Validate( hb, inv, D, dinv );    }

template<class VP2, class VI, class VA, class VE> void Unvoid( 
     const vec< pair< int, vec<int> > >& s1s2, const HyperBasevectorX& hb, 
     const vec<int>& inv, const vec<int64_t>& bci, 
     vec<Bool>& dup, vec<Bool>& bad, ReadPathVecX& paths, VE& ebcx, VA& alignsb,
     digraphE<vec<int>>& D, vec<int>& dinv, vec<vec<vec<vec<int>>>>& dlines,
     VP2& dpaths, VI& dpaths_index, Bool use_rights, ostream& gout,
     vec< vec< triple< int, pair<int,int>, pair<int,int> > > >& closures1,
     vec< vec< digraphE<vec<int>> > >& closures2, vec< vec< vec<int> > >& closures3,
     const String& DIR, const Bool SAVE_LOCAL, const int MAX_READS, 
     const Bool SINGLE, const Bool verbose, const Bool new_test )
{    
     // Set up variables.

     closures1.clear( ), closures2.clear( ), closures3.clear( );
     closures1.resize( s1s2.size( ) );
     closures2.resize( s1s2.size( ) );
     closures3.resize( s1s2.size( ) );
     double clockb1 = 0, clockb2 = 0, clockb3 = 0;
     double c1 = 0, c2a = 0, c2b = 0, c3 = 0, c4 = 0;
     double c5 = 0, c6a1 = 0, c6a2 = 0, c6a3 = 0, c6b = 0, c6c = 0, c7 = 0;
     double bclock = WallClockTime( );
     double mclock = WallClockTime( );
     String link_report;
     if (SAVE_LOCAL) Mkdir777( DIR + "/a.local" );

     // Set up ancillary data structures.

     vec<Bool> internal( paths.size( ), False );
     vec<int> to_left, to_right, mult, lens;
     #pragma omp parallel sections
     {
          #pragma omp section
          {    D.ToLeft(to_left), D.ToRight(to_right);    }
          #pragma omp section
          {    ComputeMult( hb, D, mult );    }
          #pragma omp section
          {    lens.resize( D.E( ), 0 );
               #pragma omp parallel for
               for ( int e = 0; e < D.E( ); e++ )
               {    if ( D.O(e)[0] < 0 ) continue;
                    for ( int j = 0; j < D.O(e).isize( ); j++ )
                         lens[e] += hb.Kmers( D.O(e)[j] );    }    }    }
     vec<int> tol( D.E( ), -1 );
     for ( int i = 0; i < dlines.isize( ); i++ )
     for ( int j = 0; j < dlines[i].isize( ); j++ )
     for ( int k = 0; k < dlines[i][j].isize( ); k++ )
     for ( int l = 0; l < dlines[i][j][k].isize( ); l++ )
          tol[ dlines[i][j][k][l] ] = i;

     // Start main loop.

     cout << Date( ) << ": memory in use = " << MemUsageGBString( ) << endl;
     cout << Date( ) << ": looping over " << s1s2.size( ) 
          << " elements, showing 100 dots:" << endl;
     int ndots = 0, stopped = 0;
     #pragma omp parallel for schedule( dynamic, 1 )
     for ( int i = 0; i < s1s2.isize( ); i++ )
     {    const int s1 = s1s2[i].first;
          const vec<int>& s2s = s1s2[i].second;
          #pragma omp critical
          {    clockb1 -= WallClockTime( );    }
          ostringstream mout;
          double pclock = WallClockTime( );

          // Set up variables.

          vec<int32_t> bcl;
          HyperBasevector hbl;
          HyperBasevectorX hbxl;
          vec<int> kmersl;
          vec<int> invl;
          ReadPathVec pathsl;
          digraphE<vec<int>> Dl; 
          vec<int> dinvl; 
          ReadPathVec dpathsl;
          VecULongVec dpaths_indexl;
          vec<Bool> dupl;
          MasterVec< SerfVec<triple<int,int,int> > > alignsbl;
          vec<vec<vec<vec<int>>>> dlinesl;
          String link_report;
          Bool closed;
          int d1, p1, d2, p2;
          digraphE<vec<int>> Dlp0;
          vec<int> dinvlp0;
          #pragma omp critical
          {    clockb1 += WallClockTime( );
               clockb2 -= WallClockTime( );    }

          // Build local assembly, part 1.

          Bool results_only = True;
          if (SAVE_LOCAL) results_only = False;
          BuildLocal1( s1, s2s, hb, inv, bci, dup, bad, 
               paths, alignsb, ebcx, D, dinv, dlines, dpaths, dpaths_index, 
               internal, to_left, mult, lens, use_rights, bcl, hbl, hbxl, 
               kmersl, invl, pathsl, Dl, dinvl, dpathsl, dpaths_indexl, dupl, 
               alignsbl, dlinesl, link_report, Dlp0, verbose, results_only, c1, 
               c2a, c2b, c3, c4, c5, c6a1, c6a2, c6a3, c6b, c6c, c7, MAX_READS,
               SINGLE );

          // Got through each possibility.

          #pragma omp critical
          {    clockb2 += WallClockTime( );
               clockb3 -= WallClockTime( );    }

          mout << "\n" << Date( ) << ": using S1 = " << s1 << endl;
          for ( int si = 0; si < s2s.isize( ); si++ )
          {    const int s2 = s2s[si];
               mout << Date( ) << ": ==> S2 = " << s2 << endl;
          
               // Build local assembly, part 2.
     
               digraphE<vec<int>> Dlp(Dlp0);
               vec<int> dinvlp(dinvlp0); // DOES THIS MAKE SENSE?
               BuildLocal2( s1, s2, D, dlines, hbl, invl, pathsl, Dl, 
                    dinvl, dpathsl, dpaths_indexl, dupl, dlinesl, Dlp, dinvlp, 
                    closed, d1, p1, d2, p2, verbose, results_only, SINGLE,
                    new_test );

               if (closed) 
               {    mout << Date( ) << ": closed!" << endl;
                    closures1[i].push( s2, make_pair(d1,p1), make_pair(d2,p2) );
                    closures2[i].push_back( Dlp );
                    closures3[i].push_back( dinvlp );    }
               if (SAVE_LOCAL)
               {    String DIRL = DIR + "/a.local/" 
                         + ToString(s1) + "." + ToString( D.O(s1).back( ) );
                    Mkdir777(DIRL);
                    Ofstream( out, DIRL + "/link_report" );
                    out << link_report;
                    Ofstream( bout, DIRL + "/bc" );
                    for ( int j = 0; j < bcl.isize( ); j++ )
                         bout << bcl[j] << "\n";
                    BinaryWriter::writeFile( DIRL + "/a.sup", Dl );
                    BinaryWriter::writeFile( DIRL + "/a.kmers", kmersl );
                    BinaryWriter::writeFile( DIRL + "/a.inv", invl );
                    pathsl.WriteAll( DIRL + "/a.paths" );
                    VecULongVec pathsl_index;
                    invert( pathsl, pathsl_index );
                    pathsl_index.WriteAll( DIRL + "/a.paths.inv" );
                    BinaryWriter::writeFile( DIRL + "/a.sup.inv", dinvl );
                    dpathsl.WriteAll( DIRL + "/a.dpaths" );
                    dpaths_indexl.WriteAll( DIRL + "/a.dpaths.index" );
                    alignsbl.WriteAll( DIRL + "/a.alignsb" );
                    BinaryWriter::writeFile( DIRL + "/a.dup", dupl );
                    BinaryWriter::writeFile( DIRL + "/a.hbx", hbxl );
                    BinaryWriter::writeFile( DIRL + "/a.sup.lines", dlinesl );
                    vecbasevector edges( hbxl.E( ) );
                    for ( int e = 0; e < (int) hbxl.E( ); e++ )
                         edges[e] = hbxl.O(e);
                    edges.WriteAll( DIRL + "/a.fastb" );
                    Cp2( DIRL + "/../../../genome.names", DIRL );    }    }
          
          // If there are two closures, see if one is clearly the winner.

          ChooseClosure( 
               hb, D, dlines, tol, closures1[i], closures2[i], closures3[i], mout );

          // Dump output and note progress.

          mout << Date( ) << ": done, time used = " << TimeSince(pclock) << endl;
          #pragma omp critical
          {    gout << mout.str( );
               clockb3 += WallClockTime( );
               MakeDots( stopped, ndots, s1s2.size( ) );    }    }

     // Summarize.

     DPRINT3( clockb1, clockb2, clockb3 );
     DPRINT3( c1, c2a, c2b );
     DPRINT3( c3, c4, c5 );
     DPRINT3( c6a1, c6a2, c6a3 );
     DPRINT3( c6b, c6c, c7 );
     cout << Date( ) << ": " << TimeSince(bclock) << " used capturing gaps" 
          << endl;    }

template void Unvoid( const vec< pair< int, vec<int> > >& s1s2, 
     const HyperBasevectorX& hb, const vec<int>& inv, 
     const vec<int64_t>& bci, vec<Bool>& dup, vec<Bool>& bad, 
     ReadPathVecX & paths, MasterVec<IntVec>& ebcx, 
     MasterVec< SerfVec<triple<int,int,int> > >& alignsb,
     digraphE<vec<int>>& D, vec<int>& dinv, vec<vec<vec<vec<int>>>>& dlines,
     MasterVec<ReadPath>& dpaths, 
     MasterVec<ULongVec>& dpaths_index, Bool use_rights, ostream& gout,
     vec< vec< triple< int, pair<int,int>, pair<int,int> > > >& closures1,
     vec< vec< digraphE<vec<int>> > >& closures2, vec< vec< vec<int> > >& closures3,
     const String& DIR, const Bool SAVE_LOCAL, const int MAX_READS, const Bool,
     const Bool, const Bool );

template void Unvoid( const vec< pair< int, vec<int> > >& s1s2, 
     const HyperBasevectorX& hb, const vec<int>& inv, 
     const vec<int64_t>& bci, vec<Bool>& dup, vec<Bool>& bad, 
     ReadPathVecX& paths, VirtualMasterVec<IntVec>& ebcx, 
     VirtualMasterVec< SerfVec<triple<int,int,int> > >& alignsb,
     digraphE<vec<int>>& D, vec<int>& dinv, vec<vec<vec<vec<int>>>>& dlines,
     MasterVec<ReadPath>& dpaths, 
     MasterVec<ULongVec>& dpaths_index, Bool use_rights, ostream& gout,
     vec< vec< triple< int, pair<int,int>, pair<int,int> > > >& closures1,
     vec< vec< digraphE<vec<int>> > >& closures2, vec< vec< vec<int> > >& closures3,
     const String& DIR, const Bool SAVE_LOCAL, const int MAX_READS, const Bool,
     const Bool, const Bool );

// If there are two closures, see if one is clearly the winner.

void ChooseClosure( const HyperBasevectorX& hb, const digraphE<vec<int>>& D,
     const vec<vec<vec<vec<int>>>>& dlines, const vec<int>& tol,
     vec< triple< int, pair<int,int>, pair<int,int> > >& closures1,
     vec< digraphE<vec<int>> >& closures2, vec< vec<int> >& closures3,
     ostream& mout )
{
     if ( closures1.size( ) == 2 )
     {    vec<vec<int>> dcontent(2), lcontent(2);
          vec<int> ltotal( 2, 0 );
          vec<int> l2(2);
          for ( int j = 0; j < 2; j++ )
          {    const digraphE<vec<int>>& Dl = closures2[j];
               for ( int l = 0; l < Dl.E( ); l++ )
                    if ( Dl.O(l)[0] >= 0 ) dcontent[j].append( Dl.O(l) );
               l2[j] = tol[ closures1[j].first ];
               const vec<vec<vec<int>>>& L2 = dlines[l2[j]];
               for ( int a = 0; a < L2.isize( ); a++ )
               for ( int b = 0; b < L2[a].isize( ); b++ )
               for ( int c = 0; c < L2[a][b].isize( ); c++ )
               {    int d = L2[a][b][c];
                    if ( D.O(d)[0] >= 0 ) lcontent[j].append( D.O(d) );    }
               UniqueSort( dcontent[j] ), UniqueSort( lcontent[j] );
               for ( int i = 0; i < lcontent[j].isize( ); i++ )
                    ltotal[j] += hb.Kmers( lcontent[j][i] );    }
          vec<double> frac(2);
          for ( int j = 0; j < 2; j++ )
          {    int present = 0;
               for ( int i = 0; i < lcontent[j].isize( ); i++ )
               {    int e = lcontent[j][i];
                    if ( BinMember( dcontent[1-j], e ) ) present += hb.Kmers(e);    }
               mout << "frac of kmers from L2 = " << l2[j] 
                    << " that are contained in assembly to L2 = " << l2[1-j] << " = "
                    << PERCENT_RATIO( 3, present, ltotal[j] ) << endl; 
               frac[j] = double(present) / double( ltotal[j] );    }
          double MIN_AD = 0.9;
          if ( frac[0] - frac[1] >= MIN_AD )
          {    mout << "L" << l2[0] << " wins" << endl;
               closures1.resize(1);
               closures2.resize(1);
               closures3.resize(1);    }
          if ( frac[1] - frac[0] >= MIN_AD )
          {    mout << "L" << l2[1] << " wins" << endl;
               closures1 = { closures1[1] };
               closures2 = { closures2[1] };
               closures3 = { closures3[1] };    }    }    }
