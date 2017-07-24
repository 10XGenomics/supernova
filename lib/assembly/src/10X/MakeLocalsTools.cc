// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"
#include "Equiv.h"
#include "ParallelVecUtilities.h"
#include "TokenizeString.h"
#include "feudal/PQVec.h"
#include "graph/DigraphTemplate.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "paths/long/large/GapToyTools.h"
#include "paths/long/large/Lines.h"
#include "10X/Capture.h"
#include "10X/CleanThe.h"
#include "10X/mergers/ClosuresToGraph.h"
#include "10X/Gap.h"
#include "10X/Heuristics.h"
#include "10X/LineOO.h"
#include "10X/LocalTools.h"
#include "10X/MakeLocalsTools.h"
#include "10X/Peeler.h"
#include "10X/PlaceReads.h"
#include "10X/PullApart.h"
#include "10X/Scaffold.h"
#include "10X/SecretOps.h"
#include "10X/Super.h"
#include "10X/astats/FinAlign.h"
#include "10X/astats/RefLookup.h"
#include "10X/astats/View.h"
#include "10X/mergers/GetMergers.h"
#include "10X/DfTools.h"

void FindLinesFixed( const digraphE<vec<int>>& D, const vec<int>& dinv,
     vec<vec<vec<vec<int>>>>& dlines, const int max_cell_paths,
     const int max_cell_depth, const Bool verbose )
{    FindLines( D, dinv, dlines, max_cell_paths, max_cell_depth, verbose );
     if (verbose) cout << Date( ) << ": computing ancillaries" << endl;
     vec<int> to_left, to_right;
     vec<pair<int,int>> tol;
     MakeTol2( D, dlines, tol );
     D.ToLeft(to_left), D.ToRight(to_right);
     for ( int d = 0; d < D.E( ); d++ )
     {    if ( tol[d].first >= 0 ) continue;
          int v = to_left[d], w = to_right[d];
          if ( v < 0 || w < 0 ) continue;
          if ( D.To(v).solo( ) && D.From(w).solo( ) )
          {    int f1 = D.ITo(v,0), f2 = D.IFrom(w,0);
               if ( tol[f1].first >= 0 && tol[f2].first == tol[f1].first
                    && tol[f2].second == tol[f1].second + 2 )
               {    continue;    }    }
          tol[d] = make_pair( dlines.size( ), 0 );
          dlines.push_back( {{{d}}} );
          if ( tol[ dinv[d] ].first < 0 )
          {    tol[ dinv[d] ] = make_pair( dlines.size( ), 0 );
               dlines.push_back( {{{dinv[d]}}} );    }    }    }

void LBC2GBC(const vec<vec<int>>& BCLIST, const int A, const vec<int>& bcseq){
    cout<<"Barcodes: ";
    for(auto lbc:bcseq)
        cout<<BCLIST[A][lbc]<<",";
    cout<<"-"<<endl;
}

int GE2LA(const vec<digraphE<vec<int>>>& DL, const int64_t gE){
   // compute dci
   vec<int64_t> dci;
   int64_t N = 0;
   for(int j = 0; j< DL.isize() ; j++){
       N += DL[j].E();
       dci.push_back(N);
   }
   int lA = 0;
   while(gE<dci.back()){
        if(gE>dci[lA])
            lA++;
        else
            break;
   }
   lA--;
   cout<<"lA:"<<lA<<", lE:"<<gE-dci[lA]<<endl;
   return lA; 
}

int GE2LA(const vec<int64_t>& dci, const int64_t gE){
   // compute dci
   int lA = 0;
   while(gE<dci.back()){
        if(gE>dci[lA])
            lA++;
        else
            break;
   }
   lA--;
   cout<<"lA:"<<lA<<", lE:"<<gE-dci[lA]<<endl;
   return lA; 
}

// template<class T> 
void MakeLocalAssemblies( const HyperBasevectorX& hb,
     const vec<int>& inv, const vec<int64_t>& bci, 
     // T& pathsx, 
     ReadPathVecX& pathsx,
     const vec<Bool>& dup,
     VirtualMasterVec<basevector>& bases, VirtualMasterVec<PQVec>& quals,
     const vec<vec<int>>& bs, vec< digraphE<vec<int>> >& DL,
     vec< vec<int> >& DLI, vec< vec<int> >& BCLIST,
     vec<vec<vec<int>>>& BPATHS, const int MAX_BARCODES, const int MIN_LINE_TO_WALK,
     const String& udir, const String& suffix, 
     const String& DIR, const MasterVec< SerfVec<triple<int,int,int> > >& alignsb, 
     Bool C2G, Bool STORE_LOCALS, const Bool PRINT_DETAILS, const String& INSTANCE,
     const Bool NEW_STUFF, const Bool CANON )
{

     // Delete preexisting local assemblies.

     if (STORE_LOCALS)
     {    String lldir = udir + "/" + INSTANCE + "/locas";
          if ( IsDirectory(lldir) )
          {    vec<String> all = AllFiles(lldir);
               for ( auto d : all )
               {    if ( d.IsInt( ) )
                         System( "rm -rf " + lldir + "/" + d );    }    }    }

     // some lambdas
     auto RecordMe = [&](double& clk){
         clk = WallClockTime(); };
     auto AccrueTime = [&](vec<double>& clks, vec<double>& clks2, int finalstage){
         #pragma omp critical
         {
             for(int s = 0; s < finalstage; s++)
                 clks[s] += clks2[s+1]-clks2[s];
         }
     };

     // Find all edges that appear in each barcode set.

     cout << Date( ) << ": finding all edges in barcode sets" << endl;
     double clock = WallClockTime( );
     DL.resize( bs.size( ) ), DLI.resize( bs.size( ) );
     vec<int> minv( hb.E( ) );
     #pragma omp parallel for
     for ( int e = 0; e < hb.E( ); e++ )
          minv[e] = Min( e, inv[e] );
     int nthreads = omp_get_max_threads( );
     vec<Bool> processed( bs.size( ), False );
     double iclock = WallClockTime( );

     vec<double> clks(30,0); // one for each stage as marked
     int totalstages = 0;
     // #pragma omp parallel for schedule(dynamic,1) firstprivate(paths)
     #pragma omp parallel for schedule(dynamic,1)
     for ( int bi = 0; bi < nthreads; bi++ )
     {    
          // Shared data structures.

          vec<int> es;
          vec<triple<int,int32_t,int32_t>> esb;
          vec<int> invl;
          vec< pair<int64_t,int> > idsb;
          HyperBasevector hbl;
          vec<int32_t> bcl;
          vec<int> dinvl;

          // Find barcode group to process.

          int z = 0;
          while(1)
          {  
               // Advance to next unprocessed barcode group.
               if ( z == bs.isize( ) ) break;
               if ( processed[z] )
               {    z++;
                    continue;    }
               else
               {
                    Bool next = False;
                    #pragma omp critical
                    {    if ( processed[z] )
                         {    z++;
                              next = True;    }
                         else 
                         {    processed[z] = True;
                              if ( z > 0 && z % 100 == 0 ) 
                              {    double t = WallClockTime( ) - iclock;
                                   t *= bs.size( ) / double(z);
                                   t /= 3600;
                                   cout << z << " of " << bs.size( )
                                        << ", predicted total time = "
                                        << t << " hours" << endl;    }    }    }
                    if (next) continue;    }

               // Now assemble it.  Note that the two references to pathsx
               // below are NOT appropriate the global case.

               vec<double> iclks(30,0); // one for each stage as marked
               int stage = 0;

               // ------------  STAGE 1 ------------------------------
               // Take set of barcodes in bs, and extract all reads
               // Store each base edge instance in reads 
               // Sort on base edge to line up all of its instances
               RecordMe(iclks[stage++]);
               es.clear( ), esb.clear( );
               const vec<int>& b = bs[z];
               for ( auto x : b )
               {    for ( int64_t id = bci[x]; id < bci[x+1]; id++ )
                    {    if ( dup[id/2] ) continue;
                         // const ReadPath& p = paths[id];
                         ReadPath p;
                         pathsx.unzip( p, hb, id );
                         for ( auto e : p )
                              esb.push( minv[e], x, id - bci[x] );    }    }
               Sort(esb);    

               
               // ------------  STAGE 2 ------------------------------
               // Scan through sorted base edge instances
               // for e, find number of bc associated with it and total instances
               // store e and inv(e) in es if number of bc is too low, 
               // if number of instances is too low
               // or if it has large #kmers
               // rewrite esb accordingly
               // WARN: could cause uncaps!!
               RecordMe(iclks[stage++]);
               int c = 0;
               for ( int i = 0; i < esb.isize( ); i++ )
               {    int j, nbc = 1, e = esb[i].first;
                    for ( j = i + 1; j < esb.isize( ); j++ )
                    {    if ( esb[j].first != e ) break;
                         if ( esb[j].second != esb[j-1].second ) nbc++;    }
                    int n = j - i;
                    if ( n == 1 || ( n < 4 && hb.Kmers(e) >= 1000 )
                         || ( nbc < 3 && hb.Kmers(e) >= 500 )
                         || ( nbc < 2 && hb.Kmers(e) >= 200 ) );
                    else 
                    {    for ( int k = i; k < j; k++ ) esb[c++] = esb[k];
                         es.push_back(e);
                         if ( inv[e] != e ) es.push_back( inv[e] );    }
                    i = j - 1;    }
               esb.resize(c);
               Sort(es);

               // Build superassembly, following BuildLocal1.

               hbl.Clear( );
               ReadPathVec pathsl;
               bcl.clear( );
               digraphE<vec<int>> Dl;
               dinvl.clear( );
               ReadPathVec dpathsl;

               // ------------  STAGE 3 ------------------------------
               // Form local HyperBasevector.
               // This is very slow.  The bulk of the slowness may be because
               // resizing vec<vec<int>> objects is expensive, but that's just
               // a hypothesis.  The binary search in the constructor adds to
               // that, but is a relatively small component.

               RecordMe(iclks[stage++]);
               // Creating a subgraph consisting of only edges we picked
               ( (digraphE<basevector>&) hbl ).Initialize( 
                    digraphE<basevector>::COMPLETE_SUBGRAPH_EDGES, hb, es );
               hbl.SetK( hb.K( ) );

               // ------------  STAGE 4 ------------------------------
               // Form local inversion.

               RecordMe(iclks[stage++]);
               invl.clear( );
               invl.resize( es.size( ) );
               for ( int i = 0; i < es.isize( ); i++ )
                    invl[i] = BinPosition( es, inv[ es[i] ] );

               // Make list of finally used global read ids, and their barcodes.  

               idsb.clear( );
               for ( int i = 0; i < esb.isize( ); i++ )
               {    int b = esb[i].second;
                    int64_t id = bci[b] + esb[i].third;
                    idsb.push( id, b );    }
               UniqueSort(idsb);
     
               // Form local paths of selected reads, truncating as needed.  
               // Note that offsets should be adjusted, but they're not.
               // Populate barcode used

               vec<int> x, y, lens;
               vec<int64_t> ids2;
               for ( int z = 0; z < idsb.isize( ); z++ )
               {    int64_t id = idsb[z].first;
                    int b = idsb[z].second;

                    // Check for missing partner.

                    if ( id % 2 == 1 && ( z == 0 || idsb[z-1].first != id - 1 ) )
                    {    pathsl.push_back( ReadPath( ) );
                         ids2.push_back(id-1);
                         bcl.push_back(b);    }
                    ids2.push_back(id);
                              
                    // Set x to the longest subpath in es.

                    // const ReadPath& p = paths[id];
                    ReadPath p;
                    pathsx.unzip( p, hb, id );
                    x.clear( );
                    for ( auto e : p ) x.push_back(e);
                    vec<vec<int>> xs;
                    for ( int i = 0; i < x.isize( ); i++ )
                    {    if ( !BinMember( es, x[i] ) ) continue;
                         int j;
                         for ( j = i + 1; j < x.isize( ); j++ )
                              if ( !BinMember( es, x[j] ) ) break;
                         y.clear( );
                         for ( int k = i; k < j; k++ ) y.push_back( x[k] );
                         xs.push_back(y);
                         i = j - 1;    }
                    {    lens.clear( );
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

                    // Save path.

                    ReadPath pl;
                    pl.setOffset( p.getOffset( ) );
                    SerfVec<int> y( x.size( ) );
                    for ( int j = 0; j < (int) x.isize( ); j++ )
                         y[j] = BinPosition( es, x[j] );
                    (SerfVec<int>&) pl = y;
                    pathsl.push_back(pl);
                    bcl.push_back(b);

                    // Check for missing partner.

                    if ( id % 2 == 0 && ( z == idsb.isize( ) - 1 
                         || idsb[z+1].first != id + 1 ) )
                    {    pathsl.push_back( ReadPath( ) );
                         ids2.push_back(id+1);
                         bcl.push_back(b);    }    }

               // ------------  STAGE 5 ------------------------------
               // Simplify the graph.  Start by removing unneeded vertices.
               // Note that we observe cases of circles where no unneeded 
               // vertices are removed, and there is no intersection with the
               // involution.  Just an interesting observation, nothing to do 
               // with this module per se.
               //
               // We keep the original HyperBasevector, so we can translate 
               // back to it.
               
               RecordMe(iclks[stage++]);
               HyperBasevector hbls(hbl);
               RemoveUnneededVertices2( hbl, invl, pathsl, False, True );
               CleanupCore( hbl, invl, pathsl );

               // NEW CODE BELOW.
               //
               // Identify certain "bad" local base edges.  First we declare local
               // base edges that are palindromes to be bad.  Second, in a commented
               // out block of code, we identify local base edges tha are inversion
               // crosses, much as in SnipFlipSquares.  The reason that this block
               // of code is commented out is that it may not be needed.
               //
               // We don't actually delete the bad edges because given the code
               // organization, it's not clear how to do that.  Instead, below,
               // we delete pair closures that involve any of the bad edges.
               //
               // It's not clear that palindromes should be deleted here, or
               // instead deleted from the global base graph.  Also note that
               // the originating case for this was a single kmer edge, findable
               // at R=10:63.162-63.165.

               /*
               vec<int> delsx;
               for ( int e = 0; e < hbl.E( ); e++ )
                    if ( invl[e] == e ) delsx.push_back(e);
               */

               /*
               {
               const int MAX_CAN_INS_DEL = 3;
               const int MIN_CAN_INS_RATIO = 5;
               VecULongVec paths_indexl;
               invert( pathsl, paths_indexl, hbl.E( ) );
               for ( int v = 0; v < hbl.N( ); v++ )
               {    if ( hbl.From(v).size( ) != 2 ) continue;
                    for ( int j = 0; j < 2; j++ )
                    {    int e1 = hbl.IFrom(v,j), e2 = hbl.IFrom(v,1-j);
                         int w = hbl.From(v)[j];
                         if ( !hbl.From(w).solo( ) ) continue;
                         if ( !hbl.To(v).solo( ) ) continue;
                         int x = hbl.From(w)[0], y = hbl.From(v)[1-j];
                         int h = hbl.ITo(v,0), f = hbl.IFrom(w,0);
                         if ( hbl.To(w).size( ) != 2 ) continue;
                         if ( hbl.From(y).size( ) != 2 ) continue;
                         int g = -1;
                         for ( int l = 0; l < 2; l++ )
                         {    if ( hbl.From(y)[l] == x )
                              {    g = hbl.IFrom(y,l);
                                   break;    }    }
                         if ( !IsUnique( vec<int>{ v, w, x, y } ) ) continue;
                         if ( !IsUnique( vec<int>{ e1, e2, f, g } ) ) continue;
                         if ( invl[e2] != f || invl[e1] != g ) continue;
                         int n1 = 0, n2 = 0;
                         for ( int mpass = 1; mpass <= 2; mpass++ )
                         for ( int pass = 1; pass <= 2; pass++ )
                         {    int f;
                              if ( mpass == 1 ) f = ( pass == 1 ? e1 : invl[e1] );
                              else f = ( pass == 1 ? e2 : invl[e2] );
                              for ( int j = 0; 
                                   j < (int) paths_indexl[f].size( ); j++ )
                              {    int64_t id = paths_indexl[f][j];
                                   const ReadPath& p = pathsl[id];
                                   for ( int l = 0; l < (int) p.size( ) - 1; l++ )
                                   {    if ( pass == 1 )
                                        {    if ( p[l] == h && p[l+1] == f )
                                             {    ( mpass == 1 ? n1 : n2 )++;
                                                  break;    }    }
                                        else
                                        {    if ( p[l] == f && p[l+1] == invl[h] )
                                             {    ( mpass == 1 ? n1 : n2 )++;
                                                  break;    }    }    }    }    }
                         if ( n1 > MAX_CAN_INS_DEL ) continue;
                         if ( n2 == 0 || n2 < MIN_CAN_INS_RATIO * n1 ) continue;
                         delsx.push_back( e1, invl[e1] );    }
                    UniqueSort(delsx);    }
               }
               */

               // Map the HyperBasevector edges back to the original.

               const int K = 48;
               ForceAssertEq( K, hbl.K( ) );
               vec<kmer<K>> starts( hbls.E( ) );
               for ( int e = 0; e < hbls.E( ); e++ )
                    starts[e].SetToSubOf( hbls.O(e), 0 );
               vec<int> sids( hbls.E( ), vec<int>::IDENTITY );
               SortSync( starts, sids );
               vec<vec<int>> backto( hbl.E( ) );
               vec<int> sto_right;
               hbls.ToRight(sto_right);
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
                                   break;    }    }    }    }

               // Form paths index.

               VecULongVec paths_indexl;
               invert( pathsl, paths_indexl, hbl.E( ) );

               // ------------  STAGE 6 ------------------------------
               // Form closures and convert to graph.  Note that scaffolding
               // does not use the PCR-free reads, which may not make sense.

               RecordMe(iclks[stage++]);
               vec<vec<int>> all_closuresl;
               vec<Bool> dupl( pathsl.size( )/2, False );
               vec<Bool> badl( pathsl.size( )/2, False );
               // NOTE EXPENSIVE CONVERSION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               ReadPathVecX pathslx;
               HyperBasevectorX hblx(hbl);
               pathslx.parallel_append(pathsl,hblx);
               // This sets graph properties = 200mer graph.
               // Also has fragmentation as it doesn't include other reads
               MakeClosures( hblx, invl, pathslx, paths_indexl, 
                    dupl, badl, all_closuresl, False );
               /*
               vec<Bool> cdel( all_closuresl.size( ), False );
               for ( int i = 0; i < all_closuresl.isize( ); i++ )
               for ( int j = 0; j < all_closuresl[i].isize( ); j++ )
                    if ( BinMember( delsx, all_closuresl[i][j] ) ) cdel[i] = True;
               EraseIf( all_closuresl, cdel );
               */
               Destroy(paths_indexl);

               // Make graph from closures.

               ClosuresToGraph(hblx,invl,all_closuresl,Dl,dinvl,False,"",CANON);
               // TODO:
               Cleaner( hblx, invl, pathsl, dupl, Dl, dinvl, dpathsl, False );

               // ------------  STAGE 7 ------------------------------
               RecordMe(iclks[stage++]);
               const Bool SINGLE = False; // ??????????????????????????????????
               PlaceReads( hblx, pathsl, dupl, Dl, dpathsl, False, SINGLE );
               VecIntVec ebcxl( hbl.E( ) );
               // populate esb2 with barcodes landing on e in various instances
               {    vec<vec<int>> esb2( es.size( ) );
                    for ( int i = 0; i < esb.isize( ); i++ )
                    {    int j;
                         for ( j = i + 1; j < esb.isize( ); j++ )
                              if ( esb[j].first != esb[i].first ) break;
                         int e = esb[i].first;
                         int p = BinPosition( es, e );
                         if ( p >= 0 ) // TODO: no need for this check
                         {    int rp = BinPosition( es, inv[e] );
                              for ( int k = i; k < j; k++ )
                              {    esb2[p].push_back( esb[k].second );
                                   esb2[rp].push_back( esb[k].second );    
                                        }    }
                         i = j - 1;    }
                    // At this point the hbl has been trimmed so may not contain
                    // all of es / esb
                    for ( int e = 0; e < hbl.E( ); e++ )
                    {    vec<int> x;
                        // backto translates current edge in hbl into edge(s) in es
                        // after all the trimmings to hbl
                        // then we associate uniq bc to e in hbl
                         for ( int i = 0; i < backto[e].isize( ); i++ )
                         {    int f = backto[e][i];
                              for ( auto b : esb2[f] ) x.push_back(b);    }
                         UniqueSort(x);
                         ebcxl[e].resize( x.size( ) );
                         for ( int j = 0; j < x.isize( ); j++ )
                              ebcxl[e][j] = x[j];    }    }

               // ------------  STAGE 8 ------------------------------
               // Create a vector of integers, one for each read, such that
               // "having two" nonzero elements is enough. These elements
               // indicate modified barcode ids - antiquated
               RecordMe(iclks[stage++]);

               DataSet dl;
               dl.dt = ReadDataType::BAR_10X;
               dl.start = 0;
               vec<DataSet> datasetsl = {dl};
               vec<int64_t> bidl( pathsl.size( ) );
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

               // Kill inversion bubbles.  For unknown reasons this causes
               // a crash.  Should work with a call to PlaceReads added at end,
               // could try later.

               vec<int> dels;
               /*
               ZapInversionBubbles( Dl, dinvl, dels );
               Dl.DeleteEdges(dels);
               RemoveUnneededVertices( Dl, dinvl );    
               CleanupCore( Dl, dinvl );
               */

               // ------------  STAGE 9 ------------------------------
               // Delete very weak inversion bubbles.
               // First calculate dpaths_index and dlens
               // For each vertex v, identify bubble topolgy:
               //            ---f--->
               //  --d---> v          w -i[d]->
               //            --i[f]->
               //  If the bubble is short (len[f] small)
               //  and if the number of reads containing both d and i[d] is insufficient
               //  then delete the bubble edges f and i[f]
               //  TODO: this introduces gaps and terminates lines!
               RecordMe(iclks[stage++]);

               const int MAX_INV_DIST = 200;
               const int MIN_INV_BRIDGES = 2;
               {    
                    IntIndex dpaths_indexl( dpathsl, Dl.E( ) );

                    vec<int> lens( Dl.E( ), 0 );
                    for ( int e = 0; e < Dl.E( ); e++ )
                    {    if ( Dl.O(e)[0] < 0 ) continue;
                         for ( int j = 0; j < Dl.O(e).isize( ); j++ )
                              lens[e] += hbl.Kmers( Dl.O(e)[j] );    }

                    for ( int v = 0; v < Dl.N( ); v++ )
                    {    if ( Dl.From(v).size( ) != 2 ) continue;
                         if ( Dl.From(v)[0] != Dl.From(v)[1] ) continue;
                         int w = Dl.From(v)[0];
                         if ( v == w ) continue;
                         if ( !Dl.To(v).solo( ) || !Dl.From(w).solo( ) ) continue;
                         int d1 = Dl.ITo(v,0), d2 = Dl.IFrom(w,0);
                         if ( dinvl[d1] != d2 ) continue;
                         int f1 = Dl.IFrom(v,0), f2 = Dl.IFrom(v,1);
                         if ( dinvl[f1] != f2 ) continue;
                         if ( lens[f1] > MAX_INV_DIST ) continue;
                         int nbridges = 0;
                         for ( int j = 0; j < dpaths_indexl.Count(d1); j++ )
                         {    int64_t id = dpaths_indexl.Val(d1,j);
                              const ReadPath& p = dpathsl[id];
                              for ( int l = 0; l < (int) p.size( ) - 2; l++ )
                              {    if ( p[l] == d1 && p[l+2] == d2 )
                                   {    nbridges++;
                                        break;    }    }    }
                         if ( nbridges >= MIN_INV_BRIDGES ) continue;
                         dels.push_back( f1, f2 );    }    }
               Dl.DeleteEdges(dels);
               RemoveUnneededVertices( Dl, dinvl );    
               CleanupCore( Dl, dinvl );

               // Disengage certain rare terminal loops.

               /*
               vec<int> to_left, to_right;
               Dl.ToLeft(to_left), Dl.ToRight(to_right);
               for ( int v = 0; v < Dl.N( ); v++ )
               {    if ( !Dl.From(v).solo( ) || Dl.To(v).size( ) != 2 ) continue;
                    int w = Dl.From(v)[0];
                    if ( !Dl.From(w).solo( ) || !Dl.To(w).solo( ) ) continue;
                    if ( Dl.From(w)[0] != v ) continue;
                    int d = Dl.IFrom(w,0);
                    int rd = dinvl[d];
                    if ( d == rd ) continue;
                    int N = Dl.N( );
                    Dl.AddVertices(2);
                    Dl.GiveEdgeNewToVx( d, v, N );
                    int rv = to_left[rd];
                    Dl.GiveEdgeNewFromVx( rd, rv, N+1 );    }
               RemoveUnneededVertices( Dl, dinvl );    
               CleanupCore( Dl, dinvl );
               */

               // ------------  STAGE 10 -----------------------------
               // Make scaffolds in two pass. After first scaffolding pass
               // in the second pass, delete a non cyclic edge d if it is short
               // and the l and r vertices are terminal.
               RecordMe(iclks[stage++]);

               PlaceReads( hblx, pathsl, dupl, Dl, dpathsl, False, SINGLE );
               for ( int pass = 1; pass <= 2; pass++ )
               {    
                    vec<int> to_left, to_right;
                    Dl.ToLeft(to_left), Dl.ToRight(to_right);

                    vec<int> lens( Dl.E( ), 0 );
                    for ( int e = 0; e < Dl.E( ); e++ )
                    {    if ( Dl.O(e)[0] < 0 ) continue;
                         for ( int j = 0; j < Dl.O(e).isize( ); j++ )
                              lens[e] += hbl.Kmers( Dl.O(e)[j] );    }
                    
                    if ( pass == 2 )
                    {    vec<int> dels;
                         for ( int d = 0; d < Dl.E( ); d++ )
                         {    int v = to_left[d], w = to_right[d];
                              if ( v < 0 ) continue;
                              if ( Dl.To(v).nonempty( ) || !Dl.From(v).solo( ) )
                                   continue;
                              if ( Dl.From(w).nonempty( ) || !Dl.To(w).solo( ) )
                                   continue;
                              if ( v == w ) continue;
                              if ( lens[d] >= 320 ) continue;
                              dels.push_back(d);    }
                         Dl.DeleteEdges(dels);
                         CleanupCore( Dl, dinvl );
                         RemoveUnneededVertices( Dl, dinvl );
                         PlaceReads( hblx,
                              pathsl, dupl, Dl, dpathsl, False, SINGLE );    }
                    String link_report;
                    Scaffold( hblx, invl, ebcxl, Dl, dinvl, dpathsl,
                         bidl, datasetsl, False, link_report, SINGLE );    }
               // Destroy(pathsl);
               // Destroy(dpathsl);

               // ------------  STAGE 11 -----------------------------
               // Delete components having no edge of size at least 
               // MIN_LINE_TO_WALK.
               // TODO: this is dangerous!
               RecordMe(iclks[stage++]);

               {    vec<vec<int>> comp;
                    Dl.ComponentsE(comp);

                    vec<int> lens( Dl.E( ), 0 );
                    for ( int e = 0; e < Dl.E( ); e++ )
                    {    if ( Dl.O(e)[0] < 0 ) continue;
                         for ( int j = 0; j < Dl.O(e).isize( ); j++ )
                              lens[e] += hbl.Kmers( Dl.O(e)[j] );    }

                    vec<int> dels;
                    for ( int i = 0; i < comp.isize( ); i++ )
                    {    int n = 0;
                         for ( int j = 0; j < comp[i].isize( ); j++ )
                         {    int d = comp[i][j];
                              n = Max( n, lens[d] );    }
                         if ( n < MIN_LINE_TO_WALK ) dels.append( comp[i] );    }
                    Dl.DeleteEdges(dels);
                    CleanupCore( Dl, dinvl );
                    RemoveUnneededVertices( Dl, dinvl );    }

               // ------------  STAGE 12 -----------------------------
               // Kill inversion artifacts.
               RecordMe(iclks[stage++]);
          
               PlaceReads( hblx, pathsl, dupl, Dl, dpathsl, False, SINGLE );
               dels.clear( );
               {    IntIndex dpaths_indexl( dpathsl, Dl.E( ) );
                    const int MAX_CAN_INS_DEL = 5;
                    KillInversionArtifacts( hblx, Dl, dinvl, dpathsl, 
                         dpaths_indexl, bidl, dels, MAX_CAN_INS_DEL, False );    }
               Dl.DeleteEdges(dels);
               RemoveUnneededVertices( Dl, dinvl );    
               CleanupCore( Dl, dinvl );

               // ------------  STAGE 13 -----------------------------
               // Pull apart stuff.

               /*
               PullApartInversions( Dl, dinvl );
               LineLevelPullApart( bcl, hblx, invl, pathsl, dup, Dl, dinvl );
               dels.clear( );
               PlaceReads( hblx, pathsl, dupl, Dl, dpathsl, False, SINGLE );
               PullApart( invl, Dl, dinvl, dpathsl, dels, False );
               Dl.DeleteEdges(dels);
               RemoveUnneededVertices( Dl, dinvl );    
               CleanupCore( Dl, dinvl );
               */

               RecordMe(iclks[stage++]);
               BigPull( hblx, invl, Dl, dinvl, pathsl, dupl, False );



               // =================================================================
               // =================================================================

               // New stuff.

               if (NEW_STUFF)
               {
                    // Glue dead ends that have long perfect overlaps.

                    vec< triple<int,int,int> > ends;
                    for ( int v = 0; v < Dl.N( ); v++ )
                    {    if ( Dl.From(v).empty( ) && Dl.To(v).solo( ) )
                         {    int d = Dl.ITo(v,0);
                              if ( Dl.O(d)[0] >= 0 )
                              {    int e = Dl.O(d).back( );
                                   ends.push( e, 1, d );    }    }
                         if ( Dl.To(v).empty( ) && Dl.From(v).solo( ) )
                         {    int d = Dl.IFrom(v,0);
                              if ( Dl.O(d)[0] >= 0 )
                              {    int e = Dl.O(d).front( );
                                   ends.push( e, 2, d );    }    }    }
                    Sort(ends);
                    vec<int> dels, to_left, to_right;
                    Dl.ToLeft(to_left), Dl.ToRight(to_right);
                    for ( int i = 0; i < ends.isize( ); i++ )
                    {    int j;
                         for ( j = i + 1; j < ends.isize( ); j++ )
                              if ( ends[j].first != ends[i].first ) break;
                         if ( j - i == 2 && ends[i].second == 1
                              && ends[i+1].second == 2 )
                         {    int d1 = ends[i].third, d2 = ends[i+1].third;
                              int rd1 = dinvl[d1], rd2 = dinvl[d2];
                              if ( IsUnique( vec<int>{d1,d2,rd1,rd2} )
                                   && make_pair(d1, d2) < make_pair(rd2, rd1) )
                              {    vec<int>& x = Dl.OMutable(d1);
                                   for ( int l = 1; l < Dl.O(d2).isize( ); l++ )
                                        x.push_back( Dl.O(d2)[l] );
                                   dels.push_back( d2, rd2 );
                                   int v = to_right[d1], w = to_right[d2];
                                   Dl.GiveEdgeNewToVx( d1, v, w );
                                   Dl.OMutable(rd1) = x;
                                   Dl.OMutable(rd1).ReverseMe( );
                                   for ( auto& e : Dl.OMutable(rd1) ) e = invl[e];
                                   int rv = to_left[rd1], rw = to_left[rd2];
                                   Dl.GiveEdgeNewFromVx( rd1, rv, rw );    }    }
                         i = j - 1;    }
                    Dl.DeleteEdges(dels);
                    RemoveUnneededVertices( Dl, dinvl );    
                    CleanupCore( Dl, dinvl );    

                    // Delete hanging ends.

                    HangBeGone( hblx, Dl, dinvl, 10, 1000, 0 );    }

               // =================================================================
               // =================================================================



               // ------------  STAGE 14 -----------------------------
               RecordMe(iclks[stage++]);
               // Experimental remove components with lines having low copy number
               // This should replace the subsequent block of code
               /* PlaceReads( hblx, pathsl, dupl, Dl, dpathsl, False, SINGLE ); */
               /* { */   
               /*     cout << "BEGIN CNV" << endl; */
               /*     vec<vec<pair<int,int>>> lbpl; */
               /*     MasterVec<SerfVec<pair<int,int>>> lbpxl; */
               /*     vec<double> cov; */
               /*     IntIndex dpaths_indexl1(dpathsl,Dl.E()); */
               /*     vec<vec<vec<vec<int>>>> dlinesl; */
               /*     vec<int> kmersl( hblx.E( ) ), llensl; */
               /*     for ( int e = 0; e < hblx.E( ); e++ ) */
               /*          kmersl[e] = hblx.Kmers(e); */
               /*     FindLines( Dl, dinvl, dlinesl, MAX_CELL_PATHS, MAX_CELL_DEPTH ); */
               /*     // suspicious!!!!!!!!! */
               /*     BarcodePos(bcl,hblx,Dl,dinvl,dpathsl,dpaths_indexl1,dlinesl,lbpl,0); */
               /*     for(auto x : lbpl) */
               /*     {   SerfVec<pair<int,int>> y(x.size()); */
               /*         for (int j = 0; j < x.isize(); j++) */
               /*             y[j] = x[j]; */
               /*         lbpxl.push_back(y);   } */
               /*     GetLineLengths (hblx, Dl, dlinesl, llensl); */
               /*     LineCN(kmersl,lbpxl, Dl, dlinesl, llensl, cov); */   
               /*     cout << "END CNV" << endl; */
               /* } */
                       

               // ------------  STAGE 15 -----------------------------
               // Remove components supported by very few barcodes.
               // where support is calculated by total unique bc landing on it.
               // This requires to have knowledge of copynumber and fine control on MIN_BC
               // TODO: needs inspection!
               RecordMe(iclks[stage++]);

               PlaceReads( hblx, pathsl, dupl, Dl, dpathsl, False, SINGLE );
               vec<vec<int>> comp;
               Dl.ComponentsEFast(comp);
               dels.clear( );
               const int MIN_BC = 30;
               {    IntIndex dpaths_indexl( dpathsl, Dl.E( ) );
                    for ( auto& x : comp )
                    {    vec<int> b;
                         for ( auto d : x )
                         for ( int pass = 1; pass <= 2; pass++ )
                         {    int g = ( pass == 1 ? d : dinvl[d] );
                              for ( int i = 0; i < dpaths_indexl.Count(g); i++ )
                              {    int64_t id = dpaths_indexl.Val( g, i );
                                   b.push_back( bcl[id] );    }    }
                         UniqueSort(b);
                         if ( b.isize( ) < MIN_BC )
                         {    for ( auto d : x )
                                   dels.push_back( d, dinvl[d] );    }    }    }
               Dl.DeleteEdges(dels);
               RemoveUnneededVertices( Dl, dinvl );    
               CleanupCore( Dl, dinvl );

               // ------------  STAGE 16 -----------------------------
               RecordMe(iclks[stage++]);
               // Remove lines that fail a diploid filter.
               DiploidFilter( hblx, Dl, dinvl );

               // Remove very small components.
               RemoveVerySmallComponents( hblx, Dl, dinvl, 200 );

               // ------------  STAGE 17 -----------------------------
               RecordMe(iclks[stage++]);
               //[C2G ---------------------------
               if(C2G){
                     RemoveUnneededVertices( Dl, dinvl );
                     CleanupCore( Dl, dinvl );
                     Validate( hblx, invl, Dl, dinvl );
               }
               //-------------------------------------------]
               
               // Write local assembly if requested.

               // TO DO, MAYBE:
               // basedir/a.fin.alignsb
               // builddir/a.sup.lhood
               // builddir/a.sup.lcov
               // builddir/a.sup.lbpx
               // builddir/a.sup.linelines

               if(STORE_LOCALS){
                    String zstr = String(std::to_string(z));
                    String lldir = udir + "/" + INSTANCE + "/locas";
                    Mkdir777(lldir);
                    String ldir = lldir + "/" + zstr;
                    Mkdir777(ldir);
                    String basedir = ldir + "/a.base", datadir = ldir + "/data";
                    Mkdir777( datadir ), Mkdir777( basedir );
                    String builddir = basedir + "/build";
                    Mkdir777(builddir);
                    BinaryWriter::writeFile( basedir + "/a.hbx", hblx );
                    int K = hblx.K( );
                    Ofstream( out, basedir + "/a.k" );
                    out << K << endl;
                    BinaryWriter::writeFile( basedir + "/a.inv", invl );
                    vec<int> to_left, to_right;
                    hbl.ToLeft(to_left), hbl.ToRight(to_right);
                    BinaryWriter::writeFile( basedir + "/a.to_left", to_left );
                    BinaryWriter::writeFile( basedir + "/a.to_right", to_right );
                    BinaryWriter::writeFile( basedir + "/a.dup", dupl );
                    vec<int> kmersl( hbl.E( ) );
                    for ( int e = 0; e < hbl.E( ); e++ )
                         kmersl[e] = hbl.O(e).isize( ) - K + 1;
                    BinaryWriter::writeFile( basedir + "/a.kmers", kmersl );
                    ReadPathVecX pathsX;
                    pathsX.parallel_append(pathsl,hblx);
                    pathsX.WriteAll( basedir + "/a.pathsX" );
                    VecULongVec paths_indexl;
                    invert( pathsl, paths_indexl, hbl.E( ) );
                    paths_indexl.WriteAll( basedir + "/a.paths.inv" );
                    vecbasevector tigsx( hbl.E( ) );
                    for ( int e = 0; e < hbl.E( ); e++ ) tigsx[e] = hbl.O(e);
                    tigsx.WriteAll( basedir + "/a.fastb" );
                    vec<String> opts 
                         = { "sample", "genome.fastb", "genome_amb.fastb" };
                    for ( auto s : opts )
                    {    if ( IsRegularFile( DIR + "/../" + s ) )
                              Cp2( DIR + "/../" + s, ldir );    }
                    BinaryWriter::writeFile( builddir + "/a.sup", Dl );
                    BinaryWriter::writeFile( builddir + "/a.sup.inv", dinvl );
                    PlaceReads(hblx, pathsl, dupl, Dl, dpathsl, False, SINGLE);
                    dpathsl.WriteAll( builddir + "/a.dpaths" );
                    VecULongVec dpaths_indexl;
                    invert( dpathsl, dpaths_indexl, Dl.E( ) );
                    dpaths_indexl.WriteAll( builddir + "/a.dpaths.index" );
                    vec<vec<vec<vec<int>>>> dlinesl;
                    FindLinesFixed( 
                         Dl, dinvl, dlinesl, MAX_CELL_PATHS, MAX_CELL_DEPTH );
                    BinaryWriter::writeFile(builddir + "/a.sup.lines", dlinesl);
                    BinaryWriter::writeFile( basedir + "/a.readids", ids2 );
                    ReadPathVec pathsg( ids2.size( ) );
                    for ( int i = 0; i < ids2.isize( ); i++ )
                         pathsx.unzip( pathsg[i], hb, ids2[i] );
                    pathsg.WriteAll( basedir + "/a.pathsg" );

                    vecbasevector pathsgb( ids2.size( ) );
                    for ( int i = 0; i < ids2.isize( ); i++ )
                    {    if ( pathsg[i].size( ) > 0 ) 
                         {    vec<int> x;
                              for ( auto e : pathsg[i] ) x.push_back(e);
                              pathsgb[i] = hb.Cat(x);    }    }
                    pathsgb.WriteAll( basedir + "/a.pathsgb" );

                    vecbasevector basesb( ids2.size( ) );
                    VecPQVec qualsb( ids2.size( ) );
                    #pragma omp critical
                    {    for ( int i = 0; i < ids2.isize( ); i++ )
                         {    basesb[i] = bases[ ids2[i] ];
                              qualsb[i] = quals[ ids2[i] ];    }    }
                    basesb.WriteAll( datadir + "/frag_reads_orig.fastb" );
                    qualsb.WriteAll( datadir + "/frag_reads_orig.qualp" );
                    vec<int64_t> bcil;
                    for ( int i = 0; i < bcl.isize( ); i++ )
                         if ( i == 0 || bcl[i] != bcl[i-1] ) bcil.push_back(i);
                    bcil.push_back( bcl.size( ) );
                    BinaryWriter::writeFile( 
                         datadir + "/frag_reads_orig.bci", bcil );
                    vecbasevector ac( all_closuresl.size( ) );
                    for ( int i = 0; i < all_closuresl.isize( ); i++ )
                         ac[i] = hbl.Cat( all_closuresl[i] );
                    ac.WriteAll( builddir + "/all_closures.fastb" );
                    vec< triple<int,int,int> > qeptl;
                    AllTinksCore( hblx, invl, pathsX, bcil, ebcxl, qeptl, 0 );
                    BinaryWriter::writeFile( basedir + "/a.bc_links", qeptl );
                    vec<int> llensl;
                    vec<vec<pair<int,int>>> lbpl;
                    {    IntIndex dpaths_indexl( dpathsl, Dl.E( ) );
                         BarcodePos( bcl, hblx, Dl, dinvl, 
                              dpathsl, dpaths_indexl, dlinesl, lbpl, 0, False );    }
    	            MasterVec<SerfVec<pair<int,int>>> lbpxl;
                    vec<double> covl;
                    for ( auto x : lbpl )
                    {    SerfVec<pair<int,int>> y( x.size( ) );
                         for ( int j = 0; j < x.isize( ); j++ )
                              y[j] = x[j];
                         lbpxl.push_back(y);    }
                    GetLineLengths( hblx, Dl, dlinesl, llensl );
                    BinaryWriter::writeFile( builddir + "/a.sup.llens", llensl );
                    LineCN( kmersl, lbpxl, Dl, dlinesl, llensl, covl, False );
                    lbpxl.WriteAll( builddir + "/a.sup.lbpx" );
                    BinaryWriter::writeFile( builddir + "/a.sup.lcov", covl );
                    vec< vec< pair<int,int> > > lhoodl;
                    LineProx( hblx, invl, ebcxl, Dl, dinvl, dlinesl, qeptl, lhoodl,
                         False );
                    BinaryWriter::writeFile( builddir + "/a.sup.lhood", lhoodl );
                    vec<vec< pair<int,int> >> linelocsl( kmersl.size( ) );
                    for ( int i = 0; i < dlinesl.isize( ); i++ )
                    for ( int j = 0; j < dlinesl[i].isize( ); j++ )
                    for ( int k = 0; k < dlinesl[i][j].isize( ); k++ )
                    for ( int l = 0; l < dlinesl[i][j][k].isize( ); l++ )
                    {    int d = dlinesl[i][j][k][l];
                         if ( Dl.O(d)[0] < 0 ) continue;
                         for ( int m = 0; m < Dl.O(d).isize( ); m++ )
                         {    int e = Dl.O(d)[m];
                              linelocsl[e].push( i, j );    }    }
                    if ( alignsb.size( ) > 0 )
                    {    MasterVec< SerfVec<triple<int,int,int> > > 
                              alignsbl( hbl.E( ) );
                         for ( int e = 0; e < hbl.E( ); e++ )
                         {    const vec<int>& x = backto[e];
                              int pos = 0;
                              vec<triple<int,int,int>> y;
                              for ( int j = 0; j < x.isize( ); j++ )
                              {    int f = es[ x[j] ];
                                   for ( auto m : alignsb[f] )
                                   {    y.push_back( make_triple( m.first,
                                             m.second - pos, m.third - pos ) );    }
                                   pos += hb.Kmers(f);    }
                              UniqueSort(y);
                              for ( auto t : y ) alignsbl[e].push_back(t);    }
                         alignsbl.WriteAll( basedir + "/a.alignsb" );
                         MasterVec< SerfVec< 
                              quad<int,Bool,ho_interval,ho_interval> > > 
                              viewl( dlinesl.size( ) );
                         #pragma omp parallel for
                         for ( int l = 0; l < dlinesl.isize( ); l++ )
                         {    View( l, hbl.K( ), kmersl, invl, Dl, dlinesl, 
                                   linelocsl, alignsbl, viewl[l] );    }
                         viewl.WriteAll( builddir + "/a.sup.view" );    }

               }

               // ------------  STAGE 18 -----------------------------
               // Now convert back to original coordinates.  So Dl and Dlp are
               // "identical" assemblies except for the edge space in which
               // they're defined.
               RecordMe(iclks[stage++]);
               
               DL[z] = Dl;
               digraphE<vec<int>>& Dlp = DL[z];
               for ( int e = 0; e < Dl.E( ); e++ )
               {    const vec<int>& x = Dl.O(e);
                    if ( IsCell(x) )
                    {    cell c;
                         c.CellDecode(x);
                         digraphE<vec<int>>& G = c.GMutable( );
                         for ( int d = 0; d < G.E( ); d++ )
                         {    vec<int> y;
                              for ( auto e : G.O(d) )
                              {    for ( int j = 0; 
                                        j < backto[e].isize( ); j++ )
                                   {    y.push_back( es[ backto[e][j] ] );    }
                                        }
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
               DLI[z] = dinvl;

               // ------------  STAGE 19 -----------------------------
               RecordMe(iclks[stage++]);
               // Clean up.
               // [ C2G --------------------------------------
               if(C2G==False){
                   RemoveUnneededVertices( DL[z], DLI[z] );
                   CleanupCore( DL[z], DLI[z] );
               }
               // -------------------------------------------]
               
               // Print stats.

               int64_t total_kmers = 0, base_edges = 0;
               for ( int j = 0; j < DL[z].E( ); j++ )
               {    const vec<int>& x = DL[z].O(j);
                    if ( x[0] < 0 ) continue;
                    base_edges += x.size( );
                    for ( auto e : x ) total_kmers += hb.Kmers(e);    }
               if (PRINT_DETAILS)
               {
                    #pragma omp critical
                    {    cout << "z = " << z << ", super edges = " << DL[z].E( )
                              << ", base edges = " << base_edges << ", kmers = "
                              << ToStringAddCommas(total_kmers) << endl;    }    }
               //[C2G---------------------------------------------------
               // Obtains BPATHS via placereads, read not aliased
               if(C2G){
                   BCLIST[z] = bcl;
                   UniqueSort(BCLIST[z]);
                   ReadPathVec dpl;
                   PlaceReads( hblx, pathsl, dupl, Dl, dpl, False, SINGLE );
                   // convert to bc placement
                   BPATHS[z].resize(BCLIST[z].size());
                   vec<vec<int>> dpbc(BCLIST[z].size());
                   for(int64_t r = 0; r< (int64_t) dpl.size(); r++){
                       int pos = BinPosition(BCLIST[z],bcl[r]);
                       for(auto d: dpl[r])
                           dpbc[pos].push_back(d);
                       UniqueSort(dpbc[pos]);
                   }
                   for(int b = 0; b< dpbc.isize(); b++){
                       BPATHS[z][b].reserve(dpbc[b].size());
                       for(auto d: dpbc[b])
                           BPATHS[z][b].push_back(d);
                   }
               }
               //----------------------------------------------------------------]
                   if (PRINT_DETAILS)
                   {
                        #pragma omp critical
                        {    cout << Date( ) << ": peak mem = "
                                  <<PeakMemUsageGBString( )<<
                                  ", mem = "<<MemUsageGBString( ) << endl;    }    }

               // Advance to next barcode group.
               RecordMe(iclks[stage++]);
               AccrueTime(clks,iclks,stage);
               totalstages = stage;
               }    }
     double timeused = WallClockTime() - clock;
     cout << Date( ) << ": done, time used = " << TimeSince(clock) << endl;
     cout << Date( ) << ": time spent in various stages (sorted by time): "<<endl;
     vec<int> stg( totalstages-1, vec<int>::IDENTITY );
     clks.resize(totalstages-1);
     SortSync(clks,stg);
     for(int s = 0; s < totalstages-1; s++)
         cout << "STAGE "<<stg[s]+1 <<" : "<<setprecision(4)<<fixed<<(100.0*clks[s]/nthreads/timeused)<<" %"<<endl;
     if (STORE_LOCALS)
     {    String lldir = udir + "/" + INSTANCE + "/locas";
          cout << "local assemblies: DIR = " << RealPath(lldir)
               << "/n/a.base [n = 0,...," << bs.size()-1 << "], "
               << "SUB = build" << endl;    }    }

void AddPairGaps( const HyperBasevectorX& hb, const vec<int> & inv,
     digraphE<vec<int>>& D, vec<int>& dinv, const ReadPathVec & dpaths,
     Bool verbose )
{
     vec<int> to_left, to_right;
     D.ToLeft( to_left );
     D.ToRight( to_right );

     cout << Date( ) << ": adding pair gaps" << endl;
     // Find hanging ends
     vec<int> lhangs, rhangs;
     for ( int v = 0; v < D.N( ); v++ ) {
          if ( D.To(v).size() == 0 && D.From(v).size() == 1 )
               lhangs.push_back( D.IFrom(v, 0) );
          else if ( D.From(v).size() == 0 && D.To(v).size() == 1 )
               rhangs.push_back( D.ITo(v, 0) );
     }
     
     VecULongVec dpaths_index;
     invert( dpaths, dpaths_index, D.E( ) );
     // Fix lhangs first
     for ( auto & he : lhangs ) {
          const int rhe = dinv[he];
          if ( he == rhe ) continue;
          vec <pair<int, int>> cand;
          for ( auto & rid : dpaths_index[rhe] ) {
               const int64_t prid = (rid % 2 == 0 ? rid+1: rid-1);
               auto & p = dpaths[prid];
               if ( p.size() ) {
                    const int d = p[0];
                    if ( d == he || d == rhe )
                         continue;
                    // postpone this
                    const int offset = p.getOffset();
                    if ( offset < 0 ) continue;
                    cand.push( d, -offset);
               }
          }
          UniqueSort( cand );
          vec <pair<int, int>> connects;
          for ( int i = 0; i < cand.isize(); i++ ) {
               int j = i+1;
               for ( ; j < cand.isize(); j++ ) {
                    if ( cand[j].first != cand[i].first )
                         break;
               }
               if ( j - i >= 2 ) {
                    const int e = cand[i].first;
                    const int offset = -cand[i].second;
                    const int numpairs = j - i;
                    if ( verbose ) {
                         cout << he << " -> " << e << "(offset = "
                              << offset << ", " << numpairs
                              << " pairs) "; }
                    connects.push( e, offset );
               }
               i = j;
          }
          if ( connects.size() && verbose )
               cout << endl;
          else
               cout << he << " NONE" << endl;

          // connect with pair gap
          // TODO: properly handle offsets
          for ( auto & eo : connects ) {
               const int d = eo.first, offset=eo.second;
               int E = D.E( );
               int src = to_left[d], snk = to_left[he];
               D.AddEdgeWithUpdate( src, snk, {-1}, to_left, to_right );
               const int rd = dinv[d], rhe = dinv[he];
               src = to_right[rhe];
               snk = to_right[rd];
               D.AddEdgeWithUpdate( src, snk, {-1}, to_left, to_right );
               dinv.push_back( E+1 );
               dinv.push_back( E );
          }

     }
}

void MergeBarcodeSets( vec<vec<int>>& bs, const vec<int64_t>& bci, 
     const int MAX_BARCODES )
{
     // Find barcode set intersections, then merge any two sets if one contains
     // least half of the other.

     cout << Date( ) << ": finding barcode set intersections" << endl;
     for ( int pass = 1; pass <= 2; pass++ )
     {    equiv_rel e( bs.size( ) );
          vec< pair<int,int> > bis;
          for ( int i = 0; i < bs.isize( ); i++ )
          for ( int j = 0; j < bs[i].isize( ); j++ )
               bis.push( bs[i][j], i );
          // cout << Date( ) << ": sorting bis" << endl;
          ParallelSort(bis);
          // cout << Date( ) << ": getting overs" << endl;
          vec< pair<int,int> > overs;
          for ( int64_t i = 0; i < bis.jsize( ); i++ )
          {    int64_t j;
               for ( j = i + 1; j < bis.jsize( ); j++ )
                    if ( bis[j].first != bis[i].first ) break;
               for ( int64_t k1 = i; k1 < j; k1++ )
               for ( int64_t k2 = k1 + 1; k2 < j; k2++ )
               {    overs.push( bis[k1].second, bis[k2].second );
                    overs.push( bis[k2].second, bis[k1].second );    }
               i = j - 1;    }
          // cout << Date( ) << ": sorting overs" << endl;
          ParallelSort(overs);
          for ( int64_t i = 0; i < overs.jsize( ); i++ )
          {    int64_t j = overs.NextDiff(i);
               int b1 = overs[i].first, b2 = overs[i].second;
               int n = j - i;
               int n1 = bs[b1].size( ), n2 = bs[b2].size( );
               const double MIN_OVER = 0.75;          
               double o = Max( double(n)/double(n1), double(n)/double(n2) );
               if ( o >= MIN_OVER ) e.Join( b1, b2 );
               i = j - 1;    }
          cout << Date( ) << ": found " << e.OrbitCount( ) << " orbits" << endl;
          vec<vec<int>> bs2;
          vec<int> reps;
          e.OrbitRepsAlt(reps);
          for ( int i = 0; i < reps.isize( ); i++ )
          {    vec<int> b, o;
               e.Orbit( reps[i], o );
               for ( auto x : o ) b.append( bs[x] );
               UniqueSort(b);
               if ( b.isize( ) <= MAX_BARCODES ) bs2.push_back(b);
               else for ( auto x : o ) bs2.push_back( bs[x] );    }
          bs = bs2;    }
     int64_t total = 0, total_reads = 0;
     for ( auto b : bs ) 
     {    total += b.size( );
          for ( auto x : b ) total_reads += bci[x+1] - bci[x];    }
     cout << Date( ) << ": total barcode sets = " << ToStringAddCommas( bs.size( ) )
          << endl;
     cout << Date( ) << ": in total " << ToStringAddCommas(total) << " barcodes" 
          << " and " << ToStringAddCommas(total_reads) << " reads" << endl;    }

void FillPairGaps( const HyperBasevectorX& hb, const vec<int>& inv,
     digraphE<vec<int>>& D, vec<int>& dinv, const ReadPathVecX& pathsx,
     const vec<int64_t>& bci, const vec<Bool>& dup )
{
     // Fill gaps using single reads that cross them.  There are multiple 
     // problems with this, including that it might collapse tandem 
     // duplications and that it might miss alleles.

     cout << Date( ) << ": filling pair gaps" << endl;
     vec<int> to_left, to_right;
     D.ToLeft(to_left), D.ToRight(to_right);
     vec<pair<int,int>> ends;
     vec<triple<int,int,int>> ds;
     for ( int d = 0; d < D.E( ); d++ )
     {    if ( !IsPairGap( D.O(d) ) || dinv[d] <= d ) continue;
          int v = to_left[d], w = to_right[d];
          if ( v == w ) continue;
          if ( !D.To(v).solo( ) || !D.From(w).solo( ) ) continue;
          int d1 = D.ITo(v,0), d2 = D.IFrom(w,0);
          if ( D.O(d1)[0] < 0 || D.O(d2)[0] < 0 ) continue;
          int e1 = D.O(d1).back( ), e2 = D.O(d2).front( );
          ends.push( e1, ds.size( ) );
          ds.push( d, e1, e2 );    }
     Sort(ends);
     // vec<int> b = Flatten(bs);
     double clock = WallClockTime( );
     vec<vec<vec<int>>> lexts( ds.size( ) );
     vec<int> lb( hb.E( ) + 1, -1 );
     lb.back( ) = ends.size( );
     for ( int i = ends.isize( ) - 1; i >= 0; i-- ) lb[ ends[i].first ] = i;
     for ( int i = hb.E( ) - 1; i >= 0; i-- ) if ( lb[i] < 0 ) lb[i] = lb[i+1];
     #pragma omp parallel for schedule( dynamic, 1 )
     for ( int i = 1; i < (int) bci.size( ) - 1; i++ )
     {    ReadPath p;
          vec<int> x, y;
          for ( int64_t id = bci[i]; id < bci[i+1]; id++ )
          {    if ( dup[id/2] ) continue;
               pathsx.unzip(p,hb,id);
               int n = p.size( );
               y.clear( );
               for ( auto e : p ) y.push_back(e);
               for ( int pass = 1; pass <= 2; pass++ )
               {    if ( pass == 2 )
                    {    y.ReverseMe( );
                         for ( int j = 0; j < n; j++ ) y[j] = inv[ y[j] ];    }
                    for ( int l = 0; l < n - 1; l++ )
                    {    int e = y[l];
                         int low = lb[e];
                         for ( int j = low; j < ends.isize( ); j++ )
                         {    if ( ends[j].first != e ) break;
                              int di = ends[j].second;
                              int e1 = ds[di].second, e2 = ds[di].third;
                              x.clear( );
                              Bool across = False;
                              int k;
                              for ( k = l + 1; k < n; k++ ) 
                              {    if ( y[k] == e2 )
                                   {    across = True;
                                        break;    }    }
                              if ( !across ) continue;
                              for ( int j = l + 1; j < k; j++ ) 
                                   x.push_back( y[j] );
                              #pragma omp critical
                              {    lexts[di].push_back(x);    }    }
                                   }    }    }    }
     vec<int> dels;
     for ( int i = 0; i < ds.isize( ); i++ )
     {    int d = ds[i].first, e1 = ds[i].second, e2 = ds[i].third;
          Sort(lexts[i]);
          // cout << endl;
          // PRINT4( i, d, e1, e2 );
          int count = 0;
          vec<vec<int>> alts;
          for ( int j = 0; j < lexts[i].isize( ); j++ )
          {    int k = lexts[i].NextDiff(j);
               // cout << "[" << ++count << "] " << printSeq( lexts[i][j] )
               //      << " (" << k-j << ")" << endl;
               alts.push_back( lexts[i][j] );
               j = k - 1;    }
          Bool have_empty = False;
          for ( auto& x : alts ) if ( x.empty( ) ) have_empty = True;
          if (have_empty) continue;
          if ( alts.size( ) == 1 || alts.size( ) == 2 )
          {    int v = to_left[d], w = to_right[d], rd = dinv[d];
               int rv = to_left[rd], rw = to_right[rd];
               dels.push_back( d, rd );
               for ( auto& x : alts ) 
               {    int E = D.E( );
                    D.AddEdge( v, w, x );
                    x.ReverseMe( );
                    for ( auto& e : x ) e = inv[e];
                    D.AddEdge( rv, rw, x );
                    dinv.push_back( E+1, E );    }    }    }
     D.DeleteEdges(dels);
     CleanupCore( D, dinv );    
     cout << Date( ) << ": " << TimeSince(clock) 
          << " used filling pair gaps" << endl;    }

// WARNING!  There are two versions of FormBarcodeSets, which need to be kept
// in sync.  It would be better to extract common code to reduce redundancy.

// FBSCore: nondeterministic order for outputs.

// WARNING!  Note that there is a bugfix, that is under a switch.

void FBSCore( const String& DIR, const HyperBasevectorX& hb,
     const digraphE<vec<int>>& D, const vec<int>& dinv,
     const vec<vec<vec<vec<int>>>>& dlines,
     const int MIN_LINE_TO_WALK, const int MIN_KMERS, const int MAX_BARCODES,
     const int MAX_VIEW, vec<vec<int>>& bs, vec<int>& bsl, const Bool FIX_BUG )
{
     // Compute ancillary stuff.

     cout << Date( ) << ": loading ebcx" << endl;
     VecIntVec ebcx( DIR + "/a.ebcx" );
     cout << Date( ) << ": forming barcode sets" << endl;
     vec<int> llens, linv, mult;
     GetLineLengths( hb, D, dlines, llens );
     LineInv( dlines, dinv, linv );
     ComputeMult( hb, D, mult );

     // Form barcode sets.  These are from the left end of each line l that is
     // long enough.  We get the other end via linv[l], but we only do both if
     // the line is long enough.

     cout << Date( ) << ": forming barcode sets" << endl;
     const int batch = 1000;
     bs.clear( ), bsl.clear( );
     #pragma omp parallel for
     for ( int bli = 0; bli < dlines.isize( ); bli += batch )
     {    vec<vec<int>> bsi;
          vec<int> bsli;
          for ( int li = bli; li < Min( bli + batch, dlines.isize( ) ); li++ )
          {    if ( llens[li] < MIN_LINE_TO_WALK ) continue;
               if ( linv[li] < li && llens[li] <= MAX_VIEW ) continue;
               const vec<vec<vec<int>>>& L = dlines[li];
               vec<int> b;
               int pos = 0;
               for ( int j = 0; j < L.isize( ); j++ )
               {    vec<int> lensj;
                    for ( int k = 0; k < L[j].isize( ); k++ )
                    {    int len = 0;
                         for ( int l = 0; l < L[j][k].isize( ); l++ )
                         {    int d = L[j][k][l];
                              if ( D.O(d)[0] < 0 ) continue;
                              for ( int m = 0; m < D.O(d).isize( ); m++ )
                              {    int e = D.O(d)[m];
                                   if ( mult[e] == 1 && hb.Kmers(e) >= MIN_KMERS 
                                        && (int) ebcx[e].size( ) <= MAX_BARCODES )
                                   {    for ( int l = 0; 
                                             l < (int) ebcx[e].size( ); l++ )
                                        {    b.push_back( ebcx[e][l] );    }    }
                                   len += hb.Kmers(e);
                                   if ( pos + len > MAX_VIEW ) break;    }    }
                         lensj.push_back(len);    }
                    Sort(lensj);
                    if ( lensj.nonempty( ) ) 
                    {    pos += Median(lensj);
                         if ( pos > MAX_VIEW ) break;    }    }

               // Save barcode set.

               UniqueSort(b);
               bsi.push_back(b), bsli.push_back(li);

               if (FIX_BUG)
               {    if ( llens[li] <= MAX_VIEW && li < linv[li] )
                    {    bsi.push_back(b), bsli.push_back( linv[li] );    }    }    }

          #pragma omp critical
          {    bs.append(bsi), bsl.append(bsli);    }    }    }

void FormBarcodeSets1( const String& DIR, const HyperBasevectorX& hb,
     const digraphE<vec<int>>& D, const vec<int>& dinv,
     const vec<vec<vec<vec<int>>>>& dlines, const vec<int64_t>& bci,
     const int MIN_LINE_TO_WALK, const int MIN_KMERS, const int MAX_BARCODES,
     const int MAX_VIEW, vec<vec<int>>& bs, 
     const MasterVec< SerfVec<triple<int,int,int> > >& alignsbf )
{
     double bclock = WallClockTime( );
     vec<int> bsl;
     FBSCore( DIR, hb, D, dinv, dlines, MIN_LINE_TO_WALK, MIN_KMERS, MAX_BARCODES,
          MAX_VIEW, bs, bsl );
     cout << Date( ) << ": sorting " << bs.size( ) << " barcode sets" << endl;
     ParallelSortSync( bs, bsl );
     DPRINT( bs.size( ) );
     cout << Date( ) << ": " << TimeSince(bclock)
          << " used creating barcode sets" << endl;
     MergeBarcodeSets( bs, bci, MAX_BARCODES );    }

void FormBarcodeSets2( const String& DIR, const HyperBasevectorX& hb,
     const vec<int>& inv, const digraphE<vec<int>>& D, const vec<int>& dinv,
     const vec<vec<vec<vec<int>>>>& dlines, const vec<int64_t>& bci,
     const int MIN_LINE_TO_WALK, const int MIN_KMERS, const int MAX_BARCODES,
     const int MAX_VIEW, vec<vec<int>>& bs, const vec<int>& fin,
     const String& REGION, const vec<int>& linelist,
     const MasterVec< SerfVec<triple<int,int,int> > >& alignsb,
     const MasterVec< SerfVec<triple<int,int,int> > >& alignsbf,
     const Bool PRINT_DETAILS )
{
     // Compute ancillary stuff.

     cout << Date( ) << ": compute ancillary stuff" << endl;
     int g_R = -1, start_R = -1, stop_R = -1;
     ParseRegion( REGION, g_R, start_R, stop_R );
     vec<int> llens, linv, mult;
     GetLineLengths( hb, D, dlines, llens );
     LineInv( dlines, dinv, linv );
     ComputeMult( hb, D, mult );

     // Form barcode sets.  These are from the left end of each line l that is
     // long enough.  We get the other end via linv[l], but we only do both if
     // the line is long enough.

     cout << Date( ) << ": forming barcode sets" << endl;
     double bclock = WallClockTime( );
     const int batch = 1000;
     vec<int> bsl;
     #pragma omp parallel for
     for ( int bli = 0; bli < dlines.isize( ); bli += batch )
     {    VirtualMasterVec<IntVec> ebcx( DIR + "/a.ebcx" );
          vec<vec<int>> bsi;
          vec<int> bsli;
          for ( int li = bli; li < Min( bli + batch, dlines.isize( ) ); li++ )
          {    if ( linelist.empty( ) )
               {    if ( llens[li] < MIN_LINE_TO_WALK ) continue;
                    if ( linv[li] < li && llens[li] <= MAX_VIEW ) continue;    }
               if ( linelist.nonempty( ) && !BinMember( linelist, li ) ) continue;
               const vec<vec<vec<int>>>& L = dlines[li];
               vec<int> b;

               // If FIN or REGION specified, pretest for hit, 
               // without looking at ebcx.
               
               if ( ( fin.nonempty( ) || REGION != "" ) && linelist.empty( ) )
               {    int pos = 0;
                    Bool hits_fin = False;
                    for ( int j = 0; j < L.isize( ); j++ )
                    {    vec<int> lensj;
                         for ( int k = 0; k < L[j].isize( ); k++ )
                         {    int len = 0;
                              for ( int l = 0; l < L[j][k].isize( ); l++ )
                              {    int d = L[j][k][l];
                                   if ( D.O(d)[0] < 0 ) continue;
                                   for ( int m = 0; m < D.O(d).isize( ); m++ )
                                   {    int e = D.O(d)[m];
                                        if (mult[e] == 1 && hb.Kmers(e) >= MIN_KMERS)
                                        {    for ( int pass = 1; pass <= 2; pass++ )
                                             {    int h = ( pass == 1 ? e : inv[e] );

                                                  if ( fin.nonempty( ) )
                                                  {
                                                  for ( int r = 0; 
                                                       r < (int) alignsbf[h].size( );
                                                       r++ )
                                                  {    int g = alignsbf[h][r].first;
                                                       if ( BinMember( fin, g ) ) 
                                                       {    hits_fin = True;   
                                                            goto finhit;    
                                                            }    }    
                                                  }
                                                  else
                                                  {
                                                  for ( int r = 0; 
                                                       r < (int) alignsb[h].size( );
                                                       r++ )
                                                  {    int g = alignsb[h][r].first;
                                                       int start 
                                                            = alignsb[h][r].second;
                                                       int stop 
                                                            = alignsb[h][r].third;
                                                       if ( g != g_R ) continue;
                                                       if ( IntervalOverlap(
                                                            start_R, stop_R,
                                                            start, stop ) > 0 )
                                                       {    hits_fin = True;   
                                                            goto finhit;    
                                                            }    }    
                                                  }

                                                       }    }
                                        len += hb.Kmers(e);
                                        if ( pos + len > MAX_VIEW ) break;    }    }
                              lensj.push_back(len);    }
                         Sort(lensj);
                         if ( lensj.nonempty( ) ) 
                         {    pos += Median(lensj);
                              if ( pos > MAX_VIEW ) break;    }    }
                    finhit: if ( !hits_fin ) continue; }

               // Now do full test.

               Bool hits_fin = False;
               int pos = 0;
               for ( int j = 0; j < L.isize( ); j++ )
               {    vec<int> lensj;
                    for ( int k = 0; k < L[j].isize( ); k++ )
                    {    int len = 0;
                         for ( int l = 0; l < L[j][k].isize( ); l++ )
                         {    int d = L[j][k][l];
                              if ( D.O(d)[0] < 0 ) continue;
                              for ( int m = 0; m < D.O(d).isize( ); m++ )
                              {    int e = D.O(d)[m];
                                   if ( mult[e] == 1 && hb.Kmers(e) >= MIN_KMERS 
                                        && (int) ebcx[e].size( ) <= MAX_BARCODES )
                                   {    if ( fin.size( ) > 0 )
                                        {    for ( int pass = 1; pass <= 2; pass++ )
                                             {    int h = ( pass == 1 ? e : inv[e] );

                                                  if ( fin.nonempty( ) )
                                                  {
                                                  for ( int r = 0; 
                                                       r < (int) alignsbf[h].size( );
                                                       r++ )
                                                  {    int g = alignsbf[h][r].first;
                                                       if ( BinMember( fin, g ) ) 
                                                            hits_fin = True;    }
                                                  }
                                                  else
                                                  {
                                                  for ( int r = 0; 
                                                       r < (int) alignsb[h].size( );
                                                       r++ )
                                                  {    int g = alignsb[h][r].first;
                                                       int start 
                                                            = alignsb[h][r].second;
                                                       int stop 
                                                            = alignsb[h][r].third;
                                                       if ( g != g_R ) continue;
                                                       if ( IntervalOverlap(
                                                            start_R, stop_R,
                                                            start, stop ) > 0 )
                                                       {    hits_fin = True;   
                                                            }    }    
                                                  }
                                                            }    }
                                        for ( int l = 0; 
                                             l < (int) ebcx[e].size( ); l++ )
                                        {    b.push_back( ebcx[e][l] );    }    }
                                   len += hb.Kmers(e);
                                   if ( pos + len > MAX_VIEW ) break;    }    }
                         lensj.push_back(len);    }
                    Sort(lensj);
                    if ( lensj.nonempty( ) ) 
                    {    pos += Median(lensj);
                         if ( pos > MAX_VIEW ) break;    }    }
               if ( fin.nonempty( ) && !hits_fin && linelist.empty( ) ) continue;

               // Save barcode set.

               UniqueSort(b);
               bsi.push_back(b), bsli.push_back(li);    }
          #pragma omp critical
          {    bs.append(bsi), bsl.append(bsli);    }    }
     
     cout << Date( ) << ": sorting " << bs.size( ) << " barcode sets" << endl;
     vec<int> bsls(bsl);
     Sort(bsls);
     if (PRINT_DETAILS)
          cout << Date( ) << ": using left ends of lines " << printSeq(bsls) << endl;
     ParallelSortSync( bs, bsl );
     DPRINT( bs.size( ) );
     cout << Date( ) << ": " << TimeSince(bclock)
          << " used creating barcode sets" << endl;
     MergeBarcodeSets( bs, bci, MAX_BARCODES );    }

void EvalVersusFinished( const HyperBasevectorX& hb, const vec<int>& inv,
     const digraphE<vec<int>>& D, const vec<int>& dinv, const String& FIN,
     const vec<int>& fin, const String& DIR, const String& OUTDIR, const Bool VERB1 )
{    
     cout << Date( ) << ": begin finished evaluation" << endl;
     vecbasevector G;
     String sample = FirstLineOfFile( DIR + "/../sample" );
     FetchFinished( sample, G, -1 );
     vec< vec< vec< triple< vec<int>, align, int > > > > Matches;
     vec<vec<vec<vec<int>>>> dlines;
     FindLines( D, dinv, dlines, MAX_CELL_PATHS_EVALUATION, MAX_CELL_DEPTH_EVALUATION );
     vecbasevector tigs( DIR + "/a.fastb" );
     MasterVec< SerfVec< triple< ho_interval, int, int > > > flocs;
     BinaryReader::readFile( DIR + "/a.fin.locs", &flocs );
     Bool VERBOSE = False;
     vec<vec<int>> stats;
     String suffix = FIN;
     if ( FIN.Contains( "[0,", 0 ) && FIN.Contains( ")", -1 ) )
          suffix = FIN.Between( ",", ")" );
     if ( FIN.Contains( "[0,", 0 ) && FIN.Contains( "]", -1 ) )
          suffix = ToString( FIN.Between( ",", "]" ).Int( ) + 1 );
     String result = FinAlign( suffix, fin, VERBOSE, OUTDIR, 
          hb, inv, tigs, D, dinv, dlines, G, flocs, stats );
     vec<String> res;
     Tokenize( result.After( "{" ).Before( "}" ), ',', res );
     StatLogger::log( "tot_edges",     res[0].Int(),  "Fin edges" );
     StatLogger::log( "tot_lines",     res[1].Int(),  "Fin lines" );
     StatLogger::log( "tot_cov_perc",  res[2].Double(),  "Total fin cov %" );
     StatLogger::log( "tot_errs",      res[3].Double(),  "Total errs" );
     StatLogger::log( "tot_cap",       res[4].Int(),  "Total cap gaps" );
     StatLogger::log( "tot_uncap",     res[5].Int(),  "Total uncap gaps" );
     StatLogger::setSilent(true);
     for ( int i = 0; i < fin.isize( ); i++ )
     {    String F = ToString(fin[i]);
          StatLogger::log( "edges_" + F, stats[i][0], "" );
          StatLogger::log( "lines_" + F, stats[i][1], "" );
          StatLogger::log( "cov_perc_" + F, (100.0*stats[i][2])/stats[i][3], "" );
          StatLogger::log( "errs_" + F, stats[i][4], "" );
          StatLogger::log( "cap_" + F, stats[i][5], "" );
          StatLogger::log( "uncap_" + F, stats[i][6], "" );    }
     StatLogger::setSilent(false); }

void PrintSummaryStats( const HyperBasevectorX& hb,
     const digraphE<vec<int>>& D, const vec<int>& dinv, ostream & out)
{ 
     vec<vec<int>> comp;
     D.ComponentsEFast(comp);
     int64_t nkmers = 0; 
     for ( int d = 0; d < D.E( ); d++ )
     {    if ( D.O(d)[0] >= 0 ) for ( auto e : D.O(d) ) nkmers += hb.Kmers(e);    }
     cout << "full assembly = " << D.E( ) << " edges in " << comp.size( ) 
          << " components" << ", comprising " << ToStringAddCommas(nkmers) 
          << " kmers" << endl;
     out  << "full assembly = " << D.E( ) << " edges in " << comp.size( ) 
          << " components" << ", comprising " << ToStringAddCommas(nkmers) 
          << " kmers" << endl;

     // log these stats

     StatLogger::log( "asm_edges", D.E( ), "Assembly edges" );
     StatLogger::log( "asm_comps", comp.size(), "Assembly components" );
     StatLogger::log( "asm_kmers", nkmers, "Assembly kmers" );  }   

void ExtractSubsetOfGlobal( const HyperBasevectorX& hb, const digraphE<vec<int>>& D,
     const vec<int>& dinv,
     const vec<vec<vec<vec<int>>>>& dlines, const vec<pair<int,int>>& tol,
     const MasterVec< SerfVec<triple<int,int,int> > >& alignsb,
     const MasterVec< SerfVec<triple<int,int,int> > >& alignsbf,
     const vec<int>& fin, const String& REGION, const Bool PRINT_FETCH,
     vec<Bool>& keep )
{
     keep.resize( D.E( ), False );
     const int MIN_KMERS = 5;
     int g_R = -1, start_R = -1, stop_R = -1;
     if ( REGION != "" )
     {    String chr = REGION.Before( ":" );
          if ( chr == "X" ) g_R = 22;
          else if ( chr == "Y" ) g_R = 23;
          else g_R = chr.Int( ) - 1;
          start_R = 1000000 * REGION.Between( ":", "-" ).Double( );
          stop_R = 1000000 * REGION.After( "-" ).Double( );    }
     vec<int> to_right;
     D.ToRight(to_right);
     #pragma omp parallel for
     for ( int d = 0; d < D.E( ); d++ )
     {    if ( D.O(d)[0] < 0 ) continue;
          Bool hits_fin = False;
          for ( int m = 0; m < D.O(d).isize( ); m++ )
          {    int e = D.O(d)[m];
               if ( hb.Kmers(e) >= MIN_KMERS )
               {    if ( fin.size( ) > 0 )
                    {    for ( int r = 0; r < (int) alignsbf[e].size( ); r++ )
                         {    int g = alignsbf[e][r].first;
                              if ( BinMember( fin, g ) ) hits_fin = True;    }    }
                    else if ( REGION != "" )
                    {    for ( int r = 0; r < (int) alignsb[e].size( ); r++ )
                         {    int g = alignsb[e][r].first;
                              int start = alignsb[e][r].second;
                              int stop = alignsb[e][r].third;
                              if ( g != g_R ) continue;
                              if ( IntervalOverlap(
                                   start_R, stop_R, start, stop ) > 0 )
                              {    hits_fin = True;    }    }    }    }
                                        }
          if (hits_fin)
          {    int l = tol[d].first, m = tol[d].second;
               const int FLANK = 3;
               if (PRINT_FETCH)
               {    cout << "using global L" << l << "." << m-FLANK << "-" 
                         << m+FLANK << " (half-open interval)" << endl;    }
               const vec<vec<vec<int>>>& L = dlines[l];
               int low = Max( 0, m - FLANK ), high = Min( L.isize( ), m + FLANK );
               #pragma omp critical
               {    for ( int r = low; r < high; r++ )
                    {    const vec<vec<int>>& M = L[r];
                         for ( int i = 0; i < M.isize( ); i++ )
                         for ( int j = 0; j < M[i].isize( ); j++ )
                         {    int d = M[i][j];
                              keep[d] = True;
                              keep[ dinv[d] ] = True;    }    
                         if ( r < high - 1 && r % 2 == 0 )
                         {    int d = L[r][0][0];
                              int v = to_right[d];
                              if ( D.From(v).solo( ) )
                              {    int g = D.IFrom(v,0);
                                   keep[g] = keep[ dinv[g] ] = True;
                                        }    }    }    }    }    }    }

void ParseRegion( const String& REGION, int& g_R, int& start_R, int& stop_R )
{    if ( REGION != "" )
     {    String chr = REGION.Before( ":" );
          if ( chr == "X" ) g_R = 22;
          else if ( chr == "Y" ) g_R = 23;
          else g_R = chr.Int( ) - 1;
          start_R = 1000000 * REGION.Between( ":", "-" ).Double( );
          stop_R = 1000000 * REGION.After( "-" ).Double( );    }    }

void PlaceReadsMasked2( const HyperBasevectorX& hb, const digraphE<vec<int>>& D,
     const vec<Bool>& dup, const ReadPathVecX& pathsx,
     const vec<Bool>& pmask, MasterVec<IntVec>& dpaths )
{
     // combine pmask and dup

     vec<Bool> md( dup.jsize( ) );
     if (pmask.size()) {
     for ( int64_t pid = 0; pid < dup.jsize(); pid++ )
          md[pid] = pmask[pid] || dup[pid];
     } else
          md = dup;

     int64_t to_place = 0;
     for ( auto noproc : md ) if ( !noproc ) to_place++;
     cout << Date( ) << ": placing " << ToStringAddCommas(to_place) 
          << " non-dup pairs on final assembly" << endl;
     const Bool ALIGN2 = False;
     vec<int64_t> maprr(pathsx.size());
     for(int64_t k = 0; k < (int64_t) pathsx.size(); k++)
         maprr[k] = k;
     PlaceReads2( hb, pathsx, maprr, md, D, dpaths, pathsx.size(), True, False, ALIGN2 );
     int64_t num_placed = 0;
     for ( auto & p : dpaths ) if ( p.size() > 0 ) num_placed++;
     cout << Date( ) << ": placed " << num_placed << " reads; " 
         << "peak mem = "<<PeakMemUsageGBString( )<<
         ", mem = "<<MemUsageGBString( ) << endl;
}

void PlaceReadsMasked( const HyperBasevectorX& hb, const digraphE<vec<int>>& D,
     const vec<Bool>& dup, const ReadPathVecX& pathsx,
     const vec<Bool>& pmask, ReadPathVec& dpaths )
{
     // combine pmask and dup

     vec<Bool> md( dup.jsize( ) );
     if (pmask.size()) {
     for ( int64_t pid = 0; pid < dup.jsize(); pid++ )
          md[pid] = pmask[pid] || dup[pid];
     } else
          md = dup;

     int64_t to_place = 0;
     for ( auto noproc : md ) if ( !noproc ) to_place++;
     cout << Date( ) << ": placing " << ToStringAddCommas(to_place) 
          << " non-dup pairs on final assembly" << endl;
     const Bool ALIGN2 = False;
     PlaceReads( hb, pathsx, md, D, dpaths, True, False, ALIGN2 );
     int64_t num_placed = 0;
     for ( auto & p : dpaths ) if ( p.size() > 0 ) num_placed++;
     cout << Date( ) << ": placed " << num_placed << " reads; " 
         << " peak mem = "<<PeakMemUsageGBString( )<<
         ", mem = "<<MemUsageGBString( ) << endl;   }

void PlaceLinkedReadsMasked( const HyperBasevectorX & hb, const vec<int> & inv,
     const digraphE<vec<int>> & D, const vec<int> & dinv, 
     const vec<vec<vec<vec<int>>>> & dlines, const vec<Bool> & dup,
     const vec<int64_t> & bci, const ReadPathVecX & pathsx, ReadPathVec & dpaths,
     const vec<Bool> & pmask, const int MAX_PASSES, const Bool verbose )
{
     // combine pmask and dup

     vec<Bool> md( dup.jsize( ) );
     if (pmask.size()) {
     for ( int64_t pid = 0; pid < dup.jsize(); pid++ )
          md[pid] = pmask[pid] || dup[pid];
     } else
          md = dup;

     int64_t to_place = 0;
     for ( auto noproc : md ) if ( !noproc ) to_place++;
     cout << Date( ) << ": placing " << ToStringAddCommas(to_place) 
          << " non-dup pairs on final assembly" << endl;
     const Bool ALIGN2 = False;
     
     PlaceLinkedReads( hb, inv, D, dinv, dlines, md, bci, 
          pathsx, dpaths, MAX_PASSES, verbose );

     int64_t num_placed = 0;
     for ( auto & p : dpaths ) if ( p.size() > 0 ) num_placed++;
     cout << Date( ) << ": placed " << num_placed << " reads; " 
         << " peak mem = "<<PeakMemUsageGBString( )<<
         ", mem = "<<MemUsageGBString( ) << endl;
}

