// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>

#include "CoreTools.h"
#include "ParallelVecUtilities.h"
#include "graph/DigraphTemplate.h"
#include "math/HoInterval.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "paths/long/large/Lines.h"
#include "10X/DfTools.h"
#include "10X/Gap.h"
#include "10X/Heuristics.h"
#include "10X/IntIndex.h"
#include "10X/LineOO.h"
#include "10X/PlaceReads.h"
#include "10X/PullApart.h"
#include "10X/Super.h"

void BucketLines( const vec<vec<vec<vec<int>>>>& dlines, const vec<int>& llens,
     vec<vec<int>>& buckets, const int min_len )
{    buckets.clear( );
     const int bsize = 1000000;
     int nd = dlines.size( );
     vec<int> llensx(llens), ids( nd, vec<int>::IDENTITY );
     ParallelReverseSortSync( llensx, ids );
     for ( int i = 0; i < nd; i++ )
     {    if ( llensx[i] < min_len ) continue;
          int64_t sum = llensx[i];
          int j;
          for ( j = i + 1; j < nd; j++ )
          {    if ( llensx[j] < min_len ) continue;
               if ( sum >= bsize ) break;
               sum += llensx[j];    }
          vec<int> b;
          for ( int k = i; k < j; k++ )
          {    if ( llensx[k] < min_len ) continue;
               b.push_back( ids[k] );    }
          buckets.push_back(b);
          i = j - 1;    }    }

void FindMoleculesOnLines(    vec<vec<vec<vec<int>>>>& dlines, 
                              vec <int> &llens,
                              vec<vec<triple<int32_t,int,int64_t>>>& lrpb,
                              vec<double> & fhist,
                              vec<pair<float, int>> & lr,
                              int & mol_lwml,
                              const int MAX_GAP_FRAG,
                              const int MIN_READS_FRAG,
                              const int MIN_FRAG_LEN )
{
     const int MIN_SAMPLES = 100000;
     const double BIN_WIDTH   = 1000.0; // 1kb bins for histogram
     vec<int> flens;
     // create an index that sorts lines, longest first.
     vec<int> llens_copy = llens, lindex( llens.size(), vec<int>::IDENTITY );
     ReverseSortSync(llens_copy, lindex);
     llens_copy.clear();
     int longest_scaffold=0;

     for ( int li = 0; li < dlines.isize( ); li++ ) {    
          const vec<triple<int32_t,int,int64_t>> & RPB = lrpb[lindex[li]];
          // line too short we can't really find any molecules
          if ( llens[lindex[li]] <= MIN_FRAG_LEN )
               break;

          // if line has too few reads then continue
          if ( RPB.isize() <= MIN_READS_FRAG )
               continue;
          
          for ( int i = 0; i < RPB.isize( ); i++ ) {    
               int j=i+1;
               set <int64_t> pids_seen;
               vec<int> rpos;
               pids_seen.insert(RPB[i].third);
               for ( j = i + 1; j < RPB.isize( ); j++ ) {
                    if ( RPB[j].first != RPB[i].first ) break;
                    int rs = RPB[j].second - RPB[j-1].second;
                    if ( rs > MAX_GAP_FRAG ) break;
                    if ( !pids_seen.count( RPB[j].third ) ) {
                         pids_seen.insert( RPB[j].third );
                         rpos.push_back( RPB[j].second );
                    }
               }
               // reads per molecule
               int rpm = rpos.size();
               if ( rpm >=  MIN_READS_FRAG ) {
                    int start = RPB[i].second, stop=RPB[j-1].second;
                    // diameter or max length-span of reads
                    int dia = stop - start;
                    int closest_edge = Max(Min(start, llens[li]-stop), 0);
                    int max_sep = -1;
                    for ( int k = 0; k != rpm-1; k++ ) {
                         int sep = rpos[k+1] - rpos[k];
                         if ( sep > max_sep )
                              max_sep = sep;
                    }
                    // length correction for sequence that isn't captured 
                    // in the length-span between the reads
                    int lest = int( (rpm + 1.0)/(rpm - 1.0)*dia );
                    if ( closest_edge > MAX_GAP_FRAG && lest > MIN_FRAG_LEN ) {
                         flens.push_back( lest );
                         lr.push_back( make_pair( lest/1000.0, rpm ) );
                    }
               }
               // move to start of next barcode
               i = j - 1;    
          }
          if ( flens.isize() > MIN_SAMPLES ) {
               longest_scaffold = llens[lindex[li]];
               break;
          }
     }
     Sort(flens);
     ReverseSort( lr );
     int N50f = ( flens.nonempty( ) ? N50(flens) : 0 );
     // compute molecule length histogram
     double max_frag_len = 0;
     int num_bins = 0;
     if ( !flens.empty( ) ) {
          max_frag_len = flens[flens.size()-1];
          num_bins = int(max_frag_len/BIN_WIDTH)+1;
          fhist.resize(num_bins, 0.0);
     }
     for ( auto f: flens ) {
          fhist[ int(f/BIN_WIDTH) ] += 1;
     }
     int sample = flens.size();
     double lwml=0.0, totlen=0.0;
     for ( int i = 0; i != sample; i++ ) {
          totlen += flens[i];
          lwml   += double(flens[i])*double(flens[i]);
     }
     lwml /= Max(totlen, 1.0);
     mol_lwml = int(lwml);
     // compute lwml from frag hist
     double l2tot=0.0, ltot=0.0;
     for ( int l = 0; l != num_bins; l++ ) {
          double bin_mp = (l+0.5)*BIN_WIDTH;
          l2tot += (fhist[l]*bin_mp*bin_mp);
          ltot  += (fhist[l]*bin_mp);
     }
     double lwml2 = l2tot/Max(ltot, 1.0);
     cout << Date( ) << ": N50/LWML/LWML from hist = " << ToStringAddCommas(N50f) << "/"
          << ToStringAddCommas(int(lwml)) << "/" << ToStringAddCommas(int(lwml2))
          << " [sample size = " << ToStringAddCommas( sample ) 
          << ", longest scaffold = " << ToStringAddCommas( longest_scaffold) << "]\n";
     
}

void LineCN( const vec<int>& kmers, const MasterVec<SerfVec<pair<int,int>>>& lbp,
     const digraphE<vec<int>>& D, const vec<vec<vec<vec<int>>>>& dlines, 
     const vec<int>& llens, vec<double>& COV )
{
     // Test for broken input.

     int num_lines = dlines.size();
     StatLogger::issue_alert("num_lines", num_lines);

     // Heuristics.

     const double FLANK_MUL = 1.1;
     const int FLANK_MIN = 1000;
     const int FLANK_MAX = 100000;
     const double BASE_CN = 2.0;
     const int CN_BOUND = 4;

     // Compute bounded physical coverage (approximate description).  For a given 
     // integer f > 0 and a given point p on a line of length >= 2f, let cov(p,f) be
     //  the number of barcodes that have a read placed in both the intervals 
     // [p-f,p) and (p,p+f].  Let cov(f) be the mean over all points p on the line 
     // (excluding the two ends of length f).

     cout << Date( ) << ": finding bridging" << endl;
     vec<int> flanks = {FLANK_MIN};
     while( flanks.back( ) < FLANK_MAX ) 
          flanks.push_back( flanks.back( ) * FLANK_MUL );
     int nf = flanks.size( );
     vec<int64_t> count( dlines.size( ) * nf, 0 );
     #pragma omp parallel for schedule(dynamic, 1)
     for ( int64_t l = 0; l < dlines.isize( ); l++ )
     {    if ( llens[l] <= 2*flanks[0] ) continue;
          const SerfVec< pair<int,int> >& x = lbp[l];
          for ( int i = 0; i < (int) x.size( ); i++ )
          {    int j;
               for ( j = i + 1; j < (int) x.size( ); j++ )
                    if ( x[j].first != x[i].first ) break;
               for ( int f = 0; f < nf; f++ )
               {    int flank = flanks[f];
                    for ( int k = i + 1; k < j; k++ )
                    {    int left = x[k-1].second, right = x[k].second;
                         if ( right - left >= 2*flank ) continue;
                         int extra = right - left - flank;
                         if ( extra > 0 ) { left += extra, right -= extra; }
                         left = Max( left, flank );
                         right = Min( right, llens[l] - flank );
                         if ( right <= left ) continue; // can it happen?
                         count[l*nf+f] += right - left;    }    }
               i = j - 1;    }    }    
     vec<int64_t> num( nf, 0 ), den( nf, 0 );
     #pragma omp parallel for schedule(dynamic,1)
     for ( int f = 0; f < nf; f++ )
     {    for ( int64_t l = 0; l < dlines.jsize( ); l++ ) 
          {    num[f] += count[l*nf+f];
               if ( llens[l] > 2*flanks[f] ) 
                    den[f] += llens[l] - 2*flanks[f];    }    }
     vec<double> cmean(nf);
     for ( int f = 0; f < nf; f++ ) cmean[f] = double(num[f])/double(den[f]);

     // Compute copy number of individual lines.  This is done by choosing a 
     // flank and computing the normalized bounded physical coverage relative to
     // that flank.

     cout << Date( ) << ": computing copy number" << endl;
     COV.resize_and_set( dlines.size( ), -1 );
     #pragma omp parallel for schedule(dynamic, 1000)
     for ( int64_t l = 0; l < dlines.jsize( ); l++ )
     {    if ( llens[l] <= 2*flanks[0] ) continue;
          int f;
          for ( f = 1; f < nf; f++ ) if ( llens[l] <= CN_BOUND*flanks[f] ) break;
          f--;
          double c = double(count[l*nf+f]) / double(llens[l]-2*flanks[f]);
          COV[l] = BASE_CN * c / cmean[f];    }

     // Assay behavior, informational only.

     for ( int pass = 1; pass <= 3; pass++ )
     {    int MIN_REP, MAX_REP;
          if ( pass == 1 ) { MIN_REP = 20000, MAX_REP = 50000; }
          if ( pass == 2 ) { MIN_REP = 50000, MAX_REP = 100000; }
          if ( pass == 3 ) { MIN_REP = 100000, MAX_REP = 1000000000; }
          const double CN_MULT = 1.5;
          vec<double> covs;
          int count = 0;
          for ( int i = 0; i < dlines.isize( ); i++ )
          {    if ( llens[i] < MIN_REP || llens[i] >= MAX_REP ) continue;
               if ( COV[i] > BASE_CN * CN_MULT ) continue;
               if ( COV[i] < BASE_CN / CN_MULT ) continue;
               count++;
               covs.push_back( COV[i] / BASE_CN );    }
          if ( covs.nonempty( ) )
          {    double dev = StdDev( covs, 1.0 );
               cout << Date( ) << ": copy number dev for lines ";
               if ( pass == 1 ) cout << "20-50 kb";
               if ( pass == 2 ) cout << "50-100 kb";
               if ( pass == 3 ) cout << "100+ kb";
               cout << " = " << PERCENT_RATIO( 3, dev*1000, 1000 )
                    << " [count=" << count << "]" << endl;    }    }    }


void FixMisassemblies( const HyperBasevectorX& hb, const vec<Bool>& dup,
     const vec<int64_t>& bci, const ReadPathVecX& paths,
     digraphE<vec<int>>& D, vec<int>& dinv, ReadPathVec& dpaths )
{
     cout << Date( ) << ": begin fixing misassemblies" << endl;

     // Place reads.

     PlaceReads( hb, paths, dup, D, dpaths, True, False );
     cout << Date( ) << ": making index" << endl;
     IntIndex dpaths_index( dpaths, D.E( ), False );

     // Expand barcode index.

     cout << Date( ) << ": expanding barcode index" << endl;
     vec<int32_t> bc( bci.back( ), -1 );
     #pragma omp parallel for
     for ( int b = 0; b < bci.isize( ) - 1; b++ )
     {    int64_t start = bci[b], stop = bci[b+1];
          for ( int64_t j = start; j < stop; j++ )
               bc[j] = b;    }

     // Zap weird inversion bubbles.

     vec<int> dels;
     ZapInversionBubbles( D, dinv, dels );
     D.DeleteEdges(dels);

     // Kill cells that look misassembled.

     const int BC_REQUIRE = 5000;
     const int BC_FLANK = 20000;
     const int BC_IGNORE = 2000;
     vec<vec<vec<vec<int>>>> dlines;
     vec <double> fhist;
     vec <pair<float, int>> lr;
     FindLines( D, dinv, dlines, MAX_CELL_PATHS, MAX_CELL_DEPTH );
     KillMisassembledCells( hb, dup, bci, paths, D, dinv, dpaths, dpaths_index,
          bc, dlines, dels, fhist, lr, BC_REQUIRE, BC_FLANK, BC_IGNORE, False);

     // Delete the edges.

     UniqueSort(dels);
     cout << Date( ) << ": deleting " << dels.size( )
          << " putatatively bad edges" << endl;
     D.DeleteEdges(dels);    }

void KillMisassembledCells( const HyperBasevectorX& hb, const vec<Bool>& dup,
     const vec<int64_t>& bci, const ReadPathVecX& paths,
     digraphE<vec<int>>& D, vec<int>& dinv, ReadPathVec& dpaths,
     IntIndex& dpaths_index, vec<int32_t>& bc,
     vec<vec<vec<vec<int>>>>& dlines, vec<int>& dels, vec <double>& fhist,
     vec <pair<float, int>> & lr,
     const int BC_REQUIRE, int BC_FLANK, int BC_IGNORE, const Bool verbose )
{
     cout << Date( ) << ": begin fixing misassembled cells" << endl;

     // Heuristics.

     
     const int BC_MIN = 10;
     const int BC_MAX_CELL = 1000;

     // Get superedge lengths.

     vec<int> dlens( D.E( ), 0 );
     #pragma omp parallel for
     for ( int e = 0; e < D.E( ); e++ )
     {    if ( D.O(e)[0] < 0 ) continue;
          for ( int j = 0; j < D.O(e).isize( ); j++ )
               dlens[e] += hb.Kmers( D.O(e)[j] );    }

          vec<int> to_left, to_right, llens;
     D.ToLeft(to_left), D.ToRight(to_right);
     GetLineLengths( hb, D, dlines, llens );
     
     
     
     // Estimate fragment length distribution.

     const int MAX_GAP_FRAG = 60000;
     const int MIN_READS_FRAG = 4;
     const int MIN_FRAG_LEN = 1000;
     int mol_lwml=0; 
     // Define local block so lrpb can be destroyed
     // by going out of scope
     {
          // compute read, pos on lines
          vec< vec<triple<int32_t, int, int64_t>>> lrpb;
          ReadPosLine( bc, hb, D, dinv, dpaths, dpaths_index, dlines, lrpb, 0 );

          FindMoleculesOnLines( dlines, llens, lrpb, fhist, lr, mol_lwml,
                                MAX_GAP_FRAG, MIN_READS_FRAG, MIN_FRAG_LEN);
     }
     cout << Date( ) << ": Line N50 = " << ToStringAddCommas(N50(llens)) << endl;

     // Scale down if parameters are not supported by the fragment size.

     BC_IGNORE = Min( BC_IGNORE, mol_lwml/4 );
     
     // Compute barcode positions.

     vec< vec<pair<int, int>>> lbp;
     BarcodePos( bc, hb, D, dinv, dpaths, dpaths_index, dlines, lbp, 0 );


     // Get mean number of positions in a window.

     int64_t total_bases = 0, total_pos = 0;
     for ( int i = 0; i < dlines.isize( ); i++ )
     {    if ( llens[i] < BC_FLANK ) continue;
          total_bases += llens[i];
          total_pos += lbp[i].size( );    }
     double winpos
          = double(BC_FLANK-BC_IGNORE) * double(total_pos) / double(total_bases);
     cout << Date( ) << ": " << winpos << " positions per window" << endl;

     // Define line buckets.

     vec<vec<int>> buckets;
     BucketLines( dlines, llens, buckets );

     // Loop across cells.

     cout << Date( ) << ": kill cells at misassemblies, main loop" << endl;
     #pragma omp parallel for schedule(dynamic, 1)
     for ( int bi = 0; bi < buckets.isize( ); bi++ )
     {    vec< pair<int,int> > orc;
          vec<int> lens, bleft, bright, es;
          for ( int ib = 0; ib < buckets[bi].isize( ); ib++ )
          {    int dl = buckets[bi][ib];
               const vec<vec<vec<int>>>& L = dlines[dl];
               orc.clear( );
               orc.resize( lbp[dl].size( ) );
               for ( int j = 0; j < orc.isize( ); j++ )
                    orc[j] = make_pair( lbp[dl][j].second, lbp[dl][j].first );
               Sort(orc);
               int pos = 0;
               for ( int m = 0; m < L.isize( ); m++ )
               {    const vec<vec<int>>& M = L[m];

                    // Get cell length.

                    int ncell = 0;
                    lens.clear( );
                    for ( int i = 0; i < M.isize( ); i++ )
                    {    const vec<int>& x = M[i];
                         int len = 0;
                         for ( int j = 0; j < x.isize( ); j++ )
                         {    int d = x[j];
                              if ( D.O(d)[0] >= 0 ) len += dlens[d];    }
                         lens.push_back(len);    }
                    Sort(lens);
                    if ( lens.nonempty( ) ) ncell = Median(lens);
                    if ( m % 2 == 0 )
                    {    pos += ncell;
                         continue;    }

                    // If it's a gap edge, find it's identifier.  Convoluted.

                    int g = -1;
                    if ( M.solo( ) && M[0].empty( ) && m > 0 && L[m-1].nonempty( )
                         && L[m-1][0].nonempty( ) )
                    {    int d = L[m-1][0][0];
                         int v = to_right[d];
                         if ( v >= 0 ) g = D.IFrom(v,0);    }

                    // Check for enough room.

                    int mid = pos + ncell/2;
                    if ( mid >= BC_REQUIRE && llens[dl] - mid >= BC_REQUIRE
                         && ncell <= BC_MAX_CELL )
                    {
                         // Compute number of bridging barcodes.

                         int64_t low = LowerBound1( orc, mid - BC_FLANK );
                         int64_t high = UpperBound1( orc, mid + BC_FLANK );
                         bleft.clear( ), bright.clear( );
                         for ( int64_t j = low; j < high; j++ )
                         {    int start = orc[j].first, b = orc[j].second;
                              if ( start >= mid - BC_FLANK
                                   && start <= mid - BC_IGNORE )
                              bleft.push_back(b);
                              if ( start <= mid + BC_FLANK
                                   && start >= mid + BC_IGNORE )
                              bright.push_back(b);    }
                         int nleft = bleft.size( ), nright = bright.size( );
                         int n = Min( nleft, nright );
                         UniqueSort(bleft), UniqueSort(bright);
                         int bridge = MeetSize( bleft, bright );
                         int expect = Min( 1.0, n/winpos ) * BC_MIN;

                         // Check for not enough.

                         es.clear( );
                         if ( bridge < expect || verbose )
                         {    if ( g >= 0 ) es.push_back(g);
                              for ( int i = 0; i < M.isize( ); i++ )
                                   es.append( M[i] );
                              UniqueSort(es);    }
                         if (verbose)
                         {
                              #pragma omp critical
                              {    cout << "cell {" << printSeq(es) << "} from line "
                                        << dl << " has support " << bridge
                                        << ", expect " << expect;
                                   if ( bridge < expect ) cout << " [DELETE]";
                                   cout << endl;    }    }
                         if ( bridge < BC_MIN )
                         {    int nes = es.size( );
                              for ( int i = 0; i < nes; i++ )
                                   es.push_back( dinv[ es[i] ] );
                              UniqueSort(es);
                              #pragma omp critical
                              {    dels.append(es);    }    }    }

                    // Advance position.

                    pos += ncell;    }    }    }    }

void Cleaner( const HyperBasevectorX& hb, const vec<int>& inv,
     const ReadPathVecX& paths, const vec<Bool>& dup, digraphE< vec<int> >& D,
     vec<int>& dinv, ReadPathVec& dpaths, const Bool verbose )
{
     // Place reads on the graph.

     if (verbose) cout << Date( ) << ": placing reads on supergraph" << endl;
     const Bool single = False;
     PlaceReads( hb, paths, dup, D, dpaths, verbose, single );

     // Index the paths.

     if (verbose) cout << Date( ) << ": indexing paths" << endl;
     IntIndex dpaths_index( dpaths, D.E( ), verbose );

     // Index D.

     if (verbose) cout << Date( ) << ": creating index" << endl;
     vec<vec<pair<int,int>>> nd( hb.E( ) );
     for ( int e = 0; e < D.E( ); e++ )
     {    const vec<int>& x = D.O(e);
          if ( x[0] < 0 ) continue;
          for ( int j = 0; j < x.isize( ); j++ )
               nd[ x[j] ].push( e, j );    }

     // Heuristics.  Note that MAX_END was 50, and was changed to 250.  This is
     // kind of covering up stuff, and could be losing power.
     // Ditto for MAX_TINY, which was 10.

     const int MAX_KILL = 320;
     const double MIN_RATIO = 25;
     const int MIN_EXT = 500;
     const int MAX_TINY = 150;
     const int MIN_STANDALONE = 200;

     // Find lengths.

     if (verbose) cout << Date( ) << ": computing lengths" << endl;
     vec<int> lens( D.E( ), 0 );
     #pragma omp parallel for
     for ( int e = 0; e < D.E( ); e++ )
     {    if ( D.O(e)[0] < 0 ) continue;
          for ( int j = 0; j < D.O(e).isize( ); j++ )
               lens[e] += hb.Kmers( D.O(e)[j] );    }

     // Get distances to end.

     vec<int> dfw;
     DistancesToEndArr( D, lens, MAX_KILL * MIN_RATIO, True, dfw );

     // Removing compound hanging ends.

     vec<int> dels;
     FindCompoundHangs( D, dinv, lens, dfw, dels, MAX_TINY, MIN_RATIO, verbose );
     vec<int> to_left, to_right;
     D.ToLeft(to_left), D.ToRight(to_right);

     // Look for hanging ends.

     if (verbose) cout << Date( ) << ": finding hanging ends" << endl;
     #pragma omp parallel for
     for ( int v = 0; v < D.N( ); v++ )
     {    for ( int j1 = 0; j1 < D.From(v).isize( ); j1++ )
          {    int w = D.From(v)[j1];
               if ( D.From(w).nonempty( ) || !D.To(w).solo( ) ) continue;
               int e1 = D.IFrom( v, j1 );
               const vec<int>& x1 = D.OFrom( v, j1 );
               if ( x1[0] < 0 ) continue;
               int n1 = lens[e1];
               if ( n1 > MAX_KILL ) continue;
               for ( int j2 = 0; j2 < D.From(v).isize( ); j2++ )
               {    int e2 = D.IFrom( v, j2 );
                    const vec<int>& x2 = D.OFrom( v, j2 );
                    int n2 = lens[ D.IFrom(v,j2) ] + dfw[ D.From(v)[j2] ];
                    if ( n2 < MIN_RATIO * n1 ) continue;
                    Bool contained = ( x2.Contains( x1, 0 ) );
                    int re1 = dinv[e1], support = 0;
                    support += dpaths_index.Count(e1);
                    support += dpaths_index.Count(re1);
                    Bool weak = ( support <= 1 );
                    int f = x1.front( );
                    Bool extended = False;
                    for ( int l = 0; l < nd[f].isize( ); l++ )
                    {    int e = nd[f][l].first, ep = nd[f][l].second;
                         const vec<int>& y = D.O(e);
                         int ext = dfw[ to_right[e] ];
                         for ( int j = ep + x1.isize( ); j < y.isize( ); j++ )
                              ext += hb.Kmers( y[j] );
                         if ( ext >= MIN_EXT && y.Contains( x1, ep ) )
                         {    extended = True;
                              break;    }    }
                    if ( !contained && !extended && !weak ) continue;
                    #pragma omp critical
                    {    dels.push_back( e1, dinv[e1] );    }    }    }    }
     UniqueSort(dels);
     if (verbose)
          cout << Date( ) << ": to delete " << dels.size( ) << " edges" << endl;

     // Look for weak branches.  We look for a two-way branch in which one is
     // a dead end having at most two reads supporting it, and the other has
     // at least ten times as high support.

     vec<Bool> del1( D.E( ), False );
     for ( int i = 0; i < dels.isize( ); i++ )
          del1[ dels[i] ] = True;
     int nweak = 0;
     #pragma omp parallel for
     for ( int v = 0; v < D.N( ); v++ )
     {    if ( D.From(v).size( ) != 2 ) continue;
          for ( int mm = 0; mm < 2; mm++ )
          {    int j1 = mm;
               int j2 = 1 - j1;
               int d1 = D.IFrom(v,j1), d2 = D.IFrom(v,j2);
               if ( del1[d1] || del1[d2] ) continue;
               int w1 = D.From(v)[j1], w2 = D.From(v)[j2], n1 = 0, n2 = 0;
               if ( D.From(w1).nonempty( ) || !D.To(w1).solo( ) ) continue;
               for ( int p = 1; p <= 2; p++ )
               {    int f = ( p == 1 ? d1 : dinv[d1] );
                    for ( int j = 0; j < dpaths_index.Count(f); j++ )
                    {    int64_t id = dpaths_index.Val( f, j );
                         n1++;    }    }
               for ( int p = 1; p <= 2; p++ )
               {    int f = ( p == 1 ? d2 : dinv[d2] );
                    for ( int j = 0; j < dpaths_index.Count(f); j++ )
                    {    int64_t id = dpaths_index.Val( f, j );
                         n2++;    }    }
               if ( n1 > 2 ) continue;
               if ( n2 < 10 * n1 ) continue;
               #pragma omp critical
               {    dels.push_back( d1, dinv[d1] );
                    nweak++;    }    }    }
     if (verbose)
          cout << Date( ) << ": " << nweak << " weak branches detected" << endl;

     // Look for weak branches, part two.  We look for a two-way branch in which
     // one branch has no support, and the other has support at least ten times
     // as high.

     vec<Bool> del2( D.E( ), False );
     for ( int i = 0; i < dels.isize( ); i++ )
          del2[ dels[i] ] = True;
     int nweak2 = 0;
     #pragma omp parallel for
     for ( int v = 0; v < D.N( ); v++ )
     {    if ( D.From(v).size( ) != 2 ) continue;
          for ( int mm = 0; mm < 2; mm++ )
          {    int j1 = mm;
               int j2 = 1 - j1;
               int d1 = D.IFrom(v,j1), d2 = D.IFrom(v,j2);
               if ( del2[d1] || del2[d2] ) continue;
               int n1 = 0, n2 = 0;
               for ( int p = 1; p <= 2; p++ )
               {    int f = ( p == 1 ? d1 : dinv[d1] );
                    for ( int j = 0; j < dpaths_index.Count(f); j++ )
                    {    int64_t id = dpaths_index.Val( f, j );
                         n1++;    }    }
               for ( int p = 1; p <= 2; p++ )
               {    int f = ( p == 1 ? d2 : dinv[d2] );
                    for ( int j = 0; j < dpaths_index.Count(f); j++ )
                    {    int64_t id = dpaths_index.Val( f, j );
                         n2++;    }    }
               if ( n1 > 0 || n2 < 10 ) continue;
               #pragma omp critical
               {    dels.push_back( d1, dinv[d1] );
                    del2[d1] = del2[ dinv[d1] ] = True;
                    nweak2++;    }    }    }
     if (verbose)
          cout << Date( ) << ": " << nweak2 << " weak branches detected (2)" << endl;
     D.DeleteEdges(dels);

     // Look for weak branches, part three.  We look for a two-way branch in which
     // one branch has at most one read support "from the left", and the other has
     // at least ten reads supporting it.  This caused a roughly 2% drop in N50
     // perfect stretch.

     const int MAX_WEAK3_LOSE = 1;
     const int MIN_WEAK3_WIN = 10;
     const int MIN_WEAK3_RATIO = 10;
     DelWeak3( D, dinv, dpaths, dpaths_index, dels,
          MAX_WEAK3_LOSE, MIN_WEAK3_WIN, MIN_WEAK3_RATIO, verbose );

     // Completely clean up graph and rebuild, preparatory to tidying.

     if (verbose) cout << Date( ) << ": interim cleaning" << endl;
     D.DeleteEdges(dels);
     RemoveUnneededVertices( D, dinv );
     CleanupCore( D, dinv );
     dels.clear( );
     lens.clear( );
     lens.resize( D.E( ), 0 );
     #pragma omp parallel for
     for ( int e = 0; e < D.E( ); e++ )
     {    if ( D.O(e)[0] < 0 ) continue;
          for ( int j = 0; j < D.O(e).isize( ); j++ )
               lens[e] += hb.Kmers( D.O(e)[j] );    }
     D.ToLeft(to_left), D.ToRight(to_right);
     PlaceReads( hb, paths, dup, D, dpaths, verbose, single );

     // Tidy graph by removing duff at the end of longish lines.

     RemoveDuff( hb, D, dinv, dels, verbose );

     // Delete short standalones.

     if (verbose) cout << Date( ) << ": removing short standalones" << endl;
     for ( int d = 0; d < D.E( ); d++ )
     {    int v = to_left[d], w = to_right[d];
          if ( D.To(v).nonempty( ) || !D.From(v).solo( ) ) continue;
          if ( D.From(w).nonempty( ) || !D.To(w).solo( ) ) continue;
          if ( v == w ) continue;
          if ( lens[d] >= MIN_STANDALONE ) continue;
          dels.push_back(d);    }

     // Kill components having low unique content.

     KillLowUnique( hb, D, dels, verbose );
     D.DeleteEdges(dels);

     // Look for canonical pull aparts.  There are two cases: either one edge
     // in the middle, or just a vertex.   Do three passes.

     for ( int pass = 1; pass <= 3; pass++ )
     {    PullApart( inv, D, dinv, dpaths, dels, verbose );    }

     // Clean the graph.

     if (verbose) cout << Date( ) << ": cleaning graph" << endl;
     D.DeleteEdges(dels);
     RemoveUnneededVertices( D, dinv );
     if (verbose) cout << Date( ) << ": more cleaning" << endl;
     CleanupCore( D, dinv );
     Emanate( D, dinv, verbose );
     Emanate2( inv, D, dinv, verbose );
     if (verbose)
          cout << Date( ) << ": final graph has " << D.E( ) << " edges" << endl;    }


void FlattenSomeBubbles( const HyperBasevectorX& hb, const vec<Bool>& dup,
     const ReadPathVecX& paths, digraphE<vec<int>>& D, vec<int>& dinv, 
     ReadPathVec& dpaths, const int MAX_DELTA, const double MIN_RATIO,
     const int MAX_DEL )
{
     cout << Date( ) << ": flattening bubbles" << endl;
     vec<int> dels;
     PlaceReads( hb, paths, dup, D, dpaths, True, False );
     IntIndex dpaths_index( dpaths, D.E( ), False );
     vec<int> lens( D.E( ), 0 );
     #pragma omp parallel for
     for ( int e = 0; e < D.E( ); e++ )
     {    if ( D.O(e)[0] < 0 ) continue;
          for ( int j = 0; j < D.O(e).isize( ); j++ )
               lens[e] += hb.Kmers( D.O(e)[j] );    }
     #pragma omp parallel for
     for ( int v = 0; v < D.N( ); v++ )
     {    if ( D.To(v).size( ) != 1 || D.From(v).size( ) != 2 ) continue;
          if ( D.From(v)[0] != D.From(v)[1] ) continue;
          int w = D.From(v)[0];
          if ( v == w ) continue;
          if ( D.To(w).size( ) != 2 || D.From(w).size( ) != 1 ) continue;
          if ( D.O( D.IFrom(v,0) )[0] < 0 || D.O( D.IFrom(v,1) )[0] < 0 ) continue;
          int delta = AbsDiff( lens[ D.IFrom(v,0)], lens[ D.IFrom(v,1)] );
          if ( delta == 0 || delta > MAX_DELTA ) continue;
          vec<vec<int64_t>> supp(2);
          for ( int j = 0; j < 2; j++ )
          {    int d = D.IFrom(v,j);
               int rd = dinv[d];
               for ( int i = 0; i < dpaths_index.Count(d); i++ )
                    supp[j].push_back( dpaths_index.Val(d,i) / 2 );
               for ( int i = 0; i < dpaths_index.Count(rd); i++ )
                    supp[j].push_back( dpaths_index.Val(rd,i) / 2 );
               UniqueSort( supp[j] );    }
          int m = MeetSize( supp[0], supp[1] );
          for ( int j = 0; j < 2; j++ )
          {    int d1 = D.IFrom(v,j);
               int n1 = supp[j].isize( ) - m, n2 = supp[1-j].isize( ) - m;
               if ( n1 > MAX_DEL || n2 < MIN_RATIO * Max(1, n1) ) continue;
               #pragma omp critical
               {    dels.push_back( d1, dinv[d1] );    }    }    }
     #pragma omp parallel for
     for ( int v = 0; v < D.N( ); v++ )
     {    int n = D.From(v).size( );
          if ( D.To(v).size( ) != 1 || n < 2 ) continue;
          int w = D.From(v)[0];
          Bool bad = False;
          for ( int j = 1; j < n; j++ ) if ( D.From(v)[j] != w ) bad = True;
          if ( bad || v == w ) continue;
          if ( D.To(w).isize( ) != n || D.From(w).size( ) != 1 ) continue;
          for ( int j = 0; j < n; j++ )
               if ( D.O( D.IFrom(v,j) )[0] < 0 ) bad = True;
          if (bad) continue;
          vec<vec<int64_t>> supp(n);
          for ( int j = 0; j < n; j++ )
          {    int d = D.IFrom(v,j);
               int rd = dinv[d];
               for ( int i = 0; i < dpaths_index.Count(d); i++ )
                    supp[j].push_back( dpaths_index.Val(d,i) / 2 );
               for ( int i = 0; i < dpaths_index.Count(rd); i++ )
                    supp[j].push_back( dpaths_index.Val(rd,i) / 2 );
               UniqueSort( supp[j] );    }
          int M = 0;
          for ( int j = 0; j < n; j++ ) M = Max( M, supp[j].isize( ) );
          const double MIN_RATIO2 = 10.0;
          const int MAX_KILL = 2;
          for ( int j = 0; j < n; j++ )
          {    int d = D.IFrom(v,j), m = supp[j].size( );
               if ( m <= MAX_KILL && M >= MIN_RATIO2 * Max(m,1) )
               {
                    #pragma omp critical
                    {    dels.push_back( d, dinv[d] );    }    }    }    }
     UniqueSort(dels);
     cout << Date( ) << ": deleting " << ToStringAddCommas( dels.size( ) )
          << " weak bubble edges" << endl;
     D.DeleteEdges(dels);    }

void ComputeMult( const HyperBasevectorX& hb, const digraphE<vec<int>>& D,
     vec<int>& mult )
{    mult.clear( );
     mult.resize( hb.E( ), 0 );
     for ( int d = 0; d < D.E( ); d++ )
     {    if ( D.O(d)[0] < 0 ) continue;
          for ( int i = 0; i < D.O(d).isize( ); i++ )
               mult[ D.O(d)[i] ]++;    }    }

void KillMisassembledCellsAlt( const HyperBasevectorX& hb, digraphE<vec<int>>& D,
     vec<int>& dinv, const VecIntVec& ebcx, vec<vec<vec<vec<int>>>>& dlines,
     vec<int>& dels )
{
     cout << Date( ) << ": begin fixing misassembled cells-alt" << endl;

     // Heuristics.

     const double MIN_SHARE_FRAC = 0.25;
     const int BC_MAX_CELL = 1000;
     const double SURPRISE = 4.0;

     // Get superedge lengths.

     vec<int> dlens( D.E( ), 0 );
     #pragma omp parallel for
     for ( int e = 0; e < D.E( ); e++ )
     {    if ( D.O(e)[0] < 0 ) continue;
          for ( int j = 0; j < D.O(e).isize( ); j++ )
               dlens[e] += hb.Kmers( D.O(e)[j] );    }

     // Define line buckets.

     vec<int> llens;
     GetLineLengths( hb, D, dlines, llens );
     vec<vec<int>> buckets;
     BucketLines( dlines, llens, buckets );

     // Loop across cells.

     cout << Date( ) << ": kill cells at misassemblies-alt, main loop" << endl;
     vec<int> to_right, mult;
     D.ToRight(to_right);
     ComputeMult( hb, D, mult );
     int dcells = 0;
     #pragma omp parallel for schedule(dynamic, 1)
     for ( int bi = 0; bi < buckets.isize( ); bi++ )
     {    vec< pair<int,int> > orc;
          vec<int> lens, bleft, bright, es;
          for ( int ib = 0; ib < buckets[bi].isize( ); ib++ )
          {    int dl = buckets[bi][ib];
               const vec<vec<vec<int>>>& L = dlines[dl];
               for ( int m = 1; m < L.isize( ); m += 2 )
               {    const vec<vec<int>>& M = L[m];

                    // Get cell length.

                    int ncell = 0;
                    lens.clear( );
                    for ( int i = 0; i < M.isize( ); i++ )
                    {    const vec<int>& x = M[i];
                         int len = 0;
                         for ( int j = 0; j < x.isize( ); j++ )
                         {    int d = x[j];
                              if ( D.O(d)[0] >= 0 ) len += dlens[d];    }
                         lens.push_back(len);    }
                    Sort(lens);
                    if ( lens.nonempty( ) ) ncell = Median(lens);
                    if ( ncell > BC_MAX_CELL ) continue;

                    // Compute barcodes on edges on both sides.

                    vec<int> b1, b2;
                    int d1 = L[m-1][0][0], d2 = L[m+1][0][0];
                    if ( D.O(d1)[0] < 0 || D.O(d2)[0] < 0 ) continue;
                    for ( int i = 0; i < D.O(d1).isize( ); i++ )
                    {    int e = D.O(d1)[i];
                         if ( mult[e] != 1 ) continue;
                         for ( int j = 0; j < (int) ebcx[e].size( ); j++ )
                              b1.push_back( ebcx[e][j] );    }
                    for ( int i = 0; i < D.O(d2).isize( ); i++ )
                    {    int e = D.O(d2)[i];
                         if ( mult[e] != 1 ) continue;
                         for ( int j = 0; j < (int) ebcx[e].size( ); j++ )
                              b2.push_back( ebcx[e][j] );    }
                    UniqueSort(b1), UniqueSort(b2);

                    // Check for poor overlap.

                    int n = Min( b1.size( ), b2.size( ) ), k = MeetSize( b1, b2 );
                    if ( n < 10 ) continue;
                    if ( ( k + SURPRISE * sqrt(k) ) / double(n) >= MIN_SHARE_FRAC )
                         continue;

                    // Kill the cell.

                    int g = -1;
                    if ( M.solo( ) && M[0].empty( ) && m > 0 && L[m-1].nonempty( )
                         && L[m-1][0].nonempty( ) )
                    {    int d = L[m-1][0][0];
                         int v = to_right[d];
                         if ( v >= 0 ) g = D.IFrom(v,0);    }
                    vec<int> es;
                    if ( g >= 0 ) es.push_back(g);
                    for ( int i = 0; i < M.isize( ); i++ )
                         es.append( M[i] );
                    UniqueSort(es);
                    #pragma omp critical
                    {    dels.append(es);
                         dcells++;    }    }    }    }
     cout << Date( ) << ": killed " << dcells << " cells" << endl;    }

// Splay vertices at the ends of long lines.

void Splay( digraphE<vec<int>>& D, const vec<int>& dinv,
     const vec<vec<vec<vec<int>>>>& dlines, const vec<int>& llens,
     const int MIN_SPLAY )
{
     cout << Date( ) << ": splaying vertices" << endl;
     vec<int> splays, to_left, to_right;
     D.ToLeft(to_left), D.ToRight(to_right);
     vec<int> linv;
     LineInv( dlines, dinv, linv );
     #pragma omp parallel for
     for ( int i = 0; i < dlines.isize( ); i++ )
     {    if ( llens[i] < MIN_SPLAY ) continue;
          int ip = linv[i];
          if ( ip <= i ) continue;
          for ( int pass = 1; pass <= 2; pass++ )
          {    int li = ( pass == 1 ? i : ip );

               const vec<vec<vec<int>>>& L = dlines[li];
               int v = to_left[ L.front( )[0][0] ], w = to_right[ L.back( )[0][0] ];
               if ( v < 0 || w < 0 ) continue;
               if ( D.From(v).size( ) + D.To(v).size( ) > 1 )
               {
                    #pragma omp critical
                    {    splays.push_back(v);    }    }
               if ( D.From(w).size( ) + D.To(w).size( ) > 1 )
               {
                    #pragma omp critical
                    {    splays.push_back(w);    }    }    }    }
     UniqueSort(splays);
     for ( int i = 0; i < splays.isize( ); i++ )
          D.SplayVertex( splays[i] );    }

void GetLineLengths( const HyperBasevectorX& hb, const digraphE<vec<int>>& D,
     const vec<vec<vec<vec<int>>>>& dlines, vec<int>& llens )
{
     // Compute line lengths.

     llens.clear( );
     llens.resize( dlines.size( ), 0 );
     #pragma omp parallel for
     for ( int i = 0; i < dlines.isize( ); i++ )
     {    int pos = 0;
          const vec<vec<vec<int>>>& L = dlines[i];
          for ( int j = 0; j < L.isize( ); j++ )
          {    vec<int> lensj;
               for ( int k = 0; k < L[j].isize( ); k++ )
               {    int len = 0;
                    for ( int l = 0; l < L[j][k].isize( ); l++ )
                    {    int d = L[j][k][l];
                         if ( D.O(d)[0] < 0 ) continue;
                         for ( int m = 0; m < D.O(d).isize( ); m++ )
                              len += hb.Kmers( D.O(d)[m] );    }
                    lensj.push_back(len);    }
               Sort(lensj);
               if ( lensj.nonempty( ) ) pos += Median(lensj);    }
          llens[i] = pos;    }    }

void GetLineLengths( const vec<int>& kmers, const digraphE<vec<int>>& D,
     const vec<vec<vec<vec<int>>>>& dlines, vec<int>& llens )
{
     // Compute line lengths.

     llens.clear( );
     llens.resize( dlines.size( ), 0 );
     #pragma omp parallel for schedule(dynamic, 1000)
     for ( int i = 0; i < dlines.isize( ); i++ )
     {    int pos = 0;
          const vec<vec<vec<int>>>& L = dlines[i];
          for ( int j = 0; j < L.isize( ); j++ )
          {    vec<int> lensj;
               for ( int k = 0; k < L[j].isize( ); k++ )
               {    int len = 0;
                    for ( int l = 0; l < L[j][k].isize( ); l++ )
                    {    int d = L[j][k][l];
                         if ( D.O(d)[0] < 0 ) continue;
                         for ( int m = 0; m < D.O(d).isize( ); m++ )
                              len += kmers[ D.O(d)[m] ];    }
                    lensj.push_back(len);    }
               Sort(lensj);
               if ( lensj.nonempty( ) ) pos += Median(lensj);    }
          llens[i] = pos;    }    }

void GetLineLengths( const vec<int>& dlens,
     const vec<vec<vec<vec<int>>>>& dlines, vec<int>& llens )
{    llens.clear( );
     llens.resize( dlines.size( ), 0 );
     #pragma omp parallel for schedule(dynamic, 1000)
     for ( int i = 0; i < dlines.isize( ); i++ )
     {    int pos = 0;
          const vec<vec<vec<int>>>& L = dlines[i];
          for ( int j = 0; j < L.isize( ); j++ )
          {    vec<int> lensj;
               for ( int k = 0; k < L[j].isize( ); k++ )
               {    int len = 0;
                    for ( int l = 0; l < L[j][k].isize( ); l++ )
                         len += dlens[ L[j][k][l] ];
                    lensj.push_back(len);    }
               Sort(lensj);
               if ( lensj.nonempty( ) ) pos += Median(lensj);    }
          llens[i] = pos;    }    }

void FindCompoundHangs( const digraphE<vec<int>>& D, const vec<int>& dinv,
     const vec<int>& lens, const vec<int>& dfw, vec<int>& dels,
     const int MAX_TINY, const double MIN_RATIO, const Bool verbose )
{
     if (verbose) cout << Date( ) << ": removing compound hanging ends" << endl;
     vec<int> to_left, to_right;
     D.ToLeft(to_left), D.ToRight(to_right);
     #pragma omp parallel for
     for ( int v = 0; v < D.N( ); v++ )
     {    for ( int j1 = 0; j1 < D.From(v).isize( ); j1++ )
          {    int f1 = D.IFrom(v,j1);
               if ( D.O(f1)[0] < 0 ) continue;
               int n1 = lens[f1] + dfw[ D.From(v)[j1] ];
               if ( n1 > MAX_TINY ) continue;
               for ( int j2 = 0; j2 < D.From(v).isize( ); j2++ )
               {    int f2 = D.IFrom(v,j2);
                    if ( D.O(f2)[0] < 0 ) continue;
                    int n2 = lens[f2] + dfw[ D.From(v)[j2] ];
                    if ( n2 < MIN_RATIO * n1 ) continue;
                    vec<int> es = {f1};
                    for ( int i = 0; i < es.isize( ); i++ )
                    {    int w = to_right[ es[i] ];
                         for ( int j = 0; j < D.From(w).isize( ); j++ )
                         {    int e = D.IFrom(w,j);
                              if ( Member( es, e ) ) continue; // prob. can't happen
                              es.push_back(e);    }    }
                    #pragma omp critical
                    {    for ( int i = 0; i < es.isize( ); i++ )
                              dels.push_back( es[i], dinv[ es[i] ] );    }
                                   }    }    }    }

void CleanupCore( digraphE<vec<int>>& D, vec<int>& dinv )
{    if ( D.E( ) != dinv.isize( ) )
     {    cout << "\nInternal error, involution size not equal to number of edges.\n"
               << endl;
          PRINT2( D.E( ), dinv.size( ) );
          TracebackThisProcess( );
          Scram(1);    }
     vec<Bool> used;
     D.Used(used);
     for ( int d = 0; d < D.E( ); d++ )
     {    if ( dinv[d] < 0 || dinv[d] >= D.E( ) )
          {    cout << "\nInternal error, involution value doesn't make sense.\n"
                    << endl;
               PRINT2( d, dinv[d ]);
               TracebackThisProcess( );
               Scram(1);    }
          if ( dinv[dinv[d]] != d )
          {    cout << "\nInternal error, involution is not an involution.\n"
                    << endl;
               PRINT3( d, dinv[d], dinv[dinv[d]] );
               TracebackThisProcess( );
               Scram(1);    }
          if ( used[d] != used[ dinv[d] ] )
          {    cout << "\nInternal error, involution maps used edge to unused "
                    << "edge.\n" << endl;
               PRINT(d);
               TracebackThisProcess( );
               Scram(1);    }    }
     vec<int> to_new_id( used.size( ), -1 );
     {    int count = 0;
          for ( int i = 0; i < used.isize( ); i++ )
               if ( used[i] ) to_new_id[i] = count++;    }
     vec<int> dinv2;
     for ( int i = 0; i < D.EdgeObjectCount( ); i++ )
     {    if ( !used[i] ) continue;
          if ( dinv[i] < 0 ) dinv2.push_back(-1); // SHOULD NOT HAPPEN!
          else dinv2.push_back( to_new_id[ dinv[i] ] );    }
     dinv = dinv2;
     D.RemoveDeadEdgeObjects( );
     D.RemoveEdgelessVertices( );    }

void RemoveUnneededVertices( digraphE<vec<int>>& D, vec<int>& dinv )
{
    // new algorithm
    // 1. make a list of vertices to kill
    // 2. find the "boundary" edges at the front and tail of the run to remove
    // 3. walk edges between/including boundary edges, building up new edge object,
    //    and recording edge mappings.
    // 4. add new edge object between boundary vertices;
    // 5. delete all edges marked for deletion
    //
    // this all assumes that the involution can't share edges in any
    // string of edges that we're interested in.  The only way I can see
    // this happening is a single, palindromic edge, but that would not
    // qualify for

    vec<bool> vertex_kill( D.N(), false );
    vec<int> vertex_queue, to_left, to_right;
    D.ToLeft(to_left);
    D.ToRight(to_right);

    // step 1: make a list of vertices to kill
    // o----o----o----o
    // v0   v1   v2   v3

    for ( int v = 0; v < D.N(); ++v )
        if ( D.FromSize(v) == 1 && D.ToSize(v) == 1
                && D.From(v)[0] != D.To(v)[0]
                && D.O( D.IFrom( v, 0 ) )[0] >= 0
                && D.O( D.ITo( v, 0 ) )[0] >= 0 ) {
            vertex_kill[v] = true;
            vertex_queue.push_back(v);
        }

    // step 2

    vec<pair<int,int>> bound;
    while ( vertex_queue.size() ) {
        int v = vertex_queue.back();
        vertex_queue.pop_back();
        if ( !vertex_kill[v] ) continue;     // already done it
        int eleft;
        int vleft = v;
        size_t runsize = 0;
        do {
            runsize++;
            vertex_kill[vleft] = false;
            eleft = D.EdgeObjectIndexByIndexTo(vleft,0);
            vleft = D.To(vleft)[0];
        } while ( vertex_kill[vleft] );
        int eright;
        int vright = v;
        do {
            runsize++;
            vertex_kill[vright] = false;
            eright = D.EdgeObjectIndexByIndexFrom(vright,0);
            vright = D.From(vright)[0];
        } while ( vertex_kill[vright] );
        runsize--;      // 'v' gets counted twice

        // We rely on the fact that the involution is not
        // tied up with the path here.  We decide to push on the involution
        // here, too, so that we *know* what the inv[] of the new edge is.  However,
        // this requires that we canonicalize, so we don't do this twice.  This
        // canonicalization looks odd, but is correct (I think).

        if ( eleft < dinv[eright] ) {
            // WARNING: code below relies on the fact that we're pushing on a
            // run and its involution adjacent in this list.
            bound.push(eleft,eright);
            bound.push(dinv[eright], dinv[eleft]);

        }
    }

    // steps 3 and 4

    vec<int> edge_renumber0( D.EdgeObjectCount(), vec<int>::IDENTITY );
    vec<int> new_edge_numbers;
    vec<int> to_delete;
    while ( bound.size() ) {
        auto bounds = bound.back();
        bound.pop_back();
        int new_edge_no = D.EdgeObjectCount();
        vec<int> new_edge( D.EdgeObject( bounds.first ) );
        edge_renumber0[ bounds.first ] = new_edge_no;
        to_delete.push_back( bounds.first );

        for ( int v = to_right[bounds.first]; v != to_right[bounds.second]; v = D.From(v)[0] ) {
            // for each edge...
            int edge = D.EdgeObjectIndexByIndexFrom(v,0);
            to_delete.push_back(edge);
            new_edge.SetCat( new_edge, D.EdgeObject(edge) );
            edge_renumber0[edge] = new_edge_no;
        }
        D.AddEdge(to_left[bounds.first], to_right[bounds.second], new_edge);
        new_edge_numbers.push_back(new_edge_no);
    }
    D.DeleteEdges(to_delete);

    // for each pair of newly created edges, update mInv

    dinv.resize( D.E( ) );
    for (auto itr = new_edge_numbers.begin();
            itr != new_edge_numbers.end(); advance(itr,2)) {
         dinv[itr[0]] = itr[1];
         dinv[itr[1]] = itr[0];
    }    }

// Kill components having low unique content.

void KillLowUnique( const HyperBasevectorX& hb, digraphE<vec<int>>& D,
     vec<int>& dels, const Bool verbose )
{
     const int MIN_UNIQ = 75;
     vec<int> mult;
     ComputeMult( hb, D, mult );
     if (verbose) cout << Date( ) << ": finding components" << endl;
     vec<vec<int>> comp;
     D.ComponentsEFast(comp);
     if (verbose) cout << Date( ) << ": finding unique content" << endl;
     for ( int i = 0; i < comp.isize( ); i++ )
     {    int uc = 0;
          for ( int j = 0; j < comp[i].isize( ); j++ )
          {    int d = comp[i][j];
               if ( D.O(d)[0] < 0 ) continue;
               for ( int l = 0; l < D.O(d).isize( ); l++ )
               {    int e = D.O(d)[l];
                    if ( mult[e] == 1 ) uc += hb.Kmers(e);    }    }
          if ( uc < MIN_UNIQ )
          {
               #pragma omp critical
               {    dels.append( comp[i] );    }    }    }    }

void KillLowUniqueFrac( const HyperBasevectorX& hb, digraphE<vec<int>>& D,
     vec<int>& dels, const double MIN_UNIQ_FRAC, const Bool verbose )
{
     vec<int> mult;
     ComputeMult( hb, D, mult );
     if (verbose) cout << Date( ) << ": finding components" << endl;
     vec<vec<int>> comp;
     D.ComponentsEFast(comp);
     if (verbose) cout << Date( ) << ": finding unique content" << endl;
     for ( int i = 0; i < comp.isize( ); i++ )
     {    int64_t uc = 0, total = 0;
          for ( int j = 0; j < comp[i].isize( ); j++ )
          {    int d = comp[i][j];
               if ( D.O(d)[0] < 0 ) continue;
               for ( int l = 0; l < D.O(d).isize( ); l++ )
               {    int e = D.O(d)[l];
                    total += hb.Kmers(e);
                    if ( mult[e] == 1 ) uc += hb.Kmers(e);    }    }
          if ( double(uc)/double(total) < MIN_UNIQ_FRAC )
          {
               #pragma omp critical
               {    dels.append( comp[i] );    }    }    }    }

// Tidy graph by removing duff at the end of longish lines.

void RemoveDuff( const HyperBasevectorX& hb, digraphE<vec<int>>& D,
     const vec<int>& dinv, vec<int>& dels, const Bool verbose )
{
     // Heuristics.

     const int MIN_PRE_END = 500;
     const int MAX_END = 250;
     const int MAX_KILL = 320;
     const double MIN_RATIO = 25;

     // Do it.

     vec<int> lens( D.E( ), 0 );
     #pragma omp parallel for
     for ( int e = 0; e < D.E( ); e++ )
     {    if ( D.O(e)[0] < 0 ) continue;
          for ( int j = 0; j < D.O(e).isize( ); j++ )
               lens[e] += hb.Kmers( D.O(e)[j] );    }
     vec<int> dfw;
     DistancesToEndArr( D, lens, MAX_KILL * MIN_RATIO, True, dfw );
     if (verbose) cout << Date( ) << ": finding lines" << endl;
     vec<vec<vec<vec<int>>>> dlines;
     FindLines( D, dinv, dlines, MAX_CELL_PATHS, MAX_CELL_DEPTH );

     // Get line lengths.

     if (verbose) cout << Date( ) << ": getting line lengths" << endl;
     vec<int> llens;
     GetLineLengths( hb, D, dlines, llens );

     // Remove duff.

     if (verbose) cout << Date( ) << ": removing duff" << endl;
     vec<int> to_left, to_right;
     D.ToLeft(to_left), D.ToRight(to_right);
     #pragma omp parallel for schedule( dynamic, 1000 )
     for ( int li = 0; li < dlines.isize( ); li++ )
     {    if ( llens[li] < MIN_PRE_END ) continue;
          int e = dlines[li].front( )[0][0];
          int f = dlines[li].back( )[0][0];

          // Now we've identified a line that goes from e to f.  Trim duff on the
          // right, if there's not too much.

          if ( dfw[to_right[f]] > MAX_END  ) continue;
          vec<int> vs;
          vs.push_back( to_right[f] );
          Bool have_gap = False;
          for ( int i = 0; i < vs.isize( ); i++ )
          {    int v = vs[i];
               for ( int j = 0; j < D.From(v).isize( ); j++ )
               {    int w = D.From(v)[j];
                    if ( !Member( vs, w ) ) vs.push_back(w);
                    int g = D.IFrom(v,j);
                    if ( D.O(g)[0] < 0 ) have_gap = True;    }    }
          if (have_gap) continue;
          for ( int i = 0; i < vs.isize( ); i++ )
          {    int v = vs[i];
               for ( int j = 0; j < D.From(v).isize( ); j++ )
               {    int w = D.From(v)[j];
                    if ( !Member( vs, w ) ) vs.push_back(w);
                    int g = D.IFrom(v,j);
                    #pragma omp critical
                    {    dels.push_back( g, dinv[g] );    }    }    }    }    }

// Remove duff.  Suppose that at the end of the line is a collection of edges 
// containing in total a "small" number of kmers.  Then delete these edges.  A bit 
// dangerous in how cells (gap edge type) are handled.  Also we could just liberate 
// the duff rather than delete it.

void RemoveDuff2( const HyperBasevectorX& hb, digraphE<vec<int>>& D, vec<int>& dinv )
{
     cout << Date( ) << ": removing duff 2" << endl;
     vec<vec<vec<vec<int>>>> dlines;
     FindLines( D, dinv, dlines, MAX_CELL_PATHS, MAX_CELL_DEPTH );
     vec<int> to_left, to_right;
     D.ToLeft(to_left), D.ToRight(to_right);
     vec<int> dlens( D.E( ), 0 );
     #pragma omp parallel for
     for ( int e = 0; e < D.E( ); e++ )
     {    if ( D.O(e)[0] < 0 ) continue;
          for ( int j = 0; j < D.O(e).isize( ); j++ )
               dlens[e] += hb.Kmers( D.O(e)[j] );    }
     vec<int> llens;
     GetLineLengths( hb, D, dlines, llens );
     const int MIN_LEN = 1000;
     const int MAX_DEL = 1000;
     const int MIN_RATIO = 10;
     vec<int> dels;
     for ( int l = 0; l < dlines.isize( ); l++ )
     {    if ( llens[l] < MIN_LEN ) continue;
          int d = dlines[l].back( )[0][0];
          int v = to_right[d];
          int n = 0;
          vec<int> ds;
          Bool too_big = False;
          for ( int j = 0; j < D.To(v).isize( ); j++ )
          {    int f = D.ITo(v,j);
               if ( f != d && !Member( ds, f ) ) 
               {    ds.push_back(f);
                    n += dlens[f];
                    if ( n > MAX_DEL || MIN_RATIO * n > llens[l] )
                    {    too_big = True;
                         break;    }    }    }
          for ( int j = 0; j < D.From(v).isize( ); j++ )
          {    int f = D.IFrom(v,j);
               if ( !Member( ds, f ) )
               {    ds.push_back(f);
                    n += dlens[f];
                    if ( n > MAX_DEL || MIN_RATIO * n > llens[l] )
                    {    too_big = True;
                         break;    }    }    }
          if (too_big) continue;
          for ( int i = 0; i < ds.isize( ); i++ )
          {    for ( auto x : { to_left[ds[i]], to_right[ds[i]] } )
               {    for ( int j = 0; j < D.To(x).isize( ); j++ )
                    {    int f = D.ITo(x,j);
                         if ( f != d && !Member( ds, f ) ) 
                         {    ds.push_back(f);
                              n += dlens[f];
                              if ( n > MAX_DEL || MIN_RATIO * n > llens[l] )
                              {    too_big = True;
                                   goto next;    }    }    }
                    for ( int j = 0; j < D.From(x).isize( ); j++ )
                    {    int f = D.IFrom(x,j);
                         if ( !Member( ds, f ) )
                         {    ds.push_back(f);
                              n += dlens[f];
                              if ( n > MAX_DEL || MIN_RATIO * n > llens[l] )
                              {    too_big = True;
                                   goto next;    }    }    }    }    }
          for ( int i = 0; i < ds.isize( ); i++ )
               dels.push_back( ds[i], dinv[ ds[i] ] );
          next: continue;    }
     D.DeleteEdges(dels);
     RemoveUnneededVertices( D, dinv );
     CleanupCore( D, dinv );    }

void DelWeak3( digraphE< vec<int> >& D, vec<int>& dinv, const ReadPathVec& dpaths,
     const IntIndex& dpaths_index, vec<int>& dels,
     const int MAX_WEAK3_LOSE, const int MIN_WEAK3_WIN, const int MIN_WEAK3_RATIO,
     const Bool verbose )
{
     // Look for weak branches, part three.  We look for a two-way branch in which
     // one branch has at most one read support "from the left", and the other has
     // at least ten reads supporting it.  This caused a roughly 2% drop in N50
     // perfect stretch.

     // vec<Bool> del2( D.E( ), False );
     int nweak3 = 0;
     #pragma omp parallel for
     for ( int v = 0; v < D.N( ); v++ )
     {    if ( D.From(v).size( ) != 2 ) continue;
          for ( int mm = 0; mm < 2; mm++ )
          {    int j1 = mm;
               int j2 = 1 - j1;
               int d1 = D.IFrom(v,j1), d2 = D.IFrom(v,j2);
               if ( D.O(d1)[0] < 0 || D.O(d2)[0] < 0 ) continue;
               // if ( del2[d1] || del2[d2] ) continue;
               if ( !D.To(v).solo( ) ) continue;
               int e = D.ITo(v,0);
               int n1 = 0, n2 = 0;
               for ( int mpass = 1; mpass <= 2; mpass++ )
               for ( int pass = 1; pass <= 2; pass++ )
               {    int f;
                    if ( mpass == 1 ) f = ( pass == 1 ? d1 : dinv[d1] );
                    else f = ( pass == 1 ? d2 : dinv[d2] );
                    for ( int j = 0; j < dpaths_index.Count(f); j++ )
                    {    int64_t id = dpaths_index.Val( f, j );
                         const ReadPath& p = dpaths[id];
                         for ( int l = 0; l < (int) p.size( ) - 1; l++ )
                         {    if ( pass == 1 )
                              {    if ( p[l] == e && p[l+1] == f )
                                   {    ( mpass == 1 ? n1 : n2 )++;
                                        break;    }    }
                              else
                              {    if ( p[l] == f && p[l+1] == dinv[e] )
                                   {   ( mpass == 1 ? n1 : n2 )++;
                                        break;    }    }    }   }   }
               if ( n1 > MAX_WEAK3_LOSE ) continue;
               if ( n2 < MIN_WEAK3_WIN ) continue;
               if ( n2 < MIN_WEAK3_RATIO * n1 ) continue;
               #pragma omp critical
               {    dels.push_back( d1, dinv[d1] );
                    // del2[d1] = del2[ dinv[d1] ] = True; // unstable
                    nweak3++;    }    }    }
     if (verbose)
     {    cout << Date( ) << ": " << nweak3 << " weak branches detected (3)"
               << endl;    }
     UniqueSort(dels);    }

void DelWeak4( digraphE< vec<int> >& D, vec<int>& dinv, const ReadPathVec& dpaths,
     const IntIndex& dpaths_index, vec<int>& dels,
     const int MAX_WEAK4_LOSE, const int MAX_WEAK4_LOSE_TOTAL, 
     const int MIN_WEAK4_WIN, const int MIN_WEAK4_RATIO, const Bool verbose )
{
     // Look for weak branches, part three.  We look for a two-way branch in which
     // one branch has at most one read support "from the left", and the other has
     // at least ten reads supporting it.  This caused a roughly 2% drop in N50
     // perfect stretch.

     // vec<Bool> del2( D.E( ), False );
     int nweak3 = 0;
     #pragma omp parallel for
     for ( int v = 0; v < D.N( ); v++ )
     {    if ( D.From(v).size( ) != 2 ) continue;
          for ( int mm = 0; mm < 2; mm++ )
          {    int j1 = mm;
               int j2 = 1 - j1;
               int d1 = D.IFrom(v,j1), d2 = D.IFrom(v,j2);
               if ( D.O(d1)[0] < 0 || D.O(d2)[0] < 0 ) continue;
               // if ( del2[d1] || del2[d2] ) continue;
               if ( !D.To(v).solo( ) ) continue;
               int e = D.ITo(v,0);
               int n1 = 0, n2 = 0;
               int total1 = 0, total2 = 0;
               for ( int mpass = 1; mpass <= 2; mpass++ )
               for ( int pass = 1; pass <= 2; pass++ )
               {    int f;
                    if ( mpass == 1 ) f = ( pass == 1 ? d1 : dinv[d1] );
                    else f = ( pass == 1 ? d2 : dinv[d2] );
                    for ( int j = 0; j < dpaths_index.Count(f); j++ )
                    {    int64_t id = dpaths_index.Val( f, j );
                         ( mpass == 1 ? total1 : total2 )++;
                         const ReadPath& p = dpaths[id];
                         for ( int l = 0; l < (int) p.size( ) - 1; l++ )
                         {    if ( pass == 1 )
                              {    if ( p[l] == e && p[l+1] == f )
                                   {    ( mpass == 1 ? n1 : n2 )++;
                                        break;    }    }
                              else
                              {    if ( p[l] == f && p[l+1] == dinv[e] )
                                   {   ( mpass == 1 ? n1 : n2 )++;
                                        break;    }    }    }   }   }
               if ( total1 > MAX_WEAK4_LOSE_TOTAL ) continue;
               if ( n1 > MAX_WEAK4_LOSE ) continue;
               if ( n2 < MIN_WEAK4_WIN ) continue;
               if ( n2 < MIN_WEAK4_RATIO * n1 ) continue;
               #pragma omp critical
               {    dels.push_back( d1, dinv[d1] );
                    // del2[d1] = del2[ dinv[d1] ] = True; // unstable
                    nweak3++;    }    }    }
     if (verbose)
     {    cout << Date( ) << ": " << nweak3 << " weak branches detected (3)"
               << endl;    }
     UniqueSort(dels);    }




void DelWeak5( digraphE< vec<int> >& D, vec<int>& dinv, const ReadPathVec& dpaths,
     const IntIndex& dpaths_index, vec<int>& dels,
     const int MAX_WEAK4_LOSE, const int MAX_WEAK4_LOSE_TOTAL, 
     const int MIN_WEAK4_WIN, const int MIN_WEAK4_RATIO, const Bool verbose )
{
     int nweak3 = 0;
     #pragma omp parallel for
     for ( int v = 0; v < D.N( ); v++ )
     {    if ( D.To(v).size( ) != 1 || D.From(v).size( ) != 2 ) continue;
          if ( D.From(v)[0] != D.From(v)[1] ) continue;
          int w = D.From(v)[0];
          if ( v == w ) continue;
          if ( D.To(w).size( ) != 2 || D.From(w).size( ) != 1 ) continue;
          if ( D.O( D.IFrom(v,0) )[0] < 0 || D.O( D.IFrom(v,1) )[0] < 0 ) continue;
          for ( int mm = 0; mm < 2; mm++ )
          {    int j1 = mm;
               int j2 = 1 - j1;
               int d1 = D.IFrom(v,j1), d2 = D.IFrom(v,j2);
               if ( D.O(d1)[0] < 0 || D.O(d2)[0] < 0 ) continue;
               if ( !D.To(v).solo( ) ) continue;
               int e = D.ITo(v,0);
               int n1 = 0, n2 = 0;
               int total1 = 0, total2 = 0;
               for ( int mpass = 1; mpass <= 2; mpass++ )
               for ( int pass = 1; pass <= 2; pass++ )
               {    int f;
                    if ( mpass == 1 ) f = ( pass == 1 ? d1 : dinv[d1] );
                    else f = ( pass == 1 ? d2 : dinv[d2] );
                    for ( int j = 0; j < dpaths_index.Count(f); j++ )
                    {    int64_t id = dpaths_index.Val( f, j );
                         ( mpass == 1 ? total1 : total2 )++;
                         const ReadPath& p = dpaths[id];
                         for ( int l = 0; l < (int) p.size( ) - 1; l++ )
                         {    if ( pass == 1 )
                              {    if ( p[l] == e && p[l+1] == f )
                                   {    ( mpass == 1 ? n1 : n2 )++;
                                        break;    }    }
                              else
                              {    if ( p[l] == f && p[l+1] == dinv[e] )
                                   {   ( mpass == 1 ? n1 : n2 )++;
                                        break;    }    }    }   }   }
               if ( total1 > MAX_WEAK4_LOSE_TOTAL ) continue;
               if ( n1 > MAX_WEAK4_LOSE ) continue;
               if ( n2 < MIN_WEAK4_WIN ) continue;
               if ( n2 < MIN_WEAK4_RATIO * n1 ) continue;
               #pragma omp critical
               {    dels.push_back( d1, dinv[d1] );
                    nweak3++;    }    }    }
     if (verbose)
     {    cout << Date( ) << ": " << nweak3 << " weak branches detected (5)"
               << endl;    }
     UniqueSort(dels);    }




void Cleaner( const HyperBasevectorX& hb, const vec<int>& inv,
     const ReadPathVec& paths, const vec<Bool>& dup, digraphE< vec<int> >& D,
     vec<int>& dinv, ReadPathVec& dpaths, const Bool verbose )
{
     // Place reads on the graph.

     if (verbose) cout << Date( ) << ": placing reads on supergraph" << endl;
     const Bool single = False;
     PlaceReads( hb, paths, dup, D, dpaths, verbose, single );

     // Index the paths.

     if (verbose) cout << Date( ) << ": indexing paths" << endl;
     IntIndex dpaths_index( dpaths, D.E( ), verbose );

     // Index D.

     if (verbose) cout << Date( ) << ": creating index" << endl;
     vec<vec<pair<int,int>>> nd( hb.E( ) );
     for ( int e = 0; e < D.E( ); e++ )
     {    const vec<int>& x = D.O(e);
          if ( x[0] < 0 ) continue;
          for ( int j = 0; j < x.isize( ); j++ )
               nd[ x[j] ].push( e, j );    }

     // Heuristics.  Note that MAX_END was 50, and was changed to 250.  This is
     // kind of covering up stuff, and could be losing power.
     // Ditto for MAX_TINY, which was 10.

     const int MAX_KILL = 320;
     const double MIN_RATIO = 25;
     const int MIN_EXT = 500;
     const int MAX_TINY = 150;
     const int MIN_STANDALONE = 200;

     // Find lengths.

     if (verbose) cout << Date( ) << ": computing lengths" << endl;
     vec<int> lens( D.E( ), 0 );
     #pragma omp parallel for
     for ( int e = 0; e < D.E( ); e++ )
     {    if ( D.O(e)[0] < 0 ) continue;
          for ( int j = 0; j < D.O(e).isize( ); j++ )
               lens[e] += hb.Kmers( D.O(e)[j] );    }

     // Get distances to end.

     vec<int> dfw;
     DistancesToEndArr( D, lens, MAX_KILL * MIN_RATIO, True, dfw );

     // Removing compound hanging ends.

     vec<int> dels;
     FindCompoundHangs( D, dinv, lens, dfw, dels, MAX_TINY, MIN_RATIO, verbose );
     vec<int> to_left, to_right;
     D.ToLeft(to_left), D.ToRight(to_right);

     // Look for hanging ends.

     if (verbose) cout << Date( ) << ": finding hanging ends" << endl;
     #pragma omp parallel for
     for ( int v = 0; v < D.N( ); v++ )
     {    for ( int j1 = 0; j1 < D.From(v).isize( ); j1++ )
          {    int w = D.From(v)[j1];
               if ( D.From(w).nonempty( ) || !D.To(w).solo( ) ) continue;
               int e1 = D.IFrom( v, j1 );
               const vec<int>& x1 = D.OFrom( v, j1 );
               if ( x1[0] < 0 ) continue;
               int n1 = lens[e1];
               if ( n1 > MAX_KILL ) continue;
               for ( int j2 = 0; j2 < D.From(v).isize( ); j2++ )
               {    int e2 = D.IFrom( v, j2 );
                    const vec<int>& x2 = D.OFrom( v, j2 );
                    int n2 = lens[ D.IFrom(v,j2) ] + dfw[ D.From(v)[j2] ];
                    if ( n2 < MIN_RATIO * n1 ) continue;
                    Bool contained = ( x2.Contains( x1, 0 ) );
                    int re1 = dinv[e1], support = 0;
                    support += dpaths_index.Count(e1);
                    support += dpaths_index.Count(re1);
                    Bool weak = ( support <= 1 );
                    int f = x1.front( );
                    Bool extended = False;
                    for ( int l = 0; l < nd[f].isize( ); l++ )
                    {    int e = nd[f][l].first, ep = nd[f][l].second;
                         const vec<int>& y = D.O(e);
                         int ext = dfw[ to_right[e] ];
                         for ( int j = ep + x1.isize( ); j < y.isize( ); j++ )
                              ext += hb.Kmers( y[j] );
                         if ( ext >= MIN_EXT && y.Contains( x1, ep ) )
                         {    extended = True;
                              break;    }    }
                    if ( !contained && !extended && !weak ) continue;
                    #pragma omp critical
                    {    dels.push_back( e1, dinv[e1] );    }    }    }    }
     UniqueSort(dels);
     if (verbose)
          cout << Date( ) << ": to delete " << dels.size( ) << " edges" << endl;

     // Look for weak branches.  We look for a two-way branch in which one is
     // a dead end having at most two reads supporting it, and the other has
     // at least ten times as high support.

     vec<Bool> del1( D.E( ), False );
     for ( int i = 0; i < dels.isize( ); i++ )
          del1[ dels[i] ] = True;
     int nweak = 0;
     #pragma omp parallel for
     for ( int v = 0; v < D.N( ); v++ )
     {    if ( D.From(v).size( ) != 2 ) continue;
          for ( int mm = 0; mm < 2; mm++ )
          {    int j1 = mm;
               int j2 = 1 - j1;
               int d1 = D.IFrom(v,j1), d2 = D.IFrom(v,j2);
               if ( del1[d1] || del1[d2] ) continue;
               int w1 = D.From(v)[j1], w2 = D.From(v)[j2], n1 = 0, n2 = 0;
               if ( D.From(w1).nonempty( ) || !D.To(w1).solo( ) ) continue;
               for ( int p = 1; p <= 2; p++ )
               {    int f = ( p == 1 ? d1 : dinv[d1] );
                    for ( int j = 0; j < dpaths_index.Count(f); j++ )
                    {    int64_t id = dpaths_index.Val( f, j );
                         n1++;    }    }
               for ( int p = 1; p <= 2; p++ )
               {    int f = ( p == 1 ? d2 : dinv[d2] );
                    for ( int j = 0; j < dpaths_index.Count(f); j++ )
                    {    int64_t id = dpaths_index.Val( f, j );
                         n2++;    }    }
               if ( n1 > 2 ) continue;
               if ( n2 < 10 * n1 ) continue;
               #pragma omp critical
               {    dels.push_back( d1, dinv[d1] );
                    nweak++;    }    }    }
     if (verbose)
          cout << Date( ) << ": " << nweak << " weak branches detected" << endl;

     // Look for weak branches, part two.  We look for a two-way branch in which
     // one branch has no support, and the other has support at least ten times
     // as high.

     int nweak2 = 0;
     #pragma omp parallel for
     for ( int v = 0; v < D.N( ); v++ )
     {    if ( D.From(v).size( ) != 2 ) continue;
          for ( int mm = 0; mm < 2; mm++ )
          {    int j1 = mm;
               int j2 = 1 - j1;
               int d1 = D.IFrom(v,j1), d2 = D.IFrom(v,j2);
               int n1 = 0, n2 = 0;
               for ( int p = 1; p <= 2; p++ )
               {    int f = ( p == 1 ? d1 : dinv[d1] );
                    for ( int j = 0; j < dpaths_index.Count(f); j++ )
                    {    int64_t id = dpaths_index.Val( f, j );
                         n1++;    }    }
               for ( int p = 1; p <= 2; p++ )
               {    int f = ( p == 1 ? d2 : dinv[d2] );
                    for ( int j = 0; j < dpaths_index.Count(f); j++ )
                    {    int64_t id = dpaths_index.Val( f, j );
                         n2++;    }    }
               if ( n1 > 0 || n2 < 10 ) continue;
               #pragma omp critical
               {    dels.push_back( d1, dinv[d1] );
                    nweak2++;    }    }    }
     if (verbose)
          cout << Date( ) << ": " << nweak2 << " weak branches detected (2)" << endl;
     D.DeleteEdges(dels);

     // Look for weak branches, part three.  We look for a two-way branch in which
     // one branch has at most one read support "from the left", and the other has
     // at least ten reads supporting it.  This caused a roughly 2% drop in N50
     // perfect stretch.

     const int MAX_WEAK3_LOSE = 1;
     const int MIN_WEAK3_WIN = 10;
     const int MIN_WEAK3_RATIO = 10;
     DelWeak3( D, dinv, dpaths, dpaths_index, dels,
          MAX_WEAK3_LOSE, MIN_WEAK3_WIN, MIN_WEAK3_RATIO, verbose );

     // Completely clean up graph and rebuild, preparatory to tidying.

     if (verbose) cout << Date( ) << ": interim cleaning" << endl;
     D.DeleteEdges(dels);
     RemoveUnneededVertices( D, dinv );
     CleanupCore( D, dinv );
     dels.clear( );
     lens.clear( );
     lens.resize( D.E( ), 0 );
     #pragma omp parallel for
     for ( int e = 0; e < D.E( ); e++ )
     {    if ( D.O(e)[0] < 0 ) continue;
          for ( int j = 0; j < D.O(e).isize( ); j++ )
               lens[e] += hb.Kmers( D.O(e)[j] );    }
     D.ToLeft(to_left), D.ToRight(to_right);
     PlaceReads( hb, paths, dup, D, dpaths, verbose, single );

     // Tidy graph by removing duff at the end of longish lines.

     RemoveDuff( hb, D, dinv, dels, verbose );

     // Delete short standalones.

     if (verbose) cout << Date( ) << ": removing short standalones" << endl;
     for ( int d = 0; d < D.E( ); d++ )
     {    int v = to_left[d], w = to_right[d];
          if ( D.To(v).nonempty( ) || !D.From(v).solo( ) ) continue;
          if ( D.From(w).nonempty( ) || !D.To(w).solo( ) ) continue;
          if ( v == w ) continue;
          if ( lens[d] >= MIN_STANDALONE ) continue;
          dels.push_back(d);    }

     // Kill components having low unique content.

     KillLowUnique( hb, D, dels, verbose );
     D.DeleteEdges(dels);

     // Look for canonical pull aparts.  There are two cases: either one edge
     // in the middle, or just a vertex.   Do three passes.

     for ( int pass = 1; pass <= 3; pass++ )
     {    PullApart( inv, D, dinv, dpaths, dels, verbose );    }

     // Clean the graph.

     if (verbose) cout << Date( ) << ": cleaning graph" << endl;
     D.DeleteEdges(dels);
     RemoveUnneededVertices( D, dinv );
     if (verbose) cout << Date( ) << ": more cleaning" << endl;
     CleanupCore( D, dinv );
     Emanate( D, dinv, verbose );
     Emanate2( inv, D, dinv, verbose );
     if (verbose)
          cout << Date( ) << ": final graph has " << D.E( ) << " edges" << endl;    }

void Zipper( digraphE<vec<int>>& D, vec<int>& dinv )
{
    cout << Date( ) << ": zippering" << endl;
    // for each super-edge, the vertices to the left and to the right
    vec<int> to_left, to_right;
    D.ToLeft(to_left), D.ToRight(to_right);
    // total count of zipper operations
    int allzips = 0;
    // total number of times we are iterating over the entire graph.
    int loop_count=0;
    // list of vertices to delete after each zippering loop
    vec<int> dels;
    while(1) { //we are going to iteratively zipper
        int zips = 0;
        // list of candidate vertices for zippering.
        vec<int> can;
        #pragma omp parallel for
        for ( int v = 0; v < D.N( ); v++ )  {
        // there must be two outgoing edges at least.
        if ( D.From(v).size( ) >= 2 )   {
        // for each pair of outgoing edges...
        for ( int j1 = 0; j1 < D.From(v).isize( ); j1++ )
            for ( int j2 = j1+1; j2 < D.From(v).isize( ); j2++ )    {
                int d1 = D.IFrom(v,j1), d2 = D.IFrom(v,j2);
                // if the super edge shares a common base graph (BG) edge
                if ( D.O(d1)[0] == D.O(d2)[0] ) {
                    #pragma omp critical
                    {    can.push_back(v);    }    }    }    }    }
        // sort the vertices and remove duplicates
        ParallelUniqueSort(can);
        // go through the list of candidate vertices
        for ( int cv = 0; cv < can.isize( ); cv++ ) {
            int v = can[cv];
            // go through all pairs of outgoing edges
            for ( int j1 = 0; j1 < D.From(v).isize( ); j1++ )  {
                for ( int j2 = 0; j2 < D.From(v).isize( ); j2++ )   {
                    int d1 = D.IFrom(v,j1), d2 = D.IFrom(v,j2);
                    // canonically order the edges so that d1 has fewer BG edges than d2.
                    if ( j1 == j2 || D.O(d1).isize( ) > D.O(d2).isize( ) ) continue;
                    // don't zipper if the edge is a gap.
                    if ( D.O(d1)[0] < 0 || D.O(d2)[0] < 0)
                        continue;
                    int e = D.O(d1)[0];
                    // if this pair of edges has a common BG edge
                    if ( e == D.O(d2)[0] )  {
                        // reverse complement edges.
                        int rd1 = dinv[d1], rd2 = dinv[d2];
                        int rv = to_right[rd1];
                        // don't zipper if rd1, rd2 is degenerate with d1 and d2
                        if ( !IsUnique( d1, d2, rd1, rd2 ) )
                            continue;
                        // potential pathologies with shared vertices
                        // if d1 and d2 form a bubble...
                        if ( to_right[d1] == to_right[d2] )  {
                            if (!IsUnique(vec<int> {v, to_right[d1], rv, to_left[rd1]}))
                continue; }
            // else if they are going to distinct vertices...
            else {
                            if (!IsUnique(vec<int>  {v, to_right[d1], to_right[d2],
                                             rv, to_left[rd1], to_left[rd2] }))
                continue;
                        }
                        // BG edges that will be zippered together
                        vec<int> common, rcommon;
                        // remaining BG edges on d1 and d2 and rd1 and rd2
                        vec<int> x1, x2, rx1, rx2;
                        // these are the number of edge labels that can be zippered
                        int num_common = 0;
                        int n1 = D.O(d1).isize( );
                        for (num_common = 0; num_common < n1; num_common++) {
                            if ( D.O(d1)[num_common] != D.O(d2)[num_common] ||
                                 D.O(d1)[num_common] < 0 )
                                 break;
                            common.push_back( D.O(d1)[num_common]);
                            rcommon.push_back( D.O(rd1)[n1 - num_common - 1]);
                        }
                        // reverse rcommon
                        for ( int l = 0; l < num_common/2; l++ )    {
                            int temp = rcommon[l];
                            rcommon[l] = rcommon[num_common-l-1];
                            rcommon[num_common-l-1] = temp;
                        }
                        //cout << "common elements " << num_common << endl;
                        // now find the edge labels that are not common
                        for ( int l = num_common; l < D.O(d1).isize( ); l++ )
                            x1.push_back( D.O(d1)[l] );
                        for ( int l = num_common; l < D.O(d2).isize( ); l++ )
                            x2.push_back( D.O(d2)[l] );
                        for ( int l = 0; l < D.O(rd1).isize( ) - num_common; l++ )
                            rx1.push_back( D.O(rd1)[l] );
                        for ( int l = 0; l < D.O(rd2).isize( ) - num_common; l++ )
                            rx2.push_back( D.O(rd2)[l] );
                        zips++;
                        int N = D.N( );
                        // one for the new edge, and one for the reverse complement.
                        D.AddVertices(2);
                        int d = D.AddEdgeWithUpdate( v, N, common, to_left, to_right );
                        int rd = D.AddEdgeWithUpdate( N+1, rv, rcommon, to_left, to_right );
                        // update the involution. dinv[d] = rd and dinv[rd]=d
                        dinv.push_back( rd, d );
                        // now modify the graph.
                        if ( x1.empty( ) && x2.empty( ) ) {
                            dels.push_back( d1, d2, rd1, rd2 );
                            D.TransferEdgesWithUpdate( to_right[d1], to_right[d],
                                                       to_left, to_right);
                            D.TransferEdgesWithUpdate( to_left[rd1], to_left[rd],
                                                       to_left, to_right );
                            if ( to_right[d1] != to_right[d2] ) {
                                D.TransferEdgesWithUpdate( to_right[d2], to_right[d],
                                                           to_left, to_right);
                                D.TransferEdgesWithUpdate( to_left[rd2],
                                                           to_left[rd], to_left, to_right );
                                }
                            }
                        else if ( x1.empty( ) ) {
                            dels.push_back( d1, rd1 );
                            D.TransferEdgesWithUpdate(  to_right[d1], to_right[d],
                                                        to_left, to_right);
                            D.TransferEdgesWithUpdate(  to_left[rd1], to_left[rd],
                                                        to_left, to_right );
                            D.OMutable(d2) = x2, D.OMutable(rd2) = rx2;
                            D.GiveEdgeNewFromVx( d2, v, N );
                            to_left[d2] = N;
                            D.GiveEdgeNewToVx( rd2, rv, N+1 );
                            to_right[rd2] = N+1;
                        }
                        else  {
                        // note that x2 cannot be empty since it is at least as long as x1
                        D.OMutable(d1) = x1, D.OMutable(rd1) = rx1;
                        D.OMutable(d2) = x2, D.OMutable(rd2) = rx2;
                        D.GiveEdgeNewFromVx( d1, v, N );
                        to_left[d1] = N;
                        D.GiveEdgeNewToVx( rd1, rv, N+1 );
                        to_right[rd1] = N+1;
                        D.GiveEdgeNewFromVx( d2, v, N );
                        to_left[d2] = N;
                        D.GiveEdgeNewToVx( rd2, rv, N+1 );
                        to_right[rd2] = N+1;    }
                        // we are going to zipper once, and move to next candidate vertex.
                        // since we are going to edit the graph after one zipper,
                        // we zipper once for each vertex. Break out of the double loop using goto.
                        goto next_candidate;
                    } // close of if
                }  } // close of double for loop
        next_candidate: continue;
        } //close of candidate loop
    // if we didn't zip at all and parsed the entire graph, quit.
    if ( zips == 0 )
        break;
    allzips += zips;
    D.DeleteEdgesWithUpdate( dels, to_left, to_right );
    dels.clear( );
    loop_count++;
    } // close of while (1)
    cout << Date( ) << ": made " << allzips << " zips. ";
    cout << loop_count << " loops." << endl;
    RemoveUnneededVertices( D, dinv );
    CleanupCore( D, dinv );
}

void Validate( const HyperBasevectorX& hb, const vec<int>& inv,
     const digraphE<vec<int>>& D, const vec<int>& dinv )
{
     // Check for illegal base edge ids.

     vec<int> bads;
     #pragma omp parallel for schedule( dynamic, 10000 )
     for ( int d = 0; d < D.E( ); d++ )
     {    const vec<int>& x = D.O(d);
          if ( x.empty( ) || x[0] < 0 ) continue;
          for ( auto e : x )
          {    if ( e < 0 || e >= hb.E( ) )
               {
                    #pragma omp critical
                    {    bads.push_back(d);    }    }    }    }
     if ( bads.nonempty( ) )
     {    Sort(bads);
          int d = bads[0];
          cout << "\nIllegal superedge " << d << " = "
               << printSeq( D.O(d) ) << endl;
          TracebackThisProcess( );
          Scram(1);    }

     // Check for empty edges.

     #pragma omp parallel for schedule( dynamic, 10000 )
     for ( int e = 0; e < D.E( ); e++ )
     {    const vec<int>& x = D.O(e);
          if ( x.empty( ) )
          {
               #pragma omp critical
               {    bads.push_back(e);    }    }    }
     Sort(bads);
     for ( int bi = 0; bi < bads.isize( ); bi++ )
     {    int e = bads[bi];
          const vec<int>& x = D.O(e);
          cout << "\nIllegal empty edge." << endl;
          TracebackThisProcess( );
          Scram(1);    }

     // Check for illegal adjacencies.

     #pragma omp parallel for schedule( dynamic, 10000 )
     for ( int e = 0; e < D.E( ); e++ )
     {    const vec<int>& x = D.O(e);
          if ( x[0] < 0 ) continue;
          for ( int j = 0; j < x.isize( ) - 1; j++ )
          {    if ( !hb.Abut( x[j], x[j+1] ) )
               {
                    #pragma omp critical
                    {    bads.push_back(e);    }    }    }    }
     Sort(bads);
     for ( int bi = 0; bi < bads.isize( ); bi++ )
     {    int e = bads[bi];
          const vec<int>& x = D.O(e);
          if ( x[0] < 0 ) continue;
          for ( int j = 0; j < x.isize( ) - 1; j++ )
          {    if ( !hb.Abut( x[j], x[j+1] ) )
               {    cout << "\nNon-adjacency between " << x[j] << " and "
                         << x[j+1] << " within edge " << e << "." << endl;
                    TracebackThisProcess( );
                    Scram(1);    }    }    }

     // Second check for illegal adjacencies.

     #pragma omp parallel for schedule( dynamic, 10000 )
     for ( int v = 0; v < D.N( ); v++ )
     {    for ( int j1 = 0; j1 < D.To(v).isize( ); j1++ )
          for ( int j2 = 0; j2 < D.From(v).isize( ); j2++ )
          {    const vec<int>& x1 = D.OTo( v, j1 );
               const vec<int>& x2 = D.OFrom( v, j2 );
               if ( x1[0] < 0 || x2[0] < 0 ) continue;
               if ( !hb.Abut( x1.back( ), x2.front( ) ) )
               {
                    #pragma omp critical
                    {    bads.push_back(v);    }    }    }    }
     Sort(bads);
     for ( int bi = 0; bi < bads.isize( ); bi++ )
     {    int v = bads[bi];
          for ( int j1 = 0; j1 < D.To(v).isize( ); j1++ )
          for ( int j2 = 0; j2 < D.From(v).isize( ); j2++ )
          {    const vec<int>& x1 = D.OTo( v, j1 );
               const vec<int>& x2 = D.OFrom( v, j2 );
               if ( x1[0] < 0 || x2[0] < 0 ) continue;
               if ( !hb.Abut( x1.back( ), x2.front( ) ) )
               {    cout << "\nEncountered nonadjacency." << endl;
                    cout << "Superedge " << D.ITo( v, j1 ) << " ends with base edge "
                         << x1.back( ) << "." << endl;
                    cout << "Superedge " << D.IFrom( v, j2 ) 
                         << " begins with base edge " << x2.front( ) << "." << endl;
                    cout << "But " << x1.back( ) << "," << x2.front( )
                         << " is invalid." << endl;
                    TracebackThisProcess( );
                    Scram(1);    }    }    }

     // Check involution.

     if ( D.E( ) != dinv.isize( ) )
     {    cout << "\nInvolution has wrong size." << endl;
          TracebackThisProcess( );
          Scram(1);    }
     #pragma omp parallel for schedule( dynamic, 10000 )
     for ( int e = 0; e < D.E( ); e++ )
     {    if ( dinv[e] < 0 )
          {
               #pragma omp critical
               {    bads.push_back(e);    }    }
          else if ( dinv[dinv[e]] != e )
          {
               #pragma omp critical
               {    bads.push_back(e);    }    }    }
     Sort(bads);
     for ( int bi = 0; bi < bads.isize( ); bi++ )
     {    int e = bads[bi];
          if ( dinv[e] < 0 )
          {    cout << "\nInvolution of supergraph edge " << e
                    << " is " << dinv[e] << ", which doesn't make sense." << endl;
               TracebackThisProcess( );
               Scram(1);    }
          if ( dinv[dinv[e]] != e )
          {    cout << "\nInvolution is not an involution." << endl;
               TracebackThisProcess( );
               Scram(1);    }    }

     // Check inv versus graph structure.

     vec<int> to_left, to_right;
     D.ToLeftParallel(to_left), D.ToRightParallel(to_right);
     #pragma omp parallel for schedule( dynamic, 10000 )
     for ( int v = 0; v < D.N( ); v++ )
     for ( int j1 = 0; j1 < D.To(v).isize( ); j1++ )
     for ( int j2 = 0; j2 < D.From(v).isize( ); j2++ )
     {    int e1 = D.ITo( v, j1 ), e2 = D.IFrom( v, j2 );
          int re1 = dinv[e1], re2 = dinv[e2];
          if ( to_right[re2] != to_left[re1] )
          {
               #pragma omp critical
               {    bads.push_back(v);    }    }    }
     Sort(bads);
     for ( int bi = 0; bi < bads.isize( ); bi++ )
     {    int v = bads[bi];
          for ( int j1 = 0; j1 < D.To(v).isize( ); j1++ )
          for ( int j2 = 0; j2 < D.From(v).isize( ); j2++ )
          {    int d1 = D.ITo( v, j1 ), d2 = D.IFrom( v, j2 );
               int rd1 = dinv[d1], rd2 = dinv[d2];
               if ( to_right[rd2] != to_left[rd1] )
               {    cout << "\n";
                    PRINT4( d1, d2, rd1, rd2 );
                    cout << "superedge d1 abuts superedge d2, "
                         "so rd2 should abut rd1, but it doesn't" << endl;
                    cout << "actual edges after rd2 are:";
                    int v = to_right[rd2], w = to_left[rd1];
                    for ( int l = 0; l < D.From(v).isize( ); l++ )
                         cout << " " << D.IFrom( v, l );
                    cout << "\nactual edges before rd1 are:";
                    for ( int l = 0; l < D.To(w).isize( ); l++ )
                         cout << " " << D.ITo( w, l );
                    cout << endl;
                    cout << "Involution not consistent with graph structure."
                         << endl;
                    TracebackThisProcess( );
                    Scram(1);    }    }    }
     #pragma omp parallel for schedule( dynamic, 10000 )
     for ( int d = 0; d < D.E( ); d++ )
     {    const vec<int> &x = D.O(d), &y = D.O( dinv[d] );
          if ( x.size( ) != y.size( ) )
          {
               #pragma omp critical
               {    bads.push_back(d);    }    }
          if ( x[0] < 0 || y[0] < 0 )
          {    if ( x[0] >= 0 || y[0] >= 0 )
               {
                    #pragma omp critical
                    {    bads.push_back(d);    }    }    }
          else
          {    for ( int j = 0; j < x.isize( ); j++ )
               {    if ( y[j] != inv[ x[ x.isize( ) - j - 1 ] ] )
                    {
                         #pragma omp critical
                         {    bads.push_back(d);    }    }    }    }    }
     Sort(bads);
     for ( int bi = 0; bi < bads.isize( ); bi++ )
     {    int d = bads[bi];
          const vec<int> &x = D.O(d), &y = D.O( dinv[d] );
          if ( x.size( ) != y.size( ) )
          {    cout << "\nInvoluted edge has different size." << endl;
               TracebackThisProcess( );
               Scram(1);    }
          if ( x[0] < 0 || y[0] < 0 )
          {    if ( x[0] >= 0 || y[0] >= 0 )
               {    cout << "\nInconsistent involuted gap edges." << endl;
                    TracebackThisProcess( );
                    Scram(1);    }    }
          else
          {    for ( int j = 0; j < x.isize( ); j++ )
               {    if ( y[j] != inv[ x[ x.isize( ) - j - 1 ] ] )
                    {    cout << "\nThe two involutions are inconsistent." << endl;
                         TracebackThisProcess( );
                         Scram(1);    }    }    }    }

     // Validate gap edges.
     
     ValidateGapEdges( hb, inv, D, dinv );

     // Check for zipper down.  Turned off.  This wouldn't work per se because
     // some things can't be zippered: see Zipper.

     /*
     for ( int v = 0; v < D.N( ); v++ )
     {    vec<int> x;
          if ( D.From(v).size( ) > 1 )
          {    x.clear( );
               for ( int j = 0; j < D.From(v).isize( ); j++ )
                    x.push_back( D.O( D.IFrom(v,j) )[0] );
               Sort(x);
               for ( int i = 0; i < x.isize( ); i++ )
               {    int j = x.NextDiff(i);
                    if ( j - i > 1 )
                    {    cout << "\nEdge " << x[i] << " appears multiply "
                              << "exiting the same edge." << endl;
                         TracebackThisProcess( );
                         Scram(1);    }
                    i = j - 1;    }    }
          if ( D.To(v).size( ) > 1 )
          {    x.clear( );
               for ( int j = 0; j < D.To(v).isize( ); j++ )
                    x.push_back( D.O( D.ITo(v,j) ).back( ) );
               Sort(x);
               for ( int i = 0; i < x.isize( ); i++ )
               {    int j = x.NextDiff(i);
                    if ( j - i > 1 )
                    {    cout << "\nEdge " << x[i] << " appears multiply "
                              << "entering the same edge." << endl;
                         TracebackThisProcess( );
                         Scram(1);    }
                    i = j - 1;    }    }    }
     */

          }

void Emanate( digraphE<vec<int>>& D, vec<int>& dinv, const Bool verbose )
{    if (verbose) cout << Date( ) << ": start emanate" << endl;
     vec<int> to_right;
     D.ToRight(to_right);
     int N = D.N( );
     for ( int w = 0; w < N; w++ )
     {    for ( int j1 = 0; j1 < D.To(w).isize( ); j1++ )
          for ( int j2 = j1+1; j2 < D.To(w).isize( ); j2++ )
          {    int w1 = D.To(w)[j1], w2 = D.To(w)[j2];
               if ( D.From(w1).size( ) != 1 || D.To(w1).nonempty( ) ) continue;
               if ( D.From(w2).size( ) != 1 || D.To(w2).nonempty( ) ) continue;
               int d1 = D.ITo(w,j1), d2 = D.ITo(w,j2);
               vec<int> &x1 = D.OMutable(d1), &x2 = D.OMutable(d2);
               if ( x1[0] < 0 || x2[0] < 0 ) continue;
               if ( x1[0] != x2[0] ) continue;
               if ( x1.solo( ) || x2.solo( ) ) continue; // shouldn't happen
               int rd1 = dinv[d1], rd2 = dinv[d2];
               if ( !IsUnique( d1, d2, rd1, rd2 ) ) continue;
               int com; // number of shared hb edges
               for ( com = 1; com < Min( x1.isize( ), x2.isize( ) ) - 1; com++ )
                    if ( x1[com] != x2[com] ) break;
               vec<int> y, ry;
               y.SetToSubOf( x1, 0, com );
               x1.SetToSubOf( x1, com, x1.isize( ) - com );
               x2.SetToSubOf( x2, com, x2.isize( ) - com );
               D.TransferEdges( w1, w2 );
               dinv.push_back( D.E( ) + 1 );
               dinv.push_back( D.E( ) );
               D.AddEdge( w1, w2, y );
               vec<int> &rx1 = D.OMutable(rd1), &rx2 = D.OMutable(rd2);
               ry.SetToSubOf( rx1, rx1.isize( ) - com, com );
               rx1.SetToSubOf( rx1, 0, rx1.isize( ) - com );
               rx2.SetToSubOf( rx2, 0, rx2.isize( ) - com );
               int rw1 = to_right[rd1], rw2 = to_right[rd2];
               D.TransferEdges( rw1, rw2 );
               D.AddEdge( rw2, rw1, ry );
               goto next_w;    }
          next_w: continue;    }    }

void Emanate2( const vec<int>& inv, digraphE<vec<int>>& D, vec<int>& dinv,
     const Bool verbose )
{    if (verbose) cout << Date( ) << ": start emanate2" << endl;
     int edits = 0;
     vec<int> to_right;
     D.ToRight(to_right);
     int N = D.N( );
     for ( int w = 0; w < N; w++ )
     {    if ( D.To(w).size( ) != 2 ) continue;
          int w1 = D.To(w)[0], w2 = D.To(w)[1];
          if ( D.From(w1).size( ) != 1 || D.To(w1).nonempty( ) ) continue;
          if ( D.From(w2).size( ) != 1 || D.To(w2).nonempty( ) ) continue;
          int d1 = D.ITo(w,0), d2 = D.ITo(w,1);
          vec<int> &x1 = D.OMutable(d1), &x2 = D.OMutable(d2);
          if ( x1.solo( ) || x2.solo( ) ) continue;
          if ( x1[0] < 0 || x2[0] < 0 ) continue;
          int rd1 = dinv[d1], rd2 = dinv[d2];
          if ( !IsUnique( d1, d2, rd1, rd2 ) ) continue;
          if ( x1[0] == x2[0] ) continue;
          vec<int> pos1, pos2;
          for ( int j = 0; j < x2.isize( ); j++ )
               if ( x1[0] == x2[j] ) pos2.push_back(j);
          for ( int j = 0; j < x1.isize( ); j++ )
               if ( x2[0] == x1[j] ) pos1.push_back(j);
          if ( pos1.size( ) + pos2.size( ) != 1 ) continue;
          if ( pos1.size( ) == 1 )
          {    int p1 = pos1[0];
               if ( p1 == x1.isize( ) - 1 ) continue;
               edits++;
               int com; // number of shared hb edges
               for ( com = 1; com < Min( x1.isize( ) - p1 - 1, x2.isize( ) ) - 1;
                    com++ )
               {    if ( x1[p1+com] != x2[com] ) break;    }
               vec<int> y;
               y.SetToSubOf( x1, 0, p1 + com );
               D.TransferEdges( w2, w1 );
               dinv.push_back( D.E( ) + 1 );
               dinv.push_back( D.E( ) );
               x1.SetToSubOf( x1, p1 + com, x1.isize( ) - p1 - com );
               x2.SetToSubOf( x2, com, x2.isize( ) - com );
               vec<int> &rx1 = D.OMutable(rd1), &rx2 = D.OMutable(rd2);
               int rw1 = to_right[rd1], rw2 = to_right[rd2];
               vec<int> ry = y;
               rx1 = x1, rx2 = x2;
               ry.ReverseMe( ), rx1.ReverseMe( ), rx2.ReverseMe( );
               for ( int j = 0; j < ry.isize( ); j++ )
                    ry[j] = inv[ ry[j] ];
               for ( int j = 0; j < rx1.isize( ); j++ )
                    rx1[j] = inv[ rx1[j] ];
               for ( int j = 0; j < rx2.isize( ); j++ )
                    rx2[j] = inv[ rx2[j] ];
               D.AddEdge( w2, w1, y );
               D.TransferEdges( rw2, rw1 );
               D.AddEdge( rw1, rw2, ry );     }
          if ( pos2.size( ) == 1 ) // duplicates above with 1 <--> 2
          {    int p2 = pos2[0];
               if ( p2 == x2.isize( ) - 1 ) continue;
               edits++;
               int com; // number of shared hb edges
               for ( com = 1; com < Min( x2.isize( ) - p2 - 1, x1.isize( ) ) - 1;
                    com++ )
               {    if ( x2[p2+com] != x1[com] ) break;    }
               vec<int> y;
               y.SetToSubOf( x2, 0, p2 + com );
               D.TransferEdges( w1, w2 );
               dinv.push_back( D.E( ) + 1 );
               dinv.push_back( D.E( ) );
               x2.SetToSubOf( x2, p2 + com, x2.isize( ) - p2 - com );
               x1.SetToSubOf( x1, com, x1.isize( ) - com );
               vec<int> &rx2 = D.OMutable(rd2), &rx1 = D.OMutable(rd1);
               int rw2 = to_right[rd2], rw1 = to_right[rd1];
               vec<int> ry = y;
               rx2 = x2, rx1 = x1;
               ry.ReverseMe( ), rx1.ReverseMe( ), rx2.ReverseMe( );
               for ( int j = 0; j < ry.isize( ); j++ )
                    ry[j] = inv[ ry[j] ];
               for ( int j = 0; j < rx2.isize( ); j++ )
                    rx2[j] = inv[ rx2[j] ];
               for ( int j = 0; j < rx1.isize( ); j++ )
                    rx1[j] = inv[ rx1[j] ];
               D.AddEdge( w1, w2, y );
               D.TransferEdges( rw1, rw2 );
               D.AddEdge( rw2, rw1, ry );     }    }
     if (verbose) cout << Date( ) << ": made " << edits << " edits" << endl;    }

// LineProx.  Find lines that appear to be proximate to other lines.  This only
// looks at the ends of lines, which is both a plus and a minus.

void LineProx( const HyperBasevectorX& hb, const vec<int>& inv,
     const VecIntVec& ebcx, const digraphE<vec<int>>& D, const vec<int>& dinv,
     const vec<vec<vec<vec<int>>>>& dlines,
     const vec< triple<int,int,int> >& qept, vec< vec< pair<int,int> > >& lhood )
{
     // Heuristics.

     const int LOOK_IN = 10000;

     // Stuff in parallel.

     vec<int> mult, tod( hb.E( ), -1 ), tol( D.E( ), -1 ), linv;
     #pragma omp parallel sections
     {
          // Index edges.

          #pragma omp section
          {    ComputeMult( hb, D, mult );
               for ( int d = 0; d < D.E( ); d++ )
               {    if ( D.O(d)[0] < 0 ) continue;
                    for ( int i = 0; i < D.O(d).isize( ); i++ )
                         tod[ D.O(d)[i] ] = d;    }    }

          // Index lines.

          #pragma omp section
          {    for ( int i = 0; i < dlines.isize( ); i++ )
               for ( int j = 0; j < dlines[i].isize( ); j++ )
               for ( int k = 0; k < dlines[i][j].isize( ); k++ )
               for ( int l = 0; l < dlines[i][j][k].isize( ); l++ )
               {    int e = dlines[i][j][k][l];
                    if ( e >= 0 ) tol[e] = i;    }    }

          // Get line inversion.

          #pragma omp section
          {    LineInv( dlines, dinv, linv );    }    }

     // Compute barcode sets for each line, looking only at the ends.

     cout << Date( ) << ": compute barcode sets for lines" << endl;
     vec<vec<int>> lbc( dlines.size( ) );
     #pragma omp parallel for
     for ( int i = 0; i < dlines.isize( ); i++ )
     {    const vec<vec<vec<int>>>& L = dlines[i];
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
                              if ( pos + len <= LOOK_IN && mult[e] == 1 )
                              {    for ( auto b : ebcx[e] ) lbc[i].push_back(b);    }
                              len += hb.Kmers(e);    }    }
                    lensj.push_back(len);    }
               Sort(lensj);
               if ( lensj.nonempty( ) ) pos += Median(lensj);    }
          pos = 0;
          for ( int j = L.isize( ) - 1; j >= 0; j-- )
          {    vec<int> lensj;
               for ( int k = 0; k < L[j].isize( ); k++ )
               {    int len = 0;
                    for ( int l = L[j][k].isize( ) - 1; l >= 0; l-- )
                    {    int d = L[j][k][l];
                         if ( D.O(d)[0] < 0 ) continue;
                         for ( int m = D.O(d).isize( ) - 1; m >= 0; m-- )
                         {    int e = D.O(d)[m];
                              if ( pos + len <= LOOK_IN && mult[e] == 1 )
                              {    for ( auto b : ebcx[e] ) lbc[i].push_back(b);    }
                              len += hb.Kmers(e);    }    }
                    lensj.push_back(len);    }
               Sort(lensj);
               if ( lensj.nonempty( ) ) pos += Median(lensj);    }
          UniqueSort( lbc[i] );    }

     // Index barcode links.

     cout << Date( ) << ": indexing barcode links" << endl;
     vec<int> qep_index( hb.E( ) + 1, -1 );
     qep_index.back( ) = qept.size( );
     for ( int64_t i = qept.jsize( ) - 1; i >= 0; i-- )
          qep_index[ qept[i].first ] = i;
     for ( int e = hb.E( ) - 1; e >= 0; e-- )
          if ( qep_index[e] < 0 ) qep_index[e] = qep_index[e+1];

     // Find barcode neighbors.

     cout << Date( ) << ": finding barcode neighbors" << endl;
     lhood.clear_and_resize( dlines.size( ) );
     const int MIN_LEN = 100;
     const int MIN_LINKS = 6;
     const double MIN_NHOOD_FRAC = 0.1;
     #pragma omp parallel for
     for ( int i1 = 0; i1 < dlines.isize( ); i1++ )
     {    if ( linv[i1] < i1 ) continue;
          vec<int> n;
          const vec<vec<vec<int>>>& L = dlines[i1];
          for ( int j = 0; j < L.isize( ); j++ )
          for ( int k = 0; k < L[j].isize( ); k++ )
          for ( int l = 0; l < L[j][k].isize( ); l++ )
          {    int d = L[j][k][l];
               if ( D.O(d)[0] < 0 ) continue;
               for ( int m = 0; m < D.O(d).isize( ); m++ )
               {    int e = D.O(d)[m];
                    if ( mult[e] != 1 ) continue;
                    e = Min( e, inv[e] );
                    for ( int64_t l = qep_index[e]; l < qep_index[e+1]; l++ )
                    {    int f = qept[l].second;
                         if ( hb.Kmers(f) < MIN_LEN || mult[f] != 1 ) continue;
                         int i2a = tol[tod[f]];
                         if ( i2a < 0 ) continue;
                         int i2 = Min( i2a, linv[i2a] );
                         if ( i2 != i1 ) n.push_back(i2);    }    }    }
          UniqueSort(n);
          for ( int j = 0; j < n.isize( ); j++ )
          {    int c = MeetSize( lbc[i1], lbc[ n[j] ] );
               if ( c >= MIN_LINKS ) lhood[i1].push( c, n[j] );    }
          ReverseSort( lhood[i1] );
          for ( int j = 1; j < lhood[i1].isize( ); j++ )
          {    if ( double( lhood[i1][j].first ) / lhood[i1][0].first
                    < MIN_NHOOD_FRAC )
               {    lhood[i1].resize(j);
                    break;    }    }
          if ( linv[i1] > i1 ) lhood[ linv[i1] ] = lhood[i1];    }    }

int LineN50( const HyperBasevectorX& hb, const digraphE<vec<int>>& D,
     const vec<vec<vec<vec<int>>>>& dlines )
{    vec<int> llens;
     GetLineLengths( hb, D, dlines, llens );
     return ( llens.size( ) > 0 ? N50(llens) : 0 );    }

void KillInversionArtifacts( const digraphE<vec<int>>& D, const vec<int>& dinv,
     const ReadPathVec& dpaths, const IntIndex& dpaths_index, 
     const vec<int64_t>& bid, vec<int>& dels, const int MAX_CAN_INS_DEL )
{
     // Kill low depth canonical inversions.

     cout << Date( ) << ": killing low depth canonical inversions" << endl;
     const int MIN_CAN_INS_RATIO = 5;
     int ndels = 0;
     Bool verbose = False;
     #pragma omp parallel for
     for ( int v = 0; v < D.N( ); v++ )
     {    if ( !D.To(v).solo( ) || D.From(v).size( ) != 2 ) continue;
          for ( int j1 = 0; j1 < 2; j1++ )
          {    int j2 = 1 - j1;
               int w = D.From(v)[j1];
               if ( !D.From(w).solo( ) ) continue;
               int e = D.IFrom(w,0), z = D.IFrom(v,j1), f = D.IFrom(v,j2);
               int h = D.ITo(v,0);
               if (verbose)
               {
                    #pragma omp critical
                    {    cout << "trying z = " << z << " = " << printSeq( D.O(z) ) 
                              << ", e = " << e << ", z = " << z
                              << ", f = " << f << ", h = " << h << endl;    }    }

               // Do we seem to be at an inversion?  Look for inverse match 
               // between top and bottom "lines", near this point, using an
               // arbitrary threshold for "near".

               vec<vec<int>> nhood = { {h,f}, {e} };
               if ( D.To(w).size( ) == 2 )
               {    for ( int l = 0; l < D.To(w).isize( ); l++ )
                    {    if ( D.ITo(w,l) != z ) 
                         {    nhood[1].push_back( D.ITo(w,l) );
                              int m = D.To(w)[l];
                              if ( D.To(m).size( ) == 2 && D.To(m)[0] == D.To(m)[1] )
                              {    int p = D.To(m)[0];
                                   if ( D.To(p).solo( ) )
                                   {    nhood[1].push_back( 
                                             D.ITo(p,0) );    }    }    }    }    }
          
               Bool looks_like = False;
               for ( auto x : nhood[0] )
               for ( auto y : nhood[1] )
                    if ( dinv[x] == y ) looks_like = True;
               if (verbose)
               {
                    #pragma omp critical
                    {    cout << "z = " << z << ", nhood[0] = "
                              << printSeq(nhood[0]) << ", nhood[1] = "
                              << printSeq(nhood[1]) 
                              << ", looks_like = " << int(looks_like)
                              << endl;    }    }
               if ( !looks_like ) continue;

               // Compute support.

               vec<int> j = {j1, j2};
               vec<vec<int64_t>> bids(2);
               for ( int q = 0; q < 2; q++ )
               {    int g = D.IFrom( v, j[q] );
                    for ( int l = 0; l < (int) dpaths_index.Count(g); l++ )
                    {    int64_t id = dpaths_index.Val(g,l);
                         const ReadPath& p = dpaths[id];

                         if ( verbose )
                         {
                              #pragma omp critical
                              {    cout << "fw, z = " << z << ", q = " << q
                                        << ", id = " << id << ", bid = " << bid[id]
                                        << ", p = " << printSeq(p) << endl;    }    }

                         if ( bid[id] == 0 ) continue;
                         Bool ok = False;
                         for ( int m = 0; m < (int) p.size( ) - 1; m++ )
                         {    if ( p[m] == h && p[m+1] == g )
                              {    ok = True;
                                   break;    }    }
                         if (ok) bids[q].push_back( bid[id] );    }
                    int rg = dinv[g], rh = dinv[h];
                    for ( int l = 0; l < (int) dpaths_index.Count(rg); l++ )
                    {    int64_t id = dpaths_index.Val(rg,l);
                         const ReadPath& p = dpaths[id];

                         if ( verbose )
                         {
                              #pragma omp critical
                              {    cout << "rc, z = " << z << ", q = " << q
                                        << ", id = " << id << ", bid = " << bid[id]
                                        << ", p = " << printSeq(p) << endl;    }    }

                         if ( bid[id] == 0 ) continue;
                         Bool ok = False;
                         for ( int m = 0; m < (int) p.size( ) - 1; m++ )
                         {    if ( p[m] == rg && p[m+1] == rh )
                              {    ok = True;
                                   break;    }    }
                         if (ok) bids[q].push_back( bid[id] );    }
                    UniqueSort( bids[q] );    }
               if (verbose)
               {
                    #pragma omp critical
                    {    cout << "z = " << z 
                              << ", bids[0] = " << printSeq(bids[0])
                              << ", bids[1] = " << printSeq(bids[1]) 
                              << endl;    }    }
               if ( bids[0].isize( ) > MAX_CAN_INS_DEL ) continue;
               if ( bids[1].isize( ) < MIN_CAN_INS_RATIO * bids[0].isize( ) )
                    continue;
               if ( bids[1].empty( ) ) continue;
               #pragma omp critical
               {    dels.push_back( z, dinv[z] );
                    ndels += 2;    }
               break;    }    }
     cout << Date( ) << ": identified " << ndels << " edges to delete" << endl;    }

void SimpleHangs( const HyperBasevectorX& hb, digraphE<vec<int>>& D,
     const vec<int>& dinv, vec<int>& dels, const int MAX_KILL,
     const double MIN_RATIO, const Bool verbose, const Bool single )
{
     // Find lengths.

     vec<int> lens( D.E( ), 0 );
     int nthreads = ( single ? 1 : omp_get_max_threads( ) );
     #pragma omp parallel for num_threads(nthreads)
     for ( int e = 0; e < D.E( ); e++ )
     {    if ( D.O(e)[0] < 0 ) continue;
          for ( int j = 0; j < D.O(e).isize( ); j++ )
               lens[e] += hb.Kmers( D.O(e)[j] );    }

     // Get distances to end.

     vec<int> dfw;
     DistancesToEndArr( D, lens, MAX_KILL * MIN_RATIO, True, dfw );

     // Look for hanging ends.

     #pragma omp parallel for num_threads(nthreads)
     for ( int v = 0; v < D.N( ); v++ )
     {    for ( int j1 = 0; j1 < D.From(v).isize( ); j1++ )
          {    int w = D.From(v)[j1];
               if ( D.From(w).nonempty( ) || !D.To(w).solo( ) ) continue;
               int e1 = D.IFrom( v, j1 );
               const vec<int>& x1 = D.OFrom( v, j1 );
               if ( x1[0] < 0 ) continue;
               int n1 = lens[e1];
               if ( n1 > MAX_KILL ) continue;
               for ( int j2 = 0; j2 < D.From(v).isize( ); j2++ )
               {    int e2 = D.IFrom( v, j2 );
                    const vec<int>& x2 = D.OFrom( v, j2 );
                    int n2 = lens[ D.IFrom(v,j2) ] + dfw[ D.From(v)[j2] ];
                    if ( n2 < MIN_RATIO * n1 ) continue;
                    #pragma omp critical
                    {    dels.push_back( e1, dinv[e1] );    }    }    }    }    }

void ZapInversionBubbles( const digraphE<vec<int>>& D, const vec<int>& dinv,
     vec<int>& dels )
{    cout << Date( ) << ": zapping inversion bubbles" << endl;
     vec<vec<vec<vec<int>>>> dlines;
     FindLines( D, dinv, dlines, MAX_CELL_PATHS, MAX_CELL_DEPTH );
     int nzaps = 0;
     #pragma omp parallel for schedule( dynamic, 1000 )
     for ( int dl = 0; dl < dlines.isize( ); dl++ )
     {    const vec<vec<vec<int>>>& L = dlines[dl];
          for ( int i = 0; i < L.isize( ) - 2; i += 2 )
          {    if ( L[i].size( ) != 1 || L[i+2].size( ) != 1 ) continue;
               if ( L[i][0].size( ) != 1 || L[i+2][0].size( ) != 1 ) continue;
               if ( dinv[ L[i][0][0] ] != L[i+2][0][0] ) continue;
               #pragma omp critical
               {    nzaps++;
                    for ( int j = 0; j < L[i+1].isize( ); j++ )
                    for ( int k = 0; k < L[i+1][j].isize( ); k++ )
                    {    dels.push_back( L[i+1][j][k],
                              dinv[ L[i+1][j][k] ] );    }    }    }    }
     cout << Date( ) << ": zapped " << nzaps << " inversion bubbles" << endl;    }

void ZapMegaInversionBubbles( digraphE<vec<int>>& D, const vec<int>& dinv )
{    cout << Date( ) << ": zapping mega inversion bubbles" << endl;
     vec<vec<vec<vec<int>>>> dlines;
     FindLines( D, dinv, dlines, MAX_CELL_PATHS, MAX_CELL_DEPTH );
     vec<int> to_left, to_right, tol( D.E( ), -1 ), linv, splays;
     D.ToLeft(to_left), D.ToRight(to_right);
     for ( int i = 0; i < dlines.isize( ); i++ )
     for ( int j = 0; j < dlines[i].isize( ); j++ )
     for ( int k = 0; k < dlines[i][j].isize( ); k++ )
     for ( int l = 0; l < dlines[i][j][k].isize( ); l++ )
     {    int e = dlines[i][j][k][l];
          tol[e] = i;    }
     LineInv( dlines, dinv, linv );
     #pragma omp parallel for schedule( dynamic, 1000 )
     for ( int dl = 0; dl < dlines.isize( ); dl++ )
     {    const vec<vec<vec<int>>>& L = dlines[dl];
          int v = to_left[ L.front( )[0][0] ], w = to_right[ L.back( )[0][0] ];
          if ( D.To(v).solo( ) && D.From(v).size( ) == 2 )
          {    int l1 = tol[ D.IFrom(v,0) ], l2 = tol[ D.IFrom(v,1) ];
               if ( linv[l1] == l2 )
               {    
                    #pragma omp critical
                    {    splays.push_back( v, w );    }    }    }    }
     UniqueSort(splays);
     for ( auto v : splays ) D.SplayVertex(v);
     cout << Date( ) << ": zapped " << splays.size( )/2 
          << " mega inversion bubbles" << endl;    }

int64_t CheckSum( const HyperBasevectorX& hb, const vec<int>& inv,
     const digraphE<vec<int>>& D, const vec<int>& dinv )
{    int64_t C = hb.CheckSum( ) + D.CheckSum( );
     for ( int e = 0; e < inv.isize( ); e++ )
          C += (e+1) * (e+1) * (inv[e] + 1);
     for ( int d = 0; d < D.E( ); d++ )
     for ( int j = 0; j < D.O(d).isize( ); j++ )
          C += (d+1) * (j+1) * (D.O(d)[j]+1);
     for ( int d = 0; d < dinv.isize( ); d++ )
          C += (d+1) * (d+1) * (dinv[d] + 1);
     return C;    }

// MakeMerges: M = { ( (d1,p1), (d2,p2), n ) } where 
// interval [p1,p1+n) on superedge d1 matches interval [p2,p2+n) on superedge d2.
// It looks like merge set should be symmetric with respect to the involution,
// but not checked carefully.

int MakeMerges( const vec< triple< pair<int,int>, pair<int,int>, int > >& M,
     digraphE<vec<int>>& D, vec<int>& dinv, const Bool verbose )
{
     // Segment edges.

     if (verbose)
     {    cout << Date( ) << ": segmenting based on " 
               << ToStringAddCommas( M.size( ) ) << " matches" << endl;    }
     vec<Bool> touched( D.E( ), False ), to_use( M.size( ), False );
     vec<vec<int>> segs( D.E( ) );
     for ( int i = 0; i < M.isize( ); i++ )
     {    int d1 = M[i].first.first, j1 = M[i].first.second;
          int d2 = M[i].second.first, j2 = M[i].second.second;
          int n = M[i].third;
          if ( touched[d1] || touched[d2] ) continue;
          int rd1 = dinv[d1], rd2 = dinv[d2];
          touched[d1] = touched[d2] = touched[rd1] = touched[rd2] = True;
          to_use[i] = True;
          segs[d1] = { 0, j1, j1+n, D.O(d1).isize( ) };
          segs[d2] = { 0, j2, j2+n, D.O(d2).isize( ) };
          UniqueSort( segs[d1] ), UniqueSort( segs[d2] );
          for ( int j = 0; j < segs[d1].isize( ); j++ )
               segs[rd1].push_back( D.O(d1).isize( ) - segs[d1][j] );
          for ( int j = 0; j < segs[d2].isize( ); j++ )
               segs[rd2].push_back( D.O(d2).isize( ) - segs[d2][j] );
          segs[rd1].ReverseMe( ), segs[rd2].ReverseMe( );    }

     // Merge along matches.

     if (verbose) cout << Date( ) << ": merging" << endl;
     vec<int> to_left, to_right, dels;
     D.ToLeft(to_left), D.ToRight(to_right);
     int joins = 0;
     for ( int i = 0; i < M.isize( ); i++ )
     {    if ( !to_use[i] ) continue;
          int d1 = M[i].first.first, j1 = M[i].first.second;
          int d2 = M[i].second.first, j2 = M[i].second.second;
          if ( dinv[d1] < d1 ) continue;
          int n = M[i].third;
          int rd1 = dinv[d1], rd2 = dinv[d2];
          dels.push_back( d1, d2, rd1, rd2 );
          int v1 = to_left[d1], w1 = to_right[d1];
          int v2 = to_left[d2], w2 = to_right[d2];
          int rv1 = to_left[rd1], rw1 = to_right[rd1];
          int rv2 = to_left[rd2], rw2 = to_right[rd2];
          if ( !IsUnique( vec<int>{ v1, w1, v2, w2, rv1, rw1, rv2, rw2 } ) )
               continue;
          int rj1 = D.O(d1).isize( ) - j1 - n, rj2 = D.O(d2).isize( ) - j2 - n;
          joins += 2;
          const vec<int> &segs1 = segs[d1], &segs2 = segs[d2];
          vec<vec<int>> u1( segs1.isize( ) - 1 ), u2( segs2.isize( ) - 1 );
          for ( int i = 0; i < u1.isize( ); i++ )
               u1[i].SetToSubOf( D.O(d1), segs1[i], segs1[i+1] - segs1[i] );
          digraphE<vec<int>> U1( u1, digraphE<vec<int>>::EDGES_IN_LINE );
          for ( int i = 0; i < u2.isize( ); i++ )
               u2[i].SetToSubOf( D.O(d2), segs2[i], segs2[i+1] - segs2[i] );
          digraphE<vec<int>> U2( u2, digraphE<vec<int>>::EDGES_IN_LINE );
          const vec<int> &rsegs1 = segs[rd1], &rsegs2 = segs[rd2];
          vec<vec<int>> ru1( rsegs1.isize( ) - 1 ), ru2( rsegs2.isize( ) - 1 );
          for ( int i = 0; i < ru1.isize( ); i++ )
               ru1[i].SetToSubOf( D.O(rd1), rsegs1[i], rsegs1[i+1] - rsegs1[i] );
          digraphE<vec<int>> RU1( ru1, digraphE<vec<int>>::EDGES_IN_LINE );
          for ( int i = 0; i < ru2.isize( ); i++ )
               ru2[i].SetToSubOf( D.O(rd2), rsegs2[i], rsegs2[i+1] - rsegs2[i] );
          digraphE<vec<int>> RU2( ru2, digraphE<vec<int>>::EDGES_IN_LINE );
          int N = D.N( ), E = D.E( );
          D.AppendWithUpdate( U1, to_left, to_right ); 
          D.AppendWithUpdate( RU1, to_left, to_right ); 
          int n1 = U1.E( ), n2 = U2.E( ), z1 = U1.N( ), z2 = U2.N( );
          for ( int i = 0; i < U1.E( ); i++ )
               dinv.push_back( E + U1.E( ) + n1 - i - 1 );
          for ( int i = 0; i < U1.E( ); i++ )
               dinv.push_back( E + n1 - i - 1 );
          E = D.E( );
          D.AppendWithUpdate( U2, to_left, to_right ); 
          D.AppendWithUpdate( RU2, to_left, to_right ); 
          for ( int i = 0; i < U2.E( ); i++ )
               dinv.push_back( E + U2.E( ) + n2 - i - 1 );
          for ( int i = 0; i < U2.E( ); i++ )
               dinv.push_back( E + n2 - i - 1 );
          D.TransferEdgesWithUpdate( v1, N, to_left, to_right );
          D.TransferEdgesWithUpdate( w1, N + z1 - 1, to_left, to_right );
          D.TransferEdgesWithUpdate( v2, N + 2*z1, to_left, to_right );
          D.TransferEdgesWithUpdate( w2, N + 2*z1 + z2 - 1, to_left, to_right );
          D.TransferEdgesWithUpdate( rv1, N + z1, to_left, to_right );
          D.TransferEdgesWithUpdate( rw1, N + 2*z1 - 1, to_left, to_right );
          D.TransferEdgesWithUpdate( rv2, N + 2*z1 + z2, to_left, to_right );
          D.TransferEdgesWithUpdate( rw2, N + 2*z1 + 2*z2 - 1, to_left, to_right );
          for ( int i1 = 0; i1 < segs1.isize( ); i1++ )
          {    int p1 = segs1[i1];
               if ( p1 != j1 && p1 != j1 + n ) continue;
               int p2 = p1 + j2 - j1;
               int i2 = Position( segs2, p2 );
               if ( i2 < 0 ) continue;
               D.TransferEdgesWithUpdate( 
                    N + i1, N + 2*z1 + i2, to_left, to_right );    }
          for ( int i1 = 0; i1 < rsegs1.isize( ); i1++ )
          {    int p1 = rsegs1[i1];
               if ( p1 != rj1 && p1 != rj1 + n ) continue;
               int p2 = p1 + rj2 - rj1;
               int i2 = Position( rsegs2, p2 );
               if ( i2 < 0 ) continue;
               D.TransferEdgesWithUpdate( N + z1 + i1, N + 2*z1 + z2 + i2, 
                    to_left, to_right );    }    }

     // Clean up.

     if (verbose)
     {    cout << Date( ) << ": found " << ToStringAddCommas(joins) 
               << " joins" << endl;    }
     D.DeleteEdges(dels);
     Zipper( D, dinv );
     return joins;    }

// Try to merge along relatively short overlaps that are proximate in the graph.

void MergeShortOverlaps( HyperBasevectorX& hb, vec<int>& inv, digraphE<vec<int>>& D,
     vec<int>& dinv, const int LOOK_MERGE, const Bool allow_two )
{
     cout << Date( ) << ": start micromerger" << endl;
     const int LOOK = 6;
     vec<int> to_left, to_right, lens( D.E( ), 0 );
     D.ToLeft(to_left), D.ToRight(to_right);
     #pragma omp parallel for
     for ( int e = 0; e < D.E( ); e++ )
     {    if ( D.O(e)[0] < 0 ) continue;
          for ( int j = 0; j < D.O(e).isize( ); j++ )
               lens[e] += hb.Kmers( D.O(e)[j] );    }
     cout << Date( ) << ": and begin the loop" << endl;
     vec< triple< pair<int,int>, pair<int,int>, int > > M;
     #pragma omp parallel for schedule(dynamic, 10000)
     for ( int v = 0; v < D.N( ); v++ )
     {    vec<vec<int>> paths;
          for ( int j = 0; j < D.From(v).isize( ); j++ )
               paths.push_back( { D.IFrom(v,j) } );
          for ( int m = 1; m <= LOOK - 1; m++ )
          {    vec<vec<int>> paths2;
               for ( int j = 0; j < paths.isize( ); j++ )
               {    const vec<int>& p = paths[j];
                    int w = to_right[ p.back( ) ];
                    if ( D.From(w).empty( ) ) paths2.push_back(p);
                    for ( int i = 0; i < D.From(w).isize( ); i++ )
                    {    int f = D.IFrom(w,i);
                         if ( Member( p, f ) ) continue;
                         vec<int> q(p);
                         q.push_back(f);
                         paths2.push_back(q);    }    }
               paths = paths2;    }
          for ( int u = 0; u < LOOK; u++ )
          {    vec<int> see;
               for ( int j = 0; j < paths.isize( ); j++ )
               {    int n = paths[j].size( );
                    if ( u < n ) see.push_back( paths[j][u] );
                    if ( u-1 >= 0 && u-1 < n ) see.push_back( paths[j][u-1] );
                    if ( u-2 >= 0 && u-2 < n ) see.push_back( paths[j][u-2] );    }
               UniqueSort(see);
               for ( int j1 = 0; j1 < see.isize( ); j1++ )
               for ( int j2 = j1+1; j2 < see.isize( ); j2++ )
               {    int d1 = see[j1], d2 = see[j2];
                    if ( d1 == d2 ) continue;
                    if ( lens[d1] < LOOK_MERGE || lens[d2] < LOOK_MERGE ) continue;

                    // Look for an overlap of at least LOOK_SEE kmers between
                    // edges d1 and d2.  Give up if there is more than one.

                    const vec<int> &x1 = D.O(d1), &x2 = D.O(d2);
                    vec< triple<int,int,int> > com;
                    for ( int i = 0; i < x1.isize( ); i++ )
                         com.push( x1[i], i, 1 );
                    for ( int i = 0; i < x2.isize( ); i++ )
                         com.push( x2[i], i, 2 );
                    Sort(com);
                    vec< triple<int,int,int> > overs;
                    for ( int r = 0; r < com.isize( ); r++ )
                    {    int s;
                         for ( s = r + 1; s < com.isize( ); s++ )
                              if ( com[s].first != com[r].first ) break;
                         for ( int t1 = r; t1 < s; t1++ )
                         {    if ( com[t1].third != 1 ) continue;
                              for ( int t2 = r; t2 < s; t2++ )
                              {    if ( com[t2].third != 2 ) continue;
                                   int p1 = com[t1].second, p2 = com[t2].second;
                                   if ( p1 > 0 && p2 > 0 && x1[p1-1] == x2[p2-1] )
                                        continue;
                                   int n = 0, l1;
                                   for ( l1 = p1; l1 < x1.isize( ); l1++ )
                                   {    int l2 = l1 + p2 - p1;
                                        if ( l2 == x2.isize( ) ) break;
                                        if ( x1[l1] != x2[l2] ) break;
                                        n += hb.Kmers( x1[l1] );    }
                                   if ( n < LOOK_MERGE ) continue;
                                   int z = l1 - p1;
                                   overs.push( p1, p2, z );    }    }
                         r = s - 1;    }
                    if ( overs.size( ) == 0 || overs.size( ) > 2 ) continue;
                    if ( overs.size( ) == 2 && !allow_two ) continue;
                    if ( overs.size( ) == 2 && overs[1].third > overs[0].third ) 
                         swap( overs[0], overs[1] );

                    // Test for niceness and save.

                    int rd1 = dinv[d1], rd2 = dinv[d2];
                    if ( !IsUnique( vec<int>{ d1, d2, rd1, rd2 } ) ) continue;
                    int v1 = to_left[d1], w1 = to_right[d1];
                    int v2 = to_left[d2], w2 = to_right[d2];
                    int rv1 = to_left[rd1], rw1 = to_right[rd1];
                    int rv2 = to_left[rd2], rw2 = to_right[rd2];
                    if ( !IsUnique( 
                         vec<int>{ v1, w1, v2, w2, rv1, rw1, rv2, rw2 } ) ) 
                    {    continue;    }
                    int p1 = overs[0].first, p2 = overs[0].second;
                    int n = overs[0].third;
                    int n1 = D.O(d1).size( ), n2 = D.O(d2).size( );
                    #pragma omp critical
                    {    M.push( make_pair( d1, p1 ), make_pair( d2, p2 ), n );    
                         M.push( make_pair( dinv[d1], n1-p1-n ),
                              make_pair( dinv[d2], n2-p2-n ), n );    }
                    goto next_v_candidate;    }    }
          next_v_candidate: continue;    }
     cout << Date( ) << ": sorting merges" << endl;
     ParallelUniqueSort(M);

     // Make the merges.

     int count = MakeMerges( M, D, dinv );
     cout << Date( ) << ": made " << count << " merges" << endl;    }

void LineLevelPullApart( const vec<int32_t>& bc, const HyperBasevectorX& hb,
     const vec<int>& inv, const ReadPathVec& paths, const vec<Bool>& dup,
     digraphE<vec<int>>& D, vec<int>& dinv )
{
     cout << Date( ) << ": start process of making line-level pullaparts" << endl;

     // Get lines and barcode positions on them.

     cout << Date( ) << ": get lines and barcode positions on them" << endl;
     vec<vec<vec<vec<int>>>> dlines;
     vec<int> llens, tol, linv;
     vec<vec<pair<int,int>>> lbp;
     vec<int> to_left, to_right;
     {    D.ToLeft(to_left), D.ToRight(to_right);
          FindLines( D, dinv, dlines, MAX_CELL_PATHS, MAX_CELL_DEPTH );
          GetLineLengths( hb, D, dlines, llens );
          tol.resize( D.E( ) );
          for ( int i = 0; i < dlines.isize( ); i++ )
          for ( int j = 0; j < dlines[i].isize( ); j++ )
          for ( int k = 0; k < dlines[i][j].isize( ); k++ )
          for ( int l = 0; l < dlines[i][j][k].isize( ); l++ )
               tol[ dlines[i][j][k][l] ] = i;
          LineInv( dlines, dinv, linv );
          ReadPathVec dpaths;
          PlaceReads( hb, paths, dup, D, dpaths, True, False );
          const int BC_VIEW = 50000;
          {    IntIndex dpaths_index( dpaths, D.E( ), False );
               BarcodePos( bc, hb, D, dinv, dpaths, dpaths_index, dlines, 
                    lbp, BC_VIEW );    }    }

     // Make line-level pullaparts.

     cout << Date( ) << ": finding line-level pullaparts" << endl;
     const int MIN_LONG = 2000;
     const int IGNORE = 1000;
     const int VIEW = 20000;
     const int MIN_SUPPORT_PLUS = 10;
     const int MAX_SUPPORT_MINUS = 2;
     int lpulls = 0;
     int N0 = D.N( );
     for ( int l1 = 0; l1 < dlines.isize( ); l1++ )
     {    
          // Look for {l1,l2}m{r1,r2} that should turn into l1mr1 + l2mr2.

          if ( llens[l1] < MIN_LONG ) continue;
          int d1 = dlines[l1].back( )[0][0];
          int v = to_right[d1], d2 = -1;
          if ( D.To(v).size( ) != 2 || D.From(v).size( ) != 1 ) continue;
          for ( int j : {0,1} ) if ( D.ITo(v,j) != d1 ) d2 = D.ITo(v,j);
          int l2 = tol[d2];
          if ( llens[l2] < MIN_LONG ) continue;
          int m = tol[ D.IFrom(v,0) ];
          int w = to_right[ dlines[m].back( )[0][0] ];
          if ( v == w ) continue;
          if ( D.From(w).size( ) != 2 ) continue;
          int f1 = D.IFrom(w,0), f2 = D.IFrom(w,1);
          int r1 = tol[f1], r2 = tol[f2];
          int rl1 = linv[l1], rl2 = linv[l2];
          int rm = linv[m], rr1 = linv[r1], rr2 = linv[r2];
          if ( !IsUnique( vec<int>{ l1, l2, m, r1, r2, rl1, rl2, rm, rr1, rr2 } ) ) 
               continue;
          if ( llens[r1] < MIN_LONG || llens[r2] < MIN_LONG ) continue;
          vec<int32_t> bc_l1, bc_l2, bc_r1, bc_r2;
          for ( auto x : lbp[l1] )
          {    int b = x.first, p = x.second;
               if ( llens[l1] - p >= IGNORE && llens[l1] - p <= VIEW
                    && ( bc_l1.empty( ) || b != bc_l1.back( ) ) )
                    {    bc_l1.push_back(b);    }    }
          for ( auto x : lbp[l2] )
          {    int b = x.first, p = x.second;
               if ( llens[l2] - p >= IGNORE && llens[l2] - p <= VIEW
                    && ( bc_l2.empty( ) || b != bc_l2.back( ) ) )
               {    bc_l2.push_back(b);    }    }
          for ( auto x : lbp[r1] )
          {    int b = x.first, p = x.second;
               if ( p >= IGNORE && p <= VIEW
                    && ( bc_r1.empty( ) || b != bc_r1.back( ) ) )
               {    bc_r1.push_back(b);    }    }
          for ( auto x : lbp[r2] )
          {    int b = x.first, p = x.second;
               if ( p >= IGNORE && p <= VIEW
                    && ( bc_r2.empty( ) || b != bc_r2.back( ) ) )
               {    bc_r2.push_back(b);    }    }
          int u11 = MeetSize( bc_l1, bc_r1 ), u12 = MeetSize( bc_l1, bc_r2 );
          int u21 = MeetSize( bc_l2, bc_r1 ), u22 = MeetSize( bc_l2, bc_r2 );
          if ( u12 > u11 )
          {    swap(r1,r2), swap(u11,u12), swap(u21,u22); swap(f1,f2);    }
          if ( u11 < MIN_SUPPORT_PLUS || u22 < MIN_SUPPORT_PLUS ) continue;
          if ( u12 > MAX_SUPPORT_MINUS || u21 > MAX_SUPPORT_MINUS ) continue;

          // Make the edit.

          vec<int> em = Contents( dlines[m] ), mto_left, mto_right;
          digraphE<vec<int>> M( digraphE<vec<int>>::COMPLETE_SUBGRAPH_EDGES,
               D, em, to_left, to_right );
          M.ToLeft(mto_left), M.ToRight(mto_right);
          int N = D.N( ), E = D.E( );
          int rd1 = dinv[d1], rd2 = dinv[d2], rf1 = dinv[f1], rf2 = dinv[f2];
          D.Append(M);
          int u1 = dlines[m].front( )[0][0], u2 = dlines[m].back( )[0][0];
          int p1 = BinPosition( em, u1 ), p2 = BinPosition( em, u2 );
          D.GiveEdgeNewToVx( d2, to_right[d2], N + mto_left[p1] );
          D.GiveEdgeNewFromVx( f2, to_left[f2], N + mto_right[p2] );
          vec<int> erm = Contents( dlines[rm] ), rmto_left, rmto_right;
          digraphE<vec<int>> RM( digraphE<vec<int>>::COMPLETE_SUBGRAPH_EDGES,
               D, erm, to_left, to_right );
          RM.ToLeft(rmto_left), RM.ToRight(rmto_right);
          D.Append(RM);
          int rp1 = BinPosition( erm, dinv[u1] ), rp2 = BinPosition( erm, dinv[u2] );
          D.GiveEdgeNewFromVx( rd2, to_left[rd2], N + M.N( ) + rmto_right[rp1] );
          D.GiveEdgeNewToVx( rf2, to_right[rf2], N + M.N( ) + rmto_left[rp2] );
          for ( int g = 0; g < M.E( ); g++ ) 
               dinv.push_back( E + M.E( ) + BinPosition( erm, dinv[ em[g] ] ) );
          for ( int g = 0; g < M.E( ); g++ ) 
               dinv.push_back( E + BinPosition( em, dinv[ erm[g] ] ) );
          lpulls++;    }
     cout << Date( ) << ": made " << lpulls << " line-level pullaparts" << endl;

     // Look for 2-1 splits.

     cout << Date( ) << ": looking for two-to-one splits" << endl;
     int splits = 0;
     vec<int> vs;
     vec<Bool> vs_swap;
     #pragma omp parallel for schedule(dynamic, 10000)
     for ( int v = 0; v < N0; v++ )
     {    if ( D.To(v).size( ) != 1 || D.From(v).size( ) != 2 ) continue;
          int d = D.ITo(v,0), f1 = D.IFrom(v,0), f2 = D.IFrom(v,1);
          int l = tol[d], r1 = tol[f1], r2 = tol[f2];
          int rd = dinv[d], rf1 = dinv[f1], rf2 = dinv[f2];
          if ( !IsUnique( vec<int>{ d, f1, f2, rd, rf1, rf2 } ) ) continue;
          if ( llens[l] < MIN_LONG ) continue;
          if ( llens[r1] < MIN_LONG || llens[r2] < MIN_LONG ) continue;
          vec<int32_t> bc_l, bc_r1, bc_r2;
          for ( auto x : lbp[l] )
          {    int b = x.first, p = x.second;
               if ( llens[l] - p >= IGNORE && llens[l] - p <= VIEW
                    && ( bc_l.empty( ) || b != bc_l.back( ) ) )
               {    bc_l.push_back(b);    }    }
          for ( auto x : lbp[r1] )
          {    int b = x.first, p = x.second;
               if ( p >= IGNORE && p <= VIEW
                    && ( bc_r1.empty( ) || b != bc_r1.back( ) ) )
               {    bc_r1.push_back(b);    }    }
          for ( auto x : lbp[r2] )
          {    int b = x.first, p = x.second;
               if ( p >= IGNORE && p <= VIEW
                    && ( bc_r2.empty( ) || b != bc_r2.back( ) ) )
               {    bc_r2.push_back(b);    }    }
          int u1 = MeetSize( bc_l, bc_r1 ), u2 = MeetSize( bc_l, bc_r2 );
          Bool swapped = False;
          if ( u1 > u2 )
          {    swapped = True;
               swap(r1,r2), swap(u1,u2), swap(f1,f2), swap(rf1,rf2);    }
          if ( u2 < MIN_SUPPORT_PLUS || u1 > MAX_SUPPORT_MINUS ) continue;
          #pragma omp critical
          {    vs.push_back(v);    
               vs_swap.push_back(swapped);    }    }
     SortSync( vs, vs_swap );
     for ( int i = 0; i < vs.isize( ); i++ )
     {    int v = vs[i];
          Bool swapped = vs_swap[i];
          if ( D.To(v).size( ) != 1 || D.From(v).size( ) != 2 ) continue;
          int d = D.ITo(v,0), f1 = D.IFrom(v,0), f2 = D.IFrom(v,1);
          int l = tol[d], r1 = tol[f1], r2 = tol[f2];
          int rd = dinv[d], rf1 = dinv[f1], rf2 = dinv[f2];
          if (swapped)
          {    swap(r1,r2), swap(f1,f2), swap(rf1,rf2);    }
          splits++;
          int z = to_left[f1], rz = to_right[rf1];
          int N = D.N( );
          D.AddVertices(2);
          D.GiveEdgeNewFromVx( f1, z, N );
          D.GiveEdgeNewToVx( rf1, rz, N+1 );    }
     Validate( hb, inv, D, dinv );
     cout << Date( ) << ": made " << splits << " two-to-one splits" << endl;

     // Clean up.

     RemoveUnneededVertices( D, dinv );
     CleanupCore( D, dinv );    }

void FlattenSomeBubbles( const HyperBasevectorX& hb, const vec<Bool>& dup,
     const ReadPathVec& paths, digraphE<vec<int>>& D, vec<int>& dinv, 
     ReadPathVec& dpaths, const int MAX_DELTA, const double MIN_RATIO,
     const int MAX_DEL )
{
     cout << Date( ) << ": flattening bubbles" << endl;
     vec<int> dels;
     PlaceReads( hb, paths, dup, D, dpaths, True, False );
     IntIndex dpaths_index( dpaths, D.E( ), False );
     vec<int> lens( D.E( ), 0 );
     #pragma omp parallel for
     for ( int e = 0; e < D.E( ); e++ )
     {    if ( D.O(e)[0] < 0 ) continue;
          for ( int j = 0; j < D.O(e).isize( ); j++ )
               lens[e] += hb.Kmers( D.O(e)[j] );    }
     #pragma omp parallel for
     for ( int v = 0; v < D.N( ); v++ )
     {    if ( D.To(v).size( ) != 1 || D.From(v).size( ) != 2 ) continue;
          if ( D.From(v)[0] != D.From(v)[1] ) continue;
          int w = D.From(v)[0];
          if ( v == w ) continue;
          if ( D.To(w).size( ) != 2 || D.From(w).size( ) != 1 ) continue;
          if ( D.O( D.IFrom(v,0) )[0] < 0 || D.O( D.IFrom(v,1) )[0] < 0 ) continue;
          int delta = AbsDiff( lens[ D.IFrom(v,0)], lens[ D.IFrom(v,1)] );
          if ( delta == 0 || delta > MAX_DELTA ) continue;
          vec<vec<int64_t>> supp(2);
          for ( int j = 0; j < 2; j++ )
          {    int d = D.IFrom(v,j);
               int rd = dinv[d];
               for ( int i = 0; i < dpaths_index.Count(d); i++ )
                    supp[j].push_back( dpaths_index.Val(d,i) / 2 );
               for ( int i = 0; i < dpaths_index.Count(rd); i++ )
                    supp[j].push_back( dpaths_index.Val(rd,i) / 2 );
               UniqueSort( supp[j] );    }
          int m = MeetSize( supp[0], supp[1] );
          for ( int j = 0; j < 2; j++ )
          {    int d1 = D.IFrom(v,j);
               int n1 = supp[j].isize( ) - m, n2 = supp[1-j].isize( ) - m;
               if ( n1 > MAX_DEL || n2 < MIN_RATIO * Max(1, n1) ) continue;
               #pragma omp critical
               {    dels.push_back( d1, dinv[d1] );    }    }    }
     #pragma omp parallel for
     for ( int v = 0; v < D.N( ); v++ )
     {    int n = D.From(v).size( );
          if ( D.To(v).size( ) != 1 || n < 2 ) continue;
          int w = D.From(v)[0];
          Bool bad = False;
          for ( int j = 1; j < n; j++ ) if ( D.From(v)[j] != w ) bad = True;
          if ( bad || v == w ) continue;
          if ( D.To(w).isize( ) != n || D.From(w).size( ) != 1 ) continue;
          for ( int j = 0; j < n; j++ )
               if ( D.O( D.IFrom(v,j) )[0] < 0 ) bad = True;
          if (bad) continue;
          vec<vec<int64_t>> supp(n);
          for ( int j = 0; j < n; j++ )
          {    int d = D.IFrom(v,j);
               int rd = dinv[d];
               for ( int i = 0; i < dpaths_index.Count(d); i++ )
                    supp[j].push_back( dpaths_index.Val(d,i) / 2 );
               for ( int i = 0; i < dpaths_index.Count(rd); i++ )
                    supp[j].push_back( dpaths_index.Val(rd,i) / 2 );
               UniqueSort( supp[j] );    }
          int M = 0;
          for ( int j = 0; j < n; j++ ) M = Max( M, supp[j].isize( ) );
          const double MIN_RATIO2 = 10.0;
          const int MAX_KILL = 2;
          for ( int j = 0; j < n; j++ )
          {    int d = D.IFrom(v,j), m = supp[j].size( );
               if ( m <= MAX_KILL && M >= MIN_RATIO2 * Max(m,1) )
               {
                    #pragma omp critical
                    {    dels.push_back( d, dinv[d] );    }    }    }    }
     UniqueSort(dels);
     cout << Date( ) << ": deleting " << ToStringAddCommas( dels.size( ) )
          << " weak bubble edges" << endl;
     D.DeleteEdges(dels);    }
