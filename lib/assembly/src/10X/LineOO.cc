// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"
#include "ParallelVecUtilities.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "10X/LineOO.h"
#include "10X/Super.h"

// returns (barcode, pos, pid) of reads that map to a line.
void ReadPosLine( const vec<int32_t>& bc, const HyperBasevectorX& hb,
     const digraphE<vec<int>>& D, const vec<int>& dinv,
     const ReadPathVec& dpaths, const IntIndex& dpaths_index,
     const vec<vec<vec<vec<int>>>>& dlines, vec<vec<triple<int32_t,int,int64_t>>>& lrpb,
     const int view )
{
     // Get superedge lengths.

     cout << Date( ) << ": getting superedge lengths a" << endl;
     vec<int> dlens( D.E( ), 0 );
     #pragma omp parallel for
     for ( int e = 0; e < D.E( ); e++ )
     {    if ( D.O(e)[0] < 0 ) continue;
          for ( int j = 0; j < D.O(e).isize( ); j++ )
               dlens[e] += hb.Kmers( D.O(e)[j] );    }

     // Get line lengths.

     vec<int> llens;
     GetLineLengths( hb, D, dlines, llens );

     // Define line buckets.

     cout << Date( ) << ": defining line buckets a" << endl;
     int nd = dlines.size( );
     vec<int> llensx(llens), ids( nd, vec<int>::IDENTITY );
     ParallelReverseSortSync( llensx, ids );
     vec<vec<int>> buckets;
     const int bsize = 1000000;
     for ( int i = 0; i < nd; i++ )
     {    int64_t sum = llensx[i];
          int j;
          for ( j = i + 1; j < nd; j++ )
          {    if ( sum >= bsize ) break;
               sum += llensx[j];    }
          vec<int> b;
          for ( int k = i; k < j; k++ )
               b.push_back( ids[k] );
          buckets.push_back(b);
          i = j - 1;    }

     // Compute barcode positions.

     cout << Date( ) << ": computing barcode positions on lines" << endl;
     lrpb.clear_and_resize( dlines.size( ) );
     #pragma omp parallel for schedule(dynamic, 1)
     for ( int bi = 0; bi < buckets.isize( ); bi++ )
     {    vec<int64_t> ids;
          vec<int> starts, lens;
          for ( int ib = 0; ib < buckets[bi].isize( ); ib++ )
          {    int dl = buckets[bi][ib];
               const vec<vec<vec<int>>>& L = dlines[dl];
               int pos = 0;
               ids.clear( );
               starts.clear( );
               for ( int m = 0; m < L.isize( ); m++ )
               {    const vec<vec<int>>& M = L[m];
                    int pos_add = 0;
                    {    lens.clear( );
                         for ( int k = 0; k < M.isize( ); k++ )
                         {    int len = 0;
                              for ( int l = 0; l < M[k].isize( ); l++ )
                              {    int d = M[k][l];
                                   if ( D.O(d)[0] < 0 ) continue;
                                   for ( int m = 0; m < D.O(d).isize( ); m++ )
                                        len += hb.Kmers( D.O(d)[m] );    }
                              lens.push_back(len);    }
                         Sort(lens);
                         if ( lens.nonempty( ) ) pos_add = Median(lens);    }
                    if ( view > 0 && pos > view && llens[dl] - pos - pos_add > view )
                    {    pos += pos_add;
                         continue;    }
                    for ( int i = 0; i < M.isize( ); i++ )
                    {    const vec<int>& x = M[i];
                         int len = 0;
                         for ( int j = 0; j < x.isize( ); j++ )
                         {    int d = x[j], rd = dinv[ x[j] ];
                              if ( D.O(d)[0] < 0 ) continue;
                              for ( int l = 0; l < dpaths_index.Count(d); l++ )
                              {    int64_t id = dpaths_index.Val(d,l);
                                   if ( bc[id] == 0 ) continue;
                                   const ReadPath& p = dpaths[id];
                                   int offset = p.getOffset( );
                                   for ( int j = 0; j < (int) p.size( ); j++ )
                                   {    if ( p[j] == d ) break;
                                        offset -= dlens[ p[j] ];    }
                                   int s = pos+len+offset;
                                   if ( view == 0 || s <= view 
                                        || llens[dl] - s <= view )
                                   {    ids.push_back(id);
                                   starts.push_back(s);    }    }
                              for ( int l = 0; 
                                   l < (int) dpaths_index.Count(rd); l++ )
                              {    int64_t id = dpaths_index.Val(rd,l);
                                   if ( bc[id] == 0 ) continue;
                                   const ReadPath& p = dpaths[id];
                                   int offset = p.getOffset( );
                                   for ( int j = 0; j < (int) p.size( ); j++ )
                                   {    if ( p[j] == rd ) break;
                                        offset -= dlens[ p[j] ];    }
                                   // crazy, see SecretOps.cc:
                                   offset = dlens[d] + 200 - 1 - offset;
                                   int s = pos+len+offset;
                                   if ( view == 0 || s <= view 
                                        || llens[dl] - s <= view )
                                   {    ids.push_back(id);
                                        starts.push_back(s);    }    }
                              len += dlens[d];    }    }
                    pos += pos_add;    }
               UniqueSortSync( ids, starts );
               SortSync( starts, ids );
               for ( int l = 0; l < ids.isize( ); l++ )
                    lrpb[dl].push( bc[ids[l]], starts[l], ids[l]/2 );
               Sort( lrpb[dl] );    }    }    }

// BarcodePos: find pairs (barcode,pos) on a line.
// If view > 0, only show (barcode,pos) within that distance from a line end.

void BarcodePos( const vec<int32_t>& bc, const HyperBasevectorX& hb,
     const digraphE<vec<int>>& D, const vec<int>& dinv,
     const ReadPathVec& dpaths, const IntIndex& dpaths_index,
     const vec<vec<vec<vec<int>>>>& dlines, vec<vec<pair<int,int>>>& lbp,
     const int view )
{
     // Get superedge lengths.

     cout << Date( ) << ": getting superedge lengths a" << endl;
     vec<int> dlens( D.E( ), 0 );
     #pragma omp parallel for
     for ( int e = 0; e < D.E( ); e++ )
     {    if ( D.O(e)[0] < 0 ) continue;
          for ( int j = 0; j < D.O(e).isize( ); j++ )
               dlens[e] += hb.Kmers( D.O(e)[j] );    }

     // Get line lengths.

     vec<int> llens;
     GetLineLengths( hb, D, dlines, llens );

     // Define line buckets.

     cout << Date( ) << ": defining line buckets a" << endl;
     int nd = dlines.size( );
     vec<int> llensx(llens), ids( nd, vec<int>::IDENTITY );
     ParallelReverseSortSync( llensx, ids );
     vec<vec<int>> buckets;
     const int bsize = 1000000;
     for ( int i = 0; i < nd; i++ )
     {    int64_t sum = llensx[i];
          int j;
          for ( j = i + 1; j < nd; j++ )
          {    if ( sum >= bsize ) break;
               sum += llensx[j];    }
          vec<int> b;
          for ( int k = i; k < j; k++ )
               b.push_back( ids[k] );
          buckets.push_back(b);
          i = j - 1;    }

     // Compute barcode positions.

     cout << Date( ) << ": computing barcode positions on lines" << endl;
     lbp.clear_and_resize( dlines.size( ) );
     #pragma omp parallel for schedule(dynamic, 1)
     for ( int bi = 0; bi < buckets.isize( ); bi++ )
     {    vec<int64_t> ids;
          vec<int> starts, lens;
          for ( int ib = 0; ib < buckets[bi].isize( ); ib++ )
          {    int dl = buckets[bi][ib];
               const vec<vec<vec<int>>>& L = dlines[dl];
               int pos = 0;
               ids.clear( );
               starts.clear( );
               for ( int m = 0; m < L.isize( ); m++ )
               {    const vec<vec<int>>& M = L[m];
                    int pos_add = 0;
                    {    lens.clear( );
                         for ( int k = 0; k < M.isize( ); k++ )
                         {    int len = 0;
                              for ( int l = 0; l < M[k].isize( ); l++ )
                              {    int d = M[k][l];
                                   if ( D.O(d)[0] < 0 ) continue;
                                   for ( int m = 0; m < D.O(d).isize( ); m++ )
                                        len += hb.Kmers( D.O(d)[m] );    }
                              lens.push_back(len);    }
                         Sort(lens);
                         if ( lens.nonempty( ) ) pos_add = Median(lens);    }
                    if ( view > 0 && pos > view && llens[dl] - pos - pos_add > view )
                    {    pos += pos_add;
                         continue;    }
                    for ( int i = 0; i < M.isize( ); i++ )
                    {    const vec<int>& x = M[i];
                         int len = 0;
                         for ( int j = 0; j < x.isize( ); j++ )
                         {    int d = x[j], rd = dinv[ x[j] ];
                              if ( D.O(d)[0] < 0 ) continue;
                              for ( int l = 0; l < dpaths_index.Count(d); l++ )
                              {    int64_t id = dpaths_index.Val(d,l);
                                   if ( bc[id] == 0 ) continue;
                                   const ReadPath& p = dpaths[id];
                                   int offset = p.getOffset( );
                                   for ( int j = 0; j < (int) p.size( ); j++ )
                                   {    if ( p[j] == d ) break;
                                        offset -= dlens[ p[j] ];    }
                                   int s = pos+len+offset;
                                   if ( view == 0 || s <= view 
                                        || llens[dl] - s <= view )
                                   {    ids.push_back(id);
                                   starts.push_back(s);    }    }
                              for ( int l = 0; 
                                   l < (int) dpaths_index.Count(rd); l++ )
                              {    int64_t id = dpaths_index.Val(rd,l);
                                   if ( bc[id] == 0 ) continue;
                                   const ReadPath& p = dpaths[id];
                                   int offset = p.getOffset( );
                                   for ( int j = 0; j < (int) p.size( ); j++ )
                                   {    if ( p[j] == rd ) break;
                                        offset -= dlens[ p[j] ];    }
                                   // crazy, see SecretOps.cc:
                                   offset = dlens[d] + 200 - 1 - offset;
                                   int s = pos+len+offset;
                                   if ( view == 0 || s <= view 
                                        || llens[dl] - s <= view )
                                   {    ids.push_back(id);
                                        starts.push_back(s);    }    }
                              len += dlens[d];    }    }
                    pos += pos_add;    }
               UniqueSortSync( ids, starts );
               SortSync( starts, ids );
               for ( int l = 0; l < ids.isize( ); l++ )
                    lbp[dl].push( bc[ ids[l] ], starts[l] );
               Sort( lbp[dl] );    }    }    }

double ScoreOrder( 
     const vec<int>& L,                    // list of line ids
     const vec<vec< pair<int,int> >>& lbp, // barcode positions on lines
     const vec<int>& llens,                // line lengths
     vec< triple<int,int,int> >& M )       // scratch
{
     // Create merged list.

     M.clear( );
     int pos = 0;
     for ( int i = 0; i < L.isize( ); i++ )
     {    int lid = L[i];
          const vec< pair<int,int> >& x = lbp[lid];
          for ( int z = 0; z < x.isize( ); z++ )
               M.push( x[z].first, i, pos + x[z].second );
          pos += llens[lid];    }
     Sort(M);

     // Score it.  Removing the top measurement (from one barcode)
     // does not appear to help.

     double ad = 0.0;
     vec<int> gaps;
     for ( int k = 0; k < M.isize( ); k++ )
     {    
          // Find one barcode group.

          int l;
          for ( l = k + 1; l < M.isize( ); l++ )
               if ( M[l].first != M[k].first ) break;

          // Process group for one barcode.

          double n = l - k - 1;
          double mean_gap = ( M[l-1].third - M[k].third ) / n;
          double MIN_ADD = 2.0; // increasing to 4.0 appeared worse
          for ( int z = k + 1; z < l; z++ )
          {    if ( M[z].second > M[z-1].second )
               {    double plus = 

                         // NOT SURE IF WE SHOULD HAVE THIS OR NOT
                         // Commenting out this line decreased chr error and
                         // and increased N50 scaffold size, but increased
                         // ori error.
                         // ( M[z].second - M[z-1].second ) *

                         double( M[z].third - M[z-1].third ) / mean_gap;
                    if ( plus >= MIN_ADD ) ad += plus;    }    }

          // Go on to next barcode.

          k = l - 1;    }

     return ad;    }

// DUPLICATES ABOVE:
double ScoreOrder( 
     const vec<int>& L,                    // list of line ids
     VirtualMasterVec<SerfVec<pair<int,int>>> lbpx,  // barcode positions on lines
     const vec<int>& llens,                // line lengths
     vec< triple<int,int,int> >& M )       // scratch
{
     // Create merged list.

     M.clear( );
     int pos = 0;
     for ( int i = 0; i < L.isize( ); i++ )
     {    int lid = L[i];
          const SerfVec< pair<int,int> >& x = lbpx[lid];
          for ( int z = 0; z < (int) x.size( ); z++ )
               M.push( x[z].first, i, pos + x[z].second );
          pos += llens[lid];    }
     Sort(M);

     // Score it.  Removing the top measurement (from one barcode)
     // does not appear to help.

     double ad = 0.0;
     vec<int> gaps;
     for ( int k = 0; k < M.isize( ); k++ )
     {    
          // Find one barcode group.

          int l;
          for ( l = k + 1; l < M.isize( ); l++ )
               if ( M[l].first != M[k].first ) break;

          // Process group for one barcode.

          double n = l - k - 1;
          double mean_gap = ( M[l-1].third - M[k].third ) / n;
          double MIN_ADD = 2.0; // increasing to 4.0 appeared worse
          for ( int z = k + 1; z < l; z++ )
          {    if ( M[z].second > M[z-1].second )
               {    double plus = 

                         // NOT SURE IF WE SHOULD HAVE THIS OR NOT
                         // Commenting out this line decreased chr error and
                         // and increased N50 scaffold size, but increased
                         // ori error.
                         // ( M[z].second - M[z-1].second ) *

                         double( M[z].third - M[z-1].third ) / mean_gap;
                    if ( plus >= MIN_ADD ) ad += plus;    }    }

          // Go on to next barcode.

          k = l - 1;    }

     return ad;    }

double MemoryScoreOrder( 
     const vec<int>& L,                    // list of line ids
     const vec<vec< pair<int,int> >>& lbp, // barcode positions on lines
     const vec<int>& llens,                // line lengths
     vec< triple<int,int,int> >& M,        // scratch
     map< vec<int>, double >& memory )
{ 
     double s = -1;
     #pragma omp critical
     {    if ( memory.find(L) != memory.end( ) )
               s = memory[L];    }
     if ( s < 0 )
     {    s = ScoreOrder( L, lbp, llens, M );
          #pragma omp critical
          {    memory[L] = s;    }    }
     return s;    }

void OrderN( 
     // lines to OO
     const vec<int>& L0,                    // list of line ids

     // assembly info:
     const vec<vec< pair<int,int> >>& lbp, // barcode positions on lines
     const vec<int>& llens,                // line lengths

     // working data structures
     vec< triple<int,int,int> >& M,

     // answer:
     double& advantage, vec<int>& L )
{
     // Create candidates.

     ForceAssertGe( L0.isize( ), 2 );
     int n = L0.size( );
     vec< vec<int> > Lt;
     {    Lt.push_back( { } );
          for ( int d = 0; d < n; d++ ) // extend by one at each step
          {    vec<vec<int>> L2;
               for ( int j = 0; j < Lt.isize( ); j++ ) // entry from last round
               for ( int k = 0; k < n; k++ ) // let's add L0[k]
               {    if ( Member( Lt[j], L0[k] ) ) continue;
                    vec<int> c = Lt[j];
                    c.push_back( L0[k] );
                    L2.push_back(c);    }
               Lt = L2;    }    }
     int N = Lt.size( );

     // Score candidates and report the winner.

     vec<double> ad( N, 0 );
     for ( int j = 0; j < N; j++ )
          ad[j] = ScoreOrder( Lt[j], lbp, llens, M );
     vec<int> ids( N, vec<int>::IDENTITY );
     SortSync( ad, ids );
     L = Lt[ ids[0] ];
     advantage = ad[1] - ad[0];    }

// DUPLICATES ABOVE:
void OrderN( 
     // lines to OO
     const vec<int>& L0,                    // list of line ids

     // assembly info:
     VirtualMasterVec<SerfVec<pair<int,int>>> lbpx,  // barcode positions on lines
     const vec<int>& llens,                // line lengths

     // working data structures
     vec< triple<int,int,int> >& M,

     // answer:
     double& advantage, vec<int>& L )
{
     // Create candidates.

     ForceAssertGe( L0.isize( ), 2 );
     int n = L0.size( );
     vec< vec<int> > Lt;
     {    Lt.push_back( { } );
          for ( int d = 0; d < n; d++ ) // extend by one at each step
          {    vec<vec<int>> L2;
               for ( int j = 0; j < Lt.isize( ); j++ ) // entry from last round
               for ( int k = 0; k < n; k++ ) // let's add L0[k]
               {    if ( Member( Lt[j], L0[k] ) ) continue;
                    vec<int> c = Lt[j];
                    c.push_back( L0[k] );
                    L2.push_back(c);    }
               Lt = L2;    }    }
     int N = Lt.size( );

     // Score candidates and report the winner.

     vec<double> ad( N, 0 );
     for ( int j = 0; j < N; j++ )
          ad[j] = ScoreOrder( Lt[j], lbpx, llens, M );
     vec<int> ids( N, vec<int>::IDENTITY );
     SortSync( ad, ids );
     L = Lt[ ids[0] ];
     advantage = ad[1] - ad[0];    }

void MemoryOrderN(
     // lines to OO
     const vec<int>& L0,                    // list of line ids

     // assembly info:
     const vec<vec< pair<int,int> >>& lbp, // barcode positions on lines
     const vec<int>& llens,                // line lengths

     // working data structures
     vec< triple<int,int,int> >& M,

     // memory
     map< vec<int>, double >& memory,      

     // answer:
     double& advantage, vec<int>& L )
{
     // Create candidates.

     ForceAssertGe( L0.isize( ), 2 );
     int n = L0.size( );
     vec< vec<int> > Lt;
     {    Lt.push_back( { } );
          for ( int d = 0; d < n; d++ ) // extend by one at each step
          {    vec<vec<int>> L2;
               for ( int j = 0; j < Lt.isize( ); j++ ) // entry from last round
               for ( int k = 0; k < n; k++ ) // let's add L0[k]
               {    if ( Member( Lt[j], L0[k] ) ) continue;
                    vec<int> c = Lt[j];
                    c.push_back( L0[k] );
                    L2.push_back(c);    }
               Lt = L2;    }    }
     int N = Lt.size( );

     // Score candidates and report the winner.

     vec<double> ad( N, 0 );
     for ( int j = 0; j < N; j++ )
          ad[j] = MemoryScoreOrder( Lt[j], lbp, llens, M, memory );
     vec<int> ids( N, vec<int>::IDENTITY );
     SortSync( ad, ids );
     L = Lt[ ids[0] ];
     advantage = ad[1] - ad[0];    }

void MemoryOrderAndOrientN( 
     // lines to OO
     const vec<int>& L0,                    // list of line ids

     // assembly info:
     const vec<vec< pair<int,int> >>& lbp, // barcode positions on lines
     const vec<int>& llens,                // line lengths
     const vec<int>& linv,                 // line inversion

     // working data structures
     vec< triple<int,int,int> >& M,

     // memory
     map< vec<int>, double >& memory,      

     // answer:
     double& advantage, vec<int>& L )
{
     // Create candidates.

     int n = L0.size( );
     vec< vec<int> > Lt;
     vec< vec<Bool> > fwt;
     {    Lt.push_back( { } );
          fwt.push_back( { } );
          for ( int d = 0; d < n; d++ ) // extend by one at each step
          {    vec<vec<int>> L2;
               vec<vec<Bool>> fw2;
               for ( int j = 0; j < Lt.isize( ); j++ ) // entry from last round
               for ( int k = 0; k < n; k++ ) // let's add L0[k]
               {    if ( Member( Lt[j], L0[k] ) ) continue;
                    for ( int p = 0; p < 2; p++ ) // two possible orientations
                    {    if ( k == 0 && p == 1 ) continue; // fix one orientation
                         vec<int> c = Lt[j];
                         c.push_back( L0[k] );
                         vec<Bool> b = fwt[j];
                         b.push_back( p == 0 ? True : False );
                         L2.push_back(c);
                         fw2.push_back(b);    }    }
               Lt = L2;
               fwt = fw2;    }    }
     int N = Lt.size( );

     // Score candidates and report the winner.

     vec<double> ad( N, 0 );
     for ( int j = 0; j < N; j++ )
     {    for ( int l = 0; l < Lt[j].isize( ); l++ )
               if ( !fwt[j][l] ) Lt[j][l] = linv[ Lt[j][l] ];
          ad[j] = MemoryScoreOrder( Lt[j], lbp, llens, M, memory );    }
     vec<int> ids( N, vec<int>::IDENTITY );
     SortSync( ad, ids );
     L = Lt[ ids[0] ];
     advantage = ad[1] - ad[0];    }
