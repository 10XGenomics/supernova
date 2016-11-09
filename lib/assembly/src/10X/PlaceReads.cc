// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>

#include "CoreTools.h"
#include "math/Functions.h"
#include "math/HoInterval.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "paths/long/large/Lines.h"
#include "10X/DfTools.h"
#include "10X/PlaceReads.h"

// Align a path to the supergraph.

void Align( 
     const digraphE<vec<int>>& D, // supergraph
     const vec<int>& to_left,     // edge to left
     const vec<int>& to_right,    // edge to right
     const vec<int>& x,           // input path
     ho_interval& xpos,           // aligned start/stop positions on x
     vec<int>& d,                 // path in D
     int& dstart,                 // index of start edge on first edge in d
     int& dstop                   // index+1 of stop edge on last edge in d
          )
{    
     // Extend right.

     while(1)
     {    if ( xpos.Stop( ) == x.isize( ) ) break;
          else if ( dstop < D.O( d.back( ) ).isize( ) )
          {    if ( x[ xpos.Stop( ) ] != D.O( d.back( ) )[dstop] ) break;
               else
               {    xpos.AddToStop(1);
                    dstop++;    }    }
          else
          {    int v = to_right[ d.back( ) ], dn = -1;
               for ( int j = 0; j < D.From(v).isize( ); j++ )
               {    int n = D.IFrom( v, j );
                    if ( D.O(n)[0] < 0 ) continue;
                    if ( D.O(n)[0] == x[ xpos.Stop( ) ] )
                    {    dn = n;
                         break;    }    }
               if ( dn < 0 ) break;
               xpos.AddToStop(1);
               d.push_back(dn);
               dstop = 1;    }    }

     // Extend left.

     d.ReverseMe( );
     while(1)
     {    if ( xpos.Start( ) == 0 ) break;
          else if ( dstart > 0 )
          {    if ( x[ xpos.Start( ) - 1 ] != D.O( d.back( ) )[dstart-1] ) break;
               {    xpos.AddToStart(-1);
                    dstart--;    }    }
          else
          {    int v = to_left[ d.back( ) ], dn = -1;
               for ( int j = 0; j < D.To(v).isize( ); j++ )
               {    int n = D.ITo( v, j );
                    if ( D.O(n)[0] < 0 ) continue;
                    if ( D.O(n).back( ) == x[ xpos.Start( ) - 1 ] )
                    {    dn = n;
                         break;    }    }
               if ( dn < 0 ) break;
               xpos.AddToStart(-1);
               d.push_back(dn);
               dstart = D.O(dn).isize( ) - 1;    }    }
     d.ReverseMe( );    }

// Align2: this version more correctly handles duplicate branches, as would arise
// in case zippering was incomplete.  Probably this should supplant Align at some
// point.

void Align2( 
     const digraphE<vec<int>>& D, // supergraph
     const vec<int>& to_left,     // edge to left
     const vec<int>& to_right,    // edge to right
     const vec<int>& x,           // input path
     ho_interval& xpos,           // aligned start/stop positions on x
     vec<int>& d,                 // path in D
     int& dstart,                 // index of start edge on first edge in d
     int& dstop                   // index+1 of stop edge on last edge in d
          )
{    
     // Extend right.

     while(1)
     {    if ( xpos.Stop( ) == x.isize( ) ) break;
          else if ( dstop < D.O( d.back( ) ).isize( ) )
          {    if ( x[ xpos.Stop( ) ] != D.O( d.back( ) )[dstop] ) break;
               else
               {    xpos.AddToStop(1);
                    dstop++;    }    }
          else
          {    int v = to_right[ d.back( ) ], dn = -1;
               int count = 0;
               for ( int j = 0; j < D.From(v).isize( ); j++ )
               {    int n = D.IFrom( v, j );
                    if ( D.O(n)[0] < 0 ) continue;
                    if ( D.O(n)[0] == x[ xpos.Stop( ) ] )
                    {    dn = n;
                         count++;    }    }
               if ( count == 0 ) break;
               if ( count == 1 )
               {    xpos.AddToStop(1);
                    d.push_back(dn);
                    dstop = 1;    }
               else
               {    vec<vec<int>> exts;
                    for ( int j = 0; j < D.From(v).isize( ); j++ )
                    {    int n = D.IFrom( v, j );
                         if ( D.O(n)[0] < 0 ) continue;
                         if ( D.O(n)[0] == x[ xpos.Stop( ) ] )
                         {    int v2 = to_right[n];
                              if ( D.From(v2).size( ) == 0 ) exts.push_back({n});
                              else for ( int l = 0; l < D.From(v2).isize( ); l++ )
                              {    vec<int> z = { n, D.IFrom(v2,l) };
                                   exts.push_back(z);    }    }    }
                    vec<int> matches;
                    int len = x.isize( ) - xpos.Stop( );
                    for ( int j = 0; j < exts.isize( ); j++ )
                    {    const vec<int>& z = exts[j];
                         vec<int> t;
                         for ( auto m : z )
                         {    for ( auto r : D.O(m) )
                              {    t.push_back(r);
                                   if ( t.isize( ) == len ) break;    }
                              if ( t.isize( ) == len ) break;    }
                         Bool ok = True;
                         for ( int l = 0; l < len; l++ )
                         {    if ( t[l] != x[ xpos.Stop( ) + l ] )
                              {    ok = False;
                                   break;    }    }
                         if (ok) matches.push_back(j);    }
                    if ( !matches.solo( ) ) break;
                    xpos.AddToStop(1);
                    d.push_back( exts[ matches[0] ][0] );
                    dstop = 1;    }    }    }

     // Extend left.

     d.ReverseMe( );
     while(1)
     {    if ( xpos.Start( ) == 0 ) break;
          else if ( dstart > 0 )
          {    if ( x[ xpos.Start( ) - 1 ] != D.O( d.back( ) )[dstart-1] ) break;
               {    xpos.AddToStart(-1);
                    dstart--;    }    }
          else
          {    int v = to_left[ d.back( ) ], dn = -1;
               int count = 0;
               for ( int j = 0; j < D.To(v).isize( ); j++ )
               {    int n = D.ITo( v, j );
                    if ( D.O(n)[0] < 0 ) continue;
                    if ( D.O(n).back( ) == x[ xpos.Start( ) - 1 ] )
                    {    dn = n;
                         count++;    }    }
               if ( count == 0 ) break;
               if ( count == 1 )
               {    xpos.AddToStart(-1);
                    d.push_back(dn);
                    dstart = D.O(dn).isize( ) - 1;    }
               else
               {    vec<vec<int>> exts;
                    for ( int j = 0; j < D.To(v).isize( ); j++ )
                    {    int n = D.ITo( v, j );
                         if ( D.O(n)[0] < 0 ) continue;
                         if ( D.O(n).back( ) == x[ xpos.Start( ) - 1 ] )
                         {    int v2 = to_left[n];
                              if ( D.To(v2).size( ) == 0 ) exts.push_back({n});
                              else for ( int l = 0; l < D.To(v2).isize( ); l++ )
                              {    vec<int> z = { n, D.ITo(v2,l) };
                                   exts.push_back(z);    }    }    }
                    vec<int> matches;
                    int len = xpos.Start( );
                    for ( int j = 0; j < exts.isize( ); j++ )
                    {    const vec<int>& z = exts[j];
                         vec<int> t;
                         for ( auto m : z )
                         {    for ( int k = D.O(m).isize( ) - 1; k >= 0; k-- )
                              {    int r = D.O(m)[k];
                                   t.push_back(r);
                                   if ( t.isize( ) == len ) break;    }
                              if ( t.isize( ) == len ) break;    }
                         Bool ok = True;
                         for ( int l = 0; l < len; l++ )
                         {    if ( t[l] != x[ xpos.Start( ) - l - 1 ] )
                              {    ok = False;
                                   break;    }    }
                         if (ok) matches.push_back(j);    }
                    if ( !matches.solo( ) ) break;
                    xpos.AddToStart(-1);
                    int n = exts[ matches[0] ][0];
                    d.push_back(n);
                    dstart = D.O(n).isize( ) - 1;    }    }    }
     d.ReverseMe( );    }

void PlaceReads( const HyperBasevectorX& hb, const ReadPathVec& paths, 
     const vec<Bool>& dup, const digraphE<vec<int>>& D, ReadPathVec& dpaths,
     const Bool verbose, const Bool single, const Bool align2 )
{
     // Compute dlens.  This is only used in the sanity check below.  It would be
     // good to eliminate this calculation if possible.

     vec<int> dlens( D.E( ), 0 );
     #pragma omp parallel for
     for ( int e = 0; e < D.E( ); e++ )
     {    if ( D.O(e)[0] < 0 ) continue;
          for ( int j = 0; j < D.O(e).isize( ); j++ )
               dlens[e] += hb.Kmers( D.O(e)[j] );    }

     // Index D.

     if (verbose) cout << Date( ) << ": creating index for read placement" << endl;
     vec<Bool> used;
     D.Used(used);
     vec<vec<pair<int,int>>> nd( hb.E( ) );
     for ( int e = 0; e < D.E( ); e++ )
     {    if ( !used[e] ) continue;
          const vec<int>& x = D.O(e);
          if ( x[0] < 0 ) continue;
          for ( int j = 0; j < x.isize( ); j++ )
               nd[ x[j] ].push( e, j );    }

     // Build paths.

     if ( dpaths.size( ) == 0 )
     {    if (verbose) cout << Date( ) << ": defining data structure" << endl;
          dpaths = paths;    } // TERRIBLE WAY TO RESERVE SPACE!!!!!!!!!!!!!!!!
     if (verbose) cout << Date( ) << ": creating to_left and to_right" << endl;
     vec<int> to_left, to_right;
     D.ToLeft(to_left), D.ToRight(to_right);
     if (verbose) cout << Date( ) << ": looking at reads" << endl;
     double clock = WallClockTime( );
     const int batch = 100000;
     int nthreads = ( single ? 1 : omp_get_max_threads( ) );
     #pragma omp parallel for num_threads(nthreads)
     for ( int64_t bi = 0; bi < (int64_t) paths.size( ); bi += batch )
     {    vec<int> aplace, d, x;
          for ( int64_t id = bi; 
               id < Min( bi + batch, (int64_t) paths.size( ) ); id++ )
          {    dpaths[id].resize(0); 
               if ( dup[id/2] ) continue;
               const ReadPath &p = paths[id]; 
               if ( p.size( ) == 0 ) continue;

               // Find minimum multiplicity edge in p.

               int m = 1000000000, mp = -1;
               for ( int j = 0; j < (int) p.size( ); j++ )
               {    int n = nd[ p[j] ].size( );
                    if ( n == 0 ) continue;
                    if ( n < m )
                    {    m = n;
                         mp = j;    }    }

               // Test placements.

               if ( mp < 0 || nd[ p[mp] ].empty( ) ) continue;
               int nplaces = 0, pos = -1;
               for ( int i = 0; i < nd[ p[mp] ].isize( ); i++ )
               {    int e = nd[ p[mp] ][i].first, epos = nd[ p[mp] ][i].second;
                    x.clear( );
                    for ( int j = 0; j < (int) p.size( ); j++ )
                         x.push_back( p[j] );
                    ho_interval xpos( mp, mp+1 );
                    d = {e};
                    int dstart = epos, dstop = epos + 1;
                    ( !align2 ? Align : Align2 )( 
                         D, to_left, to_right, x, xpos, d, dstart, dstop );
     
                    // For now, toss improper alignments.

                    if ( xpos.Start( ) > 0 )
                    {    if ( dstart > 0 ) continue;
                         if ( dstart == 0 )
                         {    int v = to_left[ d.front( ) ];
                              if ( D.To(v).nonempty( ) ) continue;    }    }
                    if ( xpos.Stop( ) < x.isize( ) ) 
                    {    if ( dstop < D.O( d.back( ) ).isize( ) ) continue;
                         int w = to_right[ d.back( ) ];
                         if ( D.From(w).nonempty( ) ) continue;    }

                    // Compute start position of read on d[0].

                    pos = p.getOffset( );
                    for ( int j = 0; j < dstart; j++ )
                         pos += hb.Kmers( D.O(d[0])[j] );

                    // Give up if more than one placement.

                    nplaces++;
                    if ( nplaces > 1 ) break;

                    // Save.  
                    // Just saving d, but could save ( xpos, d, dstart, dstop ).

                    aplace.clear( );
                    for ( int j = 0; j < d.isize( ); j++ )
                         aplace.push_back( d[j] );    }

               // For now require that there is a unique proper placement.

               if ( nplaces != 1 ) continue;

               // Check for broken placement.  Need to track down how this could
               // ever happen.

               int n = hb.K( ) - 1;
               for ( auto f : aplace ) n += dlens[f];
               if ( pos > n ) continue;

               // Save placement.

               dpaths[id].resize( aplace.size( ) );
               for ( int j = 0; j < aplace.isize( ); j++ )
                    dpaths[id][j] = aplace[j];
               dpaths[id].setOffset(pos);    }     }

     if (verbose) cout << Date( ) << ": done, used " << TimeSince(clock) << endl;
     int64_t placed = 0;
     for ( int64_t id = 0; id < (int64_t) paths.size( ); id++ )
          if ( dpaths[id].size( ) > 0 ) placed++;
     if (verbose)
     {    cout << Date( ) << ": " 
               << PERCENT_RATIO( 3, placed, (int64_t) paths.size( ) ) 
               << " placed" << endl;    }    }

void PlaceReadsSmart( const HyperBasevectorX& hb, const ReadPathVec& paths, 
     const vec<Bool>& dup, const digraphE<vec<int>>& D, const vec<int>& dinv,
     ReadPathVec& dpaths, const vec<vec<vec<vec<int>>>>& dlines, 
     const vec<int64_t>& bci, const Bool verbose, const int btest, 
     const Bool align2 )
{
     // Set up.

     vec<int> dlens( D.E( ), 0 ), to_left, to_right, linv;
     #pragma omp parallel for
     for ( int e = 0; e < D.E( ); e++ )
     {    if ( D.O(e)[0] < 0 ) continue;
          for ( int j = 0; j < D.O(e).isize( ); j++ )
               dlens[e] += hb.Kmers( D.O(e)[j] );    }
     D.ToLeft(to_left), D.ToRight(to_right);
     LineInv( dlines, dinv, linv );

     // To each nongap edge, assign a line, and a starting position on that line.
     // For nonpalindromic lines L, we pick just one of L, rc(L), and place both
     // edges d and dinv[d] at the same place on that line.

     vec<quad<int,int,Bool,int>> tol( D.E( ), make_quad(-1,-1,True,-1) );
     cout << Date( ) << ": defining tol" << endl;
     for ( int i = 0; i < dlines.isize( ); i++ )
     {    if ( linv[i] < i ) continue;
          const vec<vec<vec<int>>>& L = dlines[i];
          int pos = 0;
          for ( int j = 0; j < L.isize( ); j++ )
          {    vec<int> lens;
               for ( int k = 0; k < L[j].isize( ); k++ )
               {    int len = 0;
                    for ( int l = 0; l < L[j][k].isize( ); l++ )
                    {    int d = L[j][k][l];
                         if ( D.O(d)[0] < 0 ) continue;
                         if ( tol[d].first < 0 )
                         {    tol[d] = make_quad( i, pos+len, True, j );
                              tol[ dinv[d] ] 
                                   = make_quad( i, pos+len, False, j );    }
                         len += dlens[d];    }
                    lens.push_back(len);    }
               Sort(lens);
               if ( lens.nonempty( ) ) pos += Median(lens);    }    }

     // Define the territory of each barcode.  The granularity of this definition
     // is weak, as we use only start/stop positions (i,j) on lines, whereas we
     // should use base positions.

     // Go through the barcodes.

     cout << Date( ) << ": start traversal" << endl;
     const int MAX_BC_GAP = 100000;
     const int MIN_BC_GROUP = 3;
     int ndots = 0, stopped = 0;
     #pragma omp parallel for schedule(dynamic, 1000)
     for ( int b = 1; b < bci.isize( ) - 1; b++ )
     {    if ( btest >= 0 && b != btest ) continue;

          // Logging.

          ostringstream out;
          if (verbose)
          {    out << "\nlooking at barcode " << b << " [" << bci[b+1]-bci[b]
                    << " reads]" << endl;    }

          // Define the territory of each barcode.  The granularity of this
          // definition is weak, as we use only start/stop positions (i,j) on lines,
          // whereas we should use base positions.

          int total = 0;
          vec<triple<int,int,int>> lplaces; // (line, start on line, line unit)
          for ( int64_t id = bci[b]; id < bci[b+1]; id++ )
          {    const ReadPath& p = dpaths[id];
               int pos = p.getOffset( );
               for ( int j = 0; j < (int) p.size( ); j++ )
               {    if ( j == 0 || tol[ p[j] ].first != tol[ p[j-1] ].first )
                    {    triple<int,int,int> x; 
                         x.first = tol[ p[j] ].first;
                         int xpos = pos;
                         if ( !tol[ p[j] ].third ) xpos = dlens[ p[j] ] - xpos;
                         x.second = tol[ p[j] ].second + xpos;
                         x.third = tol[ p[j] ].fourth;
                         lplaces.push_back(x);    }
                    pos -= dlens[ p[j] ];    }    }
          Sort(lplaces);
          vec<triple<int,int,int>> groups;
          for ( int i = 0; i < lplaces.isize( ); i++ )
          {    int j;
               for ( j = i + 1; j < lplaces.isize( ); j++ )
               {    if ( lplaces[j].first != lplaces[i].first ) break;
                    if ( lplaces[j].second - lplaces[j-1].second > MAX_BC_GAP ) 
                         break;    }
               if ( j - i >= MIN_BC_GROUP )
               {    int len = lplaces[j-1].second - lplaces[i].second;
                    int start = lplaces[i].third, stop = lplaces[j-1].third;
                    int l = lplaces[i].first;
                    const vec<vec<vec<int>>>& L = dlines[l];

                    // Extend start and stop positions.

                    const int EXT = 25000;
                    int left = 0, right = 0;
                    while( left < EXT && start >= 1 )
                    {    start--;
                         vec<int> lens;
                         for ( int k = 0; k < L[start].isize( ); k++ )
                         {    int len = 0;
                              for ( int l = 0; l < L[start][k].isize( ); l++ )
                              {    int d = L[start][k][l];
                                   len += dlens[d];    }
                              lens.push_back(len);    }
                         Sort(lens);
                         if ( lens.nonempty( ) ) left += Median(lens);    }
                    while( right < EXT && stop < L.isize( ) - 1 )
                    {    stop++;
                         vec<int> lens;
                         for ( int k = 0; k < L[stop].isize( ); k++ )
                         {    int len = 0;
                              for ( int l = 0; l < L[stop][k].isize( ); l++ )
                              {    int d = L[stop][k][l];
                                   len += dlens[d];    }
                              lens.push_back(len);    }
                         Sort(lens);
                         if ( lens.nonempty( ) ) right += Median(lens);    }
                    
                    // Save.

                    if (verbose)
                    {    out << l << "." << lplaces[i].second << "-" 
                              << lplaces[j-1].second << " [len=" << len 
                              << ",count=" << j-i 
                              << ",start=" << start << "/" << L[start][0][0]
                              << ",stop=" << stop << "/" << L[stop][0][0]
                              << "]" << endl;    }
                    total += len;
                    groups.push( lplaces[i].first, start, stop );    }
               i = j - 1;    }    
          if (verbose) out << "total = " << ToStringAddCommas(total) << endl;

          // Define the edges in the territory.

          vec<pair<int,int>> locs;
          int nd = 0;
          for ( int i = 0; i < groups.isize( ); i++ )
          {    const vec<vec<vec<int>>>& L = dlines[ groups[i].first ];
               int start = groups[i].second, stop = groups[i].third;
               for ( int j = start; j <= stop; j++ )
               {    const vec<vec<int>>& M = L[j];
                    vec<int> ds;
                    for ( int k = 0; k < M.isize( ); k++ )
                    for ( int l = 0; l < M[k].isize( ); l++ )
                    {    int d = M[k][l];
                         if ( D.O(d)[0] >= 0 ) ds.push_back( d, dinv[d] );    }
                    UniqueSort(ds);
                    nd += ds.size( );
                    for ( int k = 0; k < ds.isize( ); k++ )
                    {    int d = ds[k];
                         vec<int> es = D.O(d);
                         UniqueSort(es);
                         for ( auto e : es ) locs.push( e, d );    }    }    }
          Sort(locs);
          if (verbose) 
          {    out << "locs = " << ToStringAddCommas( locs.size( ) ) 
                    << " base edges in " << ToStringAddCommas(nd) 
                    << " superedges" << endl;    }

          // Try to place more reads.

          int pl1 = 0, pl2 = 0;
          for ( int64_t id = bci[b]; id < bci[b+1]; id++ )
          {    if ( dpaths[id].size( ) > 0 ) pl1++;
               const ReadPath& p = paths[id];
               if ( p.size( ) == 0 || dpaths[id].size( ) > 0 || dup[id/2] ) continue;

               // Now p = paths[id] represents a nonduplicate read from the barcode
               // that is placed on the base graph.  Find supergraph edges assigned 
               // to the barcode (as defined by locs), and find the matches of the 
               // read to one those edges that extends for at least 30 kmers; also 
               // ignore matches that are less than half as long as the best match.
               // (Not sure the "30" is enforced below.)
               //
               // First let pl =
               // { (supergraph edge d, pos on p, pos on d, match length in kmers) }.

               vec<quad<int,int,int,int>> pl;
               for ( int j = 0; j < (int) p.size( ); j++ )
               {    int e = p[j];
                    int low = LowerBound1(locs, e), high = UpperBound1(locs, e);
                    for ( int m = low; m < high; m++ )
                    {    int d = locs[m].second;
                         for ( int l = 0; l < D.O(d).isize( ); l++ )
                         {    if ( D.O(d)[l] != e ) continue;
                              if ( l > 0 && j > 0 && D.O(d)[l-1] == p[j-1] ) 
                                   continue;
                              int n = 0, r;
                              for ( r = l; r < D.O(d).isize( ); r++ )
                              {    if ( r - l + j >= (int) p.size( ) ) break;
                                   if ( D.O(d)[r] != p[r-l+j] ) break;
                                   n += hb.Kmers( D.O(d)[r] );    }
                              if ( btest >= 0 ) PRINT5_TO( out, id, d, j, l, n );
                              pl.push( d, j, l, n );    }    }    } 

               // Now filter the matches as described above.

               int nmax = 0;
               for ( int i = 0; i < pl.isize( ); i++ )
                    nmax = Max( nmax, pl[i].fourth );
               vec<Bool> to_delete( pl.size( ), False );
               for ( int i = 0; i < pl.isize( ); i++ )
               {    int n = pl[i].fourth;
                    if ( nmax >= 30 && n < nmax/2 ) to_delete[i] = True;    }
               EraseIf( pl, to_delete );
               if ( btest >= 0 ) PRINT2_TO( out, id, pl.size( ) );

               // If there is a unique placement, place the read.

               if ( pl.solo( ) )
               {    int d = pl[0].first, j = pl[0].second, l = pl[0].third;
                    vec<int> x;
                    for ( auto e : p ) x.push_back(e);
                    ho_interval xpos( j, j + 1 );
                    vec<int> q = {d};
                    int qstart = l, qstop = l + 1;
                    ( !align2 ? Align : Align2 )( 
                         D, to_left, to_right, x, xpos, q, qstart, qstop );

                    // For now, toss improper alignments.

                    if ( btest >= 0 )
                    {    PRINT6_TO( out, id, xpos.Start( ), qstart, xpos.Stop( ),
                              x.size( ), qstop );    }
                    Bool bad = False;
                    if ( xpos.Start( ) > 0 )
                    {    if ( qstart > 0 ) continue;
                         if ( qstart == 0 )
                         {    int v = to_left[ q.front( ) ];
                              if ( D.To(v).nonempty( ) ) continue;    }    }
                    if ( xpos.Stop( ) < x.isize( ) )
                    {    if ( qstop < D.O( q.back( ) ).isize( ) ) continue;
                         int w = to_right[ q.back( ) ];
                         if ( D.From(w).nonempty( ) ) continue;    }

                    // Compute new offset.

                    int pos = p.getOffset( );
                    for ( int j = 0; j < qstart; j++ )
                         pos += hb.Kmers( D.O(q[0])[j] );

                    // Check for broken placements.  Not sure how this can happen.

                    int n = hb.K( ) - 1;
                    for ( auto f : q ) n += dlens[f];
                    if ( pos > n ) continue;

                    // Save the new placement.

                    dpaths[id].resize( q.size( ) );
                    for ( int j = 0; j < q.isize( ); j++ )
                         dpaths[id][j] = q[j];
                    // out << "placing read " << id << " at pos " << pos
                    //      << ", path = " << printSeq(q) << endl;
                    dpaths[id].setOffset(pos);    
                    pl2++;    }    }
          if (verbose) 
          {    out << "orig reads placed = " << pl1 << ", new reads placed = "
                    << pl2 << "\n";
               #pragma omp critical
               {    cout << out.str( );    }    }
          else
          {
               #pragma omp critical
               {    MakeDots( stopped, ndots, bci.isize( ) - 2 );    }    }    }    }


void PlaceReads( const HyperBasevectorX& hb, const ReadPathVecX& paths, 
     const vec<Bool>& dup, const digraphE<vec<int>>& D, ReadPathVec& dpaths,
     const Bool verbose, const Bool single, const Bool align2 )
{
     // Compute dlens.  This is only used in the sanity check below.  It would be
     // good to eliminate this calculation if possible.

     vec<int> dlens( D.E( ), 0 );
     #pragma omp parallel for
     for ( int e = 0; e < D.E( ); e++ )
     {    if ( D.O(e)[0] < 0 ) continue;
          for ( int j = 0; j < D.O(e).isize( ); j++ )
               dlens[e] += hb.Kmers( D.O(e)[j] );    }

     // Index D.

     if (verbose) cout << Date( ) << ": creating index for read placement" << endl;
     vec<Bool> used;
     D.Used(used);
     vec<vec<pair<int,int>>> nd( hb.E( ) );
     for ( int e = 0; e < D.E( ); e++ )
     {    if ( !used[e] ) continue;
          const vec<int>& x = D.O(e);
          if ( x[0] < 0 ) continue;
          for ( int j = 0; j < x.isize( ); j++ )
               nd[ x[j] ].push( e, j );    }

     // Build paths.

     if ( dpaths.size( ) == 0 )
     {    if (verbose) cout << Date( ) << ": defining data structure" << endl;
         // BETTER ALLOCATION NEEDED? !!!!!!!!!!!!!!!!!!!!!!!!!!
         dpaths.resize(paths.size());
         int mpl = 1;// heuristic
         for(int64_t id = 0; id<paths.size(); id++){
             dpaths[id].reserve(mpl); 
         }
     }
     if (verbose) cout << Date( ) << ": creating to_left and to_right" << endl;
     vec<int> to_left, to_right;
     D.ToLeft(to_left), D.ToRight(to_right);
     if (verbose) cout << Date( ) << ": looking at reads" << endl;
     double clock = WallClockTime( );
     const int batch = 100000;
     int nthreads = ( single ? 1 : omp_get_max_threads( ) );
     #pragma omp parallel for num_threads(nthreads)
     for ( int64_t bi = 0; bi < (int64_t) paths.size( ); bi += batch )
     {    vec<int> aplace, d, x;
          for ( int64_t id = bi; 
               id < Min( bi + batch, (int64_t) paths.size( ) ); id++ )
          {    dpaths[id].resize(0); 
               if ( dup[id/2] ) continue;
               ReadPath p;
               paths.unzip(p,hb,id);
               if ( p.size( ) == 0 ) continue;

               // Find minimum multiplicity edge in p.

               int m = 1000000000, mp = -1;
               for ( int j = 0; j < (int) p.size( ); j++ )
               {    int n = nd[ p[j] ].size( );
                    if ( n == 0 ) continue;
                    if ( n < m )
                    {    m = n;
                         mp = j;    }    }

               // Test placements.

               if ( mp < 0 || nd[ p[mp] ].empty( ) ) continue;
               int nplaces = 0, pos = -1;
               for ( int i = 0; i < nd[ p[mp] ].isize( ); i++ )
               {    int e = nd[ p[mp] ][i].first, epos = nd[ p[mp] ][i].second;
                    x.clear( );
                    for ( int j = 0; j < (int) p.size( ); j++ )
                         x.push_back( p[j] );
                    ho_interval xpos( mp, mp+1 );
                    d = {e};
                    int dstart = epos, dstop = epos + 1;
                    ( !align2 ? Align : Align2 )( 
                         D, to_left, to_right, x, xpos, d, dstart, dstop );
     
                    // For now, toss improper alignments.

                    if ( xpos.Start( ) > 0 )
                    {    if ( dstart > 0 ) continue;
                         if ( dstart == 0 )
                         {    int v = to_left[ d.front( ) ];
                              if ( D.To(v).nonempty( ) ) continue;    }    }
                    if ( xpos.Stop( ) < x.isize( ) ) 
                    {    if ( dstop < D.O( d.back( ) ).isize( ) ) continue;
                         int w = to_right[ d.back( ) ];
                         if ( D.From(w).nonempty( ) ) continue;    }

                    // Compute start position of read on d[0].

                    pos = p.getOffset( );
                    for ( int j = 0; j < dstart; j++ )
                         pos += hb.Kmers( D.O(d[0])[j] );

                    // Give up if more than one placement.

                    nplaces++;
                    if ( nplaces > 1 ) break;

                    // Save.  
                    // Just saving d, but could save ( xpos, d, dstart, dstop ).

                    aplace.clear( );
                    for ( int j = 0; j < d.isize( ); j++ )
                         aplace.push_back( d[j] );    }

               // For now require that there is a unique proper placement.

               if ( nplaces != 1 ) continue;

               // Check for broken placement.  Need to track down how this could
               // ever happen.

               int n = hb.K( ) - 1;
               for ( auto f : aplace ) n += dlens[f];
               if ( pos > n ) continue;

               // Save placement.

               dpaths[id].resize( aplace.size( ) );
               for ( int j = 0; j < aplace.isize( ); j++ )
                    dpaths[id][j] = aplace[j];
               dpaths[id].setOffset(pos);    }     }

     if (verbose) cout << Date( ) << ": done, used " << TimeSince(clock) << endl;
     int64_t placed = 0;
     for ( int64_t id = 0; id < (int64_t) paths.size( ); id++ )
          if ( dpaths[id].size( ) > 0 ) placed++;
     if (verbose)
     {    cout << Date( ) << ": " 
               << PERCENT_RATIO( 3, placed, (int64_t) paths.size( ) ) 
               << " placed" << endl;    }    }

void PlaceReadsSmart( const HyperBasevectorX& hb, const ReadPathVecX& paths, 
     const vec<Bool>& dup, const digraphE<vec<int>>& D, const vec<int>& dinv,
     ReadPathVec& dpaths, const vec<vec<vec<vec<int>>>>& dlines, 
     const vec<int64_t>& bci, const Bool verbose, const int btest,
     const Bool align2 )
{
     // Set up.

     vec<int> dlens( D.E( ), 0 ), to_left, to_right, linv;
     #pragma omp parallel for
     for ( int e = 0; e < D.E( ); e++ )
     {    if ( D.O(e)[0] < 0 ) continue;
          for ( int j = 0; j < D.O(e).isize( ); j++ )
               dlens[e] += hb.Kmers( D.O(e)[j] );    }
     D.ToLeft(to_left), D.ToRight(to_right);
     LineInv( dlines, dinv, linv );

     // To each nongap edge, assign a line, and a starting position on that line.
     // For nonpalindromic lines L, we pick just one of L, rc(L), and place both
     // edges d and dinv[d] at the same place on that line.

     vec<quad<int,int,Bool,int>> tol( D.E( ), make_quad(-1,-1,True,-1) );
     cout << Date( ) << ": defining tol" << endl;
     for ( int i = 0; i < dlines.isize( ); i++ )
     {    if ( linv[i] < i ) continue;
          const vec<vec<vec<int>>>& L = dlines[i];
          int pos = 0;
          for ( int j = 0; j < L.isize( ); j++ )
          {    vec<int> lens;
               for ( int k = 0; k < L[j].isize( ); k++ )
               {    int len = 0;
                    for ( int l = 0; l < L[j][k].isize( ); l++ )
                    {    int d = L[j][k][l];
                         if ( D.O(d)[0] < 0 ) continue;
                         if ( tol[d].first < 0 )
                         {    tol[d] = make_quad( i, pos+len, True, j );
                              tol[ dinv[d] ] 
                                   = make_quad( i, pos+len, False, j );    }
                         len += dlens[d];    }
                    lens.push_back(len);    }
               Sort(lens);
               if ( lens.nonempty( ) ) pos += Median(lens);    }    }

     // Define the territory of each barcode.  The granularity of this definition
     // is weak, as we use only start/stop positions (i,j) on lines, whereas we
     // should use base positions.

     // Go through the barcodes.

     cout << Date( ) << ": start traversal" << endl;
     const int MAX_BC_GAP = 100000;
     const int MIN_BC_GROUP = 3;
     int ndots = 0, stopped = 0;
     #pragma omp parallel for schedule(dynamic, 1000)
     for ( int b = 1; b < bci.isize( ) - 1; b++ )
     {    if ( btest >= 0 && b != btest ) continue;

          // Logging.

          ostringstream out;
          if (verbose)
          {    out << "\nlooking at barcode " << b << " [" << bci[b+1]-bci[b]
                    << " reads]" << endl;    }

          // Define the territory of each barcode.  The granularity of this
          // definition is weak, as we use only start/stop positions (i,j) on lines,
          // whereas we should use base positions.

          int total = 0;
          vec<triple<int,int,int>> lplaces; // (line, start on line, line unit)
          for ( int64_t id = bci[b]; id < bci[b+1]; id++ )
          {    const ReadPath& p = dpaths[id];
               int pos = p.getOffset( );
               for ( int j = 0; j < (int) p.size( ); j++ )
               {    if ( j == 0 || tol[ p[j] ].first != tol[ p[j-1] ].first )
                    {    triple<int,int,int> x; 
                         x.first = tol[ p[j] ].first;
                         int xpos = pos;
                         if ( !tol[ p[j] ].third ) xpos = dlens[ p[j] ] - xpos;
                         x.second = tol[ p[j] ].second + xpos;
                         x.third = tol[ p[j] ].fourth;
                         lplaces.push_back(x);    }
                    pos -= dlens[ p[j] ];    }    }
          Sort(lplaces);
          vec<triple<int,int,int>> groups;
          for ( int i = 0; i < lplaces.isize( ); i++ )
          {    int j;
               for ( j = i + 1; j < lplaces.isize( ); j++ )
               {    if ( lplaces[j].first != lplaces[i].first ) break;
                    if ( lplaces[j].second - lplaces[j-1].second > MAX_BC_GAP ) 
                         break;    }
               if ( j - i >= MIN_BC_GROUP )
               {    int len = lplaces[j-1].second - lplaces[i].second;
                    int start = lplaces[i].third, stop = lplaces[j-1].third;
                    int l = lplaces[i].first;
                    const vec<vec<vec<int>>>& L = dlines[l];

                    // Extend start and stop positions.

                    const int EXT = 25000;
                    int left = 0, right = 0;
                    while( left < EXT && start >= 1 )
                    {    start--;
                         vec<int> lens;
                         for ( int k = 0; k < L[start].isize( ); k++ )
                         {    int len = 0;
                              for ( int l = 0; l < L[start][k].isize( ); l++ )
                              {    int d = L[start][k][l];
                                   len += dlens[d];    }
                              lens.push_back(len);    }
                         Sort(lens);
                         if ( lens.nonempty( ) ) left += Median(lens);    }
                    while( right < EXT && stop < L.isize( ) - 1 )
                    {    stop++;
                         vec<int> lens;
                         for ( int k = 0; k < L[stop].isize( ); k++ )
                         {    int len = 0;
                              for ( int l = 0; l < L[stop][k].isize( ); l++ )
                              {    int d = L[stop][k][l];
                                   len += dlens[d];    }
                              lens.push_back(len);    }
                         Sort(lens);
                         if ( lens.nonempty( ) ) right += Median(lens);    }
                    
                    // Save.

                    if (verbose)
                    {    out << l << "." << lplaces[i].second << "-" 
                              << lplaces[j-1].second << " [len=" << len 
                              << ",count=" << j-i 
                              << ",start=" << start << "/" << L[start][0][0]
                              << ",stop=" << stop << "/" << L[stop][0][0]
                              << "]" << endl;    }
                    total += len;
                    groups.push( lplaces[i].first, start, stop );    }
               i = j - 1;    }    
          if (verbose) out << "total = " << ToStringAddCommas(total) << endl;

          // Define the edges in the territory.

          vec<pair<int,int>> locs;
          int nd = 0;
          for ( int i = 0; i < groups.isize( ); i++ )
          {    const vec<vec<vec<int>>>& L = dlines[ groups[i].first ];
               int start = groups[i].second, stop = groups[i].third;
               for ( int j = start; j <= stop; j++ )
               {    const vec<vec<int>>& M = L[j];
                    vec<int> ds;
                    for ( int k = 0; k < M.isize( ); k++ )
                    for ( int l = 0; l < M[k].isize( ); l++ )
                    {    int d = M[k][l];
                         if ( D.O(d)[0] >= 0 ) ds.push_back( d, dinv[d] );    }
                    UniqueSort(ds);
                    nd += ds.size( );
                    for ( int k = 0; k < ds.isize( ); k++ )
                    {    int d = ds[k];
                         vec<int> es = D.O(d);
                         UniqueSort(es);
                         for ( auto e : es ) locs.push( e, d );    }    }    }
          Sort(locs);
          if (verbose) 
          {    out << "locs = " << ToStringAddCommas( locs.size( ) ) 
                    << " base edges in " << ToStringAddCommas(nd) 
                    << " superedges" << endl;    }

          // Try to place more reads.

          int pl1 = 0, pl2 = 0;
          for ( int64_t id = bci[b]; id < bci[b+1]; id++ )
          {    if ( dpaths[id].size( ) > 0 ) pl1++;
               ReadPath p;
               paths.unzip(p,hb,id);
               if ( p.size( ) == 0 || dpaths[id].size( ) > 0 || dup[id/2] ) continue;

               // Now p = paths[id] represents a nonduplicate read from the barcode
               // that is placed on the base graph.  Find supergraph edges assigned 
               // to the barcode (as defined by locs), and find the matches of the 
               // read to one those edges that extends for at least 30 kmers; also 
               // ignore matches that are less than half as long as the best match.
               // (Not sure the "30" is enforced below.)
               //
               // First let pl =
               // { (supergraph edge d, pos on p, pos on d, match length in kmers) }.

               vec<quad<int,int,int,int>> pl;
               for ( int j = 0; j < (int) p.size( ); j++ )
               {    int e = p[j];
                    int low = LowerBound1(locs, e), high = UpperBound1(locs, e);
                    for ( int m = low; m < high; m++ )
                    {    int d = locs[m].second;
                         for ( int l = 0; l < D.O(d).isize( ); l++ )
                         {    if ( D.O(d)[l] != e ) continue;
                              if ( l > 0 && j > 0 && D.O(d)[l-1] == p[j-1] ) 
                                   continue;
                              int n = 0, r;
                              for ( r = l; r < D.O(d).isize( ); r++ )
                              {    if ( r - l + j >= (int) p.size( ) ) break;
                                   if ( D.O(d)[r] != p[r-l+j] ) break;
                                   n += hb.Kmers( D.O(d)[r] );    }
                              if ( btest >= 0 ) PRINT5_TO( out, id, d, j, l, n );
                              pl.push( d, j, l, n );    }    }    } 

               // Now filter the matches as described above.

               int nmax = 0;
               for ( int i = 0; i < pl.isize( ); i++ )
                    nmax = Max( nmax, pl[i].fourth );
               vec<Bool> to_delete( pl.size( ), False );
               for ( int i = 0; i < pl.isize( ); i++ )
               {    int n = pl[i].fourth;
                    if ( nmax >= 30 && n < nmax/2 ) to_delete[i] = True;    }
               EraseIf( pl, to_delete );
               if ( btest >= 0 ) PRINT2_TO( out, id, pl.size( ) );

               // If there is a unique placement, place the read.

               if ( pl.solo( ) )
               {    int d = pl[0].first, j = pl[0].second, l = pl[0].third;
                    vec<int> x;
                    for ( auto e : p ) x.push_back(e);
                    ho_interval xpos( j, j + 1 );
                    vec<int> q = {d};
                    int qstart = l, qstop = l + 1;
                    ( !align2 ? Align : Align2 )( 
                         D, to_left, to_right, x, xpos, q, qstart, qstop );

                    // For now, toss improper alignments.

                    if ( btest >= 0 )
                    {    PRINT6_TO( out, id, xpos.Start( ), qstart, xpos.Stop( ),
                              x.size( ), qstop );    }
                    Bool bad = False;
                    if ( xpos.Start( ) > 0 )
                    {    if ( qstart > 0 ) continue;
                         if ( qstart == 0 )
                         {    int v = to_left[ q.front( ) ];
                              if ( D.To(v).nonempty( ) ) continue;    }    }
                    if ( xpos.Stop( ) < x.isize( ) )
                    {    if ( qstop < D.O( q.back( ) ).isize( ) ) continue;
                         int w = to_right[ q.back( ) ];
                         if ( D.From(w).nonempty( ) ) continue;    }

                    // Compute new offset.

                    int pos = p.getOffset( );
                    for ( int j = 0; j < qstart; j++ )
                         pos += hb.Kmers( D.O(q[0])[j] );

                    // Check for broken placements.  Not sure how this can happen.

                    int n = hb.K( ) - 1;
                    for ( auto f : q ) n += dlens[f];
                    if ( pos > n ) continue;

                    // Save the new placement.

                    dpaths[id].resize( q.size( ) );
                    for ( int j = 0; j < q.isize( ); j++ )
                         dpaths[id][j] = q[j];
                    // out << "placing read " << id << " at pos " << pos
                    //      << ", path = " << printSeq(q) << endl;
                    dpaths[id].setOffset(pos);    
                    pl2++;    }    }
          if (verbose) 
          {    out << "orig reads placed = " << pl1 << ", new reads placed = "
                    << pl2 << "\n";
               #pragma omp critical
               {    cout << out.str( );    }    }
          else
          {
               #pragma omp critical
               {    MakeDots( stopped, ndots, bci.isize( ) - 2 );    }    }    }    }
