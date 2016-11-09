// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>

#include "CoreTools.h"
#include "Set.h"
#include "graph/DigraphTemplate.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "paths/long/large/Lines.h"
#include "10X/Capture.h"
#include "10X/DfTools.h"
#include "10X/Heuristics.h"
#include "10X/IntIndex.h"
#include "10X/Scaffold.h"
#include "10X/Super.h"

void Scaffold( const HyperBasevectorX& hb, const vec<int>& inv,
     const VecIntVec& ebcx, digraphE<vec<int>>& D, vec<int>& dinv, 
     const ReadPathVec& dpaths, const vec<int64_t>& bid,
     const vec<DataSet>& datasets, const Bool verbose, String& link_report,
     const Bool single )
{
     vec<Bool> used;
     D.Used(used);
     if (verbose)
          cout << Date( ) << ": entering Scaffold, sanity check input" << endl;
     #pragma omp parallel for
     for ( int64_t i = 0; i < (int64_t) dpaths.size( ); i++ )
     {    const ReadPath& p = dpaths[i];
          for ( int j = 0; j < (int) p.size( ); j++ )
          {    if ( p[j] >= D.E( ) )
               {    cout << "\nFound path that references out of range edge id."
                         << endl;
                    PRINT2( p[j], D.E( ) );
                    TracebackThisProcess( );
                    Scram(1);    }
               if ( !used[ p[j] ] )
               {    cout << "\nFound path through unused edge." << endl;
                    TracebackThisProcess( );
                    Scram(1);    }    }    }

     // Find lines, and track positions of edges on them.

     if (verbose) cout << Date( ) << ": finding lines" << endl;
     vec<vec<vec<vec<int>>>> lines;
     FindLines( D, dinv, lines, MAX_CELL_PATHS, MAX_CELL_DEPTH );
     // { (line_id, base_pos_on_line) }:
     vec< pair<int,int> > tol( D.E( ), make_pair(-1,-1) ); 
     vec<int> lens( lines.size( ), 0 );
     for ( int i = 0; i < lines.isize( ); i++ )
     {    int pos = 0;
          const vec<vec<vec<int>>>& L = lines[i];
          for ( int j = 0; j < L.isize( ); j++ )
          {    vec<int> lensj;
               for ( int k = 0; k < L[j].isize( ); k++ )
               {    int len = 0;
                    for ( int l = 0; l < L[j][k].isize( ); l++ )
                    {    // Note overwriting, not necessarily good:
                         tol[ L[j][k][l] ] = make_pair( i, pos + len );    
                         int d = L[j][k][l];
                         if ( D.O(d)[0] < 0 ) continue;
                         for ( int m = 0; m < D.O(d).isize( ); m++ )
                              len += hb.Kmers( D.O(d)[m] );    }
                    lensj.push_back(len);    }
               Sort(lensj);
               if ( lensj.nonempty( ) ) pos += Median(lensj);    }
          lens[i] = pos + hb.K( ) - 1;    }

     // Sanity check of lines.

     {    vec<Bool> lused( D.E( ), False );
          for ( int i = 0; i < lines.isize( ); i++ )
          for ( int j = 0; j < lines[i].isize( ); j++ )
          for ( int k = 0; k < lines[i][j].isize( ); k++ )
          for ( int l = 0; l < lines[i][j][k].isize( ); l++ )
               lused[ lines[i][j][k][l] ] = True;
          for ( int d = 0; d < D.E( ); d++ )
          {    if ( !used[d] || lused[d] ) continue;
               if ( D.O(d)[0] < 0 ) continue;
               cout << "\nEdge d = " << d << " is used but not in a line." << endl;
               cout << "d = " << printSeq( D.O(d) ) << endl << endl;
               Scram(1);    }    }

     // Find barcodes on lines, excluding an end.

     const int END_DEL = 1000;
     const int MIN_REMAINING = 1000;
     if (verbose) cout << Date( ) << ": finding multiplicities" << endl;
     vec<int> mult;
     ComputeMult( hb, D, mult );
     if (verbose) cout << Date( ) << ": finding barcodes on lines" << endl;
     vec<int> left_enuf( lines.size( ), False );
     vec<int> right_enuf( lines.size( ), False );
     vec<vec<int>> left_bc( lines.size( ) ), right_bc( lines.size( ) );
     vec<vec<int>> buckets;
     BucketLines( lines, lens, buckets );
     #pragma omp parallel for schedule(dynamic, 1)
     for ( int b = 0; b < buckets.isize( ); b++ )
     {    vec<int> x, y, lensj;
          for ( int bi = 0; bi < buckets[b].isize( ); bi++ )
          {    int i = buckets[b][bi];
               if ( lens[i] < END_DEL + MIN_REMAINING ) continue;
               const vec<vec<vec<int>>>& L = lines[i];
               x.clear( ), y.clear( );
               for ( int j = 0; j < L.isize( ); j += 2 )
               for ( int m = 0; m < L[j][0].isize( ); m++ )
               {    int d = L[j][0][m];
                    if ( D.O(d)[0] >= 0 ) x.append( D.O(d) );    }
               for ( int j = 0; j < x.isize( ); j++ )
                    y.push_back( hb.Kmers( x[j] ) );
               int sum = 0;
               for ( int j = 0; j < y.isize( ); j++ )
               {    sum += y[j];
                    if ( sum >= END_DEL )
                    {    if ( Sum(y) - sum >= MIN_REMAINING ) 
                              left_enuf[i] = True;    }    }
               sum = 0;
               for ( int j = y.isize( ) - 1; j >= 0; j-- )
               {    sum += y[j];
                    if ( sum >= END_DEL )
                    {    if ( Sum(y) - sum >= MIN_REMAINING ) 
                              right_enuf[i] = True;    }    }
               int pos = 0;
               for ( int j = 0; j < L.isize( ); j++ )
               {    lensj.clear( );
                    for ( int k = 0; k < L[j].isize( ); k++ )
                    {    int len = 0;
                         for ( int l = 0; l < L[j][k].isize( ); l++ )
                         {    int d = L[j][k][l];
                              if ( D.O(d)[0] < 0 ) continue;
                              for ( int m = 0; m < D.O(d).isize( ); m++ )
                              {    int e = D.O(d)[m];
                                   if ( pos + len >= END_DEL && mult[e] == 1 )
                                        left_bc[i].append( ebcx[e] );
                                   len += hb.Kmers( D.O(d)[m] );    }    }
                         lensj.push_back(len);    }
                    Sort(lensj);
                    if ( lensj.nonempty( ) ) pos += Median(lensj);    }
               pos = 0;
               for ( int j = L.isize( ) - 1; j >= 0; j-- )
               {    lensj.clear( );
                    for ( int k = 0; k < L[j].isize( ); k++ )
                    {    int len = 0;
                         for ( int l = L[j][k].isize( ) - 1; l >= 0; l-- )
                         {    int d = L[j][k][l];
                              if ( D.O(d)[0] < 0 ) continue;
                              for ( int m = D.O(d).isize( ) - 1; m >= 0; m-- )
                              {    int e = D.O(d)[m];
                                   if ( pos + len >= END_DEL && mult[e] == 1 )
                                        right_bc[i].append( ebcx[e] );
                                   len += hb.Kmers( D.O(d)[m] );    }    }
                         lensj.push_back(len);    }
                    Sort(lensj);
                    if ( lensj.nonempty( ) ) pos += Median(lensj);    }
               UniqueSort( left_bc[i] ), UniqueSort( right_bc[i] );    }    }

     // Translate paths into line space.

     if (verbose) cout << Date( ) << ": originating lpaths" << endl;
     ReadPathVec lpaths(dpaths); // TERRIBLE WAY TO RESERVE!!!!!!!!!!!!!!!!!!!!!!!!!
     if (verbose) cout << Date( ) << ": translating paths into line space" << endl;
     int nthreads = ( single ? 1 : omp_get_max_threads( ) );
     #pragma omp parallel for num_threads(nthreads)
     for ( int64_t id = 0; id < (int64_t) dpaths.size( ); id++ )
     {    lpaths[id].resize(0);
          const ReadPath& p = dpaths[id];
          if ( p.size( ) == 0 ) continue;
          int offset = p.getOffset( );
          int d = p[0];
          ForceAssertGe( d, 0 );
          ForceAssertGe( D.O(d)[0], 0 );
          int l = tol[d].first, lpos = tol[d].second;
          ForceAssertGe( l, 0 );
          lpaths[id].setOffset( offset + lpos );
          for ( int j = 0; j < (int) p.size( ); j++ )
          {    int l = tol[ p[j] ].first;
               if ( j == 0 || l != lpaths[id].back( ) )
                    lpaths[id].push_back(l);    }    }

     // Index line paths.

     if (verbose) cout << Date( ) << ": indexing line paths" << endl;
     IntIndex lpaths_index( lpaths, lines.size( ), verbose );

     // Find links.

     const int MAX_SEP = 1000;
     const int MIN_BC = 5;
     vec<vec< pair<int,int64_t> >> linksto( D.E( ) );
     vec<int> to_left, to_right;
     D.ToLeft(to_left), D.ToRight(to_right);
     if (verbose) cout << Date( ) << ": start link loop" << endl;
     #pragma omp parallel for num_threads(nthreads)
     for ( int l1 = 0; l1 < lines.isize( ); l1++ )
     {    for ( int i = 0; i < lpaths_index.Count(l1); i++ )
          {    
               // Get distance of first read start from right end of line l1.

               const int64_t id1 = lpaths_index.Val( l1, i );
               const ReadPath& p1 = lpaths[id1];
               int dist1 = lens[l1] - p1.getOffset( );

               // Get distance of second read start from right end of line l2.

               const int64_t id2 = ( id1 % 2 == 0 ? id1+1 : id1-1 );
               const ReadPath& p2 = lpaths[id2];
               if ( p2.size( ) == 0 ) continue;
               int l2 = p2[0];
               int dist2 = lens[l2] - p2.getOffset( );

               // Require that implied fragment size is not too large.

               if ( dist1 + dist2 > MAX_SEP ) continue;

               // Don't link to self or inverse.

               int d1 = lines[l1].back( )[0][0];
               int d2 = dinv[ lines[l2].back( )[0][0] ];
               if ( d2 == d1 || d2 == dinv[d1] ) continue;

               // Let's only go sink to source.

               int v = to_right[d1], w = to_left[d2];
               if ( !D.Sink(v) || !D.Source(w) ) continue;
               if ( D.To(v).size( ) > 1 || D.From(w).size( ) > 1 ) continue;

               // Require barcode linking if there's an opportunity.

               if ( right_enuf[l1] && left_enuf[l2]
                    && MeetSize( right_bc[l1], left_bc[l2] ) < MIN_BC )
               {    continue;    }

               // Record link.

               linksto[d1].push( d2, id1 );    }    }

     // Record links.

     vec<vec<int>> linksto2( D.E( ) );
     ostringstream lrout2;
     {    for ( int d = 0; d < D.E( ); d++ )
          {    Sort( linksto[d] );
               vec< triple< int, int, vec<int64_t> > > where;
               for ( int i = 0; i < linksto[d].isize( ); i++ )
               {    int j;
                    for ( j = i + 1; j < linksto[d].isize( ); j++ )
                         if ( linksto[d][j].first != linksto[d][i].first ) break;
                    vec<int64_t> pids;
                    vec<int> bids;
                    for ( int k = i; k < j; k++ )
                    {    pids.push_back( linksto[d][k].second/2 );
                         bids.push_back( bid[ linksto[d][k].second ] );    }
                    UniqueSort(bids);
                    if ( bids.size( ) > 1 )
                    {    linksto2[d].push_back( linksto[d][i].first );
                         where.push( linksto[d][i].first, bids.size( ), pids );    }
                    i = j - 1;    }
               if ( where.nonempty( ) )
               {    lrout2 << "\nreaching from " << d 
                         << " (inv=" << dinv[d] << "):" << endl;
                    for ( int j = 0; j < where.isize( ); j++ )
                    {    lrout2 << "[" << j+1 << "] " << where[j].first
                                << " (x" << where[j].third.size( ) << ",b" 
                                << where[j].second << ")";
                         if ( where[j].third.size( ) <= 6 )
                              lrout2 << "  pids=" << printSeq( where[j].third );
                         lrout2 << endl;    }    }    }    }
     link_report = lrout2.str( );

     // Now make the gaps.  This is the first point where we edit the graph.
     // The only action taken is to add edges.

     if (verbose) cout << Date( ) << ": making gaps" << endl;
     set< pair<int,int> > gaps;
     int ne2 = D.E( );
     for ( int d1 = 0; d1 < ne2; d1++ )
     {    if ( !linksto2[d1].solo( ) ) continue;
          int d2 = linksto2[d1][0];
          int rd1 = dinv[d1], rd2 = dinv[d2];
          if ( !linksto2[rd2].solo( ) ) continue;
          if ( linksto2[rd2][0] != rd1 ) continue;
          if ( make_pair(rd2,rd1) < make_pair(d1,d2) ) continue;

          // Assume no intertwining with inversion.

          if ( !IsUnique( d1, d2, rd1, rd2 ) ) continue;

          // If we haven't seen the gap before, insert it.

          if ( !Member( gaps, make_pair( to_right[d1], to_left[d2] ) ) )
          {    dinv.push_back( D.E( ) + 1, D.E( ) );
               D.AddEdge( to_right[d1], to_left[d2], {-1} );
               D.AddEdge( to_right[rd2], to_left[rd1], {-1} );
               gaps.insert( make_pair( to_right[d1], to_left[d2] ) );
               gaps.insert( make_pair( to_right[rd2], to_left[rd1] ) );    }    }
     if (verbose) cout << Date( ) << ": found " << gaps.size( ) << " gaps" << endl;

     // Close some gaps.  This is the second point where we edit the graph.  The
     // only action taken is local surgery on edges.

     if (verbose) cout << Date( ) << ": closing some gaps" << endl;
     int njoins = 0;
     // not actually sure this looping is needed
     while(1)
     {    D.ToLeft(to_left), D.ToRight(to_right);
          vec<Bool> touched( D.E( ), False );
          vec< pair<int,int> > joiners;
          #pragma omp parallel for schedule(dynamic, 10000), num_threads(nthreads)
          for ( int d = 0; d < D.E( ); d++ )
          {    if ( D.O(d)[0] != -1 ) continue;
               int rd = dinv[d];
               if ( rd <= d ) continue;
               if ( !D.To( to_left[d] ).solo( ) ) continue;
               if ( !D.From( to_right[d] ).solo( ) ) continue;
               int d1 = D.ITo( to_left[d], 0 ), d2 = D.IFrom( to_right[d], 0 );
               if ( d1 == d2 ) continue;
               if ( to_left[d1] == to_right[d2] ) continue;
               const vec<int> &x1 = D.O(d1), &x2 = D.O(d2);
               vec< pair<int,int> > overs;
               int count = 0;
               for ( int i1 = 0; i1 < x1.isize( ); i1++ )
               for ( int i2 = 0; i2 < x2.isize( ); i2++ )
               {    if ( x1[i1] != x2[i2] ) continue;
                    if ( i1 > 0 && i2 > 0 && x1[i1-1] == x2[i2-1] ) continue;
                    int n;
                    for ( n = 1; ; n++ )
                    {    if ( i1 + n == x1.isize( ) || i2 + n == x2.isize( ) ) break;
                         if ( x1[i1+n] != x2[i2+n] ) break;    }
                    if ( i1 + n != x1.isize( ) ) continue;
                    if ( i2 > 0 ) continue;
                    overs.push( i1, n );    }
               if ( overs.solo( ) )
               {
                    #pragma omp critical
                    {    joiners.push( d, overs[0].first );    }    }    }
          Bool made_join = False;
          Sort(joiners);
          for ( int jo = 0; jo < joiners.isize( ); jo++ )
          {    int d = joiners[jo].first, off = joiners[jo].second;
               int d1 = D.ITo( to_left[d], 0 ), d2 = D.IFrom( to_right[d], 0 );
               if ( touched[d] || touched[d1] || touched[d2] ) continue;
               // actually not sure the following lines do anything....
               if ( !D.To( to_left[d] ).solo( ) ) continue; // how possible?
               if ( !D.From( to_right[d] ).solo( ) ) continue; // how possible?
               made_join = True;
               njoins += 2;
               touched[d] = touched[d1] = touched[d2] = True;
               vec<int> y = D.O(d1);
               y.resize(off);
               y.append( D.O(d2) );
     
               int rd = dinv[d], rd1 = dinv[d1], rd2 = dinv[d2];
               touched[rd] = touched[rd1] = touched[rd2] = True;
               vec<int> ry = y;
               ry.ReverseMe( );
               for ( int i = 0; i < ry.isize( ); i++ )
                    ry[i] = inv[ ry[i] ];
     
               D.OMutable(d1) = y;
               D.DeleteEdgeFrom( to_right[d1], 0 );
               D.DeleteEdgeFrom( to_right[d], 0 );
               D.TransferEdges( to_right[d1], to_right[d2] );
               to_right[d1] = to_right[d2];
               D.OMutable(rd1) = ry;
               D.DeleteEdgeTo( to_left[rd1], 0 );
               D.DeleteEdgeTo( to_left[rd], 0 );
               D.TransferEdges( to_left[rd1], to_left[rd2] );
               to_left[rd1] = to_left[rd2];    }
          if ( !made_join ) break;    }
     if (verbose) cout << Date( ) << ": " << njoins << " gaps closed" << endl; 

     // Look for canonical and simple loops.

     vec<int> dels;
     CaptureCanonicalLoops( D, dinv, dels, verbose, single );
     CaptureSimpleLoops( D, dinv, dels, verbose, single );

     // Remove hanging ends.

     const int MAX_KILL = 320;
     const double MIN_RATIO = 25;
     SimpleHangs( hb, D, dinv, dels, MAX_KILL, MIN_RATIO, verbose, single );
     UniqueSort(dels);
     if (verbose) 
     {    cout << Date( ) << ": to delete " << dels.size( ) << " hanging ends" 
               << endl;    }

     // Delete edges.  This is the third point where we edit the graph.
     // Edges are deleted, then unneeded vertices are removed.  The latter adds
     // and deletes edges.

     if (single) D.DeleteEdges(dels);
     else D.DeleteEdgesParallel(dels);
     if (verbose) cout << Date( ) << ": del unneeded verts" << endl;
     RemoveUnneededVertices( D, dinv );    }
