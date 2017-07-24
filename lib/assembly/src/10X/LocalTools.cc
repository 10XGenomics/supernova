// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"
#include "ParallelVecUtilities.h"
#include "paths/HyperBasevector.h"
#include "paths/long/large/Lines.h"
#include "10X/Heuristics.h"
#include "10X/LocalTools.h"
#include "10X/PlaceReads.h"
#include "10X/Super.h"

void GlueAssemblies( const HyperBasevectorX& hb, digraphE<vec<int>>& D, 
     vec<int>& dinv, const vec<int> mult,
     const int MIN_LINE_TO_WALK, const int MAX_MULT, const Bool verbose )
{
     // Go through multiple join steps.

     for ( int pass = 1; ; pass++ )
     {
          // Gather long edges in the graph.

          if (verbose) cout << "\n" << Date( ) << ": pass = " << pass << endl;
          if (verbose) cout << Date( ) << ": getting edge lengths" << endl;
          vec<int> dlens( D.E( ), 0 ), dmm( D.E( ), 1000000000 );
          // for each super-edge d in D
          #pragma omp parallel for
          for ( int d = 0; d < D.E( ); d++ )
          {    
               // don't deal with bad d
               if ( D.O(d)[0] < 0 ) continue;
               // for each basegraph edge, get kmer size and minimum mutliplicity of an edge in d -- i.e. support
               for ( int j = 0; j < D.O(d).isize( ); j++ )
               {    dlens[d] += hb.Kmers( D.O(d)[j] );
                    dmm[d] = Min( dmm[d], mult[ D.O(d)[j] ] );    }    }

          if (verbose) cout << Date( ) << ": processing long edges" << endl;
          vec<triple<int,int,int>> xlocs;
          // for each edge in D
          for ( int d = 0; d < D.E( ); d++ )
          {    
               // skip short edge
               if ( dlens[d] < MIN_LINE_TO_WALK ) continue;
               // skip edges that have too many dups
               if ( dmm[d] > MAX_MULT ) continue;
               // for each edge e in d, store position if multiplicity > 1
               for ( int j = 0; j < D.O(d).isize( ); j++ )
               {    int e = D.O(d)[j];
                    if ( mult[e] > 1 ) xlocs.push( e, d, j );    }    }
          if (verbose) cout << Date( ) << ": sorting" << endl;
          // sort so that instances related to an edge appear consecutively
          ParallelSort(xlocs);
          vec<int64_t> xstarts( hb.E( ) + 1, -1 );
          xstarts[0] = 0;
          // store cumulative index for each new edge e
          for ( int64_t j = xlocs.jsize( ) - 1; j >= 0; j-- )
               xstarts[ xlocs[j].first ] = j;
          xstarts[ hb.E( ) ] = xlocs.size( );
          for ( int64_t j = hb.E( ) - 1; j >= 0; j-- )
               if ( xstarts[j] < 0 ) xstarts[j] = xstarts[j+1];

          // Find long matches.  This finds instances of two edges (d1, d2) that 
          // share an identical stretch of >= MIN_LINE_TO_WALK kmers.  Once we find 
          // an instance for an edge d1, we move on to the next d1.

          if (verbose) cout << Date( ) << ": looking for long matches" << endl;
          vec< triple< pair<int,int>, pair<int,int>, int > > M;
          // for each superedge d1 in D
          #pragma omp parallel for schedule(dynamic,10000)
          for ( int d1 = 0; d1 < D.E( ); d1++ )
          {    
               // skip if not long enough
               if ( dlens[d1] < MIN_LINE_TO_WALK ) continue;
               // slip if reps
               if ( dmm[d1] > MAX_MULT ) continue;
               // for each edge e in d1
               for ( int j1 = 0; j1 < D.O(d1).isize( ); j1++ )
               {    int e = D.O(d1)[j1], jx, s0 = 0;
                    // for edges after e in d1
                    for ( jx = j1; jx < D.O(d1).isize( ); jx++ )
                    {    
                         // accure kmer size past e
                         s0 += hb.Kmers( D.O(d1)[jx] );
                         // stop when kmers long
                         if ( s0 >= MIN_LINE_TO_WALK ) break;    }
                    // if past e we don't find long seq, skip d1
                    if ( s0 < MIN_LINE_TO_WALK ) break;
                    // otherwise if d1 passes these tests,
                    // for each instance of e in various superedges of D
                    int64_t low = xstarts[e], high = xstarts[e+1];
                    for ( int64_t i = low; i < high; i++ )
                    {    int d2 = xlocs[i].second, j2 = xlocs[i].third;
                         // skip d1
                         if ( d2 == d1 && j2 == j1 ) continue;
                         // if not start to the match, skip
                         if ( j1 > 0 && j2 > 0 && D.O(d1)[j1-1] == D.O(d2)[j2-1] ) 
                              continue;
                         // if not enough stuff to match, skip
                         if ( j2 - j1 + jx >= D.O(d2).isize( ) ) continue;
                         // get rcs
                         int rd1 = dinv[d1], rd2 = dinv[d2];
                         // skip chimeras?
                         if ( !IsUnique( vec<int>{ d1, d2, rd1, rd2 } ) ) continue;
                         Bool fail = False;
                         // scan down from jx to j1 and match corresponding part of d2
                         for ( int l1 = jx; l1 > j1; l1-- )
                         {    int l2 = j2 - j1 + l1;
                              if ( D.O(d1)[l1] != D.O(d2)[l2] )
                              {    fail = True;
                                   break;    }    }
                         // if match failed, skip d2
                         if (fail) continue;
                         // otherwise match found
                         int s = 0, l1;
                         // get size of match
                         for ( l1 = jx + 1; l1 < D.O(d1).isize( ); l1++ )
                         {    int l2 = j2 - j1 + l1;
                              if ( l2 >= D.O(d2).isize( ) ) break;
                              if ( D.O(d1)[l1] != D.O(d2)[l2] ) break;
                              s += hb.Kmers( D.O(d1)[l1] );    }
                         int min_mult = 1000000000;
                         // find min mult in region of match and skip if too big
                         for ( int m = j1; m < l1; m++ )
                              min_mult = Min( min_mult, mult[ D.O(d1)[m] ] );
                         if ( min_mult > MAX_MULT ) continue;
                         // somehow this is symmetric wrt involution
                         #pragma omp critical
                         {    M.push( 
                                   make_pair(d1,j1), make_pair(d2,j2), l1-j1 );    }
                         goto next_d1;    }    }
                    next_d1: continue;    }
          if (verbose) cout << Date( ) << ": sorting matches" << endl;
          ParallelSort(M);

          // Make merges.

          if ( MakeMerges( M, D, dinv, verbose ) == 0 ) break;    }    }

// Find general two-in two-out constructs that can be pulled apart.

void BigPull( const HyperBasevectorX& hb, const vec<int>& inv, 
     digraphE<vec<int>>& D, vec<int>& dinv, const ReadPathVec& paths, 
     const vec<Bool>& dup, const Bool verbose )
{
     // Heuristics.

     const int MIN_LONG_LINE = 1000;
     const int MAX_DIST = 400;
     const int MAX_DEPTH = 4;
     const int MIN_RATIO = 5;
     const int MAX_KILL = 4;

     // Build ancillary data.

     vec<int> to_left, to_right, llens, dlens( D.E( ), 0 );
     D.ToLeft(to_left), D.ToRight(to_right);
     vec<vec<vec<vec<int>>>> dlines;
     FindLines( D, dinv, dlines, MAX_CELL_PATHS, MAX_CELL_DEPTH );
     GetLineLengths( hb, D, dlines, llens );
     #pragma omp parallel for
     for ( int d = 0; d < D.E( ); d++ )
     {    if ( D.O(d)[0] < 0 ) continue;
          for ( int j = 0; j < D.O(d).isize( ); j++ )
               dlens[d] += hb.Kmers( D.O(d)[j] );    }
     digraphE<int> G( D, dlens );
     vec<int> tol1( D.E( ), -1 ), tol2( D.E( ), -1 );
     for ( int i = 0; i < dlines.isize( ); i++ )
     {    if ( llens[i] >= MIN_LONG_LINE )
          {    const vec<vec<vec<int>>>& L = dlines[i];
               tol1[ L.front( )[0][0] ] = i; 
               tol2[ L.back( )[0][0] ] = i;    }    }
     ReadPathVec dpaths;
     PlaceReads( hb, paths, dup, D, dpaths, verbose, False );
     IntIndex dpaths_index( dpaths, D.E( ) );

     // Find candidates.

     if (verbose) cout << Date( ) << ": finding candidates" << endl;
     vec< pair< vec<int>, vec<int> > > cans;
     for ( int l1 = 0; l1 < dlines.isize( ); l1++ )
     {    if ( llens[l1] < MIN_LONG_LINE ) continue;
          const vec<vec<vec<int>>>& L1 = dlines[l1];
          int v1 = to_right[ L1.back( )[0][0] ];

          // Find vertices near v1.

          vec<pair<int,int>> nears;
          nears.push( v1, 0 );
          for ( int p = 1; p <= MAX_DEPTH; p++ )
          {    vec<pair<int,int>> nears2 = nears;
               for ( int i = 0; i < nears.isize( ); i++ )
               {    int v = nears[i].first, n = nears[i].second;
                    for ( int j = 0; j < D.From(v).isize( ); j++ )
                    {    int v2 = D.From(v)[j];
                         int n2 = n + dlens[ D.IFrom(v,j) ];
                         if ( n2 <= MAX_DIST )
                         {    Bool found = False;
                              for ( int l = 0; l < nears2.isize( ); l++ )
                              {    if ( nears2[l].first != v2 ) continue;
                                   found = True;
                                   nears2[l].second = Min( nears2[l].second, n2 );
                                   break;    }
                              if ( !found ) nears2.push( v2, n2 );    }    }    
                    for ( int j = 0; j < D.To(v).isize( ); j++ )
                    {    int v2 = D.To(v)[j];
                         int n2 = n + dlens[ D.ITo(v,j) ];
                         if ( n2 <= MAX_DIST )
                         {    Bool found = False;
                              for ( int l = 0; l < nears2.isize( ); l++ )
                              {    if ( nears2[l].first != v2 ) continue;
                                   found = True;
                                   nears2[l].second = Min( nears2[l].second, n2 );
                                   break;    }
                              if ( !found ) nears2.push( v2, n2 );    }    }    }
               nears = nears2;    }

          // Find long lines that start or stop at one of those vertices.

          vec<int> lefts, rights;
          for ( int i = 0; i < nears.isize( ); i++ )
          {    int v = nears[i].first;
               for ( int j = 0; j < D.From(v).isize( ); j++ )
               {    int d = D.IFrom(v,j);
                    int l = tol1[d];
                    if ( l >= 0 ) rights.push_back(l);    }
               for ( int j = 0; j < D.To(v).isize( ); j++ )
               {    int d = D.ITo(v,j);
                    int l = tol2[d];
                    if ( l >= 0 ) lefts.push_back(l);    }    }
          UniqueSort(lefts), UniqueSort(rights);
          if ( lefts.size( ) == 2 && rights.size( ) == 2 )
               cans.push( lefts, rights );    }
     UniqueSort(cans);

     // Track touching of the left and right vertices of an edge.

     vec<Bool> touchedl( D.E( ), False ), touchedr( D.E( ), False );

     // Define edits.  There are two passes.  On the first pass we define
     // the edits.  On the second pass we make them.

     if (verbose) cout << Date( ) << ": finding support" << endl;
     int count = 0;
     vec<Bool> to_edit( cans.size( ), False ), swapped( cans.size( ), False );
     vec< vec<vec<vec<int>>> > PATHS( cans.size( ) );
     ostringstream out;
     for ( int zpass = 1; zpass <= 2; zpass++ )
     for ( int i = 0; i < cans.isize( ); i++ )
     {    if ( zpass == 2 && !to_edit[i] ) continue;
          int L1 = cans[i].first[0], L2 = cans[i].first[1];
          int R1 = cans[i].second[0], R2 = cans[i].second[1];
          int d1 = dlines[L1].back( )[0][0], d2 = dlines[L2].back( )[0][0];
          int f1 = dlines[R1].front( )[0][0], f2 = dlines[R2].front( )[0][0];
          int rd1 = dinv[d1], rd2 = dinv[d2], rf1 = dinv[f1], rf2 = dinv[f2];
          int v1 = to_right[d1], v2 = to_right[d2];
          int w1 = to_left[f1], w2 = to_left[f2];
          if ( zpass == 2 )
          {    if ( touchedr[d1] || touchedr[d2] || touchedl[f1] || touchedl[f2] )
                    continue;    }
          vec<int64_t> s11, s12, s21, s22;
          vec< pair< int, vec<int> > > X;
          vec<int> ds = {d1,d2};

          // Compute support.

          Bool solved = False;
          if ( zpass == 1 )
          {    for ( int pass = 1; pass <= 2; pass++ )
               for ( int di = 0; di < ds.isize( ); di++ )
               {    int d = ds[di];
                    if ( pass == 2 ) d = dinv[d];
                    for ( int j = 0; j < dpaths_index.Count(d); j++ )
                    {    int64_t id1 = dpaths_index.Val( d, j );
                         int64_t id2 = ( id1 % 2 == 0 ? id1 + 1 : id1 - 1 );
                         const ReadPath &p1 = dpaths[id1], &p2 = dpaths[id2];
                         vec<int> x1, x2;
                         for ( auto d : p1 ) x1.push_back(d);
                         for ( auto d : p2 ) x2.push_back( dinv[d] );
                         x2.ReverseMe( );
                         vec<int> x(x1);
                         x.append(x2);
                         if ( pass == 2 )
                         {    x.ReverseMe( );
                              for ( auto& g : x ) g = dinv[g];    }
                         X.push( id1/2, x );    }    }
               for ( int i = 0; i < X.isize( ); i++ )
               {    int p1 = Position( X[i].second, d1 ); 
                    int p2 = Position( X[i].second, f1 );
                    if ( p1 >= 0 && p1 < p2 ) s11.push_back( X[i].first );
                    p2 = Position( X[i].second, f2 );
                    if ( p1 >= 0 && p1 < p2 ) s12.push_back( X[i].first );
                    p1 = Position( X[i].second, d2 ); 
                    p2 = Position( X[i].second, f1 );
                    if ( p1 >= 0 && p1 < p2 ) s21.push_back( X[i].first );
                    p2 = Position( X[i].second, f2 );
                    if ( p1 >= 0 && p1 < p2 ) s22.push_back( X[i].first );    }
               UniqueSort(s11), UniqueSort(s12); 
               UniqueSort(s21), UniqueSort(s22);

               // Can be solved by support.

               if ( s12.isize( ) <= MAX_KILL && s21.isize( ) <= MAX_KILL
                    && s11.isize( ) >= MIN_RATIO * Max( 1, s12.isize( ) )
                    && s22.isize( ) >= MIN_RATIO * Max( 1, s21.isize( ) ) )
               {    solved = True;    }
               if ( s11.isize( ) <= MAX_KILL && s22.isize( ) <= MAX_KILL
                    && s12.isize( ) >= MIN_RATIO * Max( 1, s11.isize( ) )
                    && s21.isize( ) >= MIN_RATIO * Max( 1, s22.isize( ) ) )
               {    solved = True;
                    swapped[i] = True;    }

               // Can be solved by inversion.
               
               if ( !solved )
               {    vec<int> m1 = {d1,d2}, m2 = {dinv[f1],dinv[f2]};
                    Sort(m1), Sort(m2);
                    if ( m1 == m2 )
                    {    if ( d1 == dinv[f1] ) swapped[i] = True;
                         solved = True;    }    }    }

          if ( swapped[i] )
          {    swap(R1,R2), swap(s11,s12), swap(s21,s22); 
               swap(f1,f2), swap(rf1,rf2), swap(w1,w2);    }

          // Find the paths across.  If we hit a gap edge, the computation
          // may blow up, and we give up.  Should be handled differently.

          vec<vec<vec<int>>> paths(2);
          if ( zpass == 1 )
          {    
               // Print.

               out << "\n[" << ++count << "] ";
               out << "L1 = " << L1 << ", L2 = " << L2 << ", R1 = " << R1 
                    << ", R2 = " << R2 << " (" << s11.size( ) << "," 
                    << s12.size( ) << "," << s21.size( ) << "," 
                    << s22.size( ) << ")" << endl;
               if ( !solved )
               {    out << "Not solved." << endl;
                    continue;    }

               Bool fail = False;
               for ( int p = 0; p < 2; p++ )
               {    int v = ( p == 0 ? v1 : v2 ), w = ( p == 0 ? w1 : w2 );
                    const int MAX_PATHS = 10;
                    const int MAX_PARTIALS = 100;
                    Bool ok = G.AllPathsLengthRangeAlt( v, w, 0, MAX_DIST, 
                         to_right, paths[p], MAX_PATHS, 0, False, False, 
                         MAX_PARTIALS );
                    if ( p == 0 ) out << "paths from L1 to R1:\n";
                    else out << "paths from L2 to R2:\n";
                    if ( !ok )
                    {    out << "Attempt to find all paths failed." << endl;
                         fail = True;    }
                    else if ( paths[p].empty( ) )
                    {    out << "Found no paths." << endl;
                         fail = True;    }
                    else
                    {    for ( int j = 0; j < paths[p].isize( ); j++ )
                         {    out << "[" << j+1 << "] " 
                                   << printSeq( paths[p][j] ) << endl;    }    }
                    if (fail) break;    }
               if (fail) continue;

               // Explicitly test for gap edge and give up if there is one.

               Bool gap = False;
               for ( int p = 0; p < 2; p++ )
               for ( int j = 0; j < paths[p].isize( ); j++ )
               for ( int l = 0; l < paths[p][j].isize( ); l++ )
                    if ( D.O( paths[p][j][l] )[0] < 0 ) gap = True;
               if (gap) 
               {    out << "Encountered gap edge." << endl;
                    continue;    }

               // Don't allow empty edges.  We could allow them, but the code
               // below would have to be modified to avoid creating empty edges 
               // in D.
     
               Bool empty = False;
               for ( int p = 0; p < 2; p++ )
               for ( int j = 0; j < paths[p].isize( ); j++ )
                    if ( paths[p][j].empty( ) ) empty = True;
               if (empty) 
               {    out << "Found empty path." << endl;
                    continue;    }

               // OK we're good.

               to_edit[i] = True;
               PATHS[i] = paths;

               continue;    }
     
          // Now edit the graph.

          paths = PATHS[i];
          int rv1 = to_right[rf1], rv2 = to_right[rf2];
          int rw1 = to_left[rd1], rw2 = to_left[rd2];
          touchedr[d1] = touchedr[d2] = touchedl[f1] = touchedl[f2] = True;
          touchedr[rf1] = touchedr[rf2] = touchedl[rd1] = touchedl[rd2] = True;
          vec<int> verts = {v1,v2,w1,w2,rv1,rv2,rw1,rw2};
          UniqueSort(verts);
          for ( auto x : verts ) D.SplayVertexWithUpdate(x, to_left, to_right);
          vec< triple<int,int,vec<int>> > adds;
          map< triple<int,int,vec<int>>, triple<int,int,vec<int>> > rc;
          for ( int p = 0; p < 2; p++ )
          for ( int j = 0; j < paths[p].isize( ); j++ )
          {    int v = ( p == 0 ? to_right[d1] : to_right[d2] ); 
               int w = ( p == 0 ? to_left[f1] : to_left[f2] );
               vec<int> x;
               for ( auto d : paths[p][j] ) x.append( D.O(d) );
               adds.push( v, w, x );
               int rv = ( p == 0 ? to_right[rf1] : to_right[rf2] ); 
               int rw = ( p == 0 ? to_left[rd1] : to_left[rd2] );
               vec<int> rx = x;
               rx.ReverseMe( );
               for ( int l = 0; l < rx.isize( ); l++ ) rx[l] = inv[ rx[l] ];
               adds.push( rv, rw, rx );
               rc[ make_triple( v, w, x ) ] = make_triple( rv, rw, rx );
               rc[ make_triple( rv, rw, rx ) ] = make_triple( v, w, x );    }
          UniqueSort(adds);
          int E = D.E( );
          for ( int j = 0; j < adds.isize( ); j++ )
          {    int v = adds[j].first, w = adds[j].second;
               const vec<int>& y = adds[j].third;
               D.AddEdgeWithUpdate( v, w, y, to_left, to_right );
               int j2;
               for ( j2 = 0; j2 < adds.isize( ); j2++ )
                    if ( adds[j2] == rc[ adds[j] ] ) break;
               dinv.push_back( E + j2 );    }    }
     if (verbose) cout << out.str( );    }

// Find subset of read paths that should be mapped to the final union of assemblies
// grab reads that used to map to the relevant region
// and reads from barcodes used for local assembly
void FindRelevantPaths ( digraphE<vec<int>> & D, vec<int> & dinv, ReadPathVecX & pathsx,
     VecULongVec & dpaths_index, vec<int64_t> & bci, 
     vec <Bool> & keep, vec < vec<int> > & bs, vec <Bool> & pmask )
{
     // by default all read-pairs are going to be skipped
     ForceAssertGt( pathsx.size(), 0 );
     pmask.resize (pathsx.size()/2, True );
     
     // grab read ids that would have been placed on super edges
     // that we are keeping around
     {
          for ( int d = 0; d < D.E( ); d++ ) {
               if ( keep[d] ) {
                    for ( auto id : dpaths_index[d] ) {
                         pmask[id/2] = False;
                    }
                    for ( auto id : dpaths_index[dinv[d]] ) {
                         pmask[id/2] = False;
                    }
               }
          }
     }
     // all barcodes that are in bs
     for ( auto & bcs : bs ) {
          // add in barcodes from bs
          for ( auto & b : bcs ) {
               for ( int64_t id = bci[b]; id < bci[b+1]; id+=2 ) {
                    pmask[id/2]=False;
               }
          }    
     }
}
