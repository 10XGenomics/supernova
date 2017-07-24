// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

// Merge along shorter overlaps but nearby in the graph as defined by heuristics
// LOOK, SEE, LOOK_MERGE

#include "10X/mergers/ShortMergers.h"

// Given two super edge objects x1 and x2 look for matching sub-sequences of
// base vectors that have a length >= MIN_OVERLAP. If allow_two then allow for up to
// two such matches.
// return True if a unique match exists and store that information in overs.
// Otherwise, return False

inline Bool FindOverlap( const HyperBasevectorX & hb, 
     const vec<int> & x1, const vec<int> & x2, 
     vec<triple<int,int,int>> & overs, const int & MIN_OVERLAP, const Bool & allow_two )
{
     vec< triple<int,int,int> > com;
     for ( int i = 0; i < x1.isize( ); i++ )
          com.push( x1[i], i, 1 );
     for ( int i = 0; i < x2.isize( ); i++ )
          com.push( x2[i], i, 2 );
     Sort(com);
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
                    if ( n < MIN_OVERLAP ) continue;
                    int z = l1 - p1;
                    overs.push( p1, p2, z );
                    if ( overs.size( ) > 2 ) return False; }    }
          r = s - 1;    }
     if ( overs.size( ) == 0 ) return False;
     if ( overs.size( ) == 2 && !allow_two ) return False;
     if ( overs.size( ) == 2 && overs[1].third > overs[0].third ) 
          swap( overs[0], overs[1] );
     return True;   }

// Explore the graph D towards the right (using From) from a vertex v 
// up to a depth=DEPTH and store the resulting edges that are encountered in
// de (depth, edge). Note that DEPTH <= 255 (uint8_t). 

void ExploreRightToDepth( const digraphE<vec<int>> & D, const vec<int> & to_right,
     const vec<Bool> & used, const int & v, vec<pair<uint8_t, int>> & de,
     const uint8_t DEPTH )
{
     de.clear( );
     
     set<int> all = {v};
     set<int> vset = {v};
     for ( uint8_t depth = 0; depth < DEPTH; depth++ ) {
          const int start = de.size( );
          for ( auto & vp : vset ) {
               for ( int i = 0; i != D.From(vp).isize(); i++ ) {
                    const int edge = D.IFrom( vp, i );
                    if ( used[edge] ) de.push( depth, edge );
               }
          }
          vset.clear( );
          for ( int i = start; i != de.isize( ); i++ ) {
               const int w = to_right[de[i].second]; 
               if ( all.count( w ) ) continue;
               vset.insert( w );
               all.insert( w );
          }
     }
     Sort( de );
}

// Support function for FindZipperCandidates( ) that implements the 
// recursive traversal

void FindBaseEdgeOverlap( const digraphE<vec<int>> & D, const vec<int> & dinv,
     const vec<int> & to_left, const vec<int> & to_right,     
     const int & v, const int & voffset, const int & w, const int & woffset,
     set<int> & vvisited, vec<nptuple> & M, Bool & overflow, int & depth,
     const int & MAX_DEPTH )
{
     depth++;
     if ( depth >= MAX_DEPTH ) return;

     // Find edges emanating from v and w
     // that are not gaps or abutting seq gap edges with non zero trim
     vec<int> ve, we;
     for ( int pass = 1; pass <= 2; pass++ ) {
          if ( pass == 2 && v == w ) {
               we = ve;
               break;
          }
          const int vert = ( pass == 1 ? v : w );
          vec<int> & es = ( pass == 1 ? ve : we );
          for ( int i = 0; i != D.From(vert).isize(); i++ ) {
               const int d = D.IFrom(vert,i);
               if ( D.O(d)[0] < 0 ) continue;
               const Bool seqgap = DetectSeqGapsWithTrim( D, to_left, to_right, d );
               if ( seqgap ) continue;
               es.push_back( d );
          }
     }
     
     if ( ve.empty( ) || we.empty( ) ) return;

     vec <quad<int,int,int,int>> next;
     for ( auto & d1 : ve ) {
          for ( auto & d2 : we ) {
               if ( d1 == d2 ) continue;
               if ( v == w && d1 > d2 ) continue;
               const int rv = to_right[d1], rw = to_right[d2];
               
               // don't zipper if the edge is a loop
               if ( rv == v || rw == w ) continue;

               const int rd1 = dinv[d1], rd2 = dinv[d2];        
               // are we in awkward configurations with the inv edge?
               if ( !IsUnique( rd1, rd2, d1, d2 ) )
                    continue;
               
               const vec<int> & x1 = D.O(d1);
               const vec<int> & x2 = D.O(d2);
               const int n1 = x1.size(), n2 = x2.size();
               int p1 = voffset, p2 = woffset;
               while ( p1 < n1 && p2 < n2 ) {
                    if ( x1[p1] != x2[p2] ) break;
                    p1++;
                    p2++;
               }
               int len = p1 - voffset;
               if ( len == 0 ) continue;
               
               // don't create loops
               // and avoid potential pathologies with shared vertices
               
               const int rrdv = to_right[rd1];
               if ( to_right[d1] == to_right[d2] ) {
                    if (!IsUnique(vec<int> {v, to_right[d1], rrdv, to_left[rd1]}))
                         continue;
                    if ( (p1 == n1 && p2 == n2-1) || (p1 == n1-1 && p2 == n2) ) {
                         p1--;
                         p2--;
                         len--;
                    }
                    if ( len == 0 ) continue;
               } else if (!IsUnique(vec<int>  {v, to_right[d1], to_right[d2],
                         rrdv, to_left[rd1], to_left[rd2] }))
                    continue;
 
               if ( len >= numeric_limits<uint16_t>::max() ) {
                    overflow=True;
                    continue;
               }
               // we can zipper
               
               M.push( make_pair( d1, voffset ),
                    make_pair( d2, woffset ), len );
               M.push( make_pair( dinv[d1], n1-len-voffset ),
                    make_pair( dinv[d2], n2-len-woffset ), len );
               
               // if both matches stop short of end
               if ( p1 < n1 && p2 < n2 ) continue;
               
               const int rvseen = vvisited.count( rv );
               const int rwseen = vvisited.count( rw );
               // otherwise go on zippering
               if ( p1 == n1 && p2 < n2 && rvseen == 0 )
                    next.push( rv, 0, w, p2 );
               if ( p1 < n1 && p2 == n2 && rwseen == 0)
                    next.push( v, p1, rw, 0 );
               if ( p1 == n1 && p2 == n2 && rvseen == 0 && rwseen == 0)
                    next.push( rv, 0, rw, 0 );
         }
     }
     UniqueSort( next );
     
     for ( int i = 0; i != next.isize( ); i++ ) {
          const auto & n1 = next[i];
          if ( n1.second == 0 )
               vvisited.insert( n1.first );
          if ( n1.fourth == 0 )
               vvisited.insert( n1.third );
          for ( int j = 0; j != next.isize( ); j++ ) {
               const auto & n2 = next[j];
               FindBaseEdgeOverlap( D, dinv, to_left, to_right,
                    n1.first, n1.second, n2.third, n2.fourth,
                    vvisited, M, overflow, depth, MAX_DEPTH );
          }
     }
}

// Find candidates for zippering by starting at a vertex and exploring
// along zipper lines recursively.

void FindZipperCandidates( const digraphE<vec<int>> & D, const vec<int> & dinv,
     const vec<int> & to_left, const vec<int> & to_right,
     vec<nptuple> & M, const int MAX_DEPTH )
{
     Bool detect_overflow = False;
     int done = 0, prev = 0;
     #pragma omp parallel for schedule (dynamic, 1) reduction(+:done)
     for ( int v = 0; v < D.N( ); v++ ) {
          set<int> vvisited = {v};
          int off1 = 0, off2 = 0;
          int v1 = v, v2 = v;
          vec<nptuple> Mpiece;
          Bool overflow = False;
          int depth = 0;
          FindBaseEdgeOverlap( D, dinv, to_left, to_right, v1, off1, v2, off2,
               vvisited, Mpiece, overflow, depth, MAX_DEPTH );
          done++;
          if ( Mpiece.empty( ) ) continue;
          UniqueSort( Mpiece );
          #pragma omp critical
          {
               M.append( Mpiece );
               if ( !detect_overflow && overflow )
                    detect_overflow=True;
          }
     }
     if ( detect_overflow )
          cout << "WARNING: skipped paths too long for uint16_t" << endl;
}

// Given a list of matches (d1, p1), (d2, p2), overlap, verify that they
// are really matches in the super graph

void VerifyMatchList( const digraphE<vec<int>> & D, const vec<nptuple> & M )
{
     #pragma omp parallel for
     for ( size_t i = 0; i < M.size( ); i++ ) {
          auto & m = M[i];
          const int d1 = m.first.first, p1 = m.first.second;
          const int d2 = m.second.first, p2 = m.second.second;
          const int n = m.third;
          const vec<int> & x1 = D.O(d1), & x2 = D.O(d2);
          Bool mismatch = False;
          for ( int j = 0; j < n; j++ ) {
               if ( x1[p1+j] != x2[p2+j] )
                    mismatch = True;
          }
          if ( mismatch ) {
               #pragma omp critical
               {
                    cout << "Invalid match: ";
                    PRINT5(d1, p1, d2, p2, n);
                    cout << "D.O(d1): ";
                    for ( auto & e : x1 )
                         cout << e << " ";
                    cout << endl << "D.O(d2): ";
                    for ( auto & e : x2 )
                         cout << e << " ";
                    cout << endl;
               }
          }
     }
}


// Try to merge along relatively short overlaps that are proximate in the graph.
// This is a low memory version of MergeShortOverlaps

void MergeShortOverlapsLowMem( const HyperBasevectorX& hb, const vec<int>& inv, 
     digraphE<vec<int>>& D, vec<int>& dinv, const int LOOK_MERGE, const int LOOK, 
     const Bool allow_two, const Bool verbose, const Bool debug )
{
     if (verbose) cout << Date( ) << ": start micromerger" << endl;
     vec<int> to_left, to_right, lens( D.E( ), 0 );
     D.ToLeft(to_left), D.ToRight(to_right);
     vec<Bool> used;
     D.Used( used );
     #pragma omp parallel for
     for ( int e = 0; e < D.E( ); e++ )
     {    if ( D.O(e)[0] < 0 ) continue;
          for ( int j = 0; j < D.O(e).isize( ); j++ )
               lens[e] += hb.Kmers( D.O(e)[j] );    }
     if (verbose) cout << Date( ) << ": and begin the loop" << endl;
     
     // this is a triple <pair,pair,int> of ints
     // named a nested penta-tuple
     vec <nptuple> M;
     const int SEE=Min(3, LOOK);
     #pragma omp parallel for schedule(dynamic, 10000)
     for ( int v = 0; v < D.N( ); v++ ) {
          vec<pair<uint8_t, int>> de;
          ExploreRightToDepth( D, to_right, used, v, de, uint8_t(LOOK) );
          for ( int i = 0; i != de.isize(); i++ ) {
               int j = i;
               set<int> see;
               for ( ; j != de.isize(); j++ ) {
                    if ( de[j].first - de[i].first >= SEE ) break;
                    const int d = de[j].second;
                    if ( lens[d] >= LOOK_MERGE )
                         see.insert( d );
               }

               for ( auto it1 = see.begin(); it1 != see.end(); it1++ ) {
                    const int d1 = *it1;

                    // skip out if there is a seq gap edge with non-trivial trim
                    Bool seqgap = DetectSeqGapsWithTrim( D, to_left, to_right, d1 );
                    if ( seqgap ) continue;

                    for ( auto it2 = it1; it2 != see.end(); it2++ ) {
                    const int d2 = *it2;
                    if ( d1 == d2 ) continue;
                    
                    Bool seqgap = DetectSeqGapsWithTrim( D, to_left, to_right, d2 );
                    if ( seqgap ) continue;
                    // Look for an overlap of at least LOOK_SEE kmers between
                    // edges d1 and d2.  Give up if there is more than one.

                    const vec<int> &x1 = D.O(d1), &x2 = D.O(d2);
                    vec< triple<int,int,int> > overs;
                    const Bool found = FindOverlap( hb, x1, x2, overs, LOOK_MERGE,
                         allow_two ); 
                    if ( ! found ) continue;

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
               i = j-1;
          }
          next_v_candidate: continue;    }
     
     if (verbose)
          cout << Date( ) << ": sorting merges; peak mem = " << PeakMemUsageGBString( )
               << ", mem = "<<MemUsageGBString( )<< endl;
     ParallelUniqueSort(M);
     
     if ( debug ) {
          cout << Date( ) << ": validating matches" << endl;
          VerifyMatchList( D, M );
     }

     if (verbose){ cout << Date( ) << ": done; peak mem = " << PeakMemUsageGBString( )
         << ", mem = "<<MemUsageGBString( )<< endl;}
     if ( verbose ) cout << Date( ) << ": found " << ToStringAddCommas( M.size( ) ) 
          << " candidates with overlap" << endl;
     
     // Make the merges using C2G algorithm
     String debugdir = "";
     digraphE<vec<int>> D_new;
     vec<int> dinv_new;
     cout << Date( ) << ": initiate gluing" << endl;
     Bool GET_META = False;
     vec<vec<quad<int,uint16_t,int,uint16_t>>> dummy_paths;
     GlueGraphs( D, dinv, M, hb.E(), D_new, dinv_new, 
             GET_META, dummy_paths, verbose, debugdir );
     #pragma omp parallel sections
     {
     #pragma omp section
     {    D    = D_new;  }
     #pragma omp section
     {    dinv = dinv_new; }
     }
     if (verbose) 
          cout << Date( ) << ": done. " << " merges; peak mem = " 
               << PeakMemUsageGBString( ) << ", mem = " << MemUsageGBString( )<< endl;
}

// Find long zippering matches recursively up to MAX_DEPTH > 0
// and repeat the process MAX_PASSES times, or until the number
// of zips is below MIN_ZIPS. MIN_ZIPS = -1 turns this off. 
// In each iteration,
// we form a large vector of pair-wise super edge gluings
// and call GlueGraph to perform the gluings.

void ZipperRecursive( const HyperBasevectorX & hb, const vec<int> & inv, 
     digraphE<vec<int>> & D, vec<int> & dinv, const String DEBUG_DIR,
     const int MAX_DEPTH, const int MAX_PASSES, const int MIN_ZIPS,
     const Bool verbose, const Bool gg_verbose )
{
     const Bool GG=True;
     cout << Date( ) << ": making up to " << MAX_PASSES 
          << " passes to zipper" << endl;
          
     int zips = 0;
     int passes = 0;
     const Bool debug = False;
     // this is a triple <pair,pair,int> of ints
     // named a nested penta-tuple
     vec <nptuple> M;
     vec<int> to_left, to_right;
     vec<Bool> used;
     const int nthreads = omp_get_max_threads( );
     while ( passes < MAX_PASSES ) {
          passes++;
          if ( verbose ) cout << Date( ) << ": pass " << passes << endl;
          to_left.clear();
          to_right.clear();
          used.clear();
          M.clear( );
          D.ToLeft(to_left), D.ToRight(to_right);
          D.Used( used );
          
          if ( verbose ) cout << Date( ) << ": finding edge pairs to zip" 
                              << endl;

          FindZipperCandidates( D, dinv, to_left, to_right, M, MAX_DEPTH );
          
          if ( verbose ) cout << Date( ) << ": unique sorting matches" << endl;
          if ( M.empty() )
               break;
          ParallelUniqueSort(M);
          
          if ( MIN_ZIPS >= 0 && M.isize( ) < MIN_ZIPS )
               break;

          if ( debug ) {
               cout << Date( ) << ": writing data" << endl;
               VerifyMatchList( D, M );
               String base = DEBUG_DIR + "/zip_" + ToString( passes );
               Mkdir777( base );
               ofstream out( base + "/matches.out" );
               for ( auto & m : M )
                    out  << m.first.first << "," << m.first.second << " "
                         << m.second.first << "," << m.second.second << " "
                         << m.third << endl;
               cout << Date( ) << ": validate" << endl;
               Validate( hb, inv, D, dinv );
               BinaryWriter::writeFile( base + "/a.sup", D );
               BinaryWriter::writeFile( base + "/a.sup.inv", dinv );
          }
          if ( verbose ) cout << Date( ) << ": found " << ToStringAddCommas(M.size( ))
                              << " zips" << endl;
          if ( GG ) {
               // Make the merges using C2G algorithm
               String debugdir = "";
               digraphE<vec<int>> D_new;
               vec<int> dinv_new;
               if ( verbose )
                    cout << Date ( ) << ": using GG algorithm" << endl;
               zips += M.size();
               Bool GET_META = False;
               vec<vec<quad<int,uint16_t,int,uint16_t>>> dummy_paths;
               GlueGraphs<int,uint16_t>( D, dinv, M, hb.E(), D_new, dinv_new,
                    GET_META, dummy_paths, gg_verbose, debugdir );
               
               #pragma omp parallel sections
               {
                    #pragma omp section
                    {    D    = D_new;  }
                    #pragma omp section
                    {    dinv = dinv_new; }
               }
               if ( verbose ) cout << Date( ) << ": done with graph gluing" << endl;
          } else {
               if ( verbose )
                    cout << Date( ) << ": using MM algorithm" << endl;
               int joins = MakeMergesOnly<int>( M, D, dinv, verbose );
               if ( verbose ) 
                    cout << Date( ) << ": made " << joins << " joins" << endl;
               zips += joins;
          }
         
          if ( verbose ) 
               cout << Date( ) << ": graph checksum "
                    << ToStringAddCommas(D.CheckSum( )) << endl;
          else {
               if ( passes % 50 == 0 )
                    cout << endl;
               else if ( passes % 10 == 0 )
                    cout << " ";
               cout << ".";
               flush(cout);
          }

     }
     RemoveUnneededVertices( D, dinv );
     CleanupCore( D, dinv );
     
     cout << endl << Date( ) << ": made " << ToStringAddCommas(zips) << " zips in "
          << passes << " passes" << endl;
     cout << Date( ) << ": graph checksum "
          << ToStringAddCommas(D.CheckSum( )) << endl;
}

// Zipper large graphs by first calling ZipperRecursive followed by the old
// Zipper() function that zippers a graph iteratively. The combination is
// supposed to optimize performance. When the number of zips is large, 
// ZipperRecursive is fast, but doesn't do it completely, and once the number
// of zips is tractable use the old Zipper function

void ZipperCombo( const HyperBasevectorX & hb, const vec<int> & inv, 
     digraphE<vec<int>> & D, vec<int> & dinv, const String DEBUG_DIR,
     const Bool verbose )
{
     const int MAX_DEPTH = 50;
     const int MAX_PASSES = 100;
     const int MIN_ZIPS = 1000;
     cout << Date( ) << ": begin recursive zippering" << endl;
     ZipperRecursive( hb, inv, D, dinv, DEBUG_DIR, MAX_DEPTH, MAX_PASSES,
          MIN_ZIPS, verbose, verbose );
     cout << Date( ) << ": begin iterative zippering" << endl;
     const Bool SINGLE=False;
     Zipper( D, dinv, SINGLE, True );
     cout << Date( ) << ": final graph checksum " << ToStringAddCommas( D.CheckSum( ) )
          << endl;
     cout << Date( ) << ": done; peak = " << PeakMemUsageGBString( )
          << ", mem = " << MemUsageGBString( )<< endl;
}
