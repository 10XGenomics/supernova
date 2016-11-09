///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"
#include "ParallelVecUtilities.h"
#include "Qualvector.h"
#include "graph/FindCells.h"
#include "kmers/KmerRecord.h"
#include "kmers/MakeLookup.h"
#include "math/Functions.h"
#include "paths/HyperBasevector.h"
#include "paths/RemodelGapTools.h"
#include "paths/long/KmerCount.h"
#include "paths/long/ReadPath.h"
#include "paths/long/large/GapToyTools.h"

// AnalyzeBranches: note not adjusting to_right.  This is wrong.

void AnalyzeBranches( HyperBasevector& hb, vec<int>& to_right, const vec<int>& inv2, 
     ReadPathVec& paths2, const Bool ANALYZE_BRANCHES_REV, 
     const int min_ratio2, const Bool ANALYZE_BRANCHES_VERBOSE )

{    double clock0 = WallClockTime( );
     vec<int> to_left;
     hb.ToLeft(to_left);
     for ( int i = 0; i < (int) paths2.size( ); i++ )
     {    ReadPath& p = paths2[i];
          for ( int j = 0; j < (int) p.size( ); j++ )
          {    if ( p[j] >= hb.EdgeObjectCount( ) ) p[j] = -1;
               if ( j > 0 && p[j-1] >= 0 && p[j] >= 0 
                    && to_right[ p[j-1] ] != to_left[ p[j] ] )
               {    p[j] = -1;    }    }    }

     // Heuristics.

     const int max_dist = 4;
     const int min_ratio = 5;
     const int max_kill = 2;

     vec< pair<int,int> > breaks;
     vec< vec<int> > froms( hb.EdgeObjectCount( ) ), tos( hb.EdgeObjectCount( ) );
     LogTime( clock0, "analyzing branches 0" );
     double clock1 = WallClockTime( );
     for ( int pass = 1; pass <= 2; pass++ )
     {    const int batch = 10000;
          int64_t npids = paths2.size( )/2;
          #pragma omp parallel for
          for ( int64_t bi = 0; bi < npids; bi += batch )
          {    vec< pair<int,int> > PP;
               for ( int64_t pid = bi; pid < Min( bi + batch, npids ); pid++ )
               {    vec<int> x, y;
                    for ( int j = 0; j < (int) paths2[2*pid].size( ); j++ )
                         x.push_back( paths2[2*pid][j] );
                    for ( int j = 0; j < (int) paths2[2*pid+1].size( ); j++ )
                         y.push_back( paths2[2*pid+1][j] );
                    y.ReverseMe( );
                    for ( int j = 0; j < y.isize( ); j++ )
                         if ( y[j] >= 0 ) y[j] = inv2[ y[j] ];
                    if ( pass == 2 )
                    {    swap( x, y );
                         x.ReverseMe( ), y.ReverseMe( );
                         for ( int j = 0; j < x.isize( ); j++ )
                              if ( x[j] >= 0 ) x[j] = inv2[ x[j] ];
                         for ( int j = 0; j < y.isize( ); j++ )
                              if ( y[j] >= 0 ) y[j] = inv2[ y[j] ];    }
                    pair< vec<int>, vec<int> > p = make_pair( x, y );
                    vec< pair<int,int> > P;
                    for ( int j1 = 0; j1 < p.first.isize( ) - 1; j1++ )
                    {    if ( p.first[j1] >= 0 && p.first[j1+1] >= 0 )
                              P.push( p.first[j1], p.first[j1+1] );    }
                    for ( int j1 = 0; j1 < p.second.isize( ) - 1; j1++ )
                    {    if ( p.second[j1] >= 0 && p.second[j1+1] >= 0 )
                              P.push( p.second[j1], p.second[j1+1] );    }
                    for ( int j1 = 0; j1 < p.first.isize( ); j1++ )
                    {    int x1 = p.first[j1];
                         if ( x1 >= 0 )
                         {    int m = Position( p.second, x1 );
                              if ( m < 0 && p.second.nonempty( ) 
                                   && p.second[0] >= 0 ) 
                              {    P.push( x1, p.second[0] );    }    }    }
                    UniqueSort(P);
                    PP.append(P);    }
               #pragma omp critical
               {    for ( int j = 0; j < PP.isize( ); j++ )
                    {    froms[ PP[j].first ].
                              push_back( PP[j].second );
                         tos[ PP[j].second ].
                              push_back( PP[j].first );    }    }    }    }
     LogTime( clock1, "analyzing branches 1" );
     double clock1b = WallClockTime( );
     #pragma omp parallel for
     for ( int e = 0; e < hb.EdgeObjectCount( ); e++ )
     {    Sort( froms[e] );
          Sort( tos[e] );    }
     if (ANALYZE_BRANCHES_VERBOSE) cout << "\nforward reach:\n";
     LogTime( clock1b, "analyzing branches 1b" );
     double clock2 = WallClockTime( );

     /*
     for ( int e = 0; e < hb.EdgeObjectCount( ); e++ )
     {    int v = to_right[e];
          if ( hb.From(v).size( ) != 2 || hb.To(v).size( ) > 1 ) continue;
          int n1 = 0, n2 = 0;
          for ( int j = 0; j < froms[e].isize( ); j++ )
          {    if ( froms[e][j] = hb.IFrom( v, 0 ) ) n1++;
               if ( froms[e][j] = hb.IFrom( v, 1 ) ) n2++;    }
          if ( n1 >= 10 && n2 <= 1 ) breaks.push( e, hb.IFrom( v, 1 ) );
          if ( n2 >= 10 && n1 <= 1 ) breaks.push( e, hb.IFrom( v, 0 ) );    }
     */

     for ( int e = 0; e < hb.EdgeObjectCount( ); e++ )
     {    int v = to_right[e];
          if ( hb.From(v).size( ) <= 1 ) continue;
          if ( hb.To(v).size( ) > 1 ) continue;
          vec< vec<int> > follow( hb.From(v).size( ) );
          vec<int> branches;
          for ( int j = 0; j < hb.From(v).isize( ); j++ )
          {    int f = hb.EdgeObjectIndexByIndexFrom( v, j );
               branches.push_back(f);    }
          int nbranches = branches.size( );

          for ( int j = 0; j < hb.From(v).isize( ); j++ )
          {    int f = hb.EdgeObjectIndexByIndexFrom( v, j );
               int w = to_right[f];
               for ( int l = 0; l < hb.From(w).isize( ); l++ )
                    follow[j].push_back( hb.EdgeObjectIndexByIndexFrom(w, l) );    }

          for ( int dpass = 1; dpass < max_dist; dpass++ )
          {    for ( int i = 0; i < nbranches; i++ )
               {    int n = follow[i].size( );
                    for ( int j = 0; j < n; j++ )
                    {    int w = to_right[ follow[i][j] ];
                         follow[i].append( hb.FromEdgeObj(w) );    }
                    UniqueSort( follow[i] );    }    }

          vec<int> fr, count;
          for ( int i = 0; i < froms[e].isize( ); i++ )
          {    int j = froms[e].NextDiff(i);
               int c = j - i, f = froms[e][i];
               fr.push_back(f), count.push_back(c);
               i = j - 1;    }
          vec<Bool> to_delete( fr.size( ), False );
          for ( int i = 0; i < fr.isize( ); i++ )
          {    vec<int> homes;
               for ( int j = 0; j < follow.isize( ); j++ )
                    if ( Member( follow[j], fr[i] ) ) homes.push_back(j);
               if ( homes.size( ) == follow.size( ) ) count[i] = 0;
               if ( homes.solo( ) )
               {    for ( int j = 0; j < fr.isize( ); j++ )
                    {    if ( fr[j] == hb.EdgeObjectIndexByIndexFrom( v, homes[0] ) )
                         {    count[j] += count[i];
                              count[i] = 0;    }    }    }    }
          for ( int i = 0; i < fr.isize( ); i++ )
               if ( count[i] == 0 ) to_delete[i] = True;
          EraseIf( fr, to_delete ), EraseIf( count, to_delete );
          vec<int> s1 = fr, s2 = branches;
          Sort(s1), Sort(s2);
          if ( s1 == s2 && s1.size( ) == 2 )
          {    if ( count[0] < min_ratio * count[1] 
                    && count[1] < min_ratio * count[0] )
               {    continue;    }    }
          ReverseSortSync( count, fr );
          if (ANALYZE_BRANCHES_VERBOSE)
          {    cout << e << " -->";
               for ( int i = 0; i < fr.isize( ); i++ )
                    cout << " " << fr[i] << "[" << count[i] << "]";    }
          if ( count.size( ) >= 2 && count[0] >= min_ratio2 * Max( 1, count[1] )
               && count[1] <= max_kill && Member( branches, fr[0] ) )
          {    if (ANALYZE_BRANCHES_VERBOSE) cout << " -- RECOMMEND PRUNING";
               for ( int j = 0; j < branches.isize( ); j++ )
                    if ( branches[j] != fr[0] ) breaks.push( e, branches[j] );    }
          if (ANALYZE_BRANCHES_VERBOSE) cout << "\n";    }
     UniqueSort(breaks);
     for ( int i = 0; i < breaks.isize( ); i++ )
     {    int e = breaks[i].first, f = breaks[i].second;
          int n = hb.N( );
          hb.AddVertices(2);
          hb.GiveEdgeNewFromVx( f, to_right[e], n );
          to_left[f] = n;
          int re = inv2[e], rf = inv2[f];
          if ( re >= 0 && rf >= 0 ) 
          {    hb.GiveEdgeNewToVx( rf, to_right[rf], n+1 );
               to_right[rf] = n+1;    }    }

     if (ANALYZE_BRANCHES_REV)
     {
     vec< pair<int,int> > breaksr;
     if (ANALYZE_BRANCHES_VERBOSE) cout << "\nbackward reach:\n";

     /*
     for ( int e = 0; e < hb.EdgeObjectCount( ); e++ )
     {    int v = to_left[e];
          if ( hb.To(v).size( ) != 2 || hb.From(v).size( ) > 1 ) continue;
          int n1 = 0, n2 = 0;
          for ( int j = 0; j < tos[e].isize( ); j++ )
          {    if ( tos[e][j] = hb.ITo( v, 0 ) ) n1++;
               if ( tos[e][j] = hb.ITo( v, 1 ) ) n2++;    }
          if ( n1 >= 10 && n2 <= 1 ) breaksr.push( hb.ITo( v, 1 ), e );
          if ( n2 >= 10 && n1 <= 1 ) breaksr.push( hb.ITo( v, 0 ), e );    }
     */

     for ( int e = 0; e < hb.EdgeObjectCount( ); e++ )
     {    int v = to_left[e];
          if ( hb.To(v).size( ) <= 1 ) continue;
          if ( hb.From(v).size( ) > 1 ) continue;
          vec< vec<int> > preceed( hb.To(v).size( ) );
          vec<int> branches;
          for ( int j = 0; j < hb.To(v).isize( ); j++ )
          {    int f = hb.EdgeObjectIndexByIndexTo( v, j );
               branches.push_back(f);    }
          int nbranches = branches.size( );

          for ( int j = 0; j < hb.To(v).isize( ); j++ )
          {    int f = hb.EdgeObjectIndexByIndexTo( v, j );
               int w = to_left[f];
               for ( int l = 0; l < hb.To(w).isize( ); l++ )
                    preceed[j].push_back( hb.EdgeObjectIndexByIndexTo(w, l) );    }

          for ( int dpass = 1; dpass < max_dist; dpass++ )
          {    for ( int i = 0; i < nbranches; i++ )
               {    int n = preceed[i].size( );
                    for ( int j = 0; j < n; j++ )
                    {    int w = to_left[ preceed[i][j] ];
                         preceed[i].append( hb.ToEdgeObj(w) );    }
                    UniqueSort( preceed[i] );    }    }

          vec<int> fr, count;
          for ( int i = 0; i < tos[e].isize( ); i++ )
          {    int j = tos[e].NextDiff(i);
               int c = j - i, f = tos[e][i];
               if ( to_right[f] == to_left[e] )
               {   fr.push_back(f), count.push_back(c);    }
               i = j - 1;    }
          vec<Bool> to_delete( fr.size( ), False );
          for ( int i = 0; i < fr.isize( ); i++ )
          {    vec<int> homes;
               for ( int j = 0; j < preceed.isize( ); j++ )
                    if ( Member( preceed[j], fr[i] ) ) homes.push_back(j);
               if ( homes.size( ) == preceed.size( ) ) count[i] = 0;
               if ( homes.solo( ) )
               {    for ( int j = 0; j < fr.isize( ); j++ )
                    {    if ( fr[j] == hb.EdgeObjectIndexByIndexTo( v, homes[0] ) )
                         {    count[j] += count[i];
                              count[i] = 0;    }    }    }    }
          for ( int i = 0; i < fr.isize( ); i++ )
               if ( count[i] == 0 ) to_delete[i] = True;
          EraseIf( fr, to_delete ), EraseIf( count, to_delete );
          vec<int> s1 = fr, s2 = branches;
          Sort(s1), Sort(s2);
          if ( s1 == s2 && s1.size( ) == 2 )
          {    if ( count[0] < min_ratio * count[1] 
                    && count[1] < min_ratio * count[0] )
               {    continue;    }    }
          ReverseSortSync( count, fr );
          if (ANALYZE_BRANCHES_VERBOSE)
          {    cout << e << " <--";
               for ( int i = 0; i < fr.isize( ); i++ )
                    cout << " " << fr[i] << "[" << count[i] << "]";    }
          if ( count.size( ) >= 2 && count[0] >= min_ratio2 * Max( 1, count[1] )
               && count[1] <= max_kill && Member( branches, fr[0] ) )
          {    if (ANALYZE_BRANCHES_VERBOSE) cout << " -- RECOMMEND PRUNING";
               for ( int j = 0; j < branches.isize( ); j++ )
                    if ( branches[j] != fr[0] ) breaksr.push( branches[j], e );    }
          if (ANALYZE_BRANCHES_VERBOSE) cout << "\n";    }
     UniqueSort(breaksr);
     for ( int i = 0; i < breaksr.isize( ); i++ )
     {    int e = breaksr[i].first, f = breaksr[i].second;
          int n = hb.N( );
          hb.AddVertices(2);
          hb.GiveEdgeNewToVx( e, to_left[f], n );
          to_right[e] = n;
          
          int re = inv2[e], rf = inv2[f];
          if ( re >= 0 && rf >= 0 ) 
          {    hb.GiveEdgeNewFromVx( re, to_left[re], n+1 );
               to_left[re] = n+1;    }    }
     breaks.append(breaksr);
     }

     int nb = breaks.size( );
     for ( int i = 0; i < nb; i++ )
          breaks.push( inv2[ breaks[i].second ], inv2[ breaks[i].first ] );
     UniqueSort(breaks);
     #pragma omp parallel for
     for ( int64_t i = 0; i < (int64_t) paths2.size( ); i++ )
     {    ReadPath& p = paths2[i];
          Bool bad = False;
          for ( int j = 0; j < ( (int) p.size( ) ) - 1; j++ )
          {    pair<int,int> x = make_pair( p[j], p[j+1] );
               if ( BinMember( breaks, x ) ) bad = True;    }
          if (bad) p.resize(0);    }
     if (ANALYZE_BRANCHES_VERBOSE) cout << "\n";
     LogTime( clock2, "analyzing branches 2" );    }

void ExtendTerminalEdges( const HyperBasevector& hb, 
     const vec< vec<int> >& layout_pos, const vec< vec<int64_t> >& layout_id, 
     const vec< vec<Bool> >& layout_or, const vecbasevector& bases, 
     const VecPQVec& quals )
{    
     vec<int> to_left, to_right;
     hb.ToLeft(to_left), hb.ToRight(to_right);
     qvec qv;
     for ( int e = 0; e < hb.EdgeObjectCount( ); e++ )
     {    if ( !hb.Sink( to_right[e] ) ) continue;
          vec< triple<int,char,uchar> > exts;
          cout << "extensions of edge " << e << endl;
          for ( int j = 0; j < layout_pos[e].isize( ); j++ )
          {    int id = layout_id[e][j];
               quals[id].unpack(&qv);
               int stop = layout_pos[e][j] + bases[id].isize( );
               if ( layout_pos[e][j] > hb.EdgeLengthBases(e) )
               {    cout << "weird" << endl;
                    continue;    }
               int ext = stop - hb.EdgeLengthBases(e);
               if ( ext <= 0 ) continue;
               if ( layout_or[e][j] )
               {    cout << "fw ";
                    for ( int l = 0; l < ext; l++ )
                    {    int pos = bases[id].isize( ) - ext + l;
                         // if ( quals[id][pos] <= 2 ) break;
                         cout << as_base( bases[id][pos] );
                         exts.push( l, bases[id][pos], qv[pos] );    }
                    cout << "\n";    }    
               else
               {    cout << "rc ";
                    for ( int l = 0; l < ext; l++ )
                    {    int pos = bases[id].isize( ) - 1
                              - ( bases[id].isize( ) - ext + l );
                         // if ( quals[id][pos] <= 2 ) break;
                         cout << as_base( 3 - bases[id][pos] );
                         exts.push( l, 3 - bases[id][pos], qv[pos] );   }
                    cout << "\n";    }    }    
          Sort(exts);
          for ( int i = 0; i < exts.isize( ); i++ )
          {    int j;
               for ( j = i + 1; j < exts.isize( ); j++ )
                    if ( exts[j].first != exts[i].first ) break;
               vec<int> qsum( 4, 0 ), ids( 4, vec<int>::IDENTITY );
               for ( int k = i; k < j; k++ )
                    qsum[ exts[k].second ] += exts[k].third;
               ReverseSortSync( qsum, ids );
               cout << exts[i].first << ":";
               for ( int l = 0; l < 4; l++ )
               {    if ( qsum[l] == 0 ) break;
                    cout << " " << as_base( ids[l] ) << "[" 
                         << qsum[l] << "]";    }
               cout << "\n";
               i = j - 1;    }    }    }

void PlaceMore( const HyperBasevector& hb, const vecbasevector& bases, 
		const VecPQVec& quals, ReadPathVec& paths2, vec<int64_t>& placed,
		const int place_more_level, const Bool verbose )
{    
     double clock = WallClockTime( );
     cout << Date( ) << ": entering PlaceMore, peak mem = " 
          << PeakMemUsageGBString( ) << endl;
     vec<int64_t> unplaced;
     for ( int64_t id = 0; id < (int64_t) paths2.size( ); id += 2 )
     {    if ( paths2[id].size( ) == 0 && paths2[id+1].size( ) == 0 )
               unplaced.push_back( id, id+1 );    }
     cout << Date( ) << ": found " << ToStringAddCommas( unplaced.size( ) ) 
          << " unplaced reads" << endl;
     vecbasevector all;
     int nedges = hb.EdgeObjectCount( );
     all.reserve( nedges + unplaced.isize( ) );
     for ( int e = 0; e < nedges; e++ )
          all.push_back( hb.EdgeObject(e) );
     for ( int j = 0; j < unplaced.isize( ); j++ )
          all.push_back( bases[ unplaced[j] ] );
     const int M = 32;
     vec< vec< pair<int,int> > > places( unplaced.size( ) );
     vec<String> prefixes;
     if ( place_more_level == 0 ) prefixes = { "" };
     if ( place_more_level == 1 ) prefixes = { "A", "C", "G", "T" };
     if ( place_more_level == 2 )
     {    prefixes = { "AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA",
               "GC", "GG", "GT", "TA", "TC", "TG", "TT" };    }
     int npasses = IPow( 4, place_more_level );
     for ( int pass = 0; pass < npasses; pass++ )
     {    vec< triple<kmer<M>,int,int> > kmers_plus;
          // cout << Date( ) << ": pass = " << pass 
          //      << ", building kmers_plus, peak mem = " << PeakMemUsageGBString( ) 
          // << endl;
          MakeKmerLookup0Pre( all, prefixes[pass], kmers_plus );
          // cout << Date( ) << ": building places, peak mem = " 
          //      << PeakMemUsageGBString( ) << endl;
          const int64_t batches = 100;
          vec<int64_t> bstart(batches+1);
          for ( int64_t i = 0; i <= batches; i++ )
               bstart[i] = ( (int64_t) kmers_plus.size( ) * i ) / batches;
          for ( int64_t i = 1; i < batches; i++ )
          {    int64_t& s = bstart[i];
               while( s > bstart[i-1] 
                    && kmers_plus[s].first == kmers_plus[s-1].first )
               {    s--;    }    }
          // cout << Date( ) << ": start parallel for loop" << endl;
          #pragma omp parallel for schedule(dynamic, 1)
          for ( int64_t bi = 0; bi < batches; bi++ )
          for ( int64_t i = bstart[bi]; i < bstart[bi+1]; i++ )
          {    int64_t j;
               for ( j = i + 1; j < bstart[bi+1]; j++ )
                    if ( kmers_plus[j].first != kmers_plus[i].first ) break;
               vec<int64_t> g;
               for ( int64_t k = i; k < j; k++ )
                    if ( kmers_plus[k].second < nedges ) g.push_back(k);
               if ( g.solo( ) )
               {    for ( int64_t k = i; k < j; k++ )
                    {    if ( kmers_plus[k].second >= nedges )
                         {    int e = kmers_plus[ g[0] ].second;
                              int epos = kmers_plus[ g[0] ].third; 
                              int rpos = kmers_plus[k].third;
                              pair<int,int> p = make_pair( e, epos - rpos );
                              vec< pair<int,int> >& x 
                                   = places[ kmers_plus[k].second - nedges ];
                              if ( !Member( x, p ) )
                              {    
                                   #pragma omp critical
                                   {    x.push_back(p);    }    }    }    }     }

               i = j - 1;    }    }
     cout << Date( ) << ": placing reads, peak mem = " 
          << PeakMemUsageGBString( ) << endl;
     qvec qv;
     for ( int i = 0; i < unplaced.isize( ); i++ )
     {    UniqueSort( places[i] );
          int min_off = 1000000000;
          for ( int j = 0; j < places[i].isize( ); j++ )
               min_off = Min( min_off, places[i][j].second );
          if ( min_off >= 0 )
          {    vec<Bool> to_delete( places[i].size( ), False );
               int64_t id = unplaced[i];
               quals[id].unpack(&qv);
               int m = bases[id].size( );
               for ( int j = 0; j < places[i].isize( ); j++ )
               {    int e = places[i][j].first, offset = places[i][j].second;
                    m = Min( m, hb.EdgeLengthBases(e) - offset );    }
               vec<int> score( places[i].size( ), 0 );
               vec<int> ids( places[i].size( ), vec<int>::IDENTITY );
               for ( int j = 0; j < places[i].isize( ); j++ )
               {    int e = places[i][j].first, offset = places[i][j].second;
                    for ( int l = 0; l < m; l++ )
                    {    if ( bases[id][l] != hb.EdgeObject(e)[offset+l] )
                              score[j] += qv[l];    }    }
               SortSync( score, ids );
               for ( int j = 1; j < places[i].isize( ); j++ )
                    if ( score[j] > score[0] ) to_delete[ ids[j] ] = True;
               EraseIf( places[i], to_delete );    }

          for ( int j = 0; j < places[i].isize( ); j++ )
          {    int xid = i;
               int64_t id = unplaced[xid];
               int e = places[i][j].first, offset = places[i][j].second;
               if (verbose)
               {    cout << "read " << id
                         << " placed at " << e << "." << offset << "\n";     }
               if ( places[i].solo( ) )
               {    IntVec x;
                    x.push_back(e);
                    ReadPath& p = paths2[id];
                    (IntVec&) p = x;
                    p.setOffset(offset);
		    placed.push_back(id);
                    /*
                    int last_skip 
                         = hb.EdgeLengthBases(e) - bases[id].isize( ) - offset;
                    p.setLastSkip(last_skip);
                    */
                    }    }    }
     cout << TimeSince(clock) << " used in PlaceMore" << endl;    }

void ExtraPaths( const HyperBasevector& hb, const vecbasevector& bases,
     const VecPQVec& quals, ReadPathVec& paths2 )
{    cout << Date( ) << ": setup" << endl;
     vecbasevector edges( hb.EdgeObjectCount( ) );
     for ( int e = 0; e < hb.EdgeObjectCount( ); e++ )
          edges[e] = hb.EdgeObject(e);
     int64_t total_kmers = 0; // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     for ( int e = 0; e < hb.EdgeObjectCount( ); e++ ) // XXXXXXXXXXXXXXXXXXXXX
          total_kmers += hb.EdgeLengthKmers(e); // XXXXXXXXXXXXXXXXXXXXXXXXXXXX
     cout << "total kmers = " << ToStringAddCommas(total_kmers) << endl; // XXX
     const int K0 = 60;
     vec< triple<kmer<K0>,int,int> > kmers_plus;
     MakeKmerLookup0( edges, kmers_plus );
     cout << Date( ) << ": made kmer lookup table" << endl; // XXXXXXXXXXXXXXXX
     vec< kmer<K0> > kmers( kmers_plus.size( ) );
     for ( int64_t i = 0; i < kmers_plus.isize( ); i++ )
          kmers[i] = kmers_plus[i].first;
     const int qfloor = 25;
     int count = 0;
     int64_t npids = bases.size( )/2;
     cout << Date( ) << ": main loop" << endl;
     #pragma omp parallel for
     for ( int64_t pid = 0; pid < npids; pid++ )
     {    int64_t id1 = 2*pid, id2 = 2*pid + 1;
          qvec qv;
          quals[id1].unpack(&qv);
          Bool placed1 = ( paths2[id1].size( ) > 0 );
          Bool placed2 = ( paths2[id2].size( ) > 0 );
          if ( placed1 == placed2 ) continue;
          if (placed1) swap( id1, id2 );
          int stop = qv.size( );
          vec< pair<int,int> > locs;
          kmer<K0> x;
          for ( int j = 0; j <= stop - K0; j++ )
          {    Bool bad = False;
               for ( int l = 0; l < K0; l++ )
               {    if ( qv[j+l] < qfloor )
                    {    bad = True;
                         break;    }    }
               if (bad) continue;
               x.SetToSubOf( bases[id1], j );
               int64_t low = LowerBound( kmers, x ); 
               int64_t high = UpperBound( kmers, x );
               if ( low < high )
               {    for ( int64_t l = low; l < high; l++ )
                    {    locs.push( kmers_plus[l].second, 
                              kmers_plus[l].third - j );    }
                    break;    }    }
          UniqueSort(locs);
          if ( !locs.solo( ) ) continue;
          int64_t id = id1;
          for ( int l = 0; l < locs.isize( ); l++ )
          {    int edge = locs[l].first, pos = locs[l].second;
               IntVec v;
               v.push_back(edge);
               (IntVec&) paths2[id] = v;
               paths2[id].setOffset(pos);
               int64_t pid = id/2;
               // PRINT5( pid, id, edge, pos, x.ToString( ) );    
                    }
          #pragma omp critical
          {    count++;    }    }
     cout << Date( ) << ": placed partners for " << count << " pairs, "
          << PERCENT_RATIO( 3, count, npids ) << " of total" << endl;    }

void LayoutReads( const HyperBasevector& hb, const vec<int>& inv,
     const vecbasevector& bases, const ReadPathVec& paths,
     vec<vec<int>>& layout_pos, vec<vec<int64_t>>& layout_id, 
     vec<vec<Bool>>& layout_or )
{
     int nedges = hb.EdgeObjectCount( );
     layout_pos.resize(nedges), layout_id.resize(nedges), layout_or.resize(nedges);
     for ( int64_t i = 0; i < (int64_t) paths.size( ); i++ )
     {    vec<int> x;
          for ( int j = 0; j < (int) paths[i].size( ); j++ )
               x.push_back( paths[i][j] );
          if ( x.empty( ) ) continue;
          int pos = paths[i].getOffset( );
          for ( int j = 0; j < x.isize( ); j++ )
          {    if ( j > 0 && j < x.isize( ) - 1 ) continue;
               layout_pos[ x[j] ].push_back(pos);
               layout_id[ x[j] ].push_back(i);
               layout_or[ x[j] ].push_back(True);
               pos -= hb.EdgeLengthKmers( x[j] );    }
          x.ReverseMe( );
          for ( int j = 0; j < x.isize( ); j++ )
               x[j] = inv[ x[j] ];
          pos = paths[i].getOffset( ) + bases[i].isize( );
          int len = hb.EdgeLength( x[0] );
          for ( int j = 1; j < x.isize( ); j++ )
               len += hb.EdgeLengthKmers( x[j] );
          pos = len - pos;
          for ( int j = 0; j < x.isize( ); j++ )
          {    if ( j > 0 && j < x.isize( ) - 1 ) continue;
               layout_pos[ x[j] ].push_back(pos);
               layout_id[ x[j] ].push_back(i);
               layout_or[ x[j] ].push_back(False);
               pos -= hb.EdgeLengthKmers( x[j] );    }    }
     #pragma omp parallel for
     for ( int e = 0; e < nedges; e++ )
          SortSync( layout_pos[e], layout_or[e], layout_id[e] );    }

void SortBlobs( const HyperBasevector& hb,
     const vec< triple< pair<int,int>, triple<int,vec<int>,vec<int>>, vec<int> > >&
          blobber,
     vec< pair<int,int> >& blobs )
{    
     vec< vec<int> > comp, compv;
     // hb.ComponentsE(comp); // extremely slow
     hb.Components(compv);
     comp.resize( compv.size( ) );
     #pragma omp parallel for
     for ( int i = 0; i < comp.isize( ); i++ )
     {    const vec<int>& c = compv[i];
          for ( int j = 0; j < c.isize( ); j++ )
          {    int v = c[j];
               for ( int l = 0; l < hb.From(v).isize( ); l++ )
               {    comp[i].push_back( 
                         hb.EdgeObjectIndexByIndexFrom( v, l ) );    }    }
          UniqueSort( comp[i] );    }
     vec<int> compsize( comp.size( ), 0 );
     for ( int i = 0; i < comp.isize( ); i++ )
     for ( int j = 0; j < comp[i].isize( ); j++ )
          compsize[i] += hb.EdgeLengthKmers( comp[i][j] );
     vec<int> ecompsize( hb.EdgeObjectCount( ) );
     for ( int i = 0; i < comp.isize( ); i++ )
     for ( int j = 0; j < comp[i].isize( ); j++ )
          ecompsize[ comp[i][j] ] = compsize[i];
     vec<int> minlen( blobs.size( ) );
     for ( int bl = 0; bl < blobs.isize( ); bl++ )
     {    int i = blobs[bl].first, j = blobs[bl].second;
          vec<int> lefts, rights;
          for ( int k = i; k < j; k++ )
          {    lefts.append( blobber[k].second.second );
               rights.append( blobber[k].second.third );    }
          UniqueSort(lefts), UniqueSort(rights);
          minlen[bl] 
               = Min( ecompsize[ lefts[0] ], ecompsize[ rights[0] ] );    }
     ReverseSortSync( minlen, blobs );    }

void RemoveHangs( HyperBasevector& hb, vec<int>& inv, ReadPathVec& paths, 
     const int max_del )
{    double clock1 = WallClockTime( );
     const double junk_ratio = 10.0;
     vec<kmer_count> kc( hb.EdgeObjectCount( ) );
     for ( int e = 0; e < hb.EdgeObjectCount( ); e++ )
          kc[e].n = hb.EdgeObject(e).isize( ) - hb.K( ) + 1;
     digraphE<kmer_count> shb_kc( hb, kc );
     LogTime( clock1, "removing hanging ends 1" );
     double clock2 = WallClockTime( );
     RemoveHangingEnds( shb_kc, &kmer_count::N, max_del, junk_ratio );
     LogTime( clock2, "removing hanging ends 2" );
     double clock3 = WallClockTime( );
     vec<int> e_to_delete;
     vec<Bool> used;
     shb_kc.Used(used);
     for ( int e = 0; e < hb.EdgeObjectCount( ); e++ )
          if ( !used[e] ) e_to_delete.push_back(e);    
     hb.DeleteEdges(e_to_delete);
     LogTime( clock3, "removing hanging ends 3" );    }

void Insert( VecULongVec& paths2_index, const int e, const int64_t id )
{    if ( !BinMember( paths2_index[e], id ) )
     {    paths2_index[e].push_back(id);
          Sort( paths2_index[e] );    }    }

void Patch( HyperBasevector& hb, const vec< pair<int,int> >& blobs, 
     vec<HyperBasevector>& mhbp, const String& work_dir, const vec<String>& mreport, 
     vecbvec& new_stuff )
{    
     double clock = WallClockTime( );
     new_stuff.clear();
     cout << Date( ) << ": patch reserving space" << endl; // XXXXXXXXXXXXXXXXXXXXXX
     new_stuff.reserve(std::accumulate(mhbp.begin(),mhbp.end(),0ul,
                     [](size_t nnn,HyperBasevector const& hbv)
                     { nnn += hbv.E();
                       auto end = hbv.To().end();
                       auto itr2 = hbv.From().begin();
                       for ( auto itr=hbv.To().begin(); itr!=end; ++itr,++itr2 )
                           nnn += itr->size()*itr2->size();
                       return nnn; }));
     cout << Date( ) << ": memory in use = " << MemUsageGBString( ) << endl;

     int K = hb.K( );
     Ofstream( rout, work_dir + "/patch" );
     for ( int bl = 0; bl < blobs.isize( ); bl++ )
     {    rout << mreport[bl];
          PRINT_TO( rout, mhbp[bl].N( ) );
          if ( mhbp[bl].N( ) > 0 )
          {    HyperBasevector const& hbp = mhbp[bl];
               PRINT_TO( rout, hbp.EdgeObjectCount( ) );

               // Delete detritus.
     
               PRINT2_TO( rout, left, right );

               // Insert patch.

               for ( int e = 0; e < hbp.EdgeObjectCount( ); e++ )
                    new_stuff.push_back( hbp.EdgeObject(e) );
               for ( int v = 0; v < hbp.N( ); v++ )
               for ( int i1 = 0; i1 < hbp.To(v).isize( ); i1++ )
               for ( int i2 = 0; i2 < hbp.From(v).isize( ); i2++ )
               {    basevector const& e1 = hbp.EdgeObjectByIndexTo( v, i1 );
                    basevector const& e2 = hbp.EdgeObjectByIndexFrom( v, i2 );
                    new_stuff.push_back(TrimCat(K,e1,e2));    
                    // new_stuff.back( ).Print( rout, "patch." + ToString(v) 
                    //      + "." + ToString(i1) + "." + ToString(i2) );    
                         }
               rout << "Inserting patch." << endl;    }    }

     cout << TimeSince(clock) << " used patching"
          << ", peak mem usage = " << PeakMemUsageGBString( ) << endl;    }

void Clean200( HyperBasevector& hb, vec<int>& inv, ReadPathVec& paths,
     const vecbasevector& bases, const VecPQVec& quals, const int verbosity )
{
     // Start.

     double clock = WallClockTime( );

     // Heuristics.

     const int min_win = 100;
     const int max_lose = 50;
     const int min_ratio = 5;
     const int max_del = 15;
     const int max_exts = 10;
     const int npasses = 2;

     // Run two passes.

     for ( int zpass = 1; zpass <= npasses; zpass++ )
     {
     
     // ---->

     // Set up.

     vec<int> to_right;
     hb.ToRight(to_right);
     VecULongVec paths_index;
     invert( paths, paths_index, hb.EdgeObjectCount( ) );

     // Look for weak branches.

     cout << Date( ) << ": start walking" << endl;
     const int max_rl = 250;
     vec<int> to_delete;
     #pragma omp parallel for
     for ( int e = 0; e < hb.EdgeObjectCount( ); e++ )
     {    int v = to_right[e];
          if ( !hb.To(v).solo( ) || hb.From(v).size( ) <= 1 ) continue;

          // Find extensions of e.

          int n = hb.From(v).size( );
          int depth = max_rl;
          vec< vec<int> > exts;
          for ( int pass = 1; pass <= 2; pass++ )
          {    exts.clear( );
               for ( int j = 0; j < n; j++ )
               {    vec<int> x = { hb.IFrom( v, j ) };
                    exts.push_back(x);    }
               for ( int i = 0; i < exts.isize( ); i++ )
               {    if ( i >= max_exts ) break;
                    int len = 0;
                    for ( int l = 0; l < exts[i].isize( ); l++ )
                         len += hb.Kmers( exts[i][l] );
                    if ( len >= depth ) continue;
                    int w = to_right[ exts[i].back( ) ];
                    if ( hb.From(w).empty( ) ) 
                    {    depth = Min( depth, len );
                         continue;    }
                    vec<int> p = exts[i];
                    for ( int m = 0; m < hb.From(w).isize( ); m++ )
                    {    vec<int> q(p);
                         q.push_back( hb.IFrom( w, m ) );
                         if ( m == 0 ) exts[i] = q;
                         else exts.push_back(q);    }
                    i--;    }    }
          if ( exts.isize( ) > max_exts ) continue;
          int N = exts.size( );
          vec<int> ei(N);
          for ( int i = 0; i < N; i++ )
          for ( int j = 0; j < n; j++ )
               if ( exts[i][0] == hb.IFrom( v, j ) ) ei[i] = j;

          // Convert to basevectors.

          vec<basevector> bexts, rbexts;
          for ( int i = 0; i < N; i++ )
          {    const vec<int>& x = exts[i];
               basevector b = hb.Cat(x);
               bexts.push_back(b);
               b.ReverseComplement( );
               rbexts.push_back(b);    }

          // Score each read path containing e or a successor of it.

          vec<vec<int>> scores(n);
          vec< pair<int64_t,int> > pi; // {(id,start)}
          for ( int i = 0; i < (int) paths_index[e].size( ); i++ )
          {    int64_t id = paths_index[e][i];
               const ReadPath& p = paths[id];
               for ( int j = 0; j < (int) p.size( ); j++ )
               {    if ( p[j] == e ) 
                    {    int start = p.getOffset( );
                         for ( int l = 0; l <= j; l++ )
                              start -= hb.Kmers( p[l] );
                         pi.push( id, start );    }    }    }
          for ( int m = 0; m < n; m++ )
          {    int ep = hb.IFrom( v, m );
               for ( int i = 0; i < (int) paths_index[ep].size( ); i++ )
               {    int64_t id = paths_index[ep][i];
                    const ReadPath& p = paths[id];
                    for ( int j = 0; j < (int) p.size( ); j++ )
                    {    if ( p[j] == ep ) 
                         {    if ( j > 0 && p[j-1] == e ) continue;
                              int start = p.getOffset( );
                              for ( int l = 0; l < j; l++ )
                                   start -= hb.Kmers( p[l] );
                              pi.push( id, start );    }    }    }    }
          qvec qv;
          for ( int i = 0; i < pi.isize( ); i++ )
          {    int64_t id = pi[i].first; 
               int start = pi[i].second;
               quals[id].unpack(&qv);
               const ReadPath& p = paths[id];
               vec<int> q( N, 0 );
               for ( int pos = 0; pos < depth + hb.K( ) - 1; pos++ )
               {    int rpos = pos - start;
                    if ( rpos < 0 || rpos >= bases[id].isize( ) ) continue;
                    for ( int l = 0; l < N; l++ )
                    {    if ( bexts[l][pos] != bases[id][rpos] )
                              q[l] += qv[rpos];    }    }
               vec<int> qq( n, 1000000000 );
               for ( int l = 0; l < N; l++ )
                    qq[ ei[l] ] = Min( qq[ ei[l] ], q[l] );
               vec<int> idx( n, vec<int>::IDENTITY );
               SortSync( qq, idx );
               if ( qq[0] < qq[1] ) scores[ idx[0] ].push_back( qq[1] - qq[0] );    }

          // Score each read path containing rc of e or its successor.

          vec< pair<int64_t,int> > rpi; // {(id, start of read rel edge re)}
          int re = inv[e];
          for ( int i = 0; i < (int) paths_index[re].size( ); i++ )
          {    int64_t id = paths_index[re][i];
               const ReadPath& p = paths[id];
               for ( int j = 0; j < (int) p.size( ); j++ )
               {    if ( p[j] == re ) 
                    {    int start = p.getOffset( );
                         for ( int l = 0; l < j; l++ )
                              start -= hb.Kmers( p[l] );
                         rpi.push( id, start );    }    }    }
          for ( int m = 0; m < n; m++ )
          {    int rep = inv[ hb.IFrom( v, m ) ];
               for ( int i = 0; i < (int) paths_index[rep].size( ); i++ )
               {    int64_t id = paths_index[rep][i];
                    const ReadPath& p = paths[id];
                    for ( int j = 0; j < (int) p.size( ); j++ )
                    {    if ( p[j] == rep ) 
                         {    if ( j < (int) p.size( ) - 1 && p[j+1] == re ) 
                                   continue;
                              int start = p.getOffset( );
                              for ( int l = 0; l <= j; l++ )
                                   start -= hb.Kmers( p[l] );
                              rpi.push( id, start );    }    }    }    }
          for ( int i = 0; i < rpi.isize( ); i++ )
          {    int64_t id = rpi[i].first; 
               int start = rpi[i].second;
               quals[id].unpack(&qv);
               const ReadPath& p = paths[id];
               vec<int> q( N, 0 );
               for ( int pos = 0; pos < depth + hb.K( ) - 1; pos++ )
               {    int rpos = hb.K( ) - 2 - pos - start;
                    if ( rpos < 0 || rpos >= bases[id].isize( ) ) continue;
                    for ( int l = 0; l < N; l++ )
                    {    int s = bexts[l].size( );
                         if ( rbexts[l][s-pos-1] != bases[id][rpos] )
                              q[l] += qv[rpos];    }    }
               vec<int> qq( n, 1000000000 );
               for ( int l = 0; l < N; l++ )
                    qq[ ei[l] ] = Min( qq[ ei[l] ], q[l] );
               vec<int> idx( n, vec<int>::IDENTITY );
               SortSync( qq, idx );
               if ( qq[0] < qq[1] ) scores[ idx[0] ].push_back( qq[1] - qq[0] );    }

          // Analyze scores.

          for ( int j = 0; j < n; j++ )
               ReverseSort( scores[j] );
          //Bool win = False;
          for ( int d = 0; d <= max_del; d++ )
          {    vec<int> qsum(n);
               for ( int j = 0; j < n; j++ )
               {    for ( int i = 0; i < scores[j].isize( ); i++ )
                    {    if ( scores[j][i] <= d ) break;
                         qsum[j] += scores[j][i];    }    }
               vec<int> ids( n, vec<int>::IDENTITY );
               ReverseSortSync( qsum, ids );
               if ( qsum[0] >= min_win && qsum[1] <= max_lose
                    && qsum[0] >= min_ratio * qsum[1] )
               {    //win = True;
                    #pragma omp critical
                    {    for ( int j = 1; j < n; j++ )
                         {    int e2 = hb.IFrom( v, ids[j] );
                              to_delete.push_back( e2, inv[e2] );    }
                         if ( verbosity >= 1 )
                         {    cout << "\n" << e << " --> ";
                              for ( int j = 0; j < n; j++ )
                              {    int e2 = hb.IFrom( v, ids[j] );
                                   if ( j > 0 ) cout << ",";
                                   cout << e2;    }
                              cout << "\n";
                              for ( int j = 0; j < n; j++ )
                              {    cout << "e" << j+1 << ": ";
                                   cout << printSeq( scores[ ids[j] ] ) << endl;    }
                              cout << "deleting e2";
                              if ( n > 2 ) cout << "-" << n << endl;
                              cout << endl;    }
                         if ( verbosity >= 2 )
                         {    for ( int j = 1; j < n; j++ )
                              {    int e2 = hb.IFrom( v, ids[j] );
                                   hb.EdgeObject(e2).Print( cout,
                                       ToString(zpass) + "." + ToString(e2) );
                                   cout << "-----" << endl;    }    }    }

                    break;    }    }    }

     // Clean up.

     hb.DeleteEdges(to_delete);
     Cleanup( hb, inv, paths );    }
     TestInvolution( hb, inv );
     Validate( hb, paths );    
     cout << TimeSince(clock) << " used cleaning 200-mer graph" << endl;    }

template<class H> void DegloopCore( const int mode, H& hb, vec<int>& inv, 
     ReadPathVec& paths, const vecbasevector& bases, const VecPQVec& quals,
     const VecULongVec& paths_index, const int v, const int pass,
     const double min_dist, vec<int>& EDELS, const int verbosity,
     const vec<int>* ids )
{
     int K = hb.K( );
     int n = ( pass == 1 ? hb.From(v).size( ) : hb.To(v).size( ) );
     if ( n >= 2 )
     {    vec<vec<int>> qs(n);

          // Don't mess with homopolymers.

          Bool homop = False;
          for ( int i = 0; i < n; i++ )
          {    int e = ( pass == 1 ? hb.IFrom( v, i ) : hb.ITo( v, i ) );
               if ( hb.Bases(e) == 0 ) continue;
               int ne = hb.Bases(e);
               vec<char> b;
               const int hcount = 10;
               if ( pass == 1 )
               {    for ( int j = 0; j < hcount; j++ )
                         b.push_back( hb.EdgeObject(e)[ K - j - 1 ] );    }
               else
               {    for ( int j = 0; j < hcount; j++ )
                         b.push_back( hb.EdgeObject(e)[ ne - K + j ] );    }
               UniqueSort(b);
               if ( b.solo( ) ) homop = True;    }
          if (homop) return;

          // Proceed.

          int min_edge = 1000000000;
          for ( int i = 0; i < n; i++ )
          {    int e = ( pass == 1 ? hb.IFrom( v, i ) : hb.ITo( v, i ) );
               if ( hb.Bases(e) == 0 ) continue;
               min_edge = Min( min_edge, hb.Bases(e) );    }
          auto qvItr = quals.begin();
          for ( int i = 0; i < n; i++ )
          {    int e = ( pass == 1 ? hb.IFrom( v, i ) : hb.ITo( v, i ) );
               if ( hb.Bases(e) == 0 ) continue;
               int ne = hb.Bases(e), re = inv[e];
               for ( int xpass = 1; xpass <= 2; xpass++ )
               {    int x = ( xpass == 1 ? e : re );
                    for (int j = 0; j < (int) paths_index[x].size( ); j++)
                    {    int64_t id = paths_index[x][j];
                         const ReadPath& p = paths[id];
                         const basevector& b = bases[id];
                         const qualvector& q = qvItr[id];

                         // Set homopolymer base quality to min across it.

                         /*
                         vec<int> q( bases[id].size( ) );
                         for ( int j = 0; j < b.isize( ); j++ )
                              q[j] = quals[id][j];
                         for ( int j = 0; j < b.isize( ); j++ )
                         {    int k;
                              for ( k = j + 1; k < b.isize( ); k++ )
                                   if ( b[k] != b[j] ) break;
                              int m = 1000000000;
                              for ( int l = j; l < k; l++ )
                                   m = Min( m, (int) q[l] );
                              for ( int l = j; l < k; l++ )
                                   q[l] = m;
                              j = k - 1;    }
                         */

                         for ( int l = 0; l < (int) p.size( ); l++ )
                         {    if ( p[l] != x ) continue;
                              int estart = p.getOffset( ); // start of read on edge
                              for ( int m = 0; m < l; m++ )
                                   estart -= hb.Kmers( p[m] );
                              int estop = estart + b.size( ); // stop of read on edge
                              int rpos = ( xpass == 1 ^ pass == 1 ?
                                   -estart + ne - K : -estart + K - 1 );
                              if ( rpos < 0 || rpos >= b.isize( ) ) continue;

                              // Very stringent condition, probably too stringent
                              // in some cases.

                              if ( !( xpass == 1 ^ pass == 1 ) )
                              {    if ( IntervalOverlap(
                                        0, min_edge, estart, estop ) < K )
                                   {    continue;   }    }
                              else
                              {    if ( IntervalOverlap(
                                        ne - min_edge, ne, estart, estop ) < K )
                                   {    continue;   }    }

                              /*
                              if ( b[rpos] != hb.EdgeObject(x)[
                                   xpass == 1 ^ pass == 1 ? ne-K : K-1 ] )
                              {    continue;    }
                              */

                              if ( verbosity >= 3 )
                              {    ForceAssert( ids != NULL );
                                   cout << "read " << (*ids)[id] 
                                        << " supports edge " << e << " with quality "
                                        << int(q[rpos]) << endl;    }

                              qs[i].push_back( q[rpos] );    }    }    }
               ReverseSort( qs[i] );    }

          // Assess quality score distribution difference.

          vec<double> m( n, -1 );
          vec<int> k(n);
          for ( int i = 0; i < n; i++ )
          {    k[i] = qs[i].size( );
               if ( qs[i].nonempty( ) ) m[i] = Mean( qs[i] );    }
          vec<int> dels;
          vec<double> dists;
          for ( int i1 = 0; i1 < n; i1++ )
          for ( int i2 = 0; i2 < n; i2++ )
          {    if ( i1 == i2 ) continue;
               int good1 = 0;
               for ( int j = 0; j < qs[i1].isize( ); j++ )
                    if ( qs[i1][j] >= 30 ) good1++;
               int good2 = 0;
               for ( int j = 0; j < qs[i2].isize( ); j++ )
                    if ( qs[i2][j] >= 30 ) good2++;
               int e2 = ( pass == 1 ? hb.IFrom( v, i2 ) : hb.ITo( v, i2 ) );
               int ne2 = hb.Kmers(e2);

               if ( mode >= 2 && k[i2] == 0 && good1 >= 10 && ne2 <= 200 )
                    dels.push_back(i2);

               if ( k[i1] == 0 || k[i2] == 0 ) continue;
               double dist = (m[i1]-m[i2]) 
                    / sqrt( m[i1]*m[i1]/k[i1] + m[i2]*m[i2]/k[i2] );
               Bool kill2 = ( dist >= min_dist && good2 <= 1 && ne2 <= 200 );
               if (kill2) dels.push_back(i2);
               if ( kill2 || ( dist >= 0 && verbosity >= 2 ) )
                    dists.push_back(dist);    }
          UniqueSort(dels);

          if ( verbosity >= 2 )
          {    cout << "\n" << ( pass == 1 ? "from" : "to" ) << "\n";
               for ( int i = 0; i < n; i++ )
               {    int e = ( pass == 1 ? hb.IFrom( v, i ) : hb.ITo( v, i ) );
                    cout << e << ": " << printSeq( qs[i] ) << endl;    }
               cout << "dist = " << printSeq(dists) << endl;    }

          // Record edges for deletion.

          if ( dels.nonempty( ) )
          {
               #pragma omp critical
               {    vec<int> edels;
                    if ( verbosity == 1 )
                         cout << "\n" << ( pass == 1 ? "from" : "to" ) << "\n";
                    for ( int i = 0; i < n; i++ )
                    {    int e = ( pass == 1 ? hb.IFrom( v, i ) : hb.ITo( v, i ) );
                         if ( Member( dels, i ) ) edels.push_back(e);
                         if ( verbosity == 1 )
                              cout << e << ": " << printSeq( qs[i] ) << endl;    }
                    Sort(edels);
                    EDELS.append(edels);
                    if ( verbosity >= 1 ) 
                         cout << "delete edges: " << printSeq(edels) << endl;
                    if ( verbosity == 1 )
                    {    cout << "dist = " << printSeq(dists) 
                              << endl;    }    }    }    }    }

template void DegloopCore( const int mode, HyperBasevector& hb, vec<int>& inv, 
     ReadPathVec& paths, const vecbasevector& bases, const VecPQVec& quals,
     const VecULongVec& paths_index, const int v, const int pass,
     const double min_dist, vec<int>& EDELS, const int verbosity,
     const vec<int>* ids );

template void DegloopCore( const int mode, HyperBasevectorX& hb, vec<int>& inv, 
     ReadPathVec& paths, const vecbasevector& bases, const VecPQVec& quals,
     const VecULongVec& paths_index, const int v, const int pass,
     const double min_dist, vec<int>& EDELS, const int verbosity,
     const vec<int>* ids );

// Go through branch points.
// Score branches by computing quality score at Kth base.
// Uses version of quality score distribution test from DivineBubbles, but
// less sophisticated.

void Degloop( const int mode, HyperBasevector& hb, vec<int>& inv, ReadPathVec& paths,
     const vecbasevector& bases, const VecPQVec& quals, const double min_dist,
     const int verbosity )
{    cout << Date( ) << ": start degloop" << endl;
     cout << Date( ) << ": creating path index" << endl;
     VecULongVec paths_index;
     invert( paths, paths_index );

     // Main loop.

     cout << Date( ) << ": starting loop" << endl;
     int K = hb.K( );
     vec<int> EDELS;
     #pragma omp parallel for
     for ( int v = 0; v < hb.N( ); v++ )
     {    for ( int pass = 1; pass <= 2; pass++ )
          {    DegloopCore( mode, hb, inv, paths, bases, quals, paths_index, 
                    v, pass, min_dist, EDELS, verbosity );    }    }
     int ed = EDELS.size( );
     for ( int i = 0; i < ed; i++ )
          EDELS.push_back( inv[ EDELS[i] ] );
     UniqueSort(EDELS);
     hb.DeleteEdges(EDELS);
     cout << Date( ) << ": degloop complete" << endl;    }
