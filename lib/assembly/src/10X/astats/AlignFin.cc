// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

// Align to finished sequence and compute N50 perfect stretch.
// Stats shown for mrx/final:toes.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "MainTools.h"
#include "PrintAlignment.h"
#include "paths/HyperBasevector.h"
#include "random/Shuffle.h"
#include "10X/Gap.h"
#include "10X/Super.h"
#include "10X/astats/AlignFin.h"
#include "10X/astats/GenomeAlign.h"
#include "10X/astats/RefAlign.h"

void AlignFin( 

     // inputs:

     const HyperBasevectorX& hb, vec<int> inv,
     digraphE<vec<int>> D, vec<int> dinv,
     const vec<vec<vec<vec<int>>>>& dlines, 
     const vecbasevector& genome,  // finished sequence
     MasterVec< SerfVec<triple<int,int,int> > > galignsb,  // align->finished

     // outputs:

     vec< vec< vec< triple< vec<int>, align, int > > > >& Matches,
     double& errw,       // weighted error rate
     int64_t& N50PS,     // N50 perfect stretch
     String& report,

     // control:

     const Bool longest_only, // kill all but longest alignment
     const Bool stats, // compute N50 perfect stretch and weighted error rate
     const vec<int>& targets // use only these finished ids

     )

{
     // Set up.

     double clock = WallClockTime( );
     const int HBK = hb.K( );
     vecbasevector tigs = hb.Edges( );
     Matches.resize( genome.size( ) );
     vec<vec<ho_interval>> perfs( genome.size( ) );
     vec< pair<int,int> > err( genome.size( ), make_pair(0,0) );

     // Reinsert loops.

     vec<int> llens;
     GetLineLengths( hb, D, dlines, llens );
     const int MIN_LINE = 10000;
     int nd = D.E( );
     ReinsertLoops( hb, inv, D, dinv );
     vec<int> to_left, to_right;
     D.ToLeft(to_left), D.ToRight(to_right);

     // Identify edges in short lines.

     cout << Date( ) << ": finding edges in short lines" << endl;
     vec<Bool> keep( D.E( ), False );
     for ( int i = 0; i < dlines.isize( ); i++ )
     {    if ( llens[i] < MIN_LINE ) continue;
          vec<int> C = Contents( dlines[i] );
          for ( auto d : C ) keep[d] = True;    }

     // Compute line index.

     for ( int d = nd; d < D.E( ); d++ ) keep[d] = True;
     vec<pair<int,int>> tol( D.E( ), make_pair(-1,-1) );
     for ( int i = 0; i < dlines.isize( ); i++ )
     for ( int j = 0; j < dlines[i].isize( ); j++ )
     for ( int k = 0; k < dlines[i][j].isize( ); k++ )
     for ( int l = 0; l < dlines[i][j][k].isize( ); l++ )
          tol[ dlines[i][j][k][l] ] = make_pair( i, j );

     // Munch.

     cout << Date( ) << ": munching" << endl;
     int TE = tigs.size( );
     Munch( D, dinv, tigs, inv, HBK );

     // Create auxiliary aligns.

     cout << Date( ) << ": creating auxiliary aligns" << endl;
     {    vecbasevector tigs2;
          for ( int e = TE; e < (int) tigs.size( ); e++ ) tigs2.push_back( tigs[e] );
          vec< vec< pair<int,int> > > aligns;
          const int align_genome_K = 80;
          AlignToGenomeCore( tigs2, genome, aligns, align_genome_K );
          galignsb.resize( tigs.size( ) );
          const int MAX_JUMP = 100;
          for ( int e = TE; e < (int) tigs.size( ); e++ )
          for ( int i = 0; i < aligns[e-TE].isize( ); i++ )
          {    int j, g = aligns[e-TE][i].first;
               for ( j = i + 1; j < aligns[e-TE].isize( ); j++ )
               {    if ( aligns[e-TE][j].first != g ) break;
                    if ( aligns[e-TE][j].second - aligns[e-TE][j-1].second 
                         > MAX_JUMP )
                    {    break;    }    }
               galignsb[e].push_back( make_triple( g, aligns[e-TE][i].second,
                    aligns[e-TE][j-1].second + tigs[e].size( ) ) );
               i = j - 1;    }    }

     // Create second alignments.

     cout << Date( ) << ": creating second alignments" << endl;
     int64_t N50_perf = -1;
     MasterVec<SerfVec<refalign>> galigns;
     const Bool allow_second_best = True;
     RefAlign( genome, HBK, tigs, inv, D, dinv, dlines, 
          galignsb, galigns, 0, vec<int>( ), allow_second_best );

     // Delete alignments of edges in short lines.

     for ( int d = 0; d < D.E( ); d++ ) if ( !keep[d] ) galigns[d].resize(0);

     // Go through the finished sequences.

     vec<String> Report( genome.size( ) );
     #pragma omp parallel for schedule(dynamic,1)
     for ( int g = 0; g < (int) genome.size( ); g++ )
     {    if ( targets.nonempty( ) && !BinMember( targets, g ) ) continue;
          ostringstream out;
          out << "============================================================="
               << "=======================\n\n";
          out << "ALIGNMENTS TO FINISHED SEQUENCE " << g << endl << endl;
          vec< pair< vec<int>, align > > matches;

          // Find the edges that appear.

          vec<int> ds;
          for ( int d = 0; d < (int) galigns.size( ); d++ )
          for ( int j = 0; j < (int) galigns[d].size( ); j++ )
               if ( galigns[d][j].chr - 1 == g ) ds.push_back(d);
          UniqueSort(ds);

          // Now sort them by their line location.

          vec<triple<int,int,int>> lds;
          for ( int i = 0; i < ds.isize( ); i++ )
          {    int d = ds[i];
               lds.push( tol[d].first, tol[d].second, d );    } 
          Sort(lds);

          // Cell contents functions.

          auto CellContents = [&]( const int v, const int w )
          {    vec<int> ds = D.IFrom(v);
               ds.append( D.ITo(w) );
               UniqueSort(ds);
               for ( int k = 0; k < ds.isize( ); k++ )
               {    int d = ds[k];
                    int x = to_left[d], y = to_right[d];
                    for ( auto z : {x,y} )
                    {    if ( z == v || z == w ) continue;
                         {    for ( auto f : D.IFrom(z) )
                                   if ( !Member( ds, f ) ) ds.push_back(f);
                              for ( auto f : D.ITo(z) )
                              {    if ( !Member( ds, f ) )
                                   ds.push_back(f);    }    }    }    }
               return ds;    };
          auto CellHasGap = [&]( const int v, const int w )
          {    vec<int> ds = CellContents( v, w );
               for ( auto d : ds ) if ( D.O(d)[0] < 0 ) return True;
               return False;    };

          // Global variables.

          const int K = 40;
          RefAlignScratch<K> s;
          Bool verbose = False;
          String report;
          SerfVec<refalign> x;
          const basevector& GG = genome[g];
          const int mismatch = 10;
          const int gap_open = 10;
          const int gap_extend = 1;

          // Define path alignment functions.

          auto ScorePath = [&]( const vec<int>& p )
          {    vec<int> y;
               for ( auto d : p ) y.append( D.O(d) );
               RefAlignCore<K>( y, HBK, tigs, inv, D, dinv, dlines, genome, 
                    galignsb, x, verbose, report, s, g );
               basevector B = tigs[ y[0] ];
               for ( int l = 1; l < y.isize( ); l++ )
               {    B.resize( B.isize( ) - (HBK-1) );
                    B = Cat( B, tigs[ y[l] ] );    }
               int best_score = 1000000000;
               for ( auto s : x )
               {    if ( s.chr != g + 1 ) continue;
                    const align& a = s.a;
                    int score = 0; 
                    int p1 = a.pos1( ), p2 = a.pos2( );
                    for ( int u = 0; u < a.Nblocks( ); u++ )
                    {    if ( a.Gaps(u) != 0 )
                              score += gap_open + (Abs(a.Gaps(u))-1) * gap_extend;
                         if ( a.Gaps(u) > 0 ) p2 += a.Gaps(u);
                         if ( a.Gaps(u) < 0 ) p1 -= a.Gaps(u);
                         for ( int x = 0; x < a.Lengths(u); x++ )
                         {    if ( B[p1] != GG[p2] ) score += mismatch;
                              ++p1; ++p2;    }    }
                    best_score = Min( best_score, score );    }
               return best_score;    };
          auto PathAlign = [&]( const vec<int>& p, align& a )
          {    vec<int> y;
               for ( auto d : p ) y.append( D.O(d) );
               // to switch to verbose, inscrutable but potentially informative
               // const int verbosity = 1;
               const int verbosity = 0;
               RefAlignCore<K>( y, HBK, tigs, inv, D, dinv, dlines, genome, 
                    galignsb, x, verbosity, report, s, g );
               basevector B = tigs[ y[0] ];
               for ( int l = 1; l < y.isize( ); l++ )
               {    B.resize( B.isize( ) - (HBK-1) );
                    B = Cat( B, tigs[ y[l] ] );    }
               int best_score = 1000000000;
               int best_m = -1;
               for ( int m = 0; m < (int) x.size( ); m++ )
               {    const refalign& s = x[m];
                    if ( s.chr != g + 1 ) continue;
                    const align& a = s.a;
                    int score = 0; 
                    int p1 = a.pos1( ), p2 = a.pos2( );
                    for ( int u = 0; u < a.Nblocks( ); u++ )
                    {    if ( a.Gaps(u) != 0 )
                              score += gap_open + (Abs(a.Gaps(u))-1) * gap_extend;
                         if ( a.Gaps(u) > 0 ) p2 += a.Gaps(u);
                         if ( a.Gaps(u) < 0 ) p1 -= a.Gaps(u);
                         for ( int x = 0; x < a.Lengths(u); x++ )
                         {    if ( B[p1] != GG[p2] ) score += mismatch;
                              ++p1; ++p2;    }    }
                    if ( score < best_score )
                    {    best_score = score;
                         best_m = m;    }    }
               if ( best_m < 0 ) 
               {    /*
                    out << "Not aligned." << endl;
                    // PROBABLY DELETE THIS EXTRA STUFF LATER:
                    const int verbosity = 2;
                    RefAlignCore<K>( y, HBK, tigs, inv, D, dinv, dlines, 
                         genome, galignsb, x, verbosity, report, s, g );
                    out << report;
                    */
                    return False;    }
               else 
               {    a = x[best_m].a;
                    /*
                    // PROBABLY DELETE THIS EXTRA STUFF LATER:
                    if ( !Proper( a, B.size( ), GG.size( ) ) )
                    {    out << "not proper" << endl;
                         const int verbosity = 2;
                         RefAlignCore<K>( y, HBK, tigs, inv, D, dinv, dlines, 
                              genome, galignsb, x, verbosity, report, s, g );
                         out << report;    }
                    */
                    return True;    }    };

          // Form them into chains.

          int count = 0;
          for ( int i = 0; i < lds.isize( ); i++ )
          {    int lid = lds[i].first;
               if ( lid < 0 ) continue;
               const vec<vec<vec<int>>>& L = dlines[lid];

               // Don't allow starting at a cell containing a gap.
     
               {    int cid = lds[i].second; // cell id
                    if ( cid % 2 == 1 )
                    {    if ( cid == L.isize( ) - 1 ) continue; // very weird case
                         int l = cid;
                         int v = to_right[ L[l-1][0][0] ]; 
                         int w = to_left[ L[l+1][0][0] ];
                         if ( CellHasGap( v, w ) ) continue;    }    }

               // Proceed.

               int j;
               vec<int> z = { lds[i].third };
               for ( j = i + 1; j < lds.isize( ); j++ )
               {    int cid = lds[j].second; // cell id

                    // Require same line.

                    if ( lds[j].first != lds[i].first ) break;

                    // Require not skipping a long even cell or a gap cell.

                    Bool do_break = False;
                    for ( int l = lds[j-1].second + 1; l < lds[j].second; l++ )
                    {    if ( l % 2 == 0 ) 
                         {    int n = 0, d = L[l][0][0];
                              for ( int j = 0; j < D.O(d).isize( ); j++ )
                                   n += tigs[ D.O(d)[j] ].isize( ) - HBK + 1;
                              if ( n >= 1000 ) do_break = True;    }
                         else if ( L[l].solo( ) )
                         {    int v = to_right[ L[l-1][0][0] ]; 
                              if ( D.From(v).solo( ) )
                              {    if ( D.OFrom(v,0)[0] < 0 ) 
                                        do_break = True;    }    }    }
                    if (do_break) break;

                    // For now, skip any cell that contains a gap edge.

                    if ( cid % 2 == 1 )
                    {    if ( cid == L.isize( ) - 1 ) break; // very weird case
                         int l = cid;

                         // Check for gap cell.

                         int v = to_right[ L[l-1][0][0] ]; 
                         int w = to_left[ L[l+1][0][0] ];
                         if ( CellHasGap( v, w ) ) break;    }

                    // Save edge.

                    z.push_back( lds[j].third );    }
               int start = lds[i].second, stop = lds[j-1].second;
               
               // Convert chain into a path.  Even cells represent single nongap
               // edges, where odd cells represent potentially complicated stuff,
               // including gaps and cycles.

               int l = start;
               while( l <= stop )
               {    vec<int> master;
                    int l2;
                    for ( l2 = l; l2 <= stop; l2++ )
                    {    
                         // Even cells should represent single nongap edges.

                         if ( l2 % 2 == 0 ) 
                         {    master.push_back( L[l2][0][0] );
                              continue;    }

                         // Deal with odd cells.  First determine the cell contents.
     
                         int d1 = L[l2-1][0][0], d2 = L[l2+1][0][0];
                         int v = to_right[d1], w = to_left[d2];
                         vec<int> ds = CellContents( v, w );
     
                         // Find paths through the cell.  We start with allowing
                         // one-fold duplication and keep going until going further
                         // doesn't help.  Note the danger that this computation
                         // could explode.
                         // Note assumptions:
                         // 1. There is always a path across the cell.
                         // 2. Our approach will find the best path.
     
                         int BEST_SCORE = 1000000000;
                         vec<int> best_path;
                         for ( int max_copies = 1; ; max_copies++ )
                         {    
                              // out << endl; // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                              // PRINT5_TO(out,g,max_copies,l2,start,stop); // XXXXXX
                              // PRINT2_TO( out, L[l2-1][0][0], L[l2+1][0][0] );//XXX
                              vec<vec<int>> paths;
                              double clock = WallClockTime( );
                              D.EdgePathsLim( to_left, to_right, v, w, d2, 
                                   paths, max_copies, -1, -1 );
                              /*
                              if ( WallClockTime( ) - clock > 1 )
                              {    out << Date( ) << ": " << TimeSince(clock) 
                                        << " used finding edge paths"
                                        << endl;    }
                              */
                              // PRINT_TO( out, paths.size( ) ); // XXXXXXXXXXXXXXXX
                              const int MAX_PATHS = 200;
                              if ( paths.isize( ) > MAX_PATHS ) break;
                              vec<int> SCORE( paths.size( ), 1000000000 );
                              clock = WallClockTime( );
                              for ( int m = 0; m < paths.isize( ); m++ )
                              {    vec<int> p = paths[m];

                                   // Check for gap.  We're supposed to have checked
                                   // this already, so it's not clear why it's 
                                   // needed.

                                   Bool gap = False;
                                   for ( auto d : p )
                                        if ( D.O(d)[0] < 0 ) gap = True;
                                   if (gap) continue;

                                   // Expand p on the left and right if possible.
                                   // Then convert to a list y of base edges.

                                   if ( l2 > start ) p.push_front( L[l2-1][0][0]);
                                   if ( l < stop ) p.push_back( L[l2+1][0][0] );
                                   SCORE[m] = ScorePath(p);
                                   // out << "p = " << printSeq(p) // XXXXXXXXXXX
                                   //      << " ==> " << SCORE[m] << endl; // XXXX
                                        }
                              /*
                              if ( WallClockTime( ) - clock > 1 )
                              {    out << Date( ) << ": " << TimeSince(clock) 
                                        << " used scoring " << paths.size( )
                                        << " paths" << endl;    }
                              */
                              if ( Min(SCORE) == BEST_SCORE )
                              {    for ( int k = 0; k < SCORE.isize( ); k++ )
                                   {    if ( SCORE[k] == BEST_SCORE )
                                        {    best_path = paths[k];
                                             break;    }    }
                                   break;    }
                              else BEST_SCORE = Min(SCORE);    }

                         // Append to master unless failure.

                         if ( BEST_SCORE == 1000000000 ) 
                         {    
                              /*
                              out << "align at L" << lid << "." << l2 // XXXXXXXXXXX
                                   << " failed" << endl; // XXXXXXXXXXXXXXXXXXXXXXXX
                              */
                              break;
                                        } // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                         master.append(best_path);    }

                    // Save master.

                    if ( master.nonempty( ) ) 
                    {    double clock = WallClockTime( );
                         align a;
                         Bool aligned = PathAlign( master, a );

                         // Try to extend master.

                         /*
                         if (aligned)
                         {    vec<int> y;
                              for ( auto d : master ) y.append( D.O(d) );
                              basevector B = tigs[ y[0] ];
                              for ( int l = 1; l < y.isize( ); l++ )
                              {    B.resize( B.isize( ) - (HBK-1) );
                                   B = Cat( B, tigs[ y[l] ] );    }
                              Bool changed = False;
                              int v = to_left[ master.front( ) ];
                              if ( a.pos1( ) == 0 && a.pos2( ) > 0 
                                   && D.To(v).solo( ) )
                              {    int d = D.ITo(v,0);
                                   if ( D.O(d)[0] >= 0 )
                                   {    master.push_front(d);
                                        changed = True;    }    }
                              int w = to_right[ master.back( ) ];
                              if ( a.Pos1( ) == B.isize( ) 
                                   && a.Pos2( ) < genome[g].isize( )
                                   && D.From(w).solo( ) )
                              {    int d = D.IFrom(w,0);
                                   if ( D.O(d)[0] >= 0 )
                                   {    master.push_back(d);
                                        changed = True;    }    }
                              if (changed) aligned = PathAlign( master, a );    }
                         */

                         // Save.

                         if (aligned) matches.push( master, a );
                         /*
                         if ( WallClockTime( ) - clock > 1 )
                         {    out << Date( ) << ": " << TimeSince(clock) 
                                   << " used aligning master" << endl;    }    
                         */
                              }
                    l = l2 + 1;    }

               // Advance to next chain.

               i = j - 1;    }

          // Add error count to matches.

          vec< triple< vec<int>, align, int > > matchesx;
          for ( int j = 0; j < matches.isize( ); j++ )
          {    const vec<int>& p = matches[j].first;
               const align& a = matches[j].second;
               vec<int> y;
               for ( auto d : p ) y.append( D.O(d) );
               basevector B = tigs[ y[0] ];
               for ( int l = 1; l < y.isize( ); l++ )
               {    B.resize( B.isize( ) - (HBK-1) );
                    B = Cat( B, tigs[ y[l] ] );    }
               const basevector& GG = genome[g];
               const int mismatch = 10;
               const int gap_open = 10;
               const int gap_extend = 1;
               int score = 0; 
               int p1 = a.pos1( ), p2 = a.pos2( );
               for ( int u = 0; u < a.Nblocks( ); u++ )
               {    if ( a.Gaps(u) != 0 )
                         score += gap_open + (Abs(a.Gaps(u))-1) * gap_extend;
                    if ( a.Gaps(u) > 0 ) p2 += a.Gaps(u);
                    if ( a.Gaps(u) < 0 ) p1 -= a.Gaps(u);
                    for ( int x = 0; x < a.Lengths(u); x++ )
                    {    if ( B[p1] != GG[p2] ) score += mismatch;
                         ++p1; ++p2;    }    }
               matchesx.push( matches[j].first, matches[j].second, score );    }

          // Merge across line boundaries.

          vec<Bool> to_delete( matchesx.size( ), False );
          for ( int i1 = 0; i1 < matchesx.isize( ); i1++ )
          for ( int i2 = 0; i2 < matchesx.isize( ); i2++ )
          {    int d1 = matchesx[i1].first.back( );
               int d2 = matchesx[i2].first.front( );
               if ( to_right[d1] == to_left[d2] 
                    && tol[d2].first >= 0 && tol[d2].second == 0 )
               {    vec<int> p = matchesx[i1].first;
                    p.append( matchesx[i2].first );
                    int score = ScorePath(p);
                    if ( score <= matchesx[i1].third + matchesx[i2].third )
                    {    align a;
                         PathAlign( p, a );
                         matchesx.push( p, a, score );
                         to_delete.push_back(False);    
                         to_delete[i1] = to_delete[i2] = True;    }    }    }
          EraseIf( matchesx, to_delete );

          // Combine matches across gaps.  

          vec< vec< triple< vec<int>, align, int > > > matches2( matchesx.size( ) );
          for ( int i = 0; i < matchesx.isize( ); i++ )
               matches2[i].push_back( matchesx[i] );
          while(1)
          {    vec<Bool> to_delete( matches2.size( ), False );
               for ( int i1 = 0; i1 < matches2.isize( ); i1++ )
               for ( int i2 = 0; i2 < matches2.isize( ); i2++ )
               {    if ( to_delete[i1] || to_delete[i2] ) continue;
                    int d1 = matches2[i1].back( ).first.back( );
                    int d2 = matches2[i2].front( ).first.front( );
                    if ( tol[d1].first == tol[d2].first
                         && tol[d2].second == tol[d1].second + 2
                         && D.From( to_right[d1] ).solo( )
                         && D.OFrom( to_right[d1], 0 )[0] < 0 )
                    {    matches2[i1].append( matches2[i2] );
                         to_delete[i2] = True;    }    }
               if ( Sum(to_delete) == 0 ) break;
               EraseIf( matches2, to_delete );    }

          // Delete inferior placements.

          for ( int i = 0; i < matches2.isize( ); i++ )
          {    int errsi = 0;
               for ( auto x : matches2[i] ) errsi += x.third;
               Bool beaten = False;
               for ( int j = 0; j < matches2.isize( ); j++ )
               {    int errsj = 0;
                    for ( auto x : matches2[j] ) errsj += x.third;
                    if ( errsj >= errsi ) continue;
                    Bool cov = True;
                    for ( auto x : matches2[i] )
                    {    Bool c = False;
                         for ( auto y : matches2[j] )
                         {    if ( x.second.pos2( ) >= y.second.pos2( )
                                   && x.second.Pos2( ) <= y.second.Pos2( ) )
                              {    c = True;
                                   break;    }    }
                         if ( !c ) 
                         {    cov = False;
                              break;    }    }
                    if (cov)
                    {    beaten = True;
                         break;    }    }
               if ( !beaten ) Matches[g].push_back( matches2[i] );    }

          // Delete matches having no long perfect stretch.

          const int MINPERF = 400;
          vec<Bool> to_delete3( Matches[g].size( ), False );
          for ( int i = 0; i < Matches[g].isize( ); i++ )
          {    int maxperf = 0;
               for ( int j = 0; j < Matches[g][i].isize( ); j++ )
               {    const vec<int>& p = Matches[g][i][j].first;
                    const align& a = Matches[g][i][j].second;
                    vec<int> y;
                    for ( auto d : p ) y.append( D.O(d) );
                    basevector B = tigs[ y[0] ];
                    for ( int l = 1; l < y.isize( ); l++ )
                    {    B.resize( B.isize( ) - (HBK-1) );
                         B = Cat( B, tigs[ y[l] ] );    }
                    maxperf = Max( maxperf, a.MaxPerf( B, GG ) );    }
               if ( maxperf < MINPERF ) to_delete3[i] = True;    }
          EraseIf( Matches[g], to_delete3 );

          // List canned results.

          vec<vec<vec<int>>> results = {
               { {0,20261,20}, {20453,92196,156} },              // #0
               { {0,98770,31} },                                 // #1
               { {} }, // not fully aligned on either end        // #2 NOT
               { {0,51923,32} },                                 // #3
               { {0,64445,93}, {64709,80645,68},                 // #4
                 {80930,142094,95} },
               { {0,56487,0}, {56598,83589,13},                  // #5
                 {83803,92357,30} },
               { {0,5060,47}, {6729,21920,0}, {22744,35984,59},  // #6
                 {36033,57339,73}, {57490,68629,0}, 
                 {68612,100976,500} },
               { {0,17116,30}, {23215,24089,20},                 // #7
                 {23982,63460,71}, {63576,109155,90} },
               { {0,28561,41}, {28544,59231,24} },               // #8
               { {} }, // looks unsolved due to misassembly      // #9 NOT
               { {} }, // multiple placements                    // #10 NOT
               { {0,137387,98} },                                // #11
               { {0,11400,0}, {11482,57876,0} },                 // #12
               // ------------------------------------------------------------------
               { {} },  // multiple placements                   // #13 NOT
               { {0,72535,81} },                                 // #14
               { {0,71872,85} },                                 // #15
               { {0,58648,19}, {61044,77913,32} },               // #16
               { {0,126270,230} },                               // #17
               { {0,65646,117} },                                // #18
               { {0,35231,70}, {35019,53140,10},                 // #19
                 {53209,147588,71}, {147786,186044,20} },
               { {0,19158,10}, {20656,52552,0} },                // #20
               { {0,96470,31}, {97447,130607,20} },              // #21
               { {0,61240,274}, {63117,75076,0} },               // #22
               { {0,19111,187} },                                // #23
               // ------------------------------------------------------------------
               { {0,67505,140}, {67812,91690,601} },             // #24
               { {0,87206,239} },                                // #25
               { {0,119241,172} },                               // #26
               { {} },  // slighly worse alignment               // #27 NOT
               { {0,21079,0}, {21007,63660,20} },                // #28
               { {0,97037,200} },                                // #29
               { {0,2593,0}, {2635,23564,182},                   // #30
                 {29324,74562,31} },
               { {0,97259,133}, {97403,116137,20} },             // #31
               { {} },  // looks unsolved due to misassembly     // #32 NOT
               { {0,144546,109} },                               // #33
               { {0,139011,203} },                               // #34
               { {0,100948,348} },                               // #35
               { {} },  // looks unsolved due to misassembly     // #36 NOT
               // ------------------------------------------------------------------
               { {0,67033,44} },                                 // #37
               { {} },  // not solved                            // #38 NOT
               { {0,12070,0} },                                  // #39
               { {0,27003,81}, {27055,40645,41},                 // #40
                 {40748,73333,60}, {73411,82939,483} },
               { {0,2518,0}, {2326,48832,61} },                  // #41
               { {} },  // not solved                            // #42 NOT
               { {0,139563,38} },                                // #43
               { {} },  // not solved                            // #44 NOT
               { {0,6163,21}, {6466,33028,202},                  // #45
                 {33023,68243,225}, {68205,68845,20}, 
                 {69016,69701,0}, {69702,70122,0}, 
                 {70332,83674,20}, {83786,103452,93}, 
                 {103235,191308,90} },
               { {0,47330,24}, {47362,152395,356} },             // #46
               { {} },  // looks unsolved due to misassembly     // #47 NOT
               { {0,30709,40}, {30950,112066,273} },             // #48
               { {0,140692,25}, {143761,156054,0} },             // #49
               { {} },  // looks unsolved due to misassembly     // #50 NOT
               { {0,19186,31}, {19402,126990,280} },             // #51
               { {} },  // not solved                            // #52 NOT
               // ------------------------------------------------------------------
               { {0,42938,87} },                                 // #53
               { {0,145167,44} },                                // #54
               { {0,167662,51}, {167716,186882,173} },           // #55
               { {0,104594,24} },                                // #56
               { {0,4549,0} },                                   // #57
               { {} },  // not solved                            // #58 NOT
               { {0,20200,460}, {19895,61566,30} },              // #59
                                                                 // solved = 46
               };

          // Kill all but longest only if requested.

          if ( longest_only && Matches[g].size( ) > 1 )
          {    vec<int> lens( Matches[g].size( ), 0 );
               vec<int> ids( Matches[g].size( ), vec<int>::IDENTITY );
               for ( int i = 0; i < Matches[g].isize( ); i++ )
               {    for ( int j = 0; j < Matches[g][i].isize( ); j++ )
                         lens[i] += Matches[g][i][j].second.extent2( );    }
               ReverseSortSync( lens, ids );
               if ( ids[0] != 0 ) Matches[g][0] = Matches[g][ ids[0] ];
               Matches[g].resize(1);    }

          // Compute perfect match blocks and tally errors.

          for ( int i = 0; i < Matches[g].isize( ); i++ )
          {    const vec< triple< vec<int>, align, int > >& M = Matches[g][i];
               for ( int j = 0; j < M.isize( ); j++ )
               {    const vec<int>& p = M[j].first;
                    const align& a = M[j].second;
                    vec<int> y;
                    for ( auto d : p ) y.append( D.O(d) );
                    basevector B = tigs[ y[0] ];
                    for ( int l = 1; l < y.isize( ); l++ )
                    {    B.resize( B.isize( ) - (HBK-1) );
                         B = Cat( B, tigs[ y[l] ] );    }
                    vec<ho_interval> pi; 
                    a.PerfectIntervals2( B, GG, pi );
                    const int mismatch = 10;
                    const int gap_open = 10;
                    const int gap_extend = 1;
                    int p1 = a.pos1( ), p2 = a.pos2( );
                    err[g].second += a.extent1( );
                    for ( int u = 0; u < a.Nblocks( ); u++ )
                    {    if ( a.Gaps(u) != 0 )
                         {    err[g].first += gap_open 
                                   + (Abs(a.Gaps(u))-1) * gap_extend;    }
                         if ( a.Gaps(u) > 0 ) p2 += a.Gaps(u);
                         if ( a.Gaps(u) < 0 ) p1 -= a.Gaps(u);
                         for ( int x = 0; x < a.Lengths(u); x++ )
                         {    if ( B[p1] != GG[p2] ) err[g].first += mismatch;
                              ++p1; ++p2;    }    }
                    for ( int k = 0; k < pi.isize( ); k++ )
                    {    ho_interval h = pi[k];
                         int s = perfs[g].empty( ) ? 0 : perfs[g].back( ).Stop( );
                         if ( h.Start( ) > s )
                         {    for ( int m = s; m < h.Start( ); m++ )
                                   perfs[g].push( m, m+1 );    }
                         int start = Max( s, h.Start( ) ), stop = h.Stop( );
                         if ( stop > start ) 
                              perfs[g].push( start, stop );    }    }    }
          if ( perfs[g].nonempty( ) )
          {    while( perfs[g].back( ).Stop( ) < GG.isize( ) )
               {    perfs[g].push( perfs[g].back( ).Stop( ),
                         perfs[g].back( ).Stop( ) + 1 );    }    }
          else
          {    for ( int s = 0; s < GG.isize( ); s++ )
                    perfs[g].push( s, s+1 );    }
                    
          // Print.

          int mcount = 0;
          for ( int i = 0; i < Matches[g].isize( ); i++ )
          {    const vec< triple< vec<int>, align, int > >& M = Matches[g][i];
               out << "PLACEMENT " << ++mcount << endl << endl;
               out << "result = {";
               vec<vec<int>> res;
               for ( int k = 0; k < M.isize( ); k++ )
               {    const align& a = M[k].second;
                    res.push( vec<int> { a.pos2( ), a.Pos2( ), M[k].third } );
                    out << " {" << printSeq( res.back( ) ) << "}";
                    if ( k < M.isize( ) - 1 ) out << ",";    }
               out << " }";
               if ( Matches[g].solo( ) && g < results.isize( ) && res == results[g] )
                    out << " ✔";
               out << endl << endl;
               for ( int k = 0; k < M.isize( ); k++ )
               {    if ( k > 0 ) out << "(gap)\n" << endl;
                    out << "[" << g << "=";
                    for ( int s = 0; s < M[k].first.isize( ); s++ )
                    {    int d1 = M[k].first[s];
                         int lid = tol[d1].first, t;
                         for ( t = s + 1; t < M[k].first.isize( ); t++ )
                         {    int lid2 = tol[ M[k].first[t] ].first;
                              if ( lid2 >= 0 && lid2 != lid ) break;    }
                         int d2 = M[k].first[t-1];
                         int c1 = tol[d1].second, c2 = tol[d2].second;
                         if ( s > 0 ) out << ",";
                         out << "L" << lid << "." << c1 << "-" << c2;
                         s = t - 1;    }
                    const vec<int>& p = M[k].first;
                    out << "] " << "master = " << printSeq(p) << endl;
                    const align& a = M[k].second;
                    out << "\naligned to " << g << "." << a.pos2( ) << "-" 
                         << a.Pos2( );
                    const basevector& GG = genome[g];
                    if ( a.Pos2( ) == GG.isize( ) ) out << " = end";
                    out << ", errs = " << M[k].third/10.0 << endl;
                    vec<int> y;
                    for ( auto d : p ) y.append( D.O(d) );
                    basevector B = tigs[ y[0] ];
                    for ( int l = 1; l < y.isize( ); l++ )
                    {    B.resize( B.isize( ) - (HBK-1) );
                         B = Cat( B, tigs[ y[l] ] );    }
                    PrintVisualAlignmentClean( True, out, B, GG, a );    }    }
          Report[g] = out.str( );    }
     for ( int g = 0; g < (int) genome.size( ); g++ ) report += Report[g];

     if (stats)

     {    // Compute weighted error rate.

          int64_t num = 0, den = 0;
          for ( int g = 0; g < (int) genome.size( ); g++ )
          {    num += err[g].first;
               den += err[g].second;    }
          errw = (num/10.0) / den;
     
          // Compute N50, allowing for random changing of the Fosmids.
     
          const int npasses = 1000;
          vec<int> N50s(npasses);
          vec<vec<int>> R(npasses);
          for ( int pass = 0; pass < npasses; pass++ )
          {    R[pass] = vec<int>( genome.size( ), vec<int>::IDENTITY );
               Shuffle( R[pass].begin( ), R[pass].end( ) );    }
          cout << Date( ) << ": computing N50" << endl;
          #pragma omp parallel for
          for ( int pass = 0; pass < npasses; pass++ )
          {    vec<int> q;
               for ( int gg = 0; gg < (int) genome.size( ); gg++ )
               {    int g = R[pass][gg];
                    if ( gg > 0 )
                    {    q.back( ) += perfs[g].front( ).Length( );
                         for ( int j = 1; j < perfs[g].isize( ); j++ )
                              q.push_back( perfs[g][j].Length( ) );    }
                    else
                    {    for ( int j = 0; j < perfs[g].isize( ); j++ )
                              q.push_back( perfs[g][j].Length( ) );    }    }
               N50s[pass] = N50(q);    }
          N50PS = int(round( Mean(N50s)));    }

     // Done.

     cout << Date( ) << ": done, time used = " << TimeSince(clock) 
          << ", peak mem = " << PeakMemUsageGBString( ) << endl;    }
