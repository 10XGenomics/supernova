// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

// Note that much of the code here that assays assemblies using a reference
// sequence is designed for HUMAN samples.  All of it would have to be carefully
// reviewed and modified to allow for nonhuman samples.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"
#include "ParallelVecUtilities.h"
#include "VecUtilities.h"
#include "graph/DigraphTemplate.h"
#include "math/Functions.h"
#include "paths/HyperBasevector.h"
#include "paths/long/large/Lines.h"
#include "10X/DfTools.h"
#include "10X/Gap.h"
#include "10X/Heuristics.h"
#include "10X/LineLine.h"
#include "10X/Super.h"
#include "10X/astats/AlignFin.h"
#include "10X/astats/MeasureGaps.h"
#include "10X/astats/Misassembly.h"
#include "10X/astats/View.h"
#include "10X/MakeHist.h"
#include "10X/MakeLocalsTools.h"

template <class E>
String PrettyFormat( E value, const String unit="", const int MAX_EXP=9 )
{
     double x = value;
     map<int,String> suffix;
     suffix[0] = " ";
     suffix[3] = "K";
     suffix[6] = "M";
     suffix[9] = "G";
     double y = x;
     int s = 0;
     while ( y > 1000 && s < MAX_EXP ) {
          y/=1000;
          s += 3;
     }
     ostringstream out;
     out << fixed << setprecision(2) << setw(7) << right << y
         << " " << suffix[s] << (unit.empty() ? String(" ") : unit);
     
     return String(out.str());
}

template<class T>
void EraseIf( SerfVec<T>& v, const vec<Bool>& to_delete )
{    SerfVec<T> v2;
     for ( int64_t i = 0; i < (int64_t) v.size( ); i++ )
          if ( !to_delete[i] ) v2.push_back( v[i] );
     v = v2;    }

void ReportAssemblyStats( const vec<int64_t>& bci, const vecbasevector& genome,
     const vec< pair<int,ho_interval> >& ambint, const HyperBasevectorX& hb,
     const vec<int>& inv, digraphE<vec<int>> D, vec<int> dinv,
     const vec<vec<vec<vec<int>>>>& dlines, 
     MasterVec< SerfVec<triple<int,int,int> > >& alignsb, const vecbasevector& G, 
     ostream& out, const String& DIR, const String& OUTDIR )
{
     // Get checksum.

     auto checksum = CheckSum( hb, inv, D, dinv );

     // Declare the minimum line that we'll use.

     const int MIN_LINE0 = 1000;
     const int MIN_LINE = 10000;

     // Compute some auxiliary stuff.

     vec<int> kmers( hb.E( ) ), llens;
     // We have two versions of GetLineLengths, one with order D, hb
     // and one with order hb, D.  I hope these are identical.  The order
     // used below also appears in one other place in this file.
     GetLineLengths( D, hb, dlines, llens );
     #pragma omp parallel for
     for ( int e = 0; e < hb.E( ); e++ )
          kmers[e] = hb.Kmers(e);
     int K = hb.K( );

     // Get linelocs.

     vec<vec< pair<int,int> >> linelocs( kmers.size( ) );
     for ( int i = 0; i < dlines.isize( ); i++ )
     for ( int j = 0; j < dlines[i].isize( ); j++ )
     for ( int k = 0; k < dlines[i][j].isize( ); k++ )
     for ( int l = 0; l < dlines[i][j][k].isize( ); l++ )
     {    int d = dlines[i][j][k][l];
          if ( D.O(d)[0] < 0 ) continue;
          for ( int m = 0; m < D.O(d).isize( ); m++ )
          {    int e = D.O(d)[m];
               linelocs[e].push( i, j );    }    }

     // Track misassembly stats.

     int64_t total_err_dis_num = 0, total_err_dis_den = 0;
     int64_t total_err_ori_num = 0, total_err_ori_den = 0;
     int64_t total_err_ord_num = 0, total_err_ord_den = 0;

     // Compute coverage.  Horrible, done by deleting "one chromosome".

     auto per = [&]( double n, double d )
     {    ostringstream out;
          out << fixed << setprecision(2) << setw(7) << right << 100*n/d;
          return out.str( );    };
     int64_t genome_size = 0, cap_gap_total = 0, cov_total = 0;
     int64_t contig_line_N50 = 0, pairtig_N50 = 0;

     // Find lines of lines.

     cout << Date( ) << ": finding lines of lines" << endl;
     vec<vec<vec<vec<int>>>> dlines2;
     FindLineLines( D, dinv, dlines, dlines2 );
     BinaryWriter::writeFile( OUTDIR + "/a.sup.linelines", dlines2 );
     vec<int> linv;
     LineInv( dlines, dinv, linv );

     // For each line line bubble, delete one branch.

     vec<int> dels;
     for ( int i = 0; i < dlines2.isize( ); i++ )
     {    const vec<vec<vec<int>>>& L = dlines2[i];
          for ( int j = 0; j < L.isize( ); j++ )
          {    if ( L[j].size( ) == 2 && L[j][0].solo( ) && L[j][1].solo( ) )
               {    vec<int> x = Contents( L[j] );
                    if ( x.solo( ) ) continue; // should never happen
                    vec<int> y = { linv[ L[j][0][0] ], linv[ L[j][1][0] ] };
                    Sort(y);
                    if ( y < x ) continue;
                    if ( llens[ x[0] ] >= llens[ x[1] ] )
                         dels.append( Contents( dlines[ x[1] ] ) );
                    else dels.append( Contents( dlines[ x[0] ] ) );    }    }    }
     int nd = dels.size( );
     for ( int i = 0; i < nd; i++ ) dels.push_back( dinv[ dels[i] ] );
     digraphE<vec<int>> F(D);
     vec<int> finv(dinv);
     cout << Date( ) << ": deleting edges from F" << endl;
     F.DeleteEdgesParallel(dels);
     RemoveUnneededVertices( F, finv );
     CleanupCore( F, finv );

     // Find lines, and locations of base edges on them.

     vec<vec<vec<vec<int>>>> flines;
     cout << Date( ) << ": finding lines again" << endl;
     FindLines( F, finv, flines, MAX_CELL_PATHS, MAX_CELL_DEPTH, False );
     cout << Date( ) << ": computing linelocs" << endl;
     vec<vec< pair<int,int> >> linelocs2( kmers.size( ) );
     for ( int i = 0; i < flines.isize( ); i++ )
     for ( int j = 0; j < flines[i].isize( ); j++ )
     for ( int k = 0; k < flines[i][j].isize( ); k++ )
     for ( int l = 0; l < flines[i][j][k].isize( ); l++ )
     {    int d = flines[i][j][k][l];
          if ( F.O(d)[0] < 0 ) continue;
          for ( int m = 0; m < F.O(d).isize( ); m++ )
          {    int e = F.O(d)[m];
               linelocs[e].push( i, j );    }    }
     #pragma omp parallel for
     for ( int e = 0; e < hb.E( ); e++ )
          UniqueSort( linelocs[e] );
     vec<int> llens2, to_left2, to_right2;
     GetLineLengths( F, hb, flines, llens2 );
     F.ToLeft(to_left2), F.ToRight(to_right2);

     // Compute contig N50.

     cout << Date( ) << ": computing contig N50" << endl;
     for ( int pass = 1; pass <= 2; pass++ )
     {    vec<int> lens;
          for ( int i = 0; i < flines.isize( ); i++ )
          {    if ( llens2[i] < MIN_LINE ) continue;
               int pos = 0;
               const vec<vec<vec<int>>>& L = flines[i];
               for ( int j = 0; j < L.isize( ); j++ )
               {    const vec<vec<int>>& M = L[j];

                    // Break contig at cell if every path in the cell
                    // contains either a gap (pair or barcode), or a cycle.

                    Bool gap = True;
                    for ( int k = 0; k < M.isize( ); k++ )
                    {    if ( M[k].empty( ) )
                         {    int g = -1;
                              if ( j > 0 && L[j-1].nonempty( )
                                   && L[j-1][0].nonempty( ) )
                              {    int d = L[j-1][0][0];
                                   int v = to_right2[d];
                                   if ( v >= 0 ) g = F.IFrom(v,0);    }
                              if ( IsSequence( F.O(g) ) ) gap = False;
                              if ( pass == 2 )
                              {    if ( g >= 0 && IsPairGap( F.O(g) ) )
                                        gap = False;    }    }
                         else
                         {    Bool gappy = False;
                              const vec<int>& x = M[k];
                              for ( int l = 0; l < x.isize( ); l++ )
                              {    if ( pass == 2 && IsPairGap( F.O(x[l]) ) )
                                        continue;
                                   if ( IsSequence( F.O(x[l]) ) ) continue;
                                   if ( F.O(x[l])[0] < 0 ) gappy = True;    }
                              if ( !gappy ) gap = False;    }    }
                    if (gap)
                    {    if ( pos >= 1 ) lens.push_back(pos);
                         pos = 0;    }
                    vec<int> lensj;
                    for ( int k = 0; k < L[j].isize( ); k++ )
                    {    int len = 0;
                         for ( int l = 0; l < L[j][k].isize( ); l++ )
                         {    int d = L[j][k][l];
                              if ( IsSequence( F.O(d) ) )
                              {    int ltrim, rtrim;
                                   basevector x;
                                   GapToSeq( F.O(d), ltrim, rtrim, x );
                                   len += x.isize( ) - K + 1;
                                   int v = to_left2[d];
                                   if ( F.IFrom(v,0) == d )
                                        len = len - ltrim - rtrim;
                                   continue;    }
                              else if ( F.O(d)[0] < 0 ) continue;
                              for ( int m = 0; m < F.O(d).isize( ); m++ )
                                   len += hb.Kmers( F.O(d)[m] );    }
                         lensj.push_back(len);    }
                    Sort(lensj);
                    if ( lensj.nonempty( ) ) pos += Median(lensj);    }
               if ( pos >= 1 ) lens.push_back(pos);    }
          if ( lens.nonempty( ) )
          {    if ( pass == 1 ) {
                    contig_line_N50 = N50(lens);
                    // contig bin size 1 kb
                    ComputeAndWriteHist( lens, 1000, OUTDIR+"/stats",
                         "contig", "final", False );
               } else pairtig_N50 = N50(lens);    }    }

     // Proceed with flines analysis.

     if ( genome.size( ) > 0 )
     {    
          // Quantify misassemblies.

          // #pragma omp parallel for
          for ( int L = 0; L < flines.isize( ); L++ )
          {    if ( llens2[L] < MIN_LINE ) continue;
               SerfVec< quad<int,Bool,ho_interval,ho_interval> > view;
               View( L, K, kmers, inv, F, flines, linelocs2, alignsb, view );
               Misassembly( view, total_err_dis_num, total_err_dis_den,
                    total_err_ori_num, total_err_ori_den, total_err_ord_num,
                    total_err_ord_den );    }

          // Compute coverage.

          vec< quad< int, ho_interval, int, vec<String> > > all_notes;
          cout << Date( ) << ": start main loop" << endl;
          #pragma omp parallel for schedule( dynamic, 100 )
          for ( int L = 0; L < flines.isize( ); L++ )
          {    if ( llens2[L] < MIN_LINE ) continue;
               SerfVec< quad<int,Bool,ho_interval,ho_interval> > view;
               View( L, K, kmers, inv, F, flines, linelocs2, alignsb, view );

               // Delete short and rc segments.

               vec<Bool> to_delete( view.size( ), False );
               for ( int i = 0; i < (int) view.size( ); i++ )
               {    if ( !view[i].second ) to_delete[i] = True;
                    else
                    {    int gstart = view[i].third.Start( );
                         int gstop = view[i].third.Stop( );
                         if ( gstop - gstart < MIN_LINE )
                         {    to_delete[i] = True;    }    }    }
               EraseIf( view, to_delete );

               // Merge proximate segments.

               Sort(view);
               to_delete.resize_and_set( view.size( ), False );
               for ( int i = 0; i < (int) view.size( ); i++ )
               {    int j;
                    for ( j = i + 1; j < (int) view.size( ); j++ )
                    {    if ( view[j].first != view[i].first ) break;
                         if ( view[j].third.Start( ) - view[j-1].third.Stop( )
                              > 50000 )
                         {    break;    }    }
                    view[i].third.SetStop( view[j-1].third.Stop( ) );
                    view[i].fourth.SetStop( view[j-1].fourth.Stop( ) );
                    for ( int k = i + 1; k < j; k++ ) to_delete[k] = True;
                    i = j - 1;    }
               EraseIf( view, to_delete );

               // Convert to notes.

               vec< quad< int, ho_interval, int, vec<String> > > notes;
               for ( int i = 0; i < (int) view.size( ); i++ )
               {    int gstart = view[i].third.Start( );
                    int gstop = view[i].third.Stop( );
                    int g = view[i].first;
                    int astart = view[i].fourth.Start( );
                    int astop = view[i].fourth.Stop( );
                    String chr = ToString(g+1);
                    if ( g == 22 ) chr = "X";
                    if ( g == 23 ) chr = "Y";
                    ostringstream out0, out0b, out1, out2, out3, out4, out5;
                    out0 << "L" << L;
                    // out0b << ToString( COV[L], 2 );
                    out1 << chr << ":" << ToString(gstart/1000000.0,3) + "-"
                         + ToString(gstop/1000000.0,3);
                    out2 << ToStringAddCommas(gstop-gstart);
                    out3 << ToStringAddCommas(astop-astart);

                    // Compute "gap", the estimated sum of captured gaps within
                    // the unit.

                    int64_t gap = (gstop-gstart) - (astop-astart);
                    out4 << ToStringAddCommas(gap);
                    out5 << fixed << setprecision(2) << right
                         << (100.0 * gap) / (gstop-gstart) << "%";
                    notes.push( g, view[i].third, gap, vec<String> { out0.str( ),
                         // out0b.str( ), 
                         out1.str( ), out2.str( ), out3.str( ),
                         out4.str( ), out5.str( ) } );    }
               #pragma omp critical
               {    all_notes.append(notes);    }    }
          cout << Date( ) << ": sorting notes" << endl;
          Sort(all_notes);

          // Determine if chromosome Y appears to be present.  Note different
          // treatment of ambiguous bases in the reference, because there are so
          // many of them on Y.  Data for cov_Y on which Y threshold was based:
          // NA12878 (female):  0.33%
          // NA24385 (male):   44.29%
          // HGP (male):       42.91%.
          // All these are from 1200M reads.  We set a low threshold because 
          // for lower coverage, or lower quality reads, coverage for a male sample
          // might be much lower.

          int64_t genome_size_Y_not_N = genome[23].size( );
          for ( int l = 0; l < ambint.isize( ); l++ )
          {    if ( ambint[l].first == 23 ) 
                    genome_size_Y_not_N -= ambint[l].second.Length( );    }
          int64_t cap_gap_total_Y = 0, cov_total_Y = 0;
          for ( int i = 0; i < all_notes.isize( ); i++ )
          {    int j, g = all_notes[i].first;
               if ( g != 23 ) continue;
               vec<ho_interval> cov;
               for ( j = i + 1; j < all_notes.isize( ); j++ )
                    if ( all_notes[j].first != g ) break;
               for ( int k = i; k < j; k++ )
               {    cap_gap_total_Y += Max( 0, all_notes[k].third );
                    cov.push_back( all_notes[k].second );    }
               cov_total_Y += TotalCovered(cov);
               i = j - 1;    }
          double cap_frac_Y = cap_gap_total_Y / double(genome_size_Y_not_N);
          double uncap_frac_Y 
               = (genome_size_Y_not_N-cov_total_Y) / double(genome_size_Y_not_N);
          double cov_Y = 1 - cap_frac_Y - uncap_frac_Y;
          cout << Date( ) << ": nominal coverage of chromosome Y = "
               << per( cov_Y*10000, 10000 ) << "%" << endl;
          Bool male = ( cov_Y >= 0.05 );
          cout << Date( ) << ": calling this sample "
               << ( male ? "MALE" : "FEMALE" ) << endl;
          int max_chr = ( male ? 23 : 22 );

          // Compute coverage.

          for ( int g = 0; g <= max_chr; g++ )
               genome_size += genome[g].size( );
          for ( int i = 0; i < all_notes.isize( ); i++ )
          {    int j, g = all_notes[i].first;
               if ( g > max_chr ) break;
               vec<ho_interval> cov;
               for ( j = i + 1; j < all_notes.isize( ); j++ )
                    if ( all_notes[j].first != g ) break;
               for ( int k = i; k < j; k++ )
               {    cap_gap_total += Max( 0, all_notes[k].third );
                    cov.push_back( all_notes[k].second );    }
               for ( int l = 0; l < ambint.isize( ); l++ )
                    if ( ambint[l].first == g ) cov.push_back( ambint[l].second );
               cov_total += TotalCovered(cov);
               i = j - 1;    }

          // Compute genome view.

          {    vec<vec<String>> rows;
               vec<String> row = { "xline", /* "CN", */ "genome", "glength", 
                    "alength", "gap", "gapfrac" };
               rows.push_back(row);
               rows.push_back( vec<String>( ) );
               for ( int i = 0; i < all_notes.isize( ); i++ )
               {    if ( i > 0 && all_notes[i].first != all_notes[i-1].first )
                         rows.push_back( vec<String>( ) );
                    rows.push_back( all_notes[i].fourth );    }
               Ofstream( xout, OUTDIR + "/o.genome" );
               xout << "\nWARNING: line numbers below don't make sense\n\n";
               PrintTabular( xout, rows, 3, "lrlrrrr" );    }    }

     // Measure gaps.  Partially outdated....

     double cap_frac = -1;
     if ( genome.size( ) > 0 )
     {    cout << Date( ) << ": set up to measure gaps, saving to o.gaps" << endl;
          int64_t total_lens = 0;
          int total_gaps = 0;
          for ( int i = 0; i < dlines.isize( ); i++ )
          {    total_lens += llens[i] + (K-1);
               const vec<vec<vec<int>>>& L = dlines[i];
               for ( int j = 0; j < L.isize( ); j++ )
               {    if ( !L[j].solo( ) || !L[j][0].empty( ) ) continue;
                    total_gaps++;    }    }
          vec< pair<int,int> > gaps;
          Bool verbose = False;
          MeasureGaps( K, kmers, inv, D, dlines, llens, alignsb, gaps, verbose );
          int64_t gap_sum = 0;
          for ( int i = 0; i < gaps.isize( ); i++ )
               if ( gaps[i].first >= 0 ) gap_sum += gaps[i].first;
          double total_gap_lens
               = gap_sum * double(total_gaps) / double( gaps.size( ) );
          cap_frac = total_gap_lens / total_lens;
          Ofstream( gout, OUTDIR + "/o.gaps" );
          gout << "#  gap_size  gap_edge_id";
          for ( int i = 0; i < gaps.isize( ); i++ )
               gout << gaps[i].first << "  " << gaps[i].second << "\n";    }

     // Compute total edges.

     int64_t total_edges = 0;
     for ( int e = 0; e < D.E( ); e++ )
          total_edges += D.O(e).size( );

     // Find lines of lines, and their lengths, then use this to compute 
     // scaffold N50.  Also computed estimated genome size.

     cout << Date( ) << ": finding lines of lines" << endl;
     vec<vec<vec<vec<int>>>> dlines2x;
     FindLineLines( D, dinv, dlines, dlines2x );
     vec<int> lens2;
     GetLineLineLengths( llens, dlines2x, lens2 );
     ReverseSortSync( lens2, dlines2x );
     {    Ofstream( out, OUTDIR + "/o.linelines" );
          vec<int> x;
          for ( int i = 0; i < dlines2x.isize( ); i++ )
          {    x.clear( );
               const vec<vec<vec<int>>>& M = dlines2x[i];
               for ( int j = 0; j < M.isize( ); j += 2 ) // ONLY PRINTING EVEN!
                    x.push_back( M[j][0][0] );
               out << "M" << i << "[l=" << ToStringAddCommas( lens2[i] )
                    << "]: " << printSeq(x) << "\n";    }    }
     int64_t est_genome_size = 0;
     vec<int> lllens;
     int nlines0 = 0;
     for ( auto x : lens2 ) 
     {    if ( x >= MIN_LINE0 ) nlines0++;
          if ( x >= MIN_LINE ) lllens.push_back(x);
          if ( x >= MIN_LINE ) est_genome_size += x + hb.K( ) - 1;    }
     est_genome_size /= 2;
     int nlines = lllens.size( );

     int64_t scaffold_N50 = N50PL(lllens);
     int64_t scaffold_N60 = NPL(lllens, 0.6);
     ComputeAndWriteHist( lllens, 10000, OUTDIR+"/stats", "scaffold", 
          "final", False );

     cout << Date( ) << ": scaffold N50 = " << ToStringAddCommas(scaffold_N50)
          << endl;
     cout << Date( ) << ": scaffold N60 = " << ToStringAddCommas(scaffold_N60) 
          << endl;
     if (nlines) {
          cout << Date( ) << ": scaffold length-weighted mean length = "
               << ToStringAddCommas( int( round( WeightedMean(lllens) ) ) ) << endl;
     }
     cout << Date( ) << ": assembly size (sum of scaffolds >= 10 kb) = "
          << ToStringAddCommas(est_genome_size) << endl;
     
     // Print large scaffold sizes.  Stupid because it's showing the value for each
     // lineline and its involution, so you get pairwise duplications in the file.

     {    Ofstream( out, OUTDIR + "/stats/large_scaffold_sizes.txt" );
          for ( int i = 0; i < nlines; i++ )
               out << ToStringAddCommas( lllens[i] ) << endl;    }

     // Compute phasetig N50.

     cout << Date( ) << ": computing phasetig N50" << endl;
     vec<int> bblens; // bubble branch lengths
     vec<Bool> marked( D.E( ), False );
     vec<int> to_left, to_right;
     D.ToLeft(to_left), D.ToRight(to_right);
     for ( int i = 0; i < dlines2x.isize( ); i++ )
     {    const vec<vec<vec<int>>>& L = dlines2x[i];
          for ( int j = 0; j < L.isize( ); j++ )
          {    if ( L[j].size( ) != 2 ) continue;
               vec<int> blens( 2, 0 );
               for ( int k = 0; k < 2; k++ )
               {    for ( int l = 0; l < L[j][k].isize( ); l++ )
                    {    int M = L[j][k][l];
                         blens[k] += llens[M];
                         for ( auto d : Contents(dlines[M]) )
                              marked[d] = True;    }    }
               bblens.push_back( Mean(blens) );    }    }
     for ( int i = 0; i < dlines.isize( ); i++ )
     {    const vec<vec<vec<int>>>& L = dlines[i];
          for ( int j = 0; j < L.isize( ); j++ )
          {    if ( L[j].size( ) < 2 ) continue;
               int d1 = L[j][0].front( ), d2 = L[j][0].back( );
               if ( marked[d1] || marked[d2] ) continue;
               int v = to_left[d1], w = to_right[d2];
               if ( D.From(v).size( ) != 2 || D.To(w).size( ) != 2 ) continue;
               int f1 = D.IFrom(v,0), f2 = D.IFrom(v,1);
               vec<int> G1, G2;
               for ( int k = 0; k < L[j].isize( ); k++ )
               {    if ( L[j][k].front( ) == f1 ) G1.push_back( L[j][k].back( ) );
                    if ( L[j][k].front( ) == f2 ) G2.push_back( L[j][k].back( ) );
                         }
               if ( !Contents(G1).solo( ) || !Contents(G2).solo( ) ) continue;
               vec<int> blens;
               for ( int k = 0; k < L[j].isize( ); k++ )
               {    int len = 0;
                    for ( int l = 0; l < L[j][k].isize( ); l++ )
                         len += hb.Kmers( L[j][k][l] );
                    blens.push_back(len);    }
               Sort(blens);
               bblens.push_back( Median(blens) );    }    }
     int64_t phasetig_N50 = ( bblens.size( ) > 0 ? N50(bblens) : 0 );
     ComputeAndWriteHist( bblens, 1000, OUTDIR+"/stats", "phase_block",
          "final", False );

     // Compute N50 perfect stretch.
     
     int64_t N50_perf = -1;
     double errw = -1;
     cout << Date( ) << ": G has size " << G.size( ) << endl;
     if ( G.size( ) > 0 )
     {    MasterVec< SerfVec<triple<int,int,int> > > galignsb;
          galignsb.ReadAll( DIR + "/a.fin.alignsb" );
          vec< vec< vec< triple< vec<int>, align, int > > > > Matches;
          String report;
          AlignFin( hb, inv, D, dinv, dlines, G, galignsb,
               Matches, errw, N50_perf, report, True );
          Ofstream( rout, OUTDIR + "/stats/o.report" );
          rout << report;
     
          StatLogger::log("perf_N50", N50_perf, "N50 perfect stretch size" );
          // Important: the error rate must be multiplied by 1000
          // since it is measured per kb
          StatLogger::log("weighted_error_rate", 1000.0*errw, "weighted error rate" );
     }
     
     // Finished Sequence Evaluation
     
     if ( !G.empty() )
     {
          vec<int> fin (G.size(), vec<int>::IDENTITY);
          EvalVersusFinished( hb, inv, D, dinv, "all", fin, DIR, OUTDIR, False );
     }

     // Reinsert loops.

     ReinsertLoops( hb, inv, D, dinv );

     // Find N50 edge size.

     vec<int> lens;
     for ( int d = 0; d < D.E( ); d++ )
     {    if ( D.O(d)[0] < 0 ) continue;
          int len = 0;
          for ( int i = 0; i < D.O(d).isize( ); i++ )
               len += kmers[ D.O(d)[i] ];
          lens.push_back(len);    }
     ParallelSort(lens);
     int N50_edge = N50(lens);
     ComputeAndWriteHist( lens, 1000, OUTDIR+"/stats", "edge",
          "final", False );

     // Find the barcode distribution histogram
     {
          vec<int64_t> sizes;
          for ( int i = 1; i < bci.isize( ) - 1; i++ )
          {    int n = bci[i+1] - bci[i];
               if ( n > 0 ) sizes.push_back(n);    }
          Sort(sizes);
          
          ComputeAndWriteHist( sizes, int64_t(10), OUTDIR+"/stats",
               "reads_per_barcode", "final", False );
     }
     
     // ============================================================================
     // Log all the output stats
     // ============================================================================

     StatLogger::log( "scaffolds_1kb_plus", nlines0/2, "scaffolds >= 1 kb", true);
     StatLogger::log( "scaffolds_10kb_plus", nlines/2, "scaffolds >= 10 kb", true);

     // Log edge N50.

     StatLogger::log( "edge_N50", N50_edge, "N50 edge length", true);

     // Log contig N50.

     StatLogger::log( "contig_N50", contig_line_N50, "N50 contig length", true);

     // Log phase block N50.

     StatLogger::log( "phase_block_N50", phasetig_N50, "N50 phase block length", true);

     // Log scaffold N50 and N60.  Treated differently.

     StatLogger::log( 
          "scaffold_N50", scaffold_N50, "N50 scaffold length", true);
     StatLogger::log( 
          "scaffold_N60", scaffold_N60, "N60 scaffold length" );

     // Log approximate genome size.

     StatLogger::log( 
          "assembly_size", est_genome_size, "Assembly size", true);
     
     // Log misassembly stats.

     if ( genome.size( ) > 0 ) {    
          // If we have no data to evaluate these error percentages
          // we set them to 100 %
          // misasm_perc is the sum of individual error percents (when not 100%)
          // if we have no info at all then misasm_perc is 100%
          double dis_err_perc = 100.0, ori_err_perc = 100.0, ord_err_perc = 100.0;
          double mis = 0.0;
          bool no_mis_data = true;
          //
          if ( total_err_dis_den > 0 ) {
               dis_err_perc = 100.0*total_err_dis_num/total_err_dis_den;
               mis += double(total_err_dis_num)/double(total_err_dis_den);
               no_mis_data = false;
          } 
          StatLogger::log( "dis_err_perc",
               dis_err_perc,
               "pct of assembly in distant misjoin" );
          //
          if ( total_err_ori_den > 0 ) {
               ori_err_perc = 100.0*total_err_ori_num / total_err_ori_den;
               mis += double(total_err_ori_num)/double(total_err_ori_den);
               no_mis_data = false;
          }
          StatLogger::log( "ori_err_perc",
               ori_err_perc,
               "pct of assembly having wrong orientation" );
          //
          if ( total_err_ord_den > 0 ) {
               ord_err_perc = 100.0*total_err_ord_num/total_err_ord_den;
               mis += double(total_err_ord_num)/double(total_err_ord_den);
               no_mis_data = false;
          }
          StatLogger::log( "ord_err_perc",
               ord_err_perc,
               "pct of assembly out of order" );
          // if we have no data on mis assemblies then set error to 100% 
          if (no_mis_data)
               mis = 1.0;
          StatLogger::log( "misasm_perc", mis*100, "Misassembled pct" );
     }

     // Log gap stats.

     if ( genome.size( ) > 0 ){
          cap_frac = cap_gap_total / double(genome_size);
          StatLogger::log( "air_err_perc", 100*cap_frac,
                    "pct of assembly in captured gaps");
          
          double uncap_frac = (genome_size-cov_total) / double(genome_size);
          StatLogger::log( "vac_err_perc",
                    100*uncap_frac, "pct of assembly in uncaptured gaps");
          
          StatLogger::log("total_gap_perc",
                    (cap_frac+uncap_frac)*100, "Total gap pct");
     }


     // Define some printing functions.

     auto nper = []( double n, double d ) { return 100*n/d; };
     auto rat = [&]( double n, double d )
     {    ostringstream out;
          out << fixed << setprecision(2) << setw(7) << right << n/d;
          return out.str( );    };
     auto hline = [&]( )
     {    out << "----------------------------------------"
               << "----------------------------------------\n";    };
     // {    out << "========================================"
     //           << "====================================\n";    };

     // ============================================================================
     // Print final report.
     // ============================================================================

     // Print sample id and sample description
     // if there's no sample id, then just print DIR/../

     String asmdir = LineOfOutput("readlink -f " + DIR + "/../", true, true );
     asmdir = asmdir.RevAfter( "/" );
     out << "\n";
     hline( );
     out << "SUMMARY" << endl;
     hline( );
     out << "- " << Date( ) << endl;
     String CS_SAMPLE_ID = StatLogger::getStringStat( "CS_SAMPLE_ID", "" );
     String CS_SAMPLE_DESC = StatLogger::getStringStat( "CS_SAMPLE_DESC", "" );

     out << "- [" << (CS_SAMPLE_ID.empty() ? asmdir : CS_SAMPLE_ID)  << "]";
     out << "  " << CS_SAMPLE_DESC << endl;
     
     //if ( WRITE_SUB != "" ) out << "- SUB = " << WRITE_SUB.After( ":" ) << endl;
     
     // Print the code commit hash
     
     String commit_hash = StatLogger::getStringStat( "commit_hash", "" );
     String software_release = StatLogger::getStringStat("software_release","");
#ifdef CS
     out << "- software release = " << software_release << "(" << commit_hash << ")" << endl;
#else
     out << "- software release = " << software_release << endl;
     out << "- commit hash = " << commit_hash << endl;
#endif
     
     // Print final assembly checksum

     StatLogger::log( "checksum", checksum, "Assembly checksum", true);
     out << "- assembly checksum = " << ToStringAddCommas(checksum) << endl;
     
     //---------------------------------------------------------------------------------
     // Automated printing of stats
     //---------------------------------------------------------------------------------
     
     // To add a stat to output
     // add a row with
     // internal stat name, display name, description, units, max exponent suffix
     // For example,
     // {"nreads", "READS","number of reads; ideal 800M-1200M for human", "", "6"},
     // "6" means the highest suffix will be M (not G), empty => 9
     // unit is a one character unit, like b for bases, or % for percent, default:<space>
     // this line will produce the output:
     //
     // - 1200.00 M   = READS          = number of reads; ideal 800-1200 for human
     // 
     // if you get the internal stat name wrong, then you will get a message 
     // on the stdout that the stat could not be found. For example:
     // 
     // Fri Dec 16 18:19:56 2016: could not find a value for nreads
     //
     vec<vec<String>> display={{"INPUT", "", "", "", ""},
{"nreads", "READS","number of reads; ideal 800M-1200M for human", "", "6"},
{"bases_per_read", "MEAN READ LEN","mean read length after trimming; ideal 140", "b", ""},
{"effective_coverage", "EFFECTIVE COV","effective read coverage; ideal ~42 for nominal 56x cov", "x", ""},
{"q30_r2_perc", "READ TWO Q30","fraction of Q30 bases in read 2; ideal 75-85", "%", ""},
{"median_ins_sz", "MEDIAN INSERT","median insert size; ideal 0.35-0.40", "b", ""},
{"proper_pairs_perc", "PROPER PAIRS","fraction of proper read pairs; ideal >= 75", "%", ""},
{"lw_mean_mol_len", "MOLECULE LEN","weighted mean molecule size; ideal 50-100", "b", ""},
{"hetdist", "HETDIST","mean distance between heterozygous SNPs", "b", ""},
{"unbar_perc", "UNBAR","fraction of reads that are not barcoded", "%", ""},
{"rpb_N50", "BARCODE N50","N50 reads per barcode", "", ""},
{"dup_perc", "DUPS","fraction of reads that are duplicates", "%", ""},
{"placed_perc", "PHASED","nonduplicate and phased reads; ideal 45-50", "%", ""},
{"OUTPUT", "", "", "", ""},
{"scaffolds_10kb_plus", "LONG SCAFFOLDS","number of scaffolds >= 10 kb", "", "" },
{"edge_N50", "EDGE N50","N50 edge size", "b", ""},
{"contig_N50", "CONTIG N50","N50 contig size", "b", ""},
{"phase_block_N50", "PHASEBLOCK N50","N50 phase block size", "b", ""},
{"scaffold_N50", "SCAFFOLD N50","N50 scaffold size", "b", ""},
{"scaffold_N60", "SCAFFOLD N60","N60 scaffold size", "b", ""},
{"assembly_size", "ASSEMBLY SIZE","assembly size (only scaffolds >= 10 kb)", "b", ""}};
     
     // These statistics are only computed if a reference genome
     // is supplied
     // currently: human
     if ( genome.size() > 0 ) {
          display.append(vec<vec<String>>{
{"dis_err_perc", "dis error","fraction of assembly in distant misjoin", "%", ""},
{"ori_err_perc","ori error","fraction of assembly having wrong orientation", "%", ""},
{"ord_err_perc","ord error","fraction of assembly out of order", "%", ""},
{"misasm_perc", "MISASSEMBLED", "total mis-assembly fraction", "%", ""},
{"air_err_perc", "air error","fraction of assembly in captured gaps", "%", ""},
{"vac_err_perc", "vac error","fraction of assembly in uncaptured gaps", "%", ""},
{"total_gap_perc", "GAP", "total fraction of assembly in gaps", "%", ""}});
     }

     // These statistics are computed only when finished
     // sequence is available
     // currently: NA12878, HGP, CHM.
     if ( G.size() > 0 ) {
          display.append(vec<vec<String>>{
{"FINISHED SEQUENCE EVALUATION", "", "", "", ""},
{"perf_N50", "PERFECT N50", "N50 perfect stretch size", "b", ""},
{"weighted_error_rate", "ERROR RATE", "weighted error rate per kb", "", ""},
{"tot_edges", "FIN EDGES","number of edges containing finished sequence", "", ""},
{"tot_lines", "FIN LINES","number of lines containing finished sequence", "", ""},
{"tot_uncap", "UNCAPS","number of uncaptured gaps", "", ""},
{"tot_cap", "CAPS","number of captured gaps", "", ""},
{"tot_cov_perc", "FIN COV%","% of finished sequence covered", "%", ""},
{"tot_errs", "FIN ERRORS","number of errors vs finished sequence", "", ""}});
     }
     // Should we display CS stats only?
     // by default false, but is overridden if CS is set during compile time
     Bool show_cs_only = False;
     
     #ifdef CS
          show_cs_only = True;
     #endif
     
     // Go through all the stats and print them to out

     vec<int> max_width (3, 0);
     vec<vec<String>> disp_data;
     for ( auto & row : display ) {
          if ( row[1].empty() ) {
               disp_data.push_back( {row[0]} );
          } else {
               Bool cs = StatLogger::isStatCS( row[0] );
               if ( !show_cs_only || (show_cs_only && cs)) {
                    double val = StatLogger::getNumStat( row[0], -1 );
                    if ( val == -1 ) {
                         cout << Date( ) << ": could not find a value for " << row[0]
                              << endl;
                         continue;
                    }
                    int max = 9;
                    if ( !row[4].empty() )
                         max=row[4].Int();
                    // IMPORTANT: val is printed as xxxx.yy
                    vec<String> disp_row = {PrettyFormat( val, row[3], max ), 
                         row[1], row[2]};
                    disp_data.push_back( disp_row );
                    for ( int i = 0; i < 3; i++ )
                         max_width[i] = Max( max_width[i], (int)disp_row[i].size() );
               }
          }
     }
     String cat = "";
     Bool printed = False;
     for ( auto & row : disp_data ) {
          if ( row.solo() ) {
               cat = row[0];
               printed = False;
          } else {
               if ( ! printed ) {
                    hline( );
                    out << cat << endl;
                    printed = True;
               }
               out << "- ";
               out << row[0];
               out << " = ";
               out << row[1];
               for ( int j = 0; j < max_width[1]-row[1].isize(); j++ )
                    out << " ";
               out << " = ";
               out << row[2];
               out << endl;
          }
     }
     hline( );


     // write only CS facing stats
     StatLogger::dump_csv( OUTDIR + "/summary_cs.csv" );
     // and all the stats here
     StatLogger::dump_json( OUTDIR + "/all_stats.json" );
     StatLogger::dump_json( OUTDIR + "/summary.json", True );
     StatLogger::write( OUTDIR+"/a.perf_stats" );
}
