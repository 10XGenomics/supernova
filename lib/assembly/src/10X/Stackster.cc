// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

// Try to close one gap.

/*

EXAMPLES

====================================================================================

[1]

e1=91229227 e2=56965248  ==>  offset = 211, gap = 61

CAAGCACTTCATAAATTTCTTTTACTTTGGGGTCAATTACTGTAGAGTTTTTATGTTCCTTTGGTGGTATCATATTTCCT
CAAGCACTTCATAAATTTCTTTTACTTTGGGGTCAATTACTGTAGAGTTTTTATGTTCCTTTGGTGGTATCATATTTCCT


(400 matching bases)

                     *
TGTATTTTTAGTAGAGACAGGATTTC
TGTATTTTTAGTAGAGACAGGGTTTC

Agrees with DDN as the only answer.

====================================================================================

[2]

e1=22942409 e2=4067771  ==>  offset = 246, gap = 238

0fw vs 6, 0 mismatches/15 indels (of 628), from 0-628 to 82992895-82993538


AGAAAAAAAAAATAAACATACCAGCTTAAATAACACTGCTACAATCAGCTGAGTGTTTGTGGAAAATCTATTATTTATAT
AGAAAAAAAAAATAAACATACCAGCTTAAATAACACTGCTACAATCAGCTGAGTGTTTGTGGAAAATCTATTATTTATAT

(160 matching bases)
                                         |||||||||||||||
TTAATATGTTTATTAGGCTACAGGTGTCACAGGTCAATTCA               TTTTTTTTTTTTTTTTTTTTTTTT
TTAATATGTTTATTAGGCTACAGGTGTCACAGGTCAATTCATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT


(323 matching bases)

Agrees with DDN as the only answer.

1fw vs 6, 0 mismatches/14 indels (of 647), from 0-647 to 82992877-82993538

TATAAAGCTATCCTTCCAAGAAAAAAAAAATAAACATACCAGCTTAAATAACACTGCTACAATCAGCTGAGTGTTTGTGG
TATAAAGCTATCCTTCCAAGAAAAAAAAAATAAACATACCAGCTTAAATAACACTGCTACAATCAGCTGAGTGTTTGTGG

(160 matching bases)
                                                           ||||||||||||||
CAATACTAAGTACTTTTCTTAATATGTTTATTAGGCTACAGGTGTCACAGGTCAATTCA              TTTTTTT
CAATACTAAGTACTTTTCTTAATATGTTTATTAGGCTACAGGTGTCACAGGTCAATTCATTTTTTTTTTTTTTTTTTTTT

(341 matching bases)

(now have a third)

====================================================================================

[3]

e1=29459996 e2=67092327  ==>  offset = 19, gap = -30

0fw vs 6, 2 mismatches/4 indels (of 419), from 0-419 to 83003415-83003838

GGAGGCTGAGGCAAGAGATTTGCTTGAACCCAGGAGGCAGAGGTTGCGATAAGCCGAGATCACGCCATTGCACTCCAGCC
GGAGGCTGAGGCAAGAGATTTGCTTGAACCCAGGAGGCAGAGGTTGCGATAAGCCGAGATCACGCCATTGCACTCCAGCC

            *              ||||
TGGGTTACAAGAGCAAAACTGCGTCTC    AAAAAAAAAAAAAAAAAAGACAGAAAGAAAGAAATTCAAATTTAAACAA
TGGGTTACAAGACCAAAACTGCGTCTCAAAAAAAAAAAAAAAAAAAAAAGACAGAAAGAAAGAAATTCAAATTTAAACAA

(240 matching bases)

     *
TAGCAGGAGTTAATTACACTCTG
TAGCAAGAGTTAATTACACTCTG

DDN suggests this is exactly right and the only answer.

====================================================================================

[4]

e1=88136837 e2=36566940
e1=88136837 e2=52977845
e1=88136837 e2=80403271
don't work

this does
e1=88136837 e2=30621213  ==>  offset = 300, gap = 23
which uses a 28-kmer e2

0rc vs 7, 0 mismatches/0 indels (of 468), from 0-468 to 82967401-82967869

(perfect match of length 468)

====================================================================================

[5]

e1=20488822 e2=44779836  ==>  offset = 179, gap = 130

0fw vs 7, 0 mismatches/0 indels (of 425), from 0-425 to 82989853-82990278

(perfect match of length 425)

====================================================================================

[6]

e1=152656857 e2=170743295  ==>  offset = 226, gap = 10  or  offfset = 228, gap = 12

0rc vs 12, 0 mismatches/2 indels (of 626), from 0-626 to 83025074-83025702

AGATTTTTTTCTCAGACCCAAAAAGAATTTTAAAATCCATATGAAATAAGAAAGGATTATAAACATCCAAATCAATACTG
AGATTTTTTTCTCAGACCCAAAAAGAATTTTAAAATCCATATGAAATAAGAAAGGATTATAAACATCCAAATCAATACTG

(240 matching bases)
                       ||
CCCATGCAAAAAAAAAAAAAAAA  TGAGCTTGCCTACCTTATACCACACTCAAAAATAAACTTAAAATAGATTAAAGAC
CCCATGCAAAAAAAAAAAAAAAAAATGAGCTTGCCTACCTTATACCACACTCAAAAATAAACTTAAAATAGATTAAAGAC

(228 matching bases)

plus perfect match to hg19

both versions appear to match DDN assembly perfectly

====================================================================================

[7]

e1=173231299 e2=28024325  ==>  offset = 307, gap = 228

28024325rc vs 2, 0 mismatches/0 indels (of 405), from 0-405 to 71475572-71475977

(perfect match of length 405)

173231299rc vs 2, 5 mismatches/0 indels (of 578), from 0-578 to 71476205-71476783 

* ****                                                                          
AAAAAAAAAAAAAAAAAAAAAAATCAGATGCATTTTTTTTTTCATATAATCACCCCAATTTATGCTTTTATGTCTCTCTT
CATCTCAAAAAAAAAAAAAAAAATCAGATGCATTTTTTTTTTCATATAATCACCCCAATTTATGCTTTTATGTCTCTCTT

(498 matching bases)

0rc vs 2, 2 mismatches/0 indels (of 707), from 0-707 to 71475811-71476518

TCAATTCTATGCATGGTGTTCTTCCTCTCATCCTCCAACACTTCACTTAAAAATTACCTTCTTGGGAAGTCCTTCTAAAA
TCAATTCTATGCATGGTGTTCTTCCTCTCATCCTCCAACACTTCACTTAAAAATTACCTTCTTGGGAAGTCCTTCTAAAA

(240 matching bases)
                                                                   *          *
AACCTGGGAGGCAGAGGCTGCAGTGAGCCAAGATCGCACCACTGCACTCTGGCCTGGGTGACAGAGCAAGACTCCATCAC
AACCTGGGAGGCAGAGGCTGCAGTGAGCCAAGATCGCACCACTGCACTCTGGCCTGGGTGACAGAGCGAGACTCCATCTC

(307 matching bases)

The mismatches are puzzling.

====================================================================================

[8]

e1=29564916 e2=72614413  ==>  offset = -6, gap = -44

0fw vs 2, 0 mismatches/1 indels (of 394), from 0-394 to 71512432-71512827

                                                                         |
TGAACCCGGGAGGCGGAGGTTGCAGTGAGCTGAGATCACACCACTGCTGGGTGACACAGTAAGACACTGTCTC AAAAAA
TGAACCCGGGAGGCGGAGGTTGCAGTGAGCTGAGATCACACCACTGCTGGGTGACACAGTAAGACACTGTCTCAAAAAAA

(315 matching bases)

====================================================================================

[9]

e1=69686072 e2=46071345  ==>  offset = 212

0rc vs 3, 0 mismatches/0 indels (of 458), from 411-458 to 102012020-102012067

====================================================================================

*/

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"
#include "Set.h"
#include "feudal/PQVec.h"
#include "paths/RemodelGapTools.h"
#include "paths/long/ReadPath.h"
#include "paths/long/ReadStack.h"
#include "random/Bernoulli.h"
#include "10X/DfTools.h"
#include "10X/Stackster.h"
#include "feudal/SubsetMasterVec.h"

void StackWalk( const basevector& E, const vecbasevector& bases,
     const vecqualvector& quals, vec<Bool>& placed, vec<int>& start, basevector& EE )
{    
     // Heuristics.

     const int K = 48;
     const int MIN_WIN = 60;
     const double MIN_WIN_RATIO = 2.0;

     // Set up computation.

     int n = bases.size( );
     placed.resize( n, False );
     start.resize( n, -1 );
     vec< triple<kmer<K>,int,int> > kmers_plus;
     MakeKmerLookup0( bases, kmers_plus );
     EE.SetToSubOf( E, 0, K );

     // Extend.

     while(1)
     {    
          // Get first kmer.

          kmer<K> x;
          x.SetToSubOf( EE, EE.isize( ) - K );

          // Place more reads.

          int low = LowerBound1( kmers_plus, x );
          int high = UpperBound1( kmers_plus, x );
          Bool duplicate = False;
          for ( int i = low + 1; i < high; i++ )
          {    if ( kmers_plus[i].second == kmers_plus[i-1].second )
               {    duplicate = True;
                    break;    }    }
          if ( !duplicate )
          {    for ( int i = low; i < high; i++ )
               {    int id = kmers_plus[i].second;
                    if ( placed[id] ) continue;
                    int pos = kmers_plus[i].third;
                    placed[id] = True;
                    start[id] = EE.size( ) - K - pos;    }    }

          // Compute next base.

          vec<int> qsum( 4, 0 ), ids( 4, vec<int>::IDENTITY );
          for ( int id = 0; id < n; id++ )
          {    if ( !placed[id] ) continue;
               int pos = EE.isize( ) - start[id];
               if ( pos >= bases[id].isize( ) ) continue;
               qsum[ bases[id][pos] ] += quals[id][pos];    }
          ReverseSortSync( qsum, ids );
          if ( qsum[0] < MIN_WIN ) return;
          if ( double(qsum[0]) < MIN_WIN_RATIO * double(qsum[1]) ) return;
          EE.push_back( ids[0] );    }    }

template <unsigned N>
class PrecomputedBinomialSums
{
public:
    PrecomputedBinomialSums( unsigned w, double p )
    { for ( unsigned n = w; n < N; ++n )
      { double* pVal = mVal[n];
        for ( unsigned k = 0; k != n; ++k )
          *pVal++ = log10(BinomialSum(n,k,p)); } }

    double const* operator[]( unsigned n ) const
    { AssertLt(n,N); return mVal[n]; }

private:
    double mVal[N][N];
};

template< class VB, class VQ, class VP, class VPI > 
void Stackster( const int e1, const int e2, const vec<basevector>& edges,
     VB const& bases, VQ const& quals, const int K, const vec<DataSet>& datasets, 
     const vec<int>& kmers, const vec<int>& inv, const vec<Bool>& dup,
     VP const& paths, VPI const& paths_index, vec<basevector>& closures, 
     vec<int>& trim, const int VERBOSITY, const Bool ALT, const Bool EXP,
     const vec< pair<int64_t,Bool> >& idsfw2 )
{
     // Create the read set, including partners that are off the ends.

     if ( VERBOSITY >= 1 ) cout << Date( ) << ": creating read set" << endl;
     const int MAX_DIST = 800;
     const int MAX_PROX = 400;
     vec< pair<int64_t,Bool> > idsfw;
     for ( int spass = 1; spass <= 2; spass++ )
     {    int e = ( spass == 1 ? e1 : inv[e2] );
          if ( e < 0 ) break;
          for ( int pass = 1; pass <= 2; pass++ )
          {    const int f = ( pass == 1 ? e : inv[e] );
               for ( int i = 0; i < (int) paths_index[f].size( ); i++ )
               {    int64_t id = paths_index[f][i];
                    if ( dup[id/2] ) continue;
                    const ReadPath& p = paths[id];
                    int N = edges[spass-1].size( );
                    int start = p.getOffset( );
                    if ( pass == 2 ) start = N - start - bases[id].isize( );
                    int stop = start + bases[id].isize( );
                    if ( ALT && N - stop > MAX_PROX ) continue;
                    idsfw.push( id, spass == 1 ^ pass == 2 );   }   }
          for ( int i = 0; i < (int) paths_index[e].size( ); i++ )
          {    int64_t id1 = paths_index[e][i];
               ForceAssertLt(id1/2, dup.size() );
               if ( dup[id1/2] ) continue;
               int64_t id2 = ( id1 % 2 == 0 ? id1+1 : id1-1 );
               const ReadPath& p1 = paths[id1];
               for ( int j = 0; j < (int) p1.size( ); j++ )
               {    if ( p1[j] == e )
                    {    int pos1 = p1.getOffset( );
                         for ( int l = 0; l < j; l++ )
                              pos1 -= kmers[ p1[l] ];
                         if ( kmers[e] + K - 1 - pos1 > MAX_DIST ) continue;
                         idsfw.push( id2, spass == 2 );    }    }    }    }
     idsfw.append(idsfw2);
     UniqueSort(idsfw);

     // Load reads.

     if ( VERBOSITY >= 1 ) cout << Date( ) << ": loading reads" << endl;
     vec<int64_t> ids;
     vecbasevector basesx;
     VecPQVec qualsx;
     for ( int i = 0; i < idsfw.isize( ); i++ )
     {    int64_t id = idsfw[i].first;
          ids.push_back( id );
          basesx.push_back( bases[id] );
          qualsx.push_back( quals[id] );    }

     // Unpack data.

     if ( VERBOSITY >= 1 ) cout << Date( ) << ": unpacking" << endl;
     vecbasevector basesy(basesx);
     vecqualvector qualsy( ids.size( ) );
     for ( int i = 0; i < (int) basesx.size( ); i++ )
          qualsx[i].unpack( &qualsy[i] );

     // Lower quality scores in long homopolymers.

     const int HSTAR = 5;
     if ( VERBOSITY >= 1 ) cout << Date( ) << ": lowering" << endl;
     for ( int i = 0; i < (int) basesy.size( ); i++ )
     {    qualvector& q = qualsy[i];
          const basevector& b = basesy[i];
          for ( int j = 0; j < (int) q.size( ); j++ )
          {    int k;
               for ( k = j + 1; k < (int) q.size( ); k++ )
                    if ( b[k] != b[j] ) break;
               if ( k - j >= 20 )
               {    for ( int l = j; l < k; l++ )
                         q[l] = Min( (int) q[l], HSTAR );    }
               j = k - 1;    }    }

     // Build read stacks.

     if ( VERBOSITY >= 1 ) 
          cout << "\n" << Date( ) << ": building stacks" << endl;
     vec<readstack> stacks(2);
     vec<vec<Bool>> placed(2);
     for ( int i = 0; i < (int) basesx.size( ); i++ )
     {    if ( !idsfw[i].second )
          {    basesy[i].ReverseComplement( );
               qualsy[i].ReverseMe( );    }    }
     trim.clear( );
     trim.resize( 2, 0 );
     vec<basevector> EE(2);
     for ( int spass = 1; spass <= 2; spass++ )
     {    int e = ( spass == 1 ? e1 : e2 );
          if ( spass == 2 )
          {    for ( int i = 0; i < (int) basesy.size( ); i++ )
               {    basesy[i].ReverseComplement( );
                    qualsy[i].ReverseMe( );    }    }

          // Walk the stack.

          vec<int> start;
          basevector E = edges[spass-1];
          if ( spass == 2 ) E.ReverseComplement( );
          if (ALT)
          {    if ( E.isize( ) > MAX_PROX )
               {    trim[spass-1] = E.isize( ) - MAX_PROX;
                    E.SetToSubOf( E, E.isize( ) - MAX_PROX, MAX_PROX );    }    }
          StackWalk( E, basesy, qualsy, placed[spass-1], start, EE[spass-1] );
          if ( VERBOSITY >= 2 ) 
          {    if ( spass == 1 ) cout << endl;
               EE[spass-1].Print( cout, "walk_" + ToString(spass) );    }

          // Convert to stack.

          vec< triple<int64_t,int,Bool> > locs;
          vec<int> src;
          for ( int i = 0; i < (int) basesy.size( ); i++ )
          {    if ( !placed[spass-1][i] ) continue;
               src.push_back(i);
               locs.push( ids[i], start[i], idsfw[i].second ^ spass == 2 );    }
          if ( locs.empty( ) )
          {    if ( VERBOSITY >= 1 )
               {    cout << "\nNo reads placed for " << e << endl;
                    cout << "Giving up." << endl;    }
               return;    }
          int M = 0;
          for ( int i = 0; i < locs.isize( ); i++ )
          {    int64_t id = locs[i].first;
               int pos = locs[i].second;
               Bool fw = locs[i].third;
               M = Max( M, pos + basesy[ src[i] ].isize( ) );    }
          readstack& stack = stacks[spass-1];
          stack.Initialize( locs.size( ), M );
          for ( int i = 0; i < locs.isize( ); i++ )
          {    int64_t id = locs[i].first;
               const qualvector& q = qualsy[ src[i] ];
               const basevector& b = basesy[ src[i] ];
               int pos = locs[i].second;
               Bool fw = locs[i].third;
               stack.SetId( i, id );
               stack.SetLen( i, b.size( ) );
               stack.SetRc2( i, !fw );
               stack.SetOffset( i, pos );
               for ( int j = 0; j < b.isize( ); j++ )
               {    int p = pos + j;
                    if ( p >= 0 && p < stack.Cols( ) )
                    {    stack.SetBase( i, p, b[j] );
                         stack.SetQual( i, p, q[j] );    }    }    }

          // Trim stack.

          if ( !ALT )
          {    const int MAX_WIDTH = 400;
               if ( stack.Cols( ) > MAX_WIDTH )
               {    trim[spass-1] = stack.Cols( ) - MAX_WIDTH;
                    stack.Trim( 
                         stack.Cols( ) - MAX_WIDTH, stack.Cols( ) );    }    }    }

     // Orient and print stacks.

     for ( int spass = 1; spass <= 2; spass++ )
     {    readstack& stack = stacks[spass-1];
          if ( spass == 2 ) stack.Reverse( );
          if ( VERBOSITY >= 2 )
          {    vec<String> title( stack.Rows( ) );
               for ( int i = 0; i < stack.Rows( ); i++ )
               {    int64_t id = stack.Id(i);
                    int di;
                    for ( di = 0; di < datasets.isize( ); di++ )
                         if ( id < datasets[di].start ) break;
                    const ReadDataType& dtype = datasets[di-1].dt;
                    if ( dtype == ReadDataType::BAR_10X ) title[i] = "10X";
                    else if ( dtype == ReadDataType::UNBAR_10X ) title[i] = "10X";
                    else if ( dtype == ReadDataType::PCR_FREE ) title[i] = "free";
                    else if ( dtype == ReadDataType::PCR ) title[i] = "PCR";    }
               stack.Print( cout, 69, title );    }    }

     // Look for possible matches between reads in one stack and the consensus in
     // the opposite stack.

     const int MM = 8;
     basevector bb(MM);
     int n1 = stacks[0].Cols( ), n2 = stacks[1].Cols( );
     vec<Bool> allowed( n1 + n2, False ), disallowed( n1 + n2, False );
     const int QUALCAP = 100;
     basevector con1x, con2x;
     qualvector conq1x, conq2x;
     stacks[0].Consensus1( con1x, conq1x, QUALCAP );
     stacks[1].Consensus1( con2x, conq2x, QUALCAP );
     const int MIN_FLAT = 10;
     for ( int i = 0; i < con1x.isize( ); i++ )
     {    int j;
          for ( j = i + 1; j < con1x.isize( ); j++ )
               if ( con1x[j] != con1x[i] ) break;
          if ( j - i >= MIN_FLAT )
          {    for ( int k = i; k < j; k++ )
                    conq1x[k] = 0;    }
          i = j - 1;    }
     for ( int i = 0; i < con2x.isize( ); i++ )
     {    int j;
          for ( j = i + 1; j < con2x.isize( ); j++ )
               if ( con2x[j] != con2x[i] ) break;
          if ( j - i >= MIN_FLAT )
          {    for ( int k = i; k < j; k++ )
                    conq2x[k] = 0;    }
          i = j - 1;    }
     vec< pair<double,int> > results;
     for ( int s = 0; s < 2; s++ )
     {    const readstack &s1 = stacks[s], &s2 = stacks[1-s];
          const basevector& con1 = ( s == 0 ? con1x : con2x );
          set< pair<int,int> > sees;

          // Find the kmers on con1 that are acceptable.  First half can't be a
          // homopolymer and kmer must contain at least three different bases.

          vec< pair<basevector,int> > bbs;
          for ( int j = 0; j <= s1.Cols( ) - MM; j++ )
          {    for ( int m = 0; m < MM; m++ )
                    bb.Set( m, con1[j+m] );
               Bool diff = False;
               for ( int m = 1; m < MM/2; m++ )
                    if ( bb[m] != bb[m-1] ) diff = True;
               if ( !diff ) continue;
               diff = False;
               for ( int m = 1; m < MM/2; m++ )
                    if ( bb[MM/2+m] != bb[MM/2+m-1] ) diff = True;
               if ( !diff ) continue;
               vec<Bool> present(4,False);
               for ( int m = 0; m < MM; m++ )
                    present[ bb[m] ] = True;
               const int MIN_BASES = 3;
               if ( Sum(present) < MIN_BASES ) continue;
               bbs.push( bb, j );    }
          Sort(bbs);

          // Look at every read kmer.

          for ( int r = 0; r < s2.Rows( ); r++ )
          {    for ( int l = 0; l <= s2.Cols( ) - MM; l++ )
               {    
                    // Is kmer defined?

                    Bool defined = True;
                    for ( int m = 0; m < MM; m++ )
                    {    if ( !s2.Def( r, l+m ) )
                         {    defined = False;
                              break;    }    }
                    if ( !defined ) continue;

                    // Find matches of kmer to con1.

                    basevector bb2(MM);
                    for ( int m = 0; m < MM; m++ )
                         bb2.Set( m, s2.Base( r, l+m ) );
                    int low = LowerBound1( bbs, bb2 );
                    int high = UpperBound1( bbs, bb2 );
                    for ( int v = low; v < high; v++ )
                    // for ( int v = 0; v < bbs.isize( ); v++ )
                    {    int j = bbs[v].second;
                         int offset = ( s == 0 ? j - l : l - j );
                         if ( Member( sees, make_pair( r, offset ) ) ) continue;

                         // Test for inconsistent columns.

                         if ( disallowed[ offset + n2 ] ) continue;
                         if ( !allowed[ offset + n2 ] )
                         {    int conflicts = 0;
                              for ( int i1 = 0; i1 < con1x.isize( ); i1++ )
                              {    int i2 = i1 - offset;
                                   if ( i2 < 0 || i2 >= con2x.isize( ) ) continue;
                                   if ( conq1x[i1] == QUALCAP 
                                        && conq2x[i2] == QUALCAP 
                                        && con1x[i1] != con2x[i2] )
                                   {    conflicts++;    }    }
                              const int CONFLICTS = 10;
                              if ( conflicts < CONFLICTS )
                                   allowed[ offset + n2 ] = True;
                              else
                              {    disallowed[ offset + n2 ] = True;
                                   continue;    }    }

                         sees.insert( make_pair( r, offset ) );

                         int implied = offset + trim[0]
                              + stacks[1].Cols( ) - ( kmers[e2] + K - 1 - trim[1] );
                         int sep = implied - ( kmers[e1] + K - 1 );

                         // Compute bits (v1).

                         const int WID = 20;
                         const int MAX_OVERLAP = 1000;
                         static PrecomputedBinomialSums<MAX_OVERLAP> gBS(WID,0.75);
                         int offx = ( s == 0 ? offset : -offset );
                         int pos1 = ( offx >= 0 ? offx : 0 ); 
                         int pos2 = ( offx <= 0 ? -offx : 0 );
                         int overlap = IntervalOverlap( 
                              0, s1.Cols( ), offx, offx + s2.Cols( ) );
                         double minp = 0;

                         vec<int> err_sum( overlap + 1, 0 );
                         for ( int v = 1; v <= overlap; v++ )
                         {    err_sum[v] = err_sum[v-1];
                              if ( con1[pos1+v-1] != s2.Base( r, pos2+v-1 ) ) 
                                   err_sum[v]++;    }

                         for ( int start = 0; start < overlap; start++ )
                         {    for ( int n = WID; n <= overlap - start; n++ )
                              {    if ( !s2.Def( r, pos2+start )
                                        || !s2.Def( r, pos2+start+n-1 ) )
                                   {    continue;     }
                                   int k = err_sum[start+n] - err_sum[start];
                                   minp = Min( minp, gBS[n][k] );    }    }
                         double bits = -minp * 10.0 / 6.0;
                         const double MIN_BITS = 20.0;
                         if ( bits < MIN_BITS ) continue;

                         // Save offset.

                         results.push( bits, offset );

                         // Print.

                         if ( VERBOSITY >= 2 )
                         {    cout << endl << "match " << bb.ToString( ) 
                                   << " of " << s2.Id(r) << " to consensus of " 
                                   << ( s == 0 ? "1st" : "2nd" )
                                   << " stack, offset = " << offset 
                                   << ", bits = " << bits << endl;    
                              cout << "implied gap between edges = " << sep << endl;
                              const int flank = 35;
                              String q1, q2;
                              for ( int m = j-flank; m < j+MM+flank; m++ )
                              {    if ( m < 0 || m >= s1.Cols( ) ) 
                                        q1.push_back( ' ' );
                                   else q1.push_back( as_base( con1[m] ) );    }
                              for ( int m = l-flank; m < l+MM+flank; m++ )
                              {    if ( m < 0 || m >= s2.Cols( ) ) 
                                        q2.push_back( ' ' );
                                   else q2.push_back( s2.BaseVis( r, m ) );    }
                              for ( int j = 0; j < q1.isize( ); j++ )
                              {    if ( q1[j] != ' ' && q2[j] != ' ' 
                                        && q1[j] != q2[j] )
                                   {    cout << "*";     }
                                   else cout << " ";    }
                              cout << endl << q1 << endl << q2 << endl;
                                   }    }    }    }    }

     // Filter offsets.  Note that we could require that an offset occur at least
     // twice.

     ReverseSort(results);
     const double MAX_DIFF = 20.0;
     vec<int> offsets;
     for ( int i = 0; i < results.isize( ); i++ )
     {    double score = results[i].first;
          int offset = results[i].second;
          if ( offsets.empty( ) ) offsets.push_back(offset);
          else if ( Member( offsets, offset ) ) continue;
          else if ( results[0].first - score <= MAX_DIFF ) 
          {    offsets.push_back(offset);    }    }
     Sort(offsets);
     if ( VERBOSITY >= 1 ) cout << "\nusing offsets = " << printSeq(offsets) << endl;

     // Enrich the stacks by adding back reads that were originally placed
     // on the edges.  Could do this earlier.  And should understand why they're
     // not being adding in originally.

     /*
     vec<int> start(2); // start of edge on stack
     vec<Bool> started(2,False);
     for ( int s = 0; s < s; s++ )
     {    basevector con;
          qualvector conq;
          stacks[s].Consensus1( con, conq );
          int e = ( s == 0 ? e1 : e2 );
          for ( int i = 0; i <= con.isize( ) - K; i++ )
          for ( int j = 0; j <= kmers[e] + K - 1 - K; j++ )
          {    Bool match = True;
               for ( int l = 0; l < K; l++ )
               {    if ( edges[s][j+l] != con[i+l] )
                    {    match = False;
                         break;    }    }
               if (match)
               {    start[s] = i - j;
                    started[s] = True;
                    goto done;    }    }
          done: continue;    }
     for ( int s = 0; s < 2; s++ )
     {    readstack& stack = stacks[s];
          vec<int64_t> sids;
          for ( int i = 0; i < stack.Rows( ); i++ )
               sids.push_back( stack.Id(i) );
          UniqueSort(sids);
          int e = ( s == 0 ? e1 : e2 );
          for ( int pass = 1; pass <= 2; pass++ )
          {    const int f = ( pass == 1 ? e : inv[e] );
               for ( int i = 0; i < (int) paths_index[f].size( ); i++ )
               {    int64_t id = paths_index[f][i];
                    int pp = BinPosition( ids, id ); // not quite right, fw/rc both?
                    if ( pp < 0 ) continue;
                    if ( BinMember( sids, id ) ) continue;
                    const ReadPath& p = paths[id];
                    int rpos = p.getOffset( ); // start of read on edge
                    if ( pass == 2 )
                         rpos = kmers[e] + K - 1 - rpos - basesx[pp].isize( );
                    int spos = rpos + start[s]; // start of read on stack
                    if ( IntervalOverlap( 0, stacks[s].Cols( ), 
                         spos, spos + basesx[pp].isize( ) ) <= 0 ) 
                    {    continue;    }
                    cout << "looks like read " << id << " could be placed in "
                         << ( s == 0 ? "1st" : "2nd" ) << " stack at position " 
                         << spos << endl;    
                    int N = stacks[s].Rows( );
                    stacks[s].AddRows(1);
                    const basevector& b = basesy[pp];
                    const qualvector& q = qualsy[pp];
                    stacks[s].SetId( N, id );
                    stacks[s].SetLen( N, basesy[pp].size( ) );
                    stacks[s].SetRc2( N, pass == 2 );
                    stacks[s].SetOffset( N, spos );
                    for ( int j = 0; j < b.isize( ); j++ )
                    {    int p = spos + j;
                         if ( p >= 0 && p < stacks[s].Cols( ) )
                         {    stacks[s].SetBase( N, p, b[j] );
                              stacks[s].SetQual( N, p, q[j] );    }    }
                                   }    }    }
     */

     // Merge stacks.

     const int MAX_OFFSETS = 3;
     if ( offsets.isize( ) <= MAX_OFFSETS )
     {    for ( int j = 0; j < offsets.isize( ); j++ )
          {    readstack s1( stacks[0] ), s2( stacks[1] );
               int offset = offsets[j];

               // Trim stacks to prevent extension beyond the original edges.

               int o = offset;
               int right_ext = s1.Cols( ) - offset - s2.Cols( );
               if ( right_ext > 0 ) s1.Trim( 0, s1.Cols( ) - right_ext );
               int left_ext = -offset;
               if ( left_ext > 0 )
               {    s2.Trim( left_ext, s2.Cols( ) );
                    o = 0;    }

               // Now merge the stacks.

               readstack s(s1);
               s.Merge( s2, o );
               if ( VERBOSITY >= 2 )
               {    cout << "\njoint stack for offset " << offset << ":\n";
                    vec<String> title( s.Rows( ) );
                    for ( int i = 0; i < s.Rows( ); i++ )
                    {    int64_t id = s.Id(i);
                         int di;
                         for ( di = 0; di < datasets.isize( ); di++ )
                              if ( id < datasets[di].start ) break;
                         const ReadDataType& dtype = datasets[di-1].dt;
                         if ( dtype == ReadDataType::BAR_10X ) title[i] = "10X";
                         else if ( dtype == ReadDataType::UNBAR_10X ) 
                              title[i] = "10X";
                         else if ( dtype == ReadDataType::PCR_FREE ) 
                              title[i] = "free";
                         else if ( dtype == ReadDataType::PCR ) 
                              title[i] = "PCR";    }
                    s.Print( cout, 69, title );    }
               basevector f;
               qualvector fragq;
               s.Consensus1( f, fragq );
               if ( VERBOSITY >= 2 )
               {    cout << endl;
                    f.Print( cout, "consensus" );     }
               closures.push_back(f);    }    }

     // Print stack consensuses.

     if ( VERBOSITY >= 2 )
     {    cout << "\nstack consensuses:\n";
          for ( int j = 0; j < 2; j++ )
          {    basevector f;
               qualvector q;
               stacks[j].Consensus1(f,q);
               f.Print( cout, ( j == 0 ? "left" : "right" ) );    }    }

// =================================================================================

     // Look directly for bridges.

     if ( !EXP ) return;
     cout << "\n" << Date( ) << ": start experimental approach" << endl;
     EE[1].ReverseComplement( );
     for ( int i = 0; i < (int) basesy.size( ); i++ )
     {    basesy[i].ReverseComplement( );
          qualsy[i].ReverseMe( );    }
     vecbasevector all(basesy);
     all.push_back( EE[0] );
     all.push_back( EE[1] );
     vec< triple<kmer<MM>,int,int> > kmers_plus;
     MakeKmerLookup0( all, kmers_plus );
     vec< triple<int,int,int> > matches;
     for ( int i = 0; i < kmers_plus.isize( ); i++ )
     {    int j;
          for ( j = i + 1; j < kmers_plus.isize( ); j++ )
               if ( kmers_plus[j].first != kmers_plus[i].first ) break;
          vec<int> count(2,0);
          for ( int m = 0; m < 2; m++ )
          for ( int k = i; k < j; k++ )
               if ( kmers_plus[k].second == (int) basesy.size( ) + m ) count[m]++;
          if ( Sum(count) == 0 )
          {    i = j - 1;
               continue;    }
           
          // Check kmer for acceptable.  First half can't be a homopolymer 
          // and kmer must contain at least three different bases.

          basevector bb;
          kmers_plus[i].first.GetBasevector(bb);
          Bool good = False, diff = False;
          for ( int m = 1; m < MM/2; m++ ) if ( bb[m] != bb[m-1] ) diff = True;
          if (diff)
          {    diff = False;
               for ( int m = 1; m < MM/2; m++ )
                    if ( bb[MM/2+m] != bb[MM/2+m-1] ) diff = True;
               if (diff)
               {    vec<Bool> present(4,False);
                    for ( int m = 0; m < MM; m++ ) present[ bb[m] ] = True;
                    const int MIN_BASES = 3;
                    if ( Sum(present) >= MIN_BASES ) good = True;    }    }
          if ( !good )
          {    i = j - 1;   
               continue;    }
          for ( int k1 = i; k1 < j; k1++ )
          for ( int k2 = i; k2 < j; k2++ )
          {    if ( kmers_plus[k1].second >= (int) basesy.size( ) ) continue;
               int m = kmers_plus[k2].second - (int) basesy.size( );
               if ( m < 0 ) continue;
               int cpos = kmers_plus[k2].third;
               int rid = kmers_plus[k1].second, rpos = kmers_plus[k1].third;

               // Validate.  Require an overlap of size <= K, having at least K/2
               // bases in agreement, and only test the three cases where the MM-mer
               // is on the left end, in the middle, or on the right end.

               Bool valid = False;
               const basevector &r = basesy[rid], &c = EE[m];
               int match = MM;
               for ( int l = 0; l < K - MM; l++ )
               {    if ( rpos+MM+l >= r.isize( ) || cpos+MM+l >= c.isize( ) ) break;
                    if ( r[rpos+MM+l] == c[cpos+MM+l] ) match++;    }
               if ( match >= K/2 ) valid = True;
               match = MM;
               for ( int l = 1; l <= K - MM; l++ )
               {    if ( rpos-l <= 0 || cpos-l <= 0 ) break;
                    if ( r[rpos-l] == c[cpos-l] ) match++;    }
               if ( match >= K/2 ) valid = True;
               match = MM;
               for ( int l = 0; l < (K-MM)/2; l++ )
               {    if ( rpos+MM+l >= r.isize( ) || cpos+MM+l >= c.isize( ) ) break;
                    if ( r[rpos+MM+l] == c[cpos+MM+l] ) match++;    }
               for ( int l = 1; l <= (K-MM)/2; l++ )
               {    if ( rpos-l <= 0 || cpos-l <= 0 ) break;
                    if ( r[rpos-l] == c[cpos-l] ) match++;    }
               if ( match >= K/2 ) valid = True;
               if ( !valid ) continue;

               // Record overlap.

               int offset = cpos - rpos;
               matches.push( rid, m, offset );    }
          i = j - 1;    }

     // Screen matches.

     UniqueSort(matches);
     vec< quad<int,int,int,double> > matches2;
     for ( int i = 0; i < matches.isize( ); i++ )
     {    int rid = matches[i].first, m = matches[i].second; 
          int offx = matches[i].third;

          // Compute bits (v1).

          const int WID = 20;
          const int MAX_OVERLAP = 1000;
          static PrecomputedBinomialSums<MAX_OVERLAP> gBS(WID,0.75);
          int pos1 = ( offx >= 0 ? offx : 0 ), pos2 = ( offx <= 0 ? -offx : 0 );
          int overlap = IntervalOverlap( 
               0, EE[m].isize( ), offx, offx + basesy[rid].isize( ) );
          double minp = 0;
          vec<int> err_sum( overlap + 1, 0 );
          for ( int v = 1; v <= overlap; v++ )
          {    err_sum[v] = err_sum[v-1];
               if ( EE[m][pos1+v-1] != basesy[rid][pos2+v-1] ) err_sum[v]++;    }
          for ( int start = 0; start < overlap; start++ )
          {    for ( int n = WID; n <= overlap - start; n++ )
               {    int k = err_sum[start+n] - err_sum[start];
                    minp = Min( minp, gBS[n][k] );    }    }
          double bits = -minp * 10.0 / 6.0;
          const double MIN_BITS = 20.0;
          if ( bits >= MIN_BITS ) matches2.push( rid, m, offx, bits );    }

     // Find the bridges.

     PRINT( matches2.size( ) );
     int bridges = 0;
     vec< triple<int,double,double> > obb;
     for ( int i = 0; i < matches2.isize( ); i++ )
     {    int j;
          for ( j = i + 1; j < matches2.isize( ); j++ )
               if ( matches2[j].first != matches2[i].first ) break;
          for ( int k1 = i; k1 < j; k1++ )
          for ( int k2 = i; k2 < j; k2++ )
          {    if ( matches2[k1].second != 0 || matches2[k2].second != 1 ) continue;
               bridges++;
               int off = matches2[k1].third - matches2[k2].third;
               double bits1 = matches2[k1].fourth, bits2 = matches2[k2].fourth;
               int offorig 
                    = off - trim[0] - trim[1] - stacks[1].Cols( ) + EE[1].isize( );
               obb.push( offorig, bits1, bits2 );    }
          i = j - 1;    }

     // Filter.

     Sort(obb);
     double best = -1;
     for ( int i = 0; i < obb.isize( ); i++ )
     {    double b1 = obb[i].second, b2 = obb[i].third;
          best = Max( best, Min( b1, b2 ) );    }
     vec<Bool> to_delete( obb.size( ), False );
     for ( int i = 0; i < obb.isize( ); i++ )
     {    int j;
          for ( j = i + 1; j < obb.isize( ); j++ )
               if ( obb[j].first != obb[i].first ) break;
          if ( j - i == 1 )
          {    int m = Min( obb[i].second, obb[i].third );
               if ( m <= 35 && best - m >= 20 ) to_delete[i] = True;    }
          i = j - 1;    }
     EraseIf( obb, to_delete );

     for ( int i = 0; i < obb.isize( ); i++ )
     {    int off = obb[i].first;
          double b1 = obb[i].second, b2 = obb[i].third;
          PRINT3( off, b1, b2 );    }
     cout << Date( ) << ": experimental approach complete" << endl;    }

template void Stackster( const int e1, const int e2, const vec<basevector>& edges,
     VirtualMasterVec<basevector> const& bases, VirtualMasterVec<PQVec> const& quals,
     const int K, const vec<DataSet>& datasets, 
     const vec<int>& kmers, const vec<int>& inv, const vec<Bool>& dup,
     VirtualMasterVec<ReadPath> const& paths, 
     VirtualMasterVec<ULongVec> const& paths_index, vec<basevector>& closures, 
     vec<int>& trim, const int VERBOSITY, const Bool ALT, const Bool EXP,
     const vec< pair<int64_t,Bool> >& idsfw2 );

template void Stackster( const int e1, const int e2, const vec<basevector>& edges,
     MasterVec<basevector> const& bases, MasterVec<PQVec> const& quals,
     const int K, const vec<DataSet>& datasets, 
     const vec<int>& kmers, const vec<int>& inv, const vec<Bool>& dup,
     MasterVec<ReadPath> const& paths, MasterVec<ULongVec> const& paths_index,
     vec<basevector>& closures, vec<int>& trim, const int VERBOSITY, 
     const Bool ALT, const Bool EXP,
     const vec< pair<int64_t,Bool> >& idsfw2 );

template void Stackster( const int e1, const int e2, const vec<basevector>& edges,
     MasterVec<basevector> const& bases, SubsetMasterVec<PQVec> const& quals,
     const int K, const vec<DataSet>& datasets, 
     const vec<int>& kmers, const vec<int>& inv, const vec<Bool>& dup,
     SubsetMasterVec<ReadPath> const& paths, SubsetMasterVec<ULongVec> const& paths_index,
     vec<basevector>& closures, vec<int>& trim, const int VERBOSITY, 
     const Bool ALT, const Bool EXP, const vec< pair<int64_t,Bool> >& idsfw2 );
