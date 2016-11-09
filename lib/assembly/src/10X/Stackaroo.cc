// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

// Try to close gaps in supergraph assembly.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"
#include "FastIfstream.h"
#include "PrintAlignment.h"
#include "feudal/PQVec.h"
#include "math/HoInterval.h"
#include "paths/RemodelGapTools.h"
#include "paths/long/MakeKmerStuff.h"
#include "paths/long/ReadPath.h"
#include "paths/long/ReadStack.h"
#include "paths/long/large/GapToyTools.h"
#include "random/Bernoulli.h"
#include "random/Random.h"
#include "10X/DfTools.h"
#include "10X/Gap.h"
#include "10X/Stackaroo.h"
#include "10X/Super.h"

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

// #define MAKE_KMER_LOOKUP MakeKmerLookup0
#define MAKE_KMER_LOOKUP MakeKmerLookup3

void GetSupport( const HyperBasevectorX& hb, const HyperBasevectorX& hbl2,
     const vecbasevector& basesy, vec<vec<int>>& support )
{
     support.clear( );
     support.resize( hbl2.E( ) );
     {    const int K = 48;
          ForceAssertEq( K, hb.K( ) );
          vec< triple<kmer<K>,int,int> > kmers_plus;
          vecbasevector allb = basesy;
          for ( int e = 0; e < hbl2.E( ); e++ )
               allb.push_back( hbl2.O(e) );
          MAKE_KMER_LOOKUP( allb, kmers_plus );
          for ( int i = 0; i < kmers_plus.isize( ); i++ )
          {    int j;
               for ( j = i + 1; j < kmers_plus.isize( ); j++ )
                    if ( kmers_plus[j].first != kmers_plus[i].first ) break;

               // Don't use kmers ending in 10 or more of the same base.

               basevector b;
               kmers_plus[i].first.GetBasevector(b);
               Bool bad1 = True, bad2 = True;
               for ( int l = 1; l < 10; l++ )
               {    if ( b[l] != b[l-1] ) bad1 = False;
                    if ( b[ K - l - 1 ] != b[ K - l ] ) bad2 = False;    }

               if ( !bad1 && !bad2 )
               {    for ( int k1 = i; k1 < j; k1++ )
                    for ( int k2 = i; k2 < j; k2++ )
                    {    int id1 = kmers_plus[k1].second;
                         if ( id1 >= (int) basesy.size( ) ) continue;
                         int id2 = kmers_plus[k2].second - (int) basesy.size( );
                         if ( id2 < 0 ) continue;
                         support[id2].push_back(id1);    }    }
               i = j - 1;    }
          for ( int e = 0; e < hbl2.E( ); e++ )
               UniqueSort( support[e] );    }    }

void Stackaroo( 

     // inputs:

     VirtualMasterVec<basevector> bases,
     VirtualMasterVec<PQVec> quals, const HyperBasevectorX& hb,
     const vec<int>& inv, VirtualMasterVec<ReadPath> xpaths,
     VirtualMasterVec<ULongVec> xpaths_index, 

     // inputs and outputs:

     digraphE<vec<int>>& D, vec<int>& dinv, 

     // control over gap set:

     const String& S, String& R, const Bool ALL,

     // logging:

     const Bool VERBOSE, const int VERBOSITY, const HyperBasevectorX& ddn,
     const Bool VISUAL_ABBR, const Bool DIRECT, const Bool SHOW_MERGED_STACKS
     )
{
     // Start!

     double clock = WallClockTime( );

     // Parse R and S.

     int opts = 0;
     if ( S != "" ) opts++;
     if ( R != "" ) opts++;
     if (ALL) opts++;
     if ( opts != 1 )
     {    cout << "Exactly one of S and R and ALL must be specified." 
               << endl << endl;
          Scram(1);    }
     vec<int> rs, ss;
     if ( R != "" && isdigit( R[0] ) ) R = "{" + R + "}";
     ParseIntSet( R, rs, False );
     ParseIntSet( S, ss, False );

     // Create ancillary data structures.

     vec<int> to_left, to_right;
     D.ToLeft(to_left), D.ToRight(to_right);
     int K = hb.K( );
     vec<int> dkmers( D.E( ), 0 );
     #pragma omp parallel for
     for ( int d = 0; d < D.E( ); d++ )
     {    if ( D.O(d)[0] < 0 ) continue;
          for ( int j = 0; j < D.O(d).isize( ); j++ )
               dkmers[d] += hb.Kmers( D.O(d)[j] );    }

     // Find the gaps.

     vec<int> gaps;
     cout << Date( ) << ": finding gaps" << endl;
     for ( int d = 0; d < D.E( ); d++ )
     {    if ( D.O(d)[0] != -1 ) continue;
          if ( ALL && dinv[d] < d ) continue;
          int v = to_left[d], w = to_right[d];
          if ( !D.From(v).solo( ) || !D.To(v).solo( ) ) continue;
          if ( !D.From(w).solo( ) || !D.To(w).solo( ) ) continue;
          int d1 = D.ITo(v,0), d2 = D.IFrom(w,0);
          if ( D.O(d1)[0] < 0 || D.O(d2)[0] < 0 ) continue;
          gaps.push_back(d);    }
     cout << Date( ) << ": total number of pair gaps = "
          << ToStringAddCommas( gaps.size( ) ) << endl;

     // Heuristics.

     const int MAX_CLOSURES = 4;
     const int USELEN = 400;
     const int MAX_EDGES = 50;
     const int MM = 8;
     const int MAX_BRIDGES = 5000;
     const int HSTAR = 5;

     // Track basic stats.

     int ngaps = Max( rs.size( ), ss.size( ) );
     if (ALL) ngaps = gaps.size( );
     int nclosed = 0, nclosures = 0, nclosures_true = 0;
     int semiclosed = 0, ugly = 0;

     // Go through the specified gaps.

     int nthreads = ( ( DIRECT || VERBOSE ) ? 1 : omp_get_max_threads( ) );
     vec<vec<basevector>> bpaths(ngaps);
     vec<int> thegap(ngaps);
     double mclock = WallClockTime( );
     int ndots = 0, done = 0;
     #pragma omp parallel for num_threads(nthreads) \
          firstprivate( bases, quals, xpaths, xpaths_index ) schedule(dynamic, 1)
     for ( int ng = 0; ng < ngaps; ng++ )
     {
          // Find the gap.

          int R = -1, S = -1;
          if (ALL) S = gaps[ng];
          else if ( ss.nonempty( ) ) S = ss[ng];
          else
          {    R = rs[ng];
               int delta = ( ng == 0 ? rs[ng] : rs[ng] - rs[ng-1] );
               for ( int j = 0; j < delta; j++ ) 
                    S = gaps[ randomx( ) % gaps.size( ) ];    }

          ostringstream goutx;
          ostream& gout = ( DIRECT ? cout : goutx );
          gout << "\n========================================================"
               << "============================\n\n";
          if ( R >= 0 ) { DPRINT2_TO( gout, R, S ); }
          else DPRINT_TO( gout, S );
          thegap[ng] = S;

          // Get edge info.

          double clock1 = WallClockTime( );
          int v = to_left[S], w = to_right[S];
          int d1 = D.ITo(v,0), d2 = D.IFrom(w,0);
          DPRINT2_TO( gout, d1, d2 );
          vecbasevector edgesx;
          edgesx.push_back( hb.Cat( D.O(d1) ) );
          edgesx.push_back( hb.Cat( D.O(d2) ) );
          vecbasevector edgesy(edgesx);
          if ( edgesy[0].isize( ) > USELEN )
               edgesy[0].SetToSubOf( 
                    edgesx[0], edgesx[0].isize( ) - USELEN, USELEN );
          if ( edgesy[1].isize( ) > USELEN ) edgesy[1].resize(USELEN);
          vec<basevector> edges(2);
          for ( int i = 0; i < 2; i++ )
               edges[i] = edgesx[i];
          if (VERBOSE)
          {    gout << "\nBOUNDING EDGELETS" << endl;
               edgesy[0].Print( gout, "left" );
               edgesy[1].Print( gout, "right" );
               gout << endl;    }
     
          // Find patches.

          gout << Date( ) << ": start patching" << endl;
          vec<basevector> closures;
          const Bool ALT = True;
          vec<int> trim;
          vec< pair<int64_t,Bool> > idsfw2;
          const int MAX_DIST = 800;
          const int MAX_PROX = 400;
          Bool super_verbose = False;
     
          // Build read set.

          for ( int spass = 1; spass <= 2; spass++ )
          {    int d = ( spass == 1 ? d1 : dinv[d2] );
               if (super_verbose) PRINT_TO(gout,d);
               if ( d < 0 ) break;
               for ( int pass = 1; pass <= 2; pass++ )
               {    const int f = ( pass == 1 ? d : dinv[d] );
                    for ( int i = 0; i < (int) xpaths_index[f].size( ); i++ )
                    {    int64_t id = xpaths_index[f][i];
                         if (super_verbose) PRINT_TO(gout,id);
                         const ReadPath& p = xpaths[id];
                         int N = edges[spass-1].size( );
                         int offset = p.getOffset( );
                         for ( int j = 0; j < (int) p.size( ); j++ )
                         {    if ( p[j] != f )
                              {     offset -= dkmers[ p[j] ];
                                    continue;    }
                              int start;
                              if ( pass == 1 ) start = offset;
                              else start = dkmers[f] + K - 1 
                                   - offset - bases[id].isize( );
                              int stop = start + bases[id].isize( );
                              if (super_verbose) 
                                   PRINT4_TO( gout, id, pass, stop, N );
                              if ( ALT && N - stop > MAX_PROX ) continue;
                              idsfw2.push( 
                                   id, spass == 1 ^ pass == 2 );    }    }    }
               for ( int i = 0; i < (int) xpaths_index[d].size( ); i++ )
               {    int64_t id1 = xpaths_index[d][i];
                    int64_t id2 = ( id1 % 2 == 0 ? id1+1 : id1-1 );
                    const ReadPath& p1 = xpaths[id1];
                    for ( int j = 0; j < (int) p1.size( ); j++ )
                    {    if ( p1[j] == d )
                         {    int pos1 = p1.getOffset( );
                              for ( int l = 0; l < j; l++ )
                                   pos1 -= dkmers[ p1[l] ];
                              if ( dkmers[d] + K - 1 - pos1 > MAX_DIST ) continue;
                              idsfw2.push( id2, spass == 2 );    }    }    }    }
          UniqueSort(idsfw2);
          gout << Date( ) << ": " << TimeSince(clock1) << " used on segment 1" 
               << endl;

          // Create solidtigs.

          double rclock = WallClockTime( );
          vecbasevector basesx;
          vecqualvector qualsx( idsfw2.size( ) );
          for ( int i = 0; i < idsfw2.isize( ); i++ )
          {    int64_t id = idsfw2[i].first;
               basesx.push_back( bases[id] );
               quals[id].unpack( &qualsx[i] );
               if ( !idsfw2[i].second ) 
               {    basesx.back( ).ReverseComplement( );
                    qualsx.back( ).ReverseMe( );    }    }
          if ( VERBOSE && ngaps == 1 ) basesx.WriteAll( "snorgle.fastb" );
          gout << Date( ) << ": " << TimeSince(rclock) << " used fetching reads" 
               << endl;
          gout << Date( ) << ": there are " << ToStringAddCommas( basesx.size( ) )
               << " reads" << endl;

          // Lower quality scores in homopolymers.

          double clock2 = WallClockTime( );
          vecbasevector basesy(basesx);
          vecqualvector qualsy(qualsx);
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

          // Directly generate the graph.  This should only be forward, because the
          // reads are oriented.  Horrendously inefficient calculation.

          HyperBasevector hbl;
          vecbasevector moo(basesy);
          moo.Append(edgesy);
          BasesToGraph( moo, K, hbl );
          if (VERBOSE) gout << endl;
          gout << Date( ) << ": local graph has " << hbl.E( ) << " edges" << endl;
     
          // Find edge support.  Insane duplication of calculation from BasesToGraph.

          vec<vec<int>> support( hbl.E( ) );
          {    const int K = 48;
               ForceAssertEq( K, hb.K( ) );
               vec< triple<kmer<K>,int,int> > kmers_plus;
               vecbasevector allb = basesy;
               for ( int e = 0; e < hbl.E( ); e++ )
                    allb.push_back( hbl.O(e) );
               MAKE_KMER_LOOKUP( allb, kmers_plus );
               for ( int i = 0; i < kmers_plus.isize( ); i++ )
               {    int j;
                    for ( j = i + 1; j < kmers_plus.isize( ); j++ )
                         if ( kmers_plus[j].first != kmers_plus[i].first ) break;
                    for ( int k1 = i; k1 < j; k1++ )
                    for ( int k2 = i; k2 < j; k2++ )
                    {    int id1 = kmers_plus[k1].second;
                         if ( id1 >= (int) basesy.size( ) ) continue;
                         int id2 = kmers_plus[k2].second - (int) basesy.size( );
                         if ( id2 < 0 ) continue;
                         support[id2].push_back(id1);    }
                    i = j - 1;    }
               for ( int e = 0; e < hbl.E( ); e++ )
                    UniqueSort( support[e] );    }

          // Delete weak edges.

          gout << Date( ) << ": deleting weak edges" << endl;
          vec<int> dels;
          for ( int v = 0; v < hbl.N( ); v++ )
          {    if ( hbl.From(v).size( ) <= 1 ) continue;
               vec<int> s( hbl.From(v).size( ) );
               for ( int j = 0; j < hbl.From(v).isize( ); j++ )
                    s[j] = support[ hbl.IFrom(v,j) ].size( );
               vec<int> ids( hbl.From(v).size( ), vec<int>::IDENTITY );
               ReverseSortSync( s, ids );
               for ( int j = 1; j < s.isize( ); j++ )
               {    int e = hbl.IFrom( v, ids[j] );
                    if ( s[0] >= 3 && s[j] == 0 ) dels.push_back(e);
                    if ( s[0] >= 6 && s[j] == 1 ) dels.push_back(e);
                    if ( s[0] >= 8 && s[j] == 2 ) dels.push_back(e);    }    }
          for ( int v = 0; v < hbl.N( ); v++ )
          {    if ( hbl.To(v).size( ) <= 1 ) continue;
               vec<int> s( hbl.To(v).size( ) );
               for ( int j = 0; j < hbl.To(v).isize( ); j++ )
                    s[j] = support[ hbl.ITo(v,j) ].size( );
               vec<int> ids( hbl.To(v).size( ), vec<int>::IDENTITY );
               ReverseSortSync( s, ids );
               for ( int j = 1; j < s.isize( ); j++ )
               {    int e = hbl.ITo( v, ids[j] );
                    if ( s[0] >= 3 && s[j] == 0 ) dels.push_back(e);
                    if ( s[0] >= 6 && s[j] == 1 ) dels.push_back(e);
                    if ( s[0] >= 8 && s[j] == 2 ) dels.push_back(e);    }    }
          for ( int e = 0; e < hbl.E( ); e++ )
               if ( support[e].size( ) <= 2 ) dels.push_back(e);
          hbl.DeleteEdges(dels);
          hbl.RemoveUnneededVertices( );
          hbl.RemoveDeadEdgeObjects( );
          hbl.RemoveEdgelessVertices( );
          dels.clear( );
          for ( int e = 0; e < hbl.E( ); e++ )
               if ( hbl.Kmers(e) < 20 ) dels.push_back(e);
          hbl.DeleteEdges(dels);
          hbl.RemoveUnneededVertices( );
          hbl.RemoveDeadEdgeObjects( );
          hbl.RemoveEdgelessVertices( );
          // ---->

     // Again find edge support.  Insane duplication.

     gout << Date( ) << ": finding edge support" << endl;
     support.clear( );
     support.resize( hbl.E( ) );
     {    const int K = 48;
          ForceAssertEq( K, hb.K( ) );
          vec< triple<kmer<K>,int,int> > kmers_plus;
          vecbasevector allb = basesy;
          for ( int e = 0; e < hbl.E( ); e++ )
               allb.push_back( hbl.O(e) );
          MAKE_KMER_LOOKUP( allb, kmers_plus );
          for ( int i = 0; i < kmers_plus.isize( ); i++ )
          {    int j;
               for ( j = i + 1; j < kmers_plus.isize( ); j++ )
                    if ( kmers_plus[j].first != kmers_plus[i].first ) break;
               for ( int k1 = i; k1 < j; k1++ )
               for ( int k2 = i; k2 < j; k2++ )
               {    int id1 = kmers_plus[k1].second;
                    if ( id1 >= (int) basesy.size( ) ) continue;
                    int id2 = kmers_plus[k2].second - (int) basesy.size( );
                    if ( id2 < 0 ) continue;
                    support[id2].push_back(id1);    }
               i = j - 1;    }
          for ( int e = 0; e < hbl.E( ); e++ )
               UniqueSort( support[e] );    }

     // Compute read paths.  This is very much suboptimal.
     // For now, not setting offset.
     // For now, just a list of edges, ignoring dups and gaps.

     gout << Date( ) << ": computing read paths" << endl;
     ReadPathVec pathsl( basesy.size( ) );
     {    const int K = 48;
          ForceAssertEq( K, hb.K( ) );
          vec< triple<kmer<K>,int,int> > kmers_plus;
          vecbasevector allb = basesy;
          for ( int e = 0; e < hbl.E( ); e++ )
               allb.push_back( hbl.O(e) );
          MAKE_KMER_LOOKUP( allb, kmers_plus );
          vec<triple<int,int,int>> places;
          for ( int i = 0; i < kmers_plus.isize( ); i++ )
          {    int j;
               for ( j = i + 1; j < kmers_plus.isize( ); j++ )
                    if ( kmers_plus[j].first != kmers_plus[i].first ) break;
               for ( int k1 = i; k1 < j; k1++ )
               for ( int k2 = i; k2 < j; k2++ )
               {    int id = kmers_plus[k1].second;
                    if ( id >= (int) basesy.size( ) ) continue;
                    int e = kmers_plus[k2].second - (int) basesy.size( );
                    if ( e < 0 ) continue;
                    int offset = kmers_plus[k1].third - kmers_plus[k2].third;
                    // edge e is placed on read id at the given offset
                    places.push( id, offset, e );    }
               i = j - 1;    }
          Sort(places);
          for ( int i = 0; i < places.isize( ); i++ )
          {    int j;
               int id = places[i].first;
               for ( j = i + 1; j < places.isize( ); j++ )
                    if ( places[j].first != id ) break;
               pathsl[id].push_back( places[i].third );
               for ( int k = i+1; k < j; k++ )
               {    int e = places[k].third;
                    if ( e != pathsl[id].back( ) ) pathsl[id].push_back(e);    }
               i = j - 1;    }    }

     // Print graph.

     gout << Date( ) << ": now local graph has " << hbl.E( ) << " edges" << endl;
     /*
     if (VERBOSE)
     {    Ofstream( out, "whatever.dot" );
          hbl.PrintSummaryDOT0w( out, False, False, True );
          for ( int e = 0; e < hbl.E( ); e++ )
          {    hbl.O(e).Print( gout, 
                    ToString(e) + "[" + ToString( support[e].size( ) ) + "]" );    }
          gout << endl;    }
     */
     if ( hbl.E( ) > MAX_EDGES )
     {    gout << Date( ) << ": too many edges, giving up" << endl;
          gout << "\nCONCLUSION: there are 0 closures for ";
          if ( R > 0 ) gout << "R = " << R << endl;
          else gout << "S = " << S << endl;
          if ( VERBOSITY == 0 ) 
          {   
               #pragma omp critical
               {    MakeDots( done, ndots, ngaps );    }    }
          if ( !DIRECT && VERBOSITY > 0 )
          {
               #pragma omp critical
               {    cout << goutx.str( );    }    }
          continue;    }

     // Show read pair linking between edges.

     vec< pair<int,int> > links;
     for ( int e1 = 0; e1 < hbl.E( ); e1++ )
     for ( int e2 = 0; e2 < hbl.E( ); e2++ )
     for ( int j1 = 0; j1 < support[e1].isize( ); j1++ )
     for ( int j2 = 0; j2 < support[e2].isize( ); j2++ )
     {    if ( e1 == e2 ) continue;
          int id1 = support[e1][j1], id2 = support[e2][j2];
          int64_t ID1 = idsfw2[id1].first, ID2 = idsfw2[id2].first;
          if ( ID1/2 != ID2/2 ) continue;
          if ( !idsfw2[id1].second || idsfw2[id2].second ) continue;
          const ReadPath &p1 = pathsl[id1], &p2 = pathsl[id2];
          vec<int> x1, x2;
          for ( int j = 0; j < (int) p1.size( ); j++ )
               x1.push_back( p1[j] );
          for ( int j = 0; j < (int) p2.size( ); j++ )
               x2.push_back( p2[j] );
          if ( Position( x1, e1 ) >= 0 && Position( x1, e2 ) >= 0
               && Position( x1, e2 ) < Position( x1, e1 ) )
          {    continue;    }
          if ( Position( x2, e1 ) >= 0 && Position( x2, e2 ) >= 0
               && Position( x2, e2 ) < Position( x2, e1 ) )
          {    continue;    }
          links.push( e1, e2 );    }
     UniqueSort(links);
     if (VERBOSE)
     {    gout << "READ PAIR LINKS" << endl << endl;
          for ( int i = 0; i < links.isize( ); i++ )
          {    int j = links.NextDiff(i);
               gout << links[i].first << " --> " << links[i].second
                    << " [" << j-i << "]" << endl;
               i = j - 1;    }
          gout << endl;    }

     // Find connected components, and which are linked to which.

     gout << Date( ) << ": find connected components" << endl;
     vec<vec<int>> comp;
     hbl.ComponentsEFast(comp);
     vec<vec<int>> linksto( comp.size( ) );
     vec<int> tocomp( hbl.E( ) );
     for ( int j = 0; j < comp.isize( ); j++ )
     for ( int l = 0; l < comp[j].isize( ); l++ )
          tocomp[ comp[j][l] ] = j;
     for ( int i = 0; i < links.isize( ); i++ )
     {    int e1 = links[i].first, e2 = links[i].second;
          int c1 = tocomp[e1], c2 = tocomp[e2];
          if ( c1 != c2 ) linksto[c1].push_back(c2);    }
     for ( int i = 0; i < comp.isize( ); i++ )
          UniqueSort( linksto[i] );
     if (VERBOSE)
     {    gout << "COMPONENTS AND LINKING" << endl;
          for ( int i = 0; i < comp.isize( ); i++ )
          {    Sort( comp[i] );
               gout << "[" << i+1 << "] = " << printSeq( comp[i] ) << " -->";
               for ( int j = 0; j < linksto[i].isize( ); j++ )
                    gout << " [" << linksto[i][j] + 1 << "]";
               gout << endl;    }
          gout << endl;    }
     gout << Date( ) << ": " << TimeSince(clock2) << " used on segment 2" << endl;

     // Compute eightmer frequencies in the assembly.

     vec<kmer<MM>> hlocs;
     for ( int e = 0; e < hbl.E( ); e++ )
     {    const basevector& b = hbl.O(e);
          kmer<MM> x;
          for ( int j = 0; j <= b.isize( ) - K + MM; j++ )
          {    x.SetToSubOf( b, j );
               hlocs.push_back(x);    }    }
     Sort(hlocs);
     
     // Look directly for bridges.

     gout << Date( ) << ": finding bridges" << endl;
     double bclock = WallClockTime( );
     {
     if (VERBOSE) gout << "ANCHORS" << endl;
     vecbasevector all(basesy);
     for ( int e = 0; e < hbl.E( ); e++ )
          all.push_back( hbl.O(e) );
     vec< triple<kmer<MM>,int,int> > kmers_plus;
     MAKE_KMER_LOOKUP( all, kmers_plus );
     vec< triple<int,int,int> > matches;
     const int MAX_MATCHES = 1000000;
     for ( int i = 0; i < kmers_plus.isize( ); i++ )
     {    int j;
          for ( j = i + 1; j < kmers_plus.isize( ); j++ )
               if ( kmers_plus[j].first != kmers_plus[i].first ) break;
          vec<int> count( hbl.E( ), 0 );
          for ( int m = 0; m < hbl.E( ); m++ )
          for ( int k = i; k < j; k++ )
               if ( kmers_plus[k].second == (int) basesy.size( ) + m ) count[m]++;
          if ( Sum(count) == 0 )
          {    i = j - 1;
               continue;    }
           
          // Check kmer for acceptable.  First half can't be a homopolymer, nor
          // can second half, and kmer must contain at least two different bases
          // and can't be a dinuke.
          // And count can't be 6:2.

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
                    const int MIN_BASES = 2;
                    if ( Sum(present) >= MIN_BASES ) good = True;    }    }
          if (good)
          {    diff = False;
               for ( int m = 2; m < MM; m += 2 )
                    if ( bb[m] != bb[m-2] || bb[m+1] != bb[m-1] ) diff = True;
               if ( !diff ) good = False;    }

          if (good)
          {    vec<Bool> count(4,0);
               for ( int m = 0; m < MM; m++ ) count[ bb[m] ]++;
               ReverseSort(count);
               if ( count[0] == 6 && count[1] == 2 ) good = False;    }

          const int MAX_MM_FREQ = 1;
          if (good)
          {    if ( UpperBound( hlocs, kmers_plus[i].first )
                    - LowerBound( hlocs, kmers_plus[i].first ) > MAX_MM_FREQ )
               {    good = False;    }    }
     
          if ( !good )
          {    i = j - 1;   
               continue;    }

          for ( int k2 = i; k2 < j; k2++ )
          {    int m = kmers_plus[k2].second - (int) basesy.size( );
               if ( m < 0 ) continue;
               for ( int k1 = i; k1 < j; k1++ )
               {    if ( kmers_plus[k1].second >= (int) basesy.size( ) ) continue;
                    int cpos = kmers_plus[k2].third;
                    int rid = kmers_plus[k1].second, rpos = kmers_plus[k1].third;
     
                    // Validate.  Require an overlap of size <= K, having at least
                    // K/2 bases in agreement, and only test the three cases where
                    // the MM-mer is on the left end, in the middle, or on the right
                    // end.
     
                    Bool valid = False;
                    const basevector &r = basesy[rid], &c = hbl.O(m);
                    int match = MM;
                    for ( int l = 0; l < K - MM; l++ )
                    {    if ( rpos+MM+l >= r.isize( ) || cpos+MM+l >= c.isize( ) ) 
                              break;
                         if ( r[rpos+MM+l] == c[cpos+MM+l] ) match++;    }
                    if ( match >= K/2 ) valid = True;
                    if ( !valid )
                    {    match = MM;
                         for ( int l = 1; l <= K - MM; l++ )
                         {    if ( rpos-l <= 0 || cpos-l <= 0 ) break;
                              if ( r[rpos-l] == c[cpos-l] ) match++;    }
                         if ( match >= K/2 ) valid = True;
                         if ( !valid )
                         {    match = MM;
                              for ( int l = 0; l < (K-MM)/2; l++ )
                              {    if ( rpos+MM+l >= r.isize( ) 
                                        || cpos+MM+l >= c.isize( ) ) 
                                   {    break;    }
                                   if ( r[rpos+MM+l] == c[cpos+MM+l] ) match++;    }
                              for ( int l = 1; l <= (K-MM)/2; l++ )
                              {    if ( rpos-l <= 0 || cpos-l <= 0 ) break;
                                   if ( r[rpos-l] == c[cpos-l] ) match++;    }
                              if ( match >= K/2 ) valid = True;    }    }
                    if ( !valid ) continue;
     
                    // Record overlap.
     
                    int offset = cpos - rpos;
                    matches.push( rid, m, offset );
                    if ( matches.isize( ) > MAX_MATCHES ) break;    }
               if ( matches.isize( ) > MAX_MATCHES ) break;    }
          i = j - 1;    }
     gout << Date( ) << ": " << TimeSince(bclock) << " used finding matches" << endl;
     if ( matches.isize( ) > MAX_MATCHES )
     {    gout << Date( ) << ": too many matches, giving up" << endl;
          gout << "\nCONCLUSION: there are 0 closures for ";
          if ( R > 0 ) gout << "R = " << R << endl;
          else gout << "S = " << S << endl;
          if ( VERBOSITY == 0 ) 
          {   
               #pragma omp critical
               {    MakeDots( done, ndots, ngaps );    }    }
          if ( !DIRECT && VERBOSITY > 0 )
          {
               #pragma omp critical
               {    cout << goutx.str( );    }    }
          continue;    }
     gout << Date( ) << ": initially there are " 
          << ToStringAddCommas( matches.size( ) ) << " matches" << endl;

     // Screen matches.

     gout << Date( ) << ": screening matches" << endl;
     double sclock = WallClockTime( );
     UniqueSort(matches);
     vec< quad<int,int,int,double> > matches1, matches2;
     for ( int i = 0; i < matches.isize( ); i++ )
     {    int rid = matches[i].first, m = matches[i].second; 
          int offx = matches[i].third;

          // Compute bits (v1).

          const int WID = 20;
          const int MAX_OVERLAP = 1000;
          static PrecomputedBinomialSums<MAX_OVERLAP> gBS(WID,0.75);
          int pos1 = ( offx >= 0 ? offx : 0 ), pos2 = ( offx <= 0 ? -offx : 0 );
          int overlap = IntervalOverlap( 
               0, hbl.Bases(m), offx, offx + basesy[rid].isize( ) );
          double minp = 0;
          vec<int> err_sum( overlap + 1, 0 );
          for ( int v = 1; v <= overlap; v++ )
          {    err_sum[v] = err_sum[v-1];
               if ( hbl.O(m)[pos1+v-1] != basesy[rid][pos2+v-1] ) err_sum[v]++;    }
          for ( int start = 0; start < overlap; start++ )
          {    for ( int n = WID; n <= overlap - start; n++ )
               {    int k = err_sum[start+n] - err_sum[start];
                    minp = Min( minp, gBS[n][k] );    }    }
          double bits = -minp * 10.0 / 6.0;
          const double MIN_BITS = 20.0;
          if ( bits < MIN_BITS ) continue;
          matches1.push( rid, m, offx, bits );    }
     gout << Date( ) << ": " << TimeSince(sclock) << " used computing bits" << endl;
     double tclock = WallClockTime( );
     vec<Bool> to_delete1( matches1.size( ), False );
     for ( int i = 0; i < matches1.isize( ); i++ )
     {    int j;
          vec<int> ms;
          for ( j = i; j < matches1.isize( ); j++ )
          {    if ( matches1[j].first != matches1[i].first ) break;
               ms.push_back( matches1[j].second );    }
          UniqueSort(ms);
          if ( ms.solo( ) )
          {    for ( int k = i; k < j; k++ )
                    to_delete1[k] = True;    }
          i = j - 1;    }
     EraseIf( matches1, to_delete1 );
     for ( int i = 0; i < matches1.isize( ); i++ )
     {    int rid = matches1[i].first, m = matches1[i].second; 
          int offx = matches1[i].third;
          double bits = matches1[i].fourth;
          int j;
          double maxbits = bits;
          for ( j = i + 1; j < matches1.isize( ); j++ )
          {    if ( matches1[j].first != matches1[i].first ) break;
               if ( matches1[j].second != matches1[i].second ) break;
               maxbits = Max( bits, matches1[i].fourth );    }
          if ( bits >= maxbits - 10.0 )
          {
               // Show the anchor.

               if (VERBOSE)
               {    int pos1 = ( offx >= 0 ? offx : 0 ); 
                    int pos2 = ( offx <= 0 ? -offx : 0 );
                    int overlap = IntervalOverlap( 
                         0, hbl.Bases(m), offx, offx + basesy[rid].isize( ) );
                    gout << endl;
                    PRINT4_TO( gout, rid, m, offx, bits );
                    avector<int> gaps(1), lengths(1);
                    gaps(0) = 0;
                    lengths(0) = overlap;
                    align a( pos2, pos1, gaps, lengths );
                    PrintVisualAlignmentClean( 
                         True, gout, basesy[rid], hbl.O(m), a, qualsy[rid] );    }
     
               // Save.
     
               matches2.push( rid, m, offx, bits );    }

          i = j - 1;    }
     gout << Date( ) << ": " << TimeSince(tclock) << " used screening matches"
          << endl;
     gout << Date( ) << ": now have " << ToStringAddCommas( matches2.size( ) )
          << " matches" << endl;
     
     // Find the bridges.

     double bclock2 = WallClockTime( );
     gout << Date( ) << ": finding bridges" << endl;
     if (VERBOSE) gout << "\nBRIDGES" << endl << endl;
     vec< pair<int,int> > loc;
     vec<int> lleft, lright;
     hbl.ToLeft(lleft), hbl.ToRight(lright);
     vec< quad< pair<int,int>, int, pair<double,double>, int > > mogo;
     for ( int i = 0; i < matches2.isize( ); i++ )
     {    int j;
          for ( j = i + 1; j < matches2.isize( ); j++ )
               if ( matches2[j].first != matches2[i].first ) break;
          for ( int k1 = i; k1 < j; k1++ )
          for ( int k2 = i; k2 < j; k2++ )
          {    int e1 = matches2[k1].second, e2 = matches2[k2].second;
               if ( e1 == e2 ) continue;
               if ( hbl.From( lright[e1] ).nonempty( ) ) continue;
               if ( hbl.To( lleft[e2] ).nonempty( ) ) continue;
               int off = matches2[k1].third - matches2[k2].third;
               double bits1 = matches2[k1].fourth, bits2 = matches2[k2].fourth;
               mogo.push( make_pair(e1, e2), off, make_pair(bits1, bits2),
                    matches2[i].first );    }
          i = j - 1;    }
     gout << Date( ) << ": " << TimeSince(bclock2) << " used bridging" << endl;

     // Filter.

     double fclock = WallClockTime( );
     Sort(mogo);
     vec<Bool> to_delete( mogo.size( ), False );
     for ( int r = 0; r < mogo.isize( ); r++ )
     {    int s;
          for ( s = r + 1; s < mogo.isize( ); s++ )
               if ( mogo[s].first != mogo[r].first ) break;
          double best = -1;
          for ( int i = r; i < s; i++ )
          {    double b1 = mogo[i].third.first, b2 = mogo[i].third.second;
               best = Max( best, Min( b1, b2 ) );    }
          for ( int i = r; i < s; i++ )
          {    int j;
               for ( j = i + 1; j < s; j++ )
                    if ( mogo[j].second != mogo[i].second ) break;
               if ( j - i == 1 )
               {    int m = Min( mogo[i].third.first, mogo[i].third.second );
                    if ( m <= 35 && best - m >= 20 ) to_delete[i] = True;    }
               i = j - 1;    }
          r = s - 1;    }
     EraseIf( mogo, to_delete );

     // Filter the bridges more.

     vec<triple<int,int,int>> bridges;
     for ( int i = 0; i < mogo.isize( ); i++ )
     {    int e1 = mogo[i].first.first, e2 = mogo[i].first.second;
          int c1 = tocomp[e1], c2 = tocomp[e2];
          int rid = mogo[i].fourth;

          // Note below that we're counting links between e1 and e2, but not
          // checking that the offset is the same.  I don't think that makes
          // sense.

          int count = 1;
          for ( int j = i - 1; j >= 0; j-- )
          {    if ( mogo[j].first != mogo[i].first ) break;
               count++;    }
          for ( int j = i + 1; j < mogo.isize( ); j++ )
          {    if ( mogo[j].first != mogo[i].first ) break;
               count++;    }
          if ( count == 1 && !BinMember( linksto[c1], c2 ) ) continue;

          int off = mogo[i].second;
          int gap = off - hbl.Bases(e1);
          // MIGHT NEED TO INCREASE 12!  INCREASED TO 40.
          if ( gap < -(K+40) ) continue;
          double bits1 = mogo[i].third.first, bits2 = mogo[i].third.second;
          if (VERBOSE)
          {    gout << e1 << " --(" << gap << ")--> " << e2 << " [" << bits1 << "," 
                    << bits2 << ";rid=" << rid << "]" << endl;    }
          bridges.push( e1, e2, off );    }
     if (VERBOSE) gout << endl;
     gout << Date( ) << ": " << TimeSince(bclock2) << " used filtering" << endl;

     // Try to walk across each bridge.

     double clock3 = WallClockTime( );
     gout << Date( ) << ": walking " << bridges.size( ) << " bridges" << endl;
     if ( bridges.isize( ) > MAX_BRIDGES )
     {    gout << Date( ) << ": too many bridges, giving up" << endl;
          gout << "\nCONCLUSION: there are 0 closures for ";
          if ( R > 0 ) gout << "R = " << R << endl;
          else gout << "S = " << S << endl;
          if ( VERBOSITY == 0 ) 
          {   
               #pragma omp critical
               {    MakeDots( done, ndots, ngaps );    }    }
          if ( !DIRECT && VERBOSITY > 0 )
          {
               #pragma omp critical
               {    cout << goutx.str( );    }    }
          continue;    }
     vec<basevector> bseq;
     {
     ostringstream out;
     UniqueSort(bridges);
     for ( int br = 0; br < bridges.isize( ); br++ )
     {    vec<readstack> stacks(2);
          vecbasevector EE;
          int e1 = bridges[br].first, e2 = bridges[br].second;
          EE.push_back( hbl.O(e1) );
          EE.push_back( hbl.O(e2) );
          for ( int si = 0; si < 2; si++ )
          {    readstack& stack = stacks[si];
               const basevector b = EE[si];
               vecbasevector wow = basesy;
               wow.push_back(b);
               if ( si == 1 )
               {    for ( int i = 0; i < (int) wow.size( ); i++ )
                         wow[i].ReverseComplement( );    }
               const int K = 48;
               vec< triple<kmer<K>,int,int> > kmers_plus;
               MAKE_KMER_LOOKUP( wow, kmers_plus );

               vec< triple<int,int,int> > locs;
               for ( int i = 0; i < kmers_plus.isize( ); i++ )
               {    int j;
                    for ( j = i + 1; j < kmers_plus.isize( ); j++ )
                         if ( kmers_plus[j].first != kmers_plus[i].first ) break;
                    for ( int k2 = i; k2 < j; k2++ )
                    {    int id2 = kmers_plus[k2].second - (int) basesy.size( );
                         if ( id2 != 0 ) continue;
                         for ( int k1 = i; k1 < j; k1++ )
                         {    int id1 = kmers_plus[k1].second;
                              if ( id1 >= (int) basesy.size( ) ) continue;
                              int pos1 = kmers_plus[k1].third;
                              int pos2 = kmers_plus[k2].third;
                              if ( pos1 > 0 && pos2 > 0 
                                   && b[pos2-1] == basesy[id1][pos1-1] )
                              {    continue;    }
                              int n = 1;
                              while(1)
                              {    pos1++;
                                   pos2++;
                                   if ( pos2 == b.isize( ) ) break;
                                   if ( pos1 == basesy[id1].isize( ) ) break;
                                   if ( b[pos2] != basesy[id1][pos1] ) break;
                                   n++;    }
                              int off = pos2 - pos1;
                              locs.push( id1, off, n );    }    }
                    i = j - 1;    }
               Sort(locs);
               vec<triple<int,int,int>> locs2;
               for ( int i = 0; i < locs.isize( ); i++ )
               {    int j;
                    int n = locs[i].third;
                    for ( j = i + 1; j < locs.isize( ); j++ )
                    {    if ( locs[j].first != locs[i].first ) break;
                         if ( locs[j].second != locs[i].second ) break;
                         n += locs[j].third;    }
                    locs2.push( locs[i].first, n, locs[i].second );
                    i = j - 1;    }
               vec<pair<int,int>> locs3;
               for ( int i = 0; i < locs2.isize( ); i++ )
               {    int j;
                    for ( j = i + 1; j < locs2.isize( ); j++ )
                         if ( locs2[j].first != locs2[i].first ) break;
                    locs3.push( locs2[i].first, locs2[i].third );
                    i = j - 1;    }
               int M = 0;
               for ( int i = 0; i < locs3.isize( ); i++ )
               {    int id = locs3[i].first, pos = locs3[i].second;
                    M = Max( M, pos + basesy[id].isize( ) );    }
               stack.Initialize( locs3.size( ), M );
               for ( int i = 0; i < locs3.isize( ); i++ )
               {    int64_t id = locs3[i].first;
                    basevector b = basesy[id];
                    qualvector q = qualsy[id];
                    if ( si == 1 )
                    {    b.ReverseComplement( );
                         q.ReverseMe( );    }
                    Bool fw = idsfw2[id].second;
                    int pos = locs3[i].second;
                    stack.SetId( i, id );
                    stack.SetPid( i, id/2 );
                    stack.SetLen( i, b.size( ) );
                    stack.SetRc2( i, !fw );
                    stack.SetOffset( i, pos );
                    for ( int j = 0; j < b.isize( ); j++ )
                    {    int p = pos + j;
                         if ( p >= 0 && p < stack.Cols( ) )
                         {    stack.SetBase( i, p, b[j] );
                              stack.SetQual( i, p, q[j] );    }    }    }    }
          stacks[1].Reverse( );
          int off = bridges[br].third;
          {    int gap = off - EE[0].isize( );
               readstack s( stacks[0] );
               int offset = off - ( stacks[1].Cols( ) - EE[1].isize( ) );
               s.Merge( stacks[1], offset );

               // Not sure this does anything, but it seems like it should:

               s.SortByPid( -1, -1, -1 );
               s.Unique( );

               basevector f;
               qualvector fragq;
               s.Consensus1( f, fragq );
          
               // Analyze stack for minimal consistency.

               const int flank = 2;
               Bool good = True;
               for ( int i = flank; i < s.Cols( ) - flank; i++ )
               {    Bool ok = False;
                    for ( int j = 0; j < s.Rows( ); j++ )
                    {    Bool mismatch = False;
                         for ( int l = i - flank; l <= i + flank; l++ )
                         {    if ( s.Base( j, l ) != f[l] )
                              {    mismatch = True;
                                   break;    }    }
                         if ( !mismatch )
                         {    ok = True;
                              break;    }    }
                    if ( !ok )
                    {    good = False;
                         break;     }    }
               if ( !good )
               {    int ltrim = -offset;
                    if ( ltrim < 0 ) ltrim = 0;
                    int rtrim = stacks[0].Cols( ) - offset - stacks[1].Cols( );
                    if ( rtrim < 0 ) rtrim = 0;
                    if ( ltrim + rtrim < s.Cols( ) )
                         s.Trim( ltrim, s.Cols( ) - rtrim );
                    s.Consensus1( f, fragq );

                    if (VERBOSE)
                    {    good = True;
                         for ( int i = flank; i < s.Cols( ) - flank; i++ )
                         {    Bool ok = False;
                              for ( int j = 0; j < s.Rows( ); j++ )
                              {    Bool mismatch = False;
                                   for ( int l = i - flank; l <= i + flank; l++ )
                                   {    if ( s.Base( j, l ) != f[l] )
                                        {    mismatch = True;
                                             break;    }    }
                                   if ( !mismatch )
                                   {    ok = True;
                                        break;    }    }
                              if ( !ok )
                              {    good = False;
                                   break;     }    }    }    }

               // Look for consistent closures.

               if (VERBOSE)
               {
               double cclock = WallClockTime( );
               int reg = 2*flank + 1;
               vec<String> walksx, walks;
               vec<int> support;
               String x(5);
               for ( int i = 0; i < s.Rows( ); i++ )
               {    Bool defined = True;
                    for ( int j = 0; j < reg; j++ )
                         if ( !s.Def(i,j) ) defined = False;
                    if ( !defined ) continue;
                    for ( int j = 0; j < reg; j++ )
                         x[j] = s.Base(i,j);
                    walksx.push_back(x);    }
               Sort(walksx);
               int max_support = 0;
               for ( int i = 0; i < walksx.isize( ); i++ )
               {    int j = walksx.NextDiff(i);
                    max_support = Max( max_support, j - i );
                    i = j - 1;    }
               for ( int i = 0; i < walksx.isize( ); i++ )
               {    int j = walksx.NextDiff(i);
                    if ( j - i >= max_support/2 )
                    {    walks.push_back( walksx[i] );
                         support.push_back( j - i );    }
                    i = j - 1;    }
               for ( int p = 1; p <= s.Cols( ) - reg; p++ )
               {    vec<String> nexts0, nexts;
                    for ( int i = 0; i < s.Rows( ); i++ )
                    {    String n(reg);
                         Bool def = True;
                         for ( int j = 0; j < reg; j++ )
                         {    if ( !s.Def(i,p+j) )
                              {    def = False;
                                   break;    }
                              n[j] = s.Base(i,p+j);    }
                         if (def) nexts0.push_back(n);    }
                    Sort(nexts0);
                    int maxc = 0;
                    for ( int i = 0; i < nexts0.isize( ); i++ )
                    {    int j = nexts0.NextDiff(i);
                         maxc = Max( maxc, j - i );
                         i = j - 1;    }
                    for ( int i = 0; i < nexts0.isize( ); i++ )
                    {    int j = nexts0.NextDiff(i);
                         if ( j - i >= maxc/4 ) nexts.push_back( nexts0[i] );
                         i = j - 1;    }
                    vec<String> walks2;
                    vec<int> support2;
                    for ( int w = 0; w < walks.isize( ); w++ )
                    {    String x = walks[w];
                         x.resize( x.size( ) + 1 );
                         vec<char> exts;
                         for ( int i = 0; i < s.Rows( ); i++ )
                         {    Bool ext = True;
                              for ( int j = 0; j < reg - 1; j++ )
                              {    if ( s.Base(i,p+j) != x[p+j] )
                                   {    ext = False;
                                        break;    }    }
                              if ( !s.Def( i, p+reg-1 ) ) ext = False;
                              if (ext) exts.push( s.Base( i, p+reg-1 ) );    }
                         Sort(exts);
                         int max_support = 0;
                         for ( int i = 0; i < exts.isize( ); i++ )
                         {    int j = exts.NextDiff(i);
                              max_support = Max( max_support, j - i );
                              i = j - 1;    }
                         for ( int i = 0; i < exts.isize( ); i++ )
                         {    int j = exts.NextDiff(i);
                              x.back( ) = exts[i];
                              String n = x.substr( x.isize( ) - reg, reg );
                              if ( !BinMember( nexts, n ) ) continue;
                              if ( j - i >= max_support/2 ) 
                              {    walks2.push_back(x);    
                                   support2.push_back( support[w] + j - i );    }
                              i = j - 1;    }    }
                    walks = walks2;
                    PRINT2_TO( gout, p, walks.size( ) );
                    if ( walks.isize( ) <= 20 )
                    {    for ( int i = 0; i < walks.isize( ); i++ )
                         {    gout << "[" << i+1 << "] "; 
                              for ( int j = 0; j < walks[i].isize( ); j++ )
                                   gout << as_base( walks[i][j] );
                              gout << endl;    }    }
                    support = support2;    }
               gout << "found " << walks.size( ) << " consistent closures, "
                    << "support = " << printSeq(support)
                    << ", used " << TimeSince(cclock) << endl;    
               }

               // Save.

               bseq.push_back(f);
               if (VERBOSE)
               {    String title = "e1=" + ToString(e1) + ",e2=" + ToString(e2)
                         + ",gap=" + ToString(gap) + ( good ? "" : ",BAD" );
                    f.Print( out, title );
                    if (SHOW_MERGED_STACKS) 
                    {    gout << "\n" << title << endl;
                         s.Print(gout);    }    }    }    }
     /*
     if (VERBOSE)
     {    {    Ofstream( fout, "whatever.fasta" );
               fout << out.str( );    }
          SystemSucceed( "QueryLookupTable K=12 MM=12 MC=0.15 SEQS=whatever.fasta "
               "L=whatever/genome.lookup VISUAL=True NH=True "
               "QUIET=True MF=5K SMITH_WAT=True QUERY_NAMING=from_record > whatever" 
               );    
          Remove( "whatever.fasta" );
          fast_ifstream in( "whatever" );
          String line;
          while(1)
          {    getline( in, line );
               if ( in.fail( ) ) break;
               gout << line << "\n";    }
          Remove( "whatever" );    }
     */
     }
     gout << Date( ) << ": " << TimeSince(clock3) << " used on segment 3" << endl;

     // Create merged graph.

     double clock4 = WallClockTime( );
     gout << Date( ) << ": creating merged graph" << endl;
     HyperBasevector hbl2;
     vecbasevector bork;
     for ( int e = 0; e < hbl.E( ); e++ )
          bork.push_back( hbl.O(e) );
     for ( int v = 0; v < hbl.N( ); v++ )
     for ( int j1 = 0; j1 < hbl.To(v).isize( ); j1++ )
     for ( int j2 = 0; j2 < hbl.From(v).isize( ); j2++ )
     {    int e1 = hbl.ITo(v,j1), e2 = hbl.IFrom(v,j2);
          basevector b = hbl.Cat( e1, e2 );
          bork.push_back(b);    }
     for ( int i = 0; i < bseq.isize( ); i++ )
          bork.push_back( bseq[i] );
     bork.Append(edgesy);
     BasesToGraph( bork, K, hbl2 );

     // Clean the graph.

     gout << Date( ) << ": cleaning the graph" << endl;
     int b1 = -1, b2 = -1;
     String s1 = edgesy[0].ToString( ), s2 = edgesy[1].ToString( );
     s1.resize(K);
     s2 = s2.substr( s2.isize( ) - K, K );
     vec<int> left2, right2;
     hbl2.ToLeft(left2), hbl2.ToRight(right2);
     for ( int e = 0; e < hbl2.E( ); e++ )
     {    int p1 = hbl2.O(e).ToString( ).Position(s1);
          if ( p1 >= 0 )
          {    b1 = e;    
               hbl2.SplayVertex( left2[e] );
               hbl2.OMutable(e).SetToSubOf( hbl2.O(e), p1, hbl2.Bases(e) - p1 );    }
          int p2 = hbl2.O(e).ToString( ).Position(s2);
          if ( p2 >= 0 )
          {    b2 = e;    
               hbl2.SplayVertex( right2[e] );    
               if ( p2 + K < hbl2.Bases(e) )
                    hbl2.OMutable(e).resize( p2 + K );    }    }
     hbl2.ToLeft(left2), hbl2.ToRight(right2);
     if ( b1 < 0 || b2 < 0 )
     {    gout << Date( ) << ": can't find ends, giving up" << endl;
          gout << "\nCONCLUSION: there are 0 closures for ";
          if ( R > 0 ) gout << "R = " << R << endl;
          else gout << "S = " << S << endl;
          if ( VERBOSITY == 0 ) 
          {   
               #pragma omp critical
               {    MakeDots( done, ndots, ngaps );    }    }
          if ( !DIRECT && VERBOSITY > 0 )
          {
               #pragma omp critical
               {    cout << goutx.str( );    }    }
          continue;    }
     int v = right2[b1], w = left2[b2];
     vec<int> bx;
     if ( b1 == b2 ) bx = {b1};
     else
     {    v = right2[b1], w = left2[b2];
          if ( v == w )
          {    bx.push_back( b1, b2 );    }
          else
          {    bx = hbl2.EdgesSomewhereBetween( v, w );
               if ( bx.empty( ) )
               {    gout << "\nWARNING(1): COULDN'T FIND ANY EDGES BETWEEN!!" 
                         << endl << endl;    }
               else bx.push_back( b1, b2 );    }    }
     UniqueSort(bx);
     vec<int> dels2;
     for ( int e = 0; e < hbl2.E( ); e++ )
          if ( !BinMember( bx, e ) ) dels2.push_back(e);
     hbl2.DeleteEdges(dels2);
     hbl2.RemoveUnneededVertices( );
     hbl2.RemoveDeadEdgeObjects( );
     hbl2.RemoveEdgelessVertices( );
     gout << Date( ) << ": merged graph has " << hbl2.E( ) << " edges" << endl;
     /*
     if (VERBOSE)
     {    Ofstream( out, "whatever1.dot" );
          hbl2.PrintSummaryDOT0w( out, False, False, True );    }
     */

     // Again again find edge support.

     // NOTE EXPENSIVE CONVERSION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     GetSupport( hb, HyperBasevectorX(hbl2), basesy, support );
     if (VERBOSE)
     {    gout << "\nMERGED GRAPH" << endl << endl;
          for ( int e = 0; e < hbl2.E( ); e++ )
          {    hbl2.O(e).Print( gout, 
                    ToString(e) + "[" + ToString( support[e].size( ) ) + "]" );    }
          gout << endl;    }

     // Delete weak edges.  The second pass deletes some edges that match ddn,
     // and which might be true.

     for ( int pass = 1; pass <= 2; pass++ )
     {    dels.clear( );
          for ( int v = 0; v < hbl2.N( ); v++ )
          {    if ( hbl2.From(v).size( ) <= 1 ) continue;
               vec<int> s( hbl2.From(v).size( ) );
               for ( int j = 0; j < hbl2.From(v).isize( ); j++ )
                    s[j] = support[ hbl2.IFrom(v,j) ].size( );
               vec<int> ids( hbl2.From(v).size( ), vec<int>::IDENTITY );
               ReverseSortSync( s, ids );
               for ( int j = 1; j < s.isize( ); j++ )
               {    int e = hbl2.IFrom( v, ids[j] );
                    if ( s[0] >= 3 && s[j] == 0 ) dels.push_back(e);
                    if ( s[0] >= 6 && s[j] == 1 ) dels.push_back(e);
                    if ( s[0] >= 8 && s[j] == 2 ) dels.push_back(e);    }    }
          for ( int v = 0; v < hbl2.N( ); v++ )
          {    if ( hbl2.To(v).size( ) <= 1 ) continue;
               vec<int> s( hbl2.To(v).size( ) );
               for ( int j = 0; j < hbl2.To(v).isize( ); j++ )
               s[j] = support[ hbl2.ITo(v,j) ].size( );
               vec<int> ids( hbl2.To(v).size( ), vec<int>::IDENTITY );
               ReverseSortSync( s, ids );
               for ( int j = 1; j < s.isize( ); j++ )
               {    int e = hbl2.ITo( v, ids[j] );
                    if ( s[0] >= 3 && s[j] == 0 ) dels.push_back(e);
                    if ( s[0] >= 6 && s[j] == 1 ) dels.push_back(e);
                    if ( s[0] >= 8 && s[j] == 2 ) dels.push_back(e);    }    }
          hbl2.DeleteEdges(dels);
          hbl2.RemoveUnneededVertices( );
          hbl2.RemoveDeadEdgeObjects( );
          hbl2.RemoveEdgelessVertices( );
          if ( pass == 2 ) 
               // NOTE EXPENSIVE CONVERSION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               GetSupport( hb, HyperBasevectorX(hbl2), basesy, support );    }

     // And clean again.

     b1 = -1, b2 = -1;
     for ( int e = 0; e < hbl2.E( ); e++ )
     {    if ( hbl2.O(e).ToString( ).Contains(s1) ) b1 = e;
          if ( hbl2.O(e).ToString( ).Contains(s2) ) b2 = e;    }
     if ( b1 < 0 || b2 < 0 )
     {    gout << "No edges between, so nothing more to do." << endl << endl;
          gout << "\nCONCLUSION: there are 0 closures for ";
          if ( R > 0 ) gout << "R = " << R << endl;
          else gout << "S = " << S << endl;    }
     else
     {    hbl2.ToLeft(left2), hbl2.ToRight(right2);
          if ( b1 == b2 ) bx = {b1};
          else
          {    v = right2[b1], w = left2[b2];
               bx = hbl2.EdgesSomewhereBetween( v, w );
               if ( bx.empty( ) )
               {    gout << "\nWARNING(2): COULDN'T FIND ANY EDGES BETWEEN!!" 
                         << endl << endl;    }
               else bx.push_back( b1, b2 );    }
          UniqueSort(bx);
          dels2.clear( );
          for ( int e = 0; e < hbl2.E( ); e++ )
               if ( !BinMember( bx, e ) ) dels2.push_back(e);
          hbl2.DeleteEdges(dels2);
          hbl2.RemoveUnneededVertices( );
          hbl2.RemoveDeadEdgeObjects( );
          hbl2.RemoveEdgelessVertices( );
     
          // Print the graph.
     
          gout << "\nfinal local graph has " << hbl2.E( ) << " edges";
          /*
          if (VERBOSE) gout << ", below, and printing as whatever.dot";
          gout << endl;
          if (VERBOSE)
          {    Ofstream( out, "/whatever/whatever.dot" );
               hbl2.PrintSummaryDOT0w( out, False, False, True );    }
          // NOTE EXPENSIVE CONVERSION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          GetSupport( hb, HyperBasevectorX(hbl2), basesy, support );
          if (VERBOSE)
          {    for ( int e = 0; e < hbl2.E( ); e++ )
               {    hbl2.O(e).Print( gout, ToString(e) + "[" 
                         + ToString( support[e].size( ) ) + "]" );    }
               gout << endl;    }
          */

          if ( hbl2.E( ) == 0 )
          {    gout << "\nCONCLUSION: there are 0 closures for ";
               if ( R > 0 ) gout << "R = " << R << endl;
               else gout << "S = " << S << endl;    }

          // Get some stats about the graph.

          if ( hbl2.E( ) > 0 )
          {    hbl2.ToLeft(left2), hbl2.ToRight(right2);
               v = -1, w = -1;
               for ( int x = 0; x < hbl2.N( ); x++ )
               {    if ( hbl2.Source(x) ) v = x;
                    if ( hbl2.Sink(x) ) w = x;    }
               vec<vec<int>> paths;
               const int max_paths = 200;
               const int max_iterations = 1000000;
               Bool OK = hbl2.EdgePaths( left2, right2, v, w, paths, -1, max_paths,
                    max_iterations );
               if ( !OK ) gout << "path calculation blew up" << endl;
     
               gout << "CONCLUSION: there are ";
               if (OK) gout << paths.size( );
               else gout << "unknown number of";
               gout << " closures for ";
               if ( R > 0 ) gout << "R = " << R << endl;
               else gout << "S = " << S << endl;
               if ( !OK || paths.nonempty( ) )
               {    if ( OK && paths.isize( ) <= MAX_CLOSURES )
                    {    nclosed++;
                         nclosures += paths.size( );    }
                    else semiclosed++;    }

               // Create closure sequences.

               if ( OK && paths.isize( ) <= MAX_CLOSURES )
               {    for ( int i = 0; i < paths.isize( ); i++ )
                    {    basevector b = hbl2.Cat( paths[i] );
                         bpaths[ng].push_back(b);    }    }

               // Check to see if closures exactly match DDN.  This implementation
               // will not follow multiple branches that agree on the "first" base,
               // but these should be rare.

               if ( OK && paths.isize( ) <= MAX_CLOSURES && VERBOSE )
               {    gout << "searching ddn" << endl;
                    double clock = WallClockTime( );

                    const int KK = 200;
                    vec<basevector> bstarts;
                    if (OK)
                    {    for ( int i = 0; i < paths.isize( ); i++ )
                         {    basevector b = bpaths[ng][i];
                              if ( b.isize( ) >= KK )
                              {    b.resize(KK);
                                   bstarts.push_back(b);    }    }
                         UniqueSort(bstarts);    }
                    vec< pair<int,int> > starts;
                    #pragma omp parallel for schedule( dynamic, 1000 )
                    for ( int u = 0; u < ddn.E( ); u++ )
                    {    const basevector& U = ddn.O(u);
                         for ( int j = 0; j <= U.isize( ) - KK; j++ )
                         for ( int i = 0; i < bpaths[ng].isize( ); i++ )
                         {    Bool mismatch = False;
                              for ( int l = 0; l < KK; l++ )
                              {    if ( bstarts[i][l] != U[j+l] )
                                   {    mismatch = True;
                                        break;    }    }
                              if (mismatch) continue;
                              #pragma omp critical
                              {    starts.push( u, j );    }    }    }
                    vec<Bool> found( paths.size( ), False );
                    for ( int i = 0; i < paths.isize( ); i++ )
                    {    basevector b = hbl2.Cat( paths[i] );
                         if ( b.isize( ) < KK ) continue;
                         for ( int s = 0; s < starts.isize( ); s++ )
                         {    int u = starts[s].first, pos = starts[s].second;
                              Bool mismatch = False;
                              int bpos;
                              for ( bpos = 0; 
                                   bpos < Min( b.isize( ), ddn.Bases(u) - pos ); 
                                   bpos++ )
                              {    if ( b[bpos] != ddn.O(u)[pos+bpos] )
                                   {    mismatch = True;
                                        break;    }    }
                              if (mismatch) continue;
                              while( bpos < b.isize( ) && !mismatch )
                              {    int v = ddn.ToRight(u);
                                   Bool match = False;
                                   for ( int j = 0; 
                                        j < (int) ddn.From(v).size( ); j++ )
                                   {    int up = ddn.IFrom( v, j );
                                        if ( ddn.O(up)[KK-1] != b[bpos] ) continue;
                                        int upos = KK;
                                        u = up;
                                        bpos++;
                                        Bool bad = False;
                                        while( bpos < b.isize( )
                                             && upos < ddn.Bases(u) )
                                        {    if ( ddn.O(up)[upos] != b[bpos] )
                                             {    bad = True;
                                                  break;    }
                                             bpos++; upos++;    }
                                        if ( !bad )
                                        {    match = True;
                                             break;    }    }
                                   if ( !match ) 
                                   {    mismatch = True;
                                        break;    }    }
                              if ( !mismatch )
                              {    found[i] = True;
                                   break;    }    }    }
                    gout << "found " << Sum(found) << " closures in ddn" << endl;
                    if ( paths.isize( ) <= MAX_CLOSURES )
                         nclosures_true += Sum(found);
                    gout << "searching took " << TimeSince(clock) << endl;    }

               // Align the paths to the reference.
     
               /*
               if ( VERBOSE && paths.size( ) <= max_paths )
               {    {    Ofstream( out, "bpaths.fasta" );
                         for ( int i = 0; i < paths.isize( ); i++ )
                         {    ostringstream tout;
                              tout << printSeq( paths[i] );
                              hbl2.Cat( paths[i] ).Print( out, 
                                   "[" + ToString(i) + "] = " 
                                   + tout.str( ) );    }    }
                    if (VERBOSE)
                    {    gout << "\nALIGNMENT OF PATHS FOR ";
                         if ( R > 0 ) gout << "R = " << R << endl;
                         else gout << "S = " << S << endl;
                         SystemSucceed( "QueryLookupTable K=12 MM=12 MC=0.15 "
                              "SEQS=bpaths.fasta "
                              "L=whatever/genome.lookup "
                              "VISUAL=True NH=True QUIET=True MF=5K SMITH_WAT=True "
                              "QUERY_NAMING=from_record > whatever" );
                         CpAppend( "whatever", gout );
                         if ( paths.isize( ) <= MAX_CLOSURES )
                         {    fast_ifstream in( "whatever" );
                              while(1)
                              {    String line;
                                   getline( in, line );
                                   if ( in.fail( ) ) break;
                                   if ( !line.Contains( "mismatches" ) ) continue;
                                   int mismatches 
                                        = line.Between( ", ", " mismatches" ).Int( );
                                   int indels 
                                        = line.Between( "/", " indels" ).Int( );
                                   if ( mismatches + indels > 25 ) 
                                   {    gout << "see ugly closure" << endl;
                                        ugly++;    }    }    }
                         Remove( "whatever" );    }    }
               */
               if ( VERBOSE && OK )
               {    gout << "\nanalyzing closures" << endl;
                    for ( int i = 0; i < paths.isize( ); i++ )
                    {    basevector c = hbl2.Cat( paths[i] );
                         const int K = 48;
                         ForceAssertEq( K, hb.K( ) );
                         vec< triple<kmer<K>,int,int> > kmers_plus;
                         vecbasevector allc = basesy;
                         allc.push_back(c);
                         MAKE_KMER_LOOKUP( allc, kmers_plus );
                         vec< pair<int,int> > places;
                         for ( int j = 0; j < kmers_plus.isize( ); j++ )
                         {    int k;
                              for ( k = j + 1; k < kmers_plus.isize( ); k++ )
                              {    if ( kmers_plus[k].first != kmers_plus[j].first )
                                        break;    }
                              for ( int m1 = j; m1 < k; m1++ )
                              for ( int m2 = j; m2 < k; m2++ )
                              {    int id = kmers_plus[m1].second;
                                   if ( id == (int) bases.size( ) ) continue;
                                   if ( kmers_plus[m2].second != (int) bases.size() )
                                        continue;
                                   int offset = kmers_plus[m2].third
                                        - kmers_plus[m1].third;
                                   places.push( id, offset );    }
                              j = k - 1;    }
                         UniqueSort(places);
                         vec<Bool> to_delete( places.size( ), False );
                         for ( int j = 0; j < places.isize( ); j++ )
                         {    int k;
                              for ( k = j + 1; k < places.isize( ); k++ )
                                   if ( places[k].first != places[j].first ) break;
                              if ( k - j > 1 )
                              {    for ( int l = j; l < k; l++ )
                                        to_delete[l] = True;    }    }
                         EraseIf( places, to_delete );

                         // NOW ANALYZE SUPPORT
               
                         }    }    }

          // Print graph.

          if (VERBOSE)
          {    gout << "FINAL GRAPH" << endl << endl;
               for ( int e = 0; e < hbl2.E( ); e++ )
               {    hbl2.O(e).Print( cout, ToString(e) + "[" 
                         + ToString( support[e].size( ) ) + "]" );    }
               gout << endl;    }    }    

     gout << Date( ) << ": " << TimeSince(clock4) << " used on segment 4" << endl;
     if ( VERBOSITY == 0 )
     {
          #pragma omp critical
          {    MakeDots( done, ndots, ngaps );    }    }
     if ( !DIRECT && VERBOSITY > 0 )
     {
          #pragma omp critical
          {    cout << goutx.str( );    }    }

          }    }

     // Prep and insert the gap closures.

     cout << Date( ) << ": prepping gap closures" << endl;
     int nfails = 0;
     vec<int> dels;
     // #pragma omp parallel for num_threads(nthreads) schedule(dynamic, 1)
     for ( int ng = 0; ng < ngaps; ng++ )
     {    int S = thegap[ng];
          vec<basevector>& b = bpaths[ng];
          if ( b.empty( ) ) continue;
          int v = to_left[S], w = to_right[S];
          int d1 = D.ITo(v,0), d2 = D.IFrom(w,0);
          int start = Max( 0, dkmers[d1] + K - 1 - USELEN );
          vecbasevector edgesx;
          edgesx.push_back( hb.Cat( D.O(d1) ) );
          edgesx.push_back( hb.Cat( D.O(d2) ) );
          basevector E1 = edgesx[0], E2 = edgesx[1];
          int n1 = E1.isize( ), n2 = E2.isize( );
          int stop = Min( USELEN, n2 );
          int start0 = start, stop0 = stop;
          for ( ; start < n1; start++ )
          {    Bool mismatch = False;
               for ( int j = 0; j < b.isize( ); j++ )
               {    if ( start-start0 == b[j].isize( ) 
                         || b[j][start-start0] != E1[start] ) 
                    {    mismatch = True;
                         break;    }    }
               if (mismatch) break;    }
          start -= (K-1);

          // To find the match interval, take the K-1 bases beginning at start
          // on E1, and the first K-1 bases on b, after trimming start-start0
          // bases from its left.

          for ( ; stop > 0; stop-- )
          {    Bool mismatch = False;
               for ( int j = 0; j < b.isize( ); j++ )
               {    if ( stop0-stop == b[j].isize( )
                         || b[j][ b[j].isize( ) - (stop0-stop) - 1 ] != E2[stop-1] )
                    {    mismatch = True;
                         break;    }    }
               if (mismatch) break;    }
          stop += (K-1);

          // To find the match interval, take the K-1 bases ending at stop on E2,
          // and the last K-1 bases on b, after tripping stop - (K-1) bases from
          // its right.

          int ltrim = start - start0, rtrim = stop0 - stop;
          /*
          PRINT2( ltrim, rtrim );
          cout << endl;
          if ( start < start0 ) cout << "problem" << endl;
          E1.Print( cout, "d1" );
          E2.Print( cout, "d2" );
          for ( int j = 0; j < b.isize( ); j++ )
               b[j].Print( cout, j );
          */

          int trim1 = E1.isize( ) - ( start + K - 1 );
          ForceAssertGe( trim1, 0 );
          E1.resize( E1.isize( ) - trim1 );
          ForceAssertLe( stop, E2.isize( ) );
          int trim2 = stop - (K-1);
          E2.SetToSubOf( E2, trim2, E2.isize( ) - trim2 );
          ForceAssert( E1.isize( ) >= K-1 );
          ForceAssert( E2.isize( ) >= K-1 );
          // PRINT2( trim1, trim2 );

          Bool OK = True;
          for ( int j = 0; j < b.isize( ); j++ )
          {    if ( !( b[j].isize( ) - ltrim - rtrim >= K - 1 ) )
               {    OK = False;
                    break;    }
               b[j].SetToSubOf( b[j], ltrim, b[j].isize( ) - ltrim - rtrim );    }
          if ( !OK )
          {    // cout << "need different representation of patch" << endl;
               // cout << "probably there are shared base edges between left and "
               //      << "right" << endl;
               vec<int> left = D.O(d1), right = D.O(d2);
               Sort(left), Sort(right);
               vec<int> lr = Intersection( left, right );
               int share = 0;
               for ( int i = 0; i < lr.isize( ); i++ )
                    share += hb.Kmers( lr[i] );
               // cout << "total shared kmers = " << share << endl;
               #pragma omp critical
               {    nfails++;    }
               continue;    }

          for ( int j = 0; j < b.isize( ); j++ )
               ForceAssert( b[j].isize( ) >= (K-1) );

          /*
          E1.Print( cout, "d1_trim" );
          E2.Print( cout, "d2_trim" );
          for ( int j = 0; j < b.isize( ); j++ )
               b[j].Print( cout, ToString(j) + "_trim" );
          */

          for ( int j = 0; j < b.isize( ); j++ )
          for ( int i = 0; i < K - 1; i++ )
          {    if ( b[j][i] != E1[ E1.isize( ) - K + 1 + i ] )
               {    // cout << "left mismatch" << endl;
                    ForceAssert( 0 == 1 );
                    break;    }
               if ( b[j][ b[j].isize( ) - i - 1 ] != E2[ K - 2 - i ] )
               {    // cout << "right mismatch" << endl;
                    ForceAssert( 0 == 1 );
                    break;    }    }

          vec<int> n(2, 0);
          for ( int i = 0; i < 2; i++ )
          {    const int d = ( i == 0 ? d1 : d2 );
               for ( int j = 0; j < D.O(d).isize( ); j++ )
                    n[i] += hb.Kmers( D.O(d)[j] );
               int trim = ( i == 0 ? trim1 : trim2 );
               int v = to_left[d], w = to_right[d];
               int ltrim, rtrim;
               if ( D.To(v).nonempty( ) && IsSequence( D.OTo(v,0) ) )
               {    basevector x;
                    GapToSeq( D.OTo(v,0), ltrim, rtrim, x );
                    trim += rtrim;    }
               if ( D.From(w).nonempty( ) && IsSequence( D.OFrom(w,0) ) )
               {    basevector x;
                    GapToSeq( D.OFrom(w,0), ltrim, rtrim, x );
                    trim += ltrim;    }
               n[i] -= trim;    }
          if ( n[0] < 1 || n[1] < 1 )
          {
               #pragma omp critical
               {    cout << "FAIL" << endl;
                    nfails++;    }
               continue;    }

          // Patch the gap.

          int rv = to_left[ dinv[S] ], rw = to_right[ dinv[S] ];
          for ( int j = 0; j < b.isize( ); j++ )
          {    vec<int> x;
               SeqToGap( trim1, trim2, b[j], x );
               int N = D.E( );
               D.AddEdge( v, w, x );
               b[j].ReverseComplement( );
               SeqToGap( trim2, trim1, b[j], x );
               D.AddEdge( rv, rw, x );
               dinv.push_back( N+1, N );
               dels.push_back( S, dinv[S] );    }    }

     // Clean up.

     DPRINT(nfails);
     cout << Date( ) << ": cleaning up" << endl;
     D.DeleteEdges(dels);
     cout << Date( ) << ": cleaning up" << endl;
     CleanupCore( D, dinv );
     Validate( hb, inv, D, dinv );

     // Summary stats.

     if ( !ALL )
     {    cout << "\nSUMMARY STATS" << endl;
          cout << "- closed " << nclosed << " of " << ngaps << " gaps" << endl;
          cout << "- plus " << semiclosed << " unusable (more than MAX_CLOSURES = " 
               << MAX_CLOSURES << " closures)" << endl;
          cout << "- total usable closures = " << nclosures;
          if (VERBOSE)
          {    cout << " of which " << nclosures_true << " match ddn and " 
                    << ugly << " are ugly";    }
          cout << endl;
          cout << "- used " << TimeSince(clock) << endl << endl;
          cout << Date( ) << ": done, " << TimeSince(mclock)
               << " used on main loop" << endl << endl;    }
     else 
     {    cout << Date( ) << ": post-patching complete" << endl;
          cout << Date( ) << ": closed " << nclosed << " of " << ngaps 
               << " gaps" << endl;
          cout << Date( ) << ": plus " << semiclosed 
               << " unusable (more than MAX_CLOSURES = " 
               << MAX_CLOSURES << " closures)" << endl;
          cout << Date( ) << ": total usable closures = " << nclosures
               << endl;    }    }
