///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "Intvector.h"
#include "ParallelVecUtilities.h"
#include "VecUtilities.h"
#include "efasta/EfastaTools.h"
#include "math/Functions.h"
#include "graph/FindCells.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "paths/long/large/GapToyTools.h"
#include "paths/long/large/Lines.h"
#include "paths/long/large/CN1PeakFinder.h"
#include "system/SortInPlace.h"

template <class DE>
vec<int> OrderEdges( const DE& dgraph ) {
     int nobj = dgraph.EdgeObjectCount( );
     vec<int> ids( nobj, vec<int>::IDENTITY );
     std::sort(ids.begin(),ids.end(),
             [&dgraph](int i1,int i2)
             {return dgraph.EdgeObject(i1).size()>dgraph.EdgeObject(i2).size();});
     return ids;
}

template <>
vec<int> OrderEdges( const digraphE<int>& dgraph ) {
     int nobj = dgraph.EdgeObjectCount( );
     vec<int> ids( nobj, vec<int>::IDENTITY );
     return ids;
}

template <class DE> class GraphRuler;

template <>
struct GraphRuler<HyperBasevector>
{
     static bool isGap( HyperBasevector const& hb, int edge ) {
          return hb.EdgeLengthBases(edge) == 0;
     }
};

template <>
struct GraphRuler<digraphE<basevector>>
{
     static bool isGap( digraphE<basevector> const& hb, int edge ) {
          return False;
     }
};

template <>
struct GraphRuler<digraphE<vec<int>>>
{
     static bool isGap( digraphE< vec<int> > const& super, int edge ) {
          return super.O(edge)[0] < 0;
     }
};

template <>
struct GraphRuler<digraphE<int>>
{
     static bool isGap( digraphE<int> const& super, int edge ) {
          return super.O(edge) < 0;
     }
};

// WARNING: the single option to FindLines is only partially implemented.

template <class DE>
void FindLines( const DE& dgraph, const vec<int>& dinv,
     vec<vec<vec<vec<int>>>>& lines, const int max_cell_paths, const int max_depth,
     const Bool verbose, const Bool single )
{    double clock = WallClockTime( );
     vec<int> to_left, to_right;
     dgraph.ToLeft(to_left), dgraph.ToRight(to_right);

     // Heuristics.

     const int verts_mul = 2;

     // Find some cells.  These do not include standard bubbles.

     vec<vec<vec<int>>> xpaths;
     vec< pair<int,int> > bounds;
     {    int max_cell_verts = verts_mul * max_cell_paths;
          vec< pair<int,int> > bounds0;
          if (verbose) cout << Date( ) << ": finding cells" << endl;
          FindSomeCells( dgraph, max_cell_verts, max_depth, bounds0 );

          // Symmetrize cells.

          if (verbose) cout << Date( ) << ": symmetrizing cells" << endl;
          int nb = bounds0.size( );
          for ( int i = 0; i < nb; i++ )
          {    int v = bounds0[i].first, w = bounds0[i].second;
               int rv = to_right[ dinv[ dgraph.IFrom(v,0) ] ];
               int rw = to_left[ dinv[ dgraph.ITo(w,0) ] ];
               bounds0.push( rw, rv );    }
          if (single) UniqueSort(bounds0);
          else ParallelUniqueSort(bounds0);

          // Find paths across cells.

          xpaths.resize( bounds0.size( ) );
          bounds = bounds0;
          vec<Bool> xdel( xpaths.size( ), False );
          #pragma omp parallel for
          for ( int i = 0; i < bounds0.isize( ); i++ )
          {    int v = bounds0[i].first, w = bounds0[i].second;
               Bool OK = dgraph.EdgePaths( to_left, to_right, v, w, xpaths[i], -1,
                    max_cell_paths );
               int nxpaths = xpaths[i].size( );
               if ( !OK || nxpaths > max_cell_paths ) xdel[i] = True;    }
          EraseIf( xpaths, xdel ), EraseIf( bounds, xdel );    }

     // Remove subset cells.

     if (verbose) cout << Date( ) << ": removing subset cells" << endl;
     int nobj = dgraph.EdgeObjectCount( );
     vec<Bool> used;
     dgraph.Used(used);
     VecIntVec contents( bounds.size( ) );
     IntVec e;
     for ( int i = 0; i < bounds.isize( ); i++ )
     {    int v = bounds[i].first, w = bounds[i].second;
          e.clear();
          e.push_back(dgraph.IFrom(v,0)).push_back(dgraph.ITo(w,0));
          for ( int j = 0; j < xpaths[i].isize( ); j++ )
            for ( int k = 0; k < xpaths[i][j].isize( ); k++ )
               e.push_back( xpaths[i][j][k] );
          std::sort(e.begin(),e.end());
          e.erase(std::unique(e.begin(),e.end()),e.end());
          contents[i] = e;    }
     VecIntVec cell_index;
     invert(contents,cell_index,dgraph.EdgeObjectCount());
     vec<Bool> xdel2( bounds.size( ), False );
     for ( int e = 0; e < dgraph.EdgeObjectCount( ); e++ )
     for ( unsigned j1 = 0; j1 < cell_index[e].size( ); j1++ )
     {    int c1 = cell_index[e][j1];
          if ( xdel2[c1] ) continue;
          for ( unsigned j2 = 0; j2 < cell_index[e].size( ); j2++ )
          {    if ( j1 == j2 ) continue;
               int c2 = cell_index[e][j2];
               if ( xdel2[c2] ) continue;
               if ( bounds[j1].second == bounds[j2].first ) continue;
               if ( bounds[j2].second == bounds[j1].first ) continue;
               if ( contents[c1].size( ) >= contents[c2].size( ) ) continue;
               if ( BinSubset( contents[c1], contents[c2] ) )
                    xdel2[c1] = True;    }    }
     cell_index.clear();
     EraseIf( bounds, xdel2 );
     EraseIf( xpaths, xdel2 );

     // Add in gaps.

     if (verbose) cout << Date( ) << ": adding in gaps" << endl;
     //vec<vec<int>> with a single empty element.
     vec<vec<int>> x(1);
     for ( int e = 0; e < nobj; e++ )
     {    if ( !used[e] ) continue;
          int v = to_right[e];
          if ( !dgraph.To(v).solo( ) || !dgraph.From(v).solo( ) ) continue;
          int f = dgraph.IFrom( v, 0 ), w = dgraph.From(v)[0];
          if ( !GraphRuler<DE>::isGap( dgraph, f) ) continue;
          if ( !dgraph.To(w).solo( ) || !dgraph.From(w).solo( ) ) continue;
          bounds.push( v, w );
          xpaths.push_back(x);    }

     // Index bounds.

     if (verbose) cout << Date( ) << ": indexing bounds" << endl;
     if (single) SortSync( bounds, xpaths );
     else ParallelSortSync( bounds, xpaths );
     vec< vec<int> > left_ind( dgraph.N( ) ), right_ind( dgraph.N( ) );
     for ( int i = 0; i < bounds.isize( ); i++ )
     {    left_ind[ bounds[i].first ].push_back(i);
          right_ind[ bounds[i].second ].push_back(i);    }

     // Place edge indices in order by length, largest first.

     vec<int> ids = OrderEdges( dgraph );

     // Go through the edges and build lines.

     if (verbose) cout << Date( ) << ": building lines" << endl;
     vec<Bool> marked( nobj, False );
     lines.clear( );
     // #pragma omp parallel for // seems unsafe and doesn't save time
     vec<vec<vec<int>>> line;
     for ( int ie = 0; ie < nobj; ie++ )
     {    int e = ids[ie];
          if ( GraphRuler<DE>::isGap(dgraph, e) ) continue;
          if ( !used[e] ) continue;
          if ( marked[e] ) continue;

          // Build line, extending first to the left, then to the right.

          line = {{{e}}};
          Bool circle = False;
          while(1)
          {    int w = to_left[ line.back( )[0][0] ];
               if ( !dgraph.From(w).solo( ) || !right_ind[w].solo( ) ) break;
               int bid = right_ind[w][0];
               int v = bounds[bid].first;
               line.push_back( xpaths[bid] );
               int eb = dgraph.ITo(v, 0);
               line.push_back( {{eb}} );
               if ( eb == e )
               {    circle = True;
                    break;     }    }
          line.ReverseMe( );
          if ( !circle )
          {    while(1)
               {    int v = to_right[ line.back( )[0][0] ];
                    if ( !dgraph.To(v).solo( ) || !left_ind[v].solo( ) ) break;
                    int bid = left_ind[v][0];
                    int w = bounds[bid].second;
                    int eb = dgraph.IFrom(w, 0);
                    line.push_back( xpaths[bid], {{eb}} );
                    if ( eb == e ) cout << "CIRCLE!" << endl;    }    }

          // Save.

          lines.push_back(line);

          // Mark.

          for ( int i = 0; i < line.isize( ); i++ )
          for ( int j = 0; j < line[i].isize( ); j++ )
          for ( int k = 0; k < line[i][j].isize( ); k++ )
          {    int f = line[i][j][k];
               marked[f] = marked[ dinv[f] ] = True;    }

          // Generate reverse complement of line and save.

          line.ReverseMe( );
          for ( int i = 0; i < line.isize( ); i++ )
          for ( int j = 0; j < line[i].isize( ); j++ )
          {    line[i][j].ReverseMe( );
               for ( int k = 0; k < line[i][j].isize( ); k++ )
                    line[i][j][k] = dinv[ line[i][j][k] ];    }
          lines.push_back(line);    }

     // Order paths.

     if (verbose) cout << Date( ) << ": ordering paths" << endl;
     int nthreads = ( single ? 1 : omp_get_max_threads( ) );
     #pragma omp parallel for num_threads(nthreads)
     for ( int i = 0; i < lines.isize( ); i++ )
     {    for ( int j = 0; j < lines[i].isize( ); j++ )
               Sort( lines[i][j] );    }

     // Remove lines having identical content.

     if (verbose) cout << Date( ) << ": removing identical content" << endl;
     if (single) Sort(lines);
     else sortInPlaceParallel( lines.begin(), lines.end() );
     Unique(lines);
     size_t nLines = lines.size();
     contents.clear().resize(nLines);
     #pragma omp parallel for num_threads(nthreads)
     for ( size_t i = 0; i < nLines; i++ )
     {   vec<vec<vec<int>>> const& vvv = lines[i];
         size_t res = 0;
         for ( vec<vec<int>> const& vv : vvv )
             for ( vec<int> const& v : vv ) res += v.size();
         IntVec iv;
         iv.reserve(res);
         for ( vec<vec<int>> const& vv : vvv )
             for ( vec<int> const& v : vv )
                 for ( int e : v ) iv.push_back(e);
         std::sort(iv.begin(),iv.end());
         iv.erase(std::unique(iv.begin(),iv.end()),iv.end());
         contents[i] = iv;    }

     vec<int> ids2( nLines, vec<int>::IDENTITY );
     ParallelSort( ids2,
                 [&contents](int i1,int i2){return contents[i1]<contents[i2];});

     vec<Bool> to_delete1( nLines, False );
     for ( size_t i = 0; i != nLines; i++ )
     {    int m = ids2[i];
          IntVec const& probe = contents[m];
          size_t j;
          for ( j = i + 1; j != nLines; j++ )
               if ( contents[ids2[j]] != probe ) break;
          for ( size_t k = i+1; k != j; k++ )
               m = std::min( m, ids2[k] );
          for ( size_t k = i; k != j; k++ )
               if ( ids2[k] != m ) to_delete1[ ids2[k] ] = True;
          i = j - 1; }

     EraseIf( lines, to_delete1 );

     // Define line lengths in such a way that a subset line will have a smaller
     // length.

     vec< pair<int,int> > llen( lines.size( ), make_pair(0,0) );
     // OLD DEFINITION
     // FAILED ON SHOWING THAT
     // 9545651 {{}} 8965780 {{}} 8965779 {{}} 9545652
     // IS A SUBSET OF
     // 2820324 {{2820325},{2820326}} 2820327 {{2820328,11158223,2820323},
     // {9545651,11194873,8965780,11187223,8965779,11187222,9545652}} 2820324
     // THIS IS BECAUSE OF THE HETEROGENEOUS WAY WE HANDLE GAPS.
     // for ( int i = 0; i < lines.isize( ); i++ )
     // for ( int j = 0; j < lines[i].isize( ); j++ )
     //      llen[i] += lines[i][j].size( );
     // WAS IN GapToy.HG03642
     // NEW DEFINITION
     for ( int i = 0; i < lines.isize( ); i++ )
     {    llen[i].second = i;
          for ( int j = 0; j < lines[i].isize( ); j++ )
          for ( int k = 0; k < lines[i][j].isize( ); k++ )
               llen[i].first += lines[i][j][k].size( );    }

     // Remove subset lines.

     if (verbose) cout << Date( ) << ": removing subsets" << endl;
     if (single) ReverseSortSync( llen, lines );
     else ParallelReverseSortSync( llen, lines );
     vec<Bool> to_delete( lines.size( ), False );
     vec< vec<int> > lines_index(nobj);
     for ( int i = 0; i < lines.isize( ); i++ )
     for ( int j = 0; j < lines[i].isize( ); j++ )
     for ( int k = 0; k < lines[i][j].isize( ); k++ )
     for ( int l = 0; l < lines[i][j][k].isize( ); l++ )
     {    int e = lines[i][j][k][l];
          if ( lines_index[e].nonempty( ) && lines_index[e][0] == i ) continue;
          lines_index[e].push_back(i);    }
     for ( int e = 0; e < nobj; e++ )
     {    if ( lines_index[e].size( ) > 1 )
          {
               // XXX:
               if ( llen[ lines_index[e][0] ].first 
                    == llen[ lines_index[e][1] ].first )
               {    cout << "\nconflicting lines:\n";
                    for ( int j = 0; j < 2; j++ )
                    {    const vec<vec<vec<int>>>& line = lines[ lines_index[e][j] ];
                         cout << "\n";
                         for ( int i = 0; i < line.isize( ); i++ )
                         {    if ( i > 0 ) cout << " ";
                              if ( i % 2 == 0 ) cout << line[i][0][0];
                              else
                              {    cout << "{";
                                   for ( int j = 0; j < line[i].isize( ); j++ )
                                   {    if ( j > 0 ) cout << ",";
                                        cout << "{" << printSeq( line[i][j] )
                                             << "}";    }
                                   cout << "}";    }    }
                         cout << "\n";    }    }
               ForceAssertGt( llen[ lines_index[e][0] ].first, 
                    llen[ lines_index[e][1] ].first );

               for ( int j = 1; j < lines_index[e].isize( ); j++ )
                    to_delete[ lines_index[e][j] ] = True;    }    }
     EraseIf( lines, to_delete );
     if (verbose) cout << Date( ) << ": found " << lines.size( ) << " lines" << endl;
     LogTime( clock, "used finding lines" );    }

void GetTol( const HyperBasevector& hb,
     const vec<vec<vec<vec<int>>>>& lines, vec<int>& tol )
{
     tol.resize_and_set( hb.EdgeObjectCount( ), -1 );
     for ( int i = 0; i < lines.isize( ); i++ )
     for ( int j = 0; j < lines[i].isize( ); j++ )
     for ( int k = 0; k < lines[i][j].isize( ); k++ )
     for ( int l = 0; l < lines[i][j][k].isize( ); l++ )
          tol[ lines[i][j][k][l] ] = i;    }

void GetTol( const HyperBasevectorX& hb,
     const vec<vec<vec<vec<int>>>>& lines, vec<int>& tol )
{
     tol.resize_and_set( hb.E( ), -1 );
     for ( int i = 0; i < lines.isize( ); i++ )
     for ( int j = 0; j < lines[i].isize( ); j++ )
     for ( int k = 0; k < lines[i][j].isize( ); k++ )
     for ( int l = 0; l < lines[i][j][k].isize( ); l++ )
          tol[ lines[i][j][k][l] ] = i;    }

void GetLineNpairs( const HyperBasevector& hb, const vec<int>& inv,
     const ReadPathVec& paths, const vec<vec<vec<vec<int>>>>& lines,
     vec<int>& npairs )
{
     npairs.resize_and_set( lines.size( ), 0 );
     vec<int> tol;
     GetTol( hb, lines, tol );
     for ( int64_t pid = 0; pid < (int64_t) paths.size( ) / 2; pid++ )
     {    int64_t id1 = 2*pid, id2 = 2*pid+1;
          vec<int> e;
          for ( int j = 0; j < (int) paths[id1].size( ); j++ )
               e.push_back( tol[ paths[id1][j] ], tol[ inv[ paths[id1][j] ] ] );
          for ( int j = 0; j < (int) paths[id2].size( ); j++ )
               e.push_back( tol[ paths[id2][j] ], tol[ inv[ paths[id2][j] ] ] );
          UniqueSort(e);
          for ( int j = 0; j < e.isize( ); j++ )
               npairs[ e[j] ]++;    }    }

void WriteLineStats( const String& head, const vec<vec<vec<vec<int>>>>& lines,
     const vec<int>& llens, const vec<int>& npairs, const vec<vec<covcount>>& covs )
{    Ofstream( out, head + ".lines.stats" );
     vec<String> cov( lines.size( ) );
     int ns = covs.size( );
     for ( int i = 0; i < lines.isize( ); i++ )
     {    int e = lines[i][0][0][0];
          Bool defined = False;
          int ns = covs.size( );
          for ( int ss = 0; ss < ns; ss++ )
               if ( covs[ss][e].Def( ) ) defined = True;
          if (defined)
          {    for ( int ss = 0; ss < ns; ss++ )
               {    if ( ss > 0 ) cov[i] += ",";
                    if ( covs[ss][e].Def( ) )
                         cov[i] += ToString( covs[ss][e].Cov( ), 2 ) + "x";
                    else cov[i] += "?x";    }    }    }
     for ( int i = 0; i < lines.isize( ); i++ )
     {    int e1 = lines[i].front( )[0][0], e2 = lines[i].back( )[0][0];
          out << "line[" << i << "] " << e1 << ".." << e2 << " len=" << llens[i]
               << " npairs=" << npairs[i];
          if ( cov[i].nonempty( ) ) out << " cov=" << cov[i];
          out << "\n";    }    }

int64_t LineN50( const HyperBasevector& hb,
     const vec<vec<vec<vec<int>>>>& lines, const int min_len )
{    vec<int> llens, lens;
     GetLineLengths( hb, lines, llens );
     for ( int i = 0; i < lines.isize( ); i++ )
          if ( llens[i] >= min_len ) lens.push_back( llens[i] + hb.K( ) - 1 );
     if ( lens.empty( ) ) return 0;
     return N50(lens);    }

vec<int64_t> Pids( const int e, const int ss, const HyperBasevector& hb,
     const vec<int>& inv, const ReadPathVec& paths, const VecULongVec& paths_index,
     const vec<int64_t>& subsam_starts )
{    vec<int64_t> pids;
     int ns = subsam_starts.size( );
     for ( int pass = 1; pass <= 2; pass++ )
     {    int d = ( pass == 1 ? e : inv[e] );
          for ( int j = 0; j < (int) paths_index[d].size( ); j++ )
          {    int64_t pid = paths_index[d][j]/2;
               int ssx;
               for ( ssx = 0; ssx < subsam_starts.isize( ); ssx++ )
                    if ( ssx == ns - 1 || 2*pid < subsam_starts[ssx+1] ) break;
               if ( ssx == ss ) pids.push_back(pid);    }    }
     UniqueSort(pids);
     return pids;    }

vec<int64_t> Pids( const int e, const HyperBasevector& hb,
     const vec<int>& inv, const ReadPathVec& paths, const VecULongVec& paths_index )
{    vec<int64_t> pids;
     for ( int pass = 1; pass <= 2; pass++ )
     {    int d = ( pass == 1 ? e : inv[e] );
          for ( int j = 0; j < (int) paths_index[d].size( ); j++ )
          {    int64_t pid = paths_index[d][j]/2;
               pids.push_back(pid);    }    }
     UniqueSort(pids);
     return pids;    }

vec<int64_t> LinePids( const vec<vec<vec<int>>>& L, const HyperBasevector& hb,
     const vec<int>& inv, const ReadPathVec& paths, const VecULongVec& paths_index )
{    vec<int64_t> pids;
     vec<int> e;
     for ( int i1 = 0; i1 < L.isize( ); i1++ )
     for ( int i2 = 0; i2 < L[i1].isize( ); i2++ )
     for ( int i3 = 0; i3 < L[i1][i2].isize( ); i3++ )
          e.push_back( L[i1][i2][i3] );
     UniqueSort(e);
     for ( int j = 0; j < e.isize( ); j++ )
          pids.append( Pids( e[j], hb, inv, paths, paths_index ) );
     UniqueSort(pids);
     return pids;    }

double RawCoverage( const int e, const int ss, const HyperBasevector& hb,
     const vec<int>& inv, const ReadPathVec& paths, const VecULongVec& paths_index,
     const vec<int64_t>& subsam_starts )
{    vec<int64_t> pids = Pids( e, ss, hb, inv, paths, paths_index, subsam_starts );
     return pids.size( ) / double( hb.Kmers(e) );    }

double RawCoverage( const int e, const HyperBasevector& hb,
     const vec<int>& inv, const ReadPathVec& paths, const VecULongVec& paths_index )
{    vec<int64_t> pids = Pids( e, hb, inv, paths, paths_index );
     return pids.size( ) / double( hb.Kmers(e) );    }

void ComputeCoverage( const HyperBasevector& hb, const vec<int>& inv,
     const ReadPathVec& paths, const vec<vec<vec<vec<int>>>>& lines,
     const vec<int64_t>& subsam_starts, vec<vec<covcount>>& covs )
{
     // Heuristics.

     const int min_line = 1000;
     const int top_group = 50;

     // Index paths (better done outside this program).

     VecULongVec paths_index;
     invert( paths, paths_index, hb.EdgeObjectCount( ) );

     // Compute pairs touching each line.

     int ns = subsam_starts.size( );
     covs.resize(ns);
     vec<vec<int>> npairs( ns, vec<int>( lines.size( ), 0 ) );
     vec<int> tol;
     GetTol( hb, lines, tol );
     {
	 vec<int> e;
	 for ( int64_t pid = 0; pid < (int64_t) paths.size( ) / 2; pid++ ) {
	     int64_t id1 = 2*pid, id2 = 2*pid+1;
	     e.clear();
	     for ( int j = 0; j < (int) paths[id1].size( ); j++ )
		 e.push_back( tol[ paths[id1][j] ], tol[ inv[ paths[id1][j] ] ] );
	     for ( int j = 0; j < (int) paths[id2].size( ); j++ )
		 e.push_back( tol[ paths[id2][j] ], tol[ inv[ paths[id2][j] ] ] );
	     UniqueSort(e);
	     int ss;
	     for ( ss = 0; ss < ns; ss++ )
		 if ( ss == ns - 1 || 2*pid < subsam_starts[ss+1] ) break;
	     for ( int j = 0; j < e.isize( ); j++ )
		 npairs[ss][ e[j] ]++;
	 }
     }

     // Compute coverage of lines.

     vec<int> lens;
     GetLineLengths( hb, lines, lens );
     vec<vec<double>> covl( ns, vec<double>( lines.size( ), 0 ) );
     for ( int ss = 0; ss < ns; ss++ )
	 for ( int l = 0; l < lines.isize( ); l++ )
	     covl[ss][l] = double(npairs[ss][l]) / double(lens[l]);

     // Compute baseline coverage.

     vec<double> base_cov(ns);  // median coverage of 50 largest lines
     {    vec<vec<double>> top(ns);
          vec<int> lensx(lens), ids( lines.size( ), vec<int>::IDENTITY );
          ReverseSortSync( lensx, ids );
          for ( int ss = 0; ss < ns; ss++ )
	  // use largest 50 (or all in fewer) largest lines
          {    for ( int i = 0; i < Min( top_group, lensx.isize( ) ); i++ )
                    top[ss].push_back( covl[ss][ ids[i] ] );
               Sort( top[ss] );
               base_cov[ss] = Median( top[ss] );    }    }
     const int max_len = *std::max_element(lens.begin(), lens.end());  // largest contig
     const int min_len = Min(10000, max_len/2) ;  // Allow for the case that no contigs are longer than 10000
     const double radius = 0.08;
     const double cov_low = 0.5;
     const double cov_high = 1.5;
     for ( int ss = 0; ss < ns; ss++ ) {
	 vec<double> covx;  // coverage of lines longer than min_len
	 vec<int> ids;      // associated ids

	 // Ofstream(covx_out, "covx." + ToString(ss) );

	 cout << Date() << ": determining candidates" << endl;
	 for ( int i = 0; i < lines.isize( ); i++ ) {
	     if ( lens[i] >= min_len && covl[ss][i] > 0) {
		 ids.push_back(i);
		 covx.push_back( covl[ss][i] );
	     }
	 }
	 SortSync( covx, ids );
	 vec<int64_t> mass(covx.size()); // total length of lines with coverage within 'radius'

	 for ( int i = 0; i < covx.isize( ); i++ ) {
	     int64_t this_mass = 0;  // mass, total length of lines with coverage within 'radius' of covx[i]
	     this_mass += lens[ ids[i] ];
	     for ( int j = i - 1; j >= 0; j-- ) {
		 if ( covx[i] - covx[j] > radius * covx[i] ) break;
		 this_mass += lens[ ids[j] ];
	     }
	     for ( int j = i + 1; j < covx.isize( ); j++ ) {
		 if ( covx[j] - covx[i] > radius * covx[i] ) break;
		 this_mass += lens[ ids[j] ];
	     }
	     mass[i] = this_mass;
	     // covx_out << covx[i] << ", " << lens[ids[i]] << ", " << this_mass << endl;
	 }

	 // Find CN1 peak (diploid peak for diploid samples)

	 CN1PeakFinder peak_finder;
	 double cn1_cov = peak_finder.FindPeak(covx,mass);
	 if (cn1_cov == 0)
	     cout << "WARNING: Unable to find coverage peak for sample: " << ss << endl;

	 base_cov[ss] = cn1_cov;
     }

     // Set coverage based on lines.

     for ( int ss = 0; ss < ns; ss++ )
     {    covs[ss].resize( hb.E( ) );
          for ( int l = 0; l < lines.isize( ); l++ )
          {    if ( lens[l] >= min_line )
               {    for ( int j = 0; j < lines[l].isize( ); j += 2 )
                    {    int e = lines[l][j][0][0];
                         covs[ss][e].Set(
                              covl[ss][l] / base_cov[ss] );    }    }    }    }

     // Set coverage for some lines within cells.

     for ( int l = 0; l < lines.isize( ); l++ )
     for ( int j = 1; j < lines[l].isize( ); j += 2 )
     {    const vec<vec<int>>& x = lines[l][j];
          if ( x.empty( ) || x[0].empty( ) ) continue;
          vec< pair<int,int> > io;
          for ( int i = 0; i < x.isize( ); i++ )
               io.push( x[i].front( ), x[i].back( ) );
          UniqueSort(io);
          vec<int> ins, outs;
          for ( int i = 0; i < x.isize( ); i++ )
          {    ins.push_back( x[i].front( ) );
               outs.push_back( x[i].back( ) );    }
          UniqueSort(ins), UniqueSort(outs);
          if ( ins.size( ) != outs.size( ) || ins.size( ) != io.size( ) ) continue;
          if ( ins.size( ) < 2 ) continue;
          vec<vec<int>> pi( io.size( ) );
          for ( int i = 0; i < x.isize( ); i++ )
          {    pi[ BinPosition( io, make_pair( x[i].front( ), x[i].back( ) ) ) ]
                    .push_back(i);    }
          for ( int i = 0; i < io.isize( ); i++ )
          {    vec<int> lens, e;
               for ( int k = 0; k < pi[i].isize( ); k++ )
               {    int p = pi[i][k];
                    int len = 0;
                    for ( int u = 0; u < x[p].isize( ); u++ )
                         len += hb.Kmers( x[p][u] );
                    lens.push_back(len);    }
               Sort(lens);
               int len = Median(lens);
               if ( len < min_line ) continue;
               for ( int k = 0; k < pi[i].isize( ); k++ )
               {    int p = pi[i][k];
                    e.append( x[p] );    }
               UniqueSort(e);

               vec<vec<int64_t>> pids(ns);
               for ( int k = 0; k < e.isize( ); k++ )
               for ( int ss = 0; ss < ns; ss++ )
               {    vec<int64_t> x = Pids( e[k], ss, hb, inv, paths, paths_index,
                         subsam_starts );
                    pids[ss].append(x);    }

               for ( int ss = 0; ss < ns; ss++ )
               {    UniqueSort( pids[ss] );
                    double c = double( pids[ss].size( ) ) / double(len);
                    vec<vec<int>> keys;
                    for ( int k = 0; k < pi[i].isize( ); k++ )
                    {    int p = pi[i][k];
                         vec<int> y = x[p];
                         UniqueSort(y);
                         keys.push_back(y);    }
                    vec<int> z;
                    Intersection( keys, z );
                    for ( int l = 0; l < z.isize( ); l++ )
                         covs[ss][ z[l] ].Set( c / base_cov[ss] );    }    }    }

     // Set coverage for long edges.

     for ( int ss = 0; ss < ns; ss++ )
     for ( int e = 0; e < hb.E( ); e++ )
     {    if ( !covs[ss][e].Def( ) && hb.Kmers(e) >= min_line )
          {    double c =
                    RawCoverage( e, ss, hb, inv, paths, paths_index, subsam_starts );
               covs[ss][e].Set( c / base_cov[ss] );    }    }    }

// Write the line to standard out

void DisplayLine(const vec<vec<vec<int>>>& line) {
    bool first = false;
    for (const auto& paths : line) {
	if (!first)
	    cout << " -> ";
	cout << "{";
	for (const auto& path : paths) {
	    cout << " (";
	    for (int edge: path)
		cout << edge << " ";
	    cout << ") ";
	}
	cout << "}";
	first = false;
    }
    cout << " ";
}

void LineInv( const vec<vec<vec<vec<int>>>>& lines, const vec<int>& inv,
     vec<int>& linv )
{    linv.clear( );
     linv.resize( lines.size( ), -1 );
     vec< pair<int,int> > line_ends( lines.size( ) );
     #pragma omp parallel for
     for ( int i = 0; i < lines.isize( ); i++ )
     {    line_ends[i] = make_pair(
               lines[i].front( )[0][0], lines[i].back( )[0][0] );    }
     vec<int> ids( lines.size( ), vec<int>::IDENTITY );
     ParallelSortSync( line_ends, ids );
     Bool failed = False;
     #pragma omp parallel for
     for ( int i = 0; i < line_ends.isize( ); i++ )
     {    int e1 = line_ends[i].first, e2 = line_ends[i].second;
          int p = BinPosition( line_ends, make_pair( inv[e2], inv[e1] ) );
          if ( p >= 0 ) linv[ ids[i] ] = ids[p];
          else
          {
               // Check for the case where the line is a circular component and
               // is its own inverse.

               Bool ok = False;
               const vec<vec<vec<int>>>& L = lines[ ids[i] ];
               int e = inv[ L[0][0][0] ];
               for ( int j = 0; j < L.isize( ); j++ )
               for ( int k = 0; k < L[j].isize( ); k++ )
               for ( int l = 0; l < L[j][k].isize( ); l++ )
                    if ( L[j][k][l] == e ) ok = True;
               if (ok) linv[ ids[i] ] = ids[i];
               else
               {
                    #pragma omp critical
                    {    if ( !failed )
                         {    cout << "\nFailed to find rc of line " << ids[i] 
                                   << " = " << e1 << ".." << e2 
                                   << " [" << L.size( ) << "]." << endl;
                              cout << "Fatal Error." << endl;    }
                         failed = True;    }    }    }    }
     if (failed) TracebackThisProcess( );
     for ( int i = 0; i < lines.isize( ); i++ )
     {    ForceAssertGe( linv[i], 0 );
          ForceAssertEq( linv[linv[i]], i );    }    }

void SortLines( vec<vec<vec<vec<int>>>>& lines, const HyperBasevector& hb,
     const vec<int>& inv )
{    int N = lines.size( );
     vec<int> lens, ids( N, vec<int>::IDENTITY ), F(N), B(N);
     GetLineLengths( hb, lines, lens );
     for ( int i = 0; i < N; i++ )
     {    F[i] = lines[i].front( )[0][0], B[i] = lines[i].back( )[0][0];    }
     ParallelSort( ids, [&F,&B,&lens,&inv]( int i1, int i2 )
     {    return make_triple( -lens[i1], Min( F[i1], inv[ B[i1] ] ), F[i1] )
               < make_triple( -lens[i2], Min( F[i2], inv[ B[i2] ] ), F[i2] );
                    }    );
     vec<int> idsx(N);
     for ( int i = 0; i < N; i++ )
          idsx[ ids[i] ] = i;
     PermuteVec( lines, idsx );    }

void DumpLineFiles( const vec<vec<vec<vec<int>>>>& lines, const HyperBasevector& hb,
     const vec<int>& inv, const ReadPathVec& paths, const String& dir )
{
     const int gap = 100;
     const int K = hb.K( );

     vec<int> to_left, to_right;
     hb.ToLeft(to_left), hb.ToRight(to_right);
     VecULongVec paths_index;
     invert( paths, paths_index, hb.E( ) );
     Ofstream( out1, dir + "/a.lines.efasta" );
     Ofstream( out2, dir + "/a.lines.fasta" );
     for ( int i = 0; i < lines.isize( ); i++ )
     {
          // Don't print both a line and its rc.

          if ( i > 0 && lines[i-1].front( )[0][0] == inv[ lines[i].back( )[0][0] ] )
               continue;

          const vec<vec<vec<int>>>& L = lines[i];
          Bool circular1 = ( L.size( ) > 1 && L.front( )[0][0] == L.back( )[0][0] );
          Bool circular2 = ( L.solo( )
               && to_left[ L[0][0][0] ] == to_right[ L[0][0][0] ] );
          String b1, b2;
          for ( int j = 0; j < L.isize( ); j++ )
          {    if ( circular1 && j == L.isize( ) - 1 ) break;
               const vec<vec<int>>& x = L[j];
               if ( x.solo( ) && x[0].empty( ) )
               {    b1 += String( gap, 'N' ), b2 += String( gap, 'N' );    }
               else
               {
                    // Find the "most likely" path.  Note that we only consider
                    // paths entering from the left.  This asymmetry doesn't make
                    // sense.  Should do both sides.

                    int best = 0;
                    if ( j % 2 == 1 )
                    {    vec<int> cov( x.size( ), 0 );
                         int e = L[j-1][0][0];
                         for ( int l = 0; l < (int) paths_index[e].size( ); l++ )
                         {    const ReadPath& p = paths[ paths_index[e][l] ];
                              for ( int m = 0; m < (int) p.size( ); m++ )
                              {    if ( p[m] != e ) continue;
                                   vec<Bool> match( x.size( ), True );
                                   for ( int r = 0; r < x.isize( ); r++ )
                                   {    for ( int s = 0; s < x[r].isize( ); s++ )
                                        {    if ( m + 1 + s >= (int) p.size( ) )
                                                  break;
                                             if ( p[m+1+s] != x[r][s] )
                                             {    match[r] = False;
                                                  break;    }    }    }
                                   if ( Sum(match) == 1 )
                                   {    for ( int r = 0; r < x.isize( ); r++ )
                                             if ( match[r] ) cov[r]++;    }    }    }
                         int re = inv[e];
                         for ( int l = 0; l < (int) paths_index[re].size( ); l++ )
                         {    const ReadPath& q = paths[ paths_index[re][l] ];
                              vec<int> p;
                              for ( int m = q.size( ) - 1; m >= 0; m-- )
                                   p.push_back( inv[ q[m] ] );
                              for ( int m = 0; m < (int) p.size( ); m++ )
                              {    if ( p[m] != e ) continue;
                                   vec<Bool> match( x.size( ), True );
                                   for ( int r = 0; r < x.isize( ); r++ )
                                   {    for ( int s = 0; s < x[r].isize( ); s++ )
                                        {    if ( m + 1 + s >= (int) p.size( ) )
                                                  break;
                                             if ( p[m+1+s] != x[r][s] )
                                             {    match[r] = False;
                                                  break;    }    }    }
                                   if ( Sum(match) == 1 )
                                   {    for ( int r = 0; r < x.isize( ); r++ )
                                             if ( match[r] ) cov[r]++;    }    }    }
                         vec<int> ids( x.size( ), vec<int>::IDENTITY );
                         ReverseSortSync( cov, ids );
                         best = ids[0];    }

                    // Add to fasta/efasta.

                    vec<basevector> bs;
                    for ( int m = 0; m < x.isize( ); m++ )
                    {    bs.push_back( hb.Cat( x[m] ) );
                         if ( j < L.isize( ) - 1 )
                              bs.back( ).resize( bs.back( ).isize( ) - (K-1) );    }
                    b1 += efasta(bs);
                    b2 += bs[best].ToString( );    }    }
          String header = "line_" + ToString(i);
          if ( circular1 || circular2 ) header += " circular";
          efasta(b1).Print( out1, header );
          efasta(b2).Print( out2, "flattened_" + header );    }

     Ofstream( out3, dir + "/a.lines.src" );
     for ( int i = 0; i < lines.isize( ); i++ )
     {    const vec<vec<vec<int>>>& L = lines[i];
          for ( int j = 0; j < L.isize( ); j++ )
          {    if ( j > 0 ) out3 << ",";
               if ( j % 2 == 0 ) out3 << L[j][0][0];
               else
               {    out3 << "{";
                    for ( int k = 0; k < L[j].isize( ); k++ )
                    {    if ( k > 0 ) out3 << ",";
                         out3 << "{" << printSeq(L[j][k]) << "}";    }
                    out3 << "}";    }    }
          out3 << "\n";    }    }

void MakeTigs( const vec<vec<vec<int>>>& L, vec<vec<vec<vec<int>>>>& tigs )
{    tigs.clear( );
     int g = 0;
     for ( int i = 1; i < L.isize( ); i += 2 )
     {    if ( L[i].size( ) != 1 ) continue;
          if ( L[i][0].size( ) != 0 ) continue;
          vec<vec<vec<int>>> C;
          for ( int x = g; x < i; x++ )
               C.push_back( L[x] );
          tigs.push_back(C);
          g = i + 1;    }
     vec<vec<vec<int>>> C;
     for ( int x = g; x < L.isize( ); x++ )
          C.push_back( L[x] );
     if ( C.nonempty( ) ) tigs.push_back(C);    }

template void FindLines<HyperBasevector>( const HyperBasevector& dgraph, 
     const vec<int>& dinv, vec<vec<vec<vec<int>>>>& lines, 
     const int max_cell_paths, const int max_depth, const Bool verbose,
     const Bool single );
template void FindLines<digraphE<basevector>>( const digraphE<basevector>& dgraph, 
     const vec<int>& dinv, vec<vec<vec<vec<int>>>>& lines, 
     const int max_cell_paths, const int max_depth, const Bool verbose,
     const Bool single );
template void FindLines<digraphE<vec<int>>>( const digraphE<vec<int>>& dgraph, 
     const vec<int>& dinv, vec<vec<vec<vec<int>>>>& lines, 
     const int max_cell_paths, const int max_depth, const Bool verbose,
     const Bool single );
template void FindLines<digraphE<int>>( const digraphE<int>& dgraph, 
     const vec<int>& dinv, vec<vec<vec<vec<int>>>>& lines, 
     const int max_cell_paths, const int max_depth, const Bool verbose,
     const Bool single );

void Canonicalize( vec<vec<vec<vec<int>>>>& L )
{    for ( int j = 0; j < L.isize( ); j++ )
     {    vec<vec<vec<int>>>& x = L[j];
          if ( x.front( ) != x.back( ) || x.solo( ) ) continue;
          int m = x[0][0][0];
          int mi = 0;
          for ( int k = 2; k < x.isize( ); k += 2 )
          {    if ( x[k][0][0] < m )
               {    m = x[k][0][0];
                    mi = k;    }    }
          vec<int> mis;
          for ( int k = 0; k < x.isize( ); k += 2 )
               if ( x[k][0][0] == m ) mis.push_back(k);
          if ( mis.solo( ) && mis[0] == 0 ) continue;
          vec<vec<vec<vec<int>>>> y( mis.size( ) );
          for ( int i = 0; i < mis.isize( ); i++ )
          {    for ( int j = mi; j < x.isize( ); j++ )
                    y[i].push_back( x[j] );
               for ( int j = 1; j <= mi; j++ )
                    y[i].push_back( x[j] );    }
          x = Min(y);    }
     Sort(L);    }
