/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2014) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <set>

#include "CoreTools.h"
#include "graph/Digraph.h"
#include "graph/FindCells.h"

// Return false iff:
// 1) The cell contains any vertices, other than the cell's opening vertex, which
//    have predecessors outside the cell.
// 2) The cell contains any vertices, other than the cell's closing vertex, which
//    have successors outside the cell.
// Otherwise, return true.
// Assumes the cell is organized such that cell[0] is the opening vertex and
// cell.back() is the closing vertex.
bool
CellIsClosed( const vec<int> & cell, const digraph& G )
{
  
  for ( size_t i = 1; i < cell.size(); i++ ) // note that i=0 is skipped
    for ( size_t j = 0; j < G.To( cell[i] ).size(); j++ )
      if ( !Member( cell, G.To(cell[i])[j] ) )
	return false;
  
  for ( size_t i = 0; i+1 < cell.size(); i++ ) // note that i=cell.back() is skipped
    for ( size_t j = 0; j < G.From( cell[i] ).size(); j++ )
      if ( !Member( cell, G.From(cell[i])[j] ) )
	return false;
  
  return true;
}




// BuildCell: A helper function for FindCells, below.
// Attempt to find a cell going from v1 to vn, using {candidates} as a set of
// potential intermediates.
// If a cell is found, BuildCell returns a vec<int> cell, with cell.front() = v1
// and cell.back() = vn.
// If no cell is found, BuildCell returns an empty set.
//
// Note that the candidates list must be sorted.
vec<int>
BuildCell( const digraph & G, const int & v1, const int & vn,
	   const vec<int> & candidates, const size_t max_cell_size )
{
  vec<int> empty_cell(0);
  
  if ( candidates.solo() ) return empty_cell;
  
  
  // Step backward from vn and compile the set {vi} of vertices making up the
  // cell.
  vec<int> vis;
  vec<int> growth_spots( 1, vn );
  vec<int> next_growth_spots;
  
  for ( size_t i = 0; i < max_cell_size; i++ ) {
    
    // Add the current list of growth spots to the set {vi}.
    vis.append( growth_spots );
    
    // If {vi} (which does not yet include v1) is too big, give up.
    if ( vis.size() >= max_cell_size ) return empty_cell;
    
    
    // Step backward in the search, and find another iteration of vi's.
    next_growth_spots.clear();
    for ( size_t j = 0; j < growth_spots.size(); j++ ) {
      vec<int> tos = G.To( growth_spots[j] );
      for ( size_t k = 0; k < tos.size(); k++ ) {
	
	// Do not add v1 to the list {vi}.  This is how we reduce the list of
	// growth_spots to empty when we get back to the start of the cell.
	if ( tos[k] != v1 )
	  next_growth_spots.push_back( tos[k] );
      }
    }
    UniqueSort( next_growth_spots );
    
    // If any of these vi's have successors not in the candidates list, then the
    // cell is ruined.
    for ( size_t j = 0; j < next_growth_spots.size(); j++ ) {
      vec<int> froms = G.From( next_growth_spots[j] );
      for ( size_t k = 0; k < froms.size(); k++ ) {
	if ( BinPosition( candidates, froms[k] ) == -1 )
	  return empty_cell;
      }
    }
    
    growth_spots = next_growth_spots;
    if ( growth_spots.empty() ) break; // nothing left to do
  }
  
  
  // Sort the vis list, and make sure it begins with v1 and ends with vn.
  Sort( vis );
  vis.EraseValue( vn );
  vis.push_front(v1);
  vis.push_back(vn);
  
  if ( !CellIsClosed( vis, G ) ) return empty_cell;
  
  return vis;
}




// For documentation, see FindCells.h
void FindCells( const digraph& G, const size_t max_cell_size, vec< vec<int> >& cells )
{
  cells.clear();
  
  // Loop over every possible value of v1 - the starting vertex of the cell.
  for ( int v1 = 0; v1 < G.N( ); v1++ ) {
    
    // Ignore v1 if it is a self-loop.
    if ( BinPosition( G.To( v1 ), v1 ) != -1 ) continue;
    
    // Find candidate choices for vn - the endpoint of the cell.
    // This means doing a breadth-first search from v1, going as many as
    // max_cell_size steps.
    vec<int> growth_spots( 1, v1 );
    vec<int> candidates;
    for ( size_t i = 0; i < max_cell_size; i++ ) {
      
      // Step forward in the search.
      vec<int> next_growth_spots;
      for ( size_t j = 0; j < growth_spots.size(); j++ )
	next_growth_spots.append( G.From( growth_spots[j] ) );
      
      candidates.append( growth_spots );
      growth_spots = next_growth_spots;
    }
    UniqueSort( candidates );
    
    // Remove from the candidates list any vertices v which contain predecessors
    // not in the candidates list.  (By design this includes v1, which until now
    // is in the list.)  These vertices cannot possibly create cells.
    // Note that this excludes "hairs" that branch backwards.
    vec<Bool> erase( candidates.size(), False );
    for ( size_t i = 0; i < candidates.size(); i++ ) {
      vec<int> tos = G.To( candidates[i] );
      for ( size_t j = 0; j < tos.size(); j++ )
	if ( BinPosition( candidates, tos[j] ) == -1 ) {
	  //if ( G.From( tos[j] ).solo() && G.To( tos[j] ).empty() ) continue;
	  erase[i] = True;
	  break;
	}
    }
    EraseIf( candidates, erase );
    
    
    // Go through each candidate vn, and consider the possibility that it is
    // the tail end of a cell.  If we find a cell, add it to the list.
    for ( size_t i = 0; i < candidates.size( ); i++ ) {
      
      // Ignore vn if it is a self-loop.
      int vn = candidates[i];
      if ( BinPosition( G.To( vn ), vn ) != -1 ) continue;
      
      // Attempt to build a cell from v1 to vn.
      vec<int> cell = BuildCell( G, v1, vn, candidates, max_cell_size );
      if ( cell.empty() ) continue;
      
      // If the cell contains only two elements (v1 and vn), require them to
      // have multiple links between them.
      if ( cell.size() == 2 ) {
	vec<int> froms = G.From( cell[0] );
	if ( froms.CountValue( cell[1] ) < 2 ) continue;
      }
      
      // Note that v1 can only have one minimal cell coming out of it.  Hence if
      // we find a vn that completes a cell, there's no need to keep looking
      // at other vn's for this v1.
      cells.push_back( cell );
      break;
    }
  }
}

void FindSomeCells( const digraph& G, const int max_cell_size,
     const int max_depth, vec< pair<int,int> >& bounds )
{
     bounds.clear( );
     #pragma omp parallel for
     for ( int v = 0; v < G.N( ); v++ )
     {    
          // Consider only canonical cell entry vertices v.

          if ( !G.To(v).solo( ) || G.From(v).size( ) <= 1 ) continue;
          if ( Member( G.From(v), v ) ) continue;
          
          // Find vertices a bit downstream of the immediate successors of v.

          int no = G.From(v).size( );
          vec<vec<int>> down(no), downd(no);
          for ( int j = 0; j < no; j++ )
          {    down[j].push_back( G.From(v)[j] );
               downd[j].push_back(0);
               for ( int i = 0; i < down[j].isize( ); i++ )
               {    if ( downd[j][i] == max_depth ) break;
                    for ( int l = 0; l < G.From( down[j][i] ).isize( ); l++ )
                    {    int w = G.From( down[j][i] )[l], d = downd[j][i] + 1;
                         int p = Position( down[j], w );
                         if ( p < 0 || downd[j][p] > d )
                         {    down[j].push_back(w);
                              downd[j].push_back(d);    }    }    }
               UniqueSort( down[j] );    }

          // Find candidates for canonical cell exit vertices w.

          vec<int> ex;
          Intersection( down, ex );
          vec<Bool> to_del( ex.size( ), True );
          for ( int i = 0; i < ex.isize( ); i++ )
          {    int w = ex[i];
               if ( !G.From(w).solo( ) || G.To(w).size( ) <= 1 ) continue;
               if ( Member( G.To(w), w ) ) continue;
               to_del[i] = False;    }
          EraseIf( ex, to_del );

          // Test candidates.

          vec<int> ex2;
          vec<vec<int>> xs;
          for ( int i = 0; i < ex.isize( ); i++ )
          {    int w = ex[i];

               // Check for bounding of cell by v..w, and check cell size.

               vec<int> x = {v};
               Bool bad = False;
               for ( int j = 0; j < x.isize( ); j++ )
               {    if ( x.isize( ) > max_cell_size || G.From( x[j] ).empty( )
                         || G.To( x[j] ).empty( ) )
                    {    bad = True;
                         break;    }
                    if ( x[j] != w )
                    {    for ( int l = 0; l < G.From( x[j] ).isize( ); l++ )
                         {    int t = G.From( x[j] )[l];
                              if ( t == v )
                              {    bad = True;
                                   break;    }
                              if ( !Member( x, t ) ) x.push_back(t);    }    }
                    if ( x[j] != v )
                    {    for ( int l = 0; l < G.To( x[j] ).isize( ); l++ )
                         {    int t = G.To( x[j] )[l];
                              if ( t == w )
                              {    bad = True;
                                   break;    }
                              if ( !Member( x, t ) ) x.push_back(t);    }    }    }
               if ( bad || x.isize( ) > max_cell_size ) continue;

               // Check for cycles.

               for ( int j = 0; j < x.isize( ); j++ )
               {    if (bad) break;
                    if ( x[j] == w ) continue;
                    vec<int> m = { x[j] };
                    for ( int l = 0; l < m.isize( ); l++ )
                    {    if (bad) break;
                         for ( int r = 0; r < G.From( m[l] ).isize( ); r++ )
                         {    int z = G.From( m[l] )[r];
                              if ( z == x[j] )
                              {    bad = True;
                                   break;    }
                              if ( z == w ) continue;
                              if ( !Member( m, z ) ) m.push_back(z);    }    }    }
               if (bad) continue;
               xs.push_back(x);
               ex2.push_back(w);    }

          // Pick smallest.

          if ( ex2.empty( ) ) continue;
          vec<int> len( xs.size( ) ), ids( xs.size( ), vec<int>::IDENTITY );
          for ( int i = 0; i < xs.isize( ); i++ )
               len[i] = xs[i].size( );
          SortSync( len, ids );
          if ( ex2.size( ) >= 2 && len[0] == len[1] ) continue; // possible???
          int w = ex2[ ids[0] ];
          #pragma omp critical
          {    bounds.push( v, w );    }    }
     Sort(bounds);    }
