///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <map>
#include <queue>

#include "Bitvector.h"
#include "CoreTools.h"
#include "Equiv.h"
#include "Set.h"
#include "STLExtensions.h"
#include "graph/Digraph.h"
#include "math/Functions.h"

Bool digraph::HavePath( const int v, const int w ) const
{    vec<int> x, y;
     GetSuccessors1( v, x );
     GetPredecessors1( w, y );
     return Meet( x, y );    }

Bool digraph::IsSimpleLine( ) const
{    if ( N( ) == 1 && From(0).empty( ) && To(0).empty( ) ) return True;
     vec<Bool> used( N( ), False );
     for ( int v = 0; v < N( ); v++ )
     {    if ( To(v).nonempty( ) ) continue;
          used[v] = True;
          if ( !From(v).solo( ) ) return False;
          v = From(v)[0];
          used[v] = True;
          while(1)
          {    if ( !To(v).solo( ) ) return False;
               if ( From(v).empty( ) ) return N( ) == Sum(used);
               if ( !From(v).solo( ) ) return False;
               if ( v == From(v)[0] ) return False;
               v = From(v)[0];
               used[v] = True;    }
          return False;    }
     return False;    }

void digraph::FindSimpleLines( vec<vec<int>>& lines ) const
{    lines.clear( );
     vec<Bool> used( N( ), False );
     for ( int v = 0; v < N( ); v++ )
     {    vec<int> x = {v};
          if ( used[v] ) continue;
          used[v] = True;

          // Extend, first to the right, then to the left.

          int w = v;
          Bool circle = False;
          while ( From(w).size( ) == 1 )
          {    w = From(w)[0];
               if ( To(w).size( ) != 1 ) break;
               circle = Member( x, w );
               x.push_back(w);
               if (circle) break;
               used[w] = True;    }
          if ( !circle )
          {    w = v;
               while ( To(w).size( ) == 1 )
               {    w = To(w)[0];
                    if ( From(w).size( ) != 1 ) break;
                    x.push_front(w);
                    used[w] = True;    }    }
          lines.push_back(x);    }    }

void digraph::FindSimpleLines2( vec<vec<int>>& lines ) const
{    lines.clear( );
     vec<Bool> used( N( ), False );
     for ( int v = 0; v < N( ); v++ )
     {    vec<int> x = {v};
          if ( used[v] ) continue;

          // Extend, first to the right, then to the left.

          int w = v;
          Bool circle = False;
          while ( From(w).size( ) == 1 )
          {    w = From(w)[0];
               circle = Member( x, w );
               x.push_back(w);
               if ( To(w).size( ) != 1 ) break;
               if (circle) break;    }
          if ( !circle )
          {    w = v;
               while ( To(w).size( ) == 1 )
               {    w = To(w)[0];
                    x.push_front(w);
                    if ( From(w).size( ) != 1 ) break;    }    }
          if (circle)
          {    for ( int j = 0; j < x.isize( ); j++ )
                    used[ x[j] ] = True;    }
          if ( !circle && From( x[0] ).size( ) == 1 && x[0] != v ) continue;
          if ( !circle && From( x[0] ).size( ) > 1 && x[1] != v ) continue;
          lines.push_back(x);    }    }

void digraph::TransferEdges( int v, int w, const Bool enter_only )
{    ForceAssert( v != w );

     // Change edges v --> v to edges w --> w.

     if ( !enter_only )
     {
     vec<Bool> remove_from_v;
     remove_from_v.resize_and_set( from_[v].size( ), False );
     for ( int i = 0; i < from_[v].isize( ); i++ )
     {    if ( from_[v][i] == v )
          {    from_[w].push_back(w);
               to_[w].push_back(w);
               remove_from_v[i] = True;
               int j = InputFromOutputTo( v, i );
               to_[v].erase( to_[v].begin( ) + j );    }    }
     EraseIf( from_[v], remove_from_v );
     Sort( from_[w] ), Sort( to_[w] );
     }

     // Change edges u --> v to edges u --> w.
     
     for ( int i = 0; i < to_[v].isize( ); i++ )
     {    int u = to_[v][i];
          int j = InputToOutputFrom( v, i );
          from_[u][j] = w;
          Sort( from_[u] );    }

     // Change edges v --> x to edges w --> x.

     // if ( !enter_only )
     {    for ( int i = 0; i < from_[v].isize( ); i++ )
          {    int x = from_[v][i];
               int j = InputFromOutputTo( v, i );
               if ( !enter_only ) to_[x][j] = w;
               else to_[x][j] = v;
               Sort( to_[x] );    }    }

     // Do the rest.

     if ( !enter_only )
     {    from_[w].append( from_[v] );    }
     Sort( from_[w] );
     to_[w].append( to_[v] );
     Sort( to_[w] );
     to_[v].clear( );
     if ( !enter_only ) from_[v].clear( );    }

// CyclicCore successively deletes vertices and edges from the graph, without 
// actually deleting them, but tracking instead the number of edges entering and
// exiting each vertex.

void digraph::CyclicCore( vec<int>& core ) const
{    vec<int> in( N( ) ), out( N( ) ), sources, sinks;
     for ( int v = 0; v < N( ); v++ )
     {    in[v] = To(v).size( ), out[v] = From(v).size( );    }
     Sources(sources), Sinks(sinks);
     for ( int i = 0; i < sources.isize( ); i++ )
     {    int v = sources[i];
          out[v] = 0;
          for ( int j = 0; j < From(v).isize( ); j++ )
          {    int w = From(v)[j];
               in[w]--;
               if ( in[w] == 0 ) sources.push_back(w);    }    }
     for ( int i = 0; i < sinks.isize( ); i++ )
     {    int v = sinks[i];
          if ( in[v] == 0 ) continue;
          for ( int j = 0; j < To(v).isize( ); j++ )
          {    int w = To(v)[j];
               if ( in[w] == 0 ) continue;
               out[w]--;
               if ( out[w] == 0 ) sinks.push_back(w);    }    }
     core.clear( );
     for ( int v = 0; v < N( ); v++ )
          if ( in[v] > 0 && out[v] > 0 ) core.push_back(v);    }

Bool digraph::HasEdge( int64_t v, int64_t w ) const {
  return find( from_[v].begin(),  from_[v].end(), w ) != from_[v].end();
}

void digraph::Sources( vec<int>& v ) const
{    v.clear( );
     for ( int i = 0; i < N( ); i++ )
          if ( Source(i) ) v.push_back(i);    }

void digraph::Sinks( vec<int>& v ) const
{    v.clear( );
     for ( int i = 0; i < N( ); i++ )
          if ( Sink(i) ) v.push_back(i);    }

Bool digraph::LoopAt( const int v ) const 
{    set<int> A, B;
     for ( int i = 0; i < From(v).isize( ); i++ )
     {    int w = From(v)[i];
          if ( w == v ) return True;
          A.insert(w);    }
     while( A.size( ) > 0 )
     {    int x = *A.begin( );
          A.erase( A.begin( ) );
          B.insert(x);
          for ( int i = 0; i < From(x).isize( ); i++ )
          {    int y = From(x)[i];
               if ( y == v ) return True;
               if ( !Member( B, y ) ) A.insert(y);    }    }
     return False;    }

Bool digraph::HasCycle( const vec<int>& sub ) const
{    ForceAssert( sub.nonempty( ) );
     vec<int> s = sub;
     Sort(s);
     const int white = 0, red = 1, black = 2;
     for ( int i = 0; i < s.isize( ); i++ )
     {    vec<int> color( s.size( ), white );
          vec<int> reds;
          reds.push_back(i);
          color[i] = red;
          while( reds.nonempty( ) )
          {    int j = reds.back( );
               color[j] = black;
               reds.resize( reds.isize( ) - 1 );
               int v = s[j];
               for ( int t = 0; t < From(v).isize( ); t++ )
               {    int w = From(v)[t];
                    int p = BinPosition( s, w );
                    ForceAssertGe( p, 0 );
                    if ( p == i ) return True;
                    if ( color[p] == white ) 
                    {    color[p] = red;    
                         reds.push_back(p);    }    }    }    }
     return False;    }

void digraph::Initialize( const vec< vec<int> >& from, const vec< vec<int> >& to )
{    ForceAssertEq( from.size( ), to.size( ) );
     int N = from.size( );
     from_ = from, to_ = to;
     for ( int i = 0; i < N; i++ )
     {    for ( int j = 0; j < from[i].isize( ); j++ )
               CheckGoodVertex( from[i][j] );
          for ( int j = 0; j < to[i].isize( ); j++ )
               CheckGoodVertex( to[i][j] );
          ForceAssert( from[i].Ordered( ) );
          ForceAssert( to[i].Ordered( ) );    }
     for ( int i = 0; i < N; i++ )
     {    for ( int j = 0; j < from[i].isize( ); j++ )
               ForceAssert( BinMember( to[ from[i][j] ], i ) );
          for ( int j = 0; j < to[i].isize( ); j++ )
               ForceAssert( BinMember( from[ to[i][j] ], i ) );    }    }

digraph::digraph( const vec< vec<int> >& from, const vec< vec<int> >& to )
{    Initialize( from, to );    }

// Cribbed off digraphE::Initialize in DigraphTemplate.h - Josh Burton, May '08
void digraph::Initialize( const digraph& g, const vec<int>& v )
{   
  from_.resize( v.size( ) ), to_.resize( v.size( ) );
  vec<int> vsorted(v), vindex( v.size( ), vec<int>::IDENTITY );
  SortSync( vsorted, vindex );
  for ( int i = 0; i < v.isize( ); i++ )
    {    int x = v[i];
    for ( int j = 0; j < g.From(x).isize( ); j++ )
      {    int y = g.From(x)[j];
      int p2 = BinPosition( vsorted, y );
      if ( p2 < 0 ) continue;
      int i2 = vindex[p2];
      from_[i].push_back(i2);    
      to_[i2].push_back(i);
      }
    }
  for ( int i = 0; i < v.isize( ); i++ )
    { Sort( from_[i] );
      Sort( to_[i] ); }
}

digraph::digraph( const digraph& g, const vec<int>& v )
{ Initialize ( g, v ); }

digraph::digraph( const digraph& g, const int i )
{
  vec<vec<int> > components;
  g.Components( components );
  if ( i >= components.isize( ) || i < 0 )
    FatalErr("ERROR: input component ID (i = " << i << ") must be non-negative and smaller than the number of components (" << components.isize( ) << ")" );
  Initialize( g, components[i] );
}


digraph::digraph( const matrix< Bool >& A ) {
  Initialize( A );
}

void digraph::Initialize( const matrix< Bool >& A ) {
  ForceAssert( A.Square() );
  vec< vec< int > > newFrom( A.Ncols() ), newTo( A.Nrows() );
  for ( int v = 0; v < A.Ncols(); v++ ) {
    for ( int w = 0; w < A.Nrows(); w++ ) {
      if ( A( v, w ) ) {
	newFrom[ v ].push_back( w );
	newTo[ w ].push_back( v );
      }
    }
  }

  from_ = newFrom;
  to_ = newTo;
}


void PrintStandardDOTHeader( ostream& out )
{    out << "digraph G {\n\n";
     out << "rankdir=LR;\n";
     // out << "orientation=landscape;\n";
     out << "node [width=0.1,height=0.1,fontsize=10,shape=point,margin=\"0.05,0.045\"];\n";
     out << "edge [fontsize=12];\n";    }

void digraph::DOT( ostream& out ) const
{    PrintStandardDOTHeader(out);
     for ( int v = 0; v < N( ); v++ )
     {    for ( int j = 0; j < from_[v].isize( ); j++ )
               out << v << " -> " << from_[v][j] << ";\n";    }
     out << "\n}\n";    }

void digraph::DOT( ostream& out, const vec<String>& vertex_colors ) const
{    PrintStandardDOTHeader(out);
     for ( int v = 0; v < vertex_colors.isize( ); v++ )
     {    if ( vertex_colors[v] != "" )
          {    out << v << " [color=" << vertex_colors[v] << "];\n";    }    }
     for ( int v = 0; v < N( ); v++ )
     {    for ( int j = 0; j < from_[v].isize( ); j++ )
               out << v << " -> " << from_[v][j] << ";\n";    }
     out << "\n}\n";    }

void digraph::DOT( ostream& out, const vec<String>& vertex_colors,
     const vec< vec<String> >& edge_colors ) const
{    PrintStandardDOTHeader(out);
     for ( int v = 0; v < vertex_colors.isize( ); v++ )
     {    if ( vertex_colors[v] != "" )
          {    out << v << " [color=" << vertex_colors[v] << "];\n";    }    }
     for ( int v = 0; v < N( ); v++ )
     {    for ( int j = 0; j < from_[v].isize( ); j++ )
          {    const String& color = edge_colors[v][j];
               out << v << " -> " << from_[v][j];
               if ( color != "" ) out << " [color=" << color << "]";
               out << ";\n";    }    }
     out << "\n}\n";    }

void digraph::DOT( ostream& out, const vec< vec<String> >& edge_labels ) const
{    PrintStandardDOTHeader(out);
     for ( int v = 0; v < N( ); v++ )
     {    for ( int j = 0; j < from_[v].isize( ); j++ )
          {    int w = from_[v][j];
               const String& label = edge_labels[v][j];
               if ( label.size( ) > 0 )
               {    String ename = "E_" + ToString(v) + "_" + ToString(j);
                    out << "\n" << ename;
                    out << " [label=\"" << label << "\",shape=box]";
                    out << ";\n" << v << " -> " << ename << " [arrowhead=none];\n";
                    out << ename << " -> " << w << ";\n";    }
               else out << "\n" << v << " -> " << w << ";\n";    }    }
     out << "\n}\n";    }

void digraph::DOT( ostream& out, const vec< vec<String> >& edge_labels,
     const vec<String>& vertex_colors ) const
{    PrintStandardDOTHeader(out);
     for ( int v = 0; v < N( ); v++ )
     {    if ( vertex_colors[v].size( ) > 0 )
               out << v << "[color=" << vertex_colors[v] << "];\n";
          for ( int j = 0; j < from_[v].isize( ); j++ )
          {    int w = from_[v][j];
               const String& label = edge_labels[v][j];
               if ( label.size( ) > 0 )
               {    String ename = "E_" + ToString(v) + "_" + ToString(j);
                    out << "\n" << ename;
                    out << " [label=\"" << label << "\",shape=box]";
                    out << ";\n" << v << " -> " << ename << " [arrowhead=none];\n";
                    out << ename << " -> " << w << ";\n";    }
               else out << "\n" << v << " -> " << w << ";\n";    }    }
     out << "\n}\n";    }

void digraph::DOT( ostream& out, const vec< vec<String> >& edge_labels,
     const vec<String>& vertex_colors, const vec<vec<String>>& edge_colors, 
     const vec<vec<String>>& edge_attrs, const vec<vec<String>>& edge_attrs2,
     const vec<vec<String>>& legends, const vec<String>& legend_colors,
     const String& layout ) const
{    out << "digraph G {\n\n";
     out << "rankdir=LR;\n";
     out << "node [width=0.1,height=0.1,fontsize=10,shape=point,"
          << "margin=\"0.05,0.045\"];\n";
     out << "edge [fontsize=12];\n";
     if ( layout != "" ) out << "layout=" << layout << ";\n";
     for ( int v = 0; v < N( ); v++ )
     {    if ( vertex_colors[v].size( ) > 0 )
               out << v << "[color=" << vertex_colors[v] << "];\n";
          for ( int j = 0; j < from_[v].isize( ); j++ )
          {    int w = from_[v][j];
               const String& label = edge_labels[v][j];
               const String& color = edge_colors[v][j];
               const String& attr = edge_attrs[v][j];
               const String& attr2 = edge_attrs2[v][j];
               if ( label.size( ) > 0 )
               {    String ename = "E_" + ToString(v) + "_" + ToString(j);
                    out << "\n" << ename;
                    out << " [label=\"" << label << "\",shape=box";
                    if ( color != "" ) out << ",color=" << color;
                    if ( attr != "" ) out << "," << attr;
                    out << "]";
                    out << ";\n" << v << " -> " << ename << " [arrowhead=none";
                    if ( attr2 != "" ) out << "," << attr2;
                    if ( color != "" ) out << ",color=" << color;
                    out << "];\n";
                    out << ename << " -> " << w;
                    if ( color != "" || attr2 != "" )
                    {    out << "[";
                         if ( color != "" ) out << "color=" << color;
                         if ( color != "" && attr2 != "" ) out << ",";
                         if ( attr2 != "" ) out << attr2;
                         out << "]";    }
                    out << ";\n";    }
               else 
               {    out << "\n" << v << " -> " << w;
                    if ( color != "" || attr2 != "" )
                    {    out << "[";
                         if ( color != "" ) out << "color=" << color;
                         if ( color != "" && attr2 != "" ) out << ",";
                         if ( attr2 != "" ) out << attr2;
                         out << "]";    }
                    out << ";\n";    }    }    }
     PrintDotLegends( out, legends, legend_colors );
     out << "\n}\n";    }

void digraph::DOT_vl( ostream& out, const vec<String> & vertex_labels,
     const String& layout, const vec<String>& legend, const String& color,
     const String& vshape ) const
{    vec< vec<String> > legends;
     legends.push_back(legend);
     vec<String> colors;
     colors.push_back(color);
     DOT_vl( out, vertex_labels, layout, legends, colors, vec<vec<String>>( ),
          vec<String>( ), vshape );    }

String Quotex( const String& s )
{    return "\"" + s + "\"";    }

void digraph::DOT_el( ostream& out, const String& layout,
          const vec< vec<String> >& edge_labels) const
{    PrintStandardDOTHeader(out);
     if ( layout != "" ) out << "layout=" << layout << ";\n";
     for ( int v = 0; v < N( ); v++ )
     {    out << v << ";\n";
          for ( size_t j = 0; j < from_[v].size( ); j++ )
          {    out << v << " -> " << from_[v][j];
               if ( edge_labels.nonempty( ) )
                    out << " [label=" << Quotex(edge_labels[v][j]) << "]";
               out << ";\n";    }    }
     out << "\n}" << endl;    }

void PrintDotLegends( 
     ostream& out, const vec<vec<String>>& legends, const vec<String>& colors )
{    for ( int z = 0; z < legends.isize( ); z++ )
     {    const vec<String>& legend = legends[z];
          String quote = "\"";
          if ( legend.nonempty( ) )
          {    out << "node [shape=plaintext];\n";
               out << "struct" + ToString(z) + " [label=<\n";
               String color = ( colors.empty( ) ? "white"
                    : ( colors[z] == "" ? "white" : colors[z] ) );
               out << "<TABLE CELLSPACING=" << Quotex("0") << " BORDER=" 
                    << Quotex("0") << " CELLBORDER=" << Quotex("2") 
                    << " CELLPADDING=" << Quotex("0") << " BGCOLOR="
                    << Quotex(color) << ">\n";
               out << "<TR><TD>\n";
               out << "<TABLE CELLSPACING=" << Quotex("0") << " BORDER=" 
                    << Quotex("0") << ">\n";
               for ( int j = 0; j < legend.isize( ); j++ )
               {    out << "<TR><TD ALIGN=" << Quotex("left") << ">";
                    out << legend[j];
                    out << "</TD></TR>\n";    }
               out << "</TABLE>\n";
               out << "</TD></TR></TABLE>>\n";
               out << "];\n"; }    }    }

void digraph::DOT_vl( ostream& out, const vec<String> & vertex_labels,
     const String& layout, const vec< vec<String> > & legends, 
     const vec<String>& colors, const vec< vec<String> >& edge_labels,
     const vec<String>& vertex_colors, const String& vshape ) const
{    ForceAssertEq( N(), vertex_labels.isize() );
     PrintStandardDOTHeader(out);
     if ( layout != "" ) out << "layout=" << layout << ";\n";
     // Modify the standard DOT header to draw nodes as ellipses containing labels
     out << "node [shape=" << ( vshape == "" ? "ellipse" : vshape ) << "];\n";
     for ( int v = 0; v < N( );  v++ )
     {    if ( vertex_labels[v] != "" )
          {    out << v << " [label=\"" << vertex_labels[v] << "\"";
               if ( vertex_colors.nonempty( ) && vertex_colors[v] != "" )
               {    out << ",color=" << vertex_colors[v]
                        << ",fontcolor=" << vertex_colors[v];    }
               out  << "];\n";    }    }
     for ( int v = 0; v < N( ); v++ )
     {    for ( size_t j = 0; j < from_[v].size( ); j++ )
          {    out << v << " -> " << from_[v][j];
               if ( edge_labels.nonempty( ) )
                    out << " [label=" << Quotex(edge_labels[v][j]) << "]";
               out << ";\n";    }    }
     PrintDotLegends( out, legends, colors );
     out << "\n}" << endl;    }

template< > void digraphE<int>::DOT_vel( 
     ostream& out, const vec<String> & vertex_labels )
{    ForceAssertEq( N(), vertex_labels.isize() );
     PrintStandardDOTHeader(out);
     // Modify the standard DOT header to draw nodes as ellipses containing labels
     out << "node [shape=ellipse];\n";
     for ( int v = 0; v < N( );  v++ )
     {    if ( vertex_labels[v] != "" )
               out << v << " [label=\"" << vertex_labels[v] << "\"];\n";    }
     for ( int v = 0; v < N( ); v++ )
     {    for ( size_t j = 0; j < from_[v].size( ); j++ )
          {    out << v << " -> " << from_[v][j];
               out << " [label=" << Quotex(
                    ToString( EdgeObjectByIndexFrom( v, j ) ) ) << "]";
               out << ";\n";    }    }
     out << "\n}" << endl;    }

template< > void digraphE<int>::DOT0( ostream& out )
{    PrintStandardDOTHeader(out);
     for ( int v = 0; v < N( ); v++ )
     {    for ( size_t j = 0; j < from_[v].size( ); j++ )
          {    out << v << " -> " << from_[v][j];
               out << " [label=" << Quotex(
                    ToString( EdgeObjectByIndexFrom( v, j ) ) ) << "]";
               out << ";\n";    }    }
     out << "\n}" << endl;    }

template< > void digraphE<String>::DOT0( ostream& out )
{    PrintStandardDOTHeader(out);
     for ( int v = 0; v < N( ); v++ )
     {    for ( size_t j = 0; j < from_[v].size( ); j++ )
          {    out << v << " -> " << from_[v][j];
               out << " [label=" << Quotex( EdgeObjectByIndexFrom( v, j ) ) << "]";
               out << ";\n";    }    }
     out << "\n}" << endl;    }

void digraph::DOT_vl( ostream& out, const vec<String> & vertex_labels,
     vec<int> vs, const vec<String>& vertex_colors ) const
{    ForceAssertEq( N(), vertex_labels.isize() );
     Sort(vs);
     PrintStandardDOTHeader(out);
     // Modify the standard DOT header to draw nodes as ellipses containing labels
     out << "node [shape=ellipse];\n";
     for ( int v = 0; v < N();  v++ )
     {    if ( BinMember( vs, v ) )
          {    if ( vertex_labels[v] != "" )
                    out << v << " [label=\"" << vertex_labels[v] << "\"";
               if ( vertex_colors.nonempty( ) && vertex_colors[v] != "" )
                    out << ",color=\"" << vertex_colors[v] << "\"";
               out << "];\n";    }    }
     for ( int v = 0; v < N(); v++ )
     {    if ( BinMember( vs, v ) )
          {    for ( size_t j = 0; j < from_[v].size( ); j++ )
               {    int w = from_[v][j];
                    if ( BinMember( vs, w ) ) 
                         out << v << " -> " << w << ";\n";    }    }    }
     out << "\n}" << endl;    }

// Create a representation of a given graph in XML.
// http://graphml.graphdrawing.org
// The edge_labels arguments should be in bijective 
// correspondence with from_.
void digraph::WriteGraphML( ostream& out, const vec< vec<String> >& edge_labels ) const {
  out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
    "<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\"  \n"
   "  xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"  \n"
   "  xsi:schemaLocation=\"http://graphml.graphdrawing.org/xmlns \n"
   "  http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd\">  \n";
  out << "<graph id=\"G\" edgedefault=\"directed\">\n";

  out << "<key id=\"d0\" for=\"edge\"  attr.name=\"edgeName\" attr.type=\"string\" />\n";
  
  for ( int v = 0; v < N(); v++ )
    out << "  <node id=\"n" << v << "\" />\n";

  for ( int v = 0; v < N(); v++ ) {
    for ( int j = 0; j < from_[v].isize(); j++ ) {
      int w = from_[v][j];
      out << "<edge source=\"n" << v << "\" target=\"n" << w << "\" >\n";
      out << "  <data key=\"d0\">" << edge_labels[v][j] << "</data>\n";
      out << "</edge>\n";
    }
  }
  out << "</graph>\n";
  
  out << "</graphml>\n";
}

void digraph::ComponentRelation( equiv_rel& e ) const
{    e.Initialize( N( ) );
     for ( int v = 0; v < N( ); v++ )
     {    for ( int j = 0; j < From(v).isize( ); j++ )
          {    int w = From(v)[j];
               e.Join(v, w);    }    }    }

void digraph::Components( vec< vec<int> >& comp, const vec<Bool>* invisible ) const
{    comp.clear( );
     vec<Bool> used( N( ), False );
     if ( invisible != NULL ) used = *invisible;
     vec<int> C, Cnext;
     for ( int v = 0; v < N( ); v++ ) 
     {    if ( used[v] ) continue;
          C.clear( ), Cnext.clear( );
          Cnext.push_back(v);
          while( Cnext.nonempty( ) ) 
          {    int w = Cnext.back( );
               Cnext.pop_back( );
               if ( used[w] ) continue;
               used[w] = True;
               C.push_back(w);
               // Expand this component recursively
               Cnext.append( From(w) );
               Cnext.append( To(w) );    }
         Sort(C);
         comp.push_back(C);    }    }

void digraphX::Components( vec< vec<int> >& comp, const vec<Bool>* invisible ) const
{    comp.clear( );
     vec<Bool> used( N( ), False );
     if ( invisible != NULL ) used = *invisible;
     vec<int> C, Cnext;
     for ( int v = 0; v < N( ); v++ ) 
     {    if ( used[v] ) continue;
          C.clear( ), Cnext.clear( );
          Cnext.push_back(v);
          while( Cnext.nonempty( ) ) 
          {    int w = Cnext.back( );
               Cnext.pop_back( );
               if ( used[w] ) continue;
               used[w] = True;
               C.push_back(w);
               // Expand this component recursively
               Cnext.append( From(w) );
               Cnext.append( To(w) );    }
         Sort(C);
         comp.push_back(C);    }    }

void digraph::ComponentsAlt( vec< vec<int> >& comp ) const
{    comp.clear( );
     equiv_rel e( N( ) );
     for ( int v = 0; v < N( ); v++ )
     {    for ( int j = 0; j < From(v).isize( ); j++ )
               e.Join( v, From(v)[j] );    }
     for ( int x = 0; x < N( ); x++ )
     {    if ( e.Representative(x) )
          {    vec<int> o;
               e.Orbit( x, o );
               Sort(o);
               comp.push_back(o);    }    }    }


size_t digraph::NComponents() const
{
  vec< vec<int> > comps;
  this->Components( comps );
  return comps.size();
}

int digraph::ConnectedComponents( ) const
{    equiv_rel e;
     ComponentRelation(e);
     return e.OrbitCount( );    }

void digraph::CutPoints( vec<int>& cuts ) const
{    cuts.clear( );
     vec<int> reps;
     equiv_rel E;
     ComponentRelation(E);
     E.OrbitRepsAlt(reps);
     for ( int m = 0; m < reps.isize( ); m++ )
     {    vec<int> o;
          E.Orbit( reps[m], o );
          Sort(o);
          for ( int i = 0; i < o.isize( ); i++ )
          {    equiv_rel e( o.size( ) );
               for ( int v = 0; v < o.isize( ); v++ )
               {    if ( v == i ) continue;
                    for ( int j = 0; j < From( o[v] ).isize( ); j++ )
                    {    int w = From( o[v] )[j];
                         if ( w == o[i] ) continue;
                         e.Join( v, BinPosition( o, w ) );    }    }
               if ( e.OrbitCount( ) > 1 ) cuts.push_back( o[i] );    }    }
     Sort(cuts);    }

void digraph::GetPredecessors( const vec<int>& v, vec<int>& to_v ) const
{    set<int> check, tov;
     for ( int i = 0; i < v.isize( ); i++ )
     {    check.insert( v[i] );
          tov.insert( v[i] );    }
     while( !check.empty( ) )
     {    int x = *check.begin( );
          check.erase( check.begin( ) );
          for ( int i = 0; i < To(x).isize( ); i++ )
          {    int y = To(x)[i];
               if ( Member( tov, y ) ) continue;
               check.insert(y);
               tov.insert(y);    }    }
     to_v.clear( );
     for ( set<int>::iterator i = tov.begin( ); i != tov.end( ); ++i )
          to_v.push_back(*i);    }

void digraph::GetSuccessors( const vec<int>& v, vec<int>& from_v ) const
{    set<int> check, fromv;
     for ( int i = 0; i < v.isize( ); i++ )
     {    check.insert( v[i] );
          fromv.insert( v[i] );    }
     while( !check.empty( ) )
     {    int x = *check.begin( );
          check.erase( check.begin( ) );
          for ( int i = 0; i < From(x).isize( ); i++ )
          {    int y = From(x)[i];
               if ( Member( fromv, y ) ) continue;
               check.insert(y);
               fromv.insert(y);    }    }
     from_v.clear( );
     for ( set<int>::iterator i = fromv.begin( ); i != fromv.end( ); ++i )
          from_v.push_back(*i);    }

void digraph::GetPredecessors1( const int v, vec<int>& to_v ) const
{    vec<int> V;
     V.push_back(v);
     GetPredecessors( V, to_v );    }

void digraph::GetSuccessors1( const int v, vec<int>& from_v ) const
{    vec<int> V;
     V.push_back(v);
     GetSuccessors( V, from_v );    }

void digraph::SubgraphSources( const vec<int>& S, vec<int>& v ) const
{    ForceAssert( S.UniqueOrdered( ) );
     v.clear( );
     for ( int i = 0; i < S.isize( ); i++ )
     {    int x = S[i];
          Bool source = True;
          for ( int j = 0; j < To(x).isize( ); j++ )
          {    int y = To(x)[j];
               if ( BinMember( S, y ) )
               {    source = False;
                    break;    }    }
          if (source) v.push_back(x);    }    }

void digraph::SubgraphSinks( const vec<int>& S, vec<int>& v ) const
{    ForceAssert( S.UniqueOrdered( ) );
     v.clear( );
     for ( int i = 0; i < S.isize( ); i++ )
     {    int x = S[i];
          Bool sink = True;
          for ( int j = 0; j < From(x).isize( ); j++ )
          {    int y = From(x)[j];
               if ( BinMember( S, y ) )
               {    sink = False;
                    break;    }    }
          if (sink) v.push_back(x);    }    }

template<class E>
     void digraphE<E>::Distance( int v, int w, int max_dist, vec<int>& D ) const
{    set< pair<int,int> > unprocessed, processed;
     unprocessed.insert( make_pair( v, 0 ) );
     while( !unprocessed.empty( ) )
     {    int x = unprocessed.begin( )->first;
          int dx = unprocessed.begin( )->second;
          unprocessed.erase( unprocessed.begin( ) );
          processed.insert( make_pair( x, dx ) );
          for ( int j = 0; j < From(x).isize( ); j++ )
          {    int y = From(x)[j];
	        int dy = dx + EdgeObjectByIndexFrom( x, j );
               if ( dy > max_dist ) continue;
               if ( Member( processed, make_pair( y, dy ) ) ) continue;
               unprocessed.insert( make_pair( y, dy ) );    }    }
     D.clear( );
     for ( set< pair<int,int> >::iterator i = processed.begin( ); 
          i != processed.end( ); ++i )
     {    int x = i->first;
          int d = i->second;
          if ( x == w ) D.push_back(d);    }
     UniqueSort(D);    }

template 
     void digraphE<int>::Distance( int v, int w, int max_dist, vec<int>& D ) const;

Bool digraph::AllPaths( int v, int w, vec< vec<int> >& paths, int maxpaths,
     const Bool allow_self_loop, const int maxpushes ) const
{    paths.clear( );
     set< vec<int> > pathss;
     int pushes = 0;
     if ( v < 0 && w < 0 )
     {    vec<int> sinks;
          Sinks(sinks);
	  vec< vec<int> > paths0;
          for ( int j = 0; j < sinks.isize( ); j++ )
          {    Bool OK = AllPaths( -1, sinks[j], paths0, maxpaths );
               if ( !OK ) return False;
               paths.append(paths0);
               if ( maxpaths > 0 && paths.isize( ) > maxpaths ) return False;    }
          return True;    }
     set< pair< vec<int>, set<int> > > partials;
     vec<int> just;
     set<int> justs;

     // Case 1.

     if ( v >= 0 && w >= 0 )
     {
          // First find all the vertices up2 which are upstream of w.  But stop at v.

          set<int> up1, up2;
          up1.insert(w);
          while( !up1.empty( ) )
          {    int y = *up1.begin( );
               up1.erase( up1.begin( ) );
               up2.insert(y);
               if ( y == v && !allow_self_loop ) continue;
               for ( int i = 0; i < To(y).isize( ); i++ )
               {    int x = To(y)[i];
                    if ( !Member( up2, x ) ) up1.insert(x);    }    }
          if ( !Member( up2, v ) ) return False;

          // Now build the answer.

          just.push_back(v), justs.insert(v);
          partials.insert( make_pair( just, justs ) );
          while( !partials.empty( ) )
          {    pair< vec<int>, set<int> > p = *partials.begin( );
               partials.erase( partials.begin( ) );
               int x = p.first.back( );
               Bool to_process = True;
               if ( x == w ) 
               {    pathss.insert( p.first );
                    if ( maxpaths > 0 && (int) pathss.size( ) > maxpaths ) 
                         return False;
                    to_process = False;
                    if ( allow_self_loop && p.first.size( ) == 1 )
                         to_process = True;    }
               if (to_process)
               {    for ( int i = 0; i < From(x).isize( ); i++ )
                    {    int y = From(x)[i];
                         if ( allow_self_loop && v == w && y == v )
                         {    vec<int> pf = p.first;
                              pf.push_back(v);
                              pathss.insert(pf);
                              if ( maxpaths > 0 && (int) pathss.size( ) > maxpaths ) 
                                   return False;
                              continue;    }
                         if ( !Member( up2, y ) || Member( p.second, y ) ) continue;
                         pair< vec<int>, set<int> > q = p;
                         q.first.push_back(y);
                         q.second.insert(y);
                         partials.insert(q);    
                         if ( maxpushes >= 0 && ++pushes >= maxpushes ) return False;
                         if ( maxpaths > 0 && (int) partials.size( ) > maxpaths ) 
                              return False;    }    }    }    }

     // Case 2.

     if ( w < 0 )
     {    just.push_back(v), justs.insert(v);
          partials.insert( make_pair( just, justs ) );
          while( !partials.empty( ) )
          {    pair< vec<int>, set<int> > p = *partials.begin( );
               partials.erase( partials.begin( ) );
               int x = p.first.back( );
               if ( From(x).empty( ) )
               {    pathss.insert( p.first );
                    if ( maxpaths > 0 && (int) pathss.size( ) > maxpaths ) 
                         return False;    }
               else
               {    for ( int i = 0; i < From(x).isize( ); i++ )
                    {    int y = From(x)[i];
                         if ( Member( p.second, y ) ) continue;
                         pair< vec<int>, set<int> > q = p;
                         q.first.push_back(y);
                         q.second.insert(y);
                         partials.insert(q);
                         if ( maxpushes >= 0 && ++pushes >= maxpushes ) return False;
                         if ( maxpaths > 0 && (int) partials.size( ) > maxpaths ) 
                              return False;    }    }    }    }

     // Case 3.

     if ( v < 0 )
     {    just.push_back(w), justs.insert(w);
          partials.insert( make_pair( just, justs ) );
          while( !partials.empty( ) )
          {    pair< vec<int>, set<int> > p = *partials.begin( );
               partials.erase( partials.begin( ) );
               int y = p.first.front( );
               if ( To(y).empty( ) )
               {    pathss.insert( p.first );
                    if ( maxpaths > 0 && (int) pathss.size( ) > maxpaths ) 
                         return False;    }
               else
               {    for ( int i = 0; i < To(y).isize( ); i++ )
                    {    int x = To(y)[i];
                         if ( Member( p.second, x ) ) continue;
                         pair< vec<int>, set<int> > q = p;
                         q.first.push_front(x);
                         q.second.insert(x);
                         partials.insert(q);
                         if ( maxpushes >= 0 && ++pushes >= maxpushes ) return False;
                         if ( maxpaths > 0 && (int) partials.size( ) > maxpaths ) 
                              return False;    }    }    }    }

     paths.reserve( pathss.size( ) );
     for ( set< vec<int> >::iterator i = pathss.begin( ); i != pathss.end( ); ++i )
          paths.push_back(*i);
     return True;    }


template<class E> void digraphE<E>::AllPathsFixedLength( 
     int v, int w, int L, vec< vec<int> >& paths ) const
{    ForceAssertGe( L, 0 );
     paths.clear( );
     vec<int> to_right;
     ToRight(to_right);
     vec< vec<int> > partials;
     partials.push_back( vec<int>( ) );
     while( partials.nonempty( ) )
     {    vec<int> p = partials.back( );
          partials.pop_back( );
          int l = 0;
          for ( int i = 0; i < p.isize( ); i++ )
               l += EdgeObject( p[i] );
          int vn = ( p.empty( ) ? v : to_right[ p.back( ) ] );
          if ( l == L ) 
          {    if ( vn == w ) paths.push_back(p);    }
          else
          {    for ( int j = 0; j < From(vn).isize( ); j++ )
               {    int e = EdgeObjectIndexByIndexFrom( vn, j );
                    if ( l + EdgeObject(e) > L ) continue;
                    vec<int> p2 = p;
                    p2.push_back(e);
                    partials.push_back(p2);    }    }    }    }

template void digraphE<int>::AllPathsFixedLength( 
     int v, int w, int L, vec< vec<int> >& paths ) const;

template<class E> Bool digraphE<E>::AllPathsLengthRange( int v, int w, int L1, 
     int L2, const vec<int>& to_left, const vec<int>& to_right, 
     vec< vec<int> >& paths, int max_paths, int max_loops, const Bool no_dups, 
     const int max_copies ) const
{    paths.clear( );
     if ( L2 < L1 ) return True;

     // Find all edges that are to the right of v, and whose left end is at 
     // distance <= L2 from it.

     map<int,int> c1;
     vec<pair<int,int>> Q;
     for ( int i = 0; i < From(v).isize( ); i++ ) Q.push( IFrom(v,i), 0 );
     int loopcount = 0;
     while( Q.nonempty( ) )
     {    if ( max_loops > 0 && loopcount++ > max_loops ) return False;
          int e = Q.back( ).first, d = Q.back( ).second;
          Q.pop_back( );
          if ( c1.find(e) == c1.end( ) || d < c1[e] ) 
          {    c1[e] = d;
               int d2 = d + O(e);
               if ( d2 <= L2 )
               {    int x = to_right[e];
                    for ( int i = 0; i < From(x).isize( ); i++ )
                         Q.push( IFrom(x,i), d2 );    }    }    }

     // Now find all edges that are to the left of w, and whose right end is at
     // distance <= L2 from it.

     map<int,int> c2;
     Q.clear( );
     for ( int i = 0; i < To(w).isize( ); i++ ) Q.push( ITo(w,i), 0 );
     while( Q.nonempty( ) )
     {    if ( max_loops > 0 && loopcount++ > max_loops ) return False;
          int e = Q.back( ).first, d = Q.back( ).second;
          Q.pop_back( );
          if ( c2.find(e) == c2.end( ) || d < c2[e] ) 
          {    c2[e] = d;
               int d2 = d + O(e);
               if ( d2 <= L2 )
               {    int x = to_left[e];
                    for ( int i = 0; i < To(x).isize( ); i++ )
                         Q.push( ITo(x,i), d2 );    }    }    }

     // Find the edges in both sets.

     vec<int> x1, x2;
     for ( auto m = c1.begin( ); m != c1.end( ); m++ ) x1.push_back(m->first);
     for ( auto m = c2.begin( ); m != c2.end( ); m++ ) x2.push_back(m->first);
     Sort(x1), Sort(x2);
     vec<int> I = Intersection( x1, x2 );

     // Now start the main search.

     vec< vec<int> > partials;
     partials.push_back( vec<int>( ) );
     while( partials.nonempty( ) )
     {    if ( max_loops > 0 && loopcount++ > max_loops ) return False;
          vec<int> p = partials.back( );
          partials.pop_back( );
          int l = 0;
          for ( int i = 0; i < p.isize( ); i++ )
               l += EdgeObject( p[i] );
          int vn = ( p.empty( ) ? v : to_right[ p.back( ) ] );
          if ( l >= L1 && l <= L2 && vn == w ) 
          {    paths.push_back(p);
               if ( max_paths > 0 && paths.isize( ) > max_paths ) return False;    }
          if ( l < L2 )
          {    for ( int j = 0; j < From(vn).isize( ); j++ )
               {    int e = EdgeObjectIndexByIndexFrom( vn, j );
                    if ( l + EdgeObject(e) > L2 ) continue;
                    if ( no_dups && Member( p, e ) ) continue;
                    if ( !BinMember( I, e ) ) continue;
                    vec<int> p2 = p;
                    p2.push_back(e);
                    if ( max_copies > 0 )
                    {    vec<int> q(p2);
                         Sort(q);
                         Bool fail = False;
                         for ( int r = 0; r < q.isize( ); r++ )
                         {    int s = q.NextDiff(r);
                              if ( s - r > max_copies )
                              {    fail = True;
                                   break;    }
                              r = s - 1;    }
                         if (fail) continue;    }
                    partials.push_back(p2);    }    }    }
     return True;    }

template<class E> Bool digraphE<E>::AllPathsLengthRangeAlt( int v, int w, 
     int L1, int L2, const vec<int>& to_right, vec< vec<int> >& paths, 
     int max_paths, int max_loops, const Bool no_dups, const Bool eq_ok,
     const int max_partials ) const
{    ForceAssertGe( L1, 0 );
     paths.clear( );
     if ( L2 < L1 ) return True;
     vec< vec<int> > partials;
     partials.push_back( vec<int>( ) );
     int loopcount = 0;
     int first = 0;
     while( first < partials.isize( ) )
     {    if ( max_loops > 0 && loopcount++ > max_loops ) return False;
          vec<int> p = partials[first++];
          int l = 0;
          for ( int i = 0; i < p.isize( ); i++ )
               l += EdgeObject( p[i] );
          int vn = ( p.empty( ) ? v : to_right[ p.back( ) ] );
          if ( l >= L1 && l <= L2 && vn == w ) 
          {    paths.push_back(p);
               if ( eq_ok && paths.isize( ) == max_paths ) return True;
               if ( max_paths > 0 && paths.isize( ) > max_paths ) return False;    }
          if ( l < L2 )
          {    for ( int j = 0; j < From(vn).isize( ); j++ )
               {    int e = EdgeObjectIndexByIndexFrom( vn, j );
                    if ( l + EdgeObject(e) > L2 ) continue;
                    if ( no_dups && Member( p, e ) ) continue;
                    vec<int> p2 = p;
                    p2.push_back(e);
                    partials.push_back(p2);    
                    if ( max_partials > 0 && partials.isize( ) > max_partials )
                         return False;    }    }    }
     return True;    }

template Bool digraphE<int>::AllPathsLengthRange( int v, int w, int L1, 
     int L2, const vec<int>& to_left, const vec<int>& to_right, 
     vec< vec<int> >& paths, int max_paths,
     int max_loops, const Bool no_dups, const int ) const;

template Bool digraphE<int>::AllPathsLengthRangeAlt( int v, int w, int L1, 
     int L2, const vec<int>& to_right, vec< vec<int> >& paths, int max_paths,
     int max_loops, const Bool no_dups, const Bool eq_ok,
     const int max_partials ) const;

void digraph::DeleteEdgesAtVertex( int v )
{    for ( int i = 0; i < from_[v].isize( ); i++ )
     {    int w = from_[v][i];
          if ( v == w ) continue;
          for ( int j = to_[w].isize( ) - 1; j >= 0; j-- )
               if ( to_[w][j] == v ) to_[w].erase( to_[w].begin( ) + j );    }
     for ( int i = 0; i < to_[v].isize( ); i++ )
     {    int w = to_[v][i];
          if ( v == w ) continue;
          for ( int j = from_[w].isize( ) - 1; j >= 0; j-- )
               if ( from_[w][j] == v ) from_[w].erase( from_[w].begin( ) + j );    }
     from_[v].clear( ), to_[v].clear( );    }

/// Method: ReplaceWithTransitiveClosure()
/// Replace this graph with its transitive closure.  That is, if there is a
/// directed path from u to w in the original graph, there is a direct edge
/// (u,w) in the new graph.

void digraph::ReplaceWithTransitiveClosure() {
  matrix< Bool > A( N(), N() );
  for ( int v = 0; v < N(); v++ )
    for ( int w = 0; w < N(); w++ )
      A( v, w ) = HasEdge( v, w );

  int pathLen = 1;
  while ( pathLen < N() ) {
    matrix< Bool > Asq;
    mul( A, A, Asq );
    for ( int v = 0; v < N(); v++ )
      for ( int w = 0; w < N(); w++ )
	A( v, w ) = ( A( v, w ) || Asq( v, w ) );
    pathLen *= 2;
  }

  Initialize( A );
}

vec<int> digraph::ComponentOf( const int v ) const
{    vec< vec<int> > comp;
     Components(comp);
     for ( int i = 0; i < comp.isize( ); i++ )
          if ( BinMember( comp[i], v ) ) return comp[i];
     ForceAssert( 0 == 1 );
     return vec<int>( );    }

vec<int> digraph::VerticesConnectedTo( const vec<int>& v ) const
{    vec<int> x(v);
     UniqueSort(x);
     while(1)
     {    int n1 = x.size( );
          for ( int i = 0; i < n1; i++ )
          {    x.append( From( x[i] ) );
               x.append( To( x[i] ) );    }
          UniqueSort(x);
          if ( x.isize( ) == n1 ) break;    }
     return x;    }

void digraph::AddEdge( int64_t v, int64_t w )
{
  CheckGoodVertex( v );
  CheckGoodVertex( w );
       
  Assert( is_sorted( from_[v].begin(), from_[v].end() ) );
  int i = upper_bound( from_[v].begin(), from_[v].end(), w ) - from_[v].begin();
  from_[v].insert( from_[v].begin()+i, w );
  Assert( is_sorted( from_[v].begin(), from_[v].end() ) );

  Assert( is_sorted( to_[w].begin(), to_[w].end() ) );
  int j = upper_bound( to_[w].begin(), to_[w].end(), v ) - to_[w].begin();
  to_[w].insert( to_[w].begin()+j, v );
  Assert( is_sorted( to_[w].begin(), to_[w].end() ) );
  }

Bool digraph::Acyclic( ) const
{    vec< vec<int> > comp;
     Components(comp);
     Bool acyclic = True;
     for ( int i = 0; i < comp.isize( ); i++ )
          if ( HasCycle( comp[i] ) ) return False;
     return True;    }


// ribeiro 2009-09-16: 
// this is a non-recursive implementation of Tarjan's Strongly Connected Components algoritm 
// you can find a recursive pseudocode implementation of the algorithm on wikipedia

void digraph::StronglyConnectedComponents( vec< vec<int> >& SCCs ) const
{
  enum State { start, 
               run_iteration, 
               next_iteration,
               end, 
               finished };  // type defining the various states of the algorithm

  SCCs.clear(); // vector of Strongly Connected Components

  vec<int> Iv(N(), -1);  // index   of vertex v
  vec<int> LLv(N(), -1); // lowlink of vertex v

  vec<int> S;            // vec and set 
  set<int> vStk;           //  representing the vertex stack

  int iv = 0;              // the vertex ID 
  int iw = -1;             // the ID of a neighboring vertex to iv  
  int iv_neig;             // the index of a neighboring vertex of iv (0..n_neig-1) 

  vec<int> iv_stack;       // iv stack for mimicking a recursive call
  vec<int> iv_neig_stack;  // iv_neig stack for mimicking a recursive call

  State state = start;
  int I = 0;             // vertex index in order of visit

  while (state != finished) {  

    if (state == start) {

      if (iv != N()) { // still have vertices to visit
        Iv[iv] = LLv[iv] = I;
        I++;
        S.push_back(iv), vStk.insert(iv); // add iv to the stack
        
        iv_neig = 0;              // prepare iteration over 'to' neighbors
        state = run_iteration;
      }
      else {
        state = finished;
      }
    }


    if (state == run_iteration) {

      if (iv_neig != this->From(iv).isize()) { // still have 'to' neighbors to visit
        iw = this->From(iv)[iv_neig];  // vertex ID of neighbor
        
        if (Iv[iw] < 0) { // iw has not been visited yet
          iv_stack.push_back(iv);          // push to recursive stack
          iv_neig_stack.push_back(iv_neig);
          iv = iw;
          iv_neig = -1;
          iw = -1;

          state = start; // recursive call
        }
        else {
          if (Member(vStk, iw))  // vertex in the stack
            LLv[iv] = Min(LLv[iv], Iv[iw]);
          
          state = next_iteration;
        }
      }
      else {  // no more neighbors
        state = end;
      }
    }
    

    if (state == end) {

      if (LLv[iv] == Iv[iv]) { // vertex iv is the root of a SCC; remove it from the atack
        vec<int> SCC;
        int jv;
        do {
          jv = S.back();
          S.pop_back(), vStk.erase(jv);
          SCC.push_back(jv);
        } while (iv != jv);
        Sort(SCC);
        SCCs.push_back(SCC);  // add the SCC to the SCCs
      }

      if (iv_stack.size()) {   // if iv_stack is not empty this was a recursive call

        iw = iv;
        iv = iv_stack.back(); 
        iv_neig = iv_neig_stack.back();

        iv_stack.pop_back();
        iv_neig_stack.pop_back();
        
        LLv[iv] = Min(LLv[iv], LLv[iw]);

        state = next_iteration;
      }
      else {  // not a recursive call; go on to the next unvisited vertex

        while (iv < N() && Iv[iv] >= 0)
          iv++;
        
        state = start; 
      }
    }


    if (state == next_iteration) {  // next 'to' neighbor
      iv_neig ++;       
      state = run_iteration;
    }

  }
  Sort(SCCs);

}

template void digraphE<int>::DeleteEdgeFrom(int, int);
template void digraphE<int>::DeleteEdges(vec<int> const&);

template void digraphV<int>::DeleteEdgeFrom(int, int);
template void digraphV<int>::DeleteEdgeTo(int, int);

void digraph::PrettyDOT( ostream& out,
			 Bool label_contigs,
			 Bool label_vertices, 
			 const vec<int>* componentsToPrint, 
			 const vec<String> *label_contigs_extra,
			 const vec<int> *verticesToPrint ) const
{
  // Set up output.
  out << "digraph G {\n\n";
  if (label_vertices)
    out << "node [width=0.1,height=0.1,fontsize=12,shape=plaintext];\n";
  else out << "node [width=0.1,height=0.1,fontsize=10,shape=point];\n";
  out << "edge [fontsize=12];\n";
  if (label_contigs) out << "margin=1.0;\n";
  out << "rankdir=LR;\n";
  out << "labeljust=l;\n";
  
  // Define components.
  equiv_rel e;
  ComponentRelation(e);
  vec<int> reps;
  e.OrbitRepsAlt(reps);
  
  // Contig labels, and estimate space after label.
  vec<int> label_distance;
  vec<String> strContigLabel;
  if ( label_contigs ) {
    if ( label_contigs_extra )
      strContigLabel = *label_contigs_extra;
    else {
      strContigLabel.resize( reps.size( ) );
      for (size_t ii=0; ii<reps.size( ); ii++) {
	strContigLabel[ii] = "contig " + ToString( ii );
      }
    }
    
    label_distance.resize( strContigLabel.size( ), 0 );
    for (int ii=0; ii<(int)strContigLabel.size( ); ii++)
      label_distance[ii] = 1 + (int)( strContigLabel[ii].size( ) / 2 );
  }
  
  // Selected components.
  vec<int> select;
  if ( componentsToPrint ) select = *componentsToPrint;
  else {
    select.reserve( reps.size( ) );
    for (int ii=0; ii<(int)reps.size( ); ii++) select.push_back( ii );
  }
  
  // Vertices to be skipped.
  vec<bool> skip_vtx;
  if ( verticesToPrint ) {
    skip_vtx.resize( this->N( ), true );
    for (size_t ii=0; ii<verticesToPrint->size( ); ii++)
      skip_vtx[ (*verticesToPrint)[ii] ] = false;
  }

  // Print the contigs.  We put each contig in its own cluster (the
  // subgraph's name MUST start with "cluster" for this to have any effect).
  for ( int sel_id = select.isize( ) - 1; sel_id >= 0; sel_id-- ) {
    int i = select[sel_id];

    out << "\nsubgraph cluster" << i << " {\n";
    out << "color=white;\n";
    if ( label_contigs && label_contigs_extra )
      out << "label=\"" << strContigLabel[i]
	  << "\","
	  << "fontsize=18,"
	  << "fontname=\"Times-Bold\"\n";
    
    vec<int> o;
    e.Orbit( reps[i], o );
    
    // Find "leftmost" vertex.
    Sort(o);
    vec<float> pos( o.size( ) );
    vec<Bool> placed( o.size( ), False );
    pos[0] = 0.0, placed[0] = True;
    while( Sum(placed) < o.isize( ) ) {
      for ( int i1 = 0; i1 < o.isize( ); i1++ ) {
	int v = o[i1];
	for ( int j = 0; j < From(v).isize( ); j++ ) {
	  int w = From(v)[j];
	  int i2 = BinPosition( o, w );
	  if ( !( placed[i1] ^ placed[i2] ) ) continue;
	  if ( placed[i1] ) pos[i2] = pos[i1] + 1;
	  else pos[i1] = pos[i2] - 1;
	  placed[i1] = placed[i2] = True;
	}
      }
    }
    
    float left = Min(pos);
    int leftj = 0;
    for ( leftj = 0; leftj < pos.isize( ); leftj++ )
      if ( pos[leftj] == left ) break;
    int leftv = o[leftj];
    
    // Print component.
    for ( int vi = 0; vi < o.isize( ); vi++ ) {
      int v = o[vi];
      if ( verticesToPrint && skip_vtx[v] ) continue;

      if (label_vertices)
	out << v << " [label=" << "\"" << v << "\""  << ",fontcolor=black];\n";
      
      for ( int j = 0; j < From(v).isize( ); j++ ) {
	int w = From(v)[j];
	float wd = 2.0;
	String color, label;
	Bool bold = False;
        color = "black";
	out << v << " -> " << w
	    << " [minlen=" << wd << ",color=" << color;
	if (bold) out << ",style=bold";
	if ( label != "" ) out << ",label=\"" << label << "\"";
	if ( label_contigs && v == leftv && j == 0 && ! label_contigs_extra )
	  out << ",taillabel=\"" << strContigLabel[i] 
	      << "\",labelangle=180,"
	      << "weight=10000,"
	      << "labeldistance=" << label_distance[i] << ",labelfontsize=18,"
	      << "labelfontname=\"Times-Bold\"";
	out << "];\n";
      }
    }
    out << "}\n";
  }
  out << "\n}" << endl;
  
}

void digraph::readBinary( BinaryReader& reader )
{
    typedef vec<int>::const_iterator Itr;
    reader.read(&from_);
    to_.clear();
    to_.resize(from_.size());
    for ( size_t iii = 0; iii < from_.size(); ++iii )
    {
        vec<int> const& fff = from_[iii];
        for ( Itr itr(fff.begin()), end(fff.end()); itr != end; ++itr )
            to_[*itr].push_back(iii);
    }
}

void digraph::Reverse( )
{    for ( int i = 0; i < N( ); i++ )
          swap( from_[i], to_[i] );    }

Bool digraph::TestValid( const Bool exit ) const
{    if ( from_.size( ) != to_.size( ) )
          DIGRAPH_INVALID( "sizes of from_ and to_ are different", exit );
     for ( int v = 0; v < N( ); v++ )
     {    if ( !from_[v].Ordered( ) )
          {    DIGRAPH_INVALID( "from_[" << v << "] is not ordered", exit );    }
          if ( !to_[v].Ordered( ) )
          {    DIGRAPH_INVALID( "to_[" << v << "] is not ordered", exit );    }    }
     for ( int v = 0; v < N( ); v++ )
     {    if ( !to_[v].Ordered( ) )
               DIGRAPH_INVALID( "to_[" << v << "] is not sorted", exit );    }
     for ( int v = 0; v < N( ); v++ )
     {    for ( int j = 0; j < from_[v].isize( ); j++ )
          {    int w = from_[v][j];
               if ( w < 0 )
               {    DIGRAPH_INVALID( "There is an edge from " << v << " to \"vertex\" "
                         << w << ".", exit );    }    }    }
     for ( int v = 0; v < N( ); v++ )
     {    for ( int j = 0; j < to_[v].isize( ); j++ )
          {    int w = to_[v][j];
               if ( w < 0 )
               {    DIGRAPH_INVALID( "There is an edge from \"vertex\" " << w
                         << " to " << v << ".", exit );    }    }    }
     return True;    }

Bool digraph::TestValidParallel( ) const
{    if ( from_.size( ) != to_.size( ) )
          DIGRAPH_INVALIDP( "sizes of from_ and to_ are different" );
     #pragma omp parallel for
     for ( int v = 0; v < N( ); v++ )
     {    if ( !from_[v].Ordered( ) )
          {    DIGRAPH_INVALIDP( "from_[" << v << "] is not ordered" );    }
          if ( !to_[v].Ordered( ) )
          {    DIGRAPH_INVALIDP( "to_[" << v << "] is not ordered" );    }    }
     #pragma omp parallel for
     for ( int v = 0; v < N( ); v++ )
     {    if ( !to_[v].Ordered( ) )
               DIGRAPH_INVALIDP( "to_[" << v << "] is not sorted" );    }
     #pragma omp parallel for
     for ( int v = 0; v < N( ); v++ )
     {    for ( int j = 0; j < from_[v].isize( ); j++ )
          {    int w = from_[v][j];
               if ( w < 0 )
               {    DIGRAPH_INVALIDP( "There is an edge from " << v << " to \"vertex\" "
                         << w << "." );    }    }    }
     #pragma omp parallel for
     for ( int v = 0; v < N( ); v++ )
     {    for ( int j = 0; j < to_[v].isize( ); j++ )
          {    int w = to_[v][j];
               if ( w < 0 )
               {    DIGRAPH_INVALIDP( "There is an edge from \"vertex\" " << w
                         << " to " << v << "." );    }    }    }
     return True;    }

// Note on the implementation of ShortestPath.  Revision 43041 got rid of its
// quadratic behavior, but is known to have made it slower in some cases.

template<class EdgeT>
void digraphE<EdgeT>::ShortestPath(
     const int start, const int stop, vec<int>& path ) const
{    path.clear( );
     ForceAssertGe( start, 0 );
     ForceAssertGe( stop, 0 );
     ForceAssertLt( start, N( ) );
     ForceAssertLt( stop, N( ) );
     EdgeT const infinity = std::numeric_limits<EdgeT>::max();
     vec<EdgeT> dist( N( ), infinity );
     vec<int> prev( N( ), -1 );
     Bool negative = False;
     for ( int e = 0; e < EdgeObjectCount( ); e++ )
          if ( EdgeObject(e) < 0 ) negative = True;
     if ( !negative )
     {    dist[start] = 0;
          set< pair<EdgeT,int> > X;
          for ( int i = 0; i < N( ); i++ )
               X.insert( make_pair( dist[i], i ) );
          vec<Bool> Q( N( ), True );
          while( !X.empty( ) )
          {    pair<EdgeT,int> p = *X.begin( );
               EdgeT min_dist = p.first;
               int u = p.second;
               X.erase( X.begin( ) );
               Q[u] = False;
               if ( min_dist == infinity ) return;
               if ( u == stop ) break;
               for ( int j = 0; j < From(u).isize( ); j++ )
               {    int v = From(u)[j];
                    if ( Q[v] )
                    {    EdgeT alt = dist[u] + EdgeObjectByIndexFrom( u, j );
                         if ( alt < dist[v] )
                         {    X.erase( X.find( make_pair( dist[v], v ) ) );
                              dist[v] = alt;
                              X.insert( make_pair( dist[v], v ) );
                              prev[v] = u;    }    }    }    }    }
     else
     {    dist[start] = 0;
          for ( int i = 0; i < N( ); i++ )
          {    for ( int u = 0; u < N( ); u++ )
               {    for ( int j = 0; j < From(u).isize( ); j++ )
                    {    int v = From(u)[j];
                         if ( dist[u] + EdgeObjectByIndexFrom( u, j ) < dist[v] )
                         {    if ( i == N( ) - 1 )
                              {    cout << "ShortestPath: encountered negative "
                                        << "cycle." << endl;
                                   ForceAssert( 0 == 1 );    }
                              dist[v] = dist[u] + EdgeObjectByIndexFrom( u, j );
                              prev[v] = u;    }    }    }    }    }
     int u = stop;
     while(1)
     {    path.push_back(u);
          if ( prev[u] < 0 ) break;
          u = prev[u];    }
     path.ReverseMe( );    }

// Virtual ShortestPath.  Note optimization not inserted in non-virtual
// ShortestPath.  Also switch to priority_queue, not inserted.

template<class T> T digraphE_V1<T>::ShortestPath(
     const int start, const int stop, vec<int>& path ) const
{    path.clear( );
     ForceAssertGe( start, 0 );
     ForceAssertGe( stop, 0 );
     ForceAssertLt( start, N( ) );
     ForceAssertLt( stop, N( ) );
     T const infinity = std::numeric_limits<T>::max();
     vec<T> dist( N( ), infinity );
     vec<int> prev( N( ), -1 );
     dist[start] = 0;
     priority_queue< 
          pair<T,int>, vector< pair<T,int> >, std::greater< pair<T,int> > > X;
     for ( int i = 0; i < N( ); i++ )
          X.push( make_pair( dist[i], i ) );
     vec<Bool> Q( N( ), True );
     while( !X.empty( ) )
     {    const pair<T,int>& p = X.top( );
          T min_dist = p.first;
          int u = p.second;
          X.pop( );
          if ( min_dist > dist[u] ) continue;
          Q[u] = False;
          if ( min_dist == infinity ) return infinity;
          if ( u == stop ) break;
          vec< pair<int,T> > y = FromVec(u);
          for ( int j = 0; j < y.isize( ); j++ )
          {    int v = y[j].first;
               if ( Q[v] )
               {    ForceAssertGe( y[j].second, 0 );
                    T alt = dist[u] + y[j].second;
                    if ( alt < dist[v] )
                    {    dist[v] = alt;
                         X.emplace( dist[v], v );
                         prev[v] = u;    }    }    }    }
     int u = stop;
     while(1)
     {    path.push_back(u);
          if ( prev[u] < 0 ) break;
          u = prev[u];    }
     path.ReverseMe( );
     return dist[stop];    }

template double digraphE_V1<double>::ShortestPath(
     const int start, const int stop, vec<int>& path ) const;

template float digraphE_V1<float>::ShortestPath(
     const int start, const int stop, vec<int>& path ) const;

void digraph::DeleteEdgeTo( int w, int j )
{    int v = to_[w][j];
     int i = InputToOutputFrom( w, j );
     to_[w].erase( to_[w].begin( ) + j );
     from_[v].erase( from_[v].begin( ) + i );    }

void digraph::DeleteEdgeFrom( int v, int j )
{    int w = from_[v][j];
     int i = InputFromOutputTo( v, j );
     from_[v].erase( from_[v].begin( ) + j );
     to_[w].erase( to_[w].begin( ) + i );    }

int digraph::InputToOutputFrom( int w, int i ) const
{    int v = to_[w][i];
     for ( int j = 0; j < from_[v].isize( ); j++ )
          if ( from_[v][j] == w ) return j;
     ForceAssert( 0 == 1 );
     return -1;    }

int digraph::InputFromOutputTo( int w, int i ) const
{    int v = from_[w][i];
     for ( int j = 0; j < to_[v].isize( ); j++ )
          if ( to_[v][j] == w ) return j;
     ForceAssert( 0 == 1 );
     return -1;    }

// TEMPLATE INSTANTIATIONS

#include "graph/DigraphTemplate.h"

template digraphE<int>::digraphE();
template digraphE<int>::digraphE(digraph const&);
template digraphE<int>::digraphE(const vec<int>&, const ConstructorBehavior);
template digraphE<int>::digraphE(const digraphE<int>&, const equiv_rel&);
template digraphE<int>::digraphE(const vec<vec<int> >&, const vec<vec<int> >&, const vec<int>&, const vec<vec<int> >&, const vec<vec<int> >&, const Bool);
template int digraphE<int>::AddEdge(const int, const int, const int&);
template void digraphE<int>::Clear();
template const int& digraphE<int>::EdgeObject(int) const;
template int const& digraphE<int>::EdgeObjectByIndexFrom(int, int) const;
template int digraphE<int>::IFrom(int, int) const;
template int digraphE<int>::EdgeObjectCount() const;
template int digraphE<int>::EdgeObjectIndexByIndexFrom(int, int) const;
template int& digraphE<int>::EdgeObjectByIndexFromMutable(int, int);
template int const& digraphE<int>::EdgeObjectByIndexTo(int, int) const;
template int digraphE<int>::EdgeObjectIndexByIndexTo(int, int) const;
template vec<int> const& digraphE<int>::Edges() const;

template Bool digraphE<int>::EdgePaths(const int, const int, vec<vec<int> >&, 
     const int, const int, const int, const Bool ) const;

template vec<int> digraphE<int>::EdgesBetween(const int, const int) const;
template vec<vec<int> > const& digraphE<int>::FromEdgeObj() const;
template vec<int> const& digraphE<int>::FromEdgeObj(int) const;
template vec<int>& digraphE<int>::FromEdgeObjMutable(int);
template void digraphE<int>::Initialize(const vec<vec<int> >&, const vec<vec<int> >&, const vec<int>&, const vec<vec<int> >&, const vec<vec<int> >&, const Bool);
template void digraphE<int>::Initialize(digraph const&, vec<int> const&);
template void digraphE<int>::Initialize(digraphE<int> const&, equiv_rel_template<int> const&);
template int digraphE<int>::InputFromOutputTo(int, int) const;
template int digraphE<int>::InputToOutputFrom(int, int) const;
template int digraphE<int>::MaxEdge(int, int);
template int digraphE<int>::MinEdge(int, int);

template digraphE<double>::digraphE(const vec<vec<int> >&, const vec<vec<int> >&, const vec<double>&, const vec<vec<int> >&, const vec<vec<int> >&, const Bool);

template int digraphE< vec<int> >::EdgeObjectIndexByIndexTo(int, int) const;

template void digraphE<int>::PrettyDOT( ostream& out, 
     const vec<double>& lengths, const edge_label_info, Bool label_contigs, 
     Bool label_vertices, const vec<int>* componentsToPrint,
     const vec<String> *label_contigs_extra, const vec<int>* verticesToPrint,
     const vec<Bool>* dashed, const vec<Bool>* invisible,
     const vec<String>* edge_color, const vec<int>* pen_widths, const String,
     const double, const double, const double, const double ) const;

template vec<int> digraphE<int>::RemoveDeadEdgeObjects();
template void digraphE<int>::RemoveEdgelessVertices();
template void digraphE<int>::RemoveEdgelessVertices(vec<int> const&);
template void digraphE<int>::Reverse();
template void digraphE<int>::ShortestPath(int const,int const,vec<int>&) const;

template digraphE<int> digraphE<int>::Subgraph(vec<int> const& ) const;
template digraphE< vec<int> > digraphE< vec<int> >::Subgraph(vec<int> const& ) const;

template vec<vec<int> > const& digraphE<int>::ToEdgeObj() const;
template vec<int> const& digraphE<int>::ToEdgeObj(int) const;
template vec<int>& digraphE<int>::ToEdgeObjMutable(int);
template void digraphE<int>::ToLeft(vec<int>&) const;
template void digraphE<int>::ToRight(vec<int>&) const;
template void digraphE<int>::TransferEdges(int, int, unsigned char);
template void digraphE<int>::Used(vec<unsigned char>&) const;
template void digraphE<int>::readBinary(BinaryReader&);
template void digraphE<int>::writeBinary(BinaryWriter&) const;
template void digraphE<int>::JoinEdges(int,int const&);
template int& digraphE<int>::EdgeObjectMutable(int);

typedef std::pair<int,int> IntPair;
template digraphE<IntPair>::digraphE(const vec<vec<int> >&, const vec<vec<int> >&, const vec<IntPair>&, const vec<vec<int> >&, const vec<vec<int> >&, const Bool);

template void digraphE<IntPair>::DeleteEdgeFrom(int, int);
template void digraphE<IntPair>::DeleteEdges(const vec<int>&);
template IntPair const& digraphE<IntPair>::EdgeObjectByIndexFrom(int, int) const;
template int digraphE<IntPair>::EdgeObjectIndexByIndexFrom(int, int) const;
template vec<IntPair> digraphE<IntPair>::EdgeObjectsBetween(const int, const int) const;
template vec<int> digraphE<IntPair>::EdgesBetween(const int, const int) const;
template void digraphE<IntPair>::Initialize(vec<vec<int> > const&, vec<vec<int> > const&, vec<IntPair> const&, vec<vec<int> > const&, vec<vec<int> > const&, const Bool);
template int digraphE<IntPair>::InputFromOutputTo(int, int) const;

template digraphVE<int, int>::digraphVE();
template digraphVE<int,int>::digraphVE(const digraphE<int>&, const vec<int>&);
template int digraphVE<int, int>::N() const;
template const int& digraphVE<int,int>::Vert(int) const;
template int& digraphVE<int,int>::VertMutable(int);
template void digraphVE<int,int>::readBinary(BinaryReader&);
template void digraphVE<int,int>::writeBinary(BinaryWriter&) const;

template void digraphE<int>::ComponentsE( vec< vec<int> >& comp ) const;

template void digraphE<int>::AddVertices(int);

typedef pair<String,vec<int>> F2;
template digraphV<F2>::digraphV( const vec< vec<int> >&, 
     const vec< vec<int> >&, const vec<F2>& );
template const F2& digraphV<F2>::Vert( int v ) const;
template F2& digraphV<F2>::VertMutable( int );
template void digraphV<F2>::RemoveEdgelessVertices( );

template digraphE<int>::edge_label_info::edge_label_info(digraphE<int>::edge_label_info::ConstructorBehavior, unsigned char, unsigned char, vec<FeudalString<char, std::char_traits<char> > > const*);

template void 
     digraphE<int>::GetSuccessors(vec<int> const&, vec<std::pair<int, int> >&);

template int64_t digraphE<int>::CheckSum() const;

template
Bool digraphE<int>::EdgePaths( const vec<int>& to_left, const vec<int>& to_right,
     const int v, const int w, vec< vec<int> >& paths, const int max_copies,
     const int max_paths, const int max_iterations, const Bool ) const;

template void digraphE<double>::ShortestPath(int const,int const,vec<int>&) const;

/* stuff0 in src/10X/Stuff.cc */

// digraphE<vec<int>> instantiations

template void digraphE<vec<int>>::ToLeftParallel( vec<int>& to_left ) const;
template void digraphE<vec<int>>::ToRightParallel( vec<int>& to_right ) const;
template void digraphE<vec<int>>::readBinary(BinaryReader&);
template void digraphE<vec<int>>::writeBinary(BinaryWriter&) const;
template void digraphE<vec<int>>::TransferEdges(int, int, const Bool);
template void digraphE<vec<int>>::SplayVertex( const int v );
template void digraphE<vec<int>>::Used( vec<Bool>& used ) const;
template void digraphE<vec<int>>::Append( const digraphE<vec<int>>& );
template void digraphE<vec<int>>::Clear( );
template void digraphE<vec<int>>::ValidateLR( const vec<int>&, const vec<int>& ) 
     const;
template void digraphE<vec<int>>::ValidateLRParallel( const vec<int>&, 
     const vec<int>& ) const;
template void digraphE<vec<int>>::DeleteEdgesWithUpdate( 
     const vec<int>&, vec<int>&, vec<int>& );
template void digraphE<vec<int>>::TransferEdgesWithUpdate( int v, int w,
          vec<int>& to_left, vec<int>& to_right, const Bool enter_only );
template void digraphE<vec<int>>::AppendWithUpdate(
     const digraphE<vec<int>>& D, vec<int>& to_left, vec<int>& to_right );
template int digraphE<vec<int>>::AddEdgeWithUpdate( const int v, const int w, 
     const vec<int>& e, vec<int>& to_left, vec<int>& to_right );
template void digraphE<vec<int>>::Reverse( );
template bool operator==( 
     const digraphE<vec<int>>& g1, const digraphE<vec<int>>& g2 );
template bool operator!=( 
     const digraphE<vec<int>>& g1, const digraphE<vec<int>>& g2 );
template void digraphE<vec<int>>::ComponentsEFast( vec<vec<int>>& ) const;
template Bool digraphE< vec<int> >::EdgePaths( const int, const int, 
     vec< vec<int> >&, const int, const int, const int, const Bool ) const;
template Bool digraphE< vec<int> >::EdgePathsLim( const vec<int>&, const vec<int>&,
     const int, const int, const int, vec< vec<int> >&, const int, const int,
     const int ) const;
template vec<int>& digraphE<vec<int>>::EdgeObjectMutable( int i );
template void digraphE<vec<int>>::SplayVertexWithUpdate(
          const int v, vec<int>& to_left, vec<int>& to_right );
template void digraphE<int>::DeleteEdgesAtVertex(int);

template const pair<int,int>& digraphE<pair<int,int>>::EdgeObject(int) const;

template int digraphE<int>::IGoes(int, int) const;
template const int& digraphE<int>::OGoes(int, int) const;

template const int& digraphE<int>::OFrom(int, int) const;

// THIS DOESN'T BELONG HERE:

template void digraphE<BaseVec>::Initialize( 
     const digraphE<BaseVec>::ConstructorType2, 
     digraphEX<BaseVec> const&, vec<int> const& );
