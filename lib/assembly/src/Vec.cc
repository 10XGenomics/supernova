///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include <algorithm>

#include <fcntl.h>
#include <fstream>
#include <math.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "String.h"
#include "system/ErrNo.h"
#include "system/Types.h"
#include "Vec.h"

void PrettyPrint(ostream& o, const vec<int>& v, int max_items, String terminator)
{    int chars_in_line = 0;
     for ( vec<int>::size_type i = 0; i < v.size( ); i++ )
     {    chars_in_line += int(floor(log10( static_cast<float>( v[i] )))) + 2;
          if ( (int) i == max_items - 1 )
          {    if ( chars_in_line + 3 > 80 ) o << "\n" << v[i];
               else o << v[i];
               o << " ...";
               break;    }
          if ( chars_in_line > 80 )
          {    o << "\n";
               chars_in_line = int(floor(log10( static_cast<float>( v[i] )))) + 2;    }
          o << v[i] << " ";    }
     o << terminator;    }

void PrettyPrint(ostream& o, const vec<longlong>& v, int max_items, 
     String terminator)
{    int chars_in_line = 0;
     for ( vec<longlong>::size_type i = 0; i < v.size( ); i++ )
     {    chars_in_line += int(floor(log10( static_cast<float>( v[i] )))) + 2;
          if ( (int) i == max_items - 1 )
          {    if ( chars_in_line + 3 > 80 ) o << "\n" << v[i];
               else o << v[i];
               o << " ...";
               break;    }
          if ( chars_in_line > 80 )
          {    o << "\n";
               chars_in_line = int(floor(log10( static_cast<float>( v[i] )))) + 2;    }
          o << v[i] << " ";    }
     o << terminator;    }

void PrettyPrint(ostream& o, const vec<double>& v, int max_items, String terminator)
{    int chars_in_line = 0;
     for ( vec<double>::size_type i = 0; i < v.size( ); i++ )
     {    String s = ToString( v[i] );
          chars_in_line += s.size( ) + 1;
          if ( (int) i == max_items - 1 )
          {    if ( chars_in_line + 3 > 80 ) o << "\n" << v[i];
               else o << v[i];
               o << " ...";
               break;    }
          if ( chars_in_line > 80 )
          {    o << "\n";
               chars_in_line = s.size( ) + 1;    }
          o << v[i] << " ";    }
     o << terminator;    }

void PrettyPrint(ostream& o, const vec<TraceInt>& v, int max_items, String terminator)
{    int chars_in_line = 0;
     for ( vec<TraceInt>::size_type i = 0; i < v.size( ); i++ )
     {    chars_in_line += int(floor(log10( static_cast<float>( v[i] )))) + 2;
          if ( (int) i == max_items - 1 )
          {    if ( chars_in_line + 3 > 80 ) o << "\n" << v[i];
               else o << v[i];
               o << " ...";
               break;    }
          if ( chars_in_line > 80 )
          {    o << "\n";
               chars_in_line = int(floor(log10( static_cast<float>( v[i] )))) + 2;    }
          o << v[i] << " ";    }
     o << terminator;    }


istream& operator>>(istream& s, vec<String>& v)
{    int n;
     s >> n;
     v.resize(n);
     char c;
     s.get(c);
     for ( vec<String>::size_type i = 0; i < v.size( ); i++ )
          getline( s, v[i] );
     return s;    }

// istream& operator>>(istream& s, vec<String>& v)
// {    int n;
//      s >> n;
//      v.resize(n);
//      char c;
//      s.get(c);
//      for ( vec<Sting>::size_type i = 0; i < v.size( ); i++ )
//           getline( s, v[i] );
//      return s;    }

ostream& operator<<(ostream& s, const vec<double>& v)
{    s << v.size( ) << "\n";
     for ( vec<double>::size_type i = 0; i < v.size( ); i++ )
          s << v[i] << "\n";
     return s;    }

ostream& operator<<(ostream& s, const vec<unsigned short>& v)
{    s << v.size( ) << "\n";
     for ( vec<unsigned short>::size_type i = 0; i < v.size( ); i++ )
          s << v[i] << "\n";
     return s;    }

ostream& operator<<(ostream& s, const vec<int>& v)
{    s << v.size( ) << "\n";
     for ( vec<int>::size_type i = 0; i < v.size( ); i++ )
          s << v[i] << "\n";
     return s;    }

ostream& operator<<(ostream& s, const vec<longlong>& v)
{    s << v.size( ) << "\n";
     for ( vec<longlong>::size_type i = 0; i < v.size( ); i++ )
          s << v[i] << "\n";
     return s;    }

ostream& operator<<(ostream& s, const vec<float>& v)
{    s << v.size( ) << "\n";
     for ( vec<float>::size_type i = 0; i < v.size( ); i++ )
          s << v[i] << "\n";
     return s;    }

ostream& operator<<(ostream& s, const vec<String>& v)
{    s << v.size( ) << "\n";
     for ( vec<String>::size_type i = 0; i < v.size( ); i++ )
          s << v[i] << "\n";
     return s;    }

template <>
void WriteAppend( const String& f, const vec<String>& v )
{
    ForceAssert( !IsRegularFile( f + ".gz" ) );
    if ( !IsRegularFile(f) )
    {
        std::ofstream out( f.c_str() );
        out << setfill('0') << setw(15) << v.size() << setfill(' ') << '\n';
        for ( size_t i = 0; i < v.size( ); i++ )
            out << v[i] << '\n';
        out.close();
    }
    else
    {
        std::fstream out(f.c_str(),std::ios_base::in|std::ios_base::out);
        size_t n;
        out >> n;
        size_t const max_size_bound = 10000000ul * 100000000ul;
        ForceAssertLt( n+v.size(), max_size_bound );
        out.seekp( 0, std::ios_base::beg );
        out << setfill('0') << setw(15) << n+v.size() << setfill(' ') << '\n';
        out.seekp( 0, std::ios_base::end );
        for ( size_t i = 0; i < v.size( ); i++ )
            out << v[i] << '\n';
        out.close();
    }
}

void PrintTabular( ostream& out, const vec< vec<String> >& rows, int sep,
     String justify)
{    int nrows = rows.size( ), ncols = 0;
     for ( int i = 0; i < nrows; i++ )
          ncols = std::max( ncols, (int) rows[i].size( ) );
     vec<int> maxcol;
     maxcol.resize_and_set( ncols, 0 );
     for ( int i = 0; i < (int) rows.size( ); i++ )
     {    for ( int j = 0; j < (int) rows[i].size( ); j++ )
               maxcol[j] = std::max( maxcol[j], (int) rows[i][j].size( ) );    }
     for ( int i = 0; i < (int) rows.size( ); i++ )
     {    for ( int j = 0; j < (int) rows[i].size( ); j++ )
          {    if ( j < (int) justify.size( ) && justify[j] == 'r' )
               {    for ( int k = 0; k < maxcol[j] - (int) rows[i][j].size( ); k++ )
                         out << " ";
                    out << rows[i][j];
                    if ( j < rows[i].isize( ) - 1 )
                    {    for ( int k = 0; k < sep; k++ )
                              out << " ";    }    }
               else
               {    out << rows[i][j];
                    if ( j < rows[i].isize( ) - 1 )
                    {    for ( int k = 0; 
                              k < maxcol[j] - (int) rows[i][j].size( ) + sep; k++ )
                              out << " ";    }    }    }
          out << "\n";    }    }

void PrintCSV(ostream& out, const vec< vec<String> >& rows)
{
    for (unsigned int i = 0; i < rows.size(); i++)
    {
        for (unsigned int j = 0; j < rows[i].size(); j++)
        {
            out << "\"" << rows[i][j] << "\"";
            if (j != rows[i].size()-1) { out << ","; }
            else { out << "\n"; }
        }
    }
}

bool IsAsciiVec( const String &filename )
{
  ifstream in( filename.c_str() );
  char c;
  while ( in )
  {
    in.get( c );
    if ( ! isdigit( c ) && ! isspace( c ) )
      return false;

    if ( c == '\n' ) 
      break;
  }
  return true;
}

longlong AsciiVecSize( const String& filename )
{
  String ns;
  Ifstream( in, filename );
  in >> ns;  
  ForceAssert( ns.IsInt( ) );
  longlong n = ns.Int( );
  return n;
}

#include "graph/DigraphTemplate.h"

template digraphV<vec<int> >::digraphV();
template void digraphV<vec<int> >::DeleteVertex(int);
template void digraphV<vec<int> >::DeleteVertices(const vec<int>&);
template void digraphV<vec<int> >::Initialize(const vec<vec<int> >&, const vec<vec<int> >&, const vec<vec<int> >&);
template const vec<int>& digraphV<vec<int> >::Vert(int) const;

template void digraphE< vec<int> >::RemoveUnneededVertices( );
template vec<int> digraphE< vec<int> >::RemoveDeadEdgeObjects( );
template void digraphE<vec<int> >::Initialize(vec<vec<int> > const&, vec<vec<int> > const&, vec<vec<int> > const&, vec<vec<int> > const&, vec<vec<int> > const&, 
const Bool);
template const vec<int>& digraphE<vec<int> >::EdgeObject(int) const;
template void digraphE<vec<int> >::PrettyDOT(std::basic_ostream<char, std::char_traits<char> >&, vec<double> const&, digraphE<vec<int> >::edge_label_info, unsigned char, unsigned char, vec<int> const*, vec<FeudalString<char, std::char_traits<char> > > const*, vec<int> const*, vec<unsigned char> const*, const vec<Bool>*,
     const vec<String>*, const vec<int>*, const String, const double,
     const double, const double, const double ) const;
template void digraphE<vec<int> >::ToLeft(vec<int>&) const;
template void digraphE<vec<int> >::ToRight(vec<int>&) const;
template void digraphE<vec<int> >::DeleteEdges(vec<int> const&, vec<int> const&);
