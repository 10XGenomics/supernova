///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "feudal/BinaryStream.h"
#include "graph/Digraph.h"
#include "paths/HyperEfasta.h"

vec<String> SpecialEfastaTag2Graph( digraphE<basevector>& out, efasta in){
    vec<String> tagList;
    String extra;
    out.Clear();

    in.erase(std::remove(in.begin(), in.end(), ' '), in.end());
    in.erase(std::remove(in.begin(), in.end(), '\t'), in.end());
    in.erase(std::remove(in.begin(), in.end(), '\n'), in.end());

    if(in.size()<1) return tagList;
    auto nOpen  = std::count(in.begin(),in.end(),'{');
    auto nClose = std::count(in.begin(),in.end(),'}');
    auto nPOpen  = std::count(in.begin(),in.end(),'(');
    auto nPClose = std::count(in.begin(),in.end(),')');
    if( nClose!=nClose || nOpen!=1 || nClose!=1 || in.front()!='{' || in.back()!='}'
       || nPOpen != nPClose
      ){
        std::cout<<"WARNING: SpecialEfastaTag2Graph() unable to decode a tag "
                 << in.front() << "..." << in.back() << " with "
                 << nOpen << " {, " << nClose << "}, "
                 << nPOpen << " (, " << nPClose << "). "
                 << "Empty graph is returned." << std::endl;
        return tagList;
    }

    typedef std::tuple<String,int,int> entry_t;
    vec<entry_t> entries;

    const efasta::size_type end = in.size()-1;
    if(1==end||in[1]==','){
        entries.push_back(std::make_tuple("",-1,-1));
    }
    for(efasta::size_type begin=1 ; begin!=efasta::npos && begin<end; ){
        const char character = in[begin];
        if( character=='[' ){
            auto back = in.find_first_of(']',begin+1);
            if( back==efasta::npos || back>=end){
                std::cout<<"WARNING: SpecialEfastaTag2Graph() unable to decode a tag due to unmatched []. "
                         << "Empty graph is returned." << std::endl;
                tagList.clear();
                out.Clear();
                return tagList;
            }
            extra.append(in,begin,back+1-begin);
            begin=back+1;
        }
        else if( character=='('){
            auto back = in.find_first_of(')',begin+1);
            if( back==efasta::npos || back>=end){
                std::cout<<"WARNING: SpecialEfastaTag2Graph() unable to decode a tag due to unmatched (). "
                         << "Empty graph is returned." << std::endl;
                tagList.clear();
                out.Clear();
                return tagList;
            }
            auto nComma = std::count(in.begin()+begin,in.begin()+back+1,',');
            if( nComma != 2 ){
                std::cout<<"WARNING: SpecialEfastaTag2Graph() unable to decode a tag due to (...) does not contain two commas. "
                         << "Empty graph is returned." << std::endl;
                tagList.clear();
                out.Clear();
                return tagList;
            }
            String sLoc(in,begin+1,back-begin-1);
            auto src = sLoc.Before(",").Int();
            auto dst = sLoc.Between(",",",").Int();
            if( src<0 || dst <0){
                std::cout<<"WARNING: SpecialEfastaTag2Graph() unable to decode a tag due to negative vertex indices " << src << " " << dst << "."
                         << "Empty graph is returned." << std::endl;
                tagList.clear();
                out.Clear();
                return tagList;
            }
            entries.push_back(std::make_tuple(sLoc.SafeAfterLast(",") ,src,dst));
            begin=back+1;
        }
        else if(character==','){
            if( begin+1==end || in[begin+1]==','){
                entries.push_back(std::make_tuple("",-1,-1));
            }
            ++begin;
        }
        else if(  character=='A'||character=='C'||character=='G'||character=='T'
                ||character=='a'||character=='c'||character=='g'||character=='t'){
            auto back_plus_one = in.find_first_not_of("ACGTacgt",begin+1);
            if( back_plus_one==efasta::npos || back_plus_one>=end){
                entries.push_back(std::make_tuple(in.substr(begin,end-begin),-1,-1));
                begin=end;
            }
            else{
                entries.push_back(std::make_tuple(in.substr(begin,back_plus_one-begin),-1,-1));
                begin=back_plus_one;
            }
        }
        else{
            std::cout<<"WARNING: SpecialEfastaTag2Graph() unable to decode a tag due the presence of character " << character << "."
                     << "Empty graph is returned." << std::endl;
            tagList.clear();
            out.Clear();
            return tagList;
        }
    }
    int maxIndex=-1;
    for( const auto& entry:entries){
        maxIndex = std::max(maxIndex,std::get<1>(entry));
        maxIndex = std::max(maxIndex,std::get<2>(entry));
    }
    int nVertices=maxIndex+1;
    for( auto& entry:entries){
        if(std::get<1>(entry)<0){
            std::get<1>(entry)=nVertices++;
        }
        if(std::get<2>(entry)<0){
            std::get<2>(entry)=nVertices++;
        }
    }
    out.Clear();
    out.AddVertices(nVertices);
    for( const auto& entry:entries){
        out.AddEdge(std::get<1>(entry),std::get<2>(entry),basevector(std::get<0>(entry)));
    }
    if( extra.size()>0){
        tagList.push_back( "HyperEfastaTag=\""+extra+"\"");
    }
    return tagList;
}

HyperEfasta::HyperEfasta( const HyperBasevector& h  )
{    vec<efasta> edges( h.EdgeObjectCount( ) );
     for ( int i = 0; i < edges.isize( ); i++ )
          edges[i] = h.EdgeObject(i).ToString( );
     (digraphE<efasta>&)(*this) = digraphE<efasta>( h.From( ), h.To( ), edges,
          h.ToEdgeObj( ), h.FromEdgeObj( ), True );
     K_ = h.K( );    }

void HyperEfasta::Initialize( const HyperBasevector& h )
{    vec<efasta> edges( h.EdgeObjectCount( ) );
     for ( int i = 0; i < edges.isize( ); i++ )
          edges[i] = h.EdgeObject(i).ToString( );
     (digraphE<efasta>&)(*this) = digraphE<efasta>( h.From( ), h.To( ), edges,
          h.ToEdgeObj( ), h.FromEdgeObj( ) );
     K_ = h.K( );    }

void HyperEfasta::writeBinary( BinaryWriter& writer ) const
{    writer.write(K_);
     writer.write(static_cast<digraphE<efasta> const&>(*this));    }

void HyperEfasta::readBinary( BinaryReader& reader )
{    reader.read(&K_);
     reader.read(static_cast<digraphE<efasta>*>(this));    }

void HyperEfasta::PrintSummaryDOT0w( ostream& out, Bool label_contigs,
     Bool label_vertices, Bool label_edges, const vec<int>* componentsToPrint,
     const Bool edge_labels_base_alpha, const vec<String> *label_edges_extra,
     const vec<String> *label_contigs_extra, const vec<int>* verticesToPrint ) const
{    vec<double> lengths( EdgeObjectCount( ) );
     for ( int i = 0; i < EdgeObjectCount( ); i++ )
          lengths[i] = EdgeLengthKmers(i);
     PrettyDOT( out, lengths, edge_label_info( 
          edge_label_info::EXPLICIT, label_edges,
          edge_labels_base_alpha, label_edges_extra ), 
          label_contigs, label_vertices, 
          componentsToPrint, label_contigs_extra, verticesToPrint );    }

void HyperEfasta::RemoveUnneededVertices( )
{    for ( int i = 0; i < N( ); i++ )
     {    if ( From(i).size( ) == 1 && To(i).size( ) == 1 && From(i)[0] != i )
          {    efasta p = EdgeObjectByIndexTo( i, 0 );
               efasta q = EdgeObjectByIndexFrom( i, 0 );
               ForceAssertGe( p.isize( ), K( ) - 1 );
               ForceAssertGe( q.isize( ), K( ) - 1 );
               Bool right_ok = True;
               for ( int j = 0; j < K( ) - 1; j++ )
               {    if ( q[j] == '{' )
                    {    right_ok = False;
                         break;    }    }
               if (right_ok) q = String( q, K( ) - 1, q.isize( ) - ( K( ) - 1 ) );
               else
               {    for ( int j = 0; j < K( ) - 1; j++ )
                    {    if ( p[ p.isize( ) - 1 - j ] == '}' )
                              cout << "\np = " << p << "\n\n";
                         ForceAssert( p[ p.isize( ) - 1 - j ] != '}' );    }
                    p.resize( p.size( ) - K( ) + 1 );    }
               p.append(q);
               JoinEdges( i, p );    }    }
     RemoveEdgelessVertices( );    }

#include "graph/DigraphTemplate.h"

template const vec< vec<int> >& digraphE<efasta>::FromEdgeObj() const;
template const vec< vec<int> >& digraphE<efasta>::ToEdgeObj() const;
template digraphE<efasta>::digraphE();
template digraphE<efasta>::digraphE(const vec<vec<int> >&, const vec<vec<int> >&, const vec<efasta>&, const vec<vec<int> >&, const vec<vec<int> >&, const Bool );
template int digraphE<efasta>::AddEdge(const int, const int, const efasta&);
template void digraphE<efasta>::DeleteEdgeTo(int, int);
template void digraphE<efasta>::DeleteEdgeFrom(int, int);
template const efasta& digraphE<efasta>::EdgeObject(int) const;
template efasta& digraphE<efasta>::EdgeObjectMutable(int) ;
template efasta const& digraphE<efasta>::EdgeObjectByIndexFrom(int, int) const;
template efasta const& digraphE<efasta>::EdgeObjectByIndexTo(int, int) const;
template int digraphE<efasta>::EdgeObjectCount() const;
template int digraphE<efasta>::EdgeObjectIndexByIndexFrom(int, int) const;
template void digraphE<efasta>::Initialize(vec<vec<int> > const&, vec<vec<int> > const&, vec<efasta> const&, vec<vec<int> > const&, vec<vec<int> > const&, const Bool);
template int digraphE<efasta>::InputFromOutputTo(int, int) const;
template int digraphE<efasta>::InputToOutputFrom(int, int) const;
template void digraphE<efasta>::JoinEdges(int, const efasta&);

template void digraphE<efasta>::PrettyDOT( ostream& out, 
     const vec<double>& lengths, const edge_label_info, Bool label_contigs, 
     Bool label_vertices, const vec<int>* componentsToPrint,
     const vec<String> *label_contigs_extra, const vec<int>* verticesToPrint,
     const vec<Bool>* dashed, const vec<Bool>* invisible,
     const vec<String>* edge_color, const vec<int>* pen_widths, const String,
     const double, const double, const double, const double ) const;

template void digraphE<efasta>::readBinary(BinaryReader&);
template vec<int> digraphE<efasta>::RemoveDeadEdgeObjects();
template void digraphE<efasta>::RemoveEdgelessVertices();
template void digraphE<efasta>::RemoveEdgelessVertices(vec<int> const&);
template void digraphE<efasta>::Used(vec<unsigned char>&) const;
template void digraphE<efasta>::writeBinary(BinaryWriter&) const;
template void digraphE<efasta>::ToLeft(vec<int>&) const;
template void digraphE<efasta>::ToRight(vec<int>&) const;
template int digraphE<efasta>::EdgeObjectIndexByIndexTo(int, int) const;
template void digraphE<efasta>::DeleteEdges(vec<int> const&);

template const vec<int>& digraphE<efasta>::FromEdgeObj(int) const;
template const vec<int>& digraphE<efasta>::ToEdgeObj(int) const;

template digraphE<efasta>::edge_label_info::edge_label_info(digraphE<efasta>::edge_label_info::ConstructorBehavior, unsigned char, unsigned char, vec<FeudalString<char, std::char_traits<char> > > const*);

template digraphE<efasta>::edge_label_info::edge_label_info(digraphE<efasta>::edge_label_info::ConstructorBehavior, vec<FeudalString<char, std::char_traits<char> > > const*);

template efasta& digraphE<efasta>::EdgeObjectByIndexFromMutable(int, int);
template const vec<efasta>& digraphE<efasta>::Edges() const;

template Bool digraphE<efasta>::EdgePaths(int, int, vec<vec<int> >&, int, int, int, const Bool) 
const;

template Bool digraphE<efasta>::EdgePaths( const vec<int>& left, 
     const vec<int>& right, const int v, const int w, vec< vec<int> >& paths, 
     const int max_copies, const int max_paths, const int max_iterations, const Bool ) const;
