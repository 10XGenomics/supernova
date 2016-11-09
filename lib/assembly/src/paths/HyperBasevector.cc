///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <strstream>

#include "CoreTools.h"
#include "Equiv.h"
#include "ParallelVecUtilities.h"
#include "Qualvector.h"
#include "feudal/BinaryStream.h"
#include "feudal/IncrementalWriter.h"
#include "graph/Digraph.h"
#include "paths/HyperBasevector.h"
#include "paths/KmerBaseBroker.h"

HyperBasevector::HyperBasevector( const HyperBasevectorX& hbx ) :
    digraphE<basevector>( hbx.AsDigraphE() ) {
    K_ = hbx.K( );
}

void HyperBasevector::SetToDisjointUnionOf( const vec<HyperBasevector>& v, const bool clear /* = true */ )
{
     for ( size_t i = 0; i < v.size(); i++ ) {
          ForceAssertEq( K( ), v[i].K( ) );
          ForceAssertNe( &v[i], this );			// probably not a good idea to pass one's self in
     }

     if ( clear ) Clear();				// should clear N() and EdgeObjectCount() to zero

     int nvert = N(), nedge = EdgeObjectCount();	// .reserve() for total number of (v,e) that we'll have
     for ( size_t i = 0; i < v.size(); i++ )
     {    nvert += v[i].N( );
          nedge += v[i].EdgeObjectCount( );    }
     FromMutable( ).reserve(nvert), ToMutable( ).reserve(nvert);
     FromEdgeObjMutable( ).reserve(nvert), ToEdgeObjMutable( ).reserve(nvert);
     EdgesMutable( ).reserve(nedge);

     int vcount = N(), ecount = EdgeObjectCount();	// array offset for the next (v,e) to be appended
     for ( size_t i = 0; i < v.size(); i++ )
     {    const HyperBasevector& h = v[i];
          EdgesMutable( ).append( h.Edges( ) );
          vec< vec<int> > from = h.From( ), to = h.To( );
          vec< vec<int> > frome = h.FromEdgeObj( ), toe = h.ToEdgeObj( );
          for ( int j = 0; j < h.N( ); j++ )
          {    for ( size_t u = 0; u < from[j].size(); u++ )
               {    from[j][u] += vcount;
                    frome[j][u] += ecount;    }
               for ( size_t u = 0; u < to[j].size(); u++ )
               {    to[j][u] += vcount;
                    toe[j][u] += ecount;    }    }
          FromMutable( ).append(from), ToMutable( ).append(to);
          FromEdgeObjMutable( ).append(frome), ToEdgeObjMutable( ).append(toe);
          vcount += h.N( );
          ecount += h.EdgeObjectCount( );    }    }


void HyperBasevector::Reverse( )
{    for ( int i = 0; i < EdgeObjectCount( ); i++ )
          EdgeObjectMutable(i).ReverseComplement( );
     digraphE<basevector>::Reverse( );    }

void HyperBasevector::PrintSummaryDOT0w( ostream& out, Bool label_contigs,
     Bool label_vertices, Bool label_edges, const vec<int>* componentsToPrint,
     const Bool edge_labels_base_alpha, const vec<String> *label_edges_extra,
     const vec<String> *label_contigs_extra, const vec<int>* verticesToPrint,
     const vec<String> *edge_color) const
{    vec<double> lengths( EdgeObjectCount( ) );
     for ( int i = 0; i < EdgeObjectCount( ); i++ )
         lengths[i] = EdgeLengthKmers(i);
     PrettyDOT( out, lengths, digraphE<basevector>::edge_label_info(
          digraphE<basevector>::edge_label_info::EXPLICIT, label_edges,
          edge_labels_base_alpha, label_edges_extra ),
          label_contigs, label_vertices,
          componentsToPrint, label_contigs_extra, verticesToPrint, NULL,
          NULL, edge_color, NULL, "", K( ) + K( )/2 );    }

template <class KBB>
vec<basevector> Convert( const vec<KmerPath>& in, KBB const& kbb )
{    vec<basevector> out;
     for ( size_t j = 0; j < in.size(); j++ )
     {    SuperBaseVector s = kbb.ToSequence( in[j] );
          ForceAssertEq( s.size( ), 1 );
          out.push_back( s.Seq(0) );    }
     return out;    }

HyperBasevector::HyperBasevector( const String& filename )
{    BinaryReader::readFile( filename, this ); }

HyperBasevector::HyperBasevector( const HyperKmerPath& h, const KmerBaseBroker& kbb )
     : digraphE<basevector>( h.From( ), h.To( ), Convert( h.Edges( ), kbb ),
                             h.ToEdgeObj( ), h.FromEdgeObj( ) ),
       K_( h.K() )
{
  TestValid();
}

HyperBasevector::HyperBasevector( const HyperKmerPath& h, const KmerBaseBrokerBig& kbb )
     : digraphE<basevector>( h.From( ), h.To( ), Convert( h.Edges( ), kbb ),
                             h.ToEdgeObj( ), h.FromEdgeObj( ) ),
       K_( h.K() )
{
  TestValid();
}

void HyperBasevector::Initialize( const HyperKmerPath& h, const KmerBaseBroker& kbb )
{    (digraphE<basevector>&)(*this) =
          digraphE<basevector>( h.From( ), h.To( ), Convert( h.Edges( ), kbb ),
          h.ToEdgeObj( ), h.FromEdgeObj( ) );
     K_ = h.K( );
     TestValid();
}

void HyperBasevector::writeBinary( BinaryWriter& writer ) const
{
    writer.write(K_);
    writer.write(static_cast<digraphE<basevector> const&>(*this));
}

void HyperBasevector::readBinary( BinaryReader& reader )
{
    reader.read(&K_);
    reader.read(static_cast<digraphE<basevector>*>(this));
}

void HyperBasevectorX::writeBinary( BinaryWriter& writer ) const
{
    writer.write(K_);
    writer.write(static_cast<digraphEX<basevector> const&>(*this));
}

void HyperBasevectorX::readBinary( BinaryReader& reader )
{
    reader.read(&K_);
    reader.read(static_cast<digraphEX<basevector>*>(this));
}

// Check that the edge adjacencies in this HyperBasevector make sense.
// If they don't, this causes a FatalErr.

Bool HyperBasevector::TestValid( const Bool exit ) const
{

     if ( !digraphE<basevector>::TestValid(exit) ) return False;

     const String err = "This HyperBasevector is internally inconsistent.";

     for ( int e = 0; e < EdgeObjectCount( ); e++ )
     {    if ( EdgeLengthKmers(e) <= 0 )
          {    cout << "HyperBasevector edge " << e << " has length "
                    << EdgeLengthKmers(e) << " kmers." << endl;
               if ( !exit ) return False;
               FatalErr(err);    }    }

  // For every vertex, determine the (K-1)-mer adjacency implied by that vertex.
  for ( int v = 0; v < N(); v++ ) {
    const vec<int> & to = ToEdgeObj(v), from = FromEdgeObj(v);
    if ( to.empty() && from.empty() ) continue;
    basevector adj;
    if ( to.empty() ) adj.SetToSubOf( EdgeObject( from[0] ), 0, K_-1 );
    else adj.SetToSubOf( EdgeObject( to[0] ), EdgeObject( to[0] ).size() - K_ + 1, K_-1 );

    // Verify that this (K-1)-mer appears in every edge to/from this vertex.
    Bool broken = False;
    for ( size_t i = 0; i < to.size(); i++ )
      if ( adj != basevector( EdgeObject( to[i] ), EdgeObject( to[i] ).size() - K_ + 1, K_-1 ) )
	broken = True;

    for ( size_t i = 0; i < from.size(); i++ )
      if ( adj != basevector( EdgeObject( from[i] ), 0, K_-1 ) )
	broken = True;

    if (broken)
    {    cout << "\nProblem in HyperBasevector at vertex " << v
              << ", entering edges ";
         for ( int j = 0; j < To(v).isize( ); j++ )
         {    if ( j > 0 ) cout << ",";
              int e = EdgeObjectIndexByIndexTo( v, j );
              cout << e << "[l=" << EdgeLengthKmers(e) << "]";    }
         cout << " and exiting edges ";
         for ( int j = 0; j < From(v).isize( ); j++ )
         {    if ( j > 0 ) cout << ",";
              int e = EdgeObjectIndexByIndexFrom( v, j );
              cout << e << "[l=" << EdgeLengthKmers(e) << "]";    }
         cout << "." << endl;
         if ( !exit ) return False;
         FatalErr(err);    }
  }
  return True;
}


void HyperBasevector::RemoveSmallComponents( int min_kmers )
{    vec< vec<int> > comps;
     Components(comps);
     vec<int> o, keep;
     for ( size_t i = 0; i < comps.size(); i++ )
     {    const vec<int>& o = comps[i];
          int nkmers = 0;
          for ( size_t j = 0; j < o.size(); j++ )
          {    int v = o[j];
               for ( size_t t = 0; t < From(v).size(); t++ )
                    nkmers += EdgeObjectByIndexFrom( v, t ).size() - K( ) + 1;    }
          if ( nkmers < min_kmers ) continue;
          keep.append(o);    }
     HyperBasevector h( *this, keep );
     *this = h;    }

void HyperBasevector::LowerK( int newK )
{    ForceAssertLe( newK, K( ) );
     for ( int i = 0; i < EdgeObjectCount( ); i++ )
     {    ForceAssertGt( EdgeObject(i).size( ), 0u );
          EdgeObjectMutable(i).resize( EdgeObject(i).size( ) + newK - K( ) );    }
     K_ = newK;    }

void HyperBasevector::ReduceK( int newK )
{    int delta = K( ) - newK;
     ForceAssertGe( delta, 0 );
     ForceAssertEq( delta % 2, 0 );
     for ( int v = 0; v < N( ); v++ )
     {    if ( From(v).empty( ) || To(v).empty( ) );
          else if ( From(v).solo( ) && To(v).size( ) > 1 )
          {    for ( int j = 0; j < To(v).isize( ); j++ )
               {    int e = EdgeObjectIndexByIndexTo( v, j );
                    EdgeObjectMutable(e).resize(
                         EdgeLengthBases(e) - delta );    }    }
          else if ( To(v).solo( ) && From(v).size( ) > 1 )
          {    for ( int j = 0; j < From(v).isize( ); j++ )
               {    int e = EdgeObjectIndexByIndexFrom( v, j );
                    EdgeObjectMutable(e).SetToSubOf( EdgeObject(e), delta,
                         EdgeLengthBases(e) - delta );    }    }
          else
          {    for ( int j = 0; j < To(v).isize( ); j++ )
               {    int e = EdgeObjectIndexByIndexTo( v, j );
                    EdgeObjectMutable(e).resize( EdgeLengthBases(e) - delta/2 );    }
               for ( int j = 0; j < From(v).isize( ); j++ )
               {    int e = EdgeObjectIndexByIndexFrom( v, j );
                    EdgeObjectMutable(e).SetToSubOf( EdgeObject(e), delta/2,
                         EdgeLengthBases(e) - delta );    }    }    }
     K_ = newK;    }

void HyperBasevector::GenerateVecbasevector( vecbvec &bases ) const
{
  bases.clear( );

  // Reserve memory.
  longlong n_reads = 0;
  longlong n_bases = 0;
  vec< vec<int> > components;
  this->ComponentEdges( components );
  for (size_t ii=0; ii<components.size(); ii++) {
    n_reads += components[ii].size();
    for (size_t jj=0; jj<components[ii].size(); jj++)
      n_bases += this->EdgeObject( components[ii][jj] ).size( );
  }
  bases.Reserve( n_bases / 16 + n_reads, n_reads );

  // Fill bases.
  for (size_t ii=0; ii<components.size(); ii++)
    for (size_t jj=0; jj<components[ii].size(); jj++)
      bases.push_back( this->EdgeObject( components[ii][jj] ) );
}

void
HyperBasevector::DumpFasta( const String& fn, Bool edge_alpha /* = True */ ) const
{
  Ofstream( out, fn );
  vec<int> to_left, to_right;
  ToLeft(to_left), ToRight(to_right);
  for ( int e = 0; e < EdgeObjectCount( ); e++ ) {
    int v = to_left[e], w = to_right[e];
    if ( edge_alpha )
	out << ">edge_" << BaseAlpha(e);
    else
	out << ">" << e;
    out << " " << v << ":" << w << "\n";
    EdgeObject(e).Print(out);
  }
}

void
HyperBasevector::DumpFastb( const String& fn ) const
{
  IncrementalWriter<basevector> basesOut( fn.c_str( ) );
  for ( int e = 0; e < EdgeObjectCount( ); e++ ) {
    basesOut.add( EdgeObject(e) );
  }
  basesOut.close();
}

Bool operator==( const HyperBasevector& h1, const HyperBasevector& h2 )
{    if ( h1.K( ) != h2.K( ) ) return False;
     return (const digraphE<basevector>&) h1
          == (const digraphE<basevector>&) h2;    }

void HyperBasevector::RemoveUnneededVertices( )
{    for ( int i = 0; i < N( ); i++ )
     {    if ( From(i).size( ) == 1 && To(i).size( ) == 1 && From(i)[0] != i )
          {    basevector p = EdgeObjectByIndexTo( i, 0 );
               p.resize( p.size( ) - K( ) + 1 );
               p.append( EdgeObjectByIndexFrom( i, 0 ) );
               JoinEdges( i, p );    }    }
     RemoveEdgelessVertices( );    }

basevector HyperBasevector::EdgePathToBases( const vec<int>& e ) const
{    ForceAssert( e.nonempty( ) );
     return Cat(e);    }

void HyperBasevector::DeleteReverseComplementComponents( )
{    vec<int> to_rc( EdgeObjectCount( ), -1 );
     {    vec<basevector> edges;
          vec<int> ids( EdgeObjectCount( ), vec<int>::IDENTITY );
          for ( int i = 0; i < EdgeObjectCount( ); i++ )
               edges.push_back( EdgeObject(i) );
          SortSync( edges, ids );
          for ( int i = 0; i < EdgeObjectCount( ); i++ )
          {    basevector q = EdgeObject(i);
               q.ReverseComplement( );
               int pos = BinPosition( edges, q );
               if ( pos >= 0 ) to_rc[i] = ids[pos];    }    }
     vec< vec<int> > components;
     ComponentEdges(components);
     Sort(components);
     vec<int> rc_to_delete;
     for ( size_t i = 0; i < components.size(); i++ )
     {    vec<int> rc;
          for ( size_t j = 0; j < components[i].size(); j++ )
               rc.push_back( to_rc[ components[i][j] ] );
          Sort(rc);
          int p = BinPosition( components, rc );
          if ( p > (int)i )
          {    for ( size_t j = 0; j < components[p].size(); j++ )
                    rc_to_delete.push_back( components[p][j] );    }    }
     Sort(rc_to_delete);
     DeleteEdges(rc_to_delete);    }

// FindIsomorphicComponents2: virtually identically to the HyperKmerPath
// function of the same name.

void HyperBasevector::FindIsomorphicComponents2(
     equiv_rel& comp1, equiv_rel& comp2 ) const
{
  comp1.Initialize( N() );
  ComponentRelation( comp1 );
  typedef vec<int> component_signature_t;
  vec< component_signature_t > componentSignatures;
  vec<int> componentReps;
  comp1.OrbitRepsAlt( componentReps );
  vec<int> reps = componentReps;
  comp2.Initialize( componentReps.size( ) );
  for ( size_t i = 0 ; i < componentReps.size(); i++ )
  {    int thisComponentRep = componentReps[i];
       vec<int> thisComponentVerts;
       comp1.Orbit( thisComponentRep, thisComponentVerts );
       component_signature_t thisComponentSig;
       thisComponentSig.push_back( thisComponentVerts.size() );
       int thisComponentNumEdges = 0;
       for ( size_t j = 0; j < thisComponentVerts.size(); j++ )
            thisComponentNumEdges += From(thisComponentVerts[j]).size();
       thisComponentSig.push_back( thisComponentNumEdges );
       vec<int> thisComponentOutDegrees;
       for ( size_t j = 0; j < thisComponentVerts.size(); j++ )
            thisComponentOutDegrees.push_back( From(thisComponentVerts[j]).size() );
       Sort( thisComponentOutDegrees );
       thisComponentSig.append( thisComponentOutDegrees );
       vec<int> thisComponentInDegrees;
       for ( size_t j = 0; j < thisComponentVerts.size(); j++ )
           thisComponentInDegrees.push_back( To(thisComponentVerts[j]).size() );
       Sort( thisComponentInDegrees );
       thisComponentSig.append( thisComponentInDegrees );
       vec<int> thisComponentEdgeLens;
       for ( size_t j = 0; j < thisComponentVerts.size(); j++ )
       {    const vec<int>& fromThisVrtx = FromEdgeObj(thisComponentVerts[j]);
            for ( size_t k = 0; k < fromThisVrtx.size(); k++ )
	    thisComponentEdgeLens.push_back(
                 EdgeObject( fromThisVrtx[k] ).isize( ) - K( ) + 1 );    }
       Sort( thisComponentEdgeLens );
       thisComponentSig.append( thisComponentEdgeLens );
       componentSignatures.push_back( thisComponentSig );    }
  SortSync( componentSignatures, componentReps );
  for ( size_t i = 0; i+1 < componentSignatures.size(); i++ )
  {    if ( componentSignatures[i] == componentSignatures[i+1] &&
  	    ComponentsAreIsomorphic( *this, componentReps[i],
                *this, componentReps[i+1] ) )
       {    comp2.Join( Position( reps, componentReps[i] ),
                        Position( reps, componentReps[i+1] ) );    }    }    }

// ComponentsAreIsomorphic: same as for HyperKmerPath.

bool ComponentsAreIsomorphic( const HyperBasevector& hkp1, int seed_vx_1,
			      const HyperBasevector& hkp2, int seed_vx_2,
			      vec<int>* p_iso ) {
  vec<int> local_iso;
  if( p_iso == NULL ) p_iso = &local_iso;

  p_iso->resize( hkp1.N(), -1 );
  // (*p_iso)[i]=-1 means hkp1 vx i unused; j>=0 means matches hkp2 vx j.
  vec<Bool> used2( hkp2.N(), false );

  // The hypothesis <i,j> means vertices i and j ought to match
  vec< pair<int,int> > hypotheses;
  // The total number of hypotheses will be two per edge in the
  // connected component (but generally not all on the stack at once).
  hypotheses.reserve( 2 * hkp1.EdgeObjectCount() );
  // Initial hypothesis: the given seed vertices match.
  hypotheses.push_back( make_pair(seed_vx_1,seed_vx_2) );

  while( ! hypotheses.empty() ) {

    int v1 = hypotheses.back().first, v2=hypotheses.back().second;
    hypotheses.pop_back();

    if( (*p_iso)[v1]==v2 )                     // already knew this hypothesis
      continue;
    else if( used2[v2] || (*p_iso)[v1]!=-1 )   // knew a contradictory hypothesis
      return false;

    // Otherwise:
    // (1) Assume this hypothesis is true
    (*p_iso)[v1] = v2;
    used2[v2] = true;

    // (2) check that the edges in/out of these vertices match up
    vec<int> from_v1 = hkp1.From(v1), to_v1 = hkp1.To(v1);
    vec<int> from_v2 = hkp2.From(v2), to_v2 = hkp2.To(v2);

    if( from_v1.size() != from_v2.size() || to_v1.size() != to_v2.size() )
      return false;

    vec<basevector> from_v1_edges( from_v1.size() );
    vec<basevector> from_v2_edges( from_v2.size() );
    vec<basevector> to_v1_edges( to_v1.size() );
    vec<basevector> to_v2_edges( to_v2.size() );

    for(size_t i=0; i<from_v1.size(); i++)
      from_v1_edges[i] = hkp1.EdgeObjectByIndexFrom(v1,i);
    for(size_t i=0; i<from_v2.size(); i++)
      from_v2_edges[i] = hkp2.EdgeObjectByIndexFrom(v2,i);

    for(size_t j=0; j<to_v1.size(); j++)
      to_v1_edges[j] = hkp1.EdgeObjectByIndexTo(v1,j);
    for(size_t j=0; j<to_v2.size(); j++)
      to_v2_edges[j] = hkp2.EdgeObjectByIndexTo(v2,j);

    SortSync( from_v1_edges, from_v1 );
    SortSync( from_v2_edges, from_v2 );
    SortSync( to_v1_edges, to_v1 );
    SortSync( to_v2_edges, to_v2 );

    if( from_v1_edges != from_v2_edges || to_v2_edges != to_v2_edges )
      return false;

    ForceAssert( from_v1_edges.UniqueOrdered() );
    ForceAssert( from_v2_edges.UniqueOrdered() );
    ForceAssert( to_v1_edges.UniqueOrdered() );
    ForceAssert( to_v2_edges.UniqueOrdered() );

    for( size_t i=0; i < from_v1.size(); i++ )
      hypotheses.push_back( make_pair(from_v1[i],from_v2[i]) );
    for( size_t j=0; j < to_v1.size(); j++ )
      hypotheses.push_back( make_pair(to_v1[j],to_v2[j]) );


    // (3) formulate all the corollary hypotheses imposed by edge matching
    for( size_t i=0; i < from_v1.size(); i++ )
      hypotheses.push_back( make_pair(from_v1[i],from_v2[i]) );
    for( size_t j=0; j < to_v1.size(); j++ )
      hypotheses.push_back( make_pair(to_v1[j],to_v2[j]) );

  } // end of processing of hypotheses.back()

  // If no hypotheses gave rise to a contradiction,
  // then we have an isomorphism.  Hooray!
  return true;
}


#include "graph/DigraphTemplate.h"
template Bool digraphE<BaseVec>::TestValid(const Bool) const;
template digraphE<BaseVec>::digraphE();
template digraphE<BaseVec>::digraphE(vec<BaseVec> const&, const ConstructorBehavior);
template digraphE<BaseVec>::digraphE(const vec<vec<int> >&, const vec<vec<int> >&, const vec<BaseVec>&, const vec<vec<int> >&, const vec<vec<int> >&, const Bool);
template int digraphE<BaseVec>::AddEdge(int, int, BaseVec const&);
template void digraphE<BaseVec>::AddVertices(int);
template void digraphE<BaseVec>::Clear();
template void digraphE<BaseVec>::ComponentEdges( vec< vec<int> >& edges ) const;
template void digraphE<BaseVec>::ComponentsE(vec<vec<int> >&) const;
template void digraphE<BaseVec>::DeleteEdgeFrom(int, int);

template void digraphE<BaseVec>::DeleteEdges(vec<int> const&);
template void digraphE<BaseVec>::DeleteEdges( vec<int> const&, vec<int> const& );

template void digraphE<BaseVec>::DeleteEdgesAtVertex(int);
template void digraphE<BaseVec>::DumpGraphML( const String& ) const;
template int digraphE<BaseVec>::UsedCount() const;

template const BaseVec& digraphE<BaseVec>::EdgeObject(int) const;
template const BaseVec& digraphEX<BaseVec>::EdgeObject(int) const;

template void digraphEX<BaseVec>::Used(vec<Bool>&) const;
template int digraphEX<BaseVec>::UsedCount() const;

template BaseVec const& digraphE<BaseVec>::EdgeObjectByIndexFrom(int, int) const;
template BaseVec& digraphE<BaseVec>::EdgeObjectByIndexFromMutable(int, int);
template BaseVec const& digraphE<BaseVec>::EdgeObjectByIndexTo(int, int) const;
template BaseVec& digraphE<BaseVec>::EdgeObjectByIndexToMutable(int, int);
template int digraphE<BaseVec>::EdgeObjectCount() const;
template int digraphE<BaseVec>::EdgeObjectIndexByIndexFrom(int, int) const;
template vec<BaseVec> const& digraphE<BaseVec>::Edges() const;
template vec<BaseVec>& digraphE<BaseVec>::EdgesMutable();
template BaseVec& digraphE<BaseVec>::EdgeObjectMutable(int);
template int digraphE<BaseVec>::EdgeObjectIndexByIndexTo(int, int) const;
template vec<int> digraphE<BaseVec>::EdgesBetween( const int, const int) const;
template vec<int> digraphE<BaseVec>::EdgesSomewhereBetween( const int, const int) const;
template vec<int> digraphE<BaseVec>::EdgesConnectedTo(vec<int> const&) const;
template vec<vec<int> > const& digraphE<BaseVec>::FromEdgeObj() const;
template vec<int> const& digraphE<BaseVec>::FromEdgeObj(int) const;
template vec<vec<int> >& digraphE<BaseVec>::FromEdgeObjMutable();
template vec<int>& digraphE<BaseVec>::FromEdgeObjMutable(int);
template void digraphE<BaseVec>::Initialize(vec<vec<int> > const&, vec<vec<int> > const&, vec<BaseVec> const&, vec<vec<int> > const&, vec<vec<int> > const&, const Bool);
template int digraphE<BaseVec>::InputFromOutputTo(int, int) const;
template int digraphE<BaseVec>::InputToOutputFrom(int, int) const;
template void digraphE<BaseVec>::JoinEdges(int, const BaseVec&);
template void digraphE<BaseVec>::PopBubbles(const vec<int>&);
template void digraphE<BaseVec>::PopHyperBubbles(const vec<int>&);

template void digraphE<BaseVec>::PrettyDOT( ostream& out,
     const vec<double>& lengths, const edge_label_info, Bool label_contigs,
     Bool label_vertices, const vec<int>* componentsToPrint,
     const vec<String> *label_contigs_extra, const vec<int>* verticesToPrint,
     const vec<Bool>* dashed, const vec<Bool>* invisible,
     const vec<String>* edge_color, const vec<int>* pen_widths, const String,
     const double, const double, const double, const double ) const;

template void digraphEX<BaseVec>::PrettyDOT( ostream& out,
     const vec<double>& lengths, const edge_label_info, Bool label_contigs,
     Bool label_vertices, const vec<int>* componentsToPrint,
     const vec<String> *label_contigs_extra, const vec<int>* verticesToPrint,
     const vec<Bool>* dashed, const vec<Bool>* invisible,
     const vec<String>* edge_color, const vec<int>* pen_widths, const String,
     const double, const double, const double, const double ) const;

template vec<int> digraphE<BaseVec>::RemoveDeadEdgeObjects();
template void digraphE<BaseVec>::RemoveEdgelessVertices( );
template void digraphE<BaseVec>::RemoveEdgelessVertices(vec<int> const&);
template void digraphE<BaseVec>::RemoveUnneededVertices();
template void digraphE<BaseVec>::Reverse();
template void digraphE<BaseVec>::SplitEdge(int, int, BaseVec const&, BaseVec const&);
template digraphE<BaseVec> digraphE<BaseVec>::Subgraph(const vec<int>&) const;
template vec<vec<int> > const& digraphE<BaseVec>::ToEdgeObj() const;
template vec<int> const& digraphE<BaseVec>::ToEdgeObj(int) const;
template vec<vec<int> >& digraphE<BaseVec>::ToEdgeObjMutable();
template vec<int>& digraphE<BaseVec>::ToEdgeObjMutable(int);
template void digraphE<BaseVec>::ToLeft(vec<int>&) const;
template void digraphE<BaseVec>::ToRight(vec<int>&) const;
template void digraphE<BaseVec>::TransferEdges(int, int, unsigned char );
template void digraphE<BaseVec>::Used(vec<unsigned char>&) const;

template void digraphE<BaseVec>::readBinary(BinaryReader&);
template void digraphE<BaseVec>::writeBinary(BinaryWriter&) const;

template void digraphEX<BaseVec>::readBinary(BinaryReader&);
template void digraphEX<BaseVec>::writeBinary(BinaryWriter&) const;

template digraphE<BaseVec>::digraphE(digraphE<BaseVec>::ConstructorType1, digraphE<BaseVec> const&, vec<int, std::allocator<int> > const&);

template void digraphE<BaseVec>::Initialize(digraphE<BaseVec>::ConstructorType1, digraphE<BaseVec> const&, vec<int, std::allocator<int> > const&);

template bool operator==( const digraphE<BaseVec>& g1, const digraphE<BaseVec>& g2 );

template digraphE<BaseVec>::edge_label_info::edge_label_info(ConstructorBehavior, Bool, Bool, vec<String> const*);
template digraphE<BaseVec>::edge_label_info::edge_label_info(ConstructorBehavior, vec<String> const*);
template Bool digraphE<BaseVec>::EdgePaths(int, int, vec<vec<int> >&, const int, const int, const int, const Bool) const;

template void digraphE<BaseVec>::Initialize(int);

template void digraphE<BaseVec>::GiveEdgeNewFromVx(int, int, int);
template void digraphE<BaseVec>::GiveEdgeNewToVx(int, int, int);

template void digraphE<BaseVec>::ReverseComponent(int);

template void digraphE<BaseVec>::DeleteEdgesFrom(int, vec<int> const&);
template void digraphE<BaseVec>::DeleteEdgesTo(int, vec<int> const&);
template void digraphE<BaseVec>::DeleteEdgeTo(int, int);
template void digraphE<BaseVec>::SplayVertex(int);

template void digraphE<BaseVec>::InitialEdges(vec<int>&) const;
template void digraphE<BaseVec>::TerminalEdges(vec<int>&) const;

template digraphE<BaseVec>::digraphE(digraphE<BaseVec>::ConstructorType2, digraphE<BaseVec> const&, vec<int> const&, vec<int> const&, vec<int>const& );

template const basevector& digraphE<BaseVec>::O(int) const;

basevector HyperBasevector::Cat( const int e1, const int e2 ) const
{    return TrimCat( K( ), EdgeObject(e1), EdgeObject(e2) );    }
basevector HyperBasevector::Cat( const int e1, const int e2, const int e3 ) const
{    return TrimCat( K( ), EdgeObject(e1), EdgeObject(e2), EdgeObject(e3) );    }
basevector HyperBasevector::Cat( const int e1, const int e2, const int e3,
     const int e4 ) const
{    return TrimCat( K( ), EdgeObject(e1), EdgeObject(e2), EdgeObject(e3),
          EdgeObject(e4) );    }
basevector HyperBasevector::Cat( const int e1, const int e2, const int e3,
     const int e4, const int e5 ) const
{    return TrimCat( K( ), EdgeObject(e1), EdgeObject(e2), EdgeObject(e3),
          EdgeObject(e4), EdgeObject(e5) );    }
basevector HyperBasevector::Cat( const vec<int>& e ) const
{    basevector x = EdgeObject( e[0] );
     for ( int j = 1; j < e.isize( ); j++ )
          x = TrimCat( K( ), x, EdgeObject( e[j] ) );
     return x;    }

basevector HyperBasevectorX::Cat( const vec<int>& e ) const
{    basevector x = EdgeObject( e[0] );
     for ( int j = 1; j < e.isize( ); j++ )
          x = TrimCat( K( ), x, EdgeObject( e[j] ) );
     return x;    }

template <class TEdgeGraph>
static int64_t EdgeGraphCheckSum(TEdgeGraph const* pG) {
     int64_t x = 0;
     for ( int e = 0; e < pG->E( ); e++ )
     {    const basevector& b = pG->EdgeObject(e);
          for ( int j = 0; j < b.isize( ); j++ )
               x += ( b[j] + 1 ) * ( e + 1 ) * ( j + 1 );    }
     x += ((typename TEdgeGraph::Superclass&)(*pG)).CheckSum( );
     return x;
}


int64_t HyperBasevector::CheckSum( ) const
{    return EdgeGraphCheckSum( this ); }

int64_t HyperBasevectorX::CheckSum() const
{    return EdgeGraphCheckSum( this ); } 

template int64_t digraphE<BaseVec>::CheckSum() const;
template int64_t digraphEX<BaseVec>::CheckSum() const;

template vec<int>
digraphE<BaseVec>::EdgesBoundedBy(int, int, vec<int, std::allocator<int> > const&, vec<int, std::allocator<int> > const&) const;

template Bool digraphE<basevector>::EdgePaths( const vec<int>& to_left,
     const vec<int>& to_right,
     const int v, const int w, vec< vec<int> >& paths, const int max_copies,
     const int max_paths, const int max_iterations, const Bool ) const;

template int digraphE<BaseVec>::IFrom(int, int) const;
template int digraphE<BaseVec>::ITo(int, int) const;
template int digraphE<BaseVec>::E() const;

template void digraphE<BaseVec>::DeleteEdgesParallel( const vec<int>& );
template void digraphE<BaseVec>::DeleteEdgesParallel( const vec<Bool>& );

template void digraphE<BaseVec>::ComponentsEFast( vec<vec<int>>& ) const;

template digraphEX<BaseVec>::digraphEX(digraphE<BaseVec> const&);

template
void DistancesToEndArr<BaseVec>(digraphE<BaseVec> const&, vec<int, std::allocator<int> > const&, int, unsigned char, vec<int, std::allocator<int> >&);

void HyperBasevector::Involution( vec<int>& inv )
{    vecbasevector edges( EdgeObjectCount( ) );
     for ( int e = 0; e < EdgeObjectCount( ); e++ )
          edges[e] = EdgeObject(e);
     vec<int> x1( E( ), vec<int>::IDENTITY ), x2( E( ), vec<int>::IDENTITY );
     inv.resize( edges.size( ) );
     ParallelSort( x1, [&edges](int i1, int i2){ return edges[i1] < edges[i2]; } );
     #pragma omp parallel for
     for ( int i = 0; i < (int) edges.size( ); i++ )
          edges[i].ReverseComplement( );
     ParallelSort( x2, [&edges](int i1, int i2){ return edges[i1] < edges[i2]; } );
     for ( int i = 0; i < E( ); i++ )
          inv[ x1[i] ] = x2[i];    }

template void digraphE<basevector>::Initialize(
     const ConstructorType2 constructor_type, const digraphE& g, const vec<int>& ed,
     const vec<int>& to_left, const vec<int>& to_right );

template Bool digraphEX<basevector>::EdgePaths( const int, const int,
     vec< vec<int> >&, const int, const int, const int, const Bool ) const;
