///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Fastavector.h"
#include "Bitvector.h"
#include "graph/Digraph.h"
#include "system/System.h"
#include <cstring>
#include <istream>
#include <string>

// Split this into chunks, separated by gaps ('n'), and return each chunk as a
// gapless fastavector.
// TODO: generalize this into a templatized STL algorithm.
vec<fastavector>
fastavector::SplitOnGaps() const
{
    vec<fastavector> chunks;
    const_iterator stop(end());
    for ( const_iterator itr(begin()); itr != stop; ++itr )
    {
        if ( *itr != 'n' ) // found a chunk's start
        {
            const_iterator itr2(itr);
            while ( ++itr2 != stop && *itr2 != 'n' )
                ; // just keep going until we find the chunk's end

            chunks.push_back(fastavector(itr,itr2));

            if ( itr2 == stop )
                break;
            itr = itr2;
        }
    }
    return chunks;
}

// Returns a basevector, and sets a bitvector of ambiguous bases, if supplied.
basevector 
fastavector::ToBasevector( bitvector* ambiguous) const
{
  basevector bv(size());
  if (ambiguous) {
    ambiguous->resize(size());
    ambiguous->Zero();
  }

  // loop over bases
  for ( size_type idx = 0; idx < size(); ++idx ) {
    GeneralizedBase const& gb = GeneralizedBase::fromChar((*this)[idx]);
    char const* endBase = gb.end();
    char const* curBase = gb.bases();

    // set basevec base to first possible base code
    bv.set(idx, Base::char2Val(*curBase));

    // mark if this base is ambiguous
    if (endBase - curBase > 1 && ambiguous)
      ambiguous->set(idx, 1);
  }
  return bv;
}

basevector 
fastavector::ToBasevectorRandom( ) const
{
  basevector bv(size());

  // loop over bases
  for ( size_type idx = 0; idx < size(); ++idx ) {
    GeneralizedBase const& gb = GeneralizedBase::fromChar((*this)[idx]);
    char const* endBase = gb.end();
    char const* curBase = gb.bases();

    // set basevec base to first possible base code
    bv.set(idx, Base::char2Val(*curBase));

    if (endBase - curBase > 1){
      bv.Set(idx,gb.random());
    }
  }
  return bv;
}



// Returns the set of all basevectors that this fastavector could
// possibly represent.  Fails if this fastavector has any gaps.
// Returns an empty set if the number of possibilities is more than max.
vecbasevector
fastavector::AllBasevectors( size_t maxVecs ) const
{
  ForceAssert( !HasGaps() );
  vecbasevector all;
  
  size_t nVecs = 1;
  const_iterator stop(end());
  for ( const_iterator itr(begin()); itr != stop; ++itr ) {
    GeneralizedBase const& gb = GeneralizedBase::fromChar(*itr);
    nVecs *= gb.end() - gb.bases();
    if ( nVecs > maxVecs )
      return all; // too many options - return an empty vecbasevector
  }
  
  all.reserve(nVecs);
  
  // Create an initial copy -- ambiguity codes will be squashed randomly,
  // which is OK.  We'll nail 'em down later.
  all.push_back(bvec(begin(),end(),GenCharToRandomBaseMapper()));
  
  // If there is no ambiguity anywhere, we can return the single vecbasevector.
  if ( nVecs == 1 ) return all;
  
  // Find each ambiguous position, with more than one base option
  for ( size_type idx = 0; idx < size(); ++idx ) {
    
    GeneralizedBase const& gb = GeneralizedBase::fromChar((*this)[idx]);
    char const* endBase = gb.end();
    char const* curBase = gb.bases();
    if ( endBase - curBase == 1 ) continue; // only one base option here
    
    // take the current set of bvecs and set the code in the
    // ambiguous position to the first legal code
    bvec::value_type code = Base::char2Val(*curBase);
    vecbvec::iterator end(all.end());
    for ( vecbvec::iterator itr(all.begin()); itr != end; ++itr )
      itr->set(idx,code);
    
    // for the remaining legal codes in the ambiguity set at this idx
    for ( ++curBase; curBase < endBase; ++curBase ) {
      // append a copy of the current set, as it was initially,
      // for each additional legal base, and set the code in the
      // ambiguous position to the next legal value
      code = Base::char2Val(*curBase);
      
      for ( vecbvec::iterator itr(all.begin()); itr != end; ++itr ) {
	// this could theoretically invalidate the iterators,
	// but it doesn't because we've reserved enough space
	all.push_back(*itr);
	all.back().set(idx,code);
      }
    }
  }
  
  return all;
}

// Returns the set of all kmers contained in AllBasevectors (but if a kmer
// location has more than small_max possibilities, ignore it.)
// Much more likely to return a reasonably sized set than AllBasevectors.
vecbasevector
fastavector::AllKmers( fastavector::size_type K, size_t small_max ) const
{
  vecbasevector result;
  if ( size() < K ) return result; // empty result
  
  size_type nKmers = size() - K + 1;
  result.reserve(nKmers);
  
  const_iterator stop(begin(nKmers));
  const_iterator itr2(begin(K));
  
  bool unambig = false;
  
  for ( const_iterator itr(begin()); itr != stop; ++itr, ++itr2 ) {
    
    // As a time-saving step, we use the unambig flag to keep track of when
    // the fastavector contains no ambiguous bases.  In these areas we can
    // create the basevector directly from the fastavector, instead of calling
    // AllBasevectors.
    int n_ambig;
    if ( unambig ) {
      result.push_back( basevector( itr, itr2, GenCharToRandomBaseMapper() ) );
      n_ambig = 1;
    }
    else {
      vecbasevector vbv(0);
      vbv = fastavector(itr,itr2).AllBasevectors(small_max);
      result.append( vbv.begin(), vbv.end() );
      n_ambig = vbv.size();
    }
    
    // Update the unambig flag.
    unambig = false;
    if ( n_ambig == 1 )
      if ( itr2 != end() ) {
	GeneralizedBase const& gb = GeneralizedBase::fromChar(*itr2);
	if ( gb.end() - gb.bases() == 1 )
	  unambig = true;
      }
  }
  
  return result;
}

void fastavector::Print( ostream& out, const String& id ) const
{
    out << '>' << id;
    for ( size_type i = 0; i < size(); ++i )
    {
        if ( !(i % 80) ) out << '\n';
        out << (*this)[i];
    }
    out << '\n';
}

void fastavector::ReverseComplement()
{
    GeneralizedBase::reverseComplement(begin(), end());
}

void LoadFromFastaFile( const String& f, vec<fastavector>& vfv )
{
    ifstream in(f.c_str());
    std::string buf;
    buf.reserve( 8192 );
    fastavector fv;
    fv.reserve( 100000 );

    vfv.clear();

    bool dataPresent = false;
    while ( getline(in,buf) )
    {
        size_t sz = buf.size();
        if ( !sz )
            continue; // NON-STRUCTURED!  ignore blank lines

        if ( buf[0] != '>' ) // if not a FASTA comment
        {
            if ( fv.size() + sz > fv.capacity() )
                // grow rapidly to avoid a realloc per line on whole-genome
                // files (which typically have a gigantic chr1 coming first)
                fv.reserve( 4*fv.capacity() );

            // buffer the FASTA
            fv.append(buf.begin(),buf.end());
        }
        else if ( dataPresent )
        {
            vfv.push_back(fv);
            fv.clear();
        }
        dataPresent = true;
    }

    if ( !dataPresent )
    {
        if ( !in.eof() )
            FatalErr("Unable to read FASTA file " << f);
    }
    else
    {
        vfv.push_back(fv);
        if ( !in.eof() )
            FatalErr("Unable to read FASTA file " << f << " to completion");
    }
}

void LoadFromFastaFile( const String& f, vec<fastavector>& vfv, vec<String>& vnames )
{
    ifstream in(f.c_str());
    std::string buf;
    buf.reserve( 8192 );
    fastavector fv;
    fv.reserve( 100000 );

    vfv.clear();

    bool dataPresent = false;
    while ( getline(in,buf) )
    {
        size_t sz = buf.size();
        if ( !sz )
            continue; // NON-STRUCTURED! ignore blank lines

        if ( buf[0] != '>' ) // if not a FASTA comment
        {
            if ( fv.size() + sz > fv.capacity() )
                // grow rapidly to avoid a realloc per line on whole-genome
                // files (which typically have a gigantic chr1 coming first)
                fv.reserve(4 * fv.capacity());

            // buffer the FASTA
            fv.append(buf.begin(), buf.end());

            // if we're jumping right into FASTA without having seen a comment
            if ( !dataPresent )
                vnames.push_back("");
        }
        else
        {
            String name(buf.begin()+1,buf.end());
            DeleteTrailingWhiteSpace(name);
            vnames.push_back(name);
            if ( dataPresent )
            {
                vfv.push_back(fv);
                fv.clear();
            }
        }
        dataPresent = true;
    }

    if ( !dataPresent )
    {
        if ( !in.eof() )
            FatalErr("Unable to read FASTA file " << f);
    }
    else
    {
        vfv.push_back(fv);
        if ( !in.eof() )
            FatalErr("Unable to read FASTA file " << f << " to completion");
    }

    ForceAssertEq( vfv.size(), vnames.size() );
}



// Load multiple, possibly overlapping regions from an indexed fasta file
// we parse the index and then the fasta.
//
// Output is vec<fastavector> of contigs and vec<String> of contig names
// in the SAME ORDER as the specified regions vec.
//
// optional: pad specifies an amount to pad the region on either side, whenever possible
// optional: padded are pairs indicating how much padding was actually used at the start and end
//
void LoadRegionsFromIndexedFastaFile( const String& fasta, const String& index,
	     const vec<FastaRegion>& regions,
	     vec<fastavector>& v,
	     int pad,
	     vec<std::pair<int,int>>* padded
	     )
{
    ForceAssertGe(pad, 0);
    fastaindex fai( index );
    ifstream in( fasta );
    string buf;
    fastavector fv;

    if ( !in.is_open() ) FatalErr("error opening FASTA file " << fasta );

    for ( auto iter_region = regions.begin(); iter_region != regions.end(); ++iter_region ) {

	// chr:start-end
	const String& chr = std::get<0>(*iter_region);
	const int start = std::get<1>(*iter_region);
	const int end = std::get<2>(*iter_region);
	ForceAssertGe(start, 0);

	if ( !fai.has_key( chr ) )
	    FatalErr( "FASTA index " << index << " does not reference chromosome " << chr );


	if ( !in.seekg( fai[chr].start, in.beg ).good() ) {
	    FatalErr("Failure seeking in FASTA file.  Maybe needs an update?");
	}

	// pad region if necessary and possible
	int new_start = max( 0, start - pad );
	int ilen = fai[chr].len;
	int new_end = min( ilen, end + pad );

	// read file
	int so_far = 0;		// characters seen, including ones to skip
	while ( getline(in, buf) && buf.size() > 0 && buf[0] != '>' && so_far < new_end ) {
	    int bufsize = buf.size();
	    if ( so_far + bufsize > new_start ) {	// are we (almost) there yet?
		auto istart = buf.begin() + max(0, new_start - so_far );
		auto iend = min( buf.end(), istart + new_end - so_far );
		fv.append( istart, iend );
	    }
	    so_far += bufsize;
	}

	if ( static_cast<int>(fv.size()) != new_end - new_start )
	    FatalErr("short contig/chromosome " << chr << " in FASTA file " << fasta );

	// push back results
	v.push_back(fv);
	if ( padded ) {
	    padded->push_back( std::make_pair( start - new_start, new_end - end  ) );
	}
    }


}

void LoadRegionFromIndexedFastaFile( const String& fasta, const String& index,
	    const String& chr, const int start, const int end,
	     fastavector& fv )
{
    fastaindex fai( index );
    LoadRegionFromIndexedFastaFile( fasta, fai, chr, start, end, fv );
}

void LoadRegionFromIndexedFastaFile( const String& fasta, const fastaindex& fai,
	const String& chr, const int start, const int end,
	 fastavector& fv
	 )
{
    fv.clear();
    fv.reserve( end - start );
    ifstream in( fasta );
    string buf;

    if ( !in.is_open() ) FatalErr("error opening FASTA file " << fasta );

    ForceAssertGe(start, 0);

    if ( !fai.has_key( chr ) )
	FatalErr( "FASTA index " << fai.filename() << " does not reference chromosome " << chr );


    if ( !in.seekg( fai[chr].start, in.beg ).good() ) {
	FatalErr("Failure seeking in FASTA file.  Maybe needs an update?");
    }

    int ilen = fai[chr].len;

    if ( end > ilen )
	FatalErr("Requested contig end at position " << end << " but index says end is " << ilen);

    // read file
    int so_far = 0;		// characters seen, including ones to skip
    while ( getline(in, buf) && buf[0] != '>' && so_far < end ) {
	int bufsize = buf.size();
	if ( bufsize > 0 && so_far + bufsize > start ) {	// are we (almost) there yet (and not blank)?
	    auto istart = buf.begin() + max(0, start - so_far );
	    auto iend = buf.end() + min(0,  end - so_far - bufsize );
	    fv.append( istart, iend );
	}
	so_far += bufsize;
    }

    if ( static_cast<int>(fv.size()) != end - start )
	FatalErr("short contig/chromosome " << chr << " in FASTA file " << fasta );
}



void WriteScaffoldedFasta( const String &out_file,
			   const vec<fastavector> &fasta,
			   const vec<superb> &scaffolds,
			   const vec<Bool> &rc,
			   const int min_gap, 
			   const char gap_char,
			   const Bool ncbi_format)
{
  Ofstream( out, out_file );
  
  // Loop over all superbs (scaffolds).
  for ( size_t i = 0; i < scaffolds.size(); i++ ) {
    const superb& S = scaffolds[i];
    if (ncbi_format) {	  // one-based, zero padding for lexical order
      out << ">scaffold";
      out.width(5);
      out.fill('0');
      out << i+1 << "\n";
    } else {
      out << ">scaffold_" << i << "\n";
    }
    
    // Find the base sequence of this scaffold, including 'n's for
    // gaps. WARNING: gaps < min_gap will be reset at min_gap.
    vec<char> s;
    for ( int j = 0; j < S.Ntigs( ); j++ ) {
      bool need_rcing = rc.nonempty() && rc[ S.Tig(j) ];
      fastavector b_rc;
      if ( need_rcing ) {
	b_rc = fasta[S.Tig(j)];
	b_rc.ReverseComplement( );
      }
      const fastavector &b = need_rcing ? b_rc : fasta[S.Tig(j)];
      for ( unsigned int l = 0; l < b.size( ); l++ )
	s.push_back( b[l] );
      if ( j < S.Ntigs( ) - 1 )
	s.push_back_copies( gap_char, Max( min_gap, S.Gap(j) ) );
    }
    
    // Print the bases, adding line breaks where necessary.
    int printed = 1;
    for ( unsigned int j = 0; j < s.size( ); j++, printed++ ) {
      out << s[j];
      if ( printed % 80 == 0 ) out << "\n";
    }
    if ( printed == 1 || printed % 80 != 1 ) out << "\n";
  }
  
  out.close( );
}  

void RenumberAndMinimize( vec<fastavector>& f, vec<superb>& s )
{    vec<Bool> used( f.size( ), False );
     for ( size_t i = 0; i < s.size( ); i++ )
     {    for ( int j = 0; j < s[i].Ntigs( ); j++ )
          {    int m = s[i].Tig(j);
               ForceAssertLt( m, f.isize( ) );
               used[m] = True;    }    }
     vec<int> to_new( f.size( ), -1 );
     size_t count = 0;
     for ( size_t i = 0; i < f.size( ); i++ )
     {    if ( used[i] ) 
          {    f[count] = f[i];
               to_new[i] = count++;    }    }
     f.resize(count);
     for ( size_t i = 0; i < s.size( ); i++ )
     {    for ( int j = 0; j < s[i].Ntigs( ); j++ )
          {    int m = s[i].Tig(j);
               s[i].SetTig( j, to_new[m] );    }    }    }

#include "feudal/OuterVecDefs.h"
template class OuterVec<fastavector>;

#include "graph/DigraphTemplate.h"
template digraphE<fastavector>::digraphE();
template int digraphE<fastavector>::AddEdge(int, int, fastavector const&);
template void digraphE<fastavector>::AddVertices(int);
template void digraphE<fastavector>::Clear();
template void digraphE<fastavector>::ComponentEdges(vec<vec<int> >&) const;
template void digraphE<fastavector>::ComponentsE(vec<vec<int> >&) const;
template void digraphE<fastavector>::DeleteEdgeFrom(int, int);
template void digraphE<fastavector>::DeleteEdges(vec<int> const&);
template void digraphE<fastavector>::DeleteEdges(const vec<int>&, const vec<int>&);
template void digraphE<fastavector>::DeleteEdgeTo(int, int);
template const fastavector& digraphE<fastavector>::EdgeObject(int) const;
template fastavector const& digraphE<fastavector>::EdgeObjectByIndexFrom(int, int) const;
template fastavector const& digraphE<fastavector>::EdgeObjectByIndexTo(int, int) const;
template int digraphE<fastavector>::EdgeObjectCount() const;
template int digraphE<fastavector>::EdgeObjectIndexByIndexFrom(int, int) const;
template int digraphE<fastavector>::EdgeObjectIndexByIndexTo(int, int) const;
template fastavector& digraphE<fastavector>::EdgeObjectMutable(int);
template vec<fastavector> digraphE<fastavector>::EdgeObjectsBetween(int, int) const;
template vec<fastavector> const& digraphE<fastavector>::Edges() const;
template vec<int> digraphE<fastavector>::EdgesBetween(int, int) const;
template vec<fastavector>& digraphE<fastavector>::EdgesMutable();
template vec<vec<int> > const& digraphE<fastavector>::FromEdgeObj() const;
template vec<int> const& digraphE<fastavector>::FromEdgeObj(int) const;
template vec<vec<int> >& digraphE<fastavector>::FromEdgeObjMutable();
template int digraphE<fastavector>::InputFromOutputTo(int, int) const;
template int digraphE<fastavector>::InputToOutputFrom(int, int) const;
template void digraphE<fastavector>::JoinEdges(int, fastavector const&);

template void digraphE<fastavector>::PrettyDOT( ostream& out, 
     const vec<double>& lengths, const edge_label_info, Bool label_contigs, 
     Bool label_vertices, const vec<int>* componentsToPrint,
     const vec<String> *label_contigs_extra, const vec<int>* verticesToPrint,
     const vec<Bool>* dashed, const vec<Bool>* invisible,
     const vec<String>* edge_color, const vec<int>* pen_width, const String,
     const double, const double, const double, const double ) const;

template vec<int> digraphE<fastavector>::RemoveDeadEdgeObjects();
template void digraphE<fastavector>::RemoveEdgelessVertices();
template void digraphE<fastavector>::RemoveEdgelessVertices(vec<int> const&);
template void digraphE<fastavector>::Reverse();
template vec<vec<int> > const& digraphE<fastavector>::ToEdgeObj() const;
template vec<vec<int> >& digraphE<fastavector>::ToEdgeObjMutable();
template void digraphE<fastavector>::ToLeft(vec<int>&) const;
template void digraphE<fastavector>::ToRight(vec<int>&) const;
template void digraphE<fastavector>::Used(vec<unsigned char>&) const;
template void digraphE<fastavector>::readBinary(BinaryReader&);
template void digraphE<fastavector>::writeBinary(BinaryWriter&) const;

template digraphE<fastavector>::edge_label_info::edge_label_info(digraphE<fastavector>::edge_label_info::ConstructorBehavior, unsigned char, unsigned char, vec<FeudalString<char, std::char_traits<char> > > const*);

template const vec<int>& digraphE<fastavector>::ToEdgeObj(int) const;
