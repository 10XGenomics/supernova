///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// This file defines "fastavector", which stores bases in fasta format.

#ifndef FASTA_VECTOR
#define FASTA_VECTOR

#include "Basevector.h"
#include "Bitvector.h"
#include "Charvector.h"
#include "String.h"
#include "Superb.h"
#include "Vec.h"
#include "dna/Bases.h"
#include "graph/Digraph.h"
#include <cstddef>
#include <tuple>
#include <iostream>
#include <unistd.h>

typedef std::tuple<String, int, int> FastaRegion;

class fastaindex {
private:
    struct line_data {
	String info;	// extra metadata past contig name
	size_t len;	// contig length
	size_t start;	// byte offset (0-based) of contig in file
	size_t idx;
	// we're excluding the line length stuff in an faidx as we need to deal gracefully with
	// fasta files that have non-uniform line lengths.
    };
    std::map<String, line_data> _lines;
    const String _filename;		// useful for error msgs.

    static bool fexist( const string& filename ) {
	struct stat stat_buf;
	return ( !stat(filename.c_str(), &stat_buf) );
    }




public:
    static String filename_fix(const String& filename, const bool verbose = false ) {
        std::vector<String> base_names={filename};
        const std::vector<String> extensions={".fai",".faidx",".index"};
        std::vector<String> tested_names;
        auto ext_start=filename.rfind(".");
        if( ext_start != String::npos){ base_names.push_back(filename.substr(0,ext_start)); }
	for ( const auto& base: base_names ) {
	    for ( const auto& ext: extensions ) {
	        const String name = base+ext;
	        tested_names.push_back(name);
		if ( fexist( name ) )
		    return name;
	    }
	}

	if ( verbose ) {
	    std::cout << "WARNING: none of the following files have been found to serve as FASTA index:\n";
	    for (const auto& name : tested_names){ std::cout << name << "\n"; }
	    std::cout<<std::endl;
	}

	return filename;
    }
    explicit fastaindex( const String& filename ) : _filename(filename_fix(filename)) {
	ifstream in;
	in.open( _filename );
	if ( !in.is_open() ) FatalErr("error opening index file " << filename << " as " << _filename );


	std::string linebuf;

	size_t idx = 0;
	while ( getline(in, linebuf) ) {		// not terribly efficient, but the index will be short
	    if ( linebuf.size() > 0 ) {			// allow blank lines which probably shouldn't exist
		line_data tmp;

		// should be fasta comment up front until first tab
		size_t i = linebuf.find_first_of(" \t");	// stop at space or tab to pull out just the contig name
		if ( i == string::npos ) FatalErr( "non-emtpy line with no tab or space in FASTA index " << filename);

		// chromosome name
		std::string chr = linebuf.substr(0, i);

		size_t j = i;

		// comment, if there
		if ( linebuf[i] == ' ' ) {			// if we didn't hit a tab, then there's a fasta comment
		    j = linebuf.find_first_of('\t',i);		// find the tab, starting at i (avoids checking i+1 length)
		    if ( j == string::npos ) FatalErr( "non-emtpy line with no 2nd tab in FASTA index " << filename);
		    if ( j > i+1 ) tmp.info = linebuf.substr(i+1,j-i-1);	// ensure non-emtpy
		}

		// rest grab contig length and start and ignore the line length stuff at the end
		if ( j+1 >= linebuf.size() ) FatalErr("missing part of line in FASTA index " << filename);
		istringstream s( linebuf.substr(j+1) );

		size_t tmp_uint;
		if ( (s>>tmp.len).bad() ) FatalErr("failure parsing length of contig from FASTA index " << filename);

		if ( (s>>tmp.start).bad() ) FatalErr("failure parsing length of contig from FASTA index " << filename);

//		PRINT4( chr, tmp.len, tmp.start, tmp.info );

		tmp.idx = idx++;

		_lines[chr] = tmp;
	    }
	}
    }

    const line_data& operator[]( const String& chr ) const {
	ForceAssertGt(_lines.count(chr), 0U);	// use has_key to check
	return _lines.at(chr);
    }

    Bool has_key( const String& chr ) const {
	return (_lines.count(chr) > 0);
    }

    const String& filename() const { return _filename; }
};


class fastavector : public CharVec
{
    typedef CharVec::alloc_type Alloc;

public:
    fastavector() {}

    explicit fastavector( Alloc const& alloc ) : CharVec(alloc) {}

    explicit fastavector( size_type n, char val = GeneralizedBase::X.asChar(),
                          size_type capacity = 0, Alloc const& alloc = Alloc() )
    : CharVec(n,val,capacity,alloc) {}

    template <class Itr>
    fastavector( Itr first, Itr const& last, size_type capacity = 0,
              Alloc const& alloc = Alloc() )
    : CharVec(first,last,capacity,alloc) {}

    fastavector( const String & s ) { assign(s.begin(),s.end()); }

    explicit fastavector( const bvec& b ) { SetFrom(b); }

    // compiler-supplied copying and destructor are OK

    fastavector& resize( size_type n, char val = GeneralizedBase::X.asChar() )
    { CharVec::resize(n,val); return *this; }

    void SetFrom( const bvec& b )
    { clear().reserve(b.size());
    bvec::const_iterator end(b.end());
    for ( bvec::const_iterator itr(b.begin()); itr != end; ++itr )
       push_back(Base::val2Char(*itr)); }

    String ToString() const
    { String s; s.Set(&(*this)[0],size()); return s; }

    /// SetToSubOf(fv, start, len):  Set this to the length len piece of fv,
    /// starting at position start.  The case where this == &fv is allowed.
    /// If len == -1, go all the way to the end of fv.
    fastavector& SetToSubOf( fastavector const& fv, size_type start,
                             size_type len )
    { AssertLe( start, fv.size() );
      if ( len == ~0U ) len = fv.size() - start;
      AssertLe( len, fv.size()-start );
      reserve(len);
      memcpy(data(),fv.data()+start,len);
      setSize(len);
      return *this; }

    bool HasGaps() const
    { const_iterator stop(end());
      for ( const_iterator itr(begin()); itr != stop; ++itr )
          if ( *itr == 'n' ) return true;
      return false; }

    /// Split this into chunks, separated by gaps ('n'), and return each chunk
    /// as a gapless fastavector.
    vec<fastavector> SplitOnGaps() const;

    /// Returns a basevector, and sets a bitvector of ambiguous bases, if supplied
    /// possibly represent.  Fails if this fastavector has any gaps.
    /// Returns an empty set if the number of possibilities is more than max.
    basevector ToBasevector( bitvector* ambiguous = NULL) const;


    // same as ToBasevector, but randomizes ambiguous bases
    basevector ToBasevectorRandom() const;

    /// Returns the set of all basevectors that could this fastavector could
    /// possibly represent.  Fails if this fastavector has any gaps.
    /// Returns an empty set if the number of possibilities is more than max.
    vecbasevector AllBasevectors( size_t max = 1000 ) const;

    /// Returns the set of all kmers contained in AllBasevectors (but if a kmer
    /// location has more than small_max possibilities, ignore it.)
    /// Much more likely to return a reasonably sized set than AllBasevectors.
    vecbasevector AllKmers( size_type J, size_t small_max = 100 ) const;

    /// Prints in a fasta format: "><string_id>\n" followed by the full base
    /// sequence stored in the basevector; breaks
    /// long sequences nicely into 80-character lines
    void Print( ostream& out, const String& id = "" ) const;

    void ReverseComplement();

    void Append( const fastavector& f )
    { append(f.begin(),f.end()); }

    int AmbCount() const
    { int result = 0;
      const_iterator stop(end());
      for ( const_iterator itr(begin()); itr != stop; ++itr )
          if ( GeneralizedBase::fromChar(*itr).isAmbiguous() ) result += 1;
      return result; }

  int MaxWindowAmbCount(const int sz_window) const 
  {
    int max_counts = 0;
    const int nb = size();
    if (nb < sz_window) {
      for (int ib = 0; ib < nb; ib++)
        if (GeneralizedBase::fromChar((*this)[ib]).isAmbiguous())
          max_counts++;
    }
    else {
      int counts = 0;
      for (int ib = 0; ib < sz_window; ib++)
        if (GeneralizedBase::fromChar((*this)[ib]).isAmbiguous())
          counts++;
    
      max_counts = counts;

      for (int ib = sz_window; ib < nb; ib++) {
        if (GeneralizedBase::fromChar((*this)[ib            ]).isAmbiguous()) counts++;
        if (GeneralizedBase::fromChar((*this)[ib - sz_window]).isAmbiguous()) counts--;
        if (counts > max_counts)
          max_counts = counts;
      }
    }
    return max_counts;
  }
    



    fastavector& combine( fastavector const& fv )
    { if ( fv.size() > size() ) resize(fv.size());
      iterator dst(begin());
      const_iterator stop(fv.end());
      for ( const_iterator itr(fv.begin()); itr != stop; ++itr, ++dst )
          *dst = GeneralizedBase::ambiguityCode(*itr,*dst);
      return *this; }

    /// Compute the concatenation of two fastavectors
    friend fastavector Cat( const fastavector& left, const fastavector& right )
    { fastavector result(left.begin(),left.end(),left.size()+right.size());
      result.append(right.begin(),right.end());
      return result; }

    friend fastavector Cat( const fastavector& v1, const fastavector& v2,
         const fastavector& v3 )
    {    fastavector v( v1.begin( ), v1.end( ), v1.size() + v2.size() + v3.size() );
         v.append( v2.begin( ), v2.end( ) );
         v.append( v3.begin( ), v3.end( ) );
         return v;    }

    /// Combine: merge fastavectors starting at their beginnings.
    friend fastavector Combine( const vec<fastavector>& F )
    { fastavector result;
      vec<fastavector>::const_iterator end(F.end());
      for ( vec<fastavector>::const_iterator itr(F.begin()); itr != end; ++itr )
          result.combine(*itr);
      return result; }

     friend fastavector Combine( const fastavector& f1, const fastavector& f2 )
     { fastavector result(f1); result.combine(f2); return result; }

     // Load from fasta file.

     friend void LoadFromFastaFile( const String& f, vec<fastavector>& v );
  
     friend void LoadFromFastaFile( const String& f, vec<fastavector>& v, 
				    vec<String>& names );
     

     friend void LoadRegionsFromIndexedFastaFile( const String& fasta, const String& index,
	     const vec<FastaRegion>& regions,
	     vec<fastavector>& v,
	     int pad=0,
	     vec<std::pair<int,int>>* padded=nullptr
	     );

     friend void LoadRegionFromIndexedFastaFile( const String& fasta, const fastaindex& fai,
	    const String& chr, const int start, const int end,
	     fastavector& fv
	     );

     friend void LoadRegionFromIndexedFastaFile( const String& fasta, const String& index,
	    const String& chr, const int start, const int end,
	     fastavector& fv
	     );

     // Write a fasta-format file which contains the bases in the
     // input fasta scaffolded together (with gaps, etc.) as defined
     // by the scaffolds.  If supplied, rc defines which contigs are
     // reverse-complement. Gaps < min_gap will be reset at min_gap.
     friend void WriteScaffoldedFasta( const String &out_file,
				       const vec<fastavector> &fasta,
				       const vec<superb> &scaffolds,
				       const vec<Bool> &rc = vec<Bool>( ),
				       const int min_gap = 1, 
				       const char gap_char = 'N',
				       const Bool ncbi_format = False);

     // Remove contigs that don't appear in a scaffold.  

     friend void RenumberAndMinimize( vec<fastavector>& f, vec<superb>& s );
  
};


SELF_SERIALIZABLE(fastavector);
extern template class digraphE<fastavector>;

typedef MasterVec<fastavector> vecfastavector;
typedef vecfastavector vecfvec;
extern template class OuterVec<fastavector>;

typedef fastavector fvec;
typedef fastavector FastaVec;
typedef vecfastavector FastaVecVec;


#endif
