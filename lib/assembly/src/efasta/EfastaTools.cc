///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "Bitvector.h"
#include "CoreTools.h"
#include "Fastavector.h"
#include "FastIfstream.h"
#include "Superb.h"
#include "TokenizeString.h"
#include "dna/Bases.h"
#include "efasta/EfastaTools.h"
#include "math/Functions.h"
#include "math/Array.h"
#include <iostream>

#define Err(message)                                      \
{    cout << message << endl << "\nInvalid.\n" << endl;   \
     TracebackThisProcess( );    }

efasta::efasta( const fastavector& v ) { 
  (*this).reserve( v.size() );
  for ( size_t i = 0; i < v.size(); i++ )
   (*this) += ExpandAmbCode( v[i] );
}

efasta::efasta( const basevector& v ) { 
  (*this) = v.ToString() ;
}

efasta::efasta( const vec<basevector>& x )
{    if ( x.empty( ) ) return;
     if ( x.solo( ) ) 
     {    *this = x[0].ToString( );
          return;    }

     // Find the common parts of the x's on the left and on the right, and 
     // separate those out from the variant parts.

     int left_share = 0, right_share = 0;
     for ( int r = 0; r < x[0].isize( ); r++ )
     {    Bool agree = True;
          for ( int j = 0; j < x.isize( ); j++ )
          {    if ( x[j].isize( ) <= r || x[j][r] != x[0][r] )
               {    agree = False;
                    break;    }    }
          if ( !agree ) break;
          left_share++;    }
     for ( int r = 0; r < x[0].isize( ) - left_share; r++ )
     {    Bool agree = True;
          for ( int j = 0; j < x.isize( ); j++ )
          {    if ( x[j].isize( ) - left_share <= r
                    || x[j][ x[j].isize() - r - 1 ] != x[0][ x[0].isize() - r - 1 ] )
               {    agree = False;
                    break;    }    }
          if ( !agree ) break;
          right_share++;    }
     basevector Left_share( x[0], 0, left_share );
     basevector Right_share( x[0], x[0].isize( ) - right_share, right_share );
     vec<basevector> xdiff( x.size( ) );
     for ( int i = 0; i < x.isize( ); i++ )
     {    xdiff[i].SetToSubOf( x[i], left_share,
               x[i].isize( ) - left_share - right_share );    }
     *this = Left_share.ToString( ) + "{" + xdiff[0].ToString( );
     for ( int j = 1; j < xdiff.isize( ); j++ )
          *this += "," + xdiff[j].ToString( );
     *this += "}" + Right_share.ToString( );    }

efasta::efasta( const VecEFasta& x )
{    vec<basevector> e;
     for ( size_t j = 0; j != x.size( ); j++ )
     {    vec<basevector> t;
          x[j].ExpandTo(t);
          e.append(t);    }
     *this = efasta(e);    }

efasta::efasta( const basevector& seq1, const basevector& seq2, const String& marking )
{   const int D_Score=1 , I_Score=1, S_Score=1;
    int M0 = seq1.size(), N0 = seq2.size();
    // trim off same head and tail only  if matching score is zereo
    int len_head = 0;
    int len_tail = 0;
    while ( len_head < M0 && len_head < N0 
            && seq1[ len_head ] == seq2[ len_head ] ) 
        len_head++;
    while ( M0 - len_tail > len_head && N0 - len_tail > len_head
            && seq1[ M0- len_tail -1 ] == seq2[ N0 - len_tail  -1 ] ) 
        len_tail++; 
    // match the middle parts of the two sequences
    BaseVec a( seq1, len_head, M0 - len_head - len_tail );
    BaseVec b( seq2, len_head, N0 - len_head - len_tail );
    int M = a.size(), N = b.size();
    String aa = a.ToString();
    String bb = b.ToString();
    RecArray<char> s( M+1, N+1 ); // store directions, 0:match, 1:mis, 2:del, 3:ins
    vec<int> prev_row( N+1, 0 );
    vec<int> curr_row( N+1, 0 );
    for ( int j = 0; j <= N; ++j ) {
        prev_row[j] = I_Score * j;
        s[0][j] = 3;
    }
    const int MinMatchBlockSize = 10; 
    for ( int i = 1; i <= M; i++ ) {
        curr_row[0] = D_Score * i;
        s[i][0] = 2;
        for ( int j = 1; j <= N; j++ ) {
            // skip continuous matching to speed up
            if ( aa[i-1] == bb[j-1] && s[i-1][j-1] == 0 ) {
                s[i][j] = 0;
                curr_row[j] = prev_row[j-1];
                continue;
            }
            // Other cases
            vec<char> dirs(4, vec<char>::IDENTITY);
            vec<int> scores(4); 
            const int inf= std::numeric_limits<int>::max();
            if ( aa[i-1] == bb[j-1] ) {
                int len_match_block = 1;
                for( int x = 1; x < M-i+1 && x < N-j+1; x++ )
                    if ( aa[i-1+x] == bb[j-1+x] ) len_match_block++;
                    else break;
                for( int x = 1; x < i && x < j; x++ )
                    if ( aa[i-1-x] == bb[j-1-x] ) len_match_block++;
                    else break;
                if ( len_match_block >= MinMatchBlockSize )
                    scores[0] = prev_row[j-1];
                else 
                    scores[0] = inf;
            } else {
                scores[0] = inf;
            }
            scores[1] = prev_row[j-1] + S_Score;
            scores[2] = prev_row[j] + D_Score;
            scores[3] = curr_row[j-1] + I_Score;
            SortSync(scores, dirs);
            curr_row[j] = scores.front();
            s[i][j] = dirs.front();
        }
        std::swap( prev_row, curr_row );
    }
    //// Now find the edits that generate seq b from seq a
    //// Remember to add back the trimming from front
    int i = M , j = N ;
    vec< vec< vec<char> > > efasta_reverse;
    efasta_reverse.push_back( vec<vec<char> >(1) );
    while ( i > 0 || j > 0 ) {
        const int dir = s[i][j];
        if ( dir == 0 ) {
            if ( efasta_reverse.back().size() != 1 )
                efasta_reverse.push_back( vec< vec<char> >(1) );
            efasta_reverse.back()[0].push_back( aa[i-1] );
            i--, j--;
        } 
        else {
            if ( efasta_reverse.back().size() == 1 ) {
                efasta_reverse.push_back( vec< vec<char> >(2) );
            } 
            if ( dir == 1 ) {
                efasta_reverse.back()[0].push_back( aa[i-1] );
                efasta_reverse.back()[1].push_back( bb[j-1] );
                i--, j--;
            } 
            else if ( dir == 2 ) {
                efasta_reverse.back()[0].push_back( aa[i-1] );
                i--;
            } 
            else if ( dir == 3 ) {
                efasta_reverse.back()[1].push_back( bb[j-1] );
                j--;
            } 
            else {
                cout << "Error found" << endl;
                CRD::exit(1);
            }
        }
    }
    String head = BaseVec( seq1.begin(), seq1.begin() + len_head ).ToString();
    String tail = BaseVec( seq1.end() - len_tail, seq1.end() ).ToString();
    (*this).append( head.begin(), head.end() );

    efasta head_tag( "{"+marking);
    for ( vec< vec< vec<char> > >::const_reverse_iterator it = efasta_reverse.rbegin(), end = efasta_reverse.rend();
            it != end; it++ ) {
        if ( it->size() != 1 ) {
//            (*this).push_back( '{' );
            (*this).append(head_tag);
        }
        for ( size_t i = 0; i < it->size(); ++i ) {
            if ( i != 0 )
                (*this).push_back(',');
            (*this).append( (*it)[i].rbegin(), (*it)[i].rend() );
        }
        if ( it->size() != 1 ) {
            (*this).push_back( '}' );
        }
    }
    (*this).append( tail.begin(), tail.end() );
}


String AllPlus( const vec<String>& lines )
{    String x;
     for ( size_t j = 0; j < lines.size( ); j++ )
     {    x.append( lines[j] );
          if ( j != lines.size( ) - 1 ) x.push_back( '\n' );    }
     return x;    }

void ValidationError( const String& err_msg, const vec<String>& lines,
     const int i, const int j )
{    vec<char> context;
     const int flank = 40;
     int ix = i, jx = j;
     int lcount = flank, rcount = flank;
     for ( int f = 0; f < flank; f++ )
     {    jx--;
          if ( jx < 0 )
          {    ix--;
               if ( ix < 0 ) break;
               if ( lines[ix].size( ) == 0 ) break;
               jx = lines[ix].isize( ) - 1;
               context.push_front( '\n' );    }
          context.push_front( lines[ix][jx] );    }
     context.push_front( '\n' );
     for ( int l = 0; l < 3; l++ )
          context.push_front( '.' );
     context.push_back( ' ', lines[i][j], ' ' );
     ix = i, jx = j;
     for ( int f = 0; f < flank; f++ )
     {    jx++;
          if ( jx == lines[ix].isize( ) )
          {    ix++;
               if ( ix == lines.isize( ) ) break;
               jx = 0;
               if ( lines[ix].size( ) == 0 ) break;
               context.push_back( '\n' );    }
          context.push_back( lines[ix][jx] );    }
     context.push_back( '\n', '.', '.', '.' );
     cout << "\nProblem with efasta record.\n" << err_msg << ".\n" << "context =\n"; 
     cout << "\n";
     for ( int l = 0; l < context.isize( ); l++ )
          cout << context[l];
     cout << "\n\nInvalid.\n" << endl;
     TracebackThisProcess( );    }

void ValidateEfastaRecord( const vec<String>& lines, const String& msg, 
     const Bool allow_empty_record )
{    if ( lines.empty( ) ) 
     {    if (allow_empty_record) return;
          Err( "Illegal empty efasta record." );    }
     int brackets = 0, commas = 0;
     for ( size_t i = 0; i < lines.size( ); i++ )
     {    if ( lines[i].size( ) == 0 ) 
          {    if ( msg != "" ) cout << msg << endl;
               Err( "Illegal blank line in:\n" << AllPlus(lines) );    }
          for ( size_t j = 0; j < lines[i].size( ); j++ )
          {    char c = lines[i][j];
               if ( c == 'A' || c == 'C' || c == 'G' || c == 'T' ) continue;
               else if ( c == 'N' )
               {    if ( brackets != 0 )
                    {    if ( msg != "" ) cout << msg << endl;
                         ValidationError( "Illegal use of N in choose expression",
                              lines, i, j );    }    }
               else if ( c == ',' )
               {    commas++;
                    if ( brackets != 1 )
                    {    if ( msg != "" ) cout << msg << endl;
                         ValidationError( "Illegal comma", lines, i, j );    }    }
               else if ( c == '{' ) 
               {    brackets++;
                    if ( brackets != 1 )
                    {    if ( msg != "" ) cout << msg << endl;
                         ValidationError( "Illegal bracket use",
                              lines, i, j );    }    }
               else if ( c == '}' ) 
               {    brackets--;
                    if ( brackets != 0 )
                    {    if ( msg != "" ) cout << msg << endl;
                         ValidationError( "Illegal bracket use",
                              lines, i,  j );    }
                    if ( commas == 0 )
                    {    if ( msg != "" ) cout << msg << endl;
                         ValidationError( "Illegal solo choose expression",
                              lines, i, j );    }
                    commas = 0;    }
               else 
               {    if ( msg != "" ) cout << msg << endl;
                    ValidationError( "Illegal character '" + ToString(c) + "'",
                         lines, i, j );    }    }    }
     if ( brackets != 0 )
     {    if ( msg != "" ) cout << msg << endl;
          Err( "Illegal bracket use in:\n" << AllPlus(lines) );    }    }

void LoadEfastaFlat( const String& fn, vecbasevector& bases )
{    bases.clear( );
     fast_ifstream in(fn);
     String line;
     while(1)
     {    getline( in, line );
          if ( in.fail( ) ) 
          {    if ( line.size( ) == 0 ) Err( "Illegal empty file." );
               break;    }
          if ( line.size( ) == 0 ) Err( "Illegal empty record." );
          if ( line[0] != '>' )
          {    Err( "See line = '" << line << "', which was expected "
                    << "to start with >." );    }
          vec<String> lines;
          Bool eof = False;
          while(1)
          {    char c;
               in.peek(c);
               if ( in.fail( ) ) { eof = True; break; }
               if ( c == '>' ) break;
               getline( in, line );
               lines.push_back(line);    }
          String all;
          int64_t all_size = 0;
          for ( size_t i = 0; i < lines.size( ); i++ )
               all_size += lines[i].size( );
          all.reserve(all_size);
          for ( size_t i = 0; i < lines.size( ); i++ )
               all.append( lines[i] );
          ValidateEfastaRecord(lines);
          basevector b;
          efasta(all).FlattenTo(b);
          bases.push_back(b);
          if (eof) break;    }    }

void LoadEfastaFlatGaps( const String& fn, vecbitvector& gaps )
{    gaps.clear( );
     fast_ifstream in(fn);
     String line;
     while(1)
     {    getline( in, line );
          if ( in.fail( ) ) 
          {    if ( line.size( ) == 0 ) Err( "Illegal empty file." );
               break;    }
          if ( line.size( ) == 0 ) Err( "Illegal empty record." );
          if ( line[0] != '>' )
          {    Err( "See line = '" << line << "', which was expected "
                    << "to start with >." );    }
          vec<String> lines;
          Bool eof = False;
          while(1)
          {    char c;
               in.peek(c);
               if ( in.fail( ) ) { eof = True; break; }
               if ( c == '>' ) break;
               getline( in, line );
               lines.push_back(line);    }
          if ( lines.empty( ) ) Err( "Illegal record of empty length." );
          String all;
          int64_t all_size = 0;
          for ( size_t i = 0; i < lines.size( ); i++ )
               all_size += lines[i].size( );
          all.reserve(all_size);
          for ( size_t i = 0; i < lines.size( ); i++ )
               all.append( lines[i] );
          ValidateEfastaRecord(lines);
          bitvector b;
          for ( size_t i = 0; i < all.size( ); i++ )
          {    while( i < all.size( ) && all[i] != '{' ) 
               {    if ( all[i] == 'N' ) b.push_back(1);
                    else b.push_back(0);
                    i++;    }
               if ( i < all.size( ) && all[i] == '{' ) i++;
               while ( i < all.size( ) && all[i] != ',' )
               {    b.push_back(0);
                    i++;    }
               while ( i < all.size( ) && all[i] != '}' ) i++;    }
          gaps.push_back(b);
          if (eof) break;    }    }

void LoadEfastaIntoStrings( const String& fn, VecEFasta& x, vec<String>& headers,
     const Bool keep_headers, const Bool allow_empty_record )
{    x.clear( );
     fast_ifstream in(fn);
     String line;
     while(1)
     {    getline( in, line );
          if ( in.fail( ) ) 
          {    if ( line.size( ) == 0 ) Err( "Illegal empty file." );
               break;    }
          if ( line.size( ) == 0 ) 
          {    if ( !allow_empty_record ) Err( "Illegal empty record." );
               x.push_back( efasta( ) );
               continue;    }
          if ( line[0] != '>' )
          {    Err( "See line = '" << line << "', which was expected "
                    << "to start with >." );    }
          if (keep_headers) headers.push_back(line);
          vec<String> lines;
          Bool eof = False;
          while(1)
          {    char c;
               in.peek(c);
               if ( in.fail( ) ) { eof = True; break; }
               if ( c == '>' ) break;
               getline( in, line );
               lines.push_back(line);    }
          x.resize( x.size( ) + 1 );
          efasta& all = x.back( );
          int64_t all_size = 0;
          for ( size_t i = 0; i < lines.size( ); i++ )
               all_size += lines[i].size( );
          all.resize(all_size);
          size_t pos = 0;
          for ( size_t i = 0; i < lines.size( ); i++ )
          {    memcpy( &all[pos], &lines[i][0], lines[i].size( ) );
               pos += lines[i].size( );    }
          ValidateEfastaRecord( lines, "", allow_empty_record );
          if (eof) break;    }    }

void LoadEfastaIntoStrings( const String& fn, VecEFasta& x,
     const Bool allow_empty_record )
{    vec<String> headers;
     LoadEfastaIntoStrings( fn, x, headers, False, allow_empty_record );    }

void LoadEfastaIntoStrings( const String& fn, VecEFasta& x, vec<String>& headers,
     const Bool allow_empty_record )
{    LoadEfastaIntoStrings( fn, x, headers, True, allow_empty_record );    }

int efasta::Length1( const Bool count_Ns ) const
{    int count = 0;
     for ( size_t i = 0; i < size( ); i++ )
     {    while( i < size( ) && (*this)[i] != '{' ) 
          {    if ( count_Ns || (*this)[i] != 'N' ) count++;
               i++;    }
          while ( i < size( ) && (*this)[i] != ',' )
	    {  if ( (*this)[i] != '{' ) count++;
               i++;    }
          while ( i < size( ) && (*this)[i] != '}' ) i++;    }
     return count;    }

int efasta::Index1(const int n) const
{
  int count = 0;
  for ( size_t i = 0, alt_start = 0; i < size( ); i++ ) {
    while( i < size( ) && (*this)[i] != '{' )  {
      if (count >= n) return i;
      if ( (*this)[i] != 'N' ) count++;
      i++;
    }
    alt_start = i;       
    while ( i < size( ) && (*this)[i] != ',' ) {
      if ( (*this)[i] != '{' ) count++;
      i++;
    }
    while ( i < size( ) && (*this)[i] != '}' ) i++;
    if (count > n) return alt_start;
  }
  return size();
}

int efasta::Index1Alt( const int n ) const
{    int count = 0;
     for ( size_t i = 0; i < size( ); i++ ) 
     {    while( i < size( ) && (*this)[i] != '{' )  
          {    if ( count++ == n ) return i;
               i++;    }
          while ( i < size( ) && (*this)[i] != ',' ) 
          {    if ( (*this)[i] != '{' && (*this)[i] != ',' ) 
               {    if ( count++ == n ) return i;    };
               i++;    }
          while ( i < size( ) && (*this)[i] != '}' ) i++;    }
     return size( );    }


int efasta::MinLength( ) const {
  Bool amb = False;
  vec<int> aLocLens;
  int count = 0;
  for ( size_t i = 0; i < size(); i++ ){
    if ( (*this)[i] == '{' ){
      amb = True;
      aLocLens.resize(1);
      aLocLens.back() = 0;
    }else if ( (*this)[i] == '}' ){
      amb = False;
      count += Min( aLocLens );
    }else if ( (*this)[i] == ',' ){
      aLocLens.push_back(0);
    }else{
      if ( amb )
	aLocLens.back()++;
      else
	count++;
    }
  }
  return count;
}

int efasta::MaxLength( ) const {
  Bool amb = False;
  vec<int> aLocLens;
  int count = 0;
  for ( size_t i = 0; i < size(); i++ ){
    if ( (*this)[i] == '{' ){
      amb = True;
      aLocLens.resize(1);
      aLocLens.back() = 0;
    }else if ( (*this)[i] == '}' ){
      amb = False;
      count += Max( aLocLens );
    }else if ( (*this)[i] == ',' ){
      aLocLens.push_back(0);
    }else{
      if ( amb )
	aLocLens.back()++;
      else
	count++;
    }
  }
  return count;
}

void SplitEfastaIntoContigs( const VecEFasta& scaffolds,
     VecEFasta& contigs, vec<superb>& scaffold_structure )
{    contigs.clear( ), scaffold_structure.clear( );
     for ( size_t i = 0; i != scaffolds.size( ); i++ )
     {    const String& s = scaffolds[i];
          vec<int> contig_ids, gaps;
          ForceAssert( s.size( ) > 0 );
          for ( size_t j = 0; j < s.size( ); j++ )
          {    if ( s[j] == 'N' )
               {    size_t j_start = j;
                    while ( j < s.size( ) && s[j] == 'N' ) j++;
                    gaps.push_back( j_start - j );    }
               else
               {    String c;
                    while ( j < s.size( ) && s[j] != 'N' ) c.push_back( s[j++] );
                    contig_ids.push_back( contigs.size( ) );
                    contigs.push_back(c);    }
               j--;    }
          if ( gaps.size( ) >= contig_ids.size( ) )
          {    cout << "SplitEfastaIntoContigs: leading and trailing gaps are not "
                         << "allowed." << endl << "Abort." << endl;
               exit(1);    }
          superb ss;
          ss.SetNtigs( contig_ids.size( ) );
          for ( size_t j = 0; j < contig_ids.size( ); j++ )
          {    ss.SetTig( j, contig_ids[j] );
               ss.SetLen( j, contigs[ contig_ids[j] ].Length1( ) );
               if ( j < contig_ids.size( ) - 1 ) 
               {    ss.SetGap( j, gaps[j] );
                    ss.SetDev( j, 0 );    }    }
          scaffold_structure.push_back(ss);    }    }

void efasta::FlattenTo( basevector& b ) const
{    b.clear( );
     for ( size_t i = 0; i < size( ); i++ )
     {    while( i < size( ) && (*this)[i] != '{' ) 
          {    if ( (*this)[i] == 'N' ) b.push_back(0);
               else b.push_back( as_char( (*this)[i] ) );
               i++;    }
          i++;
          while ( i < size( ) && (*this)[i] != ',' )
          {    b.push_back( as_char( (*this)[i] ) );
               i++;    }
          while ( i < size( ) && (*this)[i] != '}' ) i++;    }    }


void efasta::FlattenTo( basevector& b, bitvector& gaps ) const
{    FlattenTo(b);
     gaps.resize( b.size( ), False );
     int count = 0;
     for ( size_t i = 0; i < size( ); i++ )
     {    while( i < size( ) && (*this)[i] != '{' ) 
          {    gaps.Set( count++, (*this)[i] == 'N' );
               i++;    }
          i++;
          while ( i < size( ) && (*this)[i] != ',' )
          {    gaps.Set( count++, False );
               i++;    }
          while ( i < size( ) && (*this)[i] != '}' ) i++;    }    }



char Amb(vec<String> p)
{
  Sort(p);
  if (p.size() == 2) {
    if      (p[0] == "A" && p[1] == "C") return 'M';
    else if (p[0] == "A" && p[1] == "G") return 'R';
    else if (p[0] == "A" && p[1] == "T") return 'W';
    else if (p[0] == "C" && p[1] == "G") return 'S';
    else if (p[0] == "C" && p[1] == "T") return 'Y';
    else if (p[0] == "G" && p[1] == "T") return 'K';
  }
  else if (p.size() == 3) {
    if      (p[0] == "A" && p[1] == "C" && p[2] == "G") return 'V';
    else if (p[0] == "A" && p[1] == "C" && p[2] == "T") return 'H';
    else if (p[0] == "A" && p[1] == "G" && p[2] == "T") return 'D';    
    else if (p[0] == "C" && p[1] == "G" && p[2] == "T") return 'B';    
  }
  else if (p.size() == 4 && p[0] == "A" && p[1] == "C" && p[2] == "G" && p[3] == "T") {
    return 'N';
  }
  return 0;
}

int efasta::Ambiguities( ) const
{    int comma = 0;
     for ( size_t i = 0; i < size( ); i++ )
          if ( (*this)[i] == ',' ) ++comma;
     return comma;    }

int efasta::AmbEventCount( ) const
{    int bra = 0;
     for ( size_t i = 0; i < size( ); i++ )
          if ( (*this)[i] == '{' ) ++bra;
     return bra;    }

int efasta::AmbCount( ) const
{    int count = 0;
     for ( size_t i = 0; i < size( ); i++ )
     {    while( i < size( ) && (*this)[i] != '{' ) i++;
          i++;
          int c = 0, mc = 0;
          while( i < size( ) && (*this)[i] != '}' )
          {    if ( (*this)[i] == ',' ) 
               {    mc = Max( mc, c );
                    c = 0;    }
               else c++;
               i++;    }
          mc = Max( mc, c );
          count += mc;    }
     return count;    }


Bool efasta::IndexInAmbiguity( const longlong pos ) const {
  ForceAssertLt( pos, size() );
  ForceAssertGe( pos, 0 );
  if ( (*this)[pos] == '}' || (*this)[pos] == '{' || (*this)[pos] == ',' )
    return True;
  if ( (*this).find("}",pos+1) < (*this).find("{",pos+1) )
    return True;
  return False;
}


// Return total ambiguous bases, plus count snp and indel events
int efasta::AmbCount(int& snp_count, int& indel_count) const
{
  int total_bases = 0;
  snp_count = indel_count = 0;

  vec<char> seps;
  seps.push_back(',');

  for (size_t i = 0; i < size(); i++) {
    while (i < size() && (*this)[i] != '{') {
      if ((*this)[i] == 'N') ++total_bases;
      ++i;
    }
    if (i < size()) {

      size_t j = i + 1;
      while (j < size() && (*this)[j] != '}') 
        j++;
      
      String stuff;
      for (size_t l = i + 1; l < j; l++)
        stuff.push_back((*this)[l]);

      vec<String> parts;
      TokenizeStrictly(stuff, seps, parts);
      char amb = Amb(parts);


      if (amb != 0) {
	++snp_count;
	++total_bases;
      } else {
	int maxc = 0, minc = isize();
        for (size_t ip = 0; ip < parts.size(); ip++) {
	  int c = parts[ip].isize();
	  maxc = Max(maxc, c);
	  minc = Min(minc, c);
	}
	total_bases += maxc;
	//if (maxc == minc) snp_count += maxc;
	//else ++indel_count;
	++indel_count;
      }
    }
  }
  return total_bases;
}

void efasta::FlattenTo(fastavector & f, char gap_charactor) const
{
  f.clear();
  f.reserve(size());
  vec<char> seps;
  seps.push_back(',');
  for (size_t i = 0; i < size(); i++) {
    while(i < size() && (*this)[i] != '{') {
      if ((*this)[i] == 'N') f.push_back(gap_charactor);
      else f.push_back((*this)[i]);
      i++;
    }
    if (i < size()) {

      size_t j = i + 1;
      while (j < size() && (*this)[j] != '}') 
        j++;
      
      String stuff;
      for (size_t l = i + 1; l < j; l++)
        stuff.push_back((*this)[l]);
      
      vec<String> parts;
      TokenizeStrictly(stuff, seps, parts);
      char amb = Amb(parts);
      
      if (amb != 0) 
        f.push_back(amb);
      else {
        for (size_t l = 0; l < parts[0].size(); l++)
          f.push_back(parts[0][l]);
      }
      i = j;
    }
  }
}





void efasta::FlattenTo(fastavector & f, vec<Ambiguity> & va,
		       const Bool ambiguous_base_codes) const
{
  f.clear();
  f.reserve(size());
  va.clear();
  vec<char> seps;
  seps.push_back(',');
  for (size_t i = 0; i < size(); i++) {
    while (i < size() && (*this)[i] != '{') {
      if ((*this)[i] == 'N') f.push_back('n');
      else                   f.push_back((*this)[i]);
      i++;
    }
    if (i < size()) {

      size_t j = i + 1;
      while (j < size() && (*this)[j] != '}') 
        j++;
      
      String stuff;
      for (size_t l = i + 1; l < j; l++)
        stuff.push_back((*this)[l]);

      vec<String> parts;
      TokenizeStrictly(stuff, seps, parts);
      char amb = Amb(parts);

      if (ambiguous_base_codes && amb != 0) {
        f.push_back(amb);
      }
      else {
        for (size_t ip = 1; ip < parts.size(); ip++) {
          Ambiguity a;
          a.start   = f.size();
          a.size    = parts[0].size();
          a.replace = parts[ip];     
          va.push_back(a);
        }

        for (size_t l = 0; l < parts[0].size(); l++)
          f.push_back(parts[0][l]);

      }

      i = j;
    }
  }
}


void efasta::FlattenNMaxTo(fastavector & f) const
{
  f.clear();
  f.reserve(size());
  vec<char> seps;
  seps.push_back(',');
  for (size_t i = 0; i < size(); i++) {
    while(i < size() && (*this)[i] != '{') {
      if ((*this)[i] == 'N') f.push_back('n');
      else f.push_back((*this)[i]);
      i++;
    }
    if (i < size()) {

      size_t j = i + 1;
      while (j < size() && (*this)[j] != '}') 
        j++;
      
      String stuff;
      for (size_t l = i + 1; l < j; l++)
        stuff.push_back((*this)[l]);
      
      vec<String> parts;
      TokenizeStrictly(stuff, seps, parts);
     
      int im = parts.size() -1;
      int vm = parts[im].size();
      for ( int k = 0; k < parts.isize() -1; k++ )
	if ( parts[k].isize() > vm ){
	  im = k;
	  vm = parts[k].size();
	}
      
      for (size_t l = 0; l < parts[im].size(); l++)
	f.push_back('N');
      
      i = j;
    }
  }
}

void efasta::FlattenMaxTo(fastavector & f) const
{
  f.clear();
  f.reserve(size());
  vec<char> seps;
  seps.push_back(',');
  for (size_t i = 0; i < size(); i++) {
    while(i < size() && (*this)[i] != '{') {
      if ((*this)[i] == 'N') f.push_back('n');
      else f.push_back((*this)[i]);
      i++;
    }
    if (i < size()) {

      size_t j = i + 1;
      while (j < size() && (*this)[j] != '}') 
        j++;
      
      String stuff;
      for (size_t l = i + 1; l < j; l++)
        stuff.push_back((*this)[l]);
      
      vec<String> parts;
      TokenizeStrictly(stuff, seps, parts);
      Sort( parts );
      int im = 0;
      int vm = parts[im].size();
      for ( int k = 0; k < parts.isize(); k++ )
	if ( parts[k].isize() > vm ){
	  im = k;
	  vm = parts[k].size();
	}
      
      for (size_t l = 0; l < parts[im].size(); l++)
	f.push_back( parts[im][l] );
      
      i = j;
    }
  }
}

void efasta::FlattenMinTo(fastavector & f) const
{
  f.clear();
  f.reserve(size());
  vec<char> seps;
  seps.push_back(',');
  for (size_t i = 0; i < size(); i++) {
    while(i < size() && (*this)[i] != '{') {
      if ((*this)[i] == 'N') f.push_back('n');
      else f.push_back((*this)[i]);
      i++;
    }
    if (i < size()) {

      size_t j = i + 1;
      while (j < size() && (*this)[j] != '}') 
        j++;
      
      String stuff;
      for (size_t l = i + 1; l < j; l++)
        stuff.push_back((*this)[l]);
      
      vec<String> parts;
      TokenizeStrictly(stuff, seps, parts);
      Sort( parts );
      int im = 0;
      int vm = parts[im].size();
      for ( int k = 0; k < parts.isize(); k++ )
	if ( parts[k].isize() < vm ){
	  im = k;
	  vm = parts[k].size();
	}
      
      for (size_t l = 0; l < parts[im].size(); l++)
	f.push_back( parts[im][l] );
      
      i = j;
    }
  }
}

void efasta::ReverseComplement()
{
    reverse(begin(),end());
    for ( char& chr : *this )
    {
        switch ( chr )
        {
        case ',': break;
        case '{': chr = '}'; break;
        case '}': chr = '{'; break;
        default: chr = GeneralizedBase::complementChar(chr); break;
        }
    }
}


void efasta::GetBlocks( vec< vec<String> >& blocks ) const
{
  blocks.clear();
  int ambCount = 0;
  vec<char> seps;
  seps.push_back(',');
  for ( String::const_iterator it = begin(); it != end(); it++) {
    blocks.push_back( vec<String>(1,"") );
    while( it != end() && *it != '{') {
      blocks.back().back() += *it;
      it++;
    }
    if ( *it == '{' ) ambCount++;
    if (it != end()) {
      String::const_iterator jt = it + 1;
      while (jt != end() && *jt != '}') jt++;
      String stuff;
      for (String::const_iterator lt = it + 1; lt < jt; lt++)
        stuff.push_back(*lt);
      
      vec<String> parts;
      TokenizeStrictly(stuff, seps, parts);
      blocks.push_back( parts );
      it = jt;
    }
    if ( it == end() ) it--;
  }
  //PRINT2( ambCount, blocks.size() );
}

void efasta::Erase1( const int start, const int stop ){

   vec< vec<String> > blocks;
  (*this).GetBlocks( blocks );
  

  int len1 = (*this).Length1();
  int pos = len1;
  for ( int ib = blocks.isize() -1; ib >= 0; ib-- ){
    for ( int ip = blocks[ib][0].isize() -1; ip >= 0; ip-- ){
      pos--;
      if ( pos < stop && pos >= start ){
	if ( pos - ip >= start ){
	  blocks[ib][0].erase( blocks[ib][0].begin(), blocks[ib][0].begin() + ip + 1 );
	  pos -= ip;
	  ip  -= ip;
	}else{
	  blocks[ib][0].erase( blocks[ib][0].begin() + ip - (pos - start), blocks[ib][0].begin() + ip + 1 );
	  pos -= pos - start;
	  ip  -= pos - start;
	}
      }
    }
  }
  
  efasta result;
  for ( int bi = 0; bi < blocks.isize(); bi++ ){
    if ( blocks[bi].size() > 1 ){
      result += "{";
      for ( int ii = 0; ii < blocks[bi].isize(); ii++ ){
	result += blocks[bi][ii];
	if ( ii < blocks[bi].isize() -1 )
	  result += ",";
	else result += "}";
      }
    }else
      result  += blocks[bi][0];
  }
  
  *this = result;
  return;
}

void efasta::Extract1( const int start, const int stop, efasta& result ){

  if ( stop > (*this).Length1() || stop < 0 ||
       start < 0 || start >= (*this).Length1() ){
    PRINT3( start, stop, (*this).Length1() );
    FatalErr("Inconsistent start and stop positions");
  }
  
  result = (*this);
  if ( stop < (*this).isize() )
    result.Erase1( stop, result.Length1() );
  if ( start > 0 )
    result.Erase1( 0, start );
  
  ForceAssertEq( result.Length1(), stop - start );

  return;
}



void efasta::Canonicalize( efasta& canonical )
{
  canonical.clear();
  vec< vec<String> > blocks;
  GetBlocks( blocks );

  // make sure that ambiguous blocks are bracketed by non-ambiguous ones (even if empty)
  
  if ( blocks.front().size() > 1 )
    blocks.push_front( vec<String>(1,"") );
  if ( blocks.back().size() > 1 )
    blocks.push_back( vec<String>(1,"") );
  for ( int i = 0; i < blocks.isize() -1; i++ ){
    if ( blocks[i].size() > 1 && blocks[i+1].size() > 1 ){
      blocks.insert( blocks.begin() + i + 1, vec<String>(1,"") );
    }
  }

  // see if we can compute alternative representation
  
  for ( int i = 1; i < blocks.isize()-1; i++ ){
    if ( blocks[i].size() <= 1 ) continue;
    int je = -1;
    for ( int j = 0; j < blocks[i].isize(); j++ )
      if ( blocks[i][j] == "" ){
	je = j;
	break;
      }
    if ( je == -1 ) continue;
    vec<String> block2( blocks[i] );
    vec<int> lens2;
    for ( int j = 0; j < block2.isize(); j++ ){
      if ( j == je ) block2[j] = blocks[i+1][0];
     lens2.push_back( block2[j].isize() );
    }
    int max_p = Min( lens2 );
    int se = -1;
    String sshift = "";
    for ( int s = 0; s < max_p; s++ ){
      set<char> bases;
      for ( int j = 0; j < block2.isize(); j++ )
	bases.insert( block2[j][s] );
      if ( bases.size() > 1 )
	break;
      else{ 
	se = s;
	sshift += *bases.begin();
      }
    }
    
    if ( se > -1 ){
      //cout << "found room for shifting" << endl;
      //PRINT( sshift );
      blocks[i-1][0] += sshift;
      for ( int j = 0; j < blocks[i].isize(); j++ )
	if ( j != je ) blocks[i][j].erase(0, sshift.size() );
      blocks[i+1][0].erase(0,sshift.size() );
    }
 
  }

  for ( int i = 0; i < blocks.isize(); i++ ){
    if ( blocks[i].isize() == 1 ){
      canonical += blocks[i][0];
    }else if ( blocks[i].isize() > 1 ){
      canonical += '{';
      for ( int j = 0; j < blocks[i].isize(); j++ ){
	canonical += blocks[i][j];
	if ( j < blocks[i].isize() -1 )
	  canonical += ',';
      }
      canonical += '}';
    }
  }
  // if ( *this != canonical ){
//     cout << "original efasta and canonical version different:" << endl;
//     PRINT( *this );
//     PRINT( canonical );
//   }
}


// ExpandTo.  Convert to a list of fastavectors.  This will find ambiguous
// base codes.  Return False if max_count specified and the number of 
// fastavectors would exceed it.

Bool efasta::ExpandTo( vec<fastavector>& v, const int max_count ) const
{    v.clear( );
     vec< vec<fastavector> > G;
     MakeGraph(G);
     int64_t count = 1;
     for ( int l = 0; l < G.isize( ); l++ )
     {    count *= G[l].size( );
          if ( ( max_count >= 0 && count > max_count ) || count > 1000000000 ) 
               return False;    }
     v.push_back( fastavector( ) );
     for ( int l = 0; l < G.isize( ); l++ )
     {    vec<fastavector> w;
          for ( size_t r = 0; r < v.size( ); r++ )
          {    for ( int j = 0; j < G[l].isize( ); j++ )
               {    w.push_back( Cat( v[r], G[l][j] ) );    }    }
          v = w;    }
     return True;    }

Bool efasta::ExpandTo( vec<basevector>& v, const int max_count ) const
{    v.clear( );
     vec< vec<basevector> > G;
     MakeGraph(G);
     int64_t count = 1;
     for ( int l = 0; l < G.isize( ); l++ )
     {    count *= G[l].size( );
          if ( ( max_count >= 0 && count > max_count ) || count > 1000000000 ) 
               return False;    }
     v.push_back( basevector( ) );
     for ( int l = 0; l < G.isize( ); l++ )
     {    vec<basevector> w;
          for ( size_t r = 0; r < v.size( ); r++ )
          {    for ( int j = 0; j < G[l].isize( ); j++ )
               {    w.push_back( Cat( v[r], G[l][j] ) );    }    }
          v = w;    }
     return True;    }

Bool efasta::ExpandTo( vec< vec<basevector> >& G, vec< vec<int> >& paths, const int max_count ) const
{    
  G.clear();
  MakeGraph(G);
  int64_t count = 1;
  for ( int l = 0; l < G.isize( ); l++ ){    
    count *= G[l].size( );
    if ( ( max_count >= 0 && count > max_count ) || count > 1000000000 ) 
      return False;    
  }
  
  paths.clear();
  for ( int l = 0; l < G.isize( ); l++ ){    
    vec< vec<int> > w;
    if ( paths.size() > 0 ){
      for ( size_t r = 0; r < paths.size( ); r++ ){    
	for ( int j = 0; j < G[l].isize( ); j++ ){
	  w.push_back( paths[r] );
	  w.back().push_back(j);    
	}    
      }
    }else{
      for ( int j = 0; j < G[l].isize( ); j++ )
	w.push_back( vec<int>(1,j) );
    }
    paths = w;    
  }
  return True;    
}


void efasta::MakeGraph( vec< vec<fastavector> >& G ) const
{    G.clear( );
     vec<char> seps;
     seps.push_back( ',' );
     fastavector f;
     f.reserve( size( ) );
     for ( size_t i = 0; i < size( ); i++ )
     {    while( i < size( ) && (*this)[i] != '{' ) 
          {    if ( (*this)[i] == 'N' ) f.push_back( 'n' );
               else f.push_back( (*this)[i] );
               i++;    }
          vec<fastavector> g(1);
          if ( f.size( ) > 0 )
          {    g[0] = f;
               G.push_back(g);
               f.clear( );    }
          if ( i == size( ) ) break;
          size_t j;
          for ( j = i+1; j < size( ); j++ )
               if ( (*this)[j] == '}' ) break;
          String stuff;
          for ( size_t l = i + 1; l < j; l++ )
               stuff.push_back( (*this)[l] );
          vec<String> parts;
          TokenizeStrictly( stuff, seps, parts );
          char amb = Amb(parts);
          if ( amb != 0 )
          {    g[0].resize(1);
               g[0][0] = amb;    }
          else
          {    g.resize( parts.size( ) );
               for ( size_t l = 0; l < parts.size( ); l++ )
                    g[l] = parts[l];    }
          G.push_back(g);
          i = j;    }    }

void efasta::MakeGraph( vec< vec<basevector> >& G ) const
{    G.clear( );
     vec<char> seps;
     seps.push_back( ',' );
     basevector f;
     f.reserve( size( ) );
     for ( size_t i = 0; i < size( ); i++ )
     {    while( i < size( ) && (*this)[i] != '{' ) 
          {    if ( (*this)[i] == 'N' ) f.push_back(0);
               else f.push_back( as_char( (*this)[i] ) );
               i++;    }
          vec<basevector> g(1);
          if ( f.size( ) > 0 )
          {    g[0] = f;
               G.push_back(g);
               f.clear( );    }
          if ( i == size( ) ) break;
          size_t j;
          for ( j = i+1; j < size( ); j++ )
               if ( (*this)[j] == '}' ) break;
          String stuff;
          for ( size_t l = i + 1; l < j; l++ )
               stuff.push_back( (*this)[l] );
          vec<String> parts;
          TokenizeStrictly( stuff, seps, parts );
          g.clear_and_resize( parts.size( ) );
          for ( size_t l = 0; l < parts.size( ); l++ )
          {    for ( int j = 0; j < parts[l].isize( ); j++ )
                    g[l].push_back( as_char( parts[l][j] ) );    }
          G.push_back(g);
          i = j;    }    }

void efasta::MakeFromGraph( const vec< vec<basevector> >& G )
{    clear( );
     for ( int i = 0; i < G.isize( ); i++ )
     {    ForceAssertGe( G[i].isize( ), 1 );
          if ( G[i].size( ) > 1 ) push_back( '{' );
          for ( int j = 0; j < G[i][0].isize( ); j++ )
               push_back( as_base( G[i][0][j] ) );
          for ( int l = 1; l < G[i].isize( ); l++ )
          {    push_back( ',' );
               for ( int j = 0; j < G[i][l].isize( ); j++ )
                    push_back( as_base( G[i][l][j] ) );    }
          if ( G[i].size( ) > 1 ) push_back( '}' );    }    }

class gpath {

     public:

     gpath( const vec< vec<basevector> >* G, const int start ) 
          : G(G), start(start) { }

     const vec< vec<basevector> >* G;
     int start;
     vec<int> p;

     Bool Next( )
     {    if ( p.empty( ) ) return False;
          while( p.back( ) == (*G)[ start + p.isize( ) - 1 ].isize( ) - 1 )
          {    p.pop_back( );
               if ( p.empty( ) ) return False;    }
          p.back( )++;
          return True;    }

     int Length( ) const
     {    int len = 0;
          for ( int j = 0; j < p.isize( ); j++ )
               len += (*G)[ start + j ][ p[j] ].isize( );
          return len;    }

     Bool AtEnd( ) const
     {    return G->isize( ) == start + p.isize( );    }

     void AddZero( )
     {    p.push_back(0);    }

     void AddTo( basevector& b ) const
     {    for ( int j = 0; j < p.isize( ); j++ )
               b.append( (*G)[ start + j ][ p[j] ] );    }
               
};

void efasta::GetKmers( const int K, vecbasevector& kmers, const int max_per ) const
{    ForceAssertGe( K, 1 );
     vec< vec<basevector> > G;
     MakeGraph(G);
     basevector b, b0;
     for ( int i = 0; i < G.isize( ); i++ )
     {    for ( int j = 0; j < G[i].isize( ); j++ )
          {    for ( int l = 0; l < G[i][j].isize( ); l++ )
               {    if ( G[i][j].isize( ) - l >= K )
                    {    b.SetToSubOf( G[i][j], l, K );
                         kmers.push_back_reserve(b);
                         continue;    }
                    b0.SetToSubOf( G[i][j], l, G[i][j].isize( ) - l );
                    int count = 0;
                    gpath gp( &G, i + 1 );
                    while( max_per < 0 || count < max_per )
                    {    while( gp.Length( ) + b0.isize( ) < K )
                         {    if ( gp.AtEnd( ) ) break;
                              gp.AddZero( );    }
                         if ( gp.Length( ) + b0.isize( ) >= K )
                         {    b = b0;
                              gp.AddTo(b);
                              b.resize(K);
                              kmers.push_back(b);
                              count++;    }    
                         if ( !gp.Next( ) ) break;    }    }    }    }    }

String ExpandAmbCode( char x )
{    if ( x == 'N' ) return "{A,C,G,T}";
     x = toupper(x);
     if ( x == 'A' || x == 'C' || x == 'G' || x == 'T' || x == 'N' ) 
          return String(x);
     if ( x == 'K' ) return "{G,T}";
     if ( x == 'M' ) return "{A,C}";
     if ( x == 'R' ) return "{A,G}";
     if ( x == 'Y' ) return "{C,T}";
     if ( x == 'S' ) return "{C,G}";
     if ( x == 'W' ) return "{A,T}";
     if ( x == 'B' ) return "{C,G,T}";
     if ( x == 'V' ) return "{A,C,G}";
     if ( x == 'H' ) return "{A,C,T}";
     if ( x == 'D' ) return "{A,G,T}";
     cout << "Unable to expand base '" << x << "'." << endl << "Abort." << endl;
     TracebackThisProcess( );
     return String( );    }

String ExpandAmbCode( const String& x )
{    String y;
     y.reserve( x.size( ) );
     for ( int j = 0; j < x.isize( ); j++ )
     {    char c = x[j];
          if ( c == '{' || c == '}' || c == ',' ) y += c;
          else y += ExpandAmbCode(c);    }
     return y;    }

void efasta::Print( ostream& out, const String& id ) const
{
    out << '>' << id << '\n';
    unsigned sz = size();
    char const* itr = data();
    unsigned wholeLines = sz/80;
    while ( wholeLines-- )
    {
        out.write(itr,80) << '\n';
        itr += 80;
    }
    if ( (sz %= 80) )
        out.write(itr,sz) << '\n';
}

void efasta::PrintAmbiguities( ostream& out, const String& id ) const
{
  int eventCount = -1;
  vec<char> seps;
  seps.push_back(',');
  for (size_t i = 0; i < size(); i++) {
    while(i < size() && (*this)[i] != '{') {
      i++;
    }
    if (i < size()) {

      size_t j = i + 1;
      while (j < size() && (*this)[j] != '}') 
        j++;
      
      String stuff;
      for (size_t l = i + 1; l < j; l++)
        stuff.push_back((*this)[l]);
      
      eventCount++;
      vec<String> parts;
      TokenizeStrictly(stuff, seps, parts);
      for ( int pi = 0; pi < parts.isize(); pi++ ){
	out << '>' << id << "_" << eventCount << "-" << pi;
	for ( size_t si = 0; si < parts[pi].size(); ++si ){
	  if ( !(si % 80) ) out << '\n';
	  out << parts[pi][si];
	}
	out << '\n';
      }
      
      i = j;
    }
  }
}


void WriteScaffoldedEFasta( const String &out_file,
			    const VecEFasta &fasta,
			    const vec<superb> &scaffolds,
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
    // gaps. WARNING: Negative gaps reset to 1.
    vec<char> s;
    for ( int j = 0; j < S.Ntigs( ); j++ ) {
      const efasta &b = fasta[S.Tig(j)];
      for ( unsigned int l = 0; l < b.size( ); l++ )
	s.push_back( b[l] );
      if ( j < S.Ntigs( ) - 1 )
	s.push_back_copies( 'N', Max( 1, S.Gap(j) ) );
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

String ReplaceByLengths( const efasta& e )
{    String s;
     s.reserve( e.size( ) );
     for ( int i = 0; i < e.isize( ); i++ )
     {    if ( e[i] != '{' ) s += e[i];
          else
          {    s += '{';
               int count = 0;
               vec<int> counts;
               for ( i++; i < e.isize( ); i++ )
               {    if ( e[i] == ',' || e[i] == '}' )
                    {    counts.push_back(count);
                         if ( e[i] == '}' ) break;
                         count = 0;    }
                    else count++;    }
               Sort(counts);
               for ( int j = 0; j < counts.isize( ); j++ )
               {    if ( j > 0 ) s += ",";
                    s += ToString( counts[j] );    }
               s += '}';    }    }
     return s;    }

void ValidateEfastaRecord( const String& line, const String& msg,
     const Bool allow_empty_record )
{    vec<String> lines;
     lines.push_back(line);
     ValidateEfastaRecord( lines, msg, allow_empty_record );    }

void GetShares( const vec<basevector>& patches, int& left_share, int& right_share )
{    left_share = 0, right_share = 0;
     for ( int r = 0; r < patches[0].isize( ); r++ )
     {    Bool agree = True;
          for ( int j = 0; j < patches.isize( ); j++ )
          {    if ( patches[j].isize( ) <= r || patches[j][r] != patches[0][r] )
               {    agree = False;
                    break;    }    }
          if ( !agree ) break;
          left_share++;    }
     for ( int r = 0; r < patches[0].isize( ) - left_share; r++ )
     {    Bool agree = True;
          for ( int j = 0; j < patches.isize( ); j++ )
          {    if ( patches[j].isize( ) - left_share <= r
                    || patches[j][ patches[j].isize( ) - r - 1 ]
                    != patches[0][ patches[0].isize( ) - r - 1 ] )
               {    agree = False;
                    break;    }    }
          if ( !agree ) break;
          right_share++;    }    }

#include "feudal/OuterVecDefs.h"
template class OuterVec<efasta>;
