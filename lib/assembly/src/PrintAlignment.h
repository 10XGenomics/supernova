///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef PRINTALIGNMENT
#define PRINTALIGNMENT

#include "Alignment.h"
#include "Basevector.h"
#include "CoreTools.h"
#include "Qualvector.h"

void PrintBlanks( ostream& out, int n );

template<class BASEVEC>
void PrintBases( ostream& out, const BASEVEC& rd, int from, int to );

template<class BASEVEC1, class BASEVEC2>
void PrintVisualAlignment( Bool abbreviate, ostream& out, const BASEVEC1& rd1, 
     const BASEVEC2& rd2, const align& a, 
     const qualvector& scores1 = qualvector(0), 
     const qualvector& scores2 = qualvector(0), 
     int begin = 0, Bool one_frame = false, int min_score_to_abbrev = 0,
     Bool abbeviate_poor = False, float min_fract_poor = 2.0,
     Bool abbreviate_good = False, float max_fract_good = 0.05,
     Bool print_heads_and_tails = True, const Bool CtoT_special = False,
     const int pw = 80 );

// This version deletes spurious blank lines.  Ugly as sin.

template<class BASEVEC1, class BASEVEC2>
void PrintVisualAlignmentClean( Bool abbreviate, ostream& out, const BASEVEC1& rd1, 
     const BASEVEC2& rd2, const align& a, 
     const qualvector& scores1 = qualvector(0), 
     const qualvector& scores2 = qualvector(0), 
     int begin = 0, Bool one_frame = false, int min_score_to_abbrev = 0,
     Bool abbeviate_poor = False, float min_fract_poor = 2.0,
     Bool abbreviate_good = False, float max_fract_good = 0.05,
     Bool print_heads_and_tails = True, const Bool CtoT_special = False,
     const int pw = 80 )
{
     ostringstream outx;
     PrintVisualAlignment( abbreviate, outx, rd1, rd2, a, scores1, scores2,
          begin, one_frame, min_score_to_abbrev, abbeviate_poor, min_fract_poor,
          abbreviate_good, max_fract_good, print_heads_and_tails, CtoT_special, pw );
     istringstream in( outx.str( ) );
     String line;
     vec<String> lines;
     while(1)
     {    getline( in, line );
          if ( in.fail( ) ) break;
          lines.push_back(line);    }
     if ( lines.nonempty( ) && lines[0].Contains( "(p", 0 ) )
     {    out << "\n" << lines[0] << "\n\n";
          return;    }
     int start = 0;
     if ( lines.size( ) >= 2 && WhiteSpaceFree(lines[0]).empty( ) 
          && WhiteSpaceFree(lines[1]).empty( ) ) 
     {    start = 1;    }
     ostringstream out1;
     for ( int i = start; i < lines.isize( ); i++ )
     {    out1 << lines[i] << "\n";
          int j;
          for ( j = i + 1; j < lines.isize( ); j++ )
               if ( WhiteSpaceFree( lines[j] ).size( ) > 0 ) break;
          for ( int k = 0; k < Min( 1, j - i - 1 ); k++ )
               out1 << "\n";
          i = j - 1;    }
     istringstream in2( out1.str( ) );
     lines.clear( );
     while(1)
     {    getline( in2, line );
          if ( in2.fail( ) ) break;
          lines.push_back(line);    }
     for ( int i = 0; i < lines.isize( ); i++ )
     {    if ( i > 0 && i < lines.isize( ) - 1 && WhiteSpaceFree(lines[i]).empty( ) )
          {    int n = lines[i-1].size( );
               int b;
               for ( b = 0; b < lines[i+1].isize( ); b++ )
                    if ( lines[i+1][b] != ' ' ) break;
               if ( b - n >= 10 ) continue;    }
          out << lines[i] << "\n";    }    }

template<class BASEVEC1, class BASEVEC2>
void PrintVisualAlignment( Bool rd2_is_rc, Bool abbreviate, ostream& out, 
     const BASEVEC1& rd1, BASEVEC2 rd2, const align& a, 
     const qualvector& scores1 = qualvector(0), qualvector scores2 = qualvector(0),
     int begin = 0, Bool one_frame = false, int min_score_to_abbrev = 0,
     Bool abbeviate_poor = False, float min_fract_poor = 2.0,
     Bool abbreviate_good = False, float max_fract_good = 0.05, const int pw = 80 );

template<class BASEVEC1, class BASEVEC2>
inline void PrintVisualAlignment( Bool abbreviate, ostream& out, 
     const BASEVEC1& rd1, const BASEVEC2& rd2, const alignment& a, 
     const qualvector& scores1 = qualvector(0), 
     const qualvector& scores2 = qualvector(0),
     int begin = 0, Bool one_frame = false, int min_score_to_abbrev = 0,
     Bool abbreviate_poor = False, float min_fract_poor = 2.0, const int pw = 80 )
{    PrintVisualAlignment( abbreviate, out, rd1, rd2, align(packalign(a)),
          scores1, scores2, begin, one_frame, min_score_to_abbrev,
          abbreviate_poor, min_fract_poor, pw );    }

template<class BASEVEC1, class BASEVEC2>
inline void PrintVisualAlignment( Bool rd2_is_rc, Bool abbreviate, ostream& out, 
     const BASEVEC1& rd1, BASEVEC2 rd2, const alignment& a, 
     const qualvector& scores1 = qualvector(0),
     qualvector scores2 = qualvector(0),
     int begin = 0, Bool one_frame = false, int min_score_to_abbrev = 0,
     Bool abbreviate_poor = False, float min_fract_poor = 2.0, const int pw = 80 )
{    PrintVisualAlignment( rd2_is_rc, abbreviate, out, rd1, rd2, 
          align(packalign(a)), scores1, scores2, begin, one_frame,
          min_score_to_abbrev, abbreviate_poor, min_fract_poor, pw );    }

void PrintAlignment( ostream& out, const basevector& rd1, 
     const basevector& rd2, const alignment& a );
void PrintAlignment( Bool rd2_is_rc, ostream& out, const basevector& rd1, 
     basevector rd2, const alignment& a );
void PrintReadWithScores( basevector& B, qualvector& Q, ostream& out,
     const int pw = 80 );
void PrintErrorsInAlignment(ostream& out, const basevector& rd1, 
     const basevector& rd2, const alignment& a, const qualvector& scores1,
     const qualvector& scores2 );

#endif
