///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef LIST_FOSMIDS_H
#define LIST_FOSMIDS_H

#include "CoreTools.h"
#include "ParseSet.h"

// good_fosmids: Fosmids that assembly without a gap

inline vec<int> good_fosmids( )
{    return { 1,3,4,5,6,7,8,10,11,12,13,14,15,17,18,19,21,23,25,26,27,28,30,31,32,
          39,40,41,42,43,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,63,64,
          65,66,68,69,70,71,72,73,74,75,76,77,78,80,82,84,85,86,88,89,90,91,92,93,
          94,95,96,97,98,100,101,102,103,104,105 };    }

// ok_fosmids: Fosmids that assemble with a gap but appear to have no 
// significant gap relative to the reference

inline vec<int> ok_fosmids( )
{    return { 0,2,20,44,67,79,81 };    }

inline vec<int> AllFosmids( )
{    return {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,18,19,20,21,22,
                  23,24,25,26,27,28,29,30,31,32,34,37,39,40,41,42,43,44,45,46,47,
                  48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,68,69,
                  70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,
                  91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106};    }

inline vec<int> expand_fosmids( const String& CLASS )
{    vec<int> answer;
     vec<String> classes;
     ParseStringSet( "{" + CLASS + "}", classes );
     for ( int i = 0; i < classes.isize( ); i++ )
     {    if ( classes[i] == "goods" ) answer.append( good_fosmids( ) );
          else if ( classes[i] == "oks" ) answer.append( ok_fosmids( ) );
          else FatalErr( "Illegal CLASS specification." );    }
     UniqueSort(answer);
     return answer;    }

vec<String> AllFosmidRegions( size_t padding = 0 );
String FosmidRegion( int fid, size_t padding = 0 );
void GetRegionInfo( const String& ID, int& g, int& rstart, int& rstop,
     String& gid, String& loc );


#endif
