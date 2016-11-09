///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "CoreTools.h"
#include "PackAlign.h"
#include "paths/long/AssessBestAlignCore.h"

// PushBoundaries: for tandem-repeat indels, push left/right boundaries so that
// they extend up to the ends of the repeat.

void PushBoundaries( basevector b1, basevector b2,
     int& left1, int& right1, int& left2, int& right2 )
{    
     // Swap.

     Bool swapped = False;
     if ( left1 == right1 )
     {    swap( b1, b2 );
          swap( left1, left2 );
          swap( right1, right2 );
          swapped = True;    }

     // Try each possible repeat period.

     int n = right1 - left1;
     for ( int p = 1; p <= Min( n, left1 ); p++ )
     {    if ( n % p != 0 ) continue;
          Bool repeat = True;
          for ( int j = left1 + p; j < right1; j++ )
          {    if ( b1[j] != b1[ j + (j-left1) % p ] )
               {    repeat = False;
                    break;    }    }
          if ( !repeat ) continue;
          for ( int m = right1; m < b1.isize( ); m += p )
          {    Bool match = True;
               for ( int i = 0; i < p; i++ )
               {    if ( b1[ m + i ] != b1[ left1 + i ] )
                    {    match = False;
                         break;    }    }
               if ( !match ) break;
               if ( right2 + p > b2.isize( ) ) break;
               right1 += p;
               right2 += p;    }
          for ( int m = left1 - p; left1 >= 0; m -= p )
          {    if ( m < 0 ) break;
               Bool match = True;
               for ( int i = 0; i < p; i++ )
               {    if ( b1[ m + i ] != b1[ left1 + i ] )
                    {    match = False;
                         break;    }    }
               if ( !match ) break;
               if ( left2 - p < 0 ) break;
               left1 -= p;
               left2 -= p;    }
          break;    }

     // Swap back.

     if (swapped)
     {    swap( left1, left2 );
          swap( right1, right2 );    }    }

void DecomposeAlign( const align& a, const basevector& b1, const basevector& b2,
     vec< pair<int,int> >& P1, vec< pair<int,int> >& P2 )
{
     int n1 = b1.size( ), n2 = b2.size( );
     const int sep = 30;
     int p1 = a.pos1( ), p2 = a.pos2( );
     for ( int j = 0; j < a.Nblocks( ); j++ )
     {    Bool tied = False;
          if ( P1.nonempty( ) )
          {    if ( p1 < P1.back( ).second ) tied = True;
               if ( p2 < P2.back( ).second ) tied = True;    }
          if ( a.Gaps(j) > 0 ) 
          {    int left1 = p1, right1 = p1;
               int left2 = p2, right2 = p2 + a.Gaps(j);
               PushBoundaries( b1, b2, left1, right1, left2, right2 );
               if ( !tied )
               {    P1.push( Max( 0, left1 - sep ), Min( n1, right1 + sep ) );
                    P2.push( Max( 0, left2 - sep ), 
                         Min( n2, right2 + a.Gaps(j) + sep ) );    }
               else
               {    P1.back( ).second = Min( n1, right1 + sep );
                    P2.back( ).second = Min( n2, right2 + a.Gaps(j) + sep );    }
               p2 += a.Gaps(j);    }
          if ( a.Gaps(j) < 0 ) 
          {    int left1 = p1, right1 = p1 - a.Gaps(j);
               int left2 = p2, right2 = p2;
               PushBoundaries( b1, b2, left1, right1, left2, right2 );
               if ( !tied )
               {    P1.push( Max( 0, left1 - sep ), 
                         Min( n1, right1 - a.Gaps(j) + sep ) );
                    P2.push( Max( 0, left2 - sep ), Min( n2, right2 + sep ) );    }
               else
               {    P1.back( ).second = Min( n1, right1 + sep );
                    P2.back( ).second = Min( n2, right2 + sep );    }
               p1 -= a.Gaps(j);    }
          for ( int x = 0; x < a.Lengths(j); x++ ) 
          {    if ( b1[p1] != b2[p2] )
               {    Bool tied = False;
                    if ( P1.nonempty( ) )
                    {    if ( p1 < P1.back( ).second ) tied = True;
                         if ( p2 < P2.back( ).second ) tied = True;    }
                    if ( !tied )
                    {    P1.push( Max( 0, p1 - sep ), Min( n1, p1 + sep ) );
                         P2.push( Max( 0, p2 - sep ), 
                              Min( n2, p2 + sep ) );    }
                    else
                    {    P1.back( ).second = Min( n1, p1 + sep );
                         P2.back( ).second = Min( n2, p2 + sep );    }    }
               ++p1; ++p2;    }    }    }
