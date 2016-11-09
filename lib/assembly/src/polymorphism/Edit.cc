///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "PackAlign.h"
#include "polymorphism/Edit.h"

vec<pair<int,edit0>> AlignToEdits(const align& a, const basevector& S, 
        const basevector& T) 
{
    vec<pair<int,edit0>> edits;
    int p1 = a.pos1( ), p2 = a.pos2( );
    Bool first = True;
    for ( int j = 0; j < a.Nblocks( ); j++ ) 
    {    if ( a.Gaps(j) > 0 )  
        {    if ( !( first && edits.nonempty( )
                  && p2 <= edits.back( ).first ) )
             {    edits.push( p2, edit0( DELETION, a.Gaps(j) ) );    }
             if ( first && edits.nonempty( )
                  && p2 >= edits.back( ).first )
             {    first = False;    }
             p2 += a.Gaps(j);    }
        if ( a.Gaps(j) < 0 ) 
        {    if ( !( first && edits.nonempty( )
                  && p2 <= edits.back( ).first ) )
             {    edits.push( p2, edit0( INSERTION, 
                       basevector( S, p1, -a.Gaps(j) )
                       .ToString( ) ) );    }
             if ( first && edits.nonempty( )
                  && p2 >= edits.back( ).first )
             {    first = False;    }
             p1 -= a.Gaps(j);    }
        for ( int x = 0; x < a.Lengths(j); x++ ) 
        {    if ( S[p1] != T[p2] )
             {    if ( !( first && edits.nonempty( )
                       && p2 <= edits.back( ).first ) )
                  {    edits.push( p2, edit0( SUBSTITUTION, 
                            (char) as_base(S[p1]) ) );    }
                  if ( first && edits.nonempty( )
                       && p2 >= edits.back( ).first )
                  {    first = False;    }    }
             ++p1; ++p2;    }    }    
    return edits;
}
