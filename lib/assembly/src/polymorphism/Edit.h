///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef EDIT_H
#define EDIT_H

#include "PackAlign.h"

// An "edit" represents a change to a read implied by its alignment to another
// read.  A "substitution" includes the case of no change.  For insertions, because 
// the base that follows the insertion has the same coordinate on the first read, we
// include it as part of the edit.

enum ETYPE { INSERTION = 1, DELETION = 2, SUBSTITUTION = 3 };

class edit {

     public:

     edit( ) { }
     edit( const ETYPE ins, const String& x, const int id ) 
          : etype(ins), n(0), seq(x), id(id)
     {    ForceAssertEq( (int) etype, (int) INSERTION );    }

     edit( const ETYPE del, const int n, const int id ) : etype(del), n(n), id(id)
     {    ForceAssertEq( (int) etype, (int) DELETION );    }

     edit( const ETYPE sub, const char b, const int id ) : etype(sub), n(0), id(id)
     {    ForceAssertEq( (int) etype, (int) SUBSTITUTION );
          seq.push_back(b);    }

     friend Bool operator<( const edit& e1, const edit& e2 )
     {    if ( e1.etype < e2.etype ) return True;
          if ( e1.etype > e2.etype ) return False;
          if ( e1.n < e2.n ) return True;
          if ( e1.n > e2.n ) return False;
          return e1.seq < e2.seq;    }

     friend Bool operator==( const edit& e1, const edit& e2 )
     {    return e1.etype == e2.etype && e1.n == e2.n && e1.seq == e2.seq;    }
     friend Bool operator!=( const edit& e1, const edit& e2 )
     {    return !( e1 == e2 );    }

     ETYPE etype;
     int n;
     String seq;
     int id;

};

// edit0: edit without id

class edit0 {

     public:

     edit0( ) { }
     edit0( const ETYPE ins, const String& x ) : etype(ins), n(0), seq(x)
     {    ForceAssertEq( (int) etype, (int) INSERTION );    }

     edit0( const ETYPE del, const int n ) : etype(del), n(n)
     {    ForceAssertEq( (int) etype, (int) DELETION );    }

     edit0( const ETYPE sub, const char b ) : etype(sub), n(0)
     {    ForceAssertEq( (int) etype, (int) SUBSTITUTION );
          seq.push_back(b);    }

     friend Bool operator<( const edit0& e1, const edit0& e2 )
     {    if ( e1.etype < e2.etype ) return True;
          if ( e1.etype > e2.etype ) return False;
          if ( e1.n < e2.n ) return True;
          if ( e1.n > e2.n ) return False;
          return e1.seq < e2.seq;    }

     friend Bool operator==( const edit0& e1, const edit0& e2 )
     {    return e1.etype == e2.etype && e1.n == e2.n && e1.seq == e2.seq;    }
     friend Bool operator!=( const edit0& e1, const edit0& e2 )
     {    return !( e1 == e2 );    }

     friend ostream& operator<<( ostream& out, const edit0& e )
     {    if ( e.etype == INSERTION ) out << "insertion of " << e.seq;
          else if ( e.etype == DELETION ) out << "deletion of " << e.n << " bases";
          else out << "substitution of " << e.seq;
          return out;    }

     int Dist( ) const // edit distance
     {    if ( etype == DELETION ) return n;
          else return seq.size( );    }

     ETYPE etype;
     int n;
     String seq;

};

vec<pair<int,edit0>> AlignToEdits(const align& a, const basevector& S,
        const basevector& T);

#endif
