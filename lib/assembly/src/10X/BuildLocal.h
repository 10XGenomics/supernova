// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

// Build a local assembly, trying to walk from the end of one line to another.

#ifndef TENX_BUILD_LOCALX_H
#define TENX_BUILD_LOCALX_H

#include "CoreTools.h"
#include "graph/Digraph.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "10X/paths/ReadPathVecX.h"

template<class VE> void GetBarcodes( 

     // inputs:
     const int s1, const int GRAB, const int MAX_BARCODES, const int MIN_KMERS, 
     const HyperBasevectorX& hb, const digraphE<vec<int>>& D, 
     const vec<int>& to_left, const vec<int>& lens, const vec<int>& mult, VE& ebcx, 
     
     // output:
     vec<int>& b )
{
     int total = 0, sc = s1;
     while(1)
     {    if ( D.O(sc)[0] < 0 ) break;
          total += lens[sc];
          for ( int j = 0; j < D.O(sc).isize( ); j++ )
          {    int e = D.O(sc)[j];
               if ( mult[e] != 1 ) continue;
               if ( hb.Kmers(e) < MIN_KMERS ) continue;
               if ( (int) ebcx[e].size( ) > MAX_BARCODES ) continue;
               for ( int l = 0; l < (int) ebcx[e].size( ); l++ )
                    b.push_back( ebcx[e][l] );    }
          if ( total >= GRAB ) break;
          int w = to_left[sc];

          // Skip over gaps.
          
          if ( D.To(w).solo( ) && D.From(w).solo( ) && D.OTo( w, 0 )[0] < 0 )
          {    int x = D.To(w)[0];
               if ( D.From(x).solo( ) && D.To(x).solo( ) )
               {    int sc = D.ITo(x,0);
                    continue;    }    }

          if ( D.To(w).size( ) != 2 || D.From(w).size( ) != 1 ) break;
          if ( D.To(w)[0] != D.To(w)[1] ) break;
          int v = D.To(w)[0];
          if ( D.From(v).size( ) != 2 || D.To(v).size( ) != 1 ) break;
          int d1 = D.ITo(w,0), d2 = D.ITo(w,1);
          for ( int pass = 1; pass <= 2; pass++ )
          {    int d = ( pass == 1 ? d1 : d2 );
               if ( D.O(d)[0] < 0 ) continue;
               for ( int j = 0; j < D.O(d).isize( ); j++ )
               {    int e = D.O(d)[j];
                    if ( mult[e] != 1 ) continue;
                    if ( hb.Kmers(e) < MIN_KMERS ) continue;
                    if ( (int) ebcx[e].size( ) > MAX_BARCODES ) continue;
                    for ( int l = 0; l < (int) ebcx[e].size( ); l++ )
                         b.push_back( ebcx[e][l] );    }    }
          total += lens[d1] + lens[d2];
          if ( total >= GRAB ) break;
          sc = D.ITo(v,0);    }
     UniqueSort(b);    }

// Differences of GetBarcodes2 with GetBarcodes:
// - uses uni filter instead of mult (don't really know that this matters);
// - different filtering for walking over bubbles.

template<class VE> void GetBarcodes2( 

     // inputs:
     const int s1, const int GRAB, const int MAX_BARCODES, const int MIN_KMERS, 
     const HyperBasevectorX& hb, const digraphE<vec<int>>& D, 
     const vec<int>& to_left, const vec<int>& lens, const vec<Bool>& uni, VE& ebcx, 
     
     // output:
     vec<int>& b )
{
     int total = 0;
     vec<int> ss = {s1};
     while(1)
     {    int tot = 0;
          for ( int m = 0; m < ss.isize( ); m++ )
          {    int s = ss[m];
               if ( D.O(s)[0] < 0 ) continue;
               tot += lens[s];
               for ( int j = 0; j < D.O(s).isize( ); j++ )
               {    int e = D.O(s)[j];
                    if ( !uni[e] ) continue;
                    if ( hb.Kmers(e) < MIN_KMERS ) continue;
                    if ( (int) ebcx[e].size( ) > MAX_BARCODES ) continue;
                    for ( int l = 0; l < (int) ebcx[e].size( ); l++ )
                         b.push_back( ebcx[e][l] );    }    }
          total += tot / ss.size( );
          if ( total >= GRAB ) break;
          int w = to_left[ ss[0] ];
          if ( D.To(w).size( ) == 0 || D.To(w).size( ) > 2 ) break;
          if ( D.To(w).size( ) == 2 && D.To(w)[0] != D.To(w)[1] ) break;
          ss.clear( );
          for ( int j = 0; j < D.To(w).isize( ); j++ )
               ss.push_back( D.ITo(w,j) );    }
     UniqueSort(b);    }

template<class VA, class VE> void BuildLocal1(

     // Edges to walk between:

     const int s1, const vec<int>& l2s,

     // Global inputs:

     const HyperBasevectorX& hb, const vec<int>& inv,
     const vec<int64_t>& bci, vec<Bool>& dup, vec<Bool>& bad, ReadPathVecX& paths, 
     VA& alignsb, VE& ebcx, const digraphE<vec<int>>& D, const vec<int>& dinv,
     const vec<vec<vec<vec<int>>>>& dlines, const ReadPathVec& dpaths, 
     const VecULongVec& dpaths_index, const vec<Bool>& internal,
     const vec<int>& to_left, const vec<int>& mult, const vec<int>& lens,

     // Control:

     const Bool use_rights,

     // Local outputs:

     vec<int32_t>& bcl,
     HyperBasevector& hbl,
     HyperBasevectorX& hbxl, vec<int>& kmersl, vec<int>& invl, ReadPathVec& pathsl,
     digraphE<vec<int>>& Dl, vec<int>& dinvl, 
     ReadPathVec& dpathsl, VecULongVec& dpaths_indexl,
     vec<Bool>& dupl,
     MasterVec< SerfVec<triple<int,int,int> > >& alignsbl,
     vec<vec<vec<vec<int>>>>&  dlinesl,
     String& link_report,
     digraphE<vec<int>>& Dlp,

     // Logging etc.:

     Bool verbose, const Bool results_only,
     double& c1, double& c2a, double& c2b, double& c3, double& c4, 
     double& c5, double& c6a1, double& c6a2, double& c6a3, double& c6b, 
     double& c6c, double& c7,

     // Control:

     const int MAX_READS, const Bool SINGLE = True );

void BuildLocal2(

     // Edges to walk between:

     const int s1, const int s2,

     // Global inputs:

     const digraphE<vec<int>>& D, const vec<vec<vec<vec<int>>>>& dlines,

     // Local outputs:

     HyperBasevector& hbl, vec<int>& invl, ReadPathVec& pathsl,
     digraphE<vec<int>>& Dl, vec<int>& dinvl, 
     ReadPathVec& dpathsl, VecULongVec& dpaths_indexl, vec<Bool>& dupl,
     vec<vec<vec<vec<int>>>>&  dlinesl, digraphE<vec<int>>& Dlp, vec<int>& dinvlp,

     // Status:
     // (On Dlp, match beginning of global edge s1 to position p1 on
     // local edge d1, stop of global edges s2 to position p2 on local edge d2.)

     Bool& closed, int& d1, int& p1, int& d2, int& p2,

     // Logging etc.:

     Bool verbose, const Bool results_only, const Bool SINGLE = True, 
     const Bool new_test = False );

void Surgery( const HyperBasevectorX& hb, const vec<int>& inv,
     digraphE<vec<int>>& D, vec<int>& dinv,
     vec<vec<vec<vec<int>>>>& dlines, vec< pair< int, vec<int> > >& s1s2,
     const vec< vec< triple< int, pair<int,int>, pair<int,int> > > >& closures1,
     const vec< vec< digraphE<vec<int>> > >& closures2,
     const vec< vec< vec<int> > >& closures3, const Bool require_simple,
     const Bool require_simpler, ostream& out );

template<class VP2, class VI, class VA, class VE> void Unvoid(
     const vec< pair< int, vec<int> > >& s1s2,
     const HyperBasevectorX& hb, const vec<int>& inv, const vec<int64_t>& bci,
     vec<Bool>& dup, vec<Bool>& bad, ReadPathVecX& paths, VE& ebcx, VA& alignsb,
     digraphE<vec<int>>& D, vec<int>& dinv, vec<vec<vec<vec<int>>>>& dlines,
     VP2& dpaths, VI& dpaths_index, Bool use_rights, ostream& gout,
     vec< vec< triple< int, pair<int,int>, pair<int,int> > > >& closures1,
     vec< vec< digraphE<vec<int>> > >& closures2, vec< vec< vec<int> > >& closures3,
     const String& DIR, const Bool SAVE_LOCAL, const int MAX_READS, 
     const Bool SINGLE = True, const Bool verbose = False, 
     const Bool new_test = False  );

// If there are two closures, see if one is clearly the winner.

void ChooseClosure( const HyperBasevectorX& hb, const digraphE<vec<int>>& D,
     const vec<vec<vec<vec<int>>>>& dlines, const vec<int>& tol,
     vec< triple< int, pair<int,int>, pair<int,int> > >& closures1,
     vec< digraphE<vec<int>> >& closures2, vec< vec<int> >& closures3,
     ostream& mout );

#endif
