// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

#include "CoreTools.h"
#include "feudal/SmallVecDefs.h"
#include "feudal/OuterVecDefs.h"
#include "graph/Digraph.h"
#include "math/Functions.h"
#include "math/HoInterval.h"
#include "10X/Gap.h"
#include "10X/astats/View.h"

// View: find genomic locations of a line.  
// Returns {(chr,orientation,genomic interval,assembly interval)
// where the assembly interval is the start/stop on the line, counting gaps as
// zero length.
//
// Note:
// 1. This only looks at non-bubble parts of the line.
// 2. Designed for human genomes!
// 3. It only looks at aligns to 1,...,22,X,Y.

template<class VA> void View( const int l, const int K, const vec<int>& kmers, 
     const vec<int>& inv, const digraphE<vec<int>>& D,
     const vec<vec<vec<vec<int>>>>& dlines,
     const vec<vec< pair<int,int> >>& linelocs,
     VA&& alignsb, SerfVec< quad<int,Bool,ho_interval,ho_interval> >& view )
{    
     // Heuristics.  We might wish to make some of these into arguments.

     const int MIN_KMERS_TO_SHOW = 500;
     const int MAX_CHR = 23; // 0-based numbering; human chromosome Y
     const int MAX_GAP = 50000;

     // Get view.

     view.clear( );
     vec< quad<int,Bool,ho_interval,ho_interval> > view0;
     const vec<vec<vec<int>>>& L = dlines[l];
     int pos = 0;
     for ( int j = 0; j < L.isize( ); j++ )
     {    vec<int> lens;
          for ( int k = 0; k < L[j].isize( ); k++ )
          {    int len = 0;
               for ( int l = 0; l < L[j][k].isize( ); l++ )
               {    int d = L[j][k][l];
                    if ( IsSequence( D.O(d) ) )
                    {    basevector x;
                         int ltrim, rtrim;
                         GapToSeq( D.O(d), ltrim, rtrim, x );
                         pos += x.isize( ) - K + 1 - ltrim - rtrim;    }
                    if ( D.O(d)[0] < 0 ) continue;
                    for ( int m = 0; m < D.O(d).isize( ); m++ )
                    {    int e = D.O(d)[m];
                         if ( k == 0 )
                         {
                         const vec<pair<int,int>>& locs = linelocs[e];
                         if ( kmers[e] < MIN_KMERS_TO_SHOW || locs.size( ) > 2
                              || ( locs.size( ) == 2 
                                   && locs[0].first != locs[1].first ) )
                         {    len += kmers[e];
                              continue;    }
                         int re = inv[e];
                         if ( alignsb[e].size( ) + alignsb[re].size( ) != 1 )
                         {    len += kmers[e];
                              continue;    }
                         if ( alignsb[e].size( ) == 1 )
                         {    int g = alignsb[e][0].first;
                              if ( g > MAX_CHR )
                              {    len += kmers[e];
                                   continue;    }
                              int start = alignsb[e][0].second; 
                              int stop = alignsb[e][0].third;
                              if ( start < 0 || stop <= start ) // shouldn't happen
                              {    len += kmers[e];
                                   continue;    }
                              view0.push( g, True, ho_interval(start,stop),
                                   ho_interval( 
                                        pos + len, pos + len + kmers[e] ) );    }
                         if ( alignsb[re].size( ) == 1 )
                         {    int g = alignsb[re][0].first;
                              if ( g > MAX_CHR )
                              {    len += kmers[e];
                                   continue;    }
                              int start = alignsb[re][0].second; 
                              int stop = alignsb[re][0].third;
                              if ( start < 0 || stop <= start ) // shouldn't happen
                              {    len += kmers[e];
                                   continue;    }
                              view0.push( g, False, ho_interval(start,stop),
                                   ho_interval( 
                                        pos + len, pos + len + kmers[e] ) );    }
                         }
                         len += kmers[e];    }    }
               lens.push_back(len);    }
          Sort(lens);
          if ( lens.nonempty( ) ) pos += Median(lens);    }

     // Condense.

     /*
     cout << "\ninitial view" << endl;
     for ( int i = 0; i < view0.isize( ); i++ )
     {    cout << "[" << i+1 << "] " << ( view0[i].second ? "+" : "-" )
               << view0[0].first << ":" << view0[i].third
               << "   [" << view0[i].fourth << "]" << endl;    }
     */

     int g = -1, gstart = -1, gstop = -1;
     int astart = -1, astop = -1;
     Bool fw = True;
     for ( int i = 0; i < view0.isize( ); i++ )
     {    if ( g < 0 )
          {    g = view0[i].first;
               fw = view0[i].second;
               gstart = view0[i].third.Start( ), gstop = view0[i].third.Stop( );
               astart = view0[i].fourth.Start( ), astop = view0[i].fourth.Stop( );
               continue;    }
          Bool newview = False;
          if ( view0[i].first != g || view0[i].second != fw ) newview = True;
          else if ( fw && Abs( view0[i].third.Start( ) - gstop ) > MAX_GAP ) 
               newview = True;
          else if ( !fw && Abs( gstart - view0[i].third.Stop( ) ) > MAX_GAP ) 
               newview = True;
          if (newview)
          {    view.push_back( make_quad( g, fw, 
                    ho_interval(gstart,gstop), ho_interval(astart,astop) ) );
               g = view0[i].first;
               fw = view0[i].second;
               gstart = view0[i].third.Start( ), gstop = view0[i].third.Stop( );
               astart = view0[i].fourth.Start( ), astop = view0[i].fourth.Stop( );
                    }
          else
          {    gstart = Min( gstart, view0[i].third.Start( ) );
               gstop = Max( gstop, view0[i].third.Stop( ) );    
               astart = Min( astart, view0[i].fourth.Start( ) );
               astop = Max( astop, view0[i].fourth.Stop( ) );    }    }
     if ( g >= 0 ) 
     {    view.push_back( make_quad( g, fw, ho_interval(gstart,gstop),
               ho_interval(astart,astop) ) );    }    

     /*
     cout << "\ncondensed view" << endl;
     for ( int i = 0; i < (int) view.size( ); i++ )
     {    cout << "[" << i+1 << "] " << ( view[i].second ? "+" : "-" )
               << view[0].first << ":" << view[i].third
               << "   [" << view[i].fourth << "]" << endl;    }
     */

     }

template void View( const int l, const int K, const vec<int>& kmers, 
     const vec<int>& inv, 
     const digraphE<vec<int>>& D, const vec<vec<vec<vec<int>>>>& dlines,
     const vec<vec< pair<int,int> >>& linelocs,
     VirtualMasterVec< SerfVec< triple<int,int,int> > >& alignsb,
     SerfVec< quad<int,Bool,ho_interval,ho_interval> >& view );

template void View( const int l, const int K, const vec<int>& kmers, 
     const vec<int>& inv, 
     const digraphE<vec<int>>& D, const vec<vec<vec<vec<int>>>>& dlines,
     const vec<vec< pair<int,int> >>& linelocs,
     MasterVec< SerfVec< triple<int,int,int> > >& alignsb,
     SerfVec< quad<int,Bool,ho_interval,ho_interval> >& view );

typedef quad<int,Bool,ho_interval,ho_interval> porkypie;
extern template class SmallVec<porkypie,MempoolAllocator<porkypie>>;
extern template class OuterVec< SerfVec<porkypie> >;
template class SmallVec<porkypie,MempoolAllocator<porkypie> >;
template class OuterVec< SerfVec<porkypie> >;
