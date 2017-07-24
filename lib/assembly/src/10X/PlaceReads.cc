// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>

#include "CoreTools.h"
#include "math/Functions.h"
#include "math/HoInterval.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "paths/long/large/Lines.h"
#include "10X/DfTools.h"
#include "10X/PlaceReads.h"
#include "ParallelVecUtilities.h"
#include "10X/CleanThe.h"


// Align a path to the supergraph.

void Align( 
     const digraphE<vec<int>>& D, // supergraph
     const vec<int>& to_left,     // edge to left
     const vec<int>& to_right,    // edge to right
     const vec<int>& x,           // input path
     ho_interval& xpos,           // aligned start/stop positions on x
     vec<int>& d,                 // path in D
     int& dstart,                 // index of start edge on first edge in d
     int& dstop                   // index+1 of stop edge on last edge in d
          )
{    
     // Extend right.

     while(1)
     {    if ( xpos.Stop( ) == x.isize( ) ) break;
          else if ( dstop < D.O( d.back( ) ).isize( ) )
          {    if ( x[ xpos.Stop( ) ] != D.O( d.back( ) )[dstop] ) break;
               else
               {    xpos.AddToStop(1);
                    dstop++;    }    }
          else
          {    int v = to_right[ d.back( ) ], dn = -1;
               for ( int j = 0; j < D.From(v).isize( ); j++ )
               {    int n = D.IFrom( v, j );
                    if ( D.O(n)[0] < 0 ) continue;
                    if ( D.O(n)[0] == x[ xpos.Stop( ) ] )
                    {    dn = n;
                         break;    }    }
               if ( dn < 0 ) break;
               xpos.AddToStop(1);
               d.push_back(dn);
               dstop = 1;    }    }

     // Extend left.

     d.ReverseMe( );
     while(1)
     {    if ( xpos.Start( ) == 0 ) break;
          else if ( dstart > 0 )
          {    if ( x[ xpos.Start( ) - 1 ] != D.O( d.back( ) )[dstart-1] ) break;
               {    xpos.AddToStart(-1);
                    dstart--;    }    }
          else
          {    int v = to_left[ d.back( ) ], dn = -1;
               for ( int j = 0; j < D.To(v).isize( ); j++ )
               {    int n = D.ITo( v, j );
                    if ( D.O(n)[0] < 0 ) continue;
                    if ( D.O(n).back( ) == x[ xpos.Start( ) - 1 ] )
                    {    dn = n;
                         break;    }    }
               if ( dn < 0 ) break;
               xpos.AddToStart(-1);
               d.push_back(dn);
               dstart = D.O(dn).isize( ) - 1;    }    }
     d.ReverseMe( );    }

// !!! THIS VERSION IS ACTUALLY BUGGY!
// Align2: this version more correctly handles duplicate branches, as would arise
// in case zippering was incomplete.  Probably this should supplant Align at some
// point.

void Align2( 
     const digraphE<vec<int>>& D, // supergraph
     const vec<int>& to_left,     // edge to left
     const vec<int>& to_right,    // edge to right
     const vec<int>& x,           // input path
     ho_interval& xpos,           // aligned start/stop positions on x
     vec<int>& d,                 // path in D
     int& dstart,                 // index of start edge on first edge in d
     int& dstop                   // index+1 of stop edge on last edge in d
          )
{    
     // Sanity check D, to_left and to_right.  
     // Commented out because insanely expensive.

     /*
     #ifndef NDEBUG
     {    vec<int> to_leftx, to_rightx;
          D.ToLeft(to_leftx), D.ToRight(to_rightx);
          Assert( to_left == to_leftx );
          Assert( to_right == to_rightx );
          for ( int f = 0; f < D.E( ); f++ ) Assert( D.O(f).size( ) > 0 );    }
     #endif
     */

     // Sanity check the rest of the input data.
     // Commented out because not sufficiently tested, but probably correct.

     /*
     AssertGt( x.isize( ), 0 );
     AssertGt( d.isize( ), 0 );
     AssertGe( dstart, 0 );
     AssertGt( dstop, 0 );
     AssertLe( dstop, D.O( d.back( ) ).isize( ) );
     AssertGe( xpos.Start( ), 0 );
     AssertLe( xpos.Stop( ), x.isize( ) );
     #ifndef NDEBUG
     {    for ( int i = 0; i < d.isize( ) - 1; i++ )
               AssertEq( to_right[ d[i] ], to_left[ d[i+1] ] );
          if ( d.solo( ) )
          {    AssertLt( dstart, dstop );    }
          vec<int> z1, z2;
          for ( int j = xpos.Start( ); j < xpos.Stop( ); j++ ) z1.push_back( x[j] );
          if ( d.solo( ) )
          {    for ( int j = dstart; j < dstop; j++ )
                    z2.push_back( D.O( d[0] )[j] );    }
          else
          {    for ( int j = dstart; j < D.O(d[0]).isize( ); j++ )
                    z2.push_back( D.O( d.front( ) )[j] );
               for ( int j = 1; j < d.isize( ) - 1; j++ ) z2.append( D.O(d[j]) );
               for ( int j = 0; j < dstop; j++ )
                    z2.push_back( D.O( d.back( ) )[j] );    }
          Assert( z1 == z2 );    }
     #endif
     */

     // Extend right.

     while(1)
     {    if ( xpos.Stop( ) == x.isize( ) ) break;
          else if ( dstop < D.O( d.back( ) ).isize( ) )
          {    if ( x[ xpos.Stop( ) ] != D.O( d.back( ) )[dstop] ) break;
               else
               {    xpos.AddToStop(1);
                    dstop++;    }    }
          else
          {    int v = to_right[ d.back( ) ], dn = -1;
               int count = 0;
               for ( int j = 0; j < D.From(v).isize( ); j++ )
               {    int n = D.IFrom( v, j );
                    if ( D.O(n)[0] < 0 ) continue;
                    if ( D.O(n)[0] == x[ xpos.Stop( ) ] )
                    {    dn = n;
                         count++;    }    }
               if ( count == 0 ) break;
               if ( count == 1 )
               {    xpos.AddToStop(1);
                    d.push_back(dn);
                    dstop = 1;    }
               else
               {    
                    vec<vec<int>> exts;
                    for ( int j = 0; j < D.From(v).isize( ); j++ )
                    {    int n = D.IFrom( v, j );
                         if ( D.O(n)[0] < 0 ) continue;
                         if ( D.O(n)[0] == x[ xpos.Stop( ) ] )
                         {    int v2 = to_right[n];
                              if ( D.From(v2).size( ) == 0 ) exts.push_back({n});
                              else for ( int l = 0; l < D.From(v2).isize( ); l++ )
                              {    vec<int> z = { n, D.IFrom(v2,l) };
                                   exts.push_back(z);    }    }    }
                    vec<int> matches;
                    int len = x.isize( ) - xpos.Stop( );
                    for ( int j = 0; j < exts.isize( ); j++ )
                    {    const vec<int>& z = exts[j];
                         vec<int> t;
                         for ( auto m : z )
                         {    for ( auto r : D.O(m) )
                              {    t.push_back(r);
                                   if ( t.isize( ) == len ) break;    }
                              if ( t.isize( ) == len ) break;    }
                         Bool ok = True;
                         if ( t.isize( ) < len ) ok = False;
                         for ( int l = 0; l < t.isize( ); l++ )
                         {    if ( t[l] != x[ xpos.Stop( ) + l ] )
                              {    ok = False;
                                   break;    }    }
                         if (ok) matches.push_back(j);    }
                    if ( !matches.solo( ) ) break;
                    xpos.AddToStop(1);
                    d.push_back( exts[ matches[0] ][0] );
                    dstop = 1;    }    }    }

     // Extend left.

     d.ReverseMe( );
     while(1)
     {    if ( xpos.Start( ) == 0 ) break;
          else if ( dstart > 0 )
          {    if ( x[ xpos.Start( ) - 1 ] != D.O( d.back( ) )[dstart-1] ) break;
               {    xpos.AddToStart(-1);
                    dstart--;    }    }
          else
          {    int v = to_left[ d.back( ) ], dn = -1;
               int count = 0;
               for ( int j = 0; j < D.To(v).isize( ); j++ )
               {    int n = D.ITo( v, j );
                    if ( D.O(n)[0] < 0 ) continue;
                    if ( D.O(n).back( ) == x[ xpos.Start( ) - 1 ] )
                    {    dn = n;
                         count++;    }    }
               if ( count == 0 ) break;
               if ( count == 1 )
               {    xpos.AddToStart(-1);
                    d.push_back(dn);
                    dstart = D.O(dn).isize( ) - 1;    }
               else
               {    vec<vec<int>> exts;
                    for ( int j = 0; j < D.To(v).isize( ); j++ )
                    {    int n = D.ITo( v, j );
                         if ( D.O(n)[0] < 0 ) continue;
                         if ( D.O(n).back( ) == x[ xpos.Start( ) - 1 ] )
                         {    int v2 = to_left[n];
                              if ( D.To(v2).size( ) == 0 ) exts.push_back({n});
                              else for ( int l = 0; l < D.To(v2).isize( ); l++ )
                              {    vec<int> z = { n, D.ITo(v2,l) };
                                   exts.push_back(z);    }    }    }
                    vec<int> matches;
                    int len = xpos.Start( );
                    for ( int j = 0; j < exts.isize( ); j++ )
                    {    const vec<int>& z = exts[j];
                         vec<int> t;
                         for ( auto m : z )
                         {    for ( int k = D.O(m).isize( ) - 1; k >= 0; k-- )
                              {    int r = D.O(m)[k];
                                   t.push_back(r);
                                   if ( t.isize( ) == len ) break;    }
                              if ( t.isize( ) == len ) break;    }
                         Bool ok = True;
                         if ( t.isize( ) < len ) ok = False;
                         for ( int l = 0; l < t.isize( ); l++ )
                         {    if ( t[l] != x[ xpos.Start( ) - l - 1 ] )
                              {    ok = False;
                                   break;    }    }
                         if (ok) matches.push_back(j);    }
                    if ( !matches.solo( ) ) break;
                    xpos.AddToStart(-1);
                    int n = exts[ matches[0] ][0];
                    d.push_back(n);
                    dstart = D.O(n).isize( ) - 1;    }    }    }
     d.ReverseMe( );    }

void ExtendRight(const digraphE<vec<int>>& D, int w, int edge, const vec<int>& x, ho_interval xpos, int dseed,
        vec<int>& sofar, vec<triple<int,vec<int>,ho_interval>>& hit){ 
     // artifical brakes
     if(hit.size()>20) return;

     int len = x.isize() - xpos.Stop(); 
     int seeklen = Min(D.O(edge).isize()-dseed,len);         
     Bool ok = True;
     int l = 0;
     for(l = 0; l < seeklen; l++){
         if(x[xpos.Stop()+l]!=D.O(edge)[dseed+l]){
             ok = False;
             break;
         }
     }
     dseed += l;
     sofar.push_back(edge);
     xpos.AddToStop(l);

     if(!ok || xpos.Stop()==x.isize()){
         hit.push(dseed,sofar,xpos);
         sofar.resize(sofar.size()-1);
         return;
     }
    
     int ext = 0;
     for(int k = 0; k < D.From(w).isize(); k++){
         if(D.O(D.IFrom(w,k))[0]<0) continue;
         if(x[xpos.Stop()]==D.O(D.IFrom(w,k))[0]){
             ext++;
             ExtendRight(D, D.From(w)[k],D.IFrom(w,k),x,xpos,0,sofar,hit);
         }
     }
     if(ext==0)
         hit.push(dseed,sofar,xpos);
     sofar.resize(sofar.size()-1);
}

void ExtendLeft(const digraphE<vec<int>>& D, int v, int edge, const vec<int>& x, ho_interval xpos, int dseed,
        vec<int>& sofar, vec<triple<int,vec<int>,ho_interval>>& hit){ 
     // artifical brakes
     if(hit.size()>20) return;

     int len = xpos.Start(); 
     int seeklen = Min(dseed,len);         
     Bool ok = True;
     int l = 0;
     for(l = 0; l < seeklen; l++){
         if(x[xpos.Start()-l-1]!=D.O(edge)[dseed-l-1]){
             ok = False;
             break;
         }
     }
     dseed -= l;
     sofar.push_back(edge);
     xpos.SubFromStart(l);

     if(!ok || xpos.Start()==0){
         hit.push(dseed,sofar,xpos);
         sofar.resize(sofar.size()-1);
         return;
     }

     int ext = 0;
     for(int k = 0; k < D.To(v).isize(); k++){
         if(D.O(D.ITo(v,k))[0]<0) continue;
         if(x[xpos.Start()-1]==D.O(D.ITo(v,k)).back()){
             ext++;
             ExtendLeft(D, D.To(v)[k],D.ITo(v,k),x,xpos,D.O(D.ITo(v,k)).size(),sofar,hit);
         }
     }
     if(ext==0)
         hit.push(dseed,sofar,xpos);
     sofar.resize(sofar.size()-1);
}

// This version handles unzippered vertices carefully. It first creates path extension to left and right
// from a seed edge. Then it picks the maximum length unambiguous path in either direction and 
// joins them at the seed to form the final path
void Align2_new( 
     const digraphE<vec<int>>& D, // supergraph
     const vec<int>& to_left,     // edge to left
     const vec<int>& to_right,    // edge to right
     const vec<int>& x,           // input path
     ho_interval& xpos,           // aligned start/stop positions on x
     vec<int>& d,                 // path in D
     int& dstart,                 // index of start edge on first edge in d
     int& dstop                   // index+1 of stop edge on last edge in d
          )
{    
     // process hits
     vec<triple<int,vec<int>,ho_interval>> hitR, hitL;
     triple<int,vec<int>,ho_interval> qL, qR;
     vec<int> sofar;
     int edge = d.back();
     int w = to_right[edge];
     int v = to_left[edge];

     ExtendRight(D, w, edge, x, xpos, dstop, sofar, hitR);
     ForceAssertEq(sofar.size(),0);
     ExtendLeft(D, v, edge, x, xpos, dstart, sofar, hitL);
     ForceAssertEq(sofar.size(),0);

     /* for(auto h : hitL) */
     /*     cout << h.first << " : " << printSeq(h.second) << " : " << h.third.Start() << " " << h.third.Stop() << endl; */
     /* cout << " -- " << endl; */
     /* for(auto h : hitR) */
     /*     cout << h.first << " : " << printSeq(h.second) << " : " << h.third.Start() << " " << h.third.Stop() << endl; */
     /* cout << endl; */

     vec<int> mR,mL;
     vec<Bool> delR(hitR.size(),False),delL(hitL.size(),False);
     int idx;

     for(auto& h : hitR)
         mR.push_back(h.third.Stop());
     SortSync(mR,hitR);
     idx = 0;
     while(mR.size()>1 && idx+1 < mR.isize()){
        if(mR[idx]==mR[idx+1]){
            delR[idx] = True;
            delR[idx+1] = True;
        }
        idx++;
     }
     auto keepR = hitR;
     auto kR = mR;
     EraseIf(keepR,delR);
     EraseIf(kR,delR);

     for(auto& h : hitL)
         mL.push_back(h.third.Start());
     ReverseSortSync(mL,hitL);
     idx = 0;
     while(mL.size()>1 && idx+1 < mL.isize()){
        if(mL[idx]==mL[idx+1]){
            delL[idx] = True;
            delL[idx+1] = True;
        }
        idx++;
     }
     auto keepL = hitL;
     auto kL = mL;
     EraseIf(keepL,delL);
     EraseIf(kL,delL);

     // if no single path, trim to basic
     if(kR.size()==0){
         qR = hitR[0];
         qR.second.resize(1);
         xpos.AddToStop( D.O(qR.second.back()).isize()-dstop );
         dstop = D.O(qR.second.back()).isize();
     }else{
         qR = keepR.back(); // max stop
         xpos.SetStop( qR.third.Stop() );
         dstop = qR.first;
     }

     // if no single path, trim to basic
     if(kL.size()==0){
         qL = hitL[0];
         qL.second.resize(1);
         xpos.SubFromStart( dstart );
         dstart = 0;
     }else{
         qL = keepL.back();
         xpos.SetStart( qL.third.Start() );
         dstart = qL.first;
     }
       
     // sanity
     ForceAssertGe(qL.second.size(),1);
     ForceAssertGe(qR.second.size(),1);
     ForceAssertEq(qL.third.Stop(),qR.third.Start()+1);

     // build d
     d = qL.second;
     d.ReverseMe();
     d.insert(d.end(),qR.second.begin()+1,qR.second.end());

     /* cout << "d " << printSeq(d) << " dstart "  << dstart << " dstop " << dstop << " xpos " << xpos.Start() << " " << xpos.Stop() << " : " << x.isize() << endl; */

}

// This is the correct version of FindPlaces that should be used. It uses all min mult edges as seeds.
Bool FindPlacesAlt( const ReadPath& p, vec<int>& d, vec<int>& x, 
     const HyperBasevectorX& hb, const digraphE<vec<int>>& D, const vec<int>& dlens,
     const vec<int>& to_left, const vec<int>& to_right,
     const vec<vec<pair<int,int>>>& nd, const Bool align2,
     int& nplaces, int& pos, vec<int>& aplace )
{
     // Empty path?
     if ( p.size( ) == 0 ) return False;

     // Find minimum multiplicity edge in p.
     int m = 1000000000, mpe = -1;
     for ( int j = 0; j < (int) p.size( ); j++ )
     {    int n = nd[ p[j] ].size( );
          if ( n == 0 ) continue;
          if ( n < m )
          {    m = n;
               mpe = j;    }    }

     // Find all min mult edges
     vec<int> minE;
     for( int j = 0; j < (int) p.size( ); j++)
         if(m==nd[p[j]].isize( ))
             minE.push(j);


     // Test placements.
     nplaces = 0, pos = -1;
     vec<triple<ho_interval,vec<int>,int>> accepted;
     for(auto mp : minE){
         if ( mp < 0 || nd[ p[mp] ].empty( ) ) return False; // pathological case not gonna happen
         for ( int i = 0; i < nd[ p[mp] ].isize( ); i++ )
         {    int e = nd[ p[mp] ][i].first, epos = nd[ p[mp] ][i].second;
              x.clear( );
              for ( int j = 0; j < (int) p.size( ); j++ ) x.push_back( p[j] );
              ho_interval xpos( mp, mp+1 );
              d = {e};
              int dstart = epos, dstop = epos + 1;
              ( !align2 ? Align : Align2 )(
                   D, to_left, to_right, x, xpos, d, dstart, dstop );

              // For now, toss improper alignments.
              if ( xpos.Start( ) > 0 )
              {    if ( dstart > 0 ) continue;
                   if ( dstart == 0 )
                   {    int v = to_left[ d.front( ) ];
                        if ( D.To(v).nonempty( ) ) continue;    }    }
              if ( xpos.Stop( ) < x.isize( ) )
              {    if ( dstop < D.O( d.back( ) ).isize( ) ) continue;
                   int w = to_right[ d.back( ) ];
                   if ( D.From(w).nonempty( ) ) continue;    }

              // Compute start position of read on d[0].
              pos = p.getOffset( );
              for ( int j = 0; j < dstart; j++ )
                   pos += hb.Kmers( D.O(d[0])[j] );

              // Save 
              accepted.push(xpos,d,pos);
         }
     }

     if(accepted.size()==0)
         return False;

     // mix of full and improper placements
     ParallelUniqueSort(accepted);
     vec<int> klen(accepted.size(),0);
     for(int j = 0; j < accepted.isize(); j++){
         auto xpos1 = accepted[j].first;
         for(int kk = xpos1.Start(); kk < xpos1.Stop(); kk++)
             klen[j] += hb.Kmers(x[kk]);
     }
     ReverseSortSync(klen,accepted);

     // if largest placement is not unique (be it proper or improper), return false;
     if(klen.size()>1)
         if(klen[0]==klen[1])
             return False;

     d = accepted[0].second;
     pos = accepted[0].third;
     nplaces = 1;
     aplace.clear( );
     for ( int j = 0; j < d.isize( ); j++ ) aplace.push_back( d[j] );    

     // Check for broken placement.  Need to track down how this could ever happen.
     int n = hb.K( ) - 1;
     for ( auto f : aplace ) n += dlens[f];
     if ( pos > n ) return False;
     return True;    }

Bool FindPlaces( const ReadPath& p, vec<int>& d, vec<int>& x, 
     const HyperBasevectorX& hb, const digraphE<vec<int>>& D, const vec<int>& dlens,
     const vec<int>& to_left, const vec<int>& to_right,
     const vec<vec<pair<int,int>>>& nd, const Bool align2,
     int& nplaces, int& pos, vec<int>& aplace )
{
     // Empty path?

     if ( p.size( ) == 0 ) return False;

     // Find minimum multiplicity edge in p.

     int m = 1000000000, mp = -1;
     for ( int j = 0; j < (int) p.size( ); j++ )
     {    int n = nd[ p[j] ].size( );
          if ( n == 0 ) continue;
          if ( n < m )
          {    m = n;
               mp = j;    }    }

     // Test placements.

     if ( mp < 0 || nd[ p[mp] ].empty( ) ) return False;
     nplaces = 0, pos = -1;
     for ( int i = 0; i < nd[ p[mp] ].isize( ); i++ )
     {    int e = nd[ p[mp] ][i].first, epos = nd[ p[mp] ][i].second;
          x.clear( );
          for ( int j = 0; j < (int) p.size( ); j++ ) x.push_back( p[j] );
          ho_interval xpos( mp, mp+1 );
          d = {e};
          int dstart = epos, dstop = epos + 1;
          ( !align2 ? Align : Align2 )(
               D, to_left, to_right, x, xpos, d, dstart, dstop );

          // For now, toss improper alignments.

          if ( xpos.Start( ) > 0 )
          {    if ( dstart > 0 ) continue;
               if ( dstart == 0 )
               {    int v = to_left[ d.front( ) ];
                    if ( D.To(v).nonempty( ) ) continue;    }    }
          if ( xpos.Stop( ) < x.isize( ) )
          {    if ( dstop < D.O( d.back( ) ).isize( ) ) continue;
               int w = to_right[ d.back( ) ];
               if ( D.From(w).nonempty( ) ) continue;    }

          // Compute start position of read on d[0].

          pos = p.getOffset( );
          for ( int j = 0; j < dstart; j++ )
               pos += hb.Kmers( D.O(d[0])[j] );

          // Give up if more than one placement.

          nplaces++;
          if ( nplaces > 1 ) break;

          // Save.  Just saving d, but could save ( xpos, d, dstart, dstop ).

          aplace.clear( );
          for ( int j = 0; j < d.isize( ); j++ ) aplace.push_back( d[j] );    }

     // For now require that there is a unique proper placement.

     if ( nplaces != 1 ) return False;

     // Check for broken placement.  Need to track down how this could ever happen.

     int n = hb.K( ) - 1;
     for ( auto f : aplace ) n += dlens[f];
     if ( pos > n ) return False;
     return True;    }

void PlaceReadsAlt( const HyperBasevectorX& hb, const ReadPathVec& paths, 
     const vec<Bool>& dup, const digraphE<vec<int>>& D, ReadPathVec& dpaths,
     const Bool verbose, const Bool single, const Bool align2 )
{
     // Compute dlens.  This is only used in the sanity check below.  It would be
     // good to eliminate this calculation if possible.  Then index D.

     vec<int> dlens;
     ComputeDlens( hb, D, dlens );
     vec<vec<pair<int,int>>> nd;
     MakeIndex( hb, D, nd, verbose );

     // Build paths.

     if ( dpaths.size( ) == 0 )
     {    if (verbose) cout << Date( ) << ": defining data structure" << endl;
          dpaths = paths;    } // TERRIBLE WAY TO RESERVE SPACE!!!!!!!!!!!!!!!!
     else {
         dpaths.resize(paths.size());
     }
     if (verbose) cout << Date( ) << ": creating to_left and to_right" << endl;
     vec<int> to_left, to_right;
     D.ToLeft(to_left), D.ToRight(to_right);
     if (verbose) cout << Date( ) << ": looking at reads" << endl;
     double clock = WallClockTime( );
     const int batch = 100000;
     int nthreads = ( single ? 1 : omp_get_max_threads( ) );
     #pragma omp parallel for num_threads(nthreads) schedule(dynamic,1)
     for ( int64_t bi = 0; bi < (int64_t) paths.size( ); bi += batch )
     {    vec<int> aplace, d, x;
          for ( int64_t id = bi; 
               id < Min( bi + batch, (int64_t) paths.size( ) ); id++ )
          {    dpaths[id].resize(0); 
               if ( dup[id/2] ) continue;
               const ReadPath &p = paths[id]; 
               int nplaces, pos;
               if ( !FindPlacesAlt( p, d, x, hb, D, dlens, to_left, to_right, 
                    nd, align2, nplaces, pos, aplace ) )
               {    continue;    }
               dpaths[id].resize( aplace.size( ) );
               for ( int j = 0; j < aplace.isize( ); j++ )
                    dpaths[id][j] = aplace[j];
               dpaths[id].setOffset(pos);    }     }
     if (verbose) cout << Date( ) << ": done, used " << TimeSince(clock) << endl;
     int64_t placed = 0;
     for ( int64_t id = 0; id < (int64_t) paths.size( ); id++ )
          if ( dpaths[id].size( ) > 0 ) placed++;
     if (verbose)
     {    cout << Date( ) << ": " 
               << PERCENT_RATIO( 3, placed, (int64_t) paths.size( ) ) 
               << " placed" << endl;    }    }

void ComputeDlens( const HyperBasevectorX& hb, const digraphE<vec<int>>& D,
     vec<int>& dlens )
{    dlens.resize_and_set( D.E( ), 0 );
     #pragma omp parallel for
     for ( int e = 0; e < D.E( ); e++ )
     {    if ( D.O(e)[0] < 0 ) continue;
          for ( int j = 0; j < D.O(e).isize( ); j++ )
               dlens[e] += hb.Kmers( D.O(e)[j] );    }    }

void MakeIndex( const HyperBasevectorX& hb, const digraphE<vec<int>>& D,
     vec<vec<pair<int,int>>>& nd, const Bool verbose )
{    if (verbose) cout << Date( ) << ": creating index for read placement" << endl;
     vec<Bool> used;
     D.Used(used);
     nd.clear_and_resize( hb.E( ) );
     for ( int e = 0; e < D.E( ); e++ )
     {    if ( !used[e] ) continue;
          const vec<int>& x = D.O(e);
          if ( x[0] < 0 ) continue;
          for ( int j = 0; j < x.isize( ); j++ )
               nd[ x[j] ].push( e, j );    }    }

void PlaceReads2( const HyperBasevectorX& hb, const ReadPathVec& paths, 
     vec<int64_t>& maprr, const vec<Bool>& dup, const digraphE<vec<int>>& D, 
     MasterVec<IntVec>& dpaths, int rd_set,
     const Bool verbose, const Bool single, const Bool align2 )
{
     // Compute dlens.  This is only used in the sanity check below.  It would be
     // good to eliminate this calculation if possible.  Then index D.

     vec<int> dlens;
     ComputeDlens( hb, D, dlens );
     vec<vec<pair<int,int>>> nd;
     MakeIndex( hb, D, nd, verbose );

     // Build paths.

     if ( dpaths.size( ) == 0 )
     {    if (verbose) cout << Date( ) << ": defining data structure" << endl;
          dpaths.resize(rd_set);    }
     if (verbose) cout << Date( ) << ": creating to_left and to_right" << endl;
     vec<int> to_left, to_right;
     D.ToLeft(to_left), D.ToRight(to_right);
     if (verbose) cout << Date( ) << ": looking at reads" << endl;
     double clock = WallClockTime( );
     const int batch = 100000;
     int nthreads = ( single ? 1 : omp_get_max_threads( ) );
     #pragma omp parallel for num_threads(nthreads) schedule(dynamic,1)
     for ( int64_t bi = 0; bi < (int64_t) paths.size( ); bi += batch )
     {    vec<int> aplace, d, x;
          for ( int64_t id = bi;
               id < Min( bi + batch, (int64_t) paths.size( ) ); id++ )
          {    dpaths[maprr[id]].resize(0);
               if ( dup[id/2] ) continue;
               const ReadPath &p = paths[id];
               int nplaces, pos;
               if ( !FindPlaces( p, d, x, hb, D, dlens, to_left, to_right, 
                    nd, align2, nplaces, pos, aplace ) )
               {    continue;    }
               dpaths[maprr[id]].resize( aplace.size( ) );
               for ( int j = 0; j < aplace.isize( ); j++ )
                    dpaths[maprr[id]][j] = aplace[j]; }     }

     if (verbose) cout << Date( ) << ": done, used " << TimeSince(clock) << endl;
     int64_t placed = 0;
     for ( int64_t id = 0; id < (int64_t) paths.size( ); id++ )
          if ( dpaths[maprr[id]].size( ) > 0 ) placed++;
     if (verbose)
     {    cout << Date( ) << ": " 
               << PERCENT_RATIO( 3, placed, (int64_t) paths.size( ) ) 
               << " placed" << endl;    }    }

void PlaceReads2( const HyperBasevectorX& hb, const ReadPathVecX& paths, 
     vec<int64_t>& maprr, const vec<Bool>& dup, const digraphE<vec<int>>& D, 
     MasterVec<IntVec>& dpaths, int rd_set,
     const Bool verbose, const Bool single, const Bool align2 )
{
     // Compute dlens.  This is only used in the sanity check below.  It would be
     // good to eliminate this calculation if possible.  Then index D.

     vec<int> dlens;
     ComputeDlens( hb, D, dlens );
     vec<vec<pair<int,int>>> nd;
     MakeIndex( hb, D, nd, verbose );

     // Build paths.


     int64_t to_place = 0;
     #pragma omp parallel for reduction(+:to_place)
     for ( int64_t pid = 0; pid < dup.jsize(); pid++ ) 
     {    if ( dup[pid] ) continue;
          to_place++;    }

     // Reserve space when placing a large number of paths.

     if ( to_place > 10000000 && dpaths.size( ) == 0 )
     {    if (verbose) cout << Date( ) << ": defining data structure" << endl;
          // BETTER ALLOCATION NEEDED? !!!!!!!!!!!!!!!!!!!!!!!!!!
          dpaths.resize(rd_set);
          int mpl = 1;// heuristic
          for(int64_t id = 0; id<rd_set; id++)
               dpaths[id].reserve(mpl);    }
     else{
         dpaths.resize(rd_set);
     }
     if (verbose) cout << Date( ) << ": creating to_left and to_right" << endl;
     vec<int> to_left, to_right;
     D.ToLeft(to_left), D.ToRight(to_right);
     if (verbose) cout << Date( ) << ": looking at reads" << endl;
     double clock = WallClockTime( );
     const int batch = 100000;
     int nthreads = ( single ? 1 : omp_get_max_threads( ) );
     #pragma omp parallel for num_threads(nthreads) schedule(dynamic,1)
     for ( int64_t bi = 0; bi < (int64_t) paths.size( ); bi += batch )
     {    vec<int> aplace, d, x;
          ReadPath p;
          for ( int64_t id = bi;
               id < Min( bi + batch, (int64_t) paths.size( ) ); id++ )
          {    dpaths[maprr[id]].resize(0);
               if ( dup[id/2] ) continue;
               p.clear( );
               paths.unzip(p,hb,id);
               int nplaces, pos;
               if ( !FindPlaces( p, d, x, hb, D, dlens, to_left, to_right, 
                    nd, align2, nplaces, pos, aplace ) )
               {    continue;    }
               dpaths[maprr[id]].resize( aplace.size( ) );
               for ( int j = 0; j < aplace.isize( ); j++ )
                    dpaths[maprr[id]][j] = aplace[j]; }     }
     if (verbose) cout << Date( ) << ": done, used " << TimeSince(clock) << endl;
     int64_t placed = 0;
     for ( int64_t id = 0; id < (int64_t) paths.size( ); id++ )
          if ( dpaths[maprr[id]].size( ) > 0 ) placed++;
     if (verbose)
     {    cout << Date( ) << ": " 
               << PERCENT_RATIO( 3, placed, (int64_t) paths.size( ) ) 
               << " placed" << endl;    }    }

void PlaceReads( const HyperBasevectorX& hb, const ReadPathVec& paths, 
     const vec<Bool>& dup, const digraphE<vec<int>>& D, ReadPathVec& dpaths,
     const Bool verbose, const Bool single, const Bool align2 )
{
     // Compute dlens.  This is only used in the sanity check below.  It would be
     // good to eliminate this calculation if possible.  Then index D.

     vec<int> dlens;
     ComputeDlens( hb, D, dlens );
     vec<vec<pair<int,int>>> nd;
     MakeIndex( hb, D, nd, verbose );

     // Build paths.

     if ( dpaths.size( ) == 0 )
     {    if (verbose) cout << Date( ) << ": defining data structure" << endl;
          dpaths = paths;    } // TERRIBLE WAY TO RESERVE SPACE!!!!!!!!!!!!!!!!
     else {
         dpaths.resize(paths.size());
     }
     if (verbose) cout << Date( ) << ": creating to_left and to_right" << endl;
     vec<int> to_left, to_right;
     D.ToLeft(to_left), D.ToRight(to_right);
     if (verbose) cout << Date( ) << ": looking at reads" << endl;
     double clock = WallClockTime( );
     const int batch = 100000;
     int nthreads = ( single ? 1 : omp_get_max_threads( ) );
     #pragma omp parallel for num_threads(nthreads) schedule(dynamic,1)
     for ( int64_t bi = 0; bi < (int64_t) paths.size( ); bi += batch )
     {    vec<int> aplace, d, x;
          for ( int64_t id = bi; 
               id < Min( bi + batch, (int64_t) paths.size( ) ); id++ )
          {    dpaths[id].resize(0); 
               if ( dup[id/2] ) continue;
               const ReadPath &p = paths[id]; 
               int nplaces, pos;
               if ( !FindPlaces( p, d, x, hb, D, dlens, to_left, to_right, 
                    nd, align2, nplaces, pos, aplace ) )
               {    continue;    }
               dpaths[id].resize( aplace.size( ) );
               for ( int j = 0; j < aplace.isize( ); j++ )
                    dpaths[id][j] = aplace[j];
               dpaths[id].setOffset(pos);    }     }
     if (verbose) cout << Date( ) << ": done, used " << TimeSince(clock) << endl;
     int64_t placed = 0;
     for ( int64_t id = 0; id < (int64_t) paths.size( ); id++ )
          if ( dpaths[id].size( ) > 0 ) placed++;
     if (verbose)
     {    cout << Date( ) << ": " 
               << PERCENT_RATIO( 3, placed, (int64_t) paths.size( ) ) 
               << " placed" << endl;    }    }

void PlaceReadsSmart( const HyperBasevectorX& hb, const ReadPathVec& paths, 
     const vec<Bool>& dup, const digraphE<vec<int>>& D, const vec<int>& dinv,
     ReadPathVec& dpaths, const vec<vec<vec<vec<int>>>>& dlines, 
     const vec<int64_t>& bci, const Bool verbose, const int btest, 
     const Bool align2 )
{
     // Set up.

     vec<int> dlens( D.E( ), 0 ), to_left, to_right, linv;
     #pragma omp parallel for
     for ( int e = 0; e < D.E( ); e++ )
     {    if ( D.O(e)[0] < 0 ) continue;
          for ( int j = 0; j < D.O(e).isize( ); j++ )
               dlens[e] += hb.Kmers( D.O(e)[j] );    }
     D.ToLeft(to_left), D.ToRight(to_right);
     LineInv( dlines, dinv, linv );

     // To each nongap edge, assign a line, and a starting position on that line.
     // For nonpalindromic lines L, we pick just one of L, rc(L), and place both
     // edges d and dinv[d] at the same place on that line.

     vec<quad<int,int,Bool,int>> tol( D.E( ), make_quad(-1,-1,True,-1) );
     cout << Date( ) << ": defining tol" << endl;
     for ( int i = 0; i < dlines.isize( ); i++ )
     {    if ( linv[i] < i ) continue;
          const vec<vec<vec<int>>>& L = dlines[i];
          int pos = 0;
          for ( int j = 0; j < L.isize( ); j++ )
          {    vec<int> lens;
               for ( int k = 0; k < L[j].isize( ); k++ )
               {    int len = 0;
                    for ( int l = 0; l < L[j][k].isize( ); l++ )
                    {    int d = L[j][k][l];
                         if ( D.O(d)[0] < 0 ) continue;
                         if ( tol[d].first < 0 )
                         {    tol[d] = make_quad( i, pos+len, True, j );
                              tol[ dinv[d] ] 
                                   = make_quad( i, pos+len, False, j );    }
                         len += dlens[d];    }
                    lens.push_back(len);    }
               Sort(lens);
               if ( lens.nonempty( ) ) pos += Median(lens);    }    }

     // Define the territory of each barcode.  The granularity of this definition
     // is weak, as we use only start/stop positions (i,j) on lines, whereas we
     // should use base positions.

     // Go through the barcodes.

     cout << Date( ) << ": start traversal" << endl;
     const int MAX_BC_GAP = 100000;
     const int MIN_BC_GROUP = 3;
     int ndots = 0, stopped = 0;
     #pragma omp parallel for schedule(dynamic, 1000)
     for ( int b = 1; b < bci.isize( ) - 1; b++ )
     {    if ( btest >= 0 && b != btest ) continue;

          // Logging.

          ostringstream out;
          if (verbose)
          {    out << "\nlooking at barcode " << b << " [" << bci[b+1]-bci[b]
                    << " reads]" << endl;    }

          // Define the territory of each barcode.  The granularity of this
          // definition is weak, as we use only start/stop positions (i,j) on lines,
          // whereas we should use base positions.

          int total = 0;
          vec<triple<int,int,int>> lplaces; // (line, start on line, line unit)
          for ( int64_t id = bci[b]; id < bci[b+1]; id++ )
          {    const ReadPath& p = dpaths[id];
               int pos = p.getOffset( );
               for ( int j = 0; j < (int) p.size( ); j++ )
               {    if ( j == 0 || tol[ p[j] ].first != tol[ p[j-1] ].first )
                    {    triple<int,int,int> x; 
                         x.first = tol[ p[j] ].first;
                         int xpos = pos;
                         if ( !tol[ p[j] ].third ) xpos = dlens[ p[j] ] - xpos;
                         x.second = tol[ p[j] ].second + xpos;
                         x.third = tol[ p[j] ].fourth;
                         lplaces.push_back(x);    }
                    pos -= dlens[ p[j] ];    }    }
          Sort(lplaces);
          vec<triple<int,int,int>> groups;
          for ( int i = 0; i < lplaces.isize( ); i++ )
          {    int j;
               for ( j = i + 1; j < lplaces.isize( ); j++ )
               {    if ( lplaces[j].first != lplaces[i].first ) break;
                    if ( lplaces[j].second - lplaces[j-1].second > MAX_BC_GAP ) 
                         break;    }
               if ( j - i >= MIN_BC_GROUP )
               {    int len = lplaces[j-1].second - lplaces[i].second;
                    int start = lplaces[i].third, stop = lplaces[j-1].third;
                    int l = lplaces[i].first;
                    const vec<vec<vec<int>>>& L = dlines[l];

                    // Extend start and stop positions.

                    const int EXT = 25000;
                    int left = 0, right = 0;
                    while( left < EXT && start >= 1 )
                    {    start--;
                         vec<int> lens;
                         for ( int k = 0; k < L[start].isize( ); k++ )
                         {    int len = 0;
                              for ( int l = 0; l < L[start][k].isize( ); l++ )
                              {    int d = L[start][k][l];
                                   len += dlens[d];    }
                              lens.push_back(len);    }
                         Sort(lens);
                         if ( lens.nonempty( ) ) left += Median(lens);    }
                    while( right < EXT && stop < L.isize( ) - 1 )
                    {    stop++;
                         vec<int> lens;
                         for ( int k = 0; k < L[stop].isize( ); k++ )
                         {    int len = 0;
                              for ( int l = 0; l < L[stop][k].isize( ); l++ )
                              {    int d = L[stop][k][l];
                                   len += dlens[d];    }
                              lens.push_back(len);    }
                         Sort(lens);
                         if ( lens.nonempty( ) ) right += Median(lens);    }
                    
                    // Save.

                    if (verbose)
                    {    out << l << "." << lplaces[i].second << "-" 
                              << lplaces[j-1].second << " [len=" << len 
                              << ",count=" << j-i 
                              << ",start=" << start << "/" << L[start][0][0]
                              << ",stop=" << stop << "/" << L[stop][0][0]
                              << "]" << endl;    }
                    total += len;
                    groups.push( lplaces[i].first, start, stop );    }
               i = j - 1;    }    
          if (verbose) out << "total = " << ToStringAddCommas(total) << endl;

          // Define the edges in the territory.

          vec<pair<int,int>> locs;
          int nd = 0;
          for ( int i = 0; i < groups.isize( ); i++ )
          {    const vec<vec<vec<int>>>& L = dlines[ groups[i].first ];
               int start = groups[i].second, stop = groups[i].third;
               for ( int j = start; j <= stop; j++ )
               {    const vec<vec<int>>& M = L[j];
                    vec<int> ds;
                    for ( int k = 0; k < M.isize( ); k++ )
                    for ( int l = 0; l < M[k].isize( ); l++ )
                    {    int d = M[k][l];
                         if ( D.O(d)[0] >= 0 ) ds.push_back( d, dinv[d] );    }
                    UniqueSort(ds);
                    nd += ds.size( );
                    for ( int k = 0; k < ds.isize( ); k++ )
                    {    int d = ds[k];
                         vec<int> es = D.O(d);
                         UniqueSort(es);
                         for ( auto e : es ) locs.push( e, d );    }    }    }
          Sort(locs);
          if (verbose) 
          {    out << "locs = " << ToStringAddCommas( locs.size( ) ) 
                    << " base edges in " << ToStringAddCommas(nd) 
                    << " superedges" << endl;    }

          // Try to place more reads.

          int pl1 = 0, pl2 = 0;
          for ( int64_t id = bci[b]; id < bci[b+1]; id++ )
          {    if ( dpaths[id].size( ) > 0 ) pl1++;
               const ReadPath& p = paths[id];
               if ( p.size( ) == 0 || dpaths[id].size( ) > 0 || dup[id/2] ) continue;

               // Now p = paths[id] represents a nonduplicate read from the barcode
               // that is placed on the base graph.  Find supergraph edges assigned 
               // to the barcode (as defined by locs), and find the matches of the 
               // read to one those edges that extends for at least 30 kmers; also 
               // ignore matches that are less than half as long as the best match.
               // (Not sure the "30" is enforced below.)
               //
               // First let pl =
               // { (supergraph edge d, pos on p, pos on d, match length in kmers) }.

               vec<quad<int,int,int,int>> pl;
               for ( int j = 0; j < (int) p.size( ); j++ )
               {    int e = p[j];
                    int low = LowerBound1(locs, e), high = UpperBound1(locs, e);
                    for ( int m = low; m < high; m++ )
                    {    int d = locs[m].second;
                         for ( int l = 0; l < D.O(d).isize( ); l++ )
                         {    if ( D.O(d)[l] != e ) continue;
                              if ( l > 0 && j > 0 && D.O(d)[l-1] == p[j-1] ) 
                                   continue;
                              int n = 0, r;
                              for ( r = l; r < D.O(d).isize( ); r++ )
                              {    if ( r - l + j >= (int) p.size( ) ) break;
                                   if ( D.O(d)[r] != p[r-l+j] ) break;
                                   n += hb.Kmers( D.O(d)[r] );    }
                              if ( btest >= 0 ) PRINT5_TO( out, id, d, j, l, n );
                              pl.push( d, j, l, n );    }    }    } 

               // Now filter the matches as described above.

               int nmax = 0;
               for ( int i = 0; i < pl.isize( ); i++ )
                    nmax = Max( nmax, pl[i].fourth );
               vec<Bool> to_delete( pl.size( ), False );
               for ( int i = 0; i < pl.isize( ); i++ )
               {    int n = pl[i].fourth;
                    if ( nmax >= 30 && n < nmax/2 ) to_delete[i] = True;    }
               EraseIf( pl, to_delete );
               if ( btest >= 0 ) PRINT2_TO( out, id, pl.size( ) );

               // If there is a unique placement, place the read.

               if ( pl.solo( ) )
               {    int d = pl[0].first, j = pl[0].second, l = pl[0].third;
                    vec<int> x;
                    for ( auto e : p ) x.push_back(e);
                    ho_interval xpos( j, j + 1 );
                    vec<int> q = {d};
                    int qstart = l, qstop = l + 1;
                    ( !align2 ? Align : Align2 )( 
                         D, to_left, to_right, x, xpos, q, qstart, qstop );

                    // For now, toss improper alignments.

                    if ( btest >= 0 )
                    {    PRINT6_TO( out, id, xpos.Start( ), qstart, xpos.Stop( ),
                              x.size( ), qstop );    }
                    Bool bad = False;
                    if ( xpos.Start( ) > 0 )
                    {    if ( qstart > 0 ) continue;
                         if ( qstart == 0 )
                         {    int v = to_left[ q.front( ) ];
                              if ( D.To(v).nonempty( ) ) continue;    }    }
                    if ( xpos.Stop( ) < x.isize( ) )
                    {    if ( qstop < D.O( q.back( ) ).isize( ) ) continue;
                         int w = to_right[ q.back( ) ];
                         if ( D.From(w).nonempty( ) ) continue;    }

                    // Compute new offset.

                    int pos = p.getOffset( );
                    for ( int j = 0; j < qstart; j++ )
                         pos += hb.Kmers( D.O(q[0])[j] );

                    // Check for broken placements.  Not sure how this can happen.

                    int n = hb.K( ) - 1;
                    for ( auto f : q ) n += dlens[f];
                    if ( pos > n ) continue;

                    // Save the new placement.

                    dpaths[id].resize( q.size( ) );
                    for ( int j = 0; j < q.isize( ); j++ )
                         dpaths[id][j] = q[j];
                    // out << "placing read " << id << " at pos " << pos
                    //      << ", path = " << printSeq(q) << endl;
                    dpaths[id].setOffset(pos);    
                    pl2++;    }    }
          if (verbose) 
          {    out << "orig reads placed = " << pl1 << ", new reads placed = "
                    << pl2 << "\n";
               #pragma omp critical
               {    cout << out.str( );    }    }
          else
          {
               #pragma omp critical
               {    MakeDots( stopped, ndots, bci.isize( ) - 2 );    }    }    }    }


void PlaceReads( const HyperBasevectorX& hb, const ReadPathVecX& paths, 
     const vec<Bool>& dup, const digraphE<vec<int>>& D, ReadPathVec& dpaths,
     const Bool verbose, const Bool single, const Bool align2 )
{
     // Compute dlens.  This is only used in the sanity check below.  It would be
     // good to eliminate this calculation if possible.  Then index D.

     vec<int> dlens;
     ComputeDlens( hb, D, dlens );
     vec<vec<pair<int,int>>> nd;
     MakeIndex( hb, D, nd, verbose );

     // Build paths.

     int64_t to_place = 0;
     #pragma omp parallel for reduction(+:to_place)
     for ( int64_t pid = 0; pid < dup.jsize(); pid++ ) 
     {    if ( dup[pid] ) continue;
          to_place++;    }
     // reserve space when placing a large number of paths
     if ( to_place > 10000000 && dpaths.size( ) == 0 )
     {    if (verbose) cout << Date( ) << ": defining data structure" << endl;
          // BETTER ALLOCATION NEEDED? !!!!!!!!!!!!!!!!!!!!!!!!!!
          dpaths.resize(paths.size());
          int mpl = 1;// heuristic
          for (int64_t id = 0; id<paths.size(); id++)
               dpaths[id].reserve(mpl);    }
     else{
         dpaths.resize(paths.size());
     }
     if (verbose) cout << Date( ) << ": creating to_left and to_right" << endl;
     vec<int> to_left, to_right;
     D.ToLeft(to_left), D.ToRight(to_right);
     if (verbose) cout << Date( ) << ": looking at reads" << endl;
     double clock = WallClockTime( );
     const int batch = 100000;
     int nthreads = ( single ? 1 : omp_get_max_threads( ) );
     #pragma omp parallel for num_threads(nthreads) schedule(dynamic,1)
     for ( int64_t bi = 0; bi < (int64_t) paths.size( ); bi += batch )
     {    vec<int> aplace, d, x;
          ReadPath p;
          for ( int64_t id = bi; 
               id < Min( bi + batch, (int64_t) paths.size( ) ); id++ )
          {    dpaths[id].resize(0); 
               if ( dup[id/2] ) continue;
               p.clear( );
               paths.unzip(p,hb,id);
               int nplaces, pos;
               if ( !FindPlaces( p, d, x, hb, D, dlens, to_left, to_right, 
                    nd, align2, nplaces, pos, aplace ) )
               {    continue;    }
               dpaths[id].resize( aplace.size( ) );
               for ( int j = 0; j < aplace.isize( ); j++ )
                    dpaths[id][j] = aplace[j];
               dpaths[id].setOffset(pos);    }     }
     if (verbose) cout << Date( ) << ": done, used " << TimeSince(clock) << endl;
     int64_t placed = 0;
     for ( int64_t id = 0; id < (int64_t) paths.size( ); id++ )
          if ( dpaths[id].size( ) > 0 ) placed++;
     if (verbose)
     {    cout << Date( ) << ": " 
               << PERCENT_RATIO( 3, placed, (int64_t) paths.size( ) ) 
               << " placed" << endl;    }    }

void PlaceReadsSmart( const HyperBasevectorX& hb, const ReadPathVecX& paths, 
     const vec<Bool>& dup, const digraphE<vec<int>>& D, const vec<int>& dinv,
     ReadPathVec& dpaths, const vec<vec<vec<vec<int>>>>& dlines, 
     const vec<int64_t>& bci, const Bool verbose, const int btest,
     const Bool align2 )
{
     // Set up.

     vec<int> dlens( D.E( ), 0 ), to_left, to_right, linv;
     #pragma omp parallel for
     for ( int e = 0; e < D.E( ); e++ )
     {    if ( D.O(e)[0] < 0 ) continue;
          for ( int j = 0; j < D.O(e).isize( ); j++ )
               dlens[e] += hb.Kmers( D.O(e)[j] );    }
     D.ToLeft(to_left), D.ToRight(to_right);
     LineInv( dlines, dinv, linv );

     // To each nongap edge, assign a line, and a starting position on that line.
     // For nonpalindromic lines L, we pick just one of L, rc(L), and place both
     // edges d and dinv[d] at the same place on that line.

     vec<quad<int,int,Bool,int>> tol( D.E( ), make_quad(-1,-1,True,-1) );
     cout << Date( ) << ": defining tol" << endl;
     for ( int i = 0; i < dlines.isize( ); i++ )
     {    if ( linv[i] < i ) continue;
          const vec<vec<vec<int>>>& L = dlines[i];
          int pos = 0;
          for ( int j = 0; j < L.isize( ); j++ )
          {    vec<int> lens;
               for ( int k = 0; k < L[j].isize( ); k++ )
               {    int len = 0;
                    for ( int l = 0; l < L[j][k].isize( ); l++ )
                    {    int d = L[j][k][l];
                         if ( D.O(d)[0] < 0 ) continue;
                         if ( tol[d].first < 0 )
                         {    tol[d] = make_quad( i, pos+len, True, j );
                              tol[ dinv[d] ] 
                                   = make_quad( i, pos+len, False, j );    }
                         len += dlens[d];    }
                    lens.push_back(len);    }
               Sort(lens);
               if ( lens.nonempty( ) ) pos += Median(lens);    }    }

     // Define the territory of each barcode.  The granularity of this definition
     // is weak, as we use only start/stop positions (i,j) on lines, whereas we
     // should use base positions.

     // Go through the barcodes.

     cout << Date( ) << ": start traversal" << endl;
     const int MAX_BC_GAP = 100000;
     const int MIN_BC_GROUP = 3;
     int ndots = 0, stopped = 0;
     #pragma omp parallel for schedule(dynamic, 1000)
     for ( int b = 1; b < bci.isize( ) - 1; b++ )
     {    if ( btest >= 0 && b != btest ) continue;

          // Logging.

          ostringstream out;
          if (verbose)
          {    out << "\nlooking at barcode " << b << " [" << bci[b+1]-bci[b]
                    << " reads]" << endl;    }

          // Define the territory of each barcode.  The granularity of this
          // definition is weak, as we use only start/stop positions (i,j) on lines,
          // whereas we should use base positions.

          int total = 0;
          vec<triple<int,int,int>> lplaces; // (line, start on line, line unit)
          for ( int64_t id = bci[b]; id < bci[b+1]; id++ )
          {    const ReadPath& p = dpaths[id];
               int pos = p.getOffset( );
               for ( int j = 0; j < (int) p.size( ); j++ )
               {    if ( j == 0 || tol[ p[j] ].first != tol[ p[j-1] ].first )
                    {    triple<int,int,int> x; 
                         x.first = tol[ p[j] ].first;
                         int xpos = pos;
                         if ( !tol[ p[j] ].third ) xpos = dlens[ p[j] ] - xpos;
                         x.second = tol[ p[j] ].second + xpos;
                         x.third = tol[ p[j] ].fourth;
                         lplaces.push_back(x);    }
                    pos -= dlens[ p[j] ];    }    }
          Sort(lplaces);
          vec<triple<int,int,int>> groups;
          for ( int i = 0; i < lplaces.isize( ); i++ )
          {    int j;
               for ( j = i + 1; j < lplaces.isize( ); j++ )
               {    if ( lplaces[j].first != lplaces[i].first ) break;
                    if ( lplaces[j].second - lplaces[j-1].second > MAX_BC_GAP ) 
                         break;    }
               if ( j - i >= MIN_BC_GROUP )
               {    int len = lplaces[j-1].second - lplaces[i].second;
                    int start = lplaces[i].third, stop = lplaces[j-1].third;
                    int l = lplaces[i].first;
                    const vec<vec<vec<int>>>& L = dlines[l];

                    // Extend start and stop positions.

                    const int EXT = 25000;
                    int left = 0, right = 0;
                    while( left < EXT && start >= 1 )
                    {    start--;
                         vec<int> lens;
                         for ( int k = 0; k < L[start].isize( ); k++ )
                         {    int len = 0;
                              for ( int l = 0; l < L[start][k].isize( ); l++ )
                              {    int d = L[start][k][l];
                                   len += dlens[d];    }
                              lens.push_back(len);    }
                         Sort(lens);
                         if ( lens.nonempty( ) ) left += Median(lens);    }
                    while( right < EXT && stop < L.isize( ) - 1 )
                    {    stop++;
                         vec<int> lens;
                         for ( int k = 0; k < L[stop].isize( ); k++ )
                         {    int len = 0;
                              for ( int l = 0; l < L[stop][k].isize( ); l++ )
                              {    int d = L[stop][k][l];
                                   len += dlens[d];    }
                              lens.push_back(len);    }
                         Sort(lens);
                         if ( lens.nonempty( ) ) right += Median(lens);    }
                    
                    // Save.

                    if (verbose)
                    {    out << l << "." << lplaces[i].second << "-" 
                              << lplaces[j-1].second << " [len=" << len 
                              << ",count=" << j-i 
                              << ",start=" << start << "/" << L[start][0][0]
                              << ",stop=" << stop << "/" << L[stop][0][0]
                              << "]" << endl;    }
                    total += len;
                    groups.push( lplaces[i].first, start, stop );    }
               i = j - 1;    }    
          if (verbose) out << "total = " << ToStringAddCommas(total) << endl;

          // Define the edges in the territory.

          vec<pair<int,int>> locs;
          int nd = 0;
          for ( int i = 0; i < groups.isize( ); i++ )
          {    const vec<vec<vec<int>>>& L = dlines[ groups[i].first ];
               int start = groups[i].second, stop = groups[i].third;
               for ( int j = start; j <= stop; j++ )
               {    const vec<vec<int>>& M = L[j];
                    vec<int> ds;
                    for ( int k = 0; k < M.isize( ); k++ )
                    for ( int l = 0; l < M[k].isize( ); l++ )
                    {    int d = M[k][l];
                         if ( D.O(d)[0] >= 0 ) ds.push_back( d, dinv[d] );    }
                    UniqueSort(ds);
                    nd += ds.size( );
                    for ( int k = 0; k < ds.isize( ); k++ )
                    {    int d = ds[k];
                         vec<int> es = D.O(d);
                         UniqueSort(es);
                         for ( auto e : es ) locs.push( e, d );    }    }    }
          Sort(locs);
          if (verbose) 
          {    out << "locs = " << ToStringAddCommas( locs.size( ) ) 
                    << " base edges in " << ToStringAddCommas(nd) 
                    << " superedges" << endl;    }

          // Try to place more reads.

          int pl1 = 0, pl2 = 0;
          ReadPath p;
          for ( int64_t id = bci[b]; id < bci[b+1]; id++ )
          {    if ( dpaths[id].size( ) > 0 ) pl1++;
               p.clear( );
               paths.unzip(p,hb,id);
               if ( p.size( ) == 0 || dpaths[id].size( ) > 0 || dup[id/2] ) continue;

               // Now p = paths[id] represents a nonduplicate read from the barcode
               // that is placed on the base graph.  Find supergraph edges assigned 
               // to the barcode (as defined by locs), and find the matches of the 
               // read to one those edges that extends for at least 30 kmers; also 
               // ignore matches that are less than half as long as the best match.
               // (Not sure the "30" is enforced below.)
               //
               // First let pl =
               // { (supergraph edge d, pos on p, pos on d, match length in kmers) }.

               vec<quad<int,int,int,int>> pl;
               for ( int j = 0; j < (int) p.size( ); j++ )
               {    int e = p[j];
                    int low = LowerBound1(locs, e), high = UpperBound1(locs, e);
                    for ( int m = low; m < high; m++ )
                    {    int d = locs[m].second;
                         for ( int l = 0; l < D.O(d).isize( ); l++ )
                         {    if ( D.O(d)[l] != e ) continue;
                              if ( l > 0 && j > 0 && D.O(d)[l-1] == p[j-1] ) 
                                   continue;
                              int n = 0, r;
                              for ( r = l; r < D.O(d).isize( ); r++ )
                              {    if ( r - l + j >= (int) p.size( ) ) break;
                                   if ( D.O(d)[r] != p[r-l+j] ) break;
                                   n += hb.Kmers( D.O(d)[r] );    }
                              if ( btest >= 0 ) PRINT5_TO( out, id, d, j, l, n );
                              pl.push( d, j, l, n );    }    }    } 

               // Now filter the matches as described above.

               int nmax = 0;
               for ( int i = 0; i < pl.isize( ); i++ )
                    nmax = Max( nmax, pl[i].fourth );
               vec<Bool> to_delete( pl.size( ), False );
               for ( int i = 0; i < pl.isize( ); i++ )
               {    int n = pl[i].fourth;
                    if ( nmax >= 30 && n < nmax/2 ) to_delete[i] = True;    }
               EraseIf( pl, to_delete );
               if ( btest >= 0 ) PRINT2_TO( out, id, pl.size( ) );

               // If there is a unique placement, place the read.

               if ( pl.solo( ) )
               {    int d = pl[0].first, j = pl[0].second, l = pl[0].third;
                    vec<int> x;
                    for ( auto e : p ) x.push_back(e);
                    ho_interval xpos( j, j + 1 );
                    vec<int> q = {d};
                    int qstart = l, qstop = l + 1;
                    ( !align2 ? Align : Align2 )( 
                         D, to_left, to_right, x, xpos, q, qstart, qstop );

                    // For now, toss improper alignments.

                    if ( btest >= 0 )
                    {    PRINT6_TO( out, id, xpos.Start( ), qstart, xpos.Stop( ),
                              x.size( ), qstop );    }
                    Bool bad = False;
                    if ( xpos.Start( ) > 0 )
                    {    if ( qstart > 0 ) continue;
                         if ( qstart == 0 )
                         {    int v = to_left[ q.front( ) ];
                              if ( D.To(v).nonempty( ) ) continue;    }    }
                    if ( xpos.Stop( ) < x.isize( ) )
                    {    if ( qstop < D.O( q.back( ) ).isize( ) ) continue;
                         int w = to_right[ q.back( ) ];
                         if ( D.From(w).nonempty( ) ) continue;    }

                    // Compute new offset.

                    int pos = p.getOffset( );
                    for ( int j = 0; j < qstart; j++ )
                         pos += hb.Kmers( D.O(q[0])[j] );

                    // Check for broken placements.  Not sure how this can happen.

                    int n = hb.K( ) - 1;
                    for ( auto f : q ) n += dlens[f];
                    if ( pos > n ) continue;

                    // Save the new placement.

                    dpaths[id].resize( q.size( ) );
                    for ( int j = 0; j < q.isize( ); j++ )
                         dpaths[id][j] = q[j];
                    // out << "placing read " << id << " at pos " << pos
                    //      << ", path = " << printSeq(q) << endl;
                    dpaths[id].setOffset(pos);    
                    pl2++;    }    }
          if (verbose) 
          {    out << "orig reads placed = " << pl1 << ", new reads placed = "
                    << pl2 << "\n";
               #pragma omp critical
               {    cout << out.str( );    }    }
          else
          {
               #pragma omp critical
               {    MakeDots( stopped, ndots, bci.isize( ) - 2 );    }    }    }    }

// Alternative way to make the index for read placement, it's faster
// but has a bigger footprint

void MakeIndexAlt( const HyperBasevectorX& hb, const digraphE<vec<int>>& D,
     vec<vec<pair<int,int>>> & nd, const Bool verbose )
{    if (verbose) cout << Date( ) << ": creating index for read placement" << endl;
     vec<Bool> used;
     D.Used(used);
     vec<triple<int,int,int>> edj;
     for ( int d = 0; d < D.E( ); d++ )
     {    if ( !used[d] ) continue;
          const vec<int>& x = D.O(d);
          if ( x[0] < 0 ) continue;
          for ( int j = 0; j < x.isize( ); j++ )
               edj.push( x[j], d, j );    }    
     ParallelSort( edj );
     int prev = -1;
     for ( int64_t i = 0; i != edj.jsize( ); i++ ) {
          const auto & t = edj[i];
          if ( t.first != prev ) {
               nd.append( vec<vec<pair<int,int>>> ( t.first - prev ) );
          }
          nd[t.first].push ( t.second, t.third );
          prev = t.first;
     }
     if ( prev != hb.E( )-1 )
          nd.append( vec<vec<pair<int,int>>> ( hb.E( ) - 1  - prev ) );
     ForceAssertEq( hb.E( ), nd.isize( ) );
}

void ExploreEdgesWithinDistanceDepth ( const digraphE<vec<int>> & D, 
     const vec<int> & to_left, const vec<int> & to_right,
     const vec<int> & dlens, const int & v, vec<int> & expl, 
     const int & DIST, const uint8_t & DEPTH, const Bool & RIGHT )
{
     set<int> all = {v};
     set<int> vset = {v};
     int dist = 0;
     for ( uint8_t depth = 0; depth < DEPTH; depth++ ) {
          const int start = expl.size( );
          int d = 0;
          for ( auto & vp : vset ) {
               const int imax = ( RIGHT ? D.IFrom(vp).size() : D.ITo(vp).size() );
               for ( int i = 0; i != imax; i++ ) {
                    const int edge = (RIGHT ? D.IFrom( vp, i ) : D.ITo( vp, i ));
                    expl.push_back( edge );
                    d = Max( dlens[edge], d );
               }
          }
          dist += d;
          if ( dist > DIST )
               break;
          vset.clear( );
          for ( int i = start; i != expl.isize( ); i++ ) {
               const int w = (RIGHT ? to_right[expl[i]] : to_left[expl[i]]); 
               if ( all.count( w ) ) continue;
               vset.insert( w );
               all.insert( w );
          }
     }
     UniqueSort( expl );
}

// This is the same code as FindPlaces except that it allows for multiple
// placements.

void FindAllPlacements( const ReadPath& p, const HyperBasevectorX& hb, 
     const digraphE<vec<int>>& D, const vec<int>& dlens,
     const vec<int>& to_left, const vec<int>& to_right,
     const vec<vec<pair<int,int>>>& nd, const Bool align2,
     int& nplaces, vec<pair<int,vec<int>>> & placements )
{
     nplaces=0;

     // Empty path?
     
     if ( p.size( ) == 0 ) return;

     // Find minimum multiplicity edge in p.

     int m = 1000000000, mp = -1;
     for ( int j = 0; j < (int) p.size( ); j++ )
     {    int n = nd[ p[j] ].size( );
          if ( n == 0 ) continue;
          if ( n < m )
          {    m = n;
               mp = j;    }    }

     // Test placements.

     if ( mp < 0 || nd[ p[mp] ].empty( ) ) return;
     int pos = -1;
     placements.clear( );
     vec<int> x,d;
     for ( int i = 0; i < nd[ p[mp] ].isize( ); i++ )
     {    int e = nd[ p[mp] ][i].first, epos = nd[ p[mp] ][i].second;
          x.clear( );
          for ( int j = 0; j < (int) p.size( ); j++ ) x.push_back( p[j] );
          ho_interval xpos( mp, mp+1 );
          d = {e};
          int dstart = epos, dstop = epos + 1;
          ( !align2 ? Align : Align2 )(
               D, to_left, to_right, x, xpos, d, dstart, dstop );

          // For now, toss improper alignments.

          if ( xpos.Start( ) > 0 )
          {    if ( dstart > 0 ) continue;
               if ( dstart == 0 )
               {    int v = to_left[ d.front( ) ];
                    if ( D.To(v).nonempty( ) ) continue;    }    }
          if ( xpos.Stop( ) < x.isize( ) )
          {    if ( dstop < D.O( d.back( ) ).isize( ) ) continue;
               int w = to_right[ d.back( ) ];
               if ( D.From(w).nonempty( ) ) continue;    }

          // Compute start position of read on d[0].

          pos = p.getOffset( );
          for ( int j = 0; j < dstart; j++ )
               pos += hb.Kmers( D.O(d[0])[j] );
          
          // Check for broken placement.  Need to track down how this could ever happen.

          int n = hb.K( ) - 1;
          for ( auto f : d ) n += dlens[f];
          if ( pos > n ) continue;
          
          nplaces++;
          placements.push_back( make_pair( pos, d ) );
     }
}


// Linked-read placement
// - First place reads uniquely
// - Then place read partners: if a read is uniquely placed, but its partner
// has multiple potential placements, then check if there is a unique placement
// in a 1 kb neighborhood of the uniquely placed read. This is done by
// exploring the graph up to a fixed distance/depth.
// - Within a barcode take all the uniquely placed reads ( and partners )
// and mark a zone of 50 kb around each one. If a read in that barcode
// has multiple placement, but only ONE placement inside such a zone, then
// choose that one. The 50 kb zones are half-open intervals on edges in the 
// line graph

void PlaceLinkedReads( const HyperBasevectorX & hb, const vec<int> & inv,
     const digraphE<vec<int>> & D, const vec<int> & dinv, 
     const vec<vec<vec<vec<int>>>> & dlines, const vec<Bool> & dup,
     const vec<int64_t> & bci, const ReadPathVecX & pathsx, ReadPathVec & dpaths,
     const int MAX_PASSES, const Bool verbose )
{
     // Define Heuristics
     // Which alignment code to use (we never use ALIGN2)
     const Bool ALIGN2 = False;
     // How far (in edges) to explore the graph to look for a partner
     const uint8_t DEPTH = 10;
     // How far to explore the line graph
     const int MAX_DEPTH = 10;
     // How far in bases to explore the graph
     const int DIST = 1000;
     // What is the max separation of linked reads within a molecule
     const int MAX_BC_SEP = 30000;
     
     ForceAssertGt( MAX_PASSES, 0 );
     if ( verbose ) cout << Date( ) << ": setting up data structures" << endl;
     vec<int> to_left, to_right, dlens(D.E( ), 0);
     D.ToLeft( to_left );
     D.ToRight( to_right );
     
     // Compute super edge lengths
     
     #pragma omp parallel for
     for ( int d = 0; d < D.E( ); d++ ) {
          if ( D.O(d)[0] < 0 ) continue;
          int & len = dlens[d];
          for ( auto & e : D.O(d) )
               len += hb.Kmers(e);
     }

     // Compute index for read placement
     
     vec<vec<pair<int,int>>> nd;
     MakeIndex( hb, D, nd, True );

     // Compute line lengths (using median cell path length)
     // and compute super edge location along line

     if ( verbose ) cout << Date( ) << ": computing edge locations on lines" << endl;
     vec<int> llens(dlines.size( ), 0);
     vec<pair<int,int>> dloc(D.E( ));
     for ( int i = 0; i < dlines.isize(); i++ ) {
          const auto & L = dlines[i];
          int pos = 0;
          for ( int j = 0; j < L.isize( ); j++ ) {
               vec<int> plens;
               for ( auto & path : L[j] ) {
                    int plen = 0;
                    for ( auto & d : path ) {
                         dloc[d] = make_pair( i, plen + pos );
                         plen += dlens[d];
                    }
                    plens.push_back( plen );     
               }
               pos += Median( plens );
          }
          llens[i] = pos;
     }
     
     if ( verbose ) cout << Date( ) << ": set up line graph" << endl;
     // Compute involution on lines

     vec<int> linv;
     LineInv( dlines, dinv, linv );
     
     // Build line graph
     
     digraphE<int> LG;
     BuildLineGraph( D, dinv, dlines, LG );
     
     vec<int> lg_to_left, lg_to_right;
     LG.ToLeft( lg_to_left );
     LG.ToRight( lg_to_right );

     // Clear dpaths and reserve space if first time placing reads

     size_t N = pathsx.size();
     if ( dpaths.size() != N ) {
          cout << Date( ) << ": reserving space (slow)" << endl;
          dpaths.clear( );
          dpaths.resize( N );
          for ( size_t id = 0; id != N; id++ ) {
               if ( dup[id/2] ) continue;
               dpaths[id].reserve( 2 );
          }
     } else {
          #pragma omp parallel for
          for ( size_t id = 0; id < N; id++ )
               dpaths[id].clear();
     }

     const Bool PAIR_EXPLORE_RIGHT = True;
     int64_t tot_placed = 0;
     vec<vec<pair<int,vec<int>>>> all_placements;
     if ( verbose ) cout << Date( ) << ": begin read placement" << endl;
     
     // First place barcode 0 reads
     
     if ( verbose ) cout << Date( ) << ": placing un-barcoded reads" << endl;

     #pragma omp parallel for schedule( dynamic, 1 ) reduction (+:tot_placed) \
          firstprivate( all_placements )
     for ( int64_t id1 = bci[0]; id1 < bci[1]; id1 += 2 ) {
          if ( dup[id1/2] ) continue;

          const int64_t id2 = id1 + 1;
          all_placements.clear( );
          all_placements.resize( 2 );
          // try to place the individual reads uniquely
          for ( int pi = 0; pi < 2; pi++ ) {
               const int64_t id = id1+pi;
               ReadPath & p = dpaths[id];
               
               ReadPath bp;
               pathsx.unzip( bp, hb, id );
               if ( bp.size() == 0 ) continue;

               vec<pair<int,vec<int>>> & placements = all_placements[pi];
               placements.clear( );
               int nplaces = 0;
               FindAllPlacements( bp, hb, D, dlens, to_left, to_right, nd, ALIGN2,
                    nplaces, placements );
               
               if ( nplaces == 1 ) {
                    p.setOffset( placements[0].first );
                    for ( auto  & d : placements[0].second )
                         p.push_back( d );
                    tot_placed++;
               }
          }

          // if we placed exactly one read of the pair
          // then try to place the partner

          if ( dpaths[id1].empty( ) != dpaths[id2].empty( ) ) {
               const int64_t placed   = ( dpaths[id1].empty() ? id2 : id1 );
               const int64_t unplaced = ( placed % 2 == 0 ? placed+1 : placed-1 );
               const auto & placements = all_placements[unplaced-id1];
               
               if ( placements.nonempty() ) {
                    ForceAssert( placements.size( ) != 1 );
                    // p1 has already been placed
                    // p2 is empty
                    const ReadPath & p1 = dpaths[placed];
                    ReadPath & p2 = dpaths[unplaced];
                    
                    const int v = to_right[p1.back( )];
                    vec<int> expl;
                    for ( auto & d : p1 )
                         expl.push_back(d);
                    ExploreEdgesWithinDistanceDepth( D, to_left, to_right, dlens, v, 
                         expl, DIST, DEPTH, PAIR_EXPLORE_RIGHT );

                    vec<int> hits;
                    for ( int i = 0; i != placements.isize(); i++ ) {
                         const int d = dinv[placements[i].second.back()];
                         if ( BinMember( expl, d ) )
                              hits.push_back(i);
                    }
                    if ( hits.size() == 1 ) {
                         p2.setOffset( placements[hits[0]].first );
                         for ( auto & d : placements[hits[0]].second )
                              p2.push_back( d );
                         tot_placed++;
                    } 
               }
          }
     }
     
     // Now place the barcoded reads in multiple passes

     if ( verbose ) cout << Date( ) << ": placing barcoded reads" << endl;
     all_placements.clear( );
     map<int,ho_interval> lzone, lzone2, new_lines;

     int pass = 0;
     while ( pass < MAX_PASSES ) {
     pass++;
     if ( verbose )
          cout << Date( ) << ": pass " << pass << " of " << MAX_PASSES << endl;
     all_placements.clear( );
     lzone.clear( );
     lzone2.clear( );
     #pragma omp parallel for schedule( dynamic, 1 ) reduction (+:tot_placed) \
          firstprivate( all_placements,lzone,lzone2,new_lines )
     for ( int b = 1; b < bci.isize()-1; b++ ) {
          const int64_t start = bci[b], stop = bci[b+1];
          
          all_placements.clear( );
          lzone.clear( );
          lzone2.clear( );
          all_placements.resize( stop-start );

          vec<int> expl;
          for ( int64_t id1 = start; id1 < stop; id1 += 2 ) {
               if ( dup[id1/2] ) continue;

               const int64_t id2 = id1 + 1;
               // try to place the individual reads uniquely
               for ( int pi = 0; pi < 2; pi++ ) {
                    const int64_t id = id1+pi;
                    ReadPath & p = dpaths[id];
                    
                    // if we already have a path, it means we placed the read
                    // in a previous iteration, so skip

                    if ( p.size() > 0 )
                         continue;
                         
                    ReadPath bp;
                    pathsx.unzip( bp, hb, id );
                    if ( bp.size() == 0 ) continue;

                    vec<pair<int,vec<int>>> & placements = all_placements[id-start];
                    placements.clear( );
                    int nplaces = 0;
                    FindAllPlacements( bp, hb, D, dlens, to_left, to_right, nd, ALIGN2,
                         nplaces, placements );
                    
                    if ( nplaces == 1 ) {
                         p.setOffset( placements[0].first );
                         for ( auto  & d : placements[0].second )
                              p.push_back( d );
                         tot_placed++;
                    }
               }

               // if we placed exactly one read of the pair
               // then try to place the partner

               if ( dpaths[id1].empty( ) != dpaths[id2].empty( ) ) {
                    const int64_t placed   = ( dpaths[id1].empty() ? id2 : id1 );
                    const int64_t unplaced = ( placed % 2 == 0 ? placed+1 : placed-1 );
                    const auto & placements = all_placements[unplaced-start];
                    
                    if ( placements.nonempty() ) {
                         // p1 has already been placed
                         // p2 is empty
                         ReadPath & p1 = dpaths[placed];
                         ReadPath & p2 = dpaths[unplaced];
                         
                         const int v = to_right[p1.back( )];
                         expl.clear( );
                         for ( auto & d : p1 )
                              expl.push_back(d);
                         ExploreEdgesWithinDistanceDepth( D, to_left, to_right, dlens, v, 
                              expl, DIST, DEPTH, PAIR_EXPLORE_RIGHT );
                         
                         vec<int> hits;
                         for ( int i = 0; i != placements.isize(); i++ ) {
                              const int d = dinv[placements[i].second.back()];
                              if ( BinMember( expl, d ) )
                                   hits.push_back(i);
                         }
                         if ( hits.size() == 1 ) {
                              p2.setOffset( placements[hits[0]].first );
                              for ( auto & d : placements[hits[0]].second )
                                   p2.push_back( d );
                              tot_placed++;
                         }
                    }
               }
               
               // now take the uniquely placed reads in the pair
               // and create a zone around them defined by the barcode

               for ( int pi = 0; pi < 2; pi++ ) {
                    const int64_t id = id1+pi;
                    if ( dpaths[id].size() > 0 ) {
                         const ReadPath & p = dpaths[id];
                         const int l = dloc[p[0]].first;
                         const int loffset = p.getOffset( ) + dloc[p[0]].second;
                         ho_interval bczone(  loffset-MAX_BC_SEP,
                                              loffset+MAX_BC_SEP ); 
                         if ( lzone.count( l ) == 0 )
                              lzone[l]       = bczone;
                         else {
                              auto & linfo = lzone[l];
                              linfo  = Span( linfo, bczone );
                         }
                    }
               }
          }
          
          // Extend this zone by following the line graph

          for ( int depth = 0; depth < MAX_DEPTH; depth++ ) {
               new_lines.clear( );
               int num_ext = 0;
               for ( auto & it : lzone ) {
                    const int l = it.first;
                    auto & z = it.second;
                    int & start = z.StartMutable( ), & stop = z.StopMutable( );

                    if ( start < 0 ) {
                         num_ext++;
                         int dist = -start;
                         start = 0;
                         const int v = lg_to_left[l];
                         for ( int i = 0; i != LG.To(v).isize( ); i++ ) {
                              const int lp = LG.ITo(v, i);
                              if ( lzone.count( lp ) )
                                   continue;
                              ho_interval lp_zone ( llens[lp] - dist, llens[lp] );
                              if ( new_lines.count( lp ) ) {
                                   auto & lp_info = new_lines[lp];
                                   lp_info = Span( lp_info, lp_zone );
                              } else
                                   new_lines[lp] = lp_zone;
                         }
                    }
                    if ( stop > llens[l] ) {
                         num_ext++;
                         int dist = stop-llens[l];
                         stop=llens[l];
                         const int v = lg_to_right[l];
                         for ( int i = 0; i != LG.From(v).isize( ); i++ ) {
                              const int lp = LG.IFrom(v, i);
                              if ( lzone.count(lp) )
                                   continue;
                              ho_interval lp_zone ( 0, dist );
                              if ( new_lines.count( lp ) ) {
                                   auto & lp_info = new_lines[lp];
                                   lp_info = Span( lp_info, lp_zone );
                              } else
                                   new_lines[lp] = lp_zone;
                         }
                    }
               }
               if ( num_ext == 0 )
                    break;
               // add in new lines that we found
               
               for ( auto & it : new_lines )
                    lzone[it.first] = it.second;
          }
          
          // make it closed under the involution
          
          for ( auto & it1 : lzone ) {
               const int l = it1.first;
               const int rl = linv[l];
               ho_interval & i1 = it1.second;
               ho_interval i2( llens[l] - i1.Stop( ), llens[l] - i1.Start( ) );
               if ( lzone.count( rl ) == 0 ) {
                    lzone2[l] = i1;
                    lzone2[rl] = i2;
               } else {
                    ho_interval j2 = Span(lzone[rl],i2);
                    ho_interval j1 = Span( i1, ho_interval( llens[l]-j2.Stop( ), 
                                                            llens[l]-j2.Start( ) ));
                    lzone2[l]=j1;
                    lzone2[rl]=j2;
               }
          }

          // Now place reads with multiple placements

          for ( int64_t id = start; id < stop; id++ ) {
               if ( dup[id/2] ) continue;
               if ( dpaths[id].size( ) > 0 ) continue;
               
               const vec<pair<int,vec<int>>> & placements = all_placements[id-start];
               if ( placements.empty() ) continue;
               int hit = -1;
               for ( int pli = 0; pli != placements.isize(); pli++ ) {
                    const auto & pl = placements[pli];
                    const vec<int> & dp = pl.second;
                    ForceAssertGt( dp.size(), 0 );
                    const int l = dloc[dp[0]].first;
                    const int loffset = pl.first + dloc[dp[0]].second;
                    if ( lzone2.count( l ) == 0 ) continue;
                    const auto & zone = lzone2[l];
                    if ( zone.Contains( loffset ) ) {
                         if ( hit == -1 )
                              hit = pli;
                         else {
                              hit = -2;
                              break;
                         }
                    }
               }

               if ( hit >= 0 ) {
                    const auto & pl = placements[hit];
                    const auto & dp = pl.second;
                    ForceAssertGt( dp.size(), 0 );
                    ReadPath & p = dpaths[id];
                    p.setOffset( pl.first );
                    for ( auto & d : dp )
                         p.push_back( d );
                    tot_placed++;
               }
          }
     } 
          if ( verbose )
               cout << Date( ) << ": placed " << ToStringAddCommas(tot_placed)
                    << " reads" << endl;
          if ( tot_placed == 0 )
               break;
     } // END multiple passes

     if ( verbose )
          cout << Date( ) << ": placed " 
               << PERCENT_RATIO( 3, tot_placed, (int64_t) dpaths.size( ) ) << " of reads"
               << endl;

}
