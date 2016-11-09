// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"
#include "paths/HyperBasevector.h"
#include "10X/Gap.h"
#include "10X/Super.h"

void ReinsertLoop( const int d, const HyperBasevectorX& hb, const vec<int>& inv,
     digraphE<vec<int>>& D, vec<int>& dinv, vec<int>& to_left, vec<int>& to_right )
{    int rd = dinv[d];
     const vec<int>& x = D.O(d);
     int v = to_left[d], w = to_right[d];
     int rv = to_right[rd], rw = to_left[rd];
     cell c, cr;
     c.CellDecode(x);
     const vec<int>& rx = D.O(rd);
     cr.CellDecode(rx);
     const digraphE<vec<int>> &G = c.G( ), &GR = cr.G( );
     int N = D.N( ), E = D.E( ), n = G.N( ), e = G.E( );
     D.AppendWithUpdate( G, to_left, to_right );
     D.AppendWithUpdate( GR, to_left, to_right );
     if ( rd == d )
     {    // TOTALLY NOT RIGHT FOR THIS CASE
          cout << "WRONG" << endl;
          Scram(0);
          for ( int f = 0; f < e; f++ )
               dinv.push_back( E + f );
          for ( int f = 0; f < e; f++ )
               dinv.push_back( E + e + f );    }
     else
     {    for ( int f = 0; f < e; f++ ) dinv.push_back( E + e + f );
          for ( int f = 0; f < e; f++ ) dinv.push_back( E + f );    }
     D.TransferEdgesWithUpdate( N + c.Left( ), v, to_left, to_right );
     D.TransferEdgesWithUpdate( N + n + cr.Left( ), rw, to_left, to_right );
     if ( c.Left( ) == c.Right( ) )
     {    if ( w != v )
          {    D.TransferEdgesWithUpdate( w, v, to_left, to_right );
               D.TransferEdgesWithUpdate( rv, rw, to_left, to_right );    }    }
     else
     {    D.TransferEdgesWithUpdate( N + c.Right( ), w, to_left, to_right );
          D.TransferEdgesWithUpdate( 
               N + n + cr.Right( ), rv, to_left, to_right );    }
     D.DeleteEdgesWithUpdate( {d,rd}, to_left, to_right );    }


void ReinsertLoopsMap( digraphE<vec<int>> const& D, vec<int> const& dinv, 
          std::map<int,int>& edge_map )
{
     vec<int> to_left, to_right;
     D.ToLeft(to_left), D.ToRight(to_right);
     int nd = D.E( );
     for ( int d = 0; d < D.E(); d++ )
     {    int rd = dinv[d];
          if ( rd < d ) continue;
          if ( rd == d ) continue; // punting on this case (lazy)
          if ( IsCell( D.O(d) ) ) {    
               int rd = dinv[d];
               const vec<int>& x = D.O(d);
               const vec<int>& rx = D.O(rd);
               cell c, cr;
               c.CellDecode(x);
               cr.CellDecode(rx);
               const digraphE<vec<int>> &G = c.G( ), &GR = cr.G( );
               edge_map[d] = nd;
               nd += G.E();
               edge_map[rd] = nd;
               nd += GR.E();
          }
     }
}

void ReinsertLoops( const HyperBasevectorX& hb, const vec<int>& inv,
     digraphE<vec<int>>& D, vec<int>& dinv )
{
     // please keep ReinsertLoopsMap above in sync with
     // changes below
     cout << Date( ) << ": reinserting loops" << endl;
     vec<int> to_left, to_right;
     D.ToLeft(to_left), D.ToRight(to_right);
     int nd = D.E( );
     for ( int d = 0; d < nd; d++ )
     {    int rd = dinv[d];
          if ( rd < d ) continue;
          if ( rd == d ) continue; // punting on this case (lazy)
          if ( IsCell( D.O(d) ) )
          {    ReinsertLoop( d, hb, inv, D, dinv, to_left, to_right );    }    }
     Validate( hb, inv, D, dinv );    }

// this is a copy of ReinsertLoops above that is to be used by
// ValidateMakeFasta. It enables tracking of cell index.
void ReinsertLoopsWithTracking( const HyperBasevectorX& hb, const vec<int>& inv,
     digraphE<vec<int>>& D, vec<int>& dinv , vec <int> & ctracker )
{
     // please keep ReinsertLoopsMap above in sync with
     // changes below
     cout << Date( ) << ": reinserting loops (with tracking)" << endl;
     vec<int> to_left, to_right;
     D.ToLeft(to_left), D.ToRight(to_right);
     int nd = D.E( );
     ctracker.clear();
     ctracker.resize( D.E( ), -1 );
     int cell_index = 0;
     for ( int d = 0; d < nd; d++ )
     {    int rd = dinv[d];
          if ( rd < d ) continue;
          if ( rd == d ) continue; // punting on this case (lazy)
          const int ne = D.E( );
          if ( IsCell( D.O(d) ) )
          {    ReinsertLoop( d, hb, inv, D, dinv, to_left, to_right );    }
          // track new edges
          for ( int j = 0; j < (D.E( )-ne); j++ )
               ctracker.push_back( cell_index );
          cell_index++;
          }
     Validate( hb, inv, D, dinv );    }


// Encode: Encode a vec<int> as a vec<int>.  Appends.
// Decode: undo it.

void Encode( const vec<int>& x, vec<int>& y )
{    y.push_back( x.size( ) );
     y.append(x);    }
void Decode( const vec<int>& y, int& pos, vec<int>& x )
{    int n = y[pos++];
     x.reserve(n);
     for ( int j = 0; j < n; j++ ) x.push_back( y[pos++] );    }

// Encode: Encode a vec<vec<int>> as a vec<int>.  Appends.
// Decode: undo it.

void Encode( const vec<vec<int>>& x, vec<int>& y )
{    y.push_back( x.size( ) );
     for ( int i = 0; i < x.isize( ); i++ )
     {    y.push_back( x[i].size( ) );
          y.append( x[i] );    }    }
void Decode( const vec<int>& y, int& pos, vec<vec<int>>& x )
{    int n = y[pos++];
     x.resize(n);
     for ( int j = 0; j < n; j++ )
     {    int k = y[pos++];
          x[j].reserve(k);
          for ( int i = 0; i < k; i++ ) x[j].push_back( y[pos++] );    }    }
               
// Encode a digraphE<vec<int>> as a vec<int>.  Appends.

void Encode( const digraphE<vec<int>>& D, vec<int>& x )
{    Encode( D.From( ), x );
     Encode( D.To( ), x );
     Encode( D.FromEdgeObj( ), x );
     Encode( D.ToEdgeObj( ), x );
     Encode( D.Edges( ), x );    }

// Undo encode.

void Decode( const vec<int>& x, int& pos, digraphE<vec<int>>& D )
{    D.Clear( );
     Decode( x, pos, D.FromMutable( ) );
     Decode( x, pos, D.ToMutable( ) );
     Decode( x, pos, D.FromEdgeObjMutable( ) );
     Decode( x, pos, D.ToEdgeObjMutable( ) );
     Decode( x, pos, D.EdgesMutable( ) );    }

void cell::CellEncode( vec<int>& x ) const
{    x.push_back(-4);
     Encode( G_, x );
     x.push_back( v_, w_ );    }

void cell::CellDecode( const vec<int>& x )
{    int pos = 1;
     Decode( x, pos, G_ );
     v_ = x[pos++];
     w_ = x[pos++];    }

void SeqToGap( const int ltrim, const int rtrim, const basevector& x, vec<int>& y )
{    int n = x.size( );
     ForceAssert( n > 0 );
     y.clear( );
     y.resize( 4 + (n+15)/16, 0 );
     ForceAssertGe( ltrim, 0 );
     ForceAssertGe( rtrim, 0 );
     y[0] = -3, y[1] = ltrim, y[2] = rtrim, y[3] = n;
     for ( int j = 0; j < n; j++ )
     {    y[ 4 + j/16 ] ^= x[j] << (2*(j%16));    }    }

void GapToSeq( const vec<int>& y, int& ltrim, int& rtrim, basevector& x )
{    ForceAssertGe( y.isize( ), 5 );
     ForceAssertEq( y[0], -3 );
     ltrim = y[1];
     rtrim = y[2];
     int n = y[3];
     ForceAssert( n > 0 );
     x.resize(n);
     ForceAssertEq( y.isize( ), 4 + (n+15)/16 );
     for ( int j = 0; j < n; j++ )
          x.Set( j, ( y[ 4 + j/16 ] >> (2*(j%16)) ) & 3 );    }

void cell::FindPath( vec<int>& p ) const
{    p.clear( );
     vec<int> to_left, to_right;
     G( ).ToLeft(to_left), G( ).ToRight(to_right);
     vec<vec<int>> paths;
     const int MAX_COPIES = 2;
     const int MAX_PATHS = 1000;
     const int MAX_ITERATIONS = 10000;
     Bool ok = G( ).EdgePaths( to_left, to_right, Left( ), Right( ), paths, 
          MAX_COPIES, MAX_PATHS, MAX_ITERATIONS, True );
     if ( !ok || paths.empty( ) ) return;
     vec< pair<int,int> > score( paths.size( ) );
     for ( int i = 0; i < paths.isize( ); i++ )
          score[i] = make_pair( -Contents( paths[i] ).size( ), paths[i].size( ) );
     SortSync( score, paths );
     p = paths[0];    }

void ValidateGapEdges( const HyperBasevectorX& hb, const vec<int>& inv,
     const digraphE<vec<int>>& D, const vec<int>& dinv )
{    vec<int> to_left, to_right;
     D.ToLeft(to_left), D.ToRight(to_right);
     for ( int d1 = 0; d1 < D.E( ); d1++ )
     {    int d2 = dinv[d1];
          if ( d2 < d1 ) continue;
          const vec<int> &x1 = D.O(d1), &x2 = D.O(d2);
          if ( x1[0] >= 0 ) continue;
          if ( x1[0] != x2[0] )
          {    cout << "ValidateGapEdges: gap type inequality." << endl;
               Scram(1);    }
          if ( IsSequence(x1) )
          {    int ltrim1, rtrim1, ltrim2, rtrim2;
               basevector b1, b2;
               GapToSeq( x1, ltrim1, rtrim1, b1 );
               GapToSeq( x2, ltrim2, rtrim2, b2 );
               if ( ltrim1 != rtrim2 || rtrim1 != ltrim2 )
               {    cout << "ValidateGapEdges: illegal trimming." << endl;
                    Scram(1);    }
               b2.ReverseComplement( );
               if ( b1 != b2 )
               {    cout << "ValidateGapEdges: rc failed." << endl;
                    Scram(1);    }    }
          if ( IsCell(x1) )
          {    cell c1, c2;
               c1.CellDecode(x1), c2.CellDecode(x2);
               const digraphE<vec<int>>& G1 = c1.G( );
               for ( int g = 0; g < G1.E( ); g++ )
               {    if ( IsCell( G1.O(g) ) )
                    {    cout << "ValidateGapEdges: illegal embedded cell in "
                              << "gap cell edge " << d1 << "." << endl;
                         Scram(1);    }
                    // should do sanity test on seq gap edges 
                         }
               digraphE<vec<int>> G2 = c2.G( );
               G2.Reverse( );
               for ( int g = 0; g < G2.E( ); g++ )
               {    if ( IsPairGap( G2.O(g) ) );
                    else if ( IsBarcodeOnlyGap( G2.O(g) ) );
                    else if ( IsSequence( G2.O(g) ) )
                    {    int ltrim, rtrim;
                         basevector b;
                         GapToSeq( G2.O(g), ltrim, rtrim, b );
                         b.ReverseComplement( );
                         vec<int> x;
                         SeqToGap( rtrim, ltrim, b, x );
                         G2.OMutable(g) = x;    }
                    else
                    {    G2.OMutable(g).ReverseMe( );
                         for ( int j = 0; j < G2.O(g).isize( ); j++ )
                         {    int e = G2.O(g)[j];
                              if ( e < 0 || e >= hb.E( ) )
                              {    cout << "ValidateGapEdges: illegal base edge."
                                        << endl;
                                   Scram(1);    }
                              G2.OMutable(g)[j] = inv[e];    }    }    }
               // Note we should have G1 = G2, except for the possibility that
               // vertices have been permuted.  We could test to see if G1 and G2
               // are equal (up to vertex permutation), but don't.
               if ( G1.E( ) != G2.E( ) )
               {    cout << "ValidateGapEdges: asymmetry in gap cells." << endl;
                    cout << "Edge counts differ: " << G1.E( ) << " vs "
                         << G2.E( ) << endl;
                    Scram(1);    }
               for ( int g = 0; g < G1.E( ); g++ )
               {    if ( G1.O(g) != G2.O(g) )
                    {    cout << "ValidateGapEdges: gap cell asymmetry." << endl;
                         Scram(1);    }    }    }    }    }

void Munch( digraphE<vec<int>>& D, const vec<int>& dinv, vecbasevector& tigs, 
     vec<int>& inv, const int HBK )
{    
     // Create ancillary data structures.

     vec<int> to_left, to_right, dlens( D.E( ), 0 );
     D.ToLeft(to_left), D.ToRight(to_right);
     #pragma omp parallel for
     for ( int d = 0; d < D.E( ); d++ )
     {    if ( D.O(d)[0] < 0 ) continue;
          for ( int j = 0; j < D.O(d).isize( ); j++ )
               dlens[d] += tigs[ D.O(d)[j] ].isize( ) - HBK + 1;    }

     // Find vertices to the left of sequence gaps.

     vec<int> seqverts;
     #pragma omp parallel for
     for ( int d = 0; d < D.E( ); d++ )
     {    if ( IsSequence( D.O(d) ) )
          {
               #pragma omp critical
               {    seqverts.push_back( to_left[d] );    }    }    }
     UniqueSort(seqverts);

     // Process sequence gaps.

     int nprobs = 0;
     for ( int i = 0; i < seqverts.isize( ); i++ )
     {    int v = seqverts[i];
          int ltrim, rtrim;
          basevector b;

          // Check if already processed.

          int d = D.IFrom( v, 0 );
          if ( !IsSequence( D.O(d) ) ) continue;

          // Test for legal problems.

          int w = to_right[d];
          GapToSeq( D.O(d), ltrim, rtrim, b );
          int d1 = D.ITo(v,0), d2 = D.IFrom(w,0);
          Bool problem = False;
          if ( dlens[d1] - 1 < ltrim ) problem = True;
          if ( dlens[d2] - 1 < rtrim ) problem = True;
          if ( d1 == d2 && dlens[d1] - 1 < ltrim + rtrim ) problem = True;
          if (problem)
          {    nprobs++;
               continue;    }

          // Test for illegal problems.

          for ( int j = 0; j < D.From(v).isize( ); j++ )
          {    int d = D.IFrom( v, j );
               int w = to_right[d];
               ForceAssert( D.To(v).solo( ) );
               if ( !D.From(w).solo( ) ) cout << "Problem at edge " << d << endl;
               ForceAssert( D.From(w).solo( ) );    }

          // Do the edit.

          dlens[d1] -= ltrim;
          dlens[d2] -= rtrim;
          int rd1 = dinv[d1], rd2 = dinv[d2];
          if ( rd1 != d1 ) dlens[rd1] -= ltrim;
          if ( rd2 != d2 ) dlens[rd2] -= rtrim;
          for ( int j = 0; j < D.From(v).isize( ); j++ )
          {    int d = D.IFrom( v, j );
               int w = to_right[d];
               GapToSeq( D.O(d), ltrim, rtrim, b );

               // First edit the gap edge itself.

               D.OMutable(d) = { (int) tigs.size( ) };
               int E = tigs.size( );
               tigs.push_back(b);
               int rd = dinv[d];
               if ( rd != d )
               {    b.ReverseComplement( );
                    D.OMutable(rd) = { (int) tigs.size( ) };
                    tigs.push_back(b);
                    inv.push_back( E+1, E );    }
               else inv.push_back(E);

               // Now edit the adjacent edges.

               if ( j == 0 )
               {    // int d1 = D.ITo(v,0), d2 = D.IFrom(w,0);
                    ForceAssert( d1 != rd1 );
                    ForceAssert( d2 != rd2 );
                    if ( ltrim > 0 )
                    {    while(1)
                         {    int n = tigs[ D.O(d1).back( ) ].isize( ) - HBK + 1;
                              if ( n <= ltrim )
                              {    ltrim -= n;
                                   D.OMutable(d1).pop_back( );
                                   D.OMutable(rd1).pop_front( );    }
                              else break;    }
                         if ( ltrim > 0 )
                         {    b = tigs[ D.O(d1).back( ) ];
                              b.resize( b.isize( ) - ltrim );
                              D.OMutable(d1).back( ) = tigs.size( );
                              int E = tigs.size( );
                              tigs.push_back(b);    
                              b.ReverseComplement( );
                              D.OMutable(rd1).front( ) = tigs.size( );
                              tigs.push_back(b);    
                              inv.push_back( E+1, E );    }    }
                    if ( rtrim > 0 )
                    {    while(1)
                         {    int n = tigs[ D.O(d2).front( ) ].isize( ) - HBK + 1;
                              if ( n <= rtrim )
                              {    rtrim -= n;
                                   D.OMutable(d2).pop_front( );
                                   D.OMutable(rd2).pop_back( );    }
                              else break;    }
                         if ( rtrim > 0 )
                         {    b = tigs[ D.O(d2).front( ) ];
                              b.SetToSubOf( b, rtrim, b.isize( ) - rtrim );
                              D.OMutable(d2).front( ) = tigs.size( );
                              int E = tigs.size( );
                              tigs.push_back(b);    
                              b.ReverseComplement( );
                              D.OMutable(rd2).back( ) = tigs.size( );
                              tigs.push_back(b);    
                              inv.push_back( E+1, E );    }    }    }    }    }

     // Done.

     cout << Date( ) << ": there were " << nprobs << " unconverted sequence gaps"
          << endl;    }

/*
ISSUES
- Unconverted sequence gaps need to be converted to plain gaps.
- Some of these could probably be recovered.
*/
