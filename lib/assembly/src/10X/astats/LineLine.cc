// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"
#include "paths/long/large/Lines.h"
#include "10X/Heuristics.h"
#include "10X/astats/LineLine.h"

void FindLineLines(
     // inputs:
     const digraphE<vec<int>>& D, const vec<int>& dinv,
     const vec<vec<vec<vec<int>>>>& dlines,
     // outputs:
     vec<vec<vec<vec<int>>>>& dlines2,
     vec<int>* linv2p /* = nullptr */ )
{
     // Form lines into a graph.

     vec<int> to_left, to_right;
     D.ToLeft(to_left), D.ToRight(to_right);
     vec<int> verts( 2 * dlines.size( ) );
     #pragma omp parallel for
     for ( int i = 0; i < dlines.isize( ); i++ )
     {    int v = to_left[ dlines[i].front( )[0][0] ];
          int w = to_right[ dlines[i].back( )[0][0] ];
          verts[2*i] = v, verts[2*i+1] = w;    }
     UniqueSort(verts);
     int N = verts.size( );
     vec<vec<int>> from(N), to(N), from_edge_obj(N), to_edge_obj(N);
     vec<int> edges( dlines.size( ), vec<int>::IDENTITY );
     for ( int i = 0; i < dlines.isize( ); i++ )
     {    int v = BinPosition( verts, to_left[ dlines[i].front( )[0][0] ] );
          int w = BinPosition( verts, to_right[ dlines[i].back( )[0][0] ] );
          from[v].push_back(w), to[w].push_back(v);
          from_edge_obj[v].push_back(i), to_edge_obj[w].push_back(i);    }
     #pragma omp parallel for
     for ( int v = 0; v < N; v++ )
     {    SortSync( from[v], from_edge_obj[v] );
          SortSync( to[v], to_edge_obj[v] );    }
     digraphE<int> D2( from, to, edges, to_edge_obj, from_edge_obj );
     vec<int> linv;
     LineInv( dlines, dinv, linv );

     // Find lines of lines.

     cout << Date( ) << ": forming lines of lines" << endl;
     FindLines( D2, linv, dlines2, MAX_CELL_PATHS, MAX_CELL_DEPTH, False );
     cout << Date( ) << ": found " << dlines.size( ) << " lines, and "
          << dlines2.size( ) << " lines of lines" << endl;    

     if ( linv2p )  { LineInv( dlines2, linv, *linv2p ); }
}

void GetLineLineLengths( const vec<int>& llens,
     const vec<vec<vec<vec<int>>>>& dlines2, vec<int>& lens2 )
{    lens2.resize( dlines2.size( ) );
     #pragma omp parallel for
     for ( int i = 0; i < dlines2.isize( ); i++ )
     {    const vec<vec<vec<int>>>& L = dlines2[i];
          int pos = 0;
          for ( int j = 0; j < L.isize( ); j++ )
          {    vec<int> lensj;
               for ( int k = 0; k < L[j].isize( ); k++ )
               {    int len = 0;
                    for ( int l = 0; l < L[j][k].isize( ); l++ )
                    {    int u = L[j][k][l];
                         len += llens[u];    }
                    lensj.push_back(len);    }
               Sort(lensj);
               if ( lensj.nonempty( ) ) pos += Median(lensj);    }
          lens2[i] = pos;    }    }
