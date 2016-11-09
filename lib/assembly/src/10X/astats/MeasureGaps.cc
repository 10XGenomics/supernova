// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>

#include "CoreTools.h"
#include "graph/Digraph.h"
#include "math/Functions.h"
#include "10X/astats/MeasureGaps.h"

void MeasureGaps( const int K, const vec<int>& kmers, const vec<int>& inv,
     const digraphE<vec<int>>& D, const vec<vec<vec<vec<int>>>>& dlines,
     const vec<int>& llens,
     const MasterVec< SerfVec<triple<int,int,int> > >& alignsb,
     vec< pair<int,int> >& gaps, const Bool verbose )
{
     vec<int> mult( kmers.size( ), 0 ), to_right;
     for ( int d = 0; d < D.E( ); d++ )
     {    if ( D.O(d)[0] < 0 ) continue;
          for ( int i = 0; i < D.O(d).isize( ); i++ )
               mult[ D.O(d)[i] ]++;    }
     D.ToRight(to_right);
     int nthreads = ( verbose ? 1 : omp_get_max_threads( ) );
     #pragma omp parallel for num_threads(nthreads)
     for ( int i = 0; i < dlines.isize( ); i++ )
     {    int pos = 0;
          const vec<vec<vec<int>>>& L = dlines[i];
          vec< vec< triple<int,int,Bool> > > places(L.size()); // (gpos-lpos,g,fw)
          vec<int> starts( L.size( ) );
          for ( int j = 0; j < L.isize( ); j++ )
          {    starts[j] = pos;
               vec<int> lensj;
               for ( int k = 0; k < L[j].isize( ); k++ )
               {    int len = 0;
                    for ( int l = 0; l < L[j][k].isize( ); l++ )
                    {    int d = L[j][k][l];
                         if ( D.O(d)[0] < 0 ) continue;
                         for ( int m = 0; m < D.O(d).isize( ); m++ )
                         {    int e = D.O(d)[m];
                              int re = inv[e];
                              if ( mult[e] == 1 && kmers[e] >= 50
                                   && alignsb[e].size( ) + alignsb[re].size( ) == 1 )
                              {    int g = -1, gpos = -1;
                                   Bool fw = True;
                                   if ( alignsb[e].size( ) == 1 )
                                   {    g = alignsb[e][0].first;
                                        gpos = alignsb[e][0].second;
                                        places[j].push( gpos-(pos+len), g, fw );    }
                                   else
                                   {    g = alignsb[re][0].first;
                                        fw = False;
                                        gpos = alignsb[re][0].second;
                                        places[j].push( 
                                             gpos+(pos+len)+kmers[e], g, fw );    }
                                             }
                              len += kmers[e];    }    }
                    lensj.push_back(len);    }
               Sort(lensj);
               if ( lensj.nonempty( ) ) pos += Median(lensj);    }
          for ( int j = 0; j < L.isize( ); j++ )
          {    if ( !L[j].solo( ) || !L[j][0].empty( ) ) continue;

               // Get gap edge id dgap: convoluted.

               if ( j == 0 || L[j-1].empty( ) || L[j-1][0].empty( ) ) 
               {    if (verbose) cout << "L" << i << "." << j << ": UNKNOWN" << endl;
                    continue;    }
               int d = L[j-1][0][0];
               int v = to_right[d];
               if ( !D.From(v).solo( ) ) 
               {    if (verbose) cout << "L" << i << "." << j << ": UNKNOWN" << endl;
                    continue;    }
               int dgap = D.IFrom(v,0);

               // Compute gap.

               Bool left_found = False, right_found = False;
               int left_offset = -1, right_offset = -1, left_g = -1, right_g = -1;
               Bool left_fw = True, right_fw = True;
               for ( int k = j-1; k >= 0; k-- )
               {    for ( int l = places[k].isize( ) - 1; l >= 0; l-- )
                    {    left_found = True;
                         left_offset = places[k][l].first;
                         left_g = places[k][l].second;
                         left_fw = places[k][l].third;
                         break;    }
                    if (left_found) break;    }
               for ( int k = j+1; k < L.isize( ); k++ )
               {    for ( int l = 0; l < places[k].isize( ); l++ )
                    {    right_found = True;
                         right_offset = places[k][l].first;
                         right_g = places[k][l].second;
                         right_fw = places[k][l].third;
                         break;    }
                    if (right_found) break;    }
               const int MAX_BELIEVABLE_GAP = 20000;
               if ( left_found && right_found && left_g == right_g
                    && left_fw == right_fw )
               {    int gap;
                    if (left_fw) gap = right_offset - left_offset - (K-1);
                    else         gap = left_offset - right_offset - (K-1);
                    if ( Abs(gap) <= MAX_BELIEVABLE_GAP )
                    {
                         #pragma omp critical
                         {    gaps.push( gap, dgap );    }
                         if (verbose)
                         {    cout << "L" << i << "." << j << " = " << dgap
                                   << ": gap " << gap << endl;    }    }
                    else if (verbose)
                    {    cout << "L" << i << "." << j << " = " << dgap
                              << " UNKNOWN" << endl;    }    }
               else if (verbose)
               {    cout << "L" << i << "." << j << " = " << dgap
                              << " UNKNOWN" << endl;    }    }    }
     Sort(gaps);    }
