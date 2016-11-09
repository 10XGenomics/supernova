// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "10X/Decycle.h"
#include "10X/Gap.h"
#include "10X/Super.h"

// Remove some simple cycles.

void Decycle( const HyperBasevectorX& hb, const vec<int>& inv,
     digraphE<vec<int>>& D, vec<int>& dinv, const ReadPathVecX& paths )
{
     cout << Date( ) << ": set up for decycling" << endl;
     double clock = WallClockTime( );
     vec<int> to_left, to_right;
     D.ToLeft(to_left), D.ToRight(to_right);

     // First build a mini-index.  Note duplication with code below.

     const int MAX_PATHS = 1000;
     cout << Date( ) << ": building need_count" << endl;
     vec<Bool> need_count( hb.E( ), False );
     #pragma omp parallel for schedule(dynamic,1000)
     for ( int d = 0; d < D.E( ); d++ )
     {    if ( !IsCell( D.O(d) ) ) continue;
          cell c;
          c.CellDecode( D.O(d) );
          if ( c.G( ).E( ) != 2 ) continue;
          int v = to_left[d], w = to_right[d];
          if ( !D.To(v).solo( ) || !D.From(w).solo( ) ) continue;
          if ( c.Left( ) != 0 || c.Right( ) != 1 ) continue;
          const vec<int>& X = D.O( D.ITo(v,0) ), &Y = D.O( D.IFrom(w,0) );
          if ( X[0] < 0 || Y[0] < 0 ) continue;
          int x = X.back( ), y = Y.front( );
          int rx = inv[x], ry = inv[y];
          if ( make_pair( x,y ) >= make_pair( ry, rx ) ) continue;
          vec<int> a, b;
          for ( int v = 0; v < c.G( ).N( ); v++ )
          for ( int j = 0; j < c.G( ).From(v).isize( ); j++ )
          {    int w = c.G( ).From(v)[0];
               const vec<int>& m = c.G( ).O( c.G( ).IFrom(v,j) );
               if ( m[0] >= 0 )
               {    if ( v == 0 && w == 1 ) a = m;
                    if ( v == 1 && w == 0 ) b = m;    }    }
          if ( a.empty( ) || b.empty( ) ) continue;
          vec<int> ra(a), rb(b);
          ra.ReverseMe( ), rb.ReverseMe( );
          for ( auto& e : ra ) e = inv[e];
          for ( auto& e : rb ) e = inv[e];
          for ( auto e : { x, ry, b.front( ), b.back( ), rb.front( ), rb.back( ) } )
               need_count[e] = True;    }
     cout << Date( ) << ": building mini-index" << endl;
     vec<vec<int64_t>> paths_index( hb.E( ) );
     const int64_t batch = 10000;
     int64_t N = paths.size( );
     #pragma omp parallel for schedule(dynamic,1)
     for ( int64_t bi = 0; bi < N; bi += batch )
     {    vec<pair<int,int64_t>> stuff;
          ReadPath p;
          for ( int64_t id = bi; id < Min( bi + batch, N ); id++ )
          {    paths.unzip(p,hb,id);
               for ( int j = 0; j < (int) p.size( ); j++ )
               {    int e = p[j];
                    if ( need_count[e] ) stuff.push( e, id );    }    }
          #pragma omp critical
          {    for ( int64_t i = 0; i < stuff.jsize( ); i++ )
               {    int e = stuff[i].first;
                    int64_t id = stuff[i].second;
                    if ( paths_index[e].isize( ) <= MAX_PATHS )
                         paths_index[e].push_back(id);    }    }    }
     cout << Date( ) << ": sorting" << endl;
     #pragma omp parallel for schedule(dynamic,10000)
     for ( int e = 0; e < hb.E( ); e++ )
          Sort( paths_index[e] );

     // Main loop.

     cout << Date( ) << ": start loop" << endl;
     vec<int> dels;
     vec<vec<vec<int>>> patches;
     #pragma omp parallel for schedule(dynamic,1000)
     for ( int d = 0; d < D.E( ); d++ )
     {    if ( !IsCell( D.O(d) ) ) continue;
          cell c;
          c.CellDecode( D.O(d) );
          if ( c.G( ).E( ) != 2 ) continue;
          int v = to_left[d], w = to_right[d];
          if ( !D.To(v).solo( ) || !D.From(w).solo( ) ) continue;
          if ( c.Left( ) != 0 || c.Right( ) != 1 ) continue;
          const vec<int>& X = D.O( D.ITo(v,0) ), &Y = D.O( D.IFrom(w,0) );
          if ( X[0] < 0 || Y[0] < 0 ) continue;
          int x = X.back( ), y = Y.front( );
          int rx = inv[x], ry = inv[y];
          if ( make_pair( x,y ) >= make_pair( ry, rx ) ) continue;
          vec<int> a, b;
          for ( int v = 0; v < c.G( ).N( ); v++ )
          for ( int j = 0; j < c.G( ).From(v).isize( ); j++ )
          {    int w = c.G( ).From(v)[0];
               const vec<int>& m = c.G( ).O( c.G( ).IFrom(v,j) );
               if ( m[0] >= 0 )
               {    if ( v == 0 && w == 1 ) a = m;
                    if ( v == 1 && w == 0 ) b = m;    }    }
          if ( a.empty( ) || b.empty( ) ) continue;
          vec<int> fulls;
          vec<int> ra(a), rb(b);
          ra.ReverseMe( ), rb.ReverseMe( );
          for ( auto& e : ra ) e = inv[e];
          for ( auto& e : rb ) e = inv[e];
          Bool too_many = False;
          for ( auto e : { x, ry, b.front( ), b.back( ), rb.front( ), rb.back( ) } )
               if ( (int) paths_index[e].size( ) > MAX_PATHS ) too_many = True;
          if (too_many) continue;
          ReadPath p;
          for ( int i = 0; i < (int) paths_index[x].size( ); i++ )
          {    int64_t id = paths_index[x][i];
               paths.unzip(p,hb,id);
               for ( int j = 0; j < (int) p.size( ); j++ )
               {    if ( p[j] == x )
                    {    for ( int k = j + 1; k < (int) p.size( ); k++ )
                         {    if ( p[k] == y )
                              {    fulls.push_back( 
                                        (k-j)/2 - 1 );    }    }    }    }    }
          for ( int i = 0; i < (int) paths_index[ry].size( ); i++ )
          {    int64_t id = paths_index[ry][i];
               paths.unzip(p,hb,id);
               for ( int j = 0; j < (int) p.size( ); j++ )
               {    if ( p[j] == ry )
                    {    for ( int k = j + 1; k < (int) p.size( ); k++ )
                         {    if ( p[k] == rx )
                              {    fulls.push_back( 
                                        (k-j)/2 - 1 );    }    }    }    }    }
          Sort(fulls);
          vec<int64_t> babs, bays;
          vec<int> bab = { b.back( ) }, bay = { b.back( ) };
          bab.append(a), bay.append(a);
          bab.push_back( b.front( ) );
          bay.push_back(y);
          vec<int> q;
          for ( int i = 0; i < (int) paths_index[b.back()].size( ); i++ )
          {    int64_t id = paths_index[b.back()][i];
               paths.unzip(p,hb,id);
               q.clear( );
               for ( auto e : p ) q.push_back(e);
               if ( q.Contains(bab) ) babs.push_back(id);
               if ( q.Contains(bay) ) bays.push_back(id);    }
          bab.ReverseMe( ), bay.ReverseMe( );
          for ( auto& e : bab ) e = inv[e];
          for ( auto& e : bay ) e = inv[e];
          for ( int i = 0; i < (int) paths_index[rb.back()].size( ); i++ )
          {    int64_t id = paths_index[rb.back()][i];
               paths.unzip(p,hb,id);
               q.clear( );
               for ( auto e : p ) q.push_back(e);
               if ( q.Contains(bab) ) babs.push_back(id);
               if ( q.Contains(bay) ) bays.push_back(id);    }
          UniqueSort(babs), UniqueSort(bays);
          vec<int> ufulls(fulls);
          UniqueSort(ufulls);

          /*
          vec<int> tests = { 52718452, 50354723, 50741737, 15097537 };
          if ( Member( tests, a[0] ) )
          {
               #pragma omp critical
               {    int p = Position( tests, a.front() ) + 1;
                    PRINT2( p, a );
                    cout << "kmers = " << hb.Kmers(a.front()) << "," 
                         << hb.Kmers(b.front())
                         << "; babs = " << babs.size( ) 
                         << "; bays = " << bays.size( ) 
                         << "; fulls = " << printSeq(fulls) << endl;    }    }
          */

          if ( ( bays.isize( ) >= 3 && babs.isize( ) == 0 )
               || ( bays.isize( ) >= 10 && babs.isize( ) == 1 ) )
          {    ufulls.push_back(1);
               UniqueSort(ufulls);

               #pragma omp critical
               {    dels.push_back(d);
                    /*
                    {    int p = Position( tests, a ) + 1;
                         if ( p >= 1 ) PRINT(p);    }
                    */
                    vec<vec<int>> p;
                    for ( int i = 0; i < ufulls.isize( ); i++ )
                    {    int n = ufulls[i];
                         vec<int> x = a;
                         for ( int j = 0; j < n; j++ )
                         {    x.append(b);
                              x.append(a);    }
                         p.push_back(x);    }
                    patches.push_back(p);    }    }    }

     // Make edits;

     cout << Date( ) << ": making " << dels.size( ) << " edits" << endl;
     SortSync( dels, patches );
     int nd = dels.size( );
     for ( int i = 0; i < nd; i++ )
     {    int d = dels[i];
          vec<vec<int>>& p = patches[i];
          int v = to_left[d], w = to_right[d];
          int E = D.E( );
          for ( int j = 0; j < p.isize( ); j++ ) 
              D.AddEdgeWithUpdate( v, w, p[j], to_left, to_right );
          int rd = dinv[d];
          int rv = to_left[rd], rw = to_right[rd];
          for ( int j = 0; j < p.isize( ); j++ )
          {    p[j].ReverseMe( );
               for ( auto& e : p[j] ) e = inv[e];
               D.AddEdgeWithUpdate( rv, rw, p[j], to_left, to_right );    }
          for ( int j = 0; j < p.isize( ); j++ ) dinv.push_back( E + p.isize() + j );
          for ( int j = 0; j < p.isize( ); j++ ) dinv.push_back( E + j );
          D.DeleteEdgesWithUpdate( {d,dinv[d]}, to_left, to_right );
          dels.push_back(rd);    }

     // Clean up.

     cout << Date( ) << ": replaced " << dels.size( ) << " cells" << endl;
     RemoveUnneededVertices( D, dinv );
     CleanupCore( D, dinv );
     cout << Date( ) << ": done, time used = " << TimeSince(clock) << endl;    }
