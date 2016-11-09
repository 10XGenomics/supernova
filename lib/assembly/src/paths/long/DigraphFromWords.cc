///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"
#include "Equiv.h"
#include "graph/Digraph.h"
#include "math/Functions.h"
#include "paths/long/DigraphFromWords.h"
#include "paths/long/Fix64_6.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/SupportedHyperBasevector.h"

void DigraphFromWords(const vec<vec<int> >& w, digraphE<int>& D) {
     ForceAssert( w.UniqueOrdered( ));

     // Determine if first entries in w are equal, and likewise for last entries.

     Bool equal1 = True, equal2 = True;
     if (w.nonempty()) {
	  for (int v = 0; v < w.isize(); v++) {
	       if (w[v].size() < 3 || w[v].front() != w[0].front())
		    equal1 = False;
	       if (w[v].size() < 3 || w[v].back() != w[0].back())
		    equal2 = False;
	  }
     }

     // Initialize graph as disjoint union.

     int N = 0;
     for (int i = 0; i < w.isize(); i++) {
	  N += w[i].isize() + 1;
	  if (i > 0)
	       N -= 2;
	  if (i > 0 && equal1)
	       N--;
	  if (i > 0 && equal2)
	       N--;
     }
     vec<vec<int> > from(N), to(N), from_edge_obj(N), to_edge_obj(N);
     vec<int> edges;
     int n = 0;
     for (int i = 0; i < w.isize(); i++) {
	  int start = 0, stop = w[i].isize();
	  if (i > 0 && equal1)
	       start++;
	  if (i > 0 && equal2)
	       stop--;
	  for (int j = start; j < stop; j++) {
	       int x = n, y = n + 1;
	       if (j == 0 && i > 0) {
		    n--;
		    y--;
		    x = 0;
	       }
	       if (j == 1 && equal1 && i > 0) {
		    x = 1;
		    y--;
		    n--;
	       }
	       if (j + 1 == w[i].isize() && i > 0)
		    y = w[0].isize();
	       if (j + 1 == w[i].isize() - 1 && equal2 && i > 0)
		    y = w[0].isize() - 1;
	       from[x].push_back(y), to[y].push_back(x);
	       from_edge_obj[x].push_back(edges.size());
	       to_edge_obj[y].push_back(edges.size());
	       edges.push_back(w[i][j]);
	       n++;
	  }
	  if (i == 0)
	       n++;
     }
     for (int v = 0; v < N; v++) {
	  SortSync(from[v], from_edge_obj[v]);
	  SortSync(to[v], to_edge_obj[v]);
     }
     D.Initialize(from, to, edges, to_edge_obj, from_edge_obj);
     restart: vec<int> to_left, to_right;
     D.ToLeft(to_left), D.ToRight(to_right);
     edges = D.Edges();
     vec<int> edge_ids(edges.size(), vec<int>::IDENTITY);
     SortSync(edges, edge_ids);
     for (int i = 0; i < edges.isize(); i++) {
	  int j = edges.NextDiff(i);
	  for (int k1 = i; k1 < j; k1++)
	       for (int k2 = k1 + 1; k2 < j; k2++) {
		    int e1 = edge_ids[k1], e2 = edge_ids[k2];
		    digraphE<int> D2(D);
		    int l1 = to_left[e1], r1 = to_right[e1];
		    int l2 = to_left[e2], r2 = to_right[e2];

		    vec<int> es;
		    es.push_back(e2);
		    D2.DeleteEdges(es);

		    // D2.DeleteEdge(e2);

		    if (l2 != l1)
			 D2.TransferEdges(l2, l1);
		    if (r2 != r1)
			 D2.TransferEdges(r2, r1);
		    D2.RemoveEdgelessVertices();
		    if (!D2.Acyclic())
			 continue;
		    D2.RemoveDeadEdgeObjects();
		    vec<vec<Bool> > to_delete_from(D2.N()), to_delete_to(
			      D2.N());
		    for (int v = 0; v < D2.N(); v++) {
			 to_delete_from[v].resize(D2.From(v).size(), False);
			 to_delete_to[v].resize(D2.To(v).size(), False);
		    }
		    for (int v = 0; v < D2.N(); v++) {
			 vec<pair<int, int> > eval;
			 vec<int> id(D2.From(v).size(), vec<int>::IDENTITY);
			 for (int l = 0; l < D2.From(v).isize(); l++)
			      eval.push(D2.From(v)[l],
					D2.EdgeObjectByIndexFrom(v, l));
			 SortSync(eval, id);
			 for (int l = 0; l < eval.isize(); l++) {
			      int m = eval.NextDiff(l);
			      for (int k = l + 1; k < m; k++) {
				   to_delete_from[v][id[k]] = True;
				   int w = D2.From(v)[id[k]];
				   to_delete_to[w][D2.InputFromOutputTo(v,
					     id[k])] = True;
			      }
			      l = m - 1;
			 }
		    }
		    for (int v = 0; v < D2.N(); v++) {
			 EraseIf(D2.FromMutable(v), to_delete_from[v]);
			 EraseIf(D2.FromEdgeObjMutable(v), to_delete_from[v]);
			 EraseIf(D2.ToMutable(v), to_delete_to[v]);
			 EraseIf(D2.ToEdgeObjMutable(v), to_delete_to[v]);
		    }

		    vec<vec<int> > paths;
		    for (int v = 0; v < D2.N(); v++) {
			 if (!D2.Source(v))
			      continue;
			 for (int w = 0; w < D2.N(); w++) {
			      if (!D2.Sink(w))
				   continue;
			      vec<vec<int> > pathsvw;
                              // could upgrade to faster version that uses
                              // to_left, to_right
			      D2.EdgePaths(v, w, pathsvw);
			      paths.append(pathsvw);
			 }
		    }
		    for (int m = 0; m < paths.isize(); m++)
			 for (int l = 0; l < paths[m].isize(); l++)
			      paths[m][l] = D2.EdgeObject(paths[m][l]);
		    Sort(paths);
		    if (paths == w) {
			 D = D2;
			 goto restart;
		    }
	       }
	  i = j - 1;
     }
}

void WordsToDigraphAlt( const vec< vec<int> >& w, const vec<fix64_6>& ww,
     digraphE< vec<int> >& D, vec<int>& inv, 
     const long_logging_control& log_control, const long_logging& logc )
{    
     // Validate input, and create an index to w.

     double clock = WallClockTime( );
     int wmax = -1;
     for ( int i = 0; i < w.isize( ); i++ )
     for ( int j = 0; j < w[i].isize( ); j++ )
     {    ForceAssertGe( w[i][j], 0 );
          wmax = Max( wmax, w[i][j] );    }
     vec< vec< pair<int,int> > > windex( wmax + 1 );
     for ( int i = 0; i < w.isize( ); i++ )
     for ( int j = 0; j < w[i].isize( ); j++ )
               windex[ w[i][j] ].push( i, j );
     
     // Find overlaps.

     vec< vec< triple<int,int,int> > > over( w.size( ) );
     for ( int i1 = 0; i1 < w.isize( ); i1++ )
     {    for ( int j1 = 0; j1 < w[i1].isize( ); j1++ )
          {    for ( int l = 0; l < windex[ w[i1][j1] ].isize( ); l++ )
               {    int i2 = windex[ w[i1][j1] ][l].first;
                    int j2 = windex[ w[i1][j1] ][l].second;
                    int o = j1 - j2;
                    int ov = IntervalOverlap( 
                         0, w[i1].isize( ), o, o + w[i2].isize( ) );
                    if ( i1 == i2 && o == 0 ) continue;
                    if ( Overlap( w[i1], w[i2], o ) ) 
                         over[i1].push( i2, o, ov );    }    }
          UniqueSort( over[i1] );    }

     // XXX:
     /*
     for ( int i = 0; i < w.isize( ); i++ )
     for ( int j = 0; j < over[i].isize( ); j++ )
     {    if ( inv[ w[i][0] ] < 0 ) continue;
          int ip = over[i][j].first;
          vec<int> p( w[i] );
          p.ReverseMe( );
          for ( int l = 0; l < p.isize( ); l++ )
               p[l] = inv[ p[l] ];
          int i2 = BinPosition( w, p );
          ForceAssertGe( i2, 0 );
          vec<int> pp( w[ip] );
          pp.ReverseMe( );
          for ( int l = 0; l < pp.isize( ); l++ )
               pp[l] = inv[ pp[l] ];
          int ip2 = BinPosition( w, pp );
          ForceAssertGe( ip2, 0 );
          Bool found = False;
          for ( int j2 = 0; j2 < over[i2].isize( ); j2++ )
          {    if ( over[i2][j2].first == ip2
                    && over[i2][j2].third == over[i][j].third )
               {    found = True;    }    }
          ForceAssert(found);    }
     */

     // Filter overlaps by weight.

     const fix64_6 minw = 5;
     vec< vec<Bool> > to_del( w.size( ) );
     for ( int i = 0; i < w.isize( ); i++ )
          to_del[i].resize( over[i].size( ), False );
     for ( int i = 0; i < w.isize( ); i++ )
     {    vec<Bool> right_extend( over[i].size( ), False );
          for ( int j = 0; j < over[i].isize( ); j++ )
          {    if ( over[i][j].second + w[ over[i][j].first ].isize( )
                    > w[i].isize( ) )
               {    right_extend[j] = True;    }    }
          int M = 0;
          for ( int j = 0; j < over[i].isize( ); j++ )
               if ( right_extend[j] ) M = Max( M, over[i][j].third );
          vec<fix64_6> ow( M+1, 0 );
          for ( int j = 0; j < over[i].isize( ); j++ )
          {    if ( right_extend[j] )
               {    for ( int l = 1; l <= over[i][j].third; l++ )
                         ow[l] += ww[ over[i][j].first ];    }    }
          int x;
          for ( x = M; x >= 2; x-- )
               if ( ow[x] >= minw ) break;
          for ( int j = 0; j < over[i].isize( ); j++ )
          {    if ( right_extend[j] && over[i][j].third < x ) 
               {    to_del[i][j] = True;
                    int ip = over[i][j].first;
                    for ( int l = 0; l < over[ip].isize( ); l++ )
                    {    if ( over[ip][l].first == i 
                              && over[ip][l].third == over[i][j].third )
                         {    to_del[ip][l] = True;    }    }    }    }    }

     for ( int i = 0; i < w.isize( ); i++ )
     {    vec<Bool> left_extend( over[i].size( ), False );
          for ( int j = 0; j < over[i].isize( ); j++ )
               if ( over[i][j].second < 0 ) left_extend[j] = True;
          int M = 0;
          for ( int j = 0; j < over[i].isize( ); j++ )
               if ( left_extend[j] ) M = Max( M, over[i][j].third );
          vec<fix64_6> ow( M+1, 0 );
          for ( int j = 0; j < over[i].isize( ); j++ )
          {    if ( left_extend[j] )
               {    for ( int l = 1; l <= over[i][j].third; l++ )
                         ow[l] += ww[ over[i][j].first ];    }    }
          int x;
          for ( x = M; x >= 2; x-- )
               if ( ow[x] >= minw ) break;
          for ( int j = 0; j < over[i].isize( ); j++ )
          {    if ( left_extend[j] && over[i][j].third < x ) 
               {    to_del[i][j] = True;
                    int ip = over[i][j].first;
                    for ( int l = 0; l < over[ip].isize( ); l++ )
                    {    if ( over[ip][l].first == i 
                              && over[ip][l].third == over[i][j].third )
                         {    to_del[ip][l] = True;    }    }    }    }    }

     /*
     cout << "\noverlaps:\n"; // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     for ( int i = 0; i < w.isize( ); i++ ) // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     for ( int j = 0; j < over[i].isize( ); j++ ) // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     {    int i1 = i, i2 = over[i][j].first; // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          int ov = over[i][j].third; // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          cout << printSeq( w[i1] ) << " overlaps " << printSeq( w[i2] ) // XXXXXXXX
               << " by " << ov; // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          if ( to_del[i][j] ) cout << " [delete]"; // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          cout << "\n";    } // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     */

     // Forcing symmetry.

     for ( int i = 0; i < w.isize( ); i++ )
     for ( int j = 0; j < over[i].isize( ); j++ )
     {    if ( !to_del[i][j] ) continue;
          if ( inv[ w[i][0] ] < 0 ) continue;
          int ip = over[i][j].first;
          vec<int> p( w[i] );
          p.ReverseMe( );
          for ( int l = 0; l < p.isize( ); l++ )
               p[l] = inv[ p[l] ];
          int i2 = BinPosition( w, p );
          ForceAssertGe( i2, 0 );
          vec<int> pp( w[ip] );
          pp.ReverseMe( );
          for ( int l = 0; l < pp.isize( ); l++ )
               pp[l] = inv[ pp[l] ];
          int ip2 = BinPosition( w, pp );
          ForceAssertGe( ip2, 0 );
          for ( int j2 = 0; j2 < over[i2].isize( ); j2++ )
          {    if ( over[i2][j2].first == ip2
                    && over[i2][j2].third == over[i][j].third )
               {    to_del[i2][j2] = True;    }    }    }
          
     for ( int i = 0; i < w.isize( ); i++ )
          EraseIf( over[i], to_del[i] );

     // XXX:
     /*
     for ( int i = 0; i < w.isize( ); i++ )
     for ( int j = 0; j < over[i].isize( ); j++ )
     {    if ( inv[ w[i][0] ] < 0 ) continue;
          int ip = over[i][j].first;
          vec<int> p( w[i] );
          p.ReverseMe( );
          for ( int l = 0; l < p.isize( ); l++ )
               p[l] = inv[ p[l] ];
          int i2 = BinPosition( w, p );
          ForceAssertGe( i2, 0 );
          vec<int> pp( w[ip] );
          pp.ReverseMe( );
          for ( int l = 0; l < pp.isize( ); l++ )
               pp[l] = inv[ pp[l] ];
          int ip2 = BinPosition( w, pp );
          ForceAssertGe( ip2, 0 );
          Bool found = False;
          for ( int j2 = 0; j2 < over[i2].isize( ); j2++ )
          {    if ( over[i2][j2].first == ip2
                    && over[i2][j2].third == over[i][j].third )
               {    found = True;    }    }
          ForceAssert(found);    }
     */

     // From the digraph W that is the disjoint union of the words, with each word
     // stretched out over a sequence of edges (corresponding to the letters in
     // the word).

     int nv = 0, ne = 0;
     for ( int i = 0; i < w.isize( ); i++ )
     {    ne += w[i].size( );
          nv += w[i].size( ) + 1;    }
     vec< vec<int> > from(nv), to(nv), from_edge_obj(nv), to_edge_obj(nv);
     vec<int> edges(ne);
     int iv = 0, ie = 0;
     vec<int> vstart( w.size( ) ), estart( w.size( ) );
     vec< pair<int,int> > vorigin, eorigin;
     for ( int i = 0; i < w.isize( ); i++ )
     {    vstart[i] = iv, estart[i] = ie;
          for ( int j = 0; j < w[i].isize( ); j++ )
          {    vorigin.push( i, j );
               eorigin.push( i, j );
               from[iv].push_back(iv+1), to[iv+1].push_back(iv);
               from_edge_obj[iv].push_back(ie);
               to_edge_obj[iv+1].push_back(ie);
               edges[ie++] = w[i][j];
               iv++;    }
          vorigin.push( i, w[i].size( ) );
          iv++;    }
     digraphE<int> W( from, to, edges, to_edge_obj, from_edge_obj );

     // From the overlaps, deduce equivalence relations on the vertices and edges
     // of W.

     equiv_rel evert(iv), eedge(ie);
     for ( int i1 = 0; i1 < over.isize( ); i1++ )
     for ( int j = 0; j < over[i1].isize( ); j++ )
     {    int i2 = over[i1][j].first, o = over[i1][j].second;
          for ( int l1 = 0; l1 <= w[i1].isize( ); l1++ )
          {    int l2 = l1 - o;
               if ( l2 < 0 || l2 > w[i2].isize( ) ) continue;
               if ( l1 == w[i1].isize( ) || l2 == w[i2].isize( ) )
                    evert.Join( vstart[i1] + l1, vstart[i2] + l2 );
               if ( l1 == w[i1].isize( ) || l2 == w[i2].isize( ) ) continue;
               if ( l1 == w[i1].isize( ) ) continue;
               evert.Join( vstart[i1] + l1, vstart[i2] + l2 );
               eedge.Join( estart[i1] + l1, estart[i2] + l2 );    }    }

     // In cases where the resulting graph would have sources that emanate long
     // the same edge, join them.  Ditto for sinks.

     {    vec<int> vreps, ereps;
          evert.OrbitRepsAlt(vreps), eedge.OrbitRepsAlt(ereps);
          int nvreps = vreps.size( ), nereps = ereps.size( );
          vec< vec<int> > from2(nvreps), to2(nvreps); 
          vec< vec<int> > from_edge_obj2(nvreps), to_edge_obj2(nvreps);
          vec<int> edges2(nereps);
          for ( int v1 = 0; v1 < W.N( ); v1++ )
          for ( int j = 0; j < W.From(v1).isize( ); j++ )
          {    int v2 = W.From(v1)[j];
               int e = W.EdgeObjectIndexByIndexFrom( v1, j );
               int pv1 = BinPosition( vreps, evert.ClassId(v1) );
               int pv2 = BinPosition( vreps, evert.ClassId(v2) );
               int pe = BinPosition( ereps, eedge.ClassId(e) );
               if ( Member( from_edge_obj2[pv1], pe ) ) continue;
               from2[pv1].push_back(pv2), to2[pv2].push_back(pv1);
               from_edge_obj2[pv1].push_back(pe); 
               to_edge_obj2[pv2].push_back(pe);    }
          for ( int i = 0; i < nereps; i++ )
          {    int x = w[ eorigin[ ereps[i] ].first ][ eorigin[ ereps[i] ].second ];
               edges2[i] = x;    }
          for ( int i = 0; i < nvreps; i++ )
          {    SortSync( from2[i], from_edge_obj2[i] );
               SortSync( to2[i], to_edge_obj2[i] );    }
          digraphE<int> D( from2, to2, edges2, to_edge_obj2, from_edge_obj2 );
          vec<int> sources, sinks;
          D.Sources(sources), D.Sinks(sinks);
          vec< pair<int,int> > source_edges, sink_edges;
          for ( int i = 0; i < sources.isize( ); i++ )
          {    int v = sources[i];
               for ( int j = 0; j < D.From(v).isize( ); j++ )
                    source_edges.push( D.EdgeObjectByIndexFrom( v, j ), v );    }
          for ( int i = 0; i < sinks.isize( ); i++ )
          {    int v = sinks[i];
               for ( int j = 0; j < D.To(v).isize( ); j++ )
                    sink_edges.push( D.EdgeObjectByIndexTo( v, j ), v );    }
          UniqueSort(source_edges), UniqueSort(sink_edges);
          for ( int i = 0; i < source_edges.isize( ); i++ )
          {    int j;
               for ( j = i + 1; j < source_edges.isize( ); j++ )
                    if ( source_edges[j].first != source_edges[i].first ) break;
               for ( int k = i + 1; k < j; k++ )
               {    int v1 = source_edges[i].second, v2 = source_edges[k].second;
                    evert.Join( vreps[v1], vreps[v2] );    }
               i = j - 1;    }
          for ( int i = 0; i < sink_edges.isize( ); i++ )
          {    int j;
               for ( j = i + 1; j < sink_edges.isize( ); j++ )
                    if ( sink_edges[j].first != sink_edges[i].first ) break;
               for ( int k = i + 1; k < j; k++ )
               {    int v1 = sink_edges[i].second, v2 = sink_edges[k].second;
                    evert.Join( vreps[v1], vreps[v2] );    }
               i = j - 1;    }    }

     // Zipper.

     while(1)
     {    int joins = 0;
          vec<int> vreps, ereps;
          evert.OrbitRepsAlt(vreps), eedge.OrbitRepsAlt(ereps);
          int nvreps = vreps.size( ), nereps = ereps.size( );
          vec< vec<int> > from2(nvreps), to2(nvreps); 
          vec< vec<int> > from_edge_obj2(nvreps), to_edge_obj2(nvreps);
          vec< vec<int> > edges2(nereps);
          for ( int v1 = 0; v1 < W.N( ); v1++ )
          for ( int j = 0; j < W.From(v1).isize( ); j++ )
          {    int v2 = W.From(v1)[j];
               int e = W.EdgeObjectIndexByIndexFrom( v1, j );
               int pv1 = BinPosition( vreps, evert.ClassId(v1) );
               int pv2 = BinPosition( vreps, evert.ClassId(v2) );
               int pe = BinPosition( ereps, eedge.ClassId(e) );
               if ( Member( from_edge_obj2[pv1], pe ) ) continue;
               from2[pv1].push_back(pv2), to2[pv2].push_back(pv1);
               from_edge_obj2[pv1].push_back(pe); 
               to_edge_obj2[pv2].push_back(pe);    }
          for ( int i = 0; i < nereps; i++ )
          {    vec<int> x(1);
               x[0] = w[ eorigin[ ereps[i] ].first ][ eorigin[ ereps[i] ].second ];
               edges2[i] = x;    }
          #pragma omp parallel for
          for ( int v = 0; v < nvreps; v++ )
          {    vec< triple<int,int,int> > s;
               for ( int j = 0; j < from2[v].isize( ); j++ )
               {    int ei = from_edge_obj2[v][j];
                    s.push( edges2[ei][0], from2[v][j], j );    }
               Sort(s);
               for ( int j = 0; j < s.isize( ); j++ )
               {    int k;
                    for ( k = j + 1; k < s.isize( ); k++ )
                         if ( s[k].first != s[j].first ) break;
                    for ( int l = j + 1; l < k; l++ )
                    {    
                         #pragma omp critical
                         {    evert.Join( vreps[ s[l].second ], 
                                   vreps[ s[j].second ] );
                              eedge.Join( ereps[ from_edge_obj2[v][ s[l].third ] ], 
                                   ereps[ from_edge_obj2[v][ s[j].third ] ] );
                              joins++;    }    }
                    j = k - 1;    }    }
          #pragma omp parallel for
          for ( int v = 0; v < nvreps; v++ )
          {    vec< triple<int,int,int> > s;
               for ( int j = 0; j < to2[v].isize( ); j++ )
               {    int ei = to_edge_obj2[v][j];
                    s.push( edges2[ei][0], to2[v][j], j );    }
               Sort(s);
               for ( int j = 0; j < s.isize( ); j++ )
               {    int k;
                    for ( k = j + 1; k < s.isize( ); k++ )
                         if ( s[k].first != s[j].first ) break;
                    for ( int l = j + 1; l < k; l++ )
                    {    
                         #pragma omp critical
                         {    evert.Join( vreps[ s[l].second ], 
                                   vreps[ s[j].second ] );
                              eedge.Join( ereps[ to_edge_obj2[v][ s[l].third ] ], 
                                   ereps[ to_edge_obj2[v][ s[j].third ] ] );
                              joins++;    }    }
                    j = k - 1;    }    }
          if ( joins == 0 ) break;    }

     // Apply the equivalence relation to W, yielding the digraph D.

     vec<int> vreps, ereps;
     evert.OrbitRepsAlt(vreps), eedge.OrbitRepsAlt(ereps);
     int nvreps = vreps.size( ), nereps = ereps.size( );
     vec< vec<int> > from2(nvreps), to2(nvreps); 
     vec< vec<int> > from_edge_obj2(nvreps), to_edge_obj2(nvreps);
     vec< vec<int> > edges2(nereps);
     for ( int v1 = 0; v1 < W.N( ); v1++ )
     for ( int j = 0; j < W.From(v1).isize( ); j++ )
     {    int v2 = W.From(v1)[j];
          int e = W.EdgeObjectIndexByIndexFrom( v1, j );
          int pv1 = BinPosition( vreps, evert.ClassId(v1) );
          int pv2 = BinPosition( vreps, evert.ClassId(v2) );
          int pe = BinPosition( ereps, eedge.ClassId(e) );
          if ( Member( from_edge_obj2[pv1], pe ) ) continue;
          from2[pv1].push_back(pv2), to2[pv2].push_back(pv1);
          from_edge_obj2[pv1].push_back(pe), to_edge_obj2[pv2].push_back(pe);    }
     for ( int i = 0; i < nereps; i++ )
     {    vec<int> x(1);
          x[0] = w[ eorigin[ ereps[i] ].first ][ eorigin[ ereps[i] ].second ];
          edges2[i] = x;    }
     for ( int i = 0; i < nvreps; i++ )
     {    SortSync( from2[i], from_edge_obj2[i] );
          SortSync( to2[i], to_edge_obj2[i] );    }
     D.Initialize( from2, to2, edges2, to_edge_obj2, from_edge_obj2 );

     // Compute involution.

     vec<int> inv2( D.EdgeObjectCount( ), -1 );
     for ( int e = 0; e < D.EdgeObjectCount( ); e++ )
     {    int i = eorigin[ ereps[e] ].first, pos = eorigin[ ereps[e] ].second;
          const vec<int>& y = w[i];
          vec<int> yinv;
          if ( inv[ y[0] ] < 0 ) continue;
          for ( int j = 0; j < y.isize( ); j++ )
               yinv.push_back( inv[ y[j] ] );
          yinv.ReverseMe( );
          int posinv = y.isize( ) - pos - 1;
          int iinv = Position( w, yinv );
          ForceAssertGe( iinv, 0 );
          int einv = BinPosition( ereps, eedge.ClassId( estart[iinv] + posinv ) );
          ForceAssertGe( einv, 0 );
          inv2[e] = einv;    }
     inv = inv2;
     REPORT_TIME( clock, "used in WordsToDigraph" );    }

void WordsToDigraph( const vec< vec<int> >& w, digraphE< vec<int> >& D,
     vec<int>& inv, const long_logging_control& log_control, const long_logging& logc )
{    
     // Validate input, and create an index to w.

     double clock = WallClockTime( );
     int wmax = -1;
     for ( int i = 0; i < w.isize( ); i++ )
     for ( int j = 0; j < w[i].isize( ); j++ )
     {    ForceAssertGe( w[i][j], 0 );
          wmax = Max( wmax, w[i][j] );    }
     vec< vec< pair<int,int> > > windex( wmax + 1 );
     for ( int i = 0; i < w.isize( ); i++ )
     for ( int j = 0; j < w[i].isize( ); j++ )
               windex[ w[i][j] ].push( i, j );
     
     // Find overlaps.

     vec< vec< pair<int,int> > > over( w.size( ) );
     for ( int i1 = 0; i1 < w.isize( ); i1++ )
     {    for ( int j1 = 0; j1 < w[i1].isize( ); j1++ )
          {    for ( int l = 0; l < windex[ w[i1][j1] ].isize( ); l++ )
               {    int i2 = windex[ w[i1][j1] ][l].first;
                    int j2 = windex[ w[i1][j1] ][l].second;
                    int o = j1 - j2;
                    if ( i1 == i2 && o == 0 ) continue;
                    if ( Overlap( w[i1], w[i2], o ) ) 
                         over[i1].push( i2, o );    }    }
          UniqueSort( over[i1] );    }

     // From the digraph W that is the disjoint union of the words, with each word
     // stretched out over a sequence of edges (corresponding to the letters in
     // the word).

     int nv = 0, ne = 0;
     for ( int i = 0; i < w.isize( ); i++ )
     {    ne += w[i].size( );
          nv += w[i].size( ) + 1;    }
     vec< vec<int> > from(nv), to(nv), from_edge_obj(nv), to_edge_obj(nv);
     vec<int> edges(ne);
     int iv = 0, ie = 0;
     vec<int> vstart( w.size( ) ), estart( w.size( ) );
     vec< pair<int,int> > vorigin, eorigin;
     for ( int i = 0; i < w.isize( ); i++ )
     {    vstart[i] = iv, estart[i] = ie;
          for ( int j = 0; j < w[i].isize( ); j++ )
          {    vorigin.push( i, j );
               eorigin.push( i, j );
               from[iv].push_back(iv+1), to[iv+1].push_back(iv);
               from_edge_obj[iv].push_back(ie);
               to_edge_obj[iv+1].push_back(ie);
               edges[ie++] = w[i][j];
               iv++;    }
          vorigin.push( i, w[i].size( ) );
          iv++;    }
     digraphE<int> W( from, to, edges, to_edge_obj, from_edge_obj );

     // From the overlaps, deduce equivalence relations on the vertices and edges
     // of W.

     equiv_rel evert(iv), eedge(ie);
     for ( int i1 = 0; i1 < over.isize( ); i1++ )
     for ( int j = 0; j < over[i1].isize( ); j++ )
     {    int i2 = over[i1][j].first, o = over[i1][j].second;
          for ( int l1 = 0; l1 <= w[i1].isize( ); l1++ )
          {    int l2 = l1 - o;
               if ( l2 < 0 || l2 > w[i2].isize( ) ) continue;
               if ( l1 == w[i1].isize( ) || l2 == w[i2].isize( ) )
                    evert.Join( vstart[i1] + l1, vstart[i2] + l2 );
               if ( l1 == w[i1].isize( ) || l2 == w[i2].isize( ) ) continue;
               if ( l1 == w[i1].isize( ) ) continue;
               evert.Join( vstart[i1] + l1, vstart[i2] + l2 );
               eedge.Join( estart[i1] + l1, estart[i2] + l2 );    }    }

     // In cases where the resulting graph would have sources that emanate long
     // the same edge, join them.  Ditto for sinks.

     {    vec<int> vreps, ereps;
          evert.OrbitRepsAlt(vreps), eedge.OrbitRepsAlt(ereps);
          int nvreps = vreps.size( ), nereps = ereps.size( );
          vec< vec<int> > from2(nvreps), to2(nvreps); 
          vec< vec<int> > from_edge_obj2(nvreps), to_edge_obj2(nvreps);
          vec<int> edges2(nereps);
          for ( int v1 = 0; v1 < W.N( ); v1++ )
          for ( int j = 0; j < W.From(v1).isize( ); j++ )
          {    int v2 = W.From(v1)[j];
               int e = W.EdgeObjectIndexByIndexFrom( v1, j );
               int pv1 = BinPosition( vreps, evert.ClassId(v1) );
               int pv2 = BinPosition( vreps, evert.ClassId(v2) );
               int pe = BinPosition( ereps, eedge.ClassId(e) );
               if ( Member( from_edge_obj2[pv1], pe ) ) continue;
               from2[pv1].push_back(pv2), to2[pv2].push_back(pv1);
               from_edge_obj2[pv1].push_back(pe); 
               to_edge_obj2[pv2].push_back(pe);    }
          for ( int i = 0; i < nereps; i++ )
          {    int x = w[ eorigin[ ereps[i] ].first ][ eorigin[ ereps[i] ].second ];
               edges2[i] = x;    }
          for ( int i = 0; i < nvreps; i++ )
          {    SortSync( from2[i], from_edge_obj2[i] );
               SortSync( to2[i], to_edge_obj2[i] );    }
          digraphE<int> D( from2, to2, edges2, to_edge_obj2, from_edge_obj2 );
          vec<int> sources, sinks;
          D.Sources(sources), D.Sinks(sinks);
          vec< pair<int,int> > source_edges, sink_edges;
          for ( int i = 0; i < sources.isize( ); i++ )
          {    int v = sources[i];
               for ( int j = 0; j < D.From(v).isize( ); j++ )
                    source_edges.push( D.EdgeObjectByIndexFrom( v, j ), v );    }
          for ( int i = 0; i < sinks.isize( ); i++ )
          {    int v = sinks[i];
               for ( int j = 0; j < D.To(v).isize( ); j++ )
                    sink_edges.push( D.EdgeObjectByIndexTo( v, j ), v );    }
          UniqueSort(source_edges), UniqueSort(sink_edges);
          for ( int i = 0; i < source_edges.isize( ); i++ )
          {    int j;
               for ( j = i + 1; j < source_edges.isize( ); j++ )
                    if ( source_edges[j].first != source_edges[i].first ) break;
               for ( int k = i + 1; k < j; k++ )
               {    int v1 = source_edges[i].second, v2 = source_edges[k].second;
                    evert.Join( vreps[v1], vreps[v2] );    }
               i = j - 1;    }
          for ( int i = 0; i < sink_edges.isize( ); i++ )
          {    int j;
               for ( j = i + 1; j < sink_edges.isize( ); j++ )
                    if ( sink_edges[j].first != sink_edges[i].first ) break;
               for ( int k = i + 1; k < j; k++ )
               {    int v1 = sink_edges[i].second, v2 = sink_edges[k].second;
                    evert.Join( vreps[v1], vreps[v2] );    }
               i = j - 1;    }    }

     // Zipper.

     while(1)
     {    int joins = 0;
          vec<int> vreps, ereps;
          evert.OrbitRepsAlt(vreps), eedge.OrbitRepsAlt(ereps);
          int nvreps = vreps.size( ), nereps = ereps.size( );
          vec< vec<int> > from2(nvreps), to2(nvreps); 
          vec< vec<int> > from_edge_obj2(nvreps), to_edge_obj2(nvreps);
          vec< vec<int> > edges2(nereps);
          for ( int v1 = 0; v1 < W.N( ); v1++ )
          for ( int j = 0; j < W.From(v1).isize( ); j++ )
          {    int v2 = W.From(v1)[j];
               int e = W.EdgeObjectIndexByIndexFrom( v1, j );
               int pv1 = BinPosition( vreps, evert.ClassId(v1) );
               int pv2 = BinPosition( vreps, evert.ClassId(v2) );
               int pe = BinPosition( ereps, eedge.ClassId(e) );
               if ( Member( from_edge_obj2[pv1], pe ) ) continue;
               from2[pv1].push_back(pv2), to2[pv2].push_back(pv1);
               from_edge_obj2[pv1].push_back(pe); 
               to_edge_obj2[pv2].push_back(pe);    }
          for ( int i = 0; i < nereps; i++ )
          {    vec<int> x(1);
               x[0] = w[ eorigin[ ereps[i] ].first ][ eorigin[ ereps[i] ].second ];
               edges2[i] = x;    }
          #pragma omp parallel for
          for ( int v = 0; v < nvreps; v++ )
          {    vec< triple<int,int,int> > s;
               for ( int j = 0; j < from2[v].isize( ); j++ )
               {    int ei = from_edge_obj2[v][j];
                    s.push( edges2[ei][0], from2[v][j], j );    }
               Sort(s);
               for ( int j = 0; j < s.isize( ); j++ )
               {    int k;
                    for ( k = j + 1; k < s.isize( ); k++ )
                         if ( s[k].first != s[j].first ) break;
                    for ( int l = j + 1; l < k; l++ )
                    {    
                         #pragma omp critical
                         {    evert.Join( vreps[ s[l].second ], 
                                   vreps[ s[j].second ] );
                              eedge.Join( ereps[ from_edge_obj2[v][ s[l].third ] ], 
                                   ereps[ from_edge_obj2[v][ s[j].third ] ] );
                              joins++;    }    }
                    j = k - 1;    }    }
          #pragma omp parallel for
          for ( int v = 0; v < nvreps; v++ )
          {    vec< triple<int,int,int> > s;
               for ( int j = 0; j < to2[v].isize( ); j++ )
               {    int ei = to_edge_obj2[v][j];
                    s.push( edges2[ei][0], to2[v][j], j );    }
               Sort(s);
               for ( int j = 0; j < s.isize( ); j++ )
               {    int k;
                    for ( k = j + 1; k < s.isize( ); k++ )
                         if ( s[k].first != s[j].first ) break;
                    for ( int l = j + 1; l < k; l++ )
                    {    
                         #pragma omp critical
                         {    evert.Join( vreps[ s[l].second ], 
                                   vreps[ s[j].second ] );
                              eedge.Join( ereps[ to_edge_obj2[v][ s[l].third ] ], 
                                   ereps[ to_edge_obj2[v][ s[j].third ] ] );
                              joins++;    }    }
                    j = k - 1;    }    }
          if ( joins == 0 ) break;    }

     // Apply the equivalence relation to W, yielding the digraph D.

     vec<int> vreps, ereps;
     evert.OrbitRepsAlt(vreps), eedge.OrbitRepsAlt(ereps);
     int nvreps = vreps.size( ), nereps = ereps.size( );
     vec< vec<int> > from2(nvreps), to2(nvreps); 
     vec< vec<int> > from_edge_obj2(nvreps), to_edge_obj2(nvreps);
     vec< vec<int> > edges2(nereps);
     for ( int v1 = 0; v1 < W.N( ); v1++ )
     for ( int j = 0; j < W.From(v1).isize( ); j++ )
     {    int v2 = W.From(v1)[j];
          int e = W.EdgeObjectIndexByIndexFrom( v1, j );
          int pv1 = BinPosition( vreps, evert.ClassId(v1) );
          int pv2 = BinPosition( vreps, evert.ClassId(v2) );
          int pe = BinPosition( ereps, eedge.ClassId(e) );
          if ( Member( from_edge_obj2[pv1], pe ) ) continue;
          from2[pv1].push_back(pv2), to2[pv2].push_back(pv1);
          from_edge_obj2[pv1].push_back(pe), to_edge_obj2[pv2].push_back(pe);    }
     for ( int i = 0; i < nereps; i++ )
     {    vec<int> x(1);
          x[0] = w[ eorigin[ ereps[i] ].first ][ eorigin[ ereps[i] ].second ];
          edges2[i] = x;    }
     for ( int i = 0; i < nvreps; i++ )
     {    SortSync( from2[i], from_edge_obj2[i] );
          SortSync( to2[i], to_edge_obj2[i] );    }
     D.Initialize( from2, to2, edges2, to_edge_obj2, from_edge_obj2 );

     // Compute involution.

     vec<int> inv2( D.EdgeObjectCount( ), -1 );
     for ( int e = 0; e < D.EdgeObjectCount( ); e++ )
     {    int i = eorigin[ ereps[e] ].first, pos = eorigin[ ereps[e] ].second;
          const vec<int>& y = w[i];
          vec<int> yinv;
          if ( inv[ y[0] ] < 0 ) continue;
          for ( int j = 0; j < y.isize( ); j++ )
               yinv.push_back( inv[ y[j] ] );
          yinv.ReverseMe( );
          int posinv = y.isize( ) - pos - 1;
          int iinv = Position( w, yinv );
          ForceAssertGe( iinv, 0 );
          int einv = BinPosition( ereps, eedge.ClassId( estart[iinv] + posinv ) );
          ForceAssertGe( einv, 0 );
          inv2[e] = einv;    }
     inv = inv2;
     REPORT_TIME( clock, "used in WordsToDigraph" );    }

void WordsToDigraphAlt2( const vec< vec<int> >& w, 
     const vec< vec< pair<int,int> > >& over, digraphE< vec<int> >& D,
     vec<int>& inv, const long_logging_control& log_control, const long_logging& logc )
{    
     // Validate input, and create an index to w.

     double clock = WallClockTime( );
     int wmax = -1;
     for ( int i = 0; i < w.isize( ); i++ )
     for ( int j = 0; j < w[i].isize( ); j++ )
     {    ForceAssertGe( w[i][j], 0 );
          wmax = Max( wmax, w[i][j] );    }
     vec< vec< pair<int,int> > > windex( wmax + 1 );
     for ( int i = 0; i < w.isize( ); i++ )
     for ( int j = 0; j < w[i].isize( ); j++ )
               windex[ w[i][j] ].push( i, j );
     
     // From the digraph W that is the disjoint union of the words, with each word
     // stretched out over a sequence of edges (corresponding to the letters in
     // the word).

     int nv = 0, ne = 0;
     for ( int i = 0; i < w.isize( ); i++ )
     {    ne += w[i].size( );
          nv += w[i].size( ) + 1;    }
     vec< vec<int> > from(nv), to(nv), from_edge_obj(nv), to_edge_obj(nv);
     vec<int> edges(ne);
     int iv = 0, ie = 0;
     vec<int> vstart( w.size( ) ), estart( w.size( ) );
     vec< pair<int,int> > vorigin, eorigin;
     for ( int i = 0; i < w.isize( ); i++ )
     {    vstart[i] = iv, estart[i] = ie;
          for ( int j = 0; j < w[i].isize( ); j++ )
          {    vorigin.push( i, j );
               eorigin.push( i, j );
               from[iv].push_back(iv+1), to[iv+1].push_back(iv);
               from_edge_obj[iv].push_back(ie);
               to_edge_obj[iv+1].push_back(ie);
               edges[ie++] = w[i][j];
               iv++;    }
          vorigin.push( i, w[i].size( ) );
          iv++;    }
     digraphE<int> W( from, to, edges, to_edge_obj, from_edge_obj );

     // From the overlaps, deduce equivalence relations on the vertices and edges
     // of W.

     equiv_rel evert(iv), eedge(ie);
     for ( int i1 = 0; i1 < over.isize( ); i1++ )
     for ( int j = 0; j < over[i1].isize( ); j++ )
     {    int i2 = over[i1][j].first, o = over[i1][j].second;
          for ( int l1 = 0; l1 <= w[i1].isize( ); l1++ )
          {    int l2 = l1 - o;
               if ( l2 < 0 || l2 > w[i2].isize( ) ) continue;
               if ( l1 == w[i1].isize( ) || l2 == w[i2].isize( ) )
                    evert.Join( vstart[i1] + l1, vstart[i2] + l2 );
               if ( l1 == w[i1].isize( ) || l2 == w[i2].isize( ) ) continue;
               if ( l1 == w[i1].isize( ) ) continue;
               evert.Join( vstart[i1] + l1, vstart[i2] + l2 );
               eedge.Join( estart[i1] + l1, estart[i2] + l2 );    }    }

     // In cases where the resulting graph would have sources that emanate long
     // the same edge, join them.  Ditto for sinks.

     {    vec<int> vreps, ereps;
          evert.OrbitRepsAlt(vreps), eedge.OrbitRepsAlt(ereps);
          int nvreps = vreps.size( ), nereps = ereps.size( );
          vec< vec<int> > from2(nvreps), to2(nvreps); 
          vec< vec<int> > from_edge_obj2(nvreps), to_edge_obj2(nvreps);
          vec<int> edges2(nereps);
          for ( int v1 = 0; v1 < W.N( ); v1++ )
          for ( int j = 0; j < W.From(v1).isize( ); j++ )
          {    int v2 = W.From(v1)[j];
               int e = W.EdgeObjectIndexByIndexFrom( v1, j );
               int pv1 = BinPosition( vreps, evert.ClassId(v1) );
               int pv2 = BinPosition( vreps, evert.ClassId(v2) );
               int pe = BinPosition( ereps, eedge.ClassId(e) );
               if ( Member( from_edge_obj2[pv1], pe ) ) continue;
               from2[pv1].push_back(pv2), to2[pv2].push_back(pv1);
               from_edge_obj2[pv1].push_back(pe); 
               to_edge_obj2[pv2].push_back(pe);    }
          for ( int i = 0; i < nereps; i++ )
          {    int x = w[ eorigin[ ereps[i] ].first ][ eorigin[ ereps[i] ].second ];
               edges2[i] = x;    }
          for ( int i = 0; i < nvreps; i++ )
          {    SortSync( from2[i], from_edge_obj2[i] );
               SortSync( to2[i], to_edge_obj2[i] );    }
          digraphE<int> D( from2, to2, edges2, to_edge_obj2, from_edge_obj2 );
          vec<int> sources, sinks;
          D.Sources(sources), D.Sinks(sinks);
          vec< pair<int,int> > source_edges, sink_edges;
          for ( int i = 0; i < sources.isize( ); i++ )
          {    int v = sources[i];
               for ( int j = 0; j < D.From(v).isize( ); j++ )
                    source_edges.push( D.EdgeObjectByIndexFrom( v, j ), v );    }
          for ( int i = 0; i < sinks.isize( ); i++ )
          {    int v = sinks[i];
               for ( int j = 0; j < D.To(v).isize( ); j++ )
                    sink_edges.push( D.EdgeObjectByIndexTo( v, j ), v );    }
          UniqueSort(source_edges), UniqueSort(sink_edges);
          for ( int i = 0; i < source_edges.isize( ); i++ )
          {    int j;
               for ( j = i + 1; j < source_edges.isize( ); j++ )
                    if ( source_edges[j].first != source_edges[i].first ) break;
               for ( int k = i + 1; k < j; k++ )
               {    int v1 = source_edges[i].second, v2 = source_edges[k].second;
                    evert.Join( vreps[v1], vreps[v2] );    }
               i = j - 1;    }
          for ( int i = 0; i < sink_edges.isize( ); i++ )
          {    int j;
               for ( j = i + 1; j < sink_edges.isize( ); j++ )
                    if ( sink_edges[j].first != sink_edges[i].first ) break;
               for ( int k = i + 1; k < j; k++ )
               {    int v1 = sink_edges[i].second, v2 = sink_edges[k].second;
                    evert.Join( vreps[v1], vreps[v2] );    }
               i = j - 1;    }    }

     // Zipper.
     while(1)
     {    int joins = 0;
          vec<int> vreps, ereps;
          evert.OrbitRepsAlt(vreps), eedge.OrbitRepsAlt(ereps);
          int nvreps = vreps.size( ), nereps = ereps.size( );
          vec< vec<int> > from2(nvreps), to2(nvreps); 
          vec< vec<int> > from_edge_obj2(nvreps), to_edge_obj2(nvreps);
          vec< vec<int> > edges2(nereps);
          for ( int v1 = 0; v1 < W.N( ); v1++ )
          for ( int j = 0; j < W.From(v1).isize( ); j++ )
          {    int v2 = W.From(v1)[j];
               int e = W.EdgeObjectIndexByIndexFrom( v1, j );
               int pv1 = BinPosition( vreps, evert.ClassId(v1) );
               int pv2 = BinPosition( vreps, evert.ClassId(v2) );
               int pe = BinPosition( ereps, eedge.ClassId(e) );
               if ( Member( from_edge_obj2[pv1], pe ) ) continue;
               from2[pv1].push_back(pv2), to2[pv2].push_back(pv1);
               from_edge_obj2[pv1].push_back(pe); 
               to_edge_obj2[pv2].push_back(pe);    }
          for ( int i = 0; i < nereps; i++ )
          {    vec<int> x(1);
               x[0] = w[ eorigin[ ereps[i] ].first ][ eorigin[ ereps[i] ].second ];
               edges2[i] = x;    }
// the execution order actually matters downstream, and the region is critical'ed anyways
//          #pragma omp parallel for
          for ( int v = 0; v < nvreps; v++ )
          {    vec< triple<int,int,int> > s;
               for ( int j = 0; j < from2[v].isize( ); j++ )
               {    int ei = from_edge_obj2[v][j];
                    s.push( edges2[ei][0], from2[v][j], j );    }
               Sort(s);
               for ( int j = 0; j < s.isize( ); j++ )
               {    int k;
                    for ( k = j + 1; k < s.isize( ); k++ )
                         if ( s[k].first != s[j].first ) break;
                    for ( int l = j + 1; l < k; l++ )
                    {    
//                         #pragma omp critical
                         {    evert.Join( vreps[ s[l].second ], 
                                   vreps[ s[j].second ] );
                              eedge.Join( ereps[ from_edge_obj2[v][ s[l].third ] ], 
                                   ereps[ from_edge_obj2[v][ s[j].third ] ] );
                              joins++;    }    }
                    j = k - 1;    }    }
// the execution order actually matters downstream, and the region is critical'ed anyways
//          #pragma omp parallel for
          for ( int v = 0; v < nvreps; v++ )
          {    vec< triple<int,int,int> > s;
               for ( int j = 0; j < to2[v].isize( ); j++ )
               {    int ei = to_edge_obj2[v][j];
                    s.push( edges2[ei][0], to2[v][j], j );    }
               Sort(s);
               for ( int j = 0; j < s.isize( ); j++ )
               {    int k;
                    for ( k = j + 1; k < s.isize( ); k++ )
                         if ( s[k].first != s[j].first ) break;
                    for ( int l = j + 1; l < k; l++ )
                    {    
//                         #pragma omp critical
                         {    evert.Join( vreps[ s[l].second ], 
                                   vreps[ s[j].second ] );
                              eedge.Join( ereps[ to_edge_obj2[v][ s[l].third ] ], 
                                   ereps[ to_edge_obj2[v][ s[j].third ] ] );
                              joins++;    }    }
                    j = k - 1;    }    }
          if ( joins == 0 ) break;    }
     // Apply the equivalence relation to W, yielding the digraph D.

     vec<int> vreps, ereps;
     evert.OrbitRepsAlt(vreps), eedge.OrbitRepsAlt(ereps);
     int nvreps = vreps.size( ), nereps = ereps.size( );
     vec< vec<int> > from2(nvreps), to2(nvreps); 
     vec< vec<int> > from_edge_obj2(nvreps), to_edge_obj2(nvreps);
     vec< vec<int> > edges2(nereps);
     for ( int v1 = 0; v1 < W.N( ); v1++ )
     for ( int j = 0; j < W.From(v1).isize( ); j++ )
     {    int v2 = W.From(v1)[j];
          int e = W.EdgeObjectIndexByIndexFrom( v1, j );
          int pv1 = BinPosition( vreps, evert.ClassId(v1) );
          int pv2 = BinPosition( vreps, evert.ClassId(v2) );
          int pe = BinPosition( ereps, eedge.ClassId(e) );
          if ( Member( from_edge_obj2[pv1], pe ) ) continue;
          from2[pv1].push_back(pv2), to2[pv2].push_back(pv1);
          from_edge_obj2[pv1].push_back(pe), to_edge_obj2[pv2].push_back(pe);    }
     for ( int i = 0; i < nereps; i++ )
     {    vec<int> x(1);
          x[0] = w[ eorigin[ ereps[i] ].first ][ eorigin[ ereps[i] ].second ];
          edges2[i] = x;    }
     for ( int i = 0; i < nvreps; i++ )
     {    SortSync( from2[i], from_edge_obj2[i] );
          SortSync( to2[i], to_edge_obj2[i] );    }
     D.Initialize( from2, to2, edges2, to_edge_obj2, from_edge_obj2 );

     // Compute involution.

     vec<int> inv2( D.EdgeObjectCount( ), -1 );
     for ( int e = 0; e < D.EdgeObjectCount( ); e++ )
     {    int i = eorigin[ ereps[e] ].first, pos = eorigin[ ereps[e] ].second;
          const vec<int>& y = w[i];
          vec<int> yinv;
          if ( inv[ y[0] ] < 0 ) continue;
          for ( int j = 0; j < y.isize( ); j++ )
               yinv.push_back( inv[ y[j] ] );
          yinv.ReverseMe( );
          int posinv = y.isize( ) - pos - 1;
          int iinv = Position( w, yinv );
          ForceAssertGe( iinv, 0 );
          int einv = BinPosition( ereps, eedge.ClassId( estart[iinv] + posinv ) );
          ForceAssertGe( einv, 0 );
          inv2[e] = einv;    }
     inv = inv2;
     REPORT_TIME( clock, "used in WordsToDigraph" );    }
