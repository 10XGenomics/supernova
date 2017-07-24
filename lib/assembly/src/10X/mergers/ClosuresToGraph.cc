// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>

#include "10X/mergers/ClosuresToGraph.h"
// Convert a digraphE<int> into a digraphE<vec<int>> in which unneeded vertices
// have been removed.

void Vectorify( 

     // input: pointer to digraphE<int> -- GETS DELETED MIDSTREAM!!

     digraphE<int>* D0p, 
     
     // input and output:

     vec<int>& dinv, 

     // output:

     digraphE<vec<int>>& D,

     // control:

     const Bool verbose, const Bool single, const Bool use_inv )
{    
     int nthreads = ( single ? 1 : omp_get_max_threads( ) );
     vec<int> to_left, to_right, edges0;
     vec<int> dinv2;
     {    digraphE<int>& D0 = *D0p;
          if (verbose)
          {    cout << Date( ) << ": converting, peak mem = "
                    << PeakMemUsageGBString( ) <<", mem = "<<MemUsageGBString() << endl;    }
          D0.ToLeft(to_left), D0.ToRight(to_right);
          vec<Bool> used( D0.E( ), False );
          for ( int e = 0; e < D0.E( ); e++ )
          {    Bool lverbose = False;
               if ( used[e] ) continue;
               if ( lverbose ) cout << "past continue for e=" << e << endl;
               vec<int> x = {e};
               // used[e] = True; // true but pointless
               Bool circle = False;
               while(1)
               {    int e = x.back( );
                    int v = to_right[e];
                    if ( !D0.From(v).solo( ) || !D0.To(v).solo( ) ) break;
                    int f = D0.IFrom( v, 0 );
                    if ( f == x[0] )
                    {    circle = True;
                         break;    }
                    x.push_back(f);
                    used[f] = True;    }
               if ( lverbose ) 
               {    cout << "grew forward for e=" << e << ", x=" << printSeq(x) 
                         << endl;
                    cout << "circle is " << (circle?"True":"False") << endl;    }
               if ( !circle )
               {    x.ReverseMe( );
                    while(1)
                    {    int e = x.back( );
                         int v = to_left[e];
                         if ( !D0.From(v).solo( ) || !D0.To(v).solo( ) ) break;
                         int f = D0.ITo( v, 0 );
                         if ( f == x[0] )
                         {    circle = True;
                              break;    }
                         x.push_back(f);
                         used[f] = True;    }
                    x.ReverseMe( );    }
               if ( lverbose )
               {    cout << "grew backward for e=" << e << ", x=" << printSeq(x) 
                         << endl;
                    cout << "circle is " << (circle?"True":"False") << endl;    }
               int E = D.E( );
               D.EdgesMutable( ).push_back(x);
               if (use_inv)
               {    vec<int> rx( x.size( ) );
                    for ( int i = 0; i < x.isize( ); i++ )
                    {    rx[ x.isize( ) - i - 1 ] = dinv[ x[i] ];
                         if ( lverbose ) PRINT3( i, x[i], dinv[x[i]] );    }
                    if ( lverbose ) cout << "rx is " << printSeq(rx) << endl;
                    for ( int i = 0; i < rx.isize( ); i++ )
                         used[ rx[i] ] = True;
                    int idx = E;    // index of the one we pushed on
                    if ( x != rx ) 
                    {    D.EdgesMutable().push_back(rx);
                         dinv2.push_back( idx+1 );    }
                    dinv2.push_back( idx );
                    AssertEq(D.E(), dinv2.isize() );    }    }
          if (verbose)
          {    cout << Date( ) << ": copying edges, peak mem = "
                    << PeakMemUsageGBString( ) <<", mem = "<<MemUsageGBString() << endl;    }
          // Note that the following copy seems "theoretically" unnecessary:
          edges0 = D0.Edges( );    
          delete D0p;    }
     if (use_inv) dinv = dinv2;
     if (verbose)
     {    cout << Date( ) << ": making verts, peak mem = "
               << PeakMemUsageGBString( ) <<", mem = "<<MemUsageGBString() << endl;    }
     vec<int> verts( 2 * D.E( ) );
     #pragma omp parallel for schedule(dynamic, 10000) num_threads(nthreads)
     for ( int i = 0; i < D.E( ); i++ )
     {    verts[2*i] = to_left[ D.O(i).front( ) ];
          verts[2*i+1] = to_right[ D.O(i).back( ) ];    }
     if (verbose)
     {    cout << Date( ) << ": sorting verts, peak mem = "
               << PeakMemUsageGBString( ) <<", mem = "<<MemUsageGBString() << endl;    }
     if (single) UniqueSort(verts);
     else ParallelUniqueSort(verts);
     int NV = verts.size( );
     if (verbose)
     {    cout << Date( ) << ": resizing, peak mem = " << PeakMemUsageGBString( )
               <<", mem = "<<MemUsageGBString() << endl;    }
     D.FromMutable( ).resize(NV), D.ToMutable( ).resize(NV);
     D.FromEdgeObjMutable( ).resize(NV), D.ToEdgeObjMutable( ).resize(NV);
     if (verbose)
     {    cout << Date( ) << ": filling out graph, peak mem = "
               << PeakMemUsageGBString( ) <<", mem = "<<MemUsageGBString() << endl;    }
     for ( int e = 0; e < D.E( ); e++ )
     {    int v = BinPosition( verts, to_left[ D.O(e).front( ) ] );
          int w = BinPosition( verts, to_right[ D.O(e).back( ) ] );
          D.FromMutable(v).push_back(w), D.ToMutable(w).push_back(v);
          D.FromEdgeObjMutable(v).push_back(e);
          D.ToEdgeObjMutable(w).push_back(e);    }
     Destroy(to_left), Destroy(to_right), Destroy(verts);
     if (verbose)
     {    cout << Date( ) << ": fixing edges, peak mem = "
               << PeakMemUsageGBString( ) <<", mem = "<<MemUsageGBString() << endl;    }
     #pragma omp parallel for schedule(dynamic, 10000) num_threads(nthreads)
     for ( int e = 0; e < D.E( ); e++ )
     {    for ( int j = 0; j < D.O(e).isize( ); j++ )
               D.OMutable(e)[j] = edges0[ D.O(e)[j] ];    }
     Destroy(edges0);
     if (verbose)
     {    cout << Date( ) << ": completing D construction, peak mem = "
               << PeakMemUsageGBString( ) <<", mem = "<<MemUsageGBString() << endl;    }
     #pragma omp parallel for schedule(dynamic, 10000) num_threads(nthreads)
     for ( int v = 0; v < D.N( ); v++ )
     {    SortSync( D.FromMutable(v), D.FromEdgeObjMutable(v) );
          SortSync( D.ToMutable(v), D.ToEdgeObjMutable(v) );    }

     // Done.

     if (verbose)
          cout << Date( ) << ": final graph has " << D.E( ) << " edges" << endl;    }

template <class T, class EL>
void GetMatches(vec<vec<om<EL>>>& omatch, const HyperBasevectorX& hb,
        const vec<vec<T>>& all_closures, const vec<vec<int>>& ci, 
        const vec<int64_t>& cinv, Bool verbose)
{
    int64_t N = all_closures.size();
    omatch.clear();
    omatch.resize(N); // { (c2,start1,start2,len) }
    if (verbose){
        cout << Date( ) << ": finding short overlaps; peak mem = "
           << PeakMemUsageGBString( )<< ", mem = "<<MemUsageGBString( ) << endl;  
    }

    const int MIN_OVER = 200 - (hb.K()-1);
    #pragma omp parallel for schedule( dynamic, 1000 )
    for ( int i1 = 0; i1 < N; i1++ ) 
    {    // take a closure x1
        const vec<int>& x1 = all_closures[i1];
        int nkmers = 0, b = -1, best = 1000000000;
        // find leftmost least multiplicity index b in x1 in last 200 base seq of x1
        for ( int j = x1.isize( ) - 1; j >= 0; j-- )
        {    if ( ci[ x1[j] ].isize( ) < best )
           {    best = ci[ x1[j] ].size( );
                b = j;    }
           nkmers += hb.Kmers( x1[j] );
           if ( nkmers >= MIN_OVER ) break;    }
        // get all closures sharing this least multiplicity edge
        const vec<int>& ids = ci[ x1[b] ];
        for ( int u = 0; u < ids.isize( ); u++ )
        {    int i2 = ids[u];
          // get edge list corresponding to each such closure
           const vec<int>& x2 = all_closures[i2];
           int j1 = b;
           int stop2 = Min( x2.isize( ), j1 + x2.isize( ) - x1.isize( ) ); // always latter
           for ( int j2 = 0; j2 < stop2; j2++ )
           {    // find first index in x2 that matches least multiplicity edge
                if ( x1[j1] != x2[j2] ) continue;
                int over = 0;
                om<EL> ShortMatch(i2,j1,j2,1);
                ShortMatch.template Extend<T>(x1,x2);
                if (ShortMatch.Stop1()<x1.isize())
                    continue;
                if(ShortMatch.Start2()>0 && ShortMatch.Start1()>0)
                    continue;
                for (int in=ShortMatch.Start1(); in<ShortMatch.Stop1(); in++)
                    over += hb.Kmers(x1[in]);
                if(over<MIN_OVER)
                    continue;
                omatch[i1].push(ShortMatch);
           }    
        }    
    }

    // symmetrize
    if (verbose){
        cout << Date( ) << ": symmetrizing short overlaps; peak mem = "
           << PeakMemUsageGBString( )<< ", mem = "<<MemUsageGBString( ) << endl;  }

    vec<int> oms(N);
    for ( int64_t i = 0; i < N; i++ )
      oms[i] = omatch[i].size( );
    for ( int64_t i1 = 0; i1 < N; i1++ ){
        for ( int j = 0; j < oms[i1]; j++ )
        {   int64_t i2 = omatch[i1][j].c2;
            int start1 = omatch[i1][j].start1, start2 = omatch[i1][j].start2;
            int len = omatch[i1][j].len;
            om<EL> o(i1,start2,start1,len);
            if(!Member(omatch[i2],o))
               omatch[i2].push(o);    
        }
    }

    // Add all long matches.
    if (verbose){
        cout << Date( ) << ": adding long matches + symmetrize; peak mem = "
           << PeakMemUsageGBString( ) << ", mem = "<<MemUsageGBString( ) << endl;    }

    const int batch = 10000;
    #pragma omp parallel for schedule( dynamic, 1 )
    for ( int bi = 0; bi < hb.E( ); bi += batch )
    {    vec< pair<int,int> > Q;
        // start with 10000 edges at a time
        for ( int e = bi; e < Min( bi + batch, hb.E( ) ); e++ )
        {    
            // if edge has not enough kmers or not shared, roll over
            if ( hb.Kmers(e) < MIN_OVER ) continue;
            if ( ci[e].size( ) <= 1 ) continue;

            // Find all closure,location pairs for edge at location in closure
            Q.clear( );
            for ( int j = 0; j < ci[e].isize( ); j++ )
            {   int c = ci[e][j];
                const vec<int>& x = all_closures[c];
                for ( int m = 0; m < x.isize( ); m++ )
                     if ( x[m] == e ) Q.push( c, m );    }

            // Join. Take a pair from Q
            for ( int i1 = 0; i1 < Q.isize( ); i1++ )
            {   int c1 = Q[i1].first, m1 = Q[i1].second;
               // take another pair from Q not already taken
                for ( int i2 = i1 + 1; i2 < Q.isize( ); i2++ )
                {   int c2 = Q[i2].first, m2 = Q[i2].second;
                    Bool match = False;
                    // continue to next loop if these two closures already have been matched
                    for ( int l = 0; l < omatch[c1].isize( ); l++ )
                    {    if ( omatch[c1][l].c2 != c2 ) continue;
                         if ( omatch[c1][l].Offset( ) != m1 - m2 ) continue;
                         if ( m1 < omatch[c1][l].Start1( ) ) continue;
                         if ( m1 >= omatch[c1][l].Stop1( ) ) continue;
                         match = True;
                         break;    }
                    if (match) continue;
                    for ( int l = 0; l < omatch[c2].isize( ); l++ )
                    {    if ( omatch[c2][l].c2 != c1 ) continue;
                         if ( omatch[c2][l].Offset( ) != m2 - m1 ) continue;
                         if ( m2 < omatch[c2][l].Start1( ) ) continue;
                         if ( m2 >= omatch[c2][l].Stop1( ) ) continue;
                         match = True;
                         break;    }
                    if (match) continue;
                    // otherwise we have a new distinct long-long match
                    #pragma omp critical
                    {    omatch[c1].push( c2, m1, m2, 1 );
                        // we glue these together at edge, and extend before and after
                         omatch[c1].back( ).template Extend<T>(
                              all_closures[c1], all_closures[c2] );
                         om<EL> o(c1, omatch[c1].back().Start2(), omatch[c1].back().Start1(),
                                 omatch[c1].back().Len());
                         if(!Member(omatch[c2],o))
                             omatch[c2].push(o);
                    }    
                }    
            }    
        }    
    } 

    // Force match symmetry with respect to involution.
    if (verbose){
        cout << Date( ) << ": symmetrizing wrt involution; peak mem = "
           << PeakMemUsageGBString( ) << ", mem = "<<MemUsageGBString( ) << endl;    }

    for ( int64_t i = 0; i < N; i++ )
        oms[i] = omatch[i].size( );

    int64_t iadds = 0, itotal = 0;
    double iclock = WallClockTime( );
    for ( int64_t i1 = 0; i1 < N; i1++ )
    {    itotal += oms[i1];
        for ( int j = 0; j < oms[i1]; j++ )
        {   int64_t i2 = omatch[i1][j].c2;
            int start1 = omatch[i1][j].start1, start2 = omatch[i1][j].start2;
            int len = omatch[i1][j].len;
            int64_t ip1 = cinv[i1], ip2 = cinv[i2];
            int istart1 = all_closures[i1].isize( ) - start1 - len;
            int istart2 = all_closures[i2].isize( ) - start2 - len;
            om<EL> o( ip2, istart1, istart2, len );
            if ( !Member( omatch[ip1], o ) ){
                omatch[ip1].push_back(o);
                iadds++;    
            }
        }
    }
    if (verbose)
    {   DPRINT2_TO( cout, iadds, itotal );
        cout << Date( ) << ": done, time used = " << TimeSince(iclock)
           << endl;    }
    Destroy(oms);
}

template <class T>
void ClosuresToGraph( const HyperBasevectorX& hb, const vec<int>& inv,
     vec<vec<T>>& all_closures, digraphE<vec<T>>& D, vec<int>& dinv,
     const Bool verbose, String dir, Bool Canon)
{

    int64_t hbE = hb.E();
    vec<int64_t> cinv;
    vec<vec<int >> ci;
    vec<vec<om<uint16_t>>> omatch;

    // closure indexing
    ComputeClosureIndex<T>(inv, cinv, ci, hbE, all_closures, verbose); 

    // get matches
    GetMatches<T,uint16_t>(omatch, hb, all_closures, ci, cinv, verbose);  

    // create graph
    Bool GET_META = False;
    vec<vec<quad<int,uint16_t,int,uint16_t>>> closure_paths;
    NucleateGraph<T,uint16_t>(omatch, hbE, all_closures, cinv, ci, D, dinv, 
            GET_META, closure_paths, verbose, dir, Canon);
}

template void ClosuresToGraph( const HyperBasevectorX& hb, const vec<int>& inv,
     vec<vec<int>>& all_closures, digraphE<vec<int>>& D, vec<int>& dinv,
     const Bool verbose, String dir, Bool Canon);

template <class T>
void ClosuresToGraph_HP( const HyperBasevectorX& hb, const vec<int>& inv,
     vec<vec<T>>& all_closures, digraphE<vec<T>>& D, vec<int>& dinv,
     const Bool verbose, String dir, Bool Canon)
{

    int64_t hbE = hb.E();
    vec<int64_t> cinv;
    vec<vec<int >> ci;
    vec<vec<om<uint16_t>>> omatch;

    // closure indexing
    ComputeClosureIndex<T>(inv, cinv, ci, hbE, all_closures, verbose); 

    // get matches
    GetMatches<T,uint16_t>(omatch, hb, all_closures, ci, cinv, verbose);  

    // create graph
    Bool GET_META = False;
    vec<vec<quad<int,uint16_t,int,uint16_t>>> closure_paths;
    NucleateGraph_HP<T,uint16_t>(omatch, hbE, all_closures, cinv, ci, D, dinv, 
            GET_META, closure_paths, verbose, dir, Canon);
}

template void ClosuresToGraph_HP( const HyperBasevectorX& hb, const vec<int>& inv,
     vec<vec<int>>& all_closures, digraphE<vec<int>>& D, vec<int>& dinv,
     const Bool verbose, String dir, Bool Canon);
