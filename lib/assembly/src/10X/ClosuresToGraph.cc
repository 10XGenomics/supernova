// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>

#include "CoreTools.h"
#include "ParallelVecUtilities.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "10X/ClosuresToGraph.h"
#include "10X/Super.h"
#include "10X/MaybeWriter.h"

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

static void InvertClosures(const vec<vec<int>>& all_closures, vec<vec<int> >& ci, const int64_t nE){
    ci.resize(nE);
    for ( int64_t i = 0; i < (int64_t) all_closures.size(); i++ )
    for ( int j = 0; j < all_closures[i].isize( ); j++ ){
        ci[ all_closures[i][j] ].push_back(i);
    }

    #pragma omp parallel for schedule( dynamic, 1000 )
    for ( int e = 0; e < nE; e++ )
    {   UniqueSort( ci[e] );}
}

template <class T>
static void SwapF(T& a, T& b){
    T c = a;
    a = b;
    b = c;
}

namespace Aliasing{
    // use only with primitive data types 

    void JoinClasses(vec<etuple>& e_inst, const vec<int64_t>& estart, 
            int64_t idx1, int64_t idx2){
        if(e_inst[idx1].eoID==e_inst[idx2].eoID) return;
        SwapF(e_inst[idx1].c,e_inst[idx2].c); 
        SwapF(e_inst[idx1].m,e_inst[idx2].m);
        int64_t idx = estart[e_inst[idx1].c] + e_inst[idx1].m;
        while(e_inst[idx].eoID!=e_inst[idx1].eoID){
            e_inst[idx].eoID = e_inst[idx1].eoID;
            idx = estart[e_inst[idx].c] + e_inst[idx].m;
        }
    }

    void JoinVerts(vec<vtuple>& v_inst, const int idx1, const int idx2){
        if(v_inst[idx1].voID == v_inst[idx2].voID) return; // nothing to do
        SwapF(v_inst[idx1].nv,v_inst[idx2].nv);
        int idx = v_inst[idx1].nv; 
        while(v_inst[idx].voID!=v_inst[idx1].voID){
            v_inst[idx].voID = v_inst[idx1].voID;
            idx = v_inst[idx].nv; 
        }
    }

}

namespace BuildSuperGraph{

    void IsLeftMost(vec<etuple>& e_inst, const vec<int64_t>& estart, 
            const vec<pair<int,unsigned short>>& orb, const vec<vec<int>>& all_closures){ 
        // cycle through given orbit{c.m} of e and choose c with m > 0
        int allStarting = 1;
        unsigned int pr;
        std::set<pair<int,unsigned int>> prevsInNonStarting;
        for (unsigned int p = 0; p<orb.size(); p++){
            if(orb[p].second>0){
                allStarting = 0;
                prevsInNonStarting.insert(make_pair(all_closures[orb[p].first][orb[p].second-1],
                            e_inst[estart[orb[p].first]+orb[p].second-1].eoID));
                pr = p;
            }
        }

        // dangling start means leftmost
        if(allStarting){
            for(const auto& p:orb)
                e_inst[estart[p.first]+p.second].type = 1;
            return;
        }

        if(prevsInNonStarting.size()>1){
            for(const auto& p:orb)
                e_inst[estart[p.first]+p.second].type = 1;
            return;
        }
        

        // get orbit of e-- containing (closure,instance) connecting with special e
        int64_t idx;
        vec<pair<int,unsigned short>> Q;
        Q.push(orb[pr].first,orb[pr].second-1);
        int64_t start = estart[Q.back().first]+Q.back().second;
        while(1){
            idx = estart[Q.back().first]+Q.back().second;
            if(estart[e_inst[idx].c]+e_inst[idx].m == start)
                break;
            Q.push(e_inst[idx].c,e_inst[idx].m);
        }

        // consider non-ending closures on e-- and get set of nexts
        std::set<pair<int,unsigned int>> nextsInNonEnding;
        for(pr = 0; pr<Q.size(); pr++){
            if(Q[pr].second < all_closures[Q[pr].first].size()-1)
                nextsInNonEnding.insert(make_pair(all_closures[Q[pr].first][Q[pr].second+1],
                            e_inst[estart[Q[pr].first]+Q[pr].second+1].eoID));
        }

        if(nextsInNonEnding.size()>1){
            for(const auto& p:orb)
                e_inst[estart[p.first]+p.second].type = 1;
            return;
        }

        return;
    }

    void IsRightMost(vec<etuple>& e_inst, const vec<int64_t>& estart, 
            const vec<pair<int,unsigned short>>& orb, const vec<vec<int>>& all_closures){
        // same as detecting leftmost
        int allEnding = 1;
        unsigned int pr;
        std::set<pair<int,unsigned int>> nextsInNonEnding;
        for (unsigned int p = 0; p<orb.size(); p++){
            if(orb[p].second < all_closures[orb[p].first].size()-1){
                allEnding = 0;
                nextsInNonEnding.insert(make_pair(all_closures[orb[p].first][orb[p].second+1],
                            e_inst[estart[orb[p].first]+orb[p].second+1].eoID));
                pr = p;
            }
        }

        // dangling end means rightmost
        if(allEnding){
            for(const auto& p:orb)
                e_inst[estart[p.first]+p.second].type = 
                    (e_inst[estart[p.first]+p.second].type>0)? 3:2;
            return;
        }

        if(nextsInNonEnding.size()>1){
            for(const auto& p:orb)
                e_inst[estart[p.first]+p.second].type = 
                    (e_inst[estart[p.first]+p.second].type>0)? 3:2;
            return;
        }

        // get orbit of e++ containing (closure,instance) connecting with special e
        int64_t idx;
        vec<pair<int,unsigned short>> Q;
        Q.push(orb[pr].first,orb[pr].second+1);
        int64_t start = estart[Q.back().first]+Q.back().second;
        while(1){
            idx = estart[Q.back().first]+Q.back().second;
            if(estart[e_inst[idx].c]+e_inst[idx].m == start)
                break;
            Q.push(e_inst[idx].c,e_inst[idx].m);
        }

        // consider  non-ending closures on e++ 
        std::set<pair<int,unsigned int>> prevsInNonStarting;
        for(pr = 0; pr<Q.size(); pr++){
            if(Q[pr].second > 0)
                prevsInNonStarting.insert(make_pair(all_closures[Q[pr].first][Q[pr].second-1],
                            e_inst[estart[Q[pr].first]+Q[pr].second-1].eoID));
        }

        if(prevsInNonStarting.size()>1){
            for(const auto& p:orb)
                e_inst[estart[p.first]+p.second].type = 
                    (e_inst[estart[p.first]+p.second].type>0)? 3:2;
            return;
        }

        return;
    }

    void SmartJoinR2L(vec<vtuple>& v_inst, const vec<etuple>&  e_inst, 
            const vec<int64_t>& estart, const int i,
            std::unordered_map<int64_t,pair<int,int>>& map_ei)
    {
        // get eorbit for v_inst[i]:
        vec<pair<int,int64_t>> iorb;
        int64_t ne = v_inst[i].v;
        int64_t ni = ne;
        while(1){
            ne = estart[e_inst[ni].c]+e_inst[ni].m; 
            iorb.push(e_inst[ni].c,ne);
            if(ne==v_inst[i].v) break;
            ni = ne;
        }

        // get prev edge's index in e_inst and map back to v_inst, 
        // which will be joined to i and removed from verts. this joins rights to a left.
        for(const auto& io : iorb){
            if(io.second>estart[io.first]){
                // find orbit of io.second-1
                int64_t pe = io.second-1;
                while(1){
                    pe = estart[e_inst[pe].c]+e_inst[pe].m; 
                    if(map_ei.count(pe)){
                        #pragma omp critical
                        {Aliasing::JoinVerts(v_inst, i, map_ei[pe].second);}
                        break;
                    }
                    if(pe==io.second-1) break;
                }
            }
        }
    }

    void SmartJoinL2R(vec<vtuple>& v_inst, const vec<etuple>&  e_inst, 
            const vec<int64_t>& estart, const int i,
            std::unordered_map<int64_t,pair<int,int>>& map_ei)
    {
        // get eorbit for v_inst[i]:
        vec<pair<int,int64_t>> iorb;
        int64_t ne = v_inst[i].v;
        int64_t ni = ne;
        while(1){
            ne = estart[e_inst[ni].c]+e_inst[ni].m; 
            iorb.push(e_inst[ni].c,ne);
            if(ne==v_inst[i].v) break;
            ni = ne;
        }

        // get next edge's index in e_inst and map back to v_inst, which will be 
        // joined to i and removed from verts. this joins rights to a left.
        for(const auto& io : iorb){
            if(io.second < estart[io.first+1]-1){
                // find orbit of io.second+1
                int64_t pe = io.second+1;
                while(1){
                    pe = estart[e_inst[pe].c]+e_inst[pe].m; 
                    if(map_ei.count(pe)){
                        #pragma omp critical
                        {Aliasing::JoinVerts(v_inst, i, map_ei[pe].first);}
                        break;
                    }
                    if(pe==io.second+1) break;
                }
            }
        }
    }


    int RecursiveWalk(vec<int>& unipath, const vec<etuple>& e_inst,
            const vec<int64_t>& estart, const int64_t v,
            const vec<vec<int>>& all_closures,
            std::unordered_map<int64_t,pair<int,int>>& map_ei){
        unsigned int j = 0;
        int64_t nv = v;
        while(j<all_closures[e_inst[nv].c].size()-e_inst[nv].m){
            unipath.push_back(all_closures[e_inst[nv].c][e_inst[nv].m+j]);
            // added element is rightmost if it is type 3 or 2 in einst
            if(e_inst[estart[e_inst[nv].c]+e_inst[nv].m+j].type>1){
                // get eorbit
                vec<pair<int,unsigned short>> Q;
                int64_t ne = estart[e_inst[nv].c]+e_inst[nv].m+j;
                int64_t st = ne;
                while(1){
                    ne = estart[e_inst[ne].c]+e_inst[ne].m; 
                    if(map_ei.count(ne)){
                        return map_ei[ne].second;
                        break;
                    }
                    if(ne==st) break;
                }
            }
            j++;
        }
        // otherwise reached end of closure, get aligning closure and redo above step
        // get singleton eorbit of (e_inst[nv].c,e_inst[nv].m+j-1)
        int64_t idx;
        vec<pair<int,unsigned short>> Q;
        Q.push(e_inst[nv].c,e_inst[nv].m+j-1);
        int64_t start = estart[Q.back().first]+Q.back().second;
        while(1){
            idx = estart[Q.back().first]+Q.back().second;
            if(estart[e_inst[idx].c]+e_inst[idx].m == start)
                break;
            Q.push(e_inst[idx].c,e_inst[idx].m);
        }
        ForceAssertGe(Q.size(),2); // this confirms that there is only one aligning closure
        for(const auto& q:Q){
            if((unsigned short)all_closures[q.first].size()==q.second+1)
                continue;
            nv = estart[q.first]+q.second+1; 
            break;
        }
        return RecursiveWalk(unipath, e_inst, estart, nv, all_closures, map_ei);
    }
}

void ClosuresToGraph( const HyperBasevectorX& hb, const vec<int>& inv,
     vec<vec<int>>& all_closures, digraphE<vec<int>>& D, vec<int>& dinv,
     const Bool verbose, String dir)
{

    // Index closures.
    int64_t N = all_closures.size( ), total = 0;
    for ( int64_t i = 0; i < N; i++ )
      total += all_closures[i].size( );

    if (verbose){
        cout << Date( ) << ": " << ToStringAddCommas(N) << " closures "
           << "having in total " << ToStringAddCommas(total) << " edges, with repetitions from " << hb.E( ) << " distinct edges" <<endl; }

    if (verbose){
        cout << Date( ) << ": indexing closures + uniquesorting, peak mem = "
           << PeakMemUsageGBString( )<< ", mem = "<<MemUsageGBString( ) << endl;    }

    vec<vec<int >> ci( hb.E( ) );
    vec<vec<om>> omatch;

    omatch.resize(N); // { (c2,start1,start2,len) }
    // short overlaps
    InvertClosures(all_closures,ci,hb.E());
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
                om ShortMatch(i2,j1,j2,1);
                ShortMatch.Extend(x1,x2);
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
            om o(i1,start2,start1,len);
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
                         omatch[c1].back( ).Extend(
                              all_closures[c1], all_closures[c2] );
                         om o(c1, omatch[c1].back().Start2(), omatch[c1].back().Start1(),
                                 omatch[c1].back().Len());
                         if(!Member(omatch[c2],o))
                             omatch[c2].push(o);
                    }    
                }    
            }    
        }    
    }

    // Force match symmetry.
    if (verbose){
        cout << Date( ) << ": computing involution for all_closures; peak mem = "
           << PeakMemUsageGBString( ) << ", mem = "<<MemUsageGBString( ) << endl;    }

    vec<int64_t> cinv( N, -1 );
    #pragma omp parallel for
    for ( int64_t i = 0; i < N; i++ ){
        const vec<int>& x = all_closures[i];
        vec<int> rx( x.size( ) );
        for ( int j = 0; j < x.isize( ); j++ )
            rx[ x.isize( ) - j - 1 ] = inv[ x[j] ];
        int64_t ip = BinPosition( all_closures, rx );
        ForceAssertGe( ip, 0 );
        cinv[i] = ip;
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
            om o( ip2, istart1, istart2, len );
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
    // write ci to disk and destroy
    MaybeWriter finalWriter(true);
    if(dir!="")
    {
        finalWriter.writeFile(dir+"/orig/a.ci",ci);
        Destroy(ci);
    // write omatches to disk and destroy
        finalWriter.writeFile(dir+"/orig/a.omatches",omatch);
        Destroy(omatch);
    }

    // define a struct used later
    if (verbose){
        cout << Date( ) << ": forming data structures for supergraph, peak mem = "
           << PeakMemUsageGBString( )<< ", mem = "<<MemUsageGBString( ) << endl;    }

    vec<int64_t> estart(N+1);
    estart[0] = 0;
    for ( int i = 1; i < N+1; i++ )
        estart[i] = estart[i-1]+all_closures[i-1].size();

    vec<etuple> e_inst;
    e_inst.resize(total);

    vec<unsigned int> eorbitIds(hb.E(),0);
    for(int c = 0; c < N; c++){
        // generate pairs (closure containing e, index)
        for ( int m = 0; m < all_closures[c].isize(); m++ ){
            e_inst[estart[c]+m].c = c; 
            e_inst[estart[c]+m].m = m;
            e_inst[estart[c]+m].eoID = eorbitIds[all_closures[c][m]];
            e_inst[estart[c]+m].type = 0;
            eorbitIds[all_closures[c][m]]++;
        }    
    }
    Destroy(eorbitIds);

    // write all_closures to disk and destroy
    if(dir!=""){
        finalWriter.writeFile(dir+"/orig/a.all_closures",all_closures);
        Destroy(all_closures);
        // load omatches
        BinaryReader::readFile(dir+"/orig/a.omatches",&omatch);
    }

    if (verbose)
    {   cout << Date( ) << ": defining local equivalence on edges of basegraph -- caution: this can take many hours!; peak mem = " 
           << PeakMemUsageGBString( ) << ", mem = "<<MemUsageGBString() << endl;    }

    // loop through all closures/words. 
    // This way is fast since each match is looped through once.
    // However we can't make it parallelized
    // Slower way is to do this to consider the set of instances of e : (c,m) and partition them based on omatch info
    int64_t l1, l2;
    for ( int i1 = 0; i1 < omatch.isize( ); i1++ ){
        for ( int j = 0; j < omatch[i1].isize( ); j++ )
        {    
            l1 = estart[i1]+omatch[i1][j].start1; 
            l2 = estart[omatch[i1][j].c2]+ omatch[i1][j].start2;
            for ( int l = 0; l < omatch[i1][j].len; l++ )
                Aliasing::JoinClasses(e_inst,estart, l1 + l, l2 + l);
        }    
    }    
    Destroy(omatch);
    if(dir!=""){
        // load all_closures
        BinaryReader::readFile(dir+"/orig/a.all_closures",&all_closures);
    }

    if (verbose)
    {    cout << Date( ) << ": finding edge type & marking edges to/from vertices; peak mem = " 
           << PeakMemUsageGBString( ) << ", mem = "<<MemUsageGBString() << endl;    }
    if(dir!=""){
        // load ci
        BinaryReader::readFile(dir+"/orig/a.ci",&ci);
    }

    // populate type 
    // 0 if ordinary
    // 1 if leftmost-only
    // 2 if rightmost-only
    // 3 if both left and rightmost
    #pragma omp parallel for schedule( dynamic, 1 )
    for(int64_t e = 0; e < hb.E(); e++){
        // get partitions/orbits of (c,m) wrt e
        vec<unsigned int> eoID;
        vec<vec<pair<int,unsigned short>>> orbits;
        for(int j = 0; j< ci[e].isize(); j++){
            int c = ci[e][j];
            const vec<int>& x = all_closures[c];
            for(int m = 0; m<x.isize(); m++){
                if (x[m] == e){
                    unsigned int thiseoID = e_inst[estart[c]+m].eoID;
                    unsigned int pos = std::find(eoID.begin(),eoID.end(),thiseoID)-eoID.begin();
                    // if not already taken, add to eoID and orbits
                    if( pos == eoID.size()){
                        orbits.push_back(vec<pair<int,unsigned short>>());
                        orbits.back().push(c,m);
                        eoID.push_back(thiseoID);
                    }else{ // add to existing places
                        orbits[pos].push(c,m); 
                    }
                }
            }
        }
        //for each orbit
        for(const auto& orb : orbits){
            BuildSuperGraph::IsLeftMost(e_inst,estart,orb,all_closures);
        } 
    }
    #pragma omp parallel for schedule( dynamic, 1 )
    for(int64_t e = 0; e < hb.E(); e++){
        // get partitions/orbits of (c,m) wrt e
        vec<unsigned int> eoID;
        vec<vec<pair<int,unsigned short>>> orbits;
        for(int j = 0; j< ci[e].isize(); j++){
            int c = ci[e][j];
            const vec<int>& x = all_closures[c];
            for(int m = 0; m<x.isize(); m++){
                if (x[m] == e){
                    unsigned int thiseoID = e_inst[estart[c]+m].eoID;
                    unsigned int pos = std::find(eoID.begin(),eoID.end(),thiseoID)-eoID.begin();
                    // if not already taken, add to eoID and orbits
                    if( pos == eoID.size()){
                        orbits.push_back(vec<pair<int,unsigned short>>());
                        orbits.back().push(c,m);
                        eoID.push_back(thiseoID);
                    }else{ // add to existing places
                        orbits[pos].push(c,m); 
                    }
                }
            }
        }
        //for each orbit
        for(const auto& orb : orbits){
            BuildSuperGraph::IsRightMost(e_inst,estart,orb,all_closures);
        } 
    }
    if (verbose)
    {   cout << Date( ) <<": defining global equivalence on vertices of supergraph; peak mem = " 
           << PeakMemUsageGBString( ) << ", mem = "<<MemUsageGBString() << endl;    }

    // Construct a dense index of vtuple(id in e_inst, cycle vertex, void, type) to e_inst 
    // for all etuple with type > 0. Add type==1 and type==2; type==3 is added once each as type 1 and 2. 
    vec<vtuple> v_inst;
    v_inst.reserve(2*hb.E()); // guessing the size
    std::unordered_map<int64_t,pair<int,int>> map_ei; // build inverse index, can be a vec
    map_ei.reserve(2*hb.E());
    int vorbitIds = 0;
    #pragma omp parallel for ordered schedule( dynamic, 1 )
    for(int64_t e = 0; e < hb.E(); e++){
        // get partitions/orbits of (c,m) wrt e
        vec<unsigned int> eoID;
        vec<vec<pair<int,unsigned short>>> orbits;
        for(int j = 0; j< ci[e].isize(); j++){
            int c = ci[e][j];
            const vec<int>& x = all_closures[c];
            for(int m = 0; m<x.isize(); m++){
                if (x[m] == e){
                    unsigned int thiseoID = e_inst[estart[c]+m].eoID;
                    unsigned int pos = std::find(eoID.begin(),eoID.end(),thiseoID)-eoID.begin();
                    // if not already taken, add to eoID and orbits
                    if( pos == eoID.size()){
                        orbits.push_back(vec<pair<int,unsigned short>>());
                        orbits.back().push(c,m);
                        eoID.push_back(thiseoID);
                    }else{ // add to existing places
                        orbits[pos].push(c,m); 
                    }
                }
            }
        }
        //for each orbit
       #pragma omp ordered
       {
        for(const auto& orb : orbits){
           int64_t ii = estart[orb[0].first]+orb[0].second;
                int turn = 0;
                if(e_inst[ii].type == 3){
                    map_ei.insert(make_pair(ii,make_pair(vorbitIds,vorbitIds+1)));
                    v_inst.push(ii, vorbitIds ,vorbitIds,1);
                    v_inst.push(ii, vorbitIds+1 ,vorbitIds+1,2);
                    vorbitIds+=2;
                    turn = 1;
                }
                if(e_inst[ii].type == 1){
                    map_ei.insert(make_pair(ii,make_pair(vorbitIds,-1)));
                    v_inst.push(ii, vorbitIds, vorbitIds, e_inst[ii].type);
                    vorbitIds++;
                    turn = 1;
                }
                if(e_inst[ii].type == 2){
                    map_ei.insert(make_pair(ii,make_pair(-1,vorbitIds)));
                    v_inst.push(ii, vorbitIds, vorbitIds, e_inst[ii].type);
                    vorbitIds++;
                    turn = 1;
                }
                if(turn){
                    for(const auto& q: orb)
                        e_inst[estart[q.first]+q.second].type *= -1;
                }
            }
        }
    }
    // reset edge type
    for(auto& e : e_inst)
        e.type *= -1;

    Destroy(ci);

    if (verbose)
    {   cout << Date( ) << ": creating representatives for verts; peak mem = " 
           << PeakMemUsageGBString( ) << ", mem = "<<MemUsageGBString() << endl;    }


    // get vorbit reps for left and right each and join them
    vec<int> repL, repR;
    repL.reserve(v_inst.size()/2);
    repR.reserve(v_inst.size()/2);
    for(int i = 0; i< v_inst.isize(); i++){
    if(v_inst[i].type==1)
        repL.push_back(i);
    else
        repR.push_back(i);
    }

    if (verbose)
    {  
        cout<<"number of edge_obj in supergraph: "<<v_inst.size()/2<<endl;
        cout << Date( ) << ": completing vertex aliasing -- caution: this can take many hours!; peak mem = " 
           << PeakMemUsageGBString( ) << ", mem = "<<MemUsageGBString() << endl;    }

    // v1 equiv v2 if they form a successive left-right pair along a common closure
    #pragma omp parallel for
    for(int i = 0; i< repL.isize(); i++){
        BuildSuperGraph::SmartJoinR2L(v_inst, e_inst, estart, repL[i], map_ei);
    }
    #pragma omp parallel for
    for(int i = 0; i< repR.isize(); i++){
        BuildSuperGraph::SmartJoinL2R(v_inst, e_inst, estart, repR[i], map_ei);
    } 

    if (verbose)
    {   cout << Date( ) << ": remapping vertex orbitIDs supergraph; peak mem = " 
           << PeakMemUsageGBString( ) << ", mem = "<<MemUsageGBString() << endl;    }

    // make the vertex orbit IDs to be consecutive, starting from 0
    std::unordered_map<int,int> vmap;
    int ix = 0;
    for(const auto& v : v_inst){
    if(vmap.count(v.voID)==0)
        vmap[v.voID] = ix++;
    }
      
    if (verbose)
    {   
        cout<<"number of vertices in supergraph: "<<vmap.size()<<endl;
        cout << Date( ) << ": constructing edge objects in supergraph; peak mem = " 
           << PeakMemUsageGBString( ) << ", mem = "<<MemUsageGBString() << endl;    }

    // construct edge_obj 
    vec<vec<int>> unipaths(repL.size());
    vec<int> permRepR(repR.size(),-1);
    std::unordered_map<vec<int64_t>,int> Start;
    Start.reserve(unipaths.size());
    // with each leftmost edge
    #pragma omp parallel for
    for(int iv = 0; iv < repL.isize(); iv++){
        permRepR[iv] = BinPosition(repR,BuildSuperGraph::RecursiveWalk(unipaths[iv],
                  e_inst,estart,v_inst[repL[iv]].v,all_closures, map_ei));
        // get eorbit of v_inst[repL[iv]].v as (c,m)
        vec<int64_t> Q;
        int64_t ne = v_inst[repL[iv]].v;
        while(1){
            ne = estart[e_inst[ne].c]+e_inst[ne].m; 
            Q.push(ne);
            if(ne==v_inst[repL[iv]].v) break;
        }
        Sort(Q);
        #pragma omp critical
        {Start[Q] = iv;}
    }
    int64_t SUM = 0;
    for(const auto& upath: unipaths)
        SUM += upath.size();

    if (verbose)
    {   
        cout<<"All edge objects together have "<<SUM<<" edges"<<endl;
        cout << Date( ) << ": constructing involution on the supergraph; peak mem = "
           << PeakMemUsageGBString( ) << ", mem = "<<MemUsageGBString( ) << endl;    }

    vec<int> ws(N,0);
    for ( int64_t i = 0; i < N; i++ )
        ws[i] = all_closures[i].size( );
    Destroy(all_closures); 

    dinv.resize(unipaths.size(),-1);
    #pragma omp parallel for
    for(int p = 0; p < unipaths.isize(); p++){
        // get eorbit of v_inst[repR[permRepR[p]]].v as (c,m)
        vec<pair<int,unsigned short>> Q;
        int64_t ne = v_inst[repR[permRepR[p]]].v;
        while(1){
           Q.push(e_inst[ne].c, e_inst[ne].m);
           ne = estart[e_inst[ne].c]+e_inst[ne].m; 
           if(ne==v_inst[repR[permRepR[p]]].v) break;
        }
        vec<int64_t> RCQ;
        for(const auto& q: Q)
            RCQ.push(estart[cinv[q.first]]+ws[q.first]-q.second-1);
        Sort(RCQ);
        #pragma omp critical
        { dinv[p] = Start[RCQ]; }
    }

    Destroy(ws);
    Destroy(e_inst);
    Destroy(estart);

    if (verbose)
    {   cout << Date( ) << ": constructing supergraph; peak mem = "
           << PeakMemUsageGBString( ) << ", mem = "<<MemUsageGBString( ) << endl;    }

    // construct the digraph
    D.EdgesMutable().resize(unipaths.size());
    #pragma omp parallel for
    for(unsigned int i = 0; i<unipaths.size(); i++){
        D.OMutable(i) = unipaths[i];
    }
    D.FromMutable().resize(vmap.size());
    D.ToMutable().resize(vmap.size());
    D.FromEdgeObjMutable().resize(vmap.size());
    D.ToEdgeObjMutable().resize(vmap.size());

    for(unsigned int i = 0; i<unipaths.size(); i++){
        D.FromMutable(vmap[v_inst[repL[i]].voID]).push_back(vmap[v_inst[repR[permRepR[i]]].voID]);
        D.ToMutable(vmap[v_inst[repR[permRepR[i]]].voID]).push_back(vmap[v_inst[repL[i]].voID]);
        D.FromEdgeObjMutable(vmap[v_inst[repL[i]].voID]).push_back(i);
        D.ToEdgeObjMutable(vmap[v_inst[repR[permRepR[i]]].voID]).push_back(i);
    }

    // sort syncing
    if(verbose) cout<< Date()<<": sortsyncing"<< ",mem = "<<MemUsageGBString()<<endl;
    #pragma omp parallel for
    for(unsigned int i =0; i<vmap.size(); i++){
        SortSync(D.FromMutable(i),D.FromEdgeObjMutable(i));
        SortSync(D.ToMutable(i),D.ToEdgeObjMutable(i));
    }
    // check inversion
    int unset = 0, uninv = 0;
    for ( int i = 0; i < dinv.isize( ); i++ ){
        if ( dinv[i] < 0 ) unset++;
        else if ( dinv[dinv[i]] != i ) uninv++;
    }
    ForceAssertEq( unset, 0 );
    ForceAssertEq( uninv, 0 ); 

    Destroy(v_inst);
    Destroy(repR);
    Destroy(repL);
    Destroy(permRepR);
    vmap.clear();
}
