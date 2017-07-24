// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>

#include "10X/mergers/NucleateGraph_HP.h"
#include "10X/mergers/NicePrints.h"

namespace BuildSuperGraph_HP{
    template <class T, class EL>
    void IsLeftMost(vec<etuple<EL>>& e_inst, const vec<int64_t>& estart, 
            const vec<pair<int,EL>>& orb, const vec<vec<T>>& all_closures,
            const vec<vec<vec<pair<int,EL>>>>& eorbits){ 
        // cycle through given orbit{c.m} of e and choose c with m > 0
        int allStarting = 1;
        unsigned int pr;
        std::set<pair<T,unsigned int>> prevsInNonStarting;
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
        T emm = all_closures[orb[pr].first][orb[pr].second-1];
        const vec<pair<int,EL>>& Q = eorbits[emm][e_inst[estart[orb[pr].first]+orb[pr].second-1].eoID];

        // consider non-ending closures on e-- and get set of nexts
        std::set<pair<T,unsigned int>> nextsInNonEnding;
        for(pr = 0; pr<Q.size(); pr++){
            if(Q[pr].second < all_closures[Q[pr].first].size()-1)
                nextsInNonEnding.insert(make_pair(
                            all_closures[Q[pr].first][Q[pr].second+1],
                            e_inst[estart[Q[pr].first]+Q[pr].second+1].eoID));
            if(nextsInNonEnding.size()>1)
                break;
        }

        if(nextsInNonEnding.size()>1){
            for(const auto& p:orb)
                e_inst[estart[p.first]+p.second].type = 1;
            return;
        }

        return;
    }

    template <class T, class EL>
    void IsRightMost(vec<etuple<EL>>& e_inst, const vec<int64_t>& estart, 
            const vec<pair<int,EL>>& orb, const vec<vec<T>>& all_closures,
            const vec<vec<vec<pair<int,EL>>>>& eorbits){
        // same as detecting leftmost
        int allEnding = 1;
        unsigned int pr;
        std::set<pair<T,unsigned int>> nextsInNonEnding;
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
        T epp = all_closures[orb[pr].first][orb[pr].second+1];
        const vec<pair<int,EL>>& Q = eorbits[epp][e_inst[estart[orb[pr].first]+orb[pr].second+1].eoID];

        // consider  non-ending closures on e++ 
        std::set<pair<T,unsigned int>> prevsInNonStarting;
        for(pr = 0; pr<Q.size(); pr++){
            if(Q[pr].second > 0)
                prevsInNonStarting.insert(make_pair(
                            all_closures[Q[pr].first][Q[pr].second-1],
                            e_inst[estart[Q[pr].first]+Q[pr].second-1].eoID));
            if(prevsInNonStarting.size()>1)
                break;
        }

        if(prevsInNonStarting.size()>1){
            for(const auto& p:orb)
                e_inst[estart[p.first]+p.second].type = 
                    (e_inst[estart[p.first]+p.second].type>0)? 3:2;
            return;
        }

        return;
    }
    
    template <class T, class EL>
    void SmartJoinR2L(vec<vtuple>& v_inst, const vec<etuple<EL>>&  e_inst, 
            const vec<int64_t>& estart, const int i,
            std::unordered_map<pair<T,int>,pair<int,int>>& map_ei,
            const vec<vec<T>>& all_closures,
            const vec<vec<vec<pair<int,EL>>>>& eorbits)
    {
        // get eorbits for v_inst[i]:
        int64_t ne = v_inst[i].v;
        T e = all_closures[e_inst[ne].c][e_inst[ne].m];
        const vec<pair<int,EL>>& iorb = eorbits[e][e_inst[ne].eoID];

        // get prev edge's index in e_inst and map back to v_inst, 
        // which will be joined to i and removed from verts. this joins rights to a left.
        for(const auto& io : iorb){
            if(io.second>0){
                T emm = all_closures[io.first][io.second-1];
                int cID = e_inst[estart[io.first]+io.second-1].eoID;
                pair<T,int> mp(emm,cID);
                ForceAssertEq(map_ei.count(mp),1);
                #pragma omp critical
                {GlueElements::JoinVs(v_inst, i, map_ei[mp].second);}
            }
        }
    }

    template <class T, class EL>
    void SmartJoinL2R(vec<vtuple>& v_inst, const vec<etuple<EL>>&  e_inst, 
            const vec<int64_t>& estart, const int i,
            std::unordered_map<pair<T,int>,pair<int,int>>& map_ei,
            const vec<vec<T>>& all_closures,
            const vec<vec<vec<pair<int,EL>>>>& eorbits)
    {
        // get eorbits for v_inst[i]:
        int64_t ne = v_inst[i].v;
        T e = all_closures[e_inst[ne].c][e_inst[ne].m];
        const vec<pair<int,EL>>& iorb = eorbits[e][e_inst[ne].eoID];

        // get next edge's index in e_inst and map back to v_inst, which will be 
        // joined to i and removed from verts. this joins rights to a left.
        for(const auto& io : iorb){
            if(io.second < all_closures[io.first].size()-1){
                T epp = all_closures[io.first][io.second+1];
                int cID = e_inst[estart[io.first]+io.second+1].eoID;
                pair<T,int> mp(epp,cID);
                ForceAssertEq(map_ei.count(mp),1);
                #pragma omp critical
                {GlueElements::JoinVs(v_inst, i, map_ei[mp].first);}
            }
        }
    }


    template <class T, class EL>
    int RecursiveWalk(vec<T>& unipath, const vec<etuple<EL>>& e_inst,
            const vec<int64_t>& estart, const int64_t v,
            const vec<vec<T>>& all_closures,
            std::unordered_map<pair<T,int>,pair<int,int>>& map_ei,
            const vec<vec<vec<pair<int,EL>>>>& eorbits,
            const Bool GET_META, vec<int64_t>& metadata){
        unsigned int j = 0;
        int64_t nv = v;
        while(j<all_closures[e_inst[nv].c].size()-e_inst[nv].m){
            unipath.push_back(all_closures[e_inst[nv].c][e_inst[nv].m+j]);
            if(GET_META)
                metadata.push_back(estart[e_inst[nv].c]+e_inst[nv].m+j);
            // added element is rightmost if it is type 3 or 2 in einst
            if(e_inst[estart[e_inst[nv].c]+e_inst[nv].m+j].type>1){
                // get eorbits
                T e = all_closures[e_inst[nv].c][e_inst[nv].m+j];
                int cID = e_inst[estart[e_inst[nv].c]+e_inst[nv].m+j].eoID;
                pair<T,int> mp(e,cID);
                ForceAssertEq(map_ei.count(mp),1);
                return map_ei[mp].second;
            }
            j++;
        }
        // otherwise reached end of closure, get aligning closure and redo above step
        // get singleton eorbits of (e_inst[nv].c,e_inst[nv].m+j-1)
        int64_t idx;
        T ep = all_closures[e_inst[nv].c][e_inst[nv].m+j-1];
        vec<pair<int,EL>> Q = eorbits[ep][e_inst[estart[e_inst[nv].c]+e_inst[nv].m+j-1].eoID];
        ForceAssertGe(Q.size(),2); // this confirms that there is only one aligning closure
        for(const auto& q:Q){
            if((EL)all_closures[q.first].size()==q.second+1)
                continue;
            nv = estart[q.first]+q.second+1; 
            break;
        }
        return RecursiveWalk<T>(unipath, e_inst, estart, 
                nv, all_closures, map_ei, eorbits,
                GET_META, metadata);
    }
}

template <class T, class EL>
void NucleateGraph_HP(vec<vec<om<EL>>>& omatch, T hbE, vec<vec<T>>& all_closures,
        vec<int64_t>& cinv, vec<vec<int>>& ci, digraphE<vec<T>>& D, vec<int>& dinv, 
        const Bool GET_META, vec<vec<quad<int,EL,int,EL>>>& closure_supp,
        const Bool verbose, String dir, Bool Canon)
{
    // write ci to disk and destroy
    if(dir!="")
    {
        // []
        Mkdir777(dir+"/orig");
        BinaryWriter::writeFile(dir+"/orig/a.ci",ci);
        Destroy(ci);
    // write omatches to disk and destroy
        BinaryWriter::writeFile(dir+"/orig/a.omatches",omatch);
        Destroy(omatch);
    }

    int64_t N = all_closures.size( ), total = 0;
    for ( int64_t i = 0; i < N; i++ )
      total += all_closures[i].size( );

    if (verbose){
        cout << Date( ) << ": " << ToStringAddCommas(N) << " closures "
           << "having in total " << ToStringAddCommas(total) 
           << " edges, with repetitions from " << hbE
           << " distinct edges" <<endl; }


    // define a struct used later
    if (verbose) PRINTDEETS("forming data structures for supergraph");

    vec<int64_t> estart(N+1);
    estart[0] = 0;
    for ( int i = 1; i < N+1; i++ )
        estart[i] = estart[i-1]+all_closures[i-1].size();

    vec<etuple<EL>> e_inst;
    e_inst.resize(total);

    vec<unsigned int> eorbitIds(hbE,0);
    #pragma omp parallel for
    for(int c = 0; c < N; c++){
        // generate pairs (closure containing e, index)
        for ( int m = 0; m < all_closures[c].isize(); m++ ){
            e_inst[estart[c]+m].c = c; 
            e_inst[estart[c]+m].m = m;
            e_inst[estart[c]+m].type = 0;
        }    
    }
    for(int c = 0; c < N; c++){
        for ( int m = 0; m < all_closures[c].isize(); m++ ){
            e_inst[estart[c]+m].eoID = eorbitIds[(int)all_closures[c][m]];
            eorbitIds[(int)all_closures[c][m]]++;
        }
    }
    // []
    Destroy(eorbitIds);

    // write all_closures to disk and destroy
    if(dir!=""){
        // []
        BinaryWriter::writeFile(dir+"/orig/a.all_closures",all_closures);
        Destroy(all_closures);
        // load omatches
        BinaryReader::readFile(dir+"/orig/a.omatches",&omatch);
    }

    if (verbose)
        PRINTDEETS("defining local equivalence on edges of basegraph-- caution: this can take many hours!");

    // loop through all closures/words. 
    // This way is fast since each match is looped through once.
    // However we can't make it parallelized
    // Slower way is to do this to consider the set of instances of e : (c,m) 
    // and partition them based on omatch info
    int nummatch = 0, stopped = 0, dots = 0;
    int64_t l1, l2;
    for ( int i1 = 0; i1 < omatch.isize( ); i1++ )
        nummatch += omatch[i1].size();
    for ( int i1 = 0; i1 < omatch.isize( ); i1++ ){
        for ( int j = 0; j < omatch[i1].isize( ); j++ )
        {    
            l1 = estart[i1]+omatch[i1][j].start1; 
            l2 = estart[omatch[i1][j].c2]+ omatch[i1][j].start2;
            for ( int l = 0; l < omatch[i1][j].len; l++ )
                GlueElements::JoinEs<EL>(e_inst,estart, l1 + l, l2 + l);
            if(verbose)
                MakeDots(stopped,dots,nummatch);
        }    
    }
 
    // []
    Destroy(omatch);

    if(verbose) PRINTQUANTS("number of mergers done = ",ToStringAddCommas(nummatch));
    
    if(dir!=""){
        BinaryReader::readFile(dir+"/orig/a.all_closures",&all_closures);
        BinaryReader::readFile(dir+"/orig/a.ci",&ci);
    }

    if (verbose) PRINTDEETS("precompute edge orbits and canonicalize classIds");
    // HIGH MEM
    vec<vec<vec<pair<int,EL>>>> eorbits(hbE);
    #pragma omp parallel for schedule( dynamic, 1 )
    for(T e = 0; e < hbE; e++){
        // get partitions/orbits of (c,m) wrt e
        std::unordered_map<unsigned int,unsigned int> IDpos;
        vec<vec<pair<int,EL>>>& orbits = eorbits[e];
        unsigned int pos =0;
        for(int j = 0; j< ci[e].isize(); j++){
            int c = ci[e][j];
            const vec<T>& x = all_closures[c];
            for(int m = 0; m<x.isize(); m++){
                if (x[m] == e){
                    unsigned int thiseoID = e_inst[estart[c]+m].eoID;
                    if(IDpos.count(thiseoID)==0){
                        IDpos[thiseoID] = pos++;
                        orbits.push_back(vec<pair<int,EL>>());
                        orbits.back().push(c,m);
                    }else 
                        orbits[IDpos[thiseoID]].push(c,m); 
                    
                }
            }
        }
    }
    for(T e = 0; e < hbE; e++){
        //for each orbit
        auto& orbits = eorbits[e];
        #pragma omp parallel for
        for(unsigned int k = 0; k < orbits.size(); k++){
            for(auto q: orbits[k])
                e_inst[estart[q.first]+q.second].eoID = k;
        }
    }

    if (verbose) PRINTDEETS("finding edge type & marking edges to/from vertices");

    // populate type 
    // 0 if ordinary
    // 1 if leftmost-only
    // 2 if rightmost-only
    // 3 if both left and rightmost
    #pragma omp parallel for schedule( dynamic, 1 )
    for(T e = 0; e < hbE; e++){
        // get partitions/orbits of (c,m) wrt e
        auto& orbits = eorbits[e];
        //for each orbit
        for(const auto& orb : orbits){
            BuildSuperGraph_HP::IsLeftMost<T,EL>(e_inst,estart,orb,all_closures,eorbits);
        } 
    }
    #pragma omp parallel for schedule( dynamic, 1 )
    for(T e = 0; e < hbE; e++){
        // get partitions/orbits of (c,m) wrt e
        auto& orbits = eorbits[e];
        //for each orbit
        for(const auto& orb : orbits){
            BuildSuperGraph_HP::IsRightMost<T,EL>(e_inst,estart,orb,all_closures,eorbits);
        } 
    }
    if(verbose) PRINTDEETS("defining global equivalence on vertices of supergraph");

    // Construct a dense index of vtuple(id in e_inst, cycle vertex, void, type) to e_inst 
    // for all etuple with type > 0. Add type==1 and type==2; type==3 is added once each as type 1 and 2. 
    vec<vtuple> v_inst;
    v_inst.reserve(2*hbE); // guessing the size
    std::unordered_map<pair<T,int>,pair<int,int>> map_ei; // build inverse index, can be a vec
    map_ei.reserve(2*hbE);

    if(Canon){
        int vorbitIds = 0;
        vec<int64_t> select_inst;
        #pragma omp parallel for schedule( dynamic, 1 )
        for(T e = 0; e < hbE; e++){
            // get partitions/orbits of (c,m) wrt e
            auto& orbits = eorbits[e];
            //for each orbit
            for(const auto& orb : orbits){
               int64_t ii = estart[orb[0].first]+orb[0].second;
               if(e_inst[ii].type > 0){
                  #pragma omp critical 
                  { select_inst.push_back(ii); }
                   for(const auto& q: orb)
                       e_inst[estart[q.first]+q.second].type *= -1;
               }
           }
        }
        // reset edge type
        for(auto& e : e_inst)
            e.type *= -1;

        // canonicalize stuff!!
        ParallelSort(select_inst);
        for(auto ii:select_inst){
            T ed = all_closures[e_inst[ii].c][e_inst[ii].m];
            int cID = e_inst[ii].eoID;
            if(e_inst[ii].type==3){
                map_ei[make_pair(ed,cID)] = make_pair(vorbitIds, vorbitIds+1);
                v_inst.push(ii, vorbitIds ,vorbitIds,1);
                v_inst.push(ii, vorbitIds+1 ,vorbitIds+1,2);
                vorbitIds+=2;
            }
            if(e_inst[ii].type==1){
                map_ei[make_pair(ed,cID)]=make_pair(vorbitIds,-1);
                v_inst.push(ii, vorbitIds, vorbitIds, e_inst[ii].type);
                vorbitIds++;
            }
            if(e_inst[ii].type==2){
                map_ei[make_pair(ed,cID)]=make_pair(-1,vorbitIds);
                v_inst.push(ii, vorbitIds, vorbitIds, e_inst[ii].type);
                vorbitIds++;
            }
        }
        // []
        Destroy(select_inst);
    }
    else{
        int vorbitIds = 0;
        #pragma omp parallel for ordered schedule( dynamic, 1 )
        for(T e = 0; e < hbE; e++){
            // get partitions/orbits of (c,m) wrt e
            auto& orbits = eorbits[e]; 
            //for each orbit
            #pragma omp ordered
            {  for(const auto& orb : orbits){
                   int64_t ii = estart[orb[0].first]+orb[0].second;
                   T ed = all_closures[orb[0].first][orb[0].second];
                   int cID = e_inst[ii].eoID;
                   int turn = 0;
                   if(e_inst[ii].type == 3){
                       map_ei[make_pair(ed,cID)]= make_pair(vorbitIds,vorbitIds+1);
                       v_inst.push(ii, vorbitIds ,vorbitIds,1);
                       v_inst.push(ii, vorbitIds+1 ,vorbitIds+1,2);
                       vorbitIds+=2;
                       turn = 1;
                   }
                   if(e_inst[ii].type == 1){
                       map_ei[make_pair(ed,cID)]=make_pair(vorbitIds,-1);
                       v_inst.push(ii, vorbitIds, vorbitIds, e_inst[ii].type);
                       vorbitIds++;
                       turn = 1;
                   }
                   if(e_inst[ii].type == 2){
                       map_ei[make_pair(ed,cID)]=make_pair(-1,vorbitIds);
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
    }
    // []
    Destroy(ci);

    if (verbose) PRINTDEETS("creating representatives for verts");

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
        PRINTQUANTS("#E in supergraph",ToStringAddCommas(v_inst.size()/2));
        PRINTDEETS("completing vertex aliasing --caution: this can take many hours!");
    }

    // v1 equiv v2 if they form a successive left-right pair along a common closure
    stopped = 0, dots = 0;
    #pragma omp parallel for schedule(dynamic,10)
    for(int i = 0; i< repL.isize(); i++){
        BuildSuperGraph_HP::SmartJoinR2L<T,EL>(v_inst, e_inst, estart, 
                repL[i], map_ei, all_closures, eorbits);
        if(verbose){
            #pragma omp critical
            { MakeDots(stopped,dots,repL.isize());}
        }
    }

    if (verbose) PRINTDEETS("remapping vertex orbitIDs supergraph");

    // make the vertex orbit IDs to be consecutive, starting from 0
    std::unordered_map<int,int> vmap;
    int ix = 0;
    for(const auto& v : v_inst){
    if(vmap.count(v.voID)==0)
        vmap[v.voID] = ix++;
    }
      
    if (verbose)
    {   
        PRINTQUANTS("#V in supergraph",ToStringAddCommas(vmap.size()));
        PRINTDEETS("constructing edge objects in supergraph");
    }

    // construct edge_obj 
    vec<vec<T>> unipaths(repL.size());
    vec<int> permRepR(repR.size(),-1);
    std::unordered_map<vec<int64_t>,int> Start;
    Start.reserve(unipaths.size());
    vec<vec<int64_t>> metadata;
    metadata.resize(repL.size());
    // with each leftmost edge
    #pragma omp parallel for
    for(int iv = 0; iv < repL.isize(); iv++){
        permRepR[iv] = BinPosition(repR,BuildSuperGraph_HP::RecursiveWalk<T,EL>(
                    unipaths[iv],e_inst,estart,v_inst[repL[iv]].v,
                    all_closures, map_ei, eorbits, GET_META, metadata[iv]));
        // get eorbits of v_inst[repL[iv]].v as (c,m)
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

    // []
    Destroy(eorbits);

    int64_t SUM = 0;
    for(const auto& upath: unipaths)
        SUM += upath.size();

    if (verbose)
    {   
        PRINTQUANTS("#base edges in supergraph",ToStringAddCommas(SUM));
        PRINTDEETS("constructing involution on the supergraph");
    }

    vec<int> ws(N,0);
    for ( int64_t i = 0; i < N; i++ )
        ws[i] = all_closures[i].size( );
    // []
    Destroy(all_closures); 

    dinv.resize(unipaths.size(),-1);
    #pragma omp parallel for
    for(int p = 0; p < unipaths.isize(); p++){
        // get eorbits of v_inst[repR[permRepR[p]]].v as (c,m)
        vec<pair<int,EL>> Q;
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

    if(GET_META){
        if (verbose) PRINTDEETS("get meta data; WARNING: SLOW & HIGH MEM");
        
        // for each unipath metadata
        closure_supp.clear();
        closure_supp.resize(metadata.size());
        #pragma omp parallel for
        for(int64_t kk = 0; kk < metadata.size(); kk++){
            vec<triple<int,EL,int>> match;
            match.reserve(10*metadata[kk].size()); //heuristic
            for(int idx = 0; idx<  metadata[kk].isize(); idx++){
                int64_t nv = metadata[kk][idx];
                // get nv orb: {c,m}
                int64_t ne = estart[e_inst[nv].c]+e_inst[nv].m;
                int64_t st = ne;
                while(1){
                    match.push(e_inst[ne].c,e_inst[ne].m,idx);
                    ne = estart[e_inst[ne].c]+e_inst[ne].m; 
                    if(ne==st) break;
                }
            }
            Sort(match);
            // mark boundaries and populate paths in sorted manner
            int cl = -1;
            for(int i = 0; i < match.isize(); i++){ 
                if(i>0){
                    // if another segment matches, new start!
                    if( (cl==match[i].first) &&
                        ((match[i].second!=match[i-1].second+1) || 
                        (match[i].third!=match[i-1].third+1)) )
                        closure_supp[kk].push_back(quad<int,EL,int,EL>(
                            match[i].first, match[i].second, match[i].third,0));
                }
                if(match[i].first!=cl){
                    cl = match[i].first;
                    closure_supp[kk].push_back(quad<int,EL,int,EL>(
                        match[i].first, match[i].second, match[i].third,0));
                }
                closure_supp[kk].back().fourth++;
            }
        }

        if(dir!="")
            BinaryWriter::writeFile(dir+"/orig/a.clpaths",closure_supp);
    }

    // []
    Destroy(metadata);

    // []
    Destroy(ws);
    Destroy(e_inst);
    Destroy(estart);

    if (verbose) PRINTDEETS("constructing supergraph");

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
    if(verbose) PRINTDEETS("sortsyncing");

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

    // []
    Destroy(v_inst);
    Destroy(repR);
    Destroy(repL);
    Destroy(permRepR);
    vmap.clear();
}

// for base edge sequence
template void NucleateGraph_HP(vec<vec<om<uint16_t>>>& omatch, uchar hbE, 
        vec<vec<uchar>>& all_closures, vec<int64_t>& cinv, vec<vec<int>>& ci, 
        digraphE<vec<uchar>>& D, vec<int>& dinv, 
        const Bool GET_META, 
        vec<vec<quad<int,uint16_t,int,uint16_t>>>& closure_supp,
        const Bool verbose, String dir, Bool Canon);

// standard case
template void NucleateGraph_HP(vec<vec<om<uint16_t>>>& omatch, int hbE, 
        vec<vec<int>>& all_closures, vec<int64_t>& cinv, vec<vec<int>>& ci, 
        digraphE<vec<int>>& D, vec<int>& dinv, 
        const Bool GET_META, 
        vec<vec<quad<int,uint16_t,int,uint16_t>>>& closure_supp,
        const Bool verbose, String dir, Bool Canon);
