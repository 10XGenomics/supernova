// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>

#include "10X/mergers/GlueGraphs.h"
#include "10X/mergers/NicePrints.h"

namespace ExtractGraphs{

    template <class T,class EL>
    void IsLeftMost(vec<etuple<EL>>& e_inst, const vec<int64_t>& estart, 
            const vec<pair<int,EL>>& orb, const digraphE<vec<T>>& D,
            const vec<int>& to_lefts, const vec<int>& to_rights, T legitE){ 

        if(All_Closures::first<T>(D,orb[0].first)>=legitE && 
                orb[0].second==0){
            for(const auto& p:orb)
                e_inst[estart[p.first]+p.second].type = 
                    (e_inst[estart[p.first]+p.second].type>0)? 3:1;
            return;
        }

        // cycle through given orbit{c.m} of e 
        int pr = -1;
        std::set<pair<T,unsigned int>> prevs;
        vec<int> leftds;
        for (unsigned int p = 0; p<orb.size(); p++){
            if(orb[p].second>0){
                prevs.insert(make_pair(
                            All_Closures::elem<T>(D,orb[p].first,orb[p].second-1),
                            e_inst[estart[orb[p].first]+orb[p].second-1].eoID));
                pr = p;
            }else{
                // get prevs along connecting edges to left, if present 
                vec<int> ds; LeftEdges<T>(D,orb[p].first,to_lefts,ds);
                for(auto d : ds){
                    prevs.insert(make_pair(All_Closures::last<T>(D,d),
                                e_inst[estart[d+1]-1].eoID));
                    leftds.push_back(d);
                }
            }
        }

        // if all starting OR multiple distinct e-- instances, then LM
        // if e is start of a gap edge! then LM
        if(prevs.size()!=1){
            for(const auto& p:orb)
                e_inst[estart[p.first]+p.second].type = 
                    (e_inst[estart[p.first]+p.second].type>0)? 3:1;
            return;
        }

               
        if((*prevs.begin()).first>=legitE && orb[0].second==0){
            for(const auto& p:orb)
                e_inst[estart[p.first]+p.second].type = 
                    (e_inst[estart[p.first]+p.second].type>0)? 3:1;
            return;
        }
 
        // get orbit of e-- containing (closure,instance) connecting with special e
        int64_t idx;
        vec<pair<int,EL>> Q;
        if(pr>-1)
            Q.push(orb[pr].first,orb[pr].second-1);
        else
            Q.push(leftds[0],All_Closures::size(estart,leftds[0])-1); 
        int64_t start = estart[Q.back().first]+Q.back().second;
        while(1){
            idx = estart[Q.back().first]+Q.back().second;
            if(estart[e_inst[idx].c]+e_inst[idx].m == start)
                break;
            Q.push(e_inst[idx].c,e_inst[idx].m);
        }

        // consider non-ending closures on e-- and get set of nexts
        // consider connecting closures to e-- and get set of nexts
        std::set<pair<T,unsigned int>> nextsOfprevs;
        for(pr = 0; pr<Q.isize(); pr++){
            if(Q[pr].second < All_Closures::isize(estart,Q[pr].first)-1){
                nextsOfprevs.insert(make_pair(
                            All_Closures::elem<T>(D,Q[pr].first,Q[pr].second+1),
                            e_inst[estart[Q[pr].first]+Q[pr].second+1].eoID));
            }else{
                vec<int> ds; RightEdges<T>(D,Q[pr].first,to_rights,ds);
                for(auto d: ds)
                    nextsOfprevs.insert(make_pair(
                                All_Closures::first<T>(D,d),
                                e_inst[estart[d]].eoID));
            }
            if(nextsOfprevs.size()>1)
                break;
        }

        if(nextsOfprevs.size()>1){
            for(const auto& p:orb)
                e_inst[estart[p.first]+p.second].type = 
                    (e_inst[estart[p.first]+p.second].type>0)? 3:1;
            return;
        }
        return;
    }

    template <class T, class EL>
    void IsRightMost(vec<etuple<EL>>& e_inst, const vec<int64_t>& estart, 
            const vec<pair<int,EL>>& orb, const digraphE<vec<T>>& D,
            const vec<int>& to_lefts, const vec<int>& to_rights, T legitE){

        if(All_Closures::last<T>(D,orb[0].first)>=legitE && 
                orb[0].second==All_Closures::isize(estart,orb[0].first)-1){
            for(const auto& p:orb)
                e_inst[estart[p.first]+p.second].type = 
                    (e_inst[estart[p.first]+p.second].type>0)? 3:2;
            return;
        }

        int pr = -1;
        std::set<pair<T,unsigned int>> nexts;
        vec<int> rightds;
        for (unsigned int p = 0; p<orb.size(); p++){
               
            if(orb[p].second < All_Closures::isize(estart,orb[p].first)-1){
                nexts.insert(make_pair(
                            All_Closures::elem<T>(D,orb[p].first,orb[p].second+1),
                            e_inst[estart[orb[p].first]+orb[p].second+1].eoID));
                pr = p;
            }else{
                // get nexts along connecting edges to right, if present 
                vec<int> ds; RightEdges<T>(D,orb[p].first,to_rights,ds);
                for(auto d : ds){
                    nexts.insert(make_pair(All_Closures::first<T>(D,d),
                                e_inst[estart[d]].eoID));
                    rightds.push_back(d);
                }
            }
        }

        // if ending, or lost of nexts or gap
        if(nexts.size()!=1){
            for(const auto& p:orb)
                e_inst[estart[p.first]+p.second].type = 
                    (e_inst[estart[p.first]+p.second].type>0)? 3:2;
            return;
        }  

        if((*nexts.begin()).first>=legitE && 
                orb[0].second==All_Closures::isize(estart,orb[0].first)-1){
            for(const auto& p:orb)
                e_inst[estart[p.first]+p.second].type = 
                    (e_inst[estart[p.first]+p.second].type>0)? 3:2;
            return;
        }

        // get orbit of e++ containing (closure,instance) connecting with special e
        int64_t idx;
        vec<pair<int,EL>> Q;
        if(pr>-1)
            Q.push(orb[pr].first,orb[pr].second+1);
        else
            Q.push(rightds[0],0);
        int64_t start = estart[Q.back().first]+Q.back().second;
        while(1){
            idx = estart[Q.back().first]+Q.back().second;
            if(estart[e_inst[idx].c]+e_inst[idx].m == start)
                break;
            Q.push(e_inst[idx].c,e_inst[idx].m);
        }

        // consider  non-ending closures on e++ 
        std::set<pair<T,unsigned int>> prevsOfnexts;
        for(pr = 0; pr<Q.isize(); pr++){
            if(Q[pr].second > 0){
                prevsOfnexts.insert(make_pair(
                            All_Closures::elem<T>(D,Q[pr].first,Q[pr].second-1),
                            e_inst[estart[Q[pr].first]+Q[pr].second-1].eoID));
            }else{
                vec<int> ds; LeftEdges<T>(D,Q[pr].first,to_lefts,ds);
                for(auto d: ds)
                    prevsOfnexts.insert(make_pair(
                                All_Closures::last<T>(D,d),
                                e_inst[estart[d+1]-1].eoID));

            }
            if(prevsOfnexts.size()>1)
                break;
        }

        if(prevsOfnexts.size()>1){
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
            const digraphE<vec<T>>& D,
            std::unordered_map<pair<T,int>,pair<int,int>>& map_ei,
            const vec<int>& to_lefts)
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

        // get prev edge's index in e_inst either along closure or via connection
        // and map back to v_inst to get vertex id j
        // which will be joined to i. this joins rights to a left.
        for(const auto& io : iorb){
            if(io.second>estart[io.first]){
                // find orbit of io.second-1
                int64_t pe = io.second-1;
                T emm = All_Closures::elem(D,e_inst[pe].c,e_inst[pe].m);
                int cID = e_inst[pe].eoID;
                pair<T,int> mp(emm,cID);
                #pragma omp critical
                {GlueElements::JoinVs(v_inst, i, map_ei[mp].second);}
            }else{
                // get left ds
                vec<int> ds; LeftEdges<T>(D,io.first,to_lefts,ds);
                // get ends of ds and their orbits
                // map ends of ds to vs and join with i
                for(auto d : ds){
                    int64_t pe = estart[d+1]-1;
                    T emm = All_Closures::elem(D,e_inst[pe].c,e_inst[pe].m);
                    int cID = e_inst[pe].eoID;
                    pair<T,int> mp(emm,cID);
                    #pragma omp critical
                    {GlueElements::JoinVs(v_inst, i, map_ei[mp].second);}
                }
            }
        }
    }

    template <class T,class EL>
    void SmartJoinL2R(vec<vtuple>& v_inst, const vec<etuple<EL>>&  e_inst, 
            const vec<int64_t>& estart, const int i,
            const digraphE<vec<T>>& D,
            std::unordered_map<pair<T,int>,pair<int,int>>& map_ei,
            const vec<int>& to_rights)
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
                T epp = All_Closures::elem(D,e_inst[pe].c,e_inst[pe].m);
                int cID = e_inst[pe].eoID;
                pair<T,int> mp(epp,cID);
                #pragma omp critical
                {GlueElements::JoinVs(v_inst, i, map_ei[mp].first);}

            }else{
                // get right ds
                vec<int> ds; RightEdges<T>(D,io.first,to_rights,ds);
                // get ends of ds and their orbits
                // map ends of ds to vs and join with i
                for(auto d : ds){
                    int64_t pe = estart[d];
                    T epp = All_Closures::elem(D,e_inst[pe].c,e_inst[pe].m);
                    int cID = e_inst[pe].eoID;
                    pair<T,int> mp(epp,cID);
                    #pragma omp critical
                    {GlueElements::JoinVs(v_inst, i, map_ei[mp].first);}
                }
            }
        }
    }

    template <class T, class EL>
    int Walk(vec<T>& unipath, const vec<etuple<EL>>& e_inst,
            const vec<int64_t>& estart, int64_t& v,
            const digraphE<vec<T>>& D,
            std::unordered_map<pair<T,int>,pair<int,int>>& map_ei, 
            const vec<int>& to_rights, const Bool GET_META,
            vec<int64_t>& metadata){
        unsigned int j = 0;
        int64_t nv = v;
        while(j<All_Closures::size(estart,e_inst[nv].c)-e_inst[nv].m)
        {
            unipath.push_back(All_Closures::elem<T>(D,e_inst[nv].c,e_inst[nv].m+j));
            if(GET_META)
                metadata.push_back(estart[e_inst[nv].c]+e_inst[nv].m+j);
            // added element is rightmost if it is type 3 or 2 in einst
            if(e_inst[estart[e_inst[nv].c]+e_inst[nv].m+j].type>1)
            {
                // get eorbit
                int64_t ne = estart[e_inst[nv].c]+e_inst[nv].m+j;
                T ed = All_Closures::elem(D,e_inst[ne].c,e_inst[ne].m);
                int cID = e_inst[ne].eoID;
                pair<T,int> mp(ed,cID);
                ForceAssertEq(map_ei.count(mp),1);
                return map_ei[mp].second;
            }
            j++;
        }
        // otherwise reached end of closure, get aligning closure and redo above step
        // get singleton eorbit of (e_inst[nv].c,e_inst[nv].m+j-1)
        int64_t idx;
        vec<pair<int,EL>> Q;
        Q.push(e_inst[nv].c,e_inst[nv].m+j-1);
        int64_t start = estart[Q.back().first]+Q.back().second;
        while(1){
            idx = estart[Q.back().first]+Q.back().second;
            if(estart[e_inst[idx].c]+e_inst[idx].m == start)
                break;
            Q.push(e_inst[idx].c,e_inst[idx].m);
        }

        // get non ending aligned closure
        // if all end, get a right
        // otherwise we got a problem
        vec<int> ds;
        Bool extends = False;
        for(const auto& q:Q){ 
            if((EL)All_Closures::isize(estart,q.first)==(q.second+1)){
                // get rights
                vec<int> dss;
                RightEdges<T>(D,q.first,to_rights,dss);
                ds.insert(ds.end(), dss.begin(), dss.end());
                continue;
            }
            nv = estart[q.first]+q.second+1; 
            extends = True;
            break;
        }
        if(!extends){
            // check logic
            if(Q.size()>1)
                ForceAssertGe(ds.size(),1);
            // This check fails in zippering
            // as one could have multiple edges emanating from
            // one end of closures, all waiting to be zipped!
            /* else */
            /*     ForceAssertEq(ds.size(),1); */ 
            // just pick one
            nv = estart[ds[0]];
        }
        v = nv;
        return -2;
    }

    template <class T, class EL>
    int RecursiveWalk(vec<T>& unipath, const vec<etuple<EL>>& e_inst,
            const vec<int64_t>& estart, const int64_t v,
            const digraphE<vec<T>>& D,
            std::unordered_map<pair<T,int>,pair<int,int>>& map_ei, 
            const vec<int>& to_rights, const Bool GET_META,
            vec<int64_t>& metadata){

        int64_t nv = v;
        int k = -1;
        while(k<0){
            k = Walk<T,EL>(unipath, e_inst, estart, nv, 
                    D, map_ei, to_rights, GET_META, metadata);
        }
        return k;
    }
}

template <class T, class EL>
void GlueGraphs(digraphE<vec<T>>& D, vec<int>& dinv, vec<nptuple>& Matches, 
        T hbE, digraphE<vec<T>>& D_new, vec<int>& dinv_new, 
        const Bool GET_META, vec<vec<quad<int,EL,int,EL>>>& merger_supp,
        const Bool verbose, String dir)
{
    if ( verbose ) PRINTDEETS("Gluing Graphs");

    // ------------- GATHER SOME STATS -------------------
    int64_t total = 0;
    int max_e_sz = 0;
    int64_t N = D.E();
    for(int d = 0; d<D.E(); d++){
        total += D.O(d).size();
        max_e_sz = Max(max_e_sz, D.O(d).isize());
    }
    if (verbose){
        PRINTQUANTS("total superedges = ",ToStringAddCommas(N));
        PRINTQUANTS("total baseedges = ",ToStringAddCommas(total));
        PRINTQUANTS("base assembly size = ",ToStringAddCommas((int)hbE));
        PRINTQUANTS("max size of superedges = ",ToStringAddCommas(max_e_sz));
    }

    // ------------ REPLACE GAPS -------------------------

    std::unordered_map<T,vec<T>> Gaps;
    T nE = hbE;
    int numgaps = 0; int64_t totalgapedges = 0; int nonzeros = 0;
    for(int c =0; c<D.E(); c++){
        if(D.O(c).size()==0) continue;
        nonzeros++;
        if(D.O(c)[0]<0){
            numgaps++;
            totalgapedges += D.O(c).size();
            Gaps.insert(make_pair(nE,D.O(c)));
            // modifiy gaps seq
            for(int m = 0; m< D.O(c).isize(); m++){
                D.OMutable(c)[m] = nE; nE++;
            }
        }
    }

    if (verbose){
        PRINTQUANTS("number of gap superedges = ",ToStringAddCommas(numgaps));
        PRINTQUANTS("total gap base edges = ",ToStringAddCommas(totalgapedges));
        PRINTQUANTS("number of nonzero superedges = ",ToStringAddCommas(nonzeros));
    }

    if(verbose) PRINTDEETS("replace gaps with placeholder seq");

    // ------------ INDEX CLOSURES ------------------------
    if (verbose) PRINTDEETS("indexing closures + uniquesorting");

    vec<vec<int >> ci;
    InvertClosures<T>( D,ci,nE);

    if (verbose) 
        PRINTDEETS(" added fake "<<ToStringAddCommas(nE-hbE)
                   <<" base edges corresponding to gaps");

    // []
    if(dir!=""){
        BinaryWriter::writeFile(dir+"/a.ci",ci);
        Destroy(ci);
    } 
 
    // ---------------- EDGE TYPE INHERITANCE -----------------------------
    
    if (verbose) PRINTDEETS("init structs");

    vec<int64_t> estart(N+1);
    estart[0] = 0;
    for ( int i = 1; i < N+1; i++ )
        estart[i] = estart[i-1]+All_Closures::isize<T>(D,i-1);
    ForceAssertEq(estart.back(),total);

    vec<etuple<EL>> e_inst;
    e_inst.resize(total);

    vec<unsigned int> eorbitIds(nE,0);
    #pragma omp parallel for
    for(int c = 0; c < N; c++){
        for ( int m = 0; m < All_Closures::isize(estart,c); m++ ){
            e_inst[estart[c]+m].c = c; 
            e_inst[estart[c]+m].m = m;
            e_inst[estart[c]+m].type = 0;
        }
    }
    for(int c = 0; c < N; c++){
        for ( int m = 0; m < All_Closures::isize(estart,c); m++ ){
            e_inst[estart[c]+m].eoID = eorbitIds[All_Closures::elem<T>(D,c,m)];
            eorbitIds[All_Closures::elem<T>(D,c,m)]++;
        }
    }
 
    // []
    Destroy(eorbitIds);

    // ----------------- EDGE ALIASING ----------------------------
    if (verbose) 
        PRINTDEETS("join orbits of base edge instance-- caution: this can take many hours!");

    // Can't parallelize this at all.
    // Make sure matches are symmetric wrt involution!!!
    int64_t nummerges = 0 , stopped = 0, dots = 0;
    for(int64_t i = 0; i < (int64_t) Matches.size(); i++){
        int64_t l1 =0, l2=0, l1i =0, l2i =0;
        l1 = estart[Matches[i].first.first] + Matches[i].first.second; 
        l2 = estart[Matches[i].second.first] + Matches[i].second.second; 
        /* l1i = estart[dinv[Matches[i].first.first]] + */ 
        /*     D.O(dinv[Matches[i].first.first]).size()-1-Matches[i].first.second; */
        /* l2i = estart[dinv[Matches[i].second.first]] + */ 
        /*     D.O(dinv[Matches[i].second.first]).size()-1-Matches[i].second.second; */
        for(int l = 0; l<Matches[i].third; l++){
            if(l==0)
                nummerges++;
            GlueElements::JoinEs<EL>(e_inst,estart, l1+l, l2+l);
            /* GlueElements::JoinEs<EL>(e_inst,estart, l1i-l, l2i-l); */
        }
        if(verbose)
            MakeDots(stopped,dots,(int64_t) Matches.size());
    }
 
    if (verbose) 
        PRINTDEETS("performed "<<ToStringAddCommas(nummerges)<<" mergers out of "
                   <<ToStringAddCommas(Matches.size())<<" candidates");

    // []
    Destroy(Matches);
    
    Bool diagnose = False;

    if(diagnose)
        cout << "Checksums #1: "<<checksum<EL>(e_inst) <<endl;

    // ---------------- INSERT DIAGNOSTICS1 -----------------
    // ------------------------------------------------------
    
    // --------------- IDENTIFY LEFTMOST AND RIGHTMOST ----------------
    if (verbose) PRINTDEETS("finding edge type & marking edges to/from vertices");

    if(dir!=""){
        BinaryReader::readFile(dir+"/a.ci",&ci);
    }
 
    vec<int> to_lefts, to_rights;
    D.ToLeftParallel(to_lefts);
    D.ToRightParallel(to_rights);

    // populate type 
    // 0 if ordinary
    // 1 if leftmost-only
    // 2 if rightmost-only
    // 3 if both left and rightmost
    #pragma omp parallel for schedule( dynamic, 1000 )
    for(T e = 0; e < nE; e++){
        // get partitions/orbits of (c,m) wrt e
        vec<unsigned int> eoID;
        vec<vec<pair<int,EL>>> orbits;
        for(int j = 0; j< ci[e].isize(); j++){
            int c = ci[e][j];
            const vec<T> x = All_Closures::obj<T>(D,c);
            for(int m = 0; m<x.isize(); m++){
                if (x[m] == e){
                    unsigned int thiseoID = e_inst[estart[c]+m].eoID;
                    unsigned int pos = std::find(eoID.begin(),eoID.end(),thiseoID)-eoID.begin();
                    // if not already taken, add to eoID and orbits
                    if( pos == eoID.size()){
                        orbits.push_back(vec<pair<int,EL>>());
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
            ExtractGraphs::IsLeftMost<T,EL>(e_inst,estart,orb,D,to_lefts,to_rights,hbE);
        } 
    }
    #pragma omp parallel for schedule( dynamic, 1000 )
    for(T e = 0; e < nE; e++){
        // get partitions/orbits of (c,m) wrt e
        vec<unsigned int> eoID;
        vec<vec<pair<int,EL>>> orbits;
        for(int j = 0; j< ci[e].isize(); j++){
            int c = ci[e][j];
            const vec<T> x = All_Closures::obj<T>(D,c);
            for(int m = 0; m<x.isize(); m++){
                if (x[m] == e){
                    unsigned int thiseoID = e_inst[estart[c]+m].eoID;
                    unsigned int pos = std::find(eoID.begin(),eoID.end(),thiseoID)-eoID.begin();
                    // if not already taken, add to eoID and orbits
                    if( pos == eoID.size()){
                        orbits.push_back(vec<pair<int,EL>>());
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
            ExtractGraphs::IsRightMost<T,EL>(e_inst,estart,orb,D,to_lefts,to_rights,hbE);
        } 
    }

    if(diagnose)
        cout << "Checksums #2: "<<checksum<EL>(e_inst) <<endl;

    // ---------------- INSERT DIAGNOSTICS2 -----------------
    // ------------------------------------------------------

    // --------------------- STORE VERTS -----------------------
    if (verbose) PRINTDEETS("joining orbits of edges to define vertices");

    // Construct a dense index of vtuple(id in e_inst, cycle vertex, void, type) to e_inst 
    // for all etuple with type > 0. Add type==1 and type==2; type==3 is 
    // added once each as type 1 and 2. 
    vec<vtuple> v_inst;
    v_inst.reserve(2*nE); // guessing the size
    std::unordered_map<pair<T,int>,pair<int,int>> map_ei; // build inverse index, can be a vec
    map_ei.reserve(2*nE);
    int vorbitIds = 0;
    vec<int64_t> select_inst;
    #pragma omp parallel for schedule( dynamic, 1000 )
    for(T e = 0; e < nE; e++){ 
        // get partitions/orbits of (c,m) wrt e
        vec<unsigned int> eoID;
        vec<vec<pair<int,EL>>> orbits;
        for(int j = 0; j< ci[e].isize(); j++){
            int c = ci[e][j];
            const vec<T> x = All_Closures::obj<T>(D,c);
            for(int m = 0; m<x.isize(); m++){
                if (x[m] == e){
                    unsigned int thiseoID = e_inst[estart[c]+m].eoID;
                    unsigned int pos = std::find(eoID.begin(),eoID.end(),thiseoID)-eoID.begin();
                    // if not already taken, add to eoID and orbits
                    if( pos == eoID.size()){
                        orbits.push_back(vec<pair<int,EL>>());
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

    // Canonicalize stuff!!
    ParallelSort(select_inst);
    for(auto ii:select_inst){
        T ed = All_Closures::elem(D,e_inst[ii].c,e_inst[ii].m);
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
    Destroy(ci);
    if(diagnose)
        cout << "Checksums #3: "<<checksum(v_inst) <<endl;

    // -------------------- GET VERTS REPS ---------------
    if (verbose) PRINTDEETS("getting representatives for vertices");

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
    ForceAssertEq(repL.size(),repR.size());
    if(diagnose)
        cout << "Checksums #4: "<<checksum(repR)<<","<<checksum(repL) <<endl;

    // -------------------- SMART JOIN L & R --------------
    if (verbose){
        PRINTQUANTS("Final graph: #E = ", ToStringAddCommas(v_inst.size()/2));
        PRINTDEETS("finish merging vertex orbits-- caution: this can take many hours!");
    }

    // v1 equiv v2 if they form a successive left-right pair along a common closure
    stopped = 0, dots = 0;
    #pragma omp parallel for schedule(dynamic,100)
    for(int i = 0; i< repL.isize(); i++){
        ExtractGraphs::SmartJoinR2L<T,EL>(v_inst,e_inst,estart,repL[i],D,map_ei,to_lefts);
        if(verbose){
            #pragma omp critical
            { MakeDots(stopped,dots,(int64_t) repL.size());}
        }
    }

    //[]
    Destroy(to_lefts);
    if(diagnose)
        cout << "Checksums #5: "<<checksum(v_inst) <<endl;

    // ----------------- RENUMBERING VERTS ----------------
    if (verbose) PRINTDEETS("remapping vertex orbitIDs in supergraph");

    // make the vertex orbit IDs to be consecutive, starting from 0
    std::unordered_map<int,int> vmap;
    int ix = 0;
    for(const auto& v : v_inst){
    if(vmap.count(v.voID)==0)
        vmap[v.voID] = ix++;
    }
     
    // ---------------- RECURSIVE WALK ------------
    if (verbose)
    {   
        PRINTQUANTS("Final graph : #V= ",ToStringAddCommas(vmap.size()));
        PRINTDEETS("constructing edges in supergraph");
    }

    // construct edge_obj 
    vec<vec<T>> unipaths(repL.size());
    vec<int> permRepR(repR.size(),-1);
    std::unordered_map<vec<int64_t>,int> Start;
    Start.reserve(unipaths.size());
    vec<vec<int64_t>> metadata;
    metadata.resize(unipaths.size());
    // with each leftmost edge
    #pragma omp parallel for
    for(int iv = 0; iv < repL.isize(); iv++){
        permRepR[iv] = BinPosition(repR,ExtractGraphs::RecursiveWalk<T,EL>(
                    unipaths[iv], e_inst,estart,v_inst[repL[iv]].v, 
                    D, map_ei, to_rights, GET_META, metadata[iv]));
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
    //[]
    Destroy(to_rights);

    int64_t SUM = 0;
    for(const auto& upath: unipaths)
        SUM += upath.size();

    // --------------- GET INVOLUTION ----------------
    if (verbose)
    {  
        PRINTQUANTS("Final graph: total base edges = ",ToStringAddCommas(SUM));
        PRINTDEETS("constructing involution on the supergraph");
    }

    vec<int> ws(N,0);
    for ( int64_t c = 0; c < N; c++ )
        ws[c] = All_Closures::isize<T>(D,c);

    // []
    /* destroy(D); */
    dinv_new.resize(unipaths.size(),-1);
    #pragma omp parallel for
    for(int p = 0; p < unipaths.isize(); p++){
        // get eorbit of v_inst[repR[permRepR[p]]].v as (c,m)
        vec<pair<int,EL>> Q;
        int64_t ne = v_inst[repR[permRepR[p]]].v;
        while(1){
           Q.push(e_inst[ne].c, e_inst[ne].m);
           ne = estart[e_inst[ne].c]+e_inst[ne].m; 
           if(ne==v_inst[repR[permRepR[p]]].v) break;
        }
        vec<int64_t> RCQ;
        for(const auto& q: Q)
            RCQ.push(estart[dinv[q.first]]+ws[q.first]-q.second-1);
        Sort(RCQ);
        #pragma omp critical
        { dinv_new[p] = Start[RCQ]; }
    }

    if(GET_META){
        if (verbose) PRINTDEETS("get meta data; WARNING: SLOW & HIGH MEM");

        // for each unipath metadata
        merger_supp.clear();
        merger_supp.resize(metadata.size());
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
                        merger_supp[kk].push_back(quad<int,EL,int,EL>(
                            match[i].first, match[i].second, match[i].third,0));
                }
                if(match[i].first!=cl){
                    cl = match[i].first;
                    merger_supp[kk].push_back(quad<int,EL,int,EL>(
                        match[i].first, match[i].second, match[i].third,0));
                }
                merger_supp[kk].back().fourth++;
            }
        }

        if(dir!="")
            BinaryWriter::writeFile(dir+"/a.clpaths",merger_supp);
    }
    // []
    Destroy(metadata);

    // []
    Destroy(dinv);
    Destroy(ws);
    Destroy(e_inst);
    Destroy(estart);

    // ---------------- CONSTRUCT SUPERGRAPH --------------
    if (verbose) PRINTDEETS("constructing supergraph");

    // construct the digraph
    D_new.EdgesMutable().resize(unipaths.size());
    #pragma omp parallel for
    for(unsigned int i = 0; i<unipaths.size(); i++){
        D_new.OMutable(i) = unipaths[i];
    }
    D_new.FromMutable().resize(vmap.size());
    D_new.ToMutable().resize(vmap.size());
    D_new.FromEdgeObjMutable().resize(vmap.size());
    D_new.ToEdgeObjMutable().resize(vmap.size());

    for(unsigned int i = 0; i<unipaths.size(); i++){
        D_new.FromMutable(vmap[v_inst[repL[i]].voID]).push_back(vmap[v_inst[repR[permRepR[i]]].voID]);
        D_new.ToMutable(vmap[v_inst[repR[permRepR[i]]].voID]).push_back(vmap[v_inst[repL[i]].voID]);
        D_new.FromEdgeObjMutable(vmap[v_inst[repL[i]].voID]).push_back(i);
        D_new.ToEdgeObjMutable(vmap[v_inst[repR[permRepR[i]]].voID]).push_back(i);
    }

    // sort syncing
    if(verbose) PRINTDEETS("sortsyncing");

    #pragma omp parallel for
    for(unsigned int i =0; i<vmap.size(); i++){
        SortSync(D_new.FromMutable(i),D_new.FromEdgeObjMutable(i));
        SortSync(D_new.ToMutable(i),D_new.ToEdgeObjMutable(i));
    }

    // --------------- REINTRODUCE GAPS -----------------
    if(verbose) PRINTDEETS("recover gaps");

    // loop through all d, if d[0] >= hbE, then D.OMutable(d) = Gaps[d[0]];
    for(int64_t d = 0; d<D_new.E(); d++){
        if (D_new.O(d).size() == 0) continue;
        T start = D_new.O(d)[0];
        if (start < hbE ) continue;
        D_new.OMutable(d) = Gaps[start];
    }

    // check inversion
    int unset = 0, uninv = 0;
    for ( int i = 0; i < dinv_new.isize( ); i++ ){
        if ( dinv_new[i] < 0 ) unset++;
        else if ( dinv_new[dinv_new[i]] != i ) uninv++;
    }
    ForceAssertEq( unset, 0 );
    ForceAssertEq( uninv, 0 ); 


    // []
    Gaps.clear();
    Destroy(v_inst);
    Destroy(repR);
    Destroy(repL);
    Destroy(permRepR);
    vmap.clear();
}

template void GlueGraphs(digraphE<vec<int>>& D, vec<int>& dinv, 
    vec<nptuple>& Matches, int hbE, digraphE<vec<int>>& D_new, vec<int>& dinv_new,
    const Bool GET_META, 
    vec<vec<quad<int,uint16_t,int,uint16_t>>>& merger_supp,
    const Bool verbose, String dir);

template void GlueGraphs(digraphE<vec<unsigned char>>& D, 
    vec<int>& dinv, vec<nptuple>& Matches, unsigned char hbE, 
    digraphE<vec<unsigned char>>& D_new, vec<int>& dinv_new, 
    const Bool GET_META, 
    vec<vec<quad<int,uint16_t,int,uint16_t>>>& merger_supp,
    const Bool verbose, String dir);
