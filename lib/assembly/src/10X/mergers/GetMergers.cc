// Copyright of 10x Genomics (2016)

#include "10X/mergers/GetMergers.h"
#include "10X/mergers/NicePrints.h"

nptuple make_tup(int d1, int p1, int d2, int p2, int len){
    return make_triple(make_pair(d1,p1),make_pair(d2,p2),len);
}


// can be made smarter. 
void MarkOverLaps(const vec<int>& x1, const vec<int>& x2, const vec<int>& mult, 
        const HyperBasevectorX& hb, const int MIN_GLUE, 
        const int MAX_MULT, vec<olap>& Olaps, Bool verbose){
    
    // Dynamic program
    vec<vec<int>> matrix;
    matrix.resize(x1.size()+1);
    for(int i1 = 0; i1<x1.isize()+1; i1++)
        matrix[i1].resize(x2.size()+1,0);

    // artificially raise the bar!
    vec<int> m1(x1);
    for(int i1 = 0; i1<x1.isize(); i1++)
        m1[i1] = (mult[x1[i1]] > MAX_MULT) ? 1:0;

    for(int i1 = 0; i1 <x1.isize(); i1++){
        for(int i2 = 0; i2<x2.isize(); i2++){
            if(x1[i1]==x2[i2])
                matrix[i1+1][i2+1] = (m1[i1]) ? 0 : matrix[i1][i2] + hb.Kmers(x1[i1]);
        }
    }

    // gather longest common substrings that are acceptable (min kmer, max mult) and unique
    // filter out repeats

    vec<triple<int,int,int>> Suffs;
    vec<int> lens;
    for(int i1 = 0; i1 <x1.isize(); i1++){
        for(int i2 = 0; i2<x2.isize(); i2++){
            if(matrix[i1+1][i2+1]>0){
                Suffs.push_back(make_triple(matrix[i1+1][i2+1],i1,i2));
                matrix[i1+1][i2+1] = 0; // mark as read
                int k = 1;
                while(i1+k < x1.isize() && i2+k< x2.isize() && matrix[i1+1+k][i2+1+k]>0){
                    Suffs.back() = make_triple(matrix[i1+1+k][i2+1+k],i1+k,i2+k);
                    matrix[i1+1+k][i2+1+k] = 0; // mark as read
                    k++;
                }
                if (Suffs.back().first < MIN_GLUE/2) // lenient here!!!!
                    Suffs.pop_back(); // prune by length
                else
                    lens.push_back(k);
            }
        }
    }

    // sort by kmer size of match
    SortSync(Suffs,lens);

    if(Suffs.size()){
        // TODO:choose to keep longest overlap, or longest two, 
        // or all in same O.O and > MIN_GLUE
        auto& t = Suffs.back();
        int L = lens.back();
        Olaps.push_back(olap(t.second-L+1,t.third-L+1,L));

        // pretty print all overlaps
        if(verbose){
            int max_mult = 0;
            for(auto e: x1)
                max_mult = Max(max_mult,mult[e]);
            int num = 0;
            for(int i = 0; i<Suffs.isize(); i++){
                auto& t1 = Suffs[i];
                num++;
                cout<<"Common Seq #"<<num<<" : ";
                olap O(t1.second-lens[i]+1,t1.third-lens[i]+1,lens[i]);
                for(int r = 0; r<O.len; r++)
                    cout << "\033[1;31m "<<x1[r+O.start1]<<"\033[0m ";
                cout<<endl<<"Max multiplicity of a base edge : "<<max_mult<<endl<<endl;
            }
        }
    }else{
        if(verbose){
            int max_mult = 0;
            for(auto e: x1)
                max_mult = Max(max_mult,mult[e]);
            cout<<"No large overlaps"<<endl;
            cout<<"Max multiplicity of a base edge : "<<max_mult<<endl<<endl;
        }
    }
}

namespace MergeHelper{
void DebugPairs(qtuple Pairs, const HyperBasevectorX& hbx, const vec<vec<int>>& DLENS, const vec<digraphE<vec<int>>>& DD)
{
     vec<int> hb1(DD[Pairs.a1].O(Pairs.d1));
     vec<int> hb2(DD[Pairs.a2].O(Pairs.d2));
     vec<int> k1, k2;
     for(auto e:hb1)
        k1.push_back(hbx.Kmers(e));
     for(auto e:hb2)
        k2.push_back(hbx.Kmers(e));

     vec<int> D1 = hb1, D2 = hb2;
     Sort(hb1); Sort(hb2);
     vec<int> ce,ke;
     Intersection(hb1,hb2,ce);
     hb1 = D1;
     hb2 = D2;
     for(auto e: ce)
         ke.push_back(hbx.Kmers(e));
     int sum =0;
     for(auto k:ke)
         sum += k;
     vec<int> ce1 = ce, ce2 = ce;

     cout<<"----->"<<endl;
     cout<<"D1:"<<Pairs.d1<<", A1:"<<Pairs.a1<<"; ";
     cout<<"D2:"<<Pairs.d2<<", A2:"<<Pairs.a2<<endl;
     cout<<"Base edge seq : ";
     for(auto d: D1){
         if(BinPosition(ce1,d)>-1){
             cout << "\033[1;31m "<<d<<"\033[0m ";
             hb1[Position(hb1,d)] = -1;
             ce1[BinPosition(ce1,d)] = -1;
             Sort(ce1);
         }
         else
             cout<<d<<" ";
     }
     cout<<endl;
 
     cout<<"Base edge seq : ";
     for(auto d: D2){
         if(BinPosition(ce2,d)>-1){
             cout << "\033[1;31m "<<d<<"\033[0m ";
             hb2[Position(hb2,d)] = -1;
             ce2[BinPosition(ce2,d)] = -1;
             Sort(ce2);
         }else
             cout<<d<<" ";
     }
     cout<<endl;

     cout<<"Total kmers:"<<DLENS[Pairs.a1][Pairs.d1]<<" : ";
     for(int i = 0; i< k1.isize(); i++){
         int d = k1[i];
         if(hb1[i]<0)
             cout << "\033[1;31m "<<d<<"\033[0m ";
         else
             cout<<d<<" ";
     }
     cout<<endl;

     cout<<"Total kmers:"<<DLENS[Pairs.a2][Pairs.d2]<<" : ";
     for(int i = 0; i< k2.isize(); i++){
         int d = k2[i];
         if(hb2[i]<0)
             cout << "\033[1;31m "<<d<<"\033[0m ";
         else
             cout<<d<<" ";
     }
     cout<<endl;
     cout<<"Total common kmers: "<<sum;
     cout<<endl;
     cout<<"<-----"<<endl;
     cout<<endl;
 }

template<class Iter, class T>
Iter binary_findL(Iter begin, Iter end, T val)
{
    return std::lower_bound(begin, end, val);
}

template<class Iter, class T>
Iter binary_findR(Iter begin, Iter end, T val)
{
    return std::upper_bound(begin, end, val);
}

Bool NumCommonKmers(const vec<int>& x1, const vec<int>& x2, 
        const HyperBasevectorX& hb, const int MIN_GLUE){

    int sz = 0;
    if(x2.size()==0) return False;
    if(x1.size()==0) return False;

    auto s1 = binary_findL(x1.begin(), x1.end(), x2[0]);
    auto e1 = binary_findR(s1, x1.end(), x2.back());
    auto s2 = binary_findL(x2.begin(), x2.end(), x1[0]);
    auto e2 = binary_findR(s2, x2.end(), x1.back());

    auto i1 = s1;
    auto i2 = s2;
    while(i1!=e1 && i2!=e2){
        if(*i1<*i2) i1++;
        else{
            if(*i2<*i1) i2++;
            else{
                sz += hb.Kmers(*i1);
                if(sz>=MIN_GLUE/2) return True;
                i1++;i2++;
            }
        }
    }
    return False;
}

Bool NumCommonKmers_slow(const vec<int>& x1, const vec<int>& x2, 
        const HyperBasevectorX& hb, const int MIN_GLUE){

    int i1 = 0, i2 = 0, sz = 0;
    while(i1<x1.isize() && i2<x2.isize()){
        if(x1[i1]<x2[i2]){
            i1++;
            continue;
        }
        if(x2[i2]<x1[i1]){
            i2++;
            continue;
        }
        if(x1[i1]==x2[i2]){
            sz += hb.Kmers(x1[i1]);
            i1++;i2++;
            if(sz>=MIN_GLUE/2)
                return True;
            continue;
        }
    }
    return False;
}

Bool MarkOverLapsLenient(const vec<int>& x1, const vec<int>& x2, const HyperBasevectorX& hb, 
        const int MIN_GLUE){
    
    // Dynamic program
    vec<vec<int>> matrix;
    matrix.resize(x1.size()+1);
    for(int i1 = 0; i1<x1.isize()+1; i1++)
        matrix[i1].resize(x2.size()+1,0);

    for(int i1 = 0; i1 <x1.isize(); i1++){
        for(int i2 = 0; i2<x2.isize(); i2++){
            if(x1[i1]==x2[i2])
                matrix[i1+1][i2+1] = matrix[i1][i2] + hb.Kmers(x1[i1]);
            if(matrix[i1+1][i2+1]>=MIN_GLUE/2) // linient
                return True;
        }
    }
    return False;
}

void Homogenize_DPI(vec<vec<int>>& BPATHS_INDEX, const vec<int>& dinv, const int numE){
    for(int d = 0; d<numE; d++){
        int d1 = min(d,dinv[d]);
        int d2 = max(d,dinv[d]);
        if(d1!=d2){
            if(d1==d){ // do this only when smallest
                BPATHS_INDEX[d1].insert(BPATHS_INDEX[d1].end(),
                BPATHS_INDEX[d2].begin(),
                BPATHS_INDEX[d2].end());
                UniqueSort(BPATHS_INDEX[d1]);
                BPATHS_INDEX[d2] = BPATHS_INDEX[d1];
            }
        }else
            UniqueSort(BPATHS_INDEX[d]);
    }
}

void MarkIgnores(vec<Bool>& IgnoreDs, const vec<vec<int>>& BPATHS_INDEX, 
    const vec<int>& DLENS, const digraphE<vec<int>>& DL, 
    const int MIN_COMMON_BC, const int MIN_GLUE){
    IgnoreDs.resize(DL.E(),True);
    for(int d = 0; d < DL.E(); d++){
        if(DL.O(d)[0]<0) continue;
        if(DLENS[d] < MIN_GLUE) continue;
        if(BPATHS_INDEX[d].isize() < MIN_COMMON_BC) continue;
        IgnoreDs[d] = False;
    }
}

Bool EnoughCommons(const vec<int>& x1, const vec<int>& x2, 
        const int threshold){

    int sz = 0;
    if(x2.size()==0) return False;
    if(x1.size()==0) return False;

    auto s1 = binary_findL(x1.begin(), x1.end(), x2[0]);
    auto e1 = binary_findR(s1, x1.end(), x2.back());
    auto s2 = binary_findL(x2.begin(), x2.end(), x1[0]);
    auto e2 = binary_findR(s2, x2.end(), x1.back());

    auto i1 = s1;
    auto i2 = s2;
    while(i1!=e1 && i2!=e2){
        if(*i1<*i2) i1++;
        else{
            if(*i2<*i1) i2++;
            else{
                sz++; 
                if(sz>=threshold) return True;
                i1++;i2++;
            }
        }
    }
    return False;
}

Bool EnoughCommons_slow(const vec<int>& B1, const vec<int>& B2, 
        const int threshold){
    int num_common_bcs=0;
    int ib1 = 0, ib2 = 0;
    while(ib1<B1.isize() && ib2<B2.isize()){
        if(B1[ib1]<B2[ib2]){
            ib1++;
            continue;
        }
        if(B1[ib1]>B2[ib2]){
            ib2++;
            continue;
        }
        if(B1[ib1]==B2[ib2]){
            ib1++;ib2++;
            num_common_bcs++;
            if(num_common_bcs>=threshold){
                return True;
            }
            continue;
        }
    }
    return False;
}

void GreatIntersector1(const vec<pair<int,int>>& crd, const vec<vec<vec<int>>>& BPATHS, 
        const vec<vec<int>>& DLI, vec<qtuple>& FD1D2s, const vec<Bool>& IgnoreD1s, 
        const vec<Bool>& IgnoreD2s, const vec<vec<vec<int>>>& sDL, 
        const vec<vec<int>> BPATHS_INDEX_A1, const vec<vec<int>>& BPATHS_INDEX_A2, 
        const HyperBasevectorX& hbx, const int A1, const int A2, 
        const int MIN_GLUE, const int MIN_COMMON_BC){

    std::unordered_set<pair<int,int>> LD1D2s;

    for(int ar = 0; ar < crd.isize(); ar++){

        int64_t bc1 = crd[ar].first;
        int64_t bc2 = crd[ar].second;
        vec<int> d1s, d2s;

        // filter based on ignores and used up barcodes
        for(auto d : BPATHS[A1][bc1])
            if(!IgnoreD1s[d])
                d1s.push_back(d);
        for(auto d : BPATHS[A2][bc2])
            if(!IgnoreD2s[d])
                d2s.push_back(d);

        for(int d1 = 0; d1< (int) d1s.size(); d1++){
            int D1 = d1s[d1];
            for(int d2 =0; d2< (int) d2s.size(); d2++){
                int D2 = d2s[d2];

                if(LD1D2s.count(make_pair(D1,D2))) continue;
                if(LD1D2s.count(make_pair(DLI[A1][D1],DLI[A2][D2]))) continue;
                LD1D2s.insert(make_pair(D1,D2));
                
                const vec<int>& B2 = BPATHS_INDEX_A2[D2];
                const vec<int>& B1 = BPATHS_INDEX_A1[D1];

                if(EnoughCommons(B1,B2,MIN_COMMON_BC)){
                    if(NumCommonKmers(sDL[A1][D1], sDL[A2][D2], 
                            hbx, MIN_GLUE)){
                        FD1D2s.push_back(qtuple{D1,A1,D2,A2});
                        FD1D2s.push_back(qtuple{DLI[A1][D1],A1,DLI[A2][D2],A2});
                    }
                }
            }
        }
    }
}

void ConstructMaps(const vec<Bool>& IgnoreD1s, const vec<Bool>& IgnoreD2s, 
        vec<int>& mapD1s, vec<int>& imapD1s, vec<int>& mapD2s, 
        vec<int>& imapD2s, vec<int>& LD1D2s, 
        const vec<pair<int,int>>& crd, const vec<vec<int>>& BPATHS, 
        const vec<int>& DLI, Bool constructor2 = True){

    imapD2s.reserve(IgnoreD2s.size());
    mapD2s.resize(IgnoreD2s.size(),-1);
    for(int i = 0; i<IgnoreD2s.isize(); i++)
        if(!IgnoreD2s[i])
            imapD2s.push_back(i);
    // mapd2s
    for(int i = 0; i<imapD2s.isize(); i++)
        mapD2s[imapD2s[i]]=i;

    if(constructor2){
        // at least start with A2 sized stuff
        imapD1s.reserve(400000); //heuristic
        mapD1s.resize(IgnoreD1s.size(),-1);
        for(int ar = 0; ar < crd.isize(); ar++){
            for(auto& d : BPATHS[crd[ar].first]){
                if(!IgnoreD1s[d] && mapD1s[d]==-1){
                    mapD1s[d]=1;
                    mapD1s[DLI[d]]=1;
                }
            }
        }
        int idx = 0;
        for(int i = 0; i < mapD1s.isize(); i++){
            if(mapD1s[i]==1){
                imapD1s.push_back(i);
                mapD1s[i] = idx; idx++;
            }
        }
        // size the [][]
        LD1D2s.resize(imapD1s.size()*imapD2s.size(),0);
    }else{
        // at least start with A2 sized stuff
        imapD1s.reserve(IgnoreD1s.size());
        mapD1s.resize(IgnoreD1s.size(),-1);
        for(int i = 0; i<IgnoreD1s.isize(); i++)
            if(!IgnoreD1s[i])
                imapD1s.push_back(i);
        // size the [][]
        LD1D2s.resize(imapD1s.size()*imapD2s.size(),0);
        // mapd1s
        for(int i = 0; i<imapD1s.isize(); i++){
            mapD1s[imapD1s[i]]=i;
        }
    } 
}

void GreatIntersector(const vec<pair<int,int>>& crd, const vec<vec<vec<int>>>& BPATHS, 
        const vec<vec<int>>& DLI, vec<qtuple>& FD1D2s, const vec<Bool>& IgnoreD1s, 
        const vec<Bool>& IgnoreD2s, const vec<vec<vec<int>>>& sDL, 
        const vec<vec<int>>& BPATHS_INDEX_A1, const vec<vec<int>>& BPATHS_INDEX_A2, 
        const HyperBasevectorX& hbx, const int A1, const int A2, 
        const int MIN_GLUE, const int MIN_COMMON_BC, Bool constructor2=True){

    vec<int> LD1D2s;
    vec<int> mapD1s, mapD2s, imapD2s, imapD1s;
    ConstructMaps(IgnoreD1s, IgnoreD2s, mapD1s, imapD1s, mapD2s, imapD2s,
            LD1D2s, crd, BPATHS[A1], DLI[A1], constructor2);

    vec<int> d1s, d2s;
    for(int ar = 0; ar < crd.isize(); ar++){
        d1s.clear();
        d2s.clear();
        // filter based on ignores and used up barcodes
        for(auto d : BPATHS[A1][crd[ar].first])
            if(!IgnoreD1s[d])
                d1s.push_back(d);
        for(auto d : BPATHS[A2][crd[ar].second])
            if(!IgnoreD2s[d])
                d2s.push_back(d);
        for(int d1 = 0; d1< (int) d1s.size(); d1++)
            for(int d2 =0; d2< (int) d2s.size(); d2++)
                LD1D2s[mapD1s[d1s[d1]]*imapD2s.size() + mapD2s[d2s[d2]]]++;
    }

    // process all entries
    vec<pair<int,int>> inst;
    for(int i = 0; i< (int) imapD1s.size(); i++){
        for(int j = 0; j< (int) imapD2s.size(); j++){
            if (LD1D2s[i*imapD2s.size()+j]<1) continue;
            int D1 = imapD1s[i];
            int D2 = imapD2s[j];
            int cnt = 0;
            inst.clear();
            inst.push_back(make_pair(D1,D2));
            inst.push_back(make_pair(DLI[A1][D1],D2));
            inst.push_back(make_pair(D1,DLI[A2][D2]));
            inst.push_back(make_pair(DLI[A1][D1],DLI[A2][D2]));
            UniqueSort(inst);
            for(auto& q: inst){
                if(LD1D2s[mapD1s[q.first]*imapD2s.size()+mapD2s[q.second]]>0){
                    cnt += LD1D2s[mapD1s[q.first]*imapD2s.size()+mapD2s[q.second]];
                    LD1D2s[mapD1s[q.first]*imapD2s.size()+mapD2s[q.second]] = 0;
                }else
                    q.first = -1; //mark non existent
            }
            
            set<pair<int,int>> repeats;
            if(cnt>=MIN_COMMON_BC){
                for(auto& q: inst){
                    if(q.first>-1){
                        pair<int,int> PP(q.first,q.second),
                        QQ(DLI[A1][q.first],DLI[A2][q.second]);
                        if(repeats.count(PP)) continue; // already done!
                        if(EnoughCommons(BPATHS_INDEX_A1[q.first],
                                    BPATHS_INDEX_A2[q.second],MIN_COMMON_BC)){
                            if(NumCommonKmers(sDL[A1][q.first], 
                                        sDL[A2][q.second], hbx, MIN_GLUE)){
                                FD1D2s.push_back(qtuple{q.first,A1,q.second,A2});
                                repeats.insert(PP);
                                if(repeats.count(QQ)==0){
                                    repeats.insert(QQ);
                                    FD1D2s.push_back(qtuple{DLI[A1][q.first],A1,
                                        DLI[A2][q.second],A2});
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}


}

static void TranslatePairsToMatches(const vec<digraphE<vec<int>>>& DL, 
        const vec<vec<int>> & DLI, const HyperBasevectorX& hbx, 
        vec<qtuple>& Pairs, vec<nptuple>& Matches, const int MIN_GLUE,
        const int MAX_MULT, String OUTDIR, Bool verbose ){

     if (verbose) PRINTDEETS("translating pairs to matches");

     // calculate mult
     vec<int> mult;
     mult.resize( hbx.E( ), 0 );
     for(int j = 0; j<DL.isize(); j++){
         for ( int d = 0; d < DL[j].E( ); d++ )
         {    if ( DL[j].O(d)[0] < 0 ) continue;
              for ( int i = 0; i < DL[j].O(d).isize( ); i++ )
                   mult[ DL[j].O(d)[i] ]++;
         }
     }
 
     /* // debug */
     /* if (!verbose){ */
     /*     for(int i = 0; i < Min(Pairs.isize(),1000); i++){ */
     /*        cout<<"Candidate pair :"<<i<<endl; */
     /*        vec<olap> Olaps; */
     /*        qtuple& p = Pairs[i]; */
     /*        MergeHelper::DebugPairs(p, hbx, DLENS, DL); */
     /*        MarkOverLaps(DL[p.a1].O(p.d1), DL[p.a2].O(p.d2), mult, */ 
     /*                hbx, MIN_GLUE, MAX_MULT, Olaps, True); */
     /*     } */
     /* } */

     // flatten Ds
     // get total # closures
     int64_t N = 0;
     vec<int64_t> dci; dci.push_back(0);
     vec<int> dc;
     for(int j = 0; j< DL.isize() ; j++){
         N += DL[j].E();
         dci.push_back(N);
         for(int d = 0; d< DL[j].E(); d++)
             dc.push_back(j); 
     } 
     BinaryWriter::writeFile(OUTDIR+"/a.dci",dci);

     if(verbose) 
         PRINTQUANTS("original assembly #E = ",ToStringAddCommas(DL.back().E()));

     // Translate to matches from candidate pairs
     #pragma omp parallel for schedule(dynamic,1)
     for(int64_t i = 0; i < (int64_t)Pairs.size(); i++){
         int64_t l1 =0, l2=0, l1i =0, l2i =0;
         const qtuple & p = Pairs[i];
         vec<olap> Olaps;
         MarkOverLaps(DL[p.a1].O(p.d1),DL[p.a2].O(p.d2),
                 mult,hbx,MIN_GLUE,MAX_MULT,Olaps);
         for(auto& ol: Olaps){
             // sanity check - gaps are removed from pairs!
             if( DL[p.a1].O(p.d1)[0] < 0 || DL[p.a2].O(p.d2)[0] < 0)
                 cout<<"Gaps detected in candidate pairs -- likely gonna fail!"<<endl;
             #pragma omp critical
             { 
                 Matches.push_back(make_tup(dci[p.a1]+p.d1,ol.start1,
                         dci[p.a2]+p.d2,ol.start2,ol.len));
                 Matches.push_back(make_tup(dci[p.a1]+DLI[p.a1][p.d1],
                         DL[p.a1].O(DLI[p.a1][p.d1]).size()-ol.start1-ol.len,
                         dci[p.a2]+DLI[p.a2][p.d2],
                         DL[p.a2].O(DLI[p.a2][p.d2]).size()-ol.start2-ol.len,
                         ol.len));
             }
         }
     }

     if(verbose) PRINTDEETS("sorting");
     ParallelUniqueSort(Matches);
     if(verbose) PRINTDEETS("finished translation");
}

void GetMergerCandidates(const vec<digraphE<vec<int>>>& DL, const vec<vec<int>>& DLI, 
        vec<vec<int>>& BCLIST, vec<vec<vec<int>>>& BPATHS, const vec<int>& bc, 
        const vec<int64_t>& bci, const HyperBasevectorX& hbx, vec<nptuple>& Matches, 
        int MIN_GLUE, int MAX_MULT, int MIN_BC_OLAP, 
        int MIN_COMMON_BC, String OUTDIR, Bool verbose)
{
     double clock = WallClockTime();
     int bs_size = DL.isize();
               
     // DLENS
     if(verbose) PRINTDEETS("creating dlens");
    
     vec<vec<int>> DLENS(bs_size);
     #pragma omp parallel for schedule(dynamic,1)
     for(int z =0; z<bs_size; z++){
         auto& Dl = DL[z];
         DLENS[z].resize(Dl.E(),0);
         for(int e =0; e<Dl.E(); e++){
             if(Dl.O(e)[0]<0) continue;
             for(int j =0; j<Dl.O(e).isize(); j++)
                 DLENS[z][e] += hbx.Kmers(Dl.O(e)[j]);
         }
     }
 
     // Start building other data structs
     std::set<pair<int,int>> APairSet;
     vec<pair<int,int>> APairs;

     if (verbose)
         PRINTDEETS("marking highly overlapping assemblies among "
                    <<bs_size<<" assemblies");

     #pragma omp parallel for schedule(dynamic,1)
     for(int i1 = 0; i1<bs_size; i1++){
         vec<int> common_bcs;
         
         for(int i2 = i1+1; i2<bs_size; i2++){
             common_bcs.clear();
             
            Bool enuf = False;
            if(BCLIST[i2].isize()==bci.isize()-1) // global assembly
                enuf = (BCLIST[i1].isize()>=MIN_BC_OLAP);
            else{
                enuf = MergeHelper::EnoughCommons(BCLIST[i1],BCLIST[i2],MIN_BC_OLAP);
            }
             #pragma omp critical
             { 
                 if (enuf){
                     // sort by bclist size
                     if(BCLIST[i1].size()>BCLIST[i2].size())
                         APairSet.insert(make_pair(i1,i2));
                     else
                         APairSet.insert(make_pair(i2,i1));
                 }
             }
         }
     }
     for(auto& pr: APairSet)
         APairs.push_back(pr);
     APairSet.clear();

     // build APairsI
     vec<int> APairsI;
     int idx = 0;
     int prev = -1; 
     for(auto pr : APairs){
         if(prev != pr.first){
             APairsI.push_back(idx);
             prev = pr.first;
         }
         idx++;
     }
     APairsI.push_back(idx);

     if (verbose){
         cout<<"size of pairs of overlapping assemblies with a "
             <<" threshold MIN_BC_OLAP:"<< MIN_BC_OLAP
             <<" is :"<<APairs.size()<<"/"
             <<(((int64_t) bs_size)*(bs_size-1))/2<<endl; 
     }

     int64_t SUM = 0;
     for(auto& Dl: DL)
         SUM += Dl.E();
     int stopped = 0, ndots = 0;

     if (verbose) PRINTDEETS("sort kmers in edges");

     vec<vec<vec<int>>> sDL(bs_size);
     #pragma omp parallel for schedule(dynamic,1)
     for(int a = 0; a < bs_size; a++){
         sDL[a] = DL[a].Edges();
         for(int d = 0; d< DL[a].E(); d++)
             Sort(sDL[a][d]);
     }

     if (verbose) PRINTDEETS("marking pairs of superedges for overlap calculation");

     if(APairsI.size()==1) // just one assembly
     {
         //[]
         Destroy(BPATHS);
         Destroy(BCLIST);
         return;
     }

     vec<qtuple> Pairs;

     // handle noglobals option internally
     Bool NOGLOBALS = False;
     if(BCLIST[APairs.back().first].size()!=bci.size()-1)
            NOGLOBALS = True;

     int nastybit = 0;
     if(!NOGLOBALS)
         nastybit = APairsI[APairsI.size()-1]-APairsI[APairsI.size()-2];
     int64_t APrange = (NOGLOBALS) ?  APairsI.size()-1 : APairsI.size()-2;

     double mclock=WallClockTime();
/* #define _jumpso_ */
#ifdef _jumpso_
     goto jumper;
#endif
     #pragma omp parallel for schedule(dynamic,1)
     for(int64_t ia = 0; ia < APrange; ia++){
         int A1 = APairs[APairsI[ia]].first;

         vec<vec<int>> BPATHS_INDEX_A1;
         invert(BPATHS[A1],BPATHS_INDEX_A1,DL[A1].E());
         MergeHelper::Homogenize_DPI(BPATHS_INDEX_A1,DLI[A1],DL[A1].E());
         for(int d = 0; d<DL[A1].E(); d++)
             for(int b = 0; b< (int) BPATHS_INDEX_A1[d].size(); b++)
                 BPATHS_INDEX_A1[d][b] = BCLIST[A1][BPATHS_INDEX_A1[d][b]];

         vec<pair<int,int>> crd;
         vec<vec<int>> BPATHS_INDEX_A2;

         // get ignores
         vec<Bool> IgnoreD1s;
         MergeHelper::MarkIgnores(IgnoreD1s,BPATHS_INDEX_A1, DLENS[A1], 
                 DL[A1], MIN_COMMON_BC, MIN_GLUE);
         for(int64_t ap = APairsI[ia]; ap < (int64_t) APairsI[ia+1]; ap++){
             
            int A2 = APairs[ap].second;

            BPATHS_INDEX_A2.clear();
            invert(BPATHS[A2],BPATHS_INDEX_A2,DL[A2].E());
            MergeHelper::Homogenize_DPI(BPATHS_INDEX_A2,DLI[A2],DL[A2].E());
            for(int d = 0; d<DL[A2].E(); d++)
                for(int b = 0; b< (int) BPATHS_INDEX_A2[d].size(); b++)
                    BPATHS_INDEX_A2[d][b] = BCLIST[A2][BPATHS_INDEX_A2[d][b]]; 

            // ON THE FLY
            crd.clear();
            int i1 = 0, i2 = 0;
            while(i1< BCLIST[A1].isize() && i2 < BCLIST[A2].isize()){
                if(BCLIST[A1][i1] < BCLIST[A2][i2]){
                    i1++;
                    continue;
                }
                if(BCLIST[A2][i2] < BCLIST[A1][i1]){
                    i2++;
                    continue;
                }
                if(BCLIST[A1][i1] == BCLIST[A2][i2]){
                    crd.push_back(make_pair(i1,i2));
                    i1++; i2++;
                    continue;
                }
            }

            vec<Bool> IgnoreD2s;
            MergeHelper::MarkIgnores(IgnoreD2s,BPATHS_INDEX_A2, DLENS[A2], 
                    DL[A2], MIN_COMMON_BC, MIN_GLUE);
            vec<qtuple> FD1D2s;

            MergeHelper::GreatIntersector(crd, BPATHS, DLI, FD1D2s, 
                    IgnoreD1s, IgnoreD2s, sDL, BPATHS_INDEX_A1, 
                    BPATHS_INDEX_A2, hbx, A1, A2, MIN_GLUE, MIN_COMMON_BC);

            #pragma omp critical
            {
                Pairs.insert(Pairs.end(),FD1D2s.begin(),FD1D2s.end());
                MakeDots(stopped, ndots, APairs.isize() - nastybit);
            }
         } 
     } 
#ifdef _jumpso_
jumper:
#endif
     if(NOGLOBALS)
         goto skipnasty;

     // do last one separately as it is bad
     if(verbose) PRINTDEETS("take care of nasty assembly");

     stopped = 0; ndots = 0;

     for(int64_t ia = APairsI.size()-2; ia < (int64_t) APairsI.size()-1; ia++){

         int A1 = APairs[APairsI[ia]].first;

         vec<vec<int>> BPATHS_INDEX_A1;
         invert(BPATHS[A1],BPATHS_INDEX_A1,DL[A1].E());
         MergeHelper::Homogenize_DPI(BPATHS_INDEX_A1,DLI[A1],DL[A1].E());
         for(int d = 0; d<DL[A1].E(); d++)
             for(int b = 0; b< (int) BPATHS_INDEX_A1[d].size(); b++)
                 BPATHS_INDEX_A1[d][b] = BCLIST[A1][BPATHS_INDEX_A1[d][b]];

         vec<Bool> IgnoreD1s;
         MergeHelper::MarkIgnores(IgnoreD1s, BPATHS_INDEX_A1, DLENS[A1], 
                 DL[A1], MIN_COMMON_BC, MIN_GLUE);

         // stats
         Bool printstats = False;
         if(printstats){
             int64_t sz = 0, mk = 0, dl = 0, el = 0;
             for(int d = 0; d <DL[A1].E(); d++){
                 if(IgnoreD1s[d])
                     mk++;
                 else{
                     sz += BPATHS_INDEX_A1[d].isize();
                     dl += DLENS[A1][d];
                 }
             }
             for(int bb = 0; bb<BCLIST[A1].isize(); bb++)
                for(auto d : BPATHS[A1][bb])
                    if (!IgnoreD1s[d])
                        el++;
             if(verbose){
                 cout<<Date( )<<": ignore fraction: "<<((double) mk)/DL[A1].E()
                 <<"; mean bc per valid edge: "<< ((double)sz)/(DL[A1].E()-mk)
                 <<"; mean kmer size per valid edge: "<<((double)dl)/(DL[A1].E()-mk)
                 <<"; mean valid edges per barcode: "<<((double)el)/BCLIST[A1].size()<<endl;
             }
         }

         if(verbose) PRINTDEETS("finished precomputation for big guy");
    
         #pragma omp parallel for schedule(dynamic,10)
         for(int64_t ap = APairsI[ia]; ap < (int64_t) APairsI[ia+1]; ap++){
             
            int A2 = APairs[ap].second;
            /* #pragma omp critical */
            /* { cout << "A2 "<<A2<<endl; } */
            vec<vec<int>> BPATHS_INDEX_A2;
            invert(BPATHS[A2],BPATHS_INDEX_A2,DL[A2].E());
            MergeHelper::Homogenize_DPI(BPATHS_INDEX_A2,DLI[A2],DL[A2].E());
            for(int d = 0; d<DL[A2].E(); d++)
                for(int b = 0; b< (int) BPATHS_INDEX_A2[d].size(); b++)
                    BPATHS_INDEX_A2[d][b] = BCLIST[A2][BPATHS_INDEX_A2[d][b]];

            // ON THE FLY
            vec<pair<int,int>> crd;
            int i2 = 0;
            int i1 = 0;
            // sanity check
            ForceAssertEq(BCLIST[A1].size(),bci.size()-1);
            while(i2< BCLIST[A2].isize()){
                crd.push_back(make_pair(BCLIST[A2][i2],i2));
                i2++;
            }

            vec<Bool> IgnoreD2s; 
            MergeHelper::MarkIgnores(IgnoreD2s, BPATHS_INDEX_A2, DLENS[A2], 
                    DL[A2],MIN_COMMON_BC, MIN_GLUE);
            vec<qtuple> FD1D2s;
            
            MergeHelper::GreatIntersector(crd, BPATHS, DLI, FD1D2s, 
                    IgnoreD1s, IgnoreD2s, sDL, BPATHS_INDEX_A1, 
                    BPATHS_INDEX_A2, hbx, A1, A2, MIN_GLUE, MIN_COMMON_BC);

            #pragma omp critical
            {
                Pairs.insert(Pairs.end(),FD1D2s.begin(),FD1D2s.end());
                MakeDots(stopped, ndots, nastybit);
            }
         } 
     }
skipnasty:
     if (verbose) PRINTDEETS("done, used "<<TimeSince(mclock));

     // []
     Destroy(BPATHS);
     Destroy(BCLIST); 

     if(verbose) 
         PRINTDEETS("built a dictionary of candidate d1d2 of size: " << 
                    Pairs.size() <<" in "<<SUM<<" total superedges");

     // translate pairs
     TranslatePairsToMatches(DL, DLI, hbx, Pairs, Matches,
             MIN_GLUE, MAX_MULT, OUTDIR, verbose );

     if (verbose) PRINTDEETS("done, used "<<TimeSince(clock));
     
     // []
     Destroy(DLENS);
}


