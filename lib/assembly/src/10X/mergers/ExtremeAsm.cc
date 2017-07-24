// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>

#include "10X/mergers/ExtremeAsm.h"
#include "10X/mergers/NicePrints.h"
#include "10X/mergers/NucleateGraph.h"
#include "10X/mergers/NucleateGraph_HP.h"
#include "10X/mergers/GraphEditors.h"
#include "10X/mergers/EdgeSupport.h"


// Smith-Waterman kind of aligner
void MarkAlignments(vec<unsigned char> x1, vec<unsigned char> x2, 
        const vec<int> weight, const int MIN_WT, 
        vec<olap>& Olaps, const vec<String>& controls, Bool verbose){
    
    // take care of dups COSTLY!!!
    if(x1==x2){
        Olaps.push_back(olap(0,0,x1.isize()));
        return;
    }

    // Dynamic program
    vec<vec<int>> matrix;
    matrix.resize(x1.size()+1);
    for(int i1 = 0; i1<x1.isize()+1; i1++)
        matrix[i1].resize(x2.size()+1,0);

    for(int i1 = 0; i1 <x1.isize(); i1++){
        for(int i2 = 0; i2<x2.isize(); i2++){
            if(x1[i1]==x2[i2])
                matrix[i1+1][i2+1] = matrix[i1][i2] + weight[x1[i1]];
        }
    }

    // gather longest common substrings that are acceptable 
    // (min kmer, max extreme content) 

    vec<triple<int,int,int>> Suffs;
    vec<int> Lens;
    for(int i1 = 0; i1 <x1.isize(); i1++){
        for(int i2 = 0; i2<x2.isize(); i2++){
            if(matrix[i1+1][i2+1]>0){
                Suffs.push_back(make_triple(matrix[i1+1][i2+1],i1,i2));
                matrix[i1+1][i2+1] = 0; // mark as read
                int k = 1;
                while(i1+k < x1.isize() && 
                        i2+k< x2.isize() && 
                        matrix[i1+1+k][i2+1+k]>0){
                    Suffs.back() = make_triple(matrix[i1+1+k][i2+1+k],i1+k,i2+k);
                    matrix[i1+1+k][i2+1+k] = 0; // mark as read
                    k++;
                }
                if (Suffs.back().first < MIN_WT) // lenient here!!!!
                    Suffs.pop_back(); // prune by length
                else
                    Lens.push_back(k);
            }
        }
    }

    // sort by kmer size of match
    ReverseSortSync(Suffs,Lens);

    // introduce switches as string for now
    std::set<String> switches;
    for(auto con: controls){
        if(con=="ea")
            switches.insert("EndAligns");
        if(con=="nr")
            switches.insert("NonRepeats");
        if(con=="si")
            switches.insert("Single");
    }

    if(Suffs.size()){
        vec<Bool> ignore(Suffs.size(),False);
        int L1 = x1.isize();
        int L2 = x2.isize();
        int ltrim = 5;
        int rtrim = 5;
        if(switches.count("EndAligns")){
            for(int i = 0; i < Suffs.isize(); i++){
                if(ignore[i]) continue;
                // switches
                // accept only if near ends
                // ----------
                //     ----------
                auto& t = Suffs[i];
                int L = Lens[i];
                int s1 = t.second-L+1;
                int s2 = t.third-L+1;
                int e1 = t.second+1;
                int e2 = t.third+1;
                if(e1>L1-rtrim || e2>L2-rtrim){ // 5 is heuristic
                    if(e1>L1-rtrim && e2>L2-rtrim){ // refined
                        if(s1<ltrim && s2<ltrim) continue;
                    }else{
                        if(s1<ltrim || s2<ltrim) continue;
                    }
                }
                ignore[i]=True;
            }
        }
        if(switches.count("NonRepeats")){
            for(int i = 0; i < Suffs.isize(); i++){
                if(ignore[i]) continue;
                auto& t = Suffs[i];
                int L = Lens[i];
                int s1 = t.second-L+1;
                int e1 = t.second+1;
                int A = 0, G = 0, C = 0, T = 0;
                for(int s = s1; s < e1; s++){
                    if(x1[s]==0) A++;
                    if(x1[s]==1) C++;
                    if(x1[s]==2) G++;
                    if(x1[s]==3) T++;
                }
                // ignore if A^n
                // ignore if AT^n
                // ignore if GC^n
                if(1.0*A/L > 0.8 || 1.0*T/L > 0.8 || 1.0*C/L > 0.8 || 1.0*G/L > 0.8) 
                    ignore[i]=True;
                if(1.0*(A+T)/L > 0.8 || 1.0*(G+C)/L > 0.8)
                    ignore[i]=True;
            }
        }
        if(switches.count("Single")){
            int found = -1;
            for(int i = 0; i < ignore.isize(); i++){
                if(!ignore[i] && found==-1)
                    found = 1;
                else
                    ignore[i] = True;
            }
        }

        // create olaps from what remains
        for(int i = 0; i < Suffs.isize(); i++){
            if(ignore[i]) continue;
            auto& t = Suffs[i];
            int L = Lens[i];
            Olaps.push_back(olap(t.second-L+1,t.third-L+1,L));
        }
    }

    if(verbose){
        if(Suffs.size()){
        // pretty print at most top 4 overlaps
            String BaseChar[4] = {"A","C","G","T"};
            int num = 0;
            for(int i = 0; i < Min(Suffs.isize(),4); i--){
                auto& t1 = Suffs[i];
                num++;
                cout<<"Common Seq #"<<num<<" : "<<endl;
                olap O(t1.second-Lens[i]+1,t1.third-Lens[i]+1,Lens[i]);
                for(int kr = 0; kr<O.start1; kr++)
                    cout << BaseChar[x1[kr]];
                for(int r = 0; r<O.len; r++)
                    cout << "\033[1;31m"<<BaseChar[x1[r+O.start1]]<<"\033[0m";
                for(int kr = O.start1+O.len; kr<x1.isize(); kr++)
                    cout << BaseChar[x1[kr]];
                cout<<endl;
                for(int kr = 0; kr<O.start2; kr++)
                    cout << BaseChar[x2[kr]];
                for(int r = 0; r<O.len; r++)
                    cout << "\033[1;31m"<<BaseChar[x2[r+O.start2]]<<"\033[0m";
                for(int kr = O.start2+O.len; kr<x2.isize(); kr++)
                    cout << BaseChar[x2[kr]];
                cout<<endl<<endl;
            }
        }else
            cout<<"No large overlaps"<<endl<<endl;
    }
}

void GetEvilMatches(vec<vec<om<uint16_t>>>& omatch, const int K,
        const vec<vec<unsigned char>>& all_closures, const vec<vec<int>>& ci, 
        const vec<int64_t>& cinv, const vec<String>& controls, Bool verbose)
{
    // consider all edge pairs d1,d2
    // good if d1 and d2 align >= 8 bases over a non repeat region
    // good if d1 and d2 align >= 8 bases over repeat tails and heads
    omatch.resize(all_closures.isize());
   
    // controls
    const int MIN_WT = K;
    vec<int> weight = {1,1,1,1};
    
    #pragma omp parallel for
    for(int i1 = 0; i1 < all_closures.isize(); i1++){
        for (int i2 = i1; i2<all_closures.isize(); i2++){
            const auto& x1 = all_closures[i1];
            const auto& x2 = all_closures[i2];
            vec<olap> Olaps;
            MarkAlignments(x1, x2, weight, MIN_WT, Olaps, controls, verbose);
            // store all good matches and their rcs
            for(auto ol: Olaps){
                if(ol.len == all_closures[i1].isize() && 
                        ol.len==all_closures[i2].isize()) continue;
                #pragma omp critical
                { omatch[i1].push_back(om<uint16_t>(i2,ol.start1,ol.start2,ol.len));
                  omatch[cinv[i1]].push_back(om<uint16_t>(cinv[i2],
                            x1.isize()-ol.start1-ol.len,
                            x2.isize()-ol.start2-ol.len,ol.len)); }
            }
        }
    }

    // reduce size and thus increase speed
    #pragma omp parallel for
    for(int i = 0; i < omatch.isize(); i++)
        UniqueSort(omatch[i]);
}

void CreateMiniAsm(vecbasevector& reads, HyperBasevector& hbl,
     vec<int> invl, const int K){
    BasesToGraph(reads, K, hbl);
    hbl.Involution(invl);
}

class TranslateBases{
    private:
        String BaseChar[4] = {"A","C","G","T"};
    public:
        void CharToBaseString(const vec<unsigned char> edge, String& bs){
            bs = "";
            for(auto e:edge)
                bs += BaseChar[e];
        };
        void BaseStringToChar(String& bstring, vec<unsigned char>& bsint){
            bsint.clear();
            for(auto bs: bstring){
                int bint = -1;
                if(bs=='A')
                    bint = 0;
                if(bs=='C')
                    bint = 1;
                if(bs=='G')
                    bint = 2;
                if(bs=='T')
                    bint = 3;
                bsint.push_back(bint);
            }
        };
};

template <class VBV>
void CharEdges2BaseVec(const digraphE<vec<unsigned char>>& Dp, 
        VBV& superedgeseq){
    
    TranslateBases TB;

    for(int d = 0; d < Dp.E(); d++){
        String bs;
        TB.CharToBaseString(Dp.O(d),bs);
        basevector seq(bs);
        superedgeseq.push_back(seq);
    }
}

template void CharEdges2BaseVec(const digraphE<vec<unsigned char>>& Dp,
        vecbasevector& superedgeseq);
template void CharEdges2BaseVec(const digraphE<vec<unsigned char>>& Dp,
        vec<basevector>& superedgeseq);

template <class VBV>
void Reads2CharClosures(vec<vec<unsigned char>>& all_closuresp, 
        const VBV& reads){ 
    
    TranslateBases TB;
    all_closuresp.clear();
    all_closuresp.reserve(reads.size());
    for(const auto& read: reads){
        String bstring = read.ToString();
        vec<unsigned char> bsint;
        TB.BaseStringToChar(bstring,bsint);
        all_closuresp.push_back(bsint);
    }
}

template void Reads2CharClosures(vec<vec<unsigned char>>& all_closuresp, 
        const vec<basevector>& reads);
template void Reads2CharClosures(vec<vec<unsigned char>>& all_closuresp, 
        const vecbasevector& reads);


template <class VBV>
void ExtremeAsm( VBV& bases, const vec<int64_t> & readlist, 
        digraphE<basevector>& Dl, vec<int>& dinvl, const int K, 
        const vec<String>& controls, const Bool verbose, String out_dir){
    
    vecbasevector reads;
    reads.reserve(readlist.size());
    vec<int64_t> rdlist;
    rdlist.reserve(readlist.size());
    for(const auto rd: readlist)
        rdlist.push_back(abs(rd));
    
    ParallelUniqueSort(rdlist);
    for(auto rd: rdlist)
        reads.push_back(bases[rd]);

    // add inverses
    int64_t rdsz = reads.size();
    for(int64_t i = 0; i<rdsz; i++)
       reads.push_back( ReverseComplement( reads[i] ) );

    ExtremeAsm(reads, Dl, dinvl, K, controls, verbose, out_dir);
}

template void ExtremeAsm( VirtualMasterVec<basevector>& bases, 
        const vec<int64_t>& readlist, digraphE<basevector>& Dl, 
        vec<int>& dinvl, const int K, const vec<String>& controls, 
        const Bool verbose, String out_dir);

template void ExtremeAsm( vecbasevector& bases, const vec<int64_t>& readlist,
        digraphE<basevector>& Dl, vec<int>& dinvl, const int K, 
        const vec<String>& controls, const Bool verbose, String out_dir);

void ExtremeAsm( vecbasevector& reads, digraphE<basevector>& Dl, 
        vec<int>& dinvl, const int K, const vec<String>& controls,
        const Bool verbose, String out_dir)
{
    double clock = WallClockTime();

    // ----------- Convert Reads to Closures --------------------
    if(verbose) PRINTDEETS("creating closures from reads");

    vec<vec<unsigned char>> all_closuresp,all_cl;
    Reads2CharClosures(all_closuresp, reads);
    all_cl = all_closuresp;
    

    // ----------- Do extreme assembly --------------------------
    if(verbose) PRINTDEETS("creating super graph");

    // prepare for assembly building
    unsigned char numE = 4;
    vec<int> invp = {3,2,1,0}; // A = 0, C = 1, G = 2, T = 3
    
    // closure indexing
    vec<int64_t> cinvp,cinv2; 
    vec<vec<int>> cip;
    cip.resize(numE);
    for ( int64_t i = 0; i < (int64_t) all_closuresp.size(); i++ )
        for ( int j = 0; j < all_closuresp[i].isize( ); j++ )
            cip[ all_closuresp[i][j] ].push_back(i);

    #pragma omp parallel for schedule( dynamic, 10 )
    for ( unsigned char e = 0; e < numE; e++ )
    { UniqueSort( cip[e] ); }

    int64_t N = all_closuresp.size();
    cinvp.resize(N,-1);
    for(int kk = 0; kk < N/2; kk++){
        cinvp[kk] = kk+N/2;
        cinvp[kk+N/2] = kk;
    }
    cinv2 = cinvp;

    /* ComputeClosureIndex(invp, cinvp, cip, numE, all_closuresp, verbose); */

    // get matches
    vec<vec<om<uint16_t>>> omatch;
    GetEvilMatches(omatch, K, all_closuresp, cip, cinvp, controls, False);

    // create graph
    vec<int> dinvp;
    digraphE<vec<unsigned char>> Dp;
    Bool GET_META = True;
    vec<vec<quad<int,uint16_t,int,uint16_t>>> read_supp;
    NucleateGraph_HP(omatch, numE, all_closuresp, cinvp, cip, Dp, dinvp,
            GET_META, read_supp, verbose,"",True);

    // ----------- Readpathing for intermediate stage --------------
    if(verbose) PRINTDEETS("read-pathing");
    ReadPathVec paths;
    vec<vec<int>> paths_index;
    GetReadPathsAndIndex(read_supp, paths, paths_index,reads.size());

    // ----------- Do zippering and remove hanging ends ------------
    Bool TINKER = True;
    if(TINKER>0){
        CallGraphEditors(Dp,dinvp,all_cl,cinv2,read_supp,verbose);
    }else{
        if(verbose) PRINTDEETS("clean-up graph");
        BabyZipper(Dp, dinvp, numE, 10, 10);

        vec<int> dels;
        RemoveSimpleHangs(Dp, dinvp,dels);
        ParallelUniqueSort(dels);
        Dp.DeleteEdgesParallel(dels);
        RemoveUnneededVertices(Dp,dinvp);
        CleanupCore(Dp,dinvp);

    }

    // ------------ translations ------------------------------------
    if(verbose) PRINTDEETS("translate graphs to known language");

    // convert edges from char to string to basevec 
    vec<basevector> superedgeseq;
    CharEdges2BaseVec(Dp, superedgeseq);
    auto const& new_from_edge_obj = Dp.FromEdgeObj();
    auto const& new_to_edge_obj = Dp.ToEdgeObj();
    auto const& new_from = Dp.From();
    auto const& new_to = Dp.To();
    Dl.Initialize(new_from,new_to,superedgeseq,new_to_edge_obj,new_from_edge_obj);
    dinvl = dinvp;

    // save graph at out_dir
    if(out_dir!=""){
        if(verbose) PRINTDEETS("saving graphs to disk");
        Mkdir777(out_dir);
        BinaryWriter::writeFile(out_dir+"/"+"a.sup_loc",Dl);
        BinaryWriter::writeFile(out_dir+"/"+"a.sup_loc.inv",dinvl);
    }
    
    if(verbose) { cout << Date() << ": took "<<TimeSince(clock)<<endl; }
    
}

void SanitizeMatches(vec<vec<om<uint16_t>>>& omatch, vec<int64_t>& cinvp,
        const vec<vec<unsigned char>>& all_closuresp){
    // close omatches under inv
    auto omatch2 = omatch;
    #pragma omp parallel for
    for(int kk = 0; kk < omatch2.isize(); kk++){
        for(auto& o: omatch2[kk]){
            // add rcs
            auto ro = o;
            ro.c2 = cinvp[o.c2];
            ro.start1 = all_closuresp[kk].size()-o.start1-o.len;
            ro.start2 = all_closuresp[o.c2].size()-o.start2-o.len;
            #pragma omp critical
            { omatch[cinvp[kk]].push_back(ro); }
        }
    }
    Destroy(omatch2);
    // uniqsort
    #pragma omp parallel for
    for(int i = 0; i < omatch.isize(); i++)
        UniqueSort(omatch[i]);

    // validate
    int unalign = 0;
    #pragma omp parallel for
    for(int kk = 0; kk < omatch.isize(); kk++){
        for(auto& o: omatch[kk]){
            // add rcs
            for(int nn = 0; nn < o.len; nn++){
                if(all_closuresp[kk][o.start1+nn]!=
                        all_closuresp[o.c2][o.start2+nn]){
                    #pragma omp atomic
                    unalign = unalign+1;
                }
            }
        }
    }
    ForceAssertEq(unalign,0); 
}

template <class VBV>
void ExtremeAsm( VBV& bases, const vec<int64_t> & readlist, 
        digraphE<basevector>& Dl, vec<int>& dinvl, const int K, 
        const vec<String>& controls, const vec<nptuple>& Matches, 
        const Bool verbose, String out_dir, unsigned int TINKER){
    
    double clock1 = WallClockTime( );
    vec<int64_t> rdlist;
    rdlist.reserve(readlist.size());
    for(auto r:readlist)
        if(r!=0) // protection
            rdlist.push_back(r);
    auto rdlist2 = rdlist; 

    // add inverses
    rdlist.reserve(2*rdlist.size());
    for(const auto rd: rdlist2)
        rdlist.push_back(-rd);
    Destroy(rdlist2);

    // Damn!!! 
    ParallelUniqueSort(rdlist);
    // rewrite midpoint onwards to have correct cinv!
    int hs = rdlist.isize()/2;
    #pragma omp parallel for
    for(int i = 0; i < hs; i++)
        rdlist[i+hs] = -rdlist[i];

    // get the reads
    vecbasevector reads;
    reads.reserve(rdlist.size());
    for(auto rd: rdlist){
        if(rd>0)
            reads.push_back(bases[rd]);
        else
            reads.push_back(ReverseComplement(bases[-rd]));
    }

    // create map and translate matches
    std::unordered_map<int,int> readtoindex;
    for(int i = 0; i < rdlist.isize(); i++)
        readtoindex[rdlist[i]]=i;
    
    vec<vec<om<uint16_t>>> omatch;
    omatch.resize(rdlist.size());
    for(auto np: Matches){
        if(np.first.first!=0 && np.second.first!=0){
            int id1 = readtoindex[np.first.first];
            int id2 = readtoindex[np.second.first];
            omatch[id1].push_back(
                    om<uint16_t>(id2,np.first.second,
                    np.second.second,np.third));
        }
    } 
    cout << Date( ) << ": " << TimeSince(clock1) << " used in ExtremeAsm.1" << endl;
    
    ExtremeAsm(reads, Dl, dinvl, K, controls, omatch, verbose, out_dir, TINKER);
}

template void ExtremeAsm( VirtualMasterVec<basevector>& bases, 
        const vec<int64_t>& readlist, digraphE<basevector>& Dl, 
        vec<int>& dinvl, const int K, const vec<String>& controls,
        const vec<nptuple>& Matches, const Bool verbose, String out_dir,
        unsigned int TINKER);

template void ExtremeAsm( vecbasevector& bases, const vec<int64_t>& readlist,
        digraphE<basevector>& Dl, vec<int>& dinvl, const int K, 
        const vec<String>& controls, const vec<nptuple>& Matches, 
        const Bool verbose, String out_dir, unsigned int TINKER);

void ExtremeAsm( vecbasevector& reads, digraphE<basevector>& Dl, 
        vec<int>& dinvl, const int K, const vec<String>& controls, 
        vec<vec<om<uint16_t>>>& omatch, const Bool verbose, String out_dir, 
        unsigned int TINKER)
{
    double clock = WallClockTime();
    double clock2a = WallClockTime( );

    // ----------- Convert Reads to Closures --------------------
    if(verbose) PRINTDEETS("creating closures from reads");

    vec<vec<unsigned char>> all_closuresp;
    Reads2CharClosures(all_closuresp, reads);
    auto all_cl = all_closuresp;
    

    // ----------- Do extreme assembly --------------------------
    if(verbose) PRINTDEETS("creating super graph");

    // prepare for assembly building
    unsigned char numE = 4;
    vec<int> invp = {3,2,1,0}; // A = 0, C = 1, G = 2, T = 3
    
    // closure indexing
    vec<int64_t> cinvp; 
    vec<vec<int>> cip;
    InvertClosures(all_closuresp, cip, numE);

    int64_t N = all_closuresp.size();
    cinvp.resize(N,-1);
    for(int kk = 0; kk < N/2; kk++){
        cinvp[kk] = kk+N/2;
        cinvp[kk+N/2] = kk;
    }
    auto cinv2 = cinvp;
  
    SanitizeMatches(omatch,cinvp,all_closuresp);

    cout << Date( ) << ": " << TimeSince(clock2a)
          << " used in ExtremeAsm.2a" << endl;
    double clock2b = WallClockTime( );

    // create graph
    vec<int> dinvp;
    digraphE<vec<unsigned char>> Dp;
    Bool GET_META = True;
    vec<vec<quad<int,uint16_t,int,uint16_t>>> read_supp;
    NucleateGraph_HP(omatch, numE, all_closuresp, cinvp, cip, Dp, dinvp, 
            GET_META, read_supp, verbose,"",True);
    cout << Date( ) << ": " << TimeSince(clock2b)
          << " used in ExtremeAsm.2b" << endl;
    double clock2c = WallClockTime( );

    // ---------- Do any tinkering you like here -------------------
    // 1. unsigned TINKER is 0 by default.
    // 2. the graph edges are vec<uchar> where {A,C,G,T} maps to {0u,1u,2u,3u}.
    //    so that complementary bases add to 3 and allow easy RC. 
    // 3. Input: Dp, dinvp, all_closuresp (which are vec<vec<uchar>> reads), supports
    // 4. Peer at RemoveSimpleHangs to get an idea of how to use existing machinery
    //    to do edits. 
    // 5. Any advanced machinery from before may be built on assumptions
    //    and also may need templating to work with uchars.
    // 6. Note that the first half of reads and the second half of reads 
    //    are RCs of each other.
    if(TINKER>0){
        CallGraphEditors(Dp,dinvp,all_cl,cinv2,read_supp,verbose);
    }else{
        if(verbose) PRINTDEETS("clean-up graph");
        BabyZipper(Dp, dinvp, numE, 10, 10);

        vec<int> dels;
        RemoveSimpleHangs(Dp, dinvp,dels);
        ParallelUniqueSort(dels);
        Dp.DeleteEdgesParallel(dels);
        RemoveUnneededVertices(Dp,dinvp);
        CleanupCore(Dp,dinvp);

    }

    double clock3 = WallClockTime( );
    // ------------ translations ------------------------------------
    if(verbose) PRINTDEETS("translate graphs to known language");

    // convert edges from char to string to basevec
    vec<basevector> superedgeseq;
    CharEdges2BaseVec(Dp, superedgeseq);
    auto const& new_from_edge_obj = Dp.FromEdgeObj();
    auto const& new_to_edge_obj = Dp.ToEdgeObj();
    auto const& new_from = Dp.From();
    auto const& new_to = Dp.To();
    Dl.Initialize(new_from,new_to,superedgeseq,new_to_edge_obj,new_from_edge_obj);
    dinvl = dinvp;

    // save graph at out_dir
    if(out_dir!=""){
        if(verbose) PRINTDEETS("saving graphs to disk");
        Mkdir777(out_dir);
        BinaryWriter::writeFile(out_dir+"/"+"a.sup_loc",Dl);
        BinaryWriter::writeFile(out_dir+"/"+"a.sup_loc.inv",dinvl);
    }
    cout << Date( ) << ": " << TimeSince(clock3) << " used in ExtremeAsm.3" << endl;
    
    if(verbose)
        cout << Date() << ": took "<<TimeSince(clock)<<endl;
    
}

