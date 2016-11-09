///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include <atomic>
#include <unordered_set>

#include "paths/long/EvalAssembly.h"
#include "paths/long/RefTrace.h"
#include "paths/long/RefTraceControl.h"
#include "paths/long/RefTraceTools.h"

#include "Basevector.h"
#include "CoreTools.h"
#include "Equiv.h"
#include "ParallelVecUtilities.h"
#include "PrintAlignment.h"
#include "graph/Digraph.h"
#include "kmers/KmerRecord.h"
#include "math/Functions.h"
#include "math/HoInterval.h"
#include "pairwise_aligners/SmithWatAffine.h"
#include "pairwise_aligners/SmithWatBandedA.h"
#include "paths/HyperBasevector.h"
#include "paths/RemodelGapTools.h"
#include "paths/long/AssessBestAlignCore.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/MakeKmerStuff.h"
#include "polymorphism/Edit.h"
#include "reporting/PerfStat.h"
#include "paths/long/Variants.h"
#include "paths/long/ReadOriginTracker.h"
#include "paths/long/CreateGenome.h"
#include "paths/long/Logging.h"
#include "util/NullOStream.h"

#include "graph/DigraphTemplate.h"

namespace { // open anonymous namespace


int PenaltyFunc( const int e, const int g, const int d )
{    return 3*e + 100*g + d;    }

vec<pair<int,edit0>> AlignToEdits(const align& a, const basevector& aSeq,
        const basevector& rSeq, N50Calculator* pN50Calc )
{
    vec<pair<int,edit0>> edits;
    int p1 = a.pos1( ), p2 = a.pos2( );
    if ( pN50Calc && p2 ) pN50Calc->flush();
    Bool first = True;
    for ( int j = 0; j < a.Nblocks( ); j++ ) 
    {   if ( a.Gaps(j) > 0 )
        {    if ( !( first && edits.nonempty( )
                  && p2 <= edits.back( ).first ) )
             {    edits.push( p2, edit0( DELETION, a.Gaps(j) ) );    }
             if ( first && edits.nonempty( )
                  && p2 >= edits.back( ).first )
             {    first = False;    }
             p2 += a.Gaps(j);
             if ( pN50Calc ) pN50Calc->flush();    }
        if ( a.Gaps(j) < 0 ) 
        {    if ( !( first && edits.nonempty( )
                  && p2 <= edits.back( ).first ) )
             {    edits.push( p2, edit0( INSERTION, 
                       basevector( aSeq, p1, -a.Gaps(j) )
                       .ToString( ) ) );    }
             if ( first && edits.nonempty( )
                  && p2 >= edits.back( ).first )
             {    first = False;    }
             p1 -= a.Gaps(j);
             if ( pN50Calc ) pN50Calc->flush();    }
        unsigned matchCount = 0;
        for ( int x = 0; x < a.Lengths(j); x++ ) 
        {    if ( aSeq[p1] == rSeq[p2] )
                  matchCount += 1;
             else
             {    if ( !( first && edits.nonempty( )
                       && p2 <= edits.back( ).first ) )
                  {    edits.push( p2, edit0( SUBSTITUTION, 
                            (char) as_base(aSeq[p1]) ) );    }
                  if ( first && edits.nonempty( )
                       && p2 >= edits.back( ).first )
                  {    first = False;    }
                  if ( pN50Calc )
                  {    if ( matchCount )
                       {    pN50Calc->add(matchCount); matchCount = 0;    }
                       pN50Calc->flush();    }    }
             ++p1; ++p2;    }
        if ( pN50Calc && matchCount )
             pN50Calc->add(matchCount);    }
    return edits;
}

vec< std::tuple<int64_t,int64_t,int64_t,int64_t,int64_t> >
PerfectMatch(const align&a, const basevector&S,const basevector&T){
    vec< std::tuple<int64_t,int64_t,int64_t,int64_t,int64_t> > out;
    int p1 = a.pos1( ), p2 = a.pos2( );

    Bool first = True;
    for ( int j = 0; j < a.Nblocks( ); j++ ){
        if ( a.Gaps(j) > 0 ) p2 += a.Gaps(j);
        if ( a.Gaps(j) < 0 ) p1 -= a.Gaps(j);
        int match_size=0;
        for ( int x = 0; x < a.Lengths(j); x++ ){
            if ( S[p1] != T[p2] ){
                if(match_size>0){
                    out.emplace_back(match_size,p1-match_size,p1,p2-match_size,p2);
                }
                match_size=0;
            }
            else{
                ++match_size;
            }
            ++p1; ++p2;
        }
        if(match_size>0){
//std::cout << "match " << match_size << " "
//                     << p1-match_size << " "
//                     << p1 << " "
//                     << p2-match_size <<" "
//                     << p2 << std::endl;
            out.emplace_back(match_size,p1-match_size,p1,p2-match_size,p2);
        }
    }
    return out;
}


unsigned int SmithWatAffineWrapper(const basevector & S, const basevector & T,
                            alignment & a,
                            bool penalize_left_gap, bool penalize_right_gap,
                            const int mismatch_penalty=3, const int gap_open_penalty=12, const int gap_extend_penalty=1){
    static std::atomic<int> depth;
    unsigned int out;
//    std::cout<<"SmithWatAffineWrapper at depth " << depth << "looking at " << S.size() << " " << T.size() << std::endl;

    if(depth==0 || (depth < 4 && (S.size() > 40000 || T.size() > 40000)) ){
        ++depth;
//    std::cout<<"SmithWatAffineWrapper calling super" << std::endl;
        out = SmithWatAffineSuper(S,T,a,501,penalize_left_gap,penalize_right_gap,0,mismatch_penalty,gap_open_penalty,gap_extend_penalty,SmithWatAffineWrapper);
        --depth;
    }
    else{
//    std::cout<<"SmithWatAffineWrapper calling parallel2" << std::endl;
        out = SmithWatAffineParallel2(S,T,a,penalize_left_gap,penalize_right_gap,mismatch_penalty,gap_open_penalty,gap_extend_penalty);
    }
    return out;
};
// Find the best alignment between the assembly sequence and reference.
// Right now the slow SmithWatAffineParallel is the bottle neck.
align BestAlign(const basevector& b1, const basevector& b2, ostream& out, int iAlignerK)
{
    alignment al;
    double tclock = WallClockTime( );
    switch(iAlignerK){
    case -1:
        SmithWatAffineWrapper( b1, b2, al, False, False );
//        std::cout << "using SmithWatAffineWrapper ran for " << TimeSince(tclock) << std::endl;
        break;
    case 0:
//        SmithWatAffineParallel( b1, b2, al, False, False );
        SmithWatAffineParallel2( b1, b2, al, False, False );
//        std::cout << "using SmithWatAffineParallel2 ran for " << TimeSince(tclock) << std::endl;
        break;
    default:
        SmithWatAffineSuper( b1, b2, al, iAlignerK, False, False,0 );
        // std::cout << "SmithWatAffineSuper at K=" << iAlignerK << " ran for " << TimeSince(tclock)<< std::endl;
        break;
    }
    align a(al);
    out << "Alignment of best path sequences [" << a.pos1() << ", " << a.Pos1()
        << " (len=" << b1.size() << ")"
        << ") to reference [" << a.pos2() << ", " << a.Pos2() << ")"
        << " (len=" << b2.size() << ")"<< endl;;
    PrintVisualAlignmentClean(True, out, b1, b2, a);
    return a;
}

// Find the best alignment between the assembly sequence and a circular reference.
// as of r48407, SmithWatAffineSuper fails to align doubled reference so this is a work around.
// if SmithWatAffineSuper can actually align b1 to b2 (for sequences of the form b.append(b)), a straight up call would suffice
// basically, use the information in coors_hbp to find an underlying edge sequence to "split" the alignment problem
// right in the middle. This removes duplicated kmers when the reference is doubled, and will work
// as long as the underlying aligner works with the original, undoubled reference
align BestAlignForDoubled(const basevector& b1, const basevector& b2
                         , const vec<std::tuple<int64_t,int64_t,int,int64_t,int64_t,int64_t,int64_t>>& coors_hbp
                         , ostream& out)
{
//    std::cout << "in ForDoubled aligning " << b1.size() << " vs " << b2.size() << std::endl;
    int iAlignerK=501;

    int64_t b1_2 = std::get<5>(coors_hbp.front());
    int64_t e1_2 = std::get<6>(coors_hbp.back());
    int64_t max_range = b2.size();
    ForceAssert(max_range%2==0);
    int64_t range = max_range/2;

    align output;
    if( e1_2-b1_2 < range+1000){
//        std::cout << "ForDoubled: less than half-filled, shrinking" << std::endl;
        int64_t padding = (range-e1_2+b1_2)/2;

        int64_t b2_2 = max(int64_t(0),b1_2-padding);
        int64_t e2_2 = min(max_range,e1_2+padding);
        alignment al;
        SmithWatAffineWrapper( b1, basevector(b2.begin()+b2_2,b2.begin()+e2_2), al, False, False );
        output=align(al);
        output.Setpos2(output.pos2()+b2_2);
    }
    else{
//        std::cout << "ForDoubled: more than half-filled, splitting" << std::endl;
        vec< std::pair<int64_t,size_t> > distance_index; distance_index.reserve(coors_hbp.size());
        for(size_t cc=0;cc<coors_hbp.size();++cc){
            distance_index.emplace_back(  min( abs(std::get<5>(coors_hbp[cc])-range)
                                             , abs(std::get<6>(coors_hbp[cc])-range)
                                             ) , cc);

        }
        std::sort(distance_index.begin(),distance_index.end());
        bool bFoundSeed=false;
        int64_t b_seed_1=-1,e_seed_1=-1,b_seed_2=-1,e_seed_2=-1;
        for(size_t ee=0;!bFoundSeed&&ee<distance_index.size();++ee){
            int64_t be_1 = std::get<0>(coors_hbp[distance_index[ee].second]);
            int64_t ee_1 = std::get<1>(coors_hbp[distance_index[ee].second]);

            int64_t be_2 = std::get<5>(coors_hbp[distance_index[ee].second]);
            int64_t ee_2 = std::get<6>(coors_hbp[distance_index[ee].second]);
//std::cout << "finding seed with " << std::endl;
//PRINT4(be_1,ee_1,be_2,ee_2);

            int64_t padding=500;
            int64_t minmatch=3000;
            int64_t seed_trim=100;
            int64_t b2_2 = max(int64_t(0),be_2-padding);
            int64_t e2_2 = min(max_range,ee_2+padding);
            if( e2_2-b2_2 < 2*padding || be_1 > ee_1 ) continue;

//PRINT6(be_1,ee_1,be_2,ee_2,b2_2,e2_2);
            alignment al_sub;
            const basevector b1_sub(b1.begin()+be_1,b1.begin()+ee_1);
            const basevector b2_sub(b2.begin()+b2_2,b2.begin()+e2_2);
//            SmithWatAffineSuper( b1_sub, b2_sub, al_sub, iAlignerK, False, False,0 );
            SmithWatAffineWrapper( b1_sub, b2_sub, al_sub, False, False);
            if(    al_sub.StartOnQuery()==0 && al_sub.EndOnQuery()==ee_1-be_1){
//std::cout << "fully aligned" << std::endl;
                align a_sub(al_sub);
                vec< std::tuple<int64_t,int64_t,int64_t,int64_t,int64_t> > pms=PerfectMatch(al_sub, b1_sub,b2_sub);
                std::sort(pms.rbegin(),pms.rend());
//PRINT(pms.size());
                if(pms.size()>0&&std::get<0>(pms.front())>=minmatch){
//PRINT(std::get<0>(pms.front()));
                    b_seed_1=be_1+std::get<1>(pms.front())+seed_trim;
                    e_seed_1=be_1+std::get<2>(pms.front())-seed_trim;
                    b_seed_2=b2_2+std::get<3>(pms.front())+seed_trim;
                    e_seed_2=b2_2+std::get<4>(pms.front())-seed_trim;

//                    PRINT4(be_1,ee_1,be_2,ee_2);

//                    PRINT5(std::get<0>(pms.front())
//                          ,std::get<1>(pms.front())
//                          ,std::get<2>(pms.front())
//                          ,std::get<3>(pms.front())
//                          ,std::get<4>(pms.front()));
                    ForceAssert(std::equal(b1.begin()+b_seed_1
                                          ,b1.begin()+e_seed_1
                                          ,b2.begin()+b_seed_2));
                    bFoundSeed=true;
                }
            }
        }
        bool bAligned=false;
        if(bFoundSeed){
            int64_t bA_1=0;
            int64_t eA_1=b_seed_1;
            int64_t bA_2=0;
            int64_t eA_2=b_seed_2;
            avector<int> A_gaps, A_lengths;
            int A_pos1, A_Pos1, A_pos2, A_Pos2, A_errors;
            if( eA_1 > bA_1){
                alignment al_A;
                // reference has already been cut in half, if the following call fail (n^2 wise), such
                // failure belongs completely in the scope of SmithWatAffineSuper
                SmithWatAffineSuper( basevector(b1.begin()+bA_1,b1.begin()+eA_1)
                                   , basevector(b2.begin()+bA_2,b2.begin()+eA_2)
                                   , al_A, iAlignerK, False, true,0 );
                al_A.Unpack( A_pos1, A_pos2, A_errors, A_gaps, A_lengths );
                A_Pos1=al_A.Pos1();
                A_Pos2=al_A.Pos2();
            }
            else{
                A_pos1=bA_1;
                A_Pos1=eA_1;
                A_pos2=bA_2;
                A_Pos2=eA_2;
                A_errors=0;
                A_gaps.Append(0);
                A_lengths.Append(0);
            }

            int64_t bB_1=e_seed_1;
            int64_t eB_1=b1.size();
            int64_t bB_2=e_seed_2;
            int64_t eB_2=b2.size();
            alignment al_B;

            SmithWatAffineSuper( basevector(b1.begin()+bB_1,b1.begin()+eB_1)
                               , basevector(b2.begin()+bB_2,b2.begin()+eB_2)
                               , al_B, iAlignerK, true, False,0 );

            avector<int> B_gaps, B_lengths;
            int B_pos1, B_pos2, B_errors;
            al_B.Unpack( B_pos1, B_pos2, B_errors, B_gaps, B_lengths );

//std::cout << std::endl;
//            PRINT2(A_gaps.length,B_gaps.length);
//            PRINT2(A_lengths.length,B_lengths.length);
//            PRINT2(A_gaps(0),B_gaps(0));
//std::cout << std::endl;
//            PRINT3(bA_1,eA_1,eA_1-bA_1);
//            PRINT3(bA_2,eA_2,eA_2-bA_2);
//            PRINT3(bB_1,eB_1,eB_1-bB_1);
//            PRINT3(bB_2,eB_2,eB_2-bB_2);
//std::cout << std::endl;
//            PRINT4(A_pos1,A_Pos1,A_pos2,A_Pos2);
//            PRINT4(B_pos1,al_B.Pos1(),B_pos2,al_B.Pos2());
//std::cout << std::endl;

            if(   B_gaps.length > 0 && A_gaps.length > 0
               && B_gaps(0)==0 && B_pos1==0 && B_pos2==0
               && A_Pos1 == eA_1-bA_1 && A_Pos2 ==eA_2-bA_2
              ){
                A_lengths(A_lengths.length-1) += e_seed_1-b_seed_1 + B_lengths(0);
                for( size_t ii=1 ; ii<B_lengths.length ; ++ii){
                    A_lengths.Append(B_lengths(ii));
                    A_gaps.Append(B_gaps(ii));
                }
                output.Set(A_pos1,A_pos2,A_gaps,A_lengths,A_lengths.length);
                bAligned=true;
            }
            else{
                std::cout << "WARNING: BestAlignForDoubled couldn't merge splitted alignment." << std::endl;
            }
        }
        if(!bAligned){
            alignment al;
            SmithWatAffineWrapper( b1, b2, al, False, False );
            output=align(al);
        }
    }

//    double tclock = WallClockTime( );
//    std::cout << "using SmithWatAffineWrapper ran for " << TimeSince(tclock) << std::endl;
    out << "Alignment of best path sequences [" << output.pos1() << ", " << output.Pos1()
        << " (len=" << b1.size() << ")"
        << ") to reference [" << output.pos2() << ", " << output.Pos2() << ")"
        << " (len=" << b2.size() << ")"<< endl;;
    PrintVisualAlignmentClean(True, out, b1, b2, output);
    return output;
}


// The best path may contain gaps (non-continuous edge sequences). We will
// break the path into segments at the gaps and align the segments separately.
void BreakPath(const HyperBasevector& hbp, const vec<int>& path, const vec<int>& eids, const vec<std::pair<int,int>>& limits
                                         , vec<vec<int>>& path_segs, vec<vec<int>>& eid_segs, vec<vec<std::pair<int,int>>>& limits_segs)
{
    vec<int> to_left, to_right;
    hbp.ToLeft(to_left); hbp.ToRight(to_right);
    for (size_t i = 0; i < path.size(); i++) {
        size_t j = i+1;
        while (j < path.size() && to_left[path[j]] == to_right[path[j-1]]) 
            j++;
        path_segs.push_back(vec<int>(path.begin()+i, path.begin()+j));
        eid_segs.push_back(vec<int>(eids.begin()+i, eids.begin()+j));
        limits_segs.emplace_back(limits.begin()+i, limits.begin()+j);
        i = j-1;
    }
}

// Print the path. Convert the edge id from hbp to original assembly graph.
void PrintPath(const vec<int>& path, const vec<pair<int,Bool>>& hbp_to_hb, 
        ostream& out) 
{
    for (size_t i = 0; i < path.size(); i++) {
        if (i!=0) out << " ";
        out << hbp_to_hb[path[i]].first;
    }
    out << endl;
}

// From the alignment, accumulate assembly error stats.
void AddErrors(const align& a, const basevector& aSeq, const basevector& rSeq,
                    N50Calculator* pN50Calc, AssemblyError::err_per_ref& res)
{
    for(const auto& edit: AlignToEdits(a,aSeq,rSeq,pN50Calc))
    {
        switch (edit.second.etype) {
            case INSERTION:
                res.AddIndel(edit.first,+int64_t(edit.second.seq.size()));
                break;
            case DELETION:
                res.AddIndel(edit.first,-int64_t(edit.second.n));
                break;
            case SUBSTITUTION:
                res.AddSub(edit.first,1);
                break;
        }
    }
}

class circular_tester
{
public:
    circular_tester(const HyperBasevector&in):hbp(in),bCircular(false){
        vec<vec<int>> c_e;
        hbp.ComponentsE(c_e);
        edge_component.resize(hbp.EdgeObjectCount());
        for(size_t cc=0;cc<c_e.size();++cc){
            for(const auto& entry: c_e[cc]){
                edge_component[entry]=cc;
            }
        }
        hbp.ToRight( to_right );
    };
    void LoadSeq(const vec<int>&in){ seq=in; EvalSeq();};
    void UpdateSeq(bool bPopFront, bool bPushBack, int right=std::numeric_limits<int>::min()){
        if(bPopFront){
            ForceAssert(seq.size()>0);
            for(size_t ii=0;ii+1<seq.size();++ii){ seq[ii]=seq[ii+1]; }
            seq.resize(seq.size()-1);
        }
        if(bPushBack){ seq.push_back(right); }
        if(bPopFront || bPushBack){
//            std::cout <<"circular_tester::UpdateSeq " << (int)bPopFront << " " << (int)bPushBack << " " << right << std::endl;
            EvalSeq();
        }
    };
    void EvalSeq(){
        bCircular=true;

        for(size_t ss=0; bCircular && ss<seq.size(); ++ss){
            int curr_edge = seq[ss];
            if( curr_edge < 0 ) continue;

            int next_edge=-1;
            size_t shift=1;
            for(; next_edge < 0 && shift<seq.size();++shift){
                next_edge = seq[(ss+shift)%seq.size()];
            }
            if( next_edge >= 0){
                shift-=1;
                if(shift==1 || edge_component[curr_edge] == edge_component[next_edge]){
                    int vv = to_right[curr_edge];
                    bool bFound=false;
                    for(int jj=0;!bFound && jj<hbp.FromSize(vv);++jj){
                        bFound = next_edge == hbp.EdgeObjectIndexByIndexFrom(vv,jj) ;
                    }
//std::cout << ss << "+" << shift << "|" << curr_edge << " " << next_edge << "|"
//          << int(bFound) <<  " " << int(ss+1==seq.size()) << " " << int(shift) << " " << int(curr_edge==next_edge) << std::endl;
                    bCircular=bFound || (ss+1==seq.size() && shift==1 && curr_edge==next_edge) ;
                }
            }
        }
        /*
        if(bCircular){
            std::cout << "circular: ";
        }
        else{
            std::cout << "non-circular: ";
        }

            auto oss = ostream_iterator<int>(std::cout, " ");
            std::copy(seq.begin(),seq.end(),oss);
            std::cout << std::endl;
            */
    };
    bool IsCircular()const{return bCircular;};
    const vec<int>& getSeq()const{return seq;};
private:
    circular_tester();
    const HyperBasevector& hbp;
    vec<int> seq;
    vec<size_t> edge_component;
    vec<int> to_right;
    bool bCircular;
};

void MakeHBPDot(ostream& out,const HyperBasevector&hbp,const vec<pair<int,Bool>>& hbp_to_hb){
    vec<String> vl(hbp.N());
    vec<vec<String>> el(hbp.N());
    for(int64_t vv=0;vv<hbp.N();++vv){
        vl[vv]=ToString(vv);
        el[vv].resize(hbp.FromSize(vv));
        for(int64_t ee=0;ee<hbp.FromSize(vv);++ee){
            int eid = hbp.EdgeObjectIndexByIndexFrom(vv,ee);
            el[vv][ee]=ToString(eid) + " " + ToString( hbp_to_hb[ eid ].first);
        }
    }
    hbp.DOT_vl(out,vl,"",vec<vec<String>>(),vec<String>(),el);
}


}// close anonymous namespace

void AssemblyError::PrintSummary(ostream& out)
{
    auto gaps=GetGaps();
    auto indels=GetIndels();
    const int subs=GetSubs();
    const auto viLGap=GetLeftGaps();
    const auto viRGap=GetRightGaps();

    const bool bGoodTopology=std::all_of(this->begin(),this->end()
                                        ,[](const std::pair<int,err_per_ref> & in ){return in.second.GoodTopology();} );

    Sort(gaps);
    Sort(indels);
    out << "Best paths:\n";
    for(const auto& bp: GetBestPaths()){
        std::copy(bp.begin(),bp.end(),ostream_iterator<int>(out,","));
        out << std::endl;
    }

    out << "EVAL_ASSEMBLY: ";

    bool bHasLGap = std::any_of(viLGap.begin(),viLGap.end(),[](int in){return in!=0;});
    bool bHasRGap = std::any_of(viRGap.begin(),viRGap.end(),[](int in){return in!=0;});

    if (gaps.empty() && indels.empty() && subs == 0 && !bHasLGap && !bHasRGap ){
        out << "perfect" << endl;
        if(!bGoodTopology){
            out << "circular-topology problem" << std::endl;
        }
        return;
    }

    auto oss = ostream_iterator<int>(out, " ");
    bool prev = false;

    if(bHasLGap){
        out << "left gaps=";
        copy(viLGap.begin(), viLGap.end(), oss);
        prev = true;
    }
    if(bHasRGap){
        out << "right gaps=";
        copy(viRGap.begin(), viRGap.end(), oss);
        prev = true;
    }

    if (!gaps.empty()) {
        out << "gaps=";
        copy(gaps.begin(), gaps.end(), oss);
        prev = true;
    }
    if (!indels.empty()) {
        if (prev) out << ", ";
        out << "indels=";
        copy(indels.begin(), indels.end(), oss);
        prev = true;
    } 
    if (subs > 0) {
        if (prev) out << ", ";
        out << "nsubs=" << subs;
        prev = true;
    }
    out << endl;
    if(!bGoodTopology){
        out << "circular-topology problem" << std::endl;
    }
}

namespace{
    void FixInversion( const HyperBasevector& hb, vec<int>& inv2 )
    {    inv2.resize( hb.EdgeObjectCount( ), -1 );
         vec< triple<basevector,int,int> > uni2;
         vec<Bool> used;
         hb.Used(used);
         for ( int e = 0; e < hb.EdgeObjectCount( ); e++ )
         {    if ( !used[e] ) continue;
              uni2.push( hb.EdgeObject(e), 0, e );
              basevector b = hb.EdgeObject(e);
              b.ReverseComplement( );
              uni2.push( b, 1, e );    }
         Sort(uni2);
         for ( int i = 0; i < uni2.isize( ); i++ )
         {    int j;
              for ( j = i + 1; j < uni2.isize( ); j++ )
                   if ( uni2[j].first != uni2[i].first ) break;
              if ( j - i == 2 && uni2[i].second != uni2[i+1].second )
              {    inv2[ uni2[i].third ] = uni2[i+1].third;
                   inv2[ uni2[i+1].third ] = uni2[i].third;    }
              else if ( j - i == 1 ) inv2[ uni2[i].third ] = -1;
              // else PRINT(j-i);
              i = j - 1;    }    };
    void UnrollAssemblyLoops(HyperBasevector& out, const HyperBasevector& in, vec<std::pair<int,int>>& out_to_in){
        out=in;
        out_to_in.clear_and_resize(out.EdgeObjectCount());
        for(int ee=0;ee<out_to_in.isize();++ee){ out_to_in[ee]=std::make_pair(ee,1); }
        const int threshold= 3*out.K()-2;
        for(int vv=0;vv<out.N();++vv){
            const auto& from=out.From(vv);
            vec<basevector> loops;
            vec<int> loops_idx;
            for(size_t jj=0;jj<from.size();++jj){
                if(from[jj]==vv){
                    loops.push_back( out.EdgeObjectByIndexFrom(vv,jj));
                    loops_idx.push_back( out.EdgeObjectIndexByIndexFrom(vv,jj));
                }
            }
            vec<basevector> new_loops;
            vec<std::pair<int,int>> new_loops_info;
            for(size_t ll=0;ll<loops.size();++ll){
                basevector loop=TrimCat(out.K(),loops[ll],loops[ll]);
                for(int mul=2
                   ; loop.isize() < threshold && !Member(new_loops,loop)
                   ; ++mul, loop=TrimCat(out.K(),loop,loops[ll])){
                    new_loops.push_back(loop);
                    new_loops_info.emplace_back(loops_idx[ll],mul);
                }
            }
            for(size_t ll=0;ll<new_loops.size();++ll){
                if(!Member(loops,new_loops[ll])){
                    out_to_in.emplace_back(new_loops_info[ll].first,new_loops_info[ll].second);
                    out.AddEdge(vv,vv,new_loops[ll]);
                }
            }
        }
    };

class GraphTransformer{
    template <typename ELEM_T>
    class edge{
    public:
        edge(const edge&L,const edge&R,int K):seq(L.GetSeq()),len(L.GetLen()+R.GetLen()+1-K){ seq.append(R.GetSeq()); };
        edge():seq(),len(0){};
        edge(const vec<ELEM_T>& v,int64_t l ):seq(v),len(l){};
        edge(const ELEM_T& e,int64_t l ):seq({e}),len(l){};
        friend bool operator==(const edge& L, const edge& R){
            return L.GetSeq().size()==R.GetSeq().size() && std::equal(L.GetSeq().begin(),L.GetSeq().end(),R.GetSeq().begin());
        }
        int64_t& GetLen(){return len;};
        const int64_t& GetLen()const{return len;};
        vec<ELEM_T>& GetSeq(){return seq;};
        const vec<ELEM_T>& GetSeq()const{return seq;};
    private:
        vec<ELEM_T> seq; //a sequence of original graph edge
        int64_t len; // length in units of bp
    };
public:

    //an element of the symbolic graph edge
    class elem_t{
    public:
        // edge segment flag, note that FULL==L M R, where L/M/R denotes left/middle/right segments
        enum SEG_ENUM{SEG_FULL,SEG_LM,SEG_MR,SEG_L,SEG_M,SEG_R,SEG_ERR};
        explicit elem_t( int e):eid(e),seg(SEG_ENUM::SEG_FULL){};
        elem_t( int e, SEG_ENUM s):eid(e),seg(s){};
        int& GetEid(){return eid;};
        const int& GetEid()const{return eid;};
        SEG_ENUM& GetSeg(){return seg;};
        const SEG_ENUM& GetSeg()const{return seg;};
        friend bool operator==(const elem_t& L, const elem_t& R){ return L.eid==R.eid && L.seg==R.seg; }
    private:
        int eid; //edge id
        SEG_ENUM seg; //segment
    };
    static String ToString(const elem_t& elem){
        String out = ::ToString(elem.GetEid());
        switch( elem.GetSeg() ){
            case elem_t::SEG_FULL: out+=":F"; break;
            case elem_t::SEG_LM: out+=":LM"; break;
            case elem_t::SEG_MR: out+=":MR"; break;
            case elem_t::SEG_L: out+=":L"; break;
            case elem_t::SEG_M: out+=":M"; break;
            case elem_t::SEG_R: out+=":R"; break;
            default: out+=":!!!!!!"; break;
        }
        return out;
    };
    typedef edge<elem_t> edge_t;
    typedef digraphE< edge_t > graph_t;

    GraphTransformer(const HyperBasevector& org):m_org(org)
                                                ,m_new(org,vec<edge_t>(org.EdgeObjectCount()) ){
        for(int ee=0;ee<m_org.EdgeObjectCount();++ee){
            m_new.EdgeObjectMutable(ee).GetSeq()={elem_t(ee)};
            m_new.EdgeObjectMutable(ee).GetLen()=m_org.EdgeObject(ee).size();
        }
    };
    /* for all edges neighboring some edge shorter than 'short_neighbor', break it in to 2 or 3 pieces
     * if the segment can be at least short_segment
     * This is a performance tweak to yield shorter edges (of short_segment length) that can be aligned
     * much more easily. There for short_segment length increases with specificity but decreases with
     * down stream performance.
     */
    void BreakLongEdges(int64_t short_neighbor, int64_t short_segment);

    /* if there are at least one short edge going from v->w and at least another going w->v
     * expand the graph topology such that there are 4 vertices with no loop back v1->w1 and w2->v2
     * this alows AbsorbShortEdgges to wipe out short edges
     */
    void BreakTwoVerticesLoopBack(int64_t limit,size_t max_step=100){ for( size_t ss=0; BreakTwoVerticesLoopBackStep(limit) && ss<max_step; ++ss) {}; }

    /* look for edges shorter than 'limit' between v->w, and also no other edges going w->w, v->v, w->v
     * sort them in ascending order of shortest neighboring edges to vertices-z (z->v and w->z)
     * absorb such short edges such that z->v->w becomes z->w, or v->w->z becomes v->z
     * the original vertices v and w might get deleted, depends on if there are residue edges linking them
     */
    void AbsorbShortEdges(int64_t limit);


    /* look for edges E v->v and unroll them to be e^N, until the length is longer than 'limit'
     */
    bool UnrollShortLoops(int64_t limit);

    const graph_t& getNewGraph()const{return m_new;};


    // from the symbolic transformation graph, construct a list of edge objects containing the actual bp
    void getNewEdges(vec<basevector>& out, int64_t short_segment)const{
        out.clear();
        out.resize(m_new.EdgeObjectCount());
        for(size_t ee=0;ee<out.size();++ee){
            const vec<elem_t>& edge_seq = m_new.EdgeObject(ee).GetSeq();
            out[ee]=getSequence(edge_seq.front(),short_segment);
            for(size_t ss=1;ss<edge_seq.size();++ss){
                out[ee]=TrimCat(m_org.K(),out[ee],getSequence(edge_seq[ss],short_segment));
            }
        }
    }

    // return the bp representation of the graph edge, according to L-M-R breakdown with L/R being of length short_segment
    basevector getSequence(elem_t elem, int64_t short_segment)const
    {
        basevector::const_iterator ib = m_org.EdgeObject(elem.GetEid()).begin();
        basevector::const_iterator ie = m_org.EdgeObject(elem.GetEid()).end();
        switch(elem.GetSeg()){
            case elem_t::SEG_FULL :
                break;
            case elem_t::SEG_L :
                ie=ib+short_segment;
                break;
            case elem_t::SEG_M :
                ib += short_segment - m_org.K() + 1;
                ie -= short_segment - m_org.K() + 1;
                break;
            case elem_t::SEG_R :
                ib=ie-short_segment;
                break;
            case elem_t::SEG_LM :
                ie -= short_segment - m_org.K() + 1;
                break;
            case elem_t::SEG_MR :
                ib += short_segment - m_org.K() + 1;
                break;
            default:
                FatalErr("inconsistent graph");

        };
//        std::cout << m_org.EdgeObject(elem.GetEid()).size() << " "
//                  << std::distance(ib,ie) << " "
//                  << elem.GetEid() << " "
//                  << elem.GetSeg() << " "
//                  << std::endl;
        return basevector(ib,ie);

    };

    vec<pair<int,Bool>> new_to_old( const vec<pair<int,Bool>>& new_seq){
        vec<pair<int,Bool>> old_seq; old_seq.reserve(new_seq.size());
        int last_eid=-1;
        Bool bLastFlag=True;
        elem_t::SEG_ENUM last_seg=elem_t::SEG_ERR;

        for(const auto&entry: new_seq){ // for each new edge in new coordinate
            auto new_edge_seq = m_new.EdgeObject( entry.first ).GetSeq(); // the sequences of old edges in old coordinate
            const Bool bFlag = entry.second; // if the new edge is FW or not

            if( bFlag == False){ // reverse edges if new edge is RC'ed
                std::reverse(new_edge_seq.begin(),new_edge_seq.end());
            }

            for( const elem_t& elem: new_edge_seq){
                int eid = elem.GetEid();
                elem_t::SEG_ENUM seg = elem.GetSeg();
                if ( last_eid >=0 && eid != last_eid){ // if there's an unfinished edge, as we look at new edge
                    old_seq.emplace_back(last_eid,bLastFlag);
                    last_eid=-1;
                }
                if(   seg == elem_t::SEG_FULL
                   || (bFlag  && ( seg == elem_t::SEG_R || seg == elem_t::SEG_MR) ) //FW edges end with R segment
                   || (!bFlag && ( seg == elem_t::SEG_L || seg == elem_t::SEG_LM) ) //RC edges end with L segment
                  ){ // if an edge is done
                    old_seq.emplace_back(eid,bFlag);
                    last_eid=-1;
                }
                else{ // if an edge is not done, log it
                    last_eid=eid;
                    bLastFlag=bFlag;
                }
            }
        }
        if(last_eid>=0){
            old_seq.emplace_back(last_eid,bLastFlag);
        }
        return old_seq;
    }
    void MakeDot(ostream& out){
        vec<String> vl(m_new.N());
        vec<vec<String>> el(m_new.N());
        for(int64_t vv=0;vv<m_new.N();++vv){
            vl[vv]=::ToString(vv);
            el[vv].resize(m_new.FromSize(vv));
            for(int64_t ee=0;ee<m_new.FromSize(vv);++ee){
                int eid = m_new.EdgeObjectIndexByIndexFrom(vv,ee);
                el[vv][ee]=::ToString(eid)+":";
                for(const auto& ss: m_new.EdgeObject(eid).GetSeq()){
                    el[vv][ee]+=" " + ToString(ss);
                }
            }
        }
        m_new.DOT_vl(out,vl,"",vec<vec<String>>(),vec<String>(),el);
    }
private:
    bool BreakTwoVerticesLoopBackStep(int64_t limit);
    GraphTransformer();
    const HyperBasevector& m_org;
    graph_t m_new;
};
void GraphTransformer::BreakLongEdges(const int64_t short_neighbor, const int64_t short_segment){
    if( short_neighbor >= short_segment ){
        std::cout << "WARNING: GraphTransformer::BreakLongEdges() is not being used properly" << std::endl;
        return;
    }

    const int64_t min_length_to_work = 3*short_segment-2*(m_org.K()-1) ;


    vec<int> to_left, to_right;
    m_new.ToLeft(to_left), m_new.ToRight(to_right);
    vec<int> e_delete;
    m_new.DeleteEdges(e_delete);

    for(int eid=0,nOrgEdges=m_new.EdgeObjectCount();eid<nOrgEdges;++eid){
        const int64_t edge_len = m_new.EdgeObject(eid).GetLen();
        if( edge_len < min_length_to_work ) continue;

        const int vv = to_left[eid];
        const int ww = to_right[eid];

        bool bBreakLeft=false;
        for(int jj=0;!bBreakLeft && jj<m_new.ToSize(vv);++jj){
            bBreakLeft = m_new.EdgeObjectByIndexTo(vv,jj).GetLen() <= short_neighbor ;
        }

        bool bBreakRight=false;
        for(int jj=0;!bBreakRight && jj<m_new.FromSize(ww);++jj){
            bBreakRight = m_new.EdgeObjectByIndexFrom(ww,jj).GetLen() <= short_neighbor ;
        }
        if( bBreakLeft && bBreakRight){
            vec<elem_t> mid_seq = m_new.EdgeObject(eid).GetSeq();
            int64_t len_l = short_segment;
            int64_t len_r = short_segment;
            int64_t len_m = m_new.EdgeObject(eid).GetLen() + 2 * (m_org.K() - 1) -len_l - len_r;
            ForceAssert(len_m >= short_segment);

            if( mid_seq.front().GetSeg() ==elem_t::SEG_FULL){
                mid_seq.front().GetSeg() = elem_t::SEG_MR;
            }
            else if( mid_seq.front().GetSeg() ==elem_t::SEG_LM){
                mid_seq.front().GetSeg() = elem_t::SEG_M;
            }
            else{
                continue;
            }

            if( mid_seq.back().GetSeg() ==elem_t::SEG_FULL){
                mid_seq.back().GetSeg() = elem_t::SEG_LM;
            }
            else if( mid_seq.back().GetSeg() ==elem_t::SEG_MR){
                mid_seq.back().GetSeg() = elem_t::SEG_M;
            }
            else{
                continue;
            }
            int v_inner = m_new.N();
            int w_inner = v_inner+1;
            m_new.AddVertices(2);
            m_new.AddEdge(vv,v_inner,edge_t(elem_t(mid_seq.front().GetEid(),elem_t::SEG_L),len_l));
            m_new.AddEdge(v_inner,w_inner,edge_t(mid_seq,len_m));
            m_new.AddEdge(w_inner,ww,edge_t(elem_t(mid_seq.back().GetEid(),elem_t::SEG_R),len_r));
            e_delete.push_back(eid);
        }
        else if( bBreakLeft){
            int nv = m_new.N();
            vec<elem_t> right_seq = m_new.EdgeObject(eid).GetSeq();
            int64_t len_l = short_segment;
            int64_t len_r = m_new.EdgeObject(eid).GetLen() + m_org.K() - 1 -len_l;
            ForceAssert(len_r >= short_segment);

            if( right_seq.front().GetSeg() ==elem_t::SEG_FULL){
                right_seq.front().GetSeg() = elem_t::SEG_MR;
            }
            else if( right_seq.front().GetSeg() ==elem_t::SEG_LM){
                right_seq.front().GetSeg() = elem_t::SEG_M;
            }
            else{
                continue;
            }
            m_new.AddVertices(1);
            m_new.AddEdge(vv,nv,edge_t(elem_t(right_seq.front().GetEid(),elem_t::SEG_L),len_l));
            m_new.AddEdge(nv,ww,edge_t(right_seq,len_r));
            e_delete.push_back(eid);
        }
        else if( bBreakRight){
            int nv = m_new.N();
            vec<elem_t> left_seq = m_new.EdgeObject(eid).GetSeq();
            int64_t len_r = short_segment;
            int64_t len_l = m_new.EdgeObject(eid).GetLen() + m_org.K() - 1 -len_r;
            ForceAssert(len_l >= short_segment);

            if( left_seq.back().GetSeg() ==elem_t::SEG_FULL){
                left_seq.back().GetSeg() = elem_t::SEG_LM;
            }
            else if( left_seq.back().GetSeg() ==elem_t::SEG_MR){
                left_seq.back().GetSeg() = elem_t::SEG_M;
            }
            else{
                continue;
            }
            m_new.AddVertices(1);
            m_new.AddEdge(vv,nv,edge_t(left_seq,len_l));
            m_new.AddEdge(nv,ww,edge_t(elem_t(left_seq.back().GetEid(),elem_t::SEG_R),len_r));
            e_delete.push_back(eid);
        }
    }
    m_new.DeleteEdges(e_delete);
    m_new.RemoveDeadEdgeObjects();
}

bool GraphTransformer::UnrollShortLoops(int64_t limit){
    bool bChanged=false;
    for(int vv=0;vv<m_new.N();++vv){
        const auto& from=m_new.From(vv);
        vec<edge_t> loops;
        for(size_t jj=0;jj<from.size();++jj){
            if(from[jj]==vv){
                loops.push_back( m_new.EdgeObjectByIndexFrom(vv,jj));
            }
        }
        vec<edge_t> new_loops;
        for(size_t ll=0;ll<loops.size();++ll){
            for( edge_t loc_edge(loops[ll],loops[ll],m_org.K())
               ; loc_edge.GetLen() < limit
               ; loc_edge=edge_t(loc_edge,loops[ll],m_org.K())
               ){
                if( !Member(new_loops,loc_edge))  new_loops.push_back(loc_edge);
            }
        }
        for(size_t ll=0;ll<new_loops.size();++ll){
            if(!Member(loops,new_loops[ll])){
                m_new.AddEdge(vv,vv,new_loops[ll]);
                bChanged=true;
            }
        }
    }
    return bChanged;
}
bool GraphTransformer::BreakTwoVerticesLoopBackStep(int64_t limit){
    bool bChanged=false;

    vec< std::pair<int64_t,int> > candidates;
    candidates.reserve(m_new.EdgeObjectCount());

    vec<int> to_left, to_right;
    m_new.ToLeft(to_left), m_new.ToRight(to_right);
    auto HasLoopBack = [&to_left,&to_right,this](int ee){
        const int vv = to_left[ee];
        const int ww = to_right[ee];
        return vv!=ww && std::any_of(m_new.From(ww).begin(),m_new.From(ww).end(),[&vv](int nv){return nv==vv;});
    };
    //there can be dead edges
    for(int vv=0;vv<m_new.N();++vv){
        for(int jj=0;jj<m_new.FromSize(vv);++jj){
            auto ee = m_new.EdgeObjectIndexByIndexFrom(vv,jj);
            auto len = m_new.EdgeObject(ee).GetLen();
            if( len < limit && HasLoopBack(ee)){
                candidates.emplace_back(len,ee);
            }
        }
    }
    std::sort(candidates.begin(),candidates.end());
    std::unordered_set<int> frozen_vertices;
    vec<int> e_delete;
    for(const auto& candidate:candidates){
        const int ee = candidate.second;
        const int vv = to_left[ee];
        const int ww = to_right[ee];

        if(frozen_vertices.find(vv)!=frozen_vertices.end() || frozen_vertices.find(ww)!=frozen_vertices.end()) continue;
        vec<int> v_in,v_out,v_loop, w_in,w_out,w_loop,fwd,bck;
        vec<int> v_candidates{vv,ww};
        for(int jj=0;jj<m_new.FromSize(vv);++jj){
            const int eid = m_new.EdgeObjectIndexByIndexFrom(vv,jj);
            const int nv=m_new.From(vv)[jj];
            v_candidates.push_back(nv);
            if(nv==vv)     { v_loop.push_back(eid); }
            else if(nv==ww){ fwd.push_back(eid); }
            else           { v_out.push_back(eid); }
        }
        for(int jj=0;jj<m_new.ToSize(vv);++jj){
            const int eid = m_new.EdgeObjectIndexByIndexTo(vv,jj);
            const int nv = m_new.To(vv)[jj];
            v_candidates.push_back(nv);
            if( nv!=vv && nv!=ww) v_in.push_back(eid);
        }

        for(int jj=0;jj<m_new.FromSize(ww);++jj){
            const int eid = m_new.EdgeObjectIndexByIndexFrom(ww,jj);
            const int nv = m_new.From(ww)[jj];
            v_candidates.push_back(nv);
            if( nv == vv)  { bck.push_back(eid); }
            else if(nv==ww){ w_loop.push_back(eid); }
            else           { w_out.push_back(eid); }
        }
        for(int jj=0;jj<m_new.ToSize(ww);++jj){
            const int eid = m_new.EdgeObjectIndexByIndexTo(ww,jj);
            const int nv = m_new.To(ww)[jj];
            v_candidates.push_back(nv);
            if( nv!=vv && nv!= ww) w_in.push_back(eid);
        }

        if(std::any_of(v_candidates.begin(),v_candidates.end()
                      ,[&](int vv){return frozen_vertices.find(vv)!=frozen_vertices.end();})){
            continue;
        }
        if(! ( (v_in.size()>0 && w_out.size() >0) || (v_out.size()>0 && w_out.size()>0) )  ) continue;
        for(const auto& nv:v_candidates){ frozen_vertices.insert(nv); }
        if( v_out.size()>0 && w_out.size() >0 ){
            for(const auto bid: bck){
                for(const auto fid: fwd){
                    m_new.AddEdge(ww,ww,edge_t( m_new.EdgeObject(bid)
                                              , m_new.EdgeObject(fid)
                                              , m_org.K()));
                }
            }

            const int vv_new = m_new.N();
            const int ww_new = m_new.N()+1;
            m_new.AddVertices(2);

            for(const auto eid: v_in) { m_new.AddEdge(to_left[eid],vv_new,m_new.EdgeObject(eid)); }
            for(const auto eid: v_out){ m_new.AddEdge(vv_new,to_right[eid],m_new.EdgeObject(eid)); }
            for(const auto eid: v_loop){ m_new.AddEdge(vv_new,vv_new,m_new.EdgeObject(eid)); }

            for(const auto eid: w_in) { m_new.AddEdge(to_left[eid],ww_new,m_new.EdgeObject(eid)); }
            for(const auto eid: w_out){ m_new.AddEdge(ww_new,to_right[eid],m_new.EdgeObject(eid)); }
            for(const auto eid: w_loop){ m_new.AddEdge(ww_new,ww_new,m_new.EdgeObject(eid)); }

            for(const auto eid: bck){ m_new.AddEdge(ww_new,vv_new,m_new.EdgeObject(eid)); }

            for(const auto bid: bck){
                for(const auto fid: fwd){
                    m_new.AddEdge(vv_new,vv_new,edge_t( m_new.EdgeObject(fid)
                                                      , m_new.EdgeObject(bid)
                                                      , m_org.K()));
                }
            }
            e_delete.append(bck);
//            m_new.DeleteEdges(bck);
        }
        else if( v_out.size() > 0 ){
            const int nv=m_new.N();
            m_new.AddVertices(1);
            for(const auto fid: fwd){
                m_new.AddEdge(vv,nv,m_new.EdgeObject(fid));
                for(const auto bid: bck){
                    m_new.AddEdge(vv,vv,edge_t( m_new.EdgeObject(fid)
                                              , m_new.EdgeObject(bid)
                                              , m_org.K()));
                }
            }
            e_delete.append(fwd);
//            m_new.DeleteEdges(fwd);
        }
        else if( w_out.size() > 0){
            const int nv=m_new.N();
            m_new.AddVertices(1);
            for(const auto bid: bck){
                m_new.AddEdge(ww,nv,m_new.EdgeObject(bid));
                for(const auto fid: fwd){
                    m_new.AddEdge(ww,ww,edge_t( m_new.EdgeObject(bid)
                                              , m_new.EdgeObject(fid)
                                              , m_org.K()));
                }
            }
            e_delete.append(bck);
//            m_new.DeleteEdges(bck);
        }
        else{
std::cout<<"WARNING: BreakTwoVerticesLoopBackStep() encountered problems. retreating" << std::endl;
            break;
        }
        bChanged=true;
    }
    if(bChanged){
        m_new.DeleteEdges(e_delete);
        m_new.RemoveDeadEdgeObjects();
        m_new.RemoveEdgelessVertices();
    }

    return bChanged;
}

void GraphTransformer::AbsorbShortEdges(int64_t limit){
    vec<int> to_left, to_right;
    to_left.reserve(m_new.EdgeObjectCount()*2);
    to_right.reserve(m_new.EdgeObjectCount()*2);
    m_new.ToLeft(to_left), m_new.ToRight(to_right);
    auto GetLongestNeighbor = [&to_left,&to_right,this](int ee){
        std::pair<int64_t,int> len_dir(-1,-1);
        const int vv = to_left[ee];
        const int ww = to_right[ee];
        bool bHasLoop =   std::any_of(m_new.From(ww).begin(),m_new.From(ww).end(),[&](int nv){return nv==ww || nv==vv;})
                       && std::any_of(m_new.To(vv).begin(),m_new.To(vv).end(),[&](int nv){return nv==ww || nv==vv;}) ;
        //do not do anything if there are loop backs between v and w
        if(!bHasLoop){
            int64_t max_neighbor_len_R=-1;
            for(int jj=0;jj<m_new.FromSize(ww);++jj){
                max_neighbor_len_R = std::max(m_new.EdgeObjectByIndexFrom(ww,jj).GetLen(),max_neighbor_len_R);
            }
            int64_t max_neighbor_len_L=-1;
            for(int jj=0;jj<m_new.ToSize(vv);++jj){
                max_neighbor_len_L = std::max(m_new.EdgeObjectByIndexTo(vv,jj).GetLen(),max_neighbor_len_L);
            }
            if( max_neighbor_len_R >0 && max_neighbor_len_L > 0){
                if( max_neighbor_len_R > max_neighbor_len_L){
                    len_dir=std::make_pair(max_neighbor_len_L,0);
                }
                else{
                    len_dir=std::make_pair(max_neighbor_len_R,1);
                }
            }
            else if (max_neighbor_len_R>0){
                len_dir=std::make_pair(max_neighbor_len_R,1);
            }
            else if (max_neighbor_len_L>0){
                len_dir=std::make_pair(max_neighbor_len_L,0);
            }
        }
        return len_dir;
    };

    vec< std::pair<int64_t,int> > candidates;
    candidates.reserve(m_new.EdgeObjectCount());
    std::unordered_set<int> frozen_vertices;
    vec<int> e_delete, v_delete;

    for(bool bLoop=true; bLoop;){
        bLoop=false;

        candidates.clear();
        frozen_vertices.clear();
        e_delete.clear();
        v_delete.clear();

        //there can be dead edges
        for(int vv=0;vv<m_new.N();++vv){
            for(int jj=0;jj<m_new.FromSize(vv);++jj){
                int ee = m_new.EdgeObjectIndexByIndexFrom(vv,jj);
                if( ee != m_new.From(vv)[jj] && m_new.EdgeObject(ee).GetLen() < limit){
                    auto len = std::get<0>(GetLongestNeighbor(ee));
                    if(len >= 0) candidates.emplace_back(len,ee);
                }
            }
        }
        if(candidates.size()==0){
            bLoop=false;
            break;
        };
        std::sort(candidates.begin(),candidates.end());

        for(const auto& candidate:candidates){
            int ee = candidate.second;
            const int vv = to_left[ee];
            const int ww = to_right[ee];
            ForceAssert(vv!=ww);
            // if both vertices are untouched
            if( frozen_vertices.find(vv)==frozen_vertices.end() && frozen_vertices.find(ww)==frozen_vertices.end() )
            {
                frozen_vertices.insert(vv);
                frozen_vertices.insert(ww);
                auto info = GetLongestNeighbor(ee);
                if(std::get<0>(info)>0){
                    // get all short edges going v->w
                    vec<int> short_eids;
                    for(int jj=0;jj<m_new.FromSize(vv);++jj){
                        const auto& loc_eid = m_new.EdgeObjectIndexByIndexFrom(vv,jj);
                        if( m_new.From(vv)[jj] == ww && m_new.EdgeObject(loc_eid).GetLen() < limit ){
                            short_eids.push_back(loc_eid);
                        }
                    }
                    if(short_eids.size()>0){
                        if(std::get<1>(info)==1){ // to the right
                            // if all vertices z (v->w->z) has not been touched
                            if(!std::any_of(m_new.From(ww).begin(),m_new.From(ww).end()
                                           ,[&](int nv){return frozen_vertices.find(nv)!=frozen_vertices.end();}
                                           )){
                                //delete vertex w, if only short edges are going into it
                                bool bDeleteVertex = m_new.ToSize(ww) == short_eids.isize();

                                // for each short edge going v->w, merge with all edges w->z, resulting in edges v->z
                                for(const auto& eid: short_eids){
                                    for(int jj=0;jj<m_new.FromSize(ww);++jj){
                                        int nv = m_new.From(ww)[jj];
                                        frozen_vertices.insert(nv);
                                        m_new.AddEdge(vv,nv ,edge_t(m_new.EdgeObject(eid) ,m_new.EdgeObjectByIndexFrom(ww,jj) ,m_org.K()));
                                        to_left.push_back(vv);
                                        to_right.push_back(nv);
                                    }
                                }
                                if(bDeleteVertex){ v_delete.push_back(ww); }
                                else{ e_delete.append(short_eids); }
                                bLoop=true;
                            }
                        }
                        else{//to the left
                            if(!std::any_of(m_new.To(vv).begin(),m_new.To(vv).end()
                                           ,[&](int nv){return frozen_vertices.find(nv)!=frozen_vertices.end();}
                                           )){
                                bool bDeleteVertex = m_new.FromSize(vv) == short_eids.isize();
                                for(const auto& eid: short_eids){
                                    for(int jj=0;jj<m_new.ToSize(vv);++jj){
                                        int nv = m_new.To(vv)[jj];
                                        frozen_vertices.insert(nv);
                                        m_new.AddEdge(nv,ww ,edge_t(m_new.EdgeObjectByIndexTo(vv,jj)
                                                                       ,m_new.EdgeObject(eid)
                                                                       ,m_org.K()));
                                        to_left.push_back(nv);
                                        to_right.push_back(ww);
                                    }
                                }
                                if(bDeleteVertex){ v_delete.push_back(vv); }
                                else{ e_delete.append(short_eids); }
                                bLoop=true;
                            }
                        }
                    }
                }
            }
        }

        for(auto vv: v_delete){
            for(int jj=0;jj<m_new.FromSize(vv);++jj){
                e_delete.push_back(m_new.EdgeObjectIndexByIndexFrom(vv,jj));
            }
            for(int jj=0;jj<m_new.ToSize(vv);++jj){
                e_delete.push_back(m_new.EdgeObjectIndexByIndexTo(vv,jj));
            }
        }
        m_new.DeleteEdges(e_delete);

        //uncommenting the following will yield serious slow down in graph-to-ref alignment, but gain specificity
        /*
        for(int vv=0;vv<m_new.N();++vv){
            if( m_new.ToSize(vv)==1 && m_new.FromSize(vv)==1
                && m_new.From(vv)[0]!=m_new.To(vv)[0]
                && m_new.From(vv)[0]!=vv
                && m_new.To(vv)[0]!=vv
                ) {
                int nv = m_new.To(vv).front();
                int nw = m_new.From(vv).front();
                m_new.JoinEdges(vv,edge_t( m_new.EdgeObjectByIndexTo(vv,0)
                                         , m_new.EdgeObjectByIndexFrom(vv,0)
                                         , m_org.K()));
                to_left.push_back(nv);
                to_right.push_back(nw);
                bLoop=true;
            }
        }
        */

        if(bLoop && m_new.EdgeObjectCount() > std::numeric_limits<int>::max()/4){
            m_new.RemoveDeadEdgeObjects();
            m_new.RemoveEdgelessVertices();
            m_new.ToLeft(to_left);
            m_new.ToRight(to_right);
        }
    }
    for(int vv=0;vv<m_new.N();++vv){
        if( m_new.ToSize(vv)==1 && m_new.FromSize(vv)==1
            && m_new.From(vv)[0]!=m_new.To(vv)[0]
            && m_new.From(vv)[0]!=vv
            && m_new.To(vv)[0]!=vv
            ) {
            m_new.JoinEdges(vv,edge_t( m_new.EdgeObjectByIndexTo(vv,0)
                                     , m_new.EdgeObjectByIndexFrom(vv,0)
                                     , m_org.K()));
        }
    }
    m_new.RemoveDeadEdgeObjects();
    m_new.RemoveEdgelessVertices();
}

AssemblyError EvalAssemblyCore( const ref_data& ref, const HyperBasevector& hb,
        const vec<int>& inv, ostream& out, int iSWAk, int verbosity,
        N50Calculator* pN50Calc )
{
    ForceAssertEq( hb.EdgeObjectCount( ), inv.isize( ) );
    const vec<HyperBasevector>& GH = ref.GH; 
    const vec<bool>& is_circular = ref.is_circular;

    double rclock = WallClockTime( );
    if ( verbosity >= 2 )
        out << "assembly has " << hb.EdgeObjectCount( ) << " edges" << endl;
    // Test and expand reference.
    LinearRef linear_ref(GH,is_circular);
    if ( verbosity >= 1 ) out << Date( ) << ": building hbp" << endl;
    HyperBasevector hbp; 
    vec<pair<int,Bool>> hbp_to_hb;
    CreateHBPlus(hb, inv, hbp, hbp_to_hb);
    if ( verbosity >= 2 ) {    
        out << "'doubled' assembly has " << hbp.EdgeObjectCount( ) << " edges"
            << endl;    
    }
    if( verbosity >=3){
        ofstream ofs("hbp.dot");
        MakeHBPDot(ofs,hbp,hbp_to_hb);
    }
    // Heuristics.
    const int L = 40;
    const int min_dist = -10000;
    const int max_dist = 10000;
    int best_pass;
    int best_score;
    vec<int> best_path;
    string best_output;
    EdgePlacements edge_placements(hbp, hbp_to_hb, linear_ref.Seqs());
    if ( verbosity >= 1 ) out << Date() <<":b4 AlignEdgesToRefExp<L>" << endl;
    edge_placements.AlignEdgesToRefExp<L> (verbosity, out);
    if ( verbosity >= 1 ) out << Date() <<":after AlignEdgesToRefExp<L>" << endl;
    if ( verbosity >= 1 ) out << Date() <<":b4 RemoveBadPlacements()" << endl;
    edge_placements.RemoveBadPlacements();
    if ( verbosity >= 1 ) out << Date() <<":after RemoveBadPlacements()" << endl;
    if ( verbosity >= 1 ) out << Date() <<":b4 TwiddleSmart()" << endl;
    edge_placements.TwiddleSmart();
    if ( verbosity >= 1 ) out << Date() <<":after TwiddleSmart()" << endl;
    vec< vec<int> > spaths;
    vec< triple<int,int,int> > spaths_egd;
    vec< pair<int,int> > spaths_gg_pen;
    GraphZ gz(edge_placements, &PenaltyFunc);
    gz.FindShortestPath(min_dist, max_dist, spaths, spaths_egd, spaths_gg_pen,
            out, verbosity);
    vec<vec<int>> best_path_local, eids;
    int best_g;
    gz.FindBestEdgePath(spaths_egd, spaths, best_path_local, eids, best_g);

    vec<vec<std::pair<int,int>>> best_path_local_limits(best_path_local.size());
    for(size_t gg=0;gg<best_path_local.size();++gg){
        for(auto eid: eids[gg]){
            best_path_local_limits[gg].emplace_back( gz.verts[gz.edges[eid].first].second, gz.verts[gz.edges[eid].second].second);
        }
    }
    /*
    for(size_t gg=0;gg<best_path_local.size();++gg){
        for(size_t pp=0;pp<best_path_local[gg].size();++pp){
            std::cout << best_path_local[gg][pp] << " "
                      << eids[gg][pp] << " "
                      << hbp.EdgeObject(gz.edges[eids[gg][pp]].third.first).size() << " "
                      << best_path_local_limits[gg][pp].first << " "
                      << best_path_local_limits[gg][pp].second << " "
                      << std::endl;
        }
    }
    */
    AssemblyError res;
//    {
//        ofstream ofs("bayo_ref.fasta");
//        for ( int g = 0; g < (int) edge_placements.G.size( ); g++ ) {
//            ofs << ">" << g<<"\n";
//            ofs << edge_placements.G[g].ToString() << std::endl;
//        }
//    }
    for ( int g = 0; g < (int) edge_placements.G.size( ); g++ ) {    
        out << "Full best path(g=" << g << "): ";
        PrintPath(best_path_local[g], hbp_to_hb, out);
        vec<vec<int>> path_segs, eid_segs;
        vec<vec<std::pair<int,int>>> limits_segs;
        BreakPath(hbp, best_path_local[g], eids[g], best_path_local_limits[g]
                 , path_segs, eid_segs, limits_segs);

        int iFront = std::numeric_limits<int>::max();
        int iBack = std::numeric_limits<int>::min();

        // how the sequence of best-path edges combined into a sequence
        // 0: begin on sequence
        // 1: end on sequence
        // 2: edge index
        // 3: begin of edge
        // 4: end of edge
        // 5: rough begin on reference
        // 6: rough end on reference
        vec<std::tuple<int64_t,int64_t,int,int64_t,int64_t,int64_t,int64_t>> coors_hbp;
        if(!linear_ref.IsDoubled(g)) // if reference is not doubled, due to circular reference
        {
            vec<align> aligns(path_segs.size());
            for (size_t i = 0; i < path_segs.size(); i++) {
                out << "segment " << i << ": ";
                PrintPath(path_segs[i], hbp_to_hb, out);
                decltype(coors_hbp) coors_hbp_loc;
                basevector seq_seg = edge_placements.BestSeq(path_segs[i], eid_segs[i],limits_segs[i],coors_hbp_loc);
                if ( verbosity >= 1 ) out << Date() <<":before long-to-long alignment" << endl;
                aligns[i] = BestAlign(seq_seg, edge_placements.G[g], out, iSWAk);
                if ( verbosity >= 1 ) out << Date() <<":after long-to-long alignment" << endl;
                AddErrors(aligns[i], seq_seg, edge_placements.G[g], pN50Calc, res.GetErrPerRef(g));
//                std::cout << aligns[i].pos1() << " "
//                          << aligns[i].Pos1() << " "
//                          << aligns[i].pos2() << " "
//                          << aligns[i].Pos2() << std::endl;
                for(auto& entry: coors_hbp_loc){
//                    std::cout << std::get<0>(entry) << " " <<  std::get<1>(entry)  << " --> " ;
                    std::get<0>(entry) = aligns[i].PosOn2(std::get<0>(entry));
                    std::get<1>(entry) = aligns[i].PosOn2(std::get<1>(entry)-1)+1;
//                    std::cout << std::get<0>(entry) << " " <<  std::get<1>(entry) << " ; "
//                              << std::get<2>(entry) << " " << hbp_to_hb[std::get<2>(entry)].first << " ; "
//                              << std::get<3>(entry) << " " << std::get<4>(entry) << " ; "
//                              << std::endl;
                }
                coors_hbp.append(coors_hbp_loc);

                iFront = min(iFront,aligns[i].pos2());
                iBack  = max(iBack,aligns[i].Pos2());
            }
            for (size_t i = 1; i < aligns.size(); i++)
                res.AddGap(g,aligns[i-1].Pos2(),aligns[i].pos2() - aligns[i-1].Pos2());
        }
        else{ // for doubled reference, which has a high chance of breaking the underlying "hash-based" aligner
              // back to n^2 behavior.
              // this is to code around the deficiency, and should not be needed with a more resiliant aligner
            basevector full_seq;
            for (size_t i = 0; i < path_segs.size(); i++) {
                decltype(coors_hbp) coors_hbp_loc;
                int64_t end = full_seq.size();
                full_seq.append( edge_placements.BestSeq(path_segs[i], eid_segs[i],limits_segs[i],coors_hbp_loc) );
                for(auto& entry:coors_hbp_loc){
                    std::get<0>(entry)+=end;
                    std::get<1>(entry)+=end;
                }
                coors_hbp.append(coors_hbp_loc);
            }
            align alignment = BestAlignForDoubled(full_seq, edge_placements.G[g], coors_hbp,out);
//std::cout << "new alignmet" << std::endl;
//PrintVisualAlignmentClean(True, cout, full_seq, edge_placements.G[g], alignment);
//std::cout << "new alignmet end" << std::endl;
            AddErrors(alignment, full_seq, edge_placements.G[g], pN50Calc, res.GetErrPerRef(g));
//            std::cout << alignment.pos1() << " " << alignment.Pos1() << " "
//                      << alignment.pos2() << " " << alignment.Pos2() << std::endl;
            for(auto& entry: coors_hbp){
//                std::cout << std::get<0>(entry) << " " <<  std::get<1>(entry)  << " --> " ;
                std::get<0>(entry) = alignment.PosOn2(std::get<0>(entry));
                std::get<1>(entry) = alignment.PosOn2(std::get<1>(entry)-1)+1;
//                std::cout << std::get<0>(entry) << " " <<  std::get<1>(entry) << " ; "
//                          << std::get<2>(entry) << " " << hbp_to_hb[std::get<2>(entry)].first << " ; "
//                          << std::get<3>(entry) << " " << std::get<4>(entry) << " ; "
//                          << std::endl;
            }
            iFront = alignment.pos2();
            iBack  = alignment.Pos2();
        }
        res.AddLeftGap(g,iFront);
        res.AddRightGap(g,edge_placements.G[g].isize()-iBack);

        res.GetErrPerRef(g).SetRange(linear_ref.Seq(g).size(),linear_ref.IsDoubled(g));
        res.GetErrPerRef(g).SetBestPath(best_path_local[g],hbp_to_hb);
        res.GetErrPerRef(g).Compact(hbp,coors_hbp,hbp_to_hb);
    }
    return res;
}


}

//set iSWAk to 501 to use hash-based aligner to avoid n^2 blow up
AssemblyError EvalAssembly( const ref_data& ref, const HyperBasevector& hb,
        const vec<int>& inv, ostream& out, int iSWAk, int verbosity,
        N50Calculator* pN50Calc )
{
    if ( verbosity >= 1 ) out << Date() << ":EvalAssembly" << endl;
    AssemblyError  res;
    bool has_circular = std::count(ref.is_circular.begin(),ref.is_circular.end(),true)>0;
    if(has_circular) {
        /* This is for circular reference. Basically, the graph transformation required
         * is different from those of the linear ones (test case specific)
         * at this time there is very little circular test cases so only unroll Assembly loops
         * is needed -- but in the future, one should peform the transformations as
         * done for  linear-reference
         */
        vec<std::pair<int,int>> new_to_old;
        HyperBasevector hb_loc;
        UnrollAssemblyLoops(hb_loc,hb,new_to_old);
        vec<int> inv_loc;
        FixInversion(hb_loc,inv_loc);
        res = EvalAssemblyCore( ref, hb_loc, inv_loc, out, iSWAk, verbosity, pN50Calc );

        for(auto& entry: res){
            auto& err_per_ref = entry.second;
            const auto& bp_new = err_per_ref.GetBestPath();
            vec<pair<int,Bool>> bp_old; bp_old.reserve(bp_new.size());
            for(const auto& entry: bp_new){
                bp_old.insert(bp_old.end()
                             ,new_to_old[entry.first].second
                             ,make_pair(new_to_old[entry.first].first,entry.second));
            }
            err_per_ref.SetBestPath(bp_old);
        }

    }
    else{
        /* The older assembly evaluation has significant specificity problem when the graph
         * features short edges, as well as loopbacks involving short edges.
         * The following is a series of graph transformation to "unroll" those features
         * into long edges. The goal is to have longer edges to reduce multiple-alignment
         * of an edge, which tends to lead to failure of the evaluation.
         * The price to pay includes
         * a) the (in principle) expoential increase in the number of edges, which is
         *    rather moderate due to the seletion algorithm
         * b) longer edge in costing slower alignment, there is no free lunch in gaining
         *    specificity, but breaking of long edges tend to minimize such cost
         *
         * See the member function's declaration for description. fosmid-52 seems to be a good
         * illustration of the algorithm (look at the dot)
         */
        if ( verbosity >= 1 ) out << Date() << ":calling GraphTransformer" << endl;
        const int64_t absorb_limit = hb.K()*2+1-2;
        const int64_t unroll_limit = absorb_limit + hb.K();
        const int64_t break_tgt = unroll_limit + 1000; //1000bp power
        GraphTransformer engine(hb);
        if(verbosity>=1) {
            ofstream ofs("gt_init.dot");
            engine.MakeDot(ofs);
        }
        if ( verbosity >= 1 ) out << Date() << ":b4 GraphTransformer::ble" << endl;
        engine.BreakLongEdges(absorb_limit, break_tgt);
        if ( verbosity >= 1 ) out << Date() << ":after GraphTransformer::ble" << endl;
        if(verbosity>=1) {
            ofstream ofs("gt_ble.dot");
            engine.MakeDot(ofs);
        }
        if ( verbosity >= 1 ) out << Date() << ":b4 GraphTransformer::2vlb" << endl;
        engine.BreakTwoVerticesLoopBack(absorb_limit);
        if ( verbosity >= 1 ) out << Date() << ":after GraphTransformer::2vlb" << endl;
        if(verbosity>=1) {
            ofstream ofs("gt_2loop.dot");
            engine.MakeDot(ofs);
        }
        if ( verbosity >= 1 ) out << Date() << ":b4 GraphTransformer::ase" << endl;
        engine.AbsorbShortEdges(absorb_limit);
        if ( verbosity >= 1 ) out << Date() << ":after GraphTransformer::ase" << endl;
        if(verbosity>=1) {
            ofstream ofs("gt_abs.dot");
            engine.MakeDot(ofs);
        }
        if ( verbosity >= 1 ) out << Date() << ":b4 GraphTransformer::use" << endl;
        engine.UnrollShortLoops(unroll_limit);
        if ( verbosity >= 1 ) out << Date() << ":after GraphTransformer::use" << endl;
        if(verbosity>=1) {
            ofstream ofs("gt_unroll.dot");
            engine.MakeDot(ofs);
        }
        if ( verbosity >= 1 ) out << Date() << ":getting new edges" << endl;
        vec<basevector> edges;
        engine.getNewEdges(edges,break_tgt);
        HyperBasevector hb_loc(hb.K(),engine.getNewGraph(),edges);
        vec<int> inv_loc;
        FixInversion(hb_loc,inv_loc);
        if ( verbosity >= 1 ) out << Date() << ":calling core" << endl;
        res = EvalAssemblyCore( ref, hb_loc, inv_loc, out, iSWAk, verbosity, pN50Calc );
        if ( verbosity >= 1 ) out << Date() << ":done core" << endl;
        for(auto& entry: res){
            auto& err_per_ref = entry.second;
            /*
            const auto& bp_new = err_per_ref.GetBestPath();
            vec<pair<int,Bool>> bp_old; bp_old.reserve(bp_new.size());
            for(const auto& entry: bp_new){
                const auto& seq = engine.new_to_old(entry.first);
                if( entry.second ){ //if forward
                    std::transform(seq.begin(),seq.end(),std::back_inserter(bp_old),[&](int i){return make_pair(i,entry.second);});
                }
                else{
                    std::transform(seq.rbegin(),seq.rend(),std::back_inserter(bp_old),[&](int i){return make_pair(i,entry.second);});
                }
            }
            err_per_ref.SetBestPath(bp_old);
            */
            //report the best path in the original graph's coordinate
            err_per_ref.SetBestPath( engine.new_to_old(err_per_ref.GetBestPath() )    ) ;
        }
//        res = EvalAssemblyCore( ref, hb, inv, out, iSWAk, verbosity );
    }
    res.PrintSummary(out);
    if ( verbosity >= 1 ) out << Date() << ":EvalAssembly done " << endl;
    return res;
}

void AssemblyError::err_per_ref::Compact(const HyperBasevector& hbp
                                        ,vec<std::tuple<int64_t,int64_t,int,int64_t,int64_t,int64_t,int64_t>> coors_hbp
                                        ,const vec<pair<int,Bool>>&hbp_to_hb){
    const int verbosity = 0;
    ForceAssert(range<=max_range);
    if(range == max_range){
        m_GoodTopology=true;
        return;
    }
    if(verbosity) std::cout <<"AssemblyError::err_per_ref::Compact(): compacting range." << std::endl;

    std::sort(coors_hbp.begin(),coors_hbp.end());
    bool bCompactable=    coors_hbp.size() > 0
                       && std::get<0>(coors_hbp.front()) >=0
                       && std::get<1>(coors_hbp.front()) >= std::get<0>(coors_hbp.front()) ;
//    std::cout << std::get<0>(coors_hbp.front()) << " " <<  std::get<1>(coors_hbp.front())  << std::endl;
    for(size_t ii=1;ii<coors_hbp.size();++ii){
        bCompactable =    bCompactable
                       && std::get<0>(coors_hbp[ii]) >= std::get<1>(coors_hbp[ii-1])
                       && std::get<1>(coors_hbp[ii]) >= std::get<0>(coors_hbp[ii]);
//        std::cout << std::get<0>(coors_hbp[ii]) << " " <<  std::get<1>(coors_hbp[ii])  << std::endl;
    }

    if(! bCompactable ){
        m_GoodTopology=false;
        viLGap.push_back( -1 );
        viRGap.push_back( -1 );
        if(verbosity) std::cout <<"AssemblyError::err_per_ref::Compact(): circular reference not compacted due to bad mapped coordinates." << std::endl;
        return;
    }
    vec<std::pair<int,int64_t>> base_coor(max_range,std::make_pair(-1,-1LL));
    for(int cc=0;cc<coors_hbp.isize();++cc){
        std::fill(base_coor.begin()+ std::get<0>(coors_hbp[cc])
                 ,base_coor.begin()+ std::get<1>(coors_hbp[cc])
                 ,std::make_pair(cc,-1));
    }

    vec<bool> is_sub(max_range,false);
    vec<bool> is_del(max_range,false);
    vec<int64_t> ins_count(max_range,0);

    for(size_t bb=0;bb<base_coor.size();++bb){
        if(base_coor[bb].first<0) is_del[bb]=true;
//        else
//            std::cout << bb << " check " << base_coor[bb].first << " " << base_coor[bb].second << std::endl;
    }
//    std::cout << std::count(is_del.begin(),is_del.end(),true) << " bases not covered" << std::endl;

    ForceAssert(viLGap.size()==1);
    std::fill(is_del.begin(),is_del.begin()+viLGap.front(),true);
//    std::cout << std::count(is_del.begin(),is_del.end(),true) << " bases not covered" << std::endl;

    ForceAssert(viRGap.size()==1);
    std::fill(is_del.end()-viRGap.front(),is_del.end(),true);
//    std::cout << std::count(is_del.begin(),is_del.end(),true) << " bases not covered" << std::endl;

    for(const auto&entry:indels){
        ForceAssert(entry.first>=0 && entry.first < max_range && entry.second!=0);
        if( entry.second > 0 )   { ins_count[entry.first]+=entry.second; }
        else if(entry.second < 0){ std::fill(is_del.begin()+entry.first,is_del.begin()+entry.first-entry.second,true); }
    }
//    std::cout << std::count(is_del.begin(),is_del.end(),true) << " bases not covered" << std::endl;

    for(int cc=0;cc<coors_hbp.isize();++cc){
        auto seq1 = std::get<0>(coors_hbp[cc]);
        auto seq2 = std::get<1>(coors_hbp[cc]);
        auto edge1 = std::get<3>(coors_hbp[cc]);
        auto edge2 = std::get<4>(coors_hbp[cc]);
        auto ee=edge1;
        for(auto ss=seq1;ss<seq2;++ss){
            ForceAssert(base_coor[ss].first==cc);
            if( ins_count[ss] > 0 && ss!=seq1){
                ee+=ins_count[ss];
            }
            if(!is_del[ss]){
                base_coor[ss].second=ee;
                ++ee;
            }
        }
    }
//    for(size_t bb=0;bb<base_coor.size();++bb){
//        std::cout << "> " << bb << " " << base_coor[bb].first << " " << base_coor[bb].second << std::endl;
//    }
    auto coor_is_circular=[&base_coor,&coors_hbp](int64_t left,int64_t right){
        const auto& ll = base_coor[left];
        std::pair<int,int64_t> rr = (right<0)?std::make_pair(-1,int64_t(-1)):base_coor[right];
        if( ll.first >=0 && ll.first >= 0){
//        std::cout << left << " " << right << "->"
//                  << ll.first << " " << ll.second  << ":"
//                  << std::get<0>(coors_hbp[ll.first]) << " "
//                  << std::get<1>(coors_hbp[ll.first]) << " "
//                  << std::get<2>(coors_hbp[ll.first]) << " "
//                  << std::get<3>(coors_hbp[ll.first]) << " "
//                  << std::get<4>(coors_hbp[ll.first]) << " "
//                  << "|"
//                  << rr.first << " " << rr.second << ":"
//                  << std::get<0>(coors_hbp[rr.first]) << " "
//                  << std::get<1>(coors_hbp[rr.first]) << " "
//                  << std::get<2>(coors_hbp[rr.first]) << " "
//                  << std::get<3>(coors_hbp[rr.first]) << " "
//                  << std::get<4>(coors_hbp[rr.first]) << " "
//                  << std::endl;
        }
        return    ll.first<0
               || rr.first<0
               || (std::get<2>(coors_hbp[ll.first]) == std::get<2>(coors_hbp[rr.first]) && ll.second == rr.second )
               ;
    };

    for(const auto&entry:subs){
        ForceAssert(entry.first>=0 && entry.first < max_range && entry.second==1);
        is_sub[entry.first]=true;
    }

    int64_t nError=0;
    for(int64_t bb=0;bb<range;++bb){ nError += is_sub[bb] + is_del[bb] + ins_count[bb]; }

    circular_tester tester(hbp);
    tester.UpdateSeq(false,true,(base_coor[0].first<0)?-1: std::get<2>(coors_hbp[base_coor[0].first]));
    for(int64_t bb=1;bb<range;++bb){
        if( base_coor[bb].first != base_coor[bb-1].first){
            tester.UpdateSeq(false,true,(base_coor[bb].first<0)?-1: std::get<2>(coors_hbp[base_coor[bb].first]));
        }
    }

    int64_t min_error = (tester.IsCircular()&&coor_is_circular(0,range))?nError:std::numeric_limits<int64_t>::max();
    uint64_t best_front = 0;
    vec<int> best_seq(tester.IsCircular()?tester.getSeq():vec<int>());

    int64_t min_error_linear(min_error);
    uint64_t best_front_linear(best_front);
    vec<int> best_seq_linear(best_seq);


    for( int64_t front=1,back=range ; back < max_range ; ++front, ++back){
        nError -= is_sub[front-1] + is_del[front-1] + ins_count[front-1];
        nError += is_sub[back] + is_del[back] + ins_count[back];
        tester.UpdateSeq( base_coor[front].first!=base_coor[front-1].first
                        , base_coor[back].first != base_coor[back-1].first
                        , (base_coor[back].first < 0)?-1:std::get<2>(coors_hbp[base_coor[back].first]));
        if( (tester.IsCircular()&&coor_is_circular(front,(back+1==max_range)?-1:back+1)) && nError < min_error){
            min_error = nError;
            best_front = front;
            best_seq=tester.getSeq();
        }
        if(nError < min_error_linear){
            min_error_linear=nError;
            best_front_linear = front;
            best_seq_linear=tester.getSeq();
        }
    }
    viLGap.clear();
    viRGap.clear();
    gaps.clear();
    indels.clear();
    subs.clear();

    if(verbosity) std::cout <<"AssemblyError::err_per_ref::Compact(): (min_error front)=(" << min_error << " " << best_front << ")" << std::endl;
    if(min_error==std::numeric_limits<int64_t>::max()){
        if(verbosity) std::cout <<"AssemblyError::err_per_ref::Compact(): circular reference not compacted due to bad circular mapping." << std::endl;

        viLGap.push_back( min(-1,-int(base_coor[best_front_linear].second)) );

        auto end_element = base_coor[(best_front_linear+range)%max_range];
        viRGap.push_back( min(-1,-int( std::get<2>(coors_hbp[end_element.first]) - end_element.second )) );

        best_front=best_front_linear;
        best_seq=best_seq_linear;
        m_GoodTopology=false;
    }
    else{
        m_GoodTopology=true;
    }


    int iDelCount=0;
    int64_t block_start=-1;

    for(int64_t bb=0;bb<range;++bb){
        int64_t idx = best_front+bb;
        if( iDelCount > 0 && (is_sub[idx] || ins_count[idx] || !is_del[idx]) ){
            AddIndel(block_start,-iDelCount);
            block_start=-1;
            iDelCount=0;
        }
        if(is_sub[idx]){
            AddSub(bb,1);
        }
        if(ins_count[idx]>0){
            AddIndel(bb,ins_count[idx]);
        }
        if(is_del[idx]){
            if(iDelCount==0){
                block_start=bb;
            }
            ++iDelCount;
        }
    }
    if(iDelCount>0){
        AddIndel(block_start,-iDelCount);
    }
    best_seq.erase(std::remove_if(std::begin(best_seq),std::end(best_seq),[](int64_t in){return in < 0;}),std::end(best_seq));
    SetBestPath(best_seq,hbp_to_hb);
    if(verbosity) std::cout <<"AssemblyError::err_per_ref::Compact(): reducint max range from " << max_range << " to "  << range << std::endl;
    max_range=range;
}
