///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "paths/long/VariantReadSupport.h"

#include "CoreTools.h"
#include "Map.h"
#include "paths/long/Logging.h"
#include "paths/long/SupportedHyperBasevector.h"
#include "pairwise_aligners/SmithWatAffine.h"
#include "PrintAlignment.h"
#include "util/TextTable.h"
#include "graph/Digraph.h"
#include "graph/DigraphTemplate.h"
#include "paths/long/ReadOriginTracker.h"
#include "paths/long/EvalByReads.h"

namespace {
void CalcLengthProbSimple(vec<double>&vOut,const double dProbIns_len, const double dProbDel_len, const uint64_t org_length){
    const double dProbIns=dProbIns_len;
    const double dProbDel=dProbDel_len;
    const double dProbNor=1.-dProbIns-dProbDel;
    ForceAssert(   dProbNor >=0.0 && dProbIns >=0.0 && dProbDel >=0.0
                && dProbNor <=1.0 && dProbIns <=1.0 && dProbDel <=1.0 );

    vOut.clear();
    vOut.resize( 2*org_length+1 , 0.0);
    vOut[org_length] = dProbNor;
    if( org_length >0) vOut[org_length-1] = dProbDel;
    vOut[org_length+1] = dProbIns;
}
void CalcLengthProb(vec<double>&vOut,const double dProbIns_len, const double dProbDel_len, const uint64_t org_length){
    const double dProbIns=dProbIns_len/double(org_length);
    const double dProbDel=dProbDel_len/double(org_length);
    const double dProbNor=1.-dProbIns-dProbDel;
    ForceAssert(   dProbNor >=0.0 && dProbIns >=0.0 && dProbDel >=0.0
                && dProbNor <=1.0 && dProbIns <=1.0 && dProbDel <=1.0 );

    vOut.clear();
    vOut.resize( 2*org_length+1 , 0.0);


    double dNextZeroCol = pow(dProbNor,double(org_length));
    const uint64_t nRow=org_length+1;
    const uint64_t nCol=org_length+1;
    for( uint64_t row = 0 ; row < nRow ; ++row){
        double dColValue = dNextZeroCol ;
        dNextZeroCol *= dProbIns / dProbNor * ( org_length - row ) / (row+1.0) ;
        for( uint64_t col = 0 ; col < nCol ; ++col){
            vOut[org_length + row - col] += dColValue;
            dColValue *= dProbDel / dProbNor * ( org_length - row - col) / (col+1.0) ;
        }
    }
};

void ComputeHomopolymerProb_AT_Diploid( vec<std::tuple<uint64_t,double,uint64_t>>& n_q_len){
    const int verbosity=0;
    const double minQ=2.0;
    const uint64_t nFlavors = n_q_len.size();
    const double dTrustOfGapFreeAlignmentPrior=0.5;
    const double critical_ratio=0.40;
    ForceAssert(nFlavors>0);
    if(verbosity>0) std::cout.precision(5);
    uint64_t max_n=0;
    for( auto& entry: n_q_len){
        if( get<1>(entry) < minQ ){ get<0>(entry)=0; }
        else{ max_n=max(get<0>(entry),max_n); }
    }
    if(max_n>3){
        for( auto& entry: n_q_len){
            if( get<0>(entry) == 1 ){ get<1>(entry)=min(get<1>(entry),0.); }
        }
    }
    Sort(n_q_len,std::greater<std::tuple<uint64_t,double,uint64_t>>());

    if(verbosity>0) std::cout<<"in coming, reads:"<< std::endl;
    if(verbosity>0) for(const auto& entry:n_q_len){ std:: cout << get<0>(entry) << " " << get<1>(entry) << " " << get<2>(entry) << std::endl; }

    vec<vec<double>> vvdLengthProbs(nFlavors);
    vec<double> p_from_q(nFlavors);
    uint64_t nReads=0;

    if (get<0>(n_q_len[0]) ==0) return;
    for(uint64_t ff=0;ff<nFlavors;++ff){
        uint64_t length = get<2>(n_q_len[ff]);
        uint64_t nread_loc = get<0>(n_q_len[ff]);
        nReads+=nread_loc;
        p_from_q[ff] = (nread_loc>0) ? (1.0-pow(10.0 , -0.1*get<1>(n_q_len[ff]))): 0.0;

        const double P_D = 0.01 * max( 1.0 , 0.4 * (length-10.0) + 1.0 ) ;
        const double P_I = 0.01 * max( 0.25 , 0.1 * (length-10.0) + 0.25 ) ;

        CalcLengthProb(vvdLengthProbs[ff],P_I,P_D,length);
    }

    if(verbosity>0) std::cout<<"after init, reads:"<< std::endl;
    if(verbosity>0) for(const auto& entry:n_q_len){ std:: cout << get<0>(entry) << " " << get<1>(entry) << " " << get<2>(entry) << std::endl; }
    if( nReads==0 ) return;
    vec<long double> P_o_gt(nFlavors*nFlavors,0.0);
    vec<long double> P_gt_o(nFlavors*nFlavors,0.0);
    vec<long double> eta_gt(nFlavors*nFlavors,0.0);
    uint64_t nGT = nFlavors*(nFlavors+1)/2;

    long double dTmp=0;
    for( uint64_t row=0;row<nFlavors;++row){
        for( uint64_t col=0;col<nFlavors;++col){
            long double& dLoc = eta_gt[ row*nFlavors+col ];
            dLoc=1.0;
            for(size_t ff=0;ff<nFlavors;++ff){
                if( ff==row || ff==col){ dLoc*=p_from_q[ff]; }
                else                   { dLoc*=1.0-p_from_q[ff]; }
            }
            if( col >= row ){ dTmp += dLoc; }
        }
    }
    if( dTmp < std::numeric_limits<long double>::epsilon() * 10){
        for( uint64_t row=0;row<nFlavors;++row){
            for( uint64_t col=0;col<nFlavors;++col){
                eta_gt[ row*nFlavors+col ] = 1.0/(long double)(nGT);
            }
        }
    }
    else{
        dTmp = (1.0-dTmp) / (long double)(nGT);
        for( uint64_t row=0;row<nFlavors;++row){
            for( uint64_t col=0;col<nFlavors;++col){
                long double& dLoc = eta_gt[ row*nFlavors+col ];
                dLoc+=dTmp;
                dLoc*=dTrustOfGapFreeAlignmentPrior;
                dLoc+= (1.0-dTrustOfGapFreeAlignmentPrior) / (long double)(nGT);

            }
        }
    }
    if(verbosity>0){
        std::cout << "prior: " << std::endl;
        for( uint64_t row=0;row<nFlavors;++row){
            for( uint64_t col=0;col<nFlavors;++col){
                std::cout << eta_gt[row*nFlavors+col] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
        std::cout << std::endl;
    }

    long double dFactor = 1.0;
    //prevent underflow
    for( size_t ll = get<0>(n_q_len[0])+1 ; ll <= nReads ; ++ll){ dFactor *= ll; }
    for(uint64_t ff=1;ff<nFlavors;++ff){
        for( size_t ll = 2 ; ll <= get<0>(n_q_len[ff]) ; ++ll){ dFactor /= (long double)(ll); }
    }
    if( dFactor==std::numeric_limits<long double>::infinity() || isnan(dFactor) ){ dFactor=1.0; }
    if(verbosity>0) std::cout << "dFactor " << dFactor << std::endl;

    long double dDenom=0.0;
    for( uint64_t row=0;row<nFlavors;++row){
        for( uint64_t col=row;col<nFlavors;++col){
            long double& dLoc = P_o_gt[ row*nFlavors+col ];
            dLoc = dFactor*eta_gt[row*nFlavors+col];
            if(verbosity>0) std::cout << row << " " << col << " " << dLoc << std::endl;
            for( size_t li=0;li<nFlavors;++li){
                const auto& length = get<2>(n_q_len[li]);
                const auto& n = get<0>(n_q_len[li]);
                if(n!=0){
                    if(verbosity){
                        std::cout << row << " " << col << " " << li << " (" << length << "," << n << ") "
                                  << ((length<vvdLengthProbs[row].size())?  vvdLengthProbs[row][length] : 0.0) << " "
                                  << ((length<vvdLengthProbs[col].size())?  vvdLengthProbs[col][length] : 0.0) << " "
                                  << pow(0.5*( ((length<vvdLengthProbs[row].size())?  vvdLengthProbs[row][length] : 0.0)
                                              +((length<vvdLengthProbs[col].size())?  vvdLengthProbs[col][length] : 0.0)
                                             )
                                        , double(n)
                                  )<<std::endl;;
                    }
                    dLoc *=pow(0.5*( ((length<vvdLengthProbs[row].size())?  vvdLengthProbs[row][length] : 0.0)
                                    +((length<vvdLengthProbs[col].size())?  vvdLengthProbs[col][length] : 0.0)
                                   )
                              , double(n)
                              );
                }
            }
            if(verbosity>0) std::cout << row << " " << col << " " << dLoc << std::endl;
            dDenom+=dLoc;
        }
    }
    if(verbosity>0) std::cout <<"denom: " << dDenom<< std::endl;
    if( dDenom < std::numeric_limits<long double>::epsilon() * 100 || isnan(dDenom)) return;
    ForceAssert(dDenom>0);
    if(verbosity>0){
        for( uint64_t row=0;row<nFlavors;++row){
            for( uint64_t col=0;col<nFlavors;++col){
                std::cout << P_o_gt[row*nFlavors+col] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
    double dTotal=0.0;
    for( uint64_t row=0;row<nFlavors;++row){
        P_gt_o[ row*nFlavors+row] = P_o_gt[ row*nFlavors+row] / dDenom;
        dTotal+=P_gt_o[ row*nFlavors+row];
        for( uint64_t col=row+1;col<nFlavors;++col){
            P_gt_o[ row*nFlavors+col] = P_o_gt[ row*nFlavors+col] / dDenom;
            dTotal+=P_gt_o[ row*nFlavors+col];
            P_gt_o[ col*nFlavors+row] = P_gt_o[ row*nFlavors+col];
        }
    }
    if(verbosity>0){
        for( uint64_t row=0;row<nFlavors;++row){
            for( uint64_t col=0;col<nFlavors;++col){
                std::cout << P_gt_o[row*nFlavors+col] << " ";
                ForceAssert( P_gt_o[ col*nFlavors+row] == P_gt_o[ row*nFlavors+col] );
            }
            std::cout << std::endl;
        }
        std::cout<< "dTotal = " << dTotal << std::endl;
    }
//    ForceAssert ( fabs(dTotal-1.0) < 10*std::numeric_limits<double>::epsilon());
    ForceAssert ( dTotal > -10*std::numeric_limits<float>::epsilon());
    ForceAssert ( dTotal < 1.0+10*std::numeric_limits<float>::epsilon());
    for( auto& entry:P_gt_o){ entry /= dTotal; }
    for( uint64_t ff=0;ff<nFlavors;++ff){
        const auto itr = P_gt_o.begin() + ff*nFlavors;
        double p=std::accumulate(itr,itr+nFlavors,0.0);

        const auto itr0 = eta_gt.begin() + ff*nFlavors;
        double p0=std::accumulate(itr0,itr0+nFlavors,0.0);

//        if( p<p0 ){
            double q = -10.0 * log10(1.0-p);
            if ( q < minQ) q=0;
            if ( q == std::numeric_limits<double>::infinity() ) q = std::numeric_limits<double>::max();
            get<1>(n_q_len[ff]) = q;
//        }

        if(verbosity>0){
            std::cout << get<0>(n_q_len[ff]) << " "
                      << get<1>(n_q_len[ff]) << " "
                      << get<2>(n_q_len[ff]) << " "
                      <<std::endl;
        }
    }
};

void AdjustHomopolymerProb(unsigned char baseval, std::unordered_map< uint64_t , vec<std::tuple<uint64_t,double,int>>>& length_lookup){
    vec<std::tuple<uint64_t,double,uint64_t>> n_q_len;
    uint64_t max_length=0;
    uint64_t nReads=0;
    for( const auto&len_entry: length_lookup){
        uint64_t count=0;
        double qsum=0.0;
        for(const auto& var_entry: len_entry.second){
            count+=std::get<2>(var_entry);
            qsum+=std::get<1>(var_entry);
        }
        n_q_len.push_back( make_tuple(count,qsum,len_entry.first));
        max_length = std::max(max_length,len_entry.first);
        nReads+=count;
    }
    if( baseval == Base::char2Val('A') || baseval==Base::char2Val('T')){ // if it's a A/T homopolymer
        const uint64_t critical_length = 10; // this seems to be where we want to draw the line as discussed in the meeting
        if ( max_length < critical_length || nReads ==0) return;
        ComputeHomopolymerProb_AT_Diploid( n_q_len);
    }
    else{
        return ;
    }
    for(auto entry: n_q_len){
        auto loc_length = get<2>(entry);
        auto loc_q = get<1>(entry);
        for( auto& var_entry: length_lookup[loc_length]){ get<1>(var_entry)= std::min( get<1>(var_entry), loc_q); }
    }
}
void RemoveDanglingCalls( map<Variant, vec<pair<double,double>>>& probs
                        , const vec<VariantCallGroup>& vcall_groups
                        ){
    const double dHighQ=20;
    const double dLowQ=4;
    const int64_t critical_factor=4;

    if(probs.size() ==0 ) return;
    size_t nSamples = (*probs.begin()).second.size();
    if(nSamples ==0 ) return;

    for( const VariantCallGroup& vgroup: vcall_groups){
        map<Variant, set<int>> variant_branches;
        const auto& vcalls=vgroup.GetVariantCalls();

        const size_t nBranches = vcalls.size();
        for (size_t branch = 0; branch < nBranches; branch++) {
            for (const VariantCall& x: vcalls[branch]){
                variant_branches[x.variant].insert(branch);
            }
        }

//actually it doesn't matter that they overlap since we're looking for dangling calls
        //variant must not overlap
//        int64_t back=std::numeric_limits<int64_t>::min();
//        bool bProceed=true;
//        for( const auto& entry: variant_branches){
//            if( entry.first.pos <= back){
//                bProceed=false;
//            }
//            back = max( back , int64_t(entry.first.pos+entry.first.ref.size()-1));
//        }
//        if(!bProceed ) continue;

        vec<vec<size_t>>                           appearance(nBranches,vec<size_t>(nSamples,0));
        vec<vec<size_t>>                        must_be_there(nBranches,vec<size_t>(nSamples,0));
        vec<vec<map<Variant,std::pair<bool,bool> > > > hinges(nBranches,vec<map<Variant ,std::pair<bool,bool> >> (nSamples));

        for (auto it = variant_branches.begin(); it != variant_branches.end(); ++it) { // for each variant
            const Variant& var = it->first;
            const set<int>& pathids = it->second;

            vec<int> ref_list,alt_list;
            for (size_t branch = 0; branch < nBranches; branch++) {
                if(pathids.find(branch)!=pathids.end()){
                    alt_list.push_back(branch);
                }
                else{
                    ref_list.push_back(branch);
                }
            }
            auto p_itr = probs.find(var);
            if(p_itr==probs.end())continue;
            for (size_t sample = 0; sample < nSamples; sample++){ //for each sample
                double p_ref = p_itr->second[sample].first;
                if(p_ref >= dLowQ){
                    for( auto entry: ref_list){
                        ++appearance[entry][sample];
                    }
                    if(ref_list.size()==1 && p_ref >= dHighQ){
                        ++must_be_there[ref_list.front()][sample];
                        hinges[ref_list.front()][sample][var].first=true;
                    }
                }
                double p_alt = p_itr->second[sample].second;
                if(p_itr->second[sample].second >= dLowQ){
                    for( auto entry: alt_list){
                        ++appearance[entry][sample];
                    }
                    if(alt_list.size()==1 && p_alt >= dHighQ){
                        ++must_be_there[alt_list.front()][sample];
                        hinges[alt_list.front()][sample][var].second=true;
                    }
                }
            }
        }
        for(size_t branch=0;branch<nBranches;++branch){
            for(size_t sample=0;sample<nSamples;++sample){
                if( must_be_there[branch][sample] && appearance[branch][sample] * critical_factor < variant_branches.size()){
                    for( auto itr=hinges[branch][sample].begin() ; itr!=hinges[branch][sample].end(); ++itr){
                        if( (*itr).second.first ){ probs[(*itr).first][sample].first=0; }
                        if( (*itr).second.second ){ probs[(*itr).first][sample].second=0; }
                    }
                }
            }
        }
    }
}

void FillInZeroProbAccordingToGroup( map<Variant, vec<pair<double,double>>>& probs
                                   , const vec<VariantCallGroup>& vcall_groups
                                   ){
    const double dHighQ=20;
    const double dLowQ=4;

    if(probs.size() ==0 ) return;
    size_t nSamples = (*probs.begin()).second.size();
    if(nSamples ==0 ) return;

    for( const VariantCallGroup& vgroup: vcall_groups){
        map<Variant, set<int>> variant_branches;
        const auto& vcalls=vgroup.GetVariantCalls();
        
        const size_t nBranches = vcalls.size();
        for (size_t branch = 0; branch < nBranches; branch++) {
            for (const VariantCall& x: vcalls[branch]){
                variant_branches[x.variant].insert(branch);
            }
        }
        vec<Variant> v_list;v_list.reserve(variant_branches.size());
        for(const auto&entry:variant_branches){ v_list.push_back(entry.first); }
        Sort(v_list ,[](const Variant&L,const Variant&R){ if(L.gid!=R.gid) return L.gid<R.gid;if(L.pos!=R.pos) return L.pos<R.pos; return L.ref.size()<R.ref.size();});
        ForceAssert(v_list.size()==variant_branches.size());
        
        if( v_list.size() == 0) continue;
        vec<vec<Variant>> var_clusters;
        int64_t back=std::numeric_limits<int64_t>::min();
        for( const Variant& var: v_list){
            if( back < var.pos) var_clusters.push_back(vec<Variant>());
            var_clusters.back().push_back(var);
            back = std::max(back, int64_t(var.pos+var.ref.size()-1));
        }
        
        vec<vec<double>> max_branch_q(nBranches,vec<double>(nSamples,0));
        vec<vec<size_t>> appearance(nBranches,vec<size_t>(nSamples,0));
        vec<vec<size_t>> must_be_there(nBranches,vec<size_t>(nSamples,0));
        vec<vec<int>> cluster_ref_lists;cluster_ref_lists.resize(var_clusters.size());
        
        for( const auto& cluster:var_clusters){
            set<int> cluster_alts;
            for(const auto& var: cluster){
                const set<int>& pathids = variant_branches[var];
                for(size_t branch=0;branch<nBranches;++branch){
                    if(pathids.find(branch)!=pathids.end()){
                        cluster_alts.insert(branch);
                    }
                }
            }
            cluster_ref_lists.push_back(vec<int>());
            cluster_ref_lists.back().reserve(nBranches-cluster_alts.size());
            for(size_t branch=0;branch<nBranches;++branch){ if(cluster_alts.find(branch)==cluster_alts.end()){ cluster_ref_lists.back().push_back(branch); } }
            const vec<int>& cluster_ref_list=cluster_ref_lists.back();
            
            for(const auto& var: cluster){
                auto p_itr = probs.find(var);
                if(p_itr==probs.end())continue;
                
                vec<int> alt_list;
                const set<int>& pathids = variant_branches[var];
                for(size_t branch=0;branch<nBranches;++branch){
                    if(pathids.find(branch)!=pathids.end()){
                        alt_list.push_back(branch);
                    }
                }
                for (size_t sample = 0; sample < nSamples; sample++){
                    double p_ref = p_itr->second[sample].first;
                    if(p_ref >= dLowQ){
                        for( auto entry: cluster_ref_list){
                            ++appearance[entry][sample];
                        }
                        if(cluster_ref_list.size()==1 && p_ref >= dHighQ){
                            max_branch_q[cluster_ref_list.front()][sample] = max(max_branch_q[cluster_ref_list.front()][sample],p_ref);
                            ++must_be_there[cluster_ref_list.front()][sample];
                        }
                    }
                    double p_alt = p_itr->second[sample].second;
                    if(p_itr->second[sample].second >= dLowQ){
                        for( auto entry: alt_list){
                            ++appearance[entry][sample];
                        }
                        if(alt_list.size()==1 && p_alt >= dHighQ){
                            max_branch_q[alt_list.front()][sample] = max(max_branch_q[alt_list.front()][sample],p_alt);
                            ++must_be_there[alt_list.front()][sample];
                        }
                    }
                }
            }
        }

        for( size_t cc=0 ; cc < var_clusters.size() ; ++cc){
            const auto& cluster=var_clusters[cc];
            const vec<int>& cluster_ref_list=cluster_ref_lists[cc];
            for(const auto& var: cluster){
                const set<int>& pathids = variant_branches[var];
    
                vector<bool> bvAlt(nBranches,false);
                for (size_t k = 0; k < vcalls.size(); k++) {
                    bvAlt[k] = pathids.find(k)!=pathids.end();
                }
                auto p_itr = probs.find(var);
                if(p_itr==probs.end())continue;
                
                for (size_t sample = 0; sample < nSamples; sample++){
                    for(size_t branch=0;branch<nBranches;++branch){
                        ForceAssert( appearance[branch][sample] <= variant_branches.size());
                        if(    must_be_there[branch][sample]
                            && (    //appearance[branch][sample] * 4 >= variant_branches.size() * 3
    //                              ||
                                  appearance[branch][sample] + 1 >= variant_branches.size()
                               )
                          ){
                            if(bvAlt[branch]){
                                if((*p_itr).second[sample].second<dLowQ){
                                    (*p_itr).second[sample].second = max((*p_itr).second[sample].second , max_branch_q[branch][sample]);
                                }
                            }
                            else{
                                if((*p_itr).second[sample].first<dLowQ && (cluster.size()==1||Member(cluster_ref_list,int(branch)))){
                                    (*p_itr).second[sample].first = max((*p_itr).second[sample].first , max_branch_q[branch][sample]);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void DumpGraph(const String filename, const HyperBasevector& hb) 
{
    vec<String> edge_names2(hb.EdgeObjectCount());
    vec<double> lengths2( hb.EdgeObjectCount( ) );
    for (size_t i = 0; i < edge_names2.size(); ++i) {
        edge_names2[i] = ToString(i);
        lengths2[i] = hb.EdgeLengthKmers(i);
    }
    ofstream dout(filename);
    hb.PrettyDOT( dout, lengths2, HyperBasevector::edge_label_info(
                HyperBasevector::edge_label_info::DIRECT, &edge_names2 ) );
}

// Extend the reference edge of the K=1 HyperBasevector.
void ExtendedBasevector(const HyperBasevector& hb, int edge, int ext_max,
        basevector& edge_extended, int& len_extended_left) 
{
    ForceAssertEq(hb.K(), 1);
    vec<int> to_left; hb.ToLeft(to_left);
    vec<int> to_right; hb.ToRight(to_right);

    vec<int> edges = {edge};

    int ext_left = 0;
    int node_left = to_left[edge];
    while (hb.To(node_left).size() > 0 && ext_left < ext_max) {
        int e_new = hb.EdgeObjectIndexByIndexTo(node_left, hb.To(node_left).size() -1);
        edges.push_front(e_new);
        ext_left += hb.EdgeLengthKmers(e_new);
        node_left = hb.To(node_left).back();
    }

    int ext_right = 0;
    int node_right = to_right[edge];
    while (hb.From(node_right).size() > 0 && ext_right < ext_max) {
        int e_new = hb.EdgeObjectIndexByIndexFrom(node_right, hb.From(node_right).size()-1);
        edges.push_back(e_new);
        node_right = hb.From(node_right).back();
        ext_right += hb.EdgeLength(e_new);
    }
    basevector edges_to_bases = hb.EdgePathToBases(edges);
    int trim_left = max(0, ext_left - ext_max);
    int trim_right = max(0, ext_right - ext_max);
    len_extended_left = ext_left - trim_left;
    edge_extended = basevector(edges_to_bases.begin() + trim_left,
                          edges_to_bases.end() - trim_right);
}


void PrintReadSupports(const HyperBasevector& bubble_graph,
        const vecbasevector& bases,
        const vec<int> sample_ids,
        vec<tuple<int,int,read_place>> edge_rid_place) 
{
    if (edge_rid_place.empty()) return;
    vec<int> to_left; bubble_graph.ToLeft(to_left);

    auto CompareByEdgeSampleRead = 
        [&sample_ids](const tuple<int,int,read_place>& a, const tuple<int,int,read_place>&b)
        {   
            int e, rid, e2, rid2;
            tie(e, rid, ignore) = a;
            tie(e2, rid2, ignore) = b;
            int sample = sample_ids[rid];
            int sample2 = sample_ids[rid2];
            return tie(e, sample, rid) < tie(e2, sample2, rid2);
        };
    sort(edge_rid_place.begin(), edge_rid_place.end(),
            CompareByEdgeSampleRead);

    vec<int> shifts(edge_rid_place.size(), 0);
    for (size_t i = 0; i < edge_rid_place.size(); i++) {
        int e = -1, rid = -1;
        read_place rp;
        tie(e, rid, rp) = edge_rid_place[i];
        vec<int> rshifts;
        if (!rp.Fw()) rp.Reverse(bases[rid], bubble_graph);
        rp.FindReadShift(e, bubble_graph, rshifts);
        shifts[i] = rshifts[0];
    }

    int min_shift = *min_element(shifts.begin(), shifts.end());

    TextTable tb;
    for (size_t i = 0; i < edge_rid_place.size(); i++) {
        int e = get<0>(edge_rid_place[i]);
        size_t j = i+1;
        while (j < edge_rid_place.size() && 
                get<0>(edge_rid_place[j]) == e ) j++;

        tb.SetRawLine();
        tb << " ======================================== Edge " << e 
            << " ======================================= " << EndRow<TextTable>;

        // the reference segments +/- 100
        const int VarRefExtMax = 100;
        int node_left = to_left[e];
        int e_ref = bubble_graph.EdgeObjectIndexByIndexFrom(node_left,
                bubble_graph.FromSize(node_left)-1);
        basevector ref_seg;
        int ext_len = 0;
        ExtendedBasevector(bubble_graph, e_ref, min(VarRefExtMax, -min_shift), 
                ref_seg , ext_len);
        ForceAssertLe(ext_len, -min_shift);

        tb << "REF extended" << Tab<TextTable> << String(-ext_len-min_shift, ' ') << ref_seg << EndRow<TextTable>;
        // the variant 
        tb << "ALT" << Tab<TextTable> << String(-min_shift, ' ') << bubble_graph.EdgeObject(e) << EndRow<TextTable>;
        tb << "REF" << Tab<TextTable> << String(-min_shift, ' ') << bubble_graph.EdgeObject(e_ref) << EndRow<TextTable>;
        // the reads
        for (size_t k = i; k < j; k++) {
            int rid = get<1>(edge_rid_place[k]);
            const read_place& rp = get<2>(edge_rid_place[k]);
            basevector read = bases[rid];
            int sample = sample_ids[rid];
            if (!rp.Fw()) 
                read.ReverseComplement();
            tb << sample << ":" << rid << (rp.Fw() ? "+":"-" ) 
                << " q " << rp.Qsum()/1000 << Tab<TextTable>;
            tb << String(shifts[k] - min_shift, ' ') << read << EndRow<TextTable>;
        }
        i = j - 1;
    }
    tb.Print(cout, 1, "ll");
}

// Extension of variant edges to avoid the following read placement artificats:
//             C  
// ATATATATATAT   CGCGTAGT
//             ATC
// Read ATATATATATATCCGCGTAGT.. should support two edges equally well. However
// because we seed the read with 12mers , it will not thread through the bottom
// edge. This happens for tandem repeats. It can happen on both sides.
// 
// We fix the problem by extend the two edges on both side up to the base
// they diverge.
    
void ExtendTandemRepeatEdges(HyperBasevector& hb, int verbosity = 1) 
{
    ForceAssertEq(hb.K(), 1);
    for (int node = 0; node < hb.N(); node++) {
        if (hb.FromSize(node) <= 1) continue;
        if (hb.ToSize(node) != 1) continue;
        int node2 = hb.From(node)[0];
        for (int j = 0; j < hb.FromSize(node); j++) 
            ForceAssertEq(hb.From(node)[j], node2);
        if (hb.FromSize(node2) != 1) continue;

        vec<int> edges, lens;
        vec<basevector> bases;
        for (int j = 0; j < hb.FromSize(node); j++) {
            edges.push_back(hb.EdgeObjectIndexByIndexFrom(node,j));
            bases.push_back(hb.EdgeObjectByIndexFrom(node,j));
            lens.push_back(hb.EdgeObjectByIndexFrom(node,j).size());
        }
        SortSync(lens, bases, edges);
        if (lens.front() == lens.back()) continue;

        int head_edge = hb.EdgeObjectIndexByIndexTo(node, 0);
        basevector head = hb.EdgeObject(head_edge);

        int tail_edge = hb.EdgeObjectIndexByIndexFrom(node2, 0);
        basevector tail = hb.EdgeObject(tail_edge);

        int max_ext_left = 0, max_ext_right = 0;
        for (size_t i = 0; i < edges.size(); i++) {
            basevector short_edge = bases[i];
            for (size_t j = i+1; j < lens.size(); j++) {
                if (lens[j] == lens[i]) continue;

                int ext = lens[j] - lens[i];
                if ((int)head.size() > ext) {
                    basevector ext_base(head.end()-ext, head.end());
                    basevector short_extended = Cat(ext_base, short_edge);
                    if (short_extended == bases[j]) {
                        int p1 = head.size() -1 - ext, p2 = head.size() -1;
                        // leave at least one base in head
                        while (p1 > 0 && head[p1] == head[p2]) {
                            p1--, p2--;
                        }
                        int more_ext = head.size() - 1 - p2;
                        max_ext_left = max(max_ext_left, ext + more_ext);
                    }
                }
                if ((int)tail.size() > ext) {
                    basevector ext_base(tail.begin(), tail.begin()+ext);
                    basevector short_extended = Cat(short_edge, ext_base);
                    if (short_extended == bases[j]) {
                        int p1 = ext, p2 = 0;
                        // leave at least one base in tail
                        while (p1 < (int)tail.size()-1 && tail[p1] == tail[p2]) {
                            p1++, p2++;
                        }
                        max_ext_right = max(max_ext_right, ext + p2);
                    }
                }
            }
        }
        if (max_ext_right == 0 && max_ext_left == 0) continue;

        if (verbosity >= 1) {
            cout << "Extension for all edges between " << head_edge << " and " 
                << tail_edge << " are found " << endl;
            cout << "ext_left= " << max_ext_left << " ext_right= " << max_ext_right
                << " for edges ";
            edges.Println(cout);
        }
        if (max_ext_left > 0) 
            hb.EdgeObjectMutable(head_edge) = basevector(head.begin(), head.end() - max_ext_left);
        if (max_ext_right > 0)
            hb.EdgeObjectMutable(tail_edge) = basevector(tail.begin()+max_ext_right, tail.end());
        basevector left_add = basevector(head.end()-max_ext_left, head.end());
        basevector right_add = basevector(tail.begin(), tail.begin()+max_ext_right);
        for (size_t i = 0; i < edges.size(); i++) {
            hb.EdgeObjectMutable(edges[i]) =
                Cat(left_add, bases[i], right_add);
        }
    }
}

void ConnectLastTwoVertices(HyperBasevector& graph, map<int,vec<int>>& edge_to_varid,
        const basevector& edge, const vec<int>& vids) 
{
    graph.AddEdge(graph.N()-2, graph.N()-1, edge);
    edge_to_varid[graph.EdgeObjectCount()-1] = vids;
}

}

// Build a K=1 hyperbasevector of the reference with bubbles for each variant,
// or several variants if they are in the same ref position. Then calculate the
// probability of the variants by aligning reads back to it.
// The output is a map of each variant to a vector of pairs:
// {   {prob_ref, prob_var} .. from reads in sample0
//     {prob_ref, prob_var} .. from reads in sample1
//     .... }
void FindVariantProb(const ReadOriginTracker* p_read_tracker,
        const vec<VariantCallGroup>& vcall_groups,
        const vecbasevector& Gplus, 
        const vec<int>& Gplus_ext, 
        map<Variant, vec<pair<double,double>>>& probs,
        const long_logging* logc) 
{
    int verbosity = logc->verb["REFTRACE_VARIANTS"];
    bool bSafeFindPlaces = logc->REFTRACE_VARIANTS_LIMIT_K1_EFFORT;
    bool bMinorPhaseReconstruction = logc->REFTRACE_VARIANTS_MINOR_PHASE_RECONSTRUCTION;
    bool DetectSingleEdge = logc->DETECT_VARIANT_ON_SINGLE_EDGES;

    vec<int> vids_to_show_supports;
    if (logc->SHOW_READS_ON_VARIANT != -1)
        vids_to_show_supports.push_back(logc->SHOW_READS_ON_VARIANT - 1);

    if (verbosity >= 1) 
        cout << Date() << ": Finding probabilities for each variant" << endl;

    const vecbasevector& bases = p_read_tracker->Reads();
    const vecqualvector& quals = p_read_tracker->Quals();
    int nsamples = p_read_tracker->getSampleList().size();

    map<Variant, set<pair<int,int>>> var_group_branch;
    for (size_t grpid = 0; grpid < vcall_groups.size(); grpid++) {
        const VariantCallGroup& group = vcall_groups[grpid];
        for (int branchid = 0; branchid < group.GetNBranch(); branchid++) {
            const vec<VariantCall>& vcalls = group.GetVariantCalls(branchid);
            for (const VariantCall& x: vcalls)
                var_group_branch[x.variant].insert(make_pair( grpid, branchid ));
        }
    }
    if (var_group_branch.empty()) return;

    vec<Variant> vars;
    for (auto& x: var_group_branch) vars.push_back(x.first);
    vec<Variant> vars_ext(vars); // in Gplus
    for (auto& x: vars_ext) { x.pos += Gplus_ext[x.gid]; }

    //collect homopolymers
    enum { HOMO_BASE, HOMO_FRONT, HOMO_BACK, HOMO_LEN, HOMO_VIDX};
    typedef std::tuple<unsigned char, uint64_t,uint64_t,uint64_t, uint64_t> homopolymer_t;
    vec<homopolymer_t> homo_collection;
    const uint64_t homopolymer_critical_length=3;
    for(size_t vid=0;vid<vars.size();++vid){
        const Variant& variant = vars[vid];
        const auto& ref = variant.ref;
        const auto& alt = variant.alt;
        const uint64_t v_gid  = variant.gid;
        const auto& genome = Gplus[v_gid];
        const uint64_t v_start= variant.pos + Gplus_ext[v_gid];

        if( ref.size() > alt.size() ){ // deletion
            //            ForceAssert( alt.size() == 1);
            //            ForceAssert( alt[0] == ref[0] );
            //            ForceAssert( ref[0] == Base::val2Char(genome[v_start]));
            unsigned char baseval = Base::char2Val(ref[1]);
            size_t cc=2; for( ; cc<ref.size() && ref[cc]==ref[cc-1]; ++cc) {};
            if( cc==ref.size()){
                uint64_t front=v_start+1;
                uint64_t back =v_start+1;
                for( ; back+1 < genome.size() && genome[back] == genome[back+1] ; ++back){ }
                for( ; front > 0 && genome[front] == genome[front-1] ; --front){ }
                const uint64_t nDeleted = ref.size()-alt.size();
                if( back - front +1 >= homopolymer_critical_length + nDeleted){
                    if(verbosity>0) std::cout << "deletion with homopolymer of " << Base::val2Char(baseval) << " from  " << back-front+1 << " to " << back-front + 1 - nDeleted << std::endl;
                    homo_collection.push_back( std::make_tuple(baseval,front,back,back-front+1-nDeleted,vid));
                }
            }
        }//deletion
        else if( ref.size() < alt.size() ){ //insertion
            //            ForceAssert( ref.size() == 1);
            //            ForceAssert( alt[0] == ref[0] );
            //            ForceAssert( ref[0] == Base::val2Char(genome[v_start]));
            unsigned char baseval = Base::char2Val(alt[1]);
            uint64_t nInserted=1;
            for(;nInserted+1<alt.size()&&alt[nInserted+1]==alt[nInserted];++nInserted){};
            if(nInserted+1==alt.size()){
                uint64_t front=std::numeric_limits<uint64_t>::max();
                uint64_t back=std::numeric_limits<uint64_t>::max();
                if( genome[v_start+1] == baseval){
                    front=v_start+1;
                    back=v_start+1;
                    for( ; back+1 < genome.size() && genome[back] == genome[back+1] ; ++back){ }
                }
                if( genome[v_start] == baseval){
                    front=v_start;
                    if(back==std::numeric_limits<uint64_t>::max()){
                        back=v_start;
                    }
                    for( ; front > 0 && genome[front] == genome[front-1] ; --front){ }
                }
                // std::cout << front << " " << back << std::endl;
                if(   front!=std::numeric_limits<uint64_t>::max()
                        && back!=std::numeric_limits<uint64_t>::max()
                        && back-front+1+nInserted >= homopolymer_critical_length
                  ){
                    homo_collection.push_back( std::make_tuple(baseval,front,back,back-front+1+nInserted,vid));
                    if(verbosity>0) std::cout << "insertion with homopolymer of " << Base::val2Char(baseval) << " from  " << back-front+1 << " to " << back-front + 1 + nInserted << std::endl;
                }
                else if( back == std::numeric_limits<uint64_t> :: max() && nInserted > homopolymer_critical_length){
                    homo_collection.push_back( std::make_tuple(baseval,v_start,v_start,nInserted,vid));
                    if(verbosity>0) std::cout << "insertion with homopolymer of " << Base::val2Char(baseval) << " from  " << 0 << " to " << nInserted << std::endl;
                }
            }
        }//insertion
    }//vid
    Sort(homo_collection);
    typedef vec<std::tuple<double,double,int,int>>  homopolymer_info_t;
    typedef std::unordered_map<uint64_t, homopolymer_info_t> homopolymer_log_t;
    homopolymer_log_t homopolymer_log;
    vec< vec<std::pair<uint64_t,uint64_t>> > homopolymer_groups;
    for( size_t cc=0;cc<homo_collection.size();){
        auto end = cc+1;
        homopolymer_groups.push_back(vec<std::pair<uint64_t,uint64_t> >());
        homopolymer_log[std::get<HOMO_VIDX>(homo_collection[cc])]=homopolymer_info_t();
        homopolymer_groups.back().push_back(make_pair(std::get<HOMO_VIDX>(homo_collection[cc]),cc) );
        if(verbosity>0){
            std::cout << "homopolymer of " << Base::val2Char(std::get<HOMO_BASE>(homo_collection[cc])) << " "
                << std::get<HOMO_FRONT> (homo_collection[cc]) << " "
                << std::get<HOMO_BACK> (homo_collection[cc]) - std::get<HOMO_FRONT> (homo_collection[cc]) + 1 << " "
                << std::get<HOMO_LEN> (homo_collection[cc]) << " "
                << std::get<HOMO_VIDX> (homo_collection[cc]) << " "
                << std::endl;
        }
        for(;   end < homo_collection.size()
                && std::get<HOMO_BASE>(homo_collection[end]) == std::get<HOMO_BASE>(homo_collection[cc])
                && std::get<HOMO_FRONT>(homo_collection[end]) == std::get<HOMO_FRONT>(homo_collection[cc])
                && std::get<HOMO_BACK>(homo_collection[end]) == std::get<HOMO_BACK>(homo_collection[cc])
                ;++end){
            homopolymer_log[std::get<HOMO_VIDX>(homo_collection[end])]=homopolymer_info_t();
            homopolymer_groups.back().push_back(make_pair(std::get<HOMO_VIDX>(homo_collection[end]),end) );
            if(verbosity > 0){
                std::cout << "homopolymer of " << Base::val2Char(std::get<HOMO_BASE>(homo_collection[end])) << " "
                    << std::get<HOMO_FRONT> (homo_collection[end]) << " "
                    << std::get<HOMO_BACK> (homo_collection[end]) - std::get<HOMO_FRONT> (homo_collection[end]) + 1 << " "
                    << std::get<HOMO_LEN> (homo_collection[end]) << " "
                    << std::get<HOMO_VIDX> (homo_collection[end]) << " "
                    << std::endl;
            }
        }
        if(verbosity>0)std::cout<<" ----- " << std::endl;
        cc=end;
    }

    vec<tuple<int,int,int,int,int>> overlapping_vars; 
    vec<int> vars_to_overlaps(vars_ext.size(), -1);
    for (size_t i = 0; i < vars_ext.size(); i++) {
        int g = vars_ext[i].gid;
        int start = vars_ext[i].pos;
        int end = start + vars_ext[i].ref.size();
        size_t i2 = i+1;
        while (i2 < vars_ext.size() && vars_ext[i2].gid == g && 
                vars_ext[i2].pos < end) {
            end = max(end, vars_ext[i2].pos + (int)vars_ext[i2].ref.size());
            start = min(start, vars_ext[i2].pos);
            i2++;
        }
        overlapping_vars.push(i, i2, g, start, end);
        for (size_t k = i; k < i2; k++) 
            vars_to_overlaps[k] = overlapping_vars.size()-1;
        i = i2 - 1;
    }

    vec<tuple<int,int,int,int,int>> clusters; 
    // Heuristics for cluster separation.
    // effectively MinClusterSep cannot be larger than  MaxClusterSpan
    const int MaxClusterSpan = 20;
    const int MinClusterSep = 30; 
    // Check if variants from i to j appears more than once
    auto HasDupVars = [&](int i, int j) {
        for (int k = i; k < j; k++) 
            if (var_group_branch[vars[k]].size() > 1) return true;
        return false;
    };

    for (size_t i = 0; i < overlapping_vars.size(); i++) {
        int index1, index2, g, start, end;
        tie(index1, index2, g, start, end) = overlapping_vars[i];

        bool has_dup = HasDupVars(index1, index2);
        size_t i2 = i+1;
        while (i2 < overlapping_vars.size() && !has_dup) {
            int next_index1, next_index2, next_g, next_end;
            tie(next_index1, next_index2, next_g, ignore, next_end) = overlapping_vars[i2];
            bool next_has_dup= HasDupVars(next_index1, next_index2);
            if (next_g == g && next_end < end + MinClusterSep 
                    && !next_has_dup && end - start < MaxClusterSpan) {
                index2 = next_index2;
                end = next_end;
                i2++;
            } else 
                break;
        }
        clusters.push(index1, index2, g, start, end);
        i = i2 - 1;
    }

    // regroup the variants within a cluster
    vec<vec<vec<int>>> cluster_paths(clusters.size());
    for (size_t i = 0; i < clusters.size(); i++) {
        int index1, index2, g, start, end;
        tie(index1, index2, g, start, end) = clusters[i];
        map<int,map<int,vec<int>>> grouped; // separate by groupid branchid
        for (int k = index1; k < index2; k++) 
            for (auto& x: var_group_branch[vars[k]])
                grouped[x.first][x.second].push_back(k);

        vec<vec<int>> expanded{vec<int>{}};
        for (auto it = grouped.begin(); it != grouped.end(); ++it) {
            const map<int, vec<int>>& branches = (*it).second;
            vec<vec<int>> expanded_new;
            expanded_new.reserve(2*expanded.size());
            for (auto it2 = branches.begin(); it2 != branches.end(); it2++) {
                const auto& tail=(*it2).second;
                vec<vec<int>> expanded_copy;
                expanded_copy.reserve(expanded.size()*2);
                for( const auto& head : expanded){
                    if(head.size()==0){
                        expanded_copy.push_back(tail);
                    }
                    else if(       vars_ext[head.back()].pos+int64_t(vars_ext[head.back()].ref.size())-1
                               <   vars_ext[tail.front()].pos
                            || vars_ext[head.back()].gid!=vars_ext[tail.front()].gid
                           ){
                        expanded_copy.push_back(head);
                        expanded_copy.back().append(tail);
                    }
                    else{
std::cout << "WARNING: imminent crashing/ill-defined behavior in current variant clustering approach, patching with extra bubbles." << std::endl;
//The current clustering approach assumes that one v-group's v-calls would span a locus which does not overlap with that of another v-group.
//This is not true as of Sept 4, 2013 and would lead to down-stream crashing/undefined behavior.
//Here is the point-of-no-return leading to those crashes, and the following is a patch,
// designed according to the no-change-of-results requirement established on Sept 4, 2013.
                        int64_t head_back=head.size()-1;
                        for(;   head_back>=0
                             &&    vars_ext[head[head_back]].pos+int64_t(vars_ext[head[head_back]].ref.size())-1
                                >= vars_ext[tail[0]].pos
                            ;--head_back
                           ){}
                        expanded_copy.push_back(head);
                        expanded_copy.back().resize(head_back+1);
                        expanded_copy.back().append(tail);

                        size_t tail_front=0;
                        for(;   tail_front<tail.size()
                             &&    vars_ext[head.back()].pos+int64_t(vars_ext[head.back()].ref.size())-1
                                >= vars_ext[tail[tail_front]].pos
                            ;++tail_front
                           ){}
                        expanded_copy.push_back(head);
                        expanded_copy.back().insert(expanded_copy.back().end()
                                                   ,tail.begin()+tail_front
                                                   ,tail.end());
                    }
                }//head
                expanded_new.append(expanded_copy);
            }
            expanded = expanded_new;
        }
        UniqueSort(expanded);
        cluster_paths[i] = expanded;
    }

    // create the edges
    vec<vec<basevector>> cluster_seq(clusters.size());// [cluster][branch]
    for (size_t i = 0; i < clusters.size(); i++) {
        int index1, index2, g, start, end;
        tie(index1, index2, g, start, end) = clusters[i];

        basevector ref_seg(Gplus[g], start, end - start);
        for (size_t j = 0; j < cluster_paths[i].size(); j++) {
            const vec<int>& path = cluster_paths[i][j];
            basevector seq;
            int prev_stop = 0;
            for (size_t k = 0; k < path.size(); k++) {
                int vid = path[k];
                int offset = vars_ext[vid].pos - start;
                basevector alt( vars_ext[vid].alt );
                basevector ref( vars_ext[vid].ref );
                ForceAssertGe(offset, 0);
                ForceAssertLe(prev_stop, end - start);
                // In the case an indel follows a sub, the variant start position
                // will be the same in the reference.
                if (offset < prev_stop) {
                    ForceAssertEq(offset, prev_stop-1);
                    ForceAssertNe(alt.size(), ref.size());
                    offset++;
                    alt.SetToSubOf(alt, 1, alt.size()-1);
                    ref.SetToSubOf(ref, 1, ref.size()-1);
                }
                seq.append( ref_seg.begin()+ prev_stop, ref_seg.begin()+ offset );
                seq.append( alt.begin(), alt.end() );
                prev_stop = offset + ref.size();
            }
            if (prev_stop != (int)ref_seg.size())
                seq.append(ref_seg.begin() + prev_stop, ref_seg.end());
            cluster_seq[i].push_back(seq);
        }
    }

    HyperBasevector variant_graph(1);
    map<int,vec<int>>  edge_to_varid;
    for (size_t cid = 0; cid < clusters.size(); cid++) {
        int prev_g = (cid == 0 ? -1 : get<2>(clusters[cid-1]));
        int next_g = (cid == clusters.size()-1 ? -1 : get<2>(clusters[cid+1]));
        int g, start, end;
        tie(ignore, ignore, g, start, end) = clusters[cid];

        // add preceding edge if necessary
        if (g != prev_g) {
            variant_graph.AddVertices(1);
            if (start > 0) { 
                basevector ref_edge(Gplus[g], 0, start);
                variant_graph.AddVertices(1);
                ConnectLastTwoVertices(variant_graph, edge_to_varid, ref_edge, vec<int>());
            }
        } 

        // Add the variant edge, also add the reference
        variant_graph.AddVertices(1);
        for (size_t j = 0; j < cluster_seq[cid].size(); j++) {
            ConnectLastTwoVertices(variant_graph, edge_to_varid, 
                    cluster_seq[cid][j], cluster_paths[cid][j]);
        }
        basevector ref_edge(Gplus[g], start, end - start);
        ConnectLastTwoVertices(variant_graph, edge_to_varid, ref_edge, vec<int>());

        // add the ref edge between variant clusters in the same chromosome
        if (g == next_g && get<3>(clusters[cid+1]) > end) {
            basevector ref_edge(Gplus[g],  end, get<3>(clusters[cid+1]) - end);
            variant_graph.AddVertices(1);
            ConnectLastTwoVertices(variant_graph, edge_to_varid, ref_edge, vec<int>());
        }

        // add trailing edge if necessary
        if ( g != next_g && end < (int)Gplus[g].size() ) {
            basevector ref_edge(Gplus[g].begin()+end,  Gplus[g].end());
            variant_graph.AddVertices(1);
            ConnectLastTwoVertices(variant_graph, edge_to_varid, ref_edge, vec<int>());
        }
    }
    map<int,vec<int>> edge_to_varid_plus(edge_to_varid);
    for (auto& x: edge_to_varid_plus) 
        for_each(x.second.begin(), x.second.end(), [](int& x){x++;});

    if (verbosity >= 1) 
        DumpGraph("vgroup.dot", variant_graph);

    ExtendTandemRepeatEdges(variant_graph, verbosity);

    if (verbosity >= 3)
        for (int e = 0; e < variant_graph.EdgeObjectCount(); e++) {
            cout << "e= " << e << " vid= " << edge_to_varid[e]
                << " : " << variant_graph.EdgeObject(e) << endl;
        }
    
    
    vec<int> homopolymer_edges,bayo_tmp;
    for (int e = 0; e < variant_graph.EdgeObjectCount(); e++) {
        for( auto vid: edge_to_varid[e]){
            if( homopolymer_log.find(vid)!=homopolymer_log.end()){
                homopolymer_edges.push_back(e);
                break;
            }
        }
    }
    Sort(homopolymer_edges);

    vec<int> edges_to_show_supports(vids_to_show_supports.size());
    if (!vids_to_show_supports.empty()) {
        for (size_t i = 0; i < vids_to_show_supports.size(); i++) 
            for (int e = 0; e < variant_graph.EdgeObjectCount(); e++)
                if (Member(edge_to_varid[e], vids_to_show_supports[i]))
                    edges_to_show_supports[i] = e;
        cout << "Looking for support for the following variants" << endl;
        for (size_t i = 0; i < edges_to_show_supports.size(); i++) {
            cout << "vid= " << vids_to_show_supports[i] 
                << " edge= " << edges_to_show_supports[i] << endl;
        }
        // Always add the ref edge for each variant
        vec<int> edges_to_show_supports2;
        vec<int> to_left; variant_graph.ToLeft(to_left);
        for (int e: edges_to_show_supports) {
            int node_left = to_left[e];
            edges_to_show_supports2.push_back(e);
            for (int i = 0; i < variant_graph.FromSize(node_left); i++)
                edges_to_show_supports2.push_back(variant_graph.
                        EdgeObjectIndexByIndexFrom(node_left,i));
        }
        swap(edges_to_show_supports2, edges_to_show_supports);
    }

    vec< vec< pair<int,int> > > home_index; // pair of (rid, qsum)
    vec<tuple<int,int,read_place>> edge_rid_place;
    FindReadHomesBest(bases, quals, variant_graph, &home_index,
            edges_to_show_supports, &edge_rid_place,
            verbosity,bSafeFindPlaces);

    if (!edge_rid_place.empty()) {
        vec<int> sample_ids(bases.size(), -1);
        for (size_t i = 0; i < sample_ids.size(); i++) 
            sample_ids[i] = p_read_tracker->getSampleID(i);
        PrintReadSupports(variant_graph, bases, sample_ids, edge_rid_place);
    }

    map<int, vec<pair<uint64_t,uint64_t>>> read_support_count;
    map<int, vec<tuple<int,int,int,int>>> read_support_by_strand;
    for (int n1 = 0; n1 < variant_graph.N(); n1++) {
        if (variant_graph.From(n1).size() <= 1) continue;
        vec<int> branches(variant_graph.From(n1).size());
        for (size_t j = 0; j < variant_graph.From(n1).size(); j++) 
            branches[j] = variant_graph.EdgeObjectIndexByIndexFrom(n1,j);
        
        bool bHasHomoPolymer=false;
        size_t nbranches = branches.size();
        const int qgood = 4;
        // the reference branch is the last one
        vec<vec<double>> qsum(nbranches, vec<double>(nsamples, 0.0));
        vec<vec<int>> nreads(nbranches, vec<int>(nsamples, 0));
        // for + and - read dir
        vec<vec<pair<int,int>>> nreads_by_strand(nbranches, vec<pair<int,int>>(nsamples));
        for (size_t j = 0; j < branches.size(); j++) {
            int e = branches[j];
            if( Member(homopolymer_edges,e)) bHasHomoPolymer=true;
            for (auto& x: home_index[e]) {
                int rid = x.first;
                int sample_id = p_read_tracker->getSampleID(rid);
                ReadOriginTracker::READ_DIR rdir = p_read_tracker->ReadDirection(rid);
                if (x.second/1000 >= qgood && sample_id >= 0) {
                    qsum[j][sample_id] += x.second/1000;
                    nreads[j][sample_id] ++;
                    if (rdir == ReadOriginTracker::PLUS)
                        nreads_by_strand[j][sample_id].first++;
                    else if (rdir == ReadOriginTracker::MINUS)
                        nreads_by_strand[j][sample_id].second++;
                    if (verbosity >= 2) {
                        cout << "add read_" << rid << " (" 
                            << sample_id << ") to edge " << e 
                            << " qsum= " << x.second/1000 << endl;
                    }
                }
            }
        }
        if (verbosity >= 2) {
            cout << "branch: " << endl;
            for (size_t j = 0; j < branches.size(); j++) {
                int e = branches[j];
                if (j == branches.size() -1) 
                    cout << "reference ";
                else
                    cout << "vid = " << edge_to_varid_plus[e];
                cout << " e= " << e;
                for (int sample = 0; sample < nsamples; sample++)
                    cout << " nreads= " << nreads[j][sample] << " q= " << qsum[j][sample];
                cout  << endl;
            }
        }
        // Doe not call if there are nread <=2 && q <= 40 while all other
        // branches have 5 times ore more qsum.
        if(!bHasHomoPolymer){
            const int WeakSupport = 2, WeakQSum = 40;
            const int StrongQsumFactor = 5;
            for (int sample = 0; sample < nsamples; sample++) {
                vec<int> weak;
                double max_weak_support = 0;
                for (size_t j = 0; j < branches.size(); j++) {
                    if (nreads[j][sample] > 0 && nreads[j][sample] <= WeakSupport &&
                            qsum[j][sample] <= WeakQSum) {
                        weak.push_back(j);
                        max_weak_support = max(max_weak_support, qsum[j][sample]);
                    }
                }
                if (weak.empty()) continue;
                size_t ngood = 0, nzero = 0;
                for (size_t j = 0; j < branches.size(); j++) 
                    if (qsum[j][sample] >= StrongQsumFactor * max_weak_support)
                        ngood++;
                    else if (nreads[j][sample] == 0) 
                        nzero++;
                if (weak.size() + ngood + nzero == branches.size() && ngood > 0) {
                    for (int j: weak) {
                        qsum[j][sample] = 0.;
                        if (verbosity >= 2 && nreads[j][sample] != 0)
                            cout << "Disable branch " << j << " for sample " << sample
                                << " vid = " << edge_to_varid_plus[branches[j]]
                                << " nread= " << nreads[j][sample]
                                << " qsum= " << qsum[j][sample]
                                << endl;
                    }
                }
            }
        }
        // Find the qsum of reads that support the reference for given variant
        // located in branch j.
        auto QsumRef = [&](int vid, int sample) {
            double qsum_var = 0;
            double qsum_ref = 0;
            int nreads_var = 0, nreads_ref = 0;
            int nreads_var_plus = 0, nreads_var_minus = 0,
                nreads_ref_plus = 0, nreads_ref_minus = 0;
            for (int k = 0; k < branches.isize(); k++) {
                bool overlap = false;
                bool hasvid = false;
                for (int vid2: edge_to_varid[branches[k]]) 
                    if (vid2 == vid) {
                        hasvid = true; 
                        break;
                    } else if (vars_to_overlaps[vid2] == vars_to_overlaps[vid]) {
                        overlap = true;
                        break;
                    }
                if (hasvid) {
                    qsum_var +=  qsum[k][sample];
                    nreads_var += nreads[k][sample];
                    nreads_var_plus +=  nreads_by_strand[k][sample].first;
                    nreads_var_minus +=  nreads_by_strand[k][sample].second;
                }
                else if (!overlap) {
                    qsum_ref += qsum[k][sample];
                    nreads_ref += nreads[k][sample];
                    nreads_ref_plus +=  nreads_by_strand[k][sample].first;
                    nreads_ref_minus +=  nreads_by_strand[k][sample].second;
                }
            }
            return make_tuple(qsum_ref, qsum_var, nreads_ref, nreads_var,
                    nreads_ref_plus, nreads_ref_minus,
                    nreads_var_plus, nreads_var_minus);
        };
        set<int> vids;
        for (size_t j = 0; j < branches.size()-1; j++)
            for (int vid: edge_to_varid[branches[j]])
                vids.insert(vid);
        for (int vid: vids) {
            probs[vars[vid]].resize(nsamples);
            read_support_count[vid].resize(nsamples);
            read_support_by_strand[vid].resize(nsamples);
            for (int sample = 0; sample < nsamples; sample++) {
                auto ans = QsumRef(vid, sample);
//std::cout << "bayo: vid=" <<vid << " sample="<<sample << " " << get<0>(ans) << " " << get<1>(ans) << " " << get<2>(ans) << " " << get<3>(ans) << " " << get<4>(ans) << " " << get<5>(ans) << " " << get<6>(ans) << " " << get<7>(ans) << std::endl;
                probs[vars[vid]][sample].first += get<0>(ans);
                probs[vars[vid]][sample].second += get<1>(ans);
                read_support_count[vid][sample].first += get<2>(ans);
                read_support_count[vid][sample].second += get<3>(ans);
                get<0>(read_support_by_strand[vid][sample]) += get<4>(ans);
                get<1>(read_support_by_strand[vid][sample]) += get<5>(ans);
                get<2>(read_support_by_strand[vid][sample]) += get<6>(ans);
                get<3>(read_support_by_strand[vid][sample]) += get<7>(ans);
            }
            auto itr=homopolymer_log.find(vid);
            if( itr!=homopolymer_log.end()){
                auto& entry = (*itr).second;
                entry.resize(nsamples);
                for (int sample = 0; sample < nsamples; sample++) {
                    auto ans = QsumRef(vid, sample);
                    entry[sample] = make_tuple(get<0>(ans), get<1>(ans),
                           get<2>(ans), get<3>(ans));
                }
            }
        }
    }
    for ( const auto& group: homopolymer_groups){
        if( group.size()==0) continue;
        const auto first_indices = group.front();
        const unsigned char baseval = std::get<HOMO_BASE> (homo_collection[first_indices.second]);
        const auto reference_length = std::get<HOMO_BACK> (homo_collection[first_indices.second]) - std::get<HOMO_FRONT> (homo_collection[first_indices.second]) + 1;

//        if(verbosity>0)std::cout << "group of " << std::get<HOMO_FRONT> (homo_collection[first_indices.second]) << std::endl;

        for(int sample=0;sample<nsamples;sample++){
            std::unordered_map< uint64_t , vec<std::tuple<uint64_t,double,int>>> length_lookup;


            length_lookup[reference_length].push_back( make_tuple(std::numeric_limits<uint64_t>::max()
                        ,std::get<0>(homopolymer_log[first_indices.first][sample])
                        ,std::get<2>(homopolymer_log[first_indices.first][sample])
                        ));
            for( const auto& entry: group){
                auto length = std::get<HOMO_LEN> (homo_collection[entry.second]);
                auto vid = entry.first;
                if( std::get<2>(length_lookup[reference_length][0]) >  std::get<2>(homopolymer_log[vid][sample]) ){
                    length_lookup[reference_length][0]= make_tuple(vid
                            ,std::get<0>(homopolymer_log[vid][sample])
                            ,std::get<2>(homopolymer_log[vid][sample])
                            );
                }
                length_lookup[length].push_back( make_tuple(vid
                            ,std::get<1>(homopolymer_log[vid][sample])
                            ,std::get<3>(homopolymer_log[vid][sample])
                            ));
            }

            AdjustHomopolymerProb(baseval,length_lookup);


            double ref_prob = std::get<1>(length_lookup[reference_length][0]);

            for( const auto&len_entry: length_lookup){
                for(const auto& var_entry: len_entry.second){
                    const auto& vid  = std::get<0>(var_entry);
                    const auto& prob = std::get<1>(var_entry);
                    if( vid != std::numeric_limits<uint64_t>::max()){
                        probs[vars[vid]][sample]=make_pair(ref_prob,prob);
                    }
                }
            }
//            if(verbosity>0) std::cout << "--end of sample-- " << std::endl;
        }
//        if(verbosity>0) std::cout << "--end of group-- " << std::endl;
    }

    // set q=3 for one read supported variant
    for( const auto& entry: read_support_count){
        const auto& vid = entry.first;
        const vec<pair<uint64_t,uint64_t>> & log=entry.second;
        for(size_t sample=0;sample<log.size();++sample){
            if( log[sample].first ==1){
                if( probs[vars[vid]][sample].first > 0.){
                probs[vars[vid]][sample].first = min( 3.0 , probs[vars[vid]][sample].first );
                if(verbosity>0) std::cout<<"WARNING: 1-read filter: zero'ing P-ref of sample "<<sample << " of variant with vid " << vid << " REF=" << vars[vid].ref << " ALT="<<vars[vid].alt<< std::endl;
                }
            }
            if( log[sample].second ==1){
                if( probs[vars[vid]][sample].second > 0.){
                probs[vars[vid]][sample].second = min( 3.0 , probs[vars[vid]][sample].second );
                if(verbosity>0) std::cout<<"WARNING: 1-read filter: zero'ing P-alt of sample "<<sample << " of variant with vid " << vid << " REF=" << vars[vid].ref << " ALT="<<vars[vid].alt<< std::endl;
                }
            }
        }
    }

    // print strand support 
    if (verbosity >= 1){
        // only variants on double-banch bubble are considered
        set<Variant> var_on_double_edges;
        for (size_t grpid = 0; grpid < vcall_groups.size(); grpid++) {
            const VariantCallGroup& group = vcall_groups[grpid];
            if (group.GetNBranch() != 2) continue;
            for (int branchid = 0; branchid < group.GetNBranch(); branchid++) {
                const vec<VariantCall>& vcalls = group.GetVariantCalls(branchid);
                for (const VariantCall& x: vcalls)
                    var_on_double_edges.insert(x.variant);
            }
        }
        // remove variant appearing in both branches or overlap with other
        // variants
        auto it = var_on_double_edges.begin();
        for (; it != var_on_double_edges.end(); ) {
            bool appear_both_branch = var_group_branch[*it].size() > 1;
            int vid = find(vars.begin(), vars.end(), *it) - vars.begin();
            auto overlap = overlapping_vars[vars_to_overlaps[vid]];
            bool overlap_with_other = (get<1>(overlap) - get<0>(overlap)) > 1;
            if (appear_both_branch || overlap_with_other) 
                var_on_double_edges.erase(it++);
            else
                ++it;
        }

        TextTable tb;
        tb << "VID" << Tab<TextTable> << "+Ref" << Tab<TextTable> << "-Ref" << Tab<TextTable> << "+Var" << Tab<TextTable> 
            << "-Var" << Tab<TextTable> << "Strand_Bias" << EndRow<TextTable>;
        tb << DoubleLine;
        for (size_t vid = 0; vid < vars.size(); vid++) {
            Variant var = vars[vid];
            if (var_on_double_edges.find(var) == var_on_double_edges.end())
                continue;
            tb << vid+1 << Tab<TextTable>;
            for(int sample = 0; sample < nsamples; sample++) {
                if (sample != 0) tb << Tab<TextTable>;
                int var_plus, var_minus, ref_plus, ref_minus;
                tie(ref_plus, ref_minus,var_plus, var_minus)
                    = read_support_by_strand[vid][sample];
                tb << ref_plus << Tab<TextTable>;  //a
                tb << ref_minus << Tab<TextTable>; //c
                tb << var_plus << Tab<TextTable>;  //b
                tb << var_minus << Tab<TextTable>;        //d
                double a = ref_plus, c = ref_minus, b = var_plus, d = var_minus;
                if (a+c < b+d) {
                    swap(a,b);
                    swap(c,d);
                }
                double sb = (b/(a+b) - d/(c+d))/((b+d)/(a+b+c+d));
                tb << sb << Tab<TextTable>;
            }
            tb << EndRow<TextTable>;
        }
        ofstream ofs("strand_bias.txt");
        tb.Print(ofs, 1, "rrrrr");
    }

    if (DetectSingleEdge) {
        const double QsumRefSingleEdge = 0;
        set<Variant> var_on_single_edges;
        for (size_t grpid = 0; grpid < vcall_groups.size(); grpid++) {
            const VariantCallGroup& group = vcall_groups[grpid];
            if (group.GetNBranch() > 1) continue;
            for (int branchid = 0; branchid < group.GetNBranch(); branchid++) {
                const vec<VariantCall>& vcalls = group.GetVariantCalls(branchid);
                for (const VariantCall& x: vcalls)
                    var_on_single_edges.insert(x.variant);
            }
        }
        if (DetectSingleEdge) {
            cout << "Overiding probability for " << var_on_single_edges.size()
                << " variants found on single edges" << endl;
        }
        for (auto it = probs.begin(); it != probs.end(); ++it) {
            if (var_on_single_edges.find(it->first) != var_on_single_edges.end())
                for (auto& p: it->second)
                    p.first = QsumRefSingleEdge;
        }
    }

    if( bMinorPhaseReconstruction )
    {
        RemoveDanglingCalls( probs, vcall_groups);
        FillInZeroProbAccordingToGroup( probs, vcall_groups);
//std::cout << "----------------------------------------------------------------" << std::endl;
//for (size_t ii=106 ; ii<109 ; ++ii) for( const auto& entry: probs[vars[ii]]){ std::cout << ii << " " << entry.first << " " << entry.second << std::endl; }
    }
}


// Performe gapless alignment of reads to hyperbasevector and report the
// preference of reads aligned to each edge as calculated by qualsum difference
// between best in-edge placement and best out-edge placement.
void FindReadHomesBest(const vecbasevector& bases, const vecqualvector& quals, 
        const HyperBasevector& bubble_graph,
        vec< vec< pair<int,int> > >* homes_index,
        vec<int> edges_to_show_supports,
        vec<tuple<int,int,read_place>>* edge_rid_place,
        int verbosity, const bool bSafeFindPlaces)
{
    if (verbosity >= 1)
        cout << Date() << ": Finding placements of " << bases.size() << " reads " << endl;

    ForceAssertEq(bases.size(), quals.size());
    const int L = 12;
    const int infinity = 1000000000;

    vec<int> to_left; bubble_graph.ToLeft(to_left);
    vec<int> to_right; bubble_graph.ToRight(to_right);

    HyperBasevector hb_fw(bubble_graph), hb_rc(bubble_graph);
    hb_rc.Reverse( );
    vec<int> to_right_fw, to_right_rc;
    hb_fw.ToRight(to_right_fw), hb_rc.ToRight(to_right_rc);
    vecbasevector x_fw, x_rc;
    for ( int i = 0; i < hb_fw.EdgeObjectCount( ); i++ )
        x_fw.push_back( hb_fw.EdgeObject(i) );
    for ( int i = 0; i < hb_rc.EdgeObjectCount( ); i++ )
        x_rc.push_back( hb_rc.EdgeObject(i) );
    VecIntPairVec locs_fw, locs_rc;
    CreateGlocs( x_fw, L, locs_fw );
    CreateGlocs( x_rc, L, locs_rc );

    vec< vec<read_place> > PLACES2( bases.size( ) );
    const int min_qual = 1;
    const double prox = 50 * 1000;
    uint64_t nTerminated=0;
    if( verbosity >= 1){
        if( bSafeFindPlaces ){
            std::cout << "FindPlaces being safe" << std::endl;
        }
        else{
            std::cout << "FindPlaces being unsafe" << std::endl;
        }
    }

    #pragma omp parallel for schedule (dynamic, 1) reduction(+:nTerminated)
    for (size_t id = 0; id < bases.size(); id++) {
        if ( bases[id].isize( ) < L ) continue;
        int n = KmerId( bases[id], L, 0 );
        int qual_sum = infinity;
        if( bSafeFindPlaces ){
            const uint64_t maxNSteps=5000000;
            auto tmp = SafeFindPlaces( bases[id], quals[id], n, hb_fw, hb_rc, to_right_fw,
                                     to_right_rc, locs_fw, locs_rc, PLACES2[id], qual_sum,
                                     maxNSteps,
                                     min_qual, prox );
            if(tmp>0){
                ++nTerminated;
            }
        }
        else{
            FindPlaces( bases[id], quals[id], n, hb_fw, hb_rc, to_right_fw,
               to_right_rc, locs_fw, locs_rc, PLACES2[id], qual_sum, min_qual, prox );
        }
    }

    if( nTerminated > 0){
/*
        std::cout << "\nWARNING: " + ToString(nTerminated) + "/" + ToString(bases.size())
                  + " reads have not been aligned to the K=1 bubble graph to avoid exponentially long run time.\n" << std::endl;
*/
    }

    // Find out the best placement for each edge
    if (verbosity >= 1)
        cout << Date() << ": Finding the best placement for each read " << endl;

    vec<int> best_qsums(bubble_graph.EdgeObjectCount(), infinity);
    #pragma omp parallel
    {
        vec<int> best_qsums_loc(best_qsums);
        // NOTE(tim): possible intel compiler bug - putting PLACES2.size()
        // inside the for-statement causes an "internal error: null pointer"
        // on compilation
        size_t places2size = PLACES2.size();
        #pragma omp barrier
        #pragma omp for schedule (dynamic, 1)
        for (size_t id = 0; id < places2size; id++)
            for (int j = 0; j < PLACES2[id].isize(); j++)
                for (auto e: PLACES2[id][j].E())
                    best_qsums_loc[e] = min(best_qsums_loc[e], PLACES2[id][j].Qsum());
        #pragma omp critical
        {
            for( size_t edge=0;edge<best_qsums.size();++edge){
                best_qsums[edge] = min(best_qsums_loc[edge], best_qsums[edge]);
            }
        }
    }

    if (verbosity>=2) {
        const int MaxPlacementToDisplay = 1000;
        for (size_t id = 0; id < bases.size(); id++ ) {   
            cout << "read_" << id << endl;
            for ( int j = 0; j < PLACES2[id].isize( ); j++ ) {
                if (j >= MaxPlacementToDisplay) {
                    cout << "... " << PLACES2[id].isize() - j 
                        << " more ignored ..." << endl;
                    break;
                }
                vec<int> v;
                for ( int l = 0; l < PLACES2[id][j].N( ); l++ )
                    v.push_back( PLACES2[id][j].E(l) );
                cout << "    p=" << PLACES2[id][j].P() 
                    << " q= " <<  PLACES2[id][j].Qsum() << " ";
                v.Println(cout);
            }
        }
    }

    // Find the home of the reads, and quality score difference between the
    // best placement on this edge or in other edge
    // { { (edge11, diff11), (edge12, diff12), ...} // for read1
    //   { (edge21, diff21), (edge22, diff22), ...} // for read2
    //   ...}

    if (verbosity >= 1)
        cout << Date() << ": Finding the read support and qsum difference " << endl;

    vec< vec< pair<int,double>>> homes( bases.size( ) );
    vec<int> ids(bases.size(), vec<int>::IDENTITY);
    vec<int> nplaces(bases.size(), 0);
    for (size_t i = 0; i < PLACES2.size(); i++) 
        nplaces[i] = PLACES2[i].size();
    ReverseSortSync(nplaces, ids);
    
    #pragma omp parallel for schedule (dynamic, 1)
    for (size_t idx = 0; idx < ids.size(); idx++ ) {   
        int id = ids[idx];
        vec<int> v;
        for ( int j = 0; j < PLACES2[id].isize( ); j++ )
            for ( int l = 0; l < PLACES2[id][j].N( ); l++ )
                v.push_back( PLACES2[id][j].E(l) );
        UniqueSort(v);
        for ( int i = 0; i < v.isize( ); i++ ) {
            int e = v[i];
            bool to_show_read_support = Member(edges_to_show_supports, e);
            // what are the sister branches?
            vec<int> e2s;
            int n_left = to_left[e];
            for (size_t k = 0; k < bubble_graph.From(n_left).size(); k++) {
                int e2 = bubble_graph.EdgeObjectIndexByIndexFrom(n_left, k);
                e2s.push_back(e2);
            }
            int group_best_qsum = best_qsums[e];
            for (int e2: e2s) group_best_qsum = min(group_best_qsum, best_qsums[e2]);

            int q1 = infinity, q2 = infinity;
            int best_placement = -1;
            for ( int j = 0; j < PLACES2[id].isize( ); j++ ) {
                int qsum = PLACES2[id][j].Qsum();
                for (int e2: e2s) {
                    if (qsum > group_best_qsum + prox) continue;
                    if ( Member(PLACES2[id][j].E( ), e2) ) {
                        if (e2 == e) {
                            if (qsum < q1) {
                                q1 = qsum;
                                if (to_show_read_support) 
                                    best_placement = j;
                            }
                        }
                        else
                            q2 = Min( q2, qsum );
                    }
                }
            }
            if (q1 == infinity) continue;
            if (q2 == infinity) q2 = q1 + 30000; // +30 for single placement
            if (q1 < q2) {
                //double dq = q2 - q1;
                //double prob2 = pow(10.0, -dq/1e4) + 0.25* q1/dq/bases[id].size();
                //prob2 = min(1.0, prob2);
                //double dq2 = -log(prob2)/log(10.0)*1e4;
                //homes[id].push( e, dq2 );    
                //#pragma omp critical
                //cout << "dq= " << q2 - q1 << " dq2= " << dq2 
                //    << " prob2= " << prob2 
                //    << " q1= " << q1 << " q2= " << q2 << endl;
                homes[id].push( e, q2 - q1 );    
                if (best_placement != -1) {
                    #pragma omp critical
                    edge_rid_place->push(e, id, PLACES2[id][best_placement]);
                }
            }
        }    
    }

    homes_index->clear_and_resize(bubble_graph.EdgeObjectCount());
    for (size_t id = 0; id < bases.size(); id++)
        for (int j = 0; j < homes[id].isize( ); j++)
            (*homes_index)[ homes[id][j].first ].push(id, homes[id][j].second);

    if (verbosity>=2) {
        vec<int> nreads( bubble_graph.EdgeObjectCount(), 0);
        for (size_t id = 0; id < bases.size(); id++ ) 
            for ( int j = 0; j < PLACES2[id].isize( ); j++ )
                for ( int l = 0; l < PLACES2[id][j].N( ); l++ )
                    nreads[ PLACES2[id][j].E(l) ]++;
        cout << "read support for each edges " << endl;
        for (size_t e = 0; e < nreads.size(); e++) {
            cout << "edge" << e << " nreads=" << nreads[e] << ": ";
            for (auto& x: (*homes_index)[e])
                cout << x.first << "(" << x.second << ")" << " ";
            cout << endl;
        }
    }
}

