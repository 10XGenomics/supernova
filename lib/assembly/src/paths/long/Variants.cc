///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS
#include "paths/long/Variants.h"

#include "CoreTools.h"
#include "paths/HyperEfasta.h"
#include "paths/long/CleanEfasta.h"
#include "paths/long/Logging.h"
#include "paths/long/SupportedHyperBasevector.h"
#include "util/TextTable.h"
#include "paths/long/ReadOriginTracker.h"
#include "paths/long/VariantCallTools.h"
#include "paths/long/VariantReadSupport.h"
#include "paths/long/VariantFilters.h"
#include "paths/long/VariantPostProcess.h"

int MarkVariants( HyperEfasta& he, const vec<VariantSignature>& v_signatures,
     const long_logging& logc )
{    double clock = WallClockTime( );
     const int verbosity = 0;
     int nfound = 0;
     vec<efasta> snp_bubbles_str( v_signatures.size( ) );

     typedef efasta::size_type st_t;
     typedef std::tuple<st_t,st_t,st_t,st_t> bubble_signature_t;
     std::vector<bubble_signature_t> bubble_signatures(v_signatures.size());
     vec<size_t> bubbles_activated;

     typedef vec<std::pair<size_t,size_t> >edit_t;
     vec<edit_t> edits(he.EdgeObjectCount());

//     vec<efasta> old_new_edges(he.EdgeObjectCount());

     #pragma omp parallel for schedule(dynamic,1)
     for ( size_t k = 0; k < v_signatures.size( ); ++k )
     {    snp_bubbles_str[k] = efasta( v_signatures[k].Branches() );

          st_t one,two,three;
          one = snp_bubbles_str[k].find( '{' , 0);
          ForceAssert( one!=efasta::npos );
          two = snp_bubbles_str[k].find( ',' , one+1);
          ForceAssert( two!=efasta::npos );
          three = snp_bubbles_str[k].find( '}' , two+1);
          ForceAssert( three!=efasta::npos );
          bubble_signatures[k] = std::make_tuple(one,v_signatures.size()-k,two-one,three-one);

     }
     std::sort(bubble_signatures.begin(),bubble_signatures.end(),std::greater<bubble_signature_t>());

     #pragma omp parallel
     {
     vec<size_t> bubbles_activated_loc;
     #pragma omp for schedule(dynamic,1) reduction(+:nfound)
     for ( int i = 0; i < he.EdgeObjectCount( ); ++i )
     {    const efasta& edge = he.EdgeObject(i);

          st_t start=0;
          for( st_t open_pos=edge.find('{',start); open_pos!=efasta::npos; open_pos=edge.find('{',start)){
              size_t k;
              for( k=0; k < bubble_signatures.size() ; ++k){
                  if(   open_pos >= std::get<0>(bubble_signatures[k])
                     && open_pos + std::get<2>(bubble_signatures[k]) < edge.size()
                     && edge[open_pos + std::get<2>(bubble_signatures[k])] == ','
                     && open_pos + std::get<3>(bubble_signatures[k]) < edge.size()
                     && edge[open_pos + std::get<3>(bubble_signatures[k])] == '}'
                     && edge.Contains( snp_bubbles_str[v_signatures.size()-std::get<1>(bubble_signatures[k])]
                                     , open_pos - std::get<0>(bubble_signatures[k])
                                     )
                    ){
                    break;
                  }
              }
              if( k==bubble_signatures.size()){
                  auto end_pos= open_pos+1 ;
                  edits[i].push_back(std::make_pair(start,end_pos));
                  start=end_pos;
              }
              else{
                  auto v_index = v_signatures.size()-std::get<1>(bubble_signatures[k]);
                  auto v_start = open_pos - std::get<0>(bubble_signatures[k]);
                  // this if-statement enforces behavior identical to the old code, namely,
                  // if the tail of one bubble overlaps with the head of another, the later won't be marked
                  if( v_start >= start){

                      edits[i].push_back(std::make_pair(start,v_start));
                      edits[i].push_back(std::make_pair(std::numeric_limits<size_t>::max(),v_index));
                      bubbles_activated_loc.push_back(v_index);

                      nfound++;
                      if ( verbosity >=1 ) {
                           #pragma omp critical
                           {    cout << "jump: marked variant " << v_index << " in efasta edge "
                                     << i << endl;    }    }

                      start = v_start + snp_bubbles_str[v_index].size();
                  }
                  else{
                      auto end_pos= open_pos+1 ;
                      edits[i].push_back(std::make_pair(start,end_pos));
                      start=end_pos;
                  }
              }
              if( start >= edge.size()) break;
          }
          if( start < edge.size()){
              edits[i].push_back(std::make_pair(start,edge.size()));
          }
#if 0
          efasta& new_edge=old_new_edges[i];
          for ( int pos = 0; pos < edge.isize( ); pos++ )
          {    int k;
               for ( k = 0; k < v_signatures.isize( ); ++k )
                    if ( edge.Contains( snp_bubbles_str[k], pos ) ) break;
               if ( k == v_signatures.isize( ) )
               {    new_edge.push_back( edge[pos] );
                    continue;    }
               efasta branch_combined(
                    v_signatures[k].Branch1(), v_signatures[k].Branch2() );
               efasta branch_combined_tagged;
               for( auto seq : branch_combined )
               {    if ( seq == '{' ) branch_combined_tagged.append("{[Variant]");
                    else branch_combined_tagged.push_back(seq);    }
               new_edge += branch_combined_tagged;
               pos += snp_bubbles_str[k].size( ) - 1;
//               nfound++;
               if ( verbosity >=1 )
               {
                    #pragma omp critical
                    {    cout << "marked variant " << k << " in efasta edge "
                              << i << endl;    }    }    }

//          edge = new_edge;
#endif
     }
     #pragma omp critical
     {
          bubbles_activated.append(bubbles_activated_loc);
     }
     }
     vec<efasta> snp_bubbles_combined_tagged( v_signatures.size( ) );
     UniqueSort(bubbles_activated);
     const size_t nActivated=bubbles_activated.size();
     #pragma omp parallel for schedule(dynamic,1)
     for(size_t ii=0;ii<nActivated;++ii) {
         auto v_index=bubbles_activated[ii];
         efasta loc( v_signatures[v_index].Branch1( ), v_signatures[v_index].Branch2( ),"[Variant]" );
         swap( loc , snp_bubbles_combined_tagged[v_index]);
     }

     #pragma omp parallel for schedule(dynamic,1)
     for ( int i = 0; i < he.EdgeObjectCount( ); ++i ) {
          const efasta& old_edge = he.EdgeObject(i);
          size_t edge_size=0;
          for( const auto& entry: edits[i]){
              if( entry.first== std::numeric_limits<size_t>::max()){
                   auto v_index=entry.second;
                   ForceAssert( (snp_bubbles_combined_tagged[v_index].size()==0) == (snp_bubbles_str[v_index].size()==0));
                   edge_size+=snp_bubbles_combined_tagged[v_index].size();
              }
              else{
                   ForceAssert( entry.second <= old_edge.size());
                   edge_size+=entry.second-entry.first;
              }
          }
          efasta new_edge;
          new_edge.reserve(edge_size);
          for( const auto& entry: edits[i]){
              if( entry.first== std::numeric_limits<size_t>::max()){
                   new_edge+=snp_bubbles_combined_tagged[entry.second];
              }
              else{
                   new_edge.append(old_edge,entry.first,entry.second-entry.first);
              }
          }
          swap( new_edge , he.EdgeObjectMutable(i));
//          ForceAssert(old_new_edges[i]==he.EdgeObject(i));
     }

     REPORT_TIME( clock, "used marking variants" );
     if (logc.STATUS_LOGGING)
     {    cout << Date( ) << ": marked " << nfound << " of " << v_signatures.size()
               << " bubbles" << endl;    }
     return nfound;    }

void VariantsMarkedDot(const String& head, const SupportedHyperBasevector& shb, 
       const vec<int>& varient_edge_list, const long_logging& logc )
{
     double clock = WallClockTime( );
     const int Varient_Pen_Width = 4;
     vec<Bool> invisible( shb.EdgeObjectCount( ), false );
     vec<int> varient_edges( varient_edge_list );
     Sort(varient_edges);
     // Below are copied from SupportedHyperBasevector::Dumpdot
     vec<Bool> hide;
     FlagEdgesForHiding( shb, shb.Inv( ), hide, logc );
     const Bool DOT_LABEL_CONTIGS = True;
     const Bool DOT_LABEL_VERTICES = False;
     vec<double> lengths( shb.EdgeObjectCount( ) );
     for ( int i = 0; i < shb.EdgeObjectCount( ); i++ )
         lengths[i] = shb.EdgeLengthKmers(i);
     vec<String> edge_id_names( shb.EdgeObjectCount( ) );
     for ( int i = 0; i < shb.EdgeObjectCount( ); i++ ) {    
         if ( !shb.InvDef(i) ) edge_id_names[i] = ToString(i);
         else edge_id_names[i] = ToString(i) + "=" + ToString( shb.Inv(i) ) + "'";    
     }
     vec<int> pen_widths(shb.EdgeObjectCount(),0);
     for ( int i = 0; i < shb.EdgeObjectCount( ); i++ ) {    
         if (BinMember(varient_edges,i))
             pen_widths[i] = Varient_Pen_Width;
     }
     Ofstream( dout, head + ".dot" );
     shb.PrettyDOT( dout, lengths, HyperBasevector::edge_label_info(
          HyperBasevector::edge_label_info::DIRECT, &edge_id_names ),
          DOT_LABEL_CONTIGS, DOT_LABEL_VERTICES, NULL, NULL, NULL, &hide,
          &invisible, NULL, &pen_widths );    
     REPORT_TIME( clock, "used in VariantsMarkedDot" );    
}

void ReftraceVariants(ostream& out, ostream& callsOut, const vecbasevector& G, 
        const vecbasevector& Gplus, const vec<int>& Gplus_ext, 
        const HyperBasevector& hbp, const vec<pair<int,Bool>>& hbp_to_hb,
        const vec<vec<vec<pair<int,int>>>>& ref_pos_all_s, const vec<vec<vec<int>>>& best_path_s,
        const String& sVarFile, const String& sRefFile, 
        const RefTraceControl& rtc,
        const ReadOriginTracker* p_read_tracker,
        const long_logging* logc)
{
    int verbosity = logc->verb["REFTRACE_VARIANTS"];
    const vec<size_t>& ref_index = rtc.getRefIndex();
    const vec<String>& ref_tags = rtc.getRefTags();
    const vec<size_t>& ref_start = rtc.getRefStart();
    const vec<fastavector>& ref_seqs = rtc.getRefSeqs();

    bool UseInternalEdgeID = false;
    bool UseInternalCoord  = false;
    bool bExternalSeqs = ref_seqs.size() > 0;
    if(bExternalSeqs){
        ForceAssert(ref_seqs.size()==G.size());
        for(size_t ii=0 ; ii < ref_seqs.size() ; ++ii){
            ForceAssert(ref_seqs[ii].size()==G[ii].size());
        }
    }
    bool bExternalTags = ref_tags.size() > 0;
    if(bExternalTags){
        ForceAssert(ref_tags.size()==G.size());
    }

    // Find variant by unrolling assembly on the reference, guided by best path.
    vec<VariantCallGroup> groups;
    vec<vec<align>> all_aligns(G.size());  // alignments of all edges
    vec<EdgesOnRef> unrolled_graph;
    for (size_t g = 0; g < G.size(); g++) {
        for( size_t component=0;component<best_path_s[g].size();++component){
            if(verbosity) std::cout<< "working on (g,component) (" <<  g << "," << component << ")"<< std::endl;
            unrolled_graph.push(EdgesOnRef(hbp, hbp_to_hb, G[g], Gplus[g], 
                        Gplus_ext[g], p_read_tracker));
            EdgesOnRef& graph = unrolled_graph.back();//[g];
            graph.InitFromBestPath(best_path_s[g][component], ref_pos_all_s[g][component]);
            const int lo = (component==0)?std::numeric_limits<int>::min()
                                         :ref_pos_all_s[g][component-1].front().first;
            const int hi = (component+1<best_path_s[g].size())?ref_pos_all_s[g][component+1].back().second
                                                              :std::numeric_limits<int>::max();
            graph.UnrollAll(verbosity,lo,hi);
            graph.MakeBubbleGraph(verbosity);
            if (p_read_tracker != NULL)
                graph.PathProb(p_read_tracker->Reads(), p_read_tracker->Quals(), verbosity);
            graph.CallVariantsGroupedWithProb(g, &groups, &all_aligns[g], verbosity);
        }
    }
    if (verbosity >= 1) {
        for (size_t g = 0, index=0; g < best_path_s.size(); g++) {
            for( size_t component=0;component<best_path_s[g].size();++component,++index){
                unrolled_graph[index].DumpUnrolled( "unrolled" + ToString(g) 
                                                  + "_" + ToString(component) + ".dot", 
                        (UseInternalEdgeID ? NULL : &hbp_to_hb) );
            }
        }
    }

    // remove variants on edges with multiple placement on the genome
    if (logc->IGNORE_VARIANT_ON_REPETITIVE_EDGE) {
        ForceAssertEq(G.size(), ref_tags.size());
        ForceAssertEq(G.size(), ref_start.size());
        //String ref_m100_file = "/wga/scr4/bigrefs/human19/genome.m100";
        String ref_m100_file = sRefFile + ".fasta.m100";
        RemoveRepetitiveEdges(groups, ref_index, ref_start, ref_m100_file, 
                hbp, hbp_to_hb, verbosity);
    }

    // Find variant friends.
    map<Variant, vec<pair<int,int>>> var_friending;
    FindVariantFriends(groups, all_aligns, hbp, hbp_to_hb,  &var_friending);

    // Calculate the probability of variants
    UniqueSort(groups);
    map<Variant, vec<pair<double,double>>> probs;
    if (p_read_tracker != NULL) {
        FindVariantProb(p_read_tracker, groups, Gplus, Gplus_ext, probs, logc); 
    }

    TextTable tb2;
    VCFWriter vcf;
    const vec<String> vsDefaultSample={"ALL_EDGES"};
    bool bVCFOut=true;//p_read_tracker!=NULL;
    const auto& sample_list = (p_read_tracker!=NULL)?p_read_tracker->getSampleList():vsDefaultSample;
    for( const auto& entry: sample_list){
        vcf.setSample(entry,"unphased representation of DISCOVAR variant ouput for dataset "+entry);
    }
    vcf.setFormat("GP","G","Float","phred-scaled genotype posterior probabilities");
    vcf.setFormat("REFP","1","Float","phred-scaled probability that reference exists");
    vcf.setFormat("ALTP","A","Float","phred-scaled probability that alternatives exist");
    //vcf.setSample("EMPTY","This sample should be empty. For demonstration purposes");

    ForceAssertEq(ref_tags.size(), ref_seqs.size());
    ForceAssertEq(ref_tags.size(), ref_start.size());
    for ( size_t i = 0; i < ref_tags.size(); ++i ) {
	vcf.AddMeta("DiscovarRegion", {
		{ "CHR", QuoteString(ref_tags[i]) },
		{ "START", ToString( ref_start[i] ) },
		{ "END", ToString( ref_start[i]+ref_seqs[i].size() ) } } );
    }
    vcf.AddInfo( VCFWriter::info_t("BL", "1", "Integer", "Discovar variant block" ) );

    size_t max_vindex=0;
    {
        set<Variant> tmp;
        for( const auto& group:groups){
          for( const auto& vcall_vec: group.GetVariantCalls()){
            for (const VariantCall& x: vcall_vec) {
              tmp.insert(x.variant); } } }
        max_vindex=tmp.size();
        if(max_vindex<10)max_vindex=10;
    }

    //table headers
    tb2 << "#BLOCK" << Tab<TextTable> << "ID" << Tab<TextTable> << "PATHS" << EndRow<TextTable>;
    tb2 << "#VAR" << Tab<TextTable> << "ID" << Tab<TextTable> << "LOCATION" << Tab<TextTable> << "REF" << Tab<TextTable> << "ALT" << Tab<TextTable> << "PHASE"<< Tab<TextTable> <<"INFO"<<EndRow<TextTable>;

    int var_index = 0, block_index = 0;			// indexes are updated in PrintTable
    for (size_t i = 0; i < groups.size(); i++) {
        vec<size_t> ref_shift_empty;
        groups[i].PrintTable(tb2, callsOut, var_index, block_index, hbp_to_hb, ref_tags,
                (UseInternalCoord ? ref_shift_empty: ref_start), 
                var_friending, max_vindex, probs);
        groups[i].FillVCF(vcf,  hbp_to_hb, ref_tags,
                (UseInternalCoord ? ref_shift_empty: ref_start), var_friending, probs,sample_list, block_index);
    }

    std::ofstream ofs,ofs_vcf;
    if( sVarFile!=""){
        ofs.open(sVarFile.c_str(),std::ios::out);
        if(bVCFOut) ofs_vcf.open((sVarFile+".vcf").c_str(),std::ios::out);
    }
    {
        std::ostream& effective_out = (ofs.is_open() ) ? ofs : out;
        tb2.Print(effective_out, 1, "llllllll", true);
    }

    if(ofs.is_open()){ofs.close(); }
    if(ofs_vcf.is_open()){
        vcf.sort();
        ofs_vcf << vcf;
        ofs_vcf.close();
        if(filter_vcf(sVarFile+".vcf",sVarFile+".filtered.vcf")){
            std::cout << "There has been an error filtering VCF. " << sVarFile+".filtered.vcf has been removed." << std::endl;
        }
    }

}

void VariantCallGroup::PrintTable(TextTable& tb, ostream& callsOut, int& index, int& block_id,
        const vec<pair<int,Bool>>& hbp_to_hb, 
        const vec<String> ref_tags, const vec<size_t> ref_start,
        const map<Variant, vec<pair<int,int>>>& var_friending,
        size_t max_index, const map<Variant, vec<pair<double,double>>> probs) const
{
    const bool PrintPathProb =  ( ref_weight_.size() && weights_.size() );  // print, if we have weights
    const size_t MaxSeqLen = 10;
    map<Variant, set<int>> variant_callers; 
    map<Variant, vec<double>> variant_prob;
    for (size_t i = 0; i < vcalls_.size(); i++) {
        const vec<VariantCall>& vcall_vec = vcalls_[i];
        for (const VariantCall& x: vcall_vec) 
           variant_callers[x.variant].insert(i); 
        // Add probabilities
        set<Variant> vs;
        for (const VariantCall& x: vcall_vec) 
           vs.insert(x.variant); 
        for (auto& x: vs) {
            // add new record if needed
            int nsample = weights_[i].size();
            variant_prob.insert(make_pair(x,vec<double>(nsample,0.0)));
            for (int sample=0; sample<nsample; sample++){
                variant_prob[x][sample] += weights_[i][sample];
            }
        }
    }
    if (variant_callers.empty()) return;

    const int block_id_loc=++block_id;
    int shift=0;
    for(size_t factor=10; max_index/factor>0 ; factor*=10,++shift){};
    for(size_t factor=10; block_id_loc/factor>0 ; factor*=10,--shift){};
    if(shift<0) shift=0;

    tb << "BLOCK  " << block_id_loc << String(1+shift,' ');
    callsOut << "BLOCK\t" << block_id;
    char prefix = '\t';
    for (size_t i = 0; i < branches_.size(); i++) {
        if (i != 0) tb << "|";
        const vec<int>& path = GetBranchPath(i);
        for (size_t j = 0; j < path.size(); j++) {
            if (j != 0) tb << ",";
            char sign = (hbp_to_hb[path[j]].second ? '+' : '-');
            int edge0 = hbp_to_hb[path[j]].first;
            tb << edge0 << sign;
            callsOut << prefix << edge0 << sign;
            prefix = '|';
        }
    }
    tb.SetRawLine(); // do not auto adjust this line
    tb << EndRow<TextTable>;
    callsOut << '\n';

    if (PrintPathProb) {
        tb << "PATH_QSUM ";
        for (size_t i = 0; i < weights_.size(); i++) {
            if (i != 0) tb << "|";
            for (size_t sample = 0; sample < weights_[i].size(); sample++) {
                tb << (sample>0?String(","):String("")) + ToString(weights_[i][sample]);
            }
        }
        tb << "|ref=";
        for (size_t sample = 0; sample < ref_weight_.size(); sample++) 
            tb << (sample>0?String(","):String("")) + ToString(ref_weight_[sample]);
        tb.SetRawLine(); // do not auto adjust this line
        tb << EndRow<TextTable>;
    }

    for (auto it = variant_callers.begin(); it != variant_callers.end(); ++it) {
        const Variant& var = it->first;
        const set<int>& pathids = it->second;
        tb << "VAR  " << Tab<TextTable>;
        tb << ++index << Tab<TextTable>;
        String tag = (ref_tags.empty() ? ToString(var.gid+1) : ref_tags[var.gid]);

        int64_t shift = (ref_start.empty() ? 0 : ref_start[var.gid]);
        int64_t genome_pos_one_based  = shift + var.pos + 1;
        tb << tag << ":" << genome_pos_one_based << Tab<TextTable>;

        String ref = var.ref, alt = var.alt;
        if (ref.size() > MaxSeqLen) {
            ref.resize(MaxSeqLen);
            ref.back() = '+';
        }
        if (alt.size() > MaxSeqLen) {
            alt.resize(MaxSeqLen);
            alt.back() = '+';
        }
        tb << ref << Tab<TextTable>;
        tb << alt << Tab<TextTable>;
        callsOut << "VAR\t" << index << '\t' << tag << ':'
                 << genome_pos_one_based << '\t'
                 << var.ref << '\t' << var.alt << '\t';

        for (size_t k = 0; k < branches_.size(); k++) {
            char mark = (pathids.find(k) != pathids.end() ? '+' : '-');
            if (k != 0) tb << '|';
            tb << mark;
            callsOut << mark;
        }
        tb << Tab<TextTable>;

        // Any edge not containing this variant and other overlapping variants 
        // will be considered supporting reference. 
        vec<bool> vbRefPresent(branches_.size(),false);
        vec<bool> vbAltPresent(branches_.size(),false);
        int var_start = var.pos;
        int var_end = var.pos + var.ref.size();
        for (size_t k = 0; k < branches_.size(); k++) {
            if(pathids.find(k) == pathids.end()){
                bool overlap = false;
                const vec<VariantCall>& vcall_vec = vcalls_[k];
                for (const VariantCall& x: vcall_vec) {
                    int var_start2 = x.variant.pos;
                    int var_end2 = x.variant.pos + x.variant.ref.size();
                    if (max(var_start, var_start2) < min(var_end, var_end2)) {
                        overlap = true;
                        break;
                    }
                }
                if (!overlap) vbRefPresent[k] = true;
            } else {
                vbAltPresent[k] = true;
            }
        }
        tb << Tab<TextTable>;
        if (PrintPathProb && !probs.empty()){
            tb << setprecision(2);
            tb << "PROB= ";
            auto it = probs.find(var);
            if (it != probs.end()) {
                for (size_t sample = 0; sample < it->second.size(); sample++){
                    if (sample != 0) tb << " ";
                    tb << 1 - pow(10, - it->second[sample].first/10) << "/"
                       << 1 - pow(10, - it->second[sample].second/10);
                }
            }
        }
        if (var.ref.size() > MaxSeqLen || var.alt.size() > MaxSeqLen) {
            tb << " REF=" << var.ref;
            tb << ";ALT=" << var.alt;
        }
        auto it2 = var_friending.find(var);
        if (it2 != var_friending.end()) {
            const vec<pair<int,int>>& friend_locs = it2->second;

            vec<size_t> real_friends;
            for (size_t k = 0; k < friend_locs.size(); k++) {
                int gid2 = friend_locs[k].first;
                if( gid2!=var.gid || friend_locs[k].second != var.pos ) real_friends.push_back(k);
            }
            if( real_friends.size()>0){
                tb << " FRIENDS=";
                for (auto k: real_friends){
                    if (k != 0) tb << ",";
                    int64_t shift2 = (ref_start.empty() ? 0 : ref_start[friend_locs[k].first]);
                    int gid2 = friend_locs[k].first;
                    String tag = (ref_tags.empty() ? ToString(gid2+1) : ref_tags[gid2]);
                    tb << tag << ":" <<friend_locs[k].second +1 + shift2;
                    callsOut << '\t' << friend_locs[k].first + 1 << ":"
                             << friend_locs[k].second +1 + shift2;
                }
            }
        }
        tb << EndRow<TextTable>;
        callsOut << '\n';
    }
    // Adding empty row at block end
    tb.SetRawLine(); tb << EndRow<TextTable>;
}

void VariantCallGroup::FillVCF(VCFWriter& vcf
                              ,const vec<pair<int,Bool>>& hbp_to_hb
                              ,const vec<String> ref_tags
                              ,const vec<size_t> ref_start
                              ,const map<Variant, vec<pair<int,int>>>& var_friending
                              ,const map<Variant, vec<pair<double,double>>> probs
                              ,const vec<String>& sample_list
                              ,const size_t block_no) const
{
    const bool PrintPathProb = true;
    const size_t MaxSeqLen = 10;
    size_t nSamples=0;
    map<Variant, set<int>> variant_callers;
    map<Variant, vec<double>> variant_prob;
    for (size_t i = 0; i < vcalls_.size(); i++) {
        const vec<VariantCall>& vcall_vec = vcalls_[i];
        for (const VariantCall& x: vcall_vec)
           variant_callers[x.variant].insert(i);
        // Add probabilities
        set<Variant> vs;
        for (const VariantCall& x: vcall_vec)
           vs.insert(x.variant);
        for (auto& x: vs) {
            // add new record if needed
            int nsample = weights_[i].size();
            nSamples = std::max(nSamples,size_t(nsample));
            variant_prob.insert(make_pair(x,vec<double>(nsample,0.0)));
            for (int sample=0; sample<nsample; sample++){
                variant_prob[x][sample] += weights_[i][sample];
            }
        }
    }
    if (variant_callers.empty()) return;

    enum {MERGE_GID,MERGE_FRONT,MERGE_BACK,MERGE_REF,MERGE_EDGES,MERGE_VARS};
    typedef std::tuple< int64_t, int64_t, int64_t, String,vec<String>,Variant> entry_t;
    vec<entry_t> vcf_entries;

    for (auto it = variant_callers.begin(); it != variant_callers.end(); ++it) {
        const Variant& var = it->first;
        const set<int>& pathids = it->second;
        String tag = (ref_tags.empty() ? ToString(var.gid+1) : ref_tags[var.gid]);

        int64_t shift = (ref_start.empty() ? 0 : ref_start[var.gid]);
        int64_t genome_pos_one_based  = shift + var.pos + 1;

        const int64_t back_one_based = genome_pos_one_based + var.ref.size()-1;
        vec<String> edge_sequences(branches_.size());
        vec<bool> vbRefPresent(branches_.size(),false);
        for (size_t k = 0; k < branches_.size(); k++) {
            if(pathids.find(k) != pathids.end()){
                edge_sequences[k]=var.alt;
            }
            else{
                vbRefPresent[k]=true;
                edge_sequences[k]=var.ref;
            }
        }
        vcf_entries.push_back( std::make_tuple(var.gid,genome_pos_one_based,back_one_based,var.ref,edge_sequences,var));
    }
    Sort(vcf_entries); //sort them with gid, then front, then back
    // go through the entries, merge entries with same gid, and covering overlaping region (determined from front/back) in the reference
    for( size_t ee=0;ee< vcf_entries.size();){
        const auto& start_entry = vcf_entries[ee];
        auto gid = std::get<MERGE_GID>(start_entry);
        auto range_front = std::get<MERGE_FRONT>(start_entry);
        auto range_back = std::get<MERGE_BACK>(start_entry);
        String ref = std::get<MERGE_REF>(start_entry);
        size_t ff = ee+1;
        for( ;   ff < vcf_entries.size()
              && std::get<MERGE_GID>(vcf_entries[ff]) == gid
              && std::get<MERGE_FRONT>(vcf_entries[ff]) <= range_back
             ; ++ff){
            auto loc_back= std::get<MERGE_BACK>(vcf_entries[ff]);
            if( loc_back > range_back){
                size_t extra =loc_back - range_back;
                const auto& loc_ref=std::get<MERGE_REF>(vcf_entries[ff]);
                ref.append(loc_ref,loc_ref.size()-extra,extra);
                range_back=loc_back;
            }
        }
        vec<String> edge_alts(branches_.size(),ref);
        vec< vec<Variant>> edge_vars(branches_.size());
        for(size_t edge=0;edge<branches_.size();++edge){
            for(size_t redit=0;redit+ee<ff;++redit){
                size_t edit = ff-1-redit;
                const auto& loc_entry=vcf_entries[edit];
                const String& loc_alt = std::get<MERGE_EDGES>(loc_entry)[edge];
                const String& loc_ref = std::get<MERGE_REF>(loc_entry);
                if( loc_alt != loc_ref){
                    size_t pos = std::get<MERGE_FRONT>(loc_entry) - range_front;
                    size_t len = loc_ref.size();
                    ForceAssert( ref.substr(pos,len) == loc_ref);
                    edge_alts[edge].replace(pos,len,loc_alt);
                    edge_vars[edge].push_back(std::get<MERGE_VARS>(loc_entry));
                }
            }
        }
        vec<String> sequences={ref};
        vec< vec<double> > sequence_weights;
        double dAltWeight=0.0;


        if(probs.empty()){
            sequence_weights.push_back(ref_weight_);
            for(size_t edge=0;edge<branches_.size();++edge){
                size_t ss=0;
                for(; ss < sequences.size() && edge_alts[edge]!=sequences[ss]; ++ss){ }
                if( ss!=0){
                    for( size_t si=0;si<nSamples;++si){
                        dAltWeight+=weights_[edge][si];
                    }
                }
                if( ss >= sequences.size()){
                    sequences.push_back(edge_alts[edge]);
                    sequence_weights.push_back(weights_[edge]);
                }
                else{
                    for( size_t si=0;si<nSamples;++si){
                        sequence_weights[ss][si]+=weights_[edge][si];
                    }
                }
            }
        }
        else{
            sequence_weights.push_back( vec<double>(nSamples,std::numeric_limits<double>::max()));
            for(size_t edge=0;edge<branches_.size();++edge){
                vec<double> dvAlt=vec<double>(nSamples,std::numeric_limits<double>::max());
                for( const auto& var:edge_vars[edge]){
                    const auto& it = probs.find(var);
                    if (it != probs.end()) {
                        for(size_t si=0;si<nSamples;++si){
                            sequence_weights[0][si] = std::min(sequence_weights[0][si],(*it).second[si].first);
                            dvAlt[si] = std::min(dvAlt[si],(*it).second[si].second);
                        }
                    }
                }
                size_t ss=0;
                for(; ss < sequences.size() && edge_alts[edge]!=sequences[ss]; ++ss){ }
                if( ss >= sequences.size()){
                    sequences.push_back(edge_alts[edge]);
                    sequence_weights.push_back(dvAlt);
                }
                else{
                    for( size_t si=0;si<nSamples;++si){
                        sequence_weights[ss][si]=std::min(sequence_weights[ss][si],dvAlt[si]);
                    }
                }
            }
            for( auto& a: sequence_weights){
                replace(a.begin(),a.end(),std::numeric_limits<double>::max(),0.0);
            }
            for(size_t ii=1;ii<sequence_weights.size();++ii){
                dAltWeight+=std::accumulate( sequence_weights[ii].begin(),sequence_weights[ii].end(),0.0);
            }
        }
        std::map<String,vec<std::pair<String,String>>>  gtfs;
        for(size_t sample=0;sample<nSamples;++sample){
            auto& sample_name = (sample<sample_list.size())?sample_list[sample]:ToString(sample);

            vec<double> p_exists(sequence_weights.size());
            for(size_t seq=0;seq<p_exists.size();++seq){
                p_exists[seq] = 1.0 - pow(10.0, - sequence_weights[seq][sample]/10.0);
            }

            String GP;
            vec<double> vd_q_genotype;

            for(size_t gj=0 ; gj<p_exists.size();++gj){
                for(size_t gi=0 ; gi<=gj;++gi){
                    double p_genotype=1.0;
                    for(size_t gk=0;gk<p_exists.size();++gk){
                        if( gk==gi || gk==gj){ p_genotype *= p_exists[gk]; }
                        else{ p_genotype *= 1.0-p_exists[gk]; }
                    }
                    double q_genotype = -10.0*log10(1.0-p_genotype);
                    if(q_genotype<std::numeric_limits<float>::epsilon()) q_genotype=0;
                    if(q_genotype==std::numeric_limits<double>::infinity()){
                        q_genotype=sequence_weights[gj][sample];
                        if(gi!=gj){
                            auto mima = std::minmax(sequence_weights[gi][sample],sequence_weights[gj][sample]);
                            q_genotype = mima.first -10.0*log10(1+ pow(10.,-0.1*(mima.second-mima.first)) - pow(10.,-0.1*(mima.second)));
                            if(q_genotype==std::numeric_limits<double>::infinity()){
                                q_genotype=mima.first;
                            }
                        }
                    }
                    vd_q_genotype.push_back(q_genotype);
                    if(GP.size()>0) GP+=",";
                    GP+= ToString( q_genotype);
                }
            }

            String GT;
            double max_q=0.0;
            for(size_t gj=0,count=0 ; gj<p_exists.size();++gj){
                for(size_t gi=0 ; gi<=gj;++gi,++count){
                    if( vd_q_genotype[count] >= max_q){
                        max_q = vd_q_genotype[count];
                        GT = ToString(gi)+"/"+ToString(gj);
                    }

                }
            }

            gtfs[sample_name].push_back( std::make_pair("GT",GT));
            gtfs[sample_name].push_back( std::make_pair("GP",GP));
            gtfs[sample_name].push_back( std::make_pair("REFP",ToString(sequence_weights[0][sample])));
            String altp;
            for(size_t aa=1;aa<p_exists.size();++aa){
                if(aa>1) altp+=',';
                altp+=ToString(sequence_weights[aa][sample]);
            }
            gtfs[sample_name].push_back( std::make_pair("ALTP",altp));
        }
        String tag = (ref_tags.empty() ? ToString(gid+1) : ref_tags[gid]);
        vec<std::pair<String,String> > infos = {{ "BL", ToString( block_no ) }};
        vcf.AddEntry(tag,range_front,"",ref,vec<String>(sequences.begin()+1,sequences.end()),ToString(dAltWeight),vec<String>(),infos,gtfs);
        ee=ff;
    }
}

bool operator==(const VariantCallGroup&left,const VariantCallGroup&right){
    if( left.GetGID() != right.GetGID()) return false;
    std::set< int > l_loc,r_loc;
    for( const auto& branch: left.GetVariantCalls()){
        for(const auto& entry: branch){
            l_loc.insert( entry.variant.pos );
        }
    }
    for( const auto& branch: right.GetVariantCalls()){
        for(const auto& entry: branch){
            r_loc.insert( entry.variant.pos );
        }
    }
    return l_loc==r_loc;
}
bool operator<(const VariantCallGroup&left,const VariantCallGroup&right){
    return (   left.GetGID() < right.GetGID()
            || (left.GetGID()==right.GetGID() && left.GetStartPos() < right.GetStartPos()));
}
