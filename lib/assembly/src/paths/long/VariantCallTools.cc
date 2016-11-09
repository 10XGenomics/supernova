///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "paths/long/VariantCallTools.h"

#include "CoreTools.h"
#include "Map.h"
#include "math/Functions.h"
#include "paths/RemodelGapTools.h"
#include "paths/long/Logging.h"
#include "paths/long/SupportedHyperBasevector.h"
#include "pairwise_aligners/SmithWatAffine.h"
#include "PrintAlignment.h"
#include "util/TextTable.h"
#include "graph/Digraph.h"
#include "graph/DigraphTemplate.h"
#include "paths/long/ReadOriginTracker.h"
#include "paths/long/EvalByReads.h"
#include "paths/long/VariantReadSupport.h"
#include "kmers/KMer.h"

namespace {

bool IsBubbleEdge(int eid, const HyperBasevector& hb, const vec<int>& to_left,
        const vec<int>& to_right)
{
    int node_left = to_left[eid], node_right = to_right[eid];
    bool ends_match = true;
    for (int j = 0; j < hb.FromSize(node_left); j++) {
        int edge_alt = hb.EdgeObjectIndexByIndexFrom(node_left, j);
        if (to_right[edge_alt] != node_right) {
            ends_match = false;
            break;
        }
    }
    return ends_match && hb.FromSize(node_left) > 1
        && hb.ToSize(node_left) == 1 && hb.FromSize(node_right) == 1;
}

// Find half of the edges that are most similar to the selected target
void FindSimilar(const basevector& selected, const vec<basevector>& edges,
        const vec<int>& edge_ids, vec<int>& similar_edges) 
{
    if (edges.size() < 2) return;
    const int K = 12;
    const double MinSimilarity = 0.8;
    typedef KMer<K> KmerT;
    vec<KmerT> kmers;
    for (size_t i = 0; i < selected.size() - K + 1; i++) 
        kmers.push(selected.begin() + i);
    UniqueSort(kmers);
    vec<pair<double,int>> sim_index;
    for (size_t i = 0; i < edges.size(); i++) {
        //cout << "edge " << i << ": ";
        vec<KmerT> kmers2;
        for (size_t j = 0; j < edges[i].size() - K + 1; j++) {
            kmers2.push(edges[i].begin() + j);
        }
        UniqueSort(kmers2);
        int common = 0, unique1 = 0, unique2 = 0;
        auto it1 = kmers.begin(), it2 = kmers2.begin();
        while (it1 != kmers.end() || it2 != kmers.end()) {
            if (it1 == kmers.end()) {
                unique2 += distance(it2, kmers2.end());
                it2 == kmers2.end();
                break;
            }
            if (it2 == kmers2.end()) {
                unique1 += distance(it1, kmers.end());
                it1 == kmers.end();
                break;
            }
            if ( (*it1) == (*it2) ) {
                common++;
                it1++, it2++;
            } else if (*it1 < *it2) {
                unique1++;
                it1++;
            } else {
                unique2++;
                it2++;
            }
        }
        double similarity = (double) common/ (double)(common + unique1 + unique2);
        //cout << "common= " << common << " unique1= " << unique1 
        //    << " unique2= " << unique2 << " simi= " << similarity << endl;
        sim_index.push(similarity, i);
    }
    Sort(sim_index);
    for (auto x: sim_index) {
        if (x.first > MinSimilarity) break;
        similar_edges.push_back(edge_ids[x.second]);
    }
}


// RestrictedAlign: given a proposed SmithWatAffineBanded alignment with given
// offset and bandwidth, do the alignment, but first try to choose a better offset
// and bandwidth.  The first step is to find all kmer matches between the two
// sequences (K = 40).  Each such match defines a possible offset, and we take the
// middle 90% of these offsets to define initial bounds on the offset.  Then we
// traverse the remaining offsets.  We ignore offsets for which there is another
// offset at the same position on b1 that are within the 90% bounds.  All other
// offsets are put into a pile called 'extras'.  At first we align, ignoring the
// extras.  Then we ask if use of any of the extras might possibly have improved
// the alignment (using a test that is not really airtight).  If so we extend the 
// offset bounds and try again.
//
// Note that the kmer lookup is parallelized.

void RestrictedAlign( const basevector& b1, const basevector& b2, 
     const int offset, const int bandwidth, align& a )
{
     // Control.

     const int L = 40;
     const int mul = 20;
     const int bw_add = 10;
     const int bw_mul = 10;
     const Bool compare_to_old = False;
     const Bool exit_bad = True; // not in effect
     const Bool verbose = False;

#if 0
     int smallest = min( b1.isize(), b2.isize() );
     if ( smallest > 5000 ) {
         int K = 501; //max( 60, min( smallest/100, 501 ) );
         alignment al;
         SmithWatAffineSuper( b1, b2, al, K, false, false );
         a=align(al);
         return;
     }
#endif

     // Find offsets.

     int low = Max( 0, offset - bandwidth );
     int high = Min( b2.isize( ), offset + bandwidth + b1.isize( ) );
     vecbasevector x(2);
     x[0] = b1, x[1] = basevector( b2, low, high - low );
     vec< triple<kmer<L>,int,int> > kmers_plus;
     MakeKmerLookup0( x, kmers_plus );
     vec<int> offsets, pos;
     for ( int i = 0; i < kmers_plus.isize( ); i++ )
     {    int j;
          for ( j = i + 1; j < kmers_plus.isize( ); j++ )
               if ( kmers_plus[j].first != kmers_plus[i].first ) break;
          int m;
          for ( m = i; m < j; m++ )
               if ( kmers_plus[m].second == 1 ) break;
          for ( int k = i; k < m; k++ )
          for ( int l = m; l < j; l++ )
          {    offsets.push( kmers_plus[l].third - kmers_plus[k].third + low );    
               pos.push_back( kmers_plus[k].third );    }
          i = j - 1;    }
     int o1 = -1, o2 = -1, offset2, bandwidth2;
     vec<int> extras;

     // If there are no offsets, use the original offset and bandwidth.

     if ( offsets.empty( ) )
     {    offset2 = offset;
          bandwidth2 = bandwidth;    }

     // Main case.

     else
     {    Sort(offsets);

          // Get 90% bounds.

          o1 = offsets[ offsets.size( ) / mul ];
          o2 = offsets[ ( (mul-1) * offsets.size( ) ) / mul ];
          ostringstream out;
          out << "90% = [" << o1 << "," << o2 << "]" << endl;

          // Look at offsets associated with each position, and find 'extra'
          // offsets.

          SortSync( pos, offsets );
          for ( int i = 0; i < pos.isize( ); i++ )
          {    int j = pos.NextDiff(i);
               //Bool inside = False;
               //for ( int k = i; k < j; k++ )
               //     if ( o1 <= offsets[k] && offsets[k] <= o2 ) inside = True;
               // if ( !inside )
               {    for ( int k = i; k < j; k++ )
                    {    if ( o1 <= offsets[k] && offsets[k] <= o2 ) continue;
                         //out << pos[k] << " goes to offset " << offsets[k] << endl;
                         extras.push_back( offsets[k] );    }    }
               i = j - 1;    }
          
          // Define new offset.

          offset2 = ( o2 + o1 ) / 2;
          bandwidth2 = Max( o2 - offset2, offset2 - o1 ) + bw_add;
          if (verbose)
          {    cout << "\n";
               PRINT2( b1.size( ), b2.size( ) );
               PRINT2( o1, o2 );
               PRINT4( offset2, bandwidth2, offset, bandwidth );
               cout << out.str( ) << "\n";    }    }

     // Align.

     align a2;
     if (verbose) cout << Date( ) << ": alignment starting" << endl;
     int nerrors;
     int err = SmithWatAffineBanded( b1, b2, -offset2, bandwidth2, a2, nerrors );
     if (verbose)
     {    cout << Date( ) << ": alignment complete" << endl;
          PRINT(err);    }
     /*
     if ( exit_bad && err > 10000 )
     {    PRINT(err);
          b1.Print( cout, "b1" );
          b2.Print( cout, "b2" );
          Scram(0);     }
     */
     
     // If we might possibly have done better using some of the extra
     // offsets, add them in and try again.  The test used here is not quite right.

     Sort(extras);
     Bool retry = False;
     int o1_new(o1), o2_new(o2);
     const int mismatch_penalty = 3;
     const int gap_open_penalty = 12;
     const int gap_extend_penalty = 1;
     for ( int i = 0; i < extras.isize( ); i++ )
     {    if ( extras[i] < o1_new || extras[i] > o2_new )
          {    int gap = Max( o2, extras[i] ) - Min( o1, extras[i] );
               int min_penalty = gap_open_penalty + (gap-1) * gap_extend_penalty;
               if ( min_penalty < err )
               {    retry = True;
                    o1_new = Min( o1_new, extras[i] );
                    o2_new = Max( o2_new, extras[i] );    }    }    }
     if (retry)
     {    int bandwidthx = Max( o2_new - offset2, offset2 - o1_new );
          int bandwidth3 = bandwidth2;
          while ( bandwidth3 < bandwidthx )
          {    bandwidth3 = Min( bandwidthx, bw_mul * bandwidth3 );
               if (verbose)
               {    cout << "retrying" << endl;
                    PRINT3( err, offset2, bandwidth3 );    }
               err = SmithWatAffineBanded(b1, b2, -offset2, bandwidth3, a2, nerrors);
               int gap = 2 * bandwidth3;
               int min_penalty = gap_open_penalty + (gap-1) * gap_extend_penalty;
               if ( min_penalty > err ) break;    }    }

     // Test to see if results are different than old method.

     if (compare_to_old)
     {    int nerror;
          int err_old = SmithWatAffineBanded(b1, b2, -offset, bandwidth, a, nerror);
          if ( !( a2 == a ) ) 
          {    cout << "not equal!\n";
               PRINT2( err, err_old );
               cout << "offsets in old alignment:\n";
               int p1 = a.pos1( ), p2 = a.pos2( );
               cout << p2 - p1 << endl;
               for ( int j = 0; j < a.Nblocks( ); j++ ) 
               {    if ( a.Gaps(j) > 0 )  
                    {    p2 += a.Gaps(j);
                         cout << p2 - p1 << endl;    }
                    if ( a.Gaps(j) < 0 ) 
                    {    p1 -= a.Gaps(j);
                         cout << p2 - p1 << endl;    }
                    p1 += a.Lengths(j), p2 += a.Lengths(j);    }
               // b1.Print( cout, "b1" );
               // b2.Print( cout, "b2" );
               cout << "old alignment:\n";
               PRINT3( a.pos2( ), a.Pos2( ), b2.size( ) );
               PrintVisualAlignment( True, cout, b1, b2, a );
               cout << "new alignment:\n";
               PrintVisualAlignment( True, cout, b1, b2, a2 );
               Scram(1);    }    }

     // Done.

     if (verbose) cout << "done" << endl;
     a = a2;
}

int CorrelatePositionsRangeChecked( const align& a, const int x1 )
{    
     int pos1 = a.pos1( ), pos2 = a.pos2( );
     if ( x1 < pos1 || x1 > a.Pos1()) return -1; // off the end, shouldn't happen
     if ( x1 == pos1 ) return pos2;
     int nblocks = a.Nblocks( );
     const avector<int> &gaps = a.Gaps( ), &lengths = a.Lengths( );
     for ( int j = 0; j < nblocks; j++ )
     {    if ( gaps(j) > 0 ) pos2 += gaps(j);
          if ( gaps(j) < 0 )
          {    for ( int x = 0; x < -gaps(j); x++ )
                    if ( x1 == pos1++ ) return pos2;    } // really in gap
          if ( x1 < pos1 + lengths(j) ) return pos2 + x1 - pos1;
          else
          {    pos1 += lengths(j), pos2 += lengths(j);    }    }
     return -1; // off the end, shouldn't happen
} 

template<typename GraphT>
void DumpGraph(const String filename, const GraphT& shb) {
    vec<String> edge_names2(shb.EdgeObjectCount(),"");
    vec<double> lengths2( shb.EdgeObjectCount( ) );
    for (size_t i = 0; i < edge_names2.size(); ++i) {
        edge_names2[i] = ToString(i);
        lengths2[i] = shb.EdgeLengthKmers(i);
    }
    ofstream dout(filename);
    shb.PrettyDOT( dout, lengths2, HyperBasevector::edge_label_info(
                HyperBasevector::edge_label_info::DIRECT, &edge_names2 ) );
}

// Find edits from the alignments, return pos1, pos2, edits triple, where
// pos1 is the position in source sequence (assembly), pos2 is the position
// in target (reference), and edit is the mutation from refernce to
// assembly.
// Generate insertion edits of unaligned head/tail source sequence
// if head_ins/tail_ins are set. We don't want to call those variants near the
// begin or end of the chromosome or chromosome segments.
void GetEditsFromAlign(const basevector& s, const basevector& t, const align& a, 
        vec<triple<int,int,String>>* edits, vec<pair<String,String>>* change = NULL, 
        char prev_char_s = 'X', char prev_char_t = 'X', bool head_ins = false,
        bool tail_ins = true)
{
    // prefix the two sequences by a leading base
    String ss(1, prev_char_t);
    //String ss(1, prev_char_s);
    ss += s.ToString();
    String tt(1, prev_char_t);
    tt += t.ToString();
    int p1 = a.pos1();
    int p2 = a.pos2();
    if (p1 > 0 && head_ins) { // insertion of p1 seq at beginning
        String alt = ss.substr(0, p1+1);
        String ref = tt.substr(p2,1);
        if (change != NULL)
            change->push(ref, alt);
        edits->push(p1, 0, "Ins: " + ss.substr(1,p1));
    }
    for (int j = 0; j < a.Nblocks(); ++j){
        int gap = a.Gaps(j);
        int len = a.Lengths(j);
        if ( gap > 0 ) {
            String ref = tt.substr(p2,gap+1); // equiv to p2-1 in t
            String alt = ss.substr(p1, 1);    // equiv to p1-1 in s
            alt[0] = ref[0];
            if (change != NULL)
                change->push(ref, alt);
            edits->push(p1, p2, "Del: "+ ToString(basevector(t,p2,gap)));
            p2 += gap;
        }
        if ( gap < 0 ) {
            int ndel = -gap;
            String ref = tt.substr(p2,1);
            String alt = ss.substr(p1, ndel+1);
            alt[0] = ref[0];
            if (change != NULL)
                change->push(ref, alt);
            edits->push(p1, p2, "Ins: " + ToString(basevector(s, p1, ndel)));
            p1 += ndel;
        }
        for ( int x = 0; x < len; x++ ) {
            if ( s[p1] != t[p2] ) {
                String ref = tt.substr(p2+1,1);
                String alt = ss.substr(p1+1,1);
                if (change != NULL)
                    change->push(ref, alt);
                edits->push(p1, p2, "Sub: " + ToString(basevector(s,p1,1)));
            }
            ++p1, ++p2;
        }
    }
    // insertion at the end
    if (p1 != (int)s.size() && tail_ins) {
        int nins = s.size() - p1;
        String ref = tt.substr(p2,1);
        String alt = ss.substr(p1, nins+1);
        alt[0] = ref[0];
        if (change != NULL)
            change->push(ref, alt);
        edits->push((int)s.size(), p2, "Ins: " + ss.substr(p1+1,nins));
    }
}


struct SubEdgeLoc {
    int group_id;
    int branch_id;
    int shift;
    bool dir;
};

} // end namespace


// Given the variants and callers, find the multiple placement of the caller
// edges and genome location (gid, pos) of the alternative placement.
void FindVariantFriends(const vec<VariantCallGroup>& vcall_groups, 
        const vec<vec<align>>& all_aligns, const HyperBasevector& hbp,
        const vec<pair<int,Bool>>& hbp_to_hb,
        map<Variant, vec<pair<int,int>>> *p_var_friending)
{
    map<int, vec<SubEdgeLoc>> edge_homes;
    for (size_t grpid = 0; grpid < vcall_groups.size(); grpid++) {
        const VariantCallGroup group = vcall_groups[grpid];
        for (int branch_id = 0; branch_id < group.GetNBranch(); branch_id++ ) {
            const vec<int>& path = group.GetBranchPath(branch_id);
            int shift = -group.GetTrim(branch_id).first;
            for (size_t sub_id = 0; sub_id < path.size(); sub_id++) {
                int edge0 = hbp_to_hb[path[sub_id]].first;
                bool dir  = hbp_to_hb[path[sub_id]].second;
                edge_homes[edge0].push_back(SubEdgeLoc{(int)grpid, branch_id, shift, dir});
                shift += hbp.EdgeLengthKmers(path[sub_id]);
            }
        }
    }

    // (grpid, caller) pair for each variant and (grpid, vcall) pair for each
    map<Variant, vec<pair<int,Caller>>> vars; 
    for (size_t grpid = 0; grpid < vcall_groups.size(); grpid++) {
        const VariantCallGroup group = vcall_groups[grpid];
        for (int i = 0; i < group.GetNBranch(); i++) {
            const vec<VariantCall>& vcall_vec = group.GetVariantCalls(i);
            for (const VariantCall& x: vcall_vec) {
                int edge0 = hbp_to_hb[x.caller.sub_edge_id].first;
                vars[x.variant].push(grpid, x.caller);
            }
        }
    }
    for (auto it = vars.begin(); it != vars.end(); ++it) {
        const Variant& var = it->first;
        // now assign friend locations
        for (const pair<int,Caller>& x: it->second) {
            int sub_edge = x.second.sub_edge_id;
            int sub_pos  = x.second.sub_edge_pos;
            int edge0    = hbp_to_hb[sub_edge].first;
            bool dir     = hbp_to_hb[sub_edge].second;
            int edge_len = hbp.EdgeLengthBases(sub_edge);

            set<int> grpids; // other groups this sub-edge belongs to
            for (SubEdgeLoc& loc: edge_homes[edge0]) { grpids.insert(loc.group_id); }
            if (grpids.size() <= 1) continue;

            int grpid = x.first;
            for (const SubEdgeLoc& loc: edge_homes[edge0]) { 
                if (loc.group_id == grpid) continue;
                int pos_friend_edge = (loc.dir == dir ? sub_pos : edge_len - 1 - sub_pos);
                int gpos = CorrelatePositionsRangeChecked(
                        vcall_groups[loc.group_id].GetAlign(loc.branch_id),
                        loc.shift + pos_friend_edge );
                if (gpos >= 0)
                    (*p_var_friending)[var].push(vcall_groups[loc.group_id].GetGID(), gpos);
            }
        }
    }
    for (auto it = (*p_var_friending).begin(); it != (*p_var_friending).end(); ++it) {
        UniqueSort(it->second);
    }
}


// ======================= EdgesOnRef method implementation ====================================

void EdgesOnRef::InitFromBestPath(const vec<int>& path, 
        const vec<pair<int,int>>& ref_pos) 
{
    ForceAssert(dg_.N() == 0);
    vec<int> to_right; hb_.ToRight(to_right);
    vec<int> to_left; hb_.ToLeft(to_left);
    dg_.AddVertices(1); // the first vertex
    for (size_t i = 0; i < path.size(); i++) {
        // additional vertex if two edges are not connected in assembly graph
        if (i > 0 && to_right[path[i-1]] != to_left[path[i]])
            dg_.AddVertices(1);
        dg_.AddVertices(1);
        dg_.AddEdge(dg_.N()-2, dg_.N()-1, EdgeLoc(path[i], 
                    ref_pos[i].first + gplus_ext_));
    }
}

static void make_rc_pairs( const HyperBasevector& hb, const vector<int>& list, vector< pair<int,int> >& pair_list)
{
    pair_list.clear();
    vector<int> tmp(list);

    while ( tmp.size() ) {
	int el = tmp.back();
	tmp.pop_back();
	auto hb1 = hb.EdgeObject(el);
	for ( size_t i = 0; i < tmp.size(); ++i ) {
	    auto hb2 = hb.EdgeObject(tmp[i]);
	    if ( hb1 == ReverseComplement(hb2) ) {
		pair_list.push_back( make_pair( el, tmp[i]) );
		tmp.erase(tmp.begin()+i);
	    }
	}
    }
}


// Unroll the graph, find alignments of the unrolled edges, find variants,
// and return the path of the edges.
void EdgesOnRef::UnrollAll(const int verbosity, const int iLoBound, const int iHiBound, const bool bCorrection)
{
    std::unordered_set<int> prefered_edges;
    for(int ii=0;ii<dg_.EdgeObjectCount();++ii){
        prefered_edges.insert(ii);
    }
    set<int> best_path;
    for ( int i = 0; i < dg_.EdgeObjectCount(); i++ )
	best_path.insert(GetOriginalId(i));

    if (verbosity >= 1) {
        cout << Date() << ": Start unrolling using best paths" << endl;
        cout << "    best path (sorted)= ";
        for ( auto edge : best_path )
            cout << edge << " ";
        cout << endl;
        std::cout <<"    unroll limit: " << iLoBound << " " <<iHiBound<<std::endl;
        DumpGraph("hbp.dot", hb_);
    }
    ForceAssertGt(dg_.EdgeObjectCount(), 0);
    const int LenMax = 10000;   // for graph traversal


    vec<int> to_right; hb_.ToRight(to_right);
    vec<int> to_left; hb_.ToLeft(to_left);


    hb_.ToRight(to_right);
    hb_.ToLeft(to_left);
    vec<int> backbone(dg_.EdgeObjectCount(), vec<int>::IDENTITY);
    const int MinOverlap = 80;
    const double MinMatchRate = 0.3;
    const double dRequiredAlignedFractionForExtension=0.2;
    const int64_t iAcceptableFrontExtensionLength= 2 * gplus_ext_;
    const int64_t iAcceptableBackExtensionLength= 2 * (gplus_.size()-genome_.size()-gplus_ext_);

    // record the alternative extension and scores in bubbles edges
    vec<vec<pair<double,int>>> backbone_extension_alts(backbone.size());

    if (verbosity >= 1) cout << Date() << ": Extending the path tail" << endl;
    int tail_edge = dg_.EdgeObjectCount() - 1;
    while (true) {
        int tail_edge0 = GetOriginalId(tail_edge);
        int tail_edge_end = hb_.EdgeLengthKmers(tail_edge0) + GetStart(tail_edge);
        if ( size_t(tail_edge_end) >= gplus_.size()){ break;}
        if (verbosity >= 2) 
            cout << "Trailing edge " << tail_edge << "->" << tail_edge0 << endl;
        if (hb_.From(to_right[tail_edge0]).empty()) break;

        int new_edge0 = -1;
        double best_match_rate = 0.0;
        align a;
        vec<pair<double,int>> alt_edges;
        for (int i = 0; i < hb_.FromSize(to_right[tail_edge0]); i++) {
            int eid = hb_.EdgeObjectIndexByIndexFrom(to_right[tail_edge0], i);
            align ai;
            int Bandwidth = max(1, (int)hb_.EdgeObject(eid).size()/2);
            RestrictedAlign(hb_.EdgeObject(eid), gplus_, tail_edge_end,
                    Bandwidth, ai );
            int nmatch = ai.MatchingBases(hb_.EdgeObject(eid), gplus_);
            double match_rate = (double) nmatch / 
                min(hb_.EdgeObject(eid).size(), gplus_.size()-tail_edge_end);

            if (match_rate > best_match_rate) {
                best_match_rate = match_rate;
                new_edge0 = eid;
                a = ai;
            }
            alt_edges.push(match_rate, eid);
        }
        ReverseSort(alt_edges);

        bool is_loop = any_of(backbone.begin(), backbone.end(), 
                [=](int x){return GetOriginalId(x) == new_edge0;});
        if (is_loop) {
            if (verbosity >= 2) cout << "loop encountered" << endl;
            break;
        }
        int new_edge_start = a.pos2() - a.pos1();
        if (  new_edge_start + MinOverlap > (int)gplus_.size()
                || best_match_rate < MinMatchRate
            || a.pos2() > iHiBound || a.Pos2() > iHiBound
           ) break;

        if (verbosity >= 2) {
            cout << "align at " << a.pos1() << "," << a.Pos1() << " to " << a.pos2() << "," << a.Pos2() << endl;
            if ( a.pos2() - a.pos1()+ MinOverlap < (int)gplus_.size() )
            PrintVisualAlignment(True, cout, hb_.EdgeObject(new_edge0), gplus_, a);
        }
        if(   dRequiredAlignedFractionForExtension * hb_.EdgeObject(new_edge0).size() > a.Pos1()-a.pos1()
           && hb_.EdgeObject(new_edge0).size() > iAcceptableBackExtensionLength
          ){
            break;
        }

        if (verbosity >= 2) 
            cout << "Adding tail " << new_edge0 << " start " << new_edge_start << endl;
        vec<pair<int,int>> tail_extension;
        tail_extension.push(new_edge0, new_edge_start );
        AddConnections(tail_edge, -1, tail_extension);
        tail_edge = dg_.EdgeObjectCount() - 1;
        backbone.push_back(tail_edge);
        if (!IsBubbleEdge(new_edge0, hb_, to_left, to_right)) 
            alt_edges.clear();
        backbone_extension_alts.push_back(alt_edges);
    }

    if (verbosity >= 1) cout << Date() << ": Extending the path head" << endl;
    int head_edge = 0;
    while (true) {
        int head_edge0 = GetOriginalId(head_edge);
        int start0 = GetStart(head_edge);
        if ( start0 < 0){ break; }
        if (verbosity >= 2) 
            cout << "head edge " << head_edge0 << " start at " << start0 << endl;
        if(hb_.To(to_left[head_edge0]).empty()) break;

        int new_edge0 = -1; 
        double best_match_rate = 0.0;
        align a;
        vec<pair<double,int>> alt_edges;
        for (int i = 0; i < hb_.ToSize(to_left[head_edge0]); i++) {
            int eid = hb_.EdgeObjectIndexByIndexTo(to_left[head_edge0], i);
            align ai;
            int Bandwidth = max(1, (int)hb_.EdgeObject(eid).size()/2);
            int start_prev = start0 - hb_.EdgeLengthKmers(eid);
            RestrictedAlign(hb_.EdgeObject(eid), gplus_, start_prev, Bandwidth, ai);
            int nmatch = ai.MatchingBases(hb_.EdgeObject(eid), gplus_);
            double match_rate = (double) nmatch / min(hb_.EdgeObject(eid).isize(), 
                                                      start0);
            if (match_rate > best_match_rate) {
                best_match_rate = match_rate;
                new_edge0 = eid;
                a = ai;
            }
            alt_edges.push(match_rate, eid);
        }
        ReverseSort(alt_edges);

        bool is_loop = any_of(backbone.begin(), backbone.end(), 
                [=](int x){return GetOriginalId(x) == new_edge0;});
        if (is_loop) {
            if (verbosity >= 2) 
                cout << "loop encountered" << endl;
            break;
        }
        if (a.Pos2() < MinOverlap || best_match_rate < MinMatchRate
            || a.pos2() < iLoBound || a.Pos2() < iLoBound
           ) break;
        if (verbosity >= 2) {
            cout << "align at " << a.pos1() << "," << a.Pos1() << " to " << a.pos2() << "," << a.Pos2() << endl;
            PrintVisualAlignment(True, cout, hb_.EdgeObject(new_edge0), gplus_, a);
        }
        if(   dRequiredAlignedFractionForExtension * hb_.EdgeObject(new_edge0).size() > a.Pos1()-a.pos1()
           && hb_.EdgeObject(new_edge0).size() > iAcceptableFrontExtensionLength
          ){
            break;
        }

        int new_edge_start = a.pos2() - a.pos1();
        vec<pair<int,int>> head_extension;
        head_extension.push(new_edge0, new_edge_start);
        if (verbosity >= 2) 
            cout << "Adding head " << new_edge0 << " starting at " << new_edge_start << endl;
        AddConnections(-1, head_edge, head_extension);
        head_edge = dg_.EdgeObjectCount() - 1;
        backbone.push_front(head_edge);
        if (!IsBubbleEdge(new_edge0, hb_, to_left, to_right)) 
            alt_edges.clear();
        backbone_extension_alts.push_front(alt_edges);
    }

    ///////////////////////////////////////////////////////////////////
    // START OF INVERSION CODE
    // try to remove inversion edges not on the best path
    if ( verbosity >= 1 ) cout << Date() << ": looking for inversions" << endl;

    vec<int> to_inv(hb_.EdgeObjectCount(), -1);
    for (int e = 0; e < hb_.EdgeObjectCount(); e++) {
        basevector base_inv = ReverseComplement(hb_.EdgeObject(e));
        for (int e2 = e+1; e2 < hb_.EdgeObjectCount(); e2++) {
            const basevector& base2 = hb_.EdgeObject(e2);
            if (base_inv.size() == base2.size() &&
                base_inv == base2) {
                to_inv[e] = e2;
                to_inv[e2] = e;
            }
        }
    }

    // Fix the backboen extension if theire are inversion bubbles
    {
        vec<int> backbone0 = ConvertToAssemblyPathFromUnrolled(backbone);
        for (size_t i = 0; i < backbone0.size(); i++) {
            double best_match_rate = 1.0;
            if (!backbone_extension_alts[i].empty())
                best_match_rate = backbone_extension_alts[i][0].first;
            for (size_t j = 0; j < backbone0.size(); j++) {
                if( j != i && to_inv[backbone0[i]] == backbone0[j]
                        && !backbone_extension_alts[j].empty()){
                   double  best_match_rate2 = backbone_extension_alts[j][0].first;
                   if (best_match_rate2 < best_match_rate && 
                           backbone_extension_alts[j].size() >= 2) {
                       int edge_alt =  backbone_extension_alts[j][1].second; // second best
                       EdgeLoc& edge_j = dg_.EdgeObjectMutable(backbone[j]);
                       edge_j.original_edge = edge_alt;
                       if (verbosity >= 2)
                           cout << "Replacing backbone edge " << backbone[j] << " from " 
                               << backbone0[j] << " to " << edge_alt << endl;
                       backbone0[j] = edge_alt;
                       backbone_extension_alts[j].erase(backbone_extension_alts[j].begin());
                   }
                }
            }
        }
    }

    // Divide the graph into compoments
    vec<int> component_ids(hb_.EdgeObjectCount(), -1);
    {
        vec<vec<int>> component_edges;
        hb_.ComponentEdges(component_edges);
        for (size_t i = 0; i < component_edges.size(); i++)
            for (int e: component_edges[i])
                component_ids[e] = i;

        // join compoments that are linked by backbone path
        int comp_id0 = component_ids[GetOriginalId(backbone[0])];
        for (size_t i = 1; i < backbone.size(); i++) {
            int comp_id = component_ids[GetOriginalId(backbone[i])];
            if (comp_id != comp_id0) {
                if (verbosity >= 2) {
                    cout << "joining component " << comp_id << " (" << 
                        component_edges[comp_id].size() << " edges) to " << comp_id0 << endl;
                }
                for (size_t j = 0; j < component_edges[comp_id].size(); j++) {
                    int e = component_edges[comp_id][j];
                    component_ids[e] = comp_id0;
                }
            }
            copy(component_edges[comp_id].begin(), component_edges[comp_id].end(),
                    back_inserter(component_edges[comp_id0]));
            component_edges[comp_id].clear();
        }
    }

    // get edge object lengths for quick comparison
    vector<int> edge_length( hb_.EdgeObjectCount() );
    map<int, vector<int> > length_to_edgeid;
    for ( int e = 0; e < hb_.EdgeObjectCount(); e++ ) {
        size_t edge_length = hb_.EdgeObject(e).size();
        length_to_edgeid[edge_length].push_back(e);
    }
    if ( verbosity >= 2 ) {
        cout << Date() << ": all edge lengths" << endl;
        for ( auto el : length_to_edgeid ) {
            cout << el.first << " --> ";
            for ( auto l : el.second )
        	cout << l << " ";
            cout << endl;
        }
    }
    // new algorithm
    // 1. find all pairs of RC edges
    // 2. if one is one best path and the other is not, kill it
    // 3. if one is adjacent to best path and the other is not, kill it
    // cout << "BEST PATH=";
    // for ( auto best : best_path )
    //     cout << best << " ";
    // cout << endl;
    vec<int> to_delete;
    for ( auto el : length_to_edgeid ) {
        auto list = el.second;
        vector< pair<int, int> > rc_pairs;
        make_rc_pairs(hb_, list, rc_pairs);
        for ( auto pair : rc_pairs ) {
            int p1 = pair.first;
            int p2 = pair.second;
            bool p1_best = best_path.find(p1) != best_path.end();
            bool p2_best = best_path.find(p2) != best_path.end();
            // PRINT4(p1,p1_best,p2,p2_best);
            if ( p1_best && !p2_best ) to_delete.push_back(p2);
            else if ( p2_best && !p1_best ) to_delete.push_back(p1);
            else if ( !p1_best && !p2_best ) {
        	// look for 2nd order effects
        	bool p1_adj_best = false;
        	bool p2_adj_best = false;
        	for ( int i = 0; i < hb_.FromSize( to_left[p1] ) && !p1_adj_best; i++ ) {
        	    int eid = hb_.EdgeObjectIndexByIndexFrom(to_left[p1], i);
        	    if ( best_path.find(eid) != best_path.end() ) p1_adj_best = true;
        	}
        	for ( int i = 0; i < hb_.ToSize( to_right[p1] ) && !p1_adj_best; i++ ) {
        	    int eid = hb_.EdgeObjectIndexByIndexTo(to_right[p1], i);
        	    if ( best_path.find(eid) != best_path.end() ) p1_adj_best = true;
        	}
        	for ( int i = 0; i < hb_.FromSize( to_left[p2] ) && !p2_adj_best; i++ ) {
        	    int eid = hb_.EdgeObjectIndexByIndexFrom(to_left[p2], i);
        	    if ( best_path.find(eid) != best_path.end() ) p2_adj_best = true;
        	}
        	for ( int i = 0; i < hb_.ToSize( to_right[p2] ) && !p2_adj_best; i++ ) {
        	    int eid = hb_.EdgeObjectIndexByIndexTo(to_right[p2], i);
        	    if ( best_path.find(eid) != best_path.end() ) p2_adj_best = true;
        	}

        	if ( p1_adj_best && !p2_adj_best ) to_delete.push_back(p2);
        	else if ( p2_adj_best && !p1_adj_best ) to_delete.push_back(p1);
            }
        }
    }

    vec<int> to_del2;
    // More aggressively remove inversion edges that belongs to the same graph
    // components.

    // remove backbone branches if there are inversion mixtures. However, keep
    // the branch if the sequences is similar to backbone.
    for (size_t i = 0; i < backbone.size(); i++) {
        int edge0 = GetOriginalId(backbone[i]);
        int node_left = to_left[edge0], node_right = to_right[edge0];
        bool rc_mixture = false, is_bubble_edge = true;
        for (int j = 0; j < hb_.FromSize(node_left); j++) {
            int edge_alt = hb_.EdgeObjectIndexByIndexFrom(node_left, j);
            if (component_ids[to_inv[edge_alt]] == component_ids[edge0]) {
                rc_mixture = true;
            }
            if (to_right[edge_alt] != node_right) {
                is_bubble_edge = false;
                break;
            }
        }
        if (! rc_mixture || ! is_bubble_edge) continue;
        vec<int> edge_ids;
        vec<basevector> edges;
        for (int j = 0; j < hb_.FromSize(node_left); j++) {
            int edge_id= hb_.EdgeObjectIndexByIndexFrom(node_left, j);
            edge_ids.push_back(edge_id);
            edges.push_back(hb_.EdgeObject(edge_id));
        }
        vec<int> similar_edges;
        //cout << "edge0= " << edge0 << endl;
        FindSimilar(hb_.EdgeObject(edge0), edges, edge_ids, similar_edges);
        copy(similar_edges.begin(), similar_edges.end(), back_inserter(to_del2));
    }

    // Always remove reversed backbone edges in the branch
    auto RemoveBackboneReverse = [&](const vec<int>& backbone00, vec<int>& to_del3) {
        vec<int> backbone_reverted00;
        for (auto it = backbone00.rbegin(); it != backbone00.rend(); ++it) {
            backbone_reverted00.push_back(to_inv[*it]);
        }
        if (verbosity >= 2)
        cout << "Reverted backbone: ";
        backbone_reverted00.Println(cout);
        int v = to_left[backbone_reverted00[0]];
        for (size_t i = 0; i < backbone_reverted00.size(); i++) {
            if (v == -1)
                v = to_left[backbone_reverted00[i]];
            int expect = backbone_reverted00[i];
            int next_v = -1;
            for (int j = 0; j < hb_.FromSize(v); j++) {
                int e2 = hb_.EdgeObjectIndexByIndexFrom(v, j);
                if (e2 == expect) {
                    if (!Member(backbone00, e2)) {
                        to_del3.push_back(e2);
                    }
                    next_v = hb_.From(v)[j];
                    break;
                }
            }
            v = next_v;
        }
    };
    vec<int> backbone00 = ConvertToAssemblyPathFromUnrolled(backbone);
    vec<int> to_del3;
    RemoveBackboneReverse(backbone00, to_del3);
    if (verbosity >= 2) {
        cout << "Delete additional inversion edges that belongs to the same components: edge_id = ";
        to_del2.Println(cout);
        cout << "Delete additional inversion edges that belongs to the reverted backbone: edge_id = ";
        to_del3.Println(cout);
    }
    copy(to_del2.begin(), to_del2.end(), back_inserter(to_delete));
    copy(to_del3.begin(), to_del3.end(), back_inserter(to_delete));
    //// backbone extension
    //vec<int> backbone00_ext = backbone00;
    //while (1) {
    //    int e0 = backbone00_ext[0];
    //    int vleft = to_left[e0];
    //    vec<int> prev_edges;
    //    for (int j = 0; j < hb_.ToSize(vleft); j++) {
    //        int e1 = hb_.EdgeObjectIndexByIndexTo(vleft, j);
    //        if (!Member(to_del3, e1))
    //            prev_edges.push_back(e1);
    //    }
    //    if (prev_edges.size() != 1) 
    //        break;
    //    else
    //        backbone00_ext.push_front(prev_edges[0]);
    //}
    //vec<int> to_del4;
    //RemoveBackboneReverse(backbone00_ext, to_del4);
    //copy(to_del4.begin(), to_del4.end(), back_inserter(to_delete));
    //if (verbosity >= 2) {
    //    cout << "Delete additional inversion edges that belongs to extended reverted backbone: edge_id = ";
    //    to_del4.Println(cout);
    //}
    UniqueSort(to_delete);

    // okay, delete edges
    hb_.DeleteEdges(to_delete);
    if ( verbosity >= 2 ) {
        cout << Date() << ": inversion edges to delete:" << endl;
        for ( auto edge : to_delete ) cout << edge << " ";
        cout << endl;
        DumpGraph("hbp-noinv.dot", hb_);
    }

    // END OF INVERSION CODE
    ///////////////////////////////////////////////////////////////////
 
    if (verbosity >= 2) {
        cout << "The backbone edges are: ";
        backbone.Println(cout);
        cout << "The corresponding hbp edges are: ";
        vec<int> backbone0 = ConvertToAssemblyPathFromUnrolled(backbone);
        for (size_t i = 0; i < backbone0.size(); i++) {
            if (i!=0) cout << " ";
            cout << backbone0[i];
            if (Member(backbone0, to_inv[backbone0[i]]) )
                cout << "("<< to_inv[backbone0[i]] <<")";
        }
        cout << endl;
    }
    // each edge belongs to a path between two edges [e1,.., e2], where e1
    // and e2 are all in the backbone
    map<int, pair<int,int>> bb_anchors; 
    for (auto& x: backbone) 
        bb_anchors[x] = make_pair(x,x);


    if (verbosity >= 1) 
        cout << Date() << ": Looking for alternative paths between backbone edges" << endl;
    vec<int> backbone_dup(backbone.size(), 0);
    vec<int> backbone0 = ConvertToAssemblyPathFromUnrolled(backbone);
    const int MaxLocDevAllowed = 150;
    while (true) {
        bool no_change = true; // keep searching if alternative path is found
        vec<vec<pair<int,int>>> known_starts(hb_.EdgeObjectCount());
        for (int e2 = 0; e2 < dg_.EdgeObjectCount(); e2++) 
            known_starts[GetOriginalId(e2)].push(GetStart(e2), e2);

        for (int node2 = 0; node2 < dg_.N() && no_change; node2++) {
            for (size_t j = 0; j < dg_.From(node2).size(); j++) {
                // For each edge
                int edge2 = dg_.EdgeObjectIndexByIndexFrom(node2, j);
                int edge1 = GetOriginalId(edge2);
                int node2_next = dg_.From(node2)[j];
                int node1_next = to_right[edge1];
                // start from node1_next in graph1 and find a new path from
                // known paths starting from node_next in graph2 
                set<int> known_edges1;
                for (size_t k = 0; k < dg_.From(node2_next).size(); k++) {
                    int edge2_next = dg_.EdgeObjectIndexByIndexFrom(node2_next, k);
                    known_edges1.insert(GetOriginalId(edge2_next));
                }
                for (size_t k = 0; k < hb_.From(node1_next).size() && no_change; k++) {
                    int edge1_next = hb_.EdgeObjectIndexByIndexFrom(node1_next, k);
                    int node1_nn = hb_.From(node1_next)[k];
                    if (find(known_edges1.begin(), known_edges1.end(), edge1_next)
                            != known_edges1.end()) continue;
                    if (Member(backbone0, edge1_next)) continue;

                    if (verbosity >= 2)
                        cout << "try path "<< edge2 << "(" << edge1 << ") -> (" 
                            << edge1_next << ") " << endl;

                    // Recursively follow the edge until reaching in backbone again
                    // search for possible alternative paths from edge1 to back_edge1:
                    //    (xx) --edge1--> (node1_next) --edge1_next--> (node1_nn) --back_edge1--
                    // (node2) --edge2--> (node2_next) --edge2_next--> (xx)       --back_edge2--
                    int back_edge = -1;
                    vec<pair<int,int>> alternative_path_and_start; 
                    int start_next = GetStart(edge2) + hb_.EdgeLengthKmers(edge1);
                    alternative_path_and_start.push(edge1_next, start_next);
                    if (hb_.To(node1_next).solo() && hb_.From(node1_nn).solo() &&
                            known_edges1.size() == 1 && 
                            to_right[*known_edges1.begin()] == node1_nn) { // alway add simple bubble
                        int edge1_nn = hb_.EdgeObjectIndexByIndexFrom(node1_nn, 0);
                        if (known_starts[edge1_nn].solo()) {
                            back_edge = known_starts[edge1_nn][0].second;
                        }
                    }
                    if (back_edge == -1) {
                        FindAlternativePath(node1_nn, &alternative_path_and_start, &back_edge, 
                                known_starts, LenMax, MaxLocDevAllowed, verbosity);
                    }
                    if (back_edge != -1) {
                        if (verbosity >= 2)
                        {
                            cout << "found alternative path between " << edge2 << "(" << edge1 << ")" << " and " 
                                << back_edge << "(" << GetOriginalId(back_edge) << ")" << endl;
                            for (auto&x: alternative_path_and_start) 
                                cout << x.first << "@" << x.second << " ";
                            cout << endl;
                        }
                        int new_edge_start = dg_.EdgeObjectCount();
                        AddConnections(edge2, back_edge, alternative_path_and_start);
                        int new_edge_end = dg_.EdgeObjectCount();

                        int anchor_start = bb_anchors[edge2].first;
                        int anchor_end = bb_anchors[back_edge].second;
                        for (int x = new_edge_start; x <= new_edge_end; x++) 
                            bb_anchors[x] = make_pair(anchor_start, anchor_end);
                        int pos1 = find(backbone.begin(), backbone.end(), anchor_start) - backbone.begin();
                        int pos2 = find(backbone.begin(), backbone.end(), anchor_end) - backbone.begin();
                        ForceAssertNe(pos1, backbone.isize());
                        ForceAssertNe(pos2, backbone.isize());
                        if (pos2 > pos1) {
                            for (int x = pos1+1; x < pos2; x++)
                                backbone_dup[x]++;
                        } else { // loop back
                            for (int x = pos2; x <= pos1; x++)
                                backbone_dup[x]++;
                        }
                        no_change = false;
                    }
                }
            }
        }
        if (no_change) break;
    }


    if(verbosity >= 2){
        std::cout<< "prefered edges ";
        for( const auto& entry: prefered_edges){
            std::cout << entry << " " ;
        }
        std::cout<< std::endl;
    }



for( bool bLoop=bCorrection;bLoop;)
{
    bLoop=false;
    auto dups = count_dg_dups();
    if( dups.size() > 0){
        vec<int> dg_to_right; dg_.ToRight(dg_to_right);
        vec<int> dg_to_left; dg_.ToLeft(dg_to_left);
        for( auto& entry: dups){
            std::cout << std::get<0>(entry) << " "
                      << std::get<1>(entry) << " "
                      << std::get<2>(entry) << " "
                      << std::get<3>(entry) << " "
                      << std::endl;
        }
        auto connection_counts = count_dg_prefered_edge_connection(prefered_edges,dg_to_left,dg_to_right);
//        std::cout <<"connection count: "<< std::endl;
//        for(size_t ii=0;ii<connection_counts.size();++ii){
//            std::cout << ii << " " << connection_counts[ii].first << " " << connection_counts[ii].second << std::endl;
//        }
        vec<std::tuple<size_t,int,size_t,bool,int,int,int,int>> dups_counted(dups.size());
        for(size_t ii=0;ii<dups.size();++ii){
            const auto& entry= dups[ii];
            auto vertex = std::get<0>(entry) ;
            auto vi = std::get<1>(entry) ;
            auto vj = std::get<2>(entry) ;
            auto dir = std::get<3>(entry) ;

            if( dir==0){ // to
                auto ei = dg_.EdgeObjectIndexByIndexTo(vertex,vi);
                auto ej = dg_.EdgeObjectIndexByIndexTo(vertex,vj);
                ForceAssert( dg_.EdgeObject(ei).original_edge == dg_.EdgeObject(ej).original_edge);

                size_t ci = connection_counts[ei].first;
                size_t cj = connection_counts[ej].first;
                bool bji = cj>ci;
                int loser = bji?ei:ej;
                dups_counted[ii]=std::make_tuple(max(ci,cj),loser,min(ci,cj),bji,vertex,vi,vj,dir);
            }
            else{ // from
                auto ei = dg_.EdgeObjectIndexByIndexFrom(vertex,vi);
                auto ej = dg_.EdgeObjectIndexByIndexFrom(vertex,vj);
                ForceAssert( dg_.EdgeObject(ei).original_edge == dg_.EdgeObject(ej).original_edge);

                size_t ci = connection_counts[ei].second;
                size_t cj = connection_counts[ej].second;
                bool bji = cj>ci;
                int loser = bji?ei:ej;
                dups_counted[ii]=std::make_tuple(max(ci,cj),loser,min(ci,cj),bji,vertex,vi,vj,dir);
            }
        }
        Sort(dups_counted);
        for( auto& entry: dups_counted){
            auto c1 = std::get<0>(entry) ;
            auto loser = std::get<1>(entry);
            auto c2 = std::get<2>(entry) ;
            auto bji = std::get<3>(entry) ;
            auto vertex = std::get<4>(entry) ;
            auto vi = std::get<5>(entry) ;
            auto vj = std::get<6>(entry) ;
            auto dir = std::get<7>(entry) ;

            auto ci = (bji)?c2:c1;
            auto cj = (bji)?c1:c2;

            int ei,ej;
            if( dir==0){ // to
                ei = dg_.EdgeObjectIndexByIndexTo(vertex,vi);
                ej = dg_.EdgeObjectIndexByIndexTo(vertex,vj);
            }
            else{
                ei = dg_.EdgeObjectIndexByIndexFrom(vertex,vi);
                ej = dg_.EdgeObjectIndexByIndexFrom(vertex,vj);
            }
            ForceAssert( dg_.EdgeObject(ei).original_edge == dg_.EdgeObject(ej).original_edge);
            std::cout << dg_.EdgeObject(ei).original_edge << " "
                      << dir << "|"
                      << loser << " "
                      << ei << " "
                      << ej << " "
                      << ci << " "
                      << cj << " "<< std::endl;
        }
        vec<int> to_delete;
        clear_dg_branch(to_delete, std::get<1>(dups_counted.back()),std::get<7>(dups_counted.back()),dg_to_left,dg_to_right  );
        Sort(to_delete);
        for(auto& entry: backbone){
            int count = 0;
            for(count=0;count<entry;++count){ }
            ForceAssert(count<=entry);
            entry-=count;
        }
        dg_.DeleteEdges(to_delete);
        bLoop=true;
    }
}

    if (verbosity >= 2)
    {
        cout << "Seeding edges are: " ;
        for (size_t i = 0; i < backbone.size(); i++) {
            if (backbone_dup[i] == 0) cout << backbone[i] << " ";
        }
        cout << endl;
    }
    backbone_ = backbone;
    backbone_dup_ = backbone_dup;
};


vec<std::tuple<int,int,int,int>> EdgesOnRef::count_dg_dups()const{
    vec<std::tuple<int,int,int,int>> out;
    for (int vertex = 0; vertex < dg_.N()  ; vertex++) {
        for (size_t ii = 0; ii < dg_.To(vertex).size(); ii++) {
            auto ei_i = dg_.EdgeObjectIndexByIndexTo(vertex,ii);
            for (size_t jj = ii+1; jj < dg_.To(vertex).size(); jj++) {
                auto ei_j = dg_.EdgeObjectIndexByIndexTo(vertex,jj);
                if( dg_.EdgeObject(ei_i).original_edge == dg_.EdgeObject(ei_j).original_edge){
                    std::cout <<  dg_.EdgeObject(ei_i).original_edge << std::endl;
                    out.push_back( std::make_tuple(vertex,ii,jj,0) );
                }
            }
        }
        for (size_t ii = 0; ii < dg_.From(vertex).size(); ii++) {
            auto ei_i = dg_.EdgeObjectIndexByIndexFrom(vertex,ii);
            for (size_t jj = ii+1; jj < dg_.From(vertex).size(); jj++) {
                auto ei_j = dg_.EdgeObjectIndexByIndexFrom(vertex,jj);
                if( dg_.EdgeObject(ei_i).original_edge == dg_.EdgeObject(ei_j).original_edge){
                    std::cout <<  dg_.EdgeObject(ei_i).original_edge << std::endl;
                    out.push_back( std::make_tuple(vertex,ii,jj,1) );
                }
            }
        }
    }
    return out;
}
vec<std::pair<size_t,size_t>> EdgesOnRef::count_dg_prefered_edge_connection(const std::unordered_set<int>& prefered_edges
                                                                           ,const vec<int>& to_left, const vec<int>& to_right)const{
    vec<std::pair<size_t,size_t>> out(dg_.EdgeObjectCount(),std::make_pair(0ul,0ul));
    for( auto pe: prefered_edges){
        //going right
        const size_t nBases = hb_.EdgeObject(GetOriginalId(pe)).size();
        {
            std::deque<int> to_do;
            std::set<int> discovered;
            to_do.push_back(pe);
            discovered.insert(pe);
            while(!to_do.empty()){
                auto edge=to_do.front();
                to_do.pop_front();

                out[edge].first+=nBases;
                auto w = to_right[edge];
                for( int ee=0 ; ee < dg_.FromSize(w);++ee){
                    auto ei = dg_.EdgeObjectIndexByIndexFrom(w,ee);
                    if( discovered.find(ei) == discovered.end()){
                        discovered.insert(ei);
                        to_do.push_back(ei);
                    }
                }
            }
        }
        //going left
        {
            std::deque<int> to_do;
            std::set<int> discovered;
            to_do.push_back(pe);
            discovered.insert(pe);
            while(!to_do.empty()){
                auto edge=to_do.front();
                to_do.pop_front();

                out[edge].second+=nBases;
                auto v = to_left[edge];
                for( int ee=0 ; ee < dg_.ToSize(v);++ee){
                    auto ei = dg_.EdgeObjectIndexByIndexTo(v,ee);
                    if( discovered.find(ei) == discovered.end()){
                        discovered.insert(ei);
                        to_do.push_back(ei);
                    }
                }
            }
        }
    }
    return out;
}

void EdgesOnRef::clear_dg_branch(vec<int>& to_delete, int root_edge,int dir,const vec<int>& to_left,const vec<int>&to_right){
    to_delete.clear();
    std::set<int> deleted;
    std::deque<int> to_do;
    to_do.push_back(root_edge);
    while(!to_do.empty()){
        int edge = to_do.front();
        to_do.pop_front();
        if( deleted.find(edge)== deleted.end()){
            deleted.insert(edge);
            to_delete.push_back(edge);


            size_t nParallelAway=0;
            if( dir == 0){
                int vertex = to_left[edge];
                for(int ee=0;ee<dg_.FromSize(vertex);++ee){
                    if( deleted.find(dg_.EdgeObjectIndexByIndexFrom(vertex,ee))==deleted.end()){
                        ++nParallelAway;
                    }
                }
            }
            else{
                int vertex = to_right[edge];
                for(int ee=0;ee<dg_.ToSize(vertex);++ee){
                    if( deleted.find(dg_.EdgeObjectIndexByIndexTo(vertex,ee))==deleted.end()){
                        ++nParallelAway;
                    }
                }
            }
            if(nParallelAway==0){
                if( dir == 0){
                    int vertex = to_left[edge];
                    for(int ee=0;ee<dg_.ToSize(vertex);++ee){
                        to_do.push_back( dg_.EdgeObjectIndexByIndexTo(vertex,ee));
                    }
                }
                else{
                    int vertex = to_right[edge];
                    for(int ee=0;ee<dg_.FromSize(vertex);++ee){
                        to_do.push_back( dg_.EdgeObjectIndexByIndexFrom(vertex,ee));
                    }
                }
            }
        }
    }
    Sort(to_delete,std::greater<int>());
    /*
    std::cout<<"deleting: ";
    for(auto entry:to_delete){
        std::cout << entry << " ";
    }
    std::cout<<std::endl;
    dg_.DeleteEdges(to_delete);
    */
}

// Starting from an edge, exploring all path until it lands to the unrolled graph and
// the inferred edge location is consistent with know value.
bool EdgesOnRef::FindAlternativePath(int node, vec<pair<int,int>> *p_edge_starts, int *p_back_edge, 
        const vec<vec<pair<int,int>>>& known_starts, 
        const int LenMax, const int MaxLocDevAllowed,
        int verbosity) 
{
    ForceAssertGt(p_edge_starts->size(), 0u);
    for (size_t i = 0; i < hb_.From(node).size(); i++) {
        int node_new = hb_.From(node)[i];
        int edge_new = hb_.EdgeObjectIndexByIndexFrom(node, i);
        int start_new = p_edge_starts->back().second + 
            hb_.EdgeLengthKmers(p_edge_starts->back().first) ;
        // do not allow loop
        bool is_loop = false;
        for (auto& x: *p_edge_starts) {
            if (x.second == edge_new) {
                is_loop = true;
                break;
            }
        }
        if (is_loop) continue;
        // check if it's a valid alternative path, either if the path end
        // shift agrees with known one, or the path belongs to a simple bubble
        const vec<pair<int,int>>& starts = known_starts[edge_new];
        *p_back_edge = -1; // find the tail edge?
        for (const auto& x : starts) {
            if (abs(x.first - start_new) < MaxLocDevAllowed) {
                *p_back_edge = x.second;
                break;
            }
        }
        if (*p_back_edge != -1) // path end agrees
            return true;
        // not found, track down
        if (start_new - p_edge_starts->front().second > LenMax
                || starts.size() > 0 || p_edge_starts->size() > 10) {
            if (starts.size() > 0 && verbosity >= 2) {
                cout << "skip mismatch edge  " << edge_new << " expect_start= " << start_new 
                    << " known_start= " ;
                for (const auto& x : starts) cout << x.first << " "; cout << endl;
            }
            continue;
        }
        p_edge_starts->push_back(make_pair(edge_new,start_new));
        //cout << "follow edge ->" << edge_new << endl;
        if (FindAlternativePath(node_new, p_edge_starts, p_back_edge, known_starts, 
                    LenMax, MaxLocDevAllowed, verbosity))
            return true;
    }
    // all branches failed, clear last node, track back.
    p_edge_starts->pop_back();
    return false;
}

void EdgesOnRef::AddConnections(int front_edge, int back_edge, const vec<pair<int,int>>& alternative_path_and_start)
{
    vec<int> to_right; dg_.ToRight(to_right);
    vec<int> to_left; dg_.ToLeft(to_left);
    // add the missing node if needed
    int last_node = (back_edge == -1 ? dg_.N() : to_left[back_edge]);
    if (back_edge == -1) dg_.AddVertices(1);
    int first_node = (front_edge == -1 ? dg_.N() : to_right[front_edge]);
    if (front_edge == -1) dg_.AddVertices(1);
    if (alternative_path_and_start.solo()) {
        dg_.AddEdge(first_node, last_node, EdgeLoc(alternative_path_and_start[0]));
    } else {
        int start_new_vertices = dg_.N();
        dg_.AddVertices(alternative_path_and_start.size() - 1);
        dg_.AddEdge(first_node, start_new_vertices, 
                EdgeLoc(alternative_path_and_start.front()));
        for (size_t i = 1; i < alternative_path_and_start.size() - 1; ++i) {
            dg_.AddEdge(start_new_vertices + i-1, start_new_vertices+i, 
                    EdgeLoc(alternative_path_and_start[i]));
        }
        dg_.AddEdge(dg_.N()-1, last_node, EdgeLoc(alternative_path_and_start.back()));
    }
}

// Align each edge from bubble graph to reference and call variants.
// Also report the alignment of these edges, which are used to find friend
// locations.
void EdgesOnRef::CallVariantsGroupedWithProb(int gid, vec<VariantCallGroup> *p_groups, 
        vec<align>* p_edge_aligns, int verbosity) 
{
    if (verbosity >=1) cout << Date() << ": Start calling variants" << endl;
    // identify anchoring edges from bubble graph
    vec<int> anchoring_edges;
    vec<pair<bool,bool>> free_ends;
    for (int n1 = 0; n1 < bubble_graph_.N(); n1++) {
        if (bubble_graph_.From(n1).size() == 1) {
            int n2 = bubble_graph_.From(n1)[0];
            int e = bubble_graph_.EdgeObjectIndexByIndexFrom(n1,0);
            ForceAssertEq(bubble_graph_.To(n2).size(), 1u);
            anchoring_edges.push_back(e);
            free_ends.push(bubble_graph_.To(n1).empty(),
                    bubble_graph_.From(n2).empty());
        }
    }
    int nedges = bubble_graph_.EdgeObjectCount();
    int ngroups = anchoring_edges.size() * 2 - 1;
    ForceAssert( ngroups >= 0 );
    // if ( ngroups < 0 ) ngroups = 0;

    // to be calculated 
    vec<VariantCallGroup> all_groups(ngroups);
    vec<align>   align_all_edges(nedges);
    for (size_t i = 0; i < anchoring_edges.size(); i++) {
        int e = anchoring_edges[i];
        const vec<int>& path = bubble_graph_edge_path_[e];
        vec<int> path0 = ConvertToAssemblyPathFromUnrolled(path);
        basevector edge_seq = bubble_graph_.EdgeObject(e);
        // trim K-1 base if overlaping with cells
        int trim_start = (free_ends[i].first ? 0 : hb_.K()-1);
        int trim_end = (free_ends[i].second ? 0 : hb_.K()-1); 

        char edge_prev_base = ( trim_start == 0 ? 'X' : 
                as_base(edge_seq[trim_start-1]) );
        edge_seq.SetToSubOf(edge_seq, trim_start, 
                edge_seq.size() - trim_start - trim_end);
        int path_offset = GetStart(path.front()) + trim_start;

        align a;
        int bandwidth = max(1, (int)edge_seq.size() /2 );

        RestrictedAlign( edge_seq, gplus_, path_offset, bandwidth, a );
        if (verbosity >= 2) {
            cout << "anchoring edge "  << i << " [" <<  a.pos1() << "," 
                << a.Pos1() << ") aligns to [" << a.pos2() << "," 
                << a.Pos2() << ") " << endl;
            cout << "path0= ";
            path0.Println(cout);
            if (trim_start == 0 || trim_end == 0) 
                cout << "trims= " << trim_start << "," << trim_end << endl;
            PrintVisualAlignment(True, cout, edge_seq, gplus_, a);
        }
        align_all_edges[e] = a;
        align_all_edges[e].Setpos2(a.pos2() - gplus_ext_); // related to G

        // find the edits
        vec<triple<int,int,String>> edits;
        vec<pair<String,String>> change;
        bool head_ins = (trim_start != 0);
        bool tail_ins = (trim_end != 0);
        GetEditsFromAlign(edge_seq, gplus_, a, &edits, &change, 
                edge_prev_base, 'X', head_ins, tail_ins);
        FilterAndModifyEdits(edits, change);
        vec<VariantCall> vcalls;
        for (size_t j = 0; j < edits.size(); j++) {
            vec<pair<int,int>> edge0_pos = FindEdgeHome(edits[j].first+trim_start, path0);
            for (size_t kk = 0; kk < edge0_pos.size(); kk++) {
                VariantCall vcall = { { gid, edits[j].second, edits[j].third, 
                                        change[j].first, change[j].second}, 
                                      { e, edits[j].first, 
                                        edge0_pos[kk].first, edge0_pos[kk].second} }; 
                vcalls.push_back(vcall);
            }
        }
        VariantCallGroup vcall_group(gid, a.pos2(), a.Pos2());
        a.Setpos2(a.pos2() - gplus_ext_); // related to G
        //new_hb_edge_mapping_.ConvertPathToOld(path0);
        vcall_group.AddBranch(path0, vcalls, a, make_pair(trim_start, trim_end),
                bubble_graph_edge_weight_[e]);
        vcall_group.AddRefWeight(vec<double>(bubble_graph_edge_weight_[e].size(), 0.0));
        all_groups[i*2] = vcall_group;
    }

    // removes some calls from one anchoring edge if neighboring edges overlap -- so that the same region wont get called twice
    for (size_t i = 0; i+1 < anchoring_edges.size(); i++) {
        auto& groupL = all_groups[2*i];
        auto& groupR = all_groups[2*i+2];
        if( groupL.GetEndPos() > groupR.GetStartPos()){
            int64_t l_last_v_pos = std::numeric_limits<int64_t>::min();
            for( const auto& b: groupL.GetVariantCalls() ){ for( const auto& v: b){ l_last_v_pos=std::max(int64_t(v.variant.pos),l_last_v_pos); } }
            int64_t r_first_v_pos = std::numeric_limits<int64_t>::max();
            for( const auto& b: groupR.GetVariantCalls() ){ for( const auto& v: b){ r_first_v_pos=std::min(int64_t(v.variant.pos),r_first_v_pos); } }
            if ( l_last_v_pos >= r_first_v_pos ){
                auto l_pad = abs(groupL.GetEndPos()-l_last_v_pos);
                auto r_pad = abs(groupR.GetEndPos()-r_first_v_pos);
                if( r_pad>l_pad ){
                    std::cout << "\nremoving calls from group " << 2*i << " located in its overlap with group " << 2*i+2 << "\n" <<std::endl;
                    for( auto& b: groupL.GetVariantCalls() ){
                        vec<Bool> to_delete(b.size(),False);
                        for(size_t ii=0;ii<b.size();++ii){ if( b[ii].variant.pos>=r_first_v_pos ){ to_delete[ii]=True; } }
                        EraseIf(b,to_delete);
                    }
                }
                else{
                    std::cout << "\nremoving calls from group " << 2*i+2 << " located in its overlap with group " << 2*i << "\n" <<std::endl;
                    for( auto& b: groupR.GetVariantCalls() ){
                        vec<Bool> to_delete(b.size(),False);
                        for(size_t ii=0;ii<b.size();++ii){ if( b[ii].variant.pos<=l_last_v_pos ){ to_delete[ii]=True; } }
                        EraseIf(b,to_delete);
                    }
                }
            }
        }
    }
    // find all paths between cell region
    int ncells = anchoring_edges.size() - 1;

    vec<int> to_right3; bubble_graph_.ToRight(to_right3);
    vec<int> to_left3;  bubble_graph_.ToLeft(to_left3);
    for ( int64_t i = 0; i < anchoring_edges.jsize()-1; i++ ) {
        int header_edge = anchoring_edges[i];
        int trailer_edge =  anchoring_edges[i+1];
        const basevector& header = bubble_graph_.EdgeObject(header_edge);
        int n1 = to_right3[header_edge];
        int n2 = to_left3[trailer_edge];
        int nbranches = bubble_graph_.From(n1).size();
        ForceAssertEq(bubble_graph_.From(n1).size(), bubble_graph_.To(n2).size());
        int cell_start = all_groups[2*i].GetEndPos();
        int cell_end = all_groups[2*i+2].GetStartPos();
        VariantCallGroup vcall_group(gid, cell_start, cell_end);
        if (cell_start >= cell_end) {
            cout << "   skip the overlapping region between "
                << header_edge << " and " << trailer_edge << endl;
            all_groups[2*i+1] = vcall_group;
            continue;
        }
        basevector cell_ref(gplus_, cell_start, cell_end - cell_start);

        vec<double> ref_weight;
        for (int j = 0; j < nbranches; j++) {
            int e = bubble_graph_.EdgeObjectIndexByIndexFrom(n1, j);
            vec<int> path0 =
                ConvertToAssemblyPathFromUnrolled(bubble_graph_edge_path_[e]);
            basevector branch_base = bubble_graph_.EdgeObject(e);

            // Align the edge to reference. penalize end-gaps during alignment
            int nerror; alignment al;

            /*
            int score = SmithWatAffine(branch_base, cell_ref, al, true, true);
            align a = align(al);
            */

            align a;
            int n1 = branch_base.size( ), n2 = cell_ref.size( );
            int offset = (-n1+n2)/2;
            int bandwidth = Max( offset + n1, n2 - offset );
            RestrictedAlign( branch_base, cell_ref, offset, bandwidth, a );
            vec<ho_interval> p1, p2;
            a.PerfectIntervals1( branch_base, cell_ref, p1 );
            a.PerfectIntervals2( branch_base, cell_ref, p2 );
            const int ML = 20;
            Bool left_ok = False, right_ok = False;
            for ( int i = 0; i < p1.isize( ); i++ )
            {    if ( p1[i].Start( ) == 0 && p2[i].Start( ) == 0 
                    && p1[i].Length( ) >= ML ) 
                 {    left_ok = True;    }
                 if ( p1[i].Stop( ) == branch_base.isize( )
                      && p2[i].Stop( ) == cell_ref.isize( )
                      && p1[i].Length( ) >= ML ) 
                 {    right_ok = True;    }    }
            if ( !left_ok || !right_ok )
            {    alignment al;
                 SmithWatAffine(branch_base, cell_ref, al, true, true);
                 a = align(al);    }

            if (verbosity >= 2) {
                cout << "align branch"  << j << " [" <<  a.pos1() << "," << a.Pos1() 
                    << ") aligns to [" << a.pos2() << "," << a.Pos2() << ") " << endl;
                cout << "path0= ";
                path0.Println(cout);
                PrintVisualAlignment(True, cout, branch_base, cell_ref, a);
            }
            a.Setpos2(a.pos2() + cell_start);
            align_all_edges[e] = a;
            align_all_edges[e].Setpos2(a.pos2() - gplus_ext_); // related to G

            vec<triple<int,int,String>> edits;
            vec<pair<String,String>> change;
            char prev_base = as_base(*(header.end() - bubble_graph_.K()));
            GetEditsFromAlign(branch_base, gplus_, a, &edits, &change, 
                    prev_base, 'X', true, true);
            FilterAndModifyEdits(edits, change);

            vec<VariantCall> vcalls;
            for (size_t k = 0; k < edits.size(); k++) {
                vec<pair<int,int>> edge0_pos = FindEdgeHome(edits[k].first, path0);
                for (size_t kk = 0; kk < edge0_pos.size(); kk++) {
                    VariantCall vcall = { { gid, edits[k].second, edits[k].third, 
                                            change[k].first, change[k].second }, 
                                          { e, edits[k].first, 
                                            edge0_pos[kk].first, edge0_pos[kk].second } }; 
                    vcalls.push_back(vcall);
                }
            }
            if (bubble_graph_edge_is_assembly_[e]) {
                a.Setpos2(a.pos2() - gplus_ext_); // related to G
                //new_hb_edge_mapping_.ConvertPathToOld(path0);
                vcall_group.AddBranch(path0, vcalls, a, make_pair(0, 0), bubble_graph_edge_weight_[e]);
            }
            if (bubble_graph_edge_is_ref_[e]) 
                ref_weight = bubble_graph_edge_weight_[e];
        }
        vcall_group.AddRefWeight(ref_weight);
        all_groups[2*i+1] = vcall_group;
    }
    copy(all_groups.begin(), all_groups.end(), back_inserter(*p_groups));
    copy(align_all_edges.begin(), align_all_edges.end(), back_inserter(*p_edge_aligns));

    if (verbosity >=1) cout << Date() << ": Ended calling variants" << endl;
};

// Determine the location on sub-edges given a position in the untrimmed path edge,
// The sub-edges are the edge id from hbp
vec<pair<int,int>> EdgesOnRef::FindEdgeHome(int pos_on_edge, const vec<int>& path0) const
{
    vec<pair<int,int>> edge0_start;
    int sub_start = 0, sub_end = 0;
    for (size_t i = 0; i < path0.size(); i++) {
        sub_start = (i == 0 ? 0 : sub_end - hb_.K() + 1);
        sub_end = sub_start + hb_.EdgeLengthBases(path0[i]);
        if (sub_start >= pos_on_edge) break;
        if (pos_on_edge < sub_end ) 
            edge0_start.push(path0[i], pos_on_edge - sub_start);
    }
    for (auto& p: edge0_start) {
        new_hb_edge_mapping_.ConvertEdgeIdAndOffsetToOld(p);
    }
    auto end = std::unique(edge0_start.begin(), edge0_start.end());
    edge0_start.resize(distance(edge0_start.begin(), end));
    return edge0_start;
}

void EdgesOnRef::DumpUnrolled(String filename, 
        const vec<pair<int,Bool>>* p_hbp_to_hb) const
{
    vec<String> edge_names2(dg_.EdgeObjectCount(),"");
    vec<double> lengths2( dg_.EdgeObjectCount( ) );
    for (size_t i = 0; i < edge_names2.size(); ++i) {
        int eid = GetOriginalId(i);
        if (p_hbp_to_hb != NULL)
            eid = (*p_hbp_to_hb)[eid].first;
        edge_names2[i] = ToString(i) + "->" + ToString(eid) + "@"
            + ToString(GetStart(i));
        lengths2[i] = hb_.EdgeLengthKmers(GetOriginalId(i));
    }
    ofstream dout(filename);
    dg_.PrettyDOT( dout, lengths2, digraphE<EdgeLoc>::edge_label_info(digraphE<EdgeLoc>
                ::edge_label_info::DIRECT, &edge_names2 ) );
}

// Linearize the graph to a string of bubbles connected by single-edge 
// anchoring edges. 
void EdgesOnRef::MakeBubbleGraph(int verbosity) 
{
    if (verbosity >=1)
        cout << Date() << ": Make bubble graph" << endl;
    vec<vec<int>> single_edges;
    vec<align> single_edges_aligns;
    FindAnchoringEdges(single_edges);
    AlignAnchoringEdges(single_edges, single_edges_aligns, verbosity);
    if (verbosity >=2)
        cout << "build the graph backbone" << endl;
    bubble_graph_.Clear();
    bubble_graph_.SetK(hb_.K());
    bubble_graph_.AddVertices(single_edges.size() * 2);
    for (size_t i = 0; i < single_edges.size(); i++) {
        vec<int> path0 = ConvertToAssemblyPathFromUnrolled(single_edges[i]);
        basevector edge = hb_.EdgePathToBases(path0);
        bubble_graph_.AddEdge(i*2, i*2+1, edge);
        bubble_graph_edge_path_.push_back(single_edges[i]);
        bubble_graph_edge_is_ref_.push_back(false);
        bubble_graph_edge_is_assembly_.push_back(true);
    }
    if (verbosity >=2)
        cout << "build the graph bubbles" << endl;
    for ( int64_t i = 0; i < single_edges.jsize()-1; i++ ) {
        int ref_start = single_edges_aligns[i].Pos2();
        int ref_stop = single_edges_aligns[i+1].pos2();

        int entrance_edge = single_edges[i].back();
        int exit_edge     = single_edges[i+1].front();
        const int NMaxPath = 50;
        vec<vec<int>> cell_paths;
        FindAllPathsNoLoop(dg_, entrance_edge, exit_edge, &cell_paths, NMaxPath);
        if (cell_paths.empty()) continue;
        if (ref_stop < ref_start) continue;

        basevector ref_seg;
        ref_seg.SetToSubOf(gplus_, ref_start, ref_stop-ref_start);
        vec<basevector> branches;
        for (size_t j = 0; j < cell_paths.size(); j++) {
            vec<int> path0 = ConvertToAssemblyPathFromUnrolled(cell_paths[j]);
            branches.push_back(hb_.EdgePathToBases(path0));
        }
        // remove mixed inversion edges from the graph, where structures like
        // X -> { Y, Y' } -> X' are found.
        const basevector& entrance_base = hb_.EdgeObject(GetOriginalId(entrance_edge));
        const basevector& exit_base = hb_.EdgeObject(GetOriginalId(exit_edge));
        if (entrance_base == ReverseComplement(exit_base)) {
            if (verbosity >=1) {
                cout << "Try to remove inversion between dg_ edge " <<
                    entrance_edge  << " and " << exit_edge << endl; 
            }
            RemoveInversion(cell_paths, branches);
        }

        bool has_ref = false;
        for (size_t j = 0; j < cell_paths.size(); j++) {
            const basevector& edge = branches[j];
            bubble_graph_.AddEdge(i*2+1, i*2+2, edge);
            bubble_graph_edge_path_.push_back(cell_paths[j]);
            bubble_graph_edge_is_ref_.push_back(edge == ref_seg);
            bubble_graph_edge_is_assembly_.push_back(true);
            if (edge == ref_seg) { has_ref = true; }
        }
        if (!has_ref) {
            bubble_graph_.AddEdge(i*2+1, i*2+2, ref_seg);
            bubble_graph_edge_path_.push_back(vec<int>());
            bubble_graph_edge_is_ref_.push_back(true);
            bubble_graph_edge_is_assembly_.push_back(false);
        }
    }
    ForceAssertEq(bubble_graph_edge_path_.isize(), bubble_graph_.EdgeObjectCount());
    bubble_graph_edge_weight_.assign(bubble_graph_edge_path_.size(), vec<double>());
}

void EdgesOnRef::DumpBubbleGraph(String filename) const
{
    vec<double> lengths(bubble_graph_.EdgeObjectCount());
    vec<String> edge_names(bubble_graph_.EdgeObjectCount());
    for (size_t e = 0; e < lengths.size(); e++) {
        lengths[e] = bubble_graph_.EdgeLengthKmers(e);
        edge_names[e] += ToString(e) + ": ";
        const vec<int>& path= bubble_graph_edge_path_[e];
        vec<int> path0 = ConvertToAssemblyPathFromUnrolled(path);
        for (size_t k = 0; k < path.size(); k++) {
            if (k != 0 ) edge_names[e] += ",";
            edge_names[e] += ToString(path0[k]);
        }
    }
    ofstream dout(filename);
    bubble_graph_.PrettyDOT(dout, lengths, HyperBasevector::edge_label_info(
                HyperBasevector::edge_label_info::DIRECT, &edge_names) );
}

// Calculate probability of being real for each paths.
void EdgesOnRef::PathProb(const vecbasevector& bases, const vecqualvector& quals, 
        int verbosity) 
{
    if (verbosity >= 1) 
        cout << Date() << ": Calculate probabilities of edges" << endl;
    
    // { { (rid11, diff11), (rid12, diff12), ...} // for edge1
    //   { (rid21, diff21), (rid22, diff22), ...} // for edge2
    //   ...  // }
    vec<vec<pair<int,int>>> homes_index(bubble_graph_.EdgeObjectCount());
//    FindReadHomesBest(bases, quals, bubble_graph_, &homes_index, vec<int>(), NULL, verbosity);
    vec<String> samples = p_read_tracker_->getSampleList();
    int nsamples = p_read_tracker_->getSampleList().size();
    if (verbosity >= 2) {
        for (size_t i = 0; i < samples.size(); i++) 
            cout << "sample " << i << " " << samples[i] << endl;
    }
    vec<vec<int>> supports(bubble_graph_.EdgeObjectCount(), vec<int>(nsamples,0));
    vec<vec<double>> probs(bubble_graph_.EdgeObjectCount(), vec<double>(nsamples,0.0));
    for (int n1 = 0; n1 < bubble_graph_.N(); n1++) {
        if (bubble_graph_.From(n1).size() < 1) continue;
        int n2 = bubble_graph_.From(n1)[0];
        vec<int> branches(bubble_graph_.From(n1).size());
        for (size_t j = 0; j < bubble_graph_.From(n1).size(); j++) 
            branches[j] = bubble_graph_.EdgeObjectIndexByIndexFrom(n1,j);
        const int qgood = 4;
        for (size_t j = 0; j < branches.size(); j++) {
            int e = branches[j];
            if (verbosity>=2) {
                cout << "branches " << j << " edge= " << e << endl;
                for (auto& x: homes_index[e]) {
                    int rid = x.first;
                    int sample_id = p_read_tracker_->getSampleID(rid);
                    if (x.second/1000 >= qgood && sample_id >= 0) 
                        cout << x.first << "," << x.second << " ";
                }
                cout << endl;
            }
            for (auto& x: homes_index[e]) {
                int rid = x.first;
                int sample_id = p_read_tracker_->getSampleID(rid);
                if (x.second/1000 >= qgood && sample_id >= 0) {
                    supports[e][sample_id] += x.second/1000;    
                    probs[e][sample_id] += x.second/1000;    
                }
            }
        }
    }
    bubble_graph_edge_weight_ = probs;
    bubble_graph_edge_supports_ = supports;
    if (verbosity >=1 ) 
        DumpBubbleGraph("bubble.dot");
}

// Find all paths between node n1 and n2 from graph dg_. Do not allow
// traversal of same edges twice

// auxilary data structure for FindAllPathsNoLoop
namespace {
struct PartialPath {
    vec<int> nodes;
    vec<int> edges;
};
}

template<typename GraphT>
void EdgesOnRef::FindAllPathsNoLoop(const GraphT& dg, int entrace_edge, int exit_edge, 
        vec<vec<int>>* allpaths, const int NMaxPath) const
{
    const int MaxDepth = 100;
    vec<int> to_right2; dg.ToRight(to_right2);
    vec<int> to_left2; dg.ToLeft(to_left2);
    int n1 = to_right2[entrace_edge];
    int n2 = to_left2[exit_edge];

    PartialPath start{{n1},vec<int>{}};
    stack<PartialPath> visited;
    visited.push(start);
    while (! visited.empty()) {
        if ((int)allpaths->size() > NMaxPath) break;
        PartialPath p0 = visited.top();
        visited.pop();
        int node = p0.nodes.back();
        if (node == n2) { allpaths->push_back(p0.edges); }
        if ((int)p0.edges.size() < MaxDepth) {
            for (int i = dg.From(node).size() -1; i >= 0; i--){
                int node_new = dg.From(node)[i];
                int edge_new = dg.EdgeObjectIndexByIndexFrom(node, i);
                if (!Member(p0.edges, edge_new)
                        && edge_new != exit_edge) {
                    PartialPath pnew = p0;
                    pnew.nodes.push_back(node_new);
                    pnew.edges.push_back(edge_new);
                    visited.push(pnew);
                }
            }
        }
    }
}

// We apply the following filters:
// 1. Delete edits that are outof range. 
// 2. make corrections to  position of edits as original alignment was with Gplus.
// 3. change the position of insertion and deletion to previous base.
// 4. Delete the large clumps of edits from inverted region or mis-alignments
void EdgesOnRef::FilterAndModifyEdits( vec<triple<int,int,String>>& edits, 
        vec<pair<String,String>>& change)
{
    ForceAssertEq(edits.size(), change.size());
    // only take edits in the genome_ region
    vec<Bool> todel(edits.size(), false);
    for (size_t j = 0; j < edits.size(); j++) {
        edits[j].second -= gplus_ext_;
        if (change[j].first.size() != change[j].second.size()) {
            edits[j].second -= 1;
            edits[j].first -= 1;
        }
        if ( edits[j].second < 0 || edits[j].second >= (int)genome_.size())
            todel[j] = True; 
    }
    // detect and remove large clumps of edits
    const int MinClumpSep = 30;
    const int MinEditsInClump = 30;
    const double MinFractionBaseChangeInClump = 0.5;
    for (size_t i = 0; i < edits.size(); i++) {
        int inserted_base = 0;
        int nmatch = 0;
        bool i_is_indel = (change[i].first.size() != change[i].second.size());
        if (i_is_indel) inserted_base += change[i].second.size()-1;
        size_t j = i + 1;
        while (j < edits.size() && abs(edits[j].second - edits[j-1].second 
                    - change[j-1].first.size()) < MinClumpSep) {
            nmatch += edits[j].second - edits[j-1].second - change[j-1].first.size();
            bool j_is_indel = (change[j].first.size() != change[j].second.size());
            if (j_is_indel) 
                inserted_base += change[j].second.size()-1;
            j++;
        }
        int total_bases = edits[j-1].second + change[j-1].first.size() 
            - edits[i].second + inserted_base;
        double fraction_base_change = 1 - (double)nmatch/ total_bases;
        if (int(j - i) > MinEditsInClump && 
                fraction_base_change > MinFractionBaseChangeInClump) {
            for (size_t k = i; k < j; k++)
                todel[k] = True;
        }
        /*
        cout << "i= " << i << " total_base= " << total_bases
            << " nmatch= " <<nmatch 
            << " nedits= " << j-i << endl;
        */
        i = j - 1;
    }
    EraseIf(edits, todel);
    EraseIf(change, todel);
}

vec<int> EdgesOnRef::ConvertToAssemblyPathFromUnrolled(const vec<int>& unrolled_path) const
{
    vec<int> path0;
    for (auto& x: unrolled_path) path0.push_back(GetOriginalId(x)); 
    return path0;
}

// Prepair backbone of the bubble graph. The backbone is single edged part in
// the graph.  Also remove any edge if the length is shorter than K after
// trimming the overlapping K-1 bases on both ends, or cannot be anchored in
// the genome (to be implemented). 
void EdgesOnRef::FindAnchoringEdges(vec<vec<int>>& single_edges) const
{
    single_edges.clear();
    vec<int> to_right2; dg_.ToRight(to_right2);
    vec<int> to_left2; dg_.ToLeft(to_left2);
    // backboen that are not covered by alternative path
    vec<int> backbone;
    for (size_t i = 0; i < backbone_.size(); i++) 
        if (backbone_dup_[i] == 0) backbone.push_back(backbone_[i]);
    // stitch continuous backbone edges
    for (size_t i = 0; i < backbone.size(); i++) {
        vec<int> path = {backbone[i]};
        size_t j = i + 1; // extension
        while (j < backbone.size()) {
            int node1 = to_right2[backbone[j-1]];
            int node2 = to_left2[backbone[j]];
            if (node1 != node2) {
                vec<vec<int>> allpaths;
                FindAllPathsNoLoop(dg_, backbone[j-1], backbone[j], &allpaths, 3);
                if (allpaths.size() != 1) break;
                else path.append(allpaths[0]);
            }
            path.push_back(backbone[j]);
            j++;
        }
        single_edges.push_back(path);
        i = j-1;
    }
}

// Remove short edges. Also remove edges that cannot be anchored in reference,
// which are probably due to large insertions.
void EdgesOnRef::AlignAnchoringEdges(vec<vec<int>>& single_edges,
        vec<align>& single_edge_aligns, int verbosity) 
{
    vec<int> to_right2; dg_.ToRight(to_right2);
    vec<int> to_left2; dg_.ToLeft(to_left2);
    single_edge_aligns.clear();
    single_edge_aligns.resize(single_edges.size());
    vec<Bool> todel(single_edges.size(), False);
    const int MinAlignedBase = 20;
    for (size_t i = 0; i < single_edges.size(); i++) {
        const vec<int>& path = single_edges[i];    // path in dg_
        vec<int> path0 = ConvertToAssemblyPathFromUnrolled(path);
        basevector edge = hb_.EdgePathToBases(path0);
        int trim_left = 0, trim_right = 0;
        if (dg_.To(to_left2[single_edges[i].front()]).size() != 0)
            trim_left = hb_.K() - 1;
        if (dg_.From(to_right2[single_edges[i].back()]).size() != 0)
            trim_right = hb_.K() - 1;
        if ((int)edge.size() - trim_left - trim_right < MinAlignedBase) {
            todel[i] = True;
            if (verbosity >= 2) {
                cout << "Erasing short anchoring edge " << i << " len=" 
                    << edge.size() << " trim_left= " << trim_left 
                    << " trim_right= " << trim_right << endl;
                cout << " path= ";
                path.Println(cout);
                cout << " path0= ";
                path0.Println(cout);
            }
        } else {
            // make the alignment
            basevector edge_trimmed(edge.begin()+trim_left, 
                    edge.end() - trim_right);
            align a;

            int bandwidth = max(1, (int)edge_trimmed.size() /2 );
            int path_offset = GetStart(path.front()) + trim_left;

            RestrictedAlign( edge_trimmed, gplus_,  path_offset, bandwidth, a );

            single_edge_aligns[i] = a;
            if (verbosity >=2) {
                cout << "anchoring edge " << i << " trim_left= " << trim_left
                    << " trim_right= " << trim_right  << " path0= ";
                path0.Println(cout);
                cout << "    found alignment " << " [" <<  a.pos1() << "," << a.Pos1() << ") to [" 
                    << a.pos2() << "," << a.Pos2() << ") " << endl;
            }
            if (a.Pos1() - a.pos1() < MinAlignedBase ||
                    a.Pos2() - a.pos2() < MinAlignedBase ){
                todel[i] = True;
                if (verbosity >=2) {
                    cout << "    Unable to align anchoring edge " << i << "to reference: "
                        << " alignment= " << " [" <<  a.pos1() << "," << a.Pos1() << ") to [" 
                        << a.pos2() << "," << a.Pos2() << ") " << endl;
                    cout << "   path0= ";
                    path0.Println(cout);
                    cout << "   path= ";
                    path.Println(cout);
                }
            } 
        }
    }
    if (todel.CountValue(True) > 0) {
        if (verbosity >= 2) {
            cout << "Erasing short anchoring edges:  " << endl;
            for (size_t i = 0; i < todel.size(); i++) {
                if (todel[i]) {
                    cout << " i= " << i << " path= ";
                    single_edges[i].Println(cout);
                }
            }
        }
        EraseIf(single_edges, todel);
        EraseIf(single_edge_aligns, todel);
    }
    ForceAssertEq(single_edges.size(), single_edge_aligns.size());
}

void EdgesOnRef::RemoveInversion(vec<vec<int>>& paths, vec<basevector>& edges) const 
{
    ForceAssertEq(paths.size(), edges.size());
    vec<int> backbone_counts(paths.size(), 0);
    for (size_t i = 0; i < paths.size(); i++) {
        int total = 0;
        for (auto& x: paths[i]) 
           if (Member(backbone_, x)) total++;
        backbone_counts[i] = total;
    }
    vec<Bool> todel(paths.size(), False);
    for (size_t i = 0; i < edges.size(); i++) {
        basevector rc_i = ReverseComplement(edges[i]);
        for (size_t j = i+1; j < edges.size(); j++) {
            if (todel[j]) continue;
            if (edges[j] == rc_i) 
                todel[backbone_counts[i] < backbone_counts[j] ? i:j] = True;
        }
    }
    if (todel.CountValue(True) > 0) {
        cout << "Delete " << todel.CountValue(True) << " inversion edges "  << endl;
        for (size_t i = 0; i < paths.size(); i++) {
            cout << "path= ";
            paths[i].Println(cout);
            cout << "path0= ";
            ConvertToAssemblyPathFromUnrolled(paths[i]).Println(cout);
        }
    } else {
        cout << "No inversion detected." << endl;
    }
    EraseIf(paths, todel);
    EraseIf(edges, todel);
}


// explicit instantiation of some digraph functions used locally
template int digraphE<EdgeLoc>::AddEdge( const int v, const int w, const EdgeLoc& e );
template void digraphE<EdgeLoc>::AddVertices(int nadd);
template void digraphE<EdgeLoc>::ToLeft(vec<int>& to_left) const;
template void digraphE<EdgeLoc>::ToRight(vec<int>& to_right) const;
