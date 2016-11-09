///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
#ifndef VARIANT_CALL_TOOLS_H
#define VARIANT_CALL_TOOLS_H

#include "CoreTools.h"
#include "paths/long/SupportedHyperBasevector.h"
#include <unordered_set>
#include <unordered_map>

// A class to re-align assembly to reference and call variants.
// 
// We start from the best path, which is essentially a linear graph. Then we
// add the other edges based on the connectivity in the assembly. At the end,
// every edge is assigned, and every path in the resulting graph is a possible
// chromosome.
//
// Then we find the anchoring edges, which are the single edges that can
// be unambigiously located on reference. We then enumerate the path between
// anchoring edges and make a bubble graph, which is a set of multiple-branch
// bubbles connected by the anchoring edges.
//
// The variants are then called using the bubble graph. The probability of the
// variants are calculated by realign reads back to the graph.
//

// forward declarations
class VariantCallGroup; 
class ReadOriginTracker;

struct EdgeLoc{
    int original_edge;
    int start;

    EdgeLoc(int original_edge, int start) : original_edge(original_edge), start(start) {}
    EdgeLoc(const pair<int,int>& edge_start) :original_edge(edge_start.first), start(edge_start.second) {}
};

// Track the original of the edge id and offset after long edge been broken
class EdgeTracker {
public:
   void AddMapping(int new_id, int old_id, int offset)  {
       edge_offset_mapping[new_id] = make_pair(old_id, offset);
   }
   void ConvertEdgeIdAndOffsetToOld(pair<int,int>& id_offset) const {
       auto it = edge_offset_mapping.find(id_offset.first);
       if (it != edge_offset_mapping.end()) {
           id_offset.first = it->second.first;
           id_offset.second += it->second.second;
       }
   }
   void ConvertPathToOld(vec<int>& path0) const {
       for (int& p0: path0) {
           auto it = edge_offset_mapping.find(p0);
           if (it != edge_offset_mapping.end()) p0 = it->second.first;
       }
   }
private:
    std::map<int,pair<int,int>> edge_offset_mapping;
};

class EdgesOnRef {
public:
    EdgesOnRef(const HyperBasevector& hb, const vec<pair<int,Bool>>& hbp_to_hb,
            const basevector& genome, const basevector& gplus, int gplus_ext,
            const ReadOriginTracker *p_read_tracker = NULL) 
        : hb_(hb), hbp_to_hb_(hbp_to_hb), genome_(genome), 
        gplus_(gplus), gplus_ext_(gplus_ext), p_read_tracker_(p_read_tracker)
    { genome_prev_char_ = (gplus_ext == 0 ? 'X': gplus[gplus_ext-1]); }

    // ============= Unrolled graph functions =============================
    //
    void InitFromBestPath(const vec<int>& path, const vec<pair<int,int>>& ref_pos) ;

    // Unroll the graph, find alignments of the unrolled edges, find variants,
    // and return the path of the edges.
    void UnrollAll(const int verbosity = 0, const int iLoBound=-1, const int iHiBound=std::numeric_limits<int>::max(),const bool bCorrection=false);

    // ============= Bubble graph functions =============================
    //
    // Dump the unrolled graph to a dot file.
    void DumpUnrolled(String filename, const vec<pair<int,Bool>>* p_hbp_to_hb = NULL ) const;

    // Linearize the graph to a string of bubbles connected by single-edge 
    // anchoring edges. 
    void MakeBubbleGraph(int verbosity = 0) ;

    // Dump the bubble graph to a dot file.
    void DumpBubbleGraph(String filename) const;

    // Calculate probability of being real for each paths.
    void PathProb(const vecbasevector& bases, const vecqualvector& quals, int verbosity = 0);

    // Call variants. Using bubble graph
    void CallVariantsGroupedWithProb(int gid, vec<VariantCallGroup> *p_groups, 
            vec<align>* p_edge_aligns, int verbosity = 0) ;

    // ============= Utility functions =============================
    
    // Find all paths between two edges from graph. Do not allow
    // traversal of same edges twice
    template<typename GraphT>
    void FindAllPathsNoLoop(const GraphT& dg, int entrace_edge, int exit_edge, 
            vec<vec<int>>* allpaths, const int NMaxPath) const;

    // Remove the invesion paths
    void RemoveInversion(vec<vec<int>>& paths, vec<basevector>& edges) const;

    // Delete edits that are outof range. Also make corrections to 
    // position of edits as original alignment was with Gplus.
    void FilterAndModifyEdits( vec<triple<int,int,String>>& edits, vec<pair<String,String>>& change);

    // Starting from an edge, exploring all path until it lands to the unrolled graph and
    // the inferred edge location is consistent with know value.
    bool FindAlternativePath(int node, vec<pair<int,int>> *p_edge_starts, 
            int *p_back_edge, const vec<vec<pair<int,int>>>& known_starts, 
            const int LenMax, const int MaxLocDevAllowed, int verbosity = 0);

    // Add edge in the unrolled graph given the alternative path info.
    void AddConnections(int front_edge, int back_edge, 
            const vec<pair<int,int>>& alternative_path_and_start);

    // Determine the location on sub-edges given a position in the untrimmed path edge,
    // Used for variant friend finding.
    vec<pair<int,int>> FindEdgeHome(int pos_on_edge, const vec<int>& path0) const;

    // Retern the original assembly edge id given the unrolled graph edge id;
    int GetOriginalId(int edge) const { return dg_.EdgeObject(edge).original_edge; }

    // Retern the ref location given the unrolled graph edge id;
    int GetStart(int edge) const { return dg_.EdgeObject(edge).start; }

    vec<int> ConvertToAssemblyPathFromUnrolled(const vec<int>& unrolled_path) const;

    // Prepair backbone of the bubble graph. The backbone is single edged part in
    // the graph.  Also remove any edge if the length is shorter than K after
    // trimming the overlapping K-1 bases on both ends, or cannot be anchored in
    // the genome (to be implemented). 
    void FindAnchoringEdges(vec<vec<int>>& single_edges) const;

    void AlignAnchoringEdges(vec<vec<int>>& single_edges,
            vec<align>& single_edge_aligns, int verbosity);

private:
    HyperBasevector                       hb_;
    const vec<pair<int,Bool>>           & hbp_to_hb_;
    const basevector                    & genome_;
    const basevector                    & gplus_;
    const int                             gplus_ext_;
    const ReadOriginTracker             * p_read_tracker_;
    // the previous base before the genome segment (needed for variant calling)
    // This is mostly useless now since we have gplus_ .
    char                        genome_prev_char_;  

    // unrolled graph
    digraphE<EdgeLoc>           dg_;
    vec<int>                    backbone_;
    vec<int>                    backbone_dup_;

    HyperBasevector             bubble_graph_;              //note that this hyperbasevector is not strictly a Kmer graph after fiticious reference edges, which can be shorter than K-1, have been added
    vec<vec<int>>               bubble_graph_edge_path_;    //edge_ids from dg_
    vec<vec<double>>            bubble_graph_edge_weight_;  //weight for each edge
    vec<vec<int>>               bubble_graph_edge_supports_;
    vec<bool>                   bubble_graph_edge_is_ref_;      // the edge is same as ref
    vec<bool>                   bubble_graph_edge_is_assembly_; // the edge was in assemlby

    EdgeTracker                 new_hb_edge_mapping_; //record the broken edges

    vec< std::tuple<int,int,int,int>> count_dg_dups()const;
    vec<std::pair<size_t,size_t>> count_dg_prefered_edge_connection(const std::unordered_set<int>& prefered_edges
                                                                   ,const vec<int>& to_left, const vec<int>& to_right)const;
    void clear_dg_branch(vec<int>& to_delete, int root_edge, int dir,const vec<int>& to_left,const vec<int>&to_right);

};



// Given the variants and callers, find the multiple placement of the caller
// edges and genome location (gid, pos) of the alternative placement.
void FindVariantFriends(const vec<VariantCallGroup>& vcall_groups, 
        const vec<vec<align>>& all_aligns, const HyperBasevector& hbp,
        const vec<pair<int,Bool>>& hbp_to_hb,
        map<Variant, vec<pair<int,int>>> *p_var_friending);

#endif
