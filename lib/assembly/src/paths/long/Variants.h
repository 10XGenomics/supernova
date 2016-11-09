///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef LONG_VARIANTS_H
#define LONG_VARIANTS_H

#include "CoreTools.h"
#include "efasta/EfastaTools.h"
#include "paths/HyperEfasta.h"
#include "paths/long/LongProtoTools.h"
#include "util/TextTable.h"
#include "paths/simulation/VCF.h"
// forward declarations
class SupportedHyperBasevector;
class HyperBasevector; 
class ReadOriginTracker;

// Store the signature sequence that can be used to find the variants in the assembly.
class VariantSignature {
public:
    // Only record the two branches at this moment
    VariantSignature(const basevector& branch1, const basevector& branch2)
        : head(), tail(),branches(2) {branches[0]=branch1;branches[1]=branch2;}
    const basevector& Branch1() const { return branches[0]; }
    const basevector& Branch2() const { return branches[1]; }
    const vec<basevector>& Branches() const { return branches; }
public:
    basevector head, /*branch1, branch2,*/ tail;
    vec<basevector> branches; // 0 - branch1, 1-branch2
};

// Mark the variants from the hyper efasta assembly. Return the actual number marked.
int MarkVariants( HyperEfasta& he, const vec<VariantSignature>& v_signatures,
     const long_logging& logc );

void VariantsMarkedDot(const String& head, const SupportedHyperBasevector& shb, 
       const vec<int>& varient_edge_list, const long_logging& logc );

void ReftraceVariants(ostream& out, ostream& callsOut, const vecbasevector& G, 
        const vecbasevector& Gplus, const vec<int>& Gplus_ext, 
        const HyperBasevector& hbp, const vec<pair<int,Bool>>& to_hb,
        const vec<vec<vec<pair<int,int>>>>& ref_pos_all_s, const vec<vec<vec<int>>>& best_path_s,
        const String& sVarFile="", const String& sRefFile="OUT_GENOME", 
        const RefTraceControl& ref_trace_control = RefTraceControl(),
        const ReadOriginTracker* p_read_tracker = NULL,
        const long_logging* logc = NULL) ;

// A variant describes a difference of the assembly to the genome sequence. Here
// the difference is recorded as an edit from the genome at given location.
// # All the locations are relative to the input genome, which maybe truncated.
// # The ref and alt start from pos-1, as specified by the output format.
struct Variant {
    int gid, pos;
    String edit;
    String ref, alt;
    friend bool operator< (const Variant& a, const Variant& b) {
        if (a.gid != b.gid) return a.gid < b.gid;
        if (a.pos != b.pos) return a.pos < b.pos;
        // sub is before ins/del if both at same position
        if (a.ref.size() == a.alt.size() && b.ref.size() != b.alt.size()) return true;
        if (a.ref.size() != a.alt.size() && b.ref.size() == b.alt.size()) return false;
        return a.edit < b.edit;
    }
    friend bool operator== (const Variant& a, const Variant& b) {
        return a.gid == b.gid && a.pos == b.pos && a.edit == b.edit;
    }
    friend ostream& operator<<(ostream& out, const Variant& var) {
        return out << var.gid << ":" << var.pos << var.ref << "->" << var.alt;
    }
};

// A caller is an assembly edge that is different from the genome at a
// certain position. Since the edge can be multiple edges from the assembly graph
// stitched together, we also record which sub-edge ID and position it belongs
// to.
struct Caller {
    int edge, pos;
    int sub_edge_id;            
    int sub_edge_pos;
    friend bool operator<(const Caller& a, const Caller& b){
        if (a.edge != b.edge) return a.edge < b.edge;
        if (a.pos != b.pos) return a.pos < b.pos;
        if (a.sub_edge_id != b.sub_edge_id) return a.sub_edge_id < b.sub_edge_id;
        return a.sub_edge_pos < b.sub_edge_pos;
    }
    friend ostream& operator<<(ostream& out, const Caller& var) {
        return out << var.pos << "@" << var.edge << "(" << var.sub_edge_pos 
            << "@" << var.sub_edge_id << ")";
    }
};

struct VariantCall {
    Variant variant; 
    Caller caller;  
    friend ostream& operator<< (ostream& out, const VariantCall& vc) {
        return out << vc.variant << " by " << vc.caller;
    }
};


// A VariantCallGroup contains a cell or a single branches in the unrolled graph.
// All branches in the group align to same region in the genome. 
class VariantCallGroup {
public:
    VariantCallGroup () :gid(-1), start(-1), end(-1) {};
    VariantCallGroup (int gid, int start, int end) :gid(gid), start(start), end(end){};

    void AddBranch(const vec<int>& path, const vec<VariantCall>& vcalls, 
            const align& a, const pair<int,int>& trim, vec<double> weight = vec<double>()) 
    { branches_.push_back(path); vcalls_.push_back(vcalls); 
      aligns_.push_back(a); weights_.push_back(weight); trims_.push_back(trim);}

    void AddRefWeight(vec<double> ref_weight) { ref_weight_ = ref_weight; }

    int GetGID()      const { return gid; }
    int GetStartPos() const { return start; }
    int GetEndPos()   const { return end; }

    // ================== branch accessor =========================
    int GetNBranch() const { return vcalls_.size(); }

    pair<int,int> GetTrim(int branch_id) const 
    { return trims_[branch_id];  }

    const vec<int>& GetBranchPath(int branch_id) const 
    { return branches_[branch_id];  }

    const align& GetAlign(int branch_id) const 
    { return aligns_[branch_id]; }

    const vec<VariantCall>& GetVariantCalls (int branch_id) const 
    { return vcalls_[branch_id]; }

    const vec<vec<VariantCall>>& GetVariantCalls () const 
    { return vcalls_; }

    vec<VariantCall>& GetVariantCalls (int branch_id) 
    { return vcalls_[branch_id]; }

    vec<vec<VariantCall>>& GetVariantCalls()
    { return vcalls_; }

    // ================== Output ==================================
    void PrintTable(TextTable& tb, ostream& callsOut, int& index, int& block_id,
            const vec<pair<int,Bool>>& hbp_to_hb, 
            const vec<String> ref_tags, const vec<size_t> ref_start,
            const map<Variant, vec<pair<int,int>>>& var_friending,
            size_t max_index, const map<Variant, vec<pair<double,double>>> probs
        ) const;

    void FillVCF(VCFWriter& vcf
                ,const vec<pair<int,Bool>>& hbp_to_hb
                ,const vec<String> ref_tags
                ,const vec<size_t> ref_start
                ,const map<Variant, vec<pair<int,int>>>& var_friending
                ,const map<Variant, vec<pair<double,double>>> probs
                ,const vec<String>&sample_list
                ,const size_t block_no) const;

private:
    int gid;                             // genome id
    int start, end;                      // covered region in reference genome
    vec<vec<int>> branches_;             // path (in hbp) of each branch
    vec<vec<VariantCall>> vcalls_;       // variant calls from each branch
    vec<align> aligns_;                  // alignment of the branch to reference
    vec<pair<int,int>> trims_;            // branch was trimmed 
    vec<vec<double>> weights_;           // weights for each branch
    vec<double> ref_weight_;             // weight for ref brach 
};

bool operator==(const VariantCallGroup&left, const VariantCallGroup&right);
bool operator<(const VariantCallGroup&left, const VariantCallGroup&right);

#endif
