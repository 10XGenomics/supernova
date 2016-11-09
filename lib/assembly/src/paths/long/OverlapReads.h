///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef OVERLAP_READS_H
#define OVERLAP_READS_H

#include "CoreTools.h"
#include "Basevector.h"
#include "efasta/EfastaTools.h"
// The utilities to overlay a small number of efasta reads ( < 100 ) and create
// a linear assembly. 

class ReadOverlap{
public:
    int rid1, rid2;        // read ids
    int alt1, alt2;        // which alternatives from the two reads
    int overlap;           // how much is the overlap
    ReadOverlap(int r1, int r2, int a1, int a2, int o):
        rid1(r1),rid2(r2),alt1(a1),alt2(a2),overlap(o) {};
    ReadOverlap(){ rid1=-1; rid2=-1;};
    bool operator== ( const ReadOverlap& r ) const {
        return rid1 == r.rid1 && rid2 == r.rid2 &&
            alt1 == r.alt1 && alt2 == r.alt2 && overlap == r.overlap ;
    }
    friend ostream& operator<<( ostream& out, const ReadOverlap& r);
};

inline ostream& operator<< ( ostream& out, const ReadOverlap& r) {
    return out << r.rid1 << "_" << r.alt1 << "( " <<  r.overlap <<" )" << r.rid2 << "_" << r.alt2;
}

inline ostream& operator<<  ( ostream& out, const vec<ReadOverlap>& rs) {
    out << "[";
    for ( size_t i = 0; i < rs.size(); ++i ) out << rs[i] << ", ";
    out << "]";
    return out;
}


class ReadOverlapGraph {
  public:
    // Constructor and destructors

    ReadOverlapGraph( const vec< vec<basevector> > & reads ) : reads_(reads) {};

    // facilitating methods
    
    size_t NReads() const { return reads_.size(); }
    size_t Len(int rid) const { return reads_[rid][0].size(); }
    size_t Len(int rid, int alt) const { return reads_[rid][alt].size(); }
    int PathToLen( const vec<ReadOverlap>& path ) const;

    // Read overlap calculation
    
    // If tail b1[len1-overlap:len1) is the same as head b2[0: overlap)
    static bool BaseMatch( const basevector& b1,  const basevector& b2, int overlap );
    // return the overlap between the two reads e1 and e2, so that e1[len1 - overlap: len1)
    // perfectly overlap with e2[0, overlap). 
    // Todo: multiple overlap values
    static unsigned int Overlap( const basevector& r1, const basevector& r2, unsigned int overlap_lb, unsigned int overlap_ub );
    // Find if there are any overlaps from read 1 to read 2
    void FindReadOverlaps( int rid1, int rid2, int lb, int ub, vec<ReadOverlap> *er_overlaps );
    bool IsOverlaped( int rid1, int alt1, int rid2, int alt2 );

    // Graph generation and editing

    void Init( int Min_Overlap=100, const vec< pair<int,int> > *pt_starts=0 ) ;
    // If there are overlaps A->B, B->C, and A->C, remove the redundent A->C edge
    void RemoveTransitive();
    int RemoveTransitiveFrom(int rid, int alt, vec< vec<bool> >& visited );
    // Combine expanded reads if ambiguity sites do not overlap with other reads
    void CombineAmbReads( const vec< vec<Ambiguity> >& creads_ambs ) ;
    void GetAmbsOnPath( const vec<ReadOverlap>&  path, vec<Ambiguity>& new_ambs ) const;

    // Assembly methods

    void AllPathsFrom( int rid, int alt, vec< vec<ReadOverlap> >* p_paths, vec<int>* p_lens, int max_path ) const;
    void GetAssembly( vec<basevector>& assemblies, vec< vec<ReadOverlap> >& paths, int max_outputs = 10) const;
    void LinearAssembly( vec<basevector>& assemblies, vec< vec<ReadOverlap> >& paths) const;
    void PerfectAssembly( vec<basevector>& assembly, vec< vec<ReadOverlap> >& path) const;
    bool PerfectAssemblyBackward( int i, int j, basevector& assembly, vec< vec<int> >& read_locs, 
                                  vec<ReadOverlap>& path ) const;
    bool PerfectAssemblyForward( int i, int j, basevector& assembly, vec< vec<int> >& read_locs,
                                 vec<ReadOverlap>& path ) const;

    // debugging 

    void PrintGraph( ostream& out );
  private:
    vec< vec< vec<ReadOverlap> > > tos_;   // tos_[i] are a vector of overlaps from read i
    vec< vec< vec<ReadOverlap> > > froms_; // froms_[i] are a vector of overlaps to read i
    vec< vec<bool> > disable_;             // some reads are disabled because they are contained in other reads
    vec< vec<Ambiguity> > combined_ambs_;  // some reads expansions can be combined because the ambiguous sites don't overlap
    const vec< vec<basevector> >& reads_;
};


class SomeReadsAssembler {
private:
    vec< pair<int,int> > starts_;     // optional read starts
    vec< vec<Ambiguity> > ambs_;      // optional ambiguities 
    int Min_Overlap_ ;                // default min read overlap is 100
    int Min_Assembly_Len_ ;           // default min read overlap is 100
    int Max_Assembly_Num_;
    vec<basevector> assemblies_;      // possible assemblies
    vec< vec<ReadOverlap> > paths_;   // assembly paths
public:
    ReadOverlapGraph graph_;
    SomeReadsAssembler( const vec< vec<basevector> >& reads );
    // set assembly paramters
    void AddPosConstraints( const vec<pair<int,int> >& starts) { starts_ = starts; }
    void SetMinOverlap(int overlap) { Min_Overlap_ = overlap; }
    void SetMinAssemblyLen( int len ) { Min_Assembly_Len_ = len; }
    void SetMaxAssemblyNum( int n ) { Max_Assembly_Num_ = n; }
    void SetCombineAmbReads( const vec< vec<Ambiguity> >& ambs ) { ambs_ = ambs; }
    void InitGraph() {
        if ( starts_.empty() ) graph_.Init(Min_Overlap_);
        else graph_.Init(Min_Overlap_, &starts_);
        if ( ! ambs_.empty() ) graph_.CombineAmbReads( ambs_ );
    }
    // run assembly
    void Assembly() { graph_.GetAssembly(assemblies_, paths_, Max_Assembly_Num_); RemoveShortAssembly(); }
    // Generate assembly from the overlapping read, requires :
    // 1. No loop in the assembly
    // 2. The assembly contains more than 50% of the reads
    void LinearAssembly() { graph_.LinearAssembly(assemblies_, paths_); RemoveShortAssembly(); }
    // Return only if all the reads form a perfect assembly
    void PerfectAssembly() { graph_.PerfectAssembly(assemblies_,  paths_); RemoveShortAssembly(); }
    void RemoveShortAssembly();
    // Retrieving results
    const vec<basevector>& GetAssemblies() const { return assemblies_; }
    efasta GetEfastaAssembly(int i, size_t start, size_t stop) const;
    const vec<ReadOverlap>& GetAssemblyPath( int i ) const { return paths_[i]; }
    // get the sub-path the correspond to the region [start,stop) in the assembly
    vec<ReadOverlap> GetAssemblyPathRegion( int i, int start, int stop ) const;
};

#endif
