///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef GAP_TOY_TOOLS_H
#define GAP_TOY_TOOLS_H

#include "Bitvector.h"
#include "CoreTools.h"
#include "Intvector.h"
#include "feudal/ObjectManager.h"
#include "feudal/PQVec.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "paths/long/SupportedHyperBasevector.h"
#include "paths/long/large/DiscoStats.h"
#include "paths/long/large/GapToyTools2.h"
#include "paths/long/large/GapToyTools3.h"
#include "paths/long/large/GapToyTools4.h"
#include "paths/long/large/GapToyTools5.h"
#include "paths/long/large/GapToyTools6.h"
#include "paths/long/large/Lines.h"
#include "10X/paths/ReadPathVecX.h"
#include "10X/Gap.h"
#include "system/SpinLockedData.h"
#include <memory>
#include <fstream>

class PerfStatLogger
{
public:
    static void init( String const& workDir )
    { delete gInst.mpOS;
      gInst.mpOS = new ofstream(workDir+"/statistics.txt"); }

    template <class T>
    static void log( String const& statName, T const& val, String const& gloss )
    { if ( gInst.mpOS )
       (*gInst.mpOS) << statName << '\t' << val << '\t' << gloss << std::endl; }

private:
    PerfStatLogger() : mpOS(nullptr) {}
    ~PerfStatLogger() { delete mpOS; }

    std::ofstream* mpOS;
    static PerfStatLogger gInst;
};

class GapToyResults {
     public:
     GapToyResults( ) : nedges(0) { }
     Bool Defined( ) const { return nedges > 0; }
     int nedges;
     double meanlen;
     int rev;
     int indels;
     int gaps;
     int subs;
};

void JoinPaths0( vec<int> const& inv, ReadPathVec& paths );

void JoinPaths(vec<int> const& inv, ReadPathVec& paths, HyperBasevector const& hbv,
     bool test = false );

//void JoinPathsAndHackReads(vec<int>& inv, String const& in_head, ReadPathVec& paths, HyperBasevector& hbv, bool test = false );

size_t FixPaths(HyperBasevector const& hbv, ReadPathVec& paths);

void DumpBPaths( vec<basevector> const& bp, int lroot,
          int rroot, String const& head );

void FetchPids( const String& FETCH_PIDS, const vec<int>& fosmids,
     const String& work_dir, const int mq );

String ToStringN( const vec<int>& x, const int vis = -1 );

void GapToyEvaluate( const String& SAMPLE, const String& species,
     const HyperBasevector& hb, const vecbasevector& G,
     const vec<int>& fosmids, const String& work_dir, const String& final_dir,
     int iAlignerK, const GapToyResults& res, Bool verbose );

void PrintDotMatchingGenome( const HyperBasevector& hb, const vecbasevector& G,
     vecString const& gNames, const String& work_dir );

void Dot( const int nblobs, int& nprocessed, int& dots_printed,
     const Bool ANNOUNCE, const int bl );

void FixInversion( const HyperBasevector& hb, vec<int>& inv2 );

void InsertPatch( HyperBasevector& hb, vec<int>& to_left, vec<int>& to_right,
     const HyperBasevector& hbp, const int lroot, const int rroot, const int left,
     const int right, vec<Bool>& used );

void DefinePairs( const ReadPathVec& paths, const vec<int>& inv,
     vec< pair<vec<int>,vec<int>> >& pairs, vec<int64_t>& pairs_pid,
     const String& dir );

void BasesToGraph( vecbasevector& bpathsx, const int K, HyperBasevector& hb );

void GetRoots( const HyperBasevector& hb, vec<int>& to_left, vec<int>& to_right,
     const vec<int>& lefts, const vec<int>& rights, int& lroot, int& rroot );

void MakeLocalAssembly1( const int lroot, const int rroot,
     const HyperBasevector& hb, const vecbasevector& bases,
     const VecPQVec& quals, const vec<int64_t>& pids, const String& TMP,
     ostringstream& mout, const Bool LOCAL_LAYOUT, const int K2_FLOOR,
     const String& work_dir, VecEFasta& corrected, vecbasevector& creads,
     vec<pairing_info>& cpartner, vec<int>& cid, LongProtoTmpDirManager& tmp_mgr );

void PlaceMore( const HyperBasevector& hb, const vecbasevector& bases,
     const VecPQVec& quals, ReadPathVec& paths2, vec<int64_t>& placed,
     const int place_more_level, const Bool verbose = False );

void ExtraPaths( const HyperBasevector& hb, const vecbasevector& bases,
     const VecPQVec& quals, ReadPathVec& paths2 );

void AnalyzeBranches( HyperBasevector& hb, vec<int>& to_right, const vec<int>& inv2,
     ReadPathVec& paths2, const Bool ANALYZE_BRANCHES_REV,
     const int min_ratio2, const Bool ANALYZE_BRANCHES_VERBOSE );

void ExtendTerminalEdges( const HyperBasevector& hb,
     const vec< vec<int> >& layout_pos, const vec< vec<int64_t> >& layout_id,
     const vec< vec<Bool> >& layout_or, const vecbasevector& bases,
     const VecPQVec& quals );

void SelectSpecials( const HyperBasevector& hb, vecbasevector& bases,
     VecPQVec const& quals, const ReadPathVec& paths2, const String& work_dir );

void LayoutReads( const HyperBasevector& hb, const vec<int>& inv,
     const vecbasevector& bases, const ReadPathVec& paths,
     vec<vec<int>>& layout_pos, vec<vec<int64_t>>& layout_id,
     vec<vec<Bool>>& layout_or );

void SortBlobs( const HyperBasevector& hb,
     const vec< triple< pair<int,int>, triple<int,vec<int>,vec<int>>, vec<int> > >&
          blobber,
     vec< pair<int,int> >& blobs );

void RemoveHangs( HyperBasevector& hb, vec<int>& inv, ReadPathVec& paths,
     const int max_del );

void Degloop( const int mode, HyperBasevector& hb, vec<int>& inv, ReadPathVec& paths,
     const vecbasevector& bases, const VecPQVec& quals, const double min_dist,
     const int verbosity = 0 );

// DegloopCore: H = HyperBasevector or HyperBasevectorX.

template<class H> void DegloopCore( const int mode, H& hb, vec<int>& inv,
     ReadPathVec& paths, const vecbasevector& bases, const VecPQVec& quals,
     const VecULongVec& paths_index, const int v, const int pass,
     const double min_dist, vec<int>& EDELS, const int verbosity,
     const vec<int>* ids = NULL );

void Patch( HyperBasevector& hb, const vec< pair<int,int> >& blobs,
     vec<HyperBasevector>& mhbp, const String& work_dir, const vec<String>& mreport,
     vecbvec& new_stuff );

void CleanupCore( HyperBasevector& hb, vec<int>& inv, ReadPathVec& paths );

void Cleanup( HyperBasevector& hb, vec<int>& inv, ReadPathVec& paths );

void CleanupLoops( HyperBasevector& hb, vec<int>& inv, ReadPathVec& paths );

void RemoveUnneededVertices( HyperBasevector& hb, vec<int>& inv,
     ReadPathVec& paths );

// For two-edge loops only:

void RemoveUnneededVerticesLoopsOnly( HyperBasevector& hb, vec<int>& inv,
     ReadPathVec& paths );

void RemoveUnneededVertices2( HyperBasevector& hb, vec<int>& inv,
     ReadPathVec& paths, Bool debug = false, const Bool single = False );

void RemoveSmallComponents3( HyperBasevector& hb,
     const Bool remove_small_cycles = False );

void Empty( const HyperBasevector& hb, const vec< pair<vec<int>,vec<int>> >& pairs,
     const vec<int64_t>& pairs_pid, vec<vec<int>>& left_empty,
     vec<vec<int>>& right_empty, const Bool EMPTY2 );

void AddNewStuff( vecbvec& new_stuff, HyperBasevector& hb, vec<int>& inv2,
     ReadPathVec& paths2, const vecbasevector& bases, const VecPQVec& quals,
     const int MIN_GAIN, const vec<int>& TRACE_PATHS, const String& work_dir,
     const int EXT_MODE );

void ExtendPath( ReadPath& p, const int64_t i, const HyperBasevector& hb,
     const vec<int>& to_right, const bvec& bases,
     const qvec& quals, const int min_gain, const Bool verbose, const int mode );

void ExtendPath( ReadPath& p, const int64_t i, const HyperBasevectorX& hb,
     const bvec& bases,
     const qvec& quals, const int min_gain, const Bool verbose, const int mode );

void ExtendPath2( ReadPath& p, const int64_t i, const HyperBasevector& hb,
     const vec<int>& to_left, const vec<int>& to_right, const bvec& bases,
     const qvec& quals, const int min_gain, const Bool verbose,
     const int mode );

// a class for dealing with bubbles and support
// the basic usage is
// 1) bubble_logger(graph,inv)
// 2) for each read, call log_read(read,qual,read-path)
// 3) call getData() to analyze data
class bubble_logger{
public:
    //stores the state of the bubble
    struct bubble_data_t{
        struct support_t{
            support_t(int a, int b):read_branch(a),branch_branch(b){}
            int read_branch;   // qsum between read and path
            int branch_branch; // qsum difference against the losing branch
        };
        bubble_data_t(const int branch0_edge,const int branch1_edge)
            :branch_edges({branch0_edge,branch1_edge})
            ,branch_supports(2)
            ,lock_ptr(new SpinLockedData)
            {}
        bubble_data_t(const int branch0_edge,const int branch1_edge
                   ,const int branch0_edge_rc, const int branch1_edge_rc)
            :branch_edges({branch0_edge,branch1_edge,branch0_edge_rc,branch1_edge_rc})
            ,branch_supports(4)
            ,lock_ptr(new SpinLockedData)
            {}
        bubble_data_t(bubble_data_t && other)noexcept
            :branch_edges(std::move(other.branch_edges))
            ,branch_supports(std::move(other.branch_supports))
            ,lock_ptr(std::move(other.lock_ptr)) { }

        void addSupport(size_t branch, support_t const&weight);
        vec<support_t> getSupport(size_t branch)const{
            SpinLocker locker(*lock_ptr);
            return branch_supports[branch];
        };
        vec<int>const& getEdges()const{ return branch_edges;};

    private:
        bubble_data_t();
        bubble_data_t(bubble_data_t const&other);
        const vec< int > branch_edges;
        vec< vec<support_t> > branch_supports;
        std::unique_ptr<SpinLockedData> lock_ptr;
    };

    // given the problem statement this is the only meaningful constructor
    bubble_logger(const HyperBasevector& hb, const vec<int>& inv);

    // if one or more edges in rp is part of a bubble, perform gap-free alignment on the path as the alternate path
    // and collect the result
    // returns true if at any point (qsum of alt path) < (qsum of orig path)
    bool log_read(basevector const&read, qualvector const&qual, ReadPath const&rp, bool bVerbose=false);

    //do a gap-free alignment of the read against the graph according to rp
    //returns the sum of read-quality-score at the position the sequence mismatches
    int getQ(basevector const&read, qualvector const&qual, ReadPath const&rp, const qualvector::value_type min_q=4);

    // assign weight to the bubble-branch of an edge
    void addWeight(int edge,bubble_data_t::support_t const&weight){
        const auto& b_b=edge_bubble_branch_[edge];
        int bubble_data_size = bubble_data_.size();
        ForceAssertLt(b_b.first, bubble_data_size);
        bubble_data_[ b_b.first ].addSupport(b_b.second,weight);
    }

    // if edge is part of a bubble, return the edge index corresponding to another branch, or -1 otherwise
    int alt(int edge)const{ return edge_alt_[edge];}

    // if an edge is in a bubble
    bool inBubble(int edge)const{ return edge_bubble_branch_[edge].first>=0;}

    // return the current list of bubble data
    const std::vector<bubble_data_t>& getData()const{return bubble_data_;};

private:
    bubble_logger();
    const HyperBasevector& hb_;
    vec<int> edge_alt_;                            // [edge_id] -> edge id of another branch
    vec< std::pair<int,int> > edge_bubble_branch_; // [edge_id] -> (bubble,branch)
    std::vector< bubble_data_t > bubble_data_;             // [bubble_id] -> data

};
std::ostream& operator<<(std::ostream&os, bubble_logger const& in);
std::ostream& operator<<(std::ostream& os, bubble_logger::bubble_data_t const&in);
void PopBubbles( HyperBasevector& hb , const vec<int>& inv2
               , const vecbasevector& bases, const VecPQVec& quals, const ReadPathVec& paths2);
void PrintBubbles( std::ostream& os, HyperBasevector& hb , const vec<int>& inv2
               , const vecbasevector & bases, const VecPQVec& quals, const ReadPathVec& paths2);

// HIGHLY INCOMPLETE:

void Validate( const HyperBasevector& hb, const ReadPathVec& paths );
void Validate( const HyperBasevector& hb, const ReadPathVecX& paths );
void Validate( const HyperBasevectorX& hb, const ReadPathVecX& paths );
void Validate( const HyperBasevector& hb, const HyperBasevectorX& hbx, const ReadPathVecX& paths );

void TestIndex( const HyperBasevector& hb,
        const ReadPathVec& paths, const VecULongVec& invPaths);

void TestInvolution( const HyperBasevector& hb, const vec<int>& inv );

void DeleteFunkyPathPairs( const HyperBasevector& hb, const vec<int>& inv,
     const vecbasevector& bases, ReadPathVec& paths, const Bool verbose );

void AlignToGenome( const HyperBasevector& hb, const vec<int>& inv,
     const vecbasevector& genome, vec< vec< pair<int,int> > >& hits,
     const int K2 = 500 );

// A perf_place is a perfect match between an assembly (extending across one or
// more edges) and a genome reference sequence.  The starting position on the
// genome is (g,start).  The sequence of assembly edges is e.  The start on the
// first edge is estart.  The total length of the match in bases is len.

class perf_place {

     public:

     int G( ) const { return g; }
     int Len( ) const { return len; }
     int Gstart( ) const { return gstart; }
     int Gstop( ) const { return gstart + len; }
     const vec<int>& E( ) const { return e; }

     int Estart( ) const { return estart; }

     int Estop( const HyperBasevector& hb ) const
     {    int x = estart + len;
          for ( int j = 0; j < e.isize( ) - 1; j++ )
               x -= hb.Kmers( e[j] );
          return x;    }

     int g;
     int gstart;
     int len;
     vec<int> e;
     int estart;

};

// AlignToGenomePerf: find perfect fw matches of length at least 200 between the
// assembly and the genome.
//
// perfs: { (g,gstart), (e,estart), len ) }

void AlignToGenomePerf( const HyperBasevector& hb, const vecbasevector& genome,
     vec< triple< pair<int,int>, pair<int,int>, int > >& perfs,
     vec<perf_place>& places );

void ReroutePaths( const HyperBasevector& hb, const vec<int>& inv,
     ReadPathVec& paths, const vecbasevector& bases, const VecPQVec& quals );

void Tamp( HyperBasevector& hb, vec<int>& inv, ReadPathVec& paths,
     const int max_shift );

void ExtendPairs60( HyperBasevector& hb, vec<int>& inv, ReadPathVec& paths );

void PlacePartners( const HyperBasevector& hb, const vec<int>& inv,
     ReadPathVec& paths, const vecbasevector& bases, const VecPQVec& quals );

void ReduceK( const int newK, HyperBasevector& hb, const vec<int>& inv,
     ReadPathVec& paths );

void LogTime( const double clock, const String& what, const String& work_dir = "" );

void AssayMisassemblies( const HyperBasevector& hbx, const vec<int>& inv,
     const vec< vec< pair<int,int> > >& hits,
     const vec<vec<vec<vec<int>>>>& linesx, const String& final_dir );

void Clean200( HyperBasevector& hb, vec<int>& inv, ReadPathVec& paths,
     const vecbasevector& bases, const VecPQVec& quals,
     const int verbosity = 0 // 0 or 1
     );

void MakeFinalFasta( const HyperBasevector& hbx, const vec<int>& inv2,
     const vec<vec<vec<vec<int>>>>& linesx, const vec<int>& npairsx,
     const vec<vec<covcount>>& covs, vec< vec< pair<int,int> > > hits,
     const String& final_dir, const String& work_dir, const Bool ALIGN_TO_GENOME );

String Chr( const int g );

void PartnersToEnds( const HyperBasevector& hb, ReadPathVec& paths,
                        const vecbasevector& bases, const VecPQVec& quals );
void PartnersToEndsOld( const HyperBasevector& hb, ReadPathVec& paths,
                        const vecbasevector& bases, const VecPQVec& quals );

void BuildGenomeMap( const HyperBasevector& hb, const vec<int>& inv,
     const vec<vec<vec<vec<int>>>>& lines, const vec<vec<covcount>>& covs,
     const vec< vec< pair<int,int> > >& hits, const vecbitvector& genome_amb,
     const vec<String>& genome_names, const String& final_dir );

void FragDist( const HyperBasevector& hb, const vec<int>& inv,
     const ReadPathVec& paths, const String out_file );

void UnwindThreeEdgePlasmids(HyperBasevector& hb, vec<int>& inv, ReadPathVec& paths);

// Compute fraction of edges in the assembly whose copy number is within 10% of
// the nearest integer value.

double CNIntegerFraction(const HyperBasevector& hb, const vec<vec<covcount>>& covs,
			 const double frac = 0.1, const int min_edge_size = 2000);

String PrintHits( const int e, vec<vec<pair<int,int>>>& hits,
     const HyperBasevectorX& hb, const vec<int>& inv,
     const vec<String>& genome_names );

void RemoveUnneededVerticesGeneralizedLoops( HyperBasevector& hb, vec<int>& inv,
     ReadPathVec& paths );

void BuildSeqGaps( vecbasevector& allx, const HyperBasevector& hb, 
        const digraphE<vec<int>>& D );

void BuildAll( vecbasevector& all, const HyperBasevector& hb,
     const int64_t extra = 0 );

void TranslatePaths( ReadPathVec & paths2, const HyperBasevector & hb3,
     const vec<vec<int>>& to3, const vec<int>& left3 );

void TranslatePaths( ReadPathVecX & paths2, const HyperBasevectorX & hbx3,
     const vec<vec<int>>& to3, const vec<int>& left3, const String paths_file );

int N50PerfectStretch( const HyperBasevector& hb, const vecbasevector& genome,
     const Bool concatenate );

void ReportMemory( const disco_stats& stats );

void PrintSysInfo( );

void MemoryCheck( const Bool MEMORY_CHECK, const String& work_dir );

void DefineRegions( const String& X, vec<int>& fosmids, vec<String>& regions,
     map<String,GapToyResults>& res, const int PAD, Bool& all,
     const String& SAMPLE, String& EVALUATE, const String& F );

#endif
