///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// A SupportedHyperBasevector is a HyperBasevector, together with paths through it,
// represented as sequences of edges, and for each path, a weight.  It also includes
// an 'involution', by which we mean a map i: {edges} --> {edges} such that 
// i^2 = 1, and which maps an edge to an edge that is the reverse complement of it, 
// however if an entire component has already had its reverse complement deleted, we 
// map all its edges to -1.  In addition, we include the median corrected read 
// length, and might add other data of this sort.  The paths are assumed to be 
// unique-ordered, however some methods do the unique ordering as a last, separate 
// step.
//
// A path is represented as a sequence of integers x.  The special value -1 is used 
// to represent a temporary hole that has been introduced in the path, typically 
//
// A SupportedHyperBasevector also has a set of pairs.  These are pairs of paths,
// together with {(weight,trim,lib)} each representing a weight, a trim amount, and
// a library id.

// NOTE: The .cc files for this class are split.

#ifndef SUPPORTED_HYPERBASEVECTOR_H
#define SUPPORTED_HYPERBASEVECTOR_H

#include "Basevector.h"
#include "CoreTools.h"
#include "Qualvector.h"
#include "math/IntDistribution.h"
#include "paths/HyperBasevector.h"
#include "paths/long/Fix64_6.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/SupportedHyperBasevector2.h"
#include "paths/long/SupportedHyperBasevector3.h"
#include "paths/long/SupportedHyperBasevector4.h"
#include "paths/long/SupportedHyperBasevector5.h"
#include "paths/long/SupportedHyperBasevector6.h"
#include "paths/long/SupportedHyperBasevector7.h"
#include "paths/long/SupportedHyperBasevector8.h"
#include "paths/long/PairInfo.h"
#include "paths/long/Variants.h"
#include "paths/long/RefTrace.h"

class SupportedHyperBasevectorPrivate : public HyperBasevector {

     public:

     // **** CONSTRUCT AND TEST METHODS *****

     SupportedHyperBasevectorPrivate( ) { }

     SupportedHyperBasevectorPrivate( const HyperBasevector& hb, 
          const vec<int>& inv, const vec< vec<int> >& paths, 
          const vec<fix64_6>& weights_fw, 
          const vec<fix64_6>& weights_rc, 
          const vec< vec< pair<fix64_6,int64_t> > >& weights_fw_origin,
          const vec< vec< pair<fix64_6,int64_t> > >& weights_rc_origin,
          const vec< pair< vec<int>, vec<int> > >& pairs,
          const vec< vec<pair_point> >& pair_data,
          const int64_t nreads, const IntDistribution& read_length_dist,
          const double fudge_mult )
          : HyperBasevector(hb), inv_(inv), paths_(paths), weights_fw_(weights_fw),
          weights_rc_(weights_rc), weights_fw_origin_(weights_fw_origin),
          weights_rc_origin_(weights_rc_origin), pairs_(pairs), 
          pair_data_(pair_data), nreads_(nreads), 
          read_length_dist_(read_length_dist), fudge_mult_(fudge_mult) { }

     void TestValid( const long_logging& logc, const Bool test_paths = True ) const;

     // **** ACCESSORS FOR PRIVATE DATA *****

     int NPaths( ) const { return paths_.size( ); }
     const vec<int>& Path( const int i ) const { return paths_[i]; }
     vec<int>& PathMutable( const int i ) { return paths_[i]; }
     int Path( const int i, const int j ) const { return paths_[i][j]; }
     int& PathMutable( const int i, const int j ) { return paths_[i][j]; }
     const vec< vec<int> >& Paths( ) const { return paths_; }
     vec< vec<int> >& PathsMutable( ) { return paths_; }

     fix64_6 WeightFw( int i ) const { return weights_fw_[i]; }
     const vec< pair<fix64_6,int64_t> >& WeightFwOrigin( int i ) const
     { return weights_fw_origin_[i]; }
     vec< pair<fix64_6,int64_t> >& WeightFwOriginMutable( int i )
     { return weights_fw_origin_[i]; }
     fix64_6& WeightFwMutable( const int i ) { return weights_fw_[i]; }
     const vec<fix64_6>& WeightsFw( ) const { return weights_fw_; }
     vec<fix64_6>& WeightsFwMutable( ) { return weights_fw_; }
     fix64_6 WeightRc( int i ) const { return weights_rc_[i]; }
     const vec< pair<fix64_6,int64_t> >& WeightRcOrigin( int i ) const
     { return weights_rc_origin_[i]; }
     vec< pair<fix64_6,int64_t> >& WeightRcOriginMutable( int i )
     { return weights_rc_origin_[i]; }
     fix64_6& WeightRcMutable( const int i ) { return weights_rc_[i]; }
     const vec<fix64_6>& WeightsRc( ) const { return weights_rc_; }
     vec<fix64_6>& WeightsRcMutable( ) { return weights_rc_; }
     fix64_6 Weight( int i ) const { return weights_fw_[i] + weights_rc_[i]; }
     const vec< vec< pair<fix64_6,int64_t> > >& WeightsFwOrigin( ) const
     { return weights_fw_origin_; }
     const vec< vec< pair<fix64_6,int64_t> > >& WeightsRcOrigin( ) const
     { return weights_rc_origin_; }
     vec< vec< pair<fix64_6,int64_t> > >& WeightsFwOriginMutable( )
     { return weights_fw_origin_; }
     vec< vec< pair<fix64_6,int64_t> > >& WeightsRcOriginMutable( )
     { return weights_rc_origin_; }

     int NPairs( ) const { return pairs_.size( ); }
     const pair< vec<int>, vec<int> >& 
          Pair( const int i ) const { return pairs_[i]; }
     pair< vec<int>, vec<int> >& PairMutable( const int i ) { return pairs_[i]; }
     const vec<int>& PairLeft( int i ) const { return pairs_[i].first; }
     vec<int>& PairLeftMutable( int i ) { return pairs_[i].first; }
     const vec<int>& PairRight( int i ) const { return pairs_[i].second; }
     vec<int>& PairRightMutable( int i ) { return pairs_[i].second; }
     int PairLeft( const int i, const int j ) const { return pairs_[i].first[j]; }
     int& PairLeftMutable( const int i, const int j ) { return pairs_[i].first[j]; }
     int PairRight( const int i, const int j ) const { return pairs_[i].second[j]; }
     int& PairRightMutable( const int i, const int j ) 
          { return pairs_[i].second[j]; }
     const vec< pair< vec<int>, vec<int> > >& Pairs( ) const { return pairs_; }
     vec< pair< vec<int>, vec<int> > >& PairsMutable( ) { return pairs_; }

     const vec< vec<pair_point> >& PairData( ) const { return pair_data_; }
     vec< vec<pair_point> >& PairDataMutable( ) { return pair_data_; }
     const vec<pair_point>& PairData( const int i ) const { return pair_data_[i]; }
     vec<pair_point>& PairDataMutable( const int i ) { return pair_data_[i]; }
     const pair_point& PairData( const int i, const int j ) const
     {    return pair_data_[i][j];    }
     pair_point& PairDataMutable( const int i, const int j )
     {    return pair_data_[i][j];    }

     int Inv( const int x ) const { return inv_[x]; }
     int& InvMutable( const int x ) { return inv_[x]; }
     const vec<int>& Inv( ) const { return inv_; }
     vec<int>& InvMutable( ) { return inv_; }
     Bool InvDef( const int x ) const { return inv_[x] != -1; }

     int64_t ReadCount( ) const { return nreads_; }
     int64_t& ReadCountMutable( ) { return nreads_; }
     int MedianCorrectedReadLength( ) const { return read_length_dist_.median( ); }
     int MedianCorrectedReadLengthFudge( ) const 
     { return int( round( read_length_dist_.median( ) * 1.0 ) ); }
     const IntDistribution& ReadLengthDist( ) const { return read_length_dist_; }
     IntDistribution& ReadLengthDistMutable( ) { return read_length_dist_; }
     double FudgeMult( ) const { return fudge_mult_; }
     double& FudgeMultMutable( ) { return fudge_mult_; }

     // **** PRIVATE MEMBERS *****

     private:

     vec<int> inv_;
     vec< vec<int> > paths_;
     vec<fix64_6> weights_fw_, weights_rc_;
     vec< vec< pair<fix64_6,int64_t> > > weights_fw_origin_;
     vec< vec< pair<fix64_6,int64_t> > > weights_rc_origin_;
     vec< pair< vec<int>, vec<int> > > pairs_;
     vec< vec<pair_point> > pair_data_;
     int64_t nreads_;
     IntDistribution read_length_dist_;
     double fudge_mult_;

};

class SupportedHyperBasevector : public SupportedHyperBasevectorPrivate {

     public:

     // **** BASIC UTILITIES *****

     SupportedHyperBasevector( ) : SupportedHyperBasevectorPrivate( ) { }

     SupportedHyperBasevector( const HyperBasevector& hb, 
          const vec<int>& inv, const vec< vec<int> >& paths, 
          const vec<fix64_6>& weights_fw, const vec<fix64_6>& weights_rc, 
          const vec< vec< pair<fix64_6,int64_t> > >& weights_fw_origin,
          const vec< vec< pair<fix64_6,int64_t> > >& weights_rc_origin,
          const vec< pair< vec<int>, vec<int> > >& pairs,
          const vec< vec<pair_point> >& pair_data,
          const int64_t nreads, const IntDistribution& read_length_dist,
          const double fudge_mult )
          : SupportedHyperBasevectorPrivate( hb, inv, paths, weights_fw, weights_rc,
          weights_fw_origin, weights_rc_origin, pairs, pair_data, nreads, 
          read_length_dist, fudge_mult ) { }

     void Reverse( );

     void AddTrim( const int id, const int t )
     {    for ( int l = 0; l < PairDataMutable(id).isize( ); l++ )
               PairDataMutable(id,l).TrimMutable( ) += t;    }

     // ***** COMPLEX MODIFIERS *****

     // PullApart: consider a diagram in which there is an edge r, with exactly two
     // edges x1, x2 abutting r on the left, and exactly two edges y1, y2 abutting r
     // on the right.  Suppose both x1,r,y1 and x2,r,y2 occur in the paths with
     // weight >= min_weight_split, and that x1,r,y2 and x2,r,y1 occur not at all.
     // Then replace all five edges by two new edges z1, z2.  Paths beginning or
     // ending in any of the five edges are translated appropriately.

     void PullApart( const double min_weight_split, const long_logging& logc );

     // PullApart2: a bit more general than PullApart in some ways, but only
     // implemented for inversion-free components.
     // Defined in SupportedHyperBasevector6.cc.
     void PullApart2( const double min_weight_split, const long_logging& logc );

     // UnwindAssembly.  Attempt to simplify the assembly.
     // 
     // WARNING: This is not a bona fide SupportedHyperBasevector method, because
     // it does not update the paths and weights.  Call BootStrap to fix.

     void UnwindAssembly( const long_logging_control& log_control,
          const long_logging& logc );

     // Delete very weak edges.

     void DeleteWeakEdges( const long_logging& logc );

     // ChunkEdges.  Find sequences of edges x1,...,xn satisfying the following:
     // (1) n >= 2.
     // (2) The length of x2,...,xn-1 in kmers is at most the median corrected
     //     read length.
     // (3) The sequence occurs at least min_chunk = 5 times in the paths.
     // (4) Every path containing either x1 or xn is compatible with this sequence.
     // If for a given x1, there is more than one x1,...,xn, we choose the longest.
     // Now we traverse all x1,..,xn as above, from longest to shortest.
     // Then:
     // (A) Add an edge x, the concatenation of x1,...,xn.
     // (B) Delete x1 and xn.
     // (C) Any path containing x1 or xn is altered to contain x instead.
     // (D) For each i = 2,...,n-1, if each path containing xi is compatible
     //     with x1,...,xn, delete xi and replace the path by x.
     // (E) Given a contiguous subsequence s of x2,...,xn-1, if every path containing
     //     it is compatible with x1,...,xn, if a path contains s, replace it by x.

     void ChunkEdges( const long_logging& logc );

     void WordifyAlt2( const long_logging_control& log_control,
          const long_logging& logc );

     void Hookup( const long_logging& logc );

     void KillWeakExits2( const long_heuristics& heur, const long_logging& logc );

     // PopBubbles.  Look for edges e: v --> w whose forward or reverse multiplicity 
     // is <= max_pop_del, and such that both are <= max_pop_del2.  Then find all 
     // alternative paths from v to w whose length in kmers is within 
     // min_pop_ratio * delta_kmers of len(e), and whose forward 
     // and reverse multiplicities are > max_pop_del.  Suppose that there is a
     // unique such alternative path p, that its length is within delta_kmers of 
     // len(e), and that its forward multiplicity is >= min_pop_ratio * mult_fw(e), 
     // or its reverse multiplicity is >= min_pop_ratio * mult_rc(e).  Then delete 
     // edge e and replace all occurrences of it by p.
     struct BubbleParams
     {
         BubbleParams( double maxPopDel, double maxPopDel2, double minPopRatio,
                             double deltaKmers )
         : max_pop_del(maxPopDel), max_pop_del2(maxPopDel2),
           min_pop_ratio(minPopRatio), delta_kmers(deltaKmers)
         {}

         double max_pop_del;
         double max_pop_del2;
         double min_pop_ratio;
         double delta_kmers;
     };
     void PopBubbles( BubbleParams const& params, unsigned nThreads,
                             const long_logging_control& log_control,
                             const long_logging& logc, const long_heuristics& heur );

     struct BubbleAux
     {
         vec<int> to_left;
         vec<int> to_right;
         vec< vec< pair<int,int> > > paths_index;
         vec<fix64_6> mult_fw;
         vec<fix64_6> mult_rc;
     };
     bool IsPoppable( int edgeId, BubbleParams const& parms,
                         BubbleAux const& aux,
                         long_logging_control const& log_control,
                         const long_logging& logc,
                         vec<int>* pSubstPath, const long_heuristics& heur ) const;

     // MakeLoops.  Given a pair of edges between two vertices, replace them by
     // a loop.
     //
     // * --*--> v --e1--> w --*--> *
     //            <--e2--

     void MakeLoops( const long_logging& logc );

     // UnrollLoops.  Attempt to find all instantiations of a simple loop.

     void UnrollLoops( const long_logging& logc );

     // Gulp.  Whenever we have u --e--> v --f1,f2--> w1,w2, replace the three
     // edges by two edges ef1, ef2, provided that the length of e in kmers is
     // at most 20.  And the other way too.

     void Gulp( const long_logging_control& log_control,
          const long_logging& logc );

     // Ungulp.  Reverse gulping, using the same 20 threshold.
     
     void Ungulp( const long_logging& logc );

     // TrimHangingEnds.  Remove edges that 'go nowhere'.

     void TrimHangingEnds( const int max_del, const double junk_ratio,
          const long_heuristics& heur, const long_logging& logc );

     // DeleteLowCoverage.  Whenever there is a branch, with one branch having
     // coverage <= 2, and the other branch having coverage at least 5 times higher,
     // delete the low coverage branch.

     void DeleteLowCoverage( const long_heuristics& heur,
          const long_logging_control& log_control,
          const long_logging& logc );

     // Check all bubbles and delete unlikely branches. 

     void DivineBubbles( const int L, const vecbasevector& bases, 
          const vecqualvector& quals, const long_heuristics& heur, 
          const long_logging& logc);
     
     // Check all bubbles and delete unlikely branches. 
     void DivineSingleMutationBubbles( const int L, const vecbasevector& bases, 
          const vecqualvector& quals, const long_heuristics& heur, 
          const long_logging& logc);

     // Check if their are variants in the assembly.  Record all the detected
     // SNP branches. 
     void DetectVarients( const vecbasevector& bases, const vecqualvector& quals,
             const long_logging& logc,
             vec<VariantSignature>* snp_bubbles, vec<int>* snp_edge_list) const;

     // ***** METHODS TO FIX BROKEN DATA *****

     // UniqueOrderPaths: fix intermediate state in which paths are not
     // unique-ordered.

     void UniqueOrderPaths( );

     // FixWeights.  Make sure that weights are symmetric with respect to the
     // involution.  
     //
     // NOTE THAT MOST OF THESE CALLS ARE PROBABLY NO LONGER NEEDED!

     void FixWeights( const long_logging& logc );

     // TruncatePaths.  If a path contains unused edges, replace it by the longest
     // subpath (in kmers) that does not contain unused edges.  In the case of a 
     // tie we delete the path.

     void TruncatePaths( const long_logging& logc );

     // Bootstrap.  Given a SupportedHyperBasevector shb0, map 
     // HyperBasevector(*this) onto it and then lift the paths and weights of
     // shb0 to *this.  Normally this would be applied in the case where shb0 is
     // a subgraph of a unipath graph.

     void Bootstrap( const SupportedHyperBasevector& shb0,
          const long_logging& logc );

     // DeleteUnusedPaths.  Remove paths (and weights) that have dead edges in
     // them.  This is for fixing a SupportedHyperBasevector after deleting
     // some of its components, as in 
     // HyperBasevector::DeleteReverseComplementComponents, but can also be used
     // to delete all paths containing any dead edges.

     void DeleteUnusedPaths( );

     // TransformPaths.  Suppose given edge sequences x1,...,xn, which are converted
     // into new edges y1,...,yn.  Then transform paths and weights appropriately.  
     // Note that this may not correctly handle the case where a path contains 
     // multiple xs in an overlapping fashion.  Note also that ambiguously 
     // transformable paths are killed, as are 'illegal' paths.

     void TransformPaths( const vec< vec<int> >& x, const vec<int>& y,
          vec< vec< pair<int,int> > >& paths_index );
     void TransformPaths( const vec< vec<int> >& x, const vec<int>& y )
     {    vec< vec< pair<int,int> > > paths_index;
          TransformPaths( x, y, paths_index );    }

     // ***** CLEANUP METHODS *****

     void RemoveDeadEdgeObjects0( );
     void RemoveDeadEdgeObjects( );
     void RemoveUnneededVertices0( vec< triple<int,int,int> >& merges );
     void RemoveUnneededVertices( );
     void DeleteReverseComplementComponents( const long_logging& logc, 
          const int64_t iDirectionalSortFactor=0 );
     void RemoveSmallComponents( const int K );
     void RemoveSmallMostlyAcyclicComponents( const long_logging& logc, 
          const int max_small_comp = 1500 );
     void RemoveSmallComponents2( const long_logging& logc, 
          const int max_small_comp );
     void RemovePathsWithoutReverseComplements();
     void RemovePairsWithoutReverseComplements();

     // **** WRITING METHODS *****

     // DumpFiles(HEAD): generate HEAD.{shbv,dot,support,fasta}.

     void DumpFiles( const String& head, const long_logging_control& log_control,
          const long_logging& logc ) const;
     void DumpFilesStandard( const long_logging_control& log_control,
          const long_logging& logc, const int id ) const;

     void DumpDot( const String& head, const vec<Bool>& invisible,
          const vec<String>& edge_color, const long_logging& logc,
          const Bool hide_inv = True, 
          const vec<String>& edge_names = vec<String>( ) ) const;

     void writeBinary( BinaryWriter& writer ) const;
     void readBinary( BinaryReader& reader );

};

template<> struct Serializability<SupportedHyperBasevector>
{ typedef SelfSerializable type; }; 

// Overlap.  Determine if two vectors v, w of integers overlap with the given
// offset.  Positive offset corresponds to w being to the right of v.

Bool Overlap( const vec<int>& v, const vec<int>& w, const int o );

void OrientToReference( SupportedHyperBasevector& shb, const vecbasevector& genome,
     const long_logging& logc );

void TraceEdges( const SupportedHyperBasevector& shb, const String& TRACE_EDGES,
     const vecbasevector& bases, const vecqualvector& quals );

void AnalyzeAssembly( const SupportedHyperBasevector& shb, const vecbasevector& G,
     const int LG, const VecIntPairVec& Glocs );

void AssessAssembly( const String& SAMPLE, const SupportedHyperBasevector& shb,
     const HyperEfasta& he, const vec<Bool>& hide, const String& TMP, 
     const ref_data& ref, const String& HUMAN_CONTROLS, 
     const long_logging& logc, const uint NUM_THREADS, RefTraceControl RTCtrl=RefTraceControl() );

void ReportAssemblyStats( const SupportedHyperBasevector& shb );

void DumpEfastaAssembly( const SupportedHyperBasevector& shb, const HyperEfasta& he,
     const vec<int>& inv, const vec<Bool>& hide, const String& OUT_HEAD,
     const long_logging& logc );

void CountCov( const SupportedHyperBasevector& shb, const String& TMP, 
     const int gp1 );

void RemoveNegatives( vec<int>& p );

#endif
