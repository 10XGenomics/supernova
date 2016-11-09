///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////


// RefTrace.  Assess HyperBasevector assembly hb using reference sequence G, by
// aligning each reference contig G[g] to hb.  Each such alignment consists of
// aligning segments, punctuated by gaps.  A number of choices are involved in
// defining these alignments.  End gaps are treated differently: we simply tally
// the number of unaligned bases.  Some gaps are likely not true gaps, but rather
// artifacts of the alignment process.
//
// HOW IT WORKS:
//
// 1. We build a larger HyperBasevector hbp in which every edge lives in both
//    orientations.
// 2. Each edge in hbp is aligned to the genome, allowing for multiple placements
//    (and subject to several heuristics).
// 3. To each alignment we associate its start and stop points on the genome.
//    The stop point is defined by finding the point on the genome that corresponds
//    to the point on the edge that is K-1 bases from its right end.  (This is not
//    well-behaved if there is an indel at this point.)
// 4  Because of 'random' choices in alignment, there may be small inconsistencies
//    between these points arising from an edge entering and an edge exiting a given
//    vertex.  We 'twiddle' the start and stop points slightly to get them to agree
//    in some cases (so long as cycles are not introduced).
// 5. We form a graph in which the vertices are certain triples (g,p,v) where g is
//    a reference contig id, p is a position on it, and v is a vertex id in hbp.
//    Each alignment defines two vertices and an edge between them.
// 6. We add edges to the graph to represent gaps, whenever we observe sink-source
//    'gaps' of size between -10000 and 10000 bases.  In the case of negative gaps,
//    this introduces some inaccuracy because of overcounting, although we eliminate
//    some of this below.
// 7. Then, without recomputing sources and sinks, we add 'master' sources and sinks
//    for each chromosome, and connect the master sources to the old sources and the
//    old sinks to the master sinks.
// 8. Each edge is assigned a penalty, 1*e + 100*g + 0.1*d, where
//    - e = alignment error count for the edge;
//    - g = the number of gaps;
//    - d = the number of gap bases.
// 9. We find shortest paths through the graph.  By construction there is a single
//    shortest path for each chromosome, which is our final answer.  Note that there
//    may be errors in common for each pair of edges * ----e1----> * ----e2----> * 
//    since they overlap by K-1 bases, which would result in overcounting.  Negative
//    gaps as described above could have the same effect.  This introduces some 
//    inaccuracy into the finding of shortest paths.  However having found the 
//    shortest paths, we then compensate (although not fully) by finding the total
//    nonredundant error count for a given path, in an appropriate sense.

// ISSUES:
//
// 1. As described above, there are several inconsistencies in the calculation
//    arising from indels and incompatible alignments.
// 2. Circular chromosomes are incorrectly handled.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS
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

namespace { // open anonymous namespace

int Penalty( const int e, const int g, const int d )
{    return int( round( 100 * ( e + 100*g + d/100.0 ) ) );    }

// Use this penalty function to generate longer best path, useful for variant calling.
int PenaltyForVariantCalling( const int e, const int g, const int d )
{    return 3*e + 100*g + d;    }

// Generate the sequences from the best path. The first K-1 bases will always be
// trimmed except for the first edge. This can be unwanted behavior if the path
// is continuous.
basevector BestGlobal(const vec<vec<int>>& best_path, const vec<vec<int>>& eids, 
        int g, const HyperBasevector& hbp, 
        const EdgePlacements& edge_placements ) 
{
    int left_trim = edge_placements.vedata[eids[g].front( )].fourth.first;
    int right_trim = hbp.EdgeObject( best_path[g].back( ) ).isize( )
         - edge_placements.vedata[ eids[g].back( ) ].fourth.second;
    basevector best_bpath = hbp.EdgeObject( best_path[g][0] );
    for ( int i = 1; i < best_path[g].isize( ); i++ )
    {    best_bpath = TrimCat( hbp.K( ), 
              best_bpath, hbp.EdgeObject( best_path[g][i] ) );    }
    return basevector(best_bpath, left_trim, best_bpath.isize( ) - 
            left_trim - right_trim);
}

// Compare the best_path sequences to the reference and show the difference.
void ShowReftraceEvents(const basevector& best_bpath, const basevector& ref, int g,
        int& meta_events, ostream& out)
{
    alignment al;
    const basevector& b1 = best_bpath, &b2 = ref;
    SmithWatAffineParallel( b1, b2, al );
    align a(al);
    vec< pair<int,int> > P1, P2;
    DecomposeAlign( a, b1, b2, P1, P2 );
    for ( int z = 0; z < P1.isize( ); z++ ) {    
        basevector B1( b1, P1[z].first, 
                P1[z].second - P1[z].first );
        basevector B2( b2, P2[z].first, 
                P2[z].second - P2[z].first );
        align A = a.TrimmedTo1( 
                P1[z].first, P1[z].second - P1[z].first );
        vec<int> mgg = A.MutationsGap1Gap2( B1, b2 );
        int start = P2[z].first, stop = P2[z].second;
        out << "\nMETA-EVENT " << ++meta_events << " (at " << g
            << "." << start << "-" << stop << ", subs = "
            << mgg[0] << ", dels = "
            << mgg[1] << ", ins = " << mgg[2] << ")\n";
        PrintVisualAlignmentClean(
                False, out, B1, b2, A );
    }
}

// Check if the best-path sequence covers significant portion of the reference.
// The criterion is that 30% of the reference base is covered and the largest
// gap size is no larger than 2000.
bool HasGoodCoverage(const vec<int>& best_path, vec<int>& eids,
        const EdgePlacements& edge_placements, const basevector& Gg) 
{
    const double MinBestPathCoverage = 0.3;
    const int MaxMaxGapBases = 2000;
    const HyperBasevector& hbp = edge_placements.hbp;

    int start = edge_placements.vedata[eids.front()].first.second;
    int stop = edge_placements.vedata[eids.back()].second.second;
    double cov_rate1 = (double)(stop - start)/(double)Gg.size();

    vec<int> to_left, to_right;
    hbp.ToLeft(to_left); hbp.ToRight(to_right);
    int max_gap_bases = 0;
    for (size_t i = 1; i < best_path.size(); i++) {
        if (to_left[best_path[i]] != to_right[best_path[i-1]]) {
            int ref_start = edge_placements.vedata[eids[i]].first.second;
            int last_stop = edge_placements.vedata[eids[i-1]].second.second;
            max_gap_bases = max(max_gap_bases, ref_start - last_stop);
        } 
    }
    return cov_rate1 >= MinBestPathCoverage && max_gap_bases <= MaxMaxGapBases;
}

// Given the best path, find the start and stop position of each edge on the
// reference.  It re-calculate the best path if the original one has poor
// coverage, bad design.
void FindPositionOnRef( vec<vec<vec<int>>>& best_path_s,
              vec<vec<vec<pair<int,int>>>>& ref_pos_all_s,
              vec<vec<int>> best_path2, vec<vec<int>> eids2,
              int best_g, const EdgePlacements& edge_placements,
              int min_dist, int max_dist, int verbosity,
              ostream& out)
{
    if (! HasGoodCoverage(best_path2[best_g], eids2[best_g], edge_placements, 
                edge_placements.G[best_g])) {
        out << "Best path has low quality. Try alternative "
            "penalty function that do not over-penalize gaps." << endl;
        GraphZ gz2(edge_placements, &PenaltyForVariantCalling);
        vec< vec<int> > spaths2; 
        vec< triple<int,int,int> > spaths_egd2;
        vec< pair<int,int> > spaths_gg_pen2;
        int best_g2;
        gz2.FindShortestPath(min_dist, max_dist, spaths2, spaths_egd2,
                spaths_gg_pen2, out, verbosity);
        gz2.FindBestEdgePath(spaths_egd2, spaths2, best_path2, eids2, best_g2);
    }

    int gsize = edge_placements.G.size();
    vec<vec<pair<int,int>>> ref_pos_all(gsize); 
    for (size_t g = 0; g < ref_pos_all.size(); g++) 
        ref_pos_all[g].resize(best_path2[g].size());

    for (int g = 0; g < gsize; g++) {    
        for (int i = 0; i < (int)best_path2[g].size(); i++) {
            int eid = eids2[g][i];
            int start_on_ref = edge_placements.vedata[eid].first.second;
            int start_on_edge = edge_placements.vedata[eid].fourth.first;
            ref_pos_all[g][i].first = start_on_ref - start_on_edge;
            ref_pos_all[g][i].second = ref_pos_all[g][i].first + 
                edge_placements.hbp.EdgeLength(best_path2[g][i]);
        }
    }
    // Below should be removed so that each g should have one best path
    best_path_s.resize(gsize);
    ref_pos_all_s.resize(gsize);
    for(int gg=0; gg<gsize; ++gg){
        best_path_s[gg].push_back( best_path2[gg]);
        ref_pos_all_s[gg].push_back(ref_pos_all[gg]);
    }
}

// The monster Reftrace and variant calling function. It performs multiple tasks
// depending on the input from long_logging and other switches. The main steps:
// 1. Perform reftrace to find the best path in the assembly graph that the
//    reference will thread through. The algorithm is explained on top of this
//    page. The data is wrapped into two classes of EdgePlacements that contains
//    the alignments of edges to reference, and GraphZ that contains
//    edge connections from left to right ends of the reference.
// 2. Print out the best path sequence, show reference events, or perform
//    variant calling. Note that for variant calling, we will run the reftrace
//    again in slight different parameter if the original best path coverage is
//    bad. We then call variant using the new copy of the best path, without
//    altering the original one, which is also used for assessment.
// 3. Report the assembly assessment from the reftrace result.
RefTraceResults RefTraceInternal( const ref_data& ref,
     const HyperBasevector& hb, const vec<int>& inv, 
     const int verbosity, const long_logging& logc, 
     ostream& out, ostream& callsOut, RefTraceHeuristics rth,
     const String& BEST_GLOBAL_OUT, const Bool fix_bug,
     const RefTraceControl* p_ref_trace_control, 
     const ReadOriginTracker* p_read_tracker)
{    

    const vec<HyperBasevector>& GH = ref.GH; 
    const vec<bool>& is_circular = ref.is_circular;
    const vecbasevector& Gplus = ref.G3plus;
    const vec<int>& Gplus_ext = ref.Gplus_ext;

     double rclock = WallClockTime( );
     int K = hb.K( );
     if ( verbosity >= 2 )
          out << "assembly has " << hb.EdgeObjectCount( ) << " edges" << endl;
     ForceAssertEq( hb.EdgeObjectCount( ), inv.isize( ) );

     // Test and expand reference.

     LinearRef linear_ref(GH);

     if ( verbosity >= 1 ) out << Date( ) << ": building hbp" << endl;
     HyperBasevector hbp; 
     vec<pair<int,Bool>> hbp_to_hb;
     CreateHBPlus(hb, inv, hbp, hbp_to_hb);
     if ( verbosity >= 2 )
     {    out << "'doubled' assembly has " << hbp.EdgeObjectCount( ) << " edges"
               << endl;    }

     // Heuristics.

     const int L = 40;
     const double min_cov_frac = 0.5;
     const int min_dist = -10000;
     const int max_dist = 10000;

     // Align edges of hbp to reference.
     EdgePlacements edge_placements_org(hbp, hbp_to_hb, linear_ref.Seqs());
     edge_placements_org.AlignEdgesToRef<L> (min_cov_frac, rth.max_error_rate,
             rth.max_offset_diff, rth.min_group_frac, rth.offset_add,
             rth.min_group_save, fix_bug, logc.REFTRACE_VARIANTS, verbosity, out);

     // Try to twiddle alignment ends.  Be careful not to create cycles.
     edge_placements_org.Twiddle(rth.max_twiddle);
     const bool bOldMethodFailed=edge_placements_org.vedata.size()==0;

     EdgePlacements edge_placements_new(hbp, hbp_to_hb, linear_ref.Seqs());
     if(bOldMethodFailed){
          if(verbosity){
               out << "\nWARNING: AlignEdgesToRef<L>+Twiddle failed, redoing using new code\n" << std::endl;
          }
          edge_placements_new.AlignEdgesToRefExp<L> (verbosity, out);
          edge_placements_new.RemoveBadPlacements();
          edge_placements_new.TwiddleSmart();
     }

     const EdgePlacements& edge_placements=(!bOldMethodFailed)?edge_placements_org:edge_placements_new;


     vec< vec<int> > spaths;
     vec< triple<int,int,int> > spaths_egd;
     vec< pair<int,int> > spaths_gg_pen;

     // Build graph.
     GraphZ gz(edge_placements, &Penalty);
     gz.FindShortestPath(min_dist, max_dist, spaths, spaths_egd, spaths_gg_pen,
             out, verbosity);

     REPORT_TIME( rclock, "used tracing reference - main" );

     // Print best global alignment.

     int meta_events = 0;
     if ( BEST_GLOBAL_OUT != "" || logc.SHOW_REFTRACE_EVENTS || logc.REFTRACE_VARIANTS)
     {    double gclock = WallClockTime( );
          const vecbasevector& G = linear_ref.Seqs();
          if ( spaths.empty( ) )
          {    out << "BEST_GLOBAL_OUT: there are no paths, giving up." << endl;
               _exit(1);    }

          vec< vec<int> > best_path( G.size( ) ), eids( G.size( ) );
          int best_g = -1;
          gz.FindBestEdgePath(spaths_egd, spaths, best_path, eids, best_g);

          if ( BEST_GLOBAL_OUT != "" ) {
              basevector best_bpath = BestGlobal(best_path, eids, 
                      best_g, hbp, edge_placements);
              vecbasevector best;
              best.push_back(best_bpath);
              best.push_back(G[best_g]);
              best.WriteAll(BEST_GLOBAL_OUT);
          }

          // call variants guided by the best path alignment on reference

          if (logc.REFTRACE_VARIANTS) {
              String sVarFile = p_ref_trace_control->getVariantOutFile();
              String sRefFile = p_ref_trace_control->getRefHead();

              vec<int> Gplus_ext_zero(G.size(), 0);
              const vecbasevector& Gplus_sel = (Gplus.empty() ? G : Gplus);
              const vec<int>& Gplus_ext_sel = (Gplus.empty() ? Gplus_ext_zero : Gplus_ext);

              vec<vec< vec<int> >> best_path_s( G.size( ) );
              vec<vec<vec<pair<int,int>>>> ref_pos_all_s(G.size());

              FindPositionOnRef(best_path_s, ref_pos_all_s, best_path, eids, best_g,  
                      edge_placements, min_dist, max_dist, verbosity, out);
              ReftraceVariants(out, callsOut, G, Gplus_sel, Gplus_ext_sel, hbp,
                      hbp_to_hb, ref_pos_all_s, best_path_s, sVarFile, sRefFile,
                      *p_ref_trace_control, p_read_tracker, &logc);
          }

          // Show RefTrace events.  This is part of what is in
          // auxmain/AssessBestAlign.cc.

          if (logc.SHOW_REFTRACE_EVENTS) {    
               for ( int g = 0; g < (int) G.size( ); g++ ) {    
                   basevector best_bpath = BestGlobal(best_path, eids, 
                           g, hbp, edge_placements);
                   ShowReftraceEvents(best_bpath, G[g], g, meta_events, out);
               }
          }
          REPORT_TIME( gclock, "used tracing reference - global" ); 
     }
     // Report gaps.

     double tclock = WallClockTime( );
     if ( verbosity >= 1 )
     {    out << "\ngaps:\n";
          for ( int i = 0; i < spaths.isize( ); i++ )
          {    const vec<int>& p = spaths[i];
               for ( int i = 0; i < p.isize( ) - 1; i++ )
               {    int v = p[i], w = p[i+1];
                    int min_penalty = 1000000000, best_edge = -1;
                    for ( int l = 0; l < gz.Z.From(v).isize( ); l++ )
                    {    if ( gz.Z.From(v)[l] != w ) continue;
                         int p = gz.Z.EdgeObjectByIndexFrom( v, l );
                         if ( p < min_penalty )
                         {    min_penalty = p;
                              best_edge = gz.Z.EdgeObjectIndexByIndexFrom( 
                                   v, l );    }    }
                    if ( best_edge < 0 ) continue;
                    if ( gz.egd[best_edge].second == 0 ) continue;
                    int start = gz.verts[v].second, stop = gz.verts[w].second;
                    out << gz.verts[v].first << "." << start << "-" << stop
                         << " [len = " << stop-start << "]\n";    }    }    }

     // Extract metrics from paths.

     int total_errs = 0, total_gaps = 0, total_gap_bases = 0;
     const int infinity = 100000000;
     vec<int> e_best( GH.size( ), infinity );
     vec<int> g_best( GH.size( ), infinity );
     vec<int> d_best( GH.size( ), infinity );
     for ( int i = 0; i < spaths.isize( ); i++ )
     {    int m = spaths_gg_pen[i].first;
          int x = linear_ref.Source(m);
          e_best[x] = Min( e_best[x], spaths_egd[i].first );
          g_best[x] = Min( g_best[x], spaths_egd[i].second );
          d_best[x] = Min( d_best[x], spaths_egd[i].third );    }
     for ( int x = 0; x < (int) GH.size(  ); x++ )
     {    if ( e_best[x] < infinity )
          {    total_errs += e_best[x];
               total_gaps += g_best[x];
               total_gap_bases += d_best[x];    }
          else
          {    total_gaps++;
               for ( int i = 0; i < linear_ref.N(); i++ )
               {    if ( linear_ref.Source(i) == x )
                    {    total_gap_bases += linear_ref.Seq(i).size( );
                         break;    }    }    }    }

     // Print summary statistics.

     int penalty = Penalty( total_errs, total_gaps, total_gap_bases );
     int64_t gsize = 0;
     for ( int g = 0; g < linear_ref.N(); g++ )
          gsize += linear_ref.Seq(g).size( );
     double errors_per_Mb = 1000000.0 * double(total_errs) / double(gsize);
     int free_gaps = 0;
     for ( int g = 0; g < linear_ref.N(); g++ )
          if ( !is_circular[g] ) free_gaps += 2;
     double gaps_per_Mb
          = double(total_gaps - free_gaps) / ( double(gsize)/1000000.0 );
     int gap_bases_per
          = int( round( double(total_gap_bases) / double( linear_ref.N() ) ) );

     if(logc.REFTRACE_VARIANTS_SUMMARY){
         out << penalty << " penalty (p = 1*e + 100*g + 0.1*d)" << endl;
         out << "-[" << setiosflags(ios::fixed) << setprecision(2) << errors_per_Mb
              << resetiosflags(ios::fixed) << "%%%]- errors"
              << " (e = errs = " << total_errs << ")" << endl;
         if (logc.SHOW_REFTRACE_EVENTS)
              out << "-[" << meta_events << "]- error meta-events" << endl;
         out << "-[" << setiosflags(ios::fixed) << setprecision(2) << gaps_per_Mb
              << resetiosflags(ios::fixed) << "]- excess gaps per Mb "
              << "(g = total gaps = " << total_gaps
              << ", \"excess\" = " << total_gaps - free_gaps << ")" << endl;
         out << "-[" << gap_bases_per << "]- gap bases per reference contig "
              << "(d = gap bases = " << total_gap_bases << ")\n";
     }

     if ( logc.PERF_STATS )
     {    PerfStat::log( ) << std::fixed << setprecision(0)
               << PerfStat( "gap_bases_per",
               "gap bases per reference contig", gap_bases_per );
          PerfStat::log( ) << std::fixed << setprecision(2)
               << PerfStat( "errors_per_Mb", "errors per Mb", errors_per_Mb );
          PerfStat::log( ) << std::fixed << setprecision(2)
               << PerfStat( "gaps_per_Mb", "excess gaps per Mb", gaps_per_Mb );    }
     REPORT_TIME( tclock, "used tracing reference - tail" );
     return RefTraceResults( penalty, total_gaps, meta_events );    }

}// close anonymous namespace

RefTraceResults RefTrace( const ref_data& ref,
     const HyperBasevector& hb, const vec<int>& inv, 
     const int verbosity, const long_logging& logc, 
     ostream& out, RefTraceHeuristics rth,
     const String& BEST_GLOBAL_OUT, const Bool fix_bug)
{
    NullOStream nullOS;
    return RefTraceInternal(ref, hb, inv, verbosity, logc, out, nullOS,
            rth, BEST_GLOBAL_OUT, fix_bug, NULL, NULL);
}


RefTraceResults RefTraceAndCallVaraint( const ref_data& ref,
     const HyperBasevector& hb, const vec<int>& inv, 
     const int verbosity, const long_logging& logc, 
     ostream& out, ostream& callsOut, RefTraceHeuristics rth,
     const String& BEST_GLOBAL_OUT, const Bool fix_bug,
     const RefTraceControl& ref_trace_control, 
     const ReadOriginTracker* p_read_tracker)
{
    return RefTraceInternal(ref, hb, inv, verbosity, logc, out, callsOut,
            rth, BEST_GLOBAL_OUT, fix_bug, &ref_trace_control, p_read_tracker);
}

