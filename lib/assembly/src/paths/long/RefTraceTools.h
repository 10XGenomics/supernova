///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
#ifndef REFTRACE_TOOLS_H
#define REFTRACE_TOOLS_H
// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "CoreTools.h"
#include "paths/HyperBasevector.h"
#include "pairwise_aligners/SmithWatAffine.h"
#include "pairwise_aligners/SmithWatBandedA.h"
#include "paths/long/MakeKmerStuff.h"
#include "PrintAlignment.h"
#include "paths/long/RefTrace.h"


// Create a HyperBasevector hbp that equals hb plus its reverse complement.
// However only do this for components that need it.

void CreateHBPlus(const HyperBasevector& hb, const vec<int>& inv,
        HyperBasevector& hbp, vec<pair<int,Bool>>& hbp_to_hb);


// Linearized reference sequences. Expand the paths from the source to the sink
// of the reference graph, built a vecbasevector of expended sequence, and
// record the origion of of the chromosome id.

class LinearRef {
public:
    LinearRef(const vec<HyperBasevector>& GH,const vec<bool>& c=vec<bool>());
    int   N() const                    { return G.size(); }
    int   Source(int g) const          {return G_source[g]; }
    bool  IsDoubled(int g) const          {return isDoubled[g]; }
    const basevector& Seq(int g) const { return G[g]; }
    const vecbasevector& Seqs() const  { return G; }
private:
    vec<int> G_source;
    vecbasevector G;
    vec<bool> isDoubled;
};


// Some data structures.
// The structure vedata has the following structure:
// { (genome_tig_id, start_pos_on_genome_tig, left_vertex_of_edge_in_hbp ),
//   (genome_tig_id, stop_pos_on_genome_tig-K+1, right_vertex_of_edge_in_hbp ),
//   (hbp_edge_id, error_count),
//   (start_pos_on_hbp_edge, stop_pos_on_hbp_edge) }.

class EdgePlacements {
public:
    EdgePlacements(const HyperBasevector& hbp, const vec<pair<int,Bool>>& hbp_to_hb,
            const vecbasevector& G) : hbp(hbp), hbp_to_hb(hbp_to_hb), G(G) {}

    // Align edges of hbp to reference.
    template<int L> 
    void AlignEdgesToRef( 
         // heuristics:
              const double min_cov_frac, const double max_error_rate,
              const int max_offset_diff, const double min_group_frac,
              const int offset_add, const int min_group_save, const Bool fix_bug,
         // logging:
              bool REFTRACE_VARIANTS, const int verbosity, ostream& out );

    template <int L>
    void AlignEdgesToRefExp(const int verbosity, ostream& out);
    void RemoveBadPlacements();
    void Twiddle(const int max_twiddle);
    void TwiddleSmart();

    // Generate the matching sequences from the best path.
    basevector BestSeq(const vec<int>& best_path, const vec<int>& eids
                      , const vec<std::pair<int,int>>& limits
                      , vec<std::tuple<int64_t,int64_t,int,int64_t,int64_t,int64_t,int64_t>>& coors_edge);

public:
    const HyperBasevector& hbp;
    const vec<pair<int,Bool>>& hbp_to_hb;
    const vecbasevector& G;

    vec< quad< triple<int,int,int>, triple<int,int,int>, 
        pair<int,int>, pair<int,int> > > vedata;
    vec<align> aligns;
    vec<int> xto_left, xto_right;
private:
    int CorrelatePositionsAlways(const align& a, const int x1)const;
};

class GraphZ {
public:
    typedef int (*PenaltyFuncT)(int, int, int);
    GraphZ(const EdgePlacements& ep, PenaltyFuncT pf) 
        : edge_placements(ep), hbp(ep.hbp), hbp_to_hb(ep.hbp_to_hb), G(ep.G)
    { Penalty = pf; }

    void FindShortestPath(const int min_dist, const int max_dist,
            vec< vec<int> >& spaths, vec< triple<int,int,int> >& spaths_egd,
            vec< pair<int,int> >& spaths_gg_pen,
            ostream& out, int verbosity = 0);

    // Find the corresponding best path in hbp edges. 
    void FindBestEdgePath( const vec< triple<int,int,int> >& spaths_egd,
            const vec< vec<int> >& spaths,
            vec<vec<int>>& best_path, vec<vec<int>>& eids, int& best_g) ;
public:
    const EdgePlacements& edge_placements;
    const HyperBasevector& hbp;
    const vec<pair<int,Bool>>& hbp_to_hb;
    const vecbasevector& G;
    PenaltyFuncT Penalty;
    vec< triple<int,int,int> > verts;
    vec< triple< int, int, pair<int,int> > > edges;
    vec< triple<int,int,int> > egd; 
    digraphE<int> Z;
private:
    void BuildGraph(const int verbosity, ostream& out);

    void AddGapEdges(const int min_dist, const int max_dist, const int verbosity, 
            ostream& out ,const bool bPreserveDisconnectedComponents=false);

    void AddConnectedGapEdges(const int min_dist, const int max_dist, const int verbosity,
            ostream& out ,const bool bPreserveDisconnectedComponents=false);

    void AddSpecialVerts( const int K, const vec<int>& sources, const vec<int>& sinks,
         const bool bPreserveDisconnectedComponents=false);

    void AnnouncePaths( const vec< vec<int> >& spaths, const int K, 
         const vec< triple<int,int,int> >& spaths_egd, 
         const int verbosity, ostream& out ) const;

    void FindShortestPathBetween( const int this_source, const int this_sink, 
         const digraphE<int>& ZS, const vec<int>& suc, 
         vec< vec<int> >& spaths, vec< triple<int,int,int> >& spaths_egd, 
         vec< pair<int,int> >& spaths_gg_pen,
         const int verbosity, ostream& out ) const;

    void MakeZDot(ostream& os);
};


template<int L> 
void EdgePlacements::AlignEdgesToRef( 
          const double min_cov_frac, const double max_error_rate,
          const int max_offset_diff, const double min_group_frac,
          const int offset_add, const int min_group_save, const Bool fix_bug,
     // logging:
          bool REFTRACE_VARIANTS, const int verbosity, ostream& out )
{
     // Setup for alignment.
     vecbasevector all(G);
     vec< triple<kmer<L>,int,int> > kmers_plus;
     MakeKmerLookup0( all, kmers_plus );
     vec< kmer<L> > kmers( kmers_plus.size( ) );
     for ( int64_t i = 0; i < kmers_plus.jsize( ); i++ )
          kmers[i] = kmers_plus[i].first;
     hbp.ToLeft(xto_left), hbp.ToRight(xto_right);

     // Go through the edges of the (doubled) assembly.

     #pragma omp parallel for schedule(dynamic,1)
     for ( int i = 0; i < hbp.EdgeObjectCount( ); i++ )
     {    const basevector& e = hbp.EdgeObject(i);

          // For each kmer in the edge, find its hits to the reference and find
          // the kmers having the most hits.

          int nkmers = e.isize( ) - L + 1;
          vec< triple<int64_t,int64_t,int64_t> > locs(nkmers);
          vec<int> pos( nkmers, vec<int>::IDENTITY );
          kmer<L> x;
          for ( int j = 0; j < nkmers; j++ )
          {    x.SetToSubOf( e, j );
               int64_t low = LowerBound(kmers, x), high = UpperBound(kmers, x);
               locs[j].first = high - low;
               locs[j].second = low, locs[j].third = high;    }
          if (fix_bug) ReverseSortSync( locs, pos );
          else SortSync( locs, pos );

          // Determine cutoff 'top'.

          double mcf = min_cov_frac;
          if ( REFTRACE_VARIANTS ) mcf = 0.6;
          int t = int( floor( nkmers * mcf ) ), top;
          for ( top = t + 1; top < nkmers; top++ )
               if ( locs[top].first > locs[t].first ) break;

          // Find the associated offsets.

          vec< pair<int,int> > offset;
          for ( int j = 0; j < top; j++ )
          {    for ( int64_t m = locs[j].second; m < locs[j].third; m++ )
               {    int g = kmers_plus[m].second, o = kmers_plus[m].third - pos[j];
                    offset.push( g, o );    }    }
          Sort(offset);

          // Form offsets into groups.

          vec< triple< int, int, pair<int,int>  > > og;
          int mod = max_offset_diff;
          if ( REFTRACE_VARIANTS ) mod = 500;
          for ( int j = 0; j < offset.isize( ); j++ )
          {    int k;
               for ( k = j + 1; k < offset.isize( ); k++ )
               {    if ( offset[k].first != offset[j].first ) break;
                    if ( offset[k].second - offset[k-1].second > mod )
                         break;    }
               og.push( k - j, offset[j].first, 
                    make_pair( offset[j].second, offset[k-1].second ) );
               j = k - 1;    }
          ReverseSort(og);
          if ( verbosity >= 4 )
          {    
               #pragma omp critical
               {    out << "\noriginal edge " << hbp_to_hb[i].first << ": ";
                    PRINT4_TO( out, nkmers, top, 
                         offset.size( ), og.size( ) );    
                    for ( int j = 0; j < og.isize( ); j++ )
                         PRINT2_TO( out, j, og[j].first );    }    }

          // Filter offset groups.

          double mgf = min_group_frac;
          if ( REFTRACE_VARIANTS ) mgf = 0.65;
          int gj;
          for ( gj = 0; gj < og.isize( ); gj++ )
          {    if ( og[gj].first < min_group_save
                    && og[gj].first < mgf * og[0].first ) 
               {    break;    }    }
          og.resize(gj);
          if ( verbosity >= 3 && og.nonempty( ) )
          {    
               #pragma omp critical
               {    out << "\noffsets for edge " << i << " (hb_edge=" 
                         << hbp_to_hb[i].first << ", nkmers=" << nkmers << ")" << endl;
                    for ( int j = 0; j < og.isize( ); j++ )
                    {    out << "[" << j << "] " << og[j].second << "."
                              << og[j].third.first << "-" << og[j].third.second 
                              << " (" << og[j].first << ")" << endl;    }    }    }

          // Align.  The reason for adding to the offset is that there could be in
          // indel in the first or last L bases.

          for ( int j = 0; j < og.isize( ); j++ )
          {    int g = og[j].second;
               int off_low = og[j].third.first, off_high = og[j].third.second;
               int mid_offset = ( off_low + off_high ) / 2;
               int bandwidth 
                    = Max(mid_offset - off_low, off_high - mid_offset) + offset_add;

               // Do the alignment.  This is kludgy.  If the alignment has too 
               // many errors and the edge is long, we suspect that the problem
               // might be with a big indel, so we align using a larger bandwidth.
               // Note the unfortunate us of hardcoded constants.

               align a;
               int errors;

               if ( !REFTRACE_VARIANTS )
               {
                   const int SMA_method = 1;
                   if (SMA_method == 1) {
                       SmithWatBandedA( hbp.EdgeObject(i), G[g],
                            -mid_offset, bandwidth, a, errors, 0, 1, 1 );
                   } else if(SMA_method == 2) {
                       SmithWatAffineBanded( hbp.EdgeObject(i), G[g],
                               -mid_offset, bandwidth, a, errors );
                   } else {
                       cout << "unrecognized SMA_method" << endl;
                   }
                   if ( double(errors) / double( a.extent2( ) ) > max_error_rate )
                   {
                       // So the following code (after the continue;
                       // bandwidth=5000) was taking a ton of time
                       // (0.5-1 sec per alignment).  Also in my tests
                       // it had a very low success rate <0.5% AND its
                       // removal does not seem to impact the result.
                       // We should do something clever with the
                       // alignments (super aligner?) if we end up
                       // needing it. -- nw
                       continue;
#if 0
                        const int long_edge = 5000;
                        const int max_indel = 5000;
                        if ( hbp.EdgeLengthBases(i) < long_edge ) continue;
                        SmithWatBandedA( hbp.EdgeObject(i), G[g],
                             -mid_offset, max_indel, a, errors, 0, 1, 1 );
                        if ( double(errors) / double( a.extent2( ) ) > max_error_rate )
                             continue;
#endif
                   }
               }
               else
               {    double score = SmithWatAffineBanded( hbp.EdgeObject(i), G[g],
                           -mid_offset, bandwidth, a, errors ) / 3.0;
                    if ( verbosity >= 3 )
                    {    
                         #pragma omp critical
                         {    double err_rate = score / double( a.extent2( ) );
                              int hb_edge = hbp_to_hb[i].first;
                              int offset = -mid_offset;
                              PRINT5( hb_edge, offset, bandwidth, score, err_rate );
                                   }    }
                    double var_max_error_rate = 0.3;
                    if ( score / double( a.extent2( ) ) > var_max_error_rate )
                         continue;    }

               if ( verbosity >= 3 )
               {
                    #pragma omp critical
                    {    out << "\nalignment " << j << " of edge " << i << " ("
                              << xto_left[i] << " --> " << xto_right[i] 
                              << ", hb_edge=" << hbp_to_hb[i].first << ")" << endl;
                         vec<int> errs 
                              = a.MutationsGap1Gap2( hbp.EdgeObject(i), G[g] );
                         int mismatches = errs[0];
                         int indels = errs[1] + errs[2];
                         PRINT5_TO( out, g, a.pos2( ), a.Pos2( ), 
                              mismatches, indels );
                         if ( verbosity == 4 )
                         {    PrintVisualAlignment( True, out, hbp.EdgeObject(i),
                                   G[g], a );    }
                         if ( verbosity >= 5 )
                         {    PrintVisualAlignment( False, out, hbp.EdgeObject(i),
                                   G[g], a );    }    }    }

               // Figure out where the position e.isize( ) - K + 1 should map to
               // under the alignment.  Note that because there could be an indel
               // there, this is not necessarily a meaningful answer.

               int x1 = e.isize( ) - hbp.K( ) + 1;
               int x2 = CorrelatePositionsAlways( a, x1 );

               // Save results.

               #pragma omp critical
               {    vedata.push( make_triple( g, a.pos2( ), xto_left[i] ),
                         make_triple( g, x2, xto_right[i] ),
                         make_pair( i, errors ), make_pair( a.pos1( ), a.Pos1( ) ) );
                    aligns.push_back(a);    }    }    }    
     // Sort the output to avoid the stochastic downstream behavior of BuildGraph
     // that seems depend on the input order of the alignment data.
     SortSync(vedata, aligns);
}

// An experimental version of function to align edges to reference that
// automatically adjust heuristics for best results.

template<int L> 
void EdgePlacements::AlignEdgesToRefExp(const int verbosity, ostream& out)
{
     // Setup for alignment.
     vecbasevector all(G);
     vec< triple<kmer<L>,int,int> > kmers_plus;
     MakeKmerLookup0( all, kmers_plus );
     vec< kmer<L> > kmers( kmers_plus.size( ) );
     for ( int64_t i = 0; i < kmers_plus.jsize( ); i++ )
          kmers[i] = kmers_plus[i].first;
     hbp.ToLeft(xto_left), hbp.ToRight(xto_right);

     unsigned int max_g_len = G.front().size();
     for(size_t gg=1;gg<G.size();++gg){max_g_len=max(max_g_len,G[gg].size());}

     vec<std::pair<int,int>> permutation(hbp.EdgeObjectCount());
     for(int ii=0;ii<hbp.EdgeObjectCount();++ii){ permutation[ii]=std::make_pair(hbp.EdgeObject(ii).isize(),ii);}
     std::sort(permutation.rbegin(),permutation.rend());


     //very dirty way of load balance, should be coded with a worklist.h instead.
     typedef triple< int, int, pair<int,int> >  og_type;   // the og specification from old code
     typedef std::tuple<og_type,double,int> work_type; // og_type, max_error_rate, offset_add

     vec< vec<work_type> > ee_vec_work( hbp.EdgeObjectCount() ); // edge_idx -> a list of work_type
     vec< std::pair<size_t,size_t> > unit_idx_tt_vec_work; // flattened indices of ee_vec_work

     const int np=3;//number of passes
     #pragma omp parallel
     {
     SmithWatBandedAEngine swbae(sqrt(max_g_len)*2,sqrt(max_g_len));

     #pragma omp for schedule(dynamic,1)
     for ( int ee = 0; ee < hbp.EdgeObjectCount( ); ee++ )
     {
          int i=permutation[ee].second;
          const basevector& e = hbp.EdgeObject(i);
          // For each kmer in the edge, find its hits to the reference and find
          // the kmers having the most hits.
          int nkmers = e.isize( ) - L + 1;
          vec< triple<int64_t,int64_t,int64_t> > locs(nkmers);
          vec<int> pos( nkmers, vec<int>::IDENTITY );
          kmer<L> x;
          for ( int j = 0; j < nkmers; j++ )
          {    x.SetToSubOf( e, j );
               int64_t low = LowerBound(kmers, x), high = UpperBound(kmers, x);
               locs[j].first = high - low;
               locs[j].second = low, locs[j].third = high;    }
          ReverseSortSync( locs, pos );

          // Determine cutoff 'top'.
          double min_cov_frac = 0.5;
          int t = int( floor( nkmers * min_cov_frac ) ), top;
          for ( top = t + 1; top < nkmers; top++ )
               if ( locs[top].first > locs[t].first ) break;

          // Find the associated offsets.
          vec< pair<int,int> > offset;
          for ( int j = 0; j < top; j++ )
          {    for ( int64_t m = locs[j].second; m < locs[j].third; m++ )
               {    int g = kmers_plus[m].second, o = kmers_plus[m].third - pos[j];
                    offset.push( g, o );    }    }
          Sort(offset);

          for(int pass = 0; pass < np; pass++)
          {
//              auto pt = getenv("PASS");
//              if (pt) {
//                  pass = atoi(pt);
//                  np = 1;
//                  cout << "pass= " << pass << endl;
//              }
              RefTraceHeuristics rth;
              switch (pass) {
                  case 0: 
                      //rth.max_offset_diff = 10;    // default
                      //rth.max_error_rate = 0.05;   
                      //rth.offset_add = 1;          // default
                      //rth.max_twiddle = 3;         // default    
                      rth.min_group_frac = 0.1;
                      rth.min_group_save = 200;    
                      break;
                  case 1: 
                      rth.max_offset_diff = 30;
                      rth.max_error_rate = 0.31;
                      rth.offset_add = 5;
                      rth.min_group_frac = 0.1;
                      rth.max_twiddle = 5;    
                      break;
                  case 2: 
                      rth.max_offset_diff = 350;
                      rth.max_error_rate = 0.31;
                      rth.offset_add = 5;
                      rth.min_group_frac = 0.75;
                      rth.max_twiddle = 120;    
                      break;
              }

              // Form offsets into groups.
              vec< triple< int, int, pair<int,int>  > > og;
              for ( int j = 0; j < offset.isize( ); j++ )
              {    int k;
                   for ( k = j + 1; k < offset.isize( ); k++ )
                   {    if ( offset[k].first != offset[j].first ) break;
                        if ( offset[k].second - offset[k-1].second > rth.max_offset_diff )
                             break;    }
                   og.push( k - j, offset[j].first, 
                        make_pair( offset[j].second, offset[k-1].second ) );
                   j = k - 1;    }
              ReverseSort(og);
    
              // Filter offset groups.
              int gj;
              for ( gj = 0; gj < og.isize( ); gj++ )
              {    if ( og[gj].first < rth.min_group_save
                        && og[gj].first < rth.min_group_frac * og[0].first ) 
                   {    break;    }    }
              og.resize(gj);
              for( const auto& entry: og ){
                  ee_vec_work[i].emplace_back( entry , rth.max_error_rate, rth.offset_add);
              }
          }
     }
     {
         #pragma omp barrier
     }
     #pragma omp master
     {
         const size_t n=std::accumulate(ee_vec_work.begin(),ee_vec_work.end(),size_t(0),[](size_t a,vec<work_type>const&b){return a+b.size();});
         unit_idx_tt_vec_work.reserve(n);
         for(const auto& entry: permutation){
             for(size_t ff=0;ff<ee_vec_work[entry.second].size();++ff){
                 unit_idx_tt_vec_work.emplace_back(entry.second,ff);
             }
         }
     }
     {
         #pragma omp barrier
     }
     #pragma omp for schedule(dynamic,1) nowait
     for ( size_t og_idx = 0 ; og_idx < unit_idx_tt_vec_work.size() ; ++og_idx)
     {
         {
              // Align.  The reason for adding to the offset is that there could be in
              // indel in the first or last L bases.
//              for ( int j = 0; j < og.isize( ); j++ )
              {
                   const auto& indices = unit_idx_tt_vec_work[og_idx];
                   const auto& entry = ee_vec_work[indices.first][indices.second];
//                   int g = og[j].second;
//                   int off_low = og[j].third.first, off_high = og[j].third.second;
                   int g = std::get<0>(entry).second;
                   int off_low = std::get<0>(entry).third.first, off_high = std::get<0>(entry).third.second;
                   int mid_offset = ( off_low + off_high ) / 2;
                   int bandwidth 
                        = Max(mid_offset - off_low, off_high - mid_offset) + std::get<2>(entry);// rth.offset_add;
                   // Do the alignment.  This is kludgy.  If the alignment has too 
                   // many errors and the edge is long, we suspect that the problem
                   // might be with a big indel, so we align using a larger bandwidth.
                   // Note the unfortunate us of hardcoded constants.
                   align a;
                   int errors;
                   swbae.run( hbp.EdgeObject(indices.first), G[g],
                        -mid_offset, bandwidth, a, errors, 0, 1, 1 );
                   if ( double(errors) / double( a.extent2( ) ) > std::get<1>(entry) /*rth.max_error_rate*/ ) {
                        const int long_edge = 5000;
                        const int max_indel = 5000;
                        if ( hbp.EdgeLengthBases(indices.first) < long_edge ) continue;
                        swbae.run( hbp.EdgeObject(indices.first), G[g],
                             -mid_offset, max_indel, a, errors, 0, 1, 1 );
                        if ( double(errors) / double( a.extent2( ) ) > std::get<1>(entry)/*rth.max_error_rate*/ )
                             continue;    
                   }

//                   errors += a.pos1();
//                   errors += hbp.EdgeObject(i).size()-a.Pos1();

                   // Figure out where the position e.isize( ) - K + 1 should map to
                   // under the alignment.  Note that because there could be an indel
                   // there, this is not necessarily a meaningful answer.
                   int x1 = hbp.EdgeObject(indices.first).isize( ) - hbp.K( ) + 1;
                   int x2 = CorrelatePositionsAlways( a, x1 );
                   #pragma omp critical
                   {    vedata.push( make_triple( g, a.pos2( ), xto_left[indices.first] ),
                             make_triple( g, x2, xto_right[indices.first] ),
                             make_pair( indices.first, errors ), make_pair( a.pos1( ), a.Pos1( ) ) );
                        aligns.push_back(a);    }
              }    
          }
     }    

     }//omp parallel
     // Sort the output to avoid the stochastic downstream behavior of BuildGraph
     // that seems depend on the input order of the alignment data.
     UniqueSortSync(vedata, aligns);
}


#endif

