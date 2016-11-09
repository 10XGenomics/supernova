///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "CoreTools.h"
#include "PairsManager.h"
#include "ParallelVecUtilities.h"
#include "Qualvector.h"
#include "efasta/EfastaTools.h"
#include "math/HoInterval.h"
#include "math/IntDistribution.h"
#include "paths/HyperBasevector.h"
#include "paths/KmerBaseBroker.h"
#include "paths/ReadsToPathsCoreX.h"
#include "paths/UnibaseUtils.h"
#include "paths/Unipath.h"
#include "paths/long/CreateGenome.h"
#include "paths/long/EvalByReads.h"
#include "paths/long/LongHyper.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/LongReadsToPaths.h"
#include "paths/long/PairInfo.h"
#include "paths/long/SupportedHyperBasevector.h"

Bool LongHyper( const String& READS, const VecEFasta& correctede, 
     const vec<pairing_info>& cpartner, SupportedHyperBasevector& shb, 
     const long_heuristics& heur,
     const long_logging_control& log_control, const long_logging& logc, 
     const LongProtoTmpDirManager& tmp_mgr, bool useOldLRPMethod )
{
     // Choose K2.

     double K2frac_mult = 1.0;
     if ( READS != "" && !READS.Contains( ".fastb", -1 ) && heur.CORRECT_PAIRS ) 
          K2frac_mult = 0.82;
     double K2frac_to_use = heur.K2frac * K2frac_mult;
     int K2 = SelectK2( correctede, K2frac_to_use, logc, heur );
     K2 = Max( K2, heur.K2_FLOOR );
     if ( heur.K2_FORCE >= 0 ) K2 = heur.K2_FORCE;
     double fudge_mult = K2frac_mult;

     // Setup for PERTURB_TRANSLATIONS.

     const bool bUseOrgReads= ( heur.PERTURB_TRANSLATIONS || logc.TRACE_EDGES0 != "" || heur.REMAP_READS );
     vecbasevector const& Bases = bUseOrgReads? tmp_mgr.get("frag_reads_orig").reads() : vecbasevector();
     vecqualvector const& Quals = bUseOrgReads? tmp_mgr.get("frag_reads_orig").quals() : vecqualvector();
     if(bUseOrgReads){ tmp_mgr.get("frag_reads_orig").pairs().makeCache(); }
     PairsManager const& Pairs = bUseOrgReads? tmp_mgr.get("frag_reads_orig").pairs() : PairsManager();

     // Expand correctede.  Note that this is exponential and thus totally unsound.

     double eclock = WallClockTime( );
     if (logc.STATUS_LOGGING) ReportPeakMem( );
     if (logc.STATUS_LOGGING) cout << Date( ) << ": making hyper" << endl;
     vecbasevector correctedv;
     vec< triple<int,int,int> > origin;
     for ( size_t id = 0; id < correctede.size( ); id++ )
     {    vec<basevector> b;
          correctede[id].ExpandTo(b);
          for ( int j = 0; j < b.isize( ); j++ )
          {    correctedv.push_back_reserve( b[j] );
               origin.push( id, j, b.size( ) );    }    }
     REPORT_TIME( eclock, "used in expansion" );
     int64_t read_count = correctedv.size( );
     if ( heur.INJECT_REF && log_control.G != 0 ) 
          correctedv.Append( *(log_control.G) );

     // Make paths.

     double pclock = WallClockTime( );
     HyperBasevector hb;
     HyperKmerPath h;
     vecKmerPath paths, paths_rc;
     if (logc.STATUS_LOGGING)
          cout << Date() << ": calling LongReadsToPaths" << endl;
     unsigned const COVERAGE = 50u;
     LongReadsToPaths(correctedv,K2,COVERAGE,logc.verb[ "LRP" ],useOldLRPMethod,
                          &hb,&h,&paths,&paths_rc);
     REPORT_TIME( pclock, "used in pathing" );

     // Trace reads through h.  Ignore reads that lie entirely on one edge.
     // Note that in fact each bona fide read can give rise to several 'reads'
     // via EFASTA and that each of these are traced here.  Note also that we 
     // discard reads that do not have a perfect mapping to the graph, and this 
     // could have bad consequences.

     double tclock = WallClockTime( );
     vec< vec<int> > usu;
     vec<fix64_6> count_fw, count_rc;
     vec< pair< vec<int>, vec<int> > > pairs;
     vec< vec<pair_point> > pair_data;
     vec< vec< pair<fix64_6,int64_t> > > weights_fw_origin, weights_rc_origin;
     {    if (logc.STATUS_LOGGING) cout << Date( ) << ": start tracing reads" << endl;
          vecKmerPath hpaths;
          vec<big_tagged_rpint> hpathsdb;
          for ( int e = 0; e < h.EdgeObjectCount( ); e++ )
               hpaths.push_back_reserve( h.EdgeObject(e) );
          CreateDatabase( hpaths, hpathsdb );

          // Track in parallel:
          // us = sequences in unipath graph
          // weight_fw = weight of sequence (forward orientation)
          // weight_rc = weight of sequence (reverse orientation)
          // denom_fw = reciprocal of forward weight (up to rounding)
          // denom_rc = reciprocal of reverse weight (up to rounding)
          // left = left extension of read by sequence
          // right = right extension of read by sequence
          // ids = original read id
          // fwplace = True if forward orientation

          vec< vec<int> > us;
          vec<fix64_6> weight_fw, weight_rc;
          vec<int> denom_fw, denom_rc;
          vec<int> left, right;
          vec<int> ids;
          vec<Bool> fwplace;

          // Go through the reads.

          vec<int> starts;
          for ( int id = 0; id < read_count; id++ )
          {    starts.push_back( us.size( ) );
               for ( int pass = 1; pass <= 2; pass++ )
               {    const KmerPath& p = ( pass == 1 ? paths[id] : paths_rc[id] );
                    vec< triple<ho_interval,int,ho_interval> > M, M2;
                    vec<int> u;
                    int rpos = 0;
                    for ( int j = 0; j < p.NSegments( ); j++ )
                    {    const KmerPathInterval& I = p.Segment(j);
                         vec<longlong> locs;
                         Contains( hpathsdb, I, locs );
                         for ( int l = 0; l < locs.isize( ); l++ )
                         {    const big_tagged_rpint& t = hpathsdb[ locs[l] ];
                              int hid = t.PathId( );
                              if ( hid < 0 ) continue;
                              longlong hpos = I.Start( ) - t.Start( );
                              longlong start = Max( I.Start( ), t.Start( ) );
                              longlong stop = Min( I.Stop( ), t.Stop( ) );
                              longlong hstart = start - t.Start( );
                              for ( int r = 0; r < t.PathPos( ); r++ )
                                   hstart += hpaths[hid].Segment(r).Length( );
                              longlong hstop = hstart + stop - start;
                              longlong rstart = rpos + start - I.Start( );
                              longlong rstop = rstart + stop - start;
                              M.push( ho_interval( rstart, rstop ), hid,
                                   ho_interval( hstart, hstop ) );    }
                         rpos += I.Length( );    }    
                    Bool bad = False;
                    for ( int i = 0; i < M.isize( ); i++ )
                    {    int j;
                         for ( j = i + 1; j < M.isize( ); j++ )
                         {    if ( M[j].first.Start( ) != M[j-1].first.Stop( ) + 1 ) 
                                   break;
                              if ( M[j].second != M[j-1].second ) break;
                              if ( M[j].third.Start( ) != M[j-1].third.Stop( ) + 1 ) 
                                   break;    }
                         u.push_back( M[i].second );
                         Bool incomplete = False;
                         if ( i > 0 && M[i].third.Start( ) > 0 ) incomplete = True;
                         if ( j < M.isize( ) && M[j-1].third.Stop( )
                              != hpaths[ M[i].second ].KmerCount( ) - 1 )
                         {    incomplete = True;
                              bad = True;    }
                         if ( logc.TRACE_READS && i == 0 )
                         {    cout << "\n" << ( pass == 1 ? "+" : "-" ) 
                                   << origin[id].first << "."
                                   << origin[id].second + 1 << "\n";    }
                         if ( i == 0 && j == M.isize( ) && !incomplete )
                         {    i = j - 1;
                              continue;    }
                         int last = ( i == 0 ? -1 : M2.back( ).first.Stop( ) );
                         if ( M[i].first.Start( ) > last + 1 )
                         {    if (logc.TRACE_READS)
                              {    cout << last+1 << "-" << M[i].first.Start( ) - 1
                                        << " --> MISSING\n";    }
                              bad = True;    }
                         M2.push( ho_interval( M[i].first.Start( ),
                              M[j-1].first.Stop( ) ), M[i].second,
                              ho_interval( M[i].third.Start( ),
                              M[j-1].third.Stop( ) ) );
                         if (logc.TRACE_READS)
                         {    cout << M[i].first.Start( ) << "-" 
                                   << M[j-1].first.Stop( ) << " --> " 
                                   << M[i].second << "." << M[i].third.Start( )
                                   << "-" << M[j-1].third.Stop( )
                                   << ( incomplete ? " [INCOMPLETE]" : "" ) 
                                   << "\n";    }
                         if ( j == M.isize( ) && M[j-1].first.Stop( ) 
                              < p.KmerCount( ) - 1 )
                         {    if (logc.TRACE_READS)
                              {    cout << M[j-1].first.Stop( ) + 1 << "-"
                                        << p.KmerCount( ) - 1 
                                        << " --> MISSING\n";    }
                              bad = True;    }
                         i = j - 1;    }
                    if ( !bad && u.nonempty( ) )
                    {    if (logc.TRACE_READS)
                         {    cout << "u =";
                              for ( int j = 0; j < u.isize( ); j++ )
                                   cout << " " << u[j];
                              cout << "\n";    }
                         us.push_back(u);    
                         if ( pass == 1 )
                         {    weight_fw.push_back( fix64_6( 1, origin[id].third ) );
                              weight_rc.push_back(0);
                              denom_fw.push_back( origin[id].third );
                              denom_rc.push_back( 0 );    }
                         else
                         {    weight_rc.push_back( fix64_6( 1, origin[id].third ) );
                              weight_fw.push_back(0);
                              denom_rc.push_back( origin[id].third );
                              denom_fw.push_back( 0 );    }
                         left.push_back( M.front( ).third.Start( ) );
                         right.push_back( hb.EdgeObject( u.back( ) ).isize( ) 
                              - ( M.back( ).third.Stop( ) + K2 ) );    
                         ids.push_back( origin[id].first );    
                         fwplace.push_back( pass == 1 );    }    }    }
          starts.push_back( us.size( ) );

          // Remap reads.  Perhaps not compatible with PERTURB_TRANSLATIONS.

          if ( heur.REMAP_READS )
          {    double clock = WallClockTime( );
               cout << Date( ) << ": remapping reads" << endl;
               const int infinity = 1000000000;
               const int L = 12;
               HyperBasevector hb_fw(hb), hb_rc(hb);
               hb_rc.Reverse( );
               vec<int> to_left_fw, to_left_rc;
               hb_fw.ToLeft(to_left_fw), hb_rc.ToLeft(to_left_rc);
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
               vec< vec<int> > ids_index( Bases.size( ) );
               for ( int64_t i = 0; i < ids.jsize( ); i++ )
                    ids_index[ ids[i] ].push_back(i);
               vec<int> LL;
               for ( int j = 0; j < hb.EdgeObjectCount( ); j++ )
                    LL.push_back( hb.EdgeLengthKmers(j) );
               digraphE<int> G( hb, LL );
               #pragma omp parallel for
               for ( int64_t pid = 0; pid < (int64_t) Pairs.nPairs( ); pid++ )
               {    int64_t id1 = Pairs.ID1(pid), id2 = Pairs.ID2(pid);
                    vec< vec<read_place> > places(2);
                    for ( int pass = 0; pass < 2; pass++ )
                    {    int64_t id = ( pass == 0 ? id1 : id2 );
                         int n = KmerId( Bases[id], L, 0 );
                         int qual_sum = infinity;
                         FindPlaces( Bases[id], Quals[id], n, hb_fw, hb_rc, 
                              to_right_fw, to_right_rc, locs_fw, locs_rc, 
                              places[pass], qual_sum );    }
                    if ( places[0].size( ) != 2 || places[1].size( ) != 2 ) continue;
                    {    vec< vec< vec<int> > > joins(2);
                         for ( int pass = 0; pass < 2; pass++ )
                         {    int m = ( pass == 0 ? 0 : 1 );
                              for ( int r1 = 0; r1 < 2; r1++ )
                              for ( int r2 = 0; r2 < 2; r2++ )
                              {    if ( !places[m][r1].Fw( ) ) continue;
                                   if ( places[1-m][r2].Fw( ) ) continue;
                                   vec<int> v1 = places[m][r1].E( );
                                   vec<int> v2 = places[1-m][r2].E( );
                                   v2.ReverseMe( );
                                   for ( int r1 = 0; r1 < v1.isize( ); r1++ )
                                   {    Bool match = True;
                                        for ( int r2 = 0; r2 < v2.isize( ); r2++ )
                                        {    if ( r1 + r2 == v1.isize( ) ) break;
                                             if ( v1[r1+r2] != v2[r2] )
                                             {    match = False;
                                                  break;    }    }
                                        if ( !match ) continue;
                                        vec<int> join;
                                        for ( int k = 0; k < r1; k++ )
                                             join.push_back( v1[k] );
                                        join.append(v2);
                                        joins[pass].push_back(join);    }
                                   const int max_sep = 500; // should not be hardcoded!
                                   int v = to_right_fw[ v1.back( ) ];
                                   int w = to_left_fw[ v2.front( ) ];
                                   vec< vec<int> > paths;
                                   const int max_paths = 5;
                                   Bool ok = G.AllPathsLengthRange( v, w, 0, max_sep, 
                                        to_right_fw, paths, max_paths );
                                   for ( int l = 0; l < paths.isize( ); l++ )
                                   {    vec<int> x(v1);
                                        x.append( paths[l] );
                                        x.append(v2);
                                        joins[pass].push_back(x);    }    }    }
                         if ( !joins[0].solo( ) || !joins[1].solo( ) ) continue;
                         Bool owned = False;
                         for ( int i = 0; i < ids_index[id1].isize( ); i++ )
                         {    int m = ids_index[id1][i];
                              if ( us[m] == joins[0][0] ) owned = True;    }
                         for ( int i = 0; i < ids_index[id2].isize( ); i++ )
                         {    int m = ids_index[id2][i];
                              if ( us[m] == joins[1][0] ) owned = True;    }
                         if (owned) continue;

                         Bool remap_verbose = False;
                         if (remap_verbose)
                         {
                         #pragma omp critical
                         {    cout << "\n";
                              PRINT2( id1, id2 );
                              for ( int pass = 0; pass < 2; pass++ )
                              {    for ( int l = 0; l < joins[pass].isize( ); l++ )
                                   {    cout << "new " << pass+1 << " " 
                                             << printSeq( joins[pass][l] ) 
                                             << endl;    }     }
                              for ( int i = 0; i < ids_index[id1].isize( ); i++ )
                              {    int m = ids_index[id1][i];
                                   cout << "old 1 --> " 
                                        << ( fwplace[m] ? "fw" : "rc" ) << " "
                                        << printSeq( us[m] ) << endl;    }
                              for ( int i = 0; i < ids_index[id2].isize( ); i++ )
                              {    int m = ids_index[id2][i];
                                   cout << "old 2 --> " 
                                        << ( fwplace[m] ? "fw" : "rc" ) << " "
                                        << printSeq( us[m] ) << endl;    }    }    
                         }

                         // Let's first deal with the case where there was no old
                         // placement.

                         if ( ids_index[id1].empty( ) && ids_index[id2].empty( ) )
                         {
                              #pragma omp critical
                              {    us.push_back( joins[0][0] );
                                   weight_fw.push_back(1), weight_rc.push_back(0);
                                   denom_fw.push_back(1), denom_rc.push_back(0);
                                   fwplace.push_back(True);
                                   ids.push_back(id1);
                                   int L, R;
                                   if ( places[0][0].Fw( ) ) L = places[0][0].P( );
                                   else L = places[0][1].P( );
                                   if ( !places[0][0].Fw( ) ) R = places[0][0].P( );
                                   else R = places[0][1].P( );
                                   left.push_back(L), right.push_back(R);   
                                   us.push_back( joins[1][0] );
                                   weight_fw.push_back(0), weight_rc.push_back(1);
                                   denom_fw.push_back(0), denom_rc.push_back(1);
                                   fwplace.push_back(False);
                                   ids.push_back(id1);
                                   left.push_back(R);
                                   right.push_back(L);    }    }    }    }
               if (logc.STATUS_LOGGING)
                    cout << "\n" << Date( ) << ": done" << endl;    }

          // Look up pairs.

          for ( int id1 = 0; id1 < cpartner.isize( ); id1++ )
          {    const pairing_info& p = cpartner[id1];
               if ( p.Status( ) != 1 ) continue;
               // cout << "\npair\n";
               int id2 = p.Partner( ), lib = p.LibId( );
               for ( int m1 = starts[id1]; m1 < starts[id1+1]; m1++ )
               for ( int m2 = starts[id2]; m2 < starts[id2+1]; m2++ )
               {    if ( fwplace[m1] == fwplace[m2] ) continue;
                    int n1(m1), n2(m2);
                    if ( !fwplace[m1] ) swap( n1, n2 );
                    fix64_6 w( 1, ( denom_fw[n1] + denom_rc[n1] ) 
                         * ( denom_fw[n2] + denom_rc[n2] ) );
                    int trim = right[n1] + left[n2];
                    if (logc.PRINT_INITIAL_PAIRS)
                    {    cout << printSeq( us[n1] ) << " [w=" << weight_fw[n1] 
                              << "+" << weight_rc[n1] << ",trim=" << right[n1] << "]"
                              << " --" << lib << "--> " << printSeq( us[n2] )
                              << " [w=" << weight_fw[n2] << "+" << weight_rc[n2]
                              << ",trim=" << left[n2] << "]"
                              << " (trim=" << trim << ", weight=" << w << ")\n";    }
                    pairs.push( us[n1], us[n2] );
                    vec<pair_point> p;
                    p.push( trim, w, lib );
                    pair_data.push(p);    }    }
          SortSync( pairs, pair_data );
          vec< pair< vec<int>, vec<int> > > pairs2;
          vec< vec<pair_point> > pair_data2;
          for ( int i = 0; i < pairs.isize( ); i++ )
          {    int j = pairs.NextDiff(i);
               pairs2.push_back( pairs[i] ); 
               vec<pair_point> p;
               for ( int k = i; k < j; k++ )
                    p.push_back( pair_data[k][0] );
               Sort(p);
               pair_data2.push_back(p);
               i = j - 1;    }
          pairs = pairs2;
          pair_data = pair_data2;

          // Experimental code to attempt to perturb read translations.

          if ( heur.PERTURB_TRANSLATIONS && Bases.size( ) > 0 )
          {    
               Bool verbose = logc.verb[ "PERTURB_TRANSLATIONS" ] >= 1;
               vec<Bool> udel( us.size( ), False );

               const int L = 12;
               vecbasevector x;
               for ( int i = 0; i < hb.EdgeObjectCount( ); i++ )
                    x.push_back( hb.EdgeObject(i) );
               VecIntPairVec locs;
               CreateGlocs( x, L, locs );

               vec<int> inv;
               vecbasevector hbx;
               for ( int e = 0; e < hb.EdgeObjectCount( ); e++ )
                    hbx.push_back_reserve( hb.EdgeObject(e) );
               UnibaseInvolution( hbx, inv );
     
               vec<int> to_left, to_right;
               hb.ToLeft(to_left), hb.ToRight(to_right);
     
               vec< vec<int> > usu(us);
               ParallelUniqueSort(usu);
               vec< vec< pair<int,int> > > usu_index( hb.EdgeObjectCount( ) );
               for ( int i = 0; i < usu.isize( ); i++ )
               for ( int j = 0; j < usu[i].isize( ); j++ )
                    usu_index[ usu[i][j] ].push( i, j );
               int nus = us.size( );

               // Make a copy of us.  This is a workaround for what seems to be
               // a bug in OpenMP (no it is not).  Otherwise, 'target' (below) 
               // would sometimes
               // end up with garbage in it, and this is presumably associated with
               // the pushback onto us and presumed resizing of that.  Interesting.

               vec< vec<int> > us_copy(us);

               #pragma omp parallel for
               for ( int i = 0; i < nus; i++ )
               {    
                    vec<int> target = us_copy[i];
                    if ( !fwplace[i] )
                    {    target.ReverseMe( );
                         for ( int j = 0; j < target.isize( ); j++ )
                              target[j] = inv[ target[j] ];    }

                    int64_t id1 = ids[i], id2 = Pairs.getPartnerID( ids[i] );

                    int n1 = KmerId( Bases[id1], L, 0 );
                    int n2 = KmerId( Bases[id2], L, 0 );

                    vec< triple<int,ho_interval,ho_interval> > segs1x, segs2x;
                    vec<int> val1x, val2x;
     
                    vec< triple< vec<int>, int, int > > where;
                    const int64_t max_locs = 1000;
                    if ( 1l*locs[n1].size()*locs[n2].size( ) > max_locs ) continue;
                    for ( unsigned j1 = 0; j1 < locs[n1].size( ); j1++ )
                    for ( unsigned j2 = 0; j2 < locs[n2].size( ); j2++ )
                    {    int u1 = locs[n1][j1].first;
                         int u2 = inv[ locs[n2][j2].first ];
                         int p1 = locs[n1][j1].second, p2 = locs[n2][j2].second;
                         p2 = hb.EdgeObject(u2).isize( ) - p2 - L;
                         for ( int l1 = 0; l1 < usu_index[u1].isize( ); l1++ )
                         for ( int l2 = 0; l2 < usu_index[u2].isize( ); l2++ )
                         {    int pid = usu_index[u1][l1].first;
                              if ( usu_index[u2][l2].first != pid ) continue;
                              if ( usu_index[u1][l1].second > 
                                   usu_index[u2][l2].second ) 
                              {    continue;    }
                              vec<int> w;
                              for ( int m = usu_index[u1][l1].second;
                                   m <= usu_index[u2][l2].second; m++ )
                              {    w.push_back( usu[pid][m] );    }
                              where.push( w, p1, p2 + L );    }    }
                    UniqueSort(where);
                    vec<int> Q( where.size( ), -1 );
                    for ( int j = 0; j < where.isize( ); j++ )
                    {    int u1 = where[j].first.front( ); 
                         int un = where[j].first.back( );
                         int p1 = where[j].second, p2 = where[j].third;
                         if ( hb.To( to_left[un] ).nonempty( ) && p2 < K2 ) continue;
                         if ( hb.From( to_right[u1] ).nonempty( ) 
                              && p1 > hb.EdgeObject(u1).isize( ) - K2 )
                         {    continue;    }

                         // Map out where the first read goes.

                         vec< triple<int,ho_interval,ho_interval> > segs1;
                         Bool fail = False;
                         {    int u = u1, p = p1, i = 0, l = 0;
                              while(1)
                              {    int dist = Min( Bases[id1].isize( ) - l,
                                        hb.EdgeObject(u).isize( ) - p );
                                   int q = p + dist;
                                   segs1.push( u, ho_interval(p,q),
                                        ho_interval( l, l + dist ) );
                                   i++;
                                   l += dist;
                                   if ( i == where[j].first.isize( ) )
                                   {    if ( q > p2 || l != Bases[id1].isize( ) ) 
                                             fail = True;
                                        break;   }
                                   if ( l == Bases[id1].isize( ) ) break;
                                   p = K2 - 1;
                                   u = where[j].first[i];    }    }
                         if (fail) continue;

                         // Compute the corresponding qsum.

                         int qsum = 0;
                         for ( int s = 0; s < segs1.isize( ); s++ )
                         {    Bool found = False;
                              for ( int m = 0; m < segs1x.isize( ); m++ )
                              {    if ( segs1[s] == segs1x[m] )
                                   {    found = True;
                                        qsum += val1x[m];
                                        break;    }    }
                              if (found) continue;
                              int u = segs1[s].first;
                              int q = 0;
                              for ( int l = segs1[s].third.Start( );
                                   l < segs1[s].third.Stop( ); l++ )
                              {    int p = segs1[s].second.Start( )
                                        + l - segs1[s].third.Start( );
                                   if ( Bases[id1][l] != hb.EdgeObject(u)[p] )
                                        q += Quals[id1][l];    }
                              qsum += q;
                              segs1x.push_back( segs1[s] );
                              val1x.push_back(q);    }

                         // Map out where the second read goes.

                         vec< triple<int,ho_interval,ho_interval> > segs2;
                         {    int u = un, p = p2; 
                              int i = where[j].first.isize( ) - 1, l = 0;
                              while(1)
                              {    int dist = Min( Bases[id2].isize( ) - l, p );
                                   int q = p - dist;
                                   segs2.push( u, ho_interval(q,p),
                                        ho_interval( l, l + dist ) );
                                   i--;
                                   l += dist;
                                   if ( i < 0 )
                                   {    if ( q < p1 || l != Bases[id2].isize( ) ) 
                                             fail = True;
                                        break;   }
                                   if ( l == Bases[id2].isize( ) ) break;
                                   u = where[j].first[i];
                                   p = hb.EdgeObject(u).isize( ) - K2 + 1;    }    }
                         if (fail) continue;

                         // Compute the corresponding qsum.

                         for ( int s = 0; s < segs2.isize( ); s++ )
                         {    Bool found = False;
                              for ( int m = 0; m < segs2x.isize( ); m++ )
                              {    if ( segs2[s] == segs2x[m] )
                                   {    found = True;
                                        qsum += val2x[m];
                                        break;    }    }
                              if (found) continue;
                              int u = segs2[s].first;
                              int q = 0;
                              for ( int l = segs2[s].third.Start( );
                                   l < segs2[s].third.Stop( ); l++ )
                              {    int p = segs2[s].second.Stop( ) - 1
                                        - ( l - segs2[s].third.Start( ) );
                                   if ( 3 - Bases[id2][l] != hb.EdgeObject(u)[p] )
                                        q += Quals[id2][l];    }
                              qsum += q;
                              segs2x.push_back( segs2[s] );
                              val2x.push_back(q);    }
                         Q[j] = qsum;    }
     
                    SortSync( Q, where );
                    vec<Bool> to_delete( Q.size( ), False );
                    for ( int j = 0; j < Q.isize( ); j++ )
                         if ( Q[j] < 0 ) to_delete[j] = True;
                    EraseIf( Q, to_delete );
                    EraseIf( where, to_delete );
                    const int max_diff = 10;
                    for ( int j = 1; j < Q.isize( ); j++ )
                    {    if ( Q[j] - Q[0] > max_diff )
                         {    Q.resize(j);
                              where.resize(j);
                              break;    }    }
     
                    if ( where.empty( ) ) continue;
                    fix64_6 weight = weight_fw[i] + weight_rc[i];
                    if ( weight < 1 ) continue;
                    if ( where.size( ) > 1 || where[0].first != target )
                    {    vec<long double> ps;
                         for ( int j = 0; j < where.isize( ); j++ )
                              ps.push_back( exp10l( -Q[j]/10.0 ) );
                         Sort(ps);
                         long double psum = Sum(ps);

                         if (verbose)
                         {    
                              #pragma omp critical
                              {    cout << "\nlooking at pair " << id1 << "/" << id2 
                                        << endl;
                                   cout << "weight = " << weight << endl;
                                   PRINT( int(fwplace[i]) );
                                   cout << "translates to " << printSeq(target) 
                                        << endl;
                                   cout << "start1:\n";
                                   for ( unsigned j1 = 0; j1 < locs[n1].size( ); j1++ )
                                   {    int u1 = locs[n1][j1].first; 
                                        int p1 = locs[n1][j1].second;
                                        cout << u1 << "." << p1 << "\n";    }
                                   cout << "start2:\n";
                                   for ( unsigned j2 = 0; j2 < locs[n2].size( ); j2++ )
                                   {    int u2 = inv[ locs[n2][j2].first ]; 
                                        int p2 = locs[n2][j2].second;
                                        p2 = hb.EdgeObject(u2).isize( ) - p2 - L;
                                        cout << u2 << "." << p2 << "\n";    }
                                   cout << "locs:\n";
                                   for ( int j = 0; j < where.isize( ); j++ )
                                   {    int u1 = where[j].first.front( ); 
                                        int un = where[j].first.back( );
                                        int p1 = where[j].second; 
                                        int p2 = where[j].third;
                                        cout << u1 << "." << p1;
                                        for ( int l = 1; l < where[j].first.isize( );
                                             l++ )
                                        {    cout << "," << where[j].first[l];    }
                                        long double p = exp10l( -Q[j]/10.0 ) / psum;
                                        fix64_6 weight( p * 1000000, 1000000 );
                                        cout << "." << where[j].third 
                                             << " [Q=" << Q[j] << ",w=" << weight 
                                             << "]";
                                        if ( where[j].first == target )
                                             cout << " ***";
                                   cout << "\n";    }    }    }
                              
                         #pragma omp critical
                         {    udel[i] = True;
                              for ( int j = 0; j < where.isize( ); j++ )
                              {    long double p = exp10l( -Q[j]/10.0 ) / psum;
                                   fix64_6 weight( p * 1000000, 1000000 );
                                   if ( !fwplace[i] )
                                   {    where[j].first.ReverseMe( );
                                        for ( int l = 0; l < where[j].first.isize( );
                                             l++ )
                                        {    where[j].first[l] = inv[ 
                                                  where[j].first[l] ];    }    }
                                   us.push_back( where[j].first );
                                   ids.push_back(id1);
                                   udel.push_back(False);
                                   if ( weight_fw[i] > 0 )
                                   {    weight_fw.push_back( weight/2 );
                                        weight_rc.push_back(0);    }
                                   else
                                   {    weight_rc.push_back( weight/2 );
                                        weight_fw.push_back(0);    
                                        }    }    }    }    }
               EraseIf( us, udel );
               EraseIf( weight_fw, udel );
               EraseIf( weight_rc, udel );
               EraseIf( ids, udel );    }

          // Generate main words.

          SortSync( us, weight_fw, weight_rc, ids );    
          for ( int i = 0; i < us.isize( ); i++ )
          {    int j = us.NextDiff(i);
               usu.push_back( us[i] );
               fix64_6 cfw = 0, crc = 0;
               for ( int k = i; k < j; k++ )
               {    cfw += weight_fw[k];
                    crc += weight_rc[k];    }
               count_fw.push_back(cfw);
               count_rc.push_back(crc);
               vec< pair<fix64_6,int64_t> > wfw, wrc;
               for ( int k = i; k < j; k++ )
               {    if ( weight_fw[k] > 0 ) wfw.push( weight_fw[k], ids[k] );
                    if ( weight_rc[k] > 0 ) wrc.push( weight_rc[k], ids[k] );    }
               weights_fw_origin.push_back(wfw);
               weights_rc_origin.push_back(wrc);
               i = j - 1;    }
          if (logc.TRACE_READS) 
          {    cout << "\ntraces:\n";
               for ( int i = 0; i < usu.isize( ); i++ )
               {    cout << "[" << i+1 << "," << setiosflags(ios::fixed) 
                         << setprecision(1) << count_fw[i] + count_rc[i]
                         << resetiosflags(ios::fixed) << "x]";
                    for ( int j = 0; j < usu[i].isize( ); j++ )
                         cout << " " << usu[i][j];
                    cout << "\n";    }    
               cout << "\n";    }
          if (logc.STATUS_LOGGING)
               cout << Date( ) << ": read tracing complete" << endl;    }

     // Get median read length and create SupportedHyperBasevector.

     vec<int> len;
     for ( size_t i = 0; i < correctede.size( ); i++ )
     {    int n = correctede[i].Length1( );
          if ( n == 0 ) continue;
          len.push_back(n);    }
     Sort(len);
     IntDistribution read_length_dist;
     read_length_dist.from_hits(len);
     vec<int> inv;
     vecbasevector hbx;
     for ( int e = 0; e < hb.EdgeObjectCount( ); e++ )
          hbx.push_back_reserve( hb.EdgeObject(e) );
     REPORT_TIME( tclock, "used in tracing" );
     double iclock = WallClockTime( );
     UnibaseInvolution( hbx, inv );
     REPORT_TIME( iclock, "used finding involution" );
     double z1clock = WallClockTime( );
     if ( usu.empty( ) ) return False;
     shb = SupportedHyperBasevector( hb, inv, usu, count_fw, count_rc, 
          weights_fw_origin, weights_rc_origin, pairs, pair_data, 
          correctede.size( ), read_length_dist, fudge_mult );
     REPORT_TIME( z1clock, "used after tracing" );
     shb.FixWeights(logc);
     shb.TestValid(logc);
     shb.DumpFilesStandard( log_control, logc, 1 );
     if (logc.STATUS_LOGGING) ReportPeakMem( );

     // Trace edges.

     if ( logc.TRACE_EDGES0 != "" )
          TraceEdges( shb, logc.TRACE_EDGES0, Bases, Quals );

     return True;    }
