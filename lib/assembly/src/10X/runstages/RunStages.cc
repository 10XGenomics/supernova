// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

#include "10X/runstages/RunStages.h"
#include "feudal/SubsetMasterVec.h"
#include "10X/Super.h"
#include "10X/Rescue.h"
#include "10X/PathsIndex.h"
#include "10X/ClosuresToGraph.h"
#include "10X/paths/ReadPathVecX.h"

void StageClosures(HyperBasevectorX& hb, vec<int>& inv, ReadPathVecX& paths,
          String paths_index_file, vec<Bool>& dup, vec<Bool>& bad,
          vec<vec<int>>& all_closures, digraphE<vec<int>>& D,
          vec<int>& dinv, String fin_dir)
{
     STAGE(Closures);
     // Make closures and the graph from them.
     cout << Date( ) << ": start making closures" << endl;
     
     VecULongVec paths_index;
     MakeClosures( hb, inv, paths, paths_index, dup, bad, all_closures, True, paths_index_file );
     MEM(before_closures_to_graph);
     ClosuresToGraph( hb, inv, all_closures, D, dinv, True, fin_dir );
     Validate( hb, inv, D, dinv );
     MEM(after_closures_to_graph);
}


void StageEBC(const ReadPathVecX & paths,const HyperBasevectorX& hb, 
        const vec <int32_t> & bc, const vec<int> & inv, 
        const vec<int64_t> & bci, VecIntVec & ebcx )
{
     STAGE(EBC);
     computeEdgeToBarcodeX(paths, hb, bc, inv, bci, ebcx, false );

}

void StageDups( HyperBasevectorX& hb, vecbasevector& bases,
          ObjectManager<VecPQVec>& quals_om, ReadPathVecX& paths,
          vec<int32_t>& bc, vec<Bool>& dup, vec<Bool>& bad, double& interdup_rate)

{
     STAGE(Dups);
//     VecPQVec& quals = quals_om.load_mutable( );
     quals_om.unload();
     VirtualMasterVec<PQVec> vmv( quals_om.filename() );
//     MarkDups( bases, quals, paths, bc, dup, interdup_rate );
     MarkDups( bases, vmv, paths, hb, bc, dup, interdup_rate );
//     MarkBads( hb, bases, quals, paths, bad );
     MarkBads( hb, bases, vmv, paths, bad );
//     quals_om.unload( );
}

void StageTrim( HyperBasevectorX& hb, HyperBasevector& hbv, vec<int>& inv,
          vecbasevector& bases, ObjectManager<VecPQVec>& quals_om,
          ReadPathVecX& paths, String paths_index_file, String countsb_file, 
          vec<int32_t>& bc)
{
     STAGE(Trim);
     quals_om.unload();
     MEM(trim_unload_quals_om);
     VirtualMasterVec<PQVec> vmv( quals_om.filename() );
     vec<int> dels;
     vec<Bool> dup;
     double interdup;
     MarkDups( bases, vmv, paths, hb, bc, dup, interdup );
//     VecPQVec& quals = quals_om.load_mutable( );
     MowLawn( hb, inv, bases, quals_om, paths, paths_index_file, bc, dup, interdup, dels );
     quals_om.unload();

     // NOTE EXPENSIVE CONVERSION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     cout << Date( ) << ": hbx -> hbv" << endl;
     hbv = HyperBasevector(hb);
     cout << Date( ) << ": delete edges" << endl;
     hbv.DeleteEdgesParallel(dels);
     cout << Date( ) << ": pathsX -> paths" << endl;
     ReadPathVec rpaths;
     paths.parallel_unzip(rpaths,hb);
     cout << Date( ) << ": destroy pathsX" << endl;
     Destroy(paths);
     cout << Date( ) << ": cleanup" << endl;
     Cleanup( hbv, inv, rpaths );
     // NOTE EXPENSIVE CONVERSION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     cout << Date( ) << ": hbv -> hbx" << endl;
     hb = HyperBasevectorX(hbv);
    
     dels.clear();

     // Compute edge support
     // TODO make this a function and disentangle from paths index
     cout << Date( ) << ": computing edge support, ";
     vec <int> rs (hb.E( ), 0);
     { // Block to kill rsb
     // compute number of threads based on available memory
     const int T = Max(int(1), Min( int(omp_get_max_threads()), int(MemAvailable( 0.9 ) / (4*hb.E( )) - 1) ) );
     cout << "using " << T << " thread(s)" << endl;
     vec <vec<int>> rsb( T, vec<int>( hb.E( ), 0 ) );
     const int64_t rbatch = rpaths.size()/T+1;
     #pragma omp parallel for num_threads(T)
     for ( int64_t start = 0; start < (int64_t) rpaths.size(); start+=rbatch ) {
          const int64_t stop = Min( start+rbatch, (int64_t) rpaths.size() );
          const int t = start/rbatch;
          for ( int64_t id = start; id < stop; id++ ) {
               ReadPath & p = rpaths[id];
               for ( auto & e : p ) {
                    rsb[t][e]++;
                    rsb[t][inv[e]]++;
               }
          }
     }
     #pragma omp parallel for
     for ( int e = 0; e < hb.E( ); e++ ) {
          for ( int t = 0; t < T; t++ )
               rs[e] += rsb[t][e];
     }
     }
     cout << Date( ) << ": finding hanging ends" << endl;
     // Let's kill some hanging ends.
     // Warning:: if you want to use paths_index here, you have to 
     // recompute it since the graph was modified above.
     // Only rs (read support) is up to date.
     for ( int e = 0; e < hb.E( ); e++ )
     {    int re = inv[e];
          //if ( in_dels[e] ) continue; //already deleted in MowLawn
          
          if ( rs[e] > 0 )
               continue;
          
          if ( hb.Kmers(e) > 200 ) continue; // not thought out!!
          
          // Bool alt = False;
          int v = hb.ToLeft(e), w = hb.ToRight(e);
          if ( hb.To(v).size( ) == 0 && hb.From(v).size( ) == 1 )
          {    dels.push_back( e, re );    }
          else if ( hb.From(w).size( ) == 0 && hb.To(w).size( ) == 1 )
          {    dels.push_back( e, re );    }    }

     UniqueSort(dels);
     cout << Date( ) << ": deleting " << dels.size( ) << " total edges"
          << endl;
     hbv.DeleteEdgesParallel(dels);
     cout << Date( ) << ": cleaning up" << endl;
     Cleanup( hbv, inv, rpaths );
     // NOTE EXPENSIVE CONVERSION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     hb = HyperBasevectorX(hbv);
     paths.parallel_append(rpaths,hb);
     // IF YOU WANT OLD READPATHS, THIS IS THE EARLIEST YOU CAN GET IT!
     // NEED TO PASS FIN_DIR AS ARGUMENT
     /* rpaths.WriteAll(fin_dir + "/a.paths"); */
}


void StageExtension( HyperBasevector& hbv, vec<int>& inv,
          vecbasevector& bases, ObjectManager<VecPQVec>& quals_om,
          ReadPathVecX& paths, Bool const BACK_EXTEND )
{
     STAGE(Extension);
//     VecPQVec& qualsx = quals_om.load_mutable( );
     quals_om.unload();
     VirtualMasterVec<PQVec> vquals(quals_om.filename());
     cout << "Stage Extension, quals loaded, now current mem = " << MemUsageGBString( ) << endl;
     cout << "paths size before ExtendPathsNew " << paths.SizeSum() << endl;
     ExtendPathsNew( hbv, inv, bases, vquals, paths, BACK_EXTEND );
     cout << "Stage Extension, paths extended, now current mem = " << MemUsageGBString( ) << endl;
     cout << "paths size after ExtendPathsNew " << paths.SizeSum() << endl;
     quals_om.unload();
     cout << "after quals unload now current mem = " << MemUsageGBString( ) << endl;
}

void StagePatch(String const& dir, int const K, vecbasevector& bases, 
          ObjectManager<VecPQVec>& quals_om, HyperBasevector& hbv, 
          HyperBasevectorX& hb, ReadPathVec& paths, String pi_file,
          vec<int>& inv, const vec<Bool>& dup, vec<Bool>& bad, 
          vec<DataSet>& datasets, vec<int32_t>& bc, const int max_width, 
          Bool ONE_GOOD, vec<basevector>& closures, vec<pair<int,int>>& pairs, 
          Bool const CG2, const Bool STACKSTER,
          const Bool STACKSTER_ALT, const Bool RESCUE )
{
     // Get started.

     STAGE(Patch);
     ForceAssertEq( dup.size()*2, bases.size() );
     MEM(patch_start);
     // NOTE EXPENSIVE CONVERSION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     hb = HyperBasevectorX(hbv);
     MEM(hbv_copy);

     // Rescue kmers.

     if (RESCUE) Rescue( 0, bases, bc, datasets, hb, inv, paths, closures );

     // Proceed with the rest of the patching.

     VirtualMasterVec<PQVec> vmv( quals_om.filename() );
     quals_om.unload();
//     VecPQVec& quals = quals_om.load_mutable( );
//     MarkBads( hb, bases, quals, paths, bad );
     MarkBads( hb, bases, vmv, paths, bad );

     // Look for edge pairs that we should try to close.
     cout << Date( ) << ": finding edge pairs" << endl;
     FindEdgePairs(
          hb, inv, paths, pi_file, bad, pairs, datasets, bc, ONE_GOOD );
     MEM(edge_pairs);

     // Determine which read ids will be used

     cout << Date() << ": determining subset of interest" << endl;
     vec<uint64_t> read_ids;
#if 1
     vec<bool> in_pairs( hbv.E( ), false);
     for ( auto & p : pairs ) {
          in_pairs[p.first]=true;
          in_pairs[p.second]=true;
          in_pairs[inv[p.first]]=true;
          in_pairs[inv[p.second]]=true;
     }
     
     // edges of interest that are in pairs
     vec<size_type> eoi;
     for ( int e = 0; e < hbv.E( ); e++ ) {
          if (in_pairs[e])
               eoi.push_back(e);
     }

     for ( uint64_t id = 0; id < paths.size(); id++ ) {
          for ( auto e : paths[id] ) {
               if ( in_pairs[e] ) {
                    read_ids.push_back( id );
                    if ( id % 2 == 0 ) {
                         id++; //no need to go to the pair again
                         read_ids.push_back( id );
                    } else
                         read_ids.push_back( id-1 );
                    break;
               }
          }
     }
     ParallelSort( read_ids ); 
#else
     for (size_t i = 0; i < paths.size(); ++i )   // make this a NO-OP for testing
          read_ids.push_back(i);
#endif
     cout << Date() << ": subset is size " << read_ids.size() << endl;

     cout << Date() << ": Destroying paths and populating subsets from disk" << endl;
     Destroy(paths);
     MEM(destroy_paths);

     SubsetMasterVec<PQVec> quals_subset( quals_om.filename(), read_ids );
     SubsetMasterVec<ReadPath> paths_subset( dir+"/a.paths", read_ids );
     SubsetMasterVec<ULongVec> paths_index_subset( pi_file, eoi );

     MEM(after_populate_subset);

#if 0
     // TEMP - check that the NO-OP worked
     cout << Date() << ": checking paths and quals" << endl;
     ForceAssertEq(quals_subset.size(), quals.size() );
     ForceAssertEq(paths_subset.size(), paths.size() );
     for ( size_t i = 0; i < quals_subset.size(); ++i ) {
          if ( paths_subset[i] !=  paths[i] ) FatalErr( "paths bad" );
          for ( size_t j = 0; j < quals_subset[i].size(); ++j )  {
               qualvector q1, q2;
               quals_subset[i].unpack(&q1);
               quals[i].unpack(&q2);
               if ( q1 != q2 ) FatalErr("quals bad" );
          }
     }
     cout << Date() << ": done checking paths and quals" << endl;
#endif


     // Traverse pairs.

     cout << Date( ) << ": start traversing " << ToStringAddCommas( pairs.size( ) )
          << " pairs" << endl;
     double pclock = WallClockTime( );
     const int batches = 1000;
     vec< vec<basevector> > closuresi(batches);
     int64_t NP = pairs.size( );
     vec<int> kmers( hb.E( ) );
     #pragma omp parallel for
     for ( int e = 0; e < hb.E( ); e++ )
          kmers[e] = hb.Kmers(e);
//     auto& quals = quals_om.load_mutable();       // TO-DO: REMOVE ME 
     #pragma omp parallel for
     for ( int bi = 0; bi < batches; bi++ )
     {    for ( int pi = bi * NP / batches; pi < (bi+1) * NP / batches; pi++ )
          {    int e1 = pairs[pi].first, e2 = pairs[pi].second;
               if ( STACKSTER )
               {    vec<basevector> edges = { hb.O(e1), hb.O(e2) };
                    const int VERBOSITY = 0;
                    vec<basevector> f;
                    const Bool EXP = False;
                    vec<int> trim;
                    Stackster( e1, e2, edges, bases, quals_subset, K, datasets, 
                         kmers, inv, dup, paths_subset, paths_index_subset, f, trim,
                         VERBOSITY, STACKSTER_ALT, EXP );
                    closuresi[bi].append(f);    }
               Bool verbose = False;
               vec<basevector> f;
               if ( CG2 ) CloseGap2( hb, inv, bases, quals_subset, paths_subset, paths_index_subset,
                    e1, e2, f, verbose, pi, max_width );
               else CloseGap( hb, inv, bases, quals_subset, paths_subset, paths_index_subset,
                    e1, e2, f, verbose, pi, max_width );
               {    closuresi[bi].append(f);    }    }    }

     MEM(traverse_pairs);
     Destroy(paths_index_subset);
     MEM(destroy_paths_index_subset);
     quals_om.unload( );
     MEM(quals_unload);
     Destroy(paths_subset);
     MEM(destroyed_paths_subset);
     Destroy(quals_subset);
     MEM(destroyed_quals_subset);
     Destroy(bases);
     MEM(destroyed_bases);

     for ( int bi = 0; bi < batches; bi++ )
          closures.append( closuresi[bi] );
     cout << Date( ) << ": found " << closures.size( ) << " closures" << endl;
     cout << TimeSince(pclock) << " used closing pairs" << endl;
     // cout << Date( ) << ": sorting closures" << endl;
     // ParallelSort(closures);

     // Insert gap patches into assembly.

     ForceAssertEq( K, hbv.K( ) );
     HyperBasevector hb3;
     vec<vec<int>> to3( hbv.E( ) );
     vec<int> left3( hbv.E( ) );

     MEM(start_build);
     {    ReadPathVec allx_paths;
          vecbasevector allx;
          BuildAll( allx, hbv, closures.size( ) );
          MEM(build_all);
          allx.append( closures.begin( ), closures.end( ) );
          cout << Date( ) << ": building hb2" << endl;
          cout << "memory in use now = "
               << ToStringAddCommas( MemUsageBytes( ) ) << endl;
          double clock2 = WallClockTime( );
          const int coverage = 4;
          MEM(start_build_big);
          buildBigKHBVFromReads( K, allx, coverage, &hb3, &allx_paths );
          MEM(build_big);
          cout << Date( ) << ": back from buildBigKHBVFromReads" << endl;

          // build to3 and left3 from allx_paths

          for ( int i = 0; i < hb.E( ); ++i )
          {    for ( auto const& p : allx_paths[i] )
                    to3[i].push_back( p );
               left3[i] = allx_paths[i].getFirstSkip();    }

          cout << TimeSince(clock2) << " used in new stuff 2 test" << endl;    }

     cout << "peak mem usage = " << PeakMemUsageGBString( ) << endl;
     MEM(before_load_paths);
     paths.ReadAll( dir + "/a.paths" );
     MEM(load_paths);
     double clock3 = WallClockTime( );
     TranslatePaths( paths, hb3, to3, left3 );
     Destroy(to3), Destroy(left3);
     hbv = hb3;
     hbv.Involution(inv);
     Validate( hbv, paths );

}




void StageBuildGraph( int const K, vecbasevector& bases, ObjectManager<VecPQVec>& quals_om,
          int MIN_QUAL, int MIN_FREQ, int MIN_BC, vec<int32_t> const& bc, int64_t bc_start,
          std::string const GRAPH, double const GRAPHMEM,
          String const& work_dir, String const& read_head, HyperBasevector& hbv, ReadPathVec& paths, vec<int>& inv)
{
     STAGE(BuildGraph);

     cout << Date() << ": barcoded *datatypes* start at " << bc_start << endl;

     MEM(before_graph_creation);
     if ( K != 48 ) FatalErr("remember that we have changes for K=48 not propagated to other K values.");
     buildReadQGraph48(work_dir, read_head, GRAPH, bases, quals_om, False, False, MIN_QUAL, MIN_FREQ, bc_start, MIN_BC, &bc,
                              .75, 0, "", True, False, &hbv, &paths, GRAPHMEM, False );
     cout << Date( ) << ": back from buildReadQGraph" << endl;
     cout << Date( ) << ": memory in use = " << MemUsageGBString( )
          << ", peak = " << PeakMemUsageGBString( ) << endl;
     hbv.Involution(inv);

     // Report assembly statistic.

     if ( hbv.E( ) == 0 )
          Martian::exit("The initial assembly graph built from kmers is empty and execution will be terminated.");
     vec<int> len( hbv.E( ) );
     for ( int e = 0; e < hbv.E( ); e++ )
          len[e] = hbv.Kmers(e);
     Sort(len);
     cout << "N50 edge length = " << N50(len) << endl;

     // Write files.

#if 0          // should be able to retire FixPaths now
     auto bads = FixPaths( hbv, paths ); // needed? - let's find out... see below
     if ( bads ) cout << Date() << ": FixPaths truncated " << bads << " BAD paths" << endl;
#endif

}


