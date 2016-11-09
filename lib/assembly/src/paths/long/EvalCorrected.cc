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
#include "ParallelVecUtilities.h"
#include "PrintAlignment.h"
#include "efasta/EfastaTools.h"
#include "kmers/KmerRecord.h"
#include "lookup/LookAlign.h"
#include "pairwise_aligners/ClusterAligner.h"
#include "paths/BigMapTools.h"
#include "paths/RemodelGapTools.h"
#include "paths/long/CreateGenome.h"
#include "paths/long/EvalCorrected.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/MakeKmerStuff.h"
#include "random/Random.h"

void EvalCorrected( 
     const VecEFasta& corrected0,             // corrected reads
     const vec<int>& cid,                     // ids of corrected reads
     const ref_data& ref,
     const long_logging_control& log_control,
     const long_logging& logc )
{
     const vecbasevector& G = ref.G;
     const int LG = ref.LG;
     const VecIntPairVec& Glocs = ref.Glocs;

     if ( corrected0.empty( ) )
     {    cout << "\nThere are no reads to evaluate!\n";
          return;    }
     double clock = WallClockTime( );
     vecbasevector all;
     all.Append(G);

     // Expanded corrected reads.  This is exponential and thus totally unsound.

     vec< vec<basevector> > corrected( corrected0.size( ) );
     for ( size_t id = 0; id < corrected0.size( ); id++ )
          corrected0[id].ExpandTo( corrected[id] );
     vec<basevector> correctedv1;
     for ( int id = 0; id < corrected.isize( ); id++ )
     {    for ( int j = 0; j < corrected[id].isize( ); j++ )
               correctedv1.push_back( corrected[id][j] );    }

     // Compute ambiguity rate.

     double ambiguities = 0.0, amb_denom = 0.0;
     for ( int i = 0; i < corrected.isize( ); i++ )
     {    if ( corrected[i].size( ) > 1 ) 
               ambiguities += log2( corrected[i].size( ) );
          double d = 0.0;
          for ( int j = 0; j < corrected[i].isize( ); j++ )
               d += corrected[i][j].size( );
          amb_denom += d / corrected[i].size( );    }
     double amb_rate = ambiguities / amb_denom;

     // We use both the reads and their reverse complements.  This isn't quite 
     // right -- we should only be using one, depending on which matches the
     // reference better.

     const int M = 100;
     vec< triple<kmer<M>,int,int> > kmers_plus;
     MakeKmerLookup0( all, kmers_plus );
     vec< kmer<M> > kmers( kmers_plus.size( ) );
     for ( size_t i = 0; i < kmers.size( ); i++ )
          kmers[i] = kmers_plus[i].first;
     ParallelUniqueSort(kmers);
     vec<Bool> found( kmers.size( ), False );
     cout << Date( ) << ": start looking for reference kmers" << endl;
     #pragma omp parallel for
     for ( int j = 0; j < correctedv1.isize( ); j++ )
     {    basevector b = correctedv1[j];
          kmer<M> x;
          for ( int pass = 1; pass <= 2; pass++ )
          {    if ( pass == 2 ) b.ReverseComplement( );
               for ( int l = 0; l <= b.isize( ) - M; l++ )
               {    x.SetToSubOf( b, l );
                    int64_t p = BinPosition( kmers, x );
                    if ( p >= 0 ) found[p] = True;    }    }    }

     // Print missing kmers.

     if ( logc.SHOW_MISSING > 0 )
     {    cout << "\n";
          vec< kmer<M> > missing;
          for ( int64_t i = 0; i < (int64_t) found.size( ); i++ )
               if ( !found[i] ) missing.push_back( kmers[i] );
          Sort(missing);
          vec< pair<int,int> > mlocs;
          for ( int64_t j = 0; j < (int64_t) kmers_plus.size( ); j++ )
          {    if ( BinMember( missing, kmers_plus[j].first ) )
                    mlocs.push( kmers_plus[j].second, kmers_plus[j].third );    }
          Sort(mlocs);
          int count = 0;
          for ( int i = 0; i < mlocs.isize( ); i++ )
          {    int gid = mlocs[i].first;
               int j;
               for ( j = i + 1; j < mlocs.isize( ); j++ )
               {    if ( mlocs[j].first != gid ) break;
                    if ( mlocs[j].second != mlocs[j-1].second + 1 ) break;    }
               int gstart = mlocs[i].second, gstop = mlocs[j-1].second + M;
               basevector b( G[gid], gstart, gstop - gstart );
               b.Print( cout, "missing_" + ToString(++count) + ", " + ToString(j-i) 
                    + " kmers, " + ToString(gid) + "." + ToString(gstart) + "-"
                    + ToString(gstop) );
               if ( count == logc.SHOW_MISSING ) break;
               i = j - 1;    }
          cout << "\n";    }

     // Select reads for alignment.

     int unaligned = 0;
     vec<double> err_rate;
     vec<int> ids;
     if ( (int) corrected.size( ) <= logc.EVAL_SAMPLE )
     {    for ( int i = 0; i < (int) corrected.size( ); i++ )
               ids.push_back(i);    }
     else
     {    set<int> idsx;
          while( (int) idsx.size( ) < logc.EVAL_SAMPLE )
          {    int x = randomx( ) % corrected.size( );
               idsx.insert(x);    }
          for ( set<int>::iterator i = idsx.begin( ); i != idsx.end( ); i++ )
               ids.push_back(*i);    }

     // Compute perfect alignments.

     vecbasevector correctedv0;
     vec< pair<int,int> > origin0;
     for ( int i = 0; i < ids.isize( ); i++ )
     {    int id = ids[i];
          for ( int j = 0; j < corrected[id].isize( ); j++ )
          {    correctedv0.push_back_reserve( corrected[id][j] );
               origin0.push( i, j );    }    }
     cout << Date( ) << ": starting perfect alignment of "
          << ToStringAddCommas( correctedv0.size( ) ) << " read instantiations" 
          << endl;
     vec<Bool> perfect( ids.size( ), False );
     #pragma omp parallel for
     for ( size_t r = 0; r < correctedv0.size( ); r++ )
     {    int i = origin0[r].first, j = origin0[r].second;
          if ( perfect[i] ) continue;
          vec<placementy> p = FindGenomicPlacementsY( 
               0, corrected[ ids[i] ][j], LG, G, Glocs );
          if ( p.nonempty( ) )
          {    
               #pragma omp critical
               {    perfect[i] = True;    }    }    }

     // Set up for imperfect alignment.

     vecbasevector correctedv;
     vec< pair<int,int> > origin;
     int todo = 0;
     for ( int i = 0; i < ids.isize( ); i++ )
     {    int id = ids[i];
          if ( !perfect[i] )
          {    for ( int j = 0; j < corrected[id].isize( ); j++ )
               {    correctedv.push_back_reserve( corrected[id][j] );
                    origin.push( todo, i );    }
               todo++;    }    }
     cout << Date( ) << ": done, aligning " 
          << ToStringAddCommas( correctedv.size( ) ) << " read instantiations" 
          << endl;

     // Print initial results.

     if ( logc.MIN_ERRS_TO_PRINT >= 0 )
     {    cout << "\n==========================================================="
               << "====================\n";
          cout << "alignments of reads having at least " << logc.MIN_ERRS_TO_PRINT
               << " errors" << endl;
          cout << "==========================================================="
               << "====================\n";    }
     for ( int i = 0; i < ids.isize( ); i++ )
     {    if ( perfect[i] > 0 )
          {    err_rate.push_back(0);
               if ( logc.MIN_ERRS_TO_PRINT == 0 )
               {    cout << "\nalignment of read " << cid[ ids[i] ] 
                         << " is perfect\n";    }    }    }

     // Start main alignment loop.

     vec< vec<double> > er(todo);
     vec< vec<String> > reports(todo);
     #pragma omp parallel for schedule(dynamic, 1)
     for ( size_t z = 0; z < correctedv.size( ); z++ )
     {    int t = origin[z].first, i = origin[z].second;
          int id = ids[i];
          vec<look_align> xaligns;
          ClusterAligner( correctedv[z], G, LG, Glocs, xaligns );
          if ( xaligns.nonempty( ) )
          {    int e = xaligns[0].Errors( );
               ostringstream out;
               if ( logc.MIN_ERRS_TO_PRINT >= 0 && e >= logc.MIN_ERRS_TO_PRINT )
               {    basevector b = correctedv[z];
                    if ( xaligns[0].Rc1( ) ) b.ReverseComplement( );
                    out << "\nalignment of read " << cid[id] 
                         << " is " << ( xaligns[0].Rc1( ) ? "rc" : "fw" )
                         << " and has " << e << " errors\n";
                    out << "aligned at " << xaligns[0].target_id
                         << "." << xaligns[0].a.pos2( ) << "-"
                         << xaligns[0].a.Pos2( ) << ":"
                         << ( xaligns[0].Fw1( ) ? "fw" : "rc" ) << "\n";
                    PrintVisualAlignment( True, out, b,
                         G[ xaligns[0].target_id ], xaligns[0].a );    }
               #pragma omp critical
               {    er[t].push_back( xaligns[0].ErrorRate( ) );    
                    reports[t].push_back( out.str( ) );    }    }    }

     // Summarize results.

     for ( int t = 0; t < todo; t++ )
     {    if ( er[t].empty( ) ) unaligned++;
          else
          {    SortSync( er[t], reports[t] );
               err_rate.push_back( er[t][0] );
               cout << reports[t][0];    }    }
     Sort(err_rate);
     cout << "\n" << PERCENT_RATIO( 4, Sum(found), (int64_t) kmers.size( ) ) 
          << " sensitivity" 
          << " (fraction of genomic 100-mers present in corrected reads)" << endl;
     cout << PERCENT_RATIO( 3, unaligned, ids.isize( ) )
          << " of corrected reads are unaligned" << endl;
     cout << setiosflags(ios::fixed) << setprecision(2) << 10000.0 * Mean(err_rate)
          << "%% error rate for aligned reads" << endl;
     cout << setprecision(2) << 10000.0 * amb_rate << resetiosflags(ios::fixed) 
          << "%% ambig rate for aligned reads" << endl;
     REPORT_TIME( clock, "in EvalCorrected" );    }
