///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Correct1.  Attempt to correct errors in reads.
//
// Inputs and outputs:
// bases: reads
// quals: quality scores associated to reads
//
// Inputs:
// to_edit: which reads to edit
//
// Outputs:
// trim_to: recommended length to trim read to.
//
// Method:
//
// 1. Find all gap-free alignments to a given read, seeded on perfect K-mer matches.
//
// 2. Delete some putatively false friends.
//
// 3. Build matrix and proceed through its columns.
//
// 4. Form the quality score sums for each column.  For the 'losers', delete their
//    top score.
//
// 5. Declare a base to be 'solved' if the winning score is at least 100, >= 10 times
//    the second best, and the second best is at most 100.  Change the base if 
//    needed.
//
// 6. If a base is not solved, truncate the read there and stop.
// 
// Ideas and notes:
//
// 1. Based on an examination of three Fosmids, about 1/400 human reads has an
//    isolated indel, as defined by having ten matching bases on each side.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "CoreTools.h"
#include "PairsManager.h"
#include "ParallelVecUtilities.h"
#include "PrintAlignment.h"
#include "Qualvector.h"
#include "lookup/LookAlign.h"
#include "math/Functions.h"
#include "pairwise_aligners/ClusterAligner.h"
#include "paths/BigMapTools.h"
#include "paths/long/Correct1.h"
#include "paths/long/CorrectByStack.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/FriendAligns.h"
#include "paths/long/MakeKmerStuff.h"
#include "paths/long/ReadStack.h"


void Correct1( String const& tmpDir, const int K, const int max_freq, vecbasevector& bases,
     vecqualvector& quals, const PairsManager& pairs, const vec<Bool>& to_edit, 
     vec<int>& trim_to, const vec<int>& trace_ids, 
     const long_logging_control& log_control, const long_logging& logc,
     const long_heuristics& heur )
{    double clock = WallClockTime( );
     cout << Date( ) << ": entering Correct1 with K = " << K << endl;

     // Build alignments.
     double time = WallClockTime();
     FriendAligner faligns(bases,quals,to_edit,
                             tmpDir+"/correct1."+ToString(K)+".friends",
                             heur.FF_MAKE_ALIGN_IMPL, K,
                             heur.FF_MIN_FREQ,heur.FF_MAX_FREQ,
                             heur.FF_MIN_QUAL,heur.FF_COVERAGE,
                             heur.FF_DOWN_SAMPLE,heur.FF_VERBOSITY);
     cout << "Time use in FriendAligner " << TimeSince( time ) << endl;

     //// Define read starts.

     //vec<int64_t> id1_start( bases.size( ) + 1, -1 );
     //id1_start[ bases.size( ) ] = aligns_all.size( );
     //for ( int64_t j = aligns_all.jsize( ) - 1; j >= 0; j-- )
     //     id1_start[ aligns_all[j].id1 ] = j;
     //for ( int id = (int) bases.size( ) - 1; id >= 0; id-- )
     //     if ( id1_start[id] < 0 ) id1_start[id] = id1_start[id+1];
     //cout << Date( ) << ": read starts defined" << endl;

     // Go through the reads.

     double rclock = WallClockTime( );
     trim_to.resize( bases.size( ) );

     vec<int64_t> use;
     for ( int64_t id1 = 0; id1 < (int64_t) bases.size( ); id1++ )
          if ( to_edit[id1] ) use.push_back(id1);

     static int count_closed(0);

     vecbasevector bases_new(bases);
     vecqualvector quals_new(quals);

     vecbasevector bases_new_special(bases);
     vecqualvector quals_new_special(quals);
     vec<Bool> solved_special( bases.size( ), False );

     int closed = 0, closed_uniquely = 0;

     cout << Date( ) << ": starting main correction loop" << endl;
     #pragma omp parallel for schedule(dynamic, 100)
     for ( int64_t id1x = 0; id1x < (int64_t) use.size( ); id1x++ )
     {    int64_t id1 = use[id1x];
          if ( bases[id1].empty( ) ) continue;
          trim_to[id1] = bases[id1].size( );

          // Fetch 'alignments' to other reads.

          basevector b1 = bases[id1];
          qualvector q1 = quals[id1];
          vec< triple<int,int,Bool> > offset_id_rc2;
          offset_id_rc2.push( 0, id1, False );
          Friends aligns;
          faligns.getAligns( id1, &aligns );
          for ( size_t l = 0; l < aligns.size(); l++ )
          {    offset_id_rc2.push( aligns[l].offset(),
                    aligns[l].readId(), aligns[l].isRC() );    }
          //for ( int64_t l = id1_start[id1]; l < id1_start[id1+1]; l++ )
          //{    offset_id_rc2.push( aligns_all[l].offset,
          //          aligns_all[l].id2, aligns_all[l].rc2 );    }

          // Generate matrix.

          int n = offset_id_rc2.size( );
          int k = bases[id1].size( );
          StackBaseVecVec call( n, StackBaseVec( k, ' ' ) );
          StackQualVecVec callq( n, StackQualVec( k, -1 ) );
          for ( int j = 0; j < n; j++ )
          {    int id2 = offset_id_rc2[j].second, offset = offset_id_rc2[j].first;
               Bool rc2 = offset_id_rc2[j].third;
               const basevector& b2 = bases[id2];
               const qualvector& q2 = quals[id2];
               for ( int p2 = 0; p2 < b2.isize( ); p2++ )
               {    int p1 = p2 + offset;
                    if ( p1 < 0 || p1 >= k ) continue;
                    if ( !rc2 )
                    {    call[j][p1] = b2[p2];
                         callq[j][p1] = q2[p2];    }    
                    else
                    {    call[j][p1] = 3 - b2[ b2.isize( ) - p2 - 1 ];
                         callq[j][p1] = q2[ b2.isize( ) - p2 - 1 ];    }    }    }

          // Perform corrections.

          CorrectByStack( call, callq, bases_new[id1], quals_new[id1],
               trim_to[id1], offset_id_rc2, BinMember(trace_ids, id1) );    

          // Align to reference.

          if ( BinMember( trace_ids, id1 ) && log_control.G != 0 )
          {    basevector b = bases_new[id1];
               #pragma omp critical
               {    for ( int pass = 1; pass <= 2; pass++ )
                    {    if ( pass == 2 ) b.resize( trim_to[id1] );
                         cout << "\nalignments of " 
                              << ( pass == 1 ? "" : "trimmed " )
                              << "trace read " << id1 << " to reference:\n\n";
                         vec<look_align> aligns;
                         ClusterAligner( b, *(log_control.G), log_control.LG, 
                              *(log_control.Glocs), aligns );
                         for ( int j = 0; j < aligns.isize( ); j++ )
                         {    look_align& la = aligns[j];
                              la.query_id = id1;
                              int g = la.target_id;
                              la.PrintReadableBrief( cout, "read_" + ToString(id1),
                                   ToString(g) );
                              la.PrintVisual( cout, b, quals_new[id1], 
                                   (*(log_control.G))[g], 
                                   0 );    }    }    }    }    }

     // Save results.

     bases = bases_new;
     quals = quals_new;
     cout << "\n" << TimeSince(clock) << " used in Correct1, of which "
          << TimeSince(rclock) << " was in going through reads\n" << endl;    }
