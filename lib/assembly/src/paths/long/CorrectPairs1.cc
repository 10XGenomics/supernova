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
#include "PrintAlignment.h"
#include "Qualvector.h"
#include "efasta/EfastaTools.h"
#include "lookup/LookAlign.h"
#include "math/Functions.h"
#include "pairwise_aligners/ClusterAligner.h"
#include "pairwise_aligners/SmithWatBandedA.h"
#include "paths/BigMapTools.h"
#include "paths/long/CorrectPairs1.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/FriendAligns.h"
#include "paths/long/MakeKmerStuff.h"
#include "paths/long/ReadStack.h"
#include "random/Bernoulli.h"

namespace { // open anonymous namespace

Bool cmp_ho_start_stop( const ho_interval& h1, const ho_interval& h2 )
{    if ( h1.Start( ) < h2.Start( ) ) return True;
     if ( h1.Start( ) > h2.Start( ) ) return False;
     return h1.Stop( ) > h2.Stop( );    }

} // close anonymous namespace

void CorrectPairs1( String const& tmpDir, const int K, const int max_freq, vecbasevector& bases,
     vecqualvector& quals, const PairsManager& pairs, const vec<Bool>& to_edit, 
     const vec<int>& trace_ids, const long_heuristics& heur, 
     const long_logging_control& log_control, const long_logging& logc, 
     VecEFasta& corrected )
{    double clock = WallClockTime( );
     if (logc.STATUS_LOGGING)
          cout << Date( ) << ": entering CorrectPairs1 with K = " << K << endl;

     // Build alignments.

     FriendAligner faligns(bases,quals,to_edit,
                             tmpDir+"/correctPairs1."+ToString(K)+".friends",
                             heur.FF_MAKE_ALIGN_IMPL, K,
                             heur.FF_MIN_FREQ,heur.FF_MAX_FREQ,
                             heur.FF_MIN_QUAL,heur.FF_COVERAGE,
                             heur.FF_DOWN_SAMPLE,heur.FF_VERBOSITY);

     //// Define read starts.

     //vec<int64_t> id1_start( bases.size( ) + 1, -1 );
     //id1_start[ bases.size( ) ] = aligns_all.size( );
     //for ( int64_t j = aligns_all.jsize( ) - 1; j >= 0; j-- )
     //     id1_start[ aligns_all[j].id1 ] = j;
     //for ( int id = (int) bases.size( ) - 1; id >= 0; id-- )
     //     if ( id1_start[id] < 0 ) id1_start[id] = id1_start[id+1];
     //cout << Date( ) << ": read starts defined" << endl;

     // Go through the reads.

     vec<int64_t> use;
     for ( int64_t id1 = 0; id1 < (int64_t) bases.size( ); id1++ )
     {    int64_t id2 = pairs.getPartnerID(id1);
          if ( to_edit[id1] && to_edit[id2] && bases[id1].size( ) > 0 && id2 < id1 )
               use.push_back(id1);    }

     static int count_closed(0);
     /*
     {    cout << Date( ) << ": begin raise" << endl;
          vecqualvector quals2(quals);
          double timer1 = 0, timer2 = 0;
          #pragma omp parallel for schedule(dynamic, 100)
          for ( int64_t id1x = 0; id1x < (int64_t) use.size( ); id1x++ )
          {    int64_t id1 = use[id1x];
               if ( bases[id1].empty( ) ) continue;
               double clock1 = WallClockTime( );
               readstack stack( id1, aligns_all, id1_start[id1], id1_start[id1+1],
                    readstack::strict, bases, quals, pairs );
               timer1 += WallClockTime( ) - clock1;
               double clock2 = WallClockTime( );
               stack.Raise1( 0, 21, True );
               timer2 += WallClockTime( ) - clock2;
               for ( int i = 0; i < bases[id1].isize( ); i++ )
                    quals2[id1][i] = stack.Qual( 0, i );    }
          quals = quals2;
          cout << timer1 << " seconds spent constructing" << endl;
          cout << timer2 << " seconds spent raising" << endl;
          if (logc.STATUS_LOGGING) cout << Date( ) << ": done" << endl;    }
     */

     int closed = 0, closed_uniquely = 0;

     REPORT_TIME( clock, "used in pair correction head" );
     if (logc.STATUS_LOGGING)
     {    ReportPeakMem( "start main correction, " + ToString( use.size( ) )
               + " reads" );    }
     int batch = (int64_t) use.size( ) / omp_get_max_threads( );
     batch = Min( 100, Max( 1, batch ) );
     #pragma omp parallel for schedule(dynamic, batch)
     for ( int64_t id1x = 0; id1x < (int64_t) use.size( ); id1x++ )
     {    double aclock = WallClockTime( );
          int64_t id1 = use[id1x];
          // Build stacks.

          ostringstream out;
          int64_t id1p = pairs.getPartnerID(id1);
          Bool trace = BinMember(trace_ids, id1) || BinMember(trace_ids, id1p);
          // Get aligns for this read
          Friends aligns;
          faligns.getAligns( id1, &aligns );
          readstack stack1( id1, aligns, 0, aligns.size(),
               readstack::right_extended, bases, quals, pairs );
          Friends aligns_p;
          faligns.getAligns( id1p, &aligns_p );
          readstack stack2( id1p, aligns_p, 0, aligns_p.size(),
               readstack::right_extended, bases, quals, pairs );
          REPORT_TIMEX( aclock, "used making stacks for pair correction", out );
          if ( stack1.Rows( ) > heur.MAX_STACK || stack2.Rows( ) > heur.MAX_STACK )
               continue;
          double bclock = WallClockTime( );

          if (trace)
          {
               #pragma omp critical
               {    cout << "\ntracing reads id1 = " << id1 
                         << ", id1p = " << id1p << endl;
                    cout << "\ninitial stack1:\n";
                    stack1.Print(cout);
                    cout << "initial stack2:\n";
                    stack2.Print(cout);    }    }

          // Filter out low-quality reads.

          int total_bases = 0, total_qual = 0;
          for ( int i = 0; i < stack1.Cols( ); i++ )
          {    if ( stack1.Def(0,i) ) total_bases++;
               if ( stack1.Qual(0,i) >= 2 ) total_qual += stack1.Qual(0,i);    }
          for ( int i = 0; i < stack2.Cols( ); i++ )
          {    if ( stack2.Def(0,i) ) total_bases++;
               if ( stack2.Qual(0,i) >= 2 ) total_qual += stack2.Qual(0,i);    }
          double this_qual = double(total_qual)/double(total_bases);

          int bases_all = 0, total_all = 0;
          vec<int> ids_all;
          for ( int pass = 1; pass <= 2; pass++ )
          {    const readstack& s = ( pass == 1 ? stack1 : stack2 );
               for ( int j = 0; j < s.Rows( ); j++ )
                    ids_all.push_back( s.Id(j) );    }
          UniqueSort(ids_all);
          for ( int m = 0; m < ids_all.isize( ); m++ )
          {    int id = ids_all[m];
               for ( int j = 0; j < (int) quals[id].size( ); j++ )
               {    bases_all++;
                    if ( quals[id][j] >= 2 ) total_all += quals[id][j];    }    }
          if ( bases_all == 0 ) bases_all++;
          double all_qual = double(total_all)/double(bases_all);

          if (trace)
          {
               #pragma omp critical
               {    cout << "\ntracing reads id1 = " << id1 
                         << ", id1p = " << id1p << endl;
                    PRINT2( this_qual, all_qual );    }    }
          REPORT_TIMEX( bclock, "used qual filtering in pair correction", out );

          if ( all_qual - this_qual > heur.CP_MAX_QDIFF ) 
          {    if ( omp_get_thread_num( ) == 0 ) cout << out.str( );
               continue;    }

          // Remove friends having inadequate glue to the founder.

          double dclock = WallClockTime( );
          vec<Bool> suspect;
          stack1.FlagNoise(suspect);
          stack1.Erase(suspect);
          stack2.FlagNoise(suspect);
          stack2.Erase(suspect);

          // Proceed.

          const int q_solid = 30;
          stack1.Raise1(0);
          stack1.MotifDiff(1,suspect);
          stack1.Erase(suspect);

          if (trace)
          {
               #pragma omp critical
               {    cout << "\ntracing reads id1 = " << id1 
                         << ", id1p = " << id1p << endl;
                    cout << "\npost-motif stack1:\n";
                    stack1.Print(cout);    }    }

          stack1.HighQualDiff( q_solid, 1, suspect );
          stack1.Erase(suspect);

          if (trace)
          {
               #pragma omp critical
               {    cout << "\ntracing reads id1 = " << id1 
                         << ", id1p = " << id1p << endl;
                    cout << "\npost-hq stack1:\n";
                    stack1.Print(cout);    }    }

          stack2.Raise1(0);
          stack2.MotifDiff(1,suspect);
          stack2.Erase(suspect);
          stack2.HighQualDiff( q_solid, 1, suspect );
          stack2.Erase(suspect);

          if (trace)
          {
               #pragma omp critical
               {    cout << "\ntracing reads id1 = " << id1 
                         << ", id1p = " << id1p << endl;
                    cout << "\npost-hq stack2:\n";
                    stack2.Print(cout);    }    }

          // PairWeak2( stack1, stack2, suspect1, suspect2 );
          // stack1.Erase(suspect1);
          // stack2.Erase(suspect2);

          // Identify reads having a Q30 base that is unsupported.

          /*
          Bool dam = False;
          for ( int pass = 1; pass <= 2; pass++ )
          {    const int64_t id = ( pass == 1 ? id1 : id1p );
               const basevector& b1 = bases[id];
               const readstack& s = ( pass == 1 ? stack1 : stack2 );
               for ( int j = 0; j < b1.isize( ); j++ )
               {    if ( s.Qual(0,j) < 30 ) continue;
                    int count30 = 0;
                    for ( int i = 1; i < s.Rows( ); i++ )
                    {    if ( s.Qual(i,j) >= 30 && s.Base(i,j) == s.Base(0,j)
                              && s.Pid(i) != s.Pid(0) )
                         {    count30++;    }    }
                    if ( count30 == 0 )
                    {    dam = True;
                         break;    }    }    }
          if (dam)
          {
               #pragma omp critical
               {    damaged++;    }
               continue;    }
          */

          // Try to recruit other reads.  Computationally unsound.

          /*
          stack1.Recruit( K, aligns_all, id1_start, bases, quals, pairs );
          stack1.Raise1(0);
          stack1.MotifDiff(1,suspect);
          stack1.Erase(suspect);
          stack1.HighQualDiff( q_solid, 1, suspect );
          stack1.Erase(suspect);
          stack2.Recruit( K, aligns_all, id1_start, bases, quals, pairs );
          stack2.Raise1(0);
          stack2.MotifDiff(1,suspect);
          stack2.Erase(suspect);
          stack2.HighQualDiff( q_solid, 1, suspect );
          stack2.Erase(suspect);
          */
          REPORT_TIMEX( dclock, "used in filtering pair correction", out );
          double d2clock = WallClockTime( );

          stack2.Reverse( );
          basevector con1, con2;
          qualvector conq1, conq2;
          stack1.Consensus1( con1, conq1 );
          stack2.Consensus1( con2, conq2 );
          const int L = 20;

          // Check for single-base indel (second method).

          /*
          if ( offsets.size( ) == 2 && offsets[1] - offsets[0] == 1 )
          {    align a;
               int errors;
               SmithWatBandedA( con1, con2, offsets[0], 2, a, errors, 0, 1, 1 );
               if ( errors == 1 )
               {
                    #pragma omp critical
                    {    indel++;    }
                    continue;    }    }
          */
          REPORT_TIMEX( d2clock, "used in consensus in pair correction", out );
          double d3clock = WallClockTime( );

          vec<int> offsets = GetOffsets1( stack1, stack2, 0, heur.DELTA_MIS );
          vec<int> final_offsets = offsets;

          // Trace.

          if (trace)
          {
               #pragma omp critical
               {    cout << "\ntracing reads id1 = " << id1 
                         << ", id1p = " << id1p << endl;
                    cout << "final offsets = " << printSeq(final_offsets) 
                         << endl;    
                    vec<int> moffsets = GetOffsets1( 
                         stack1, stack2, 2, heur.DELTA_MIS );    }    }

          // For each offset, create the merged stack associated to it.

          vec<basevector> closures;
          vec<qualvector> closuresq;
          vec<readstack> stacks;
          vec<int> closureso;
          REPORT_TIMEX( d3clock, "used in offset creation in pair correction", out );
          double eclock = WallClockTime( );
          for ( int oj = 0; oj < offsets.isize( ); oj++ )
          {    int minq_floor = ( offsets.size( ) > 1 ? heur.CP_MINQ_FLOOR 
                    : 5 /* Min( heur.CP_MINQ_FLOOR, 5 ) */ );
               int min_glue_floor = ( offsets.size( ) > 1 ? heur.CP_MIN_GLUE 
                    : Min( heur.CP_MIN_GLUE, 20 ) );
               readstack stack(stack1);
               stack.Merge( stack2, offsets[oj] );
               stack.SortByPid( pairs.getPairID(id1), 0, stack1.Rows( ) );
               stack.Unique( );
               stack.Raise1(0), stack.Raise1(1);

               if (trace)
               {    
                    #pragma omp critical
                    {    cout << "\ninitial merged stack for " << stack.Id(0) 
                              << endl;
                         stack.Print(cout);    }    }

               /*
               stack.MotifDiff(2,suspect);
               if ( suspect[0] || suspect[1] ) continue;
               stack.Erase(suspect);
               */

               stack.HighQualDiff( q_solid, 2, suspect );
               if ( suspect[0] || suspect[1] ) continue;
               stack.Erase(suspect);

               if (trace)
               {    
                    #pragma omp critical
                    {    cout << "\npost hq-diff merged stack for " << stack.Id(0) 
                              << endl;
                         stack.Print(cout);    }    }

               // stack.AddPartners( 32, 2, bases, quals, pairs );

               stack.PairWeak1( suspect );
               if ( suspect[0] || suspect[1] ) continue;
               stack.Erase(suspect);

               // stack.CleanColumns(2,suspect);
               // if ( !suspect[0] && !suspect[1] ) stack.Erase(suspect);

               int start, stop;
               for ( start = 0; start < stack.Cols( ); start++ )
                    if ( stack.Def(0,start) ) break;
               for ( stop = stack.Cols( ) - 1; stop >= 0; stop-- )
                    if ( stack.Def(1,stop) ) break;
               stop++;
               Bool weird = ( start >= stop );
               if ( !weird ) stack.Trim( start, stop );

               // Create and edit consensus.

               basevector con;
               qualvector conq;
               stack.StrongConsensus1( con, conq, heur.CP_RAISE_ZERO );

               const int protected_bases = 10;
               const int q_to_protect = 20;
               for ( int j = 0; j < protected_bases; j++ )
               {    if ( j >= stack.Cols( ) ) break;
                    if ( stack.Def(0,j) && stack.Base(0,j) != con[j]
                         && stack.Qual(0,j) >= q_to_protect )
                    {    con.Set( j, stack.Base(0,j) );
                         conq[j] = stack.Qual(0,j);    }    }
               for ( int j = 0; j < protected_bases; j++ )
               {    int jr = stack.Cols( ) - j - 1;
                    if ( jr < 0 ) break;
                    if ( stack.Def(1,jr) && stack.Base(1,jr) != con[jr]
                         && stack.Qual(1,jr) >= q_to_protect )
                    {    con.Set( jr, stack.Base(1,jr) );
                         conq[jr] = stack.Qual(1,jr);    }    }

               for ( int j = 0; j < con.isize( ); j++ )
               {    if ( stack.Qual(0,j) >= 30 && stack.Base(0,j) != con[j] )
                         conq[j] = 0;
                    if ( stack.Qual(1,j) >= 30 && stack.Base(1,j) != con[j] )
                         conq[j] = 0;    }

               // Test for suspicious inconsistencies between the founder
               // and the consensus.

               for ( int j = 0; j < con.isize( ); j++ )
               for ( int m = 0; m < 2; m++ )
               {    if ( !stack.Def(m,j) || stack.Base(m,j) == con[j] ) 
                         continue;
                    const int flank = 5;
                    const int min_mult = 3;
                    if ( j < flank || j + flank >= con.isize( ) ) continue;
                    Bool mismatch = False;
                    for ( int l = 0; l < flank; l++ )
                    {    if ( stack.Base(m,j-l-1) != con[j-l-1] )
                              mismatch = True;
                         if ( stack.Base(m,j+l+1) != con[j+l+1] )
                              mismatch = True;    }
                    if (mismatch) continue;
                    int mult = 0;
                    for ( int r = 2; r < stack.Rows( ); r++ )
                    {    Bool mismatch = False;
                         for ( int p = j - flank; p <= j + flank; p++ )
                         {    if ( stack.Base(r,p) != stack.Base(m,p) )
                              {    mismatch = True;
                                   break;    }    }
                         if (mismatch) continue;
                         if ( ++mult == min_mult ) break;    }
                    if ( mult == min_mult ) conq[j] = 0;    }

               // Attempt to recover conflicted columns.

               vec<Bool> to_del( stack.Rows( ), False );
               for ( int j = 0; j < stack.Cols( ); j++ )
               {    if ( conq[j] < minq_floor )
                    {    const int qmin = 2;
                         const int qdelta = 10;
                         if ( stack.Qual(0,j) < qmin && stack.Qual(1,j) < qmin )
                              continue;
                         if ( stack.Qual(0,j) >= qmin && stack.Qual(1,j) >= qmin
                              && stack.Base(0,j) != stack.Base(1,j)
                              && Abs( stack.Qual(0,j) - stack.Qual(1,j) ) < qdelta )
                         {    continue;    }
                         char b;
                         if ( stack.Qual(0,j) >= qmin 
                              && stack.Qual(0,j) >= stack.Qual(1,j) ) 
                         {    b = stack.Base(0,j);    }
                         else b = stack.Base(1,j);
                         for ( int i = 2; i < stack.Rows( ); i++ )
                         {    if ( stack.Qual(i,j) >= qmin && stack.Base(i,j) != b )
                              {    if (trace)
                                   {    cout << "conflict in column " << j
                                             << ", so deleting read "
                                             << stack.Id(i) << " from merged stack "
                                             << " for read " << stack.Id(0)
                                             << endl;    }
                                   to_del[i] = True;    }    }    }    }
               stack.Erase(to_del);
               stack.StrongConsensus1( con, conq, heur.CP_RAISE_ZERO );

               for ( int j = 0; j < protected_bases; j++ )
               {    if ( j >= stack.Cols( ) ) break;
                    if ( stack.Def(0,j) && stack.Base(0,j) != con[j]
                         && stack.Qual(0,j) >= q_to_protect )
                    {    con.Set( j, stack.Base(0,j) );
                         conq[j] = stack.Qual(0,j);    }    }
               for ( int j = 0; j < protected_bases; j++ )
               {    int jr = stack.Cols( ) - j - 1;
                    if ( jr < 0 ) break;
                    if ( stack.Def(1,jr) && stack.Base(1,jr) != con[jr]
                         && stack.Qual(1,jr) >= q_to_protect )
                    {    con.Set( jr, stack.Base(1,jr) );
                         conq[jr] = stack.Qual(1,jr);    }    }

               const int qfloor = 20;
               vec<placementy> p;
                    
               // Probably doesn't do anything:

               Bool yes1 = False, yes2 = False;
               for ( int j = 0; j < stack.Cols( ); j++ )
               {    if ( stack.Def(0,j) ) yes1 = True;
                    if ( stack.Def(1,j) ) yes2 = True;    }
               if ( !yes1 || !yes2 ) continue;

               // Determine minimum consensus quality.

               int minq = 1000000000;
               for ( int j = 0; j < con.isize( ); j++ )
                    minq = Min( minq, (int) conq[j] );

               // Check for glue.

               vec<ho_interval> agree;
               for ( int i = 0; i < stack.Rows( ); i++ )
               for ( int j = 0; j < stack.Cols( ); j++ )
               {    if ( stack.Base(i,j) != con[j] ) continue;
                    int k;
                    for ( k = j + 1; k < stack.Cols( ); k++ )
                         if ( stack.Base(i,k) != con[k] ) break;
                    if ( k - j >= 40 ) agree.push( j, k );
                    j = k - 1;    }
               sort( agree.begin( ), agree.end( ), cmp_ho_start_stop );
               vec<Bool> to_delete( agree.size( ), False );
               for ( int i = 0; i < agree.isize( ); i++ )
               {    int j;
                    for ( j = i + 1; j < agree.isize( ); j++ )
                    {    if ( agree[j].Stop( ) > agree[i].Stop( ) ) break;    }
                    for ( int k = i + 1; k < j; k++ )
                         to_delete[k] = True;
                    i = j - 1;    }
               EraseIf( agree, to_delete );
               int min_glue;
               if ( agree.empty( ) || agree[0].Start( ) > 0 ) min_glue = 0;
               else
               {    min_glue = agree[0].Length( );
                    int stop = agree[0].Stop( );
                    for ( int i = 1; i < agree.isize( ); i++ )
                    {    if ( agree[i].Stop( ) > stop )
                         {    min_glue = Min( min_glue, stop - agree[i].Start( ) );
                              stop = agree[i].Stop( );    }    }
                    if ( stop < con.isize( ) ) min_glue = 0;    }

               // Trace.

               if (trace)
               {
                    #pragma omp critical
                    {    cout << "\ntracing reads id1 = " << id1 << ", id1p = "
                              << id1p << "\noffset = " << offsets[oj] << endl;
                         PRINT4( minq, minq_floor, min_glue, min_glue_floor );
                         if ( minq < minq_floor )
                         {    cout << "low-quality positions:";
                              int count = 0;
                              for ( int j = 0; j < (int) conq.size( ); j++ )
                              {    if ( conq[j] == minq )
                                   {    if ( ++count <= 10 ) cout << " " << j;
                                        else cout << " ...";    }    }
                              cout << "\n";    }
                         cout << "\nalignment of closure:\n";
                         vec<look_align> aligns;
                         ClusterAligner( con, *(log_control.G), log_control.LG, 
                              *(log_control.Glocs), aligns );
                         for ( int j = 0; j < aligns.isize( ); j++ )
                         {    int g = aligns[j].target_id;
                              aligns[j].PrintReadableBrief( cout,
                                   "con", ToString(g) );
                              aligns[j].PrintVisual( cout, con, conq,
                                   (*(log_control.G))[g], False );    }
                         /*
                         vec<basevector> cs = stack.Consensuses2( 30, 2 );
                         cout << "found " << cs.size( ) << " consensuses\n";
                         const int max_cs_print = 100;
                         for ( int l = 0; 
                              l < Min( max_cs_print, cs.isize( ) ); l++ )
                         {    cout << "alignment of consensus " << l << endl;
                              ClusterAligner( cs[l], *(log_control.G), 
                                   log_control.LG, *(log_control.Glocs), 
                                   aligns );
                              for ( int j = 0; j < aligns.isize( ); j++ )
                              {    int g = aligns[j].target_id;
                                   aligns[j].PrintReadableBrief( cout,
                                        "cs[l]", ToString(g) );
                                   aligns[j].PrintVisual( cout, cs[l],
                                        (*(log_control.G))[g], False );    }
                                   }
                         */

                         cout << "stack1:\n";
                         stack1.Print(cout);
                         cout << "stack2:\n";
                         stack2.Print(cout);
                         cout << "stack:\n";
                         stack.Print(cout);    }    }

               if ( minq >= minq_floor && min_glue >= min_glue_floor )
               {    if (trace)
                    {
                         #pragma omp critical
                         {    cout << "\ntracing reads id1 = " << id1 << ", id1p = "
                                   << id1p << "\noffset = " << offsets[oj] << endl;
                              cout << "saving closure\n";    }    }
                    closures.push_back(con), closuresq.push_back(conq);    
                    closureso.push_back( offsets[oj] );
                    stacks.push_back(stack);    }    }
          REPORT_TIMEX( eclock, "used processing offsets in pair corr", out );
          double fclock = WallClockTime( );

          // Save closures.  Currently we only save a longest one.
          // (NO, CHANGED.)

          if ( closures.nonempty( ) )
          {    
               if (heur.CP_CONDENSE_HOMOPOLYMERS)
               {    efasta eclosure(closures);
                    if ( eclosure.AmbEventCount( ) == 1 )
                    {    String mid = eclosure.Between( "{", "}" );
                         mid.GlobalReplaceBy( ",", "" );
                         Bool homopolymer = True;
                         for ( int j = 1; j < mid.isize( ); j++ )
                              if ( mid[j] != mid[0] ) homopolymer = False;
                         if (homopolymer)
                         {    corrected[id1] = eclosure;
                              eclosure.ReverseComplement( );
                              corrected[id1p] = eclosure;
                              #pragma omp critical
                              {    count_closed++;    }
                              REPORT_TIMEX( fclock, 
                                   "used in closures, pair corr", out );
                              if ( omp_get_thread_num( ) == 0 ) cout << out.str( );
                              continue;    }    }    }

               int mc = 1000000000;
               for ( int i = 0; i < closures.isize( ); i++ )
                    mc = Min( mc, closures[i].isize( ) );
               basevector left(mc), right(mc);

               #pragma omp critical
               {    closed++;
                    if ( closures.solo( ) ) closed_uniquely++;    }
               
               for ( int j = 0; j < mc; j++ )
               {    vec<char> seen;
                    for ( int i = 0; i < closures.isize( ); i++ )
                         seen.push_back( closures[i][j] );
                    UniqueSort(seen);
                    if ( seen.solo( ) ) left.Set( j, seen[0] );
                    else
                    {    left.resize(j);
                         break;    }    }    
               for ( int j = 0; j < mc; j++ )
               {    vec<char> seen;
                    for ( int i = 0; i < closures.isize( ); i++ )
                    {    int jb = closures[i].isize( ) - j - 1;
                         seen.push_back( closures[i][jb] );    }
                    UniqueSort(seen);
                    if ( seen.solo( ) ) right.Set( mc - j - 1, seen[0] );
                    else
                    {    right.SetToSubOf( right, mc - j, j );
                         break;    }    }

               if ( logc.verb[ "CORRECT1" ] >= 1 )
               {    vec<placementy> p1, p2;
                    if ( left.size( ) > 0 )
                         p1 = log_control.FindGenomicPlacements(left);
                    if ( right.size( ) > 0 )
                         p2 = log_control.FindGenomicPlacements(right);
                    if ( ( left.size( ) > 0 && p1.empty( ) ) 
                         || ( right.size( ) > 0 && p2.empty( ) ) )
                    {
                         #pragma omp critical
                         {    cout << "\n";
                              PRINT2( id1, id1p );
                              PRINT2( con1.size( ), con2.size( ) );
                              PRINT2( left.size( ), right.size( ) );
                              cout << "final offsets = " << printSeq(final_offsets) 
                                   << endl;
     
                              vec<int> perfs;
                              for ( int i = 0; i < closures.isize( ); i++ )
                              {    cout << "\nclosure " << i;
                                   cout << "\noffset = " << closureso[i] << endl;
                                   vec<placementy> p = log_control.
                                        FindGenomicPlacements(closures[i]);
                                   if ( p.nonempty( ) ) cout << "perfect!" << endl;
                                   else cout << "imperfect" << endl;
                                   if ( p.nonempty( ) ) perfs.push_back(i);    
                                   const readstack& s = stacks[i];
                                   PRINT( s.Cols( ) );
                                   int totalq = 0;
                                   for ( int j = 0; j < s.Rows( ); j++ )
                                   for ( int l = 0; l < s.Cols( ); l++ )
                                   {    if ( s.Qual(j,l) > 2 ) 
                                        {    if ( s.Base(j,l) == closures[i][l] )
                                                  totalq += s.Qual(j,l);
                                             else totalq -= s.Qual(j,l);    }    }
                                   cout << "mean column qsum = "
                                        << double(totalq) / double( s.Cols( ) ) 
                                        << endl;    }
                              cout << "perfect ids: " << printSeq(perfs) << endl;
                              cout << "imperfect closure for " << id1 << ":" << endl;
                              static int count(0);
                              if ( count < 100 )
                              {    vec<look_align> aligns;
                                   ClusterAligner(left, *(log_control.G), 
                                        log_control.LG, *(log_control.Glocs), aligns);
                                   for ( int j = 0; j < aligns.isize( ); j++ )
                                   {    int g = aligns[j].target_id;
                                        aligns[j].PrintReadableBrief( cout, 
                                             ToString(id1), ToString(g) );
                                        aligns[j].PrintVisual( cout, left,
                                             (*(log_control.G))[g] );    }
                                   ClusterAligner(right, *(log_control.G), 
                                        log_control.LG, *(log_control.Glocs), aligns);
                                   for ( int j = 0; j < aligns.isize( ); j++ )
                                   {    int g = aligns[j].target_id;
                                        aligns[j].PrintReadableBrief( cout, 
                                             ToString(id1p), ToString(g) );
                                        aligns[j].PrintVisual( cout, right,
                                             (*(log_control.G))[g] );    }
                                   cout << "\nstack 1:\n";
                                   stack1.Print(cout);
                                   cout << "\nstack 2:\n";
                                   stack2.Print(cout);
                                   cout << "\nzero stack:\n";
                                   stacks[0].Print(cout);
                                   count++;    }    }    }    }

               corrected[id1] = left;
               #pragma omp critical
               {    count_closed++;    }
               REPORT_TIMEX( fclock, "used in closures, pair corr", out );
               #pragma omp critical
               {    cout << out.str( );    }
               if ( left != right )
               {    right.ReverseComplement( );
                    corrected[id1p] = right;    }    }    }

     // Done.

     if (logc.STATUS_LOGGING)
     {    cout << Date( ) << ": " << PERCENT_RATIO( 3, count_closed, use.isize( ) )
               << " of pairs attempted were closed" << endl;
          if ( closed > 0 )
          {    cout << Date( ) << ": " << PERCENT_RATIO( 3, closed_uniquely, closed )
                    << " of closed pairs had a unique closure" << endl;    }    }    }
