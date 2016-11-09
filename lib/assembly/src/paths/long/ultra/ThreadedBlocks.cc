///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "CoreTools.h"
#include "PrintAlignment.h"
#include "VecUtilities.h"
#include "efasta/EfastaTools.h"
#include "pairwise_aligners/SmithWatFree.h"
#include "paths/BigMapTools.h"
#include "paths/LongReadTools.h"
#include "paths/long/ultra/ConsensusScoreModel.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/ultra/MultipleAligner.h"
#include "paths/long/ultra/ThreadedBlocks.h"
#include "polymorphism/Edit.h"

void threaded_blocks::PrintThreadSummary( const int g, ostream& out,
     const vecbasevector& genome, const int LG, 
     const VecIntPairVec& Glocs,
     const vec< vec<basevector> >& gap_truth ) const
{
     // Summarize threads.

     Bool found_truth = False;
     vec<String> threads;
     for ( int id = 0; id < NReads( ); id++ )
     {    if ( !Alive(id) ) continue;
          if ( !Member( ThreadRange(id), g ) ) continue;
          threads.push_back( Thread( id, g ).ToString( ) );    }
     Sort(threads);
     int n = 0;
     for ( int i = 0; i < threads.isize( ); i++ )
          n += threads[i].size( );
     out << "\nthreads for gap " << g+1 << " (" << threads.size( );
     if ( threads.nonempty( ) ) out << " of mean size " << n / threads.size( );
     out << "):\n";
     for ( int i = 0; i < threads.isize( ); i++ )
     {    int j = threads.NextDiff(i);
          out << j-i << " x " << threads[i];
          if ( genome.size( ) > 0 && gap_truth[g].solo( ) 
               && gap_truth[g][0].ToString( ) == threads[i] )
          {    out << " (TRUE)";
               found_truth = True;    }
          out << "\n";
          i = j - 1;    }    

     // Print truth.

     if ( genome.size( ) > 0 && !found_truth )
     {    if ( gap_truth[g].empty( ) ) 
               out << "failed to find true sequence for gap" << endl;
          else if ( !gap_truth[g].solo( ) ) 
          {    out << "found multiple sequences for gap" << endl;    }
          else if ( gap_truth[g][0].size( ) > 1000 )
          {    out << "true gap size appears to be " << gap_truth[g][0].size( )
                    << ", which seems fishy" << endl;    }
          else
          {    out << "true sequence for gap =\n";
               gap_truth[g][0].Print(out);    }    }    }

void FindEdits( const vec<basevector>& threads, const basevector& t,
     vec< vec<edit0> >& edits )
{    edits.clear_and_resize( t.size( ) + 1 );
     for ( int i = 1; i < threads.isize( ); i++ )
     {    if ( threads[i].size( ) == 0 && t.size( ) > 0 )
               edits[0].push( DELETION, t.isize( ) );
          if ( threads[i].size( ) > 0 && t.size( ) == 0 )
               edits[0].push( INSERTION, threads[i].ToString( ) );
          if ( threads[i].size( ) == 0 || t.size( ) == 0 ) continue;
          align a;

          /*
          int errors;
          const int bw_add = 5;
          int offset = ( threads[i].isize( ) - t.isize( ) ) / 2;
          int bandwidth = Abs(offset) + bw_add;
          SmithWatBandedA( threads[i], t, offset, bandwidth, a, errors, 0, 1, 1 );
          */

          // Note that in place of the following, we might want a banded
          // Smith-Waterman that has end penalties.  

          SmithWatFreeSym( threads[i], t, a, true, true, 1, 1, 1 );

          if ( a.pos1( ) > 0 )
          {    basevector b( threads[i], 0, a.pos1( ) );
               edits[0].push( INSERTION, b.ToString( ) );    }
          if ( a.pos2( ) > 0 ) edits[0].push( DELETION, a.pos2( ) );
          int p1 = a.pos1( ), p2 = a.pos2( );
          for ( int j = 0; j < a.Nblocks( ); j++ ) 
          {    if ( a.Gaps(j) > 0 )
               {    edits[p2].push( DELETION, a.Gaps(j) );
                    p2 += a.Gaps(j);    }
               if ( a.Gaps(j) < 0 )
               {    basevector b( threads[i], p1, -a.Gaps(j) );
                    edits[p2].push( INSERTION, b.ToString( ) );
                    p1 -= a.Gaps(j);    }
               for ( int x = 0; x < a.Lengths(j); x++ ) 
               {    if ( threads[i][p1] != t[p2] )
                    {    edits[p2].push( SUBSTITUTION, 
                              (char) as_base( threads[i][p1] ) );    }
                    p1++, p2++;    }    }    
          if ( a.Pos1( ) < threads[i].isize( ) )
          {    basevector b( threads[i], a.Pos1( ), 
                    threads[i].isize( ) - a.Pos1( ) );
               edits[ t.size( ) ].push( INSERTION, b.ToString( ) );    }
          if ( a.Pos2( ) < t.isize( ) )
               edits[ a.Pos2( ) ].push( DELETION, t.isize( ) - a.Pos2( ) );    }
     vec< pair<int,edit0> > keepers;
     for ( int p = 0; p <= t.isize( ); p++ )
          Sort( edits[p] );
     return;    }

void MakeEdits( const vec<basevector>& threads, basevector& t, 
     const vec< vec<edit0> >& edits )
{    vec< pair<int,edit0> > keepers;
     for ( int p = 0; p <= t.isize( ); p++ )
     {    for ( int i = 0; i < edits[p].isize( ); i++ )
          {    int j = edits[p].NextDiff(i);
               if ( 2*(j-i) > threads.isize( ) ) keepers.push( p, edits[p][i] );
               i = j - 1;    }    }
     for ( int j = keepers.isize( ) - 1; j >= 0; j-- )
     {    if ( j < keepers.isize( ) - 1 && keepers[j].first == keepers[j+1].first )
               continue;
          int p = keepers[j].first;
          const edit0& e = keepers[j].second;
          if ( e.etype == INSERTION )
          {    basevector b1( t, 0, p );
               basevector b2( e.seq );
               basevector b3( t, p, t.size( ) - p );
               t = Cat( b1, b2, b3 );    }
          if ( e.etype == DELETION )
          {    basevector b1( t, 0, p );
               basevector b2( t, p + e.n, t.isize( ) - ( p + e.n ) );
               t = Cat( b1, b2 );    }
          if ( e.etype == SUBSTITUTION ) t.Set( p, as_char( e.seq[0] ) );    }    }


// Improve the consensus guided by a scoring function that measures the (log)
// probability of generate the observed threads from the consensus. 
// Method: 
// 1. Score the consensus, and generate edits sugguested from alignment to other threads.
// 2. Starting from the most popular edits, iterative check all edits until find score improvements. 
// 3. Stop when not improves are found. Otherwise update the new consensus, and repeat from step 1.
void ImproveConsensus( const ConsensusScoreModel& model, const vec<basevector>& threads, basevector& t )
{    
     bool Use_Fast_Scorer = true; 
     const unsigned int Min_Edit_Votes = 2;
     int NIter = 10; // maximum iteration
     for( int iter=0; iter < NIter; iter++ ) {
         vec< pair<int, basevector> > candidates; // (score, new_consensus) pair
         // Use Fast Scoer for long sequence
         if ( Use_Fast_Scorer && t.size() >= 50 ) {
             FastScorer fs( model, t, threads );
             int score_new = fs.GetScore();
             vec< pair<int, edit0> > loc_edits;
             vec<unsigned int> counts;
             fs.GetLocEdits( &loc_edits, &counts, Min_Edit_Votes );
             int score = score_new;
             for ( size_t i = 0; i < loc_edits.size(); ++i ) {
                 basevector t2 = FastScorer::NewSeq( t, loc_edits[i] );
                 int score2 = fs.ScoreEdit( loc_edits[i] );
                 if ( score2 < score ) {
                     candidates.push( score2, t2 );
                     break;
                 } 
             }
         }
         if ( t.size() < 50 ) {
             vec< pair<int, edit0> > loc_edits;
             int score = model.Score( t, threads, Min_Edit_Votes,  &loc_edits );
             for ( size_t i = 0; i < loc_edits.size(); ++i ) {
                 basevector t2 = FastScorer::NewSeq( t, loc_edits[i] );
                 int score2 = model.Score( t2, threads );
                 if ( score2 < score ) {
                     candidates.push( score2, t2 );
                     break;
                 } 
             }
         }
         if ( candidates.empty() ) break;
         Sort( candidates );
         t = candidates.front().second;
     }
}

// Derived from MakeEdit, but allow votes from neighbor sites to be combined.
void MakeEdits2( const vec<basevector>& threads, basevector& t, 
     const vec< vec<edit0> >& edits )
{    vec< pair<int,edit0> > keepers;
     for ( int p = 0; p <= t.isize( ); p++ ) {    
          // the number of votes also includes the neighboring bases
          int n_votes = edits[p].isize( );
          if ( p > 0 ) n_votes += edits[p-1].size();
          if ( p < t.isize() ) n_votes += edits[p+1].size();
          if ( n_votes * 2 <= threads.isize() ) continue;
          int max_vote_i = -1, max_vote_p = -1, max_vote = 0;
          for ( int pp = p-1; pp <= p+1; pp++ ) {    
               if ( pp < 0 || pp >= edits.isize() ) continue;
               for ( int i = 0; i < edits[pp].isize( ); i++ )
               {    int j = edits[pp].NextDiff(i);
                    if ( (j-i) > max_vote ) {
                         max_vote = j - i;
                         max_vote_p = pp;
                         max_vote_i = i;
                    }
                    i = j - 1;
               }    
          }
          // Only if the most popular edit is at base p and the base be not modified previously.
          if ( max_vote_p == p && ( keepers.empty() || keepers.back().first != p ) )
               keepers.push( p, edits[p][max_vote_i] );
     }
     for ( int j = keepers.isize( ) - 1; j >= 0; j-- )
     {    if ( j < keepers.isize( ) - 1 && keepers[j].first == keepers[j+1].first )
               continue;
          int p = keepers[j].first;
          const edit0& e = keepers[j].second;
          if ( e.etype == INSERTION )
          {    basevector b1( t, 0, p );
               basevector b2( e.seq );
               basevector b3( t, p, t.size( ) - p );
               t = Cat( b1, b2, b3 );    }
          if ( e.etype == DELETION )
          {    basevector b1( t, 0, p );
               if ( p + e.n <= t.isize( ) )
               {    basevector b2( t, p + e.n, t.isize( ) - ( p + e.n ) );
                    t = Cat( b1, b2 );    }    }
          if ( e.etype == SUBSTITUTION ) t.Set( p, as_char( e.seq[0] ) );    }    }

// Check the homopolymers in the consensus sequence and determine the most likely length. 
// If cannot determine, output several possible sequences.
void ImproveHomopolymer( const ConsensusScoreModel& model, const vec<basevector>& threads, 
          const basevector& left_ext, const basevector& right_ext,
          vec< basevector>& consensus ) 
{
     const double prob_ratio_cutoff = 0.1; // any sequence with probability less than the most likely sequence will be deleted
     vec<double> logps( consensus.size() );
     for ( size_t i = 0; i < consensus.size(); ++i ) {
          // extent
          BaseVec extended_con = Cat( left_ext, consensus[i], right_ext );
          double logp = 0; // the log likelihood of this consensus sequence.
          for ( size_t j = 0; j < threads.size(); ++j ) {
               BaseVec extended_thread = Cat( left_ext, threads[j], right_ext );
               double p = model.Probability( extended_con, extended_thread );
               logp += log(p);
          }
          logps[i] = logp;
     }
     SortSync( logps, consensus );
     vec<Bool> to_del( consensus.size(), False );
     double logp_cutoff = logps.back() + log( prob_ratio_cutoff );  
     for ( size_t i = 0; i < logps.size(); ++i ) 
          if ( logps[i] < logp_cutoff ) to_del[i] = True;
     EraseIf( consensus, to_del );
     ForceAssert( ! consensus.empty() );
}

void threaded_blocks::ThreadConsensus( const int g, const Scorer& scorer,
     const vec<basevector>& gap_truth, vec<basevector>& consensus, 
     String& report, const ConsensusScoreModel& error_model, 
     ostream& xout, const long_heuristics& heur,
     const long_logging_control& log_control, const long_logging& logc ) const
{    
     double clock1 = WallClockTime( );
     vec<basevector> threads;
     vec<int> thread_read_ids; // which read the thread is from
     for ( int id = 0; id < NReads( ); id++ )
     {    if ( !Alive(id) ) continue;
          if ( !Member( ThreadRange(id), g ) ) continue;
          threads.push_back( Thread( id, g ) );    
          thread_read_ids.push_back( id );    }
     if ( threads.empty( ) )
     {    if (logc.PRINT_EDITS) report = "No threads, giving up.\n";
          return;    }
     ostringstream out;
     vec< vec<edit0> > edits;
     FindEdits( threads, threads[0], edits );
     basevector t = threads[0];
     MakeEdits( threads, t, edits );
     REPORT_TIMEX( clock1, "used in consensus 1", xout );

     // Second iteration, allowing neighbor votes to be combined.

     double clock2 = WallClockTime( );
     FindEdits( threads, t, edits );
     MakeEdits2( threads, t, edits );
     REPORT_TIMEX( clock2, "used in consensus 2", xout );

     // Third iteration, improve the consensus guided by a scoring 
     // function that measures the (log) probability of generate the 
     // observed threads from the consensus:
     // 1. Score the consensus, and generate edits suggested from alignment to other 
     //    threads.
     // 2. Starting from the most popular edits, iterative check all edits until 
     //    find score improvements. 
     // 3. Stop when not improves are found. Otherwise update the new consensus, 
     //    and repeat from step 1.
     // 4. Optionally output the alternatives edits and read ids that supports the 
     //    alternatives   

     if ( ! error_model.Empty() ) 
     {    if ( !heur.USE_MULTIPLE_ALIGNER ) {
               double iclock = WallClockTime( );
               ImproveConsensus( error_model, threads, t );
               REPORT_TIMEX( iclock, "used in ImproveConsensus", xout );
          }
          else
          {    double mclock = WallClockTime( );
               MultipleAligner ma( scorer, t );
               ma.addReads(threads.begin(),threads.end());
               t = ma.getConsensus();    
               REPORT_TIMEX( mclock, "used in mult aligner", xout );    }    }

     double clock3b = WallClockTime( );
     if ( !heur.LENGTHEN_HOMOPOLYMERS ) consensus.push_back(t);
     else
     {    vec<basevector> c(1);
          for ( int i = 0; i < t.isize( ); i++ )
          {    int j;
               for ( j = i + 1; j < t.isize( ); j++ )
                    if ( t[j] != t[i] ) break;
               const int trig = 20;
               const int trig2 = 30;
               const basevector &b1 = Block(g), &b2 = Block(g+1);

               // First look for long dinucleotide runs.

               if ( heur.LENGTHEN_DINUKES && j - i == 1 && j < t.isize( ) )
               {    int jsave = j;
                    for ( j = i + 2; j < t.isize( ) - 1; j += 2 )
                         if ( t[j] != t[i] || t[j+1] != t[i+1] ) break;
                    int count = j - i;
                    if ( i == 0 )
                    {    for ( int r = b1.isize( ) - 2; r >= 0; r -= 2 )
                         {    if ( b1[r] != t[i] || b1[r+1] != t[i+1] ) break;
                              count += 2;    }    }
                    if ( j == t.isize( ) )
                    {    for ( int r = 0; r < b2.isize( ) - 1; r += 2 )
                         {    if ( b2[r] != t[i] || b2[r+1] != t[i+1] ) break;
                              count += 2;    }    }
                    if ( j == t.isize( ) - 1 && t.back( ) == t[i] )
                    {    for ( int r = 0; r < b2.isize( ) - 1; r += 2 )
                         {    if ( b2[r] != t[i+1] || b2[r+1] != t[i] ) break;
                              count += 2;    }    }
                    if ( count >= trig2 )
                    {    for ( int l = 0; l < c.isize( ); l++ )
                         for ( int k = i; k < j; k++ )
                              c[l].push_back( t[k] );
                         // xout << "lengthening dinuke, count = " << count << "\n";
                         vec<basevector> cnew(c);
                         for ( int l = 0; l < c.isize( ); l++ )
                         {    c[l].push_back( t[i] );
                              c[l].push_back( t[i+1] );    }
                         cnew.append(c);
                         c = cnew;
                         i = j - 1;
                         continue;    }
                    j = jsave;    }

               // Now look for mononucleotide runs.

               for ( int l = 0; l < c.isize( ); l++ )
               for ( int k = i; k < j; k++ )
                    c[l].push_back( t[i] );
               int count = j - i;
               if ( i == 0 )
               {    for ( int r = b1.isize( ) - 1; r >= 0; r-- )
                    {    if ( b1[r] != t[i] ) break;
                         count++;    }    }
               if ( j == t.isize( ) )
               {    for ( int r = 0; r < b2.isize( ); r++ )
                    {    if ( b2[r] != t[i] ) break;
                         count++;    }    }
               if ( count >= trig )
               {    vec<basevector> cnew(c);
                    // xout << "lengthening homopolymer, count = " << count << "\n";
                    for ( int l = 0; l < c.isize( ); l++ )
                         c[l].push_back( t[i] );
                    cnew.append(c);
                    c = cnew;    }
               if ( heur.LENGTHEN_HOMOPOLYMERS2 && count >= 2*trig )
               {    vec<basevector> cnew(c);
                    for ( int l = 0; l < c.isize( ); l++ )
                         c[l].push_back( t[i] );
                    cnew.append(c);
                    c = cnew;    }
               i = j - 1;    }
          consensus.append(c);    }
     REPORT_TIMEX( clock3b, "used in consensus 3b", xout );

     double clock4 = WallClockTime( );
     if ( heur.IMPROVE_HOMOPOLYMERS && consensus.size() > 1 ) {
          // find the edge
          const basevector &b1 = Block(g), &b2 = Block(g+1);
          const basevector c =  consensus[0];
          int len = c.size();
          basevector right_ext;
          // mononucleotide extention
          for ( size_t i = 0; i < b2.size(); i++ ) 
               if ( b2[i] == c.back() ) right_ext.push_back( b2[i] );
               else break;
          // dinucleotide extention
          if ( right_ext.empty() && c.size() >=2 ) {
               for ( size_t i = 0; i < b2.size() - 1; i++ ) 
                    if ( b2[i] == c[len-2] && b2[i+1] == c[len-1] ){
                         right_ext.push_back( b2[i] );
                         right_ext.push_back( b2[i+1] );
                    }
                    else break;
          }
          basevector left_ext;
          for ( int i = b1.size()-1; i >=0;  i-- ) 
               if ( b1[i] == c.front() ) left_ext.push_back( b1[i] );
               else break;
          if ( left_ext.empty() && c.size() >=2 ) {
               int n = b1.size();
               int j = 0;
               for ( ; j < n -1; j += 2)
                    if ( ! ( b1[n - j - 2] == c[0] && b1[n - j -1] == c[1] ) ) break;
               if ( j > 0 ) left_ext.SetToSubOf( b1, b1.size() - j, j );
          }
          ImproveHomopolymer( error_model, threads, left_ext, right_ext, consensus );
     }
     REPORT_TIMEX( clock4, "used in consensus 4", xout );

     if (logc.PRINT_EDITS)
     {    //out << "consensus = " << t.ToString( );
          Bool simple = ( gap_truth.solo( ) && t == gap_truth[0] );
          if ( !simple ) 
          {    out << "\nedits for gap " << g+1 << " (" 
                    << threads.size( ) << " threads)" << endl;
               for ( int p = 0; p <= threads[0].isize( ); p++ )
               {    
                    // Work around for undiagnosed problem.

                    if ( p >= edits.isize( ) ) continue;

                    for ( int i = 0; i < edits[p].isize( ); i++ )
                    {    int j = edits[p].NextDiff(i);
                         const int min_to_print = 3;
                         if ( j-i >= min_to_print || 2*(j-i) > threads.isize( ) )
                         {    out << p << " --> " << j-i << " x ";
                              if ( edits[p][i].etype == INSERTION )
                                   out << "insertion of " << edits[p][i].seq;
                              if ( edits[p][i].etype == DELETION )
                                   out << "deletion of " << edits[p][i].n;
                              if ( edits[p][i].etype == SUBSTITUTION )
                                   out << "substitution of " << edits[p][i].seq;
                              if ( 2*(j-i) > threads.isize( ) ) out << " ***";
                              out << endl;    }
                         i = j - 1;    }    }
               out << "initial = " << threads[0].ToString( ) << endl;
               if (  !Alive(0) || !Member( ThreadRange(0), g ) )
               {    out << "don't have thread for read 0, "
                         << "so this is a bit fishy" << endl;    }
               out << "consensus = " << t.ToString( ) << endl;
               if ( gap_truth.empty( ) ) out << "No true sequence for gap.\n";
               if ( gap_truth.solo( ) )
               {    const basevector& gt = gap_truth[0];
                    out << "truth     = " << gt.ToString( ) << "\n\n";    
                    if ( t.size( ) > 0 && gt.size( ) > 0 )
                    {    align a;

                         // Print the scores.

                         out << "alignment of consensus to truth (scores = "
                              << error_model.Score( t, threads ) << ", "
                                   << error_model.Score( gt, threads ) 
                                   << ", respectively)" << endl;
                         SmithWatFreeSym( t, gt, a, true, true, 1, 1, 1 );
                         PrintVisualAlignment( True, out, t, gt, a );    
                         if ( a.pos1( ) > 0 )
                         {    out << "insertion of " << a.pos1( )
                                   << " on left" << endl;    }
                         if ( a.pos2( ) > 0 ) 
                         {    out << "deletion of " << a.pos1( )
                                   << " on left" << endl;    }
                         if ( a.Pos1( ) < t.isize( ) )
                         {    out << "insertion of " << t.isize( ) - a.Pos1( )
                                   << " on left" << endl;    }
                         if ( a.Pos2( ) < gt.isize( ) )
                         {    out << "insertion of " << gt.isize( ) - a.Pos2( )
                                   << " on left" << endl;    }    }    }
               if ( gap_truth.size( ) > 1 )
               {    out << "multiple truth values\n";
                    if ( gap_truth.size( ) <= 10 )
                    {    for ( int j = 0; j < gap_truth.isize( ); j++ )
                         {    const basevector& gt = gap_truth[j];
                              gt.Print( out, "truth." 
                                   + ToString(j+1) );    }    }    }    }
          report = out.str( );    }    }

efasta threaded_blocks::MakeCorrectedRead( const ConsensusScoreModel& error_model, 
     ostream& out, const long_heuristics& heur, 
     const long_logging_control& log_control, const long_logging& logc ) const
{    efasta r;
     if ( NBlocks( ) == 0 ) return r;
     double clock1 = WallClockTime( );
     vec< vec<basevector> > gap_truth;
     if (logc.PRINT_EDITS) 
     {    GetGapTruth( *log_control.G, log_control.LG, *log_control.Glocs, 
               gap_truth );    }
     vec< vec<basevector> > consensus( NGaps( ) );
     vec<String> reports( NGaps( ) );
     REPORT_TIMEX( clock1, "used in correction setup", out );

     // The model for scoring the consensus sequence. 

     double sclock = WallClockTime( );
     Scorer scorer( error_model.GetSubRate( ), error_model.GetDelRate( ), 
          error_model.GetInsRate( ) );
     REPORT_TIMEX( sclock, "used creating scorer", out );

     // Go through the gaps.

     for ( int g = 0; g < NGaps( ); g++ )
     {    vec<basevector> gt;
          if (logc.PRINT_EDITS) gt = gap_truth[g];
          ThreadConsensus( g, scorer, gt, consensus[g], reports[g], 
               error_model, out, heur, log_control, logc );    }

     // Finish.

     double fclock = WallClockTime( );
     if (logc.PRINT_EDITS)
     {    for ( int j = 0; j < reports.isize( ); j++ )
               out << reports[j];    }
     for ( int j = 0; j < NBlocks( ); j++ )
     {    r += Block(j).ToString( );
          if ( j < NGaps( ) )
          {    if ( consensus[j].empty( ) ) break;
               r += efasta( consensus[j] );    }    }
     REPORT_TIMEX( fclock, "used finishing consensus", out );
     return r;    }

void threaded_blocks::GetGapTruth( const vecbasevector& genome, const int LG, 
     const VecIntPairVec& Glocs,
     vec< vec<basevector> >& gap_truth ) const
{    if ( NBlocks() == 0 ) { gap_truth.clear(); return; }
     gap_truth.clear_and_resize( NGaps( ) );
     for ( int g = 0; g < NGaps( ); g++ )
     {    const basevector &b1 = Block(g), &b2 = Block(g+1);
          vec<placementy> places1 = FindGenomicPlacementsY(0, b1, LG, genome, Glocs);
          vec<placementy> places2 = FindGenomicPlacementsY(0, b2, LG, genome, Glocs);
          for ( int j1 = 0; j1 < places1.isize( ); j1++ )
          for ( int j2 = 0; j2 < places2.isize( ); j2++ )
          {    const placementy &p1 = places1[j1], &p2 = places2[j2];
               if ( p1.g != p2.g || p1.fw != p2.fw ) continue;
               basevector x;
               if ( p1.fw )
               {    if ( p1.Pos > p2.pos ) continue;
                    x.SetToSubOf( genome[p1.g], p1.Pos, p2.pos - p1.Pos );    }
               else
               {    if ( p2.Pos > p1.pos ) continue;
                    x.SetToSubOf( genome[p1.g], p2.Pos, p1.pos - p2.Pos );
                    x.ReverseComplement( );    }
               gap_truth[g].push_back(x);    }
          UniqueSort( gap_truth[g] );    }    }
