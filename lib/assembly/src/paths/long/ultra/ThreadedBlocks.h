///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// This defines a class for a "system of threaded blocks".
//
// blocks_       = the b sequences that are the blocks in the n reads
// threads_      = an n x (b-1) matrix containing the threads
// alive_        = is read alive?
// thread_range_ = start and stop positions for threads for a given live read

#ifndef THREADED_BLOCKS_H
#define THREADED_BLOCKS_H

#include "CoreTools.h"
#include "IntPairVec.h"
#include "efasta/EfastaTools.h"
#include "math/HoInterval.h"
#include "paths/long/ultra/ConsensusScoreModel.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/ultra/MultipleAligner.h"

class threaded_blocks {

     public:

     threaded_blocks( ) { }
     threaded_blocks( const vec<basevector>& blocks, 
          const vec< vec<basevector> >& threads, const vec<Bool>& alive, 
          const vec<ho_interval>& thread_range ) : blocks_(blocks), 
               threads_(threads), alive_(alive), thread_range_(thread_range) { }

     int NBlocks( ) const { return blocks_.size( ); }
     int NGaps( ) const { return blocks_.isize( ) - 1; }
     int NReads( ) const { return threads_.size( ); }
     const basevector& Block( int b ) const { return blocks_[b]; }
     Bool Alive( const int id ) const { return alive_[id]; }
     const ho_interval& ThreadRange( int id ) const { return thread_range_[id]; }
     const basevector& Thread( const int id, const int b ) const
     {    Assert( Alive(id) );
          Assert( Member( ThreadRange(id), b ) );
          return threads_[id][b];    }
     basevector& Thread( const int id, const int b )
     {    Assert( Alive(id) );
          Assert( Member( ThreadRange(id), b ) );
          return threads_[id][b];    }

     void ThreadConsensus( const int g, const Scorer& scorer,
          const vec<basevector>& gap_truth, vec<basevector>& consensus, 
          String& report, const ConsensusScoreModel& error_model, ostream& xout, 
          const long_heuristics& heur, const long_logging_control& log_control,
          const long_logging& logc ) const;

     // Build a 'corrected' read by choosing the most common threads, or by
     // computing consensus, in one of two ways.

     efasta MakeCorrectedRead( const ConsensusScoreModel& error_model, ostream& out,
          const long_heuristics& heur, const long_logging_control& log_control,
          const long_logging& logc ) const ;

     // Get true sequence for gaps.

     void GetGapTruth( const vecbasevector& genome, const int LG,
          const VecIntPairVec& Glocs,
          vec< vec<basevector> >& gap_truth ) const;

     // Print out the threads, collecting identical ones together in a single
     // line.  If genome provided, attempt to print true value for the sequence
     // that should go between the blocks.

     void PrintThreadSummary( const int g, ostream& out,
          const vecbasevector& genome, const int LG, 
          const VecIntPairVec& Glocs,
          const vec< vec<basevector> >& gap_truth ) const;

     void PrintEdits( ) const;

     private:

     vec<basevector> blocks_;
     vec< vec<basevector> > threads_;
     vec<Bool> alive_;
     vec<ho_interval> thread_range_;

};

#endif
