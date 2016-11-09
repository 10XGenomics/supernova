///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Given 'reads', find the friends of read ec, and kmer alignments of the gant
// to that read.

#ifndef GET_FRIENDS_AND_ALIGNS_INITIAL_H
#define GET_FRIENDS_AND_ALIGNS_INITIAL_H

#include "Basevector.h"
#include "CoreTools.h"
#include "kmers/KmerRecord.h"
#include "paths/long/ultra/ConsensusScoreModel.h"
#include "paths/long/Friends.h"
#include "paths/long/Logging.h"
#include "paths/long/LongProtoTools.h"

template<int K> void GetFriendsAndAlignsInitial( 
     const IAndOsVec& F,                  // given friends,
     const vecbasevector& reads,           // reads
     int ec,                               // where founder was
     vecbasevector& gang,                  // the gang; read 0 = founder
     vec<int> * p_friend_ids,               // optionally output friend ids in the gang
     vec< vec< pair<int,int> > >& a,       // kmer aligns of gang to read 0
     ostream& out,                         // for logging
     const ConsensusScoreModel& error_model,
     const long_heuristics& heur,
     const long_logging_control& log_control,
     const long_logging& logc );

#endif
