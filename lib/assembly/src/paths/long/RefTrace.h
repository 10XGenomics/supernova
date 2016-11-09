///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef REFTRACE_H
#define REFTRACE_H

#include "CoreTools.h"
#include "paths/HyperBasevector.h"

class RefTraceResults {

     public:

     RefTraceResults( ) : penalty(0), gaps(0), meta_events(0) { }
     RefTraceResults( const int penalty, const int gaps, const int meta_events ) 
          : penalty(penalty), gaps(gaps), meta_events(meta_events) { }

     int penalty;
     int gaps;
     int meta_events;

};

class RefTraceHeuristics {

     public:

     RefTraceHeuristics( ) : max_offset_diff(10), max_error_rate(0.05), 
          offset_add(1), max_twiddle(3), min_group_frac(0.75), 
          min_group_save(1000000000) { }

     RefTraceHeuristics( const int max_offset_diffx, const double max_error_ratex,
          const int offset_addx, const int max_twiddlex,
          const double min_group_fracx, const int min_group_savex )
          : max_offset_diff(max_offset_diffx), max_error_rate(max_error_ratex),
          offset_add(offset_addx), max_twiddle(max_twiddlex),
          min_group_frac(min_group_fracx), min_group_save(min_group_savex)
     {    RefTraceHeuristics rth;
          if ( max_offset_diff < 0 ) max_offset_diff = rth.max_offset_diff;
          if ( max_error_rate < 0 ) max_error_rate = rth.max_error_rate;
          if ( offset_add < 0 ) offset_add = rth.offset_add;
          if ( max_twiddle < 0 ) max_twiddle = rth.max_twiddle;
          if ( min_group_frac < 0 ) min_group_frac = rth.min_group_frac;     
          if ( min_group_save < 0 ) min_group_save = rth.min_group_save;    }

     int max_offset_diff;
     double max_error_rate;
     int offset_add;
     int max_twiddle;
     double min_group_frac;
     int min_group_save;

};


// BEST_GLOBAL_OUT: if provided, name of fastb file to write two basevectors to:
// first the best assembly sequence, then the best reference sequence.

class SupportedHyperBasevector; 
class ReadOriginTracker;
class ref_data;
class RefTraceControl;
class long_logging;

RefTraceResults RefTrace( const ref_data& ref,
     const HyperBasevector& hb, const vec<int>& inv, 
     const int verbosity, const long_logging& logc, 
     ostream& out, RefTraceHeuristics rth,
     const String& BEST_GLOBAL_OUT, const Bool fix_bug);

RefTraceResults RefTraceAndCallVaraint( const ref_data& ref,
     const HyperBasevector& hb, const vec<int>& inv, 
     const int verbosity, const long_logging& logc, 
     ostream& out, ostream& callsOut, RefTraceHeuristics rth,
     const String& BEST_GLOBAL_OUT, const Bool fix_bug,
     const RefTraceControl& ref_trace_control, 
     const ReadOriginTracker* p_read_tracker = NULL);

#endif
