///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef DISCO_STATS_H
#define DISCO_STATS_H

#include "CoreTools.h"

class disco_stats {

     public:

     double mean_read; // mean read length
     double mean_qual; // mean base quality
     double total_bases; // total bases in reads

     void Compute( const String& work_dir );

};

TRIVIALLY_SERIALIZABLE(disco_stats);

#endif
