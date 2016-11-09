///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "CoreTools.h"
#include "feudal/PQVec.h"
#include "math/Functions.h"
#include "paths/long/large/DiscoStats.h"

void disco_stats::Compute( const String& work_dir )
{    int64_t nreads =
          MastervecFileObjectCount( work_dir + "/data/frag_reads_orig.fastb" );
     vec<int64_t> sample;
     srandomx(54321237);
     for ( int i = 0; i < Min( (int64_t) 10000, nreads ); i++ )
          sample.push_back( big_random( ) % nreads );
     vecbasevector some_reads;
     VecPQVec some_quals;
     some_reads.Read( work_dir + "/data/frag_reads_orig.fastb", sample );
     some_quals.Read( work_dir + "/data/frag_reads_orig.qualp", sample );
     some_reads.WriteAll( work_dir + "/sample_reads.fastb" );
     some_quals.WriteAll( work_dir + "/sample_reads.qualp" );
     BinaryWriter::writeFile( work_dir + "/sample_reads.ids", sample );
     vec<int> rlens, rquals;
     for ( int i = 0; i < (int) some_reads.size( ); i++ )
     {    rlens.push_back( some_reads[i].size( ) );
          qvec q;
          some_quals[i].unpack(&q);
          for ( int j = 0; j < some_reads[i].isize( ); j++ )
               rquals.push_back( q[j] );    }
     mean_read = Mean(rlens);
     mean_qual = Mean(rquals);
     total_bases = mean_read * nreads;    }
