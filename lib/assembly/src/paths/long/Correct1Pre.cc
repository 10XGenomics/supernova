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
#include "Qualvector.h"
#include "paths/BigMapTools.h"
#include "paths/long/Correct1Pre.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/FriendAligns.h"
#include "paths/long/MakeKmerStuff.h"
#include "paths/long/ReadStack.h"
//#include "paths/long/ColorCompare.h"


namespace {
template <class veclike>
void append_read( const std::string& filename, size_t serial, const veclike& r) 
{
     ofstream out( filename.c_str(), ios_base::app );
     ForceAssert(out.good());

     out << ">" << serial << endl;
     for ( size_t i = 0; i < r.size(); ++i ) {
	    if ( i && (i%80)==0 )
		       out << endl;
		out << Base::val2Char(r[i]);
		  }
     out << endl;
       ForceAssert(out.good());
	 out.close();
     }
}

void Correct1Pre( String const& tmpDir, const int K, const int max_freq, 
     vecbasevector& bases, vecqualvector& quals, const PairsManager& pairs, 
     const vec<Bool>& to_edit, vec<int>& trim_to, const vec<int>& trace_ids, 
     const long_logging& logc, const long_heuristics& heur )
{
     double clock1 = WallClockTime( );
     if (logc.STATUS_LOGGING) ReportPeakMem( "begin Correct1Pre" );

     // Build alignments.

     FriendAligner faligns(bases,quals,to_edit,
                             tmpDir+"/correct1pre."+ToString(K)+".friends",
                             heur.FF_MAKE_ALIGN_IMPL, K,
                             heur.FF_MIN_FREQ,heur.FF_MAX_FREQ,
                             heur.FF_MIN_QUAL,heur.FF_COVERAGE,
                             heur.FF_DOWN_SAMPLE,heur.FF_VERBOSITY);
     if (logc.STATUS_LOGGING) ReportPeakMem( "alignment data created" );

     //// Define read starts.

     //vec<int64_t> id1_start( bases.size( ) + 1, -1 );
     //id1_start[ bases.size( ) ] = aligns_all.size( );
     //for ( int64_t j = aligns_all.jsize( ) - 1; j >= 0; j-- )
     //     id1_start[ aligns_all[j].id1 ] = j;
     //for ( int id = (int) bases.size( ) - 1; id >= 0; id-- )
     //     if ( id1_start[id] < 0 ) id1_start[id] = id1_start[id+1];
     //cout << Date( ) << ": read starts defined" << endl;

     // Go through the reads.

     trim_to.resize( bases.size( ) );

     vec<int64_t> use;
     for ( int64_t id1 = 0; id1 < (int64_t) bases.size( ); id1++ )
          if ( to_edit[id1] ) use.push_back(id1);

//     vecbasevector bases_new(bases);
//     vecqualvector quals_new(quals);

     //output buffer, of size use.size(), since base are 2 bits, 4 values,
     // there are not invalid value, so might as well use arbitrary but trackable
     // bases[0] and quals[0] for initial value
     vecbasevector bases_loc(use.size(),bases[0]);
     vecqualvector quals_loc(use.size(),quals[0]);


     // Do the corrections.

     if (logc.STATUS_LOGGING)
     {    cout << Date( ) << ": going through " << ToStringAddCommas( use.size( ) )
               << " reads" << endl;    }
     const int batch = 100;
     REPORT_TIME( clock1, "used in Correct1Pre head" );
     double clock2 = WallClockTime( );
     #pragma omp parallel for schedule(dynamic, 1)
     for ( int64_t id1a = 0; id1a < (int64_t) use.size( ); id1a += batch )
     {    readstack stack;
          vec<Bool> suspect;
          for ( int64_t id1x = id1a; 
               id1x < Min( id1a + batch, (int64_t) use.size( ) ); id1x++ )
          {    int64_t id1 = use[id1x];
               BaseVec& currentBaseVec = bases_loc[id1x];
               UCharVec& currentCharVec = quals_loc[id1x];
               currentBaseVec = bases[id1];
               currentCharVec = quals[id1];

               trim_to[id1] = bases[id1].size( );
               if ( bases[id1].empty( ) ) continue;
               Friends aligns;
               faligns.getAligns( id1, &aligns );

               // we could win more here by not returning > MAX_STACK aligns
               if ( aligns.size() > static_cast<unsigned>(heur.MAX_STACK) )
               {    //bases_new[id1] = bases[id1];
                    //quals_new[id1] = quals[id1];
                    continue; }

               stack.Initialize( id1, aligns, 0, aligns.size(),
                    readstack::strict, bases, quals, pairs );

               if ( BinMember( trace_ids, id1 ) )
               {    vec<basevector> cons;
                    cons.push_back( currentBaseVec );
                    #pragma omp critical
                    {    cout << "\nCorrect1Pre, K = " << K 
                              << ", initial stack, tracing read " << id1 << endl;
                         cout << "stack:\n";
                         stack.Print( cout, cons );    }    }

#if 0
               if ( stack.Rows( ) > heur.MAX_STACK )
               { //   bases_new[id1] = bases[id1];
                 //   quals_new[id1] = quals[id1];
                    continue; }                 // this continue was missing
#endif

               stack.HighQualDiff( 30, 1, suspect );
               stack.Erase(suspect);
               if ( heur.HQ_DIFF_WINDOW )
               {    stack.HighQualDiffWindow( suspect );
                    stack.Erase(suspect);    }
               if ( ! heur.USE_EM )
               {    Bool verbose = False;
                    stack.CorrectAll( 
                         currentBaseVec, currentCharVec, trim_to[id1], 
                         verbose );    }
               if ( BinMember( trace_ids, id1 ) )
               {    vec<basevector> cons;
                    cons.push_back( currentBaseVec );
                    #pragma omp critical
                    {    cout << "\nCorrect1Pre, K = " << K 
                              << ", final stack, tracing read " << id1 << endl;
                         cout << "stack:\n";
                         stack.Print( cout, cons );    }    }    }    }

     // Save results.

     #pragma omp parallel for schedule(dynamic)
     for( size_t id1x = 0 ; id1x < use.size() ; ++id1x){
          const int64_t& id1 = use[id1x];
          bases[id1] = bases_loc[id1x];
          quals[id1] = quals_loc[id1x];
     }
     REPORT_TIME( clock2, "used in Correct1Pre main" );    }
