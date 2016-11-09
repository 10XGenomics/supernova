///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <sys/utsname.h>
#include <sys/wait.h>

#include "CoreTools.h"
#include "FastIfstream.h"
#include "paths/long/fosmid/Fosmids.h"
#include "paths/long/large/DiscoStats.h"
#include "paths/long/large/GapToyTools.h"

void ReportMemory( const disco_stats& stats )
{
     // Report memory.

     uint64_t phys_mem = physicalMemory( );
     cout << Date( ) << ": see total physical memory of "
          << ToStringAddCommas(phys_mem) << " bytes" << endl;
     rlimit max_data;
     if ( getrlimit( RLIMIT_DATA, &max_data ) == 0 )
     {    uint64_t max_mem = max_data.rlim_max;
          if ( max_mem < phys_mem )
          {    cout << "\nWell this is very interesting.  Apparently "
                    << "your memory usage is capped at\n" 
                    << ToStringAddCommas(max_mem)
                    << ".  This is less than the physical memory on\n" 
                    << "your machine, and may result in a crash.  "
                    << "Please let us know\nif this happens, as we can "
                    << "make our code respect the memory cap.\n\n";    }    }
     uint64_t allowed_mem = GetMaxMemory( );
     if ( allowed_mem > 0 && allowed_mem < phys_mem )
     {    cout << Date( ) << ": see user-imposed limit on memory of "
               << ToStringAddCommas(allowed_mem) << " bytes" << endl;    }

     // Report bytes per base and issue warning if appropriate.

     int64_t total_bytes = GetMaxMemory( );
     double bytes_per_base = total_bytes / stats.total_bases;
     cout << Date( ) << ": " << setiosflags(ios::fixed) 
          << setprecision(2) << bytes_per_base << resetiosflags(ios::fixed)
          << " bytes per read base, assuming max memory available" << endl;
     // 2.25 OK in 51452.YRI
     // 1.85 failed in 51454.F3
     if ( bytes_per_base < 1.8 )
     {    cout << "\nWARNING: generally 2.0 bytes per read base of memory are "
               << "needed.  You have\nsubstantially less than this, so the "
               << "odds of your assembly completing are low.\n" << endl;    }
     if ( bytes_per_base < 2.0 )
     {    cout << "\nWARNING: generally 2.0 bytes per read base of memory are "
               << "needed.  You have\nsomewhat less than this, so it is "
               << "possible that your assembly will not finish.\n" << endl;    }
     if ( bytes_per_base < 2.4 )
     {    cout << "\nWARNING: generally about 2.0 bytes per read base is enough, "
               << "but we have seen\ncases with up to 2.32 bytes per read base "
               << "where memory is exhausted.\n" << endl;    }    }

void PrintSysInfo( )
{
     // Print system info.

     cout << "SYSTEM INFO" << endl;
     utsname un;
     int ok = uname( &un );
     if ( ok != 0 )
     {    cout << "- Warning: attempt to determine system name through uname "
               << "failed." << endl;
          cout << "- This is weird but not a problem." << endl;    }
     else
     {    cout << "- OS: " << un.sysname << " :: " << un.release
               << " :: " << un.version << endl;
          cout << "- node name: " << un.nodename << endl;
          cout << "- hardware type: " << un.machine << endl;    }
     if ( IsRegularFile( "/proc/cpuinfo" ) )
     {    fast_ifstream in( "/proc/cpuinfo" );
          String line;
          vec<String> lines;
          while(1)
          {    getline( in, line );
               if ( in.fail( ) ) break;
               line.GlobalReplaceBy( "\t", "" );
               lines.push_back(line);    }
          UniqueSort(lines);
          for ( int i = 0; i < lines.isize( ); i++ )
          {    const String& L = lines[i];
               if ( L.Contains("model name", 0) ) cout << "- cpu " << L << endl;
               if ( L.Contains("cpu MHz", 0) || L.Contains("cache size", 0) )
                    cout << "- " << L << endl;    }    }
     cout << "- physical memory: " <<  PhysicalMemGBString( ) << endl << endl;    }

void MemoryCheck( const Bool MEMORY_CHECK, const String& work_dir )
{    if (MEMORY_CHECK)
     {    cout << "MEMORY CHECK (typically takes several minutes; could cause "
               << "\nmachine to become sluggish or result in this job being killed)" 
               << endl;
          Remove( work_dir + "/available_mem" );
          pid_t pid = fork( );
          vec<char*> ptr(100);
          if ( pid == -1 )
          {    cout << "- Warning: tried to fork process to test available memory, "
                    << "but the fork failed." << endl;    }
          else if ( pid == 0 )
          {    int64_t P = GetMaxMemory( ), d;
               int64_t M = P/100;
               for ( d = 1; d <= 100; d++ )
               {    // probably OK if executed before RunTime( ):
                    // ptr[d-1] = new (std::nothrow) char[M];
                    ptr[d-1] = (char*) malloc(M);
                    if ( ptr[d-1] == NULL )
                    {    d--;
                         cout << "- Apparently able to allocate " << d << "% of "
                              << "nominally available memory (and not more)." 
                              << endl;
                         cout << "- Can access " << (d*M) / (1024*1024*1024) 
                              << " GB." << endl;
                         if ( d < 90 )
                         {   cout << "- WARNING: you may have less available memory "
                                  << "than you think!" << endl;    }
                         else 
                         {    cout << "- But you're pretty close, so things should "
                                   << "be OK." << endl;    }
                         Ofstream( out, work_dir + "/available_mem" );
                         out << d*M << endl;
                         break;    }
                    #pragma omp parallel for
                    for ( int64_t i = 0; i < M/1024; i++ )
                         ptr[d-1][1024*i] = 0;    }
               if ( d > 100 )
               {    cout << "- Apparently able to allocate 100% of nominally "
                         << "available memory." << endl;
                    cout << "- Can access at least " << P / (1024*1024*1024) 
                         << " GB." << endl;    }
               _exit(0);    }
          else 
          {    int status;
               waitpid( pid, &status, 0 );    }
          if ( IsRegularFile( work_dir + "/available_mem" ) )
          {    Ifstream( memin, work_dir + "/available_mem" );
               int64_t mema;
               memin >> mema;
               SetMaxMemory(mema);
               cout << "- Lowering max memory to " << ToStringAddCommas(mema)
                    << " bytes." << endl;    }
          cout << endl;    }
     else
     {    cout << "Omitting memory check.  If you run into problems with memory,\n"
               << "you might try rerunning with MEMORY_CHECK=True." 
               << endl << endl;    }    }

void DefineRegions( const String& X, vec<int>& fosmids, vec<String>& regions,
     map<String,GapToyResults>& res, const int PAD, Bool& all, const String& SAMPLE,
     String& EVALUATE, const String& F )
{
     if ( X == "tiny" )
     {    fosmids = {1};
          regions.push_back( "1:24.8M-24.9M" );      // F1
          }
     else if ( X.size() > 3 && X.StartsWith("fos") ) // X=fosY for fosmid Y
     {    int fosno = X.After("fos").Int();
          fosmids = {fosno};
          regions.push_back( FosmidRegion(fosno, PAD) );
          cout << "assembling fosmid id = " << fosno << endl;
          }
     else if ( X == "small" ) // 3.45 Mb
     {    for ( int n = 1; n <= 39; n++ )
          {    if ( n == 16 || n == 33 || n == 35 || n == 36 || n == 38 ) continue;
               fosmids.push_back(n);    }
          res[X].rev = 49066;
          res[X].nedges = 13963, res[X].meanlen = 566.536;
          res[X].gaps = 26, res[X].indels = 43, res[X].subs = 109;
          // Notes for particular Fosmids, run with CUTTING_EDGE = True.  
          // * = good candidate for investigation
          regions.push_back( "1:24.8M-24.9M" );      // F1
          // F2. There is a small gap around 43100.  There is a gap in the same 
          // place in the unassisted LongProto assembly.
          regions.push_back( "1:54.75M-54.85M" );    // F2
          regions.push_back( "1:164.9M-165M" );      // F3
          regions.push_back( "2:106.75M-106.85M" );  // F4
          // F5.  Two gaps.  The first gap is not a real gap but the alternative
          // has a 102-base indel.  The second gap is 294 bases of low-complexity
          // sequence.  It is closed in the unassisted LongProto assembly.
          regions.push_back( "2:239.3M-239.4M" );    // F5
          regions.push_back( "3:11M-11.1M" );        // F6
          regions.push_back( "3:61.5M-61.6M" );      // F7
          regions.push_back( "5:111M-111.1M" );      // F8
          regions.push_back( "5:177.65M-177.75M" );  // F9
          // F9.  Two gaps.  One gap can be eliminated with CONSERVATIVE_KEEP=True
          // and K2_FLOOR=60.
          regions.push_back( "5:179.3M-179.4M" );    // F10
          // F11.  One gap.  The same gap is in the unassisted LongProto assembly.
          regions.push_back( "5:179.2M-179.3M" );    // F11
          // F12.  One gap, new.
          regions.push_back( "6:19.6M-19.7M" );      // F12
          regions.push_back( "6:92.55M-92.65M" );    // F13
          regions.push_back( "7:3.85M-3.95M" );      // F14
          regions.push_back( "7:38.7M-38.8M" );      // F15
          // F17.  Two gaps.
          regions.push_back( "8:23.2M-23.3M" );      // F17
          regions.push_back( "8:30.75M-30.85M" );    // F18
          // F19.  One gap, new.
          regions.push_back( "8:72.75M-72.85M" );    // F19
          // *F20.  One gap.  New to r48420.  This involves a huge homopolymer-like
          // cell that contains the true path but which DeleteLowCoverage wrecks.
          // Turning off DeleteLowCoverage also results in a gap because because
          // the cell contains a cycle.
          regions.push_back( "8:128.75M-128.85M" );  // F20
          regions.push_back( "10:30.85M-30.95M" );   // F21
          // *F22.  One gap.  New to r48420.  This gap can be closed by *deleting*
          // two pids: 6295, 6708.  This suggests that we may need a better local
          // assembly process.
          regions.push_back( "11:44.9M-45M" );       // F22
          regions.push_back( "11:45.5M-45.6M" );     // F23
          // F24.  Two large gaps.  The first one is around 10300 and the second one
          // is around 33300.  If you assemble the region using unassisted 
          // LongProto, the first gap is filled but associated with a 655 base 
          // deletion, and the second gap is a gap in the LongProto assembly.
          regions.push_back( "11:64.95M-65.05M" );   // F24
          regions.push_back( "11:67.7M-67.8M" );     // F25
          regions.push_back( "11:75.45M-75.55M" );   // F26
          regions.push_back( "11:111.8M-111.9M" );   // F27
          regions.push_back( "12:3.1M-3.2M" );       // F28
          // F29.  Two gaps.  One is a large insertion so we don't expect to get it.
          // The other gap is unclosed in unassisted LongProto, but looks like it
          // might be closable.
          regions.push_back( "12:7M-7.1M" );         // F29
          regions.push_back( "12:14.8M-14.9M" );     // F30
          // *F31.  One gap.  Unassisted LongProto assembly does not have a gap.
          // Single Fosmid assembly has the same gap.
          regions.push_back( "12:57.55M-57.7M" );    // F31
          regions.push_back( "12:113.95M-114.05M" ); // F32
          regions.push_back( "14:104M-104.1M" );     // F34
          // F37.  Three gaps.  Not examined yet.
          regions.push_back( "15:30.4M-30.5M" );     // F37
          regions.push_back( "15:74.9M-75.0M" );     // F39
               }
     else if ( X == "medium" ) // 6 Mb
     {    fosmids = {1,2,3,4,5,6,7,8,9,10,11,12};
          res[X].rev = 48385;
          res[X].nedges = 31518, res[X].meanlen = 450.469;
          res[X].gaps = 9, res[X].indels = 14, res[X].subs = 24;
          regions.push_back( "1:24.6M-25.1M" );      // F1
          regions.push_back( "1:54.55M-55.05M" );    // F2
          regions.push_back( "1:164.7M-165.2M" );    // F3
          regions.push_back( "2:106.55M-107.05M" );  // F4
          regions.push_back( "2:239.1M-239.6M" );    // F5
          regions.push_back( "3:10.8M-11.3M" );      // F6
          regions.push_back( "3:61.3M-61.8M" );      // F7
          regions.push_back( "5:110.8M-111.3M" );    // F8
          regions.push_back( "5:177.45M-177.95M" );  // F9
          regions.push_back( "5:178.8M-179.8M" );    // F10,11
          regions.push_back( "6:19.4M-19.9M" );      // F12
               }
     else if ( X == "mediump" ) // 11.5 Mb
     {    fosmids = {1,2,3,4,5,6,7,8,9,10,11,12};
          res[X].rev = 48385;
          res[X].nedges = 62015, res[X].meanlen = 439.684;
          res[X].gaps = 11, res[X].indels = 17, res[X].subs = 36;
          regions.push_back( "1:24.4M-25.4M" );      // F1
          regions.push_back( "1:54.35M-55.35M" );    // F2
          regions.push_back( "1:164.5M-165.5M" );    // F3
          regions.push_back( "2:106.35M-107.35M" );  // F4
          regions.push_back( "2:238.9M-239.9M" );    // F5
          regions.push_back( "3:10.6M-11.6M" );      // F6
          regions.push_back( "3:61.1M-62.1M" );      // F7
          regions.push_back( "5:110.6M-111.6M" );    // F8
          regions.push_back( "5:177.25M-178.25M" );  // F9
          regions.push_back( "5:178.6M-180.1M" );    // F10,11
          regions.push_back( "6:19.2M-20.2M" );      // F12
               }
     else if ( X == "chr11" )                        // chromosome 11
     {    fosmids = {22,23,24,25,26,27,88,89,90,91};
          regions.push_back("11");    }
     else if ( X == "large" ) // the entire aligned genome
     {    fosmids = {1,2,3,4,5,6,7,8,9,10,11,12};
          for ( int c = 1; c <= 22; c++ )
               regions.push_back( ToString(c) );
          regions.push_back( "X" );    }
     else if ( X == "allf" ) 
     {    fosmids = AllFosmids( );
	  regions = AllFosmidRegions(PAD);    } // all fosmids with padding
     else if ( X == "all" ) // the entire genome
     {    if ( SAMPLE == "NA12878" ) 
          {    fosmids = AllFosmids( );
               if ( EVALUATE == "" ) EVALUATE = "True";    }
          all = True;    }
     else if ( X.Contains( ":" ) && X.After( ":" ).Contains( "-" ) )
     {    regions.push_back(X);    }
     else
     {    cout << "Unknown X argument." << endl;
          Scram(1);    }
     if ( F != "default" ) 
     {    vec<int> ftmp;
          ParseIntSet( F, ftmp );
          if ( !Subset(ftmp, fosmids) ) 
          {    cout << "WARNING: your eval fosmids may not intersect with the "
                    << "loaded data" << endl;    }
          fosmids = ftmp;    }    }
