// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

// MakeDepend: library JEMALLOC
// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS
//
#include <omp.h>
#include "10X/runstages/RunStages.h"
#include "system/jemalloc-hooks.h"
#include "10X/astats/GenomeAlign.h"
#include "10X/astats/RefLookup.h"
#include "10X/DfTools.h"
#include "10X/PathsIndex.h"
#include "system/System.h"

void FirstLoadData( String const& work_dir, String const& read_head, vec<String>& cleanupFiles,
          Bool const PREDUP, String const& KEEP, String const& R, vec<String> const& lr,
          vec<double> const& LR_SELECT_FRAC,
          vec<int16_t>& lens, vec<vec<vec<int64_t>>>& qhist,
          vecbvec& bases, ObjectManager<VecPQVec>& quals_om, vec<int64_t>& bci,
          vec<String>& subsam_names, vec<int64_t>& subsam_starts, vec<DataSet>& datasets, 
          int & max_read_length )
{
     auto& quals = quals_om.create();
     LoadData( work_dir, R, lr, LR_SELECT_FRAC, bases, quals_om, bci, subsam_names, subsam_starts, datasets );

     cout << "datatypes loaded:" << endl;
     for ( auto const& d : datasets ) cout << "\t" << d << endl;

     cleanupFiles.push_back( work_dir + "/data/frag_reads_orig.fastb" );
     cleanupFiles.push_back( work_dir + "/data/frag_reads_orig.qualp" );
     cleanupFiles.push_back( work_dir + "/data/frag_reads_orig.bci" );

     if ( PREDUP ) {
          auto nbases = bases.size();
          Predup( bases, quals, bci );
          cout << Date() << ": Predup removed " <<
               ( nbases - bases.size() )  * 10000 / nbases / 100.0 << "% of reads" << endl;
          cout << Date() << ": writing Predup versions of bases, quals" << endl;
          bases.WriteAll( work_dir + "/data/frag_reads_predup.fastb" );
          quals_om.newFile( work_dir + "/data/frag_reads_predup.qualp" );
          quals_om.store();
          BinaryWriter::writeFile( work_dir + "/data/frag_reads_predup.bci", bci );

          cleanupFiles.push_back( work_dir + "/data/frag_reads_predup.fastb" );
          cleanupFiles.push_back( work_dir + "/data/frag_reads_predup.qualp" );
          cleanupFiles.push_back( work_dir + "/data/frag_reads_predup.bci" );
     }

     max_read_length=0;
     lens.resize( bases.size( ) );
     for ( int64_t i = 0; i < (int64_t) bases.size( ); i++ ) {
          lens[i] = bases[i].size( );
          if ( max_read_length < lens[i] )
               max_read_length=lens[i];
     }

     GetQualStats( quals, qhist, max_read_length );

     if ( KEEP != "none" ) {
          // this writes to either pre-dup or orig as needed
          cout << Date() << ": writing " << work_dir << read_head << ".{lens,qhist}" << endl;
          BinaryWriter::writeFile(
               work_dir + read_head + ".lens", lens );
          BinaryWriter::writeFile(
               work_dir + read_head + ".qhist", qhist );
          BinaryWriter::writeFile(
               work_dir + read_head + ".dti", datasets );
     }

     quals_om.unload();
}

int main( int argc, char *argv[] )
{    RunTime( );
     double all_clock = WallClockTime( );
     String start_time = Date( );

     BeginCommandArguments;
     CommandArgument_Bool_OrDefault_Doc(TRACK_SOME_MEMORY, False,
          "track some memory allocations, for debugging");

     // product development command line arguments
     // cleaner would be "PD" versions of CommandArguments.  
     // Please follow the lead and keep these grouped.
#if !defined(CS)
     CommandArgument_String_OrDefault_Doc(SAMPLE, "NA12878", 
          "sample name; human or NA12878 or HGP or unknown");
     CommandArgument_String_OrDefault_Doc(R, "",
          "comma-separated list of plain-read input files" );
     CommandArgument_String_OrDefault_Valid_Doc(PIPELINE, "pd", "{cs,pd}",
               "are we running a CS or PD pipeline" );
     CommandArgument_Bool_OrDefault_Doc(JE_STATS, False,
               "print JEMALLOC stats at exit");
     CommandArgument_String_OrDefault_Doc(DATASET, "1",
          "1, X1, X2, X3 or X4, and only meaningful for SAMPLE=NA12878");
     CommandArgument_String_OrDefault_Doc(INSTANCE, "1",
          "to allow multiple concurrent runs");
     CommandArgument_Bool_OrDefault_Doc(ALIGN, True,
          "align to entire genome");
     CommandArgument_Bool_OrDefault_Doc(FINAL_ALIGN, False, "generate alignsb just for the a.base stage");
     CommandArgument_Bool_OrDefault_Doc(CG2, True, "use CloseGap2");
     CommandArgument_Bool_OrDefault_Doc(STACKSTER, True, "use Stackster");
     CommandArgument_Bool_OrDefault_Doc(STACKSTER_ALT, False, "passed to Stackster");
     CommandArgument_Double_OrDefault_Doc(GRAPHMEM, 0.9,
		 "fraction of memory to allow ReadQGrapher to grab" );
     CommandArgument_Bool_OrDefault_Doc(EXIT_LOAD, False,
          "exit after loading and writing data");
     CommandArgument_String_OrDefault_Doc(CHR, "", "if set, only align to this "
          "chromosome, which should be an integer");
     CommandArgument_String_OrDefault_Doc(FINAL, "a.base",
          "name of final directory");
     CommandArgument_String_OrDefault_Doc(GRAPH, "", "if set, use MSP graph csv file ");
     CommandArgument_Bool_OrDefault_Doc(PREDUP, False,
               "use Predup for cleaning input reads");
     CommandArgument_String_OrDefault_Doc(REF, "hg19", 
          "reference sequence, either hg19 or fos100");
     CommandArgument_Bool_OrDefault_Doc(RESCUE, False, "rescue kmers");
     // NON-ALGORITHMIC OPTIONS

     CommandArgument_String_OrDefault_Doc(USER, "",
          "use this instead of Getenv(USER)");
     CommandArgument_String_OrDefault_Doc(ROOT, "/mnt/assembly",
          "root for output directory location");
     CommandArgument_String_OrDefault_Doc(REPORT, "",
               "output summary report to a file" );
     // ALGORITHMIC HEURISTICS

     CommandArgument_Int_OrDefault_Doc(K, 48, "K value for base graph");
     CommandArgument_Int_OrDefault_Doc(MIN_FREQ, 3, "passed to ReadQGrapher");
     CommandArgument_Int_OrDefault_Doc(MIN_BC, 2, "passed to ReadQGrapher");
     CommandArgument_Int_OrDefault_Doc(MIN_QUAL, 7, "passed to ReadQGrapher");
     CommandArgument_Int_OrDefault_Doc(MIN_LINKS, 4, "minimum links to patch a gap");
     CommandArgument_Bool_OrDefault_Doc(ONE_GOOD, False,
          "find links from dead end to anywhere");
     CommandArgument_Bool_OrDefault_Doc(BACK_EXTEND, False,
          "try to extend paths backward");
     CommandArgument_String_OrDefault_Doc(START, "",
          "start at a particular point in the assembly process:\n"
          "(shown below in order)\n"
          "loaded: right after data is loaded\n"
          "patch: at gap patching\n"
          "trim: at trimming\n"
          "dups: at duplicate marking\n"
          "alltinks: at barcode link finding\n"
               );
#else
     // Use these defaults for the CS pipeline
     String SAMPLE       = "unknown";
     String R            = "";
     String PIPELINE     = "cs";
     Bool JE_STATS       = False;
     String DATASET      = "1";
     String INSTANCE     = "1";
     Bool ALIGN          = False;
     Bool FINAL_ALIGN    = False;
     Bool CG2            = True;
     Bool STACKSTER      = True;
     Bool STACKSTER_ALT  = False;
     double GRAPHMEM     = 0.9;
     Bool EXIT_LOAD      = False;
     String CHR          = "";
     String FINAL        = "a.base";
     String GRAPH        = "";
     Bool PREDUP         = False;
     String REF          = "hg19";
     Bool RESCUE         = False;
     String USER         = "";
     String ROOT         = "/mnt/assembly";
     String REPORT       = "";
     int K               = 48;
     int MIN_FREQ        = 3;
     int MIN_BC          = 2;
     int MIN_QUAL        = 7;
     int MIN_LINKS       = 4;
     Bool ONE_GOOD       = False;
     Bool BACK_EXTEND    = False;
     String START        = "";
#endif

     // CS pipeline uses these arguments
     CommandArgument_String_OrDefault_Doc(LR, "",
          "comma-separated list of linked-read input files, each ending in fastb, "
          "each obtained from ParseBarcodedBam");
     CommandArgument_DoubleSet_OrDefault_Doc(LR_SELECT_FRAC, "1.0", "subsample LR data" );
     CommandArgument_String_OrDefault_Doc(OUT_DIR, "", "name of output directory"); 
     
     // PERFORMANCE OPTIONS - THREADS AND MEMORY
     CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS,0,
          "Number of threads.  By default, the number of processors online.");
     CommandArgument_Double_OrDefault_Doc(MAX_MEM_GB, 0,
          "if specified, maximum *suggested* RAM use in GB; in some cases may be "
          "exceeded by our code");
     
     EndCommandArguments;

     if (TRACK_SOME_MEMORY) DeclareThatWeAreTrackingSomeMemory( );
     
     // retired or fixed or massaged command-line arguments
     String KEEP = "all";

     if ( PIPELINE == "cs" ) SAMPLE="unknown";

     if ( SAMPLE == "unknown" ) ALIGN=False;

     // Check args.
     RunStages start( START, {"", "loaded", "patch", "trim", "dups", "alltinks"} );

     ForceAssert( K == 40 || K == 48 || K == 60 );
     if ( KEEP == "none" && START == "loaded" )
     {    cout << "I'm not sure you really want to do this, since it will\n"
               << "delete your starting files.  So I'm going to quit." << endl;
          Scram(1);    }
     if ( R == "" && LR == "" && START == "" )
     {    cout << "I'm not sure you really want to do this, since it may\n"
               << "delete your starting files.  So I'm going to quit." << endl;
          Scram(1);    }

     // Set up directories, etc.

     String work_dir = ROOT + "/GapToy/" + INSTANCE;
     if ( OUT_DIR != "" ) work_dir = OUT_DIR;
     Mkpath(work_dir);
     String fin_dir = work_dir + "/" + FINAL;
     if ( KEEP != "none" ) Mkdir777(fin_dir);

     // Initialize martian alert logger
     // after OUT_DIR has been created
     Martian::init(work_dir);
     
     // Initialize logger
     StatLogger::init("", work_dir + "/alerts.list");


     MaybeWriter finalWriter( KEEP != "none" );        // conditions determine whether these actually emit files or not
     MaybeWriter nonFinalWriter( KEEP == "all" );      // usage: finalWriter.writeFiles( filename, obj ) same as
                                                       //        BinaryWriter::writeFiles
     vec<String> cleanupFiles;  // files to be written, but deleted if KEEP=none

     // Set computational limits, etc.

     SetThreads( NUM_THREADS, False );
     SetMaxMemoryGBCheck(MAX_MEM_GB);
     
     const int PI_CHUNKS = 30;
     // Parse linked read data request.

     vec<String> lr;
     ParseStringSet( "{" + LR + "}", lr );

     for ( size_t i = 0; i < lr.size(); ++i )
     {    String head = lr[i].Before(".fastb");
          if ( !IsRegularFile( head + ".fastb" )
               || !IsRegularFile( head + ".qualp" )
               || !IsRegularFile( head + ".bci" ) )
          {    cout << "\nCan't file your LR input files " << head << ".*." << endl;
               cout << "Giving up." << endl;
               Scram(1);    }    }

     // Handle samples and subsamples.  A given sample may be
     // divided into multiple subsamples.

     vec<String> subsam_names = { "C" };
     BinaryWriter::writeFile( work_dir + "/subsam.names", subsam_names );
     vec<int64_t> subsam_starts( subsam_names.size( ), 0 );

     // Define heuristics.

     Remove( work_dir + "/clock.log" );
     LogTime( 0, "", work_dir );
     {    string hostname = getHostName();
          hostname = hostname.substr( 0, hostname.find('.') );
          OfstreamMode( xout, work_dir + "/the_command", ios::app );
          xout << "\n" << hostname << ": " << command.TheCommand( ) << endl;    
          // write commit hash
          Ofstream( yout, work_dir + "/.commit_hash" );
          yout << command.GetCommitHash( );
          }
     String run_head = work_dir + "/a";
     String tmp_dir1 = work_dir + "/data";
     Mkdir777(tmp_dir1);
     String log_dir = work_dir + "/logs";
     Mkdir777(log_dir);
     SystemSucceed( "/bin/rm -rf " + work_dir + "/loc/*" );
     Mkdir777( work_dir + "/loc" ), Mkdir777( work_dir + "/special" );

     // Load data

     vecbvec bases;
     ObjectManager<VecPQVec> quals_om(work_dir + "/data/frag_reads_orig.qualp");
     vec<int64_t> bci;
     vec<vec<vec<int64_t>>> qhist;
     vec<int16_t> lens;
     int max_read_length=0;
     const char* read_head = ( PREDUP ? "/data/frag_reads_predup" : "/data/frag_reads_orig" );
     vec<DataSet> datasets;


     if ( START == "" ) {
          FirstLoadData( work_dir, read_head, cleanupFiles, PREDUP, KEEP, R, lr,
                    LR_SELECT_FRAC, lens, qhist, bases, quals_om, bci, 
                    subsam_names, subsam_starts, datasets, max_read_length );
     } else {
          /* this is awful */
          quals_om.newFile( work_dir + read_head + ".qualp" );
     }


     if ( START == "loaded" )
     {    if ( !IsRegularFile( work_dir + read_head + ".lens" ) )
          {    vecbasevector bases( work_dir + read_head + ".fastb" );
               lens.resize( bases.size( ) );
               for ( int64_t i = 0; i < (int64_t) bases.size( ); i++ )
                    lens[i] = bases[i].size( );
               BinaryWriter::writeFile(
                    work_dir + read_head + ".lens", lens );    }
          if ( !IsRegularFile( work_dir + read_head + ".qhist" ) )
          {    VecPQVec quals( work_dir + read_head + ".qualp" );
               BinaryReader::readFile(work_dir + read_head + ".lens", &lens);
               max_read_length=0;
               for (int64_t i = 0; i != int64_t(lens.size()); i++) {
                    if (max_read_length < lens[i]) 
                         max_read_length=lens[i];
               }
               GetQualStats( quals, qhist, max_read_length );
               BinaryWriter::writeFile(
                    work_dir + read_head + ".qhist", qhist ); }    }
     if ( START != "" )
     {    cout << Date( ) << ": loading data" << endl;
          if ( START != "alltinks" )
          {    bases.ReadAll( work_dir + read_head + ".fastb" ); }
          BinaryReader::readFile( work_dir + read_head + ".bci", &bci );
          cout << Date( ) << ": total reads = "
               << ToStringAddCommas( bases.size( ) ) << endl;
          BinaryReader::readFile( work_dir + read_head + ".lens", &lens );
          BinaryReader::readFile(
                    work_dir + read_head + ".qhist", &qhist );


          BinaryReader::readFile( work_dir + read_head + ".dti", &datasets);
     }

     ForceAssertGt(datasets.size(), 0u);
     int64_t bc_start=0;
     for ( auto const& d : datasets )
          if ( d.dt == ReadDataType::BAR_10X || d.dt == ReadDataType::UNBAR_10X ) {     // R data must come first
               bc_start = d.start;
               break;
          }
     // Block to evaluate Q-score related metrics and warn user
     // bad cycle detection
     // % Q30 on R2
     {
          // Compute whether we have a cycle failure
          // defined as 50 % or greater of the reads 
          // having a base with Q <= 2 at a fixed cycle
          vec <int> bad_cycles(2, 0);
          double max_low_q_base_fraction=0.0;
          int MIN_READS = 1000; // only alert if we have at least 1000 reads 
          for ( int ri = 0; ri != 2; ri++ ) {
               for ( int pos = 0; pos != int(qhist[ri].size()); pos++ ) {
                    float low_q=0.0, all_q=0.0;
                    for ( int q = 0; q != int(qhist[ri][pos].size()); q++ ) {
                         if ( q <= 2 )
                              low_q += qhist[ri][pos][q];
                         all_q += qhist[ri][pos][q];
                    }
                    if ( all_q > MIN_READS ) {
                         if ( max_low_q_base_fraction < low_q/all_q ) 
                              max_low_q_base_fraction = low_q/all_q;
                         if ( low_q/all_q > 0.5 ) {
                              bad_cycles[ri]++;
                         }
                    }
               }
          }
          
          // construct format string
          ostringstream msg;
          if ( bad_cycles[0] > 0 )
               msg << bad_cycles[0] << " cycle(s) on read one";
          if ( bad_cycles[0] > 0 && bad_cycles[1] > 0 )
               msg << ", ";
          if ( bad_cycles[1] > 0 )
               msg << bad_cycles[1] << " cycle(s) on read two";
          String format_string = msg.str();
          StatLogger::log("max_low_q_base_frac", max_low_q_base_fraction, \
               "Max low Q base frac per cycle");
          StatLogger::log("bad_cycles_r1", bad_cycles[0], \
               "# bad cycles on R1");
          StatLogger::log("bad_cycles_r2", bad_cycles[1], \
               "# bad cycles on R2");
          StatLogger::issue_alert("max_low_q_base_frac", max_low_q_base_fraction, format_string);

          // % Q30 warning
          double total = 0, total30 = 0;
          int max_read_length = qhist[0].size();
          int max_qual        = qhist[0][0].size();
          for ( int pos = 0; pos < max_read_length; pos++ ) {
               for ( int qv = 0; qv != max_qual; qv++ ) {    
                    total += qhist[1][pos][qv];
                    if ( qv >= 30 ) total30 += qhist[1][pos][qv];    }
          }
          double q30_r2_perc = 0;
          if ( total > 0 )
               q30_r2_perc = 100.0 * double(total30) / double(total);
          StatLogger::issue_alert("q30_r2_perc", q30_r2_perc);
          StatLogger::log("q30_r2_perc", q30_r2_perc, \
               "Percent of bases on read 2 with Q-score >= 30", true );
     }
     // test number of reads that have white-listed barcodes
     // Issue a Martian::exit or Martian::alarm accordingly
     {
          int64_t unbar_start = 0;
          int64_t bar_start   = 0;
          int64_t total_reads = bases.size();
          for ( auto const & d: datasets ) {
               if ( d.dt == ReadDataType::BAR_10X )
                    bar_start = d.start;
               if ( d.dt == ReadDataType::UNBAR_10X )
                    unbar_start = d.start;
          }
          // this is the fraction of 10X reads that have whitelisted barcodes 
          double valid_bc_perc = 0;
          if (total_reads > unbar_start )
               valid_bc_perc = 100*double( total_reads - bar_start )/ double(total_reads - unbar_start);
          cout << Date( ) << ": " << valid_bc_perc 
               << "%% of linked reads have valid barcodes" << endl;
          StatLogger::log( "valid_bc_perc", valid_bc_perc, \
               "% reads with valid barcodes", true);
          StatLogger::issue_alert( "valid_bc_perc", valid_bc_perc );
     }
     cout << Date() << ": " << ToStringAddCommas(bases.size()) << " total reads loaded (R+LR)." << endl;
     MEM(reads_loaded);
     PRINT2(subsam_names.size(), subsam_starts.size() );

     //ForceAssert( quals_om.filename().StartsWith( work_dir + read_head ) );
     String quals_head = quals_om.filename().SafeAfterLast("/").Before(".");
     ForceAssert( String(read_head).EndsWith( quals_head ) );

     // Sanity check barcode counts.
     SanityCheckBarcodeCounts(bci);

     // Expand barcode index;

     cout << Date( ) << ": expanding barcode index" << endl;
     vec<int32_t> bc( bci.back( ), -1 );
     #pragma omp parallel for
     for ( int b = 0; b < bci.isize( ) - 1; b++ )
     {    int64_t start = bci[b], stop = bci[b+1];
          for ( int64_t j = start; j < stop; j++ )
               bc[j] = b;    }

     // Report some data stats.

     if ( START != "alltinks" )
     {    int min_read = 1000000000, max_read = 0;
          int64_t total_bases = 0, total_kmers = 0;
          for ( int64_t i = 0; i < (int64_t) bases.size( ); i++ )
          {    min_read = Min( min_read, bases[i].isize( ) );
               max_read = Max( max_read, bases[i].isize( ) );
               total_bases += bases[i].size( );
               if ( bases[i].isize( ) >= K )
                    total_kmers += bases[i].isize( ) - K + 1;    }
          cout << "read sizes: ";
          if ( min_read == max_read ) cout << min_read;
          else cout << min_read << "-" << max_read;
          cout << endl;

          int64_t gsize = 3200 * (int64_t) 1000000;
          double cov = double(total_bases) / gsize;
          double cov_K = double(total_kmers) / gsize;
          PRINT2( cov, cov_K );    }

     // Force subsam stuff.

     subsam_names.resize(1);
     subsam_names[0] = "C";
     subsam_starts.resize(1);
     subsam_starts[0] = 0;
     BinaryWriter::writeFile( work_dir + "/subsam.starts", subsam_starts );
     BinaryWriter::writeFile( work_dir + "/subsam.names", subsam_names );
     if (EXIT_LOAD) Scram(0);

     // evil copying of reference
     // TODO: make this a parameter
     if ( START == "" && ALIGN )
     {
          ForceAssert( REF == "hg19" || REF == "fos100" );
          if ( REF == "hg19" )
          {
               Cp2( "/mnt/opt/meowmix_git/assembly/refs/hg19/genome.fastb", work_dir );
               Cp2( "/mnt/opt/meowmix_git/assembly/refs/hg19/genome.ambint", work_dir );
               Cp2( "/mnt/opt/meowmix_git/assembly/refs/hg19/genome.names", work_dir );
          }
          if ( REF == "fos100" )
          {
               Cp2( "/mnt/opt/meowmix_git/assembly/refs/NA12878/fos100.fastb", work_dir + "/genome.fastb" );
               Cp2( "/mnt/opt/meowmix_git/assembly/refs/NA12878/fos100.ambint", work_dir + "/genome.ambint" );
               Cp2( "/mnt/opt/meowmix_git/assembly/refs/NA12878/fos100.names", work_dir + "/genome.names" );
          }
     }

     cleanupFiles.push_back( work_dir + "/genome.fastb" );
     cleanupFiles.push_back( work_dir + "/genome.ambint" );
     cleanupFiles.push_back( work_dir + "/genome.names" );

     // Set up assembly data structures.

     HyperBasevector hbv;
     HyperBasevectorX hb;
     vec<int> inv;
     ReadPathVec paths;
     VecULongVec paths_index;
     vecbasevector genome;
     vec< pair<int,ho_interval> > ambint;
     
     if ( ALIGN ) {
          genome.ReadAll( work_dir + "/genome.fastb" );
          if ( CHR != "" && CHR.IsInt( ) )
          {    int g = CHR.Int( ) - 1;
               for ( int j = 0; j < (int) genome.size( ); j++ )
                    if ( j != g ) genome[j].resize(0);    }
          BinaryReader::readFile( work_dir + "/genome.ambint", &ambint );
     }
     vecbasevector G;

     if ( IsRegularFile(work_dir+"/sample") ) 
          Remove(work_dir+"/sample");
     if (SAMPLE != "unknown")
     {    FetchFinished( SAMPLE, G );
          Ofstream( out, work_dir + "/sample" );
          out << SAMPLE << endl;    
     } 
     MasterVec< SerfVec<triple<int,int,int> > > alignsb;

     // Now make the hbv.
     vec<Bool> dup;
     if ( start.at_or_before("loaded") ) {
         // TODO: WAIT
          StageBuildGraph( K, bases, quals_om, MIN_QUAL, MIN_FREQ, MIN_BC,
                    bc, bc_start, GRAPH, GRAPHMEM, work_dir, read_head,
                    hbv, paths, inv);
          
          // write files here
          String dir = work_dir + "/a." + ToString(K);
          cout << Date( ) << ": inverting paths index, mem usage = "
               << MemUsageGBString( ) << endl;
          Mkdir777( dir );
          writePathsIndex( paths, inv, dir, "a.paths.inv", "a.countsb", 
                           PI_CHUNKS, false ); // no verbose
          
          if ( KEEP == "all" ) {
               WriteAssemblyFiles( hbv, inv, paths, paths_index, bci, ALIGN,
                    genome, dir, alignsb );
          }

          VirtualMasterVec<PQVec> vmv( quals_om.filename() );
          double interdup_rate;
          MarkDups( bases, vmv, paths, bc, dup, interdup_rate );
          nonFinalWriter.writeFile( work_dir + "/a." + ToString(K) + "/a.dup", dup);
     }

     // Start to patch gaps.
     vec<basevector> closures;
     vec<Bool> bad;
     ReadPathVecX pathsX;
     
     if ( start.at("patch") ) {
          cout << Date()  << ": starting at patch... reading data" << endl;
          quals_om.unload();
          String dir = work_dir + "/a." + ToString(K);
          BinaryReader::readFile( dir + "/a.hbv", &hbv );
          BinaryReader::readFile( dir + "/a.hbx", &hb );
          BinaryReader::readFile( dir + "/a.inv", &inv );
          paths.ReadAll( dir + "/a.paths" );
          if ( !IsRegularFile( dir + "/a.dup" ) ) {
               cout << Date()  << ": dups file missing, recreating it" << endl;
               VirtualMasterVec<PQVec> vmv( quals_om.filename() );
               double interdup_rate;
               MarkDups( bases, vmv, paths, bc, dup, interdup_rate );
               nonFinalWriter.writeFile( dir + "/a.dup", dup);
          } else
               BinaryReader::readFile( dir + "/a.dup", &dup );
          //paths_index.ReadAll( dir + "/a.paths.inv" );
     }
    
     if ( start.at_or_before("patch") ) {
          String dir = work_dir + "/a." + ToString(K);
          const int max_width = 400;
          vec< pair<int,int> > pairs;
          String paths_index_file = dir + "/a.paths.inv";
          StagePatch(dir, K, bases, quals_om, hbv, hb, paths, paths_index_file,
                    inv, dup,
                    bad, datasets, bc, max_width, ONE_GOOD, closures, pairs, CG2,
                    STACKSTER, STACKSTER_ALT, RESCUE );
          Validate( hbv, paths );
          cout << Date( ) << ": hbv has checksum " << hbv.CheckSum( ) << endl;

          cout << Date() << ": re-read bases :(" << endl;
          bases.ReadAll( work_dir + read_head + ".fastb" );
          cout << Date() << ": re-read bases :( done :)" << endl;
          nonFinalWriter.writeFile( work_dir + "/a." + ToString(K) + "/a.hops", pairs );
          nonFinalWriter.writeFile( work_dir + "/a." + ToString(K) + "/a.bad", bad );
          BinaryWriter::writeFile( work_dir + "/closures.fastb", closures );

          // Extend paths.
          //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CONVERT TO PATHSX HERE !!!!!!!!!!
          pathsX.parallel_append(paths,hb);
          Destroy(paths);
          //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          double clock5 = WallClockTime( );
          cout << Date( ) << ": reloading quals, current mem = "
               << MemUsageGBString( ) << ", peak = "
               << PeakMemUsageGBString( ) << endl;
          // Note that presumably we could instead bring in quals as a
          // VirtualMasterVec.
          StageExtension( hbv, inv, bases, quals_om, pathsX, BACK_EXTEND );
          cout << TimeSince(clock5) << " used in new stuff 5" << endl;
          cout << "now current mem = " << MemUsageGBString( ) << endl;
          // NOTE EXPENSIVE CONVERSION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          hb = HyperBasevectorX(hbv);

          // Write files.
          quals_om.unload();
          MEM(patched_write_quals_unload);
          Destroy(bases);
          MEM(patched_write_bases_unload);

          cout << Date( ) << ": inverting paths index, mem usage = "
               << MemUsageGBString( ) << endl;
          Mkdir777( work_dir + "/a.patched" );
          writePathsIndex( pathsX, hb, inv, work_dir + "/a.patched", 
                           "a.paths.inv", "a.countsb", PI_CHUNKS, false );
          if ( KEEP == "all" ) WriteAssemblyFiles( hbv, inv, pathsX, paths_index,
               bci, ALIGN, genome, work_dir + "/a.patched", alignsb );

          // Report assembly statistic.

          {    vec<int> len( hbv.E( ) );
               for ( int e = 0; e < hbv.E( ); e++ )
                    len[e] = hbv.Kmers(e);
               Sort(len);
               cout << "N50 edge length = " << N50(len) << endl;     }

          bases.ReadAll( work_dir + read_head + ".fastb" );
          MEM(patched_write_reload_bases);
     }

     if ( start.at("trim") ) {
          cout << Date( ) << ": loading assembly" << endl;
          #pragma omp parallel sections
          {
               #pragma omp section
               { BinaryReader::readFile( work_dir + "/a.patched/a.hbv", &hbv ); }
               #pragma omp section
               { 
                 BinaryReader::readFile( work_dir + "/a.patched/a.hbx", &hb ); 
                 pathsX.ReadAll( work_dir +"/a.patched/a.pathsX");
               }
               #pragma omp section
               { BinaryReader::readFile( work_dir + "/a.patched/a.inv", &inv ); }
          }
     }


     // Trim.
     vec<int64_t> icount;

     if ( start.at_or_before("trim") ) {
          String paths_index_file = work_dir + "/a.patched/a.paths.inv";
          String countsb_file     = work_dir + "/a.patched/a.countsb";
          // TODO: PATHX CLEANUP( ) IN STAGETRIM
          StageTrim( hb, hbv, inv, bases, quals_om, pathsX,
                     paths_index_file, countsb_file, bc);
          cout << Date( ) << ": hbv has checksum " << hbv.CheckSum( ) << endl;
          hb = HyperBasevectorX(hbv);
          Validate( hbv, hb, pathsX );
          // update paths_index after trimming
          writePathsIndex( pathsX, hb, inv, fin_dir, "a.paths.inv", "a.countsb", 
                           PI_CHUNKS, false );

          // Write files.
          if ( hb.E( ) == 0 )
               Martian::exit("The assembly graph contains zero edges.");
          if ( KEEP != "none" )
          {    
               WriteAssemblyFiles(hbv, inv, pathsX, paths_index,
                    bci, (ALIGN || FINAL_ALIGN), genome, fin_dir, alignsb );

               if ( SAMPLE != "unknown" && ( ALIGN || FINAL_ALIGN ) )
               {    MasterVec< SerfVec<triple<int,int,int> > > alignsb;
                    const int align_genome_K = 60;
                    GenomeAlign<align_genome_K>( hb, inv, G, alignsb );
                    alignsb.WriteAll( fin_dir + "/a.fin.alignsb" );    }
               else Remove( fin_dir + "/a.fin.alignsb" );

               FragDist( hb, inv, pathsX, icount );
               // Warn if insert is too short
               {
                    int64_t NI = Sum(icount), isum = 0;
                    int m;
                    for ( m = 0; m < icount.isize( ); m++ )
                    {    isum += icount[m];
                         if ( isum >= NI/2 ) break;    }
                    StatLogger::log("median_ins_sz", m, "median insert size", true);
                    StatLogger::issue_alert( "median_ins_sz", m );
               }

               BinaryWriter::writeFile( fin_dir + "/a.ins_dist", icount );    }
          else
               Mkdir777(fin_dir);

     }

     if ( start.at_or_after("dups" ) ) {
          cout << Date( ) << ": loading assembly" << endl;
          #pragma omp parallel sections
          {
          #pragma omp section
          {    BinaryReader::readFile( fin_dir + "/a.hbv", &hbv ); }
          #pragma omp section
          {    
               BinaryReader::readFile( fin_dir + "/a.hbx", &hb ); 
               pathsX.ReadAll( fin_dir +"/a.pathsX");
          }
          #pragma omp section
          {    BinaryReader::readFile( fin_dir + "/a.inv", &inv ); }
          #pragma omp section
          {    BinaryReader::readFile( work_dir + "/closures.fastb", &closures ); }
          #pragma omp section
          {    if ( ALIGN ) alignsb.ReadAll( fin_dir + "/a.alignsb" ); }
          #pragma omp section
          {    BinaryReader::readFile( fin_dir + "/a.ins_dist", &icount ); }
          }
     }

     // Mark duplicates and bads.
     if ( start.at_or_before("dups") ) {
          double interdup_rate;
          StageDups( hb, bases, quals_om, pathsX, bc, dup, bad, interdup_rate);
          MEM(before_destroy_bases);
          Destroy(bases);
          MEM(before_unload_again);
          quals_om.unload();
          MEM(after_unload_again);
          finalWriter.writeFile( fin_dir + "/a.dup", dup );
          Ofstream( out, fin_dir + "/a.interdup" );
          out << interdup_rate << endl;
          finalWriter.writeFile( fin_dir + "/a.bad", bad );


          // Map the closures back to the assembly.  Note that because the assembly
          // has been trimmed, not all the closures will map back.  Note also that
          // clop is not sync'ed to closures, in general clop will have less entries
          // than closures.  And some entries in clop will be zero-length readpaths.

          /*
          if ( K == 40 ) MapClosures<40>( hb, closures, clop );
          if ( K == 48 ) MapClosures<48>( hb, closures, clop );
          if ( K == 60 ) MapClosures<60>( hb, closures, clop );
          */

          // Compute edge to barcode map.
          {
          VecIntVec ebcx;
          StageEBC( pathsX, hb, bc, inv, bci, ebcx );

          if (KEEP != "none") ebcx.WriteAll( fin_dir + "/a.ebcx" );
          }

          // Make closures and the graph from them.

          {
          vec<vec<int>> all_closures;
          digraphE<vec<int>> D;
          vec<int> dinv;
          String paths_index_file = fin_dir + "/a.paths.inv";
          Mkdir777( fin_dir + "/orig" );
          StageClosures( hb, inv, pathsX, paths_index_file, dup, bad, all_closures,
                         D, dinv, fin_dir );
          const int64_t base_checksum = hbv.CheckSum( );
          cout << Date( ) << ": hbv has checksum " << base_checksum << endl;
          StatLogger::log("base_checksum", base_checksum, "Base checksum");
          cout << Date( ) << ": D has checksum " << D.CheckSum( ) << endl;
          finalWriter.writeFile( fin_dir + "/orig/a.sup", D );
          finalWriter.writeFile( fin_dir + "/orig/a.sup.inv" , dinv );
          // finalWriter.writeFile( fin_dir + "/orig/a.all_closures", all_closures );
          }
     }

     // Find barcode links.

     vec< triple<int,int,int> > qept;
     VecIntVec ebcx;

     if ( start.at("alltinks") ) {
          BinaryReader::readFile( fin_dir + "/a.dup", &dup );
          BinaryReader::readFile( fin_dir + "/a.bad", &bad );
     }
     if ( start.at_or_before("alltinks") ) {
          STAGE(AllTinks);
          ebcx.ReadAll( fin_dir + "/a.ebcx" );
          AllTinksCore( hb, inv, pathsX, bci, ebcx, qept );
          finalWriter.writeFile( fin_dir + "/a.bc_links", qept );
     }

     cout << "\n" << Date( ) << ": done, " << TimeSince(all_clock)
          << " used, peak mem = " << PeakMemUsageGBString( ) << endl << endl;

     if (KEEP=="none")
     {
          cout << Date() << ": Removing files because KEEP=none:" << endl;
          for ( auto const& fn : cleanupFiles ) {
               cout << "\t" << fn << endl;
               Remove(fn);
          }
     }
     // Record etime, peak mem
     double etime = (WallClockTime( ) - all_clock)/(60.*60.);
     StatLogger::log( "etime_df_h", etime, "Elapsed time DF" );
     StatLogger::log( "mem_peak_df_gb", PeakMemUsageGB(), "Mem peak DF (gb)" );

     // Done.
     StatLogger::write( fin_dir + "/a.perf_stats" );
#ifdef JEMALLOC_HOOKS
     if ( JE_STATS ) (void) je_malloc_stats_print(nullptr, nullptr, nullptr);
#endif

     Scram(0);
}

