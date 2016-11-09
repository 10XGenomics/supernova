// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

// CP.  Original naming from "Close Pairs", doesn't make a lot of sense now.
// This is a major Supernova module, run after DF in the Supernova pipeline.

// MakeDepend: library JEMALLOC
// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "system/jemalloc-hooks.h"

#include "FetchReads.h"
#include "Intvector.h"
#include "MainTools.h"
#include "ParallelVecUtilities.h"
#include "graph/DigraphTemplate.h"
#include "pairwise_aligners/SmithWatAffine.h"
#include "paths/HyperBasevector.h"
#include "paths/long/DiscovarTools.h"
#include "paths/long/ReadPath.h"
#include "paths/long/large/Lines.h"
#include "system/HostName.h"
#include "10X/BuildLocal.h"
#include "10X/Capture.h"
#include "10X/Decycle.h"
#include "10X/DfTools.h"
#include "10X/Flipper.h"
#include "10X/Gap.h"
#include "10X/Gaprika.h"
#include "10X/Heuristics.h"
#include "10X/IntIndex.h"
#include "10X/InvFix.h"
#include "10X/LineOO.h"
#include "10X/LocalTools.h"
#include "10X/PlaceReads.h"
#include "10X/PullApart.h"
#include "10X/Scaffold.h"
#include "10X/Splat.h"
#include "10X/Stackaroo.h"
#include "10X/Star.h"
#include "10X/Super.h"
#include "10X/SuperFiles.h"
#include "10X/astats/AlignFin.h"
#include "10X/astats/AssemblyStats.h"
#include "10X/astats/RefAlign.h"
#include "10X/astats/RefLookup.h"
#include "10X/astats/View.h"
#include "10X/paths/ReadPathVecX.h"
#include "10X/MakeHist.h"

int main( int argc, char *argv[] )
{    RunTime( );
     double clock = WallClockTime( );

     BeginCommandArguments;
     CommandArgument_Bool_OrDefault_Doc(TRACK_SOME_MEMORY, False,
          "track some memory allocations, for debugging");

     // Define arguments that are not in the customer pipeline.

     #if !defined(CS)
     CommandArgument_String_OrDefault_Doc(SAMPLE, "NA12878",
          "sample name; human or NA12878 or HGP or unknown");
     CommandArgument_Bool_OrDefault_Doc(WRITE, True,
          "write output files; not fully respected");
     CommandArgument_String_OrDefault_Doc(READ_SUB, "",
          "sub for input supergraph file names");
     CommandArgument_String_OrDefault_Doc(WRITE_SUB, "",
          "sub for output supergraph file names");
     CommandArgument_Bool_OrDefault_Doc(BUILD_ONLY, False, "stop after build");
     CommandArgument_Bool_OrDefault_Doc(SAVE_LOCAL1, False,
          "save local assemblies at first stage");
     CommandArgument_Bool_OrDefault_Doc(SAVE_LOCAL2, False,
          "save local assemblies at second stage");
     CommandArgument_Int_OrDefault_Doc(S1, -1, "for testing using just this value");
     CommandArgument_Bool_OrDefault_Doc(REQUIRE_SIMPLER, True,
          "only allow particularly simple patches");
     CommandArgument_Double_OrDefault_Doc(MIN_ADVANTAGE, 60.0, "passed to Star");
     CommandArgument_Bool_OrDefault_Doc(OO, False, "passed to Star");
     CommandArgument_Bool_OrDefault_Doc(SURGERY_VERBOSE2, False,
          "extra logging for round 2 surgery");
     CommandArgument_Int_OrDefault_Doc(MIN_SPLAY1, 3500, "minimum line to splay");
     CommandArgument_Bool_OrDefault_Doc(NEW_TEST, True, "passed to BuildLocal2");
     CommandArgument_Bool_OrDefault_Doc(JE_STATS, False, 
          "print JEMALLOC stats at exit");
     CommandArgument_String_OrDefault_Doc(START, "", "can be path or scaffold or "
          "build or surgery or star or starstar or patch or phase or post or "
          "fix or fase or tidy or canon or report");
     CommandArgument_Bool_OrDefault_Doc(DEBUG, False, 
          "turn on extra output; intended for R&D assemblies");
     CommandArgument_Bool_OrDefault_Doc(REWRITE, False, 
          "special option: read assembly and write more files, "
          "intended for use with Probe");

     // EXPERIMENTAL OPTIONS IN PROGRESS.

     CommandArgument_Bool_OrDefault_Doc(ULTRA, False, 
          "do things to make longer scaffolds");
     CommandArgument_Bool_OrDefault_Doc(GLUE, False, "glue after surgery");
     CommandArgument_Int_OrDefault_Doc(MAX_GAP, 100000, "max gap for GAPRIKA");

     // NOTE TWICE!!!

     // Define defaults for command line parameters that are not in customer 
     // pipeline.

     #else
     String SAMPLE            = "unknown";
     Bool WRITE               = True;
     String READ_SUB          = "";
     String WRITE_SUB         = "";
     Bool BUILD_ONLY          = False;
     Bool SAVE_LOCAL1         = False;
     Bool SAVE_LOCAL2         = False;
     int S1                   = -1;
     Bool REQUIRE_SIMPLER     = True;
     double MIN_ADVANTAGE     = 60.0;
     Bool OO                  = False;
     Bool SURGERY_VERBOSE2    = False;
     int MIN_SPLAY1           = 3500;
     Bool NEW_TEST            = True;
     Bool JE_STATS            = False;
     String START             = "";
     Bool DEBUG               = False;
     Bool ULTRA               = False;
     Bool GLUE                = False;
     int MAX_GAP              = 100000;
     Bool REWRITE             = False;
     #endif

     // These parameters are exposed in the customer pipeline.

     CommandArgument_String_OrDefault_Doc(DIR, ".", "assembly directory");

     CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS,0,
         "Number of threads.  By default, the number of processors online.");
     CommandArgument_Double_OrDefault_Doc(MAX_MEM_GB, 0,
          "if specified, maximum *suggested* RAM use in GB; in some cases may be "
          "exceeded by our code");

     CommandArgument_String_OrDefault_Doc(CS_SAMPLE_ID, "", "customer sample id");
     CommandArgument_String_OrDefault_Doc(CS_SAMPLE_DESC, "", 
          "customer sample desc");
     EndCommandArguments;

     if (TRACK_SOME_MEMORY) DeclareThatWeAreTrackingSomeMemory( );
     
     // Initialize martian alert logger.
     
     String BASE_DIR = DIR + "/../";
     Martian::init(BASE_DIR);
     
     // Initialize StatLogger
     String df_stats="", alerts_file="";
     if ( IsRegularFile(DIR + "/a.perf_stats") )
          df_stats=DIR + "/a.perf_stats";
     if ( IsRegularFile(BASE_DIR + "/alerts.list") )
          alerts_file = BASE_DIR + "/alerts.list";
     StatLogger::init( df_stats, alerts_file );

     // Log the command that was run.

     {    string hostname = getHostName();
          hostname = hostname.substr( 0, hostname.find('.') );
          OfstreamMode( xout, DIR + "/../the_command", ios::app );
          xout << "\n" << hostname << ": " << command.TheCommand( ) << endl;    }

     // Setups and sanity checks.

     SetThreads(NUM_THREADS, False);
     SetMaxMemoryGBCheck(MAX_MEM_GB);

     vec<String> stagelist = { "", "path", "scaffold", "build", "surgery", "star", 
          "starstar", "patch", "phase", "post", "fix", "fase", "tidy", "canon", 
          "report", "perfect" };
     int startpos = Position( stagelist, START );
     if ( startpos < 0 )
     {    cout << "Illegal START option." << endl;
          Scram(1);    }
     if ( IsRegularFile( DIR + "/../sample" ) ) {    
          String sample = FirstLineOfFile( DIR + "/../sample" );
          if ( sample != SAMPLE && SAMPLE != "unknown" )
          {    cout << "Inconsistent sample." << endl;
               Scram(1);    }    
     } else 
          SAMPLE="unknown";
     Bool report = ( START == "report" );
     vec<int64_t> bci;
     BinaryReader::readFile( DIR + "/../data/frag_reads_orig.bci", &bci );
     SanityCheckBarcodeCounts(bci);
     Bool ALIGN = ( SAMPLE != "unknown" );
     if (ALIGN)
     {    
          if ( !IsRegularFile( DIR + "/a.alignsb" ) )
          {    cout << "The file: a.alignsb does not exist.\n" <<endl;
               Scram(1);   }
          if ( LastModified( DIR + "/a.hbx" ) > LastModified( DIR + "/a.alignsb" ) )
          {    cout << "Bad time stamp for a.alignsb, suggests something went "
                    << "wrong when running DF.\n" << endl;
               Scram(1);    }    }
     if ( READ_SUB != "" ) READ_SUB = ":" + READ_SUB;
     if ( WRITE_SUB != "" ) WRITE_SUB = ":" + WRITE_SUB;
     Mkdir777( DIR + "/final" + WRITE_SUB );
     Mkdir777( DIR + "/final" + WRITE_SUB + "/stats" );

     // Define and carry forward stage-specific stats.  Also functions to write the
     // stats, along with a commented-out generic lambda that we could sub in if we 
     // were using gcc >= 4.9.  

     vec<pair<String,String>> stagestats;
     stagestats.push( "patch", "histogram_molecules.json" );
     stagestats.push( "patch", "s.lr" );
     stagestats.push( "fase", "s.hetdist" );
     for ( int i = 0; i < stagestats.isize( ); i++ )
     {    String stage = stagestats[i].first, fn = stagestats[i].second;
          if ( startpos > Position( stagelist, stage ) )
          {    String INDIR = DIR + "/" + "final" + READ_SUB;
               String OUTDIR = DIR + "/" + "final" + WRITE_SUB;
               if ( INDIR != OUTDIR && IsRegularFile( INDIR + "/stats/" + fn ) )
               {    Mkdir777(OUTDIR);
                    Mkdir777( OUTDIR + "/stats" );
                    Cp2( INDIR + "/stats/" + fn, 
                         OUTDIR + "/stats/" + fn );    }    }    }
     auto WriteFileStatsInt = [&]( int& x, const String& fn )
     {    String OUTDIR = DIR + "/" + "final" + WRITE_SUB;
          Mkdir777(OUTDIR);
          BinaryWriter::writeFile( OUTDIR + "/stats/" + fn, x );
          String stats_dir = DIR + "/../stats" + WRITE_SUB;
          Mkdir777(stats_dir);
          Cp2( OUTDIR + "/stats/" + fn, stats_dir + "/" + fn );    };
     auto WriteFileStatsVecPairFloatInt = [&]( vec<pair<float, int>>& x, const String& fn )
     {    String OUTDIR = DIR + "/" + "final" + WRITE_SUB;
          Mkdir777(OUTDIR);
          BinaryWriter::writeFile( OUTDIR + "/stats/" + fn, x );
          String stats_dir = DIR + "/../stats" + WRITE_SUB;
          Mkdir777(stats_dir);
          Cp2( OUTDIR + "/stats/" + fn, stats_dir + "/" + fn );    };


     // Declare some variables.

     vec<DataSet> datasets;
     HyperBasevectorX hb;
     vec<int> inv, kmers;
     vec<Bool> dup, bad;
     vecbasevector G;
     VecIntVec ebcx;
     ReadPathVecX paths;
     vec< triple<int,int,int> > qept;
     digraphE<vec<int>> D;
     vec<int> dinv;
     vec<vec<vec<vec<int>>>> dlines;
     vec<int32_t> bc;
     vec<int64_t> bid;
     MasterVec< SerfVec<triple<int,int,int> > > alignsb;

     // Run "perfect" stage, which is special.  What it does: reinsert loops, then
     // compute N50 perfect stretch, and write report files.

     if ( START == "perfect" )
     {    String INDIR = DIR + "/final" + READ_SUB;
          String OUTDIR = DIR + "/final" + WRITE_SUB;
          ForceAssert( SAMPLE != "unknown" );
          BinaryReader::readFile( DIR + "/a.hbx", &hb );
          BinaryReader::readFile( DIR + "/a.inv", &inv );
          BinaryReader::readFile( INDIR + "/a.sup", &D );
          BinaryReader::readFile( INDIR + "/a.sup.inv", &dinv );
          BinaryReader::readFile( INDIR + "/a.sup.lines", &dlines );
          FetchFinished( SAMPLE, G );    
          int64_t N50_perf = -1;
          double errw = -1;
          if ( G.size( ) > 0 )
          {    MasterVec< SerfVec<triple<int,int,int> > > galignsb;
               galignsb.ReadAll( DIR + "/a.fin.alignsb" );
               vec< vec< vec< triple< vec<int>, align, int > > > > Matches;
               String report;
               AlignFin( hb, inv, D, dinv, dlines, G, galignsb,
                    Matches, errw, N50_perf, report, True );
               cout << Date( ) << ": N50 perfect stretch = "
                    << ToStringAddCommas(N50_perf) << endl;
               double errwf = int( round( ( errw*1000 ) * 100 ) ) / 100.0;
               cout << Date( ) << ": error rate = " << errwf
                    << " per kb" << endl;
               Ofstream( rout, OUTDIR + "/o.report" );
               rout << report;    }
          cout << Date( ) << ": done, time used = " << TimeSince(clock) << endl;
          Scram(0);    }

     // Load data.

     cout << Date( ) << ": loading assembly" << endl;
     #pragma omp parallel sections
     {
          #pragma omp section
          {    if ( !report )
               {    BinaryReader::readFile( DIR + "/../data/frag_reads_orig.dti",
                         &datasets );    }    }
          #pragma omp section
          {    BinaryReader::readFile( DIR + "/a.inv", &inv );    }
          #pragma omp section
          {    BinaryReader::readFile( DIR + "/a.kmers", &kmers );    }
          #pragma omp section
          {    if ( !report) BinaryReader::readFile( DIR + "/a.dup", &dup );    }
          #pragma omp section
          {    if ( !report) BinaryReader::readFile( DIR + "/a.bad", &bad );    }
          #pragma omp section
          {    BinaryReader::readFile( DIR + "/a.hbx", &hb );
               if ( !report ) paths.readBinary(DIR + "/a.pathsX");  }
          #pragma omp section
          if ( SAMPLE != "unknown" ) FetchFinished( SAMPLE, G );
          #pragma omp section
          {    if ( !report ) ebcx.ReadAll( DIR + "/a.ebcx" );    }
          #pragma omp section
          {    if ( !report)
                    BinaryReader::readFile( DIR + "/a.bc_links", &qept );    }    }
     int K = hb.K( );
     vecbasevector genome;
     vec< pair<int,ho_interval> > ambint;
     if ( SAMPLE == "human" || SAMPLE == "NA12878" || SAMPLE == "HGP" )
     {    genome.ReadAll( "/mnt/opt/meowmix_git/assembly/refs/hg19/genome.fastb" );
          BinaryReader::readFile(
               "/mnt/opt/meowmix_git/assembly/refs/hg19/genome.ambint", 
               &ambint );    }
     vec<double> COV;

     // Expand barcode index.

     cout << Date( ) << ": expanding barcode index" << endl;
     if ( START != "report" )
     {    bc.resize( bci.back( ), -1 );
          #pragma omp parallel for
          for ( int b = 0; b < bci.isize( ) - 1; b++ )
          {    int64_t start = bci[b], stop = bci[b+1];
               for ( int64_t j = start; j < stop; j++ )
                    bc[j] = b;    }    }
     
     // Compute read two percent proper metric
     // TODO we should be computing this stat in DF. When we allow
     // for StatLogger to be written (SELF_SERIALIZABLE) we can move this
     // computation there. Right now the variable perc_proper_pairs is being
     // passed to ReportAssemblyStats

     double perc_proper_pairs = -1;
     {
          // temporary nonsense until we move things
          /* VirtualMasterVec<ReadPath> vmv(DIR+"/a.paths"); */
          // IMPORTANT: this function returns a percentage
          // not a fraction!!
          ReadTwoPctProper( hb, inv, paths, perc_proper_pairs );
          StatLogger::log("proper_pairs_perc", perc_proper_pairs, "% proper pairs", true);
          StatLogger::issue_alert("proper_pairs_perc", perc_proper_pairs );
     }
     
     // Create a vector of integers, one for each read, such that "having two"
     // nonzero elements is enough.  Note nested parallel for loops, not really
     // what we want.

     cout << Date( ) << ": creating bid" << endl;
     if ( START != "report" )
     {    bid.resize( paths.size( ) );
          #pragma omp parallel for schedule( dynamic, 1 )
          for ( int b = 0; b < bci.isize( ) - 1; b++ )
          {    int64_t start = bci[b], stop = bci[b+1];
               #pragma omp parallel for
               for ( int64_t id = start; id < stop; id++ )
               {    int di;
                    for ( di = 0; di < datasets.isize( ); di++ )
                         if ( id < datasets[di].start ) break;
                    const ReadDataType& dtype = datasets[di-1].dt;
                    if ( dtype == ReadDataType::BAR_10X )
                         bid[id] = (int64_t) paths.size( ) + b + 1;
                    else if ( dtype == ReadDataType::UNBAR_10X ) bid[id] = 0;
                    else if ( dtype == ReadDataType::PCR_FREE ) bid[id] = id + 1;
     
                    // Probably not what we want:
     
                    else if ( dtype == ReadDataType::PCR ) 
                         bid[id] = id + 1;    }    }    }

     // Define writing functions.

     auto BasicWrite = [&]( const String& stage )
     {    if (WRITE)
          {    cout << Date( ) << ": writing files before stage " << stage << endl;
               String OUTDIR = DIR + "/" + stage + WRITE_SUB;
               Mkdir777(OUTDIR);
               BinaryWriter::writeFile( OUTDIR + "/a.sup", D );
               BinaryWriter::writeFile( OUTDIR + "/a.sup.inv", dinv );
               Remove( OUTDIR + "/a.sup.lines" );
               Remove( OUTDIR + "/a.dpaths" );
               Remove( OUTDIR + "/a.dpaths.index" );
               Remove( OUTDIR + "/a.sup.galigns" );    }
          cout << Date( ) << ": ===== you can now use START=" << stage << " ====="
               << endl;    };
     auto ProbeWrite = [&]( const String& stage )
     {    if (WRITE)
          {    cout << Date( ) << ": writing files before stage " << stage << endl;
               String OUTDIR = DIR + "/" + stage + WRITE_SUB;
               Mkdir777(OUTDIR);
               BinaryWriter::writeFile( OUTDIR + "/a.sup", D );
               BinaryWriter::writeFile( OUTDIR + "/a.sup.inv", dinv );
               FindLines( D, dinv, dlines, MAX_CELL_PATHS, MAX_CELL_DEPTH );
               BinaryWriter::writeFile( OUTDIR + "/a.sup.lines", dlines );
               ReadPathVec dpaths;
               PlaceReads( hb, paths, dup, D, dpaths, True, False );
               PlaceReadsSmart(hb, paths, dup, D, dinv, dpaths, dlines, bci, False);
               dpaths.WriteAll( OUTDIR + "/a.dpaths" );
               {    VecULongVec dpaths_index;
                    invert( dpaths, dpaths_index, D.E( ) );
                    dpaths_index.WriteAll( OUTDIR + "/a.dpaths.index" );    }
               if (ALIGN)
               {    MasterVec< SerfVec<triple<int,int,int> > > alignsb;
                    alignsb.ReadAll( DIR + "/a.alignsb" );
                    MasterVec<SerfVec<refalign>> galigns;
                    RefAlign( genome, hb.K( ), hb.Edges( ), inv, D, dinv, dlines, 
                         alignsb, galigns, False, vec<int>( ) );
                    galigns.WriteAll( OUTDIR + "/a.sup.galigns" );    }    }
          cout << Date( ) << ": ===== you can now use START=" << stage << " ====="
               << endl;    };
     auto PlusWrite = [&]( const String& stage )
     {    if (WRITE)
          {    cout << Date( ) << ": writing files before stage " << stage << endl;
               String OUTDIR = DIR + "/" + stage + WRITE_SUB;
               Mkdir777(OUTDIR);
               BinaryWriter::writeFile( OUTDIR + "/a.sup" , D );
               BinaryWriter::writeFile( OUTDIR + "/a.sup.inv", dinv );    
               FindLines( D, dinv, dlines, MAX_CELL_PATHS, MAX_CELL_DEPTH );
               BinaryWriter::writeFile( OUTDIR + "/a.sup.lines", dlines );
               vec<int> llens;
               GetLineLengths( hb, D, dlines, llens );
               BinaryWriter::writeFile( OUTDIR + "/a.sup.llens", llens );
               vec< vec< pair<int,int> > > lhood;
               LineProx( hb, inv, ebcx, D, dinv, dlines, qept, lhood );
               BinaryWriter::writeFile( OUTDIR + "/a.sup.lhood", lhood );
               {    ReadPathVec dpaths;
                    PlaceReads( hb, paths, dup, D, dpaths, True, False );
                    dpaths.WriteAll( OUTDIR + "/a.dpaths" ); 
                    {    VecULongVec dpaths_index;
                         invert( dpaths, dpaths_index, D.E( ) );
                         dpaths_index.WriteAll( OUTDIR + "/a.dpaths.index" );    }
                    IntIndex dpaths_index( dpaths, D.E( ) );
                    vec<vec<pair<int,int>>> lbp;
                    BarcodePos( 
                         bc, hb, D, dinv, dpaths, dpaths_index, dlines, lbp, 0 );
                    BinaryWriter::writeFile( OUTDIR + "/a.sup.lbp", lbp );
                    MasterVec<SerfVec<pair<int,int>>> lbpx;
                    for ( auto x : lbp )
                    {    SerfVec<pair<int,int>> y( x.size( ) );
                         for ( int j = 0; j < x.isize( ); j++ )
                              y[j] = x[j];
                         lbpx.push_back(y);    }
                    lbpx.WriteAll( OUTDIR + "/a.sup.lbpx" );
                    vec<int> kmers( hb.E( ) );
                    #pragma omp parallel for
                    for ( int e = 0; e < hb.E( ); e++ )
                         kmers[e] = hb.Kmers(e);
                    vec<double> COV;
                    LineCN( kmers, lbpx, D, dlines, llens, COV );
                    BinaryWriter::writeFile( OUTDIR + "/a.sup.lcov", COV );    }
               if (ALIGN)
               {    MasterVec< SerfVec<triple<int,int,int> > > alignsb;
                    alignsb.ReadAll( DIR + "/a.alignsb" );
                    vec<vec< pair<int,int> >> linelocs( kmers.size( ) );
                    for ( int i = 0; i < dlines.isize( ); i++ )
                    for ( int j = 0; j < dlines[i].isize( ); j++ )
                    for ( int k = 0; k < dlines[i][j].isize( ); k++ )
                    for ( int l = 0; l < dlines[i][j][k].isize( ); l++ )
                    {    int d = dlines[i][j][k][l];
                         if ( D.O(d)[0] < 0 ) continue;
                         for ( int m = 0; m < D.O(d).isize( ); m++ )
                         {    int e = D.O(d)[m];
                              linelocs[e].push( i, j );    }    }
                    MasterVec< SerfVec< quad<int,Bool,ho_interval,ho_interval> > > 
                         view( dlines.size( ) );
                    #pragma omp parallel for
                    for ( int l = 0; l < dlines.isize( ); l++ )
                    {    View( l, hb.K( ), kmers, inv, D, dlines, 
                              linelocs, alignsb, view[l] );    }
                    view.WriteAll( OUTDIR + "/a.sup.view" );    }
               if (ALIGN)
               {    MasterVec< SerfVec<triple<int,int,int> > > alignsb;
                    alignsb.ReadAll( DIR + "/a.alignsb" );
                    MasterVec<SerfVec<refalign>> galigns;
                    RefAlign( genome, hb.K( ), hb.Edges( ), inv, D, dinv, dlines, 
                         alignsb, galigns, False, vec<int>( ) );
                    galigns.WriteAll( OUTDIR + "/a.sup.galigns" );    }    }
          cout << Date( ) << ": ===== you can now use START=" << stage << " ====="
               << endl;    };

     // ============================================================================

     // Convert closures to graph. (THIS COMMENT DOESN'T MAKE SENSE.)

     ReadPathVec dpaths;
     vec< vec< pair<int,int> > > lhood;
     const int MIN_LEN = 0;
     const Bool single = False;
     cout << Date( ) << ": reading supergraph" << endl;
     auto load = [&]( const String& start )
     {    String INDIR = DIR + "/" + start + READ_SUB;
          BinaryReader::readFile( INDIR + "/a.sup", &D );
          BinaryReader::readFile( INDIR + "/a.sup.inv", &dinv );    };
     if (REWRITE)
     {    load(START);
          PlusWrite(START);
          Scram(0);    }
     if ( startpos >= Position( stagelist, String("star") )
          && startpos <= Position( stagelist, String("tidy") ) )
     {    load(START);    }
     else if ( START == "report" )
     {    load( "final" );
          String fdir = DIR + "/final" + READ_SUB;
          BinaryReader::readFile( fdir + "/a.sup.lines", &dlines );
          BinaryReader::readFile( fdir + "/a.sup.lcov", &COV );
          if (ALIGN) alignsb.ReadAll( DIR + "/a.alignsb" );    }
     else if ( START == "build" || START == "surgery" )
     {    String INDIR = DIR + "/" + "build" + READ_SUB;
          BinaryReader::readFile( INDIR + "/a.sup", &D );
          BinaryReader::readFile( INDIR + "/a.sup.inv", &dinv );
          BinaryReader::readFile( INDIR + "/a.sup.lines", &dlines );
          BinaryReader::readFile( INDIR + "/a.sup.lhood", &lhood );    }
     else
     {    if ( START == "path" || START == "scaffold" ) load( "pre" );
          else
          {    load( "orig" );

               // Clean the graph.  Zippering is temporarily turned off because the
               // existing zippering doesn't respect the involution.  When we fix
               // this and turn zippering back on, we should also turn on the
               // zipper testing in Validate.

               cout << Date( ) << ": cleaning graph" << endl;
               Cleaner( hb, inv, paths, dup, D, dinv, dpaths, True );
               // Zipper( D, inv );
               String OUTDIR = DIR + "/pre" + WRITE_SUB;
               if (WRITE)
               {    cout << Date( ) << ": writing assembly" << endl;
                    Mkdir777(OUTDIR);
                    #pragma omp parallel sections
                    { 
                         #pragma omp section
                         {    BinaryWriter::writeFile( OUTDIR + "/a.sup", D );    }
                         #pragma omp section
                         {    BinaryWriter::writeFile( 
                                   OUTDIR + "/a.sup.inv", dinv );    }    }    }
               Validate( hb, inv, D, dinv );
               if (WRITE)
               {    HyperBasevectorX hbd;
                    SuperToSeqGraph( hb, D, hbd );
                    hbd.Edges( ).WriteAll( OUTDIR + "/a.sup.fastb" );    }

               // Make digraphEX version.

               /*
               cout << Date( ) << ": converting" << endl;
               vec<SerfVec<int>> edges( D.E( ) );
               for ( int e = 0; e < D.E( ); e++ )
               {    SerfVec<int>& x = edges[e];
                    x.resize( D.O(e).size( ) );
                    for ( int j = 0; j < D.O(e).isize( ); j++ )
                         x[j] = D.O(e)[j];    }
               digraphE<SerfVec<int>> DX( D, edges );
               digraphEX<SerfVec<int>> DY(DX);
               if (WRITE)
                    BinaryWriter::writeFile( OUTDIR + "/a.sup", DY );
               */

               cout << Date( ) << ": ===== you can now use START=scaffold ====="
                    << endl;    }

          // Place reads on the graph.

          PlaceReads( hb, paths, dup, D, dpaths, True, single );

          // Make scaffolds.

          String OUTDIR = DIR + "/" + "build" + WRITE_SUB;
          Mkdir777(OUTDIR);
          for ( int pass = 1; pass <= 2; pass++ )
          {    vec<int> lens( D.E( ), 0 );
               #pragma omp parallel for
               for ( int e = 0; e < D.E( ); e++ )
               {    if ( D.O(e)[0] < 0 ) continue;
                    for ( int j = 0; j < D.O(e).isize( ); j++ )
                         lens[e] += hb.Kmers( D.O(e)[j] );    }
               if ( pass == 2 )
               {    vec<int> dels;
                    for ( int v = 0; v < D.N( ); v++ )
                    for ( int j = 0; j < D.From(v).isize( ); j++ )
                    {    int w = D.From(v)[j], d = D.IFrom(v,j);
                         if ( D.To(v).nonempty( ) || !D.From(v).solo( ) ) continue;
                         if ( D.From(w).nonempty( ) || !D.To(w).solo( ) ) continue;
                         if ( v == w ) continue;
                         if ( lens[d] >= 320 ) continue;
                         dels.push_back(d);    }
                    RemoveDuff( hb, D, dinv, dels, True );
                    D.DeleteEdgesParallel(dels);
                    CleanupCore( D, dinv );
                    PlaceReads( hb, paths, dup, D, dpaths, True, single );    }
               String link_report;
               Scaffold( hb, inv, ebcx, D, dinv, dpaths, bid,
                    datasets, True, link_report, single );
               Ofstream( out, OUTDIR + "/link_report" );
               out << link_report;    }

          // Kill inversion artifacts.

          cout << Date( ) << ": recompute paths" << endl;
          PlaceReads( hb, paths, dup, D, dpaths, True, single );
          {    IntIndex dpaths_index( dpaths, D.E( ) );
               vec<int> dels;
               const int MAX_CAN_INS_DEL = 4;
               KillInversionArtifacts( 
                    D, dinv, dpaths, dpaths_index, bid, dels, MAX_CAN_INS_DEL );
               D.DeleteEdges(dels);    }
          RemoveUnneededVertices( D, dinv );
          CleanupCore( D, dinv );
            
          // Find lines.

          cout << Date( ) << ": finding lines" << endl;
          FindLines( D, dinv, dlines, MAX_CELL_PATHS, MAX_CELL_DEPTH );

          // Find proximate lines.

          cout << Date( ) << ": finding proximate lines" << endl;
          LineProx( hb, inv, ebcx, D, dinv, dlines, qept, lhood );

          // Compute line lengths.

          vec<int> llens;
          GetLineLengths( hb, D, dlines, llens );

          // Splay vertices at the ends of long lines.

          // PlusWrite( "presplay" ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          Splay( D, dinv, dlines, llens, MIN_SPLAY1 );

          // Write files, to permit reentry after this point.

          if (WRITE)
          {    cout << Date( ) << ": writing files before build" << endl;
               String OUTDIR = DIR + "/" + "build" + WRITE_SUB;
               Mkdir777(OUTDIR);
               BinaryWriter::writeFile( OUTDIR + "/a.sup", D );
               BinaryWriter::writeFile( OUTDIR + "/a.sup.inv", dinv );
               BinaryWriter::writeFile( OUTDIR + "/a.sup.lines", dlines );
               vec<int> llens;
               GetLineLengths( hb, D, dlines, llens );
               BinaryWriter::writeFile( OUTDIR + "/a.sup.llens", llens );
               BinaryWriter::writeFile( OUTDIR + "/a.sup.lhood", lhood );    }
          cout << Date( ) << ": ===== you can now use START=build ====="
               << endl;    
          // PlusWrite( "build" ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
               }
     Validate( hb, inv, D, dinv );

     // ============================================================================

     // Output data structures for surgery.

     vec< pair< int, vec<int> > > s1s2;
     vec< vec< triple< int, pair<int,int>, pair<int,int> > > > closures1;
     vec< vec< digraphE<vec<int>> > > closures2;
     vec< vec< vec<int> > > closures3;

     // Begin building to capture gaps.

     if ( startpos <= Position( stagelist, String("build") ) )
     {
          // Compute ancillary data structures.

          cout << Date( ) << ": computing ancillaries" << endl;
          vec<int> llens, linv, to_right;
          GetLineLengths( hb, D, dlines, llens );
          LineInv( dlines, dinv, linv );
          D.ToRight(to_right);

          // Place reads again, and index.

          PlaceReads( hb, paths, dup, D, dpaths, True, single );
          cout << Date( ) << ": indexing paths" << endl;
          VecULongVec dpaths_index;
          invert( dpaths, dpaths_index, D.E( ) );

          // Find reads that are internal to a line.  This is turned off because it
          // makes results worse and results in only a modest speedup.  Perhaps
          // it's because of misplaced reads.

          vec<Bool> internal( dpaths.size( ), False );
          /*
          cout << Date( ) << ": finding internals" << endl;
          const int END_SIZE = 5000;
          #pragma omp parallel for
          for ( int i = 0; i < dlines.isize( ); i++ )
          {    const vec<vec<vec<int>>>& L = dlines[i];
               int pos = 0;
               for ( int j = 0; j < L.isize( ); j++ )
               {    vec<int> lensj;
                    for ( int k = 0; k < L[j].isize( ); k++ )
                    {    int len = 0;
                         for ( int l = 0; l < L[j][k].isize( ); l++ )
                         {    int d = L[j][k][l];
                              if ( D.O(d)[0] < 0 ) continue;
                              if ( pos + len > END_SIZE
                                   && llens[i] - pos - len > END_SIZE )
                              {    for ( int z = 0;
                                        z < (int) dpaths_index[d].size( ); z++ )
                                   {    int64_t id1 = dpaths_index[d][z];
                                        const ReadPath& p = dpaths[id1];
                                        int offset = p.getOffset( );
                                        for ( int m = 0; m < (int) p.size( ); m++ )
                                        {    if ( p[m] == d ) break;
                                             for ( auto e : D.O( p[m] ) )
                                                  offset -= hb.Kmers(e);    }
                                        if ( llens[i] - pos - len - p.getOffset( )
                                             <= END_SIZE )
                                        {    continue;    }
                                        internal[id1] = True;
                                        int64_t id2
                                             = ( id1 % 2 == 0 ? id1+1 : id1-1 );
                                        internal[id2] = True;    }    }
                              for ( int m = 0; m < D.O(d).isize( ); m++ )
                              {    int e = D.O(d)[m];
                                   len += hb.Kmers(e);    }    }
                         lensj.push_back(len);    }
                    Sort(lensj);
                    if ( lensj.nonempty( ) ) pos += Median(lensj);    }    }
          */

          // Heuristics for closing.

          const int MIN_LINE_TO_WALK = 1000;
          const int NHOOD_DEPTH = 3;

          // Define the walk set s1s2.

          {    vec<int> l1s;
               for ( int l1 = 0; l1 < dlines.isize( ); l1++ )
               {    if ( llens[l1] < MIN_LINE_TO_WALK ) continue;
                    // test for very weird thing that might happen:
                    if ( D.O( dlines[l1].back( )[0][0] )[0] < 0 ) continue;
                    int d = dlines[l1].back( )[0][0];
                    int v = to_right[d];
                    if ( D.From(v).nonempty( ) || D.To(v).size( ) > 1 ) continue;
                    l1s.push_back(l1);    }
               for ( int li = 0; li < l1s.isize( ); li++ )
               {    int l1 = l1s[li];

                    // Define possible values for L2.

                    vec<int> l2s;
                    for ( int i = 0;
                         i < Min( NHOOD_DEPTH, lhood[l1].isize( ) ); i++ )
                    {    int l2 = lhood[l1][i].second;
                         int s2 = dlines[l2].front( )[0][0];
                         l2s.push_back( l2 );
                         int rs2 = dlines[ linv[l2] ].front( )[0][0];
                         l2s.push_back( linv[l2] );
                         UniqueSort(l2s);    }

                    // Check for overlaps.

                    Bool overlaps = False;
                    int s1 = dlines[l1].back( )[0][0];
                    for ( int i = 0; i < l2s.isize( ); i++ )
                    {    int s2 = dlines[ l2s[i] ].front( )[0][0];
                         if ( !IsUnique( vec<int>{ s1, s2, dinv[s1], dinv[s2] } ) )
                              overlaps = True;    }
                    if (overlaps)
                    {    // cout << Date( ) << ": overlaps, giving up" << endl;
                         continue;    }

                    // Save.

                    vec<int> s2s;
                    for ( int i = 0; i < l2s.isize( ); i++ )
                    {    int s2 = dlines[ l2s[i] ].front( )[0][0];
                         s2s.push_back(s2);    }
                    s1s2.push( s1, s2s );    }    }

          // Try to capture gaps.

          cout << Date( ) << ": start gap capturing loop" << endl;
          String OUTDIR = DIR + "/" + "star" + WRITE_SUB;
          Mkdir777(OUTDIR);
          Ofstream( gout, OUTDIR + "/capture.log" );
          cout << Date( ) << ": logging to " << OUTDIR << "/capture.log" << endl;
          const Bool use_rights = False;
          {    MasterVec< SerfVec<triple<int,int,int> > > alignsb;
               if (ALIGN) alignsb.ReadAll( DIR + "/a.alignsb" );
               Unvoid( s1s2, hb, inv, bci, dup, bad, paths, ebcx, alignsb,
                    D, dinv, dlines, dpaths, dpaths_index, use_rights, gout,
                    closures1, closures2, closures3, DIR, SAVE_LOCAL1, -1 );    }
          BinaryWriter::writeFile( OUTDIR + "/a.closures1", closures1 );
          BinaryWriter::writeFile( OUTDIR + "/a.closures2", closures2 );
          BinaryWriter::writeFile( OUTDIR + "/a.closures3", closures3 );
          BinaryWriter::writeFile( OUTDIR + "/a.s1s2", s1s2 );
          cout << Date( ) << ": ===== you can now use START=surgery =====" << endl;
          if (BUILD_ONLY) Scram(0);    }

     // ============================================================================

     // Carry out surgery.

     if ( startpos <= Position( stagelist, String("surgery") ) )
     {    if ( START == "surgery" )
          {    String OUTDIR = DIR + "/" + "star" + WRITE_SUB;
               if ( START == "surgery" ) OUTDIR = DIR + "/" + "star" + READ_SUB;
               BinaryReader::readFile( OUTDIR + "/a.closures1", &closures1 );
               BinaryReader::readFile( OUTDIR + "/a.closures2" , &closures2 );
               BinaryReader::readFile( OUTDIR + "/a.closures3", &closures3 );
               BinaryReader::readFile( OUTDIR + "/a.s1s2", &s1s2 );    }
          String OUTDIR = DIR + "/" + "star" + WRITE_SUB;
          Mkdir777(OUTDIR);
          {    Ofstream( out, OUTDIR + "/surgery.log" );
               Surgery( hb, inv, D, dinv, dlines, s1s2, closures1, closures2, 
                    closures3, False, False, out );    }
          Zipper( D, dinv );
          Validate( hb, inv, D, dinv );

          // Merge some unmerged stuff.  There are uncollapsed long perfect repeats
          // in the assembly at this point.  One solution (below) is to merge them
          // now.  Other possible solutions:
          // 1. Call Surgery with one of the require_simp* options on.
          // 2. Trace the duplication to its source and eliminate it.

          if (GLUE)
          {    vec<int> mult;
               ComputeMult( hb, D, mult );
               const int MIN_LINE_TO_WALK = 1000;
               const int MAX_MULT = 5;
               GlueAssemblies( hb, D, dinv, mult, MIN_LINE_TO_WALK, MAX_MULT );    }

          // Clean up.

          vec<int> dels;
          cout << Date( ) << ": placing reads again" << endl;
          PlaceReads( hb, paths, dup, D, dpaths, True, single );
          const int MAX_KILL = 350;
          const double MIN_RATIO = 25;
          Bool verbose = True;
          {    KillLowUnique( hb, D, dels, True );
               KillLowUniqueFrac( hb, D, dels );
               cout << Date( ) << ": deleting " << ToStringAddCommas( dels.size( ) )
                    << " edges by low unique" << endl;
               SimpleHangs( 
                    hb, D, dinv, dels, MAX_KILL, MIN_RATIO, verbose, single );
               IntIndex dpaths_index( dpaths, D.E( ) );
               const int MAX_CAN_INS_DEL = 4;
               KillInversionArtifacts( 
                    D, dinv, dpaths, dpaths_index, bid, dels, MAX_CAN_INS_DEL );    }

          // Remove compound hangs.

          const int MAX_TINY = 350;
          if (verbose) cout << Date( ) << ": computing lengths" << endl;
          vec<int> lens( D.E( ), 0 );
          #pragma omp parallel for
          for ( int e = 0; e < D.E( ); e++ )
          {    if ( D.O(e)[0] < 0 ) continue;
               for ( int j = 0; j < D.O(e).isize( ); j++ )
                    lens[e] += hb.Kmers( D.O(e)[j] );    }
          vec<int> dfw;
          DistancesToEndArr( D, lens, MAX_TINY * MIN_RATIO, True, dfw );
          FindCompoundHangs(
               D, dinv, lens, dfw, dels, MAX_TINY, MIN_RATIO, verbose );

          // Capture loops.

          CaptureCanonicalLoops( D, dinv, dels, verbose, single );
          CaptureSimpleLoops( D, dinv, dels, verbose, single );
          //clean up necessary after this step above
          UniqueSort(dels);
          D.DeleteEdgesParallel(dels);
          dels.clear();

          // Pull apart, then capture messy loops.

          PullApart( inv, D, dinv, dpaths, dels, verbose );
          UniqueSort(dels);
          D.DeleteEdgesParallel(dels);
          dels.clear( );
          CaptureMessyLoops( hb, inv, D, dinv, dels );

          // Complete clean up.

          UniqueSort(dels);
          D.DeleteEdgesParallel(dels);
          RemoveUnneededVertices( D, dinv );
          CleanupCore( D, dinv );
          Validate( hb, inv, D, dinv );
          FindLines( D, dinv, dlines, MAX_CELL_PATHS, MAX_CELL_DEPTH );
          LineProx( hb, inv, ebcx, D, dinv, dlines, qept, lhood );
          BasicWrite( "star" );    }

     // ============================================================================

     // Splay vertices at ends of long lines, then fix misassemblies.

     const int MIN_SPLAY2 = 5000;
     if ( startpos <= Position( stagelist, String("star") ) )
     {    cout << Date( ) << ": resplaying" << endl;
          FindLines( D, dinv, dlines, MAX_CELL_PATHS, MAX_CELL_DEPTH );
          vec<int> llens;
          GetLineLengths( hb, D, dlines, llens );
          Splay( D, dinv, dlines, llens, MIN_SPLAY2 );
          if (GLUE)
          {    vec<int> dels;
               KillLowUnique( hb, D, dels, True );
               cout << Date( ) << ": post splay cleanup, deleting "
                    << ToStringAddCommas( dels.size( ) ) << " edges" << endl;
               D.DeleteEdges(dels);
               RemoveUnneededVertices( D, dinv );
               CleanupCore( D, dinv );    }
          Validate( hb, inv, D, dinv );
          cout << Date( ) << ": placing reads" << endl;
          PlaceReads( hb, paths, dup, D, dpaths, True, False );
          cout << Date( ) << ": fixing misassemblies" << endl;
          FixMisassemblies( hb, dup, bci, paths, D, dinv, dpaths );
          CleanupCore( D, dinv );
          if (DEBUG) PlusWrite( "starstar" );
          else BasicWrite( "starstar" );    }

     // ============================================================================

     // Introduce barcode-only gaps.

     if ( startpos <= Position( stagelist, String("starstar") ) )
     {    const Bool DJANGO = False;
          Star( hb, inv, dup, bc, paths, D, dinv, dpaths, ebcx, qept,
               MIN_ADVANTAGE, OO, DJANGO );

          // Kill some gaps that look misassembled.

          cout << Date( ) << ": begin fixing misassemblies" << endl;
          PlaceReads( hb, paths, dup, D, dpaths, True, False );
          cout << Date( ) << ": making index" << endl;
          vec<int> dels;
          FindLines( D, dinv, dlines, MAX_CELL_PATHS, MAX_CELL_DEPTH );
          const int BC_REQUIRE = 5000;
          const int BC_FLANK = 20000;
          const int BC_IGNORE = 2000;
          vec <double> fhist;
          vec <pair<float, int>> lr;
          {    IntIndex dpaths_index( dpaths, D.E( ), False );
               KillMisassembledCells( hb, dup, bci, paths, D, dinv, dpaths, 
                    dpaths_index, bc, dlines, dels, fhist, lr, BC_REQUIRE, 
                    BC_FLANK, BC_IGNORE, False);    }
          // warn if the mean length-weighted molecule length is
          // less than 20 kb.
          {
               const double BIN_WIDTH=1000;
               double lw_mean_mol_length = 0.0;
               double totlen = 0.0;
               for ( int fbin = 0 ; fbin != int(fhist.size()); fbin++ ) {
                    double bin_len = (fbin + 0.5)*BIN_WIDTH;
                    totlen += fhist[fbin] * bin_len;
                    lw_mean_mol_length   += fhist[fbin] * bin_len * bin_len;
               }
               if ( fhist.size() > 0 )
                    lw_mean_mol_length /= totlen;
               // warn if fragments are too short.
               StatLogger::log( "lw_mean_mol_len", lw_mean_mol_length, 
                    "Length-weighted mean molecule length", true );
               StatLogger::issue_alert( "lw_mean_mol_len", lw_mean_mol_length );
               String OUTDIR = DIR + "/" + "final" + WRITE_SUB;
               Mkdir777(OUTDIR);
               WriteHistToJson( fhist, 0.0, BIN_WIDTH*fhist.size(), BIN_WIDTH,
                    OUTDIR+"/stats", "molecules", "patch" );
          }
          WriteFileStatsVecPairFloatInt( lr, "s.lr" );
          // Kill misassemblies.

          KillMisassembledCellsAlt( hb, D, dinv, ebcx, dlines, dels );
          UniqueSort(dels);
          cout << Date( ) << ": deleting " << dels.size( )
               << " putatatively bad edges" << endl;
          D.DeleteEdges(dels);
          CleanupCore( D, dinv );
          BasicWrite( "patch" );    }

     // ============================================================================

     // Patch bc-only gaps.

     if ( startpos <= Position( stagelist, String("patch") ) )
     {    cout << Date( ) << ": start stage to patch bc-only gaps" << endl;
          PlaceReads( hb, paths, dup, D, dpaths, True, single );
          s1s2.clear( );
          for ( int v = 0; v < D.N( ); v++ )
          {    if ( !D.From(v).solo( ) ) continue;
               int d = D.IFrom(v,0);
               if ( !IsBarcodeOnlyGap( D.O(d) ) ) continue;
               int w = D.From(v)[0];
               if ( !D.To(v).solo( ) || !D.From(w).solo( ) || !D.To(w).solo( ) )
                    continue;
               int d1 = D.ITo(v,0), d2 = D.IFrom(w,0);
               if ( D.O(d1)[0] < 0 || D.O(d2)[0] < 0 ) continue;
               s1s2.push( d1, vec<int>{d2} );    }
          cout << Date( ) << ": start gap filling loop" << endl;
          String OUTDIR = DIR + "/" + "phase" + WRITE_SUB;
          Mkdir777(OUTDIR);
          Ofstream( gout, OUTDIR + "/fill.log" );
          cout << Date( ) << ": logging to a.base/fill.log" << endl;
          // SHOULD MAKE TRUE? // **************************************************
          const Bool use_rights = False; // ****************************************
          const int MAX_READS = 1000000;
          FindLines( D, dinv, dlines, MAX_CELL_PATHS, MAX_CELL_DEPTH );
          cout << Date( ) << ": N50 line = " 
               << ToStringAddCommas( LineN50( hb, D, dlines ) ) << endl;
          {    MasterVec< SerfVec<triple<int,int,int> > > alignsb;
               if (ALIGN) alignsb.ReadAll( DIR + "/a.alignsb" );
               VecULongVec dpaths_index2;
               invert( dpaths, dpaths_index2, D.E( ) );
               Unvoid( s1s2, hb, inv, bci, dup, bad, paths, ebcx, alignsb, D,
                    dinv, dlines, dpaths, dpaths_index2, use_rights, gout,
                    closures1, closures2, closures3, DIR, SAVE_LOCAL2,
                    MAX_READS );    }
          {    Ofstream( out, OUTDIR + "/surgery.log" );
               Surgery( hb, inv, D, dinv, dlines, s1s2, closures1, closures2, 
                    closures3, True, REQUIRE_SIMPLER, out );    }
          Zipper( D, dinv );
          {    vec<int> dels;

               // Using crazy high threshold for hanging ends here.  What we
               // should be doing is folding these back.

               const int MAX_KILL = 15000;
               const double MIN_RATIO = 5;
               Bool verbose = True;
               SimpleHangs(
                    hb, D, dinv, dels, MAX_KILL, MIN_RATIO, verbose, single );
               D.DeleteEdges(dels);
               RemoveUnneededVertices( D, dinv );
               CleanupCore( D, dinv );
               Validate( hb, inv, D, dinv );

               // Fix misassemblies.

               PlaceReads( hb, paths, dup, D, dpaths, True, False );
               cout << Date( ) << ": making index" << endl;
               {    IntIndex dpaths_index( dpaths, D.E( ), False );
                    cout << Date( ) << ": fixing misassemblies" << endl;
                    dels.clear( );
                    FindLines( D, dinv, dlines, MAX_CELL_PATHS, MAX_CELL_DEPTH );
                    cout << Date( ) << ": N50 line = " 
                         << ToStringAddCommas( LineN50( hb, D, dlines ) ) << endl;
                    const int BC_REQUIRE = 5000;
                    const int BC_FLANK = 20000;
                    const int BC_IGNORE = 2000;
                    vec<double> fhist;
                    vec<pair<float,int>> lr;
                    KillMisassembledCells( hb, dup, bci, paths, D, dinv, dpaths,
                         dpaths_index, bc, dlines, dels, fhist, lr,
                         BC_REQUIRE, BC_FLANK, BC_IGNORE, False);    }
               cout << Date( ) << ": cleaning up" << endl;
               D.DeleteEdges(dels);
               CleanupCore( D, dinv );

               // Introduce barcode-only gaps.

               PlaceReads( hb, paths, dup, D, dpaths, True, False );
               const Bool DJANGO = True;
               Star( hb, inv, dup, bc, paths, D, dinv, dpaths, ebcx,
                    qept, MIN_ADVANTAGE, OO, DJANGO );
               Validate( hb, inv, D, dinv );    }

          // Add a final pass of breaking.

          cout << Date( ) << ": start final breaking" << endl;
          FindLines( D, dinv, dlines, MAX_CELL_PATHS, MAX_CELL_DEPTH );
          cout << Date( ) << ": N50 line = " 
               << ToStringAddCommas( LineN50( hb, D, dlines ) ) << endl;
          PlaceReads( hb, paths, dup, D, dpaths, True, False );
          IntIndex dpaths_index( dpaths, D.E( ), False );
          vec<int> dels;
          const int BC_REQUIRE2 = 5000;
          const int BC_FLANK2 = 20000;
          const int BC_IGNORE2 = 5000;
          vec <double> fhist;
          vec <pair<float, int>> lr;
          KillMisassembledCells( hb, dup, bci, paths, D, dinv, dpaths, dpaths_index,
               bc, dlines, dels, fhist, lr, BC_REQUIRE2, BC_FLANK2, BC_IGNORE2, False );
          D.DeleteEdges(dels);
          CleanupCore( D, dinv );
          Validate( hb, inv, D, dinv );

          // Another round of breaking.

          cout << Date( ) << ": start final final breaking" << endl;
          FindLines( D, dinv, dlines, MAX_CELL_PATHS, MAX_CELL_DEPTH );
          cout << Date( ) << ": N50 line = " 
               << ToStringAddCommas( LineN50( hb, D, dlines ) ) << endl;
          PlaceReads( hb, paths, dup, D, dpaths, True, False );
          dpaths_index.Initialize( dpaths, D.E( ), False );
          const int BC_REQUIRE3 = 25000;
          const int BC_FLANK3 = 40000;
          const int BC_IGNORE3 = 20000;
          dels.clear( );
          KillMisassembledCells( hb, dup, bci, paths, D, dinv, dpaths, dpaths_index,
               bc, dlines, dels, fhist, lr, BC_REQUIRE3, BC_FLANK3, BC_IGNORE3, False );
          D.DeleteEdges(dels);
          CleanupCore( D, dinv );
          Validate( hb, inv, D, dinv );
          BasicWrite( "phase" );    }

     // ============================================================================

     // Another round of filling.  If you want to test this round on a single gap,
     // set START = phase, SAVE_LOCAL2 = True, and S1 = the left edge.  This 
     // modifies some log files so you should set WRITE_SUB for such a test.

     if ( startpos <= Position( stagelist, String("phase") ) )
     {    cout << Date( ) << ": start second stage to patch bc-only gaps" << endl;
          String OUTDIR = DIR + "/" + "post" + WRITE_SUB;
          Mkdir777(OUTDIR);

          // Running two passes because in the first pass, if there are barcode
          // gaps on the left and the right of an edge, both gaps will not be
          // closed.  We try to avoid duplicating too many calculations.

          vec<vec<int>> D1;
          for ( int pass = 1; pass <= 2; pass++ )
          {    cout << Date( ) << ": start pass " << pass << endl;
               if ( pass == 1 ) 
               {    D1 = D.Edges( );
                    ParallelUniqueSort(D1);    }
               PlaceReads( hb, paths, dup, D, dpaths, True, single );
               s1s2.clear( );
               for ( int v = 0; v < D.N( ); v++ )
               {    if ( !D.From(v).solo( ) ) continue;
                    int d = D.IFrom(v,0);
                    if ( !IsBarcodeOnlyGap( D.O(d) ) ) continue;
                    int w = D.From(v)[0];
                    if ( !D.To(v).solo( ) || !D.From(w).solo( ) || !D.To(w).solo( ) )
                         continue;
                    int d1 = D.ITo(v,0), d2 = D.IFrom(w,0);
                    if ( D.O(d1)[0] < 0 || D.O(d2)[0] < 0 ) continue;
                    if ( S1 >= 0 && d1 != S1 ) continue;
                    if ( pass == 2 && BinMember( D1, D.O(d1) ) 
                         && BinMember( D1, D.O(d2) ) )
                    {    continue;    }
                    s1s2.push( d1, vec<int>{d2} );    }
               cout << Date( ) << ": start gap filling loop two" << endl;
               String logfile = OUTDIR + "/fill.log" + ToString(pass);
               Ofstream( gout, logfile );
               cout << Date( ) << ": logging to " << logfile << endl;
               const Bool use_rights = True;
               const int MAX_READS = 1000000;
               FindLines( D, dinv, dlines, MAX_CELL_PATHS, MAX_CELL_DEPTH );
               cout << Date( ) << ": N50 line = " 
                    << ToStringAddCommas( LineN50( hb, D, dlines ) ) << endl;
               {    MasterVec< SerfVec<triple<int,int,int> > > alignsb;
                    if (ALIGN) alignsb.ReadAll( DIR + "/a.alignsb" );
                    VecULongVec dpaths_index2;
                    invert( dpaths, dpaths_index2, D.E( ) );
                    Bool unvoid_verbose = SAVE_LOCAL2;
                    Unvoid( s1s2, hb, inv, bci, dup, bad, paths, ebcx, alignsb, D,
                         dinv, dlines, dpaths, dpaths_index2, use_rights, gout,
                         closures1, closures2, closures3, DIR, SAVE_LOCAL2,
                         MAX_READS, True, unvoid_verbose, NEW_TEST );    }
               {    Ofstream( out, OUTDIR + "/surgery.log" + ToString(pass) );
                    Surgery( hb, inv, D, dinv, dlines, s1s2, closures1, closures2, 
                         closures3, True, REQUIRE_SIMPLER, out );    }
               if (SAVE_LOCAL2)
               {    cout << "\n" << Date( ) << ": done, time used = " 
                         << TimeSince(clock) << ", peak mem = " 
                         << PeakMemUsageGBString( ) << endl;
                    Scram(0);    }
               Zipper( D, dinv );
               vec<int> dels;

               // Using crazy high threshold for hanging ends here.  What we
               // should be doing is folding these back.

               const int MAX_KILL = 15000;
               const double MIN_RATIO = 5;
               Bool verbose = True;
               SimpleHangs( 
                    hb, D, dinv, dels, MAX_KILL, MIN_RATIO, verbose, single );
               D.DeleteEdges(dels);
               RemoveUnneededVertices( D, dinv );
               CleanupCore( D, dinv );
               Validate( hb, inv, D, dinv );    }

          // Flatten some bubbles.

          const int MAX_DELTA = 2;
          const double MIN_RATIO = 6.0;
          const int MAX_DEL = 4;
          FlattenSomeBubbles( hb, dup, paths, D, dinv, dpaths, MAX_DELTA,
               MIN_RATIO, MAX_DEL );

          // Delete low unique stuff.

          vec<int> dels;
          KillLowUnique( hb, D, dels, True );
          KillLowUniqueFrac( hb, D, dels );
          UniqueSort(dels);
          cout << Date( ) << ": deleting " << ToStringAddCommas( dels.size( ) )
               << " low-uniq edges (" << PERCENT_RATIO( 3, dels.isize( ), D.E( ) )
               << ") before patching" << endl;
          D.DeleteEdgesParallel(dels);
          CleanupCore( D, dinv );

          // Patch in original gap closures.

          {    vecbasevector closures;
               cout << Date() << ": about to release paths, "
                    << "memory = " << MemUsageGBString() << endl;
               Destroy(paths);
               cout << Date() << ": released paths, memory = " <<
                    MemUsageGBString() << endl;
               BinaryReader::readFile( DIR + "/../closures.fastb", &closures );
               Splat( hb, inv, closures, D, dinv );
               cout << Date() << ": about to reload paths" << endl;
               double clock = WallClockTime();
               paths.readBinary( DIR+"/a.pathsX" );
               cout << Date() << ": reloaded paths in " 
                    << TimeSince(clock) << endl;    }
          BasicWrite( "post" );    }

     // ============================================================================

     // Next stage.

     if ( startpos <= Position( stagelist, String("post") ) )
     {    
          // Convert some barcode-only gaps to pair gaps.

          {    cout << Date( ) << ": converting gaps" << endl;
               PlaceReads( hb, paths, dup, D, dpaths, True, False );
               IntIndex dpaths_index( dpaths, D.E( ), False );
               vec<int> to_left, to_right;
               D.ToLeft(to_left), D.ToRight(to_right);
               int count = 0;
               #pragma omp parallel for schedule(dynamic,10)
               for ( int d = 0; d < D.E( ); d++ )
               {    if ( dinv[d] <= d || !IsBarcodeOnlyGap( D.O(d) ) ) continue;
                    int v = to_left[d], w = to_right[d];
                    if ( !D.To(v).solo( ) || !D.From(w).solo( ) ) continue;
                    int f = D.ITo(v,0), g = D.IFrom(w,0);
                    Bool linked = False;
                    for ( int i = 0; i < dpaths_index.Count(f); i++ )
                    {    int64_t id = dpaths_index.Val(f,i);
                         int64_t idp = ( id % 2 == 0 ? id+1 : id-1 );
                         for ( int j = 0; j < (int) dpaths[idp].size( ); j++ )
                         {    if ( dpaths[idp][j] == dinv[g] )
                              {    linked = True;
                                   break;    }    }
                         if (linked) break;    }
                    if (linked)
                    {    D.OMutable(d) = D.OMutable( dinv[d] ) = {-1};    
                         #pragma omp critical
                         {    count += 2;    }    }    }
               cout << Date( ) << ": converted " << count << " gaps" << endl;    }

          // Smart read placement. 

          cout << Date( ) << ": setting up for smart read placement" << endl;
          String OUTDIR = DIR + "/" + "fix" + WRITE_SUB;
          Mkdir777(OUTDIR);
          FindLines( D, dinv, dlines, MAX_CELL_PATHS, MAX_CELL_DEPTH );
          PlaceReads( hb, paths, dup, D, dpaths, True, False );
          PlaceReadsSmart( hb, paths, dup, D, dinv, dpaths, dlines, bci, False );
          dpaths.WriteAll( OUTDIR + "/a.xpaths" );
          {    VecULongVec dpaths_index;
               invert( dpaths, dpaths_index, D.E( ) );
               dpaths_index.WriteAll( OUTDIR + "/a.xpaths.index" );    }

          // Close gaps.

          cout << Date( ) << ": post patching" << endl;
          HyperBasevectorX ddn;
          VirtualMasterVec<basevector> bases(
               DIR + "/../data/frag_reads_orig.fastb" );
          VirtualMasterVec<PQVec> quals( DIR + "/../data/frag_reads_orig.qualp" );
          VirtualMasterVec<ReadPath> xpaths( OUTDIR + "/a.xpaths" );
          VirtualMasterVec<ULongVec> xpaths_index( OUTDIR + "/a.xpaths.index" );
          String R, S;
          Stackaroo( bases, quals, hb, inv, xpaths, xpaths_index, D, dinv, R, S,
               True, False, 0, ddn, False, False, False );
          RemoveDuff2( hb, D, dinv );
          BasicWrite( "fix" );    }

     // ============================================================================

     // Fix some inversions and do other stuff.

     if ( startpos <= Position( stagelist, String("fix") ) )
     {    
          // Now try to join lines.

          cout << Date( ) << ": trying to join lines" << endl;
          if (ULTRA)
          {    FindLines( D, dinv, dlines, MAX_CELL_PATHS, MAX_CELL_DEPTH );
               vec<int> llens;
               GetLineLengths( hb, D, dlines, llens );
               Splay( D, dinv, dlines, llens, MIN_SPLAY2 );    }
          PlaceReads( hb, paths, dup, D, dpaths, True, False );
          const Bool DJANGO = True;
          const double MIN_ADVANTAGE = 60.0;
          Star( hb, inv, dup, bc, paths, D, dinv, dpaths, ebcx,
               qept, MIN_ADVANTAGE, False, DJANGO );
          cout << Date( ) << ": begin fixing misassemblies" << endl;
          PlaceReads( hb, paths, dup, D, dpaths, True, False );

          // Try scaffolding.  This is commented out because it is expensive
          // and does not seem to change much.  The expensive part is at
          // "finding barcodes on lines" in Scaffold.cc.

          /*
          {    String link_report;
               Scaffold( hb, inv, ebcx, D, dinv, dpaths, bid,
                    datasets, True, link_report, single );    }
          RemoveDuff2( hb, D, dinv );
          PlaceReads( hb, paths, dup, D, dpaths, True, False );
          */

          // Kill misassemblies.

          cout << Date( ) << ": making index" << endl;
          FindLines( D, dinv, dlines, MAX_CELL_PATHS, MAX_CELL_DEPTH );
          const int BC_REQUIRE = 25000;
          const int BC_FLANK = 40000;
          const int BC_IGNORE = 20000;
          int mol_lwml = - 1;
          vec<int> dels;
          {    IntIndex dpaths_index( dpaths, D.E( ), False );
               vec<double> fhist;
               vec<pair<float,int>> lr;
               KillMisassembledCells( hb, dup, bci, paths, D, dinv, dpaths,
                    dpaths_index, bc, dlines, dels, fhist, lr, BC_REQUIRE, BC_FLANK,
                    BC_IGNORE, False );    }
          D.DeleteEdges(dels);
          RemoveUnneededVertices( D, dinv );
          CleanupCore( D, dinv );

          // Remove some cycles.

          {    Decycle( hb, inv, D, dinv, paths );
               const Bool align2 = True;
               PlaceReads( hb, paths, dup, D, dpaths, True, False, align2 );
               cout << Date( ) << ": deleting weak branches after decycle" << endl;
               const int MAX_WEAK3_LOSE = 2;
               const int MAX_WEAK3_LOSE_TOTAL = 2;
               const int MIN_WEAK3_WIN = 5;
               const int MIN_WEAK3_RATIO = 5;
               vec<int> dels;
               {    IntIndex dpaths_index( dpaths, D.E( ), False );
                    DelWeak5( D, dinv, dpaths, dpaths_index, dels,
                         MAX_WEAK3_LOSE, MAX_WEAK3_LOSE_TOTAL, MIN_WEAK3_WIN, 
                         MIN_WEAK3_RATIO, True );    }
               D.DeleteEdges(dels);
               RemoveUnneededVertices( D, dinv );
               CleanupCore( D, dinv );    }
          
          // Zap weird inversion bubbles.

          dels.clear( );
          ZapInversionBubbles( D, dinv, dels );
          D.DeleteEdges(dels);
          RemoveUnneededVertices( D, dinv );
          CleanupCore( D, dinv );    
          ZapMegaInversionBubbles( D, dinv );

          // Fix some inversions.

          {    FindLines( D, dinv, dlines, MAX_CELL_PATHS, MAX_CELL_DEPTH );
               PlaceReads( hb, paths, dup, D, dpaths, True, False );
               vec<vec<pair<int,int>>> lbp;
               {    IntIndex dpaths_index( dpaths, D.E( ) );
                    BarcodePos( bc, hb, D, dinv, 
                         dpaths, dpaths_index, dlines, lbp, 0 );    }
               MasterVec<SerfVec<pair<int,int>>> lbpx;
               for ( auto x : lbp )
               {    SerfVec<pair<int,int>> y( x.size( ) );
                    for ( int j = 0; j < x.isize( ); j++ )
                         y[j] = x[j];
                    lbpx.push_back(y);    }
               MasterVec<SerfVec<refalign>> galigns;
               if (ALIGN)
               {    MasterVec< SerfVec<triple<int,int,int> > > alignsb;
                    alignsb.ReadAll( DIR + "/a.alignsb" );
                    RefAlign( genome, hb.K( ), hb.Edges( ), inv, D, dinv, dlines, 
                         alignsb, galigns, False, vec<int>( ) );    }
               InvFix( hb, inv, kmers, D, dinv, dlines, lbpx, galigns );    }

          // Save.

          if (WRITE)
          {    cout << Date( ) << ": writing files before actual phasing" << endl;
               String OUTDIR = DIR + "/" + "fase" + WRITE_SUB;
               Mkdir777(OUTDIR);
               BinaryWriter::writeFile( OUTDIR + "/a.sup" , D );
               BinaryWriter::writeFile( OUTDIR + "/a.sup.inv", dinv );    
               FindLines( D, dinv, dlines, MAX_CELL_PATHS, MAX_CELL_DEPTH );
               BinaryWriter::writeFile( OUTDIR + "/a.sup.lines", dlines );
               vec<int> llens;
               GetLineLengths( hb, D, dlines, llens );
               BinaryWriter::writeFile( OUTDIR + "/a.sup.llens", llens );
               vec< vec< pair<int,int> > > lhood;
               LineProx( hb, inv, ebcx, D, dinv, dlines, qept, lhood );
               BinaryWriter::writeFile( OUTDIR + "/a.sup.lhood", lhood );
               {    PlaceReads( hb, paths, dup, D, dpaths, True, False );
                    dpaths.WriteAll( OUTDIR + "/a.dpaths" ); 
                    {    VecULongVec dpaths_index;
                         invert( dpaths, dpaths_index, D.E( ) );
                         dpaths_index.WriteAll( OUTDIR + "/a.dpaths.index" );    }
                    IntIndex dpaths_index( dpaths, D.E( ) );
                    vec<vec<pair<int,int>>> lbp;
                    BarcodePos( 
                         bc, hb, D, dinv, dpaths, dpaths_index, dlines, lbp, 0 );
                    BinaryWriter::writeFile( OUTDIR + "/a.sup.lbp", lbp );
                    MasterVec<SerfVec<pair<int,int>>> lbpx;
                    for ( auto x : lbp )
                    {    SerfVec<pair<int,int>> y( x.size( ) );
                         for ( int j = 0; j < x.isize( ); j++ )
                              y[j] = x[j];
                         lbpx.push_back(y);    }
                    lbpx.WriteAll( OUTDIR + "/a.sup.lbpx" );
                    vec<int> kmers( hb.E( ) );
                    #pragma omp parallel for
                    for ( int e = 0; e < hb.E( ); e++ )
                         kmers[e] = hb.Kmers(e);
                    vec<double> COV;
                    LineCN( kmers, lbpx, D, dlines, llens, COV );
                    BinaryWriter::writeFile( OUTDIR + "/a.sup.lcov", COV );    }
               if (ALIGN)
               {    MasterVec< SerfVec<triple<int,int,int> > > alignsb;
                    alignsb.ReadAll( DIR + "/a.alignsb" );
                    MasterVec<SerfVec<refalign>> galigns;
                    RefAlign( genome, hb.K( ), hb.Edges( ), inv, D, dinv, dlines, 
                         alignsb, galigns, False, vec<int>( ) );
                    galigns.WriteAll( OUTDIR + "/a.sup.galigns" );
                    vec<vec< pair<int,int> >> linelocs( kmers.size( ) );
                    for ( int i = 0; i < dlines.isize( ); i++ )
                    for ( int j = 0; j < dlines[i].isize( ); j++ )
                    for ( int k = 0; k < dlines[i][j].isize( ); k++ )
                    for ( int l = 0; l < dlines[i][j][k].isize( ); l++ )
                    {    int d = dlines[i][j][k][l];
                         if ( D.O(d)[0] < 0 ) continue;
                         for ( int m = 0; m < D.O(d).isize( ); m++ )
                         {    int e = D.O(d)[m];
                              linelocs[e].push( i, j );    }    }
                    MasterVec< SerfVec< quad<int,Bool,ho_interval,ho_interval> > > 
                         view( dlines.size( ) );
                    #pragma omp parallel for
                    for ( int l = 0; l < dlines.isize( ); l++ )
                    {    View( l, hb.K( ), kmers, inv, D, dlines, 
                              linelocs, alignsb, view[l] );    }
                    view.WriteAll( OUTDIR + "/a.sup.view" );    }    

               // Estimate heterozygosity rate.
               // NOTE.  This looks for bubbles in the unphased assembly.  We should
               // be restricting attention to bubbles having good coverage on both 
               // branches.

               // Heuristics.

               const int MIN_LINE = 100000;
               const int MIN_COV = 10;
               const int UNSAMPLE = 100;

               // Compute edge coverage.

               vec<int> cov( D.E( ), 0 );
               for ( auto p : dpaths )
               for ( auto d : p ) cov[d]++;

               // Now actually estimate heterozygosity rate.

               cout << Date( ) << ": computing, sampling 1/" << UNSAMPLE 
                    << " bubbles" << endl;
               int64_t total_length = 0, total_snps = 0;
               vec<pair<int,int>> all_bubbles;
               #pragma omp parallel for schedule(dynamic, 10000)
               for ( int i = 0; i < dlines.isize( ); i++ )
               {    const vec<vec<vec<int>>>& L = dlines[i];
                    int pos = 0;
                    vec<pair<int,int>> bubbles;
                    for ( int j = 0; j < L.isize( ); j++ )
                    {    if ( L[j].size( ) == 2 && L[j][0].solo( ) && L[j][1].solo( )
                              && D.O( L[j][0][0] )[0] >= 0 
                              && D.O( L[j][1][0] )[0] >= 0 )
                         {    int d1 = L[j][0][0], d2 = L[j][1][0];
                              if ( cov[d1] + cov[dinv[d1]] >= MIN_COV 
                                   && cov[d2] + cov[dinv[d2]] >= MIN_COV )
                              {    if ( ( i + j + d1 + d2 ) % UNSAMPLE == 0 )
                                        bubbles.push( d1, d2 );    }    }
                         vec<int> lensj;
                         for ( int k = 0; k < L[j].isize( ); k++ )
                         {    int len = 0;
                              for ( int l = 0; l < L[j][k].isize( ); l++ )
                              {    int d = L[j][k][l];
                                   if ( D.O(d)[0] < 0 ) continue;
                                   for ( int m = 0; m < D.O(d).isize( ); m++ )
                                        len += hb.Kmers( D.O(d)[m] );    }
                              lensj.push_back(len);    }
                         Sort(lensj);
                         if ( lensj.nonempty( ) ) pos += Median(lensj);    }
                    if ( pos < MIN_LINE ) continue;
                    #pragma omp critical
                    {    total_length += pos;
                         all_bubbles.append(bubbles);    }    }
               cout << Date( ) << ": there are " 
                    << ToStringAddCommas( all_bubbles.size( ) ) 
                    << " bubbles" << endl;
               const int batch = 1000;
               #pragma omp parallel for
               for ( int bi = 0; bi < all_bubbles.isize( ); bi += batch )
               {    int snps = 0;
                    for ( int i = bi; 
                         i < Min( bi + batch, all_bubbles.isize( ) ); i++ )
                    {    int d1 = all_bubbles[i].first, d2 = all_bubbles[i].second;
                         basevector x1 = hb.Cat( D.O(d1) ), x2 = hb.Cat( D.O(d2) );
                         alignment al;
                         SmithWatAffine( x1, x2, al );
                         align a(al);
                         vec<int> mgg = a.MutationsGap1Gap2( x1, x2 );
                         snps += mgg[0];    }
                    #pragma omp critical
                    {    total_snps += snps;    }    }
               int hetdist = 0;
               if ( total_snps != 0 )
                    hetdist = int( round( total_length / (total_snps*UNSAMPLE) ) );
               cout << Date( ) << ": estimated heterozygosity rate = 1 / "
                    << ToStringAddCommas(hetdist) << endl;
               WriteFileStatsInt( hetdist, "s.hetdist" );    }
          cout << Date( ) << ": ===== you can now use START=fase ====="
               << endl;    }

     // ============================================================================

     // Phase the assembly.

     if ( startpos <= Position( stagelist, String("fase") ) )
     {    
          // Estimate sizes of barcode-only gaps.

          FindLines( D, dinv, dlines, MAX_CELL_PATHS, MAX_CELL_DEPTH );
          {    MasterVec<SerfVec<refalign>> galigns;
               Bool EVAL = False;
               MasterVec<SerfVec<pair<int,int>>> lbpx;
               String XDIR;
               if ( startpos < Position( stagelist, String("fase") ) )
                    XDIR = DIR + "/" + "fase" + WRITE_SUB;
               else XDIR = DIR + "/" + "fase" + READ_SUB;
               lbpx.ReadAll( XDIR + "/a.sup.lbpx" );
               Gaprika( kmers, D, dinv, dlines, lbpx, MAX_GAP,
                    galigns, 0, EVAL );    }

          // Now phase.

          cout << Date( ) << ": start phasing" << endl;
          FindLines( D, dinv, dlines, MAX_CELL_PATHS, MAX_CELL_DEPTH );
          PlaceReads( hb, paths, dup, D, dpaths, True, False );
          // does not seem to help:
          // PlaceReadsSmart( hb, paths, dup, D, dinv, dpaths, dlines, bci, False );
          {    IntIndex dpaths_index( dpaths, D.E( ), False );
               vec<String> report;
               Flipper( 
                    hb, inv, dup, bad, bc, D, dinv, dlines, dpaths_index, report );
               Mkdir777( DIR + "/final" );
               // NOT THE RIGHT PLACE TO PUT THIS
               Ofstream( out, DIR + "/final/o.phase" );
               for ( int l = 0; l < report.isize( ); l++ )
                    out << report[l];    }
          BasicWrite( "tidy" );    }

     // ============================================================================

     // Tidy the assembly.

     if ( startpos <= Position( stagelist, String("tidy") ) )
     {
          // Some cleanup.  This is only needed because of artifacts introduced
          // by Local.  First set up.

          vec<int> to_left, to_right;
          D.ToLeft(to_left), D.ToRight(to_right);
          FindLines( D, dinv, dlines, MAX_CELL_PATHS, MAX_CELL_DEPTH );
          vec< pair<int,int> > tol( D.E( ), make_pair(-1,-1) );
          for ( int i = 0; i < dlines.isize( ); i++ )
          for ( int j = 0; j < dlines[i].isize( ); j++ )
          for ( int k = 0; k < dlines[i][j].isize( ); k++ )
          for ( int l = 0; l < dlines[i][j][k].isize( ); l++ )
          {    int e = dlines[i][j][k][l];
               if ( e >= 0 ) tol[e] = make_pair( i, j );    }
          vec<int> dlens( D.E( ), 0 );
          #pragma omp parallel for
          for ( int e = 0; e < D.E( ); e++ )
          {    if ( D.O(e)[0] < 0 ) continue;
               for ( int j = 0; j < D.O(e).isize( ); j++ )
                    dlens[e] += hb.Kmers( D.O(e)[j] );    }

          // Suppose a line ends in a barcode gap.  Try to fix.
          // If barcode and pair gap are adjacent, merge.

          cout << Date( ) << ": identify edits" << endl;
          Validate( hb, inv, D, dinv );
          const int MAX_TO_DIE = 500;
          vec<int> dels;
          int countp = 0;
          for ( int v = 0; v < D.N( ); v++ )
          {    if ( !D.To(v).solo( ) || !D.From(v).solo( ) ) continue;
               int d1 = D.ITo(v,0), d2 = D.IFrom(v,0);
               int i1 = tol[d1].first, j1 = tol[d1].second;
               int i2 = tol[d2].first, j2 = tol[d2].second;
               if ( i1 < 0 || i2 < 0 ) continue;
               if ( i1 == i2 ) continue;
               if ( D.O(d1)[0] >= 0 || D.O(d2)[0] >= 0 ) continue;
               int rd1 = dinv[d1], rd2 = dinv[d2];
               int ri1 = tol[rd1].first, rj1 = tol[rd1].second;
               int ri2 = tol[rd2].first, rj2 = tol[rd2].second;
               if ( ri1 < 0 || ri2 < 0 ) continue;
               if ( ri1 == ri2 ) continue;
               if ( D.O(rd1)[0] >= 0 || D.O(rd2)[0] >= 0 ) continue;
               if ( IsBarcodeOnlyGap( D.O(d1) ) && IsPairGap( D.O(d2) )
                    && IsBarcodeOnlyGap( D.O(rd2) ) && IsPairGap( D.O(rd1) ) )
               {    countp++;
                    int v = to_left[d2], w = to_right[d2];
                    D.TransferEdgesWithUpdate( w, v, to_left, to_right );
                    dels.push_back(d2);    }
               else if ( IsBarcodeOnlyGap( D.O(d2) ) && IsPairGap( D.O(d1) )
                    && IsBarcodeOnlyGap( D.O(rd1) ) && IsPairGap( D.O(rd2) ) )
               {    countp++;
                    int v = to_left[d1], w = to_right[d1];
                    D.TransferEdgesWithUpdate( v, w, to_left, to_right );
                    dels.push_back(d1);    }
               else if ( IsBarcodeOnlyGap( D.O(d1) ) && j1 % 2 == 0 && j1 >= 2
                    && IsBarcodeOnlyGap( D.O(rd2) ) && dlines[ri2].size( ) >= 2 )
               {    vec<int> ds = Contents( dlines[i1][j1-1] );
                    int n = 0;
                    for ( auto d : ds ) n += dlens[d];
                    if ( n <= MAX_TO_DIE )
                    {    countp++;
                         dels.append(ds);
                         int f = dlines[i1][j1-2][0][0];
                         int v = to_right[f], w = to_left[d1];
                         D.TransferEdgesWithUpdate( 
                              v, w, to_left, to_right );    }    }
               else if ( IsBarcodeOnlyGap( D.O(d2) ) && dlines[i2].size( ) >= 2
                    && IsBarcodeOnlyGap( D.O(rd1) ) && rj1 % 2 == 0 && rj1 >= 2 )
               {    vec<int> ds = Contents( dlines[i2][1] );
                    int n = 0;
                    for ( auto d : ds ) n += dlens[d];
                    if ( n <= MAX_TO_DIE )
                    {    countp++;
                         dels.append(ds);
                         int f = dlines[i2][2][0][0];
                         int v = to_left[f], w = to_right[d2];
                         D.TransferEdgesWithUpdate( 
                              v, w, to_left, to_right );    }    }    }
          cout << Date( ) << ": made " << countp << " edits" << endl;
          cout << Date( ) << ": deleting " << dels.size( ) << " edges" << endl;
          D.DeleteEdges(dels);
          CleanupCore( D, dinv );
          Validate( hb, inv, D, dinv );

          // Delete weak edges.  When two passes were ran, GAP increased by 0.6%,
          // which doesn't really make sense.  To investigate.

          for ( int pass = 1; pass <= 1; pass++ )
          {    cout << Date( ) << ": deleting weak branches" << endl;
               const int MAX_WEAK3_LOSE = 2;
               const int MAX_WEAK3_LOSE_TOTAL = 2;
               const int MIN_WEAK3_WIN = 5;
               const int MIN_WEAK3_RATIO = 5;
               FindLines( D, dinv, dlines, MAX_CELL_PATHS, MAX_CELL_DEPTH );
               PlaceReads( hb, paths, dup, D, dpaths, True, False );
               PlaceReadsSmart(hb, paths, dup, D, dinv, dpaths, dlines, bci, False);
               {    IntIndex dpaths_index( dpaths, D.E( ), False );
                    dels.clear( );
                    DelWeak4( D, dinv, dpaths, dpaths_index, dels,
                         MAX_WEAK3_LOSE, MAX_WEAK3_LOSE_TOTAL, MIN_WEAK3_WIN, 
                         MIN_WEAK3_RATIO, True );    }
               vec<Bool> to_delete( D.E( ), False );
               for ( int i = 0; i < dels.isize( ); i++ )
                    to_delete[ dels[i] ] = True;
               vec<int> dels2;
               for ( int i = 0; i < dlines.isize( ); i++ )
               {    const vec<vec<vec<int>>>& L = dlines[i];
                    for ( int j = 1; j < L.isize( ) - 1; j += 2 )
                    {    vec<int> all, keep;
                         for ( int k = 0; k < L[j].isize( ); k++ )
                              all.append( L[j][k] );
                         UniqueSort(all);
                         vec<Bool> del( L[j].size( ), False );
                         for ( int k = 0; k < L[j].isize( ); k++ )
                         {    for ( int l = 0; l < L[j][k].isize( ); l++ )
                              {    if ( to_delete[ L[j][k][l] ] )
                                   {    del[k] = True;
                                        break;    }    }    }
                         Bool good_path = False;
                         for ( int k = 0; k < L[j].isize( ); k++ )
                         {    if ( del[k] ) continue;
                              Bool gap = False;
                              for ( int l = 0; l < L[j][k].isize( ); l++ )
                                   if ( D.O( L[j][k][l] )[0] < 0 ) gap = True;
                              if ( !gap ) good_path = True;    }
                         if ( !good_path ) continue;
                         for ( int k = 0; k < L[j].isize( ); k++ )
                              if ( !del[k] ) keep.append( L[j][k] );
                         UniqueSort(keep);
                         for ( int k = 0; k < all.isize( ); k++ )
                         {    if ( !BinMember( keep, all[k] ) )
                              {    dels2.push_back( 
                                        all[k], dinv[ all[k] ] );    }    }    }    }
               cout << Date( ) << ": actually deleting " << dels2.size( ) << " edges"
                    << endl;
               D.DeleteEdges(dels2);
               RemoveUnneededVertices( D, dinv );
               CleanupCore( D, dinv );    }

          // Delete 3:0 bubbles.

          cout << Date( ) << ": deleting 3:0 bubbles" << endl;
          // const int MAX_DELTA = 2;
          const double MIN_RATIO = 3.0;
          const int MAX_DEL = 0;
          FindLines( D, dinv, dlines, MAX_CELL_PATHS, MAX_CELL_DEPTH );
          PlaceReads( hb, paths, dup, D, dpaths, True, False );
          PlaceReadsSmart(hb, paths, dup, D, dinv, dpaths, dlines, bci, False);
          dels.clear( );
          {    IntIndex dpaths_index( dpaths, D.E( ), False );
               vec<int> lens( D.E( ), 0 );
               #pragma omp parallel for
               for ( int e = 0; e < D.E( ); e++ )
               {    if ( D.O(e)[0] < 0 ) continue;
                    for ( int j = 0; j < D.O(e).isize( ); j++ )
                         lens[e] += hb.Kmers( D.O(e)[j] );    }
               #pragma omp parallel for
               for ( int v = 0; v < D.N( ); v++ )
               {    if ( D.To(v).size( ) != 1 || D.From(v).size( ) != 2 ) continue;
                    if ( D.From(v)[0] != D.From(v)[1] ) continue;
                    int w = D.From(v)[0];
                    if ( v == w ) continue;
                    if ( D.To(w).size( ) != 2 || D.From(w).size( ) != 1 ) continue;
                    if ( D.O( D.IFrom(v,0) )[0] < 0 || D.O( D.IFrom(v,1) )[0] < 0 ) 
                         continue;
                         int delta = 
                              AbsDiff( lens[ D.IFrom(v,0)], lens[ D.IFrom(v,1)] );
                    // if ( delta == 0 || delta > MAX_DELTA ) continue;
                    vec<vec<int64_t>> supp(2);
                    for ( int j = 0; j < 2; j++ )
                    {    int d = D.IFrom(v,j);
                         int rd = dinv[d];
                         for ( int i = 0; i < dpaths_index.Count(d); i++ )
                              supp[j].push_back( dpaths_index.Val(d,i) / 2 );
                         for ( int i = 0; i < dpaths_index.Count(rd); i++ )
                              supp[j].push_back( dpaths_index.Val(rd,i) / 2 );
                         UniqueSort( supp[j] );    }
                    int m = MeetSize( supp[0], supp[1] );
                    for ( int j = 0; j < 2; j++ )
                    {    int d1 = D.IFrom(v,j);
                         int n1 = supp[j].isize( ) - m, n2 = supp[1-j].isize( ) - m;
                         if ( n1 > MAX_DEL || n2 < MIN_RATIO * Max(1, n1) ) continue;
                         #pragma omp critical
                         {    dels.push_back( d1, dinv[d1] );    }    }    }    }
          UniqueSort(dels);
          cout << Date( ) << ": deleting " << ToStringAddCommas( dels.size( ) )
               << " weak bubble edges" << endl;
          D.DeleteEdges(dels);
          RemoveUnneededVertices( D, dinv );
          CleanupCore( D, dinv );

          // Delete terminal sequence gap edges.  These should be very rare, but
          // will make the downstream code crash.  
          // WE MIGHT WANT TO PUT THIS IN RELEASE 1.0.1.

          dels.clear( );
          D.ToLeft(to_left), D.ToRight(to_right);
          #pragma omp parallel for
          for ( int d = 0; d < D.E( ); d++ )
          {    if ( !IsSequence( D.O(d) ) ) continue;
               int v = to_left[d], w = to_right[d];
               if ( D.To(v).empty( ) || D.From(w).empty( ) )
               {
                    #pragma omp critical
                    {    dels.push_back( d, dinv[d] );    }    }    }
          D.DeleteEdges(dels);
          RemoveUnneededVertices( D, dinv );
          CleanupCore( D, dinv );
          BasicWrite( "canon" );    }

     // ============================================================================

     // Canonicalize cells.

     if ( startpos <= Position( stagelist, String("canon") ) )
     {    FindLines( D, dinv, dlines, MAX_CELL_PATHS, MAX_CELL_DEPTH );
          cout << Date( ) << ": canonicalizing cells" << endl;
          const int XMAX_CANON = 4;
          vec<int> to_left, to_right;
          D.ToLeft(to_left), D.ToRight(to_right);
          vec<int> all_dels;
          for ( int i = 0; i < dlines.isize( ); i++ )
          {    const vec<vec<vec<int>>>& L = dlines[i];
               for ( int j = 1; j < L.isize( ) - 1; j += 2 )
               {    if ( L[j].isize( ) <= 2 || L[j].isize( ) > XMAX_CANON ) continue;
                    if ( L[j-1].empty( ) || L[j+1].empty( ) ) continue;
                    if ( L[j-1][0].empty( ) || L[j+1][0].empty( ) ) continue;
                    Bool have_gap = False;
                    vec<int> dels = Contents( L[j] );
                    for ( auto d : dels ) if ( D.O(d)[0] < 0 ) have_gap = True;
                    if (have_gap) continue;
                    int d1 = L[j-1][0][0], d2 = L[j+1][0][0];
                    int v = to_right[d1], w = to_left[d2];
                    int rd1 = dinv[d2], rd2 = dinv[d1];
                    if ( !IsUnique( vec<int>{ d1, d2, rd1, rd2 } ) ) continue;
                    if ( make_pair(rd1,rd2) < make_pair(d1,d2) ) continue;
                    int rv = to_right[rd1], rw = to_left[rd2];
                    int nd = dels.size( ), E = D.E( ), n = L[j].size( );
                    for ( int i = 0; i < nd; i++ ) dels.push_back( dinv[ dels[i] ] );
                    all_dels.append(dels);
                    vec<vec<int>> new_edges(n);
                    for ( int k = 0; k < n; k++ )
                    for ( int l = 0; l < L[j][k].isize( ); l++ )
                         new_edges[k].append( D.O( L[j][k][l] ) );
                    for ( int k = 0; k < n; k++ ) D.AddEdge( v, w, new_edges[k] );
                    for ( int i = 0; i < n; i++ )
                    {    new_edges[i].ReverseMe( );
                         for ( int k = 0; k < new_edges[i].isize( ); k++ )
                              new_edges[i][k] = inv[ new_edges[i][k] ];    }
                    for ( int k = 0; k < n; k++ ) D.AddEdge( rv, rw, new_edges[k] );
                    for ( int k = 0; k < n; k++ ) dinv.push_back( E + n + k );
                    for ( int k = 0; k < n; k++ ) dinv.push_back( E + k );    }    }
          D.DeleteEdges(all_dels);
          RemoveUnneededVertices( D, dinv );
          CleanupCore( D, dinv );
          Validate( hb, inv, D, dinv );

          // Write final assembly.  At this point we assume only that D and dinv
          // are updated.

          if (ALIGN) alignsb.ReadAll( DIR + "/a.alignsb" );
          SuperFiles( hb, inv, bci, bc, paths, dup, D, dinv, dpaths, 
               dlines, COV, ebcx, qept, genome, alignsb, DIR, WRITE_SUB );    }

     // Report some stats.

     cout << Date( ) << ": reporting assembly stats" << endl;
     Destroy(paths);
     {    ostringstream xout;

          // Reload dpaths.  This is terrible, just to compute fraction of reads
          // placed, should have computed this earlier and saved.

          String OUTDIR = DIR + "/" + "final" + WRITE_SUB;
          if ( dpaths.size( ) == 0 ) dpaths.ReadAll( OUTDIR + "/a.dpaths" );

          // OK now report the stats.

          ReportAssemblyStats( bci, genome, ambint, hb, inv, D, dinv, dlines, dpaths,
               COV, alignsb, G, perc_proper_pairs, xout, DIR, WRITE_SUB, 
               CS_SAMPLE_ID, CS_SAMPLE_DESC );
          cout << xout.str( );
          Mkdir777(OUTDIR);
          Ofstream( yout, OUTDIR + "/stats/summary.txt" );
          yout << xout.str( );    }

     // Copy over all stats into one place

     String stats_dir = DIR + "/../stats" + WRITE_SUB;
     Mkdir777(stats_dir);
     // q-score histogram
     if ( IsRegularFile( DIR + "/../data/frag_reads_orig.qhist" ) ) // TEMP!
          Cp2( DIR + "/../data/frag_reads_orig.qhist", stats_dir );
     // insert size distribution
     if ( IsRegularFile( DIR + "/a.ins_dist" ) ) // TEMP!
          Cp2( DIR + "/a.ins_dist", stats_dir);

     // all stats 
     // first copy over histograms
     vec<String> histograms ({"contig", "edge", "phase_block", "scaffold", 
          "reads_per_barcode", "molecules"});
     for ( auto & name : histograms ) {
          String fn = DIR + "/final" + WRITE_SUB + "/stats/histogram_" + name+ ".json";
          if ( IsRegularFile ( fn ) )
               Cp2( fn, stats_dir);
     }
     Cp2( DIR + "/final" + WRITE_SUB + "/all_stats.json", stats_dir);
     Mv( DIR + "/final" + WRITE_SUB + "/summary.json", stats_dir + "/summary.json");         // ligo is saddened by the copy
     Cp2( DIR + "/final" + WRITE_SUB + "/summary_cs.csv", stats_dir);
     Cp2( DIR + "/final" + WRITE_SUB + "/stats/summary.txt", stats_dir);
     
     // Record etime, peak mem
     double etime_cp = (WallClockTime( ) - clock)/(60.*60.);
     StatLogger::log( "etime_cp_h", etime_cp, "Elapsed time CP" );
     StatLogger::log( "mem_peak_cp_gb", PeakMemUsageGB(), "Mem peak CP (gb)" );
     
     double etime_df = StatLogger::getNumStat("etime_df_h");
     double etime_h = etime_df + etime_cp;
     StatLogger::log( "etime_h", etime_h, "DF+CP time" );

     StatLogger::dump_text( DIR + "/final" + WRITE_SUB + "/stats/statistics.txt");
     StatLogger::dump_text( DIR + "/../stats/perf.stats" );

     // Done.

     cout << "\n" << Date( ) << ": done, time used = " << TimeSince(clock)
          << ", peak mem = " << PeakMemUsageGBString( ) << endl;

#ifdef JEMALLOC_HOOKS
     if ( JE_STATS ) (void) je_malloc_stats_print(nullptr, nullptr, nullptr);
#endif

     Scram(0);    }
