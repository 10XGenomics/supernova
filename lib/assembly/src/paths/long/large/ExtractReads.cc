///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "Bitvector.h"
#include "CoreTools.h"
#include "FastIfstream.h"
#include "FetchReads.h"
#include "PairsManager.h"
#include "ParseSet.h"
#include "Qualvector.h"
#include "TokenizeString.h"
#include "bam/ReadBAM.h"
#include "feudal/ObjectManager.h"
#include "feudal/PQVec.h"
#include "math/HoInterval.h"
#include "paths/long/LoadCorrectCore.h"
#include "paths/long/large/ExtractReads.h"
#include "paths/long/large/ReadNameLookup.h"
#include "random/Shuffle.h"

class rs_meta { // read set meta info
     public:
     String type;
     String sample;
     String lib;
     double frac;
     rs_meta( ) : type("frag"), sample("C"), lib("1"), frac(1) { }
     friend ostream& operator<<( ostream& out, const rs_meta& m )
     {    return out << "type=" << m.type << ",sample=" << m.sample << ",lib=" 
               << m.lib << ",frac=" << m.frac;    }
};

void GetAmbInt( const vecbitvector& amb, vec< pair<int,ho_interval> >& amb_int );

void GetCannedReferenceSequences( const String& sample, const String& species,
     const String& work_dir );

void ExtractReads( const String& sample, const String& species, String reads,
     String& SELECT_FRAC, const int READS_TO_USE, const vec<String>& regions, 
     const String& tmp_dir1, const String& work_dir, const Bool all, 
     const Bool USE_PF_ONLY, const Bool KEEP_NAMES, vec<String>& subsam_names, 
     vec<int64_t>& subsam_starts, vecbvec* pReads, ObjectManager<VecPQVec>& quals )
{
     double lclock = WallClockTime( );

     // this needs to be done to remove any symlinks created in a previous "all"
     // run, because a later non-all run will then over-write the symlinked data,
     // rather than the symlink itself.

     Remove( tmp_dir1 + "/frag_reads_orig.fastb" );
     Remove( tmp_dir1 + "/frag_reads_orig.qualp" );

     // Get reference sequence.

     if (all) GetCannedReferenceSequences( sample, species, work_dir );

     // Parse, check and load files.  This is the main extraction path.

     if ( reads != "" )
     {    cout << Date( ) << ": finding input files" << endl;
          reads.GlobalReplaceBy( " ", "" );
          vec< vec< String > > infiles;
          vec< vec< String > > infiles_rn;
          vecbasevector& xbases = (*pReads);
          VecPQVec& xquals = quals.create( );
          vecString xnames;
          vec< vec< pair<int,int> > > infiles_pairs;
          vec<rs_meta> infiles_meta;

          // Define read groups.
     
          vec<String> groups;
          String line;
          if ( reads.Contains( "@", 0 ) )
          {    String fn = reads.After( "@" );
               if ( !IsRegularFile(fn) )
               {    cout << "\nCan't find the file " << fn << " given by your "
                         << "READS argument.\n" << endl;
                    Scram(1);    }
               fast_ifstream in(fn);
               while(1)
               {    getline( in, line );
                    if ( in.fail( ) ) break;
                    groups.push_back(line);    }    }
          else Tokenize( reads, '+', groups );

          // Extract metainfo.

          infiles_meta.resize( groups.size( ) );
          for ( int g = 0; g < groups.isize( ); g++ )
          {    String meta;
               if ( groups[g].Contains( "::" ) )
               {    meta = groups[g].Before( "::" );
                    groups[g] = groups[g].After( "::" );    }
               vec<String> parts;
               Tokenize( meta, ',', parts );
               rs_meta& x = infiles_meta[g];
               for ( int i = 0; i < parts.isize( ); i++ )
               {    if ( !parts[i].Contains( ":" ) )
                    {    cout << "\nIllegal metainfo specification " << meta
                              << ".\nEach metainfo field must be of the form "
                              << "arg:value.\n" << endl;
                         Scram(1);    }
                    String arg = parts[i].Before( ":" ), val = parts[i].After( ":" );
                    if ( arg == "type" )
                    {    if ( val == "frag" || val == "long" || val == "jump" ) 
                              x.type = val;
                         else
                         {    cout << "\nUnrecognized type " << val 
                                   << " in metainfo specification " << meta
                                   << ".\n" << endl;
                              Scram(1);    }
                         if ( val != "frag" )
                         {    cout << "\nCurrently only type frag is implemented."
                                   << "\n" << endl;
                              Scram(1);    }    }
                    else if ( arg == "sample" ) x.sample = val;
                    else if ( arg == "lib" ) x.lib = val;
                    else if ( arg == "frac" ) 
                    {    if ( !val.IsDouble( ) || val.Double( ) <= 0
                              || val.Double( ) > 1 )
                         {    cout << "\nIllegal value " << val << " for frac in "
                                   << "metainfo specification " << meta
                                   << ".\n" << endl;
                              Scram(1);    }
                         x.frac = val.Double( );    }
                    else 
                    {    cout << "\nIllegal argument " << arg << " in "
                              << "metainfo specification " << meta
                              << ".\n" << endl;
                         Scram(1);    }    }    }
          // SHOULD ELIMINATE ******************************************************
          if ( SELECT_FRAC.IsDouble( ) && SELECT_FRAC.Double( ) < 1 ) 
          {    for ( int g = 0; g < groups.isize( ); g++ )
                    infiles_meta[g].frac = SELECT_FRAC.Double( );    }
     
          // Sort groups by sample.

          subsam_names.clear( );
          vec<String> ssn;
          for ( int g = 0; g < groups.isize( ); g++ )
          {    subsam_names.push_back( infiles_meta[g].sample );
               if ( !Member( ssn, subsam_names.back( ) ) )
                    ssn.push_back( subsam_names.back( ) );    }
          vec<int> ssi( subsam_names.size( ) );
          for ( int i = 0; i < subsam_names.isize( ); i++ )
               ssi[i] = Position( ssn, subsam_names[i] );
          SortSync( ssi, subsam_names, groups, infiles_meta );
          subsam_starts.resize_and_set( subsam_names.size( ), 0 );

          // Create file list.

          infiles.resize( groups.size( ) );
          infiles_rn.resize( groups.size( ) );
          infiles_pairs.resize( groups.size( ) );
          for ( int g = 0; g < groups.isize( ); g++ )
          {    String s = groups[g];

               // Parse data.

               vec<String> fns;
               int bcount = 0;
               String b;
               for ( int i = 0; i < s.isize( ); i++ )
               {    if ( s[i] == '{' ) 
                    {    b.push_back( s[i] );
                         bcount++;    }
                    else if ( s[i] == '}' ) 
                    {    b.push_back( s[i] );
                         bcount--;    }
                    else if ( s[i] == ',' && bcount == 0 )
                    {    fns.push_back(b);
                         b.clear( );    }
                    else b.push_back( s[i] );    }
               fns.push_back(b);
               for ( int i = 0; i < fns.isize( ); i++ )
               {    String f = fns[i];
                    vec<String> fs;
                    int ok = Glob( f, fs );
                    if ( ok != 0 )
                    {    cout << "\nFailed to glob " << f << ".\n"
                              << "This means that it does not correspond to a "
                              << "file or files according to the\n"
                              << "rules for globbing.  "
                              << "Please see the DISCOVAR de novo manual.\n"
                              << endl;
                         Scram(1);    }
                    infiles[g].append(fs);    }    }

          // Check that files are OK.

          for ( int g = 0; g < groups.isize( ); g++ )
          for ( int j = 0; j < infiles[g].isize( ); j++ )
          {    String fn = infiles[g][j];
               if ( !IsRegularFile(fn) )
               {    cout << "\nFailed to find file " << fn << ".\n" << endl;
                    Scram(1);    }
               vec<String> suf 
                    = { ".bam", ".fastq", ".fastq.gz", ".fq", ".fq.gz", ".fastb" };
               Bool ok = False;
               for ( auto s : suf )
                    if ( fn.Contains( s, -1 ) ) ok = True;
               if ( !ok ) 
               {    cout << "\nInput file " << fn << " has unsupported type.\n" 
                         << endl;
                    Scram(1);    }    }
          vec<String> flat;
          for ( int g = 0; g < groups.isize( ); g++ )
               flat.append( infiles[g] );
          Sort(flat);
          for ( int i = 1; i < flat.isize( ); i++ )
          {    if ( flat[i] == flat[i-1] || flat[i] == flat[i-1] + ".gz"
                    || flat[i] + ".gz" == flat[i-1] )
               {    cout << "\nFile " << flat[i] 
                         << " appears more than once in your "
                         << "READS specification.\n" << endl;
                    Scram(1);    }    }

          // Get first readnames for fastq files, and sort by them.

          for ( int g = 0; g < groups.isize( ); g++ )
          {    infiles_rn[g].resize( infiles[g].size( ) );
               for ( int j = 0; j < infiles[g].isize( ); j++ )
               {    String fn = infiles[g][j];
                    vec<String> suf = { ".fastq", ".fastq.gz", ".fq", ".fq.gz" };
                    Bool fq = False;
                    for ( auto s : suf )
                         if ( fn.Contains( s, -1 ) ) fq = True;
                    if ( !fq ) continue;
                    String command = "cat " + fn + " | head -1";
                    if ( fn.Contains( ".gz", -1 ) ) command = "z" + command;
                    fast_pipe_ifstream in(command);
                    getline( in, line );
                    if ( !line.Contains( "@", 0 ) || line.size( ) == 1
                         || ( line[1] == ' ' || line[1] == '/' ) )
                    {    cout << "\nSomething is wrong with the first line of your "
                              << "fastq file " << fn << "\n" << endl;
                         Scram(1);    }
                    int p = 0;
                    for ( p = 0; p < line.isize( ); p++ )
                         if ( line[p] == ' ' || line[p] == '/' ) break;
                    infiles_rn[g][j] = line.substr( 1, p - 1 );    }
               SortSync( infiles_rn[g], infiles[g] );
               for ( int j = 0; j < infiles_rn[g].isize( ); j++ )
               {    int k = infiles_rn[g].NextDiff(j);
                    if ( k - j > 2 && infiles_rn[g][j] != "" )
                    {    cout << "\nThere are more than two fastq files that "
                              << "start with the read name " << infiles_rn[g][j]
                              << ":\n";
                         for ( int l = j; l < k; l++ )
                              cout << "[" << l-j+1 << "] " << infiles[g][l] << endl;
                         cout << "Therefore it's not clear how to pair the "
                              << "files.\n\n";
                         Scram(1);    }    }    }

          // Check whether use of KEEP_NAMES makes sense.

          if (KEEP_NAMES)
          {    Bool have_bam = False, have_non_bam = False;
               for ( int g = 0; g < groups.isize( ); g++ )
               for ( int j = 0; j < infiles[g].isize( ); j++ )
               {    if ( infiles[g][j].Contains( ".bam", -1 ) ) have_bam = True;
                    else have_non_bam = True;    }
               if ( have_bam && have_non_bam )
               {    cout << "\nFor KEEP_NAMES, you can't supply a mixture of "
                         << "bam and non-bam files.\n" << endl;
                    Scram(1);    }    }

          // Precheck fastb files.

          for ( int g = 0; g < groups.isize( ); g++ )
          {    if ( infiles_meta[g].type != "long" )
               {    for ( int j = 0; j < infiles[g].isize( ); j++ )
                    {    String fn = infiles[g][j];
                         if ( !fn.Contains( ".fastb", -1 ) ) continue;
                         int64_t B = MastervecFileObjectCount(fn), Q;
                         if ( B % 2 != 0 )
                         {    cout << "\nThe file " << fn << " should be interlaced "
                                   << "and hence have an even number\n"
                                   << "of entries.  It does not.\n" << endl;
                              Scram(1);    }
                         String fn2b = fn.RevBefore( ".fastb" ) + ".qualb";
                         String fn2p = fn.RevBefore( ".fastb" ) + ".qualp";
                         String both;
                         if ( IsRegularFile(fn2b) )
                         {    Q = MastervecFileObjectCount(fn2b);
                              both = fn.Before( ".fastb" ) + ".{fastb,qualb}";    }
                         else if ( IsRegularFile(fn2p) )
                         {    Q = MastervecFileObjectCount(fn2p);
                              both = fn.Before( ".fastb" ) + ".{fastb,qualp}";    }
                         else
                         {    cout << "\nFor the file " << fn << ", there is no "
                                   << "matching qualb or qualp file.\n" << endl;
                              Scram(1);    }
                         if ( B != Q )
                         {    cout << "\nSee unequal file size for " << both
                                   << ".\n" << endl;
                              Scram(1);    }    }    }    }

          // Read the files.

          int nfiles = 0;
          for ( int g = 0; g < groups.isize( ); g++ )
               nfiles += infiles[g].size( );
          cout << Date( ) << ": reading " << nfiles 
               << " files (which may take a while)" << endl;
          for ( int g = 0; g < groups.isize( ); g++ )
          {    if ( g > 0 && subsam_names[g] != subsam_names[g-1] )
                    subsam_starts[g] = xbases.size( );
               for ( int j = 0; j < infiles[g].isize( ); j++ )
               {    String fn = infiles[g][j];
                    int64_t N0 = xbases.size( );
               
                    // Parse bam files.

                    if ( fn.Contains( ".bam", -1 ) )
                    {    bool const UNIQUIFY_NAMES = true;
                         vecString* pxnames = ( KEEP_NAMES ? &xnames : 0 );
                         BAMReader bamReader( USE_PF_ONLY, UNIQUIFY_NAMES,
                              infiles_meta[g].frac, long(READS_TO_USE) );
                              bamReader.readBAM( 
                                   fn, &xbases, &xquals, pxnames );    }

                    // Parse fastb/qualb/qualp files.

                    else if ( fn.Contains( ".fastb", -1 ) )
                    {    String fn2b = fn.RevBefore( ".fastb" ) + ".qualb";
                         String fn2p = fn.RevBefore( ".fastb" ) + ".qualp";
                         if ( IsRegularFile(fn2b) )
                         {    xbases.ReadAll( fn, True );
                              vecqualvector q;
                              q.ReadAll(fn2b);
                              convertAppendParallel( q.begin( ), q.end( ), xquals );
                              infiles[g][j] 
                                   = fn.Before( ".fastb" ) + ".{fastb,qualb}";    }
                         else if ( IsRegularFile(fn2p) )
                         {    xbases.ReadAll( fn, True );
                              xquals.ReadAll( fn2p, True );
                              infiles[g][j] = fn.Before( ".fastb" ) 
                                   + ".{fastb,qualp}";    }
                         double frac = infiles_meta[g].frac;
                         if ( frac < 1 )
                         {    int64_t total = 0, taken = 0;
                              Bool skip_next = False;
                              int64_t pos = N0;
                              for (int64_t i = N0; i < (int64_t) xbases.size( ); i++)
                              {    total++;
                                   if (skip_next)
                                   {    skip_next = False;
                                        continue;    }
                                   if ( total % 2 == 1 
                                        && double(taken)/double(total) > frac ) 
                                   {    skip_next = True;
                                        continue;    }
                                   taken++;
                                   if ( pos < i )
                                   {    xbases[pos] = xbases[i];
                                        xquals[pos] = xquals[i];    }
                                   pos++;    }
                              xbases.resize(pos), xquals.resize(pos);    }    }

                    // Parse paired fastq files.

                    else if ( infiles_rn[g][j] != "" 
                         && j < infiles_rn[g].isize( ) - 1
                         && infiles_rn[g][j] == infiles_rn[g][j+1] )
                    {    infiles_pairs[g].push( j, j+1 );
                         const String &fn1 = infiles[g][j], &fn2 = infiles[g][j+1];
                         String command1 = "cat " + fn1, command2 = "cat " + fn2;
                         if ( fn1.Contains( ".gz", -1 ) ) command1 = "z" + command1;
                         if ( fn2.Contains( ".gz", -1 ) ) command2 = "z" + command2;
                         fast_pipe_ifstream in1(command1), in2(command2);
                         String line1, line2;
                         int64_t total = 0, taken = 0;
                         double frac = infiles_meta[g].frac;

                         // Buffer for quality score compression in batches.
                         
                         const int qbmax = 10000000;
                         vec<qvec> qualsbuf;
                         MempoolOwner<char> alloc;
                         for ( int i = 0; i < qbmax; i++ )
                              qualsbuf.emplace_back(alloc);
                         int qbcount = 0;

                         // Go through the input files.

                         basevector b1, b2;
                         while(1)
                         {    getline(in1,line1), getline(in2,line2);
                              if ( in1.fail( ) && in2.fail( ) ) break;
                              if ( ( in1.fail( ) && !in2.fail( ) )
                                   || ( in2.fail( ) && !in1.fail( ) ) )
                              {    cout << "\nThe files " << fn1 << " and " << fn2 
                                        << " appear to be paired, yet have "
                                        << "different numbers of records.\n" << endl;
                                   Scram(1);    }
     
                              // Fetch bases.  Turn Ns into As.
     
                              getline( in1, line1 ), getline( in2, line2 );
                              if ( in1.fail( ) || in2.fail( ) )
                              {    cout << "\nSee incomplete record in " << fn1
                                        << " or " << fn2 << ".\n" << endl;
                                   Scram(1);    }
                              for ( int i = 0; i < line1.isize( ); i++ )
                                   if ( line1[i] == 'N' ) line1[i] = 'A';
                              for ( int i = 0; i < line2.isize( ); i++ )
                                   if ( line2[i] == 'N' ) line2[i] = 'A';
                              b1.SetFromString(line1);
                              b2.SetFromString(line2);
     
                              // Skip line.

                              getline( in1, line1 ), getline( in2, line2 );
                              if ( in1.fail( ) || in2.fail( ) )
                              {    cout << "\nSee incomplete record in " << fn1
                                        << " or " << fn2 << ".\n" << endl;
                                   Scram(1);    }
     
                              // Fetch quals.
     
                              getline( in1, line1 ), getline( in2, line2 );
                              if ( in1.fail( ) || in2.fail( ) )
                              {    cout << "\nSee incomplete record in " << fn1
                                        << " or " << fn2 << ".\n" << endl;
                                   Scram(1);    }
                              if ( b1.size( ) != line1.size( ) 
                                   || b2.size( ) != line2.size( ) )
                              {    cout << "\n1: " << b1.size( ) << " bases "
                                        << ", " << line1.size( ) << " quals" << endl;
                                   cout << "2: " << b2.size( ) << " bases "
                                        << ", " << line2.size( ) << " quals" << endl;
                                   cout << "See inconsistent base/quality lengths "
                                        << "in " << fn1 << " or " << fn2 << endl;
                                   Scram(1);    }

                              // Check frac.

                              if ( frac < 1 )
                              {    total++;
                                   if ( double(taken)/double(total) > frac ) 
                                        continue;
                                   taken++;    }

                              // Save.

                              qvec& q1 = qualsbuf[qbcount++];
                              qvec& q2 = qualsbuf[qbcount++];
                              q1.resize( line1.size( ) ), q2.resize( line2.size( ) );
                              if ( qbcount == qbmax )
                              {    convertAppendParallel( qualsbuf.begin( ), 
                                        qualsbuf.begin( ) + qbcount, xquals );
                                   qbcount = 0;     }
                              for ( int i = 0; i < line1.isize( ); i++ )
                                   q1[i] = line1[i] - 33;
                              for ( int i = 0; i < line2.isize( ); i++ )
                                   q2[i] = line2[i] - 33;
                              xbases.push_back(b1), xbases.push_back(b2);    }
                         convertAppendParallel( qualsbuf.begin( ), 
                              qualsbuf.begin( ) + qbcount, xquals );
                         j++;    }
               
                    // Parse unpaired fastq files.

                    else if ( infiles_rn[g][j] != "" )
                    {    vecqualvector Q;
                         const String& fn = infiles[g][j];
                         String command = "cat " + fn;
                         if ( fn.Contains( ".gz", -1 ) ) command = "z" + command;
                         fast_pipe_ifstream in(command);
                         int64_t total = 0, taken = 0;
                         double frac = infiles_meta[g].frac;
                         Bool skip_next = False;
                         while(1)
                         {    getline( in, line );
                              if ( in.fail( ) ) break;
          
                              // Fetch bases.  Turn Ns into As.
          
                              getline( in, line );
                              if ( in.fail( ) )
                              {    cout << "\nSee incomplete record in " << fn 
                                        << ".\n" << endl;
                                   Scram(1);    }
                              for ( int i = 0; i < line.isize( ); i++ )
                                   if ( line[i] == 'N' ) line[i] = 'A';
                              basevector b(line);
     
                              // Skip line.
          
                              getline( in, line );
                              if ( in.fail( ) )
                              {    cout << "\nSee incomplete record in " << fn 
                                        << ".\n" << endl;
                                   Scram(1);    }
     
                              // Fetch quals.
     
                              getline( in, line );
                              if ( in.fail( ) )
                              {    cout << "\nSee incomplete record in " << fn 
                                        << ".\n" << endl;
                                   Scram(1);    }
                              if ( b.size( ) != line.size( ) )
                              {    cout << "\nSee " << b.size( ) << " bases "
                                        << ", " << line.size( ) << " quals" << endl;
                                   cout << "See inconsistent base/quality lengths "
                                        << "in " << fn << ".\n" << endl;
                                   Scram(1);    }

                              // Check frac.

                              if ( frac < 1 )
                              {    total++;
                                   if (skip_next)
                                   {    skip_next = False;
                                        continue;    }
                                   if ( total % 2 == 1 
                                        && double(taken)/double(total) > frac ) 
                                   {    skip_next = True;
                                        continue;    }
                                   taken++;    }

                              // Save.

                              qualvector q( line.size( ) );
                              for ( int i = 0; i < line.isize( ); i++ )
                                   q[i] = line[i] - 33;
                              xbases.push_back(b);
                              Q.push_back(q);    }
     
                         // Check sanity and compress.

                         if ( infiles_meta[g].type != "long" && Q.size( ) % 2 != 0 )
                         {    cout << "\nThe file\n" << fn 
                                   << "\nshould be interlaced "
                                   << "and hence have an even number of entries."
                                   << "  It does not.\n" << endl;
                              Scram(1);    }
                         convertAppendParallel( 
                              Q.begin( ), Q.end( ), xquals );    }    }    }

          // Generate file list.

          vec<String> f1;
          vec< pair<String,String> > f2;
          vec<rs_meta> f1_meta, f2_meta;
          for ( int g = 0; g < infiles.isize( ); g++ )
          {    vec<Bool> used( infiles[g].size( ), False );
               for ( int j = 0; j < infiles_pairs[g].isize( ); j++ )
               {    int x1 = infiles_pairs[g][j].first; 
                    int x2 = infiles_pairs[g][j].second;
                    f2.push( infiles[g][x1], infiles[g][x2] );
                    f2_meta.push_back( infiles_meta[g] );
                    used[x1] = used[x2] = True;    }
               for ( int j = 0; j < infiles[g].isize( ); j++ )
               {    if ( !used[j] ) 
                    {    f1.push_back( infiles[g][j] );
                         f1_meta.push_back( infiles_meta[g] );    }    }    }
          cout << "\nINPUT FILES:\n";
          ostringstream iout;
          for ( int j = 0; j < f1.isize( ); j++ )
               iout << "[" << j+1 << "," << f1_meta[j] << "]  " << f1[j] << endl;
          for ( int j = 0; j < f2.isize( ); j++ )
          {    iout << "[" << f1.isize( ) + j + 1 << "a" 
                    << "," << f2_meta[j] << "] " << f2[j].first << endl;
               iout << "[" << f1.isize( ) + j + 1 << "b" 
                    << "," << f2_meta[j] << "] " << f2[j].second << endl;    }
          {    Ofstream( ioutx, work_dir + "/input_files" );
               ioutx << iout.str( );    }
          cout << iout.str( ) << "\n";

          // Fix subsam.

          int scount = 0;
          for ( int i = 0; i < subsam_names.isize( ); i++ )
          {    if ( i > 0 )
               {    if ( subsam_names[i] != subsam_names[i-1] )
                    {    subsam_names[++scount] = subsam_names[i];
                         subsam_starts[scount] = subsam_starts[i];    }    }    }
          subsam_names.resize(scount+1), subsam_starts.resize(scount+1);
          cout << Date( ) << ": found " << subsam_names.size( )
               << " samples" << endl;
          cout << Date( ) << ": starts = " << printSeq(subsam_starts) << endl;
          for ( int i = 0; i < subsam_starts.isize( ); i++ )
          {    if ( ( i < subsam_starts.isize( ) - 1 
                         && subsam_starts[i] == subsam_starts[i+1] )
                    || ( i == subsam_starts.isize( ) - 1
                         && subsam_starts[i] == (int64_t) xbases.size( ) ) )
               {    cout << "\nWARNING: It looks like you've got zero reads "
                         << "for sample " << subsam_names[i] << ".\n";
                    cout << "One way this could happen, for example, is if you "
                         << "provided a BAM file\nin which none of the reads were "
                         << "paired." << endl;    }    }

          // Work around bug.

          /*
          for ( int64_t i = 0; i < (int64_t) xbases.size( ); i++ )
          {    if ( xbases[i].size( ) < 60 )
               {    cout << "See read of length " << xbases[i].size( ) << ".  ";
                    cout << "For now, all reads must have length >= 60.\n" << endl;
                    Scram(1);    }    }
          */

          // Check for no reads.

          if ( xbases.size( ) == 0 )
          {    cout << "\nThere are no reads." << endl;
               cout << "Assembly cannot proceed, but please have an A1 day!"
                    << endl << endl;
               Scram(1);    }

          // Save files.
     
          xbases.WriteAll( work_dir + "/data/frag_reads_orig.fastb" );
          quals.store( );
          if ( xnames.size( ) > 0 )
          {    xnames.WriteAll( work_dir + "/data/frag_reads_orig.names" );
               readname_lookup look(xnames);
               BinaryWriter::writeFile( 
                    work_dir + "/data/frag_reads_orig.names.idx", look );    }

          // Report stats.

          int64_t all_reads =
               MastervecFileObjectCount( tmp_dir1 + "/frag_reads_orig.fastb" );
          cout << Date( ) << ": using " << ToStringAddCommas(all_reads)
               << " reads" << endl;
          cout << Date( ) << ": data extraction complete" 
               << ", peak mem = " << PeakMemUsageGBString( ) << endl;
          cout << TimeSince(lclock) << " used extracting reads" << endl;
          return;    }

     // Internal processing when reads are not specified.

     if (all)
     {    
          // Handle the fastb/qualb case and fastb/qualp cases.

          if ( sample == "maize" )
          {    String dir = "/wga/scr4/wg_projects/Z.mays/CML247";
               SystemSucceed( "ln -s " + dir + "/maize.fastb "
                     + tmp_dir1 + "/frag_reads_orig.fastb" );
               SystemSucceed( "ln -s " + dir + "/maize.qualp "
                     + tmp_dir1 + "/frag_reads_orig.qualp" );    }
          else if ( sample == "HCC1143+BL" || sample == "HCC1954+BL" 
               || sample == "F1" || sample == "F2" || sample == "F3"
               || sample == "CEPH" || sample == "YRI" )
          {    vec<String> sfss;
               vec<double> sfs;
               Tokenize( SELECT_FRAC, '+', sfss );
               sfs.resize( sfss.size( ) );
               for ( int j = 0; j < sfss.isize( ); j++ )
                    sfs[j] = sfss[j].Double( );
               int ns = sfs.size( );
               vec<String> dir(ns);
               String base = "/wga/scr4/wg_projects/H.sapien/";
               if ( sample == "HCC1143+BL" )
               {    dir[0] = base + "HCC1143";
                    dir[1] = base + "HCC1143BL";    }
               if ( sample == "HCC1954+BL" )
               {    dir[0] = base + "HCC1954";
                    dir[1] = base + "HCC1954BL";    }
               if ( sample == "F1" )
               {    dir[0] = base + "17E_PD";
                    dir[1] = base + "16E_MD";
                    dir[2] = base + "15E_DD";    }
               if ( sample == "F2" )
               {    dir[0] = base + "F2.1";
                    dir[1] = base + "F2.2";
                    dir[2] = base + "F2.3";    }
               if ( sample == "F3" )
               {    dir[0] = base + "F3.1";
                    dir[1] = base + "F3.2";
                    dir[2] = base + "F3.3";
                    dir[3] = base + "F3.4";
                    dir[4] = base + "F3.5";    }
               if ( sample == "CEPH" )
               {    dir[0] = base + "NA12878_CEPH";
                    dir[1] = base + "NA12892_CEPH";
                    dir[2] = base + "NA12891_CEPH";    }
               if ( sample == "YRI" )
               {    dir[0] = base + "YRI.1";
                    dir[1] = base + "YRI.2";
                    dir[2] = base + "YRI.3";    }
               vec<int64_t> all_reads(ns), using_reads(ns);
               for ( int j = 0; j < ns; j++ )
               {    all_reads[j] = MastervecFileObjectCount(
                          dir[j] + "/frag_reads_orig.fastb" );
                    using_reads[j] = int64_t(round( sfs[j] * all_reads[j] ));
                    using_reads[j] = 2 * ( using_reads[j] / 2 );    }
               for ( int j = 1; j < ns; j++ )
                    subsam_starts[j] = subsam_starts[j-1] + using_reads[j-1];
               cout << Date( ) << ": building joint qualp" << endl;
               {    VecPQVec& Q = quals.create( );
                    for ( int j = 0; j < ns; j++ )
                    {    Q.ReadRange( dir[j] + "/frag_reads_orig.qualp",
                               0, using_reads[j] );    }
                    cout << Date( ) << ": writing joint qualp" << endl;
                    quals.store( );    }
               cout << Date( ) << ": building joint fastb" << endl;
               {    vecbasevector& bases = *pReads;
                    for ( int j = 0; j < ns; j++ )
                    {    bases.ReadRange( dir[j] + "/frag_reads_orig.fastb",
                               0, using_reads[j] );    }
                    cout << Date( ) << ": writing joint fastb" << endl;
                    bases.WriteAll( tmp_dir1 + "/frag_reads_orig.fastb" );    }
               SELECT_FRAC = "1.0";    }
          else
          {    String dir;
               if ( IsDirectory( "/wga/scr4/wg_projects/H.sapien/" + sample ) )
                    dir = "/wga/scr4/wg_projects/H.sapien/" + sample;
               if ( sample == "rhino" )
                    dir = "/wga/scr4/wg_projects/C.simum/H77KJADXX";
               if ( sample == "aardvark" )
                    dir = "/wga/scr4/wg_projects/O.afer/H7CRNADXX";
               if ( sample == "lolium" ) dir = "/wga/scr4/wg_projects/L.perenne";
               SystemSucceed( "ln -s " + dir + "/frag_reads_orig.fastb "
                     + tmp_dir1 + "/frag_reads_orig.fastb" );
               pReads->ReadAll(tmp_dir1 + "/frag_reads_orig.fastb");
               SystemSucceed( "ln -s " + dir + "/frag_reads_orig.qualp "
                     + tmp_dir1 + "/frag_reads_orig.qualp" );
               quals.load();    }

          // Report stats.

          int64_t all_reads =
               MastervecFileObjectCount( tmp_dir1 + "/frag_reads_orig.fastb" );
          int64_t using_reads = all_reads;
          if ( reads == "" && SELECT_FRAC.Double( ) < 1.0 )
          {    using_reads = int64_t(round(SELECT_FRAC.Double()*all_reads));
               using_reads = 2*(using_reads/2);    }
          cout << Date( ) << ": using " << ToStringAddCommas(using_reads)
               << " reads" << endl;

          cout << TimeSince(lclock) << " used extracting reads" << endl;
          return;   }

     // Regional processing.

     vecbasevector genome( regions.size( ) );
     vec<String> bams;
     if ( reads != "" ) ParseStringSet( "{" + reads + "}", bams );
     else
     {    String base 
               = "/wga/scr4/picard/H01UJADXX/C1-508_2012-11-01_2012-11-04";
          bams.push_back( base + "/1/Solexa-125532/"
               + "H01UJADXX.1.aligned.duplicates_marked.bam" );
          bams.push_back( base + "/2/Solexa-125532/"
               + "H01UJADXX.2.aligned.duplicates_marked.bam" );    }
     #pragma omp parallel for
     for ( int i = 0; i < regions.isize( ); i++ )
     {    vecbasevector bases;
          vecqualvector quals;
          String reg = regions[i];
          if ( reg.Contains( ":" ) ) 
          {    int start = reg.Between( ":", "-" ).Int( );
               int stop = reg.After( "-" ).Int( );
               reg = reg.Before( ":" ) 
                    + ":" + ToString(start) + "-" + ToString(stop);    }
          String diri = work_dir + "/data/" + ToString(i);
          Mkdir777(diri);
          for ( int z = 0; z < bams.isize( ); z++ )
          {    String getsam = "samtools view -h " + bams[z] + " " + reg;
               SamIAm( z, getsam, diri, false, "", USE_PF_ONLY );
               PairsManager pairs( diri + "/" + ToString(z) + ".pairs" );
               int64_t N = bases.size( ), n;
               {    vecqualvector qualsi( diri + "/" + ToString(z) + ".qualb" );
                    n = pairs.nPairs( ) * 2;
                    quals.resize( N + n );
                    for ( int64_t pid = 0; pid < (int64_t) pairs.nPairs( ); pid++ )
                    {    quals[N + 2*pid] = qualsi[ pairs.ID1(pid) ];
                         quals[N + 2*pid+1] 
                              = qualsi[ pairs.ID2(pid) ];    }    }
               vecbasevector basesi( diri + "/" + ToString(z) + ".fastb" );
               bases.resize( N + n );
               for ( int64_t pid = 0; pid < (int64_t) pairs.nPairs( ); pid++ )
               {    bases[N + 2*pid] = basesi[ pairs.ID1(pid) ];
                    bases[N + 2*pid+1] = basesi[ pairs.ID2(pid) ];     }    }
          bases.WriteAll( diri + "/frag_reads_orig.fastb" );
          quals.WriteAll( diri + "/frag_reads_orig.qualb" );
          int id;
          String sid;                         
          if ( regions[i].Contains( ":" ) ) sid = regions[i].Before(":");
          else sid = regions[i];
          if ( sid == "X" ) id = 22;
          else if ( sid == "Y" ) id = 23;
          else id = sid.Int( ) - 1;
          vecbasevector genomei;
          String genome_fastb = "/wga/dev/references/Homo_sapiens/genome.fastb";
          genomei.ReadOne( genome_fastb, id );
          if ( regions[i].Contains( ":" ) ) 
          {    int start = regions[i].Between( ":", "-" ).Int( );
               int stop = regions[i].After( "-" ).Int( );
               genomei[0].SetToSubOf( genomei[0], start, stop - start );    }
          genome[i] = genomei[0];    }
     genome.WriteAll( work_dir + "/genome.fastb" );
     int64_t genome_size = 0;
     for ( int g = 0; g < (int) genome.size( ); g++ )
          genome_size += genome[g].size( );
     cout << "genome has size " << ToStringAddCommas(genome_size) << endl;
     int64_t nbases = 0;
     for ( int i = 0; i < regions.isize( ); i++ )
     {    nbases += MastervecFileObjectCount(
               tmp_dir1 + "/" + ToString(i) + "/frag_reads_orig.fastb" );    }
     int SEED = 666;
     vec<uint64_t> shuffled;
     Shuffle64( nbases/2, shuffled, (uint64_t) SEED );
     int64_t nusing = nbases;
     if ( SELECT_FRAC.Double( ) < 1.0 ) 
     {    nusing = int( round( SELECT_FRAC.Double( ) * nbases ) );
          if ( nusing % 2 == 1 ) nusing++;    }
     {    vecqualvector qqq;
          for ( int i = 0; i < regions.isize( ); i++ )
          {    qqq.ReadAll( tmp_dir1 + "/" + ToString(i)
                    + "/frag_reads_orig.qualb", True );    }
          vecqualvector quals2(nbases);
          #pragma omp parallel for
          for ( int64_t pid = 0; pid < nbases/2; pid++ )
          {    quals2[2*pid] = qqq[ 2*shuffled[pid] ];
               quals2[2*pid+1] = qqq[ 2*shuffled[pid]+1 ];    }
          convertAssignParallel(quals2.begin(),quals2.begin(nusing),
                                  quals.create());
          quals.store();    }
     {    vecbasevector bases;
          for ( int i = 0; i < regions.isize( ); i++ )
          {    bases.ReadAll( tmp_dir1 + "/" + ToString(i) 
                    + "/frag_reads_orig.fastb", True );    }
          cout << Date( ) << ": using " << nusing << " reads" << endl;
          vecbasevector& bases2 = *pReads;
          bases2.clear().resize(nbases);
          #pragma omp parallel for
          for ( int64_t pid = 0; pid < nbases/2; pid++ )
          {    bases2[2*pid] = bases[ 2*shuffled[pid] ];
               bases2[2*pid+1] = bases[ 2*shuffled[pid]+1 ];    }
          bases2.resize(nusing);
          bases2.WriteAll( tmp_dir1 + "/frag_reads_orig.fastb" );    }
     cout << TimeSince(lclock) << " used extracting reads" << endl;    }

void GetAmbInt( const vecbitvector& amb, vec< pair<int,ho_interval> >& amb_int )
{    for ( int g = 0; g < (int) amb.size( ); g++ )
     {    for ( int i = 0; i < (int) amb[g].size( ); i++ )
          {    if ( !amb[g][i] ) continue;
               int j;
               for ( j = i + 1; j < (int) amb[g].size( ); j++ )
                    if ( !amb[g][j] ) break;
               amb_int.push( g, ho_interval( i, j ) );
               i = j - 1;    }    }    }

void GetCannedReferenceSequences( const String& sample, const String& species,
     const String& work_dir )
{
     // Get reference sequence.

     String genome_fasta, genome_fastb, genome_amb, genome_names;
     String genome_fasta_alt, genome_fastb_alt, genome_amb_alt;
     String genome_names_alt;
     String href_dir = "/wga/scr4/bigrefs";
     if ( species == "human" ) 
     {    genome_fastb = href_dir + "/grch38/genome.fastb";
          genome_names = href_dir + "/grch38/genome.names";
          genome_amb = href_dir + "/grch38/genome.fastamb";
          genome_fastb_alt = href_dir + "/human19/genome.fastb";
          genome_names_alt = href_dir + "/human19/genome.names";
          genome_amb_alt = href_dir + "/human19/genome.fastamb";    }
     if ( species == "mouse" )
     {    genome_fastb = href_dir + "/mouse/genome.fastb";
          genome_names = href_dir + "/mouse/genome.names";
          genome_amb = href_dir + "/mouse/genome.fastamb";    }
     String dref = "/wga/dev/references";
     if ( sample == "rhino" )
          genome_fastb = "/wga/scr4/wg_projects/C.simum/genome.fastb";
     if ( sample == "aardvark" )
          genome_fastb = "/wga/scr4/wg_projects/O.afer/genome.fastb";
     if ( sample == "ecoli12" )
          genome_fastb = dref + "/Escherichia_coli/genome.fastb";
     if ( sample == "rhody" )
          genome_fastb = dref + "/Rhodobacter_sphaeroides/genome.fastb";
     if ( sample == "plasmo" )
     {    genome_fastb = dref + "/Plasmodium_falciparum/genome.fastb";
          genome_names = dref + "/Plasmodium_falciparum/genome.names";    }
     if ( sample == "tb148" )
     {    genome_fasta = "/seq/references/Mycobacterium_tuberculosis_W148/v0/"
               "Mycobacterium_tuberculosis_W148.fasta";    }
     if ( sample == "tbHaarlem" )
          genome_fasta = dref + "/Mycobacterium_tuberculosis/genome.fasta";

     // Build genome files.

     if ( genome_fastb != "" || genome_fasta != "" )
     {    vecbasevector genome;
          if ( genome_fastb != "" )
          {    Cp2( genome_fastb, work_dir );
               genome.ReadAll( work_dir + "/genome.fastb" );    }
          else
          {    FetchReads( genome, 0, genome_fasta );
               genome.WriteAll( work_dir + "/genome.fastb" );    }
          if ( genome_names != "" ) Cp2( genome_names, work_dir );
          else
          {    Ofstream( nout, work_dir + "/genome.names" );
               for ( int g = 0; g < (int) genome.size( ); g++ )
                    nout << g << "\n";    }
          if ( genome_amb != "" ) 
          {    Cp2( genome_amb, work_dir );
               vec< pair<int,ho_interval> > amb_int;
               vecbitvector amb(genome_amb);
               GetAmbInt( amb, amb_int );
               BinaryWriter::writeFile( work_dir + "/genome.ambint", amb_int );    }
          int64_t genome_size = 0;
          for ( int g = 0; g < (int) genome.size( ); g++ )
               genome_size += genome[g].size( );
          cout << "genome has size " << ToStringAddCommas(genome_size) << endl;    }

     // Build genome alt files.

     if ( genome_fastb_alt != "" || genome_fasta_alt != "" )
     {    vecbasevector genome;
          if ( genome_fastb_alt != "" )
          {    Cp2( genome_fastb_alt, work_dir + "/genome.fastb_alt" );
               genome.ReadAll( work_dir + "/genome.fastb_alt" );    }
          else
          {    FetchReads( genome, 0, genome_fasta_alt );
               genome.WriteAll( work_dir + "/genome.fastb_alt" );    }
          if ( genome_names_alt != "" ) 
               Cp2( genome_names_alt, work_dir + "/genome.names_alt" );
          else
          {    Ofstream( nout, work_dir + "/genome.names_alt" );
               for ( int g = 0; g < (int) genome.size( ); g++ )
                    nout << g << "\n";    }
          if ( genome_amb_alt != "" ) 
          {    Cp2( genome_amb_alt, work_dir + "/genome.fastamb_alt" );
               vec< pair<int,ho_interval> > amb_int;
               vecbitvector amb(genome_amb_alt);
               GetAmbInt( amb, amb_int );
               BinaryWriter::writeFile( 
                    work_dir + "/genome.ambint_alt", amb_int );    }
          int64_t genome_size = 0;
          for ( int g = 0; g < (int) genome.size( ); g++ )
               genome_size += genome[g].size( );
          cout << "genome alt has size " << ToStringAddCommas(genome_size) 
               << endl;    }    }
