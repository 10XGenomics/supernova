///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <sys/time.h>

#include "Basevector.h"
#include "CoreTools.h"
#include "FastIfstream.h"
#include "FetchReads.h"
#include "ParallelVecUtilities.h"
#include "ParseSet.h"
#include "TokenizeString.h"
#include "math/Functions.h"
#include "paths/HyperBasevector.h"
#include "paths/LongReadTools.h"
#include "paths/long/CreateGenome.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/fosmid/FosmidPool.h"
#include "paths/long/fosmid/Fosmids.h"
#include "random/Random.h"
#include "system/Worklist.h"
#include "paths/long/RefTrace.h"
#include "paths/long/RefTraceControl.h"

void ParseRegions( const String& X, vec<String>& regions )
{    regions.clear( );
     vec<int> tigs;
     int fail = 0;
     ParseIntSet( X, tigs, fail );
     if ( fail == 0 )
     {    for ( int j = 0; j < tigs.isize( ); j++ )
               regions.push_back( ToString( tigs[j] ) );    }
     else if ( !X.Contains( "random", 0 ) ) ParseStringSet( X, regions );    }

void GetRegion( const String& SAMPLE, const String& region, 
     const vecbasevector& genome, const vecbasevector& genomeg, vecbasevector& G, 
     vecbasevector& Gplus, vec<int>& Gplus_ext, vec<HyperBasevector>& GH, 
     vec<double>& ploidy, const int ext, int& found )
{    
     if ( SAMPLE == "human" )
     {    String rdir = "/wga/dev/references/Homo_sapiens/NA12878_regions";
          vec<String> knowns = AllFiles(rdir);
          int p = BinPosition( knowns, region + ".vhbv" );
          if ( p >= 0 )
          {    vec<HyperBasevector> x;
               BinaryReader::readFile( rdir + "/" + region + ".vhbv", &x );
               ForceAssertEq( x.isize( ), 2 );
               vec<int> sources, sinks;
               for ( int m = 0; m < x.isize( ); m++ )
               {    x[m].Sources(sources); 
                    x[m].Sinks(sinks);
                    if ( !sources.solo() || !sinks.solo() )
                    {    FatalErr("\nEach reference sequence must have a "
                                   "unique source and a unique sink.\n"
                                   "Abort.");    }
                    vec< vec<int> > paths;
                    // note could upgrade to faster version that uses to_left,right
                    x[m].EdgePaths( sources[0], sinks[0], paths );
                    basevector z = x[m].EdgeObject( paths[0][0] );
                    for ( int i = 1; i < paths[0].isize( ); i++ )
                         z = TrimCat( x[m].K( ), z, x[m].EdgeObject( paths[0][i] ) );
                    G.push_back(z);
                    Gplus.push_back(z);
                    Gplus_ext.push_back(0);
                    ploidy.push_back(0.5);    }
               GH.append(x);
               found++;
               return;    }    }
     int GSTART = region.Between( ":", "-" ).Int( );
     int GSTOP = region.After( "-" ).Int( );
     ForceAssert( GSTOP >= 0 );
     ForceAssertLe( GSTART, GSTOP );
     if ( GSTOP > genome[0].isize( ) )
     {    cout << "\n";
          PRINT2( GSTOP, genome[0].size( ) );
          FAIL_MSG( "Region specified by X goes beyond chromosome end." );    }
     basevector g = genome[0];
     g.SetToSubOf( g, GSTART, GSTOP - GSTART );
     G.push_back_reserve(g);
     basevector gplus = genome[0];
     int GSTART_ext = Max( 0, GSTART - ext );
     int GSTOP_ext = Min( genome[0].isize( ), GSTOP + ext );
     gplus.SetToSubOf( gplus, GSTART_ext, GSTOP_ext - GSTART_ext );
     Gplus.push_back_reserve(gplus);    
     Gplus_ext.push_back( GSTART - GSTART_ext );    
     ploidy.push_back(1.0);    
     vec<basevector> x;
     x.push_back(g);
     const int K = 100;
     GH.push( K, x );    }

Bool TestAmb( const String& region, const String& genome_head, const int id, 
     int start, int stop, const Bool extend )
{    if ( IsRegularFile( genome_head + ".fastamb" ) )
     {    const int ext = 50000;
          if (extend)
          {    start -= ext;
               stop += ext;    }
          vecbitvector amb;
          vec<int> amblen;
          amb.ReadOne( genome_head + ".fastamb", id );
          for ( int j = start; j < stop; j++ )
          {    if ( j < 0 || j >= (int) amb[0].size( ) ) continue;
               if ( !amb[0][j] ) continue;
               int k;
               for ( k = j + 1; k < stop; k++ )
                    if ( !amb[0][k] ) break;
               amblen.push_back( k - j );
               j = k - 1;    }
          ReverseSort(amblen);
          const int min_amb_call = 100;
          if ( amblen.nonempty( ) && amblen[0] >= min_amb_call )
          {    cout << Date( ) << ": Warning, ";
               if (extend) cout << "extension of ";
               cout << region << " has " << "seq of ambiguous bases of len "
                    << amblen[0] << "." << endl;
               return True;    }    }
     return False;    }

// note that the SAMPLE="human" IN_GENOME="" path of
// this function is duplicated and modded to
// CreateGenome_Discovar, as a quick fix of memory/speed-related issue.
// Any changes here should also be made there

void CreateGenome( const String& IN_GENOME, const String& SAMPLE, const String& X, 
     const String& HUMAN_CONTROLS, String& X_actual, ref_data& ref,
     const double GENOME_SUB_PERCENT, RefTraceControl& RTCtrl )
{    
     vecbasevector& G = ref.G;
     vecbasevector& Gplus = ref.G3plus;
     vec<int>& Gplus_ext = ref.Gplus_ext;
     vec<HyperBasevector>& GH = ref.GH;
     vec<bool>& is_circular = ref.is_circular;
     vec<double>& ploidy = ref.ploidy;
     if ( IN_GENOME != "" )
     {    if ( IN_GENOME.Contains( ".fastb", -1 ) ) G.ReadAll(IN_GENOME);
          if ( IN_GENOME.Contains( ".fasta", -1 ) ) 
          {    FetchReads(G, 0, IN_GENOME);
               vec<basevector> p;
               for ( int i = 0; i < (int) G.size( ); i++ )
                    p.push_back( G[i] );
               const int K = 80; // no justification for this
               GH.push_back( HyperBasevector( K, p ) );    }
          /* Note not setting GH. */
          if ( IN_GENOME.Contains( ".hbv", -1 ) ) 
          {    HyperBasevector hb;
               BinaryReader::readFile( IN_GENOME, &hb );
               for ( int e = 0; e < hb.EdgeObjectCount( ); e++ )
                    G.push_back( hb.EdgeObject(e) );
               GH.push_back(hb);    }    
          Gplus = G;
          Gplus_ext = vec<int>( G.size( ), 0 );
          is_circular.resize( G.size( ), false );
          ploidy.resize( G.size( ), 1 );
          return;    }
     
     if ( !( SAMPLE == "plasmo" || SAMPLE == "ecoli" || SAMPLE == "ecoli11"
          || SAMPLE == "bcereus" || SAMPLE == "entero" || SAMPLE == "tb"
          || SAMPLE == "bifi" || SAMPLE == "scardovia" || SAMPLE == "tb148"
          || SAMPLE == "rhody" || SAMPLE.Contains( "random:", 0 ) 
          || SAMPLE == "mouse" || SAMPLE == "human" || SAMPLE == "hpool1" 
          || SAMPLE == "hpool2" || SAMPLE == "hpool3" || SAMPLE == "human.hpool2"
          || SAMPLE == "human.hpool3" || SAMPLE == "neiss" || SAMPLE == "arabidopsis"
          || SAMPLE == "ecoli_scs" || SAMPLE == "rhino"
          || SAMPLE == "ecoli12" || SAMPLE.Contains( ".fastb", -1 ) ) )
     {    FAIL_MSG( "Illegal specification of SAMPLE." );    }

     // Parse regions.

     vec<String> regions;
     ParseRegions( X, regions );
     const int ext = 500;

     // If human controls are to be used, set the reference sequence appropriately.

     if ( HUMAN_CONTROLS != "" )
     {    String dir 
               = "/wga/dev/references/Homo_sapiens/NA12878_Fosmid_Pool.regions.fin";
          vec<String> all = GetFinishedFosmidFiles( );
          vec<int> hc;
          if ( HUMAN_CONTROLS != "all" )
          {    if ( isdigit( HUMAN_CONTROLS[0] ) ) ParseIntSet( HUMAN_CONTROLS, hc );
               else hc = expand_fosmids(HUMAN_CONTROLS);    }
          int hcount = 0;
          for ( const String& f : all )
          {    int id = f.Between( ".", "." ).Int( );
               if ( HUMAN_CONTROLS != "all" && !BinMember( hc, id ) ) continue;
               hcount++;
               vecbasevector b;
               FetchReads( b, 0, dir + "/" + f );
               G.Append(b);    }
          if ( HUMAN_CONTROLS != "all" && hcount != hc.isize( ) )
          {    FatalErr( "It would appear that you may have specified indices "
                    "of Fosmids for which there is no reference sequence.\n"
                    << "Abort." );    }
          Gplus = G;
          Gplus_ext.resize( G.size( ), 0 );
          GH.resize( G.size( ) );
          for ( int g = 0; g < (int) G.size( ); g++ )
          {    vec<basevector> x;
               x.push_back( G[g] );
               const int K = 100;
               GH[g] = HyperBasevector( K, x );    }
          is_circular.resize( G.size( ), False );
          ploidy.resize( G.size( ), 1.0 );
          double true_size = 0;
          for ( int g = 0; g < (int) G.size( ); g++ )
               true_size += G[g].size( ) * ploidy[g];
          cout << Date( ) << ": genome/region has size " 
               << ToStringAddCommas( int64_t(round(true_size)) ) << endl;
          return;    }

     // Handle the case of a random genome.

     if ( SAMPLE.Contains( "random:", 0 ) )
     {    vec<String> chrs;
          vec<char> comma, vert;
          comma.push_back( ',' ), vert.push_back( '|' );
          TokenizeStrictly( SAMPLE.After( "random:" ), comma, chrs );
          map<String,basevector> defs;
          for ( int g = 0; g < chrs.isize( ); g++ )
          {    vec<String> x;
               TokenizeStrictly( chrs[g], vert, x );
               basevector B;
               bool circ = false;
               for ( int j = 0; j < x.isize( ); j++ )
               {    String c = x[j];
                    if ( c == "" && j == x.isize( ) - 1 ) circ = true;
                    else if ( c.Contains( "(" ) && c.After( "(" ).Contains( ")" ) )
                    {    String head = c.Before( "(" ), tail = c.After( ")" );
                         String sz = c.Between( "(", ")" );
                         int64_t mult = 1;
                         if ( sz.Contains( "K", -1 ) )
                         {    mult = 1000;
                              sz = sz.RevBefore( "K" );    }
                         if ( defs.find(head) != defs.end( ) || head.Contains( "'" )
                              || tail != "" || !sz.IsInt( ) )
                         {    FAIL_MSG( "Illegal random genome spec." );    }
                         int64_t s = sz.Int( ) * mult;
                         basevector b(s);
                         for ( int64_t i = 0; i < s; i++ )
                              b.Set( i, randomx( ) % 4 );
                         defs[head] = b;
                         B = Cat( B, b );    }
                    else
                    {    Bool rc = False;
                         if ( c.Contains( "'", -1 ) )
                         {    c = c.RevBefore( "'" );
                              rc = True;    }
                         if ( defs.find(c) == defs.end( ) )
                              FAIL_MSG( "Illegal random genome spec." );
                         basevector b = defs[c];
                         if (rc) b.ReverseComplement( );
                         B = Cat( B, b );    }    }
               G.push_back_reserve(B);
               ploidy.push_back(1.0);
               is_circular.push_back(circ);    }
          Gplus = G;
          Gplus_ext.resize( G.size( ), 0 );    }

     // Handle the case where the genome is nonrandom.

     String genome_head = "/wga/dev/references/";
     vecbasevector genome;
     if ( SAMPLE == "ecoli" ) genome_head += "Escherichia_coli/genome";
     if ( SAMPLE == "ecoli11" ) genome_head = "/seq/references/"
          "Escherichia_coli_EC11-9941/v0/Escherichia_coli_EC11-9941";
     if ( SAMPLE == "ecoli12" ) 
          genome_head = "/wga/dev/references/Escherichia_coli/genome";
     if ( SAMPLE == "ecoli_scs" )
          genome_head = "/wga/dev/references/Escherichia_coli/SCS_110/genome";
     if ( SAMPLE == "bcereus" ) genome_head = "/seq/references/"
          "Bacillus_cereus_VD022/v0/Bacillus_cereus_VD022";
     if ( SAMPLE == "rhino" ) 
     {    genome_head = "/seq/references/Ceratotherium_simum_simum/v0/"
               "Ceratotherium_simum_simum";    }
     if ( SAMPLE == "entero" ) 
          genome_head = "/wga/dev/references/Enterococcus_casseliflavus/genome";
     if ( SAMPLE == "tb148" ) genome_head = "/seq/references/"
               "Mycobacterium_tuberculosis_W148/v0/Mycobacterium_tuberculosis_W148";
     if ( SAMPLE == "bifi" )
          genome_head = "/wga/dev/references/Bifidobacterium_bifidum/genome";
     if ( SAMPLE == "scardovia" )
          genome_head = "/wga/dev/references/Scardovia_wiggsiae/genome";
     if ( SAMPLE == "tb" )
          genome_head = "/wga/dev/references/Mycobacterium_tuberculosis/genome";
     if ( SAMPLE == "rhody" ) genome_head += "Rhodobacter_sphaeroides/genome";
     if ( SAMPLE == "plasmo" ) genome_head += "Plasmodium_falciparum/genome";
     if ( SAMPLE == "neiss" ) genome_head += "Neisseria_meningitidis/genome";
     if ( SAMPLE == "mouse" ) genome_head += "Mus_musculus/genome_fixed_break1000";
     if ( SAMPLE == "human" || SAMPLE == "human.hpool2" || SAMPLE == "human.hpool3" )
          genome_head += "Homo_sapiens/genome";
     if ( SAMPLE == "hpool1" ) genome_head += "Homo_sapiens/WIBR_Fosmid_Pool";
     if ( SAMPLE == "hpool2" || SAMPLE == "hpool3" ) 
     {    // genome_head += "Homo_sapiens/NA12878_Fosmid_Pool";
          genome_head += "Homo_sapiens/genome";    }
     if ( SAMPLE == "arabidopsis" ) genome_head = "/seq/references/"
          "Arabidopsis_thaliana_TAIR10/v1/Arabidopsis_thaliana_TAIR10";
     if ( SAMPLE.Contains( ".fastb", -1 ) ) genome_head = SAMPLE.RevBefore(".fastb");
     if ( !SAMPLE.Contains( "random:", 0 ) )
     {    String genome_fastb = genome_head + ".fastb";
          if ( X.Contains( "random", 0 ) )
          {    int n = X.Between( "random(", ")" ).Int( );
               if ( IsRegularFile(genome_fastb) ) genome.ReadAll(genome_fastb);
               else FetchReads( genome, 0, genome_head + ".fasta" );
               Bool possible = False;
               for ( int g = 0; g < (int) genome.size( ); g++ )
                    if ( genome[g].isize( ) >= n ) possible = True;
               if ( !possible )
                    FAIL_MSG( "Genome doesn't have a region of size " << n << "." );
               timeval t;
               gettimeofday( &t, NULL );
               unsigned int seed = t.tv_sec + t.tv_usec;
               srandomx(seed);
               int64_t N = genome.SizeSum( );
               G.resize(1);
               Gplus.resize(1);
               Gplus_ext.resize(1);
               ploidy.push_back(1.0);

               RTCtrl.getRefHead() = genome_head;
               RTCtrl.getRefTags().resize(1);
               RTCtrl.getRefStart().resize(1);
               RTCtrl.getRefSeqs().resize(1);
               while(1)
               {    int64_t start = big_random( ) % N;
                    for ( int g = 0; g < (int) genome.size( ); g++ )
                    {    if ( start + n <= genome[g].isize( ) )
                         {    int64_t stop = start + n;
                              G[0].SetToSubOf( genome[g], start, stop - start );

                              if(SAMPLE=="human"){
                                  switch(g){
                                      case 22: RTCtrl.getRefTags()[0]="X"; break;
                                      case 23: RTCtrl.getRefTags()[0]="Y"; break;
                                      default: RTCtrl.getRefTags()[0]=ToString(g+1);break; } }
                              else{ RTCtrl.getRefTags()[0]=ToString(g); }
                              RTCtrl.getRefIndex()[0]= g;
                              RTCtrl.getRefStart()[0]=start;
                              RTCtrl.getRefSeqs()[0]=fastavector(G[0]);

                              int64_t start_ext = Max( (int64_t) 0, start - ext );
                              int64_t stop_ext 
                                   = Min( (int64_t) genome[g].size( ), stop + ext );
                              Gplus[0].SetToSubOf( 
                                   genome[g], start_ext, stop_ext - start_ext );
                              Gplus_ext[0] = start - start_ext;
                              X_actual = ToString(g) + ":" + ToString(start)
                                   + "-" + ToString(start+n);
                              cout << Date( ) << ": using X=" << g << ":"
                                   << start << "-" << start+n << endl;
                              break;    }
                         start -= genome[g].isize( );
                         if ( start < 0 ) break;    }
                    if ( G[0].size( ) > 0 ) break;    }    }
          else
          {    
               // Translate regions in hpool2 and hpool3 cases.

               if ( SAMPLE == "hpool2" || SAMPLE == "hpool3" )
               {    vec<String> regions0, regionsx, regions2;
                    vec< vec< pair<String,String> > > junctions, breaks, edits;
                    ParseFosmidPoolMetainfo( regions0, junctions, breaks, edits );
                    for ( int i = 0; i < regions0.isize( ); i++ )
                    {    String id = regions0[i].Before( ":" );
                         int nid = ( id == "X" ? 23 : id.Int( ) ) - 1;
                         String range = regions0[i].After( ":" );
                         if ( range.Contains( " " ) ) range = range.Before( " " );
                         regionsx.push_back( ToString(nid) + ":" + range );    }
                    if ( regions.empty( ) ) 
                    {    if ( SAMPLE == "hpool2" )
                         {    for ( int id = 0; id <= 55; id++ )
                                   regions.push_back( regionsx[id] );    }
                         else
                         {    for ( int id = 56; id <= 106; id++ )
                                   regions.push_back( regionsx[id] );    }    }
                    else
                    {    for ( int i = 0; i < regions.isize( ); i++ )
                         {    ForceAssert( regions[i].IsInt( ) );
                              int id = regions[i].Int( );
                              ForceAssert( id >= 0 && id < regionsx.isize( ) );
                              regions2.push_back( regionsx[id] );    }
                         regions = regions2;    }    }

               // Main case.

               String gdir = "/wga/dev/references/Homo_sapiens";
               if ( regions.empty( ) ) 
               {    if ( IsRegularFile(genome_fastb) ) genome.ReadAll(genome_fastb);
                    else FetchReads( genome, 0, genome_head + ".fasta" );
                    G = genome;
                    for ( int j = 0; j < (int) genome.size( ); j++ )
                         ploidy.push_back( 1.0 );
                    Gplus = genome;
                    Gplus_ext.resize( G.size( ), 0 );    }
               else 
               {    int npasses = 1;
                    int found = 0;
                    for ( int pass = 0; pass < npasses; pass++ )
                    for ( int j = 0; j < regions.isize( ); j++ )
                    {    int id;
                         String sid;                         
                         if ( regions[j].Contains( ":" ) ) 
                              sid = regions[j].Before(":");
                         else sid = regions[j];
                         if ( SAMPLE == "human" )
                         {    if ( sid == "X" ) id = 22;
                              else if ( sid == "Y" ) id = 23;
                              else id = sid.Int( ) - 1;    }
                         else id = sid.Int( );    
                         if ( !regions[j].Contains( ":" ) )
                         {    genome.clear( );
                              if ( IsRegularFile(genome_fastb) ) 
                                   genome.ReadOne( genome_fastb, id );
                              else 
                              {    vecbasevector genomea;
                                   FetchReads( genomea, 0, genome_head + ".fasta" );
                                   genome.resize(1);
                                   genome[0] = genomea[id];    }
                              if ( IsRegularFile(genome_fastb) ) 
                              {    TestAmb( regions[j], genome_head, id, 0, 
                                        genome[0].size( ), False );    }
                              G.push_back_reserve( genome[0] );
                              RTCtrl.getRefHead() = genome_head;
                              RTCtrl.getRefIndex().push_back(id);
                              RTCtrl.getRefTags().push_back(ToString(sid));
                              RTCtrl.getRefStart().push_back(0);
                              RTCtrl.getRefSeqs().push_back(fastavector(G.back()));
                              Gplus.push_back_reserve( genome[0] );
                              Gplus_ext.push_back(0);
                              ploidy.push_back( 1.0 );    }
                         else
                         {    genome.clear( );
                              if ( IsRegularFile(genome_fastb) ) 
                                   genome.ReadOne( genome_fastb, id );
                              else 
                              {    vecbasevector genomea;
                                   FetchReads( genomea, 0, genome_head + ".fasta" );
                                   genome.resize(1);
                                   genome[0] = genomea[id];    }
                              int start = regions[j].Between( ":", "-" ).Int( );
                              int stop = regions[j].After( "-" ).Int( );
                              if ( IsRegularFile(genome_fastb) )
                              {    if ( !TestAmb( regions[j], genome_head, id, 
                                        start, stop, False ) )
                                   {    TestAmb( regions[j], genome_head, id, 
                                             start, stop, True );    }    }
                              vecbasevector genomeg;
                              GetRegion( SAMPLE, regions[j], genome, genomeg, G, 
                                   Gplus, Gplus_ext, GH, ploidy, ext, found );    
                              RTCtrl.getRefHead() = genome_head;
                              RTCtrl.getRefIndex().push_back(id);
                              RTCtrl.getRefTags().push_back(ToString(sid));
                              RTCtrl.getRefStart().push_back(start);
                              RTCtrl.getRefSeqs().push_back(fastavector(G.back()));
                                   }    }
                    if ( found > 0 )
                    {    cout << Date( ) << ": found known assemblies for " << found 
                              << " of " << regions.size( ) << " regions" 
                              << endl;    }    }    }    }
     double true_size = 0;
     for ( int g = 0; g < (int) G.size( ); g++ )
          true_size += G[g].size( ) * ploidy[g];
     cout << Date( ) << ": genome/region has size " 
          << ToStringAddCommas( int64_t(round(true_size)) ) << endl;
     if ( !SAMPLE.Contains( "random:", 0 ) ) is_circular.resize( G.size( ), False );
     if ( SAMPLE == "rhody" || SAMPLE == "neiss" || SAMPLE == "ecoli11" 
          || SAMPLE == "bcereus" || SAMPLE == "tb" || SAMPLE == "entero"
          || SAMPLE == "bifi" || SAMPLE == "scardovia" || SAMPLE == "tb148"
          || SAMPLE == "ecoli12" || SAMPLE == "ecoli_scs" )
     {    for ( int i = 0; i < (int) G.size( ); i++ )
               if ( X == "" || regions[i].IsInt( ) ) is_circular[i] = True;    }

     // Create diploid genome.

     if ( GENOME_SUB_PERCENT > 0.0 )
     {    vecbasevector G2(G), G2plus(Gplus);
          G2.Append(G);
          is_circular.append(is_circular);
          RTCtrl.getRefIndex().append(RTCtrl.getRefIndex());
          RTCtrl.getRefTags().append(RTCtrl.getRefTags());
          RTCtrl.getRefStart().append(RTCtrl.getRefStart());
          RTCtrl.getRefSeqs().append(RTCtrl.getRefSeqs());
          vec<double> ploidy2;
          for ( int j = 0; j < ploidy.isize( ); j++ )
               ploidy2.push_back( ploidy[j]/2, ploidy[j]/2 );
          ploidy = ploidy2;
          srandomx(1234567);
          for ( int g = (int) G.size( ); g < (int) G2.size( ); g++ )
          {    basevector& b = G2[g];
               for ( int j = 0; j < b.isize( ); j++ )
               {    if ( randomx( ) % 1000000 < GENOME_SUB_PERCENT * 10000.0 )
                    {    int add = 1 + ( randomx( ) % 3 );
                         b.Set( j, ( b[j] + add ) % 4 );    }    }    }
          G = G2;    
          Gplus = G2;    
          Gplus_ext.append(Gplus_ext);    }

     // Create HyperBasevector reference.

     if ( GH.empty( ) )
     {    GH.resize( G.size( ) );
          for ( int g = 0; g < (int) G.size( ); g++ )
          {    vec<basevector> x;
               x.push_back( G[g] );
               const int K = 100;
               GH[g] = HyperBasevector( K, x );    }    }    }


// note that this function is duplicated and modded from the SAMPLE="human" IN_GENOME="" path of
// CreateGenome, as a quick fix of memory/speed-related issue.
// Any changes there should also be made here
void CreateGenome_Discovar( ref_data& ref, const DiscovarTools::DiscovarRefTraceControl& vControl )
{
     vecbasevector& G = ref.G;
     vecbasevector& Gplus = ref.G3plus;
     vec<int>& Gplus_ext = ref.Gplus_ext;
     vec<HyperBasevector>& GH = ref.GH;
     vec<bool>& is_circular = ref.is_circular;
     vec<double>& ploidy = ref.ploidy;
     //this is mimic the SAMPLE="human" IN_GENOME="" execution path of CreateGenome called by LongProto
     {
         const size_t nSequence=vControl.getRefSeqs().size();
         srandomx(1643295582);
         for(size_t ii=0;ii<nSequence;++ii){
             Gplus.push_back_reserve(vControl.getRefSeqsExt()[ii].ToBasevectorRandom());

             basevector tmp( Gplus.back(), vControl.getRefSeqsExtExt()[ii], vControl.getRefSeqs()[ii].size());

             G.push_back_reserve(tmp);

             Gplus_ext.push_back( vControl.getRefSeqsExtExt()[ii] );
             ploidy.push_back(1.0);
             vec<basevector> x;
             x.push_back(tmp);
             const int K = 100;
             GH.push( K, x );

         }
     }
     double true_size = 0;
     for ( int g = 0; g < (int) G.size( ); g++ )
          true_size += G[g].size( ) * ploidy[g];
     cout << Date( ) << ": genome/region has size "
          << ToStringAddCommas( int64_t(round(true_size)) ) << endl;
     is_circular.resize( G.size( ), False );

     // Create HyperBasevector reference.

     if ( GH.empty( ) )
     {    GH.resize( G.size( ) );
          for ( int g = 0; g < (int) G.size( ); g++ )
          {    vec<basevector> x;
               x.push_back( G[g] );
               const int K = 100;
               GH[g] = HyperBasevector( K, x );    }    }    }


namespace
{

struct WorkItem
{
    WorkItem() : mRefId(0), mStart(0), mEnd(0), mDest(0) {}
    WorkItem( size_t refId, unsigned start, unsigned end, size_t dest )
    : mRefId(refId), mStart(start), mEnd(end), mDest(dest) {}

    // default copying, moving, and destructor are OK

    size_t mRefId;
    unsigned mStart;
    unsigned mEnd;
    size_t mDest;
};

typedef triple<int,int,int> IntTriple;

class KmerizationProc
{
public:
    KmerizationProc( vecbvec const& ref, int K, vec<IntTriple>& out )
    : mpRef(&ref), mpOut(&out), mK(K) {}

    // default copying, moving, and destructor are OK

    void operator()( WorkItem const& item ) const
    { bvec const& tig = (*mpRef)[item.mRefId];
      auto oItr = mpOut->begin()+item.mDest;
      for ( unsigned idx = item.mStart; idx != item.mEnd; ++idx )
      { *oItr = IntTriple(KmerId(tig,mK,idx),item.mRefId,idx); ++oItr; } }

private:
    vecbvec const* mpRef;
    vec<triple<int,int,int>>* mpOut;
    int mK;
};


} // end of anonymous namespace

void CreateGlocs( const vecbasevector& G, unsigned const LG, VecIntPairVec& Glocs )
{
    size_t nKmers = 0;
    for ( auto itr=G.begin(),end=G.end(); itr != end; ++itr )
        if ( itr->size() >= LG )
            nKmers += itr->size()-LG+1;
    vec<IntTriple> out(nKmers);
    if ( true )
    {
        KmerizationProc proc(G,LG,out);
        Worklist<WorkItem,KmerizationProc> wl(proc);
        size_t nTigs = G.size();
        size_t dest = 0;
        unsigned const BATCH_SIZE = 100000;
        ForceAssertGe(BATCH_SIZE,LG);
        for ( size_t idx = 0; idx != nTigs; ++idx )
        {
            unsigned off = 0;
            unsigned sz = G[idx].size();
            if ( sz < LG )
                continue;
            sz = sz - LG + 1;
            while ( sz >= BATCH_SIZE )
            {
                wl.add(WorkItem(idx,off,off+BATCH_SIZE,dest));
                off += BATCH_SIZE;
                dest += BATCH_SIZE;
                sz -= BATCH_SIZE;
            }
            if ( sz )
            {
                wl.add(WorkItem(idx,off,off+sz,dest));
                dest += sz;
            }
        }
    }
    ParallelSort(out);
    Glocs.clear().resize( 1 << (2*LG) );
    auto beg = out.begin(), end = out.end();
    while ( beg != end )
    {
        auto itr = beg;
        int kmerId = beg->first;
        while ( ++itr != end )
            if ( kmerId != itr->first )
                break;
        IntPairVec& locs = Glocs[kmerId];
        locs.reserve(itr-beg);
        while ( beg != itr )
        {
            locs.push_back(std::pair<int,int>(beg->second,beg->third));
            ++beg;
        }
    }
}
