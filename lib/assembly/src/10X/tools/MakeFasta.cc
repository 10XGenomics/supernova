// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

// MakeFasta.  Generate Supernova output

#include "MainTools.h"
#include "graph/Digraph.h"
#include "graph/DigraphTemplate.h"
#include "10X/astats/AssemblyStats.h"
#include "paths/long/large/Lines.h"
#include "VecString.h"
#include "10X/writestuff/ScafLinePrinter.h"
#include "10X/SuperFiles.h"
#include "system/System.h"


const String base_version="1.7";

void TXT2GZ(const String& fn){
     ifstream in;
     std::string data_fn;
     in.open(fn);
     if(in){
         if (!in.seekg(0, std::ios::end)){
             cout<<"Error seeking in file " << fn <<endl;
             Scram(1);
         }
         data_fn.reserve(in.tellg());
         if(!in.seekg(0, std::ios::beg)){
             cout<<"Error seeking in file " << fn <<endl;
             Scram(1);
         }
         data_fn.assign((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
         in.close();
     }else{
         cout<<"Error opening " << fn << " for reading"<<endl;
         Scram(1);
     }
     FstreamWriteGZ(fn, data_fn.c_str());
     cout<<"writing gzipped object"<<endl;
}

int main( int argc, char *argv[] )
{
     RunTime( );

     String sflavors = "{raw,megabubbles,pseudohap,pseudohap2}";
     vec<String> vflavors;
     ParseStringSet(sflavors, vflavors);

     BeginCommandArguments;
     CommandArgument_String_OrDefault_Doc(DIR, ".", "/final assembly directory");
     CommandArgument_String_Doc(OUT_HEAD, "base filename for the output files" );
     CommandArgument_Bool_OrDefault_Doc(ABBREV, True, "abbreviate edge ids" );
     CommandArgument_Bool_OrDefault_Doc(ZIP, True, "output in gzip format" );
     CommandArgument_UnsignedInt_OrDefault_Doc(MINSIZE, 1000, "minimum FASTA record size (except for 'raw')" );
     CommandArgument_String_OrDefault_Valid_Doc(FLAVOR, "megabubbles", sflavors,
               "output flavor, please see documentation at http://software.10xgenomics.com/supernova" );

#if !defined(CS)
     // These arguments are hidden from user
     CommandArgument_Bool_OrDefault_Doc(SKIP_GAPS, False, "when GAPS=separate, skip gap records altogether");
     CommandArgument_String_OrDefault_Doc(REINS_DIR, "", "write reinserted graph to a directory");
     CommandArgument_String_OrDefault_Doc(LOG,"","write verbose information to log file");
#else
     Bool SKIP_GAPS=False;
     String REINS_DIR="";
     String LOG="";
#endif
     EndCommandArguments;



     int nflavor = Position( vflavors, FLAVOR );
     if ( nflavor < 0 ) FatalErr( "unknown output flavor: " + FLAVOR );

     String version = base_version + " style=" + ToString(nflavor+1);

     cout << "version=" << version << endl;

     HyperBasevectorX hbd;
     digraphE<vec<int>> D;
     vec<int> dinv;
     LineVec dlines;
     HyperBasevectorX hb;
     vec<int> inv;

     auto start=WallClockTime();

     Bool const SEQ_IDS = True;

     cout << Date() << ": loading assembly files" << endl;
     BinaryReader::readFile( DIR + "/a.sup.hbx", &hbd);
     BinaryReader::readFile( DIR + "/a.sup", &D);
     BinaryReader::readFile( DIR + "/a.sup.inv", &dinv);
     BinaryReader::readFile( DIR + "/a.sup.lines", &dlines );
     BinaryReader::readFile( DIR + "/../a.hbx", &hb);
     BinaryReader::readFile( DIR + "/../a.inv", &inv);

     cout << Date() << ": cell reinsertion" << endl;
     digraphE<vec<int>> D2(D);
     vec<int> dinv2(dinv);
     ReinsertLoops(hb,inv,D2,dinv2);

     HyperBasevectorX hbd2;
     SuperToSeqGraph( hb, D2, hbd2 );




     auto dpaths_fn = DIR + "/a.dpaths.counts";
     vec<uint16_t> dpaths_counts;
     if ( !IsRegularFile(dpaths_fn) ) {
          cout << Date() << ": missing dpaths counts file.. generating" << endl;
          ReadPathVec dpaths;
          VecULongVec dpaths_index;
          cout << Date() << ":\tloading dpaths..." << endl;
          dpaths.ReadAll( DIR+"/a.dpaths" );
          cout << Date() << ":\tloading dpaths index..." << endl;
          dpaths_index.ReadAll( DIR+"/a.dpaths.index" );
          cout << Date() << ":\tcalculating counts" << endl;
          GenerateDpathsCounts(dpaths, dpaths_index, dinv, dpaths_counts);
          BinaryWriter::writeFile( dpaths_fn, dpaths_counts );
          size_t zero_count=0;
          for ( size_t i = 0; i < dpaths_index.size(); ++i )
               if ( dpaths_index[i].size() == 0 ) zero_count++;
          cout << 100.0 * zero_count / dpaths_index.size() <<
               "% of edges in D are not in read paths" << endl;
     } else {
          BinaryReader::readFile( dpaths_fn, &dpaths_counts );
     }

     if ( REINS_DIR != "" ) { 
          cout << Date() << ": writing reinserted graph to: " << REINS_DIR << endl;
          if (!IsDirectory(REINS_DIR)) Mkdir777(REINS_DIR);
          BinaryWriter::writeFile( REINS_DIR + "/a.sup", D2);
          BinaryWriter::writeFile( REINS_DIR + "/a.sup.inv", dinv2);
          BinaryWriter::writeFile( REINS_DIR+ "/a.sup.hbx", hbd2 );
          hbd2.Edges( ).WriteAll( REINS_DIR +"/a.sup.fastb"  );
     }

     if ( FLAVOR == "raw" ) {
          String fasta_fn = OUT_HEAD+".fasta";

          // TODO: gap absorb needs to be pulled out of ScaflinePrinter if it's going to be here...
          FastaEdgeWriter fWriter(fasta_fn, version, hbd2.K(), False /* abbrev */, True /* absorb gaps */, LOG, False /* skip gaps */, SEQ_IDS, 0 /* minsize */);

          vec<Bool> used;
          hbd2.Used(used);

          int iter=0, ndots=0;
          cout << "writing 100 dots in 10 groups of 10:" << endl;
          for ( int i = 0; i < hbd2.E(); ++i ) {
               int vleft = hbd2.ToLeft(i);
               int vright = hbd2.ToRight(i);
               basevector const& obj = hbd2.O(i);
               if ( used[i] ) {
                    if ( obj.size() == 0 )
                         fWriter.AddGapEdge( vleft, vright, i);
                    else
                         fWriter.AddEdge( vleft, vright, i, hbd2.O(i) );
                    fWriter.Break();
               }
               MakeDots(iter,ndots,hbd2.E());
          }

          cout << endl << "completed conversion in " << TimeSince(start) << endl;
          if(ZIP) TXT2GZ(fasta_fn);
          return 0;
     }


     ScafLinePrinter slp( hb, hbd, D, dinv, dlines, dpaths_counts, D2, hbd2 );
     slp.SetKeepRc( False );
     slp.SetSeparateGaps( False );
     slp.SetBreakBubbles( False );

     if ( FLAVOR == "megabubbles" ) {
          String fasta_fn = OUT_HEAD+".fasta";
          slp.SetMashMegaBubbles( False );
          cout << Date() << ": writing the FASTA portion of output to " << fasta_fn << endl;
          FastaEdgeWriter fWriter(fasta_fn, version, hbd.K(), ABBREV, True /* absorb gaps */, LOG, False /* skip gaps */, SEQ_IDS, MINSIZE);
          slp.WalkScaffoldLines( fWriter );
          if(ZIP) TXT2GZ(fasta_fn);
     } else if ( FLAVOR == "pseudohap" ) {
          slp.SetMashMegaBubbles( True );
          String fasta_fn = OUT_HEAD+".fasta";
          cout << Date() << ": writing the FASTA portion of output to " << fasta_fn << endl;
          // TODO: gap absorb needs to be pulled out of ScaflinePrinter if it's going to be here...
          FastaEdgeWriter fWriter(fasta_fn, version, hbd.K(), ABBREV, True /* absorb gaps */, LOG, False /* skip gaps */, SEQ_IDS, MINSIZE);
          slp.WalkScaffoldLines( fWriter );
          if(ZIP) TXT2GZ(fasta_fn);
     } else if ( FLAVOR == "pseudohap2" ) {
          slp.SetMashMegaBubbles( True );
          for ( int allele = 0; allele < 2; ++allele ) {
               String fasta_fn = OUT_HEAD+"."+ToString(allele+1)+".fasta";
               cout << Date() << ": writing the FASTA portion of output to " << fasta_fn << endl;
               // TODO: gap absorb needs to be pulled out of ScaflinePrinter if it's going to be here...
               String this_log = LOG.SafeBeforeLast(".") + "." + ToString(allele) + ".log";
               FastaEdgeWriter fWriter(fasta_fn, version, hbd.K(), ABBREV, True /* absorb gaps */, this_log, False /* skip gaps */, SEQ_IDS, MINSIZE);
               slp.WalkScaffoldLines( fWriter, allele );
               if(ZIP) TXT2GZ(fasta_fn);
          }
     } else FatalErr("BUG: unknown FLAVOR late in code: " + FLAVOR );

     cout << "completed conversion in " << TimeSince(start) << endl;

     return 0;
}
