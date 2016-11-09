// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
//
// Parse special barcoded fastq files.  These are not regulation fastqs.
//
// FASTQ lines:
// - read1, qual1, read2, qual2, barcode, barcode_qual, sample index, sample index qual
// - no plus (+) signs.  Line ending signifies the end of each of the above entries in a
//   record.
// - PHRED offset TBD (+33, +64)

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "FastIfstream.h"
#include "MainTools.h"
#include "ParallelVecUtilities.h"
#include "Qualvector.h"
#include "feudal/PQVec.h"
#include "kmers/KmerRecord.h"
#include "feudal/BinaryStream.h"
#include "feudal/IncrementalWriter.h"
#include "feudal/VirtualMasterVec.h"
#include <fstream>
#include <string>
#include <array>

namespace {

struct ParseError : std::runtime_error {
     ParseError(String const& s) : std::runtime_error(s) {};
     ParseError(string const& s) : std::runtime_error(s) {};
     ParseError(char const * s) : std::runtime_error(s) {};
};
class EndOfFile {};

void convertPhred( String const& buf, qualvector& qtmp, const uchar phred_offset=33 )
{
     qtmp.clear();
     for ( char const c : buf  ) {
          if ( c == '\n' || c == '\r' ) continue;
          qtmp.push_back( c - phred_offset );
     }
}


bool hasGemGroup( String const& s )
{
     // presence of a dash -- searching from end for efficiency
     auto itr=s.end();
     while ( itr-- != s.begin() ) if ( *itr == '-' ) return true;
     return false;
}


void newUnpackBarcodeSortedFastq( const String& file, vecbasevector& b0_bases,
          VecPQVec& b0_quals, vecbasevector& bc_bases, VecPQVec& bc_quals,
          vec<uint32_t>& bcs )
{
     // TODO: work for non-gzipped
     ForceAssert(file.EndsWith(".fastq.gz"));

     cout << Date() << ": " << file << endl;

     String buf;
     fast_pipe_ifstream input( "zcat " + file );
     enum { START, READ1, QUAL1, READ2, QUAL2, BARC, QUALBARC, INDEX, QUALINDEX, END };
     array<basevector,2> btmp;
     array<qualvector,2> qtmp;
     uint32_t barc = bcs.size() ? bcs.back() : 0u;
     String lastb;
     size_t line = 0;
     std::map<String,size_t> seen;

     while ( 1 ) {
          try {
               // name
               getline(input, buf); line++;
               if ( input.fail() ) throw EndOfFile();
               if ( !buf.StartsWith("@") ) throw ParseError(buf);

               // read a block
               for ( size_t blocki = START+1; blocki < END; ++blocki  ) {
                    getline(input, buf); line++;
                    if ( input.fail() ) throw ParseError(buf);

                    switch( blocki ) {
                         case READ1:
                              btmp[0].SetFromStringWithNs( buf );
                              break;
                         case READ2:
                              btmp[1].SetFromStringWithNs( buf );
                              break;
                         case QUAL1:
                              convertPhred( buf, qtmp[0] );
                              break;
                         case QUAL2:
                              convertPhred( buf, qtmp[1] );
                              break;
                         case BARC:
                              // lack of a gem group indicates failing barcode
                              // presence of a gem group, but lack of a barcode indicates no barcode read
                              // AND no whitelist.  Treat them the same.
                              if ( hasGemGroup( buf ) && buf[0] != '-' ) {
                                   if ( lastb != buf ) {    // should trigger 1st time
                                        barc++;             // ...so 1st bc is 1
                                        lastb = buf;
                                        if ( seen.count(buf) > 0 ) {
                                             cout  << "barcode " << buf << " at line " << line << " already seen at line "  << seen[buf] << endl;
                                        }
                                        seen[buf]=line;
//                                        PRINT3( barc, seen.size(), buf );
                                   }
                                   bcs.push_back(barc);
                                   bcs.push_back(barc);
                                   bc_bases.push_back( btmp[0] );
                                   bc_bases.push_back( btmp[1] );
                                   bc_quals.push_back( PQVec( qtmp[0] ) );
                                   bc_quals.push_back( PQVec( qtmp[1] ) );
                              } else {
                                   b0_bases.push_back( btmp[0] );
                                   b0_bases.push_back( btmp[1] );
                                   b0_quals.push_back( PQVec( qtmp[0] ) );
                                   b0_quals.push_back( PQVec( qtmp[1] ) );
                              }

                              break;
                         case QUALBARC:
                         case INDEX:
                         case QUALINDEX:
                         default:
                              break;
                    }
               }

          } catch ( ParseError const& e ) {
               FatalErr( "out of sync reading line: " + string(e.what()) );
          } catch ( EndOfFile const& e ) {
               break;
          }
     }
}


void readIndirect( String const& filename, vec<String>& v )
{
     Ifstream(INPUT, filename);
     String buf;
     while (1)  {
          getline(INPUT, buf);
          if ( INPUT.fail() ) break;
          if ( buf.size() > 0 ) {
               buf.resize( buf.size()-1 );
               v.push_back(buf);
          }
     }
}

void mergeBarcodedReadFiles( vec<String> const& MERGE_HEADS, String const& OUT_HEAD )
{
     // parse input, including any indirect files
     vec<String> merge_heads;
     for ( auto const& merge_head : MERGE_HEADS ) {
          ForceAssertGt(merge_head.size(), 0u);
          if ( merge_head[0] == '@' ) readIndirect( merge_head.After("@"), merge_heads );
          else merge_heads.push_back( merge_head );
     }

     // create output files
     BinaryIteratingWriter<vec<int64_t>>     bci_out( OUT_HEAD + ".bci" );
     IncrementalWriter<basevector>           bases_out( OUT_HEAD + ".fastb" );
     IncrementalWriter<PQVec>                quals_out( OUT_HEAD + ".qualp" );

     // two passes -- do zero bc first because they must be aggregated.
     // all other files are assumed to have unique barcodes
     size_t count = 0;
     uint32_t bc_bias = 0;
     bci_out.write(0u);

     for (unsigned pass = 0; pass <=1; ++pass ) {
          for ( std::string merge_head : merge_heads ) {

               cout << Date() << ": pass " << pass+1 << " (of 2): " << merge_head << endl;

               VirtualMasterVec<basevector> bases( merge_head + ".fastb" );
               VirtualMasterVec<PQVec> quals( merge_head + ".qualp" );
               vec<int64_t> bci; BinaryReader::readFile( merge_head + ".bci", &bci );

               ForceAssertEq( bases.size(), quals.size() );
               ForceAssertEq( bases.size(), bci.back() );

               if ( bases.size() == 0 ) continue;

               for ( size_t i = 1; i < bci.size(); ++i )
                    if ( bci[i] < bci[i-1] ) {
                         cout << "barcode indices are out of order" << endl;
                         PRINT2(bci[i], bci[i-1]);
                    }

               if ( pass == 0 ) {
                    size_t this_count = bci[1];
                    for ( size_t i = 0; i < this_count; ++i ) {
                         bases_out.add( bases[i] );
                         quals_out.add( quals[i] );
                    }

                    count += this_count;
               } else {
                    size_t zero_count = bci[1];

                    // increment all of the offsets
                    ForceAssertGe( count, zero_count );

                    for ( size_t i = 1; i < bci.size(); ++i )
                         bci[i] = bci[i] - zero_count + count;

                    // write all reads and barcodes past bc0
                    for ( size_t i = zero_count; i < bases.size(); ++i )  {
                         bases_out.add( bases[i] );
                         quals_out.add( quals[i] );
                    }

                    for ( size_t i=1; i < bci.size()-1; ++i )
                         bci_out.write( bci[i] );

                    count += bases.size() - zero_count;
               }
          }

          if ( pass == 1 ) bci_out.write( count );
     }
}


};

int main(int argc, char *argv[])
{
     double clock = WallClockTime();

     RunTime( );

     BeginCommandArguments;
     CommandArgument_StringSet_OrDefault_Doc(FASTQS, "", "list of barcode-fastq files with reads\n"
     "will consider a read to be 'unbarcoded' if there is no appended gem group (e.g. BBBBBBBBBBBBBB-1)");
     CommandArgument_String_Doc(OUT_HEAD, "basename of output files {.fastb, .qualp, .bc, .bci}");
     CommandArgument_StringSet_OrDefault_Doc(MERGE_HEADS, "",
               "list of heads for fastb/qualp/bc/bci files to merge, rather than produce" );
     EndCommandArguments;

     // Define data structures.

     if ( MERGE_HEADS.size() ) {

          // MERGE path
          if ( FASTQS.size() ) FatalErr( "must not supply FASTQ/FASTQ0 with MERGE_HEADS" );

          mergeBarcodedReadFiles( MERGE_HEADS, OUT_HEAD );

     } else {
          // FASTQ->FASTB/QUALP path

          if ( !FASTQS.size() ) FatalErr( "no inputs specified FASTQ or FASTQ0" );

          vec<String> fastqs;

          for ( auto const& file : FASTQS )
               if ( file[0] == '@' ) readIndirect(file.After("@"), fastqs );
               else fastqs.push_back(file);


          // we're assuming that barcodes are not split across files
          vecbasevector bc_bases,b0_bases;
          VecPQVec bc_quals,b0_quals;
          vec<uint32_t> bcs;
          for ( auto const& file : FASTQS )
               newUnpackBarcodeSortedFastq( file, b0_bases, b0_quals, bc_bases, bc_quals, bcs );


          vec<int64_t> bci;
          bci.push_back(0);
          int64_t b0_size = b0_bases.size();
          int64_t cur = 0;
          for (size_t i = 0; i < bcs.size(); ++i )
               if ( bcs[i] != cur ) {        // bcs should start != 0, so we cover the end of 0 here
                    cur = bcs[i];
                    bci.push_back(i+b0_size);
               }
          bci.push_back(b0_size+bcs.size());


          cout << Date() << ": writing output to " + OUT_HEAD + " .fastb,.qualp,.bci " << endl;
          IncrementalWriter<basevector> bases_out( OUT_HEAD + ".fastb" );
          IncrementalWriter<PQVec>      quals_out( OUT_HEAD + ".qualp" );
          bases_out.add( b0_bases.begin(), b0_bases.end() );
          bases_out.add( bc_bases.begin(), bc_bases.end() );
          quals_out.add( b0_quals.begin(), b0_quals.end() );
          quals_out.add( bc_quals.begin(), bc_quals.end() );
          BinaryWriter::writeFile( OUT_HEAD + ".bci", bci ) ;

          // barcode stats
          size_t count = 0;
          for ( size_t i = 1; i < bci.size()-1; ++i )  {
               if ( bci[i+1] - bci[i] > 10 ) count++;
          }
          cout << Date() << ": " << count << " barcodes have more than 10 reads, out of a total of "
               << bci.size()-2 << " barcodes" << endl;

          // Done.
     }

     cout << TimeSince(clock) << " used, peak mem = " << PeakMemUsageGBString( )
          << endl;
     cout << Date( ) << ": done" << endl;
     Scram(0);
}
