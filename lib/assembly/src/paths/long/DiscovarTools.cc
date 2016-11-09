///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Discovar messages are tagged "DISCOVAR MESSAGE".

#include <omp.h>

#include "CoreTools.h"
#include "FastIfstream.h"
#include "ParseSet.h"
#include "VecUtilities.h"
#include "lookup/SAM.h"
#include "paths/long/DiscovarTools.h"

namespace DiscovarTools{

void DiscovarUnhappy( )
{    // DISCOVAR MESSAGE
     cout << "Sorry, Discovar cannot proceed.\n" << endl;
     Scram(1);    }

void ExitAssemblyEmpty( )
{    // DISCOVAR MESSAGE
     cout << "\nDear user, we are sad to report that your assembly is now "
          << "empty.\nThis could be because your coverage is too low, or "
          << "that something\nelse is wrong with your data.  Assembly of very "
          << "small regions\n(e.g. ~1 kb) can also cause empty assemblies.\n"
          << "Discovar complete.  Have a nice day.\n\n";
     Scram( );    }

void ExitPathsEmpty( )
{    // DISCOVAR MESSAGE
     cout << "\nDear user, we are sad to report the current assembly is not supported by reads.\n"
          << "This could happen if the reads do not align well with the reference (if supplied), or "
          << "that something\nelse is wrong with your data.\n"
          << "Discovar complete.  Have a nice day.\n\n";
     Scram( );    }

void ExitNoReads( )
{    // DISCOVAR MESSAGE
     cout << "\nDear user, we are sad to report that no reads were found.\n"
          << "Further assembly is thus impossible.\n"
          << "Discovar complete.  Have a nice day.\n\n";
     Scram( );    }

void ExitNoCorrectedReads( )
{    // DISCOVAR MESSAGE
     cout << "\nDear user, we are sad to report that after error correction, "
          << "no reads remain.\n"
          << "Further assembly is thus impossible.\n"
          << "Discovar complete.  Have a nice day.\n\n";
     Scram( );    }

void ExitSamtoolsFailed( )
{    // DISCOVAR MESSAGE
     cout << "\nSamtools failed, Discovar cannot proceed.\n" << endl;
     Scram(1);    }

void ExitShortReads(const String& additional_info )
{
    // DISCOVAR MESSAGE
    cout << "\nDiscovar has found that all reads in your input are too short";
    if ( additional_info != "" )
        cout << ":\n" << additional_info << "\n";
    cout << "\nDiscovar complete.  Have a nice day.\n\n";
    Scram( );    }

void CheckDiscovarSystemRequirements( )
{
     // Check samtools existence and version.

     if ( System( "which samtools > /dev/null 2>&1" ) != 0 )
     {    // DISCOVAR MESSAGE
          cout << "I can't find the samtools executable.  It may be that you need\n"
               << "to install it on your system.  Or perhaps you just need to add\n"
               << "it to your path.\n";
          DiscovarUnhappy( );    }
     String samvs = StringOfOutput( "samtools 2>&1 | grep Version", 2 );

     // Parse samvs as v1.v2.v3... where v1-3 are nonnegative integers and ... is 
     // anything.  Or read as v1.v2.

     int samv1 = 0, samv2 = 0, samv3 = 0;
     const triple<int,int,int> min_samv = {0,1,18};
     Bool confused = False;

     // First deal with the v1.v2 case.

     if ( samvs.Contains( "." ) && samvs.Before( "." ).IsInt( )
          && samvs.After( "." ).IsInt( ) && !samvs.After( "." ).Contains( "." ) )
     {    samv1 = samvs.Before( "." ).Int( ), samv2 = samvs.After( "." ).Int( );    }

     // Now the other cases.

     else
     {    if ( !samvs.Contains( "." ) || !samvs.After( "." ).Contains( "." )
               || !samvs.Before( "." ).IsInt( ) 
               || !samvs.Between( ".", "." ).IsInt( ) )
          {    confused = True;    }
          String samt;
          if ( !confused )
          {    samt = samvs.After( "." ).After( "." );
               int n;
               for ( n = 0; n < samt.isize( ); n++ )
                    if ( !isdigit( samt[n] ) ) break;
               samt.resize(n);
               if ( samt.empty( ) || !samt.IsInt( ) ) confused = True;   }
          if (confused)
          {    cout << "I can't determine which version of samtools you're using.\n";
               cout << "Version read as: " << samvs << "." << endl;
               cout << "Perhaps you have an old version." << endl;
               cout << "Please use at least version " << min_samv.first << "."
                    << min_samv.second << "." << min_samv.third << "." << endl;
               DiscovarUnhappy( );    }
          samv1 = samvs.Before( "." ).Int( );
          samv2 = samvs.Between( ".", "." ).Int( );
          samv3 = samt.Int( );    }
     triple<int,int,int> samv = {samv1,samv2,samv3};
     if ( samv < min_samv )
     {    // DISCOVAR MESSAGE
          cout << "You appear to have version " << samvs << " of samtools.\n"
               << "Please use at least version " << min_samv.first << "."
               << min_samv.second << "." << min_samv.third << "." << endl;
          DiscovarUnhappy( );    }    }

void CheckDiscovarBams( const vec<String>& bams )
{    for ( const String& b : bams )
     {    if ( !b.Contains( ".bam", -1 ) )
          {    // DISCOVAR MESSAGE
               cout << "The file " << b << " that you provided as part of READS "
                    << "does not end in .bam.\n";
               DiscovarUnhappy( );    }
          if ( !IsRegularFile(b) )
          {    // DISCOVAR MESSAGE
               cout << "I can't find the file " << b << " that you provided as "
                    << "part of READS.\n";
               DiscovarUnhappy( );    }
          if ( System( "samtools view -H " + b + " > /dev/null 2>&1" ) != 0 )
          {    // DISCOVAR MESSAGE
               cout << "When I try to view " << b << " with samtools, something\n"
                    << "goes wrong.  Please check to see that your input files\n"
                    << "have been created correctly." << endl;
               DiscovarUnhappy( );    }    }    }

void CheckDiscovarOutHead( const String& OUT_HEAD ){
    if( OUT_HEAD.Contains("/")){
      String directory = OUT_HEAD.SafeBeforeLast("/");
      if( ! IsDirectory(directory)){
        std::cout << "The directory " << directory << " in OUT_HEAD=" << OUT_HEAD << " does not exist."<< std::endl;
        DiscovarUnhappy();
      }
    }
}

void CheckDiscovarTmp( const String& TMP ){
    if ( IsRegularFile(TMP) ){
        std::cout << "Temporarily directory TMP="<<TMP << " already exists as a regular file."<< std::endl;
        DiscovarUnhappy();
    }
}

void CheckDiscovarReads( const String& READS ){
    if ( READS.Contains("~") ){
        std::cout << "The READS argument cannot contain '~'. Please use full paths instead." << std::endl;
        DiscovarUnhappy();
    }
}

void CheckDiscovarRegions( const String& REGIONS )
{    vec<String> regions;
     size_t region_size=0;
     ParseStringSet( "{"+REGIONS +"}", regions );
     for ( auto r : regions )
     {    if ( r.Contains( ":" ) )
          {    String range = r.After( ":" );
               Bool bad = False;
               if ( !range.Contains( "-" ) ) bad = True;
               if ( !bad )
               {    if ( !range.Before( "-" ).IsInt( ) 
                    || !range.After( "-" ).IsInt( ) )
                    bad = True;    }
               if ( !bad )
               {    int start = range.Before( "-" ).Int( );
                    int stop = range.After( "-" ).Int( );
                    if (  start > stop  || stop <0 || start<0){
                        bad = True;
                    }
                    else{
                        region_size += stop - start;
                    }
               }
               if ( !bad )
                    for ( auto c : range ) if ( isalpha(c) ) bad = True;
               if (bad)
               {    // DISCOVAR MESSAGE
                    cout << "REGIONS argument " << r << " doesn't make sense: "
                         << "there needs to be a\nstart-stop integer range after "
                         << "the colon, with start <= stop." << endl;
                    DiscovarUnhappy( );    }    }
     }
     CheckRegionSize(region_size);
}
               
void TestDiscovarRegionsBamsCompatibility( const String& REGIONS, 
     const vec<String>& bams )
{    vec<String> regions;
     ParseStringSet( "{" + REGIONS + "}", regions );
     size_t totalFileSize=0;
     for ( auto b : bams )
     {    SAM::BAMFile bf( b, "", true );
          vec<RefDesc> rd = bf.getRefDescs( );
          vec<String> refnames;
          vec<int> reflens;
          for ( int i = 0; i < rd.isize( ); i++ )
          {    refnames.push_back( rd[i].getName( ) );
               reflens.push_back( rd[i].getLength( ) );    }
          UniqueSortSync( refnames, reflens );

          size_t region_size=0;
          for ( auto region : regions )
          {    String ref = region;
               if ( ref.Contains( ":" ) ) ref = region.Before( ":" );
               int p = BinPosition( refnames, ref );
               if ( p < 0 )
               {    // DISCOVAR MESSAGE
                    cout << "There appears to be an incompatibility between your "
                         << "READS and REGIONS arguments.\nSpecifically, the "
                         << "region " << region << " refers to reference record\n"
                         << ref << ", which is not declared in the header for bam "
                         << "file\n\n" << b << "." << "\n" << endl;
                    DiscovarUnhappy( );    }
               if ( region.Contains( ":" ) )
               {    String range = region.After( ":" );
                    int64_t stop = range.After( "-" ).Int( );
                    if ( stop > reflens[p] )
                    {    // DISCOVAR MESSAGE
                         cout << "It would appear that the stop position for region "
                              << region << "\nextends beyond the end of the "
                              << "reference sequence." << endl;
                         DiscovarUnhappy( );    }
                    region_size += stop-region.After(":").Before("-").Int();
               }
               else{
                    region_size += reflens[p];
               }
          }

          {
              ifstream file(b.c_str(),ios::in|ios::binary|ios::ate);
              totalFileSize+=file.tellg();
              file.close();
          }
          CheckRegionSize(region_size);
     }
     CheckBAMSize(totalFileSize,REGIONS);
}

          /*
          // fast_pipe_ifstream in( "samtools view -H " + b );
          // String line;
          // while(1)
          // {    getline( in, line );
          //      if ( in.
          */

               // }    }

//size_t DiscovarRefTraceControl::MIN_LENGTH = 1000;

void CheckReferenceInput( const String& REFERENCE, const String& OUT_HEAD ){
    if(REFERENCE=="") return;
    if ( !IsRegularFile(REFERENCE) ){
        std::cout << "REFERENCE=" << REFERENCE << " is not a regular file."<< std::endl;
        DiscovarUnhappy();
    }
    String REF_BARE = REFERENCE.SafeAfterLast("/");
    String OUT_BARE = OUT_HEAD.SafeAfterLast("/");
    if(REF_BARE.Contains(".")){
        if( REF_BARE.SafeBeforeLast(".") == OUT_BARE){
            std::cout << "Please refrain from naming REFERENCE with a prefix identical to OUT_HEAD, as it might get overwritten.\n" << std::endl;
            DiscovarUnhappy( );
        }
    }
};

void DiscovarRefTraceControl::CheckConsistency() const{
    bool bConsistent =   getRefSeqs().size()==getRefTags().size()
                      && ref_seqs_ext.size() == getRefSeqs().size()
                      && ref_seqs_ext_ext.size() == ref_seqs_ext.size()
                      && getRefStart().size() == ref_seqs_ext_ext.size()
                      && ref_length.size() == getRefStart().size()
                      && getRefIndex().size() == ref_length.size();
    ForceAssert(bConsistent);
    for(size_t ii=0; ii<getRefSeqs().size(); ++ii){
        ForceAssert( ref_seqs_ext[ii].size() >= getRefSeqs()[ii].size() );
        if( ref_seqs_ext[ii].size() > getRefSeqs()[ii].size() + 2*EXT_LENGTH){
            cout << "WARNING: extended reference sequence is too long." << std::endl;
        }
    }
};

DiscovarRefTraceControl::DiscovarRefTraceControl(const String& REF_FASTA, const String& REGIONS, const long_logging& logc
                                                , const String& sVariantOutFile_
                                                ):base_t(sVariantOutFile_)
{
    bRefTraceOn = REF_FASTA != "";
    if(bRefTraceOn) {
        if( REF_FASTA == fastaindex::filename_fix(REF_FASTA)){
            std::cout << "Please create a FASTA index for the reference genome." << std::endl;
            std::cout << "    e.g.    samtools faidx " << REF_FASTA << std::endl << std::endl;
            DiscovarUnhappy( );
        }
        CheckDiscovarRegions(REGIONS);
        vec<String> regions;
        ParseStringSet( "{" + REGIONS + "}", regions );

        fastaindex fai( REF_FASTA );


        if (logc.STATUS_LOGGING)
        {    std::cout << Date() << ": RefTrace is on. Reading " << REF_FASTA 
                  << " for REGIONS=" << REGIONS << std::endl;    }

        size_t region_size=0;

        for(const auto& entry: regions){
            String chr;
	    int start=-1,end=-1;
            if( entry.Contains(":")){
        	chr = entry.Before(":");
        	//note that this is using using base-1 coordinate as base-0 coordinate.
        	//This is a decision made by the group to achieve compliance with LongProto
                start = entry.After(":").Before("-").Int();
                end = entry.After(":").After("-").Int();
            }
            else{
        	chr = entry;
                start=0;
                end = -1;
            }
            if ( !fai.has_key(chr) ) {    // DISCOVAR MESSAGE
                cout << "There appears to be an incompatibility between your "
                     << "REF_FASTA and REGIONS arguments.\nSpecifically, the "
                     << "region " << entry << " refers to a reference record,\n"
                     << "which is not found in the FASTA "
                     << "file\n\n" << REF_FASTA << "." << "\n" << endl;
                DiscovarUnhappy( );
            }
            int ilen = fai[chr].len;

            if( end!=-1){
                //note that this is a hack for using base-1 coordinate as base-0 coordinate.
                //This is a decision made by the group to achieve compliance with LongProto
                if (ilen==end+1) {--end;};
            }
            else{
        	end = ilen;
            }
            if( end > ilen ) {
                cout << "There appears to be an incompatibility between your "
                     << "REF_FASTA and REGIONS arguments.\nSpecifically, the "
                     << "region " << entry << " refers to a sequence,\n"
                     << "which is longer than that specified in the FASTA "
                     << "file\n\n" << REF_FASTA << "." << "\n" << endl;
                DiscovarUnhappy( );
            }


            if( end-start < MIN_LENGTH){
                std::cout << "Reference sequence shorter than " << MIN_LENGTH << "-base is not allowed." << std::endl;
                DiscovarUnhappy( );

            }

            getRefTags().push_back(chr);

            if (logc.STATUS_LOGGING)
            {
            cout << Date() << ": about to load from indexed Fasta: " << getRefTags().back() << ":" <<
        	    start << "-" << end << endl;
            }
            int ext_start = std::max(0,start-EXT_LENGTH);
            int ext_end = std::min(ilen,end+EXT_LENGTH);

            fastavector fv;
            LoadRegionFromIndexedFastaFile(REF_FASTA,fai, chr, ext_start, ext_end, fv);
            ref_seqs_ext.push_back(fv);
            ref_seqs_ext_ext.push_back(start-ext_start);

            getRefSeqs().push_back(fastavector());
            getRefSeqs().back().SetToSubOf(ref_seqs_ext.back(),start-ext_start,end-start);

            region_size+=end-start;
            getRefIndex().push_back(fai[chr].idx);
            getRefStart().push_back(start);
            ref_length.push_back(end-start);
        }

        // if region not specified
        if( getRefSeqs().size()==0){
	    vec<fastavector> loc_seqs;
	    vec<String> loc_tags;
	    LoadFromFastaFile(REF_FASTA,loc_seqs,loc_tags);

            getRefStart().resize(loc_tags.size(),0);
            ref_length.resize(loc_tags.size());
            getRefIndex().resize(loc_seqs.size());
            for(size_t ii=0;ii<loc_seqs.size();++ii){
                getRefIndex()[ii]=ii;
                getRefStart()[ii]=0;
                ref_length[ii]=loc_seqs[ii].size();
                region_size+=loc_seqs[ii].size();
            }
            ref_seqs_ext = loc_seqs;
            ref_seqs_ext_ext.assign( ref_seqs_ext.size(),0);
            swap(loc_seqs,getRefSeqs());
            swap(loc_tags,getRefTags());
        }
        CheckRegionSize(region_size,"When reference tracing is on, ");
    }
    CheckConsistency();
};

void CheckRegionSize(size_t size, const String& sPrefix){
    if( size > DiscovarTools::MaxRegionSize){
        std::cout << sPrefix
                  << "REGIONS totaled more than " << DiscovarTools::MaxRegionSize << "\n"
                  << " bases, which is not officially supported at this time."<< std::endl;
        DiscovarUnhappy();
    }
}

void CheckDiscovarRegionsInput(const String& REGIONS){

    if( REGIONS==""){
        std::cout << "REGIONS cannot be an empty string." << std::endl;
        DiscovarUnhappy();
    }

}

void CheckBAMSize(size_t size, const String& REGIONS){
    if( size > DiscovarTools::MaxBAMSize && (REGIONS=="" || REGIONS=="all")){
        std::cout << "Using no REGIONS restriction when the BAM files total more than " << DiscovarTools::MaxBAMSize << "\n"
                  << "bytes, which is not officially supported at this time."<< std::endl;
        DiscovarUnhappy();
    }
}

void DiscovarRefTraceControl::AssertEqual(const DiscovarRefTraceControl&other){
    ForceAssert(bRefTraceOn==other.bRefTraceOn);
    ForceAssert(getRefTags()==other.getRefTags());
    ForceAssert(getRefSeqs()==other.getRefSeqs());
    ForceAssert(getRefSeqs().size()==other.getRefSeqs().size());
    ForceAssert(ref_seqs_ext==other.ref_seqs_ext);
    ForceAssert(ref_seqs_ext_ext==other.ref_seqs_ext_ext);
    ForceAssert(getRefStart()==other.getRefStart());
    ForceAssert(ref_length==other.ref_length);
    ForceAssert(getRefIndex()==other.getRefIndex());
};
}


