///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Various things go here, including all abnormal termination error messages for
// Discovar (not including internal errors).  This allows us to maintain some
// consistency between these messages.

#ifndef DISCOVAR_TOOLS_H
#define DISCOVAR_TOOLS_H

#include "CoreTools.h"
#include "Fastavector.h"
#include "paths/HyperBasevector.h"
#include "paths/long/Logging.h"
#include "paths/long/RefTrace.h"
#include "paths/long/RefTraceControl.h"

namespace DiscovarTools{

    // Exit 0 upon finding that the assembly is empty, or there are no reads.

    void ExitAssemblyEmpty( );
    void ExitPathsEmpty( );
    void ExitNoReads( );
    void ExitNoCorrectedReads( );
    void ExitShortReads( const String& additional_info="" );

    // Final message and exit 1 when something really wrong.

    void DiscovarUnhappy( );

    void CheckDiscovarSystemRequirements( );

    void CheckDiscovarBams( const vec<String>& bams );

    //check if REGIONS input is valid or not
    void CheckDiscovarRegionsInput(const String& REGIONS);

    //check if parsed REGIONS input is valid or not
    void CheckDiscovarRegions( const String& REGIONS );

    //check if parsed OUT_HEAD input is valid or not
    void CheckDiscovarOutHead( const String& OUT_HEAD );

    //check if parsed TMP input is valid or not
    void CheckDiscovarTmp( const String& TMP );

    //check if parsed READS input is valid or not
    void CheckDiscovarReads( const String& READS );

    void TestDiscovarRegionsBamsCompatibility( const String& REGIONS,
         const vec<String>& bams );

    void ExitSamtoolsFailed( );

    void CheckReferenceInput( const String& REFERENCE, const String& OUT_HEAD );

    // Forbids region longer than certain length
    const size_t MaxRegionSize=50000000; // Maximum allowable region size
    void CheckRegionSize(size_t size,const String&sPrefix=""); // helper function: check if 'size' is longer than MaxRegionSize

    // Forbids total BAM file size over certain number of bytes, typically called when REGIONS="all", which is then translated to REGIONS==""
    const size_t MaxBAMSize=10ULL*1024*1024*1024; // Maximum allowable BAM size
    void CheckBAMSize(size_t size,const String& REGIONS); // helper function: check if 'size' is larger than MaxBAMSize if REGIONS=="all" or ""


    //control data structure for ref trace/ variant calling
    class DiscovarRefTraceControl: public RefTraceControl{
    public:
        typedef RefTraceControl base_t;
        DiscovarRefTraceControl():bRefTraceOn(false){};
        DiscovarRefTraceControl( const String& REF_FASTA
                               , const String& REGIONS
                               , const long_logging& logc
                               , const String& sVariantOutFile_);

        //ForceAssert that data size are consistent
        void CheckConsistency() const;

        bool IsRefTraceOn()const {return bRefTraceOn;}
        const vec<fastavector>& getRefSeqsExt()const{return ref_seqs_ext;}
        const vec<size_t>& getRefSeqsExtExt()const{return ref_seqs_ext_ext;}
        const vec<size_t>& getRefLength()const{return ref_length;}
        const vec<size_t>& getChromIdx()const{return RefTraceControl::getRefIndex();}

        void AssertEqual(const DiscovarRefTraceControl&other);

    private:
        bool bRefTraceOn;              //if Ref Trace functionality is on or not

        // one entry for each region entry, or if region is not supplied, one entry for one entry in REFERENCE file
        vec<fastavector> ref_seqs_ext; // each entry is the sequence (trimmed according to REGION, plus up to EXT_LENGTH more bases on both ends)
        vec<size_t> ref_seqs_ext_ext; // the start of ref_seqs relative to ref_seqs_ext
        vec<size_t> ref_length;     // length of trimmed region in reference sequence

        static const int MIN_LENGTH = 1000; // minimum allowable reference length
        static const int EXT_LENGTH = 500; // try to expand region by EXT_LENGTH on both ends for ref_seqs_ext, for LongProto compatibility
    };
}//namespace DiscovarTools

#endif
