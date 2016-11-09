///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
#ifndef REFTRACE_CONTROL_H
#define REFTRACE_CONTROL_H

#include "PairsManager.h"
#include "Basevector.h"
#include "Qualvector.h"
#include "PairsManager.h"
#include "Fastavector.h"

class RefTraceControl{
public:
    RefTraceControl():sVariantOutFile(""),read_head(""),bPairsManagerInitialized(false){};
    RefTraceControl(const String& sOutFile):sVariantOutFile(sOutFile),read_head(""),bPairsManagerInitialized(false){};
    RefTraceControl(const String& sOutFile,const String& read_head) :sVariantOutFile(sOutFile),read_head(read_head),bPairsManagerInitialized(false){};
    virtual ~RefTraceControl(){};

    void setReadHead(const String& in){
        if(bPairsManagerInitialized){
            pm.clear();
            bPairsManagerInitialized=false;
        }
        read_head=in;
    };

    const String& getVariantOutFile()const{return sVariantOutFile;};

    const String& getRefHead() const { return ref_head; }
    const vec<String>& getRefTags()const{ return ref_tags; }
    const vec<size_t>& getRefIndex()const{ return ref_index; }
    const vec<size_t>& getRefStart()const{return ref_start;}
    const vec<fastavector>& getRefSeqs()const{return ref_seqs;}

    String& getRefHead() { return ref_head; }
    vec<size_t>& getRefIndex() { return ref_index; }
    vec<String>& getRefTags(){ return ref_tags; }
    vec<fastavector>& getRefSeqs(){return ref_seqs;}
    vec<size_t>& getRefStart(){return ref_start;}

    //this is NOT thread safe because PairsManager is not thread safe
    void ReadySampleLookup();                   //ready the pair manager
    // one can easily consolidate multiple library into a single sample, to be worked out later
    vec<String> getSampleList()const;           //get the full list of sample names
    int getSampleID(size_t read_ID)const ;      //get sample index from a read id
    String getSampleName(size_t read_ID)const ; //get sample name from a read id
    const vecbasevector& Reads() const { return reads; }
    const vecqualvector& Quals() const { return quals; }
    const PairsManager& GetPM() const { return pm; }

private:
    String sVariantOutFile;

    // one entry for each region entry
    String ref_head;               // reference genome sequence file (without .fasta or .fastb extension)
    vec<size_t> ref_index;         // which contig index in the reference fasta file 
    vec<String> ref_tags;          // each entry is the sequence tag in REFERENCE
    vec<fastavector> ref_seqs;     // each entry is the sequence (trimmed according to REGION)
    vec<size_t> ref_start;         // start of the trimmed region in reference sequence

    String read_head;
    bool bPairsManagerInitialized;
    PairsManager pm;
    vecbasevector reads;
    vecqualvector quals;
};

#endif
