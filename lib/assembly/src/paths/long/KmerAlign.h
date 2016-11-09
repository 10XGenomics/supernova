///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
#ifndef KMER_ALIGN_H
#define KMER_ALIGN_H

//  Kmer alignment between two reads. 
//  Input 
//        Kmer-matching offset list, which has the positions of the kmers in
//        read1 (pos1), and the offsets of the matching kmer position in read2 (
//        pos1 - pos2 )
//  Output:
//        kmer alignment pairs (pos2, pos1).
//
//  Method:
//    First we devide all the offsets into different groups. Within each group
//    are kmer-matchings do not have pos1 and offset difference larger than a
//    threold.  
//    Then we order the groups by the size, biggest group first. We
//    iterate through all the groups to add kmer-matchings one by one as long
//    as it does not conflict with previously added kmer-matching.  Two
//    alignments a1 and a2 are non-conflict if a1.pos1 > a2.pos1 && a1.pos2 >
//    a2.pos2 or a1.pos1 < a2.pos1 && a1.pos2 < a2.pos2.
//    
// Todo:
//   If we need more than one plausible answer for output.
void KmerAlign(
        const vec< pair<int,int> > & offset, 
        vec< pair<int,int> > & aligns, int VERBOSITY = 0 );


void KmerAlign2(
        const vec< pair<int,int> > & offset, 
        vec< pair<int,int> > & aligns, int VERBOSITY = 0 );

void KmerAlign3(
        const vec< pair<int,int> > & offset, 
        vec< pair<int,int> > & aligns, int VERBOSITY = 0 );

void GroupAndSelect( 
        const vec< pair<int,int> > & offset_in, 
        vec< pair<int,int> > & offset_out,
        int Max_Pos_Diff, double Max_Offset_Diff ,
        int VERBOSITY);

void LargestOffsetCluster( 
        const vec< pair<int,int> > & offset_in, 
        vec< pair<int,int> > & offset_out,
        int Max_Pos_Diff, double Max_Offset_Diff);


// Assessment by comparison between two kmer alignment methods
void KmerAlignCompare( 
        const vec< pair<int,int> > & aligns, 
        const vec< pair<int,int> > & ref_align ) ;

#endif
