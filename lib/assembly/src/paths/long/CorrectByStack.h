///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// CorrectByStack takes as input a matrix whose first row is a 'founder' read
// and whose subsequent rows are 'friend' reads, or more precisely the bases on
// those friends that align to bases on the founder read.  These alignments are
// assumed to be gap-free.
//
// Output: a corrected read (both bases and quality scores), together with a
// recommended 'trim point'.

#ifndef CORRECT_BY_STACK_H
#define CORRECT_BY_STACK_H

#include "Basevector.h"
#include "CoreTools.h"
#include "Qualvector.h"
#include "paths/long/ReadStack.h"

void CorrectByStack(

     // The matrices call and callq provide bases and quality scores for the
     // founder (row 0) and its friends, for the positions aligned to the founder.  
     // Bases are 0, 1, 2, 3 (not '0',...), or ' ' (blank) in the case for a 
     // position on a friend that is before or after its end.  In the latter
     // case the quality score is shown as -1.
     //
     // For convenience, call and callq may be modified, but this is not treated
     // as a return value.

     StackBaseVecVec& call,
     StackQualVecVec& callq,

     // The new consensus is returned via bases_new and quals_new.  Their lengths
     // should not be changed.  It is assumed that upon input bases_new and 
     // quals_new are the founder.

     basevector& bases_new,
     qualvector& quals_new,
    
     // The recommended trimmed length of the founder is returned as trim_to.
     // Set it to the length of the founder if it is not to be trimmed at all.
     // Note that bases_new and quals_new should not themselves be trimmed.

     int& trim_to,

     // Triples (offset,id,rc2?) are provided, and may be modified as with
     // call and callq.  This is mostly for diagnostic purposes.

     vec< triple<int,int,Bool> >& offset_id_rc2,

     // Provide logging?

     const Bool verbose );

// TraceRead: provide view of what's happening for a particular read.

void TraceRead( StackBaseVecVec& call, StackQualVecVec& callq,
     vec< triple<int,int,Bool> >& offset_id_rc2 );

#endif
