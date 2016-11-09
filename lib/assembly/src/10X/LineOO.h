// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

#ifndef TENX_LINE_OO_H
#define TENX_LINE_OO_H

#include "CoreTools.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "10X/IntIndex.h"
#include "10X/LineOO.h"

void ReadPosLine( const vec<int32_t>& bc, const HyperBasevectorX& hb,
     const digraphE<vec<int>>& D, const vec<int>& dinv,
     const ReadPathVec& dpaths, const IntIndex& dpaths_index,
     const vec<vec<vec<vec<int>>>>& dlines, vec<vec<triple<int32_t,int,int64_t>>>& lrpb,
     const int view );

void BarcodePos( const vec<int32_t>& bc, const HyperBasevectorX& hb,
     const digraphE<vec<int>>& D, const vec<int>& dinv,
     const ReadPathVec& dpaths, const IntIndex& dpaths_index,
     const vec<vec<vec<vec<int>>>>& dlines, vec<vec<pair<int,int>>>& lbp,
     const int view );

double ScoreOrder( 
     const vec<int>& L,                    // list of line ids
     const vec<vec< pair<int,int> >>& lbp, // barcode positions on lines
     const vec<int>& llens,                // line lengths
     vec< triple<int,int,int> >& M );      // scratch

double ScoreOrder( 
     const vec<int>& L,                    // list of line ids
     VirtualMasterVec<SerfVec<pair<int,int>>> lbpx,  // barcode positions on lines
     const vec<int>& llens,                // line lengths
     vec< triple<int,int,int> >& M );      // scratch

double MemoryScoreOrder(
     const vec<int>& L,                    // list of line ids
     const vec<vec< pair<int,int> >>& lbp, // barcode positions on lines
     const vec<int>& llens,                // line lengths
     vec< triple<int,int,int> >& M,        // scratch
     map< vec<int>, double >& memory );

void OrderN(
     // lines to OO
     const vec<int>& L0,                    // list of line ids

     // assembly info:
     const vec<vec< pair<int,int> >>& lbp, // barcode positions on lines
     const vec<int>& llens,                // line lengths

     // working data structures
     vec< triple<int,int,int> >& M,

     // answer:
     double& advantage, vec<int>& L );

void OrderN(
     // lines to OO
     const vec<int>& L0,                    // list of line ids

     // assembly info:
     VirtualMasterVec<SerfVec<pair<int,int>>> lbpx,  // barcode positions on lines
     const vec<int>& llens,                // line lengths

     // working data structures
     vec< triple<int,int,int> >& M,

     // answer:
     double& advantage, vec<int>& L );

void MemoryOrderN(
     // lines to OO
     const vec<int>& L0,                    // list of line ids

     // assembly info:
     const vec<vec< pair<int,int> >>& lbp, // barcode positions on lines
     const vec<int>& llens,                // line lengths

     // working data structures
     vec< triple<int,int,int> >& M,

     // memory
     map< vec<int>, double >& memory,

     // answer:
     double& advantage, vec<int>& L );

void MemoryOrderAndOrientN( 
     // lines to OO
     const vec<int>& L0,                    // list of line ids

     // assembly info:
     const vec<vec< pair<int,int> >>& lbp, // barcode positions on lines
     const vec<int>& llens,                // line lengths
     const vec<int>& linv,                 // line inversion

     // working data structures
     vec< triple<int,int,int> >& M,

     // memory
     map< vec<int>, double >& memory,      

     // answer:
     double& advantage, vec<int>& L );

#endif
