// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

#ifndef TENX_GLUE_GRAPHS_HP_H
#define TENX_GLUE_GRAPHS_HP_H

#include "10X/mergers/MergerStructs.h"
#include "10X/mergers/MergerHelper.h"

// This function takes as Input
// a superassembly in <T> 
// a set of "Matches" ie. overlap segments between various edges of input graph
// 
// and gives Output
// a superassembly in <T> and its involution by gluing along overlaps
// also optionally gives you merger paths if GET_META = True
//
// memory saved by supplying a dir to store intermediate calculations
// set verbose True to help debug
//
// KNOWN ISSUES/ LIMITS:
// the type T is typically in {int32,uchar}
// the type EL is the smallest datatype that 
// can store the largest length among all closures.
// Typical usage : <T,EL> = <int,uint16_t>
// This function is used in Zippering and Gluing along barcode anchored overlaps
// Typically does 100s of millions of gluing operations
template <class T, class EL>
void GlueGraphs_HP( digraphE<vec<T>>& D, vec<int>& dinv, 
        vec<nptuple>& Matches, T hbE, 
        digraphE<vec<T>>& D_new, vec<int>& dinv_new, 
        const Bool GET_META, vec<vec<quad<int,EL,int,EL>>>& merger_supp,
        const Bool verbose = False, String dir="" );


#endif
