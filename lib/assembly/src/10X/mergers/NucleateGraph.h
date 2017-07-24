// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

#ifndef TENX_NUCLEAT_GRAPH_H
#define TENX_NUCLEAT_GRAPH_H

#include "10X/mergers/MergerStructs.h"
#include "10X/mergers/MergerHelper.h"

// This function takes as Input
// a set of "closures" ie. sequences in <T>
// a set of "omatch" ie. overlap segments between various closures
// an involution on closures "cinv"
// an index into closures "ci"
// 
// and gives Output
// a superassembly in <T> and its involution by gluing along overlaps
// also gives optional meta data output in the form of closure_supp
//
// memory saved by supplying a dir to store intermediate calculations
// set verbose True to help debug
//
// KNOWN ISSUES/ LIMITS:
// the type T is typically in {int32,uchar}
// the type EL is the smallest datatype that 
// can store the largest length among all closures.
// Typical usage : <T,EL> = <int,uint16_t>, <uchar,uint16_t>
// The algorithm could fail if the closures aren't complete under involution
// and if there are duplicates in the closures.
// Typically does 100s of millions of gluing operations
template <class T,class EL>
void NucleateGraph(vec<vec<om<EL>>>& omatch, T hbE, vec<vec<T>>& all_closures,
        vec<int64_t>& cinv, vec<vec<int>>& ci, digraphE<vec<T>>& D, vec<int>& dinv, 
        const Bool GET_META,
        vec<vec<quad<int,EL,int,EL>>>& closure_supp,
        const Bool verbose, String dir="", Bool Canon = True);

#endif
