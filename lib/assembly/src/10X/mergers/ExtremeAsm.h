// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

#ifndef TENX_EXTREME_ASM_H
#define TENX_EXTREME_ASM_H

#include "10X/mergers/MergerStructs.h"
#include "10X/mergers/MergerHelper.h"

void ExtremeAsm( vecbasevector& reads, 
     digraphE<basevector>& Dl, vec<int>& dinvl, const int K,
     const vec<String>& controls, const Bool verbose, String out_dir="");

template <class VBV>
void ExtremeAsm( VBV& bases, const vec<int64_t>& readlist,
        digraphE<basevector>& Dl, vec<int>& dinvl, const int K,
        const vec<String>& controls, const Bool verbose, String out_dir="");

template <class VBV>
void ExtremeAsm( VBV& bases, const vec<int64_t>& readlist,
        digraphE<basevector>& Dl, vec<int>& dinvl, const int K, 
        const vec<String>& controls, const vec<nptuple>& Matches, 
        const Bool verbose, String out_dir="", unsigned int TINKER = 0);

void ExtremeAsm( vecbasevector& reads, digraphE<basevector>& Dl, 
        vec<int>& dinvl, const int K, const vec<String>& controls, 
        vec<vec<om<uint16_t>>>& omatch, const Bool verbose, String out_dir="",
        unsigned int TINKER = 0);

template <class VBV>
void Reads2CharClosures(vec<vec<unsigned char>>& all_closuresp,
        const VBV& reads);

void GetEvilMatches(vec<vec<om<uint16_t>>>& omatch, const int K,
        const vec<vec<unsigned char>>& all_closures, const vec<vec<int>>& ci,
        const vec<int64_t>& cinv, const vec<String>& controls, Bool verbose);

template <class VBV>
void CharEdges2BaseVec(const digraphE<vec<unsigned char>>& Dp, VBV& superedgeseq);

void SanitizeMatches(vec<vec<om<uint16_t>>>& omatch, vec<int64_t>& cinvp,
        const vec<vec<unsigned char>>& all_closuresp);

#endif
