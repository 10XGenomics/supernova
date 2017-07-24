#ifndef TENX_GET_MERGERS_H
#define TENX_GET_MERGERS_H

#include "10X/mergers/MergerStructs.h"

void DebugPairs(qtuple Pairs, const HyperBasevectorX& hbx, const vec<vec<int>>& DLENS, 
        const vec<digraphE<vec<int>>>& DD);

void MarkOverLaps(const vec<int>& x1, const vec<int>& x2, const vec<int>& mult, 
        const HyperBasevectorX& hb, const int MIN_OLAP_KMERS, 
        const int MAX_MULT, vec<olap>& Olaps, Bool verbose=False);

void GetMergerCandidates(const vec<digraphE<vec<int>>>& DL, const vec<vec<int>>& DLI, 
        vec<vec<int>>& BCLIST, vec<vec<vec<int>>>& BPATHS, const vec<int>& bc, 
        const vec<int64_t>& bci, const HyperBasevectorX& hbx, 
        vec<triple<pair<int,int>, pair<int,int>, int>>& Matches, 
        int MIN_GLUE, int MAX_MULT, int MIN_BC_OLAP, int MIN_COMMON_BC, 
        String OUTDIR, Bool verbose );

#endif
