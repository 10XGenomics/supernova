#ifndef _TENX_DISP_G_H
#define _TENX_DISP_G_H

#include "Basevector.h"
#include "CoreTools.h"
#include "graph/Digraph.h"
#include "paths/long/ReadPath.h"
#include "graph/DigraphTemplate.h"
#include "10X/Super.h"
#include "10X/Glue.h"

template <class F>
void PrintGraphPaths(const digraphE<F>& ghb, const vec<vec<int>>& gold,
        String fname, vec<String> color);

template <class F>
void PrintGraphPaths(const digraphE<F>& ghb, const vec<vec<int>>& gold,
        const vec<int>& meta, String fname, vec<String> color);

template <class F>
void DisplayGraphPaths(const digraphE<F>& dhb, const ReadPathVec& dpaths,
                 const vec<vec<int>>& dpaths_index, const vec<int> IDs, Bool verbose, String text, String fname);

template <class F>
void DisplayGraphPaths(const digraphE<F>& dhb, const vec<vec<int>>& Gold, const vec<String>& colors,
                 const vec<vec<int>>& dpaths_index, Bool verbose, String fname);

void DisplayRefGraphPath(const digraphE<bvec>& dhb, const vec<bvec>& ref, 
        const vec<vec<int>>& dpaths_index, Bool verbose, String fname);

template <class E, class F>
void DisplaySCC(const digraphE<E>& D, const digraphE<F>& dhb, 
        const vec<vec<int>>& dpaths_index, Bool verbose, String fname);
#endif
