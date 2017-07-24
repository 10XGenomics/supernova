// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

#ifndef TENX_GLUE_H
#define TENX_GLUE_H

#include "Basevector.h"
#include "CoreTools.h"
#include "graph/Digraph.h"

// Glue: return True if glue rule known.

Bool Glue( VirtualMasterVec<basevector>& bases, const vec<int64_t>& ids, 
     const String& glue_rules, digraphE<basevector>& ghb, vec<int>& ginv );

void SeqSelfAlignments(String& x1,const int MIN_WT, 
        const int top, ostringstream& out);

void SeqAlignments(String& x1, String& x2, const int MIN_WT, 
        const int top, ostringstream& out);

// Functions for finding complete perfect representations of a given sequence in a 
// gluegraph.  The function IsPath tells you whether there is a path.  The 
// function GetPaths will find all the paths, but be careful what you wish for!
// The number of paths can grow exponentially, and even if it doesn't, GetPaths can
// be very slow.  The function GetPathCov finds the union of all the edges in all
// the paths (without actually finding the paths).  This is much better behaved.

Bool IsPath( const digraphE<basevector>& ghb, const basevector& b );

void GetPaths( const digraphE<basevector>& ghb, const basevector& b,
     vec<vec<int>>& p );

void GetPathCov( const digraphE<basevector>& ghb, const basevector& b,
     vec<int>& cov );

#endif
