// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

// This file is For assembly operations that will appear as part of the assembly 
// pipeline.  This operations are often templatized over MasterVec/VirtualMasterVec 
// to permit both global (full) calculations and local (tiny rapid) calculations.

#ifndef TENX_SECRET_OPS_H
#define TENX_SECRET_OPS_H

#include "CoreTools.h"
#include "feudal/PQVec.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "10X/paths/ReadPathVecX.h"

void MakeClosures( const HyperBasevectorX& hb, const vec<int>& inv,
     ReadPathVecX& paths, VecULongVec& paths_index,
     const vec<Bool>& dup, const vec<Bool>& bad,
     vec<vec<int>>& all_closures, const Bool verbose, String paths_index_file = "" );

// MarkBads.  Mark read pairs as 'bad' if at least one of the reads has > 5
// high-quality (Q30+) mismatches with the assembly.  Note that we should be
// trimming off adapter before doing this.

template< class VB, class VQ, class VP > void MarkBads( 
     // inputs:
     const HyperBasevectorX& hb, VB& bases, VQ& quals, VP& paths, 
     // output:
     vec<Bool>& bad );

template< class VB, class VQ > void MarkBads( 
     // inputs:
     const HyperBasevectorX& hb, VB& bases, VQ& quals, ReadPathVecX& paths, 
     // output:
     vec<Bool>& bad );

// IntSpaceExts.  Given an edge e, find a set of paths emanating from e that
// in principle have the property that every molecule containing e in the sample
// traverses one of the paths.  This uses only read paths containing e and their
// partners.  (The exact method needs to be written down here.)

class IntSpaceExtsWorkspace {
     public:
     vec<int> n, x, y, z;
     vec<Bool> edel, to_delete2, to_delete3;
     vec<vec<int>> ext2, ext3, pp;
     vec< pair< int64_t, vec<int> > > epo;
};

template< class VP, class VPI > void IntSpaceExts( 
     // inputs:
     const int e, const HyperBasevectorX& hb, const vec<int>& kmers, 
     const vec<int>& inv, VP& paths, VPI& paths_index, const vec<int32_t>& bc,
     const vec<Bool>& dup, const vec<Bool>& bad,
     // logging:
     const Bool GLOBAL, const Bool PRINT_RIGHT_EXTS, const int mu, 
     // output:
     vec<vec<int>>& ext, 
     // logging:
     const Bool verbose,
     // workspace: 
     IntSpaceExtsWorkspace& work );

// MarkDups.  Code to mark read PAIRS as duplicates.  If read placements share the
// same first edge and offset, and first 5 bases on the partner read, call them
// duplicates of each other.  Choose the one having the highest quality score sum
// across the pair and mark the rest as duplicates.  In case of a tie keep the pair
// having lowest id.
//
// Also flag reads that are putative artifactual (informatic) duplicates.  These 
// are reads that have identical bases and identical quality scores.

template< class VB, class VQ, class VP > void MarkDups( 
     // inputs:
     VB& bases, VQ& quals, VP& paths, const vec<int32_t>& bc,
     // outputs:
     vec<Bool>& dup, double& interdup_rate,
     // logging:
     const Bool verbose = False );

template< class VB, class VQ> void MarkDups( 
     // inputs:
     VB& bases, VQ& quals, ReadPathVecX& paths, const HyperBasevectorX& hb,
     const vec<int32_t>& bc,
     // outputs:
     vec<Bool>& dup, double& interdup_rate,
     // logging:
     const Bool verbose = False );


// AllTinksCore: returns triples
// ( e1, e2, #barcodes shared by e1 and e2 )
// but *** WARNING *** e1 is always the minimum of e1 and inv[e1].

void AllTinksCore( const HyperBasevectorX& hb, const vec<int>& inv,     
     const ReadPathVecX& paths, const vec<int64_t>& bci, const VecIntVec & ebcx,
     vec< triple<int,int,int> >& qept, const int verbosity = 1 );

// Compute barcode positions on chain chains.

double SinglesRate( const vec<int>& kmers, const vec<int>& inv, 
     const vec<Bool>& dup, const VecULongVec& paths_index, 
     const vec<int64_t>& bci );

#endif
