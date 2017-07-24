// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

// RefAlign.  Align a supernova assembly to a reference sequence.  Intended for 
// aligning a human assembly to hg19.

#ifndef TENX_ASTATS_REF_ALIGN_H
#define TENX_ASTATS_REF_ALIGN_H

#include "CoreTools.h"
#include "PackAlign.h"
#include "feudal/SmallVecDefs.h"
#include "feudal/OuterVecDefs.h"
#include "kmers/KmerRecord.h"
#include "paths/HyperBasevector.h"

class refalign {

     public:

     refalign( ) { }
     refalign( const int chr, const align& a ) : chr(chr), a(a) { }
     int chr;   // e.g. +1 = chr1, -1 = cr of chr1
     align a;

     void writeBinary(BinaryWriter& writer ) const
     {    writer.write(chr);
          writer.write(a);    }

     void readBinary(BinaryReader& reader)
     {    reader.read(&chr);
          reader.read(&a);    }

     static size_t externalSizeof( ) { return 0; }

};

SELF_SERIALIZABLE(refalign);

extern template class SmallVec<refalign,MempoolAllocator<refalign> >;
extern template class OuterVec< SerfVec<refalign> >;
template class SmallVec<refalign,MempoolAllocator<refalign> >;
template class OuterVec< SerfVec<refalign> >;

template<int K> class RefAlignScratch {

     public:

     RefAlignScratch( ) 
     {    locs.resize(2);
          blob.resize(2);    }

     vec<vec<triple<int,int,int>>> locs;
     vec< triple< pair<int,int>, triple< int, Bool, int >, String > > valigns;
     vec< pair<ho_interval,ho_interval> > M;
     vec< triple<kmer<K>,int,int> > kmers_plus;
     basevector B;
     basevector G;
     vecbasevector blob;
     vec<Bool> to_delete;
     vec<int> path;

};

void RefAlign(
     // inputs:
     const vecbasevector& genome, const int HBK, const vecbasevector& tigs,
     const vec<int>& inv, const digraphE<vec<int>>& D, 
     const vec<int>& dinv, const vec<vec<vec<vec<int>>>>& dlines,
     const MasterVec< SerfVec<triple<int,int,int> > >& alignsb,
     // output:
     MasterVec<SerfVec<refalign>>& galigns,
     // control:
     const int verbosity,
     const vec<int>& test_edges,
     const Bool allow_second_best = False
     );

template<int K> void RefAlignCore(

     // Edge to be aligned.

     const vec<int>& X,

     // Assembly info.

     const int HBK, 
     const vecbasevector& tigs,
     const vec<int>& inv,
     const digraphE<vec<int>>& D,
     const vec<int>& dinv,
     const vec<vec<vec<vec<int>>>>& dlines,
     const vecbasevector& genome,
     const MasterVec< SerfVec<triple<int,int,int> > >& alignsb,

     // Output and output control.

     SerfVec<refalign>& galigns,
     const int verbosity,
     String& report,

     // Scratch data carried around to reduce memory allocation overhead.

     RefAlignScratch<K>& s,

     // Use only this g.

     const int g_use = -1,

     // More control.
     
     const Bool allow_second_best = False,
     const Bool optb = False

);

#endif
