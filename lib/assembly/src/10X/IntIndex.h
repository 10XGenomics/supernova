// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

// Inverter.  Effectively, transpose a vec<vec<int>>.
// However, applied to a ReadPathVec.

#ifndef TENX_INT_INDEX_H
#define TENX_INT_INDEX_H

#include "CoreTools.h"
#include "paths/long/ReadPath.h"

// An IntIndex is a feudal-type data structure that can dynamically manage
// integers of size four or eight bytes, depending on how big they are.

class IntIndex {

     public:

     IntIndex( const ReadPathVec& paths, const int NE, const Bool verbose = False );
     void Initialize( 
          const ReadPathVec& paths, const int NE, const Bool verbose = False );

     int64_t N( ) const { return index_.size( ) - 1; }

     int64_t Count( const int e ) const
     {    return index_[e+1] - index_[e];    }

     int64_t Val( const int e, const int i ) const
     {    if ( !big_ ) return core_[ index_[e] + i ];
          else return core_big_[ index_[e] + i ];    }

     // Need:

     // 1. read all
     // 2. write all
     // 3. read selected
     // 4. virtual read

     private:

     vec<uint32_t> core_;
     vec<uint64_t> core_big_;
     vec<int64_t> index_;
     Bool big_;

};

#endif
