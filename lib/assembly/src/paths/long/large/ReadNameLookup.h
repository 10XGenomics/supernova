///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Class readname_lookup is a class for keeping track of Illumina readnames.
//
// Read name assumptions:
//
// 1. Read name consists of colon-separated fields.
// 2. Read name suffixed by .1 or .2 to indicate first/second read of pair.
// 2. After removal of the suffix, all but one of the fields are nonnegative
//    integers.
// 3. Numeric fields together can be crammed into eight bytes.
// 4. Names for a given set all have the same layout.
// 5. At most 2^33 readnames.

#ifndef READNAME_LOOKUP_H
#define READNAME_LOOKUP_H

#include "CoreTools.h"

class readname_lookup {

     public:

     readname_lookup( ) { fcpos_ = 0; }
     readname_lookup( const vecString& names );

     uint64_t GetReadId( const String& n );

     void writeBinary( BinaryWriter& writer ) const;
     void readBinary( BinaryReader& reader );

     private:

     int fcpos_;
     vec<uint64_t> top_;
     vec<String> fcnames_;
     vec<uint64_t> keys_;
     vec<uint32_t> pids_;

     uint64_t KeyFromName( String name );

};

SELF_SERIALIZABLE(readname_lookup);

#endif
