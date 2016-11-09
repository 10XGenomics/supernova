/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "feudal/BaseVec.h"
#include "kmers/KmerRecord.h"
#include "Vec.h"

int const kmer_with_count_base::max_count;

template<int K> void kmer<K>::SetToSubOf( const basevector& source, 
     const size_type start )
{   size_t len = K;
    AssertLe( start, source.size() );
    AssertLe( len, source.size()-start );
    size_t end = len & ~15;
    int32_t* dst = (int32_t*) &data_;
    for ( size_t idx = 0; idx < end; idx += 16)
         *dst++ = source.extractKmer(start+idx, 16);
    if ( end < len ) *dst = source.extractKmer(start+end, len-end);   }

template void kmer<4>::SetToSubOf( const basevector&, const size_type );
template void kmer<6>::SetToSubOf( const basevector&, const size_type );
template void kmer<8>::SetToSubOf( const basevector&, const size_type );
template void kmer<12>::SetToSubOf( const basevector&, const size_type );
template void kmer<16>::SetToSubOf( const basevector&, const size_type );
template void kmer<20>::SetToSubOf( const basevector&, const size_type );
template void kmer<24>::SetToSubOf( const basevector&, const size_type );
template void kmer<28>::SetToSubOf( const basevector&, const size_type );
template void kmer<32>::SetToSubOf( const basevector&, const size_type );
template void kmer<36>::SetToSubOf( const basevector&, const size_type );
template void kmer<40>::SetToSubOf( const basevector&, const size_type );
template void kmer<48>::SetToSubOf( const basevector&, const size_type );
template void kmer<60>::SetToSubOf( const basevector&, const size_type );
template void kmer<61>::SetToSubOf( const basevector&, const size_type );
template void kmer<64>::SetToSubOf( const basevector&, const size_type );
template void kmer<72>::SetToSubOf( const basevector&, const size_type );
template void kmer<80>::SetToSubOf( const basevector&, const size_type );
template void kmer<84>::SetToSubOf( const basevector&, const size_type );
template void kmer<88>::SetToSubOf( const basevector&, const size_type );
template void kmer<96>::SetToSubOf( const basevector&, const size_type );
template void kmer<100>::SetToSubOf( const basevector&, const size_type );
template void kmer<128>::SetToSubOf( const basevector&, const size_type );
template void kmer<144>::SetToSubOf( const basevector&, const size_type );
template void kmer<172>::SetToSubOf( const basevector&, const size_type );
template void kmer<192>::SetToSubOf( const basevector&, const size_type );
template void kmer<200>::SetToSubOf( const basevector&, const size_type );
template void kmer<320>::SetToSubOf( const basevector&, const size_type );
template void kmer<368>::SetToSubOf( const basevector&, const size_type );
template void kmer<400>::SetToSubOf( const basevector&, const size_type );
template void kmer<440>::SetToSubOf( const basevector&, const size_type );
template void kmer<460>::SetToSubOf( const basevector&, const size_type );
template void kmer<500>::SetToSubOf( const basevector&, const size_type );
template void kmer<544>::SetToSubOf( const basevector&, const size_type );
template void kmer<640>::SetToSubOf( const basevector&, const size_type );
template void kmer<720>::SetToSubOf( const basevector&, const size_type );
template void kmer<1000>::SetToSubOf( const basevector&, const size_type );
template void kmer<1200>::SetToSubOf( const basevector&, const size_type );
template void kmer<1600>::SetToSubOf( const basevector&, const size_type );
template void kmer<2000>::SetToSubOf( const basevector&, const size_type );
template void kmer<4000>::SetToSubOf( const basevector&, const size_type );
template void kmer<10000>::SetToSubOf( const basevector&, const size_type );
