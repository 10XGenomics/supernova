///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// author: Filipe Ribeiro 2012

#ifndef KMERS__NAIF_KMER__KMER_FUNCTIONS_H
#define KMERS__NAIF_KMER__KMER_FUNCTIONS_H





// ------ comparison functions ------


template <class KMER_t>
inline
bool operator<=(const KMER_t & a, const KMER_t & b) { return !(a > b); }

template <class KMER_t>
inline
bool operator>=(const KMER_t & a, const KMER_t & b) { return !(a < b); }

template <class KMER_t>
inline
bool operator!=(const KMER_t & a, const KMER_t & b) { return !(a == b); }




// ---- printing functions ----


inline
String stringify_base(const unsigned base, const char * ascii = "ACGT")
{
  ForceAssert(base < 4u);
  return ToString(ascii[base]);
}


template<class KMER_t>
inline 
String stringify_kmer(const KMER_t & kmer, const char * ascii = "ACGT")
{
  String s = "";
  const unsigned K = kmer.size();
  for (unsigned i = 0; i < K; i++)
    s += stringify_base(kmer[i], ascii);
  return s;
}


inline
String hieroglyph(const unsigned base) { return stringify_base(base, "^(-."); }

template<class KMER_t>
inline
String hieroglyphs(const KMER_t & kmer) { return stringify_kmer(kmer, "^(-."); }







// ---- conversion functions


template<class KMER_t>
inline 
KMER_t reverse_complement(const KMER_t & kmerFW) 
{
  const unsigned K = kmerFW.size();
  KMER_t kmerRC = kmerFW;          // copy, because KMER_t might have other stuff 
  for (unsigned i = 0; i != K; i++) 
    kmerRC.set(K - i - 1, 3u ^ kmerFW[i]);   // set bases from the begining
  
  return kmerRC;
}



template<class KMER_t>
inline
KMER_t canonical(const KMER_t & kmerFW) 
{
  const unsigned K = kmerFW.K();
  const KMER_t kmerRC = reverse_complement(kmerFW); // copy, because KMER_t might have other stuff 
  return (kmerFW < kmerRC) ? kmerFW : kmerRC;
}




// ---- other general utility functions ----



template<class KMER_t>
vec<unsigned> acgt_content(const KMER_t & kmer)
{
  vec<unsigned> result(4, 0);
  const unsigned K = kmer.size();
  for (unsigned i = 0; i < K; i++)
    result[kmer[i]]++;
  return result;
}







#endif
