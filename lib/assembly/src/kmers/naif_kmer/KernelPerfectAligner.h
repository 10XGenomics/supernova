///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////



#ifndef KMERS__NAIF_KMER__KERNEL_PERFECT_ALIGNER_H
#define KMERS__NAIF_KMER__KERNEL_PERFECT_ALIGNER_H



#include "system/LockedData.h"
#include "kmers/naif_kmer/Kmers.h"
#include "kmers/naif_kmer/NaifKmerizer.h"
#include "ParallelVecUtilities.h"
#include "Vec.h"

typedef vec<triple<int64_t, int64_t, int> > PerfectAligns;



// -------- REC_t class --------

template<class KMER_t>
class KmerAndLocPA : public KMER_t
{
public:
  uint64_t i_read         : 32;   // up to 2^32 =  4G
  uint64_t i_base         : 30;   // up to 2^30 =  1G 
  uint64_t is_rc          :  1;   // for kmers, whether it is RC of FW
  uint64_t is_palindrome  :  1;   // for kmers, whether it is a palindrome

  explicit KmerAndLocPA(const unsigned K = 0) :
    KMER_t(K),
    i_read(0),
    i_base(0),
    is_rc(0),
    is_palindrome(0)
  {}

};



// -------- KERNEL_t class --------

template<class KMER_t>
class KernelPerfectAligner
{
private:
  static const unsigned tag_A = 0; // these tags, upon sort, put the set A before set B
  static const unsigned tag_B = 1;

  const BaseVecVec & _bases;      // all bv's, first from set A followed by set B
  const size_t       _nbvA;       // number of bv's in set A
  const size_t       _K;
  const size_t       _max_placements;

  PerfectAligns      _pa_tmp;     // only used in clones
  PerfectAligns    * _p_pa;       // final results go here; only merge() adds to it

  const bool         _ignore_palindromes;       // whether to skip palindromes
  const bool         _ignore_rel_rc;            // only fw<->fw or rc<->rc

  LockedData         _lock;       // lock for merge()


public:
  typedef KmerAndLocPA<KMER_t>   rec_type;

  KernelPerfectAligner(const BaseVecVec & bases, 
                       const size_t       nbvA, 
                       const size_t       K,
                       const size_t       max_placements,
                       PerfectAligns    * p_pa,
                       const bool         ignore_palindromes = false,
                       const bool         ignore_rel_rc = false
                       ):
    _bases(bases), 
    _nbvA(nbvA), 
    _K(K), 
    _max_placements(max_placements),
    _pa_tmp(),
    _p_pa(p_pa), 
    _ignore_palindromes( ignore_palindromes ),
    _ignore_rel_rc( ignore_rel_rc ),
    _lock() 
  {}

  // copy constructor for temporary kernels
  explicit KernelPerfectAligner(const KernelPerfectAligner & that) :
    _bases(that._bases),
    _nbvA(that._nbvA), 
    _K(that._K),
    _max_placements(that._max_placements),
    _pa_tmp(),
    _p_pa(0),
    _ignore_palindromes( that._ignore_palindromes ),
    _ignore_rel_rc( that._ignore_rel_rc ),
    _lock()
  {}

  // ==== interface function needed by naif_kmerize() ====
  size_t K() const { return _K; }

  // ==== interface function needed by naif_kmerize() ====
  const BaseVecVec & bases() const { return _bases; }




  // ==== interface function needed by naif_kmerize() ====
  void parse_base_vec(ParcelBuffer<rec_type> * p_buf, 
                      const size_t ibv)
  {
    ParcelBuffer<rec_type> & parcels = *p_buf;
    SubKmers<BaseVec, rec_type> kmer_cur(_K, _bases[ibv]);
    unsigned const tag = (ibv < _nbvA) ?  tag_A : tag_B; 
    
    if (tag == tag_A) { // store only first kmer for the read set A
      if (kmer_cur.not_done()) {  // checks if there are any kmers in read
        rec_type & kmer = kmer_cur.canonical();
        if (parcels.in_one_parcel(kmer)) {
          kmer.tag_set(tag); 
          kmer.i_read        = ibv;
          kmer.i_base        = kmer_cur.index_start();
          kmer.is_rc         = kmer_cur.is_canonical_rc();
          kmer.is_palindrome = kmer_cur.is_palindrome();
          parcels.add(kmer);
          if (0 && ibv == 34056714) {
            cout << "debug: fw " << hieroglyphs(kmer_cur.fw()) << endl;
            cout << "debug: bw " << hieroglyphs(kmer_cur.rc()) << endl;
            cout << "debug: palidrome " << kmer.is_palindrome << endl;
          }
        }
      }
    }
    else {  // (tag == tag_B) store all kmers from read set B
      while (kmer_cur.not_done()) {
        rec_type & kmer = kmer_cur.canonical();
        if (parcels.in_one_parcel(kmer)) {
          kmer.tag_set(tag); 
          kmer.i_read        = ibv;
          kmer.i_base        = kmer_cur.index_start();
          kmer.is_rc         = kmer_cur.is_canonical_rc();
          kmer.is_palindrome = kmer_cur.is_palindrome();
          parcels.add(kmer);
        }
        kmer_cur.next();          
      }
    }
  }



  // ==== interface function needed by naif_kmerize() ====
  // 
  // here the kmer alignments are extended
  //
  void summarize(const vec<rec_type> & kmer_recs, 
                 const size_t i0,
                 const size_t i1)
  {
    size_t nA = 0; // set A only the 1st kmer was stored and the sort puts them first
    size_t i = i0;
    while (i != i1 && kmer_recs[i].tag() == tag_A) {
      nA++;
      i++;
    }
    const size_t nB = i1 - i0 - nA;

    if (nA > 0 && nB > 0) {
      
      const bool is_palindrome = kmer_recs[i0].is_palindrome;
	
      for (size_t jA = 0; jA != nA &&
           ( !_ignore_palindromes || !is_palindrome ); jA++) {

	const size_t       iA = i0 + jA;
	const size_t  i_readA =   kmer_recs[iA].i_read;
	const bool     is_fwA = ! kmer_recs[iA].is_rc;
	const BaseVec & readA = _bases[i_readA];
	const size_t      nbA = readA.size();


	// ---- read same length as K, found perfect alignment

	if (_bases[i_readA].size() == _K) {  
	  
	  if ((!is_palindrome && _max_placements >=     nB) || 
              ( is_palindrome && _max_placements >= 2 * nB) || 
              _max_placements == 0) {
	    
            for (size_t iB = i0 + nA; iB != i1; iB++) {
	      const size_t   i_readB =   kmer_recs[iB].i_read;
	      const size_t   i_baseB =   kmer_recs[iB].i_base;
	      const bool      is_fwB = ! kmer_recs[iB].is_rc;
	      const int         posB = (is_fwB == is_fwA) ? i_baseB : - signed(i_baseB) - 1;
	      if ( !_ignore_rel_rc || is_fwB == is_fwA ) {
                   _pa_tmp.push(i_readA, i_readB - _nbvA,  posB);
                   if (is_palindrome)
                     _pa_tmp.push(i_readA, i_readB - _nbvA,  -posB - 1);
	      }
	    }
	  }
        }

        // ---- read longer than K, must extend to check perfect alignment
      
	else { 

	  PerfectAligns pa_local;
          pa_local.reserve(nB);
          
	  for (size_t iB = i0 + nA; iB != i1; iB++) {
	    if (_max_placements >= pa_local.size() || 
		_max_placements == 0) { 

	      const size_t   i_readB =   kmer_recs[iB].i_read;
	      const size_t   i_baseB =   kmer_recs[iB].i_base;
	      const bool      is_fwB = ! kmer_recs[iB].is_rc;
	      const BaseVec  & readB = _bases[i_readB];
	      const size_t       nbB = readB.size();
	      
	      if (is_fwB == is_fwA || is_palindrome) { // both FW or both RC alignment or palindrome
		if (i_baseB + nbA <= nbB) {  // read A fully inside read B
		  bool matches = true;
		  size_t ibA = _K;
                  size_t ibB = i_baseB + _K;
		  while (matches && ibA != nbA) {
		    if (readA[ibA] != readB[ibB]) 
		      matches = false;
		    ibA++;
                    ibB++;
		  }
		  if (matches) 
		    pa_local.push(i_readA, i_readB - _nbvA,  int(i_baseB));
		}
	      }
	      if ((!_ignore_rel_rc && is_fwB != is_fwA) || is_palindrome) { // one FW and one RC alignment or palindrome
		
		if (nbA <= i_baseB + _K) { // read A fully inside read B
		  bool matches = true;
		  size_t ibA = _K;
                  size_t ibB = i_baseB - 1;
		  while (matches && ibA != nbA) {
		    if (readA[ibA] != 3 - readB[ibB]) 
		      matches = false;
		    ibA++;
                    ibB--;
		  }
		  if (matches) 
		    pa_local.push(i_readA, i_readB - _nbvA,  -int(i_baseB + _K - nbA) - 1);
		}
	      }

	    }  
	  } // for (size_t iB = i0 + nA; iB != i1; iB++)
	  
	  if (_max_placements >= pa_local.size() || 
	      _max_placements == 0) 
	    _pa_tmp.insert(_pa_tmp.end(), pa_local.begin(), pa_local.end());
	}
        
	
      } // for (size_t jA = 0; jA != nA; jA++)
      
    } // if (nA > 0 && nB > 0)
  
  }
  

  
  // ==== interface function needed by naif_kmerize() ====
  // 
  //  only called on the main kernel to merge the results of the temp kernels 
  // 
  void merge(const KernelPerfectAligner & kpa_tmp,
	     const size_t i_parcel)
  {
    Locker lock(_lock);
    
    _p_pa->insert(_p_pa->end(), kpa_tmp._pa_tmp.begin(), kpa_tmp._pa_tmp.end());
  }

};

template <class KMER_t> unsigned const KernelPerfectAligner<KMER_t>::tag_A;
template <class KMER_t> unsigned const KernelPerfectAligner<KMER_t>::tag_B;


// -------- KERNEL_t class --------
typedef vec<triple<int64_t, int, int64_t> > KmerAligns;

template<class KMER_t>
class KernelAllKmerAligns
{
private:
  static const unsigned tag_A = 0; // these tags, upon sort, put the set A before set B
  static const unsigned tag_B = 1;

  const BaseVecVec & _bases;      // all bv's, first from set A followed by set B
  const size_t       _nbvA;       // number of bv's in set A
  const size_t       _K;

  KmerAligns      _pa_tmp;     // only used in clones
  KmerAligns    * _p_pa;       // final results go here; only merge() adds to it

  const bool         _ignore_palindromes;       // whether to skip palindromes
  const bool         _ignore_rel_rc;            // only fw<->fw or rc<->rc

  LockedData         _lock;       // lock for merge()


public:
  typedef KmerAndLocPA<KMER_t>   rec_type;

  KernelAllKmerAligns(const BaseVecVec & bases,
                       const size_t       nbvA,
                       const size_t       K,
                       KmerAligns    * p_pa,
                       const bool         ignore_palindromes = false,
                       const bool         ignore_rel_rc = false
                       ):
    _bases(bases),
    _nbvA(nbvA),
    _K(K),
    _pa_tmp(),
    _p_pa(p_pa),
    _ignore_palindromes( ignore_palindromes ),
    _ignore_rel_rc( ignore_rel_rc ),
    _lock()
  {}

  // copy constructor for temporary kernels
  explicit KernelAllKmerAligns(const KernelAllKmerAligns & that) :
    _bases(that._bases),
    _nbvA(that._nbvA),
    _K(that._K),
    _pa_tmp(),
    _p_pa(0),
    _ignore_palindromes( that._ignore_palindromes ),
    _ignore_rel_rc( that._ignore_rel_rc ),
    _lock()
  {}

  // ==== interface function needed by naif_kmerize() ====
  size_t K() const { return _K; }

  // ==== interface function needed by naif_kmerize() ====
  const BaseVecVec & bases() const { return _bases; }




  // ==== interface function needed by naif_kmerize() ====
  void parse_base_vec(ParcelBuffer<rec_type> * p_buf,
                      const size_t ibv)
  {
    ParcelBuffer<rec_type> & parcels = *p_buf;
    SubKmers<BaseVec, rec_type> kmer_cur(_K, _bases[ibv]);
    unsigned const tag = (ibv < _nbvA) ?  tag_A : tag_B;

    while (kmer_cur.not_done()) {
         rec_type & kmer = kmer_cur.canonical();
         if (parcels.in_one_parcel(kmer)) {
              kmer.tag_set(tag);
              kmer.i_read        = ibv;
              kmer.i_base        = kmer_cur.index_start();
              kmer.is_rc         = kmer_cur.is_canonical_rc();
              kmer.is_palindrome = kmer_cur.is_palindrome();
              parcels.add(kmer);
         }
         kmer_cur.next();
    }
  }



  // ==== interface function needed by naif_kmerize() ====
  //
  // here the kmer alignments are extended
  //
  void summarize(const vec<rec_type> & kmer_recs,
                 const size_t i0,
                 const size_t i1)
  {
    size_t nA = 0; // set A only the 1st kmer was stored and the sort puts them first
    size_t i = i0;
    while (i != i1 && kmer_recs[i].tag() == tag_A) {
      nA++;
      i++;
    }
    const size_t nB = i1 - i0 - nA;

    if (nA > 0 && nB > 0) {

         const bool is_palindrome = kmer_recs[i0].is_palindrome;

         for (size_t jA = 0; jA != nA &&
              ( !_ignore_palindromes || !is_palindrome ); jA++) {

              const size_t       iA = i0 + jA;
              const size_t  i_readA =   kmer_recs[iA].i_read;
              const size_t  i_baseA =   kmer_recs[iA].i_base;
              const bool     is_fwA = ! kmer_recs[iA].is_rc;
              const BaseVec & readA = _bases[i_readA];
              const size_t      nbA = readA.size();
              const int        posA = i_baseA;


                   for (size_t iB = i0 + nA; iB != i1; iB++) {
                        const size_t   i_readB =   kmer_recs[iB].i_read;
                        const size_t   i_baseB =   kmer_recs[iB].i_base;
                        const bool      is_fwB = ! kmer_recs[iB].is_rc;
                        const int         posB = (is_fwB == is_fwA) ? i_baseB : - signed(i_baseB) - 1;
                        if ( !_ignore_rel_rc || is_fwB == is_fwA ) {
                             _pa_tmp.push(i_readA, posA - posB, i_readB - _nbvA);
                             if (is_palindrome)
                                  _pa_tmp.push(i_readA,  posA - (-posB - 1), i_readB - _nbvA);
                        }
                   }
              }

    } // if (nA > 0 && nB > 0)

  }


  // ==== interface function needed by naif_kmerize() ====
  //
  //  only called on the main kernel to merge the results of the temp kernels
  //
  void merge(KernelAllKmerAligns & kpa_tmp,
	     const size_t i_parcel)
  {

    KmerAligns& tmp = kpa_tmp._pa_tmp;

    std::sort(tmp.begin(),tmp.end());
    tmp.resize( std::unique( tmp.begin(), tmp.end() ) - tmp.begin());

    Locker lock(_lock);

    size_t norig = _p_pa->size();

    _p_pa->insert(_p_pa->end(), tmp.begin(), tmp.end() );
    //    std::inplace_merge(_p_pa->begin(), _p_pa->begin() + norig, _p_pa->end() );
    //    _p_pa->resize(std::unique( _p_pa->begin(), _p_pa->end() ) - _p_pa->begin() );


//    _p_pa->insert(_p_pa->end(), kpa_tmp._pa_tmp.begin(), kpa_tmp._pa_tmp.end());

    // _p_pa->insert(_p_pa->end(), tmp.begin(), tmp.end() );
    // std::sort( _p_pa->begin() + norig, _p_pa->end() );
    // auto end = std::unique( _p_pa->begin() + norig, _p_pa->end() );
    // _p_pa->resize( end - _p_pa->begin() );
    // std::inplace_merge(_p_pa->begin(), _p_pa->begin()+norig, _p_pa->end() );
//    auto end = std::unique( _p_pa->begin(), _p_pa->end() );
//    _p_pa->resize( end - _p_pa->begin() );
  }

};


template <class KMER_t> unsigned const KernelAllKmerAligns<KMER_t>::tag_A;
template <class KMER_t> unsigned const KernelAllKmerAligns<KMER_t>::tag_B;








/*


template<class KMER_t>
class KernelPerfectAligner2
{
public:
  typedef KmerAndLocPA<KMER_t>   rec_type;

private:
  const BaseVecVec & _bvsA;      // the queries
  const BaseVecVec & _bvsB;      // the targets
  const size_t       _K;
  const size_t       _n_threads;
  const size_t       _max_placements;
  
  vec<rec_type>      _kmer_recs;
  vec<rec_type>    * _p_kmer_recs_tmp;
  PerfectAligns    * _p_pa;       // final results go here; only process_merged() adds to it
  LockedData         _lock;       // lock for merge()


public:

  KernelPerfectAligner2(const BaseVecVec & bvsA, 
			const BaseVecVec & bvsB, 
			const size_t       K,
			const size_t       n_threads,
			const size_t       max_placements,
			PerfectAligns    * p_pa) :
    _bvsA(bvsA), 
    _bvsB(bvsB), 
    _K(K), 
    _n_threads(n_threads),
    _max_placements(max_placements),
    _kmer_recs(0),
    _p_pa(p_pa), 
    _lock() 
  {}

  // copy constructor for temporary kernels
  explicit KernelPerfectAligner2(const KernelPerfectAligner2 & that) :
    _bvsA(that._bvsA),
    _bvsB(that._bvsB),
    _K(that._K),
    _n_threads(that._n_threads),
    _max_placements(that._max_placements),
    _kmer_recs(0),
    _p_pa(0),
    _lock()
  {}

  // ==== interface function needed by naif_kmerize() ====
  size_t K() const { return _K; }

  // ==== interface function needed by naif_kmerize() ====
  const BaseVecVec & bases() const { return _bvsB; }



 
  // ==== interface function needed by naif_kmerize() ====
  void parse_base_vec(ParcelBuffer<rec_type> * p_buf, 
                      const size_t ibv)
  {
    ParcelBuffer<rec_type> & parcels = *p_buf;
    SubKmers<BaseVec, rec_type> kmer_cur(_K, _bvsB[ibv]);
    while (kmer_cur.not_done()) {
      rec_type & kmer = kmer_cur.canonical();
      if (parcels.in_one_parcel(kmer)) {
        kmer.tag_set(tag); 
        kmer.i_read        = ibv;
        kmer.i_base        = kmer_cur.index_start();
        kmer.is_rc         = kmer_cur.is_canonical_rc();
        kmer.is_palindrome = kmer_cur.is_palindrome();
        parcels.add(kmer);
      }
      kmer_cur.next();          
    }
  }


  // ==== interface function needed by naif_kmerize() ====
  // 
  void summarize(const vec<rec_type> & kmer_recs, 
                 const size_t i0,
                 const size_t i1)
  {
    // just need to store the pointer to the outside kmer recs vector
    if (i0 == 0) _p_kmer_recs_tmp = & kmer_recs;
  }

  

  // ==== interface function needed by naif_kmerize() ====
  // 
  //  only called on the main kernel to merge the results of the temp kernels 
  // 
  void merge(const KernelPerfectAligner2 & kernel_tmp,
	     const size_t i_parcel)
  {
    // merge all the parcels for this pass into one
    Locker lock(_lock);
    _kmer_recs.insert(_kmer_recs.begin(), 
                      _p_kmer_recs_tmp->begin(),
                      _p_kmer_recs_tmp->end());
  }


  // ==== interface function needed by naif_kmerize() ====
  // 
  //  only called on the main kernel to process merged results
  // 
  void process_merged(const vec<rec_type> & parcel,
                      const size_t i_parcel)
  {
    {


    const size_t nbvA = _bvsA.size();
    const size_t i0bvA = (nbvA *  i_parcel     ) / _n_threads;
    const size_t i1bvA = (nbvA * (i_parcel + 1)) / _n_threads;

    for (size_t ibv = i0bvA; ibv != i1bvA; ibv++) {

      const size_t i0 = lower_bound();
    }

      const size_t nB = i1 - i0 - nA;

    if (nA > 0 && nB > 0) {
      
      const bool is_palindrome = kmer_recs[i0].is_palindrome;
	
      for (size_t jA = 0; jA != nA; jA++) {

	const size_t       iA = i0 + jA;
	const size_t  i_readA =   kmer_recs[iA].i_read;
	const bool     is_fwA = ! kmer_recs[iA].is_rc;
	const BaseVec & readA = _bases[i_readA];
	const size_t      nbA = readA.size();


	// ---- read same length as K, found perfect alignment

	if (_bases[i_readA].size() == _K) {  
	  
	  if ((!is_palindrome && _max_placements >=     nB) || 
              ( is_palindrome && _max_placements >= 2 * nB) || 
              _max_placements == 0) {
	    
            for (size_t iB = i0 + nA; iB != i1; iB++) {
	      const size_t   i_readB =   kmer_recs[iB].i_read;
	      const size_t   i_baseB =   kmer_recs[iB].i_base;
	      const bool      is_fwB = ! kmer_recs[iB].is_rc;
	      const int         posB = (is_fwB == is_fwA) ? i_baseB : - signed(i_baseB) - 1;
	      _pa_tmp.push(i_readA, i_readB - _nbvA,  posB);
	      if (is_palindrome)
		_pa_tmp.push(i_readA, i_readB - _nbvA,  -posB - 1);
	    }
	  }
        }

        // ---- read longer than K, must extend to check perfect alignment
      
	else { 

	  PerfectAligns pa_local;
          pa_local.reserve(nB);
          
	  for (size_t iB = i0 + nA; iB != i1; iB++) {
	    if (_max_placements >= pa_local.size() || 
		_max_placements == 0) { 

	      const size_t   i_readB =   kmer_recs[iB].i_read;
	      const size_t   i_baseB =   kmer_recs[iB].i_base;
	      const bool      is_fwB = ! kmer_recs[iB].is_rc;
	      const BaseVec  & readB = _bases[i_readB];
	      const size_t       nbB = readB.size();
	      
	      if (is_fwB == is_fwA || is_palindrome) { // both FW or both RC alignment or palindrome
		if (i_baseB + nbA <= nbB) {  // read A fully inside read B
		  bool matches = true;
		  size_t ibA = _K;
                  size_t ibB = i_baseB + _K;
		  while (matches && ibA != nbA) {
		    if (readA[ibA] != readB[ibB]) 
		      matches = false;
		    ibA++;
                    ibB++;
		  }
		  if (matches) 
		    pa_local.push(i_readA, i_readB - _nbvA,  int(i_baseB));
		}
	      }
	      if (is_fwB != is_fwA || is_palindrome) { // one FW and one RC alignment or palindrome
		
		if (nbA <= i_baseB + _K) { // read A fully inside read B
		  bool matches = true;
		  size_t ibA = _K;
                  size_t ibB = i_baseB - 1;
		  while (matches && ibA != nbA) {
		    if (readA[ibA] != 3 - readB[ibB]) 
		      matches = false;
		    ibA++;
                    ibB--;
		  }
		  if (matches) 
		    pa_local.push(i_readA, i_readB - _nbvA,  -int(i_baseB + _K - nbA) - 1);
		}
	      }

	    }  
	  } // for (size_t iB = i0 + nA; iB != i1; iB++)
	  
	  if (_max_placements >= pa_local.size() || 
	      _max_placements == 0) 
	    _pa_tmp.insert(_pa_tmp.end(), pa_local.begin(), pa_local.end());
	}
        
	
      } // for (size_t jA = 0; jA != nA; jA++)
      
    } // if (nA > 0 && nB > 0)
  
  }
  








  }





};

*/



#endif
