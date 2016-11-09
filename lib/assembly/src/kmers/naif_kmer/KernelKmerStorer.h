///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

//   
//  A kmer storer kernel generates a vector of kmer records.
//
//  Each record contains the kmer and at least the frequency. 
//
//
//

#ifndef KMERS__NAIF_KMER__KERNEL_KMER_STORER_H
#define KMERS__NAIF_KMER__KERNEL_KMER_STORER_H

#include "kmers/naif_kmer/Kmers.h"
#include "kmers/naif_kmer/NaifKmerizer.h"

// -------- KERNEL_t class --------
//
//  ELEM_t must provide the method     ELEM_t::set_freq(freq)
//  

template<class ELEM_t>
class KernelKmerStorer
{
private:
  const BaseVecVec   & _bvv;
  const size_t         _K;

  vec<ELEM_t>          _kvec_tmp; // only used in temporaries
  vec<ELEM_t>        * _p_kvec;     // final results go here; only merge() adds to it
  LockedData           _lock;     // lock for merge()

  const Validator    * _p_validator;
  
public:
  typedef typename ELEM_t::kmer_type  rec_type;

  KernelKmerStorer(const BaseVecVec & bvv, 
                   const size_t       K,
                   vec<ELEM_t>      * p_kvec,
                   const Validator  * p_validator = 0) :
    _bvv(bvv), 
    _K(K), 
    _kvec_tmp(), 
    _p_kvec(p_kvec), 
    _lock(),
    _p_validator(p_validator)
  {}

  // copy constructor for temporary kernels
  explicit KernelKmerStorer(const KernelKmerStorer<ELEM_t> & that) :
    _bvv(that._bvv),
    _K(that._K),
    _kvec_tmp(),
    _p_kvec(0),
    _lock(),
    _p_validator(that._p_validator)
  {}

  // interface function needed by naif_kmerize()
  size_t K() const { return _K; }

  // interface function needed by naif_kmerize()
  const BaseVecVec & bases() const { return _bvv; }
 

  // interface function needed by naif_kmerize()
  void parse_base_vec(ParcelBuffer<rec_type> * p_buf, 
                      const size_t ibv)
  {
    ParcelBuffer<rec_type> & parcel_buf = *p_buf;
    SubKmers<BaseVec, rec_type> kmer_cur(_K, _bvv[ibv]);
    while (kmer_cur.not_done()) {
      rec_type & kmer = kmer_cur.canonical();
      if (parcel_buf.in_one_parcel(kmer))
        parcel_buf.add(kmer);
      kmer_cur.next();          
    }
  }


  // interface function needed by naif_kmerize()
  void summarize(const vec<rec_type> & krecs,  
                 const size_t i0k,
                 const size_t i1k)
  {
    const size_t kf = i1k - i0k;

    if (!_p_validator || (*_p_validator)(kf)) {
      ELEM_t krec(krecs[i0k]);
      krec.set_freq(kf);
      _kvec_tmp.push_back(krec);
    }
  }
 

  // interface function needed by naif_kmerize()
  void merge(const KernelKmerStorer<ELEM_t> & kernel_tmp,
	     const size_t i_parcel)
  {
    Locker lock(_lock);
    _p_kvec->insert(_p_kvec->end(), kernel_tmp._kvec_tmp.begin(), kernel_tmp._kvec_tmp.end());
  }

};




// -------- KERNEL_t class --------
//
//  ELEM_t must provide the method     TBA -- methods for adding
//  locations
//

template<class ELEM_t>
class KernelKmerBVLocStorer
{
private:
  const BaseVecVec   & _bvv;
  const size_t         _K;

  vec<ELEM_t>          _kvec_tmp; // only used in temporaries
  vec<ELEM_t>        * _p_kvec;     // final results go here; only merge() adds to it
  LockedData           _lock;     // lock for merge()

  const Validator    * _p_validator;

  typedef typename ELEM_t::kmer_type kmer_type;
public:
  typedef KmerBVLoc<kmer_type> rec_type;

  KernelKmerBVLocStorer(const BaseVecVec & bvv,
                   const size_t       K,
                   vec<ELEM_t>      * p_kvec,
                   const Validator  * p_validator = 0) :
    _bvv(bvv),
    _K(K),
    _kvec_tmp(),
    _p_kvec(p_kvec),
    _lock(),
    _p_validator(p_validator)
  {}

  // copy constructor for temporary kernels
  explicit KernelKmerBVLocStorer(const KernelKmerBVLocStorer<ELEM_t> & that) :
    _bvv(that._bvv),
    _K(that._K),
    _kvec_tmp(),
    _p_kvec(0),
    _lock(),
    _p_validator(that._p_validator)
  {}

  // interface function needed by naif_kmerize()
  size_t K() const { return _K; }

  // interface function needed by naif_kmerize()
  const BaseVecVec & bases() const { return _bvv; }


  // interface function needed by naif_kmerize()
  void parse_base_vec(ParcelBuffer<rec_type> * p_buf,
                      const size_t ibv)
  {
    ParcelBuffer<rec_type> & parcel_buf = *p_buf;
    SubKmers<BaseVec, rec_type> kmer_cur(_K, _bvv[ibv]);
    while (kmer_cur.not_done()) {
      rec_type & kmer = kmer_cur.canonical();
      if (parcel_buf.in_one_parcel(kmer)) {
        kmer.set_ibv(ibv);
        kmer.set_ib(kmer_cur.index_start());

        if (kmer_cur.is_palindrome())  kmer.set_palindrome(true);
        else                           kmer.set_fw(kmer_cur.is_canonical_fw());

        parcel_buf.add(kmer);
      }
      kmer_cur.next();
    }
  }


  // interface function needed by naif_kmerize()
  void summarize(const vec<rec_type> & krecs,
                 const size_t i0k,
                 const size_t i1k)
  {
    const size_t kf = i1k - i0k;

    if (!_p_validator || (*_p_validator)(kf)) {
      ELEM_t krec(krecs[i0k]);
      for ( size_t i = i0k; i < i1k; ++i ) {
          krec.locs().push_back( krecs[i] );     // should peel off the BVLocs 
      }
      _kvec_tmp.push_back(krec);
    }
  }


  // interface function needed by naif_kmerize()
  void merge(const KernelKmerBVLocStorer<ELEM_t> & kernel_tmp,
	     const size_t i_parcel)
  {
    Locker lock(_lock);
    _p_kvec->insert(_p_kvec->end(), kernel_tmp._kvec_tmp.begin(), kernel_tmp._kvec_tmp.end());
  }

};












// -------- KERNEL_t class --------
//
//  Stores the frequency and the affixes (prefixes and suffixes) of each kmer
//
//  ELEM_t must provide the methods
//
//     ELEM_t::set_freq(freq)
//     ELEM_t::set_prefix(base)
//     ELEM_t::set_suffix(base)
//


template<class ELEM_t>
class KernelKmerAffixesStorer
{
private:
  const BaseVecVec   & _bvv;
  const size_t         _K;
  const size_t         _Ks;

  vec<ELEM_t>        * _p_kvec;   // final results go here; only merge() adds to it
  vec<ELEM_t>          _kvec_tmp; // only used in temporaries
  LockedData           _lock;     // lock for merge()

  const Validator    * _p_validator;

  bool _is_valid_kf(const size_t kf) const 
  {
    return (kf > 0 && (!_p_validator || (*_p_validator)(kf)));
  }

  typedef typename ELEM_t::kmer_type kmer_type;
public:
  typedef KmerBVLoc<kmer_type> rec_type;

  KernelKmerAffixesStorer(const BaseVecVec & bvv, 
                          const size_t       K,
                          vec<ELEM_t>      * p_kvec,
                          const Validator  * p_validator = 0) :
    _bvv(bvv), 
    _K(K), 
    _Ks(K - 2),
    _p_kvec(p_kvec), 
    _kvec_tmp(), 
    _lock(),
    _p_validator(p_validator)
  { 
    //ForceAssert(K % 2 == 1);
  }

  // copy constructor for temporary kernels
  explicit KernelKmerAffixesStorer(const KernelKmerAffixesStorer<ELEM_t> & that) :
    _bvv(that._bvv),
    _K(that._K),
    _Ks(that._Ks),
    _p_kvec(0),
    _kvec_tmp(),
    _lock(),
    _p_validator(that._p_validator)
  {}

  // interface function needed by naif_kmerize()
  size_t K() const { return _K; }

  // interface function needed by naif_kmerize()
  const BaseVecVec & bases() const { return _bvv; }
 

  // interface function needed by naif_kmerize()
  // Here we look for (k-2)mers and store their locations
  void parse_base_vec(ParcelBuffer<rec_type> * p_parcels, 
                      const size_t ibv)
  {
    ParcelBuffer<rec_type> & parcels = *p_parcels;
    SubKmers<BaseVec, rec_type> kmer_cur(_Ks, _bvv[ibv]);
    while (kmer_cur.not_done()) {
      rec_type & kmer = kmer_cur.canonical();
      if (parcels.in_one_parcel(kmer)) {
        kmer.set_ibv(ibv);
        kmer.set_ib(kmer_cur.index_start());
        
        if (kmer_cur.is_palindrome())  kmer.set_palindrome(true);
        else                           kmer.set_fw(kmer_cur.is_canonical_fw());

        parcels.add(kmer);
      }
      kmer_cur.next();          
    }
  }


  // interface function needed by naif_kmerize()
  void summarize(const vec<rec_type> & krecs,  
                 const size_t i0k,
                 const size_t i1k)
  {
    const bool is_palindrome = krecs[i0k].is_palindrome();
    
    // compute frequencies for various kmers around central (k-2)mer

    vec<vec<size_t> > kf_pre_suf(4, vec<size_t>(4, 0));
    vec<vec<size_t> > kf_suf_suf(4, vec<size_t>(4, 0));
    vec<vec<size_t> > kf_pre_pre(4, vec<size_t>(4, 0));

    for (size_t ik = i0k; ik < i1k; ik++) {
      const rec_type & krec    = krecs[ik];
      const BaseVec  & bv      = _bvv[krec.ibv()];
      const unsigned   nb      = bv.size();
      const unsigned   ib_kmer = krec.ib();

      if (krec.is_fw()) {         // (k-2)mer has FW orientation in read

        ForceAssert(!is_palindrome); // make sure it's not a palindrome
	
	const bool good_pre2 = (ib_kmer >= 2);
	const bool good_pre1 = (ib_kmer >= 1);
	const bool good_suf1 = (ib_kmer + _Ks < nb);
	const bool good_suf2 = (ib_kmer + _Ks + 1 < nb);
	
	const unsigned pre2 = good_pre2 ? bv[ib_kmer - 2] : 0;
	const unsigned pre1 = good_pre1 ? bv[ib_kmer - 1] : 0;
	const unsigned suf1 = good_suf1 ? bv[ib_kmer + _Ks] : 0;
	const unsigned suf2 = good_suf2 ? bv[ib_kmer + _Ks + 1] : 0;
	
	if (good_pre2)              kf_pre_pre[pre2][pre1]++;
	if (good_pre1 && good_suf1) kf_pre_suf[pre1][suf1]++;
	if (good_suf2)              kf_suf_suf[suf1][suf2]++;
      }         

      if (krec.is_rc()) {         // (k-2)mer has RC orientation in read

        ForceAssert(!is_palindrome); // make sure it's not a palindrome

	const bool good_suf2 = (ib_kmer >= 2);
	const bool good_suf1 = (ib_kmer >= 1);
	const bool good_pre1 = (ib_kmer + _Ks < nb);
	const bool good_pre2 = (ib_kmer + _Ks + 1 < nb);
	
	const unsigned suf2 = good_suf2 ? (3u ^ bv[ib_kmer - 2]) : 0;
	const unsigned suf1 = good_suf1 ? (3u ^ bv[ib_kmer - 1]) : 0;
	const unsigned pre1 = good_pre1 ? (3u ^ bv[ib_kmer + _Ks]) : 0;
	const unsigned pre2 = good_pre2 ? (3u ^ bv[ib_kmer + _Ks + 1]) : 0;
	
	if (good_pre2)              kf_pre_pre[pre2][pre1]++;
	if (good_pre1 && good_suf1) kf_pre_suf[pre1][suf1]++;
	if (good_suf2)              kf_suf_suf[suf1][suf2]++;
      }         

      if (is_palindrome) {         // (k-2)mer is a palindrome

	const bool good_pre2 = (ib_kmer >= 2);
	const bool good_pre1 = (ib_kmer >= 1);
	const bool good_suf1 = (ib_kmer + _Ks < nb);
	const bool good_suf2 = (ib_kmer + _Ks + 1 < nb);
	
	const unsigned pre2_fw = good_pre2 ? bv[ib_kmer - 2] : 0;
	const unsigned pre1_fw = good_pre1 ? bv[ib_kmer - 1] : 0;
	const unsigned suf1_fw = good_suf1 ? bv[ib_kmer + _Ks] : 0;
	const unsigned suf2_fw = good_suf2 ? bv[ib_kmer + _Ks + 1] : 0;

	const unsigned pre2_rc = 3u ^ suf2_fw;
	const unsigned pre1_rc = 3u ^ suf1_fw;
	const unsigned suf1_rc = 3u ^ pre1_fw;
	const unsigned suf2_rc = 3u ^ pre2_fw;

	if (good_pre1 && good_suf1) {
	  // pick canonical k-mer
	  if (pre1_fw < pre1_rc) kf_pre_suf[pre1_fw][suf1_fw]++;
	  else                   kf_pre_suf[pre1_rc][suf1_rc]++;
	}

	if (good_pre2) {
          kf_pre_pre[pre2_fw][pre1_fw]++;
	  kf_suf_suf[suf1_rc][suf2_rc]++;
        }

	if (good_suf2) {
          kf_suf_suf[suf1_fw][suf2_fw]++;
	  kf_pre_pre[pre2_rc][pre1_rc]++;
        }
      }         

    }

    // summarize kmers and store them with number of affixes and kmer frequencies

    
    const ELEM_t kmer0(krecs[i0k]);


    for (unsigned pre1 = 0; pre1 < 4; pre1++) {
      for (unsigned suf1 = 0; suf1 < 4; suf1++) {

        const size_t & kf = kf_pre_suf[pre1][suf1];
        
        if (_is_valid_kf(kf)) {
          
          // build kmer with affix info
          
          KmerFWRC<ELEM_t> kmerFR(kmer0); // keeps both FW and RC

          kmerFR.add_left(pre1);    // add prefix
          kmerFR.add_right(suf1);   // add suffix

	  kmerFR.fw().set_freq(kf);
	  kmerFR.rc().set_freq(kf);

	  {
	    KmerFWRC<ELEM_t> kmerFRpre(kmerFR); // keeps both FW and RC
	    kmerFRpre.push_left(0u);

	    for (unsigned pre2 = 0; pre2 < 4; pre2++) {
	      kmerFRpre.set(0, pre2);
	      const bool is_pal = kmerFRpre.is_palindrome();
	      const size_t & kfpp = kf_pre_pre[pre2][pre1];
	      if (_is_valid_kf(is_pal ? kfpp/2 : kfpp)) {
		kmerFR.fw().set_prefix(pre2);
		kmerFR.rc().set_suffix(3u ^ pre2);
	      }
	    }
	  }
	  {
	    KmerFWRC<ELEM_t> kmerFRsuf(kmerFR); // keeps both FW and RC
	    kmerFRsuf.push_right(0u);
	    
	    for (unsigned suf2 = 0; suf2 < 4; suf2++) {
	      kmerFRsuf.set(_K - 1, suf2);
	      const bool is_pal = kmerFRsuf.is_palindrome();
	      const size_t & kfss = kf_suf_suf[suf1][suf2];
	      if (_is_valid_kf(is_pal ? kfss/2 : kfss)) {
		kmerFR.fw().set_suffix(suf2);
		kmerFR.rc().set_prefix(3u ^ suf2);
	      }
	    }
	  }

	  _kvec_tmp.push_back(kmerFR.canonical());
        }
      }
    }
    
  }
    

  // interface function needed by naif_kmerize()
  void merge(const KernelKmerAffixesStorer<ELEM_t> & kernel_tmp,
	     const size_t i_parcel)
  {
    Locker lock(_lock);
    _p_kvec->insert(_p_kvec->end(), kernel_tmp._kvec_tmp.begin(), kernel_tmp._kvec_tmp.end());
  }

};






// -------- KERNEL_t class --------
//
//  Stores two frequencies 
//
//  ELEM_t must provide the methods
//
//     ELEM_t::set_freq_A(freq)
//     ELEM_t::set_freq_B(freq)
//  
//

template<class ELEM_t>
class KernelKmerBiStorer
{
private:
  const BaseVecVec   & _bvv;        // all bv's, first from set A followed by set B
  const size_t         _nbvA;       // number of bv's in set A
  const size_t         _K;

  vec<ELEM_t>          _kvec_tmp; // only used in temporaries
  vec<ELEM_t>        * _p_kvec;     // final results go here; only merge() adds to it
  LockedData           _lock;     // lock for merge()

  const BiValidator  * _p_bi_validator;
  
public:
  typedef typename ELEM_t::kmer_type  rec_type;

  KernelKmerBiStorer(const BaseVecVec  & bvv, 
                     const size_t        nbvA, 
                     const size_t        K,
                     vec<ELEM_t>       * p_kvec,
                     const BiValidator * p_bi_validator = 0) :
    _bvv(bvv), 
    _nbvA(nbvA),
    _K(K), 
    _kvec_tmp(), 
    _p_kvec(p_kvec), 
    _lock(),
    _p_bi_validator(p_bi_validator)
  {}

  // copy constructor for temporary kernels
  explicit KernelKmerBiStorer(const KernelKmerBiStorer<ELEM_t> & that) :
    _bvv(that._bvv),
    _nbvA(that._nbvA), 
    _K(that._K),
    _kvec_tmp(),
    _p_kvec(0),
    _lock(),
    _p_bi_validator(that._p_bi_validator)
  {}

  // interface function needed by naif_kmerize()
  size_t K() const { return _K; }

  // interface function needed by naif_kmerize()
  const BaseVecVec & bases() const { return _bvv; }
 

  // interface function needed by naif_kmerize()
  void parse_base_vec(ParcelBuffer<rec_type> * p_buf, 
                      const size_t ibv)
  {
    ParcelBuffer<rec_type> & parcel_buf = *p_buf;
    SubKmers<BaseVec, rec_type> kmer_cur(_K, _bvv[ibv]);
    const bool is_A = (ibv < _nbvA); 

    while (kmer_cur.not_done()) {
      rec_type & kmer = kmer_cur.canonical();
      if (parcel_buf.in_one_parcel(kmer)) {
        kmer.tag_set(is_A); // tag kmer as A or B
        parcel_buf.add(kmer);
      }
      kmer_cur.next();          
    }
  }


  // interface function needed by naif_kmerize()
  void summarize(const vec<rec_type> & krecs,  
                 const size_t i0k,
                 const size_t i1k)
  {
    const size_t kf = i1k - i0k;
    
    size_t kfA = 0; 
    for (size_t ik = i0k; ik < i1k; ik++) 
      if (krecs[ik].tag())    // tag is 1 for As
        kfA++;

    const size_t kfB = kf - kfA;

    if (!_p_bi_validator || (*_p_bi_validator)(kfA, kfB)) {
      ELEM_t krec(krecs[i0k]);
      krec.set_freq_A(kfA);
      krec.set_freq_B(kfB);
      _kvec_tmp.push_back(krec);
    }
  }
 

  // interface function needed by naif_kmerize()
  void merge(const KernelKmerBiStorer<ELEM_t> & kernel_tmp,
	     const size_t i_parcel)
  {
    Locker lock(_lock);
    _p_kvec->insert(_p_kvec->end(), kernel_tmp._kvec_tmp.begin(), kernel_tmp._kvec_tmp.end());
  }

};













#endif
