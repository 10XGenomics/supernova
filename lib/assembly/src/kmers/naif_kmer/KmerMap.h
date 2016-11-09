///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////



#ifndef KMERS__NAIF_KMER__KMER_MAP_H
#define KMERS__NAIF_KMER__KMER_MAP_H

#include "feudal/BinaryStream.h"


// ---- this hash table contains a vec<REC_t> so it more than duplicates 
//      memory usage upon construction.
//      
//      KmerMap implemented as a chain hash table.
//
//      Deleting records cannot be implemented in an efficient way in a chain hash.
//

template<class REC_t>
class KmerMap
{
  typedef typename REC_t::kmer_type Kmer_t;

  vec<REC_t>     _hash;
  size_t         _n_rec;
  
  size_t _ih_next(size_t & ih) const
  { 
    ih++;
    if (ih >= _hash.size()) ih -= _hash.size();
    return ih;
  }

  size_t _ih0_from_kmer(const Kmer_t & kmer) const 
  { 
    return kmer.hash_64bits() % _hash.size(); 
  }
public:

  typedef REC_t rec_type;

  KmerMap() : _hash(), _n_rec(0) {}

  KmerMap(const vec<REC_t> & kvec, 
          const float ratio = 1.5) : 
    _hash(),
    _n_rec(0)
  { from_kmer_vec(kvec, ratio); }


  // ---- build the hash from a kmer record vector

  void from_kmer_vec(const vec<REC_t> & kvec, 
                     const float ratio = 1.5,
		     const unsigned verbosity = 1)
  {
    ForceAssertGt(ratio, 1.0);
    _n_rec = kvec.size();

    _hash.clear();
    // compute hash size --if we have zero recs, we still want a hash 
    // -- it will just have one empty (!valid) kmer.  This is so we
    // avoid calculating a (hash mod zero) later on.
    // TS: actually, i guess we better make sure that we have at least one
    // empty slot.
    size_t hash_size = std::max<size_t>( ratio*_n_rec, _n_rec+1 );
    _hash.resize(hash_size , REC_t(0));
    // if (verbosity) cout << "nh= " << _hash.size() << endl;

    // ---- set up the hash for each value in _kvec
    // for (size_t i = 0; i != _n_rec; dots_pct(i++, _n_rec, verbosity > 0)) {
    for (size_t i = 0; i != _n_rec; i++) {
      const REC_t & rec = kvec[i];
      size_t ih = _ih0_from_kmer(rec);
      while (_hash[ih].is_valid_kmer())
        _ih_next(ih);
      _hash[ih] = rec;
    }
  }


  // ---- kmer record index if kmer record exists; otherwise, index of empty record

  size_t ih_of_kmer(const Kmer_t & kmer) const
  {
    size_t ih = _ih0_from_kmer(kmer);
    while (true) {
      const REC_t & rec = _hash[ih];
      if (rec.is_valid_kmer() && ! kmer.match(rec))
        _ih_next(ih);
      else 
        return ih;
    }
  }


  // ---- insert a new record

  void insert_rec(const REC_t & rec_in)
  {
    ForceAssertLt(num_recs() + 1, size_hash());
    const size_t ih = ih_of_kmer(rec_in);

    if (!_hash[ih].is_valid_kmer()) // inserting rather than replacing
      _n_rec++;
    _hash[ih] = rec_in;  
  }



  REC_t   operator()(const Kmer_t & kmer) const { return _hash[ih_of_kmer(kmer)]; }
  REC_t & operator()(const Kmer_t & kmer)       { return _hash[ih_of_kmer(kmer)]; }

  REC_t   operator[](const size_t ih) const { return _hash[ih]; }
  REC_t & operator[](const size_t ih)       { return _hash[ih]; }

  size_t  size_hash()                  const { return _hash.size(); }
  size_t  num_recs()                   const { return _n_rec; }
  float   ratio()                      const { return float(size_hash()) / float(num_recs()); }


  
  // ---- outputs the frequencies of number of steps to find a record 
  void report() const
  {
    vec<size_t> freqs(1, 0);
    size_t not_found = 0;
    const size_t nh = _hash.size();
    for (size_t ih = 0; ih < nh; ih++) {
      const REC_t & rec0 = _hash[ih];
      if (rec0) {
        size_t jh = _ih0_from_kmer(rec0);
        REC_t rec = _hash[jh]; 
        size_t n_seeks = 1;
        while (rec && !rec0.match(rec)) {
          rec = _hash[_ih_next(jh)];
          n_seeks++;
        }
        
        
        
        if (_hash[jh]) {
          if (n_seeks >= freqs.size()) freqs.resize(n_seeks, 0);
          freqs[n_seeks]++;
        }
        else {
          not_found++;
        }
      }
      else {
        freqs[0]++;
      }
    }
    for (size_t i = 0; i != freqs.size(); i++) 
      cout << setw(10) << i << " " << setw(10) << freqs[i] << endl;
    cout << setw(10) << not_found << " not found." <<endl;
  }



  void write_binary(const String & fn) const
  {
    BinaryWriter::writeFile(fn.c_str(), _hash);
  }
  
  size_t read_binary(const String & fn) 
  {
    BinaryReader::readFile(fn.c_str(), &_hash);
    _n_rec = 0;
    const size_t nh = _hash.size();
    for (size_t ih = 0; ih < nh; ih++)
      if (_hash[ih].is_valid_kmer()) _n_rec++;
    return _hash.size();
  }
  


};

















#endif

