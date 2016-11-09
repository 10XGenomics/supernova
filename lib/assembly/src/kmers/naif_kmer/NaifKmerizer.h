///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

//
//   
//  NaifKmerizer.h 
//
//  This header file provides a templatized framework for efficiently 
//  generating and processing kmers in parallel.
//
//
//  Here only the general algorithm is defined.  All the specifics of how the 
//  kmers are generated and what to do with them (count them, store them, 
//  etc...) is done by a kernel object that must be provided.  
//
//  Use by simply calling the function:
// 
//     naif_kmerize(& kernel, num_threads, mem_mean);  
//
//
//  The kernel MUST provide the following:
//      
//     A copy constructor that generates temporary kernels (see examples)
//
//     kernel::rec_type
//     kernel.K()
//     kernel.bases()     // a reference to the vector of base vectors
//     kernel.parse_base_vec(& sub_parcel_buffer, i_read)
//     kernel.summarize(parcel, i, j)
//     kernel.merge(tmp_kernel)
//
//  
//  The algorithm is the following: 
//
//    1. Determine the number of parcels and passes so that the whole process 
//       fits in memory.
//
//    2. On each pass and on each thread, parse a block of reads (there are as 
//       many blocks as threads) for kmers belonging to a subset of parcels 
//       (as many subsets as threads).  Stores them in buffers, one for each 
//       parcel.  Total of n_threads^2 buffers.
//
//    3. Wait for all threads to finish parsing their read blocks. 
// 
//    3. Combine all the kmer data collected in the buffers into a single kmer
//       parcel per thread.  
//
//    4. Sort all the parcels.  Each thread sorts their own.
//  
//    5. Parse parcels for consecutive kmers.  Give each kmer batch to a 
//       temporary kernel copy for processing. 
//
//    6. Merge results stored in the temporary kernels into the main kernel.  
//       All thread locking is left to the kernel. 
//
//    7. Wait for all threads to be done with this pass before moving to the 
//       next pass.
//  
  



#ifndef KMERS__KMER_DATA__NAIF_KMERIZER__H
#define KMERS__KMER_DATA__NAIF_KMERIZER__H

#include "system/SysConf.h"  // processorsOnline()
#include "system/WorklistN.h"


static inline 
String TagNK(String S = "NK") { return Date() + " (" + S + "): "; } 







template<class REC_t>
class ParcelBuffer : public vec< vec<REC_t> >
{
  const size_t _n_parcels;
  size_t       _i0_parcel;
  size_t       _i1_parcel;
  
public: 
  ParcelBuffer(const size_t i0_parcel, 
               const size_t i1_parcel, 
               const size_t n_parcels) :
    vec< vec<REC_t> >(i1_parcel - i0_parcel),
    _n_parcels(n_parcels),
    _i0_parcel(i0_parcel),
    _i1_parcel(i1_parcel)
  {}

  size_t i_parcel_compute(const REC_t & rec) const 
  {
    // only look at first 20 bits
    return  ((rec.hash_64bits() & 0x000fffff) * _n_parcels) >> 20;
  }

  bool in_one_parcel(const REC_t & rec) const 
  {
    const size_t i_parcel = i_parcel_compute(rec);
    return (_i0_parcel <= i_parcel && i_parcel < _i1_parcel);
  }

  void add(const REC_t & rec)
  {
    const size_t n_sub_parcels = this->size();
    const size_t i_sub_parcel = i_parcel_compute(rec) % n_sub_parcels;
    (*this)[i_sub_parcel].push_back(rec);
  }
};


// -------- a class to keep track of various times

class TimerNK
{
public:
  double parse;
  double collect;
  double sort;
  double summarize;
  double merge;
  double wait;
  double idle;

  TimerNK() :
    parse(0),
    collect(0),
    sort(0),
    summarize(0),
    merge(0),
    wait(0),
    idle(0)
  {}

  TimerNK(const vec<TimerNK> & timers) :
    parse(0),
    collect(0),
    sort(0),
    summarize(0),
    merge(0),
    wait(0),
    idle(0)
  {
    const size_t n = timers.size();
    double max = 0;
    for (size_t i = 0; i < n; i++) {
      parse     += timers[i].parse     / n;
      collect   += timers[i].collect   / n;
      sort      += timers[i].sort      / n;
      summarize += timers[i].summarize / n; 
      merge     += timers[i].merge     / n;
      wait      += timers[i].wait      / n;
      const double total = timers[i].total();
      if (total > max) max = total;
    }
    for (size_t i = 0; i < n; i++) 
      idle += max - timers[i].total();
    idle /= n;
  }

  double total() const { return parse + collect + sort + summarize + merge + wait; }

  void print() const 
  {
    double tot = total();
    cout << fixed << setprecision(0);
    cout << TagNK() << "  time total       = " 
         << setw(6) << tot       << "s (100%)" << endl;
    cout << TagNK() << "  time parsing     = " 
         << setw(6) << parse     << "s (" << setw(3) << 100 * parse / tot << "%)" << endl;
    cout << TagNK() << "  time collecting  = " 
         << setw(6) << collect   << "s (" << setw(3) << 100 * collect / tot << "%)" << endl;
    cout << TagNK() << "  time sorting     = " 
         << setw(6) << sort      << "s (" << setw(3) << 100 * sort / tot << "%)" << endl;
    cout << TagNK() << "  time summarizing = " 
         << setw(6) << summarize << "s (" << setw(3) << 100 * summarize / tot << "%)" << endl;
    cout << TagNK() << "  time merging     = " 
         << setw(6) << merge     << "s (" << setw(3) << 100 * merge / tot << "%)" << endl;
    cout << TagNK() << "  time waiting     = " 
         << setw(6) << wait      << "s (" << setw(3) << 100 * wait / tot << "%)" << endl;
    cout << TagNK() << "  time idling      = " 
         << setw(6) << idle      << "s (" << setw(3) << 100 * idle / tot << "%)" << endl;
  }
};


// ----------------------------------------
//  ParcelProc<KERNEL_t>
// ----------------------------------------



template<class KERNEL_t>
class ParcelProc 
{
private:
  typedef typename KERNEL_t::rec_type rec_t;

  const size_t                 _n_threads;
  const size_t                 _n_passes;
  KERNEL_t                   * _p_kernel_main;

  const vec<size_t>          & _i1_bv;

  // one ParcelBuffer<rec_t> for each thread
  vec<ParcelBuffer<rec_t> *> & _parcel_bufs_pt; 
  vec<size_t>                & _sync_states; 
  vec<TimerNK>               & _timers;
  unsigned                     _verbosity;


public:
  ParcelProc(const size_t                  n_threads,
             const size_t                  n_passes,
             KERNEL_t                    * p_kernel_main,
             const vec<size_t>           & i1_bv,
             vec<ParcelBuffer<rec_t> * > & parcel_bufs_pt,
             vec<size_t>                 & sync_states,
             vec<TimerNK>                & timers, 
             const unsigned                verbosity = 1) :
    _n_threads(n_threads),
    _n_passes(n_passes),
    _p_kernel_main(p_kernel_main),
    _i1_bv(i1_bv),
    _parcel_bufs_pt(parcel_bufs_pt),
    _sync_states(sync_states),
    _timers(timers),
    _verbosity(verbosity)
  {}

  // copy constructor
  ParcelProc(const ParcelProc<KERNEL_t> & that) :
    _n_threads(that._n_threads),
    _n_passes(that._n_passes),
    _p_kernel_main(that._p_kernel_main),
    _i1_bv(that._i1_bv),
    _parcel_bufs_pt(that._parcel_bufs_pt),
    _sync_states(that._sync_states),
    _timers(that._timers),
    _verbosity(that._verbosity)
  {}


private:
  // ---- this is just a dumb way of synchronizing the threads
  void _naif_sync(const size_t i_thread)
  {
    _sync_states[i_thread]++;
    {
      bool done = false;
      while (!done) {
        done = true;
        for (size_t i = 0; i < _n_threads; i++)
          if (_sync_states[i] < _sync_states[i_thread]) 
            done = false;
      }
    } 
  }
  
  
  // ---- add all the blocks into block 0 to build the full parcel
  vec<rec_t> & _build_parcel(const size_t i_sub_parcel)
  {
    vec<rec_t> & parcel = (*(_parcel_bufs_pt[0]))[i_sub_parcel];
    const size_t n_blks = _n_threads;
    
    size_t n_recs = 0;
    for (size_t i_blk = 0; i_blk < n_blks; i_blk++)
      n_recs += (*(_parcel_bufs_pt[i_blk]))[i_sub_parcel].size();
    parcel.reserve(n_recs);

    if (_verbosity >= 3) 
      cout << TagNK() << "i_sub_parcel= " << i_sub_parcel 
           << ": reserved space for " << n_recs << " records (n_blks= " << n_blks << ")." << endl;
          
    // start at 1 because 0 is the full parcel
    for (size_t i_blk = 1; i_blk < n_blks; i_blk++) {
      vec<rec_t> & parcel_blk = (*(_parcel_bufs_pt[i_blk]))[i_sub_parcel];
      parcel.insert(parcel.end(), parcel_blk.begin(), parcel_blk.end());
    }
    return parcel;
  }


  
  size_t _consecutive_records_find_end(const vec<rec_t> & recs, 
                                       const size_t i0)
  {
    // **** can optimize with binary forward search!!!! ****
    // actually, not worth the complexity.  This is not a bottleneck. far from it.
    
    const size_t n_recs = recs.size();
    const rec_t & rec = recs[i0];
    size_t i1 = i0 + 1;
    while (i1 < n_recs && recs[i1].match(rec))
      i1++;
    return i1;
  }




public:
  // ---- called by the Worklist by each thread
  void operator() (const size_t i_thread) 
  {
    const String str_thread = (i_thread < 10 ? " [thread  " : " [thread ") + ToString(i_thread) + "]";
    
    TimerNK & timer = _timers[i_thread];

    // ---- some declarations
    const size_t n_parcels     = _n_passes * _n_threads;
    const size_t n_sub_parcels = _n_threads;
    const size_t i_sub_parcel  = i_thread;
    const size_t i_blk         = i_thread;

    
    for (size_t i_pass = 0; i_pass < _n_passes; i_pass++) {
      
      const String str_pass = (i_pass < 10 ? "[pass  " : "[pass ") + ToString(i_pass) + "]";

      if (i_thread == 0) {
        if (_verbosity == 1) 
          cout << TagNK() << str_pass << " parse" << flush;
        if (_verbosity >= 2) 
          cout << TagNK() << str_pass << " Parsing reads for kmers." << endl;
      }

      // ---- define the starting and ending indices for the parcel subset

      const size_t i0_parcel = i_pass * n_sub_parcels;
      const size_t i1_parcel = i0_parcel + n_sub_parcels;


      // ---- declare parcel buffers for the current parcel subset

      ParcelBuffer<rec_t> sub_parcels(i0_parcel, i1_parcel, n_parcels);
      _parcel_bufs_pt[i_sub_parcel] = &sub_parcels; // update "global" pointers

      
      // ---- define the starting and ending base vector indices

      const size_t i0_bv = (i_blk == 0) ? 0 : _i1_bv[i_blk - 1];
      const size_t i1_bv = _i1_bv[i_blk];

      if (_verbosity >= 3) 
	cout << TagNK() << str_pass << str_thread 
	     << " bvs: [" << i0_bv << " : " << i1_bv << "] " 
	     << i1_bv - i0_bv << " to parse." << endl;


      // ---- parse base vectors block for kmers in parcel subset
      //
      //      'parsing' is where the user's kernel parses all the base vectors (reads)
      //      for a sub-group of kmers and other associated metadata of the user's choosing.
      //
      //      Each thread only looks at a subset of reads but searches for 
      //      for various subsets of kmers. 

      timer.parse -= WallClockTime();
      for (size_t i_bv = i0_bv; i_bv < i1_bv; i_bv++)
        _p_kernel_main->parse_base_vec(&sub_parcels, i_bv);
      timer.parse += WallClockTime();
        
      if (_verbosity >= 3) 
	cout << TagNK() << str_pass << str_thread << " Done parsing." << endl;


      // ---- a thread barrier to make sure all threads are here


      timer.wait -= WallClockTime();
      _naif_sync(i_thread);
      timer.wait += WallClockTime();


      // ---- collect this thread's parcel from all the sub parcels
      //
      //      Now, all the various kmer subsets, dispersed over all the threads 
      //      buffers, are brought together.
      
      const size_t i_parcel = i0_parcel + i_sub_parcel;

      if (_verbosity >= 3) 
	cout << TagNK() << str_pass << str_thread 
	     << " Collect parcel " << i_parcel << "." << endl;
      else if (i_thread == 0) {
        if (_verbosity >= 2) cout << TagNK() << str_pass << " Collecting parcels." << endl;
        if (_verbosity == 1) cout << " collect" << flush;
      }

        
      timer.collect -= WallClockTime();
      vec<rec_t> & parcel = _build_parcel(i_sub_parcel);
      timer.collect += WallClockTime();
      

      // ---- sort kmer records in this parcel

      if (_verbosity >= 3) 
	cout << TagNK() << str_pass << str_thread 
	     << " Done collecting parcel of size " << parcel.size() << ". Sorting now." << endl;
      else if (i_thread == 0) {
        if (_verbosity >= 2) cout << TagNK() << str_pass << " Sorting parcels." << endl;
        if (_verbosity == 1) cout << " sort" << flush;
      }
      
      timer.sort -= WallClockTime();
      sort(parcel.begin(), parcel.end());
      timer.sort += WallClockTime();

      
      // ---- summarize parcel into a temporary kernel
      // 
      //      'summarization' is when the user's temporary kernel actually processes 
      //      a kmer parcel.
      //
      //      A temporary kernel is created because it might have local data structures
      //      to be merged later in the 'merging' step.
      
      KERNEL_t kernel_tmp(*_p_kernel_main);
      const size_t n_recs = parcel.size();

      if (_verbosity >= 3) 
	cout << TagNK() << str_pass << str_thread 
	     << " Done sorting. Summarizing " << n_recs << " records now." << endl;
      else if (i_thread == 0) {
        if (_verbosity >= 2) cout << TagNK() << str_pass << " Summarizing consecutive kmers." << endl;
        if (_verbosity == 1) cout << " summarize" << flush;
      }
	      
      timer.summarize -= WallClockTime();
      size_t i0 = 0;
      while (i0 < n_recs) {
	const size_t i1 = _consecutive_records_find_end(parcel, i0);
	kernel_tmp.summarize(parcel, i0, i1);
	i0 = i1;
      }
      timer.summarize += WallClockTime();

      
      // ---- merge summarized parcel in temporary kernel into global kernel
      //
      //      'merging' is when the user's main kernel brings together 
      //      the results summarized by each temporary kernel.

      if (_verbosity >= 3) 
	cout << TagNK() << str_pass << str_thread 
	     << " Done summarizing. Merging now." << endl;
      else if (i_thread == 0) {
        if (_verbosity >= 2) cout << TagNK() << str_pass << " Merging results." << endl;
        if (_verbosity == 1) cout << " merge" << flush;
      }

      timer.merge -= WallClockTime();
      _p_kernel_main->merge(kernel_tmp, i_parcel);
      timer.merge += WallClockTime();
      

      // ---- a barrier to make sure all threads are here before parcel destructors kick in

      timer.wait -= WallClockTime();
      _naif_sync(i_thread);
      timer.wait += WallClockTime();
      
      if (i_thread == 0) {
        if (_verbosity >= 2) cout << TagNK() << str_pass << " Done." << endl;
        if (_verbosity == 1) cout << endl;
      }
    }

  }  // for (size_t i_pass = 0; i_pass < _n_passes; i_pass++)


};










// ------------------------
//   main routine to call
// ------------------------

template<class KERNEL_t> 
void naif_kmerize(KERNEL_t     * p_kernel_main,
                  const size_t   n_threads_arg,
                  const size_t   verbosity = 3,
                  const size_t   mem_mean_ceil_user = 0)
{
  typedef typename KERNEL_t::rec_type KmerRec_t;

  if (verbosity > 0) cout << TagNK() << "Starting kmerization." << endl;
  const double time_start = WallClockTime();
  
  const size_t n_threads = boundNumThreads(n_threads_arg);
 
  const BaseVecVec & bvv = p_kernel_main->bases();

  const size_t K    = (*p_kernel_main).K();
  const size_t n_bv = bvv.size();
 
  ForceAssertGt(n_bv, 0u);
    
  // ---- compute number of kmers 

  size_t n_k_total = 0;
  for (size_t i_bv = 0; i_bv < n_bv; i_bv++) {
    const size_t n_b = bvv[i_bv].size();
    if (n_b >= K) 
      n_k_total += n_b - K + 1;  // add the number of kmers on base vector 
  }
  if (n_k_total == 0) {
    cout << TagNK() << "No kmers to process. Something is wrong with the data. Aborting." << endl;
    CRD::exit(1);
  }

  // ---- compute ending indices for read blocks

  const size_t n_blks = n_threads;
  size_t n_k_per_blk = n_k_total / n_blks;
  vec<size_t> i1_bv(n_blks, n_bv);
  size_t n_k = 0;
  for (size_t i_bv = 0; i_bv < n_bv; i_bv++) {
    const size_t n_b = bvv[i_bv].size();
    if (n_b >= K) 
      n_k += n_b - K + 1;

    const size_t i_blk = (n_k == 0) ? 0 : (n_k - 1) * n_blks / n_k_total;
    ForceAssertLt(i_blk, n_blks);
    i1_bv[i_blk] = i_bv + 1;
  }
  
  sort(i1_bv.begin(), i1_bv.end()); // because some entries are skipped and left with 0
  ForceAssertEq(n_bv, i1_bv.back());

  // ---- compute number of parcels that fit in memory

  // need factor of 2 because kmers are buffered from different threads before
  // the kmer vec is built, so each kmer ends up being stored twice. 
  const size_t mem_needed = 2 * n_k_total * sizeof(KmerRec_t);
  const size_t mem_total = GetMaxMemory();
  const size_t mem_used  = MemUsageBytes();
  const size_t mem_avail = MemAvailable();
  ForceAssertGt(MemAvailable(), 0ul);

  // use, on average, 1/4 of available memory.
  size_t mem_mean = mem_avail / 4; 

  // if mem_mean > mem_mean_ceil_user, cap it.
  if (mem_mean_ceil_user > 0 && mem_mean > mem_mean_ceil_user)
    mem_mean = mem_mean_ceil_user;

  size_t n_passes  = 1 + (mem_needed - 1) / mem_mean;

  // run at least 10 passes, so that in large memory meachines and small genomes 
  // it doesn't look like it's hogging memory 

  if (n_passes < 10) n_passes = 10;  

  const size_t n_parcels = n_threads * n_passes;

  const size_t mem_to_use = mem_needed / n_passes;
    
  if (verbosity > 0) {
    cout << TagNK() << "  n_threads     = " << setw(18) << ToStringAddCommas(n_threads) << endl;
    cout << TagNK() << "  K             = " << setw(18) << ToStringAddCommas(K) << endl;
    cout << TagNK() << "  sizeof_krec   = " << setw(18) << ToStringAddCommas(sizeof(KmerRec_t)) << endl;
    cout << TagNK() << "  n_bv          = " << setw(18) << ToStringAddCommas(n_bv) << endl;
    cout << TagNK() << "  n_kmers       = " << setw(18) << ToStringAddCommas(n_k_total) << endl;
    cout << TagNK() << "  mem_total     = " << setw(18) << ToStringAddCommas(mem_total) << endl;
    cout << TagNK() << "  mem_used      = " << setw(18) << ToStringAddCommas(mem_used) << endl;
    cout << TagNK() << "  mem_avail     = " << setw(18) << ToStringAddCommas(mem_avail) << endl;
    cout << TagNK() << "  mem_mean      = " << setw(18) << ToStringAddCommas(mem_mean) << endl;
    cout << TagNK() << "  mem_mean_user = " << setw(18) << ToStringAddCommas(mem_mean_ceil_user) << endl;
    cout << TagNK() << "  mem_needed    = " << setw(18) << ToStringAddCommas(mem_needed) << endl;
    cout << TagNK() << "  mem_to_use    = " << setw(18) << ToStringAddCommas(mem_to_use) << endl;
    cout << TagNK() << "  n_parcels     = " << setw(18) << ToStringAddCommas(n_parcels) << endl;
    cout << TagNK() << "  n_passes      = " << setw(18) << ToStringAddCommas(n_passes) << endl;
  }

  // ---- parcel buffers pointers
  //      This is a bit awkward.
  //      First I tried declaring buffers once and reusing them at each pass.
  //      This leads to a continuous memory increase as each pass statistically can 
  //      need slightly larger buffers than the previous one. 
  //      To keep memory usage under control we must create and destroy
  //      the buffers at each pass.  
  //      Therefore, here we declare a "global" vector of pointers to those buffers.
  //      This vector must be updated at each pass upon buffer declaration.

  vec<ParcelBuffer<KmerRec_t> *> parcel_bufs_pt(n_threads);
    
  vec<size_t> sync_states(n_threads, 0); // ---- for syncing the threads 
  
  vec<TimerNK> timers(n_threads);

  // ---- run parallel

  if (true) {
    ParcelProc<KERNEL_t> proc(n_threads, 
                              n_passes,
                              p_kernel_main,
                              i1_bv,
                              parcel_bufs_pt,
                              sync_states,
                              timers,
			      verbosity);
    
    if (n_threads <= 1) { // no need for Worklist
      proc(0);
    }
    else {
      parallelFor(0ul,n_threads,proc,n_threads);
    }
  }

  if (verbosity > 0) {
    TimerNK timer(timers);
    timer.print();
    cout << TagNK() << "Done with kmerization. Took " << TimeSince(time_start) << "." << endl;
  }
}










#endif
