// Copyright (c) 2003 Broad Institute/Massachusetts Institute of Technology

#ifndef SHUFFLE_H
#define SHUFFLE_H

#include "Vec.h"
#include "random/RNGen.h"


/**
 * Shuffle
 *
 * Given an integer N>0, it will fill shuffled with the shuffled integers
 * between 0 and N-1 (included).
 *
 * It uses std::random_shuffle and a function object which uses
 * drand48_r. Thus, it is thread safe and multiprocessing safe, in the
 * sense that the same seed will always produce the same sequence.
 */
void Shuffle( int N, vec<int> &shuffled, int seed = 0 );
void Shuffle64( uint64_t N, vec<uint64_t> &shuffled, uint64_t seed = 0 );

/**
 * Shuffle that actually is independent of parallelism.
 * The functions above were thread-safe, but don't actually guarantee the same
 * order in a multi-threaded context, because they use a shared RNG.
 * This one uses an independent RNG, which, if seeded in a thread-independent
 * way, will produce stable results in a multi-threaded context.
 */
template <class Itr> // Itr is a random access iterator
void Shuffle( Itr const& beg, Itr const& end, unsigned seed = RNGen::random() )
{ RNGen rng(seed);
  if ( beg !=  end )
  { long idx = 1; using std::iter_swap;
    Itr itr(beg);
    while ( ++itr != end )
      iter_swap(itr, beg+rng.next()%++idx); } }

#endif
