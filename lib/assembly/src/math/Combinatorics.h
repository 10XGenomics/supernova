/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef COMBINATORICS_H
#define COMBINATORICS_H

#include "Vec.h"

/// Step through every way of turning on bits_on out of num_bits.
/// Designed for the case where you need to look through all
/// of the objects each time, and do something different based
/// on whether they are chosen or not.  If you only need to touch
/// the chosen things, this is wrong wrong wrong.
///
/// Sample usage:
///
/// vec<Bool> foo;
/// FirstCombination(foo, n, k);
/// do {
///   ...
/// } while ( NextCombination(foo) );

template<class T>
void FirstCombination( vec<T>& c, int num_bits, int bits_on ) {
  c.clear();
  c.resize(bits_on, 1);
  c.resize(num_bits, 0);
}

/// Knuths "Algorithm L": push up the first bit you can;
/// all earlier bits reset to its minimum possible position.
/// Returns false when there are no more combinations.
template<class T>
bool NextCombination( vec<T>& c ) {
  unsigned int i=0, j=0;
  for(; i<c.size() && !c[i]; i++)  ;
  for(; i<c.size()-1 && c[i+1]; i++) { c[i]=0; c[j++]=1; }
  if( i < c.size()-1 ) {c[i]=0; c[i+1]=1; return true;}
  return false;
}

// Number of combinations of n objects taken k at a time, following
// Knuth's definitions to extend the classical case when n is positive
// and 0<=k<=n.
double Choices(int n, int k);
  
// Binomial probability: probability that in n trials with prob. p of
// success we see exactly k successes.
double BinomialProb(double p, int n, int k);

// Compute the probability of getting k or fewer successes in 
// n trials with probability of success p.
//
// See also BinomialSum in random/Bernoulli.h.

double BinomialProbCum( double p, int n, int k );

// Probability of k or more failures in n trials with probability of failure p.
inline double BinomialProbCumFailure( double p, int n, int k )
{ return BinomialProbCum(1.-p,n,n-k); }

#endif
