/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "math/Combinatorics.h"
#include "CoreTools.h"
#include <cmath>
#include <limits>

double Choices(int n, int k)
{
  if (k<0)
    return 0.0;
  double prod = 1.0;
  for ( int j = 1; j <= k; ++j )
    prod *= double(n + 1 - j)/j;
  return prod;
}
  
double BinomialProb(double p, int n, int k)
{
  return Choices(n,k) * pow(p, k) * pow(1.0-p, n-k);
}

namespace
{

inline double minTerm( double val )
{
    double const MIN_TERM = 2.*std::numeric_limits<double>::min();
    return ( fabs(val) < MIN_TERM ) ? MIN_TERM : val;
}

double beta_cont_frac( double a, double b, double x )
{
    int const N_ITERATIONS = 300;
    double const EPS = 2.*std::numeric_limits<double>::epsilon();
    int nnn = 0;

    double num = 1.;
    double denom = 1./minTerm(1. - (a+b)*x/(a+1.));
    double cf = denom;
    double aLess1 = a - 1.;
    double aPlus1 = a + 1.;
    double aPlusb = a + b;

    while ( ++nnn <= N_ITERATIONS )
    {
        double coeff = nnn*(b - nnn)*x/((aLess1 + 2*nnn)*(a + 2*nnn));
        num = minTerm(1. + coeff/num);
        denom = 1./minTerm(1. + coeff*denom);
        cf *= num*denom;

        coeff = -(a + nnn)*(aPlusb + nnn)*x/((a + 2*nnn)*(aPlus1 + 2*nnn));
        num = minTerm(1. + coeff/num);
        denom = 1./minTerm(1. + coeff*denom);

        double delta = denom*num;
        cf *= delta;
        if ( fabs(delta - 1.) < EPS )
            break;
    }

    ForceAssertLe(nnn,N_ITERATIONS);
    return cf;
}

inline double betaI( double a, double b, double x )
{
    ForceAssertGt(a,0.);
    ForceAssertGt(b,0.);
    ForceAssertGe(x,0.);
    ForceAssertLe(x,1.);

    if ( x == 0. )
        return 0.;
    if ( x == 1. )
        return 1.;

    double ln_beta = lgamma(a) + lgamma(b) - lgamma(a + b);
    double val = exp(a*log(x) + b*log1p(-x) - ln_beta);

    if ( x < (a + 1.)/(a + b + 2.) )
        return val*beta_cont_frac(a,b,x)/a;

    return 1. - val*beta_cont_frac(b,a,1.-x)/b;
}

}

double BinomialProbCum( double p, int n, int k )
{
    return betaI(n-k,k+1,1.-p);
}
