///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "math/Array.h"
#include "CoreTools.h"
#include "Vec.h"

vec<double> GaussianKernal( int delta, int width_n_delta ) {
    // prepare the kernal function for smoothing
    vec<double> kernal( 2 * width_n_delta * delta + 1 );
    int center = width_n_delta * delta; // center position of the gaussian
    for ( int i = 0; i < kernal.isize(); i++ ) {
        int x = i - center;
        kernal[i] = exp( - x * x / double( 2 * delta * delta ) );
    }
    double total = Sum( kernal );
    for ( int i = 0; i < kernal.isize(); i++ ) 
        kernal[i] /= total;
    return kernal;
}

void SmoothArrayGaussian ( vec<double> & distr, int delta )
{
    vec<double> kernal = GaussianKernal( delta, 4 );
    int center = (kernal.size() - 1) /2;
    // smoothed array
    vec<double> distr_s( distr.size(), 0 );
    for ( int i = 0; i < distr.isize(); i++ ) 
        for ( int k = 0; k < kernal.isize(); k++ ) {
            int x = i + k - center;
            if ( x < 0 || x >= distr.isize() ) continue;
            distr_s[x] += distr[i] * kernal[k];
        }
    swap(distr, distr_s);
}

void SmoothArrayGaussian ( StdMap<int, double> & distr, int delta )
{
    vec<double> kernal = GaussianKernal( delta );
    int center = (kernal.size() - 1) /2;
     // smoothed array
     StdMap<int, double> distr_s;
     for ( map<int, double>::iterator it = distr.begin(); 
               it != distr.end(); it++ )
          for ( int k = 0; k < kernal.isize(); k++ ) {
               int x = it->first  + k - center;
               distr_s[x] += it->second * kernal[k];
          }
     swap(distr, distr_s);
}

