///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef CN1_PEAK_FINDER_H
#define CN1_PEAK_FINDER_H

#include "math/Functions.h"
#include "Vec.h"

class CN1PeakFinder {

public:

    vec<size_t> candidates;  // list of candidate peaks    

    vec<size_t> cn_peaks;  // list of peaks
    vec<int> cn_values;    // list of peak CNs
    double cn1_coverage; // CN1 coverage value - not necessarily exact value of first peak
    bool diploid = false;  // is the organism diploid


    uint64_t high_CN_prefilter = 5;  // filter peaks with presumed CN greater than this
    double max_peak_tolerance = 0.1;  // in scoring, max tolorance when matching peaks
    bool debug = false;

    double FindPeak (const vec<double>& coverage, const vec<int64_t>& mass);
    
private:
    bool MatchPeak(const vec<double>& coverage, vec<int>& used, double base, double multiplier = 1);

    // find the largest peak, assume this is CN1 or near to it
    size_t MaxPeak(const vec<int64_t>& mass);

    void PrefilterHighCNPeaks (const vec<double>& coverage, const vec<int64_t>& mass);

};

#endif
