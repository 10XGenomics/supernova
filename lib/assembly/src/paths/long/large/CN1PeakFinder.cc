///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "util/PeakFinder.h"

#include "paths/long/large/CN1PeakFinder.h"
#include "math/Functions.h"
#include "Vec.h"

void CN1PeakFinder::PrefilterHighCNPeaks (const vec<double>& coverage, const vec<int64_t>& mass) {
    size_t peak_count = candidates.size();

    if (peak_count < 2)
	return;
    
    // find the largest peak, assume this is CN1 or near to it
    double max_coverage = coverage[candidates[MaxPeak(mass)]];
    
    if (debug)
	cout << "Prefilter assumed CN1 peak coverage: " << max_coverage << endl;
    
    // Eliminate peaks more than 5 times larger than assumed CN1 peak
    size_t max_candidate = 0;
    while (max_candidate < peak_count && 
	   (coverage[candidates[max_candidate]] <= high_CN_prefilter * max_coverage))
	max_candidate++;
    candidates.resize(max_candidate);
    
    if (debug)
	cout << "Prefilter removed last " << peak_count - max_candidate << " peaks" << endl;
    
    return;
}


double CN1PeakFinder::FindPeak (const vec<double>& coverage, const vec<int64_t>& mass) {

    // Special case if no data is supplied.
    if (mass.size() == 0) {
	if (debug) cout << "No valid data, cannot find peak" << endl;
	return 0;
    }

    PeakFinder<double,int64_t> peak_finder;
    candidates = peak_finder.FindPeaks(coverage,mass);
    
    PrefilterHighCNPeaks(coverage,mass);

    size_t peak_count = candidates.size();
	
    if (peak_count == 1 ) {
	// Single peak found
	cn_peaks.push_back(candidates[0]);
	cn_values.push_back(1);
	if (debug) cout << "Single peak" << endl;
    } else if (peak_count == 0) {
	// No peak found - use max mass instead
	cn_peaks.push_back(std::max_element(mass.begin(), mass.end()) - mass.begin());
	cn_values.push_back(1);
	if (debug) cout << "No peaks - using max value" << endl;
    } else {
	
	// Score all possible CN1 peaks
	
	size_t max_peak = MaxPeak(mass);

	int best_score = 0;
	vec<int> best_used;
	for (size_t i = 0; i < peak_count; i++) {
	    int score = 0;
		
	    double base_cov = coverage[candidates[i]];

	    vec<int> used(peak_count, 0);
	    used[i] = 1;
	    
	    // Handle special diploid case
	    if (i > 0)
		MatchPeak(coverage, used, base_cov, 0.5);
	    
	    // Look for CN2,3,4,5... peaks
	    for (size_t multiplier = 2; multiplier <= high_CN_prefilter; multiplier++) {
		MatchPeak(coverage, used, base_cov, multiplier);
	    }
	    
	    // Compute the score
	    score = std::count_if(used.begin(), used.end(),  [] (int val) {return val != 0;});
	    if (debug)
		cout << "SCORE: " << i << " " << score << " " << ToStringBool(diploid) << endl;
		
	    // Only consider if we used the largest peak in some way
	    if (used[max_peak] != 0 ) {
		if (score == best_score) {
		    // prefer diploid scores if the 1st peak is 10x smaller than the 2nd
		    size_t dip = std::find(used.begin(), used.end(), -2) - used.begin();
		    if (dip != used.size() && mass[candidates[dip]] * 10 < mass[candidates[i]]) {
			best_score = score;
			best_used = used;
		    }
		}  else if (score > best_score) {
		    best_score = score;
		    best_used = used;
		}
	    }
	}
	

	// Determine used peaks
	for (size_t i = 0; i < best_used.size(); i++) {
	    if (best_used[i] != 0) {
		cn_peaks.push_back(candidates[i]);
		cn_values.push_back(best_used[i]);
	    }
	}
    }

    // Compute coverage using largest of the first two peaks to improve accuracy.
    
    if (cn_peaks.size() > 1 && mass[cn_peaks[0]] < mass[cn_peaks[1]])
	cn1_coverage = coverage[cn_peaks[1]] / 2.0;
    else 
	cn1_coverage = coverage[cn_peaks[0]];
    
    diploid = (cn_values[0] == -2);
    
    return cn1_coverage;
}

    
bool CN1PeakFinder::MatchPeak(const vec<double>& coverage, vec<int>& used, double base, double multiplier) {
    double target = base * multiplier;
    for (size_t i = 0; i < used.size(); i++) {
	if (used[i] == 0 && Abs(target - coverage[candidates[i]]) < max_peak_tolerance * target) {
	    used[i] = (multiplier >= 1 ? multiplier : -1.0/multiplier );
	    return true;
	}
    }
    return false;
}


// find the largest peak, assume this is CN1 or near to it
size_t CN1PeakFinder::MaxPeak(const vec<int64_t>& mass) {
    size_t max_peak = 0;
    for (size_t i = 0; i < candidates.size(); i++)  {
	if (mass[candidates[i]] > mass[candidates[max_peak]]) 
	    max_peak = i;
    }
    return max_peak;
}
