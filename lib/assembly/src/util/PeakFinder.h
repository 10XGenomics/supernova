///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef PEAK_FINDER_H
#define PEAK_FINDER_H

#include "math/Functions.h"
#include "Vec.h"

template<class X, class Y> class PeakFinder {

public:

    double window = 0.05;  // fractional window size
    size_t min_shoulder = 10;  // minimum number of points on either side of a valid peak
    int64_t min_peak_ratio = 10000; // filter out peaks smaller than max_peak/min_peak_ratio
    double min_peak_height_ratio = 1.2; // filter out peaks whose troughs are within this ratio
    bool debug = false;


    vec<size_t> FindPeaks(const vec<Y>& y_val) {

	Y global_peak = *std::max_element(y_val.begin(), y_val.end());

	vec<size_t> candidates;
	
	if (y_val.size() > min_shoulder * 2 )
	    for (auto it = y_val.begin() + min_shoulder; it < y_val.end() - min_shoulder; it++) {
		auto peak = std::max_element(it - min_shoulder, it + min_shoulder + 1 );
		if (peak ==  it) {
		    // background noise filter
		    if (*it < global_peak/min_peak_ratio)
			continue;

		    candidates.push_back(it - y_val.begin());
		}
	    }

	return candidates;
    }

   
    vec<size_t> FindPeaks(const vec<X>& x_val, const vec<Y>& y_val) {
	
	ForceAssertEq(x_val.size(), y_val.size());

	vec<size_t> candidates;

	// No data, no peaks.
	if (x_val.size() == 0)
	    return candidates;
	
	// Find all peaks with minimum shoulder - ignoring the x-range window size.
	vec<size_t> simple_peaks = FindPeaks(y_val);

	// Filter out peaks that aren't the largest point within the x-range window size.

	for (auto i : simple_peaks) {
	    X center_x = x_val[i];
	    Y center_y = y_val[i];

	    // Determine left (inclusive) and right (exclusive) indices
	    
	    // left = x  > (1 - window) * center_x (range start - inclusive)
	    auto left = std::upper_bound(x_val.begin(), x_val.end(), center_x * (1.0 - window));
	    // right = x > (1 + window) * center_x (range end - exclusive)
	    auto right = std::upper_bound(x_val.begin(), x_val.end(), center_x * (1.0 + window));
	    size_t left_index = left - x_val.begin();  // inclusive
	    size_t right_index = right - x_val.begin();  // exclusive
	    
	    // edge of data filter
	    if (left == x_val.begin() || right == x_val.end())
		continue;
	    
	    // sparse data filter
	    if (i - left_index < min_shoulder)
		continue;
	    if (right_index - i - 1 < min_shoulder)
		continue;
	    
	    // Find the peak position in the window
	    auto peak = std::max_element(y_val.begin() + left_index, y_val.begin() + right_index);
	    
	    if (peak == y_val.begin() + i) {
		
		candidates.push_back(i);
		if (debug) {
		    cout << "CANDIDATES :" << i << "  " << center_x << "  " << center_y << "  " << right_index - i - 1 << "  " << i - left_index << endl;
		}
	    }
	    
	}

	size_t peak_count = candidates.size();

	// filter peaks without deep enough troughs on either side
	vec<Bool> peaks_to_erase(peak_count, false);
	for (size_t i = 0; i < peak_count; i++) {
	    size_t left_peak = (i == 0 ? 0 : candidates[i - 1]);
	    size_t this_peak = candidates[i];
	    size_t right_peak = (i == peak_count - 1 ? x_val.size() : candidates[i + 1]);
	    Y left_min = *std::min_element(y_val.begin() + left_peak, y_val.begin() + this_peak);
	    Y right_min = *std::min_element(y_val.begin() + this_peak, y_val.begin() + right_peak);
	    if (Max(left_min, right_min) * min_peak_height_ratio > y_val[this_peak]) {
		peaks_to_erase[i] = true;
		if (debug)
		    cout << "SHALLOW: " << i << "  " << x_val[this_peak] << " " << y_val[this_peak] << " " << 1.0 * y_val[this_peak]/left_min << " " << 1.0 * y_val[this_peak]/right_min << endl;
	    }
	}
	EraseIf(candidates, peaks_to_erase);

	peak_count = candidates.size();

	// Centralize peaks - very basic
	for (size_t i = 0; i < peak_count; i++) {
	    size_t start = candidates[i];
	    size_t end = start + 1;
	    while (end < y_val.size() && y_val[end] == y_val[start])
		end++;
	    candidates[i] = start + (end - start - 1) / 2;
	    if (debug)
		cout << "CENTRED:" << i << "  " << candidates[i] << "  " << x_val[candidates[i]] << " " << y_val[candidates[i]] << endl;
	}
	
	return candidates;
    }
};

#endif
