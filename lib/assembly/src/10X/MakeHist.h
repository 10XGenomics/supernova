// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

#ifndef _makehist_h_
#define _makehist_h_

#include "MainTools.h"
#include "VecString.h"
#include "10X/SuperFiles.h"
#include "system/System.h"


template <class T, class VOut>
void ComputeHist(const vec<T>& in, VOut& out, T& minVal, T& maxVal, T& binsize);

template <class VVIn, class T, class VOut>
void ComputeHistSizes(const VVIn& in, VOut& out, T& minVal, T& maxVal, T& binsize);

template <class T, class VOut>
void WriteHistToJson(const VOut& out, const T minVal, const T maxVal, const T binsize, String dir, String desc, String stage="");

template<class T>
void ComputeAndWriteHist(const vec<T> VIn , T binsize, String dir, String desc, String stage="none", Bool verbose=False);

#endif
