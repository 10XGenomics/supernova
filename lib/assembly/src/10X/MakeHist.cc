// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

// Tool to compute histogram for an arbitrary vec<T> and write to JSON
#include "10X/MakeHist.h"

// #define TESTHIST

// compute histogram of a vector
template <class T, class VOut>
void ComputeHist(const vec<T>& in, VOut& out, T& minVal, T& maxVal, T& binsize){

    if ( in.size() == 0 ) {
        maxVal = T(0);
        minVal = T(0);
        out.resize(0);
        return;
    }

    maxVal= in.back();
    minVal = maxVal;
    for(auto & val : in)
            maxVal = std::max(maxVal,val), minVal = std::min(minVal,val);
    ForceAssertGt(binsize,T(0));

    // heuristic
    if((double)abs(maxVal-minVal)<1e-9){
        out.push_back(in.size());
        return;
    }

    // recompute bins
    size_t nBins = size_t((maxVal-minVal)/binsize);
    while(nBins == 0){ 
        binsize = binsize/10;
        nBins = size_t((maxVal-minVal)/binsize);
    }

    out.clear();
    // add an extra bin to handle outlier case
    out.resize(nBins+1,0);
    for(auto const& val : in)
        out[size_t((val-minVal)/binsize)]++;
    // remove extra bin if empty
    if (out.back()==0)
       out.resize(out.size()-1);
}

template void ComputeHist(const vec<int>& in, vec<int64_t>& out, int& minVal, int& maxVal, int& binsize);

template void ComputeHist(const vec<int64_t>& in, vec<int64_t>& out, int64_t& minVal, int64_t& maxVal, int64_t& binsize);

template void ComputeHist(const vec<double>& in, vec<int64_t>& out, double& minVal, double& maxVal, double& binsize);

// compute histogram of sizes
template <class VVIn, class T, class VOut>
void ComputeHistSizes(const VVIn& in, VOut& out, T& minVal, T& maxVal, T& binsize){

    vec<T> inSize; inSize.reserve(in.size());
    for(auto const& vec : in)
        inSize.push_back(vec.size());

    ComputeHist<T,VOut>(inSize,out,minVal,maxVal,binsize);
}


// write histogram to json
template <class T, class VOut>
void WriteHistToJson(const VOut& out, const T minVal, const T maxVal, const T binsize, String dir, String desc, String stage){
    
    ofstream OS;
    String filename = dir+"/histogram_"+desc+".json";
    OS.open(filename);
    if(OS){
        OS << "{" << endl;
        OS << "\t" << "\"description\": " << "\"" << desc << "\"," << endl;
        OS << "\t" << "\"stage\": " << "\"" << stage << "\"," << endl;
        OS << "\t" << "\"binsize\": " << binsize << "," << endl;
        OS << "\t" << "\"min\": " << minVal << "," << endl;
        OS << "\t" << "\"max\": " << maxVal << "," << endl;
        OS << "\t" << "\"numbins\": " << out.size() << "," << endl;
        OS << "\t" << "\"vals\": [" ;
        for(size_t i = 0; i<out.size(); i++){
            OS << out[i];
            if(i!=out.size()-1)
                OS << "," ;
        }
        OS << "]" <<endl;
        OS << "}" <<endl;

        OS.close();
    }
}

template void WriteHistToJson(const vec<double>& out, const double minVal, const double maxVal, const double binsize, String dir, String desc, String stage);

template void WriteHistToJson(const vec<int64_t>& out, const int minVal, const int maxVal, const int binsize, String dir, String desc, String stage);

template void WriteHistToJson(const vec<int64_t>& out, const int64_t minVal, const int64_t maxVal, const int64_t binsize, String dir, String desc, String stage);

template void WriteHistToJson(const vec<int64_t>& out, const double minVal, const double maxVal, const double binsize, String dir, String desc, String stage);

template<class T>
void ComputeAndWriteHist(const vec<T> VIn , T binsize, String dir, String desc, String stage, Bool verbose){

     auto start=WallClockTime();

     // input
     T maxV= T(0), minV=T(0);

     vec<int64_t> hist;

     // compute histogram
     if (verbose)
          cout << Date() << ": computing histogram" << endl;
     ComputeHist<T,vec<int64_t>>(VIn,hist,minV,maxV,binsize);

     // write json file
     if (verbose)
          cout << Date() << ": writing json file" << endl;
     WriteHistToJson<T,vec<int64_t>>(hist,minV,maxV,binsize,dir,desc,stage);

     // checksum
     size_t sum = 0ul;
     for(auto & h: hist)
         sum+=h;
     if(sum!=VIn.size())
         FatalErr("Incorrect histogram size\n");
     if (verbose)
          cout << "done, in " << TimeSince(start) << endl;
}

template void ComputeAndWriteHist(const vec<int> VIn , int binsize, String dir, String desc, String stage, Bool verbose);

template void ComputeAndWriteHist(const vec<int64_t> VIn , int64_t binsize, String dir, String desc, String stage, Bool verbose);

template void ComputeAndWriteHist(const vec<double> VIn , double binsize, String dir, String desc, String stage, Bool verbose);

