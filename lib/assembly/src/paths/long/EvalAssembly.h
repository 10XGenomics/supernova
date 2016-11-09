///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
#ifndef EVAL_ASSEMBLY_H
#define EVAL_ASSEMBLY_H

#include "CoreTools.h"
#include <fstream>
#include <map>

class ref_data;
class HyperBasevector;
class align;

class AssemblyError {
public:
    typedef std::pair<int64_t,int64_t>p_c_t;
    //core data structure
    class err_per_ref{
    public:
        err_per_ref():range(-1),max_range(-1),m_GoodTopology(true){};
        void SetRange(uint64_t len,bool isDoubled){ range=len; max_range=len; if(isDoubled) range/=2; };
        void AddIndel(int64_t p, int64_t indel_size){ indels.emplace_back(p,indel_size); }
        void AddSub(int64_t p, int64_t n)           { subs.emplace_back(p,n); }
        void AddGap(int64_t p, int64_t gap_size)    { gaps.emplace_back(p,gap_size); }
        void AddLeftGap(int64_t in)                 { viLGap.push_back(in);};
        void AddRightGap(int64_t in)                { viRGap.push_back(in);};

        void SetBestPath(const vec<int>& hbp,const vec<pair<int,Bool>>& hbp_to_hb){
            best_path.clear_and_resize(hbp.size());
            std::transform(hbp.begin(),hbp.end(),best_path.begin(),[&hbp_to_hb](int p){return hbp_to_hb[p];});
        };
        void SetBestPath(const vec<pair<int,Bool>>& in) {best_path=in;};
        const vec<pair<int,Bool>>& GetBestPath()const{return best_path;};
        template<typename _OIter>
        void LogBestPath(_OIter itr)const{
            std::transform(best_path.begin(),best_path.end(),itr,[](const pair<int,Bool>&in){return in.first;});
        }


        bool GoodTopology()const {return m_GoodTopology && range==max_range;};
        template<typename _OIter>
        void LogGaps(_OIter itr)const { std::transform(gaps.begin(),gaps.end(),itr,[](const p_c_t& in){return in.second;}); }

        template<typename _OIter>
        void LogIndels(_OIter itr)const { std::transform(indels.begin(),indels.end(),itr,[](const p_c_t& in){return abs(in.second);}); }

        int CountSubs()const { return std::accumulate(subs.begin(),subs.end(),0,[](int i,const p_c_t& in){return i+in.second;}); }
        vec<int> viLGap;
        vec<int> viRGap;

        void Compact(const HyperBasevector& hbp, vec<std::tuple<int64_t,int64_t,int,int64_t,int64_t,int64_t,int64_t>> coors_hbp, const vec<pair<int,Bool>>&hbp_to_hb);
    private:
        vec<p_c_t> gaps;
        vec<p_c_t> indels;
        vec<p_c_t> subs;
        vec<pair<int,Bool>> best_path;
        int64_t range,max_range;
        bool m_GoodTopology;
    };

    //Modification functions
    void AddIndel(int g, int64_t p, int64_t indel_size) { g_err[g].AddIndel(p,indel_size); }
    void AddSub(int g, int64_t p, int64_t n)            { g_err[g].AddSub(p,n); }
    void AddGap(int g, int64_t p, int64_t gap_size)     { g_err[g].AddGap(p,gap_size); }
    void AddLeftGap(int g, int64_t in)                  { g_err[g].AddLeftGap(in);};
    void AddRightGap(int g, int64_t in)                 { g_err[g].AddRightGap(in);};
    err_per_ref& GetErrPerRef(int g){return g_err[g];};
    std::map<int,err_per_ref>::iterator begin(){return g_err.begin();};
    std::map<int,err_per_ref>::iterator end(){return g_err.end();};


    //Access-only functions
    void PrintSummary(ostream& out);

    vec<int> GetGaps()const{
        vec<int> out;
        for(auto&entry:g_err) entry.second.LogGaps(std::back_inserter(out));
        return out;
    }
    vec<int> GetIndels()const{
        vec<int> out;
        for(auto&entry:g_err) entry.second.LogIndels(std::back_inserter(out));
        return out;
    }
    int GetSubs()const{
        int out=0;
        for(auto&entry:g_err) out+=entry.second.CountSubs();
        return out;
    }

    vec<int> GetLeftGaps()const{
        vec<int> out;
        for(auto&entry:g_err) std::copy_if(entry.second.viLGap.begin(),entry.second.viLGap.end(),std::back_inserter(out) , [](int i){return i!=0;});
        return out;
    };
    vec<int> GetRightGaps()const{
        vec<int> out;
        for(auto&entry:g_err) std::copy_if(entry.second.viRGap.begin(),entry.second.viRGap.end(),std::back_inserter(out) , [](int i){return i!=0;});
        return out;
    };
    //returns a list of best paths, each correspond to a reference sequence, with which the errors are calculated
    vec<vec<int>> GetBestPaths()const{
        vec<vec<int>> out; out.reserve(g_err.size());
        for(const auto& ref: g_err){
            out.push_back(vec<int>());
            ref.second.LogBestPath(std::back_inserter(out.back()));
        }
        return out;
    };

private:
    std::map<int,err_per_ref> g_err;
};

class N50Calculator
{
public:
    N50Calculator( String const& logFile )
    : mPending(0), mOS(logFile) { mVals.reserve(10000); }

    void start( int gID )
    { gID = ~gID;
      mOS.write(reinterpret_cast<char const*>(&gID),sizeof(gID)); }
    void add( unsigned val )
    { if ( !val ) return;
      mPending += val;
      mOS.write(reinterpret_cast<char const*>(&val),sizeof(val)); }
    void flush()
    { if ( mPending ) {mVals.push_back(mPending); mPending = 0;}
      int val = 0;
      mOS.write(reinterpret_cast<char const*>(&val),sizeof(val)); }

    unsigned getN50()
    { flush();
      if ( mVals.empty() ) return 0;
      std::sort(mVals.begin(),mVals.end());
      return N50(mVals); }

private:
    vec<unsigned> mVals;
    unsigned mPending;
    std::ofstream mOS;
};

//set iSWAk to 501 to use hash-based aligner to avoid n^2 blow up
AssemblyError EvalAssembly( const ref_data& ref, const HyperBasevector& hb, 
        const vec<int>& inv, ostream& out, int iSWAk, int verbosity,
        N50Calculator* pN50Calc = nullptr );

#endif

