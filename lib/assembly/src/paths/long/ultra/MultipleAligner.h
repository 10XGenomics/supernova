///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file MultipleAligner.h
 * \author tsharpe
 * \date May 7, 2012
 *
 * \brief
 */
#ifndef PATHS_LONG_MULTIPLEALIGNER_H_
#define PATHS_LONG_MULTIPLEALIGNER_H_

#include "Basevector.h"
#include "Charvector.h"
#include "system/Assert.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <iterator>
#include <numeric>
#include <ostream>
#include <vector>

/// Adds extra '-' code to allow gapped bvecs.
class gvec
{
public:
    gvec() {}

    gvec( bvec const& bv ) : mVals(bv.begin(),bv.end()) {}

    template <class Itr> // Itr returns base codes or gvec codes
    gvec( Itr const& beg, Itr const& end ) : mVals(beg,end) {}

    gvec( gvec const& gv ) : mVals(gv.mVals) {}

    gvec& operator=( gvec const& gv ) { mVals = gv.mVals; return *this; }

    typedef std::vector<unsigned char> container;
    typedef container::value_type value_type;
    typedef container::size_type size_type;
    typedef container::difference_type difference_type;
    typedef container::const_iterator const_iterator;
    typedef container::const_reverse_iterator const_reverse_iterator;

    const_iterator begin() const { return mVals.begin(); }
    const_iterator cbegin() const { return mVals.begin(); }
    const_iterator end() const { return mVals.end(); }
    const_iterator cend() const { return mVals.end(); }
    const_reverse_iterator rbegin() const { return mVals.rbegin(); }
    const_reverse_iterator crbegin() const { return mVals.rbegin(); }
    const_reverse_iterator rend() const { return mVals.rend(); }
    const_reverse_iterator crend() const { return mVals.rend(); }

    // note: there is no non-const operator[], because we can't let people
    // shove arbitrary values into slots.  instead, you have to use set so that
    // we can check for legal values.

    value_type const& operator[]( size_type idx ) const
    { return mVals[idx]; }

    value_type at( size_type idx ) const
    { ForceAssertLt(idx,mVals.size()); return mVals[idx]; }

    gvec& reserve( size_type nnn ) { mVals.reserve(nnn); return *this; }

    gvec& set( size_type idx, value_type val )
    { AssertLt(idx,mVals.size()); ForceAssertLt(val,sizeof(gCharSet));
      mVals[idx] = val; return *this; }

    gvec& push_back( value_type val )
    { ForceAssertLt(val,sizeof(gCharSet));
      mVals.push_back(val); return *this; }

    gvec& insert( size_type idx, value_type val )
    { AssertLe(idx,mVals.size()); ForceAssertLt(val,sizeof(gCharSet));
      mVals.insert(mVals.begin()+idx,val); return *this; }

    gvec& reverse()
    { std::reverse(mVals.begin(),mVals.end()); return *this; }

    size_type size() const { return mVals.size(); }

    static char charFor( gvec::value_type val ) { return gCharSet[val]; }

    friend bool operator==( gvec const& gv1, gvec const& gv2 )
    { return gv1.mVals == gv2.mVals; }
    friend bool operator!=( gvec const& gv1, gvec const& gv2 )
    { return gv1.mVals != gv2.mVals; }
    friend bool operator<( gvec const& gv1, gvec const& gv2 )
    { return gv1.mVals < gv2.mVals; }

    friend std::ostream& operator<<( std::ostream& os, gvec const& gv )
    { std::ostream_iterator<char> osi(os);
      std::transform(gv.begin(),gv.end(),osi,Val2CharMapper());
      return os; }

    struct Val2CharMapper
    { char operator()( value_type val ) const { return gCharSet[val]; } };

    enum VAL { A, C, G, T, GAP };
    static gvec::value_type const NVALS = 5;

private:
    container mVals;

    static char const gCharSet[NVALS];
};


class Scorer
{
public:
    class Builder
    {
    public:
        Builder() { memset(mBins,0,sizeof(mBins)); }

        size_t count( gvec::value_type testVal,
                        gvec::value_type consVal ) const
        { AssertLt(testVal,NVALS); AssertLt(consVal,NVALS);
          return mBins[NVALS*consVal+testVal]; }

        size_t consTot( gvec::value_type consVal ) const
        { AssertLt(consVal,NVALS);
          size_t const* ppp = mBins+NVALS*consVal;
          return std::accumulate(ppp,ppp+NVALS,0ul); }

        void match( gvec::value_type testVal, gvec::value_type consVal )
        { AssertLt(testVal,NVALS); AssertLt(consVal,NVALS);
          mBins[NVALS*consVal+testVal] += 1; }
        void insertion( gvec::value_type testVal )
        { match(testVal,gvec::GAP); }
        void deletion( gvec::value_type consVal )
        { match(gvec::GAP,consVal); }

    private:
        static gvec::value_type const NVALS = gvec::NVALS;
        size_t mBins[NVALS*NVALS];
    };

    // make an edit-distance scorer
    Scorer()
    { std::fill(mPhred,mPhred+NVALS*NVALS,1u);
      for ( size_t idx = 0; idx != NVALS; ++idx )
          mPhred[NVALS*idx+idx] = 0; }

    Scorer( double substRate, double delRate, double insRate )
    { ForceAssertGe(substRate,0.);
      ForceAssertLe(substRate,1.);
      ForceAssertGe(delRate,0.);
      ForceAssertLe(delRate,1.);
      ForceAssertGe(insRate,0.);
      ForceAssertLe(insRate,1.);
      double mult = 10.;
      if ( insRate > 0. )
        while ( !phred(1.-insRate,mult) ) mult *= 10;
      if ( delRate > 0. || substRate > 0. )
        while ( !phred(1.-delRate-substRate,mult) ) mult *= 10;
      unsigned delVal = phred(delRate,mult);
      unsigned substVal = phred(substRate/3.,mult);
      unsigned matchVal = phred(1.-delRate-substRate,mult);
      unsigned insVal = phred(insRate/4.,mult);
      unsigned noInsVal = phred(1.-insRate,mult);
      unsigned* pPhred = mPhred;
      for ( gvec::value_type consVal = 0; consVal != NVALS; ++consVal )
      { if ( consVal == gvec::GAP )
          for ( gvec::value_type testVal = 0; testVal != NVALS; ++testVal )
            if ( testVal == gvec::GAP ) *pPhred++ = noInsVal;
            else *pPhred++ = insVal;
        else
          for ( gvec::value_type testVal = 0; testVal != NVALS; ++testVal )
            if ( testVal == consVal ) *pPhred++ = matchVal;
            else if ( testVal == gvec::GAP ) *pPhred++ = delVal;
            else *pPhred++ = substVal; } }

    Scorer( Builder const& builder )
    { unsigned* ppp = mPhred;
      for ( gvec::value_type consVal = 0; consVal != NVALS; ++consVal )
      { double tot = builder.consTot(consVal);
        ForceAssertGt(tot,0.);
        double mult = 1.;
        do
        { mult *= 10.;
          for ( gvec::value_type testVal = 0; testVal != NVALS; ++testVal )
          { unsigned score = phred(builder.count(testVal,consVal)/tot,mult);
            if ( !score ) break;
            *ppp++ = score; }
        } while ( ppp != mPhred+NVALS*NVALS ); } }

    unsigned operator()( gvec::value_type testVal,
                                gvec::value_type consVal ) const
    { AssertLt(testVal,NVALS);  AssertLt(consVal,NVALS);
      return mPhred[NVALS*consVal+testVal]; }

private:
    static unsigned phred( double rate, double mult )
    { return rate > 0. ? std::min(6*mult,-mult*log10(rate))+.5 : 6*mult; }

    static gvec::value_type const NVALS = gvec::NVALS;
    unsigned mPhred[NVALS*NVALS];
};


class MultipleAligner
{
public:
    MultipleAligner( Scorer const& scorer, gvec const& consensus )
    : mNW(scorer)
    { mScores.reserve(2*consensus.size());
      typedef gvec::const_iterator Itr;
      Itr end(consensus.end());
      for ( Itr itr(consensus.begin()); itr != end; ++itr )
        mScores.push_back(ConsScore(*itr)); }

    MultipleAligner( Scorer const& scorer, size_t widthEstimate = 0 )
    : mNW(scorer) { mScores.reserve(widthEstimate); }

    template <class Itr> // Itr returns read-like objects
    void addReads( Itr itr, Itr const& end )
    { while ( itr != end ) { addRead(itr->begin(),itr->end()); ++itr; } }

    template <class Itr> // Itr returns base codes
    void addRead( Itr const& beg, Itr const& end )
    { mNW.align(mScores,beg,end);
      mNW.tracebackAdjustScores(mScores,beg); }

    template <class Itr> // Itr returns base codes
    gvec alignRead( Itr const& beg, Itr const& end )
    { mNW.align(mScores,beg,end);
      return mNW.traceback(mScores,beg); }

    Scorer const& getScorer() const { return mNW.getScorer(); }

    bvec getConsensus() const
    { bvec result; result.reserve(mScores.size());
      for ( CSItr itr(mScores.begin()), end(mScores.end()); itr != end; ++itr )
        if ( itr->bestVal() != gvec::GAP )
          result.push_back(itr->bestVal());
      return result; }

    size_t getConsensusLen() const { return mScores.size(); }


    class ConsScore
    {
    public:
        ConsScore( gvec::value_type val )
        : mBestVal(val)
        { memset(mCounts,0,sizeof(mCounts)); }

        ConsScore( gvec::value_type val, size_t nGaps, Scorer const& scorer )
        : mBestVal(gvec::GAP)
        { memset(mCounts,0,sizeof(mCounts));
          mCounts[gvec::GAP] = nGaps;
          add(val,scorer); }

        gvec::value_type bestVal() const { return mBestVal; }

        size_t nVals() const
        { return std::accumulate(mCounts,mCounts+NVALS,0ul); }

        unsigned getCount( gvec::value_type code ) const
        { AssertLt(code,NVALS); return mCounts[code]; }

        unsigned getBestCount() const { return mCounts[mBestVal]; }

        unsigned getErrCount() const{ return nVals()-mCounts[mBestVal]; }

        void add( gvec::value_type val, Scorer const& scorer )
        { mCounts[val] += 1;
          if ( val != mBestVal )
          { if ( eval(val,scorer) < eval(mBestVal,scorer) )
              mBestVal = val; } }

        void print( ostream& os, Scorer const& scorer ) const
        { os << gvec::charFor(mBestVal); char sep = '[';
          for ( gvec::value_type consVal = 0; consVal!=NVALS; ++consVal )
          { os << sep; sep = ','; os << eval(consVal,scorer); }
          os << ']'; }

    private:
        unsigned eval( gvec::value_type consVal, Scorer const& scorer )const
        { unsigned result = 0;
          for ( gvec::value_type testVal = 0; testVal!=NVALS; ++testVal )
            result += mCounts[testVal]*scorer(testVal,consVal);
          return result; }

        static gvec::value_type const NVALS = gvec::NVALS;
        gvec::value_type mBestVal;
        unsigned mCounts[NVALS];
    };
    typedef std::vector<ConsScore> VecConsScore;
    typedef VecConsScore::const_iterator CSItr;

    VecConsScore const& getScores() const { return mScores; }

    template <class Itr> // Itr returns base codes
    gvec getAlignment( Itr const& beg, Itr const& end ) const
    { mNW.align(mScores,beg,end);
      return mNW.tracebackGVec(beg); }

    template <class Itr> // Itr returns read-like objects
    void printMultipleAlignment( std::ostream& os,
                                        Itr itr, Itr const& end ) const
    { printGappedConsensus(os);
      while ( itr != end )
      { printAlignment(os,itr->begin(),itr->end()); ++itr; } }

    friend ostream& operator<<( ostream& os, MultipleAligner const& ma )
    { VecConsScore const& scores = ma.getScores();
      for ( CSItr itr(scores.begin()), end(scores.end()); itr != end; ++itr )
        itr->print(os,ma.getScorer());
      return os; }

private:
    class NeedlemanWunscher
    {
    public:
        NeedlemanWunscher( Scorer const& scorer )
        : mScorer(scorer) {}

        Scorer const& getScorer() const { return mScorer; }

        template <class Itr>
        unsigned align( VecConsScore const& scores, Itr itr, Itr const& end )
        { using std::distance;
          mBackPtrs.clear().reserve(distance(itr,end)+1);
          mRow[0].clear(); mRow[0].reserve(2*scores.size());
          mRow[1].clear(); mRow[1].reserve(2*scores.size());
          firstRow(scores);
          size_t whichRow = 0;
          for ( ; itr != end; ++itr )
            nextRow(scores,whichRow ^= 1,*itr);
          return mRow[whichRow].back(); }

        template <class Itr> // Itr returns base codes
        void tracebackAdjustScores(VecConsScore& scores, Itr const& read) const
        { size_t consIdx = scores.size();
          size_t testIdx = mBackPtrs.size()-1;
          size_t nReads = scores.size() ? scores[0].nVals() : 0;
          while ( true )
          { switch ( mBackPtrs[testIdx][consIdx] )
            {case DONE:
              AssertEq(testIdx,0ul); AssertEq(consIdx,0ul);
              return;
             case MATCH:
              scores[--consIdx].add(read[--testIdx],mScorer);
              break;
             case DEL:
              scores[--consIdx].add(gvec::GAP,mScorer);
              break;
             case INS:
              ConsScore cs(read[--testIdx],nReads,mScorer);
              scores.insert(scores.begin()+consIdx,cs);
              break; } } }

        template <class Itr> // Itr returns base codes
        gvec traceback( VecConsScore& scores, Itr const& read ) const
        { size_t consIdx = mBackPtrs[0].size()-1;
          size_t testIdx = mBackPtrs.size()-1;
          size_t nReads = scores.size() ? scores[0].nVals() : 0;
          gvec result; result.reserve(consIdx);
          bool done = false;
          while ( !done )
          { bvec::value_type baseCode;
            switch ( mBackPtrs[testIdx][consIdx] )
            {case DONE:
              AssertEq(testIdx,0ul); AssertEq(consIdx,0ul);
              result.reverse();
              done = true;
              break;
             case MATCH:
              baseCode = read[--testIdx];
              scores[--consIdx].add(baseCode,mScorer);
              result.push_back(baseCode);
              break;
             case DEL:
              scores[--consIdx].add(gvec::GAP,mScorer);
              result.push_back(gvec::GAP);
              break;
             case INS:
              baseCode = read[--testIdx];
              ConsScore cs(baseCode,nReads,mScorer);
              scores.insert(scores.begin()+consIdx,cs);
              result.push_back(baseCode);
              break; } }
          return result; }

        template <class Itr> // Itr returns base codes
        gvec tracebackGVec( Itr const& read ) const
        { size_t consIdx = mBackPtrs[0].size()-1;
          size_t testIdx = mBackPtrs.size()-1;
          gvec result; result.reserve(consIdx);
          bool done = false;
          while ( !done )
          { switch ( mBackPtrs[testIdx][consIdx] )
            {case DONE:
              AssertEq(testIdx,0ul); AssertEq(consIdx,0ul);
              result.reverse();
              done = true;
              break;
             case MATCH:
              result.push_back(read[--testIdx]);
              consIdx -= 1;
              break;
             case DEL:
              result.push_back(gvec::GAP);
              consIdx -= 1;
              break;
             case INS:
              result.push_back(read[--testIdx]);
              break; } }
          return result; }

    private:
        enum BACKPTR { DONE, INS, DEL, MATCH };

        typedef std::vector<unsigned> DPRow;

        void firstRow( VecConsScore const& scores )
        { CharVec& vecBackPtr = mBackPtrs.resize(1).back();
          vecBackPtr.resize(scores.size()+1,DEL);
          DPRow& row = mRow[0];
          row.clear();
          row.push_back(0u);
          vecBackPtr[0] = DONE;
          for ( CSItr itr(scores.begin()), end(scores.end()); itr!=end; ++itr )
            row.push_back(row.back()+mScorer(gvec::GAP,itr->bestVal())); }

        void nextRow( VecConsScore const& scores, size_t whichRow,
                            gvec::value_type testVal )
        { CharVec& vecBackPtr = mBackPtrs.resize(mBackPtrs.size()+1).back();
          vecBackPtr.reserve(scores.size()+1);
          DPRow& row = mRow[whichRow];
          row.clear();
          DPRow::const_iterator pPrevRow = mRow[whichRow ^ 1].begin();
          row.push_back(*pPrevRow+mScorer(testVal,gvec::GAP));
          vecBackPtr.push_back(INS);
          for ( CSItr itr(scores.begin()), end(scores.end()); itr!=end; ++itr )
          { gvec::value_type consVal = itr->bestVal();
            unsigned score = *pPrevRow+mScorer(testVal,consVal);
            BACKPTR backPtr = MATCH;
            unsigned score2 = row.back()+mScorer(gvec::GAP,consVal);
            if ( score2 < score )
            { score = score2; backPtr = DEL; }
            score2 = *++pPrevRow+mScorer(testVal,gvec::GAP);
            if ( score2 < score )
            { score = score2; backPtr = INS; }
            row.push_back(score);
            vecBackPtr.push_back(backPtr); } }

        Scorer const& mScorer;
        VecCharVec mBackPtrs; // a double-array of backptrs
        DPRow mRow[2]; // a pair of rows for keeping the DP scores
    };

    void printGappedConsensus( std::ostream& os ) const
    { for ( CSItr itr(mScores.begin()), end(mScores.end()); itr != end; ++itr )
        os << gvec::charFor(itr->bestVal());
      os << '\n'; }

    template <class Itr> // Itr returns base codes
    void printAlignment( std::ostream& os,
                                Itr const& beg, Itr const& end ) const
    { unsigned score = mNW.align(mScores,beg,end);
      os << mNW.tracebackGVec(beg) << ' ' << score << '\n'; }

    mutable NeedlemanWunscher mNW;
    VecConsScore mScores; // a score for each consensus base
};

#endif /* PATHS_LONG_MULTIPLEALIGNER_H_ */
