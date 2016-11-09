///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// A simplest context-free error model to score concensus sequence. 
//
// The required parameters are deletion, insertion, and substitution rate. The score
// is calcualted by sum over all the log likelyhood of the concensus to the observed 
// threads. 
//
// Currently all log likelyhood are converted to integers with two precision digits to
// speedup. The implementation if scoring method is very inefficient.
//

#ifndef CONSENSUS_SCORE_MODEL_H
#define CONSENSUS_SCORE_MODEL_H

#include "CoreTools.h"
#include "Vec.h"
#include "feudal/BaseVec.h"
#include "polymorphism/Edit.h"
#include "feudal/TrackingAllocator.h"

class ConsensusScoreModel {
    double del_rate, ins_rate, sub_rate; // save the input rate
    int U_Score; //= 4;    (the reference value are for del_rate = 3%, ins_rate = 0.2%, sub_rate = 0.8%)
    int D_Score; //= 350;  
    int I_Score; //= 760;
    int S_Score; //= 593;
    bool empty;  // if the model is uninitialized
    bool m_ignore_matching_score; // if ignore matching bases
    bool m_score_fast;            // if use fast scoring method

    public:
    ConsensusScoreModel() { empty = true; }
    ConsensusScoreModel( double del_rate, double ins_rate, double sub_rate, bool ignore_matching_score = true, bool score_fast = true ); 

    // Initialized the model with parameters
    // if ignore_matching_score is set, all the base-matching score will be ignored.
    void Init( double del_rate, double ins_rate, double sub_rate, bool ignore_matching_score=true, bool score_fast=true ) ;
    bool Empty() const { return empty; }
    double GetDelRate() const { return del_rate; }
    double GetInsRate() const { return ins_rate; }
    double GetSubRate() const { return sub_rate; }

    // The score of generating sequence b from sequence a. Lower score corresponds
    // to higher probability. (  = - logP * 100 )
    // Also optionally return the alignment strings.
    int Score( const BaseVec& a, const BaseVec& b, String* pt_align_a=0, String* pt_align_b=0 ) const ;
    int ScoreFast( const BaseVec& a, const BaseVec& b, String* pt_align_a=0, String* pt_align_b=0 ) const ;
    int ScoreFull( const BaseVec& a, const BaseVec& b, String* pt_align_a=0, String* pt_align_b=0 ) const ;
    // The score of generating the observed set of threads from sequence a. The
    // score is a summation over scores from a to each thread.
    // Also optionally return the alignment strings.
    int Score( const BaseVec& a, const vec<BaseVec>& threads, unsigned int min_vote = 0, 
            vec< pair<int,edit0> >* edits = 0 ) const ;

    // Calculate the probability of generating sequence b from sequence a
    // Multiple paths are considered
    double Probability( const BaseVec& a, const BaseVec& b) const;
 
    // Obtain the edits from the alignment strings of two sequences
    // !!! The new edits are ADDED to the existing vector.
    void GetEdits( const String & align_a, const String & align_b, vec< pair<int,edit0> > *loc_eidts ) const;

    void PrintScoreMatrix(  const BaseVec& a0, const BaseVec& b0 ) const ;

    friend class FastScorer;
};


class Row; // forward declaration
// A class  for error model-based alignment and scoring.  Similar to
// ConsensusScoreModel, but with speed improvement for long reads by caching
// the DP matrix and performing local re-calculation if only a small 
// edit is introduced into the founder read.
// Example:
// Suppose a substitution edit is introduced at position i in the consensus sequence a.
// Name the new sequence aa, we need to find the optimal alignment between aa and b.
// We know that all the optimal alignment between aa[0:i-1] and b[0:any j] is not affected.
// In other words, the rows [0:i] in the forward scoring matrix will be the same.
// Counting from backward, optimal alignments between aa[i+1,M-1] and b[any j:N-1] 
// are also the same.
// The optimal alignment between aa and b is the minimum of 
// optimal alignment aa[0:i] to b[0:j] ) plus alignment aa[i+1:M-1] to b[j+1,N-1] ),
// for j in [0,N). 
// Only need to calculate one row of the matrix for optimal alignment from aa[0:i] to b[0:j], which
// can be derived from aa[0:i-1] row. 
//
class FastScorer {
public:
    // constructor
    FastScorer( const ConsensusScoreModel& model, const BaseVec& read0, const BaseVec& read1 ) ;
    FastScorer( const ConsensusScoreModel& model_, const BaseVec& read0, const vec<BaseVec>& threads );
    // accessors
    int NThreads( ) const { return threads_.size(); }
    int GetScore( ) const;
    int GetScore(int ithread ) const;
    void PrintScoreMatrix( int i ) const;

    int ScoreEdit( const pair<int, edit0>& loc_edit ) ;
    int ScoreEdit1( int ithread, const pair<int, edit0> & loc_edit ) ;

    const BaseVec& FounderSeq() const { return read0_; }
    void GetLocEdit1( int ithread, vec< pair<int,edit0> > *p_loc_edits ) const;
    void GetLocEdits( vec< pair<int, edit0> > *p_loc_edits, vec<unsigned int> *p_counts, unsigned int mis_vote=0 ) const;
    void GetAlignString( int ithread, String *pt_align_a, String *pt_align_b ) const;

    // heuristics 
    static const int inf;                      // numerical limit of all scores.
    // command functions
    static void AlignToLocEdits( const String & align_a, const String & align_b, vec< pair<int,edit0> >* loc_edits );
    static BaseVec NewSeq( const BaseVec& read0,  const pair<int, edit0>& loc_edit );

private:
    void Init();
    void InitBackward();

private:
    const ConsensusScoreModel& model_;
    BaseVec read0_;
    vec<BaseVec> threads_;
    typedef vec< Row > ScoreMatrixType;
    vec< bool > valid_;
    vec< ScoreMatrixType > s_forwards_; // scoring matrix
    vec< ScoreMatrixType > s_backwards_; // scoring matrix from the sequence tail
    bool backward_set_;                  // if the backward scoring matrix is set
};

// A row in a matrix
class Row : private vec<int> {
public:
    typedef vec<int>::allocator_type allocator_type;
    int lb, ub;
    //-------------------------
    Row(int N, int val) : vec<int>(N, val), lb(0), ub(N) {};
    Row(int lbi, int ubi, int val ) : vec<int>(ubi-lbi, val), lb(lbi), ub(ubi) { };
    void SetCol(int j, int val)    { //ForceAssertLe( lb, j ); 
                                     //ForceAssertLt( j, ub );
                                     vec<int>::operator[](j-lb) = Min(val, FastScorer::inf); }
    int  GetCol(int j)   const     { if ( j < lb || j >=  ub ) return FastScorer::inf; 
                                     else return vec<int>::operator[]( j-lb ); }
    //int  operator[](int j) const   { return GetCol(j); }
    int  back() const              { return vec<int>::back(); }
};

#endif
