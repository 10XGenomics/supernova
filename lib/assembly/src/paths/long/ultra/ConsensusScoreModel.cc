///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "paths/long/ultra/ConsensusScoreModel.h"
#include "CoreTools.h"
#include "Vec.h"
#include "VecUtilities.h"
#include <map>

ConsensusScoreModel::ConsensusScoreModel( double del_rate, double ins_rate, double sub_rate,
       bool ignore_matching_score, bool score_fast )
{
    Init(del_rate, ins_rate, sub_rate, ignore_matching_score, score_fast); 
    empty = false; 
}

// Initialized the model with parameters
void ConsensusScoreModel::Init( double del_rate_, double ins_rate_, double sub_rate_, 
        bool ignore_matching_score, bool score_fast ) 
{
    // save the input rate for future  reference
    del_rate = del_rate_;
    ins_rate = ins_rate_;
    sub_rate = sub_rate_;
    m_ignore_matching_score = ignore_matching_score;
    m_score_fast = score_fast;

    // If any of the error rates are zero, make them a small positive number
    // instead.  This is obviously not the correct approach, but it is better
    // than crashing.

    if ( del_rate == 0 ) del_rate = 0.000001;
    if ( ins_rate == 0 ) ins_rate = 0.000001;
    if ( sub_rate == 0 ) sub_rate = 0.000001;

    // sanity check
    ForceAssert( del_rate > 0 && del_rate < 1 );
    ForceAssert( ins_rate > 0 && ins_rate < 1 );
    ForceAssert( sub_rate > 0 && sub_rate < 1 );
    // convert log likelihood integers for fast calculation
    D_Score = - log( del_rate ) * 100;
    I_Score = - log( ins_rate/4 ) * 100; // 4 different choices for insertion
    S_Score = - log( sub_rate/3 ) * 100; // 3 different choices for substituttion
    U_Score = - log( 1 - del_rate - ins_rate - sub_rate ) * 100;
    // set U_Score = 0 to speed up
    if ( ignore_matching_score ) {
        D_Score -= U_Score;
        I_Score -= U_Score;
        S_Score -= U_Score;
        U_Score  = 0;
        ForceAssert( D_Score > 0 && I_Score > 0 && S_Score > 0 ); // The function will not work in this rate range
    }
    empty = false;
}



void ConsensusScoreModel::PrintScoreMatrix(  const BaseVec& a0, const BaseVec& b0 ) const {
    int M0 = a0.size(), N0 = b0.size();
    // trim off same head and tail.
    int len_head = 0;
    while ( len_head < M0 && len_head < N0 && a0[len_head] == b0[len_head] ) len_head++;
    int len_tail = 0;
    while ( M0 - len_tail - len_head > 0 && N0 - len_tail - len_head > 0
            && a0[M0-len_tail -1] == b0[N0- len_tail  -1] ) 
        len_tail++; 
    // match the middle parts of the two sequences
    BaseVec a( a0, len_head, M0 - len_head - len_tail );
    BaseVec b( b0, len_head, N0 - len_head - len_tail );
    int M = a.size(), N = b.size();
    const int infinity = 100000000;
    vec< vec<int> > s(M+1, vec<int>(N+1, infinity) );
    for ( int i = 0; i <= M; ++i ) s[i][0] = D_Score * i;
    for ( int j = 0; j <= N; ++j ) s[0][j] = I_Score * j;
    for ( int i = 1; i <= M; i++ ) {
        int j_lb = i * N * 0.8 / M - 1 ;
        if ( j_lb < 1 ) j_lb = 1;
        int j_ub = i * N * 1.2 / M  + 1;
        if ( j_ub > N ) j_ub = N;
        for ( int j = j_lb; j <= j_ub; j++ ) {
            int sub_score = s[i-1][j-1];
            if ( a[i-1] != b[j-1] ) sub_score += S_Score;
            int del_score = s[i-1][j] + D_Score;
            int ins_score = s[i][j-1] + I_Score;
            s[i][j] = std::min( {sub_score, del_score, ins_score} );
        }
    }
    cout << "len_head= " << len_head << endl;
    cout << "len_tail= " << len_tail << endl;
    for ( size_t x = 0; x < s.size(); ++x ) {
        cout << "row " << x << ": ";
        for( int j = 0; j <= N; j++ ) {
            String sf = ( s[x][j] >= infinity ? "***" : ToString(s[x][j]) );
            cout << sf << " ";
        }
        cout << endl;
    }
}


// The score of genrating sequence b from sequence a. Lower score corresponds
// to higher probability. (  = - logP * 100 )
// Also optionally return the alignment strings.
// 
// Matching sequences at head and tail are trimmed before actual alignment.
// Improve speed by ignoring unecessary DP calculatios
int ConsensusScoreModel::ScoreFast( const BaseVec& a0, const BaseVec& b0, String* pt_align_a, String* pt_align_b ) const 
{
    int M0 = a0.size(), N0 = b0.size();
    // trim off same head and tail.
    int len_head = 0;
    while ( len_head < M0 && len_head < N0 && a0[len_head] == b0[len_head] ) len_head++;
    int len_tail = 0;
    while ( M0 - len_tail - len_head > 0 && N0 - len_tail - len_head > 0
            && a0[M0-len_tail -1] == b0[N0- len_tail  -1] ) 
        len_tail++; 
    // match the middle parts of the two sequences
    BaseVec a( a0, len_head, M0 - len_head - len_tail );
    BaseVec b( b0, len_head, N0 - len_head - len_tail );
    int M = a.size(), N = b.size();
    const int infinity = 100000000;
    vec< vec<int> > s(M+1, vec<int>(N+1, infinity) );
    for ( int i = 0; i <= M; ++i ) s[i][0] = D_Score * i;
    for ( int j = 0; j <= N; ++j ) s[0][j] = I_Score * j;
    for ( int i = 1; i <= M; i++ ) {
        int j_lb = i * N * 0.8 / M - 1 ;
        if ( j_lb < 1 ) j_lb = 1;
        int j_ub = i * N * 1.2 / M  + 1;
        if ( j_ub > N ) j_ub = N;
        for ( int j = j_lb; j <= j_ub; j++ ) {
            int sub_score = s[i-1][j-1];
            if ( a[i-1] != b[j-1] ) sub_score += S_Score;
            int del_score = s[i-1][j] + D_Score;
            int ins_score = s[i][j-1] + I_Score;
            s[i][j] = std::min( {sub_score, del_score, ins_score} );
        }
    }
    //for ( size_t x = 0; x < s.size(); ++x ) {
    //    cout << "row " << x << ": ";
    //    copy( s[x].begin(), s[x].end(), ostream_iterator<int>(cout," ") );
    //    cout << endl;
    //}
    // alignment strings from trace back to the matrix
    if ( pt_align_a != 0 && pt_align_b != 0 ) {
        // sanity check
        if ( s[M][N] >= infinity ) return s[M][N];

        String & align_a = * pt_align_a; // shortcuts
        String & align_b = * pt_align_b;
        align_a.clear();
        align_b.clear();
        String aa = a.ToString();
        String bb = b.ToString();
        int i = M , j = N ;
        while ( i > 0 && j > 0 ){
            if ( s[i][j] == s[i-1][j-1] + ( a[i-1] == b[j-1] ? U_Score : S_Score ) ) {
                align_a = aa[i-1] + align_a;
                align_b = bb[j-1] + align_b;
                i--, j--;
            }
            else if ( s[i][j] == s[i-1][j] + D_Score ) {
                align_a = aa[i-1] + align_a;
                align_b = '-'  + align_b;
                i--;
            }
            else if ( s[i][j] == s[i][j-1] + I_Score ) {
                align_a = '-'  + align_a;
                align_b = bb[j-1] + align_b;
                j--;
            }
            else { 
                cout << "Error at position " << i << " " << j << endl; 
                cout << "Score= " << s[M][N] << endl;
                cout << "The two sequences are: " << endl; 
                cout << aa << endl;
                cout << bb << endl;
                cout << "The matrix: " << endl;
                for ( int ii = 0; ii <= i; ++ii ) {
                    cout << "row " << ii << ": ";
                    for ( int jj = 0; jj <= j; ++jj ) 
                        if ( s[ii][jj] == infinity )
                            cout << -1 << "\t";
                        else
                            cout << s[ii][jj] << "\t";
                    cout << endl;
                }
                return s[M][N];
                //CRD::exit(1); 
            }
        }
        while ( i > 0 ) {
            align_a = aa[i-1] + align_a;
            align_b = '-' + align_b;
            i--;
        }
        while ( j > 0 ) {
            align_b = bb[j-1] + align_b;
            align_a = '-' + align_a;
            j--;
        }
        // Add back the previously trimmed head and tail
        String head = "";
        String tail = "";
        if ( len_head > 0 )
            head = BaseVec(a0, 0, len_head).ToString();
        if ( len_tail > 0 )
            tail = BaseVec(a0, M0 - len_tail, len_tail ).ToString();
        align_a = head + align_a + tail;
        align_b = head + align_b + tail;
    }
    return s[M][N];
}


// The score of genrating sequence b from sequence a. Lower score corresponds
// to higher probability. (  = - logP * 100 )
// Also optionally return the alignment strings.
// 
// Matching sequences at head and tail are trimmed before actual alignment.
int ConsensusScoreModel::ScoreFull( const BaseVec& a0, const BaseVec& b0, String* pt_align_a, String* pt_align_b ) const 
{
    int M0 = a0.size(), N0 = b0.size();
    // trim off same head and tail only  if matching score is zereo
    int len_head = 0;
    int len_tail = 0;
    if ( U_Score == 0 ) { // 
        while ( len_head < M0 && len_head < N0 && a0[len_head] == b0[len_head] ) len_head++;
        while ( M0 - len_tail - len_head > 0 && N0 - len_tail - len_head > 0
                && a0[M0-len_tail -1] == b0[N0- len_tail  -1] ) 
            len_tail++; 
    }
    // match the middle parts of the two sequences
    BaseVec a( a0, len_head, M0 - len_head - len_tail );
    BaseVec b( b0, len_head, N0 - len_head - len_tail );
    int M = a.size(), N = b.size();
    vec< vec<int> > s(M+1, vec<int>(N+1,0) );
    for ( int i = 0; i <= M; ++i ) s[i][0] = D_Score * i;
    for ( int j = 0; j <= N; ++j ) s[0][j] = I_Score * j;
    for ( int i = 1; i <= M; i++ )
        for ( int j = 1; j <= N; j++ ) {
            int sub_score = s[i-1][j-1];
            if ( a[i-1] != b[j-1] ) sub_score += S_Score;
            else sub_score += U_Score;
            int del_score = s[i-1][j] + D_Score;
            int ins_score = s[i][j-1] + I_Score;
            s[i][j] = std::min( {sub_score, del_score, ins_score} );
        }
    // alignment strings from trace back to the matrix
    if ( pt_align_a != 0 && pt_align_b != 0 ) {
        String & align_a = * pt_align_a; // shortcuts
        String & align_b = * pt_align_b;
        align_a.clear();
        align_b.clear();
        String aa = a.ToString();
        String bb = b.ToString();
        int i = M , j = N ;
        while ( i > 0 && j > 0 ){
            int score = s[i][j];
            if ( score == s[i-1][j-1] + ( a[i-1] == b[j-1] ? U_Score : S_Score ) ) {
                align_a = aa[i-1] + align_a;
                align_b = bb[j-1] + align_b;
                i--, j--;
            }
            else if ( score == s[i-1][j] + D_Score ) {
                align_a = aa[i-1] + align_a;
                align_b = '-'  + align_b;
                i--;
            }
            else if ( score == s[i][j-1] + I_Score ) {
                align_a = '-'  + align_a;
                align_b = bb[j-1] + align_b;
                j--;
            }
            else { cout << "Error at position " << i << " " << j << endl; CRD::exit(1); }
        }
        while ( i > 0 ) {
            align_a = aa[i-1] + align_a;
            align_b = '-' + align_b;
            i--;
        }
        while ( j > 0 ) {
            align_b = bb[j-1] + align_b;
            align_a = '-' + align_a;
            j--;
        }
        // Add back the previously trimmed head and tail
        String head = "";
        String tail = "";
        if ( len_head > 0 )
            head = BaseVec(a0, 0, len_head).ToString();
        if ( len_tail > 0 )
            tail = BaseVec(a0, M0 - len_tail, len_tail ).ToString();
        align_a = head + align_a + tail;
        align_b = head + align_b + tail;
    }
    return s[M][N];
}

// Wrapper of different scoring methods.
int ConsensusScoreModel::Score( const BaseVec& a0, const BaseVec& b0, String* pt_align_a, String* pt_align_b ) const 
{
    if ( m_score_fast )
        return ScoreFast( a0, b0, pt_align_a, pt_align_b );
    else
        return ScoreFull( a0, b0, pt_align_a, pt_align_b );
}

// The score of genrating the observed set of threads from sequence a. The
// score is a summation over scores from a to each thread.
// Also optionally return the edits to t suggested by alignments to other threads.
int ConsensusScoreModel::Score( const BaseVec& a, const vec<BaseVec>& threads, unsigned int min_vote, 
        vec< pair<int,edit0> >* loc_edits ) const 
{
    int score = 0;
    if ( loc_edits == 0 ) { // no edits required
        for ( size_t i = 0; i < threads.size(); ++i ) 
            score += Score( a, threads[i] );
        return score;
    }
    loc_edits->clear();
    // otherwise calculate score and suggested edits
    map< pair<int,edit0>, vec<int> > vote_pool; // record the votes
    for ( size_t i = 0; i < threads.size(); ++i ) {
        String align_a, align_b;
        score += Score( a, threads[i], &align_a, &align_b );
        vec< pair<int,edit0> > edits;
        GetEdits( align_a, align_b, &edits);
        for( size_t j = 0; j < edits.size(); j++ )
            vote_pool[ edits[j] ].push_back(i);
    }
    vec <int> counts;
    for( map< pair<int,edit0>, vec<int> >::iterator it = vote_pool.begin();
            it != vote_pool.end(); it++ ) {
        if ( it->second.size() >= min_vote ) {
            loc_edits->push_back( it->first );
            counts.push_back( it->second.size() );
        }
    }
    ReverseSortSync( counts, *loc_edits );
    return score;
}

// The probability of generating sequence b from sequence a. 
// Degeneracy of generating from alternative alignments are considered.
// 
double ConsensusScoreModel::Probability( const BaseVec& a, const BaseVec& b ) const 
{
    int M = a.size(), N = b.size();
    vec< vec<double> > s(M+1, vec<double>(N+1,0) );
    double pd = del_rate;
    double pi = ins_rate / 4;
    double ps = sub_rate / 3;
    double pu = 1 - del_rate - ins_rate - sub_rate;
    s[0][0] = 1;
    for ( int i = 1; i <= M; ++i ) s[i][0] = pd * s[i-1][0];
    for ( int j = 1; j <= N; ++j ) s[0][j] = pi * s[0][j-1];
    for ( int i = 1; i <= M; i++ )
        for ( int j = 1; j <= N; j++ ) {
            double sub_score = s[i-1][j-1];
            if ( a[i-1] != b[j-1] ) sub_score *= ps;
            else sub_score *= pu;
            double del_score = s[i-1][j] * pd;
            double ins_score = s[i][j-1] * pi;
            s[i][j] = sub_score + del_score + ins_score;
        }
    return s[M][N];
}


// Obtain the edits from the alignment strings of two sequences. 
// !!! The new edits are ADDED to the existing vector.
void ConsensusScoreModel::GetEdits( const String & align_a, const String & align_b, vec< pair<int,edit0> >* loc_edits ) const 
{
    ForceAssertEq( align_a.size(), align_b.size() );
    int pos = 0;
    int actual_pos = 0; // actual editing position in sequence a, effectively a counter of none '-' 
    while ( pos < align_a.isize() ) {
        // Deal with insertion first. Preceed until firt none '-' position
        if ( align_a[pos] == '-' ) {
            size_t next_pos = pos;
            while ( next_pos < align_a.size() && align_a[next_pos] == '-' ) next_pos++;
            String seq = align_b.substr( pos, next_pos - pos );
            loc_edits->push( actual_pos, edit0( INSERTION, seq) );
            pos = next_pos;
            if ( pos >= align_a.isize() ) break;
        }
        // other cases
        if ( align_a[pos] != align_b[pos] ) 
            if ( align_b[pos] == '-' ) 
                loc_edits->push( actual_pos,  edit0( DELETION, 1 ) );
            else 
                loc_edits->push( actual_pos, edit0( SUBSTITUTION, align_b[pos] ) );
        pos++;
        actual_pos++;
    }
}

// Constructor
FastScorer::FastScorer( const ConsensusScoreModel& model, const BaseVec& read0, const vec<BaseVec>& threads ) : 
    model_(model), read0_(read0), threads_(threads)
{   Init();   }

// Constructor With only one thread, for debug purpose
FastScorer::FastScorer( const ConsensusScoreModel& model, const BaseVec& read0, const BaseVec& read1 ) : 
    model_(model), read0_(read0)
{   threads_.push_back( read1 );
    Init();   }

void FastScorer::Init() 
{   // heuristics
    const int MaxShift = 5;                   // To speedup, elements too far from diagnal will not be scored,
    const double ShiftRatio = 0.2;            // The max shift from diagnal is ( ShiftRatio * i + MaxShift );
    s_forwards_.assign( NThreads(), ScoreMatrixType() );
    valid_.assign( NThreads(), true );
    const int M = read0_.size();
    for( size_t ithread = 0; ithread < threads_.size(); ithread++ ) {
        const BaseVec& read1 = threads_[ithread];
        const int N = read1.size();
        if ( abs(M - N) > int( Max(M,N) * ShiftRatio ) ) { // size mis-smatch
            valid_[ithread] = false;
            continue;
        }
        ScoreMatrixType& s_forward = s_forwards_[ithread]; // best score for aligning read0[0:i] to read1[0:j]
        s_forward.clear();
        // calculate the boundary
        for ( int i = 0; i <= M; i++ ) {
            int j_avg = i * N / M; // the expected diagnal position
            int shift = ShiftRatio * Min( i, M-i ) * N / M + MaxShift;
            int j_lb = Max( 0, j_avg - shift );
            int j_ub = Min( N+1, j_avg + shift );
            ForceAssertLt( j_lb, j_ub );
            s_forward.push_back( Row( j_lb, j_ub, inf ) ); 
        }
        // first row
        for ( int j = s_forward[0].lb; j < s_forward[0].ub; ++j ) 
            s_forward[0].SetCol( j , model_.I_Score * j );
        for ( int i = 1; i <= M; i++ ) {
            int j_lb = s_forward[i].lb;
            int j_ub = s_forward[i].ub;
            if ( j_lb == 0 ) {
                s_forward[i].SetCol(0, model_.D_Score * i );
                ++j_lb;
            }
            for ( int j = j_lb; j < j_ub; j++ ) {
                int sub_score = s_forward[i-1].GetCol(j-1);
                if ( read0_[i-1] != read1[j-1] ) sub_score += model_.S_Score;
                int del_score = s_forward[i-1].GetCol(j)  + model_.D_Score;
                int ins_score = s_forward[i].GetCol( j-1 ) + model_.I_Score;
                s_forward[i].SetCol( j, std::min( {sub_score, del_score, ins_score} ) );
            }
        }
        // disable "bad" threads 
        if ( s_forward[M].GetCol(N) >= inf ) valid_[ithread] = false;
    }
    // Remember that we haven't don't the backward scoring matrix yet
    backward_set_ = false;
}

void FastScorer::InitBackward() 
{
    s_backwards_.resize( NThreads() );
    const int M = read0_.size();
    for( size_t ithread = 0; ithread < threads_.size(); ithread++ ) {
        if ( !valid_[ithread] ) continue;
        const BaseVec& read1 = threads_[ithread];
        const int N = read1.size();
        ScoreMatrixType& s_backward = s_backwards_[ithread];
        s_backward.clear();
        // calculate the boundary
        for ( int i = 0; i <= M; i++ ) {
            int j_lb = s_forwards_[ithread][i].lb;
            int j_ub = s_forwards_[ithread][i].ub;
            s_backward.push_back( Row( j_lb, j_ub, inf ) ); 
        }
        // last row
        for ( int j = s_backward[M].lb; j < s_backward[M].ub; ++j ) 
            s_backward[M].SetCol( j, model_.I_Score * (N-j) );
        // other rows
        for ( int i = M-1; i >=0; i-- ) {
            int j_lb = s_backward[i].lb;
            int j_ub = s_backward[i].ub;
            if ( j_ub == N+1 ) { // boundary
                s_backward[i].SetCol( N, model_.D_Score * (M-i) );
                --j_ub;
            }
            for ( int j = j_ub-1; j >= j_lb; j-- ) {
                int sub_score = s_backward[i+1].GetCol( j+1 );
                if ( read0_[i] != read1[j] ) sub_score += model_.S_Score;
                int del_score = s_backward[i+1].GetCol( j ) + model_.D_Score;
                int ins_score = s_backward[i].GetCol( j+1 ) + model_.I_Score;
                s_backward[i].SetCol( j, std::min( {sub_score, del_score, ins_score} ) );
            }
        }
    }
    backward_set_ = true;
}

int FastScorer::GetScore(int ithread) const {
    if ( ! valid_[ithread] ) return inf;
    return Min( inf, s_forwards_[ithread].back().back() );
}

int FastScorer::GetScore() const {
    int score = 0;
    for( size_t ithread = 0; ithread < threads_.size(); ithread++ ) 
        if ( valid_[ithread] ) score += GetScore(ithread);
    return score;
}

// print the scoring matrix for debugging
void FastScorer::PrintScoreMatrix(int ithread) const {
    const ScoreMatrixType& s_f = s_forwards_[ithread];
    const ScoreMatrixType& s_b = s_backwards_[ithread];
    int M = read0_.size();
    int N = threads_[ithread].size();
    for( int i = 0; i <= M; i++ ) {
        cout << "row " << i << ": ";
        for( int j = 0; j <= N; j++ ) {
            String sf = ( s_f[i].GetCol( j ) == inf ? "***" : ToString( s_f[i].GetCol( j )) );
            String sb = ( s_b[i].GetCol( j ) == inf ? "***" : ToString( s_b[i].GetCol( j )) );
            cout << sb << "/" << sf << " ";
        }
        cout << endl;
    }
}

// Re-score the consensus sequence if an edit is introduced 
int FastScorer::ScoreEdit( const pair<int, edit0> & loc_edit )  
{
    if ( ! backward_set_ ) InitBackward(); // lazy evaluation
    int score = 0;
    for( size_t ithread = 0; ithread < threads_.size(); ithread++ ) {
        if ( valid_[ithread] ) score += ScoreEdit1( ithread, loc_edit );
        if ( score >= inf )    return inf;
    }
    return score;
}

// Re-score the consensus sequence if an edit is introduced 
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
int FastScorer::ScoreEdit1( int ithread, const pair<int, edit0> & loc_edit )  
{
    if ( ! backward_set_ ) InitBackward();
    if ( ! valid_[ithread] ) return inf; 
    const BaseVec& read1 = threads_[ithread];
    const ScoreMatrixType& s_forward = s_forwards_[ithread];
    const ScoreMatrixType& s_backward = s_backwards_[ithread];
    int M = read0_.size(), N = read1.size();
    int pos = loc_edit.first;
    const edit0& e = loc_edit.second;
    int prev_row = -1, next_row = -1; 
    BaseVec read2;
    // Delete n sequences a position pos, forward scoring read0[0:pos-1], and backward scoring 
    // for read0[pos+n:M-1] are not affected. The inserted sequence is empty.
    if ( e.etype == DELETION ) {
        ForceAssertGe(pos, 0);
        ForceAssertLt(pos, M);
        prev_row = pos;
        next_row = pos + e.n;
    } 
    // Insert a sequence a position pos, forward scoring read0[0:pos-1], and backward scoring 
    // for read0[pos:M-1] are not affected. The inserted sequence is re-scored.
    else if ( e.etype == INSERTION ) {
        ForceAssertGe(pos, 0);
        ForceAssertLt(pos, M+1);
        prev_row = pos;
        next_row = pos;
        read2 = BaseVec( e.seq );
    }
    // substitute a sequence a position pos, forward scoring read0[0:pos-1], and backward scoring 
    // for read0[pos+1:M-1] are not affected. The modified sequence is re-scored.
    else if ( e.etype == SUBSTITUTION ) {
        ForceAssertGe(pos, 0);
        ForceAssertLt(pos, M);
        prev_row = pos;
        next_row = pos + 1;
        read2 = BaseVec( e.seq );
    } 
    else {
        cout << "Unknown edit type " << (int) e.etype << " at " << pos<< endl;
        CRD::exit(1);
    }
    int M2 = read2.size();
    vec< Row > s2( M2+1, s_forward[prev_row] );
    // i = 0 row is just a copy from s_forward[prev_row]
    for ( int i = 1; i <= M2; i++ ) {
        int j_lb = s2[i].lb;
        int j_ub = s2[i].ub;
        if ( j_lb == 0 ) {
            s2[i].SetCol( 0, (prev_row+i)* model_.D_Score ); 
            ++j_lb;
        }
        for ( int j = j_lb; j < j_ub; j++ ) {
            int sub_score = s2[i-1].GetCol( j-1 );
            if ( read2[i-1] != read1[j-1] ) sub_score += model_.S_Score;
            int del_score = s2[i-1].GetCol( j ) + model_.D_Score;
            int ins_score = s2[i].GetCol( j-1 ) + model_.I_Score;
            s2[i].SetCol(j, std::min( {sub_score, del_score, ins_score} ) );
        }
    }
    int j_lb = s_forward[prev_row].lb;
    int j_ub = s_forward[prev_row].ub;
    int min_score = inf;
    for( int j = j_lb; j < j_ub; j++ ) 
        min_score = Min( min_score,  s_backward[next_row].GetCol( j ) + s2.back().GetCol(j) );
    return min_score;
}

void FastScorer::GetAlignString( int ithread, String *pa, String *pb ) const
{
    String aa = read0_.ToString();
    String bb = threads_[ithread].ToString();
    const ScoreMatrixType& s = s_forwards_[ithread];
    String align_a( aa.size() + bb.size() );
    String align_b( aa.size() + bb.size() );
    String::iterator ia = align_a.end();
    String::iterator ib = align_b.end();
    int i = aa.size() , j = bb.size() ;
    while ( i > 0 && j > 0 ){
        int score = s[i].GetCol( j );
        if ( score >= inf ) { return; }
        if ( score == s[i-1].GetCol( j-1 ) + ( aa[i-1] == bb[j-1] ? model_.U_Score : model_.S_Score ) ) {
            *(--ia) = aa[--i];
            *(--ib) = bb[--j];
        }
        else if ( score == s[i-1].GetCol( j ) + model_.D_Score ) {
            *(--ia) = aa[--i];
            *(--ib) = '-';
        }
        else if ( score == s[i].GetCol( j-1 ) + model_.I_Score ) {
            *(--ia) = '-';
            *(--ib) = bb[--j];
        }
        else { cout << "Error at position " << i << " " << j << endl; CRD::exit(1); }
    }
    while ( i > 0 ) {
        *(--ia) = aa[--i];
        *(--ib) = '-'; 
    }
    while ( j > 0 ) {
        *(--ia) = '-';
        *(--ib) = bb[--j];
    }
    pa->append( ia, align_a.end() );
    pb->append( ib, align_b.end() );
}

void FastScorer::GetLocEdit1( int ithread, vec< pair<int,edit0> > *p_loc_edits ) const
{
    if ( ! valid_[ithread] ) return;
    String a1, a2;
    GetAlignString( ithread, &a1, &a2 );
    AlignToLocEdits( a1, a2, p_loc_edits );
}

void FastScorer::GetLocEdits( vec< pair<int, edit0> > *p_loc_edits, vec<unsigned int> * p_counts, unsigned int min_vote ) const
{
    ForceAssert( p_loc_edits != 0 );
    ForceAssert( p_counts != 0 );
    p_loc_edits->clear();
    p_counts->clear();
    for ( size_t ithread = 0; ithread < threads_.size(); ++ithread ) {
        if ( !valid_[ithread] ) continue;
        vec< pair<int, edit0> > new_loc_edits;
        GetLocEdit1( ithread, &new_loc_edits);
        // sanity check
        for ( size_t i = 0; i < new_loc_edits.size(); ++i ) {
            size_t p = new_loc_edits[i].first;
            const edit0&  e = new_loc_edits[i].second;
            if ( p >= read0_.size() && e.etype != INSERTION ) {
                cout << "Edit out of range " << endl;
                String str = "sub";
                if ( e.etype == INSERTION ) str = "ins";
                if ( e.etype == DELETION ) str = "del";
                cout << p << " " << str << " " 
                    << e.n << " " << e.seq << endl;
                cout << "read0= " << read0_.ToString() << endl;
                cout << "read1= " << threads_[i].ToString() << endl;
                CRD::exit(1);
            }
        }
        p_loc_edits->append( new_loc_edits );
    }
    UniqueSortAndCount( *p_loc_edits, *p_counts );
    ReverseSortSync( *p_counts, *p_loc_edits );
    // min_vote requirement
    size_t n = 0;
    while( n < p_counts->size() && (*p_counts)[n] >= min_vote ) n++;
    p_counts->resize(n);
    p_loc_edits->resize(n);
}


// ---------------------------- static constant and methods -----------------------------------------------
const int FastScorer::inf = 100*1000*1000; 

void FastScorer::AlignToLocEdits( const String & align_a, const String & align_b, vec< pair<int,edit0> >* loc_edits )
{
    ForceAssertEq( align_a.size(), align_b.size() );
    ForceAssert( loc_edits != 0 );
    loc_edits->clear();
    size_t pos = 0;
    int actual_pos = 0; // actual editing position in sequence a, effectively a counter of none '-' 
    while ( pos < align_a.size() ) {
        // Deal with insertion first. Preceed until firt none '-' position
        if ( align_a[pos] == '-' ) {
            size_t next_pos = pos;
            while ( next_pos < align_a.size() && align_a[next_pos] == '-' ) next_pos++;
            String seq = align_b.substr( pos, next_pos - pos );
            loc_edits->push_back( make_pair(actual_pos, edit0( INSERTION, seq)) );
            pos = next_pos;
            if ( pos >= align_a.size() ) break;
        }
        // other cases
        if ( align_a[pos] != align_b[pos] ) 
            if ( align_b[pos] == '-' ) 
                loc_edits->push( make_pair(actual_pos,  edit0( DELETION, 1 )) );
            else 
                loc_edits->push( make_pair(actual_pos, edit0( SUBSTITUTION, align_b[pos] )) );
        pos++;
        actual_pos++;
    }
}

BaseVec FastScorer::NewSeq( const BaseVec& t, const pair<int, edit0>& loc_edit ) 
{
    size_t p = loc_edit.first;
    const edit0& e = loc_edit.second;
    BaseVec t2;
    if ( e.etype == INSERTION ) {
        // note that p == t.size() can happen if insertion happens at the end of the sequence
        ForceAssertLe( p, t.size() );
        BaseVec b1( t, 0, p );
        BaseVec b2( e.seq );
        BaseVec b3;
        if ( (size_t) p < t.size() )
            b3.SetToSubOf( t, p, t.size( ) - p );
        t2 = Cat( b1, b2, b3 );    
    }
    else if ( e.etype == DELETION ) {
        ForceAssertLe( p + e.n , t.size( ) );
        BaseVec b1( t, 0, p );
        BaseVec b2( t, p + e.n, t.size( ) - ( p + e.n ) );
        t2 = Cat( b1, b2 );    
    }
    else if ( e.etype == SUBSTITUTION ) {
        ForceAssertLt( p, t.size() );
        t2 = t;
        t2.Set( p, as_char( e.seq[0] ) );    
    }
    else {
        cout << "Unknown edit type " << (int) e.etype << " at " << p << endl;
        CRD::exit(1);
    }
    return t2;
}
