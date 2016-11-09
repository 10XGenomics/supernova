///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// A readstack represents a square matrix.  Each row of the matrix corresponds
// to a contiguous segment of a read, truncated only at the sides of the matrix.
// Each entry in the matrix is a base and a quality score, or undefined in the
// case one is off the end of the read.
//
// Undefined entries are shown as having a blank base ' ' and quality -1.

#ifndef READ_STACK_H
#define READ_STACK_H

#include "Basevector.h"
#include "CoreTools.h"
#include "PairsManager.h"
#include "Qualvector.h"
#include "feudal/PQVec.h"
#include "paths/long/FriendAligns.h"
#include "paths/long/MakeAlignments.h"

 typedef SerfVec<char>			StackBaseVec;
 typedef MasterVec<StackBaseVec > 	StackBaseVecVec;
 typedef SerfVec<char>			StackQualVec;
 typedef MasterVec<StackQualVec > 	StackQualVecVec;

class readstack {

     public:

     // ======================== CONSTRUCTORS ======================================

     // Constructor: initialize to undefined.

     readstack()=default;
     readstack( const int nrows, const int ncols ) { Initialize(nrows,ncols); }

     // Constructor from aligns [start,stop).  In the 'strict' case the columns
     // are defined by the bases of read id1, whereas in the 'right_extended' case
     // the column start at the leftmost base of read id1 but can extend beyond its
     // right end.
     enum con_type { strict, right_extended };
     readstack( const int64_t id1, Friends const& aligns,
          const int64_t start, const int64_t stop, con_type ctype,
          const vecbasevector& bases, const vecqualvector& quals, 
          const PairsManager& pairs, const Bool use_pairs = True )
     { Initialize(id1,aligns,start,stop,ctype,bases,quals,pairs,use_pairs); }

     void Initialize( const int nrows, const int ncols );
     void Initialize( const int64_t id1, Friends const& aligns,
          const int64_t start, const int64_t stop, con_type ctype,
          const vecbasevector& bases, const vecqualvector& quals, 
          const PairsManager& pairs, const Bool use_pairs = True );
     void Initialize( const int64_t id1, Friends const& aligns,
          const int64_t start, const int64_t stop, con_type ctype,
          const vecbasevector& bases, const VecPQVec& quals,
          const PairsManager& pairs, const Bool use_pairs = True );

     void AddRows( int n )
     {    int N = bases_.size( );
          bases_.resize( N + n );
          for ( int j = N; j < N + n; j++ )
               bases_[j].assign( Cols( ), ' ' );
          quals_.resize( N + n );
          for ( int j = N; j < N + n; j++ )
               quals_[j].assign( Cols( ), -1 );
          id_.resize( N + n, -1 );
          rc2_.resize( N + n, False );
          pid_.resize( N + n, -1 );
          pair_pos_.resize( N + n, -1 );
          offset_.resize( N + n, -1 );
          len_.resize( N + n, -1 );    }

     // RefStack: Given a chunk of reference sequence, build the stack consisting
     // of all reads having a perfect K-base overlap with it.

     template<int K> void RefStack( const basevector& ref,
          const vecbasevector& bases, const vecqualvector& quals,
          const PairsManager& pairs );

     int CheckSum( ) const
     {    int x = 0;
          for ( int i = 0; i < Rows( ); i++ )
          for ( int j = 0; j < Cols( ); j++ )
          {    x += Qual(i,j) + (int) Base(i,j);
               x *= 12345;    }
          return x;    }

     // ========================= ACCESSORS ========================================

     // Rows, Cols: return number of rows and columns.

     int Rows() const { return bases_.size( ); }
     int Cols() const { return cols_; }

     // Def: tell whether entry is defined.

     Bool Def( const int i, const int j ) const { return quals_[i][j] >= 0; }

     // Base: return given base, represented as 0, 1, 2, 3 or ' ' if undefined.

     char Base( const int i, const int j ) const { return bases_[i][j]; }

     char BaseVis( const int i, const int j ) const
     {    if ( bases_[i][j] == ' ' ) return ' ';
          else return as_base( bases_[i][j] );    }

     // Qual: return a given quality score (0-255), or -1 if undefined.

     int Qual( const int i, const int j ) const { return quals_[i][j]; }

     // SetBase, SetQual: set base or quality score.

     void SetBase( const int i, const int j, char b ) { bases_[i][j] = b; }
     void SetQual( const int i, const int j, int q ) { quals_[i][j] = q; }

     // Offset: implied start position of read or its rc relative to first column.

     int Offset( const int i ) const { return offset_[i]; }

     // SetOffset: set offset value.

     void SetOffset( const int i, const int offset ) { offset_[i] = offset; }

     // Len, SetLen: return or set read length.

     int Len( const int i ) const { return len_[i]; }
     void SetLen( const int i, const int len ) { len_[i] = len; }

     // Id, Rc2: return read id and orientation.

     int64_t Id( const int i ) const { return id_[i]; }
     Bool Rc2( const int i ) const { return rc2_[i]; }

     // Pid, PairPos: return pair id, and position (0 or 1) within pair.

     int64_t Pid( const int i ) const { return pid_[i]; }
     int PairPos( const int i ) const { return pair_pos_[i]; }

     // SetId, SetPid, SetRc2, SetPairPos: set these members.

     void SetId( const int i, const int64_t id ) { id_[i] = id; }
     void SetPid( const int i, const int64_t pid ) { pid_[i] = pid; }
     void SetRc2( const int i, const Bool rc2 ) { rc2_[i] = rc2; }
     void SetPairPos( const int i, const int pp ) { pair_pos_[i] = pp; }

     // Functions to return the entire matrix.

     const StackBaseVecVec& Bases( ) const { return bases_; }
     const StackQualVecVec& Quals( ) const { return quals_; }
     const vec<int64_t>& Id( ) const { return id_; }
     const vec<Bool>& Rc2( ) const { return rc2_; }
     const vec<int64_t>& Pid( ) const { return pid_; }
     const vec<int>& PairPos( ) const { return pair_pos_; }
     const vec<int>& Offset( ) const { return offset_; }
     const vec<int>& Len( ) const { return len_; }

     // Erase: remove the given rows.

     void Erase( const vec<Bool>& to_remove );

     // Add columns on right.

     void AddColumnsOnRight( const int n );

     // ========================= OPERATORS ========================================
     
     // SortByPid: place rows in a particular order.  First group them by pid.
     // Put pid1 first.  Sort the groups by their minimum offsets.  Sort within
     // a given group by orientation (fw, then rc), and then by offset.  Note
     // that for this, the offset would not be used unless a read appeared multiple 
     // times with the same orientation (but different offsets).
     // 
     // Make sure that first two rows are i1 then i2.

     void SortByPid( const int64_t pid1, const int i1, const int i2 );

     // Unique: delete adjacent and duplicate rows.  

     void Unique( );

     // Reverse: reverse complement the entire matrix.

     void Reverse( );

     // Merge: merge another readstack into *this, placing it at the given offset
     // relative to *this.  The entire readstack is placed after *this.

     void Merge( const readstack& s, const int offset );

     // AddToStack: add in reads at given offsets.

     void AddToStack( const vec< triple<int,int64_t,Bool> >& offset_id_rc2,
          const vecbasevector& bases, const vecqualvector& quals,
          const PairsManager& pairs );

     template< class VB, class VQ > void AddToStack2( 
          const vec< triple<int,int64_t,Bool> >& offset_id_rc2,
          VB& bases, VQ& quals );

     // AddPartners: for reads placed in the stack but whose partners are not,
     // recruit the partners based on perfect K-base overlaps.

     void AddPartners( const int K, const int top, const vecbasevector& bases,
          const vecqualvector& quals, const PairsManager& pairs );

     // Recruit: bring in reads that align to a read in the stack and overlap
     // the stack by at least K.

     void Recruit( const int K, const vec<simple_align_data>& aligns,
          const vec<int64_t>& id1_start, const vecbasevector& bases,
          const vecqualvector& quals, const PairsManager& pairs );

     // Trim: keep only the columns in [start,stop).

     void Trim( const int start, const int stop );

     // ========================= ALGORITHMS =======================================

     // ColumnConsensus1: return consensus base for column.

     char ColumnConsensus1( const int i ) const;

     // Consensus1, StrongConsensus1: algorithms for computing consensus.

     basevector Consensus1( ) const;
     void Consensus1( basevector& con, qualvector& conq, 
          const int qualcap = 50 ) const;
     void StrongConsensus1( basevector& con, qualvector& conq,
          const Bool raise_zero ) const;

     // Consensuses1: an algorithm for computing multiple consensuses.

     vec<basevector> Consensuses1( ) const;

     // Consensuses2: return all consensuses that can be traverse via K-base
     // perfect matches.

     vec<basevector> Consensuses2( const int K, const int top ) const;

     // Raise1: raise quality scores on row id based on agreement with it.

     void Raise1( const int id, const int rwindow = 11,
          const Bool require_unedited = False );

     // HighQualDiff: return the rows having a Qn diff with the first 'top' rows
     // (typically 1 or 2).

     void HighQualDiff( const int n, const int top, vec<Bool>& suspect ) const;

     void HighQualDiffWindow( vec<Bool>& suspect ) const;

     // CleanColumns: another criterion for killing rows.

     void CleanColumns( const int top, vec<Bool>& suspect ) const;

     // MotifDiff: find motifs and use them to define suspect rows.

     void MotifDiff( const int top, vec<Bool>& suspect ) const;

     // PairWeak1: First find the quality score sums for each column, using only
     // reads whose partner is placed in the stack.  If the quality score sum
     // for the winning base is at least 100, and > 10 times the sum for the 
     // second best, and the second best score is < 100, then return the
     // rows having such a q30 position.  (This needs to be explained better.)

     void PairWeak1( vec<Bool>& suspect ) const;

     // PairWeak2: Just like PairWeak1, but operates on two separate stacks.

     friend void PairWeak2( const readstack& s1, const readstack& s2,
          vec<Bool>& suspect1, vec<Bool>& suspect2 );

     // CorrectAll: correct all undisputed columns, and suggest a safe trim point.

     void CorrectAll( basevector& b, qualvector& q, int& trim_to,
          const Bool verbose = False ) const;

     // GetOffsets1: compute predicted offsets.

     friend vec<int> GetOffsets1( const readstack& stack1, const readstack& stack2, 
          const int verbosity, const int delta_mis,
          const basevector* con1 = NULL, const basevector* con2 = NULL );

     // FlagNoise: identify friends  having inadequate glue to the founder.

     void FlagNoise( vec<Bool>& ) const;

     // IdentifyShifters: identify friends that align much better to the founder
     // if shifted at a homopolymer.

     void IdentifyShifters( vec<Bool>& suspect ) const;
     
     // Defenestrate: look for windows that provide evidence of false friends, then
     // toss them through it.

     void Defenestrate( vec<Bool>& suspect ) const;

     // ========================== OUTPUT ==========================================

     // Print: pretty print the matrix, with w columns per window.  Each
     // window also has about 10 columns on the right containing source info.
     // There are two forms, both of which display consensus.  The first form
     // builds the consensus inside Print, whereas for the second form, the
     // consensus is provided as a second argument.  Multiple consensuses may
     // be provided, and are rendered in a collapsed form.

     void Print( ostream& out, const int w = 70, 
          const vec<String>& title = vec<String>( ) ) const;
     void Print( ostream& out, const vec<basevector>& cons, const int w = 70, 
          const vec<String>& title = vec<String>( ) ) const;
     void Print( ostream& out, const vec<vec<char>>& con, const int w = 70, 
          const vec<String>& title = vec<String>( ) ) const;

     // ========================= PRIVATE ==========================================


     private:

     int cols_;
     StackBaseVecVec bases_;
     StackQualVecVec quals_;
     vec<int64_t> id_;
     vec<Bool> rc2_;
     vec<int64_t> pid_;
     vec<int> pair_pos_;
     vec<int> offset_;
     vec<int> len_;

};

#endif
