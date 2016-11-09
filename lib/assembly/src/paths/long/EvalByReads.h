///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef EVAL_BY_READS_H
#define EVAL_BY_READS_H

#include "Basevector.h"
#include "CoreTools.h"
#include "IntPairVec.h"
#include "Qualvector.h"
#include "paths/HyperBasevector.h"

// A read_place defines a gap-free alignment of a read to a HyperBasevector.
// In the forward case, the read is placed, starting at a given base P( ) on a given
// edge E(0) of the HyperBasevctor.  In the reverse case, the read is instead placed
// on the reverse complement of the HyperBasevector.

class read_place {

     public:

     read_place( ) { }

     read_place( const vec<int>& e, const int p, const Bool fw, const int qsum )
          : e_(e), p_(p), fw_(fw), qsum_(qsum) { }

     int N( ) const { return e_.size( ); }
     int E( const int i ) const { return e_[i]; }
     const vec<int>& E( ) const { return e_; }
     vec<int>& EMutable( ) { return e_; }
     int P( ) const { return p_; }
     void AddEdge( const int e ) { e_.push_back(e); }
     Bool Fw( ) const { return fw_; }
     int Qsum( ) const { return qsum_; }
     void SetQsum( int qsum ) { qsum_ = qsum; }

     // Reverse: define the placement one would get if the read were instead
     // reverse complemented.

     void Reverse( const basevector& b, const HyperBasevector& hb )
     {    p_ = hb.EdgeLengthBases( E(0) ) - ( P( ) + b.isize( ) );
          for ( int j = 1; j < N( ); j++ )
               p_ += hb.EdgeLengthKmers( E(j) );
          EMutable( ).ReverseMe( );
          fw_ = !fw_;    }

     friend ostream& operator<<( ostream& out, const read_place& x )
     {    out << ( x.Fw( ) ? "fw" : "rc" ) << " " << x.E(0) << "." << x.P( );
          for ( int j = 1; j < x.N( ); j++ )
               out << "," << x.E(j);
          return out << " [" << setiosflags(ios::fixed) << setprecision(1)
               << x.Qsum( ) / 1000.0 << resetiosflags(ios::fixed) << "]";    }

     void ComputeQsum( const basevector& b, const qualvector& q, 
          const HyperBasevector& hb, const int min_qual )
     {    int ei = 0, pos = P( );
          qsum_ = 0;
          for ( int l = 0; l < b.isize( ); l++ )
          {    if ( b[l] != hb.EdgeObject( E(ei) )[pos] )
               {    if ( q[l] >= min_qual ) qsum_ += q[l] * 1000;
                    else qsum_ += q[l];    }
               pos++;
               if ( pos == hb.EdgeObject( E(ei) ).isize( ) )
               {    ei++;
                    if ( ei == N( ) ) break;
                    pos = hb.K( ) - 1;    }    }    }
     
     void ComputeQsum125( const basevector& b, const qualvector& q, 
          const HyperBasevector& hb, const int min_qual )
     {    int ei = 0, pos = P( );
          qsum_ = 0;
          for ( int l = 0; l < b.isize( ); l++ )
          {    if ( b[l] != hb.EdgeObject( E(ei) )[pos] )
               {    if ( q[l] >= min_qual ) qsum_ += q[l] * 1000;
                    else qsum_ += 1250;    }
               pos++;
               if ( pos == hb.EdgeObject( E(ei) ).isize( ) )
               {    ei++;
                    if ( ei == N( ) ) break;
                    pos = hb.K( ) - 1;    }    }    }

     int MaxQ( const basevector& b, const qualvector& q, 
          const HyperBasevector& hb ) const
     {    int maxq = 0;
          int ei = 0, pos = P( );
          for ( int l = 0; l < b.isize( ); l++ )
          {    if ( b[l] != hb.EdgeObject( E(ei) )[pos] ) 
                    maxq = Max( maxq, (int) q[l] );
               pos++;
               if ( pos == hb.EdgeObject( E(ei) ).isize( ) )
               {    ei++;
                    if ( ei == N( ) ) break;
                    pos = hb.K( ) - 1;    }    }
          return maxq;    }

     void AddEdge( const int e, const basevector& b, const qualvector& q, 
          const HyperBasevector& hb, const int min_qual )
     {    int ei = 0, l = 0;
          for ( ; ei < N( ); ei++ )
          {    l += hb.EdgeObject( E(ei) ).size( );
               if ( ei == 0 ) l -= P( );
               else l -= ( hb.K( ) - 1 );    }
          e_.push_back(e);
          int pos = hb.K( ) - 1;
          for ( ; l < b.isize( ); l++ )
          {    if ( b[l] != hb.EdgeObject( E(ei) )[pos] )
               {    if ( q[l] >= min_qual ) qsum_ += q[l] * 1000;
                    else qsum_ += q[l];    }
               pos++;
               if ( pos == hb.EdgeObject( E(ei) ).isize( ) ) break;    }    }
     
     void AddEdge125( const int e, const basevector& b, const qualvector& q, 
          const HyperBasevector& hb, const int min_qual )
     {    int ei = 0, l = 0;
          for ( ; ei < N( ); ei++ )
          {    l += hb.EdgeObject( E(ei) ).size( );
               if ( ei == 0 ) l -= P( );
               else l -= ( hb.K( ) - 1 );    }
          e_.push_back(e);
          int pos = hb.K( ) - 1;
          for ( ; l < b.isize( ); l++ )
          {    if ( b[l] != hb.EdgeObject( E(ei) )[pos] )
               {    if ( q[l] >= min_qual ) qsum_ += q[l] * 1000;
                    else qsum_ += 1250;    }
               pos++;
               if ( pos == hb.EdgeObject( E(ei) ).isize( ) ) break;    }    }

     // FindReadShift: Find the relative shift of the read with respect 
     // to the given edge. Return multiple shifts if the edge repeats.

     void FindReadShift(int edge, const HyperBasevector& hb, vec<int>& shifts)
     {    shifts.clear();
          int shift = P();
          if (E(0) == edge) shifts.push_back(shift);
          for ( int j = 1; j < N( ); j++ ) {
               shift -= hb.EdgeLengthKmers(E(j-1));
               if (E(j) == edge) shifts.push_back(shift);    }    }

     void writeBinary( BinaryWriter& writer ) const;
     void readBinary( BinaryReader& reader );

     private:

     vec<int> e_;   // edge sequence
     int p_;        // start position on first edge
     Bool fw_;      // if placement of read rather than its reverse complement
     int qsum_;     // sum of quality scores at mismatch

};

template<> struct Serializability<read_place>
{ typedef SelfSerializable type; };

// WARNING - extra args at end of FindPlaces are not necessarily handled 
// correctly.  See commented out code at end of FindPlaces in the .cc file.

void FindPlaces( const basevector& b, const qualvector& q, const int n, 
     const HyperBasevector& hb_fw, const HyperBasevector& hb_rc, 
     const vec<int>& to_right_fw, const vec<int>& to_right_rc, 
     const VecIntPairVec& locs_fw,
     const VecIntPairVec& locs_rc,
     vec<read_place>& places, int& qual_sum, const int min_qual = 3, 
     const double prox = 0, const int max_diff = 1000000000 );

uint64_t SafeFindPlaces( const basevector& b, const qualvector& q, const int n,
     const HyperBasevector& hb_fw, const HyperBasevector& hb_rc,
     const vec<int>& to_right_fw, const vec<int>& to_right_rc,
     const VecIntPairVec& locs_fw,
     const VecIntPairVec& locs_rc,
     vec<read_place>& places, int& qual_sum, uint64_t maxNSteps, const int min_qual = 3,
     const double prox = 0, const int max_diff = 1000000000 );

void EvalByReads( const HyperBasevector& hb_A, const HyperBasevector& hb_R,
     const vecbasevector& bases, vecqualvector quals, int& assembly_count,
     int& reference_count, const Bool print_a, const Bool print_r );

void ExtendPlacement( const read_place& p, const basevector& b, const qualvector& q,
     const HyperBasevector& hb_fw, const HyperBasevector& hb_rc,
     const vec<int>& to_right_fw, const vec<int>& to_right_rc,
     vec<read_place>& places, int& qual_sum, const int min_qual );

#endif
