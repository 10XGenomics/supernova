///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "CoreTools.h"
#include "PackAlign.h"
#include "PairsManager.h"
#include "PrintAlignment.h"
#include "Qualvector.h"
#include "VecUtilities.h"
#include "math/Functions.h"
#include "paths/long/MakeKmerStuff.h"
#include "paths/long/ReadStack.h"
#include "random/Bernoulli.h"
#include <cstring>
#include <algorithm>

template <class MetricType>
class BaseMetrics
{
    typedef std::pair<MetricType,int> ValT;
public:
    BaseMetrics()
    { for ( unsigned idx = 0u; idx != 4u; ++idx ) mVals[idx].second = idx; }

    MetricType& val( unsigned idx )
    { AssertLt(idx,4u); return mVals[idx].first; }

    int id( unsigned idx )
    { AssertLt(idx,4u); return mVals[idx].second; }

    void sort() { std::sort(mVals,mVals+4); }
    void reverseSort() { std::sort(mVals,mVals+4,std::greater<ValT>()); }

private:
    ValT mVals[4];
};

template <unsigned N>
class PrecomputedBinomialSums
{
public:
    PrecomputedBinomialSums( unsigned w, double p )
    { for ( unsigned n = w; n < N; ++n )
      { double* pVal = mVal[n];
        for ( unsigned k = 0; k != n; ++k )
          *pVal++ = log10(BinomialSum(n,k,p)); } }

    double const* operator[]( unsigned n ) const
    { AssertLt(n,N); return mVal[n]; }

private:
    double mVal[N][N];
};

void readstack::Initialize( const int nrows, const int ncols )
{    cols_ = ncols;
     bases_.clear().resize(nrows);
     for ( auto& read : bases_ ) read.assign(ncols,' ');
     quals_.clear().resize(nrows);
     for ( auto& qual : quals_ ) qual.assign(ncols,-1);
     id_.resize_and_set( nrows, -1 );
     rc2_.resize_and_set( nrows, False ); 
     pid_.resize_and_set( nrows, -1 ); 
     pair_pos_.resize_and_set( nrows, -1 );
     offset_.resize_and_set( nrows, -1 ); 
     len_.resize_and_set( nrows, -1 );    }

void readstack::Initialize( const int64_t id1, Friends const& aligns,
     const int64_t start, const int64_t stop, con_type ctype,
     const vecbasevector& bases, const vecqualvector& quals,
     const PairsManager& pairs, const Bool use_pairs )
{    
     int n = stop - start + 1;
     int k = bases[id1].size( );
     if ( ctype == right_extended )
     {    for ( int j = 1; j < n; j++ )
          {    int id2 = aligns[start+j-1].readId();
               int offset = aligns[start+j-1].offset();
               k = Max( k, offset + bases[id2].isize( ) );    }    }
     Initialize( n, k );
     for ( int j = 0; j < n; j++ )
     {    int id2, offset;
          Bool rc2;
          if ( j == 0 )
          {    id2 = id1;
               offset = 0;
               rc2 = False;    }
          else
          {    id2 = aligns[start+j-1].readId();
               offset = aligns[start+j-1].offset();
               rc2 = aligns[start+j-1].isRC();    }
          const basevector& b2 = bases[id2];
          const qualvector& q2 = quals[id2];
          for ( int p2 = 0; p2 < b2.isize( ); p2++ )
          {    int p1 = p2 + offset;
               if ( p1 < 0 ) continue;
               if( ctype == strict && p1 >= k ) continue;
               if ( !rc2 )
               {    SetBase( j, p1, b2[p2] );
                    SetQual( j, p1, q2[p2] );    }
               else
               {    SetBase( j, p1, 3 - b2[ b2.isize( ) - p2 - 1 ] );
                    SetQual( j, p1, q2[ b2.isize( ) - p2 - 1 ] );    }    }
          SetOffset( j, offset );
          SetLen( j, b2.size( ) );
          SetId( j, id2 );
          SetRc2( j, rc2 );
          if (use_pairs)
          {    int64_t pid = pairs.getPairID(id2);
               SetPid( j, pid );
               SetPairPos( j, ( pairs.ID1(pid) == id2 ? 0 : 1 ) );    }    }    }

void readstack::Initialize( const int64_t id1, Friends const& aligns,
     const int64_t start, const int64_t stop, con_type ctype,
     const vecbasevector& bases, const VecPQVec& quals,
     const PairsManager& pairs, const Bool use_pairs )
{
     int n = stop - start + 1;
     int k = bases[id1].size( );
     if ( ctype == right_extended )
     {    for ( int j = 1; j < n; j++ )
          {    int id2 = aligns[start+j-1].readId();
               int offset = aligns[start+j-1].offset();
               k = Max( k, offset + bases[id2].isize( ) );    }    }
     Initialize( n, k );
     qvec q2;
     for ( int j = 0; j < n; j++ )
     {    int id2, offset;
          Bool rc2;
          if ( j == 0 )
          {    id2 = id1;
               offset = 0;
               rc2 = False;    }
          else
          {    id2 = aligns[start+j-1].readId();
               offset = aligns[start+j-1].offset();
               rc2 = aligns[start+j-1].isRC();    }
          const basevector& b2 = bases[id2];
          quals[id2].unpack(&q2);
          for ( int p2 = 0; p2 < b2.isize( ); p2++ )
          {    int p1 = p2 + offset;
               if ( p1 < 0 ) continue;
               if( ctype == strict && p1 >= k ) continue;
               if ( !rc2 )
               {    SetBase( j, p1, b2[p2] );
                    SetQual( j, p1, q2[p2] );    }
               else
               {    SetBase( j, p1, 3 - b2[ b2.isize( ) - p2 - 1 ] );
                    SetQual( j, p1, q2[ b2.isize( ) - p2 - 1 ] );    }    }
          SetOffset( j, offset );
          SetLen( j, b2.size( ) );
          SetId( j, id2 );
          SetRc2( j, rc2 );
          if (use_pairs)
          {    int64_t pid = pairs.getPairID(id2);
               SetPid( j, pid );
               SetPairPos( j, ( pairs.ID1(pid) == id2 ? 0 : 1 ) );    }    }    }

void readstack::AddToStack( const vec< triple<int,int64_t,Bool> >& offset_id_rc2,
     const vecbasevector& bases, const vecqualvector& quals,
     const PairsManager& pairs )
{    int oldlen = bases_.size( );
     int newlen = bases_.size( ) + offset_id_rc2.size( );
     bases_.resize(newlen);
     quals_.resize(newlen);
     id_.resize(newlen);
     rc2_.resize(newlen);
     pid_.resize(newlen);
     pair_pos_.resize(newlen);
     offset_.resize(newlen);
     len_.resize(newlen);
     for ( int m = 0; m < offset_id_rc2.isize( ); m++ )
     {    int offset = offset_id_rc2[m].first;
          int id = offset_id_rc2[m].second;
          Bool rc2 = offset_id_rc2[m].third;
          int nm = oldlen + m;
          const basevector& b = bases[id];
          const qualvector& q = quals[id];
          bases_[nm].resize( Cols( ) );
          quals_[nm].resize( Cols( ) );
          for ( int j = 0; j < Cols( ); j++ )
          {    bases_[nm][j] = ' ';
               quals_[nm][j] = -1;    }
          for ( int p2 = 0; p2 < b.isize( ); p2++ )
          {    int p1 = p2 + offset;
               if ( p1 < 0 ) continue;
               if( p1 < 0 || p1 >= Cols( ) ) continue;
               if ( !rc2 )
               {    SetBase( nm, p1, b[p2] );
                    SetQual( nm, p1, q[p2] );    }
               else
               {    SetBase( nm, p1, 3 - b[ b.isize( ) - p2 - 1 ] );
                    SetQual( nm, p1, q[ b.isize( ) - p2 - 1 ] );    }    }
          SetOffset( nm, offset );
          SetLen( nm, b.size( ) );
          SetId( nm, id );
          int64_t pid = pairs.getPairID(id);
          SetPid( nm, pid );
          SetRc2( nm, rc2 );
          SetPairPos( nm, ( pairs.ID1(pid) == id ? 0 : 1 ) );    }    }

template< class VB, class VQ > void readstack::AddToStack2( 
     const vec< triple<int,int64_t,Bool> >& offset_id_rc2, VB& bases, VQ& quals )
{    int oldlen = bases_.size( );
     int newlen = bases_.size( ) + offset_id_rc2.size( );
     bases_.resize(newlen);
     quals_.resize(newlen);
     id_.resize(newlen);
     rc2_.resize(newlen);
     pid_.resize(newlen);
     pair_pos_.resize(newlen);
     offset_.resize(newlen);
     len_.resize(newlen);
     for ( int m = 0; m < offset_id_rc2.isize( ); m++ )
     {    int offset = offset_id_rc2[m].first;
          int id = offset_id_rc2[m].second;
          Bool rc2 = offset_id_rc2[m].third;
          int nm = oldlen + m;
          const basevector& b = bases[id];
          qualvector q;
          quals[id].unpack(&q);
          bases_[nm].resize( Cols( ) );
          quals_[nm].resize( Cols( ) );
          for ( int j = 0; j < Cols( ); j++ )
          {    bases_[nm][j] = ' ';
               quals_[nm][j] = -1;    }
          for ( int p2 = 0; p2 < b.isize( ); p2++ )
          {    int p1 = p2 + offset;
               if ( p1 < 0 ) continue;
               if( p1 < 0 || p1 >= Cols( ) ) continue;
               if ( !rc2 )
               {    SetBase( nm, p1, b[p2] );
                    SetQual( nm, p1, q[p2] );    }
               else
               {    SetBase( nm, p1, 3 - b[ b.isize( ) - p2 - 1 ] );
                    SetQual( nm, p1, q[ b.isize( ) - p2 - 1 ] );    }    }
          SetOffset( nm, offset );
          SetLen( nm, b.size( ) );
          SetId( nm, id );
          SetRc2( nm, rc2 );    
          SetPid( nm, id/2 );
          SetPairPos( nm, id % 2 );    }     }

template void readstack::AddToStack2( 
     const vec< triple<int,int64_t,Bool> >& offset_id_rc2, 
     MasterVec<basevector>& bases, MasterVec<PQVec>& quals );

template void readstack::AddToStack2( 
     const vec< triple<int,int64_t,Bool> >& offset_id_rc2, 
     VirtualMasterVec<basevector>& bases, VirtualMasterVec<PQVec>& quals );

void readstack::Erase( const vec<Bool>& to_remove )
{    bases_.EraseIf( to_remove );
     quals_.EraseIf( to_remove );
     EraseIf( id_, to_remove );
     EraseIf( rc2_, to_remove );
     EraseIf( pid_, to_remove );
     EraseIf( pair_pos_, to_remove );
     EraseIf( offset_, to_remove );
     EraseIf( len_, to_remove );    }

void readstack::Unique( ) // only works with SortByPid
{    vec<Bool> to_remove( Rows( ), False );
     for ( int i = 0; i < Rows( ); i++ )
     {    int j;
          for ( j = i + 1; j < Rows( ); j++ )
          {    if ( id_[j] != id_[i] ) break;
               if ( rc2_[j] != rc2_[i] ) break;
               if ( pid_[j] != pid_[i] ) break;
               if ( pair_pos_[j] != pair_pos_[i] ) break;
               if ( offset_[j] != offset_[i] ) break;
               if ( len_[j] != len_[i] ) break;    }
          for ( int k = i + 1; k < j; k++ )
          {    for ( int c = 0; c < Cols( ); c++ )
               {    if ( !Def(i,c) && Def(k,c) )
                    {    bases_[i][c] = bases_[k][c];
                         quals_[i][c] = quals_[k][c];    }    }
               to_remove[k] = True;    }
          i = j - 1;    }
     for ( int i = 0; i < 2; i++ )
     for ( int j = 2; j < Rows( ); j++ )
     {    if ( id_[j] != id_[i] ) continue;
          if ( rc2_[j] != rc2_[i] ) continue;
          if ( pid_[j] != pid_[i] ) continue;
          if ( pair_pos_[j] != pair_pos_[i] ) continue;
          if ( offset_[j] != offset_[i] ) continue;
          if ( len_[j] != len_[i] ) continue;
          to_remove[j] = True;
          for ( int c = 0; c < Cols( ); c++ )
          {    if ( !Def(i,c) && Def(j,c) )
               {    bases_[i][c] = bases_[j][c];
                    quals_[i][c] = quals_[j][c];    }    }    }
     Erase(to_remove);    }

void readstack::SortByPid( const int64_t pid1, const int i1, const int i2 )
{    ForceAssertGe( Rows( ), 2 );
     vec<int> ident( Rows( ), vec<int>::IDENTITY ), identr( Rows( ) );
     vec<int64_t> pid( Rows( ) ), offsetp( Rows( ) );
     for ( int i = 0; i < Rows( ); i++ )
          pid[i] = Pid(i);
     SortSync( pid, ident );
     for ( int i = 0; i < Rows( ); i++ )
     {    int j = pid.NextDiff(i);
          int off = Offset( ident[i] );
          for ( int k = i+1; k < j; k++ )
               off = Min( off, Offset( ident[k] ) );
          if ( pid[i] == pid1 ) off = -1000000000;
          for ( int k = i; k < j; k++ )
               offsetp[k] = off;
          i = j - 1;    }
     vec< triple< int64_t, int64_t, pair<Bool,int> > > data( Rows( ) );
     for ( int i = 0; i < Rows( ); i++ )
          data[i] = make_triple( offsetp[i], Pid( ident[i] ), 
               make_pair( Rc2( ident[i] ), Offset( ident[i] ) ) );
     SortSync( data, ident );
     if ( ident[0] != i1 )
     {    for ( int j = 0; j < Rows( ); j++ )
          {    if ( ident[j] == i1 ) 
               {    swap( ident[0], ident[j] );
                    break;    }    }    }
     if ( ident[1] != i2 )
     {    for ( int j = 0; j < Rows( ); j++ )
          {    if ( ident[j] == i2 ) 
               {    swap( ident[1], ident[j] );
                    break;    }    }    }
     for ( int i = 0; i < Rows( ); i++ )
          identr[ ident[i] ] = i;
     PermuteVec( bases_, identr );
     PermuteVec( quals_, identr );
     PermuteVec( id_, identr );
     PermuteVec( rc2_, identr );
     PermuteVec( pid_, identr );
     PermuteVec( pair_pos_, identr );
     PermuteVec( offset_, identr );
     PermuteVec( len_, identr );    }

void readstack::Reverse( )
{    for ( int i = 0; i < Rows( ); i++ )
     {    bases_[i].ReverseMe( );
          for ( int j = 0; j < Cols( ); j++ )
               if ( bases_[i][j] != ' ' ) bases_[i][j] = 3 - bases_[i][j];
          quals_[i].ReverseMe( );
          rc2_[i] = !rc2_[i];
          offset_[i] = - ( offset_[i] + len_[i] - Cols( ) );    }    }

void readstack::Merge( const readstack& s, const int offset )
{    int rows1 = Rows( ), rows2 = s.Rows( ), cols1 = Cols( ), cols2 = s.Cols( );
     int left_ext1 = Max( 0, -offset );
     int right_ext1 = Max( 0, offset + cols2 - cols1 );
     for ( int i = 0; i < rows1; i++ )
     {    if ( left_ext1 > 0 )
          {    bases_[i].resize( left_ext1 + cols1 );
               quals_[i].resize( left_ext1 + cols1 );
               for ( int j = left_ext1 + cols1 - 1; j >= left_ext1; j-- )
               {    bases_[i][j] = bases_[i][j-left_ext1];
                    quals_[i][j] = quals_[i][j-left_ext1];    }
               for ( int j = 0; j < left_ext1; j++ )
               {    bases_[i][j] = ' ';
                    quals_[i][j] = -1;    }    }
          bases_[i].resize( bases_[i].size( ) + right_ext1, ' ' );
          quals_[i].resize( quals_[i].size( ) + right_ext1, -1 );    }
     bases_.append( s.Bases( ).begin(), s.Bases( ).end()  );
     quals_.append( s.Quals( ).begin(), s.Quals( ).end() );
     int left_ext2 = Max( 0, offset );
     int right_ext2 = Max( 0, cols1 - ( offset + cols2 ) );
     for ( int i = 0; i < rows2; i++ )
     {    int ip = rows1 + i;
          if ( left_ext2 > 0 )
          {    bases_[ip].resize( left_ext2 + cols2 );
               quals_[ip].resize( left_ext2 + cols2 );
               for ( int j = left_ext2 + cols2 - 1; j >= left_ext2; j-- )
               {    bases_[ip][j] = bases_[ip][j-left_ext2];
                    quals_[ip][j] = quals_[ip][j-left_ext2];    }
               for ( int j = 0; j < left_ext2; j++ )
               {    bases_[ip][j] = ' ';
                    quals_[ip][j] = -1;    }    }
          bases_[ip].resize( bases_[ip].size( ) + right_ext2, ' ' );
          quals_[ip].resize( quals_[ip].size( ) + right_ext2, -1 );    }
     offset_.append( s.Offset( ) );
     for ( int i = 0; i < rows1; i++ )
          SetOffset( i, Offset(i) + left_ext1 );
     for ( int i = rows1; i < rows1 + rows2; i++ )
          SetOffset( i, Offset(i) + ( offset > 0 ? offset : 0 ) );
     id_.append( s.Id( ) );
     rc2_.append( s.Rc2( ) );
     pid_.append( s.Pid( ) );
     pair_pos_.append( s.PairPos( ) );
     len_.append( s.Len( ) );
     cols_ = bases_[0].size( );    }

basevector readstack::Consensus1( ) const
{    basevector con( Cols( ) );
     for ( int i = 0; i < Cols( ); i++ )
          con.Set( i, ColumnConsensus1(i) );
     return con;    }

void readstack::Consensus1( basevector& con, qualvector& conq, const int qualcap ) 
     const
{    con.resize( Cols( ) ), conq.resize( Cols( ) );
     for ( int i = 0; i < Cols( ); i++ )
     {    
          // Compute quality score sum for each base.  Count Q0 as 0.1, Q1 as 0.2,
          // and Q2 as 0.2.
          BaseMetrics<double> mx;
          for ( int j = 0; j < Rows( ); j++ )
          {    double q = Qual(j,i);
               if ( q <= 2 ) q = Min( q, 0.2 );
               if ( q == 0 ) q = 0.1;
               if ( Qual(j,i) >= 0 ) mx.val(Base(j,i)) += q;    }
          mx.reverseSort();
          con.Set( i, mx.id(0) );
          conq[i] = Min( qualcap, int(round(mx.val(0)-mx.val(1))) );
          const int max_qcomp = 100;
          if ( mx.val(1) > max_qcomp )
          {    int badcount = 0;
               for ( int j = 0; j < Rows( ); j++ )
                    if ( Qual(j,i) >= 30 && Base(j,i) == mx.id(1) ) badcount++;
               if ( badcount >= 2 ) conq[i] = 0;    }    }    }

void readstack::StrongConsensus1( basevector& con, qualvector& conq,
     const Bool raise_zero ) const
{    con.resize( Cols( ) ), conq.resize( Cols( ) );
     for ( int i = 0; i < Cols( ); i++ )
          con.Set( i, ColumnConsensus1(i) );

     const int min_window = 41;
     const double qfudge = 0.5;
     // TODO: this is a mistake -- should be a double metric value, but I don't
     // want to change results.  --Ted 5/16/14
     vec<BaseMetrics<int>> sum(Cols());
     vec<double> q;
     for ( int j = 0; j < Rows( ); j++ )
     {    auto const& qs = quals_[j];
          q.assign(qs.begin(),qs.end());
          auto const& bs = bases_[j];
          auto beg=bs.begin();
          auto end=bs.end();
          auto cItr = con.begin();
          for ( auto itr=beg; itr != end; ++itr,++cItr )
          {
              auto itrPair = std::mismatch(itr,end,cItr);
              if ( itrPair.first-itr >= min_window )
              {
                  int i1 = itr-beg;
                  int i2 = itrPair.first-beg;
                  int i3 = i2-min_window/2;
                  for ( int l = i1+min_window/2; l <= i3; ++l )
                  {
                      int dist = std::min(l-i1,i2-l-1);
                      if ( 2*dist < min_window )
                          continue;
                      if ( !raise_zero && q[l] == 0 )
                          continue;
                      q[l] = std::max(q[l], 10.0*log10(2*dist)*qfudge);
                  }
              }
              if ( itrPair.first == end ) break;
              itr = itrPair.first;
              cItr = itrPair.second;
          }
          for ( int i = 0; i < Cols( ); i++ )
          {    if ( qs[i] >= 0 )
               {    double p = q[i];
                    if ( p <= 2 ) p = std::min( p, 0.2 );
                    if ( p == 0 ) p = 0.1;
                    sum[i].val(bs[i]) += p;    }    }    }

     for ( int i = 0; i < Cols( ); i++ )
     {    BaseMetrics<int>& mx = sum[i];
          mx.reverseSort();
          const int qual_cap = 50;
          conq[i] = Min( qual_cap, int(round(mx.val(0)-mx.val(1))) );
          const int max_qcomp = 100;
          if ( mx.val(1) > max_qcomp )
          {    int badcount = 0;
               for ( int j = 0; j < Rows( ); j++ )
                    if ( Qual(j,i) >= 30 && Base(j,i) == mx.id(1) ) badcount++;
               if ( badcount >= 2 ) conq[i] = 0;    }    }    }

void readstack::HighQualDiff( const int n, const int top, vec<Bool>& suspect ) const
{    suspect.assign( Rows( ), False );
     for ( int t = 0; t < top; t++ )
     for ( int j = top; j < Rows( ); j++ )
     {    for ( int c = 0; c < Cols( ); c++ )
          {    if ( Base(j,c) != Base(t,c) && Qual(j,c) >= n && Qual(t,c) >= n )
               {    suspect[j] = True;
                    break;    }    }    }    }

void readstack::CleanColumns( const int top, vec<Bool>& suspect ) const
{    suspect.assign( Rows( ), False );
     for ( int c = 0; c < Cols( ); c++ )
     {    const int min_q = 20;
          const int min_count = 3;
          vec<int> count(4, 0);
          for ( int j = 0; j < Rows( ); j++ )
               if ( Qual(j,c) >= min_q ) count[ Base(j,c) ]++;
          int called = 0;
          for ( int j = 0; j < 4; j++ )
               if ( count[j] >= min_count ) called++;
          if ( called >= 2 )
          {    for ( int t = 0; t < top; t++ )
               for ( int j = top; j < Rows( ); j++ )
               {    if ( Base(j,c) != Base(t,c) && Qual(j,c) >= min_q 
                         && Qual(t,c) >= min_q && count[ Base(t,c) ] >= min_count )
                    {    suspect[j] = True;    }    }    }    }    }

void readstack::Print( ostream& out, const vec<vec<char>>& con, const int w,
     const vec<String>& title ) const
{
     // Determine which reads have their partner placed.

     vec< pair<int64_t,int> > have;
     for ( int i = 0; i < Rows( ); i++ )
          have.push( Pid(i), PairPos(i) );
     UniqueSort(have);
     vec<int64_t> have2;
     for ( int i = 0; i < have.isize( ) - 1; i++ )
     {    if ( have[i].first == have[i+1].first )
               have2.push_back( have[i].first );    }

     // Set up for quality score display.

     int n = Rows( ), k = Cols( );
     vec< vec<char> > callq1( n, vec<char>( k, ' ' ) );
     vec< vec<char> > callq2( n, vec<char>( k, ' ' ) );
     for ( int j = 0; j < n; j++ )
     {    for ( int c = 0; c < k; c++ )
          {    int qa = Qual(j,c) / 10, qb = Qual(j,c) % 10;
               if ( Qual(j,c) >= 0 )
               {    callq1[j][c] = ( qa == 0 ? ' ' : '0' + qa );
                    callq2[j][c] = '0' + qb;    }    }    }

     // Display windows.

     for ( int cstart = 0; cstart < k; cstart += w )
     {    int c1 = cstart, c2 = Min( cstart + w, k );
          int wid = cstart/w + 1;
          out << "\n=================================================="
               << "==================================" << endl;
          out << "\n\n" << "[" << wid << "] window from " << c1 
               << " to " << c2 << "\n\n\n";

          int cdepth = 0;
          for ( int c = c1; c < c2; c++ )
               cdepth = Max( cdepth, con[c].isize( ) );
          for ( int j = 0; j < cdepth; j++ )
          {    for ( int c = c1; c < c2; c++ )
               {    if ( j < con[c].isize( ) ) out << as_base( con[c][j] );
                    else out << " ";    }
               out << "   consensus\n";    }
          out << "\n\n";
               
          for ( int j = 0; j < n; j++ )
          {    
               vec<String> tag(5);
               tag[0] = "   id=" + ToString( Id(j) );
               tag[1] = String("   ") + ( Rc2(j) ? "rc" : "fw" );
               tag[2] = "   " + ToString( Offset(j) );
               tag[3] = String("   ") 
                    + ( BinMember( have2, Pid(j) ) ? "pp!" 
                    : ( !Rc2(j) ? "un" : "" ) );
               tag[4] = "   window [" + ToString(wid) + "]";
               int tt = 0;

               Bool present = False, nonblank = False;
               for ( int c = c1; c < c2; c++ )
                    if ( Base(j,c) != ' ' ) present = True;
               if ( !present ) continue;

               for ( int c = c1; c < c2; c++ )
                    out << ( Base(j,c) == ' ' ? ' ' : '-' );
               if ( title.size( ) > 0 ) out << "   " << title[j];
               out << "\n";

               for ( int c = c1; c < c2; c++ )
                    if ( callq1[j][c] != ' ' ) nonblank = True;
               if (nonblank)
               {    for ( int c = c1; c < c2; c++ )
                         out << callq1[j][c];
                    if ( tt < 5 ) out << tag[tt++];
                    out << "\n";    }
               for ( int c = c1; c < c2; c++ )
                    out << callq2[j][c];
               if ( tt < 5 ) out << tag[tt++];
               out << "\n";
               for ( int c = c1; c < c2; c++ )
                    out << ( Base(j,c) == ' ' ? ' ' : ':' );
               if ( tt < 5 ) out << tag[tt++];
               out << "\n";

               int ms = 0;
               for ( int pass = 1; pass <= 2; pass++ )
               {    if ( pass == 2 && ms == 0 ) break;
                    for ( int c = c1; c < c2; c++ )
                    {    char m;  
                         if ( Base(j,c) == ' ' ) m = ' ';
                         else if ( Member( con[c], Base(j,c) ) ) m = ' ';
                         else if ( con[c].empty( ) ) m = '-';
                         else m = '*';
                         if ( m == '*' || m == '-' ) ms++;    
                         if ( pass == 2 ) out << m;    }
                    if ( pass == 2 ) 
                    {    if ( tt < 5 ) out << tag[tt++];
                         out << "\n";    }    }
               for ( int c = c1; c < c2; c++ )
               {    if ( Base(j,c) == ' ' ) out << ' '; // ?????????????
                    else out << as_base( Base(j,c) );    }
               if ( tt < 5 ) out << tag[tt++];
               out << "\n";
               for ( int c = c1; c < c2; c++ )
                    out << ( Base(j,c) == ' ' ? ' ' : '-' );
               if ( tt < 5 ) out << tag[tt++];
               out << "\n\n\n";    }    }    }

void readstack::Print( ostream& out, const vec<basevector>& cons, const int w,
     const vec<String>& title ) const
{    vec< vec<char> > con( Cols( ) );
     for ( int i = 0; i < Cols( ); i++ )
     {    vec<Bool> have( 4, False );
          for ( int j = 0; j < cons.isize( ); j++ )
               have[ cons[j][i] ] = True;
          for ( int j = 0; j < 4; j++ )
               if ( have[j] ) con[i].push_back(j);    }
     Print( out, con, w, title );    }

void readstack::Print( ostream& out, const int w,
     const vec<String>& title ) const
{    int n = Rows( ), k = Cols( );
     basevector con;
     qualvector conq;
     Consensus1( con, conq );
     vec<basevector> cons;
     cons.push_back(con);
     Print( out, cons, w, title );    }

void readstack::Raise1( const int id, const int rwindow, 
     const Bool require_unedited )
{
     // The following documentation is not entirely accurate.
     //
     // Raise quality scores.  For each 11 base window, suppose first that the
     // founder has not been edited on the window, except possibly for raising 
     // quality scores, and has score < 30 at the middle base.  Suppose in addition 
     // that there are at least 3 reads that agree with the founder, have not been 
     // edited on the window, except possibly for editing quality scores, have score 
     // > 30 at the middle base, and have no Q30 difference with the founder.
     // Then change the score of the middle base on the given read to 30.  
     //
     // Amelioration.  Don't raise q score if there is a viable alternate 
     // hypothesis.  This is a sequence in the window that agrees with the founder,
     // except for the middle base, and which occurs at least 3 times in reads 
     // that have no Q30 difference with the founder.
     //
     // Further amelioration.  These alts must be Q30 at the middle base.

     const int min_agree = 3;
     const int critical_q = 30;
     vec<int> raises;
     for ( int c = 0; c <= Cols( ) - rwindow; c++ )
     {    if ( Qual( id, c + rwindow/2 ) >= critical_q ) continue;

          if (require_unedited)
          {    Bool edited = False;
               for ( int l = 0; l < rwindow; l++ )
                    if ( Qual( id, c + l ) == 0 ) edited = True;
               if (edited) continue;    }

          if ( Qual( id, c + rwindow/2 ) == 0 ) continue;

          Bool def = True;
          for ( int l = 0; l < rwindow; l++ )
               if ( !Def( id, c+l ) ) def = False;
          if ( !def ) continue;

          int support = 0;
          for ( int j = 0; j < Rows( ); j++ )
          {    if ( j == id ) continue;
               if ( Qual( j, c + rwindow/2 ) < critical_q ) continue;
               Bool bad = False;
               for ( int l = 0; l < rwindow; l++ )
               {    if ( Qual( j, c+l ) <= 0 || Base( j, c+l ) != Base( id, c+l ) )
                    {    bad = True;
                         break;    }    }
               if (bad) continue;
               if ( ++support == min_agree ) break;    }
          if ( support < min_agree ) continue;

          int alts[4]{0,0,0,0};
          for ( int j = 0; j < Rows( ); j++ )
          {    if ( j == id ) continue;
               Bool bad = False;
               for ( int l = 0; l < rwindow; l++ )
               {    if ( l == rwindow/2 ) continue;
                    if ( Qual( j, c+l ) <= 0 || Base( j, c+l ) != Base( id, c+l ) )
                    {    bad = True;
                         break;    }    }
               if (bad) continue;
               if ( Qual( j, c + rwindow/2 ) < critical_q ) continue;
               char x = Base( j, c + rwindow/2 );
               if ( x != Base( id, c + rwindow/2 ) ) alts[int(x)]++;    }
          if ( *std::max_element(alts,alts+4) >= min_agree ) continue;

          SetQual( id, c + rwindow/2, critical_q );    }    }

void readstack::Trim( const int start, const int stop )
{    vec<Bool> to_remove( Rows( ), True );
     for ( int i = 0; i < Rows( ); i++ )
     {    for ( int j = start; j < stop; j++ )
          {    if ( Def(i,j) ) 
               {    to_remove[i] = False;
                    break;    }    }
          bases_[i].SetToSubOf( bases_[i], start, stop - start );
          quals_[i].SetToSubOf( quals_[i], start, stop - start );
          offset_[i] -= start;    }
     cols_ = stop - start;
     Erase(to_remove);    }

void readstack::PairWeak1( vec<Bool>& suspect ) const
{    suspect.assign( Rows( ), False );
     vec< pair<int64_t,int> > have;
     for ( int i = 0; i < Rows( ); i++ )
          have.push( Pid(i), PairPos(i) );
     UniqueSort(have);
     vec<int64_t> have2;
     for ( int i = 0; i < have.isize( ) - 1; i++ )
          if ( have[i].first == have[i+1].first ) have2.push_back( have[i].first );
     vec<Bool> paired( Rows( ) );
     for ( int j = 0; j < Rows( ); j++ )
          paired[j] = BinMember( have2, Pid(j) );
     for ( int i = 0; i < Cols( ); i++ )
     {    BaseMetrics<int> mx;
          for ( int j = 0; j < Rows( ); j++ )
          {    if ( !paired[j] || !Def(j,i) ) continue;
               mx.val(Base(j,i)) += Qual(j,i);    }
          mx.reverseSort();
          if ( mx.val(0) >= 100 && mx.val(0) > 10*mx.val(1) && mx.val(1) < 100 )
          for ( int j = 0; j < Rows( ); j++ )
          {    if ( Def(j,i) && Base(j,i) != mx.id(0) && Qual(j,i) >= 30 )
                    suspect[j] = True;    }    }    }

void PairWeak2( const readstack& s1, const readstack& s2,
     vec<Bool>& suspect1, vec<Bool>& suspect2 )
{    suspect1.resize_and_set( s1.Rows( ), False );
     suspect2.resize_and_set( s2.Rows( ), False );
     vec< pair<int64_t,int> > have;
     for ( int i = 0; i < s1.Rows( ); i++ )
          have.push( s1.Pid(i), s1.PairPos(i) );
     for ( int i = 0; i < s2.Rows( ); i++ )
          have.push( s2.Pid(i), s2.PairPos(i) );
     UniqueSort(have);
     vec<int64_t> have2;
     for ( int i = 0; i < have.isize( ) - 1; i++ )
          if ( have[i].first == have[i+1].first ) have2.push_back( have[i].first );
     for ( int pass = 1; pass <= 2; pass++ )
     {    const readstack& s = ( pass == 1 ? s1 : s2 );
          vec<Bool>& suspect = ( pass == 1 ? suspect1 : suspect2 );
          for ( int i = 0; i < s.Cols( ); i++ )
          {    BaseMetrics<int> mx;
               for ( int j = 0; j < s.Rows( ); j++ )
               {    if ( !s.Def(j,i) || !BinMember( have2, s.Pid(j) ) ) continue;
                    mx.val(s.Base(j,i)) += s.Qual(j,i);    }
               mx.reverseSort();
               if ( mx.val(0) >= 100 && mx.val(0) > 10*mx.val(1) && mx.val(1) < 100 )
               for ( int j = 0; j < s.Rows( ); j++ )
               {    if ( s.Def(j,i) && s.Base(j,i) != mx.id(0) && s.Qual(j,i) >= 30 )
                         suspect[j] = True;    }    }    }    }

namespace
{

class MotifLoc
{
public:
    static int const WIDTH = 10;

    MotifLoc( char const* bases, int id ) : mpBases(bases), mId(id) {}

    friend bool operator==( MotifLoc const& ml1, MotifLoc const& ml2 )
    { return std::equal(ml1.mpBases,ml1.mpBases+WIDTH,ml2.mpBases); }
    friend bool operator!=( MotifLoc const& ml1, MotifLoc const& ml2 )
    { return !(ml1==ml2); }
    friend bool operator<( MotifLoc const& ml1, MotifLoc const& ml2 )
    { return std::lexicographical_compare(ml1.mpBases,ml1.mpBases+WIDTH,
                                          ml2.mpBases,ml2.mpBases+WIDTH); }

    char const* mpBases;
    int mId;
};

}
void readstack::MotifDiff( const int top, vec<Bool>& to_delete ) const
{    int n = Rows( ), k = Cols( );
     to_delete.assign( n, False );
     int const WIDTH = MotifLoc::WIDTH;
     int const MIN_MULT = 10;
     // for each WIDTH-wide vertical stripe of columns
     vec<MotifLoc> X;
     vec<int> bigs;
     for ( int i = 0; i <= k - WIDTH; i += WIDTH )
     {    X.clear();
          // accumulate all the fully-defined WIDTH-mers
          for ( int j = 0; j < n; j++ )
          {    char const* kBeg = &bases_[j][i];
               char const* kEnd = kBeg+WIDTH;
               if ( std::find(kBeg,kEnd,' ') != kEnd ) continue;
               X.emplace_back(kBeg,j);    }
          std::sort(X.begin(),X.end());
          Bool found = False;
          bigs.clear();
          int this_one = -1;
          // for each horizontal stripe of identical WIDTH-mers
          for ( int r = 0; r < X.isize( ); r++ )
          {    int s = X.NextDiff(r);
               if ( s - r >= MIN_MULT )
               {    bool agree = false;
                    for ( int m = 0; !agree && m < top; m++ )
                         agree = MotifLoc(&bases_[m][i],m) == X[r];
                    if (agree) this_one = bigs.size( );
                    bigs.push_back(r);    }
               r = s - 1;    }
          if ( this_one >= 0 )
          {    char const* quals = &quals_[0][i];
               for ( int b = 0; b < bigs.isize( ); b++ )
               {    if ( b == this_one ) continue;
                    int that_one = bigs[b];
                    Bool hq_diff = False;
                    char const* those_bases = X[that_one].mpBases;
                    char const* these_bases = X[this_one].mpBases;
                    for ( int l = 0; l < WIDTH; l++ )
                    {    if ( those_bases[l] != these_bases[l] )
                              if ( quals[l] >= 20 ) hq_diff = True;    }
                    if ( !hq_diff ) continue;
                    for ( int d = that_one; d < X.isize( ); d++ )
                    {    if ( X[d] != X[that_one] ) break;
                         to_delete[ X[d].mId ] = True;    }    }    }    }   }

vec<basevector> readstack::Consensuses1( ) const
{    vec< vec<char> > cons, cons_term;
     vec<char> init;
     for ( int i = 0; i < Rows( ); i++ )
          if ( Qual(i,0) >= 30 ) init.push_back( Base(i,0) );
     UniqueSort(init);
     for ( int i = 0; i < init.isize( ); i++ )
     {    vec<char> x;
          x.push_back( init[i] );
          cons.push_back(x);    }
     for ( int j = 1; j < Cols( ); j++ )
     {    vec< vec<char> > cons2;
          for ( int i = 0; i < cons.isize( ); i++ )
          {    const vec<char>& x = cons[i];
               vec<char> nexts;
               for ( int i = 0; i < Rows( ); i++ )
               {    if ( Qual(i,j) < 30 ) continue;
                    Bool mismatch = False;
                    int support = 0;
                    for ( int l = 0; l < j; l++ )
                    {    if ( j-l <= 40 )
                         {    if ( Base(i,l) != x[l] )
                              {    mismatch = True;
                                   break;    }    }
                         if ( Qual(i,l) >= 30 )
                         {    support++;
                              if ( Base(i,l) != x[l] )
                              {    mismatch = True;
                                   break;    }    }    }
                    if (mismatch) continue;
                    // if ( support < Min( j/2, 20 ) ) continue;
                    nexts.push_back( Base(i,j) );    }
               Sort(nexts);
               vec<char> nextsx;
               for ( int m1 = 0; m1 < nexts.isize( ); m1++ )
               {    int m2 = nexts.NextDiff(m1);
                    if ( m2 - m1 >= 2 ) nextsx.push_back( nexts[m1] );
                    m1 = m2 - 1;    }
               nexts = nextsx;
               for ( int l = 0; l < nexts.isize( ); l++ )
               {    vec<char> y = x;
                    y.push_back( nexts[l] );
                    cons2.push_back(y);    }
               if ( nexts.empty( ) ) cons_term.push_back(x);    }
          cons = cons2;    }
     cons.append(cons_term);
     vec<basevector> x;
     for ( int i = 0; i < cons.isize( ); i++ )
     {    basevector b( cons[i].size( ) );
          for ( int j = 0; j < b.isize( ); j++ )
               b.Set( j, cons[i][j] );
          x.push_back(b);    }
     return x;    }

vec<basevector> readstack::Consensuses2( const int K, const int top ) const
{    const int min_mult = 2;
     const int hq = 30;
     vec<basevector> con;
     if ( K+1 > Cols( ) ) return con;
     basevector b(K+1);
     for ( int i = 0; i < Rows( ); i++ )
     {    Bool def = True;
          for ( int j = 0; j <= K; j++ )
          {    if ( !Def(i,j) ) 
               {    def = False;
                    break;    }    }
          if ( !def ) continue;
          Bool bad = False;
          for ( int j = 0; j <= K; j++ )
          {    for ( int t = 0; t < top; t++ )
               {    if ( Qual(t,j) >= hq && Base(i,j) != Base(t,j) )
                    {    bad = True;
                         break;    }    }
               if (bad) break;
               b.Set( j, Base(i,j) );    }
          if ( !bad ) con.push_back(b);    }
     Sort(con);
     vec<Bool> to_delete( con.size( ), False );
     for ( int i = 0; i < con.isize( ); i++ )
     {    int j = con.NextDiff(i);
          if ( j - i < min_mult )
          {    for ( int k = i; k < j; k++ )
                    to_delete[k] = True;    }
          else
          {    for ( int k = i + 1; k < j; k++ )
                    to_delete[k] = True;    }
          i = j - 1;    }
     EraseIf( con, to_delete );
     for ( int r = K+1; r < Cols( ); r++ )
     {    vec<basevector> con2, nexts;
          for ( int i = 0; i < Rows( ); i++ )
          {    Bool def = True;
               for ( int j = 0; j <= K; j++ )
               {    if ( !Def(i,j+r-K) ) 
                    {    def = False;
                         break;    }    }
               if ( !def ) continue;
               Bool bad = False;
               for ( int j = 0; j <= K; j++ )
               {    for ( int t = 0; t < top; t++ )
                    {    if ( Qual(t,j+r-K) >= hq && Base(i,j+r-K) != Base(t,j+r-K) )
                         {    bad = True;
                              break;    }    }
                    if (bad) break;
                    b.Set( j, Base(i,j+r-K) );    }
               if ( !bad ) nexts.push_back(b);    }
          Sort(nexts);
          vec<Bool> to_delete( nexts.size( ), False );
          for ( int i = 0; i < nexts.isize( ); i++ )
          {    int j = nexts.NextDiff(i);
               if ( j - i < min_mult )
               {    for ( int k = i; k < j; k++ )
                         to_delete[k] = True;    }
               else
               {    for ( int k = i + 1; k < j; k++ )
                         to_delete[k] = True;    }
               i = j - 1;    }
          EraseIf( nexts, to_delete );
          for ( int r = 0; r < con.isize( ); r++ )
          for ( int s = 0; s < nexts.isize( ); s++ )
          {    Bool agree = True;
               for ( int j = 0; j < K; j++ )
               {    if ( nexts[s][j] != con[r][ con[r].isize( ) - K + j ] )
                    {    agree = False;
                         break;    }    }
               if ( !agree ) continue;
               basevector c( con[r] );
               c.resize( c.size( ) + 1 );
               c.Set( con[r].size( ), nexts[s][ nexts[s].isize( ) - 1 ] );
               con2.push_back(c);    }
          con = con2;    }

     // For each consensus, find the reads that do not have a Q30 diff with it.

     int survivors = 0;
     vec<Bool> cdel( con.size( ), False );
     for ( int ci = 0; ci < con.isize( ); ci++ )
     {    const basevector& c = con[ci];
          vec<int> goods;
          vec<ho_interval> cov;
          for ( int i = 0; i < Rows( ); i++ )
          {    Bool bad = False;
               int qsum = 0;
               for ( int j = 0; j < Cols( ); j++ )
               {    if ( Base(i,j) != c[j] )
                    {    if ( Qual(i,j) >= 30 )
                         {    bad = True;
                              break;    }
                         if ( Qual(i,j) > 2 ) qsum += Qual(i,j);    }    }
               if ( qsum > 50 ) bad = True;
               if (bad) continue;
               goods.push_back(i);
               for ( int j = 0; j < Cols( ); j++ )
               {    if ( Base(i,j) != c[j] ) continue;
                    int k;
                    for ( k = j + 1; k < Cols( ); k++ )
                         if ( Base(i,k) != c[k] ) break;
                    if ( k - j > K ) 
                    {    int start = j, stop = k;
                         if ( stop < Cols( ) ) stop -= K;
                         cov.push( start, stop );    }
                    j = k - 1;    }    }

          if ( !BinMember( goods, 0 ) || !BinMember( goods, 1 ) )
          {    cdel[ci] = True;
               continue;    }

          ExtractGivenCoverage( Cols( ), min_mult, cov, cov );
          if ( TotalCovered(cov) < Cols( ) )
          {    cdel[ci] = True;
               continue;    }

          vec<double> support( Cols( ), 0 );
          for ( int gi = 0; gi < goods.isize( ); gi++ )
          {    int i = goods[gi];
               for ( int j = 0; j < Cols( ); j++ )
               {    if ( Base(i,j) != c[j] ) continue;
                    double q = Qual(i,j);
                    if ( q <= 2.0 ) q = 0.2;
                    support[j] += q;    }    }
          const double min_support = 10.0;
          Bool bad = False;
          for ( int j = 0; j < Cols( ); j++ )
          {    if ( support[j] < min_support )
               {    bad = True;
                    break;    }    }
          if (bad)
          {    cdel[ci] = True;
               continue;    }

          /*
          cout << "\nFound a consensus; here is the stack:\n";
          readstack stack(*this);
          vec<Bool> to_delete( Cols( ), True );
          for ( int i = 0; i < goods.isize( ); i++ )
               to_delete[ goods[i] ] = False;
          stack.Erase(to_delete);
          vec<basevector> cs;
          cs.push_back(c);
          stack.Print( cout, cs );
          */

          survivors++;    }

     PRINT2( con.size( ), survivors );
     EraseIf( con, cdel );

     return con;    }

void readstack::CorrectAll( basevector& b, qualvector& q, int& trim_to,
     const Bool verbose ) const
{    
     // Preset.

     b.resize( Cols( ) );
     q.resize( Cols( ) );
     for ( int i = 0; i < Cols( ); i++ )
     {    b.Set( i, Base(0,i) );
          q[i] = Qual(0,i);    }
     trim_to = -1;

     // Heuristics.

     const int min_win = 50;
     const int min_win_ratio = 10;
     const int max_lose = 100;

     // Go through the columns.

     for ( int i = 0; i < Cols( ); i++ )
     {    
          // Compute quality score sum for each base.  We count Q2 bases as
          // next to nothing.
          BaseMetrics<double> mx;
          int top[4]{0,0,0,0};
          for ( int j = 0; j < Rows( ); j++ )
          {    if ( Qual(j,i) >= 0 )
               {    double q = Qual(j,i);
                    if ( q <= 2 ) q = Min( q, 0.2 );
                    mx.val(Base(j,i)) += q;
                    int b = Base(j,i);
                    top[b] = Max( top[b], Qual(j,i) );    }    }
          mx.reverseSort();
          int winner = mx.id(0);

          // Drop the top score for competitors.

          vec< pair<int,int> > competitors;
          for ( int j = 1; j < 4; j++ )
              mx.val(j) -= top[ mx.id(j) ];

          Bool OK = False;
          if ( mx.val(0) >= min_win && mx.val(0) >= min_win_ratio * mx.val(1)
               && mx.val(1) <= max_lose )
          {    OK = True;    }
          if (verbose) 
          {    
               #pragma omp critical
               {    cout << "column " << i << ", read base = " << as_base(Base(0,i))
                         << ", sum[0] = " << mx.val(0) << ", sum[1] = " << mx.val(1)
                         << ", OK = " << ( OK ? "True" : "False" ) << endl;    }    }
          if (OK) 
          {    if ( Base(0,i) != winner )
               {    b.Set( i, winner );
                    q[i] = 0;    }    
               /*
               else if ( sum[0] >= 1000 && sum[1] <= 10 ) 
               {    q[i] = Max( (int) q[i], 30 );    }     
               */
                    }
          else if ( trim_to < 0 ) trim_to = i;    }
     if ( trim_to < 0 ) trim_to = Cols( );    }

void readstack::Recruit( const int K, const vec<simple_align_data>& aligns,
     const vec<int64_t>& id1_start, const vecbasevector& bases, 
     const vecqualvector& quals, const PairsManager& pairs )
{    vec< triple<int,int64_t,Bool> > offset_id_rc2_orig;
     vec<int64_t> expect;
     for ( int i = 0; i < Rows( ); i++ )
     {    int64_t id2 = Id(i);
          int offset2 = Offset(i);
          Bool rc2 = Rc2(i);
          offset_id_rc2_orig.push( offset2, id2, rc2 );
          expect.push( Pid(i) );    }
     Sort(offset_id_rc2_orig);
     UniqueSort(expect);
     vec< triple<int,int64_t,Bool> > offset_id_rc2;
     for ( int i = 0; i < Rows( ); i++ )
     {    int64_t id2 = Id(i);
          int offset2 = Offset(i);
          Bool rc2 = Rc2(i);
          int64_t start = id1_start[id2], stop = id1_start[id2+1];
          for ( int64_t l = start; l < stop; l++ )
          {    int64_t id3 = aligns[l].id2;
               int offset23 = aligns[l].offset;
               Bool rc23 = aligns[l].rc2;
               Bool rc3 = rc2 ^ rc23;
               int offset3 =
                    ( !rc2 ? offset2 + offset23
                           : offset2 - ( offset23 + bases[id3].isize( )
                                 - bases[id2].isize( ) ) );
               if ( offset3 < -bases[id3].isize( ) + K ) continue;
               if ( offset3 > Cols( ) - K ) continue;
               if ( !BinMember( expect, pairs.getPairID(id3) ) ) continue;
               if ( BinMember( offset_id_rc2_orig,
                    make_triple( offset3, id3, rc3  ) ) )
               {    continue;    }
               offset_id_rc2.push( offset3, id3, rc3 );    }    }
     UniqueSort(offset_id_rc2);
     AddToStack( offset_id_rc2, bases, quals, pairs );    }

class offset_datum {

     public:

     offset_datum( ) { }
     offset_datum( const int offset, const int overlap, const int fraglen,
          const int mismatch, const double bits ) : offset(offset), overlap(overlap),
          fraglen(fraglen), mismatch(mismatch), bits(bits) { }

     void Print(const basevector& con1, const basevector& con2, const int verbosity)
     {    if ( verbosity >= 1 )
          {    if ( verbosity >= 2 ) cout << "\n";
               cout << "offset = " << offset << ", overlap = " << overlap
                    << ", fraglen = " << fraglen
                    << ", mismatches = " << mismatch << ", rate = " 
                    << PERCENT_RATIO( 3, mismatch, overlap )
                    << ", bits = " << bits << endl;
               if ( verbosity >= 2 )
               {    align a;
                    a.SetNblocks(1);
                    a.SetGap( 0, 0 );
                    a.SetLength( 0, overlap );
                    int pos1 = ( offset < 0 ? 0 : offset );
                    int pos2 = ( offset < 0 ? -offset : 0 );
                    a.Setpos1(pos1), a.Setpos2(pos2);
                    PrintVisualAlignment( True, cout, con1, con2, a );    }    }    }

     int offset;
     int overlap;
     int fraglen;
     int mismatch;
     double bits;

};

vec<int> GetOffsets1( const readstack& stack1, const readstack& stack2,
     const int verbosity, const int delta_mis,
     const basevector* con1p, const basevector* con2p )
{
     // Announce.

     if ( verbosity >= 1 ) 
     {    cout << "\nconsensus shiftigram for reads " << stack1.Id(0) << ", " 
               << stack2.Id(0) << ", extended lengths " << stack1.Cols( )
               << ", " << stack2.Cols( ) << ":\n";    }

     // Define heuristics.

     const int min_stretch = 8;
     const int w = 20;
     const double min_bits = 25.0;
     const double min_bits_save = 40.0;
     const int wx = 40;
     const int max_ewx = 20;

     // Precomputed log binomial sums.  Note buggy hardcoding of max_overlap
     const int max_overlap = 1000;
     static PrecomputedBinomialSums<max_overlap> gBS(w,0.75);

     // Get consensus;

     double clock = WallClockTime( );
     basevector con1, con2;
     if ( con1p != NULL ) con1 = *con1p;
     else con1 = stack1.Consensus1( );
     if ( con2p != NULL ) con2 = *con2p;
     con2 = stack2.Consensus1( );

     // Check whether consensuses (which define max overlap) fit in 
     // hardcoded binomial sum size above.  But let's not actually assert.

     // ForceAssertLt( Max(con1.isize(), con2.isize()), max_overlap);              
     if ( Max( con1.isize( ), con2.isize( ) ) >= max_overlap ) return vec<int>( );

     // Look for initial offsets.

     vecbasevector con12;
     con12.push_back(con1);
     con12.push_back(con2);
     vec< triple<kmer<min_stretch>,int,int> > kmers_plus;
     MakeKmerLookup3( con12, kmers_plus );                              // sorted kmer list <kmer, con12 index, basevec offset>
     vec<int> doffsets;
     for ( int i = 0; i < kmers_plus.isize( ); i++ )
     {    int j;
          for ( j = i + 1; j < kmers_plus.isize( ); j++ )
          {    if ( kmers_plus[j].first != kmers_plus[i].first ) 
                    break;    }                                         // [i,j) bracket a range of the same kmer
          for ( int l1 = i; l1 < j; l1++ )
          {
               if ( kmers_plus[l1].second != 0 ) continue;              // is l1 kmer on first read?
               for ( int l2 = i; l2 < j; l2++ )
               {    
                    if ( kmers_plus[l2].second != 1 ) continue;         // is l2 kmer on second read? 
                    int o = kmers_plus[l1].third - kmers_plus[l2].third;
                    doffsets.push_back(o);}}                            // great, push out the offset
                    
               i = j - 1;    }
     UniqueSort(doffsets);

     // Go through initial offsets.

     vec< triple<int,int,int> > oom;                                    // this could be more efficient, but not important now
     size_t max_overlap_seen = 0;
     for ( int oi = 0; oi < doffsets.isize( ); oi++ )
     {    int o = doffsets[oi];
          int overlap = 0, mismatch = 0;
          for ( int p1 = 0; p1 < con1.isize( ); p1++ )
          {    int p2 = p1 - o;
               if ( p2 < 0 || p2 >= con2.isize( ) ) continue;
               overlap++;
               if ( con1[p1] != con2[p2] ) mismatch++;    }
          max_overlap_seen = Max( max_overlap_seen, static_cast<size_t>(overlap) );
          oom.push( o, overlap, mismatch );    }

     vec<offset_datum> offsets;
     vec<Bool> bad_window( max_overlap_seen );          // sliding window of regions with suprathreshold errors
     vec<int> sum_errors( max_overlap_seen+1 );         // cumulative sum of error counts along positions in the overlap
     if ( verbosity >= 3 ) cout << "\noffsets passing bits test:\n\n";
     for ( int j = 0; j < oom.isize( ); j++ )
     {    
          int offset = oom[j].first, overlap = oom[j].second;
          int mismatch = oom[j].third;
          int pos1 = ( offset < 0 ? 0 : offset );
          int pos2 = ( offset < 0 ? -offset : 0 );

          // Define crappy subwindows.
          
          // init and fill "sum_errors"
          ForceAssertGt(overlap, 0);
          sum_errors[0] = 0;    // always the number of errors ***prior*** to this position in the overlap
          for ( int i = 1; i <= overlap; ++i ) {
               sum_errors[i] = sum_errors[i-1] + ( ( con1[pos1+i-1] != con2[pos2+i-1] ) ? 1 : 0 );
          }

          // populate "bad_window"
          int errs = 0;                                                

          ForceAssertEq( max_overlap_seen, bad_window.size() );
          for ( size_t i = 0; i < max_overlap_seen; ++i ) {
               bad_window[i] = False;
          }
          
          for ( int m = 0; m <= overlap - wx; m++ )                     
          {    if ( con1[ pos1 + m ] != con2[ pos2 +  m ] ) errs++;
               if ( m >= wx && con1[ pos1 + m - wx ] != con2[ pos2 + m - wx ] )
                    errs--;
               if ( errs >= max_ewx ) bad_window[ Max( 0, m-wx ) ] = True;    
          }
     
          double minp = 0.;        // replaced 1 with 0 = log10(1) as log is now built into b.s. table.
          for ( int start = 0; start < overlap; start++ )
          {
               Bool bad = False;
               for ( int m = 0; !bad && m < w-wx; m++ )
                    if ( bad_window[start+m] ) 
                         bad = True;
//               if ( bad ) continue;

               for ( int n = w; n <= overlap - start; n++ )
               {    
                    if ( n - wx >= 0 && bad_window[ start + n - wx ] ) 
                         bad = True;

                    int k = sum_errors[ start + n ] - sum_errors[ start ];

                    if (bad) break;                          // okay, so if correct, we dump out at first bad window

                    minp = Min( minp, gBS[n][k] );
               } 
          }
          double bitsj = -minp * 10.0 / 6.0;
          if ( bitsj >= min_bits ) 
          {    int fraglen = offset + con2.isize( );
               offsets.push( offset, overlap, fraglen, mismatch, bitsj );
               if ( verbosity >= 3 )
                    offsets.back( ).Print( con1, con2, verbosity );    }    }
     if ( verbosity >= 2 ) cout << TimeSince(clock) << " used\n";

     // Exclude offsets that imply a Q30 mismatch between the founder reads.

     /*
     vec<Bool> to_delete0( offsets.size( ), True );
     for ( int i = 0; i < offsets.isize( ); i++ )
     {    int o = offsets[i].offset;
          Bool bad = False;
          for ( int p1 = 0; p1 < stack1.Cols( ); p1++ )
          {    int p2 = p1 - o;
               if ( p2 >= 0 && p2 < stack2.Cols( ) )
               {    if ( stack1.Qual(0,p1) >= 30 && stack2.Qual(0,p2) >= 30
                         && stack1.Base(0,p1) != stack2.Base(0,p2) )
                    {    if ( verbosity >= 3 )
                         {    cout << "excluding offset " << o << " because it "
                                   << "creates Q30 mismatch "
                                   << "with founders" << endl;    }
                         bad = True;    
                         break;    }    }    }
          if ( !bad ) to_delete0[i] = False;    }
     EraseIf( offsets, to_delete0 );
     */

     // Call a consensus position 'strong' if the founder agrees with it there
     // (and at flanking positions).

     /*
     const int sflank = 5;
     vec<Bool> fval1( con1.size( ), False ), fval2( con2.size( ), False );
     for ( int i = 0; i < stack1.Cols( ); i++ )
     {    if ( stack1.Base(0,i) != con1[i] ) continue;
          int j;
          for ( j = i + 1; j < stack1.Cols( ); j++ )
               if ( stack1.Base(0,j) != con1[j] ) break;
          for ( int k = i + sflank; k < j - sflank; k++ )
               fval1[k] = True;
          i = j - 1;    }
     for ( int i = 0; i < stack2.Cols( ); i++ )
     {    if ( stack2.Base(0,i) != con2[i] ) continue;
          int j;
          for ( j = i + 1; j < stack2.Cols( ); j++ )
               if ( stack2.Base(0,j) != con2[j] ) break;
          for ( int k = i + sflank; k < j - sflank; k++ )
               fval2[k] = True;
          i = j - 1;    }
     */

     // Exclude offsets that create inconsistent columns.

     /*
     vec<int> offsetsw;
     vec<double> bitsw;
     for ( int i = 0; i < offsets.isize( ); i++ )
     {    int o = offsets[i];
          Bool bad = False;
          for ( int p1 = 0; p1 < stack1.Cols( ); p1++ )
          {    int p2 = p1 - o;
               if ( p2 >= 0 && p2 < stack2.Cols( ) )
               {    vec<int> count1(4,0), count2(4,0);
                    for ( int j = 0; j < stack1.Rows( ); j++ )
                    {    if ( stack1.Qual(j,p1) >= 30 )
                              count1[ stack1.Base(j,p1) ]++;    }
                    for ( int j = 0; j < stack2.Rows( ); j++ )
                    {    if ( stack2.Qual(j,p2) >= 30 )
                              count2[ stack2.Base(j,p2) ]++;    }
                    Bool consistent = False;
                    for ( int l = 0; l < 4; l++ )
                         if ( count1[l] > 0 && count2[l] > 0 ) consistent = True;
                    const int min_count = 5;
                    if ( Max(count1) < min_count || Max(count2) < min_count )
                         consistent = True;
                    if ( !consistent ) 
                    {    if ( verbosity >= 3 )
                         {    cout << "dump offset " << o << " for reads "
                                   << stack1.Id(0) << "/" << stack2.Id(0) 
                                   << ", creates inconsistency in "
                                   << "columns " << p1 << "/" << p2 << endl;    }
                         bad = True;
                         break;    }    }    }
          if ( !bad ) 
          {    offsetsw.push_back( offsets[i] );
               bitsw.push_back( bits[i] );    }    }
     offsets = offsetsw;
     bits = bitsw;
     */

     // See if some offsets invalidate other offsets.

     vec< vec<Bool> > val1( offsets.size( ) ), val2( offsets.size( ) );
     for ( int i = 0; i < offsets.isize( ); i++ )
     {    
          /*
          val1[i] = fval1;
          val2[i] = fval2;
          */
          
          val1[i].resize( stack1.Cols( ), False );
          val2[i].resize( stack2.Cols( ), False );

          int o = offsets[i].offset;
          const int flank = 10;
          for ( int p1 = 0; p1 < con1.isize( ); p1++ )
          {    int p2 = p1 - o;
               if ( p2 < 0 || p2 >= con2.isize( ) ) continue;
               if ( con1[p1] != con2[p2] ) continue;
               int low1 = p1, low2 = p2;
               for ( ; p1 < con1.isize( ); p1++ )
               {    p2 = p1 - o;
                    if ( p2 >= con2.isize( ) ) break;
                    if ( con1[p1] != con2[p2] ) break;    }
               int high1 = p1, high2 = p2;

               // At this point we have intervals [low1,high1) on con1 and
               // [low2,high2) on con2 that agree.  We exclude the first and last
               // 'flank' positions, then call the rest validated.

               for ( int q1 = low1 + flank; q1 < high1 - flank; q1++ )
               {    int q2 = q1 - o;
                    val1[i][q1] = val2[i][q2] = True;    }    }    }
     vec< vec<Bool> > invalidates( 
          offsets.size( ), vec<Bool>( offsets.size( ), False ) );
     for ( int i = 0; i < offsets.isize( ); i++ )
     {    int o = offsets[i].offset;
          for ( int p1 = 0; p1 < con1.isize( ); p1++ )
          {    int p2 = p1 - o;
               if ( p2 < 0 || p2 >= con2.isize( ) ) continue;
               if ( con1[p1] == con2[p2] ) continue;
               for ( int j = 0; j < offsets.isize( ); j++ )
               {    
                    // Offset j invalidates offset i if offset j says that 
                    // p1 and p2 are both validated.

                    if ( val1[j][p1] && val2[j][p2] )
                    {    
                         /*
                         if ( verbosity >= 2 )
                         {
                         cout << "offset " << offsets[j].offset << " invalidates "
                              << "offset " << offsets[i].offset << " via positions "
                              << "(" << p1 << "," << p2 << ")\n";
                         }
                         */
                         invalidates[j][i] = True;    }    }    }    }
     vec<Bool> to_delete( offsets.size( ), False );
     for ( int i = 0; i < offsets.isize( ); i++ )
     {    Bool good = True;
          for ( int j = 0; j < offsets.isize( ); j++ )
               if ( invalidates[j][i] ) good = False;
          if ( !good ) continue;
          for ( int j = 0; j < offsets.isize( ); j++ )
               if ( invalidates[i][j] ) to_delete[j] = True;    }
     EraseIf( offsets, to_delete );

     // Big near small: delete small.

     const double min_slope = 2.0;
     const double min_add = 10.0;
     vec<Bool> to_delete2( offsets.size( ), False );
     for ( int i1 = 0; i1 < offsets.isize( ); i1++ )
     for ( int i2 = 0; i2 < offsets.isize( ); i2++ )
     {    if ( to_delete2[i1] ) continue;
          if ( offsets[i2].bits >= min_bits_save ) continue;
          double delta_bits = offsets[i1].bits - offsets[i2].bits;
          if ( delta_bits < min_add ) continue;
          if ( delta_bits / Abs( offsets[i1].offset - offsets[i2].offset ) 
               < min_slope ) 
          {    continue;    }
          to_delete2[i2] = True;    }
     EraseIf( offsets, to_delete2 );

     // Done.

     vec<int> off;
     for ( int i = 0; i < offsets.isize( ); i++ )
          off.push_back( offsets[i].offset );
     if ( verbosity >= 1 ) cout << "\nsurviving offsets: " << printSeq(off) << endl;
     for ( int i = 0; i < offsets.isize( ); i++ )
          offsets[i].Print( con1, con2, verbosity );
     return off;    }

template<int K> void readstack::RefStack( const basevector& ref,
     const vecbasevector& bases, const vecqualvector& quals, 
     const PairsManager& pairs )
{
     vec< triple<kmer<K>,int,int> > kmers_plus;
     vecbasevector all;
     all.push_back(ref);
     all.Append(bases);
     MakeKmerLookup2( all, kmers_plus );
     vec< triple<int,int64_t,Bool> > offset_id2_rc2;
     for ( int64_t i = 0; i < kmers_plus.jsize( ); i++ )
     {    int64_t j, m;
          for ( j = i + 1; j < kmers_plus.jsize( ); j++ )
               if ( kmers_plus[j].first != kmers_plus[i].first ) break;
          for ( m = i; m < j; m++ )
               if ( kmers_plus[m].second > 0 ) break;
          for ( int64_t l1 = i; l1 < m; l1++ )
          for ( int64_t l2 = m; l2 < j; l2++ )
          {    int pos1 = kmers_plus[l1].third, pos2 = kmers_plus[l2].third;
               int64_t id2 = kmers_plus[l2].second - 1; 
               Bool rc2 = ( pos1 < 0 ^ pos2 < 0 );
               if ( pos1 < 0 ) pos1 = -pos1-1;
               if ( pos2 < 0 ) pos2 = -pos2-1;
               if (rc2) pos2 = bases[id2].isize( ) - pos2 - K;
               offset_id2_rc2.push( pos1 - pos2, id2, rc2 );    }
          i = j - 1;    }
     UniqueSort(offset_id2_rc2);
     AddToStack( offset_id2_rc2, bases, quals, pairs );    }

template void readstack::RefStack<20>( const basevector& ref,
     const vecbasevector& bases, const vecqualvector& quals,
     const PairsManager& pairs );
template void readstack::RefStack<40>( const basevector& ref,
     const vecbasevector& bases, const vecqualvector& quals,
     const PairsManager& pairs );

Bool cmp23( const triple<int,int64_t,Bool>& x1, const triple<int,int64_t,Bool>& x2 )
{    if ( x1.second < x2.second ) return True;
     if ( x1.second > x2.second ) return False;
     return x1.third < x2.third;    }

template<int K> void AddPartnersCore( const int top, readstack& stack, 
     const vecbasevector& bases, const vecqualvector& quals, 
     const PairsManager& pairs )
{    vec< pair<int64_t,int> > pids;
     for ( int i = 0; i < stack.Rows( ); i++ )
          pids.push( stack.Pid(i), stack.PairPos(i) );
     UniqueSort(pids);
     vec<int64_t> pid2;
     for ( int i = 0; i < pids.isize( ); i++ )
     {    int j;
          for ( j = i + 1; j < pids.isize( ); j++ )
               if ( pids[j].first != pids[i].first ) break;
          if ( j - i == 1 ) pid2.push_back( pids[i].first );
          i = j - 1;    }
     vec< pair<int64_t,Bool> > candidates;
     for ( int i = 0; i < stack.Rows( ); i++ )
     {    if ( !BinMember( pid2, stack.Pid(i) ) ) continue;
          candidates.push( pairs.getPartnerID( stack.Id(i) ), !stack.Rc2(i) );
          candidates.push( stack.Id(i), stack.Rc2(i) );    }
     UniqueSort(candidates);
     basevector con = stack.Consensus1( );
     vecbasevector all( candidates.size( ) + 1 );
     all[0] = con;
     for ( int i = 0; i < candidates.isize( ); i++ )
     {    all[i+1] = bases[ candidates[i].first ];
          if ( candidates[i].second ) all[i+1].ReverseComplement( );    }
     vec<int64_t> starts;
     starts.push_back(0);
     for ( size_t i = 0; i < all.size( ); i++ )
     {    const basevector& u = all[i];
          starts.push_back( starts.back( ) + Max( 0, u.isize( ) - K + 1 ) );    }
     vec< triple<kmer<K>,int,int> > kmers_plus( starts.back( ) );
     for ( size_t i = 0; i < all.size( ); i++ )
     {    const basevector& u = all[i];
          for ( int j = 0; j <= u.isize( ) - K; j++ )
          {    int64_t r = starts[i] + j;
               kmers_plus[r].first.SetToSubOf( u, j );
               int64_t id2 = ( i == 0 ? 0 : candidates[i-1].first );
               kmers_plus[r].second = id2;
               Bool rc2 = ( i == 0 ? False : candidates[i-1].second );
               kmers_plus[r].third = ( i == 0 ? j : ( !rc2 ? j 
                    : - ( bases[id2].isize( ) - j - K ) - 1 ) );    }    }
     Sort(kmers_plus);
     vec< triple<int,int64_t,Bool> > offset_id2_rc2;
     for ( int64_t i = 0; i < kmers_plus.jsize( ); i++ )
     {    int64_t j, m;
          for ( j = i + 1; j < kmers_plus.jsize( ); j++ )
               if ( kmers_plus[j].first != kmers_plus[i].first ) break;
          for ( m = i; m < j; m++ )
               if ( kmers_plus[m].second > 0 ) break;
          for ( int64_t l1 = i; l1 < m; l1++ )
          for ( int64_t l2 = m; l2 < j; l2++ )
          {    int pos1 = kmers_plus[l1].third, pos2 = kmers_plus[l2].third;
               int64_t id2 = kmers_plus[l2].second; 
               Bool rc2 = ( pos2 < 0 );
               if (rc2) pos2 = -pos2-1;
               if (rc2) pos2 = bases[id2].isize( ) - pos2 - K;
               offset_id2_rc2.push( pos1 - pos2, id2, rc2 );    }
          i = j - 1;    }
     UniqueSort(offset_id2_rc2);
     vec<Bool> o_to_delete( offset_id2_rc2.size( ), False );
     for ( int i = 0; i < stack.Rows( ); i++ )
     {    int p = BinPosition( offset_id2_rc2, make_triple( stack.Offset(i),
               stack.Id(i), stack.Rc2(i) ) );
          if ( p >= 0 ) o_to_delete[p] = True;    }
     EraseIf( offset_id2_rc2, o_to_delete );
     Sort( offset_id2_rc2, cmp23 );
     int nrows_old = stack.Rows( );
     stack.AddToStack( offset_id2_rc2, bases, quals, pairs );
     vec<Bool> suspect( stack.Rows( ), False );

     /*
     readstack stackx(stack);
     vec<int> ids;
     for ( int j = 0; j < stackx.Rows( ); j++ )
          ids.push_back( stackx.Id(j) );
     UniqueSort(ids);
     vec< vec<Bool> > confirmed( ids.size( ) );
     for ( int j = 0; j < ids.isize( ); j++ )
          confirmed[j].resize( bases[ ids[j] ].size( ), False );
     const int flank = 5;
     for ( int j = 0; j < stackx.Rows( ); j++ )
     {    int id = stackx.Id(j);
          int idp = BinPosition( ids, id );
          for ( int i = flank; i < stackx.Cols( ) - flank; i++ )
          {    Bool def = True;
               for ( int l = i - flank; l <= i + flank; l++ )
               {    if ( !stackx.Def(j,l) )
                    {    def = False;
                         break;    }    }
               if ( !def ) continue;
               for ( int k = 0; k < stackx.Rows( ); k++ )
               {    if ( k == j ) continue;
                    Bool agree = True;
                    for ( int l = i - flank; l <= i + flank; l++ )
                    {    if ( stackx.Base(k,l) != stackx.Base(j,l) )
                         {    agree = False;
                              break;    }    }
                    if ( !agree ) continue;
                    int pos = i - stackx.Offset(j);
                    if ( stackx.Rc2(j) ) pos = bases[id].isize( ) - pos - 1;
                    confirmed[idp][pos] = True;    
                    break;    }    }    }
     for ( int j = 0; j < stackx.Rows( ); j++ )
     for ( int i = 0; i < stackx.Cols( ); i++ )
     {    if ( !stackx.Def(j,i) ) continue;
          int id = stackx.Id(j);
          int idp = BinPosition( ids, id );
          int pos = i - stackx.Offset(j);
          if ( stackx.Rc2(j) ) pos = bases[id].isize( ) - pos - 1;
          if ( confirmed[idp][pos] ) stackx.SetQual( j, i, 30 );    }

     const int n = 30;
     for ( int t = 0; t < top; t++ )
     for ( int j = nrows_old; j < stackx.Rows( ); j++ )
     {    for ( int c = 0; c < stackx.Cols( ); c++ )
          {    if ( stackx.Base(j,c) != stackx.Base(t,c) 
                    && stackx.Qual(j,c) >= n && stackx.Qual(t,c) >= n )
               {    suspect[j] = True;
                    break;    }    }
          int64_t id1 = pairs.getPartnerID( stackx.Id(j) );
          vec<int> locs;
          for ( int m = 0; m < nrows_old; m++ )
               if ( stackx.Id(m) == id1 ) locs.push_back(m);
          if ( locs.solo( ) )
          {    int t = locs[0];
               for ( int c = 0; c < stackx.Cols( ); c++ )
               {    if ( stackx.Base(j,c) != stackx.Base(t,c) 
                         && stackx.Qual(j,c) >= n && stackx.Qual(t,c) >= n )
                    {    suspect[j] = True;
                         break;    }    }    }    }
     */

     const int n = 30;
     for ( int t = 0; t < top; t++ )
     for ( int j = nrows_old; j < stack.Rows( ); j++ )
     {    for ( int c = 0; c < stack.Cols( ); c++ )
          {    if ( stack.Base(j,c) != stack.Base(t,c) 
                    && stack.Qual(j,c) >= n && stack.Qual(t,c) >= n )
               {    suspect[j] = True;
                    break;    }    }    }
     
     stack.Erase(suspect);

     if ( stack.Rows( ) == nrows_old ) return;
     AddPartnersCore<K>( top, stack, bases, quals, pairs );    }

void readstack::AddPartners( const int K, const int top, const vecbasevector& bases, 
     const vecqualvector& quals, const PairsManager& pairs )
{    
     if ( K == 32 ) AddPartnersCore<32>( top, *this, bases, quals, pairs );
     else if ( K == 40 ) AddPartnersCore<40>( top, *this, bases, quals, pairs );
     else
     {    FatalErr("AddPartners not implemented for K = " << K);    }    }

namespace
{
    int cappedLength( StackBaseVec::const_iterator itr,
                      StackBaseVec::const_iterator end,
                      long const homopolymerLengthCap )
    {
        int len = 0;
        while ( itr != end )
        {
            char curBase = *itr;
            auto tmp = itr+1;
            while ( tmp != end && *tmp == curBase )
                ++tmp;
            using std::distance;
            len += std::min(homopolymerLengthCap,distance(itr,tmp));
            itr = tmp;
        }
        return len;
    }
}
void readstack::FlagNoise( vec<Bool>& to_delete ) const
{    const int min_glue = 20;
     const int max_homopolymer_in_glue = 10;
     to_delete.assign( Rows( ), False );
     if ( Rows() < 2 ) return; // EARLY RETURN!

     StackBaseVec const& founder = bases_[0];
     size_t nRows = bases_.size();
     for ( size_t idx = 1; idx != nRows; ++idx )
     {
         StackBaseVec const& read = bases_[idx];
         bool readOK = false;
         auto end = read.end();
         auto itr = read.begin();
         auto fItr = founder.begin();
         while ( !readOK && itr != end )
         {
             auto mismatchItrs = std::mismatch(itr,end,fItr);
             auto tmp = mismatchItrs.first;
             if ( cappedLength(itr,tmp,max_homopolymer_in_glue) >= min_glue )
             {
                 readOK = true;
                 break;
             }
             if ( tmp == end )
                 break;
             itr = ++tmp;
             fItr = ++mismatchItrs.second;
         }
         if ( !readOK )
             to_delete[idx] = True;
     }
}

void readstack::IdentifyShifters( vec<Bool>& to_delete ) const
{    const int min_run = 15;
     const int min_err_diff = 5;
     to_delete.assign( Rows( ), False );
     for ( int p1 = 0; p1 < Cols( ); p1++ )
     {    if ( !Def(0,p1) ) break;
          int p2;
          for ( p2 = p1 + 1; p2 < Cols( ); p2++ )
               if ( Base(0,p2) != Base(0,p1) ) break;
          if ( p2 - p1 >= min_run )
          {    for ( int i = 1; i < Rows( ); i++ )
               {    int errs = 0, errsp = 0, errsm = 0;
                    for ( int j = p1; j < Cols( ); j++ )
                    {    if ( !Def(0,j) || !Def(i,j) ) break;
                         if ( Base(0,j) != Base(i,j) ) errs++;    }
                    for ( int j = p1; j < Cols( ) - 1; j++ )
                    {    if ( !Def(0,j) || !Def(i,j+1) ) break;
                         if ( Base(0,j) != Base(i,j+1) ) errsp++;    }
                    for ( int j = p1; j < Cols( ) - 1; j++ )
                    {    if ( j-1 < 0 ) continue;
                         if ( !Def(0,j) || !Def(i,j-1) ) break;
                         if ( Base(0,j) != Base(i,j-1) ) errsp++;    }
                    if ( Max( errs-errsp, errs-errsm ) >= min_err_diff )
                         to_delete[i] = True;    }
               break;    }    }    }

void readstack::Defenestrate( vec<Bool>& to_delete ) const
{    const int width = 10;
     const int min_mult = 2;
     const int min_diffs = 3;
     const int min_comp = 3;

     int n = Rows( ), k = Cols( );
     to_delete.assign( n, False );

     for ( int i = 0; i <= k - width; i += width )
     {    if ( i + width > k ) break;
          vec< vec<char> > X;
          vec<int> ids;
          for ( int j = 0; j < n; j++ )
          {    Bool bad = False;
               for ( int l = i; l < i + width; l++ )
                    if ( Base(j,l) == ' ' ) bad = True;
               if (bad) continue;
               vec<char> x(width);
               for ( int l = i; l < i + width; l++ )
                    x[l-i] = Base(j,l);
               X.push_back(x);
               ids.push_back(j);    }
          SortSync( X, ids );

          int founder = -1;
          for ( int r = 0; r < X.isize( ); r++ )
          {    int s = X.NextDiff(r);
               if ( s - r >= min_mult ) 
               {    int comp = 1;
                    for ( int j = 1; j < width; j++ )
                         if ( X[r][j] != X[r][j-1] ) comp++;
                    if ( comp >= min_comp )
                    {    founder = r;
                         break;    }    }    }
          if ( founder < 0 ) continue;

          for ( int r = 0; r < X.isize( ); r++ )
          {    int s = X.NextDiff(r);
               if ( r != founder && s - r >= min_mult )
               {    int comp = 1;
                    for ( int j = 1; j < width; j++ )
                         if ( X[r][j] != X[r][j-1] ) comp++;
                    if ( comp >= min_comp )
                    {    int diffs = 0;
                         for ( int j = 0; j < width; j++ )
                              if ( X[r][j] != X[founder][j] ) diffs++;
                         if ( diffs >= min_diffs )
                         {    for ( int l = r; l < s; l++ )
                                   to_delete[ ids[l] ] = True;    }    }    }    }    }    }

char readstack::ColumnConsensus1( const int i ) const
{
     // Compute quality score sum for each base.  Count Q0 as 0.1, Q1 as 0.2,
     // and Q2 as 0.2.

     double sum[4]{0.,0.,0.,0.};
     auto qItr = quals_.begin();
     auto bEnd = bases_.end();
     for ( auto bItr=bases_.begin(); bItr != bEnd; ++bItr,++qItr )
     {
         int q = (*qItr)[i];
         if ( q < 0 ) continue;
         double val;
         switch ( q )
         { case 0: val = .1; break;
           case 1: val = .2; break;
           case 2: val = .2; break;
           default: val = q; break; }
         sum[int((*bItr)[i])] += val;
     }
     return std::max_element(sum,sum+4)-sum; }

void readstack::HighQualDiffWindow( vec<Bool>& to_delete ) const
{    const int w = 10;
     const int min_diffs = 3;
     const int min_qsum = 30;
     const int min_qual = 10;
     to_delete.assign( Rows( ), False );
     auto fqItr = quals_[0].begin();
     auto fbItr = bases_[0].begin();
     for ( int c = 0; c <= Cols( ) - w; c++,++fqItr,++fbItr )
     {    if ( fqItr[0] < 0 || fqItr[w-1] < 0 ) continue;
          auto fbEnd = fbItr + w;
          Bool confirmed = False;
          for ( int j = 1; j < Rows( ); j++ )
          {    auto bItr = bases_[j].begin()+c;
               auto qItr = quals_[j].begin()+c;
               auto qEnd = qItr + w;
               if ( std::equal(fbItr,fbEnd,bItr) &&
                       std::find_if(qItr,qEnd,
                               [min_qual](int q){return q<min_qual;})==qEnd )
               {
                   confirmed = True;
                   break;
               }
          }
          if ( !confirmed ) continue;
          for ( int j = 1; j < Rows( ); j++ )
          {    if ( !Def(j,c) || !Def(j,c+w-1) ) continue;
               int diffs = 0, qsum = 0;
               for ( int k = 0; k < w; k++ )
               {    if ( Base(j,c+k) != Base(0,c+k) )
                    {    diffs++;
                         qsum += Qual(j,c+k);    }    }
               if ( diffs >= min_diffs && qsum >= min_qsum )
                    to_delete[j] = True;    }    }    }

void readstack::AddColumnsOnRight( const int n )
{    cols_ += n;
     for ( int i = 0; i < Rows( ); i++ )
     {    bases_[i].resize( cols_, ' ' );
          quals_[i].resize( cols_, -1 );    }    }
