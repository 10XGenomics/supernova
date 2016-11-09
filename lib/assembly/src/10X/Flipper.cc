// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

// Flipper.  Phase bubbles on lines.
//
// Partially phase.  For each line, find sets of bubbles, and the orientations of
// the bubbles in each set, defining a phasing of those bubbles.
//
// First we find the bubbles.
//
// Then we define a set of molecules.  A barcode can correspond to a set of molecules
// if there are large gaps.
//
// Initially, we try to use the molecules to phase all the bubbles.  Each phasing
// may be represented as a choice of state (fw or not) for each bubble.  We start
// with all fw, which of course is almost never true.  Relative to a given phasing,
// each molecule appears as a sequence of pluses and minuses.  We score the molecule
// as Max(pluses,minuses) - Min(pluses,minuses).  The score of a phasing is the sum
// of the scores of each molecule.
//
// We iteratively perturb the phasing in various ways in an attempt to improve its
// score:
// 1. Force each molecule to be in phase.
// 2. Try flipping each bubble.
// 3. Try pivoting at each point, flipping all the orientations to one side.
//
// After the initial phasing, we drop bubbles that phase ambiguously, and break
// the phasing at weak points (uncertain pivots).
//
// Finally, we form phase groups (bounded by uncertain pivots) into giant bubbles.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"
#include "graph/DigraphTemplate.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "paths/long/large/Lines.h"
#include "10X/DfTools.h"
#include "10X/Flipper.h"
#include "10X/LineOO.h"
#include "10X/Super.h"


class BandedMatrix {
public:
     BandedMatrix() : _rows(0), _cols(0) {} // Constructor
     
     ~BandedMatrix() {} //Default destructor

     int rows () { return _rows; }

     int cols () { return _cols; }

     void add_row ( vec<int8_t> row ); // add a row 

     int8_t entry( int i, int j ); // access elements

     void modify( int i, int j, int8_t new_val); // modify elements

     void sort (); // sort by start, end position
     
     vec < pair <int, int> > rse ( ) { return _rse; } // access the start/end
private:
     vec < vec <int8_t> > _M;
     int _rows, _cols;
     vec <pair<int, int>> _rse; // start and end of non-zero elements in each row
};

void BandedMatrix::modify( int i, int j, int8_t new_val ) {
     // Make sure that we are modifying an element that is already non-zero
     ForceAssertGe ( j, 0 );
     ForceAssertLt ( j, _cols );
     if ( j >= _rse[i].first && j < _rse[i].second )
          _M[i][j - _rse[i].first ] = new_val;
}
void BandedMatrix::sort () {
     SortSync( _rse, _M );
}
int8_t BandedMatrix::entry (int i, int j) {
     ForceAssertGe(i, 0);
     ForceAssertGe(j, 0);
     ForceAssertLt(i, _rows);
     ForceAssertLt(j, _cols);
     if ( j < _rse[i].first || j >= _rse[i].second )
          return int8_t(0);
     return _M[i][ j - _rse[i].first ];
}

void BandedMatrix::add_row ( vec <int8_t> row ) {
     int rs = row.size();
     if ( _cols > 0 )
          ForceAssertEq(rs, _cols);
     _cols = rs; // initialize columns if necessary
     _rows++;  // increment rows
     int start = -1, end = -1;
     for ( int i = 0; i != rs; i++ ) {
          if ( row[i] != 0 ) {
               start = i;
               break;
          }
     }
     for ( int i = rs - 1; i >= 0; i-- ) {
          if ( row[i] != 0 ) {
               end = i;
               break;
          }
     }
     ForceAssertEq ( start == -1, end == -1 );
     end++; // increment end to make it non-exclusive
     if ( start == -1 ) {
          start = 0;
          end   = 0;
     }
     _rse.push_back(make_pair(start, end));
     vec <int8_t> data;
     for ( int i = start; i != end; i++ ) {
          data.push_back(row[i]);
     }
     _M.push_back(data);
}

void FixColumns( const String& title, int nmol, int ncol, vec<int>& plus,
     vec<int>& minus, int& goods, int& bads, BandedMatrix & QS, vec<Bool>& fw,
     ostream& out, vec<pair<int,int>> & flm )
{
     for ( int pass = 1; ; pass++ )
     {    out << Date( ) << ": [" << title << "] fixing individual columns, pass "
               << pass << endl;
          int nfix = 0;
          for ( int i = 0; i < ncol; i++ )
          {    int delta_good = 0, delta_bad = 0;
			int first = flm[i].first, last = flm[i].second;	
               for ( int j = first; j < last; j++ )
               {    int p = plus[j], m = minus[j];
                    int g = Max( p, m ), b = Min( p, m );
                    int8_t Qji = QS.entry(j,i);
                    if ( Qji > 0 ) { p--, m++; }
                    if ( Qji < 0 ) { m--, p++; }
                    int g2 = Max( p, m ), b2 = Min( p, m );
                    delta_good += g2 - g;
                    delta_bad += b2 - b;    }
               if ( delta_good > 0 && delta_bad <= delta_good/10 )
               {    fw[i] = !fw[i];
                    nfix++;
                    for ( int j = first; j < last; j++ )
                    {    int8_t Qji = QS.entry (j,i);
                         if ( Qji == 0 ) continue;
                         if ( Qji > 0 )
                         {    plus[j]--;
                              minus[j]++;    }
                         else
                         {    minus[j]--;
                              plus[j]++;     }
                         QS.modify(j, i, -Qji);    }
                    goods += delta_good;
                    bads += delta_bad;    }    }
          if ( nfix == 0 ) break;
          out << Date( ) << ": [" << title << "] fixed " << nfix
               << " columns" << endl;    }
     PRINT2_TO( out, goods, bads );    }

void ForceAssertMatrixEq ( vec<vec<int8_t>> & Q, BandedMatrix & QS ) 
{
     int rows = QS.rows ();
     int cols = QS.cols ();
     ForceAssertEq ( rows, Q.size() );
     for ( int i = 0; i != rows; i++ ) {
          ForceAssertEq ( Q[i].size(), cols );
          for ( int j = 0; j != cols; j++ )
               ForceAssertEq ( Q[i][j], QS.entry(i, j) );
     }
}
void Flipper( const HyperBasevectorX& hb, const vec<int>& inv, const vec<Bool>& dup,
     const vec<Bool>& bad, const vec<int32_t>& bc, digraphE<vec<int>>& D,
     vec<int>& dinv, const vec<vec<vec<vec<int>>>>& dlines,
     const IntIndex& dpaths_index, vec<String>& report, const int test_line )
{
     // Compute ancillary data structures.

     vec<int> Lens( D.E( ), 0 ), linv;
     for ( int d = 0; d < D.E( ); d++ )
     {    if ( D.O(d)[0] < 0 ) continue;
          for ( int j = 0; j < D.O(d).isize( ); j++ )
               Lens[d] += hb.Kmers( D.O(d)[j] );    }
     LineInv( dlines, dinv, linv );

     // Find lines to be processed.

     cout << Date( ) << ": identifying lines" << endl;
     vec<int> llens;
     GetLineLengths( hb, D, dlines, llens );
     vec<int> proc, plens;
     for ( int l = 0; l < dlines.isize( ); l++ )
     {    if ( linv[l] < l ) continue;
          if ( test_line >= 0 && l != test_line ) continue;
          if ( dlines[l].front( )[0][0] == dlines[l].back( )[0][0] ) continue;
          int bubbles = 0;
          const vec<vec<vec<int>>>& L = dlines[l];
          for ( int i = 0; i < L.isize( ); i++ )
          {    if ( i % 2 != 1 || L[i].size( ) != 2 ) continue;
               if ( !L[i][0].solo( ) || !L[i][1].solo( ) ) continue;
               int d1 = L[i][0][0], d2 = L[i][1][0];
               if ( D.O(d1)[0] < 0 || D.O(d2)[0] < 0 ) continue;
               bubbles++;    }
          if ( bubbles >= 2 )
          {    proc.push_back(l);
               plens.push_back( llens[l] );    }    }
     ReverseSortSync( plens, proc );

     // Go through the lines.

     report.resize( dlines.size( ) );
     vec<vec<vec< pair<int,Bool> >>> groups( dlines.size( ) );
     cout << Date( ) << ": main loop, looping over " << proc.size( )
          << " lines, showing 100 dots:" << endl;
     int ndots = 0, stopped = 0;
     #pragma omp parallel for schedule(dynamic, 1)
     for ( int pl = 0; pl < proc.isize( ); pl++ )
     {    int l = proc[pl];
          ostringstream out;

          // Get start positions of line units.

          const vec<vec<vec<int>>> & L = dlines[l];
          vec<int> starts = {0};
          int pos = 0;
          for ( int i = 0; i < L.isize( ); i++ )
          {    vec<int> lens;
               for ( int k = 0; k < L[i].isize( ); k++ )
               {    int len = 0;
                    for ( int l = 0; l < L[i][k].isize( ); l++ )
                    {    int d = L[i][k][l];
                         if ( D.O(d)[0] < 0 ) continue;
                         for ( int m = 0; m < D.O(d).isize( ); m++ )
                              len += hb.Kmers( D.O(d)[m] );    }
                    lens.push_back(len);    }
               Sort(lens);
               int add = 0;
               if ( lens.nonempty( ) ) add = Median(lens);
               pos += add;
               starts.push_back(pos);    }

          // Find bubbles and their positions.  We make X and Y into sync'ed
          // vectors, with one entry per bubble, where X contains the bubble
          // branches and Y contains the line unit indices.

          const int MAX_DIST = 100000;
          vec<int> X, Y, POS;
          vec<Bool> indel;
          pos = 0;
          out << endl << "LINE " << l << endl << endl;
          for ( int i = 0; i < L.isize( ); i++ )
          {    vec<int> lens;
               for ( int k = 0; k < L[i].isize( ); k++ )
               {    int len = 0;
                    for ( int l = 0; l < L[i][k].isize( ); l++ )
                    {    int d = L[i][k][l];
                         if ( D.O(d)[0] < 0 ) continue;
                         for ( int m = 0; m < D.O(d).isize( ); m++ )
                              len += hb.Kmers( D.O(d)[m] );    }
                    lens.push_back(len);    }
               Sort(lens);
               int add = 0;
               if ( lens.nonempty( ) ) add = Median(lens);
               int pos0 = pos;
               pos += add;
               if ( i % 2 != 1 || L[i].size( ) != 2 ) continue;
               if ( !L[i][0].solo( ) || !L[i][1].solo( ) ) continue;
               int d1 = L[i][0][0], d2 = L[i][1][0];
               if ( D.O(d1)[0] < 0 || D.O(d2)[0] < 0 ) continue;
               out << "bubble " << Y.size( ) << " = " << d1 << ", " << d2 << endl;
               X.push_back( d1, d2 );
               Y.push_back(i);
               POS.push_back( pos0, pos0 );
               indel.push_back( Lens[d1] != Lens[d2] );    }

          // Create the data matrix Q.  This consists of one row for each putative
          // molecule.  At a given node, each molecule reports +1 (supports first)
          // or -1 (supports second) or 0 (supports neither or both).

          // BandedMatrix implementation
          BandedMatrix QS;

          vec<int> bars;
          vec< pair<int,int> > M; // { (barcode,node) }
          for ( int i = 0; i < X.isize( ); i++ )
          {    int d = X[i], m = i;
               for ( int xpass = 1; xpass <= 2; xpass++ )
               {    int g = ( xpass == 1 ? d : dinv[d] );
                    for ( int j = 0; j < dpaths_index.Count(g); j++ )
                    {    int64_t id = dpaths_index.Val( g, j );
                         if ( bad[id/2] || bc[id] == 0 ) continue;
                         M.push( bc[id], i );    }    }    }
          UniqueSort(M);
          int ncol = X.size( ) / 2;
          out << "\nhave " << ncol << " bubbles" << endl;
          for ( int i = 0; i < M.isize( ); i++ )
          {    int j, b = M[i].first;
               for ( j = i + 1; j < M.isize( ); j++ )
                    if ( M[j].first != b ) break;
               vec<int> x;
               for ( int k = i; k < j; k++ )
               {    if ( k > i && M[k].second/2 == M[k-1].second/2 ) continue;
                    if ( k < j - 1 && M[k].second/2 == M[k+1].second/2 ) continue;
                    int pos1 = 0;
                    if ( x.nonempty( ) ) pos1 = POS[ x.back( ) ];
                    int pos2 = POS[ M[k].second ];
                    if ( x.nonempty( ) && pos2 - pos1 > MAX_DIST )
                    {    if ( x.size( ) > 1 )
                         {    vec<int8_t> y( ncol, 0 );
                              for ( int m = 0; m < x.isize( ); m++ )
                              {    if ( x[m] % 2 == 0 ) y[ x[m]/2 ] = +1;
                                   else y[ x[m]/2 ] = -1;    }
                              QS.add_row(y);
                              bars.push_back(b);    }
                         x.clear( );    }
                    x.push_back( M[k].second );    }
               if ( x.size( ) > 1 )
               {    vec<int8_t> y( ncol, 0 );
                    for ( int m = 0; m < x.isize( ); m++ )
                    {    if ( x[m] % 2 == 0 ) y[ x[m]/2 ] = +1;
                         else y[ x[m]/2 ] = -1;    }
                    QS.add_row(y);
                    bars.push_back(b);     }
               i = j - 1;    }
          
          int nmol = QS.rows( );
          out << "observe " << nmol << " molecules" << endl;

          // Sort molecules by start, end position
          QS.sort();

          // Keep track of start/end of each molecule (in bubble coordinates)
          // molecule extends from start to end-1 both included.
		vec <pair<int, int>>  mse = QS.rse ( ); //molecule-start-end = mse
          

		// For each bubble, first and last molecule that cover it
		vec <pair <int, int>> flm (ncol); //first-last-molecule
		for ( int i = 0; i != ncol; i++ ) {
			int first = -1, last=-1;
			for ( int j = 0; j != nmol; j++ ) {
				int8_t Qji = QS.entry(j,i);
                    if ( Qji != 0 ) {
					first=j;
					break;
				}
			}
			for ( int j = nmol-1; j >= 0; j-- ) {
				int8_t Qji = QS.entry(j,i);
                    if (Qji != 0 ) {
					last=j+1;
					break;
				}
			}
			if ( first == -1 || last == -1 ) {
                    first = 0;
                    last  = 0;
               }
               ForceAssertGe(first, 0);
			ForceAssertGe(last, 0);
			flm[i] = make_pair(first, last);
		}

          // Define initial state.

          vec<Bool> fw( ncol, True );

          // Track + and - for each molecule.

          vec<int> plus( nmol, 0 ), minus( nmol, 0 );
          for ( int i = 0; i < nmol; i++ ) {
			int start = mse[i].first, end = mse[i].second;
			for ( int j = start; j < end; j++ )
			{    int8_t Qij = QS.entry(i,j);
                    if ( Qij == 0 ) continue;
				if ( Qij > 0 ) plus[i]++;
				else minus[i]++;    }
		}
          // Compute initial goods and bads.

          int goods = 0, bads = 0;
          for ( int i = 0; i < nmol; i++ )
          {    goods += Max( plus[i], minus[i] );
               bads += Min( plus[i], minus[i] );    }
          PRINT2_TO( out, goods, bads );

          // Go through the molecules and try to rectify them.

          out << Date( ) << ": start rectification" << endl;
          for ( int i = 0; i < nmol; i++ )
          {	int start = mse[i].first, end = mse[i].second;    
			vec<int> flips;
               if ( plus[i] >= minus[i] ) // try changing all minuses to pluses
               {    for ( int j = start; j < end; j++ ) {
                         int8_t Qij = QS.entry(i,j);
                         if ( Qij < 0 ) flips.push_back(j); }    }
               else // try changing all pluses to minuses
               {    for ( int j = start; j < end; j++ ) {
                         int8_t Qij = QS.entry(i, j);
                         if ( Qij > 0 ) flips.push_back(j); }    }
               if ( flips.empty( ) ) continue;
               int goodsx = goods, badsx = bads;
               int first=flm[flips[0]].first, last=flm[flips[0]].second;
			for ( int k = 1; k < flips.isize(); k++ ) {
				int new_first = flm[flips[k]].first, new_last = flm[flips[k]].second;
				if (first > new_first)
					first = new_first;
				if (last < new_last)
					last = new_last;
			}
               for ( int j = first; j < last; j++ )
               {    
                    int plusx = plus[j], minusx = minus[j];
                    for ( int k = 0; k < flips.isize( ); k++ )
                    {    int m = flips[k];
                         int8_t Qjm = QS.entry(j, m);
                         if ( Qjm == 0 ) continue;
                         if ( Qjm > 0 )
                         {    plusx--;
                              minusx++;    }
                         else
                         {    minusx--;
                              plusx++;    }    }
                    goodsx += Max( plusx, minusx ) - Max( plus[j], minus[j] );
                    badsx += Min( plusx, minusx ) - Min( plus[j], minus[j] );    }
               if ( goodsx - badsx > goods - bads )
               {    goods = goodsx, bads = badsx;
                    for ( int k = 0; k < flips.isize( ); k++ )
                    {    int m = flips[k];
                         int first_index=flm[m].first, last_index=flm[m].second;
					for ( int j = first_index; j < last_index; j++ )
                         {    int8_t Qjm = QS.entry(j, m);
                              if ( Qjm == 0 ) continue;
                              if ( Qjm > 0 )
                              {    plus[j]--;
                                   minus[j]++;   }
                              else
                              {    minus[j]--;
                                   plus[j]++; } 
                              QS.modify(j, m, -Qjm); }
                         fw[m] = !fw[m];    }    }    }
          out << Date( ) << ": done" << endl;
          PRINT2_TO( out, goods, bads );

          // Look for "pivots", where flipping all the bubbles to one side ups the
          // score.

          out << Date( ) << ": pivoting" << endl;
          vec<int> plus_left( nmol, 0 ), minus_left( nmol, 0 );
          vec<int> plus_right( nmol, 0 ), minus_right( nmol, 0 );
          for ( int i = 0; i < ncol - 1; i++ )
          {    int first_index=flm[i].first, last_index=flm[i].second;
			// keep it as 0, nmol for now.
			// because some data structures are being initialized.
			for ( int j = 0; j < nmol; j++ ) 
			{    int8_t Qji = QS.entry(j, i);
                    if ( Qji > 0 ) plus_left[j]++;
                    if ( Qji < 0 ) minus_left[j]++;
                    plus_right[j] = plus[j] - plus_left[j];
                    minus_right[j] = minus[j] - minus_left[j];    }

               // Now for each molecule j, (plus_left[j], minus_left[j]) are the
               // number of pluses and minuses up through column i, and
               // (plus_right[j], minus_right[j]) are the number that some after it.
               //
               // Compute the number of goods and bads, were we to pivot between
               // columns i and i+1.

               int goodp = 0, badp = 0;
               for ( int j = 0; j < nmol; j++ )
               {    int plusj = plus[j] - plus_left[j] + minus_left[j];
                    int minusj = minus[j] + plus_left[j] - minus_left[j];
                    goodp += Max( plusj, minusj );
                    badp += Min( plusj, minusj );    }

               // And is it better to flip?

               int advantage = (goodp-badp) - (goods-bads);
               if ( advantage > 0 )
               {    out << "flipping columns 0 through " << i << endl;
                    for ( int k = 0; k <= i; k++ ) fw[k] = !fw[k];
                    for ( int j = 0; j < nmol; j++ )
                    {    for ( int k = 0; k <= i; k++ )
                              QS.modify(j, k, -QS.entry(j, k));
                         plus[j] = plus[j] - plus_left[j] + minus_left[j];
                         minus[j] = minus[j] + plus_left[j] - minus_left[j];
                         swap( plus_left[j], minus_left[j] );    }
                    // DPRINT2_TO( out, i, advantage );
                    goods = goodp, bads = badp;    }    }
          PRINT2_TO( out, goods, bads );

          // Go through individual columns and try to fix them.

          FixColumns( "alpha", nmol, ncol, plus, minus, goods, bads, QS, fw, out, flm );

          // Try reverse rectification.

          out << Date( ) << ": start -rectification" << endl;
          for ( int i = 0; i < nmol; i++ )
          {    int start = mse[i].first, end = mse[i].second;
			vec<int> flips;
               if ( plus[i] <= minus[i] ) // try changing all minuses to pluses
               {    for ( int j = start; j < end; j++ )
                         if ( QS.entry(i, j) < 0 ) flips.push_back(j);    }
               else // try changing all pluses to minuses
               {    for ( int j = start; j < end; j++ )
                         if ( QS.entry(i, j) > 0 ) flips.push_back(j);    }
               if ( flips.empty( ) ) continue;
               int goodsx = goods, badsx = bads;
			int first=flm[flips[0]].first, last=flm[flips[0]].second;
			for ( int k = 1; k < flips.isize(); k++ ) {
				int new_first = flm[flips[k]].first, new_last = flm[flips[k]].second;
				if (first > new_first)
					first = new_first;
				if (last < new_last)
					last = new_last;
			}
               for ( int j = first; j < last; j++ )
               {    int plusx = plus[j], minusx = minus[j];
                    for ( int k = 0; k < flips.isize( ); k++ )
                    {    int m = flips[k];
                         int8_t Qjm = QS.entry(j, m);
                         if ( Qjm == 0 ) continue;
                         if ( Qjm > 0 )
                         {    plusx--;
                              minusx++;    }
                         else
                         {    minusx--;
                              plusx++;    }    }
                    goodsx += Max( plusx, minusx ) - Max( plus[j], minus[j] );
                    badsx += Min( plusx, minusx ) - Min( plus[j], minus[j] );    }
               if ( goodsx - badsx > goods - bads )
               {    goods = goodsx, bads = badsx;
                    out << "yikes!" << endl;
                    for ( int k = 0; k < flips.isize( ); k++ )
                    {    int m = flips[k];
					int first_m = flm[m].first, last_m = flm[m].second;
                         for ( int j = first_m; j < last_m; j++ )
                         {    int8_t Qjm = QS.entry(j, m);
                              if ( Qjm == 0 ) continue;
                              if ( Qjm > 0 )
                              {    plus[j]--;
                                   minus[j]++;   }
                              else
                              {    minus[j]--;
                                   plus[j]++;  }
                              QS.modify(j, m, -Qjm);   }
                         fw[m] = !fw[m];    }    }    }
          out << Date( ) << ": done" << endl;
          PRINT2_TO( out, goods, bads );

          // Go through individual columns and try to fix them (again).

          FixColumns( "beta", nmol, ncol, plus, minus, goods, bads, QS, fw, out, flm );

          // Now find ugly bubbles, marked in "ugly".

          out << Date( ) << ": finding ugly bubbles" << endl;
          Bool verbose = False;
          const double MIN_GOOD_BAD_RATIO = 4.0;
          vec<Bool> ugly( ncol, False );
          for ( int i = 0; i < ncol; i++ )
          {    int goodsi = 0, badsi = 0;
			int first = flm[i].first, last = flm[i].second;
               for ( int j = first; j < last; j++ )
               {    int8_t Qji = QS.entry(j, i);
                    if ( Qji == 0 ) continue;
                    if ( plus[j] >= minus[j] )
                    {    if ( Qji < 0 ) badsi++;
                         else  goodsi++;    }
                    else
                    {    if ( Qji > 0 ) badsi++;
                         else  goodsi++;    }    }
               if ( double(goodsi) / double( Max( 1, badsi ) ) < MIN_GOOD_BAD_RATIO )
               {    ugly[i] = True;
                    int d1 = X[2*i], d2 = X[2*i+1];
                    if (verbose)
                         PRINT6_TO( out, i, int(indel[i]), goodsi, badsi, d1, d2 );
                    for ( int j = first; j < last; j++ )
                    {    int8_t Qji = QS.entry(j, i);
                         if ( Qji == 0 ) continue;
                         if ( Qji > 0 ) plus[j]--;
                         else minus[j]--;
                         QS.modify(j, i, 0);    }    }    }
          out << Date( ) << ": found " << Sum(ugly) << " ugly bubbles" << endl;

          // Recompute goods and bads.

          goods = 0, bads = 0;
          for ( int i = 0; i < nmol; i++ )
          {    goods += Max( plus[i], minus[i] );
               bads += Min( plus[i], minus[i] );    }
          PRINT2_TO( out, goods, bads );
          out << "bad fraction = " << PERCENT_RATIO( 3, bads, goods + bads ) << endl;

          // Go through individual columns and try to fix them (again again).

          FixColumns( "gamma", nmol, ncol, plus, minus, goods, bads, QS, fw, out, flm );

          // Recompute goods and bads.

          goods = 0, bads = 0;
          for ( int i = 0; i < nmol; i++ )
          {    goods += Max( plus[i], minus[i] );
               bads += Min( plus[i], minus[i] );    }
          PRINT2_TO( out, goods, bads );
          out << "bad fraction = " << PERCENT_RATIO( 3, bads, goods + bads ) << endl;


          // Define weak pivot points that we should break at.  This info goes
          // into "breakafter".

          vec<int> breakafter;
          const int MAX_PIVOT_OK = -20;
          out << Date( ) << ": break at weak pivot points" << endl;
          plus_left.resize_and_set( nmol, 0 );
          minus_left.resize_and_set( nmol, 0 );
          plus_right.resize_and_set( nmol, 0 );
          minus_right.resize_and_set( nmol, 0 );
          out << "\nweak pivots:\n";
          for ( int i = 0; i < ncol - 1; i++ )
          {    // keep the limits of this loop intact.
	       for ( int j = 0; j < nmol; j++ )
               {    int8_t Qji = QS.entry(j, i);
                    if ( Qji > 0 ) plus_left[j]++;
                    if ( Qji < 0 ) minus_left[j]++;
                    plus_right[j] = plus[j] - plus_left[j];
                    minus_right[j] = minus[j] - minus_left[j];    }

               // Now for each molecule j, (plus_left[j], minus_left[j]) are the
               // number of pluses and minuses up through column i, and
               // (plus_right[j], minus_right[j]) are the number that some after it.
               //
               // Compute the number of goods and bads, were we to pivot between
               // columns i and i+1.

               int goodp = 0, badp = 0;
               for ( int j = 0; j < nmol; j++ )
               {    int plusj = plus[j] - plus_left[j] + minus_left[j];
                    int minusj = minus[j] + plus_left[j] - minus_left[j];
                    goodp += Max( plusj, minusj );
                    badp += Min( plusj, minusj );    }

               // Is there uncertainty about whether we should flip?

               int advantage = (goodp-badp) - (goods-bads);
               if ( advantage > 2 * MAX_PIVOT_OK ) PRINT2_TO( out, i, advantage );
               if ( advantage > MAX_PIVOT_OK )
               {    out << "breaking" << endl;
                    breakafter.push_back(i);    }    }

          // Display the final matrix.

          out << "\nfinal matrix:\n";
          int ndigits = int(ceil(log10(ncol)));
          for ( int j = 0; j < nmol; j++ )
          {    out << "\nmolecule " << j << " = barcode " << bars[j];
               int low, high;
               for ( low = 0; low < ncol; low++ )
                    if ( QS.entry(j, low) != 0 ) break;
               low = 20 * (low/20);
               for ( high = ncol - 1; high >= 0; high-- )
                    if ( QS.entry(j, high) != 0 ) break;
               high++;
               high = Min( 20 * ((high+20-1)/20), ncol );
               for ( int i = low; i < high; i++ )
               {    if ( i % 20 == 0 )
                    {    out << "\n[" << i << "]";
                         int d = ( i == 0 ? 1 : ceil(log10(i+1)) );
                         for ( int l = d; l < ndigits; l++ )
                              out << " ";    }
                    if ( i % 10 == 0 ) out << "    ";
                    out << " ";
                    if ( ugly[i] ) out << "-";
                    else
                    {    int8_t Qji = QS.entry(j, i);
                         if ( Qji == +1 ) out << "O";
                         if ( Qji ==  0 ) out << ".";
                         if ( Qji == -1 ) out << "";    }    }
               out << "\n";    }

          // Define bubble merge groups.

          for ( int i = 0; i < ncol; i++ )
          {    if ( ugly[i] ) continue;
               vec< pair<int,Bool> > g = { make_pair( Y[i], fw[i] ) };
               int j;
               for ( j = i + 1; j < ncol; j++ )
               {
                    /*
                    Bool at_gap = False;
                    for ( int z = Y[j-1] + 2; z <= Y[j] - 2; z += 2 )
                    {    if ( L[z].solo( ) && L[z][0].empty( ) )
                         {    at_gap = True;
                              break;    }
                         for ( int r = 0; r < L[z].isize( ); r++ )
                         for ( int s = 0; s < L[z][r].isize( ); s++ )
                              if ( D.O( L[z][r][s] )[0] < 0 ) at_gap = True;    }
                    if (at_gap) break;
                    */

                    if ( Member( breakafter, j-1 ) ) break;
                    if ( !ugly[j] ) g.push_back( make_pair( Y[j], fw[j] ) );    }
               if ( g.size( ) > 1 ) groups[l].push_back(g);
               i = j - 1;    }

          // Save report.

          report[l] = out.str( );

          #pragma omp critical
          {    MakeDots( stopped, ndots, proc.size( ) );    }    }

     // Edit the graph.

     cout << Date( ) << ": editing graph" << endl;
     vec<int> DELS, to_left, to_right;
     D.ToLeft(to_left), D.ToRight(to_right);
     for ( int l = 0; l < dlines.isize( ); l++ )
     {    const vec<vec<vec<int>>>& L = dlines[l];
          for ( int j = 0; j < groups[l].isize( ); j++ )
          {    const vec< pair<int,Bool> >& g = groups[l][j];

               // Find the edges in the line chunk.

               vec<int> edges;
               int i1 = g.front( ).first, i2 = g.back( ).first;
               for ( int i = i1; i <= i2; i++ )
               {    for ( int r = 0; r < L[i].isize( ); r++ )
                         edges.append( L[i][r] );
                    if ( i < i2 && L[i].nonempty( ) && L[i][0].nonempty( ) )
                    {    int f = L[i][0][0];
                         int v = to_right[f];
                         if ( D.From(v).nonempty( ) )
                              edges.push_back( D.IFrom(v,0) );    }    }
               UniqueSort(edges);

               // Clone the line chunk.  Ditto for rc.

               digraphE<vec<int>> G( digraphE<vec<int>>::COMPLETE_SUBGRAPH_EDGES,
                    D, edges, to_left, to_right );
               vec<int> redges(edges);
               // for ( auto& d : redges ) d = dinv[d];
               for ( int i = 0; i < edges.isize( ); i++ )
                    redges[i] = dinv[ edges[i] ];
               digraphE<vec<int>> RG( digraphE<vec<int>>::COMPLETE_SUBGRAPH_EDGES,
                    D, redges, to_left, to_right );

               vec<int> gto_left, gto_right;
               G.ToLeft(gto_left), G.ToRight(gto_right);

               // Add to the graph and connect the ends.

               int N = D.N( ), E = D.E( );
               int E0 = E;
               D.Append(G);
               int d1 = L[i1][0][0], d2 = L[i2][0][0];
               int p1 = BinPosition( edges, d1 ), p2 = BinPosition( edges, d2 );
               int v1 = gto_left[p1], v2 = gto_right[p2];
               D.TransferEdges( v1 + N, to_left[d1] );
               D.TransferEdges( v2 + N, to_right[d2] );

               // Delete the edges that are "phased out".

               for ( int j = 0; j < 2; j++ )
               {    for ( int m = 0; m < g.isize( ); m++ )
                    {    int p = g[m].first, q = ( j == 0 ^ g[m].second ) ? 0 : 1;
                         int d = L[p][q][0];
                         if ( j == 0 ) DELS.push_back(d);
                         else
                         {    int p = BinPosition( edges, d );
                              DELS.push_back( E + p );    }    }    }

               // Now repeat all of the above for the rc.

               vec<int> rgto_left, rgto_right;
               RG.ToLeft(rgto_left), RG.ToRight(rgto_right);
               N = D.N( ), E = D.E( );
               D.Append(RG);
               int rd1 = dinv[d1], rd2 = dinv[d2];
               int rp1 = Position( redges, rd1 ), rp2 = Position( redges, rd2 );
               int rv1 = rgto_right[rp1], rv2 = rgto_left[rp2];
               D.TransferEdges( rv1 + N, to_right[rd1] );
               D.TransferEdges( rv2 + N, to_left[rd2] );
               for ( int j = 0; j < 2; j++ )
               {    for ( int m = 0; m < g.isize( ); m++ )
                    {    int p = g[m].first, q = ( j == 0 ^ g[m].second ) ? 0 : 1;
                         int d = L[p][q][0];
                         if ( j == 0 ) DELS.push_back( dinv[d] );
                         else
                         {    int p = BinPosition( edges, d );
                              DELS.push_back( E + p );    }    }    }

               // Update the involution.

               for ( int i = 0; i < edges.isize( ); i++ )
                    dinv.push_back( E + i );
               for ( int i = 0; i < edges.isize( ); i++ )
                    dinv.push_back( E0 + i );    }    }

     // Clean up.

     cout << Date( ) << ": clean up" << endl;
     D.DeleteEdgesParallel(DELS);
     RemoveUnneededVertices( D, dinv );
     CleanupCore( D, dinv );
     Validate( hb, inv, D, dinv );    }
