// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
//
//

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS
//

#ifndef _10X_LINEGRAPHOPS_H
#define _10X_LINEGRAPHOPS_H

#include "Vec.h"
#include "Equiv.h"
#include "paths/HyperBasevector.h"
#include "10X/LineOO.h"
#include "10X/CleanThe.h"
#include "10X/PlaceReads.h"
#include "10X/Super.h"
#include "paths/long/large/Lines.h"
#include "10X/ThreadedLogger.h"
#include "10X/Heuristics.h"
#include "10X/MakeLocalsTools.h"


class LineGraphOps {
public:
     LineGraphOps( 
               /* required, static, input-only */
               HyperBasevectorX const& hb, vec<int> const& inv, 
               const ReadPathVecX & pathsx, vec<Bool> const& dup, 
               vec<int64_t> const& bci, vec< triple<int,int,int>> const& qept, 
               VecIntVec const& ebcx, 
               /* required, may be modified */
               digraphE<vec<int>>& D, vec<int>& dinv,
               vec<vec<vec<vec<int>>>>& dlines, ReadPathVec& dpaths, 
               VecULongVec& dpaths_index,
               /* optional, will be built if needed */
               vec<double>* pcov = nullptr, vec<vec<pair<int,int>>>* plhood = nullptr,
               MasterVec<SerfVec<pair<int,int>>>* plbpx = nullptr,
               /* optional */
               ThreadedLogger* plogger = nullptr )
                : _hb(hb), _inv(inv), _pathsx(pathsx), _dup(dup), 
                  _bci(bci), _qept(qept), _ebcx(ebcx), 
                  _D(D), _dinv(dinv), _dlines(dlines), 
                  _dpaths(dpaths), _dpaths_index(dpaths_index), 
                  _cov(manage(pcov)), _lhood(manage(plhood)), _lbpx(manage(plbpx)),
                  _plog(plogger), _debug(false), _we_own_plog(false), 
                  _ple(nullptr), _plreps(nullptr) {
          if ( _plog == nullptr ) {
               _we_own_plog = true;
               _plog = new ThreadedLogger;
          }
          Build();
     }

     template <class T>
     T& manage( T* p ) {
          if ( p == nullptr ) {
               p = new T;
               auto sp = static_cast<void*>(p);
               _we_own[sp] = true;
          } else {
               auto sp = static_cast<void*>(p);
               _we_own[sp] = false;
          }

          return *p;
     }

     template <class T>
     void unmanage( T* p ) {
          if ( p != nullptr && _we_own[static_cast<void*>(p)] ) delete p;
     }

     ~LineGraphOps() { 
          delete _ple; 
          delete _plreps; 
          if (_we_own_plog) delete _plog;  
          unmanage(&_cov);
          unmanage(&_lhood);
          unmanage(&_lbpx);
     }

     void Build() {
          cout << Date() << ": LineGraphOps Build" << endl;
          BuildLineGraph( _D, _dinv, _dlines, _LG );
          _LG.ToLeft( _lg_to_left ), _LG.ToRight( _lg_to_right );
          ComputeDlens( _hb, _D, _dlens );
          GetLineLengths( _dlens, _dlines, _llens );
          MakeTol( _D, _dlines, _tol );
          _D.ToRight(_to_right);
          _D.ToLeft(_to_left);
          LineInv(_dlines, _dinv, _linv);
          cout << Date() << ": LineGraphOps Build done" << endl;
     }

     void RebuildGraphChanged(
               const int max_cell_paths = MAX_CELL_PATHS_EVALUATION,
                     const int max_cell_depth = MAX_CELL_DEPTH_EVALUATION,
                        const vec<Bool> pmask = vec<Bool>()  ) {      
          cout << Date() << ": LineGraphOps RebuildGraphChanged" << endl;
          _lbpx.clear();
          _lhood.clear();
          _cov.clear();
          FindLinesFixed( _D, _dinv, _dlines, max_cell_paths, max_cell_depth );
          PlaceReadsMasked( _hb, _D, _dup, _pathsx, pmask, _dpaths );
          invert(_dpaths, _dpaths_index, _D.E() );
          cout << Date() << ": LineGraphOps Done" << endl;
          Build();
     }

     void NeedCov() { 
          #pragma omp critical
          if (_cov.size() == 0 ) this->UpdateCov(); 
     }
     void NeedLhood() { 
          #pragma omp critical
          if ( _lhood.size() == 0 ) this->UpdateLhood(); 
     }
     void NeedLbpx() { 
          #pragma omp critical
          if ( _lbpx.size() == 0 ) this->UpdateLbpx(); 
     }

     vec<double> const& GetCov() { NeedCov(); return _cov; }
     vec<vec<pair<int,int>>> const& GetLhood() { NeedLhood(); return _lhood; }
     MasterVec<SerfVec<pair<int,int>>> const& GetLbpx() { NeedLbpx(); return  _lbpx; }
     vec<int> const& GetLinv() { return _linv; }


     void SetDebug( bool debug ) { _debug = debug; }
     void Debug() { SetDebug(true); }
     void NoDebug() { SetDebug(false); }
     void DebugSync() { _plog->join(); }     // caution -- don't use while threads are writing

     std::ostringstream& Log() { return _plog->get(); }


     void UpdateLhood() {
          cout << Date() << ": \tconstruct lhood" << endl;
          LineProx( _hb, _inv, _ebcx, _D, _dinv, _dlines, _qept, _lhood );
     }

     void UpdateLbpx() {
          cout << Date() << ": \tconstruct lbpx" << endl;
          vec<vec<pair<int,int>>> lbp;
          vec<int> bc;
          BciBc( _bci, bc );
          IntIndex dpaths_index( _dpaths, _D.E( ) );
          BarcodePos( bc, _hb, _D, _dinv, _dpaths, dpaths_index, _dlines, lbp, 50000 );
          _lbpx.clear();
          for ( auto x : lbp ) {    
               SerfVec<pair<int,int>> y( x.size( ) );
               for ( int j = 0; j < x.isize( ); j++ )
                    y[j] = x[j];
               _lbpx.push_back(y);    
          }
     }


     void UpdateCov() {
          vec<vec<pair<int,int>>> lbp;
          vec<int> bc;
          BciBc( _bci, bc );
          IntIndex dpaths_index( _dpaths, _D.E( ) );
          BarcodePos( bc, _hb, _D, _dinv, _dpaths, dpaths_index, _dlines, lbp, 0 );
          MasterVec<SerfVec<pair<int,int>>> lbpx;  // barcode positions on lines
          for ( auto x : lbp ) {    
               SerfVec<pair<int,int>> y( x.size( ) );
               for ( int j = 0; j < x.isize( ); j++ )
                    y[j] = x[j];
               lbpx.push_back(y);    
          }
          cout << Date() << ": \tcopy number statistics" << endl;
          vec<int> kmers(_hb.E());
          #pragma omp parallel for
          for ( int e = 0; e < kmers.isize( ); e++ )
               kmers[e] = _hb.Kmers(e);
          LineCN( kmers, lbpx, _D, _dlines, _llens, _cov );
     }

     void BciBc( vec<int64_t> const& bci, vec<int32_t>& bc )
     {
          bc.resize( bci.back() );
          int32_t j = 0;
          for ( int64_t i = 0; i < bc.isize(); ++i )  {
               if ( i >= bci[j+1] ) j++;
               bc[i] = j;
          }
     }

     size_t InitTangles( const int MAX_SHORT = 500, const int MIN_LONG = 4000 );
     vec<int> const& GetTangleReps() { return *_plreps; };
     bool FindTangle( const int L, vec<int>& ins, vec<int>& outs, vec<int>& tangle );

     int GetTangleRep( int index ) { 
          ForceAssert( _plreps != nullptr );
          ForceAssertLt( index, _plreps->size() );
          return (*_plreps)[index]; 
     }

     struct OOScore {
          int line;
          int len;
          int bc;
          double cn;
          int side;
          double score;
          double scaled;

          bool operator<( OOScore const& b ) const {
               auto const& a=*this;
               if ( a.side < b.side ) return true;
               else if ( a.side > b.side ) return false;
               else if ( a.score < b.score ) return true;
               return false;
          }

          friend ostream& operator<<( ostream& os, const OOScore& s ) {
               auto flags = os.setf(ios::showpoint|ios::fixed);
               os << setw(2) << "L " << left << setw(8) << s.line << setw(1) << " " 
                   << setw(8) << " len kb " << right << setprecision(3) << setw(9) << ToString(s.len/1000.) << left 
                   << setw(6) << "  bc " << right << setw(5) << s.bc <<  setw(6) << "  cn " << setw(5) 
                   << right << setprecision(2) << s.cn << left << setw(10) 
                   << ( s.side == 0 ? "  center" : (s.side < 0 ? "  left" : "  right") )
                   << setw(6) << "  raw " << right << setw(10) << s.score << left << setw(8) << "    adj " 
                   << right << setw(10) << s.scaled;
               os.setf(flags);
               return os;
          }

     };

     vec<OOScore> ProbeLhood( int L, double S1 = 40, bool single = false );

     vec<int> const& LgToLeft() { return _lg_to_left; }
     vec<int> const& LgToRight() { return _lg_to_right; }

     digraphE<int> const& Lg() { return _LG; }


     bool Reach( int L, vec<int>& reach );
     int LeftmostRightNeighbor( int L, vec<int>* preach = nullptr );
     bool IsInvTangle( vec<int> const& ins, vec<int> const& outs, vec<int> const& tangle );
     void ScoreInvAndResolve(vec<pair<int,int>> const& pivots);
     int BarcodePairSupport( pair<int,int> const& lpair );

     vec<pair<int,int>> FindInvPivots(); 
     void MakeRipEdit( int L1, int L2, vec<int> lsr );

     vec<int> GetLongLines( int MIN_LONG = -1) {
          if ( MIN_LONG == -1 ) MIN_LONG = _min_long;
          ForceAssertGt(_llens.size(), 0 );
          ForceAssertGt(_dlines.size(), 0 );
          vec<int> longs;
          int N = _dlines.size();
          for ( int L = 0; L < N; ++L ) {
               if ( _llens[L] >= MIN_LONG )
                    longs.push_back(L);
          }
          return longs;
     }

private:
     HyperBasevectorX const&            _hb;
     vec<int> const&                    _inv;
     ReadPathVecX const&                _pathsx;
     vec<Bool> const&                   _dup;
     vec<int64_t> const&                _bci;
     vec< triple<int,int,int>> const&   _qept;
     VecIntVec const&                   _ebcx;
     digraphE<vec<int>>&                _D;
     vec<int>&                          _dinv;
     vec<vec<vec<vec<int>>>>&           _dlines;
     ReadPathVec&                       _dpaths;
     VecULongVec&                       _dpaths_index;

     // built by Rebuild()
     digraphE<int>                      _LG;
     vec<int>                           _lg_to_left;
     vec<int>                           _lg_to_right;
     vec<int>                           _dlens;
     vec<int>                           _llens;
     vec<int>                           _tol;
     vec<int>                           _to_left;
     vec<int>                           _to_right;
     vec<int>                           _linv;

     // dynamically allocated, if not provided
     map<void*, bool>                   _we_own;  // this must come before any variable that is used with manage() above
     vec<double>                        _cov;
     vec<vec<pair<int,int>>>            _lhood;
     MasterVec<SerfVec<pair<int,int>>>  _lbpx;

     ThreadedLogger*                    _plog;
     bool                               _debug;
     bool                               _we_own_plog;

     // optional Tangle stuff
     equiv_rel*                         _ple;
     vec<int>*                          _plreps;
     int                                _max_short;
     int                                _min_long;
};


#endif // _10X_LINEGRAPHOPS_H
