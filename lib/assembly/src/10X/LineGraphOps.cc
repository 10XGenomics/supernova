// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
//
//
// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS
//

#include "10X/LineGraphOps.h"
#include "10X/Heuristics.h"
#include <omp.h>

size_t LineGraphOps::InitTangles( const int MAX_SHORT /* = 500 */, const int MIN_LONG /* = 4000 */ ) {
     _max_short = MAX_SHORT;
     _min_long = MIN_LONG;
     if ( _ple ) delete _ple;
     _ple = new equiv_rel;
     _ple->Initialize( _dlines.size() );

     for ( int v = 0; v < _LG.N(); v++ ) {
          for ( auto & l : _LG.ITo(v) ) {
               if ( _llens[l] > MAX_SHORT ) continue;
               for ( auto & r : _LG.IFrom(v) ) {
                    if ( _llens[r] > MAX_SHORT ) continue;
                    _ple->Join( l, r );
               }
          }
     }

     if ( _plreps ) delete _plreps;
     _plreps = new vec<int>;
     _ple->OrbitRepsAlt( *_plreps );

     return _plreps->size();
}

bool LineGraphOps::FindTangle( const int L, vec<int>& ins, vec<int>& outs, vec<int>& tangle ) {
     ForceAssertLt( L, _llens.size() );
     if ( !_ple  ) FatalErr( "FindTangle called before InitTangles" );

     int rep = _ple->ClassId(L);

     if ( !BinMember( *_plreps, rep ) ) return false;
     if ( _llens[rep] > _max_short )  {
          return false;
     }

     _ple->Orbit( rep, tangle );
     Sort( tangle );

     for ( auto & t : tangle ) {
          for ( int pass = 1; pass <= 2; pass++ ) {
               const int v = ( pass == 1 ? _lg_to_left[t] : _lg_to_right[t] );
               for ( auto & x : _LG.ITo(v) )
                    if ( !BinMember( tangle, x ) )
                         ins.push_back(x);
               for ( auto & x : _LG.IFrom(v) )
                    if ( !BinMember( tangle, x ) )
                         outs.push_back(x);
          }
     }

     UniqueSort( ins ), UniqueSort( outs );

     Bool LONG = True;
     for ( auto & x : ins )
          if ( _llens[x] < _min_long )
               LONG=False;
     for ( auto & x : outs )
          if ( _llens[x] < _min_long )
               LONG=False; 

     if ( !LONG ) {
          ins.clear();
          outs.clear();
          tangle.clear();
          return false;
     }

     return true;
}


bool LineGraphOps::Reach( int L, vec<int>& reach ) {
     if ( !_ple  ) FatalErr( "Reach called before InitTangles" );
     reach.clear();

     ForceAssertLt(L, _dlines.size() );
     ForceAssertLt(_dlines[L].back().back().back(), _to_right.size() );

     // find right reach of line
     int w = _to_right[_dlines[L].back().back().back()];

     vec<int> vvv({w});
     int INTER_DEPTH=10;

     while ( vvv.size() && INTER_DEPTH-- ) { 
          w=vvv.back(); vvv.pop_back();
          for ( int j = 0; j < _D.From(w).isize(); ++j ) {
               int L2 = _tol[_D.IFrom(w,j)];
               if ( _llens[L2] >= _min_long ) reach.push_back( L2 );
               else if ( _llens[L2] <= _max_short ) {
                    vec<int> ins, outs, tangle;
                    auto good = this->FindTangle( L2, ins, outs, tangle );
                    if ( good ) reach.append(outs);
               } else {
                    w=_to_right[_dlines[L2].back().back().back()];
                    vvv.push_back(w);
               }
          }
     }

     if ( vvv.size() ) {
          reach.clear();
          return false;
     }

     UniqueSort(reach);
     return true;
}


int LineGraphOps::LeftmostRightNeighbor( int L, vec<int>* preach ) {

     NeedLhood();
     NeedLbpx();
     NeedCov();

     vec<pair<int,int>> LH;
     const int MIN_LEN = _min_long;
     const int MIN_LEFT_IGNORE = 100;
     const int MIN_BIG = 12000;
     const double MAX_CN_DIFF = 0.25;
     const int MIN_ADVANTAGE = 100;

     vec<int> tmp_reach;

     if ( !preach ) {
          Reach(L, tmp_reach);
          preach = &tmp_reach;
     }

     vec<int>& reach(*preach);

     if ( _debug ) _plog->get() << Date() << ": LGO reach for " \
          << L << ": " << printSeq(reach) << ", lhood:" << endl;

     for ( int i = 0; i < _lhood[L].isize( ); i++ )
     {    
          int L2 = _lhood[L][i].second;
          if ( L2 == L || L2 == _linv[L] ) continue;
          if ( _llens[L2] >= MIN_LEN ) {
               LH.push_back( _lhood[L][i] );
               if ( _debug ) _plog->get() << "\t(" << _lhood[L][i].first << "," << _lhood[L][i].second << ") " << endl;
          } 
     }

     Bool confused = False;
     vec<int> X;
     vec<Bool> good;
     for ( int i = 0; i < LH.isize( ); i++ ) {
          int L2 = LH[i].second;
          if ( _llens[L2] < MIN_LEN ) continue;
          vec<double> scores;
          vec< triple<int,int,int> > M;
          scores.push_back( ScoreOrder( { L2, L }, _lbpx, _llens, M ) );
          scores.push_back( ScoreOrder( { _linv[L2], L }, _lbpx, _llens, M ) );
          scores.push_back( ScoreOrder( { L, L2 }, _lbpx, _llens, M ) );
          scores.push_back( ScoreOrder( { L, _linv[L2] }, _lbpx, _llens, M ) );
          double left_adv =
               Min( scores[2], scores[3] ) - Min( scores[0], scores[1] );
          if ( reach.nonempty( ) && !BinMember( reach, L2 ) 
               && left_adv < MIN_LEFT_IGNORE && left_adv > -MIN_LEFT_IGNORE ) {    
               continue;    
          }
          if ( left_adv >= MIN_LEFT_IGNORE ) continue;
          if ( left_adv > -MIN_LEFT_IGNORE ) {    
               confused = True;
               break;    
          }
          vec<int> ids( 4, vec<int>::IDENTITY );
          SortSync( scores, ids );
          double win = scores[1] - scores[0];
          X.push_back( ids[0]==2 ? L2 : _linv[L2] );
          good.push_back( win >= MIN_ADVANTAGE );
     }
     if ( confused ) { return -2; }
     if ( X.empty() ) { return -1; }

     if ( X.size( ) > 1 ) {
          for ( int j2 = 0; j2 < X.isize( ); j2++ ) {
               int L2 = X[j2];
               Bool confused = False;
               for ( int j3 = 0; j3 < X.isize( ); j3++ ) {    
                    int L3 = X[j3];
                    if ( L3 == L2 ) continue;
                    vec<double> scores;
                    vec< triple<int,int,int> > M;
                    scores.push_back( ScoreOrder( {L3, L2}, _lbpx, _llens, M ) );
                    scores.push_back( ScoreOrder(
                         { _linv[L3], L2 }, _lbpx, _llens, M ) );
                    scores.push_back( ScoreOrder( {L2, L3}, _lbpx, _llens, M ) );
                    scores.push_back(
                         ScoreOrder( {L2, _linv[L3]}, _lbpx, _llens, M ) );
                    double left_adv = Min(
                         scores[2], scores[3] ) - Min( scores[0], scores[1] );

                    // May not be as in BarcodeJoin:

                    if ( reach.nonempty( ) && !BinMember( reach, L3 ) &&
                         left_adv >= -MIN_LEFT_IGNORE && left_adv <= 0 )
                    {    continue;    }

                    if ( left_adv >= -MIN_LEFT_IGNORE ) {
                         confused = True;
                         break;    
                    }    
               }
               if ( !confused ) {
                    if ( good[j2] ) X = {L2};
                    break;    
               }    
          }   
     }
     if ( X.solo( ) ) {    
          int L2 = X[0];
          if ( _llens[L2] >= MIN_BIG && Abs( _cov[L] - _cov[L2] ) < MAX_CN_DIFF )
               return L2;
     }

     if ( _debug ) _plog->get() << "returning without a solution, X=" << printSeq(X) << endl;

     return -3;
}


vec<LineGraphOps::OOScore> LineGraphOps::ProbeLhood( int L, double S1 /* = 40 */, bool single /* = false */ ) 
{    
     NeedLhood(); 
     NeedLbpx(); 
     NeedCov();

     vec<OOScore> scores;
     ForceAssertLt( L, _dlines.size() );

     // Define constants.

     const double S2 = 1000000000.; // high score to print
     const int MIN_LEN = 2000;  // min length of line in nhood
     const int MIN_MASS = 1000; // min mass to show align

     // Compute neighborhood.
     vec<vec<quad<int,double,int,int>>> nhoods(omp_get_max_threads() ); // ( left/right, -score, L2, b )
     const double MIN_ADVANTAGE = 40.0;
     vec<pair<int,int>> LH;
     for ( int i = 0; i < _lhood[L].isize( ); i++ )
     {    int L2 = _lhood[L][i].second;
          if ( _llens[L2] >= MIN_LEN ) LH.push_back( _lhood[L][i] );    }
     vec< triple<int,int,int> > M;
     #pragma omp parallel for schedule(dynamic,1) private(M) num_threads( single ? 1 : getConfiguredNumThreads() )
     for ( int i = 0; i < LH.isize( ); i++ ) {    
          auto tno = omp_get_thread_num();
          int L2 = LH[i].second;
          vec<double> scores;
          double start=WallClockTime();
          scores.push_back( ScoreOrder( { L2, L }, _lbpx, _llens, M ) );
          scores.push_back( ScoreOrder( { _linv[L2], L }, _lbpx, _llens, M ) );
          scores.push_back( ScoreOrder( { L, L2 }, _lbpx, _llens, M ) );
          scores.push_back( ScoreOrder( { L, _linv[L2] }, _lbpx, _llens, M ) );
          vec<int> ids( 4, vec<int>::IDENTITY );
          SortSync( scores, ids );
          double win = scores[1] - scores[0];
          if ( win < S1 || win > S2 ) continue;
          nhoods[tno].push( ( ids[0]==0||ids[0]==1 ? -1 : +1 ),
                      -win, ( ids[0]==0||ids[0]==2 ? L2 : _linv[L2] ),
                      LH[i].first );
     }
     auto nhood = Flatten( nhoods );
     Sort(nhood);

     // Generate results entries for center line.
     scores.push_back( 
               LineGraphOps::OOScore{L, _llens[L], 0, 
               std::max(0., _cov[L]), 0, 0., 0.});

     // Generate results entries for nhood.

     for ( int i = 0; i < nhood.isize( ); i++ ) {
          int L2 = nhood[i].third;
          scores.push_back(
               LineGraphOps::OOScore{ L2, _llens[L2], nhood[i].fourth,
               std::max(0., _cov[L2]), nhood[i].first, -nhood[i].second, -nhood[i].second } );
     }

    return scores; 
}

bool LineGraphOps::IsInvTangle( vec<int> const& ins, vec<int> const& outs, vec<int> const& tangle ) {
     vec<int> iouts;
     for ( auto const in : ins ) iouts.push_back( _dinv[in] );


     string const& name = "IsInvTangle";
     UniqueSort(iouts);

     bool good = false;
     if ( iouts.size() == outs.size() ) good = std::equal(outs.begin(), outs.end(), iouts.begin() );

     if ( _debug ) {
          Log() << name << " inputs and outputs are " << 
               (good ? "" : "not ") << "inversions of one another: " << endl
               << "ins: " << printSeq(ins) << endl << "outs: " << printSeq(outs) << endl;
     }

     return good;
}


vec<pair<int,int>> LineGraphOps::FindInvPivots() {
     vec<pair<int,int>> results;

     vec<int> is_v1( _LG.N(), -1 );

     #pragma omp parallel for
     for ( int i = 0; i < _LG.N(); ++i ) {
          int L0=-1, L1=-1;
          if ( _LG.FromSize(i) == 1 ) {
               L0 = L1 = _LG.EdgeObjectByIndexFrom(i,0);
          } else if ( _LG.FromSize(i) == 2 ) {
               L0 = _LG.EdgeObjectByIndexFrom(i,0);
               L1 = _LG.EdgeObjectByIndexFrom(i,1);
          } else if ( _LG.FromSize(i) == 3 ) {
               // okay, this is lame, but we'll generalize once we see if it works.
               int tL0 = _LG.EdgeObjectByIndexFrom(i,0);
               int tL1 = _LG.EdgeObjectByIndexFrom(i,1);
               int tL2 = _LG.EdgeObjectByIndexFrom(i,2);
               if ( _LG.From(i)[0] == _LG.From(i)[1] && tL0 == _linv[tL1] ) {
                    L0=tL0; L1=tL1;
               } else if ( _LG.From(i)[0] == _LG.From(i)[2] && tL0 == _linv[tL2] ) {
                    L0=tL0; L1=tL2;
               } else if ( _LG.From(i)[1] == _LG.From(i)[2] && tL1 == _linv[tL2] ) {
                    L0=tL1; L1=tL2;
               }
          }
          if ( L0 != -1 && L1 != -1 ) {
               if ( L0 == _linv[L1] && _lg_to_right[L0] == _lg_to_right[L1] ) {
                    is_v1[i] = _lg_to_right[L1];
                    if ( _debug) {
                         Log() << "InvPivot: " << i << " o " << L0 << 
                              "," << L1 << " o " << is_v1[i] << endl;
                    }
               }
          }
     }

     for ( int i = 0; i < is_v1.isize(); ++i )  {
          if ( is_v1[i] >= 0 ) {
               results.push_back( make_pair( i, is_v1[i] ) );
          }
     }

     return results;
}


void LineGraphOps::ScoreInvAndResolve(vec<pair<int,int>> const& pivots) {
     NeedCov();

     const double MAX_CN_DIFF = 0.25;
     size_t nedits=0;
     for ( auto const& pivot : pivots ) {
          int v = pivot.first;
          int w = pivot.second;

          // first case - peel away
          auto ins = _LG.ToEdgeObj(v);
          auto outs = _LG.FromEdgeObj(w);

          int L = _LG.OGoes(v,w);

          if ( ins.size() == 0 || outs.size() == 0 ) continue;

          // first stupid 2x2
          if ( ins.size() == 2 && outs.size() == 2 ) {
               if ( _debug ) Log() << "assess pivot v=" << v << ", L=" << L << endl;
               int L0,L1;
               if ( ins[0] == _linv[outs[0]] && ins[1] == _linv[outs[1]] ) {
                    L0 = ins[0];
                    L1 = outs[1];
               } else if ( ins[0] == _linv[outs[1]] && ins[1] == _linv[outs[0]] ) {
                    L0 = ins[0];
                    L1 = outs[0];
               } else continue;
               if ( _debug ) Log() << "... L0=" << L0 << ", L1=" << L1 << endl;

               if ( _llens[L0] < _min_long || _llens[L1] < _min_long ) 
                    continue;

               if ( _debug ) Log() << "... long enough " << endl;

               auto s=BarcodePairSupport( make_pair(L0,L1) );
               auto cndiff=Abs(_cov[L0]-_cov[L1]);
               if ( s < 100 || cndiff > MAX_CN_DIFF ) {
                    if ( _debug ) Log() << "... failed scoring L0=" << L0 << ", L1=" << L1 << 
                         ", score=" << s << ", cndiff=" << cndiff << endl;
                    continue;
               }

               MakeRipEdit( L0, L1, { L } );
               nedits++;
          } else if ( ins.size() == 1 && outs.size() == 1 && 
                    ins[0] == _linv[outs[0]] && outs[0] == _linv[ins[0]] ) {
               // above is over determined, but we're double checking
               vec<int> lll;
               for ( int j = 0; j < _LG.FromSize(v); ++j ) {
                    if ( _LG.From(v)[j] == w ) lll.push_back( _LG.OFrom(v, j) );
               }
               if (lll.size() == 1 && _linv[lll[0]] == lll[0] ) {
                    lll.push_back( lll[0] );
               }
               if ( lll.size() == 2 && _linv[lll[0]] == lll[1] ) {
                    int N = _D.N();
                    _D.AddVertices(2);
                    int d1 = _dlines[lll[0]].front().front().front();
                    int d2 = _dlines[lll[0]].back().back().back();
                    _D.GiveEdgeNewFromVxWithUpdate( d1, _to_left[d1], N, _to_left );
                    _D.GiveEdgeNewToVxWithUpdate( d2, _to_right[d2], N+1, _to_right );
                    if ( lll.size() == 2 ) {
                         d1 = _dlines[lll[1]].front().front().front();
                         d2 = _dlines[lll[1]].back().back().back();
                         _D.GiveEdgeNewFromVxWithUpdate( d1, _to_left[d1], N, _to_left );
                         _D.GiveEdgeNewToVxWithUpdate( d2, _to_right[d2], N+1, _to_right );
                    }
                    if ( _debug ) Log() << "removed lines " << printSeq(lll) << endl;
                    nedits++;
               }
          } else {
               if ( _debug ) Log() << "NOT assessing pivot v=" << v << ", L=" << L << endl;
          }
     }
     DebugSync();
     RemoveUnneededVertices( _D, _dinv );
     CleanupCore( _D, _dinv );
     Validate( _hb, _inv, _D, _dinv );    
     cout << Date() << ": made " << nedits << " inversion edits" << endl;
}


int LineGraphOps::BarcodePairSupport( pair<int,int> const& lpair )
{
     NeedLhood();
     int supp1 = 0, supp2 = 0;

     for ( auto const& lh : _lhood[lpair.first] )
          if ( lh.second == lpair.second ) {
               supp1 = lh.first;
               break;
          }

     for ( auto const& lh : _lhood[lpair.second] )
          if ( lh.second == lpair.first ) {
               supp2 = lh.first;
               break;
          }

     if ( _debug ) {
          Log() << "scoring " << lpair.first << "," << lpair.second << endl;
          Log() << "\tlhood " << lpair.first << ": ";
          for ( auto pair : _lhood[lpair.first] ) 
               Log() << "(" << pair.first << "," << pair.second << ") ";
          Log() << endl;
          Log() << "\tlhood " << lpair.second << ": ";
          for ( auto pair : _lhood[lpair.second] ) 
               Log() << "(" << pair.first << "," << pair.second << ") ";
          Log() << "\tsupp1=" << supp1 << ", supp2=" << supp2 << endl;
     }

     return supp1+supp2;
}


void LineGraphOps::MakeRipEdit( int L1, int L2, vec<int> lsr ) {
     if ( _debug ) Log() << "MakeRipEdit L1=" << L1 << 
          ", L2=" << L2 << ", lsr=" << lsr << endl;

     int d1 = _dlines[L1].back().back().back(); 
     int d2 = _dlines[L2].front().front().front();
     PRINT2(d1,d2);

     vec<int> em, emr;
     for ( auto L : lsr ) {    
          const vec<vec<vec<int>>>& X = _dlines[L];
          em.append( Contents(X) );
          for ( int j = 1; j < X.isize( ); j += 2 ) {
               if ( X[j].solo( ) && X[j][0].empty( ) ) {
                    int v = _to_right[ X[j-1][0][0] ];
                    em.push_back( _D.IFrom(v,0) );    
               }    
          }    
     }
     UniqueSort(em); 

     for ( auto d : em ) emr.push_back( _dinv[d] );
     int E = _D.E( ), n = em.size( );
     digraphE<vec<int>> M( digraphE<vec<int>>::COMPLETE_SUBGRAPH_EDGES, _D, em, _to_left, _to_right );
     _D.AppendWithUpdate( M, _to_left, _to_right );
     digraphE<vec<int>> RM( digraphE<vec<int>>::COMPLETE_SUBGRAPH_EDGES, _D, emr, _to_left, _to_right );
     _D.AppendWithUpdate( RM, _to_left, _to_right );

     cout << "em=" << printSeq(em) << endl;
     cout << "emr=" << printSeq(emr) << endl;

     _dinv.resize( E + 2*n );
     for ( int i = 0; i < n; i++ )
     {    _dinv[ E + i ] = E + n + i, _dinv[ E + n + i ] = E + i;    }

     int rd1 = _dinv[d1], rd2 = _dinv[d2];
     if ( n == 0 )
     {    int N = _D.N( );
          _D.AddVertices(2);
          _D.GiveEdgeNewToVxWithUpdate( d1, _to_right[d1], N, _to_right );
          _D.GiveEdgeNewFromVxWithUpdate( d2, _to_left[d2], N, _to_left );
          _D.GiveEdgeNewToVxWithUpdate( rd2, _to_right[rd2], N+1, _to_right );
          _D.GiveEdgeNewFromVxWithUpdate( rd1, _to_left[rd1], N+1, _to_left );    }
     else
     {    int v1 = -1, v2 = -1, rv1 = -1, rv2 = -1;
          for ( int i = 0; i < em.isize( ); i++ )
          {    int d = em[i];
               if ( _to_left[d] == _to_right[d1] )
               {    v1 = _to_left[E+i];
                    break;    }    }
          for ( int i = 0; i < em.isize( ); i++ )
          {    int d = em[i];
               if ( _to_right[d] == _to_left[d2] )
               {    v2 = _to_right[E+i];
                    break;    }    }
          for ( int i = 0; i < emr.isize( ); i++ )
          {    int d = emr[i];
               if ( _to_left[d] == _to_right[rd2] )
               {    rv2 = _to_left[E+n+i];
                    break;    }    }
          for ( int i = 0; i < emr.isize( ); i++ )
          {    int d = emr[i];
               if ( _to_right[d] == _to_left[rd1] )
               {    rv1 = _to_right[E+n+i];
                    break;    }    }
          PRINT4(v1,v2,rv1,rv2);
          PRINT2(rd1,rd2);
          PRINT(_to_right[d1]);
          PRINT(_to_left[d2]);
          PRINT(_to_right[rd2]);
          PRINT(_to_left[rd1]);
          _D.GiveEdgeNewToVxWithUpdate( d1, _to_right[d1], v1, _to_right );
          _D.GiveEdgeNewFromVxWithUpdate( d2, _to_left[d2], v2, _to_left );
          _D.GiveEdgeNewToVxWithUpdate( rd2, _to_right[rd2], rv2, _to_right );
          _D.GiveEdgeNewFromVxWithUpdate( rd1, _to_left[rd1], rv1, _to_left );    
     }    
}

