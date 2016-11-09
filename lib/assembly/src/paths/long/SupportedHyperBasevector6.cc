///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * SupportedHyperBasevector6.cc
 *
 *  Created on: Apr 3, 2013
 *      Author: tsharpe
 */

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "CoreTools.h"
#include "IteratorRange.h"
#include "ParallelVecUtilities.h"
#include "Qualvector.h"
#include "paths/long/CreateGenome.h"
#include "paths/long/EvalByReads.h"
#include "paths/long/KmerCount.h"
#include "paths/long/SupportedHyperBasevector.h"
#include "random/Bernoulli.h"
#include "system/WorklistN.h"
#include <algorithm>
#include <set>

namespace { // open anonymous namespace

class BubbleProc
{
public:
    BubbleProc( SupportedHyperBasevector const& shbv,
                SupportedHyperBasevector::BubbleParams const& parms,
                SupportedHyperBasevector::BubbleAux const& aux,
                long_logging_control const& log_control,
                long_logging const& logc,
                vec< std::pair<int,vec<int> > >* pSubs,
                const long_heuristics& heur )
    : mSHBV(shbv), mParms(parms), mAux(aux), mLogControl(log_control),
      mLogC(logc), mpSubs(pSubs), mheur(&heur) {}

    void operator()( size_t edgeId )
    { vec<int> substPath;
      if ( mSHBV.IsPoppable(edgeId,mParms,mAux,mLogControl,mLogC,&substPath,*mheur) )
      { SpinLocker lock(gLockedData); mpSubs->push(edgeId,substPath); } }

private:
    SupportedHyperBasevector const& mSHBV;
    SupportedHyperBasevector::BubbleParams const& mParms;
    SupportedHyperBasevector::BubbleAux const& mAux;
    long_logging_control const& mLogControl;
    long_logging const& mLogC;
    vec< std::pair<int,vec<int> > >* mpSubs;
    static SpinLockedData gLockedData;
    const long_heuristics* mheur;
};

SpinLockedData BubbleProc::gLockedData;

} // close anonymous namespace

void TruncateMe( vec<int>& p, Bool& del, int& start, int& stop, 
     const vec<Bool>& used, const SupportedHyperBasevector& shb )
{    
     del = False;
     vec< vec<int> > subs, subs_origin;
     vec<int> s, s_origin;
     for ( int j = 0; j <= p.isize( ); j++ )
     {    if ( j == p.isize( ) || ( p[j] >= 0 && !used[ p[j] ] ) )
          {    if ( s.nonempty( ) ) 
               {    if ( s.back( ) < 0 ) 
                    {    s.pop_back( );
                         s_origin.pop_back( );    }
                    subs.push_back(s);
                    subs_origin.push_back(s_origin);    }
               s.clear( );
               s_origin.clear( );    }
          else 
          {    s.push_back( p[j] );
               s_origin.push_back(j);    }    }
     vec<int> nkmers( subs.size( ), 0 ), ids( subs.size( ), vec<int>::IDENTITY );
     for ( int j = 0; j < subs.isize( ); j++ )
     {    for ( int l = 0; l < subs[j].isize( ); l++ )
               nkmers[j] += shb.EdgeLengthKmers( subs[j][l] );    }
     ReverseSortSync( nkmers, ids );
     if ( nkmers.empty( ) || ( nkmers.size( ) >= 2 && nkmers[0] == nkmers[1] ) )
          del = True;
     else
     {    p = subs[ ids[0] ];
          start = subs_origin[ ids[0] ].front( );
          stop = subs_origin[ ids[0] ].back( );    }    }

void SupportedHyperBasevector::TruncatePaths( const long_logging& logc )
{    double clock = WallClockTime( );
     vec<Bool> used, to_delete( NPaths( ), False );
     Used(used);
     for ( int i = 0; i < NPaths( ); i++ )
     {    Bool del;
          int start, stop;
          TruncateMe( PathMutable(i), del, start, stop, used, *this );
          if (del) to_delete[i] = True;    }
     EraseIf( PathsMutable( ), to_delete );
     EraseIf( WeightsFwMutable( ), to_delete );
     EraseIf( WeightsRcMutable( ), to_delete );
     // Temporary if.
     if ( WeightsFwOrigin( ).size( ) == WeightsFw( ).size( ) )
     {    EraseIf( WeightsFwOriginMutable( ), to_delete );
          EraseIf( WeightsRcOriginMutable( ), to_delete );    }
     to_delete.resize_and_set( NPairs( ), False );
     for ( int i = 0; i < NPairs( ); i++ )
     {    Bool del1, del2;
          int start1, stop1, start2, stop2;
          vec<int> left = PairLeftMutable(i), right = PairRightMutable(i);
          TruncateMe( PairLeftMutable(i), del1, start1, stop1, used, *this );
          TruncateMe( PairRightMutable(i), del2, start2, stop2, used, *this );
          if ( del1 || del2 ) to_delete[i] = True;
          else
          {    int sub = 0;
               for ( int k = stop1; k < left.isize( ); k++ )
                    sub += EdgeLengthKmers( left[k] );
               for ( int k = 0; k < start2; k++ )
                    sub += EdgeLengthKmers( right[k] );
               AddTrim( i, -sub );    }    }
     EraseIf( PairsMutable( ), to_delete );
     EraseIf( PairDataMutable( ), to_delete );
     UniqueOrderPaths( );
     REPORT_TIME( clock, "used in TruncatePaths" );    }

void SupportedHyperBasevector::PopBubbles( BubbleParams const& parms,
        unsigned nThreads, const long_logging_control& log_control, 
        const long_logging& logc, const long_heuristics& heur )
{
     double clock1 = WallClockTime( );
     BubbleAux aux;
     ToLeft(aux.to_left), ToRight(aux.to_right);
     aux.paths_index.resize( EdgeObjectCount() );
     aux.mult_fw.resize( EdgeObjectCount( ), 0 );
     aux.mult_rc.resize( EdgeObjectCount( ), 0 );
     for ( int i = 0; i < Paths( ).isize( ); i++ )
     for ( int j = 0; j < Path(i).isize( ); j++ )
     {    if ( Path(i,j) >= 0 )
          {    aux.mult_fw[ Path(i,j) ] += WeightFw(i);
               aux.mult_rc[ Path(i,j) ] += WeightRc(i);
               aux.paths_index[ Path(i,j) ].push( i, j );    }    }
     vec< pair< int, vec<int> > > subs;
     if ( true )
     {
         BubbleProc proc(*this,parms,aux,log_control,logc,&subs,heur);
         parallelFor(0,EdgeObjectCount(),proc,nThreads);
     }
     REPORT_TIME( clock1, "used popping 1" );
     double clock2 = WallClockTime( );
     vec< vec< pair<int,int> > > subs_index( EdgeObjectCount( ) );
     for ( int i = 0; i < subs.isize( ); i++ )
     for ( int j = 0; j < subs[i].second.isize( ); j++ )
          subs_index[ subs[i].second[j] ].push( i, j );
     for ( int i = 0; i < subs.isize( ); i++ )
     {    int e = subs[i].first;
          vec<int> p = subs[i].second;
          for ( int j = 0; j < aux.paths_index[e].isize( ); j++ )
          {    int id = aux.paths_index[e][j].first;
               int pos = aux.paths_index[e][j].second;
               vec<int>& q = PathsMutable( )[id];
               for ( int l = pos; l < q.isize( ); l++ )
               {    int r = BinPosition( aux.paths_index[ q[l] ], make_pair(id,l) );
                    if ( r >= 0 ) 
                    {    aux.paths_index[ q[l] ].erase(
                              aux.paths_index[ q[l] ].begin( ) + r );    }    }
               vec<int> qnew(q);
               qnew.resize(pos);
               qnew.append(p);
               for ( int l = pos+1; l < q.isize( ); l++ )
                    qnew.push_back( q[l] );
               q = qnew;
               for ( int l = pos; l < q.isize( ); l++ )
               {    aux.paths_index[ q[l] ].push_back( make_pair(id,l) );
                    Sort( aux.paths_index[ q[l] ] );    }
               aux.paths_index[e].clear( );    }
          for ( int j = 0; j < subs_index[e].isize( ); j++ )
          {    int id = subs_index[e][j].first, pos = subs_index[e][j].second;
               vec<int>& q = subs[id].second;
               for ( int l = pos; l < q.isize( ); l++ )
               {    int r = BinPosition( subs_index[ q[l] ], make_pair(id,l) );
                    if ( r >= 0 ) 
                    {    subs_index[ q[l] ].erase( 
                              subs_index[ q[l] ].begin( ) + r );    }    }
               vec<int> qnew(q);
               qnew.resize(pos);
               qnew.append(p);
               for ( int l = pos+1; l < q.isize( ); l++ )
                    qnew.push_back( q[l] );
               q = qnew;
               for ( int l = pos; l < q.isize( ); l++ )
               {    subs_index[ q[l] ].push_back( make_pair(id,l) );
                    Sort( subs_index[ q[l] ] );    }
               subs_index[e].clear( );    }    }

     vec<int> to_delete;
     to_delete.reserve(subs.size());
     for ( size_t idx = 0; idx != subs.size(); ++idx )
         to_delete.push_back(subs[idx].first);
     DeleteEdges(to_delete);
     RemoveEdgelessVertices( );
     RemoveUnneededVertices( );
     REPORT_TIME( clock2, "used popping 2" );
     RemoveDeadEdgeObjects( );
     double clock3 = WallClockTime( );
     UniqueOrderPaths( );

     // Force symmetry of paths and weights.  I suppose we might instead try to 
     // figure out how they got asymmetric.

     vec< vec<int> > new_paths;
     vec<fix64_6> new_weights_fw, new_weights_rc;
     vec< vec< pair<fix64_6,int64_t> > > new_weights_fw_origin;
     vec< vec< pair<fix64_6,int64_t> > > new_weights_rc_origin;
     for ( int i1 = 0; i1 < NPaths( ); i1++ )
     {    const vec<int>& p1 = Path(i1);
          if ( !InvDef( p1[0] ) ) continue;
          vec<int> p2;
          for ( int j = 0; j < p1.isize( ); j++ )
          {    if ( p1[j] < 0 ) p2.push_back( p1[j] );
               else p2.push_back( Inv( p1[j] ) );    }
          p2.ReverseMe( );
          int i2 = BinPosition( Paths( ), p2 );
          if ( i2 < 0 )
          {    new_paths.push_back(p2);
               new_weights_fw.push_back( WeightFw(i1) );
               new_weights_rc.push_back( WeightRc(i1) );    
               // Temporary if.
               if ( WeightsFwOrigin( ).size( ) == WeightsFw( ).size( ) )
               {    new_weights_fw_origin.push_back( WeightFwOrigin(i1) );
                    new_weights_rc_origin.push_back( WeightRcOrigin(i1) );    }
                    }    }
     PathsMutable( ).append(new_paths);
     WeightsFwMutable( ).append(new_weights_fw);
     WeightsRcMutable( ).append(new_weights_rc);
     WeightsFwOriginMutable( ).append(new_weights_fw_origin);
     WeightsRcOriginMutable( ).append(new_weights_rc_origin);
     // Temporary if.
     if ( WeightsFwOrigin( ).size( ) == WeightsFw( ).size( ) )
     {    SortSync( PathsMutable( ), WeightsFwMutable( ), WeightsRcMutable( ),
               WeightsFwOriginMutable( ), WeightsRcOriginMutable( ) );    }
     else SortSync( PathsMutable( ), WeightsFwMutable( ), WeightsRcMutable( ) );

     // Maybe we're doing something here that is fundamentally wrong.
     // (NOW TURNED OFF.)

     int npairs = NPairs( );
     /*
     vec<Bool> todel( npairs, False );
     for ( int i1 = 0; i1 < NPairs( ); i1++ )
     {    const vec<int> &p1 = PairLeft(i1), &q1 = PairRight(i1);
          if ( InvDef( p1[0] ) != InvDef( q1[0] ) ) todel[i1] = True;    }
     EraseIf( PairsMutable( ), todel );
     EraseIf( PairDataMutable( ), todel );
     */

     for ( int i1 = 0; i1 < npairs; i1++ )
     {    const vec<int> &p1 = PairLeft(i1), &q1 = PairRight(i1);
          if ( !InvDef( p1[0] ) || !InvDef( q1[0] ) ) continue;
          vec<int> p2, q2;
          for ( int j = 0; j < p1.isize( ); j++ )
          {    if ( p1[j] < 0 ) p2.push_back( p1[j] );
               else p2.push_back( Inv( p1[j] ) );    }
          for ( int j = 0; j < q1.isize( ); j++ )
          {    if ( q1[j] < 0 ) q2.push_back( q1[j] );
               else q2.push_back( Inv( q1[j] ) );    }
          p2.ReverseMe( ), q2.ReverseMe( );
          int i2 = BinPosition( Pairs( ), make_pair( q2, p2 ) );
          if ( i2 < 0 )
          {    PairsMutable( ).push_back( make_pair( q2, p2 ) );
               PairDataMutable( ).push_back( PairData(i1) );    
               SortSync( PairsMutable( ), PairDataMutable( ) );    }    }

     REPORT_TIME( clock3, "used popping 3" );
     FixWeights(logc);
     TestValid(logc);    }

bool SupportedHyperBasevector::IsPoppable( int e, BubbleParams const& parms,
                                        BubbleAux const& aux,
                                        long_logging_control const& log_control,
                                        long_logging const& logc,
                                        vec<int>* pSubstPath,
     const long_heuristics& heur ) const
{
     int v = aux.to_left[e];
     int w = aux.to_right[e];

     // don't consider certain edges with strong evidence
     if ( aux.mult_fw[e] > parms.max_pop_del &&
             aux.mult_rc[e] > parms.max_pop_del ) return false;
     if ( Max( aux.mult_fw[e], aux.mult_rc[e] ) > parms.max_pop_del2 ) return false;

     // establish a min and max path length for alternate paths
     int ne = EdgeObject(e).isize( ) - K( ) + 1;
     int low = ne - int( floor( parms.min_pop_ratio * parms.delta_kmers ) );
     int high = ne + int( floor( parms.min_pop_ratio * parms.delta_kmers ) );

     vec< vec<int> > paths, paths_final;
     vec<int> pathlens, pathlens_final;
     vec<fix64_6> mult_fw_final, mult_rc_final;

     // beginning of depth-first recursion from the starting vertex
     // to look for alternate paths
     // NOTE: it would probably be a little faster to go breadth-first
     // instead.  you might more quickly discover multiple alternatives
     // without fully exploring tiny cycles.
     for ( int l = 0; l < From(v).isize( ); l++ )
     {    int f = EdgeObjectIndexByIndexFrom( v, l );
          if ( f != e )
          { paths.push_back(vec<int>(1,f));
            pathlens.push_back( EdgeObject(f).isize( ) - K( ) + 1 );   } }

     // while there are alternatives to consider
     while( paths.nonempty( ) )
     {    vec<int> p = paths.back( );
          int n = pathlens.back( );
          paths.pop_back( );
          // cout << "found path " << printSeq(p) << endl; // XXXXXXXXXXXXXXXXX
          pathlens.pop_back( );
          int x = aux.to_right[ p.back( ) ];

          // if the path under consideration has the right length, and
          // ends up at the same vertex as the subject edge
          if ( n >= low && n <= high && x == w )
          {    fix64_6 pmult_fw = 0, pmult_rc = 0;
               for ( int l = 0; l < aux.paths_index[ p[0] ].isize( ); l++ )
               {    int id = aux.paths_index[ p[0] ][l].first;
                    int pos = aux.paths_index[ p[0] ][l].second;
                    if ( Path(id).Contains( p, pos ) )
                    {    pmult_fw += WeightFw(id);
                         pmult_rc += WeightRc(id);    }    }
               // if it has appropriate weight, it's a possibility
               if ( pmult_fw <= parms.max_pop_del || pmult_rc <= parms.max_pop_del )
                   continue; // if weight insufficient, no need to extend further
               paths_final.push_back(p);
               // cout << "accepting path " << printSeq(p) << endl; // XXX
               pathlens_final.push_back(n);
               mult_fw_final.push_back(pmult_fw);
               mult_rc_final.push_back(pmult_rc);
               // we're only looking for situations where there's a
               // single valid alternative path, so quit early when
               // we've discovered two alternatives
               if ( paths_final.size() > 1 ) break;    }

          // if the path is already too long, don't extend it further
          if ( n > high ) continue;

          // extend path and continue recursive exploration
          for ( int l = 0; l < From(x).isize( ); l++ )
          {    int f = EdgeObjectIndexByIndexFrom( x, l );
               if ( f != e )
               { paths.push_back(p);
                 paths.back().push_back(f);
                 pathlens.push_back(
                    n + EdgeObject(f).isize( ) - K( ) + 1 );    }   }   }

     if ( !paths_final.solo( ) ) return false;
     if ( Abs( pathlens_final[0] - ne ) > parms.delta_kmers ) return false;
     /*
     if ( ( mult_fw_final[0] < parms.min_pop_ratio * aux.mult_fw[e] )
          && ( mult_rc_final[0] < parms.min_pop_ratio * aux.mult_rc[e] ) )
     {    continue;    }
     */

     // Proper statistical test - note missing from other documentation.
     // Compare similar code in ImproveLongHyper.cc.

     const double max_asym_rarity = 0.00001;
     double f1 = mult_fw_final[0].ToDouble(), r1 = mult_rc_final[0].ToDouble();;
     double f2 = aux.mult_fw[e].ToDouble(), r2 = aux.mult_rc[e].ToDouble();
     if ( f2 > r2 || ( f2 == r2 && f1 > r1 ) )
     {    swap( f1, r1 );
          swap( f2, r2 );    }
     int n = int(floor(f1+r1+f2+r2));
     if ( f1 == 0 || n == 0 || n > 10000 ) return false;

     long double p;
     double q;
     if ( !heur.FIX_ASYMMETRY_BUG )
     {    p = Min( 0.5, f1/(f1+r1) ) / 2;
          q = BinomialSum( n, int(ceil(f2)), p );    }
     else
     {    p = 0.5;
          q = BinomialSum( n, int(ceil( Min(f1+r1,f2+r2) )), p );    }
     if ( q >= max_asym_rarity ) return false;

     // Report result.

     if ( logc.verb[ "POP_BUBBLES" ] >= 1 )
     {    static SpinLockedData gLockedData;
          SpinLocker lock(gLockedData);
          cout << "\nPopBubbles: replace e = " << e << " by "
               << printSeq( paths_final[0] ) << endl;
          PRINT2( aux.mult_fw[e], aux.mult_rc[e] );
          PRINT2( mult_fw_final[0], mult_rc_final[0] );    }
     *pSubstPath = paths_final[0];
     return true; }

namespace
{
typedef int EdgeId;
typedef int VertexId;
typedef unsigned SegId;

// An ordered sequence of EdgeIds that describe, e.g., a traversal of a graph
typedef vec<EdgeId> SHBVPath;

// A collection of edgeIds describing, e.g., the set of edges departing a vertex
typedef vec<EdgeId> EdgeCollection;

// A pair of edgeIds
typedef std::pair<EdgeId,EdgeId> EdgePair;

// An element of some read's path on the graph
struct ReadSegment
{
    ReadSegment()=default;
    ReadSegment( EdgeId readId, SegId segId )
    : mReadId(readId), mSegId(segId) {}

    EdgeId mReadId;
    SegId mSegId;
};

// an operation we think is a valid simplification of the graph
typedef triple<EdgeCollection,EdgeCollection,EdgeCollection> Join;

// instructions about what pull-aparts to do
struct Recommendations : private SpinLockedData
{
    void addRecommendation( std::set<SHBVPath> const& bridges11,
                            std::set<SHBVPath> const& bridges22,
                            EdgeCollection const& reach );

    vec<Join> mJoins;
    vec<vec<SHBVPath>> mPaths1, mPaths2;
};

// adding a recommendation is thread safe
void Recommendations::addRecommendation( std::set<SHBVPath> const& bridges11,
                                         std::set<SHBVPath> const& bridges22,
                                         EdgeCollection const& reach )
{
    ForceAssert(bridges11.size());
    ForceAssert(bridges22.size());
    EdgeId e1 = bridges11.begin()->front(), e2 = bridges22.begin()->front();
    EdgeId f1 = bridges11.begin()->back(), f2 = bridges22.begin()->back();
    EdgeCollection eee{e1,e2};
    EdgeCollection fff{f1,f2};
    Join join(eee,reach,fff);
    SpinLocker lock(*this);
    mJoins.push_back(join);
    mPaths1.push(bridges11.begin(),bridges11.end());
    mPaths2.push(bridges22.begin(),bridges22.end());
}

// paths and cum weights for all paths from some edge e to edges f1 or f2
struct PathWeights
{
    PathWeights( SupportedHyperBasevector const& hbv,
                    vec<ReadSegment> const& segs,
                    EdgeId e, EdgeId f1, EdgeId f2 );

    // get the total length of a portion of a path in kmers
    static int getSpan( SupportedHyperBasevector const& hbv,
                        SHBVPath::const_iterator itr,
                        SHBVPath::const_iterator end )
    { int span = 0;
      while ( itr != end ) { span += hbv.EdgeLengthKmers(*itr); ++itr; }
      return span; }

    fix64_6 weight1, weight2; // cum weights of paths to f1 and f2
    std::set<SHBVPath> bridges1, bridges2; // distinct paths to f1 and f2
    int mSpan; // max distance in kmers on any path from the edge to edge f1 or f2
};

PathWeights::PathWeights( SupportedHyperBasevector const& hbv,
                            vec<ReadSegment> const& segs,
                            EdgeId e, EdgeId f1, EdgeId f2 )
  : weight1(0), weight2(0), mSpan(0)
{
    auto segsEnd = segs.end();
    for ( auto segsItr=segs.begin(); segsItr != segsEnd; ++segsItr )
    {
        SHBVPath const& path = hbv.Path(segsItr->mReadId);
        fix64_6 const& weight = hbv.Weight(segsItr->mReadId);
        auto from = path.begin() + segsItr->mSegId;
        auto beg = from + 1;
        auto end = path.end();
        auto itr = std::find(beg,end,f1);
        if ( itr != end )
        {
            if ( std::find(beg,itr,e) != itr ) // skip it if there's a later
                continue;                      // instance of e
            weight1 += weight;
            bridges1.insert(SHBVPath(from,itr+1));
        }
        else if ( (itr=std::find(beg,end,f2)) != end )
        {
            if ( std::find(beg,itr,e) != itr ) // skip if there's a later
                continue;                      // instance of e
            weight2 += weight;
            bridges2.insert(SHBVPath(from,itr+1));
        }
        else
            continue;
        mSpan = std::max(mSpan,getSpan(hbv,beg,itr));
    }
}

// a thing that examines a vertex to decide whether it can be pulled apart
class JoinDiscoverer
{
public:
    JoinDiscoverer( SupportedHyperBasevector const& hbv,
                        vec<int> const& toLeft, vec<int> const& toRight,
                        int verbosity, double minWeightSplit )
    : mHBV(hbv), mToLeft(toLeft), mToRight(toRight),
      mPathsIndex(hbv.EdgeObjectCount()),
      mReadLenFudge(hbv.MedianCorrectedReadLength()),
      mVerbosity(verbosity),
      mMinWeightSplit(minWeightSplit), mMinWeightSplitLow(2)
    { EdgeId nEdges = hbv.Paths().isize();
      for ( EdgeId edgeId = 0; edgeId != nEdges; edgeId++ )
      { vec<EdgeId> const& path = hbv.Path(edgeId);
        SegId nSegs = path.size();
        for ( SegId segId = 0; segId != nSegs; segId++ )
          mPathsIndex[path[segId]].push_back(ReadSegment(edgeId,segId)); } }

    void processVertex( int vertexId, Recommendations* ) const;

private:
    static size_t const MAX_AFFIX_LEN = 5;

    // scan read paths to find unique prefixes of paths that enter a vertex
    vec<SHBVPath> findPrefixes( EdgeCollection const& leftAlts ) const
    {
        vec<SHBVPath> prefixes;
        auto altsEnd = leftAlts.end();
        for ( auto altsItr=leftAlts.begin(); altsItr != altsEnd; ++altsItr )
        {
            vec<ReadSegment> const& segs = mPathsIndex[*altsItr];
            auto segsEnd = segs.end();
            for ( auto segsItr=segs.begin(); segsItr != segsEnd; ++segsItr )
            {
                SHBVPath const& path = mHBV.Path(segsItr->mReadId);
                size_t nextEle = segsItr->mSegId + 1;
                auto pathEnd = path.begin()+nextEle;
                auto pathBeg = pathEnd - std::min(nextEle,MAX_AFFIX_LEN);
                bool found = false;
                auto end = prefixes.end();
                for ( auto itr=prefixes.begin(); itr != end; ++itr )
                {
                    SHBVPath const& prefix = *itr;
                    size_t len = std::min(nextEle,prefix.size());
                    if ( std::equal(pathEnd-len,pathEnd,prefix.end()-len) )
                    {
                        if ( len < MAX_AFFIX_LEN && nextEle > len )
                            *itr = SHBVPath(pathBeg,pathEnd);
                        found = true;
                        break;
                    }
                }
                if ( !found )
                    prefixes.push_back(SHBVPath(pathBeg,pathEnd));
            }
        }
        return prefixes;
    }

    // scan read paths to find unique suffixes of paths that leave a vertex
    vec<SHBVPath> findSuffixes( EdgeCollection const& rightAlts ) const
    {
        vec<SHBVPath> suffixes;
        auto altsEnd = rightAlts.end();
        for ( auto altsItr=rightAlts.begin(); altsItr != altsEnd; ++altsItr )
        {
            vec<ReadSegment> const& segs = mPathsIndex[*altsItr];
            auto segsEnd = segs.end();
            for ( auto segsItr=segs.begin(); segsItr != segsEnd; ++segsItr )
            {
                SHBVPath const& path = mHBV.Path(segsItr->mReadId);
                auto pathBeg = path.begin()+segsItr->mSegId;
                size_t pathLen = std::min(MAX_AFFIX_LEN,
                                                path.size()-segsItr->mSegId);
                bool found = false;
                auto end = suffixes.end();
                for ( auto itr=suffixes.begin(); itr != end; ++itr )
                {
                    SHBVPath const& suffix = *itr;
                    size_t len = std::min(pathLen,suffix.size());
                    if ( std::equal(pathBeg,pathBeg+len,suffix.begin()) )
                    {
                        if ( pathLen > len )
                            *itr = SHBVPath(pathBeg,pathBeg+pathLen);
                        found = true;
                        break;
                    }
                }
                if ( !found )
                    suffixes.push_back(SHBVPath(pathBeg,pathBeg+pathLen));
            }
        }
        return suffixes;
    }

    // proper edge pairs tile the set of path prefixes or suffixes, i.e.,
    // every path contains one or the other edge, and no path contains both
    static vec<EdgePair> getProperEdgePairs( vec<SHBVPath> const& paths )
    {
        std::map<EdgeId,BitVec> map;
        size_t nPaths = paths.size();
        for ( size_t idx = 0; idx != nPaths; ++idx )
        {
            SHBVPath const& path = paths[idx];
            for ( auto itr=path.begin(),end=path.end(); itr != end; ++itr )
            {
                BitVec& bits = map[*itr];
                bits.resize(nPaths,false);
                bits.set(idx,true);
            }
        }
        vec<EdgePair> pairs;
        for ( auto itr=map.begin(),end=map.end(); itr != end; ++itr )
        {
            auto itr2 = itr;
            while ( ++itr2 != end )
            {
                if ( (itr->second & itr2->second).isUniformlyFalse() &&
                        nor(itr->second,itr2->second).isUniformlyFalse() )
                    pairs.push_back(EdgePair(itr->first,itr2->first));
            }
        }
        return pairs;
    }

    // gets a sorted list of EdgeIds that lie between e1 or e2 and f1 or f2.
    void getReach( PathWeights const& pw1, PathWeights const& pw2,
                    EdgeCollection* pReach ) const
    {
        pReach->clear();
        for ( SHBVPath const& path : pw1.bridges1 )
        {
            ForceAssertGe(path.size(),2u);
            pReach->insert(pReach->end(),path.begin()+1,path.end()-1);
        }
        for ( SHBVPath const& path : pw1.bridges2 )
        {
            ForceAssertGe(path.size(),2u);
            pReach->insert(pReach->end(),path.begin()+1,path.end()-1);
        }
        for ( SHBVPath const& path : pw2.bridges1 )
        {
            ForceAssertGe(path.size(),2u);
            pReach->insert(pReach->end(),path.begin()+1,path.end()-1);
        }
        for ( SHBVPath const& path : pw2.bridges1 )
        {
            ForceAssertGe(path.size(),2u);
            pReach->insert(pReach->end(),path.begin()+1,path.end()-1);
        }

        UniqueSort(*pReach);
    }

    // returns false if the reach is not a complete sub-graph delimited by
    // e1, e2, f1, and f2 (i.e., if non-reach, non-boundary edges enter or exit
    // from reach edges).
    bool checkReach( EdgeCollection const& reach,
                        EdgeId e1, EdgeId e2, EdgeId f1, EdgeId f2 ) const
    {
        vec<int> reachVertices;
        reachVertices.reserve(reach.size()+2);
        for ( EdgeId e : reach )
        {
            reachVertices.push_back(mToLeft[e]);
        }
        reachVertices.push_back(mToLeft[f1]);
        reachVertices.push_back(mToLeft[f2]);
        UniqueSort(reachVertices);

        auto beg = reach.begin(), end = reach.end();
        for ( VertexId v : reachVertices )
        {
            for ( EdgeId e : mHBV.ToEdgeObj(v) )
                if ( e != e1 && e != e2 && !std::binary_search(beg,end,e) )
                    return false;
            for ( EdgeId e : mHBV.FromEdgeObj(v) )
                if ( e != f1 && e != f2  && !std::binary_search(beg,end,e) )
                    return false;
        }
        return true;
    }

    // do the weights allow a pull apart?
    bool pullApart( fix64_6 const& weight_11, fix64_6 const& weight_12,
                    fix64_6 const& weight_21, fix64_6 const& weight_22,
                    int span ) const
    { bool result = false;
      if ( weight_11 >= mMinWeightSplitLow && weight_22 >= mMinWeightSplitLow &&
            weight_12 == 0 && weight_21 == 0 )
        result = true;
      if ( weight_11 >= mMinWeightSplit && weight_22 >= mMinWeightSplit &&
            weight_12 + weight_21 <= 2 && span <= mReadLenFudge )
        result = true;
      if ( weight_11 >= mMinWeightSplit/2 && weight_22 >= mMinWeightSplit/2 &&
            weight_12 + weight_21 < 2 && span <= mReadLenFudge )
        result = true;
      return result; }

    // thread-safe logging of a join
    void log( EdgeId e1, EdgeId e2, EdgeId f1, EdgeId f2,
              EdgeCollection const& reach,
              std::set<SHBVPath> const& bridge11,
              std::set<SHBVPath> const& bridge22,
              fix64_6 const& weight_11, fix64_6 const& weight_12,
              fix64_6 const& weight_21, fix64_6 const& weight_22,
              int span, bool result ) const
    { SpinLocker lock(gLockCOUT);
      std::cout << "e1=" << e1
                << " e2=" << e2
                << " f1=" << f1
                << " f2=" << f2
                << "\nreach: " << rangePrinter(reach.begin(),reach.end(),",")
                << "\nbridge11:";
      int idx = 0;
      for ( auto itr=bridge11.begin(),end=bridge11.end(); itr != end; ++itr )
          std::cout << "\n[" << ++idx << "] " << rangePrinter(itr->begin(),itr->end(),",");
      std::cout << "\nbridge22:";
      idx = 0;
      for ( auto itr=bridge22.begin(),end=bridge22.end(); itr != end; ++itr )
          std::cout << "\n[" << ++idx << "] " << rangePrinter(itr->begin(),itr->end(),",");
      std::cout << "\nw11=" << weight_11
                << " w22=" << weight_22
                << " w12=" << weight_12
                << " w21=" << weight_21
                << " span=" << span << '\n';
      if ( !result ) std::cout << "rejected\n";
      std::cout << std::endl; }

    SupportedHyperBasevector const& mHBV;
    vec<int> const& mToLeft;
    vec<int> const& mToRight;
    vec<vec<ReadSegment>> mPathsIndex;
    int mReadLenFudge;
    int mVerbosity;
    double mMinWeightSplit;
    fix64_6 mMinWeightSplitLow;
    static SpinLockedData gLockCOUT; // logging lock
};

SpinLockedData JoinDiscoverer::gLockCOUT;
size_t const JoinDiscoverer::MAX_AFFIX_LEN;

// can this vertex be pulled apart?
void JoinDiscoverer::processVertex( int vertexId, Recommendations* pRecs ) const
{
    EdgeCollection const& leftAlts = mHBV.ToEdgeObj(vertexId);
    size_t nLeftAlts = leftAlts.size();
    // only process nodes with some left-diversity to cut down on re-processing
    if ( nLeftAlts <= 1 )
        return;

    EdgeCollection const& rightAlts = mHBV.FromEdgeObj(vertexId);
    size_t nRightAlts = rightAlts.size();
    if ( !nRightAlts ) // if its a dead-end, skip it
        return;

    // get all distinct paths that lead to this vertex
    vec<SHBVPath> prefixes = findPrefixes(leftAlts);
    size_t nPrefixes = prefixes.size();
    if ( nPrefixes <= 1 )
        return;

    // get all distinct paths that lead from this vertex
    vec<SHBVPath> suffixes = findSuffixes(rightAlts);
    size_t nSuffixes = suffixes.size();
    if ( nSuffixes <= 1 )
        return;

    // get candidate pairs of "e" edges upstream of vertex
    vec<EdgePair> ePairs = getProperEdgePairs(prefixes);
    if ( ePairs.empty() )
        return;

    // get candidate pairs of "f" edges downstream of vertex
    vec<EdgePair> fPairs = getProperEdgePairs(suffixes);
    if ( fPairs.empty() )
        return;

    // for each of the valid pairs of e's and f's, find all the reads that
    // go from one of the e's to one of the f's, and accumulate the read weights
    EdgeCollection reach;
    for ( auto eItr=ePairs.begin(),eEnd=ePairs.end(); eItr != eEnd; ++eItr )
    {
        EdgeId e1 = eItr->first, e2 = eItr->second;
        for ( auto fItr=fPairs.begin(),fEnd=fPairs.end(); fItr != fEnd; ++fItr )
        {
            EdgeId f1 = fItr->first, f2 = fItr->second;
	    if ( e1 == f1 || e1 == f2 || e2 == f1 || e2 == f2 )
                continue;

            PathWeights pw1(mHBV,mPathsIndex[e1],e1,f1,f2);
            PathWeights pw2(mHBV,mPathsIndex[e2],e2,f1,f2);

            getReach(pw1,pw2,&reach);
            // make sure {e1,e2,f1,f2} bound a complete sub-graph
            if ( !checkReach(reach,e1,e2,f1,f2) )
                continue;

            int span = std::max(pw1.mSpan,pw2.mSpan);
            if ( pullApart(pw1.weight1,pw1.weight2,
                            pw2.weight1,pw2.weight2,span) )
            {
                pRecs->addRecommendation(pw1.bridges1,pw2.bridges2,reach);
                if ( mVerbosity >= 2 )
                    log(e1,e2,f1,f2,reach,pw1.bridges1,pw2.bridges2,
                        pw1.weight1,pw1.weight2,pw2.weight1,pw2.weight2,
                        span,true);
            }
            else if ( pullApart(pw1.weight2,pw1.weight1,
                                    pw2.weight2,pw2.weight1,span) )
            {
                pRecs->addRecommendation(pw1.bridges2,pw2.bridges1,reach);
                if ( mVerbosity >= 2 )
                    log(e1,e2,f2,f1,reach,pw1.bridges2,pw2.bridges1,
                        pw1.weight2,pw1.weight1,pw2.weight2,pw2.weight1,
                        span,true);
            }
            else if ( mVerbosity >= 3 &&
                    !pw1.bridges1.empty() && !pw2.bridges2.empty() )
            {
                log(e1,e2,f1,f2,reach,pw1.bridges1,pw2.bridges2,
                        pw1.weight1,pw1.weight2,pw2.weight1,pw2.weight2,
                        span,false);
            }
            else if ( mVerbosity >= 3 &&
                    !pw1.bridges2.empty() && !pw2.bridges1.empty() )
            {
                log(e1,e2,f2,f1,reach,pw1.bridges2,pw2.bridges1,
                        pw1.weight2,pw1.weight1,pw2.weight2,pw2.weight1,
                        span,false);
            }
        }
    }
}

// run the vertex scanning in parallel
void findJoins( SupportedHyperBasevector const& hbv, vec<int> const& toLeft,
                    vec<int> const& toRight, int verbosity,
                    double minWeightSplit, Recommendations* pRecs )
{
    JoinDiscoverer jd(hbv,toLeft,toRight,verbosity,minWeightSplit);
    unsigned nThreads = getConfiguredNumThreads();
    size_t nVertices = hbv.N();
    size_t batchSize = ((nVertices+nThreads-1)/nThreads+9)/10;
    parallelForBatch(0ul,nVertices,batchSize,
                        [jd,pRecs]( int vertexId )
                        { jd.processVertex(vertexId,pRecs); },
                        nThreads,verbosity==1);
}

} // end of anonymous namespace

void SupportedHyperBasevector::PullApart2( const double min_weight_split,
        const long_logging& logc )
{
    // Logging stuff.
    if (logc.STATUS_LOGGING) cout << Date() << ": starting PullApart2" << endl;
    int verbosity = logc.verb["PULL_APART2"];
    double clock = WallClockTime();

    // Iterate until no improvement.
    size_t iterationCount = 0;
    bool progress = true;
    while ( progress )
    {
        progress = false;

        if ( verbosity )
            cout << Date() << ": starting iteration "
                    << ++iterationCount << std::endl;
#if 0
        if ( true )
        {
            String dotName = "iteration"+ToString(iterationCount)+".dot";
            std::ofstream out(dotName.c_str());
            PrintSummaryDOT0w(out,false,false,true);
            out.close();
        }
#endif
        vec<int> to_left, to_right;
        ToLeft(to_left); ToRight(to_right);

        Recommendations recs;
        findJoins(*this,to_left,to_right,verbosity,min_weight_split,&recs);
        vec<Join>& joins = recs.mJoins;
        if ( joins.empty() )
            break;

        // Process joins.
        cout << Date() << ": process joins" << std::endl;
        vec<vec<SHBVPath>>& paths1 = recs.mPaths1;
        vec<vec<SHBVPath>>& paths2 = recs.mPaths2;
        ParallelSortSync(joins, paths1, paths2);
        if (logc.STATUS_LOGGING)
        {    cout << Date() << ": processing " << joins.size() << " potential joins"
                  << endl;    }
        if ( verbosity >= 2 )
            cout << "\n";
        vec<Bool> touched(EdgeObjectCount(), False);
        for ( int i = 0; i < joins.isize(); i++ )
        {
            Bool overlap = False;
            for ( int j = 0; j < joins[i].first.isize(); j++ )
                if ( touched[joins[i].first[j]] )
                    overlap = True;
            for ( int j = 0; j < joins[i].second.isize(); j++ )
                if ( touched[joins[i].second[j]] )
                    overlap = True;
            for ( int j = 0; j < joins[i].third.isize(); j++ )
                if ( touched[joins[i].third[j]] )
                    overlap = True;
            if ( overlap )
                continue;

            vec<triple<vec<int>, vec<int>, vec<int> > > proc;
            vec<vec<vec<int> > > proc1, proc2;
            proc.push_back(joins[i]);
            proc1.push_back(paths1[i]), proc2.push_back(paths2[i]);

            if ( Inv(joins[i].first[0]) >= 0 )
            {
                vec<int> a, b, c;
                a.push_back(Inv(joins[i].third[0]), Inv(joins[i].third[1]));
                for ( int j = 0; j < joins[i].second.isize(); j++ )
                    b.push_back(Inv(joins[i].second[j]));
                Sort(b);
                c.push_back(Inv(joins[i].first[0]), Inv(joins[i].first[1]));
                vec<int> all1, all2;
                all1.append(joins[i].first);
                all1.append(joins[i].second);
                all1.append(joins[i].third);
                all2.append(a), all2.append(b), all2.append(c);
                Sort(all1), Sort(all2);
                if ( Meet(all1, all2) )
                    continue;
                for ( int j = 0; j < all2.isize(); j++ )
                    if ( touched[all2[j]] )
                        overlap = True;
                if ( overlap )
                    continue;
                proc.push(a, b, c);
                vec<vec<int> > p1 = paths1[i], p2 = paths2[i];
                for ( int j = 0; j < p1.isize(); j++ )
                {
                    p1[j].ReverseMe();
                    for ( int l = 0; l < p1[j].isize(); l++ )
                        p1[j][l] = Inv(p1[j][l]);
                }
                for ( int j = 0; j < p2.isize(); j++ )
                {
                    p2[j].ReverseMe();
                    for ( int l = 0; l < p2[j].isize(); l++ )
                        p2[j][l] = Inv(p2[j][l]);
                }
                proc1.push_back(p1), proc2.push_back(p2);
            }

            vec<fix64_6> weight(EdgeObjectCount(), 0);
            for ( int i = 0; i < NPaths(); i++ )
                for ( int j = 0; j < Path(i).isize(); j++ )
                    weight[Path(i)[j]] += Weight(i);
            vec<int> count(2, 0);
            for ( int p = 0; p < proc.isize(); p++ )
            {
                if ( verbosity >= 2 )
                {
                    int e1 = proc[p].first[0], e2 = proc[p].first[1];
                    int f1 = proc[p].third[0], f2 = proc[p].third[1];
                    cout << "joining: ";
                    PRINT4(e1, e2, f1, f2);
                }
                for ( int j = 0; j < proc[p].first.isize(); j++ )
                    touched[proc[p].first[j]] = True;
                for ( int j = 0; j < proc[p].second.isize(); j++ )
                    touched[proc[p].second[j]] = True;
                for ( int j = 0; j < proc[p].third.isize(); j++ )
                    touched[proc[p].third[j]] = True;
                vec<int> dels;
                dels.append(proc[p].first);
                dels.append(proc[p].second);
                dels.append(proc[p].third);
                const int max_del_weight = 4;
                Bool bad = False;
                for ( int i = 0; i < dels.isize(); i++ )
                {
                    Bool used = False;
                    for ( int j = 0; j < proc1[p].isize(); j++ )
                        if ( Member(proc1[p][j], dels[i]) )
                            used = True;
                    for ( int j = 0; j < proc2[p].isize(); j++ )
                        if ( Member(proc2[p][j], dels[i]) )
                            used = True;
                    if ( used )
                        continue;
                    if ( weight[dels[i]] > max_del_weight )
                        bad = True;
                }
                if ( bad )
                {
                    if ( verbosity >= 2 )
                        cout << "aborting join" << endl;
                    continue;
                }
                progress = true;
                int v1 = to_left[proc[p].first[0]];
                int v2 = to_left[proc[p].first[1]];
                int w1 = to_right[proc[p].third[0]];
                int w2 = to_right[proc[p].third[1]];
                vec<vec<int> > e(proc1[p]);
                e.append(proc2[p]);
                vec<int> f;
                for ( int j = 0; j < proc1[p].isize(); j++ )
                {
                    f.push_back(EdgeObjectCount());
                    int rid = -1;
                    if ( p == 1 )
                        rid = EdgeObjectCount() - count[0];
                    InvMutable().push_back(rid);
                    if ( p == 1 )
                    {
                        InvMutable(EdgeObjectCount() - count[0]) =
                                EdgeObjectCount();
                    }
                    AddEdge(v1, w1, Cat(proc1[p][j]));
                    count[p]++;
                }
                for ( int j = 0; j < proc2[p].isize(); j++ )
                {
                    f.push_back(EdgeObjectCount());
                    int rid = -1;
                    if ( p == 1 )
                        rid = EdgeObjectCount() - count[0];
                    InvMutable().push_back(rid);
                    if ( p == 1 )
                    {
                        InvMutable(EdgeObjectCount() - count[0]) =
                                EdgeObjectCount();
                    }
                    AddEdge(v2, w2, Cat(proc2[p][j]));
                    count[p]++;
                }
                DeleteEdges(dels);
                TransformPaths(e, f);
            }
        }
        if ( verbosity >= 2 )
            cout << "\n";

        // Clean up.
        cout << Date() << ": clean up" << std::endl;

        UniqueOrderPaths();
        RemoveEdgelessVertices();
        REPORT_TIME( clock, "used in PullApart2");
        RemoveDeadEdgeObjects();
        RemovePathsWithoutReverseComplements();
        RemovePairsWithoutReverseComplements();
        FixWeights(logc);
        TestValid(logc);
        DeleteReverseComplementComponents(logc);
    }
}

void SupportedHyperBasevector::Gulp( const long_logging_control& log_control, 
     const long_logging& logc )
{    double clock1 = WallClockTime( );
     if (logc.STATUS_LOGGING) cout << Date( ) << ": gulping edges" << endl;
     vec<int> to_left, to_right;
     ToLeft(to_left), ToRight(to_right);
     vec< vec< pair<int,int> > > paths_index( EdgeObjectCount( ) );
     for (int i = 0; i < Paths( ).isize(); i++)    
     for (int j = 0; j < Path(i).isize(); j++)
          paths_index[ Path(i,j) ].push(i, j);
     const int max_gulp = 20;
     while(1)
     {    if ( logc.verb[ "GULP" ] >= 1 ) 
               cout << Date( ) << ": starting gulp iteration" << endl;
          Bool changed = False;
          for ( int v = 0; v < N( ); v++ )
          {    vec<int> dels, rdels, out, rout;
               vec< vec<int> > in(2), rin(2);
               if ( To(v).size( ) == 1 && From(v).size( ) == 2 )
               {    int u = To(v)[0], w1 = From(v)[0], w2 = From(v)[1];
                    int e = EdgeObjectIndexByIndexTo( v, 0 );
                    if ( EdgeLengthKmers(e) > max_gulp ) continue;
                    int f1 = EdgeObjectIndexByIndexFrom( v, 0 );
                    int f2 = EdgeObjectIndexByIndexFrom( v, 1 );
                    // u --e--> v --f1,f2--> w1,w2
                    dels.push_back( e, f1, f2 );
                    int re = Inv(e), rf1 = Inv(f1), rf2 = Inv(f2);
                    rdels.push_back( re, rf1, rf2 );
                    UniqueSort(dels), UniqueSort(rdels);
                    if ( dels.size( ) != 3 ) continue;
                    if ( InvDef(e) && Meet( dels, rdels ) ) continue;
                    int ef1 = EdgeObjectCount( ), ef2 = EdgeObjectCount( ) + 1;
                    changed = True;
                    AddEdge( u, w1, Cat(e,f1) ), AddEdge( u, w2, Cat(e,f2) );
                    DeleteEdges(dels);
                    in[0].push_back(e, f1), in[1].push_back(e, f2);
                    out.push_back( ef1, ef2 );
                    TransformPaths( in, out /*, paths_index */ );
                    to_left.push_back(u, u), to_right.push_back(w1, w2);
                    if ( !InvDef(e) ) InvMutable( ).push_back( -1, -1 );
                    else
                    {    int rf1e = EdgeObjectCount( ); 
                         int rf2e = EdgeObjectCount( ) + 1;
                         int ru = to_right[re]; 
                         int rw1 = to_left[rf1], rw2 = to_left[rf2];
                         AddEdge( rw1, ru, Cat(rf1,re) ); 
                         AddEdge( rw2, ru, Cat(rf2,re) );
                         DeleteEdges(rdels);
                         in[0].clear( ), in[1].clear( ), out.clear( );
                         in[0].push_back(rf1, re), in[1].push_back(rf2, re);
                         out.push_back( rf1e, rf2e );
                         TransformPaths( in, out /*, paths_index */ );
                         to_left.push_back(rw1, rw2), to_right.push_back(ru, ru);
                         InvMutable( ).push_back( rf1e, rf2e, ef1, ef2 );    }    }
               if ( To(v).size( ) == 2 & From(v).size( ) == 1 )
               {    int u1 = To(v)[0], u2 = To(v)[1], w = From(v)[0];
                    int e1 = EdgeObjectIndexByIndexTo( v, 0 );
                    int e2 = EdgeObjectIndexByIndexTo( v, 1 );
                    int f = EdgeObjectIndexByIndexFrom( v, 0 );
                    if ( EdgeLengthKmers(f) > max_gulp ) continue;
                    // u1,u2 --e1,e2--> v --f--> w
                    dels.push_back( e1, e2, f );
                    int re1 = Inv(e1), re2 = Inv(e2), rf = Inv(f);
                    rdels.push_back( re1, re2, rf );
                    UniqueSort(dels), UniqueSort(rdels);
                    if ( dels.size( ) != 3 ) continue;
                    if ( InvDef(f) && Meet( dels, rdels ) ) continue;
                    int e1f = EdgeObjectCount( ), e2f = EdgeObjectCount( ) + 1;
                    changed = True;
                    AddEdge( u1, w, Cat(e1,f) ), AddEdge( u2, w, Cat(e2,f) );
                    DeleteEdges(dels);
                    in[0].push_back(e1, f), in[1].push_back(e2, f);
                    out.push_back( e1f, e2f );
                    TransformPaths( in, out /*, paths_index */ );
                    to_left.push_back(u1, u2), to_right.push_back(w, w);
                    if ( !InvDef(f) ) InvMutable( ).push_back( -1, -1 );
                    else
                    {    int rfe1 = EdgeObjectCount( ); 
                         int rfe2 = EdgeObjectCount( ) + 1;
                         DeleteEdges(rdels);
                         int ru1 = to_right[re1], ru2 = to_right[re2];
                         int rw = to_left[rf];
                         AddEdge( rw, ru1, Cat(rf,re1) );
                         AddEdge( rw, ru2, Cat(rf,re2) );
                         in[0].clear( ), in[1].clear( ), out.clear( );
                         in[0].push_back(rf, re1), in[1].push_back(rf, re2);
                         out.push_back( rfe1, rfe2 );
                         TransformPaths( in, out /*, paths_index */ );
                         to_left.push_back(rw, rw), to_right.push_back(ru1, ru2);
                         InvMutable( ).push_back( 
                              rfe1, rfe2, e1f, e2f );    }    }    }

          if ( !changed ) break;    }
     UniqueOrderPaths( );
     RemoveEdgelessVertices( );
     REPORT_TIME( clock1, "used in gulping 1" );
     RemoveDeadEdgeObjects( );
     FixWeights(logc);
     if (logc.STATUS_LOGGING)
     {    cout << Date( ) << ": after gulping edges there are "
               << EdgeObjectCount( ) << " edges" << endl;    }
     TestValid(logc);    }

void SupportedHyperBasevector::Ungulp( const long_logging& logc )
{    double clock = WallClockTime( );
     if (logc.STATUS_LOGGING) cout << Date( ) << ": ungulping edges" << endl;
     const int max_gulp = 20;
     vec<int> to_left, to_right;
     ToLeft(to_left), ToRight(to_right);
     int nv = N( );
     vec<int> ins;
     vec< vec<int> > outs;
     for ( int v = 0; v < nv; v++ )
     {    if ( From(v).size( ) != 2 ) continue;
          int w1 = From(v)[0], w2 = From(v)[1];
          if ( w1 == v || w2 == v ) continue;
          int e1 = EdgeObjectIndexByIndexFrom( v, 0 );
          int e2 = EdgeObjectIndexByIndexFrom( v, 1 );
          const basevector &E1 = EdgeObject(e1), &E2 = EdgeObject(e2);
          int re1 = Inv(e1), re2 = Inv(e2);
          int rv, rw1, rw2;
          if ( re1 >= 0 )
          {    rv = to_right[re1], rw1 = to_left[re1], rw2 = to_left[re2];
               vec<int> alle;
               alle.push_back( e1, e2, re1, re2 );
               UniqueSort(alle);
               if ( alle.size( ) != 4 ) continue;
               vec<int> vert1( {v,w1,w2} ), vert2( {rv,rw1,rw2} );
               UniqueSort(vert1), UniqueSort(vert2);
               if ( Meet( vert1, vert2 ) ) continue;   }
          int n;
          for ( n = 0; n < Min( E1.isize( ), E2.isize( ) ); n++ )
               if ( E1[n] != E2[n] ) break;
          if ( n == E1.isize( ) || n == E2.isize( ) ) n--;
          if ( n < K( ) + max_gulp - 1 ) continue;
          // PRINT4( v, E1.size( ), E2.size( ), n ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          basevector c( E1, 0, n );
          basevector f1( E1, n - K( ) + 1, E1.isize( ) - ( n - K( ) + 1 ) );
          basevector f2( E2, n - K( ) + 1, E2.isize( ) - ( n - K( ) + 1 ) );
          // PRINT2( f1.size( ), f2.size( ) ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          // e1 = f1, e2 = f2;
          int x = N( );
          // PRINT(x); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          vec<int> us = {0,1};
          DeleteEdgesFrom( v, us );
          AddVertices(1);
          AddEdge( v, x, c );
          to_left.push_back(v), to_right.push_back(x);
          AddEdge( x, w1, f1 );
          to_left.push_back(x), to_right.push_back(w1);
          AddEdge( x, w2, f2 );
          to_left.push_back(x), to_right.push_back(w2);
          ins.push_back( e1, e2 );
          {    int N = EdgeObjectCount( );
               outs.push_back( {N-3,N-2} );
               outs.push_back( {N-3,N-1} );    }
          // HyperBasevector::TestValid( ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          if ( re1 < 0 ) InvMutable( ).push_back( -1, -1, -1 );
          else
          {    basevector rc(c), rf1(f1), rf2(f2);
               rc.ReverseComplement( );
               rf1.ReverseComplement( ), rf2.ReverseComplement( );
               int rx = N( );
               DeleteEdgesTo( rv, us );
               AddVertices(1);
               AddEdge( rx, rv, rc );
               to_left.push_back(rx), to_right.push_back(rv);
               AddEdge( rw1, rx, rf1 );
               to_left.push_back(rw1), to_right.push_back(rx);
               AddEdge( rw2, rx, rf2 );
               to_left.push_back(rw2), to_right.push_back(rx);
               int N = EdgeObjectCount( );
               InvMutable( ).push_back( N-3, N-2, N-1, N-6, N-5, N-4 );
               ins.push_back( re1, re2 );
               outs.push_back( {N-1,N-3} );
               outs.push_back( {N-2,N-3} );
               // HyperBasevector::TestValid( ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                    }    }
     for ( int w = 0; w < nv; w++ )
     {    if ( To(w).size( ) != 2 ) continue;
          int v1 = To(w)[0], v2 = To(w)[1];
          if ( v1 == w || v2 == w ) continue;
          int e1 = EdgeObjectIndexByIndexTo( w, 0 );
          if ( InvDef(e1) ) continue;
          int e2 = EdgeObjectIndexByIndexTo( w, 1 );
          const basevector &E1 = EdgeObject(e1), &E2 = EdgeObject(e2);
          int n1 = E1.size( ), n2 = E2.size( );
          int n;
          for ( n = 0; n < Min( E1.isize( ), E2.isize( ) ); n++ )
               if ( E1[ n1 - n - 1 ] != E2[ n2 - n - 1 ] ) break;
          if ( n == E1.isize( ) || n == E2.isize( ) ) n--;
          if ( n < K( ) + max_gulp - 1 ) continue;
          // PRINT4( w, E1.size( ), E2.size( ), n ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          basevector c( E1, n1 - n, n );
          basevector f1( E1, 0, n1 - n + K( ) - 1 );
          basevector f2( E2, 0, n2 - n + K( ) - 1 );
          // PRINT2( f1.size( ), f2.size( ) ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          int x = N( );
          // PRINT(x); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          vec<int> us;
          us.push_back( 0, 1 );
          DeleteEdgesTo( w, us );
          AddVertices(1);
          AddEdge( x, w, c );
          AddEdge( v1, x, f1 );
          AddEdge( v2, x, f2 );
          int N = EdgeObjectCount( );
          ins.push_back(e1, e2);
          outs.push_back( {N-1,N-3} );
          outs.push_back( {N-2,N-3} );
          InvMutable( ).push_back( -1, -1, -1 );
          // HyperBasevector::TestValid( ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
               }

     SortSync( ins, outs );
     for ( int i = 0; i < NPaths( ); i++ )
     for ( int j = 0; j < Path(i).isize( ); j++ )
     {    int z = BinPosition( ins, Path(i)[j] );
          if ( z < 0 ) continue;
     
          vec<int> p;
          for ( int k = 0; k < j; k++ )
               p.push_back( Path(i)[k] );
          p.append( outs[z] );
          for ( int k = j + 1; k < Path(i).isize( ); k++ )
               p.push_back( Path(i)[k] );    }
     // Note not fixing pairs.

     // TestValid( logc, False ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     REPORT_TIME( clock, "used in ungulping" );
     RemoveUnneededVertices( );
     RemoveDeadEdgeObjects( );
     if (logc.STATUS_LOGGING) 
     {    cout << Date( ) << ": after ungulping there are " << EdgeObjectCount( ) 
               << " edges" << endl;    }
     TestValid(logc);    }

void SupportedHyperBasevector::TrimHangingEnds( const int max_del,
     const double junk_ratio, const long_heuristics& heur, 
     const long_logging& logc )
{
     double xclock = WallClockTime( );
     if (logc.STATUS_LOGGING)
     {    cout << Date( ) << ": see " << EdgeObjectCount( ) 
               << " edges, removing hanging ends" << endl;    }
     vec<kmer_count> kc( EdgeObjectCount( ) );
     for ( int e = 0; e < EdgeObjectCount( ); e++ )
          kc[e].n = EdgeObject(e).isize( ) - K( ) + 1;
     digraphE<kmer_count> shb_kc( *this, kc );
     REPORT_TIME( xclock, "used removing hanging ends head" );
     double pclock = WallClockTime( );
     const int max_paths = 100;
     RemoveHangingEnds3( shb_kc, &kmer_count::N, max_del, junk_ratio, max_paths );
     REPORT_TIME( pclock, "used removing hanging ends" );
     double q1clock = WallClockTime( );
     vec<int> e_to_delete;
     vec<Bool> used;
     shb_kc.Used(used);
     for ( int e = 0; e < EdgeObjectCount( ); e++ )
     {    if ( !used[e] )
          {    e_to_delete.push_back(e);
               if ( Inv(e) >= 0 && used[ Inv(e) ] )
                    e_to_delete.push_back( Inv(e) );    }    }
     if ( logc.verb[ "TRIM" ] >= 1 )
     {    cout << "\nTrimHangingEnds, deleting edges " << printSeq(e_to_delete)
               << endl;    }
     DeleteEdges(e_to_delete);
     REPORT_TIME( q1clock, "used removing hanging ends tail 1" );
     TruncatePaths(logc);
     if ( heur.EXTRA_STEPS )
     {    double rclock = WallClockTime( );
          RemoveUnneededVertices( );
          REPORT_TIME( rclock, "using rh tail 2" );    }
     RemoveUnneededVertices( );
     RemoveDeadEdgeObjects( );
     FixWeights(logc);
     TestValid(logc);

     // Remove terminal loops of length <= 50 kmers.

     double rclock2 = WallClockTime( );
     vec<int> ldels;
     const int maxl = 50;
     for ( int v = 0; v < N( ); v++ )
     {    if ( To(v).size( ) != 2 || From(v).size( ) != 1 ) continue;
          if ( From(v)[0] != v ) continue;
          int e = EdgeObjectIndexByIndexFrom( v, 0 );
          if ( EdgeLengthKmers(e) > maxl ) continue;
          ldels.push_back(e);    }
     for ( int v = 0; v < N( ); v++ )
     {    if ( From(v).size( ) != 2 || To(v).size( ) != 1 ) continue;
          if ( To(v)[0] != v ) continue;
          int e = EdgeObjectIndexByIndexTo( v, 0 );
          if ( EdgeLengthKmers(e) > maxl ) continue;
          ldels.push_back(e);    }
     DeleteEdges(ldels);
     TruncatePaths(logc);
     REPORT_TIME( rclock2, "used removing hanging ends tail 2" );
     RemoveDeadEdgeObjects( );
     FixWeights(logc);
     TestValid(logc);    }

void SupportedHyperBasevector::DeleteLowCoverage( const long_heuristics& heur,
     const long_logging_control& log_control, const long_logging& logc )
{
     double mclock = WallClockTime( );
     vec<int> to_left, to_right, to_delete;
     ToLeft(to_left), ToRight(to_right);
     const double low_cov = 2.0;
     vec<fix64_6> cov( EdgeObjectCount( ), 0.0 );
     for ( int i = 0; i < NPaths( ); i++ )
     {    for ( int j = 0; j < Path(i).isize( ); j++ )
          cov[ Path(i)[j] ] += Weight(i);    }

     const double pceil = 0.99;
     const double minp = 0.0001;
     vec< vec<fix64_6> > weights( EdgeObjectCount( ) );

     vec< triple<int64_t,int,fix64_6> > wid;
     vec<double> p( EdgeObjectCount( ), 1.0 );
     if ( heur.NEW_LC_FILT )
     {    for ( int i = 0; i < NPaths( ); i++ )
          {    for ( int k = 0; k < WeightFwOrigin(i).isize( ); k++ )
               {    wid.push( WeightFwOrigin(i)[k].second, i, 
                         WeightFwOrigin(i)[k].first );    }
               for ( int k = 0; k < WeightRcOrigin(i).isize( ); k++ )
               {    wid.push( WeightRcOrigin(i)[k].second, i, 
                         WeightRcOrigin(i)[k].first );    }    }
          ParallelSort(wid);
          for ( int64_t i = 0; i < wid.jsize( ); i++ )
          {    int64_t j;
               for ( j = i + 1; j < wid.jsize( ); j++ )
                    if ( wid[j].first != wid[i].first ) break;
               vec< pair<int,fix64_6> > w;
               for ( int64_t k = i; k < j; k++ )
               {    vec<int> p = Path( wid[k].second );
                    UniqueSort(p);
                    for ( int l = 0; l < p.isize( ); l++ )
                         w.push( p[l], wid[k].third );    }
               Sort(w);
               for ( int r = 0; r < w.isize( ); r++ )
               {    int s;
                    for ( s = r + 1; s < w.isize( ); s++ )
                         if ( w[s].first != w[r].first ) break;
                    fix64_6 x = 0;
                    for ( int t = r; t < s; t++ )
                         x += w[t].second;
                    weights[ w[r].first ].push_back(x);
                    r = s - 1;    }
               i = j - 1;    }
          for ( int e = 0; e < EdgeObjectCount( ); e++ )
          {    Sort( weights[e] );
               for ( int j = 0; j < weights[e].isize( ); j++ )
               {    if ( p[e] < minp ) break;
                    p[e] *= 1.0 
                         - Min( pceil, weights[e][j].ToDouble( ) );    }    }    }

     const int min_mult = 5;
     for ( int e = 0; e < EdgeObjectCount( ); e++ )
     {    int re = Inv(e), v = to_left[e], w = to_right[e];
          fix64_6 c = cov[e];
          fix64_6 rc = ( re < 0 ? 1000000000 : cov[re] );
          fix64_6 alt_c = 0, alt_rc = 0;
          for ( int j = 0; j < From(v).isize( ); j++ )
               alt_c = Max( alt_c, cov[ EdgeObjectIndexByIndexFrom( v, j ) ] );
          for ( int j = 0; j < To(w).isize( ); j++ )
               alt_c = Max( alt_c, cov[ EdgeObjectIndexByIndexTo( w, j ) ] );
          if ( re >= 0 )
          {    int v = to_left[re], w = to_right[re];
               for ( int j = 0; j < From(v).isize( ); j++ )
               {    alt_rc = Max( alt_rc, cov[ 
                         EdgeObjectIndexByIndexFrom( v, j ) ] );    }
               for ( int j = 0; j < To(w).isize( ); j++ )
               {    alt_rc = Max( alt_rc, cov[ 
                         EdgeObjectIndexByIndexTo( w, j ) ] );    }    }
          if ( heur.NEW_LC_FILT )
          {    if ( ( ( c <= low_cov || p[e] >= minp ) 
                         && alt_c >= min_mult * c ) 
                    || ( ( rc <= low_cov || ( re >= 0 && p[re] >= minp ) ) 
                         && alt_rc >= min_mult * rc ) )
               {    if ( logc.verb[ "LOW_COV" ] >= 1 )
                    {    cout << "deleting low-coverage edge " << e 
                              << " (c = " << c << ", alt_c = " << alt_c
                              << ", p[e] = " << p[e] << ")" << endl;    }
                    to_delete.push_back(e);    }    }
          else
          {    if ( heur.LC_CAREFUL && alt_c < low_cov ) continue;
               if ( ( c <= low_cov && alt_c >= min_mult * c )
                    || ( rc <= low_cov && alt_rc >= min_mult * rc ) )
               {    if ( logc.verb[ "LOW_COV" ] >= 1 )
                    {    vec<int> comp;
                         for ( int j = 0; j < From(v).isize( ); j++ )
                         {    int e = EdgeObjectIndexByIndexFrom( v, j );
                              if ( cov[e] == alt_c ) comp.push_back(e);    }
                         for ( int j = 0; j < To(w).isize( ); j++ )
                         {    int e = EdgeObjectIndexByIndexTo( w, j );
                              if ( cov[e] == alt_c ) comp.push_back(e);    }
                         UniqueSort(comp);
                         cout << "deleting low-coverage edge " << e 
                              << " using competing edge(s) " << printSeq(comp)
                              << ", c = " << c << ", alt_c = " << alt_c << endl;    }
                    to_delete.push_back(e);    }    }    } 
     DeleteEdges(to_delete);
     vec<Bool> p_to_delete( NPaths( ), False );
     for ( int i = 0; i < NPaths( ); i++ )
     {    Bool OK = True;
          for ( int j = 0; j < Path(i).isize( ); j++ )
               if ( BinMember( to_delete, Path(i)[j] ) ) OK = False;
          if ( !OK ) p_to_delete[i] = True;    }
     EraseIf( PathsMutable( ), p_to_delete );
     // Temporary if.
     if ( WeightsFwOrigin( ).size( ) == WeightsFw( ).size( ) )
     {    EraseIf( WeightsFwOriginMutable( ), p_to_delete );
          EraseIf( WeightsRcOriginMutable( ), p_to_delete );    }
     EraseIf( WeightsFwMutable( ), p_to_delete );
     EraseIf( WeightsRcMutable( ), p_to_delete );
     p_to_delete.resize_and_set( NPairs( ), False );
     for ( int i = 0; i < NPairs( ); i++ )
     {    Bool OK = True;
          for ( int j = 0; j < PairLeft(i).isize( ); j++ )
               if ( BinMember( to_delete, PairLeft(i)[j] ) ) OK = False;
          if ( !OK ) p_to_delete[i] = True;    
          for ( int j = 0; j < PairRight(i).isize( ); j++ )
               if ( BinMember( to_delete, PairRight(i)[j] ) ) OK = False;
          if ( !OK ) p_to_delete[i] = True;    }
     EraseIf( PairsMutable( ), p_to_delete );
     EraseIf( PairDataMutable( ), p_to_delete );
     RemoveEdgelessVertices( );
     RemoveUnneededVertices( );
     REPORT_TIME( mclock, "used deleting low-coverage edges" );
     RemoveDeadEdgeObjects( );    }
