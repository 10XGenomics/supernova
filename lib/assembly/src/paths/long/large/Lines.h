///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef GAPTOY_LINES_H
#define GAPTOY_LINES_H

#include "CoreTools.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"

// Description of lines data structure:
//
// vec^4<int> = vec<vec<vec<vec<int>>>> lines.
//
// Each element is a line, a vec^3<int>.
//
// Each element in a line is a vec^2<int>, which is a list
// of paths through the graph.  Alternating entries, starting
// with the first have just ONE path.
//
// So this looks like a sequence of bubbles:
//
// ----======-----=======----=====-----
//
// in which the number of branches (or paths) in a bubble is
// typically two but can be larger, or zero for a gap.
//
// Note that circles have a funny representation.  If the circle is a single
// edge, then the line has that entry alone.  However if the circle contains
// more edges (e.g. three, with two in a bubble), then the line has the same
// starting and ending edges (i.e. the same edge appears twice).  The fasta/efasta
// representation of lines in DumpLineFiles takes account of this.

typedef vec<int> LinePath; // elements are edgeIds
typedef vec<LinePath> LineSegment; // aka, a Cell
typedef vec<LineSegment> Line;
typedef vec<Line> LineVec;

// constraint: Lines have an odd number of segments
// constraint: even numbered segments have only 1 Path
// currently, nobody enforces these constraints, but code assumes them

// Counter for coverage.

class covcount {

     public:

     covcount( ) { cov_ = -1; }
     covcount( const float x ) { cov_ = x; }
     void Set( const float x ) { cov_ = x; }
     float Cov( ) const { return cov_; }
     Bool Def( ) const { return cov_ >= 0; }

     private:

     float cov_;

};

TRIVIALLY_SERIALIZABLE(covcount);

template <class DE>
void FindLines( const DE& dgraph, const vec<int>& dinv,
     vec<vec<vec<vec<int>>>>& lines, const int max_cell_paths, const int max_depth,
     const Bool verbose = False );

int64_t LineN50( const HyperBasevector& hb,
     const vec<vec<vec<vec<int>>>>& lines, const int min_len );

template <class EdgeRuler> // given an edge ID, returns its length
int GetPathLength( LinePath const& path, EdgeRuler ruler, int upTo = -1 ) {
    int sum = 0;
    for ( int edgeId : path ) {
        if ( upTo != -1 && edgeId == upTo ) break;
        sum += ruler(edgeId);
    }
    return sum;
}

template <class EdgeRuler>
int GetSegmentLength( LineSegment const& seg, EdgeRuler ruler )
{ size_t nPaths = seg.size(); ForceAssertGt(nPaths,0ul);
  if ( nPaths == 1ul ) return GetPathLength(seg.front(),ruler);
  if ( nPaths == 2ul )
    return (GetPathLength(seg.front(),ruler)+GetPathLength(seg.back(),ruler))/2;
  vec<int> pathLens; pathLens.reserve(seg.size());
  std::transform(seg.begin(),seg.end(),std::back_inserter(pathLens),
          [ruler]( LinePath const& path ){ return GetPathLength(path,ruler); });
  // return the median segment path length
  std::sort(pathLens.begin(),pathLens.end());
  if ( nPaths&1ul ) return pathLens[nPaths/2];
  return (pathLens[nPaths/2]+pathLens[nPaths/2-1])/2; }

template <class EdgeRuler>
int GetLineLength( Line const& line, EdgeRuler ruler )
{ int sum = 0;
  for ( LineSegment const& seg : line ) sum += GetSegmentLength(seg,ruler);
  return sum; }

template <class EdgeRuler>
void GetLineLengths( LineVec const& lines, vec<int>& llens, EdgeRuler ruler )
{ llens.clear(); llens.reserve(lines.size());
  std::transform(lines.begin(),lines.end(),std::back_inserter(llens),
          [ruler]( Line const& line ){ return GetLineLength(line,ruler); }); }

struct HBVKmerRuler
{ HBVKmerRuler(HyperBasevector const& hbv) : mHBV(hbv) {}
  int operator()( int edgeId ) { return mHBV.Kmers(edgeId); }
  HyperBasevector const& mHBV; };

struct HBVXKmerRuler
{ HBVXKmerRuler(HyperBasevectorX const& hbvx ) : mHBVX(hbvx) {}
  int operator()( int edgeId ) { return mHBVX.Kmers(edgeId); }
  HyperBasevectorX const& mHBVX; };

inline void GetLineLengths( HyperBasevectorX const& hb, LineVec const& lines,
                                vec<int>& llens )
{ GetLineLengths(lines,llens,HBVXKmerRuler(hb)); }

struct SupergraphKmerRuler
{
     SupergraphKmerRuler( digraphE<vec<int>> const& de, HyperBasevectorX const& hbx) : mDE(de), mKmerRuler(hbx) {};
     int operator()( int edgeId ) {
          if ( mDE.O(edgeId)[0] < 0 )
               return 0;
          else {
               int sum = 0;
               for ( auto const edge : mDE.O(edgeId) )
                    sum += mKmerRuler( edge );
               return sum;
          }
          /* NOTREACHED */
     }

     digraphE<vec<int>> const& mDE;
     HBVXKmerRuler mKmerRuler;
};

inline void GetLineLengths( HyperBasevector const& hb, LineVec const& lines,
                                vec<int>& llens )
{ GetLineLengths(lines,llens,HBVKmerRuler(hb)); }

inline void GetLineLengths( digraphE<vec<int>> const& de, HyperBasevectorX const& hbx,
          LineVec const& lines, vec<int>& llens ) {
     GetLineLengths( lines, llens, SupergraphKmerRuler( de, hbx ) );
}

void GetTol( const HyperBasevector& hb,
     const vec<vec<vec<vec<int>>>>& lines, vec<int>& tol );

void GetTol( const HyperBasevectorX& hb,
     const vec<vec<vec<vec<int>>>>& lines, vec<int>& tol );

void GetLineNpairs( const HyperBasevector& hb, const vec<int>& inv,
     const ReadPathVec& paths, const vec<vec<vec<vec<int>>>>& lines,
     vec<int>& npairs );

void WriteLineStats( const String& head, const vec<vec<vec<vec<int>>>>& lines,
     const vec<int>& llens, const vec<int>& npairs, const vec<vec<covcount>>& covs );

// Write the line to standard out

void DisplayLine(const vec<vec<vec<int>>>& line);

vec<int64_t> Pids( const int e, const int ss, const HyperBasevector& hb,
     const vec<int>& inv, const ReadPathVec& paths, const VecULongVec& paths_index,
     const vec<int64_t>& subsam_starts );

vec<int64_t> Pids( const int e, const HyperBasevector& hb,
     const vec<int>& inv, const ReadPathVec& paths, const VecULongVec& paths_index );

vec<int64_t> LinePids( const vec<vec<vec<int>>>& L, const HyperBasevector& hb,
     const vec<int>& inv, const ReadPathVec& paths, const VecULongVec& paths_index );

double RawCoverage( const int e, const int ss, const HyperBasevector& hb,
     const vec<int>& inv, const ReadPathVec& paths, const VecULongVec& paths_index,
     const vec<int64_t>& subsam_starts );

double RawCoverage( const int e, const HyperBasevector& hb,
     const vec<int>& inv, const ReadPathVec& paths, const VecULongVec& paths_index );

void ComputeCoverage( const HyperBasevector& hb, const vec<int>& inv,
     const ReadPathVec& paths, const vec<vec<vec<vec<int>>>>& lines,
     const vec<int64_t>& subsam_starts, vec<vec<covcount>>& covs );

void LineInv( const vec<vec<vec<vec<int>>>>& lines, const vec<int>& inv,
     vec<int>& linv );

// Sort lines so that they are in reverse order by length, and each line is
// adjacent to its reverse complement.

void SortLines( vec<vec<vec<vec<int>>>>& lines, const HyperBasevector& hb,
     const vec<int>& inv );

void DumpLineFiles( const vec<vec<vec<vec<int>>>>& lines, const HyperBasevector& hb,
     const vec<int>& inv, const ReadPathVec& paths, const String& dir );

// Split a line into contigs.

void MakeTigs( const vec<vec<vec<int>>>& L, vec<vec<vec<vec<int>>>>& tigs );

void Canonicalize( vec<vec<vec<vec<int>>>>& L );

#endif
