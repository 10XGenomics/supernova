///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "paths/long/ExtendReadPath.h"



namespace {     // ANONYMOUS NAMESPACE

unsigned scoreRightOverlap( bvec const& bases, qvec const& quals,
						   size_t start, bvec const& edge, int K,
						   double pDecay,
						   qual_t mapQ2,
						   qual_t pLeftOver)
{
    //      cout << "hanging bit of read: " << endl;
    //      bases.PrintBases(cout, bases.size() - start, start, false, 80 );
    //      cout << "possibly matching edge: " <<  endl;
    //      int offset = edge.size() - start - K + 1;
    ////    if ( offset > 0 ) {
    //          edge.PrintBases(cout, K-1, K-1+start, false, 80 );
    ////    } else {
    //          cout << "WHOLE EDGE:" << endl;
    //          edge.Print(cout );
    ////    }
    ForceAssertLe(pDecay,1.0);
    ForceAssertGe(pDecay,0.0);

    ForceAssertLt( start, bases.size() );
    ForceAssertEq( bases.size(), quals.size() );

    auto bitr = bases.end() - start;
    auto qitr = quals.end() - start;
    auto eitr = edge.begin() + (K-1);

    unsigned qSum = 0;
    unsigned penalty = 0;
    while ( bitr != bases.end() && qitr != quals.end() && eitr != edge.end() ) {
	if ( *bitr != *eitr ) {
	    // transform Q2 -> Q20
	    auto qscore = ( *qitr == 2 ) ? mapQ2 : *qitr;
	    penalty += qscore;
	    qSum += penalty;
	} else if ( penalty > 0 ) {
	    penalty -= (pDecay*penalty);
	}
	//          PRINT5(Base::val2Char(*bitr), Base::val2Char(*eitr), (int)*qitr, qSum, penalty);
	bitr++;
	qitr++;
	eitr++;
    }

    // penalize left over bases on the read
    while ( bitr++ != bases.end() && qitr++ != quals.end() )
	qSum += pLeftOver;

    return qSum;
}


unsigned scoreLeftOverlap( bvec const& bases, qvec const& quals,
						  size_t start, bvec const& edge, int K,
						   double pDecay,
						   qual_t mapQ2,
						   qual_t pLeftOver)
{
    //      cout << "hanging bit of read: " << endl;
    //      bases.PrintBases(cout, 0, start, false, 80 );
    //      cout << "possibly matching edge: " <<  endl;
    //      int offset = edge.size() - start - K + 1;
    //      if ( offset > 0 ) {
    //          PRINT3(offset, start, K-1);
    //          edge.PrintBases(cout, offset, start, false, 80 );
    //      } else {
    //          cout << "WHOLE EDGE:" << endl;
    //          edge.Print(cout );
    //      }
    ForceAssertLe(pDecay,1.0);
    ForceAssertGe(pDecay,0.0);

    ForceAssertLt( start, bases.size() );
    ForceAssertEq( bases.size(), quals.size() );

    auto bitr = bases.rend() - start;
    auto qitr = quals.rend() - start;
    auto eitr = edge.rbegin() + (K-1);

    unsigned qSum = 0;
    unsigned penalty = 0;
    while ( bitr != bases.rend() && qitr != quals.rend() && eitr != edge.rend() ) {
	if ( *bitr != *eitr ) {
	    // transform Q2 -> Q20
	    auto qscore = ( *qitr == 2 ) ? mapQ2 : *qitr;
	    penalty += qscore;
	    qSum += penalty;
	} else if ( penalty > 0 ) {
	    penalty -= (pDecay*penalty);
	}
	//          PRINT5(Base::val2Char(*bitr), Base::val2Char(*eitr), (int)*qitr, qSum, penalty);
	bitr++;
	qitr++;
	eitr++;
    }

    // penalize left over bases on the read
    while ( bitr++ != bases.rend() && qitr++ != quals.rend() )
	qSum += pLeftOver;

    return qSum;
}

};      /// END OF ANONYMOUS NAMESPACE

/////////// CLASS METHODS //////////////

void ExtendReadPath::attemptLeftRightExtension(  ReadPath& path, basevector const& bases, qualvector const& quals )
{
     ForceAssertEq(quals.size(), bases.size());

    while( attemptLeftwardExtension( path, bases, quals) );
     ForceAssertEq(quals.size(), bases.size());
    while( attemptRightwardExtension( path, bases, quals) );
     ForceAssertEq(quals.size(), bases.size());
}



bool ExtendReadPath::attemptLeftwardExtension(  ReadPath& path,
                        basevector const& bases, qualvector const& quals )
{
    if ( !path.size() ) return false;

    if (mDebug)
	cout << "Leftward from: " << path.getOffset() << ":" << printSeq( path ) << endl;

    if ( path.getOffset() >= 0 ) return false;

    size_t lastGap = -path.getOffset();

    if ( lastGap < 10 ) return false;

    ForceAssertLt(path.getOffset(), 0);

    // cast to be compatible with the rest of the code, now that we
    // know that we're non-negative

    // ********* Should be path.front() ***************
    size_t hbv_edge_id = path.front();
    vec<int> const& to_left = getToLeft();        // ensure that mpToLeft exists
    size_t vleft = to_left[hbv_edge_id];

    auto const& edges = mHBV.ToEdgeObj(vleft);
    auto const& vdest = mHBV.To(vleft);

    int K = mHBV.K();
    vec<int> short_dest;
    vec<bool> ehanging(edges.size(),false);
    vec<bool> elong(edges.size(), false);

    // classify outgoing edges into hanging or not and short or long
    // a short edge is one not capable of "finishing" the read and we ignore
    // short, hanging edges.
    for ( size_t i = 0; i < edges.size(); ++i ) {
	if ( mHBV.ToSize(vdest[i]) == 0 && mHBV.FromSize(vdest[i]) == 1  )
	    ehanging[i] = true;

	if ( mHBV.EdgeObject(edges[i]).size() - (K-1) >= lastGap ) elong[i] = true;

	// push back the destination vertex of a short edge
	if ( !elong[i] && !ehanging[i] )
	    short_dest.push_back( vdest[i] );
    }

    // we only process "short" edges in the following cases:
    // 1. There is only one edge (short or long), or
    //   a. there are no "long" edges mixed in, and
    //   b. all short edges lead to the same vertex, and
    //   c. there is only one outgoing edge from the vertex to which all
    //      short edges lead.
    if ( !edges.solo() ) {
	size_t nlong = Sum(elong);
	if ( short_dest.size() > 0 ) {
	    if ( nlong > 0 ) return false;
	    UniqueSort(short_dest);
	    if ( !short_dest.solo() ) return false;
	    if ( mHBV.ToSize( short_dest.back() ) != 1 ) return false;
	}
    }


    // okay, if we've made it here, we're going to score all
    // edges except hanging ends, unless all edges ARE hanging ends
    int least_edge = -1;
    unsigned least = std::numeric_limits<unsigned>::max();
    for ( size_t i = 0; i < edges.size(); ++i ) {
	//            if ( nlong == edges.size() || !ehanging[i] )
	if ( !ehanging[i] || edges.solo() ) {
	    auto score = scoreLeftOverlap( bases, quals, lastGap, mHBV.EdgeObject(edges[i]), mHBV.K(),
	            mPenaltyDecay, mMapQ2, mLeftOverPenalty);
	    if ( score < least ) {
		least_edge = edges[i];
		least = score;
	    }
	}
    }

    if (mDebug) {
        cout << "least_edge = " << least_edge << ", score of " << least <<
                  string( ( least <= lastGap*10 ) ? " is" : " is NOT" ) <<
                            " less than lastGap*10 = " <<
                            lastGap*10 << endl;
    }

    if ( least_edge == -1 || least > lastGap*10 )
	return false;

    //        cout << "SUCCESFULL LEFT EXTENSION: " << endl;
    //        cout << "WAS: " <<  path << endl;

    ReadPath tmp;
    tmp.reserve( path.size() + 1 );
    int edge_size = mHBV.EdgeLengthKmers(least_edge);
    int offset = path.getOffset() + edge_size;
    tmp.setOffset(offset);
    tmp.push_back(least_edge);
    std::copy(path.begin(), path.end(), std::back_inserter(tmp) );
    path = tmp;

    //        cout << "NOW: " <<  path << endl;

    if ( mDebug ) PRINT(lastGap);

    return true;
}


bool ExtendReadPath::attemptRightwardExtension( ReadPath& path, basevector const& bases,
						       qualvector const& quals )
{
     ForceAssertEq(quals.size(), bases.size());
    if ( !path.size() ) return false;

    if (mDebug)
	cout << "Rightward from: " << path.getOffset() << ":" << printSeq( path ) << endl;

    // compute how many bases at the end of the read extend past the
    // last edge
    int iLastGap = bases.size();
    iLastGap += path.getOffset();
    for ( auto const& p : path ) iLastGap -= mHBV.EdgeLengthKmers(p);
    iLastGap -= (mHBV.K()-1);

    if (mDebug) PRINT(iLastGap);
    if ( iLastGap < 10 ) return false;

    // cast to be compatible with the rest of the code, now that we
    // know that we're non-negative
    size_t lastGap = iLastGap;

    size_t hbv_edge_id = path.back();
    vec<int>const& to_right = getToRight();
    size_t vright = to_right[hbv_edge_id];

    auto const& edges = mHBV.FromEdgeObj(vright);
    auto const& vdest = mHBV.From(vright);

    int K = mHBV.K();
    vec<int> short_dest;
    vec<bool> ehanging(edges.size(),false);
    vec<bool> elong(edges.size(), false);

    // classify outgoing edges into hanging or not and short or long
    // a short edge is one not capable of "finishing" the read and we ignore
    // short, hanging edges.
    for ( size_t i = 0; i < edges.size(); ++i ) {
	if ( mHBV.FromSize(vdest[i]) == 0 && mHBV.ToSize(vdest[i]) == 1  )
	    ehanging[i] = true;

	if ( mHBV.EdgeObject(edges[i]).size() - (K-1) >= lastGap ) elong[i] = true;

	// push back the destination vertex of a short edge
	if ( !elong[i] && !ehanging[i] )
	    short_dest.push_back( vdest[i] );
    }

    // we only process "short" edges in the following cases:
    // 1. There is only one edge (short or long), or
    //   a. there are no "long" edges mixed in, and
    //   b. all short edges lead to the same vertex, and
    //   c. there is only one outgoing edge from the vertex to which all
    //      short edges lead.
    if ( !edges.solo() ) {
	size_t nlong = Sum(elong);
	if ( short_dest.size() > 0 ) {
	    if ( nlong > 0 ) return false;
	    UniqueSort(short_dest);
	    if ( !short_dest.solo() ) return false;
	    if ( mHBV.FromSize( short_dest.back() ) != 1 ) return false;
	}
    }

    // okay, if we've made it here, we're going to score all
    // edges except hanging ends, unless all edges ARE hanging ends
    int least_edge = -1;
    unsigned least = std::numeric_limits<unsigned>::max();
    vec<int> scores, candidates;
    for ( size_t i = 0; i < edges.size(); ++i ) {
	//            if ( nlong == edges.size() || !ehanging[i] )
	if ( !ehanging[i] || edges.solo() ) {
	    auto score = scoreRightOverlap( bases, quals, lastGap, mHBV.EdgeObject(edges[i]), mHBV.K(),
	            mPenaltyDecay, mMapQ2, mLeftOverPenalty);

	    scores.push_back( score );
	    candidates.push_back( edges[i] );

	    if ( score < least ) {
		least_edge = edges[i];
		least = score;
	    }
	}
    }

    if (mDebug) {
	PRINT(least_edge);
	cout << "original edges: " << printSeq(edges) << endl;
	cout << "considered edges: " << printSeq(candidates) << endl;
	cout << "scores: " << printSeq(scores) << endl;

        cout << "least_edge = " << least_edge << ", score of " << least <<
                  string( ( least <= lastGap*10 ) ? " is" : " is NOT" ) <<
                            " less than lastGap*10 = " <<
                            lastGap*10 << endl;
    }

    if ( least_edge == -1 || least > lastGap*10 )
	return false;

    // this is a relic from when we passed back lastGap...
    size_t edge_size = mHBV.EdgeObject(least_edge).size();
    path.push_back(least_edge);
    if ( edge_size - (K-1) >= lastGap ) {
	/* path.setLastSkip(edge_size - (K-1) - lastGap); */
	lastGap = 0;
    }
    else {
	/* path.setLastSkip(0); */
	lastGap -= (edge_size - (K-1));
    }

    if ( mDebug ) PRINT(lastGap);

    return true;
}
