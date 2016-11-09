///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef DIGRAPH_FROM_WORDS_H
#define DIGRAPH_FROM_WORDS_H

#include "CoreTools.h"
#include "graph/Digraph.h"
#include "paths/long/Fix64_6.h"
#include "paths/long/LongProtoTools.h"

// DigraphFromWords.  Given a unique-sorted collection of words w, form an acyclic 
// digraph by squeezing the words together, whose "all paths" equals w.  For example 
// if the words are {AXB,AYB} then these can be represented as a digraph (without
// squeezing) having six edges:
//
//        --A--> --X--> --B-->
//        --A--> --Y--> --B-->
//
// or after squeezing, four edges (poorly rendered here):
//
//        --A--> --X--> --B-->
//               --Y-->       
//
// The solution here is insanely suboptimal.  See "DFA minimization" in Wikipedia
// for a discussion that may be relevant, and which suggests that there is a 
// canonical solution to the problem.
//
// If w is nonempty then there will be a unique source and a unique sink.
//
// Suppose every entry in w has length at least three.  Then if the first entries
// in each w are equal, they will be presented as a single edge, and if the last
// entries in w are equal, they will be presented as a single edge.  

void DigraphFromWords( const vec< vec<int> >& w, digraphE<int>& D );

// WordsToDigraph( w, D, inv ): given a collection of words w, find the overlaps 
// between the words and then use them to define a digraph D whose edges are labeled
// by words.  For example, given words abc, bcx, bcy, one would get a graph with one 
// edge labelled abc, and two following edges labelled x and y respectively.

void WordsToDigraph( const vec< vec<int> >& w, digraphE< vec<int> >& D,
     vec<int>& inv, const long_logging_control& log_control, 
     const long_logging& logc );

void WordsToDigraphAlt( const vec< vec<int> >& w, const vec<fix64_6>& ww,
     digraphE< vec<int> >& D, vec<int>& inv, 
     const long_logging_control& log_control, const long_logging& logc );

void WordsToDigraphAlt2( const vec< vec<int> >& w, 
     const vec< vec< pair<int,int> > >& over, digraphE< vec<int> >& D,
     vec<int>& inv, const long_logging_control& log_control, 
     const long_logging& logc );

#endif
