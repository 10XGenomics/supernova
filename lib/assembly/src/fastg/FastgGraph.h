///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * FastgGraph.h
 *
 *  Created on: May 21, 2013
 *      Author: blau
 */

#ifndef FASTGGRAPH_H_
#define FASTGGRAPH_H_

#include "feudal/CharString.h"
#include "paths/HyperBasevector.h"
#include "Basevector.h"
#include "graph/Digraph.h"
#include "paths/HyperEfasta.h"
#include <list>

namespace fastg

{
//Write FASTG expression of a hyperbasevector to the specified file
void
WriteFastg(const String&filename, const HyperBasevector& graph);
//Write FASTG expression of a hyper efasta to the specified file
void
WriteFastg(const String&filename, const HyperEfasta& graph);

void FastgHeader(ostream& os);
void FastgFooter(ostream& os);

// a FASTG sequence is essentially a collection of elements,
// each containing a sequence graph and a possible sequence of N's.
// Each element contains enough information to generate a
// (possibly empty) tag plus a canonical sequence

class fastg_element // digraphE does _not_ have virtual destructor and is _dangerous_ to be serving as a base class, since we have many of these here
{
public:
    typedef digraphE<basevector> graph_t;
    static const String::size_type nCanonicalN = 20; //if the average gap is non-zero, put in this number of N as a canonical sequence
    static const int64_t invalidGap = std::numeric_limits<int64_t>::max(); //if the average gap is non-zero, put in this number of N as a canonical sequence

    fastg_element() : graph(), avgGap(invalidGap), minGap(invalidGap), maxGap( invalidGap) { } ;
    ~fastg_element() { } ;
    fastg_element(const fastg_element&in) : graph(in.graph), avgGap(in.avgGap), minGap(in.minGap), maxGap( in.maxGap) { } ;
    fastg_element& operator=(const fastg_element&in) { graph = in.graph; avgGap = in.avgGap; minGap = in.minGap; maxGap = in.maxGap; return *this; } ;
    fastg_element& operator=(const basevector&in) { graph.Clear(); addSequence(in); avgGap = invalidGap; minGap = invalidGap; maxGap = invalidGap; return *this; } ;
    void clear() { graph.Clear(); avgGap = invalidGap; minGap = invalidGap; maxGap = invalidGap; } ;

    const graph_t& getGraph() const { return graph; } ;
    graph_t& getGraph() { return graph; } ;

    // add a sequence as a disconnected component in the graph
    void addSequence(const basevector&in);
    void appendSequence(const basevector&in);
    // a non-zero avg specifies that this element is a gap
    void setGap(int64_t avg = invalidGap, int64_t minimum = invalidGap, int64_t maximum = invalidGap);
    std::tuple<int64_t,int64_t,int64_t> getGap()const{return std::make_tuple(avgGap,minGap,maxGap);};


    bool isEqual(const fastg_element&other)const;
    bool isSequence()const{ return avgGap==invalidGap && graph.EdgeObjectCount()<2 && graph.N()==graph.EdgeObjectCount()*2;};
    bool isEmpty()const{ return isSequence() && (graph.EdgeObjectCount()==0 || graph.EdgeObject(0).size()==0 );};

    //return a string
    String ToString() const;
    //append the full fastg entry to the string
    void AppendTo(String& out) const { AppendTagTo(out); AppendCanonicalTo(out); } ;
    //append the FASTG tag to a string
    void AppendTagTo(String& out) const;
    //append the FASTG canonical sequence to a string
    void AppendCanonicalTo(String& out) const;

    //the following considers all disconnected components
    size_t nFreeL() const;   // number of bases that are free to be moved from the left
    size_t nFreeR() const;  // number of bases that are free to be moved from the right
    basevector TrimL(size_t n=std::numeric_limits<size_t>::max()); // remove up to n bases from the left and return such bases
    basevector TrimR(size_t n=std::numeric_limits<size_t>::max()); // remove up to n bases from the right and return such bases
private:
    graph_t graph;
    int64_t avgGap;
    int64_t minGap;
    int64_t maxGap;
    void AppendGraphTo(String& out) const;
    vec<int> getCanonicalDigraphPath() const;
    String getFastgProperties(const vec<int>&path=vec<int>())const;
};

bool operator==(const fastg_element&left,const fastg_element&right);
bool operator!=(const fastg_element&left,const fastg_element&right);

// a fastg edge contains a sequence of fastg entries, this is useful for organization purposes
class fastg_sequence
{
public:
    fastg_sequence() : entries() { } ;
    fastg_sequence(const fastg_sequence& in):entries(in.entries){};
    explicit fastg_sequence(const basevector&in) : entries(1) { entries.back().addSequence(in); } ;
    explicit fastg_sequence(const efasta&in) {assign(in);};
    fastg_sequence& operator=(const basevector&in) { entries.resize(1); entries.back() = in; return *this; } ;

    fastg_sequence& operator=(const efasta&in) { return assign(in);};
    fastg_sequence& assign(const efasta&in);
    ~fastg_sequence() { } ;

    const std::list<fastg_element>& getEntries()const{return entries;};
    std::list<fastg_element>& getEntries(){return entries;};

    bool isEqual(const fastg_sequence&)const;

    //return a string
    String ToString() const;
    //append the full fastg entry to the string
    void AppendTo(String&out) const;
    //append the FASTG canonical sequence to a string
    void AppendCanonicalTo(String&out) const;

    void Compact();
    //the following considers all disconnected components
    size_t nFreeL() const;   // number of bases that are free to be moved from the left
    size_t nFreeR() const;  // number of bases that are free to be moved from the right
    basevector TrimL(size_t n=std::numeric_limits<size_t>::max()); // remove up to n bases from the left and return such bases
    basevector TrimR(size_t n=std::numeric_limits<size_t>::max()); // remove up to n bases from the right and return such bases
private:
    std::list<fastg_element> entries;
};
bool operator==(const fastg_sequence&left,const fastg_sequence&right);
bool operator!=(const fastg_sequence&left,const fastg_sequence&right);

// a fastg vertex stores the effective K value in the hyper kmer context; that is
// the last K-1 bases of an incoming edge overlaps with the first K-1 bases of an out-going edge
typedef size_t fastg_k;

//this is a wrapper adding axilary function to digraphVE, no additional non-static data members
class fastg_graph : public digraphVE<fastg_k, fastg_sequence>
{
public:
    typedef digraphVE<fastg_k, fastg_sequence> base_t;
    explicit fastg_graph(const HyperEfasta& in);
    explicit fastg_graph(const HyperBasevector& in);
    //reduce K->1, forming sequence graph
    void MakeK1();
    // simplify graph by merging multiple sequence edges into a fastg sequence
    void Compact();
    //interpret current graph as K=1 and write to a fastg file

    void Write(const String& filename) const;
    void Write(ostream& os) const;

    void WriteHeader(ostream& os) const;
    void WriteFooter(ostream& os) const;
};

class fastg_vertex_graph : public digraphVE<fastg_sequence, fastg_k>
{
public:
    typedef digraphVE<fastg_sequence, fastg_k> base_t;
    template<typename E> explicit fastg_vertex_graph( const digraphE<E>& graph, size_t k);
//    fastg_vertex_graph(const HyperEfasta& in);
//    fastg_vertex_graph(const HyperBasevector& in);
    //reduce K->1, forming sequence graph
    void MakeK1();
    // simplify graph by merging multiple sequence edges into a fastg sequence
    void Compact();
    void RemoveDuplicatedVertices();
    //interpret current graph as K=1 and write to a fastg file

    void Write(const String& filename) const;
    void Write(ostream& os) const;

    void WriteHeader(ostream& os) const;
    void WriteFooter(ostream& os) const;
private:
    void ConnectedVertices(std::set<int>& left, std::set<int>& edges, std::set<int>&right,int root, bool bOnRight);
};

/*
 template<typename edge_t>
 void WriteFastg( ostream&os, const digraphE<edge_t>&graph){
 for(int vv=0;vv<graph.N();++vv){
 const int nEdges = graph.From(vv).size();
 for(int ee=0;ee<nEdges;++ee){
 auto edge_index = graph.EdgeObjectIndexByIndexFrom(vv,ee);
 os << ">" << edge_index;
 auto ww = graph.From(vv)[ee];
 if( graph.From(ww).size() > 0){
 os << ":" << graph.EdgeObjectIndexByIndexFrom(ww,0);
 for(size_t ff=1;ff<graph.From(ww).size();++ff){
 os<<"," << graph.EdgeObjectIndexByIndexFrom(ww,ff);
 }
 }
 os << ";\n";
 os << graph.EdgeObjectByIndexFrom(vv,ee).ToString() << "\n";
 }
 }
 }
 */

} /* namespace fastg */

#endif /* FASTGGRAPH_H_ */
