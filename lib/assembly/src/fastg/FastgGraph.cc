///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * FastgGraph.cc
 *
 *  Created on: May 21, 2013
 *      Author: blau
 */

#include "fastg/FastgGraph.h"

#include <fstream>
namespace fastg
{
void
FastgHeader(ostream&os) {
    os << "#FASTG:begin;\n" << "#FASTG:version=1.1;\n";
}
void
FastgFooter(ostream&os) {
    os << "#FASTG:end;\n";
}

void
fastg_element::addSequence(const basevector& in) {
    int old_size = graph.N();
    graph.AddVertices(2);
    graph.AddEdge(old_size, old_size + 1, in);
}

void
fastg_element::appendSequence(const basevector& in) {
    if( graph.EdgeObjectCount()==0){
        addSequence(in);
    }
    else{
        for(int vv=0;vv<graph.N();++vv){
            if( graph.Sink(vv)){
                for(int ee=0;ee<graph.ToSize(vv);++ee){
                    graph.EdgeObjectByIndexToMutable(vv,ee).append(in.begin(),in.end());
                }
            }
        }
    }
}

void
fastg_element::setGap(int64_t avg, int64_t minimum, int64_t maximum) {
    avgGap = avg;
    minGap = minimum;
    maxGap = maximum;
}
String
fastg_element::ToString() const {
    String out;
    AppendTo(out);
    return out;
}
//append the FASTG tag to a string
void
fastg_element::AppendTagTo(String& out) const {
    if (avgGap != invalidGap) {
        out.append(
                "[" + ::ToString(nCanonicalN) + ":gap:size=("
                        + ::ToString(avgGap));
        if (minGap != invalidGap && maxGap != invalidGap ) {
            out.append("," + ::ToString(minGap) + ".." + ::ToString(maxGap));
        }
        out.append(")");
        if (graph.N() != 0) {
            auto path = getCanonicalDigraphPath();
            String property = getFastgProperties(path);
            if( property.size() >0){ out.append(","+property); }
            out.append("|");
            AppendGraphTo(out);
        }
        out.append("]");
    }
    else if (graph.EdgeObjectCount() < 2) {
        return;
    }
    else if (graph.N() == 2 * graph.EdgeObjectCount()) {
        auto canonical_size = graph.EdgeObject(0).size();
        out.append("[" + ::ToString(canonical_size) + ":alt|");
        out.append(graph.EdgeObject(0).ToString());
        for (decltype(graph.EdgeObjectCount()) ii = 1;
                ii < graph.EdgeObjectCount(); ++ii) {
            out.append(",");
            out.append(graph.EdgeObject(ii).ToString());
        }
        out.append("]");
    }
    else {
        auto path = getCanonicalDigraphPath();
        size_t canonical_size = 0;
        for (auto entry : path) {
            canonical_size += graph.EdgeObject(entry).size();
        }
        out.append("[" + ::ToString(canonical_size) + ":digraph");
        String property = getFastgProperties(path);
        if( property.size() >0){ out.append(":"+property); }
        out.append("|");
        AppendGraphTo(out);
        out.append("]");
    }
}
;

String
fastg_element::getFastgProperties(const vec<int>&path)const{
    String property;
    if( path.size() > 0){
        property.append("path=(");
        property.append(::ToString(path[0]));
        for(size_t ii=1;ii<path.size();++ii){
            property.append(","+::ToString(path[ii]));
        }
        property.append("),begin=(");
        property.append(::ToString(path.front()));
        property.append("),end=(");
        property.append(::ToString(path.back()));
        property.append(")");
    }
    return property;
}
;
vec<int>
fastg_element::getCanonicalDigraphPath() const {
    int start = -1;
    int end = -1;
    int last = -1;
    for (decltype(graph.N()) ii = 0; ii < graph.N(); ++ii) {
        if (graph.Source(ii)) {
            start = ii;
            break;
        }
    }
    if (start < 0)
        start = 0;
    std::deque<int> todo;
    vec<int> parent(graph.N(), -1);
    todo.push_back(start);
    while (!todo.empty() && end < 0) {
        auto src = todo.front();
        last = src;
        todo.pop_front();
        for (auto entry : graph.From(src)) {
            if (parent[entry] == -1) {
                parent[entry] = src;
                todo.push_back(entry);
                if (graph.Sink(entry)) {
                    end = entry;
                    break;
                }
            }
        }
    }
    if (end < 0)
        end = last;
    vec<int> path;
    for (auto vv = end; vv != start; vv = parent[vv]) {
        path.push_back(vv);
    }
    path.push_back(start);
    vec<int> out;
    for (int ii = path.isize() - 1; ii > 0; --ii) {
        auto src = path[ii];
        auto dst = path[ii - 1];
        for (int jj = 0; jj < graph.FromSize(src); ++jj) {
            if (graph.From(src)[jj] == dst) {
                out.push_back(graph.EdgeObjectIndexByIndexFrom(src, jj));
                continue;
            }
        }
    }
    return out;
}
void
fastg_element::AppendGraphTo(String& out) const {
    for (int vv = 0; vv < graph.N(); ++vv) {
        const int nEdges = graph.From(vv).size();
        for (int ee = 0; ee < nEdges; ++ee) {
            auto edge_index = graph.EdgeObjectIndexByIndexFrom(vv, ee);
            out.append(">" + ::ToString(edge_index));
            auto ww = graph.From(vv)[ee];
            if (graph.From(ww).size() > 0) {
                out.append(
                        ":"
                                + ::ToString(
                                        graph.EdgeObjectIndexByIndexFrom(ww,
                                                0)));
                for (size_t ff = 1; ff < graph.From(ww).size(); ++ff) {
                    out.append(
                            ","
                                    + ::ToString(
                                            graph.EdgeObjectIndexByIndexFrom(ww,
                                                    ff)));
                }
            }
            out.append(";");
            out.append(graph.EdgeObjectByIndexFrom(vv, ee).ToString());
        }
    }
}
//append the FASTG canonical sequence to a string
void
fastg_element::AppendCanonicalTo(String& out) const {
    if (avgGap != invalidGap) {
        out.append(nCanonicalN, 'N');
    }
    else if (graph.EdgeObjectCount() == 0) {
      return;
    }
    else if (graph.EdgeObjectCount() == 1) {
        out += graph.EdgeObject(0).ToString();
    }
    else if (graph.N() == 2 * graph.EdgeObjectCount()) {
        out.append(graph.EdgeObject(0).ToString());
    }
    else {
        const auto canonical_path = getCanonicalDigraphPath();
        for (auto entry : canonical_path) {
            out.append(graph.EdgeObject(entry).ToString());
        }
    }
}
;

size_t
fastg_element::nFreeR() const {
    size_t max_overlap = std::numeric_limits<size_t>::max();
    vec< std::pair<int,int> > edge_indices;

    for (int vv = 0; vv < graph.N(); ++vv) {
        if (graph.Sink(vv)) {
            for (int ee = 0; ee < graph.ToSize(vv); ++ee) {
                edge_indices.push_back( std::make_pair(vv,ee));
                max_overlap = std::min(max_overlap,
                        (size_t) graph.EdgeObjectByIndexTo(vv, ee).size());
            }
        }
    }
    if( max_overlap == std::numeric_limits<size_t>::max()){
        max_overlap=0;
    }
    else{
        for(size_t ii = 0 ;ii<edge_indices.size();++ii){
            const auto& edge_i=graph.EdgeObjectByIndexTo(edge_indices[ii].first,edge_indices[ii].second);
            for(size_t jj = ii+1 ;jj<edge_indices.size();++jj){
                const auto& edge_j=graph.EdgeObjectByIndexTo(edge_indices[jj].first,edge_indices[jj].second);
                size_t overlap=0;
                for(overlap=0;overlap<max_overlap && edge_i[edge_i.size()-1-overlap]==edge_j[edge_j.size()-1-overlap];++overlap){}
                max_overlap=std::min(max_overlap,overlap);
            }
        }
    }
    return max_overlap;
}
size_t
fastg_element::nFreeL() const {
    size_t max_overlap = std::numeric_limits<size_t>::max();
    vec< std::pair<int,int> > edge_indices;

    for (int vv = 0; vv < graph.N(); ++vv) {
        if (graph.Source(vv)) {
            for (int ee = 0; ee < graph.FromSize(vv); ++ee) {
                edge_indices.push_back( std::make_pair(vv,ee));
                max_overlap = std::min(max_overlap,
                        (size_t) graph.EdgeObjectByIndexFrom(vv, ee).size());
            }
        }
    }
    if( max_overlap == std::numeric_limits<size_t>::max()){
        max_overlap=0;
    }
    else{
        for(size_t ii = 0 ;ii<edge_indices.size();++ii){
            const auto& edge_i=graph.EdgeObjectByIndexFrom(edge_indices[ii].first,edge_indices[ii].second);
            for(size_t jj = ii+1 ;jj<edge_indices.size();++jj){
                const auto& edge_j=graph.EdgeObjectByIndexFrom(edge_indices[jj].first,edge_indices[jj].second);
                size_t overlap=0;
                for(overlap=0;overlap<max_overlap && edge_i[overlap]==edge_j[overlap];++overlap){}
                max_overlap=std::min(max_overlap,overlap);
            }
        }
    }
    return max_overlap;
}

basevector
fastg_element::TrimR(size_t n) {
    size_t to_trim = min(nFreeR(), n);
    basevector out;
    for (int vv = 0; vv < graph.N(); ++vv) {
        if (graph.Sink(vv)) {
            for (int ee = 0; ee < graph.ToSize(vv); ++ee) {
                basevector& edge = graph.EdgeObjectByIndexToMutable(vv, ee);
                if (to_trim > edge.size()) {
                    std::cout
                            << "WARNING: fastg_element::TrimR ordered to trim beyond edge.size()"
                            << std::endl;
                    continue;
                }
                if (out.size() < to_trim) {
                    out.SetToSubOf(edge, edge.size() - to_trim, to_trim);
                }
                edge.resize(edge.size() - to_trim);
            }
        }
    }
    return out;
}



basevector
fastg_element::TrimL(size_t n) {
    size_t to_trim = min(nFreeL(), n);
    basevector out;
    for (int vv = 0; vv < graph.N(); ++vv) {
        if (graph.Source(vv)) {
            for (int ee = 0; ee < graph.FromSize(vv); ++ee) {
                basevector& edge = graph.EdgeObjectByIndexFromMutable(vv, ee);
                if (to_trim > edge.size()) {
                    std::cout
                            << "WARNING: fastg_element::TrimL ordered to trim beyond edge.size()"
                            << std::endl;
                    continue;
                }
                if (out.size() < to_trim) {
                    out.SetToSubOf(edge, 0, to_trim);
                }
                edge.SetToSubOf(edge, to_trim, edge.size()-to_trim);
            }
        }
    }
    return out;
}

bool
fastg_element::isEqual(const fastg_element&other)const{
    return getGap()==other.getGap() && EqualExceptEdgeObjectOrder(graph,other.getGraph());
}
bool
operator==(const fastg_element&left,const fastg_element&right){
    return left.isEqual(right);
}
bool
operator!=(const fastg_element&left,const fastg_element&right){
    return !(left==right);
}

fastg_sequence&
fastg_sequence::assign(const efasta&in) {
    entries.clear();

    for (efasta::size_type begin = 0;
            begin != efasta::npos && begin < in.size();) {
        auto end = in.find_first_not_of("ACGTacgt", begin);
        if (end == efasta::npos) {
            entries.push_back(fastg_element());
            entries.back().addSequence(basevector(in.substr(begin)));
            begin = end;
        }
        else {
            if (begin < end) {
                entries.push_back(fastg_element());
                entries.back().addSequence(
                        basevector(in.substr(begin, end - begin)));
            }
            begin = end;
            if (in[begin] != '{') {
                std::cout << "WARNING: fastg_element::operator= sees character "
                        << in[begin]
                        << " outside of controlled efasta region. Terminating conversion."
                        << std::endl;
                return *this;
            }
            else {
                end = in.find_first_of("}", begin + 1);
                if (end == efasta::npos) {
                    std::cout
                            << "WARNING: fastg_element::operator= can\'t locate matching }. Terminating conversion."
                            << std::endl;
                    return *this;
                }
                entries.push_back(fastg_element());
                SpecialEfastaTag2Graph(entries.back().getGraph(),
                        in.substr(begin, end + 1 - begin));
//                if( entries.back().getGraph().EdgeObjectCount() > 1 && entries.back().getGraph().EdgeObjectCount()*2 != entries.back().getGraph().N()){
//                    entries.back().setGap(10,-2,100);
//                }
                begin = end + 1;
            }
        }
    }
    return *this;
}

void
fastg_sequence::Compact(){
    if( entries.size() < 2 ){
        return;
    }
    for( auto itr=entries.begin(); itr!=entries.end() && std::next(itr) != entries.end();){
        auto next=std::next(itr);
        if( (*itr).isSequence() && (*next).isSequence()){
            if( (*next).getGraph().EdgeObjectCount() == 1){
                (*itr).appendSequence((*next).getGraph().EdgeObject(0));
            }
            entries.erase(next);
        }
        else{
            ++itr;
        }
    }

}

String
fastg_sequence::ToString() const {
    String out;
    AppendTo(out);
    return out;
}

void
fastg_sequence::AppendTo(String&out) const {
    for (const auto& entry : entries) {
        entry.AppendTo(out);
    }
}

void
fastg_sequence::AppendCanonicalTo(String&out) const {
    for (const auto& entry : entries) {
        entry.AppendCanonicalTo(out);
    }
}

/*
 size_t fastg_sequence::nFreeL()const{
 return entries.front().nFreeL();
 }
 */
size_t
fastg_sequence::nFreeL() const {
    return entries.front().nFreeL();
}
size_t
fastg_sequence::nFreeR() const {
    return entries.back().nFreeR();
}
 basevector fastg_sequence::TrimL(size_t n){
 return entries.front().TrimL(n);
 }

basevector
fastg_sequence::TrimR(size_t n) {
    return entries.back().TrimR(n);
}
bool
fastg_sequence::isEqual(const fastg_sequence&that)const{
    auto itr_this = entries.begin();
    auto itr_that = that.entries.begin();
    auto end_this = entries.end();
    auto end_that = that.entries.end();

    while( itr_this!=end_this && itr_that!=end_that){
        if( (*itr_this).isEmpty() ){
            ++itr_this;
        }
        else if( (*itr_that).isEmpty() ){
            ++itr_that;
        }
        else{
            if( *itr_this!=*itr_that ) return false;
            ++itr_this;
            ++itr_that;
        }
    }
    for(;itr_this!=end_this;++itr_this){ if( !(*itr_this).isEmpty() ) return false; }
    for(;itr_that!=end_that;++itr_that){ if( !(*itr_that).isEmpty() ) return false; }
    return true;
}
bool
operator==(const fastg_sequence&left,const fastg_sequence&right){
    return left.isEqual(right);
}
bool
operator!=(const fastg_sequence&left,const fastg_sequence&right){
    return !(left==right);
}

fastg_graph::fastg_graph(const HyperBasevector&in) {
    const size_t nEdges = in.EdgeObjectCount();
    vec<fastg_sequence> edges(nEdges);
    for (size_t ii = 0; ii < nEdges; ++ii) {
        edges[ii] = in.EdgeObject(ii);
    }
    base_t::Initialize(in.From(), in.To(), vec<fastg_k>(in.N(), in.K()),
            edges, in.ToEdgeObj(), in.FromEdgeObj());
    MakeK1();
}

fastg_graph::fastg_graph(const HyperEfasta&in) {
    const size_t nEdges = in.EdgeObjectCount();
    vec<fastg_sequence> edges(nEdges);
    for (size_t ii = 0; ii < nEdges; ++ii) {
        edges[ii] = in.EdgeObject(ii);
    }

    base_t::Initialize(in.From(), in.To(), vec<fastg_k>(in.N(), in.K()),
            edges, in.ToEdgeObj(), in.FromEdgeObj());
    MakeK1();
}
void
fastg_graph::Write(const String&filename) const {
    ofstream out(filename.c_str(), ios::out);
    Write(out);
    out.close();
}

void
fastg_graph::MakeK1() {
    basevector tmp;
    for (int vv = 0; vv < N(); ++vv) {
        if (Vert(vv) < 2)
            continue;
        const auto k_minus_1 = Vert(vv) - 1;
        size_t to_trim = k_minus_1;
        if (FromSize(vv) > 0 && ToSize(vv) > 0) {
            for (int ii = 0; ii < ToSize(vv); ++ii) {
                to_trim = std::min(to_trim,EdgeObjectByIndexToMutable(vv, ii).nFreeR());
            }
            for (int ii = 0; ii < ToSize(vv); ++ii) {
                tmp = EdgeObjectByIndexToMutable(vv, ii).TrimR(to_trim);
                if (tmp.size() != to_trim) {
                    std::cout << "WARNING: fastg_graph::MakeK1 trimmed "
                            << tmp.size() << " instead of " << k_minus_1
                            << std::endl;
                }
            }
        }
        VertMutable(vv) -= to_trim;
    }
    for (int vv = 0; vv < N(); ++vv) {
        if (Vert(vv) < 2)
            continue;
        const auto k_minus_1 = Vert(vv) - 1;
        size_t to_trim = k_minus_1;
        if (FromSize(vv) > 0 && ToSize(vv) > 0) {
            for (int ii = 0; ii < FromSize(vv); ++ii) {
                to_trim = std::min(to_trim,EdgeObjectByIndexFromMutable(vv, ii).nFreeL());
            }
            for (int ii = 0; ii < FromSize(vv); ++ii) {
                tmp = EdgeObjectByIndexFromMutable(vv, ii).TrimL(to_trim);
                if (tmp.size() != to_trim) {
                    std::cout << "WARNING: fastg_graph::MakeK1 trimmed "
                            << tmp.size() << " instead of " << k_minus_1
                            << std::endl;
                }
            }
        }
        VertMutable(vv) -= to_trim;
    }
    for (int vv = 0; vv < N(); ++vv) {
        if( Vert(vv)!=1){
            std::cout << "WARNING: fastg_graph::MakeK1: vertex " << vv << " has K="<<Vert(vv) << std::endl;;
        }
    }
}
void
fastg_graph::Compact(){
    std::cout << "compacting" << std::endl;
    for( auto& entry: EdgesMutable()){
        entry.Compact();
    }
}

void
fastg_graph::Write(ostream&os) const {
    FastgHeader(os);
    vec< std::pair<int,int> > edge_indices(EdgeObjectCount(),std::make_pair(-1,-1));
    for (int vv = 0; vv < N(); ++vv) {
        const int nEdges = From(vv).size();
        for (int ee = 0; ee < nEdges; ++ee) {
            auto edge_index = EdgeObjectIndexByIndexFrom(vv, ee);
            auto& entry = edge_indices[edge_index];
            if( entry.first!=-1 || entry.second!=-1) std::cout<< "WARNING: fastg_graph::Write - edge discovered more than once"<< std::endl;
            entry.first=vv;
            entry.second=ee;
        }
    }
    String sBuffer;
    for(decltype(EdgeObjectCount()) ii=0;ii<EdgeObjectCount();++ii){
        auto vv = edge_indices[ii].first;
        auto ee   = edge_indices[ii].second;
        if( vv==-1 || ee==-1 || ii!=EdgeObjectIndexByIndexFrom(vv,ee)){
            std::cout<< "WARNING: fastg_graph::Write - edge logging not consistent. Edge " << ii << " not written."<< std::endl;
            continue;
        }
        os << ">" << ii;
        auto ww = From(vv)[ee];
        if (From(ww).size() > 0) {
            os << ":" << EdgeObjectIndexByIndexFrom(ww, 0);
            for (size_t ff = 1; ff < From(ww).size(); ++ff) {
                os << "," << EdgeObjectIndexByIndexFrom(ww, ff);
            }
        }
        os << ";\n";
        sBuffer.clear();
        EdgeObjectByIndexFrom(vv, ee).AppendTo(sBuffer);
        for (size_t ii = 0; ii < sBuffer.size(); ++ii) {
            os << sBuffer[ii];
            if (ii % 80 == 79) {
                os << '\n';
            }
        }
        if (sBuffer.size() % 80 != 0) os << '\n';
    }
    /*
    size_t nDup=0;
    for(decltype(EdgeObjectCount()) ii=0;ii<EdgeObjectCount();++ii){
        for(decltype(EdgeObjectCount()) jj=ii+1;jj<EdgeObjectCount();++jj){
            if( EdgeObject(ii) == EdgeObject(jj)){
                std::cout << ii << " " << jj << std::endl;
                ++nDup;
            }
        }
    }
    std::cout << "dups: " << nDup << std::endl;
    */

/*
    for (int vv = 0; vv < N(); ++vv) {
        const int nEdges = From(vv).size();
        for (int ee = 0; ee < nEdges; ++ee) {
            auto edge_index = EdgeObjectIndexByIndexFrom(vv, ee);
            os << ">" << edge_index;
            auto ww = From(vv)[ee];
            if (From(ww).size() > 0) {
                os << ":" << EdgeObjectIndexByIndexFrom(ww, 0);
                for (size_t ff = 1; ff < From(ww).size(); ++ff) {
                    os << "," << EdgeObjectIndexByIndexFrom(ww, ff);
                }
            }
            os << ";\n";
            sBuffer.clear();
            EdgeObjectByIndexFrom(vv, ee).AppendTo(sBuffer);
            for (size_t ii = 0; ii < sBuffer.size(); ++ii) {
                os << sBuffer[ii];
                if (ii % 80 == 79) {
                    os << '\n';
                }
            }
            if (sBuffer.size() % 80 != 0)
                os << '\n';
        }
    }
*/
    FastgFooter(os);
}

template<typename F>
fastg_vertex_graph::fastg_vertex_graph(const digraphE<F>&in, size_t k) {
    const size_t nVertices = in.EdgeObjectCount();

    vec<fastg_sequence> vertices(nVertices);
    for (size_t ii = 0; ii < nVertices; ++ii) {
        AddVertex(fastg_sequence(in.EdgeObject(ii)));
    }
    for(decltype(in.N()) vv=0;vv<in.N();++vv){
        for( decltype(in.FromSize(vv)) ee=0;ee<in.FromSize(vv);++ee){
            auto edge_index = in.EdgeObjectIndexByIndexFrom(vv,ee);
            auto ww = in.From(vv)[ee];
            for( decltype(in.FromSize(ww)) ff=0;ff<in.FromSize(ww);++ff){
                auto other_edge_index = in.EdgeObjectIndexByIndexFrom(ww,ff);
                AddEdge(edge_index,other_edge_index,k);
            }
        }
    }
    MakeK1();
}

void
fastg_vertex_graph::Write(const String&filename) const {
    ofstream out(filename.c_str(), ios::out);
    Write(out);
    out.close();
}
void
fastg_vertex_graph::ConnectedVertices(std::set<int>& left, std::set<int>& edges, std::set<int>&right,int root, bool bOnRight){
    left.clear();
    edges.clear();
    right.clear();

    std::deque<int> todoL,todoR;

    if(bOnRight){
        right.insert(root);
        todoR.push_back(root);
    }
    else{
        left.insert(root);
        todoL.push_back(root);
    }

    for(;!todoL.empty()||!todoR.empty();){
        while( !todoL.empty()){
            auto vv = todoL.front();
            todoL.pop_front();
            for(int ee=0;ee<FromSize(vv);++ee){
                auto ww = From(vv)[ee];
                edges.insert(EdgeObjectIndexByIndexFrom(vv,ee));
                if( right.find(ww)==right.end()){
                    todoR.push_back(ww);
                    right.insert(ww);
                }
            }
        }
        while( !todoR.empty()){
            auto ww = todoR.front();
            todoR.pop_front();
            for(int ee=0;ee<ToSize(ww);++ee){
                auto vv = To(ww)[ee];
                edges.insert(EdgeObjectIndexByIndexTo(ww,ee));
                if( left.find(vv)==left.end()){
                    todoL.push_back(vv);
                    left.insert(vv);
                }
            }
        }
    }
}


void
fastg_vertex_graph::MakeK1() {
    basevector tmp;
    std::set<int> left,edges,right;
    for( int vv=0;vv<N();++vv){
        ConnectedVertices(left,edges,right,vv,0);
        if( left.size()>0 && right.size()>0){
            size_t lock=std::numeric_limits<size_t>::max();
            for( auto entry: edges){
                if( lock == std::numeric_limits<size_t>::max() ) lock=EdgeObject(entry);
                if( lock!= EdgeObject(entry)) std::cout << "WARNING: fastg_vertex_graph::MakeK1() - inconsistent K detected" << std::endl;
            }
            if( lock >1 ){
                const auto k_minus_1 = lock-1;
                size_t to_trim = k_minus_1;
                for(auto entry:left){
                    to_trim = std::min(to_trim,Vert(entry).nFreeR());
                }
                for(auto entry:left){
                    tmp=VertMutable(entry).TrimR(to_trim);
                    if (tmp.size() != to_trim) {
                        std::cout << "WARNING: fastg_graph::MakeK1 trimmed "
                                << tmp.size() << " instead of " << k_minus_1
                                << std::endl;
                    }
                }
                for( auto entry: edges){
                    EdgeObjectMutable(entry)-=to_trim;
                }
            }
        }
    }
    for( int vv=0;vv<N();++vv){
        ConnectedVertices(left,edges,right,vv,1);
        if( left.size()>0 && right.size()>0){
            size_t lock=std::numeric_limits<size_t>::max();
            for( auto entry: edges){
                if( lock == std::numeric_limits<size_t>::max() ) lock=EdgeObject(entry);
                if( lock!= EdgeObject(entry)) std::cout << "WARNING: fastg_vertex_graph::MakeK1() - inconsistent K detected" << std::endl;
            }
            if( lock >1 ){
                const auto k_minus_1 = lock-1;
                size_t to_trim = k_minus_1;
                for(auto entry:right){
                    to_trim = std::min(to_trim,Vert(entry).nFreeL());
                }
                for(auto entry:right){
                    tmp=VertMutable(entry).TrimL(to_trim);
                    if (tmp.size() != to_trim) {
                        std::cout << "WARNING: fastg_graph::MakeK1 trimmed "
                                << tmp.size() << " instead of " << k_minus_1
                                << std::endl;
                    }
                }
                for( auto entry: edges){
                    EdgeObjectMutable(entry)-=to_trim;
                }
            }
        }
    }
    for (int ee = 0; ee < EdgeObjectCount(); ++ee) {
        if( EdgeObject(ee)!=1){
            std::cout << "WARNING: fastg_graph::MakeK1: vertex " << ee << " has K="<<EdgeObject(ee) << std::endl;;
        }
    }
}

void
fastg_vertex_graph::Compact(){
    RemoveDuplicatedVertices();
    /*
    for( auto& entry: VertMutable()){
        entry.Compact();
    }
    */
}
void
fastg_vertex_graph::RemoveDuplicatedVertices(){
    std::set< decltype(N())> to_delete;
    for(decltype(N()) vv=0; vv<N(); ++vv){
        if( to_delete.find(vv)==to_delete.end()){
            const auto& from_list = From(vv);
            for( size_t ee=0;ee<from_list.size();++ee){
                const auto& ww = from_list[ee];
                if( vv==ww || to_delete.find(ww)!=to_delete.end()) continue;
                for( size_t ff=ee+1;ff<from_list.size();++ff){
                    const auto xx = from_list[ff];
                    if( vv==xx || ww==xx || to_delete.find(xx)!=to_delete.end()) continue;

                    if( Vert(ww)==Vert(xx)){
                        TransferEdges(xx,ww);
                        to_delete.insert(xx);
                    }
                }
            }

        }
    }
    for(decltype(N()) ww=0; ww<N(); ++ww){
        if( to_delete.find(ww)==to_delete.end()){
            const auto& to_list = To(ww);
            for( size_t ee=0;ee<to_list.size();++ee){
                const auto& vv = to_list[ee];
                if( vv==ww || to_delete.find(vv)!=to_delete.end()) continue;
                for( size_t ff=ee+1;ff<to_list.size();++ff){
                    const auto uu = to_list[ff];
                    if( vv==uu || ww==uu || to_delete.find(uu)!=to_delete.end()) continue;

                    if( Vert(ww)==Vert(uu)){
                        TransferEdges(uu,vv);
                        to_delete.insert(uu);
                    }
                }
            }

        }
    }
    RemoveDuplicateEdges( );
    vec< decltype(N()) > vDel;
    vDel.reserve(to_delete.size());
    for( auto entry: to_delete){
        vDel.push_back(entry);
    }
    RemoveVertices(vDel);
}


void
fastg_vertex_graph::Write(ostream&os) const {
    FastgHeader(os);
    String sBuffer;
    for(decltype(EdgeObjectCount()) vv=0;vv<N();++vv){
        os << ">" << vv;
        if (From(vv).size() > 0) {
            os << ":" << From(vv)[0];
            for (size_t wi = 1; wi < From(vv).size(); ++wi) {
                os << "," << From(vv)[wi];
            }
        }
        os << ";\n";
        sBuffer.clear();
        Vert(vv).AppendTo(sBuffer);
        for (size_t ii = 0; ii < sBuffer.size(); ++ii) {
            os << sBuffer[ii];
            if (ii % 80 == 79) {
                os << '\n';
            }
        }
        if (sBuffer.size() % 80 != 0) os << '\n';
    }
    FastgFooter(os);
}
void
WriteFastg(const String&filename, const HyperBasevector& graph) {
//    fastg::fastg_graph fgg(graph);
//    fgg.Write(filename);
    {
        fastg::fastg_vertex_graph fgg(graph,graph.K());
        fgg.Compact();
        fgg.Write(filename);
        vec< vec<String>> edge_labels(fgg.N());
        vec< String> vertex_labels(fgg.N());
        for(size_t ii=0;ii<edge_labels.size();++ii){
            vertex_labels[ii]=ToString(ii);
            int nn = fgg.FromSize(ii);
            edge_labels[ii].resize(nn);
            for(int jj = 0 ; jj < nn ; ++jj){
                edge_labels[ii][jj] = ToString(fgg.EdgeObjectIndexByIndexFrom(ii,jj));
            }
        }
        Ofstream( dot_out, filename+".dot" );
        fgg.DOT_vl( dot_out, vertex_labels,
             "",
             vec< vec<String> >( ),
             vec<String>( ),
             vec<vec<String>>(),
             vec<String>( ) ) ;
        dot_out.close();
    }
}
void
WriteFastg(const String&filename, const HyperEfasta& graph) {
    {
        fastg::fastg_vertex_graph fgg(graph,graph.K());
        fgg.Compact();
        fgg.Write(filename);
        vec< vec<String>> edge_labels(fgg.N());
        vec< String> vertex_labels(fgg.N());
        for(size_t ii=0;ii<edge_labels.size();++ii){
            vertex_labels[ii]=ToString(ii);
            int nn = fgg.FromSize(ii);
            edge_labels[ii].resize(nn);
            for(int jj = 0 ; jj < nn ; ++jj){
                edge_labels[ii][jj] = ToString(fgg.EdgeObjectIndexByIndexFrom(ii,jj));
            }
        }
        Ofstream( dot_out, filename+".dot" );
        fgg.DOT_vl( dot_out, vertex_labels,
             "",
             vec< vec<String> >( ),
             vec<String>( ),
             vec<vec<String>>(),
             vec<String>( ) ) ;
        dot_out.close();
    }
#if 0
    fastg::fastg_graph fgg(graph);
//    fgg.Compact();
    fgg.Write(filename);
    vec< vec<String>> edge_labels(fgg.N());
    vec< String> vertex_labels(fgg.N());
    for(size_t ii=0;ii<edge_labels.size();++ii){
        vertex_labels[ii]=ToString(ii);
        int nn = fgg.FromSize(ii);
        edge_labels[ii].resize(nn);
        for(int jj = 0 ; jj < nn ; ++jj){
            edge_labels[ii][jj] = ToString(fgg.EdgeObjectIndexByIndexFrom(ii,jj));
        }
    }
    Ofstream( dot_out, filename+".dot" );
    /*
    fgg.DOT_vl( dot_out, vertex_labels,
         "",
         vec< vec<String> >( ),
         vec<String>( ),
         edge_labels,
         vec<String>( ) ) ;
    */
    fgg.DOT_el(dot_out,"", edge_labels);
    dot_out.close();
#endif
}
} /* namespace fastg */
#include "graph/DigraphTemplate.h"
template void
digraphVE<fastg::fastg_k, fastg::fastg_sequence>::Initialize(
        const vec<vec<int> >& from, const vec<vec<int> >& to,
        const vec<fastg::fastg_k>& verts,
        const vec<fastg::fastg_sequence>& edges, const vec<vec<int> >& to_edge_obj,
        const vec<vec<int> >& from_edge_obj);
template void
digraphVE<fastg::fastg_sequence,fastg::fastg_k>::AddVertex( const fastg::fastg_sequence& v );

