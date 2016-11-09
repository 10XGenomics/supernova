///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "10X/paths/ReadPathVecX.h"
#include <omp.h>

// Create a simple ruler of length size that looks like:
// 0.........10........20........30........40......
String MakeRuler(const size_t size) {
    String ruler(size, '.');
    String mark = "";
    size_t mark_pos = 0;
    for (size_t i = 0; i < size; ++i) {
	if (i%10 == 0) {
	    mark = ToString(i);
	    mark_pos = 0;
	}
	if (mark_pos < mark.size()) 
	    ruler[i] = mark[mark_pos++];
    }
    return ruler;
}


bool ValidateReadPath(const HyperBasevector& hbv, const vec<int>& to_left,
		      const vec<int>& to_right, const int offset,
		      const vec<int>& edge_list, String& message, 
		      const int read_length = 0) {

    int edge_count = hbv.EdgeObjectCount();
    uint K = hbv.K();

    // Check edge ids are in bounds
    for (auto edge_id : edge_list)
	if (edge_id >= edge_count) {
	    message = "ERROR - Invalid edge ID: " + ToString(edge_id);
	    return false;
	}
    
    // Validate path through graph
    for (uint i = 0 ; i < edge_list.size(); i++) {
	int edge = edge_list[i];
	int right_v = to_right[edge];
	if ((i < edge_list.size() - 1) && (to_left[edge_list[i+1]] != right_v)) {
	    message = "ERROR - no connection between edges : " + ToString(edge) 
		+ " and " + ToString(edge_list[i+1]);
	    return false;
	}
    }
    
    // Validate negative offset
    if (read_length != 0 && offset < 0 && -offset >= static_cast<int>(read_length)) {
	message = "ERROR - negative offset exceeds read length of: " 
	    + ToString(offset);
	return false;
    }
    
    // Validate positive offset
    int first_edge_length = hbv.EdgeLengthBases(edge_list[0]);
    if (offset < 0 && offset >= first_edge_length) {
	message = "ERROR - postive offset exceeds first edge length of: " 
	    + ToString(first_edge_length);
	return false;
    }

    // Check path through graph (only possible if we know the read length)
    if (read_length != 0) {
	int used_bases =  (offset < 0 ? -offset : 0);
	int edge_start_pos = (offset < 0 ? 0 : offset);
	
	for (size_t i = 0 ; i < edge_list.size(); i++) {

	    int remaining_bases = read_length - used_bases;
	    if ( remaining_bases == 0 ) {
		message = "ERROR - path extends beyond read onto edge: " 
		    + ToString(edge_list[i]);
		return false;
	    }

	    int edge_size = hbv.EdgeLengthBases(edge_list[i]);
	    int trimmed_edge_size = edge_size - ( i == edge_list.size() - 1 ? 0 : (K - 1)); 
	    int edge_end_pos = Min(edge_start_pos + remaining_bases, trimmed_edge_size); 
	    
	    if (edge_start_pos >= edge_end_pos) 
		edge_end_pos = edge_start_pos;

	    used_bases += (edge_end_pos - edge_start_pos);
	    
	    if (edge_start_pos == edge_end_pos)
		edge_start_pos = ( K - 1) - (edge_size - edge_start_pos);
	    else 
		edge_start_pos = 0;
	}
    }

    // Add test for ambigious overlap?

    message = "";
    return true;
}

bool ValidateAllReadPaths(const HyperBasevector& hbv, const HyperBasevectorX& hb,
        const ReadPathVecX& readpaths ) {
    // Compute left and right indices
    vec<int> to_left, to_right;
    hbv.ToLeft(to_left);
    hbv.ToRight(to_right);

    String message;
    bool found_invalid = false;
    const int max_prints = 10;
    int count = 0;
    #pragma omp parallel for schedule(dynamic,1000)
    for (size_t i = 0; i < (size_t) readpaths.size();  ++i) {
        if(count<=max_prints){
            ReadPath path; readpaths.unzip(path,hb,i);
            if (!path.empty() 
                && ValidateReadPath(hbv, to_left, to_right, path.getOffset(), 
                        vec<int>(path.begin(), path.end()), message ) == false) {
                found_invalid = true;
                #pragma omp critical
                { count++;
                  if(count==1) cout<<"Unordered output of contentious paths..."<<endl;
                  cout << "Path " << i << " of " << readpaths.size( ) << " = " 
                         << path.getOffset() << ":" 
                  << printSeq( path ) << "  " << message <<endl;
                }
            }
        }
    }
    if(count>max_prints)
        cout<<"..."<<endl;
    return !found_invalid;
}

bool ValidateAllReadPaths(const HyperBasevector& hbv, const ReadPathVecX& readpaths ) {
    // Compute left and right indices
    vec<int> to_left, to_right;
    hbv.ToLeft(to_left);
    hbv.ToRight(to_right);
    HyperBasevectorX hb;
    hb = HyperBasevectorX(hbv);

    String message;
    bool found_invalid = false;
    const int max_prints = 10;
    int count = 0;
    #pragma omp parallel for schedule(dynamic,1000)
    for (size_t i = 0; i < (size_t) readpaths.size();  ++i) {
        if(count<=max_prints){
            ReadPath path; readpaths.unzip(path,hb,i);
            if (!path.empty() 
                && ValidateReadPath(hbv, to_left, to_right, path.getOffset(), 
                        vec<int>(path.begin(), path.end()), message ) == false) {
                found_invalid = true;
                #pragma omp critical
                { count++;
                  if(count==1) cout<<"Unordered output of contentious paths..."<<endl;
                  cout << "Path " << i << " of " << readpaths.size( ) << " = " 
                         << path.getOffset() << ":" 
                  << printSeq( path ) << "  " << message <<endl;
                }
            }
        }
    }
    if(count>max_prints)
        cout<<"..."<<endl;
    return !found_invalid;
}

bool ValidateAllReadPaths(const HyperBasevector& hbv, const ReadPathVec& readpaths ) {
    // Compute left and right indices
    vec<int> to_left, to_right;
    hbv.ToLeft(to_left);
    hbv.ToRight(to_right);

    String message;
    bool found_invalid = false;
    const int max_prints = 10;
    int count = 0;

    #pragma omp parallel for schedule(dynamic,1000)
    for (size_t i = 0; i < (size_t) readpaths.size();  ++i) {
        if(count<=max_prints){
            const ReadPath& path = readpaths[i];
            if (!path.empty() 
                && ValidateReadPath(hbv, to_left, to_right, path.getOffset(), 
                        vec<int>(path.begin(), path.end()), message ) == false) {
                found_invalid = true;
                #pragma omp critical
                { count++;
                  if(count==1) cout<<"Unordered output of contentious paths..."<<endl;
                  cout << "Path " << i << " of " << readpaths.size( ) << " = " 
                         << path.getOffset() << ":" 
                  << printSeq( path ) << "  " << message <<endl;
                }
            }
        }
    }
    if(count>max_prints)
        cout<<"..."<<endl;
    return !found_invalid;
}

void DisplayReadPath(ostream &out, const vecbvec& graph_edges,
                     const HyperBasevector& hbv,
                     const vec<int>& to_left, const vec<int>& to_right,
		     const int offset, const vec<int>& edge_list,
		     const BaseVec& seq, const QualVec& quals,
		     bool show_edges, bool show_ruler, bool show_seq,
		     bool show_quals, bool show_mismatch, bool show_overlap) {
    
    show_seq = show_seq && seq.size() != 0;
    show_quals = show_quals && quals.size() != 0;
    show_mismatch = show_seq && show_mismatch;

    // Can we navigate the graph cheaply?
    bool edge_lookup = !(to_left.empty() || to_right.empty());
    bool have_hbv = (hbv.EdgeObjectCount() != 0);
    bool have_edges = !graph_edges.empty();
    
    uint K = hbv.K();
    uint read_length = seq.size();

    // Validate offset
    if (offset < 0 && -offset >= static_cast<int>(read_length)) {
	cout << "ERROR - Invalid path offset: " << offset << endl;
	return;
    }

    vec<uint> overlap_count(read_length, 0);
    BaseVec edge_path;
    vec<String> edges;

    uint used_bases = 0;
    uint edge_start_pos = (offset < 0 ? 0 : offset);
  
    // Add padding to left side of alignment
    uint padding = (offset < 0 ? -offset : 0);
    String paddingStr = String(padding, ' ');
    if (padding > 0) {
	used_bases = padding;
    }

    bool start_cap = (edge_start_pos == 0);
    bool end_cap = false;

    bool left_overlap = false;
    int left_overlap_size = 0;
    bool right_overlap = false;

    bool bad_path = false;

    // pos up the path, one edge at a time
    for (uint i = 0 ; i < edge_list.size(); i++) {

	uint remaining_bases = read_length - used_bases;
	if ( remaining_bases == 0 ) {
	    bad_path = true;
	    break;
	}
	const BaseVec& edge = (have_edges ? graph_edges[edge_list[i]] : hbv.EdgeObject(edge_list[i]));

	uint trimmed_edge_size = edge.size() - ( i == edge_list.size() - 1 ? 0 : (K - 1)); 
	uint full_edge_end_pos = Min(edge_start_pos + remaining_bases, edge.size()); // without K-1 trimming
	uint edge_end_pos = Min(edge_start_pos + remaining_bases, trimmed_edge_size);  // with K-1 trimming, except for last edge

	if ( edge_start_pos >= edge.size() ) {
	    bad_path = true;
	    break;
	}


	if (edge_start_pos >= edge_end_pos) {
	    left_overlap = true;
	    edge_end_pos = edge_start_pos;
	}

	if (remaining_bases < K) {
	    right_overlap = true;
	}

	//	PRINT4(i, edge.size(), remaining_bases, used_bases);
	//	PRINT4(edge_start_pos, edge_end_pos, trimmed_edge_size, full_edge_end_pos);

	// Build path basevector
	edge_path.append(edge.begin() + edge_start_pos, edge.begin() + edge_end_pos);
    
	// build vector of overlapping edges (trimmed to match the read)
	if (show_edges) {
	    BaseVec whole_edge(edge.begin() + edge_start_pos, edge.begin() + full_edge_end_pos);
	    String this_edge(used_bases, ' ');
	    this_edge.append(ToString(whole_edge));
	    edges.push_back(this_edge);
	}

	// update read overlap counter
	if (show_overlap)
	    for (uint j = 0; j < (full_edge_end_pos - edge_start_pos); j++) 
		overlap_count[used_bases + j]++;
  
	used_bases += (edge_end_pos - edge_start_pos);

	// Add end cap
	if (i == edge_list.size() - 1 && edge_end_pos == trimmed_edge_size )
	    end_cap = true;

	if (edge_start_pos == edge_end_pos)
	    edge_start_pos = ( K - 1) - (edge.size() - edge_start_pos);
	else 
	    edge_start_pos = 0;
    }
   
    // Build mismatch string
    String mismatches(read_length);
    if (show_mismatch) {
	for (int read_pos = 0; read_pos < static_cast<int>(read_length); read_pos++ ) {
	    int path_pos = read_pos - padding;
	    if (path_pos >= 0 && path_pos < static_cast<int>(edge_path.size()) 
		&& seq[read_pos] != edge_path[path_pos])
		mismatches[read_pos] = '*';
	    else 
		mismatches[read_pos] = ' ';
	}
    }

    // Build overlap string
    String overlaps(read_length);
    if (show_overlap) {
	for (uint i = 0; i < read_length; i++ )
	    if (overlap_count[i] == 0)
		overlaps[i] = ' ';
	    else if (overlap_count[i] == 1)
		overlaps[i] = '-';
	    else if (overlap_count[i] < 10)
		overlaps[i] = ToString(overlap_count[i])[0];
	    else  // >9 overlapping edges
		overlaps[i] = '#'; 
    }

    // Add edge caps
    if (start_cap && edge_lookup && have_hbv)
	if (hbv.Source(to_left[edge_list.front()]))
	    overlaps[padding] = '[';
	else
	    overlaps[padding] = '+';
    if (end_cap && edge_lookup & have_hbv)
	if (hbv.Sink(to_right[edge_list.back()]))
	    overlaps[used_bases - 1] = ']';
	else
	    overlaps[used_bases - 1] = '+';

    if (show_quals)    PrintStacked(out, quals);
    if (show_mismatch) out << mismatches << endl;
    if (show_seq)        out << seq << endl;
    out << paddingStr << edge_path << endl;
    if (show_overlap)  out << overlaps << endl;
    if (show_edges)   for (String& line : edges)  out << line << endl;
    if (show_ruler)  out << MakeRuler(read_length) << endl;
    if (edge_lookup && have_hbv) {
	if (left_overlap && (hbv.ToSize(to_left[edge_list[1]]) > 1) )
	    cout << "Warning: path extends ambgiously to the left >-" << endl;
	if (right_overlap && (hbv.FromSize(to_right[edge_list[edge_list.size() - 2]]) > 1) )
	    cout << "Warning: path extends ambgiously to the right -<" << endl;
    }
    if ( bad_path )
	cout << "Warning: a path was encountered where the read finished before the path" << endl;
}



void DisplayReadPath(ostream &out, const HyperBasevector& hbv,
		     const vec<int>& to_left, const vec<int>& to_right,
		     const int offset, const vec<int>& edge_list,
		     const BaseVec& seq, const QualVec& quals,
		     bool show_edges, bool show_ruler, bool show_seq,
		     bool show_quals, bool show_mismatch, bool show_overlap) {

    DisplayReadPath(out, vecbvec(), hbv, to_left, to_right,
		    offset,  edge_list,
		    seq,  quals,
		    show_edges, show_ruler, show_seq, show_quals,
		    show_mismatch, show_overlap);
}

void DisplayReadPath(ostream &out, const vecbvec& graph_edges, int K,
		     const int offset, const vec<int>& edge_list,
		     const BaseVec& seq, const QualVec& quals,
		     bool show_edges, bool show_ruler, bool show_seq,
		     bool show_quals, bool show_mismatch, bool show_overlap) {
  
    DisplayReadPath(out, graph_edges, HyperBasevector(K), vec<int>(), vec<int>(),
		    offset,  edge_list,
		    seq,  quals,
		    show_edges, show_ruler, show_seq, show_quals,
		    show_mismatch, show_overlap);
}
