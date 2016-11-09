///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
//

#include "Intvector.h"
#include "Qualvector.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "paths/long/ReadPathTools.h"
#include "kmers/naif_kmer/KernelPerfectAligner.h"
#include "kmers/naif_kmer/NaifKmerizer.h"
#include "kmers/naif_kmer/Kmers.h"
#include "paths/long/ExtendReadPath.h"
#include "paths/long/ShortKmerReadPather.h"

struct Alignment {
    int score;
    size_t edge;
    int offset;

    Alignment() : score(0), edge(0), offset(0) {};

    Alignment(int offset, size_t edge, int score) : score(score), edge(edge), offset(offset) {};
};


std::tuple<unsigned, unsigned, unsigned>  // score, mismatch, length
ShortKmerReadPather::ScoreEdge( bvec const& bases, qvec const& quals,
		    int offset, bvec const& edge) {
	
    auto bitr = bases.begin() + (offset > 0 ? 0 : -offset);
    auto qitr = quals.begin() + (offset > 0 ? 0 : -offset);
    auto eitr = edge.begin() + (offset > 0 ? offset : 0);
	
    unsigned qSum = 0;
    unsigned penalty = 0;
    unsigned mismatch = 0;
    unsigned length = 0;
    while ( bitr != bases.end() && qitr != quals.end() && eitr != edge.end() ) {
	// cout << (bitr - bases.begin()) << " " << (eitr - edge.begin() ) << " ";
	// cout << Base::val2Char(*bitr) << " " << Base::val2Char(*eitr) << endl;
	if ( *bitr != *eitr ) {
	    // transform Q2 -> Q20
	    auto qscore = ( *qitr == 2 ) ? 20 : *qitr;
	    penalty += qscore;
	    qSum += penalty;
	    mismatch++;
	} else if ( penalty > 0 ) {
	    penalty -= (0.2*penalty);
	}
	//          PRINT5(Base::val2Char(*bitr), Base::val2Char(*eitr), (int)*qitr, qSum, penalty);
	bitr++;
	qitr++;
	eitr++;

	length++;
    }

    // penalize left over bases on the read
    while ( bitr++ != bases.end() && qitr++ != quals.end() )
	qSum += 10;
	
    return std::make_tuple(qSum, mismatch, length);
}
    

void ShortKmerReadPather::FindPaths(const vecbasevector& bases,
                        const VecPQVec& quals,
                        const HyperBasevector& hbv, ReadPathVec& paths,
                        const int K, const int num_threads,
                        const vec<size_t>& ids, bool debug ) {

    vec<int> to_left, to_right;
    hbv.ToLeft(to_left);
    hbv.ToRight(to_right);
	
    bool have_ids = !ids.empty();
    size_t edge_count = hbv.EdgeObjectCount();
    size_t read_count = (have_ids ? ids.size() : bases.size());

    // Prepare bases for kmerization. First the reads
    vecbasevector seqs;
    if ( have_ids ) {
	seqs.reserve(ids.size() + edge_count);
	for ( auto i : ids) 
	    seqs.push_back(bases[i]);
    } else {
	seqs.reserve(bases.size() + edge_count);
	seqs.Append(bases);
    }
    // Now the graph edges
    seqs.append(hbv.Edges().begin(), hbv.Edges().end());

    if (debug) {
	cout << "HBV contains " << edge_count << " edges" << endl;
	cout << "Processing " << read_count << " of " << bases.size() << " reads" << endl;
    }

    double xclock = WallClockTime( );
    
    ForceAssertLe(K, 60);  // ensure that K is smaller than Kmer_t
    typedef Kmer60 Kmer_t;    // must be >= K
    typedef triple<int64_t, int, int64_t>  KmerAlign;

    KmerAligns p;  // Kmer aligmnets: read id, offset on edge, edge id
    // Kmerize!
    KernelAllKmerAligns<Kmer_t> kpa( seqs, read_count, K, &p,
				     true /* ignore palindromes */, true /* ignore relative RC */ );
    naif_kmerize(&kpa, num_threads, (debug ? 1 : 0));

    cout << TimeSince(xclock) << " used kmerizing" << endl;
    xclock = WallClockTime( );


    // Remove redundant information and sort by read, offset then edge
    UniqueSort(p);

    // Define less that only looks at the path id
    auto path_less   = [] (const KmerAlign&  one, const KmerAlign& two) 
	{ return (one.first < two.first); };
    // Define less that looks at the path id then offset only (ignore target id)
    auto offset_less = [] (const KmerAlign&  one, const KmerAlign& two)
	{ return (one.second < two.second); };

    // Process the kmer aligments
    qvec qv;
    size_t success_count = 0;
    auto next_pos = p.begin();
    while (next_pos != p.end()) {
	
	// Find next range with matching read ids
	decltype (next_pos)  read_start, read_end;
	tie (read_start, read_end) =  std::equal_range (next_pos, p.end(), *next_pos, path_less); 

	// Skip cases with too many placements
	if (read_end - read_start > 100 ) {
	    next_pos = read_end;
	    continue;
	}

	size_t read_id = (have_ids ? ids[next_pos->first] : next_pos->first);
	quals[read_id].unpack(&qv);
	// Examine all kmer alignments for this read
	vec<Alignment> aligns;
     //	vec<int> edge_list, offset_list, score_list;
	while (next_pos != read_end) {
	
	    // Find range with matching offsets and pick the best one
	    decltype (next_pos)  offset_start, offset_end;
	    tie (offset_start, offset_end) =  std::equal_range (next_pos, read_end, *next_pos, offset_less); 

	    auto best = p.end();
	    bool tied = false; // don't use edges with best matching scores
	    unsigned best_score = std::numeric_limits<int>::max(); // Lower is better
	    for (auto entry = offset_start; entry < offset_end; entry++) {
			
		int offset = -(entry->second);
		basevector edge_bv = hbv.EdgeObject(entry->third);
		unsigned mismatch, length, score;
		tie (score, mismatch, length) =  ScoreEdge( bases[read_id], qv, offset, edge_bv);

		if ( (static_cast<double>(mismatch) / (length - K) < 0.4) /* mismatch percentage smaller than 50% */
		     && (score < 10000) ) {
		    if ( score < best_score) { 
			best = entry;
			best_score = score;
			tied = false;
		    } else if ( score == best_score)
			tied = true;
		}
		if (debug) 
		    cout << "read=" << read_id << " offset=" << offset << " edge=" << entry->third << " size=" 
			 << edge_bv.size() << " score=" << score << " mismatch=" <<  mismatch 
			 << " length=" << length << (best == entry ? " *" : "") 
			 << (tied ? " Tied" : "") << endl;
	    }

	    // add best matching edge to the list
	    if (best != p.end() && !tied)
		aligns.emplace_back(-(next_pos->second), best->third, best_score);
	    
	    next_pos = offset_end;
	}
	    
	// Validate and extend the read path
		
	String message;
	if (!aligns.empty()) {
		
	    // Deal with unequal length bubbles, or branching in
	    vec<int> to_delete;
	    for (size_t i = 0 ; i < aligns.size() - 1; i++) {
		size_t edge = aligns[i].edge;
		size_t next_edge = aligns[i+1].edge;
		if (to_right[edge] == to_right[next_edge] || to_left[edge] == to_left[next_edge]) {
		    to_delete.push_back(aligns[i].score < aligns[i+1].score ? i+1 : i);
		    if (debug)  cout << "Deleting edge " << aligns[to_delete.back()].edge  << " from " 
				     << edge << "," << next_edge << endl;
		}
	    }
	    EraseTheseIndices( aligns, to_delete );

	    int initial_offset = aligns.front().offset;
	    vec<int> edge_list(aligns.size());
	    std::transform(aligns.begin(), aligns.end(), edge_list.begin(), [] (Alignment i) {return i.edge;});
		
	    if (debug) cout << "Path=" << initial_offset << ":" << printSeq(edge_list) << endl;

	    ReadPath read_path(initial_offset, edge_list);

	    // Validate and if the read path is invalid, pick the highest scoring edge
	    if (false == ValidateReadPath(hbv, to_left, to_right, initial_offset, edge_list, message)) {
		
		if (debug)
		    cout << "INVALID hbv_path= " << initial_offset << ":"<< printSeq( edge_list ) << " : " 
			 << message << endl;
		continue;

		// if ( edge_list.size() < 5) {
		//     size_t min_entry = std::min_element(score_list.begin(), score_list.end()) - score_list.begin();
		//     read_path = ReadPath(offset_list[min_entry], vec<int>(1,edge_list[min_entry]));
		//     if (debug) cout << "Picking edge :" << edge_list[min_entry] << " instead" << endl;
		// } else
		//     continue;  // To many edges to pick from
	    }

	    ExtendReadPath::attemptLeftRightExtension( read_path, bases[read_id], qv,
							   hbv, to_left, to_right, false);

	    if (debug) {
		vec<int> new_edge_list(read_path.begin(), read_path.end());
		int new_offset = read_path.getOffset();
		
		cout << read_id << "  hbv_path= " << new_offset << ":"<< printSeq( new_edge_list ) << endl;
		DisplayReadPath( cout, hbv, to_left, to_right,
				 new_offset, new_edge_list,
				 bases[read_id], qv);
		cout << endl << endl;
	    }
	    
	    paths[read_id].swap(read_path);
	    success_count++;
	    
	}
	next_pos = read_end;
    }

    cout << TimeSince(xclock) << " used exploring paths" << endl;
    if (debug) 
	cout << "Found " << success_count << " paths" << endl;
}

