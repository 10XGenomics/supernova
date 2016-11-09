///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
#include "paths/long/OverlapReads.h"
#include "VecUtilities.h"
#include "FeudalMimic.h"
#include <queue>

// ================================ static methods =============================

// If tail b1[len1-overlap:len1) is the same as head b2[0: overlap)
bool ReadOverlapGraph::BaseMatch( const basevector& b1,  const basevector& b2, int overlap )
{
    basevector::const_iterator it1 = b1.end() - overlap;
    basevector::const_iterator it2 = b2.begin();
    if ( it1 == b1.end() ) return true;
    if ( it2 == b2.end() ) return false;
    if ( *it1 != *it2 ) return false;
    while ( true ) {
        if ( ++it1 == b1.end() ) 
            if ( ++it2 == b2.end() ) 
                return false;
            else
                return true;
        else 
            if ( ++it2 == b2.end() ) 
                return false;
            else
                if ( *it1 != *it2 ) return false;
    }
    cout << " Something is wrong " << endl;
    cout << b1.ToString() << endl;
    cout << b2.ToString() << endl;
    CRD::exit(1);
}


unsigned int ReadOverlapGraph::Overlap( const basevector& r1, const basevector& r2, unsigned int overlap_lb, unsigned int overlap_ub )
{   
    //cout << "overlap between " << r1.ToString() << endl;
    //cout << "          and   " << r2.ToString() << endl;
    //PRINT2( overlap_lb, overlap_ub );
    unsigned int overlap = min( r1.size(), overlap_ub) -1;
    for( ; overlap >= overlap_lb; overlap-- ) 
        if ( BaseMatch( r1, r2, overlap ) ) return overlap;
    return 0;
}


// ==============================

void ReadOverlapGraph::Init( int Min_Overlap, const vec< pair<int,int> > *pt_starts ) 
{
    disable_.clear(); 
    Mimic( reads_, disable_, false );
    vec <ReadOverlap> er_overlaps;
    if ( pt_starts == 0 ) { // not range estimate provided
        for( size_t rid1 = 0; rid1 < NReads(); rid1++ ) 
            for( size_t rid2 = rid1+1; rid2 < NReads(); rid2++ ) {
                FindReadOverlaps( rid1, rid2, Min_Overlap, Len(rid1), &er_overlaps );
                FindReadOverlaps( rid2, rid1, Min_Overlap, Len(rid2), &er_overlaps );
            }
    }
    else {
        const vec< pair<int,int> >& starts = *pt_starts;
        for( size_t rid1 = 0; rid1 < NReads(); rid1++ ) {
            int start1_lb = starts[rid1].first;
            int start1_ub = starts[rid1].second;
            for( size_t rid2 = rid1+1; rid2 < NReads(); rid2++ ) {
                int start2_lb = starts[rid2].first;
                int start2_ub = starts[rid2].second;
                int overlap12_lb = start1_lb + Len(rid1) - start2_ub;
                int overlap12_ub = start1_ub + Len(rid1) - start2_lb;
                if ( overlap12_ub > Min_Overlap && overlap12_lb < (int)Len(rid1) )
                    FindReadOverlaps( rid1, rid2, max( Min_Overlap, overlap12_lb ), 
                            overlap12_ub, &er_overlaps );
                int overlap21_lb = start2_lb + Len(rid2) - start1_ub;
                int overlap21_ub = start2_ub + Len(rid2) - start1_lb;
                if ( overlap21_ub > Min_Overlap && overlap21_lb < (int)Len(rid2) )
                    FindReadOverlaps( rid2, rid1, max( Min_Overlap, overlap21_lb ), 
                            overlap21_ub, &er_overlaps );
            }
        }
    }
    froms_.clear(); 
    Mimic( reads_, froms_ );
    tos_.clear();
    Mimic( reads_, tos_ );
    for ( size_t i = 0; i < er_overlaps.size(); ++i ) {
        //cout << "ovelap " << i << ": ";
        //cout << er_overlaps[i].rid1 << " " << er_overlaps[i].rid2 << " " << er_overlaps[i].alt1
        //    << " " << er_overlaps[i].alt2 << " " << er_overlaps[i].overlap << endl;
        if ( disable_[ er_overlaps[i].rid1 ][ er_overlaps[i].alt1 ]
           ||disable_[ er_overlaps[i].rid2 ][ er_overlaps[i].alt2 ] ) 
            continue;
        froms_[ er_overlaps[i].rid1 ][ er_overlaps[i].alt1 ].push_back( er_overlaps[i] );
        tos_[ er_overlaps[i].rid2 ][ er_overlaps[i].alt2 ].push_back( er_overlaps[i] );
    }
    RemoveTransitive();
    //CombineAmbReads();
}


// If there are overlaps A->B, B->C, and A->C, remove the redundent A->C edge
void ReadOverlapGraph::RemoveTransitive( )
{
    int count = 0;
    vec< vec<bool> > visited;
    Mimic( reads_, visited, false );
    for ( size_t i = 0; i < froms_.size(); ++i ) {
        for ( size_t j = 0; j < froms_[i].size(); ++j ) {
            if (  ! tos_[i][j].empty() ) continue; // only source
            count += RemoveTransitiveFrom(i,j, visited);
        }
    }
    //cout << "Removed " << count << " edges." << endl;
    // update the tos
    tos_.clear();
    Mimic( reads_, tos_ );
    for ( size_t i = 0; i < froms_.size(); ++i ) {
        for ( size_t j = 0; j < froms_[i].size(); ++j ) {
            for ( size_t k = 0; k < froms_[i][j].size(); ++k ) {
                tos_[ froms_[i][j][k].rid2 ][ froms_[i][j][k].alt2 ].push_back
                    ( froms_[i][j][k] );
            }
        }
    }
}

// Remove the transit edge incident from node ( rid, alt ). 
// Mark the visited nodes, return the total number of edges removed.
int ReadOverlapGraph::RemoveTransitiveFrom(int rid, int alt, vec< vec<bool> >& visited )
{
    int count = 0;
    queue< pair<int,int> > que; 
    que.push( make_pair(rid,alt) );
    while ( ! que.empty() ) {
        pair<int,int> node = que.front();
        que.pop();
        int i = node.first, j = node.second;
        if ( visited[i][j] ) continue;
        visited[i][j] = true;
        // delete the transit edge out from this node 
        vec<Bool> to_del(  froms_[i][j].size(), False );
        for ( size_t k1 = 0; k1 < froms_[i][j].size(); ++k1 ) {
            const ReadOverlap& ro1 = froms_[i][j][k1];
            for ( size_t k2 = k1+1; k2 < froms_[i][j].size(); ++k2 ) {
                const ReadOverlap& ro2 = froms_[i][j][k2];
                if ( IsOverlaped(ro1.rid2, ro1.alt2, ro2.rid2, ro2.alt2) )
                    to_del[k2] = True;
                else if ( IsOverlaped(ro2.rid2, ro2.alt2, ro1.rid2, ro1.alt2) )
                    to_del[k1] = True;
            }
        }
        count += to_del.CountValue(True);
        EraseIf( froms_[i][j], to_del );
        for ( size_t k = 0; k < froms_[i][j].size(); ++k )
            que.push( make_pair( froms_[i][j][k].rid2, froms_[i][j][k].alt2 ) );
    }
    return count;
}


// Combine expanded reads if ambiguity sites do not overlap with other reads
void ReadOverlapGraph::CombineAmbReads( const vec< vec<Ambiguity> >& creads_ambs ) 
{
    ForceAssertEq( creads_ambs.size(), reads_.size() );
    combined_ambs_.assign( creads_ambs.size(), vec<Ambiguity>()  ); // empty by default
    // quick check
    bool has_amb = false;
    for ( size_t i = 0; i < reads_.size(); ++i ) 
        if ( reads_[i].size() > 1 ) { has_amb = true; break; }
    if ( !has_amb) return;
    for ( size_t i = 0; i < reads_.size(); ++i ) {
        bool mismatch = false;
        if ( reads_[i].size() < 2 ) continue;
        //cout << "check read " << i << " alts= " << reads_[i].size() << endl;
        for ( size_t j1 = 0; j1 < reads_[i].size(); ++j1 ) {
            for ( size_t j2 = j1+1; j2 < reads_[i].size(); ++j2 ) {
                if ( froms_[i][j1].size() != froms_[i][j2].size()
                   ||tos_[i][j1].size() != tos_[i][j2].size()  ) {
                    mismatch=true;
                    //cout << "number mismatch from " << j1 << " to " << j2 << endl;
                    break;
                }
                { // check the froms
                    vec< triple<int,int,int> > edges1( froms_[i][j1].size() );
                    for ( size_t k = 0; k < froms_[i][j1].size(); ++k ) 
                        edges1.push( froms_[i][j1][k].rid2, froms_[i][j1][k].alt2, froms_[i][j1][k].overlap );
                    Sort(edges1);
                    vec< triple<int,int,int> > edges2( froms_[i][j1].size() );
                    for ( size_t k = 0; k < froms_[i][j2].size(); ++k ) 
                        edges2.push( froms_[i][j2][k].rid2, froms_[i][j2][k].alt2, froms_[i][j2][k].overlap );
                    Sort(edges2);
                    if ( edges1 != edges2 ) {
                        //cout << "froms mismatch from " << j1 << " to " << j2 << endl;
                        mismatch = true;
                        break;
                    }
                }
                { // check the tos
                    vec< triple<int,int,int> > edges1( tos_[i][j1].size() );
                    for ( size_t k = 0; k < tos_[i][j1].size(); ++k ) 
                        edges1.push( tos_[i][j1][k].rid1, tos_[i][j1][k].alt1, tos_[i][j1][k].overlap );
                    Sort(edges1);
                    vec< triple<int,int,int> > edges2( tos_[i][j1].size() );
                    for ( size_t k = 0; k < tos_[i][j2].size(); ++k ) 
                        edges2.push( tos_[i][j2][k].rid1, tos_[i][j2][k].alt1, tos_[i][j2][k].overlap );
                    Sort(edges2);
                    if ( edges1 != edges2 ) {
                        //cout << "tos  mismatch from " << j1 << " to " << j2 << endl;
                        mismatch = true;
                        break;
                    }
                }
            }
            if ( mismatch ) break;
        }
        if ( ! mismatch ) {
            combined_ambs_[i] = creads_ambs[i];
            //reads_[i].resize(1);
            tos_[i].resize(1);
            froms_[i].resize(1);
            for ( size_t j = 0; j < froms_[i][0].size(); ++j ) {
                int rid2 = froms_[i][0][j].rid2;
                int alt2 = froms_[i][0][j].alt2;
                vec<Bool> to_del( tos_[rid2][alt2].size(), False );
                for ( size_t k = 0; k < tos_[rid2][alt2].size(); ++k ) 
                    if ( tos_[rid2][alt2][k].rid1 == (int)i && tos_[rid2][alt2][k].alt1 != 0 )
                        to_del[k] = True;
                EraseIf( tos_[rid2][alt2], to_del );
            }
            for ( size_t j = 0; j < tos_[i][0].size(); ++j ) {
                int rid1 = tos_[i][0][j].rid1;
                int alt1 = tos_[i][0][j].alt1;
                vec<Bool> to_del( froms_[rid1][alt1].size(), False );
                for ( size_t k = 0; k < froms_[rid1][alt1].size(); ++k ) 
                    if ( froms_[rid1][alt1][k].rid2 == (int)i && froms_[rid1][alt1][k].alt2 != 0 )
                        to_del[k] = True;
                EraseIf( froms_[rid1][alt1], to_del );
            }
        }
    }
    //cout << "Combined " << combined_.CountValue(true) << " amb reads " << endl;
}

// Stitch the efastas together given the estimated overlaps
void ReadOverlapGraph::GetAssembly(vec<basevector>& assemblies, 
        vec< vec<ReadOverlap> >& paths, int max_path ) const
{
    // find all paths
    assemblies.clear();
    paths.clear();
    vec<int> lens;
    for ( size_t i = 0; i < froms_.size(); ++i ) {
        for ( size_t j = 0; j < froms_[i].size(); ++j ) {
            if ( ! tos_[i][j].empty() ) continue; // only source nodes
            if ( froms_[i][j].empty() ) continue; // no paths
            AllPathsFrom( i, j, &paths, &lens, max_path );
        }
    }
    for ( size_t i = 0; i < paths.size(); ++i ) {
        const vec<ReadOverlap>& path = paths[i];
        const basevector& read0 = reads_[ path[0].rid1 ][ path[0].alt1 ];
        //// from end to end
        //basevector assembly( read0 );
        //for ( size_t i = 0; i < path.size(); ++i ) {
        //    basevector::const_iterator it = reads_[ path[i].rid2 ][ path[i].alt2 ].begin() 
        //                                  + path[i].overlap;
        //    if ( it <  reads_[ path[i].rid2 ][ path[i].alt2 ].end() )
        //        assembly.append( it, reads_[ path[i].rid2 ][ path[i].alt2 ].end() );
        //}
        // Exclude the non-overlapping region of the end reads
        basevector assembly( read0.begin() + read0.size() - path[0].overlap,
                read0.end() );
        for ( size_t i = 0; i < path.size()-1; ++i ) {
            basevector::const_iterator it = reads_[ path[i].rid2 ][ path[i].alt2 ].begin() 
                                          + path[i].overlap;
            if ( it <  reads_[ path[i].rid2 ][ path[i].alt2 ].end() )
                assembly.append( it, reads_[ path[i].rid2 ][ path[i].alt2 ].end() );
        }
        assemblies.push_back( assembly );
    }
}

void ReadOverlapGraph::GetAmbsOnPath( const vec<ReadOverlap>&  path, vec<Ambiguity>& new_ambs ) const
{
    new_ambs.clear();
    int shift = 0;
    int rid1 = path[0].rid1;
    // Exclude the beginning and the end of end reads
    const basevector& read0 = reads_[ path[0].rid1 ][ path[0].alt1 ];
    shift -= (int)read0.size() - path[0].overlap;
    for ( size_t j = 0; j < combined_ambs_[rid1].size(); ++j ) {
        if ( (int)combined_ambs_[rid1][j].start + shift < 0 ) continue;
        new_ambs.push_back( combined_ambs_[rid1][j] );
        new_ambs.back().start += shift;
    }
    for ( size_t i = 0; i < path.size(); ++i ) {
        int rid1 = path[i].rid1;
        int alt1 = path[i].alt1;
        int rid2 = path[i].rid2;
        shift += (int)reads_[rid1][alt1].size() - path[i].overlap;
        for ( size_t j = 0; j < combined_ambs_[rid2].size(); ++j ) {
            new_ambs.push_back( combined_ambs_[rid2][j] );
            new_ambs.back().start += shift;
        }
    }
}

// check if all the reads can perfectly form a linear assembly and record the starting
// position of each read on the final assembly
void ReadOverlapGraph::PerfectAssembly( vec<basevector>& assemblies, vec< vec<ReadOverlap> >& paths ) const
{
    const unsigned int Max_Assembly = 10;
    assemblies.clear();
    paths.clear();
    vec< vec<int> > read_locs;
    Mimic( reads_, read_locs, -1 );
    while (true) {
        if ( assemblies.size() >= Max_Assembly ) break;
        int rid = -1, alt = -1; // start from which read?
        for ( size_t i = 0; i < read_locs.size(); ++i ) 
            for ( size_t j = 0; j < read_locs[i].size(); ++j ) 
                if ( read_locs[i][j] == -1 ) {
                    rid = i;
                    alt = j;
                }
        if ( rid == -1 || alt == -1 ) break;
        read_locs[rid][alt] = 0;
        basevector assembly = reads_[rid][alt];
        vec<ReadOverlap> path;
        if ( PerfectAssemblyForward(rid,alt,assembly,read_locs,path) && 
                PerfectAssemblyBackward(rid,alt,assembly,read_locs,path) ) {
            if ( ! path.empty() ) {
                assemblies.push_back( assembly );
                paths.push_back( path );
            }
        }
    }
    // sort
    if ( assemblies.empty() ) return;
    vec<int> lens( assemblies.size() );
    for ( size_t i = 0; i < lens.size(); ++i ) lens[i] = assemblies[i].size();
    ReverseSortSync( lens, assemblies, paths );
}

bool ReadOverlapGraph::PerfectAssemblyForward( int i, int j, 
        basevector& assembly, vec< vec<int> >& read_locs, vec<ReadOverlap>& path ) const 
{
    if ( froms_[i][j].empty() ) return true;
    int loc0 = read_locs[i][j];
    //cout << "node " << i << ", " << j << " forward connects to: ";
    //for ( size_t k = 0; k < froms_[i][j].size(); ++k ) {
    //    const ReadOverlap& ro = froms_[i][j][k];
    //    cout << ro.rid2 << " ";
    //}
    //cout << endl;
    for ( size_t k = 0; k < froms_[i][j].size(); ++k ) {
        const ReadOverlap& ro = froms_[i][j][k];
        if ( read_locs[ro.rid2][ro.alt2] != -1 ) continue; // visited
        read_locs[ro.rid2][ro.alt2] = loc0 + reads_[i][j].size() - ro.overlap;
        // check consistency
        basevector::const_iterator it = reads_[ro.rid2][ro.alt2].begin() + ro.overlap;
        basevector::const_iterator end = reads_[ro.rid2][ro.alt2].end();
        basevector::iterator it0 = assembly.begin() + read_locs[ro.rid2][ro.alt2] + ro.overlap;
        basevector::iterator end0 = assembly.end();
        while ( true ) {
            if ( it == end ) { // read end
                break;
            }
            if ( it0 == end0 ) { // assembly end
                assembly.append( it, end );
                path.push_back( ro );
                break;
            }
            if ( *it != *it0 ) { // mismatch
                return false;
            }
            ++it; ++it0;
        }
        if ( ! PerfectAssemblyForward( ro.rid2, ro.alt2, assembly, read_locs, path ) )
            return false;
        if ( ! PerfectAssemblyBackward( ro.rid2, ro.alt2, assembly, read_locs, path ) )
            return false;
    }
    return true;
}

bool ReadOverlapGraph::PerfectAssemblyBackward( int i, int j, 
        basevector& assembly, vec< vec<int> >& read_locs, vec<ReadOverlap>& path ) const 
{
    if ( tos_[i][j].empty() ) return true;
    int loc0 = read_locs[i][j];
    //cout << "node " << i << ", " << j << " back connects to: ";
    //for ( size_t k = 0; k < tos_[i][j].size(); ++k ) {
    //    const ReadOverlap& ro = tos_[i][j][k];
    //    cout << ro.rid1 << " ";
    //}
    //cout << endl;
    for ( size_t k = 0; k < tos_[i][j].size(); ++k ) {
        const ReadOverlap& ro = tos_[i][j][k];
        if ( read_locs[ro.rid1][ro.alt1] != -1 ) continue; // visited
        read_locs[ro.rid1][ro.alt1] = loc0 - reads_[ro.rid1][ro.alt1].size() + ro.overlap;
        basevector::const_iterator it = reads_[ro.rid1][ro.alt1].end() - ro.overlap;
        basevector::const_iterator begin = reads_[ro.rid1][ro.alt1].begin();
        basevector::iterator it0 = assembly.begin() + loc0;
        basevector::iterator begin0 = assembly.begin();
        while ( true ) {
            if ( it == begin ) { // read end
                break;
            }
            if ( it0 == begin0 ) { // assembly end
                path.insert( path.begin(), ro );
                basevector ext( begin, it);
                ext.append( assembly );
                swap( ext, assembly );
                int shift = distance( begin, it );
                for ( size_t i = 0; i < read_locs.size(); ++i ) 
                    for ( size_t j = 0; j < read_locs[i].size(); ++j ) 
                        read_locs[i][j] += shift;
                break;
            }
            --it; --it0;
            if ( *it != *it0 ) { // mismatch
                return false;
            }
        }
        if ( ! PerfectAssemblyBackward( ro.rid1, ro.alt1, assembly, read_locs, path ) )
            return false;
        if ( ! PerfectAssemblyForward( ro.rid1, ro.alt1, assembly, read_locs, path ) )
            return false;
    }
    return true;
}


void ReadOverlapGraph::AllPathsFrom( int rid, int alt, vec< vec<ReadOverlap> >* p_paths, 
                                     vec<int> *p_lens, int max_path ) const
{
    unsigned int Max_Que_Size = 1000000;  // prevent memory from exploding
    if ( froms_[rid][alt].empty() ) return;
    //cout << "allpaths from " << rid << " " << alt << endl;
    queue< vec<ReadOverlap> > que;
    for ( size_t i = 0; i < froms_[rid][alt].size(); ++i ) 
        que.push( MkVec(froms_[rid][alt][i]) );
    //cout << que.size() << endl;
    while ( ! que.empty() ) {
        if ( que.size() > Max_Que_Size ) break;
        vec<ReadOverlap> q = que.front();
        //cout << "get  " << q << endl;
        que.pop();
        ReadOverlap ov = q.back();
        if ( froms_[ ov.rid2 ][ ov.alt2 ].empty() ) {
            //cout << "find  " << q << endl;
            if ( p_paths->isize() < max_path ) {
                p_paths->push_back( q );
                p_lens->push_back( PathToLen(q) );
                ReverseSortSync( *p_lens, *p_paths );
            } 
            else {
                int len = PathToLen( q );
                if ( len > p_lens->back() ) {
                    p_lens->back() = len;
                    p_paths->back() = q;
                    ReverseSortSync( *p_lens, *p_paths );
                }
            }
        }
        for ( size_t i = 0; i < froms_[ov.rid2][ov.alt2].size(); ++i ) {
            const ReadOverlap& next = froms_[ov.rid2][ov.alt2][i];
            vec<ReadOverlap> qi = q;
            if ( find( qi.begin(), qi.end(), next ) != qi.end() ) { //loop
                //cout << "find  " << q << endl;
                if ( p_paths->isize() < max_path ) {
                    p_paths->push_back( q );
                    p_lens->push_back( PathToLen(q) );
                    ReverseSortSync( *p_lens, *p_paths );
                } 
                else {
                    int len = PathToLen( q );
                    if ( len > p_lens->back() ) {
                        p_lens->back() = len;
                        p_paths->back() = q;
                        ReverseSortSync( *p_lens, *p_paths );
                    }
                }
                continue;
            }
            qi.push_back( next );
            que.push(qi);
            //cout << "add  " << qi << endl;
        }
    }
}

void ReadOverlapGraph::FindReadOverlaps( int rid1, int rid2, int lb, int ub, vec<ReadOverlap> *er_overlaps ) 
{
    for( size_t i = 0; i < reads_[rid1].size(); i++ ) 
        for( size_t j = 0; j < reads_[rid2].size(); j++ ) {
            if ( disable_[rid1][i] ) continue;
            if ( disable_[rid2][j] ) continue;
            // disable duplicate reads
            if ( reads_[rid2][j].size() <= reads_[rid1][i].size() 
                    && search( reads_[rid1][i].begin(), reads_[rid1][i].end(),
                        reads_[rid2][j].begin(), reads_[rid2][j].end() )
                    != reads_[rid1][i].end() ) {
                disable_[rid2][j] = true;
                continue;
            }
            int overlap = Overlap( reads_[rid1][i], reads_[rid2][j], lb, ub );
            if ( overlap > 0 )
                er_overlaps->push_back( ReadOverlap(rid1,rid2,i,j,overlap) );
        }
}


bool ReadOverlapGraph::IsOverlaped( int rid1, int alt1, int rid2, int alt2 ) {
    for ( size_t k1 = 0; k1 < froms_[rid1][alt1].size(); ++k1 ) 
        if ( froms_[rid1][alt1][k1].rid2 == rid2 &&  froms_[rid1][alt1][k1].alt2 == alt2 )
            return true;
    return false;   
}

int ReadOverlapGraph::PathToLen( const vec<ReadOverlap>& path ) const {
    int len = Len(path[0].rid1);
    for ( size_t j = 0; j < path.size(); ++j ) 
        len += path[j].overlap;
    return len;
}

void ReadOverlapGraph::PrintGraph( ostream& out ) 
{
    for ( size_t i = 0; i < tos_.size(); ++i ) {
        out << "Read " << i << endl;
        for ( size_t j = 0; j < reads_[i].size(); ++j ) 
            out << "    " << reads_[i][j].ToString() << endl;
        out << "  Froms " << endl;
        for ( size_t j = 0; j < froms_[i].size(); ++j ) 
            for ( size_t k = 0; k < froms_[i][j].size(); ++k ) 
                out << "    " << froms_[i][j][k] << endl;
        out << "  Tos " << endl;
        for ( size_t j = 0; j < tos_[i].size(); ++j ) 
            for ( size_t k = 0; k < tos_[i][j].size(); ++k ) 
                out << "    " << tos_[i][j][k] << endl;
    }
}

// =======================================================
//  SomeReadsAssembler 
// =======================================================

SomeReadsAssembler::SomeReadsAssembler( const vec< vec<basevector> >& reads ) :graph_(reads) 
{ 
    Min_Overlap_ = 100 ; Min_Assembly_Len_ = 500; 
    Max_Assembly_Num_ = 10; 
}


vec<ReadOverlap> SomeReadsAssembler::GetAssemblyPathRegion( int i, int start, int stop ) const{
    vec<ReadOverlap> path = GetAssemblyPath(i);
    vec<Bool> to_del( path.size(), False );
    int len = graph_.Len( path.front().rid1, path.front().alt1 );
    for ( size_t j = 0; j < path.size(); ++j ) {
        len += graph_.Len( path[j].rid2, path[j].alt2 ) - path[j].overlap;
        if ( len < start ) to_del[j] = true;
        if ( len > stop ) { 
            for ( size_t k = j+1; k < path.size(); k++ )
                to_del[k] = true; 
            break;
        }
    }
    EraseIf( path, to_del);
    return path;
}

bool operator< ( const Ambiguity& lhs, const Ambiguity& rhs) { return lhs.start < rhs.start; }

efasta SomeReadsAssembler::GetEfastaAssembly( int i, size_t start, size_t stop ) const
{
    const basevector& assembly = GetAssemblies()[i];
    const vec<ReadOverlap>& path = GetAssemblyPath(i);
    const String assembly_str = assembly.ToString();
    vec<Ambiguity> new_ambs;
    graph_.GetAmbsOnPath( path, new_ambs );
    stop = Min( stop, (size_t) assembly.size() );
    if ( new_ambs.empty() ) 
        return efasta( assembly_str.substr(start, stop - start) ) ;
    // add the ambiguity sites back
    Sort(new_ambs);

    //jcout << assembly.ToString() << endl;
    ///cout << path << endl;
    //for ( size_t i = 0; i < new_ambs.size(); ++i ) 
    //    cout << new_ambs[i].to_annotation() << endl;

    vec<Bool> todel( new_ambs.size(), False );
    for ( size_t i = 0; i < new_ambs.size(); ++i )
        if ( new_ambs[i].start < start 
           ||new_ambs[i].start >= stop ) 
            todel[i] = True;
    EraseIf( new_ambs, todel );

    if ( new_ambs.empty() ) 
        return efasta( assembly_str.substr(start, stop - start) ) ;

    efasta seq( assembly_str.substr( start, new_ambs[0].start - start ) );
    size_t len = new_ambs[0].start - start;
    for ( size_t i = 0; i < new_ambs.size(); ++i ) {
        size_t j = i+1;
        while( new_ambs[j].start == new_ambs[i].start ) 
        { ForceAssertEq( new_ambs[j].size, new_ambs[i].size ); j++; }
        seq.push_back( '{' );
        seq.append( assembly_str,  new_ambs[i].start, new_ambs[i].size );
        for ( size_t k = i; k < j; ++k ) {
            seq.push_back(',');
            seq.append( new_ambs[k].replace );
        }
        seq.push_back( '}' );
        size_t amb_stop = new_ambs[i].start + new_ambs[i].size;
        if ( j != new_ambs.size() )
            seq.append( assembly_str.substr( 
                  amb_stop, new_ambs[j].start - amb_stop ) );
        else
            seq.append( assembly_str.substr( 
                  amb_stop, stop ) );
        i = j - 1;
    }
    return seq;
}

void SomeReadsAssembler::RemoveShortAssembly() {
    vec<Bool> todel( assemblies_.size(), False );
    for ( size_t i = 0; i < assemblies_.size(); ++i ) 
        if ( (int) assemblies_[i].size() < Min_Assembly_Len_ ) todel[i] = True;
    EraseIf( assemblies_, todel );
    EraseIf( paths_, todel );
}

