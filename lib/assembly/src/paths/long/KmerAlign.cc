///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "CoreTools.h"
#include "ParallelVecUtilities.h"
#include "Equiv.h"
#include "VecUtilities.h"

void KmerAlign( 
        const vec< pair<int,int> > & offset, 
        vec< pair<int,int> > & aligns,
        int VERBOSITY )
{
    aligns.clear();
    // Heuristics
    const int Max_Pos_Diff = 200; // max position difference within the group
    const double  Max_Offset_Diff = 0.2; // allow 10% offset different
    const int Min_Group_Size = 10; // offset groups less than Min_Group_Size will not be considered
    // diagnostic output
    if ( VERBOSITY >= 1 )
        cout << "    offset.size()= " << offset.size() << endl;
    if ( VERBOSITY >= 2 ) {
        cout << "Offset= ";
        for ( size_t i = 0; i < offset.size(); ++i ) 
            cout << offset[i].first << "," << offset[i].second << " ";
        cout << endl;
    }

    // Group the kmer matching, so that within the group the offsets are less
    // than 20% and pos1 differences are less than 200 bases
    equiv_rel_template <int> equiv( offset.size() );
    for ( size_t i = 0; i < offset.size(); ++i ) 
        for ( size_t j = i+1; j < offset.size(); ++j ) {
            if ( equiv.Equiv(i,j) ) continue; // already in the same group
            int pos_diff =  abs( offset[i].second - offset[j].second );
            int offset_diff = abs( offset[i].first - offset[j].first );
            if ( pos_diff < Max_Pos_Diff
                    &&  offset_diff < int( Max_Offset_Diff * pos_diff ) + 1 )
                equiv.Join( i, j );
        }
    // Sort the groups, largest first
    vec <int> reps;
    equiv.OrbitRepsAlt( reps );
    vec <int> lens;
    for ( size_t i = 0; i < reps.size(); ++i ) 
        lens.push_back(  equiv.OrbitSize( reps[i] ) );
    ReverseSortSync( lens, reps );
    if ( VERBOSITY >= 1 ) {
        String status = ( reps.solo() || equiv.OrbitSize( reps[0] ) > equiv.OrbitSize( reps[1] ) * 2 ) ? "Good":"Bad"; 
        cout << "    " << equiv.OrbitSize( reps[0] ) << "(" << equiv.OrbitSize( reps[0] ) * 100 / offset.size() << "%)" 
            << " in first group, total number of groups= " << equiv.OrbitCount() << ", status "
            << status << endl;
        for ( size_t i = 0; i < reps.size() && equiv.OrbitSize( reps[i] ) > 10  ; ++i ) 
            cout << "        Group " << i << ": " << reps[i] << ", size " << equiv.OrbitSize( reps[i] )
                << " shift~ " << offset[ reps[i] ].first << " pos1~ " << offset[ reps[i] ].second << endl;
        cout << endl;
    }

    // Starting from the largest group, deletete inconsitent groups. 
    vec< vec< pair<int,int> > > good_groups; // good roups of alignments
    for ( size_t i = 0; i < reps.size() && equiv.OrbitSize( reps[i] ) >= Min_Group_Size; ++i ) {
        vec< pair<int,int> > the_group;
        vec<int> orbit;
        equiv.Orbit( reps[i], orbit );
        for ( size_t j = 0; j < orbit.size(); ++j ) {
            int id = orbit[j];
            the_group.push( offset[id].second - offset[id].first, offset[id].second ); 
        }
        Sort( the_group );
        int pos2_left =  the_group.front().first;
        int pos2_right = the_group.back().first;
        int pos1_left =  the_group.front().second;
        int pos1_right = the_group.back().second;
        // Decide if this group of alignments overlap with previously selected ones
        bool overlap = false;
        for ( size_t j = 0; j < good_groups.size(); ++j ) 
            if ( ! ( pos2_right < good_groups[j].front().first && pos1_right < good_groups[j].front().second
                  || pos2_left  > good_groups[j].back().first  && pos1_left  > good_groups[j].back().second  ) ){
                overlap = true;
                break;
            }    
        if ( ! overlap ) good_groups.push_back( the_group ); 
    }
    vec< pair<int,int> > good_aligns;
    for ( size_t i = 0; i < good_groups.size(); ++i ) 
        good_aligns.append( good_groups[i] );
    Sort(good_aligns);
    if ( VERBOSITY >= 1 ) cout << "good aligns= " << good_aligns.size() << endl;

    // remove nonunique kmers 
    map<int,int> kmer1, kmer2;
    for ( size_t i = 0; i < good_aligns.size(); ++i ) {
        kmer1[ good_aligns[i].first ]++;
        kmer2[ good_aligns[i].second ]++;
    }
    vec<Bool> todel(good_aligns.size(),False);
    for ( size_t i = 0; i < good_aligns.size(); ++i ) 
        if ( kmer1[ good_aligns[i].first ] > 1 ||
             kmer2[ good_aligns[i].second] > 1 )
            todel[i] = True;
    EraseIf( good_aligns, todel );
    if ( VERBOSITY >= 1 ) cout << "uniq aligns= " << good_aligns.size() << endl;
    
    // Re-order and select continueously
    vec <int> orders;
    for ( size_t i = 0; i < good_aligns.size(); ++i ) 
        orders.push_back( good_aligns[i].first + good_aligns[i].second );
    SortSync( orders, good_aligns );
    aligns.clear();
    for ( size_t i = 0; i < good_aligns.size(); ++i ) 
        if ( aligns.empty() ||
             (   good_aligns[i].first > aligns.back().first &&
                 good_aligns[i].second > aligns.back().second ) )
             aligns.push_back( good_aligns[i] );

    if ( VERBOSITY >= 1 ) {
        cout << "aligns.size()= " << aligns.size() << endl;
        cout << endl;
    }
}



// Find the largest cluster from the input kmer offsets. An offset is a pair of (pos1-pos2, pos1)
// where pos1 and pos2 are the positions of the matching kmers on the two reads.

void LargestOffsetCluster( 
        const vec< pair<int,int> > & offset_in, 
        vec< pair<int,int> > & offset_out,
        int Max_Pos_Diff, double Max_Offset_Diff)
{
    offset_out.clear();
    if ( offset_in.empty() ) return;
    // Group the kmer matching, so that within the group the offsets are less
    // than 20% and pos1 differences are less than 200 bases
    equiv_rel_template <int> equiv( offset_in.size() );
    for ( size_t i = 0; i < offset_in.size(); ++i ) 
        for ( size_t j = i+1; j < offset_in.size(); ++j ) {
            if ( equiv.Equiv(i,j) ) continue; // already in the same group
            int pos_diff =  abs( offset_in[i].second - offset_in[j].second );
            int offset_in_diff = abs( offset_in[i].first - offset_in[j].first );
            if ( pos_diff < Max_Pos_Diff
                    &&  offset_in_diff < int( Max_Offset_Diff * pos_diff ) + 1 )
                equiv.Join( i, j );
        }
    // Sort the groups, largest first
    vec <int> reps;
    equiv.OrbitRepsAlt( reps );
    vec <int> lens;
    for ( size_t i = 0; i < reps.size(); ++i ) 
        lens.push_back(  equiv.OrbitSize( reps[i] ) );
    ReverseSortSync( lens, reps );
    // output the largest group
    vec<int> orbit;
    equiv.Orbit( reps[0], orbit );
    for ( size_t j = 0; j < orbit.size(); ++j ) 
        offset_out.push_back( offset_in[ orbit[j] ] );
}

// Reusable model to remove offset outliers. Used in KmerAlign3 method
void GroupAndSelect( 
        const vec< pair<int,int> > & offset_in, 
        vec< pair<int,int> > & offset_out,
        int Max_Pos_Diff, double Max_Offset_Diff ,
        int VERBOSITY)
{
    offset_out.clear();
    // Heuristics
    if ( VERBOSITY >= 1 )
        cout << "offset_in.size()= " << offset_in.size() << endl;
    if ( VERBOSITY >= 2 ) {
        cout << "offset_in= ";
        for ( size_t i = 0; i < offset_in.size(); ++i ) 
            cout << offset_in[i].first << "," << offset_in[i].second << " ";
        cout << endl;
    }
    // Group the kmer matching, so that within the group the offsets are less
    // than 20% and pos1 differences are less than 200 bases
    equiv_rel_template <int> equiv( offset_in.size() );
    for ( size_t i = 0; i < offset_in.size(); ++i ) 
        for ( size_t j = i+1; j < offset_in.size(); ++j ) {
            if ( equiv.Equiv(i,j) ) continue; // already in the same group
            int pos_diff =  abs( offset_in[i].second - offset_in[j].second );
            int offset_in_diff = abs( offset_in[i].first - offset_in[j].first );
            if ( pos_diff < Max_Pos_Diff
                    &&  offset_in_diff < int( Max_Offset_Diff * pos_diff ) + 1 )
                equiv.Join( i, j );
        }
    // Sort the groups, largest first
    vec <int> reps;
    equiv.OrbitRepsAlt( reps );
    vec <int> lens;
    for ( size_t i = 0; i < reps.size(); ++i ) 
        lens.push_back(  equiv.OrbitSize( reps[i] ) );
    ReverseSortSync( lens, reps );
    if ( VERBOSITY >= 1 ) {
        String status = ( reps.solo() || equiv.OrbitSize( reps[0] ) > equiv.OrbitSize( reps[1] ) * 2 ) ? "Good":"Bad"; 
        cout << "  " << equiv.OrbitSize( reps[0] ) << "(" << equiv.OrbitSize( reps[0] ) * 100 / offset_in.size() << "%)" 
            << " in first group, total number of groups= " << equiv.OrbitCount() << ", status "
            << status << endl;
        for ( size_t i = 0; i < reps.size() && equiv.OrbitSize( reps[i] ) > 10  ; ++i ) 
            cout << "        Group " << i << ": " << reps[i] << ", size " << equiv.OrbitSize( reps[i] )
                << " shift~ " << offset_in[ reps[i] ].first << " pos1~ " << offset_in[ reps[i] ].second << endl;
        cout << endl;
    }
    // Starting from the largest group, deletete inconsitent groups. ( Find all the good groups )
    offset_out.clear();
    vec< int > mean_offsets;
    for ( size_t i = 0; i < reps.size(); ++i ) {
        vec< pair<int,int> > the_group;
        int mean_offset = 0;
        vec<int> orbit;
        equiv.Orbit( reps[i], orbit );
        for ( size_t j = 0; j < orbit.size(); ++j ) {
            int id = orbit[j];
            the_group.push( offset_in[id].first, offset_in[id].second ); 
            mean_offset +=  offset_in[id].first;
        }
        mean_offset /= (int) orbit.size();
        bool valid = mean_offsets.empty() ? true : false;
        for ( size_t j = 0; j < mean_offsets.size(); ++j ) {
            if ( abs( mean_offset - mean_offsets[j] )
                 < abs( mean_offsets[j] ) * Max_Offset_Diff + 1 ) {
                valid = true;
                break;
            }    
        }
        if ( valid ) {
            offset_out.append( the_group ); 
            mean_offsets.push_back( mean_offset );
        }
        if ( VERBOSITY >=2 ) {
            if ( valid ) 
                cout << "    take " << i << " " << mean_offset<< endl;
            else 
                cout << "    skip " << i << " " << mean_offset << endl;
        }
    }
}


void KmerAlign3(
        const vec< pair<int,int> > & offset, 
        vec< pair<int,int> > & aligns, int VERBOSITY = 0 )
{
    vec< pair<int,int> >  offset1, offset2, offset3; 
    if ( VERBOSITY >= 1 )
        cout << "step1 " << endl;
    GroupAndSelect( offset, offset1, 200, 0.2, VERBOSITY );
    if ( VERBOSITY >= 1 )
        cout << "step2 " << endl;
    GroupAndSelect( offset1, offset2, 50, 0.1, VERBOSITY );
    if ( VERBOSITY >= 1 )
        cout << "step3 " << endl;
    GroupAndSelect( offset2, offset3, 5, 0.05, VERBOSITY );
    vec< pair<int,int> > aligns_in;
    for ( size_t i = 0; i < offset3.size(); ++i ) 
        aligns_in.push( offset3[i].second -  offset3[i].first,
                 offset3[i].second );
    Sort( aligns_in );
    aligns.clear();
    for ( size_t i = 0; i < aligns_in.size(); ++i ) 
        if ( aligns.empty() ||
             aligns_in[i].first > aligns.back().first &&
             aligns_in[i].second > aligns.back().second )
         aligns.push_back(aligns_in[i] );   
}

void KmerAlignCompare( const vec< pair<int,int> >& aligns1, const vec< pair<int,int> >& aligns2 ) 
{
    cout << "aligns.size()= " << aligns1.size() << " " << aligns2.size() << endl;
    cout << "align1 ";
    for ( size_t i = 0; i < aligns1.size(); ++i ) 
        cout <<  aligns1[i].first << "," << aligns1[i].second << " ";
    cout << endl;
    cout << "align2 ";
    for ( size_t i = 0; i < aligns2.size(); ++i ) 
        cout << aligns2[i].first << "," << aligns2[i].second << " ";
    cout << endl;
}
