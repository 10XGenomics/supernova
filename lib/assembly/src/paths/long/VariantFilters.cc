///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
#include "paths/long/Variants.h"
#include "paths/long/VariantFilters.h"
#include "paths/long/CreateGenome.h"
#include "paths/HyperBasevector.h"
#include "PackAlign.h"
#include "Basevector.h"

void RemoveRepetitiveEdges(vec<VariantCallGroup>& groups, 
        const vec<size_t>& ref_index, 
        const vec<size_t>& ref_shift, 
        const String& ref_m100_file,
        const HyperBasevector& hbp,
        const vec<pair<int,Bool>>& hbp_to_hb,
        int verbosity)
{
    if (verbosity >= 1)
        cout << Date() << ": Removing variants on repetitive edges using genome mask file "
            << ref_m100_file  << endl;
    const int maxq = 600;
    vecbitvector amb(ref_m100_file);
    for (VariantCallGroup& group: groups) {
        int gid = group.GetGID();
        size_t chr_id = ref_index[gid];
        int start = ref_shift[gid];
        for (int i = 0; i < group.GetNBranch(); i++) {
            const vec<int>& path0 = group.GetBranchPath(i);
            int edge_len = hbp.EdgePathToBases(path0).size();
            //if (path0.size() != 1) continue;
            if (edge_len > maxq) continue;

            const align& a = group.GetAlign(i);
            pair<int,int> trims = group.GetTrim(i);
            int estart = start + a.pos2( ) - trims.first;
            int estop = start + a.Pos2( ) + trims.second;

            vec<int> assembly_edge_ids(path0.size());
            for (size_t i = 0; i < path0.size(); i++) {
                assembly_edge_ids[i] = hbp_to_hb[path0[i]].first;
            }
            if (verbosity >= 2) {
                cout << "edge "; assembly_edge_ids.Print(cout);
                cout << " path0= "; path0.Print(cout);
                cout << " estart= " << estart << " estop= " << estop << endl;
            }

            Bool unique = False;
            const int K = 100;
            for ( int j = estart; j <= estop - K; j++ )
                if ( !amb[chr_id][j] ) unique = True;
            if (unique) continue;

            if (verbosity >= 1) {
                cout << "    Removing " << group.GetVariantCalls(i).size() 
                    << " variants from edge ";
                assembly_edge_ids.Println(cout);
            }
            group.GetVariantCalls(i).clear();
        }
    }
}

