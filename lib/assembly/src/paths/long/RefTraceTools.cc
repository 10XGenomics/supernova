///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
#include "paths/long/RefTraceTools.h"
#include "polymorphism/Edit.h"
#include "VecUtilities.h"
#include "math/HoInterval.h"

void CreateHBPlus(const HyperBasevector& hb, const vec<int>& inv,
        HyperBasevector& hbp, vec<pair<int,Bool>>& hbp_to_hb) 
{
     vec< vec<int> > compx;
     hb.Components(compx);
     vec<HyperBasevector> HB;
     vec< pair<int,Bool> > to_hb;
     for ( int i = 0; i < compx.isize( ); i++ )
     {    HB.push( hb, compx[i] );
          for ( int j = 0; j < compx[i].isize( ); j++ )
          {    int v = compx[i][j];
               for ( int l = 0; l < hb.From(v).isize( ); l++ )
                    to_hb.push( hb.EdgeObjectIndexByIndexFrom(v, l), True );    }
          if ( HB.back( ).EdgeObjectCount( ) > 0 
               && inv[ to_hb.back( ).first ] < 0 )
          {    HyperBasevector rc( HB.back( ) );
               rc.Reverse( );
               HB.push_back(rc);    
               for ( int j = 0; j < compx[i].isize( ); j++ )
               {    int v = compx[i][j];
                    for ( int l = 0; l < hb.From(v).isize( ); l++ )
                    {    to_hb.push( hb.EdgeObjectIndexByIndexFrom(v, l), 
                              False );    }    }    }    }
     hbp = HyperBasevector( hb.K(), HB );
     
     // Map edges of hbp back to hb.  An edge in hbp maps to (hb_edge, fw?).
     hbp_to_hb.clear_and_resize( hbp.EdgeObjectCount( ) );
     for ( int v = 0; v < hbp.N( ); v++ )
     {    for ( int j = 0; j < hbp.From(v).isize( ); j++ )
          {    int e = hbp.EdgeObjectIndexByIndexFrom( v, j );
               int ez = to_hb[e].first;
               Bool fw = to_hb[e].second;
               hbp_to_hb[e] = make_pair( ez, fw );    }    }
}

LinearRef::LinearRef(const vec<HyperBasevector>& GH,const vec<bool>& c):isDoubled(c)
{
    if(isDoubled.size()==0) isDoubled.assign(GH.size(),false);
    ForceAssert(GH.size()==isDoubled.size());
    for ( int g = 0; g < (int) GH.size( ); g++ )
    {    const HyperBasevector& x = GH[g];
        if ( !x.Acyclic( ) )
            FatalErr("\nReference sequence must be acyclic.\nAbort.");
        vec<int> sources, sinks;
        x.Sources(sources), x.Sinks(sinks);
        if ( !sources.solo( ) || !sinks.solo( ) )
        {    FatalErr("\nEach reference sequence must have a unique source "
                "and a unique sink.\n"
                "Abort.");    }
        vec< vec<int> > paths;
        // could upgrade to faster version that uses to_left, to_right
        x.EdgePaths( sources[0], sinks[0], paths );
        for ( int j = 0; j < paths.isize( ); j++ )
        {    basevector z = x.EdgeObject( paths[j][0] );
            for ( int i = 1; i < paths[j].isize( ); i++ )
                z = TrimCat( x.K( ), z, x.EdgeObject( paths[j][i] ) );
            G.push_back(z);
            G_source.push_back(g);
            if( isDoubled[g] ){
                if( paths.isize() == 1){
                    G.back().append(z);
                }
                else{
                    FatalErr("Reference "+ToString(g)+" was specified to be circular, but it has more than one component leading to ambiguous definition of assembly error.");
                }
            }
        }
    }
}

void EdgePlacements::Twiddle(const int max_twiddle )
{
     vec< vec<int> > left( hbp.N( ) ), right( hbp.N( ) );
     for ( int i = 0; i < vedata.isize( ); i++ )
     {    int  e = vedata[i].third.first;
          left[ xto_left[e] ].push_back(i), right[ xto_right[e] ].push_back(i);    }
     for ( int v = 0; v < hbp.N( ); v++ )
     {    while(1)
          {    Bool changed = False;
               for ( int j1 = 0; j1 < right[v].isize( ); j1++ )
               for ( int j2 = 0; j2 < left[v].isize( ); j2++ )
               {    int i1 = right[v][j1], i2 = left[v][j2];
                    if ( vedata[i1].first.first != vedata[i2].first.first ) continue;
                    int& x1 = vedata[i1].second.second; 
                    int& x2 = vedata[i2].first.second;
                    if ( x1 == x2 || Abs( x1 - x2 ) > max_twiddle ) continue;

                    /*
                    int e1 = vedata[i1].third.first, e2 = vedata[i2].third.first;
                    cout << "\ninconsistent alignments of edges " << e1
                         << " and " << e2 << "?\n\n";
                    const align &a1 = aligns[i1], &a2 = aligns[i2];
                    int gg = vedata[i1].first.first;
                    PRINT2( a1.pos2( ), a1.Pos2( ) );
                    PrintVisualAlignment( True, cout, hbp.EdgeObject(e1),
                         G[gg], a1 );
                    cout << "\n";
                    PRINT2( a2.pos2( ), a2.Pos2( ) );
                    PrintVisualAlignment( True, cout, hbp.EdgeObject(e2),
                         G[gg], a2 );
                    */

                    int M = Max( x1, x2 );

                    // Avoid cycles.

                    if ( !( vedata[i1].first.second < M ) ) continue;
                    if ( !( M < vedata[i2].second.second ) ) continue;

                    x1 = x2 = M;
                    changed = True;    }
               if ( !changed ) break;    }    }    
}


// For each edge with multiple placements, remove the bad ones.
void EdgePlacements::RemoveBadPlacements()
{
    vec<vec<int>> edge_placements(hbp.EdgeObjectCount());
    for ( int i = 0; i < vedata.isize( ); i++ )
        edge_placements[vedata[i].third.first].push_back(i);
    vec<Bool> todel(vedata.size(), False);
    for (size_t i = 0; i < edge_placements.size(); i++) {
        vec<int>& places = edge_placements[i];
        if (places.solo()) continue;
        auto CompErr = [this](int i, int j)
                    { return vedata[i].third.second < vedata[j].third.second; };
        sort(places.begin(), places.end(), CompErr);
        const int ErrDiffCutOff = 5;
        for (size_t j = 1; j < places.size(); j++) {
            int err_curr = vedata[places[j]].third.second;
            int err_prev = vedata[places[j-1]].third.second;
            if (err_curr > err_prev + ErrDiffCutOff){
                for (size_t k = j; k < places.size(); k++)
                    todel[places[k]] = True;
                break;
            }
        }
    }
    EraseIf(vedata, todel);
    EraseIf(aligns, todel);
}

// Twiddle if the distance is smaller than strong_twiddle, or smaller than
// weak_twiddle && no other alignments ends in middle.
void EdgePlacements::TwiddleSmart()
{
    // location of the vertices
    vec<vec<pair<int,int>>> vlocs(hbp.N());  
    for ( int i = 0; i < vedata.isize( ); i++ ) {    
        int  e = vedata[i].third.first;
        int err = vedata[i].third.second;
        const align& a = aligns[i];
        int g = vedata[i].first.first;
        int head_err = 0, tail_err = 0;
        {
            const int Cutoff = 100;
            int len = hbp.EdgeLength(e);
            int p1 = a.pos1();
            head_err += abs(p1);
            for (int j = 0; j < a.Nblocks(); j++) {
                if (a.Gaps(j) < 0)
                    p1 -= a.Gaps(j);
                if (p1 < Cutoff)
                    head_err += abs(a.Gaps(j));
                if (p1 > len - Cutoff)
                    tail_err += abs(a.Gaps(j));
                p1 += a.Lengths(j);
            }
            tail_err += abs(len - a.Pos1());
        }
        int start = vedata[i].first.second;
        int stop =  vedata[i].second.second;
        vlocs[ xto_left[e] ].push(start - head_err, start + head_err); 
        vlocs[ xto_right[e] ].push(stop - tail_err, stop + tail_err);
    }
    // merge the vertex locations
    for (int v = 0; v < hbp.N(); v++) {
        vec<pair<int,int>>& locs = vlocs[v];
        Sort(locs);
        vec<pair<int,int>> locs_merged;
        for (size_t i = 0; i < locs.size(); i++) {
            int start = locs[i].first;
            int stop = locs[i].second;
            size_t j = i+1;
            while (j < locs.size() && locs[j].first <= stop) {
                stop = max(stop, locs[j].second);
                j++;
            }
            locs_merged.push(start, stop);
            i = j-1;
        }
        swap(locs_merged, locs);
    }
    // Adjust the end locations in vedata to match to vertex locations
    vec< vec<int> > left( hbp.N( ) ), right( hbp.N( ) );
    for ( int i = 0; i < vedata.isize( ); i++ ) {    
        int  e = vedata[i].third.first;
        left[ xto_left[e] ].push_back(i), right[ xto_right[e] ].push_back(i);    
    }
    for ( int v = 0; v < hbp.N( ); v++ ) {
        auto TwiddleLoc = [&](int& pos) {
            for (auto x: vlocs[v]) {
                if (x.first <= pos && pos <= x.second) {
                    pos = (x.first + x.second)/2;
                    return;
                }
            } };
        for ( int j1 = 0; j1 < right[v].isize( ); j1++ ) 
            TwiddleLoc(vedata[right[v][j1]].second.second); 
        for ( int j2 = 0; j2 < left[v].isize( ); j2++ )  
            TwiddleLoc(vedata[left[v][j2]].first.second);
    }
}

basevector EdgePlacements::BestSeq(const vec<int>& best_path, const vec<int>& eids
                                  , const vec<std::pair<int,int>>& limits
                                  , vec<std::tuple<int64_t,int64_t,int,int64_t,int64_t,int64_t,int64_t>>& coors_edge)
{
    ForceAssertGt(best_path.size(), 0u);
    ForceAssertEq(eids.size(), best_path.size());
    int left_trim = vedata[eids.front( )].fourth.first;
    int right_trim = hbp.EdgeObject( best_path.back( ) ).isize( )
                     - vedata[ eids.back( ) ].fourth.second;
    coors_edge.resize(best_path.size());
    size_t end=0;
    coors_edge[0]=std::make_tuple( end,end+hbp.EdgeObject(best_path[0]).size()
                                 ,best_path[0]
                                 ,size_t(0),hbp.EdgeObject(best_path[0]).size()
                                 ,limits[0].first,limits[0].second
                                 );
    end += hbp.EdgeObject(best_path[0]).size();
    basevector best_bpath = hbp.EdgeObject(best_path[0]);
    for ( int i = 1; i < best_path.isize( ); i++ )
    {    best_bpath = TrimCat( hbp.K( ), best_bpath, hbp.EdgeObject( best_path[i] ) );
         std::get<1>(coors_edge[i-1])-=hbp.K()-1;
         std::get<4>(coors_edge[i-1])-=hbp.K()-1;
         end-=hbp.K()-1;
         coors_edge[i]=std::make_tuple(end,end+hbp.EdgeObject(best_path[i]).size()
                                      ,best_path[i]
                                      ,size_t(0),hbp.EdgeObject(best_path[i]).size()
                                      ,limits[i].first,limits[i].second
                                      );
         end += hbp.EdgeObject(best_path[i]).size();
    }

    for(bool bLoop=true; bLoop ;){
        bLoop=false;
        // edge-to-ref alignment can give "negative-direction" alignment
        for(auto&entry: coors_edge){
            const auto& b2 = get<5>(entry);
            auto& e2 = get<6>(entry);
            if(e2<b2){
                e2 = b2+ get<1>(entry) - get<0>(entry);
                bLoop=true;
            }
        }
        // adjacent edges can have overlaping "rough" positions on reference
        for(size_t cc=1;cc<coors_edge.size();++cc){
            auto& b2 = get<5>(coors_edge[cc]);
            const auto& e2 = get<6>(coors_edge[cc-1]);
            if( b2 < e2){
                b2=e2;
                bLoop=true;
            }
        }
    }
    // repair coors_edge according to the left_trim and right_trim
    int64_t nDelete=0;
    for(auto& entry: coors_edge){
        std::get<0>(entry)-=left_trim;
        std::get<1>(entry)-=left_trim;
        if(std::get<1>(entry)<0)     { ++nDelete; }
        else if(std::get<0>(entry)<0){
            std::get<3>(entry)+= -std::get<0>(entry);
            std::get<0>(entry)=0;
        }
    }
    for(size_t ii=0;ii+nDelete<coors_edge.size();++ii){
        coors_edge[ii]=coors_edge[ii+nDelete];
    }
    coors_edge.resize( coors_edge.size()-nDelete);
    int64_t to_trim=right_trim;
    nDelete=0;
    for(auto ritr =coors_edge.rbegin() ;  to_trim > 0 && ritr!=coors_edge.rend(); ++ritr){
        auto loc_size = std::get<1>(*ritr)- std::get<0>(*ritr);
        if(loc_size > to_trim){
            std::get<1>(*ritr) -= to_trim;
            std::get<4>(*ritr) -= to_trim;
            to_trim=0;
        }
        else{
            ++nDelete;
            to_trim-=loc_size;
        }
    }
    coors_edge.resize( coors_edge.size()-nDelete);
    return basevector(best_bpath.begin()+left_trim, best_bpath.end()-right_trim); 
}

int EdgePlacements::CorrelatePositionsAlways( const align& a, const int x1 ) const
{    int pos1 = a.pos1( ), pos2 = a.pos2( );
     if ( x1 < pos1 ) return 0; // off the end, shouldn't happen
     if ( x1 == pos1 ) return pos2;
     int nblocks = a.Nblocks( );
     const avector<int> &gaps = a.Gaps( ), &lengths = a.Lengths( );
     for ( int j = 0; j < nblocks; j++ )
     {    if ( gaps(j) > 0 ) pos2 += gaps(j);
          if ( gaps(j) < 0 )
          {    for ( int x = 0; x < -gaps(j); x++ )
                    if ( x1 == pos1++ ) return pos2;    } // really in gap
          if ( x1 < pos1 + lengths(j) ) return pos2 + x1 - pos1;
          else
          {    pos1 += lengths(j), pos2 += lengths(j);    }    }
     return a.Pos2( );    // off the end, shouldn't happen
} 

void GraphZ::BuildGraph(const int verbosity, ostream& out)
{
     const vec< quad< triple<int,int,int>, triple<int,int,int>, 
          pair<int,int>, pair<int,int> > >& vedata = edge_placements.vedata;
     const vec<align>& aligns = edge_placements.aligns;

     for ( int i = 0; i < vedata.isize( ); i++ )
     {    edges.push( verts.size( ), verts.size( ) + 1, vedata[i].third );
          verts.push_back( vedata[i].first, vedata[i].second );    }
     if ( verbosity >= 2 ) PRINT2_TO( out, verts.size( ), edges.size( ) );
     vec<int> vid( verts.size( ), vec<int>::IDENTITY ), vidi( verts.size( ) );
     SortSync( verts, vid );
     for ( int j = 0; j < verts.isize( ); j++ )
          vidi[ vid[j] ] = j;
     vec<int> v_to_new( verts.size( ), 0 );
     vec<Bool> v_to_delete( verts.size( ), False );
     int vcount = -1;
     for ( int j = 0; j < verts.isize( ); j++ )
     {    if ( j == 0 || verts[j] != verts[j-1] ) vcount++;
          else v_to_delete[j] = True;
          v_to_new[j] = vcount;    }
     EraseIf( verts, v_to_delete );
     for ( int j = 0; j < edges.isize( ); j++ )
     {    edges[j].first = v_to_new[ vidi[ edges[j].first ] ];
          edges[j].second = v_to_new[ vidi[ edges[j].second ] ];    }
     if ( verbosity >= 2 )
     {    out << "\nvertices:\n";
          for ( int j = 0; j < verts.isize( ); j++ )
          {    out << "[" << j << "] " << verts[j].third << " : " 
                    << verts[j].first << "." << verts[j].second << endl;    }
          out << "\nedges:\n";
          for ( int j = 0; j < edges.isize( ); j++ )
          {    int g = verts[ edges[j].first ].first;
               int start = verts[ edges[j].first ].second;
               int stop = verts[ edges[j].second ].second;
               out << "[" << j << "] " << edges[j].first << " --> "
                    << edges[j].second << ", ref: " << g << "." << start << "-" 
                    << stop << ", len: " << stop - start << ", errs: " 
                    << edges[j].third.second << ", assembly edge: " 
                    << hbp_to_hb[ edges[j].third.first ].first << endl;    }    }
     int NN = verts.size( );
     vec< vec<int> > from(NN), to(NN), from_edge_obj(NN), to_edge_obj(NN);
     vec<int> penalty;
     for ( int j = 0; j < edges.isize( ); j++ )
     {    int v = edges[j].first, w = edges[j].second, errs = edges[j].third.second;
          from[v].push_back(w), to[w].push_back(v);
          from_edge_obj[v].push_back( penalty.size( ) );
          to_edge_obj[w].push_back( penalty.size( ) );
          penalty.push_back( Penalty( errs, 0, 0 ) );
          egd.push( errs, 0, 0 );    }
     for ( int v = 0; v < NN; v++ )
     {    SortSync( from[v], from_edge_obj[v] );
          SortSync( to[v], to_edge_obj[v] );    }
     Z.Initialize( from, to, penalty, to_edge_obj, from_edge_obj );    
}

void GraphZ::AddConnectedGapEdges(const int min_dist, const int max_dist, const int verbosity,
        ostream& out ,const bool bPreserveDisconnectedComponents)
{
     vec<size_t> vertex_group(Z.N(),0);
     if(bPreserveDisconnectedComponents){
          vec<vec<int>> vertex_groups;
          Z.Components(vertex_groups);
          for(size_t gg=0;gg<vertex_groups.size();++gg){
              for( const auto& member: vertex_groups[gg]){
                  vertex_group[member]=gg;
              }
          }
     }
     vec<int> sources, sinks;
     Z.Sources(sources), Z.Sinks(sinks);
     if ( verbosity >= 3 )
     {    out << "\ninitial sources:\n";
          for ( int j = 0; j < sources.isize( ); j++ )
          {    int v = sources[j];
               out << "[" << j << "] " << verts[v].first << "."
                    << verts[v].second << " (" << verts[v].third << ")\n";    }
          out << "\ninitial sinks:\n";
          for ( int j = 0; j < sinks.isize( ); j++ )
          {    int v = sinks[j];
               out << "[" << j << "] " << verts[v].first << "."
                    << verts[v].second << " (" << verts[v].third << ")\n";    }    }

     vec< int > z_to_hbp(Z.EdgeObjectCount());
     vec< vec<int> > hbp_to_z(hbp.EdgeObjectCount());
     for(size_t ze=0;ze<edges.size();++ze){
         z_to_hbp[ze] = edges[ze].third.first;
         hbp_to_z[z_to_hbp[ze]].push_back(ze);
     }
     vec<int> hbp_to_left, hbp_to_right;
     hbp.ToLeft(hbp_to_left); hbp.ToRight(hbp_to_right);
     vec<int> z_to_left, z_to_right;
     Z.ToLeft(z_to_left); Z.ToRight(z_to_right);

     for(const auto& z_src:sources){
         for(int jj=0;jj<Z.FromSize(z_src);++jj){
             const int z_src_eid = Z.EdgeObjectIndexByIndexFrom(z_src,jj);
             const int hbp_src_eid = z_to_hbp[z_src_eid];
             const int hbp_v = hbp_to_left[hbp_src_eid];
             for(int kk=0;kk<hbp.ToSize(hbp_v);++kk){
                 const int hbp_prev_eid = hbp.EdgeObjectIndexByIndexTo(hbp_v,kk);
                 for(const auto& z_prev_eid: hbp_to_z[hbp_prev_eid]){
                     int z_prev = z_to_right[z_prev_eid];
                     const bool bConnected = std::any_of(Z.From(z_prev).begin(),Z.From(z_prev).end(),[&](int v){return v==z_src;});
                     if( !bConnected && z_prev != z_src && Z.ToSize(z_prev) > 0){
                        int dist = verts[z_src].second - verts[z_prev].second;
                        if ( dist < min_dist || dist > max_dist ) continue;
                        if ( dist < 0 ) dist = 0;
                        if (!bPreserveDisconnectedComponents || vertex_group[z_prev]==vertex_group[z_src]){
                            Z.AddEdge( z_prev, z_src, Penalty( 0, 1, dist ) );
                            egd.push( 0, 1, dist );
                        }
                     }
                 }
             }
         }
     }
     for(const auto& z_sin:sinks){
         for(int jj=0;jj<Z.ToSize(z_sin);++jj){
             const int z_sin_eid = Z.EdgeObjectIndexByIndexTo(z_sin,jj);
             const int hbp_sin_eid = z_to_hbp[z_sin_eid];
             const int hbp_w = hbp_to_right[hbp_sin_eid];
             for(int kk=0;kk<hbp.FromSize(hbp_w);++kk){
                 const int hbp_next_eid = hbp.EdgeObjectIndexByIndexFrom(hbp_w,kk);
                 for(const auto& z_next_eid: hbp_to_z[hbp_next_eid]){
                     int z_next = z_to_left[z_next_eid];
                     const bool bConnected = std::any_of(Z.From(z_sin).begin(),Z.From(z_sin).end(),[&](int v){return v==z_next;});
                     if( !bConnected && z_next != z_sin && Z.FromSize(z_next) > 0){
                        int dist = verts[z_next].second - verts[z_sin].second;
                        if ( dist < min_dist || dist > max_dist ) continue;
                        if ( dist < 0 ) dist = 0;
                        if (!bPreserveDisconnectedComponents || vertex_group[z_next]==vertex_group[z_sin]){
                            Z.AddEdge( z_sin, z_next, Penalty( 0, 1, dist ) );
                            egd.push( 0, 1, dist );
                        }
                     }
                 }
             }
         }
     }
}

void GraphZ::AddGapEdges(const int min_dist, const int max_dist, const int verbosity,
        ostream& out ,const bool bPreserveDisconnectedComponents)
{
     vec<size_t> vertex_group(Z.N(),0);
     if(bPreserveDisconnectedComponents){
          vec<vec<int>> vertex_groups;
          Z.Components(vertex_groups);
          for(size_t gg=0;gg<vertex_groups.size();++gg){
              for( const auto& member: vertex_groups[gg]){
                  vertex_group[member]=gg;
              }
          }
     }
     vec<int> sources, sinks;
     Z.Sources(sources), Z.Sinks(sinks);
     if ( verbosity >= 3 )
     {    out << "\ninitial sources:\n";
          for ( int j = 0; j < sources.isize( ); j++ )
          {    int v = sources[j];
               out << "[" << j << "] " << verts[v].first << "."
                    << verts[v].second << " (" << verts[v].third << ")\n";    }
          out << "\ninitial sinks:\n";
          for ( int j = 0; j < sinks.isize( ); j++ )
          {    int v = sinks[j];
               out << "[" << j << "] " << verts[v].first << "."
                    << verts[v].second << " (" << verts[v].third << ")\n";    }    }
     vec< triple<int,int,int> > ss_data;
     for ( int j1 = 0; j1 < sources.isize( ); j1++ )
     {    int v = sources[j1];
          ss_data.push( verts[v].first, j1, 1 );    }
     for ( int j2 = 0; j2 < sinks.isize( ); j2++ )
     {    int w = sinks[j2];
          ss_data.push( verts[w].first, j2, 2 );    }
     Sort(ss_data);
     for ( int u = 0; u < ss_data.isize( ); u++ )
     {    int v;
          for ( v = u + 1; v < ss_data.isize( ); v++ )
               if ( ss_data[v].first != ss_data[u].first ) break;
          for ( int x1 = u; x1 < v; x1++ )
          {    if ( ss_data[x1].third != 1 ) continue;
               int j1 = ss_data[x1].second;
               for ( int x2 = u; x2 < v; x2++ )
               {    if ( ss_data[x2].third != 2 ) continue;
                    int j2 = ss_data[x2].second;
                    int v = sinks[j2], w = sources[j1];
                    if ( v == w ) continue;
                    int dist = verts[w].second - verts[v].second;
                    if ( dist < min_dist || dist > max_dist ) continue;
                    if ( dist < 0 ) dist = 0;
                    if (!bPreserveDisconnectedComponents || vertex_group[v]==vertex_group[w]){
                        Z.AddEdge( v, w, Penalty( 0, 1, dist ) );
                        egd.push( 0, 1, dist );
                    }
               }    }
          u = v - 1;    }
}

void GraphZ::AddSpecialVerts( const int K, const vec<int>& sources, const
        vec<int>& sinks, const bool bPreserveDisconnectedComponents)
{
     vec<size_t> vertex_group(Z.N(),0);
     vec<vec<int>> vertex_groups;
     if(bPreserveDisconnectedComponents){
          Z.Components(vertex_groups);
          for(size_t gg=0;gg<vertex_groups.size();++gg){
              for( const auto& member: vertex_groups[gg]){
                  vertex_group[member]=gg;
              }
          }
     }
     else{
         vertex_groups.push_back( vec<int>(Z.N()));
         for(size_t ii=0;ii<vertex_groups.back().size();++ii){ vertex_groups.back()[ii]=ii; }
     }
     int nz = Z.N( );
     const size_t nGroups = vertex_groups.size();
     ForceAssert(nGroups>0);

     Z.AddVertices( 2 * G.size( ) * nGroups);
     for( size_t group=0;group<nGroups;++group){
         for ( int g = 0; g < (int) G.size( ); g++ ){
              vertex_group.push_back(group);
              vertex_groups[group].push_back(verts.size());
              verts.push( g, 0, -1 );
         }
         for ( int g = 0; g < (int) G.size( ); g++ ){
              vertex_group.push_back(group);
              vertex_groups[group].push_back(verts.size());
              verts.push( g, G[g].size( ), -1 );
         }
     }
     for ( int j1 = 0; j1 < sources.isize( ); j1++ )
     {    int gg = verts[ sources[j1] ].first;
          int d = verts[ sources[j1] ].second;
          int g = ( d > 0 ? 1 : 0 );
          int group = vertex_group[sources[j1]];
          if(bPreserveDisconnectedComponents){
              Z.AddEdge( nz + group*2*G.size() + gg, sources[j1], Penalty( d, g, 0 ) );
              egd.push( d, g, 0 );
          }
          else{
              Z.AddEdge( nz + group*2*G.size() + gg, sources[j1], Penalty( 0, g, d ) );
              egd.push( 0, g, d );
          }
     }
     for ( int j2 = 0; j2 < sinks.isize( ); j2++ )
     {    int gg = verts[ sinks[j2] ].first;
          int d = Max( 0, G[gg].isize( ) - verts[ sinks[j2] ].second - (K-1) );
          int g = ( d > 0 ? 1 : 0 );
          int group = vertex_group[sinks[j2]];
          if(bPreserveDisconnectedComponents){
              Z.AddEdge( sinks[j2], nz + group*2*G.size()+(int) G.size( ) + gg, Penalty( d, g, 0 ) );
              egd.push( d, g, 0 );
          }
          else{
              Z.AddEdge( sinks[j2], nz + group*2*G.size()+(int) G.size( ) + gg, Penalty( 0, g, d ) );
              egd.push( 0, g, d );
          }
      }    
}

void GraphZ::FindBestEdgePath( const vec< triple<int,int,int> >& spaths_egd,
        const vec< vec<int> >& spaths,
        vec<vec<int>>& best_path, vec<vec<int>>& eids, int& best_g) 
{
    best_path.clear_and_resize(G.size());
    eids.clear_and_resize( G.size( ) );
    best_g = -1;
    {    vec<int> best_errs( G.size( ), 1000000000 ), best_i( G.size( ), -1 );
         for ( int i = 0; i < spaths.isize( ); i++ )
         {    const vec<int>& p = spaths[i];
              int g = verts[ p.front( ) ].first;
              if ( spaths_egd[i].first < best_errs[g] )
              {    best_errs[g] = spaths_egd[i].first;
                   best_i[g] = i;   
                   best_g = g;    }    }
         for ( int g = 0; g < (int) G.size( ); g++ )
         {    if ( best_errs[g] == Min(best_errs) )
              {    best_g = g;
                   break;    }    }
         for ( int g = 0; g < (int) G.size( ); g++ )
         {    if ( best_i[g] < 0 ) continue;
              const vec<int>& p = spaths[ best_i[g] ];
              for ( int j = 0; j < p.isize( ) - 1; j++ )
              {    int v = p[j], w = p[j+1];
                   int p = 1000000000, best_edge = -1;
                   for ( int l = 0; l < Z.From(v).isize( ); l++ )
                   {    if ( Z.From(v)[l] != w ) continue;
                        int eid = Z.EdgeObjectIndexByIndexFrom( v, l );
                        int el = egd[eid].first, gl = egd[eid].second; 
                        int dl = egd[eid].third;
                        int pl = Penalty( el, gl, dl );
                        if ( pl < p )
                        {    best_edge = l;
                             p = pl;    }    }
                   if ( best_edge < 0 ) continue;
                   int eid = Z.EdgeObjectIndexByIndexFrom( v, best_edge );
                   if ( eid < edges.isize( ) )
                   {    best_path[g].push_back( edges[eid].third.first );
                        eids[g].push_back(eid);    }    }    }    }
}

void GraphZ::FindShortestPath( const int min_dist, const int max_dist,
        vec< vec<int> >& spaths, vec< triple<int,int,int> >& spaths_egd,
        vec< pair<int,int> >& spaths_gg_pen, ostream& out, int verbosity) 
{
    if ( verbosity >= 1 ) out << Date() <<":FindShortestPath()" << endl;
    int K = hbp.K();
    spaths.clear();
    spaths_egd.clear();
    spaths_gg_pen.clear();
    if ( verbosity >= 1 ) out << Date() <<":b4 BuildGraph()" << endl;
    BuildGraph(verbosity, out);
    if ( verbosity >= 1 ) out << Date() <<":after BuildGraph()" << endl;
    if(verbosity>=3){
        std::ofstream ofs("z_init.dot");
        MakeZDot(ofs);
    }

    // Add gap edges.  Then add special vertices for the ends of the chromosomes, 
    // and add gap edges from and to them.  Note that entries in verts refer to 
    // vertex "-1" of hbp.  Note that we deliberately do not recompute sources and 
    // sinks before this step.  That's because the previous step can introduce 
    // cycles.
    vec<int> sources, sinks;
    Z.Sources(sources), Z.Sinks(sinks);
    if ( verbosity >= 1 ) out << Date() <<":b4 AddConnectedGapEdges()" << endl;
    AddConnectedGapEdges(min_dist, max_dist, verbosity, out );
    if(verbosity>=3){
        std::ofstream ofs("z_AddConnectedGapEdges.dot");
        MakeZDot(ofs);
    }
    if ( verbosity >= 1 ) out << Date() <<":after AddConnectedGapEdges()" << endl;

    if ( verbosity >= 1 ) out << Date() <<":b4 AddGapEdges()" << endl;
    AddGapEdges(min_dist, max_dist, verbosity, out );
    if(verbosity>=3){
        std::ofstream ofs("z_AddGapEdges.dot");
        MakeZDot(ofs);
    }
    if ( verbosity >= 1 ) out << Date() <<":after AddGapEdges()" << endl;
    if ( verbosity >= 1 ) out << Date() <<":b4 AddSpecialVerts()" << endl;
    AddSpecialVerts(K, sources, sinks);
    if(verbosity>=3){
        std::ofstream ofs("z_AddSpecialVerts.dot");
        MakeZDot(ofs);
    }
    if ( verbosity >= 1 ) out << Date() <<":after AddSpecialVerts()" << endl;

    // Find shortest paths through the graph and score them.
    if ( verbosity >= 3 ) out << "\nPATHS:\n";
    Z.Sources(sources), Z.Sinks(sinks);
    #pragma omp parallel for
    for ( int j1 = 0; j1 < sources.isize( ); j1++ )
    {    vec<int> suc;
         Z.GetSuccessors1( sources[j1], suc );
         vec<int> sinksx = Intersection( sinks, suc );
         if ( verbosity >= 3 )
         {
              #pragma omp critical
              {    out << Date( ) << ": start source " << j1 << " of "
                        << sources.size( ) << ", successors = " << suc.size( )
                        << " of " << Z.N( ) << endl;    }    }
         digraphE<int> ZS = Z.Subgraph(suc);
         for ( int j2 = 0; j2 < sinksx.isize( ); j2++ )
         {    if ( sources[j1] == sinksx[j2] ) continue;
              FindShortestPathBetween(sources[j1], sinksx[j2], ZS, suc,
                      spaths, spaths_egd, spaths_gg_pen, verbosity, out );    }     }
             
    // Announce paths.
    SortSync( spaths_gg_pen, spaths, spaths_egd );
    AnnouncePaths( spaths, K, spaths_egd, verbosity, out );
    if ( verbosity >= 1 ) out << Date() <<":FindShortestPath() done" << endl;
}

void GraphZ::FindShortestPathBetween( const int this_source, const int this_sink, 
     const digraphE<int>& ZS, const vec<int>& suc, vec< vec<int> >& spaths,
     vec< triple<int,int,int> >& spaths_egd, vec< pair<int,int> >& spaths_gg_pen,
     const int verbosity, ostream& out ) const
{
     vec<int> spath;
     double sclock = WallClockTime( );
     int s1 = BinPosition( suc, this_source ), s2 = BinPosition( suc, this_sink );
     vec<int> pre;
     ZS.GetPredecessors1( s2, pre );
     digraphE<int> ZSP = ZS.Subgraph(pre);
     int p1 = BinPosition(pre, s1), p2 = BinPosition(pre, s2);
     ZSP.ShortestPath( p1, p2, spath );
     for ( int l = 0; l < spath.isize( ); l++ )
          spath[l] = suc[ pre[ spath[l] ] ];
     String stime = TimeSince(sclock);

     // Compute total errors.  Note that this will not always provide the exact
     // correct answer.  An effort is made to prevent overcounting, but this can
     // result in undercounting.

     int e = 0, g = 0, d = 0;
     vec< pair<int,edit0> > edits;
     for ( int i = 0; i < spath.isize( ) - 1; i++ )
     {    int v = spath[i], w = spath[i+1], p = 1000000000;

          // Find the best edge from v to w.

          int best_edge = -1;
          for ( int l = 0; l < Z.From(v).isize( ); l++ )
          {    if ( Z.From(v)[l] != w ) continue;
               int eid = Z.EdgeObjectIndexByIndexFrom( v, l );
               int el = egd[eid].first, gl = egd[eid].second, dl = egd[eid].third;
               int pl = Penalty( el, gl, dl );
               if ( pl < p )
               {    best_edge = l;
                    p = pl;    }    } 
          int eid = Z.EdgeObjectIndexByIndexFrom( v, best_edge );
          if ( eid < edges.isize( ) )
          {    int gg = verts[v].first;
               const align& a = edge_placements.aligns[eid];
               int edge_id = edges[eid].third.first;
               int p1 = a.pos1( ), p2 = a.pos2( );
               const basevector& E = hbp.EdgeObject(edge_id);

               // Traverse the alignment of the edge to locate its errors.

               Bool first = True;
               for ( int j = 0; j < a.Nblocks( ); j++ ) 
               {    if ( a.Gaps(j) > 0 )  
                    {    if ( !( first && edits.nonempty( )
                              && p2 <= edits.back( ).first ) )
                         {    edits.push( p2, edit0( DELETION, a.Gaps(j) ) );    }
                         if ( first && edits.nonempty( )
                              && p2 >= edits.back( ).first )
                         {    first = False;    }
                         p2 += a.Gaps(j);    }
                    if ( a.Gaps(j) < 0 ) 
                    {    if ( !( first && edits.nonempty( )
                              && p2 <= edits.back( ).first ) )
                         {    edits.push( p2, edit0( INSERTION, 
                                   basevector( E, p1, -a.Gaps(j) )
                                   .ToString( ) ) );    }
                         if ( first && edits.nonempty( )
                              && p2 >= edits.back( ).first )
                         {    first = False;    }
                         p1 -= a.Gaps(j);    }
                    for ( int x = 0; x < a.Lengths(j); x++ ) 
                    {    if ( E[p1] != G[gg][p2] )
                         {    if ( !( first && edits.nonempty( )
                                   && p2 <= edits.back( ).first ) )
                              {    edits.push( p2, edit0( SUBSTITUTION, 
                                        (char) as_base(E[p1]) ) );    }
                              if ( first && edits.nonempty( )
                                   && p2 >= edits.back( ).first )
                              {    first = False;    }    }
                         ++p1; ++p2;    }    }    }
          g += egd[eid].second; d += egd[eid].third;    }
     int gg = verts[ spath.front( ) ].first;
     int gstart = verts[ spath.front( ) ].second;
     int gstop = verts[ spath.back( ) ].second;
     int M = -1;
     for ( int i = 0; i < edits.isize( ); i++ )
     {    if ( edits[i].first < M ) continue;
          const edit0& x = edits[i].second;
          M = edits[i].first;
          if ( x.etype == SUBSTITUTION ) e++;
          else if ( x.etype == DELETION ) e += x.n;
          else e += x.seq.size( );    }
     int pen = Penalty( e, g, d );
     #pragma omp critical
     {    if ( verbosity >= 3 ) 
               PRINT8_TO( out, stime, gg, gstart, gstop, pen, e, g, d );
          spaths.push_back(spath);
          spaths_egd.push( e, g, d );
          spaths_gg_pen.push( gg, pen );    }    
}

void GraphZ::AnnouncePaths( const vec< vec<int> >& spaths, const int K, 
         const vec< triple<int,int,int> >& spaths_egd, 
         const int verbosity, ostream& out ) const
{
     if ( verbosity >= 1 )
     {    out << "\nThere are " << spaths.size( ) << " paths:\n";
          for ( int i = 0; i < spaths.isize( ); i++ )
          {    const vec<int>& p = spaths[i];
               int gg = verts[ p.front( ) ].first;
               out << "[" << i+1 << "] " << gg << "." << verts[ p.front( ) ].second
                    << " (" << verts[ p.front( ) ].third << ") --> " << gg << "." 
                    << verts[ p.back( ) ].second << " (" 
                    << verts[ p.back( ) ].third << ")" << " [errs=" 
                    << spaths_egd[i].first << ",gaps=" << spaths_egd[i].second 
                    << ",gap_bases=" << spaths_egd[i].third << "]" << endl;    }    }
     if ( verbosity >= 3 )
     {    out << "\nDetailed path info:\n";
          for ( int i = 0; i < spaths.isize( ); i++ )
          {    const vec<int>& p = spaths[i];
               int gg = verts[ p.front( ) ].first;
               out << "\n[" << i+1 << "] " << gg << "." 
                    << verts[ p.front( ) ].second << " (" 
                    << verts[ p.front( ) ].third << ") --> " << gg << "." 
                    << verts[ p.back( ) ].second << " (" 
                    << verts[ p.back( ) ].third << ")" << " [errs=" 
                    << spaths_egd[i].first << ",gaps=" << spaths_egd[i].second 
                    << ",gap_bases=" << spaths_egd[i].third << "]\n\n";
               for ( int j = 0; j < p.isize( ) - 1; j++ )
               {    int v = p[j], w = p[j+1];
                    out << "   " << gg << "." << verts[ p[j] ].second << " (" 
                         << verts[ p[j] ].third << ") --> " << gg << "." 
                         << verts[ p[j+1] ].second
                              + ( j == p.isize( ) - 2 ? (K-1) : 0 )
                         << " (" << verts[ p[j+1] ].third << ")";
                    int p = 1000000000, best_edge = -1;
                    for ( int l = 0; l < Z.From(v).isize( ); l++ )
                    {    if ( Z.From(v)[l] != w ) continue;
                         int eid = Z.EdgeObjectIndexByIndexFrom( v, l );
                         int el = egd[eid].first, gl = egd[eid].second; 
                         int dl = egd[eid].third;
                         int pl = Penalty( el, gl, dl );
                         if ( pl < p )
                         {    best_edge = l;
                              p = pl;    }    }
                    if ( best_edge < 0 )
                    {    out << " -- DISCONTINUITY\n";
                         continue;    }
                    int eid = Z.EdgeObjectIndexByIndexFrom( v, best_edge );
                    String hbp_edge_id = "*", hb_edge_id = "*";
                    if ( eid < edges.isize( ) )
                    {    int x = edges[eid].third.first;
                         hbp_edge_id = ToString(x);
                         hb_edge_id = ToString( hbp_to_hb[x].first )
                              + ( hbp_to_hb[x].second ? "fw" : "rc" );    }
                    out << " [hbp_edge=" << hbp_edge_id << ",hb_edge=" 
                         << hb_edge_id << ",e=" << egd[eid].first << ",g=" 
                         << egd[eid].second << ",d=" << egd[eid].third << "]" 
                         << endl;    }    }    }    
}
void GraphZ::MakeZDot(ostream& out){
    vec<String> vl(Z.N());
    vec<vec<String>> el(Z.N());
    for(int64_t vv=0;vv<Z.N();++vv){
        vl[vv]=ToString(vv);
        if(vv<verts.isize()){ vl[vv]+=": " + ToString(verts[vv].second); }
        else                { vl[vv]+=": *";}
        el[vv].resize(Z.FromSize(vv));
        for(int64_t ee=0;ee<Z.FromSize(vv);++ee){
            int eid = Z.EdgeObjectIndexByIndexFrom(vv,ee);
            el[vv][ee]=ToString(eid);
            if(eid < edges.isize()){
                int hbp_e = edges[ eid ].third.first ;
                el[vv][ee]+=": "+ToString(hbp_e) + " " + ToString( hbp_to_hb[ hbp_e ].first);
            }
            else                   { el[vv][ee]+=": *";};
            el[vv][ee] += "("+ToString(Z.EdgeObjectByIndexFrom(vv,ee))+")";
        }
    }
    Z.DOT_vl(out,vl,"",vec<vec<String>>(),vec<String>(),el);
}
