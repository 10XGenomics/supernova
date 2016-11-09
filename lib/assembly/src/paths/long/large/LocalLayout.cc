///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "paths/long/large/LocalLayout.h"
#include "Basevector.h"
#include "feudal/BaseVec.h"
#include "CoreTools.h"
#include "graph/Digraph.h"
#include "paths/HyperBasevector.h"
#include "paths/KmerBaseBroker.h"
#include "paths/long/LongReadsToPaths.h"
#include "paths/long/MakeKmerStuff.h"
#include "CoreTools.h"



namespace {

template <class Itr>
int ScoreKmer( Itr begin, Itr end )
{
    int qsum = 0, qavg = 0;
    size_t count = 0;
//    for ( auto itr = begin; itr != end; ++itr )
//	if ( *itr > 2 ) { qavg += *itr; count++; }
//    if ( count > 0 ) qavg /= count;
//    else qavg = 20;
    qavg = 20;
    for ( auto itr = begin; itr != end; ++itr )
	qsum += ( *itr ? *itr : qavg );

    return qsum;
}

longlong KmerQScoreSum(int K, kmer_id_t kid, vecbasevector const& bases,
	vecqualvector const& quals, vecKmerPath const& paths,
	vecKmerPath const& paths_rc, vec<big_tagged_rpint> const& pathsdb, ostream& out )
{
    longlong qsumsum = 0;
    bool verbose = true;		// write logging to out

    vec<longlong> path_ids;
    Contains(pathsdb, kid, path_ids);

    for ( auto const& path_id : path_ids ) {
	ForceAssertGe(path_id,0);

	big_tagged_rpint match = pathsdb[path_id];
	longlong read_id = match.ReadId();
	bool rc = match.Rc();
	const KmerPath& the_path = ( rc ? paths_rc[read_id] : paths[read_id] );

	int offset = kid - match.Start();
	for(int i=0; i<match.PathPos(); i++)
	    offset += the_path.Length(i);


	longlong qsum;
	if ( rc ) {
	    auto rcbegin = quals[read_id].rbegin(offset);
	    auto rcend = rcbegin + K;
	    qsum = ScoreKmer( rcbegin, rcend  );
	} else {
	    auto fwbegin = quals[read_id].begin(offset);
	    auto fwend = fwbegin + K;
	    qsum = ScoreKmer( fwbegin, fwend );
	}

	// begin temporary
	qualvector kmer_quals;
	if(rc) {
	    kmer_quals.SetToSubOf( quals[read_id], quals[read_id].size()-offset-K , K );
	    kmer_quals.ReverseMe();
	} else {
	    kmer_quals.SetToSubOf( quals[read_id], offset, K);
	}
	longlong qsum0 = ScoreKmer( kmer_quals.begin(), kmer_quals.end() );
	ForceAssertEq(qsum,qsum0);
	// end temporary

	if ( verbose ) {
	    basevector kmer_bases;
	    if ( rc ) {
		kmer_bases.SetToSubOf( bases[read_id], bases[read_id].size()-offset-K , K );
		kmer_bases.ReverseComplement();
	    } else {
		kmer_bases.SetToSubOf( bases[read_id], offset, K);
	    }

	    out << "    ";
	    kmer_bases.Print(out);
            out << "    ";
            for ( auto const& qual : kmer_quals )
                out << static_cast<int>( qual / 10 );
            out << endl << "    ";
            for ( auto const& qual : kmer_quals )
                out << static_cast<int>( qual % 10 );
            out << endl;


	    out << "    [" << qsum << "] ";
	    out << endl;
	}

	qsumsum += qsum;

    }
    out << "  [TOTAL " << qsumsum << " AVG " << static_cast<double>(qsumsum) / path_ids.size() << "]" << endl;

    return qsumsum;
}

};	// end of anonymous namespace



void LocalLayout( const int lroot, const int rroot, const LongProtoTmpDirManager& dir_mgr,
     const String& work_dir )
{    
//     vecbasevector cbases( TMP + "/frag_reads_mod0.fastb" );
//     vecqualvector cquals( TMP + "/frag_reads_mod0.qualb" );
     vecbasevector const& cbases=dir_mgr.get("frag_reads_mod0").reads();
     vecqualvector const& cquals=dir_mgr.get("frag_reads_mod0").quals();

     const int K = 60;

     HyperBasevector hb;
     HyperKmerPath h;
     vecKmerPath paths, paths_rc;
     vec<big_tagged_rpint> pathsdb;
     unsigned const COVERAGE = 50u; 
     LongReadsToPaths(cbases,K,COVERAGE,0,False, &hb,&h,&paths,&paths_rc,&pathsdb);

     String head = work_dir + "/loc/uni." + ToString(lroot) + "." + ToString(rroot);
     Ofstream( bout, head + ".fasta" );
     for ( int e = 0; e < hb.EdgeObjectCount( ); e++ )
          hb.EdgeObject(e).Print(bout,e);
     BinaryWriter::writeFile( head + ".hbv", hb );

     Ofstream( out2, head + ".reads" );
     for ( int id = 0; id < h.EdgeObjectCount(); ++id ) {
	 auto const& kmer_path = h.EdgeObject(id);
	 size_t edge_len = hb.EdgeObject(id).size();
	 ForceAssertEq(kmer_path.NSegments(),1);
	 ForceAssert(kmer_path.Segment(0).isSeq());
	 auto const& kmer_path_interval = kmer_path.Segment(0);
	 vec<longlong> path_ids;
	 Contains(pathsdb, kmer_path_interval, path_ids );
	 out2 << "ID " << id << "(" << edge_len << "): ";
	 for ( auto const& path_id : path_ids )
	     out2 << path_id << " ";
	 out2 << endl;
     }

     Ofstream( out3, head + ".junctions" );
     vec<Bool> delete_me_from( h.EdgeObjectCount(), False );
     vec<Bool> delete_me_to( h.EdgeObjectCount(), False );
     map<int,int> delete_from;
     map<int,int> delete_to;
     for ( int vid = 0; vid < h.N(); ++vid ) {
	 out3 << "FROM " << vid << ": " << endl;
	 vec<double> scores;
	 for ( int fid = 0; fid < h.FromSize(vid); ++fid ) {
	     int edge_id = h.FromEdgeObj(vid)[fid];
	     kmer_id_t kid = h.EdgeObjectByIndexFrom(vid,fid).Start();
	     out3 << "  EDGE " << edge_id << ":" << endl;
	     scores.push_back( KmerQScoreSum(K, kid, cbases, cquals, paths,
			 paths_rc, pathsdb, out3 ) );
	 }
	 double mmax = scores.size() == 0 ? 0. : Max(scores);
	 double mavg = scores.size() == 0 ? 0. : (
		 std::accumulate( scores.begin(), scores.end(), 0.) /
		 static_cast<double>(scores.size()) );
	 for ( int fid = 0; fid < h.FromSize(vid); ++fid ) {
	     int edge_id = h.FromEdgeObj(vid)[fid];
	     bool test1 = scores[fid] < mmax*0.2;
	     bool test2 = scores[fid] < mavg*0.6;

	     if ( test1 || test2 ) {
		 out3 << "  REMOVE EDGE " << edge_id << " test1=" << (test1?"T":"F")
			 << " test2=" << (test2?"T":"F") << endl;
		 delete_me_from[edge_id] = True;
		 delete_from.insert(make_pair(edge_id,vid));
	     }
	 }
	 out3 << endl;

	 out3 << "TO " << vid << ": " << endl;
	 scores.clear();
	 for ( int fid = 0; fid < h.ToSize(vid); ++fid ) {
	     int edge_id = h.ToEdgeObj(vid)[fid];
	     kmer_id_t kid = h.EdgeObjectByIndexTo(vid,fid).Stop();
	     out3 << "  EDGE " << edge_id << ":" << endl;
	     scores.push_back( KmerQScoreSum(K, kid, cbases, cquals, paths,
			 paths_rc, pathsdb, out3 ) );
	 }
	 mmax = scores.size() == 0 ? 0. : Max(scores);
	 mavg = scores.size() == 0 ? 0. : (
		 std::accumulate( scores.begin(), scores.end(), 0.) /
		 static_cast<double>(scores.size()) );
	 for ( int fid = 0; fid < h.ToSize(vid); ++ fid ) {
	     int edge_id = h.ToEdgeObj(vid)[fid];
	     bool test1 = scores[fid] < mmax*0.2;
	     bool test2 = scores[fid] < mavg*0.6;

	     if ( test1 || test2 ) {
		 out3 << "  REMOVE EDGE " << edge_id << " test1=" << (test1?"T":"F")
			 << " test2=" << (test2?"T":"F") << endl;
		 delete_me_to[edge_id] = True;
		 delete_to.insert( make_pair(edge_id,vid) );
	     }
	 }
	 out3 << endl;
     }

     vec<String> graph_colors( delete_me_from.size() );
     for ( size_t eid = 0; eid < graph_colors.size(); ++eid )
	 if ( delete_me_from[eid] || delete_me_to[eid] ) graph_colors[eid] = "red";

     Ofstream( out4, head + ".dot" );
     hb.PrintSummaryDOT0w( out4, True, True, True, NULL, False, NULL,
	     NULL, NULL, &graph_colors );

//     vec<int> to_delete0;
//     for ( size_t eid = 0; eid < delete_me.size(); ++eid )
//	 if ( delete_me[eid] ) to_delete0.push_back(eid);
//     hb.DeleteEdges(to_delete0);

     vec<int> to_delete_full;
     {
     vec<int> to_left, to_right;
     hb.ToLeft( to_left );
     hb.ToRight( to_right );
     for ( int eid = 0; eid < hb.EdgeObjectCount(); ++eid )
	 if (  delete_me_from[eid] && delete_me_to[eid] ||
		 delete_me_from[eid] && hb.FromSize(to_right[eid]) == 0 ||
		 delete_me_to[eid] && hb.ToSize(to_left[eid]) == 0 ) {
	     to_delete_full.push_back( eid );
	     delete_from.erase(eid);
	     delete_to.erase(eid);
	 }
     }

     hb.DeleteEdges(to_delete_full);

     size_t nadd = delete_from.size() + delete_to.size();
     size_t new_vertex = hb.N();
     hb.AddVertices(nadd);

     for ( auto from : delete_from )
	 hb.GiveEdgeNewFromVx(from.first, from.second, new_vertex++);

     for ( auto to : delete_to )
	 hb.GiveEdgeNewToVx(to.first, to.second, new_vertex++);

     hb.RemoveEdgelessVertices();		// don't remove dead edge objects yet to avoid renumbering

     Ofstream( out5, work_dir + "/loc/uni."
          + ToString(lroot) + "." + ToString(rroot) + ".edited.dot" );
     hb.PrintSummaryDOT0w( out5, True, True, True );


     int norig = cbases.size( );
     // cbases.Append(injections);
     vecbasevector cbasesrc(cbases);
     for ( int i = 0; i < (int) cbasesrc.size( ); i++ )
          cbasesrc[i].ReverseComplement( );
     vecbasevector cbases2(cbases);
     cbases2.Append(cbasesrc);
     vec< triple<kmer<K>,int,int> > kmers_plus;
     MakeKmerLookup3( cbases2, kmers_plus );
     vec<Bool> aligned( cbases.size( ), False );
     vec< triple<int,int,int> > aligns;
     for ( int i = 0; i < kmers_plus.isize( ); i++ )
     {    int j;
          for ( j = i + 1; j < kmers_plus.isize( ); j++ )
               if ( kmers_plus[j].first != kmers_plus[i].first ) break;

          /*
          Bool twice = False;
          for ( int k = i; k < j; k++ )
          {    int l;
               for ( l = k + 1; l < j; l++ )
                    if ( kmers_plus[l].second != kmers_plus[k].second ) break;
               if ( l - k >= 2 ) twice = True;
               k = l - 1;    }
          if (twice)
          {    i = j - 1;
               continue;    }
          */

          for ( int k1 = i; k1 < j; k1++ )
          for ( int k2 = i; k2 < j; k2++ )
          {    if ( k2 == k1 ) continue;
               int id1 = kmers_plus[k1].second, id2 = kmers_plus[k2].second;
               int pos1 = kmers_plus[k1].third, pos2 = kmers_plus[k2].third;
               int offset = pos1 - pos2;
               if ( offset < 0 ) continue;
               if ( offset == 0 && id1 > id2 ) continue;
               int id1x(id1), id2x(id2);
               if ( id1x >= (int) cbases.size( ) ) id1x -= cbases.size( );
               if ( id2x >= (int) cbases.size( ) ) id2x -= cbases.size( );
               aligned[id1x] = aligned[id2x] = True;
               aligns.push( id1, id2, offset );    }
          i = j - 1;    }



     {
     const int K = 40;
     vec< triple<kmer<K>,int,int> > kmers_plus;
     MakeKmerLookup3( cbases2, kmers_plus );
     for ( int i = 0; i < kmers_plus.isize( ); i++ )
     {    int j;
          for ( j = i + 1; j < kmers_plus.isize( ); j++ )
               if ( kmers_plus[j].first != kmers_plus[i].first ) break;

          /*
          Bool twice = False;
          for ( int k = i; k < j; k++ )
          {    int l;
               for ( l = k + 1; l < j; l++ )
                    if ( kmers_plus[l].second != kmers_plus[k].second ) break;
               if ( l - k >= 2 ) twice = True;
               k = l - 1;    }
          if (twice)
          {    i = j - 1;
               continue;    }
          */

          for ( int k1 = i; k1 < j; k1++ )
          for ( int k2 = i; k2 < j; k2++ )
          {    if ( k2 == k1 ) continue;
               int id1 = kmers_plus[k1].second, id2 = kmers_plus[k2].second;
               int pos1 = kmers_plus[k1].third, pos2 = kmers_plus[k2].third;
               int offset = pos1 - pos2;
               if ( offset < 0 ) continue;
               if ( offset == 0 && id1 > id2 ) continue;
               int id1x(id1), id2x(id2);
               if ( id1x >= (int) cbases.size( ) ) id1x -= cbases.size( );
               if ( id2x >= (int) cbases.size( ) ) id2x -= cbases.size( );
               if ( aligned[id1x] && aligned[id2x] ) continue;
               aligns.push( id1, id2, offset );    }
          i = j - 1;    }
     }




     UniqueSort(aligns);
     int nb2 = 2*cbases.size( );
     vec<vec<int>> from(nb2), to(nb2), from_edge_obj(nb2), to_edge_obj(nb2);
     vec<int> edges;
     for ( int i = 0; i < aligns.isize( ); i++ )
     {    int id1 = aligns[i].first, id2 = aligns[i].second;
          int offset = aligns[i].third;
          from[id1].push_back(id2), to[id2].push_back(id1);
          from_edge_obj[id1].push_back( edges.size( ) );
          to_edge_obj[id2].push_back( edges.size( ) );
          edges.push_back(offset);    }
     for ( int v = 0; v < nb2; v++ )
     {    SortSync( from[v], from_edge_obj[v] );
          SortSync( to[v], to_edge_obj[v] );    }
     digraphE<int> G( from, to, edges, to_edge_obj, from_edge_obj );
     vec<Bool> to_delete( G.EdgeObjectCount( ), False );
     
     for ( int x = 0; x < G.N( ); x++ )
     for ( int i = 0; i < G.From(x).isize( ); i++ )
     {    int y = G.From(x)[i];
          int e = G.EdgeObjectIndexByIndexFrom( x, i );
          if ( to_delete[e] ) continue;
          int ne = G.EdgeObject(e);
          if ( ne < 0 ) continue;
          for ( int j = 0; j < G.From(y).isize( ); j++ )
          {    int z = G.From(y)[j];
               int f = G.EdgeObjectIndexByIndexFrom( y, j );
               if ( to_delete[f] ) continue;
               int nf = G.EdgeObject(f);
               if ( nf < 0 ) continue;
               for ( int k = 0; k < G.From(x).isize( ); k++ )
               {    if ( G.From(x)[k] != z ) continue;
                    int g = G.EdgeObjectIndexByIndexFrom( x, k );
                    if ( to_delete[g] ) continue;
                    int ng = G.EdgeObject(g);
                    if ( ne + nf != ng ) continue;
                    to_delete[g] = True;    }    }    }

     for ( int x = 0; x < G.N( ); x++ )
     for ( int i = 0; i < G.From(x).isize( ); i++ )
     {    int y = G.From(x)[i];
          int e = G.EdgeObjectIndexByIndexFrom( x, i );
          if ( to_delete[e] ) continue;
          int ne = G.EdgeObject(e);
          if ( ne < 0 ) continue;
          for ( int j = 0; j < G.From(y).isize( ); j++ )
          {    int z = G.From(y)[j];
               int f = G.EdgeObjectIndexByIndexFrom( y, j );
               if ( to_delete[f] ) continue;
               int nf = G.EdgeObject(f);
               if ( nf < 0 ) continue;
               for ( int k = 0; k < G.From(z).isize( ); k++ )
               {    int w = G.From(z)[k];
                    int g = G.EdgeObjectIndexByIndexFrom( z, k );
                    if ( to_delete[g] ) continue;
                    int ng = G.EdgeObject(g);
                    if ( ng < 0 ) continue;
                    for ( int l = 0; l < G.From(x).isize( ); l++ )
                    {    if ( G.From(x)[l] != w ) continue;
                         int h = G.EdgeObjectIndexByIndexFrom( x, l );
                         if ( to_delete[h] ) continue;
                         int nh = G.EdgeObject(h);
                         if ( ne + nf + ng != nh ) continue;
                         to_delete[h] = True;    }    }    }    }

     const int max_max_depth = 30;
     const int max_chains = 100;
     vec<int> to_left, to_right;
     G.ToLeft(to_left), G.ToRight(to_right);
     for ( int max_depth = 3; max_depth <= max_max_depth; max_depth++ )
     {    for ( int x = 0; x < G.N( ); x++ )
          for ( int i = 0; i < G.From(x).isize( ); i++ )
          {    int y = G.From(x)[i];
               int e = G.EdgeObjectIndexByIndexFrom( x, i ); // goal is to delete e
               if ( to_delete[e] ) continue;
               int ne = G.EdgeObject(e);
               vec<vec<int>> chains;
               vec<int> p;
               chains.push_back(p);
               for ( int d = 0; d < max_depth; d++ )
               {    if ( chains.isize( ) > max_chains ) break;
                    vec<vec<int>> chainsx;
                    for ( int c = 0; c < chains.isize( ); c++ )
                    {    vec<int> q = chains[c];
                         int v;
                         if ( chains[c].empty( ) ) v = x;
                         else v = to_right[ chains[c].back( ) ];
                         for ( int j = 0; j < G.From(v).isize( ); j++ )
                         {    vec<int> q = chains[c];
                              int f = G.EdgeObjectIndexByIndexFrom( v, j );
                              if ( f == e || to_delete[f] ) continue;
                              q.push_back(f);
                              if ( to_right[f] == y )
                              {    int ext = 0;
                                   for ( int l = 0; l < q.isize( ); l++ )
                                        ext += G.EdgeObject( q[l] );
                                   if ( ext == ne )
                                   {    to_delete[e] = True;
                                        goto tail;    }    }    
                              chainsx.push_back(q);    }    }
                    chains = chainsx;    }
               tail: continue;    }    }

     vec<int> dels;
     for ( int e = 0; e < G.EdgeObjectCount( ); e++ )
          if ( to_delete[e] ) dels.push_back(e);
     G.DeleteEdges(dels);

     int bads = 0;
     for ( int v = 0; v < G.N( ); v++ )
          if ( !G.From(v).solo( ) || !G.To(v).solo( ) ) bads++;
     bads = (bads-4)/2;
          
     Ofstream( out, work_dir + "/loc/layout."
          + ToString(lroot) + "." + ToString(rroot) + ".dot" );
     out << "// using " << dir_mgr.dir() << "\n";
     out << "// badness = " << bads << "\n";
     vec<String> vertex_labels(nb2);
     for ( int i = 0; i < (int) cbases.size( ); i++ )
     {    if ( i >= norig )
          {    vertex_labels[i] = "+I" + ToString(i-norig);
               vertex_labels[ i + (int) cbases.size( ) ]     
                    = "-I" + ToString(i-norig);    }
          else
          {    vertex_labels[i] = "+" + ToString(i);
               vertex_labels[ i + (int) cbases.size( ) ] = "-" + ToString(i);    }    }
     G.DOT_vel( out, vertex_labels );
}
