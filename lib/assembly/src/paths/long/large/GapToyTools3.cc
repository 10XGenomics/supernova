///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"
#include "ParallelVecUtilities.h"
#include "PrintAlignment.h"
#include "Qualvector.h"
#include "graphics/BasicGraphics.h"
#include "kmers/LongReadPather.h"
#include "paths/HyperBasevector.h"
#include "paths/ReadsToPathsCoreX.h"
#include "paths/RemodelGapTools.h"
#include "paths/Unipath.h"
#include "paths/long/Correct1Pre.h"
#include "paths/long/FriendAligns.h"
#include "paths/long/HBVFromEdges.h"
#include "paths/long/LongReadsToPaths.h"
#include "paths/long/ReadPath.h"
#include "paths/long/ReadPathTools.h"
#include "paths/long/ReadStack.h"
#include "paths/long/large/GapToyTools.h"

void SelectSpecials( const HyperBasevector& hb, vecbasevector& bases,
      VecPQVec const& quals, const ReadPathVec& paths2, const String& work_dir )
{    
     // Heuristics.

     const int min_qsum = 1000;
     const int min_q = 3;
     const int M = 32;
     const int min_adv = 20;
     const int max_to_call_good = 100;

     vec<int> to_left, to_right;
     hb.ToLeft(to_left), hb.ToRight(to_right);
     vecbasevector all;
     int nedges = hb.EdgeObjectCount( );
     all.reserve(nedges);
     for ( int e = 0; e < nedges; e++ )
          all.push_back( hb.EdgeObject(e) );
     vec< triple<kmer<M>,int,int> > kmers_plus;
     MakeKmerLookup0( all, kmers_plus );
     vec< kmer<M> > kmers( kmers_plus.size( ) );
     for ( int i = 0; i < (int) kmers_plus.size( ); i++ )
          kmers[i] = kmers_plus[i].first;
     int printed = 0;
     vec<int> select;
     auto qvItr = quals.begin();
     for ( int64_t id = 0; id < (int64_t) paths2.size( ); id++,++qvItr )
     {    Bool placed = False;
          int start = 0;
          const ReadPath& p = paths2[id];
          if ( paths2[id].size( ) > 0 )
          {    start = hb.EdgeLengthBases( p[0] );
               for ( int j = 1; j < (int) p.size( ); j++ )
                    start += hb.EdgeLengthKmers( p[j] );    }
          int qsum = 0;
          qvec const& qv = *qvItr;
          for ( int j = start; j < (int) qv.size( ); j++ )
               if ( qv[j] >= min_q ) qsum += qv[j];
          if ( qsum < min_qsum ) continue;

          if ( p.size( ) == 0 )
          {    vec< pair<int,int> > places;
               for ( int j = 0; j <= bases[id].isize( ) - M; j++ )
               {    kmer<M> x;
                    x.SetToSubOf( bases[id], j );
                    int low = LowerBound( kmers, x );
                    int high = UpperBound( kmers, x );
                    for ( int k = low; k < high; k++ )
                    {    int e = kmers_plus[k].second;
                         int offset = kmers_plus[k].third - j;
                         places.push( e, offset );    }    }
               UniqueSort(places);

               vec< pair<int,int> > starts;
               vec< pair<int,int> > left_termini;
               for ( int j = 0; j < places.isize( ); j++ )
               {    int e = places[j].first, offset = places[j].second;
                    int v = to_left[e];
                    if ( offset >= 0 ) starts.push_back( places[j] );
                    else if ( hb.To(v).empty( ) )
                    {    left_termini.push_back( places[j] );    }
                    else
                    {    for ( int l = 0; l < hb.To(v).isize( ); l++ )
                         {    int ep = hb.EdgeObjectIndexByIndexTo( v, l );
                              int offsetp = offset + hb.EdgeLengthKmers(ep);
                              if ( !Member( places, make_pair( ep, offsetp ) ) )
                                   places.push( ep, offsetp );    }    }    }

               vec<Bool> to_delete( starts.size( ), False );
               for ( int j = 0; j < starts.isize( ); j++ )
               {    int e = starts[j].first, offset = starts[j].second;
                    if ( hb.EdgeLengthBases(e) - offset <= hb.K( ) - 1 )
                    {    if ( hb.From( to_right[e] ).nonempty( ) ) 
                              to_delete[j] = True;    }    }
               EraseIf( starts, to_delete );

               int m = bases[id].size( );
               for ( int j = 0; j < starts.isize( ); j++ )
               {    int e = starts[j].first, offset = starts[j].second;
                    m = Min( m, hb.EdgeLengthBases(e) - offset );    }

               vec<int> qsum_starts( starts.size( ), 0 );
               for ( int j = 0; j < starts.isize( ); j++ )
               {    int e = starts[j].first, offset = starts[j].second;
                    for ( int l = 0; l < m; l++ )
                    {    if ( bases[id][l] != hb.EdgeObject(e)[offset+l] )
                              qsum_starts[j] += qv[l];    }    }
               SortSync( qsum_starts, starts );
               for ( int j = 0; j < starts.isize( ); j++ )
               {    if ( j > 0 && qsum_starts[j] - qsum_starts[j-1] >= min_adv )
                    {    starts.resize(j);
                         qsum_starts.resize(j);    }    }

               vec< vec< pair<int,int> > > full_paths;
               vec< vec< pair<int,int> > > right_termini;
               for ( int j = 0; j < starts.isize( ); j++ )
               {    vec< vec< pair<int,int> > > partials;
                    vec< pair<int,int> > ini = { starts[j] };
                    partials.push_back(ini);
                    while( partials.nonempty( ) )
                    {    vec< vec< pair<int,int> > > partials2;
                         for ( int l = 0; l < partials.isize( ); l++ )
                         {    vec< pair<int,int> > p = partials[l];
                              int e = p.back( ).first, offset = p.back( ).second;
                              int stop = offset + bases[id].isize( );
                              int w = to_right[e];
                              if ( stop <= hb.EdgeLengthBases(e) )
                                   full_paths.push_back(p);
                              else if ( hb.From(w).empty( ) )
                              {    right_termini.push_back(p);    }
                              else
                              {    for ( int i = 0; i < hb.From(w).isize( ); i++ )
                                   {    int ep = hb.EdgeObjectIndexByIndexFrom(w, i);
                                        int offsetp = offset - hb.EdgeLengthKmers(e);
                                        vec< pair<int,int> > pp(p);
                                        pp.push( ep, offsetp );
                                        partials2.push_back(pp);    }    }    }
                         partials = partials2;    }    }

               vec<int> qsum_full_paths( full_paths.size( ), 0 );
               vec<int> qsum2_full_paths( full_paths.size( ), 0 );
               for ( int i = 0; i < full_paths.isize( ); i++ )
               {    const vec< pair<int,int> >& p = full_paths[i];
                    int pos = 0;
                    int estart = 0;
                    for ( int l = 0; l < p.isize( ); l++ )
                    {    if ( pos == bases[id].isize( ) ) break;
                         int e = p[l].first, offset = p[l].second;
                         for ( int j = estart; j < hb.EdgeLengthBases(e) - offset; 
                              j++ )
                         {    if ( bases[id][pos] != hb.EdgeObject(e)[offset+j] )
                              {    qsum_full_paths[i] += qv[pos];
                                   if ( qv[pos] > 2 )
                                        qsum2_full_paths[i] += qv[pos];    }
                              pos++;
                              if ( pos == bases[id].isize( ) ) break;    }
                         estart = hb.K( ) - 1;    }    }

               // Test for good enough placement.

               if ( full_paths.nonempty( ) 
                    && Min(qsum_full_paths) <= max_to_call_good )
               {    placed = True;
                    continue;    }
               select.push_back(id);

               // Print report.

               cout << "\nread " << id << "\n";
               PRINT(m);
               cout << "starts:";
               for ( int j = 0; j < starts.isize( ); j++ )
               {    cout << " " << starts[j].first << "." << starts[j].second 
                         << "[qsum=" << qsum_starts[j] << "]";     }
               cout << "\n";
               for ( int i = 0; i < full_paths.isize( ); i++ )
               {    cout << "full path:";
                    for ( int j = 0; j < full_paths[i].isize( ); j++ )
                    {    cout << " " << full_paths[i][j].first << "."
                              << full_paths[i][j].second;    }
                    cout << " [qsum=" << qsum_full_paths[i] 
                         << ",qsum2=" << qsum2_full_paths[i] << "]\n";    }
               for ( int i = 0; i < right_termini.isize( ); i++ )
               {    cout << "right termini:";
                    for ( int j = 0; j < right_termini[i].isize( ); j++ )
                    {    cout << " " << right_termini[i][j].first << "."
                              << right_termini[i][j].second;    }
                    cout << "\n";    }
               cout << "left termini:";
               for ( int j = 0; j < left_termini.isize( ); j++ )
               {    cout << " " << left_termini[j].first << "." 
                         << left_termini[j].second;     }
               cout << "\n";    }

          if (placed) continue;
          if ( p.size( ) > 0 ) cout << "\n";
          bases[id].Print( 
               cout, "id=" + ToString(id) + ",qsum=" + ToString(qsum) );
          printed++;    }

     // Print summary.

     cout << "\n";
     cout << "selected = " << printed << " = " 
          << PERCENT_RATIO( 3, printed, (int) bases.size( ) ) << endl;
     cout << endl;

     // Find friends of selected reads.

     const int K = 24;
     vecbasevector bases2(bases), basesrc(bases);
     for ( int64_t id = 0; id < (int64_t) bases.size( ); id++ )
          basesrc[id].ReverseComplement( );
     bases2.Append(basesrc);
     vec< triple<kmer<K>,int,int> > kmers_plus2;
     vec< SerfVec<Friend> > friends( select.size( ) );
     MakeKmerLookup0( bases2, kmers_plus2 );
     for ( int64_t i = 0; i < kmers_plus2.jsize( ); i++ )
     {    int64_t j;
          for ( j = i + 1; j < kmers_plus2.jsize( ); j++ )
               if ( kmers_plus2[j].first != kmers_plus2[i].first ) break;

          if ( j - i == 2 )
          {    i = j - 1;
               continue;    }

          vec<int> ids;
          for ( int64_t k = i; k < j; k++ )
               ids.push_back( kmers_plus2[k].second );
          Sort(ids);
          int max_mult = 0;
          for ( int l = 0; l < ids.isize( ); l++ )
          {    int m = ids.NextDiff(l);
               max_mult = Max( max_mult, m - l );
               l = m - 1;    }

          if ( max_mult > 1 )
          {    i = j - 1;   
               continue;    }

          for ( int64_t k1 = i; k1 < j; k1++ )
          {    int id1 = kmers_plus2[k1].second, pos1 = kmers_plus2[k1].third;
               if ( id1 >= (int64_t) bases.size( ) ) continue;
               int id1x = BinPosition( select, id1 );
               if ( id1x < 0 ) continue;
               for ( int k2 = i; k2 < j; k2++ )
               {    if ( k2 == k1 ) continue;
                    int id2 = kmers_plus2[k2].second, pos2 = kmers_plus2[k2].third;
                    Bool rc2 = False;
                    if ( id2 >= (int) bases.size( ) )
                    {    rc2 = True;
                         id2 = id2 - (int) bases.size( );    }
                    friends[id1x].push_back( 
                         Friend( id2, pos1 - pos2, rc2 ) );    }    }
          i = j - 1;    }
     for ( int i = 0; i < select.isize( ); i++ )
     {    vec<Friend> f;
          for ( int j = 0; j < (int) friends[i].size( ); j++ )
               f.push_back( friends[i][j] );
          UniqueSort(f);
          friends[i].resize(0);
          for ( int j = 0; j < f.isize( ); j++ )
               friends[i].push_back( f[j] );    }

     for ( int i = 0; i < select.isize( ); i++ )
     {    int id1 = select[i];
          bvec const& b1 = bases[id1];
          qvec q1 = quals.begin()[id1];
          cout << "\nfriends of read " << id1 << endl;
          for ( int j = 0; j < (int) friends[i].size( ); j++ )
          {    const Friend& f = friends[i][j];
               int id2 = f.readId( );
               Bool rc2 = f.isRC( );
               int offset = f.offset( );
               bvec b2 = bases[id2];
               qvec q2 = quals.begin()[id2];
               if (rc2) 
               {    b2.ReverseComplement( );
                    q2.ReverseMe( );    }
               cout << "\n" << ( rc2 ? "rc" : "fw" ) << id2 << "." << offset << "\n";
               int pos1, pos2;
               if ( offset >= 0 )
               {    pos1 = offset;
                    pos2 = 0;    }
               else
               {    pos1 = 0;
                    pos2 = -offset;    }
               avector<int> gaps(1), lengths(1);
               gaps(0) = 0;
               lengths(0) 
                    = IntervalOverlap( 0, b1.isize( ), offset, offset + b2.isize( ) );
               align a( pos1, pos2, gaps, lengths );
               PrintVisualAlignment( True, cout, b1, b2, a, q1, q2 );    }    }

     int64_t total_friends = 0;
     for ( int i = 0; i < select.isize( ); i++ )
          total_friends += friends[i].size( );
     cout << double(total_friends) / double( select.size( ) )
          << " friends per selected read" << endl;
     for ( int i = 0; i < select.isize( ); i++ )
     {    cout << "\nfriends of read " << select[i] << " =";
          for ( int j = 0; j < (int) friends[i].size( ); j++ )
          {    const Friend& f = friends[i][j];
               cout << " " << ( f.isRC( ) ? "rc" : "fw" ) << f.readId( ) << "."
                    << f.offset( );    }
          cout << "\n";    }

     // Correct selected reads.

     vecbasevector xbases;
     xbases.reserve( select.size( ) );
     for ( int i = 0; i < select.isize( ); i++ )
          xbases.push_back( bases[ select[i] ] );
     const int SEP = 0;
     const int STDEV = 100;
     const String LIB = "woof";
     const size_t nreads = select.size( );
     PairsManager pairs( bases.size( ) );
     pairs.addLibrary( SEP, STDEV, LIB );
     size_t npairs = bases.size( ) / 2;
     for ( size_t pi = 0; pi < npairs; pi++ )
          pairs.addPairToLib( 2 * pi, 2 * pi + 1, 0 );
     pairs.makeCache( );
     vec<Bool> to_edit( xbases.size( ), True );
     vec<int> trim_to, trace_ids;
     long_logging logc( "", "" );
     ref_data ref;
     vec<ref_loc> readlocs;
     long_logging_control log_control( ref, &readlocs, "", "" );
     long_heuristics heur( "" );
     String tmp_dir = "tmp.xxx";
     vec<int> KS = {24,40};
     // vec<int> KS = {24};
     vecbasevector new_stuff;
     for ( int i = 0; i < KS.isize( ); i++ )
     {    int K = KS[i];
          trim_to.resize( bases.size( ) );
          vecbasevector bases_loc( select.size( ) );
          vecqualvector quals_loc( select.size( ) );
          const int batch = 100;
          #pragma omp parallel for schedule(dynamic, 1)
          for ( int64_t id1a = 0; id1a < (int64_t) select.size( ); id1a += batch )
          {    readstack stack;
               vec<Bool> suspect;
               for ( int64_t id1x = id1a; 
                    id1x < Min( id1a + batch, (int64_t) select.size( ) ); id1x++ )
               {    int64_t id1 = select[id1x];
                    BaseVec& currentBaseVec = bases_loc[id1x];
                    UCharVec& currentCharVec = quals_loc[id1x];
                    currentBaseVec = bases[id1];
                    currentCharVec = quals.begin()[id1];
                    trim_to[id1] = bases[id1].size( );
                    if ( bases[id1].empty( ) ) continue;
                    Friends aligns = friends[id1x];
                    if ( aligns.size() > static_cast<unsigned>(heur.MAX_STACK) ) 
                         continue;
                    stack.Initialize( id1, aligns, 0, aligns.size(),
                         readstack::strict, bases, quals, pairs );
                    stack.HighQualDiff( 30, 1, suspect );
                    stack.Erase(suspect);
                    if ( heur.HQ_DIFF_WINDOW )
                    {    stack.HighQualDiffWindow( suspect );
                         stack.Erase(suspect);    }
                    Bool verbose = False;
                    stack.CorrectAll( currentBaseVec, currentCharVec, trim_to[id1], 
                         verbose );    

                    // if ( BinMember( trace_ids, id1 ) )
                    /*
                    {    vec<basevector> cons;
                         cons.push_back( currentBaseVec );
                         #pragma omp critical
                         {    cout << "\nCorrect1Pre, K = " << K
                                   << ", final stack, tracing read " << id1 << endl;
                              cout << "stack:\n";
                              stack.Print( cout, cons );    }    }
                    */

                         }    }

          /*
          #pragma omp parallel for schedule(dynamic)
          for ( size_t id1x = 0 ; id1x < select.size() ; ++id1x )
          {    const int64_t& id1 = select[id1x];
               bases[id1] = bases_loc[id1x];
               quals[id1] = quals_loc[id1x];    }    }
          */

          if ( i == KS.isize( ) - 1 )
          {    cout << "\ncorrected reads:\n";
               for ( int i = 0; i < select.isize( ); i++ )
               {    int id = select[i];
                    basevector b = bases_loc[i];
                    // b.resize( trim_to[id] );
                    new_stuff.push_back(b);
                    b.Print( cout, "id=" + ToString(id) );    }    }

          }


     /*
     cout << "\ncorrected reads:\n";
     for ( int i = 0; i < select.isize( ); i++ )
     {    int id = select[i];
          bases[id].Print( cout, "id=" + ToString(id) );    }    
     */


     const int K0 = 60;
     unsigned const COVERAGE = 50u;
     HyperBasevector hbn;
     LongReadsToPaths( new_stuff, K0, COVERAGE, False, False, &hbn );

     cout << "hbn has " << hbn.EdgeObjectCount( ) << " edges" << endl;

     Ofstream( out, work_dir + "/special/xxx.dot" );
     hbn.PrintSummaryDOT0w( out, True, False, True );
     Ofstream( eout, work_dir + "/special/xxx.fasta" );
     for ( int e = 0; e < hbn.EdgeObjectCount( ); e++ )
          hbn.EdgeObject(e).Print( eout, e );
     BinaryWriter::writeFile( work_dir + "/special/xxx.hbv", hbn );    }

void FixPath( ReadPath& p, const vec< triple<int,int,int> >& merges,
     const vec< vec<int> >& merges_index, const HyperBasevector& hb_orig )
{
     // mark path edges no longer in the hbv as id=-1
     for ( int l = 0; l < (int) p.size( ); l++ )
          if ( p[l] >= hb_orig.EdgeObjectCount( ) ) p[l] = -1;
     while(1)
     {    Bool changed = False;
          // for each valid path edge (not -1), store and sort a
          // list of merge ids involving the edge
          vec<int> mids;
          for ( int l = 0; l < (int) p.size( ); l++ )
               if ( p[l] >= 0 ) mids.append( merges_index[ p[l] ] );
          UniqueSort(mids);
          // for each merge id, iteratively update the path to the
          // new (merged) edges
          for ( int mj = 0; mj < mids.isize( ); mj++ )
          {    int j = mids[mj];
               int e1 = merges[j].first, e2 = merges[j].second;
               ForceAssert( e1 != e2 );
               int enew = merges[j].third;
               for ( int l = 0; l < (int) p.size( ); l++ )
               {
                    // case 1: adjacent path edges cover a merged pair
                    //           p[l] == e1, p[l+1] == e2 -> p[l]=p[l+1]=enew
                    if ( l < ( (int) p.size( ) ) - 1 && p[l] == e1 && p[l+1] == e2 )
                    {    p[l] = enew;
                         for ( int m = l+2; m < (int) p.size( ); m++ )
                              p[m-1] = p[m];
                         p.pop_back( );
                         changed = True;
                    }
                    // case 2: last path edge was start of a merged pair
                    else if ( l == ( (int) p.size( ) ) - 1 && p[l] == e1 )
                    {
/*                         size_t lastSkip = p.getLastSkip();
                         p.setLastSkip(lastSkip + hb_orig.EdgeLengthKmers(e2)); */
                         p[l] = enew;
                         changed = True;
                    }
                    // case 3: first path edge was end of a merged pair
                    else if ( l == 0 && p[l] == e2 )
                    {    int offset = p.getOffset();
                         p.setOffset( offset + hb_orig.EdgeLengthKmers(e1));
                         p[l] = enew;
                         changed = True;    }
                    if (changed) break;    }
               if (changed) break;    }

          if ( !changed ) break;   // dump out if we made a clean pass

     }    // while(1)
}

void RemoveUnneededVertices2( HyperBasevector& hbv, vec<int>& inv, 
     ReadPathVec& paths, Bool debug, const Bool single )
{
    static int debug_serial = 0;
    ++debug_serial;
    std::ostringstream debug_fnam_head;
    debug_fnam_head << "graph." << std::setprecision(3) << debug_serial;

    if ( debug )
            BinaryWriter::writeFile( debug_fnam_head.str() + ".BEFORE.hbv", hbv );

    // new algorithm
    // 1. make a list of vertices to kill
    // 2. find the "boundary" edges at the front and tail of the run to remove
    // 3. walk edges between/including boundary edges, building up new edge object,
    //    and recording edge mappings.
    // 4. add new edge object between boundary vertices;
    // 5. delete all edges marked for deletion
    //
    // this all assumes that the involution can't share edges in any
    // string of edges that we're interested in.  The only way I can see
    // this happening is a single, palindromic edge, but that would not
    // qualify for

    vec<bool> vertex_kill( hbv.N(), false );
    vec<int> vertex_queue;
    vec<int> to_left, to_right;
    hbv.ToLeft(to_left);
    hbv.ToRight(to_right);

    // step 1: make a list of vertices to kill
    // o----o----o----o
    // v0   v1   v2   v3
    for ( int v = 0; v < hbv.N(); ++v )
        if ( hbv.FromSize(v) == 1 && hbv.ToSize(v) == 1
                && hbv.From(v)[0] != hbv.To(v)[0] 
                && hbv.Bases( hbv.IFrom( v, 0 ) ) > 0
                && hbv.Bases( hbv.ITo( v, 0 ) ) > 0 ) {
            vertex_kill[v] = true;
            vertex_queue.push_back(v);
        }
    if ( debug ) {
        PRINT( vertex_queue.size() );
        cout << "serialno: " << debug_serial << endl;
        cout << "vertices (edge left,right): " << endl;
        for ( auto v : vertex_queue )
            cout << v << " (" <<
            hbv.EdgeObjectIndexByIndexTo(v,0) << ","
            <<  hbv.EdgeObjectIndexByIndexFrom(v,0) << ")" << endl;
        cout << "---" << endl;
    }

    // step 2
    vec<pair<int,int>> bound;
    while ( vertex_queue.size() ) {
        int v = vertex_queue.back();
        vertex_queue.pop_back();
        if ( !vertex_kill[v] ) continue;     // already done it
        int eleft;
        int vleft = v;
        size_t runsize = 0;
        do {
            if ( debug ) cout << "vleft = " << vleft << endl;
            runsize++;
            vertex_kill[vleft] = false;
            eleft = hbv.EdgeObjectIndexByIndexTo(vleft,0);
            vleft = hbv.To(vleft)[0];
        } while ( vertex_kill[vleft] );
        int eright;
        int vright = v;
        do {
            if ( debug ) cout << "vright = " << vright << endl;
            runsize++;
            vertex_kill[vright] = false;
            eright = hbv.EdgeObjectIndexByIndexFrom(vright,0);
            vright = hbv.From(vright)[0];
        } while ( vertex_kill[vright] );
        runsize--;      // 'v' gets counted twice

        // We rely on the fact that the involution is not
        // tied up with the path here.  We decide to push on the involution
        // here, too, so that we *know* what the inv[] of the new edge is.  However,
        // this requires that we canonicalize, so we don't do this twice.  This
        // canonicalization looks odd, but is correct (I think).
        if ( eleft < inv[eright] ) {
            // WARNING: code below relies on the fact that we're pushing on a
            // run and its involution adjacent in this list.
            bound.push(eleft,eright);
            bound.push(inv[eright], inv[eleft]);

            if ( debug ) {
                cout << "eleft = " << eleft << ", eright = " << eright << endl;
                cout << "inv eleft = " << inv[eleft] <<
                        ", inv eright = " << inv[eright] << endl;
                cout << "runsize = " << runsize << endl;
                cout << "===" << endl;
            }
        }
    }
    if ( debug ) PRINT(bound.size());

    // validate step 2
    if ( debug ) {
        for ( auto const& b : bound ) {
            int vleft = to_right[b.first];
            int vright = to_left[b.second];

            for ( int v = vleft; v != vright;  ) {
                ForceAssertEq( hbv.FromSize(v), 1 );
                v = hbv.From(v)[0];
                ForceAssertEq( hbv.ToSize(v), 1 );
            }
        }
    }

    // steps 3 and 4
    vec<int> edge_renumber0( hbv.EdgeObjectCount(), vec<int>::IDENTITY );
    vec<int> offsets(hbv.EdgeObjectCount(),0);
    vec<int> new_edge_numbers;
    vec<int> to_delete;
    while ( bound.size() ) {
        auto bounds = bound.back();
        bound.pop_back();

        int new_edge_no = hbv.EdgeObjectCount();

        basevector new_edge( hbv.EdgeObject( bounds.first ) );
        int off = hbv.EdgeLengthKmers( bounds.first );
        edge_renumber0[ bounds.first ] = new_edge_no;
        to_delete.push_back( bounds.first );

        for ( int v = to_right[bounds.first]; v != to_right[bounds.second]; v = hbv.From(v)[0] ) {
            // for each edge...
            int edge = hbv.EdgeObjectIndexByIndexFrom(v,0);
            to_delete.push_back(edge);
            new_edge = TrimCat(hbv.K(), new_edge, hbv.EdgeObject(edge));
            offsets[edge] = off;
            edge_renumber0[edge] = new_edge_no;
            off += hbv.EdgeLengthKmers(edge);
        }
        hbv.AddEdge(to_left[bounds.first], to_right[bounds.second], new_edge);
        if ( debug ) {
            cout << "run from edge " << bounds.first << " to " << bounds.second <<
                    " replaced by edge " << new_edge_no << endl;
        }
        new_edge_numbers.push_back(new_edge_no);
    }
    hbv.DeleteEdges(to_delete);

    if ( debug )
        BinaryWriter::writeFile( debug_fnam_head.str() + ".AFTER.hbv", hbv );


    // for each pair of newly created edges, update mInv
    inv.resize(hbv.EdgeObjectCount() );
    for (auto itr = new_edge_numbers.begin();
            itr != new_edge_numbers.end(); advance(itr,2)) {
         inv[itr[0]] = itr[1];
         inv[itr[1]] = itr[0];
    }

    // update the read paths for the newly created edges

     int nthreads = ( single ? 1 : omp_get_max_threads( ) );
     #pragma omp parallel for num_threads(nthreads)

    for ( size_t i = 0; i < paths.size(); ++i ) {
        auto& path = paths[i];
        if ( path.size() ) {
            SerfVec<int> oldPath = path;
            path.clear();
            path.setOffset( path.getOffset() + offsets[oldPath[0]]);
            path.push_back( edge_renumber0[ oldPath[0] ] );
            for ( auto itr = oldPath.begin()+1; itr != oldPath.end(); ++itr )
                if ( edge_renumber0[ *itr ] != path.back() )
                    path.push_back( edge_renumber0[ *itr ] );
        }
    }
    Validate( hbv, paths );
    TestInvolution(hbv, inv);
}

void RemoveUnneededVertices( HyperBasevector& hb, vec<int>& inv,
     ReadPathVec& paths )
{    double clock1 = WallClockTime( );
     vec< triple<int,int,int> > merges;
     vec<int> to_left, to_right;
     hb.ToLeft(to_left), hb.ToRight(to_right);
     for ( int i = 0; i < hb.N( ); i++ )
     {    if ( hb.From(i).size( ) == 1 && hb.To(i).size( ) == 1 
               && hb.From(i)[0] != i )
          {    int e1 = hb.EdgeObjectIndexByIndexTo( i, 0 );
               int e2 = hb.EdgeObjectIndexByIndexFrom( i, 0 );
               basevector p = hb.Cat( e1, e2 );
               int enew = hb.EdgeObjectCount( );
               int re1 = inv[e1], re2 = inv[e2];
               if ( re1 < 0 || re2 < 0 ) continue; // ******************************
               int v = hb.To(i)[0], w = hb.From(i)[0];
               Bool loop = ( v == w && hb.From(v).solo( ) && hb.To(v).solo( ) );
               // v --e1--> i --e2--> w
               merges.push( e1, e2, enew );
               hb.JoinEdges( i, p );
               to_left.push_back(v), to_right.push_back(w);
               // Case 1.
               if ( re1 == e2 && re2 == e1 )
               {    // * --e1=re2--> * --e2=re1--> *
                    inv.push_back(enew);    }
               // Case 2a.
               else if ( re2 == e2 && re1 != e1 )
               {    // * --e1--> * --e2=re2--> * --re1--> *
                    int enew2 = hb.EdgeObjectCount( );
                    merges.push( enew, re1, enew2 );
                    basevector p2 = TrimCat( hb.K( ), p, hb.EdgeObject(re1) );
                    to_left.push_back(v), to_right.push_back( to_right[re1] );
                    hb.JoinEdges( w, p2 );
                    cout << "RUN case 2a" << endl; // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                    inv.push_back( -1, enew2 );    }
               // Case 2b.
               else if ( re1 == e1 && re2 != e2 )
               {    // * --re2--> * --e1=re1--> * --e2--> *
                    int enew2 = hb.EdgeObjectCount( );
                    merges.push( re2, enew, enew2 );
                    basevector p2 = TrimCat( hb.K( ), hb.EdgeObject(re2), p );
                    to_left.push_back( to_left[re2] ), to_right.push_back(w);
                    hb.JoinEdges( v, p2 );
                    cout << "RUN case 2b" << endl; // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                    inv.push_back( -1, enew2 );    }
               // Case 3.
               else if ( re1 == e1 && re2 == e2 )
               {    cout << "RUN case 3" << endl; // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                    if (loop) inv.push_back(-1);
                    else
                    {    // not sure if this can happen
                         ForceAssert( 0 == 1 );    }    }
               // Case 4.
               else
               {    // e1, e2, re1, re2 all different
                    int renew = hb.EdgeObjectCount( );
                    basevector rp = hb.Cat( re2, re1 );
                    merges.push( re2, re1, renew );
                    int ri = to_right[re2];
                    hb.JoinEdges( ri, rp );
                    int rv = to_left[re2], rw = to_right[re1];
                    to_left.push_back(rv), to_right.push_back(rw);
                    inv.push_back(renew, enew);    }    }    }
     vec< vec<int> > merges_index( hb.EdgeObjectCount( ) );
     for ( int i = 0; i < merges.isize( ); i++ )
     {    merges_index[ merges[i].first ].push_back(i);
          merges_index[ merges[i].second ].push_back(i);    }
     LogTime( clock1, "removing unneeded vertices 1" );
     double clock2 = WallClockTime( );
     #pragma omp parallel for
     for ( int64_t i = 0; i < (int64_t) paths.size( ); i++ )
          FixPath( paths[i], merges, merges_index, hb );
     LogTime( clock2, "removing unneeded vertices 2" );
     double clock3 = WallClockTime( );
     hb.RemoveEdgelessVertices( );
     LogTime( clock3, "removing unneeded vertices 3" );    }

void RemoveUnneededVerticesLoopsOnly( HyperBasevector& hb, vec<int>& inv,
     ReadPathVec& paths )
{    double clock1 = WallClockTime( );
     vec< triple<int,int,int> > merges;
     vec<int> to_left, to_right;
     hb.ToLeft(to_left), hb.ToRight(to_right);
     for ( int i = 0; i < hb.N( ); i++ )
     {    if ( hb.From(i).size( ) == 1 && hb.To(i).size( ) == 1 
               && hb.From(i)[0] != i )
          {    int e1 = hb.EdgeObjectIndexByIndexTo( i, 0 );
               int e2 = hb.EdgeObjectIndexByIndexFrom( i, 0 );
               if ( hb.EdgeLengthBases(e1) == 0 || hb.EdgeLengthBases(e2) == 0 )
                    continue;
               int v = hb.To(i)[0], w = hb.From(i)[0];
               Bool loop = ( v == w && hb.From(v).solo( ) && hb.To(v).solo( ) );
               if ( !loop ) continue;
               basevector p = hb.Cat( e1, e2 );
               int enew = hb.EdgeObjectCount( );
               int re1 = inv[e1], re2 = inv[e2];

               // v --e1--> i --e2--> w

               // Test for nasty condition.  The archetypal case here is 
               // e1 = (AT)^100, e2 = (TA)^100.  It is not possible to merge the 
               // edges and maintain consistent data structures.

               if ( re1 == e1 && re2 == e2 ) continue;

               // Start the merge.

               merges.push( e1, e2, enew );
               hb.JoinEdges( i, p );
               to_left.push_back(v), to_right.push_back(w);

               // Case 1.

               if ( re1 == e2 && re2 == e1 )
               {    // * --e1=re2--> * --e2=re1--> *

                    inv.push_back(enew);    }

               // Case 3.
               else
               {    // e1, e2, re1, re2 all different
                    int renew = hb.EdgeObjectCount( );
                    basevector rp = hb.Cat( re2, re1 );
                    merges.push( re2, re1, renew );
                    int ri = to_right[re2];
                    hb.JoinEdges( ri, rp );
                    int rv = to_left[re2], rw = to_right[re1];
                    to_left.push_back(rv), to_right.push_back(rw);
                    inv.push_back(renew, enew);    }    }    }
     vec< vec<int> > merges_index( hb.EdgeObjectCount( ) );
     for ( int i = 0; i < merges.isize( ); i++ )
     {    merges_index[ merges[i].first ].push_back(i);
          merges_index[ merges[i].second ].push_back(i);    }
     LogTime( clock1, "removing unneeded vertices 1" );
     double clock2 = WallClockTime( );
     #pragma omp parallel for
     for ( int64_t i = 0; i < (int64_t) paths.size( ); i++ )
          FixPath( paths[i], merges, merges_index, hb );
     LogTime( clock2, "removing unneeded vertices 2" );
     double clock3 = WallClockTime( );
     hb.RemoveEdgelessVertices( );
     LogTime( clock3, "removing unneeded vertices 3" );    }

void RemoveUnneededVerticesGeneralizedLoops( HyperBasevector& hb, vec<int>& inv,
     ReadPathVec& paths )
{    double clock1 = WallClockTime( );
     vec<int> to_left, to_right;
     hb.ToLeft(to_left), hb.ToRight(to_right);
     VecULongVec paths_index;
     invert( paths, paths_index, hb.EdgeObjectCount( ) );
     vec<Bool> processed( hb.N( ), False );
     vec<int> dels;
     for ( int i = 0; i < hb.N( ); i++ )
     {    if ( !processed[i] && hb.From(i).size( ) == 1 && hb.To(i).size( ) == 1 
               && hb.From(i)[0] != i )
          {    vec<int> chain;
               int v = i;
               Bool fail = False;
               while(1)
               {    chain.push_back(v);
                    v = hb.From(v)[0];
                    if ( hb.From(v).size( ) != 1 || hb.To(v).size( ) != 1 
                         || hb.From(v)[0] == v )
                    {    fail = True;
                         break;    }
                    if ( Member( chain, v ) ) break;    }
               if (fail) continue;
               vec<int> echain, rechain;
               for ( int j = 0; j < chain.isize( ); j++ )
                    echain.push_back( hb.IFrom( chain[j], 0 ) );
               for ( int j = echain.isize( ) - 1; j >= 0; j-- )
                    rechain.push_back( inv[ echain[j] ] );
               if ( Meet2( echain, rechain ) ) continue;
               basevector x = hb.Cat(echain), rx = hb.Cat(rechain);
               dels.append(echain), dels.append(rechain);
               for ( int j = 0; j < chain.isize( ); j++ )
               {    processed[to_left[chain[j]]] = True;
                    processed[to_right[chain[j]]] = True;
                    processed[to_left[rechain[j]]] = True;
                    processed[to_right[rechain[j]]] = True;    }
               hb.AddEdge( to_left[echain[0]], to_left[echain[0]], x );
               hb.AddEdge( to_left[rechain[0]], to_left[rechain[0]], rx );
               inv.push_back( hb.E( ) - 1, hb.E( ) - 2 );
               for ( int pass = 1; pass <= 2; pass++ )
               {    const vec<int>& c = ( pass == 1 ? echain : rechain );
                    for ( int j = 0; j < c.isize( ); j++ )
                    {    int e = c[j];
                         for ( int u = 0; u < (int) paths_index[e].size( ); u++ )
                         {    ReadPath& p = paths[ paths_index[e][u] ];
                              if ( p[0] != e ) continue;
                              int offset = p.getOffset( );
                              for ( int l = 0; l < j; l++ )
                                   offset += hb.Kmers( c[l] );
                              p.setOffset(offset);
                              p.resize(1);
                              if ( pass == 1 ) p[0] = hb.E( ) - 2;    
                              else p[0] = hb.E( ) - 1;    }    }    }    }    }
     hb.DeleteEdges(dels);
     CleanupCore( hb, inv, paths );    }

void RemoveSmallComponents3( HyperBasevector& hb, const Bool remove_small_cycles )
{    double clock1 = WallClockTime( );
     const int max_small_comp = 1000;
     const int min_circle = 200;
     vec<int> e_to_delete;
     vec< vec<int> > comps;
     hb.Components(comps);
     LogTime( clock1, "removing small components 1" );
     double clock2 = WallClockTime( );
     #pragma omp parallel for
     for ( size_t i = 0; i < comps.size( ); i++ )
     {    const vec<int>& o = comps[i];

          int max_edge = 0;
          for ( int j = 0; j < o.isize( ); j++ )
          for ( int l = 0; l < hb.From( o[j] ).isize( ); l++ )
          {    int e = hb.EdgeObjectIndexByIndexFrom( o[j], l );
               max_edge = max( max_edge, hb.EdgeLengthKmers(e) );    }
          if ( max_edge > max_small_comp ) continue;

          // Remove small cycles.

          int total_kmers = 0;
          for ( int j = 0; j < o.isize( ); j++ )
          for ( int l = 0; l < hb.From( o[j] ).isize( ); l++ )
          {    int e = hb.EdgeObjectIndexByIndexFrom( o[j], l );
               total_kmers += hb.EdgeLengthKmers(e);    }
          if ( total_kmers < min_circle && remove_small_cycles )
          {    for ( size_t j = 0; j < o.size( ); j++ )
               {    int v = o[j];
                    #pragma omp critical
                    {    for ( size_t t = 0; t < hb.From(v).size( ); t++ )
                         {    e_to_delete.push_back( hb.EdgeObjectIndexByIndexFrom( 
                                   v, t ) );    }    }    }
               continue;    }
          if ( hb.HasCycle(o) ) continue;

          int no = o.size( );
          digraphE<basevector> GX( digraphE<basevector>::COMPLETE_SUBGRAPH, hb, o );
          vec<int> L, sources, sinks, p;
          for ( int i = 0; i < GX.EdgeObjectCount( ); i++ )
               L.push_back( GX.EdgeObject(i).size( ) - hb.K( ) + 1 );
          digraphE<int> G( GX, L );
          G.Sources(sources), G.Sinks(sinks);
          G.AddVertices(2);
          for ( int j = 0; j < sources.isize( ); j++ )
               G.AddEdge( no, sources[j], 0 );
          for ( int j = 0; j < sinks.isize( ); j++ )
               G.AddEdge( sinks[j], no+1, 0 );
          for ( int e = 0; e < G.EdgeObjectCount( ); e++ )
               G.EdgeObjectMutable(e) = -G.EdgeObject(e);
          G.ShortestPath( no, no+1, p );
          int max_path = 0;
          for ( int j = 0; j < p.isize( ) - 1; j++ )
          {    int v1 = p[j], v2 = p[j+1];
               int m = 1000000000;
               for ( int l = 0; l < G.From(v1).isize( ); l++ )
               {    if ( G.From(v1)[l] == v2 )
                         m = Min( m, G.EdgeObjectByIndexFrom( v1, l ) );    }
               max_path -= m;    }

          if ( max_path <= max_small_comp )
          {    for ( size_t j = 0; j < o.size( ); j++ )
               {    int v = o[j];
                    #pragma omp critical
                    {    for ( size_t t = 0; t < hb.From(v).size( ); t++ )
                         {    e_to_delete.push_back( hb.EdgeObjectIndexByIndexFrom( 
                                   v, t ) );    }    }    }    }    }
     LogTime( clock2, "removing small components 2" );
     double clock3 = WallClockTime( );
     hb.DeleteEdges(e_to_delete);
     LogTime( clock3, "removing small components 3" );    }

void Empty( const HyperBasevector& hb, const vec< pair<vec<int>,vec<int>> >& pairs, 
     const vec<int64_t>& pairs_pid, vec<vec<int>>& left_empty, 
     vec<vec<int>>& right_empty, const Bool EMPTY2 )
{
     int nedges = hb.EdgeObjectCount( );
     left_empty.resize(nedges), right_empty.resize(nedges);
     if (EMPTY2)
     {    for ( int l = 0; l < pairs.isize( ); l++ )
          {    if ( pairs[l].first.nonempty( ) && pairs[l].second.empty( ) )
               {    left_empty[ pairs[l].first.back( ) ].push_back( pairs_pid[l] );
                    if ( pairs[l].first.front( ) != pairs[l].first.back( ) )
                    {    left_empty[ pairs[l].first.front( ) ].
                              push_back( pairs_pid[l] );    }    }
               if ( pairs[l].second.nonempty( ) && pairs[l].first.empty( ) )
               {    right_empty[ pairs[l].second.front( ) ]
                         .push_back( pairs_pid[l] );    
                    if ( pairs[l].second.front( ) != pairs[l].second.back( ) )
                    {    right_empty[ pairs[l].second.back( ) ]
                              .push_back( pairs_pid[l] );    }    }    }    }
     else
     {    for ( int l = 0; l < pairs.isize( ); l++ )
          {    if ( pairs[l].first.nonempty( ) && pairs[l].second.empty( ) )
                    left_empty[ pairs[l].first.back( ) ].push_back( pairs_pid[l] );
               if ( pairs[l].second.nonempty( ) && pairs[l].first.empty( ) )
               {    right_empty[ pairs[l].second.front( ) ]
                         .push_back( pairs_pid[l] );    }    }    }    }

void Validate( const HyperBasevector& hb, const ReadPathVec& paths ) {
    if (ValidateAllReadPaths(hb, paths) == false)
	TracebackThisProcess();
}

void Validate( const HyperBasevector& hb, const HyperBasevectorX& hbx, const ReadPathVecX& paths ) {
    if (ValidateAllReadPaths(hb, hbx, paths) == false)
	TracebackThisProcess();
}
void Validate( const HyperBasevector& hb, const ReadPathVecX& paths ) {
    if (ValidateAllReadPaths(hb, paths) == false)
	TracebackThisProcess();
}

void TestIndex( const HyperBasevector& hb,
        const ReadPathVec& paths, const VecULongVec& invPaths)
{
    // horribly inefficient, just meant to be quick to code
    bool good = true;

    // is each edge entry reflected in a path?
    for ( size_t edge = 0; edge < invPaths.size(); ++edge )
        for ( auto const readid : invPaths[edge] )
            if ( std::find( paths[readid].begin(),
                    paths[readid].end(), edge ) == paths[readid].end() ) {
                cout << "the index for edge " << edge << " names readid "
                        << readid << " but the readpath says "
                        << printSeq( paths[readid] ) << endl;
                good = false;
            }

    // is each path entry reflected in an edge index
    for ( size_t pathid = 0; pathid < paths.size(); ++pathid )
        for ( auto const edge : paths[pathid] )
            if ( std::find( invPaths[edge].begin(),
                    invPaths[edge].end(), pathid ) == invPaths[edge].end() ) {
                cout << "the path for read " << pathid << " names edge " <<
                edge << " but the index for that edge says " <<
                printSeq(invPaths[edge]) << endl;
                good = false;
            }

    if ( !good ) TracebackThisProcess();
}


void TestInvolution( const HyperBasevector& hb, const vec<int>& inv )
{
     vec<int> to_left, to_right;
     vec<Bool> used;
     hb.ToLeft(to_left), hb.ToRight(to_right);
     hb.Used(used);
     if ( inv.isize( ) != hb.EdgeObjectCount( ) )
     {    cout << "\n";
          PRINT2( hb.EdgeObjectCount( ), inv.size( ) );
          cout << "Involution has wrong size.\n" << "Abort." << endl;
          TracebackThisProcess( );    }
     for ( int e = 0; e < hb.EdgeObjectCount( ); e++ )
     {    if ( !used[e] ) continue;
          if ( inv[e] < 0 || inv[e] >= hb.EdgeObjectCount( ) )
          {    cout << "\n";
               PRINT3( e, inv[e], hb.EdgeObjectCount( ) );
               cout << "Illegal involution value.\n" << "Abort." << endl;
               TracebackThisProcess( );    }
          basevector b = hb.EdgeObject(e);
          b.ReverseComplement( );
          if ( b != hb.EdgeObject( inv[e] ) )
          {    cout << "\n";
               int re = inv[e];
               PRINT4( e, re, b.size( ), hb.EdgeObject(re).size( ) );
               cout << "Involution value not rc.\n" << "Abort." << endl;
               TracebackThisProcess( );    }
          if ( inv[inv[e]] != e )
          {    cout << "Involution is not an involution.\n" << "Abort." << endl;
               TracebackThisProcess( );    }    }
     for ( int v = 0; v < hb.N( ); v++ )
     {    for ( int i1 = 0; i1 < hb.To(v).isize( ); i1++ )
          for ( int i2 = 0; i2 < hb.From(v).isize( ); i2++ )
          {    int e1 = hb.EdgeObjectIndexByIndexTo( v, i1 );
               int e2 = hb.EdgeObjectIndexByIndexFrom( v, i2 );
               int re1 = inv[e1], re2 = inv[e2];
               if ( to_right[re2] != to_left[re1] )
               {    cout << "Involution does not preserve graph structure.\n";
                    PRINT3( e1, to_left[e1], to_right[e1] );
                    PRINT3( e2, to_left[e2], to_right[e2] );
                    PRINT3( re1, to_left[re1], to_right[re1] );
                    PRINT3( re2, to_left[re2], to_right[re2] );
                    cout << "Abort." << endl;
                    TracebackThisProcess( );    }    }    }
     // Note sure the following can happen:
     for ( int v = 0; v < hb.N( ); v++ )
     {    for ( int i = 0; i < hb.From(v).isize( ); i++ )
          {    int w = hb.From(v)[i];
               int e = hb.EdgeObjectIndexByIndexFrom( v, i );
               int re = inv[e];
               int rv = to_left[re], rw = to_right[re];
               if ( hb.From(rv).size( ) != hb.To(w).size( )
                    || hb.To(rw).size( ) != hb.From(v).size( ) )
               {    cout << "Graph structure is asymmetric.\n" << "Abort." << endl;
                    TracebackThisProcess( );    }    }    }     }

void FragDist( const HyperBasevector& hb, const vec<int>& inv,
     const ReadPathVec& paths, const String out_file )
{    
     // Generate fragment distribution.

     const int width = 10;
     const int max_sep = 1000;
     const int min_edge = 10000;
     vec<double> count( max_sep/width, 0 );
     for ( int64_t id1 = 0; id1 < (int64_t) paths.size( ); id1 += 2 )
     {    int64_t id2 = id1 + 1;
          if ( paths[id1].size( ) == 0 || paths[id2].size( ) == 0 ) continue;
          int e1 = paths[id1][0], e2 = inv[ paths[id2][0] ];
          int epos1 = paths[id1].getOffset( );
          if ( e1 != e2 ) continue;
          if ( hb.EdgeLengthBases(e1) < min_edge ) continue;
          int epos2 = hb.EdgeLengthBases(e2) - paths[id2].getOffset( );
          int len = epos2 - epos1;
          if ( len < 0 || len >= max_sep ) continue;
          count[ len/width ]++;    }

     // Output distribution.

     double total = Sum(count);
     {    Ofstream( out, out_file );
          out << "# fragment library size distribution" << endl;
          out << "# bins have diameter 10" << endl;
          out << "# line format:\n";
          out << "# bin_center mass" << endl;
          for ( int j = 0; j < count.isize( ); j++ )
               out << j * width + (width/2) << " " << count[j]/total << endl;    }

     // Check for abject failure.

     if ( total == 0 )
     {    Ofstream( out, out_file + ".png.FAIL" );
          Remove( out_file + ".png" );
          out << "Could not generate frags.dist.png because there was not\n"
               << "enough assembly to compute the distribution." << endl;
          return;    }

     // Check for missing executables.

     vec<String> missing;
     vec<String> ex = { "ps2epsi", "pstopnm", "pnmcrop", "pnmpad", "pnmtopng" };
     for ( auto executable : ex )
     {    if ( System( executable + " --version > /dev/null 2>&1" ) != 0 )
               missing.push_back(executable);    }
     if ( missing.nonempty( ) )
     {    Ofstream( out, out_file + ".png.FAIL" );
          out << "Could not generate frags.dist.png because the following "
               << "executables were not found:\n" << printSeq(missing) 
               << "." << endl;
          Remove( out_file + ".png" );
          return;    }

     // Make plot.
          
     String TITLE = "Fragment library size distribution";
     vec<graphics_primitive> points, lines;
     lines.push_back( SetLineWidth(1.0) );
     double xm = 0.0, ym = 0.0;
     for ( int j = 0; j < count.isize( ); j++ )
          points.push_back( Point( j * width + (width/2),  count[j]/total, 1.0 ) );
     points.push_back( SetColor(black) );
     points.append(lines);
     double x = 0, X = MaxX(points), y = 0, Y = MaxY(points);
     vec<graphics_primitive> x_axis = AxisX( x, X, 1.0, True, "", 0.0, False );
     x_axis.push_front( SetLineWidth(1.0) );
     vec<graphics_primitive> y_axis = AxisY( x, X, y, Y, True, "", 0.0, False );
     vec<graphics_primitive> points2(y_axis);
     points2.append(points);
     points = points2;
     vec< vec<graphics_primitive> > stack;
     vec<double> heights;
     stack.push_back(x_axis),    heights.push_back(0);
     stack.push_back(points),    heights.push_back(200);
     vec<graphics_primitive> title;
     title.push_back( TextCenter( TITLE, (X+x)/2.0, 0, 0, TimesBold(15) ) );
     stack.push_back(title);
     heights.push_back(20);
     String fail_msg = RenderGraphics( out_file + ".png", stack, heights, 1.0, 200, 
          50, True, 1.0, True );
     Remove( out_file + ".eps" );
     if ( fail_msg != "" )
     {    Remove( out_file + ".png" );
          Ofstream( out, out_file + ".png.FAIL" );
          out << "Could not generate frags.dist.png because something went "
               << "wrong, see below:\n\n" << fail_msg;
          return;    }    }

// UnwindThreeEdgePlasmids
//
// 1. Component has two vertices v, w.
//
// 2. Three edges: two from v to w (e1, e2), one from w to v (f)
//    (and if not, swap v and w).
//
// 3. e1, e2 well-supported (at least 10 pids each).
//
// 4. coverage of e1 and e2 within 25% of each other.
//
// Then replace by e1,f,e2,f loop.
//
// Note that this drops read placements that become ambiguous.
//
// Note that this will make a mistake if there are actually two plasmids present
// at nearly equal molarity.

void UnwindThreeEdgePlasmids(HyperBasevector& hb, vec<int>& inv, ReadPathVec& paths)
{
     // Heuristics.

     const int min_cov = 10;
     const double fudge = 0.25;
     const int min_links = 2;

     // Set up indices.

     VecULongVec paths_index;
     invert( paths, paths_index, hb.EdgeObjectCount( ) );
     vec<int> to_right;
     hb.ToRight(to_right);

     // Look for the special components.

     vec< vec<int> > comps;
     hb.Components(comps);
     vec<int> dels;
     cout << Date( ) << ": unwinding three-edge plasmids" << endl;
     for ( size_t i = 0; i < comps.size( ); i++ )
     {    const vec<int>& o = comps[i];
          if ( o.size( ) != 2 ) continue;
          int v = o[0], w = o[1];
          // PRINT2( v, w ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          if ( hb.From(v).size( ) != 2 ) swap( v, w );
          if ( hb.From(v).size( ) != 2 || hb.From(w).size( ) != 1 ) continue;
          if ( hb.From(v)[0] != w || hb.From(v)[1] != w || hb.From(w)[0] != v )
               continue;
          int e1 = hb.IFrom( v, 0 ), e2 = hb.IFrom( v, 1 ), f = hb.IFrom( w, 0 );
          // PRINT3( e1, e2, f ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          int re1 = inv[e1], re2 = inv[e2], rf = inv[f];
          vec<int> v1 = {e1,e2,f}, v2 = {re1,re2,rf};
          if ( Meet2( v1, v2 ) || Min(v2) < Min(v1) ) continue;

          // Make sure e1 and e2 are linked.

          vec<int64_t> pids;
          for ( int j = 0; j < (int) paths_index[e1].size( ); j++ )
               pids.push_back( paths_index[e1][j]/2 );
          for ( int j = 0; j < (int) paths_index[e2].size( ); j++ )
               pids.push_back( paths_index[e2][j]/2 );
          for ( int j = 0; j < (int) paths_index[f].size( ); j++ )
               pids.push_back( paths_index[f][j]/2 );
          for ( int j = 0; j < (int) paths_index[re1].size( ); j++ )
               pids.push_back( paths_index[re1][j]/2 );
          for ( int j = 0; j < (int) paths_index[re2].size( ); j++ )
               pids.push_back( paths_index[re2][j]/2 );
          for ( int j = 0; j < (int) paths_index[rf].size( ); j++ )
               pids.push_back( paths_index[rf][j]/2 );
          UniqueSort(pids);
          int links = 0;
          for ( int j = 0; j < pids.isize( ); j++ )
          {    int64_t id1 = 2*pids[j], id2 = 2*pids[j] + 1;
               vec<int> es;
               for ( int l = 0; l < (int) paths[id1].size( ); l++ )
                    es.push_back( paths[id1][l], inv[ paths[id1][l] ] );
               for ( int l = 0; l < (int) paths[id2].size( ); l++ )
                    es.push_back( paths[id2][l], inv[ paths[id2][l] ] );
               UniqueSort(es);
               if ( BinMember( es, e1 ) && BinMember( es, e2 ) ) links++;    }
          // PRINT(links); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          if ( links < min_links ) continue;

          // Make sure that e1 and e2 have enough support.  Otherwise we would 
          // worry that one of them is an error branch.

          int ne1 = Pids( e1, hb, inv, paths, paths_index ).size( );
          int ne2 = Pids( e2, hb, inv, paths, paths_index ).size( );
          if ( ne1 < min_cov || ne2 < min_cov ) continue;

          // Calculate a bunch of interesting stuff that we're not actually using.

          /*
          int nf = Pids( f, hb, inv, paths, paths_index ).size( ); // XXXXXXXXXXXXXX
          PRINT3( ne1, ne2, nf ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          PRINT3( hb.Kmers(e1), hb.Kmers(e2), hb.Kmers(f) ); // XXXXXXXXXXXXXXXXXXXX
          double ce1 = double(ne1) / ( hb.Kmers(e1) + hb.K( ) - 1 - 60 );
          double ce2 = double(ne2) / ( hb.Kmers(e2) + hb.K( ) - 1 - 60 );
          double r = Max(ce1,ce2)/Min(ce1,ce2);
          PRINT3( ce1, ce2, r ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          // double ce1 = RawCoverage( e1, hb, inv, paths, paths_index );
          // double ce2 = RawCoverage( e2, hb, inv, paths, paths_index );
          double cf = RawCoverage( f, hb, inv, paths, paths_index );
          if ( Max(ce1,ce2)/Min(ce1,ce2) - 1 > fudge ) continue;
          double r1 = (cf/ce1)/2.0, r2 = (cf/ce2)/2.0;
          PRINT2( r1, r2 ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          double low = 1 - fudge, high = 1 + fudge;
          */

          {    
               // Start assembly edit.

               // cout << "editing" << endl; // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
               vec<int> x = {e1,f,e2,f}, rx = {rf,re2,rf,re1};
               basevector b = hb.Cat(x), rb = hb.Cat(rx);
               dels.push_back( e1, e2, f, re1, re2, rf );
               int m = hb.AddEdge( v, v, b );
               int rv = to_right[re1];
               int rm = hb.AddEdge( rv, rv, rb );    
               inv.push_back( rm, m );

               // Get the relevant read pairs.

               vec<int64_t> pids;
               for ( int e : x ) 
                    pids.append( Pids( e, hb, inv, paths, paths_index ) );
               UniqueSort(pids);

               // Remap the reads.  It would be smarter to take into account
               // pairing.

               vec<int64_t> ids;
               for ( int j = 0; j < pids.isize( ); j++ )
                    ids.push_back( 2*pids[j], 2*pids[j] + 1 );
               for ( int j = 0; j < ids.isize( ); j++ )
               {    ReadPath& p = paths[ ids[j] ];
                    if ( p.size( ) == 0 ) continue;
                    Bool fixed = False;
                    for ( int l = 0; l < (int) p.size( ); l++ )
                    {    int x = p[l];
                         int pre = 0;
                         for ( int r = 0; r < l; r++ )
                              pre += hb.Kmers( p[r] );
                         fixed = True;
                         if ( x == e1 )
                         {    p[0] = m;
                              p.addOffset(-pre);    }
                         else if ( x == e2 )
                         {    p[0] = m;
                              p.addOffset( -pre + hb.Kmers(e1) + hb.Kmers(f) );    }
                         else if ( x == re2 )
                         {    p[0] = rm;
                              p.addOffset( -pre + hb.Kmers(f) );
                              fixed = True;    }
                         else if ( x == re1 )
                         {    p[0] = rm;
                              p.addOffset( -pre + 2*hb.Kmers(f) + hb.Kmers(e2) );   }
                         else fixed = False;
                         if (fixed) 
                         {    p.resize(1);
                              break;    }    }
                    if ( !fixed ) p.resize(0);    }    }    }

     // Finish assembly edit.

     hb.DeleteEdges(dels);
     CleanupCore( hb, inv, paths );    }
