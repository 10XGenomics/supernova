///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "MainTools.h"
#include "fastg/FastgGraph.h"
#include "paths/HyperBasevector.h"
#include "paths/long/CleanEfasta.h"
#include "paths/long/LargeKDispatcher.h"
#include "paths/long/MakeKmerStuff.h"
#include "paths/long/SupportedHyperBasevector.h"

namespace { // open anonymous namespace

template<int K2> void MapAssembly( const HyperBasevector& hb, 
     const HyperBasevector& hb_orig, vec< vec<int> >& content, const int verbosity )
{    if ( verbosity >= 1 ) cout << Date( ) << ": mapping hb onto hb_orig" << endl;
     vecbasevector B;
     for ( int e = 0; e < hb_orig.EdgeObjectCount( ); e++ )
          B.push_back_reserve( hb_orig.EdgeObject(e) );
     vec< triple<kmer<K2>,int,int> > kmers_plus;
     MakeKmerLookup1( B, kmers_plus );
     vec< kmer<K2> > kmers( kmers_plus.size( ) );
     for ( size_t i = 0; i < kmers.size( ); i++ )
          kmers[i] = kmers_plus[i].first;
     content.resize( hb.EdgeObjectCount( ) );
     for ( int e = 0; e < hb.EdgeObjectCount( ); e++ )
     {    const basevector& b = hb.EdgeObject(e);
          kmer<K2> x;
          int pos = 0;
          while( pos < b.isize( ) - (K2-1) )
          {    x.SetToSubOf( b, pos );
               int64_t low = LowerBound( kmers, x );
               int64_t high = UpperBound( kmers, x );

               // A kmer in hb cannot occur more than once in hb_orig.  However,
               // note that it is possible (but rare, and perhaps never occurred) 
               // to have a kmer in hb that is not present in hb_orig.  The reason 
               // for this is that at some point two edges in hb_orig might have 
               // been joined along a (K-1)-base overlap, resulting in the creation 
               // of new kmers.  It is not entirely clear what the implications of 
               // this are here.

               if ( high == low )
               {    pos++;
                    continue;    }
               ForceAssertEq( high, low + 1 );

               int id = kmers_plus[low].second;
               content[e].push_back(id);
               if ( verbosity >= 1 ) PRINT2( e, id );

               // There used to be an assert here, but it doesn't work when there
               // are instances of high = low, as discussed above.

               // ForceAssertEq( kmers_plus[low].third, 0 );

               pos += B[id].isize( ) - (K2-1) - kmers_plus[low].third;    }    }    }

template <int K>
struct MapAsmFunctor
{
    void operator()( const HyperBasevector& hb,  const HyperBasevector& hb_orig,
                        vec< vec<int> >& content, const int verbosity )
    { MapAssembly<K>(hb,hb_orig,content,verbosity); }
};

void MapOne( const vec<int>& u, const vec< vec< pair<int,int> > >& cindex,
     const vec< vec<int> >& content, const HyperBasevector& hb,
     const vec<int>& to_right, vec< vec<int> >& fulls, vec<int>& fulls_pos )
{
     fulls.clear( );
     fulls_pos.clear( );
     vec< triple< vec<int>, int, int > > partials;
     const vec< pair<int,int> >& locs = cindex[ u[0] ];
     for ( int j = 0; j < locs.isize( ); j++ )
     {    int x = locs[j].first, y = locs[j].second;
          Bool mismatch = False;
          int l;
          for ( l = 0; l < u.isize( ); l++ )
          {    if ( y+l == content[x].isize( ) ) break;
               if ( content[x][y+l] != u[l] ) 
               {    mismatch = True;
                    break;    }    }
          if (mismatch) continue;
          vec<int> p(1);
          p[0] = x;
          partials.push( p, l, y );    }
     while( partials.nonempty( ) )
     {    vec<int> p = partials.back( ).first;
          int pos = partials.back( ).second;
          int y = partials.back( ).third;
          partials.pop_back( );
          if ( pos == u.isize( ) )
          {    fulls.push_back(p);
               fulls_pos.push_back(y);
               continue;    }
          int v = to_right[ p.back( ) ];
          for ( int j = 0; j < hb.From(v).isize( ); j++ )
          {    int l, e = hb.EdgeObjectIndexByIndexFrom( v, j ); 
               Bool mismatch = False;
               for ( l = pos; l < u.isize( ); l++ )
               {    if ( l-pos == content[e].isize( ) ) break;
                    if ( content[e][l-pos] != u[l] ) 
                    {    mismatch = True;
                         break;    }    }
               if (mismatch) continue;
               vec<int> q(p);
               q.push_back(e);
               partials.push( q, l, y );    }    }    }

} // close anonymous namespace

// TransformPaths: replace each occurrence of x[i] by y[i].

void SupportedHyperBasevector::TransformPaths( const vec< vec<int> >& x, 
     const vec<int>& y, vec< vec< pair<int,int> > >& paths_index )
{    if ( paths_index.nonempty( ) ) paths_index.resize( EdgeObjectCount( ) );
     vec<int> all;
     for ( int j = 0; j < x.isize( ); j++ )
          all.append( x[j] );
     UniqueSort(all);

     vec<int> used_paths;
     if ( paths_index.nonempty( ) )
     {    for ( int i = 0; i < all.isize( ); i++ )
          {    int e = all[i];
               for ( int j = 0; j < paths_index[e].isize( ); j++ )
                    used_paths.push_back( paths_index[e][j].first );    }
          UniqueSort(used_paths);    }
     int to_use = ( paths_index.empty( ) ? NPaths( ) : used_paths.isize( ) );

     const int infty = 1000000000;
     #pragma omp parallel for
     for ( int pi = 0; pi < to_use; pi++ )
     {    int i = ( paths_index.empty( ) ? pi : used_paths[pi] );
          vec<int>& p = PathMutable(i);
          vec<int> p_old = p;
          for ( int j = 0; j < p.isize( ); j++ )
          {    if ( !BinMember( all, p[j] ) ) continue;

               // Find all the proper overlaps between p and some x[l], that use
               // the jth entry of p.

               vec< pair<int,int> > pos;
               for ( int l = 0; l < x.isize( ); l++ )
               {    for ( int m = 0; m < x[l].isize( ); m++ )
                    {    if ( x[l][m] != p[j] ) continue;
                         if ( Overlap( p, x[l], j - m ) ) pos.push(l,m);    }    }
               if ( !pos.solo( ) )
               {    p.clear( );
                    WeightFwMutable(i) = 0.0;
                    WeightRcMutable(i) = 0.0;
                    break;    }
               int l = pos[0].first, m = pos[0].second;
               int start = j;
               int stop = Min( p.isize( ), x[l].isize( ) + (j-m) );
               for ( int r = start; r < stop - 1; r++ )
                    p[r] = -1;
               p[stop-1] = y[l];    }
          RemoveNegatives(p);

          if ( paths_index.nonempty( ) && p != p_old )
          {    for ( int l = 0; l < p_old.isize( ); l++ )
               {    int e = p_old[l];
                    #pragma omp critical
                    {    auto low = LowerBound( paths_index[e], {i,0} );
                         auto high = UpperBound( paths_index[e], {i,infty} );
                         paths_index[e].erase( 
                              paths_index[e].begin( ) + low,
                              paths_index[e].begin( ) + high );    }    }
               #pragma omp critical
               {    for ( int l = 0; l < p.isize( ); l++ )
                    {    int e = p[l];
                         paths_index[e].push( i, l );
                         Sort( paths_index[e] );    }    }    }    }

     #pragma omp parallel for
     for ( int i = 0; i < NPairs( ); i++ )
     {    vec<int> &p1 = PairLeftMutable(i), &p2 = PairRightMutable(i);
          for ( int pass = 1; pass <= 2; pass++ )
          {    vec<int>& p = ( pass == 1 ? p1 : p2 );
               for ( int j = 0; j < p.isize( ); j++ )
               {    if ( !BinMember( all, p[j] ) ) continue;
                    vec< pair<int,int> > pos;
                    for ( int l = 0; l < x.isize( ); l++ )
                    {    for ( int m = 0; m < x[l].isize( ); m++ )
                         {    if ( x[l][m] != p[j] ) continue;
                              if ( Overlap( p, x[l], j - m ) ) 
                                   pos.push(l,m);    }    }
                    if ( !pos.solo( ) )
                    {    p1.clear( ), p2.clear( );
                         PairDataMutable(i).clear( );
                         break;    }
                    int l = pos[0].first, m = pos[0].second;
                    int start = j, stop = Min( p.isize( ), x[l].isize( ) + (j-m) );
                    int add = 0;
                    if ( pass == 1 && j - m + x[l].isize( ) > p.isize( ) )
                    {    for ( int r = p.isize( ); r < x[l].isize( ) + j-m; r++ )
                              add += EdgeLengthKmers( x[l][ r - (j-m) ] );    }
                    if ( pass == 2 && j - m < 0 )
                    {    for ( int r = 0; r < -(j-m); r++ )
                              add += EdgeLengthKmers( x[l][r] );    }
                    AddTrim( i, add );
                    for ( int r = start; r < stop - 1; r++ )
                         p[r] = -1;
                    p[stop-1] = y[l];    }
               RemoveNegatives(p);    }    }    }

void SupportedHyperBasevector::Bootstrap( const SupportedHyperBasevector& shb0,
     const long_logging& logc )
{    
     double clock = WallClockTime( );
     if (logc.STATUS_LOGGING) cout << Date( ) << ": begin bootstrap" << endl;
     const HyperBasevector& hb_orig = shb0;
     SupportedHyperBasevector& shb = *this;
     HyperBasevector& hb = shb;
     const int K2 = hb.K( );

     // Map reads.
     //
     // Map hb onto hb_orig.  For e an edge in hb, we set content[e] to
     // the sequence of edges in hb_orig that it maps to.

     vec< vec<int> > content;
     int verbosity = logc.verb[ "BOOTSTRAP" ];
     BigK::dispatch<MapAsmFunctor>(K2,hb,hb_orig,content,verbosity);

     // Index map.

     vec< vec< pair<int,int> > > cindex( hb_orig.EdgeObjectCount( ) );
     for ( int i = 0; i < content.isize( ); i++ )
     for ( int j = 0; j < content[i].isize( ); j++ )
          cindex[ content[i][j] ].push( i, j );

     // Map reads onto assembly.

     if (logc.STATUS_LOGGING)
          cout << Date( ) << ": mapping reads onto assembly" << endl;
     vec<int> to_left, to_right;
     hb.ToLeft(to_left), hb.ToRight(to_right);
     vec< vec<int> >& usux = shb.PathsMutable( );
     vec<fix64_6>& countx_fw = shb.WeightsFwMutable( );
     vec<fix64_6>& countx_rc = shb.WeightsRcMutable( );
     usux.clear( ), countx_fw.clear( ), countx_rc.clear( );
     vec< triple< vec<int>, fix64_6, fix64_6 > > all_fulls;
     for ( int i = 0; i < shb0.Paths( ).isize( ); i++ )
     {    vec< vec<int> > fulls;
          vec<int> fulls_pos;
          MapOne( shb0.Path(i), cindex, content, hb, to_right, fulls, fulls_pos );
          for ( int j = 0; j < fulls.isize( ); j++ )
               all_fulls.push( fulls[j], shb0.WeightFw(i) / fulls.size( ),
                    shb0.WeightRc(i) / fulls.size( ) );    }
     Sort(all_fulls);
     for ( int i = 0; i < all_fulls.isize( ); i++ )
     {    int j;
          for ( j = i + 1; j < all_fulls.isize( ); j++ )
               if ( all_fulls[j].first != all_fulls[i].first ) break;
          fix64_6 cfw = 0, crc = 0;
          for ( int k = i; k < j; k++ )
          {    cfw += all_fulls[k].second;
               crc += all_fulls[k].third;    }
          usux.push_back( all_fulls[i].first );
          countx_fw.push_back(cfw);
          countx_rc.push_back(crc);
          if ( verbosity >= 1 )
          {    cout << "[" << i+1 << "," << setiosflags(ios::fixed) 
                    << setprecision(1) << cfw << "+" << crc 
                    << resetiosflags(ios::fixed)
                    << "x] " << printSeq( usux.back( ) ) << "\n";    }
          i = j - 1;    }    

     vec< pair< vec<int>, vec<int> > >& usuy = shb.PairsMutable( );
     vec< vec<pair_point> >& county = shb.PairDataMutable( );
     usuy.clear( ), county.clear( );
     vec< pair< pair< vec<int>, vec<int> >, vec<pair_point> > > all_fullsp;
     for ( int i = 0; i < shb0.Pairs( ).isize( ); i++ )
     {    vec< vec<int> > fulls1, fulls2;
          vec<int> fulls_pos1, fulls_pos2;
          MapOne( shb0.PairLeft(i), cindex, 
               content, hb, to_right, fulls1, fulls_pos1 );
          MapOne( shb0.PairRight(i), cindex, 
               content, hb, to_right, fulls2, fulls_pos2 );

          // Get the start and stop positions of shb.PairLeft(i) and 
          // shb.PairRight(i) on the fulls.

          vec<int> left_start( fulls1.size( ), 0 ), right_start( fulls2.size( ), 0 );
          for ( int j1 = 0; j1 < fulls1.isize( ); j1++ )
          {    for ( int l = 0; l < fulls_pos1[j1]; l++ )
               {    left_start[j1] += shb0.EdgeLengthKmers( 
                         content[ fulls1[j1][0] ] [l] );    }    }
          for ( int j2 = 0; j2 < fulls2.isize( ); j2++ )
          {    for ( int l = 0; l < fulls_pos2[j2]; l++ )
               {    right_start[j2] += shb0.EdgeLengthKmers( 
                         content[ fulls2[j2][0] ] [l] );    }    }
          vec<int> left_stop( fulls1.size( ) ), right_stop( fulls2.size( ) );
          for ( int j1 = 0; j1 < fulls1.isize( ); j1++ )
          {    left_stop[j1] = left_start[j1];
               for ( int l = 0; l < shb0.PairLeft(i).isize( ); l++ )
               {    left_stop[j1] 
                         += shb0.EdgeLengthKmers( shb0.PairLeft(i)[l] );    }    }
          for ( int j2 = 0; j2 < fulls2.isize( ); j2++ )
          {    right_stop[j2] = right_start[j2];
               for ( int l = 0; l < shb0.PairRight(i).isize( ); l++ )
               {    right_stop[j2] 
                         += shb0.EdgeLengthKmers( shb0.PairRight(i)[l] );    }    }

          int nplaces = fulls1.size( ) * fulls2.size( );
          for ( int j1 = 0; j1 < fulls1.isize( ); j1++ )
          for ( int j2 = 0; j2 < fulls2.isize( ); j2++ )
          {    vec<pair_point> p = shb0.PairData(i);
               for ( int l = 0; l < p.isize( ); l++ )
               {    p[l].WeightMutable( ) /= nplaces;
                    int add = right_start[j2];
                    for ( int m = 0; m < fulls1[j1].isize( ); m++ )
                         add += EdgeLengthKmers( fulls1[j1][m] );
                    add -= left_stop[j1];
                    p[l].TrimMutable( ) += add;    }
               all_fullsp.push( make_pair( fulls1[j1], fulls2[j2] ), p );    }    }
     Sort(all_fullsp);
     for ( int i = 0; i < all_fullsp.isize( ); i++ )
     {    int j;
          for ( j = i + 1; j < all_fullsp.isize( ); j++ )
               if ( all_fullsp[j].first != all_fullsp[i].first ) break;
          usuy.push_back( all_fullsp[i].first );
          vec<pair_point> c;
          for ( int k = i; k < j; k++ )
               c.append( all_fullsp[k].second );
          Sort(c);
          county.push_back(c);
          i = j - 1;    }
     UniqueSortSync( usuy, county );    
     FixWeights(logc);
     REPORT_TIME( clock, "used bootstrapping" );    }

void SupportedHyperBasevector::writeBinary( BinaryWriter& writer ) const
{    HyperBasevector::writeBinary(writer);
     writer.write( Inv( ) );
     writer.write( Paths( ) );
     writer.write( WeightsFw( ) );
     writer.write( WeightsRc( ) );
     writer.write( WeightsFwOrigin( ) );
     writer.write( WeightsRcOrigin( ) );
     writer.write( Pairs( ) );
     writer.write( PairData( ) );
     writer.write( ReadCount( ) );
     writer.write( ReadLengthDist( ) );
     writer.write( FudgeMult( ) );    }

void SupportedHyperBasevector::readBinary( BinaryReader& reader )
{    HyperBasevector::readBinary(reader);
     reader.read( &InvMutable( ) );
     reader.read( &PathsMutable( ) );
     reader.read( &WeightsFwMutable( ) );
     reader.read( &WeightsRcMutable( ) );
     reader.read( &WeightsFwOriginMutable( ) );
     reader.read( &WeightsRcOriginMutable( ) );
     reader.read( &PairsMutable( ) );
     reader.read( &PairDataMutable( ) );
     reader.read( &ReadCountMutable( ) );
     reader.read( &ReadLengthDistMutable( ) );
     reader.read( &FudgeMultMutable( ) );    }

void SupportedHyperBasevector::DumpDot( const String& head,
     const vec<Bool>& invisible, const vec<String>& edge_color,
     const long_logging& logc, const Bool hide_inv, 
     const vec<String>& edge_names ) const
{    vec<Bool> hide;
     if (hide_inv) FlagEdgesForHiding( *this, Inv( ), hide, logc );
     else hide.resize( EdgeObjectCount( ), False );
     const Bool DOT_LABEL_CONTIGS = True;
     const Bool DOT_LABEL_VERTICES = False;
     vec<double> lengths( EdgeObjectCount( ) );
     for ( int i = 0; i < EdgeObjectCount( ); i++ )
          lengths[i] = EdgeLengthKmers(i);
     vec<String> edge_id_names( EdgeObjectCount( ) );
     if ( edge_names.size( ) > 0 ) edge_id_names = edge_names;
     else
     {    for ( int i = 0; i < EdgeObjectCount( ); i++ )
          {    if ( !InvDef(i) ) edge_id_names[i] = ToString(i);
               else 
               {    edge_id_names[i] = ToString(i) 
                         + "=" + ToString( Inv(i) ) + "'";    }    }    }
     Ofstream( dout, head + ".dot" );
     PrettyDOT( dout, lengths, HyperBasevector::edge_label_info(
          HyperBasevector::edge_label_info::DIRECT, &edge_id_names ),
          DOT_LABEL_CONTIGS, DOT_LABEL_VERTICES, NULL, NULL, NULL, &hide,
          &invisible, &edge_color, NULL, logc.LAYOUT );    }

void SupportedHyperBasevector::DumpFilesStandard( 
     const long_logging_control& log_control, const long_logging& logc,
     const int id ) const
{    if ( log_control.OUT_INT_HEAD != "" )
     {    DumpFiles( log_control.OUT_INT_HEAD + "." + ToString(id), 
               log_control, logc );    }    }

void SupportedHyperBasevector::DumpFiles( const String& head,
     const long_logging_control& log_control, const long_logging& logc ) const
{    double clock = WallClockTime( );
     
     // Output .shbv and .dot.

     BinaryWriter::writeFile( head + ".shbv", *this );
     vec<Bool> invisible( EdgeObjectCount( ), False );
     vec<String> edge_color( EdgeObjectCount( ), "" );
     DumpDot( head, invisible, edge_color, logc );

     // Output .support.

     Ofstream( sout, head + ".support" );
     for ( int i = 0; i < NPaths( ); i++ )
     {    int kmers = 0;
          for ( int j = 0; j < Path(i).isize( ); j++ )
               kmers += EdgeLengthKmers( Path(i)[j] );
          sout << "[" << i << "," << setiosflags(ios::fixed) << setprecision(1) 
               << Weight(i) << resetiosflags(ios::fixed) << "x] "
               << printSeq( Path(i) ) << " [" << kmers << " kmers]";
          // Temporary if.
          if ( WeightsFwOrigin( ).size( ) == WeightsFw( ).size( ) )
          {    vec<fix64_6> support;
               for ( int j = 0; j < WeightFwOrigin(i).isize( ); j++ )
                    support.push_back( WeightFwOrigin(i)[j].first );
               for ( int j = 0; j < WeightRcOrigin(i).isize( ); j++ )
                    support.push_back( WeightRcOrigin(i)[j].first );
               ReverseSort(support);
               sout << "; weights=";
               int count = 0;
               for ( int j = 0; j < support.isize( ); j++ )
               {    if ( j > 0 ) sout << ",";
                    if ( count++ > 5 ) 
                    {    sout << "...";
                         break;    }
                    int k = support.NextDiff(j);
                    sout << support[j] << "[" << k-j << "x]";
                    j = k - 1;    }    }
          sout << "\n";    }

     // Output .pair.

     Ofstream( pout, head + ".pairs" );
     for ( int i = 0; i < NPairs( ); i++ )
     {    vec<int> trims;
          for ( int j = 0; j < PairData(i).isize( ); j++ )
               trims.push_back( PairData(i,j).Trim( ) );
          Sort(trims);
          int Q1 = trims[ trims.size( ) / 4 ];
          int Q2 = trims[ trims.size( ) / 2 ];
          int Q3 = trims[ ( 3 * trims.size( ) ) / 4 ];
          pout << "\n[" << i << "] " << printSeqExp( PairLeft(i) ) << " --> "
               << printSeqExp( PairRight(i) ) << " (count=" << PairData(i).size( ) 
               << ", trim Q123 = " << Q1 << "," << Q2 << "," << Q3 << ")\n";    }
     pout << "\n";

     // Output .fasta.

     Ofstream( fout, head + ".fasta" );
     vec<int> to_left, to_right;
     ToLeft(to_left), ToRight(to_right);
     for ( int e = 0; e < EdgeObjectCount( ); e++ ) 
     {    int v = to_left[e], w = to_right[e];
          fout << ">edge_" << e << " " << v << ":" << w 
               << " K=" << K( ) << " kmers=" << EdgeLengthKmers(e) << "\n";
          EdgeObject(e).Print(fout);    }    

     if(logc.MAKE_FASTG) fastg::WriteFastg(head+".fastg", *this);

     // Output .glocs and .gpaths.

     if ( logc.USE_GENOME_FOR_DUMP && log_control.G->size( ) > 0 )
     {    Ofstream( out, head + ".glocs" );
          Ofstream( pout, head + ".gpaths" );

          // Compute placements of edges.

          vec< vec<placementy> > places( EdgeObjectCount( ) );
          #pragma omp parallel for
          for ( int e = 0; e < EdgeObjectCount( ); e++ ) 
               places[e] = log_control.FindGenomicPlacements( EdgeObject(e) );
          for ( int e = 0; e < EdgeObjectCount( ); e++ ) 
          {    int v = to_left[e], w = to_right[e];
               out << "edge_" << e << "[v=" << v << ",w=" << w << ",l=" 
                    << EdgeObject(e).size( ) << "]:";
               for ( int j = 0; j < places[e].isize( ); j++ )
               {    placementy& p = places[e][j];
                    int n = (*log_control.G)[p.g].size( );
                    if ( p.Pos < K( ) - 1 ) p.Pos += n;
                    out << " " << ( p.fw ? "fw" : "rc" ) << p.g << ".";
                    if (p.fw) out << p.pos << "-" << p.Pos;
                    else out << n - p.Pos << "-" << n - p.pos;    }
               out << "\n";    }    

          // Define a data structure P that assigns to each oriented edge the
          // vector of its genomic placements.

          int n = EdgeObjectCount( );
          vec< vec< triple<int,int,int> > > P(2*n);
          const vec<bool>& is_circular = *log_control.is_circular;
          for ( int e = 0; e < n; e++ )
          {    for ( int j = 0; j < places[e].isize( ); j++ )
               {    const placementy& p = places[e][j];
                    int start = p.pos;
                    int stop = start + EdgeObject(e).isize( );
                    int n = (*log_control.G)[p.g].size( );
                    if ( is_circular[p.g] && stop >= n + K( ) - 1 ) stop -= n;
                    stop -= ( K( ) - 1 );
                    P[ 2*e + ( p.fw ? 0 : 1 ) ].push( p.g, start, stop );    }    }

          // Build graph.

          // int N = 0;
          vec< triple<int,int,int> > V;
          vec<int> PS;
          PS.push_back(0);
          for ( int j = 0; j < P.isize( ); j++ )
          {    // N += P[j].size( );
               PS.push_back( PS.back( ) + P[j].isize( ) );
               for ( int l = 0; l < P[j].isize( ); l++ )
                    V.push( j, P[j][l].first, P[j][l].second );    }
          int N = PS.back( );
          vec< vec<int> > from(N), to(N);
          const int batch = 100;
          #pragma omp parallel for
          for ( int bz1 = 0; bz1 < this->N( ); bz1 += batch )
          {    vec< pair<int,int> > edges;
               for ( int z1 = bz1; z1 < Min( bz1 + batch, (int) this->N( ) ); z1++ )
               {    for ( int j1 = 0; j1 < From(z1).isize( ); j1++ )
                    {    int z2 = From(z1)[j1];
                         int v1 = EdgeObjectIndexByIndexFrom( z1, j1 );
                         if ( P[2*v1].empty( ) && P[1+2*v1].empty( ) ) continue;
                         for ( int j2 = 0; j2 < From(z2).isize( ); j2++ )
                         {    int z3 = From(z2)[j2];
                              int v2 = EdgeObjectIndexByIndexFrom( z2, j2 );
                              if ( P[2*v2].empty( ) && P[1+2*v2].empty( ) ) continue;
                              int x1, x2;
                              for ( int pass = 1; pass <= 2; pass++ )
                              {    if ( pass == 1 )
                                   {    x1 = 2*v1, x2 = 2*v2;    }
                                   else
                                   {    x2 = 2*v1 + 1, x1 = 2*v2 + 1;    }
                                   vec< triple< pair<int,int>, int, int > > X;
                                   for ( int l1 = 0; l1 < P[x1].isize( ); l1++ )
                                   {    X.push( make_pair( P[x1][l1].first,
                                             P[x1][l1].third ), 0, l1 );    }
                                   for ( int l2 = 0; l2 < P[x2].isize( ); l2++ )
                                   {    X.push( make_pair( P[x2][l2].first,
                                             P[x2][l2].second ), 1, l2 );    }
                                   Sort(X);
                                   for ( int i = 0; i < X.isize( ); i++ )
                                   {    int j, mid;
                                        for ( j = i + 1; j < X.isize( ); j++ )
                                             if ( X[j].first != X[i].first ) break;
                                        for ( mid = i + 1; mid < j; mid++ )
                                             if ( X[mid].second == 1 ) break;
                                        for ( int u1 = i; u1 < mid; u1++ )
                                        for ( int u2 = mid; u2 < j; u2++ )
                                        {    int l1 = X[u1].third, l2 = X[u2].third;
                                             int m1 = l1, m2 = l2;
                                             m1 += PS[x1];
                                             m2 += PS[x2];
                                             edges.push( m1, m2 );    }    
                                        i = j - 1;    }    }    }    }    }    
               #pragma omp critical
               {    for ( int i = 0; i < edges.isize( ); i++ )
                    {    from[ edges[i].first ].push_back( edges[i].second );
                         to[ edges[i].second ].push_back( 
                              edges[i].first );    }    }   }
          #pragma omp parallel for
          for ( int i = 0; i < N; i++ )
          {    Sort(from[i]), Sort(to[i]);    }
          digraph H( from, to );
          vec< triple< int, pair<int,int>, vec<int> > > matches;
          vec< vec<int> > paths;
          H.AllPaths( -1, -1, paths );

          if ( Sum(is_circular) > 0 ) // wrong, need granular test
          {    vec< vec<int> > all_loops;
               for ( int v = 0; v < N; v++ )
               {    vec< vec<int> > loops;
                    H.AllPaths( v, v, loops, -1, True );
                    for ( int j = 0; j < loops.isize( ); j++ )
                    {    vec<int> x = loops[j];
                         if ( x.solo( ) ) continue;
                         int m = Min(x);
                         int p;
                         for ( p = 0; p < x.isize( ); p++ )
                              if ( x[p] == m ) break;
                         vec<int> y;
                         for ( int i = p; i < x.isize( ) - 1; i++ )
                              y.push_back( x[i] );
                         for ( int i = 0; i < p; i++ )
                              y.push_back( x[i] );
                         if ( !Member( all_loops, y ) ) 
                              all_loops.push_back(y);    }    }
               paths.append(all_loops);    }

          for ( int j = 0; j < paths.isize( ); j++ )
          {    const vec<int>& p = paths[j];
               int g = V[ p.front( ) ].second, pos = V[ p.front( ) ].third;
               int Pos = V[ p.back( ) ].third
                    + EdgeObject( V[ p.back( ) ].first / 2 ).isize( );
               vec<int> q;
               for ( int l = 0; l < paths[j].isize( ); l++ )
               {    int e = V[ paths[j][l] ].first;
                    if ( e % 2 == 0 ) q.push_back(e/2);
                    else q.push_back( -e/2-1 );    }
               matches.push( g, make_pair( pos, Pos ), q );    }
          Sort(matches);
          for ( int j = 0; j < matches.isize( ); j++ )
          {    int g = matches[j].first;
               int n = (*log_control.G)[g].isize( );
               int start = matches[j].second.first, stop = matches[j].second.second;
               if ( is_circular[g] && stop >= n + K( ) - 1 ) stop -= n;
               pout << g << "." << start << "-" << stop;
               if ( stop - start == K( ) - 1 ) pout << " (perfect circle)";
               pout << " :: ";
               const vec<int>& v = matches[j].third;
               for ( int l = 0; l < v.isize( ); l++ )
               {    if ( l > 0 ) pout << ",";
                    if ( v[l] >= 0 ) pout << v[l] << "fw";
                    else pout << -v[l]-1 << "rc";    }
               pout << "\n";    }    }

     REPORT_TIME( clock, "used dumping files" );    }

void SupportedHyperBasevector::UnrollLoops( const long_logging& logc )
{    double clock = WallClockTime( );
     if (logc.STATUS_LOGGING) cout << Date( ) << ": unrolling loops" << endl;
     vec< pair< triple<int,int,int>, vec<int> > > edits;
     vec<int> to_left, to_right;
     ToLeft(to_left), ToRight(to_right);

     // Define proposed edits.

     vec< pair< triple<int,int,int>, vec<int> > > plans;
     for ( int v = 0; v < N( ); v++ )
     {    if ( To(v).size( ) != 2 || From(v).size( ) != 2 ) continue;
          int i1 = 0, i2 = 1, j1 = 0, j2 = 1;
          if ( To(v)[i1] != v ) swap(i1, i2);
          if ( To(v)[i1] != v ) continue;
          if ( From(v)[j1] != v ) swap(j1, j2);
          if ( From(v)[j1] != v ) continue;
          int u = To(v)[i2], w = From(v)[j2];
          if ( u == v || u == w || v == w ) continue;
          int e = EdgeObjectIndexByIndexTo( v, i2 );
          int f = EdgeObjectIndexByIndexTo( v, i1 );
          int g = EdgeObjectIndexByIndexFrom( v, j2 );

          // Now have loop f sandwiched between e on the left and g on the right.

          if ( EdgeLengthKmers(f) > MedianCorrectedReadLengthFudge( ) ) continue;
          vec<int> counts, counts2;
          for ( int i = 0; i < NPaths( ); i++ )
          {    const vec<int>& p = Path(i);
               for ( int j = 0; j < p.isize( ); j++ )
               {    if ( p[j] == e )
                    {    for ( int k = j + 1; k < p.isize( ); k++ )
                         {    if ( p[k] == g )
                              {    counts.push_back( k - j - 1 );
                                   break;    }    }    }    }
               for ( int j = 0; j < p.isize( ); j++ )
               {    if ( p[j] == f )
                    {    int k;
                         for ( k = j + 1; k < p.isize( ); k++ )
                              if ( p[k] != f ) break;
                         counts2.push_back( k - j );
                         j = k - 1;    }    }    }
          UniqueSort(counts), UniqueSort(counts2);
          if ( counts.empty( ) ) continue;
          if ( counts2.nonempty( ) && counts2.back( ) > counts.back( ) ) continue;
          plans.push( make_triple(e,f,g), counts );    }

     // Process proposed edits.

     vec<int> used( EdgeObjectCount( ), False );
     for ( int j = 0; j < plans.isize( ); j++ )
     {    int e = plans[j].first.first, f = plans[j].first.second;
          int g = plans[j].first.third;
          if ( used[e] || used[f] || used[g] ) continue;

          // Don't handle the case where {e,f,g} meets {re,rf,rg}.

          int re = Inv(g), rf = Inv(f), rg = Inv(e);
          if ( Meet2( vec<int>({e,f,g}), vec<int>({re,rf,rg}) ) ) continue;

          used[e] = used[f] = used[g] = True;
          if ( re >= 0 ) used[re] = used[rf] = used[rg] = True;
          const vec<int>& counts = plans[j].second;
          int npasses = ( re < 0 ? 1 : 2 );
          for ( int pass = 1; pass <= npasses; pass++ )
          {    int ex, fx, gx;
               if ( pass == 1 ) { ex = e, fx = f, gx = g; }
               else { ex = re, fx = rf, gx = rg; }
               int ux = to_left[ex], wx = to_right[gx];
               if ( logc.verb[ "UNROLL_LOOPS" ] >= 1 )
               {    cout << "replacing loop " << fx << " by " << printSeq(counts)
                         << " copies of it" << endl;    }
               vec<int> new_ids, dels( {ex,fx,gx} );
               DeleteEdges(dels);
               vec< vec<int> > x;
               for ( int i = 0; i < counts.isize( ); i++ )
               {    new_ids.push_back( EdgeObjectCount( ) );
                    vec<int> q( {ex} );
                    for ( int j = 0; j < counts[i]; j++ )
                         q.push_back(fx);
                    q.push_back(gx);
                    x.push_back(q);
                    AddEdge( ux, wx, Cat(q) );    }
               edits.push( make_triple(ex,fx,gx), new_ids );
               TransformPaths( x, new_ids );    }    }

     // Clean up.

     InvMutable( ).resize( EdgeObjectCount( ), -1 );
     Sort(edits);
     vec< triple<int,int,int> > edits1( edits.size( ) );
     for ( int i = 0; i < edits.isize( ); i++ )
          edits1[i] = edits[i].first;
     for ( int i = 0; i < edits.isize( ); i++ )
     {    int e = edits[i].first.first, f = edits[i].first.second; 
          int g = edits[i].first.third;
          int re = Inv(g), rf = Inv(f), rg = Inv(e);
          if ( re < 0 ) continue;
          int p = BinPosition( edits1, make_triple( re, rf, rg ) );
          ForceAssertGe( p, 0 );
          for ( int j = 0; j < edits[i].second.isize( ); j++ )
               InvMutable( edits[i].second[j] ) = edits[p].second[j];    }
     RemoveUnneededVertices( );
     RemoveEdgelessVertices( );
     REPORT_TIME( clock, "used unrolling loops" );
     RemoveDeadEdgeObjects( );
     if (logc.STATUS_LOGGING)
     {    cout << Date( ) << ": after unrolling loops there are "
               << EdgeObjectCount( ) << " edges" << endl;    }
     TestValid(logc);    }

void SupportedHyperBasevector::DeleteWeakEdges( const long_logging& logc )
{    double clock = WallClockTime( );
     int verbosity = logc.verb[ "DELETE_WEAK" ];
     if ( verbosity >= 1 )
     {    cout << Date( ) << ": begin deleting weak edges" << endl;
          int64_t npe = 0;
          for ( int i = 0; i < NPaths( ); i++ )
               npe += Path(i).size( );
          cout << Date( ) << ": edges = " << ToStringAddCommas( EdgeObjectCount( ) )
               << "; paths = " << ToStringAddCommas( NPaths( ) )
               << "; segs = " << ToStringAddCommas(npe) << endl;    }
 
     // Define heuristics.

     const int min_mult = 100;
     const int max_kill = 10;

     // Execute multiple passes.

     while(1)
     {
          // Define adjacencies in paths.

          if ( verbosity >= 1 ) cout << Date( ) << ": defining adjacencies" << endl;
          vec< triple<int,int,fix64_6> > adj, adjb;
          int64_t ne = EdgeObjectCount( ), np = NPaths( );
          vec<int64_t> ind( ne + 1, -1 ), indb( ne + 1, -1 );
          {    const int64_t batches = 100;
               vec< vec< triple<int,int,fix64_6> > > adjx(batches);
               #pragma omp parallel for
               for ( int64_t b = 0; b < batches; b++ )
               {    for ( int64_t i = (b*np)/batches; i < ((b+1)*np)/batches; i++ )
                    {    const vec<int>& p = Path(i);
                         for ( int j = 1; j < p.isize( ); j++ )
                              adjx[b].push( p[j-1], p[j], Weight(i) );    }    }
               if ( verbosity >= 1 ) cout << Date( ) << ": appending" << endl;
               vec<int64_t> start( batches + 1, 0 );
               for ( int b = 0; b < batches; b++ )
                    start[b+1] = start[b] + adjx[b].jsize( );
               adj.resize( start.back( ) );
               #pragma omp parallel for
               for ( int64_t b = 0; b < batches; b++ )
               {    if ( adjx[b].nonempty( ) )
                    {    memcpy( &( adj[ start[b] ] ), &( adjx[b][0] ), 
                              adjx[b].jsize( ) 
                                   * sizeof( triple<int,int,fix64_6> ) );    }    }
               if ( verbosity >= 1 )
                    cout << Date( ) << ": adding trivial adjacencies" << endl;
               vec<int> to_right;
               ToRight(to_right);
               for ( int e = 0; e < ne; e++ )
               {    int w = to_right[e];
                    for ( int j = 0; j < From(w).isize( ); j++ )
                         adj.push( e, EdgeObjectIndexByIndexFrom( w, j ), 0 );    }
               ParallelSort(adj);
               vec<Bool> to_delete( adj.size( ), False );
               for ( int64_t i = 0; i < adj.jsize( ); i++ )
               {    int64_t j;
                    for ( j = i + 1; j < adj.jsize( ); j++ )
                    {    if ( adj[j].first != adj[i].first 
                              || adj[j].second != adj[i].second )
                         {    break;    }    }
                    for ( int64_t k = i + 1; k < j; k++ )
                    {    adj[i].third += adj[k].third;
                         to_delete[k] = True;    }
                    i = j - 1;    }
               EraseIf( adj, to_delete );
               adjb = adj;
               for ( int64_t i = 0; i < adjb.jsize( ); i++ )
                    swap( adjb[i].first, adjb[i].second );
               ParallelSort(adjb);
               ind[ne] = indb[ne] = adj.size( );
               for ( int64_t i = 0; i < adj.jsize( ); i++ )
               {    ind[ adj[i].first ] = i;
                    indb[ adjb[i].first ] = i;    }
               for ( int64_t i = ne - 1; i >= 0; i-- )
               {    if ( ind[i] < 0 ) ind[i] = ind[i+1];
                    if ( indb[i] < 0 ) indb[i] = indb[i+1];    }    }
     
          // Look for edges to delete.
     
          vec<int> dels;
          for ( int e = 0; e < ne; e++ )
          {    fix64_6 me = 0, meb = 0;
               for ( int64_t i = ind[e]; i < ind[e+1]; i++ )
                    me = Max( me, adj[i].third );
               for ( int64_t i = ind[e]; i < ind[e+1]; i++ )
               {    if ( min_mult * Max( (fix64_6) 1, adj[i].third ) > me ) continue;
                    if ( adj[i].third > max_kill ) continue;
                    int f = adj[i].second;
                    fix64_6 mf = 0;
                    for ( int64_t j = ind[f]; j < ind[f+1]; j++ )
                         mf = Max( mf, adj[j].third );
                    if ( mf > max_kill ) continue;
                    fix64_6 mfb = 0;
                    for ( int64_t j = indb[f]; j < indb[f+1]; j++ )
                         mfb = Max( mfb, adjb[j].third );
                    if ( mfb > max_kill ) continue;
                    dels.push_back(f);    }
               for ( int64_t i = indb[e]; i < indb[e+1]; i++ )
                    meb = Max( meb, adjb[i].third );
               for ( int64_t i = indb[e]; i < indb[e+1]; i++ )
               {    if ( min_mult * Max( (fix64_6) 1, adjb[i].third ) > me ) continue;
                    if ( adjb[i].third > max_kill ) continue;
                    int f = adjb[i].second;
                    fix64_6 mf = 0;
                    for ( int64_t j = ind[f]; j < ind[f+1]; j++ )
                         mf = Max( mf, adj[j].third );
                    if ( mf > max_kill ) continue;
                    fix64_6 mfb = 0;
                    for ( int64_t j = indb[f]; j < indb[f+1]; j++ )
                         mfb = Max( mfb, adjb[j].third );
                    if ( mfb > max_kill ) continue;
                    dels.push_back(f);    }    }

          // Delete edges and clean up.

          int ndels = dels.size( );
          if ( ndels == 0 ) break;
          for ( int i = 0; i < ndels; i++ )
               if ( Inv( dels[i] ) >= 0 ) dels.push_back( Inv( dels[i] ) );
          UniqueSort(dels);
          if ( verbosity >= 1 )
          {    cout << Date( ) << ": deleting " << ToStringAddCommas(ndels) 
                    << " edges" << endl;    }
          DeleteEdges(dels);

          // Clean up.

          vec<Bool> delsx( ne, False );
          for ( int i = 0; i < dels.isize( ); i++ )
               delsx[ dels[i] ] = True;
          vec<Bool> p_to_delete( NPaths( ), False );
          for ( int i = 0; i < NPaths( ); i++ )
          {    for ( int j = 0; j < Path(i).isize( ); j++ )
               {    if ( delsx[ Path(i)[j] ] ) 
                    {    p_to_delete[i] = True;
                         break;    }    }    }
          EraseIf( PathsMutable( ), p_to_delete );
          // Temporary if.
          if ( WeightsFwOrigin( ).size( ) == WeightsFw( ).size( ) )
          {    EraseIf( WeightsFwOriginMutable( ), p_to_delete );
               EraseIf( WeightsRcOriginMutable( ), p_to_delete );    }
          EraseIf( WeightsFwMutable( ), p_to_delete );
          EraseIf( WeightsRcMutable( ), p_to_delete );
          p_to_delete.resize_and_set( NPairs( ), False );
          for ( int i = 0; i < NPairs( ); i++ )
          {    Bool OK = True;
               for ( int j = 0; j < PairLeft(i).isize( ); j++ )
                    if ( BinMember( dels, PairLeft(i)[j] ) ) OK = False;
               if ( !OK ) p_to_delete[i] = True;    
               for ( int j = 0; j < PairRight(i).isize( ); j++ )
                    if ( BinMember( dels, PairRight(i)[j] ) ) OK = False;
               if ( !OK ) p_to_delete[i] = True;    }
          EraseIf( PairsMutable( ), p_to_delete );
          EraseIf( PairDataMutable( ), p_to_delete );
          RemoveEdgelessVertices( );
          if ( verbosity >= 1 ) cout << Date( ) << ": removing unneeded verts" << endl;
          RemoveUnneededVertices( );
          if ( verbosity >= 1 ) cout << Date( ) << ": tail" << endl;
          RemoveDeadEdgeObjects( );
          if ( verbosity >= 1 )
          {    int64_t npe = 0;
               for ( int i = 0; i < NPaths( ); i++ )
                    npe += Path(i).size( );
               cout << Date( ) << ": edges = " 
                    << ToStringAddCommas( EdgeObjectCount( ) ) << "; paths = " 
                    << ToStringAddCommas( NPaths( ) )
                    << "; segs = " << ToStringAddCommas(npe) << endl;    }    }
     REPORT_TIME( clock, "used deleting weak edges" ); // probably double counting
     TestValid(logc);    }
