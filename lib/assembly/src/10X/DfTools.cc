// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"
#include "Intvector.h"
#include "ParallelVecUtilities.h"
#include "feudal/ObjectManager.h"
#include "feudal/PQVec.h"
#include "kmers/KmerRecord.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "paths/long/large/ExtractReads.h"
#include "10X/paths/ReadPathVecX.h"
#include "10X/DfTools.h"
#include "10X/Gap.h"
#include "10X/mergers/NicePrints.h"

size_t StatLogger::externalSizeof() { return 0; }
StatLogger StatLogger::gInst;


void InitializePathsXFromPaths( ReadPathVecX & pathsX, const HyperBasevectorX & hb,
     const String paths_file, const size_t BATCH, const Bool verbose )
{
     double clock = WallClockTime();
     double iotime = 0;
     
     pathsX.clear();

     if ( verbose ) cout << Date( ) << ": multi-threaded in pieces" << endl;
     size_t N = 0;
     { 
          VirtualMasterVec<ReadPath> vpaths ( paths_file );
          N = vpaths.size();
     }
     {
          int batches = N/BATCH+1;
          ReadPathVec piece;
          const int nthreads = omp_get_max_threads();
          ReadPathVecX xpiece;
          int batch = 0;
          size_t lenSum = 0;
          for ( size_t id = 0; id < N; id+=BATCH ) {
               batch++;
               if ( verbose )
                    cout << Date( ) << ": reading from disk (" << batch << " of " 
                         << batches << ")" << endl;
               const size_t stop = Min( id + BATCH, N );
               piece.clear();
               double ioclock = WallClockTime();
               piece.ReadRange( paths_file, id, stop );
               iotime += (WallClockTime() - ioclock);
               if ( verbose )
                    cout << Date( ) <<": parallel append" << endl;
               const int npaths = piece.size()/nthreads + 1;
               #pragma omp parallel for ordered num_threads(nthreads) \
                    firstprivate( xpiece )
               for ( int t = 0; t < nthreads; t++ ) {
                    xpiece.clear();
                    const int64_t start = t*npaths;
                    const int64_t num_reads = Min( (size_t)(t+1)*npaths, piece.size() )-start;
                    xpiece.append( piece, hb, start, num_reads );
                    #pragma omp ordered
                    {    pathsX.append( xpiece ); }
               }
               lenSum = std::accumulate( piece.begin(), piece.end(), lenSum, [](size_t b, ReadPath const& a){ return b+a.size(); } );
          }

          cout << "Total readpath length = " << lenSum << ", average=" << static_cast<long double>(lenSum)/N << ", count=" << N << endl;
     }
     int ioper = iotime/(WallClockTime()-clock ) * 100;
     cout << Date( ) << ": finished reading pathsX, used " << TimeSince(clock)
          << "; " << ioper << " \% I/O" << endl;

     // calculate average path length
}

// Create a vector of integers, one for each read, such that "having two"
// nonzero elements is enough.  Note nested parallel for loops, not really
// what we want.

void Computebid( const vec<DataSet> & datasets, const vec<int64_t> & bci,
     vec<int64_t> & bid )
{
     cout << Date( ) << ": creating bid" << endl;
     int64_t N = bci.back();
     bid.resize( N );
     #pragma omp parallel for schedule( dynamic, 1 )
     for ( int b = 0; b < bci.isize( ) - 1; b++ )
     {    int64_t start = bci[b], stop = bci[b+1];
          #pragma omp parallel for
          for ( int64_t id = start; id < stop; id++ )
          {    int di;
               for ( di = 0; di < datasets.isize( ); di++ )
                    if ( id < datasets[di].start ) break;
               const ReadDataType& dtype = datasets[di-1].dt;
               if ( dtype == ReadDataType::BAR_10X )
                    bid[id] = N + b + 1;
               else if ( dtype == ReadDataType::UNBAR_10X ) bid[id] = 0;
               else if ( dtype == ReadDataType::PCR_FREE ) bid[id] = id + 1;

               // Probably not what we want:

               else if ( dtype == ReadDataType::PCR ) 
                    bid[id] = id + 1;    }    }
}

int64_t bidc( const vec<DataSet> & datasets, const vec<int> & bc,
     const int64_t & id )
{
     const int64_t N = bc.size( );
     // First test for 10X barcoded read
     if ( bc[id] ) return (N + bc[id]);
     
     // this is barcode zero:

     // do we only have 10X data
     if ( datasets[0].dt == ReadDataType::UNBAR_10X || 
          datasets[0].dt == ReadDataType::BAR_10X )
          return N;
     
     // if we have PCR free data
     if ( id < datasets[1].start )
          return id;
     else
          return N;
}

vec<int> GetBarcodes( const int e, const vec<int>& inv,
     const VecULongVec& paths_index, const vec<int>& bc )
{    vec<int> bs;
     for ( int x : { e, inv[e] } )
     {    for ( int j = 0; j < (int) paths_index[x].size( ); j++ )
          {    int64_t id = paths_index[x][j];
               if ( bc[id] > 0 ) bs.push_back( bc[id] );    }    }
     UniqueSort(bs);
     return bs;    }

void WriteSubSample( vecbasevector const& bases, ObjectManager<VecPQVec>& quals_om,
          size_t pair_count, String const& head )
{
     // clamp pair_count to available data
     pair_count = std::min( bases.size()/2, pair_count);
     auto orig_pair_count = pair_count;

     double frac = static_cast<double>(pair_count)/(bases.size()/2);
     auto decider = [frac]() {
          return (1. * randomx() / RNGen::RNGEN_RAND_MAX) <= frac; };

     vecbasevector sbases;
     VecPQVec squals;
     sbases.reserve(pair_count*2);
     squals.reserve(pair_count*2);

     auto const& quals = quals_om.load();


     for ( size_t i = 0; i < bases.size()-1 && pair_count; i+=2 ) {
          // pull randomly until there are only enough left to fill
          // just fudge so that we end up with precisely the right num
          // (or less if simply not enough data)
          if ( decider() || pair_count*2 >= bases.size()-i ) {
               sbases.push_back( bases[i] );
               sbases.push_back( bases[i+1] );
               squals.push_back( quals[i] );
               squals.push_back( quals[i+1] );
               pair_count--;
          }
     }

     ForceAssertEq( orig_pair_count*2, sbases.size() );
     sbases.WriteAll( head + ".fastb" );
     squals.WriteAll( head + ".qualp" );
}

void LoadData( const String& work_dir, const String& R, const vec<String>& lr,
     const vec<double>& LR_SELECT_FRAC, vecbasevector& bases,
     ObjectManager<VecPQVec>& quals_om, vec<int64_t>& bci,
     vec<String>& subsam_names, vec<int64_t>& subsam_starts, vec<DataSet>& datasets )
{
     String SAMPLE, species;
     String tmp_dir1 = work_dir + "/data";
     auto& quals = quals_om.create();

     // Load ordinary read data.

     if ( R.size( ) > 0 ) {
          FatalErr("regular unbarcoded data not accepted");
          // NW - I broke this when streaming to disk.  I can resurrect, but 
          // one idea would be the stream bases and quals to disk via bases_out and quals_out and
          // then continue onward (with nbases set correctly). 
#if 0
          vec<String> regions;
          String SELECT_FRAC = "1";
          ExtractReads( SAMPLE, species, R, SELECT_FRAC, -1,
               regions, tmp_dir1, work_dir, False, False, False,
               subsam_names, subsam_starts, &bases, quals_om );
          // for now, we're assuming R is a single, PCR-free dataset
          datasets.push_back( { ReadDataType::PCR_FREE, 0 } );
#endif
     }
    
     PRINTDEETS("load reads");
     String const& OUT_HEAD = work_dir + "/data/frag_reads_orig";
     
     if ( lr.size() == 1 ) {
        // if it's only one file, just read it
        String head = lr[0].Before(".fastb");
        PRINTDEETS("short pipeline path, reading from " << head);
        bases.ReadAll(head+".fastb");
        quals.ReadAll(head+".qualp");
        BinaryReader::readFile(head+".bci", &bci);

        datasets.push_back( { ReadDataType::UNBAR_10X, bci[0] } );   // UNBAR is [ bci[0], bci[1] )
        if ( bci.size() > 2 ) {                                      // BAR is [ bci[1], bci[2] )
             // something rather wrong if we don't get here...
             datasets.push_back( { ReadDataType::BAR_10X, bci[1] } );
        }
        
        // and write it
        PRINTDEETS("writing reads/quals/bci");
        bases.WriteAll(OUT_HEAD + ".fastb");
        ForceAssertEq(quals_om.filename(), OUT_HEAD+".qualp");
        quals_om.store();
        BinaryWriter::writeFile(OUT_HEAD+".bci", bci);
     } else {
         
         BinaryIteratingWriter<vec<int64_t>>     bci_out( OUT_HEAD+".bci" );
         IncrementalWriter<basevector>           bases_out( OUT_HEAD+".fastb" );
         IncrementalWriter<PQVec>                quals_out( OUT_HEAD+".qualp" );

         bci_out.write(0u);	// null barcode always at the front


         // load linked-read (LR) data in two passes to avoid duplicating
         // reads in memory -- first unbarcoded, then barcoded
         //
         if ( lr.size() ) ForceAssertEq(lr.size(), LR_SELECT_FRAC.size() );
         cout << Date() << ": reading in linked read data" << endl;
         int64_t nbases=0;
         enum passes:int { PASS_UNBARCODED, PASS_BARCODED, PASS_LAST };
         for ( int pass = PASS_UNBARCODED; pass < PASS_LAST; ++pass ) {

              for ( size_t i = 0; i < lr.size(); ++i ) {

                   String head = lr[i].Before(".fastb");
                   VirtualMasterVec<basevector> basesi( head + ".fastb" );
                   VirtualMasterVec<PQVec> qualsi( head + ".qualp" );

                   vec<int64_t> bcii;
                   BinaryReader::readFile( head + ".bci" , &bcii );

                   // sanity check .bci to avoid old code
                   if ( bcii[0] != 0 )
                        FatalErr("barcode 0 is unbarcoded data and must start at 0");

                   auto frac = LR_SELECT_FRAC[i];
                   auto decider = [frac]() {
                        return (1. * randomx() / RNGen::RNGEN_RAND_MAX) <= frac;
                   };

                   ForceAssertEq( bcii[1] % 2, 0 );
                   ForceAssertLe( bcii[1], basesi.size() );
                   ForceAssertLe( bcii[1], qualsi.size() );

                   // TODO: these two cases below can now be merged easily
                   if ( pass == PASS_UNBARCODED ) {
                        datasets.push_back( { ReadDataType::UNBAR_10X, nbases } );

                        // append bases, quals from [ 0, bcii[1] )
                        // assume now that we only have pairs
                        for ( size_t i = 0; i < (size_t) bcii[1]; i+=2 ) {
                             if ( decider() ) {
                                  bases_out.add(basesi[i]);
                                  bases_out.add(basesi[i+1]);
                                  quals_out.add(qualsi[i]);
                                  quals_out.add(qualsi[i+1]);
                                  nbases += 2;
                             }
                        }

                   } else if ( pass == PASS_BARCODED ) {
                        datasets.push_back( { ReadDataType::BAR_10X, nbases } );

                        for ( size_t bc = 1; bc < bcii.size()-1; ++bc ) {
                             // for each barcode
                             bci_out.write( nbases );
                             for ( size_t i = (size_t) bcii[bc]; i < (size_t) bcii[bc+1]; i+=2 ) {
                                  // for each pair of this barcode
                                  if ( decider() ) {
                                       bases_out.add(basesi[i]);
                                       bases_out.add(basesi[i+1]);
                                       quals_out.add(qualsi[i]);
                                       quals_out.add(qualsi[i+1]);
                                       nbases += 2;
                                  }
                             }
                        }

                   } else FatalErr("Bug - barcodes pass");

              }	// for lr...
         } // for pass...

         bci_out.write( nbases );

         bci_out.close();
         bases_out.close();
         quals_out.close();
         
         bases.ReadAll( OUT_HEAD + ".fastb" );
         ForceAssertEq(quals_om.filename(), OUT_HEAD+".qualp");
         quals_om.unload(); quals_om.load();
         ForceAssertEq( bases.size( ), quals.size( ) );
         BinaryReader::readFile( OUT_HEAD + ".bci", &bci );
     }



     PRINTDEETS("done loading reads");
     WriteSubSample( bases, quals_om, 500, work_dir + "/data/frag_reads_orig.1000" );
}

void GetQualStats( const VecPQVec& quals, vec<vec<vec<int64_t>>>& hist, 
                   int & max_read_length )
{ 
     const int64_t N = quals.size( );
     const int T = omp_get_max_threads( );
     const int64_t batch = N/T + 1;

     // compute the histogram
     // hist is a 2 x max_read_length x 256 histogram
     cout << Date( ) << ": computing quality histogram" << endl; 
     
     // create the data structure
     hist.resize( 2, vec<vec<int64_t>>( max_read_length, vec<int64_t>(256, 0) ) );
     
     { // open block to kill temp data structure

     // create local copies of structure for each thread
     vec<vec<vec<vec<int64_t>>>> hist_part( T, hist );
     #pragma omp parallel for num_threads(T)
     for ( int t = 0; t < T; t++ ) {
          const int64_t start = batch*t;
          const int64_t stop  = Min( batch*(t+1), N );
          for ( int64_t id = start; id < stop; id++ ) {
               qualvector q;
               quals[id].unpack( &q );
               int pos = 0;
               for ( auto x : q ) {
                    hist_part[t][id%2][pos][x]++;
                    pos++;
               }
          }
     }
     
     // combine parts
     for ( int t = 0; t < T; t++ ) {
          #pragma omp parallel for collapse(3)
          for ( int read = 0; read < 2; read++ ) {
               for ( int pos = 0; pos < max_read_length; pos++ ) {
                    for ( int q = 0; q < 256; q++ ) {
                         hist[read][pos][q] += hist_part[t][read][pos][q];
                    }
               }
          }
     }

     } // end of computation block

     // find max quality in data structure.
     int max_qual;
     bool found_max=false;
     for (max_qual = 255; max_qual >= 0; max_qual--) {
          for (int pos = 0; pos != max_read_length; pos++) {
               if ( hist[0][pos][max_qual] > 0 || hist[1][pos][max_qual] > 0 ) {
                    found_max = true;
                    break;
               }
          }
          if ( found_max )
               break;
     }
     // resize data structure
     for ( int pos = 0; pos != max_read_length; pos++ ) {
          hist[0][pos].resize(max_qual+1);
          hist[1][pos].resize(max_qual+1);
     }
     
}

void FragDist( const HyperBasevectorX& hb, const vec<int>& inv,
     const ReadPathVec& paths, vec<int64_t>& count )
{    const int max_sep = 1000;
     count.resize( max_sep + 1, 0 );
     #pragma omp parallel for
     for ( int64_t id1 = 0; id1 < (int64_t) paths.size( ); id1 += 2 )
     {    int64_t id2 = id1 + 1;
          if ( paths[id1].size() == 0 || paths[id2].size() == 0 ) continue;
          const ReadPath& p1 = paths[id1];
          const ReadPath& p2 = paths[id2];
          int e1 = p1[0], e2 = inv[ p2[0] ];
          int epos1 = paths[id1].getOffset( );
          if ( e1 != e2 ) continue;
          int n = hb.Bases(e1);
          if ( epos1 + max_sep > n ) continue;
          int epos2 = hb.Bases(e2) - paths[id2].getOffset( );
          int len = epos2 - epos1;
          if ( len < 0 || len > max_sep ) continue;
          #pragma omp critical
          {    count[len]++;    }    }    }

void FragDist( const HyperBasevectorX& hb, const vec<int>& inv,
     const ReadPathVecX& paths, vec<int64_t>& count )
{    const int max_sep = 1000;
     count.resize( max_sep + 1, 0 );
     #pragma omp parallel for
     for ( int64_t id1 = 0; id1 < (int64_t) paths.size( ); id1 += 2 )
     {    int64_t id2 = id1 + 1;
          if ( paths.getNumEdges(id1) == 0 || paths.getNumEdges(id2) == 0 ) continue;
          ReadPath p1, p2; paths.unzip(p1,hb,id1); paths.unzip(p2,hb,id2);
          int e1 = p1[0], e2 = inv[ p2[0] ];
          int epos1 = paths.getOffset( id1 );
          if ( e1 != e2 ) continue;
          int n = hb.Bases(e1);
          if ( epos1 + max_sep > n ) continue;
          int epos2 = hb.Bases(e2) - paths.getOffset( id2 );
          int len = epos2 - epos1;
          if ( len < 0 || len > max_sep ) continue;
          #pragma omp critical
          {    count[len]++;    }    }    }

void ReadTwoPctProper( const HyperBasevectorX& hb, const vec<int>& inv,
     ReadPathVecX const& paths, double & r2_pct_proper )
{
     cout << Date( ) << ": computing read two pct proper";
     const int max_sep = 1000;
     double r2_map = 0.0;
     int64_t sample = 0;
     ReadPath p1, p2;
     #pragma omp parallel for reduction(+:r2_map,sample) private(p1,p2)
     for ( int64_t id1 = 0; id1 < (int64_t) paths.size( ); id1 += 2 ) {
          int64_t id2 = id1 + 1;
          // new stuff here
          // R1 must be mapped to edge e1
          if ( paths.getNumEdges( id1 ) == 0 ) continue;
          paths.unzip(p1,hb,id1);
          int e1 = p1[0];
          int le1 = hb.Bases(e1);
          int epos1 = paths.getOffset( id1 );
          
          // is the edge long enough for R2 to map here
          if ( le1 - epos1 >= max_sep ) {
               // did R2 map?
               bool mapped = false;
               if ( paths.getNumEdges( id2 ) > 0 ) 
               {
                    paths.unzip(p2,hb,id2);
                    int e2 = inv[ p2[0] ];
                    if ( e1 == e2 ) {
                         //R1,R2 mapped to the same edge
                         int epos2 = hb.Bases(e2) - paths.getOffset( id2 );
                         int len = epos2 - epos1;
                         // make sure the insert isn't too long
                         if ( len >= 0 && len <= max_sep )
                              mapped=true;
                    }
               }
               if (mapped)
                    r2_map++;
               sample++;
          }
     }
     cout << " [sample size " << sample << "]" << endl;
     if (sample > 0)
          r2_pct_proper = (r2_map/sample)*100; // This number is a percentage!
     else
          r2_pct_proper = 0;

}
template<int K> void MapClosures(
     const HyperBasevectorX& hb, const vec<basevector>& closures, ReadPathVec& clop )
{    cout << Date( ) << ": building closure lookup" << endl;

     // Remove short closures.  This should have been done elsewhere!

     vec<Bool> to_delete( closures.size( ), False );
     vec<basevector> closures2(closures);
     for ( int i = 0; i < closures.isize( ); i++ )
          if ( closures[i].isize( ) < K ) to_delete[i] = True;
     EraseIf( closures2, to_delete );

     clop.resize( closures2.size( ) );
     vec< pair<kmer<K>,int> > X( closures2.size( ) );
     #pragma omp parallel for
     for ( int i = 0; i < closures2.isize( ); i++ )
     {    const basevector& c = closures2[i];
          kmer<K> x;
          x.SetToSubOf( c, 0 );
          X[i] = make_pair( x, i );    }
     ParallelSort(X);
     cout << Date( ) << ": mapping closures back to assembly" << endl;
     #pragma omp parallel for schedule( dynamic, 10000 )
     for ( int g = 0; g < hb.E( ); g++ )
     {    const basevector& E = hb.O(g);
          kmer<K> x;
          for ( int r = 0; r < hb.Kmers(g); r++ )
          {    x.SetToSubOf( E, r );
               int64_t low = LowerBound1( X, x ), high = UpperBound1( X, x );
               for ( int64_t m = low; m < high; m++ )
               {    int c = X[m].second;
                    const basevector& C = closures2[c];
                    int e = g, epos = r;

                    // closure c starts at position epos on on edge e
                    // trace it through the assembly

                    int offset = epos;
                    vec<int> p = {e};
                    int cpos = 0;
                    while(1)
                    {    cpos += hb.Bases(e) - epos;
                         if ( cpos >= C.isize( ) ) break;
                         int v = hb.ToRight(e);
                         Bool extended = False;
                         for ( int j = 0; j < (int) hb.From(v).size( ); j++ )
                         {    int f = hb.IFrom( v, j );
                              const basevector& F = hb.O(f);
                              if ( F[K-1] == C[cpos] )
                              {    extended = True;
                                   p.push_back(f);
                                   e = f;
                                   epos = K-1;
                                   break;    }    }
                         if ( !extended ) break;    }
                    ForceAssertGe( cpos, C.isize( ) );
                    clop[c].setOffset(offset);
                    for ( int j = 0; j < p.isize( ); j++ )
                         clop[c].push_back( p[j] );    }    }    }    }

template void MapClosures<40>( const HyperBasevectorX&, const vec<basevector>&,
     ReadPathVec& );
template void MapClosures<48>( const HyperBasevectorX&, const vec<basevector>&,
     ReadPathVec& );
template void MapClosures<60>( const HyperBasevectorX&, const vec<basevector>&,
     ReadPathVec& );

int64_t EstimateGEMCount( const vec<int64_t> & bci, const int64_t total_diversity = 0)
{
     // if this is zero then we turn off the estimation
     if (total_diversity == 0)
          return 0;
     vec <int64_t> rpb ( bci.size() - 1 ); // reads per barcode
     #pragma omp parallel for
     for ( int b = 1; b < bci.isize( ) - 1; b++ )
     {    int64_t start = bci[b], stop = bci[b+1];
          ForceAssertGe( stop, start );
          rpb[b] = stop-start+1; 
     }

     ParallelSort(rpb);
     
     // Define good barcodes for GEM estimation
     // Use N99 AND
     // Use number of reads per barcode >= 4
     const double perc_threshold = 0.99;
     vec <int> cum_sum( rpb.size(), 0 );

     cum_sum[0] = rpb[0];
     for ( int i = 1; i < rpb.isize(); i++ )
          cum_sum[i] = cum_sum[i-1]+rpb[i];
     int64_t per_5_index = -1, nreads = cum_sum[cum_sum.size()-1];
     for ( per_5_index = 0; per_5_index < cum_sum.isize(); per_5_index++ ) {
          if ( cum_sum[per_5_index] >= (1-perc_threshold)*nreads && rpb[per_5_index] >= 4 )
               break;
     }
     int64_t bcs = rpb.size() - per_5_index;
     if ( bcs == 0 )
          FatalErr ("No barcodes for GEM estimation.");
     
     double mean_gems_per_bc = 0.0;
     // In this case we cannot estimate the number of GEMs
     if ( bcs >= total_diversity )
          return 0;
     else {
          double p_occupied = double(bcs)/double(total_diversity);
          mean_gems_per_bc = -std::log( 1 - p_occupied );
     }
     int64_t num_gems = mean_gems_per_bc*total_diversity;
     return num_gems;
}

void SanityCheckBarcodeCounts( const vec<int64_t>& bci )
{    if ( bci.size( ) < 2 ) return;
     // if you change big, then change it in the tenkit alert too!
     const int big = 50000;
     int64_t big_total = 0;
     vec <int> rpb;
     for ( int j = 2; j < bci.isize( ); j++ )
     {    int64_t n = bci[j] - bci[j-1];
          if ( n >= big ) big_total += n;
          if (n > 0) rpb.push_back(int(n));  }
     // First issue any rpb alerts
     int rpb_n50 = N50 ( rpb );
     StatLogger::log( "rpb_N50", rpb_n50, "N50 reads per barcode", true );
     StatLogger::issue_alert( "rpb_N50", rpb_n50 );
     // next do the big barcode alert
     int64_t total = bci.back( ) - bci[1];
     if ( total > 0 )
     {    double big_bc_perc = double(big_total) / double(total);
          StatLogger::log( "big_bc_perc", big_bc_perc, "Pct big barcodes");
          StatLogger::issue_alert( "big_bc_perc", big_bc_perc ); }    }

template <class T>
void MakeDots( T& done, T& ndots, const T total )
{    if ( done % ( (total+99) / 100 ) == 0 )
     {    if ( ndots < 100 )
          {    cout << ".";
               ndots++;
               if ( ndots > 0 && ndots % 50 == 0 ) cout << "\n";
               else if ( ndots > 0 && ndots % 10 == 0 ) cout << " ";
               flush(cout);    }    }
     if ( done == total - 1 )
     {    while ( ndots < 100 )
          {    cout << ".";
               ndots++;
               if ( ndots > 0 && ndots % 50 == 0 ) cout << "\n";
               else if ( ndots > 0 && ndots % 10 == 0 ) cout << " ";
               flush(cout);    }    }
     done++;    }

template void MakeDots( int& done, int& ndots, const int total );
template void MakeDots( int64_t& done, int64_t& ndots, const int64_t total );

// Note grossly inefficient conversion below.
// Note that we might want to factor this through Munch.

void SuperToSeqGraph( const HyperBasevectorX& hb, const digraphE<vec<int>>& D,
     HyperBasevectorX& hbd )
{    vec<basevector> edges( D.E( ) );
     vec<int> to_left, to_right;
     D.ToLeft(to_left), D.ToRight(to_right);
     vec<int> seqverts;

     // Convert most edges.

     #pragma omp parallel for
     for ( int d = 0; d < D.E( ); d++ )
     {    if ( D.O(d)[0] >= 0 ) edges[d] = hb.Cat( D.O(d) );
          else if ( IsSequence( D.O(d) ) ) 
          {
               #pragma omp critical
               {    seqverts.push_back( to_left[d] );    }    }    }
     UniqueSort(seqverts);

     // Go through the vertices to left of sequence gaps.
     // Not necessarily symmetric.

     int nprobs = 0;
     for ( int i = 0; i < seqverts.isize( ); i++ )
     {    int v = seqverts[i];
          int ltrim, rtrim;

          // Test for problem.

          int d = D.IFrom( v, 0 );
          int w = to_right[d];
          GapToSeq( D.O(d), ltrim, rtrim, edges[d] );    
          if ( D.To(v).empty( ) || D.From(w).empty( ) )
          {    cout << "\nPROBLEM:\n";
               PRINT3( d, D.To(v).size( ), D.From(w).size( ) );
               Scram(0);    }
          int d1 = D.ITo(v,0), d2 = D.IFrom(w,0);
          Bool problem = False;
          if ( edges[d1].isize( ) - ltrim < hb.K( ) )
          {    // cout << "\nPROBLEM:\n";
               // PRINT4( d, d1, ltrim, edges[d1].size( ) );
               problem = True;    }
          if ( edges[d2].isize( ) - rtrim < hb.K( ) )
          {    // cout << "\nPROBLEM:\n";
               // PRINT4( d, d2, rtrim, edges[d2].size( ) );
               problem = True;    }
          if ( d1 == d2 && edges[d1].isize( ) - ltrim - rtrim < hb.K( ) )
               problem = True;
          if (problem)
          {    for ( int j = 0; j < D.From(v).isize( ); j++ )
               {    int d = D.IFrom( v, j );
                    edges[d].clear( );    }
               nprobs++;
               continue;    }

          // Proceed.

          for ( int j = 0; j < D.From(v).isize( ); j++ )
          {    int d = D.IFrom( v, j );
               int w = to_right[d];
               if ( !D.To(v).solo( ) ) cout << "Problem 1 at edge " << d << endl;
               ForceAssert( D.To(v).solo( ) );
               if ( !D.From(w).solo( ) ) cout << "Problem 2 at edge " << d << endl;
               ForceAssert( D.From(w).solo( ) );
               GapToSeq( D.O(d), ltrim, rtrim, edges[d] );    
               if ( j == 0 )
               {    int d1 = D.ITo(v,0), d2 = D.IFrom(w,0);
                    if ( ltrim > 0 ) edges[d1].resize( edges[d1].isize( ) - ltrim );
                    if ( rtrim > 0 ) 
                    {    if ( edges[d2].isize( ) - rtrim < hb.K( ) )
                         {    cout << "\nPROBLEM:\n";
                              PRINT4( d, d2, rtrim, edges[d2].size( ) );    }
                         edges[d2].SetToSubOf( edges[d2], rtrim, 
                              edges[d2].isize( ) - rtrim );    }    }    }    }
#ifndef CS
     cout << Date( ) << ": problems converting " << nprobs << " of "
          << seqverts.size( ) << " gaps" << endl;
#endif
     HyperBasevector H( hb.K( ), D, edges );
     // NOTE EXPENSIVE CONVERSION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     hbd = HyperBasevectorX(H);    }
