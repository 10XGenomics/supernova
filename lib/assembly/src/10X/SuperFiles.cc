// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

// Generate full set of super-assembly files.  Starts with only the super-assembly
// data structures D and dinv.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"
#include "STLExtensions.h"
#include "feudal/SmallVec.h"
#include "math/HoInterval.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "paths/long/large/Lines.h"
#include "10X/DfTools.h"
#include "10X/Gap.h"
#include "10X/Heuristics.h"
#include "10X/LineOO.h"
#include "10X/PlaceReads.h"
#include "10X/Super.h"
#include "10X/SuperFiles.h"
#include "10X/astats/RefAlign.h"
#include "10X/astats/View.h"

void GenerateDpathsCounts(ReadPathVec const& dpaths, VecULongVec const& dpaths_index, 
          vec<int> const& dinv, vec<uint16_t>& dpaths_counts)
{
     dpaths_counts.resize( dpaths_index.size() );
#pragma omp parallel for
     for ( int i = 0; i < dpaths_counts.isize(); ++i ) {
          auto const& read_list = dpaths_index[i];
          auto const& inv_read_list  = dpaths_index[dinv[i]];
          vec<unsigned long> combined;
          combined.reserve( read_list.size() + inv_read_list.size() );

          copy( read_list.begin(), read_list.end(), 
                    std::back_inserter(combined) );

          if ( dinv[i] != i ) {
               copy( inv_read_list.begin(), inv_read_list.end(), 
                         std::back_inserter(combined) );
          }

          UniqueSort(combined);
          if ( combined.size() > ( 1 << 16 ) - 1 ) {
               dpaths_counts[i]=~0;
          } else {
               dpaths_counts[i] = combined.size();
          }
     }
}

void SuperFiles( const HyperBasevectorX& hb, const vec<int>& inv,
     const vec<int64_t>& bci, const vec<int32_t>& bc,
     const ReadPathVecX& paths, const vec<Bool>& dup,
     const digraphE<vec<int>>& D, const vec<int>& dinv,
     ReadPathVec& dpaths,              // recomputed!
     vec<vec<vec<vec<int>>>>& dlines,  // recomputed!
     vec<double>& COV,                 // recomputed!
     const VecIntVec& ebcx, const vec< triple<int,int,int> >& qept,
     const vecbasevector& genome,
     MasterVec< SerfVec<triple<int,int,int> > >& alignsb,
     const String& DIR, const String& WRITE_SUB )
{    
     // Place reads and find lines.

     cout << Date( ) << ": start final read placement" << endl;
     const Bool single = False;
     PlaceReads( hb, paths, dup, D, dpaths, True, single );
     FindLines( D, dinv, dlines, MAX_CELL_PATHS, MAX_CELL_DEPTH );
     PlaceReadsSmart( hb, paths, dup, D, dinv, dpaths, dlines, bci, False );
     cout << Date( ) << ": making index" << endl;
     vec<int> mult;
     cout << Date( ) << ": writing scaffolded assembly" << endl;
     String OUTDIR = DIR + "/final" + WRITE_SUB;
     Mkdir777(OUTDIR);
     auto FN = [&]( const String& x ){ return OUTDIR + "/" + x; };

     // First write some key output files, so that in case something goes wrong
     // below, we can try to trace the problem with Probe.

     #pragma omp parallel sections
     {
          #pragma omp section
          {    BinaryWriter::writeFile( FN( "/a.sup" ), D );    }
          #pragma omp section
          {    BinaryWriter::writeFile( FN( "/a.sup.inv" ), dinv );    }
          #pragma omp section
          {    BinaryWriter::writeFile( FN( "/a.sup.lines" ), dlines );    }    }
     if ( genome.size( ) > 0 )
     {    MasterVec<SerfVec<refalign>> galigns;
          RefAlign( genome, hb.K( ), hb.Edges( ), inv, D, dinv, dlines, alignsb, 
               galigns, False, vec<int>( ) );
          galigns.WriteAll( FN( "/a.sup.galigns" ) );    }

     // Now write the rest of the files.

     #pragma omp parallel sections
     {
          #pragma omp section
          {    dpaths.WriteAll( FN( "/a.dpaths" ) );    }
          #pragma omp section
          {    VecULongVec dpaths_index;
               invert( dpaths, dpaths_index, D.E( ) );
               dpaths_index.WriteAll( FN( "/a.dpaths.index" ) );    
               vec<uint16_t> dpaths_counts;
               GenerateDpathsCounts(dpaths, dpaths_index, dinv, dpaths_counts);
               BinaryWriter::writeFile( FN( "/a.dpaths.counts"), dpaths_counts  );
          }
          #pragma omp section
          {    HyperBasevectorX hbd;
               SuperToSeqGraph( hb, D, hbd );
               BinaryWriter::writeFile( FN( "/a.sup.hbx"), hbd );
               hbd.Edges( ).WriteAll( FN( "/a.sup.fastb" ) );    }
          #pragma omp section
          {
               // Find barcodes on supergraph edges.

               ComputeMult( hb, D, mult );
               vec<vec<int>> debc( D.E( ) );
               #pragma omp parallel for
               for ( int d = 0; d < D.E( ); d++ )
               {    vec<int> bc;
                    const vec<int>& x = D.O(d);
                    if ( x[0] < 0 ) continue;
                    for ( int j = 0; j < x.isize( ); j++ )
                         if ( mult[ x[j] ] == 1 ) bc.append( ebcx[ x[j] ] );
                    UniqueSort(bc);
                    debc[d] = bc;    }
               VecIntVec debcx( debc.size( ) );
               for ( int64_t i = 0; i < debc.jsize( ); i++ )
               {    debcx[i].resize( debc[i].size( ) );
                    for ( int j = 0; j < debc[i].isize( ); j++ )
                         debcx[i][j] = debc[i][j];    }
               debcx.WriteAll( FN( "/a.sup.ebc" ) );    }
          #pragma omp section
          {
               // Write line ancillary files.

               vec<int> llens;
               GetLineLengths( hb, D, dlines, llens );
               BinaryWriter::writeFile( FN( "/a.sup.llens" ), llens );
               vec< vec< pair<int,int> > > lhood;
               LineProx( hb, inv, ebcx, D, dinv, dlines, qept, lhood );
               BinaryWriter::writeFile( FN( "/a.sup.lhood" ), lhood );
               vec<int> kmers( hb.E( ) );
               #pragma omp parallel for
               for ( int e = 0; e < hb.E( ); e++ )
                    kmers[e] = hb.Kmers(e);
               vec<vec<pair<int,int>>> lbp;
               IntIndex dpaths_index( dpaths, D.E( ) );
               BarcodePos( bc, hb, D, dinv, dpaths, dpaths_index, dlines, lbp, 0 );
               BinaryWriter::writeFile( FN( "/a.sup.lbp" ), lbp );
               MasterVec<SerfVec<pair<int,int>>> lbpx;
               for ( auto x : lbp )
               {    SerfVec<pair<int,int>> y( x.size( ) );
                    for ( int j = 0; j < x.isize( ); j++ )
                         y[j] = x[j];
                    lbpx.push_back(y);    }
               lbpx.WriteAll( FN( "/a.sup.lbpx" ) );
               LineCN( kmers, lbpx, D, dlines, llens, COV );
               BinaryWriter::writeFile( FN( "/a.sup.lcov" ), COV );
               if ( alignsb.size( ) > 0 )
               {    vec<vec< pair<int,int> >> linelocs( kmers.size( ) );
                    for ( int i = 0; i < dlines.isize( ); i++ )
                    for ( int j = 0; j < dlines[i].isize( ); j++ )
                    for ( int k = 0; k < dlines[i][j].isize( ); k++ )
                    for ( int l = 0; l < dlines[i][j][k].isize( ); l++ )
                    {    int d = dlines[i][j][k][l];
                         if ( D.O(d)[0] < 0 ) continue;
                         for ( int m = 0; m < D.O(d).isize( ); m++ )
                         {    int e = D.O(d)[m];
                              linelocs[e].push( i, j );    }    }
                    MasterVec< SerfVec< quad<int,Bool,ho_interval,ho_interval> > > 
                         view( dlines.size( ) );
                    #pragma omp parallel for
                    for ( int l = 0; l < dlines.isize( ); l++ )
                    {    View( l, hb.K( ), kmers, inv, D, dlines, 
                              linelocs, alignsb, view[l] );    }
                    view.WriteAll( FN( "a.sup.view" ) );    }    }    }    }
