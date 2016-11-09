// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

// Inverter.  Effectively, transpose a vec<vec<int>>.
// However, applied to a ReadPathVec.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>

#include "MainTools.h"
#include "ParallelVecUtilities.h"
#include "paths/long/ReadPath.h"
#include "10X/IntIndex.h"

void IntIndex::Initialize( 
     const ReadPathVec& paths, const int NE, const Bool verbose )
{
     // Clear.

     core_.clear( );
     core_big_.clear( );
     index_.clear( );

     // Decide if big.

     big_ = ( paths.size( ) > 4294967296l );

     // Count elements.

     if (verbose) cout << Date( ) << ": creating an index" << endl;
     double clock = WallClockTime( );
     int64_t N = paths.size( );
     vec<int> count( NE, 0 );
     {    const int64_t batches = 10;
          vec<vec<int>> counti( batches, vec<int>( NE, 0 ) );
          #pragma omp parallel for schedule( dynamic, 1 )
          for ( int b = 0; b < batches; b++ )
          for ( int64_t id = (b*N)/batches; id < ((b+1)*N)/batches; id++ )
          {    const ReadPath& p = paths[id];
               for ( auto e : p ) counti[b][e]++;    }
          #pragma omp parallel for
          for ( int e = 0; e < NE; e++ )
          for ( int t = 0; t < batches; t++ )
               count[e] += counti[t][e];    }

     // Set up.

     const int64_t sbatch = 20000;
     vec<vec<int>> start;

     // Run in parallel.

     #pragma omp parallel sections
     {
          // Find estarts.

          #pragma omp section
          {    index_.resize( NE + 1 );
               index_[0] = 0;
               for ( int e = 0; e < NE; e++ )
                    index_[e+1] = index_[e] + count[e];    }
     
          // Allocate space.
     
          #pragma omp section
          {    int64_t total = BigSum(count);
               if ( !big_ ) core_.resize(total);
               else core_big_.resize(total);    }
     
          // Find starts.
     
          #pragma omp section
          {    int sbatches = 0;
               for ( int64_t bi = 0; bi < (int64_t) paths.size( ); bi += sbatch )
                    sbatches++;
               start.resize(sbatches);
               #pragma omp parallel for
               for ( int64_t bi = 0; bi < (int64_t) paths.size( ); bi += sbatch )
               {    int bid = bi/sbatch;
                    int nb = Min( sbatch, (int64_t) paths.size( ) - bi );
                    start[bid].resize( nb+1, 0 );
                    for ( int j = 0; j < nb; j++ )
                    {    start[bid][j+1] = start[bid][j] 
                              + paths[bi+j].size( );    }    }    }    }

     // Shovel data.

     vec<int> pos( NE, 0 );
     for ( int64_t bi = 0; bi < (int64_t) paths.size( ); bi += sbatch )
     {    int bid = bi/sbatch, nb = Min( sbatch, (int64_t) paths.size( ) - bi );
          vec< pair<int,int64_t> > D( start[bid].back( ) );
          #pragma omp parallel for
          for ( int j = 0; j < nb; j++ )
          {    const ReadPath& p = paths[bi+j];
               for ( int l = 0; l < (int) p.size( ); l++ )
               {    D[ start[bid][j] + l ] = make_pair( p[l], bi+j );    }    }
          ParallelSort(D);    
          int nd = 2 * omp_get_max_threads( );
          vec<int64_t> dstarts(nd+1);
          for ( int64_t j = 0; j <= nd; j++ )
               dstarts[j] = ( j * D.jsize( ) ) / nd;
          for ( int j = 0; j < nd; j++ )
          {    while( dstarts[j] > 0 
                    && D[ dstarts[j] ].first == D[ dstarts[j] - 1 ].first )
               {     dstarts[j]--;    }    }

          // Move the data.

          if ( !big_ )
          {
               #pragma omp parallel for schedule( dynamic, 1 )
               for ( int b = 0; b < nd; b++ )
               {    for ( int64_t m = dstarts[b]; m < dstarts[b+1]; m++ )
                    {    int e = D[m].first;
                         int64_t id = D[m].second;
                         core_[ index_[e] + pos[e] ] = id;
                         pos[e]++;    }    }    }
          else
          {    
               #pragma omp parallel for schedule( dynamic, 1 )
               for ( int b = 0; b < nd; b++ )
               {    for ( int64_t m = dstarts[b]; m < dstarts[b+1]; m++ )
                    {    int e = D[m].first;
                         int64_t id = D[m].second;
                         core_big_[ index_[e] + pos[e] ] = id;
                         pos[e]++;    }    }    }    }

     if (verbose)
     {    cout << Date( ) << ": " << TimeSince(clock) << " used indexing" 
               << endl;    }    }

IntIndex::IntIndex( const ReadPathVec& paths, const int NE, const Bool verbose )
{    Initialize( paths, NE, verbose );    }
