// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "MainTools.h"
#include "ParallelVecUtilities.h"
#include "paths/HyperBasevector.h"
#include "paths/long/large/GapToyTools.h"
#include "10X/astats/GenomeAlign.h"

void AlignToGenomeCore( const vecbasevector& tigs, const vecbasevector& genome, 
     vec< vec< pair<int,int> > >& hits, const int K2, const int max_gmult )
{
     // Heuristics.

     const int K = 60;

     // Computational performance heuristics.

     const int64_t batches = 100;

     // Go through four passes.

     int nobj = tigs.size( ), ngenome = genome.size( );
     hits.clear( );
     hits.resize(nobj);
     for ( int pi = 0; pi < 4; pi++ )
     {    
          // Build lookup table.

          double clock2 = WallClockTime( );
          vec< triple<kmer<K>,int,int> > kmers_plus;
          vec<int64_t> starts;
          starts.reserve( ngenome + nobj );
          starts.push_back(0);
          {    vec<int> counts( ngenome + nobj, 0 );
               #pragma omp parallel for
               for ( int i = 0; i < ngenome; i++ )
               {    const basevector& u = genome[i];
                    counts[i] = 0;
                    for ( int j = 0; j <= u.isize( ) - K2; j++ )
                         if ( u[j] == pi ) counts[i]++;    }
               #pragma omp parallel for
               for ( int i = 0; i < nobj; i++ )
               {    const basevector& u = tigs[i];
                    counts[ genome.size( ) + i ] = 0;
                    for ( int j = 0; j <= u.isize( ) - K2; j++ )
                         if ( u[j] == pi ) counts[ ngenome + i ]++;    }
               for ( int i = 0; i < ngenome + nobj; i++ )
                    starts.push_back( starts.back( ) + counts[i] );    }
          kmers_plus.resize( starts.back( ) );
          #pragma omp parallel for
          for ( int i = 0; i < ngenome; i++ )
          {    const basevector& u = genome[i];
               int count = 0;
               kmer<K> x;
               for ( int j = 0; j <= u.isize( ) - K2; j++ )
               {    if ( u[j] == pi )
                    {    int64_t r = starts[i] + count;
                         x.SetToSubOf( u, j ); 
                         kmers_plus[r].first = x;
                         kmers_plus[r].second = i, kmers_plus[r].third = j;    
                         count++;    }    }    }
          #pragma omp parallel for
          for ( int i = 0; i < nobj; i++ )
          {    const basevector& u = tigs[i];
               int count = 0;
               kmer<K> x;
               for ( int j = 0; j <= u.isize( ) - K2; j++ )
               {    if ( u[j] == pi )
                    {    int64_t r = starts[ngenome+i] + count;
                         x.SetToSubOf( u, j ); 
                         kmers_plus[r].first = x;
                         kmers_plus[r].second = ngenome + i, kmers_plus[r].third = j;
                         count++;    }    }    }
          ParallelSort(kmers_plus);

          // Traverse the kmers.

          vec<int64_t> bstart(batches+1);
          for ( int64_t i = 0; i <= batches; i++ )
               bstart[i] = ( (int64_t) kmers_plus.size( ) * i ) / batches;
          #pragma omp parallel for schedule(dynamic, 1)
          for ( int64_t i = 1; i < batches; i++ )
          {    int64_t& s = bstart[i];
               while( s > 0 && kmers_plus[s].first == kmers_plus[s-1].first )
               {    s--;    }    }
          #pragma omp parallel for schedule(dynamic, 1)
          for ( int64_t bi = 0; bi < batches; bi++ )
          for ( int64_t i = bstart[bi]; i < bstart[bi+1]; i++ )
          {    int64_t j, k;
               for ( j = i + 1; j < bstart[bi+1]; j++ )
                    if ( kmers_plus[j].first != kmers_plus[i].first ) break;
               for ( k = i; k < j; k++ )
                    if ( kmers_plus[k].second >= ngenome ) break;
               if ( k > i && j - k >= 1 && k - i <= max_gmult )
               {    for ( int64_t r = k; r < j; r++ )
                    {    int e = kmers_plus[r].second - ngenome;
                         int estart = kmers_plus[r].third; 
                         for ( int64_t m = i; m < k; m++ )
                         {    int g = kmers_plus[m].second;
                              int gstart = kmers_plus[m].third;
                              int offset = gstart - estart;
                              if ( !Member( hits[e], make_pair( g, offset ) ) )
                              {    Bool mismatch = False;
                                   for ( int l = K2 - 1; l >= K; l-- )
                                   {    if ( tigs[e][estart+l] 
                                             != genome[g][gstart+l] )
                                        {    mismatch = True;
                                             break;    }    }
                                   if ( !mismatch ) 
                                   {
                                        #pragma omp critical
                                        {    hits[e].push( g, offset );    
                                                  }    }    }    }    }    }
               i = j - 1;    }    }    }

// AlignToGenomeX.  Find alignments of assembly edges to the genome.  Currently this
// maps edges to pairs (g,p) consisting of a genome contig g and an inferred start
// position p of the edge on g.
//
// First we find 60-mer matches between the assembly and the genome.  We require
// that the K2-mer occur at most max_mult times in the genome.  

void AlignToGenomeX( const HyperBasevectorX& hb, const vec<int>& inv,
     const vecbasevector& genome, vec< vec< pair<int,int> > >& hits,
     const int K2, const int max_gmult, const Bool adjudicate, const Bool bubble )
{    double clock = WallClockTime( );

     // Heuristics.

     const int K = 60;
     const int max_gap = 20000;

     // Initial alignment.

     AlignToGenomeCore( hb.Edges( ), genome, hits, K2, max_gmult );

     // Copy alignments from one bubble branch to another, if the other is naked.
     // Copy alignments into the bubble if it makes sense.

     if (bubble)
     {
     for ( int v = 0; v < hb.N( ); v++ )
     {    if ( hb.To(v).size( ) != 1 || hb.From(v).size( ) != 2 ) continue;
          if ( hb.From(v)[0] != hb.From(v)[1] ) continue;
          int w = hb.From(v)[0];
          if ( hb.To(w).size( ) != 2 || hb.From(w).size( ) != 1 ) continue;
          int e1 = hb.IFrom( v, 0 ), e2 = hb.IFrom( v, 1 );
          int d = hb.ITo( v, 0 );
          int f = hb.IFrom( w, 0 );
          if ( hits[e1].empty( ) && hits[e2].empty( ) && hits[d].solo( )
               && hits[f].solo( ) && hits[d][0].first == hits[f][0].first )
          {    int pos = hits[d][0].second + hb.Bases(d) - ( hb.K( ) - 1 );
               int Pos = hits[f][0].second + ( hb.K( ) - 1 );
               if ( Pos - pos >= hb.K( ) )
               {    hits[e1].push( hits[d][0].first, pos );
                    hits[e2].push( hits[d][0].first, pos );    }    }
          else
          {    if ( hits[e1].empty( ) ) swap( e1, e2 );
               if ( hits[e1].empty( ) || hits[e2].nonempty( ) ) continue;
               hits[e2] = hits[e1];    }    }
     }

     // Adjudicate some alignments.

     if (adjudicate)
     {
     for ( int e = 0; e < hb.E( ); e++ )
     {    int re = inv[e];
          if ( hits[e].size( ) + hits[re].size( ) < 2 || e == re ) continue;
          vec<int> errors( hits[e].size( ) + hits[re].size( ), 0 );
          for ( int pass = 1; pass <= 2; pass++ )
          {    int f = ( pass == 1 ? e : re );
               int start = ( pass == 1 ? 0 : hits[e].size( ) );
               for ( int j = 0; j < hits[f].isize( ); j++ )
               {    int g = hits[f][j].first;
                    int gpos = hits[f][j].second;
                    for ( int l = 0; l < hb.Bases(f); l++ )
                    {    if ( gpos+l < 0 || gpos+l >= genome[g].isize( ) ) 
                              errors[start+j]++;
                         else if ( hb.O(f)[l] != genome[g][gpos+l] ) 
                         {    errors[start+j]++;    }    }    }
               vec<int> ids( errors.size( ), vec<int>::IDENTITY );
               SortSync( errors, ids );
               if ( errors[0] == 0 && errors[1] > 0 )
               {    int id = ids[0];
                    if ( id < hits[e].isize( ) )
                    {    hits[e] = { hits[e][id] };
                         hits[re].clear( );    }
                    else
                    {    hits[re] = { hits[re][ id - hits[e].isize( ) ] };
                         hits[e].clear( );    }    }    }    }
     }
     
     LogTime( clock, "aligning to genome" );    }

template<int K> void GenomeAlign( const HyperBasevectorX& hb, const vec<int>& inv,
     const vecbasevector& genome, 
     MasterVec< SerfVec<triple<int,int,int> > >& alignsb, const int max_gmult,
     const Bool adjudicate, const Bool bubble )
{
     vec< vec< pair<int,int> > > aligns;
     cout << Date( ) << ": aligning to genome, mem = "
          << MemUsageGBString() << endl;
     AlignToGenomeX( hb, inv, genome, aligns, K, max_gmult, adjudicate, bubble );

     // Create mastervec version, with nearby alignments merged,
     // entries presented as (chr,start,stop).

     const int MAX_JUMP = 100;
     alignsb.clear( );
     alignsb.resize( hb.E( ) );
     for ( int e = 0; e < hb.E( ); e++ )
     {    for ( int i = 0; i < aligns[e].isize( ); i++ )
          {    int j, g = aligns[e][i].first;
               for ( j = i + 1; j < aligns[e].isize( ); j++ )
               {    if ( aligns[e][j].first != g ) break;
                    if ( aligns[e][j].second - aligns[e][j-1].second > MAX_JUMP )
                         break;    }
               alignsb[e].push_back( make_triple( g, aligns[e][i].second,
                    aligns[e][j-1].second + hb.Bases(e) ) );
               i = j - 1;    }    }    }


template void GenomeAlign<60>( const HyperBasevectorX&, const vec<int>&,
     const vecbasevector&, MasterVec< SerfVec<triple<int,int,int> > >&, const int,
     const Bool, const Bool );
template void GenomeAlign<80>( const HyperBasevectorX&, const vec<int>&,
     const vecbasevector&, MasterVec< SerfVec<triple<int,int,int> > >&, const int,
     const Bool, const Bool );
