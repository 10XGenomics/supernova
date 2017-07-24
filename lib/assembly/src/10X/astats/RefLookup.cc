// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"
#include "FetchReads.h"
#include "paths/HyperBasevector.h"
#include "paths/RemodelGapTools.h"
#include "random/Shuffle.h"
#include "10X/astats/RefLookup.h"


String SampleGenome( const String& SAMPLE )
{
     if ( SAMPLE == "NA12878" ) 
          return "/mnt/opt/meowmix_git/assembly/refs/fos100.fasta";
     if ( SAMPLE == "NA12891" )
          return "/mnt/opt/meowmix_git/assembly/refs/NA12891/child.fasta";
     if ( SAMPLE == "NA12892" )
          return "/mnt/opt/meowmix_git/assembly/refs/NA12892/child.fasta";
     else if ( SAMPLE == "HGP" )
          return "/mnt/opt/meowmix_git/assembly/refs/mrx/genome.fastb";
     else if ( SAMPLE == "CHM" )
          return "/mnt/opt/meowmix_git/assembly/refs/chm/mhc_chm1_chm13.fastb";

     return "";
}


void FetchFinished( const String& SAMPLE, vecbasevector& G, int N, const int seed )
{    
     String genome = SampleGenome( SAMPLE );
     
     if ( SAMPLE == "NA12878" )
          FetchReads( G, 0, genome );
     else if ( SAMPLE == "NA12891" )
          FetchReads( G, 0, genome );
     else if ( SAMPLE == "NA12892" )
          FetchReads( G, 0, genome );
     else if ( SAMPLE == "HGP" )
     {
          int n = MastervecFileObjectCount(genome);
          vec<int> x;
          if ( seed != 0 ) Shuffle( n, x, seed );
          else x = vec<int>( n, vec<int>::IDENTITY );
          if ( N == -1 ) N = n;
          else if ( N == 0 ) N = 400;
          ForceAssertLe( N, n );
          x.resize(N);
          G.Read( genome, x );    }
     else if ( SAMPLE == "CHM" )
     {
          int n = MastervecFileObjectCount(genome);
          vec<int> x( n, vec<int>::IDENTITY );
          G.Read( genome, x );    }         
}
