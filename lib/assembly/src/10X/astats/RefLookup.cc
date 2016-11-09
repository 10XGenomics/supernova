// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"
#include "FetchReads.h"
#include "paths/HyperBasevector.h"
#include "paths/RemodelGapTools.h"
#include "random/Shuffle.h"
#include "10X/astats/RefLookup.h"

void FetchFinished( const String& SAMPLE, vecbasevector& G )
{    if ( SAMPLE == "NA12878" )
          FetchReads( G, 0, "/mnt/opt/meowmix_git/assembly/refs/fos100.fasta" );
     else if ( SAMPLE == "HGP" )
     {    String genome = "/mnt/opt/meowmix_git/assembly/refs/mrx/genome.fastb";
          int n = MastervecFileObjectCount(genome);
          vec<int> x;
          Shuffle( n, x );
          x.resize(400);
          G.Read( genome, x );    }    }
