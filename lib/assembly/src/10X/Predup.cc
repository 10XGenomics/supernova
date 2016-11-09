// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

// Pre-deduplication.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "CoreTools.h"
#include "VecUtilities.h"
#include "feudal/PQVec.h"
#include "kmers/KmerRecord.h"
#include "10X/Predup.h"

void Predup( vecbasevector& bases, VecPQVec& quals, vec<int64_t>& bci )
{
     // First mark dups.

     cout << Date( ) << ": start predup" << endl;
     vec<int64_t> sizes( bci.isize( ) - 1 );
     sizes[0] = bci[1];
     vec<Bool> to_delete( bases.size( ), False );
     cout << Date( ) << ": marking duplicates" << endl;
     #pragma omp parallel for
     for ( int bi = 1; bi < bci.isize( ) - 1; bi++ )
     {    int64_t start = bci[bi], stop = bci[bi+1];

          // Form first 40 bases of first read + first 40 bases of second read.
          // If another pair has the same value, declare them duplicates of each
          // other.

          const int K = 40;
          vec< triple< kmer<K>, kmer<K>, int64_t > > sig( (stop-start)/2 );
          for ( int64_t pid = start/2; pid < stop/2; pid++ )
          {    sig[pid - start/2].first.SetToSubOf( bases[2*pid], 0 );
               sig[pid - start/2].second.SetToSubOf( bases[2*pid+1], 0 );
               sig[pid - start/2].third = pid;    }
          Sort(sig);
          vec<int> qsum;
          qualvector q;
          for ( int64_t j = 0; j < sig.jsize( ); j++ )
          {    int64_t k;
               for ( k = j + 1; k < sig.jsize( ); k++ )
               {    if ( sig[k].first != sig[j].first ) break;
                    if ( sig[k].second != sig[j].second ) break;    }
               if ( k - j > 1 )
               {    qsum.resize_and_set( k - j, 0 );
                    for ( int64_t l = j; l < k; l++ )
                    {    int64_t pid = sig[l].third;
                         quals[2*pid].unpack(&q);
                         for ( int m = 0; m < (int) q.size( ); m++ )
                              qsum[l-j] += q[m];
                         quals[2*pid+1].unpack(&q);
                         for ( int m = 0; m < (int) q.size( ); m++ )
                              qsum[l-j] += q[m];    }
                    vec<int> ids( k - j, vec<int>::IDENTITY );
                    ReverseSortSync( qsum, ids );
                    for ( int64_t l = j; l < k; l++ )
                    {    if ( l - j == ids[0] ) continue;
                         int64_t pid = sig[l].third;
                         to_delete[ 2*pid ] = True;
                         to_delete[ 2*pid + 1 ] = True;    }    }
               sizes[bi] += 2;
               j = k - 1;    }    }

     // Now delete dups.

     cout << Date( ) << ": deleting dups" << endl;
     bases.EraseIf(to_delete), quals.EraseIf(to_delete);
     for ( int bi = 1; bi < bci.isize( ) - 1; bi++ )
          bci[bi+1] = bci[bi] + sizes[bi];    }
