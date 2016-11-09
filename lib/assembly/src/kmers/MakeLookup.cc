///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "CoreTools.h"
#include "ParallelVecUtilities.h"
#include "kmers/KmerRecord.h"
#include "kmers/MakeLookup.h"

template<int K> void MakeKmerLookup0Pre( const vecbasevector& unibases,
     const String& prefix, vec< triple<kmer<K>,int,int> >& kmers_plus )
{    vec<int64_t> starts;
     starts.push_back(0);
     {    vec<int> counts( unibases.size( ), 0 );
          #pragma omp parallel for
          for ( size_t i = 0; i < unibases.size( ); i++ )
          {    const basevector& u = unibases[i];
               String s = u.ToString( );
               counts[i] = 0;
               for ( int j = 0; j <= u.isize( ) - K; j++ )
                    if ( s.Contains( prefix, j ) ) counts[i]++;    }
          for ( size_t i = 0; i < unibases.size( ); i++ )
               starts.push_back( starts.back( ) + counts[i] );    }
     kmers_plus.resize( starts.back( ) );
     #pragma omp parallel for
     for ( size_t i = 0; i < unibases.size( ); i++ )
     {    const basevector& u = unibases[i];
          String s = u.ToString( );
          kmer<K> x;
          int count = 0;
          for ( int j = 0; j <= u.isize( ) - K; j++ )
          {    if ( s.Contains( prefix, j ) )
               {    int64_t r = starts[i] + count;
                    x.SetToSubOf( u, j ); 
                    kmers_plus[r].first = x;
                    kmers_plus[r].second = i; 
                    kmers_plus[r].third = j;    
                    count++;    }    }    }
     // cout << Date( ) << ": pre -- sorting" << endl;
     ParallelSort(kmers_plus);    
     // cout << Date( ) << ": pre -- sorting complete" << endl;    
          }

template void MakeKmerLookup0Pre( const vecbasevector& unibases,
     const String& prefix, vec< triple<kmer<32>,int,int> >& kmers_plus );

template void MakeKmerLookup0Pre( const vecbasevector& unibases,
     const String& prefix, vec< triple<kmer<40>,int,int> >& kmers_plus );

template void MakeKmerLookup0Pre( const vecbasevector& unibases,
     const String& prefix, vec< triple<kmer<60>,int,int> >& kmers_plus );
