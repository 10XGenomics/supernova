///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
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
#include "paths/long/MakeAlignments.h"
#include "paths/long/FriendAligns.h"

template<int K> class gurgle {

     public:

     gurgle( ) { }
     gurgle( const kmer<K>& mer, const signed char prev,
          const int id, const int pos ) : mer(mer), prev(prev), id(id), pos(pos) { }

     friend Bool operator<( const gurgle& x, const gurgle& y )
     {    if ( x.mer < y.mer ) return True;
          if ( x.mer > y.mer ) return False;
          if ( x.prev < y.prev ) return True;
          if ( x.prev > y.prev ) return False;
          if ( x.id < y.id ) return True;
          if ( x.id > y.id ) return False;
          return x.pos < y.pos;    }

     kmer<K> mer;
     signed char prev;
     int id;
     int pos;

};

template<int K> void MakeAlignments( const int max_freq, const vecbasevector& bases,
     const vec<Bool>& is_target, vec<simple_align_data>& aligns_all, int verbosity )
{
     // Build kmer lookup table 'kmers_plus' from the reads and their reverse
     // complements.  An entry in this table is (kmer,prev,id,pos), where prev is
     // the base before the kmer, or -1 if at the beginning of the read.

     vecbasevector bases_rc(bases);
     for ( int64_t id = 0; id < (int64_t) bases.size( ); id++ )
          bases_rc[id].ReverseComplement( );
     vecbasevector bases2(bases);
     bases2.Append(bases_rc);
     if ( verbosity ) cout << Date( ) << ": building starts" << endl;
     vec< gurgle<K> > kmers_plus;
     vec<int64_t> starts;
     starts.push_back(0);
     for ( size_t i = 0; i < bases2.size( ); i++ )
     {    const basevector& u = bases2[i];
          starts.push_back( starts.back( ) + Max( 0, u.isize( ) - K + 1 ) );    }
     if ( verbosity ) cout << Date( ) << ": making kmers_plus" << endl;
     kmers_plus.resize( starts.back( ) );
     #pragma omp parallel for
     for ( size_t i = 0; i < bases2.size( ); i++ )
     {    const basevector& u = bases2[i];
          for ( int j = 0; j <= u.isize( ) - K; j++ )
          {    int64_t r = starts[i] + j;
               kmers_plus[r].mer.SetToSubOf( u, j );
               kmers_plus[r].prev = ( j == 0 ? -1 : bases2[i][j-1] );
               kmers_plus[r].id = i;
               kmers_plus[r].pos = j;    }    }
     if ( verbosity ) cout << Date( ) << ": sorting kmers_plus, size = "
          << ToStringAddCommas( kmers_plus.size( ) ) << endl;
     ParallelSort(kmers_plus);

     // Define batches for alignment.

     if ( verbosity ) cout << Date( ) << ": defining bstart" << endl;
     int nproc = 48;
     vec<int64_t> bstart(nproc+1);
     for ( int i = 0; i <= nproc; i++ )
          bstart[i] = ( (int64_t) kmers_plus.size( ) * i ) / nproc;
     for ( int i = 1; i < nproc; i++ )
     {    int64_t& s = bstart[i];
          while( s > bstart[i-1] && kmers_plus[s].mer == kmers_plus[s-1].mer)
          {    s--;    }    }

     // Align.

     if ( verbosity ) cout << Date( ) << ": start main alignment loop" << endl;
     vec< vec<simple_align_data> > aligns( bstart.size( ) - 1 );
     #pragma omp parallel for schedule(dynamic, 1)
     for ( int i = 0; i < bstart.isize( ) - 1; i++ )
     {    int64_t start = bstart[i], stop = bstart[i+1];
          while( start < stop )
          {    vec< pair<int64_t,int64_t> > locs( 5, make_pair(-1,-1) ); // [start,stop) pairs for each preceeding bases 
               int64_t start0 = start;
               while(1)
               {    if ( locs[ kmers_plus[start].prev+1 ].first == -1 )
                       locs[ kmers_plus[start].prev+1 ].first = start;
                    locs[ kmers_plus[start].prev+1 ].second = start+1;
                    start++;
                    if ( start == stop 
                         || kmers_plus[start].mer != kmers_plus[start-1].mer )
                    {    break;    }    }
               if ( start - start0 > max_freq ) continue; // This will not do exactly what we wanted to do. To be fixed.
               // Find the alignments. To avoid duplication only aligned using the 
               // kmers where the previous kmer differs or they are all at the beginning
               // of reads.
               for ( int i1 = 0; i1 < 5; i1++ )
               for ( int i2 = 0; i2 < 5; i2++ )
               {    if ( i2 == i1 && i1 > 0 && i2 > 0 ) continue;
                    for ( int64_t j1 = locs[i1].first; j1 < locs[i1].second; j1++ )
                    {    int id1 = kmers_plus[j1].id; 
                         Bool rc1 = ( id1 >= (int) bases.size( ) );
                         if ( rc1 || !is_target[id1] ) continue;
                         int pos1 = kmers_plus[j1].pos;
                         for ( int64_t j2 = locs[i2].first; j2 < locs[i2].second; j2++ )
                         {    int id2 = kmers_plus[j2].id; 
                              int pos2 = kmers_plus[j2].pos;
                              Bool rc2 = ( id2 >= (int) bases.size( ) );
                              if (rc2) id2 -= (int) bases.size( );
                              if ( id2 == id1 ) continue;  // Don't align the read to itself
                              aligns[i].push( id1, id2, 
                                   pos1 - pos2, rc2 );    }    }    }    }    }

     // Merge and sort alignments.

     if ( verbosity ) cout << Date( ) << ": merging" << endl;
     aligns_all.clear( );
     for ( int i = 0; i < bstart.isize( ) - 1; i++ )
          aligns_all.append( aligns[i] );
     if ( verbosity ) cout << Date( ) << ": before unique-sorting there are "
          << ToStringAddCommas( aligns_all.size( ) ) << " alignments" << endl;
     if ( verbosity ) cout << Date( ) << ": sorting" << endl;
     ParallelUniqueSort(aligns_all);
     if ( verbosity ) cout << Date( ) << ": after unique-sorting there are "
          << ToStringAddCommas( aligns_all.size( ) ) << " alignments" << endl;    }

void MakeAlignments( const int K, const int max_freq, const vecbasevector& bases,
     const vec<Bool>& is_target, vec<simple_align_data>& aligns_all, int verbosity )
{    if ( K == 8 ) MakeAlignments<8>(max_freq, bases, is_target, aligns_all, verbosity);
     else if ( K == 12 ) MakeAlignments<12>(max_freq, bases, is_target, aligns_all, verbosity);
     else if ( K == 16 ) MakeAlignments<16>(max_freq, bases, is_target, aligns_all, verbosity);
     else if ( K == 20 ) MakeAlignments<20>(max_freq, bases, is_target, aligns_all, verbosity);
     else if ( K == 24 ) MakeAlignments<24>(max_freq, bases, is_target, aligns_all, verbosity);
     else if ( K == 28 ) MakeAlignments<28>(max_freq, bases, is_target, aligns_all, verbosity);
     else if ( K == 40 ) MakeAlignments<40>(max_freq, bases, is_target, aligns_all, verbosity);
     else if ( K == 60 ) MakeAlignments<60>(max_freq, bases, is_target, aligns_all, verbosity);
     else if ( K == 80 ) MakeAlignments<80>(max_freq, bases, is_target, aligns_all, verbosity);
     else
     {    FatalErr("\nIllegal K value for MakeAlignments.");    }    }
