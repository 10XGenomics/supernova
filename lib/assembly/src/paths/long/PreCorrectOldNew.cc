///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// PreCorrectOldNew.
// 
// Based on: 
// svn cat https://
// svn.broadinstitute.org/comprd/trunk/BroadCRD/paths/PreCorrectOld.cc/@32148.
//
// This is likely less efficient than pre_correct_parallel and moreover there may have
// there may have been bugfixes in the interim, in particular (r38538) the last kmer 
// might not be processed.
//
// Algorithm (initial version).
// * Find all 25-mers in the reads.  Think of each 25-mer as 12-1-12, and
// rewrite as 12-12-1, build a kmer record, and sort.  We ignore 25-mers whose
// first 12-mer equals the reverse complement of its last 12-mer.
// * Go through each stack of 25-mers for which the first 24 bases agree.
// Stacks of size < 6 are ignored.  Examine the last bases (i.e. the original
// midle bases) and their quality scores.  Find the quality score sum (qss) for
// each of A, C, G, T.  To edit a loser base, we require that the winner qss is
// at least 60 and that the loser qss is < (winner qss)/4.  We also require
// that the loser is not supported by two bases of quality 20+.
// * Do not make changes that are within 12 bases of other changes.
// * Change the base and lower the quality score to zero.
//
// Algorithm (revised version, in addition to the above).
// * For each stack, form the full stack of the reads, and look for high-quality
//   differences between winning and losing reads.  The requirement is Q >= 30 on
//   both reads.
// * If there are at least two columns, and for each such column, at least two losing
//   and at least two winning reads that exhibit such a high-quality difference, do
//   not make the change.
// !!!!!!!!! CHANGED !!!!!!!!!
//
// Problems.
// * Broken for large read ids.
// * Computational performance needs to be improved.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "CoreTools.h"
#include "ParallelVecUtilities.h"
#include "ParseSet.h"
#include "Qualvector.h"
#include "kmers/KmerRecord.h"
#include "paths/long/PreCorrectOldNew.h"
#include "Set.h"

template<int L> void PreCorrectCore( const int flank, const int u, const int l, 
     const vecbasevector& bases, const vecqualvector& quals,
     const vec< kmer_record<L,1> >& recs, ulonglong i, ulonglong j, 
     const vec<int>& bs, const vec<int>& total, vecbasevector& bases_new, 
     const vec<int>& trace_ids )
{
     Bool print = False;
     for ( ulonglong k = i; k < j; k++ )
     {    int id = recs[k].GetId( );
          if ( BinMember( trace_ids, id ) ) print = True;    }

     if (print) 
     {    cout << "\n";
          PRINT2( u, total[u] );    }
     int nrows = 0, nrows1 = 0, left_ext = 0, right_ext = 0;
     for ( ulonglong k = i; k < j; k++ )
     {    int id = recs[k].GetId( ), pos = recs[k].GetPos( );
          char middle;
          if ( pos >= 0 ) middle = bases[id][ pos + flank ];
          else middle = 3 - bases[id][ -pos - 1 + flank ];
          if ( middle != bs[u] && middle != bs[0] ) continue;
          nrows++;
          int qpos = ( pos >= 0 ? pos : 
               bases[id].isize( ) - (-pos-1) + 1 - 2*flank-2 );
          left_ext = Max( left_ext, qpos );
          right_ext = Max( right_ext, 
               bases[id].isize( ) - qpos - (2*flank+1) );    }
     int ncols = left_ext + 2*flank+1 + right_ext;

     vec< vec<char> > pbases( nrows, vec<char>( ncols, ' ' ) );
     vec< vec<signed char> > qbases( nrows, vec<signed char>( ncols, -1 ) );
     vec<int> ids;

     int row = 0;
     for ( int pass = 1; pass <= 2; pass++ )
     {    if ( pass == 2 ) nrows1 = row;
          for ( ulonglong k = i; k < j; k++ )
          {    int id = recs[k].GetId( ), pos = recs[k].GetPos( );
               char middle;
               if ( pos >= 0 ) middle = bases[id][ pos + flank ];
               else middle = 3 - bases[id][ -pos - 1 + flank ];
               if ( pass == 1 && middle != bs[u] ) continue;
               if ( pass == 2 && middle != bs[0] ) continue;
               ids.push_back(id);
               int qpos = ( pos >= 0 ? pos : 
                    bases[id].isize( ) - (-pos-1) + 1 - 2*flank-2 );
               if ( pos >= 0 )
               {    for ( int z = 0; z < bases[id].isize( ); z++ )
                    {    pbases[row][ z + left_ext - qpos ] 
                              = as_base( bases[id][z] );
                         qbases[row][ z + left_ext - qpos ] 
                              = quals[id][z];    }    }
               else
               {    basevector b = bases[id];
                    b.ReverseComplement( );
                    qualvector q = quals[id];
                    q.ReverseMe( );
                    for ( int z = 0; z < b.isize( ); z++ )
                    {    pbases[row][ z + left_ext - qpos ] = as_base( b[z] );
                         qbases[row][ z + left_ext - qpos ] = q[z];    }    }
               row++;    }    }

     if (print)
     {    cout << String( left_ext, ' ' ) << String( flank, '*' ) << " "
               << String( flank, '*' ) << "\n";
          for ( int r = 0; r < nrows; r++ )
          {    for ( int c = 0; c < ncols; c++ )
                    cout << pbases[r][c];
               cout << "\n";     }    }

     if (print) cout << "\nhigh-quality differences:\n";
     const int hq = 30;
     vec<int> badc;
     const int max_bad = 10;
     vec< StdSet<int> > diffs1(ncols), diffs2(ncols);
     vec<int> del( nrows, 0 );
     for ( int r1 = 0; r1 < nrows1; r1++ )
     {    // if ( badc.isize( ) > max_bad ) break;
          for ( int r2 = nrows1; r2 < nrows; r2++ )
          {    // if ( badc.isize( ) > max_bad ) break;
               for ( int c = 0; c < ncols; c++ )
               {    if ( pbases[r1][c] != pbases[r2][c] 
                         && qbases[r1][c] >= hq && qbases[r2][c] >= hq )
                    {    if ( c != left_ext + flank ) del[r2]++;
                         diffs1[c].insert(r1);
                         diffs2[c].insert(r2);
                         if ( !Member( badc, c ) )
                         {    badc.push_back(c);
                              if ( badc.isize( ) == max_bad + 1 )
                              {    if (print) cout << "...\n";
                                   // break;    
                                        }    }
                         if ( print && badc.isize( ) <= max_bad )
                              PRINT5( r1, ids[r1], r2, ids[r2], c );    
                              }    }    }    }
     int cbad = 0;
     for ( int c = 0; c < ncols; c++ )
          if ( diffs1[c].size( ) >= 2 && diffs2[c].size( ) >= 2 ) cbad++;
     if (print) PRINT2( nrows1, cbad );
     
     // if ( nrows1 == 1 && badc.size( ) >= 2 ) return;
     if ( cbad >= 10 ) return;

     /*
     int loser_qss = 0, winner_qss = 0;
     for ( int r1 = 0; r1 < nrows1; r1++ )
          loser_qss += qbases[r1][left_ext+flank];
     for ( int r2 = nrows1; r2 < nrows; r2++ )
          if ( del[r2] < 2 ) winner_qss += qbases[r2][left_ext+flank];
     if ( winner_qss < 60 ) return;
     if ( loser_qss >= winner_qss/4 ) return;
     */

     if (print) cout << "\n";

     for ( ulonglong k = i; k < j; k++ )
     {    int id = recs[k].GetId( ), pos = recs[k].GetPos( );

          int qmiddle;
          if ( pos >= 0 ) qmiddle = quals[id][ pos + flank ];
          else qmiddle = quals[id][ -pos - 1 + flank ];

          char middle;
          if ( pos >= 0 ) middle = bases[id][ pos + flank ];
          else middle = 3 - bases[id][ -pos - 1 + flank ];
          if ( middle != bs[u] ) continue;

          if (print) PRINT5( l, id, pos, as_base(middle), qmiddle );

          // #pragma omp critical
          {
          if ( pos >= 0 )
               bases_new[id].Set( pos + flank, bs[0] );
          else
               bases_new[id].Set( -pos-1+flank, 3 - bs[0] );    }    }    }

void PreCorrectOldNew( vecbasevector* bases, vecqualvector const& quals0,
     const vec<int>& trace_ids )
{
     // Define constants.

     const int flank = 12;
     const int L = 2*flank;
     const uint min_stack = 6;
     const int min_qss_winner = 60;
     const double qss_ratio_cutoff = 0.25;
     const int qgood = 20;
     const int qgoods = 2;

     // Execute the algorithm.

     int K = 1 + 2*flank;
     vecbasevector bases_new(*bases);
     longlong nreads = (*bases).size( );
     vec< kmer_record<2*flank,1> > recs;
     if ( !(*bases).empty( ) )
          recs.reserve( nreads * (longlong) ( (*bases)[0].isize( ) - K + 1 ) );
     basevector bi(2*flank), bir(2*flank);
     kmer_record<2*flank,1> rec;
     for ( longlong id = 0; id < nreads; id++ )
     {   const basevector& R = (*bases)[id];
         for ( int j = 0; j <= R.isize( ) - K; j++ )
         {   for ( int l = 0; l < flank; l++ )
             {    bi.Set( l, R[ j + l ] );
                  bi.Set( l + flank, R[ j + l + flank + 1 ] );     }
             bir.ReverseComplement(bi);
             if ( bi == bir ) continue;
             if ( bi < bir ) rec.Set( bi, id, j );
             else rec.Set( bir, id, -j-1 );
             recs.push_back(rec);    }    }
     ParallelSort( recs, kmer_record<L,1>::cmp_basevector );
     vec<ulonglong> start, stop;
     for ( ulonglong i = 0; i < recs.size( ); i++ )
     {    ulonglong j;
          for ( j = i + 1; j < recs.size( ); j++ )
              if ( !recs[j].EqualKmers( recs[i] ) ) break;
          start.push_back(i), stop.push_back(j);
          i = j - 1;    }
     #pragma omp parallel for
     for ( ulonglong l = 0; l < start.size( ); l++ )
     {    ulonglong i = start[l], j = stop[l];
          if ( j - i < min_stack )
          {    i = j - 1;
               continue;    }
          vec<int> qsum(4), bs(4), goods(4), total(4);
          for ( int u = 0; u < 4; u++ )
          {    qsum[u] = goods[u] = 0;
               bs[u] = u;    }
          for ( ulonglong k = i; k < j; k++ )
          {    int id = recs[k].GetId( ), pos = recs[k].GetPos( );
               char middle;
               int middle_q;
               if ( pos >= 0 )
               {    middle = (*bases)[id][ pos + flank ];
                    middle_q = quals0[id][ pos + flank ];    }
               else
               {    middle = 3 - (*bases)[id][ -pos - 1 + flank ];
                    middle_q = quals0[id][ -pos - 1 + flank ];    }
               total[middle]++;
               qsum[middle] += middle_q;
               if ( middle_q >= qgood ) ++goods[middle];    }
           ReverseSortSync( qsum, bs, goods, total );
           if ( qsum[0] >= min_qss_winner )
           {    
                #pragma omp critical
                for ( int u = 1; u < 4; u++ )
                {    if ( goods[u] < qgoods 
                          && total[u] > 0
                          && qsum[u]
                          < int( floor( qss_ratio_cutoff * double(qsum[0]) ) ) )
                     {    
                          PreCorrectCore( flank, u, l, *bases, quals0, recs,
                                          i, j, bs, total, bases_new,
                                          trace_ids );    }    }    }

           i = j - 1;    }
      vec<int> diffs;
      for ( longlong id = 0; id < nreads; id++ )
      {    basevector &R = (*bases)[id], &Rnew = bases_new[id];
           diffs.clear( );
           for ( int j = 0; j < R.isize( ); j++ )
                if ( R[j] != Rnew[j] ) diffs.push_back(j);
           for ( int i = 0; i < diffs.isize( ); i++ )
           {    if ( i > 0 && diffs[i] - diffs[i-1] <= flank ) continue;
                if ( i < diffs.isize( ) - 1 && diffs[i+1] - diffs[i] <= flank )
                     continue;
                /*
                cout << id << '\t' << diffs[i] << '\t' <<
                        Base::val2Char(R[diffs[i]]) <<
                        '\t' << Base::val2Char(Rnew[diffs[i]]) << '\n';
                */
                R.Set( diffs[i], Rnew[ diffs[i] ] );    }    }
}
