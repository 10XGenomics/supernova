// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

// Lawnmower.  (See also testone/MowOne.cc.)
//
// Look for triples of edges e1, e2, f such that
// * e1 and e2 are followed by f (and f alone);
// * nothing else is before f;
// * nothing is before e2;
// * length(e1) >= length(e2);
// * no edge other then e2 emanates from the source of e2;
// * at most two reads are assigned to e2.
//
// Then examine the reads whose placement includes e1.f or e2.f.
// (Ignore such reads having >= 5 Q30 mismatches in their canonical placement.)
// Compute the relative quality score sum for these reads.  Ignore duplicates,
// as determined by start position.
//
// For each of e1 and e2, drop the must supportive read.
//
// Then sum the relative quality score sum for each of e1, e2, capping the
// contribution of each read at 50.  Suppose that the total for the reads supporting
// e1 is at least 300 and for e2, it is zero.  Then delete e2 and its reverse
// complement.
//
// CHANGED: now ignores dups.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"
#include "feudal/PQVec.h"
#include "math/Functions.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "paths/long/large/GapToyTools.h"
#include "10X/Lawnmower.h"
#include "10X/SecretOps.h"
#include "ParallelVecUtilities.h"
#include "feudal/SubsetMasterVec.h"
#include "feudal/ObjectManager.h"

template< class VB, class VQ, class VP > void MowLawn(
     HyperBasevectorX& hb, vec<int>& inv, VB& bases, VQ& quals,
     VP& paths, String paths_index_file, const vec<int32_t>& bc,
     vec<Bool>& dup, double interdup, vec<int>& dels,
     const vec<int>& tests )
{
     // Heuristics.

     const int MAX_READS2 = 2;
     const int MAX_QS = 50;
     const int MIN_QS_RIGHT = 300;
     const int MAX_QS_WRONG = 30;
     const int MAX_HQ_DIFFS = 5;
     const int MIN_HQ = 30;

/*
     // compute edge support by non-dup reads
     const int batch = 100000;
     vec<int> npe( hb.E( ), -1 );
     for ( int start = 0; start < hb.E( ); start += batch ) {
          const int stop = Min( start+batch, hb.E( ) );
          VecULongVec paths_index_part;
          paths_index_part.ReadRange( paths_index_file, start, stop );
          #pragma omp parallel for
          for ( int e = start; e < stop; e++ ) {
               int tot = 0;
               for ( auto & id : paths_index_part[e-start] )
                    if ( !dup[ id/2 ] )
                         tot++;
               npe[e] = tot;
          }
     }

     // Delete this check below:
     // checks that all the entries were initialized
     for ( auto & n : npe )
          ForceAssertGe ( n, 0 );
     
     // Add in inverse edge support
     #pragma omp parallel for
     for ( int e = 0; e < hb.E( ); e++ ) {
          const int re = inv[e];
          if ( e < re ) {
               const int tot = npe[e] + npe[re];
               npe[e]  = tot;
               npe[re] = tot;
          }
     }
*/     
     // Find the edges whose paths_index we are going to be needing
     // and use to define a SubsetMasterVec
     vec<size_type> eoi;
     for ( int v = 0; v < hb.N( ); v++ )
     {
          if ( hb.To(v).size( ) != 2 ) continue;
          for ( int jj = 0; jj < 2; jj++ )
          {    int e1 = hb.ITo(v,jj), e2 = hb.ITo(v,1-jj);
               if ( tests.nonempty( ) && !Member( tests, e2 ) ) continue;
               int n = hb.Kmers(e2);
               if ( n > hb.Kmers(e1) ) continue;
               if ( hb.To( hb.ToLeft(e2) ).size( ) > 0 ) continue;
               if ( hb.From( hb.ToLeft(e2) ).size( ) > 1 ) continue;
               if ( hb.From(v).size( ) != 1 ) continue;
               int f = hb.IFrom( v, 0 ), re1 = inv[e1], re2 = inv[e2];
               eoi.push_back( e1 );
               eoi.push_back( e2 );
               eoi.push_back( re1 );
               eoi.push_back( re2 );
          }
     }

     ParallelUniqueSort( eoi );
     // define paths_index only over the edges of interest.
     SubsetMasterVec< ULongVec > paths_index_subset( paths_index_file, eoi );
    
     // Go for it.

     cout << Date( ) << ": start loop, mem = " << MemUsageGBString( ) << endl;
     double clock = WallClockTime( );
     dels.clear( );
     int count = 0;
     if ( tests.nonempty( ) ) omp_set_num_threads(1);
     #pragma omp parallel for schedule( dynamic, 100000 )
     for ( int v = 0; v < hb.N( ); v++ )
     {
          /*
          if ( tests.empty( ) && v % 10000000 == 0 )
          {
               #pragma omp critical
               {    count++;
                    cout << Date( ) << ": starting " << count
                         << " of " << hb.E( )/10000000 << endl;    }    }
          */

          if ( hb.To(v).size( ) != 2 ) continue;
          for ( int jj = 0; jj < 2; jj++ )
          {    int e1 = hb.ITo(v,jj), e2 = hb.ITo(v,1-jj);
               if ( tests.nonempty( ) && !Member( tests, e2 ) ) continue;
               int n = hb.Kmers(e2);
               if ( n > hb.Kmers(e1) ) continue;
               if ( hb.To( hb.ToLeft(e2) ).size( ) > 0 ) continue;
               if ( hb.From( hb.ToLeft(e2) ).size( ) > 1 ) continue;
               if ( hb.From(v).size( ) != 1 ) continue;
               int f = hb.IFrom( v, 0 ), re1 = inv[e1], re2 = inv[e2];

               // Make sure we don't have too many reads.

               int npe = 0;
               ULongVec & pi_e2  = paths_index_subset[e2];
               ULongVec & pi_re2 = paths_index_subset[re2];
               for ( int j = 0; j < (int) pi_e2.size( ); j++ )
                    if ( !dup[ pi_e2[j]/2 ] ) npe++;
               for ( int j = 0; j < (int) pi_re2.size( ); j++ )
                    if ( !dup[ pi_re2[j]/2 ] ) npe++;
               if ( npe > MAX_READS2 ) continue;

               basevector x1 = hb.O(e1), x2 = hb.O(e2);
               int del = hb.Kmers(e1) - n;
               x1.SetToSubOf( x1, del, x1.isize( ) - del );
               x1.resize(n), x2.resize(n);
               vec< triple<Bool, int,int> > oq; // (fw, offset, qual sum difference)
               for ( int pass = 1; pass <= 2; pass++ )
               {    int e = ( pass == 1 ? e1 : e2 );
                    int re = inv[e];
                    ULongVec & pi_e  = paths_index_subset[e];
                    ULongVec & pi_re = paths_index_subset[re];

                    for ( int i = 0; i < (int) pi_e.size( ); i++ )
                    {    int64_t id = pi_e[i];
                         if ( dup[id/2] ) continue;
                         const ReadPath& p = paths[id];
                         const basevector& b = bases[id];
                         qualvector q;
                         quals[id].unpack(&q);
                         int offset = p.getOffset( );

                         // Exclude any read that in its standard placement has
                         // >= 5 Q30 mismatches.

                         vec<int> x;
                         for ( int j = 0; j < (int) p.size( ); j++ )
                              x.push_back( p[j] );
                         basevector u = hb.Cat(x);
                         int q30 = 0;
                         for ( int j = 0; j < b.isize( ); j++ )
                         {    int upos = p.getOffset( ) + j;
                              if ( upos < 0 || upos >= u.isize( ) ) continue;
                              if ( q[j] < MIN_HQ ) continue;
                              if ( b[j] != u[upos] ) q30++;    }
                         if ( q30 >= MAX_HQ_DIFFS ) continue;

                         for ( int l = 0; l < (int) p.size( ) - 1; l++ )
                         {    if ( p[l] == e )
                              {    int q1 = 0, q2 = 0;
                                   for ( int m = 0; m < b.isize( ); m++ )
                                   {    int epos = offset + m;
                                        if ( pass == 1 ) epos -= del;
                                        if ( epos < 0 || epos >= n ) continue;
                                        if ( b[m] != x1[epos] ) q1 += q[m];
                                        if ( b[m] != x2[epos] ) q2 += q[m];    }
                                   if ( q1 != q2 )
                                        oq.push( True, offset, q1-q2 );    }
                              offset -= hb.Kmers( p[l] );    }    }
                    for ( int i = 0; i < (int) pi_re.size( ); i++ )
                    {    int64_t id = pi_re[i];
                         if ( dup[id/2] ) continue;
                         const ReadPath& p = paths[id];
                         basevector b = bases[id];
                         qualvector q;
                         quals[id].unpack(&q);

                         // Exclude any read that in its standard placement has
                         // >= 5 Q30 mismatches.

                         vec<int> x;
                         for ( int j = 0; j < (int) p.size( ); j++ )
                              x.push_back( p[j] );
                         basevector u = hb.Cat(x);
                         int q30 = 0;
                         for ( int j = 0; j < b.isize( ); j++ )
                         {    int upos = p.getOffset( ) + j;
                              if ( upos < 0 || upos >= u.isize( ) ) continue;
                              if ( q[j] < MIN_HQ ) continue;
                              if ( b[j] != u[upos] ) q30++;    }
                         if ( q30 >= MAX_HQ_DIFFS ) continue;

                         b.ReverseComplement( );
                         q.ReverseMe( );
                         int offset = p.getOffset( ) - hb.Kmers( p[0] );
                         for ( int l = 1; l < (int) p.size( ); l++ )
                         {    if ( p[l] == re )
                              {    int off = hb.Bases(e) - offset - b.isize( );
                                   int q1 = 0, q2 = 0;
                                   for ( int m = 0; m < b.isize( ); m++ )
                                   {    int epos = off + m;
                                        if ( pass == 1 ) epos -= del;
                                        if ( epos < 0 || epos >= n ) continue;
                                        if ( b[m] != x1[epos] ) q1 += q[m];
                                        if ( b[m] != x2[epos] ) q2 += q[m];    }
                                   if ( q1 != q2 )
                                        oq.push( False, offset, q1-q2 );    }
                              offset -= hb.Kmers( p[l] );    }    }    }
               Sort(oq);
               vec<int> z1, z2;
               for ( int i = 0; i < oq.isize( ); i++ )
               {    int j;
                    for ( j = i + 1; j < oq.isize( ); j++ )
                    {    if ( oq[j].first != oq[i].first
                              || oq[j].second != oq[i].second )
                         {    break;    }    }
                    if ( double(oq[i].third)/double(oq[j-1].third) > 0 )
                    {    if ( oq[i].third > 0 )
                              z2.push_back( Min( oq[j-1].third, MAX_QS ) );
                         else z1.push_back( Min( -oq[i].third, MAX_QS ) );    }
                    i = j - 1;    }
               Sort(z1), Sort(z2);
               if ( z1.nonempty( ) ) z1.pop_back( );
               if ( z2.nonempty( ) ) z2.pop_back( );
               int qss1 = Sum(z1), qss2 = Sum(z2);
               if ( qss1 >= MIN_QS_RIGHT && qss2 <= MAX_QS_WRONG )
               {
                    #pragma omp critical
                    {    dels.push_back( e2, inv[e2] );    }    }
               if ( tests.nonempty( ) ) PRINT4( e1, e2, qss1, qss2 );
                    }    }

     // Could kill large data structures now.

     cout << Date( ) << ": done with main loop, time = " << TimeSince(clock) << endl;
     if ( tests.empty( ) )
     {    cout << Date( ) << ": peak mem usage = " << PeakMemUsageGBString( )
               << endl;
          cout << Date( ) << ": sorting" << endl;    }
     UniqueSort(dels);
     // compute hash
/*     std::string dels_str;
     for ( auto & d : dels )
          dels_str += std::to_string( d );
     int64_t hash_dels = std::hash<std::string>()( dels_str );
     PRINT( hash_dels );*/
     cout << Date( ) << ": " << dels.size( ) << " edges scheduled for deletion"
          << endl;    }

template void MowLawn( HyperBasevectorX&, vec<int>&, VirtualMasterVec<basevector>&,
     VirtualMasterVec<PQVec>&, VirtualMasterVec<ReadPath>&,
     String, const vec<int32_t>&, vec<Bool>& dup, double interdup,
     vec<int>&, const vec<int>& );

template void MowLawn( HyperBasevectorX&, vec<int>&, MasterVec<basevector>&,
     MasterVec<PQVec>&, MasterVec<ReadPath>&, String,
     const vec<int32_t>&, vec<Bool>& dup, double interdup,
     vec<int>&, const vec<int>& );

template void MowLawn( HyperBasevectorX&, vec<int>&, MasterVec<basevector>&,
     VirtualMasterVec<PQVec>&, MasterVec<ReadPath>&, String,
     const vec<int32_t>&, vec<Bool>& dup, double interdup,
     vec<int>&, const vec<int>& );


template< class VB, class VQ> void MowLawn(
     HyperBasevectorX& hb, vec<int>& inv, VB& bases, ObjectManager<VQ>& quals_om,
     ReadPathVecX& paths, String paths_index_file, const vec<int32_t>& bc,
     vec<Bool>& dup, double interdup, vec<int>& dels,
     const vec<int>& tests )
{
     // Heuristics.

     const int MAX_READS2 = 2;
     const int MAX_QS = 50;
     const int MIN_QS_RIGHT = 300;
     const int MAX_QS_WRONG = 30;
     const int MAX_HQ_DIFFS = 5;
     const int MIN_HQ = 30;

/*
     // compute edge support by non-dup reads
     const int batch = 100000;
     vec<int> npe( hb.E( ), -1 );
     for ( int start = 0; start < hb.E( ); start += batch ) {
          const int stop = Min( start+batch, hb.E( ) );
          VecULongVec paths_index_part;
          paths_index_part.ReadRange( paths_index_file, start, stop );
          #pragma omp parallel for
          for ( int e = start; e < stop; e++ ) {
               int tot = 0;
               for ( auto & id : paths_index_part[e-start] )
                    if ( !dup[ id/2 ] )
                         tot++;
               npe[e] = tot;
          }
     }

     // Delete this check below:
     // checks that all the entries were initialized
     for ( auto & n : npe )
          ForceAssertGe ( n, 0 );
     
     // Add in inverse edge support
     #pragma omp parallel for
     for ( int e = 0; e < hb.E( ); e++ ) {
          const int re = inv[e];
          if ( e < re ) {
               const int tot = npe[e] + npe[re];
               npe[e]  = tot;
               npe[re] = tot;
          }
     }
*/     
     // Find the edges whose paths_index we are going to be needing
     // and use to define a SubsetMasterVec
     vec<size_type> eoi;
     for ( int v = 0; v < hb.N( ); v++ )
     {
          if ( hb.To(v).size( ) != 2 ) continue;
          for ( int jj = 0; jj < 2; jj++ )
          {    int e1 = hb.ITo(v,jj), e2 = hb.ITo(v,1-jj);
               if ( tests.nonempty( ) && !Member( tests, e2 ) ) continue;
               int n = hb.Kmers(e2);
               if ( n > hb.Kmers(e1) ) continue;
               if ( hb.To( hb.ToLeft(e2) ).size( ) > 0 ) continue;
               if ( hb.From( hb.ToLeft(e2) ).size( ) > 1 ) continue;
               if ( hb.From(v).size( ) != 1 ) continue;
               int f = hb.IFrom( v, 0 ), re1 = inv[e1], re2 = inv[e2];
               eoi.push_back( e1 );
               eoi.push_back( e2 );
               eoi.push_back( re1 );
               eoi.push_back( re2 );
          }
     }

     ParallelUniqueSort( eoi );
     // define paths_index only over the edges of interest.
     SubsetMasterVec< ULongVec > paths_index_subset( paths_index_file, eoi );

     // Find readpaths of interest based on paths_index_subset
     vec<size_type> poi;
     poi.reserve(paths_index_subset.SizeSum());
     for(const auto& ip: paths_index_subset){
         for(auto& rp: ip){
             poi.push_back(rp);
         }
     }
     ParallelUniqueSort(poi);
     SubsetMasterVec<PQVec> quals_subset(quals_om.filename(), poi);
     Destroy(poi);
    
     // Go for it.

     cout << Date( ) << ": start loop, mem = " << MemUsageGBString( ) << endl;
     double clock = WallClockTime( );
     dels.clear( );
     int count = 0;
     if ( tests.nonempty( ) ) omp_set_num_threads(1);
     #pragma omp parallel for schedule( dynamic, 100000 )
     for ( int v = 0; v < hb.N( ); v++ )
     {
          /*
          if ( tests.empty( ) && v % 10000000 == 0 )
          {
               #pragma omp critical
               {    count++;
                    cout << Date( ) << ": starting " << count
                         << " of " << hb.E( )/10000000 << endl;    }    }
          */

          if ( hb.To(v).size( ) != 2 ) continue;
          for ( int jj = 0; jj < 2; jj++ )
          {    int e1 = hb.ITo(v,jj), e2 = hb.ITo(v,1-jj);
               if ( tests.nonempty( ) && !Member( tests, e2 ) ) continue;
               int n = hb.Kmers(e2);
               if ( n > hb.Kmers(e1) ) continue;
               if ( hb.To( hb.ToLeft(e2) ).size( ) > 0 ) continue;
               if ( hb.From( hb.ToLeft(e2) ).size( ) > 1 ) continue;
               if ( hb.From(v).size( ) != 1 ) continue;
               int f = hb.IFrom( v, 0 ), re1 = inv[e1], re2 = inv[e2];

               // Make sure we don't have too many reads.

               int npe = 0;
               ULongVec & pi_e2  = paths_index_subset[e2];
               ULongVec & pi_re2 = paths_index_subset[re2];
               for ( int j = 0; j < (int) pi_e2.size( ); j++ )
                    if ( !dup[ pi_e2[j]/2 ] ) npe++;
               for ( int j = 0; j < (int) pi_re2.size( ); j++ )
                    if ( !dup[ pi_re2[j]/2 ] ) npe++;
               if ( npe > MAX_READS2 ) continue;

               basevector x1 = hb.O(e1), x2 = hb.O(e2);
               int del = hb.Kmers(e1) - n;
               x1.SetToSubOf( x1, del, x1.isize( ) - del );
               x1.resize(n), x2.resize(n);
               vec< triple<Bool, int,int> > oq; // (fw, offset, qual sum difference)
               for ( int pass = 1; pass <= 2; pass++ )
               {    int e = ( pass == 1 ? e1 : e2 );
                    int re = inv[e];
                    ULongVec & pi_e  = paths_index_subset[e];
                    ULongVec & pi_re = paths_index_subset[re];

                    for ( int i = 0; i < (int) pi_e.size( ); i++ )
                    {    int64_t id = pi_e[i];
                         if ( dup[id/2] ) continue;
                         ReadPath p; paths.unzip(p,hb,id);
                         const basevector& b = bases[id];
                         qualvector q;
                         quals_subset[id].unpack(&q);
                         int offset = p.getOffset( );

                         // Exclude any read that in its standard placement has
                         // >= 5 Q30 mismatches.

                         vec<int> x;
                         for ( int j = 0; j < (int) p.size( ); j++ )
                              x.push_back( p[j] );
                         basevector u = hb.Cat(x);
                         int q30 = 0;
                         for ( int j = 0; j < b.isize( ); j++ )
                         {    int upos = p.getOffset( ) + j;
                              if ( upos < 0 || upos >= u.isize( ) ) continue;
                              if ( q[j] < MIN_HQ ) continue;
                              if ( b[j] != u[upos] ) q30++;    }
                         if ( q30 >= MAX_HQ_DIFFS ) continue;

                         for ( int l = 0; l < (int) p.size( ) - 1; l++ )
                         {    if ( p[l] == e )
                              {    int q1 = 0, q2 = 0;
                                   for ( int m = 0; m < b.isize( ); m++ )
                                   {    int epos = offset + m;
                                        if ( pass == 1 ) epos -= del;
                                        if ( epos < 0 || epos >= n ) continue;
                                        if ( b[m] != x1[epos] ) q1 += q[m];
                                        if ( b[m] != x2[epos] ) q2 += q[m];    }
                                   if ( q1 != q2 )
                                        oq.push( True, offset, q1-q2 );    }
                              offset -= hb.Kmers( p[l] );    }    }
                    for ( int i = 0; i < (int) pi_re.size( ); i++ )
                    {    int64_t id = pi_re[i];
                         if ( dup[id/2] ) continue;
                         ReadPath p; paths.unzip(p,hb,id);
                         basevector b = bases[id];
                         qualvector q;
                         quals_subset[id].unpack(&q);

                         // Exclude any read that in its standard placement has
                         // >= 5 Q30 mismatches.

                         vec<int> x;
                         for ( int j = 0; j < (int) p.size( ); j++ )
                              x.push_back( p[j] );
                         basevector u = hb.Cat(x);
                         int q30 = 0;
                         for ( int j = 0; j < b.isize( ); j++ )
                         {    int upos = p.getOffset( ) + j;
                              if ( upos < 0 || upos >= u.isize( ) ) continue;
                              if ( q[j] < MIN_HQ ) continue;
                              if ( b[j] != u[upos] ) q30++;    }
                         if ( q30 >= MAX_HQ_DIFFS ) continue;

                         b.ReverseComplement( );
                         q.ReverseMe( );
                         int offset = p.getOffset( ) - hb.Kmers( p[0] );
                         for ( int l = 1; l < (int) p.size( ); l++ )
                         {    if ( p[l] == re )
                              {    int off = hb.Bases(e) - offset - b.isize( );
                                   int q1 = 0, q2 = 0;
                                   for ( int m = 0; m < b.isize( ); m++ )
                                   {    int epos = off + m;
                                        if ( pass == 1 ) epos -= del;
                                        if ( epos < 0 || epos >= n ) continue;
                                        if ( b[m] != x1[epos] ) q1 += q[m];
                                        if ( b[m] != x2[epos] ) q2 += q[m];    }
                                   if ( q1 != q2 )
                                        oq.push( False, offset, q1-q2 );    }
                              offset -= hb.Kmers( p[l] );    }    }    }
               Sort(oq);
               vec<int> z1, z2;
               for ( int i = 0; i < oq.isize( ); i++ )
               {    int j;
                    for ( j = i + 1; j < oq.isize( ); j++ )
                    {    if ( oq[j].first != oq[i].first
                              || oq[j].second != oq[i].second )
                         {    break;    }    }
                    if ( double(oq[i].third)/double(oq[j-1].third) > 0 )
                    {    if ( oq[i].third > 0 )
                              z2.push_back( Min( oq[j-1].third, MAX_QS ) );
                         else z1.push_back( Min( -oq[i].third, MAX_QS ) );    }
                    i = j - 1;    }
               Sort(z1), Sort(z2);
               if ( z1.nonempty( ) ) z1.pop_back( );
               if ( z2.nonempty( ) ) z2.pop_back( );
               int qss1 = Sum(z1), qss2 = Sum(z2);
               if ( qss1 >= MIN_QS_RIGHT && qss2 <= MAX_QS_WRONG )
               {
                    #pragma omp critical
                    {    dels.push_back( e2, inv[e2] );    }    }
               if ( tests.nonempty( ) ) PRINT4( e1, e2, qss1, qss2 );
                    }    }

     // Could kill large data structures now.

     cout << Date( ) << ": done with main loop, time = " << TimeSince(clock) << endl;
     if ( tests.empty( ) )
     {    cout << Date( ) << ": peak mem usage = " << PeakMemUsageGBString( )
               << endl;
          cout << Date( ) << ": sorting" << endl;    }
     UniqueSort(dels);
     // compute hash
/*     std::string dels_str;
     for ( auto & d : dels )
          dels_str += std::to_string( d );
     int64_t hash_dels = std::hash<std::string>()( dels_str );
     PRINT( hash_dels );*/
     cout << Date( ) << ": " << dels.size( ) << " edges scheduled for deletion"
          << endl;    }

template void MowLawn( HyperBasevectorX&, vec<int>&, VirtualMasterVec<basevector>&,
     ObjectManager<VirtualMasterVec<PQVec> >&, ReadPathVecX&,
     String, const vec<int32_t>&, vec<Bool>& dup, double interdup,
     vec<int>&, const vec<int>& );

template void MowLawn( HyperBasevectorX&, vec<int>&, MasterVec<basevector>&,
     ObjectManager<MasterVec<PQVec> >&, ReadPathVecX&, String,
     const vec<int32_t>&, vec<Bool>& dup, double interdup,
     vec<int>&, const vec<int>& );

template void MowLawn( HyperBasevectorX&, vec<int>&, MasterVec<basevector>&,
     ObjectManager<VirtualMasterVec<PQVec> >&, ReadPathVecX&, String,
     const vec<int32_t>&, vec<Bool>& dup, double interdup,
     vec<int>&, const vec<int>& );
