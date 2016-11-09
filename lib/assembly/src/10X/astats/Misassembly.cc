// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"
#include "VecUtilities.h"
#include "math/HoInterval.h"
#include "10X/astats/Misassembly.h"

void Misassembly( const SerfVec< quad<int,Bool,ho_interval,ho_interval> >& view,
     int64_t& total_err_dis_num, int64_t& total_err_dis_den,
     int64_t& total_err_ori_num, int64_t& total_err_ori_den,
     int64_t& total_err_ord_num, int64_t& total_err_ord_den )
{
     // Start to compute the core error statistics.  First define the chromosomal
     // interval chr:start-stop that we'll call the best home for the scaffold.

     const int TOO_FAR = 300000;
     vec< triple<int,int,int> > homer; // {(g,start,len)}
     for ( int i = 0; i < (int) view.size( ); i++ )
     {    homer.push( view[i].first,
               view[i].third.Start( ), view[i].third.Length( ) );    }
     Sort(homer);
     vec< quad<int,int,int,int> > inters;
     for ( int i = 0; i < homer.isize( ); i++ )
     {    int chr = homer[i].first, start = homer[i].second;
          int j, stop = homer[i].second + homer[i].third;
          int mass = homer[i].third;
          for ( j = i + 1; j < homer.isize( ); j++ )
          {    if ( homer[j].first != chr ) break;
               if ( homer[j].second - stop > TOO_FAR ) break;
               stop = Max( stop, homer[j].second + homer[j].third );
               mass += homer[j].third;    }
          inters.push( mass, chr, start, stop );
          i = j - 1;    }
     ReverseSort(inters);
     int chr = -1, start = -1, stop = -1;
     int gxlen = -1, gxlen_kmers = -1;
     if ( inters.nonempty( ) )
     {    chr = inters[0].second;
          start = inters[0].third, stop = inters[0].fourth;
          if ( stop > start )
          {    gxlen = stop - start;
               /*
               #pragma omp critical
               {    glens.push_back(gxlen);    }
               */
                     }    }

     // Now compute the distant misjoin error statistics.

     int64_t err_dis_num = 0, err_dis_den = 0;
     for ( int i = 0; i < (int) view.size( ); i++ )
     {    err_dis_den += view[i].third.Length( );
          if ( view[i].first != chr ) err_dis_num += view[i].third.Length( );
          else
          {    if ( view[i].third.Start( ) < start
                    || view[i].third.Start( ) + view[i].third.Length( ) > stop )
               {    err_dis_num += view[i].third.Length( );    }
               /*
               else gxlen_kmers += view[i].third.Length( );
               */
                    }    }
     #pragma omp critical
     {    total_err_dis_num += err_dis_num, total_err_dis_den += err_dis_den;    }

     // Next deal with orientation.

     int64_t fwn = 0, rcn = 0;
     for ( int i = 0; i < (int) view.size( ); i++ )
     {    if ( view[i].first != chr ) continue;
          if ( view[i].third.Start( ) < start || view[i].third.Stop( ) > stop )
               continue;
          if ( view[i].second ) fwn += view[i].third.Length( );
          else rcn += view[i].third.Length( );    }
     Bool is_fw = ( fwn >= rcn );
     int64_t err_ori_num = ( is_fw ? rcn : fwn ), err_ori_den = fwn + rcn;
     #pragma omp critical
     {    total_err_ori_num += err_ori_num, total_err_ori_den += err_ori_den;    }

     // Finally deal with order.  In principle, the misordered set should be the
     // smallest set that can be removed, yields a perfectly ordered set.  Here
     // smallest is measured by kmer count.  Order is only tested on the parts
     // that survive the dis and ori tests.
     //
     // First create the list of all edge alignments under consideration.

     int64_t err_ord_num = 0, err_ord_den = 0;
     vec< triple<int,int,int> > vord; // (pos,kmers,id)
     for ( int i = 0; i < (int) view.size( ); i++ )
     {    if ( view[i].first != chr ) continue;
          if ( view[i].second != is_fw ) continue;
          if ( view[i].third.Start( ) < start ) continue;
          if ( view[i].third.Stop( ) > stop ) continue;
          int start = view[i].third.Start( );
          if ( !is_fw ) start = -start;
          vord.push( start, view[i].third.Length( ), vord.size( ) );    }

     // Now define blocks to be sets whose order will not change.

     for ( int i = 0; i < vord.isize( ); i++ )
          err_ord_den += vord[i].second;
     vec< triple<int,int,int> > vords(vord);
     Sort(vords);
     vec< quad< pair<int,int>, int, int, vec<int> > > blocks;
     for ( int i = 0; i < vords.isize( ); i++ )
     {    int j, nkmers = vords[i].second;
          vec<int> es = { vords[i].third };
          for ( j = i + 1; j < vords.isize( ); j++ )
          {    if ( vords[j].third != vords[j-1].third + 1 ) break;
               nkmers += vords[j].second;
               es.push_back( vords[j].third );    }
          blocks.push( make_pair( vords[i].third, vords[j-1].third ),
               vords[i].first, nkmers, es );
          i = j - 1;    }
     Sort(blocks);
     // PRINT_TO( cout, blocks.size( ) );

     // Rank each block by the total kmers that are out of order with it.
     // Kill the worst block.  Iterate.

     while( blocks.nonempty( ) )
     {    vec<int> mis( blocks.size( ), 0 );
          vec<int> ids( blocks.size( ), vec<int>::IDENTITY );
          for ( int i = 0; i < blocks.isize( ); i++ )
          {    for ( int j = 0; j < i; j++ )
               {    if ( blocks[j].second > blocks[i].second )
                         mis[i] += blocks[j].third;    }
               for ( int j = i + 1; j < blocks.isize( ); j++ )
               {    if ( blocks[j].second < blocks[i].second )
                         mis[i] += blocks[j].third;    }    }
          ReverseSortSync( mis, ids );
          if ( mis[0] == 0 ) break;
          err_ord_num += blocks[ ids[0] ].third;
          blocks.erase( blocks.begin( ) + ids[0] );    }
     #pragma omp critical
     {    total_err_ord_num += err_ord_num, total_err_ord_den += err_ord_den;    }
          }
