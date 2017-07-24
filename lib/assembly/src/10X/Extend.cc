// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS


#include "CoreTools.h"
#include "10X/Extend.h"
#include "feudal/PQVec.h"
#include "graph/DigraphTemplate.h"
#include "paths/HyperBasevector.h"
#include "paths/long/large/GapToyTools.h"

template <class VQ>
void ExtendPathsNew( const HyperBasevectorX& hb, const vec<int>& inv,
     const vecbasevector& bases, const VQ& qualsx, ReadPathVecX& paths,
     const Bool BACK_EXTEND )
{
     // Ancillary stuff.

     int K = hb.K( );

     // Extend path backwards if possible.

     if (BACK_EXTEND)
     {    cout << Date( ) << ": extending paths backward" << endl;
          ReadPathVecX npaths;
          npaths.reserve(paths.size());
          #pragma omp parallel for ordered schedule(dynamic, 1)
          for ( int64_t id = 0; id < (int64_t) paths.size( ); id++ )
          {    ReadPath p;  paths.unzip(p,hb,id);
               // recompute p
               if(p.size()!=0){
                   int offset = p.getOffset( );
                   while( offset < 0 && hb.To( hb.ToLeft( p[0] ) ).nonempty( ) )
                   {    int v = hb.ToLeft( p[0] );
                        Bool extended = False;
                        for ( int j = 0; j < (int) hb.To(v).size( ); j++ )
                        {    int e = hb.ITo( v, j );
                             Bool mismatch = False;
                             int n = hb.Kmers(e);
                             for ( int l = 0; l < hb.Bases(e); l++ )
                             {    if ( l - offset - n < 0 ) continue;
                                  if ( hb.O(e)[l] != bases[id][ l - offset - n ] )
                                  {    mismatch = True;
                                       break;    }    }
                             if (mismatch) continue;
                             extended = True;
                             p.push_front(e);
                             p.setOffset( p.getOffset( ) + n );
                             offset += n;
                             break;    }
                        if ( !extended ) break;    
                   }
               }
               #pragma omp ordered
               npaths.append(p,hb);
          }    
          Destroy(paths);
          paths = npaths;
     }

     // Extend paths 1 (to right).

     const int MIN_GAIN = 5;
     const int EXT_MODE = 1;
     cout << Date( ) << ": extending paths" << endl;
     auto qualsx_clone = qualsx;
     int64_t batch = paths.size()/(omp_get_max_threads()*100)+100;
     vec<int64_t> starts;
     vec<ReadPathVecX> blocks;
     #pragma omp parallel for firstprivate(qualsx_clone) schedule(dynamic,1)
     for ( int64_t start = 0; start < (int64_t) paths.size( ); start+=batch )
     {  
          int64_t stop = Min(start+batch,paths.size());
          ReadPath path;
          ReadPathVecX subRPVX;
          subRPVX.reserve(stop-start);
          for(int64_t i = start; i<stop; i++){
              paths.unzip(path,hb,i);
              if ( path.size( ) > 0 ) path.resize(1);
              ExtendPath( path, i, hb, bases[i], qualsx_clone[i],
                  MIN_GAIN, False, EXT_MODE );    
              subRPVX.append(path,hb);
         }
         #pragma omp critical
         {   starts.push_back(start);
             blocks.push_back(subRPVX); }
     }
     SortSync(starts,blocks);
     Destroy(paths);
     paths.append(blocks);
     Destroy(blocks);
     Destroy(starts);
     Validate( hb, paths );

     // Extend paths 2.
     cout << Date( ) << ": starting secondary path extension, mem = " 
          << MemUsageGBString( ) << endl;
     #pragma omp parallel for schedule(dynamic,1)
     for ( int64_t start = 0; start<(int64_t)paths.size(); start+=batch)
     {    
          int64_t stop = Min(start+batch,paths.size());
          ReadPath p;
          ReadPathVecX subRPVX;
          subRPVX.reserve(stop-start);
          for(int64_t id = start; id<stop; id++){
              paths.unzip(p,hb,id);
              if ( p.size( ) != 0 ){
                  // Do last K bases of path match perfectly?
                  int pp = 0;
                  int e = p[pp];
                  int rpos = 0, epos = p.getOffset( );
                  if ( epos < 0 )
                  {    rpos = -epos;
                       epos = 0;    }
                  int mismatch = -1;
                  const basevector& b = bases[id];
                  while(1)
                  {    if ( b[rpos] != hb.O(e)[epos] )
                       {    mismatch = rpos;
                            break;   }
                       rpos++;
                       if ( rpos == b.isize( ) ) break;
                       epos++;
                       if ( epos == hb.O(e).isize( ) )
                       {    pp++;
                            if ( pp == (int) p.size( ) ) break;
                            e = p[pp];
                            epos = K - 1;    }    }
                  if (!( mismatch >= 0 && mismatch > rpos - K ) && !(rpos == b.isize( ))){
                      // Extend.
                      while(1)
                      {    if ( epos < hb.O(e).isize( ) )
                           {    if ( b[rpos] != hb.O(e)[epos] ) break;    }
                           else
                           {    int v = hb.ToRight(e);
                                Bool ext = False;
                                for ( int j = 0; j < (int) hb.From(v).size( ); j++ )
                                {    int f = hb.IFrom( v, j );
                                     if ( b[rpos] == hb.O(f)[K-1] )
                                     {    e = f;
                                          epos = K-1;
                                          p.push_back(e);
                                          ext = True;    }    }
                                if ( !ext ) break;    }
                           rpos++;
                           epos++;
                           if ( rpos == b.isize( ) ) break;    }

                      /*
                      if ( (int) p.size( ) > np )
                      {
                           #pragma omp critical
                           {    cout << "\nid = " << id << ", p = " << printSeq(p) << endl;
                                cout << "inv(p) = ";
                                for ( int j = (int) p.size( ) - 1; j >= 0; j-- )
                                {    if ( j < (int) p.size( ) - 1 ) cout << ",";
                                     cout << inv[ p[j] ];    }
                                cout << endl;    }    }
                      */
                  }
              }
              subRPVX.append(p,hb);
          }
          #pragma omp critical 
          {   starts.push_back(start);
              blocks.push_back(subRPVX); }
     }
     cout << Date( ) << ": setting paths, mem = "
          << MemUsageGBString( ) << endl;
     SortSync(starts,blocks);
     Destroy(paths);
     paths.append(blocks);
     Destroy(blocks);
     Destroy(starts);
     Validate( hb, paths );    }

template void ExtendPathsNew( const HyperBasevectorX& hb, const vec<int>& inv,
     const vecbasevector& bases, const VecPQVec& qualsx, ReadPathVecX& paths,
     const Bool BACK_EXTEND );

template void ExtendPathsNew( const HyperBasevectorX& hb, const vec<int>& inv,
     const vecbasevector& bases, const VirtualMasterVec<PQVec>& qualsx, ReadPathVecX& paths,
     const Bool BACK_EXTEND );
