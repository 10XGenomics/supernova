// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>

#include "CoreTools.h"
#include "Set.h"
#include "graph/DigraphTemplate.h"
#include "10X/Super.h"
#include "10X/mergers/ExtremeAsm.h"
#include "10X/mergers/NicePrints.h"
#include "10X/mergers/GlueGraphs.h"
#include "10X/mergers/GlueGraphs_HP.h"
#include "10X/mergers/EdgeSupport.h"
#include "10X/micro/DisplayGraphs.h"

// Templated Zipper based on c2g-gluing instead of old MakeMerges
// potentially SLOW if called on large graphs
template <class T>
void BabyZipper(digraphE<vec<T>>& D, vec<int>& dinv, T nE, const int max_pass, const int min_zips ){
    digraphE<vec<T>> D_old;
    vec<int> dinv_old;
    vec<int> to_right;
    vec<int> to_left;

    int npass = 0;
    while(npass<max_pass){
        digraphE<vec<T>> DD;
        vec<int> ddinv;
        D_old = D;
        dinv_old = dinv;
        D_old.ToRight(to_right);
        D_old.ToLeft(to_left);
        vec<nptuple> Matches;
        // get baby zips
        // schedule(dynamic,100) doesn't help much
        #pragma omp parallel for
        for (int v = 0; v < D_old.N(); v++){
            if (D_old.From(v).size()<2) continue;
            // consider all pairs of edges and possible overlaps
            auto& edges = D_old.FromEdgeObj(v);

            for (int d1 = 0; d1<edges.isize(); d1++){
                for (int d2 = 0; d2<edges.isize(); d2++){
                    int i1 = edges[d1];
                    int i2 = edges[d2];
                    if(d1==d2 || D_old.O(i1).size()>D_old.O(i2).size()) continue;
                    int ri1 = dinv_old[i1];
                    int ri2 = dinv_old[i2];
                    if(!IsUnique(i1,i2,ri1,ri2))
                        continue;

                    // pathologies of inversions
                    /* int rv = to_right[ri1]; */
                    /* // avoid pathologies for bubbles */
                    /* if(to_right[i1]==to_right[i2]){ */
                    /*     if(!IsUnique(vec<int>{v,to_right[i1],rv,to_left[ri1]}))*/
                    /*         continue; */
                    /* }else{ */
                    /*     if(!IsUnique(vec<int>{v,to_right[i1],to_right[i2], */
                    /*                 rv,to_left[ri1],to_left[ri2]})) */
                    /*         continue; */
                    /* } */

                    int n1 = D_old.O(i1).isize();
                    int n2 = D_old.O(i2).isize();

                    // get zip
                    int len = 0;
                    while(len<n1 && len<n2){
                        if(D_old.O(i1)[len]==D_old.O(i2)[len]) len++;
                        else break;
                    }

                    // don't allow some full zips
                    if(to_right[i1]==to_right[i2])
                        if((len==n1 && len < n2)||(len==n2 && len < n1))
                            len--;

                    if (len>0){
                        #pragma omp critical
                        {   // PRINT2( D.To(v).size( ), D.From(v).size( ) );
                            Matches.push(make_pair(i1,0), make_pair(i2,0),len); 
                            Matches.push(make_pair(ri1,n1-len),
                                    make_pair(ri2,n2-len),len);
                        }
                    }
                }
            }
        }
        // uniq sort
        ParallelUniqueSort(Matches);
        if(Matches.size() < min_zips)
            break;
        int nmatches = Matches.size( );
        // glue graphs
        Bool META = False;
        vec<vec<quad<int,uint16_t,int,uint16_t>>> dummy_paths;
        GlueGraphs_HP<T,uint16_t>(D_old,dinv_old,Matches,nE,DD,ddinv,META,dummy_paths);
        D = DD;
        dinv = ddinv;
        npass++;
        /*
        cout<<Date( )<<": " << nmatches << " matches, made "
            <<npass<<" zippering passes"<<endl;
        */
    }
}

template void BabyZipper(digraphE<vec<unsigned char>>& D, vec<int>& dinv, 
        unsigned char nE, const int max_pass, const int min_zips);
template void BabyZipper(digraphE<vec<int>>& D, vec<int>& dinv, 
        int nE, const int max_pass, const int min_zips);

// this serves as a good template to do some graph editing ops based on old tech
template<class T>
void RemoveSimpleHangs(digraphE<vec<T>>& D, vec<int>& dinv, vec<int>& dels ){
    #pragma omp parallel for 
    for (int v = 0; v < D.N(); v++){
        if (D.From(v).size()<2) continue;
        // now at a vertex with outs >= 2
        vec<int> delidx;
        for(int j = 0; j < D.From(v).isize(); j++){
            if(D.From(D.From(v)[j]).empty()) // reached a dead-end
                delidx.push_back(j);
        }
        // if all branches are dead-ends, don't delete anything
        if (delidx.size()==D.From(v).size()) continue;
        // otherwise mark for deletion
        #pragma omp critical
        { for(auto dj: delidx)
            dels.push_back(D.IFrom(v)[dj]); } 
    }

    // symmetrize wrt involution
    auto dels2 = dels;
    for(int k = 0; k < dels2.isize(); k++)
        dels.push_back(dinv[dels2[k]]);
    // standard procedure
}

template void RemoveSimpleHangs(digraphE<vec<unsigned char>>& D, vec<int>& dinv,
     vec<int>& dels );

void Validate( const digraphE<vec<uchar>>& D, const vec<int>& dinv )
{    if ( D.E( ) != dinv.isize( ) )
     {    cout << "\nInvolution has wrong size." << endl;
          TracebackThisProcess( );
          Scram(1);    }
     vec<int> bads;
     for ( int e = 0; e < D.E( ); e++ )
     {    for ( int j = 0; j < D.O(e).isize( ); j++ )
          {    if ( D.O(e)[j] > 3 )
               {    cout << "Illegal entry in edge " << e << "." << endl;
                    TracebackThisProcess( );
                    Scram(1);    }    }    }
     #pragma omp parallel for schedule( dynamic, 10000 )
     for ( int e = 0; e < D.E( ); e++ )
     {    if ( dinv[e] < 0 )
          {
               #pragma omp critical
               {    bads.push_back(e);    }    }
          else if ( dinv[dinv[e]] != e )
          {
               #pragma omp critical
               {    bads.push_back(e);    }    }    }
     Sort(bads);
     for ( int bi = 0; bi < bads.isize( ); bi++ )
     {    int e = bads[bi];
          if ( dinv[e] < 0 )
          {    cout << "\nInvolution of edge " << e
                    << " is " << dinv[e] << ", which doesn't make sense." << endl;
               TracebackThisProcess( );
               Scram(1);    }
          if ( dinv[dinv[e]] != e )
          {    cout << "\nInvolution is not an involution." << endl;
               TracebackThisProcess( );
               Scram(1);    }    }
     vec<int> to_left, to_right;
     D.ToLeftParallel(to_left), D.ToRightParallel(to_right);
     #pragma omp parallel for schedule( dynamic, 10000 )
     for ( int v = 0; v < D.N( ); v++ )
     for ( int j1 = 0; j1 < D.To(v).isize( ); j1++ )
     for ( int j2 = 0; j2 < D.From(v).isize( ); j2++ )
     {    int e1 = D.ITo( v, j1 ), e2 = D.IFrom( v, j2 );
          int re1 = dinv[e1], re2 = dinv[e2];
          if ( to_right[re2] != to_left[re1] )
          {
               #pragma omp critical
               {    bads.push_back(v);    }    }    }
     Sort(bads);
     for ( int bi = 0; bi < bads.isize( ); bi++ )
     {    int v = bads[bi];
          for ( int j1 = 0; j1 < D.To(v).isize( ); j1++ )
          for ( int j2 = 0; j2 < D.From(v).isize( ); j2++ )
          {    int d1 = D.ITo( v, j1 ), d2 = D.IFrom( v, j2 );
               int rd1 = dinv[d1], rd2 = dinv[d2];
               if ( to_right[rd2] != to_left[rd1] )
               {    cout << "\n";
                    PRINT4( d1, d2, rd1, rd2 );
                    cout << "edge d1 abuts edge d2, "
                         "so rd2 should abut rd1, but it doesn't" << endl;
                    cout << "actual edges after rd2 are:";
                    int v = to_right[rd2], w = to_left[rd1];
                    for ( int l = 0; l < D.From(v).isize( ); l++ )
                         cout << " " << D.IFrom( v, l );
                    cout << "\nactual edges before rd1 are:";
                    for ( int l = 0; l < D.To(w).isize( ); l++ )
                         cout << " " << D.ITo( w, l );
                    cout << endl;
                    cout << "Involution not consistent with graph structure."
                         << endl;
                    TracebackThisProcess( );
                    Scram(1);    }    }    }    }


void CallGraphEditors(digraphE<vec<unsigned char>>& D, vec<int>& dinv, 
        vec<vec<unsigned char>>& all_closures, const vec<int64_t>& cinv,
        vec<vec<quad<int,uint16_t,int,uint16_t>>>& r_supp, const Bool verbose)
{

     if(verbose) PRINTDEETS("get read paths and index");
     ReadPathVec paths;
     vec<vec<int>> paths_index;
     GetReadPathsAndIndex(r_supp, paths, paths_index, all_closures.size());

     if(verbose) PRINTDEETS("clean-up graph");


     // ============================================================================
     // CODE USING PATHS
     // ============================================================================


     // Delete some very weak bubble edges.

     cout << Date( ) << ": deleting very weak bubble edges" << endl;
     vec<int> dels;
     const int MIN_WIN = 8;
     const int MAX_LOSE = 2;
     const int MAX_EDGE = 5;
     for ( int v = 0; v < D.N( ); v++ )
     {    int s = D.From(v).size( );
          if ( !D.To(v).solo( ) || s <= 1 ) continue;
          int w = D.From(v)[0];
          Bool bad = False;
          for ( int j = 0; j < s; j++ )
          {    if ( D.From(v)[j] != w ) bad = True;
               int d = D.IFrom(v,j);
               if ( D.O(d).isize( ) > MAX_EDGE ) bad = True;    }
          if ( bad || s != D.To(w).isize( ) ) continue;
          if ( v == w || !D.From(w).solo( ) ) continue;
          vec<int> d(s), n(s);
          for ( int j = 0; j < s; j++ ) 
          {    d[j] = D.IFrom(v,j);
               n[j] = paths_index[ d[j] ].size( );    }
          ReverseSortSync( n, d );
          Bool ok = True;
          if ( n[0] < MIN_WIN ) ok = False;
          if ( s >= 2 && n[1] > MAX_LOSE ) ok = False;
          if ( s >= 2 && n[0] < MIN_WIN * n[1] ) ok = False;
          if ( !ok )
          {    cout << "bubble @ " << printSeq(n) << endl;
               continue;    }
          for ( int j = 1; j < s; j++ ) dels.push_back( d[j], dinv[d[j]] );    }
     // [dels]
     UniqueSort(dels);
     cout << Date( ) << ": deleting " << dels.size( ) << " edges" << endl;

     // Kill really weak forks.

     cout << Date( ) << ": deleting really weak forks" << endl;
     vec<int> dels2;
     vec<int> to_left, to_right;
     D.ToLeft(to_left), D.ToRight(to_right);
     for ( int d = 0; d < D.E( ); d++ )
     {    if ( paths_index[d].size( ) > 2 ) continue;
          int v = to_left[d];
          if ( !D.From(v).duo( ) ) continue;
          int f = -1, g = -1;
          for ( int j = 0; j < 2; j++ ) if ( D.IFrom(v,j) != d ) f = D.IFrom(v,j);
          if ( paths_index[f].size( ) == 0 ) continue;
          if ( paths_index[f].size( ) < 6 * paths_index[d].size( ) ) continue;
          dels2.push_back( d, dinv[d] );    }
     UniqueSort(dels2);
     cout << Date( ) << ": deleting " << dels2.size( ) << " edges" << endl;
     dels.append(dels2);
     // [dels]
     UniqueSort(dels);

     // More weak fork deletion.

     /*
     cout << Date( ) << ": new really weak fork deletion" << endl;
     dels2.clear( );
     for ( int v = 0; v < D.N( ); v++ )
     {    if ( !D.From(v).duo( ) ) continue;
          if ( !D.To(v).solo( ) ) continue;
          int n1 = 0, n2 = 0;
          int e = D.ITo(v,0);
          int f1 = D.IFrom(v,0), f2 = D.IFrom(v,1);
          for ( int j = 0; j < paths_index[e].isize( ); j++ )
          {    int id = paths_index[e][j];
               const ReadPath& p = paths[id];
               for ( int l = 0; l < (int) p.size( ) - 1; l++ )
               {    if ( p[l] == e && p[l+1] == f1 ) n1++;
                    if ( p[l] == e && p[l+1] == f2 ) n2++;    }
               // Bad: 2, 6
               // Bad: 1, 6 -- fails at 129974098-129974775
               // 2, 8: good for locus 1, bad for locus 2
               // 1, 8: good for locus 1, bad for locus 2
               // 1, 10: good for locus 1, bad for locus 2:
               // 122384187-122384265: no
               // 122384734-122384812: no
               // 122384968-122385046: no

? print fasta(ref(9,122384107,122384187)) yes
GGTATATACCCAGTAATGGGATGGCTGGGTCAAATGGTATTTCTAGTTCTAGATCCCTGAGGAATCGCCACACTGACTTC
contained in 29 reads

? print fasta(ref(9,122384187,122384265)) no
CACAATGGTTGAACTATTTTACAGTCCCACCAACAGTGTAAAAGTGTTCCTATTTCTCCACATCCTCTCCAGCACCTG
contained in 6 reads

? print fasta(ref(9,122384265,122384345))
TTGTTTCCTGACTTTTTAATGATTGCCATTCTAACTGGTGTGAGATGGTATCTCATAGTGGTTTTGATTTGCATTTCTCT

                                           GATGGTATCTCATAGTGGTTTTGATTTGCATTTCTCT 
                                           ==> 13 reads

                                AACTGGTGTGAGATGGTATCTCATAGTGGTTTTGATTTGCATTTCTCT  
                                ==> 9 reads

   TTTCCTGACTTTTTAATGATTGCCATTCTAACTGGTGTGAGATGGTATCTCATAGTGGTTTTGATTTGCATTTCTCT
==> -24764,-106890

 TGTTTCCTGACTTTTTAATGATTGCCATTCTAACTGGTGTGAGATGGTATCTCATAGTGGTTTTGATTTGCATTTCTCT
==> -24764

TTGTTTCCTGACTTTTTAATGATTGCCATTCTAACTGGTGTGA ==> 67 reads

TTGTTTCCTGACTTTTTAATGATTGCCATTCTAACTGGTGTGAGATGGTATCTCATAG ==> 1560,-13874,-24764

+24765fw vs 8, 0 mismatches/0 indels (of 151), from 0-151 to 122384004-122384155

(perfect match of length 151)

-24764fw vs 8, 0 mismatches/0 indels (of 128), from 0-128 to 122384229-122384357

(perfect match of length 128)

-106890fw vs 8, 0 mismatches/0 indels (of 128), from 0-128 to 122384267-122384395

(perfect match of length 128)

+106891fw vs 8, 0 mismatches/0 indels (of 151), from 0-151 to 122384067-122384218

(perfect match of length 151)

no:
345 (440 --> 327) = CTTC
CACAA
1480 (327 --> 365) = TGGTTGAACTA
283 (365 --> 366) = T
284 (366 --> 367) = TTTACAGTCCCACCAACAGTGTAAA 
1316 (367 --> 556) = A
446 (556 --> 557) = GTGTTCCTATTTCTCC
1736 (557 --> 405) = ACATCCTCTCCAGCACCTG
TTGTTTCCTGACTTTTTAATGAT

1334,1681,1360,1361,344,345,1480,283,284,
1316,446,1736,750,1268,752,90,91,652,1368,1646,653,33,110,111,1264,112,181,1300,1301,1296,1297,1980,1980,1999,1294,709,996,997,998,1467,1468,168,169,170,170,170,170,1342,144,144,145,1697,171,1688,141,1341,1858,1859,473,1349,414,415,416,68,69,70,427,71,1279,72,73,1801,165,166,1353,1523,172,1355,1497,1498,1494,1495,1283,1655,51,51,51,51,51,51,1285,1409,1410,1412,1822,1823,484,481,482,993,605,1239,39,1274,1215,1216,1275,1276,41,1652,1504,748,1554,158,1222,1223,1223,1223,1223,1223,1224,1645,562,563,1254,1143,1733,1885,1793,1256,21,615,616,617,363,1535,521,1202,1609,1417,1345,929,200,930,1543,840,390,391,391,391,391,391,391,1260,549,1258,1257,550,551,1263,1720,734,1352

               if ( n1 <= 1 && n2 >= 10 * n1 && n2 >= 10 )
                    dels2.push_back( f1, dinv[f1] );    }    }
     UniqueSort(dels2);
     cout << Date( ) << ": deleting " << dels2.size( ) << " edges" << endl;
     dels.append(dels2);
     UniqueSort(dels);
     */

     // Eliminate some loops.

     cout << Date( ) << ": eliminate some loops" << endl;
     int nloops = 0;
     for ( int v = 0; v < D.N( ); v++ )
     {    int n = 0;
          for ( int i = 0; i < D.From(v).isize( ); i++ )
               if ( D.From(v)[i] == v ) n++;
          if ( n == 0 ) continue;
          vec<int> loops;
          int x = -1, y = -1, f = -1, g = -1;
          for ( int i = 0; i < D.From(v).isize( ); i++ )
          {    if ( D.From(v)[i] == v ) loops.push_back( D.IFrom(v,i) );
               else
               {    g = D.IFrom(v,i);
                    y = D.From(v)[i];    }    }
          for ( int i = 0; i < D.To(v).isize( ); i++ )
          {    if ( D.To(v)[i] != v )
               {    f = D.ITo(v,i);
                    x = D.To(v)[i];    }    }
          if ( make_pair( dinv[g], dinv[f] ) < make_pair( f, g ) ) continue;
          Bool deleted = False;
          for ( auto d : loops ) if ( BinMember( dels, d ) ) deleted = True;
          if (deleted) 
          {    cout << "\nsee " << n << " loops but some already deleted" << endl;
               continue;    }
          if ( D.From(v).isize( ) != n + 1 || D.To(v).isize( ) != n + 1 ) 
          {    cout << "\nsee " << n << " loops but multiple entry/exit" << endl;
               continue;    }
          if ( BinMember( dels, f ) || BinMember( dels, g ) )
          {    cout << "\nsee " << n << " loops but left or right deleted" << endl;
               continue;    }
          if ( !IsUnique( vec<int>{ x, v, y } ) ) 
          {    cout << "\nsee " << n << " loops but nonunique1" << endl;
               continue;    }
          if ( !IsUnique( vec<int>{ f, g, dinv[f], dinv[g] } ) )
          {    cout << "\nsee " << n << " loops but nonunique2" << endl;
               continue;    }
          vec<int> s1, s2;
          for ( auto id : paths_index[f] ) s1.push_back(id);
          for ( auto id : paths_index[g] ) s2.push_back(id);
          UniqueSort(s1), UniqueSort(s2);
          vec<int> s = Intersection( s1, s2 );
          vec<vec<int>> mids;
          for ( auto id : s )
          {    const ReadPath& p = paths[id];
               for ( int l = 0; l < (int) p.size( ); l++ )
               {    if ( p[l] == f )
                    {    vec<int> m;
                         while(1)
                         {    l++;
                              if ( l == (int) p.size( ) ) break;
                              if ( p[l] == g ) break;
                              else m.push_back( p[l] );    }
                         if ( l < (int) p.size( ) ) mids.push_back(m);    }    }    }
          Sort(mids);
          vec< pair< int, vec<int> > > xc;
          for ( int i = 0; i < mids.isize( ); i++ )
          {    int j = mids.NextDiff(i);
               xc.push( j - i, mids[i] );
               i = j - 1;    }
          ReverseSort(xc);
          cout << "\nsee " << n << " loops\n";
          cout << "have " << s.size( ) << " reads that probably bridge" << endl;
          cout << "have " << xc.size( ) << " bridgers:" << endl;
          for ( int i = 0; i < xc.isize( ); i++ )
          {    cout << "[" << xc[i].first << "] "
                    << printSeq( xc[i].second ) << endl;    }
          if ( xc.empty( ) ) continue;
          int n1 = xc[0].first;
          int n2 = ( xc.size( ) >= 2 ? xc[1].first : 0 );
          if ( n2 <= 2 && n1 >= 10 )
          {    cout << "removing loops" << endl;
               nloops += 2*n;
               const vec<int>& x = xc[0].second;
               for ( auto e : x ) D.OMutable(f).append( D.O(e) );
               // [pather]
               for ( auto e : x ) AppendSupport(r_supp[f],r_supp[e]);
               Bool verbose = False;
               if (verbose)
               {    cout << "fnew = ";
                    for ( int l = 0; l < D.O(f).isize( ); l++ )
                         cout << as_base( D.O(f)[l] );
                    cout << endl;
                    cout << "g = ";
                    for ( int l = 0; l < D.O(g).isize( ); l++ )
                         cout << as_base( D.O(g)[l] );
                    cout << endl;    }
               vec<uchar> y = D.O(f);
               y.ReverseMe( );
               for ( auto& m : y ) m = 3 - m;
               D.OMutable( dinv[f] ) = y;
               // [pather]
               RCSupport(r_supp[f],r_supp[dinv[f]],cinv,all_closures,D.O(f).size());
               for ( auto d : loops ) dels.push_back( d, dinv[d] );    }    }
     cout << Date( ) << ": eliminated " << nloops << " loops" << endl;
     UniqueSort(dels);

     // Flatten some canonical loops.

     cout << Date( ) << ": flattening canonical loops" << endl;
     vec<int> delsx;
     for ( int v = 0; v < D.N( ); v++ )
     {    if ( !D.To(v).duo( ) || !D.From(v).solo( ) ) continue;
          int w = D.From(v)[0];
          if ( v == w || !D.From(w).duo( ) || !D.To(w).solo( ) ) continue;
          int a = D.IFrom(v,0), b = -1, x = -1, y = -1, r = -1, s = -1;
          for ( int j = 0; j < D.From(w).isize( ); j++ )
          {    if ( D.From(w)[j] == v ) b = D.IFrom(w,j);
               else y = D.IFrom(w,j), s = D.From(w)[j];    }
          if ( b < 0 ) continue;
          for ( int j = 0; j < D.To(v).isize( ); j++ )
               if ( D.ITo(v)[j] != b ) x = D.ITo(v,j), r = D.To(v)[j];
          if ( BinMember( dels, a ) || BinMember( dels, b ) ) continue;
          if ( BinMember( dels, x ) || BinMember( dels, y ) ) continue;
          int ra = dinv[a], rb = dinv[b], rx = dinv[x], ry = dinv[y];
          if ( ra < a ) continue;
          if ( !IsUnique( vec<int>{ r, v, w, s } ) ) continue;
          if ( !IsUnique( vec<int>{ a, b, x, y, ra, rb, rx, ry } ) ) continue;
          int bab = 0, bay = 0;
          vec<int> M;
          for ( auto id : paths_index[b] ) 
          {    const ReadPath& p = paths[id];
               for ( int j = 0; j <= (int) p.size( ) - 3; j++ )
               {    if ( p[j] == b && p[j+1] == a && p[j+2] == b ) bab++;
                    if ( p[j] == b && p[j+1] == a && p[j+2] == y ) bay++;    }    }
          PRINT2( bab, bay );
          if ( bab == 0 && bay >= 5 )
          {    D.OMutable(x).append( D.O(a) );
               D.OMutable(x).append( D.O(b) );
               // [pather]
               AppendSupport(r_supp[x],r_supp[a]);
               AppendSupport(r_supp[x],r_supp[b]);
               delsx.push_back( b, dinv[b] );
               D.OMutable(rx) = D.O(x);
               D.OMutable(rx).ReverseMe( );
               for ( auto& e : D.OMutable(rx) ) e = 3 - e;    
               // [pather]
               RCSupport(r_supp[x],r_supp[rx],cinv,all_closures,D.O(x).size());
          }    }
     dels.append(delsx);

     // Clean up.

     D.DeleteEdges(dels);
     /* RemoveUnneededVertices( D, dinv); */
     /* CleanupCore( D, dinv); */
     // [pather]
     DeleteEdges(r_supp,dels);
     RemoveUnneededVertices( D, dinv, r_supp);
     CleanupCore( D, dinv, r_supp);
     Validate( D, dinv );

     // Remove some hanging ends.  Dropping the call to RemoveSimpleHangs increased
     // the number of edges in the final assembly (when tested).

     dels.clear( );
     RemoveSimpleHangs(D,dinv,dels);
     ParallelUniqueSort(dels);
     D.DeleteEdgesParallel(dels);
     DeleteEdges(r_supp,dels);
     RemoveUnneededVertices(D,dinv,r_supp);
     CleanupCore(D,dinv,r_supp);
     /* RemoveUnneededVertices( D, dinv ); // weirdly, significantly affects results */

     // Trim the graph.

     cout << Date( ) << ": trimming" << endl;
     vec<int> lens( D.E( ), 0 ), dfw;
     for ( int e = 0; e < D.E( ); e++ ) lens[e] = D.O(e).size( );
     const int MIN_RATIO = 10;
     const int MAX_KILL = 600;
     // DistancesToEndArr( D, lens, MAX_KILL * MIN_RATIO, True, dfw );
     MaxDistanceToEndArr( D, lens, MAX_KILL * MIN_RATIO, True, dfw );
     dels.clear( );
     D.ToLeft(to_left), D.ToRight(to_right);
     for ( int v = 0; v < D.N( ); v++ )
     {    for ( int j1 = 0; j1 < D.From(v).isize( ); j1++ )
          {    int f1 = D.IFrom(v,j1);
               int n1 = lens[f1] + dfw[ D.From(v)[j1] ];
               if ( n1 > MAX_KILL ) continue;
               for ( int j2 = 0; j2 < D.From(v).isize( ); j2++ )
               {    int f2 = D.IFrom(v,j2);
                    int n2 = lens[f2] + dfw[ D.From(v)[j2] ];
                    if ( n2 < MIN_RATIO * n1 ) continue;
                    vec<int> es = {f1};
                    for ( int i = 0; i < es.isize( ); i++ )
                    {    int w = to_right[ es[i] ];
                         for ( int j = 0; j < D.From(w).isize( ); j++ )
                         {    int e = D.IFrom(w,j);
                              if ( Member( es, e ) ) continue; // prob. can't happen
                              es.push_back(e);    }    }
                    {    for ( int i = 0; i < es.isize( ); i++ )
                              dels.push_back( es[i], dinv[ es[i] ] );    }
                                   }    }    }
     cout << Date( ) << ": deleting " << dels.size( ) << " edges" << endl;
     D.DeleteEdges(dels);
     DeleteEdges(r_supp,dels);
     CleanupCore( D, dinv, r_supp );

     // Kill really weak forks (again).

     /*
     cout << Date( ) << ": deleting really weak forks (again)" << endl;
     GetReadPathsAndIndex(r_supp, paths, paths_index, all_closures.size());
     dels.clear( );
     D.ToLeft(to_left), D.ToRight(to_right);
     for ( int d = 0; d < D.E( ); d++ )
     {    if ( paths_index[d].size( ) > 2 ) continue;
          int v = to_left[d];
          if ( !D.From(v).duo( ) ) continue;
          int f = -1, g = -1;
          for ( int j = 0; j < 2; j++ ) if ( D.IFrom(v,j) != d ) f = D.IFrom(v,j);
          if ( paths_index[f].size( ) == 0 ) continue;
          if ( paths_index[f].size( ) < 6 * paths_index[d].size( ) ) continue;
          dels.push_back( d, dinv[d] );    }
     UniqueSort(dels);
     cout << Date( ) << ": deleting " << dels.size( ) << " edges" << endl;
     D.DeleteEdges(dels);
     DeleteEdges(r_supp,dels);
     CleanupCore( D, dinv, r_supp );
     */


     // ============================================================================
     // CODE NOT USING PATHS AT PRESENT
     // ============================================================================


     // Prezipper. (NOT YET, JUST COUNING.)
     // Looks for simple forks.

     double pclock = WallClockTime( );
     cout << Date( ) << ": prezipper" << endl;
     int zcount = 0;
     for ( int v = 0; v < D.N( ); v++ )
     {    if ( !D.To(v).solo( ) || !D.From(v).duo( ) ) continue;
          int f = D.ITo(v,0);
          int d1 = D.IFrom(v,0), d2 = D.IFrom(v,1);
          if ( D.O(d1).solo( ) || D.O(d2).solo( ) ) continue;
          if ( D.O(d1)[0] != D.O(d2)[0] ) continue;
          int rf = dinv[f], rd1 = dinv[d1], rd2 = dinv[d2];
          if ( !IsUnique( vec<int>{ f, d1, d2, rf, rd1, rd2 } ) ) continue;
          zcount += 2;    }
     PRINT(zcount);

     cout << Date( ) << ": " << TimeSince(pclock) << " used prezippering" << endl;

     // Zipper.  Turning off BabyZipper increased the number of edges in the 
     // final assembly (when tested).
     // 
     // BabyZipper has some verbose logging, can be processed like this:
     // cat whatever | grep size | grep From | less | sort -n -k 4 -n -k 8
     //              | uniq -c | sort -n -r -k1 | less
     // And maybe we shold turn off the verbose logging.

     cout << Date( ) << ": start baby zipper" << endl;
     double zclock = WallClockTime( );
     BabyZipper(D,dinv,(unsigned char)4,10,10);
     cout << Date( ) << ": " << TimeSince(zclock) << " used in baby zipper" << endl;

     // Assay for multiloops.

     cout << "\nmultiloops before removing redundant loops\n";
     for ( int v = 0; v < D.N( ); v++ )
     {    int n = 0;
          for ( int i = 0; i < D.From(v).isize( ); i++ )
               if ( D.From(v)[i] == v ) n++;
          if ( n >= 2 ) PRINT(n);    }
     cout << endl;

     // Assay for multiloops.

     cout << "\nmultiloops after removing redundant loops\n";
     for ( int v = 0; v < D.N( ); v++ )
     {    int n = 0;
          for ( int i = 0; i < D.From(v).isize( ); i++ )
               if ( D.From(v)[i] == v ) n++;
          if ( n >= 2 ) PRINT(n);    }
          
     // Remove small components.

     dels.clear( );
     vec<vec<int>> comp;
     D.ComponentsE(comp);
     for ( int i = 0; i < comp.isize( ); i++ )
     {    int n = 0;
          for ( auto e : comp[i] ) n += D.O(e).size( );
          if ( n < 300 ) dels.append( comp[i] );    }

     cout << Date( ) << ": deleting " << dels.size( )
          << " edges in small components" << endl;
     D.DeleteEdges(dels);
     RemoveUnneededVertices( (digraphE<vec<uchar>>&) D, dinv );
     CleanupCore( D, dinv );
     Validate( D, dinv );

     // Disengage some terminal loops.

     D.ToLeft(to_left), D.ToRight(to_right);
     for ( int v = 0; v < D.N( ); v++ )
     {    if ( D.From(v).solo( ) && D.From(v)[0] == v && D.To(v).size( ) > 1 )
          {    int d = D.IFrom(v,0);
               int rd = dinv[d];
               if ( d == rd ) continue;
               int N = D.N( );
               D.AddVertices(2);
               D.GiveEdgeNewToVx( d, v, N );
               int rv = to_left[rd];
               D.GiveEdgeNewFromVx( rd, rv, N+1 );    }    }
     RemoveUnneededVertices( (digraphE<vec<uchar>>&) D, dinv );
     CleanupCore( D, dinv );
     Validate( D, dinv );

     // Turn canonical loops into single loops.

     cout << Date( ) << ": looping canonical loops" << endl;
     dels.clear( );
     D.ToLeft(to_left), D.ToRight(to_right);
     for ( int v = 0; v < D.N( ); v++ )
     {    if ( !D.To(v).duo( ) || !D.From(v).solo( ) ) continue;
          int w = D.From(v)[0];
          if ( v == w || !D.From(w).duo( ) || !D.To(w).solo( ) ) continue;
          int a = D.IFrom(v,0), b = -1, x = -1, y = -1, r = -1, s = -1;
          for ( int j = 0; j < D.From(w).isize( ); j++ )
          {    if ( D.From(w)[j] == v ) b = D.IFrom(w,j);
               else y = D.IFrom(w,j), s = D.From(w)[j];    }
          if ( b < 0 ) continue;
          for ( int j = 0; j < D.To(v).isize( ); j++ )
               if ( D.ITo(v)[j] != b ) x = D.ITo(v,j), r = D.To(v)[j];
          int ra = dinv[a], rb = dinv[b], rx = dinv[x], ry = dinv[y];
          if ( ra < a ) continue;
          if ( !IsUnique( vec<int>{ r, v, w, s } ) ) continue;
          if ( !IsUnique( vec<int>{ a, b, x, y, ra, rb, rx, ry } ) ) continue;
          dels.push_back( b, dinv[b] );
          vec<uchar> ba = D.O(b);
          ba.append( D.O(a) );
          int E = D.E( );
          D.AddEdge( w, w, ba );
          int rw = to_right[ry];
          ba.ReverseMe( );
          for ( auto& e : ba ) e = 3 - e;
          D.AddEdge( rw, rw, ba );
          dinv.push_back( E+1, E );    }
     D.DeleteEdges(dels);
     RemoveUnneededVertices( (digraphE<vec<uchar>>&) D, dinv );
     CleanupCore( D, dinv );
     Validate( D, dinv );

     // Remove some more unneeded vertices.

     cout << Date( ) << ": removing more unneeded vertices" << endl;
     dels.clear( );
     set<int> delss;
     D.ToLeft(to_left), D.ToRight(to_right);
     for ( int v = 0; v < D.N( ); v++ )
     {    if ( !D.To(v).solo( ) || !D.From(v).solo( ) ) continue;
          int f = D.ITo(v,0), g = D.IFrom(v,0);
          int rf = dinv[f], rg = dinv[g];
          if ( !IsUnique( vec<int>{ f, g, rf, rg } ) ) continue;
          if ( Member( delss, f ) || Member( delss, g ) ) continue;
          if ( Member( delss, rf ) || Member( delss, rg ) ) continue;
          dels.push_back( f, g, rf, rg );
          delss.insert(f), delss.insert(g), delss.insert(rf), delss.insert(rg);
          int x = D.To(v)[0], y = D.From(v)[0];
          vec<uchar> h = D.O(f);
          h.append( D.O(g) );
          int E = D.E( );
          D.AddEdge( x, y, h );
          int rx = to_left[ dinv[g] ], ry = to_right[ dinv[f] ];
          h.ReverseMe( );
          for ( auto& e : h ) e = 3 - e;
          D.AddEdge( rx, ry, h );
          dinv.push_back( E+1, E );    }
     D.DeleteEdges(dels);
     CleanupCore( D, dinv );
     Validate( D, dinv );

     // Remove some redundant loops.  This looks for vertices that already have
     // A, G, T and T loops, and kills other loops there.

     cout << Date( ) << ": deleting redundant loops" << endl;
     dels.clear( );
     for ( int v = 0; v < D.N( ); v++ )
     {    int nloops = 0;
          for ( int i = 0; i < D.From(v).isize( ); i++ )
               if ( D.From(v)[i] == v ) nloops++;
          if ( nloops <= 1 ) continue;
          vec<pair<vec<uchar>,int>> loops;
          for ( int i = 0; i < D.From(v).isize( ); i++ )
          {    if ( D.From(v)[i] == v ) 
               {    int e = D.IFrom(v,i);
                    loops.push( D.O(e), e );    }    }
          Sort(loops);
          vec<Bool> have( 4, False );
          for ( auto& x : loops ) 
               if ( x.first.size( ) == 1 ) have[ x.first[0] ] = True;
          if ( Sum(have) < 4 ) continue;
          for ( int i = 0; i < loops.isize( ); i++ )
          {    if ( loops[i].first.size( ) > 1 )
               {    int d = loops[i].second;
                    if ( dinv[d] >= d ) 
                    {    dels.push_back( d, dinv[d] );
                         cout << "deleting loop from cluster of " << nloops << endl;
                              }    }    }    }
     cout << Date( ) << ": deleting " << dels.size( ) << " redundant loops" << endl;
     D.DeleteEdges(dels);
     RemoveUnneededVertices( (digraphE<vec<uchar>>&) D, dinv );
     CleanupCore( D, dinv );

     // Actually kill the ACGT loops.  They make evaluation essentially impossible.

     cout << Date( ) << ": deleting ACGT loops" << endl;
     dels.clear( );
     for ( int v = 0; v < D.N( ); v++ )
     {    int nloops = 0;
          for ( int i = 0; i < D.From(v).isize( ); i++ )
               if ( D.From(v)[i] == v ) nloops++;
          if ( nloops < 4 ) continue;
          vec<pair<vec<uchar>,int>> loops;
          for ( int i = 0; i < D.From(v).isize( ); i++ )
          {    if ( D.From(v)[i] == v ) 
               {    int e = D.IFrom(v,i);
                    loops.push( D.O(e), e );    }    }
          Sort(loops);
          vec<Bool> have( 4, False );
          for ( auto& x : loops ) 
               if ( x.first.size( ) == 1 ) have[ x.first[0] ] = True;
          if ( Sum(have) < 4 ) continue;
          for ( int i = 0; i < loops.isize( ); i++ )
          {    if ( loops[i].first.size( ) == 1 ) 
                    dels.push_back( loops[i].second );    }    }
     cout << Date( ) << ": deleting " << dels.size( ) << " ACGT loops" << endl;
     D.DeleteEdges(dels);
     RemoveUnneededVertices( (digraphE<vec<uchar>>&) D, dinv );
     CleanupCore( D, dinv );    

     // Kill some inversion squares.

     D.ToLeft(to_left), D.ToRight(to_right);
     dels.clear( );
     for ( int v = 0; v < D.N( ); v++ )
     {    if ( !D.To(v).solo( ) || !D.From(v).duo( ) ) continue;
          int j = 0;
          if ( !D.From( D.From(v)[j] ).duo( ) ) j = 1 - j;
          int f = D.IFrom(v)[j], g = D.IFrom(v)[1-j];
          int w = D.From(v)[j], x = D.From(v)[1-j];
          if ( !D.From(w).duo( ) || !D.From(x).solo( ) ) continue;
          if ( !D.To(x).duo( ) ) continue;
          int y = D.From(x)[0], rf = D.IFrom(x,0);
          if ( !D.To(y).duo( ) || !D.From(y).solo( ) ) continue;
          int rg = -1;
          for ( int l = 0; l < 2; l++ ) if ( D.From(w)[l] == y ) rg = D.IFrom(w,l);
          if ( rf != dinv[f] || rg != dinv[g] ) continue;
          dels.push_back( g, dinv[g] );    }
     D.DeleteEdges(dels);
     RemoveUnneededVertices( (digraphE<vec<uchar>>&) D, dinv );
     CleanupCore( D, dinv );

     // More zippering.

     BabyZipper(D,dinv,(unsigned char)4,100,1);
     RemoveUnneededVertices( (digraphE<vec<uchar>>&) D, dinv );
     CleanupCore( D, dinv );    }
