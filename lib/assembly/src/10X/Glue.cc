// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

#include "Basevector.h"
#include "CoreTools.h"
#include "Set.h"
#include "paths/HyperBasevector.h"
// #include "paths/long/large/GapToyTools.h"
#include "10X/mergers/ExtremeAsm.h"
#include "10X/Glue.h"

Bool Glue( VirtualMasterVec<basevector>& bases, const vec<int64_t>& ids,
     const String& glue_rules, digraphE<basevector>& ghb, vec<int>& ginv )
{
     // tokenize glue rules with K, controls
     vec<String> tokens, controls;
     int K = -1;
     Tokenize(glue_rules,',',tokens); 
     if(tokens.size())
         if(tokens[0].IsInt()){
             K = tokens[0].Int();
             for(int i = 1; i < tokens.isize(); i++)
                 controls.push_back(tokens[i]);
         }

     // if not a bad rule
     if ( K >= 8 ) 
     {
        ExtremeAsm(bases,ids, ghb, ginv, K, controls, False, "");
        return True;
     }
     else{
         cout<<"Please input correct glue rules in place of "<<glue_rules
             <<", ensure that K >= 8"<<endl;
         return False;    
     }

     return False;    
}



/************************************************************
 *             ALIGNERS                                     *
 ***********************************************************/

// default output top 4 matches of repeats >= 8bases.
void SeqAlignments(String& x1, String& x2, const int MIN_WT, 
        const int top, ostringstream& out){ 

    // controls
    std::map<char,int> weight;
    weight.insert(make_pair('A',1));
    weight.insert(make_pair('T',1));
    weight.insert(make_pair('G',1));
    weight.insert(make_pair('C',1));

    // Dynamic program
    vec<vec<int>> matrix;
    matrix.resize(x1.size()+1);
    for(int i1 = 0; i1<x1.isize()+1; i1++)
        matrix[i1].resize(x2.size()+1,0);

    for(int i1 = 0; i1 <x1.isize(); i1++){
        for(int i2 = 0; i2<x2.isize(); i2++){
            if(x1[i1]==x2[i2])
                matrix[i1+1][i2+1] = matrix[i1][i2] + weight[x1[i1]];
        }
    }

    // gather longest common substrings that are acceptable 
    // (min kmer, max extreme content) 

    vec<triple<int,int,int>> Suffs;
    vec<int> lens;
    for(int i1 = 0; i1 <x1.isize(); i1++){
        for(int i2 = 0; i2<x2.isize(); i2++){
            if(matrix[i1+1][i2+1]>0){
                Suffs.push_back(make_triple(matrix[i1+1][i2+1],i1,i2));
                matrix[i1+1][i2+1] = 0; // mark as read
                int k = 1;
                while(i1+k < x1.isize() && 
                        i2+k< x2.isize() && 
                        matrix[i1+1+k][i2+1+k]>0){
                    Suffs.back() = make_triple(matrix[i1+1+k][i2+1+k],i1+k,i2+k);
                    matrix[i1+1+k][i2+1+k] = 0; // mark as read
                    k++;
                }
                if (Suffs.back().first < MIN_WT) // lenient here!!!!
                    Suffs.pop_back(); // prune by length
                else
                    lens.push_back(k);
            }
        }
    }

    // sort by kmer size of match
    ReverseSortSync(Suffs,lens);

    if(Suffs.size()){
    // pretty print at most top 4 overlaps
        int num = 0;
        // skip the first as it is useless
        for(int i = 1; i < Min(top+1,Suffs.isize()); i++){
            auto& t1 = Suffs[i];
            num++;
            out<<"Alignment #"<<num<<" : "<<lens[i]<<" bases"<<endl;
            olap O(t1.second-lens[i]+1,t1.third-lens[i]+1,lens[i]);
            for(int kr = 0; kr<O.start1; kr++)
                out << x1[kr];
            for(int r = 0; r<O.len; r++)
                out << "\033[1;31m"<<x1[r+O.start1]<<"\033[0m";
            for(int kr = O.start1+O.len; kr<x1.isize(); kr++)
                out << x1[kr];
            out<<endl;
            for(int kr = 0; kr<O.start2; kr++)
                out << x2[kr];
            for(int r = 0; r<O.len; r++)
                out << "\033[1;31m"<<x2[r+O.start2]<<"\033[0m";
            for(int kr = O.start2+O.len; kr<x2.isize(); kr++)
                out << x2[kr];
            if(i==top || i==Suffs.isize()-1)
                out<<endl;
            else
                out<<endl<<endl;
        }
    }else
        out<<"No large duplications"<<endl<<endl;
}

void SeqSelfAlignments(String& x1,const int MIN_WT, 
        const int top, ostringstream& out){
    String x2 = x1;
    SeqAlignments(x1, x2, MIN_WT, top, out);
}

void GetPaths( const digraphE<basevector>& ghb, const basevector& b,
     vec<vec<int>>& p )
{    p.clear( );
     vec<int> gto_right;
     ghb.ToRight(gto_right);
     String m = b.ToString( );
     vec<String> M( ghb.E( ) );
     #pragma omp parallel for
     for ( int e = 0; e < ghb.E( ); e++ ) M[e] = ghb.O(e).ToString( );
     vec<pair<int,int>> epos;
     #pragma omp parallel for schedule(dynamic,1)
     for ( int e = 0; e < ghb.E( ); e++ )
     {    const basevector& g = ghb.O(e);
          for ( int pos = 0; pos < g.isize( ); pos++ )
          {    Bool mismatch = False;
               for ( int j = 0; ; j++ )
               {    if ( pos+j == g.isize( ) || j == b.isize( ) ) break;
                    if ( g[pos+j] != b[j] )
                    {    mismatch = True;
                         break;    }    }
               if (mismatch) continue;
               #pragma omp critical
               {    epos.push( e, pos );    }    }    }
     Sort(epos);
     #pragma omp parallel for schedule(dynamic,1)
     for ( int u = 0; u < epos.isize( ); u++ )
     {    int e = epos[u].first, pos = epos[u].second;
          const basevector& g = ghb.O(e);
          vec<int> lens;
          String C;
          String s;
     
          // COMPLETELY STUPID IMPLEMENTATION BELOW.

          vec<vec<int>> paths;
          vec<String> PATHS;
          paths.push_back( {e} );
          PATHS.push_back( M[e] );
          lens.clear( );
          lens.push_back( g.isize( ) - pos );
          if ( lens.back( ) >= b.isize( ) )
          {    String c;
               for ( auto e : paths.back( ) ) c.append( M[e] );
               if ( c.Contains(m) ) 
               {
                    #pragma omp critical
                    {    p.push_back( paths.back( ) );    }    }    }
          s.clear( );
          for ( int i = 0; i < paths.isize( ); i++ )
          {    if ( lens[i] >= b.isize( ) ) continue;
               int v = gto_right[ paths[i].back( ) ];
               for ( int j = 0; j < (int) ghb.From(v).size( ); j++ )
               {    vec<int> z = paths[i];
                    int f = ghb.IFrom(v,j);
                    z.push_back(f);
                    int len = lens[i] + ghb.O(f).size( );
                    C = PATHS[i];
                    C.append( M[f] );
                    s = m;
                    if ( len < s.isize( ) ) s.resize(len);
                    if ( !C.Contains( s, pos ) ) continue;
                    if ( len >= b.isize( ) )
                    {
                         #pragma omp critical
                         {    p.push_back(z);    }    }
                    paths.push_back(z), 
                    PATHS.push_back(C);
                    lens.push_back(len);    }    }    }
     Sort(p);    }

Bool IsPath( const digraphE<basevector>& ghb, const basevector& b )
{    vec<Bool> used;
     ghb.Used(used);
     vec<int> to_right;
     ghb.ToRight(to_right);
     vec<pair<int,int>> vpos;
     Bool yes = False;
     #pragma omp parallel for schedule(dynamic,1)
     for ( int e = 0; e < ghb.E( ); e++ )
     {    if ( !used[e] ) continue;
          const basevector& g = ghb.O(e);
          for ( int pos = 0; pos < g.isize( ); pos++ )
          {    
               // See if starting at position pos on g, g matches the beginning of b.

               Bool mismatch = False;
               for ( int j = 0; ; j++ )
               {    if ( j == b.isize( ) ) 
                    {    yes = True;
                         break;    }
                    if ( pos+j == g.isize( ) ) break;
                    if ( g[pos+j] != b[j] )
                    {    mismatch = True;
                         break;    }    }
               if (mismatch) continue;
               #pragma omp critical
               {    vpos.push( to_right[e], g.isize( ) - pos );    }    }    }
     if (yes) return True;
     Sort(vpos);
     set<pair<int,int>> done;
     for ( auto& x : vpos ) done.insert(x);
     while( vpos.nonempty( ) )
     {    int v = vpos.back( ).first, p = vpos.back( ).second;
          vpos.pop_back( );
          for ( int j = 0; j < ghb.From(v).isize( ); j++ )
          {    int v2 = ghb.From(v)[j], i;
               const basevector& g = ghb.OFrom(v,j);
               if ( g.empty( ) ) continue; // GAP CODE
               for ( i = 0; i < g.isize( ); i++ )
               {    if ( p+i == b.isize( ) ) return True;
                    if ( b[p+i] != g[i] ) break;    }
               if ( i == g.isize( ) )
               {    int p2 = p + g.isize( );
                    if ( Member( done, make_pair( v2, p2 ) ) ) continue;
                    done.insert( make_pair( v2, p2 ) );
                    vpos.push( v2, p2 );    }    }    }
     return False;    }

void GetPathCov( const digraphE<basevector>& ghb, const basevector& b,
     vec<int>& cov )
{    cov.clear( );
     vec<Bool> used;
     ghb.Used(used);
     vec<int> to_left, to_right;
     ghb.ToLeft(to_left), ghb.ToRight(to_right);
     vec<triple<int,int,int>> vpos;
     Bool yes = False;
     #pragma omp parallel for schedule(dynamic,1)
     for ( int e = 0; e < ghb.E( ); e++ )
     {    if ( !used[e] ) continue;
          const basevector& g = ghb.O(e);
          for ( int pos = 0; pos < g.isize( ); pos++ )
          {    Bool mismatch = False;
               for ( int j = 0; ; j++ )
               {    if ( j == b.isize( ) ) 
                    {    
                         #pragma omp critical
                         {    cov.push_back(e);    }
                         goto next;    }
                    if ( pos+j == g.isize( ) ) break;
                    if ( g[pos+j] != b[j] )
                    {    mismatch = True;
                         break;    }    }
               if (mismatch) continue;
               #pragma omp critical
               {    vpos.push( to_right[e], g.isize( ) - pos, e );    }    }
          next: continue;    }
     Sort(vpos);
     set<triple<int,int,int>> done;
     for ( auto& x : vpos ) done.insert(x);
     set<pair<int,int>> eright;
     vec<pair<int,int>> erightv;
     vec<vec<int>> pos( ghb.E( ) );
     while( vpos.nonempty( ) )
     {    int v = vpos.back( ).first, p = vpos.back( ).second;
          int f = vpos.back( ).third;
          vpos.pop_back( );
          for ( int j = 0; j < ghb.From(v).isize( ); j++ )
          {    int v2 = ghb.From(v)[j], i;
               int e = ghb.IFrom(v,j);
               const basevector& g = ghb.O(e);
               for ( i = 0; i < g.isize( ); i++ )
               {    if ( p+i == b.isize( ) ) 
                    {    eright.insert( make_pair(e,p) );
                         erightv.push(e,p);
                         cov.push_back( e, f );
                         break;    }
                    if ( b[p+i] != g[i] ) break;    }
               if ( i == g.isize( ) )
               {    int p2 = p + g.isize( );
                    if ( Member( done, make_triple( v2, p2, f ) ) ) continue;
                    done.insert( make_triple( v2, p2, f ) );
                    vpos.push( v2, p2, f );    
                    pos[e].push_back(p);    }    }    }
     UniqueSort(erightv);
     for ( int i = 0; i < erightv.isize( ); i++ )
     {    int e = erightv[i].first, pe = erightv[i].second;
          int v = to_left[e];
          for ( int j = 0; j < ghb.To(v).isize( ); j++ )
          {    int f = ghb.ITo(v,j);
               for ( int k = 0; k < pos[f].isize( ); k++ )
               {    int pf = pos[f][k];
                    if ( pe - pf != ghb.O(f).isize( ) ) continue;
                    if ( Member( eright, make_pair(f,pf) ) ) continue;
                    eright.insert( make_pair(f,pf) );
                    erightv.push(f,pf);
                    cov.push_back(f);    }    }    }
     UniqueSort(cov);     }
