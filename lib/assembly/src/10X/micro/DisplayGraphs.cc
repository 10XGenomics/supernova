// Copyright 2016, 10x Genomics, Inc.

#include "10X/micro/DisplayGraphs.h"

template <class F>
void PrintGraphPaths(const digraphE<F>& ghb, const vec<vec<int>>& gold,
        String fname, vec<String> color)
{
    vec<int> meta;
    PrintGraphPaths(ghb, gold, meta, fname, color);
}

template void PrintGraphPaths(const digraphE<vec<uchar>>& ghb, const vec<vec<int>>& gold,
        String fname, vec<String> color);
template void PrintGraphPaths(const digraphE<basevector>& ghb, const vec<vec<int>>& gold,
        String fname, vec<String> color);
template void PrintGraphPaths(const digraphE<vec<int>>& ghb, const vec<vec<int>>& gold,
        String fname, vec<String> color);

template <class F>
void PrintGraphPaths(const digraphE<F>& ghb, const vec<vec<int>>& gold,
        const vec<int>& meta, String fname, vec<String> color)
{
      Bool use_meta = False;
      if(meta.size()==ghb.E())
          use_meta = True;
      vec<vec<String>> legends;
      vec<String> legend_colors;
      String layout;
      int N = ghb.N( );
      vec<vec<String>> edge_labels(N), edge_colors(N);
      vec<vec<String>> edge_attrs(N), edge_attrs2(N);
      vec<String> vertex_colors(N);
      for ( int v = 0; v < N; v++ )
      {    edge_labels[v].resize( ghb.From(v).size( ) );
           edge_colors[v].resize( ghb.From(v).size( ) );
           edge_attrs[v].resize( ghb.From(v).size( ) );
           edge_attrs2[v].resize( ghb.From(v).size( ) );    }
      for ( int v = 0; v < N; v++ )
      for ( int j = 0; j < (int) ghb.From(v).size( ); j++ )
      {    int f = ghb.IFrom(v,j);
           edge_labels[v][j] = "#" + ToString(f) + " * " + ToString( ghb.O(f).size( ) ) 
                               + "b" + ((use_meta)? " * " + ToString(meta[f]) + "r" : "");
           for(int k = 0; k < gold.isize(); k++){
               if ( BinMember( gold[k], f ) ){
                   edge_colors[v][j] = color[k];    
               }
           }
      }
      Ofstream( dout, fname );
      ghb.DOT( dout, edge_labels, vertex_colors, edge_colors, 
           edge_attrs, edge_attrs2, legends, legend_colors, layout );
}

template void PrintGraphPaths(const digraphE<vec<uchar>>& ghb, const vec<vec<int>>& gold,
        const vec<int>& meta, String fname, vec<String> color);
template void PrintGraphPaths(const digraphE<basevector>& ghb, const vec<vec<int>>& gold,
        const vec<int>& meta, String fname, vec<String> color);
template void PrintGraphPaths(const digraphE<vec<int>>& ghb, const vec<vec<int>>& gold,
        const vec<int>& meta, String fname, vec<String> color);

template <class F>
void DisplayGraphPaths(const digraphE<F>& dhb, const ReadPathVec& dpaths,
                 const vec<vec<int>>& dpaths_index, const vec<int> IDs, Bool verbose, String text, String fname)
{
     vec<vec<int>> Gold;
     vec<String> colors;
     Gold.resize(1);
     int placed = 0;
     for(auto d : IDs){
         if(dpaths[d].size())
             placed++;
         for(auto e: dpaths[d])
             Gold.back().push_back(e);
     }
     Sort(Gold.back());
     if(verbose)
         cout << Date() << ": placed " << (100.0*placed)/IDs.size() << "% of "<< text << endl;
     colors.push_back("red");
     vec<int> meta(dpaths_index.size());
     for(int d = 0; d < dpaths_index.isize(); d++)
         meta[d] = dpaths_index[d].size();
     PrintGraphPaths(dhb,Gold,meta,fname,colors);
}

template void DisplayGraphPaths(const digraphE<vec<uchar>>& dhb, const ReadPathVec& dpaths,
                 const vec<vec<int>>& dpaths_index, const vec<int> IDs, Bool verbose, String text, String fname);
template void DisplayGraphPaths(const digraphE<basevector>& dhb, const ReadPathVec& dpaths,
                 const vec<vec<int>>& dpaths_index, const vec<int> IDs, Bool verbose, String text, String fname);
template void DisplayGraphPaths(const digraphE<vec<int>>& dhb, const ReadPathVec& dpaths,
                 const vec<vec<int>>& dpaths_index, const vec<int> IDs, Bool verbose, String text, String fname);

template <class F>
void DisplayGraphPaths(const digraphE<F>& dhb, const vec<vec<int>>& Gold, const vec<String>& colors,
                 const vec<vec<int>>& dpaths_index, Bool verbose, String fname)
{
     vec<int> meta(dpaths_index.size());
     for(int d = 0; d < dpaths_index.isize(); d++)
         meta[d] = dpaths_index[d].size();
     PrintGraphPaths(dhb,Gold,meta,fname,colors);
}

template void DisplayGraphPaths(const digraphE<vec<uchar>>& dhb, const vec<vec<int>>& Gold, const vec<String>& colors,
                 const vec<vec<int>>& dpaths_index, Bool verbose, String fname);
template void DisplayGraphPaths(const digraphE<basevector>& dhb, const vec<vec<int>>& Gold, const vec<String>& colors,
                 const vec<vec<int>>& dpaths_index, Bool verbose, String fname);
template void DisplayGraphPaths(const digraphE<vec<int>>& dhb, const vec<vec<int>>& Gold, const vec<String>& colors,
                 const vec<vec<int>>& dpaths_index, Bool verbose, String fname);

void DisplayRefGraphPath(const digraphE<bvec>& dhb, const vec<bvec>& ref, 
        const vec<vec<int>>& dpaths_index, Bool verbose, String fname)
{
     vec<vec<int>> gpaths;
     vec<vec<int>> Gold;
     vec<String> colors;
     Gold.resize(1);
     colors.push_back("red");
     for(auto r: ref){
         GetPaths(dhb,r,gpaths);
         Gold.back().append(Contents(gpaths));
     }
     ParallelUniqueSort(Gold.back());
     ForceAssertEq(dhb.E(),dpaths_index.size());
     DisplayGraphPaths(dhb,Gold,colors,dpaths_index,verbose,fname);
}

template <class E, class F>
void DisplaySCC(const digraphE<E>& D, const digraphE<F>& dhb, 
        const vec<vec<int>>& dpaths_index, Bool verbose, String fname)
{
     vec<vec<int>> Gold;
     vec<String> colors;
     vec<int> to_right,to_left;
     D.ToRight(to_right);
     D.ToLeft(to_left);
     vec<int> loopEs, loopVs;
     for(int d = 0; d < D.E(); d++)
     {
         if(to_right[d]==to_left[d]){
             loopEs.push_back(d);
             loopVs.push_back(to_right[d]);
         }
     }
     UniqueSortSync(loopVs,loopEs);// lost double loops here
     vec<vec<int>> SCC;
     D.StronglyConnectedComponents(SCC);
     vec<int> scce;
     vec<int> sccse;
     for(int s = 0; s < SCC.isize(); s++){
         if (SCC[s].size()==1){ // avoid singletons unless loops
             int ps = BinPosition(loopVs, SCC[s][0]);
             if(ps>=0){
                 vec<int> df = D.IFrom(loopVs[ps]);
                 vec<int> dt = D.ITo(loopVs[ps]);
                 Sort(df); Sort(dt);
                 vec<int> inter;
                 Intersection(df,dt,inter);
                 sccse.append(inter);
             }
             continue;
         }
         Sort(SCC[s]);
         // get all edges only inside SCC
         vec<pair<int,int>> ew;
         for(int i = 0; i < SCC[s].isize(); i++){
             for(int j = 0; j < D.From(SCC[s][i]).isize(); j++)
                 ew.push(D.From(SCC[s][i])[j],D.IFrom(SCC[s][i],j));
             for(int j = 0; j < D.To(SCC[s][i]).isize(); j++)
                 ew.push(D.To(SCC[s][i])[j],D.ITo(SCC[s][i],j));
         }
         UniqueSort(ew);
         for(auto pr : ew)
             if(BinPosition(SCC[s],pr.first)>=0)
                 scce.push_back(pr.second);
     }
     UniqueSort(scce);
     UniqueSort(sccse);
     vec<int> scccom;
     Intersection(scce,sccse,scccom);
     Sort(scccom);
     Gold.clear();
     Gold.push_back(sccse);
     Gold.push_back(scce);
     Gold.push_back(scccom);
     colors.clear(); colors.push_back("blue"); colors.push_back("red"); colors.push_back("yellow");
     DisplayGraphPaths(dhb, Gold, colors, dpaths_index, verbose, fname);
}

template void DisplaySCC(const digraphE<vec<int>>& D, const digraphE<basevector>& dhb, 
        const vec<vec<int>>& dpaths_index, Bool verbose, String fname);
template void DisplaySCC(const digraphE<int>& D, const digraphE<vec<int>>& dhb, 
        const vec<vec<int>>& dpaths_index, Bool verbose, String fname);
