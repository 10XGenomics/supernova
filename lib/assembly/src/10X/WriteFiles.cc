// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"
#include "Intvector.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "paths/long/large/GapToyTools.h"
#include "10X/DfTools.h"
#include "10X/WriteFiles.h"
#include "10X/astats/GenomeAlign.h"

void WriteAssemblyFiles( const HyperBasevector& hb, const vec<int>& inv,
     ReadPathVecX& paths, VecULongVec& paths_index, const vec<int64_t>& bci, 
     const Bool ALIGN, const vecbasevector& genome, const String& dir,
     MasterVec< SerfVec<triple<int,int,int> > >& alignsb ) // RPVX
{
     cout << Date( ) << ": writing files" << endl;
     cout << Date( ) << ": hb has checksum " << hb.CheckSum( ) << endl;
     Mkdir777(dir);
     if ( paths.size() > 0 ) {
          cout << Date() << ": writing paths" << endl;
          cout << Date() << ": Before removing pathsX, mem = " 
               << MemUsageGBString() << endl;
          paths.WriteAll( dir + "/a.pathsX" );
          Destroy(paths); //RPVX
          cout << Date() << ": Destroyed pathsX, mem = " 
               << MemUsageGBString() << endl;
     }
     Echo( ToString( hb.K( ) ), dir + "/a.k" );
     BinaryWriter::writeFile( dir + "/a.hbv", hb );
     vec<int> to_left, to_right;
     hb.ToLeft(to_left), hb.ToRight(to_right);
     BinaryWriter::writeFile( dir + "/a.to_left", to_left );
     BinaryWriter::writeFile( dir + "/a.to_right", to_right );
     BinaryWriter::writeFile( dir + "/a.inv", inv );
     HyperBasevectorX hbx(hb);
     BinaryWriter::writeFile( dir + "/a.hbx", hbx );
     {
          // Write edges and kmers.

          cout << Date( ) << ": making and writing edges" << endl;
          {    vecbvec edges( hb.Edges( ).begin( ), hb.Edges( ).end( ) );
               edges.WriteAll( dir + "/a.fastb" );    }
          {    vec<int> kmers( hb.E( ) );
               for ( int e = 0; e < hb.E( ); e++ )
                    kmers[e] = hb.Kmers(e);
               BinaryWriter::writeFile( dir + "/a.kmers", kmers );    }

          // Align to genome.

          if (ALIGN)
          {    const int align_genome_K = 80;
               GenomeAlign<align_genome_K>( hbx, inv, genome, alignsb );
               alignsb.WriteAll( dir + "/a.alignsb" );    }
          else Remove( dir + "/a.alignsb" );    }

     cout << Date() << ": reading pathsX" << endl;
     if (!IsRegularFile( dir + "/a.pathsX")){
        cout << "a.pathsX does not exist!\n\n"<<endl;
        Scram(1);
     }
     paths.ReadAll( dir + "/a.pathsX");
}


void WriteAssemblyFiles( const HyperBasevector& hb, const vec<int>& inv,
     ReadPathVec& paths, VecULongVec& paths_index, const vec<int64_t>& bci, 
     const Bool ALIGN, const vecbasevector& genome, const String& dir,
     MasterVec< SerfVec<triple<int,int,int> > >& alignsb ) // RPVX
{
     cout << Date( ) << ": writing files" << endl;
     cout << Date( ) << ": hb has checksum " << hb.CheckSum( ) << endl;
     Mkdir777(dir);
     if ( paths.size() > 0 ) {
          cout << Date() << ": writing paths" << endl;
          paths.WriteAll( dir + "/a.paths" );
          paths.resize(0);
     }
     Echo( ToString( hb.K( ) ), dir + "/a.k" );
     BinaryWriter::writeFile( dir + "/a.hbv", hb );
     vec<int> to_left, to_right;
     hb.ToLeft(to_left), hb.ToRight(to_right);
     BinaryWriter::writeFile( dir + "/a.to_left", to_left );
     BinaryWriter::writeFile( dir + "/a.to_right", to_right );
     BinaryWriter::writeFile( dir + "/a.inv", inv );
     HyperBasevectorX hbx(hb);
     BinaryWriter::writeFile( dir + "/a.hbx", hbx );
     {
          // Write edges and kmers.

          cout << Date( ) << ": making and writing edges" << endl;
          {    vecbvec edges( hb.Edges( ).begin( ), hb.Edges( ).end( ) );
               edges.WriteAll( dir + "/a.fastb" );    }
          {    vec<int> kmers( hb.E( ) );
               for ( int e = 0; e < hb.E( ); e++ )
                    kmers[e] = hb.Kmers(e);
               BinaryWriter::writeFile( dir + "/a.kmers", kmers );    }

          // Align to genome.

          if (ALIGN)
          {    const int align_genome_K = 80;
               GenomeAlign<align_genome_K>( hbx, inv, genome, alignsb );
               alignsb.WriteAll( dir + "/a.alignsb" );    }
          else Remove( dir + "/a.alignsb" );    }
     cout << Date() << ": reading paths" << endl;
     paths.ReadAll( dir + "/a.paths");
}
