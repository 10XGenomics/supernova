// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
//MakeDepend: library OMP
//MakeDepend: cflags OMP_FLAGS

#include "10X/PathsIndex.h"
#include <omp.h>
#include "10X/DfTools.h"

/* Approach to inversion:
 * Split step:
 * Split the set of edges in hbv into chunks.
 * Go through paths and create a temporary file for each chunk
 * that stores {(edge, read id) for each edge in chunk}
 * Process chunk:
 * Read in each temporary file. Sort the {(edge, read id)} vector.
 * Use a Parallel sort. Read id's that belong to the same edge
 * are now consecutive entries. 
 * write out the final paths_index file.
 */
// NEW: Also write a counts vector that measure the total read support


void writePathsIndex( ReadPathVec & paths, vec<int> & inv,
                      String base_dir, String pi_file = "a.paths.inv",
                      String counts_file = "a.countsb",
                      const int chunks = 15, bool verbose=false )
{
     cout << Date( ) << ": inverting paths" << endl;
     const int num_edges = inv.size();
     ForceAssertGt( num_edges, 0 );
     ForceAssertGt( chunks, 0 );
     const int chunk_size = num_edges/chunks + 1;
     ForceAssertGt(chunk_size, 0);

     if ( !IsDirectory( base_dir ) ) {
          FatalErr(base_dir + " is not a valid directory");
     }
     // Split into chunks 
     //
     if (verbose)
          cout << Date( ) << ": writing chunks to disk" << endl;
     vec<BinaryWriter *> out_chunks;
     vec<String> part_fn;
     for ( int i = 0; i != chunks; i++ ) {
          String fn = base_dir + "/tmp.part."+ToString(i);
          part_fn.push_back(fn);
          BinaryWriter * writer = new BinaryWriter(fn);
          out_chunks.push_back( writer );
     }
     
     for ( uint64_t id = 0; id < uint64_t(paths.size()); id++ ) {
          const ReadPath & p = paths[id];
          for ( auto e : p ) {
               const int chunk = e / chunk_size;
               out_chunks[chunk]->write( make_pair(e, id) );
          }
     }

     for ( auto * fptr : out_chunks ) {
          fptr->close();
          delete fptr;
     }
     if (verbose)
          cout << Date( ) << ": processing " << chunks << " chunks" << endl;
     
     // Process each chunk separately
     // Do not parallelize because then all of paths_index will
     // end up in memory

     // Which read ids are un-mapped to the graph?
     vec<bool> empty (num_edges, true);
     int numdots=0, done=0;
     // counts vector
     const int ns = 1;
     vec<vec<int>> countsb( ns, vec<int>(num_edges, -1) );
     IncrementalWriter<ULongVec> piw ( base_dir + "/" + pi_file );
     for ( int i = 0; i < chunks; i++ ) {
          vec<pair<int, unsigned long>> eid;
          uint64_t size = FileSize( part_fn[i] )/sizeof(pair<int, unsigned long>);
          eid.resize( size );
          BinaryReader reader (part_fn[i]);
          int64_t pos=0;
          while ( ! reader.atEOF() ) {
               pair<int, unsigned long> p;
               reader.read( &p );
               eid[pos] = p;
               pos++;
          }
          // Sort eid to put paths_index entries in consecutive order
          ParallelSort(eid);
          for ( auto & p : eid ) {
               if ( empty[ p.first ] )
                    empty[ p.first ] = false;
          }

          unsigned long j = 0;
          // Write read ids ONLY into paths_index part file
          const int emax = Min( (i+1)*chunk_size, num_edges );
          for ( int e = i*(chunk_size); e < emax; e++ ) {
               ULongVec entry;
               if ( empty[e] )
                    piw.add( entry );
               else {
                    while( j < eid.size() && eid[j].first == e ) {
                         entry.push_back( eid[j].second );
                         j++;
                    }
                    piw.add( entry );
               }
               countsb[0][e] = entry.size();
          }
          if (verbose)
               MakeDots( done, numdots, chunks );
     }
     piw.close();
     // check that all countsb elements are initialized
     for ( int n = 0; n < ns; n++ ) {
          for ( auto c : countsb[n] )
               ForceAssertGe ( c, 0 );
     }
     
     // add in the inverse edges
     #pragma omp parallel for collapse(2)
     for ( int n = 0; n < ns; n++ ) {
          for ( int e = 0; e < num_edges; e++ ) {
               const int re = inv[e];
               if ( e < re ) {
                    const int tot = countsb[n][e] + countsb[n][re];
                    countsb[n][e]  = tot;
                    countsb[n][re] = tot;
               }
          }
     }
     // write countsb
     BinaryWriter::writeFile( base_dir + "/" + counts_file, countsb );

     // Delete all the initially created tmp files
     if (verbose)
          cout << Date( ) <<": deleting tmp files" << endl;
     for ( auto fn : part_fn )
          SystemSucceed( "rm " + fn );

     cout << Date( ) << ": done" << endl;
     
}

void writePathsIndex( ReadPathVecX & paths, const HyperBasevectorX& hb, vec<int> & inv,
                      String base_dir, String pi_file = "a.paths.inv",
                      String counts_file = "a.countsb",
                      const int chunks = 15, bool verbose=false )
{
     cout << Date( ) << ": inverting pathsX" << endl;
     const int num_edges = inv.size();
     ForceAssertGt( num_edges, 0 );
     ForceAssertGt( chunks, 0 );
     const int chunk_size = num_edges/chunks + 1;
     ForceAssertGt(chunk_size, 0);

     if ( !IsDirectory( base_dir ) ) {
          FatalErr(base_dir + " is not a valid directory");
     }
     // Split into chunks 
     //
     if (verbose)
          cout << Date( ) << ": writing chunks to disk" << endl;
     vec<BinaryWriter *> out_chunks;
     vec<String> part_fn;
     for ( int i = 0; i != chunks; i++ ) {
          String fn = base_dir + "/tmp.part."+ToString(i);
          part_fn.push_back(fn);
          BinaryWriter * writer = new BinaryWriter(fn);
          out_chunks.push_back( writer );
     }
     const int pbatch = omp_get_max_threads()*1000;
     int pbatches = paths.size()/pbatch + 1;
     int numdots = 0, done=0;
     ReadPathVec bpaths;
     for ( uint64_t start = 0; start < (uint64_t) (paths.size()); start+=pbatch ) {
          const uint64_t stop = Min( uint64_t(paths.size()), start+pbatch );
          bpaths.resize(stop-start);
          if (start == 0) {
               for ( auto & p : bpaths )
                    p.reserve(2);
          }
          #pragma omp parallel for
          for ( uint64_t id = start; id < stop; id++ ) {
               ReadPath & p = bpaths[id-start];
               paths.unzip(p,hb,id);
          }
          for ( uint64_t id = start; id < stop; id++ ) {
               for ( auto & e: bpaths[id-start] ) {
                    const int chunk = e / chunk_size;
                    { out_chunks[chunk]->write( make_pair(e, id) ); }
               }
          }
          if (verbose ) {
               MakeDots( done, numdots, pbatches );
          }
     }
     Destroy(bpaths);

     for ( auto * fptr : out_chunks ) {
          fptr->close();
          delete fptr;
     }
     if (verbose)
          cout << Date( ) << ": processing " << chunks << " chunks" << endl;
     
     // Process each chunk separately
     // Do not parallelize because then all of paths_index will
     // end up in memory

     // Which read ids are un-mapped to the graph?
     vec<bool> empty (num_edges, true);
     // counts vector
     const int ns = 1;
     vec<vec<int>> countsb( ns, vec<int>(num_edges, -1) );
     numdots=0;
     done=0;
     IncrementalWriter<ULongVec> piw ( base_dir + "/" + pi_file );
     for ( int i = 0; i < chunks; i++ ) {
          vec<pair<int, unsigned long>> eid;
          uint64_t size = FileSize( part_fn[i] )/sizeof(pair<int, unsigned long>);
          eid.resize( size );
          BinaryReader reader (part_fn[i]);
          int64_t pos=0;
          while ( ! reader.atEOF() ) {
               pair<int, unsigned long> p;
               reader.read( &p );
               eid[pos] = p;
               pos++;
          }
          // Sort eid to put paths_index entries in consecutive order
          ParallelSort(eid);
          for ( auto & p : eid ) {
               if ( empty[ p.first ] )
                    empty[ p.first ] = false;
          }

          unsigned long j = 0;
          // Write read ids ONLY into paths_index part file
          const int emax = Min( (i+1)*chunk_size, num_edges );
          for ( int e = i*(chunk_size); e < emax; e++ ) {
               ULongVec entry;
               if ( empty[e] )
                    piw.add( entry );
               else {
                    while( j < eid.size() && eid[j].first == e ) {
                         entry.push_back( eid[j].second );
                         j++;
                    }
                    piw.add( entry );
               }
               countsb[0][e] = entry.size();
          }
          if (verbose)
               MakeDots( done, numdots, chunks );
     }
     piw.close();
     
     // check that all countsb elements are initialized
     for ( int n = 0; n < ns; n++ ) {
          for ( auto c : countsb[n] )
               ForceAssertGe ( c, 0 );
     }
     
     // add in the inverse edges
     #pragma omp parallel for collapse(2)
     for ( int n = 0; n < ns; n++ ) {
          for ( int e = 0; e < num_edges; e++ ) {
               const int re = inv[e];
               if ( e < re ) {
                    const int tot = countsb[n][e] + countsb[n][re];
                    countsb[n][e]  = tot;
                    countsb[n][re] = tot;
               }
          }
     }
     // write countsb
     BinaryWriter::writeFile( base_dir + "/" + counts_file, countsb );

     // Delete all the initially created tmp files
     if (verbose)
          cout << Date( ) <<": deleting tmp files" << endl;
     for ( auto fn : part_fn )
          SystemSucceed( "rm " + fn );
     
     cout << Date( ) << ": done" << endl;
     
}


/* Compute the base edge --> barcode vector map
 * Stored as VecIntVec ebcx
 */

void computeEdgeToBarcodeX( const ReadPathVecX & paths, const HyperBasevectorX& hb, 
        const vec <int32_t> & bc, const vec<int> & inv, 
        const vec<int64_t> & bci, VecIntVec & ebcx, bool verbose=false )
{
     cout << Date( ) << ": constructing ebcx" << endl;
     const int num_batches = bci.size()-1;
     const int num_edges   = inv.size();
     // for each batch of reads (with the same barcode)
     // the set of base assembly edges that they map to
     vec<vec<int>> bce(num_batches); 
     #pragma omp parallel for
     for ( int bi = 0; bi < num_batches; bi++ ) {
          const int64_t start = bci[bi];
          const int64_t stop  = bci[bi+1];
          if ( bc[start] <= 0 )
               continue;
          vec<int> edges;
          for ( int64_t id = start; id < stop; id++ ) {
              ReadPath p; paths.unzip(p,hb,id);
               for (const auto & e : p ) {
                    edges.push_back( e );
                    edges.push_back( inv[e] );
               }
          }
          UniqueSort(edges);
          bce[bi] = edges;
     }
     if (verbose)
          cout << Date( ) << ": computed bce" << endl;
     // Count edges
     vec <int> num_bc (num_edges, 0);
     for ( int bi = 0; bi < num_batches; bi++ ) {
          for ( auto e : bce[bi] ) {
               num_bc[e]++;
          }
     }
     
     // allocate memory
     ebcx.resize( num_edges );
     for ( int e = 0; e < num_edges; e++ )
          ebcx[e].resize(num_bc[e]);

     vec <int> pos ( num_edges, 0 );
     
     if (verbose)
          cout << Date( ) << ": filling up ebcx" << endl;
     // Construct the index
     for ( int bi = 0; bi < num_batches; bi++ ) {
          const int32_t b = bc[bci[bi]];
          if (b > 0) {
               vec<int> & edges = bce[bi];
               #pragma omp parallel for
               for ( int ei = 0; ei < edges.isize(); ei++ ) {
                    const int e = edges[ei];
                    ebcx[e][pos[e]] = b ;
                    pos[e]++;
               }
          }
     }
     if (verbose)
          cout << Date( ) << ": done" << endl;
}



