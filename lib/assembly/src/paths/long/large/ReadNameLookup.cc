///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"
#include "ParallelVecUtilities.h"
#include "TokenizeString.h"
#include "paths/long/large/ReadNameLookup.h"

uint64_t readname_lookup::KeyFromName( String name )
{    if ( name.Contains( ".1", -1 ) ) name = name.RevBefore( ".1" );
     else 
     {    if ( !name.Contains( ".2", -1 ) )
          {    cout << "Illegal readname " << name << "." << endl;
               cout << "Readnames are required to end with .1 or .2." << endl;
               Scram(1);    }
          name = name.RevBefore( ".2" );    }
     uint64_t id = 0, M = 1;
     vec<String> n;
     Tokenize( name, ':', n );
     if ( n.size( ) != top_.size( ) )
     {    cout << "Illegal readname " << name << endl;
          cout << "Number of fields doesn't match." << endl;
          Scram(1);    }
     for ( int l = 0; l < top_.isize( ); l++ )
     {    uint64_t x;
          if ( l == fcpos_ ) 
          {    int p = Position( fcnames_, n[l] );
               if ( p < 0 )
               {    cout << "Illegal readname " << name << "." << endl;
                    cout << "Flowcell name doesn't match." << endl;
                    Scram(1);    }
               x = p;    }
          else 
          {    if ( !n[l].IsInt( ) )
               {    cout << "Illegal readname " << name << "." << endl;
                    cout << "Non-integer field in unexpected position." << endl;
                    Scram(1);    }
               x = n[l].Int( );
               if ( x > top_[l] )
               {    cout << "Illegal readname " << name << "." << endl;
                    cout << "Field value exceeds top." << endl;
                    Scram(1);    }    }
          id += M * x;
          M *= top_[l] + 1;    }
     return id;    }

uint64_t readname_lookup::GetReadId( const String& n )
{    uint64_t key = KeyFromName(n);
     int64_t x = BinPosition( keys_, key );
     ForceAssert( x >= 0 );
     return ( 2 * (int64_t) pids_[x] ) + ( n.Contains( ".1", -1 ) ? 0 : 1 );    }

readname_lookup::readname_lookup( const vecString& names )
{
     cout << Date( ) << ": entering readname_lookup constructor" << endl;
     ForceAssert( names.size( ) > 0 );
     ForceAssertLe( (uint64_t) names.size( ), 2 * UINT32_MAX );

     // Check pairing structure.

     #pragma omp parallel for
     for ( int64_t i = 0; i < (int64_t) names.size( ); i++ )
     {    if ( i % 2 == 0 ) 
          {    if ( !names[i].Contains( ".1", -1 ) )
               {    cout << "Names not following .1, .2 convention." << endl;
                    Scram(1);    }    }
          else 
          {    if ( !names[i].Contains( ".2", -1 ) )
               {    cout << "Names not following .1, .2 convention." << endl;
                    Scram(1);    }
               ForceAssertEq( names[i].RevBefore( ".2" ), 
                    names[i-1].RevBefore( ".1" ) );    }    }

     // From first name, determine layout.

     cout << Date( ) << ": determining layout" << endl;
     vec<int> nonints;
     vec<String> n1;
     Tokenize( names[0].RevBefore( ".1" ) , ':', n1 );
     for ( int j = 0; j < n1.isize( ); j++ )
     {    String x = n1[j];
          ForceAssert( x.size( ) > 0 );
          if ( x[0] == 0 ) ForceAssert( x.size( ) == 1 );
          Bool digits = True;
          for ( int l = 0; l < x.isize( ); l++ )
               if ( !isdigit( x[l] ) ) digits = False;
          if ( !digits ) nonints.push_back(j);    }
     ForceAssert( nonints.solo( ) );
     fcpos_ = nonints[0];

     // Check entire file and determine tops.

     cout << Date( ) << ": checking all" << endl;
     int len = n1.size( );
     top_.resize( len, 0 );
     set<String> fcnames_s;
     const int64_t nbatches = 100;
     int64_t npids = names.size( ) / 2;
     vec<vec<uint64_t>> topb( nbatches, vec<uint64_t>(len,0) );
     vec< set<String> > fcnames_sb( nbatches );
     #pragma omp parallel for
     for ( int64_t bi = 0; bi < nbatches; bi++ )
     {    
          #pragma omp critical
          {    cout << "starting batch " << bi+1 << endl;    }
          vec<String> n;
          vec<int> nonints;
          for ( int64_t pi = (npids*bi)/nbatches; 
                    pi < (npids*(bi+1))/nbatches; pi++ )
          {    int64_t i = pi*2;
               nonints.clear( );
               Tokenize( names[i].RevBefore( ".1" ), ':', n );
               ForceAssertEq( n.isize( ), len );
               fcnames_sb[bi].insert( n[fcpos_] );
               for ( int j = 0; j < n.isize( ); j++ )
               {    const String& x = n[j];
                    ForceAssert( x.size( ) > 0 );
                    if ( x[0] == 0 ) ForceAssert( x.size( ) == 1 );
                    Bool digits = True;
                    for ( int l = 0; l < x.isize( ); l++ )
                         if ( !isdigit( x[l] ) ) digits = False;
                    if ( !digits ) nonints.push_back(j);    }
               if ( !nonints.solo( ) )
               {    cout << "Wrong number of noninteger fields: " 
                         << nonints.size( ) << "." << endl;
                    cout << "From: " << names[i] << endl;
                    Scram(1);    }
               if ( nonints[0] != fcpos_ )
               {    cout << "Noninteger field in wrong position." << endl;
                    cout << "From: " << names[i] << endl;
                    Scram(1);    }
               for ( int j = 0; j < len; j++ )
               {    if ( j != fcpos_ )
                    {    ForceAssert( n[j].IsInt( ) );
                         topb[bi][j] = Max( topb[bi][j], 
                              (uint64_t) n[j].Int( ) );    }    }    }    }
     for ( int64_t bi = 0; bi < nbatches; bi++ )
     {    for ( set<String>::iterator i = fcnames_sb[bi].begin( );
               i != fcnames_sb[bi].end( ); i++ )
          {    fcnames_s.insert(*i);    }
          for ( int j = 0; j < len; j++ )
               if ( j != fcpos_ ) top_[j] = Max( top_[j], topb[bi][j] );    }
     for ( set<String>::iterator i = fcnames_s.begin( ); i != fcnames_s.end( ); i++ )
          fcnames_.push_back(*i);
     top_[fcpos_] = fcnames_.size( ) - 1;

     // Check for eight-byte fit.

     uint64_t prod = 1;
     for ( int l = 0; l < len; l++ )
     {    uint64_t tp = top_[l] + 1;
          ForceAssertLt( prod, UINT64_MAX/tp );
          prod *= tp;    }

     // Translate readnames.

     cout << Date( ) << ": translating readnames" << endl;
     keys_.resize( names.size( ) / 2 );
     #pragma omp parallel for
     for ( int64_t i = 0; i < (int64_t) names.size( ); i += 2 )
          keys_[i/2] = KeyFromName( names[i] );

     // Sort.

     cout << Date( ) << ": sorting" << endl; 
     pids_ = vec<uint32_t>( names.size( )/2, vec<uint32_t>::IDENTITY );
     ParallelSortSync( keys_, pids_ );
     cout << Date( ) << ": checking" << endl;
     for ( int64_t i = 1; i < (int64_t) keys_.size( ); i++ )
     {    if ( keys_[i] == keys_[i-1] )
          {    cout << "Found duplicate key = " << keys_[i] << ",\n"
                    << "from pids " << pids_[i] << " and " << pids_[i-1] << ",\n"
                    << "associated to pairnames\n"
                    << names[ 2 * pids_[i] ] << " and\n"
                    << names[ 2 * pids_[i-1] ] << "." << endl;
               Scram(1);    }    }
     cout << Date( ) << ": done\n";    }

void readname_lookup::writeBinary( BinaryWriter& writer ) const
{    writer.write(fcpos_);
     writer.write(top_);
     writer.write(fcnames_);
     writer.write(keys_);
     writer.write(pids_);    }

void readname_lookup::readBinary( BinaryReader& reader )
{    reader.read(&fcpos_);
     reader.read(&top_);
     reader.read(&fcnames_);
     reader.read(&keys_);
     reader.read(&pids_);    }
