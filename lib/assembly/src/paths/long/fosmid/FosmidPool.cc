///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "FastIfstream.h"
#include "TokenizeString.h"
#include "VecUtilities.h"
#include "paths/long/fosmid/FosmidPool.h"

void ParseFosmidPoolMetainfo( vec<String>& regions, 
     vec< vec< pair<String,String> > >& junctions,
     vec< vec< pair<String,String> > >& breaks,
     vec< vec< pair<String,String> > >& edits )
{    fast_ifstream in( 
          "/wga/dev/references/Homo_sapiens/NA12878_Fosmid_Pool.regions" );
     String line;
     while(1)
     {    getline( in, line );
          if ( in.fail( ) ) break;

          // Delete C++-style /*...*/ comments.

          String line2;
          for ( int j = 0; j < line.isize( ); j++ )
          {    if ( j < line.isize( ) - 1 && line[j] == '/' && line[j+1] == '*' )
               {    j += 2;
                    Bool found = False;
                    for ( ; j < line.isize( ) - 1; j++ )
                    {    if ( line[j] == '*' && line[j+1] == '/' )
                         {    found = True;
                              break;    }    }
                    ForceAssert(found);
                    j += 2;    }
               if ( j == line.isize( ) ) break;
               line2.push_back( line[j] );    }
          line = line2;

          vec<String> fields, changes;
          vec< pair<String,String> > j1, b1, e1;
          Tokenize( line, {' '}, fields );
          ForceAssert( fields.isize( ) >= 2 && fields.isize( ) <= 3 );
          regions.push_back( fields[1] );
          if ( fields.isize( ) == 3 )
          {    Tokenize( fields[2], {','}, changes );
               for ( int j = 0; j < changes.isize( ); j++ )
               {    if ( changes[j] == "DELETED" ) continue;
                    else if ( changes[j] == "LOW_COVERAGE" ) continue;
                    else if ( changes[j].Contains( "|" ) ) 
                    {    j1.push( changes[j].Before( "|" ), 
                              changes[j].After( "|" ) );    }
                    else if ( changes[j].Contains( "$" ) ) 
                    {    b1.push( changes[j].Before( "$" ), 
                              changes[j].After( "$" ) );    }
                    else if ( changes[j].Contains( "-->" ) )
                    {    e1.push( changes[j].Before( "-->" ),
                              changes[j].After( "-->" ) );    }    }    }
          junctions.push_back(j1), breaks.push_back(b1); 
          edits.push_back(e1);    }    }

vec<String> GetFinishedFosmidFiles( )
{    String dir = "/wga/dev/references/Homo_sapiens/NA12878_Fosmid_Pool.regions.fin";
     vec<String> all = AllFiles(dir);
     vec<Bool> to_delete( all.size( ), False );
     for ( int i = 0; i < all.isize( ); i++ )
          if ( !all[i].Contains( ".fasta", -1 ) ) to_delete[i] = True;
     EraseIf( all, to_delete );
     vec<int> ids;
     for ( int i = 0; i < all.isize( ); i++ )
          ids.push_back( all[i].Between( ".", "." ).Int( ) );
     SortSync( ids, all );
     return all;    }
