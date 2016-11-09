///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * Fosmids.cc
 *
 *  Created on: Jul 15, 2014
 *      Author: tsharpe
 */

#include "paths/long/fosmid/Fosmids.h"
#include "math/HoInterval.h"

vec<int> ChrSizes(void)
{
     return {249250621,243199373,198022430,191154276,180915260,171115067,
          159138663,146364022,141213431,135534747,135006516,133851895,115169878,
          107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,
          51304566,155270560};
}

void GetRegionInfo( const String& ID, int& g, int& rstart, int& rstop,
     String& gid, String& loc )
{    Ifstream( hin, "/wga/dev/references/Homo_sapiens/NA12878_Fosmid_Pool.regions" );
     String line, fid;
     while(1)
     {    getline( hin, line );
          if ( hin.fail( ) )
          {    cout << "failed to find Fosmid " << ID << endl;
               _exit(1);    }
          istringstream iline( line.c_str( ) );
          iline >> fid >> loc;
          if ( fid != ID ) continue;
          gid = loc.Before( ":" );
          if ( gid == "X" ) g = 22;
          else g = gid.Int( ) - 1;
          rstart = loc.Between( ":", "-" ).Int( );
          rstop = loc.After( "-" ).Int( );
          break;    }    }

vec<String> AllFosmidRegions( size_t padding )
{    vec<String> regions;
     vec<int> const sizes = ChrSizes();
     vec< pair<String,ho_interval> > regs;
     for ( auto const fid : AllFosmids() )
     {    int g, rstart, rstop;
          String gid, loc;
          GetRegionInfo(ToString(fid), g, rstart, rstop, gid, loc );
          rstart -= padding; rstart = max( rstart, 0 );
          rstop += padding; rstop = min( rstop, sizes[g] );
          regs.push( gid, ho_interval(rstart,rstop) );   }
     Sort(regs);
     for ( int i = 0; i < regs.isize( ); i++ )
     {    int j, M = 0;
          vec<ho_interval> in, out;
          for ( j = i; j < regs.isize( ); j++ )
          {    if ( regs[j].first != regs[i].first ) break;
               in.push_back( regs[j].second );
               M = Max( M, regs[j].second.Stop( ) );    }
          ExtractGivenCoverage( M, 1, in, out );
          for ( int l = 0; l < out.isize( ); l++ )
          {    ostringstream s;
               s << regs[i].first << ":" << out[l];
               regions.push_back( s.str( ) );    }
          i = j - 1;    }
    return regions;    }

String FosmidRegion( int fid, size_t padding )
{
     vec<int> const sizes = ChrSizes();
     int g, rstart, rstop;
     String gid, loc;
     GetRegionInfo(ToString(fid), g, rstart, rstop, gid, loc );
     rstart -= padding; rstart = max( rstart, 0 );
     rstop += padding; rstop = min( rstop, sizes[g] );
     ostringstream s;
     s << gid << ":" << rstart << "-" << rstop;
     return s.str();
}
