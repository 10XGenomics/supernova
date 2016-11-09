/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2014) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include <ctype.h>

#include "graphics/BasicGraphics.h"
#include "CoreTools.h"
#include "FastIfstream.h"

Bool graphics_primitive::Valid( ) const
{    int count = 0;
     for ( int i = 0; i < (int) ps_.size( ) - 1; i++ )
     {    if ( ps_[i] == '$' && ps_[i+1] >= '1' && ps_[i+1] <= '9' )
          {    int j, arg = ps_[i+1] - '0';
               for ( j = i + 2; j < (int) ps_.size( ); j++ )
               {    if ( ps_[j] == ' ' )
                    {    ++count;
                         if ( arg != count ) return False;    }
                    else if ( isdigit( ps_[j] ) )
                    {    arg = ( 10 * arg ) + ( ps_[j] - '0' );    }
                    else break;    }    }    }
     return count == NCoords( ) * 2;    }

String graphics_primitive::Render( ) const
{    String answer = ps_;
     int count = 0, shift = 0;
     for ( int i = 0; i < (int) ps_.size( ) - 1; i++ )
     {    if ( ps_[i] == '$' && ps_[i+1] >= '1' && ps_[i+1] <= '9' )
          {    int j, arg = ps_[i+1] - '0';
               for ( j = i + 2; j < (int) ps_.size( ); j++ )
               {    if ( ps_[j] == ' ' )
                    {    ForceAssertLt( count, NCoords( ) * 2 );
                         double v;
                         if ( count % 2 == 0 ) v = Coord(count/2).x;
                         else v = Coord((count-1)/2).y;
                         ++count;
                         ForceAssertEq( arg, count );
                         answer.replace( i + shift, j - i, Ps(v) );
                         shift += Ps(v).size( ) + i - j;    }
                    else if ( isdigit( ps_[j] ) )
                    {    arg = ( 10 * arg ) + ( ps_[j] - '0' );    }
                    else break;    }    }    }
     ForceAssertEq( count, NCoords( ) * 2 );
     return answer;    }

double MinX( const vec<graphics_primitive>& g )
{    static vec<double> xs;
     xs.clear( );
     for ( int i = 0; i < g.isize( ); i++ )
     {    for ( int j = 0; j < g[i].NCoords( ); j++ )
          {    xs.push_back( g[i].Coord(j).x );    }    }
     ForceAssertGt( xs.size( ), 0u );
     return Min(xs);    }

double MaxX( const vec<graphics_primitive>& g )
{    static vec<double> xs;
     xs.clear( );
     for ( int i = 0; i < g.isize( ); i++ )
     {    for ( int j = 0; j < g[i].NCoords( ); j++ )
          {    xs.push_back( g[i].Coord(j).x );    }    }
     ForceAssertGt( xs.size( ), 0u );
     return Max(xs);    }

double MinY( const vec<graphics_primitive>& g )
{    static vec<double> ys;
     ys.clear( );
     for ( int i = 0; i < g.isize( ); i++ )
     {    for ( int j = 0; j < g[i].NCoords( ); j++ )
          {    ys.push_back( g[i].Coord(j).y );    }    }
     ForceAssertGt( ys.size( ), 0u );
     return Min(ys);    }

double MaxY( const vec<graphics_primitive>& g )
{    static vec<double> ys;
     ys.clear( );
     for ( int i = 0; i < g.isize( ); i++ )
     {    for ( int j = 0; j < g[i].NCoords( ); j++ )
          {    ys.push_back( g[i].Coord(j).y );    }    }
     ForceAssertGt( ys.size( ), 0u );
     return Max(ys);    }

int TotalCoords( const vec<graphics_primitive>& g )
{    int count = 0;
     for ( int i = 0; i < g.isize( ); i++ )
          count += g[i].NCoords( );
     return count;    }

graphics_primitive Segment( double x1, double y1, double x2, double y2 )
{    vec<coordinate> coords(2);
     coords[0] = coordinate(x1, y1), coords[1] = coordinate(x2, y2);
     String ps = "newpath $1 $2 moveto $3 $4 lineto stroke";
     return graphics_primitive( ps, coords );    }

graphics_primitive Segment( double x1, double y1, double x2, double y2, color c )
{    vec<coordinate> coords(2);
     coords[0] = coordinate(x1, y1), coords[1] = coordinate(x2, y2);
     String ps = "currentrgbcolor "
          + Ps(c.R()) + " " + Ps(c.G()) + " " + Ps(c.B()) + " setrgbcolor "
          + "newpath $1 $2 moveto $3 $4 lineto stroke setrgbcolor";
     return graphics_primitive( ps, coords );    }

graphics_primitive DottedSegment( double x1, double y1, double x2, double y2 )
{    vec<coordinate> coords(2);
     coords[0] = coordinate(x1, y1), coords[1] = coordinate(x2, y2);
     String ps = "newpath $1 $2 moveto $3 $4 lineto "
          "[2 2] 6 setdash stroke [] 0 setdash";
     return graphics_primitive( ps, coords );    }

graphics_primitive DottedSegment( double x1, double y1, double x2, double y2,
     color c )
{    vec<coordinate> coords(2);
     coords[0] = coordinate(x1, y1), coords[1] = coordinate(x2, y2);
     String ps = "currentrgbcolor "
          + Ps(c.R()) + " " + Ps(c.G()) + " " + Ps(c.B()) + " setrgbcolor "
          + "newpath $1 $2 moveto $3 $4 lineto "
          + "[2 2] 6 setdash stroke [] 0 setdash setrgbcolor";
     return graphics_primitive( ps, coords );    }

graphics_primitive TextToRight( const String& s, double x, double y, double kern )
{    ForceAssertGt( s.size( ), 0u );
     vec<coordinate> coords(1);
     coords[0] = coordinate(x, y);
     String ps = PsStringHeight( s.substr(0,1) ) + " newpath " + Ps(kern) 
          + " $1 add exch neg $2 add moveto " + "(" + s + ") show";
     return graphics_primitive( ps, coords );    }

graphics_primitive 
     TextToRight( const String& s, double x, double y, double kern, const font& f )
{    ForceAssertGt( s.size( ), 0u );
     vec<coordinate> coords(1);
     coords[0] = coordinate(x, y);
     String ps = "currentfont " + f.ps( ) + " setfont "
          + PsStringHeight( s.substr(0,1) ) + " newpath " + Ps(kern) 
          + " $1 add exch neg $2 add moveto " + "(" + s + ") show setfont";
     return graphics_primitive( ps, coords );    }

String PsStringHeight( const String& s )
{    return "gsave newpath 0 0 moveto (" + s + ") false charpath flattenpath "
          + "pathbbox grestore 3 -1 roll add 2 div 3 -1 roll pop exch pop";    }

vec<double> DefineAxis( double a, double b )
{    
     if ( !( a < b ) ) FatalErr( "DefineAxis failed." );

     // Let top = b-a, rounded down to the nearest integer power of 10.

     double top = pow( 10, floor(log10(b-a)) );

     // Write b-a as tics * atop + lower, where:
     // - tics is an integer, 2 <= tics <= 19;
     // - atop is an integer power of 10;
     // - lower is a real number, 0 <= lower < atop.

     double atop = top;
     int first_digit = int((b-a)/top);
     if ( first_digit == 1 )
     {    first_digit = 10;
          atop /= 10;    }
     int tics = first_digit;
     while( tics * atop <= b-a ) ++tics;

     // The first tic is the smallest integer multiple of atop that is >= a.

     double first_tic = ceil(a/atop) * atop;

     // Generate the tics.

     static vec<double> answer;
     answer.clear( );
     for ( int i = 0; ; i++ )
     {    double tic = first_tic + i * atop;
          if ( tic > b ) break;
          answer.push_back(tic);    }
     return answer;    }

void MakeTicLabels( const vec<double>& tics, vec<String>& labels )
{    labels.clear_and_resize( tics.size( ) );
     for ( int i = 0; i < tics.isize( ); i++ )
     {    labels[i] = ToString( tics[i], 10 );
          while( labels[i].Contains( "." ) && labels[i].Contains( "0", -1 ) )
               labels[i].resize( labels[i].size( ) - 1 );    }
     Bool alldotend = True, alldot = True;
     for ( int i = 0; i < tics.isize( ); i++ )
     {    if ( !labels[i].Contains( ".", -1 ) ) alldotend = False;
          if ( !labels[i].Contains( "." ) ) alldot = False;    }
     if (alldotend)
     {    for ( int i = 0; i < tics.isize( ); i++ )
               labels[i].resize( labels[i].size( ) - 1 );    }
     else if (alldot)
     {    int digits_to_right = 0;
          for ( int i = 0; i < tics.isize( ); i++ )
          {    digits_to_right = Max( digits_to_right, 
                    labels[i].isize( ) - labels[i].Position( "." ) - 1 );    }
          for ( int i = 0; i < tics.isize( ); i++ )
          {    int right = labels[i].isize( ) - labels[i].Position( "." ) - 1;
               for ( int j = right; j < digits_to_right; j++ )
                    labels[i] += "0";    }    }    }

vec<graphics_primitive> AxisX( double a, double b, double scale,
     Bool use_minor_tics, const String& label_suffix, const double extend,
     const Bool naked )
{    vec<graphics_primitive> answer;
     double low = a - extend * (b-a), high = b + extend * (b-a);
     answer.push_back( Segment( low, 0, high, 0 ) );
     vec<double> tics = DefineAxis( low, high );
     vec<String> labels;
     vec<double> scaled_tics(tics);
     scaled_tics *= scale;
     MakeTicLabels( scaled_tics, labels );
     for ( int j = 0; j < labels.isize( ); j++ )
          labels[j] += label_suffix;
     double fontsize = 8.0;
     double vkern = 1.5;
     if ( !naked ) answer.push_back( SetTimesRoman(fontsize) );
     for ( int i = 0; i < tics.isize( ); i++ )
     {    vec<coordinate> tic_loc(1);
          tic_loc[0].x = tics[i], tic_loc[0].y = 0;
          if ( !naked ) answer.push_back( Seg( tics[i], 0, 10, -90 ) );
          String labelid = labels[i];
          graphics_primitive 
               label( "newpath $1 $2 moveto " 
                    + PsStringWidth( labelid.substr(0,1) ) + " 2 div neg "
                    + " -10 " + Ps(fontsize) + " sub " + Ps(vkern) 
                    + " sub rmoveto (" + labelid + ") show", tic_loc );
          if ( !naked ) answer.push_back(label);    }
     tics.push_back( tics.back( ) + tics[1] - tics[0] );
     tics.insert( tics.begin( ), tics.front( ) - tics[1] + tics[0] );
     for ( int i = 0; i < tics.isize( ) - 1; i++ )
     {    for ( double j = 1.0; j <= 9.0; j += 1.0 )
          {    double start = tics[i] + j/10.0 * (tics[i+1]-tics[i]);
               if ( start < low || start > high ) continue;
               double ticheight = ( j == 5.0 ? 20.0/3.0 : 10.0/3.0 );
               if ( !naked ) 
                    answer.push_back( Seg( start, 0, ticheight, -90 ) );    }    }
     return answer;    }

graphics_primitive SetTimesRoman( double height )
{    vec<coordinate> coords;
     String ps = "/Times-Roman findfont " + Ps(height) + " scalefont setfont";
     return graphics_primitive( ps, coords );    }

graphics_primitive Seg( double x, double y, double l, double theta )
{    vec<coordinate> coords(1);
     coords[0] = coordinate(x, y);
     String ps = "newpath $1 $2 moveto " + Ps(theta) + " rotate " + Ps(l)
          + " 0 rlineto " + Ps(-theta) + " rotate stroke";
     return graphics_primitive( ps, coords );    }

graphics_primitive Seg( double x, double y, double l, double theta, color c )
{    vec<coordinate> coords(1);
     coords[0] = coordinate(x, y);
     String ps = "currentrgbcolor "
          + Ps(c.R()) + " " + Ps(c.G()) + " " + Ps(c.B()) + " setrgbcolor "
          + "newpath $1 $2 moveto " + Ps(theta) + " rotate " + Ps(l) 
          + " 0 rlineto " + Ps(-theta) + " rotate stroke setrgbcolor";
     return graphics_primitive( ps, coords );    }

graphics_primitive RectangleBaseCentered( double x, double y, double h,
     double w, color fill, Bool outline )
{    vec<coordinate> coords(4);
     coords[0] = coordinate(x - w/2, y);
     coords[1] = coordinate(x + w/2, y);
     coords[2] = coordinate(x + w/2, y + h);
     coords[3] = coordinate(x - w/2, y + h);
     String ps = "currentrgbcolor "
          "0 0 0 setrgbcolor "
          "newpath /Point1 { $1 $2 } def Point1 moveto $3 $4 lineto $5 $6 lineto "
          "$7 $8 lineto Point1 lineto closepath gsave "
          + String( outline ? "stroke " : "" ) + "grestore "
          + Ps(fill.R()) + " " + Ps(fill.G()) + " " + Ps(fill.B()) + " setrgbcolor "
          + "fill setrgbcolor";
     return graphics_primitive( ps, coords );    }

vec<graphics_primitive> DotPath( const vec<coordinate> p, double r )
{    vec<graphics_primitive> answer;
     for ( int i = 0; i < p.isize( ); i++ )
     {    answer.push_back( Point( p[i].x, p[i].y, r ) );
          if ( i < p.isize( ) - 1 )
               answer.push_back(
                    Segment( p[i].x, p[i].y, p[i+1].x, p[i+1].y ) );    }
     return answer;    }

graphics_primitive TextToLeft( const String& s, double x, double y, double kern )
{    ForceAssertGt( s.size( ), 0u );
     vec<coordinate> coords(1);
     coords[0] = coordinate(x, y);
     String ps = "newpath $1 $2 moveto " + PsStringWidth(s) + " neg " + Ps(kern)
          + " sub " + PsStringHeight( s.substr( s.size( ) - 1, 1 ) )
          // + " 2 div neg rmoveto (" + s + ") show";
          + " neg rmoveto (" + s + ") show";
     return graphics_primitive( ps, coords );    }

graphics_primitive 
     TextCenter( const String& s, double x, double y, double kern, const font& f )
{    ForceAssertGt( s.size( ), 0u );
     vec<coordinate> coords(1);
     coords[0] = coordinate(x, y);
     String ps = "currentfont " + f.ps( ) + " setfont "
          + "newpath $1 $2 moveto " + PsStringWidth(s) + " 2 div neg "
          + Ps(kern) + " sub " + PsStringHeight( s.substr( s.size( ) - 1, 1 ) )
          + " neg rmoveto (" + s + ") show setfont";
     return graphics_primitive( ps, coords );    }

graphics_primitive RainbowTextCenter( const vec<String>& s, const vec<color>& c, 
     double x, double y, double kern, const font& f )
{    ForceAssertGt( s.size( ), 0u );
     ForceAssertEq( s.size( ), c.size( ) );
     String sall;
     for ( int i = 0; i < s.isize( ); i++ )
     {    ForceAssertGt( s[i].size( ), 0u );
          sall += s[i];    }
     vec<coordinate> coords(1);
     coords[0] = coordinate(x, y);
     String ps = "currentfont " + f.ps( ) + " setfont "
          + "newpath $1 $2 moveto " + PsStringWidth(sall) + " 2 div neg "
          + Ps(kern) + " sub " + PsStringHeight( sall.substr( sall.size( ) - 1, 1 ) )
          + " neg rmoveto"; 
     for ( int i = 0; i < s.isize( ); i++ )
     {    ps += " " + Ps(c[i].R()) + " " + Ps(c[i].G()) + " " + Ps(c[i].B()) 
               + " setrgbcolor (" + s[i] + ") show";    }
     ps += " setfont";
     return graphics_primitive( ps, coords );    }

vec<graphics_primitive> AxisY( double x1, double x2, double y1, double y2,
     Bool use_minor_tics, const String& label_suffix, const double extend,
     const Bool naked )
{    vec<graphics_primitive> answer;
     double low = y1 - extend * (y2-y1), high = y2 + extend * (y2-y1);
     answer.push_back( Segment( x1, low, x1, high ) );
     vec<double> tics = DefineAxis( y1, y2 );
     vec<String> labels;
     MakeTicLabels( tics, labels );
     for ( int j = 0; j < labels.isize( ); j++ )
          labels[j] += label_suffix;
     double fontsize = 8.0;
     double hkern = 1.5;
     if ( !naked ) answer.push_back( SetTimesRoman(fontsize) );
     for ( int i = 0; i < tics.isize( ); i++ )
     {    if ( !naked ) answer.push_back( Seg( x1, tics[i], 10, 180 ) );
          if ( !naked ) answer.push_back( 
               TextToLeft( labels[i], x1, tics[i], 10 + hkern ) );    }
     tics.push_back( tics.back( ) + tics[1] - tics[0] );
     tics.insert( tics.begin( ), tics.front( ) - tics[1] + tics[0] );
     for ( int i = 0; i < tics.isize( ) - 1; i++ )
     {    for ( double j = 1.0; j <= 9.0; j += 1.0 )
          {    if ( j != 5.0 && !use_minor_tics ) continue;
               double start = tics[i] + j/10.0 * (tics[i+1]-tics[i]);
               if ( start < low || start > high ) continue;
               double ticheight = ( j == 5.0 ? 20.0/3.0 : 10.0/3.0 );
               if ( !naked ) 
                    answer.push_back( Seg( x1, start, ticheight, 180 ) );    }    }
     return answer;    }

String RenderGraphics( const String& OUT, 
     const vec< vec<graphics_primitive> >& stack,
     const vec<double>& heights, const double SCALE, const double XTRANS, 
     const double YTRANS, const Bool SHOW, const double POSTSCALE,
     const Bool allow_fail )
{    Bool to_png = OUT.Contains( ".png", -1 );
     Bool to_ps = OUT.Contains( ".ps", -1 );
     Bool to_pdf = OUT.Contains( ".pdf", -1 );
     ForceAssert( to_png || to_ps || to_pdf );
     if ( to_png && !SHOW )
     {    cout << "For a png file, you must specify SHOW=True.\n";
          exit(-1);    }
     String outhead;
     if (to_ps) outhead = OUT.Before( ".ps" );
     else if (to_pdf) outhead = OUT.Before( ".pdf" );
     else outhead = OUT.Before( ".png" );
     double width = 400;
     {    Ofstream( out, outhead + ".ps" );
          VerticalDisplay( stack, heights, width, out, SCALE, XTRANS, 
               YTRANS, SHOW );    }
     if (to_pdf)
     {    System( "ps2pdf " + outhead + ".ps" );
          Remove( outhead + ".ps" );    }
     if (to_png)
     {    String command1 = "ps2epsi " + outhead + ".ps " + outhead + ".eps";
          String command2 = "set -o pipefail; "
               "( pstopnm -portrait -xmax 8000 -ymax 8000 -stdout "
               + outhead + ".eps | pnmcrop | "
               + "pnmpad -white -left=25 -right=25 -top=25 -bottom=25 | ";
          if ( POSTSCALE != 1.0 ) 
               command2 += "pamscale " + ToString(POSTSCALE,3) + " | ";
          command2 += "pnmtopng > " + outhead + ".png )";
          for ( int pass = 1; pass <= 2; pass++ )
          {    String command = ( pass == 1 ? command1 : command2 );
               if ( !allow_fail )
               {    int status = System(command);
                    if ( status == 0 ) Remove( outhead + ".ps" );
                    else 
                    {    cout << "failed to run:\n" << command << "\n";     
                         cout << "Abort.\n";
                         TracebackThisProcess( );    }    }
               else
               {    int status = System( command + " 2> " + outhead + ".log" );
                    if ( status == 0 ) 
                    {    Remove( outhead + ".ps" ), Remove( outhead + ".log" );    }
                    else
                    {    String fail_msg = "failed to run:\n" + command + "\n";
                         if ( IsRegularFile( outhead + ".log" ) )
                         {    fast_ifstream in( outhead + ".log" );
                              String line;
                              while(1)
                              {    getline( in, line );
                                   if ( in.fail( ) ) break;
                                   fail_msg += line + "\n";    }
                              Remove( outhead + ".log" );    }
                         else fail_msg += "Oddly, no logfile was generated.\n";
                         Remove( outhead + ".ps" ); 
                         return fail_msg;    }    }    }    }
     return "";    }

void CheckValidGraphicsSuffix( const String& OUT, const Bool SHOW )
{    Bool to_png = OUT.Contains( ".png", -1 );
     Bool to_ps = OUT.Contains( ".ps", -1 );
     Bool to_pdf = OUT.Contains( ".pdf", -1 );
     ForceAssert( to_png || to_ps || to_pdf );
     ForceAssert( to_png || to_ps || to_pdf );
     if ( to_png && !SHOW )
     {    cout << "For a png file, you must specify SHOW=True.\n";
          exit(-1);    }    }
