/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2014) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// This is an experimental toolkit for generating postscript graphics, such as
// those encountered when displaying data.

// Rough explanation of usage: first you generate a vec of graphics_primitive
// objects (or vec of vec), using functions such as Segment, TextToRight, etc.
// Then you call Display or VerticalDisplay.

// Caution.  Some commands which yield text (such as TextToRight) will generate
// invisible characters if you have not first declared a font (using e.g.
// SetTimesRoman).

// EXAMPLE: See InsertPlotter.cc.  This is the only example at present!

#ifndef BASIC_GRAPHICS_H
#define BASIC_GRAPHICS_H

#include "graphics/Color.h"
#include "CoreTools.h"
#include "math/Functions.h"

// Class coordinate: keep track of cartesian coordinates.

class coordinate {
     public:
     double x, y;
     coordinate( ) { }
     coordinate( double x_arg,
		 double y_arg )
       : x(x_arg), y(y_arg)
     { } 
};

// Ps defines the "standard" translation from floating point numbers to postscript
// strings.

inline String Ps( double x ) { return ToString( x, 6 ); }

// Class font: represent a font.

class font {
     public:
     font( const String& ps ) : ps_(ps) { }
     String ps( ) const { return ps_; }
     private:
     String ps_;
};

// TimesRoman(h): represent Times Roman font at height h.

inline font TimesRoman( double height )
{    return font( "/Times-Roman findfont " + Ps(height) + " scalefont" );    }

inline font TimesBold( double height )
{    return font( "/Times-Bold findfont " + Ps(height) + " scalefont" );    }

inline font ArialBold( double height )
{    return font( "/Arial-Bold findfont " + Ps(height) + " scalefont" );    }

// A graphics_primitive object represents a postscript string, with space for
// coordinates, which are stored outside the string.  It has two components:
// the string, and the coordinates vector.  For example, a line from (0,1) to
// (5,7) could be represented as the string
//             "newpath $1 $2 moveto $3 $4 lineto stroke"
// and the coordinates vector with two entries
//             (0,1)
//             (5,7).
// As implemented, if the coordinates vector has n entries, then $1, $2, ..., $2n
// must appear in the string, in that order.
//
// Although there are no restrictions placed on the postscript which may appear in 
// the string, there are many possibilities that are inappropriate, given the way
// in which graphics_primitive objects are used.  For example, you shouldn't use
// translate, or otherwise change the coordinate system.
//
// Dimensions in the string are absolute (measured in points), whereas the
// dimensions in the coordinates vector may be externally scaled.

class graphics_primitive {

     public:

     graphics_primitive( ) { }
     
     graphics_primitive( const String& ps, const vec<coordinate>& coords )
          : ps_(ps), coords_(coords) { }

     int NCoords( ) const { return coords_.size( ); }

     coordinate Coord( int i ) const { return coords_[i]; }

     // Valid: check to see if the tokens $1, etc. are right.

     Bool Valid( ) const;

     // AffineMap( a1, b1, a2, b2 ): modify coords_ via
     // x |--> a1x+b1, y |--> a2x + b2.

     void AffineMap( double a1, double b1, double a2, double b2 ) 
     {    for ( int i = 0; i < NCoords( ); i++ )
          {    coords_[i].x = a1 * Coord(i).x + b1;
               coords_[i].y = a2 * Coord(i).y + b2;    }    }

     // Render: Substitute coords_ into ps_, generating (hopefully) valid 
     // postscript.

     String Render( ) const;

     private:

     String ps_;
     vec<coordinate> coords_;

};

// Point( x, y, r ): create a point having center (x,y) and radius r.

inline graphics_primitive Point( double x, double y, double r )
{    vec<coordinate> coords(1);
     coords[0] = coordinate(x, y);
     String ps = "$1 $2 " + Ps(r) + " 0.2 setflat newpath 0 360 arc fill 1 setflat";
     return graphics_primitive( ps, coords );    }

// SetColor( color c ): set color to c.

inline graphics_primitive SetColor( const color& c )
{    vec<coordinate> coords;
     String ps = Ps(c.R()) + " " + Ps(c.G()) + " " + Ps(c.B()) + " setrgbcolor";
     return graphics_primitive( ps, coords );    }

inline graphics_primitive SetLineWidth( float l )
{    vec<coordinate> coords;
     String ps = ToString(l, 6) + " setlinewidth";
     return graphics_primitive( ps, coords );    }

// Segment( x1, y1, x2, y2 ): generate a line segment from (x1,y1) to (x2,y2).
// Optional fifth argument: color of segment.

graphics_primitive Segment( double x1, double y1, double x2, double y2 );
graphics_primitive Segment( double x1, double y1, double x2, double y2, color c );

// DottedSegment( x1, y1, x2, y2 ): generate a dotted line segment from 
// (x1,y1) to (x2,y2).  Optional fifth argument: color.

graphics_primitive DottedSegment( double x1, double y1, double x2, double y2 );
graphics_primitive DottedSegment( double x1, double y1, double x2, double y2, 
     color c );

// Seg( x, y, l, theta ): generate a line segment from (x,y) to the point
// l absolute units away, in the direction theta.

graphics_primitive Seg( double x, double y, double l, double theta );
graphics_primitive Seg( double x, double y, double l, double theta, color c );

// RectangleBaseCentered.  Generate a rectangle of height h and width w, whose 
// base is centered at (x,y), filled with a given color, and outlined in black,
// unless outline = False is specified.

graphics_primitive RectangleBaseCentered( double x, double y, double h, 
     double w, color fill, Bool outline = True );

// DotPath: generate a sequence of line segments between the given points,
// using points of radius r.

vec<graphics_primitive> DotPath( const vec<coordinate> p, double r );

// SetTimesRoman: set Times-Roman font at the given size.

graphics_primitive SetTimesRoman( double height );

// PsStringWidth: return postscript which computes the width of a given string.

inline String PsStringWidth( const String& s )
{    return "(" + s + ") stringwidth pop";    }

// PsStringHeight: return postscript which computes the "height" of a given string.
// What is actually computed is the average of the top and bottom positions of the
// string relative to the baseline (see "add 2 div").

String PsStringHeight( const String& s );

// TextToRight( s, x, y, kern ): print string s at (x,y) + (kern,-height(s[0])/2).
// Optional last argument: font.

graphics_primitive TextToRight( const String& s, double x, double y, double kern );
graphics_primitive 
     TextToRight( const String& s, double x, double y, double kern, const font& f );

// TextToLeft( s, x, y, kern ): print string s at (x,y) + (-kern,-height(s.back)/2).

graphics_primitive TextToLeft( const String& s, double x, double y, double kern );

graphics_primitive 
     TextCenter( const String& s, double x, double y, double kern, const font& f );

graphics_primitive RainbowTextCenter( const vec<String>& s, const vec<color>& c,
     double x, double y, double kern, const font& f );

// DefineAxis: given real numbers ranging from a to b, define appropriate places
// for tic marks on an axis.

vec<double> DefineAxis( double a, double b );

// AxisX: define an x-axis appropriate for values ranging from a to b.
// This yields major tics, semimajor tics (equally spaced between major tics),
// and minor tics (at tenth spacing).
//
// Warning: use_minor_tics not implemented.

vec<graphics_primitive> AxisX( double a, double b, double scale = 1.0,
     Bool use_minor_tics = True, const String& label_suffix = "",
     const double extend = 0.05, const Bool naked = False );

// AxisY: define a y-axis appropriate for x values ranging from x1 to x2 and
// y values ranging from y1 to y2.  This yields major tics, semimajor tics (equally 
// spaced between major tics), and minor tics (at tenth spacing).  The latter may
// be turned off.

vec<graphics_primitive> AxisY( double x1, double x2, double y1, double y2,
     Bool use_minor_tics = True, const String& label_suffix = "",
     const double extend = 0.05, const Bool naked = False );

// MinX, MaxX, MinY, MaxY: return lower and upper bounds on coordinates in a
// vec of graphics_primitive objects.

double MinX( const vec<graphics_primitive>& g );
double MaxX( const vec<graphics_primitive>& g );
double MinY( const vec<graphics_primitive>& g );
double MaxY( const vec<graphics_primitive>& g );

// TotalCoords: return total number of coordinates appearing in a vec of
// graphics_primitive objects.

int TotalCoords( const vec<graphics_primitive>& g );

inline void AffineMap( vec<graphics_primitive>& g, double a1, double b1, 
     double a2, double b2 )
{    for ( int i = 0; i < g.isize( ); i++ )
          g[i].AffineMap( a1, b1, a2, b2 );    }

inline void Render( const vec<graphics_primitive>& g, ostream& out )
{    for ( int i = 0; i < g.isize( ); i++ )
          out << g[i].Render( ) << "\n";    }

// Display: Let g be a vec<graphics_primitive> object.  Scale g so it has dimension 
// w x h.  If all coordinates in g have the same y value, then g is centered.

inline void Display( vec<graphics_primitive> g, double w, double h, ostream& out )
{    if ( TotalCoords(g) == 0 ) FatalErr( "Display: no coordinates passed." );
     double x = MinX(g), X = MaxX(g), y = MinY(g), Y = MaxY(g);
     if ( x == X ) FatalErr( "Display failed." );
     
     // Scale to land in [0,w] x [0,h].

     if ( Y > y )
     {    double mx = w/(X-x), my = h/(Y-y);
          AffineMap( g, mx, -mx * x, my, -my * y );    }
     else
     {    double mx = w/(X-x);
          AffineMap( g, mx, -mx * x, 1, h/2 - y );    }

     out << "200 50 translate\n";
     Render( g, out );
     out << "showpage\n";    }

// VerticalDisplay: Let g[0], g[1], ... be vec<graphics_primitive> objects, which
// share a common x-coordinate system.  Arrange these on a page with g[0] on the
// bottom, g[1] above it, and so forth.  Scale g[i] so it has dimension w x h[i].  
// If all coordinates in g[i] have the same y value, then g[i] is centered.

inline void VerticalDisplay( vec< vec<graphics_primitive> > g, const vec<double>& h, 
     double w, ostream& out, double scale = 1.0, double xtrans = 200.0,
     double ytrans = 50.0, Bool show = True )
{    int n = g.size( );
     ForceAssertGt( n, 0 );
     for ( int i = 0; i < n; i++ )
     {    if ( TotalCoords( g[i] ) == 0 )
          {    FatalErr( "VerticalDisplay: no coordinates present in g[" 
                    << i << "]." );    }    }
     vec<double> xi(n), Xi(n), yi(n), Yi(n);
     for ( int i = 0; i < n; i++ )
     {    xi[i] = MinX( g[i] );
          Xi[i] = MaxX( g[i] );     
          yi[i] = MinY( g[i] );
          Yi[i] = MaxY( g[i] );    }
     double x = Min(xi), X = Max(Xi);
     if ( x == X ) FatalErr( "VerticalDisplay failed." );
     double floor = 0.0;
     for ( int i = 0; i < n; i++ )
     {    if ( i > 0 ) floor += h[i-1];
          if ( Yi[i] > yi[i] )
          {    double mx = w/(X-x), myi = h[i]/(Yi[i]-yi[i]);
               AffineMap( g[i], mx, -mx * x, myi, floor + (-myi * yi[i]) );    }
          else
          {    double mx = w/(X-x);
               AffineMap( g[i], mx, -mx * x, 1, floor + h[i]/2 - yi[i] );    }    }
     out << xtrans << " " << ytrans << " translate\n";
     out << scale << " " << scale << " scale\n";
     for ( int i = 0; i < n; i++ )
          Render( g[i], out );
     if (show) out << "showpage\n";    }

// RenderGraphics: if allow_fail = True, certain failure modes will be caught,
// and in that case a string explaining the failure is returned.  Otherwise the
// empty string is returned.

String RenderGraphics( const String& OUT, 
     const vec< vec<graphics_primitive> >& stack,
     const vec<double>& heights, const double SCALE = 1.0, 
     const double XTRANS = 200, const double YTRANS = 50, const Bool SHOW = True,
     const double POSTSCALE = 1.0, const Bool allow_fail = False );

void CheckValidGraphicsSuffix( const String& OUT, const Bool SHOW = True );

#endif
