///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// fix64_6: class to store a 64-bit fixed-precision number having 6 digits
// to the right of the decimal point.  Addition is associative.

#ifndef FIX64_6
#define FIX64_6

#include "CoreTools.h"

class fix64_6 {

     public:

     fix64_6( ) { }

     fix64_6( const int64_t n )
     {    val_ = n * 1000000;    }

     fix64_6( const int64_t n, const int64_t d )
     {    val_ = ( n * 1000000 ) / d;    }

     void fromRaw(const int64_t raw)
     { val_=raw; };

     friend Bool operator==( const fix64_6& x, const fix64_6& y )
     {    return x.val_ == y.val_;    }

     friend Bool operator!=( const fix64_6& x, const fix64_6& y )
     {    return x.val_ != y.val_;    }

     friend Bool operator>=( const fix64_6& x, const fix64_6& y )
     {    return x.val_ >= y.val_;    }

     friend Bool operator>( const fix64_6& x, const fix64_6& y )
     {    return x.val_ > y.val_;    }

     friend Bool operator<=( const fix64_6& x, const fix64_6& y )
     {    return x.val_ <= y.val_;    }

     friend Bool operator<( const fix64_6& x, const fix64_6& y )
     {    return x.val_ < y.val_;    }

     friend Bool operator>=( const fix64_6& x, const int64_t& y )
     {    return x.val_ >= y * 1000000;    }

     friend Bool operator>( const fix64_6& x, const int64_t& y )
     {    return x.val_ > y * 1000000;    }

     friend Bool operator<=( const fix64_6& x, const int64_t& y )
     {    return x.val_ <= y * 1000000;    }

     friend Bool operator<( const fix64_6& x, const int64_t& y )
     {    return x.val_ < y * 1000000;    }

     friend fix64_6 operator+( const fix64_6& x, const fix64_6& y )
     {    fix64_6 z;
          z.val_ = x.val_ + y.val_;
          return z;    }

     friend fix64_6 operator-( const fix64_6& x, const fix64_6& y )
     {    fix64_6 z;
          z.val_ = x.val_ - y.val_;
          return z;    }

     fix64_6& operator+=(const fix64_6& x)
     {    val_ += x.val_;
          return *this;    }

     friend fix64_6 operator*( const int64_t& x, const fix64_6& y )
     {    fix64_6 z;
          z.val_ = x * y.val_;
          return z;    }

     friend fix64_6 operator*( const fix64_6& y, const int64_t& x )
     {    fix64_6 z;
          z.val_ = x * y.val_;
          return z;    }

     friend fix64_6 operator/( const fix64_6& x, const int64_t& y )
     {    fix64_6 z;
          z.val_ = x.val_ / y;
          return z;    }

     fix64_6& operator/=(const int64_t& x)
     {    val_ /= x;
          return *this;    }

     friend ostream& operator<<( ostream& out, const fix64_6& x )
     {    int64_t z = x.val_;
          if ( z < 0 )
          {    z = -z;
               out << "-";    }
          int64_t left = z / 1000000, right = z % 1000000;
          out << left;
          if ( right > 0 ) 
          {    out << ".";
               String r = ToString(right);
               if ( right == 0 ) r = "0";
               else if ( right < 10 ) r = "00000" + r;
               else if ( right < 100 ) r = "0000" + r;
               else if ( right < 1000 ) r = "000" + r;
               else if ( right < 10000 ) r = "00" + r;
               else if ( right < 100000 ) r = "0" + r;
               while( r.size( ) > 1 && r.back( ) == '0' ) r.resize( r.size( ) - 1 );
               out << r;    }
          return out;    }

     double ToDouble( ) const
     {    return double(val_)/1000000.0;    }

     private:

     int64_t val_;

};

TRIVIALLY_SERIALIZABLE(fix64_6);

#endif
