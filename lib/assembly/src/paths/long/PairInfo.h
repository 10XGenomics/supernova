///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef PAIR_INFO
#define PAIR_INFO

#include "CoreTools.h"
#include "paths/long/Fix64_6.h"

class pair_point {

     public:

     pair_point( ) { }
     pair_point( const int trim, const fix64_6& weight, const int lib )
          : trim_(trim), weight_(weight), lib_(lib) { }

     int Trim( ) const { return trim_; }
     int& TrimMutable( ) { return trim_; }
     const fix64_6& Weight( ) const { return weight_; }
     fix64_6& WeightMutable( ) { return weight_; }
     int Lib( ) const { return lib_; }
     int& LibMutable( ) { return lib_; }

     friend Bool operator<( const pair_point& p1, const pair_point& p2 )
     {    if ( p1.trim_ < p2.trim_ ) return True;
          if ( p1.trim_ > p2.trim_ ) return False;
          if ( p1.weight_ < p2.weight_ ) return True;
          if ( p1.weight_ > p2.weight_ ) return False;
          return p1.lib_ < p2.lib_;    }

     friend Bool operator==( const pair_point& p1, const pair_point& p2 )
     {    return p1.trim_ == p2.trim_ && p1.weight_ == p2.weight_
               && p1.lib_ == p2.lib_;    }

     private:

     int trim_;
     fix64_6 weight_;
     int lib_;

};

template <> struct Serializability<pair_point> 
{ typedef TriviallySerializable type; };

// A pairing_info is intended to represent information associated to a single
// read.  It is intended to evolve.  The current version consists of:
// * status: either 0 (unpaired), 1 (1st read of pair), or 2 (2nd read of pair);
// * partner_id: id of partner read or -1 if unpaired;
// * lib_id: library id or -1 if unpaired.

class pairing_info {

     public:

     pairing_info( ) { }
     pairing_info( const int status, const int partner, const int lib_id )
          : status_(status), partner_(partner), lib_id_(lib_id) { }
     
     int Status( ) const { return status_; }
     int Partner( ) const { return partner_; }
     Bool Paired( ) const { return status_ >= 0; }
     int LibId( ) const { return lib_id_; }

     private:

     int status_;
     int partner_;
     int lib_id_;

};

// A gap_info is intended to represent information about the gap defined by
// a 'bundle' of paired reads.  This is intended to evolve.  It could be 
// anything between the simplest mean +/- dev description, and a full list of
// constituent measurements, including library references.

class gap_info {

     public:

     gap_info( ) { }

     gap_info( const double mean, const double dev ) : mean_(mean), dev_(dev) { }

     gap_info( const gap_info& gap1, const gap_info& gap2 )
     {    double g1 = gap1.Mean( ), g2 = gap2.Dev( );
          double d1 = gap1.Dev( ), d2 = gap2.Dev( );
          double w1 = 1.0/(d1*d1), w2 = 1.0/(d2*d2);
          double g = ( w1*g1 + w2*g2 ) / (w1+w2);
          double d = 1.0 / sqrt( w1+w2 );
          mean_ = g;
          dev_ = d;    }

     double Mean( ) const
     {    return mean_;    }
     double Dev( ) const
     {    return dev_;    }

     private:

     double mean_;
     double dev_;

};

template <> struct Serializability<gap_info> { typedef TriviallySerializable type; };

#endif
