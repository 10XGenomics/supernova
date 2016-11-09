// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

#ifndef TENX_DF_TOOLS_H
#define TENX_DF_TOOLS_H

#if defined(SHOWMEM)
#define MEM(X) { cout << #X << ": mem = " << MemUsageGBString( ) << ", peak = " << PeakMemUsageGBString( ) << endl; }
#else
#define MEM(X)
#endif

#include "CoreTools.h"
#include "Intvector.h"
#include "feudal/ObjectManager.h"
#include "feudal/PQVec.h"
#include "feudal/BinaryStream.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "10X/LogEntry.h"
#include "10X/Martian.h"
#include "10X/paths/ReadPathVecX.h"

enum class ReadDataType : uint8_t {
     PCR_FREE,            /* PCR-free conventional data */
     PCR,                 /* PCR-amplified conventional data */
     UNBAR_10X,           /* 10X nominally-barcoded data WITHOUT a passing barcode */
     BAR_10X              /* 10X barcoded data WITH a passing barcode */
};

struct DataSet {
     ReadDataType dt;
     int64_t start;

     friend ostream& operator<<( ostream& out, DataSet const& val ) {
          if ( val.dt == ReadDataType::PCR_FREE ) out << "PCR_FREE";
          else if ( val.dt == ReadDataType::PCR ) out << "PCR";
          else if ( val.dt == ReadDataType::UNBAR_10X ) out << "UNBAR_10X";
          else if ( val.dt == ReadDataType::BAR_10X ) out << "BAR_10X";
          else FatalErr("bad value for ReadDataType");
          out << " starts at " << val.start;
          return out;
     }
};

TRIVIALLY_SERIALIZABLE(DataSet);

class StatLogger
{
public:
     static void init(String const logfile = "", String const alertfile = "") {
          gInst.initThis();
          if (logfile != "" ) {
               BinaryReader::readFile(logfile, &gInst);
          }
          if (alertfile != "" ) {
               ifstream alerts_handle(alertfile);
               std::string str;
               while (std::getline(alerts_handle, str)) {
                    if (str == "#") {
                         vec <String> text;
                         double threshold=-1;
                         for ( int i = 0; i != 5; i++ ) {
                              std::string field;
                              std::getline(alerts_handle, field);
                              if ( i == 1 )
                                   threshold= String(field).Double();
                              else
                                   text.push_back( field );
                         }
                         ForceAssertGe(threshold, 0);
                         //(metric, compare, threshold, text, action)
                         gInst.mAlerts.push_back( _alert( text[0], text[1], threshold, text[3], text[2]) );
                    } else
                         FatalErr ("Alerts file " + alertfile + " incorrectly formatted");
               }
          }
     }

     static void setSilent( bool silent ) { gInst.mSilent = silent; }

     template <class T>
     static void log( String statName, T val, String gloss, bool cs = false ) {
          LogEntryBase * entry = new LogEntry<T>(statName, val, gloss, cs);
          bool is_int = is_integral<T>::value;
          bool is_double = is_floating_point<T>::value;
          bool is_string = std::is_assignable<String, T>::value;
          if (is_int)
               entry->type = LogEntryBase::TYPE::INT;
          else if (is_double)
               entry->type = LogEntryBase::TYPE::DOUBLE;
          else if (is_string)
               entry->type = LogEntryBase::TYPE::STRING;
          else
               entry->type = LogEntryBase::TYPE::INVALID;
          gInst.mpLines[statName] = entry;
          if ( !gInst.mSilent ) {
               cout << statName << "=" << val << "\t" << gloss << endl;
          }
     }
     
     // retrieve a stored stat as a double
     static double getNumStat( String statName, double default_value = 0.0 ) {
          double d = default_value;
          auto it = gInst.mpLines.find( statName );
          if ( it != gInst.mpLines.end() )
               it->second->to_double( d );
          return d;
     }

     // issue martian alerts
     template <class T>
     static void issue_alert( String metric_name, T value, String format_string = "" ) {
          vec <String> level, msg;
          for ( auto a : gInst.mAlerts ) {
               if (a._metric == metric_name) {
                    if ( a._is_greater ^ value < a._threshold ) {
                         level.push_back(a._action);
                         if ( format_string.size() > 0 )
                              msg.push_back(a.message(format_string));
                         else
                              msg.push_back(a.message(value));
                    }
               }
          }
          vec <String> uniq_level = level;
          UniqueSort(uniq_level);
          if ( uniq_level.size() != level.size() )
               FatalErr ("Multiple alerts of the same type triggered for " + metric_name + " = " + ToString(value) );
          // did we issue an exit, in which case don't display any other alert.
          String exit_msg = "";
          for ( int i = 0; i != level.isize(); i++ ) {
               if ( level[i] == "exit" )
                    exit_msg = msg[i];
          }
          if (exit_msg.size() > 0 )
               Martian::exit (exit_msg);
          else {
               for ( int i = 0; i != level.isize(); i++ )
                    Martian::log( level[i], msg[i] );
          }

     }


     // this function dumps out all statistics into a txt file
     // independent of the value of the flag cs
     static void dump_text( String const& filename ) {
          if ( gInst.mpLines.size() ) {
               ofstream out( filename );
               if ( out ) {
                    for ( auto entry : gInst.mpLines ) {
                         if ( entry.second->type == LogEntryBase::TYPE::INVALID )
                              continue;
                         std::string data_string;
                         entry.second->to_string(data_string);
                         out << entry.second->name << "\t" << data_string << "\t" << entry.second->gloss << endl;
                    }
               } else FatalErr("failed writing output to "+filename );
          }
     }
     
     // this function is used to define a user-facing csv
     // cs_only = true means only CS facing stats are written
     static void dump_csv( String const& filename, String const sep=",", bool cs_only = true ) {
          if ( gInst.mpLines.size() ) {
               ofstream out( filename );
               if ( out ) {
                    const int num_stats = gInst.mpLines.size();
                    bool first_entry=true;
                    for ( auto entry: gInst.mpLines ) {
                         if ( entry.second->type == LogEntryBase::TYPE::INVALID )
                              continue;
                         if ( cs_only && !entry.second->cs )
                              continue;
                         if ( !first_entry ) {
                              out << sep;
                         }
                         first_entry=false;
                         out << entry.second->name;
                    }
                    out <<endl;
                    first_entry=true;
                    for ( auto entry : gInst.mpLines ) {
                         if ( entry.second->type == LogEntryBase::TYPE::INVALID )
                              continue;
                         if ( cs_only && !entry.second->cs )
                              continue;
                         if ( !first_entry ) {
                              out << sep;
                         }
                         first_entry=false;
                         std::string data_string;
                         entry.second->to_string(data_string);
                         out << data_string;
                    }
                    out <<endl;
               } else FatalErr("failed writing output to "+filename );
          }
     }
     
     static void dump_json( String const& filename, Bool valsOnly = False ) {
          if ( gInst.mpLines.size() ) {
               ofstream out( filename );
               if ( out ) {
                    out << "{" << endl;
                    for ( auto itr = gInst.mpLines.begin(); itr != gInst.mpLines.end();  ) {
                         LogEntryBase * entry = itr->second;
                         if (entry->type == LogEntryBase::TYPE::INVALID) {
                              ++itr;
                              continue;
                         }
                         // add quotes if it is a string
                         std::string data_string;
                         entry->to_string( data_string );
                         if ( entry->type == LogEntryBase::TYPE::STRING )
                              data_string = "\"" + data_string + "\"";
                         if(valsOnly){
                             out  << '\t' << '"' << entry->name << '"' << ": " 
                                  << data_string; 
                         }else{
                             out  << '\t' << '"' << entry->name << '"' << ": { \"value\": " 
                                  << data_string << ", \"descr\": \"" << entry->gloss 
                                  << "\" }";
                         }
                         if ( ++itr != gInst.mpLines.end() ) out << "," << endl;
                    }
                    out << endl << "}" << endl;
               } else FatalErr("failed writing output to " + filename );
          }
     }

     static size_t externalSizeof();
     
     void writeBinary( BinaryWriter & writer ) const {
          writer.write(mSilent);
          int size = 0;
          for ( auto entry : mpLines ) {
               if ( entry.second->type != LogEntryBase::TYPE::INVALID )
                    size++;
          }
          writer.write(size);
          for ( auto itr : mpLines ) {
               if ( itr.second->type == LogEntryBase::TYPE::INVALID )
                    continue;
               LogEntryBase * entry = itr.second;
               writer.write(entry->type);
               writer.write(entry->cs);
               writer.write(entry->name);
               writer.write(entry->gloss);
               // read data into String
               std::string data_string;
               entry->to_string(data_string);
               switch ( entry->type ) {
                    case LogEntryBase::TYPE::INT: {
                         long long int i = atoll(data_string.c_str());
                         writer.write( i );
                         continue;
                    }
                    case LogEntryBase::TYPE::DOUBLE: {
                         double d = String(data_string).Double();
                         writer.write( d );
                         continue;
                    }
                    case LogEntryBase::TYPE::STRING: {
                         writer.write( String(data_string) );
                         continue;
                    }
                    case LogEntryBase::TYPE::INVALID: {
                         continue;
                    }
                    default: {
                         FatalErr("Unrecognized type");
                    }
               }
          }
          cout << Date( ) << ": successfully wrote " << size
               << " object(s) to log file" << endl;
     }
     
     void readBinary( BinaryReader & reader ) {
          reader.read( &mSilent );
          int size;
          reader.read(&size);
          cout << Date( ) << ": reading " << size << " logged object(s)" << endl;
          for ( auto a: mpLines)
               delete a.second;
          mpLines.clear();
          for ( int i = 0; i != size; i++ ) {
               String name, gloss;
               bool cs;
               LogEntryBase::TYPE type;
               reader.read( &type);
               reader.read( &cs);
               reader.read( &name);
               reader.read( &gloss);
               switch( type ) {
                    case LogEntryBase::TYPE::INT: {
                         long long int i;
                         reader.read( &i );
                         LogEntryBase *entry = new LogEntry<long long int>( name, i, gloss, cs);
                         entry->type = type;
                         mpLines[ name ] = entry;
                         continue;
                    }
                    case LogEntryBase::TYPE::DOUBLE: {
                         double d;
                         reader.read( &d );
                         LogEntryBase *entry = new LogEntry<double>( name, d, gloss, cs);
                         entry->type = type;
                         mpLines[ name ] = entry;
                         continue;
                    }
                    case LogEntryBase::TYPE::STRING: {
                         String s;
                         reader.read( &s );
                         LogEntryBase *entry = new LogEntry<String>( name, s, gloss, cs);
                         entry->type = type;
                         mpLines[ name ] = entry;
                         continue;
                    }
                    default: {
                         FatalErr("Unrecognized type");
                    }
               }
          }
     }

     static void write( const char* filename) {
          write(String(filename));
     }

     static void write( String const& filename ) {
          BinaryWriter::writeFile(filename, gInst);
     }


private:
     StatLogger() {};
     ~StatLogger() {
          for ( auto entry : mpLines )
               delete entry.second;
     }
     StatLogger(StatLogger const&) = delete;
     void operator=(StatLogger const&) = delete;

     void initThis() {
          gInst.mSilent = True;
          for ( auto a : gInst.mpLines )
               delete a.second;
          gInst.mpLines.clear();
     }
     struct _alert {
          _alert(String metric, String compare, double threshold, String text, String action) {
               _metric     = metric;
               _threshold  = threshold;
               _text       = text;
               _is_greater = (compare == ">");
               _action     = action;
          }
          _alert() {}
          ~_alert() {}
          template <class T>
          String message( T value ) {
               int pos = _text.find("{}");
               ostringstream valstream;
               valstream << std::fixed << std::setprecision(2);
               valstream << value;
               String s = _text;
               if ( pos >= 0 && pos < _text.isize() )
                    s = _text.replace(pos, 2, valstream.str());
               return s;
          }
          
          String _metric;
          bool   _is_greater;
          double _threshold;
          String _text;
          String _action;
     };


     bool mSilent;                 // is the logging silent or verbose
     std::map<String, LogEntryBase *> mpLines;  // holds the logging data
     vec <_alert> mAlerts;         // holds the complete set of alerts
     static StatLogger gInst;
};

SELF_SERIALIZABLE(StatLogger);

template<int K> void MapClosures( const HyperBasevectorX& hb,
     const vec<basevector>& closures, ReadPathVec& clop );

vec<int> GetBarcodes( const int e, const vec<int>& inv,
     const VecULongVec& paths_index, const vec<int>& bc );

void LoadData( const String& work_dir, const String& R, const vec<String>& lr,
     const vec<double>& LR_SELECT_FRAC, vecbasevector& bases,
     ObjectManager<VecPQVec>& quals_om, vec<int64_t>& bci,
     vec<String>& subsam_names, vec<int64_t>& subsam_starts, vec<DataSet>& datasets);

void GetQualStats( 
     const VecPQVec& quals, vec<vec<vec<int64_t>>> & hist, int & max_read_length);

void FragDist( const HyperBasevectorX& hb, const vec<int>& inv,
     const ReadPathVec& paths, vec<int64_t>& count );

void FragDist( const HyperBasevectorX& hb, const vec<int>& inv,
     const ReadPathVecX& paths, vec<int64_t>& count );

void ReadTwoPctProper( const HyperBasevectorX& hb, const vec<int>& inv,
     ReadPathVecX const& paths, double & r2_pct_proper );

typedef triple<int,int,int> trip;
extern template class SmallVec<trip,MempoolAllocator<trip> >;
extern template class OuterVec< SerfVec<trip> >;
#include "feudal/SmallVecDefs.h"
#include "feudal/OuterVecDefs.h"
template class SmallVec<trip,MempoolAllocator<trip> >;
template class OuterVec< SerfVec<trip> >;

int64_t EstimateGEMCount( const vec<int64_t> & bci, const int64_t TOTAL_DIVERSITY);
void SanityCheckBarcodeCounts( const vec<int64_t>& bci );

void MakeDots( int& done, int& ndots, const int total );

void SuperToSeqGraph( const HyperBasevectorX& hb, const digraphE<vec<int>>& D,
     HyperBasevectorX& hbd );

#endif
