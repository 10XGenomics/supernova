// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
//
#ifndef MARTIAN_H_INCLUDED
#define MARTIAN_H_INCLUDED

#include <iostream>
#include "system/System.h"
#include "Vec.h"

class Martian
{
public:
     static const int cMagic = 185;     // date of the first supernova observation

     static void init( std::string& dir ) {
          gInst.mInited = true;
          gInst.mFile = dir + "/" + gInst._stage_alerts;
          gInst.rFile = dir + "/" + gInst._rollup_alerts;
          if ( IsRegularFile( gInst.mFile ) )
               Remove( gInst.mFile );
     };

     static void init( const char* dir ) { 
          std::string sdir(dir);
          init(sdir);
     }

     static void init( String & dir ) {
          std::string sdir = std::string(dir);
          init(sdir);
     }

     static void log( std::string const& e, std::string const& msg ) {
          if ( !is_alert(e) ) FatalErr( "unknown alert tag: "+e );
          cout << "[" << e << "]" << " " << msg << endl;
          if ( gInst.mInited ) {
               gInst.mAlerts.push_back( std::make_pair( e, msg ) );
               // Update rollup alerts if not an exit
               if ( e != "exit" ) {
                    ofstream rollup;
                    rollup.open(gInst.rFile, ios::app);
                    rollup << msg << endl;
                    rollup.close();
               }
          }
     }

     static void exit( std::string const& msg ) {
          log( "exit", msg );
          Scram(cMagic);
     }

     static void alarm( std::string const& msg ) { 
          log( "alarm",  msg ); 
     }

     static void log_info( std::string const& msg ) { 
          log( "log_info", msg );
     }

     static void log_warn( std::string const& msg ) {
          log( "log_warn", msg );
     }


private:
     Martian() {};
     ~Martian() { 
          if ( gInst.mInited ) 
               gInst.dump_json(gInst.mFile);
     }

     void esc_json( std::ostream& o, std::string const& s ) {
          for ( auto c : s ) {
               if ( c == '"' ) o << "\\";
               o << c;
          }
     }

     void dump_json( std::string const& json_file ) {
          Ofstream( JSON, json_file );
          JSON << "{" << endl;
          for ( auto const& alert_tag : gAlertTypes ) {
               bool already = false;
               JSON << '"' << alert_tag << '"' << ": [";
               for ( auto const& alert : mAlerts ) {
                    if ( alert.first == alert_tag ) { 
                         if ( already ) JSON << ",";
                         else already = true;
                         JSON << '"';
                         esc_json( JSON, alert.second ); 
                         JSON << '"';
                    }
               }
               JSON << "]";
               if ( alert_tag != gAlertTypes.back() ) JSON << ",";
               JSON << endl;
          }
          JSON << "}" << endl;
     }

     static bool is_alert( std::string const& e ) {
          for ( auto const& alert : gAlertTypes ) {
               if ( alert == e ) return true;
          }
          return false;
     }
     
     std::string _stage_alerts  = "martian_alerts.json";
     std::string _rollup_alerts = "alerts_rollup.txt";
     std::string mFile;
     std::string rFile;
     bool mInited;
     std::vector<std::pair<std::string, std::string>> mAlerts;
     static const vec<std::string> gAlertTypes;
     static Martian gInst;
};

#endif // MARTIAN_H_INCLUDED
