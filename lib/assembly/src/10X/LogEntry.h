// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

#include <iostream>
#include <string>
#include <vector>
#include <typeinfo>
#include <memory>
#include "system/System.h"
#include "Vec.h"
#include <type_traits>
#include <typeinfo>
#include <typeindex>

// Base-Derived class structure to hold *any* data type
// For use in StatLogger
// Base class stores (type, cs, name, gloss)
// Derived class stores data. Can be of any C++ type.
//
// data is accessible from Base class using a virtual function to_string
// that converts data to std::string object. This is useful to read/write
// from a binary file. This currently only works for
// INT, DOUBLE, STRING types.
//
struct LogEntryBase {
     
     virtual ~LogEntryBase() { }
     virtual void print() {}
     enum class TYPE:int8_t { INT=0, DOUBLE=1, STRING=2, INVALID=3 };
     TYPE   type;        // type of object
     virtual void to_string( std::string & data_string ) {}
     virtual void to_double( double & d ) {}

     bool   cs;          // is this stat CS?
     String name;        // name of the metric
     String gloss;       // description of metric
};

TRIVIALLY_SERIALIZABLE(LogEntryBase::TYPE);
TRIVIALLY_SERIALIZABLE(long long int);

template <class T>
struct LogEntry : public LogEntryBase {
     
     LogEntry() {}
     LogEntry( String name1, T val1, String gloss1, bool cs1 ) {
          type = TYPE::INVALID;   // default type is invalid
          cs   = cs1;
          name = name1;
          data = val1;
          gloss = gloss1;
     }
     ~LogEntry() {}

     void print () override {
          if (cs)
               cout << "CS ";
          else
               cout << "PD ";
          cout << "metric: " << name; 
          if ( gloss.size() > 0 )
               cout << " (" << gloss << ")"
                         << endl;
          PRINT(data);
     }

     void to_string( std::string & data_string ) override {
          ostringstream ostr;
          ostr << data;
          data_string = ostr.str();
     }
     
     void to_double( double & d ) override {
          if ( type != TYPE::INT && type != TYPE::DOUBLE )
               FatalErr ("Cannot convert " + name + " to double");
          stringstream ss;
          ss << data;
          ss >> d;
     }
     
     T data;
};

