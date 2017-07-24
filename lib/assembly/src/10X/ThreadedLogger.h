// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
//
// MakeDepend: cflags OMP_FLAGS
// 
#ifndef _THREADEDLOGGER_H
#define _THREADEDLOGGER_H

#include <stdexcept>
#include <sstream>
#include <vector>
#include <fstream>
#include <iostream>
#include <omp.h>
#include "system/SysConf.h"

class ThreadedLogger
{
public:
     ThreadedLogger( std::string filename, size_t nthreads = getConfiguredNumThreads() ) : mBuffers(nthreads), 
               mFile(filename, std::ios::out | std::ios::trunc ) {
          if ( !mFile.good() ) {
               throw std::runtime_error(std::string("error opening log file") + filename + " for writing");
          }
     }

     ThreadedLogger( char const* filename = "/dev/stdout", size_t nthreads = getConfiguredNumThreads() ) : mBuffers(nthreads), 
               mFile(filename, std::ios::out | std::ios::trunc ) {
          if ( !mFile.good() ) {
               throw std::runtime_error(std::string("error opening log file") + filename + " for writing");
          }
     }

     virtual ~ThreadedLogger() { join(); }

     std::ostringstream& get(size_t thread = omp_get_thread_num() ) { 
          if ( !mFile.is_open() ) {
               mBuffers[0].str("");
               return mBuffers[0];
          }
          if ( thread < mBuffers.size() ) {
               return mBuffers[thread];
          }
          else throw std::runtime_error("ThreadLogger received request for non-existent thread log.");
     }

     void join() {
          if ( !mFile.is_open() ) { 
               mBuffers[0].str("");
               return;
          }
          for ( auto& s : mBuffers )  {
               mFile << s.str();
               if ( mFile.bad() ) std::cerr << s.str();
               s.str("");
          }
     }

private:
     std::vector<std::ostringstream> mBuffers;
     std::ofstream mFile;
};

#endif // _THREADEDLOGGER_H
