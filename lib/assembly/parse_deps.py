#!/usr/bin/env python

import os
import argparse
import sys
from collections import defaultdict
import glob

parser = argparse.ArgumentParser()
parser.add_argument( "deps", help="Makefile_deps file" )
parser.add_argument( "--clean_src", help ="lib/assembly folder of some repo to clean (WARNING: this deletes files!!)" )
args = parser.parse_args()

deps_file = args.deps

deps = {}
lcount = 0
for line in open( deps_file, "r" ):
     data = line.strip()
     if len(data) == 0:
          continue
     if data.count(": ") != 1:
          continue
     if data.count(":=") > 0:
          continue
     target, need = data.split(": ")
     if target in deps:
          continue
     deps[target] = need.split()
     lcount+=1

executables = ["DF", "CP", "MakeFasta", "ParseBarcodedFastqs", "FastFastbCount"]

def traverse( keys, dep_dict ):
     top_keys    = []
     lookup_keys = []
     
     for k in keys:
          if k in dep_dict:
               lookup_keys.extend( dep_dict[k] )
          else:
               top_keys.append( k )
     if len( lookup_keys ) > 0:
          return (top_keys + traverse( lookup_keys, dep_dict ))
     else:
          return top_keys

keep = set(traverse( executables, deps ))

## clean up repo

if args.clean_src:
     repo = args.clean_src
     if not os.path.exists( repo ):
          print "invalid path"
          sys.exit(1)
     os.chdir( args.clean_src )
     ## keep all top-level files
     for fn in glob.glob("*"):
          if os.path.isfile( fn ):
               keep.add( fn )
     keep.add ("src/system/MakeDepend.cc")
     keep.add ("src/LinkTimestamp.cc")
     
     print "in dir", os.path.realpath(".")
     for base, _, files in os.walk("."):
          for f in files:
               key = os.path.join(base, f)[2:]
               keybase = ".".join(key.split(".")[:-1])
               print key, 
               if (key in keep) or ((keybase + ".h") in keep ) or ( (keybase + ".cc") in keep):
                    print
               else:
                    print "DELETE"
                    os.remove( "./"+key )
     ## remove any empty dirs
     while True:
          empty_dirs = []
          for base, dirnames, files in os.walk("."):
               if len(files) == 0 and len(dirnames) == 0:
                    empty_dirs.append( base )
          if len(empty_dirs) == 0:
               break
          for d in empty_dirs:
               print "dir", d, "is empty, deleting"
               os.rmdir( d )

