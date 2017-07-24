// Copyright (c) 2016 10X Genomics, Inc. All rights reserved.

#ifndef READPATHDEFS_H_
#define READPATHDEFS_H_

#include "Vec.h"
#include "paths/long/ReadPath.h"
#include "paths/HyperBasevector.h"
#include "feudal/MasterVec.h"
#include "feudal/BinaryStream.h"

enum CODING_SCHEME{FIXED_WIDTH, HUFFMAN};

typedef vec<unsigned char> uCharArray;

typedef uCharArray::iterator uchararr_iterator;
typedef uCharArray::const_iterator const_uchararr_iterator;

typedef VirtualMasterVec<ReadPath> VReadPathVec;

#endif /* READPATHDEFS_H_ */
