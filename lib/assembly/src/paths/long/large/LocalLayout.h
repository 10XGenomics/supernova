///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef LOCAL_LAYOUT_H
#define LOCAL_LAYOUT_H
#include "feudal/CharString.h"
#include "paths/long/LongProtoTools.h"

//void LocalLayout( const int lroot, const int rroot, const String& TMP, const String& work_dir );
void LocalLayout( const int lroot, const int rroot, const LongProtoTmpDirManager&, const String& work_dir );

#endif
