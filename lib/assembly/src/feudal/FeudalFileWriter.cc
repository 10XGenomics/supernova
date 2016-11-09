///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * \file FeudalFileWriter.cc
 * \author tsharpe
 * \date Mar 13, 2009
 *
 * \brief Utility for writing Feudal Files incrementally.
 */
#include "feudal/FeudalFileWriter.h"
#include "system/Exit.h"
#include "system/System.h"
#include <iostream>
using std::cout;
using std::endl;
using std::string;

FeudalControlBlock FeudalFileWriter::gInvalidFCB(0,0,0,0,0,0);

FeudalFileWriter::FeudalFileWriter( char const* filename,
                                    FeudalFileWriter::size_type vecSize,
                                    FeudalFileWriter::size_type eltSize,
                                    FeudalFileWriter::size_type fixedLenDataLen,
                                    unsigned long estimatedNElements )
: mWriter(filename,false),
  mWriterOffsets(string(filename)+".off", false),
  mWriterFixed(string(filename)+".fix",false),
  mVecSize(vecSize),
  mEltSize(eltSize),
  mFixedLenDataLen(fixedLenDataLen),
  mElements(0)
{
    mWriter.write(gInvalidFCB);
    mWriterOffsets.write( mWriter.tell() );
}

FeudalFileWriter::~FeudalFileWriter()
{
    if ( mWriter.isOpen() )
        close();
}

void FeudalFileWriter::addElement( void const* fixedLenData )
{
    mWriterOffsets.write(mWriter.tell());
    mElements++;

    if ( mFixedLenDataLen )
    {
        char const* fStart = reinterpret_cast<char const*>(fixedLenData);
        mWriterFixed.write( fStart, fStart+mFixedLenDataLen );
    }
}

void FeudalFileWriter::checkPoint()
{
    if ( !mWriter.isOpen() )
    {
        cout << "You can't checkpoint a FeudalFileWriter after you've closed "
                "it." << endl;
        CRD::exit(1);
    }
    size_t pos = mWriter.tell();
    finish(pos, false);
    mWriter.seek(pos);
}

void FeudalFileWriter::close()
{
    if ( !mWriter.isOpen() )
        cout << "Warning:  closing feudal file " << mWriter.getFilename()
             << " that is already closed." << endl;
    else
    {
        finish(mWriter.tell());
        mWriter.close();
    }
}

namespace {
void byte_copy( std::string const& filename, BinaryWriter& dest )
{
     BinaryReader source(filename, false);
     unsigned char val;
     while ( !source.atEOF() ) {
          source.read(&val);
          dest.write(val);
     }
}


}

void FeudalFileWriter::finish( size_t pos, bool delete_tmp /* = true */ )
{
    mWriter.flush();
    mWriterOffsets.flush();
    mWriterFixed.flush();

    byte_copy( mWriterOffsets.getFilename(), mWriter );
    byte_copy( mWriterFixed.getFilename(), mWriter );

    if ( delete_tmp ) {
         mWriterOffsets.close();
         mWriterFixed.close();
         Remove( mWriterOffsets.getFilename() );
         Remove( mWriterFixed.getFilename() );
    }

    mWriter.seek(0ul);
    FeudalControlBlock fcb(mElements,
                           pos-sizeof(FeudalControlBlock),
                           mFixedLenDataLen, mVecSize, mEltSize);
    mWriter.write(fcb);
}
