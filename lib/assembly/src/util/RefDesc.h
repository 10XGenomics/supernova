/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////
/*
 * \file RefDesc.h
 * \author tsharpe
 * \date Dec 11, 2008
 *
 * \brief Metadata for reference sequence.
 */
#ifndef REFDESC_H_
#define REFDESC_H_

#include "Vec.h"
#include <string>
#include <iostream>

/**
 * \class RefDesc
 * \brief A packet of metadata describing reference sequence.
 */
class RefDesc
{
public:
    typedef unsigned long ulong;
    typedef std::string string;
    typedef vec<RefDesc> RefDescVec;
    typedef std::istream istream;
    typedef std::ostream ostream;

    /// Parses a SAM-header dictionary line to produce a RefDesc
    RefDesc( string const& samSQLine );

    RefDesc( string const& name, unsigned id, string const& uri, ulong length,
                string const& hash, string const& species, string const& asmID )
     : mName(name), mId(id), mURI(uri), mLength(length), mHash(hash),
       mSpecies(species), mAsmID(asmID)
    {}

    RefDesc() : mId(~0u), mLength(0) {}

    /// The name of the sequence.  E.g., chr18, or "whole genome".
    string const& getName() const
    { return mName; }

    unsigned getId() const
    { return mId; }

    string const& getURI() const
    { return mURI; }

    ulong getLength() const
    { return mLength; }

    string const& getHash() const
    { return mHash; }

    string const& getSpecies() const
    { return mSpecies; }

    string const& getAssemblyID() const
    { return mAsmID; }

    /// Read a FASTA file to build the dictionary
    static RefDescVec readFASTA( string const& fastaFileName );

    /// Read a dict file to build the dictionary
    static RefDescVec readDict( string const& dictFileName );

    /// Read a dictionary from a stream
    static RefDescVec readDict( istream& is );

    /// Write a dictionary to a file
    static bool writeDict( string const& dictFileName,
                                RefDescVec const& refDescs,
                                bool noHD = true );

    /// Write a dictionary to a stream
    static void writeDict( ostream& os, RefDescVec const& refDescs, bool noHD = false );

private:
    static void readDict( istream& is, RefDescVec& vec );

    string mName;
    unsigned mId;
    string mURI;
    ulong mLength;
    string mHash;
    string mSpecies;
    string mAsmID;
};

typedef RefDesc::RefDescVec RefDict;

#endif /* REFDESC_H_ */
