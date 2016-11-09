/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// TextTable:
//      An easier way to print pretty tables

// example:
//    TextTable table;
//    table << DoubleLine ;
//    table << "Name" << Tab << "Age" << Tab << "Weight" << EndRow;
//    table << CharLine('+') ;
//    table << "Jason" << Tab << 3 << Tab << 25.3<< EndRow;
//    table << "Marvin" << Tab << 1 << EndRow;
//    table << DoubleLine ;
//    table.Print( std::cout, 5, "r" );
// output: 
//    =========================
//      Name     Age     Weight
//    +++++++++++++++++++++++++
//     Jason     3       25.3  
//    Marvin     1  
//    =========================
#ifndef UTIL_TEXT_TABLE_H
#define UTIL_TEXT_TABLE_H

#include "String.h"
#include "Vec.h"
#include "feudal/TrackingAllocator.h"

class TextTable {
public:
    // == Special types to draw line and manipulate tables ==
    // A row filled with same char to draw lines
    struct CharLineT { char c; };               
    // Table manipulator
    typedef TextTable& (*TableMan) ( TextTable& ); 

    // Default constructor
    TextTable() { raw_line = false; };

    // == Basic table operations ==
    template< typename T >
    TextTable& operator<< ( const T& t ) 
    {   os << t; return *this;  }

    TextTable& operator<< ( const CharLineT& t ) 
    {   AddCharLine(t.c) ; return *this; }

    TextTable& operator<< ( const TableMan& f ) 
    {   f(*this); return *this; }

    TextTable& AddColumn ()     
    {   line.push_back( os.str() ); os.str(""); return *this;   }

    TextTable& AddRow ()        
    {   AddColumn(); lines.push_back( make_pair( line, raw_line ? '\0': ' ' ) ); 
        line.clear(); raw_line = false; return *this; }

    TextTable& SetRawLine ()        
    {   raw_line = true; return *this; }

    TextTable& AddCharLine (char c) 
    {   lines.push_back( make_pair( std::vector<std::string>(), c ) ); return *this; }

    // Get the vec<vec<String> > copy for external prettyfier
    vec<vec<String> > GetTable() const; 

    void Print( ostream& out, int padding= 1, const string& alignments="", const bool ragged_right = false );

    // Sort -- sort lines by column (zero-based), with a numeric or lexicographical sort, and optionally
    // skip some number of leading lines (e.g. for a header)
    void Sort( size_t column, bool numeric = false, size_t header_skip_size = 0);

private:
    // diable copy and assign constructors
    TextTable(const TextTable&);
    TextTable& operator= (const TextTable&);

    enum RowType { Normal, Line, DoubleLine };
    vec< pair< std::vector<std::string>, char > > lines;
    std::ostringstream os;
    std::vector<std::string> line;
    bool raw_line;
};

// Table manipulators 
template< typename TableT >
TableT& Tab( TableT& tb ) { 
    return tb.AddColumn();
}

template< typename TableT >
TableT& EndRow( TableT& tb ) { 
    return tb.AddRow();
}

// Horizontal lines

inline TextTable::CharLineT CharLine( char c ) {
    TextTable::CharLineT line = { c };
    return line;
}

static const TextTable::CharLineT SingleLine = {'-'};
static const TextTable::CharLineT DoubleLine = {'='};

#endif
