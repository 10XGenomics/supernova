/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////
#include "util/TextTable.h"
#include "String.h"
#include "Vec.h"
#include <limits>
#include <sstream>
#include <iterator>

vec<vec<String> > TextTable::GetTable( ) const {
    vec<vec<String> > table;
    for ( size_t i = 0; i < lines.size(); ++i ) {
        table.push_back( vec<String>() );
        for ( size_t j = 0; j < lines[i].first.size(); ++j ) 
            table.back().push_back( lines[i].first[j] );
    }
    return table;
}


void TextTable::Print( ostream& out, int padding, const string& alignments, const bool ragged_right ) {
    // determine the actual width for each column
    size_t ncols = 0;
    for ( size_t i = 0; i < lines.size(); ++i ) 
        if (lines[i].second != '\0')
            ncols = max( ncols, lines[i].first.size() );
    vec<size_t> widths( ncols );
    for ( size_t i = 0; i < lines.size(); ++i ) 
        if (lines[i].second != '\0')
            for ( size_t j = 0; j < lines[i].first.size(); ++j ) 
                widths[j] = max( widths[j], lines[i].first[j].size() );
    // total width of the table
    int total_width = 0;
    for ( size_t i = 0; i < widths.size(); ++i ) 
        total_width += widths[i] + padding;
    total_width -= padding;
    for ( size_t i = 0; i < lines.size(); ++i ) {
        if ( lines[i].second == '\0' ) {
            for (size_t j = 0; j < lines[i].first.size(); j++) 
                out << lines[i].first[j];
            out << endl;
        }
        else if ( lines[i].second != ' ' ) {
            out << string( total_width, lines[i].second ) << endl;
        }
        else {
            for ( size_t j = 0; j < lines[i].first.size(); ++j ) {
                char align = 'l';
                if ( j < alignments.size() ) align = alignments[j];
                if ( align == 'l' ) {
                    out << lines[i].first[j];
                    if ( j < lines[i].first.size()-1 || !ragged_right )		// don't pad out the last col if ragged_right
                	out << string( widths[j] - lines[i].first[j].size(), ' ' );
                }
                else if ( align == 'r' ) {
                    out << string( widths[j] - lines[i].first[j].size(), ' ' ) << lines[i].first[j] ;
                }
                if ( j < lines[i].first.size() -1 &&  padding > 0 )
                    out << string(padding,' ');
            }
            out << endl;
        }
    }
}


// Sort -- sort lines by column (zero-based), with a numeric or lexicographical sort, and
// optionally skip some number of leading lines (e.g. for a header)
//
// only true columnar lines (.second == ' ') with a value in the specified column are touched.
//
// inefficiency is that we duplicate the affected lines to sort -- we could just permute
// the lines in-place, but this is (slightly) complicated by the fact that some lines are
// pinned in place.  This works for now.
//
void TextTable::Sort( size_t col, bool numeric, size_t header_skip_size )
{
    typedef std::vector<std::string> SortLine;
    std::vector<SortLine> sort_lines;

    // populate the vector to sort
    for ( size_t i = header_skip_size; i < lines.size(); ++i ) {
	if ( lines[i].second == ' ' && col < lines[i].first.size() ) {
	    sort_lines.push_back( lines[i].first );
	}
    }

    // sort the vector
    if ( numeric )
	std::stable_sort( sort_lines.begin(), sort_lines.end(),
		[col] (SortLine const x, SortLine const y)  {
		    std::istringstream sx(x[col]), sy(y[col]);
		    long double dx, dy;
		    if ( !(sx >> dx) ) dx = -std::numeric_limits<long double>::max();
		    if ( !(sy >> dy) ) dy = -std::numeric_limits<long double>::max();
		    return dy < dx;		// note: reverse sort
		});
    else
	std::stable_sort( sort_lines.begin(), sort_lines.end(),
		[col] (SortLine const x, SortLine const y)  {
		    return y[col] < x[col];	// note: reverse sort
		});

    // replace the true vector
    for ( size_t i = header_skip_size; i < lines.size(); ++i ) {
	if ( lines[i].second == ' ' && col < lines[i].first.size() ) {
	    std::swap(lines[i].first, sort_lines.back());
	    sort_lines.pop_back();
	}
    }
}
