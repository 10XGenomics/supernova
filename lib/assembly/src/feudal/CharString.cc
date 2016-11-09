///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * \file CharString.cc
 * \author ghall
 */
#include "feudal/CharString.h"
#include "Charvector.h"
#include "system/ErrNo.h"
#include <cmath>
#include <cstring>

using std::cout;
using std::endl;

/*
 * External Broad Requirements
 */

// ToStringAbbrev
String ToStringAbbrev(const long x)
{ if (abs(x) < 1000) return ToString(x);
  double xDouble = x;
  xDouble /= 1000.;
  String suffix = "K";
  unsigned int precision = 0;
  if (fabs(xDouble) > 1000.)
  { xDouble /= 1000.;
    suffix = "M";
    ++precision; }
  if (fabs(xDouble) > 1000.)
  { xDouble /= 1000.;
    suffix = "G";
    ++precision; }
  return ToString(xDouble, precision) + suffix; }

// ToStringAddCommas

String ToStringAddCommas(const int64_t x)
{ String tmps = x < 0 ? ToString(-x) : ToString(x);
  String result;
  result.reserve(tmps.size()+tmps.size()/3);
  String::const_reverse_iterator end(tmps.crend());
  unsigned int phase = 0;
  for ( String::const_reverse_iterator itr(tmps.crbegin()); itr != end; ++itr )
  {
      if ( phase == 3 ) { result.push_back(','); phase = 0; }
      phase += 1;
      result.push_back(*itr);
  }
  if ( x < 0 ) result.push_back( '-' );
  using std::reverse;
  reverse(result.begin(),result.end());
  return result; }

String ToStringAddCommas(const uint64_t x)
{ String tmps = ToString(x);
  String result;
  result.reserve(tmps.size()+tmps.size()/3);
  String::const_reverse_iterator end(tmps.crend());
  unsigned int phase = 0;
  for ( String::const_reverse_iterator itr(tmps.crbegin()); itr != end; ++itr )
  {
      if ( phase == 3 ) { result.push_back(','); phase = 0; }
      phase += 1;
      result.push_back(*itr);
  }
  using std::reverse;
  reverse(result.begin(),result.end());
  return result; }

String ToStringAddCommas(const int32_t x)
{ String tmps = x < 0 ? ToString(-x) : ToString(x);
  String result;
  result.reserve(tmps.size()+tmps.size()/3);
  String::const_reverse_iterator end(tmps.crend());
  unsigned int phase = 0;
  for ( String::const_reverse_iterator itr(tmps.crbegin()); itr != end; ++itr )
  {
      if ( phase == 3 ) { result.push_back(','); phase = 0; }
      phase += 1;
      result.push_back(*itr);
  }
  if ( x < 0 ) result.push_back( '-' );
  using std::reverse;
  reverse(result.begin(),result.end());
  return result; }

String ToStringAddCommas(const uint32_t x)
{ String tmps = ToString(x);
  String result;
  result.reserve(tmps.size()+tmps.size()/3);
  String::const_reverse_iterator end(tmps.crend());
  unsigned int phase = 0;
  for ( String::const_reverse_iterator itr(tmps.crbegin()); itr != end; ++itr )
  {
      if ( phase == 3 ) { result.push_back(','); phase = 0; }
      phase += 1;
      result.push_back(*itr);
  }
  using std::reverse;
  reverse(result.begin(),result.end());
  return result; }

// RemoveCommas - useful as the sort-of inverse of above (returns a String).
String RemoveCommas(const String& s)
{
    String result;
    char const comma = ',';
    auto not_comma = [](char c)->bool { return c != comma; };
    std::back_insert_iterator<String> back_itr(result);
    std::copy_if( s.begin(), s.end(), back_itr, not_comma );
    return result;
}

// QuoteString
String QuoteString(const String& s)
{
    String quoted;
    quoted.push_back('"');
    for ( auto ch: s ) {
	if ( ch == '"' ) quoted.push_back('\\');
	quoted.push_back(ch);
    }
    quoted.push_back('"');
    return quoted;
}

// ToLower
String ToLower(const String &s)
{ String slower(s); slower.ToLower(); return slower; }

// ToUpper
String ToUpper(const String & s)
{ String supper(s); supper.ToUpper(); return supper; }

// DeleteLeadingWhiteSpace
void DeleteLeadingWhiteSpace(String& s)
{ size_type lastpos = 0;
  while((lastpos < s.size()) && isspace(s[lastpos]))
    ++lastpos;
  s.erase(0, lastpos); }

// DeleteTrailingWhiteSpace
void DeleteTrailingWhiteSpace(String& s)
{ size_type firstpos = s.size();
  while((firstpos > 0) && isspace(s[firstpos - 1]))
    --firstpos;
  s.erase(firstpos); }

// WhiteSpaceFree
String WhiteSpaceFree(const String& s)
{ String answer;
  answer.reserve(s.size());
  String::const_iterator end(s.end());
  for ( String::const_iterator itr(s.begin()); itr != end; ++itr )
      if ( !isspace(*itr) ) answer.push_back(*itr);
  return answer; }

// cmp_numeric
bool cmp_numeric(const String& s1, const String& s2)
{ for (String::size_type i = 0; i < s1.size( ); i++)
  { if (i == s2.size()) return false;
    if (!isdigit(s1[i]) || !isdigit(s2[i]))
    { if (s1[i] < s2[i]) return true;
      if (s1[i] > s2[i]) return false;
      continue; }
    String::size_type j1, j2;
    for (j1 = i + 1; j1 < s1.size(); j1++)
      if (!isdigit(s1[j1])) break;
    for (j2 = i + 1; j2 < s2.size(); j2++)
      if (!isdigit(s2[j2])) break;
    long n1 = 0, n2 = 0;
    for (String::size_type k = i; k < j1; k++)
      n1 = (10 * n1) + (s1[k] - '0');
    for (String::size_type k = i; k < j2; k++)
      n2 = (10 * n2) + (s2[k] - '0');
    if (n1 < n2) return true;
    if (n1 > n2) return false;
    if (j1 < j2) return false;
    if (j1 > j2) return true;
    i = j1 - 1; }
    return s1.size() < s2.size(); }

String BaseAlpha(unsigned int n)
{ String result;
  result.reserve(7);
  do
  {
      result.push_back('A' + n%26);
  }
  while ( n /= 26 );
  using std::reverse;
  reverse(result.begin(),result.end());
  return result; }

unsigned int UnBaseAlpha(const String& s)
{ unsigned int answer = 0;
  String::const_iterator end(s.end());
  for ( String::const_iterator itr(s.begin()); itr != end; ++itr )
      answer = 26*answer + (*itr - 'A');
  return answer; }

String BaseAlNum(unsigned int n)
{    String result;
     result.reserve(7);
     do
     {    int k = n % 36;
          if ( k < 26 ) result.push_back( 'A' + k );
          else result.push_back( '0' + (k-26) );    }
     while ( n /= 36 );
     reverse( result.begin( ),result.end( ) );
     return result;    }

std::istream& operator>> (std::istream& in, String& s)
{ std::string buf; in >> buf; s = buf.c_str(); return in; }

#include "feudal/FeudalStringDefs.h"
template class FeudalString< char >;
