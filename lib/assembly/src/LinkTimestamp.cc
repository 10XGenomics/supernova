///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
//
// MakeDepend: private

#define STRINGIZE(X) #X
#define STRINGIZER(X) STRINGIZE(X)

char const* ARACHNE_RELEASE = "v3.0";
char const* OS_RELEASE = STRINGIZER(UNAME_RELEASE);
char const* LINK_TIMESTAMP = __DATE__ " " __TIME__;
char const* SVN_REVISION = STRINGIZER(SVN_VERSION);

char const* ThisCanBeReadUsingTheStringsProgram =
 "OS Release=" STRINGIZER(UNAME_RELEASE) \
 ", Link Timestamp=" __DATE__ " " __TIME__ \
 ", SVN Revision=" STRINGIZER(SVN_VERSION) \
 ", GCC Version=" __VERSION__;
