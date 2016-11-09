////////////////////////////////////////////////////////////////////////////
//                  SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//      This software and its documentation are copyright (2009) by the   //
//  Broad Institute.  All rights are reserved.  This software is supplied //
//  without any warranty or guaranteed support whatsoever. The Broad      //
//  Institute is not responsible for its use, misuse, or functionality.   //
////////////////////////////////////////////////////////////////////////////
/*
 * \file MakeDepend.cc
 * \author tsharpe
 * \date Sep 15, 2009
 *
 * \brief Examines source dependencies, and writes a makefile.
 */
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <cerrno>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <dirent.h>
#include <sys/stat.h>

using std::cerr;
using std::cout;
using std::endl;
using std::ios;
using std::istream;
using std::ifstream;
using std::map;
using std::ofstream;
using std::ostream;
using std::set;
using std::string;
using std::vector;

char const* DEPENDENCY_DB = "dependencies.db";
char const* DEPENDENCY_DB_TMP = "dependencies.db.tmp";
char const* MAKEFILE = "Makefile_deps";
char const* MAKEFILE_TMP = "Makefile_deps.tmp";

class SourceFile;
class SourceMap;

typedef set<string> StrSet;
typedef set<SourceFile const*> DepsSet;
typedef map<string,SourceFile> SrcMap;
typedef map<string,SourceFile const*> ExeMap;
typedef map<string,StrSet> DirMap;

template <class T>
inline void write( ostream& os, T const& val )
{
    os.write( reinterpret_cast<char const*>(&val), sizeof(T) );
}

inline void write( ostream& os, string const& str )
{
    size_t len = str.size();
    write( os, len );
    if ( len )
        os.write( &str[0], len );
}

inline void write( ostream& os, StrSet const& strSet )
{
    size_t len = strSet.size();
    write( os, len );
    StrSet::const_iterator end(strSet.end());
    for ( StrSet::const_iterator itr(strSet.begin()); itr != end; ++itr )
        write( os, *itr );
}

template <class T>
inline void read( T* pVal, istream& is )
{
    is.read( reinterpret_cast<char*>(pVal), sizeof(T) );
}

inline void read( string* pStr, istream& is )
{
    size_t len = 0;
    read( &len, is );
    pStr->resize(len);
    if ( len )
        is.read( &(*pStr)[0], len );
}

void read( StrSet* pStrList, istream& is )
{
    pStrList->clear();
    StrSet::iterator end(pStrList->end());

    StrSet::size_type sz = 0;
    read( &sz, is );

    string val;

    while ( sz-- )
    {
        read( &val, is );
        end = pStrList->insert(end,val);
    }
}

class SourceFile
{
public:

    SourceFile( string const& path = "", long mTime = 0L );
    SourceFile( SourceFile const& that );
    SourceFile& operator=( SourceFile const& that );
   ~SourceFile();

    string const& getPath() const { return mPath; }
    string getDir() const { return mDirLen ? mPath.substr(0,mDirLen-1) : "."; }
    string asExecutable() const { return mPath.substr(mDirLen,mNameLen); }
    string getSourceRoot() const { return "src/"; }
    string asObjectFile() const
    { auto offset = getSourceRoot().length();
        return mPath.substr(offset,mDirLen+mNameLen-offset) + ".o"; }

    SourceFile const* getAssociatedCC( SourceMap const& ) const;

    DepsSet const& getCompilerDependents( SourceMap const& ) const;
    DepsSet const& getLinkerDependents( SourceMap const& ) const;
    StrSet const& getExecutableDependents( SourceMap const& map ) const
    { if ( !mpExeDeps ) flattenDependents(map); return *mpExeDeps; }
    StrSet const& getLibraryDependents( SourceMap const& map ) const
    { if ( !mpLibDeps ) flattenDependents(map); return *mpLibDeps; }
    StrSet const& getCFlagDependents( SourceMap const& map ) const
    { if ( !mpCFlagDeps ) getCompilerDependents(map); return *mpCFlagDeps; }

    long getMTime() const { return mMTime; }
    bool isHeader() const { return mFlags&IS_HEADER; }
    bool isArchived() const { return mFlags&IS_ARCHIVED; }
    bool isPrivate() const { return mFlags&IS_PRIVATE; }
    bool isExecutable() const { return makesObject() && mFlags&IS_EXECUTABLE; }
    bool makesObject() const { return !(mFlags&(IS_HEADER|IS_ARCHIVED)); }

    bool sameDeps( SourceFile const& that )
    { return mFlags == that.mFlags && mCFlags == that.mCFlags &&
        mLibs == that.mLibs && mExes == that.mExes &&
        mIncludes == that.mIncludes; }

    void parse(); // scan file for its dependencies

    friend void read( SourceFile* pSourceFile, istream& is )
    { read( &pSourceFile->mPath, is );
      pSourceFile->initLengths();
      read( &pSourceFile->mFlags, is );
      read( &pSourceFile->mMTime, is );
      read( &pSourceFile->mCFlags, is );
      read( &pSourceFile->mLibs, is );
      read( &pSourceFile->mExes, is );
      read( &pSourceFile->mIncludes, is );
      pSourceFile->mpAssociatedCC = 0;
      delete pSourceFile->mpCompilerDeps;
      pSourceFile->mpCompilerDeps = 0;
      delete pSourceFile->mpLinkerDeps;
      pSourceFile->mpLinkerDeps = 0; }

    friend void write( ostream& os, SourceFile const& sourceFile )
    { write( os, sourceFile.mPath );
      write( os, sourceFile.mFlags );
      write( os, sourceFile.mMTime );
      write( os, sourceFile.mCFlags );
      write( os, sourceFile.mLibs );
      write( os, sourceFile.mExes );
      write( os, sourceFile.mIncludes ); }

private:
    void initLengths();
    void parseInclude( char* buf );
    void parseDepend( char* buf );
    void flattenDependents( SourceMap const& ) const;
    void appendExeAndLibDependents( SourceFile const* pDep,
                                    SourceMap const& map ) const;

    enum Flags
    {
        IS_NOTHING = 0,
        IS_HEADER = 1,
        IS_EXECUTABLE = 2,
        IS_ARCHIVED = 4,
        IS_PRIVATE = 8
    };

    string mPath; // path relative to work area
    unsigned mDirLen; // length of the directory part of the path
    unsigned mNameLen; // length of the name part of the path
    unsigned mExtLen; // length of the extension part of the path
    mutable unsigned mFlags; // bit-wise from Flags

    long mMTime; // time of last modification


    StrSet mCFlags; // from MakeDepend: cflags
    StrSet mLibs; // libs this depends on from MakeDepend: library
    StrSet mExes; // executables this depends on from MakeDepend: dependency
    StrSet mIncludes; // sources directly included by this one

    mutable SourceFile const* mpAssociatedCC;
    mutable DepsSet* mpCompilerDeps;
    mutable DepsSet* mpLinkerDeps;
    mutable StrSet* mpExeDeps;
    mutable StrSet* mpLibDeps;
    mutable StrSet* mpCFlagDeps;
};

class SourceMap
{
public:
    /// make one by scanning the specified source directory or database
    SourceMap( string const& path );

    // compiled-supplied copying, and dtor are OK

    bool isDirty() const { return mIsDirty; }

    SourceFile const* lookup( string const& path ) const
    { SrcMap::const_iterator itr = mMap.find(path);
      return itr == mMap.end() ? 0 : &itr->second; }

    SrcMap::iterator begin() { return mMap.begin(); }
    SrcMap::iterator end() { return mMap.end(); }
    SrcMap::const_iterator begin() const { return mMap.begin(); }
    SrcMap::const_iterator end() const { return mMap.end(); }
    SrcMap::const_iterator cbegin() const { return mMap.begin(); }
    SrcMap::const_iterator cend() const { return mMap.end(); }

    void append( SourceFile const& sourceFile )
    { mMap[sourceFile.getPath()] = sourceFile; mIsDirty = true;
      updateSubMaps(sourceFile.getPath()); }

    SrcMap::iterator insert( SourceFile const& sourceFile )
    { append(sourceFile); return ++mMap.find(sourceFile.getPath()); }

    SrcMap::iterator erase( SrcMap::iterator& itr );

    ExeMap const& getExes() const { return mExes; }
    DirMap const& getDirs() const { return mDirs; }

    void update( SrcMap::iterator const& itr, SourceFile const& sourceFile )
    { itr->second = sourceFile; mIsDirty = true; }

    friend void read( SourceMap* pMap, istream& is )
    { pMap->mMap.clear();
      SourceFile source;
      SrcMap::size_type sz = 0;
      read( &sz, is );
      while ( sz-- )
      { read( &source, is );
        pMap->mMap[source.getPath()] = source;
        pMap->updateSubMaps(source.getPath()); } }

    friend void write( ostream& os, SourceMap const& map )
    { SrcMap::size_type sz = map.mMap.size();
      write( os, sz );
      SrcMap::const_iterator end(map.mMap.end());
      for ( SrcMap::const_iterator itr(map.mMap.begin()); itr != end; ++itr )
         write( os, itr->second ); }

private:
    void scanDir( string const& path );
    void updateSubMaps( string const& path );

    SrcMap mMap;
    ExeMap mExes; // executable name to source file that has its main()
    DirMap mDirs; // directory names to list of executable names
    bool mIsDirty;
};

SourceFile::SourceFile( string const& path, time_t mTime )
: mPath(path), mFlags(IS_NOTHING), mMTime(mTime),
  mpAssociatedCC(0), mpCompilerDeps(0), mpLinkerDeps(0), mpExeDeps(0),
  mpLibDeps(0), mpCFlagDeps(0)
{
    initLengths();
}

SourceFile::SourceFile( SourceFile const& that )
: mPath(that.mPath), mDirLen(that.mDirLen), mNameLen(that.mNameLen),
  mExtLen(that.mExtLen), mFlags(that.mFlags), mMTime(that.mMTime),
  mCFlags(that.mCFlags), mLibs(that.mLibs), mExes(that.mExes),
  mIncludes(that.mIncludes),
  mpAssociatedCC(0), mpCompilerDeps(0), mpLinkerDeps(0), mpExeDeps(0),
  mpLibDeps(0), mpCFlagDeps(0)
{
}

SourceFile& SourceFile::operator=( SourceFile const& that )
{
    if ( this != &that )
    {
        mPath = that.mPath;
        mDirLen = that.mDirLen;
        mNameLen = that.mNameLen;
        mExtLen = that.mExtLen;
        mFlags = that.mFlags;
        mMTime = that.mMTime;
        mCFlags = that.mCFlags;
        mLibs = that.mLibs;
        mExes = that.mExes;
        mIncludes = that.mIncludes;
        mpAssociatedCC = 0;
        delete mpCompilerDeps;
        mpCompilerDeps = 0;
        delete mpLinkerDeps;
        mpLinkerDeps = 0;
        delete mpExeDeps;
        mpExeDeps = 0;
        delete mpLibDeps;
        mpLibDeps = 0;
        delete mpCFlagDeps;
        mpCFlagDeps = 0;
    }
    return *this;
}

SourceFile::~SourceFile()
{
    delete mpCompilerDeps;
    delete mpLinkerDeps;
    delete mpExeDeps;
    delete mpLibDeps;
    delete mpCFlagDeps;
}

SourceFile const* SourceFile::getAssociatedCC( SourceMap const& map ) const
{
    if ( !mpAssociatedCC && isHeader() )
    {
        string ccFile = mPath.substr(0,mNameLen+mDirLen) + ".cc";
        mpAssociatedCC = map.lookup(ccFile);
    }
    return mpAssociatedCC;
}

DepsSet const& SourceFile::getCompilerDependents( SourceMap const& map ) const
{
    if ( !mpCompilerDeps )
    {
        mpCompilerDeps = new DepsSet();
        mpCFlagDeps = new StrSet(mCFlags);
        vector<StrSet const*> toScan1;
        toScan1.push_back(&mIncludes);
        while ( !toScan1.empty() )
        {
            vector<StrSet const*> toScan2;
            vector<StrSet const*>::iterator end(toScan1.end());
            vector<StrSet const*>::iterator itr(toScan1.begin());
            for ( ; itr != end; ++itr )
            {
                StrSet const& includes = **itr;
                StrSet::const_iterator inclEnd(includes.end());
                StrSet::const_iterator inclItr(includes.begin());
                for ( ; inclItr != inclEnd; ++inclItr )
                {
                    SourceFile const* pDep = map.lookup(getSourceRoot()+*inclItr);
                    if ( !pDep )
                        cerr << "MakeDepend warning: can't find included file "
                             << *inclItr << endl;
                    else if ( pDep != this &&
                            mpCompilerDeps->insert(pDep).second )
                    {
                        mpCFlagDeps->insert(pDep->mCFlags.begin(),
                                            pDep->mCFlags.end());
                        if ( !pDep->mIncludes.empty() )
                            toScan2.push_back(&pDep->mIncludes);
                    }
                }
            }
            toScan1 = toScan2;
        }
    }
    return *mpCompilerDeps;
}

DepsSet const& SourceFile::getLinkerDependents( SourceMap const& map ) const
{
    if ( !mpLinkerDeps )
    {
        mpLinkerDeps = new DepsSet();
        vector<SourceFile const*> toScan1;
        toScan1.push_back(this);
        while ( toScan1.size() )
        {
            vector<SourceFile const*> toScan2;
            vector<SourceFile const*>::iterator end(toScan1.end());
            vector<SourceFile const*>::iterator itr(toScan1.begin());
            for ( ; itr != end; ++itr )
            {
                DepsSet const& deps = (*itr)->getCompilerDependents(map);
                DepsSet::const_iterator depsEnd(deps.end());
                DepsSet::const_iterator depsItr(deps.begin());
                for ( ; depsItr != depsEnd; ++depsItr )
                {
                    SourceFile const* pDep = (*depsItr)->getAssociatedCC(map);
                    if ( pDep && pDep != this &&
                            mpLinkerDeps->insert(pDep).second )
                        toScan2.push_back(pDep);
                }
            }
            toScan1 = toScan2;
        }
    }
    return *mpLinkerDeps;
}

void SourceFile::parse()
{
    ifstream in(mPath.c_str(),ios::in|ios::binary);
    char buf[4096];
    while ( in.getline(buf,sizeof(buf)) )
    {
        if ( !memcmp(buf,"int main(",9) )
            mFlags |= IS_EXECUTABLE;
        else if ( !memcmp(buf,"#include",8) )
            parseInclude(buf);
        else if ( !memcmp(buf,"// MakeDepend:",14) )
            parseDepend(buf);
    }
    if ( !in.eof() )
        cerr << "MakeDepend warning:  failed to parse source file " <<
                mPath << endl;
    in.close();
}

void SourceFile::initLengths()
{
    char const* ppp = strrchr(mPath.c_str(), '/');
    mDirLen = ppp ? ppp - mPath.c_str() + 1 : 0;
    ppp = strrchr(mPath.c_str() + mDirLen, '.');
    mExtLen = ppp ? mPath.c_str() + mPath.size() - ppp : 0;
    mNameLen = mPath.size() - mDirLen - mExtLen;
    if ( ppp && memcmp(ppp,".cc",4) )
        mFlags = IS_HEADER;
}

void SourceFile::parseInclude( char* buf )
{
    bool ok = false;
    char const* originalBuf = buf;
    buf += 8;

    while ( isblank(*buf) ) ++buf;

    if ( *buf == '<' )
    {
        ok = true;
    }
    else if ( *buf == '"' )
    {
        char* ppp = strchr(++buf,'"');
        if ( ppp )
        {
            ok = true;
            *ppp = 0;
            mIncludes.insert(buf);
        }
    }

    if ( !ok )
    {
        cerr << "MakeDepend warning: failed to comprehend the line '" <<
                originalBuf << "'" << endl;
    }
}

void SourceFile::parseDepend( char* buf )
{
    char const* originalBuf = buf;
    buf += 14;

    while ( isblank(*buf) ) ++buf;

    if ( !memcmp(buf,"archived",8) )
        mFlags |= IS_ARCHIVED;
    else if ( !memcmp(buf,"private",7) )
        mFlags |= IS_PRIVATE;
    else if ( !memcmp(buf,"cflags ",7) )
        mCFlags.insert(buf+7);
    else if ( !memcmp(buf,"library ",8) )
        mLibs.insert(buf+8);
    else if ( !memcmp(buf,"dependency ",11) )
    {
        buf += 11;
        while ( true )
        {
            char chr;
            while ( isblank((chr = *buf)) ) ++buf;
            if ( !chr )
                break;

            char* ppp = buf + 1;
            while ( (chr = *ppp) && !isblank(chr) ) ++ppp;
            if ( !chr )
            {
                mExes.insert(buf);
                break;
            }
            *ppp++ = 0;
            mExes.insert(buf);
            buf = ppp;
        }
    }
    else if ( !memcmp(buf,"shared",6) )
    {
        // ignoring this feature
    }
    else
    {
        cerr << "MakeDepend warning: failed to comprehend the line '" <<
                originalBuf << "'" << endl;
    }
}

void SourceFile::flattenDependents( SourceMap const& map ) const
{
    mpExeDeps = new StrSet();
    mpLibDeps = new StrSet();

    appendExeAndLibDependents(this,map);

    DepsSet const& cDeps = getLinkerDependents(map);
    DepsSet::const_iterator end(cDeps.end());
    for ( DepsSet::const_iterator itr(cDeps.begin()); itr != end; ++itr )
        appendExeAndLibDependents(*itr,map);
}

void SourceFile::appendExeAndLibDependents( SourceFile const* pDep,
                                            SourceMap const& map ) const
{
    mpExeDeps->insert(pDep->mExes.begin(),pDep->mExes.end());
    mpLibDeps->insert(pDep->mLibs.begin(),pDep->mLibs.end());

    DepsSet const& hDeps = pDep->getCompilerDependents(map);
    DepsSet::const_iterator end(hDeps.end());
    for ( DepsSet::const_iterator itr(hDeps.begin()); itr != end; ++itr )
    {
        pDep = *itr;
        mpExeDeps->insert(pDep->mExes.begin(),pDep->mExes.end());
        mpLibDeps->insert(pDep->mLibs.begin(),pDep->mLibs.end());
    }
}

SourceMap::SourceMap( string const& path )
: mIsDirty(false)
{
    struct stat statbuf;
    if ( !stat(path.c_str(),&statbuf) )
    {
        if ( S_ISDIR(statbuf.st_mode) )
            scanDir(path);
        else
        {
            ifstream is(path.c_str());
            size_t hdr = 0;
            read( &hdr, is );
            read( this, is );
            if ( !is )
            {
                cerr << "MakeDepend error: Unable to read dependency db."
                     << endl;
                ::exit(1);
            }
        }
    }
    // else can't stat the file -- must be a new db
}

void SourceMap::scanDir( string const& path )
{
    DIR* dirStream = opendir(path.c_str());
    if ( !dirStream )
    {
        cerr << "MakeDepend error: Unable to open the directory " <<
                path << " to scan for sources.  Exiting." << endl;
        ::exit(1);
    }

    string dirPath = path + "/";
    if ( dirPath == "./" )
        dirPath.clear();

    struct dirent* pDirEnt;
    while ( (pDirEnt = readdir(dirStream)) )
    {
        string filename(dirPath+pDirEnt->d_name);
        struct stat statbuf;
        int result = lstat(filename.c_str(),&statbuf);
        if ( result )
        {
            int err = errno;
            cerr << "MakeDepend warning: Unable to stat "
                    << filename << ": " << strerror(err) << " [errno=" << err
                    << "]." << endl;
        }
        else if ( S_ISDIR(statbuf.st_mode) )
        {
            if ( pDirEnt->d_name[0] != '.' )
                scanDir( filename );
        }
        else if ( S_ISREG(statbuf.st_mode) )
        {
            size_t len = strlen(pDirEnt->d_name);
            char const* end = pDirEnt->d_name + len;
            if ( (len >= 2 && !memcmp(end-2,".h",2)) ||
                 (len >= 3 && !memcmp(end-3,".cc",3)) )
            {
                SourceFile sourceFile(filename,statbuf.st_mtime);
                mMap[sourceFile.getPath()] = sourceFile;
            }
        }
    }
    closedir(dirStream);
}

SrcMap::iterator SourceMap::erase( SrcMap::iterator& itr )
{
    SourceFile& file = itr->second;
    if ( file.isExecutable() )
    {
        string exeName = file.asExecutable();
        mExes.erase(exeName);
        string dirName = file.getDir();
        StrSet& exesInDir = mDirs[dirName];
        exesInDir.erase(exeName);
        if ( exesInDir.empty() )
            mDirs.erase(dirName);
    }
    mMap.erase(itr++);
    mIsDirty = true;
    return itr;
}

void SourceMap::updateSubMaps( string const& path )
{
    SourceFile const& file = mMap[path];
    if ( file.isExecutable() )
    {
        string exeName = file.asExecutable();

        ExeMap::iterator itr = mExes.find(exeName);
        if ( itr == mExes.end() )
        {
            mExes[exeName] = &file;
            mDirs[file.getDir()].insert(exeName);
        }
        else
        {
            cerr << "MakeDepend error: sources " << itr->second->getPath()
                    << " and " << file.getPath()
                    << " both produce executables with the name " << exeName
                    << endl;
        }
    }
}

void generateExeBits( ostream& os, SourceMap const& db, SourceFile const& file )
{
    // linker dependencies
    string exeName = file.asExecutable();
    os << "\n$(BIN)/" << exeName << ": $(OBJ)/" << file.asObjectFile();
    DepsSet const& lDeps = file.getLinkerDependents(db);
    DepsSet::const_iterator lEnd(lDeps.end());
    for ( DepsSet::const_iterator lItr(lDeps.begin()); lItr != lEnd; ++lItr )
        os << " $(OBJ)/" << (*lItr)->asObjectFile();
    os << " | $(LIBNAM)\n\t$(LINK)\n";

    StrSet const& libDeps = file.getLibraryDependents(db);
    if ( libDeps.size() )
    {
        os << "\n$(BIN)/" << exeName << ": LDLIBS := $(LDLIBS)";
        StrSet::const_iterator end(libDeps.end());
        for ( StrSet::const_iterator itr(libDeps.begin()); itr != end; ++itr )
            os << " $(" << *itr << ')';
        os << '\n';
    }

    // inter-program dependencies
    os << "\n.PHONY: " << exeName << '\n' << exeName << ": $(BIN)/" << exeName;
    StrSet const& exes = file.getExecutableDependents(db);
    if ( exes.size() )
    {
        StrSet::const_iterator end(exes.end());
        for ( StrSet::const_iterator itr(exes.begin()); itr != end; ++itr )
            os << ' ' << *itr;
    }
    os << '\n';
}

bool isSubDir( string const& sub, string const& dir )
{
    string::size_type len = dir.size();
    return sub.size() > len && sub.substr(0,len)==dir && sub[len]=='/';
}

void generateMakefile( ostream& os, SourceMap const& db )
{
    os << "LIBNAM := $(OBJ)/lib10X.a\n";

    DepsSet libObjects;

    // for each compilation unit (*.cc)
    SrcMap::const_iterator sEnd(db.end());
    for ( SrcMap::const_iterator sItr(db.begin()); sItr != sEnd; ++sItr )
    {
        SourceFile const& file = sItr->second;
        if ( file.makesObject() )
        {
            string objName = "$(OBJ)/" + file.asObjectFile();

            // compiler dependencies
            os << '\n' << objName << ": " << file.getPath();
            DepsSet const& cDeps = file.getCompilerDependents(db);
            typedef DepsSet::const_iterator CItr;
            CItr cEnd(cDeps.end());
            for ( CItr cItr(cDeps.begin()); cItr != cEnd; ++cItr )
                os << ' ' << (*cItr)->getPath();
            os << "\n\t$(COMPILE)\n";

            // special compilation flags
            StrSet const& cflagDeps = file.getCFlagDependents(db);
            if ( !cflagDeps.empty() )
            {
                os << '\n' << objName << ": CFLAGS := $(CFLAGS)";
                typedef StrSet::const_iterator SItr;
                SItr fEnd(cflagDeps.end());
                for ( SItr fItr(cflagDeps.begin()); fItr != fEnd; ++fItr )
                    os << " $(" << *fItr << ')';
                os << '\n';
            }

            // for cc files with a main
            if ( file.isExecutable() )
            {
                generateExeBits(os,db,file);
                DepsSet const& deps = file.getLinkerDependents(db);
                libObjects.insert(deps.begin(),deps.end());
            } else {
            	// WARNING!!!!
            	// EVERYTHING that's not an executable, unless marked private, NOW gets tossed
            	// into the lib10X.a
            	//
            	// that's to facilitate the inclusion of code that will ONLY
            	// be used externally -- a big change from previous behavior
            	if (!file.isPrivate()) libObjects.insert(&file);
            }
        }
    }

    // for each dir with an executable in it
    DirMap::const_iterator dEnd(db.getDirs().end());
    DirMap::const_iterator dItr(db.getDirs().begin());
    for ( ; dItr != dEnd; ++dItr )
    {
        os << "\n.PHONY: " << dItr->first << '\n' << dItr->first << ':';
        StrSet::iterator dnEnd(dItr->second.end());
        StrSet::iterator dnItr(dItr->second.begin());
        for ( ; dnItr != dnEnd; ++dnItr )
            os << ' ' << *dnItr;
        auto dItr2 = dItr;
        while ( ++dItr2 != dEnd && isSubDir(dItr2->first,dItr->first) )
            os << ' ' << dItr2->first;
        os << '\n';
    }

    os << "\nLIBOBJECTS :=";
    DepsSet::iterator oEnd(libObjects.end());
    for ( DepsSet::iterator oItr(libObjects.begin()); oItr != oEnd; ++oItr ) {
        os << " $(OBJ)/" << (*oItr)->asObjectFile();
    }

    os << "\n$(LIBNAM):\n\t$(ARCHIVE)\n";

    os << "\nMAINOBJECTS :=";
    ExeMap const& exeMap = db.getExes();
    ExeMap::const_iterator eEnd(exeMap.end());
    for ( ExeMap::const_iterator eItr(exeMap.begin()); eItr != eEnd; ++eItr )
        os << " $(OBJ)/" << eItr->second->asObjectFile();

    os << "\n.PHONY: objects\nobjects: $(LIBOBJECTS) $(MAINOBJECTS)\n.PHONY: all\nall:";
    for ( ExeMap::const_iterator eItr(exeMap.begin()); eItr != eEnd; ++eItr )
        os << " $(BIN)/" << eItr->first;
    os << '\n';
}

void doUpdate( SourceMap& db, bool writeDB )
{
    SourceMap src(".");
    SrcMap::iterator dbItr(db.begin());
    SrcMap::iterator dbEnd(db.end());
    SrcMap::iterator srcItr(src.begin());
    SrcMap::iterator srcEnd(src.end());
    bool writeMakefile = false;
    while ( srcItr != srcEnd )
    {
        SourceFile& source = srcItr->second;
        if ( dbItr == dbEnd )
        {
            source.parse();
            db.append(source);
            ++srcItr;
            writeMakefile = true;
        }
        else
        {
            SourceFile& dbSource = dbItr->second;
            if ( source.getPath() == dbSource.getPath() )
            {
                if ( source.getMTime() > dbSource.getMTime() )
                {
                    source.parse();
                    if ( !source.sameDeps(dbSource) )
                        writeMakefile = true;
                    db.update(dbItr,source);
                }
                ++dbItr;
                ++srcItr;
            }
            else if ( source.getPath() < dbSource.getPath() )
            {
                source.parse();
                dbItr = db.insert(source);
                ++srcItr;
                writeMakefile = true;
            }
            else
            {
                dbItr = db.erase(dbItr);
                writeMakefile = true;
            }
        }
    }

    if ( writeDB && db.isDirty() )
    {
        ofstream os(DEPENDENCY_DB_TMP,ios::binary|ios::out|ios::trunc);
        os.write("BINWRITE",8);
        write( os, db );
        os.close();
        if ( !os )
        {
            cerr << "MakeDepend error: Unable to write dependency db.  Exiting."
                 << endl;
            ::exit(1);
        }
        rename(DEPENDENCY_DB_TMP,DEPENDENCY_DB);
    }

    struct stat sb;
    if ( writeMakefile || stat(MAKEFILE,&sb) )
    {
        ofstream os(MAKEFILE_TMP,ios::out|ios::trunc);
        generateMakefile(os,db);
        os.close();
        if ( !os )
        {
            cerr << "MakeDepend error: Unable to write makefile.  Exiting."
                 << endl;
            ::exit(1);
        }
        rename(MAKEFILE_TMP,MAKEFILE);
    }
}

void checkUnused( SourceMap const& db )
{
    typedef DepsSet::const_iterator DepsItr;
    DepsSet allSources;
    SrcMap::const_iterator dEnd(db.end());
    for ( SrcMap::const_iterator dItr(db.begin()); dItr != dEnd; ++dItr )
        allSources.insert(&dItr->second);

    DepsSet usedSources;
    ExeMap const& exeMap = db.getExes();
    ExeMap::const_iterator eEnd(exeMap.end());
    for ( ExeMap::const_iterator eItr(exeMap.begin()); eItr != eEnd; ++eItr )
    {
        SourceFile const* pSrc = eItr->second;
        usedSources.insert(pSrc);
        DepsSet const& cDeps = pSrc->getCompilerDependents(db);
        usedSources.insert(cDeps.begin(),cDeps.end());
        DepsSet const& lDeps = pSrc->getLinkerDependents(db);
        for ( DepsItr lItr(lDeps.begin()), lEnd(lDeps.end());lItr!=lEnd;++lItr )
        {
            pSrc = *lItr;
            usedSources.insert(pSrc);
            DepsSet const& clDeps = pSrc->getCompilerDependents(db);
            usedSources.insert(clDeps.begin(),clDeps.end());
        }
    }

    StrSet unusedSources;
    DepsItr uEnd(usedSources.end());
    DepsItr aEnd(allSources.end());
    for ( DepsItr aItr(allSources.begin()); aItr != aEnd; ++aItr )
        if ( usedSources.find(*aItr) == uEnd )
            unusedSources.insert((*aItr)->getPath());

    if ( !unusedSources.empty() )
    {
        cerr << "MakeDepend warning:  The following source files are unused:"
                << endl;
        typedef StrSet::iterator StrsItr;
        StrsItr uuEnd(unusedSources.end());
        for ( StrsItr uuItr(unusedSources.begin()); uuItr != uuEnd; ++uuItr )
            cerr << '\t' << *uuItr << endl;
    }
}

void dumpModules( SourceMap const& db )
{
    typedef ExeMap::const_iterator Itr;
    ExeMap const& exes = db.getExes();
    for ( Itr itr(exes.begin()), end(exes.end()); itr != end; ++itr )
    {
        std::cout << itr->first << ":\n";
        StrSet list;
        list.insert("LinkTimestamp");
        list.insert(itr->first);
        DepsSet const& deps = itr->second->getLinkerDependents(db);
        typedef DepsSet::const_iterator DItr;
        for ( DItr dItr(deps.begin()), dEnd(deps.end()); dItr != dEnd; ++dItr )
            list.insert((*dItr)->asExecutable());
        typedef StrSet::iterator SItr;
        for ( SItr sItr(list.begin()), sEnd(list.end()); sItr != sEnd; ++sItr )
            std::cout << *sItr << std::endl;
    }
}

int main( int argc, char** argv )
{
    bool writeDB = true;
    if ( argc > 1 && !strcmp(argv[1],"-n") )
    {
        --argc; ++argv;
        writeDB = false;
    }

    SourceMap db(DEPENDENCY_DB);
    doUpdate(db,writeDB);

    if ( argc > 1 && !strcmp(argv[1],"-u") )
    {
        --argc; ++argv;
        checkUnused(db);
    }
    if ( argc > 1 && !strcmp(argv[1],"-x") )
    {
        --argc; ++argv;
        dumpModules(db);
    }

    ExeMap const& exeMap = db.getExes();
    StrSet exes;
    for ( int idx = 1; idx < argc; ++idx )
    {
        ExeMap::const_iterator itr = exeMap.find(argv[idx]);
        if ( itr == exeMap.end() )
            cerr << "MakeDepend warning: Ignoring parts-list request for " <<
                    argv[idx] << " that I don't know how to build." << endl;
        else
        {
            exes.insert(argv[idx]);
            StrSet const& subs = itr->second->getExecutableDependents(db);
            exes.insert(subs.begin(), subs.end());
        }
    }

    DepsSet deps;
    StrSet::iterator eEnd(exes.end());
    for ( StrSet::iterator eItr(exes.begin()); eItr != eEnd; ++eItr )
    {
        ExeMap::const_iterator itr = exeMap.find(*eItr);
        if ( itr == exeMap.end() )
        {
            cerr << "MakeDepend warning: Ignoring parts-list component " <<
                    *eItr << " that I don't know how to build." << endl;
            continue;
        }
        SourceFile const* pFile = itr->second;
        deps.insert(pFile);
        DepsSet const& lDeps = pFile->getLinkerDependents(db);
        deps.insert(lDeps.begin(), lDeps.end());
        DepsSet::const_iterator lEnd(lDeps.end());
        DepsSet::const_iterator lItr(lDeps.begin());
        for ( ; lItr != lEnd; ++lItr )
        {
            DepsSet const& cDeps = (*lItr)->getCompilerDependents(db);
            deps.insert(cDeps.begin(), cDeps.end());
        }
    }

    int err = 0;
    DepsSet::iterator end(deps.end());
    for ( DepsSet::iterator itr(deps.begin()); itr != end; ++itr )
    {
        cout << (*itr)->getPath();
        if ( (*itr)->isArchived() )
        {
            cout << " ***ARCHIVED***";
            err = 1;
        }
        if ( (*itr)->isPrivate() )
        {
            cout << " ***PRIVATE***";
            err = 1;
        }
        cout << '\n';
    }
    return err;
}
