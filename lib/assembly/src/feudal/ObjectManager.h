///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * ObjectManager.h
 *
 *  Created on: Sep 16, 2014
 *      Author: tsharpe
 *
 * ObjectManager pairs an object with an on-disk representation and manages
 *     Feudal versus plain binary stream storage.
 *
 * load()   - returns a pointer to the desired object, loading off disk if necessary
 * unload() - frees the in-memory representation, it will be reloaded from disk at
 *            the next load()
 * store()  - push data out to disk
 * create() - use to allocate an object when no on-disk representation exists (yet)
 *
 */
#ifndef FEUDAL_OBJECTMANAGER_H_
#define FEUDAL_OBJECTMANAGER_H_

#include "feudal/BinaryStream.h"
#include "feudal/MasterVec.h"
#include "system/Assert.h"
#include "system/file/File.h"

template <class T>
class ObjectManagerImpl
{
protected:
    void loadImpl( File const& file, T* pObj )
    { BinaryReader::readFile(file,pObj); }
    void storeImpl( File const& file, T const& obj )
    { BinaryWriter::writeFile(file,obj); }
};

template <class T>
class ObjectManagerImpl<MasterVec<T>>
{
protected:
    void loadImpl( File const& file, MasterVec<T>* pObj )
    { pObj->ReadAll(String(file.c_str())); }
    void storeImpl( File const& file, MasterVec<T> const& obj )
    { obj.WriteAll(String(file.c_str())); }
};

template <class T>
class ObjectManager : ObjectManagerImpl<T>
{
public:
    ObjectManager( String const& fileName ) : mFile(fileName), mpObj(nullptr) {}
    ObjectManager( ObjectManager const& )=delete;
    ObjectManager& operator=( ObjectManager const& )=delete;
    ~ObjectManager() { delete mpObj; }

    T& create() { delete mpObj; mpObj = new T; return *mpObj; }
    T& object() { return *mpObj; }

    T const& load() { return this->load_mutable(); }

    T & load_mutable() {
         if ( !mpObj ) {
              cout << "load: loading object from " << filename() << endl;
              mpObj=new T; this->loadImpl(mFile,mpObj);
         }
         return *mpObj;
    }

    void unload() { delete mpObj; mpObj = nullptr; }
    void store() { ForceAssert(mpObj); this->storeImpl(mFile,*mpObj); }
    bool onDisk() { return mFile.exists(true); }
    void remove() { mFile.remove(); }

    void newFile( String const& fileName ) { mFile = File(fileName); }
    String filename() { return mFile.toString(); }

private:
    File mFile;
    T* mpObj;
};

#endif /* FEUDAL_OBJECTMANAGER_H_ */
