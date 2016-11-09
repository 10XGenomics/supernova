///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef MAP_H
#define MAP_H


/** Map.h File with helper methods for std::maps.
\file Map.h

Tests for these methods are in testing/MapTest.cc. Please add appropriate 
tests if you add more methods.
*/

#include "feudal/TrackingAllocator.h"
#include "Vec.h"
#include <map>
#include <unordered_map>
#include <ext/hash_map>
using namespace __gnu_cxx;


template <class K, class V, class C=std::less<K>>
using StdMap = std::map<K,V,C,typename DefaultAllocator<std::pair<K const,V>>::type>;

template <class K, class V, class H=std::hash<K>, class P=std::equal_to<K>>
using StdUnorderedMap = std::unordered_map<K,V,H,P,typename DefaultAllocator<std::pair<K const,V>>::type>;

///Return true if k is a key in map m
template<class Key, class Value, typename Cmp> 
bool IsKey( const map<Key, Value, Cmp> & m, const Key & k ){    
  return m.find(k) != m.end( );    
}  

///Return true if k is a key in map m
template<class Key, class Value, typename Hash> 
bool IsKey( const hash_map<Key, Value, Hash> & m, const Key & k ){    
  return m.find(k) != m.end( );    
}  

template<class Key, class Value, typename Cmp> 
std::ostream & operator<<( std::ostream & os, const map<Key, Value, Cmp> & m){
  os << " map size : " << m.size() << endl;
  for (typename std::map<Key,Value,Cmp>::const_iterator iter = m.begin(); 
       iter != m.end(); 
       ++iter) {
    os << iter->first << "\t" << iter->second << "\n";
  }
  return os;
}
  
template<class Key, class Value, typename Hash> 
std::ostream & operator<<( std::ostream & os, 
                           const hash_map<Key, Value, Hash> & m){
  os << " map size : " << m.size() << endl;
  for (typename std::map<Key,Value,Hash>::const_iterator iter = m.begin(); 
       iter != m.end(); 
       ++iter) {
    os << iter->first << "\t" << iter->second << "\n";
  }
  return os;
}
  

///Get all keys in a map (simple expensive version).
template<class Key, class Value, typename Cmp> 
vec<Key> keys(const std::map<Key,Value,Cmp> &  m) {
  vec<Key> ret(m.size());
  int i = 0;
  for (typename std::map<Key,Value,Cmp>::const_iterator iter = m.begin(); 
       iter != m.end(); 
       ++i, ++iter) {
    ret[i] = iter->first;
  }
  return ret;
}

///Get all keys in a hash_map (simple expensive version).
template<class Key, class Value, typename Hash> 
vec<Key> keys(const hash_map<Key,Value,Hash> &  m) {
  vec<Key> ret(m.size());
  int i = 0;
  for (typename hash_map<Key,Value,Hash>::const_iterator iter = m.begin(); 
       iter != m.end(); 
       ++i, ++iter) {
    ret[i] = iter->first;
  }
  return ret;
}

///Get all values in a map (simple expensive version).
template<class Key, class Value, typename Cmp> 
vec<Value> values(const std::map<Key,Value,Cmp> &  m) {
  vec<Value> ret(m.size());
  int i = 0;
  for (typename std::map<Key,Value,Cmp>::const_iterator iter = m.begin();
       iter != m.end();
       ++i, ++iter) {
    ret[i] = iter->second;
  }
  return ret;
}

///Get all keys in a map (cheap but more complex version).
template<class Key, class Value, typename Cmp> 
void keyPointers(const std::map<Key,Value,Cmp> &  m, vec<const Key *> & v) {
  v.resize(m.size());
  int i = 0;
  for (typename std::map<Key,Value,Cmp>::const_iterator iter = m.begin();
       iter != m.end();
       ++i, ++iter) {
    v[i] = &(iter->first);
  }
}

///Get all values in a map (cheap but more complex version).
template<class Key, class Value, typename Cmp> 
void valuePointers(const std::map<Key,Value,Cmp> &  m,
                                 vec<const Value *> & v) {
  v.resize(m.size());
  int i = 0;
  for (typename std::map<Key,Value,Cmp>::const_iterator iter = m.begin();
       iter != m.end();
       ++i, ++iter) {
    v[i] = &(iter->second);
  }
}


/// Adds or updates entry in a map. Use for expensive to construct value objects
template <typename MapType, typename KeyArgType, typename ValueArgType>
typename MapType::iterator efficientAddOrUpdate( MapType& m,
						 const KeyArgType& k,
						 const ValueArgType& v) {
  typename MapType::iterator lb = m.lower_bound(k);

  if (lb != m.end() && !(m.key_comp()(k, lb->first))) {
    lb->second = v;
    return lb;
  } else {
    typedef typename MapType::value_type MVT;
    return m.insert(lb, MVT(k,v));
  }
}

#endif // MAP_H

/*
$Log: not supported by cvs2svn $
Revision 1.6  2007/05/11 13:17:29  palvarez
removed ifdef for gcc 2.*

Revision 1.5  2006/12/22 20:23:30  palvarez
replace an endl with "\n" for possible speed.

Revision 1.4  2004/12/27 18:27:41  palvarez
Added IsKey(), operator<<, and keys() for hash_maps

Revision 1.3  2004/11/05 22:06:29  palvarez
Added operator<<

Revision 1.2  2004/11/04 22:29:13  palvarez
Added method IsKey

Revision 1.1  2004/11/01 21:06:23  palvarez
Added helper methods keys() and values() for a std::map (two versions for
each, a simple expensive one and a cheap more complex one). Put them in new
file Map.h, by analogy to Set.h.

CV: ----------------------------------------------------------------------
  
*/  
    
  
