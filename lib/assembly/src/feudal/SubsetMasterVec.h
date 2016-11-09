// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

#ifndef FEUDAL_SUBSET_MASTERVEC_H
#define FEUDAL_SUBSET_MASTERVEC_H

#include "feudal/VirtualMasterVec.h"
#include <unordered_map>

template <class T>
class SubsetMasterVec : public MasterVec<T>
{
public:
     typedef MasterVec<T> base_type;
     typedef typename base_type::size_type size_type;
     typedef typename base_type::reference reference;
     typedef typename base_type::const_reference const_reference;
     typedef std::unordered_map< size_type, size_type> map_type;

     template <class source_type>
     SubsetMasterVec( source_type const& source, vec<size_type>& ids ) { build(source, ids); };

     SubsetMasterVec( String const& filename, vec<size_type>& ids )
          { build(VirtualMasterVec<T>(filename), ids); };

     SubsetMasterVec( char const* filename, vec<size_type>& ids )
          { build(VirtualMasterVec<T>(filename), ids); };

     template <class source_type>
     void build( source_type const& src, vec<size_t>& ids ) {
          // we'll uniquesort your ids - caveat instantiator
          if ( ! is_sorted_strict(ids.begin(), ids.end() ) )
               UniqueSort(ids);
          for ( auto const id : ids ) {
               AssertLt( id, src.size() );
               mMap[id] = this->size();
               this->push_back( src.at( id ) );
          }
     }

     inline reference operator[]( size_type idx ) {
          try { return base_type::operator[]( mMap.at(idx) ); }
          catch (std::out_of_range const& ) {
               FatalErr( "attempt to access non-existent subset value" ); }
     }

     inline const_reference operator[]( size_type idx ) const {
          try { return base_type::operator[]( mMap.at(idx) ); }
          catch (std::out_of_range const& ) {
               FatalErr( "attempt to access non-existent subset value" ); }
     }

     inline reference at( size_type idx ) {
          return base_type::at( mMap.at(idx) );
     }

     inline const_reference at(size_type idx ) const {
          return base_type::at( mMap.at(idx) );
     }


private:
     map_type mMap;
};

#endif // FEUDAL_SUBSET_MASTERVEC_H
