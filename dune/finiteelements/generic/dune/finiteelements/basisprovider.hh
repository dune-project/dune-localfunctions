// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_BASISPROVIDER_HH
#define DUNE_BASISPROVIDER_HH

#include <map>
#include <dune/grid/genericgeometry/topologytypes.hh>

namespace Dune
{

  // BasisProvider
  // -------------

  template< class BasisCreator >
  struct BasisProvider
  {
    typedef typename BasisCreator::StorageField StorageField;
    static const unsigned int dimension = BasisCreator::dimension;

    typedef typename BasisCreator::Key Key;
    typedef typename BasisCreator::Basis Basis;

    template< class Topology >
    static const Basis &basis ( const Key &key )
    {
      return instance().template getBasis< Topology >( key );
    }

    static const Basis &basis ( const unsigned int topologyId, const Key &key )
    {
      return instance().getBasis( topologyId, key );
    }

    static void release ( const Basis &basis )
    {}

  private:
    static const unsigned int numTopologies = (1 << dimension);

    typedef FieldVector< const Basis *, numTopologies > BasisArray;
    typedef std::map< Key, BasisArray > Storage;

    BasisProvider ()
    {}

    ~BasisProvider ()
    {
      const typename Storage::iterator end = basis_.end();
      for( typename Storage::iterator it = basis_.begin(); it != end; ++it )
      {
        for( unsigned int topologyId = 0; topologyId < numTopologies; ++topologyId )
        {
          const Basis *&basis = it->second[ topologyId ];
          if( basis != 0 )
            BasisCreator::release( *basis );
          basis = 0;
        }
      }
    }

    static BasisProvider &instance ()
    {
      static BasisProvider instance;
      return instance;
    }

    const Basis &getBasis ( const unsigned int topologyId, const Key &key )
    {
      typename Storage::iterator it = basis_.find( key );
      if( it == basis_.end() )
        it = basis_.insert( std::make_pair( key, BasisArray( 0 ) ) ).first;
      const Basis *&basis = it->second[ topologyId ];
      if( basis == 0 )
        GenericGeometry::IfTopology< Maker, dimension >::apply( topologyId, key, basis );
      return *basis;
    }

    template< class Topology >
    const Basis &getBasis ( const Key &key )
    {
      const unsigned int topologyId = Topology::id;
      typename Storage::iterator it = basis_.find( key );
      if( it == basis_.end() )
        it = basis_.insert( std::make_pair( key, BasisArray( 0 ) ) ).first;
      const Basis *&basis = it->second[ topologyId ];
      if( basis == 0 )
        basis = &BasisCreator::template basis< Topology >( key );
      return *basis;
    }

    template< class Topology >
    struct Maker
    {
      static void apply ( const Key &key, const Basis *&basis )
      {
        basis = &BasisCreator::template basis< Topology >( key );
      };
    };

    Storage basis_;
  };

}

#endif // DUNE_BASISPROVIDER_HH
