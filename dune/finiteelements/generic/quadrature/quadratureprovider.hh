// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_QUADRATUREPROVIDER_HH
#define DUNE_QUADRATUREPROVIDER_HH

#include <map>
#include <dune/grid/genericgeometry/topologytypes.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    // QuadratureProvider
    // ------------------

    template< class QuadratureCreator >
    struct QuadratureProvider
    {
      typedef typename QuadratureCreator::Field Field;
      static const unsigned int dimension = QuadratureCreator::dimension;

      typedef typename QuadratureCreator::Key Key;
      typedef typename QuadratureCreator::Quadrature Quadrature;

      template< class Topology >
      static const Quadrature &quadrature ( const Key &key )
      {
        return instance().template getQuadrature< Topology >( key );
      }

      static const Quadrature &
      quadrature ( const unsigned int topologyId, const Key &key )
      {
        return instance().getQuadrature( topologyId, key );
      }

      static void release ( const Quadrature &quadrature )
      {}

    private:
      static const unsigned int numTopologies = (1 << dimension);

      typedef std::map< Key, const Quadrature * > Storage;

      QuadratureProvider ()
      {}

      ~QuadratureProvider ()
      {
        for( unsigned int i = 0; i < numTopologies; ++i )
        {
          Storage &storage = storage_[ i ];
          const typename Storage::iterator end = storage.end();
          for( typename Storage::iterator it = storage.begin(); it != end; ++it )
          {
            const Quadrature *&quadrature = it->second;
            if( quadrature != 0 )
              QuadratureCreator::release( *quadrature );
            quadrature = 0;
          }
        }
      }

      static QuadratureProvider &instance ()
      {
        static QuadratureProvider instance;
        return instance;
      }

      const Quadrature &getQuadrature ( const unsigned int topologyId, const Key &key )
      {
        Storage &storage = storage_[ topologyId ];
        typename Storage::iterator it = storage.find( key );
        if( it == storage.end() )
        {
          const Quadrature *quadrature;
          GenericGeometry::IfTopology< Maker, dimension >::apply( topologyId, key, quadrature );
          it = storage.insert( std::make_pair( key, quadrature ) ).first;
        }
        return *(it->second);
      }

      template< class Topology >
      const Quadrature &getQuadrature ( const Key &key )
      {
        Storage &storage = storage_[ Topology::id ];
        typename Storage::iterator it = storage.find( key );
        if( it == storage.end() )
        {
          const Quadrature *quadrature = &(QuadratureCreator::template quadrature< Topology >( key ));
          it = storage.insert( std::make_pair( key, quadrature ) ).first;
        }
        return *(it->second);
      }

      template< class Topology >
      struct Maker
      {
        static void apply ( const Key &key, const Quadrature *&quadrature )
        {
          quadrature = &QuadratureCreator::template quadrature< Topology >( key );
        };
      };

      Storage storage_[ numTopologies ];
    };

  }

}

#endif // DUNE_BASISPROVIDER_HH
