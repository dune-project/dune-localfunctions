// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_TOPOLOGYFACTORY_HH
#define DUNE_TOPOLOGYFACTORY_HH

#include <vector>
#include <map>

#include <dune/common/fvector.hh>
#include <dune/grid/genericgeometry/topologytypes.hh>

namespace Dune
{
  template <class Traits>
  struct TopologyFactory
  {
    static const unsigned int dimension = Traits::dimension;
    typedef typename Traits::Key Key;
    typedef typename Traits::Object Object;
    typedef typename Traits::Factory Factory;
    static Object *create(unsigned int topologyId, const Key &key)
    {
      Object *object;
      GenericGeometry::IfTopology< Maker, dimension >::apply( topologyId, key, object );
      return object;
    }
    template <class Topology>
    static Object *create(const Key &key)
    {
      return Factory::template createObject<Topology> ( key );
    }
    static void release( Object *object)
    {
      delete object;
    }
  private:
    template< class Topology >
    struct Maker
    {
      static void apply ( const Key &key, Object *&object )
      {
        object = create<Topology>( key );
      };
    };
  };


  // Singleton wrapper
  template <class Factory>
  struct TopologySingletonFactory
  {
    static const unsigned int dimension = Factory::dimension;
    typedef typename Factory::Key Key;
    typedef const typename Factory::Object Object;
    template< class Topology >
    static Object *create ( const Key &key )
    {
      static_assert( (Topology::dimension == dimension),
                     "Topology with incompatible dimension used" );
      return instance().template getObject< Topology >( key );
    }
    static Object *create ( const unsigned int topologyId, const Key &key )
    {
      assert( topologyId < numTopologies );
      return instance().getObject( topologyId, key );
    }
    static void release ( Object *object )
    {}
  private:
    static TopologySingletonFactory &instance ()
    {
      static TopologySingletonFactory instance;
      return instance;
    }

    static const unsigned int numTopologies = (1 << dimension);
    typedef FieldVector< Object *, numTopologies > Array;
    typedef std::map< Key, Array > Storage;

    TopologySingletonFactory ()
    {}
    ~TopologySingletonFactory ()
    {
      const typename Storage::iterator end = storage_.end();
      for( typename Storage::iterator it = storage_.begin(); it != end; ++it )
      {
        for( unsigned int topologyId = 0; topologyId < numTopologies; ++topologyId )
        {
          Object *&object = it->second[ topologyId ];
          if( object != 0 )
            Factory::release( object );
          object = 0;
        }
      }
    }

    Object *&find( const unsigned int topologyId, const Key &key )
    {
      typename Storage::iterator it = storage_.find( key );
      if( it == storage_.end() )
        it = storage_.insert( std::make_pair( key, Array( 0 ) ) ).first;
      return it->second[ topologyId ];
    }

    Object *getObject ( const unsigned int topologyId, const Key &key )
    {
      Object *&object = find(topologyId,key);
      if( object == 0 )
        object = Factory::create( topologyId, key );
      return object;
    }
    template< class Topology >
    Object *getObject ( const Key &key )
    {
      Object *&object = find(Topology::id,key);
      if( object == 0 )
        object = Factory::template create< Topology >( key );
      return object;
    }
    Storage storage_;
  };

}
#endif
