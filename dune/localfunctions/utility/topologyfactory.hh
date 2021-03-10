// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOMETRY_TOPOLOGYFACTORY_HH
#define DUNE_GEOMETRY_TOPOLOGYFACTORY_HH

#include <cassert>

#include <array>
#include <map>
#include <memory>
#include <type_traits>
#include <vector>

#include <dune/geometry/type.hh>

namespace Dune
{

  /**
   * @brief Provide a factory over the generic topologies
   *
   * This class can be used to dynamically create objects
   * statically bound by their generic topologies.
   * The method create() returns a pointer to an object depending
   * on the topology id and a key; the dimension corresponding
   * to the topology id is static and is provided by the
   * Traits class. A static method (taking the GeometryType as template
   * argument) is also provided.
   * The Traits class must provide the space dimension,
   * the types for the key (Key),
   * the objects returned (Object), and the underlying factory
   * (Factory). This class must have a template method
   * createObject taking a key and returning a pointer to
   * the newly create Object - for destruction call the release
   * method.
   **/
  template <class Traits>
  struct TopologyFactory
  {
    // extract types from Traits class
    static const unsigned int dimension = Traits::dimension;
    typedef typename Traits::Key Key;
    typedef typename Traits::Object Object;
    typedef typename Traits::Factory Factory;

    //! dynamically create objects
    static Object *create ( const Dune::GeometryType &gt, const Key &key )
    {
      return Impl::IfGeometryType< Maker, dimension >::apply( gt.id(), key );
    }

    //! statically create objects
    template< GeometryType::Id geometryId >
    static Object *create ( const Key &key )
    {
      return Factory::template createObject< geometryId >( key );
    }

    //! release the object returned by the create methods
    static void release( Object *object ) { delete object; }

  private:
    // Internal maker class used in ifTopology helper
    template< GeometryType::Id geometryId >
    struct Maker
    {
      static Object *apply ( const Key &key )
      {
        return create< geometryId >( key );
      };
    };
  };



  /** @brief A wrapper for a TopologyFactory providing
   *         singleton storage. Same usage as TopologyFactory
   *         but with empty release method an internal storage.
   **/
  template <class Factory>
  struct TopologySingletonFactory
  {
    static const unsigned int dimension = Factory::dimension;
    typedef typename Factory::Key Key;
    typedef const typename Factory::Object Object;

    //! @copydoc TopologyFactory::create(const Dune::GeometryType &gt,const Key &key)
    static Object *create ( const Dune::GeometryType &gt, const Key &key )
    {
      assert( gt.id() < numTopologies );
      return instance().getObject( gt, key );
    }

    //! @copydoc TopologyFactory::create(const Key &key)
    template< GeometryType::Id geometryId >
    static auto create ( const Key &key )
      -> std::enable_if_t< static_cast<GeometryType>(geometryId).dim() == dimension, Object * >
    {
      return instance().template getObject< geometryId >( key );
    }

    //! @copydoc TopologyFactory::release
    static void release ( Object *object )
    {}

  private:
    struct ObjectDeleter
    {
      void operator() ( Object *ptr ) const { Factory::release( ptr ); }
    };

    static TopologySingletonFactory &instance ()
    {
      static TopologySingletonFactory instance;
      return instance;
    }

    static const unsigned int numTopologies = (1 << dimension);
    typedef std::array< std::unique_ptr< Object, ObjectDeleter >, numTopologies > Array;
    typedef std::map< Key, Array > Storage;

    TopologySingletonFactory () = default;

    std::unique_ptr< Object, ObjectDeleter > &find ( const unsigned int topologyId, const Key &key )
    {
      return storage_[ key ][ topologyId ];
    }

    Object *getObject ( const Dune::GeometryType &gt, const Key &key )
    {
      auto &object = find( gt.id(), key );
      if( !object )
        object.reset( Factory::create( gt, key ) );
      return object.get();
    }

    template< GeometryType::Id geometryId >
    Object *getObject ( const Key &key )
    {
      static constexpr GeometryType geometry = geometryId;
      auto &object = find( geometry.id(), key );
      if( !object )
        object.reset( Factory::template create< geometry >( key ) );
      return object.get();
    }

    Storage storage_;
  };

}

#endif // #ifndef DUNE_GEOMETRY_TOPOLOGYFACTORY_HH
