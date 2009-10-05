// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FINITEELEMENTS_DOFMAPPER_HH
#define DUNE_FINITEELEMENTS_DOFMAPPER_HH

#include <dune/grid/common/indexidset.hh>

#include <dune/grid/genericgeometry/conversion.hh>
#include <dune/grid/genericgeometry/referencetopology.hh>

namespace Dune
{

  // DofMapper
  // ---------

  template< class IndexSet, class Creator >
  class DofMapper
  {
    typedef DofMapper< IndexSet, Creator > This;

    struct SubEntityInfo
    {
      unsigned int codim;
      unsigned int subEntity;
      unsigned int offset;
      unsigned int numDofs;
    };

    struct MapInfo
    {
      unsigned int numDofs;
      std::vector< unsigned int > localDof;
      std::vector< SubEntityInfo > subEntityInfo;
    };

    struct IndexInfo
    {
      unsigned int offset;
      unsigned int size;
      GeometryType type;

      IndexInfo () : size( 0 ) {}
    };

    template< int topologyId >
    struct Build;

  public:
    typedef typename Creator::Key Key;

    static const unsigned int dimension = Creator::dimension;
    static const unsigned int numTopologies = (1 << dimension);

    DofMapper ( const IndexSet &indexSet, const Key &key )
      : indexSet_( indexSet )
    {
      ForLoop< Build, 0, numTopologies-1 >::apply( *this, key );
      update();
    }

    template< class Entity >
    void map ( const Entity &entity, std::vector< unsigned int > indices ) const;

    unsigned int size ( const unsigned int topologyId ) const
    {
      assert( topologyId < numTopologies );
      return mapInfo_[ topologyId ].numDofs;
    }

    template< class Entity >
    unsigned int size ( const Entity &entity ) const
    {
      return size( topologyId( entity.type() ) );
    }

    unsigned int size () const
    {
      return size_;
    }

    void update ();

  private:
    template< class Topology >
    void build ();

    void setupSize ( const unsigned int codim, const unsigned int topologyId, const unsigned int size );

    const IndexSet &indexSet_;
    MapInfo mapInfo_[ numTopologies ];
    std::vector< IndexInfo > indexInfo_[ dimension+1 ];
    unsigned int size_;
  };


  template< class IndexSet, class Creator >
  template< class Entity >
  inline void DofMapper< IndexSet, Creator >
  ::map ( const Entity &entity, std::vector< unsigned int > indices ) const
  {
    typedef typename std::vector< SubEntityInfo >::const_iterator Iterator;
    const MapInfo &mapInfo = mapInfo_[ topologyId( entity.type() ) ];

    indices.resize( mapInfo.numDofs );
    const unsigned int *localDof = &(mapInfo.localDof[ 0 ]);

    const Iterator end = mapInfo.subEntityInfo.end();
    for( Iterator it = mapInfo.subEntityInfo.begin(); it != end; ++it )
    {
      const unsigned int index = indexSet_.index( entity, it->codim, it->subEntity );
      const unsigned int offset = it->offset + index*it->numDofs;
      for( unsigned int j = 0; j < it->numDofs; ++j )
        indices[ *(localDof++) ] = offset + j;
    }
  }


  template< class IndexSet, class Creator >
  inline void DofMapper< IndexSet, Creator >::update ()
  {
    size_ = 0;
    for( unsigned int codim = 0; codim <= dimension; ++codim )
    {
      typedef typename std::vector< IndexInfo >::iterator Iterator;
      const Iterator end = indexInfo_[ codim ].end();
      for( Iterator it = indexInfo_[ codim ].begin(); it != end; ++it )
      {
        it->offset = size_;
        size_ += indexSet_.size( it->type );
      }
    }

    for( unsigned int topologyId = 0; topologyId < numTopologies; ++topologyId )
    {
      typedef typename std::vector< SubEntityInfo >::iterator Iterator;
      MapInfo &mapInfo = mapInfo_[ topologyId ];
      const Iterator end = mapInfo.subEntityInfo.end();
      for( Iterator it = mapInfo.subEntityInfo.begin(); it != end; ++it )
        it->offset = indexInfo_[ it->codim ][ it->subEntity ].offset;
    }
  }


  template< class IndexSet, class Creator >
  template< class Topology >
  inline void DofMapper< IndexSet, Creator >::build ( const Key &key )
  {
    typedef typename Creator::LocalCoefficients LocalCoefficients;

    MapInfo &mapInfo = mapInfo_[ Topology::id ];

    const LocalCoefficients &localCoefficients
      = Creator::typename localCoefficients< Topology >( key );
    mapInfo.numDofs = localCoefficients.size();
    mapInfo.localDof = new unsigned int[ mapinfo.numDofs ];

    const GenericGeometry::ReferenceTopology< dimension > &refTopology
      = GenericGeometry::ReferenceTopologies< dimension >::get( Topology::id );

    SubTopologyMapper< Topology > mapper;
    std::vector< unsigned int > counts( mapper.size(), unsigned int( 0 ) );

    for( unsigned int i = 0; i < mapInfo.numDofs; ++i )
    {
      LocalKey key = localCoefficients.localKey( i );
      ++counts[ mapper( key.codim(), key.subEntity() ) ];
    }

    for( unsigned int codim = 0; codim <= dimension; ++codim )
    {
      const unsigned int codimSize = refTopology::size( codim );

      IndexInfo &indexInfo = indexInfo_[ codim ];
      indexInfo.resize( codimSize );

      for( unsigned int subEntity = 0; subEntity < codimSize; ++subEntity )
      {
        const unsigned int count = counts[ mapper( codim, subEntity ) ];
        if( count == 0 )
          continue;

        SubEntityInfo subEntityInfo;
        subEntityInfo.codim = codim;
        subEntityInfo.subEntity = subEntity;
        subEntityInfo.numDofs = count;
        mapInfo.subEntityInfo.push_back( subEntityInfo );

        setupSize( codim, refTopology.topologyId( codim, subEntity ), count );
      }
    }

    for( unsigned int i = 0; i < mapInfo.numDofs; ++i )
    {
      typedef typename std::vector< SubEntityInfo >::iterator Iterator;
      LocalKey key = localCoefficients.localKey( i );

      unsigned int *localDof = &(mapInfo.localDof[ 0 ]);
      const Iterator end = mapInfo.subEntityInfo.end();
      for( Iterator it = mapInfo.subEntityInfo.begin(); it != end; ++it )
      {
        if( (it->codim == key.codim()) && (it->subEntity = key.subEntity()) )
        {
          *(localDof + key.index()) = i;
          break;
        }
        localDof += it->numDofs;
      }
    }

    Creator::release( localCoefficients );
  }


  template< class IndexSet, class Creator >
  inline void DofMapper< IndexSet, Creator >
  ::setupSize ( const unsigned int codim, const unsigned int topologyId, const unsigned int size )
  {
    if( size == 0 )
      return;
    unsigned int &topologySize = indexInfo_[ codim ][ topologyId ].size;
    if( topologySize > 0 )
    {
      if( topologySize != size )
        DUNE_THROW( InvalidStateException, "Inconsistent LocalCoefficients." );
    }
    else
      topologySize = size;
  }


  template< class IndexSet, class Creator >
  template< int topologyId >
  struct DofMapper< IndexSet, Creator >::Build
  {
    static void apply ( This &dofMapper, const Key &key )
    {
      typedef typename GenericGeometry::Topology< topologyId, dimension >::type Topology;
      dofMapper.template build< Topology >( key );
    }
  };

}

#endif // #ifndef DUNE_FINITEELEMENTS_DOFMAPPER_HH
