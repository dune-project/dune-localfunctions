// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FINITEELEMENTS_DOFMAPPER_HH
#define DUNE_FINITEELEMENTS_DOFMAPPER_HH

#include <dune/grid/common/indexidset.hh>

#include <dune/grid/genericgeometry/conversion.hh>
#include <dune/grid/genericgeometry/referencetopologies.hh>

#include <dune/finiteelements/common/localcoefficients.hh>

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
      unsigned int topologyId;
      unsigned int numDofs;
      unsigned int offset;
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
      unsigned int codim;
      unsigned int topologyId;

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
      GenericGeometry::ForLoop< Build, 0, numTopologies-1 >::apply( *this, key );
      update();
    }

    template< class Entity >
    void map ( const Entity &entity, std::vector< unsigned int > &indices ) const;

    unsigned int size ( const unsigned int topologyId ) const
    {
      assert( topologyId < numTopologies );
      return mapInfo_[ topologyId ].numDofs;
    }

    template< class Entity >
    unsigned int size ( const Entity &entity ) const
    {
      return size( GenericGeometry::topologyId( entity.type() ) );
    }

    unsigned int size () const
    {
      return size_;
    }

    void update ();

  private:
    template< class Topology >
    void build ( const Key &key );

    const IndexSet &indexSet_;
    MapInfo mapInfo_[ numTopologies ];
    std::vector< IndexInfo > indexInfo_[ dimension+1 ];
    unsigned int size_;
  };


  template< class IndexSet, class Creator >
  template< class Entity >
  inline void DofMapper< IndexSet, Creator >
  ::map ( const Entity &entity, std::vector< unsigned int > &indices ) const
  {
    typedef typename std::vector< SubEntityInfo >::const_iterator Iterator;
    const MapInfo &mapInfo = mapInfo_[ GenericGeometry::topologyId( entity.type() ) ];

    indices.resize( mapInfo.numDofs );
    const unsigned int *localDof = &(mapInfo.localDof[ 0 ]);

    const Iterator end = mapInfo.subEntityInfo.end();
    for( Iterator it = mapInfo.subEntityInfo.begin(); it != end; ++it )
    {
      const unsigned int index = indexSet_.subIndex( entity, it->subEntity, it->codim );
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
        const GeometryType type = GenericGeometry::geometryType( it->topologyId, dimension-(it->codim) );
        size_ += (it->size > 0 ? indexSet_.size( type ) * it->size : 0);
      }
    }

    for( unsigned int topologyId = 0; topologyId < numTopologies; ++topologyId )
    {
      typedef typename std::vector< SubEntityInfo >::iterator Iterator;
      MapInfo &mapInfo = mapInfo_[ topologyId ];
      const Iterator end = mapInfo.subEntityInfo.end();
      for( Iterator it = mapInfo.subEntityInfo.begin(); it != end; ++it )
        it->offset = indexInfo_[ it->codim ][ it->topologyId >> 1 ].offset;
    }
  }


  template< class IndexSet, class Creator >
  template< class Topology >
  inline void DofMapper< IndexSet, Creator >::build ( const Key &key )
  {
    typedef typename Creator::LocalCoefficients LocalCoefficients;

    MapInfo &mapInfo = mapInfo_[ Topology::id ];

    const LocalCoefficients &localCoefficients
      = Creator::template localCoefficients< Topology >( key );
    mapInfo.numDofs = localCoefficients.size();
    mapInfo.localDof.resize( mapInfo.numDofs );

    const GenericGeometry::ReferenceTopology< dimension > &refTopology
      = GenericGeometry::ReferenceTopologies< dimension >::get( Topology::id );

    GenericGeometry::SubTopologyMapper< Topology > mapper;
    std::vector< unsigned int > counts( mapper.size(), (unsigned int)0 );

    for( unsigned int i = 0; i < mapInfo.numDofs; ++i )
    {
      const LocalKey &key = localCoefficients.localKey( i );
      ++counts[ mapper( key.codim(), key.subEntity() ) ];
    }

    for( unsigned int codim = 0; codim <= dimension; ++codim )
    {
      const unsigned int subdimension = dimension-codim;
      indexInfo_[ codim ].resize( subdimension > 0 ? 1 << (subdimension-1) : 1 );

      const unsigned int codimSize = refTopology.size( codim );
      for( unsigned int subEntity = 0; subEntity < codimSize; ++subEntity )
      {
        const unsigned int topologyId = refTopology.topologyId( codim, subEntity );
        IndexInfo &indexInfo = indexInfo_[ codim ][ topologyId >> 1 ];
        indexInfo.codim = codim;
        indexInfo.topologyId = topologyId;

        const unsigned int count = counts[ mapper( codim, subEntity ) ];
        if( count == 0 )
          continue;

        SubEntityInfo subEntityInfo;
        subEntityInfo.codim = codim;
        subEntityInfo.subEntity = subEntity;
        subEntityInfo.topologyId = topologyId;
        subEntityInfo.numDofs = count;
        mapInfo.subEntityInfo.push_back( subEntityInfo );

        if( indexInfo.size > 0 )
        {
          if( indexInfo.size != count )
            DUNE_THROW( InvalidStateException, "Inconsistent LocalCoefficients." );
        }
        else
          indexInfo.size = count;
      }
    }

    for( unsigned int i = 0; i < mapInfo.numDofs; ++i )
    {
      typedef typename std::vector< SubEntityInfo >::iterator Iterator;
      const LocalKey &key = localCoefficients.localKey( i );

      unsigned int *localDof = &(mapInfo.localDof[ 0 ]);
      const Iterator end = mapInfo.subEntityInfo.end();
      for( Iterator it = mapInfo.subEntityInfo.begin(); true; ++it )
      {
        if( it == end )
        {
          std::cerr << "Error: (subEntity = " << key.subEntity()
                    << ", codim = " << key.codim()
                    << ") not found in subEntityInfo" << std::endl;
          std::cerr << "SubEntityInfo contains:" << std::endl;
          for( it = mapInfo.subEntityInfo.begin(); it != end; ++it )
          {
            std::cerr << "  (subEntity = " << it->subEntity
                      << ", codim = " << it->codim << ")" << std::endl;
          }
          abort();
        }

        if( (it->codim == key.codim()) && (it->subEntity == key.subEntity()) )
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
