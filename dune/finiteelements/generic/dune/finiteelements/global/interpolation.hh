// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FINITEELEMENTS_INTERPOLATION_HH
#define DUNE_FINITEELEMENTS_INTERPOLATION_HH

namespace Dune
{

  // LocalFunction
  // -------------

  template< class Function, class Entity >
  struct LocalFunction
  {
    typedef typename Function::RangeVector RangeVector;

    typedef typename Entity::Geometry Geometry;
    typedef FieldVector< typename Geometry::ctype, Geometry::dimension > DomainVector;

    LocalFunction ( const Function &function, const Entity &entity )
      : function_( function ),
        geometry_( entity.geometry_ )
    {}

    RangeVector operator() ( const DomainVector &x ) const
    {
      return function_( geometry_.global( x ) );
    }

  private:
    const Function &function_;
    const Geometry &geometry_;
  };



  // Interpolation
  // -------------

  template< class GridView, class DofMapper, class Creator >
  class Interpolation
  {
    typedef Interpolation< GridView, DofMapper, Creator > This;

    template< int topologyId >
    struct Build;

  public:
    typedef typename Creator::Key Key;
    typedef typename Creator::Field Field;

    typedef typename Creator::LocalInterpolation LocalInterpolation;

    static const unsigned int dimension = Creator::dimension;
    static const unsigned int numTopologies = (1 << dimension);

    Interpolation ( const GridView &gridView, const DofMapper &dofMapper, const Key &key )
      : gridView_( gridView ),
        dofMapper_( dofMapper )
    {
      GenericGeometry::ForLoop< Build, 0, numTopologies-1 >::apply( key, localInterpolation_ );
    }

    ~Interpolation ()
    {
      for( unsigned int topologyId = 0; topologyId < numTopologies; ++topologyId )
        Creator::release( localInterpolation_[ topologyId ] );
    }

    template< class Function >
    void operator() ( const Function &function, std::vector< Field > &dofVector ) const
    {
      typedef typename GridView::template Codim< 0 >::Iterator Iterator;
      typedef typename GridView::template Codim< 0 >::Entity Entity;

      assert( dofVector.size() == dofMapper_.size() );
      std::vector< unsigned int > indices;
      std::vector< Field > localDofs;

      const Iterator end = gridView_.end();
      for( Iterator it = gridView_.begin(); it != end; ++it )
      {
        const Entity &entity = *it;
        dofMapper_.map( entity, indices );

        const LocalInterpolation &interpolation = localInterpolation( entity );
        LocalFunction< Function, Entity > localFunction( function, entity );
        interpolation.interpolate( localFunction, localDofs );

        const unsigned int size = indices.size();
        assert( size == localDofs.size() );
        for( unsigned int i = 0; i < size; ++i )
          dofVector[ indices[ i ] ] = localDofs[ i ];
      }
    }

    const LocalInterpolation &
    localInterpolation ( const typename GridView::template Codim< 0 >::Entity &entity ) const
    {
      return localInterpolation_[ topologyId( entity.type() ) ];
    }

  private:
    const GridView &gridView_;
    const DofMapper &dofMapper_;
    const LocalInterpolation *localInterpolation_[ numTopologies ];
  };



  template< class GridView, class DofMapper, class Creator >
  template< int topologyId >
  struct Interpolation< GridView, DofMapper, Creator >::Build
  {
    static void apply ( const Key &key,
                        const LocalInterpolation *(&localInterpolation)[ numTopologies ] )
    {
      typedef typename GenericGeometry::Topology< topologyId, dimension >::type Topology;
      localInterpolation[ topologyId ] = &Creator::template localInterpolation< Topology >( key );
    }
  };
}

#endif // #ifndef DUNE_FINITEELEMENTS_INTERPOLATION_HH
