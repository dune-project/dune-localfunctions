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
    typedef FieldVector< typename Geometry::ctype, Geometry::mydimension > DomainVector;

    LocalFunction ( const Function &function, const Entity &entity )
      : function_( function ),
        geometry_( entity.geometry() )
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

  template< class DFS >
  class Interpolation
  {
    typedef Interpolation< DFS > This;

  public:
    typedef DFS DiscreteFunctionSpace;
    typedef typename DiscreteFunctionSpace::RangeField RangeField;
    typedef typename DiscreteFunctionSpace::GridView GridView;
    typedef typename DiscreteFunctionSpace::DofMapper DofMapper;

    typedef typename DiscreteFunctionSpace::LocalInterpolation LocalInterpolation;

    explicit Interpolation ( const DiscreteFunctionSpace &dfSpace )
      : dfSpace_( dfSpace )
    {}

    template< class Function >
    void operator() ( const Function &function, std::vector< RangeField > &dofs ) const
    {
      typedef typename GridView::template Codim< 0 >::Iterator Iterator;
      typedef typename GridView::template Codim< 0 >::Entity Entity;

      const GridView &gridView = dfSpace_.gridView();
      const DofMapper &dofMapper = dfSpace_.dofMapper();
      dofs.resize( dofMapper.size() );

      const Iterator end = gridView.template end< 0 >();
      for( Iterator it = gridView.template begin< 0 >(); it != end; ++it )
      {
        const Entity &entity = *it;
        static std::vector< unsigned int > indices;
        dofMapper.map( entity, indices );

        static std::vector< RangeField > localDofs;
        LocalFunction< Function, Entity > localFunction( function, entity );
        dfSpace_.interpolation( entity ).interpolate( localFunction, localDofs );

        const unsigned int size = indices.size();
        assert( size == localDofs.size() );
        for( unsigned int i = 0; i < size; ++i )
          dofs[ indices[ i ] ] = localDofs[ i ];
      }
    }

  private:
    const DiscreteFunctionSpace &dfSpace_;
  };

}

#endif // #ifndef DUNE_FINITEELEMENTS_INTERPOLATION_HH
