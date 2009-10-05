// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LAGRANGEBASIS_SPACE_HH
#define DUNE_LAGRANGEBASIS_SPACE_HH

#include <dune/finiteelements/lagrangebasis/lagrangebasis.hh>
#include <dune/finiteelements/global/dofmapper.hh>

namespace Dune
{

  template< class GV, class SF, class CF = typename ComputeField< SF, 512 >::Type >
  class LagrangeSpace
  {
    typedef LagrangeSpace< GV, SF, CF > This;

    template< int topologyId >
    struct Build;

  public:
    typedef GV GridView;
    static const unsigned int dimension = GridView::dimension;

    typedef LagrangeBasisProvider< dimension, SF, CF > BasisCreator;
    typedef LagrangePointsCreator< SF, dimension > LocalCoefficientsCreator;
    typedef LocalLagrangeInterpolationCreator< SF, dimension > LocalInterpolationCreator;

    typedef unsigned int Key;

    typedef typename BasisCreator::Basis Basis;
    typedef typename LocalInterpolationCreator::LocalInterpolation LocalInterpolation;

    typedef Dune::DofMapper< typename GridView::IndexSet, LocalCoefficientsCreator > DofMapper;

    explicit LagrangeSpace ( const GridView &gridView, const Key &order )
      : gridView_( gridView ),
        dofMapper_( gridView_.indexSet(), order )
    {
      GenericGeometry::ForLoop< Build, 0, numTopologies-1 >::apply( order, basis_, localInterpolation_ );
    }

    ~LagrangeSpace ()
    {
      for( unsigned int i = 0; i < numTopologies; ++i )
      {
        BasisProvider::release( basis_[ i ] );
        LocalInterpolationCreator::release( localInterpolation_[ i ] );
      }
    }

    const DofMapper &dofMapper () const
    {
      return dofMapper_;
    }

    const Basis &basis ( const typename GridView::template Codim< 0 >::Entity &entity )
    {
      const unsigned int topologyId = Dune::GenericGeometry::topologyId( entity.type() );
      return basis_[ topologyId ];
    }

  private:
    static const unsigned int numTopologies = (1 << dimension);

    GridView gridView_;
    DofMapper dofMapper_;
    const Basis *basis_[ numTopologies ];
    const LocalInterpolation *localInterpolation_[ numTopologies ];
  };



  template< class GV, class SF, class CF >
  template< int topologyId >
  struct LagrangeSpace< GV, SF, CF >::Build
  {
    static void apply ( const Key &order, const Basis *(&basis)[ numTopologies ],
                        const LocalInterpolation *(&localInterpolation)[ numTopologies ] )
    {
      typedef typename GenericGeometry::Topology< topologyId, dimension >::type Topology;
      basis[ topologyId ] = &BasisCreator::template basis< Topology >( order );
      localInterpolation[ topologyId ] = &LocalInterpolationCreator::template localInterpolation< Topology >( order );
    }
  };

}

#endif // #ifndef DUNE_LAGRANGEBASIS_SPACE_HH
