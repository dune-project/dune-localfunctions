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

    typedef typename GridView::Grid::ctype DomainField;
    static const unsigned int dimDomain = GridView::dimension;
    typedef FieldVector< DomainField, dimDomain > DomainVector;

    typedef SF RangeField;
    static const unsigned int dimRange = 1;
    typedef FieldVector< RangeField, dimRange > RangeVector;

    typedef LagrangeBasisProvider< dimDomain, RangeField, CF > BasisCreator;
    typedef Dune::LagrangePointsCreator< RangeField, dimDomain > LagrangePointsCreator;
    typedef LagrangePointsCreator LocalCoefficientsCreator;
    typedef LocalLagrangeInterpolationCreator< LagrangePointsCreator > LocalInterpolationCreator;

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
      for( unsigned int topologyId = 0; topologyId < numTopologies; ++topologyId )
      {
        BasisCreator::release( *(basis_[ topologyId ]) );
        LocalInterpolationCreator::release( *(localInterpolation_[ topologyId ]) );
      }
    }

    const DofMapper &dofMapper () const
    {
      return dofMapper_;
    }

    const Basis &basis ( const typename GridView::template Codim< 0 >::Entity &entity ) const
    {
      const unsigned int topologyId = Dune::GenericGeometry::topologyId( entity.type() );
      return *(basis_[ topologyId ]);
    }

  private:
    static const unsigned int numTopologies = (1 << dimDomain);

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
      typedef typename GenericGeometry::Topology< topologyId, dimDomain >::type Topology;
      basis[ topologyId ] = &BasisCreator::template basis< Topology >( order );
      localInterpolation[ topologyId ] = &LocalInterpolationCreator::template localInterpolation< Topology >( order );
    }
  };

}

#endif // #ifndef DUNE_LAGRANGEBASIS_SPACE_HH
