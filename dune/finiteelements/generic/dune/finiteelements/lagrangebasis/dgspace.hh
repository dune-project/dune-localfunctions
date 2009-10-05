// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LAGRANGEBASIS_DGSPACE_HH
#define DUNE_LAGRANGEBASIS_DGSPACE_HH

#include <dune/finiteelements/dglocalcoefficients.hh>
#include <dune/finiteelements/orthonormalbasis/orthonormalbasis.hh>
#include <dune/finiteelements/orthonormalbasis/l2interpolation.hh>
#include <dune/finiteelements/global/dofmapper.hh>
#include <dune/finiteelements/global/basisproxy.hh>
#include <dune/finiteelements/lagrangebasis.hh>

namespace Dune
{

  template< class GV, class SF, class CF = typename ComputeField< SF, 512 >::Type >
  class LagrangeDGSpace
  {
    typedef LagrangeDGSpace< GV, SF, CF > This;

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
    typedef DGLocalCoefficientsCreator< BasisCreator > LocalCoefficientsCreator;
    typedef LocalL2InterpolationCreator< BasisCreator > LocalInterpolationCreator;

    typedef unsigned int Key;

    typedef typename BasisCreator::Basis LocalBasis;
    typedef typename LocalInterpolationCreator::LocalInterpolation LocalInterpolation;

    typedef BasisProxy< LocalBasis, typename GridView::template Codim< 0 >::Geometry > Basis;

    typedef Dune::DofMapper< typename GridView::IndexSet, LocalCoefficientsCreator > DofMapper;

    explicit LagrangeDGSpace ( const GridView &gridView, const Key &order )
      : gridView_( gridView ),
        dofMapper_( gridView_.indexSet(), order )
    {
      GenericGeometry::ForLoop< Build, 0, numTopologies-1 >::apply( order, basis_, interpolation_ );
    }

    ~LagrangeDGSpace ()
    {
      for( unsigned int topologyId = 0; topologyId < numTopologies; ++topologyId )
      {
        BasisCreator::release( *(basis_[ topologyId ]) );
        LocalInterpolationCreator::release( *(interpolation_[ topologyId ]) );
      }
    }

    const GridView &gridView () const
    {
      return gridView_;
    }

    const DofMapper &dofMapper () const
    {
      return dofMapper_;
    }

    Basis basis ( const typename GridView::template Codim< 0 >::Entity &entity ) const
    {
      const unsigned int topologyId = Dune::GenericGeometry::topologyId( entity.type() );
      return Basis( *(basis_[ topologyId ]), entity.geometry() );
    }

    const LocalInterpolation &interpolation ( const typename GridView::template Codim< 0 >::Entity &entity ) const
    {
      const unsigned int topologyId = Dune::GenericGeometry::topologyId( entity.type() );
      return *(interpolation_[ topologyId ]);
    }

  private:
    static const unsigned int numTopologies = (1 << dimDomain);

    GridView gridView_;
    DofMapper dofMapper_;
    const LocalBasis *basis_[ numTopologies ];
    const LocalInterpolation *interpolation_[ numTopologies ];
  };



  template< class GV, class SF, class CF >
  template< int topologyId >
  struct LagrangeDGSpace< GV, SF, CF >::Build
  {
    static void apply ( const Key &order, const LocalBasis *(&basis)[ numTopologies ],
                        const LocalInterpolation *(&interpolation)[ numTopologies ] )
    {
      typedef typename GenericGeometry::Topology< topologyId, dimDomain >::type Topology;
      basis[ topologyId ] = &BasisCreator::template basis< Topology >( order );
      interpolation[ topologyId ] = &LocalInterpolationCreator::template localInterpolation< Topology >( order );
    }
  };

}

#endif // #ifndef DUNE_LAGRANGEBASIS_DGSPACE_HH
