// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LAGRANGEFINITEELEMENT_HH
#define DUNE_LAGRANGEFINITEELEMENT_HH

#include <dune/finiteelements/generic/localfiniteelement.hh>
#include <dune/finiteelements/lagrangebasis/lagrangepoints.hh>
#include <dune/finiteelements/lagrangebasis/lagrangebasis.hh>

namespace Dune
{

  template< unsigned int dimDomain, class D, class R, class CF=D >
  class LagrangeLocalFiniteElement
    : public GenericLocalFiniteElement< LagrangeBasisProvider< dimDomain, R, CF >,
          LagrangePointsCreator< R, dimDomain >,
          LocalLagrangeInterpolationCreator< LagrangePointsCreator< R, dimDomain > >,
          dimDomain,D,R,CF>
  {
    typedef GenericLocalFiniteElement< LagrangeBasisProvider< dimDomain, D, CF >,
        LagrangePointsCreator< D, dimDomain >,
        LocalLagrangeInterpolationCreator< LagrangePointsCreator< D, dimDomain > >,
        dimDomain,D,R,CF> Base;
    using Base::FECreator;
    using Base::FiniteElement;
  public:
    using Base::Traits;

    /** \todo Please doc me !
     */
    LagrangeLocalFiniteElement ( unsigned int topologyId,
                                 unsigned int order )
      : Base(topologyId,order)
    {}
  };

}

#endif
