// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LAGRANGEFINITEELEMENT_HH
#define DUNE_LAGRANGEFINITEELEMENT_HH

#include <dune/finiteelements/generic/localfiniteelement.hh>
#include <dune/finiteelements/lagrangebasis/lagrangecoefficients.hh>
#include <dune/finiteelements/lagrangebasis/interpolation.hh>
#include <dune/finiteelements/lagrangebasis/lagrangebasis.hh>

namespace Dune
{
  template< template <class,unsigned int> class LC,
      unsigned int dimDomain, class D, class R,
      class SF=R, class CF=SF >
  class LagrangeLocalFiniteElement
    : public GenericLocalFiniteElement< LagrangeBasisFactory< LC, dimDomain, SF, CF >,
          LagrangeCoefficientsFactory<LC, dimDomain, SF >,
          LagrangeInterpolationFactory< LC, dimDomain, SF >,
          dimDomain,D,R>
  {
    typedef GenericLocalFiniteElement< LagrangeBasisFactory< LC, dimDomain, SF, CF >,
        LagrangeCoefficientsFactory<LC, dimDomain, SF >,
        LagrangeInterpolationFactory< LC, dimDomain, SF >,
        dimDomain,D,R> Base;
  public:
    typedef typename Base::Traits Traits;

    /** \todo Please doc me !
     */
    LagrangeLocalFiniteElement ( unsigned int topologyId, unsigned int order )
      : Base( topologyId, order )
    {}
  };

}

#endif // #ifndef DUNE_LAGRANGEFINITEELEMENT_HH
