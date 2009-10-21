// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LAGRANGEFINITEELEMENT_HH
#define DUNE_LAGRANGEFINITEELEMENT_HH

#include <dune/finiteelements/generic/localfiniteelement.hh>
#include <dune/finiteelements/lagrangebasis/lagrangepoints.hh>
#include <dune/finiteelements/lagrangebasis/lagrangebasis.hh>

namespace Dune
{

  template< unsigned int dimD, class D, class R, class SF=R, class CF=SF >
  class LagrangeLocalFiniteElement
    : public GenericLocalFiniteElement< LagrangeBasisFactory< dimD, LagrangeCoefficientsFactory, SF, CF >,
          LagrangeCoefficientsFactory< SF, dimD >,
          LagrangeInterpolationFactory< LagrangeCoefficientsFactory< SF, dimD > >, dimD, D, R >
  {
    typedef GenericLocalFiniteElement< LagrangeBasisFactory< dimD, LagrangeCoefficientsFactory, SF, CF >,
        LagrangeCoefficientsFactory< SF, dimD >,
        LagrangeInterpolationFactory< LagrangeCoefficientsFactory< SF, dimD > >, dimD, D, R >
    Base;

  public:
    /** \todo Please doc me !
     */
    LagrangeLocalFiniteElement ( unsigned int topologyId, unsigned int order )
      : Base( topologyId, order )
    {}
  };

}

#endif // #ifndef DUNE_LAGRANGEFINITEELEMENT_HH
