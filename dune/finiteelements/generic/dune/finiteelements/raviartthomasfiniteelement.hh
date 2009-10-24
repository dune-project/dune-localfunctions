// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_RAVIARTTHOMASFINITEELEMENT_HH
#define DUNE_RAVIARTTHOMASFINITEELEMENT_HH

#include <dune/finiteelements/generic/localfiniteelement.hh>
#include <dune/finiteelements/raviartthomasbasis/raviartthomasbasis.hh>

namespace Dune
{

  template< unsigned int dimDomain, class D, class R,
      class SF=R, class CF=SF >
  class RaviartThomasLocalFiniteElement
    : public GenericLocalFiniteElement< RaviartThomasBasisFactory< dimDomain, SF, CF >,
          RaviartThomasCoefficientsFactory< dimDomain >,
          RaviartThomasL2InterpolationFactory< dimDomain, SF >,
          dimDomain,D,R>
  {
    typedef GenericLocalFiniteElement< RaviartThomasBasisFactory< dimDomain, SF, CF >,
        RaviartThomasCoefficientsFactory< dimDomain >,
        RaviartThomasL2InterpolationFactory< dimDomain, SF >,
        dimDomain,D,R> Base;
  public:
    using Base::Traits;

    /** \todo Please doc me !
     */
    RaviartThomasLocalFiniteElement ( unsigned int topologyId,
                                      unsigned int order )
      : Base(topologyId,order)
    {}
  };

}

#endif
