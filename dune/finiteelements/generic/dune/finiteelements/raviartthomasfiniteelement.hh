// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_RAVIARTTHOMASFINITEELEMENT_HH
#define DUNE_RAVIARTTHOMASFINITEELEMENT_HH

#include <dune/finiteelements/generic/localfiniteelement.hh>
#include <dune/finiteelements/raviartthomas/raviartthomasbasis.hh>

namespace Dune
{

  template< unsigned int dimDomain, class D, class R, class CF=D >
  class RaviartThomasLocalFiniteElement
    : public GenericLocalFiniteElement< RaviartThomasBasisCreator< dimDomain, R, CF >,
          RaviartThomasBasisCreator< dimDomain, R, CF >,
          RaviartThomasBasisCreator< dimDomain, R, CF >,
          dimDomain,D,R,CF>
  {
    typedef GenericLocalFiniteElement< RaviartThomasBasisCreator< dimDomain, R, CF >,
        RaviartThomasBasisCreator< dimDomain, R, CF >,
        RaviartThomasBasisCreator< dimDomain, R, CF >,
        dimDomain,D,R,CF> Base;
    using Base::FECreator;
    using Base::FiniteElement;
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
