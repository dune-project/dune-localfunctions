// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOBATTOFINITEELEMENT_HH
#define DUNE_LOBATTOFINITEELEMENT_HH

#if HAVE_ALGLIB

#include <dune/finiteelements/generic/localfiniteelement.hh>
#include <dune/finiteelements/lagrangebasis/lobattopoints.hh>
#include <dune/finiteelements/lagrangebasis/lagrangebasis.hh>

namespace Dune
{

  template< unsigned int dimDomain, class D, class R,
      class SF=R, class CF=R >
  class LobattoLocalFiniteElement
    : public GenericLocalFiniteElement< LagrangeBasisProvider< dimDomain, SF, CF >,
          LobattoPointsCreator< SF, dimDomain >,
          LocalLagrangeInterpolationCreator< LobattoPointsCreator< SF, dimDomain > >,
          dimDomain,D,R>
  {
    typedef GenericLocalFiniteElement< LagrangeBasisProvider< dimDomain, SF, CF >,
        LobattoPointsCreator< SF, dimDomain >,
        LocalLagrangeInterpolationCreator< LobattoPointsCreator< SF, dimDomain > >,
        dimDomain,D,R> Base;
    using Base::FECreator;
    using Base::FiniteElement;
  public:
    using Base::Traits;

    /** \todo Please doc me !
     */
    LobattoLocalFiniteElement ( unsigned int topologyId,
                                unsigned int order )
      : Base(topologyId,order)
    {}
  };

}

#else
#warning LobattoLocalFiniteElement only available with ALGLib
#endif

#endif
