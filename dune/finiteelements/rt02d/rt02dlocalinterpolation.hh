// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_RT02DLOCALINTERPOLATION_HH
#define DUNE_RT02DLOCALINTERPOLATION_HH

#include <dune/common/exceptions.hh>

#include "../common/localinterpolation.hh"

namespace Dune
{
  template<class LB>
  class RT02DLocalInterpolation
    : public LocalInterpolationInterface<RT02DLocalInterpolation<LB> >
  {
  public:

    template<typename F, typename C>
    void interpolate (const F& f, std::vector<C>& out) const
    {
      typename LB::Traits::DomainType x;
      typename LB::Traits::RangeType y;

      out.resize(3);
      DUNE_THROW(Dune::Exception,"interpolate for RT0 not implemented");
    }
  };
}

#endif
