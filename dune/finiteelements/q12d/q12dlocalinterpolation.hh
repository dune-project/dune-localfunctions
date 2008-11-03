// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_Q12DLOCALINTERPOLATION_HH
#define DUNE_Q12DLOCALINTERPOLATION_HH

#include "../common/localinterpolation.hh"

namespace Dune
{
  template<class LB>
  class Q12DLocalInterpolation
    : public LocalInterpolationInterface<Q12DLocalInterpolation<LB> >
  {
  public:

    //! \brief Local interpolation of a function
    template<typename F, typename C>
    void interpolate (const F& f, std::vector<C>& out) const
    {
      typename LB::Traits::DomainType x;
      typename LB::Traits::RangeType y;

      out.resize(4);
      x[0] = 0.0; x[1] = 0.0; f.evaluate(x,y); out[0] = y;
      x[0] = 1.0; x[1] = 0.0; f.evaluate(x,y); out[1] = y;
      x[0] = 0.0; x[1] = 1.0; f.evaluate(x,y); out[2] = y;
      x[0] = 1.0; x[1] = 1.0; f.evaluate(x,y); out[3] = y;
    }
  };
}

#endif
