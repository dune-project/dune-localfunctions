// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_P12DLOCALINTERPOLATION_HH
#define DUNE_P12DLOCALINTERPOLATION_HH

#include "../common/localinterpolation.hh"

namespace Dune
{
  template<class LB>
  class P12DLocalInterpolation
    : public LocalInterpolationInterface<P12DLocalInterpolation<LB> >
  {
  public:
    P12DLocalInterpolation ()
    {
      x[0][0] = 0.0; x[0][1] = 0.0;
      x[0][0] = 1.0; x[0][1] = 0.0;
      x[0][0] = 0.0; x[0][1] = 1.0;
    }

    //! \brief Local interpolation of a function
    template<typename F, typename C>
    void interpolate (const F& f, std::vector<C>& out) const
    {
      typename LB::Traits::RangeType y;

      out.resize(3);
      f.evaluate(x[0],y); out[0] = y;
      f.evaluate(x[1],y); out[1] = y;
      f.evaluate(x[2],y); out[2] = y;
    }
  private:
    typename LB::Traits::DomainType x[3];
  };
}

#endif
