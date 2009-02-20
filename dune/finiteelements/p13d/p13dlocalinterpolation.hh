// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_P13DLOCALINTERPOLATION_HH
#define DUNE_P13DLOCALINTERPOLATION_HH

#include "../common/localinterpolation.hh"

namespace Dune
{
  template<class LB>
  class P13DLocalInterpolation
    : public LocalInterpolationInterface<P13DLocalInterpolation<LB> >
  {
  public:
    P13DLocalInterpolation ()
    {
      x[0][0] = 0.0; x[0][1] = 0.0; x[0][2] = 0.0;
      x[0][0] = 1.0; x[0][1] = 0.0; x[0][2] = 0.0;
      x[0][0] = 0.0; x[0][1] = 1.0; x[0][2] = 0.0;
      x[0][0] = 0.0; x[0][1] = 0.0; x[0][2] = 1.0;
    }

    //! \brief Local interpolation of a function
    template<typename F, typename C>
    void interpolate (const F& f, std::vector<C>& out) const
    {
      typename LB::Traits::RangeType y;

      out.resize(4);
      f.evaluate(x[0],y); out[0] = y;
      f.evaluate(x[1],y); out[1] = y;
      f.evaluate(x[2],y); out[2] = y;
      f.evaluate(x[3],y); out[3] = y;
    }
  private:
    typename LB::Traits::DomainType x[4];
  };
}

#endif
