// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_Q22DLOCALINTERPOLATION_HH
#define DUNE_Q22DLOCALINTERPOLATION_HH

#include <vector>

namespace Dune
{
  template<class LB>
  class Q22DLocalInterpolation
  {
  public:

    //! \brief Local interpolation of a function
    template<typename F, typename C>
    void interpolate (const F& f, std::vector<C>& out) const
    {
      typename LB::Traits::DomainType x;
      typename LB::Traits::RangeType y;

      out.resize(9);
      x[0] = 0.0; x[1] = 0.0; f.evaluate(x,y); out[0] = y;
      x[0] = 1.0; x[1] = 0.0; f.evaluate(x,y); out[1] = y;
      x[0] = 0.0; x[1] = 1.0; f.evaluate(x,y); out[2] = y;
      x[0] = 1.0; x[1] = 1.0; f.evaluate(x,y); out[3] = y;
      x[0] = 0.0; x[1] = 0.5; f.evaluate(x,y); out[4] = y;
      x[0] = 1.0; x[1] = 0.5; f.evaluate(x,y); out[5] = y;
      x[0] = 0.5; x[1] = 0.0; f.evaluate(x,y); out[6] = y;
      x[0] = 0.5; x[1] = 1.0; f.evaluate(x,y); out[7] = y;
      x[0] = 0.5; x[1] = 0.5; f.evaluate(x,y); out[8] = y;

    }
  };
}

#endif
