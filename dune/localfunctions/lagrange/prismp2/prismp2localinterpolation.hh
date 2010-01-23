// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PRISM_P2_LOCALINTERPOLATION_HH
#define DUNE_PRISM_P2_LOCALINTERPOLATION_HH

#include <vector>

namespace Dune
{
  template<class LB>
  class PrismP2LocalInterpolation
  {
  public:

    //! \brief Local interpolation of a function
    template<typename F, typename C>
    void interpolate (const F& f, std::vector<C>& out) const
    {
      typename LB::Traits::RangeType y;

      out.resize(18);
      typename LB::Traits::DomainType x;

      x[0] = 0.0;   x[1] = 0.0;   x[2] = 0.0;
      f.evaluate(x,y); out[0] = y;

      x[0] = 1.0;   x[1] = 0.0;   x[2] = 0.0;
      f.evaluate(x,y); out[1] = y;

      x[0] = 0.0;   x[1] = 1.0;   x[2] = 0.0;
      f.evaluate(x,y); out[2] = y;

      x[0] = 0.0;   x[1] = 0.0;   x[2] = 1.0;
      f.evaluate(x,y); out[3] = y;

      x[0] = 1.0;   x[1] = 0.0;   x[2] = 1.0;
      f.evaluate(x,y); out[4] = y;

      x[0] = 0.0;   x[1] = 1.0;   x[2] = 1.0;
      f.evaluate(x,y); out[5] = y;

      x[0] = 0.0;   x[1] = 0.0;   x[2] = 0.5;
      f.evaluate(x,y); out[6] = y;

      x[0] = 1.0;   x[1] = 0.0;   x[2] = 0.5;
      f.evaluate(x,y); out[7] = y;

      x[0] = 0;   x[1] = 1.0;   x[2] = 0.5;
      f.evaluate(x,y); out[8] = y;

      x[0] = 0.5;   x[1] = 0.0;   x[2] = 0.0;
      f.evaluate(x,y); out[9] = y;

      x[0] = 0.0;   x[1] = 0.5;   x[2] = 0.0;
      f.evaluate(x,y); out[10] = y;

      x[0] = 0.5;   x[1] = 0.5;   x[2] = 0.0;
      f.evaluate(x,y); out[11] = y;

      x[0] = 0.5;   x[1] = 0.0;   x[2] = 1.0;
      f.evaluate(x,y); out[12] = y;

      x[0] = 0.0;   x[1] = 0.5;   x[2] = 1.0;
      f.evaluate(x,y); out[13] = y;

      x[0] = 0.5;   x[1] = 0.5;   x[2] = 1.0;
      f.evaluate(x,y); out[14] = y;

      x[0] = 0.5;   x[1] = 0.0;   x[2] = 0.5;
      f.evaluate(x,y); out[15] = y;

      x[0] = 0.0;   x[1] = 0.5;   x[2] = 0.5;
      f.evaluate(x,y); out[16] = y;

      x[0] = 0.5;   x[1] = 0.5;   x[2] = 0.5;
      f.evaluate(x,y); out[17] = y;
    }

  };
}

#endif
