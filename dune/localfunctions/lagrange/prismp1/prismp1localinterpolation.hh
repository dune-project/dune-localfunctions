// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PRISM_P1_LOCALINTERPOLATION_HH
#define DUNE_PRISM_P1_LOCALINTERPOLATION_HH

#include <vector>

namespace Dune
{
  template<class LB>
  class PrismP1LocalInterpolation
  {
  public:
    PrismP1LocalInterpolation ()
    {
      x[0][0] = 0.0; x[0][1] = 0.0; x[0][2] = 0.0;
      x[1][0] = 1.0; x[1][1] = 0.0; x[1][2] = 0.0;
      x[2][0] = 0.0; x[2][1] = 1.0; x[2][2] = 0.0;
      x[3][0] = 0.0; x[3][1] = 0.0; x[3][2] = 1.0;
      x[4][0] = 1.0; x[4][1] = 0.0; x[4][2] = 1.0;
      x[5][0] = 0.0; x[5][1] = 1.0; x[5][2] = 1.0;
    }

    //! \brief Local interpolation of a function
    template<typename F, typename C>
    void interpolate (const F& f, std::vector<C>& out) const
    {
      typename LB::Traits::RangeType y;

      out.resize(6);
      f.evaluate(x[0],y); out[0] = y;
      f.evaluate(x[1],y); out[1] = y;
      f.evaluate(x[2],y); out[2] = y;
      f.evaluate(x[3],y); out[3] = y;
      f.evaluate(x[4],y); out[4] = y;
      f.evaluate(x[5],y); out[5] = y;
    }

  private:
    typename LB::Traits::DomainType x[6];
  };
}

#endif
