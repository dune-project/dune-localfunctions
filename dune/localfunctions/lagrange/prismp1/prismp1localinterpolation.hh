// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PRISM_P1_LOCALINTERPOLATION_HH
#define DUNE_PRISM_P1_LOCALINTERPOLATION_HH

#include <vector>
#include <dune/localfunctions/common/localinterpolation.hh>

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
    void interpolate (const F& ff, std::vector<C>& out) const
    {
      auto&& f = Impl::makeFunctionWithCallOperator<typename LB::Traits::DomainType>(ff);

      out.resize(6);
      out[0] = f(x[0]);
      out[1] = f(x[1]);
      out[2] = f(x[2]);
      out[3] = f(x[3]);
      out[4] = f(x[4]);
      out[5] = f(x[5]);
    }

  private:
    typename LB::Traits::DomainType x[6];
  };
}

#endif
