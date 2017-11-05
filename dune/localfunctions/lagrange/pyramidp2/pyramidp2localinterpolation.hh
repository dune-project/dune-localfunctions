// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PYRAMID_P2_LOCALINTERPOLATION_HH
#define DUNE_PYRAMID_P2_LOCALINTERPOLATION_HH

#include <vector>
#include <dune/localfunctions/common/localinterpolation.hh>

namespace Dune
{
  template<class LB>
  class PyramidP2LocalInterpolation
  {
  public:

    //! \brief Local interpolation of a function
    template<typename F, typename C>
    void interpolate (const F& ff, std::vector<C>& out) const
    {
      out.resize(14);
      typename LB::Traits::DomainType x;

      auto&& f = Impl::makeFunctionWithCallOperator<decltype(x)>(ff);

      x[0] = 0.0;   x[1] = 0.0;   x[2] = 0.0;
      out[0] = f(x);

      x[0] = 1.0;   x[1] = 0.0;   x[2] = 0.0;
      out[1] = f(x);

      x[0] = 0.0;   x[1] = 1.0;   x[2] = 0.0;
      out[2] = f(x);

      x[0] = 1.0;   x[1] = 1.0;   x[2] = 0.0;
      out[3] = f(x);

      x[0] = 0.0;   x[1] = 0.0;   x[2] = 1.0;
      out[4] = f(x);

      x[0] = 0.0;   x[1] = 0.5;   x[2] = 0.0;
      out[5] = f(x);

      x[0] = 1.0;   x[1] = 0.5;   x[2] = 0.0;
      out[6] = f(x);

      x[0] = 0.5;   x[1] = 0.0;   x[2] = 0.0;
      out[7] = f(x);

      x[0] = 0.5;   x[1] = 1.0;   x[2] = 0.0;
      out[8] = f(x);

      x[0] = 0.0;   x[1] = 0.0;   x[2] = 0.5;
      out[9] = f(x);

      x[0] = 0.5;   x[1] = 0.0;   x[2] = 0.5;
      out[10] = f(x);

      x[0] = 0.0;   x[1] = 0.5;   x[2] = 0.5;
      out[11] = f(x);

      x[0] = 0.5;   x[1] = 0.5;   x[2] = 0.5;
      out[12] = f(x);

      x[0] = 0.5;   x[1] = 0.5;   x[2] = 0.0;
      out[13] = f(x);
    }

  };
}

#endif
