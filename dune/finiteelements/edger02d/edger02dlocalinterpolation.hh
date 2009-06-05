// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_EDGER02DLOCALINTERPOLATION_HH
#define DUNE_EDGER02DLOCALINTERPOLATION_HH

#include "../common/localinterpolation.hh"

namespace Dune
{

  /** \todo Please doc me! */
  template<class LB>
  class EdgeR02DLocalInterpolation
    : public LocalInterpolationInterface<EdgeR02DLocalInterpolation<LB> >
  {
  public:

    //! \brief Local interpolation of a function
    template<typename F, typename C>
    void interpolate (const F& f, std::vector<C>& out) const
    {
      typename LB::Traits::DomainType x;
      typename LB::Traits::RangeType y;

      out.resize(4);

      // Evaluate f at the center of the edge corresponding to the given
      // coefficient.  This is the Gauß quadrature rule:
      //
      //   out[i] = \int_{\Gamma_i} f(x)*b_i(x) dx
      //          = f(x^c_i) * b_i(x^c_i)
      //          = f(x^c_i) * t_i
      //
      // where x^c_i is the position of center of the edge number i and t_i is
      // the tangential unit vector at edge i pointing in the positive
      // direction along that axis.

      x[0] = 0.5; x[1] = 0.0; f.evaluate(x,y); out[0] = y[0];
      x[0] = 0.5; x[1] = 1.0; f.evaluate(x,y); out[1] = y[0];
      x[0] = 0.0; x[1] = 0.5; f.evaluate(x,y); out[2] = y[1];
      x[0] = 1.0; x[1] = 0.5; f.evaluate(x,y); out[3] = y[1];
    }
  };
}

#endif // DUNE_EDGER02DLOCALINTERPOLATION_HH
