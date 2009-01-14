// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_EDGER12DLOCALINTERPOLATION_HH
#define DUNE_EDGER12DLOCALINTERPOLATION_HH

#include "../common/localinterpolation.hh"

namespace Dune
{

  /** \todo Please doc me! */
  template<class LB>
  class EdgeR12DLocalInterpolation
    : public LocalInterpolationInterface<EdgeR12DLocalInterpolation<LB> >
  {
  public:

    //! \brief Local interpolation of a function
    template<typename F, typename C>
    void interpolate (const F& f, std::vector<C>& out) const
    {
      typename LB::Traits::DomainType x;
      typename LB::Traits::RangeType y;

      // Mass matrix M is
      // / 1/3 1/6  0   0  \
      // | 1/6 1/3  0   0  |
      // |  0   0  1/3 1/6 |
      // \  0   0  1/6 1/3 /
      //
      // Inverse Mass M^{-1} matrix is
      // /  4 -2  0  0 \
      // | -2  4  0  0 |
      // |  0  0  4 -2 |
      // \  0  0 -2  4 /

      out.resize(4);

      // Evaluate f at the center of the edge corresponding to the given
      // coefficient.  This is the Gauß quadrature rule:
      //
      //   out[i] = \sum_j M^{-1}_{ij} * \int_{\Gamma_i} f(x)*b_j(x) dx
      //
      // This approach has the advantage that the corresponding
      // basisfunction in the neighboring element will get the same
      // coefficient.  The correct thing AFAICT however would probably be full
      // a full fledged integral:
      //
      //   out[i] = \sum_j M^{-1}_{ij} * \int f(x)*b_j(x) dx
      //
      // where b_j is basisfunction j.

      x[0] = 0.5; x[1] = 0.0; f.evaluate(x,y); out[0] = 4*y[0];
      x[0] = 0.5; x[1] = 1.0; f.evaluate(x,y); out[1] = 4*y[0];
      x[0] = 0.0; x[1] = 0.5; f.evaluate(x,y); out[2] = 4*y[1];
      x[0] = 1.0; x[1] = 0.5; f.evaluate(x,y); out[3] = 4*y[1];

      dune_static_assert(0, "This method has yet to be verified.");
    }
  };
}

#endif // DUNE_EDGER12DLOCALINTERPOLATION_HH
