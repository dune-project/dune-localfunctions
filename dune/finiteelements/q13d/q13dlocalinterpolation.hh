// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_Q13DLOCALINTERPOLATION_HH
#define DUNE_Q13DLOCALINTERPOLATION_HH

#include "../common/localinterpolation.hh"

namespace Dune
{

  /** \todo Please doc me! */
  template<class LB>
  class Q13DLocalInterpolation
    : public LocalInterpolationInterface<Q13DLocalInterpolation<LB> >
  {
  public:

    //! \brief Local interpolation of a function
    template<typename F, typename C>
    void interpolate (const F& f, std::vector<C>& out) const
    {
      typename LB::Traits::DomainType x;
      typename LB::Traits::RangeType y;

      out.resize(8);
      x[0] = 0.0; x[1] = 0.0; x[2] = 0.0; f.evaluate(x,y); out[0] = y;
      x[0] = 1.0; x[1] = 0.0; x[2] = 0.0; f.evaluate(x,y); out[1] = y;
      x[0] = 0.0; x[1] = 1.0; x[2] = 0.0; f.evaluate(x,y); out[2] = y;
      x[0] = 1.0; x[1] = 1.0; x[2] = 0.0; f.evaluate(x,y); out[3] = y;
      x[0] = 0.0; x[1] = 0.0; x[2] = 1.0; f.evaluate(x,y); out[4] = y;
      x[0] = 1.0; x[1] = 0.0; x[2] = 1.0; f.evaluate(x,y); out[5] = y;
      x[0] = 0.0; x[1] = 1.0; x[2] = 1.0; f.evaluate(x,y); out[6] = y;
      x[0] = 1.0; x[1] = 1.0; x[2] = 1.0; f.evaluate(x,y); out[7] = y;
    }
  };
}

#endif
