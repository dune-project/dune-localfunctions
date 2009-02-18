// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_RT02DLOCALINTERPOLATION_HH
#define DUNE_RT02DLOCALINTERPOLATION_HH

#include <dune/common/exceptions.hh>

#include "../common/localinterpolation.hh"

namespace Dune
{
  template<class LB>
  class RT02DLocalInterpolation
    : public LocalInterpolationInterface<RT02DLocalInterpolation<LB> >
  {
  public:

    //! \brief Standard constructor
    RT02DLocalInterpolation ()
    {}

    //! \brief Make set numer s, where 0<=s<8
    RT02DLocalInterpolation (unsigned int s)
    {
      sign0 = sign1 = sign2 = 1.0;
      if (s&1) sign0 *= -1.0;
      if (s&2) sign1 *= -1.0;
      if (s&4) sign2 *= -1.0;
      m[0][0] = 0.5; m[0][1] = 0.5;
      m[1][0] = 0.0; m[1][1] = 0.5;
      m[2][0] = 0.5; m[2][1] = 0.0;
    }

    template<typename F, typename C>
    void interpolate (const F& f, std::vector<C>& out) const
    {
      // f gives v*outer normal at a point on the edge!
      typename F::Traits::RangeType y;

      out.resize(3);

      f.evaluate(m[0],y); out[0] = y*sign0;
      f.evaluate(m[1],y); out[1] = y*sign1;
      f.evaluate(m[2],y); out[2] = y*sign2;
    }

  private:
    typename LB::Traits::RangeFieldType sign0,sign1,sign2;
    typename LB::Traits::DomainType m[3];
  };
}

#endif
