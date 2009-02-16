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
      nu[0][0] = 1.0/sqrt(2.0); nu[0][1] = 1.0/sqrt(2.0);
      nu[1][0] = -1.0; nu[1][1] = 0.0;
      nu[2][0] = 0.0; nu[2][1] = -1.0;
      if (s&1) nu[0] *= -1.0;
      if (s&2) nu[1] *= -1.0;
      if (s&4) nu[2] *= -1.0;
      m[0][0] = 0.5; m[0][1] = 0.5;
      m[1][0] = 0.0; m[1][1] = 0.5;
      m[2][0] = 0.5; m[2][1] = 0.0;
    }

    template<typename F, typename C>
    void interpolate (const F& f, std::vector<C>& out) const
    {
      typename LB::Traits::DomainType x;
      typename LB::Traits::RangeType y;

      out.resize(3);
      f.evaluate(m[0],y); out[0] = y*nu[0];
      f.evaluate(m[1],y); out[1] = y*nu[1];
      f.evaluate(m[2],y); out[2] = y*nu[2];
    }

  private:
    typename LB::Traits::RangeType nu[3];
    typename LB::Traits::DomainType m[3];
  };
}

#endif
