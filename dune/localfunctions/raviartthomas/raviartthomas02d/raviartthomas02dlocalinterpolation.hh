// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_RT02DLOCALINTERPOLATION_HH
#define DUNE_RT02DLOCALINTERPOLATION_HH

#include <cmath>
#include <vector>
#include <dune/common/exceptions.hh>

namespace Dune
{
  template<class LB>
  class RT02DLocalInterpolation
  {
  public:

    //! \brief Standard constructor
    RT02DLocalInterpolation ()
    {
      sign0 = sign1 = sign2 = 1.0;
    }

    //! \brief Make set numer s, where 0<=s<8
    RT02DLocalInterpolation (unsigned int s)
    {
      sign0 = sign1 = sign2 = 1.0;
      if (s&1) sign0 *= -1.0;
      if (s&2) sign1 *= -1.0;
      if (s&4) sign2 *= -1.0;
      m0[0] = 0.5; m0[1] = 0.0;
      m1[0] = 0.0; m1[1] = 0.5;
      m2[0] = 0.5; m2[1] = 0.5;
      n0[0] = 0.0;           n0[1] = -1.0;
      n1[0] = -1.0;          n1[1] = 0.0;
      n2[0] = 1.0/sqrt(2.0); n2[1] = 1.0/sqrt(2.0);
      c0 = ( 0.5*n0[0] - 1.0*n0[1]);
      c1 = (-1.0*n1[0] + 0.5*n1[1]);
      c2 = ( 0.5*n2[0] + 0.5*n2[1]);
    }

    template<typename F, typename C>
    void interpolate (const F& f, std::vector<C>& out) const
    {
      // f gives v*outer normal at a point on the edge!
      typename F::Traits::RangeType y;

      out.resize(3);

      f.evaluate(m0,y); out[0] = (y[0]*n0[0]+y[1]*n0[1])*sign0/c0;
      f.evaluate(m1,y); out[1] = (y[0]*n1[0]+y[1]*n1[1])*sign1/c1;
      f.evaluate(m2,y); out[2] = (y[0]*n2[0]+y[1]*n2[1])*sign2/c2;
    }

  private:
    typename LB::Traits::RangeFieldType sign0,sign1,sign2;
    typename LB::Traits::DomainType m0,m1,m2;
    typename LB::Traits::DomainType n0,n1,n2;
    typename LB::Traits::RangeFieldType c0,c1,c2;
  };
}

#endif
