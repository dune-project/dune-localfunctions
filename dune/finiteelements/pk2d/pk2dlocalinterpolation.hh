// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_Pk2DLOCALINTERPOLATION_HH
#define DUNE_Pk2DLOCALINTERPOLATION_HH

#include "../common/localinterpolation.hh"

namespace Dune
{
  template<class LB>
  class Pk2DLocalInterpolation
    : public LocalInterpolationInterface<Pk2DLocalInterpolation<LB> >
  {
    enum {N = LB::N};
    enum {k = LB::O};
  public:

    template<typename F, typename C>
    void interpolate (const F& f, std::vector<C>& out) const
    {
      typename LB::Traits::DomainType x;
      typename LB::Traits::RangeType y;
      typedef typename LB::Traits::DomainType D;
      out.resize(N);
      int n=0;
      for (int j=0; j<=k; j++)
        for (int i=0; i<=k-j; i++)
        {
          x[0] = ((D)i)/((D)k); x[1] = ((D)j)/((D)k);
          f.evaluate(x,y);
          out[n] = y;
          n++;
        }
    }

  };
}

#endif
