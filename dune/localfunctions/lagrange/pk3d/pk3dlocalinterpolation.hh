// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_Pk3DLOCALINTERPOLATION_HH
#define DUNE_Pk3DLOCALINTERPOLATION_HH

#include <vector>

namespace Dune
{
  template<class LB>
  class Pk3DLocalInterpolation
  {
    enum {N = LB::N};
    enum {k = LB::O};

  private:
    static const int kdiv = (k == 0 ? 1 : k);
  public:

    template<typename F, typename C>
    void interpolate (const F& f, std::vector<C>& out) const
    {
      typename LB::Traits::DomainType x;
      typename LB::Traits::RangeType y;
      typedef typename LB::Traits::DomainFieldType D;
      out.resize(N);
      int n=0;
      for (int i2 = 0; i2 <= k; i2++)
        for (int i1 = 0; i1 <= k-i2; i1++)
          for (int i0 = 0; i0 <= k-i1-i2; i0++)
          {
            x[0] = ((D)i0)/((D)kdiv);
            x[1] = ((D)i1)/((D)kdiv);
            x[2] = ((D)i2)/((D)kdiv);
            f.evaluate(x,y);
            out[n] = y;
            n++;
          }
    }

  };
}

#endif
