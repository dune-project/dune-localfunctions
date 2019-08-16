// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_Pk1DLOCALINTERPOLATION_HH
#define DUNE_Pk1DLOCALINTERPOLATION_HH

#include <vector>
#include <dune/localfunctions/common/localinterpolation.hh>

#if 0
namespace Dune
{
  template<class LB>
  class Pk1DLocalInterpolation
  {
    /** \brief The number of degrees of freedom */
    enum {N = LB::N};

    /** \brief Export the element order */
    enum {k = LB::O};

  private:
    static const int kdiv = (k == 0 ? 1 : k);

  public:

    template<typename F, typename C>
    void interpolate (const F& ff, std::vector<C>& out) const
    {
      typename LB::Traits::DomainType x;
      typedef typename LB::Traits::DomainFieldType D;

      auto&& f = Impl::makeFunctionWithCallOperator<decltype(x)>(ff);

      out.resize(N);

      for (int i=0; i<N; i++)
      {
        x[0] = ((D)i)/((D)kdiv);
        out[i] = f(x);
      }
    }

  };
}
#endif

#endif
