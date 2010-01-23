// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_P1_LOCALINTERPOLATION_HH
#define DUNE_P1_LOCALINTERPOLATION_HH

#include <vector>

namespace Dune
{
  template<int dim, class LB>
  class P1LocalInterpolation
  {
  public:
    //! \brief Local interpolation of a function
    template<typename F, typename C>
    void interpolate (const F& f, std::vector<C>& out) const
    {
      typename LB::Traits::RangeType y;
      typename LB::Traits::DomainType x;

      out.resize(dim+1);

      // vertex 0
      for (int i=0; i<dim; i++)
        x[i] = 0;
      f.evaluate(x,y); out[0] = y;

      // remaining vertices
      for (int i=0; i<dim; i++) {
        for (int j=0; j<dim; j++)
          x[j] = (i==j);

        f.evaluate(x,y); out[i+1] = y;

      }

    }

  };
}

#endif
