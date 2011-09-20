// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_Q2_LOCALINTERPOLATION_HH
#define DUNE_Q2_LOCALINTERPOLATION_HH

#include <vector>

namespace Dune
{
  template<class LB>
  class Q2LocalInterpolation
  {
  public:

    //! \brief Local interpolation of a function
    template<typename F, typename C>
    void interpolate (const F& f, std::vector<C>& out) const
    {
      typename LB::Traits::DomainType x;
      typename LB::Traits::RangeType y;
      static const int dim = LB::Traits::dimDomain;

      // Compute number of Lagrange points
      size_t size = 1;
      for (int i=0; i<dim; i++)
        size *= 3;

      out.resize(size);

      for (size_t i=0; i<size; i++) {

        // Construct the i-th Lagrange point
        size_t ternary = i;
        for (int j=0; j<dim; j++) {

          int digit = ternary%3;
          ternary /= 3;

          x[j] = digit*0.5;

        }

        // Evaluate the function at this point
        f.evaluate(x,y);
        out[i] = y;

      }

    }
  };
}

#endif
