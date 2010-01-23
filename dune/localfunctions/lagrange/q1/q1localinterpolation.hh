// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_Q1_LOCALINTERPOLATION_HH
#define DUNE_Q1_LOCALINTERPOLATION_HH

#include <vector>

namespace Dune
{

  /** \todo Please doc me! */
  template<int dim, class LB>
  class Q1LocalInterpolation
  {
  public:

    //! \brief Local interpolation of a function
    template<typename F, typename C>
    void interpolate (const F& f, std::vector<C>& out) const
    {
      typename LB::Traits::DomainType x;
      typename LB::Traits::RangeType y;

      out.resize(1<<dim);

      for (int i=0; i< (1<<dim); i++) {

        // Generate coordinate of the i-th corner of the reference cube
        // We could use the ReferenceElement for this as well, but it is
        // still not clear how dune-localfunctions should have access to them.
        for (int j=0; j<dim; j++)
          x[j] = (i & (1<<j)) ? 1.0 : 0.0;

        f.evaluate(x,y); out[i] = y;

      }
    }

  };
}

#endif
