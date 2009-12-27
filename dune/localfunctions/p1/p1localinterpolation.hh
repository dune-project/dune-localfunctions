// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_P1_LOCALINTERPOLATION_HH
#define DUNE_P1_LOCALINTERPOLATION_HH

#include "../common/localinterpolation.hh"

namespace Dune
{
  template<int dim, class LB>
  class P1LocalInterpolation
#ifdef DUNE_VIRTUAL_SHAPEFUNCTIONS
    : public LocalInterpolationInterface<typename LB::Traits::DomainType, typename LB::Traits::RangeType>
#endif
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

#if DUNE_VIRTUAL_SHAPEFUNCTIONS
    typedef LocalInterpolationInterface<typename LB::Traits::DomainType, typename LB::Traits::RangeType> Base;

    /** \brief Interpolate a function given as an abstract base class */
    void interpolate(const typename Base::FunctionType& f, typename std::vector<typename Base::CoefficientType>& out) const
    {
      // Call the static implementation with a reference to the base class
      interpolate<typename Base::FunctionType, typename Base::CoefficientType>(f, out);
    }
#endif

  };
}

#endif
