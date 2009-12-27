// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALINTERPOLATION_HH
#define DUNE_LOCALINTERPOLATION_HH

#include <vector>

#include <dune/common/function.hh>

namespace Dune
{

  //! Interface class for interpolating a local basis
#ifdef DUNE_VIRTUAL_SHAPEFUNCTIONS
  template<class DomainType, class RangeType>
  class LocalInterpolationInterface
  {
  public:

    typedef VirtualFunction<DomainType, RangeType> FunctionType;
    typedef typename RangeType::field_type CoefficientType;

    //! determine coefficients interpolating a given function
    /**
     * \param[in]  f   Function instance used to interpolate.
     * \param[out] out Resulting coefficients vector.
     */
    virtual void interpolate (const FunctionType& f, std::vector<CoefficientType>& out) const = 0;

  };
#endif

}

#endif
