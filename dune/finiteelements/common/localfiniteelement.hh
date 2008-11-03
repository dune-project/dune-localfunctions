// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALFINITEELEMENT_HH
#define DUNE_LOCALFINITEELEMENT_HH

#include <iostream>
#include <vector>

#include <dune/common/geometrytype.hh>

namespace Dune {

  //! traits helper struct
  template<class LB, class LC, class LI>
  struct LocalFiniteElementTraits
  {
    typedef LB LocalBasisType;
    typedef LC LocalCoefficientType;
    typedef LI LocalInterpolationType;
  };

  //! interface for a finite element
  template<class T, class Imp>
  class LocalFiniteElementInterface
  {
  public:

    typedef T Traits;

    const typename T::LocalBasisType& localBasis () const
    {
      return asImp().localBasis();
    }

    const typename T::LocalCoefficientsType& localCoefficients () const
    {
      return asImp().localCoefficients();
    }

    const typename T::LocalInterpolationType& localInterpolation () const
    {
      return asImp().localInterpolation();
    }

    GeometryType type () const
    {
      return asImp().type();
    }

  private:
    Imp& asImp () {return static_cast<Imp &> (*this);}
    const Imp& asImp () const {return static_cast<const Imp &>(*this);}
  };

}

#endif
