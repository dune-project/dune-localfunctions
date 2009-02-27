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
    /** \todo Please doc me !
     */
    typedef LB LocalBasisType;

    /** \todo Please doc me !
     */
    typedef LC LocalCoefficientsType;

    /** \todo Please doc me !
     */
    typedef LI LocalInterpolationType;
  };

  //! interface for a finite element

#if DUNE_VIRTUAL_SHAPEFUNCTIONS
  template<class T>
#else
  template<class T, class Imp>
#endif
  class LocalFiniteElementInterface
  {
  public:

    /** \todo Please doc me !
     */
    typedef T Traits;

    /** \todo Please doc me !
     */
#if DUNE_VIRTUAL_SHAPEFUNCTIONS
    const virtual typename T::LocalBasisType& localBasis () const = 0;
#else
    const typename T::LocalBasisType& localBasis () const
    {
      return asImp().localBasis();
    }
#endif

    /** \todo Please doc me !
     */
#if DUNE_VIRTUAL_SHAPEFUNCTIONS
    const virtual typename T::LocalCoefficientsType& localCoefficients () const = 0;
#else
    const typename T::LocalCoefficientsType& localCoefficients () const
    {
      return asImp().localCoefficients();
    }
#endif

    /** \todo Please doc me !
     */
#if DUNE_VIRTUAL_SHAPEFUNCTIONS
    const virtual typename T::LocalInterpolationType& localInterpolation () const = 0;
#else
    const typename T::LocalInterpolationType& localInterpolation () const
    {
      return asImp().localInterpolation();
    }
#endif

    /** \todo Please doc me !
     */
#if DUNE_VIRTUAL_SHAPEFUNCTIONS
    virtual GeometryType type () const = 0;
#else
    GeometryType type () const
    {
      return asImp().type();
    }
#endif

#ifndef DUNE_VIRTUAL_SHAPEFUNCTIONS
  private:
    Imp& asImp () {return static_cast<Imp &> (*this);}
    const Imp& asImp () const {return static_cast<const Imp &>(*this);}
#endif
  };

}

#endif
