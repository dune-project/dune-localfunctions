// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FINITEELEMENT_HH
#define DUNE_FINITEELEMENT_HH

#include <vector>

#include <dune/common/geometrytype.hh>

namespace Dune {

  //! Traits for the FiniteElement
  /**
   *  Each finite element should either have its own specialization of this
   *  class or provide a traits class conforming to this interface in some
   *  other way.
   *
   *  \note Defining the traits inside the finite element class is impossible
   *        if the finite element is derived from FiniteElementInterface since
   *        that requires the traits class as an (implicit or explicit)
   *        template parameter.
   *
   *  \tparam Imp The type of the basis implementation.
   */
  template<typename Imp>
  struct FiniteElementTraits
  {
    //! The type of the basis
    /**
     *  This is usually something derived from BasisInterface.
     */
    typedef typename Imp::Traits::BasisType BasisType;

    //! The type of the coefficients
    /**
     *  This is usually something derived from LocalCoefficientsInterface
     *
     *  \note The Coefficients class is the same for the local and the
     *        "global" interface, so in contrast to the basis and the
     *        interpolation it is not wrapped or reimplemented.
     */
    typedef typename Imp::Traits::CoefficientsType CoefficientsType;

    //! The type of the interpolation
    /**
     *  This is usually something derived from InterpolationInterface.
     */
    typedef typename Imp::Traits::InterpolationType InterpolationType;
  };

  //! interface for a finite element
  template<typename Imp, typename T = FiniteElementTraits<Imp> >
  class FiniteElementInterface
  {
  public:
    //! Export the traits class
    typedef T Traits;

    //! get a reference to the basis object
    /**
     *  The reference is guaranteed to be valid for as long as this finite
     *  element object exists.
     */
    const typename Traits::BasisType& basis () const
    {
      return asImp().basis();
    }

    //! get a reference to the coefficients object
    /**
     *  The reference is guaranteed to be valid for as long as this finite
     *  element object exists.
     */
    const typename Traits::CoefficientsType& coefficients () const
    {
      return asImp().coefficients();
    }

    //! get a reference to the interpolation object
    /**
     *  The reference is guaranteed to be valid for as long as this finite
     *  element object exists.
     */
    const typename Traits::InterpolationType& interpolation () const
    {
      return asImp().interpolation();
    }

    //! get the geometry type of this finite element object
    GeometryType type () const
    {
      return asImp().type();
    }

  private:
    Imp& asImp () {return static_cast<Imp &> (*this);}
    const Imp& asImp () const {return static_cast<const Imp &>(*this);}
  };

} // namespace Dune

#endif // DUNE_FINITEELEMENT_HH
