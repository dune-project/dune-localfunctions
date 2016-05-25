// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_P13DLOCALFINITEELEMENT_HH
#define DUNE_P13DLOCALFINITEELEMENT_HH

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>

#include "prismp1/prismp1localbasis.hh"
#include "prismp1/prismp1localcoefficients.hh"
#include "prismp1/prismp1localinterpolation.hh"

namespace Dune
{

  /** \brief First-order Lagrangian finite element on a prism
   */
  template<class D, class R>
  class PrismP1LocalFiniteElement
  {
  public:
    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<PrismP1LocalBasis<D,R>,PrismP1LocalCoefficients,
        PrismP1LocalInterpolation<PrismP1LocalBasis<D,R> > > Traits;

    /** \todo Please doc me !
     */
    PrismP1LocalFiniteElement ()
    {
      gt.makePrism();
    }


    /** \todo Please doc me !
     */
    const typename Traits::LocalBasisType& localBasis () const
    {
      return basis;
    }

    /** \todo Please doc me !
     */
    const typename Traits::LocalCoefficientsType& localCoefficients () const
    {
      return coefficients;
    }

    /** \todo Please doc me !
     */
    const typename Traits::LocalInterpolationType& localInterpolation () const
    {
      return interpolation;
    }

    /** \brief Number of shape functions in this finite element */
    unsigned int size () const
    {
      return basis.size();
    }

    /** \todo Please doc me !
     */
    GeometryType type () const
    {
      return gt;
    }

    PrismP1LocalFiniteElement* clone () const
    {
      return new PrismP1LocalFiniteElement(*this);
    }

  private:
    PrismP1LocalBasis<D,R> basis;
    PrismP1LocalCoefficients coefficients;
    PrismP1LocalInterpolation<PrismP1LocalBasis<D,R> > interpolation;
    GeometryType gt;
  };

}

#endif
