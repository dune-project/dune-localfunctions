// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PRISM2_3DLOCALFINITEELEMENT_HH
#define DUNE_PRISM2_3DLOCALFINITEELEMENT_HH

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include "prismp2/prismp2localbasis.hh"
#include "prismp2/prismp2localcoefficients.hh"
#include "prismp2/prismp2localinterpolation.hh"

namespace Dune
{

  /** \todo Please doc me !
   */
  template<class D, class R>
  class PrismP2LocalFiniteElement
  {
  public:
    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<PrismP2LocalBasis<D,R>,
        PrismP2LocalCoefficients,
        PrismP2LocalInterpolation<PrismP2LocalBasis<D,R> > > Traits;

    /** \todo Please doc me !
     */
    PrismP2LocalFiniteElement ()
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

    PrismP2LocalFiniteElement* clone () const
    {
      return new PrismP2LocalFiniteElement(*this);
    }

  private:
    PrismP2LocalBasis<D,R> basis;
    PrismP2LocalCoefficients coefficients;
    PrismP2LocalInterpolation<PrismP2LocalBasis<D,R> > interpolation;
    GeometryType gt;
  };

}

#endif
