// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PYRAMIDP2_3DLOCALFINITEELEMENT_HH
#define DUNE_PYRAMIDP2_3DLOCALFINITEELEMENT_HH

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include "pyramidp2/pyramidp2localbasis.hh"
#include "pyramidp2/pyramidp2localcoefficients.hh"
#include "pyramidp2/pyramidp2localinterpolation.hh"

namespace Dune
{

  /** \todo Please doc me !
   */
  template<class D, class R>
  class PyramidP2LocalFiniteElement
  {
  public:
    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<PyramidP2LocalBasis<D,R>,
        PyramidP2LocalCoefficients,
        PyramidP2LocalInterpolation<PyramidP2LocalBasis<D,R> > > Traits;

    /** \todo Please doc me !
     */
    PyramidP2LocalFiniteElement ()
    {
      gt.makePyramid();
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

    PyramidP2LocalFiniteElement* clone () const
    {
      return new PyramidP2LocalFiniteElement(*this);
    }

  private:
    PyramidP2LocalBasis<D,R> basis;
    PyramidP2LocalCoefficients coefficients;
    PyramidP2LocalInterpolation<PyramidP2LocalBasis<D,R> > interpolation;
    GeometryType gt;
  };

}

#endif
