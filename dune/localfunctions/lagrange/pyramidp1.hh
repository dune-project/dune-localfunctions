// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PYRAMID_P1_LOCALFINITEELEMENT_HH
#define DUNE_PYRAMID_P1_LOCALFINITEELEMENT_HH

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>

#include "pyramidp1/pyramidp1localbasis.hh"
#include "pyramidp1/pyramidp1localcoefficients.hh"
#include "pyramidp1/pyramidp1localinterpolation.hh"

namespace Dune
{

  /** \brief First-order Lagrangian finite element on a prism
   */
  template<class D, class R>
  class PyramidP1LocalFiniteElement
  {
  public:
    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<PyramidP1LocalBasis<D,R>,PyramidP1LocalCoefficients,
        PyramidP1LocalInterpolation<PyramidP1LocalBasis<D,R> > > Traits;



    /** \todo Please doc me !
     */
    PyramidP1LocalFiniteElement ()
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

    PyramidP1LocalFiniteElement* clone () const
    {
      return new PyramidP1LocalFiniteElement(*this);
    }

  private:
    PyramidP1LocalBasis<D,R> basis;
    PyramidP1LocalCoefficients coefficients;
    PyramidP1LocalInterpolation<PyramidP1LocalBasis<D,R> > interpolation;
    GeometryType gt;
  };

}

#endif
