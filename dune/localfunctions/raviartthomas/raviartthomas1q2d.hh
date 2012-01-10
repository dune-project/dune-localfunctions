// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS1Q2DLOCALFINITEELEMENT_HH
#define DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS1Q2DLOCALFINITEELEMENT_HH

#include <dune/geometry/type.hh>

#include "../common/localfiniteelementtraits.hh"
#include "raviartthomas1q2d/raviartthomas1q2dlocalbasis.hh"
#include "raviartthomas1q2d/raviartthomas1q2dlocalinterpolation.hh"
#include "raviartthomas1q2d/raviartthomas1q2dlocalcoefficients.hh"

namespace Dune
{

  /**
   * \brief First order Raviart-Thomas shape functions on quadrilaterals.
   *
   * \tparam D Type to represent the field in the domain.
   * \tparam R Type to represent the field in the range.
   */
  template<class D, class R>
  class RT1Q2DLocalFiniteElement
  {

  public:
    typedef LocalFiniteElementTraits<RT1Q2DLocalBasis<D,R>,RT1Q2DLocalCoefficients,
        RT1Q2DLocalInterpolation<RT1Q2DLocalBasis<D,R> > > Traits;

    //! \brief Standard constructor
    RT1Q2DLocalFiniteElement ()
    {
      gt.makeQuadrilateral();
    }

    /**
     * \brief Make set number s, where 0 <= s < 16
     *
     * \param s Edge orientation indicator
     */
    RT1Q2DLocalFiniteElement (int s) : basis(s), interpolation(s)
    {
      gt.makeQuadrilateral();
    }

    const typename Traits::LocalBasisType& localBasis () const
    {
      return basis;
    }

    const typename Traits::LocalCoefficientsType& localCoefficients () const
    {
      return coefficients;
    }

    const typename Traits::LocalInterpolationType& localInterpolation () const
    {
      return interpolation;
    }

    GeometryType type () const
    {
      return gt;
    }

  private:
    RT1Q2DLocalBasis<D,R> basis;
    RT1Q2DLocalCoefficients coefficients;
    RT1Q2DLocalInterpolation<RT1Q2DLocalBasis<D,R> > interpolation;
    GeometryType gt;
  };
}
#endif // DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS1Q2DLOCALFINITEELEMENT_HH
