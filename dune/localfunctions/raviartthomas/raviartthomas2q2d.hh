// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS2Q2DLOCALFINITEELEMENT_HH
#define DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS2Q2DLOCALFINITEELEMENT_HH

#include <dune/geometry/type.hh>

#include "../common/localfiniteelementtraits.hh"
#include "raviartthomas2q2d/raviartthomas2q2dlocalbasis.hh"
#include "raviartthomas2q2d/raviartthomas2q2dlocalcoefficients.hh"
#include "raviartthomas2q2d/raviartthomas2q2dlocalinterpolation.hh"

namespace Dune
{
  /**
   * \brief Second order Raviart-Thomas shape functions on cubes.
   *
   * \tparam D Type to represent the field in the domain.
   * \tparam R Type to represent the field in the range.
   */
  template<class D, class R>
  class RT2Q2DLocalFiniteElement
  {

  public:
    typedef LocalFiniteElementTraits<RT2Q2DLocalBasis<D,R>,RT2Q2DLocalCoefficients,
        RT2Q2DLocalInterpolation<RT2Q2DLocalBasis<D,R> > > Traits;

    //! \brief Standard constructor
    RT2Q2DLocalFiniteElement ()
    {
      gt.makeQuadrilateral();
    }

    /**
     * \brief Make set number s, where 0 <= s < 16
     *
     * \param s Edge orientation indicator
     */
    RT2Q2DLocalFiniteElement (int s) : basis(s), interpolation(s)
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
    RT2Q2DLocalBasis<D,R> basis;
    RT2Q2DLocalCoefficients coefficients;
    RT2Q2DLocalInterpolation<RT2Q2DLocalBasis<D,R> > interpolation;
    GeometryType gt;
  };
}
#endif // DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS2Q2DLOCALFINITEELEMENT_HH
