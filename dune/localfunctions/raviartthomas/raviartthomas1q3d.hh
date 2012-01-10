// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS1Q3DLOCALFINITEELEMENT_HH
#define DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS1Q3DLOCALFINITEELEMENT_HH

#include <dune/geometry/type.hh>

#include "../common/localfiniteelementtraits.hh"
#include "raviartthomas1q3d/raviartthomas1q3dlocalbasis.hh"
#include "raviartthomas1q3d/raviartthomas1q3dlocalcoefficients.hh"
#include "raviartthomas1q3d/raviartthomas1q3dlocalinterpolation.hh"

namespace Dune
{
  /**
   * \brief First order Raviart-Thomas shape functions on cubes.
   *
   * \tparam D Type to represent the field in the domain.
   * \tparam R Type to represent the field in the range.
   */
  template<class D, class R>
  class RT1Q3DLocalFiniteElement
  {

  public:
    typedef LocalFiniteElementTraits<RT1Q3DLocalBasis<D,R>,RT1Q3DLocalCoefficients,
        RT1Q3DLocalInterpolation<RT1Q3DLocalBasis<D,R> > > Traits;

    //! \brief Standard constructor
    RT1Q3DLocalFiniteElement ()
    {
      gt.makeHexahedron();
    }

    /**
     * \brief Make set number s, where 0 <= s < 64
     *
     * \param s Edge orientation indicator
     */
    RT1Q3DLocalFiniteElement (int s) : basis(s), interpolation(s)
    {
      gt.makeHexahedron();
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
    RT1Q3DLocalBasis<D,R> basis;
    RT1Q3DLocalCoefficients coefficients;
    RT1Q3DLocalInterpolation<RT1Q3DLocalBasis<D,R> > interpolation;
    GeometryType gt;
  };
}
#endif // DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS1Q3DLOCALFINITEELEMENT_HH
