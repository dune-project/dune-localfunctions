// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS3_CUBE2D_LOCALFINITEELEMENT_HH
#define DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS3_CUBE2D_LOCALFINITEELEMENT_HH

#include <dune/geometry/type.hh>

#include "../common/localfiniteelementtraits.hh"
#include "raviartthomas3cube2d/raviartthomas3cube2dlocalbasis.hh"
#include "raviartthomas3cube2d/raviartthomas3cube2dlocalcoefficients.hh"
#include "raviartthomas3cube2d/raviartthomas3cube2dlocalinterpolation.hh"

namespace Dune
{
  /**
   * \brief Second order Raviart-Thomas shape functions on cubes.
   *
   * \note The Jacobian is not implemented.
   *
   * \ingroup RaviartThomas
   *
   * \tparam D Type to represent the field in the domain.
   * \tparam R Type to represent the field in the range.
   */
  template<class D, class R>
  class RT3Cube2DLocalFiniteElement
  {

  public:
    typedef LocalFiniteElementTraits<
        RT3Cube2DLocalBasis<D,R>,
        RT3Cube2DLocalCoefficients,
        RT3Cube2DLocalInterpolation<RT3Cube2DLocalBasis<D,R> > > Traits;

    //! \brief Standard constructor
    RT3Cube2DLocalFiniteElement ()
    {}

    /**
     * \brief Make set number s, where 0 <= s < 16
     *
     * \param s Edge orientation indicator
     */
    RT3Cube2DLocalFiniteElement (int s) :
      basis(s),
      interpolation(s)
    {}

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

    /** \brief Number of shape functions in this finite element */
    unsigned int size () const
    {
      return basis.size();
    }

    static constexpr GeometryType type ()
    {
      return GeometryTypes::quadrilateral;
    }

  private:
    RT3Cube2DLocalBasis<D,R> basis;
    RT3Cube2DLocalCoefficients coefficients;
    RT3Cube2DLocalInterpolation<RT3Cube2DLocalBasis<D,R> > interpolation;
  };
}
#endif // DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS3_CUBE2D_LOCALFINITEELEMENT_HH
