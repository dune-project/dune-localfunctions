// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS1_CUBE3D_LOCALFINITEELEMENT_HH
#define DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS1_CUBE3D_LOCALFINITEELEMENT_HH

#include <dune/geometry/type.hh>

#include "../common/localfiniteelementtraits.hh"
#include "raviartthomas1cube3d/raviartthomas1cube3dlocalbasis.hh"
#include "raviartthomas1cube3d/raviartthomas1cube3dlocalcoefficients.hh"
#include "raviartthomas1cube3d/raviartthomas1cube3dlocalinterpolation.hh"

namespace Dune
{
  /**
   * \brief First order Raviart-Thomas shape functions on cubes.
   *
   * \ingroup RaviartThomas
   *
   * \tparam D Type to represent the field in the domain.
   * \tparam R Type to represent the field in the range.
   */
  template<class D, class R>
  class RT1Cube3DLocalFiniteElement
  {

  public:
    typedef LocalFiniteElementTraits<
        RT1Cube3DLocalBasis<D,R>,
        RT1Cube3DLocalCoefficients,
        RT1Cube3DLocalInterpolation<RT1Cube3DLocalBasis<D,R> > > Traits;

    //! \brief Standard constructor
    RT1Cube3DLocalFiniteElement ()
    {}

    /**
     * \brief Make set number s, where 0 <= s < 64
     *
     * \param s Edge orientation indicator
     */
    RT1Cube3DLocalFiniteElement (int s) :
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
      return GeometryTypes::hexahedron;
    }

  private:
    RT1Cube3DLocalBasis<D,R> basis;
    RT1Cube3DLocalCoefficients coefficients;
    RT1Cube3DLocalInterpolation<RT1Cube3DLocalBasis<D,R> > interpolation;
  };
}
#endif // DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS1_CUBE3D_LOCALFINITEELEMENT_HH
