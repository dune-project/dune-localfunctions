// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS0_CUBE2D_LOCALFINITEELEMENT_HH
#define DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS0_CUBE2D_LOCALFINITEELEMENT_HH

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include "raviartthomas0cube2d/raviartthomas0cube2dall.hh"

namespace Dune
{
  /**
   * \brief Zero order Raviart-Thomas shape functions on rectangles.
   *
   * \ingroup RaviartThomas
   *
   * \tparam D Type to represent the field in the domain.
   * \tparam R Type to represent the field in the range.
   */
  template<class D, class R>
  class RT0Cube2DLocalFiniteElement
  {
  public:
    typedef LocalFiniteElementTraits<
        RT0Cube2DLocalBasis<D,R>,
        RT0Cube2DLocalCoefficients,
        RT0Cube2DLocalInterpolation<RT0Cube2DLocalBasis<D,R> > > Traits;

    RT0Cube2DLocalFiniteElement ()
    {}

    RT0Cube2DLocalFiniteElement (int s) :
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
    RT0Cube2DLocalBasis<D,R> basis;
    RT0Cube2DLocalCoefficients coefficients;
    RT0Cube2DLocalInterpolation<RT0Cube2DLocalBasis<D,R> > interpolation;
  };
}
#endif // DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS0_CUBE2D_LOCALFINITEELEMENT_HH
