// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS0_CUBE3D_LOCALFINITEELEMENT_HH
#define DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS0_CUBE3D_LOCALFINITEELEMENT_HH

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include "raviartthomas0cube3d/raviartthomas0cube3dall.hh"

namespace Dune
{
  /**
   * \brief Zero order Raviart-Thomas shape functions on cubes.
   *
   * \ingroup RaviartThomas
   *
   * \tparam D Type to represent the field in the domain.
   * \tparam R Type to represent the field in the range.
   */
  template<class D, class R>
  class RT0Cube3DLocalFiniteElement
  {
  public:
    typedef LocalFiniteElementTraits<
        RT0Cube3DLocalBasis<D,R>,
        RT0Cube3DLocalCoefficients,
        RT0Cube3DLocalInterpolation<RT0Cube3DLocalBasis<D,R> > > Traits;

    RT0Cube3DLocalFiniteElement ()
    {}

    RT0Cube3DLocalFiniteElement (int s) :
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
    RT0Cube3DLocalBasis<D,R> basis;
    RT0Cube3DLocalCoefficients coefficients;
    RT0Cube3DLocalInterpolation<RT0Cube3DLocalBasis<D,R> > interpolation;
  };
}
#endif // DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS0_CUBE3D_LOCALFINITEELEMENT_HH
