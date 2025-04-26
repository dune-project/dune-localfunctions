// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS0_PYRAMID_HH
#define DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS0_PYRAMID_HH

#include <dune/geometry/type.hh>

#include "../common/localfiniteelementtraits.hh"
#include "raviartthomas0pyramid/raviartthomas0pyramidlocalbasis.hh"
#include "raviartthomas0pyramid/raviartthomas0pyramidlocalcoefficients.hh"
#include "raviartthomas0pyramid/raviartthomas0pyramidlocalinterpolation.hh"

namespace Dune
{
  /**
   * \brief First order Raviart-Thomas shape functions on pyramids.
   *
   * This implements the composite element proposed by Ainsworth/Fu in
   * "A lowest-order composite finite element exact sequence on pyramids"
   * Computer Methods in Applied Mechanics and Engineering,
   * Volume 324, 1 September 2017, Pages 110-127.
   * https://arxiv.org/abs/1705.00064
   *
   * Notice that this is a composite element with a discontinuity across
   * the plane with x[0]==x[1]. It contains on DOF per face and one internal
   * DOF located at the center of the intersections of the two sub-elements.
   *
   * \ingroup RaviartThomas
   *
   * \tparam D Type to represent the field in the domain.
   * \tparam R Type to represent the field in the range.
   */
  template<class D, class R>
  class RT0PyramidLocalFiniteElement
  {

  public:
    typedef LocalFiniteElementTraits<
        RT0PyramidLocalBasis<D,R>,
        RT0PyramidLocalCoefficients,
        RT0PyramidLocalInterpolation<RT0PyramidLocalBasis<D,R> > > Traits;

    //! \brief Standard constructor
    RT0PyramidLocalFiniteElement ()
    {}

    /**
     * \brief Make set number s, where 0 <= s < 32
     *
     * \param s Face orientation indicator
     */
    RT0PyramidLocalFiniteElement (int s) :
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
      return GeometryTypes::pyramid;
    }

  private:
    RT0PyramidLocalBasis<D,R> basis;
    RT0PyramidLocalCoefficients coefficients;
    RT0PyramidLocalInterpolation<RT0PyramidLocalBasis<D,R> > interpolation;
  };
}
#endif // DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS0_PYRAMID_HH
