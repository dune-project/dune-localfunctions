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
