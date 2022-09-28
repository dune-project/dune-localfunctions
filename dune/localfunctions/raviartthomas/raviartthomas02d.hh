// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_RAVIARTTHOMAS02DLOCALFINITEELEMENT_HH
#define DUNE_RAVIARTTHOMAS02DLOCALFINITEELEMENT_HH

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include "raviartthomas02d/raviartthomas02dlocalbasis.hh"
#include "raviartthomas02d/raviartthomas02dlocalcoefficients.hh"
#include "raviartthomas02d/raviartthomas02dlocalinterpolation.hh"

namespace Dune
{

  /**
   * \brief Zero order Raviart-Thomas shape functions on triangles.
   *
   * \ingroup RaviartThomas
   *
   * \tparam D Type to represent the field in the domain.
   * \tparam R Type to represent the field in the range.
   */
  template<class D, class R>
  class
  RT02DLocalFiniteElement
  {
  public:
    typedef LocalFiniteElementTraits<RT02DLocalBasis<D,R>,RT02DLocalCoefficients,
        RT02DLocalInterpolation<RT02DLocalBasis<D,R> > > Traits;

    //! \brief Standard constructor
    RT02DLocalFiniteElement ()
    {}

    /**
     * \brief Constructor with explicitly given edge orientations
     *
     * \param s Edge orientation indicator
     */
    RT02DLocalFiniteElement (std::bitset<3> s) :
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

    unsigned int size () const
    {
      return 3;
    }

    static constexpr GeometryType type ()
    {
      return GeometryTypes::triangle;
    }

  private:
    RT02DLocalBasis<D,R> basis;
    RT02DLocalCoefficients coefficients;
    RT02DLocalInterpolation<RT02DLocalBasis<D,R> > interpolation;
  };

}

#endif
