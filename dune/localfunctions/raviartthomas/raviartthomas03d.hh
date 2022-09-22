// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS_RAVIARTTHOMAS03D_HH
#define DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS_RAVIARTTHOMAS03D_HH

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include "raviartthomas03d/raviartthomas03dlocalbasis.hh"
#include "raviartthomas03d/raviartthomas03dlocalcoefficients.hh"
#include "raviartthomas03d/raviartthomas03dlocalinterpolation.hh"

namespace Dune
{

  /**
   * \brief Zero order Raviart-Thomas shape functions on tetrahedra.
   *
   * \ingroup RaviartThomas
   *
   * \tparam D Type to represent the field in the domain.
   * \tparam R Type to represent the field in the range.
   */
  template<class D, class R>
  class
  RT03DLocalFiniteElement
  {
  public:
    typedef LocalFiniteElementTraits<RT03DLocalBasis<D,R>,RT03DLocalCoefficients,
        RT03DLocalInterpolation<RT03DLocalBasis<D,R> > > Traits;

    //! \brief Standard constructor
    RT03DLocalFiniteElement ()
    {}

    /**
     * \brief Constructor with explicitly given face orientations
     *
     * \param s Face orientation indicator
     */
    RT03DLocalFiniteElement (std::bitset<4> s) :
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
      return 4;
    }

    static constexpr GeometryType type ()
    {
      return GeometryTypes::tetrahedron;
    }

  private:
    RT03DLocalBasis<D,R> basis;
    RT03DLocalCoefficients coefficients;
    RT03DLocalInterpolation<RT03DLocalBasis<D,R> > interpolation;
  };

}

#endif
