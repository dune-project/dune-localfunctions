// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS12DLOCALFINITEELEMENT_HH
#define DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS12DLOCALFINITEELEMENT_HH

#include <dune/geometry/type.hh>

#include "../common/localfiniteelementtraits.hh"
#include "raviartthomas12d/raviartthomas12dlocalbasis.hh"
#include "raviartthomas12d/raviartthomas12dlocalcoefficients.hh"
#include "raviartthomas12d/raviartthomas12dlocalinterpolation.hh"

namespace Dune
{

  /**
   * \brief First order Raviart-Thomas shape functions on triangles.
   *
   * \ingroup RaviartThomas
   *
   * \tparam D Type to represent the field in the domain.
   * \tparam R Type to represent the field in the range.
   */
  template<class D, class R>
  class
  RT12DLocalFiniteElement
  {

  public:
    typedef LocalFiniteElementTraits<RT12DLocalBasis<D,R>,RT12DLocalCoefficients,
        RT12DLocalInterpolation<RT12DLocalBasis<D,R> > > Traits;

    //! \brief Standard constructor
    RT12DLocalFiniteElement ()
    {}

    /**
     * \brief Make set number s, where 0 <= s < 8
     *
     * \param s Edge orientation indicator
     */
    RT12DLocalFiniteElement (int s) :
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
      return GeometryTypes::triangle;
    }

  private:
    RT12DLocalBasis<D,R> basis;
    RT12DLocalCoefficients coefficients;
    RT12DLocalInterpolation<RT12DLocalBasis<D,R> > interpolation;
  };
}

#endif // DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS12DLOCALFINITEELEMENT_HH
