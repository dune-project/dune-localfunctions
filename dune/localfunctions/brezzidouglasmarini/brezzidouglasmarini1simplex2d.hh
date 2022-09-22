// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI1_SIMPLEX2D_LOCALFINITEELEMENT_HH
#define DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI1_SIMPLEX2D_LOCALFINITEELEMENT_HH

#include <dune/geometry/type.hh>

#include "../common/localfiniteelementtraits.hh"
#include "brezzidouglasmarini1simplex2d/brezzidouglasmarini1simplex2dlocalbasis.hh"
#include "brezzidouglasmarini1simplex2d/brezzidouglasmarini1simplex2dlocalcoefficients.hh"
#include "brezzidouglasmarini1simplex2d/brezzidouglasmarini1simplex2dlocalinterpolation.hh"

namespace Dune
{

  /**
   * \brief First order Brezzi-Douglas-Marini shape functions on triangles.
   *
   * \ingroup BrezziDouglasMarini
   *
   * \tparam D Type to represent the field in the domain.
   * \tparam R Type to represent the field in the range.
   */
  template<class D, class R>
  class BDM1Simplex2DLocalFiniteElement
  {

  public:
    typedef LocalFiniteElementTraits<
        BDM1Simplex2DLocalBasis<D,R>,
        BDM1Simplex2DLocalCoefficients,
        BDM1Simplex2DLocalInterpolation<BDM1Simplex2DLocalBasis<D,R> > > Traits;

    //! \brief Standard constructor
    BDM1Simplex2DLocalFiniteElement ()
    {}

    /**
     * \brief Make set number s, where 0 <= s < 8
     *
     * \param s Edge orientation indicator
     */
    BDM1Simplex2DLocalFiniteElement (int s) :
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
    BDM1Simplex2DLocalBasis<D,R> basis;
    BDM1Simplex2DLocalCoefficients coefficients;
    BDM1Simplex2DLocalInterpolation<BDM1Simplex2DLocalBasis<D,R> > interpolation;
  };
}
#endif // DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI1_SIMPLEX2D_LOCALFINITEELEMENT_HH
