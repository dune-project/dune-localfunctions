// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI2_QUBE2D_LOCALFINITEELEMENT_HH
#define DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI2_QUBE2D_LOCALFINITEELEMENT_HH

#include <dune/geometry/type.hh>

#include "../common/localfiniteelementtraits.hh"
#include "brezzidouglasmarini2cube2d/brezzidouglasmarini2cube2dlocalbasis.hh"
#include "brezzidouglasmarini2cube2d/brezzidouglasmarini2cube2dlocalcoefficients.hh"
#include "brezzidouglasmarini2cube2d/brezzidouglasmarini2cube2dlocalinterpolation.hh"

namespace Dune
{
  /**
   * \brief Second order Brezzi-Douglas-Marini shape functions on quadrilaterals.
   *
   * \ingroup BrezziDouglasMarini
   *
   * \tparam D Type to represent the field in the domain.
   * \tparam R Type to represent the field in the range.
   */
  template<class D, class R>
  class BDM2Cube2DLocalFiniteElement
  {

  public:
    typedef LocalFiniteElementTraits<
        BDM2Cube2DLocalBasis<D,R>,
        BDM2Cube2DLocalCoefficients,
        BDM2Cube2DLocalInterpolation<BDM2Cube2DLocalBasis<D,R> > > Traits;

    //! \brief Standard constructor
    BDM2Cube2DLocalFiniteElement ()
    {}

    /**
     * \brief Make set number s, where 0 <= s < ??
     *
     * \param s Edge orientation indicator
     */
    BDM2Cube2DLocalFiniteElement (int s) :
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
    BDM2Cube2DLocalBasis<D,R> basis;
    BDM2Cube2DLocalCoefficients coefficients;
    BDM2Cube2DLocalInterpolation<BDM2Cube2DLocalBasis<D,R> > interpolation;
  };
}
#endif // DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI2_QUBE2D_LOCALFINITEELEMENT_HH
