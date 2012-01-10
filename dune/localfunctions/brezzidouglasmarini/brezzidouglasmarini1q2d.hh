// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI1Q2DLOCALFINITEELEMENT_HH
#define DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI1Q2DLOCALFINITEELEMENT_HH

#include <dune/geometry/type.hh>

#include "../common/localfiniteelementtraits.hh"
#include "brezzidouglasmarini1q2d/brezzidouglasmarini1q2dlocalbasis.hh"
#include "brezzidouglasmarini1q2d/brezzidouglasmarini1q2dlocalcoefficients.hh"
#include "brezzidouglasmarini1q2d/brezzidouglasmarini1q2dlocalinterpolation.hh"

namespace Dune
{
  /**
   * \brief First order Brezzi-Douglas-Marini shape functions on quadrilaterals.
   *
   * \tparam D Type to represent the field in the domain.
   * \tparam R Type to represent the field in the range.
   */
  template<class D, class R>
  class BDM1Q2DLocalFiniteElement
  {

  public:
    typedef LocalFiniteElementTraits<BDM1Q2DLocalBasis<D,R>,BDM1Q2DLocalCoefficients,
        BDM1Q2DLocalInterpolation<BDM1Q2DLocalBasis<D,R> > > Traits;

    //! \brief Standard constructor
    BDM1Q2DLocalFiniteElement ()
    {
      gt.makeQuadrilateral();
    }

    /**
     * \brief Make set number s, where 0 <= s < 16
     *
     * \param s Edge orientation indicator
     */
    BDM1Q2DLocalFiniteElement (int s) : basis(s), interpolation(s)
    {
      gt.makeQuadrilateral();
    }

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

    GeometryType type () const
    {
      return gt;
    }

  private:
    BDM1Q2DLocalBasis<D,R> basis;
    BDM1Q2DLocalCoefficients coefficients;
    BDM1Q2DLocalInterpolation<BDM1Q2DLocalBasis<D,R> > interpolation;
    GeometryType gt;
  };
}
#endif // DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI1Q2DLOCALFINITEELEMENT_HH
