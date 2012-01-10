// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI12DLOCALFINITEELEMENT_HH
#define DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI12DLOCALFINITEELEMENT_HH

#include <dune/geometry/type.hh>

#include "../common/localfiniteelementtraits.hh"
#include "brezzidouglasmarini12d/brezzidouglasmarini12dlocalbasis.hh"
#include "brezzidouglasmarini12d/brezzidouglasmarini12dlocalcoefficients.hh"
#include "brezzidouglasmarini12d/brezzidouglasmarini12dlocalinterpolation.hh"

namespace Dune
{

  /**
   * \brief First order Brezzi-Douglas-Marini shape functions on triangles.
   *
   * \tparam D Type to represent the field in the domain.
   * \tparam R Type to represent the field in the range.
   */
  template<class D, class R>
  class BDM12DLocalFiniteElement
  {

  public:
    typedef LocalFiniteElementTraits<BDM12DLocalBasis<D,R>,BDM12DLocalCoefficients,
        BDM12DLocalInterpolation<BDM12DLocalBasis<D,R> > > Traits;

    //! \brief Standard constructor
    BDM12DLocalFiniteElement ()
    {
      gt.makeTriangle();
    }

    /**
     * \brief Make set number s, where 0 <= s < 8
     *
     * \param s Edge orientation indicator
     */
    BDM12DLocalFiniteElement (int s) : basis(s), interpolation(s)
    {
      gt.makeTriangle();
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
    BDM12DLocalBasis<D,R> basis;
    BDM12DLocalCoefficients coefficients;
    BDM12DLocalInterpolation<BDM12DLocalBasis<D,R> > interpolation;
    GeometryType gt;
  };
}
#endif // DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI12DLOCALFINITEELEMENT_HH
