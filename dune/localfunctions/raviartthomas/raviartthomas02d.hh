// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
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
    {
      gt.makeTriangle();
    }

    /**
     * \brief Make set number s, where 0 <= s < 8
     *
     * \param s Edge orientation indicator
     */
    RT02DLocalFiniteElement (int s) : basis(s), interpolation(s)
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
    RT02DLocalBasis<D,R> basis;
    RT02DLocalCoefficients coefficients;
    RT02DLocalInterpolation<RT02DLocalBasis<D,R> > interpolation;
    GeometryType gt;
  };

}

#endif
