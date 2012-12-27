// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_RAVIARTTHOMAS02DLOCALFINITEELEMENT_HH
#define DUNE_RAVIARTTHOMAS02DLOCALFINITEELEMENT_HH

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include "raviartthomas02d/raviartthomas02dlocalbasis.hh"
#include "raviartthomas02d/raviartthomas02dlocalcoefficients.hh"
#include "raviartthomas02d/raviartthomas02dlocalinterpolation.hh"

#warning This header is deprecated, please use\
  dune/localfunctions/raviartthomas/raviartthomassimplex.hh instead

namespace Dune
{

  /**
   * \deprecated This class is deprecated and will be removed after Dune 2.3.
   *             Use RaviartThomasSimplexLocalFiniteElement<2,D,R> instead.
   */
  template<class D, class R>
  class
  DUNE_DEPRECATED_MSG("Use RaviartThomasSimplexLocalFiniteElement<2,D,R> instead")
  RT02DLocalFiniteElement
  {
  public:
    typedef LocalFiniteElementTraits<RT02DLocalBasis<D,R>,RT02DLocalCoefficients,
        RT02DLocalInterpolation<RT02DLocalBasis<D,R> > > Traits;

    RT02DLocalFiniteElement ()
    {
      gt.makeTriangle();
    }

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
