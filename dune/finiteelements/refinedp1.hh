// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_REFINED_P1_LOCALFINITEELEMENT_HH
#define DUNE_REFINED_P1_LOCALFINITEELEMENT_HH

#include <dune/common/geometrytype.hh>

#include "common/localfiniteelement.hh"
#include "refinedp1/refinedp1localbasis.hh"
#include "pk2d/pk2dlocalcoefficients.hh"
#include "pk2d/pk2dlocalinterpolation.hh"

namespace Dune
{

  template<class D, class R>
  class RefinedP1LocalFiniteElement : LocalFiniteElementInterface<
                                          LocalFiniteElementTraits<RefinedP1LocalBasis<D,R>,
                                              Pk2DLocalCoefficients<2>,
                                              Pk2DLocalInterpolation<Pk2DLocalBasis<D,R,2> > >
#ifndef DUNE_VIRTUAL_SHAPEFUNCTIONS
                                          , RefinedP1LocalFiniteElement<D,R>
#endif
                                          >
  {
  public:
    typedef LocalFiniteElementTraits<RefinedP1LocalBasis<D,R>,
        Pk2DLocalCoefficients<2>,
        Pk2DLocalInterpolation<Pk2DLocalBasis<D,R,2> > > Traits;
    RefinedP1LocalFiniteElement ()
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
    RefinedP1LocalBasis<D,R> basis;
    Pk2DLocalCoefficients<2> coefficients;
    Pk2DLocalInterpolation<Pk2DLocalBasis<D,R,2> > interpolation;
    GeometryType gt;
  };

}

#endif
