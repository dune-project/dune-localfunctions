// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_P13DLOCALFINITEELEMENT_HH
#define DUNE_P13DLOCALFINITEELEMENT_HH

#include <dune/common/geometrytype.hh>

#include "common/localfiniteelement.hh"
#include "p13d/p13dlocalbasis.hh"
#include "p13d/p13dlocalcoefficients.hh"
#include "p13d/p13dlocalinterpolation.hh"

namespace Dune
{

  template<class D, class R>
  class P13DLocalFiniteElement : LocalFiniteElementInterface<
                                     LocalFiniteElementTraits<P13DLocalBasis<D,R>,P13DLocalCoefficients,
                                         P13DLocalInterpolation<P13DLocalBasis<D,R> > >,
                                     P13DLocalFiniteElement<D,R> >
  {
  public:
    typedef LocalFiniteElementTraits<P13DLocalBasis<D,R>,P13DLocalCoefficients,
        P13DLocalInterpolation<P13DLocalBasis<D,R> > > Traits;
    P13DLocalFiniteElement ()
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
    P13DLocalBasis<D,R> basis;
    P13DLocalCoefficients coefficients;
    P13DLocalInterpolation<P13DLocalBasis<D,R> > interpolation;
    GeometryType gt;
  };

}

#endif
