// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_Q22DLOCALFINITEELEMENT_HH
#define DUNE_Q22DLOCALFINITEELEMENT_HH

#include <dune/common/geometrytype.hh>

#include "common/localfiniteelement.hh"
#include "q22d/q22dlocalbasis.hh"
#include "q22d/q22dlocalcoefficients.hh"
#include "q22d/q22dlocalinterpolation.hh"

namespace Dune
{

  template<class D, class R>
  class Q22DLocalFiniteElement : LocalFiniteElementInterface<
                                     LocalFiniteElementTraits<Q22DLocalBasis<D,R>,Q22DLocalCoefficients,
                                         Q22DLocalInterpolation<Q22DLocalBasis<D,R> > >,
                                     Q22DLocalFiniteElement<D,R> >
  {
  public:
    typedef LocalFiniteElementTraits<Q22DLocalBasis<D,R>,Q22DLocalCoefficients,
        Q22DLocalInterpolation<Q22DLocalBasis<D,R> > > Traits;

    Q22DLocalFiniteElement ()
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
    Q22DLocalBasis<D,R> basis;
    Q22DLocalCoefficients coefficients;
    Q22DLocalInterpolation<Q22DLocalBasis<D,R> > interpolation;
    GeometryType gt;
  };

}

#endif
