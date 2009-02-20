// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_Q13DLOCALFINITEELEMENT_HH
#define DUNE_Q13DLOCALFINITEELEMENT_HH

#include <dune/common/geometrytype.hh>

#include "common/localfiniteelement.hh"
#include "q13d/q13dlocalbasis.hh"
#include "q13d/q13dlocalcoefficients.hh"
#include "q13d/q13dlocalinterpolation.hh"

namespace Dune
{

  template<class D, class R>
  class Q13DLocalFiniteElement : LocalFiniteElementInterface<
                                     LocalFiniteElementTraits<Q13DLocalBasis<D,R>,Q13DLocalCoefficients,
                                         Q13DLocalInterpolation<Q13DLocalBasis<D,R> > >,
                                     Q13DLocalFiniteElement<D,R> >
  {
  public:
    typedef LocalFiniteElementTraits<Q13DLocalBasis<D,R>,Q13DLocalCoefficients,
        Q13DLocalInterpolation<Q13DLocalBasis<D,R> > > Traits;

    Q13DLocalFiniteElement ()
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
    Q13DLocalBasis<D,R> basis;
    Q13DLocalCoefficients coefficients;
    Q13DLocalInterpolation<Q13DLocalBasis<D,R> > interpolation;
    GeometryType gt;
  };

}

#endif
