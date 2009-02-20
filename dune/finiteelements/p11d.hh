// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_P11DLOCALFINITEELEMENT_HH
#define DUNE_P11DLOCALFINITEELEMENT_HH

#include <dune/common/geometrytype.hh>

#include "common/localfiniteelement.hh"
#include "p11d/p11dlocalbasis.hh"
#include "p11d/p11dlocalcoefficients.hh"
#include "p11d/p11dlocalinterpolation.hh"

namespace Dune
{

  template<class D, class R>
  class P11DLocalFiniteElement : LocalFiniteElementInterface<
                                     LocalFiniteElementTraits<P11DLocalBasis<D,R>,P11DLocalCoefficients,
                                         P11DLocalInterpolation<P11DLocalBasis<D,R> > >,
                                     P11DLocalFiniteElement<D,R> >
  {
  public:
    typedef LocalFiniteElementTraits<P11DLocalBasis<D,R>,P11DLocalCoefficients,
        P11DLocalInterpolation<P11DLocalBasis<D,R> > > Traits;
    P11DLocalFiniteElement ()
    {
      gt.makeLine();
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
    P11DLocalBasis<D,R> basis;
    P11DLocalCoefficients coefficients;
    P11DLocalInterpolation<P11DLocalBasis<D,R> > interpolation;
    GeometryType gt;
  };

}

#endif
