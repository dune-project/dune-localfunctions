// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_Q12DLOCALFINITEELEMENT_HH
#define DUNE_Q12DLOCALFINITEELEMENT_HH

#include <dune/common/geometrytype.hh>

#include "common/localfiniteelement.hh"
#include "q12d/q12dlocalbasis.hh"
#include "q12d/q12dlocalcoefficients.hh"
#include "q12d/q12dlocalinterpolation.hh"

namespace Dune
{

  template<class D, class R>
  class Q12DLocalFiniteElement : LocalFiniteElementInterface<
                                     LocalFiniteElementTraits<Q12DLocalBasis<D,R>,Q12DLocalCoefficients,
                                         Q12DLocalInterpolation<Q12DLocalBasis<D,R> > >
#ifndef DUNE_VIRTUAL_SHAPEFUNCTIONS
                                     ,Q12DLocalFiniteElement<D,R>
#endif
                                     >
  {
  public:
    typedef LocalFiniteElementTraits<Q12DLocalBasis<D,R>,Q12DLocalCoefficients,
        Q12DLocalInterpolation<Q12DLocalBasis<D,R> > > Traits;

    Q12DLocalFiniteElement ()
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
    Q12DLocalBasis<D,R> basis;
    Q12DLocalCoefficients coefficients;
    Q12DLocalInterpolation<Q12DLocalBasis<D,R> > interpolation;
    GeometryType gt;
  };

}

#endif
