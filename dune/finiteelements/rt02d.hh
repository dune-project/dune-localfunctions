// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_RT02DLOCALFINITEELEMENT_HH
#define DUNE_RT02DLOCALFINITEELEMENT_HH

#include <dune/common/geometrytype.hh>

#include "common/localfiniteelement.hh"
#include "rt02d/rt02dlocalbasis.hh"
#include "rt02d/rt02dlocalcoefficients.hh"
#include "rt02d/rt02dlocalinterpolation.hh"

namespace Dune
{

  template<class D, class R>
  class RT02DLocalFiniteElement : LocalFiniteElementInterface<
                                      LocalFiniteElementTraits<RT02DLocalBasis<D,R>,RT02DLocalCoefficients,
                                          RT02DLocalInterpolation<RT02DLocalBasis<D,R> > >
#ifndef DUNE_VIRTUAL_SHAPEFUNCTIONS
                                      , RT02DLocalFiniteElement<D,R>
#endif
                                      >
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
