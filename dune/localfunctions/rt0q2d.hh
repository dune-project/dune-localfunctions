// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_RT0Q2DLOCALFINITEELEMENT_HH
#define DUNE_RT0Q2DLOCALFINITEELEMENT_HH

#include <dune/common/geometrytype.hh>

#include "common/localfiniteelement.hh"
#include "rt0q2d/rt0q2dall.hh"

namespace Dune
{
  template<class D, class R>
  class RT0Q2DLocalFiniteElement : LocalFiniteElementInterface<
                                       LocalFiniteElementTraits<RT0Q2DLocalBasis<D,R>,RT0Q2DLocalCoefficients,
                                           RT0Q2DLocalInterpolation<RT0Q2DLocalBasis<D,R> > >
#ifndef DUNE_VIRTUAL_SHAPEFUNCTIONS
                                       , RT0Q2DLocalFiniteElement<D,R>
#endif
                                       >
  {
  public:
    typedef LocalFiniteElementTraits<RT0Q2DLocalBasis<D,R>,RT0Q2DLocalCoefficients,
        RT0Q2DLocalInterpolation<RT0Q2DLocalBasis<D,R> > > Traits;

    RT0Q2DLocalFiniteElement ()
    {
      gt.makeQuadrilateral();
    }

    RT0Q2DLocalFiniteElement (int s) : basis(s), interpolation(s)
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
    RT0Q2DLocalBasis<D,R> basis;
    RT0Q2DLocalCoefficients coefficients;
    RT0Q2DLocalInterpolation<RT0Q2DLocalBasis<D,R> > interpolation;
    GeometryType gt;
  };
}
#endif
