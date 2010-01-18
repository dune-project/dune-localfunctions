// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_RAVIARTTHOMAS0Q3DLOCALFINITEELEMENT_HH
#define DUNE_RAVIARTTHOMAS0Q3DLOCALFINITEELEMENT_HH

#include <dune/common/geometrytype.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include "raviartthomas0q3d/raviartthomas0q3dall.hh"

namespace Dune
{
  template<class D, class R>
  class RT0Q3DLocalFiniteElement
  {
  public:
    typedef LocalFiniteElementTraits<RT0Q3DLocalBasis<D,R>,RT0Q3DLocalCoefficients,
        RT0Q3DLocalInterpolation<RT0Q3DLocalBasis<D,R> > > Traits;

    RT0Q3DLocalFiniteElement ()
    {
      gt.makeHexahedron();
    }

    RT0Q3DLocalFiniteElement (int s) : basis(s), interpolation(s)
    {
      gt.makeHexahedron();
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
    RT0Q3DLocalBasis<D,R> basis;
    RT0Q3DLocalCoefficients coefficients;
    RT0Q3DLocalInterpolation<RT0Q3DLocalBasis<D,R> > interpolation;
    GeometryType gt;
  };
}
#endif
