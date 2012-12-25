// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS0_CUBE3D_LOCALFINITEELEMENT_HH
#define DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS0_CUBE3D_LOCALFINITEELEMENT_HH

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include "raviartthomas0cube3d/raviartthomas0cube3dall.hh"

namespace Dune
{
  template<class D, class R>
  class RT0Cube3DLocalFiniteElement
  {
  public:
    typedef LocalFiniteElementTraits<
        RT0Cube3DLocalBasis<D,R>,
        RT0Cube3DLocalCoefficients,
        RT0Cube3DLocalInterpolation<RT0Cube3DLocalBasis<D,R> > > Traits;

    RT0Cube3DLocalFiniteElement ()
    {
      gt.makeHexahedron();
    }

    RT0Cube3DLocalFiniteElement (int s) : basis(s), interpolation(s)
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
    RT0Cube3DLocalBasis<D,R> basis;
    RT0Cube3DLocalCoefficients coefficients;
    RT0Cube3DLocalInterpolation<RT0Cube3DLocalBasis<D,R> > interpolation;
    GeometryType gt;
  };
}
#endif // DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS0_CUBE3D_LOCALFINITEELEMENT_HH
