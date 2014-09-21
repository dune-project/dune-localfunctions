// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS0_CUBE2D_LOCALFINITEELEMENT_HH
#define DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS0_CUBE2D_LOCALFINITEELEMENT_HH

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include "raviartthomas0cube2d/raviartthomas0cube2dall.hh"

namespace Dune
{
  template<class D, class R>
  class RT0Cube2DLocalFiniteElement
  {
  public:
    typedef LocalFiniteElementTraits<
        RT0Cube2DLocalBasis<D,R>,
        RT0Cube2DLocalCoefficients,
        RT0Cube2DLocalInterpolation<RT0Cube2DLocalBasis<D,R> > > Traits;

    RT0Cube2DLocalFiniteElement ()
    {
      gt.makeQuadrilateral();
    }

    RT0Cube2DLocalFiniteElement (int s) : basis(s), interpolation(s)
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

    /** \brief Number of shape functions in this finite element */
    uint size () const
    {
      return basis.size();
    }

    GeometryType type () const
    {
      return gt;
    }

  private:
    RT0Cube2DLocalBasis<D,R> basis;
    RT0Cube2DLocalCoefficients coefficients;
    RT0Cube2DLocalInterpolation<RT0Cube2DLocalBasis<D,R> > interpolation;
    GeometryType gt;
  };
}
#endif // DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS0_CUBE2D_LOCALFINITEELEMENT_HH
