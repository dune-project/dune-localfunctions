// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_EDGER02DLOCALFINITEELEMENT_HH
#define DUNE_EDGER02DLOCALFINITEELEMENT_HH

#include <dune/common/geometrytype.hh>

#include "common/localfiniteelementtraits.hh"
#include "edger02d/edger02dlocalbasis.hh"
#include "edger02d/edger02dlocalcoefficients.hh"
#include "edger02d/edger02dlocalinterpolation.hh"

namespace Dune
{

  template<class D, class R>
  class EdgeR02DLocalFiniteElement
  {
  public:
    typedef LocalFiniteElementTraits<
        EdgeR02DLocalBasis<D,R>,
        EdgeR02DLocalCoefficients,
        EdgeR02DLocalInterpolation<EdgeR02DLocalBasis<D,R> >
        > Traits;

    EdgeR02DLocalFiniteElement ()
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
    EdgeR02DLocalBasis<D,R> basis;
    EdgeR02DLocalCoefficients coefficients;
    EdgeR02DLocalInterpolation<EdgeR02DLocalBasis<D,R> > interpolation;
    GeometryType gt;
  };

}

#endif //DUNE_EDGER02DLOCALFINITEELEMENT_HH
