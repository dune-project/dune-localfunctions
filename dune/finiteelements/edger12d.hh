// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_EDGER12DLOCALFINITEELEMENT_HH
#define DUNE_EDGER12DLOCALFINITEELEMENT_HH

#include <dune/common/geometrytype.hh>

#include "common/localfiniteelement.hh"
#include "edger12d/edger12dlocalbasis.hh"
#include "edger12d/edger12dlocalcoefficients.hh"
#include "edger12d/edger12dlocalinterpolation.hh"

namespace Dune
{

  template<class D, class R>
  class EdgeR12DLocalFiniteElement
    : LocalFiniteElementInterface<
          LocalFiniteElementTraits<
              EdgeR12DLocalBasis<D,R>,
              EdgeR12DLocalCoefficients,
              EdgeR12DLocalInterpolation<EdgeR12DLocalBasis<D,R> >
              >,
          EdgeR12DLocalFiniteElement<D,R>
          >
  {
  public:
    typedef LocalFiniteElementTraits<
        EdgeR12DLocalBasis<D,R>,
        EdgeR12DLocalCoefficients,
        EdgeR12DLocalInterpolation<EdgeR12DLocalBasis<D,R> >
        > Traits;

    EdgeR12DLocalFiniteElement ()
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
    EdgeR12DLocalBasis<D,R> basis;
    EdgeR12DLocalCoefficients coefficients;
    EdgeR12DLocalInterpolation<EdgeR12DLocalBasis<D,R> > interpolation;
    GeometryType gt;
  };

}

#endif //DUNE_EDGER12DLOCALFINITEELEMENT_HH
