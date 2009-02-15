// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_P0LOCALFINITEELEMENT_HH
#define DUNE_P0LOCALFINITEELEMENT_HH

#include <dune/common/geometrytype.hh>

#include "common/localfiniteelement.hh"
#include "p0/p0localbasis.hh"
#include "p0/p0localcoefficients.hh"
#include "p0/p0localinterpolation.hh"

namespace Dune
{

  template<class D, class R, int d>
  class P0LocalFiniteElement : LocalFiniteElementInterface<
                                   LocalFiniteElementTraits<P0LocalBasis<D,R,d>,P0LocalCoefficients,
                                       P0LocalInterpolation<P0LocalBasis<D,R,d> > >,
                                   P0LocalFiniteElement<D,R,d> >
  {
  public:
    typedef LocalFiniteElementTraits<P0LocalBasis<D,R,d>, P0LocalCoefficients,
        P0LocalInterpolation<P0LocalBasis<D,R,d> > > Traits;

    P0LocalFiniteElement (GeometryType::BasicType basicType)
      : interpolation(basicType,d), gt(basicType,d)
    {}

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
    P0LocalBasis<D,R,d> basis;
    P0LocalCoefficients coefficients;
    P0LocalInterpolation<P0LocalBasis<D,R,d> > interpolation;
    GeometryType gt;
  };

}

#endif
