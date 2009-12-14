// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MIMETICLOCALFINITEELEMENT_HH
#define DUNE_MIMETICLOCALFINITEELEMENT_HH

#include <dune/common/geometrytype.hh>

#include "common/localfiniteelement.hh"
#include "mimetic/mimeticall.hh"

namespace Dune
{
  template<class D, class R, int dim>
  class MimeticLocalFiniteElement
    : public Dune::LocalFiniteElementInterface<
          Dune::LocalFiniteElementTraits<MimeticLocalBasis<D,R,dim>,
              MimeticLocalCoefficients,
              MimeticLocalInterpolation<MimeticLocalBasis<D,R,dim> > >,
          MimeticLocalFiniteElement<D,R,dim> >
  {
    Dune::GeometryType gt;
    MimeticLocalBasis<D,R,dim> basis;
    MimeticLocalCoefficients coefficients;
    MimeticLocalInterpolation<MimeticLocalBasis<D,R,dim> > interpolation;

  public:
    typedef Dune::LocalFiniteElementTraits<MimeticLocalBasis<D,R,dim>,
        MimeticLocalCoefficients,
        MimeticLocalInterpolation<MimeticLocalBasis<D,R,dim> > > Traits;

    MimeticLocalFiniteElement ()
    {}

    MimeticLocalFiniteElement (Dune::GeometryType::BasicType basicType)
      : gt(basicType,dim)
    {}

    MimeticLocalFiniteElement (Dune::GeometryType::BasicType basicType, unsigned int variant)
      : gt(basicType,dim), basis(variant), coefficients(variant)
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

    Dune::GeometryType type () const { return gt; }
  };
}

#endif
