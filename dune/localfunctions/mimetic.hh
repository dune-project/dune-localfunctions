// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_MIMETICLOCALFINITEELEMENT_HH
#define DUNE_MIMETICLOCALFINITEELEMENT_HH

#include <dune/geometry/type.hh>

#include "common/localfiniteelementtraits.hh"
#include "mimetic/mimeticall.hh"

namespace Dune
{
  /**
   * \ingroup Mimetic
   *
   * \tparam D Type to represent the field in the domain.
   * \tparam R Type to represent the field in the range.
   * \tparam dim Domain dimension
   */
  template<class D, class R, int dim>
  class MimeticLocalFiniteElement
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
