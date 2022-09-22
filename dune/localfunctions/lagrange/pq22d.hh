// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_PQ22DLOCALFINITEELEMENT_HH
#define DUNE_PQ22DLOCALFINITEELEMENT_HH

#include <dune/common/fmatrix.hh>

#include <dune/localfunctions/common/localfiniteelementvariant.hh>

#include <dune/localfunctions/lagrange/lagrangesimplex.hh>
#include <dune/localfunctions/lagrange/lagrangecube.hh>

namespace Dune
{
  template<class D, class R>
  class PQ22DLocalFiniteElement
  {
    using LFEVariant = LocalFiniteElementVariant<LagrangeSimplexLocalFiniteElement<D,R,2,2>,
                                                 LagrangeCubeLocalFiniteElement<D,R,2,2> >;
  public:
    using Traits = typename LFEVariant::Traits;

    PQ22DLocalFiniteElement ( const GeometryType &gt )
    {
      if ( gt.isTriangle() )
        lfeVariant_ = LagrangeSimplexLocalFiniteElement<D,R,2,2>();
      else if ( gt.isQuadrilateral() )
        lfeVariant_ = LagrangeCubeLocalFiniteElement<D,R,2,2>();
    }

    PQ22DLocalFiniteElement ( const GeometryType &gt, const std::vector<unsigned int> vertexmap )
    {
      if ( gt.isTriangle() )
        lfeVariant_ = LagrangeSimplexLocalFiniteElement<D,R,2,2>(vertexmap);
      else if ( gt.isQuadrilateral() )
        lfeVariant_ = LagrangeCubeLocalFiniteElement<D,R,2,2>();
    }

    const typename Traits::LocalBasisType& localBasis () const
    {
      return lfeVariant_.localBasis();
    }

    const typename Traits::LocalCoefficientsType& localCoefficients () const
    {
      return lfeVariant_.localCoefficients();
    }

    const typename Traits::LocalInterpolationType& localInterpolation () const
    {
      return lfeVariant_.localInterpolation();
    }

    /** \brief Number of shape functions in this finite element */
    unsigned int size () const
    {
      return lfeVariant_.size();
    }

    GeometryType type () const
    {
      return lfeVariant_.type();
    }

  private:

    LFEVariant lfeVariant_;
  };

}

#endif
