// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_HIERARCHICAL_P2_LOCALFINITEELEMENT_HH
#define DUNE_HIERARCHICAL_P2_LOCALFINITEELEMENT_HH

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/lagrange/lagrangesimplex.hh>

#include "hierarchicalp2/hierarchicalsimplexp2localbasis.hh"
#include "hierarchicalp2/hierarchicalsimplexp2localinterpolation.hh"

namespace Dune
{

  /** \todo Please doc me !
   */
  template<class D, class R, int dim>
  class HierarchicalP2LocalFiniteElement
  {

    static_assert(1 <= dim && dim <= 3,
                  "HierarchicalP2LocalFiniteElement only implemented for dim==1, 2, 3.");

  public:
    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<
        HierarchicalSimplexP2LocalBasis<D,R,dim>,
        typename LagrangeSimplexLocalFiniteElement<D,R,dim,2>::Traits::LocalCoefficientsType,
        HierarchicalSimplexP2LocalInterpolation<HierarchicalSimplexP2LocalBasis<D,R,dim> > > Traits;

    /** \todo Please doc me !
     */
    HierarchicalP2LocalFiniteElement ()
    {}

    /** \todo Please doc me !
     */
    const typename Traits::LocalBasisType& localBasis () const
    {
      return basis;
    }

    /** \todo Please doc me !
     */
    const typename Traits::LocalCoefficientsType& localCoefficients () const
    {
      return coefficients;
    }

    /** \todo Please doc me !
     */
    const typename Traits::LocalInterpolationType& localInterpolation () const
    {
      return interpolation;
    }

    /** \brief Number of shape functions in this finite element */
    unsigned int size () const
    {
      return basis.size();
    }

    /** \todo Please doc me !
     */
    static constexpr GeometryType type ()
    {
      return GeometryTypes::simplex(dim);
    }

  private:
    HierarchicalSimplexP2LocalBasis<D,R,dim> basis;

    typename Traits::LocalCoefficientsType coefficients;

    HierarchicalSimplexP2LocalInterpolation<HierarchicalSimplexP2LocalBasis<D,R,dim> > interpolation;
  };

}

#endif
