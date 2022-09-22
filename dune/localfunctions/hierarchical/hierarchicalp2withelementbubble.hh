// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_HIERARCHICAL_P2_WITH_ELEMENTBUBBLE_LOCALFINITEELEMENT_HH
#define DUNE_HIERARCHICAL_P2_WITH_ELEMENTBUBBLE_LOCALFINITEELEMENT_HH

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/hierarchical/hierarchicalp2withelementbubble/hierarchicalsimplexp2withelementbubble.hh>


namespace Dune
{

  /** \todo Please doc me !
   */
  template<class D, class R, int dim>
  class HierarchicalP2WithElementBubbleLocalFiniteElement
  {

    static_assert(dim==2, "HierarchicalP2WithElementBubbleLocalFiniteElement only implemented for dim==2.");

  public:
    /** \todo Please doc me !
     */
    typedef LocalFiniteElementTraits<HierarchicalSimplexP2WithElementBubbleLocalBasis<D,R,dim>,
        HierarchicalSimplexP2WithElementBubbleLocalCoefficients<dim>,
        HierarchicalSimplexP2WithElementBubbleLocalInterpolation<HierarchicalSimplexP2WithElementBubbleLocalBasis<D,R,dim> > > Traits;

    /** \todo Please doc me !
     */
    HierarchicalP2WithElementBubbleLocalFiniteElement ()
    {}

    /** \todo Please doc me !
     */
    const typename Traits::LocalBasisType& localBasis () const
    {
      return basis_;
    }

    /** \todo Please doc me !
     */
    const typename Traits::LocalCoefficientsType& localCoefficients () const
    {
      return coefficients_;
    }

    /** \todo Please doc me !
     */
    const typename Traits::LocalInterpolationType& localInterpolation () const
    {
      return interpolation_;
    }

    /** \brief Number of shape functions in this finite element */
    unsigned int size () const
    {
      return basis_.size();
    }

    /** \todo Please doc me !
     */
    static constexpr GeometryType type ()
    {
      return GeometryTypes::triangle;
    }

  private:
    HierarchicalSimplexP2WithElementBubbleLocalBasis<D,R,dim> basis_;

    HierarchicalSimplexP2WithElementBubbleLocalCoefficients<dim> coefficients_;

    HierarchicalSimplexP2WithElementBubbleLocalInterpolation<HierarchicalSimplexP2WithElementBubbleLocalBasis<D,R,dim> > interpolation_;
  };

}

#endif
