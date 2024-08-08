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
  /**
   * \brief Linear Lagrange functions enriched with quadratic edge bubble functions
   *        and an element bubble function.
   *
   * The set of basis functions contains the classical Lagrange basis functions
   * of order 1, i.e., the barycentric coordinates, edge bubble functions as the
   * product of two linear functions, and a single element bubble function that
   * vanishes on all faces of the element. The element bubble function is the
   * product of all linear basis functions and thus has polynomial order `dim+1`.
   *
   * \note The implementation here is restricted to simplex elements.
   *
   * \tparam D    Type to represent the field in the domain.
   * \tparam R    Type to represent the field in the range.
   * \tparam dim  Dimension of the domain.
   **/
  template<class D, class R, int dim>
  class HierarchicalP2WithElementBubbleLocalFiniteElement
  {
  public:
    //! Type of the local basis
    using LocalBasisType = HierarchicalSimplexP2WithElementBubbleLocalBasis<D,R,dim>;

    //! Type of the local coefficients
    using LocalCoefficientsType = HierarchicalSimplexP2WithElementBubbleLocalCoefficients<dim>;

    //! Type of the local interpolation
    using LocalInterpolationType = HierarchicalSimplexP2WithElementBubbleLocalInterpolation<LocalBasisType,dim>;

    //! Traits type that specifies the local basis, coefficients, and interpolation type.
    using Traits = LocalFiniteElementTraits<LocalBasisType,LocalCoefficientsType,LocalInterpolationType>;


    //! Returns the local basis, i.e., the set of shape functions
    const LocalBasisType& localBasis () const
    {
      return basis_;
    }

    //! Returns the assignment of the degrees of freedom to the element subentities
    const LocalCoefficientsType& localCoefficients () const
    {
      return coefficients_;
    }

    //! Returns object that evaluates degrees of freedom
    const LocalInterpolationType& localInterpolation () const
    {
      return interpolation_;
    }

    //! Returns the number of shape functions in this finite-element
    static constexpr std::size_t size () noexcept
    {
      return LocalBasisType::size();
    }

    //! Returns the type of the geometry the finite-element is attached to
    static constexpr GeometryType type () noexcept
    {
      return GeometryTypes::simplex(dim);
    }

  private:
    LocalCoefficientsType coefficients_{};
    [[no_unique_address]] LocalBasisType basis_{};
    [[no_unique_address]] LocalInterpolationType interpolation_{};
  };

}

#endif
