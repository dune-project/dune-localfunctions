// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_ENRICHED_SIMPLEXP1BUBBLE_HH
#define DUNE_LOCALFUNCTIONS_ENRICHED_SIMPLEXP1BUBBLE_HH

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/enriched/simplexp1bubble/localbasis.hh>
#include <dune/localfunctions/enriched/simplexp1bubble/localcoefficients.hh>
#include <dune/localfunctions/enriched/simplexp1bubble/localinterpolation.hh>


namespace Dune
{
  /**
   * \brief Linear Lagrange functions enriched with an element bubble function.
   *
   * The set of basis functions contains the classical Lagrange basis functions
   * of order 1, i.e., the barycentric coordinates, and a single element "bubble"
   * function that vanishes on all faces of the element. The bubble function is
   * simply defined as the product of all linear basis functions and thus has
   * polynomial order `dim+1`.
   *
   * A classical example where this kind of basis is used in the discretization
   * of the Stokes equation with the stable mixed-element called MINI element,
   * see
   *
   *   Arnold, D.N., Brezzi, F. and Fortin, M. A stable finite element for the
   *   Stokes equations. Calcolo 21, 337-344 (1984). doi: 10.1007/BF02576171
   *
   * The velocity field is discretized with continuous piecewise linear
   * functions enriched by a bubble function.
   *
   * \note The implementation here is restricted to simplex elements.
   *
   * \tparam D    Type to represent the field in the domain.
   * \tparam R    Type to represent the field in the range.
   * \tparam dim  Dimension of the domain.
   **/
  template<class D, class R, int dim>
  class SimplexP1BubbleLocalFiniteElement
  {
  public:
    //! Type of the local basis
    using LocalBasisType = SimplexP1BubbleLocalBasis<D,R,dim>;

    //! Type of the local coefficients
    using LocalCoefficientsType = SimplexP1BubbleLocalCoefficients<dim>;

    //! Type of the local interpolation
    using LocalInterpolationType = SimplexP1BubbleLocalInterpolation<LocalBasisType>;

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

} // end namespace Dune

#endif // DUNE_LOCALFUNCTIONS_ENRICHED_SIMPLEXP1BUBBLE_HH
