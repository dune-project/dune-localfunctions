// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_REFINED_REFINEDP1_HH
#define DUNE_LOCALFUNCTIONS_REFINED_REFINEDP1_HH

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/lagrange/p0.hh>

#include <dune/localfunctions/lagrange/lagrangesimplex.hh>
#include <dune/localfunctions/refined/refinedp1/refinedp1localbasis.hh>

namespace Dune
{

  /** \brief Piecewise linear continuous Lagrange functions on a uniformly refined simplex element
   *
   * \tparam D Number type used for domain coordinates
   * \tparam R Number type used for shape function values
   * \tparam dim Dimension of the domain
   */
  template<class D, class R, int dim>
  class RefinedP1LocalFiniteElement
  {
  public:
    /** \brief Export all types used by this implementation
     */
    typedef LocalFiniteElementTraits<RefinedP1LocalBasis<D,R,dim>,
                                     Impl::LagrangeSimplexLocalCoefficients<dim,2>,
                                     Impl::LagrangeSimplexLocalInterpolation<Impl::LagrangeSimplexLocalBasis<D,R,dim,2> > > Traits;

    /** \brief Default constructor
     */
    RefinedP1LocalFiniteElement ()
    {}

    /** \brief The set of shape functions
     */
    const typename Traits::LocalBasisType& localBasis () const
    {
      return basis_;
    }

    /** \brief Produces the assignments of the degrees of freedom to the element subentities
     */
    const typename Traits::LocalCoefficientsType& localCoefficients () const
    {
      return coefficients_;
    }

    /** \brief Evaluates all degrees of freedom for a given function
     */
    const typename Traits::LocalInterpolationType& localInterpolation () const
    {
      return interpolation_;
    }

    /** \brief Number of shape functions of this finite element */
    unsigned int size () const
    {
      return basis_.size();
    }

    /** \brief The element type that this finite element is defined on
     */
    static constexpr GeometryType type ()
    {
      return GeometryTypes::simplex(dim);
    }

  private:
    RefinedP1LocalBasis<D,R,dim> basis_;
    Impl::LagrangeSimplexLocalCoefficients<dim,2> coefficients_;
    // Yes, the template argument here really is LagrangeSimplexLocalBasis, even though this is not
    // the local basis of the refined locale finite element:  The reason is that LagrangeSimplexLocalInterpolation
    // uses this argument to determine the polynomial order, and RefinedP1LocalBasis returns order 1
    // whereas order 2 is needed here.
    Impl::LagrangeSimplexLocalInterpolation<Impl::LagrangeSimplexLocalBasis<D,R,dim,2> > interpolation_;
  };

}

#endif   // DUNE_LOCALFUNCTIONS_REFINED_REFINEDP1_HH
