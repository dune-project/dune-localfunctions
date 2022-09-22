// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_CROUZEIXRAVIART_HH
#define DUNE_LOCALFUNCTIONS_CROUZEIXRAVIART_HH

#include <array>
#include <numeric>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localinterpolation.hh>
#include <dune/localfunctions/common/localkey.hh>

namespace Dune { namespace Impl
{
   /** \brief Crouzeix-Raviart element on the reference simplex

     \tparam D Type to represent the field in the domain
     \tparam R Type to represent the field in the range
     \tparam dim Dimension of the domain simplex
   */
  template<class D, class R, unsigned int dim>
  class CrouzeixRaviartLocalBasis
  {
  public:
    using Traits = LocalBasisTraits<D,dim,FieldVector<D,dim>,R,1,FieldVector<R,1>,FieldMatrix<R,1,dim> >;

    /** \brief Number of shape functions -- one for each simplex facet
     */
    static constexpr unsigned int size ()
    {
      return dim+1;
    }

    //! \brief Evaluate all shape functions
    void evaluateFunction(const typename Traits::DomainType& x,
                          std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(size());

      std::fill(out.begin(), out.end()-1, 1.0);
      out.back() = 1.0-dim;

      for (unsigned int i=0; i<dim; i++)
      {
        out[i] -= dim * x[dim-i-1];
        out.back() += dim*x[i];
      }
    }

    /** \brief Evaluate Jacobians of all shape functions
     *
     * \param x Point in the reference simplex where to evaluation the Jacobians
     * \param[out] out The Jacobians of all shape functions at the point x
     */
    void evaluateJacobian(const typename Traits::DomainType& x,
                          std::vector<typename Traits::JacobianType>& out) const
    {
      out.resize(size());

      for (unsigned i=0; i<dim; i++)
        for (unsigned j=0; j<dim; j++)
          out[i][0][j] = (i==(dim-1-j)) ? -(double)dim : 0;

      std::fill(out.back()[0].begin(), out.back()[0].end(), dim);
    }

    /** \brief Evaluate partial derivatives of any order of all shape functions
     *
     * \param order Order of the partial derivatives, in the classic multi-index notation
     * \param in Position where to evaluate the derivatives
     * \param[out] out The desired partial derivatives
     */
    void partial(const std::array<unsigned int,dim>& order,
                 const typename Traits::DomainType& in,
                 std::vector<typename Traits::RangeType>& out) const
    {
      auto totalOrder = std::accumulate(order.begin(), order.end(), 0);

      out.resize(size());

      if (totalOrder == 0) {
        evaluateFunction(in, out);
        return;
      }

      if (totalOrder==1)
      {
        auto direction = std::find(order.begin(), order.end(), 1)-order.begin();

        for (unsigned int i=0; i<dim; i++)
          out[i] = (i==(dim-1-direction)) ? -(double)dim : 0.0;

        out.back()[0] = dim;
      }
      else  // all higher order derivatives are zero
        std::fill(out.begin(), out.end(), 0);
    }

    //! \brief Polynomial order of the shape functions
    static constexpr unsigned int order ()
    {
      return 1;
    }
  };

  /** \brief Associations of the Crouzeix-Raviart degrees of freedom to the facets of the reference simplex
   *
   * \tparam dim Dimension of the reference simplex
   */
  template<unsigned int dim>
  class CrouzeixRaviartLocalCoefficients
  {
  public:
    //! \brief Default constructor
    CrouzeixRaviartLocalCoefficients ()
    : localKeys_(size())
    {
      for (std::size_t i=0; i<size(); i++)
        localKeys_[i] = LocalKey(i,dim-1,0);
    }

    //! Number of coefficients
    static constexpr std::size_t size ()
    {
      return dim+1;
    }

    //! Get i-th local key
    const LocalKey& localKey (std::size_t i) const
    {
      return localKeys_[i];
    }

  private:
    std::vector<LocalKey> localKeys_;
  };

  /** \brief Evaluate the degrees of freedom of a Crouzeix-Raviart basis
   *
   * \tparam LocalBasis The corresponding set of shape functions
   */
  template<class LocalBasis>
  class CrouzeixRaviartLocalInterpolation
  {
  public:

    /** \brief Evaluate a given function at the facet midpoints
     *
     * \tparam F Type of function to evaluate
     * \tparam C Type used for the values of the function
     * \param[in] ff Function to evaluate
     * \param[out] out Array of function values
     */
    template<typename F, typename C>
    void interpolate (const F& ff, std::vector<C>& out) const
    {
      constexpr auto dim = LocalBasis::Traits::dimDomain;

      auto&& f = Impl::makeFunctionWithCallOperator<typename LocalBasis::Traits::DomainType>(ff);

      out.resize(LocalBasis::size());

      // Evaluate at the facet midpoints
      auto refSimplex = ReferenceElements<typename LocalBasis::Traits::DomainFieldType,dim>::simplex();
      for (int i=0; i<refSimplex.size(1); i++)
        out[i] = f(refSimplex.template geometry<1>(i).center());
    }
  };

} }   // namespace Dune::Impl

namespace Dune
{
  /** \brief Crouzeix-Raviart finite element
   *
   * \tparam D type used for domain coordinates
   * \tparam R type used for function values
   * \tparam dim dimension of the reference element
   */
  template<class D, class R, int dim>
  class CrouzeixRaviartLocalFiniteElement
  {
  public:
    /** \brief Export number types, dimensions, etc.
     */
    using Traits = LocalFiniteElementTraits<Impl::CrouzeixRaviartLocalBasis<D,R,dim>,
                                            Impl::CrouzeixRaviartLocalCoefficients<dim>,
                                            Impl::CrouzeixRaviartLocalInterpolation<Impl::CrouzeixRaviartLocalBasis<D,R,dim> > >;

    /** \brief Returns the local basis, i.e., the set of shape functions
     */
    const typename Traits::LocalBasisType& localBasis() const
    {
      return basis_;
    }

    /** \brief Returns the assignment of the degrees of freedom to the element subentities
     */
    const typename Traits::LocalCoefficientsType& localCoefficients() const
    {
      return coefficients_;
    }

    /** \brief Returns object that evaluates degrees of freedom
     */
    const typename Traits::LocalInterpolationType& localInterpolation() const
    {
      return interpolation_;
    }

    /** \brief The number of shape functions */
    static constexpr std::size_t size()
    {
      return dim+1;
    }

    /** \brief The reference element that the local finite element is defined on
     */
    static constexpr GeometryType type()
    {
      return GeometryTypes::simplex(dim);
    }

  private:
    Impl::CrouzeixRaviartLocalBasis<D,R,dim> basis_;
    Impl::CrouzeixRaviartLocalCoefficients<dim> coefficients_;
    Impl::CrouzeixRaviartLocalInterpolation<Impl::CrouzeixRaviartLocalBasis<D,R,dim> > interpolation_;
  };

}        // namespace Dune

#endif   // DUNE_LOCALFUNCTIONS_CROUZEIXRAVIART_HH
