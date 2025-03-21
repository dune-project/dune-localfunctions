// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_LAGRANGE_LAGRANGESIMPLEX_HH
#define DUNE_LOCALFUNCTIONS_LAGRANGE_LAGRANGESIMPLEX_HH

#include <array>
#include <numeric>
#include <algorithm>

#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/math.hh>
#include <dune/common/rangeutilities.hh>

#include <dune/geometry/referenceelements.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localkey.hh>

namespace Dune { namespace Impl
{
   /** \brief Lagrange shape functions of arbitrary order on the reference simplex

     Lagrange shape functions of arbitrary order have the property that
     \f$\hat\phi^i(x^{(j)}) = \delta_{i,j}\f$ for certain points \f$x^{(j)}\f$.

     \tparam D Type to represent the field in the domain
     \tparam R Type to represent the field in the range
     \tparam dim Dimension of the domain simplex
     \tparam k Polynomial order
   */
  template<class D, class R, unsigned int dim, unsigned int k>
  class LagrangeSimplexLocalBasis
  {

    // Compute the rescaled barycentric coordinates of x.
    // We rescale the simplex by k and then compute the
    // barycentric coordinates with respect to the points
    // p_i = e_i (for i=0,...,dim-1) and p_dim=0.
    // Notice that then the Lagrange points have the barycentric
    // coordinates (i_0,...,i_d) where i_j are all non-negative
    // integers satisfying the constraint sum i_j = k.
    constexpr auto barycentric(const auto& x) const
    {
      auto b = std::array<R,dim+1>{};
      b[dim] = k;
      for(auto i : Dune::range(dim))
      {
        b[i] = k*x[i];
        b[dim] -= b[i];
      }
      return b;
    }

    // Evaluate the univariate Lagrange polynomials L_i(t) for i=0,...,k where
    //
    // L_i(t) = (t-0)/(i-0) * ... * (t-(i-1))/(i-(i-1))
    //        = (t-0)*...*(t-(i-1))/(i!);
    static constexpr void evaluateLagrangePolynomials(const R& t, auto& L)
    {
      L[0] = 1;
      for (auto i : Dune::range(k))
        L[i+1] = L[i] * (t - i) / (i+1);
    }

    // Evaluate the univariate Lagrange polynomial derivatives L_i(t) for i=0,...,k
    // up to given maxDerivativeOrder.
    static constexpr void evaluateLagrangePolynomialDerivative(const R& t, auto& LL, unsigned int maxDerivativeOrder)
    {
      auto& L = LL[0];
      L[0] = 1;
      for (auto i : Dune::range(k))
        L[i+1] = L[i] * (t - i) / (i+1);
      for(auto j : Dune::range(maxDerivativeOrder))
      {
        auto& F = LL[j];
        auto& DF = LL[j+1];
        DF[0] = 0;
        for (auto i : Dune::range(k))
          DF[i+1] = (DF[i] * (t - i) + (j+1)*F[i]) / (i+1);
      }
    }

    using BarycentricMultiIndex = std::array<unsigned int,dim+1>;


    // This computed the required partial derivative given by the multi-index
    // beta of a product of a function given as a product of dim+1 derivatives
    // of univariate Lagrange polynomials of the dim+1 barycentric coordinates.
    // The polynomials in the product are specified as follows:
    //
    // The table L contains all required derivatives of all univariate
    // polynomials evaluated at all barycentric coordinates.  The two
    // multi-indices i and alpha encode that the polynomial for the
    // j-th barycentric coordinate is the alpha-j-th derivative of
    // the i_j-the Lagrange polynomial.
    //
    // Hence this method computes D_beta f(x) where f(x) is the product
    // \f$f(x) = \prod_{j=0}^{d} L_{i_j}^{(alpha_j)}(x_j) \f$.
    static constexpr R barycentricDerivative(
        BarycentricMultiIndex beta,
        const auto&L,
        const BarycentricMultiIndex& i,
        const BarycentricMultiIndex& alpha = {})
    {
      // If there are unprocessed derivatives left we search the first unprocessed
      // partial derivative direction j and compute it using the product and chain rule.
      // The remaining partial derivatives are forwarded to the recursion.
      // Notice that the partial derivative of the last barycentric component
      // wrt to the j-th component is -1 which is responsible for the second
      // term in the sum. Furthermore we get the factor k due to the scaling of
      // the simplex.
      for(auto j : Dune::range(dim))
      {
        if (beta[j] > 0)
        {
          auto leftDerivatives = alpha;
          auto rightDerivatives = alpha;
          leftDerivatives[j]++;
          rightDerivatives.back()++;
          beta[j]--;
          return (barycentricDerivative(beta, L, i, leftDerivatives) - barycentricDerivative(beta, L, i, rightDerivatives))*k;
        }
      }
      // If there are no unprocessed derivatives we can simply evaluate
      // the product of the derivatives of the Lagrange polynomials with
      // given indices i_j and derivative orders alpha_j
      // Evaluate the product of the univariate Lagrange polynomials
      // with given indices and orders.
      auto y = R(1);
      for(auto j : Dune::range(dim+1))
        y *= L[j][alpha[j]][i[j]];
      return y;
    }


  public:
    using Traits = LocalBasisTraits<D,dim,FieldVector<D,dim>,R,1,FieldVector<R,1>,FieldMatrix<R,1,dim> >;

    /** \brief Number of shape functions
     *
     * See https://en.wikipedia.org/wiki/Figurate_number for an explanation of the formula
     */
    static constexpr unsigned int size ()
    {
      return binomial(k+dim,dim);
    }

    //! \brief Evaluate all shape functions
    void evaluateFunction(const typename Traits::DomainType& x,
                          std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(size());

      // Specialization for zero-order case
      if (k==0)
      {
        out[0] = 1;
        return;
      }

      // Specialization for first-order case
      if (k==1)
      {
        out[0] = 1.0;
        for (size_t i=0; i<dim; i++)
        {
          out[0]  -= x[i];
          out[i+1] = x[i];
        }
        return;
      }

      // Compute rescaled barycentric coordinates of x
      auto z = barycentric(x);

      auto L = std::array<std::array<R,k+1>, dim+1>();
      for (auto j : Dune::range(dim+1))
        evaluateLagrangePolynomials(z[j], L[j]);

      if (dim==1)
      {
        unsigned int n = 0;
        for (auto i0 : Dune::range(k + 1))
          for (auto i1 : std::array{k - i0})
            out[n++] = L[0][i0] * L[1][i1];
        return;
      }
      if (dim==2)
      {
        unsigned int n=0;
        for (auto i1 : Dune::range(k + 1))
          for (auto i0 : Dune::range(k - i1 + 1))
            for (auto i2 : std::array{k - i1 - i0})
              out[n++] = L[0][i0] * L[1][i1] * L[2][i2];
        return;
      }
      if (dim==3)
      {
        unsigned int n = 0;
        for (auto i2 : Dune::range(k + 1))
          for (auto i1 : Dune::range(k - i2 + 1))
            for (auto i0 : Dune::range(k - i2 - i1 + 1))
              for (auto i3 : std::array{k - i2 - i1 - i0})
                out[n++] = L[0][i0] * L[1][i1]  * L[2][i2] * L[3][i3];
        return;
      }

      DUNE_THROW(NotImplemented, "LagrangeSimplexLocalBasis for k>=2 only implemented for dim<=3");
    }

    /** \brief Evaluate Jacobian of all shape functions
     *
     * \param x Point in the reference simplex where to evaluation the Jacobians
     * \param[out] out The Jacobians of all shape functions at the point x
     */
    void evaluateJacobian(const typename Traits::DomainType& x,
                          std::vector<typename Traits::JacobianType>& out) const
    {
      out.resize(size());

      // Specialization for k==0
      if (k==0)
      {
        std::fill(out[0][0].begin(), out[0][0].end(), 0);
        return;
      }

      // Specialization for k==1
      if (k==1)
      {
        std::fill(out[0][0].begin(), out[0][0].end(), -1);

        for (unsigned int i=0; i<dim; i++)
          for (unsigned int j=0; j<dim; j++)
            out[i+1][0][j] = (i==j);

        return;
      }

      // Compute rescaled barycentric coordinates of x
      auto z = barycentric(x);

      // L[j][m][i] is the m-th derivative of the i-th Lagrange polynomial at z[j]
      auto L = std::array<std::array<std::array<R,k+1>, 2>, dim+1>();
      for (auto j : Dune::range(dim+1))
        evaluateLagrangePolynomialDerivative(z[j], L[j], 1);

      if (dim==1)
      {
        unsigned int n = 0;
        for (auto i0 : Dune::range(k + 1))
        {
          for (auto i1 : std::array{k-i0})
          {
            out[n][0][0] = (L[0][1][i0] * L[1][0][i1] - L[0][0][i0] * L[1][1][i1])*k;
            ++n;
          }
        }
        return;
      }
      if (dim==2)
      {
        unsigned int n=0;
        for (auto i1 : Dune::range(k + 1))
        {
          for (auto i0 : Dune::range(k - i1 + 1))
          {
            for (auto i2 : std::array{k - i1 - i0})
            {
              out[n][0][0] = (L[0][1][i0] * L[1][0][i1] * L[2][0][i2] - L[0][0][i0] * L[1][0][i1] * L[2][1][i2])*k;
              out[n][0][1] = (L[0][0][i0] * L[1][1][i1] * L[2][0][i2] - L[0][0][i0] * L[1][0][i1] * L[2][1][i2])*k;
              ++n;
            }
          }
        }
        return;
      }
      if (dim==3)
      {
        unsigned int n = 0;
        for (auto i2 : Dune::range(k + 1))
        {
          for (auto i1 : Dune::range(k - i2 + 1))
          {
            for (auto i0 : Dune::range(k - i2 - i1 + 1))
            {
              for (auto i3 : std::array{k - i2 - i1 - i0})
              {
                out[n][0][0] = (L[0][1][i0] * L[1][0][i1] * L[2][0][i2] * L[3][0][i3] - L[0][0][i0] * L[1][0][i1] * L[2][0][i2] * L[3][1][i3])*k;
                out[n][0][1] = (L[0][0][i0] * L[1][1][i1] * L[2][0][i2] * L[3][0][i3] - L[0][0][i0] * L[1][0][i1] * L[2][0][i2] * L[3][1][i3])*k;
                out[n][0][2] = (L[0][0][i0] * L[1][0][i1] * L[2][1][i2] * L[3][0][i3] - L[0][0][i0] * L[1][0][i1] * L[2][0][i2] * L[3][1][i3])*k;
                ++n;
              }
            }
          }
        }

        return;
      }

      DUNE_THROW(NotImplemented, "LagrangeSimplexLocalBasis for k>=2 only implemented for dim<=3");
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
      auto totalOrder = std::accumulate(order.begin(), order.end(), 0u);

      out.resize(size());

      // Derivative order zero corresponds to the function evaluation.
      if (totalOrder == 0)
      {
        evaluateFunction(in, out);
        return;
      }

      // Derivatives of order >k are all zero.
      if (totalOrder > k)
      {
        for(auto& out_i : out)
          out_i = 0;
        return;
      }

      // It remains to cover the cases 0 < totalOrder<= k.

      if (k==1)
      {
        if (totalOrder==1)
        {
          auto direction = std::find(order.begin(), order.end(), 1);
          out[0] = -1;
          for (unsigned int i=0; i<dim; i++)
            out[i+1] = (i==(direction-order.begin()));
        }
        return;
      }

      // Since the required stack storage depends on the dynamic total order,
      // we need to do a dynamic to static dispatch by enumerating all supported
      // static orders.
      auto supportedStaticOrders = Dune::range(Dune::index_constant<1>{}, Dune::index_constant<k+1>{});
      return Dune::Hybrid::switchCases(supportedStaticOrders, totalOrder, [&](auto staticTotalOrder) {

        // Compute rescaled barycentric coordinates of x
        auto z = barycentric(in);

        // L[j][m][i] is the m-th derivative of the i-th Lagrange polynomial at z[j]
        auto L = std::array<std::array<std::array<R, k+1>, staticTotalOrder+1>, dim+1>();
        for (auto j : Dune::range(dim))
          evaluateLagrangePolynomialDerivative(z[j], L[j], order[j]);
        evaluateLagrangePolynomialDerivative(z[dim], L[dim], totalOrder);

        auto barycentricOrder = BarycentricMultiIndex{};
        for (auto j : Dune::range(dim))
          barycentricOrder[j] = order[j];
        barycentricOrder[dim] = 0;

        if constexpr (dim==1)
        {
          unsigned int n = 0;
          for (auto i0 : Dune::range(k + 1))
            for (auto i1 : std::array{k - i0})
              out[n++] = barycentricDerivative(barycentricOrder, L, BarycentricMultiIndex{i0, i1});
        }
        if constexpr (dim==2)
        {
          unsigned int n=0;
          for (auto i1 : Dune::range(k + 1))
            for (auto i0 : Dune::range(k - i1 + 1))
              for (auto i2 : std::array{k - i1 - i0})
                out[n++] = barycentricDerivative(barycentricOrder, L, BarycentricMultiIndex{i0, i1, i2});
        }
        if constexpr (dim==3)
        {
          unsigned int n = 0;
          for (auto i2 : Dune::range(k + 1))
            for (auto i1 : Dune::range(k - i2 + 1))
              for (auto i0 : Dune::range(k - i2 - i1 + 1))
                for (auto i3 : std::array{k - i2 - i1 - i0})
                  out[n++] = barycentricDerivative(barycentricOrder, L, BarycentricMultiIndex{i0, i1, i2, i3});
        }
      });
    }

    //! \brief Polynomial order of the shape functions
    static constexpr unsigned int order ()
    {
      return k;
    }
  };

  /** \brief Associations of the Lagrange degrees of freedom to subentities of the reference simplex
   *
   * \tparam dim Dimension of the reference simplex
   * \tparam k Polynomial order of the Lagrange space
   */
  template<unsigned int dim, unsigned int k>
  class LagrangeSimplexLocalCoefficients
  {
  public:
    //! \brief Default constructor
    LagrangeSimplexLocalCoefficients ()
    : localKeys_(size())
    {
      if (k==0)
      {
        localKeys_[0] = LocalKey(0,0,0);
        return;
      }

      if (k==1)
      {
        for (std::size_t i=0; i<size(); i++)
          localKeys_[i] = LocalKey(i,dim,0);
        return;
      }

      if (dim==1)
      {
        // Order is at least 2 here
        localKeys_[0] = LocalKey(0,1,0);          // vertex dof
        for (unsigned int i=1; i<k; i++)
          localKeys_[i] = LocalKey(0,0,i-1);      // element dofs
        localKeys_.back() = LocalKey(1,1,0);      // vertex dof
        return;
      }

      if (dim==2)
      {
        int n=0;
        int c=0;
        for (unsigned int j=0; j<=k; j++)
          for (unsigned int i=0; i<=k-j; i++)
          {
            if (i==0 && j==0)
            {
              localKeys_[n++] = LocalKey(0,2,0);
              continue;
            }
            if (i==k && j==0)
            {
              localKeys_[n++] = LocalKey(1,2,0);
              continue;
            }
            if (i==0 && j==k)
            {
              localKeys_[n++] = LocalKey(2,2,0);
              continue;
            }
            if (j==0)
            {
              localKeys_[n++] = LocalKey(0,1,i-1);
              continue;
            }
            if (i==0)
            {
              localKeys_[n++] = LocalKey(1,1,j-1);
              continue;
            }
            if (i+j==k)
            {
              localKeys_[n++] = LocalKey(2,1,j-1);
              continue;
            }
            localKeys_[n++] = LocalKey(0,0,c++);
          }
        return;
      }

      if (dim==3)
      {
        std::array<unsigned int, dim+1> vertexMap;
        for (unsigned int i=0; i<=dim; i++)
          vertexMap[i] = i;
        generateLocalKeys(vertexMap);
        return;
      }
      DUNE_THROW(NotImplemented, "LagrangeSimplexLocalCoefficients only implemented for k<=1 or dim<=3!");
    }

    /** Constructor for variants with permuted vertices
     *
     * \param vertexmap The permutation of the vertices.  This
     *   can for instance be generated from the global indices of
     *   the vertices by reducing those to the integers 0...dim
     */
    LagrangeSimplexLocalCoefficients (const std::array<unsigned int, dim+1> vertexMap)
    : localKeys_(size())
    {
      if (dim!=2 && dim!=3)
        DUNE_THROW(NotImplemented, "LagrangeSimplexLocalCoefficients only implemented for dim==2 and dim==3!");

      generateLocalKeys(vertexMap);
    }


    template<class VertexMap>
    LagrangeSimplexLocalCoefficients(const VertexMap &vertexmap)
    : localKeys_(size())
    {
      if (dim!=2 && dim!=3)
        DUNE_THROW(NotImplemented, "LagrangeSimplexLocalCoefficients only implemented for dim==2 and dim==3!");

      std::array<unsigned int, dim+1> vertexmap_array;
      std::copy(vertexmap, vertexmap + dim + 1, vertexmap_array.begin());
      generateLocalKeys(vertexmap_array);
    }

    //! number of coefficients
    static constexpr std::size_t size ()
    {
      return binomial(k+dim,dim);
    }

    //! get i'th index
    const LocalKey& localKey (std::size_t i) const
    {
      return localKeys_[i];
    }

  private:
    std::vector<LocalKey> localKeys_;

    void generateLocalKeys(const std::array<unsigned int, dim+1> vertexMap)
    {
      if (k==0)
      {
        localKeys_[0] = LocalKey(0,0,0);
        return;
      }

      if (dim==2)
      {
        // Create default assignment
        int n=0;
        int c=0;
        for (unsigned int j=0; j<=k; j++)
          for (unsigned int i=0; i<=k-j; i++)
          {
            if (i==0 && j==0)
            {
              localKeys_[n++] = LocalKey(0,2,0);
              continue;
            }
            if (i==k && j==0)
            {
              localKeys_[n++] = LocalKey(1,2,0);
              continue;
            }
            if (i==0 && j==k)
            {
              localKeys_[n++] = LocalKey(2,2,0);
              continue;
            }
            if (j==0)
            {
              localKeys_[n++] = LocalKey(0,1,i-1);
              continue;
            }
            if (i==0)
            {
              localKeys_[n++] = LocalKey(1,1,j-1);
              continue;
            }
            if (i+j==k)
            {
              localKeys_[n++] = LocalKey(2,1,j-1);
              continue;
            }
            localKeys_[n++] = LocalKey(0,0,c++);
          }

        // Flip edge orientations, if requested
        bool flip[3];
        flip[0] = vertexMap[0] > vertexMap[1];
        flip[1] = vertexMap[0] > vertexMap[2];
        flip[2] = vertexMap[1] > vertexMap[2];
        for (std::size_t i=0; i<size(); i++)
          if (localKeys_[i].codim()==1 && flip[localKeys_[i].subEntity()])
            localKeys_[i].index(k-2-localKeys_[i].index());

        return;
      }

      if (dim!=3)
        DUNE_THROW(NotImplemented, "LagrangeSimplexLocalCoefficients only implemented for dim==3!");

      unsigned int subindex[16];
      unsigned int codim_count[4] = {0};
      for (unsigned int m = 1; m < 16; ++m)
      {
        unsigned int codim = !(m&1) + !(m&2) + !(m&4) + !(m&8);
        subindex[m] = codim_count[codim]++;
      }

      int a1 = (3*k + 12)*k + 11;
      int a2 = -3*k - 6;
      unsigned int dof_count[16] = {0};
      unsigned int i[4];
      for (i[3] = 0; i[3] <= k; ++i[3])
        for (i[2] = 0; i[2] <= k - i[3]; ++i[2])
          for (i[1] = 0; i[1] <= k - i[2] - i[3]; ++i[1])
          {
            i[0] = k - i[1] - i[2] - i[3];
            unsigned int j[4];
            unsigned int entity = 0;
            unsigned int codim = 0;
            for (unsigned int m = 0; m < 4; ++m)
            {
              j[m] = i[vertexMap[m]];
              entity += !!j[m] << m;
              codim += !j[m];
            }
            int local_index = j[3]*(a1 + (a2 + j[3])*j[3])/6
                              + j[2]*(2*(k - j[3]) + 3 - j[2])/2 + j[1];
            localKeys_[local_index] = LocalKey(subindex[entity], codim, dof_count[entity]++);
          }
    }
  };

  /** \brief Evaluate the degrees of freedom of a Lagrange basis
   *
   * \tparam LocalBasis The corresponding set of shape functions
   */
  template<class LocalBasis>
  class LagrangeSimplexLocalInterpolation
  {
    static const int kdiv = (LocalBasis::order() == 0 ? 1 : LocalBasis::order());
  public:

    /** \brief Evaluate a given function at the Lagrange nodes
     *
     * \tparam F Type of function to evaluate
     * \tparam C Type used for the values of the function
     * \param[in] f Function to evaluate
     * \param[out] out Array of function values
     */
    template<typename F, typename C>
    void interpolate (const F& f, std::vector<C>& out) const
    {
      constexpr auto dim = LocalBasis::Traits::dimDomain;
      constexpr auto k = LocalBasis::order();
      using D = typename LocalBasis::Traits::DomainFieldType;

      typename LocalBasis::Traits::DomainType x;

      out.resize(LocalBasis::size());

      // Specialization for zero-order case
      if (k==0)
      {
        auto center = ReferenceElements<D,dim>::simplex().position(0,0);
        out[0] = f(center);
        return;
      }

      // Specialization for first-order case
      if (k==1)
      {
        // vertex 0
        std::fill(x.begin(), x.end(), 0);
        out[0] = f(x);

        // remaining vertices
        for (int i=0; i<dim; i++)
        {
          for (int j=0; j<dim; j++)
            x[j] = (i==j);

          out[i+1] = f(x);
        }
        return;
      }

      if (dim==1)
      {
        for (unsigned int i=0; i<k+1; i++)
        {
          x[0] = ((D)i)/k;
          out[i] = f(x);
        }
        return;
      }

      if (dim==2)
      {
        int n=0;
        for (unsigned int j=0; j<=k; j++)
          for (unsigned int i=0; i<=k-j; i++)
          {
            x = { ((D)i)/k, ((D)j)/k };
            out[n] = f(x);
            n++;
          }
        return;
      }

      if (dim!=3)
        DUNE_THROW(NotImplemented, "LagrangeSimplexLocalInterpolation only implemented for dim<=3!");

      int n=0;
      for (int i2 = 0; i2 <= (int)k; i2++)
        for (int i1 = 0; i1 <= (int)k-i2; i1++)
          for (int i0 = 0; i0 <= (int)k-i1-i2; i0++)
          {
            x[0] = ((D)i0)/((D)kdiv);
            x[1] = ((D)i1)/((D)kdiv);
            x[2] = ((D)i2)/((D)kdiv);
            out[n] = f(x);
            n++;
          }
    }

  };

} }    // namespace Dune::Impl

namespace Dune
{
  /** \brief Lagrange finite element for simplices with arbitrary compile-time dimension and polynomial order
   *
   * \tparam D type used for domain coordinates
   * \tparam R type used for function values
   * \tparam d dimension of the reference element
   * \tparam k polynomial order
   *
   * The Lagrange basis functions \f$\phi_i\f$ of order \f$k>1\f$ on the unit simplex
   * \f$G = \{ x \in [0,1]^{d} \,|\, \sum_{j=1}^d x_j \leq 1\}\f$ are implemented as
   * \f$\phi_i(x) = \hat{\phi}_i(kx)\f$ where \f$\hat{\phi}_i\f$ is the \f$i\f$-th
   * basis function on the scaled simplex
   * \f$\hat{G} = kG = \{ x \in [0,k]^{d} \,|\, \sum_{j=1}^d x_i \leq k\}\f$.
   * On \f$\hat{G}\f$ the uniform Lagrange points of order \f$k\f$ are
   * exactly the points \f$(i_1,\dots,i_d) \in \mathbb{N}_0^d \cap \hat{G}\f$
   * with non-negative integer coordinates in \f$\hat{G}\f$.
   * Using the lexicographic enumeration
   * of those points we can identify each \f$i=(i_1,...,i_d)\f$ with the
   * flat index \f$i\f$ of the Lagrange node and associated basis function.

   * Since we can map any point \f$\hat{x} \in \hat{G}\f$ to its barycentric coordinates
   * \f$(\hat{x}_1, \dots, \hat{x}_d, k-\sum_{j=1}^d \hat{x}_j)\f$ we define for any such point
   * its auxiliary \f$(d+1)\f$-th coordinate \f$\hat{x}_{d+1} = k-\sum_{j=1}^d \hat{x}_j\f$
   * and use this in particular for Lagrange points \f$i=(i_1,...,i_d)\f$.
   * Then the Lagrange basis function on \f$\hat{G}\f$ associated to the
   * Lagrange point \f$i=(i_1,\dots,i_d)\f$ is given by
   * \f[
   *   \hat{\phi}_i(\hat{x}) = \prod_{j=1}^{d+1} L_{i_j}(\hat{x}_j)
   * \f]
   * where we used the barycentric coordinates of \f$\hat{x}\f$ and \f$i\f$ and the
   * univariate Lagrange polynomials
   * \f[
   *   L_n(t) = \prod_{m=0}^{n-1} \frac{t-m}{n-m} = \frac{1}{n!}\prod_{m=0}^{n-1}(t-m).
   * \f]
   * This factorization can be interpreted geometrically. To this end note that
   * any component \f$1,\dots,d+1\f$ of the barycentric coordinates is associated
   * to a vertex of the simplex and thus to the opposite facet. Furthermore the Lagrange
   * nodes are all located on hyperplanes parallel to those facets. In the factorized
   * form the term \f$L_{i_j}(\hat{x}_j)\f$ guarantees that \f$\hat{\phi}_i\f$ is zero
   * on any of these hyperplane that is parallel to the facet associated with \f$j\f$
   * and lying between the Lagrange node \f$i\f$ and this facet
   * (excluding \f$i\f$ and including the facet).
   *
   * Using the factorized form, evaluating all basis functions \f$\hat{\phi}_i\f$ at a point \f$\hat{x}\f$
   * requires to first evaluate \f$\alpha_{m,j}=L_m(\hat{x}_j)\f$ for all \f$j\in \{1,\dots,d+1\}\f$
   * and \f$m \in \{0,\dots,k\}\f$ and then computing the products
   * \f$\hat{\phi}_i(\hat{x})=\prod_{j=1}^{d+1} \alpha_{i_j,j}\f$ for all admissible multi indices
   * (or, equivalently, Lagrange nodes in barycentric coordinates)
   * \f$(i_1,\dots,i_{d+1})\f$ with \f$\sum_{j=1}^{d+1} i_j = k\f$.
   * The evaluation of the univariate Lagrange polynomials can be done efficiently using the recursion
   * \f[
   *   L_0(t) = 1, \qquad L_{n+1}(t) = L_n(t)\frac{t-n}{n+1} \qquad n\geq 0.
   * \f]
   */
  template<class D, class R, int d, int k>
  class LagrangeSimplexLocalFiniteElement
  {
  public:
    /** \brief Export number types, dimensions, etc.
     */
    using Traits = LocalFiniteElementTraits<Impl::LagrangeSimplexLocalBasis<D,R,d,k>,
                                            Impl::LagrangeSimplexLocalCoefficients<d,k>,
                                            Impl::LagrangeSimplexLocalInterpolation<Impl::LagrangeSimplexLocalBasis<D,R,d,k> > >;

    /** Default-construct the finite element */
    LagrangeSimplexLocalFiniteElement() {}

    /** \brief Constructs a finite element given a vertex reordering
     * */
    template<typename VertexMap>
    LagrangeSimplexLocalFiniteElement(const VertexMap& vertexmap)
      : coefficients_(vertexmap)
    {}

    /** \brief Returns the local basis, i.e., the set of shape functions
     */
    const typename Traits::LocalBasisType& localBasis () const
    {
      return basis_;
    }

    /** \brief Returns the assignment of the degrees of freedom to the element subentities
     */
    const typename Traits::LocalCoefficientsType& localCoefficients () const
    {
      return coefficients_;
    }

    /** \brief Returns object that evaluates degrees of freedom
     */
    const typename Traits::LocalInterpolationType& localInterpolation () const
    {
      return interpolation_;
    }

    /** \brief The number of shape functions */
    static constexpr std::size_t size ()
    {
      return Impl::LagrangeSimplexLocalBasis<D,R,d,k>::size();
    }

    /** \brief The reference element that the local finite element is defined on
     */
    static constexpr GeometryType type ()
    {
      return GeometryTypes::simplex(d);
    }

  private:
    Impl::LagrangeSimplexLocalBasis<D,R,d,k> basis_;
    Impl::LagrangeSimplexLocalCoefficients<d,k> coefficients_;
    Impl::LagrangeSimplexLocalInterpolation<Impl::LagrangeSimplexLocalBasis<D,R,d,k> > interpolation_;
  };

}        // namespace Dune

#endif   // DUNE_LOCALFUNCTIONS_LAGRANGE_LAGRANGESIMPLEX_HH
