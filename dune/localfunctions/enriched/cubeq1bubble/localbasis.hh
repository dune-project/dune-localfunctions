// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_ENRICHED_CUBEQ1BUBBLE_LOCALBASIS_HH
#define DUNE_LOCALFUNCTIONS_ENRICHED_CUBEQ1BUBBLE_LOCALBASIS_HH

#include <numeric>
#include <stdexcept>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/math.hh>

#include <dune/localfunctions/common/localbasis.hh>

namespace Dune
{
  /**
   * \ingroup LocalBasisImplementation
   * \brief Q1 basis in dim-d enriched by an (order 2) element bubble function.
   *
   * The shape functions are associated to the vertices and to the barycenter. The
   * bubble function is the product of the edge bubble basis functions scaled by
   * the factor 2^(2^dim). The edge bubble functions are the product of the
   * corresponding vertex basis functions.
   *
   * \tparam D Type to represent the field in the domain.
   * \tparam R Type to represent the field in the range.
   *
   * \nosubgrouping
   **/
  template<class D, class R, int dim>
  class CubeQ1BubbleLocalBasis
  {
    template<class> friend class CubeQ1BubbleLocalInterpolation;

    //! Type of the local coordinates in the domain
    using DomainType = FieldVector<D,dim>;

    //! Range type of the basis functions
    using RangeType = FieldVector<R,1>;

    //! Type of the jacobians of the basis functions
    using JacobianType = FieldMatrix<R,1,dim>;

    static constexpr int dimension = dim;
    static constexpr int numVertices = power(2, dim);

    // scaling factor of bubble basis function for normalization
    static constexpr int scaling = power(2, 2*dim);

  public:
    //! Type traits for function signature
    using Traits = LocalBasisTraits<D,dim,DomainType,R,1,RangeType,JacobianType>;

    //! Returns number of shape functions
    static constexpr std::size_t size () noexcept
    {
      return numVertices+1;
    }

    /// Evaluate all shape functions
    static constexpr void evaluateFunction (const DomainType& x, std::vector<RangeType>& out)
    {
      out.resize(size());

      for (int i = 0; i < numVertices; ++i) {
        out[i] = 1;
        for (int j = 0; j < dimension; ++j) {
          // if j-th bit of i is set multiply with x[j], else with 1-x[j]
          out[i] *= (i & (1<<j)) ? x[j] :  1-x[j];
        }
      }

      out.back() = scaling;
      for (int j = 0; j < dimension; ++j) {
        out.back() *= (1-x[j]) * x[j];
      }
    }

    /// Evaluate Jacobian of all shape functions
    static constexpr void evaluateJacobian (const DomainType& x, std::vector<JacobianType>& out)
    {
      out.resize(size());

      // Loop over all linear shape functions
      for (int i = 0; i < numVertices; ++i) {
        // Loop over all coordinate directions
        for (int j = 0; j < dimension; ++j) {
          // Initialize: the overall expression as a product
          // if j-th bit of i is set to 1, else -11
          out[i][0][j] = (i & (1<<j)) ? 1 : -1;

          for (int k = 0; k < dimension; ++k) {
            if (j != k)
              // if k-th bit of i is set multiply with x[k], else with 1-x[k]
              out[i][0][j] *= (i & (1<<k)) ? x[k] :  1-x[k];
          }
        }
      }

      for (int j = 0; j < dimension; ++j) {
        out.back()[0][j] = scaling;
        for (int k = 0; k < dimension; ++k)
          out.back()[0][j] *= (j == k) ? (1-2*x[k]) : (1-x[k])*x[k];
      }
    }

    //! Evaluate partial derivatives of all shape functions
    static constexpr void partial (const std::array<unsigned int, dim>& order,
                                   const DomainType& x,
                                   std::vector<RangeType>& out)
    {
      unsigned int totalOrder = 0;
      for (int i = 0; i < dimension; ++i)
        totalOrder += order[i];

      switch (totalOrder) {
      case 0:
        evaluateFunction(x,out);
        break;
      case 1: {
        out.resize(size());
        int d = 0; // the direction of differentiation
        for (int i = 0; i < dimension; ++i)
          d += i * order[i];

        if (d >= dimension)
          throw std::invalid_argument("Direction of partial derivative not found!");

        // Loop over all shape functions
        for (int i = 0; i < numVertices; ++i) {
          // Initialize: the overall expression is a product
          // if j-th bit of i is set to 1, otherwise to -1
          out[i] = (i & (1<<d)) ? 1 : -1;

          for (int j = 0; j < dimension; ++j) {
            if (d != j)
              // if j-th bit of i is set multiply with in[j], else with 1-in[j]
              out[i] *= (i & (1<<j)) ? x[j] :  1-x[j];
          }
        }

        out.back() = scaling;
        for (int k = 0; k < dimension; ++k)
          out.back() *= (d == k) ? (1-2*x[k]) : (1-x[k])*x[k];

      } break;
      default:
        throw std::runtime_error("Desired derivative order is not implemented");
      }
    }

    //! Returns maximal polynomial order of the basis functions
    static constexpr unsigned int order () noexcept
    {
      return 2;
    }
  };

} // end namespace Dune

#endif // DUNE_LOCALFUNCTIONS_ENRICHED_CUBEQ1BUBBLE_LOCALBASIS_HH
