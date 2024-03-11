// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_ENRICHED_SIMPLEXP1BUBBLE_LOCALBASIS_HH
#define DUNE_LOCALFUNCTIONS_ENRICHED_SIMPLEXP1BUBBLE_LOCALBASIS_HH

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
   * \brief P1 basis in dim-d enriched by an (order dim+1) element bubble function.
   *
   * The shape functions are associated to the vertices and to the barycenter. The
   * bubble function is the product of the linear vertex basis functions scaled by
   * the factor (dim+1)^(dim+1).
   *
   * \tparam D Type to represent the field in the domain.
   * \tparam R Type to represent the field in the range.
   *
   * \nosubgrouping
   **/
  template<class D, class R, int dim>
  class SimplexP1BubbleLocalBasis
  {
    template<class> friend class SimplexP1BubbleLocalInterpolation;

    //! Type of the local coordinates in the domain
    using DomainType = FieldVector<D,dim>;

    //! Range type of the basis functions
    using RangeType = FieldVector<R,1>;

    //! Type of the jacobians of the basis functions
    using JacobianType = FieldMatrix<R,1,dim>;

    static constexpr int dimension = dim;
    static constexpr int numVertices = dim+1;

  public:
    //! Type traits for function signature
    using Traits = LocalBasisTraits<D,dim,DomainType,R,1,RangeType,JacobianType>;

    //! Returns number of shape functions
    static constexpr std::size_t size () noexcept
    {
      return numVertices+1;
    }

    //! Evaluate all shape functions
    static constexpr void evaluateFunction (const DomainType& in,
                                            std::vector<RangeType>& out)
    {
      out.resize(size());
      out[0] = 1;
      out.back() = power(dim+1, dim+1); // normalization of the bubble function
      for (int i = 0; i < dim; ++i) {
        out[0]     -= in[i];
        out[i+1]    = in[i];
        out.back() *= in[i];
      }
      out.back() *= out[0];
    }

    //! Evaluate Jacobian of all shape functions
    static constexpr void evaluateJacobian (const DomainType& in,
                                            std::vector<JacobianType>& out)
    {
      out.resize(size());
      RangeType tmp = 1;
      for (int i = 0; i < dim; ++i) {
        out[0][0][i] = -1;
        for (int j = 0; j < dim; ++j)
          out[j+1][0][i] = (i == j);
        tmp -= in[i];
      }

      for (int i = 0; i < dim; ++i) {
        out.back()[0][i] = power(dim+1, dim+1) * (tmp - in[i]);
        for (int j = i+1; j < dim+i; ++j)
          out.back()[0][i] *= in[j % dim];
      }
    }

    //! Evaluate partial derivatives of all shape functions
    static constexpr void partial (const std::array<unsigned int, dim>& order,
                                   const DomainType& in,
                                   std::vector<RangeType>& out)
    {
      unsigned int totalOrder = 0;
      for (int i = 0; i < dim; ++i)
        totalOrder += order[i];

      switch (totalOrder) {
      case 0:
        evaluateFunction(in,out);
        break;
      case 1: {
        out.resize(size());
        int d = 0; // the direction of differentiation
        for (int i = 0; i < dim; ++i)
          d += i * order[i];

        out[0] = -1;
        RangeType tmp = 1;
        for (int j = 0; j < dim; ++j) {
          out[j+1] = (d == j);
          tmp -= in[j];
        }
        out.back() = power(dim+1, dim+1) * (tmp - in[d]);
        for (int j = d+1; j < dim+d; ++j)
          out.back() *= in[j % dim];
      } break;
      default:
        throw std::runtime_error("Desired derivative order is not implemented");
      }
    }

    //! Returns maximal polynomial order of the basis functions
    static constexpr unsigned int order () noexcept
    {
      return dim+1;
    }
  };

} // end namespace Dune

#endif // DUNE_LOCALFUNCTIONS_ENRICHED_SIMPLEXP1BUBBLE_LOCALBASIS_HH
