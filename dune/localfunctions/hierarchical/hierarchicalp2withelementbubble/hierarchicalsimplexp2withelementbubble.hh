// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_HIERARCHICAL_SIMPLEX_P2_WITH_ELEMENT_BUBBLE_LOCALBASIS_HH
#define DUNE_HIERARCHICAL_SIMPLEX_P2_WITH_ELEMENT_BUBBLE_LOCALBASIS_HH

/** \file
    \brief Hierarchical p2 shape functions for the simplex
 */

#include <array>
#include <cassert>
#include <numeric>
#include <stdexcept>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/math.hh>

#include <dune/geometry/referenceelement.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localkey.hh>

namespace Dune
{
  /**
   * \ingroup LocalBasisImplementation
   * \brief P1 basis in dim-d enriched by quadratic edge bubble functions and an
   *        element bubble function of order dim+1.
   *
   * The shape functions are associated to the vertices, the edges centers, and to the element barycenter.
   * The bubble function is the product of the linear vertex basis functions. On the edges the basis
   * functions are scaled by the factor 4 and inside the element by the factor (dim+1)^(dim+1).
   *
   * \tparam D    Type to represent the field in the domain.
   * \tparam R    Type to represent the field in the range.
   * \tparam dim  Dimension of the domain.
   *
   * \nosubgrouping
   **/
  template<class D, class R, int dim>
  class HierarchicalSimplexP2WithElementBubbleLocalBasis
  {
    template<class,int> friend class HierarchicalSimplexP2WithElementBubbleLocalInterpolation;

    //! Type of the local coordinates in the domain
    using DomainType = FieldVector<D,dim>;

    //! Range type of the basis functions
    using RangeType = FieldVector<R,1>;

    //! Type of the jacobians of the basis functions
    using JacobianType = FieldMatrix<R,1,dim>;

    // Number of vertices
    static constexpr int numVertices = dim+1;

    // Number of edges (or zero for dim==1)
    static constexpr int numEdges = (dim > 1 ? ((dim+1)*dim / 2) : 0);

    // helper function to evaluate the vertex basis functions
    template <class It>
    static constexpr It evaluateVertexFunctions (const DomainType& in, It outIt)
    {
      *outIt = 1;
      for (int i = 0; i < dim; ++i)
        *outIt -= in[i];
      ++outIt;
      for (int i = 0; i < dim; ++i)
        *outIt++ = in[i];
      return outIt;
    }

    // helper function to evaluate the basis functions
    template <class It>
    static constexpr It evaluateAllFunctions (const DomainType& in, It outIt)
    {
      It vertexValues = outIt;
      outIt = evaluateVertexFunctions(in, outIt);

      if constexpr(dim > 1) {
        auto refElem = referenceElement<D,dim>(GeometryTypes::simplex(dim));
        for (int i = 0; i < numEdges; ++i) {
          const int v0 = refElem.subEntity(i,dim-1,0,dim);
          const int v1 = refElem.subEntity(i,dim-1,1,dim);
          *outIt++ = 4 * vertexValues[v0] * vertexValues[v1];
        }
      }

      // element bubble function
      *outIt = power(dim+1, dim+1);
      for (int i = 0; i < numVertices; ++i)
        *outIt *= vertexValues[i];
      return outIt;
    }

  public:
    //! Type traits for function signature
    using Traits = LocalBasisTraits<D,dim,DomainType,R,1,RangeType,JacobianType>;

    //! Returns number of shape functions
    static constexpr std::size_t size () noexcept
    {
      return numVertices + numEdges + 1;
    }

    //! Evaluate all shape functions
    static constexpr void evaluateFunction (const DomainType& in,
                                            std::vector<RangeType>& out)
    {
      out.resize(size());
      evaluateAllFunctions(in,out.begin());
    }

    //! Evaluate Jacobian of all shape functions
    static constexpr void evaluateJacobian (const DomainType& in,
                                            std::vector<JacobianType>& out)
    {
      out.resize(size());

      // vertex basis functions
      RangeType tmp = 1;
      for (int i = 0; i < dim; ++i) {
        out[0][0][i] = -1;
        for (int j = 0; j < dim; ++j)
          out[j+1][0][i] = (i == j);
        tmp -= in[i];
      }

      int n = numVertices;
      std::array<RangeType,numVertices> shapeValues;
      evaluateVertexFunctions(in, shapeValues.begin());

      // edge basis functions
      if constexpr(dim > 1) {
        auto refElem = referenceElement<D,dim>(GeometryTypes::simplex(dim));
        for (int i = 0; i < numEdges; ++i,++n) {
          const int v0 = refElem.subEntity(i,dim-1,0,dim);
          const int v1 = refElem.subEntity(i,dim-1,1,dim);
          for (int j = 0; j < dim; ++j)
            out[n][0][j] = 4 * (out[v0][0][j] * shapeValues[v1] + shapeValues[v0] * out[v1][0][j]);
        }
      }

      // element bubble function
      for (int i = 0; i < dim; ++i) {
        out[n][0][i] = power(dim+1, dim+1) * (tmp - in[i]);
        for (int j = i+1; j < dim+i; ++j)
          out[n][0][i] *= in[j % dim];
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

        // vertex basis functions
        RangeType tmp = 1;
        for (int i = 0; i < dim; ++i) {
          out[0] = -1;
          for (int j = 0; j < dim; ++j)
            out[j+1] = (dim == j);
          tmp -= in[i];
        }

        int n = numVertices;
        std::array<RangeType,numVertices> shapeValues;
        evaluateVertexFunctions(in, shapeValues.begin());

        // edge basis functions
        if constexpr(dim > 1) {
          auto refElem = referenceElement<D,dim>(GeometryTypes::simplex(dim));
          for (int i = 0; i < numEdges; ++i,++n) {
            const int v0 = refElem.subEntity(i,dim-1,0,dim);
            const int v1 = refElem.subEntity(i,dim-1,1,dim);
            out[n] = 4 * (out[v0] * shapeValues[v1] + shapeValues[v0] * out[v1]);
          }
        }

        // element bubble function
        out[n] = power(dim+1, dim+1) * (tmp - in[d]);
        for (int j = d+1; j < dim+d; ++j)
          out[n] *= in[j % dim];
      } break;
      default:
        throw std::runtime_error("Desired derivative order is not implemented");
      }
    }

    //! Polynomial order of the shape functions (4 in this case)
    static constexpr unsigned int order () noexcept
    {
      return dim+1;
    }
  };


  /**
   * \ingroup LocalBasisImplementation
   * \brief The local keys of the hierarchical basis functions with element bubble.
   *
   * This shape function set consists of three parts:
   * - Linear shape functions associated to the element vertices
   * - Quadratic edge bubbles
   * - An element bubble of order dim+1
   *
   * \nosubgrouping
   **/
  template <int dim>
  class HierarchicalSimplexP2WithElementBubbleLocalCoefficients
  {
    // Number of vertices
    static constexpr int numVertices = dim+1;

    // Number of edges (or zero for dim==1)
    static constexpr int numEdges = (dim > 1 ? ((dim+1)*dim / 2) : 0);

  public:
    //! Default constructor, initializes the local keys
    HierarchicalSimplexP2WithElementBubbleLocalCoefficients () noexcept
    {
      int n = 0;
      for (int i = 0; i < numVertices; ++i)
        li_[n++] = LocalKey(i,dim,0);       // Vertex
      if constexpr(dim > 1) {
        for (int i = 0; i < numEdges; ++i)
          li_[n++] = LocalKey(i,dim-1,0);   // Edges
      }
      li_[n++] = LocalKey(0,0,0);           // Element
    }

    //! Returns number of coefficients
    static constexpr std::size_t size () noexcept
    {
      return numVertices + numEdges + 1;
    }

    //! Returns the i'th local key
    const LocalKey& localKey (std::size_t i) const noexcept
    {
      return li_[i];
    }

  private:
    std::array<LocalKey, numVertices+numEdges+1> li_;
  };

  /**
   * \ingroup LocalFunctionsImpl
   */
  template<class LB, int dim>
  class HierarchicalSimplexP2WithElementBubbleLocalInterpolation
  {
    using LocalBasis = LB;
    using DomainType = typename LB::Traits::DomainType;
    using RangeType = typename LB::Traits::RangeType;

    // Number of vertices
    static constexpr int numVertices = dim+1;

    // Number of edges (or zero for dim==1)
    static constexpr int numEdges = (dim > 1 ? ((dim+1)*dim / 2) : 0);

  public:
    /**
     * \brief Local interpolation of the function `f`.
     * \param[in] f  A callable `f:D -> R` with domain `D=DomainType` and range
     *               `R` convertible into the coefficient type `C`.
     * \param[out] out  The interpolation coefficients `{f_i...,f_b}` are stored
     *                  in this output vector.
     **/
    template<class F, class C,
      class R = std::invoke_result_t<F, DomainType>,
      std::enable_if_t<std::is_convertible_v<R, C>, int> = 0>
    static constexpr void interpolate (const F& f, std::vector<C>& out)
    {
      auto refElem = referenceElement<typename LB::Traits::DomainFieldType,dim>(GeometryTypes::simplex(dim));

      out.resize(LB::size());
      int n = 0;

      // vertices
      assert(numVertices == refElem.size(dim));
      for (int i = 0; i < numVertices; ++i)
        out[n++] = f(refElem.position(i,dim));

      std::array<RangeType,LB::size()> shapeValues;

      // edge bubbles
      if constexpr(dim > 1) {
        assert(numEdges == refElem.size(dim-1));
        for (int i = 0; i < numEdges; ++i) {
          R y = f(refElem.position(i,dim-1));
          LB::evaluateVertexFunctions(refElem.position(i,dim-1), shapeValues.begin());
          for (int j = 0; j < numVertices; ++j)
            y -= out[j]*shapeValues[j];
          out[n++] = y;
        }
      }

      // element bubble
      R y = f(refElem.position(0,0));
      LB::evaluateAllFunctions(refElem.position(0,0), shapeValues.begin());
      for (int j = 0; j < numVertices+numEdges; ++j)
        y -= out[j]*shapeValues[j];
      out[n++] = y;
    }
  };

} // end namespace Dune

#endif // DUNE_HIERARCHICAL_SIMPLEX_P2_WITH_ELEMENT_BUBBLE_LOCALBASIS_HH
