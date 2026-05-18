// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_LAGRANGE_LAGRANGEPRISM_HH
#define DUNE_LOCALFUNCTIONS_LAGRANGE_LAGRANGEPRISM_HH

#include <array>
#include <numeric>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/math.hh>

#include <dune/geometry/referenceelements.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localkey.hh>

namespace Dune { namespace Impl
{

  // The traits provide static or dynamic order and size information
  template<int compileTimeOrder>
  struct LagrangePrismOrderTraits;

  // The traits provide static order and size information
  template<int compileTimeOrder>
  requires (compileTimeOrder >= 0)
  struct LagrangePrismOrderTraits<compileTimeOrder>
  {
    static constexpr bool is_static_order = true;

    constexpr LagrangePrismOrderTraits(int /*runTimeOrder*/ = compileTimeOrder)
    {
      static_assert((0 <= compileTimeOrder) and (compileTimeOrder <= 2), "LagrangePrism: Only order 0,1, and 2 are supported");
    }

    /**
     * \brief Polynomial order of the shape functions
     */
    static constexpr unsigned int order() {
      return compileTimeOrder;
    }

    /**
     * \brief Number of shape functions
     */
    static constexpr std::size_t size() {
      return binomial(compileTimeOrder+2,2) * (compileTimeOrder+1);
    }
  };


  // this specialization is for the dynamic order case
  template<>
  struct LagrangePrismOrderTraits<-1>
  {
    unsigned int runTimeOrder_;
    unsigned int size_;

    static constexpr bool is_static_order = false;

    constexpr explicit LagrangePrismOrderTraits(int runTimeOrder)
      : runTimeOrder_(runTimeOrder >= 0 ? (unsigned int)(runTimeOrder) : 0u)
      , size_(binomial(runTimeOrder+2,2) * (runTimeOrder+1))
    {
      if ((runTimeOrder < 0) or (runTimeOrder > 2))
        DUNE_THROW(Dune::InvalidStateException, "LagrangePrism: Only order 0,1, and 2 are supported");
    }

    /**
     * \brief Polynomial order of the shape functions
     */
    constexpr unsigned int order() const
    {
      return runTimeOrder_;
    }

    /**
     * \brief Number of shape functions
     */
    constexpr std::size_t size() const
    {
      return size_;
    }
  };

   /** \brief Lagrange shape functions of arbitrary order on the three-dimensional reference prism

     Lagrange shape functions of arbitrary order have the property that
     \f$\hat\phi^i(x_j) = \delta_{i,j}\f$ for certain points \f$x_j\f$.

     \tparam D Type to represent the field in the domain
     \tparam R Type to represent the field in the range
     \tparam compileTimeOrder Polynomial order of the Lagrange space in one direction (or -1 for dynamic order)
   */
  template<class D, class R, int compileTimeOrder>
  class LagrangePrismLocalBasis
    : public LagrangePrismOrderTraits<compileTimeOrder>
  {
    using OrderTraits = LagrangePrismOrderTraits<compileTimeOrder>;
    static constexpr std::size_t dim = 3;
  public:
    using Traits = LocalBasisTraits<D,dim,FieldVector<D,dim>,R,1,FieldVector<R,1>,FieldMatrix<R,1,dim> >;

    constexpr LagrangePrismLocalBasis(OrderTraits orderTraits)
      : OrderTraits(orderTraits)
    {}

    using OrderTraits::order;
    using OrderTraits::size;

    //! \brief Evaluate all shape functions
    void evaluateFunction(const typename Traits::DomainType& in,
                          std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(size());

      // Specialization for zero-order case
      if (order()==0)
      {
        out[0] = 1;
        return;
      }

      if (order()==1)
      {
        out[0] = (1.0-in[0]-in[1])*(1.0-in[2]);
        out[1] = in[0]*(1-in[2]);
        out[2] = in[1]*(1-in[2]);
        out[3] = in[2]*(1.0-in[0]-in[1]);
        out[4] = in[0]*in[2];
        out[5] = in[1]*in[2];

        return;
      }

      if (order()==2)
      {
        FieldVector<R,3> segmentShapeFunction;
        segmentShapeFunction[0] = 1 + in[2] * (-3 + 2*in[2]);
        segmentShapeFunction[1] =    in[2] * (4 - 4*in[2]);
        segmentShapeFunction[2] =    in[2] * (-1 + 2*in[2]);

        FieldVector<R, 6> triangleShapeFunction;
        triangleShapeFunction[0] = 2 * (1 - in[0] - in[1]) * (0.5 - in[0] - in[1]);
        triangleShapeFunction[1] = 2 * in[0] * (-0.5 + in[0]);
        triangleShapeFunction[2] = 2 * in[1] * (-0.5 + in[1]);
        triangleShapeFunction[3] = 4*in[0] * (1 - in[0] - in[1]);
        triangleShapeFunction[4] = 4*in[1] * (1 - in[0] - in[1]);
        triangleShapeFunction[5] = 4*in[0]*in[1];

        // lower triangle:
        out[0] = triangleShapeFunction[0] * segmentShapeFunction[0];
        out[1] = triangleShapeFunction[1] * segmentShapeFunction[0];
        out[2] = triangleShapeFunction[2] * segmentShapeFunction[0];

        //upper triangle
        out[3] = triangleShapeFunction[0] * segmentShapeFunction[2];
        out[4] = triangleShapeFunction[1] * segmentShapeFunction[2];
        out[5] = triangleShapeFunction[2] * segmentShapeFunction[2];

        // vertical edges
        out[6] = triangleShapeFunction[0] * segmentShapeFunction[1];
        out[7] = triangleShapeFunction[1] * segmentShapeFunction[1];
        out[8] = triangleShapeFunction[2] * segmentShapeFunction[1];

        // lower triangle edges
        out[9] = triangleShapeFunction[3] * segmentShapeFunction[0];
        out[10] = triangleShapeFunction[4] * segmentShapeFunction[0];
        out[11] = triangleShapeFunction[5] * segmentShapeFunction[0];

        // upper triangle edges
        out[12] = triangleShapeFunction[3] * segmentShapeFunction[2];
        out[13] = triangleShapeFunction[4] * segmentShapeFunction[2];
        out[14] = triangleShapeFunction[5] * segmentShapeFunction[2];

        // quadrilateral sides
        out[15] = triangleShapeFunction[3] * segmentShapeFunction[1];
        out[16] = triangleShapeFunction[4] * segmentShapeFunction[1];
        out[17] = triangleShapeFunction[5] * segmentShapeFunction[1];

        return;
      }

      DUNE_THROW(NotImplemented, "LagrangePrismLocalBasis::evaluateFunction for order " << order());
    }

    /** \brief Evaluate Jacobian of all shape functions
     *
     * \param x Point in the reference prism where to evaluation the Jacobians
     * \param[out] out The Jacobians of all shape functions at the point x
     */
    void evaluateJacobian(const typename Traits::DomainType& in,
                          std::vector<typename Traits::JacobianType>& out) const
    {
      out.resize(size());

      // Specialization for k==0
      if (order()==0)
      {
        std::fill(out[0][0].begin(), out[0][0].end(), 0);
        return;
      }

      if (order()==1)
      {
        out[0][0] = {in[2]-1, in[2]-1, in[0]+in[1]-1};
        out[1][0] = {1-in[2],       0,        -in[0]};
        out[2][0] = {      0, 1-in[2],        -in[1]};
        out[3][0] = { -in[2],  -in[2], 1-in[0]-in[1]};
        out[4][0] = {  in[2],       0,         in[0]};
        out[5][0] = {      0,   in[2],         in[1]};

        return;
      }

      if (order()==2)
      {
        // Second-order shape functions on a triangle, and the first derivatives
        FieldVector<R, 6> triangleShapeFunction;
        triangleShapeFunction[0] = 2 * (1 - in[0] - in[1]) * (0.5 - in[0] - in[1]);
        triangleShapeFunction[1] = 2 * in[0] * (-0.5 + in[0]);
        triangleShapeFunction[2] = 2 * in[1] * (-0.5 + in[1]);
        triangleShapeFunction[3] = 4*in[0] * (1 - in[0] - in[1]);
        triangleShapeFunction[4] = 4*in[1] * (1 - in[0] - in[1]);
        triangleShapeFunction[5] = 4*in[0]*in[1];

        std::array<std::array<R,2>,6> triangleShapeFunctionDer;
        triangleShapeFunctionDer[0] = {-3 + 4*(in[0] + in[1]), -3 + 4*(in[0] + in[1])};
        triangleShapeFunctionDer[1] = {          -1 + 4*in[0],                      0};
        triangleShapeFunctionDer[2] = {                     0,           -1 + 4*in[1]};
        triangleShapeFunctionDer[3] = { 4 - 8*in[0] - 4*in[1],               -4*in[0]};
        triangleShapeFunctionDer[4] = {              -4*in[1],  4 - 4*in[0] - 8*in[1]};
        triangleShapeFunctionDer[5] = {               4*in[1],                4*in[0]};

        // Second-order shape functions on a line, and the first derivatives
        FieldVector<R,3> segmentShapeFunction;
        segmentShapeFunction[0] = 1 + in[2] * (-3 + 2*in[2]);
        segmentShapeFunction[1] =     in[2] * ( 4 - 4*in[2]);
        segmentShapeFunction[2] =     in[2] * (-1 + 2*in[2]);

        FieldVector<R,3> segmentShapeFunctionDer;
        segmentShapeFunctionDer[0] = -3 + 4*in[2];
        segmentShapeFunctionDer[1] =  4 - 8*in[2];
        segmentShapeFunctionDer[2] = -1 + 4*in[2];

        // lower triangle:
        out[0][0][0] = triangleShapeFunctionDer[0][0] * segmentShapeFunction[0];
        out[0][0][1] = triangleShapeFunctionDer[0][1] * segmentShapeFunction[0];
        out[0][0][2] = triangleShapeFunction[0]       * segmentShapeFunctionDer[0];

        out[1][0][0] = triangleShapeFunctionDer[1][0] * segmentShapeFunction[0];
        out[1][0][1] = triangleShapeFunctionDer[1][1] * segmentShapeFunction[0];
        out[1][0][2] = triangleShapeFunction[1]       * segmentShapeFunctionDer[0];

        out[2][0][0] = triangleShapeFunctionDer[2][0] * segmentShapeFunction[0];
        out[2][0][1] = triangleShapeFunctionDer[2][1] * segmentShapeFunction[0];
        out[2][0][2] = triangleShapeFunction[2]       * segmentShapeFunctionDer[0];

        //upper triangle
        out[3][0][0] = triangleShapeFunctionDer[0][0] * segmentShapeFunction[2];
        out[3][0][1] = triangleShapeFunctionDer[0][1] * segmentShapeFunction[2];
        out[3][0][2] = triangleShapeFunction[0]       * segmentShapeFunctionDer[2];

        out[4][0][0] = triangleShapeFunctionDer[1][0] * segmentShapeFunction[2];
        out[4][0][1] = triangleShapeFunctionDer[1][1] * segmentShapeFunction[2];
        out[4][0][2] = triangleShapeFunction[1]       * segmentShapeFunctionDer[2];

        out[5][0][0] = triangleShapeFunctionDer[2][0] * segmentShapeFunction[2];
        out[5][0][1] = triangleShapeFunctionDer[2][1] * segmentShapeFunction[2];
        out[5][0][2] = triangleShapeFunction[2]       * segmentShapeFunctionDer[2];

        // vertical edges
        out[6][0][0] = triangleShapeFunctionDer[0][0] * segmentShapeFunction[1];
        out[6][0][1] = triangleShapeFunctionDer[0][1] * segmentShapeFunction[1];
        out[6][0][2] = triangleShapeFunction[0]       * segmentShapeFunctionDer[1];

        out[7][0][0] = triangleShapeFunctionDer[1][0] * segmentShapeFunction[1];
        out[7][0][1] = triangleShapeFunctionDer[1][1] * segmentShapeFunction[1];
        out[7][0][2] = triangleShapeFunction[1]       * segmentShapeFunctionDer[1];

        out[8][0][0] = triangleShapeFunctionDer[2][0] * segmentShapeFunction[1];
        out[8][0][1] = triangleShapeFunctionDer[2][1] * segmentShapeFunction[1];
        out[8][0][2] = triangleShapeFunction[2]       * segmentShapeFunctionDer[1];

        // lower triangle edges
        out[9][0][0] = triangleShapeFunctionDer[3][0] * segmentShapeFunction[0];
        out[9][0][1] = triangleShapeFunctionDer[3][1] * segmentShapeFunction[0];
        out[9][0][2] = triangleShapeFunction[3]       * segmentShapeFunctionDer[0];

        out[10][0][0] = triangleShapeFunctionDer[4][0] * segmentShapeFunction[0];
        out[10][0][1] = triangleShapeFunctionDer[4][1] * segmentShapeFunction[0];
        out[10][0][2] = triangleShapeFunction[4]       * segmentShapeFunctionDer[0];

        out[11][0][0] = triangleShapeFunctionDer[5][0] * segmentShapeFunction[0];
        out[11][0][1] = triangleShapeFunctionDer[5][1] * segmentShapeFunction[0];
        out[11][0][2] = triangleShapeFunction[5]       * segmentShapeFunctionDer[0];

        // upper triangle edges
        out[12][0][0] = triangleShapeFunctionDer[3][0] * segmentShapeFunction[2];
        out[12][0][1] = triangleShapeFunctionDer[3][1] * segmentShapeFunction[2];
        out[12][0][2] = triangleShapeFunction[3]       * segmentShapeFunctionDer[2];

        out[13][0][0] = triangleShapeFunctionDer[4][0] * segmentShapeFunction[2];
        out[13][0][1] = triangleShapeFunctionDer[4][1] * segmentShapeFunction[2];
        out[13][0][2] = triangleShapeFunction[4]       * segmentShapeFunctionDer[2];

        out[14][0][0] = triangleShapeFunctionDer[5][0] * segmentShapeFunction[2];
        out[14][0][1] = triangleShapeFunctionDer[5][1] * segmentShapeFunction[2];
        out[14][0][2] = triangleShapeFunction[5]       * segmentShapeFunctionDer[2];

        // quadrilateral sides
        out[15][0][0] = triangleShapeFunctionDer[3][0] * segmentShapeFunction[1];
        out[15][0][1] = triangleShapeFunctionDer[3][1] * segmentShapeFunction[1];
        out[15][0][2] = triangleShapeFunction[3]       * segmentShapeFunctionDer[1];

        out[16][0][0] = triangleShapeFunctionDer[4][0] * segmentShapeFunction[1];
        out[16][0][1] = triangleShapeFunctionDer[4][1] * segmentShapeFunction[1];
        out[16][0][2] = triangleShapeFunction[4]       * segmentShapeFunctionDer[1];

        out[17][0][0] = triangleShapeFunctionDer[5][0] * segmentShapeFunction[1];
        out[17][0][1] = triangleShapeFunctionDer[5][1] * segmentShapeFunction[1];
        out[17][0][2] = triangleShapeFunction[5]       * segmentShapeFunctionDer[1];

        return;
      }

      DUNE_THROW(NotImplemented, "LagrangePrismLocalBasis::evaluateJacobian for order " << order());
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

      if (totalOrder == 0)
      {
        evaluateFunction(in, out);
        return;
      }

      // Specialization for zero-order finite elements
      if (OrderTraits::order()==0)
      {
        out[0] = 0;
        return;
      }

      // Specialization for first-order finite elements
      if (OrderTraits::order()==1)
      {
        if (totalOrder == 1)
        {
          auto direction = std::distance(order.begin(), std::find(order.begin(), order.end(), 1));

          switch (direction) {
            case 0:
              out[0] = in[2]-1;
              out[1] = 1-in[2];
              out[2] = 0;
              out[3] = -in[2];
              out[4] = in[2];
              out[5] = 0;
              break;
            case 1:
              out[0] = in[2]-1;
              out[1] = 0;
              out[2] = 1-in[2];
              out[3] = -in[2];
              out[4] = 0;
              out[5] = in[2];
              break;
            case 2:
              out[0] = in[0]+in[1]-1;
              out[1] = -in[0];
              out[2] = -in[1];
              out[3] = 1-in[0]-in[1];
              out[4] = in[0];
              out[5] = in[1];
              break;
            default:
              DUNE_THROW(RangeError, "Component out of range.");
          }
        } else if (totalOrder == 2) {
          out.resize(size());
          if (order[0] == 1 && order[2] == 1) {
            out[0] = 1;
            out[1] =-1;
            out[2] = 0;
            out[3] =-1;
            out[4] = 1;
            out[5] = 0;
          } else if (order[1] == 1 && order[2] == 1) {
            out[0] = 1;
            out[1] = 0;
            out[2] =-1;
            out[3] =-1;
            out[4] = 0;
            out[5] = 1;
          } else {
            for (std::size_t i = 0; i < size(); ++i)
              out[i] = 0;
          }
        } else {
          out.resize(size());
          std::fill(out.begin(), out.end(), 0.0);
        }

        return;
      }

      // Specialization for second-order finite elements
      if (OrderTraits::order()==2)
      {
        if (totalOrder == 1)
        {
          auto const direction = std::distance(order.begin(), std::find(order.begin(), order.end(), 1));
          switch (direction)
          {
          case 0:
          {
            FieldVector<R,6> triangleShapeFunctionDerX;
            triangleShapeFunctionDerX[0] = -3 + 4*(in[0] +   in[1]);
            triangleShapeFunctionDerX[1] = -1 + 4* in[0];
            triangleShapeFunctionDerX[2] = 0;
            triangleShapeFunctionDerX[3] =  4 - 8* in[0] - 4*in[1];
            triangleShapeFunctionDerX[4] =                -4*in[1];
            triangleShapeFunctionDerX[5] =                 4*in[1];

            FieldVector<R,3> segmentShapeFunction;
            segmentShapeFunction[0] = 1 + in[2] * (-3 + 2*in[2]);
            segmentShapeFunction[1] =     in[2] * ( 4 - 4*in[2]);
            segmentShapeFunction[2] =     in[2] * (-1 + 2*in[2]);

            out[0]  = triangleShapeFunctionDerX[0] * segmentShapeFunction[0];
            out[1]  = triangleShapeFunctionDerX[1] * segmentShapeFunction[0];
            out[2]  = triangleShapeFunctionDerX[2] * segmentShapeFunction[0];
            out[3]  = triangleShapeFunctionDerX[0] * segmentShapeFunction[2];
            out[4]  = triangleShapeFunctionDerX[1] * segmentShapeFunction[2];
            out[5]  = triangleShapeFunctionDerX[2] * segmentShapeFunction[2];
            out[6]  = triangleShapeFunctionDerX[0] * segmentShapeFunction[1];
            out[7]  = triangleShapeFunctionDerX[1] * segmentShapeFunction[1];
            out[8]  = triangleShapeFunctionDerX[2] * segmentShapeFunction[1];
            out[9]  = triangleShapeFunctionDerX[3] * segmentShapeFunction[0];
            out[10] = triangleShapeFunctionDerX[4] * segmentShapeFunction[0];
            out[11] = triangleShapeFunctionDerX[5] * segmentShapeFunction[0];
            out[12] = triangleShapeFunctionDerX[3] * segmentShapeFunction[2];
            out[13] = triangleShapeFunctionDerX[4] * segmentShapeFunction[2];
            out[14] = triangleShapeFunctionDerX[5] * segmentShapeFunction[2];
            out[15] = triangleShapeFunctionDerX[3] * segmentShapeFunction[1];
            out[16] = triangleShapeFunctionDerX[4] * segmentShapeFunction[1];
            out[17] = triangleShapeFunctionDerX[5] * segmentShapeFunction[1];
            break;
          }
          case 1:
          {
            FieldVector<R,6> triangleShapeFunctionDerY;
            triangleShapeFunctionDerY[0] = -3 + 4*(in[0] +   in[1]);
            triangleShapeFunctionDerY[1] = 0;
            triangleShapeFunctionDerY[2] = -1 + 4*           in[1];
            triangleShapeFunctionDerY[3] =     -4* in[0];
            triangleShapeFunctionDerY[4] =  4 - 4* in[0] - 8*in[1];
            triangleShapeFunctionDerY[5] =      4* in[0];

            FieldVector<R,3> segmentShapeFunction;
            segmentShapeFunction[0] = 1 + in[2] * (-3 + 2*in[2]);
            segmentShapeFunction[1] =     in[2] * ( 4 - 4*in[2]);
            segmentShapeFunction[2] =     in[2] * (-1 + 2*in[2]);

            out[0]  = triangleShapeFunctionDerY[0] * segmentShapeFunction[0];
            out[1]  = triangleShapeFunctionDerY[1] * segmentShapeFunction[0];
            out[2]  = triangleShapeFunctionDerY[2] * segmentShapeFunction[0];
            out[3]  = triangleShapeFunctionDerY[0] * segmentShapeFunction[2];
            out[4]  = triangleShapeFunctionDerY[1] * segmentShapeFunction[2];
            out[5]  = triangleShapeFunctionDerY[2] * segmentShapeFunction[2];
            out[6]  = triangleShapeFunctionDerY[0] * segmentShapeFunction[1];
            out[7]  = triangleShapeFunctionDerY[1] * segmentShapeFunction[1];
            out[8]  = triangleShapeFunctionDerY[2] * segmentShapeFunction[1];
            out[9]  = triangleShapeFunctionDerY[3] * segmentShapeFunction[0];
            out[10] = triangleShapeFunctionDerY[4] * segmentShapeFunction[0];
            out[11] = triangleShapeFunctionDerY[5] * segmentShapeFunction[0];
            out[12] = triangleShapeFunctionDerY[3] * segmentShapeFunction[2];
            out[13] = triangleShapeFunctionDerY[4] * segmentShapeFunction[2];
            out[14] = triangleShapeFunctionDerY[5] * segmentShapeFunction[2];
            out[15] = triangleShapeFunctionDerY[3] * segmentShapeFunction[1];
            out[16] = triangleShapeFunctionDerY[4] * segmentShapeFunction[1];
            out[17] = triangleShapeFunctionDerY[5] * segmentShapeFunction[1];
            break;
          }
          case 2:
          {
            FieldVector<R, 6> triangleShapeFunction;
            triangleShapeFunction[0] = 2 * (1 - in[0] - in[1]) * (0.5 - in[0] - in[1]);
            triangleShapeFunction[1] = 2 * in[0] * (-0.5 + in[0]);
            triangleShapeFunction[2] = 2 * in[1] * (-0.5 + in[1]);
            triangleShapeFunction[3] = 4*in[0] * (1 - in[0] - in[1]);
            triangleShapeFunction[4] = 4*in[1] * (1 - in[0] - in[1]);
            triangleShapeFunction[5] = 4*in[0]*in[1];

            FieldVector<R,3> segmentShapeFunctionDer;
            segmentShapeFunctionDer[0] = -3 + 4*in[2];
            segmentShapeFunctionDer[1] =  4 - 8*in[2];
            segmentShapeFunctionDer[2] = -1 + 4*in[2];

            out[0]  = triangleShapeFunction[0] * segmentShapeFunctionDer[0];
            out[1]  = triangleShapeFunction[1] * segmentShapeFunctionDer[0];
            out[2]  = triangleShapeFunction[2] * segmentShapeFunctionDer[0];
            out[3]  = triangleShapeFunction[0] * segmentShapeFunctionDer[2];
            out[4]  = triangleShapeFunction[1] * segmentShapeFunctionDer[2];
            out[5]  = triangleShapeFunction[2] * segmentShapeFunctionDer[2];
            out[6]  = triangleShapeFunction[0] * segmentShapeFunctionDer[1];
            out[7]  = triangleShapeFunction[1] * segmentShapeFunctionDer[1];
            out[8]  = triangleShapeFunction[2] * segmentShapeFunctionDer[1];
            out[9]  = triangleShapeFunction[3] * segmentShapeFunctionDer[0];
            out[10] = triangleShapeFunction[4] * segmentShapeFunctionDer[0];
            out[11] = triangleShapeFunction[5] * segmentShapeFunctionDer[0];
            out[12] = triangleShapeFunction[3] * segmentShapeFunctionDer[2];
            out[13] = triangleShapeFunction[4] * segmentShapeFunctionDer[2];
            out[14] = triangleShapeFunction[5] * segmentShapeFunctionDer[2];
            out[15] = triangleShapeFunction[3] * segmentShapeFunctionDer[1];
            out[16] = triangleShapeFunction[4] * segmentShapeFunctionDer[1];
            out[17] = triangleShapeFunction[5] * segmentShapeFunctionDer[1];
            break;
          }
          default:
              DUNE_THROW(RangeError, "Component out of range.");
          }
        } else {
          DUNE_THROW(NotImplemented, "Desired derivative order is not implemented");
        }

        return;
      }

      DUNE_THROW(NotImplemented, "LagrangePrismLocalBasis::partial not implemented for order " << OrderTraits::order());
    }

  };

  /** \brief Associations of the Lagrange degrees of freedom to subentities of the reference prism
   *
   * \tparam compileTimeOrder Polynomial order of the Lagrange space in one direction (or -1 for dynamic order)
   */
  template<int compileTimeOrder>
  class LagrangePrismLocalCoefficients
    : public LagrangePrismOrderTraits<compileTimeOrder>
  {
    using OrderTraits = LagrangePrismOrderTraits<compileTimeOrder>;
  public:

    using OrderTraits::order;
    using OrderTraits::size;

    //! \brief Default constructor
    LagrangePrismLocalCoefficients (OrderTraits orderTraits)
      : OrderTraits(orderTraits)
      , localKeys_(size())
    {
      if (order()==0)
      {
        localKeys_[0] = LocalKey(0,0,0);
        return;
      }

      if (order()==1)
      {
        for (std::size_t i=0; i<size(); i++)
          localKeys_[i] = LocalKey(i,3,0);
        return;
      }

      if (order()==2)
      {
        // Vertex shape functions
        localKeys_[0] = LocalKey(0,3,0);
        localKeys_[1] = LocalKey(1,3,0);
        localKeys_[2] = LocalKey(2,3,0);
        localKeys_[3] = LocalKey(3,3,0);
        localKeys_[4] = LocalKey(4,3,0);
        localKeys_[5] = LocalKey(5,3,0);

        // Edge shape functions
        localKeys_[6] = LocalKey(0,2,0);
        localKeys_[7] = LocalKey(1,2,0);
        localKeys_[8] = LocalKey(2,2,0);
        localKeys_[9] = LocalKey(3,2,0);
        localKeys_[10] = LocalKey(4,2,0);
        localKeys_[11] = LocalKey(5,2,0);
        localKeys_[12] = LocalKey(6,2,0);
        localKeys_[13] = LocalKey(7,2,0);
        localKeys_[14] = LocalKey(8,2,0);

        // Quadrilateral sides shape functions
        localKeys_[15] = LocalKey(0,1,0);
        localKeys_[16] = LocalKey(1,1,0);
        localKeys_[17] = LocalKey(2,1,0);

        return;
      }

      // Now: the general case
      DUNE_THROW(NotImplemented, "LagrangePrismLocalCoefficients not implemented for order " << order());

    }

    LagrangePrismLocalCoefficients ()
    requires(OrderTraits::is_static_order)
      : LagrangePrismLocalCoefficients(OrderTraits())
    {}


    //! get i-th index
    const LocalKey& localKey (std::size_t i) const
    {
      return localKeys_[i];
    }

  private:
    std::vector<LocalKey> localKeys_;
  };

  /** \brief Evaluate the degrees of freedom of a Lagrange basis
   *
   * \tparam compileTimeOrder Polynomial order of the Lagrange space in one direction (or -1 for dynamic order)
   */
  template<class D, class R, int compileTimeOrder>
  class LagrangePrismLocalInterpolation
    : public LagrangePrismOrderTraits<compileTimeOrder>
  {
    using Traits = typename LagrangePrismLocalBasis<D,R,compileTimeOrder>::Traits;
    using OrderTraits = LagrangePrismOrderTraits<compileTimeOrder>;
  public:

    constexpr LagrangePrismLocalInterpolation(OrderTraits orderTraits)
      : OrderTraits(orderTraits)
    {}

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
      constexpr auto dim = Traits::dimDomain;
      auto k = OrderTraits::order();
      using Domain = typename Traits::DomainType;
      using DomainField = typename Traits::DomainFieldType;

      out.resize(OrderTraits::size());

      // Specialization for zero-order case
      if (k==0)
      {
        auto center = ReferenceElements<DomainField,dim>::general(GeometryTypes::prism).position(0,0);
        out[0] = f(center);
        return;
      }

      // Specialization for first-order case
      if (k==1)
      {
        for (unsigned int i=0; i<OrderTraits::size(); i++)
        {
          auto vertex = ReferenceElements<DomainField,3>::general(GeometryTypes::prism).position(i,3);
          out[i] = f(vertex);
        }
        return;
      }

      if (k==2)
      {
        out[0]  = f( Domain( {0.0, 0.0, 0.0} ) );
        out[1]  = f( Domain( {1.0, 0.0, 0.0} ) );
        out[2]  = f( Domain( {0.0, 1.0, 0.0} ) );
        out[3]  = f( Domain( {0.0, 0.0, 1.0} ) );
        out[4]  = f( Domain( {1.0, 0.0, 1.0} ) );
        out[5]  = f( Domain( {0.0, 1.0, 1.0} ) );
        out[6]  = f( Domain( {0.0, 0.0, 0.5} ) );
        out[7]  = f( Domain( {1.0, 0.0, 0.5} ) );
        out[8]  = f( Domain( {0.0, 1.0, 0.5} ) );
        out[9]  = f( Domain( {0.5, 0.0, 0.0} ) );
        out[10] = f( Domain( {0.0, 0.5, 0.0} ) );
        out[11] = f( Domain( {0.5, 0.5, 0.0} ) );
        out[12] = f( Domain( {0.5, 0.0, 1.0} ) );
        out[13] = f( Domain( {0.0, 0.5, 1.0} ) );
        out[14] = f( Domain( {0.5, 0.5, 1.0} ) );
        out[15] = f( Domain( {0.5, 0.0, 0.5} ) );
        out[16] = f( Domain( {0.0, 0.5, 0.5} ) );
        out[17] = f( Domain( {0.5, 0.5, 0.5} ) );

        return;
      }

      DUNE_THROW(NotImplemented, "LagrangePrismLocalInterpolation not implemented for order " << k);
    }

  };

} }    // namespace Dune::Impl

namespace Dune
{
  /** \brief Lagrange finite element for 3d prisms with arbitrary compile-time polynomial order
   *
   * \tparam D Type used for domain coordinates
   * \tparam R Type used for function values
   * \tparam compileTimeOrder Polynomial order of the Lagrange space in one direction (or -1 for dynamic order)
   */
  template<class D, class R, int compileTimeOrder = -1>
  class LagrangePrismLocalFiniteElement
    : public Impl::LagrangePrismOrderTraits<compileTimeOrder>
  {
    using OrderTraits = Impl::LagrangePrismOrderTraits<compileTimeOrder>;
  public:

    /** \brief Export number types, dimensions, etc.
     */
    using Traits = LocalFiniteElementTraits<Impl::LagrangePrismLocalBasis<D,R,compileTimeOrder>,
                                            Impl::LagrangePrismLocalCoefficients<compileTimeOrder>,
                                            Impl::LagrangePrismLocalInterpolation<D,R,compileTimeOrder> >;

    //! \brief Constructor for compile-time order
    constexpr LagrangePrismLocalFiniteElement()
      : OrderTraits()
      , basis_(*this)
      , coefficients_(*this)
      , interpolation_(*this)
    {
      static_assert(OrderTraits::is_static_order, "Default constructor only allowed for compile-time >= 0");
      static_assert((0 <= compileTimeOrder) and (compileTimeOrder <= 2), "LagrangePrism: Only order 0,1, and 2 are supported");
    }

    //! \brief Constructor for run-time order
    explicit constexpr LagrangePrismLocalFiniteElement(int runTimeOrder)
      : OrderTraits(runTimeOrder)
      , basis_(*this)
      , coefficients_(*this)
      , interpolation_(*this)
    {
      static_assert(not OrderTraits::is_static_order, "Passing a run-time order only allowed for compile-time = -1");
      if ((runTimeOrder < 0) or (runTimeOrder > 2))
        DUNE_THROW(Dune::InvalidStateException, "LagrangePrism: Only run-time order 0,1, and 2 are supported");
    }

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

    /** \brief The reference element that the local finite element is defined on
     */
    static constexpr GeometryType type ()
    {
      return GeometryTypes::prism;
    }

  private:
    Impl::LagrangePrismLocalBasis<D,R,compileTimeOrder> basis_;
    Impl::LagrangePrismLocalCoefficients<compileTimeOrder> coefficients_;
    Impl::LagrangePrismLocalInterpolation<D,R,compileTimeOrder> interpolation_;
  };

}        // namespace Dune

#endif   // DUNE_LOCALFUNCTIONS_LAGRANGE_LAGRANGEPRISM_HH
