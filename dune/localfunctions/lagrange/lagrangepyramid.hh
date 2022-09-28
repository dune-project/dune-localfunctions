// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_LAGRANGE_LAGRANGEPYRAMID_HH
#define DUNE_LOCALFUNCTIONS_LAGRANGE_LAGRANGEPYRAMID_HH

#include <array>
#include <numeric>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/math.hh>

#include <dune/geometry/referenceelements.hh>

#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localinterpolation.hh>
#include <dune/localfunctions/common/localkey.hh>

namespace Dune { namespace Impl
{
   /** \brief Lagrange shape functions of arbitrary order on the three-dimensional reference pyramid

     Lagrange shape functions of arbitrary order have the property that
     \f$\hat\phi^i(x_j) = \delta_{i,j}\f$ for certain points \f$x_j\f$.

     \tparam D Type to represent the field in the domain
     \tparam R Type to represent the field in the range
     \tparam k Polynomial order
   */
  template<class D, class R, unsigned int k>
  class LagrangePyramidLocalBasis
  {
  public:
    using Traits = LocalBasisTraits<D,3,FieldVector<D,3>,R,1,FieldVector<R,1>,FieldMatrix<R,1,3> >;

    /** \brief Number of shape functions
     */
    static constexpr std::size_t size ()
    {
      std::size_t result = 0;
      for (unsigned int i=0; i<=k; i++)
        result += power(i+1,2);
      return result;
    }

    //! \brief Evaluate all shape functions
    void evaluateFunction(const typename Traits::DomainType& in,
                          std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(size());

      // Specialization for zero-order case
      if (k==0)
      {
        out[0] = 1;
        return;
      }

      if (k==1)
      {
        if(in[0] > in[1])
        {
          out[0] = (1-in[0])*(1-in[1])-in[2]*(1-in[1]);
          out[1] = in[0]*(1-in[1])-in[2]*in[1];
          out[2] = (1-in[0])*in[1]-in[2]*in[1];
          out[3] = in[0]*in[1]+in[2]*in[1];
        }
        else
        {
          out[0] = (1-in[0])*(1-in[1])-in[2]*(1-in[0]);
          out[1] = in[0]*(1-in[1])-in[2]*in[0];
          out[2] = (1-in[0])*in[1]-in[2]*in[0];
          out[3] = in[0]*in[1]+in[2]*in[0];
        }

        out[4] = in[2];

        return;
      }

      if (k==2)
      {
        // transform to reference element with base [-1,1]^2
        const R x = 2.0*in[0] + in[2] - 1.0;
        const R y = 2.0*in[1] + in[2] - 1.0;
        const R z = in[2];

        if (x > y)
        {
          // vertices
          out[0] = 0.25*(x + z)*(x + z - 1)*(y - z - 1)*(y - z);
          out[1] = -0.25*(x + z)*(y - z)*((x + z + 1)*(-y + z + 1) - 4*z) - z*(x - y);
          out[2] = 0.25*(x + z)*(y - z)*(y - z + 1)*(x + z - 1);
          out[3] = 0.25*(y - z)*(x + z)*(y - z + 1)*(x + z + 1);
          out[4] = z*(2*z - 1);

          // lower edges
          out[5] = -0.5*(y - z + 1)*(x + z - 1)*(y - 1)*x;
          out[6] = -0.5*(y - z + 1)*(((x + z + 1)*(y - 1)*x - z) + z*(2*y + 1));
          out[7] = -0.5*(x + z - 1)*(((y - z - 1)*(x + 1)*y - z) + z*(2*x + 1));
          out[8] = -0.5*(y - z + 1)*(x + z - 1)*(x + 1)*y;

          // upper edges
          out[9] = z*(x + z - 1)*(y - z - 1);
          out[10] = -z*((x + z + 1)*(y - z - 1) + 4*z);
          out[11] = -z*(y - z + 1)*(x + z - 1);
          out[12] = z*(y - z + 1)*(x + z + 1);

          // base face
          out[13] = (y - z + 1)*(x + z - 1)*((y - 1)*(x + 1) + z*(x - y + z + 1));
        }
        else
        {
          // vertices
          out[0] = 0.25*(y + z)*(y + z - 1)*(x - z - 1)*(x - z);
          out[1] = -0.25*(x - z)*(y + z)*(x - z + 1)*(-y - z + 1);
          out[2] = 0.25*(x - z)*(y + z)*((x - z - 1)*(y + z + 1) + 4*z) + z*(x - y);
          out[3] = 0.25*(y + z)*(x - z)*(x - z + 1)*(y + z + 1);
          out[4] = z*(2*z - 1);

          // lower edges
          out[5] = -0.5*(y + z - 1)*(((x - z - 1)*(y + 1)*x - z) + z*(2*y + 1));
          out[6] = -0.5*(x - z + 1)*(y + z - 1)*(y + 1)*x;
          out[7] = -0.5*(x - z + 1)*(y + z - 1)*(x - 1)*y;
          out[8] = -0.5*(x - z + 1)*(((y + z + 1)*(x - 1)*y - z) + z*(2*x + 1));

          // upper edges
          out[9] = z*(y + z - 1)*(x - z - 1);
          out[10] = -z*(x - z + 1)*(y + z - 1);
          out[11] = -z*((y + z + 1)*(x - z - 1) + 4*z);
          out[12] = z*(x - z + 1)*(y + z + 1);

          // base face
          out[13] = (x - z + 1)*(y + z - 1)*((y + 1)*(x - 1) - z*(x - y - z - 1));
        }

        return;
      }

      DUNE_THROW(NotImplemented, "LagrangePyramidLocalBasis::evaluateFunction for order " << k);
    }

    /** \brief Evaluate Jacobian of all shape functions
     *
     * \param x Point in the reference pyramid where to evaluation the Jacobians
     * \param[out] out The Jacobians of all shape functions at the point x
     */
    void evaluateJacobian(const typename Traits::DomainType& in,
                          std::vector<typename Traits::JacobianType>& out) const
    {
      out.resize(size());

      // Specialization for k==0
      if (k==0)
      {
        std::fill(out[0][0].begin(), out[0][0].end(), 0);
        return;
      }

      if (k==1)
      {
        if(in[0] > in[1])
        {
          out[0][0] = {-1 + in[1], -1 + in[0] + in[2], -1 + in[1]};
          out[1][0] = { 1 - in[1],     -in[0] - in[2],     -in[1]};
          out[2][0] = {    -in[1],  1 - in[0] - in[2],     -in[1]};
          out[3][0] = {     in[1],      in[0] + in[2],      in[1]};
        }
        else
        {
          out[0][0] = {-1 + in[1] + in[2], -1 + in[0], -1 + in[0]};
          out[1][0] = { 1 - in[1] - in[2],     -in[0],     -in[0]};
          out[2][0] = {    -in[1] - in[2],  1 - in[0],     -in[0]};
          out[3][0] = {     in[1] + in[2],      in[0],      in[0]};
        }

        out[4][0] = {0, 0, 1};
        return;
      }

      if (k==2)
      {
        // transform to reference element with base [-1,1]^2
        const R x = 2.0*in[0] + in[2] - 1.0;
        const R y = 2.0*in[1] + in[2] - 1.0;
        const R z = in[2];

        // transformation of the gradient leads to a multiplication
        // with the Jacobian [2 0 0; 0 2 0; 1 1 1]
        if (x > y)
        {
          // vertices
          out[0][0][0] = 0.5*(y - z - 1)*(y - z)*(2*x + 2*z - 1);
          out[0][0][1] = 0.5*(x + z)*(x + z - 1)*(2*y - 2*z - 1);
          out[0][0][2] = 0.5*(out[0][0][0] + out[0][0][1])
                         + 0.25*((2*x + 2*z - 1)*(y - z - 1)*(y - z)
                                 + (x + z)*(x + z - 1)*(-2*y + 2*z + 1));

          out[1][0][0] = 2*(-0.25*((y - z)*((x + z + 1)*(-y + z + 1) - 4*z)
                                   + (x + z)*(y - z)*(-y + z + 1)) - z);
          out[1][0][1] = 2*(-0.25*((x + z)*((x + z + 1)*(-y + z + 1) - 4*z)
                                   + (x + z)*(y - z)*(-(x + z + 1))) + z);
          out[1][0][2] = 0.5*(out[1][0][0] + out[1][0][1])
                         - 0.25*((y - z)*((x + z + 1)*(-y + z + 1) - 4*z)
                                 - (x + z)*((x + z + 1)*(-y + z + 1) - 4*z)
                                 + (x + z)*(y - z)*(x - y + 2*z - 2))
                         - (x - y);

          out[2][0][0] = 0.5*(y - z)*(y - z + 1)*(2*x + 2*z - 1);
          out[2][0][1] = 0.5*(x + z)*(2*y - 2*z + 1)*(x + z - 1);
          out[2][0][2] = 0.5*(out[2][0][0] + out[2][0][1])
                         + 0.25*((y - x - 2*z)*(y - z + 1)*(x + z - 1)
                                 + (x + z)*(y - z)*(y - x - 2*z + 2));

          out[3][0][0] = 0.5*(y - z)*(2*x + 2*z + 1)*(y - z + 1);
          out[3][0][1] = 0.5*(2*y - 2*z + 1)*(x + z)*(x + z + 1);
          out[3][0][2] = 0.5*(out[3][0][0] + out[3][0][1])
                         + 0.25*((y - x - 2*z)*(y - z + 1)*(x + z + 1)
                                 + (y - z)*(x + z)*(y - x - 2*z));

          out[4][0][0] = 0;
          out[4][0][1] = 0;
          out[4][0][2] = 4*z - 1;

          // lower edges
          out[5][0][0] = -(y - z + 1)*(y - 1)*(2*x + z - 1);
          out[5][0][1] = -(x + z - 1)*(y - 1)*x - (y - z + 1)*(x + z - 1)*x;
          out[5][0][2] = 0.5*(out[5][0][0] + out[5][0][1])
                         + 0.5*(x + z - 1)*(y - 1)*x - 0.5*(y - z + 1)*(y - 1)*x;

          out[6][0][0] = -(y - z + 1)*(2*x + z + 1)*(y - 1);
          out[6][0][1] = -(((x + z + 1)*(y - 1)*x - z) + z*(2*y + 1)
                           + (y - z + 1)*((x + z + 1)*x + 2*z));
          out[6][0][2] = 0.5*(out[6][0][0] + out[6][0][1])
                         - 0.5*(-(((x + z + 1)*(y - 1)*x - z) + z*(2*y + 1))
                                + (y - z + 1)*(((y - 1)*x - 1) + 2*y + 1));

          out[7][0][0] = -(((y - z - 1)*(x + 1)*y - z) + z*(2*x + 1)
                           + (x + z - 1)*((y - z - 1)*y + 2*z));
          out[7][0][1] = -(x + z - 1)*(2*y - z - 1)*(x + 1);
          out[7][0][2] = 0.5*(out[7][0][0] + out[7][0][1])
                         - 0.5*(((y - z - 1)*(x + 1)*y - z) + z*(2*x + 1)
                                + (x + z - 1)*((-(x + 1)*y - 1) + 2*x + 1));

          out[8][0][0] = -(y - z + 1)*(2*x + z)*y;
          out[8][0][1] = -(2*y - z + 1)*(x + z - 1)*(x + 1);
          out[8][0][2] = 0.5*(out[8][0][0] + out[8][0][1])
                         - 0.5*(-x + y - 2*z + 2)*(x + 1)*y;

          // upper edges
          out[9][0][0] = 2*z*(y - z - 1);
          out[9][0][1] = 2*z*(x + z - 1);
          out[9][0][2] = 0.5*(out[9][0][0] + out[9][0][1])
                         + (x + z - 1)*(y - z - 1) + z*(-x + y - 2*z);

          out[10][0][0] = -2*z*(y - z - 1);
          out[10][0][1] = -2*z*(x + z + 1);
          out[10][0][2] = 0.5*(out[10][0][0] + out[10][0][1])
                          - ((x + z + 1)*(y - z - 1) + 4*z)
                          - z*(-x + y - 2*z + 2);

          out[11][0][0] = -2*z*(y - z + 1);
          out[11][0][1] = -2*z*(x + z - 1);
          out[11][0][2] = 0.5*(out[11][0][0] + out[11][0][1])
                          - (y - z + 1)*(x + z - 1) - z*(-x + y - 2*z + 2);

          out[12][0][0] = 2*z*(y - z + 1);
          out[12][0][1] = 2*z*(x + z + 1);
          out[12][0][2] = 0.5*(out[12][0][0] + out[12][0][1])
                          + (y - z + 1)*(x + z + 1) + z*(-x + y - 2*z);

          // base face
          out[13][0][0] = 2*((y - z + 1)*((y - 1)*(x + 1) + z*(x - y + z + 1))
                             + (y - z + 1)*(x + z - 1)*(y - 1 + z));
          out[13][0][1] = 2*((x + z - 1)*((y - 1)*(x + 1) + z*(x - y + z + 1))
                             + (y - z + 1)*(x + z - 1)*(x + 1 - z));
          out[13][0][2] = 0.5*(out[13][0][0] + out[13][0][1])
                          + ((-x + y - 2*z + 2)*((y - 1)*(x + 1) + z*(x - y + z + 1))
                             + (y - z + 1)*(x + z - 1)*(x - y + 2*z + 1));
        }
        else
        {
          // vertices
          out[0][0][0] = 0.5*(y + z)*(y + z - 1)*(2*x - 2*z - 1);
          out[0][0][1] = 0.5*(2*y + 2*z - 1)*(x - z - 1)*(x - z);
          out[0][0][2] = 0.5*(out[0][0][0] + out[0][0][1])
                         + 0.25*((2*y + 2*z - 1)*(x - z - 1)*(x - z)
                                 + (y + z)*(y + z - 1)*(-2*x + 2*z + 1));

          out[1][0][0] = -0.5*(y + z)*(2*x - 2*z + 1)*(-y - z + 1);
          out[1][0][1] = -0.5*(x - z)*(x - z + 1)*(-2*y - 2*z + 1);
          out[1][0][2] = 0.5*(out[1][0][0] + out[1][0][1])
                         - 0.25*((x - y - 2*z)*(x - z + 1)*(-y - z + 1)
                                 + (x - z)*(y + z)*(-x + y + 2*z - 2));

          out[2][0][0] = 0.5*((y + z)*((x - z - 1)*(y + z + 1) + 4*z)
                              + (x - z)*(y + z)*(y + z + 1) + 4*z);
          out[2][0][1] = 0.5*((x - z)*((x - z - 1)*(y + z + 1) + 4*z)
                              + (x - z)*(y + z)*(x - z - 1) - 4*z);
          out[2][0][2] = 0.5*(out[2][0][0] + out[2][0][1])
                         + 0.25*((x - y - 2*z)*((x - z - 1)*(y + z + 1) + 4*z)
                                 + (x - z)*(y + z)*(x - y - 2*z + 2) + 4*(x - y));

          out[3][0][0] = 0.5*(y + z)*(2*x - 2*z + 1)*(y + z + 1);
          out[3][0][1] = 0.5*(x - z)*(x - z + 1)*(2*y + 2*z + 1);
          out[3][0][2] = 0.5*(out[3][0][0] + out[3][0][1])
                         + 0.25*((x - y - 2*z)*(x - z + 1)*(y + z + 1)
                                 + (y + z)*(x - z)*(x - y - 2*z));

          out[4][0][0] = 0;
          out[4][0][1] = 0;
          out[4][0][2] = 4*z - 1;

          // lower edges
          out[5][0][0] = -(y + z - 1)*(2*x - z - 1)*(y + 1);
          out[5][0][1] = -(((x - z - 1)*(y + 1)*x - z) + z*(2*y + 1)
                           + (y + z - 1)*((x - z - 1)*x + 2*z));
          out[5][0][2] = 0.5*(out[5][0][0] + out[5][0][1])
                         - 0.5*((((x - z - 1)*(y + 1)*x - z) + z*(2*y + 1))
                                + (y + z - 1)*((-(y + 1)*x - 1) + 2*y + 1));

          out[6][0][0] = -(2*x - z + 1)*(y + z - 1)*(y + 1);
          out[6][0][1] = -(x - z + 1)*(2*y + z)*x;
          out[6][0][2] = 0.5*(out[6][0][0] + out[6][0][1])
                         - 0.5*(x - y - 2*z + 2)*(y + 1)*x;

          out[7][0][0] = -(2*x - z)*(y + z - 1)*y;
          out[7][0][1] = -(x - z + 1)*(2*y + z - 1)*(x - 1);
          out[7][0][2] = 0.5*(out[7][0][0] + out[7][0][1])
                         - 0.5*(x - y - 2*z + 2)*(x - 1)*y;

          out[8][0][0] = -(((y + z + 1)*(x - 1)*y - z) + z*(2*x + 1)
                           + (x - z + 1)*((y + z + 1)*y + 2*z));
          out[8][0][1] = -(x - z + 1)*(2*y + z + 1)*(x - 1);
          out[8][0][2] = 0.5*(out[8][0][0] + out[8][0][1])
                         - 0.5*(-(((y + z + 1)*(x - 1)*y - z) + z*(2*x + 1))
                                + (x - z + 1)*(((x - 1)*y - 1) + 2*x + 1));

          // upper edges
          out[9][0][0] = 2*z*(y + z - 1);
          out[9][0][1] = 2*z*(x - z - 1);
          out[9][0][2] = 0.5*(out[9][0][0] + out[9][0][1])
                         + (y + z - 1)*(x - z - 1) + z*(x - y - 2*z);

          out[10][0][0] = -2*z*(y + z - 1);
          out[10][0][1] = -2*z*(x - z + 1);
          out[10][0][2] = 0.5*(out[10][0][0] + out[10][0][1])
                          - (x - z + 1)*(y + z - 1) - z*(x - y - 2*z + 2);

          out[11][0][0] = -2*z*(y + z + 1);
          out[11][0][1] = -2*z*(x - z - 1);
          out[11][0][2] = 0.5*(out[11][0][0] + out[11][0][1])
                          - ((y + z + 1)*(x - z - 1) + 4*z) - z*(x - y - 2*z + 2);

          out[12][0][0] = 2*z*(y + z + 1);
          out[12][0][1] = 2*z*(x - z + 1);
          out[12][0][2] = 0.5*(out[12][0][0] + out[12][0][1])
                          + (x - z + 1)*(y + z + 1) + z*(x - y - 2*z);

          // base face
          out[13][0][0] = 2*((y + z - 1)*((y + 1)*(x - 1) - z*(x - y - z - 1))
                             + (x - z + 1)*(y + z - 1)*(y + 1 - z));
          out[13][0][1] = 2*((x - z + 1)*((y + 1)*(x - 1) - z*(x - y - z - 1))
                             + (x - z + 1)*(y + z - 1)*(x - 1 + z));
          out[13][0][2] = 0.5*(out[13][0][0] + out[13][0][1])
                          + (x - y - 2*z + 2)*((y + 1)*(x - 1) - z*(x - y - z - 1))
                          + (x - z + 1)*(y + z - 1)*(-(x - y - 2*z - 1));
        }

        return;
      }

      DUNE_THROW(NotImplemented, "LagrangePyramidLocalBasis::evaluateJacobian for order " << k);
    }

    /** \brief Evaluate partial derivatives of any order of all shape functions
     *
     * \param order Order of the partial derivatives, in the classic multi-index notation
     * \param in Position where to evaluate the derivatives
     * \param[out] out The desired partial derivatives
     */
    void partial(const std::array<unsigned int,3>& order,
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

      if (k==0)
      {
        out[0] = 0;
        return;
      }

      if (k==1)
      {
        if (totalOrder == 1)
        {
          auto const direction = std::distance(order.begin(), std::find(order.begin(), order.end(), 1));
          if (in[0] > in[1])
          {
            switch (direction)
            {
              case 0:
                out[0] = -1 + in[1];
                out[1] = 1  - in[1];
                out[2] = -in[1];
                out[3] = in[1];
                out[4] = 0;
                break;
              case 1:
                out[0] = -1 + in[0] + in[2];
                out[1] = -in[0] - in[2];
                out[2] = 1 - in[0] - in[2];
                out[3] = in[0]+in[2];
                out[4] = 0;
                break;
              case 2:
                out[0] = -1 + in[1];
                out[1] = -in[1];
                out[2] = -in[1];
                out[3] = in[1];
                out[4] = 1;
                break;
              default:
                DUNE_THROW(RangeError, "Component out of range.");
            }
          }
          else /* (in[0] <= in[1]) */
          {
            switch (direction)
            {
              case 0:
                out[0] = -1 + in[1] + in[2];
                out[1] = 1 - in[1] - in[2];
                out[2] = -in[1] - in[2];
                out[3] = in[1] + in[2];
                out[4] = 0;
                break;
              case 1:
                out[0] = -1 + in[0];
                out[1] = -in[0];
                out[2] = 1 - in[0];
                out[3] = in[0];
                out[4] = 0;
                break;
              case 2:
                out[0] = -1 + in[0];
                out[1] = -in[0];
                out[2] = -in[0];
                out[3] = in[0];
                out[4] = 1;
                break;
              default:
                DUNE_THROW(RangeError, "Component out of range.");
            }
          }
        } else if (totalOrder == 2)
        {
          if ((order[0] == 1 && order[1] == 1) ||
              (order[1] == 1 && order[2] == 1 && in[0] > in[1]) ||
              (order[0] == 1 && order[2] == 1 && in[0] <=in[1]))
          {
            out = {1, -1, -1, 1, 0};
          } else
          {
            out = {0, 0, 0, 0, 0};
          }

        } else
        {
          out = {0, 0, 0, 0, 0};
        }

        return;
      }

      if (k==2)
      {
        if (totalOrder == 1)
        {
          // transform to reference element with base [-1,1]^2
          const R x = 2.0*in[0] + in[2] - 1.0;
          const R y = 2.0*in[1] + in[2] - 1.0;
          const R z = in[2];

          auto const direction = std::distance(order.begin(), std::find(order.begin(), order.end(), 1));

          // transformation of the gradient leads to a multiplication
          // with the Jacobian [2 0 0; 0 2 0; 1 1 1]
          if (x > y)
          {
            switch (direction)
            {
            case 0:
              out[0] = 0.5*(y - z - 1)*(y - z)*(2*x + 2*z - 1);
              out[1] = 2*(-0.25*((y - z)*((x + z + 1)*(-y + z + 1) - 4*z) + (x + z)*(y - z)*(-y + z + 1)) - z);
              out[2] = 0.5*(y - z)*(y - z + 1)*(2*x + 2*z - 1);
              out[3] = 0.5*(y - z)*(2*x + 2*z + 1)*(y - z + 1);
              out[4] = 0;
              out[5] = -(y - z + 1)*(2*x + z - 1)*(y - 1);
              out[6] = -(y - z + 1)*(2*x + z + 1)*(y - 1);
              out[7] = -(((y - z - 1)*(x + 1)*y - z) + z*(2*x + 1) + (x + z - 1)*((y - z - 1)*y + 2*z));
              out[8] = -(y - z + 1)*(2*x + z)*y;
              out[9] = 2*z*(y - z - 1);
              out[10] = -2*z*(y - z - 1);
              out[11] = -2*z*(y - z + 1);
              out[12] = 2*z*(y - z + 1);
              out[13] = 2*((y - z + 1)*((y - 1)*(x + 1) + z*(x - y + z + 1)) + (y - z + 1)*(x + z - 1)*(y - 1 + z));
              break;
            case 1:
              out[0] = 0.5*(x + z)*(x + z - 1)*(2*y - 2*z - 1);
              out[1] = 2*(-0.25*((x + z)*((x + z + 1)*(-y + z + 1) - 4*z) + (x + z)*(y - z)*(-(x + z + 1))) + z);
              out[2] = 0.5*(x + z)*(2*y - 2*z + 1)*(x + z - 1);
              out[3] = 0.5*(2*y - 2*z + 1)*(x + z)*(x + z + 1);
              out[4] = 0;
              out[5] = -(x + z - 1)*(y - 1)*x - (y - z + 1)*(x + z - 1)*x;
              out[6] = -(((x + z + 1)*(y - 1)*x - z) + z*(2*y + 1) + (y - z + 1)*((x + z + 1)*x + 2*z));
              out[7] = -(x + z - 1)*(2*y - z - 1)*(x + 1);
              out[8] = -(2*y - z + 1)*(x + z - 1)*(x + 1);
              out[9] = 2*z*(x + z - 1);
              out[10] = -2*z*(x + z + 1);
              out[11] = -2*z*(x + z - 1);
              out[12] = 2*z*(x + z + 1);
              out[13] = 2*((x + z - 1)*((y - 1)*(x + 1) + z*(x - y + z + 1)) + (y - z + 1)*(x + z - 1)*(x + 1 - z));
              break;
            case 2:
              out[0] = -((y - z)*(2*x + 2*z - 1)*(z - y + 1))/2;
              out[1] = ((y - z + 1)*(y - 2*x + z + 2*x*y - 2*x*z + 2*y*z - 2*z*z))/2;
              out[2] = ((y - z)*(2*x + 2*z - 1)*(y - z + 1))/2;
              out[3] = ((y - z)*(2*x + 2*z + 1)*(y - z + 1))/2;
              out[4] = 4*z - 1;
              out[5] = (-(y - z + 1)*(2*x + z - 1)*(y - 1) - (x + z - 1)*(y - 1)*x - (y - z + 1)*(x + z - 1)*x + (x + z - 1)*(y - 1)*x - (y - z + 1)*(y - 1)*x)/2;
              out[6] = -((y - z + 1)*(3*y - 2*x + z + 3*x*y + x*z + y*z + x*x - 1))/2;
              out[7] = z - z*(2*x + 1) - ((2*z - y*(z - y + 1))*(x + z - 1))/2 - ((2*x - y*(x + 1))*(x + z - 1))/2 + ((x + 1)*(x + z - 1)*(z - 2*y + 1))/2 + y*(x + 1)*(z - y + 1);
              out[8] = -((y - z + 1)*(y + z + 3*x*y + x*z + y*z + x*x - 1))/2;
              out[9] = -(x + 3*z - 1)*(z - y + 1);
              out[10] = (x + z + 1)*(z - y + 1) - 2*y*z - 6*z + 2*z*z;
              out[11] = -(x + 3*z - 1)*(y - z + 1);
              out[12] = (x + 3*z + 1)*(y - z + 1);
              out[13] = (y - z + 1)*(2*y - 3*x + z + 2*x*y + 6*x*z - 2*y*z + 2*x*x + 4*z*z - 3);
              break;
            default:
              DUNE_THROW(RangeError, "Component out of range.");
            }
          }
          else // x <= y
          {
            switch (direction)
            {
            case 0:
              out[0] = -((y + z)*(2*z - 2*x + 1)*(y + z - 1))/2;
              out[1] = ((y + z)*(2*x - 2*z + 1)*(y + z - 1))/2;
              out[2] = -((y + z + 1)*(y - 3*z - 2*x*y - 2*x*z + 2*y*z + 2*z*z))/2;
              out[3] = ((y + z)*(2*x - 2*z + 1)*(y + z + 1))/2;
              out[4] = 0;
              out[5] = (y + 1)*(y + z - 1)*(z - 2*x + 1);
              out[6] = -(y + 1)*(2*x - z + 1)*(y + z - 1);
              out[7] = -y*(2*x - z)*(y + z - 1);
              out[8] = z - z*(2*x + 1) - (2*z + y*(y + z + 1))*(x - z + 1) - y*(x - 1)*(y + z + 1);
              out[9] = 2*z*(y + z - 1);
              out[10] = -2*z*(y + z - 1);
              out[11] = -2*z*(y + z + 1);
              out[12] = 2*z*(y + z + 1);
              out[13] = 2*(y + z - 1)*(2*x - z + 2*x*y - 2*x*z + 2*z*z);
              break;
            case 1:
              out[0] = -(x - z)*(y + z - 0.5)*(z - x + 1);
              out[1] = ((x - z)*(2*y + 2*z - 1)*(x - z + 1))/2;
              out[2] = -((z - x + 1)*(x + 3*z + 2*x*y + 2*x*z - 2*y*z - 2*z*z))/2;
              out[3] = ((x - z)*(2*y + 2*z + 1)*(x - z + 1))/2;
              out[4] = 0;
              out[5] = z - z*(2*y + 1) - (2*z - x*(z - x + 1))*(y + z - 1) + x*(y + 1)*(z - x + 1);
              out[6] = -x*(2*y + z)*(x - z + 1);
              out[7] = -(x - 1)*(x - z + 1)*(2*y + z - 1);
              out[8] = -(x - 1)*(x - z + 1)*(2*y + z + 1);
              out[9] = -2*z*(z - x + 1);
              out[10] = -2*z*(x - z + 1);
              out[11] = 2*z*(z - x + 1);
              out[12] = 2*z*(x - z + 1);
              out[13] = 2*(x - z + 1)*(2*x*y - z - 2*y + 2*y*z + 2*z*z);
              break;
            case 2:
              out[0] = -((x - z)*(2*y + 2*z - 1)*(z - x + 1))/2;
              out[1] = ((x - z)*(2*y + 2*z - 1)*(x - z + 1))/2;
              out[2] = ((x - z + 1)*(x - 2*y + z + 2*x*y + 2*x*z - 2*y*z - 2*z*z))/2;
              out[3] = ((x - z)*(2*y + 2*z + 1)*(x - z + 1))/2;
              out[4] = 4*z - 1;
              out[5] = z - z*(2*y + 1) - ((2*z - x*(z - x + 1))*(y + z - 1))/2 - ((2*y - x*(y + 1))*(y + z - 1))/2 + ((y + 1)*(y + z - 1)*(z - 2*x + 1))/2 + x*(y + 1)*(z - x + 1);
              out[6] = -((x - z + 1)*(x + z + 3*x*y + x*z + y*z + y*y - 1))/2;
              out[7] = -((x - z + 1)*(3*x*y - 4*y - z - x + x*z + y*z + y*y + 1))/2;
              out[8] = -((x - z + 1)*(3*x - 2*y + z + 3*x*y + x*z + y*z + y*y - 1))/2;
              out[9] = -(z - x + 1)*(y + 3*z - 1);
              out[10] = -(x - z + 1)*(y + 3*z - 1);
              out[11] = (y + z + 1)*(z - x + 1) - 2*x*z - 6*z + 2*z*z;
              out[12] = (x - z + 1)*(y + 3*z + 1);
              out[13] = (x - z + 1)*(2*x - 3*y + z + 2*x*y - 2*x*z + 6*y*z + 2*y*y + 4*z*z - 3);
              break;
            default:
              DUNE_THROW(RangeError, "Component out of range.");
            }
          }
        } else {
          DUNE_THROW(NotImplemented, "Desired derivative order is not implemented");
        }

        return;
      }

      DUNE_THROW(NotImplemented, "LagrangePyramidLocalBasis::partial for order " << k);
    }

    //! \brief Polynomial order of the shape functions
    static constexpr unsigned int order ()
    {
      return k;
    }
  };

  /** \brief Associations of the Lagrange degrees of freedom to subentities of the reference pyramid
   *
   * \tparam k Polynomial order of the Lagrange space in one direction
   */
  template<unsigned int k>
  class LagrangePyramidLocalCoefficients
  {
  public:
    //! \brief Default constructor
    LagrangePyramidLocalCoefficients ()
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
          localKeys_[i] = LocalKey(i,3,0);
        return;
      }

      if (k==2)
      {
        // Vertex shape functions
        localKeys_[0] = LocalKey(0,3,0);
        localKeys_[1] = LocalKey(1,3,0);
        localKeys_[2] = LocalKey(2,3,0);
        localKeys_[3] = LocalKey(3,3,0);
        localKeys_[4] = LocalKey(4,3,0);

        // Edge shape functions
        localKeys_[5] = LocalKey(0,2,0);
        localKeys_[6] = LocalKey(1,2,0);
        localKeys_[7] = LocalKey(2,2,0);
        localKeys_[8] = LocalKey(3,2,0);
        localKeys_[9] = LocalKey(4,2,0);
        localKeys_[10] = LocalKey(5,2,0);
        localKeys_[11] = LocalKey(6,2,0);
        localKeys_[12] = LocalKey(7,2,0);

        // base face shape function
        localKeys_[13] = LocalKey(0,1,0);

        return;
      }

      // No general case
      DUNE_THROW(NotImplemented, "LagrangePyramidLocalCoefficients for order " << k);

    }

    //! number of coefficients
    static constexpr std::size_t size ()
    {
      std::size_t result = 0;
      for (unsigned int i=0; i<=k; i++)
        result += power(i+1,2);
      return result;
    }

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
   * \tparam LocalBasis The corresponding set of shape functions
   */
  template<class LocalBasis>
  class LagrangePyramidLocalInterpolation
  {
  public:

    /** \brief Evaluate a given function at the Lagrange nodes
     *
     * \tparam F Type of function to evaluate
     * \tparam C Type used for the values of the function
     * \param[in] ff Function to evaluate
     * \param[out] out Array of function values
     */
    template<typename F, typename C>
    void interpolate (const F& ff, std::vector<C>& out) const
    {
      constexpr auto k = LocalBasis::order();
      using D = typename LocalBasis::Traits::DomainType;
      using DF = typename LocalBasis::Traits::DomainFieldType;

      auto&& f = Impl::makeFunctionWithCallOperator<D>(ff);

      out.resize(LocalBasis::size());

      // Specialization for zero-order case
      if (k==0)
      {
        auto center = ReferenceElements<DF,3>::general(GeometryTypes::pyramid).position(0,0);
        out[0] = f(center);
        return;
      }

      // Specialization for first-order case
      if (k==1)
      {
        for (unsigned int i=0; i<LocalBasis::size(); i++)
        {
          auto vertex = ReferenceElements<DF,3>::general(GeometryTypes::pyramid).position(i,3);
          out[i] = f(vertex);
        }
        return;
      }

      // Specialization for second-order case
      if (k==2)
      {
        out[0]  = f( D( {0.0, 0.0, 0.0} ) );
        out[1]  = f( D( {1.0, 0.0, 0.0} ) );
        out[2]  = f( D( {0.0, 1.0, 0.0} ) );
        out[3]  = f( D( {1.0, 1.0, 0.0} ) );
        out[4]  = f( D( {0.0, 0.0, 1.0} ) );
        out[5]  = f( D( {0.0, 0.5, 0.0} ) );
        out[6]  = f( D( {1.0, 0.5, 0.0} ) );
        out[7]  = f( D( {0.5, 0.0, 0.0} ) );
        out[8]  = f( D( {0.5, 1.0, 0.0} ) );
        out[9]  = f( D( {0.0, 0.0, 0.5} ) );
        out[10] = f( D( {0.5, 0.0, 0.5} ) );
        out[11] = f( D( {0.0, 0.5, 0.5} ) );
        out[12] = f( D( {0.5, 0.5, 0.5} ) );
        out[13] = f( D( {0.5, 0.5, 0.0} ) );

        return;
      }

      DUNE_THROW(NotImplemented, "LagrangePyramidLocalInterpolation not implemented for order " << k);
    }

  };

} }    // namespace Dune::Impl

namespace Dune
{
  /** \brief Lagrange finite element for 3d pyramids with compile-time polynomial order
   *
   * \tparam D Type used for domain coordinates
   * \tparam R Type used for shape function values
   * \tparam k Polynomial order, only orders 1 and 2 are actually implemented
   *
   * Lagrange shape functions are tricky.  In the paper mentioned below, Christian Wieners states
   * "There exists no continuously differentiable conforming shape function for the pyramid
   *  which is linear, resp. bilinear on the faces."
   * The same holds, mutatis mutandis, for second-order Lagrange functions.
   * The usual remedy, employed here, is to use shape functions that are continuous,
   * but only piecewise differentiable.
   *
   * More specifically, the implementations in this file are taken from the following papers:
   * * First order: C. Wieners, Conforming discretizations on tetrahedrons, pyramids, prisms and hexahedrons (1997)
   * * Second order: L. Liu et.al. On Higher Order Pyramidal Finite Elements, (2011) [anisotropic variant]
   *
   * These are piecewise trilinear/triquadratic basis functions with the following properties:
   * * Shape functions with Lagrange property
   * * Each basis function is bilinear/biquadratic on the rectangular face
   * * Each basis function is linear/quadratic on all triangular faces
   * * Each basis function is continuous in the interelement boundary
   *
   * As the derivatives are not continuous, numerical quadrature for expressions involving
   * shape function derivatives should employ a composite rule composed of two tetrahedral parts.
   *
   * \warning The shape functions currently do not sum up to 1, even though my understanding
   *  of the Liu et al. paper is that they should.
   *
   */
  template<class D, class R, int k>
  class LagrangePyramidLocalFiniteElement
  {
  public:
    /** \brief Export number types, dimensions, etc.
     */
    using Traits = LocalFiniteElementTraits<Impl::LagrangePyramidLocalBasis<D,R,k>,
                                            Impl::LagrangePyramidLocalCoefficients<k>,
                                            Impl::LagrangePyramidLocalInterpolation<Impl::LagrangePyramidLocalBasis<D,R,k> > >;

    /** \brief Default constructor
     *
     * \deprecated This explicit implementation only exists to work around a bug in clang 3.8
     *   which disappeared in clang 6
     */
    LagrangePyramidLocalFiniteElement() {}

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
      return Impl::LagrangePyramidLocalBasis<D,R,k>::size();
    }

    /** \brief The reference element that the local finite element is defined on
     */
    static constexpr GeometryType type ()
    {
      return GeometryTypes::pyramid;
    }

  private:
    Impl::LagrangePyramidLocalBasis<D,R,k> basis_;
    Impl::LagrangePyramidLocalCoefficients<k> coefficients_;
    Impl::LagrangePyramidLocalInterpolation<Impl::LagrangePyramidLocalBasis<D,R,k> > interpolation_;
  };

}        // namespace Dune

#endif   // DUNE_LOCALFUNCTIONS_LAGRANGE_LAGRANGEPYRAMID_HH
