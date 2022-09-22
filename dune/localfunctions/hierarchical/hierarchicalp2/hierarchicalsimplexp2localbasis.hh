// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_HIERARCHICAL_SIMPLEX_P2_LOCALBASIS_HH
#define DUNE_HIERARCHICAL_SIMPLEX_P2_LOCALBASIS_HH

/** \file
    \brief Hierarchical p2 shape functions for the simplex
 */

#include <numeric>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dune/localfunctions/common/localbasis.hh>

namespace Dune
{
  template<class D, class R, int dim>
  class HierarchicalSimplexP2LocalBasis
  {
  public:
    HierarchicalSimplexP2LocalBasis()
    {
      DUNE_THROW(Dune::NotImplemented,"HierarchicalSimplexP2LocalBasis not implemented for dim > 3.");
    }
  };

  /**@ingroup LocalBasisImplementation
     \brief Hierarchical P2 basis in 1d.

     The shape functions are associated to the following points:

     - \f$ f_0 \f$ ~ (0.0)   // linear function
     - \f$ f_1 \f$ ~ (1.0)   // linear function
     - \f$ f_2 \f$ ~ (0.5)   // quadratic bubble

     \tparam D Type to represent the field in the domain.
     \tparam R Type to represent the field in the range.

     \nosubgrouping
   */
  template<class D, class R>
  class HierarchicalSimplexP2LocalBasis<D,R,1>
  {
  public:
    //! \brief export type traits for function signature
    typedef LocalBasisTraits<D,1,Dune::FieldVector<D,1>,R,1,Dune::FieldVector<R,1>,
        Dune::FieldMatrix<R,1,1> > Traits;

    //! \brief number of shape functions
    unsigned int size () const
    {
      return 3;
    }

    //! \brief Evaluate all shape functions
    inline void evaluateFunction (const typename Traits::DomainType& in,
                                  std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(3);

      out[0] = 1-in[0];
      out[1] = 1-4*(in[0]-0.5)*(in[0]-0.5);
      out[2] = in[0];
    }

    //! \brief Evaluate Jacobian of all shape functions
    inline void
    evaluateJacobian (const typename Traits::DomainType& in,         // position
                      std::vector<typename Traits::JacobianType>& out) const      // return value
    {
      out.resize(3);

      out[0][0][0] = -1;
      out[1][0][0] = 4-8*in[0];
      out[2][0][0] =  1;
    }

    //! \brief Evaluate partial derivatives of all shape functions
    void partial (const std::array<unsigned int, 1>& order,
                  const typename Traits::DomainType& in,         // position
                  std::vector<typename Traits::RangeType>& out) const      // return value
    {
      auto totalOrder = order[0];
      if (totalOrder == 0) {
        evaluateFunction(in, out);
      } else if (totalOrder == 1) {
        out.resize(size());
        out[0] = -1;
        out[1] = 4-8*in[0];
        out[2] =  1;
      } else if (totalOrder == 2) {
        out.resize(size());
        out[0] =  0;
        out[1] = -8;
        out[2] =  0;
      } else {
        out.resize(size());
        out[0] = out[1] = out[2] = 0;
      }
    }

    /** \brief Polynomial order of the shape functions  (2, in this case)
     */
    unsigned int order () const
    {
      return 2;
    }

  };

  /**@ingroup LocalBasisImplementation
     \brief Hierarchical P2 basis in 2d.

     The shape functions are associated to the following points:

     The functions are associated to points by:

     - \f$ f_0 \f$ ~ (0.0, 0.0)
     - \f$ f_1 \f$ ~ (0.5, 0.0)
     - \f$ f_2 \f$ ~ (1.0, 0.0)
     - \f$ f_3 \f$ ~ (0.0, 0.5)
     - \f$ f_4 \f$ ~ (0.5, 0.5)
     - \f$ f_5 \f$ ~ (0.0, 1.0)

     \tparam D Type to represent the field in the domain.
     \tparam R Type to represent the field in the range.

     \nosubgrouping
   */
  template<class D, class R>
  class HierarchicalSimplexP2LocalBasis<D,R,2>
  {
  public:
    //! \brief export type traits for function signature
    typedef LocalBasisTraits<D,2,Dune::FieldVector<D,2>,R,1,Dune::FieldVector<R,1>,
        Dune::FieldMatrix<R,1,2> > Traits;

    //! \brief number of shape functions
    unsigned int size () const
    {
      return 6;
    }

    //! \brief Evaluate all shape functions
    inline void evaluateFunction (const typename Traits::DomainType& in,
                                  std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(6);

      out[0] = 1 - in[0] - in[1];
      out[1] = 4*in[0]*(1-in[0]-in[1]);
      out[2] = in[0];
      out[3] = 4*in[1]*(1-in[0]-in[1]);
      out[4] = 4*in[0]*in[1];
      out[5] = in[1];

    }

    //! \brief Evaluate Jacobian of all shape functions
    inline void
    evaluateJacobian (const typename Traits::DomainType& in,         // position
                      std::vector<typename Traits::JacobianType>& out) const      // return value
    {
      out.resize(6);

      out[0][0][0] = -1;                    out[0][0][1] = -1;
      out[1][0][0] =  4-8*in[0]-4*in[1];    out[1][0][1] = -4*in[0];
      out[2][0][0] =  1;                    out[2][0][1] =  0;
      out[3][0][0] = -4*in[1];              out[3][0][1] =  4-4*in[0]-8*in[1];
      out[4][0][0] =  4*in[1];              out[4][0][1] =  4*in[0];
      out[5][0][0] =  0;                    out[5][0][1] =  1;
    }

    //! \brief Evaluate partial derivatives of all shape functions
    void partial (const std::array<unsigned int, 2>& order,
                  const typename Traits::DomainType& in,         // position
                  std::vector<typename Traits::RangeType>& out) const      // return value
    {
      auto totalOrder = std::accumulate(order.begin(), order.end(), 0);
      if (totalOrder == 0) {
        evaluateFunction(in, out);
      } else if (totalOrder == 1) {
        out.resize(size());
        auto const direction = std::distance(order.begin(), std::find(order.begin(), order.end(), 1));

        switch (direction) {
        case 0:
          out[0] = -1;
          out[1] =  4-8*in[0]-4*in[1];
          out[2] =  1;
          out[3] = -4*in[1];
          out[4] =  4*in[1];
          out[5] =  0;
          break;
        case 1:
          out[0] = -1;
          out[1] = -4*in[0];
          out[2] =  0;
          out[3] =  4-4*in[0]-8*in[1];
          out[4] =  4*in[0];
          out[5] =  1;
          break;
        default:
          DUNE_THROW(RangeError, "Component out of range.");
        }
      } else {
        DUNE_THROW(NotImplemented, "Desired derivative order is not implemented");
      }
    }

    /** \brief Polynomial order of the shape functions  (2 in this case)
     */
    unsigned int order () const
    {
      return 2;
    }

  };

  /**@ingroup LocalBasisImplementation
     \brief Hierarchical P2 basis in 3d.

     The shape functions are associated to the following points:

     The functions are associated to points by:

     - \f$ f_0 \f$ ~ (0.0, 0.0, 0.0)
     - \f$ f_1 \f$ ~ (0.5, 0.0, 0.0)
     - \f$ f_2 \f$ ~ (1.0, 0.0, 0.0)
     - \f$ f_3 \f$ ~ (0.0, 0.5, 0.0)
     - \f$ f_4 \f$ ~ (0.5, 0.5, 0.0)
     - \f$ f_5 \f$ ~ (0.0, 1.0, 0.0)
     - \f$ f_6 \f$ ~ (0.0, 0.0, 0.5)
     - \f$ f_7 \f$ ~ (0.5, 0.0, 0.5)
     - \f$ f_8 \f$ ~ (0.0, 0.5, 0.5)
     - \f$ f_9 \f$ ~ (0.0, 0.0, 1.0)

     \tparam D Type to represent the field in the domain.
     \tparam R Type to represent the field in the range.

     \nosubgrouping
   */
  template<class D, class R>
  class HierarchicalSimplexP2LocalBasis<D,R,3>
  {
  public:
    //! \brief export type traits for function signature
    typedef LocalBasisTraits<D,3,Dune::FieldVector<D,3>,R,1,Dune::FieldVector<R,1>,
        Dune::FieldMatrix<R,1,3> > Traits;

    //! \brief number of shape functions
    unsigned int size () const
    {
      return 10;
    }

    //! \brief Evaluate all shape functions
    void evaluateFunction (const typename Traits::DomainType& in,
                           std::vector<typename Traits::RangeType>& out) const
    {
      out.resize(10);

      out[0] = 1 - in[0] - in[1] - in[2];
      out[1] = 4 * in[0] * (1 - in[0] - in[1] - in[2]);
      out[2] = in[0];
      out[3] = 4 * in[1] * (1 - in[0] - in[1] - in[2]);
      out[4] = 4 * in[0] * in[1];
      out[5] = in[1];
      out[6] = 4 * in[2] * (1 - in[0] - in[1] - in[2]);
      out[7] = 4 * in[0] * in[2];
      out[8] = 4 * in[1] * in[2];
      out[9] = in[2];
    }

    //! \brief Evaluate Jacobian of all shape functions
    void evaluateJacobian (const typename Traits::DomainType& in,         // position
                           std::vector<typename Traits::JacobianType>& out) const      // return value
    {
      out.resize(10);

      out[0][0][0] = -1;                           out[0][0][1] = -1;                            out[0][0][2] = -1;
      out[1][0][0] =  4-8*in[0]-4*in[1]-4*in[2];   out[1][0][1] = -4*in[0];                      out[1][0][2] = -4*in[0];
      out[2][0][0] =  1;                           out[2][0][1] =  0;                            out[2][0][2] =  0;
      out[3][0][0] = -4*in[1];                     out[3][0][1] =  4-4*in[0]-8*in[1]-4*in[2];    out[3][0][2] = -4*in[1];
      out[4][0][0] =  4*in[1];                     out[4][0][1] =  4*in[0];                      out[4][0][2] =  0;
      out[5][0][0] =  0;                           out[5][0][1] =  1;                            out[5][0][2] =  0;
      out[6][0][0] = -4*in[2];                     out[6][0][1] = -4*in[2];                      out[6][0][2] =  4-4*in[0]-4*in[1]-8*in[2];
      out[7][0][0] =  4*in[2];                     out[7][0][1] =  0;                            out[7][0][2] =  4*in[0];
      out[8][0][0] =  0;                           out[8][0][1] =  4*in[2];                      out[8][0][2] =  4*in[1];
      out[9][0][0] =  0;                           out[9][0][1] =  0;                            out[9][0][2] =  1;
    }

    //! \brief Evaluate partial derivatives of all shape functions
    void partial (const std::array<unsigned int, 3>& order,
                  const typename Traits::DomainType& in,         // position
                  std::vector<typename Traits::RangeType>& out) const      // return value
    {
      auto totalOrder = std::accumulate(order.begin(), order.end(), 0);
      if (totalOrder == 0) {
        evaluateFunction(in, out);
      } else if (totalOrder == 1) {
        out.resize(size());
        auto const direction = std::distance(order.begin(), std::find(order.begin(), order.end(), 1));

        switch (direction) {
        case 0:
          out[0] = -1;
          out[1] =  4-8*in[0]-4*in[1]-4*in[2];
          out[2] =  1;
          out[3] = -4*in[1];
          out[4] =  4*in[1];
          out[5] =  0;
          out[6] = -4*in[2];
          out[7] =  4*in[2];
          out[8] =  0;
          out[9] =  0;
          break;
        case 1:
          out[0] = -1;
          out[1] = -4*in[0];
          out[2] =  0;
          out[3] =  4-4*in[0]-8*in[1]-4*in[2];
          out[4] =  4*in[0];
          out[5] =  1;
          out[6] = -4*in[2];
          out[7] =  0;
          out[8] =  4*in[2];
          out[9] =  0;
          break;
        case 2:
          out[0] = -1;
          out[1] = -4*in[0];
          out[2] =  0;
          out[3] = -4*in[1];
          out[4] =  0;
          out[5] =  0;
          out[6] =  4-4*in[0]-4*in[1]-8*in[2];
          out[7] =  4*in[0];
          out[8] =  4*in[1];
          out[9] =  1;
          break;
        default:
          DUNE_THROW(RangeError, "Component out of range.");
        }
      } else {
        DUNE_THROW(NotImplemented, "Desired derivative order is not implemented");
      }
    }

    /** \brief Polynomial order of the shape functions (2 in this case)
     */
    unsigned int order () const
    {
      return 2;
    }

  };
}
#endif
